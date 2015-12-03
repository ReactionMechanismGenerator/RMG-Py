#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import numpy as np
from rmgpy.scoop_framework.util import logger as logging

CLOSE_TO_ZERO = 1E-20

def compute_reaction_rate(rxn_j, forward, T, P, coreSpeciesConcentrations): 
    """

    Computes reaction rate r as follows:

    r = k * Product(Ci^nu_ij, for all j)
    with:

    k = rate coefficient for rxn_j,
    Cij = the concentration for molecule i ,
    nu_ij = the stoichiometric coefficient for molecule i in reaction j.

    ...
    """

    species_list = rxn_j.reactants if forward == 'reactants' else rxn_j.products

    totconc = 1.0
    for spc_i in species_list:
        ci = coreSpeciesConcentrations[spc_i.label]
        if abs(ci) < CLOSE_TO_ZERO:
            return 0.
        nu_i = rxn_j.get_stoichiometric_coefficient(spc_i, forward)
        conc = ci**nu_i

        totconc *= conc

    k = rxn_j.getRateCoefficient(T,P) if forward == 'reactants' else rxn_j.getReverseRateCoefficient(T,P)
    r = k * totconc

    return r


def calc_rij(rxn_j, spc_i, isReactant, T, P, coreSpeciesConcentrations):
    """
    This function computes the rate of formation of species i
    through the reaction j.

    This function multiplies:
    - nu(i): stoichiometric coefficient of spc_i in rxn_j
    - r(rxn_j): reaction rate of rxn_j

    Returns a reaction rate

    Units: mol / m^3 s
    """
   
    nu_i = rxn_j.get_stoichiometric_coefficient(spc_i, isReactant)
    sign = -1 if isReactant else 1

    forward = isReactant

    r_j = compute_reaction_rate(rxn_j, forward, T, P, coreSpeciesConcentrations)

    rij = nu_i * sign * r_j
    return rij


def calc_Rf(spc_i, reactions, reactant_or_product, T, P, coreSpeciesConcentrations, formation_or_consumption):
    """
    Calculates the total rate of formation/consumption of species i.

    Computes the sum of the rates of formation/consumption of spc_i for all of 
    the reactions in which spc_i is a product. 

    if formation_or_consumption == 'formation', spc_i will be compared to the 
    products of the reaction. Else, spc_i will be compared to the reactants of
    the reaction.

    units of rate: mol/(m^3.s)
    """
    rate = 0.0

    for reaction in reactions:
        molecules = reaction.products if formation_or_consumption == 'formation:' else reaction.reactants
        labels = [mol.label for mol in molecules]
        if spc_i.label in labels:
            rij = calc_rij(reaction, spc_i,  reactant_or_product, T, P, coreSpeciesConcentrations)
            rate = rate + rij

    logging.debug('Rf: {rate}'.format(**locals()))

    return rate
    
def calc_Rf_closure(spc_i, reactions, reactant_or_product, T, P, coreSpeciesConcentrations):
    """
    Closure to avoid replicating function calls to calc_Rf.
    """
    def myfilter(formation_or_consumption):
        return calc_Rf(spc_i, reactions, reactant_or_product, T, P, coreSpeciesConcentrations, formation_or_consumption)
    
    return myfilter

def calc_Ri(spc_i,rij, reactions, reactant_or_product, T, P, coreSpeciesConcentrations):
    """

    Checks whether the sign of rij to decide to compute the
    total rate of formation or consumption of spc_i.

    units of rate: mol/(m^3.s)
    """

    f_closure = calc_Rf_closure(spc_i, reactions, reactant_or_product, T, P, coreSpeciesConcentrations)

    if rij > 0:
        return f_closure('formation')
    elif rij < 0:
        return f_closure('consumption') 
    elif np.absolute([rij]) < CLOSE_TO_ZERO:
        """Pick the largest value so that the ratio rij / RX remains small."""
        Rf = f_closure('formation')
        Rb = f_closure('consumption')

        """What happens when Rf ~ Rb <<< 1?"""
        return max(abs(Rf),abs(Rb))

def isImportant(rxn, species_i, reactions, reactant_or_product, tolerance, T, P, coreSpeciesConcentrations):
    """
    This function computes:
    - Ri = R(species_i)
    - rij = r(rxn)
    - alpha = ratio of rij / Ri
    
    Range of values of alpha:
    0 <= alpha <= 1

    This function also compares alpha to a user-defined tolerance TOLERANCE.
    if alpha >= tolerance:
        this reaction is important for this species.
    else:
        this reaction is unimportant for this species.

    Returns whether or not rxn is important for species_i.
    keep = True
    remove = False
    """


    rij = calc_rij(rxn, species_i, reactant_or_product, T, P, coreSpeciesConcentrations) 
    Ri = calc_Ri(species_i, rij, reactions, reactant_or_product, T, P, coreSpeciesConcentrations)

    logging.debug("rij: {rij}, Ri: {Ri}, rxn: {rxn}, species: {species_i}, reactant: {reactant_or_product}, tol: {tolerance}"\
    .format(**locals()))

    if np.any(np.absolute([rij, Ri]) < CLOSE_TO_ZERO):
       return False

    else:
        assert Ri != 0, "rij: {0}, Ri: {1}, rxn: {2}, species: {3}, reactant: {4}".format(rij, Ri, rxn, species_i, reactant_or_product)
        alpha = rij / Ri
        if alpha < 0: return False


    if alpha > tolerance :
        """
        If both values are very close to 1, then the comparison of alpha and the tolerance
        might sometimes return an unexpected value.

        When we set the tolerance to a value of 1, we want all the reactions to be unimportant,
        regardless of the value of alpha.

        """
        if np.allclose([tolerance, alpha], [1.0, 1.0]):
            return False
            
        return True
        #where tolerance is user specified tolerance
 
    elif alpha <= tolerance:
        return False