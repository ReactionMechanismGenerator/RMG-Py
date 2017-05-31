#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
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

def computeReactionRate(rxn, forward, T, P, coreSpeciesConcentrations): 
    """

    Computes reaction rate r as follows:

    r = k * Product(Ci^nuij, for all j)
    with:

    k = rate coefficient for rxn,
    Cij = the concentration for molecule i ,
    nuij = the stoichiometric coefficient for molecule i in reaction j.

    ...
    """

    speciesList = rxn.reactants if forward == 'reactants' else rxn.products

    totconc = 1.0
    for spc in speciesList:
        ci = coreSpeciesConcentrations[spc.label]
        if abs(ci) < CLOSE_TO_ZERO:
            return 0.
        nui = rxn.getStoichiometricCoefficient(spc, forward)
        conc = ci**nui

        totconc *= conc

    k = rxn.getRateCoefficient(T,P) if forward == 'reactants' else rxn.getReverseRateCoefficient(T,P)
    r = k * totconc

    return r


def calcRij(rxn, spc, isReactant, T, P, coreSpeciesConcentrations):
    """
    This function computes the rate of formation of species i
    through the reaction j.

    This function multiplies:
    - nu(i): stoichiometric coefficient of spc in rxn
    - r(rxn): reaction rate of rxn

    Returns a reaction rate

    Units: mol / m^3 s
    """
   
    nui = rxn.getStoichiometricCoefficient(spc, isReactant)
    sign = -1 if isReactant else 1

    forward = isReactant

    rj = computeReactionRate(rxn, forward, T, P, coreSpeciesConcentrations)

    rij = nui * sign * rj
    return rij


def calcRf(spc, reactions, reactantOrProduct, T, P, coreSpeciesConcentrations, formationOrConsumption):
    """
    Calculates the total rate of formation/consumption of species i.

    Computes the sum of the rates of formation/consumption of spc for all of 
    the reactions in which spc is a product. 

    if formationOrConsumption == 'formation', spc will be compared to the 
    products of the reaction. Else, spc will be compared to the reactants of
    the reaction.

    units of rate: mol/(m^3.s)
    """
    rate = 0.0

    for reaction in reactions:
        molecules = reaction.products if formationOrConsumption == 'formation:' else reaction.reactants
        labels = [mol.label for mol in molecules]
        if spc.label in labels:
            rij = calcRij(reaction, spc,  reactantOrProduct, T, P, coreSpeciesConcentrations)
            rate = rate + rij

    logging.debug('Rf: {rate}'.format(**locals()))

    return rate
    
def calcRfClosure(spc, reactions, reactantOrProduct, T, P, coreSpeciesConcentrations):
    """
    Closure to avoid replicating function calls to calcRf.
    """
    def myfilter(formationOrConsumption):
        return calcRf(spc, reactions, reactantOrProduct, T, P, coreSpeciesConcentrations, formationOrConsumption)
    
    return myfilter

def calcRi(spc,rij, reactions, reactantOrProduct, T, P, coreSpeciesConcentrations):
    """

    Checks whether the sign of rij to decide to compute the
    total rate of formation or consumption of spc.

    units of rate: mol/(m^3.s)
    """

    closure = calcRfClosure(spc, reactions, reactantOrProduct, T, P, coreSpeciesConcentrations)

    if rij > 0:
        return closure('formation')
    elif rij < 0:
        return closure('consumption') 
    elif np.absolute([rij]) < CLOSE_TO_ZERO:
        """Pick the largest value so that the ratio rij / RX remains small."""
        Rf = closure('formation')
        Rb = closure('consumption')

        """What happens when Rf ~ Rb <<< 1?"""
        return max(abs(Rf),abs(Rb))

def isImportant(rxn, spc, reactions, reactantOrProduct, tolerance, T, P, coreSpeciesConcentrations):
    """
    This function computes:
    - Ri = R(spc)
    - rij = r(rxn)
    - alpha = ratio of rij / Ri
    
    Range of values of alpha:
    0 <= alpha <= 1

    This function also compares alpha to a user-defined tolerance TOLERANCE.
    if alpha >= tolerance:
        this reaction is important for this species.
    else:
        this reaction is unimportant for this species.

    Returns whether or not rxn is important for spc.
    keep = True
    remove = False
    """


    rij = calcRij(rxn, spc, reactantOrProduct, T, P, coreSpeciesConcentrations) 
    Ri = calcRi(spc, rij, reactions, reactantOrProduct, T, P, coreSpeciesConcentrations)

    logging.debug("rij: {rij}, Ri: {Ri}, rxn: {rxn}, species: {spc}, reactant: {reactantOrProduct}, tol: {tolerance}"\
    .format(**locals()))

    if np.any(np.absolute([rij, Ri]) < CLOSE_TO_ZERO):
       return False

    else:
        assert Ri != 0, "rij: {0}, Ri: {1}, rxn: {2}, species: {3}, reactant: {4}".format(rij, Ri, rxn, spc, reactantOrProduct)
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
