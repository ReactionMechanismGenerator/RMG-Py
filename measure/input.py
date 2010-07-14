#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   MEASURE - Master Equation Automatic Solver for Unimolecular REactions
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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

import logging
import quantities
quantities.set_default_units('si')
quantities.UnitQuantity('kilocalorie', 1000.0*quantities.cal, symbol='kcal')

from chempy.species import Species, TransitionState
from chempy.reaction import Reaction
from chempy.species import LennardJones as LennardJonesModel
from chempy.thermo import WilhoitModel
from chempy.states import *
from chempy.kinetics import ArrheniusModel

from network import Network

################################################################################

# The current network
network = None

# A dict of the species in the input file, enabling lookup by label
speciesDict = {}

################################################################################

def processQuantity(quantity):
    if isinstance(quantity, tuple) or isinstance(quantity, list):
        value, units = quantity
    else:
        value = quantity; units = ''
    newUnits = str(quantities.Quantity(1.0, units).simplified.units).split()[1]
    if isinstance(value, tuple) or isinstance(value, list):
        return [float(quantities.Quantity(v, units).simplified) for v in value], newUnits
    else:
        return float(quantities.Quantity(value, units).simplified), newUnits

################################################################################

def species(label='', thermo=None, E0=None, states=None, lennardJones=None):
    global speciesDict
    if E0 is not None: E0 = processQuantity(E0)[0]
    else: E0 = 0.0
    spec = Species(label=label, thermo=thermo, states=states, E0=E0, lennardJones=lennardJones)
    speciesDict[label] = spec
    logging.debug('Found species "%s"' % spec)
    
def Wilhoit(cp0, cpInf, a0, a1, a2, a3, B, H0, S0):
    cp0, units = processQuantity(cp0)
    cpInf, units = processQuantity(cpInf)
    a0, units = processQuantity(a0)
    a1, units = processQuantity(a1)
    a2, units = processQuantity(a2)
    a3, units = processQuantity(a3)
    B, units = processQuantity(B)
    H0, units = processQuantity(H0)
    S0, units = processQuantity(S0)
    return WilhoitModel(cp0, cpInf, a0, a1, a2, a3, B, H0, S0)

def States(rotationalConstants=None, symmetry=1, frequencies=None, 
  frequencyScaleFactor=1.0, hinderedRotors=None, spinMultiplicity=1):
    modes = []
    if rotationalConstants is not None:
        inertia, units = processQuantity(rotationalConstants)
        linear = len(inertia)==1
        modes.append(RigidRotor(linear, inertia, symmetry))
    if frequencies is not None: 
        frequencies, units = processQuantity(frequencies)
        frequencies = [f / 100.0 for f in frequencies]
        modes.append(HarmonicOscillator(frequencies))
    if hinderedRotors is not None: 
        for inertia, barrier, symmetry in hinderedRotors:
            inertia, units = processQuantity(inertia)
            barrier, units = processQuantity(barrier)
            modes.append(HinderedRotor(inertia, barrier, symmetry))
    return StatesModel(modes, spinMultiplicity)

def LennardJones(sigma, epsilon):
    sigma, units = processQuantity(sigma)
    epsilon, units = processQuantity(epsilon)
    return LennardJonesModel(sigma, epsilon)

def reaction(reactants, products, kinetics=None, reversible=True, transitionState=None):
    global network
    try:
        rxn = Reaction(
            reactants = [speciesDict[label] for label in reactants],
            products = [speciesDict[label] for label in products],
            reversible = reversible,
            kinetics=kinetics,
            transitionState=transitionState,
        )
    except KeyError, e:
        raise NameError('A reaction was encountered with species "%s", but that species was not found in the input file.' % e.args[0])
        
    network.pathReactions.append(rxn)
    logging.debug('Found reaction "%s"' % rxn)

def Arrhenius(A, n, Ea):
    A, units = processQuantity(A)
    n, units = processQuantity(n)
    Ea, units = processQuantity(Ea)
    return ArrheniusModel(A=A, n=n, Ea=Ea)

def TS(E0=None):
    if E0 is not None: E0 = processQuantity(E0)[0]
    return TransitionState(E0=E0)

def collisionModel(type, parameters):
    parameters = [processQuantity(p)[0] for p in parameters] 
    if type.lower() == 'single exponential down':
        logging.debug('Collision model set to single exponential down (alpha = %g kJ/mol)' % (parameters[0] / 1000.0))
    else:
        raise NameError('Invalid collision model type "%s".' % type)
    
################################################################################

def readInput(path):

    global speciesDict, network
    
    try:
        f = open(path)
    except IOError, e:
        logging.error('The input file "%s" could not be opened.' % path)
        logging.info('Check that the file exists and that you have read access.')
        return
    
    # Clear any existing loaded species
    speciesDict = {}
    # Create new network object
    network = Network()
    
    logging.info('Reading input file "%s"...' % path)

    global_context = { '__builtins__': None }
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'species': species,
        'Wilhoit': Wilhoit,
        'States': States,
        'LennardJones': LennardJones,
        'reaction': reaction,
        'Arrhenius': Arrhenius,
        'TransitionState': TS,
        'collisionModel': collisionModel,
    }
    
    try:
        exec f in global_context, local_context
    except (NameError, TypeError, SyntaxError), e:
        logging.error('The input file "%s" was invalid:' % path)
        logging.exception(e)
    finally:
        f.close()
    
    # Figure out which configurations are isomers, reactant channels, and product channels
    for rxn in network.pathReactions:
        # Sort bimolecular configurations so that we always encounter them in the
        # same order
        # The actual order doesn't matter, as long as it is consistent
        rxn.reactants.sort()
        rxn.products.sort()
        # Reactants:
        # - All unimolecular configurations are automatically isomers
        # - All bimolecular configurations are automatically reactant channels
        if len(rxn.reactants) == 1 and rxn.reactants[0] not in network.isomers:
            network.isomers.append(rxn.reactants[0])
        elif len(rxn.reactants) > 1 and rxn.reactants not in network.isomers:
            network.reactants.append(rxn.reactants)
        # Products:
        # - If reversible, the same actions are taken as for the reactants
        # - If irreversible, configurations are treated as products
        if rxn.reversible:
            if len(rxn.products) == 1 and rxn.products[0] not in network.isomers:
                network.isomers.append(rxn.products[0])
            elif len(rxn.products) > 1 and rxn.products not in network.isomers:
                network.reactants.append(rxn.products)
        elif rxn.products not in network.products:
            network.products.append(rxn.products)
    
    return network
