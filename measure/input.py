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
from chempy.states import *
from chempy.kinetics import ArrheniusModel

from network import Network
from collision import SingleExponentialDownModel

################################################################################

# The current network
network = None

# A dict of the species in the input file, enabling lookup by label
speciesDict = {}

# The temperatures and pressures to consider
Tlist = None
Plist = None

# The energy grains to use
Elist = None

# The method to use
method = ''

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

def species(label='', E0=None, states=None, lennardJones=None, molecularWeight=0.0):
    global speciesDict
    if E0 is not None: E0 = processQuantity(E0)[0]
    else: E0 = 0.0
    spec = Species(label=label, states=states, E0=E0, lennardJones=lennardJones)
    spec.molecularWeight = processQuantity(molecularWeight)[0]
    speciesDict[label] = spec
    logging.debug('Found species "%s"' % spec)
    
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
        network.collisionModel = SingleExponentialDownModel(alpha=parameters[0])
        logging.debug('Collision model set to single exponential down (alpha = %g kJ/mol)' % (parameters[0] / 1000.0))
    else:
        raise NameError('Invalid collision model type "%s".' % type)

def bathGas(label):
    global network, speciesDict
    network.bathGas = speciesDict[label]

def temperatures(Tlist0=None, Tmin=None, Tmax=None, count=None):
    global Tlist
    if Tlist0 is not None: 
        Tlist0 = processQuantity(Tlist0)[0]
        Tlist = Tlist0
    elif Tmin is not None and Tmax is not None and count is not None:
        # Distribute temperatures evenly on a T^-1 domain
        Tmin = processQuantity(Tmin)[0]
        Tmax = processQuantity(Tmax)[0]
        Tlist = 1.0/numpy.linspace(1.0/Tmax, 1.0/Tmin, count)
    else:
        raise SyntaxError('Must specify either a list of temperatures or Tmin, Tmax, and count.')

def pressures(Plist0=None, Pmin=None, Pmax=None, count=None):
    global Plist
    if Plist0 is not None: 
        Plist0 = processQuantity(Plist0)[0]
        Plist = Plist0
    elif Pmin is not None and Pmax is not None and count is not None:
        # Distribute pressures evenly on a log domain
        Pmin = processQuantity(Pmin)[0]
        Pmax = processQuantity(Pmax)[0]
        Plist = 10.0 ** numpy.linspace(math.log10(Pmin), math.log10(Pmax), count)
    else:
        raise SyntaxError('Must specify either a list of pressures or Pmin, Pmax, and count.')

def energies(Emin=None, Emax=None, dE=None, count=None):
    global Elist, network
    if dE is not None or count is not None:
        dE = processQuantity(dE)[0]
        if dE is None: dE = 0.0
        if count is None: count = 0
        if Emin is not None and Emax is not None:
            Emin = processQuantity(Emin)[0]
            Emax = processQuantity(Emax)[0]
            Elist = network.getEnergyGrains(Emin, Emax, dE, count)
        else:
            Elist = (dE, count)
    else:
        raise SyntaxError('Must specify either dE or count.')

def _method(name):
    global method
    method = name

################################################################################

def readInput(path):

    global speciesDict, network, Tlist, Plist, Elist
    
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
        'States': States,
        'LennardJones': LennardJones,
        'reaction': reaction,
        'Arrhenius': Arrhenius,
        'TransitionState': TS,
        'collisionModel': collisionModel,
        'bathGas': bathGas,
        'temperatures': temperatures,
        'pressures': pressures,
        'energies': energies,
        'method': _method,
    }
    
    try:
        exec f in global_context, local_context
    except (NameError, TypeError, SyntaxError), e:
        logging.error('The input file "%s" was invalid:' % path)
        logging.exception(e)
        network = None
    finally:
        f.close()
    
    # If loading of the input file was unsuccessful for any reason,
    # then return None for everything so the program can terminate
    if network is None: return None, None, None, None
    
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
    
    # Print lots of information about the loaded network
    # In particular, we want to give all of the energies on the PES
    # This will help the user decide if the range of energies selected is 
    # appropriate when viewing the log file
    logging.debug('')
    logging.debug('========================================================================')
    logging.debug('Network Information')
    logging.debug('-------------------')
    logging.debug('Isomers:')
    for isomer in network.isomers:
        logging.debug('    {0:<48s} {1:12g} kJ/mol'.format(str(isomer), isomer.E0 / 1000.0))
    logging.debug('Reactant channels:')
    for reactants in network.reactants:
        logging.debug('    {0:<48s} {1:12g} kJ/mol'.format(' + '.join([str(spec) for spec in reactants]), sum([spec.E0 for spec in reactants])))
    logging.debug('Product channels:')
    for products in network.products:
        logging.debug('    {0:<48s} {1:12g} kJ/mol'.format(' + '.join([str(spec) for spec in products]), sum([spec.E0 for spec in products])))
    logging.debug('Path reactions:')
    for rxn in network.pathReactions:
        logging.debug('    {0:<48s} {1:12g} kJ/mol'.format(rxn, rxn.transitionState.E0))
    logging.debug('========================================================================')
    logging.debug('')
    
    # If there are no isomers, then there's nothing to do
    if len(network.isomers) == 0:
        logging.info('Could not find any unimolecular isomers based on this network, so there is nothing to do.')
        return None, None, None, None
      
    
    return network, Tlist, Plist, Elist, method
