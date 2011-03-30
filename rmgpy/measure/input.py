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

"""
Contains the :meth:`readInput()` function, used to read MEASURE input files.
The MEASURE input file format is based on Python syntax, and is in fact a
valid Python file. This allows us to easily load the file by defining functions
and connecting them to the constructs in the input file. If the reading of the
input file is unsuccessful, an :class:`InputError` will be raised.
"""

import logging
import quantities
quantities.set_default_units('si')
quantities.UnitQuantity('kilocalorie', 1000.0*quantities.cal, symbol='kcal')
quantities.UnitQuantity('kilojoule', 1000.0*quantities.J, symbol='kJ')

from rmgpy.chem.species import Species, TransitionState
from rmgpy.chem.reaction import Reaction
from rmgpy.chem.species import LennardJones
from rmgpy.chem.molecule import Molecule
from rmgpy.chem.states import *
from rmgpy.chem.kinetics import Arrhenius
from rmgpy.chem.thermo import *

from network import Network
from collision import SingleExponentialDownModel

################################################################################

# The current network
network = None

# A dict of the species in the input file, enabling lookup by label
speciesDict = {}

# The temperatures and pressures to consider
Tlist = None; Tparams = None
Plist = None; Pparams = None

# The energy grains to use
Elist = None

# The method to use
method = ''

# The interpolation model to use
model = ['']

################################################################################

class InputError(Exception):
    """
    An exception that is raised when a MEASURE input file is invalid for some
    reason. Pass a string detailing the circumstances of the exceptional
    behavior.
    """
    pass

################################################################################

def species(label='', E0=None, states=None, thermo=None, lennardJones=None, molecularWeight=0.0, SMILES='', InChI=''):
    global speciesDict
    if label == '': raise InputError('Missing "label" attribute in species() block.')
    spec = Species(label=label, states=states, thermo=thermo, E0=E0, molecularWeight=molecularWeight, lennardJones=lennardJones)
    if InChI != '':
        spec.molecule = [Molecule(InChI=InChI)]
    elif SMILES != '':
        spec.molecule = [Molecule(SMILES=SMILES)]
    speciesDict[label] = spec
    logging.debug('Found species "%s"' % spec)
    # If the molecular weight was not specified but the structure was, then
    # get the molecular weight from the structure
    if spec.molecularWeight == 0.0 and spec.molecule is not None and len(spec.molecule) > 0:
        spec.molecularWeight = spec.molecule[0].getMolecularWeight()

def States(rotationalConstants=None, symmetry=1, frequencies=None, 
  frequencyScaleFactor=1.0, hinderedRotors=None, spinMultiplicity=1):
    modes = []
    if rotationalConstants is not None:
        inertia = rotationalConstants
        linear = len(inertia)==1
        modes.append(RigidRotor(linear, inertia, symmetry))
    if frequencies is not None: 
        modes.append(HarmonicOscillator(frequencies))
    if hinderedRotors is not None: 
        for inertia, barrier, symmetry in hinderedRotors:
            modes.append(HinderedRotor(inertia, barrier, symmetry))
    return StatesModel(modes, spinMultiplicity)

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
        raise InputError('A reaction was encountered with species "%s", but that species was not found in the input file.' % e.args[0])
        
    network.pathReactions.append(rxn)
    logging.debug('Found reaction "%s"' % rxn)

def collisionModel(type, parameters, bathGas):
    global network, speciesDict
    if type.lower() == 'single exponential down':

        # Process parameters, making sure we have a valid set
        if len(parameters) == 1:
            if 'alpha' not in parameters:
                raise InputError('Must specify either "alpha" or ("alpha0","T0","n") as parameters for SingleExponentialDownModel.')
            alpha0 = constants.Quantity(parameters['alpha']).value
            T0 = 1000.0
            n = 0.0
        elif len(parameters) == 3:
            if 'alpha0' not in parameters or 'T0' not in parameters or 'n' not in parameters:
                raise InputError('Must specify either "alpha" or ("alpha0","T0","n") as parameters for SingleExponentialDownModel.')
            alpha0 = constants.Quantity(parameters['alpha0']).value
            T0 = constants.Quantity(parameters['T0']).value
            n = constants.Quantity(parameters['n']).value
        else:
            raise InputError('Must specify either "alpha" or ("alpha0","T0","n") as parameters for SingleExponentialDownModel.')

        # Create the collision model object
        network.collisionModel = SingleExponentialDownModel(alpha0=alpha0, T0=T0, n=n)
        logging.debug('Collision model set to single exponential down')
        
    else:
        raise NameError('Invalid collision model type "%s".' % type)
    # Set bath gas composition
    network.bathGas = {}
    for key, value in bathGas.iteritems():
        network.bathGas[speciesDict[key]] = float(value)
    # Normalize bath gas composition
    for key in network.bathGas:
        network.bathGas[key] /= sum(network.bathGas.values())
    
def temperatures(Tlist0=None, Tmin=None, Tmax=None, count=None):
    global Tlist, Tparams
    if Tlist0 is not None:
        # We've been provided a list of specific temperatures to use
        Tlist = constants.processQuantity(Tlist0)[0].value
        Tparams = None
    elif Tmin is not None and Tmax is not None and count is not None:
        # We've been provided a temperature range and number of temperatures to use
        # We defer choosing the actual temperatures because they depend on the
        # choice of interpolation model
        Tlist = None
        Tparams = [constants.Quantity(Tmin).value, constants.Quantity(Tmax).value, count]
    else:
        raise SyntaxError('Must specify either a list of temperatures or Tmin, Tmax, and count.')

def pressures(Plist0=None, Pmin=None, Pmax=None, count=None):
    global Plist, Pparams
    if Plist0 is not None:
        # We've been provided a list of specific pressures to use
        Plist = constants.processQuantity(Plist0).value
        Pparams = None
    elif Pmin is not None and Pmax is not None and count is not None:
        # We've been provided a pressures range and number of pressures to use
        # We defer choosing the actual pressures because they depend on the
        # choice of interpolation model
        Plist = None
        Pparams = [constants.Quantity(Pmin).value, constants.Quantity(Pmax).value, count]
    else:
        raise SyntaxError('Must specify either a list of pressures or Pmin, Pmax, and count.')

def energies(Emin=None, Emax=None, dE=None, count=None):
    global Elist, network
    if dE is not None or count is not None:
        dE = constants.Quantity(dE).value
        if dE is None: dE = 0.0
        if count is None: count = 0
        if Emin is not None and Emax is not None:
            Emin = constants.Quantity(Emin).value
            Emax = constants.Quantity(Emax).value
            Elist = network.getEnergyGrains(Emin, Emax, dE, count)
        else:
            Elist = (dE, count)
    else:
        raise InputError('Must specify either dE or count in energies() block.')

def _method(name):
    global method
    if name.lower() not in ['modified strong collision', 'reservoir state', 'chemically-significant eigenvalues', 'branching ratios']:
        raise InputError('Invalid method "%s"; see documentation for available methods.' % name)
    method = name

def interpolationModel(name, *args):
    global model
    model = [name]
    model.extend(list(args))

################################################################################

def generateThermoFromStates(species):
    """
    For a given :class:`Species` object `species` with molecular degrees of
    freedom data in its ``states`` attribute, generate a corresponding thermo
    model. A :class:`ThermoData` thermodynamics model is stored in the
    species ``thermo`` attribute, and nothing is returned.
    """
    # Do nothing if the species already has thermo data or if it does not have
    # states data
    if species.thermo is not None or species.states is None: return
    # States data must have external rotational modes
    if not any([isinstance(mode, RigidRotor) for mode in species.states.modes]):
        raise InputError('For species "%s", must specify external rotational constants to generate thermo model from states data.' % species)
    
    # Must use ThermoData because we can't rely on knowing
    # anything about the structure of the species
    from rmgpy.chem.thermo import ThermoData
    Tdata = numpy.linspace(numpy.min(Tlist), numpy.max(Tlist), 20.0)
    Cpdata = species.states.getHeatCapacities(Tdata)
    H298 = species.E0 + species.states.getEnthalpy(298)
    S298 = species.states.getEntropy(298)

    # Add in heat capacities for translational modes if missing
    if not any([isinstance(mode, Translation) for mode in species.states.modes]):
        if species.molecularWeight == 0:
            raise InputError('Molecular weight required for species "%s".' % species)
        trans = Translation(species.molecularWeight)
        Cpdata += trans.getHeatCapacities(Tdata)
        H298 += trans.getEnthalpy(298)
        S298 += trans.getEntropy(298)

    species.thermo = ThermoData(Tdata=(Tdata,"K"), Cpdata=(Cpdata,"J/(mol*K)"), H298=(H298/1000.,"kJ/mol"), S298=(S298,"J/(mol*K)"))

def getTemperaturesForModel(model, Tmin, Tmax, Tcount):
    """
    Returns an array of temperatures based on the interpolation `model`,
    minimum and maximum temperatures `Tmin` and `Tmax` in K, and the number of
    temperatures `Tcount`. For Chebyshev polynomials a Gauss-Chebyshev
    distribution is used; for all others a linear distribution on an inverse
    temperature domain is used. Note that the Gauss-Chebyshev grid does *not*
    place `Tmin` and `Tmax` at the endpoints, yet the interpolation is still
    valid up to these values.
    """
    if model[0].lower() == 'chebyshev':
        # Distribute temperatures on a Gauss-Chebyshev grid
        Tlist = numpy.zeros(Tcount, numpy.float64)
        for i in range(Tcount):
            T = -math.cos((2*i+1) * math.pi / (2*Tcount))
            T = 2.0 / ((1.0/Tmax - 1.0/Tmin) * T + 1.0/Tmax + 1.0/Tmin)
            Tlist[i] = T
    else:
        # Distribute temperatures evenly on a T^-1 domain
        Tlist = 1.0/numpy.linspace(1.0/Tmax, 1.0/Tmin, Tcount)
    return Tlist

def getPressuresForModel(model, Pmin, Pmax, Pcount):
    """
    Returns an array of pressures based on the interpolation `model`,
    minimum and maximum pressures `Pmin` and `Pmax` in Pa, and the number of
    pressures `Pcount`. For Chebyshev polynomials a Gauss-Chebyshev
    distribution is used; for all others a linear distribution on an logarithmic
    pressure domain is used. Note that the Gauss-Chebyshev grid does *not*
    place `Pmin` and `Pmax` at the endpoints, yet the interpolation is still
    valid up to these values.
    """
    if model[0].lower() == 'chebyshev':
        # Distribute pressures on a Gauss-Chebyshev grid
        Plist = numpy.zeros(Pcount, numpy.float64)
        for i in range(Pcount):
            P = -math.cos((2*i+1) * math.pi / (2*Pcount))
            P = 10**(0.5 * ((math.log10(Pmax) - math.log10(Pmin)) * P + math.log10(Pmax) + math.log10(Pmin)))
            Plist[i] = P
    else:
        # Distribute pressures evenly on a log domain
        Plist = 10.0 ** numpy.linspace(math.log10(Pmin), math.log10(Pmax), Pcount)
    return Plist

def readInput(path):
    """
    Reads a MEASURE input file from location `path` on disk. The input file
    format is described in the :ref:`measureusersguide`. Returns a number of
    quantities:

    * The :class:`Network` object representing the unimolecular reaction network

    * The list of temperatures in K to be used in the master equation
      calculation

    * The list of pressures in Pa to be used in the master equation
      calculation

    * A tuple containing the maximum energy grain size in J/mol and the
      minimum number of energy grains to use in the master equation calculation;
      whichever of these results in more energy grains

    * The approximate method to use to estimate the phenomenological rate
      coefficients :math:`k(T,P)`

    * The interpolation model to fit the estimated :math:`k(T,P)` values to

    * The minimum temperature in K at which the fitted interpolation model is
      valid; this is *not* necessarily equal to ``min(Tlist)``

    * The maximum temperature in K at which the fitted interpolation model is
      valid; this is *not* necessarily equal to ``max(Tlist)``

    * The minimum temperature in Pa at which the fitted interpolation model is
      valid; this is *not* necessarily equal to ``min(Plist)``

    * The maximum temperature in Pa at which the fitted interpolation model is
      valid; this is *not* necessarily equal to ``max(Plist)``
    
    """

    global speciesDict, network, Tlist, Tparams, Plist, Pparams, Elist, method, model
    
    try:
        f = open(path)
    except IOError, e:
        logging.error('The input file "%s" could not be opened.' % path)
        logging.info('Check that the file exists and that you have read access.')
        return None
    
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
        'TransitionState': TransitionState,
        'collisionModel': collisionModel,
        'temperatures': temperatures,
        'pressures': pressures,
        'energies': energies,
        'method': _method,
        'interpolationModel': interpolationModel,
        'ThermoData': ThermoData,
        'Wilhoit': Wilhoit,
        'MultiNASA': MultiNASA,
        'NASA': NASA,
    }
    
    try:
        exec f in global_context, local_context
    
    
        # If loading of the input file was unsuccessful for any reason,
        # then return None so the program can terminate
        if network is None: return None

        # Determine temperature grid if not yet known
        if Tparams is not None and Tlist is None:
            Tmin, Tmax, Tcount = Tparams
            Tlist = getTemperaturesForModel(model, Tmin, Tmax, Tcount)
        elif Tparams is None and Tlist is not None:
            Tmin = min(Tlist); Tmax = max(Tlist); Tcount = len(Tlist)
        else:
            raise InputError('No temperature() block found.')
        Tlist = numpy.array(Tlist, numpy.float64)

        # Determine pressure grid if not yet known
        if Pparams is not None and Plist is None:
            Pmin, Pmax, Pcount = Pparams
            Plist = getPressuresForModel(model, Pmin, Pmax, Pcount)
        elif Pparams is None and Plist is not None:
            Pmin = min(Plist); Pmax = max(Plist); Pcount = len(Plist)
        else:
            raise InputError('No pressure() block found.')
        Plist = numpy.array(Plist, numpy.float64)

        # Check that we have energy grain information
        if Elist is None:
            raise InputError('No energies() block found.')

        # Check that a method was specified
        if method == '':
            raise InputError('No method() block found.')

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
                network.isomers.insert(0, rxn.reactants[0])
            elif len(rxn.reactants) > 1 and rxn.reactants not in network.reactants:
                network.reactants.append(rxn.reactants)
            # Products:
            # - If reversible, the same actions are taken as for the reactants
            # - If irreversible, configurations are treated as products
            if rxn.reversible:
                if len(rxn.products) == 1 and rxn.products[0] not in network.isomers:
                    network.isomers.insert(0, rxn.products[0])
                elif len(rxn.products) > 1 and rxn.products not in network.reactants:
                    network.reactants.append(rxn.products)
            elif rxn.products not in network.products:
                network.products.append(rxn.products)
        
        # For each configuration with states data but not thermo data,
        # calculate the thermo data
        for isomer in network.isomers:
            if isomer.thermo is None and isomer.states is not None:
                generateThermoFromStates(isomer)
        for reactants in network.reactants:
            for spec in reactants:
                if spec.thermo is None and spec.states is not None:
                    generateThermoFromStates(spec)
        for products in network.products:
            for spec in products:
                if spec.thermo is None and spec.states is not None:
                    generateThermoFromStates(spec)

        # Check that we have the right data for each configuration
        # Use a string to store the errors so that the user can see all of
        # the mistakes at once
        errorString = ''
        
        for isomer in network.isomers:
            errorString0 = ''
            # All isomers must have states data
            if isomer.states is None:
                errorString0 += '    Required molecular degree of freedom data was not provided.\n'
            # All isomers must have collision parameters
            if isomer.lennardJones is None:
                errorString0 += '    Required Lennard-Jones parameters were not provided.\n'
            if isomer.molecularWeight == 0:
                errorString0 += '    Required molecular weight was not provided.\n'
            if errorString0 != '':
                errorString = 'For unimolecular isomer "%s":\n%s' % (isomer, errorString0)

        for rxn in network.pathReactions:
            errorString0 = ''
            # All reactions must have either high-pressure-limit kinetics or transition state data
            if rxn.kinetics is None and rxn.transitionState.states is None:
                errorString0 += '    Unable to determine microcanonical rate k(E); you must specify either the high-P kinetics or transition state molecular degrees of freedom.\n' % spec

            if rxn.reversible:
                # All reversible reactions must have thermo for both reactants and products
                for spec in rxn.reactants:
                    if spec.thermo is None:
                        errorString0 += '    Unable to determine thermo data for reactant "%s"; you must specify a thermodynamics model.\n' % spec
                for spec in rxn.products:
                    if spec.thermo is None:
                        errorString0 += '    Unable to determine thermo data for product "%s"; you must specify a thermodynamics model.\n' % spec
            else:
                # All irreversible reactions must have states data for reactants
                for spec in rxn.reactants:
                    if spec.states is None:
                        errorString0 += '    Required molecular degree of freedom data for reactant "%s" was not provided.\n' % spec

            if errorString0 != '':
                errorString += 'For path reaction "%s":\n%s' % (rxn, errorString0)

        if errorString != '':
            raise InputError(errorString)

    except InputError, e:
        logging.error('The input file "%s" was invalid:' % path)
        logging.info(e)
        return None

    # Print lots of information about the loaded network
    # In particular, we want to give all of the energies on the PES
    # This will help the user decide if the range of energies selected is 
    # appropriate when viewing the log file
    network.printSummary(level=logging.INFO)
    
    # If there are no isomers, then there's nothing to do
    if len(network.isomers) == 0:
        logging.info('Could not find any unimolecular isomers based on this network, so there is nothing to do.')
        return None
      
    return network, Tlist, Plist, Elist, method, model, Tmin, Tmax, Pmin, Pmax

