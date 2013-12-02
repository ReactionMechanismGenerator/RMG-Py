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

from rmgpy.quantity import Quantity
from rmgpy.species import Species, TransitionState
from rmgpy.transport import TransportData
from rmgpy.reaction import Reaction
from rmgpy.molecule import Molecule
from rmgpy.statmech import *
from rmgpy.kinetics import Arrhenius, PDepArrhenius, Chebyshev
from rmgpy.thermo import *

from network import Network
from collision import SingleExponentialDown

################################################################################

# The current MEASURE job object
measure = None

# The current network
network = None

# A dict of the species in the input file, enabling lookup by label
speciesDict = {}

################################################################################

class InputError(Exception):
    """
    An exception that is raised when a MEASURE input file is invalid for some
    reason. Pass a string detailing the circumstances of the exceptional
    behavior.
    """
    pass

################################################################################

def species(label='', E0=None, states=None, thermo=None, transportData=None, molecularWeight=0.0, collisionModel=None, SMILES='', InChI=''):
    global speciesDict
    if label == '': raise InputError('Missing "label" attribute in species() block.')
    spec = Species(label=label, states=states, thermo=thermo, E0=E0, molecularWeight=molecularWeight, transportData=transportData, collisionModel=collisionModel)
    if InChI != '':
        spec.molecule = [Molecule(InChI=InChI)]
    elif SMILES != '':
        spec.molecule = [Molecule(SMILES=SMILES)]
    speciesDict[label] = spec
    logging.debug('Found species "{0}"'.format(spec))
    # If the molecular weight was not specified but the structure was, then
    # get the molecular weight from the structure
    if spec.molecularWeight.value_si == 0.0 and spec.molecule is not None and len(spec.molecule) > 0:
        spec.molecularWeight = Quantity(spec.molecule[0].getMolecularWeight(),"kg/mol")
    
def States(translation=None, rotations=None, vibrations=None, torsions=None, frequencyScaleFactor=1.0, spinMultiplicity=1):
    modes = []
    if translation is not None:
        if isinstance(translation, list) or isinstance(translation, tuple):
            modes.extend(translation)
        else:
            modes.append(translation)
    if rotations is not None:
        if isinstance(rotations, list) or isinstance(rotations, tuple):
            modes.extend(rotations)
        else:
            modes.append(rotations)
    if vibrations is not None:
        if isinstance(vibrations, list) or isinstance(vibrations, tuple):
            modes.extend(vibrations)
        else:
            modes.append(vibrations)
    if torsions is not None:
        if isinstance(torsions, list) or isinstance(torsions, tuple):
            modes.extend(torsions)
        else:
            modes.append(torsions)
    for mode in modes:
        if isinstance(mode, HarmonicOscillator):
            mode.frequencies.value_si *= frequencyScaleFactor
    return StatesModel(modes, spinMultiplicity)

def addIsomer(species):
    global network
    network.isomers.append(species)

def addReactants(speciesA, speciesB):
    global network
    network.reactants.append([speciesA, speciesB])

def reaction(reactants, products, kinetics=None, transitionState=None):
    global network
    try:
        rxn = Reaction(
            reactants = reactants,
            products = products,
            kinetics=kinetics,
            transitionState=transitionState,
        )
    except KeyError, e:
        raise InputError('A reaction was encountered with species "{0}", but that species was not found in the input file.'.format(e.args[0]))
        
    network.pathReactions.append(rxn)
    logging.debug('Found reaction "{0}"'.format(rxn))

def pdepreaction(reactants, products, kinetics=None):
    global network
    try:
        rxn = Reaction(
            reactants = reactants,
            products = products,
            kinetics=kinetics,
        )
    except KeyError, e:
        raise InputError('A reaction was encountered with species "{0}", but that species was not found in the input file.'.format(e.args[0]))
        
    network.netReactions.append(rxn)
    logging.debug('Found pdepreaction "{0}"'.format(rxn))

def temperatures(Tlist=None, Tmin=None, Tmax=None, count=None):
    global measure
    if Tlist is not None:
        # We've been provided a list of specific temperatures to use
        measure.Tlist = Quantity(Tlist)
        measure.Tmin = Quantity(numpy.min(measure.Tlist.value_si),"K")
        measure.Tmax = Quantity(numpy.max(measure.Tlist.value_si),"K")
    elif Tmin is not None and Tmax is not None and count is not None:
        # We've been provided a temperature range and number of temperatures to use
        # We defer choosing the actual temperatures because they depend on the
        # choice of interpolation model
        measure.Tmin = Quantity(Tmin)
        measure.Tmax = Quantity(Tmax)
        measure.Tcount = count
    else:
        raise SyntaxError('Must specify either a list of temperatures or Tmin, Tmax, and count.')

def pressures(Plist=None, Pmin=None, Pmax=None, count=None):
    global measure
    if Plist is not None:
        # We've been provided a list of specific pressures to use
        measure.Plist = Quantity(Plist)
        measure.Pmin = Quantity(numpy.min(measure.Plist.value_si),"Pa")
        measure.Pmax = Quantity(numpy.max(measure.Plist.value_si),"Pa")
    elif Pmin is not None and Pmax is not None and count is not None:
        # We've been provided a pressures range and number of pressures to use
        # We defer choosing the actual pressures because they depend on the
        # choice of interpolation model
        measure.Pmin = Quantity(Pmin)
        measure.Pmax = Quantity(Pmax)
        measure.Pcount = count
    else:
        raise SyntaxError('Must specify either a list of pressures or Pmin, Pmax, and count.')

def energies(Emin=None, Emax=None, dE=None, count=None):
    global measure
    if dE is not None or count is not None:
        if dE is not None:
            measure.grainSize = Quantity(dE)
        if count is not None:
            measure.grainCount = count
        if Emin is not None and Emax is not None:
            measure.Emin = Quantity(Emin)
            measure.Emax = Quantity(Emax)
    else:
        raise InputError('Must specify either dE or count in energies() block.')

def _method(name):
    global measure
    if name.lower() not in ['modified strong collision', 'reservoir state', 'chemically-significant eigenvalues', 'branching ratios']:
        raise InputError('Invalid method "{0}"; see documentation for available methods.'.format(name))
    measure.method = name

def interpolationModel(name, *args):
    global measure
    measure.model = [name]
    measure.model.extend(list(args))

################################################################################

def generateThermoFromStates(species, Tlist):
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
    if not any([isinstance(mode, RigidRotor) for mode in species.states.modes]) and len(species.molecule[0].atoms) > 1:
        raise InputError('For species "{0}", must specify external rotational constants to generate thermo model from states data.'.format(species))
    
    # Must use ThermoData because we can't rely on knowing
    # anything about the structure of the species
    Tdata = numpy.linspace(numpy.min(Tlist), numpy.max(Tlist), 20.0)
    Cpdata = species.states.getHeatCapacities(Tdata)
    H298 = species.E0.value_si + species.states.getEnthalpy(298)
    S298 = species.states.getEntropy(298)

    # Add in heat capacities for translational modes if missing
    if not any([isinstance(mode, Translation) for mode in species.states.modes]):
        if species.molecularWeight == 0:
            raise InputError('Molecular weight required for species "{0}".'.format(species))
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

def readFile(path, measure0):
    """
    Reads a MEASURE input or output file from location `path` on disk. The file
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

    global measure, speciesDict, network, Tlist, Tparams, Plist, Pparams, Elist, method, model
    
    try:
        f = open(path)
    except IOError, e:
        logging.error('The input file "{0}" could not be opened.'.format(path))
        logging.info('Check that the file exists and that you have read access.')
        return None
    
    # Clear the job object
    measure = measure0
    measure.clear()
    # Clear any existing loaded species
    speciesDict = {}
    # Create new network object
    network = Network()
    
    title = 'Untitled'
    description = ''
    
    logging.info('Reading input file "{0}"...'.format(path))

    global_context = { '__builtins__': None }
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'species': species,
        'States': States,
        'Translation': Translation,
        'RigidRotor': RigidRotor,
        'HarmonicOscillator': HarmonicOscillator,
        'HinderedRotor': HinderedRotor,
        'TransportData': TransportData,
        'isomer': addIsomer,
        'reactants': addReactants,
        'reaction': reaction,
        'Arrhenius': Arrhenius,
        'TransitionState': TransitionState,
        'pdepreaction': pdepreaction,
        'Chebyshev': Chebyshev,
        'PDepArrhenius': PDepArrhenius,
        'SingleExponentialDown': SingleExponentialDown,
        'temperatures': temperatures,
        'pressures': pressures,
        'energies': energies,
        'method': _method,
        'interpolationModel': interpolationModel,
        'ThermoData': ThermoData,
        'Wilhoit': Wilhoit,
        'MultiNASA': MultiNASA,
        'NASA': NASA,
        'title': title,
        'description': description,
    }
    
    try:
        exec f in global_context, local_context
    
        # If loading of the input file was unsuccessful for any reason,
        # then return so the program can terminate
        if network is None: return

        # Set title and description of network
        network.title = local_context['title']
        network.description = local_context['description'].strip()
        # Log title and description
        logging.info('')
        logging.info('Network title: {0}'.format(network.title))
        logging.info('Network description: ')
        logging.info(network.description)
        
        measure.network = network
        
        # Determine temperature grid if not yet known
        if measure.Tlist is None and measure.Tmin is not None and measure.Tmax is not None and measure.Tcount is not None:
            measure.Tlist = Quantity(getTemperaturesForModel(measure.model, measure.Tmin.value_si, measure.Tmax.value_si, measure.Tcount),"K")
        elif measure.Tmin is not None and measure.Tmax is not None and measure.Tcount is not None:
            pass
        else:
            raise InputError('No temperature() block found.')
        
        # Determine pressure grid if not yet known
        if measure.Plist is None and measure.Pmin is not None and measure.Pmax is not None and measure.Pcount is not None:
            measure.Plist = Quantity(getPressuresForModel(measure.model, measure.Pmin.value_si, measure.Pmax.value_si, measure.Pcount),"Pa")
        elif measure.Pmin is not None and measure.Pmax is not None and measure.Pcount is not None:
            pass
        else:
            raise InputError('No pressure() block found.')
        
        # Check that we have energy grain information
        if measure.grainSize is None and measure.grainCount is None:
            raise InputError('No energies() block found.')

        # Check that a method was specified
        if measure.method is None:
            raise InputError('No method() block found.')

        # Convert string labels to Species objects
        network.isomers = [speciesDict[label] for label in network.isomers]
        network.reactants = [[speciesDict[label] for label in reactants] for reactants in network.reactants]
        
        # Set bath gas composition
        network.bathGas = {}
        for label, value in local_context['bathGas'].iteritems():
            network.bathGas[speciesDict[label]] = float(value)
        # Normalize bath gas composition
        for key in network.bathGas:
            network.bathGas[key] /= sum(network.bathGas.values())
        
        for reactants in network.reactants:
            reactants.sort()
        for rxn in network.pathReactions:
            rxn.reactants = [speciesDict[label] for label in rxn.reactants]
            rxn.reactants.sort()
            rxn.products = [speciesDict[label] for label in rxn.products]
            rxn.products.sort()
        for rxn in network.netReactions:
            rxn.reactants = [speciesDict[label] for label in rxn.reactants]
            rxn.reactants.sort()
            rxn.products = [speciesDict[label] for label in rxn.products]
            rxn.products.sort()
            
        # Figure out which configurations are isomers, reactant channels, and product channels
        for rxn in network.pathReactions:
            # Sort bimolecular configurations so that we always encounter them in the
            # same order
            # The actual order doesn't matter, as long as it is consistent
            rxn.reactants.sort()
            rxn.products.sort()
            # All reactant configurations not already defined as reactants or 
            # isomers are assumed to be product channels
            if len(rxn.reactants) == 1 and rxn.reactants[0] not in network.isomers and rxn.reactants not in network.products:
                network.products.append(rxn.reactants)
            elif len(rxn.reactants) > 1 and rxn.reactants not in network.reactants and rxn.reactants not in network.products:
                network.products.append(rxn.reactants)
            # All product configurations not already defined as reactants or 
            # isomers are assumed to be product channels
            if len(rxn.products) == 1 and rxn.products[0] not in network.isomers and rxn.products not in network.products:
                network.products.append(rxn.products)
            elif len(rxn.products) > 1 and rxn.products not in network.reactants and rxn.products not in network.products:
                network.products.append(rxn.products)
        
        # For each configuration with states data but not thermo data,
        # calculate the thermo data
        for isomer in network.isomers:
            if isomer.thermo is None and isomer.states is not None:
                generateThermoFromStates(isomer, measure.Tlist.value_si)
        for reactants in network.reactants:
            for spec in reactants:
                if spec.thermo is None and spec.states is not None:
                    generateThermoFromStates(spec, measure.Tlist.value_si)
        for products in network.products:
            for spec in products:
                if spec.thermo is None and spec.states is not None:
                    generateThermoFromStates(spec, measure.Tlist.value_si)

        # Check that we have the right data for each configuration
        # Use a string to store the errors so that the user can see all of
        # the mistakes at once
        errorString = ''
        
        for isomer in network.isomers:
            errorString0 = ''
            # All isomers must have states data
            if isomer.states is None:
                errorString0 += '* Required molecular degree of freedom data was not provided.\n'
            # All isomers must have collision parameters
            if isomer.transportData is None:
                errorString0 += '* Required Lennard-Jones parameters were not provided.\n'
            if isomer.molecularWeight == 0:
                errorString0 += '* Required molecular weight was not provided.\n'
            # Either the isomer or the bath gas (or both) must have collision parameters
            if isomer.collisionModel is None and not any([spec.collisionModel is not None for spec in network.bathGas]):
                errorString0 += '* Collisional energy transfer model not provided.\n'
            if errorString0 != '':
                errorString = 'For unimolecular isomer "{0}":\n{1}\n'.format(isomer, errorString0)

        for rxn in network.pathReactions:
            errorString0 = ''
            
            # If both the reactants and the products are product channels, then the path reaction is invalid
            if rxn.reactants in network.products and rxn.products in network.products:
                errorString0 += '* Both the reactants and the products are product channels.\n'
            
            # Reactions of the form A + B -> C + D are not allowed
            elif len(rxn.reactants) > 1 and len(rxn.products) > 1:
                errorString0 += '* Both the reactants and the products are bimolecular.\n'
                
            # The remaining errors are only meaningful if the path reaction is allowed
            else:
                
                # All reactions must have either high-pressure-limit kinetics or transition state data
                if rxn.kinetics is None and rxn.transitionState.states is None:
                    errorString0 += '* Unable to determine microcanonical rate k(E). you must specify either\n  the high-P kinetics or transition state molecular degrees of freedom.\n'

                # Certain reactions require that both reactants and products have thermo data
                thermoRequired = False
                if rxn.reactants not in network.products and rxn.products not in network.products:
                    # Both reactants and products are isomers or reactant channels
                    # The reaction must have thermo for both reactants and products
                    thermoRequired = True
                elif rxn.reactants in network.products and rxn.kinetics is not None and rxn.transitionState.states is None:
                    # The reactants are product channels, and we will be using the ILT method to get k(E)
                    # We will need the thermodynamics to generate the reverse k(T) and k(E)
                    thermoRequired = True
                if thermoRequired:
                    for spec in rxn.reactants:
                        if spec.thermo is None:
                            errorString0 += '* Unable to determine thermodynamics data for reactant "{0}".\n  You must specify molecular degrees of freedom or a thermodynamics model.\n'.format(spec)
                    for spec in rxn.products:
                        if spec.thermo is None:
                            errorString0 += '* Unable to determine thermodynamics data for product "{0}".\n  You must specify molecular degrees of freedom or a thermodynamics model.\n'.format(spec)
                
            if errorString0 != '':
                errorString += 'For path reaction "{0}":\n{1}\n'.format(rxn, errorString0)

        if errorString != '':
            raise InputError(errorString)

    except InputError, e:
        network.errorString = str(e)
        logging.error('The input file "{0}" was invalid:'.format(path))
        logging.info(e)
    
    # Print lots of information about the loaded network
    # In particular, we want to give all of the energies on the PES
    # This will help the user decide if the range of energies selected is 
    # appropriate when viewing the log file
    network.printSummary(level=logging.INFO)
    
    # If there are no isomers, then there's nothing to do
    if len(network.isomers) == 0:
        message = 'Could not find any unimolecular isomers based on this network, so there is nothing to do.'
        network.errorString = message
        logging.info(message)
