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

from chempy.species import Species, TransitionState
from chempy.reaction import Reaction
from chempy.species import LennardJones as LennardJonesModel
from chempy.molecule import Molecule
from chempy.states import *
from chempy.kinetics import ArrheniusModel
from chempy.thermo import *

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

def processQuantity(quantity):
    """
    Processes a `quantity` from the input file. The quantity can be a number of
    things:

    * A :data:`list` or :data:`tuple` containing a number and a string with its
      units (e.g ``[100,'kJ/mol']``)

    * A :data:`list` or :data:`tuple` containing a list or tuple of numbers and
      a string with their units (e.g ``[(1,2,3), 'm']``)

    * A :data:`list` or :data:`tuple` of numbers with no units (e.g.
      ``[1, 2, 3]``)

    * A single number with no units

    If the quantity has units, then the associated number(s) are converted to
    SI units using the ``quantities`` package.  If the quantity does not have
    units, then it is assumed to be either a dimensionless quantity or a
    quantity that is already in SI units.

    This function returns a :data:`tuple` containing the number or list of
    numbers and a string with the units of those numbers.
    """
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

def species(label='', E0=None, states=None, thermo=None, lennardJones=None, molecularWeight=0.0, SMILES='', InChI=''):
    global speciesDict
    if E0 is not None: E0 = processQuantity(E0)[0]
    else: E0 = 0.0
    spec = Species(label=label, states=states, thermo=thermo, E0=E0, lennardJones=lennardJones)
    if InChI != '':
        spec.molecule = [Molecule(InChI=InChI)]
    elif SMILES != '':
        spec.molecule = [Molecule(SMILES=SMILES)]
    spec.molecularWeight = processQuantity(molecularWeight)[0]
    speciesDict[label] = spec
    logging.debug('Found species "%s"' % spec)

def thermoGAModel(Tdata, Cpdata, H298, S298, Tmin=0.0, Tmax=99999.9, comment=''):
    return ThermoGAModel(
        Tdata=processQuantity(Tdata)[0],
        Cpdata=processQuantity(Cpdata)[0],
        H298=processQuantity(H298)[0],
        S298=processQuantity(S298)[0],
        Tmin=processQuantity(Tmin)[0],
        Tmax=processQuantity(Tmax)[0],
        comment=comment,
    )

def thermoWilhoitModel(cp0, cpInf, a0, a1, a2, a3, H0, S0, B, comment=''):
    return WilhoitModel(
        cp0=processQuantity(cp0)[0],
        cpInf=processQuantity(cpInf)[0],
        a0=a0,
        a1=a1,
        a2=a2,
        a3=a3,
        H0=H0,
        S0=S0,
        B=processQuantity(B)[0],
        comment=comment,
    )

def thermoNASAModel(polynomials=None, Tmin=0.0, Tmax=0.0, comment=''):
    return NASAModel(
        polynomials=polynomials,
        Tmin=processQuantity(Tmin)[0],
        Tmax=processQuantity(Tmax)[0],
        comment=comment,
    )

def thermoNASAPolynomial(Tmin=0.0, Tmax=0.0, coeffs=None, comment=''):
    return NASAPolynomial(
        Tmin=processQuantity(Tmin)[0],
        Tmax=processQuantity(Tmax)[0],
        coeffs=coeffs,
        comment=comment,
    )


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

def Arrhenius(A, n, Ea, T0=1.0):
    A, units = processQuantity(A)
    n, units = processQuantity(n)
    Ea, units = processQuantity(Ea)
    T0, units = processQuantity(T0)
    return ArrheniusModel(A=A, n=n, Ea=Ea, T0=T0)

def TS(E0=None, states=None, frequency=0.0):
    if E0 is not None: E0 = processQuantity(E0)[0]
    frequency = processQuantity(frequency)[0]
    return TransitionState(E0=E0, states=states, frequency=frequency)

def collisionModel(type, parameters, bathGas):
    global network, speciesDict
    if type.lower() == 'single exponential down':

        # Process parameters, making sure we have a valid set
        if len(parameters) == 1:
            if 'alpha' not in parameters:
                raise InputError('Must specify either "alpha" or ("alpha0","T0","n") as parameters for SingleExponentialDownModel.')
            alpha0 = processQuantity(parameters['alpha'])[0]
            T0 = 1000.0
            n = 0.0
        elif len(parameters) == 3:
            if 'alpha0' not in parameters or 'T0' not in parameters or 'n' not in parameters:
                raise InputError('Must specify either "alpha" or ("alpha0","T0","n") as parameters for SingleExponentialDownModel.')
            alpha0 = processQuantity(parameters['alpha0'])[0]
            T0 = processQuantity(parameters['T0'])[0]
            n = processQuantity(parameters['n'])[0]
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
        Tlist = processQuantity(Tlist0)[0]
        Tparams = None
    elif Tmin is not None and Tmax is not None and count is not None:
        # We've been provided a temperature range and number of temperatures to use
        # We defer choosing the actual temperatures because they depend on the
        # choice of interpolation model
        Tlist = None
        Tparams = [processQuantity(Tmin)[0], processQuantity(Tmax)[0], count]
    else:
        raise SyntaxError('Must specify either a list of temperatures or Tmin, Tmax, and count.')

def pressures(Plist0=None, Pmin=None, Pmax=None, count=None):
    global Plist, Pparams
    if Plist0 is not None:
        # We've been provided a list of specific pressures to use
        Plist = processQuantity(Plist0)[0]
        Pparams = None
    elif Pmin is not None and Pmax is not None and count is not None:
        # We've been provided a pressures range and number of pressures to use
        # We defer choosing the actual pressures because they depend on the
        # choice of interpolation model
        Plist = None
        Pparams = [processQuantity(Pmin)[0], processQuantity(Pmax)[0], count]
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

def interpolationModel(name, *args):
    global model
    model = [name]
    model.extend(list(args))

################################################################################

def generateThermoFromStates(species):
    """
    For a given :class:`Species` object `species` with molecular degrees of
    freedom data in its ``states`` attribute, generate a corresponding thermo
    model. A :class:`ThermoGAModel` thermodynamics model is stored in the
    species ``thermo`` attribute, and nothing is returned.
    """
    # Do nothing if the species already has thermo data or if it does not have
    # states data
    if species.thermo is not None or species.states is None: return

    # Must use ThermoGAModel because we can't rely on knowing
    # anything about the structure of the species
    from chempy.thermo import ThermoGAModel
    Tdata = numpy.linspace(numpy.min(Tlist), numpy.max(Tlist), 20.0)
    Cpdata = species.states.getHeatCapacities(Tdata)
    H298 = species.E0 + species.states.getEnthalpy(298)
    S298 = species.states.getEntropy(298)
    species.thermo = ThermoGAModel(Tdata=Tdata, Cpdata=Cpdata, H298=H298, S298=S298)

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
    format is described in the :ref:`usersguide`. Returns a number of
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
        'temperatures': temperatures,
        'pressures': pressures,
        'energies': energies,
        'method': _method,
        'interpolationModel': interpolationModel,
        'ThermoGAModel': thermoGAModel,
        'WilhoitModel': thermoWilhoitModel,
        'NASAModel': thermoNASAModel,
        'NASAPolynomial': thermoNASAPolynomial,
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
    if network is None: return None, None, None, None, None, None, None, None, None, None
    
    # Determine temperature grid if not yet known
    if Tparams is not None and Tlist is None:
        Tmin, Tmax, Tcount = Tparams
        Tlist = getTemperaturesForModel(model, Tmin, Tmax, Tcount)
    else:
        Tmin = min(Tlist); Tmax = max(Tlist); Tcount = len(Tlist)
    Tlist = numpy.array(Tlist, numpy.float64)

    # Determine pressure grid if not yet known
    if Pparams is not None and Plist is None:
        Pmin, Pmax, Pcount = Pparams
        Plist = getPressuresForModel(model, Pmin, Pmax, Pcount)
    else:
        Pmin = min(Plist); Pmax = max(Plist); Pcount = len(Plist)
    Plist = numpy.array(Plist, numpy.float64)

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

    # Print lots of information about the loaded network
    # In particular, we want to give all of the energies on the PES
    # This will help the user decide if the range of energies selected is 
    # appropriate when viewing the log file
    network.printSummary(level=logging.INFO)
    
    # If there are no isomers, then there's nothing to do
    if len(network.isomers) == 0:
        logging.info('Could not find any unimolecular isomers based on this network, so there is nothing to do.')
        return None, None, None, None, None, None, None, None, None, None
      
    return network, Tlist, Plist, Elist, method, model, Tmin, Tmax, Pmin, Pmax

