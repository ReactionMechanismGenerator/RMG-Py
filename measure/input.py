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
    global Tlist, Pparams
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
    For a given :class:`Species` object `species`, use its states model to
    generate a corresponding thermo model.
    """
    # Do nothing if the species already has thermo data or if it does not have
    # states data
    if species.thermo is not None or species.states is None: return

    # Must use ThermoGAModel because we can't rely on knowing
    # anything about the structure of the species
    from chempy.thermo import ThermoGAModel
    Tdata = numpy.linspace(numpy.min(Tlist), numpy.max(Tlist), 20.0)
    Cpdata = species.states.getHeatCapacities(Tdata)
    H298 = species.states.getEnthalpy(298)
    S298 = species.states.getEntropy(298)
    species.thermo = ThermoGAModel(Tdata=Tdata, Cpdata=Cpdata, H298=H298, S298=S298)

def readInput(path):

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
        'bathGas': bathGas,
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
    else:
        Tmin = min(Tlist); Tmax = max(Tlist); Tcount = len(Tlist)
    
    # Determine pressure grid if not yet known
    if Pparams is not None and Plist is None:
        Pmin, Pmax, Pcount = Pparams
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
    else:
        Pmin = min(Plist); Pmax = max(Plist); Pcount = len(Plist)
    
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
        return None, None, None, None
      
    return network, Tlist, Plist, Elist, method, model, Tmin, Tmax, Pmin, Pmax

################################################################################

def writeInput(path, network, Tlist, Plist, Elist, method):
    """
    Write a MEASURE input file to `path` on disk.
    """

    def writeStatesData(f, states, prefix=''):
        f.write(prefix + 'states=States(\n')
        for mode in states.modes:
            if isinstance(mode, RigidRotor):
                f.write(prefix + '    rotationalConstants=([%s], "amu/angstrom^2"),\n' % (', '.join(['%g' % (inertia*constants.Na*1e23) for inertia in mode.inertia])))
                f.write(prefix + '    symmetry=%g,\n' % (mode.symmetry))
            elif isinstance(mode, HarmonicOscillator):
                f.write(prefix + '    frequencies=([%s], "cm^-1"),\n' % (', '.join([('%g' % freq) for freq in mode.frequencies])))
                f.write(prefix + '    frequencyScaleFactor=1.0,\n')
        if any([isinstance(mode, HinderedRotor) for mode in states.modes]):
            f.write(prefix + '    hinderedRotors=[\n')
            for mode in spec.states.modes:
                if isinstance(mode, HinderedRotor):
                    f.write(prefix + '        ((%g,"amu*angstrom^2"), (%g,"kJ/mol"), %i),\n' % (mode.inertia*constants.Na*1e23, mode.barrier/1000.0, mode.symmetry))
            f.write(prefix + '    ],\n')
        f.write(prefix + '),\n')

    # Get list of all species in network
    speciesList = []
    for isomer in network.isomers:
        if isomer not in speciesList: speciesList.append(isomer)
    for reactants in network.reactants:
        for spec in reactants:
            if spec not in speciesList: speciesList.append(spec)
    for products in network.products:
        for spec in products:
            if spec not in speciesList: speciesList.append(spec)
    if network.bathGas not in speciesList: speciesList.append(network.bathGas)

    f = open(path, 'w')

    # Write each species to file
    for spec in speciesList:
        f.write('species(\n')
        f.write('    label="%s",\n' % (spec.label))
        if len(spec.molecule) > 0:
            f.write('    SMILES="%s",\n' % (spec.molecule[0].toSMILES()))
        f.write('    E0=(%g,"kJ/mol"),\n' % (spec.E0 / 1000))
        if spec.states is not None:
            writeStatesData(f, spec.states, prefix='    ')
        if spec.thermo is not None:
            if isinstance(spec.thermo, ThermoGAModel):
                f.write('    thermo=ThermoGAModel(\n')
                f.write('        Tdata=([%s], "K"),\n' % (', '.join([('%g' % T) for T in spec.thermo.Tdata])))
                f.write('        Cpdata=([%s], "J/(mol*K)"),\n' % (', '.join([('%g' % Cp) for Cp in spec.thermo.Cpdata])))
                f.write('        H298=(%g, "kJ/mol"),\n' % (spec.thermo.H298/1000.0))
                f.write('        S298=(%g, "J/(mol*K)"),\n' % (spec.thermo.S298))
                f.write('    ),\n')
            else:
                f.write('    thermo=%r,\n' % spec.thermo)
        if spec.lennardJones is not None:
            f.write('    lennardJones=LennardJones(sigma=(%g,"m"), epsilon=(%g,"J")),\n' % (spec.lennardJones.sigma, spec.lennardJones.epsilon))
        if spec.molecularWeight != 0.0:
            f.write('    molecularWeight=(%g,"g/mol"),\n' % (spec.molecularWeight*1000))
        f.write(')\n\n')

    f.write('################################################################################\n\n')

    # Write each path reaction to file
    for rxn in network.pathReactions:
        f.write('reaction(\n')
        f.write('    reactants=[%s],\n' % (', '.join([('"%s"' % spec.label) for spec in rxn.reactants])))
        f.write('    products=[%s],\n' % (', '.join([('"%s"' % spec.label) for spec in rxn.products])))
        f.write('    reversible=%s,\n' % ('True' if rxn.reversible else 'False'))
        if rxn.kinetics is not None:
            if isinstance(rxn.kinetics, ArrheniusModel):
                f.write('    kinetics=Arrhenius(\n')
                if len(rxn.reactants) == 1: units = 's^-1'
                else: units = 'm^%g/(mol^%g*s)' % (3*(len(rxn.reactants)-1), len(rxn.reactants)-1)
                f.write('        A=(%g,"%s"),\n' % (rxn.kinetics.A, units))
                f.write('        n=%g,\n' % (rxn.kinetics.n))
                f.write('        Ea=(%g,"kJ/mol"),\n' % (rxn.kinetics.Ea / 1000.0))
                f.write('        T0=(%g,"K"),\n' % (rxn.kinetics.T0))
                f.write('    ),\n')
            else:
                f.write('    kinetics=%r,\n' % rxn.kinetics)
        if rxn.transitionState is not None:
            f.write('    transitionState=TransitionState(\n')
            f.write('        E0=(%g,"kJ/mol"),\n' % (rxn.transitionState.E0/1000))
            if rxn.transitionState.states is not None:
                writeStatesData(f, rxn.transitionState.states, prefix='        ')
            f.write('    ),\n')
        f.write(')\n\n')

    f.write('################################################################################\n\n')

    f.write('collisionModel(\n')
    if isinstance(network.collisionModel, SingleExponentialDownModel):
        f.write('    type="single exponential down",\n')
        f.write('    parameters=[\n')
        f.write('        (%g,"kJ/mol"),\n' % (network.collisionModel.alpha/1000.0))
        f.write('    ],\n')
    f.write(')\n\n')

    f.write('bathGas("%s")\n\n' % (network.bathGas.label))

    f.write('temperatures(([%s],"K"))\n' % (', '.join([('%g' % T) for T in Tlist])))
    f.write('pressures(([%s],"bar"))\n' % (', '.join([('%g' % (P/1e5)) for P in Plist])))
    dE, count = Elist
    f.write('energies(dE=(%g,"kJ/mol"), count=%i)\n' % (dE/1000.0, count))
    f.write('method("%s")\n' % (method))

    f.close()
