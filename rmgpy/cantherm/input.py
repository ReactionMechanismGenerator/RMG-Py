#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains functionality for parsing CanTherm input files.
"""

import os.path
import logging

from rmgpy.species import Species, TransitionState

from rmgpy.statmech.translation import Translation, IdealGasTranslation
from rmgpy.statmech.rotation import Rotation, LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor
from rmgpy.statmech.vibration import Vibration, HarmonicOscillator
from rmgpy.statmech.torsion import Torsion, HinderedRotor, FreeRotor
from rmgpy.statmech.conformer import Conformer

from rmgpy.thermo.thermodata import ThermoData
from rmgpy.thermo.nasa import NASAPolynomial, NASA
from rmgpy.thermo.wilhoit import Wilhoit

from rmgpy.kinetics.arrhenius import Arrhenius, ArrheniusEP, PDepArrhenius, MultiArrhenius, MultiPDepArrhenius 
from rmgpy.kinetics.chebyshev import Chebyshev
from rmgpy.kinetics.falloff import ThirdBody, Lindemann, Troe
from rmgpy.kinetics.kineticsdata import KineticsData, PDepKineticsData
from rmgpy.kinetics.tunneling import Wigner, Eckart

from rmgpy.pdep.configuration import Configuration
from rmgpy.pdep.network import Network  
from rmgpy.pdep.collision import SingleExponentialDown

from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.transport import TransportData

from rmgpy.cantherm.kinetics import KineticsJob
from rmgpy.cantherm.statmech import StatMechJob
from rmgpy.cantherm.thermo import ThermoJob
from rmgpy.cantherm.pdep import PressureDependenceJob

################################################################################

speciesDict = {}
transitionStateDict = {}
reactionDict = {}
networkDict = {}
jobList = []

################################################################################

def species(label, *args, **kwargs):
    global speciesDict, jobList
    if label in speciesDict:
        raise ValueError('Multiple occurrences of species with label {0!r}.'.format(label))
    logging.info('Loading species {0}...'.format(label))
    
    spec = Species(label=label)
    speciesDict[label] = spec
    
    if len(args) == 1:
        # The argument is a path to a conformer input file
        path = args[0]
        job = StatMechJob(species=spec, path=path)
        jobList.append(job)
    
    if len(kwargs) > 0:
        # The species parameters are given explicitly
        structure = None
        E0 = None
        modes = []
        spinMultiplicity = 1
        opticalIsomers = 1
        molecularWeight = None
        collisionModel = None
        energyTransferModel = None
        thermo = None
        for key, value in kwargs.items():
            if key == 'structure':
                structure = value
            elif key == 'E0':
                E0 = value
            elif key == 'modes':
                modes = value
            elif key == 'spinMultiplicity':
                spinMultiplicity = value
            elif key == 'opticalIsomers':
                opticalIsomers = value
            elif key == 'molecularWeight':
                molecularWeight = value
            elif key == 'collisionModel':
                collisionModel = value
            elif key == 'energyTransferModel':
                energyTransferModel = value
            elif key == 'thermo':
                thermo = value
            else:
                raise TypeError('species() got an unexpected keyword argument {0!r}.'.format(key))
            
        if structure: spec.molecule = [structure]
        spec.conformer = Conformer(E0=E0, modes=modes, spinMultiplicity=spinMultiplicity, opticalIsomers=opticalIsomers)  
        spec.molecularWeight = molecularWeight
        spec.transportData = collisionModel
        spec.energyTransferModel = energyTransferModel
        spec.thermo = thermo
        
    return spec

def transitionState(label, *args, **kwargs):
    global transitionStateDict
    if label in transitionStateDict:
        raise ValueError('Multiple occurrences of transition state with label {0!r}.'.format(label))
    logging.info('Loading transition state {0}...'.format(label))
    ts = TransitionState(label=label)
    transitionStateDict[label] = ts
    
    if len(args) == 1 and len(kwargs) == 0:
        # The argument is a path to a conformer input file
        path = args[0]
        job = StatMechJob(species=ts, path=path)
        jobList.append(job)
    
    elif len(args) == 0 and len(kwargs) > 0:
        # The species parameters are given explicitly
        E0 = None
        modes = []
        spinMultiplicity = 1
        opticalIsomers = 1
        frequency = None
        for key, value in kwargs.items():
            if key == 'E0':
                E0 = value
            elif key == 'modes':
                modes = value
            elif key == 'spinMultiplicity':
                spinMultiplicity = value
            elif key == 'opticalIsomers':
                opticalIsomers = value
            elif key == 'frequency':
                frequency = value
            else:
                raise TypeError('transitionState() got an unexpected keyword argument {0!r}.'.format(key))
        
        ts.conformer = Conformer(E0=E0, modes=modes, spinMultiplicity=spinMultiplicity, opticalIsomers=opticalIsomers)  
        ts.frequency = frequency
        
    return ts

def reaction(label, reactants, products, transitionState, kinetics=None, tunneling=''):
    global reactionDict, speciesDict, transitionStateDict
    #label = 'reaction'+transitionState
    if label in reactionDict:
        label = label+transitionState
        if label in reactionDict:
            raise ValueError('Multiple occurrences of reaction with label {0!r}.'.format(label))
    logging.info('Loading reaction {0}...'.format(label))
    reactants = sorted([speciesDict[spec] for spec in reactants])
    products = sorted([speciesDict[spec] for spec in products])
    transitionState = transitionStateDict[transitionState]
    if tunneling.lower() == 'wigner':
        transitionState.tunneling = Wigner(frequency=None)
    elif tunneling.lower() == 'eckart':
        transitionState.tunneling = Eckart(frequency=None, E0_reac=None, E0_TS=None, E0_prod=None)
    elif tunneling == '' or tunneling is None:
        transitionState.tunneling = None
    elif not isinstance(tunneling, TunnelingModel):
        raise ValueError('Unknown tunneling model {0!r}.'.format(tunneling))
    rxn = Reaction(label=label, reactants=reactants, products=products, transitionState=transitionState, kinetics=kinetics)
    reactionDict[label] = rxn
    
    return rxn

def network(label, isomers=None, reactants=None, products=None, pathReactions=None, bathGas=None):
    global networkDict, speciesDict, reactionDict
    logging.info('Loading network {0}...'.format(label))
    isomers0 = isomers or []; isomers = []
    for isomer in isomers0:
        if isinstance(isomer, (list,tuple)):
            raise ValueError('Only one species can be present in a unimolecular isomer.')
        isomers.append(speciesDict[isomer])
    
    reactants0 = reactants or []; reactants = []
    for reactant in reactants0:
        if not isinstance(reactant, (list,tuple)):
            reactant = [reactant]
        reactants.append(sorted([speciesDict[spec] for spec in reactant]))
    
    if pathReactions is None:
        # If not explicitly given, use all reactions in input file
        pathReactions = reactionDict.values()
    else:
        pathReactions0 = pathReactions; pathReactions = []
        for rxn in pathReactions0:
            pathReactions.append(reactionDict[rxn])
    
    if products is None:
        # Figure out which configurations are isomers, reactant channels, and product channels
        products = []
        for rxn in pathReactions:
            # Sort bimolecular configurations so that we always encounter them in the
            # same order
            # The actual order doesn't matter, as long as it is consistent
            rxn.reactants.sort()
            rxn.products.sort()
            # All reactant configurations not already defined as reactants or 
            # isomers are assumed to be product channels
            if len(rxn.reactants) == 1 and rxn.reactants[0] not in isomers and rxn.reactants not in products:
                products.append(rxn.reactants)
            elif len(rxn.reactants) > 1 and rxn.reactants not in reactants and rxn.reactants not in products:
                products.append(rxn.reactants)
            # All product configurations not already defined as reactants or 
            # isomers are assumed to be product channels
            if len(rxn.products) == 1 and rxn.products[0] not in isomers and rxn.products not in products:
                products.append(rxn.products)
            elif len(rxn.products) > 1 and rxn.products not in reactants and rxn.products not in products:
                products.append(rxn.products)
    else:
        products0 = products or []; products = []
        for product in products0:
            if not isinstance(product, (list,tuple)):
                product = [product]
            products.append(sorted([speciesDict[spec] for spec in product]))

    isomers = [Configuration(species) for species in isomers]
    reactants = [Configuration(*species) for species in reactants]
    products = [Configuration(*species) for species in products]
        
    bathGas0 = bathGas or {}; bathGas = {}
    for spec, fraction in bathGas0.items():
        bathGas[speciesDict[spec]] = fraction
    
    network = Network(
        label = label, 
        isomers = isomers, 
        reactants = reactants, 
        products = products, 
        pathReactions = pathReactions, 
        bathGas = bathGas,
    )
    networkDict[label] = network

def kinetics(label,Tmin=None, Tmax=None,Tlist=None, Tcount=0):
    global jobList, reactionDict
    try:
        rxn = reactionDict[label]
    except KeyError:
        raise ValueError('Unknown reaction label {0!r} for kinetics() job.'.format(label))
    job = KineticsJob(reaction=rxn,Tmin=Tmin, Tmax=Tmax,Tcount = Tcount,Tlist=Tlist)
    jobList.append(job)

def statmech(label):
    global jobList, speciesDict, transitionStateDict
    if label in speciesDict or label in transitionStateDict:
        for job in jobList:
            if job.species.label == label:
                break
        else:
            raise ValueError('Could not create StatMechJob for {0!r}; no path specified.'.format(label))
    else:
        raise ValueError('Unknown species or transition state label {0!r} for statmech() job.'.format(label))

def thermo(label, thermoClass):
    global jobList, speciesDict
    try:
        spec = speciesDict[label]
    except KeyError:
        raise ValueError('Unknown species label {0!r} for thermo() job.'.format(label))
    job = ThermoJob(species=spec, thermoClass=thermoClass)
    jobList.append(job)

def pressureDependence(label, 
                       Tmin=None, Tmax=None, Tcount=0, Tlist=None,
                       Pmin=None, Pmax=None, Pcount=0, Plist=None,
                       maximumGrainSize=None, minimumGrainCount=0,
                       method=None, interpolationModel=None,
                       activeKRotor=True, activeJRotor=True, rmgmode=False):
    global jobList, networkDict
    if isinstance(interpolationModel, str):
        interpolationModel = (interpolationModel,)
    job = PressureDependenceJob(network = networkDict[label],
        Tmin=Tmin, Tmax=Tmax, Tcount=Tcount, Tlist=Tlist,
        Pmin=Pmin, Pmax=Pmax, Pcount=Pcount, Plist=Plist,
        maximumGrainSize=maximumGrainSize, minimumGrainCount=minimumGrainCount,
        method=method, interpolationModel=interpolationModel,
        activeKRotor=activeKRotor, activeJRotor=activeJRotor,
        rmgmode=rmgmode,
    )
    jobList.append(job)

def SMILES(smiles):
    return Molecule().fromSMILES(smiles)

def InChI(inchi):
    return Molecule().fromInChI(inchi)

################################################################################

def loadInputFile(path):
    """
    Load the CanTherm input file located at `path` on disk, and return a list of
    the jobs defined in that file.
    """
    global speciesDict, transitionStateDict, reactionDict, networkDict, jobList
    
    # Clear module-level variables
    speciesDict = {}
    transitionStateDict = {}
    reactionDict = {}
    networkDict = {}
    jobList = []
    
    global_context = { '__builtins__': None }
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'range': range,
        # Collision
        'TransportData': TransportData,
        'SingleExponentialDown': SingleExponentialDown,
        # Kinetics
        'Arrhenius': Arrhenius,
        # Statistical mechanics
        'IdealGasTranslation': IdealGasTranslation,
        'LinearRotor': LinearRotor,
        'NonlinearRotor': NonlinearRotor,
        'KRotor': KRotor,
        'SphericalTopRotor': SphericalTopRotor,
        'HarmonicOscillator': HarmonicOscillator,
        'HinderedRotor': HinderedRotor,
        'FreeRotor':FreeRotor,
        # Thermo
        'ThermoData': ThermoData,
        'Wilhoit': Wilhoit,
        'NASA': NASA,
        'NASAPolynomial': NASAPolynomial,
        # Functions
        'reaction': reaction,
        'species': species,
        'transitionState': transitionState,
        'network': network,
        # Jobs
        'kinetics': kinetics,
        'statmech': statmech,
        'thermo': thermo,
        'pressureDependence': pressureDependence,
        # Miscellaneous
        'SMILES': SMILES,
        'InChI': InChI,
    }

    with open(path, 'r') as f:
        try:
            exec f in global_context, local_context
        except (NameError, TypeError, SyntaxError), e:
            logging.error('The input file {0!r} was invalid:'.format(path))
            raise

    modelChemistry = local_context.get('modelChemistry', '')
    if 'frequencyScaleFactor' not in local_context:
        logging.warning('No frequency scale factor specified in input file; assuming a value of unity.')
    frequencyScaleFactor = local_context.get('frequencyScaleFactor', 1.0)
    useHinderedRotors = local_context.get('useHinderedRotors', True)
    useBondCorrections = local_context.get('useBondCorrections', False)
    
    directory = os.path.dirname(path)
    
    for job in jobList:
        if isinstance(job, StatMechJob):
            job.path = os.path.join(directory, job.path)
            job.modelChemistry = modelChemistry
            job.frequencyScaleFactor = frequencyScaleFactor
            job.includeHinderedRotors = useHinderedRotors
            job.applyBondEnergyCorrections = useBondCorrections
    
    return jobList
