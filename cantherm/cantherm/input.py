#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   CanTherm - 
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

import chempy.constants as constants
from chempy.states import HinderedRotor, HarmonicOscillator
from chempy.species import Species, TransitionState
from chempy.kinetics import ArrheniusModel
from chempy.reaction import Reaction

from gaussian import GaussianLog
from states import projectRotors, applyEnergyCorrections

################################################################################

# The model chemistry used
# The energies for each species and transition state will be automatically 
# adjusted using standard reference energies for that method
modelChemistry = ''

# A dictionary associating species identifiers with species objects
speciesDict = {}
# A dictionary associating transition state identifiers with transition state objects
transitionStateDict = {}
# A dictionary associating reaction identifiers with reaction objects
reactionDict = {}

# The file to save the output to
outputFile = ''

################################################################################

def setOutputFile(path):
    global outputFile
    outputFile = path
    f = open(path, 'w')
    f.close()
    
def setModelChemistry(method):
    """
    Set the model chemistry used in this quantum chemisty calculation to
    `method`.
    """
    global modelChemistry
    modelChemistry = method

################################################################################

def hinderedRotor(scanLog, pivots, top, symmetry):
    return [scanLog, pivots, top, symmetry]

################################################################################

def loadConfiguration(geomLog, statesLog, extSymmetry, freqScaleFactor, linear, rotors, atoms, bonds, E0=None, TS=False):
    
    logging.debug('    Reading optimized geometry...')
    log = GaussianLog(geomLog)
    geom = log.loadGeometry()
    
    logging.debug('    Reading energy...')
    if E0 is None:
        E0 = log.loadEnergy()
    else:
        E0 *= 4.35974394e-18 * constants.Na     # Hartree/particle to J/mol
    E0 = applyEnergyCorrections(E0, modelChemistry, atoms, bonds)
    logging.debug('         E0 (0 K) = %g kcal/mol' % (E0 / 4184))
    
    logging.debug('    Reading molecular degrees of freedom...')
    log = GaussianLog(statesLog)
    states = log.loadStates(symmetry=extSymmetry)

    F = log.loadForceConstantMatrix()
    
    if F is not None and len(geom.mass) > 1 and len(rotors) > 0:
        
        logging.debug('    Fitting %i hindered rotors...' % len(rotors))
        for scanLog, pivots, top, symmetry in rotors:
            log = GaussianLog(scanLog)
            fourier = log.fitFourierSeriesPotential()
            inertia = geom.getInternalReducedMomentOfInertia(pivots, top)
            rotor = HinderedRotor(inertia=inertia, symmetry=symmetry, fourier=fourier)
            states.modes.append(rotor)
            
            #import numpy
            #import pylab
            #import math
            #Vlist = log.loadScanEnergies()
            #Vlist = Vlist[:-1]
            #angle = numpy.arange(0.0, 2*math.pi+0.00001, 2*math.pi/(len(Vlist)-1), numpy.float64)
            #phi = numpy.arange(0, 6.3, 0.02, numpy.float64)
            #pylab.plot(angle, Vlist / 4184, 'ok')
            #pylab.plot(phi, rotor.getPotential(phi) / 4184, '-k')
        #pylab.show()
        
        logging.debug('    Determining frequencies from reduced force constant matrix...')
        frequencies = list(projectRotors(geom, F, rotors, linear, TS))
        
    elif len(states.modes) > 2:
        frequencies = states.modes[2].frequencies
        rotors = []
    else:
        frequencies = []
        rotors = []

    for mode in states.modes:
        if isinstance(mode, HarmonicOscillator):
            mode.frequencies = [f * freqScaleFactor for f in frequencies]

    return E0, geom, states

def loadSpecies(label, geomLog, statesLog, extSymmetry, freqScaleFactor, linear, rotors, atoms, bonds, E0=None):
    global modelChemistry
    logging.info('Loading species %s...' % label)
    E0, geom, states = loadConfiguration(geomLog, statesLog, extSymmetry, freqScaleFactor, linear, rotors, atoms, bonds, E0, TS=False)
    speciesDict[label] = Species(label=label, thermo=None, states=states, geometry=geom, E0=E0)

def loadTransitionState(label, geomLog, statesLog, extSymmetry, freqScaleFactor, linear, rotors, atoms, bonds, E0=None):
    global modelChemistry
    logging.info('Loading transition state %s...' % label)
    E0, geom, states = loadConfiguration(geomLog, statesLog, extSymmetry, freqScaleFactor, linear, rotors, atoms, bonds, E0, TS=True)
    log = GaussianLog(statesLog)
    frequency = log.loadNegativeFrequency()
    transitionStateDict[label] = TransitionState(label=label, states=states, geometry=geom, frequency=frequency, E0=E0)
    
################################################################################

def loadReaction(label, reactants, products, transitionState, degeneracy=1):
    global speciesDict, transitionStateDict, reactionDict
    logging.info('Loading reaction %s...' % label)
    rxn = Reaction(
        reactants=[speciesDict[s] for s in reactants],
        products=[speciesDict[s] for s in products],
        transitionState=transitionStateDict[transitionState],
    )
    rxn.degeneracy = degeneracy
    reactionDict[label] = rxn

################################################################################

def generateStates(label):
    global outputFile, speciesDict, transitionStateDict
    from states import saveStates
    if label in speciesDict:
        saveStates(speciesDict[label], label, outputFile)
    elif label in transitionStateDict:
        saveStates(transitionStateDict[label], label, outputFile)
    
def generateThermo(label, model, plot=False):
    global outputFile, speciesDict
    from thermo import generateThermoModel, saveThermo
    generateThermoModel(speciesDict[label], model, plot)
    saveThermo(speciesDict[label], label, outputFile)
    
def generateKinetics(label, tunneling='', plot=False):
    global outputFile, reactionDict
    from kinetics import generateKineticsModel, saveKinetics
    generateKineticsModel(reactionDict[label], tunneling, plot)
    saveKinetics(reactionDict[label], label, outputFile)
