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

################################################################################

def setModelChemistry(method):
    """
    Set the model chemistry used in this quantum chemisty calculation to
    `method`.
    """
    global modelChemistry
    modelChemistry = method

################################################################################

def loadConfiguration(geomLog, statesLog, extSymmetry, freqScaleFactor, linear, rotorPivots, rotorTops, rotorScans, rotorSymmetry, atoms, bonds):
    
    logging.debug('    Reading optimized geometry...')
    log = GaussianLog(geomLog)
    geom = log.loadGeometry()
    
    logging.debug('    Reading energy...')
    E0 = applyEnergyCorrections(log.loadEnergy(), modelChemistry, atoms, bonds)
    logging.debug('         E0 (0 K) = %g kcal/mol' % (E0 / 4184))
    
    logging.debug('    Reading molecular degrees of freedom...')
    log = GaussianLog(statesLog)
    states = log.loadStates(symmetry=extSymmetry)

    F = log.loadForceConstantMatrix()
    if F is not None and len(geom.mass) > 1 and len(rotorScans) > 0:
        logging.debug('    Fitting %i hindered rotors...' % len(rotorScans))
        rotors = loadHinderedRotors(geom, rotorPivots, rotorTops, rotorScans, rotorSymmetry)

        if len(rotors) > 0:
            logging.debug('    Determining frequencies from reduced force constant matrix...')
            frequencies = projectRotors(geom, F, rotors, rotorPivots, rotorTops, linear)
            for mode in states.modes:
                if isinstance(mode, HarmonicOscillator):
                    mode.frequencies = list(frequencies * freqScaleFactor)
            states.modes.extend(rotors)
        elif len(states.modes) > 2:
            frequencies = states.modes[2].frequencies
            rotors = []
        else:
            frequencies = []
            rotors = []
    else:
        frequencies = []
        rotors = []
    
    return E0, geom, states

def loadSpecies(label, geomLog, statesLog, extSymmetry, freqScaleFactor, linear, rotorPivots, rotorTops, rotorScans, rotorSymmetry, atoms, bonds):
    global modelChemistry
    logging.info('Loading species %s...' % label)
    E0, geom, states = loadConfiguration(geomLog, statesLog, extSymmetry, freqScaleFactor, linear, rotorPivots, rotorTops, rotorScans, rotorSymmetry, atoms, bonds)
    speciesDict[label] = Species(label=label, thermo=None, states=states, geometry=geom, E0=E0)

def loadTransitionState(label, geomLog, statesLog, extSymmetry, freqScaleFactor, linear, rotorPivots, rotorTops, rotorScans, rotorSymmetry, atoms, bonds):
    global modelChemistry
    logging.info('Loading transition state %s...' % label)
    E0, geom, states = loadConfiguration(geomLog, statesLog, extSymmetry, freqScaleFactor, linear, rotorPivots, rotorTops, rotorScans, rotorSymmetry, atoms, bonds)
    frequency = log.loadNegativeFrequency()
    transitionStateDict[label] = TransitionState(label=label, states=states, geometry=geom, frequency=frequency)
    
################################################################################

def loadReaction(label, reactants, products, transitionState):
    global speciesDict, transitionStateDict, reactionDict
    logging.info('Loading reaction %s...' % label)
    rxn = Reaction(
        reactants=[speciesDict[s] for s in reactants],
        products=[speciesDict[s] for s in products],
        transitionState=transitionStateDict[transitionState],
    )
    reactionDict[label] = rxn

################################################################################

def loadHinderedRotors(geom, pivots, top1, scans, symmetry):
    rotors = []
    for i in range(len(scans)):
        log = GaussianLog(scans[i])
        fourier = log.fitFourierSeriesPotential()
        inertia = geom.getInternalReducedMomentOfInertia(pivots[i], top1[i])
        rotor = HinderedRotor(inertia=inertia, symmetry=symmetry[i], fourier=fourier)
        rotors.append(rotor)
    return rotors

################################################################################

def generateThermo(label, plot=False):
    global speciesDict
    from thermo import generateThermoModel
    generateThermoModel(speciesDict[label], plot)

def generateKinetics(label, plot=False):
    global reactionDict
    from kinetics import generateKineticsModel
    generateKineticsModel(reactionDict[label], plot)
