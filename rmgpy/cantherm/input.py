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

import os.path
import logging
import math
import numpy
import matplotlib
matplotlib.rc('mathtext', default='regular')
            
from rmgpy.quantity import constants
from rmgpy.statmech import HinderedRotor, HarmonicOscillator
from rmgpy.species import Species, TransitionState
from rmgpy.kinetics import Arrhenius
from rmgpy.reaction import Reaction

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
# A dictionary associated species and transition state identifiers with geometry objects
geometryDict = {}

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
    pivots = [p-1 for p in pivots]
    top = [t-1 for t in top]
    return [scanLog, pivots, top, symmetry]

################################################################################

def loadConfiguration(energyLog, geomLog, statesLog, extSymmetry, spinMultiplicity, freqScaleFactor, linear, rotors, atoms, bonds, E0=None, TS=False):
    
    logging.debug('    Reading optimized geometry...')
    log = GaussianLog(geomLog)
    geom = log.loadGeometry()
    
    logging.debug('    Reading energy...')
    if E0 is None:
        if energyLog is not None: log = GaussianLog(energyLog)
        E0 = log.loadEnergy()
    else:
        E0 *= 4.35974394e-18 * constants.Na     # Hartree/particle to J/mol
    E0 = applyEnergyCorrections(E0, modelChemistry, atoms, bonds)
    logging.debug('         E0 (0 K) = %g kcal/mol' % (E0 / 4184))
    
    logging.debug('    Reading molecular degrees of freedom...')
    log = GaussianLog(statesLog)
    states = log.loadStates(symmetry=extSymmetry)
    states.spinMultiplicity = spinMultiplicity
    
    F = log.loadForceConstantMatrix()
    
    if F is not None and len(geom.mass) > 1 and len(rotors) > 0:
        
        logging.debug('    Fitting %i hindered rotors...' % len(rotors))
        for scanLog, pivots, top, symmetry in rotors:
            log = GaussianLog(scanLog)
            
            Vlist, angle = log.loadScanEnergies()
            
            inertia = geom.getInternalReducedMomentOfInertia(pivots, top)
            
            barr, symm = log.fitCosinePotential()
            cosineRotor = HinderedRotor(inertia=(inertia*constants.Na*1e23,"amu*angstrom^2"), symmetry=symm, barrier=(barr/4184.,"kcal/mol"))
            fourier = log.fitFourierSeriesPotential()
            fourierRotor = HinderedRotor(inertia=(inertia*constants.Na*1e23,"amu*angstrom^2"), symmetry=symmetry, fourier=(fourier,"J/mol"))
                
            Vlist_cosine = cosineRotor.getPotential(angle)
            Vlist_fourier = fourierRotor.getPotential(angle)
            
            rms_cosine = numpy.sqrt(numpy.sum((Vlist_cosine - Vlist) * (Vlist_cosine - Vlist)) / (len(Vlist) - 1)) / 4184.
            rms_fourier = numpy.sqrt(numpy.sum((Vlist_fourier - Vlist) * (Vlist_fourier - Vlist))/ (len(Vlist) - 1)) / 4184.
            print rms_cosine, rms_fourier, symm, symmetry
            
            # Keep the rotor with the most accurate potential
            rotor = cosineRotor if rms_cosine < rms_fourier else fourierRotor
            # However, keep the cosine rotor if it is accurate enough, the
            # fourier rotor is not significantly more accurate, and the cosine
            # rotor has the correct symmetry 
            if rms_cosine < 0.05 and rms_cosine / rms_fourier > 0.25 and rms_cosine / rms_fourier < 4.0 and symmetry == symm:
                rotor = cosineRotor
            
            states.modes.append(rotor)
            
            import pylab
            phi = numpy.arange(0, 6.3, 0.02, numpy.float64)
            fig = pylab.figure()
            pylab.plot(angle, Vlist / 4184, 'ok')
            linespec = '-r' if rotor is cosineRotor else '--r'
            pylab.plot(phi, cosineRotor.getPotential(phi) / 4184, linespec)
            linespec = '-b' if rotor is fourierRotor else '--b'
            pylab.plot(phi, fourierRotor.getPotential(phi) / 4184, linespec)
            pylab.legend(['scan', 'cosine', 'fourier'], loc=1)
            pylab.xlim(0, 2*math.pi)
            
            axes = fig.get_axes()[0]
            axes.set_xticks([float(j*math.pi/4) for j in range(0,9)])
            axes.set_xticks([float(j*math.pi/8) for j in range(0,17)], minor=True)
            axes.set_xticklabels(['$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$', '$5\pi/4$', '$3\pi/2$', '$7\pi/4$', '$2\pi$'])

            
        pylab.show()
        
        logging.debug('    Determining frequencies from reduced force constant matrix...')
        frequencies = list(projectRotors(geom, F, rotors, linear, TS))
        
    elif len(states.modes) > 2:
        frequencies = states.modes[2].frequencies.values
        rotors = []
    else:
        frequencies = []
        rotors = []

    for mode in states.modes:
        if isinstance(mode, HarmonicOscillator):
            mode.frequencies.values = numpy.array(frequencies, numpy.float) * freqScaleFactor

    return E0, geom, states

def loadSpecies(label, geomLog, statesLog, extSymmetry, spinMultiplicity, freqScaleFactor, linear, rotors, atoms, bonds, directory=None, E0=None, energyLog=None):
    global modelChemistry
    logging.info('Loading species %s...' % label)
    if directory:
        geomLog = os.path.join(directory, geomLog)
        statesLog = os.path.join(directory, statesLog)
        if energyLog: energyLog = os.path.join(directory, energyLog)
        for rotor in rotors:
            rotor[0] = os.path.join(directory, rotor[0])
    E0, geom, states = loadConfiguration(energyLog, geomLog, statesLog, extSymmetry, spinMultiplicity, freqScaleFactor, linear, rotors, atoms, bonds, E0, TS=False)
    speciesDict[label] = Species(label=label, thermo=None, states=states, E0=(E0/1000.,"kJ/mol"))
    geometryDict[label] = geom

def loadTransitionState(label, geomLog, statesLog, extSymmetry, spinMultiplicity, freqScaleFactor, linear, rotors, atoms, bonds, directory=None, E0=None, energyLog=None):
    global modelChemistry
    logging.info('Loading transition state %s...' % label)
    if directory:
        geomLog = os.path.join(directory, geomLog)
        statesLog = os.path.join(directory, statesLog)
        if energyLog: energyLog = os.path.join(directory, energyLog)
        for rotor in rotors:
            rotor[0] = os.path.join(directory, rotor[0])
    E0, geom, states = loadConfiguration(energyLog, geomLog, statesLog, extSymmetry, spinMultiplicity, freqScaleFactor, linear, rotors, atoms, bonds, E0, TS=True)
    log = GaussianLog(statesLog)
    frequency = log.loadNegativeFrequency()
    transitionStateDict[label] = TransitionState(label=label, states=states, frequency=(frequency,"cm^-1"), E0=(E0/1000.,"kJ/mol"))
    geometryDict[label] = geom
    
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
        saveStates(speciesDict[label], geometryDict[label], label, outputFile)
    elif label in transitionStateDict:
        saveStates(transitionStateDict[label], geometryDict[label], label, outputFile)
    
def generateThermo(label, model, plot=False):
    global outputFile, speciesDict
    from thermo import generateThermoModel, saveThermo
    generateThermoModel(speciesDict[label], model, plot)
    saveThermo(speciesDict[label], label, outputFile)
    
def generateKinetics(label, tunneling='', plot=False):
    global outputFile, reactionDict
    from kinetics import generateKineticsModel, saveKinetics
    generateKineticsModel(reactionDict[label], tunneling, plot)
    saveKinetics(reactionDict[label], tunneling, label, outputFile)
