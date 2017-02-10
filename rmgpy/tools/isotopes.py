#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
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
This module contains functionality for generating mechanisms with isotopes.
"""

import os
import os.path
import logging
import numpy as np
import itertools
from copy import copy
import pandas as pd
import shutil
import math

import rmgpy.constants as constants
from rmgpy.molecule import Molecule
from rmgpy.molecule.element import getElement
from rmgpy.tools.loader import loadRMGJob
from rmgpy.chemkin import ChemkinWriter
from rmgpy.rmg.input import getInput
from rmgpy.rmg.main import RMG, initializeLog
from rmgpy.rmg.model import Species
from rmgpy.rmg.listener import SimulationProfileWriter
from rmgpy.thermo.thermoengine import processThermoData
from rmgpy.data.thermo import findCp0andCpInf

def initializeIsotopeModel(rmg, isotopes):
    """
    Initialize the RMG object by using the parameter species list
    as initial species instead of the species from the RMG input file.

    """
    # Read input file
    rmg.loadInput(rmg.inputFile)

    # Check input file 
    rmg.checkInput()

    # Load databases
    rmg.loadDatabase()

    for spc in isotopes:
        spec, isNew = rmg.reactionModel.makeNewSpecies(spc)
        spec.thermo = spc.thermo        
        if isNew:
            rmg.reactionModel.addSpeciesToEdge(spec)
            rmg.initialSpecies.append(spec)

    for spec in rmg.initialSpecies:
        spec.thermo = processThermoData(spec, spec.thermo)
        if not spec.reactive:
            rmg.reactionModel.enlarge(spec)
    for spec in rmg.initialSpecies:
        if spec.reactive:
            rmg.reactionModel.enlarge(spec)

    rmg.initializeReactionThresholdAndReactFlags()
    rmg.reactionModel.initializeIndexSpeciesDict()


def generateIsotopeModel(outputDirectory, rmg0, isotopes):
    """
    Replace the core species of the rmg model with the parameter list
    of species.

    Generate all reactions between new list of core species.

    Returns created RMG object.
    """
    rmg = RMG(inputFile=rmg0.inputFile, outputDirectory=outputDirectory)
    rmg.attach(ChemkinWriter(outputDirectory))

    initializeIsotopeModel(rmg, isotopes)

    rmg.reactionModel.enlarge(reactEdge=True,
        unimolecularReact=rmg.unimolecularReact,
        bimolecularReact=rmg.bimolecularReact)

    rmg.saveEverything()

    rmg.finish()     

    return rmg   

def generateIsotopomers(spc, N=1):
    """
    Generate all isotopomers of the parameter species by adding max. N carbon isotopes to the
    atoms of the species.
    """

    mol = spc.molecule[0]
    isotope = getElement(6, 13)
    carbons = filter(lambda at: at.symbol == isotope.symbol, mol.atoms)
    
    mols = []
    addIsotope(0, N, mol, mols, isotope)

    spcs = []
    for isomol in mols:
        isotopomer = Species(molecule=[isomol], thermo=spc.thermo, transportData=spc.transportData, reactive=spc.reactive)
        isotopomer.generateResonanceIsomers()
        spcs.append(isotopomer)

    # do not retain identical species:
    filtered = []
    while spcs:
        candidate = spcs.pop()
        unique = True
        for isotopomer in filtered:
            if isotopomer.isIsomorphic(candidate):
                unique = False
                break
        if unique: filtered.append(candidate)

    for isotopomer in filtered:
        correctEntropy(isotopomer, spc)

    return filtered

def addIsotope(i, N, mol, mols, element):
    """
    
    Iterate over the atoms of the molecule, and changes the element object
    of the atom by the provide parameter element object. Add the newly created
    isotopomer to the list of Molecule objects. For each created isotopomer,
    recursively call the method, until the maximum number of isotopes per molecule
    (N) is reached.

    """
    if i == N: return
    else:
        atoms = filter(lambda at: at.symbol == element.symbol, mol.atoms)
        for at in atoms:
            if at.element == element: continue
            else:
                isotopomer = mol.copy(deep=True)
                isotopomer.atoms[mol.atoms.index(at)].element = element
                mols.append(isotopomer)
                addIsotope(i+1, N, isotopomer, mols, element)


def solve(rmg):
    """
    Solve the reaction system, read the simulation 
    profiles, and return them into a pandas dataframe.
    """

    solverDir = os.path.join(rmg.outputDirectory, 'solver')
    try:
        shutil.rmtree(solverDir)
    except OSError, e:
        pass
    
    os.mkdir(solverDir)

    reactionSysIndex = 0
    listener = SimulationProfileWriter(rmg.outputDirectory, reactionSysIndex, rmg.reactionModel.core.species)
    
    reactionSystem = rmg.reactionSystems[0]
    reactionSystem.attach(listener)

    reactionModel = rmg.reactionModel

    # run simulation:
    terminated, obj = reactionSystem.simulate(
        coreSpecies = reactionModel.core.species,
        coreReactions = reactionModel.core.reactions,
        edgeSpecies = reactionModel.edge.species,
        edgeReactions = reactionModel.edge.reactions,
        toleranceKeepInEdge = 0,
        toleranceMoveToCore = 1,
        toleranceInterruptSimulation = 1,
    ) 

    simCSV = os.path.join(rmg.outputDirectory, 'solver/simulation_{}_{}.csv'.format(reactionSysIndex + 1, len(reactionModel.core.species)))
    spcdata = pd.read_csv(simCSV)
    
    return spcdata

def cluster(spcList):
    """
    Creates subcollections of isotopomers that belong together.
    """

    unclustered = copy(spcList)

    # [[list of Species objs]]
    clusters = []

    while unclustered:
        candidate = unclustered.pop()
        for cluster in clusters:
            if any([removeIsotope(spc).isIsomorphic(removeIsotope(candidate)) for spc in cluster]):
                cluster.append(candidate)
                break
        else:
            clusters.append([candidate])

    return clusters

def removeIsotope(spc):
    """
    Create a deep copy of the first molecule of the species object and replace
    non-normal Element objects (of special isotopes) by the 
    expected isotope.
    """
    stripped = spc.copy(deep=True)

    for atom in stripped.molecule[0].atoms:
        if atom.element.isotope != -1:
            atom.element = getElement(atom.element.symbol)

    # only do it for the first molecule, generate the other resonance isomers.
    stripped.molecule = [stripped.molecule[0]]
    stripped.generateResonanceIsomers()

    return stripped

def retrieveConcentrations(spcdata, clusters):
    """
    Iterate over the species in the list of clustered species
    and return a dataframe, but filled with 
    concentration columns of the corresponding species.
    """

    concs = []

    for cluster in clusters:
        df = pd.DataFrame()
        for spc in cluster:
            try:
                header = '{}({})'.format(spc.label, spc.index)
                df[header] = spcdata[header]
            except KeyError, e:
                header = '{}'.format(spc.label)
                try:
                    df[header] = spcdata[header]
                except KeyError, e:
                    raise e
            
        concs.append(df)

    return concs

def computeProbabilities(df):
    """
    Compute the isotope probabilities by dividing the 
    species concentration by the sum of the clustered species concentrations.
    """
    probs = []
    sumConcs = df.sum(axis=1)
    for column in df:
        df[column] = df[column] / sumConcs

    return df

def generateRMGModel(inputFile, outputDirectory):
    """
    Generate the RMG-Py model NOT containing any non-normal isotopomers.

    Returns created RMG object.
    """
    initializeLog(logging.INFO, os.path.join(outputDirectory, 'RMG.log'))
    # generate mechanism:
    rmg = RMG(inputFile = os.path.abspath(inputFile),
            outputDirectory = os.path.abspath(outputDirectory)
        )
    rmg.execute()

    return rmg

def correctEntropy(isotopomer, isotopeless):
    """
    Correct the entropy of the isotopomer by the following correction for symmetry:

    S(corrected) = S(original) + R*ln(sigma(isotopeless)) - R*ln(sigma(isotopomer))
    """

    # calculate -R ln (sigma) in SI units (J/K/mol)
    Sisotopeless = - constants.R * math.log(isotopeless.getSymmetryNumber())
    Sisotopomer = - constants.R * math.log(isotopomer.getSymmetryNumber())

    # convert species thermo to ThermoData object:
    nasa = isotopomer.thermo
    thermo = nasa.toThermoData()

    # apply correction to entropy at 298K
    thermo.S298.value_si -= Sisotopeless
    thermo.S298.value_si += Sisotopomer

    # put the corrected thermo back as a species attribute:
    isotopomer.thermo = thermo

def run(inputFile, isotopeInputFile, outputDir, original=None, isotopeLoc=None):
    """
    Accepts two input files, one input file with the RMG-Py model to generate, NOT
    containing any non-normal isotopomers, and one input file for the model to be 
    generated starting from the isotopomers generated from the core species of the first model.

    Firstly, generates the RMG model for the first input file. Takes the core species of that mechanism
    and generates all isotopomers of those core species. Next, generates all reactions between the
    generated pool of isotopomers, and writes it to file. Next, reads in that newly generated mechanism
    and runs a simulation with the given conditions of the second input file, writes species mole fractions
    to file. Next, clusters the isotopomers together again, so that isotopomer probabilities can be computed.

    Returns:
    - a pandas data frame with isotopomer probabilities as a function of reaction time 
    - a list of species objects corresponding to the isotopomers of the generated model
    - a pandas data frame with isotopomer speciation data as a function of reaction time
    """

    if not original:
        outputdirRMG = os.path.join(outputDir, 'rmg')
        os.mkdir(outputdirRMG)

        rmg = generateRMGModel(inputFile, outputdirRMG)
    else:
        outputdirRMG = original
        chemkinFile = os.path.join(outputdirRMG, 'chemkin', 'chem_annotated.inp')
        dictFile = os.path.join(outputdirRMG, 'chemkin', 'species_dictionary.txt')
        rmg = loadRMGJob(inputFile, chemkinFile, dictFile, generateImages=False, useChemkinNames=True)

    isotopeInputFile = os.path.abspath(isotopeInputFile)
    rmgIso = RMG(inputFile = isotopeInputFile)
    rmgIso.loadInput(rmgIso.inputFile)

    if not isotopeLoc:
        print('Generating isotopomers for the core species in {}'.format(outputdirRMG))
        isotopes = []

        try:
            speciesConstraints = getInput('speciesConstraints')
        except Exception, e:
            logging.debug('Species constraints could not be found.')
            raise e

        try:
            maxIsotopes = speciesConstraints['maximumIsotopicAtoms']
        except KeyError, e:
            print ('Could not find the maxIsotopicAtoms constraint in the input file. Exiting...')
            raise e

        for spc in rmg.reactionModel.core.species:
            findCp0andCpInf(spc, spc.thermo)
            isotopes.extend(generateIsotopomers(spc, maxIsotopes))

        # add the original unlabeled species:
        isotopes.extend(rmg.reactionModel.core.species)
        print('Number of isotopomers: {}'.format(len(isotopes)))

        outputdirIso = os.path.join(outputDir, 'iso')
        os.mkdir(outputdirIso)

        print('Generating RMG isotope model in {}'.format(outputdirIso))
        rmgIso = generateIsotopeModel(outputdirIso, rmg, isotopes)
    else:
        outputdirIso= isotopeLoc

    chemkinFileIso = os.path.join(outputdirIso, 'chemkin', 'chem_annotated.inp')
    dictFileIso = os.path.join(outputdirIso, 'chemkin', 'species_dictionary.txt')

    print('Loading isotope chemkin model.\nInput file: {}\nChemkin file: {}\nDict file: {}'\
        .format(isotopeInputFile, chemkinFileIso, dictFileIso))

    rmgIso = loadRMGJob(isotopeInputFile, chemkinFileIso, dictFileIso, generateImages=True, useChemkinNames=True)
    rmgIso.outputDirectory = outputdirIso

    print('Clustering isotopomers...')
    clusters = cluster(rmgIso.reactionModel.core.species)

    print('Solving the Isotope model with the isotope input file...')
    spcData = solve(rmgIso)    
    
    print('Generating concentrations for lumped species...')
    concs = retrieveConcentrations(spcData, clusters)

    print('Computing isotopomer probabilities...')
    probs = []
    for df in concs:
        df = computeProbabilities(df.copy())
        probs.append(df)

    return probs, rmgIso.reactionModel.core.species, spcData
