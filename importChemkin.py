#!/usr/bin/env python
# encoding: utf-8

"""
This script enables the conversion of Chemkin files
into RMG-Py style thermo library file.
Simply pass the paths of the Chemkin files on the 
command-line, e.g.

    $ python convertThermo.py --species /path/to/chem1.inp  --thermo /path/to/therm.dat

If you supply a --species file (containing a SPECIES block) this is used to limit
the species converted.
The resulting file is saved next to the thermo input file.

If running in QTconsole, it draws pictures of the species.
"""

import os.path
import argparse
import logging
import re

import cherrypy
import json
import threading
import urllib

import rmgpy
import rmgpy.rmg
from rmgpy.display import display

from rmgpy.chemkin import loadChemkinFile, readSpeciesBlock, readThermoBlock, readReactionsBlock, removeCommentFromLine
from rmgpy.reaction import ReactionModel

from rmgpy.data.thermo import Entry, saveEntry
from rmgpy.molecule import Molecule
from rmgpy.rmg.model import Species  # you need this one, not the one in rmgpy.species!

from rmgpy.rmg.main import RMG, initializeLog
from rmgpy.molecule.draw import MoleculeDrawer

import time
import sys
# Put the RMG-database project at the start of the python path, so we use that importOldDatabase script!
databaseDirectory = rmgpy.settings['database.directory']
databaseProjectDirectory = os.path.abspath(os.path.join(databaseDirectory, '..'))
sys.path.insert(0, databaseProjectDirectory)
from importOldDatabase import getUsername
################################################################################

class MagicSpeciesDict(dict):
    """
    A dictionary that always has the species you're looking for!
    
    (If they key isn't there when you ask, it adds a blank Species() object and returns that)
    """
    def __init__(self, dictionary=None):
        self.dictionary = dictionary or {}
    def __getitem__(self, key):
        if key not in self.dictionary:
            self.dictionary[key] = Species()
        return dict.__getitem__(self.dictionary, key)

def convertFormula(formulaDict):
    """
    Given a formula in dict form {'c':2, 'h':6, 'o':0}
    return a canonical formula string "C2H6"
    
    For comparison reasons, this must be the same algorithm as used in
    rmgpy.molecule.Molecule class.
    """

    elements = {e.capitalize(): n for e, n in formulaDict.iteritems() if n > 0}
    hasCarbon = 'C' in elements
    hasHydrogen = 'H' in elements
    # Use the Hill system to generate the formula
    formula = ''
    # Carbon and hydrogen always come first if carbon is present
    if hasCarbon:
        count = elements['C']
        formula += 'C{0:d}'.format(count) if count > 1 else 'C'
        del elements['C']
        if hasHydrogen:
            count = elements['H']
            formula += 'H{0:d}'.format(count) if count > 1 else 'H'
            del elements['H']
    # Other atoms are in alphabetical order
    # (This includes hydrogen if carbon is not present)
    keys = elements.keys()
    keys.sort()
    for key in keys:
        count = elements[key]
        formula += '{0}{1:d}'.format(key, count) if count > 1 else key
    return formula

def parseCommandLineArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--species', metavar='FILE', type=str, nargs='?', default=None,
        help='the Chemkin file containing the list of species')
    parser.add_argument('--reactions', metavar='FILE', type=str, nargs='?', default=None,
        help='the Chemkin file containing the list of reactions')
    parser.add_argument('--thermo', metavar='FILE', type=str,
        help='the Chemkin files containing the thermo')
    parser.add_argument('-o', '--output-directory', type=str, nargs=1, default='',
        metavar='DIR', help='use DIR as output directory')
    parser.add_argument('-s', '--scratch-directory', type=str, nargs=1, default='',
        metavar='DIR', help='use DIR as scratch directory')

    # Options for controlling the amount of information printed to the console
    # By default a moderate level of information is printed; you can either
    # ask for less (quiet), more (verbose), or much more (debug)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-q', '--quiet', action='store_true', help='only print warnings and errors')
    group.add_argument('-v', '--verbose', action='store_true', help='print more verbose output')
    group.add_argument('-d', '--debug', action='store_true', help='print debug information')

    args = parser.parse_args()

    # add some args that RMG will want

    args.walltime = '0'
    args.restart = False

    inputDirectory = os.path.abspath(os.path.dirname(args.thermo))
    if args.output_directory == '':
        args.output_directory = os.path.join(inputDirectory, 'RMG-Py-output')
    if args.scratch_directory == '':
        args.scratch_directory = args.output_directory

    return args

def reactionMatchesFormulas(reaction, reactantFormulas):
    """
    Does the reaction match the reactantFormulas, in either direction?
    
    The reaction must contain real species with identified molecules.
    reactantFormulas is an iterable of Strings like ('CH3','H').
    Returns 'forward', 'backward', or False.
    """
    reactantFormulas = sorted(reactantFormulas)
    if reactantFormulas == sorted([s.getFormula() for s in reaction.reactants]):
        return 'forward'
    elif reactantFormulas == sorted([s.getFormula() for s in reaction.products]):
        return 'backward'
    else:
        return False



class ModelMatcher():
    """
    For identifying species in an imported model
    """
    def __init__(self, args=None):
        self.args = args
        self.speciesDict = {}
        self.thermoDict = {}
        self.formulaDict = {}
        self.smilesDict = {}
        self.identified_labels = []
        self.identified_unprocessed_labels = []
        self.speciesList = None
        self.speciesDict_rmg = {}
        self.chemkinReactions = []
        self.chemkinReactionsUnmatched = []
        self.suggestedMatches = {}
        self.votes = {}
        self.prunedVotes = {}
        self.manualMatchesToProcess = []

    def loadModel(self, species_file, reactions_file, thermo_file):
        print 'Loading model...'

        speciesAliases = {}
        speciesDict = {}
        if species_file:
            logging.info("Reading species list...")
            speciesList = []
            with open(species_file) as f:
                line0 = f.readline()
                while line0 != '':
                    line = removeCommentFromLine(line0)[0]
                    tokens_upper = line.upper().split()
                    if tokens_upper and tokens_upper[0] in ('SPECIES', 'SPEC'):
                        # Unread the line (we'll re-read it in readReactionBlock())
                        f.seek(-len(line0), 1)
                        readSpeciesBlock(f, speciesDict, speciesAliases, speciesList)
                    line0 = f.readline()
        else:
            logging.info("No species file to limit species. Will read everything in thermo file")
            speciesList = None
            speciesDict = MagicSpeciesDict(speciesDict)

        logging.info("Reading thermo...")
        with open(thermo_file) as f:
            line0 = f.readline()
            while line0 != '':
                line = removeCommentFromLine(line0)[0]
                tokens_upper = line.upper().split()
                if tokens_upper and tokens_upper[0] in ('THERMO', 'THER'):
                    # Unread the line (we'll re-read it in readThermoBlock())
                    f.seek(-len(line0), 1)
                    formulaDict = readThermoBlock(f, speciesDict)
                line0 = f.readline()

        # Save the formulaDict, converting from {'c':1,'h':4} into "CH4" in the process.
        self.formulaDict = {label: convertFormula(formula) for label, formula in formulaDict.iteritems()}

        # thermoDict contains original thermo as read from chemkin thermo file
        self.thermoDict = {s.label: s.thermo for s in speciesDict.values() }

        self.speciesList = speciesList
        self.speciesDict = speciesDict

    def loadReactions(self, reactions_file):
        logging.info("Reading reactions...")
        with open(reactions_file) as f:
            reactionList = readReactionsBlock(f, self.speciesDict, readComments=True)
        logging.info("Read {0} reactions from chemkin file.".format(len(reactionList)))
        self.chemkinReactions = reactionList
        self.chemkinReactionsUnmatched = self.chemkinReactions[:]  # make a copy


    def initializeRMG(self, args):
        """
        Creata an RMG object, store it in self.rmg_object, and set it up.
        
        This loads the database, makes some settings, etc.
        `args` should have attributes `output_directory` and `scratch_directory`.
        Also needs a global variable databseDirectory
        """
        logging.info("Loading RMG database...")
        rmg = RMG()
        rmg.outputDirectory = args.output_directory
        rmg.scratchDirectory = args.scratch_directory
        rmg.makeOutputSubdirectory('species')
        rmg.databaseDirectory = databaseDirectory
        rmg.thermoLibraries = ['primaryThermoLibrary', 'KlippensteinH2O2', 'DFT_QCI_thermo', 'CBS_QB3_1dHR']
        rmg.kineticsFamilies = ['!Substitution_O']
        rmg.reactionLibraries = [('Glarborg/HighP', False)]
        rmg.loadDatabase()
        logging.info("Loaded database.")

        rmg.reactionModel = rmgpy.rmg.model.CoreEdgeReactionModel()
        rmg.reactionModel.kineticsEstimator = 'rate rules'
        rmg.initialSpecies = []
        rmg.reactionSystems = []

        self.rmg_object = rmg
        return rmg


    def speciesMatch(self, rmg_species, chemkin_species):
        """
        Return True if the species might match, else False.
        
        i.e. if chemkin_species has been identified, it must be the rmg_species,
        but if it hasn't it must at least have the same formula.
        If it matches based only on formula, the match it is added to the self.suggestedMatches dictionary.
        """
        chemkin_label = chemkin_species.label
        identified_labels = self.identified_labels
        if chemkin_label in identified_labels:
            return self.speciesDict_rmg[chemkin_label] is rmg_species
        elif rmg_species.label in identified_labels:
            return False
        else:
            if self.formulaDict[chemkin_label] == rmg_species.molecule[0].getFormula():
                self.suggestedMatches[chemkin_label] = rmg_species
                return True
            else:
                return False

    def reactionsMatch(self, rmg_reaction, chemkin_reaction, eitherDirection=True):
        """
        This is based on the rmg.reaction.Reaction.isIsomorphic method
 
        Return ``True`` if rmg_reaction is the same as the chemkin_reaction reaction,
        or ``False`` if they are different. 
        If `eitherDirection=False` then the directions must match.
        """
        speciesMatch = self.speciesMatch
        # Compare reactants to reactants

        # get things we refer to a lot into the local namespace, to reduce lookups
        rmg_reactants = rmg_reaction.reactants
        ck_reactants = chemkin_reaction.reactants
        len_rmg_reactants = len(rmg_reactants)
        len_ck_reactants = len(ck_reactants)

        forwardReactantsMatch = False
        if len_rmg_reactants == 1 and len_ck_reactants == 1:
            if speciesMatch(rmg_reactants[0], ck_reactants[0]):
                forwardReactantsMatch = True
        elif len_rmg_reactants == 2 and len_ck_reactants == 2:
            if speciesMatch(rmg_reactants[0], ck_reactants[0]) and speciesMatch(rmg_reactants[1], ck_reactants[1]):
                forwardReactantsMatch = True
            elif speciesMatch(rmg_reactants[0], ck_reactants[1]) and speciesMatch(rmg_reactants[1], ck_reactants[0]):
                forwardReactantsMatch = True
        elif len_rmg_reactants == 3 and len_ck_reactants == 3:
            if (speciesMatch(rmg_reactants[0], ck_reactants[0]) and
                    speciesMatch(rmg_reactants[1], ck_reactants[1]) and
                    speciesMatch(rmg_reactants[2], ck_reactants[2])):
                forwardReactantsMatch = True
            elif (speciesMatch(rmg_reactants[0], ck_reactants[0]) and
                    speciesMatch(rmg_reactants[1], ck_reactants[2]) and
                    speciesMatch(rmg_reactants[2], ck_reactants[1])):
                forwardReactantsMatch = True
            elif (speciesMatch(rmg_reactants[0], ck_reactants[1]) and
                    speciesMatch(rmg_reactants[1], ck_reactants[0]) and
                    speciesMatch(rmg_reactants[2], ck_reactants[2])):
                forwardReactantsMatch = True
            elif (speciesMatch(rmg_reactants[0], ck_reactants[2]) and
                    speciesMatch(rmg_reactants[1], ck_reactants[0]) and
                    speciesMatch(rmg_reactants[2], ck_reactants[1])):
                forwardReactantsMatch = True
            elif (speciesMatch(rmg_reactants[0], ck_reactants[1]) and
                    speciesMatch(rmg_reactants[1], ck_reactants[2]) and
                    speciesMatch(rmg_reactants[2], ck_reactants[0])):
                forwardReactantsMatch = True
            elif (speciesMatch(rmg_reactants[0], ck_reactants[2]) and
                    speciesMatch(rmg_reactants[1], ck_reactants[1]) and
                    speciesMatch(rmg_reactants[2], ck_reactants[0])):
                forwardReactantsMatch = True
        elif len_rmg_reactants == len_ck_reactants:
            raise NotImplementedError("Can't check isomorphism of reactions with {0} reactants".format(len_rmg_reactants))


        rmg_products = rmg_reaction.products
        ck_products = chemkin_reaction.products
        len_rmg_products = len(rmg_products)
        len_ck_products = len(ck_products)
        # Compare products to products
        forwardProductsMatch = False
        if len_rmg_products == 1 and len_ck_products == 1:
            if speciesMatch(rmg_products[0], ck_products[0]):
                forwardProductsMatch = True
        elif len_rmg_products == 2 and len_ck_products == 2:
            if speciesMatch(rmg_products[0], ck_products[0]) and speciesMatch(rmg_products[1], ck_products[1]):
                forwardProductsMatch = True
            elif speciesMatch(rmg_products[0], ck_products[1]) and speciesMatch(rmg_products[1], ck_products[0]):
                forwardProductsMatch = True
        elif len_rmg_products == 3 and len_ck_products == 3:
            if (speciesMatch(rmg_products[0], ck_products[0]) and
                    speciesMatch(rmg_products[1], ck_products[1]) and
                    speciesMatch(rmg_products[2], ck_products[2])):
                forwardProductsMatch = True
            elif (speciesMatch(rmg_products[0], ck_products[0]) and
                    speciesMatch(rmg_products[1], ck_products[2]) and
                    speciesMatch(rmg_products[2], ck_products[1])):
                forwardProductsMatch = True
            elif (speciesMatch(rmg_products[0], ck_products[1]) and
                    speciesMatch(rmg_products[1], ck_products[0]) and
                    speciesMatch(rmg_products[2], ck_products[2])):
                forwardProductsMatch = True
            elif (speciesMatch(rmg_products[0], ck_products[2]) and
                    speciesMatch(rmg_products[1], ck_products[0]) and
                    speciesMatch(rmg_products[2], ck_products[1])):
                forwardProductsMatch = True
            elif (speciesMatch(rmg_products[0], ck_products[1]) and
                    speciesMatch(rmg_products[1], ck_products[2]) and
                    speciesMatch(rmg_products[2], ck_products[0])):
                forwardProductsMatch = True
            elif (speciesMatch(rmg_products[0], ck_products[2]) and
                    speciesMatch(rmg_products[1], ck_products[1]) and
                    speciesMatch(rmg_products[2], ck_products[0])):
                forwardProductsMatch = True
        elif len_rmg_products == len_ck_products:
            raise NotImplementedError("Can't check isomorphism of reactions with {0} products".format(len_rmg_products))

        # Return now, if we can
        if (forwardReactantsMatch and forwardProductsMatch):
            return True
        if not eitherDirection:
            return False

        # Compare reactants to products
        reverseReactantsMatch = False
        if len_rmg_reactants == 1 and len_ck_products == 1:
            if speciesMatch(rmg_reactants[0], ck_products[0]):
                reverseReactantsMatch = True
        elif len_rmg_reactants == 2 and len_ck_products == 2:
            if speciesMatch(rmg_reactants[0], ck_products[0]) and speciesMatch(rmg_reactants[1], ck_products[1]):
                reverseReactantsMatch = True
            elif speciesMatch(rmg_reactants[0], ck_products[1]) and speciesMatch(rmg_reactants[1], ck_products[0]):
                reverseReactantsMatch = True
        elif len_rmg_reactants == 3 and len_ck_products == 3:
            if (speciesMatch(rmg_reactants[0], ck_products[0]) and
                    speciesMatch(rmg_reactants[1], ck_products[1]) and
                    speciesMatch(rmg_reactants[2], ck_products[2])):
                reverseReactantsMatch = True
            elif (speciesMatch(rmg_reactants[0], ck_products[0]) and
                    speciesMatch(rmg_reactants[1], ck_products[2]) and
                    speciesMatch(rmg_reactants[2], ck_products[1])):
                reverseReactantsMatch = True
            elif (speciesMatch(rmg_reactants[0], ck_products[1]) and
                    speciesMatch(rmg_reactants[1], ck_products[0]) and
                    speciesMatch(rmg_reactants[2], ck_products[2])):
                reverseReactantsMatch = True
            elif (speciesMatch(rmg_reactants[0], ck_products[2]) and
                    speciesMatch(rmg_reactants[1], ck_products[0]) and
                    speciesMatch(rmg_reactants[2], ck_products[1])):
                reverseReactantsMatch = True
            elif (speciesMatch(rmg_reactants[0], ck_products[1]) and
                    speciesMatch(rmg_reactants[1], ck_products[2]) and
                    speciesMatch(rmg_reactants[2], ck_products[0])):
                reverseReactantsMatch = True
            elif (speciesMatch(rmg_reactants[0], ck_products[2]) and
                    speciesMatch(rmg_reactants[1], ck_products[1]) and
                    speciesMatch(rmg_reactants[2], ck_products[0])):
                reverseReactantsMatch = True
        elif len_rmg_reactants == len_ck_products:
            raise NotImplementedError("Can't check isomorphism of reactions with {0} reactants".format(len_rmg_reactants))

        # Compare products to reactants
        reverseProductsMatch = False
        if len_rmg_products == 1 and len_ck_reactants == 1:
            if speciesMatch(rmg_products[0], ck_reactants[0]):
                reverseProductsMatch = True
        elif len_rmg_products == 2 and len_ck_reactants == 2:
            if speciesMatch(rmg_products[0], ck_reactants[0]) and speciesMatch(rmg_products[1], ck_reactants[1]):
                reverseProductsMatch = True
            elif speciesMatch(rmg_products[0], ck_reactants[1]) and speciesMatch(rmg_products[1], ck_reactants[0]):
                reverseProductsMatch = True
        elif len_rmg_products == 3 and len_ck_reactants == 3:
            if (speciesMatch(rmg_products[0], ck_reactants[0]) and
                    speciesMatch(rmg_products[1], ck_reactants[1]) and
                    speciesMatch(rmg_products[2], ck_reactants[2])):
                reverseProductsMatch = True
            elif (speciesMatch(rmg_products[0], ck_reactants[0]) and
                    speciesMatch(rmg_products[1], ck_reactants[2]) and
                    speciesMatch(rmg_products[2], ck_reactants[1])):
                reverseProductsMatch = True
            elif (speciesMatch(rmg_products[0], ck_reactants[1]) and
                    speciesMatch(rmg_products[1], ck_reactants[0]) and
                    speciesMatch(rmg_products[2], ck_reactants[2])):
                reverseProductsMatch = True
            elif (speciesMatch(rmg_products[0], ck_reactants[2]) and
                    speciesMatch(rmg_products[1], ck_reactants[0]) and
                    speciesMatch(rmg_products[2], ck_reactants[1])):
                reverseProductsMatch = True
            elif (speciesMatch(rmg_products[0], ck_reactants[1]) and
                    speciesMatch(rmg_products[1], ck_reactants[2]) and
                    speciesMatch(rmg_products[2], ck_reactants[0])):
                reverseProductsMatch = True
            elif (speciesMatch(rmg_products[0], ck_reactants[2]) and
                    speciesMatch(rmg_products[1], ck_reactants[1]) and
                    speciesMatch(rmg_products[2], ck_reactants[0])):
                reverseProductsMatch = True
        elif len_rmg_products == len_ck_reactants:
            raise NotImplementedError("Can't check isomorphism of reactions with {0} products".format(len_rmg_products))

        # should have already returned if it matches forwards, or we're not allowed to match backwards
        return  (reverseReactantsMatch and reverseProductsMatch)

    def identifySmallMolecules(self):
        """Identify anything little that is uniquely determined by its chemical formula"""

        known_formulas = {
             'CH4': 'C',
             'CH3': '[CH3]',
             'CO2': 'O=C=O',
             'H2O': 'O',
             'HO': '[OH]',
             'C2H6': 'CC',
             'C2H5': 'C[CH2]',
             'C2H4': 'C=C',
             'O2': '[O][O]',
             'H2': '[H][H]',
             'H2O2': 'OO',
             'O': '[O]',
             'N2': 'N#N',
             'CO': '[C]=O',
             'HO2': '[O]O',
             'C3H8': 'CCC',
             'CH': '[CH]',
             }
        identified_labels = []
        known_labels = {
                        'mb': 'CCCC(=O)OC',
                        'mb3d': 'C=CCC(=O)OC',
                        'c2h5coch3': 'CCC(=O)C',
                        'c3h4-a': 'C=C=C',
                        'mb4j': '[CH2]CCC(=O)OC',
                        'mb3j': 'C[CH]CC(=O)OC',
                        'mbmj': 'CCCC(=O)O[CH2]',
                        'c3h5-a': 'C(=C)[CH2]',
                        'c2h3cho': 'C(=C)C=O',
                        'ch2cho': '[CH2]C=O',
                        'c2h5co': 'CC[C]=O',
'ac3h4coch3': '[CH2]C=CC(=O)C',
'ac3h5co': 'C=CC[C]=O',
'baoj': 'CCCC(=O)[O]',
'baoj2*o': 'CCC(=O)C(=O)[O]',
'c2h3co': 'C=C[C]=O',
'c2h5chco': 'CCC=C=O',
'c2h5coch2': 'CCC(=O)[CH2]',
'c2h5coco': 'CCC(=O)[C]=O',
'c2h5o2': 'CCO[O]',
'c3h5o': '[O]CC=C',
'c3h6cho-1': '[CH2]CCC=O',
'c3h6cho-3': 'CC[CH]C=O',
'c3h6coc2h5-1': '[CH2]CCC(=O)CC',
'c3h6coc2h5-3': 'CC[CH]C(=O)CC',
'c3h6coch3-1': '[CH2]CCC(=O)C',
'c3h6coch3-2': 'C[CH]CC(=O)C',
'c3h6o1-3': 'C1CCO1',
'c3h6ooh2-1': 'OOC([CH2])C',
'c5h7o2': 'COC(=O)[CH]C=C',
'ch2cch2oh': 'OC[C]=C',
'ch2ch2cho': '[CH2]CC=O',
'ch2ch2coch3': '[CH2]CC(=O)C',
'ch2choohcoch3': '[CH2]C(C(=O)C)OO',
'ch3chch2co': 'C[CH]CC=O',
'ch3chchco': 'C[CH]C=C=O',
'ch3chcho': 'C[CH]C=O',
'ch3coco': 'O=[C]C(=O)C',
'ch3cooch2': '[CH2]OC(=O)C',
'ic3h5cho': 'CC(=C)C=O',
'ic3h5co': 'CC(=C)[C]=O',
'ic3h5coc2h4p': '[CH2]CC(=O)C(=C)C',
'ic3h5coch2': 'CC=CC(=O)[CH2]',
'ic3h5coch3': 'CC=CC(=O)C',
'ic3h6coc2h5': 'CCC(=O)C(C)[CH2]',
'ic3h7co': 'O=[C]C(C)C',
'ic3h7coc2h4s': 'C[CH]C(=O)C(C)C',
'ic3h7coch2': '[CH2]C(=O)C(C)C',
'ic3h7coch3': 'CC(=O)C(C)C',
'ic3h7o': 'CC([O])C',
'ic3h7o2': '[O]OC(C)C',
'mb2o': 'CCC(C(=O)OC)[O]',
'mb2oo': 'CCC(C(=O)OC)O[O]',
'mb2ooh3j': 'C[CH]C(C(=O)OC)OO',
'mb2ooh4j': '[CH2]CC(C(=O)OC)OO',
'mb2oohmj': 'CCC(C(=O)O[CH2])OO',
'mb3o': 'COC(=O)CC([O])C',
'mb3oo': '[O]OC(CC(=O)OC)C',
'mb3ooh2j': 'OOC([CH]C(=O)OC)C',
'mb3ooh4j': 'OOC(CC(=O)OC)[CH2]',
'mb4j*o': 'COC(=O)CC[C]=O',
'mb4o': 'COC(=O)CCC[O]',
'mb4ooh2j': 'OOCC[CH]C(=O)OC',
'mb4ooh3j': 'OOC[CH]CC(=O)OC',
'mb4ooh3oo': 'OOCC(CC(=O)OC)O[O]',
'mbmj*o': 'CCCC(=O)O[C]=O',
'mbmooh2j': 'CC[CH]C(=O)OCOO',
'me2*omj*o': 'O=[C]OC(=O)C=O',
'me2j*o': 'COC(=O)[C]=O',
'mp2d2j': 'COC(=O)[C]=C',
'mp2d3j': 'COC(=O)C=[CH]',
'mp2dmj': '[CH2]OC(=O)C=C',
'mp3j*o': 'COC(=O)C[C]=O',
'mp3j*o2*o': 'COC(=O)C(=O)[C]=O',
'nc3h7coch2': 'CCCC(=O)[CH2]',
'nc3h7o': 'CCC[O]',
'pc2h4coc2h3': '[CH2]CC(=O)C=C',
'sc2h4coc2h3': 'C[CH]C(=O)C=C',
'tc3h6cho': 'O=C[C](C)C',
'tc3h6coch3': 'CC(=O)[C](C)C',
        }
       # known_labels.clear()  # for debugging
        # use speciesList if it is not None or empty, else the formulaDict keys.
        for species_label in [s.label for s in self.speciesList or []] or self.formulaDict.keys():
            formula = self.formulaDict[species_label]
            if formula in known_formulas:
                known_smiles = known_formulas[formula]
                logging.info("I think {0} is {1} based on its formula".format(species_label, known_smiles))
                smiles = known_smiles
            elif species_label in known_labels:
                known_smiles = known_labels[species_label]
                logging.info("I think {0} is {1} based on its label".format(species_label, known_smiles))
                smiles = known_smiles
            else:
                continue
            self.smilesDict[species_label] = smiles
            while formula != Molecule(SMILES=smiles).getFormula():
                smiles = raw_input("SMILES {0} has formula {1} not required formula {2}. Try again:\n".format(smiles, Molecule(SMILES=smiles).getFormula(), formula))
            species = self.speciesDict[species_label]
            species.molecule = [Molecule(SMILES=smiles)]
            species.generateResonanceIsomers()
            identified_labels.append(species_label)

        logging.info("Identified {0} species:".format(len(identified_labels)))
        for species_label in identified_labels:
            logging.info("   {0}".format(species_label))

        self.identified_labels.extend(identified_labels)

    def askForMatchSMILES(self, chemkinSpecies):
            species_label = chemkinSpecies.label
            formula = self.formulaDict[species_label]
            print "Species {species} has formula {formula}".format(species=species_label, formula=formula)
            smiles = raw_input('What is its SMILES?\n')
            while formula != Molecule(SMILES=smiles).getFormula():
                smiles = raw_input("SMILES {0} has formula {1} not required formula {2}. Try again:\n".format(smiles, Molecule(SMILES=smiles).getFormula(), formula))
            self.smilesDict[species_label] = smiles
            species = self.speciesDict[species_label]
            species.molecule = [Molecule(SMILES=smiles)]
            species.generateResonanceIsomers()
            self.identified_labels.append(species_label)

    def askForMatchID(self, speciesLabel, possibleMatches):
            """
            Ask user for a match for a given speciesLabel, choosing from the iterable possibleMatches
            """
            # print "Species {species} has formula {formula}".format(species=species_label, formula=formula)
            matchesDict = {species.index: species for species in possibleMatches }
            possibleIndicesStr = [str(i) for i in sorted(matchesDict.keys())]
            print "Species {0} could be one of:".format(speciesLabel)
            for index in sorted(matchesDict.keys()):
                rmgSpecies = matchesDict[index]
                dH = self.getEnthalpyDiscrepancy(speciesLabel, rmgSpecies)
                Nvotes = len(self.votes[speciesLabel][rmgSpecies])
                allPossibleChemkinSpecies = [ck for ck, matches in self.votes.iteritems() if rmgSpecies in matches]
                print "{0:6d} {1:18s} {2:8.1f} kJ/mol  ({3} votes) {4!s}".format(index, rmgSpecies.label, dH, Nvotes, allPossibleChemkinSpecies)
            chosenID = raw_input('What is it? (see voting info above)\n')
            while chosenID not in possibleIndicesStr:
                chosenID = raw_input("That wasn't one of {0}. Try again:\n".format(','.join(possibleIndicesStr)))

            rmgSpecies = matchesDict[int(chosenID)]
            logging.info("Based on user input, matching {0} with {1!s}".format(speciesLabel, rmgSpecies))
            self.setMatch(speciesLabel, rmgSpecies)
            return speciesLabel, rmgSpecies

    def edgeReactionsMatching(self, chemkinReaction):
        """A generator giving edge reactions that match the given chemkin reaction"""
        reactionsMatch = self.reactionsMatch
        for edgeReaction in self.rmg_object.reactionModel.edge.reactions:
            if reactionsMatch(edgeReaction, chemkinReaction):
                yield edgeReaction

    def chemkinReactionsMatching(self, rmgReaction):
        """A generator giving chemkin reactions that match the given RMG reaction"""
        reactionsMatch = self.reactionsMatch
        for chemkinReaction in self.chemkinReactions:
            if reactionsMatch(rmgReaction, chemkinReaction):
                yield chemkinReaction

    def drawSpecies(self, rmg_species):
        "Draw a species, saved in 'species' directory named after its RMG name (label and id)."
        # Draw molecules if necessary
        fstr = os.path.join(self.rmg_object.outputDirectory, 'species', '{0!s}.png'.format(rmg_species))
        if not os.path.exists(fstr):
            MoleculeDrawer().draw(rmg_species.molecule[0], 'png', fstr)

    def drawAllCandidateSpecies(self):
        """Draws all the species that are in self.prunedVotes"""
        candidateSpecies = set()
        for possibleMatches in self.votes.itervalues():
            candidateSpecies.update(possibleMatches.keys())
        for rmg_species in candidateSpecies:
            self.drawSpecies(rmg_species)

    def moveSpeciesDrawing(self, rmg_species):
        "Move a species drawing from 'species' directory to 'species/MATCHED' directory."
        source = os.path.join(self.rmg_object.outputDirectory, 'species', '{0!s}.png'.format(rmg_species))
        destination = os.path.join(self.rmg_object.outputDirectory, 'species', 'MATCHED', '{0!s}.png'.format(rmg_species))
        if os.path.exists(source):
            os.renames(source, destination)

    def getEnthalpyDiscrepancy(self, chemkinLabel, rmgSpecies):
        """
        Return the difference in enthalpy at 800K in kJ/mol
        """
        return (self.thermoDict[chemkinLabel].getEnthalpy(800.) - rmgSpecies.thermo.getEnthalpy(800.)) / 1000.

    def setMatch(self, chemkinLabel, rmgSpecies):
        """Store a match, once you've identified it"""
        self.identified_labels.append(chemkinLabel)
        self.identified_unprocessed_labels.append(chemkinLabel)
        self.speciesDict_rmg[chemkinLabel] = rmgSpecies
        enthalpyDiscrepancy = self.getEnthalpyDiscrepancy(chemkinLabel, rmgSpecies)
        logging.info("Storing match: {0} = {1!s}".format(chemkinLabel, rmgSpecies))
        logging.info("  On match, Enthalpies at 800K differ by {0:.1f} kJ/mol".format(enthalpyDiscrepancy))
        display(rmgSpecies)
        self.moveSpeciesDrawing(rmgSpecies)
        rmgSpecies.label = chemkinLabel
        with open(self.dictionaryFile, 'a') as f:
            f.write("{0}\t{1}\t{2:.1f}\n".format(chemkinLabel, rmgSpecies.molecule[0].toSMILES(), enthalpyDiscrepancy))
        self.drawSpecies(rmgSpecies)

        # For kinetics purposes, we convert the thermo to Wilhoit
        # This allows us to extrapolating H to 298 to find deltaH rxn
        # for ArrheniusEP kinetics,
        # and to 0K so we can do barrier height checks with E0.
        Cp0 = rmgSpecies.calculateCp0()
        CpInf = rmgSpecies.calculateCpInf()
        thermo = self.thermoDict[chemkinLabel]
        # pretend it was valid down to 298 K
        oldLowT = thermo.Tmin.value_si
        thermo.selectPolynomial(thermo.Tmin.value_si).Tmin.value_si = min(298.0, thermo.Tmin.value_si)
        newThermo = thermo.toWilhoit(Cp0=Cp0, CpInf=CpInf)
        thermo.comment += "\nLow T polynomial Tmin changed from {0} to {1} K when importing to RMG".format(oldLowT, 298.0)
        # thermo.selectPolynomial(thermo.Tmin.value_si).Tmin.value_si = oldLowT  # put it back
        self.thermoDict[chemkinLabel].E0 = newThermo.E0

    def getInvalidatedReactionsAndRemoveVotes(self, chemkinLabel, rmgSpecies):
        """Remove the votes, and return the list of voting reactions."""
        # Remove both chemkinLabel and rmgSpecies from the voting dictionaries
        reactionsToReCheck = set()
        possibles = self.votes.pop(chemkinLabel, {})
        for rxns in possibles.itervalues():
            for rxn in rxns:
                reactionsToReCheck.add(rxn[1])
        for ck, possibles in self.votes.iteritems():
            for rxn in possibles.pop(rmgSpecies, {}):
                reactionsToReCheck.add(rxn[1])
        return reactionsToReCheck

    def checkReactionsForMatches(self, reactionsToCheck):
        """
        Checks the given list of edge reactions for corresponding
        chemkin reactions in the self.chemkinReactionsUnmatched list.
        Updates the votes in self.votes,
        and removes fully matched things from self.chemkinReactionsUnmatched
        also prunes reactions from self.rmg_object.edge if they can't match anything ever.
        """
        chemkinReactionsUnmatched = self.chemkinReactionsUnmatched
        reactionsMatch = self.reactionsMatch
        votes = self.votes

        reactionsToPrune = set()
        for edgeReaction in reactionsToCheck:
            edgeReactionMatchesSomething = False
            for chemkinReaction in chemkinReactionsUnmatched:
                self.suggestedMatches = {}
                if reactionsMatch(edgeReaction, chemkinReaction):
                    edgeReactionMatchesSomething = True
                    logging.info("Chemkin reaction     {0}\n matches RMG reaction  {1}".format(chemkinReaction, edgeReaction))
                    if self.suggestedMatches:
                        logging.info(" suggesting new species match: {0!r}".format({l:str(s) for l, s in self.suggestedMatches.iteritems()}))
                    else:
                        logging.info(" suggesting no new species matches.")
                        for reagents in (chemkinReaction.reactants, chemkinReaction.products):
                            for reagent in reagents:
                                if reagent.label not in self.identified_labels:
                                    break
                            else:  # didn't break inner loop so these reagents have all been identified
                                continue
                            break  # did break inner loop, so break outer loop as there's an unidentified species
                        else:  # didn't break outer loop, so all species have been identified
                            # remove it from the list of useful reactions.
                            chemkinReactionsUnmatched.remove(chemkinReaction)
                    for chemkinLabel, rmgSpecies in self.suggestedMatches.iteritems():
                        if chemkinLabel not in votes:
                            votes[chemkinLabel] = {rmgSpecies: set([(chemkinReaction, edgeReaction)])}
                        else:
                            if rmgSpecies not in votes[chemkinLabel]:
                                votes[chemkinLabel][rmgSpecies] = set([(chemkinReaction, edgeReaction)])
                            else:
                                votes[chemkinLabel][rmgSpecies].add((chemkinReaction, edgeReaction))
                        # now votes is a dict of dicts of lists {'ch3':{<Species CH3>: [ voting_reactions ]}}
            if not edgeReactionMatchesSomething:
                reactionsToPrune.add(edgeReaction)
        # remove those reactions
        logging.info("Removing {0} edge reactions that didn't match anything.".format(len(reactionsToPrune)))
        prune = self.rmg_object.reactionModel.edge.reactions.remove
        for rxn in reactionsToPrune:
            try:
                prune(rxn)
            except ValueError:
                logging.info("Reaction {0!s} was not in edge! Could not remove it.".format(rxn))

    def pruneVoting(self):
        """
        Return a voting matrix with only significant (unique) votes.
        
        If the same reaction is voting for several species, remove it.
        If a match has a large enthalpy discrepancy, remove it.
        """
        votes = self.votes
        prunedVotes = {}

        # votes matrix containing sets with only the chemkin reactions, not the corresponding RMG reactions
        ckVotes = {chemkinLabel:
                   {matchingSpecies:
                    set([r[0] for r in votingReactions]) for matchingSpecies, votingReactions in possibleMatches.iteritems()
                   } for chemkinLabel, possibleMatches in votes.iteritems() }

        for chemkinLabel, possibleMatches in ckVotes.iteritems():
            for rmgSpecies in possibleMatches.keys():
                dH = self.getEnthalpyDiscrepancy(chemkinLabel, rmgSpecies)
                if abs(dH) > 150:
                    logging.info("Removing possible match {0} : {1!s}  because enthalpy discrepancy is {2:.1f} kJ/mol".format(chemkinLabel, rmgSpecies, dH))
                    del(possibleMatches[rmgSpecies])

        for chemkinLabel, possibleMatches in ckVotes.iteritems():
            if len(possibleMatches) == 0:
                logging.info("No remaining matches for {0}".format(chemkinLabel))
                continue
            if len(possibleMatches) == 1:
                prunedVotes[chemkinLabel] = possibleMatches
                continue
            commonVotes = None
            mostVotes = 0
            for matchingSpecies, votingReactions in possibleMatches.iteritems():
                mostVotes = max(mostVotes, len(votingReactions))
                if commonVotes is None:
                    commonVotes = set(votingReactions)  # make a copy!!
                else:
                    commonVotes.intersection_update(votingReactions)
            if len(commonVotes) < mostVotes:
                logging.info("Removing {0} voting reactions that are common to all {1} matches for {2}".format(
                                len(commonVotes), len(possibleMatches), chemkinLabel))
                prunedVotes[chemkinLabel] = { matchingSpecies: votingReactions.difference(commonVotes) for matchingSpecies, votingReactions in possibleMatches.iteritems() if votingReactions.difference(commonVotes)}
            else:
                prunedVotes[chemkinLabel] = possibleMatches
        self.prunedVotes = prunedVotes
        return prunedVotes

    def printVoting(self, votes):
        """
        Log the passed in voting matrix
        """
        logging.info("Current voting:::")
        chemkinControversy = {label: 0 for label in votes.iterkeys()}
        rmgControversy = {}
        flatVotes = {}
        for chemkinLabel, possibleMatches in votes.iteritems():
            for matchingSpecies, votingReactions in possibleMatches.iteritems():
                self.drawSpecies(matchingSpecies)
                flatVotes[(chemkinLabel, matchingSpecies)] = votingReactions
                chemkinControversy[chemkinLabel] += len(votingReactions)
                rmgControversy[matchingSpecies] = rmgControversy.get(matchingSpecies, 0) + len(votingReactions)

        for chemkinLabel in sorted(chemkinControversy.keys(), key=lambda label:-chemkinControversy[label]):
            possibleMatches = votes[chemkinLabel]
            logging.info("{0} matches {1} RMG species:".format(chemkinLabel, len(possibleMatches)))
            for matchingSpecies in sorted(possibleMatches.iterkeys(), key=lambda species:-len(possibleMatches[species])) :
                votingReactions = possibleMatches[matchingSpecies]
                logging.info("  {0}  matches  {1!s}  according to {2} reactions:".format(chemkinLabel, matchingSpecies, len(votingReactions)))
                logging.info("  Enthalpies at 800K differ by {0:.1f} kJ/mol".format((self.thermoDict[chemkinLabel].getEnthalpy(800) - matchingSpecies.thermo.getEnthalpy(800)) / 1000.))
                display(matchingSpecies)
                for rxn in votingReactions:
                    if isinstance(rxn, tuple):
                        logging.info("    {0!s}     //    {1!s}".format(rxn[0], rxn[1]))
                    else:
                        logging.info("    {0!s}".format(rxn))

    def main(self):
        """This is the main matcher function that does the whole thing"""
        args = self.args
        species_file = args.species
        reactions_file = args.reactions or species_file
        thermo_file = args.thermo

        outputThermoFile = os.path.splitext(thermo_file)[0] + '.thermo.py'

        self.loadModel(species_file, reactions_file, thermo_file)

        logging.info("Initializing RMG")
        self.initializeRMG(args)
        rm = self.rmg_object.reactionModel
        self.dictionaryFile = os.path.join(args.output_directory, 'MatchedSpeciesDictionary.txt')

        with open(self.dictionaryFile, 'w') as f:
            f.write("Species name\tSMILES\tEnthaply discrepancy at 800K\n")

        self.identifySmallMolecules()

        logging.info("Importing identified species into RMG model")
        # Add identified species to the reaction model complete species list
        newSpeciesDict = {}
        for species_label in self.identified_labels:
            old_species = self.speciesDict[species_label]
            logging.info(species_label)
            rmg_species, wasNew = rm.makeNewSpecies(old_species, label=old_species.label)
            assert wasNew, "Species with structure of '{0}' already created with label '{1}'".format(species_label, rmg_species.label)
            newSpeciesDict[species_label] = rmg_species
            if self.formulaDict[species_label] in {'N2'}:
                rmg_species.reactive = False
            rmg_species.generateThermoData(self.rmg_object.database)
        # Set match using the function to get all the side-effects.
        labelsToProcess = self.identified_labels
        self.identified_labels = []
        for chemkinLabel in labelsToProcess:
            # this adds it back into self.identified_labels
            self.setMatch(chemkinLabel, newSpeciesDict[chemkinLabel])

        chemkinFormulas = set(self.formulaDict.values())

        self.loadReactions(reactions_file)
        chemkinReactionsUnmatched = self.chemkinReactionsUnmatched
        votes = self.votes

        reactionsToCheck = set()
        while self.identified_unprocessed_labels:
            labelToProcess = self.identified_unprocessed_labels.pop(0)
            logging.info("Processing species {0}...".format(labelToProcess))
            if self.formulaDict[labelToProcess] in ('N2', 'Ar'):
                logging.info("Not processing {0} because I can't react it in RMG".format(labelToProcess))
                continue

            # Add species to RMG core.
            rm.enlarge(self.speciesDict_rmg[labelToProcess])

            # do a partial prune of new reactions that definitely aren't going to be useful
            reactionsToPrune = set()
            for newSpecies in rm.newSpeciesList:
                if newSpecies.molecule[0].getFormula() in chemkinFormulas:
                    continue
                # else it's not useful to us
                # identify any reactions it's involved in
                for rxn in rm.newReactionList:
                    if newSpecies in rxn.reactants or newSpecies in rxn.products:
                        reactionsToPrune.add(rxn)
            logging.info("Removing {0} edge reactions that aren't useful".format(len(reactionsToPrune)))
            # remove those reactions
            for rxn in reactionsToPrune:
                rm.edge.reactions.remove(rxn)
                rm.newReactionList.remove(rxn)
            reactionsToPrune.clear()

            logging.info("Adding {0} new RMG reactions to be checked.".format(len(rm.newReactionList)))
            reactionsToCheck.update(rm.newReactionList)
            logging.info("In total will check {0} edge reactions".format(len(reactionsToCheck)))
            logging.info("against {0} unmatched chemkin reactions.".format(len(chemkinReactionsUnmatched)))

            if len(self.identified_unprocessed_labels) == 0:
                logging.info("** Running out of things to process!")

            while reactionsToCheck:
                self.checkReactionsForMatches(reactionsToCheck)
                # Have just checked all those reactions, so clear the reactionsToCheck,
                # ready to start adding to it again based on new matches.
                reactionsToCheck.clear()

                # self.printVoting(votes)
                prunedVotes = self.pruneVoting()
                # self.printVoting(prunedVotes)

                self.drawAllCandidateSpecies()

                newMatches = []
                for chemkinLabel, possibleMatches in prunedVotes.iteritems():
                    if len(possibleMatches) == 1:
                        matchingSpecies, votingReactions = possibleMatches.items()[0]
                        logging.info("\nOnly one suggested match for {0}: {1!s}".format(chemkinLabel, matchingSpecies))
                        display(matchingSpecies)
                        logging.info("With {0} unique voting reactions:".format(len(votingReactions)))
                        for reaction in votingReactions:
                            logging.info("  {0!s}".format(reaction))
                        allPossibleChemkinSpecies = [ck for ck, matches in prunedVotes.iteritems() if matchingSpecies in matches]
                        if len(allPossibleChemkinSpecies) == 1:
                            logging.info("Only one chemkin species has this match (after pruning).")
                            self.setMatch(chemkinLabel, matchingSpecies)
                            newMatches.append((chemkinLabel, matchingSpecies))
                        else:
                            logging.info("Other Chemkin species that also match {0} (after pruning) are {1!r}".format(matchingSpecies.label, allPossibleChemkinSpecies))
                            logging.info("Will not make match at this time.")

                for chemkinLabel, matchingSpecies in newMatches:
                    invalidatedReactions = self.getInvalidatedReactionsAndRemoveVotes(chemkinLabel, matchingSpecies)
                    reactionsToCheck.update(invalidatedReactions)
                logging.info("After making {0} matches, will have to re-check {1} edge reactions".format(len(newMatches), len(reactionsToCheck)))


            logging.info("Finished processing species {0}!".format(labelToProcess))
            logging.info("Have now identified {0} of {1} species ({2:.1%}).".format(len(self.identified_labels), len(self.speciesList), float(len(self.identified_labels)) / len(self.speciesList)))
            logging.info("And fully identified {0} of {1} reactions ({2:.1%}).".format(len(self.chemkinReactions) - len(self.chemkinReactionsUnmatched), len(self.chemkinReactions), 1 - float(len(self.chemkinReactionsUnmatched)) / len(self.chemkinReactions)))
            logging.info("Still to process {0} matches: {1!r}".format(len(self.identified_unprocessed_labels), self.identified_unprocessed_labels))

            if len(self.identified_unprocessed_labels) == 0 and self.prunedVotes and not self.manualMatchesToProcess:
                logging.info("Waiting for input from the web front end..")
                while not self.manualMatchesToProcess:
                    time.sleep(1)

            while self.manualMatchesToProcess:
                chemkinLabel, matchingSpecies = self.manualMatchesToProcess.pop(0)
                logging.info("There is a manual match to process: {0} is {1!s}".format(chemkinLabel, matchingSpecies))
                self.setMatch(chemkinLabel, matchingSpecies)
                invalidatedReactions = self.getInvalidatedReactionsAndRemoveVotes(chemkinLabel, matchingSpecies)
                reactionsToCheck.update(invalidatedReactions)
                logging.info("After making that match, will have to re-check {0} edge reactions".format(len(reactionsToCheck)))

            if len(self.identified_unprocessed_labels) == 0 and self.votes:
                self.printVoting(prunedVotes)
                logging.info("Run out of options. Asking for help!")
                speciesLabel = raw_input('Which label would you like to identify? (see voting info above)\n')
                while True:
                    if speciesLabel not in self.formulaDict:
                        print("That's not a known species label")
                    elif speciesLabel in self.identified_labels:
                        print("That's already been identified")
                    elif speciesLabel not in votes:
                        print("We have no candidate matches for that label.")
                    else:  # label is valid, break out of while loop.
                        break
                    speciesLabel = raw_input("Try again:\n")
                possibleMatches = votes[speciesLabel].keys()
                chemkinLabel, matchingSpecies = self.askForMatchID(speciesLabel, possibleMatches)
                invalidatedReactions = self.getInvalidatedReactionsAndRemoveVotes(chemkinLabel, matchingSpecies)
                reactionsToCheck.update(invalidatedReactions)
                logging.info("After making that match, will have to re-check {0} edge reactions".format(len(reactionsToCheck)))


        print "Finished reading"
        with open(outputThermoFile, 'w') as f:
            counter = 0
            for species in self.speciesList:
                counter += 1
                print counter, species,
                if species.label not in self.speciesDict_rmg:
                    print ""
                    continue  # don't save unidentified species
                print "\t IDENTIFIED"
                entry = Entry()
                entry.index = counter
                entry.label = species.label
                # molecule = Molecule(SMILES=self.smilesDict.get(species.label, 'C'))
                molecule = self.speciesDict_rmg.get(species.label, Species().fromSMILES('C')).molecule[0]
                entry.item = molecule
                entry.data = self.thermoDict[species.label]
                entry.longDesc = getattr(species.thermo, 'comment', '') + 'Imported from {source}'.format(source=thermo_file)
                user = getUsername()
                event = [time.asctime(), user, 'action', '{user} imported this entry from {source}'.format(user=user, source=thermo_file)]
                entry.history = [event]
                saveEntry(f, entry)

        print "done"

    def _img(self, species):
        """Get the html tag for the image of a species"""
        imagesPath = 'img'  # to serve via cherryPy
        #imagesPath = 'file://'+os.path.abspath(os.path.join(self.args.output_directory,'species')) # to get from disk
        return "<img src='{1}/{0!s}.png' title='{0}'>".format(species, imagesPath)

    @cherrypy.expose
    def index(self):
        return self.html_head + """
<h1>Mechanism importer</h1>
<ul>
<li><a href="identified.html">List of identified species.</li>
<li><a href="votes.html">Voting reactions.</li>
</ul>
        """ + self.html_tail

    @cherrypy.expose
    def identified_html(self):
        img = self._img
        return (self.html_head + '<h1>Identified Species</h1><table style="width:500px"><tr>' +
                "</tr>\n<tr>".join(["<td>{number}</td><td>{label}</td><td>{img}</td>".format(img=img(self.speciesDict_rmg[lab]), label=lab, number=n + 1) for n, lab in enumerate(self.identified_labels)]) +
                '</tr></table>' + self.html_tail)
    @cherrypy.expose
    def identified_json(self):
        return json.dumps(self.identified_labels)

    @cherrypy.expose
    def votes_html(self):
        votes = self.votes.copy()
        img = self._img
        chemkinControversy = {label: 0 for label in votes.iterkeys()}
        rmgControversy = {}
        flatVotes = {}

        labelsWaitingToProcess = [item[0] for item in self.manualMatchesToProcess]
        speciesWaitingToProcess = [item[1] for item in self.manualMatchesToProcess]
        # to turn reactions into pictures
        searcher = re.compile('(\S+\(\d+\))\s')
        def replacer(match):
            return self._img(match.group(1))

        for chemkinLabel, possibleMatches in votes.iteritems():
            for matchingSpecies, votingReactions in possibleMatches.iteritems():
                flatVotes[(chemkinLabel, matchingSpecies)] = votingReactions
                chemkinControversy[chemkinLabel] += len(votingReactions)
                rmgControversy[matchingSpecies] = rmgControversy.get(matchingSpecies, 0) + len(votingReactions)
        output = [self.html_head]
        output.append("<h1>Votes</h1>")
        for chemkinLabel in sorted(chemkinControversy.keys(), key=lambda label:-chemkinControversy[label]):
            if chemkinLabel in labelsWaitingToProcess:
                output.append("<hr><h2>{0} has just been identified but not yet processed.</h2>".format(chemkinLabel))
                continue
            possibleMatches = votes[chemkinLabel]
            output.append("<hr><h2>{0} matches {1} RMG species</h2>".format(chemkinLabel, len(possibleMatches)))
            for matchingSpecies in sorted(possibleMatches.iterkeys(), key=lambda species:-len(possibleMatches[species])) :
                if matchingSpecies in speciesWaitingToProcess:
                    output.append("{img} which has just been identified but not yet processed.<br>".format(img=img(matchingSpecies)))
                    continue
                votingReactions = possibleMatches[matchingSpecies]
                output.append("<a href='/match.html?ckLabel={ckl}&rmgLabel={rmgl}'>{img}</a>  according to {n} reactions. ".format(ckl=urllib.quote_plus(chemkinLabel), rmgl=urllib.quote_plus(str(matchingSpecies)), img=img(matchingSpecies), n=len(votingReactions)))
                output.append("  Enthalpies at 800K differ by {0:.1f} kJ/mol<br>".format((self.thermoDict[chemkinLabel].getEnthalpy(800) - matchingSpecies.thermo.getEnthalpy(800)) / 1000.))
                output.append('<table  style="width:800px">')
                for n, rxn in enumerate(votingReactions):
                    if isinstance(rxn, tuple):
                        rmgrxn = str(rxn[1])
                        rmgRxnPics = searcher.sub(replacer, rmgrxn + ' ')
                        output.append("<tr><td>{0}</td><td> {1!s}   </td><td>  {2!s} </td></tr>".format(n + 1, rxn[0], rmgRxnPics))
                    else:
                        output.append("<tr><td>{0}</td><td> {1!s}</td></tr>".format(n + 1, rxn))
                output.append("</table>")
        output.append(self.html_tail)
        return '\n'.join(output)

    @cherrypy.expose
    def match_html(self, ckLabel=None, rmgLabel=None):
        if ckLabel not in self.votes:
            return "ckLabel not valid"
        for rmgSpecies in self.votes[ckLabel].iterkeys():
            if str(rmgSpecies) == rmgLabel:
                self.manualMatchesToProcess.append((str(ckLabel), rmgSpecies))
                break
        else:
            return "rmgLabel not a candidate for that ckLabel"
        ## Wait for it to be processed:
        #while self.manualMatchesToProcess:
        #    time.sleep(1)
        raise cherrypy.HTTPRedirect("/votes.html")

    @cherrypy.expose
    def progress_json(self):
        total = len(self.speciesList)
        identified = len(self.identified_labels) + len(self.manualMatchesToProcess)
        unprocessed = len(self.identified_unprocessed_labels) + len(self.manualMatchesToProcess)
        answer = {'processed': identified - unprocessed,
                  'unprocessed': unprocessed,
                  'unidentified': total - identified
        }
        return json.dumps(answer)

    html_head = """
<html>
<head>
<script src="http://ajax.googleapis.com/ajax/libs/jquery/1.9.1/jquery.min.js"></script>
<script>
function updateStats() {
    $.getJSON( "progress.json", function( json ) {
            console.log('Updating stats.. Unidentified now ' + json.unidentified );
            $('#processed').html(json.processed).width(json.processed);
            $('#unprocessed').html(json.unprocessed).width(json.unprocessed);
            $('#unidentified').html(json.unidentified).width(json.unidentified);
            repeater = setTimeout(updateStats, 10000); // do again in 10 seconds
        }).fail(function( jqxhr, textStatus, error ) {
              var err = textStatus + ', ' + error;
              console.log( "Request Failed: " + err);
        });
}
$( document ).ready(function() {
    updateStats();
});
</script>
<style>
#processed {background-color: #33ff33;}
#unprocessed {background-color: yellow;}
#unidentified {background-color: red;}
</style>    
</head>

<body>
    <table width=100%><tr>
        <td id="processed"></td>
        <td id="unprocessed"></td>
        <td id="unidentified"></td>
    </tr>
    </table>
    """
    html_tail = """
    </body></html>
    """


if __name__ == '__main__':

    # Parse the command-line arguments (requires the argparse module)
    args = parseCommandLineArguments()
    os.path.exists(args.output_directory) or os.makedirs(args.output_directory)

    # Initialize the logging system (resets the RMG.log file)
    level = logging.INFO
    if args.debug: level = 0
    elif args.verbose: level = logging.DEBUG
    elif args.quiet: level = logging.WARNING
    initializeLog(level, os.path.join(args.output_directory, 'RMG.log'))

    mm = ModelMatcher(args)

    t = threading.Thread(target=mm.main)
    t.daemon = True
    t.start()
    import webbrowser
    webbrowser.open('http://127.0.0.1:8080')
    cherrypy.server.socket_host = '0.0.0.0'
    cherrypy.config.update({'environment': 'production',
                            'log.error_file': 'site.log',
                            'log.screen': False})

    conf = {'/img': {'tools.staticdir.on': True,
                      'tools.staticdir.dir': os.path.join(args.output_directory, 'species'),
            }}
    cherrypy.quickstart(mm, '/', config=conf)



    # mm.main(args)
