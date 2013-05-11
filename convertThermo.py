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
"""

import os.path
import argparse
import logging

import rmgpy
import rmgpy.rmg

from rmgpy.chemkin import loadChemkinFile, readSpeciesBlock, readThermoBlock, readReactionsBlock, removeCommentFromLine
from rmgpy.reaction import ReactionModel

from rmgpy.data.thermo import Entry, saveEntry
from rmgpy.molecule import Molecule
from rmgpy.rmg.model import Species  # you need this one, not the one in rmgpy.species!

from rmgpy.rmg.main import RMG, initializeLog

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
    def __init__(self):
        self.speciesDict = {}
        self.thermoDict = {}
        self.formulaDict = {}
        self.smilesDict = {}
        self.identified = []
        self.identified_unprocessed = []
        self.speciesList = None

    def loadModel(self, species_file, reactions_file, thermo_file):
        print 'Loading model...'
        model = ReactionModel()

        speciesAliases = {}
        speciesDict = {}
        if species_file:
            logging.info("Reading species list...")
            speciesList = []
            with open(species_file) as f:
                line0 = f.readline()
                while line0 != '':
                    line = removeCommentFromLine(line0)[0]
                    tokens = line.split()
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

        with open(thermo_file) as f:
            logging.info("Reading thermo...")
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
        rmg.thermoLibraries = []
        rmg.reactionLibraries = [('Glarborg/HighP', False)]
        rmg.loadDatabase()
        logging.info("Loaded database.")

        rmg.reactionModel = rmgpy.rmg.model.CoreEdgeReactionModel()
        rmg.reactionModel.kineticsEstimator = 'rate rules'
        rmg.initialSpecies = []
        rmg.reactionSystems = []

        self.rmg_object = rmg
        return rmg

    def main(self, args):
        """This is the main matcher function that does the whole thing"""
        species_file = args.species
        reactions_file = args.reactions or species_file
        thermo_file = args.thermo

        outputThermoFile = os.path.splitext(thermo_file)[0] + '.thermo.py'

        self.loadModel(species_file, reactions_file, thermo_file)

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
             'CO': '[C]=[O]',
             'HO2': '[O]O',
             'C3H8': 'CCC',
             'CH': '[CH]',
             }

        identified = self.identified
        identified_unprocessed = self.identified_unprocessed
        # use speciesList if it is not None or empty, else the formulaDict keys.
        for species_label in [s.label for s in self.speciesList or []] or self.formulaDict.keys():
            formula = self.formulaDict[species_label]
            # print "Species {species} has formula {formula}".format(species=species_label, formula=formulaString)
            if formula in known_formulas:
                known_smiles = known_formulas[formula]
                logging.info("I think its SMILES is {0}".format(known_smiles))
                smiles = known_smiles
                # print "Hit Enter to confirm, or type the new smiles if wrong\n"
                # smiles = raw_input() or known_smiles
            else:
                continue  # Remove this line to input all SMILES strings
                smiles = raw_input('What is its SMILES?\n')
            self.smilesDict[species_label] = smiles
            while formula != Molecule(SMILES=smiles).getFormula():
                smiles = raw_input("SMILES {0} has formula {1} not required formula {2}. Try again:\n".format(smiles, Molecule(SMILES=smiles).getFormula(), formula))
            species = self.speciesDict[species_label]
            species.molecule = [Molecule(SMILES=smiles)]
            species.generateResonanceIsomers()
            identified.append(species_label)
            identified_unprocessed.append(species_label)

        logging.info("Identified {0} species:".format(len(identified)))
        for species_label in identified:
            logging.info("   {0}".format(species_label))

        logging.info("Reading reactions.")
        with open(reactions_file) as f:
            reactionList = readReactionsBlock(f, self.speciesDict, readComments=True)
        logging.info("Read {0} reactions from chemkin file.".format(len(reactionList)))

        logging.info("Initializing RMG")
        self.initializeRMG(args)
        rm = self.rmg_object.reactionModel

        logging.info("Importing identified species into RMG model")
        # Add identified species to the reaction model complete species list
        newSpeciesDict = {}
        for species_label in identified:
            old_species = self.speciesDict[species_label]
            new_species, wasNew = rm.makeNewSpecies(old_species, label=old_species.label)
            assert wasNew, "Species with structure of '{0}' already created with label '{1}'".format(species_label, new_species.label)
            # For kinetics purposes, we convert the thermo to Wilhoit
            # This allows us to extrapolating H to 298 to find deltaH rxn
            # for ArrheniusEP kinetics,
            # and to 0K so we can do barrier height checks with E0.
            thermo = old_species.thermo
            # pretend it was valid down to 298 K
            thermo.selectPolynomial(thermo.Tmin.value_si).Tmin.value_si = 298
            Cp0 = new_species.calculateCp0()
            CpInf = new_species.calculateCpInf()
            new_species.thermo = old_species.thermo.toWilhoit(Cp0=Cp0, CpInf=CpInf)
            newSpeciesDict[species_label] = new_species

        # : Dictionary of species that are in RMG core/edge reaction model, with structures
        self.speciesDict_rmg = newSpeciesDict
        rm.enlarge(newSpeciesDict['ch4'])
        rm.edge.reactions
        import ipdb; ipdb.set_trace()

        print "Finished reading"
        with open(outputThermoFile, 'w') as f:
            counter = 0
            for species in self.speciesList:
                counter += 1
                print counter, species
                entry = Entry()
                entry.index = counter
                entry.label = species.label
                molecule = Molecule(SMILES=self.smilesDict.get(species.label, 'C'))
                entry.item = molecule
                entry.data = self.thermoDict[species.label]
                entry.longDesc = getattr(species.thermo, 'comment', '') + 'Imported from {source}'.format(source=thermo_file)
                user = getUsername()
                event = [time.asctime(), user, 'action', '{user} imported this entry from {source}'.format(user=user, source=thermo_file)]
                entry.history = [event]
                saveEntry(f, entry)

        print "done"

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

    mm = ModelMatcher()
    mm.main(args)


