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

from rmgpy.chemkin import loadChemkinFile, readSpeciesBlock, readThermoBlock, readReactionsBlock, removeCommentFromLine
from rmgpy.reaction import ReactionModel

from rmgpy.data.thermo import Entry, saveEntry
from rmgpy.molecule import Molecule
from rmgpy.species import Species

import time
import sys
# Put the RMG-database project at the start of the python path, so we use that importOldDatabase script!
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'RMG-database')))
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
    given a formula in dict form {'c':2, 'h':6, 'o':0}
    return a string "C2H6"
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

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--species', metavar='FILE', type=str, nargs='?', default=None,
        help='the Chemkin file containing the list of species')
    parser.add_argument('--thermo', metavar='FILE', type=str,
        help='the Chemkin files containing the thermo')
    args = parser.parse_args()

    species_file = args.species
    thermo_file = args.thermo

    outputThermoFile = os.path.splitext(thermo_file)[0] + '.thermo.py'

    print 'Loading model...'
    model = ReactionModel()

    speciesAliases = {}
    speciesDict = {}
    if species_file:
        speciesList = []
        with open(species_file) as f:
            line0 = f.readline()
            while line0 != '':
                line = removeCommentFromLine(line0)[0]
                line = line.strip()
                tokens = line.split()
                tokens_upper = line.upper().split()

                if 'SPECIES' in line.upper():
                    print "Reading species..."
                    # Unread the line (we'll re-read it in readReactionBlock())
                    f.seek(-len(line0), 1)
                    readSpeciesBlock(f, speciesDict, speciesAliases, speciesList)
                line0 = f.readline()
    else:
        logging.info("No species file to limit species. Will read everything in thermo")
        speciesList = None
        speciesDict = MagicSpeciesDict(speciesDict)

    with open(thermo_file) as f:
        line0 = f.readline()
        while line0 != '':
            line = removeCommentFromLine(line0)[0]
            line = line.strip()
            tokens = line.split()
            tokens_upper = line.upper().split()
            if 'THERM' in line.upper():
                print "Reading thermo..."
                # Unread the line (we'll re-read it in readThermoBlock())
                f.seek(-len(line0), 1)
                formulaDict = readThermoBlock(f, speciesDict)
            line0 = f.readline()

    known = {
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
    smilesDict = {}
    for species in [s.label for s in speciesList] or formulaDict.keys():
        speciesUpper = species.upper()
        formula = formulaDict[speciesUpper]
        formulaString = convertFormula(formula)
        print "Species {species} has formula {formula}".format(species=species, formula=formulaString)
        if formulaString in known:
            known_smiles = known[formulaString]
            print "I think its SMILES is {0}".format(known_smiles)
            print "Hit Enter to confirm, or type the new smiles if wrong\n"
            smiles = raw_input() or known_smiles
        else:
            continue  # Remove this line to input all SMILES strings
            smiles = raw_input('What is its SMILES?\n')
        smilesDict[species] = smiles
        while formulaString != Molecule(SMILES=smiles).getFormula():
            smiles = raw_input("SMILES {0} has formula {1} not required formula {2}. Try again:\n".format(smiles, Molecule(SMILES=smiles).getFormula(), formulaString))

    with open(outputThermoFile, 'w') as f:
        counter = 0
        for species in speciesList:
            counter += 1
            print counter, species
            entry = Entry()
            entry.index = counter
            entry.label = species.label
            molecule = Molecule(SMILES=smilesDict.get(species, 'C'))
            entry.item = molecule
            entry.data = species.thermo
            entry.longDesc = getattr(species.thermo, 'comment', '') + 'Imported from {source}'.format(source=thermo_file)
            user = getUsername()
            event = [time.asctime(), user, 'action', '{user} imported this entry from {source}'.format(user=user, source=thermo_file)]
            entry.history = [event]
            saveEntry(f, entry)

    print "done"
