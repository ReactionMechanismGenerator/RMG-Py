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

import os.path
import logging
import numpy

import chempy.constants as constants
from chempy.states import *

from base import *

################################################################################

class GroupFrequency:
    """
    Represent a set of characteristic frequencies for a group in the frequency
    database. These frequencies are stored in the `frequencies` attribute, which
    is a ``list`` of ``tuples``, where each ``tuple`` defines a lower bound,
    upper bound, and degeneracy. Each group also has a `symmetry` correction.
    """

    def __init__(self, index=-1, frequencies=None, symmetry=1):
        self.index = index
        self.frequencies = frequencies or []
        self.symmetry = symmetry
    
    def generateFrequencies(self, count=1):
        """
        Generate a set of frequencies. For each characteristic frequency group,
        the number of frequencies returned is degeneracy * count, and these are
        distributed linearly between the lower and upper bounds.
        """
        frequencies = []
        for lower, upper, degeneracy in self.frequencies:
            number = self.degeneracy * count
            if number == 1:
                frequencies.append((self.lower + self.upper) / 2.0)
            else:
                frequencies.extend(list(numpy.linspace(self.lower, self.upper, number, endpoint=True)))
        return frequencies

################################################################################

class FrequencyDatabase:
    """
    An RMG group frequency database, which associates various functional groups
    with ranges of characteristic vibrational frequencies.
    """

    def __init__(self, path=''):
        if path != '':
            self.load(path)
        else:
            self.groupDatabase = None

    def __loadDatabase(self, dictstr, treestr, libstr):
        """
        Load a thermodynamics group additivity database. The database is stored
        in three files: `dictstr` is the path to the dictionary, `treestr` to
        the tree, and `libstr` to the library. The tree is optional, and should
        be set to '' if not desired.
        """

        # Load dictionary, library, and (optionally) tree
        database = Database()
        database.load(dictstr, treestr, libstr)

        # Convert data in library to ThermoGAModel objects or lists of
        # [link, comment] pairs
        for label, item in database.library.iteritems():
            
            if item is None:
                pass
            elif not item.__class__ is tuple:
                raise InvalidDatabaseError('Thermo library should be tuple at this point. Instead got %r'%data)
            else:
                
                index, items = item

                try:
                    comment = ''
                    items = items.split()
                    database.library[label] = self.__convertLibraryEntry(list(items), comment)
                    database.library[label].index = index

                except (ValueError, IndexError), e:
                    # Split data into link string and comment string; store
                    # as list of length 2
                    link = items[0]
                    comment = item[len(link)+1:].strip()
                    database.library[label] = [link, comment]

        # Check for well-formedness
        if not database.isWellFormed():
            raise InvalidDatabaseError('Database at "%s" is not well-formed.' % (dictstr))

        #database.library.removeLinks()

        return database

    def __convertLibraryEntry(self, freqData, comment):
        """
        Load a group frequency database. The database is stored
        in three files: `dictstr` is the path to the dictionary, `treestr` to
        the tree, and `libstr` to the library. The tree is optional, and should
        be set to '' if not desired.
        """

        gf = GroupFrequency()

        # First item is the symmetry correction
        gf.symmetry = int(freqData.pop(0))

        # Remaining items should be a multiple of three (no comments allowed at the moment)
        if len(freqData) % 3 != 0:
            raise InvalidDatabaseError('Unexpected number of items encountered in group frequencies library.')

        # Convert list of data into a list of characteristic frequencies
        count = len(freqData) / 3
        for i in range(count):
            lower = float(freqData[3*i+1])
            upper = float(freqData[3*i+2])
            degeneracy = int(freqData[3*i])
            gf.frequencies.append((lower, upper, degeneracy))

        return gf

    def __getDTLPaths(self, path, prefix):
        """
        Return a tuple of dictionary, tree, and library paths for a given
        prefix.
        """
        dict_path = os.path.join(path, '%s_Dictionary.txt' % prefix)
        tree_path = os.path.join(path, '%s_Tree.txt' % prefix)
        libr_path = os.path.join(path, '%s_Library.txt' % prefix)
        return dict_path, tree_path, libr_path

    def load(self, path):
        """
        Load a set of thermodynamics group additivity databases from the general
        database specified at `path`.
        """
        path = os.path.abspath(path)
        
        logging.info('Loading functional group frequency database from %s...' % path)
        dict_path = os.path.join(path, 'Dictionary.txt')
        tree_path = os.path.join(path, 'Tree.txt')
        libr_path = os.path.join(path, 'Library.txt')
        self.groupDatabase = self.__loadDatabase(dict_path, tree_path, libr_path)

    def getFrequencies(self, molecule):
        """
        Return the set of characteristic group frequencies corresponding to the
        speficied `molecule`. This is done by searching the molecule for
        certain functional groups for which characteristic frequencies are
        known, and using those frequencies.
        """

        frequencies = []

        # No spectral data for single atoms
        if len(molecule.atoms) > 1:
            
            # For each group in library, find all subgraph isomorphisms
            groupCount = {}
            for node, data in self.groupDatabase.library.iteritems():
                ismatch, mappings = molecule.findSubgraphIsomorphisms(self.groupDatabase.dictionary[node])
                count = len(mappings) if ismatch else 0
                if count % data[0] != 0:
                    raise Exception('Incorrect number of matches of node "%s" while estimating frequencies of %s; expected a multiple of %s, got %s.' % (node, struct, data[0], count))
                groupCount[node] = count / data[0]

            # Get characteristic frequencies
            for node, count in groupCount.iteritems():
                for charFreq in self.groupDatabase.library[node][1:]:
                    frequencies.extend(charFreq.generateFrequencies(count))

        return frequencies

################################################################################

frequencyDatabase = None

def loadFrequencyDatabase(dstr):
    """
    Load the RMG frequency database located at `dstr` into the global variable
    :data:`frequencyDatabase`.
    """
    global frequencyDatabase

    path = os.path.join(dstr,'frequencies_groups')

    # Create and load frequency databases
    frequencyDatabase = FrequencyDatabase()
    logging.debug('\tFrequencies database from '+path)
    frequencyDatabase.load(path)

    return frequencyDatabase

################################################################################

def generateFrequencyData(molecule):
    """
    Use the previously-loaded frequency database to generate a set of
    characteristic group frequencies corresponding to the speficied `molecule`.
    """
    return frequencyDatabase.getFrequencies(molecule)

################################################################################
