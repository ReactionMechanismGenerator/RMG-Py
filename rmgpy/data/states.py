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

import rmgpy.chem.constants as constants
from rmgpy.chem.states import *
from rmgpy.chem.pattern import MoleculePattern

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
            number = degeneracy * count
            if number == 1:
                frequencies.append((lower + upper) / 2.0)
            else:
                frequencies.extend(list(numpy.linspace(lower, upper, number, endpoint=True)))
        return frequencies

################################################################################

class FrequencyGroupDatabase(Database):
    """
    An RMG group frequency database, which associates various functional groups
    with ranges of characteristic vibrational frequencies.
    """

    def __init__(self, path='', old=False):
        Database.__init__(self)
        if path != '':
            self.load(path, old)

    def __loadDatabase(self, dictstr, treestr, libstr):
        """
        Load a thermodynamics group additivity database. The database is stored
        in three files: `dictstr` is the path to the dictionary, `treestr` to
        the tree, and `libstr` to the library. The tree is optional, and should
        be set to '' if not desired.
        """

        # Load dictionary, library, and (optionally) tree
        Database.load(self, dictstr, treestr, libstr)

        # Convert data in library to ThermoGAModel objects or lists of
        # [link, comment] pairs
        for label, item in self.library.iteritems():
            
            if item is None:
                pass
            elif not item.__class__ is tuple:
                raise InvalidDatabaseError('Thermo library should be tuple at this point. Instead got %r'%data)
            else:
                
                index, items = item

                try:
                    comment = ''
                    items = items.split()
                    self.library[label] = self.__convertLibraryEntry(list(items), comment)
                    self.library[label].index = index

                except (ValueError, IndexError), e:
                    # Split data into link string and comment string; store
                    # as list of length 2
                    link = items[0]
                    comment = item[len(link)+1:].strip()
                    self.library[label] = [link, comment]

        # Check for well-formedness
        if not self.isWellFormed():
            raise InvalidDatabaseError('Database at "%s" is not well-formed.' % (dictstr))

        #self.library.removeLinks()

        return self

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

    def load(self, path, old=False):
        """
        Load a set of thermodynamics group additivity databases from the general
        database specified at `path`.
        """
        path = os.path.abspath(path)
        logging.info('Loading functional group frequency database from %s...' % path)

        if old:
            dict_path = os.path.join(path, 'Dictionary.txt')
            tree_path = os.path.join(path, 'Tree.txt')
            libr_path = os.path.join(path, 'Library.txt')
            self.__loadDatabase(dict_path, tree_path, libr_path)

        logging.info('')

    def __countMatchesToNode(self, molecule, node):
        
        count = 0
        if isinstance(self.dictionary[node], MoleculePattern):
            ismatch, mappings = molecule.findSubgraphIsomorphisms(self.dictionary[node])
            count = len(mappings) if ismatch else 0
        elif isinstance(self.dictionary[node], LogicOr):
            for child in self.dictionary[node].components:
                count += self.__countMatchesToNode(molecule, child)
        return count

    def __getNode(self, molecule, atom):
        """
        Determine the group additivity thermodynamic data for the atom `atom`
        in the structure `structure`.
        """

        node0 = self.descendTree(molecule, atom, None)

        if node0 is None:
            raise KeyError('Node not found in database.')

        # It's possible (and allowed) that items in the tree may not be in the
        # library, in which case we need to fall up the tree until we find an
        # ancestor that has an entry in the library
        node = node0
        while node not in self.library and node is not None:
            node = self.tree.parent[node]
        
        return node

    def getFrequencyGroups(self, molecule):
        """
        Return the set of characteristic group frequencies corresponding to the
        speficied `molecule`. This is done by searching the molecule for
        certain functional groups for which characteristic frequencies are
        known, and using those frequencies.
        """

        frequencies = []
        groupCount = {}

        # Generate estimate of thermodynamics
        for atom in molecule.atoms:
            # Iterate over heavy (non-hydrogen) atoms
            if atom.isNonHydrogen():
                node = self.__getNode(molecule, {'*':atom})
                if node is not None:
                    try:
                        groupCount[node] += 1
                    except KeyError:
                        groupCount[node] = 1
                    
        return groupCount

################################################################################

frequencyDatabases = []

def loadFrequencyDatabase(dstr, group=True, old=False):
    """
    Load the RMG frequency database located at `dstr` into the global variable
    :data:`frequencyDatabase`.
    """
    global frequencyDatabases

    if group:
        frequencyDatabase = FrequencyGroupDatabase()
    else:
        frequencyDatabase = None
    frequencyDatabase.load(path=dstr, old=old)
    frequencyDatabases.append(frequencyDatabase)

    return frequencyDatabase

################################################################################

def generateFrequencyData(molecule, thermoModel):
    """
    Use the previously-loaded frequency database to generate a set of
    characteristic group frequencies corresponding to the speficied `molecule`.
    """

    if len(molecule.atoms) < 2:
        return None

    linear = molecule.isLinear()
    numRotors = molecule.countInternalRotors()
    numVibrations = 3 * len(molecule.atoms) - (5 if linear else 6) - numRotors

    frequencyDatabase = frequencyDatabases[-1]

    # Get characteristic frequency groups and the associated frequencies
    groupCount = frequencyDatabase.getFrequencyGroups(molecule)
    frequencies = []
    for node, count in groupCount.iteritems():
        if count != 0: frequencies.extend(frequencyDatabase.library[node].generateFrequencies(count))

    # Check that we have the right number of degrees of freedom specified
    if len(frequencies) > numVibrations:
        # We have too many vibrational modes
        difference = len(frequencies) - numVibrations
        # First try to remove hindered rotor modes until the proper number of modes remain
        if numRotors > difference:
            numRotors -= difference
            numVibrations = len(frequencies)
            logging.warning('For %s, more characteristic frequencies were generated than vibrational modes allowed. Removed %i internal rotors to compensate.' % (molecule, difference))
        # If that won't work, turn off functional groups until the problem is underspecified again
        else:
            groupsRemoved = 0
            freqsRemoved = 0
            freqCount = len(frequencies)
            while freqCount > numVibrations:
                minDegeneracy, minNode = min([(frequencyDatabase.library[node].symmetry, node) for node in groupCount if groupCount[node] > 0])
                if groupCount[minNode] > 1:
                    groupCount[minNode] -= 1
                else:
                    del groupCount[minNode]
                groupsRemoved += 1
                freqsRemoved += minDegeneracy
                freqCount -= minDegeneracy
            # Log warning
            logging.warning('For %s, more characteristic frequencies were generated than vibrational modes allowed. Removed %i groups (%i frequencies) to compensate.' % (molecule, groupsRemoved, freqsRemoved))
            # Regenerate characteristic frequencies
            frequencies = []
            for node, count in groupCount.iteritems():
                if count != 0: frequencies.extend(frequencyDatabase.library[node].generateFrequencies(count))

    # Subtract out contributions to heat capacity from the group frequencies
    Tlist = numpy.arange(300.0, 1501.0, 100.0, numpy.float64)
    Cv = numpy.array([thermoModel.getHeatCapacity(T) / constants.R for T in Tlist], numpy.float64)
    ho = HarmonicOscillator(frequencies=frequencies)
    Cv -= ho.getHeatCapacities(Tlist) / constants.R
    # Subtract out translational modes
    Cv -= 1.5
    # Subtract out external rotational modes
    Cv -= (1.5 if not linear else 1.0)
    # Subtract out PV term (Cp -> Cv)
    Cv -= 1.0
    # Check that all Cv values are still positive
    # We allow a small amount of negative values to allow for the approximate
    # nature of the heat capacity data being used
    assert all([C > -0.05 for C in Cv]), "For %s, reduced Cv data is %s" % (molecule, Cv)

    # Fit remaining frequencies and hindered rotors to the heat capacity data
    from statesfit import fitStatesToHeatCapacity
    modes = fitStatesToHeatCapacity(Tlist, Cv, numVibrations - len(frequencies), numRotors, molecule)
    statesModel = StatesModel(modes=modes)
    
    return statesModel

################################################################################
