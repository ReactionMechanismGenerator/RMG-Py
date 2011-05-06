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

def saveEntry(f, entry):
    """
    Write a Pythonic string representation of the given `entry` in the thermo
    database to the file object `f`.
    """
    f.write('entry(\n')
    f.write('    index = %i,\n' % (entry.index))
    f.write('    label = "%s",\n' % (entry.label))

    if isinstance(entry.item, Molecule):
        f.write('    molecule = \n')
        f.write('"""\n')
        f.write(entry.item.toAdjacencyList(removeH=True))
        f.write('""",\n')
    elif isinstance(entry.item, MoleculePattern):
        f.write('    group = \n')
        f.write('"""\n')
        f.write(entry.item.toAdjacencyList())
        f.write('""",\n')
    else:
        f.write('    group = "%s",\n' % (entry.item))

    if isinstance(entry.data, StatesModel):
        f.write('    states = States(\n')
        f.write('        modes = [\n')
        for mode in entry.data.modes:
            f.write('            %r,\n' % mode)
        f.write('        ],\n')
        f.write('        spinMultiplicity = %r,\n' % (entry.data.spinMultiplicity))
        f.write('    ),\n')
    elif isinstance(entry.data, GroupFrequencies):
        f.write('    states = GroupFrequencies(\n')
        f.write('        frequencies = [\n')
        for lower, upper, degeneracy in entry.data.frequencies:
            f.write('            (%g, %g, %d),\n' % (lower, upper, degeneracy))
        f.write('        ],\n')
        f.write('        symmetry = %d,\n' % (entry.data.symmetry))
        f.write('    ),\n')
    else:
        f.write('        data = %r,\n' % (entry.data))

    if entry.reference is not None: f.write('    reference = %r,\n' % (entry.reference))
    if entry.referenceType != "": f.write('    referenceType = "%s",\n' % (entry.referenceType))
    f.write('    shortDesc = """%s""",\n' % (entry.shortDesc))
    f.write('    longDesc = \n')
    f.write('"""\n')
    f.write(entry.longDesc)
    f.write('\n""",\n')

    f.write('    history = [\n')
    for time, user, action, description in entry.history:
        f.write('        ("%s","%s","%s","""%s"""),\n' % (time, user, action, description))
    f.write('    ],\n')

    f.write(')\n\n')

def generateOldLibraryEntry(data):
    """
    Return a list of values used to save entries to the old-style RMG
    thermo database based on the thermodynamics object `data`.
    """
    items = []
    if isinstance(data, StatesModel):
        pass
    elif isinstance(data, list):
        for mode in data:
            items.extend([mode.degeneracy, mode.lower, mode.upper])
    else:
        raise ValueError('data parameter must be in StatesModel format or a list; got %s instead' % (data.__class__))
    return items

def processOldLibraryEntry(data, format):
    """
    Process a list of parameters `data` as read from an old-style RMG
    thermo database, returning the corresponding thermodynamics object.
    """
    if format == 'library':
        pass
    elif format == 'groups':
        frequencies = []
        for degeneracy, lower, upper in zip(data[1::3], data[2::3], data[3::3]):
            frequencies.append((float(lower), float(upper), int(degeneracy)))
        return GroupFrequencies(frequencies=frequencies, symmetry=int(data[0]))
    else:
        raise ValueError('format parameter must be "library" or "groups"; got "%s" instead' % (format))

################################################################################

class StatesDepository(Database):
    """
    A class for working with the RMG states (frequencies) depository.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self, index, label, molecule, states, reference=None, referenceType='', shortDesc='', longDesc='', history=None):
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = Molecule().fromAdjacencyList(molecule),
            data = states,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            history = history or [],
        )

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return saveEntry(f, entry)

################################################################################

class StatesLibrary(Database):
    """
    A class for working with a RMG states (frequencies) library.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self, index, label, molecule, states, reference=None, referenceType='', shortDesc='', longDesc='', history=None):
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = Molecule().fromAdjacencyList(molecule),
            data = states,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            history = history or [],
        )

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return saveEntry(f, entry)

    def generateOldLibraryEntry(self, data):
        """
        Return a list of values used to save entries to the old-style RMG
        thermo database based on the thermodynamics object `data`.
        """
        return generateOldLibraryEntry(data)

    def processOldLibraryEntry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        thermo database, returning the corresponding thermodynamics object.
        """
        return processOldLibraryEntry(data, "library")

################################################################################

class StatesGroups(Database):
    """
    A class for working with an RMG states (frequencies) group database.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self, index, label, group, states, reference=None, referenceType='', shortDesc='', longDesc='', history=None):
        if group[0:3].upper() == 'OR{' or group[0:4].upper() == 'AND{' or group[0:7].upper() == 'NOT OR{' or group[0:8].upper() == 'NOT AND{':
            item = makeLogicNode(group)
        else:
            item = MoleculePattern().fromAdjacencyList(group)
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = item,
            data = states,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            history = history or [],
        )

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return saveEntry(f, entry)

    def generateOldLibraryEntry(self, data):
        """
        Return a list of values used to save entries to the old-style RMG
        thermo database based on the thermodynamics object `data`.
        """
        return generateOldLibraryEntry(data)

    def processOldLibraryEntry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        states database, returning the corresponding thermodynamics object.
        """
        return processOldLibraryEntry(data, "groups")

    def __countMatchesToNode(self, molecule, node):
        """
        Count the number of matches in the given :class:`Molecule` object
        `molecule` for the given `node` in the group database.
        """
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
        For a given :class:`Molecule` object `molecule` with central `atom`,
        determine the most specific functional group that describes that atom
        center and has characteristic frequencies associated with it.
        """

        node0 = self.descendTree(molecule, atom, None)

        if node0 is None:
            raise KeyError('Node not found in database.')

        # It's possible (and allowed) that items in the tree may not be in the
        # library, in which case we need to fall up the tree until we find an
        # ancestor that has an entry in the library
        node = node0
        while node.data is None and node.parent is not None:
            node = node.parent

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

    def getStatesData(self, molecule, thermoModel):
        """
        Use the previously-loaded frequency database to generate a set of
        characteristic group frequencies corresponding to the speficied
        `molecule`. The provided thermo data in `thermoModel` is used to fit
        some frequencies and all hindered rotors to heat capacity data.
        """

        # No need to determine rotational and vibrational modes for single atoms
        if len(molecule.atoms) < 2:
            return states

        linear = molecule.isLinear()
        numRotors = molecule.countInternalRotors()
        numVibrations = 3 * len(molecule.atoms) - (5 if linear else 6) - numRotors

        # Get characteristic frequency groups and the associated frequencies
        groupCount = self.getFrequencyGroups(molecule)
        frequencies = []
        for entry, count in groupCount.iteritems():
            if count != 0 and entry.data is not None: frequencies.extend(entry.data.generateFrequencies(count))

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
        
        # Fit remaining frequencies and hindered rotors to the heat capacity data
        from statesfit import fitStatesToHeatCapacity
        modes = fitStatesToHeatCapacity(Tlist, Cv, numVibrations - len(frequencies), numRotors, molecule)
        for mode in modes:
            if isinstance(mode, HarmonicOscillator):
                frequencies.extend(mode.frequencies)
                mode.frequencies = frequencies
                break
        else:
            modes.insert(0, HarmonicOscillator(frequencies=frequencies))

        statesModel = StatesModel(modes=modes)

        return (statesModel, None, None)

################################################################################

class StatesDatabase:
    """
    A class for working with the RMG states (frequencies) database.
    """

    def __init__(self):
        self.depository = {}
        self.libraries = {}
        self.groups = {}
        self.libraryOrder = []
        self.local_context = {
            'HarmonicOscillator': HarmonicOscillator,
            'RigidRotor': RigidRotor,
            'HinderedRotor': HinderedRotor,
            'Translation': Translation,
            'GroupFrequencies': GroupFrequencies,
            'States': StatesModel,
        }
        self.global_context = {}

    def load(self, path, libraries=None, depository=True):
        """
        Load the states database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        if depository:
            self.loadDepository(os.path.join(path, 'depository'))
        else:
            self.depository = {}
        self.loadLibraries(os.path.join(path, 'libraries'), libraries)
        self.loadGroups(os.path.join(path, 'groups'))

    def loadDepository(self, path):
        """
        Load the states database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        self.depository = StatesDepository().load(os.path.join(path, 'depository.py'), self.local_context, self.global_context)

    def loadLibraries(self, path, libraries=None):
        """
        Load the states database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        self.libraries = {}; self.libraryOrder = []
        for (root, dirs, files) in os.walk(os.path.join(path)):
            for f in files:
                name, ext = os.path.splitext(f)
                if ext.lower() == '.py' and (libraries is None or name in libraries):
                    logging.info('Loading frequencies library from %s in %s...' % (f, root))
                    library = StatesLibrary()
                    library.load(os.path.join(root, f), self.local_context, self.global_context)
                    library.label = os.path.splitext(f)[0]
                    self.libraries[library.label] = library
                    self.libraryOrder.append(library.label)
        if libraries is not None:
            self.libraryOrder = libraries

    def loadGroups(self, path):
        """
        Load the states database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        logging.info('Loading frequencies group database from %s...' % (path))
        self.groups = StatesGroups().load(os.path.join(path, 'groups.py' ), self.local_context, self.global_context)

    def save(self, path):
        """
        Save the states database to the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        path = os.path.abspath(path)
        if not os.path.exists(path): os.mkdir(path)

        depositoryPath = os.path.join(path, 'depository')
        if not os.path.exists(depositoryPath): os.mkdir(depositoryPath)
        self.depository.save(os.path.join(depositoryPath, 'depository.py'))

        librariesPath = os.path.join(path, 'libraries')
        if not os.path.exists(librariesPath): os.mkdir(librariesPath)
        for library in self.libraries.values():
            library.save(os.path.join(librariesPath, '%s.py' % (library.label)))

        groupsPath = os.path.join(path, 'groups')
        if not os.path.exists(groupsPath): os.mkdir(groupsPath)
        self.groups.save(os.path.join(groupsPath, 'groups.py'))

    def loadOld(self, path):
        """
        Load the old RMG thermo database from the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        # The old database does not have a depository, so create an empty one
        self.depository = StatesDepository(label='depository', name='States Depository')

        for (root, dirs, files) in os.walk(os.path.join(path, 'frequencies_libraries')):
            if os.path.exists(os.path.join(root, 'Dictionary.txt')) and os.path.exists(os.path.join(root, 'Library.txt')):
                library = StatesLibrary(label=os.path.basename(root), name=os.path.basename(root))
                library.loadOld(
                    dictstr = os.path.join(root, 'Dictionary.txt'),
                    treestr = '',
                    libstr = os.path.join(root, 'Library.txt'),
                    numParameters = -1,
                    numLabels = 1,
                    pattern = False,
                )
                library.label = os.path.basename(root)
                self.libraries[library.label] = library

        self.groups = StatesGroups(label='group', name='Functional Group Values').loadOld(
            dictstr = os.path.join(path, 'frequencies_groups', 'Dictionary.txt'),
            treestr = os.path.join(path, 'frequencies_groups', 'Tree.txt'),
            libstr = os.path.join(path, 'frequencies_groups', 'Library.txt'),
            numParameters = -1,
            numLabels = 1,
            pattern = True,
        )

    def saveOld(self, path):
        """
        Save the old RMG thermo database to the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """

        # Depository not used in old database, so it is not saved

        librariesPath = os.path.join(path, 'frequencies_libraries')
        if not os.path.exists(librariesPath): os.mkdir(librariesPath)
        for library in self.libraries.values():
            libraryPath = os.path.join(librariesPath, library.label)
            if not os.path.exists(libraryPath): os.mkdir(libraryPath)
            library.saveOld(
                dictstr = os.path.join(libraryPath, 'Dictionary.txt'),
                treestr = '',
                libstr = os.path.join(libraryPath, 'Library.txt'),
            )

        groupsPath = os.path.join(path, 'frequencies_groups')
        if not os.path.exists(groupsPath): os.mkdir(groupsPath)
        self.groups.saveOld(
            dictstr = os.path.join(groupsPath, 'Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Tree.txt'),
            libstr = os.path.join(groupsPath, 'Library.txt'),
        )

    def getStatesData(self, molecule, thermoModel=None):
        """
        Return the thermodynamic parameters for a given :class:`Molecule`
        object `molecule`. This function first searches the loaded libraries
        in order, returning the first match found, before falling back to
        estimation via group additivity.
        """
        statesModel = None
        # Check the libraries in order first; return the first successful match
        for label in self.libraryOrder:
            statesModel = self.getStatesDataFromLibrary(molecule, self.libraries[label])
            if statesModel: break
        else:
            # Thermo not found in any loaded libraries, so estimate
            statesModel = self.getStatesDataFromGroups(molecule, thermoModel)
        return statesModel[0]

    def getStatesDataFromDepository(self, molecule):
        """
        Return states data for the given :class:`Molecule` object `molecule`
        by searching the entries in the depository.
        """
        items = []
        for label, entry in self.depository.entries.iteritems():
            if molecule.isIsomorphic(entry.item):
                items.append((entry.data, self.depository['stable'], entry))
        return items

    def getStatesDataFromLibrary(self, molecule, library):
        """
        Return states data for the given :class:`Molecule` object `molecule`
        by searching the entries in the specified :class:`StatesLibrary` object
        `library`. Returns ``None`` if no data was found.
        """
        for label, entry in library.entries.iteritems():
            if molecule.isIsomorphic(entry.item):
                return (entry.data, library, entry)
        return None

    def getStatesDataFromGroups(self, molecule, thermoModel):
        """
        Return states data for the given :class:`Molecule` object `molecule`
        by estimating using characteristic group frequencies and fitting the
        remaining internal modes to heat capacity data from the given thermo
        model `thermoModel`. This always returns valid degrees of freedom data.
        """
        return self.groups.getStatesData(molecule, thermoModel)
        
################################################################################

class GroupFrequencies:
    """
    Represent a set of characteristic frequencies for a group in the frequency
    database. These frequencies are stored in the `frequencies` attribute, which
    is a ``list`` of ``tuples``, where each ``tuple`` defines a lower bound,
    upper bound, and degeneracy. Each group also has a `symmetry` correction.
    """

    def __init__(self, frequencies=None, symmetry=1):
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

        # Convert data in library to ThermoData objects or lists of
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
