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

import rmgpy.constants as constants
from rmgpy.statmech import Conformer, HarmonicOscillator, LinearRotor, NonlinearRotor, HinderedRotor, IdealGasTranslation
from rmgpy.molecule import Molecule, Group, InvalidAdjacencyListError

from base import Database, Entry, LogicNode, LogicOr, LogicAnd, makeLogicNode

################################################################################

def saveEntry(f, entry):
    """
    Write a Pythonic string representation of the given `entry` in the thermo
    database to the file object `f`.
    """
    f.write('entry(\n')
    f.write('    index = {0:d},\n'.format(entry.index))
    f.write('    label = "{0}",\n'.format(entry.label))

    if isinstance(entry.item, Molecule):
        f.write('    molecule = \n')
        f.write('"""\n')
        f.write(entry.item.toAdjacencyList(removeH=False))
        f.write('""",\n')
    elif isinstance(entry.item, Group):
        f.write('    group = \n')
        f.write('"""\n')
        f.write(entry.item.toAdjacencyList())
        f.write('""",\n')
    else:
        f.write('    group = "{0}",\n'.format(entry.item))

    if isinstance(entry.data, GroupFrequencies):
        f.write('    statmech = GroupFrequencies(\n')
        f.write('        frequencies = [\n')
        for lower, upper, degeneracy in entry.data.frequencies:
            f.write('            ({0:g}, {1:g}, {2:d}),\n'.format(lower, upper, degeneracy))
        f.write('        ],\n')
        f.write('        symmetry = {0:d},\n'.format(entry.data.symmetry))
        f.write('    ),\n')
    else:
        f.write('    statmech = {0!r},\n'.format(entry.data))

    if entry.reference is not None: f.write('    reference = {0!r},\n'.format(entry.reference))
    if entry.referenceType != "": f.write('    referenceType = "{0}",\n'.format(entry.referenceType))
    f.write('    shortDesc = u"""')
    f.write(entry.shortDesc)
    f.write('""",\n')
    f.write('    longDesc = \n')
    f.write('u"""\n')
    f.write(entry.longDesc.strip() + "\n")
    f.write('\n""",\n')

    f.write(')\n\n')

def generateOldLibraryEntry(data):
    """
    Return a list of values used to save entries to the old-style RMG
    thermo database based on the thermodynamics object `data`.
    """
    items = '{0:3d}'.format(data.symmetry)
    for lower, upper, degeneracy in data.frequencies:
        items += '     {0:3d} {1:9.1f} {2:9.1f}'.format(degeneracy,lower,upper)
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
        raise ValueError('format parameter must be "library" or "groups"; got "{0}" instead'.format(format))

################################################################################

class StatmechDepository(Database):
    """
    A class for working with the RMG statistical mechanics (frequencies)
    depository.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self,
                  index,
                  label,
                  molecule,
                  statmech,
                  reference=None,
                  referenceType='',
                  shortDesc='',
                  longDesc='',
                  ):
        
        item = Molecule().fromAdjacencyList(molecule)
        
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = item,
            data = statmech,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
        )

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return saveEntry(f, entry)

################################################################################

class StatmechLibrary(Database):
    """
    A class for working with a RMG statistical mechanics (frequencies) library.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self,
                  index,
                  label,
                  molecule,
                  statmech,
                  reference=None,
                  referenceType='',
                  shortDesc='',
                  longDesc='',
                  ):
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = Molecule().fromAdjacencyList(molecule),
            data = statmech,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
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

class StatmechGroups(Database):
    """
    A class for working with an RMG statistical mechanics (frequencies) group
    database.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self,
                  index,
                  label,
                  group,
                  statmech,
                  reference=None,
                  referenceType='',
                  shortDesc='',
                  longDesc='',
                  ):
        if ( group[0:3].upper() == 'OR{' or
             group[0:4].upper() == 'AND{' or
             group[0:7].upper() == 'NOT OR{' or
             group[0:8].upper() == 'NOT AND{'
            ):
            item = makeLogicNode(group)
        else:
            item = Group().fromAdjacencyList(group)
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = item,
            data = statmech,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
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
        statmech database, returning the corresponding thermodynamics object.
        """
        return processOldLibraryEntry(data, "groups")

    def __countMatchesToNode(self, molecule, node):
        """
        Count the number of matches in the given :class:`Molecule` object
        `molecule` for the given `node` in the group database.
        """
        count = 0
        if isinstance(self.dictionary[node], Group):
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
        if node.data is None:
            logging.warning('Statmech node {0!r} and all its parents have data=None'.format(node0))
            return None
            raise KeyError('Statmech node {0!r} and all its parents have data=None'.format(node0))
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

        # This is an additional hardcoded functional group for C-H with C in a ring
        # It is hardcoded because the adjacency list format isn't well-equipped
        # to handle this sort of functional group
        ringCH = Entry(
            label = 'ringCH',
            item = None,
            data = GroupFrequencies([(2750., 3150., 1), (900., 1100., 1)]),
        )
        
        # Generate estimate of thermodynamics
        for atom in molecule.atoms:
            # Iterate over heavy (non-hydrogen) atoms
            if atom.isHydrogen(): continue
            if molecule.isAtomInCycle(atom):
                # Atom is in cycle
                # Add each C-H bond to the ringCH group
                # This is hardcoding of functional groups!
                if atom.isCarbon():
                    for atom2 in atom.edges:
                        if atom2.isHydrogen():
                            try:
                                groupCount[ringCH] += 1
                            except KeyError:
                                groupCount[ringCH] = 1
            else:
                # Atom is not in cycle, so find a group for it
                node = self.__getNode(molecule, {'*':atom})
                if node is not None:
                    try:
                        groupCount[node] += 1
                    except KeyError:
                        groupCount[node] = 1

        return groupCount

    def getStatmechData(self, molecule, thermoModel):
        """
        Use the previously-loaded frequency database to generate a set of
        characteristic group frequencies corresponding to the speficied
        `molecule`. The provided thermo data in `thermoModel` is used to fit
        some frequencies and all hindered rotors to heat capacity data.
        """
        conformer = Conformer()
        
        # Compute spin multiplicity
        # For closed-shell molecule the spin multiplicity is 1
        # For monoradicals the spin multiplicity is 2
        # For higher-order radicals the highest allowed spin multiplicity is assumed
        conformer.spinMultiplicity = molecule.getRadicalCount() + 1
        
        # No need to determine rotational and vibrational modes for single atoms
        if len(molecule.atoms) < 2:
            return (conformer, None, None)

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
                logging.warning('For {0}, more characteristic frequencies were generated than vibrational modes allowed. Removed {1:d} internal rotors to compensate.'.format(molecule, difference))
            # If that won't work, turn off functional groups until the problem is underspecified again
            else:
                groupsRemoved = 0
                freqsRemoved = 0
                freqCount = len(frequencies)
                while freqCount > numVibrations:
                    minDegeneracy, minEntry = min([(entry.data.symmetry, entry) for entry in groupCount if groupCount[entry] > 0])
                    if groupCount[minEntry] > 1:
                        groupCount[minEntry] -= 1
                    else:
                        del groupCount[minEntry]
                    groupsRemoved += 1
                    freqsRemoved += minDegeneracy
                    freqCount -= minDegeneracy
                # Log warning
                logging.warning('For {0}, more characteristic frequencies were generated than vibrational modes allowed. Removed {1:d} groups ({2:d} frequencies) to compensate.'.format(molecule, groupsRemoved, freqsRemoved))
                # Regenerate characteristic frequencies
                frequencies = []
                for entry, count in groupCount.iteritems():
                    if count != 0: frequencies.extend(entry.data.generateFrequencies(count))

        # Subtract out contributions to heat capacity from the group frequencies
        Tlist = numpy.arange(300.0, 1501.0, 100.0, numpy.float64)
        Cv = numpy.array([thermoModel.getHeatCapacity(T) / constants.R for T in Tlist], numpy.float64)
        ho = HarmonicOscillator(frequencies=(frequencies,"cm^-1"))
        for i in range(Tlist.shape[0]):
            Cv[i] -= ho.getHeatCapacity(Tlist[i]) / constants.R
        # Subtract out translational modes
        Cv -= 1.5
        # Subtract out external rotational modes
        Cv -= (1.5 if not linear else 1.0)
        # Subtract out PV term (Cp -> Cv)
        Cv -= 1.0
        
        # Fit remaining frequencies and hindered rotors to the heat capacity data
        from statmechfit import fitStatmechToHeatCapacity
        modes = fitStatmechToHeatCapacity(Tlist, Cv, numVibrations - len(frequencies), numRotors, molecule)
        for mode in modes:
            if isinstance(mode, HarmonicOscillator):
                uncertainties = [0 for f in frequencies] # probably shouldn't be zero
                frequencies.extend(mode.frequencies.value_si)
                uncertainties.extend(mode.frequencies.uncertainty)
                mode.frequencies.value_si = numpy.array(frequencies, numpy.float)
                mode.frequencies.uncertainty = numpy.array(uncertainties, numpy.float)
                break
        else:
            modes.insert(0, HarmonicOscillator(frequencies=(frequencies,"cm^-1")))

        conformer.modes = modes

        return (conformer, None, None)

################################################################################

class StatmechDatabase(object):
    """
    A class for working with the RMG statistical mechanics (frequencies) database.
    """

    def __init__(self):
        self.depository = {}
        self.libraries = {}
        self.groups = {}
        self.libraryOrder = []
        self.local_context = {
            'HarmonicOscillator': HarmonicOscillator,
            'LinearRotor': LinearRotor,
            'NonlinearRotor': NonlinearRotor,
            'HinderedRotor': HinderedRotor,
            'IdealGasTranslation': IdealGasTranslation,
            'GroupFrequencies': GroupFrequencies,
        }
        self.global_context = {}

    def __reduce__(self):
        """
        A helper function used when pickling a StatmechDatabase object.
        """
        d = {
            'depository': self.depository,
            'libraries': self.libraries,
            'groups': self.groups,
            'libraryOrder': self.libraryOrder,
        }
        return (StatmechDatabase, (), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling a StatmechDatabase object.
        """
        self.depository = d['depository']
        self.libraries = d['libraries']
        self.groups = d['groups']
        self.libraryOrder = d['libraryOrder']

    def load(self, path, libraries=None, depository=True):
        """
        Load the statmech database from the given `path` on disk, where `path`
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
        Load the statmech database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        self.depository = {}
        self.depository['depository']  = StatmechDepository().load(os.path.join(path, 'depository.py'), self.local_context, self.global_context)

    def loadLibraries(self, path, libraries=None):
        """
        Load the statmech database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        self.libraries = {}; self.libraryOrder = []
        for (root, dirs, files) in os.walk(os.path.join(path)):
            for f in files:
                name, ext = os.path.splitext(f)
                if ext.lower() == '.py' and (libraries is None or name in libraries):
                    logging.info('Loading frequencies library from {0} in {1}...'.format(f, root))
                    library = StatmechLibrary()
                    library.load(os.path.join(root, f), self.local_context, self.global_context)
                    library.label = os.path.splitext(f)[0]
                    self.libraries[library.label] = library
                    self.libraryOrder.append(library.label)
        if libraries is not None:
            self.libraryOrder = libraries

    def loadGroups(self, path):
        """
        Load the statmech database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        logging.info('Loading frequencies group database from {0}...'.format(path))
        self.groups = {}
        self.groups['groups'] = StatmechGroups().load(os.path.join(path, 'groups.py' ), self.local_context, self.global_context)

    def save(self, path):
        """
        Save the statmech database to the given `path` on disk, where `path`
        points to the top-level folder of the statmech database.
        """
        path = os.path.abspath(path)
        if not os.path.exists(path): os.mkdir(path)
        self.saveDepository(os.path.join(path, 'depository'))
        self.saveLibraries(os.path.join(path, 'libraries'))
        self.saveGroups(os.path.join(path, 'groups'))

    def saveDepository(self, path):
        """
        Save the statmech depository to the given `path` on disk, where `path`
        points to the top-level folder of the statmech depository.
        """
        if not os.path.exists(path): os.mkdir(path)
        for name, depository in self.depository.iteritems():
            depository.save(os.path.join(path, name + '.py'))

    def saveLibraries(self, path):
        """
        Save the statmech libraries to the given `path` on disk, where `path`
        points to the top-level folder of the statmech libraries.
        """
        if not os.path.exists(path): os.mkdir(path)
        for library in self.libraries.values():
            library.save(os.path.join(path, '{0}.py'.format(library.label)))

    def saveGroups(self, path):
        """
        Save the statmech groups to the given `path` on disk, where `path`
        points to the top-level folder of the statmech groups.
        """
        if not os.path.exists(path): os.mkdir(path)
        for name, groups in self.groups.iteritems():
            groups.save(os.path.join(path, name + '.py'))

    def loadOld(self, path):
        """
        Load the old RMG thermo database from the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        # The old database does not have a depository, so create an empty one
        self.depository = {}
        self.depository['depository'] = StatmechDepository(label='depository', name='Statmech Depository')

        for (root, dirs, files) in os.walk(os.path.join(path, 'frequencies_libraries')):
            if os.path.exists(os.path.join(root, 'Dictionary.txt')) and os.path.exists(os.path.join(root, 'Library.txt')):
                library = StatmechLibrary(label=os.path.basename(root), name=os.path.basename(root))
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

        self.groups = StatmechGroups(label='group', name='Functional Group Values').loadOld(
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
        for library in self.libraries.values():
            if not os.path.exists(librariesPath): os.mkdir(librariesPath)
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

    def getStatmechData(self, molecule, thermoModel=None):
        """
        Return the thermodynamic parameters for a given :class:`Molecule`
        object `molecule`. This function first searches the loaded libraries
        in order, returning the first match found, before falling back to
        estimation via group additivity.
        """
        statmechModel = None
        # Check the libraries in order first; return the first successful match
        for label in self.libraryOrder:
            statmechModel = self.getStatmechDataFromLibrary(molecule, self.libraries[label])
            if statmechModel: break
        else:
            # Thermo not found in any loaded libraries, so estimate
            statmechModel = self.getStatmechDataFromGroups(molecule, thermoModel)
        return statmechModel[0]

    def getStatmechDataFromDepository(self, molecule):
        """
        Return statmech data for the given :class:`Molecule` object `molecule`
        by searching the entries in the depository.
        Returns a list of tuples  (statmechData, depository, entry).
        """
        items = []
        for name, depository in self.depository.iteritems():
            for label, entry in depository.entries.iteritems():
                if molecule.isIsomorphic(entry.item):
                    items.append((entry.data, self.depository[name], entry))
        return items

    def getStatmechDataFromLibrary(self, molecule, library):
        """
        Return statmech data for the given :class:`Molecule` object `molecule`
        by searching the entries in the specified :class:`StatmechLibrary` object
        `library`. Returns ``None`` if no data was found.
        """
        for label, entry in library.entries.iteritems():
            if molecule.isIsomorphic(entry.item):
                return (entry.data, library, entry)
        return None

    def getStatmechDataFromGroups(self, molecule, thermoModel):
        """
        Return statmech data for the given :class:`Molecule` object `molecule`
        by estimating using characteristic group frequencies and fitting the
        remaining internal modes to heat capacity data from the given thermo
        model `thermoModel`. This always returns valid degrees of freedom data.
        """
        return self.groups['groups'].getStatmechData(molecule, thermoModel)
        
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
