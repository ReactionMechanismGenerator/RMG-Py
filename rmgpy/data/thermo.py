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

"""

import os
import os.path
import math
import logging
import numpy
from copy import copy, deepcopy

from base import Database, Entry, makeLogicNode

from rmgpy.quantity import constants
from rmgpy.thermo import *
from rmgpy.molecule import Molecule, Atom, Bond
from rmgpy.group import Group

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
        f.write(entry.item.toAdjacencyList(removeH=True))
        f.write('""",\n')
    elif isinstance(entry.item, Group):
        f.write('    group = \n')
        f.write('"""\n')
        f.write(entry.item.toAdjacencyList())
        f.write('""",\n')
    else:
        f.write('    group = "{0}",\n'.format(entry.item))

    if isinstance(entry.data, ThermoData):
        f.write('    thermo = ThermoData(\n')
        f.write('        Tdata = {0!r},\n'.format(entry.data.Tdata))
        f.write('        Cpdata = {0!r},\n'.format(entry.data.Cpdata))
        f.write('        H298 = {0!r},\n'.format(entry.data.H298))
        f.write('        S298 = {0!r},\n'.format(entry.data.S298))
        if entry.data.Tmin is not None: f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None: f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        f.write('    ),\n')
    elif isinstance(entry.data, Wilhoit):
        f.write('    thermo = Wilhoit(\n')
        f.write('        cp0 = {0!r},\n'.format(entry.data.cp0))
        f.write('        cpInf = {0!r},\n'.format(entry.data.cpInf))
        f.write('        a0 = {0:g},\n'.format(entry.data.a0))
        f.write('        a1 = {0:g},\n'.format(entry.data.a1))
        f.write('        a2 = {0:g},\n'.format(entry.data.a2))
        f.write('        a3 = {0:g},\n'.format(entry.data.a3))
        f.write('        B = {0!r},\n'.format(entry.data.B))
        f.write('        H0 = {0!r},\n'.format(entry.data.H0))
        f.write('        S0 = {0!r},\n'.format(entry.data.S0))
        if entry.data.Tmin is not None: f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None: f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        f.write('    ),\n')
    elif isinstance(entry.data, MultiNASA):
        f.write('    thermo = MultiNASA(\n')
        f.write('        polynomials = [\n')
        for poly in entry.data.polynomials:
            f.write('            {0!r},\n'.format(poly))
        f.write('        ],\n')
        if entry.data.Tmin is not None: f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None: f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        f.write('    ),\n')
    else:
        f.write('    thermo = {0!r},\n'.format(entry.data))

    if entry.reference is not None: f.write('    reference = {0!r},\n'.format(entry.reference))
    if entry.referenceType != "": f.write('    referenceType = "{0}",\n'.format(entry.referenceType))
    f.write('    shortDesc = u"""{0}""",\n'.format(entry.shortDesc))
    f.write('    longDesc = \n')
    f.write('u"""\n')
    f.write(entry.longDesc.strip() + "\n")
    f.write('""",\n')

    f.write('    history = [\n')
    for time, user, action, description in entry.history:
        f.write('        ("{0}","{1}","{2}","""{3}"""),\n'.format(time, user, action, description))
    f.write('    ],\n')

    f.write(')\n\n')

def generateOldLibraryEntry(data):
    """
    Return a list of values used to save entries to the old-style RMG
    thermo database based on the thermodynamics object `data`.
    """
    if isinstance(data, ThermoData):
        return '{0:9g} {1:9g} {2:9g} {3:9g} {4:9g} {5:9g} {6:9g} {7:9g} {8:9g} {9:9g} {10:9g} {11:9g}'.format(
            data.H298.value/4184.,
            data.S298.value/4.184,
            data.Cpdata.values[0]/4.184,
            data.Cpdata.values[1]/4.184,
            data.Cpdata.values[2]/4.184,
            data.Cpdata.values[3]/4.184,
            data.Cpdata.values[4]/4.184,
            data.Cpdata.values[5]/4.184,
            data.Cpdata.values[6]/4.184,
            data.H298.uncertainty/4184.,
            data.S298.uncertainty/4.184,
            data.Cpdata.uncertainty/4.184,
        )
    else:
        return '{0:9g} {1:9g} {2:9g} {3:9g} {4:9g} {5:9g} {6:9g} {7:9g} {8:9g} {9:9g} {10:9g} {11:9g}'.format(
            data.getEnthalpy(298)/4184.,
            data.getEntropy(298)/4.184,
            data.getHeatCapacity(300)/4.184,
            data.getHeatCapacity(400)/4.184,
            data.getHeatCapacity(500)/4.184,
            data.getHeatCapacity(600)/4.184,
            data.getHeatCapacity(800)/4.184,
            data.getHeatCapacity(1000)/4.184,
            data.getHeatCapacity(1500)/4.184,
            0,
            0,
            0,
        )

def processOldLibraryEntry(data):
    """
    Process a list of parameters `data` as read from an old-style RMG
    thermo database, returning the corresponding thermodynamics object.
    """
    return ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],"K"),
        Cpdata = ([float(d) for d in data[2:9]],"cal/(mol*K)","+|-",float(data[11])),
        H298 = (float(data[0]),"kcal/mol","+|-",float(data[9])),
        S298 = (float(data[1]),"cal/(mol*K)","+|-",float(data[10])),
    )


################################################################################

class ThermoDepository(Database):
    """
    A class for working with the RMG thermodynamics depository.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self, index, label, molecule, thermo, reference=None, referenceType='', shortDesc='', longDesc='', history=None):
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = Molecule().fromAdjacencyList(molecule),
            data = thermo,
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

class ThermoLibrary(Database):
    """
    A class for working with a RMG thermodynamics library.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self, index, label, molecule, thermo, reference=None, referenceType='', shortDesc='', longDesc='', history=None):
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = Molecule().fromAdjacencyList(molecule),
            data = thermo,
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
        return processOldLibraryEntry(data)

################################################################################

class ThermoGroups(Database):
    """
    A class for working with an RMG thermodynamics group additivity database.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self, index, label, group, thermo, reference=None, referenceType='', shortDesc='', longDesc='', history=None):
        if group[0:3].upper() == 'OR{' or group[0:4].upper() == 'AND{' or group[0:7].upper() == 'NOT OR{' or group[0:8].upper() == 'NOT AND{':
            item = makeLogicNode(group)
        else:
            item = Group().fromAdjacencyList(group)
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = item,
            data = thermo,
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
        return processOldLibraryEntry(data)

################################################################################

class ThermoDatabase:
    """
    A class for working with the RMG thermodynamics database.
    """

    def __init__(self):
        self.depository = {}
        self.libraries = {}
        self.groups = {}
        self.libraryOrder = []
        self.local_context = {
            'ThermoData': ThermoData,
            'Wilhoit': Wilhoit,
            'NASA': NASA,
            'MultiNASA': MultiNASA,
        }
        self.global_context = {}

    def load(self, path, libraries=None, depository=True):
        """
        Load the thermo database from the given `path` on disk, where `path`
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
        Load the thermo database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        self.depository = {}
        self.depository['stable']  = ThermoDepository().load(os.path.join(path, 'stable.py'), self.local_context, self.global_context)
        self.depository['radical'] = ThermoDepository().load(os.path.join(path, 'radical.py'), self.local_context, self.global_context)

    def loadLibraries(self, path, libraries=None):
        """
        Load the thermo database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        self.libraries = {}; self.libraryOrder = []
        for (root, dirs, files) in os.walk(os.path.join(path)):
            for f in files:
                name, ext = os.path.splitext(f)
                if ext.lower() == '.py' and (libraries is None or name in libraries):
                    logging.info('Loading thermodynamics library from {0} in {1}...'.format(f, root))
                    library = ThermoLibrary()
                    library.load(os.path.join(root, f), self.local_context, self.global_context)
                    library.label = os.path.splitext(f)[0]
                    self.libraries[library.label] = library
                    self.libraryOrder.append(library.label)
        if libraries is not None:
            self.libraryOrder = libraries

    def loadGroups(self, path):
        """
        Load the thermo database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        logging.info('Loading thermodynamics group database from {0}...'.format(path))
        self.groups = {}
        self.groups['group']   =   ThermoGroups(label='group').load(os.path.join(path, 'group.py'  ), self.local_context, self.global_context)
        self.groups['gauche']  =  ThermoGroups(label='gauche').load(os.path.join(path, 'gauche.py' ), self.local_context, self.global_context)
        self.groups['int15']   =   ThermoGroups(label='int15').load(os.path.join(path, 'int15.py'  ), self.local_context, self.global_context)
        self.groups['ring']    =    ThermoGroups(label='ring').load(os.path.join(path, 'ring.py'   ), self.local_context, self.global_context)
        self.groups['radical'] = ThermoGroups(label='radical').load(os.path.join(path, 'radical.py'), self.local_context, self.global_context)
        self.groups['other']   =   ThermoGroups(label='other').load(os.path.join(path, 'other.py'  ), self.local_context, self.global_context)

    def save(self, path):
        """
        Save the thermo database to the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        path = os.path.abspath(path)
        if not os.path.exists(path): os.mkdir(path)
        self.saveDepository(os.path.join(path, 'depository'))
        self.saveLibraries(os.path.join(path, 'libraries'))
        self.saveGroups(os.path.join(path, 'groups'))

    def saveDepository(self, path):
        """
        Save the thermo depository to the given `path` on disk, where `path`
        points to the top-level folder of the thermo depository.
        """
        if not os.path.exists(path): os.mkdir(path)
        self.depository['stable'].save(os.path.join(path, 'stable.py'))
        self.depository['radical'].save(os.path.join(path, 'radical.py'))

    def saveLibraries(self, path):
        """
        Save the thermo libraries to the given `path` on disk, where `path`
        points to the top-level folder of the thermo libraries.
        """
        if not os.path.exists(path): os.mkdir(path)
        for library in self.libraries.values():
            library.save(os.path.join(path, '{0}.py'.format(library.label)))

    def saveGroups(self, path):
        """
        Save the thermo groups to the given `path` on disk, where `path`
        points to the top-level folder of the thermo groups.
        """
        if not os.path.exists(path): os.mkdir(path)
        self.groups['group'].save(os.path.join(path, 'group.py'))
        self.groups['gauche'].save(os.path.join(path, 'gauche.py'))
        self.groups['int15'].save(os.path.join(path, 'int15.py'))
        self.groups['ring'].save(os.path.join(path, 'ring.py'))
        self.groups['radical'].save(os.path.join(path, 'radical.py'))
        self.groups['other'].save(os.path.join(path, 'other.py'))

    def loadOld(self, path):
        """
        Load the old RMG thermo database from the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        # The old database does not have a depository, so create an empty one
        self.depository = {}
        self.depository['stable']  = ThermoDepository(label='stable', name='Stable Molecules')
        self.depository['radical'] = ThermoDepository(label='radical', name='Radical Molecules')
        
        for (root, dirs, files) in os.walk(os.path.join(path, 'thermo_libraries')):
            if os.path.exists(os.path.join(root, 'Dictionary.txt')) and os.path.exists(os.path.join(root, 'Library.txt')):
                library = ThermoLibrary(label=os.path.basename(root), name=os.path.basename(root))
                library.loadOld(
                    dictstr = os.path.join(root, 'Dictionary.txt'),
                    treestr = '',
                    libstr = os.path.join(root, 'Library.txt'),
                    numParameters = 12,
                    numLabels = 1,
                    pattern = False,
                )
                library.label = os.path.basename(root)
                self.libraries[library.label] = library

        self.groups = {}
        self.groups['group'] = ThermoGroups(label='group', name='Functional Group Additivity Values').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', 'Group_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', 'Group_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', 'Group_Library.txt'),
            numParameters = 12,
            numLabels = 1,
            pattern = True,
        )
        self.groups['gauche'] = ThermoGroups(label='gauche', name='Gauche Interaction Corrections').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', 'Gauche_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', 'Gauche_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', 'Gauche_Library.txt'),
            numParameters = 12,
            numLabels = 1,
            pattern = True,
        )
        self.groups['int15'] = ThermoGroups(label='int15', name='1,5-Interaction Corrections').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', '15_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', '15_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', '15_Library.txt'),
            numParameters = 12,
            numLabels = 1,
            pattern = True,
        )
        self.groups['radical'] = ThermoGroups(label='radical', name='Radical Corrections').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', 'Radical_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', 'Radical_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', 'Radical_Library.txt'),
            numParameters = 12,
            numLabels = 1,
            pattern = True,
        )
        self.groups['ring'] = ThermoGroups(label='ring', name='Ring Corrections').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', 'Ring_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', 'Ring_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', 'Ring_Library.txt'),
            numParameters = 12,
            numLabels = 1,
            pattern = True,
        )
        self.groups['other'] = ThermoGroups(label='other', name='Other Corrections').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', 'Other_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', 'Other_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', 'Other_Library.txt'),
            numParameters = 12,
            numLabels = 1,
            pattern = True,
        )

    def saveOld(self, path):
        """
        Save the old RMG thermo database to the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """

        # Depository not used in old database, so it is not saved

        librariesPath = os.path.join(path, 'thermo_libraries')
        if not os.path.exists(librariesPath): os.mkdir(librariesPath)
        for library in self.libraries.values():
            libraryPath = os.path.join(librariesPath, library.label)
            if not os.path.exists(libraryPath): os.mkdir(libraryPath)
            library.saveOld(
                dictstr = os.path.join(libraryPath, 'Dictionary.txt'),
                treestr = '',
                libstr = os.path.join(libraryPath, 'Library.txt'),
            )

        groupsPath = os.path.join(path, 'thermo_groups')
        if not os.path.exists(groupsPath): os.mkdir(groupsPath)
        self.groups['group'].saveOld(
            dictstr = os.path.join(groupsPath, 'Group_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Group_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Group_Library.txt'),
        )
        self.groups['gauche'].saveOld(
            dictstr = os.path.join(groupsPath, 'Gauche_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Gauche_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Gauche_Library.txt'),
        )
        self.groups['int15'].saveOld(
            dictstr = os.path.join(groupsPath, '15_Dictionary.txt'),
            treestr = os.path.join(groupsPath, '15_Tree.txt'),
            libstr = os.path.join(groupsPath, '15_Library.txt'),
        )
        self.groups['radical'].saveOld(
            dictstr = os.path.join(groupsPath, 'Radical_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Radical_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Radical_Library.txt'),
        )
        self.groups['ring'].saveOld(
            dictstr = os.path.join(groupsPath, 'Ring_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Ring_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Ring_Library.txt'),
        )
        self.groups['other'].saveOld(
            dictstr = os.path.join(groupsPath, 'Other_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Other_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Other_Library.txt'),
        )

    def getThermoData(self, species):
        """
        Return the thermodynamic parameters for a given :class:`Species`
        object `species`. This function first searches the loaded libraries
        in order, returning the first match found, before falling back to
        estimation via group additivity.
        """
        thermoData = None
        # Check the libraries in order first; return the first successful match
        for label in self.libraryOrder:
            thermoData = self.getThermoDataFromLibrary(species, self.libraries[label])
            if thermoData is not None: break
        else:
            # Thermo not found in any loaded libraries, so estimate
            thermoData = self.getThermoDataFromGroups(species)
        return thermoData[0]

    def getAllThermoData(self, species):
        """
        Return all possible sets of thermodynamic parameters for a given
        :class:`Species` object `species`. The hits from the depository come
        first, then the libraries (in order), and then the group additivity
        estimate. This method is useful for a generic search job.
        """
        thermoData = []
        # Data from depository comes first
        thermoData.extend(self.getThermoDataFromDepository(species))
        # Data from libraries comes second
        for label in self.libraryOrder:
            data = self.getThermoDataFromLibrary(species, self.libraries[label])
            if data: 
                thermoData.append(data)
        # Last entry is always the estimate from group additivity
        thermoData.append(self.getThermoDataFromGroups(species))
        return thermoData

    def getThermoDataFromDepository(self, species):
        """
        Return all possible sets of thermodynamic parameters for a given
        :class:`Species` object `species` from the depository. If no
        depository is loaded, a :class:`DatabaseError` is raised.
        """
        items = []
        for label, entry in self.depository['stable'].entries.iteritems():
            for molecule in species.molecule:
                if molecule.isIsomorphic(entry.item):
                    items.append((deepcopy(entry.data), self.depository['stable'], entry))
                    break
        for label, entry in self.depository['radical'].entries.iteritems():
            for molecule in species.molecule:
                if molecule.isIsomorphic(entry.item):
                    items.append((deepcopy(entry.data), self.depository['radical'], entry))
                    break
        return items

    def getThermoDataFromLibrary(self, species, library):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Species` object `species` from the specified thermodynamics
        `library`. If `library` is a string, the list of libraries is searched
        for a library with that name. If no match is found in that library,
        ``None`` is returned. If no corresponding library is found, a
        :class:`DatabaseError` is raised.
        """
        for label, entry in library.entries.iteritems():
            for molecule in species.molecule:
                if molecule.isIsomorphic(entry.item) and entry.data is not None:
                    return (deepcopy(entry.data), library, entry)
        return None

    def getThermoDataFromGroups(self, species):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Species` object `species` by estimation using the group
        additivity values. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        """
        # Ensure molecules are using explicit hydrogens
        implicitH = [mol.implicitHydrogens for mol in species.molecule]
        for molecule in species.molecule:
            molecule.makeHydrogensExplicit()
        
        thermo = []
        for molecule in species.molecule:
            molecule.clearLabeledAtoms()
            molecule.updateAtomTypes()
            tdata = self.estimateThermoViaGroupAdditivity(molecule)
            thermo.append(tdata)

        H298 = numpy.array([t.getEnthalpy(298.) for t in thermo])
        indices = H298.argsort()
        
        species.molecule = [species.molecule[ind] for ind in indices]
        implicitH = [implicitH[ind] for ind in indices]
        
        # Restore implicit hydrogens if necessary
        for implicit, molecule in zip(implicitH, species.molecule):
            if implicit: molecule.makeHydrogensImplicit()
        
        return (thermo[indices[0]], None, None)
        
    def estimateThermoViaGroupAdditivity(self, molecule):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Molecule` object `molecule` by estimation using the group
        additivity values. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        """
        implicitH = molecule.implicitHydrogens
        molecule.makeHydrogensExplicit()

        # For thermo estimation we need the atoms to already be sorted because we
        # iterate over them; if the order changes during the iteration then we
        # will probably not visit the right atoms, and so will get the thermo wrong
        molecule.sortVertices()

        thermoData = None

        if sum([atom.radicalElectrons for atom in molecule.atoms]) > 0: # radical species

            # Make a copy of the structure so we don't change the original
            saturatedStruct = molecule.copy(deep=True)

            # Saturate structure by replacing all radicals with bonds to
            # hydrogen atoms
            added = {}
            for atom in saturatedStruct.atoms:
                for i in range(atom.radicalElectrons):
                    H = Atom('H')
                    bond = Bond('S')
                    saturatedStruct.addAtom(H)
                    saturatedStruct.addBond(atom, H, bond)
                    if atom not in added:
                        added[atom] = []
                    added[atom].append([H, bond])
                    atom.decrementRadical()

            # Update the atom types of the saturated structure (not sure why
            # this is necessary, because saturating with H shouldn't be
            # changing atom types, but it doesn't hurt anything and is not
            # very expensive, so will do it anyway)
            saturatedStruct.updateConnectivityValues()
            saturatedStruct.sortVertices()
            saturatedStruct.updateAtomTypes()

            # Get thermo estimate for saturated form of structure
            thermoData = self.estimateThermoViaGroupAdditivity(saturatedStruct)
            assert thermoData is not None, "Thermo data of saturated {0} of molecule {1} is None!".format(saturatedStruct, molecule)
            # Undo symmetry number correction for saturated structure
            thermoData.S298.value += constants.R * math.log(saturatedStruct.symmetryNumber)

            # For each radical site, get radical correction
            # Only one radical site should be considered at a time; all others
            # should be saturated with hydrogen atoms
            for atom in added:

                # Remove the added hydrogen atoms and bond and restore the radical
                for H, bond in added[atom]:
                    saturatedStruct.removeBond(atom, H)
                    saturatedStruct.removeAtom(H)
                    atom.incrementRadical()

                saturatedStruct.updateConnectivityValues()
                
                try:
                    thermoData = self.__addThermoData(thermoData, self.__getGroupThermoData(self.groups['radical'], saturatedStruct, {'*':atom}))
                except KeyError:
                    logging.error("Couldn't find in radical thermo database:")
                    logging.error(molecule)
                    logging.error(molecule.toAdjacencyList())
                    raise
                        
                # Re-saturate
                for H, bond in added[atom]:
                    saturatedStruct.addAtom(H)
                    saturatedStruct.addBond(atom, H, bond)
                    atom.decrementRadical()

                # Subtract the enthalpy of the added hydrogens
                for H, bond in added[atom]:
                    thermoData.H298.value -= 52.103 * 4184

            # Correct the entropy for the symmetry number

        else: # non-radical species
            # Generate estimate of thermodynamics
            for atom in molecule.atoms:
                # Iterate over heavy (non-hydrogen) atoms
                if atom.isNonHydrogen():
                    # Get initial thermo estimate from main group database
                    try:
                        if thermoData is None:
                            thermoData = self.__getGroupThermoData(self.groups['group'], molecule, {'*':atom})
                        else:
                            thermoData = self.__addThermoData(thermoData, self.__getGroupThermoData(self.groups['group'], molecule, {'*':atom}))
                    except KeyError:
                        logging.error("Couldn't find in main thermo database:")
                        logging.error(molecule)
                        logging.error(molecule.toAdjacencyList())
                        raise
                    # Correct for gauche and 1,5- interactions
                    try:
                        thermoData = self.__addThermoData(thermoData, self.__getGroupThermoData(self.groups['gauche'], molecule, {'*':atom}))
                    except KeyError: pass
                    try:
                        thermoData = self.__addThermoData(thermoData, self.__getGroupThermoData(self.groups['int15'], molecule, {'*':atom}))
                    except KeyError: pass
                    try:
                        thermoData = self.__addThermoData(thermoData, self.__getGroupThermoData(self.groups['other'], molecule, {'*':atom}))
                    except KeyError: pass

            # Do ring corrections separately because we only want to match
            # each ring one time; this doesn't work yet
            rings = molecule.getSmallestSetOfSmallestRings()
            for ring in rings:

                # Make a temporary structure containing only the atoms in the ring
                ringStructure = Molecule()
                for atom in ring: ringStructure.addAtom(atom)
                for atom1 in ring:
                    for atom2 in ring:
                        if molecule.hasBond(atom1, atom2):
                            ringStructure.addBond(atom1, atom2, molecule.getBond(atom1, atom2))

                # Get thermo correction for this ring
                try:
                    thermoData = self.__addThermoData(thermoData, self.__getGroupThermoData(self.groups['ring'], ringStructure, {}))
                except KeyError:
                    logging.error("Couldn't find in ring database:")
                    logging.error(ringStructure)
                    logging.error(ringStructure.toAdjacencyList())
                    raise
                
        # Correct entropy for symmetry number
        molecule.calculateSymmetryNumber()
        thermoData.S298.value -= constants.R * math.log(molecule.symmetryNumber)

        if implicitH: molecule.makeHydrogensImplicit()

        return thermoData

    def __addThermoData(self, thermoData1, thermoData2):
        """
        Add two :class:`ThermoData` objects `thermoData1` and `thermoData2`
        together, returning their sum as a new :class:`ThermoData` object.
        """
        if len(thermoData1.Tdata.values) != len(thermoData2.Tdata.values) or any([T1 != T2 for T1, T2 in zip(thermoData1.Tdata.values, thermoData2.Tdata.values)]):
            raise ThermoError('Cannot add these ThermoData objects due to their having different temperature points.')
        new = ThermoData(
            Tdata = (thermoData1.Tdata.values, thermoData1.Tdata.units),
            Cpdata = (thermoData1.Cpdata.values + thermoData2.Cpdata.values, thermoData1.Tdata.units),
            H298 = (thermoData1.H298.value + thermoData2.H298.value, thermoData1.Tdata.units),
            S298 = (thermoData1.S298.value + thermoData2.S298.value, thermoData1.Tdata.units),
        )
        if thermoData1.comment == '': new.comment = thermoData2.comment
        elif thermoData2.comment == '': new.comment = thermoData1.comment
        else: new.comment = thermoData1.comment + ' + ' + thermoData2.comment
        return new
    
    def __getGroupThermoData(self, database, molecule, atom):
        """
        Determine the group additivity thermodynamic data for the atom `atom`
        in the structure `structure`.
        """

        node0 = database.descendTree(molecule, atom, None)

        if node0 is None:
            raise KeyError('Node not found in database.')

        # It's possible (and allowed) that items in the tree may not be in the
        # library, in which case we need to fall up the tree until we find an
        # ancestor that has an entry in the library
        node = node0
        while node.data is None and node is not None:
            node = node.parent
        if node is None:
            raise InvalidDatabaseError('Unable to determine thermo parameters for {0}: no library entries for {1} or any of its ancestors.'.format(molecule, node0) )

        data = node.data; comment = node.label
        while isinstance(data, str) and data is not None:
            for entry in database.entries.values():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    break
        data = deepcopy(data)
        data.comment = '{0}({1})'.format(database.label, comment)

        # This code prints the hierarchy of the found node; useful for debugging
        #result = ''
        #while node is not None:
        #   result = ' -> ' + node + result
        #   node = database.tree.parent[node]
        #print result[4:]

        return data
