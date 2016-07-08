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

import os.path
import re
import math
import logging
import numpy
from copy import deepcopy

from base import Database, Entry, makeLogicNode, DatabaseError

import rmgpy.constants as constants
from rmgpy.thermo import NASAPolynomial, NASA, ThermoData, Wilhoit
from rmgpy.molecule import Molecule, Atom, Bond, Group
import rmgpy.molecule
from rmgpy.species import Species

from rmgpy.scoop_framework.util import get

#: This dictionary is used to add multiplicity to species label
_multiplicity_labels = {1:'S',2:'D',3:'T',4:'Q',5:'V',}

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
    elif isinstance(entry.data, NASA):
        f.write('    thermo = NASA(\n')
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
    f.write('    shortDesc = u"""')
    try:
        f.write(entry.shortDesc.encode('utf-8'))
    except (UnicodeEncodeError, UnicodeDecodeError):
        f.write(entry.shortDesc.strip().encode('ascii', 'replace'))
    f.write('""",\n')
    f.write('    longDesc = \n')
    f.write('u"""\n')
    try:
        f.write(entry.longDesc.strip().encode('utf-8') + "\n")    
    except (UnicodeEncodeError, UnicodeDecodeError):
        f.write(entry.longDesc.strip().encode('ascii', 'replace') + "\n")
    f.write('""",\n')
    if entry.rank:
        f.write("    rank = {0},\n".format(entry.rank))

    f.write(')\n\n')

def generateOldLibraryEntry(data):
    """
    Return a list of values used to save entries to the old-style RMG
    thermo database based on the thermodynamics object `data`.
    """
    if isinstance(data, ThermoData):
        return '{0:9g} {1:9g} {2:9g} {3:9g} {4:9g} {5:9g} {6:9g} {7:9g} {8:9g} {9:9g} {10:9g} {11:9g}'.format(
            data.H298.value_si/4184.,
            data.S298.value_si/4.184,
            data.Cpdata.value_si[0]/4.184,
            data.Cpdata.value_si[1]/4.184,
            data.Cpdata.value_si[2]/4.184,
            data.Cpdata.value_si[3]/4.184,
            data.Cpdata.value_si[4]/4.184,
            data.Cpdata.value_si[5]/4.184,
            data.Cpdata.value_si[6]/4.184,
            data.H298.uncertainty/4184.,
            data.S298.uncertainty/4.184,
            max(data.Cpdata.uncertainty)/4.184,
        )
    elif isinstance(data, basestring):
        return data
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


def addThermoData(thermoData1, thermoData2, groupAdditivity=False):
        """
        Add the thermodynamic data `thermoData2` to the data `thermoData1`,
        and return `thermoData1`.
        
        If `groupAdditivity` is True, append comments related to group additivity estimation
        """
        if len(thermoData1.Tdata.value_si) != len(thermoData2.Tdata.value_si) or any([T1 != T2 for T1, T2 in zip(thermoData1.Tdata.value_si, thermoData2.Tdata.value_si)]):
            raise Exception('Cannot add these ThermoData objects due to their having different temperature points.')
        
        for i in range(thermoData1.Tdata.value_si.shape[0]):
            thermoData1.Cpdata.value_si[i] += thermoData2.Cpdata.value_si[i]
        thermoData1.H298.value_si += thermoData2.H298.value_si
        thermoData1.S298.value_si += thermoData2.S298.value_si

        if groupAdditivity:
            if thermoData1.comment:
                thermoData1.comment += ' + {0}'.format(thermoData2.comment)
            else:
                thermoData1.comment = 'Thermo group additivity estimation: ' + thermoData2.comment
            
        return thermoData1
    
def removeThermoData(thermoData1, thermoData2):
    """
    Remove the thermodynamic data `thermoData2` from the data `thermoData1`,
    and return `thermoData1`.
    """
    if len(thermoData1.Tdata.value_si) != len(thermoData2.Tdata.value_si) or any([T1 != T2 for T1, T2 in zip(thermoData1.Tdata.value_si, thermoData2.Tdata.value_si)]):
        raise Exception('Cannot take the difference between these ThermoData objects due to their having different temperature points.')

    for i in range(thermoData1.Tdata.value_si.shape[0]):
        thermoData1.Cpdata.value_si[i] -= thermoData2.Cpdata.value_si[i]
    thermoData1.H298.value_si -= thermoData2.H298.value_si
    thermoData1.S298.value_si -= thermoData2.S298.value_si

    return thermoData1

def averageThermoData(thermoDataList=None):
    """
    Average a list of thermoData values together.
    Sets uncertainty values to be the approximately the 95% confidence interval, equivalent to
    2 standard deviations calculated using the sample standard variance:
    
    Uncertainty = 2s
    s = sqrt( sum(abs(x - x.mean())^2) / N - 1) where N is the number of values averaged
    
    Note that uncertainties are only computed when number of values is greater than 1.
    """
    if thermoDataList is None:
        thermoDataList = []
        
    import copy
    numValues = len(thermoDataList)
        
    if numValues == 0:
        raise Exception('No thermo data values were inputted to be averaged.')
    else:
        print 'Averaging thermo data over {0} value(s).'.format(numValues)
        
        if numValues == 1:
            return copy.deepcopy(thermoDataList[0])
        
        else:
            averagedThermoData = copy.deepcopy(thermoDataList[0])
            for thermoData in thermoDataList[1:]:
                averagedThermoData = addThermoData(averagedThermoData, thermoData)


            for i in range(averagedThermoData.Tdata.value_si.shape[0]):
                averagedThermoData.Cpdata.value_si[i] /= numValues

                cpData = [thermoData.Cpdata.value_si[i] for thermoData in thermoDataList]
                averagedThermoData.Cpdata.uncertainty[i] = 2*numpy.std(cpData, ddof=1)

            HData = [thermoData.H298.value_si for thermoData in thermoDataList]
            averagedThermoData.H298.value_si /= numValues
            averagedThermoData.H298.uncertainty_si = 2*numpy.std(HData, ddof=1)

            SData = [thermoData.S298.value_si for thermoData in thermoDataList]
            averagedThermoData.S298.value_si /= numValues
            averagedThermoData.S298.uncertainty_si = 2*numpy.std(SData, ddof=1)
            return averagedThermoData

################################################################################

class ThermoDepository(Database):
    """
    A class for working with the RMG thermodynamics depository.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self, index, label, molecule, thermo, reference=None, referenceType='', shortDesc='', longDesc='', rank=None):
        entry = Entry(
            index = index,
            label = label,
            item = Molecule().fromAdjacencyList(molecule),
            data = thermo,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            rank = rank,
        )
        self.entries[label] = entry
        return entry

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

    def __init__(self, label='', name='',solvent=None, shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self,
                  index,
                  label,
                  molecule,
                  thermo,
                  reference=None,
                  referenceType='',
                  shortDesc='',
                  longDesc='',
                  rank=None,
                  ):
        
        molecule = Molecule().fromAdjacencyList(molecule)
        
        # Internal checks for adding entry to the thermo library
        if label in self.entries.keys():
            raise DatabaseError('Found a duplicate molecule with label {0} in the thermo library {1}.  Please correct your library.'.format(label, self.name))
        
        for entry in self.entries.values():
            if molecule.isIsomorphic(entry.item):
                if molecule.multiplicity == entry.item.multiplicity:
                    raise DatabaseError('Adjacency list and multiplicity of {0} matches that of existing molecule {1} in thermo library {2}.  Please correct your library.'.format(label, entry.label, self.name))
        
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = molecule,
            data = thermo,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            rank = rank,
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

    def loadEntry(self,
                  index,
                  label,
                  group,
                  thermo,
                  reference=None,
                  referenceType='',
                  shortDesc='',
                  longDesc='',
                  rank=None,
                  ):

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
            rank = rank,
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

    def copyData(self, source, destination):
        """
        This method copys the ThermoData object and all meta data
        from source to destination
        Args:
            source: The entry for which data is being copied
            destination: The entry for which data is being overwritten

        """
        destination.data = source.data
        destination.reference = source.reference
        destination._longDesc = source._longDesc
        destination._shortDesc = source._shortDesc
        destination._longDesc = source.longDesc
        destination._shortDesc = source.shortDesc
        destination.rank = source.rank
        destination.referenceType = source.referenceType    


    def removeGroup(self, groupToRemove):
        """
        Removes a group that is in a tree from the database. For thermo
        groups we also, need to re-point any unicode thermoData that may
        have pointed to the entry.

        Returns the removed group
        """

        #First call base class method
        Database.removeGroup(self, groupToRemove)

        parentR = groupToRemove.parent

        #look for other pointers that point toward entry
        for entryName, entry in self.entries.iteritems():
            if isinstance(entry.data, basestring):
                if entry.data == groupToRemove.label:
                    #if the entryToRemove.data is also a pointer, then copy
                    if isinstance(groupToRemove.data, basestring):
                        entry.data = groupToRemove.data
                    #if the parent points toward entry and the data is
                    #not a base string, we need to copy the data to the parent
                    elif entry is parentR:
                        self.copyData(groupToRemove, parentR)
                    #otherwise, point toward entryToRemove's parent
                    else:
                        entry.data = unicode(parentR.label)

        return groupToRemove
################################################################################

class ThermoDatabase(object):
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
            'NASAPolynomial': NASAPolynomial,
            'NASA': NASA,
        }
        self.global_context = {}

    def __reduce__(self):
        """
        A helper function used when pickling a ThermoDatabase object.
        """
        d = {
            'depository': self.depository,
            'libraries': self.libraries,
            'groups': self.groups,
            'libraryOrder': self.libraryOrder,
        }
        return (ThermoDatabase, (), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling a ThermoDatabase object.
        """
        self.depository = d['depository']
        self.libraries = d['libraries']
        self.groups = d['groups']
        self.libraryOrder = d['libraryOrder']

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
        self.groups['polycyclic'] = ThermoGroups(label='polycyclic').load(os.path.join(path, 'polycyclic.py'), self.local_context, self.global_context)
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
        for depo in self.depository.keys():
            self.depository[depo].save(os.path.join(path, depo+'.py'))

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
        for group in self.groups.keys():
            self.groups[group].save(os.path.join(path, group+'.py'))
            
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
        self.groups['polycyclic'] = ThermoGroups(label='other', name='Polycyclic Ring Corrections').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', 'Polycyclic_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', 'Polycyclic_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', 'Polycyclic_Library.txt'),
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
        
    def pruneHeteroatoms(self, allowed=['C','H','O','S']):
        """
        Remove all species from thermo libraries that contain atoms other than those allowed.
        
        This is useful before saving the database for use in RMG-Java
        """
        allowedElements = [rmgpy.molecule.element.getElement(label) for label in allowed]
        for library in self.libraries.values():
            logging.info("Removing hetoroatoms from thermo library '{0}'".format(library.name))
            toDelete = []
            for entry in library.entries.values():
                for atom in entry.item.atoms:
                    if atom.element not in allowedElements:
                        toDelete.append(entry.label)
                        break
            for label in toDelete:
                logging.info(" {0}".format(label))
                library.entries.pop(label)

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
        self.groups['polycyclic'].saveOld(
            dictstr = os.path.join(groupsPath, 'Polycyclic_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Polycyclic_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Polycyclic_Library.txt'),
        )
        self.groups['other'].saveOld(
            dictstr = os.path.join(groupsPath, 'Other_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Other_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Other_Library.txt'),
        )


    def getThermoData(self, species, trainingSet=None):
        """
        Return the thermodynamic parameters for a given :class:`Species`
        object `species`. This function first searches the loaded libraries
        in order, returning the first match found, before falling back to
        estimation via group additivity.
        
        Returns: ThermoData
        """
        from rmgpy.rmg.input import getInput
        
        thermo0 = None
        
        thermo0 = self.getThermoDataFromLibraries(species)
        try:
            quantumMechanics = getInput('quantumMechanics')
        except Exception, e:
            logging.debug('Quantum Mechanics DB could not be found.')
            quantumMechanics = None

        if thermo0 is not None:
            logging.info("Found thermo for {0} in {1}".format(species.label,thermo0[0].comment.lower()))
            assert len(thermo0) == 3, "thermo0 should be a tuple at this point: (thermoData, library, entry)"
            thermo0 = thermo0[0]
            
        elif quantumMechanics:
            original_molecule = species.molecule[0]
            if quantumMechanics.settings.onlyCyclics and not original_molecule.isCyclic():
                pass
            else: # try a QM calculation
                if original_molecule.getRadicalCount() > quantumMechanics.settings.maxRadicalNumber:
                    # Too many radicals for direct calculation: use HBI.
                    logging.info("{0} radicals on {1} exceeds limit of {2}. Using HBI method.".format(
                        original_molecule.getRadicalCount(),
                        species.label,
                        quantumMechanics.settings.maxRadicalNumber,
                        ))
                    
                    # Need to estimate thermo via each resonance isomer
                    thermo = []
                    for molecule in species.molecule:
                        molecule.clearLabeledAtoms()
                        # Try to see if the saturated molecule can be found in the libraries
                        tdata = self.estimateRadicalThermoViaHBI(molecule, self.getThermoDataFromLibraries)
                        priority = 1
                        if tdata is None:
                            # Then attempt quantum mechanics job on the saturated molecule
                            tdata = self.estimateRadicalThermoViaHBI(molecule, quantumMechanics.getThermoData)
                            priority = 2
                        if tdata is None:
                            # Fall back to group additivity
                            tdata = self.estimateThermoViaGroupAdditivity(molecule)
                            priority = 3
                        
                        thermo.append((priority, tdata.getEnthalpy(298.), molecule, tdata))
                    
                    if len(thermo) > 1:
                        # Sort thermo first by the priority, then by the most stable H298 value
                        thermo = sorted(thermo, key=lambda x: (x[0], x[1])) 
                        for i in range(len(thermo)): 
                            logging.info("Resonance isomer {0} {1} gives H298={2:.0f} J/mol".format(i+1, thermo[i][2].toSMILES(), thermo[i][1]))
                        # Save resonance isomers reordered by their thermo
                        species.molecule = [item[2] for item in thermo]
                        original_molecule = species.molecule[0]
                    thermo0 = thermo[0][3] 
                    
                    # If priority == 2
                    if thermo[0][0] == 2:
                        # Write the QM molecule thermo to a library so that can be used in future RMG jobs.  (Do this only if it came from a QM calculation)
                        quantumMechanics.database.loadEntry(index = len(quantumMechanics.database.entries) + 1,
                                                        label = original_molecule.toSMILES() + '_({0})'.format(_multiplicity_labels[original_molecule.multiplicity]),
                                                        molecule = original_molecule.toAdjacencyList(),
                                                        thermo = thermo0,
                                                        shortDesc = thermo0.comment
                                                        
                                                        )                    
#                    # For writing thermodata HBI check for QM molecules
#                    with open('thermoHBIcheck.txt','a') as f:
#                        f.write('// {0!r}\n'.format(thermo0).replace('),','),\n//           '))
#                        f.write('{0}\n'.format(original_molecule.toSMILES()))
#                        f.write('{0}\n\n'.format(original_molecule.toAdjacencyList(removeH=False)))

                else: # Not too many radicals: do a direct calculation.
                    thermo0 = quantumMechanics.getThermoData(original_molecule) # returns None if it fails
                
                    if thermo0 is not None:
                        # Write the QM molecule thermo to a library so that can be used in future RMG jobs.
                        quantumMechanics.database.loadEntry(index = len(quantumMechanics.database.entries) + 1,
                                                        label = original_molecule.toSMILES() + '_({0})'.format(_multiplicity_labels[original_molecule.multiplicity]),
                                                        molecule = original_molecule.toAdjacencyList(),
                                                        thermo = thermo0,
                                                        shortDesc = thermo0.comment
                                                        )                    
        if thermo0 is None:
            # Use group additivity methods to determine thermo for molecule (or if QM fails completely)
            original_molecule = species.molecule[0]
            if original_molecule.getRadicalCount() > 0:
                # Molecule is a radical, use the HBI method
                thermo = []
                for molecule in species.molecule:
                    molecule.clearLabeledAtoms()
                    # First see if the saturated molecule is in the libaries
                    tdata = self.estimateRadicalThermoViaHBI(molecule, self.getThermoDataFromLibraries)
                    priority = 1
                    if tdata is None:
                        # Otherwise use normal group additivity to obtain the thermo for the molecule
                        tdata = self.estimateThermoViaGroupAdditivity(molecule)
                        priority = 2
                    thermo.append((priority, tdata.getEnthalpy(298.), molecule, tdata))
                
                if len(thermo) > 1:
                    # Sort thermo first by the priority, then by the most stable H298 value
                    thermo = sorted(thermo, key=lambda x: (x[0], x[1]))
                    for i in range(len(thermo)): 
                        logging.info("Resonance isomer {0} {1} gives H298={2:.0f} J/mol".format(i+1, thermo[i][2].toSMILES(), thermo[i][1]))
                    # Save resonance isomers reordered by their thermo
                    species.molecule = [item[2] for item in thermo]
                thermo0 = thermo[0][3] 
            else:
                # Saturated molecule, does not need HBI method
                thermo0 = self.getThermoDataFromGroups(species)
                
        # Make sure to calculate Cp0 and CpInf if it wasn't done already
        self.findCp0andCpInf(species, thermo0)

        # Return the resulting thermo parameters
        return thermo0
    
        
    def getThermoDataFromLibraries(self, species, trainingSet=None):
        """
        Return the thermodynamic parameters for a given :class:`Species`
        object `species`. This function first searches the loaded libraries
        in order, returning the first match found, before failing and returning None.
        `trainingSet` is used to identify if function is called during training set or not.
        During training set calculation we want to use gas phase thermo to not affect reverse
        rate calculation.
        
        Returns: ThermoData or None
        """
        import rmgpy.rmg.main
        thermoData = None
        
        #chatelak 11/15/14: modification to introduce liquid phase thermo libraries
        libraryList=deepcopy(self.libraryOrder) #copy the value to not affect initial object

        if rmgpy.rmg.main.solvent is not None:
            liqLibraries=[]
            #Liquid phase simulation part: 
            #This bloc "for": Identify liquid phase libraries and store them in liqLibraries
            for iterLib in libraryList:
                if self.libraries[iterLib].solvent:
                    liqLibraries.append(iterLib)
            #Check in liqLibraries if thermo for species exists and return the first match. Only if function not called by trainingSet
            if liqLibraries and trainingSet is None:
                for label in liqLibraries:
                    thermoData = self.getThermoDataFromLibrary(species, self.libraries[label])
                    if thermoData is not None:
                        assert len(thermoData) == 3, "thermoData should be a tuple at this point"
                        #Watch out comments changed: this is used later to apply solvation or not on species matching thermo. If required, Modify this carefully.
                        thermoData[0].comment += 'Liquid thermo library: ' + label
                        return thermoData
            #Remove liqLibraries from libraryList if: called by training set (trainingSet=True) or if no thermo found in liqLibrairies
            #if no liquid library found this does nothing.   
            for libIter in liqLibraries:
                libraryList.remove(libIter)

        # Condition to execute this part: gas phase simulation or training set or liquid phase simulation with : noliquid libraries found or no matching species found in liquid libraries       
        # If gas phase simulation libraryList = self.libraryOrder (just like before modifications) and they are all gas phase, already checked by checkLibrairies function in database.load()
        # Check the libraries in order; return the first successful match
        for label in libraryList:
            thermoData = self.getThermoDataFromLibrary(species, self.libraries[label])
            if thermoData is not None:
                assert len(thermoData) == 3, "thermoData should be a tuple at this point"
                if rmgpy.rmg.main.solvent is not None and trainingSet is None:
                    thermoData[0].comment += 'Thermo library corrected for liquid phase: ' + label
                else:
                    thermoData[0].comment += 'Thermo library: ' + label
                return thermoData

        return None
    
    def findCp0andCpInf(self, species, thermoData):
        """
        Calculate the Cp0 and CpInf values, and add them to the thermoData object.
        
        Modifies thermoData in place and doesn't return anything
        """
        if not isinstance(thermoData,ThermoData):
            return # Just skip it
            raise Exception("Trying to add Cp0 to something that's not a ThermoData: {0!r}".format(thermoData))
        if thermoData.Cp0 is None:
            Cp0 = species.calculateCp0()
            thermoData.Cp0 = (Cp0,"J/(mol*K)")
        if thermoData.CpInf is None:
            CpInf = species.calculateCpInf()  
            thermoData.CpInf = (CpInf,"J/(mol*K)")
                
                
    def getAllThermoData(self, species):
        """
        Return all possible sets of thermodynamic parameters for a given
        :class:`Species` object `species`. The hits from the depository come
        first, then the libraries (in order), and then the group additivity
        estimate. This method is useful for a generic search job.
        
        Returns: a list of tuples (ThermoData, source, entry) 
        (Source is a library or depository, or None)
        """
        thermoDataList = []
        # Data from depository comes first
        thermoDataList.extend(self.getThermoDataFromDepository(species))
        # Data from libraries comes second
        for label in self.libraryOrder:
            data = self.getThermoDataFromLibrary(species, self.libraries[label])
            if data: 
                assert len(data) == 3, "thermoData should be a tuple at this point"
                data[0].comment += label
                thermoDataList.append(data)
        # Last entry is always the estimate from group additivity
        # Make it a tuple
        data = (self.getThermoDataFromGroups(species), None, None)
        thermoDataList.append(data)

        # Return all of the resulting thermo parameters
        return thermoDataList

    def getThermoDataFromDepository(self, species):
        """
        Return all possible sets of thermodynamic parameters for a given
        :class:`Species` object `species` from the depository. If no
        depository is loaded, a :class:`DatabaseError` is raised.
        
        Returns: a list of tuples (thermoData, depository, entry) without any Cp0 or CpInf data.
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
        
        Returns a tuple: (ThermoData, library, entry)  or None.
        """
        for label, entry in library.entries.iteritems():
            for molecule in species.molecule:
                if molecule.isIsomorphic(entry.item) and entry.data is not None:
                    thermoData = deepcopy(entry.data)
                    self.findCp0andCpInf(species, thermoData)
                    return (thermoData, library, entry)
        return None

    def getThermoDataFromGroups(self, species):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Species` object `species` by estimation using the group
        additivity values. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        
        The resonance isomer (molecule) with the lowest H298 is used, and as a side-effect
        the resonance isomers (items in `species.molecule` list) are sorted in ascending order.
        
        Returns: ThermoData
        """       
        thermo = []
        for molecule in species.molecule:
            molecule.clearLabeledAtoms()
            molecule.updateAtomTypes()
            tdata = self.estimateThermoViaGroupAdditivity(molecule)
            thermo.append(tdata)

        indices = self.prioritizeThermo(species, thermo)
        
        species.molecule = [species.molecule[ind] for ind in indices]
        
        thermoData = thermo[indices[0]]
        self.findCp0andCpInf(species, thermoData)
        return thermoData

    def prioritizeThermo(self, species, thermoDataList):
        """
        Use some metrics to reorder a list of thermo data from best to worst.
        Return a list of indices with the desired order associated with the index of thermo from the data list.
        """
        if len(species.molecule) > 1:
            # Go further only if there is more than one isomer
            if species.molecule[0].isCyclic():
                # Special treatment for cyclic compounds
                entries = []
                for thermo in thermoDataList:
                    ringGroups, polycyclicGroups = self.getRingGroupsFromComments(thermo)
                    
                    # Use rank as a metric for prioritizing thermo. 
                    # The smaller the rank, the better.
                    sumRank = numpy.sum([3 if entry.rank is None else entry.rank for entry in ringGroups + polycyclicGroups])
                    entries.append((thermo, sumRank))
                
                # Sort first by rank, then by enthalpy at 298 K
                entries = sorted(entries, key=lambda entry: (entry[1], entry[0].getEnthalpy(298.)))
                indices = [thermoDataList.index(entry[0]) for entry in entries]
                
            else:
                # For noncyclics, default to original algorithm of ordering thermo based on the most stable enthalpy
                H298 = numpy.array([t.getEnthalpy(298.) for t in thermoDataList])
                indices = H298.argsort()
        else:
            indices = [0]

        return indices

    def estimateRadicalThermoViaHBI(self, molecule, stableThermoEstimator ):
        """
        Estimate the thermodynamics of a radical by saturating it,
        applying the provided stableThermoEstimator method on the saturated species,
        then applying hydrogen bond increment corrections for the radical
        site(s) and correcting for the symmetry.
        
        """
        
        assert molecule.isRadical(), "Method only valid for radicals."
        saturatedStruct = molecule.copy(deep=True)
        added = saturatedStruct.saturate()
        saturatedStruct.props['saturated'] = True
        
        # Get thermo estimate for saturated form of structure
        if stableThermoEstimator == self.getThermoDataFromLibraries:
            # Get data from libraries
            saturatedSpec = Species(molecule=[saturatedStruct])
            thermoData_sat = stableThermoEstimator(saturatedSpec)
            if thermoData_sat:
                assert len(thermoData_sat) == 3, "thermoData should be a tuple at this point: (thermoData, library, entry)"
                thermoData_sat = thermoData_sat[0]
        else:
            thermoData_sat = stableThermoEstimator(saturatedStruct)
        if thermoData_sat is None:
            # logging.info("Thermo data of saturated {0} of molecule {1} is None.".format(saturatedStruct, molecule))
            return None
        assert thermoData_sat is not None, "Thermo data of saturated {0} of molecule {1} is None!".format(saturatedStruct, molecule)
        
        # Convert to ThermoData object if necessary in order to add and subtract from enthalpy and entropy values
        if not isinstance(thermoData_sat, ThermoData):
            thermoData_sat = thermoData_sat.toThermoData()
        
        
        if not stableThermoEstimator == self.computeGroupAdditivityThermo:
            #remove the symmetry contribution to the entropy of the saturated molecule
            ##assumes that the thermo data comes from QMTP or from a thermolibrary
            thermoData_sat.S298.value_si += constants.R * math.log(saturatedStruct.getSymmetryNumber())
        
        thermoData = thermoData_sat
        
        # Correct entropy for symmetry number of radical structure
        thermoData.S298.value_si -= constants.R * math.log(molecule.getSymmetryNumber())
        
        # For each radical site, get radical correction
        # Only one radical site should be considered at a time; all others
        # should be saturated with hydrogen atoms
        for atom in added:
            # Remove the added hydrogen atoms and bond and restore the radical
            for H, bond in added[atom]:
                saturatedStruct.removeBond(bond)
                saturatedStruct.removeAtom(H)
                atom.incrementRadical()
            saturatedStruct.update()
            try:
                self.__addGroupThermoData(thermoData, self.groups['radical'], saturatedStruct, {'*':atom})
            except KeyError:
                logging.error("Couldn't find in radical thermo database:")
                logging.error(molecule)
                logging.error(molecule.toAdjacencyList())
                raise
            # Re-saturate
            for H, bond in added[atom]:
                saturatedStruct.addAtom(H)
                saturatedStruct.addBond(bond)
                atom.decrementRadical()
            # Subtract the enthalpy of the added hydrogens
            for H, bond in added[atom]:
                thermoData.H298.value_si -= 52.103 * 4184

        return thermoData
        
        
    def estimateThermoViaGroupAdditivity(self, molecule):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Molecule` object `molecule` by estimation using the group
        additivity values. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        """
        # For thermo estimation we need the atoms to already be sorted because we
        # iterate over them; if the order changes during the iteration then we
        # will probably not visit the right atoms, and so will get the thermo wrong
        molecule.sortAtoms()

        if molecule.isRadical(): # radical species
            thermoData = self.estimateRadicalThermoViaHBI(molecule, self.computeGroupAdditivityThermo)
            return thermoData

        else: # non-radical species
            thermoData = self.computeGroupAdditivityThermo(molecule)
            # Correct entropy for symmetry number
            if not 'saturated' in molecule.props: 
                thermoData.S298.value_si -= constants.R * math.log(molecule.getSymmetryNumber())
            return thermoData


    def computeGroupAdditivityThermo(self, molecule):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Molecule` object `molecule` by estimation using the group
        additivity values. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        
        The entropy is not corrected for the symmetry of the molecule.
        This should be done later by the calling function.
        """

        assert not molecule.isRadical(), "This method is only for saturated non-radical species."
        # For thermo estimation we need the atoms to already be sorted because we
        # iterate over them; if the order changes during the iteration then we
        # will probably not visit the right atoms, and so will get the thermo wrong
        molecule.sortAtoms()

        # Create the ThermoData object
        thermoData = ThermoData(
            Tdata = ([300,400,500,600,800,1000,1500],"K"),
            Cpdata = ([0.0,0.0,0.0,0.0,0.0,0.0,0.0],"J/(mol*K)"),
            H298 = (0.0,"kJ/mol"),
            S298 = (0.0,"J/(mol*K)"),
        )

        cyclic = molecule.isCyclic()
        # Generate estimate of thermodynamics
        for atom in molecule.atoms:
            # Iterate over heavy (non-hydrogen) atoms
            if atom.isNonHydrogen():
                # Get initial thermo estimate from main group database
                try:
                    self.__addGroupThermoData(thermoData, self.groups['group'], molecule, {'*':atom})
                except KeyError:
                    logging.error("Couldn't find in main thermo database:")
                    logging.error(molecule)
                    logging.error(molecule.toAdjacencyList())
                    raise
                # Correct for gauche and 1,5- interactions
                if not cyclic:
                    try:
                        self.__addGroupThermoData(thermoData, self.groups['gauche'], molecule, {'*':atom})
                    except KeyError: pass
                try:
                    self.__addGroupThermoData(thermoData, self.groups['int15'], molecule, {'*':atom})
                except KeyError: pass
                try:
                    self.__addGroupThermoData(thermoData, self.groups['other'], molecule, {'*':atom})
                except KeyError: pass

        # Do ring corrections separately because we only want to match
        # each ring one time
        
        if cyclic:                
            monorings, polyrings = molecule.getDisparateRings()
            for ring in monorings:
                # Make a temporary structure containing only the atoms in the ring
                # NB. if any of the ring corrections depend on ligands not in the ring, they will not be found!
                try:
                    self.__addRingCorrectionThermoData(thermoData, self.groups['ring'], molecule, ring)
                except KeyError:
                    logging.error("Couldn't find a match in the monocyclic ring database even though monocyclic rings were found.")
                    logging.error(molecule)
                    logging.error(molecule.toAdjacencyList())
                    raise
            for ring in polyrings:
                # Make a temporary structure containing only the atoms in the ring
                # NB. if any of the ring corrections depend on ligands not in the ring, they will not be found!
                try:
                    self.__addRingCorrectionThermoData(thermoData, self.groups['polycyclic'], molecule, ring)
                except KeyError:
                    logging.error("Couldn't find a match in the polycyclic ring database even though polycyclic rings were found.")
                    logging.error(molecule)
                    logging.error(molecule.toAdjacencyList())
                    raise

        return thermoData

    def __addRingCorrectionThermoData(self, thermoData, ring_database, molecule, ring):
        """
        Determine the ring correction group additivity thermodynamic data for the given
         `ring` in the `molecule`, and add it to the existing thermo data
        `thermoData`.
        """
        matchedRingEntries = []
        # label each atom in the ring individually to try to match the group
        # for each ring, save only the ring that is matches the most specific leaf in the tree.
        for atom in ring:
            atoms = {'*':atom}
            entry = ring_database.descendTree(molecule, atoms)
            matchedRingEntries.append(entry)
        
        if matchedRingEntries is []:
            raise KeyError('Node not found in database.')
        # Decide which group to keep
        depthList = [len(ring_database.ancestors(entry)) for entry in matchedRingEntries]
        mostSpecificMatchIndices = [i for i, x in enumerate(depthList) if x == max(depthList)]
        
        mostSpecificMatchedEntries = [matchedRingEntries[idx] for idx in mostSpecificMatchIndices]
        if len(set(mostSpecificMatchedEntries)) != 1:
            logging.warning('More than one type of node was found to be most specific for this ring.')
            logging.warning('This is either due to a database error in the ring or polycyclic groups, or a partial match between the group and the full ring.')
            logging.warning(mostSpecificMatchedEntries)
            
        # Condense the number of most specific groups down to one
        mostSpecificMatchedEntry = matchedRingEntries[mostSpecificMatchIndices[0]]
        
        node = mostSpecificMatchedEntry
        while node.data is None and node is not None:
            node = node.parent
        if node is None:
            raise DatabaseError('Unable to determine thermo parameters for {0}: no data for {1} or any of its ancestors.'.format(molecule, mostSpecificGroup) )

        data = node.data; comment = node.label
        while isinstance(data, basestring) and data is not None:
            for entry in ring_database.entries.values():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    break
        data.comment = '{0}({1})'.format(ring_database.label, comment)
        
        if thermoData is None:
            return data
        else:
            return addThermoData(thermoData, data, groupAdditivity=True)

    def __addGroupThermoData(self, thermoData, database, molecule, atom):
        """
        Determine the group additivity thermodynamic data for the atom `atom`
        in the structure `structure`, and add it to the existing thermo data
        `thermoData`.
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
            raise DatabaseError('Unable to determine thermo parameters for {0}: no data for node {1} or any of its ancestors.'.format(molecule, node0) )

        data = node.data; comment = node.label
        while isinstance(data, basestring) and data is not None:
            for entry in database.entries.values():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    break
        data.comment = '{0}({1})'.format(database.label, comment)

        # This code prints the hierarchy of the found node; useful for debugging
#        result = ''
#        while node is not None:
#           result = ' -> ' + node.label + result
#           node = node.parent
#        print result[4:]
        
        if thermoData is None:
            return data
        else:
            return addThermoData(thermoData, data, groupAdditivity=True)

    def getRingGroupsFromComments(self, thermoData):
        """
        Takes a string of comments from group additivity estimation, and extracts the ring and polycyclic ring groups
        from them, returning them as lists.
        """
        tokens = thermoData.comment.split()
        ringGroups = []
        polycyclicGroups = []
        for token in tokens:
            regex = "\((.*)\)" #only hit outermost parentheses
            if 'ring' in token:
                splitTokens = re.split(regex, token)
                assert len(splitTokens) == 3, 'token: {}'.format(token)
                groupLabel = splitTokens[1]
                ringGroups.append(self.groups['ring'].entries[groupLabel])
            if 'polycyclic' in token:
                splitTokens = re.split(regex, token)
                assert len(splitTokens) == 3, 'token: {}'.format(token)
                groupLabel = splitTokens[1]
                polycyclicGroups.append(self.groups['polycyclic'].entries[groupLabel])
        return ringGroups, polycyclicGroups
