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
A class for returning and estimating the transport properties of a species

"""
import os
import os.path
import math
import logging
import numpy
from copy import copy, deepcopy

from base import Database, Entry, makeLogicNode

import rmgpy.constants as constants
from rmgpy.molecule import Molecule, Atom, Bond, Group
from rmgpy.transport import TransportData


def saveEntry(f, entry):
    """
    Write a Pythonic string representation of the given `entry` in the transport
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
        
def generateOldLibraryEntry(data):
    """
    Return a list of values used to save entries to the old-style RMG
    transport database based on the transport object `data`.
    """
    
def processOldLibraryEntry(data):
    """
    Process a list of parameters `data` as read from an old-style RMG
    transport database, returning the corresponding transport object.
    """
    
class TransportLibrary(Database):
    """
    A class for working with a RMG transport library.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self,
                  index,
                  label,
                  molecule,
                  transport,
                  reference=None,
                  referenceType='',
                  shortDesc='',
                  longDesc='',
                  history=None
                  ):
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = Molecule().fromAdjacencyList(molecule),
            data = transport,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            history = history or [],
        )
    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the transport database to the file object `f`.
        """
        return saveEntry(f, entry)

    def generateOldLibraryEntry(self, data):
        """
        Return a list of values used to save entries to the old-style RMG
        transport database based on the transport object `data`.
        """
        return generateOldLibraryEntry(data)

    def processOldLibraryEntry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        transport database, returning the corresponding transport object.
        """
        return processOldLibraryEntry(data)
    
class TransportGroups(Database):
    """
    A class for working with an RMG transport group additivity database.
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
                  history=None
                  ):
        if group[0:3].upper() == 'OR{' or group[0:4].upper() == 'AND{' or group[0:7].upper() == 'NOT OR{' or group[0:8].upper() == 'NOT AND{':
            item = makeLogicNode(group)
        else:
            item = Group().fromAdjacencyList(group)
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = item,
            data = transport,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            history = history or [],
        )
    
    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the transport database to the file object `f`.
        """
        return saveEntry(f, entry)

    def generateOldLibraryEntry(self, data):
        """
        Return a list of values used to save entries to the old-style RMG
        transport database based on the transport object `data`.
        """
        
        return generateOldLibraryEntry(data)

    def processOldLibraryEntry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        transport database, returning the corresponding transport object.
        """
        return processOldLibraryEntry(data)
    
class TransportDatabase(object):
    """
    A class for working with the RMG transport database.
    """
    
    def __init__(self):
        self.depository = {}
        self.libraries = {}
        self.groups = {}
        self.libraryOrder = []
        
    def __reduce__(self):
        """
        A helper function used when pickling a TransportDatabase object.
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
        A helper function used when unpickling a TransportDatabase object.
        """
        self.depository = d['depository']
        self.libraries = d['libraries']
        self.groups = d['groups']
        self.libraryOrder = d['libraryOrder']
    
    def getTransportProperties(self, molecule):
        """
        Return the transport properties for a given :class:`Species`
        object `species`. This function first searches the loaded libraries
        in order, returning the first match found, before falling back to
        estimation via group additivity.
        """
        transport = None
        
        for label in self.libraryOrder:
            transport = self.getTransportPropertiesFromLibrary(species, self.libraries[label])
            if transport is not None:
                transport[0].comment = label
                break
        else:
            #Transport not found in any loaded libraries, so estimate
            transport = self.getTransportPropertiesViaGroupEstimates(species)
        data, library, entry = transport
        
        return data
    
    def getAllTransportProperties(self, species):
        """
        Return all possible sets of transport parameters for a given
        :class:`Species` object `species`. The hits from the libraries (in order) come first, and then the group additivity
        estimate. This method is useful for a generic search job.
        """
        transport = []
        
        # Data from libraries comes first
        for label in self.libraryOrder:
            data = self.getTransportPropertiesFromLibrary(species, self.libraries[label])
            if data: 
                data[0].comment = label
                transport.append(data)
        # Last entry is always the estimate from group additivity
        transport.append(self.getTransportPropertiesViaGroupEstimates(species))
            
        # Return all of the resulting thermo parameters
        return transport
        
    def getTransportPropertiesFromLibrary(self, species, library):
        """
        Return the set of transport properties corresponding to a given
        :class:`Species` object `species` from the specified transport
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
    
    def getTransportPropertiesViaGroupEstimates(self,molecule):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Species` object `species` by estimation using the group
        additivity values. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        """
        #Kb boltzmans constant
        Kb = 1.3806503*10^(-23)
        groupData = []
        counter = 0
        for molecule in species.molecule:
            molecule.clearLabeledAtoms()
            molecule.updateAtomTypes()
            criticalPointContribution = self.estimateTransportViaGroupAdditivity(molecule)
            groupData.Tc += criticalPointContribution.Tc
            groupData.Pc += criticalPointContribution.Pc
            groupData.Vc += criticalPointContribution.Vc
            groupData.Tb += criticalPointContribution.Tb
            groupData.structureIndex += criticalPointContribution.structureIndex
            counter += 1
            
        groupData.Tb = 198.18 + groupData.Tb
        groupData.Vc = 17.5 + groupData.Vc
        groupData.Tc = groupData.Tb/(.584 + .965(groupData.Tc) - (groupData.Tc)^2)
        groupData.Pc = 1/(.113 + .0032*counter + groupData.Pc)^2
    
        transport = transport(
                     shapeIndex = 0,
                     epsilon = .77*groupData.Tc*Kb,
                     sigma = 2.44*(groupData.Tc/groupData.Pc)^(1/3),
                     dipoleMoment = 0,
                     polarizability = 0,
                     rotrelaxcollnum = 0,
                     comment = 'group estimate',
                     )
        return transport
        
    def estimateTransportViaGroupAdditivity(self, molecule):
        """
        Return the set of transport parameters corresponding to a given
        :class:`Molecule` object `molecule` by estimation using the group
        additivity values. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        """
        # For transport estimation we need the atoms to already be sorted because we
        # iterate over them; if the order changes during the iteration then we
        # will probably not visit the right atoms, and so will get the transport wrong
        molecule.sortVertices()
        
        criticalPointContribution = CriticalPointGroupContribution(
            Tc = 0,
            Pc = 0,
            Vc = 0,
            Tb = 0,
            structureIndex = 0,
            )
        
        if sum([atom.radicalElectrons for atom in molecule.atoms]) > 0: # radical species

            # Make a copy of the structure so we don't change the original
            saturatedStruct = molecule.copy(deep=True)

            # Saturate structure by replacing all radicals with bonds to
            # hydrogen atoms
            added = {}
            for atom in saturatedStruct.atoms:
                for i in range(atom.radicalElectrons):
                    H = Atom('H')
                    bond = Bond(atom, H, 'S')
                    saturatedStruct.addAtom(H)
                    saturatedStruct.addBond(bond)
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

            # Get critical point contribution estimates for saturated form of structure
            criticalPointContribution = self.estimateSoluteViaGroupAdditivity(saturatedStruct)
            assert soluteData is not None, "Solute data of saturated {0} of molecule {1} is None!".format(saturatedStruct, molecule)
            
            # For each radical site, get radical correction
            # Only one radical site should be considered at a time; all others
            # should be saturated with hydrogen atoms
            for atom in added:

                # Remove the added hydrogen atoms and bond and restore the radical
                for H, bond in added[atom]:
                    saturatedStruct.removeBond(bond)
                    saturatedStruct.removeAtom(H)
                    atom.incrementRadical()

                saturatedStruct.updateConnectivityValues()
                
        else: # non-radical species
            # Generate estimate of critical point contribution data
            for atom in molecule.atoms:
                # Iterate over heavy (non-hydrogen) atoms
                if atom.isNonHydrogen():
                    # Get initial critical point contribution from main group database
                    try:
                        self.__addCriticalPointContribution(criticalPointContribution, self.groups['nonring'], molecule, {'*':atom})
                    
                    except KeyError:
                        logging.error("Couldn't find in nonring database:")
                        logging.error(molecule)
                        logging.error(molecule.toAdjacencyList())
                        raise
                    # Get critical point contribution for non ring groups
                    try:
                        self.__addCriticalPointContribution(criticalPointContribution, self.groups['ring'], molecule, {'*':atom})
                    except KeyError: pass
                    
    def __addCriticalPointContribution(self, criticalPointContribution, database, molecule, atom):
        """
        Determine the critical point contribution values for the atom `atom`
        in the structure `structure`, and add it to the existing criticalPointContribution
        `criticalPointContribution`.
        """
        
        node0 = database.descendTree(molecule, atom, None)

        if node0 is None:
            raise KeyError('Node not found in database.')

        # It's possible (and allowed) that items in the tree may not be in the
        # library, in which case we need to fall up the tree until we find an
        # ancestor that has an entry in the library
        node = node0
        
        while node is not None and node.data is None:
            node = node.parent
        if node is None:
            raise KeyError('Node has no parent with data in database.')
        data = node.data
        comment = node.label
        while isinstance(data, basestring) and data is not None:
            for entry in database.entries.values():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    break
        comment = '{0}({1})'.format(database.label, comment)
        
        criticalPointContribution.Tc = data.Tc
        criticalPointContribution.Pc = data.Pc
        criticalPointContribution.Vc = data.Vc
        criticalPointContribution.Tb = data.Tb
        criticalPOintContribution.structureIndex = data.structureIndex
        
        return criticalPointContribution
    
class CriticalPointGroupContribution:
    """Joback group contribution to estimate critical properties"""
    def __init__(self, Tc=None, Pc=None, Vc=None, Tb=None, structureIndex=None):
        self.Tc = Tc
        self.Pc = Pc
        self.Vc = Vc
        self.Tb = Tb
        self.structureIndex = structureIndex
        
    def _repr_(self):
        """
        Return a string representation that can be used to reconstruct the
        CriticalPointGroupContribution object
        """
        string = 'CriticalPointGroupContribution(Tc={0!r}, Pc={1!r}, Vc={2!r}, Tb={3!r}, structureIndex={4!r}'.format(self.Tc, self.Pc, self.Vc, self.Tb, self.structureIndex)
        string += ')'
        return string
    
    def __getTc__(self):
        """Returns the value of the critical temperature of the transport group."""
        return self._Tc
    def _setTc_(self,value):
        """Sets the value of the critical temperature of the transport group."""
        self.Tc = value 
        
    def __getPc__(self):
        """Returns the value of the critical pressure of the transport group."""
        return self._Pc
    def _setPc_(self,value):
        """Sets the value of the critical pressure of the transport group"""
        self.Pc = value 
            
    def __getVc__(self):
        """Returns the value of the critical volume of the transport group."""
        return self._Vc
    def _setVc_(self,value):
        """Sets the value of the critical volume of the transport group."""
        self.Vc = value
             
    def __getTb__(self):
        """Returns the value of the boiling point of the transport group."""
        return self._Tb
    def _setTb_(self,value):
        """Sets the value of the boiling point of the transport group."""
        self.Tb = value
             
    def __getstructureIndex__(self):
        """Returns the value of the structure index of the transport group."""
        return self._structureIndex
    def _setstructureIndex_(self,value):
        """Sets the value of the structure index of the transport group."""
        self.structureIndex = value