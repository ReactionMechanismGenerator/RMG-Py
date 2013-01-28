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

import rmgpy.constants as constants
#from rmgpy.data.thermo import *
from rmgpy.molecule import Molecule, Atom, Bond, Group

################################################################################

def saveEntry(f, entry):
    """
    Write a Pythonic string representation of the given `entry` in the thermo
    database to the file object `f`.
    """
    raise NotImplementedError()

def generateOldLibraryEntry(data):
    """
    Return a list of values used to save entries to the old-style RMG
    thermo database based on the thermodynamics object `data`.
    """
    raise NotImplementedError()
    
def processOldLibraryEntry(data):
    """
    Process a list of parameters `data` as read from an old-style RMG
    thermo database, returning the corresponding thermodynamics object.
    """
    raise NotImplementedError()


class SoluteData():
    """
    Stores Abraham parameters to characterize a solute
    """
    def __init__(self, S=None, B=None, E=None, L=None, A=None, comment=""):
        #: :math:`\pi_2^H`
        self.S = S
        self.B = B
        self.E = E
        self.L = L
        self.A = A
        self.comment = comment
    def __repr__(self):
        return "SoluteData(S={0},B={1},E={2},L={3},A={4},comment={5!r})".format(self.S, self.B, self.E, self.L, self.A, self.comment)
        

################################################################################

class SoluteLibrary(Database):
    """
    A class for working with a RMG solute library.
    """
    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self,
                  index,
                  label,
                  molecule,
                  solute,
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
            data = solute,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            history = history or [],
        )

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the solute database to the file object `f`.
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

class SoluteGroups(Database):
    """
    A class for working with an RMG solute group additivity database.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self,
                  index,
                  label,
                  group,
                  solute,
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
            data = solute,
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

class SoluteDatabase(object):
    """
    A class for working with the RMG solute database.
    """

    def __init__(self):
        #self.depository = {}
        self.libraries = {}
        self.groups = {}
        self.libraryOrder = []
        self.local_context = {
            'SoluteData': SoluteData,
        }
        self.global_context = {}

    def __reduce__(self):
        """
        A helper function used when pickling a SoluteDatabase object.
        """
        d = {
            'libraries': self.libraries,
            'groups': self.groups,
            'libraryOrder': self.libraryOrder,
            }
        return (SoluteDatabase, (), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling a SoluteDatabase object.
        """
        #self.depository = d['depository']
        self.libraries = d['libraries']
        self.groups = d['groups']
        self.libraryOrder = d['libraryOrder']

    def load(self, path, libraries=None, depository=True):
        """
        Load the solute database from the given `path` on disk, where `path`
        points to the top-level folder of the solute database.
        """
        # if depository:
            # self.loadDepository(os.path.join(path, 'depository'))
        # else:
            # self.depository = {}
        #no solute library right now...
        #self.loadLibraries(os.path.join(path, 'libraries'), libraries)
        self.loadGroups(os.path.join(path, 'groups'))
        
    def loadDepository(self, path):
        """
        Load the thermo database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        raise NotImplementedError()

    def loadLibraries(self, path, libraries=None):
        """
        Load the solute database from the given `path` on disk, where `path`
        points to the top-level folder of the aolute database.
        """
        self.libraries = {}; self.libraryOrder = []
        for (root, dirs, files) in os.walk(os.path.join(path)):
            for f in files:
                name, ext = os.path.splitext(f)
                if ext.lower() == '.py' and (libraries is None or name in libraries):
                    logging.info('Loading solute library from {0} in {1}...'.format(f, root))
                    library = 	SoluteLibrary()
                    library.load(os.path.join(root, f), self.local_context, self.global_context)
                    library.label = os.path.splitext(f)[0]
                    self.libraries[library.label] = library
                    self.libraryOrder.append(library.label)
        if libraries is not None:
            self.libraryOrder = libraries

    def loadGroups(self, path):
        """
        Load the solute database from the given `path` on disk, where `path`
        points to the top-level folder of the solute database.
        """
        logging.info('Loading Platts additivity group database from {0}...'.format(path))
        self.groups = {}
        self.groups['abraham']   =   SoluteGroups(label='abraham').load(os.path.join(path, 'abraham.py'  ), self.local_context, self.global_context)
        # self.groups['gauche']  =  ThermoGroups(label='gauche').load(os.path.join(path, 'gauche.py' ), self.local_context, self.global_context)
        # self.groups['int15']   =   ThermoGroups(label='int15').load(os.path.join(path, 'int15.py'  ), self.local_context, self.global_context)
        # self.groups['ring']    =    ThermoGroups(label='ring').load(os.path.join(path, 'ring.py'   ), self.local_context, self.global_context)
        # self.groups['radical'] = ThermoGroups(label='radical').load(os.path.join(path, 'radical.py'), self.local_context, self.global_context)
        # self.groups['other']   =   ThermoGroups(label='other').load(os.path.join(path, 'other.py'  ), self.local_context, self.global_context)

    def save(self, path):
        """
        Save the solvation database to the given `path` on disk, where `path`
        points to the top-level folder of the solvation database.
        """
        path = os.path.abspath(path)
        if not os.path.exists(path): os.mkdir(path)
        #self.saveDepository(os.path.join(path, 'depository'))
        self.saveLibraries(os.path.join(path, 'libraries'))
        self.saveGroups(os.path.join(path, 'groups'))

    def saveDepository(self, path):
        """
        Save the thermo depository to the given `path` on disk, where `path`
        points to the top-level folder of the thermo depository.
        """
        raise NotImplementedError()

    def saveLibraries(self, path):
        """
        Save the solute libraries to the given `path` on disk, where `path`
        points to the top-level folder of the solute libraries.
        """
        if not os.path.exists(path): os.mkdir(path)
        for library in self.libraries.values():
            library.save(os.path.join(path, '{0}.py'.format(library.label)))

    def saveGroups(self, path):
        """
        Save the solute groups to the given `path` on disk, where `path`
        points to the top-level folder of the solute groups.
        """
        if not os.path.exists(path): os.mkdir(path)
        self.groups['abraham'].save(os.path.join(path, 'abraham.py'))
        # self.groups['gauche'].save(os.path.join(path, 'gauche.py'))
        # self.groups['int15'].save(os.path.join(path, 'int15.py'))
        # self.groups['ring'].save(os.path.join(path, 'ring.py'))
        # self.groups['radical'].save(os.path.join(path, 'radical.py'))
        # self.groups['other'].save(os.path.join(path, 'other.py'))

    def loadOld(self, path):
        """
        Load the old RMG solute database from the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        # The old database does not have a depository, so create an empty one
        # self.depository = {}
        # self.depository['stable']  = ThermoDepository(label='stable', name='Stable Molecules')
        # self.depository['radical'] = ThermoDepository(label='radical', name='Radical Molecules')
        
        for (root, dirs, files) in os.walk(os.path.join(path, 'thermo_libraries')):
            if os.path.exists(os.path.join(root, 'Dictionary.txt')) and os.path.exists(os.path.join(root, 'Library.txt')):
                library = SoluteLibrary(label=os.path.basename(root), name=os.path.basename(root))
                library.loadOld(
                    dictstr = os.path.join(root, 'Dictionary.txt'),
                    treestr = '',
                    libstr = os.path.join(root, 'Library.txt'),
                    numParameters = 5,
                    numLabels = 1,
                    pattern = False,
                )
                library.label = os.path.basename(root)
                self.libraries[library.label] = library

        self.groups = {}
        self.groups['abraham'] = SoluteGroups(label='abraham', name='Platts Group Additivity Values for Abraham Solute Descriptors').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', 'Abraham_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', 'Abraham_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', 'Abraham_Library.txt'),
            numParameters = 5,
            numLabels = 1,
            pattern = True,
        )
        # self.groups['gauche'] = ThermoGroups(label='gauche', name='Gauche Interaction Corrections').loadOld(
            # dictstr = os.path.join(path, 'thermo_groups', 'Gauche_Dictionary.txt'),
            # treestr = os.path.join(path, 'thermo_groups', 'Gauche_Tree.txt'),
            # libstr = os.path.join(path, 'thermo_groups', 'Gauche_Library.txt'),
            # numParameters = 12,
            # numLabels = 1,
            # pattern = True,
        # )
        # self.groups['int15'] = ThermoGroups(label='int15', name='1,5-Interaction Corrections').loadOld(
            # dictstr = os.path.join(path, 'thermo_groups', '15_Dictionary.txt'),
            # treestr = os.path.join(path, 'thermo_groups', '15_Tree.txt'),
            # libstr = os.path.join(path, 'thermo_groups', '15_Library.txt'),
            # numParameters = 12,
            # numLabels = 1,
            # pattern = True,
        # )
        # self.groups['radical'] = ThermoGroups(label='radical', name='Radical Corrections').loadOld(
            # dictstr = os.path.join(path, 'thermo_groups', 'Radical_Dictionary.txt'),
            # treestr = os.path.join(path, 'thermo_groups', 'Radical_Tree.txt'),
            # libstr = os.path.join(path, 'thermo_groups', 'Radical_Library.txt'),
            # numParameters = 12,
            # numLabels = 1,
            # pattern = True,
        # )
        # self.groups['ring'] = ThermoGroups(label='ring', name='Ring Corrections').loadOld(
            # dictstr = os.path.join(path, 'thermo_groups', 'Ring_Dictionary.txt'),
            # treestr = os.path.join(path, 'thermo_groups', 'Ring_Tree.txt'),
            # libstr = os.path.join(path, 'thermo_groups', 'Ring_Library.txt'),
            # numParameters = 12,
            # numLabels = 1,
            # pattern = True,
        # )
        # self.groups['other'] = ThermoGroups(label='other', name='Other Corrections').loadOld(
            # dictstr = os.path.join(path, 'thermo_groups', 'Other_Dictionary.txt'),
            # treestr = os.path.join(path, 'thermo_groups', 'Other_Tree.txt'),
            # libstr = os.path.join(path, 'thermo_groups', 'Other_Library.txt'),
            # numParameters = 12,
            # numLabels = 1,
            # pattern = True,
        # )

    def saveOld(self, path):
        """
        Save the old RMG Abraham database to the given `path` on disk, where
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
        self.groups['abraham'].saveOld(
            dictstr = os.path.join(groupsPath, 'Abraham_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Abraham_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Abraham_Library.txt'),
        )
        # self.groups['gauche'].saveOld(
            # dictstr = os.path.join(groupsPath, 'Gauche_Dictionary.txt'),
            # treestr = os.path.join(groupsPath, 'Gauche_Tree.txt'),
            # libstr = os.path.join(groupsPath, 'Gauche_Library.txt'),
        # )
        # self.groups['int15'].saveOld(
            # dictstr = os.path.join(groupsPath, '15_Dictionary.txt'),
            # treestr = os.path.join(groupsPath, '15_Tree.txt'),
            # libstr = os.path.join(groupsPath, '15_Library.txt'),
        # )
        # self.groups['radical'].saveOld(
            # dictstr = os.path.join(groupsPath, 'Radical_Dictionary.txt'),
            # treestr = os.path.join(groupsPath, 'Radical_Tree.txt'),
            # libstr = os.path.join(groupsPath, 'Radical_Library.txt'),
        # )
        # self.groups['ring'].saveOld(
            # dictstr = os.path.join(groupsPath, 'Ring_Dictionary.txt'),
            # treestr = os.path.join(groupsPath, 'Ring_Tree.txt'),
            # libstr = os.path.join(groupsPath, 'Ring_Library.txt'),
        # )
        # self.groups['other'].saveOld(
            # dictstr = os.path.join(groupsPath, 'Other_Dictionary.txt'),
            # treestr = os.path.join(groupsPath, 'Other_Tree.txt'),
            # libstr = os.path.join(groupsPath, 'Other_Library.txt'),
        # )

    def getSoluteData(self, species):
        """
        Return the solute descriptors for a given :class:`Species`
        object `species`. This function first searches the loaded libraries
        in order, returning the first match found, before falling back to
        estimation via Platts group additivity.
        """
        soluteData = None
        # Check the libraries in order first; return the first successful match
        for label in self.libraryOrder:
            soluteData = self.getSoluteDataFromLibrary(species, self.libraries[label])
            if soluteData is not None: 
                soluteData[0].comment = label
                break
        else:
            # Solute not found in any loaded libraries, so estimate
            soluteData = self.getSoluteDataFromGroups(species)
        # Add Cp0 and CpInf values
         # Cp0 = species.calculateCp0()
         # CpInf = species.calculateCpInf()
        data, library, entry = soluteData
         # if isinstance(data,SoluteData):
            # data.Cp0 = (Cp0,"J/(mol*K)")
            # data.CpInf = (CpInf,"J/(mol*K)")
        # Return the resulting solute parameters
        return data

    def getAllSoluteData(self, species):
        """
        Return all possible sets of Abraham solute descriptors for a given
        :class:`Species` object `species`. The hits from the library come
        first, then the group additivity  estimate. This method is useful 
		 for a generic search job.
        """
        thermoData = []
        # Data from depository comes first
        # thermoData.extend(self.getThermoDataFromDepository(species))
        # Data from libraries comes second
        for label in self.libraryOrder:
            data = self.getSoluteDataFromLibrary(species, self.libraries[label])
            if data: 
                data[0].comment = label
                soluteData.append(data)
        # Last entry is always the estimate from group additivity
        soluteData.append(self.getSoluteDataFromGroups(species))
		
        # Add Cp0 and CpInf values
        # Cp0 = species.calculateCp0()
        # CpInf = species.calculateCpInf()
        # for data, library, entry in thermoData:
            # if isinstance(data,ThermoData):
                # data.Cp0 = (Cp0,"J/(mol*K)")
                # data.CpInf = (CpInf,"J/(mol*K)")
				
        # Return all of the resulting thermo parameters
        return thermoData

    def getThermoDataFromDepository(self, species):
        """
        Return all possible sets of thermodynamic parameters for a given
        :class:`Species` object `species` from the depository. If no
        depository is loaded, a :class:`DatabaseError` is raised.
        """
        raise NotImplementedError()

    def getSoluteDataFromLibrary(self, species, library):
        """
        Return the set of Abraham solute descriptors corresponding to a given
        :class:`Species` object `species` from the specified solute
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

    def getSoluteDataFromGroups(self, species):
        """
        Return the set of Abraham solute parameters corresponding to a given
        :class:`Species` object `species` by estimation using the Platts group
        additivity method. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        
        It averages (linearly) over the desciptors for each Molecule (resonance isomer)
        in the Species.
        """       

        soluteData = SoluteData(0.0,0.0,0.0,0.0,0.0)
        count = 0
        comments = []
        for molecule in species.molecule:
            molecule.clearLabeledAtoms()
            molecule.updateAtomTypes()
            sdata = self.estimateSoluteViaGroupAdditivity(molecule)

            soluteData.S += sdata.S
            soluteData.B += sdata.B
            soluteData.E += sdata.E
            soluteData.L += sdata.L
            soluteData.A += sdata.A
            count += 1
            comments.append(sdata.comment)
        
        soluteData.S /= count
        soluteData.B /= count
        soluteData.E /= count
        soluteData.L /= count
        soluteData.A /= count
        soluteData.comment = "Average of {0}".format(" and ".join(comments))

        return soluteData, None, None
        
    def estimateSoluteViaGroupAdditivity(self, molecule):
        """
        Return the set of Abraham solute parameters corresponding to a given
        :class:`Molecule` object `molecule` by estimation using the Platts' group
        additivity method. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        """
        # For thermo estimation we need the atoms to already be sorted because we
        # iterate over them; if the order changes during the iteration then we
        # will probably not visit the right atoms, and so will get the thermo wrong
        molecule.sortVertices()

        # Create the SoluteData object
        soluteData = SoluteData(
            S = 0.277,
            B = 0.071,
            E = 0.248,
            L = 0.13,
            A = 0.003
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

            # Get solute descriptor estimates for saturated form of structure
            soluteData = self.estimateSoluteViaGroupAdditivity(saturatedStruct)
            assert soluteData is not None, "Solute data of saturated {0} of molecule {1} is None!".format(saturatedStruct, molecule)
            # Undo symmetry number correction for saturated structure
            # thermoData.S298.value_si += constants.R * math.log(saturatedStruct.symmetryNumber)

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
                
                
                        
                # Re-saturate
                

                # Subtract the enthalpy of the added hydrogens
            

            # Correct the entropy for the symmetry number

        else: # non-radical species
            # Generate estimate of solute data
            for atom in molecule.atoms:
                # Iterate over heavy (non-hydrogen) atoms
                if atom.isNonHydrogen():
                    # Get initial solute data from main group database
                    try:
                        self.__addGroupSoluteData(soluteData, self.groups['abraham'], molecule, {'*':atom})
                    except KeyError:
                        logging.error("Couldn't find in main abraham database:")
                        logging.error(molecule)
                        logging.error(molecule.toAdjacencyList())
                        raise
                        
                    # Correct for gauche and 1,5- interactions
                    

            # Do ring corrections separately because we only want to match
            # each ring one time; this doesn't work yet
            

                # Get thermo correction for this ring
                
                
        # Correct entropy for symmetry number
        

        return soluteData

    def __addGroupSoluteData(self, soluteData, database, molecule, atom):
        """
        Determine the Platts group additivity solute data for the atom `atom`
        in the structure `structure`, and add it to the existing solute data
        `soluteData`.
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
            raise InvalidDatabaseError('Unable to determine solute parameters for {0}: no library entries for {1} or any of its ancestors.'.format(molecule, node0) )

        data = node.data
        comment = node.label
        while isinstance(data, basestring) and data is not None:
            for entry in database.entries.values():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    break
        comment = '{0}({1})'.format(database.label, comment)

        # This code prints the hierarchy of the found node; useful for debugging
        #result = ''
        #while node is not None:
        #   result = ' -> ' + node + result
        #   node = database.tree.parent[node]
        #print result[4:]

        #if len(thermoData.Tdata.value_si) != len(data.Tdata.value_si) or any([T1 != T2 for T1, T2 in zip(thermoData.Tdata.value_si, data.Tdata.value_si)]):
            #raise ThermoError('Cannot add these ThermoData objects due to their having different temperature points.')
        
        #for i in range(7):
            #thermoData.Cpdata.value_si[i] += data.Cpdata.value_si[i]
        soluteData.S += data.S
        soluteData.B += data.B
        soluteData.E += data.E
        soluteData.L += data.L
        soluteData.A += data.A
        soluteData.comment += comment + "+"
        
        return soluteData
