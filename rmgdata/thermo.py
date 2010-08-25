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

import os
import math
import logging
import numpy

from base import *

import chempy.constants as constants
from chempy.thermo import *
from chempy.molecule import Atom, Bond, Molecule

################################################################################

class ThermoEntry(DataEntry):
    """
    A single entry in the thermodynamics database. Each entry either contains
    a thermodynamics `model` or a string label of a `node` to look at for
    thermodynamics information.
    """

    def __init__(self, model=None, node='', index=0, label='', shortComment='', longComment='', history=None):
        DataEntry.__init__(self, index, label, shortComment, longComment, history)
        self.model = model
        self.node = node

################################################################################

class ThermoDatabase:
    """
    A base class for thermodynamics databases.
    """

    def loadDatabase(self, path):

        global currentDatabase
        currentDatabase = Database()

        global_context = { '__builtins__': None }
        local_context = {
            '__builtins__': None,
            'thermo': loadThermo,
            'tree': loadTree,
            'ThermoGAModel': loadThermoGAModel,
        }
        f = open(path)
        try:
            exec f in global_context, local_context
        except (NameError, TypeError, SyntaxError), e:
            logging.error('The input file "%s" was invalid:' % path)
            logging.exception(e)
            network = None
        finally:
            f.close()

        return currentDatabase

    def getOldDTLPaths(self, path, prefix):
        """
        Return a tuple of dictionary, tree, and library paths for a given
        prefix.
        """
        dict_path = os.path.join(path, '%s_Dictionary.txt' % prefix)
        tree_path = os.path.join(path, '%s_Tree.txt' % prefix)
        libr_path = os.path.join(path, '%s_Library.txt' % prefix)
        return dict_path, tree_path, libr_path

    def loadOldDatabase(self, dictstr, treestr, libstr, pattern=True):
        """
        Load a thermodynamics group additivity database. The database is stored
        in three files: `dictstr` is the path to the dictionary, `treestr` to
        the tree, and `libstr` to the library. The tree is optional, and should
        be set to '' if not desired.
        """

        # Load dictionary, library, and (optionally) tree
        database = Database()
        database.load(dictstr, treestr, libstr, pattern)

        # Convert data in library to ThermoGAModel objects or lists of
        # [link, comment] pairs
        for label, item in database.library.iteritems():

            if item is None:
                pass
            elif not item.__class__ is tuple:
                raise InvalidDatabaseError('Thermo library should be tuple at this point. Instead got %r'%data)
            else:
                index,item = item # break apart tuple, recover the 'index' - the beginning of the line in the library file.
                # Is't it dangerous having a local variable with the same name as a module?
                # what if we want to raise another data.InvalidDatabaseError() ?
                if not ( item.__class__ is str or item.__class__ is unicode) :
                    raise InvalidDatabaseError('Thermo library data format is unrecognized.')

                items = item.split()
                try:
                    thermoData = []; comment = ''
                    # First 12 entries are thermo data
                    for i in range(12):
                        thermoData.append(float(items[i]))
                    # Remaining entries are comment
                    for i in range(12, len(items)):
                        comment += items[i] + ' '
                    comment = comment.replace('"', '').strip()

                    database.library[label] = ThermoEntry(
                        model=self.__convertOldLibraryEntry(thermoData, comment),
                        index=int(index),
                        label=label,
                        shortComment=comment,
                    )

                except (ValueError, IndexError), e:
                    # Data represents a link to a different node that contains
                    # the data to use
                    database.library[label] = ThermoEntry(
                        node=items[0],
                        index=int(index),
                        label=label,
                        shortComment=item[len(items[0])+1:].replace('"', '').strip(),
                    )

        # Check for well-formedness
        if not database.isWellFormed():
            raise InvalidDatabaseError('Database at "%s" is not well-formed.' % (dictstr))

        #database.library.removeLinks()

        return database

    def __convertOldLibraryEntry(self, data, comment):
        """
        Process a list of numbers `data` and associated description `comment`
        generated while reading from a thermodynamic database.
        """

        if len(data) != 12:
            raise Exception('Invalid list of thermo data; should be a list of numbers of length 12.')

        H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500, dH, dS, dCp = data

        H298 = float(pq.Quantity(H298, 'kcal/mol').simplified)
        S298 = float(pq.Quantity(S298, 'cal/(mol*K)').simplified)
        Tdata = numpy.array([300, 400, 500, 600, 800, 1000, 1500], numpy.float64)
        Cpdata = list(pq.Quantity([Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500], 'cal/(mol*K)').simplified)
        Cpdata = numpy.array([float(Cp) for Cp in Cpdata], numpy.float64)
        
        return ThermoGAModel(H298=H298, S298=S298, Tdata=Tdata, Cpdata=Cpdata, comment=comment)

################################################################################

class ThermoGroupDatabase(ThermoDatabase):
    """
    A set of thermodynamics group additivity databases, consisting of a primary
    database of functional groups and a number of secondary databases to provide
    corrections for 1,5-interactions, gauche interactions, radicals, rings,
    and other functionality. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `groupDatabase`     :class:`Database`   Functional group additivity values
    `radicalDatabase`   :class:`Database`   Corrections for radical species
    `ringDatabase`      :class:`Database`   Corrections for cyclic and aromatic species
    `int15Database`     :class:`Database`   Corrections for 1,5-interactions
    `gaucheDatabase`    :class:`Database`   Corrections for gauche (1,4) interactions
    `otherDatabase`     :class:`Database`   Other corrections
    =================== =================== ====================================

    """

    def __init__(self, path=''):
        if path != '':
            self.load(path)
        else:
            self.groupDatabase = None
            self.int15Database = None
            self.gaucheDatabase = None
            self.otherDatabase = None
            self.radicalDatabase = None
            self.ringDatabase = None

    def load(self, path, old=False):
        """
        Load a set of thermodynamics group additivity databases from the general
        database specified at `datapath`.
        """

        path = os.path.abspath(path)

        logging.info('Loading group thermodynamics databases from %s...' % path)
        if old:
            self.groupDatabase = self.loadOldDatabase(*self.getOldDTLPaths(path, 'Group')) # the '*' unpacks the tuple into three separate arguments
            self.int15Database = self.loadOldDatabase(*self.getOldDTLPaths(path, '15'))
            self.gaucheDatabase = self.loadOldDatabase(*self.getOldDTLPaths(path, 'Gauche'))
            self.radicalDatabase = self.loadOldDatabase(*self.getOldDTLPaths(path, 'Radical'))
            self.ringDatabase = self.loadOldDatabase(*self.getOldDTLPaths(path, 'Ring'))
            self.otherDatabase = self.loadOldDatabase(*self.getOldDTLPaths(path, 'Other'))
        else:
            self.groupDatabase = self.loadDatabase(os.path.join(path, 'group.py'))
            self.int15Database = self.loadDatabase(os.path.join(path, 'int15.py'))
            self.gaucheDatabase = self.loadDatabase(os.path.join(path, 'gauche.py'))
            self.radicalDatabase = self.loadDatabase(os.path.join(path, 'radical.py'))
            self.ringDatabase = self.loadDatabase(os.path.join(path, 'ring.py'))
            self.otherDatabase = self.loadDatabase(os.path.join(path, 'other.py'))

        logging.info('')

    def generateThermoData(self, molecule):
        """
        Determine the group additivity thermodynamic data for the given
        `molecule`.
        """

        thermoData = None

        if sum([atom.radicalElectrons for atom in molecule.atoms]) > 0:

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
            thermoData = self.generateThermoData(saturatedStruct)
            assert thermoData is not None, "Thermo data of saturated %s of molecule %s is None!" % (saturatedStruct, molecule)
            
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
            
                thermoData += self.__getThermoData(self.radicalDatabase, saturatedStruct, {'*':atom})

                # Re-saturate
                for H, bond in added[atom]:
                    saturatedStruct.addAtom(H)
                    saturatedStruct.addBond(atom, H, bond)
                    atom.decrementRadical()

                # Subtract the enthalpy of the added hydrogens
                for H, bond in added[atom]:
                    thermoData.H298 -= 52.103 * 4184

            # Correct the entropy for the symmetry number

        else:
            # Generate estimate of thermodynamics
            for atom in molecule.atoms:
                # Iterate over heavy (non-hydrogen) atoms
                if atom.isNonHydrogen():
                    # Get initial thermo estimate from main group database
                    try:
                        if thermoData is None:
                            thermoData = self.__getThermoData(self.groupDatabase, molecule, {'*':atom})
                        else:
                            thermoData += self.__getThermoData(self.groupDatabase, molecule, {'*':atom})
                    except KeyError:
                        print molecule
                        print molecule.toAdjacencyList()
                        raise
                    # Correct for gauche and 1,5- interactions
                    try:
                        thermoData += self.__getThermoData(self.gaucheDatabase, molecule, {'*':atom})
                    except KeyError: pass
                    try:
                        thermoData += self.__getThermoData(self.int15Database, molecule, {'*':atom})
                    except KeyError: pass
                    try:
                        thermoData += self.__getThermoData(self.otherDatabase, molecule, {'*':atom})
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
                thermoData += self.__getThermoData(self.ringDatabase, ringStructure, {})

        return thermoData

    def __getThermoData(self, database, molecule, atom):
        """
        Determine the group additivity thermodynamic data for the atom `atom`
        in the structure `structure`.
        """

        node = database.descendTree(molecule, atom, None)

        if node is None:
            raise KeyError('Node not found in database.')
        else:
            data = database.library[node]

        while data.model is None and data.node is not None:
            data = database.library[data.node]

        # This code prints the hierarchy of the found node; useful for debugging
        #result = ''
        #while node is not None:
        #	result = ' -> ' + node + result
        #	node = database.tree.parent[node]
        #print result[4:]

        return data.model

################################################################################

class ThermoPrimaryDatabase(ThermoDatabase):
    """
    A primary thermodynamics databases, consisting of a dictionary of species
    and a library of corresponding thermodynamic data. (No tree is utilized in
    this database.) The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `database`          :class:`Database`   Thermodynamic data for individual species
    =================== =================== ====================================

    """

    def __init__(self, path=''):
        if path != '':
            self.load(path)
        else:
            self.database = None

    def load(self, path, old=False):
        """
        Load a primary thermodynamics database from the location `path`.
        """
        path = os.path.abspath(path)
        logging.info('Loading primary thermodynamics database from %s...' % path)
        if old:
            self.database = self.loadOldDatabase(os.path.join(path,'Dictionary.txt'), '', os.path.join(path,'Library.txt'), pattern=False)
        else:
            self.database = self.loadDatabase(os.path.join(path, 'database.py'))
        logging.info('')

    def generateThermoData(self, molecule):
        """
        Determine the group additivity thermodynamic data for the given
        `molecule`.
        """
        for node, struct in self.database.dictionary.iteritems():
            if molecule.isIsomorphic(struct):
                return self.database.library[node].model
        return None

################################################################################

# A module-level variable that stores the currently-loading database
currentDatabase = None

def loadThermo(label, index, group=None, species=None, model=None, node='', short_comment='', long_comment='', history=None):
    global currentDatabase
    if species is not None:
        currentDatabase.dictionary[label] = currentDatabase.dictionary.toStructure(species, pattern=False)
    elif group is not None:
        currentDatabase.dictionary[label] = currentDatabase.dictionary.toStructure('\n' + group, pattern=True)
    else:
        raise InvalidDatabaseError('Node "%s" must contain either a group or a species.' % label)
    currentDatabase.library[label] = ThermoEntry(model, node, index, label, short_comment, long_comment, history)
    if model is not None: model.comment = short_comment

def loadTree(string):
    global currentDatabase
    currentDatabase.tree.loadString(string)

def loadThermoGAModel(Tdata, Cpdata, H298, S298, Tmin=(0.0,"K"), Tmax=(99999.9,"K")):

    # Convert all data to SI units
    Tdata = numpy.array([float(T) for T in pq.Quantity(*Tdata).simplified], numpy.float64)
    Cpdata = numpy.array([float(C) for C in pq.Quantity(*Cpdata).simplified], numpy.float64)
    H298 = float(pq.Quantity(*H298).simplified)
    S298 = float(pq.Quantity(*S298).simplified)
    Tmin = float(pq.Quantity(*Tmin).simplified)
    Tmax = float(pq.Quantity(*Tmax).simplified)

    # Create and return the ThermoGAModel object
    return ThermoGAModel(Tdata, Cpdata, H298, S298, Tmin, Tmax)

################################################################################

thermoDatabases = []

forbiddenStructures = None

def loadThermoDatabase(dstr, group, old=False):
    """
    Load the RMG thermo database located at `dstr` into the global variable
    :data:`thermoDatabase`. Also loads the forbidden structures into
    :data:`forbiddenStructures`.
    """
    global thermoDatabases
    global forbiddenStructures

    if group:
        thermoDatabase = ThermoGroupDatabase()
    else:
        thermoDatabase = ThermoPrimaryDatabase()
    thermoDatabase.load(path=dstr, old=old)
    thermoDatabases.append(thermoDatabase)

    return thermoDatabase

################################################################################

def generateThermoData(molecule):
    """
    Get the thermodynamic data associated with `molecule` by looking in the
    loaded thermodynamic database. The parameter `thermoClass` is the class of
    thermo object you want returning; default is :class:`NASAModel`.
    """

    implicitH = molecule.implicitHydrogens
    molecule.makeHydrogensExplicit()

    # For thermo estimation we need the atoms to already be sorted because we
    # iterate over them; if the order changes during the iteration then we
    # will probably not visit the right atoms, and so will get the thermo wrong
    molecule.sortVertices()

    for thermoDatabase in thermoDatabases:
        GAthermoData = thermoDatabase.generateThermoData(molecule)
        if GAthermoData is not None and isinstance(thermoDatabase, ThermoGroupDatabase):
            # Correct entropy for symmetry number
            molecule.calculateSymmetryNumber()
            GAthermoData.S298 -= constants.R * math.log(molecule.symmetryNumber)
            break
        elif GAthermoData is not None and isinstance(thermoDatabase, ThermoPrimaryDatabase):
            break
            
    if implicitH: molecule.makeHydrogensImplicit()

    return GAthermoData

def convertThermoData(thermoData, molecule, thermoClass=NASAModel):
    """
    Convert a given set of `thermoData` to the class specified by `thermoClass`.
    Raises a :class:`TypeError` if this is not possible.
    """

    from chempy.ext.thermo_converter import convertGAtoWilhoit, convertWilhoitToNASA

    # Nothing to do if we already have the right thermo model
    if isinstance(thermoData, thermoClass):
        return thermoData

    thermoData0 = thermoData

    # Convert to WilhoitModel
    if isinstance(thermoData, ThermoGAModel) and (thermoClass == WilhoitModel or thermoClass == NASAModel):
        rotors = molecule.countInternalRotors()
        atoms = len(molecule.atoms)
        linear = molecule.isLinear()
        thermoData = convertGAtoWilhoit(thermoData, atoms, rotors, linear)

    # Convert to NASAModel
    if isinstance(thermoData, WilhoitModel) and thermoClass == NASAModel:
        thermoData = convertWilhoitToNASA(thermoData, Tmin=298.0, Tmax=6000.0, Tint=1000.0)

    # Make sure we have the right class
    if not isinstance(thermoData, thermoClass):
        raise TypeError('Unable to convert thermo model of type %s to type %s.' % (thermoData0.__class__, thermoClass))

    return thermoData
