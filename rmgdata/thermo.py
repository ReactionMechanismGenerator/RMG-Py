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

################################################################################

class ThermoDatabase:
    """
    A set of thermodynamics group additivity databases, consisting of a primary
    database of functional groups and a number of secondary databases to provide
    corrections for 1,5-interactions, gauche interactions, radicals, rings,
    and other functionality.
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

    def load(self, path):
        """
        Load a set of thermodynamics group additivity databases from the general
        database specified at `datapath`.
        """

        path = os.path.abspath(os.path.join(path,'thermo_groups'))
        
        logging.info('Loading thermodynamics databases from %s...' % path)
        self.groupDatabase = self.__loadDatabase(*self.__getDTLPaths(path, 'Group')) # the '*' unpacks the tuple into three separate arguments
        self.int15Database = self.__loadDatabase(*self.__getDTLPaths(path, '15'))
        self.gaucheDatabase = self.__loadDatabase(*self.__getDTLPaths(path, 'Gauche'))
        self.radicalDatabase = self.__loadDatabase(*self.__getDTLPaths(path, 'Radical'))
        self.ringDatabase = self.__loadDatabase(*self.__getDTLPaths(path, 'Ring'))
        self.otherDatabase = self.__loadDatabase(*self.__getDTLPaths(path, 'Other'))

    def __getDTLPaths(self, path, prefix):
        """
        Return a tuple of dictionary, tree, and library paths for a given
        prefix.
        """
        dict_path = os.path.join(path, '%s_Dictionary.txt' % prefix)
        tree_path = os.path.join(path, '%s_Tree.txt' % prefix)
        libr_path = os.path.join(path, '%s_Library.txt' % prefix)
        return dict_path, tree_path, libr_path

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

                    database.library[label] = self.__convertLibraryEntry(thermoData, comment)
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

    def __convertLibraryEntry(self, data, comment):
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

    def generateThermoData(self, molecule):
        """
        Determine the group additivity thermodynamic data for the given
        `molecule`.
        """

        thermoData = None

        if sum([atom.radicalElectrons for atom in molecule.atoms]) > 0:

            # Make a copy of the structure so we don't change the original
            saturatedStruct = molecule.copy()

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
                atom.radicalElectrons = 0

            # Update the atom types of the saturated structure (not sure why
            # this is necessary, because saturating with H shouldn't be
            # changing atom types, but it doesn't hurt anything and is not
            # very expensive, so will do it anyway)
            saturatedStruct.updateAtomTypes()

            # Get thermo estimate for saturated form of structure
            thermoData = self.generateThermoData(saturatedStruct)

            # For each radical site, get radical correction
            # Only one radical site should be considered at a time; all others
            # should be saturated with hydrogen atoms
            for atom in added:

                # Remove the added hydrogen atoms and bond and restore the radical
                for H, bond in added[atom]:
                    saturatedStruct.removeBond(atom, H)
                    saturatedStruct.removeAtom(H)
                    atom.radicalElectrons += 1

                thermoData += self.__getThermoData(self.saturatedStruct, molecule, {'*':atom})

                # Re-saturate
                for H, bond in added[atom]:
                    saturatedStruct.addAtom(H)
                    saturatedStruct.addBond(atom, H, bond)
                    atom.radicalElectrons -= 1

            # Subtract the enthalpy of the added hydrogens
            thermoData_H = self.primaryDatabase.library['H']
            for bond in added[atom]:
                thermoData.H298 -= thermoData_H.H298
                #thermoData.S298 -= thermoData_H.S298

            # Correct the entropy for the symmetry number

        else:
            # Generate estimate of thermodynamics
            for atom in molecule.atoms:
                # Iterate over heavy (non-hydrogen) atoms
                if atom.isNonHydrogen():
                    # Get initial thermo estimate from main group database
                    if thermoData is None:
                        thermoData = self.__getThermoData(self.groupDatabase, molecule, {'*':atom})
                    else:
                        thermoData += self.__getThermoData(self.groupDatabase, molecule, {'*':atom})
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
        
        while data.__class__ != ThermoGAModel and data is not None:
            if data[0].__class__ == str or data[0].__class__ == unicode:
                data = database.library[data[0]]

        # This code prints the hierarchy of the found node; useful for debugging
        #result = ''
        #while node is not None:
        #	result = ' -> ' + node + result
        #	node = database.tree.parent[node]
        #print result[4:]
        
        return data

################################################################################

thermoDatabase = None

forbiddenStructures = None

def loadThermoDatabase(dstr):
    """
    Load the RMG thermo database located at `dstr` into the global variable
    `rmg.species.thermoDatabase`. Also loads the forbidden structures into
    `rmg.thermo.forbiddenStructures`.
    """
    global thermoDatabase
    global forbiddenStructures

    thermoDatabase = ThermoDatabase(path=dstr)
    
    forbiddenStructures = Dictionary()
    forbiddenStructures.load(os.path.join(dstr, 'ForbiddenStructures.txt'))

    return thermoDatabase

################################################################################

def generateThermoData(molecule, thermoClass=NASAModel):
    """
    Get the thermodynamic data associated with `molecule` by looking in the
    loaded thermodynamic database. The parameter `thermoClass` is the class of
    thermo object you want returning; default is :class:`NASAModel`.
    """

    GAthermoData = thermoDatabase.generateThermoData(molecule)

    # Correct entropy for symmetry number
    struct.calculateSymmetryNumber()
    GAthermoData.S298 -= constants.R * math.log(struct.symmetryNumber)

    logging.debug('Group-additivity thermo data: %s' % GAthermoData)

    if thermoClass == ThermoGAModel:
        return GAthermoData  # return here because Wilhoit conversion not wanted

    # Convert to Wilhoit
    rotors = struct.calculateNumberOfRotors()
    atoms = len(struct.atoms())
    linear = struct.isLinear()
    WilhoitData = convertGAtoWilhoit(GAthermoData,atoms,rotors,linear)

    logging.debug('Wilhoit thermo data: %s' % WilhoitData)

    if thermoClass == WilhoitModel:
        return WilhoitData

    # Convert to NASA
    NASAthermoData = convertWilhoitToNASA(WilhoitData)

    logging.debug('NASA thermo data: %s' % NASAthermoData)

    # compute the error for the entire conversion, printing it as info or warning (if it is sufficiently high)
    rmsErr = NASAthermoData.rmsErr(GAthermoData)
    if(rmsErr > 0.35):
        logging.warning("Poor overall GA-to-NASA fit: Overall RMS error in heat capacity fit = %.3f*R." % (rmsErr))
    else:
        logging.debug("Overall RMS error in heat capacity fit = %.3f*R" % (rmsErr))

    if thermoClass == NASAModel:
        return NASAthermoData

    # Still not returned?
    raise Exception("Cannot convert themo data into class %r"%(required_class))

################################################################################

def isStructureForbidden(struct):
    """
    Return :data:`True` if the structure `struct` contains any of the subgraphs
    listed in the forbidden structures database, or :data:`False` if not.
    """

    if forbiddenStructures:
        for lbl, s in forbiddenStructures.iteritems():
            if struct.isSubgraphIsomorphic(s):
                return True

    return False
