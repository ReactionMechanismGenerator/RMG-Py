#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
This module contains the :class:`RMGDatabase` class, which is the primary class
for working with the RMG database.
"""

import os.path
import logging

from base import ForbiddenStructures
from thermo import ThermoDatabase
from transport import TransportDatabase
from rmgpy.data.kinetics.database import KineticsDatabase
from statmech import StatmechDatabase
from solvation import SolvationDatabase
from rmgpy.exceptions import DatabaseError
from rmgpy.scoop_framework.util import get, broadcast

# Module-level variable to store the (only) instance of RMGDatabase in use.
database = None

################################################################################

class RMGDatabase:
    """
    The primary class for working with the RMG database.
    """

    def __init__(self):
        self.thermo = None
        self.transport = None
        self.forbiddenStructures = None
        self.kinetics = None
        self.statmech = None
        self.solvation = None
        
        # Store the newly created database in the module.
        global database
#        assert database is None, "Should only make one instance of RMGDatabase because it's stored as a module-level variable."
        if database is None:
            database = self
        else:
            logging.warning("Should only make one instance of RMGDatabase because it's stored as a module-level variable!")
            logging.warning("Unexpected behaviour may result!")

    def load(self,
             path,
             thermoLibraries=None,
             transportLibraries=None,
             reactionLibraries=None,
             seedMechanisms=None,
             kineticsFamilies=None,
             kineticsDepositories=None,
             statmechLibraries=None,
             depository=True,
             solvation=True,
             testing = False):
        """
        Load the RMG database from the given `path` on disk, where `path`
        points to the top-level folder of the RMG database. If none of the
        optional arguments are provided, then the entire database will be
        loaded. You can use the optional arguments to specify that only certain
        components of the database be loaded.

        Argument testing will load a lighter version of the database used for unit-tests
        """
        self.loadThermo(os.path.join(path, 'thermo'), thermoLibraries, depository)
        if not testing:
            self.loadTransport(os.path.join(path, 'transport'), transportLibraries)
            self.loadForbiddenStructures(os.path.join(path, 'forbiddenStructures.py'))
        self.loadKinetics(os.path.join(path, 'kinetics'),
                          reactionLibraries,
                          seedMechanisms,
                          kineticsFamilies,
                          kineticsDepositories
                          )
        if not testing:
            self.loadStatmech(os.path.join(path, 'statmech'), statmechLibraries, depository)
        
        if solvation:
            self.loadSolvation(os.path.join(path, 'solvation'))

    def loadThermo(self, path, thermoLibraries=None, depository=True):
        """
        Load the RMG thermo database from the given `path` on disk, where
        `path` points to the top-level folder of the RMG thermo database.
        """
        self.thermo = ThermoDatabase()
        self.thermo.load(path, thermoLibraries, depository)
        broadcast(self.thermo, 'thermo')

    def loadTransport(self, path, transportLibraries=None):
        """
        Load the RMG transport database from the given 'path' on disk, where 
        'path' points to the top-level folder of the RMG transport database.
        """
        self.transport = TransportDatabase()
        self.transport.load(path, transportLibraries)
        broadcast(self.transport, 'transport')
        
    def loadForbiddenStructures(self, path = None):
        """
        Load the RMG forbidden structures from the given `path` on disk, where
        `path` points to the forbidden structures file.

        If no path is given, a blank forbidden structures object is created.
        """
        self.forbiddenStructures = ForbiddenStructures()
        if path is not None:
            self.forbiddenStructures.load(path)
        broadcast(self.forbiddenStructures, 'forbidden')

    def loadKinetics(self,
                     path,
                     reactionLibraries=None,
                     seedMechanisms=None,
                     kineticsFamilies=None,
                     kineticsDepositories=None
                     ):
        """
        Load the RMG kinetics database from the given `path` on disk, where
        `path` points to the top-level folder of the RMG kinetics database.
        """
        kineticsLibraries = []
        libraryOrder = []
        if seedMechanisms is None and reactionLibraries is None:
            kineticsLibraries = None
        if seedMechanisms is not None:
            for library in seedMechanisms:
                kineticsLibraries.append(library)
                libraryOrder.append((library,'Seed'))
        if reactionLibraries is not None:
            for library in reactionLibraries:
                kineticsLibraries.append(library)
                libraryOrder.append((library,'Reaction Library'))
        
        self.kinetics = KineticsDatabase()
        self.kinetics.libraryOrder = libraryOrder
        self.kinetics.load(path,
                           families=kineticsFamilies,
                           libraries=kineticsLibraries,
                           depositories=kineticsDepositories
                           )

        broadcast(self.kinetics, 'kinetics')

    def loadSolvation(self, path):
        """
        Load the RMG solvation database from the given `path` on disk, where
        `path` points to the top-level folder of the RMG solvation database.
        """
        self.solvation = SolvationDatabase()
        self.solvation.load(path)
        broadcast(self.solvation, 'solvation')
        
    def loadStatmech(self, path, statmechLibraries=None, depository=True):
        """
        Load the RMG statmech database from the given `path` on disk, where
        `path` points to the top-level folder of the RMG statmech database.
        """
        self.statmech = StatmechDatabase()
        self.statmech.load(path, statmechLibraries, depository)
        broadcast(self.statmech, 'statmech')

    def loadOld(self, path):
        """
        Load the old RMG database from the given `path` on disk, where `path`
        points to the top-level folder of the old RMG database.
        """
        self.thermo = ThermoDatabase()
        self.thermo.loadOld(path)
        self.transport = TransportDatabase()
        #self.transport.loadOld(path)   #  Currently no loadOld import function available for transport groups
        self.forbiddenStructures = ForbiddenStructures() 
        self.forbiddenStructures.loadOld(os.path.join(path, 'ForbiddenStructures.txt'))
        self.kinetics = KineticsDatabase()
        self.kinetics.loadOld(path)
        self.statmech = StatmechDatabase()
        self.statmech.loadOld(path)        
        self.solvation = SolvationDatabase()
        # Not completely implemented
        # self.solvation.loadOld(path)    

    def save(self, path):
        """
        Save the RMG database to the given `path` on disk.
        """
        if not os.path.exists(path): os.makedirs(path)
        self.forbiddenStructures.save(os.path.join(path, 'forbiddenStructures.py'))
        self.thermo.save(os.path.join(path, 'thermo'))
#         self.transport.save(os.path.join(path, 'transport')) #Currently no function for saving transport groups
        self.kinetics.save(os.path.join(path, 'kinetics'))
        self.statmech.save(os.path.join(path, 'statmech'))
        self.solvation.save(os.path.join(path, 'solvation'))
        self.transport.save(os.path.join(path, 'transport'))

    def saveOld(self, path):
        """
        Save the old RMG database to the given `path` on disk.
        """
        if not os.path.exists(path): os.makedirs(path)
        self.thermo.saveOld(path)
        self.transport.saveOld(path)
        self.forbiddenStructures.saveOld(os.path.join(path, 'ForbiddenStructures.txt'))
        self.kinetics.saveOld(path)
        self.statmech.saveOld(path)

def getDB(name=''):
    """
    Returns the RMG database object that corresponds
    to the parameter name.

    First, the module level is queried. If this variable
    is empty, the broadcasted variables are queried.
    """
    global database

    if database:
        if name == '':
            return database
        elif name == 'kinetics':
            return database.kinetics
        elif name == 'thermo':
            return database.thermo
        elif name == 'transport':
            return database.transport
        elif name == 'solvation':
            return database.solvation
        elif name == 'statmech':
            return database.statmech
        elif name == 'forbidden':
            return database.forbiddenStructures
        else:
            raise Exception('Unrecognized database keyword: {}'.format(name))
    else:
        try:
            db = get(name)
            if db:
                return db
            else:
                raise DatabaseError
        except DatabaseError, e:
            logging.debug("Did not find a way to obtain the broadcasted database for {}.".format(name))
            raise e

    raise DatabaseError('Could not get database with name: {}'.format(name))
