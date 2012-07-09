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
This module contains the :class:`RMGDatabase` class, which is the primary class
for working with the RMG database.
"""

import os
import os.path

from base import ForbiddenStructures
from thermo import ThermoDatabase
from kinetics import KineticsDatabase
from states import StatesDatabase

# Module-level variable to store the (only) instance of RMGDatabase in use.
database = None

################################################################################

class RMGDatabase:
    """
    The primary class for working with the RMG database.
    """

    def __init__(self):
        self.thermo = None
        self.forbiddenStructures = None
        self.kinetics = None
        self.states = None
        
        # Store the newly created database in the module.
        global database
        assert database is None, "Should only make one instance of RMGDatabase because it's stored as a module-level variable."
        database = self

    def load(self,
             path,
             thermoLibraries=None,
             reactionLibraries=None,
             seedMechanisms=None,
             kineticsFamilies=None,
             kineticsDepositories=None,
             statesLibraries=None,
             depository=True
             ):
        """
        Load the RMG database from the given `path` on disk, where `path`
        points to the top-level folder of the RMG database. If none of the
        optional arguments are provided, then the entire database will be
        loaded. You can use the optional arguments to specify that only certain
        components of the database be loaded.
        """
        self.loadThermo(os.path.join(path, 'thermo'), thermoLibraries, depository)
        self.loadForbiddenStructures(os.path.join(path, 'forbiddenStructures.py'))
        self.loadKinetics(os.path.join(path, 'kinetics'),
                          reactionLibraries,
                          seedMechanisms,
                          kineticsFamilies,
                          kineticsDepositories
                          )
        self.loadStates(os.path.join(path, 'states'), statesLibraries, depository)

    def loadThermo(self, path, thermoLibraries=None, depository=True):
        """
        Load the RMG thermo database from the given `path` on disk, where
        `path` points to the top-level folder of the RMG thermo database.
        """
        self.thermo = ThermoDatabase()
        self.thermo.load(path, thermoLibraries, depository)

    def loadForbiddenStructures(self, path):
        """
        Load the RMG forbidden structures from the given `path` on disk, where
        `path` points to the forbidden structures file.
        """
        self.forbiddenStructures = ForbiddenStructures()
        self.forbiddenStructures.load(path)

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
        if reactionLibraries is not None and seedMechanisms is not None:
            kineticsLibraries.extend(seedMechanisms)
            kineticsLibraries.extend(reactionLibraries)
        elif reactionLibraries is not None and seedMechanisms is None:
            kineticsLibraries.extend(reactionLibraries)
        elif reactionLibraries is None and seedMechanisms is not None:
            kineticsLibraries.extend(seedMechanisms)
        else:
            kineticsLibraries = None
            
        self.kinetics = KineticsDatabase()
        self.kinetics.load(path,
                           families=kineticsFamilies,
                           libraries=kineticsLibraries,
                           depositories=kineticsDepositories
                           )

    def loadStates(self, path, statesLibraries=None, depository=True):
        """
        Load the RMG states database from the given `path` on disk, where
        `path` points to the top-level folder of the RMG states database.
        """
        self.states = StatesDatabase()
        self.states.load(path, statesLibraries, depository)

    def loadOld(self, path):
        """
        Load the old RMG database from the given `path` on disk, where `path`
        points to the top-level folder of the old RMG database.
        """
        self.thermo = ThermoDatabase()
        self.thermo.loadOld(path)
        self.forbiddenStructures = ForbiddenStructures()
        self.forbiddenStructures.loadOld(os.path.join(path, 'ForbiddenStructures.txt'))
        self.kinetics = KineticsDatabase()
        self.kinetics.loadOld(path)
        self.states = StatesDatabase()
        self.states.loadOld(path)

    def save(self, path):
        """
        Save the RMG database to the given `path` on disk.
        """
        if not os.path.exists(path): os.mkdir(path)
        self.forbiddenStructures.save(os.path.join(path, 'forbiddenStructures.py'))
        self.thermo.save(os.path.join(path, 'thermo'))
        self.kinetics.save(os.path.join(path, 'kinetics'))
        self.states.save(os.path.join(path, 'states'))

    def saveOld(self, path):
        """
        Save the old RMG database to the given `path` on disk.
        """
        if not os.path.exists(path): os.mkdir(path)
        self.thermo.saveOld(path)
        self.forbiddenStructures.saveOld(os.path.join(path, 'ForbiddenStructures.txt'))
        self.kinetics.saveOld(path)
        self.states.saveOld(path)
