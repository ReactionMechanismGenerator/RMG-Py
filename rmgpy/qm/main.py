#!/usr/bin/python
# -*- coding: utf-8 -*-
################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2012 Prof. Richard H. West (r.west@neu.edu),
#                           Prof. William H. Green (whgreen@mit.edu)
#                           and the RMG Team (rmg_dev@mit.edu)
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

import logging

import rmgpy.qm.mopac
import rmgpy.qm.gaussian
from rmgpy.data.thermo import ThermoLibrary

class QMSettings():
    """
    A minimal class to store settings related to quantum mechanics calculations.
    
    =================== ======================= ====================================
    Attribute           Type                    Description
    =================== ======================= ====================================
    `fileStore`         ``str``                 The path to the QMfiles directory
    `scratchDirectory`  ``str``                 The path to the scratch directory
    =================== ======================= ====================================
    
    """
    def __init__(self):
        self.software = None
        self.method = None
        self.fileStore = None
        self.scratchDirectory = None
        self.onlyCyclics = None
        self.maxRadicalNumber = None
        
        RMGpy_path = os.getenv('RMGpy') or os.path.normpath(os.path.join(rmgpy.getPath(),'..'))
        self.RMG_bin_path = os.path.join(RMGpy_path, 'bin')
    
    def checkAllSet(self):
        """
        Check that all the required settings are set.
        """
        from types import BooleanType, IntType
        assert self.fileStore
        #assert self.scratchDirectory
        assert self.software
        assert self.method
        assert self.onlyCyclics is not None # but it can be False
        assert type(self.onlyCyclics) is BooleanType
        assert self.maxRadicalNumber is not None # but it can be 0
        assert type(self.maxRadicalNumber) is IntType

class QMCalculator():
    """
    A Quantum Mechanics calculator object, to store settings. 
    
    The attributes are:

    =================== ======================= ====================================
    Attribute           Type                    Description
    =================== ======================= ====================================
    `settings`          :class:`QMSettings`     Settings for QM calculations
    `database`          :class:`ThermoLibrary`  Database containing QM calculations
    =================== ======================= ====================================

    """
    
    def __init__(self,
                 fileStore = None,
                 scratchDirectory = None,
                 ):
        self.settings = QMSettings()
        self.settings.fileStore = fileStore
        self.settings.scratchDirectory = scratchDirectory
        self.database = ThermoLibrary(name='QM Thermo Library')
        
    def setDefaultOutputDirectory(self, outputDirectory):
        """
        IF the fileStore or scratchDirectory are not already set, put them in here.
        """
        if not self.settings.fileStore:
            self.settings.fileStore = os.path.join(outputDirectory, 'QMfiles')
            logging.info("Setting the quantum mechanics fileStore to {0}".format(self.settings.fileStore))
        if not self.settings.scratchDirectory:
            self.settings.scratchDirectory = os.path.join(outputDirectory, 'QMscratch')
            logging.info("Setting the quantum mechanics scratchDirectory to {0}".format(self.settings.scratchDirectory))
    
    def initialize(self):
        """
        Do any startup tasks.
        """
        self.checkReady()

    def checkReady(self):
        """
        Check that it's ready to run calculations.
        """
        self.settings.checkAllSet()
        self.checkPaths()

    def checkPaths(self):
        """
        Check the paths in the settings are OK. Make folders as necessary.
        """
        self.settings.fileStore = os.path.expandvars(self.settings.fileStore) # to allow things like $HOME or $RMGpy
        self.settings.scratchDirectory = os.path.expandvars(self.settings.scratchDirectory)
        for path in [self.settings.fileStore, self.settings.scratchDirectory]:
            if not os.path.exists(path):
                logging.info("Creating directory %s for QM files."%os.path.abspath(path))
                os.makedirs(path)

        if not os.path.exists(self.settings.RMG_bin_path):
            raise Exception("RMG-Py 'bin' directory {0} does not exist.".format(self.settings.RMG_bin_path))
        if not os.path.isdir(self.settings.RMG_bin_path):
            raise Exception("RMG-Py 'bin' directory {0} is not a directory.".format(self.settings.RMG_bin_path))
            
        
    def getThermoData(self, molecule):
        """
        Generate thermo data for the given :class:`Molecule` via a quantum mechanics calculation.
        
        Ignores the settings onlyCyclics and maxRadicalNumber and does the calculation anyway if asked.
        (I.e. the code that chooses whether to call this method should consider those settings).
        """
        if self.settings.software == 'mopac':
            if self.settings.method == 'pm3':
                qm_molecule_calculator = rmgpy.qm.mopac.MopacMolPM3(molecule, self.settings)
            elif self.settings.method == 'pm6':
                qm_molecule_calculator = rmgpy.qm.mopac.MopacMolPM6(molecule, self.settings)
            elif self.settings.method == 'pm7':
                qm_molecule_calculator = rmgpy.qm.mopac.MopacMolPM7(molecule, self.settings)
            else:
                raise Exception("Unknown QM method '{0}' for mopac".format(self.settings.method))
            thermo0 = qm_molecule_calculator.generateThermoData()
        elif self.settings.software == 'gaussian':
            if self.settings.method == 'pm3':
                qm_molecule_calculator = rmgpy.qm.gaussian.GaussianMolPM3(molecule, self.settings)
            elif self.settings.method == 'pm6':
                qm_molecule_calculator = rmgpy.qm.gaussian.GaussianMolPM6(molecule, self.settings)
            else:
                raise Exception("Unknown QM method '{0}' for gaussian".format(self.settings.method))
            thermo0 = qm_molecule_calculator.generateThermoData()
        else:
            raise Exception("Unknown QM software '{0}'".format(self.settings.software))
        return thermo0
    
        