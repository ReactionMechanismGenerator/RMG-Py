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
    `software`          ``str``                 Quantum chemical package name in common letters
    `method`            ``str``                 Semi-empirical method
    `fileStore`         ``str``                 The path to the QMfiles directory
    `scratchDirectory`  ``str``                 The path to the scratch directory
    `onlyCyclics`       ``bool``                ``True`` if to run QM only on ringed species
    `maxRadicalNumber`  ``int``                 Radicals larger than this are saturated before applying HBI
    =================== ======================= ====================================
    
    """
    def __init__(self,
                 software = None,
                 method = 'pm3',
                 fileStore = None,
                 scratchDirectory = None,
                 onlyCyclics = True,
                 maxRadicalNumber = 0,
                 ):
        self.software = software
        self.method = method
        if fileStore:
            self.fileStore = os.path.join(fileStore, method)
        else:
            self.fileStore = None
        if scratchDirectory:
            self.scratchDirectory = os.path.join(scratchDirectory, method)
        else:
            self.scratchDirectory = None
        self.onlyCyclics = onlyCyclics
        self.maxRadicalNumber = maxRadicalNumber
        
        if os.sys.platform == 'win32':
            symmetryPath = os.path.join(rmgpy.getPath(),'..', 'bin', 'symmetry.exe')
            # If symmetry is not installed in the bin folder, assume it is available on the path somewhere
            if not os.path.exists(symmetryPath):
                symmetryPath = 'symmetry.exe' 
        else:
            symmetryPath = os.path.join(rmgpy.getPath(),'..', 'bin', 'symmetry')
            if not os.path.exists(symmetryPath):
                symmetryPath = 'symmetry'
        self.symmetryPath = symmetryPath

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
                 software = None,
                 method = 'pm3',
                 fileStore = None,
                 scratchDirectory = None,
                 onlyCyclics = True,
                 maxRadicalNumber = 0,
                 ):
                 
        self.settings = QMSettings(software = software,
                                   method = method,
                                   fileStore = fileStore,
                                   scratchDirectory = scratchDirectory,
                                   onlyCyclics = onlyCyclics,
                                   maxRadicalNumber = maxRadicalNumber,
                                   )
        self.database = ThermoLibrary(name='QM Thermo Library')
        
    def setDefaultOutputDirectory(self, outputDirectory):
        """
        IF the fileStore or scratchDirectory are not already set, put them in here.
        """
        
        if not self.settings.fileStore:
            self.settings.fileStore = os.path.abspath(os.path.join(outputDirectory, 'QMfiles', self.settings.method))
            logging.info("Setting the quantum mechanics fileStore to {0}".format(self.settings.fileStore))
        if not self.settings.scratchDirectory:
            self.settings.scratchDirectory = os.path.abspath(os.path.join(outputDirectory, 'QMscratch', self.settings.method))
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
                # This try/except should be redundant, but some networked file systems
                # seem to be slow or buggy or respond strangely causing problems
                # between checking the path exists and trying to create it.
                try:
                    os.makedirs(path)
                except OSError as e:
                    logging.warning("Error creating directory {0}: {1!r}".format(path, e))
                    logging.warning("Checking it already exists...")
                    assert os.path.exists(path), "Path {0} still doesn't exist?".format(path)

    def getThermoData(self, molecule):
        """
        Generate thermo data for the given :class:`Molecule` via a quantum mechanics calculation.
        
        Ignores the settings onlyCyclics and maxRadicalNumber and does the calculation anyway if asked.
        (I.e. the code that chooses whether to call this method should consider those settings).
        """
        self.initialize()
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
    
        
