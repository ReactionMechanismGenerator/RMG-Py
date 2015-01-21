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

from rmgpy.qm.molecule import QMMolecule
from rmgpy.qm.reaction import QMReaction
import rmgpy.qm.mopac
import rmgpy.qm.gaussian
import rmgpy.qm.nwchem
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
                 method = None,
                 fileStore = None,
                 scratchDirectory = None,
                 onlyCyclics = True,
                 maxRadicalNumber = 0,
                 ):
        self.software = software
        self.method = method
        self.fileStore = fileStore
        self.scratchDirectory = scratchDirectory
        self.onlyCyclics = onlyCyclics
        self.maxRadicalNumber = maxRadicalNumber
        
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
                 software = None,
                 method = None,
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
        elif self.settings.software == 'gaussian':
            if self.settings.method == 'pm3':
                qm_molecule_calculator = rmgpy.qm.gaussian.GaussianMolPM3(molecule, self.settings)
            elif self.settings.method == 'pm6':
                qm_molecule_calculator = rmgpy.qm.gaussian.GaussianMolPM6(molecule, self.settings)
            elif self.settings.method == 'b3lyp':
                qm_molecule_calculator = rmgpy.qm.gaussian.GaussianMolB3LYP(molecule, self.settings)
            else:
                raise Exception("Unknown QM method '{0}' for gaussian".format(self.settings.method))
        elif self.settings.software == 'nwchem':
            if self.settings.method =='hf':
                qm_molecule_calculator = rmgpy.qm.nwchem.NWChemMolHF(molecule, self.settings)
            else:
                raise Exception("Unknown QM method '{0}' for nwchem".format(self.settings.method))
        else:
            raise Exception("Unknown QM software '{0}'".format(self.settings.software))
        thermo0 = qm_molecule_calculator.generateThermoData()
        return thermo0
    
    def getKineticData(self, reaction, tsDatabase):
        """
        Generate thermo data for the given :class:`Molecule` via a quantum mechanics calculation.
        
        Ignores the settings onlyCyclics and maxRadicalNumber and does the calculation anyway if asked.
        (I.e. the code that chooses whether to call this method should consider those settings).
        """
        if self.settings.software == 'mopac':
            if self.settings.method == 'pm3':
                qm_reaction_calculator = rmgpy.qm.mopac.MopacTSPM3(reaction, self.settings, tsDatabase)
            elif self.settings.method == 'pm6':
                qm_reaction_calculator = rmgpy.qm.mopac.MopacTSPM6(reaction, self.settings, tsDatabase)
            elif self.settings.method == 'pm7':
                qm_reaction_calculator = rmgpy.qm.mopac.MopacTSPM7(reaction, self.settings, tsDatabase)
            else:
                raise Exception("Unknown QM method '{0}' for mopac".format(self.settings.method))
        if self.settings.software == 'gaussian':
            if self.settings.method == 'pm6':
                qm_reaction_calculator = rmgpy.qm.gaussian.GaussianTSPM6(reaction, self.settings, tsDatabase)
            elif self.settings.method == 'b3lyp':
                qm_reaction_calculator = rmgpy.qm.gaussian.GaussianTSB3LYP(reaction, self.settings, tsDatabase)
            elif self.settings.method == 'm062x':
                qm_reaction_calculator = rmgpy.qm.gaussian.GaussianTSM062X(reaction, self.settings, tsDatabase)
            else:
                raise Exception("Unknown QM method '{0}' for gaussian".format(self.settings.method, tsDatabase))
        elif self.settings.software == 'nwchem':
            if self.settings.method =='hf':
                qm_reaction_calculator = rmgpy.qm.nwchem.NWChemTSHF(reaction, self.settings, tsDatabase)
            else:
                raise Exception("Unknown QM method '{0}' for nwchem".format(self.settings.method))
        else:
            raise Exception("Unknown QM software '{0}'".format(self.settings.software))
        
        kinetics0 = qm_reaction_calculator.generateKineticData()
        return kinetics0
    
        