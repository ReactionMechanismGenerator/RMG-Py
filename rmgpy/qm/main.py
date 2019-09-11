#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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

import logging
import os
from multiprocessing import Pool

import rmgpy.qm.gaussian
import rmgpy.qm.mopac
from rmgpy.data.thermo import ThermoLibrary


class QMSettings(object):
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
                 software=None,
                 method='pm3',
                 fileStore=None,
                 scratchDirectory=None,
                 onlyCyclics=True,
                 maxRadicalNumber=0,
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
            symmetryPath = os.path.join(rmgpy.getPath(), '..', 'bin', 'symmetry.exe')
            # If symmetry is not installed in the bin folder, assume it is available on the path somewhere
            if not os.path.exists(symmetryPath):
                symmetryPath = 'symmetry.exe'
        else:
            symmetryPath = os.path.join(rmgpy.getPath(), '..', 'bin', 'symmetry')
            if not os.path.exists(symmetryPath):
                symmetryPath = 'symmetry'
        self.symmetryPath = symmetryPath

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (QMSettings, (
            self.software,
            self.method,
            self.fileStore,
            self.scratchDirectory,
            self.onlyCyclics,
            self.maxRadicalNumber,
            self.symmetryPath
        )
                )

    def checkAllSet(self):
        """
        Check that all the required settings are set.
        """
        assert self.fileStore
        assert self.software
        assert self.method
        assert self.onlyCyclics is not None  # but it can be False
        assert isinstance(self.onlyCyclics, bool)
        assert self.maxRadicalNumber is not None  # but it can be 0
        assert isinstance(self.maxRadicalNumber, int)


class QMCalculator(object):
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
                 software=None,
                 method='pm3',
                 fileStore=None,
                 scratchDirectory=None,
                 onlyCyclics=True,
                 maxRadicalNumber=0,
                 ):

        self.settings = QMSettings(software=software,
                                   method=method,
                                   fileStore=fileStore,
                                   scratchDirectory=scratchDirectory,
                                   onlyCyclics=onlyCyclics,
                                   maxRadicalNumber=maxRadicalNumber,
                                   )
        self.database = ThermoLibrary(name='QM Thermo Library')

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (QMCalculator, (self.settings, self.database))

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
        self.settings.fileStore = os.path.expandvars(self.settings.fileStore)  # to allow things like $HOME or $RMGpy
        self.settings.scratchDirectory = os.path.expandvars(self.settings.scratchDirectory)
        for path in [self.settings.fileStore, self.settings.scratchDirectory]:
            if not os.path.exists(path):
                logging.info("Creating directory %s for QM files." % os.path.abspath(path))
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

    def runJobs(self, spc_list, procnum=1):
        """
        Run QM jobs for the provided species list (in parallel if requested).
        """
        mol_list = []
        for spc in spc_list:
            if spc.molecule[0].getRadicalCount() > self.settings.maxRadicalNumber:
                for molecule in spc.molecule:
                    if self.settings.onlyCyclics and molecule.isCyclic():
                        saturated_mol = molecule.copy(deep=True)
                        saturated_mol.saturate_radicals()
                        if saturated_mol not in mol_list:
                            mol_list.append(saturated_mol)
            else:
                if self.settings.onlyCyclics and spc.molecule[0].isCyclic():
                    if spc.molecule[0] not in mol_list:
                        mol_list.append(spc.molecule[0])
        if mol_list:
            # Zip arguments for use in map.
            qm_arg_list = [(self, mol) for mol in mol_list]

            if procnum == 1:
                logging.info('Writing QM files with {0} process.'.format(procnum))
                for qm_arg in qm_arg_list:
                    _write_QMfiles_star(qm_arg)
            elif procnum > 1:
                logging.info('Writing QM files with {0} processes.'.format(procnum))
                p = Pool(processes=procnum)
                p.map(_write_QMfiles_star, qm_arg_list)
                p.close()
                p.join()


def _write_QMfiles_star(args):
    """Wrapper to unpack zipped arguments for use with map"""
    return _write_QMfiles(*args)


def _write_QMfiles(quantumMechanics, mol):
    """
    If quantumMechanics is turned on thermo is calculated in parallel here.
    """
    quantumMechanics.getThermoData(mol)


def save(rmg):
    # Save the QM thermo to a library if QM was turned on
    if rmg.quantumMechanics:
        logging.info('Saving the QM generated thermo to qmThermoLibrary.py ...')
        rmg.quantumMechanics.database.save(os.path.join(rmg.outputDirectory, 'qmThermoLibrary.py'))


class QMDatabaseWriter(object):
    """
    This class listens to a RMG subject
    and saves the thermochemistry of species computed via the 
    QMTPmethods.


    A new instance of the class can be appended to a subject as follows:
    
    rmg = ...
    listener = QMDatabaseWriter()
    rmg.attach(listener)

    Whenever the subject calls the .notify() method, the
    .update() method of the listener will be called.

    To stop listening to the subject, the class can be detached
    from its subject:

    rmg.detach(listener)
    
    """

    def __init__(self):
        super(QMDatabaseWriter, self).__init__()

    def update(self, rmg):
        save(rmg)
