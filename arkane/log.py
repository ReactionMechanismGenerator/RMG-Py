#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
A general class for parsing quantum mechanical log files
"""

class Log(object):
    """
    Represent a general log file.
    The attribute `path` refers to the location on disk of the log file of interest.
    """
    def __init__(self, path):
        self.path = path

    def getNumberOfAtoms(self):
        """
        Return the number of atoms in the molecular configuration used in
        the MolPro log file.
        """
        raise NotImplementedError("loadGeometry is not implemented for the Log class")

    def loadForceConstantMatrix(self):
        """
        Return the force constant matrix (in Cartesian coordinates) from the
        QChem log file. If multiple such matrices are identified,
        only the last is returned. The units of the returned force constants
        are J/m^2. If no force constant matrix can be found in the log file,
        ``None`` is returned.
        """
        raise NotImplementedError("loadGeometry is not implemented for the Log class")

    def loadGeometry(self):
        """
        Return the optimum geometry of the molecular configuration from the
        log file. If multiple such geometries are identified, only the
        last is returned.
        """
        raise NotImplementedError("loadGeometry is not implemented for the Log class")

    def loadConformer(self, symmetry=None, spinMultiplicity=0, opticalIsomers=None, symfromlog=None, label=''):
        """
        Load the molecular degree of freedom data from an output file created as the result of a
        QChem "Freq" calculation. As QChem's guess of the external symmetry number is not always correct,
        you can use the `symmetry` parameter to substitute your own value;
        if not provided, the value in the QChem output file will be adopted.
        """
        raise NotImplementedError("loadGeometry is not implemented for the Log class")

    def loadEnergy(self, frequencyScaleFactor=1.):
        """
        Load the energy in J/mol from a QChem log file. Only the last energy
        in the file is returned. The zero-point energy is *not* included in
        the returned value.
        """
        raise NotImplementedError("loadGeometry is not implemented for the Log class")

    def loadZeroPointEnergy(self,frequencyScaleFactor=1.):
        """
        Load the unscaled zero-point energy in J/mol from a QChem output file.
        """
        raise NotImplementedError("loadGeometry is not implemented for the Log class")

    def loadScanEnergies(self):
        """
        Extract the optimized energies in J/mol from a QChem log file, e.g. the
        result of a QChem "PES Scan" quantum chemistry calculation.
        """
        raise NotImplementedError("loadGeometry is not implemented for the Log class")

    def loadNegativeFrequency(self):
        """
        Return the imaginary frequency from a transition state frequency
        calculation in cm^-1.
        """
        raise NotImplementedError("loadGeometry is not implemented for the Log class")

