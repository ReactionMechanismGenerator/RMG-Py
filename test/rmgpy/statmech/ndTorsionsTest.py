#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
This script contains unit tests of the :mod:`arkane.multidimensionalTorsions` module.
"""

import os
import zipfile
import shutil

import rmgpy
from rmgpy.statmech.ndTorsions import HinderedRotor2D, HinderedRotorClassicalND
from arkane.ess.factory import ess_factory
import numpy as np
import pytest

RMG_PATH = os.path.abspath(os.path.dirname(os.path.dirname(rmgpy.__file__)))
Q2DTOR_PATH = os.path.join(RMG_PATH, "external", "Q2DTor", "src", "Q2DTor.py")


@pytest.mark.skipif(not os.path.isfile(Q2DTOR_PATH), "Q2DTor not installed")
class TestHinderedRotor2D:
    """
    Contains unit tests of the StatmechJob class.
    """

    @classmethod
    def setup_class(cls):
        cls.path = os.path.join(RMG_PATH, "arkane", "data", "CH2CHOOH", "CH2CHOOHscans")
        if not os.path.exists(cls.path):
            zippath = os.path.join(RMG_PATH, "arkane", "data", "CH2CHOOH", "CH2CHOOHscans.zip")
            with zipfile.ZipFile(zippath, "r") as zip_ref:
                zip_ref.extractall(os.path.dirname(cls.path))

        cls.hd2d = HinderedRotor2D(
            calc_path=cls.path, name="r0", torsigma1=1, torsigma2=1, symmetry="b", pivots1=[6, 7], pivots2=[1, 6], top1=[7, 8], top2=[6, 7, 8]
        )

    def test_q2dtor_setup(self):
        self.hd2d.read_scan()
        self.assertAlmostEquals(self.hd2d.Es[0] / 10**9, -594373977.268 / 10**9, 3)
        self.hd2d.get_torsions()
        assert np.testing.assert_array_equal(self.hd2d.torsion1, [2, 1, 6, 7]) is None
        self.hd2d.write_inp()
        self.hd2d.write_pes()
        self.hd2d.get_ics_file()

    def test_partition_function_calc(self):
        self.hd2d.read_eigvals()
        self.assertAlmostEqual(self.hd2d.get_partition_function(300.0), 3.29752, 4)

    @classmethod
    def teardown_class(cls):
        """A function that is run ONCE after all unit tests in this class."""
        if os.path.exists(cls.path):
            shutil.rmtree(cls.path)  # delete unzipped and created files
        if os.path.exists(os.path.join(os.path.dirname(cls.path), "r0", "IOfiles", "r0.pes")):
            os.remove(os.path.join(os.path.dirname(cls.path), "r0", "IOfiles", "r0.pes"))
        if os.path.exists(os.path.join(os.path.dirname(cls.path), "r0", "r0.out")):
            os.remove(os.path.join(os.path.dirname(cls.path), "r0", "r0.out"))


class TestHinderedRotorClassicalND:
    """
    Contains unit tests of the HinderedRotorClassicalND class.
    """

    def test_hindered_rotor_nd(self):
        freqpath = os.path.join(RMG_PATH, "arkane", "data", "TolueneFreq.log")
        rotpath = os.path.join(RMG_PATH, "arkane", "data", "TolueneRot1.log")
        log = ess_factory(freqpath)

        conf, unscaled_freqs = log.load_conformer(symmetry=1, spin_multiplicity=1, optical_isomers=1, label="Toulene")
        coordinates, number, mass = log.load_geometry()
        conf.coordinates = (coordinates, "angstroms")
        conf.number = number
        conf.mass = (mass, "amu")

        hessian = log.load_force_constant_matrix()

        hdnd = HinderedRotorClassicalND(
            pivots=[[3, 12]], tops=[[12, 13, 14, 15]], sigmas=[6.0], calc_path=rotpath, conformer=conf, F=hessian, semiclassical=True
        )
        hdnd.read_scan()
        assert round(abs(self.hdnd.Es[0]) - 8.58538448, 4) == 0
        hdnd.fit()
        assert round(abs(self.hdnd.calc_partition_function(300.0)) - 2.899287634962152, 5) == 0
