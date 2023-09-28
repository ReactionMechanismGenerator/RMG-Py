#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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
This script contains unit tests of the :mod:`chempy.constants` module.
"""


import math

import rmgpy.constants as constants


class TestConstants:
    """
    Contains unit tests that ensure that the physical constants are visible in
    both pure Python and compiled Cython modes, and in both cases are set to
    the appropriate values.
    """

    def test_avogadro_constant(self):
        """
        Test the value of the Avogadro constant.
        """
        Na = 6.02214179e23
        assert round(abs(constants.Na / Na - 1.0), 6) == 0, "{0} != {1}".format(constants.Na, Na)

    def test_boltzmann_constant(self):
        """
        Test the value of the Boltzmann constant.
        """
        kB = 1.3806504e-23
        assert round(abs(constants.kB / kB - 1.0), 6) == 0, "{0} != {1}".format(constants.kB, kB)

    def test_elementary_charge(self):
        """
        Test the value of the elementary charge constant.
        """
        e = 1.602176565e-19
        assert round(abs(constants.e / e - 1.0), 6) == 0, "{0} != {1}".format(constants.e, e)

    def test_gas_law_constant(self):
        """
        Test the value of the gas law constant.
        """
        R = 8.314472
        assert round(abs(constants.R / R - 1.0), 6) == 0, "{0} != {1}".format(constants.R, R)

    def test_planck_constant(self):
        """
        Test the value of the Planck constant.
        """
        h = 6.62606896e-34
        assert round(abs(constants.h / h - 1.0), 6) == 0, "{0} != {1}".format(constants.h, h)

    def test_reduced_planck_constant(self):
        """
        Test the value of the reduced Planck constant.
        """
        hbar = 1.054571726e-34
        assert round(abs(constants.hbar / hbar - 1.0), 6) == 0, "{0} != {1}".format(constants.hbar, hbar)

    def test_pi(self):
        """
        Test the value of pi.
        """
        assert round(abs(constants.pi / math.pi - 1.0), 6) == 0, "{0} != {1}".format(constants.pi, math.pi)

    def test_speed_of_light(self):
        """
        Test the value of the speed of light in a vacuum.
        """
        c = 299792458
        assert constants.c == c

    def test_electron_mass(self):
        """
        Test the value of the electron rest mass.
        """
        m_e = 9.10938291e-31
        assert round(abs(constants.m_e / m_e - 1.0), 6) == 0, "{0} != {1}".format(constants.m_e, m_e)

    def test_proton_mass(self):
        """
        Test the value of the proton rest mass.
        """
        m_p = 1.672621777e-27
        assert round(abs(constants.m_p / m_p - 1.0), 6) == 0, "{0} != {1}".format(constants.m_p, m_p)

    def test_neutron_mass(self):
        """
        Test the value of the neutron rest mass.
        """
        m_n = 1.674927351e-27
        assert round(abs(constants.m_n / m_n - 1.0), 6) == 0, "{0} != {1}".format(constants.m_n, m_n)

    def test_atomic_mass_unit(self):
        """
        Test the value of the atomic mass unit.
        """
        amu = 1.660538921e-27
        assert round(abs(constants.amu / amu - 1.0), 6) == 0, "{0} != {1}".format(constants.amu, amu)

    def test_bohr_radius(self):
        """
        Test the value of the Bohr radius.
        """
        a0 = 5.2917721092e-11
        assert round(abs(constants.a0 / a0 - 1.0), 6) == 0, "{0} != {1}".format(constants.a0, a0)

    def test_hartree_energy(self):
        """
        Test the value of the Hartree energy.
        """
        E_h = 4.35974434e-18
        assert round(abs(constants.E_h / E_h - 1.0), 6) == 0, "{0} != {1}".format(constants.E_h, E_h)
