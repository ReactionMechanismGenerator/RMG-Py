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
This script contains unit tests for the :mod:`arkane.modelchem` module.
"""

import unittest
from dataclasses import FrozenInstanceError

from arkane.modelchem import (LOT, LevelOfTheory, CompositeLevelOfTheory,
                              model_chem_to_lot, str_to_lot, get_software_id)

# Instances for use in tests
FREQ = LevelOfTheory(
    method='wB97X-D',
    basis='def2-TZVP',
    software='Gaussian 16',
    args='very-tight'
)
ENERGY = LevelOfTheory(
    method='DLPNO-CCSD(T)-F12',
    basis='def2-TZVP',
    software='Orca'
)
COMPOSITE = CompositeLevelOfTheory(
    freq=FREQ,
    energy=ENERGY
)

# Representations corresponding to instances
FREQ_REPR = "LevelOfTheory(method='wb97xd',basis='def2tzvp',software='gaussian',args=('verytight',))"
ENERGY_REPR = "LevelOfTheory(method='dlpnoccsd(t)f12',basis='def2tzvp',software='orca')"
COMPOSITE_REPR = f"CompositeLevelOfTheory(freq={FREQ_REPR},energy={ENERGY_REPR})"

# Dictionaries corresponding to instances
FREQ_DICT = {
    'class': 'LevelOfTheory',
    'method': 'wb97xd',
    'basis': 'def2tzvp',
    'software': 'gaussian',
    'args': ['verytight']  # This is a list instead of tuple because that's what YAML files expect
}
ENERGY_DICT = {
    'class': 'LevelOfTheory',
    'method': 'dlpnoccsd(t)f12',
    'basis': 'def2tzvp',
    'software': 'orca',
}
COMPOSITE_DICT = {
    'class': 'CompositeLevelOfTheory',
    'freq': FREQ_DICT,
    'energy': ENERGY_DICT
}

# Model chemistries corresponding to instances
FREQ_MODELCHEM = 'wb97xd/def2tzvp'
ENERGY_MODELCHEM = 'dlpnoccsd(t)f12/def2tzvp'
COMPOSITE_MODELCHEM = f'{ENERGY_MODELCHEM}//{FREQ_MODELCHEM}'


class TestLevelOfTheory(unittest.TestCase):
    """
    A class for testing that the LevelOfTheory class functions properly.
    """

    def test_attrs(self):
        """
        Test that instance behaves correctly.
        """
        self.assertEqual(FREQ.method, 'wb97xd')
        self.assertEqual(FREQ.basis, 'def2tzvp')
        self.assertEqual(FREQ.software, 'gaussian')
        self.assertTupleEqual(FREQ.args, ('verytight',))
        with self.assertRaises(FrozenInstanceError):
            FREQ.method = ''

        self.assertEqual(repr(FREQ), FREQ_REPR)
        self.assertEqual(repr(ENERGY), ENERGY_REPR)

        with self.assertRaises(ValueError):
            _ = LevelOfTheory(method=FREQ.method)
        lot = LevelOfTheory(method=FREQ.method, software=FREQ.software)
        self.assertIsNone(lot.basis)
        self.assertIsNone(lot.auxiliary_basis)
        self.assertIsNone(lot.cabs)
        self.assertIsNone(lot.software_version)
        self.assertIsNone(lot.solvent)
        self.assertIsNone(lot.solvation_method)
        self.assertIsNone(lot.args)

        self.assertIsInstance(FREQ, LOT)

    def test_comparison(self):
        """
        Test comparisons between instances.
        """
        self.assertIsInstance(hash(FREQ), int)
        self.assertNotEqual(FREQ, ENERGY)
        with self.assertRaises(TypeError):
            _ = ENERGY > FREQ

        # Test args in different order
        lot1 = LevelOfTheory('method', args=('arg1', 'arg2'))
        lot2 = LevelOfTheory('method', args=('arg2', 'arg1'))
        self.assertEqual(lot1, lot2)

    def test_simple(self):
        """
        Test that simple level of theory can be obtained.
        """
        lot = FREQ.simple()
        self.assertIsNot(lot, FREQ)
        self.assertEqual(lot.method, FREQ.method)
        self.assertEqual(lot.basis, FREQ.basis)
        self.assertEqual(lot.software, FREQ.software)
        for attr, val in lot.__dict__.items():
            if attr not in {'method', 'basis', 'software'}:
                self.assertIsNone(val)

    def test_to_model_chem(self):
        """
        Test conversion to model chemistry.
        """
        self.assertEqual(FREQ.to_model_chem(), FREQ_MODELCHEM)
        self.assertEqual(ENERGY.to_model_chem(), ENERGY_MODELCHEM)

        lot = LevelOfTheory(
            method='CBS-QB3',
            software='g16'
        )
        self.assertEqual(lot.to_model_chem(), 'cbsqb3')

    def test_update(self):
        """
        Test updating attributes.
        """
        lot = FREQ.update(software='Q-Chem')
        self.assertIsNot(lot, FREQ)
        self.assertEqual(lot.software, 'qchem')
        with self.assertRaises(TypeError):
            FREQ.update(test='test')

    def test_as_dict(self):
        """
        Test conversion to dictionary.
        """
        self.assertDictEqual(FREQ.as_dict(), FREQ_DICT)
        self.assertDictEqual(ENERGY.as_dict(), ENERGY_DICT)


class TestCompositeLevelOfTheory(unittest.TestCase):
    """
    A class for testing that the CompositeLevelOfTheory class functions properly.
    """

    def test_attrs(self):
        """
        Test that instance behaves correctly.
        """
        self.assertIs(COMPOSITE.freq, FREQ)
        self.assertIs(COMPOSITE.energy, ENERGY)
        self.assertEqual(repr(COMPOSITE), COMPOSITE_REPR)
        with self.assertRaises(FrozenInstanceError):
            COMPOSITE.energy = ''

        self.assertIsInstance(COMPOSITE, LOT)

    def test_comparison(self):
        """
        Test comparisons between instances.
        """
        other = CompositeLevelOfTheory(freq=ENERGY, energy=FREQ)
        self.assertIsInstance(hash(COMPOSITE), int)
        self.assertNotEqual(COMPOSITE, other)
        with self.assertRaises(TypeError):
            _ = COMPOSITE > other

    def test_simple(self):
        """
        Test that simple level of theory can be obtained.
        """
        lot = COMPOSITE.simple()
        self.assertIsNot(lot, COMPOSITE)
        self.assertEqual(lot.freq.method, COMPOSITE.freq.method)
        self.assertEqual(lot.freq.basis, COMPOSITE.freq.basis)
        self.assertEqual(lot.freq.software, COMPOSITE.freq.software)
        self.assertEqual(lot.energy.method, COMPOSITE.energy.method)
        self.assertEqual(lot.energy.basis, COMPOSITE.energy.basis)
        for attr, val in lot.freq.__dict__.items():
            if attr not in {'method', 'basis', 'software'}:
                self.assertIsNone(val)
        for attr, val in lot.energy.__dict__.items():
            if attr not in {'method', 'basis'}:
                self.assertIsNone(val)

    def test_to_model_chem(self):
        """
        Test conversion to model chemistry.
        """
        self.assertEqual(COMPOSITE.to_model_chem(), COMPOSITE_MODELCHEM)

    def test_as_dict(self):
        """
        Test conversion to dictionary.
        """
        self.assertDictEqual(COMPOSITE.as_dict(), COMPOSITE_DICT)


class TestFuncs(unittest.TestCase):
    """
    A class for testing that the functions in the modelchem module work.
    """

    def test_model_chem_to_lot(self):
        """
        Test model chemistry to quantum calculation settings conversion.
        """
        self.assertEqual(
            model_chem_to_lot(FREQ_MODELCHEM, software='gaussian', args='verytight'),
            FREQ
        )
        self.assertEqual(
            model_chem_to_lot(FREQ_MODELCHEM,
                              freq_settings={'software': 'gaussian', 'args': 'verytight'}),
            FREQ
        )
        self.assertEqual(
            model_chem_to_lot(FREQ_MODELCHEM,
                              freq_settings={'software': 'gaussian', 'args': 'verytight'},
                              energy_settings={'unused setting': None}),
            FREQ
        )
        self.assertEqual(
            model_chem_to_lot(ENERGY_MODELCHEM, energy_settings={'software': 'orca'}),
            ENERGY
        )
        self.assertEqual(
            model_chem_to_lot(COMPOSITE_MODELCHEM,
                              freq_settings={'software': 'gaussian', 'args': 'verytight'},
                              energy_settings={'software': 'orca'}),
            COMPOSITE
        )

    def test_str_to_lot(self):
        """
        Test key to quantum calculation settings conversion.
        """
        self.assertEqual(str_to_lot(FREQ_REPR), FREQ)
        self.assertEqual(str_to_lot(ENERGY_REPR), ENERGY)
        self.assertEqual(str_to_lot(COMPOSITE_REPR), COMPOSITE)

    def test_get_software_id(self):
        """
        Test standardized software identifiers.
        """
        test_names = ['gaussian', 'Gaussian 09', 'g-16', 'Gau  03']
        for name in test_names:
            self.assertEqual(get_software_id(name), 'gaussian')

        test_names = ['qchem', 'QChem', 'Q-Chem']
        for name in test_names:
            self.assertEqual(get_software_id(name), 'qchem')

        test_names = ['molpro', 'Molpro', 'MOLPRO']
        for name in test_names:
            self.assertEqual(get_software_id(name), 'molpro')

        test_names = ['orca', 'Orca', 'ORCA']
        for name in test_names:
            self.assertEqual(get_software_id(name), 'orca')

        test_names = ['terachem', 'Terachem', 'TeraChem', 'Tera-Chem', 'Tera Chem']
        for name in test_names:
            self.assertEqual(get_software_id(name), 'terachem')

        with self.assertRaises(ValueError):
            get_software_id('g')


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
