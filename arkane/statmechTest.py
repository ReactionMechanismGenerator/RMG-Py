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
This module contains unit tests of the :mod:`arkane.statmech` module.
"""

import os
import unittest

import numpy as np

from rmgpy.species import Species
from rmgpy.exceptions import InputError

from arkane import Arkane
from arkane.ess.qchem import QChemLog
from arkane.modelchem import LevelOfTheory
from arkane.statmech import ScanLog, StatMechJob, determine_rotor_symmetry, is_linear

################################################################################


class TestStatmech(unittest.TestCase):
    """
    Contains unit tests of the StatmechJob class.
    """

    @classmethod
    def setUp(cls):
        """A method that is run before each unit test in this class"""
        arkane = Arkane()
        cls.job_list = arkane.load_input_file(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                           'data', 'Benzyl', 'input.py'))

    def test_gaussian_log_file_error(self):
        """Test that the proper error is raised if gaussian geometry and frequency file paths are not the same"""
        job = self.job_list[-2]
        self.assertTrue(isinstance(job, StatMechJob))
        with self.assertRaises(InputError):
            job.load()

    def test_rotor_symmetry_determination(self):
        """
        Test that the correct symmetry number is determined for rotor potential scans.
        """
        path1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'NCC_NRotor.out')
        path2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'NCC_CRotor.out')
        scan_log1 = QChemLog(path1)
        scan_log2 = QChemLog(path2)
        v_list1, angle = scan_log1.load_scan_energies()
        v_list2, angle = scan_log2.load_scan_energies()
        symmetry1 = determine_rotor_symmetry(energies=v_list1, label='NCC', pivots=[])
        symmetry2 = determine_rotor_symmetry(energies=v_list2, label='NCC', pivots=[])
        self.assertEqual(symmetry1, 1)
        self.assertEqual(symmetry2, 3)

    def test_is_linear(self):
        """Test that we can determine the linearity of a molecule from it's coordinates"""
        xyz1 = np.array([
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.159076],
            [0.000000, 0.000000, -1.159076]])  # a trivial case
        xyz2 = np.array([
            [-0.06618943, -0.12360663, -0.07631983],
            [-0.79539707, 0.86755487, 1.02675668],
            [-0.68919931, 0.25421823, -1.34830853],
            [0.01546439, -1.54297548, 0.44580391],
            [1.94428095, 0.40772394, 1.03719428],
            [2.20318015, -0.14715186, -0.64755729],
            [1.59252246, 1.51178950, -0.33908352],
            [-0.87856890, -2.02453514, 0.38494433],
            [-1.34135876, 1.49608206, 0.53295071]])  # a non-linear multi-atom molecule
        xyz3 = np.array([
            [0.0000000000, 0.0000000000, 0.3146069129],
            [-1.0906813653, 0.0000000000, -0.1376405244],
            [1.0906813653, 0.0000000000, -0.1376405244]])  # NO2, a non-linear 3-atom molecule
        xyz4 = np.array([
            [0.0000000000, 0.0000000000, 0.1413439534],
            [-0.8031792912, 0.0000000000, -0.4947038368],
            [0.8031792912, 0.0000000000, -0.4947038368]])  # NH2, a non-linear 3-atom molecule
        xyz5 = np.array([
            [-0.5417345330, 0.8208150346, 0.0000000000],
            [0.9206183692, 1.6432038228, 0.0000000000],
            [-1.2739176462, 1.9692549926, 0.0000000000]])  # HSO, a non-linear 3-atom molecule
        xyz6 = np.array([
            [1.18784533, 0.98526702, 0.00000000],
            [0.04124533, 0.98526702, 0.00000000],
            [-1.02875467, 0.98526702, 0.00000000]])  # HCN, a linear 3-atom molecule
        xyz7 = np.array([
            [-4.02394116, 0.56169428, 0.00000000],
            [-5.09394116, 0.56169428, 0.00000000],
            [-2.82274116, 0.56169428, 0.00000000],
            [-1.75274116, 0.56169428, 0.00000000]])  # C2H2, a linear 4-atom molecule
        xyz8 = np.array([
            [-1.02600933, 2.12845307, 0.00000000],
            [-0.77966935, 0.95278385, 0.00000000],
            [-1.23666197, 3.17751246, 0.00000000],
            [-0.56023545, -0.09447399, 0.00000000]])  # C2H2, just 0.5 degree off from linearity, so NOT linear
        xyz9 = np.array([
            [-1.1998, 0.1610, 0.0275],
            [-1.4021, 0.6223, -0.8489],
            [-1.48302, 0.80682, -1.19946]])  # just 3 points in space on a straight line (not a physical molecule)
        xyz10 = np.array([
            [-1.1998, 0.1610, 0.0275]])  # mono-atomic species, non-linear
        xyz11 = np.array([
            [1.06026500, -0.07706800, 0.03372800],
            [3.37340700, -0.07706800, 0.03372800],
            [2.21683600, -0.07706800, 0.03372800]])  # CO2 at wb97xd/6-311+g(d,p), linear
        xyz12 = np.array([
            [1.05503600, -0.00335000, 0.09823600],
            [2.42816800, -0.00335000, 0.09823600],
            [-0.14726400, -0.00335000, 0.09823600],
            [3.63046800, -0.00335000, 0.09823600],
            [-1.21103500, -0.00335000, 0.09823600],
            [4.69423900, -0.00335000, 0.09823600]])  # C#CC#C at wb97xd/6-311+g(d,p), linear

        self.assertTrue(is_linear(xyz1))
        self.assertTrue(is_linear(xyz6))
        self.assertTrue(is_linear(xyz7))
        self.assertTrue(is_linear(xyz9))
        self.assertTrue(is_linear(xyz11))
        self.assertTrue(is_linear(xyz12))
        self.assertFalse(is_linear(xyz2))
        self.assertFalse(is_linear(xyz3))
        self.assertFalse(is_linear(xyz4))
        self.assertFalse(is_linear(xyz5))
        self.assertFalse(is_linear(xyz8))
        self.assertFalse(is_linear(xyz10))

    def test_specifying_absolute_file_paths(self):
        """Test specifying absolute file paths of statmech files"""
        h2o2_input = """#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {{'H-O': 2, 'O-O': 1}}

externalSymmetry = 2

spinMultiplicity = 1

opticalIsomers = 1

energy = {{'b3lyp/6-311+g(3df,2p)': Log('{energy}')}}

geometry = Log('{freq}')

frequencies = Log('{freq}')

rotors = [HinderedRotor(scanLog=Log('{scan}'), pivots=[1, 2], top=[1, 3], symmetry=1, fit='fourier')]

"""
        abs_arkane_path = os.path.abspath(os.path.dirname(__file__))  # this is the absolute path to `.../RMG-Py/arkane`
        energy_path = os.path.join('arkane', 'data', 'H2O2', 'sp_a19032.out')
        freq_path = os.path.join('arkane', 'data', 'H2O2', 'freq_a19031.out')
        scan_path = os.path.join('arkane', 'data', 'H2O2', 'scan_a19034.out')
        h2o2_input = h2o2_input.format(energy=energy_path, freq=freq_path, scan=scan_path)
        h2o2_path = os.path.join(abs_arkane_path, 'data', 'H2O2', 'H2O2.py')
        if not os.path.exists(os.path.dirname(h2o2_path)):
            os.makedirs(os.path.dirname(h2o2_path))
        with open(h2o2_path, 'w') as f:
            f.write(h2o2_input)
        h2o2 = Species(label='H2O2', smiles='OO')
        self.assertIsNone(h2o2.conformer)
        statmech_job = StatMechJob(species=h2o2, path=h2o2_path)
        statmech_job.level_of_theory = LevelOfTheory('b3lyp', '6-311+g(3df,2p)')
        statmech_job.load(pdep=False, plot=False)
        self.assertAlmostEqual(h2o2.conformer.E0.value_si, -146031.49933673252)
        os.remove(h2o2_path)

    def test_hinder_rotor_from_1d_array(self):
        """Test assigning hindered rotor 1D PES profile directly to HinderedRotor1DArray"""
        h2o2_input = """#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {{'H-O': 2, 'O-O': 1}}

externalSymmetry = 2

spinMultiplicity = 1

opticalIsomers = 1

energy = {{'b3lyp/6-311+g(3df,2p)': Log('{energy}')}}

geometry = Log('{freq}')

frequencies = Log('{freq}')

rotors = [HinderedRotor1DArray(
    angles=[0.        , 0.17453293, 0.34906585, 0.52359878, 0.6981317 ,
            0.87266463, 1.04719755, 1.22173048, 1.3962634 , 1.57079633,
            1.74532925, 1.91986218, 2.0943951 , 2.26892803, 2.44346095,
            2.61799388, 2.7925268 , 2.96705973, 3.14159265, 3.31612558,
            3.4906585 , 3.66519143, 3.83972435, 4.01425728, 4.1887902 ,
            4.36332313, 4.53785606, 4.71238898, 4.88692191, 5.06145483,
            5.23598776, 5.41052068, 5.58505361, 5.75958653, 5.93411946,
            6.10865238, 6.28318531],
    energies=[0.00000000e+00, 3.09449290e+02, 1.07459871e+03, 2.05925305e+03,
            3.02877926e+03, 3.79724994e+03, 4.23486826e+03, 4.26190303e+03,
            3.88196432e+03, 3.15173930e+03, 2.20016363e+03, 1.20431941e+03,
            3.94499732e+02, 7.23850312e+00, 2.77854025e+02, 1.40711827e+03,
            3.50375319e+03, 6.57899330e+03, 1.05208190e+04, 1.50847596e+04,
            1.99269611e+04, 2.46164740e+04, 2.86972097e+04, 3.17430074e+04,
            3.34148312e+04, 3.35267510e+04, 3.20643922e+04, 2.91936786e+04,
            2.52325029e+04, 2.06007483e+04, 1.57531541e+04, 1.11268684e+04,
            7.08120679e+03, 3.87554760e+03, 1.63995547e+03, 3.80256396e+02,
            6.14367036e-01],
    pivots=[1, 2], top=[1, 3], symmetry=1, fit='fourier')]
"""
        angles = np.array([0.        , 0.17453293, 0.34906585, 0.52359878, 0.6981317 ,
            0.87266463, 1.04719755, 1.22173048, 1.3962634 , 1.57079633,
            1.74532925, 1.91986218, 2.0943951 , 2.26892803, 2.44346095,
            2.61799388, 2.7925268 , 2.96705973, 3.14159265, 3.31612558,
            3.4906585 , 3.66519143, 3.83972435, 4.01425728, 4.1887902 ,
            4.36332313, 4.53785606, 4.71238898, 4.88692191, 5.06145483,
            5.23598776, 5.41052068, 5.58505361, 5.75958653, 5.93411946,
            6.10865238, 6.28318531])
        energies = np.array([0.00000000e+00, 3.09449290e+02, 1.07459871e+03, 2.05925305e+03,
            3.02877926e+03, 3.79724994e+03, 4.23486826e+03, 4.26190303e+03,
            3.88196432e+03, 3.15173930e+03, 2.20016363e+03, 1.20431941e+03,
            3.94499732e+02, 7.23850312e+00, 2.77854025e+02, 1.40711827e+03,
            3.50375319e+03, 6.57899330e+03, 1.05208190e+04, 1.50847596e+04,
            1.99269611e+04, 2.46164740e+04, 2.86972097e+04, 3.17430074e+04,
            3.34148312e+04, 3.35267510e+04, 3.20643922e+04, 2.91936786e+04,
            2.52325029e+04, 2.06007483e+04, 1.57531541e+04, 1.11268684e+04,
            7.08120679e+03, 3.87554760e+03, 1.63995547e+03, 3.80256396e+02,
            6.14367036e-01])
        abs_arkane_path = os.path.abspath(os.path.dirname(__file__))  # this is the absolute path to `.../RMG-Py/arkane`
        energy_path = os.path.join('arkane', 'data', 'H2O2', 'sp_a19032.out')
        freq_path = os.path.join('arkane', 'data', 'H2O2', 'freq_a19031.out')
        h2o2_input = h2o2_input.format(energy=energy_path, freq=freq_path, angles=angles, energies=energies)
        h2o2_path = os.path.join(abs_arkane_path, 'data', 'H2O2', 'H2O2_PES.py')
        os.makedirs(os.path.dirname(h2o2_path), exist_ok=True)
        with open(h2o2_path, 'w') as f:
            f.write(h2o2_input)
        h2o2 = Species(label='H2O2', smiles='OO')
        self.assertIsNone(h2o2.conformer)
        statmech_job = StatMechJob(species=h2o2, path=h2o2_path)
        statmech_job.level_of_theory = LevelOfTheory('b3lyp', '6-311+g(3df,2p)')
        statmech_job.load(pdep=False, plot=False)
        self.assertEqual(len(statmech_job.raw_hindered_rotor_data), 1)
        self.assertEqual(statmech_job.raw_hindered_rotor_data[0][2], 1)
        self.assertTrue(np.allclose(statmech_job.raw_hindered_rotor_data[0][3], angles, atol=1e-6))
        self.assertTrue(np.allclose(statmech_job.raw_hindered_rotor_data[0][4], energies, atol=1e-6))
        self.assertAlmostEqual(h2o2.conformer.E0.value_si, -146031.49933673252)
        os.remove(h2o2_path)

    def test_scanlog_class(self):
        """
        Test scanlog works for various input format and returns the correct PES profiles.
        """
        angles = np.array([0.        , 0.17453293, 0.34906585, 0.52359878, 0.6981317 ,
            0.87266463, 1.04719755, 1.22173048, 1.3962634 , 1.57079633,
            1.74532925, 1.91986218, 2.0943951 , 2.26892803, 2.44346095,
            2.61799388, 2.7925268 , 2.96705973, 3.14159265, 3.31612558,
            3.4906585 , 3.66519143, 3.83972435, 4.01425728, 4.1887902 ,
            4.36332313, 4.53785606, 4.71238898, 4.88692191, 5.06145483,
            5.23598776, 5.41052068, 5.58505361, 5.75958653, 5.93411946,
            6.10865238, 6.28318531])
        energies = np.array([0.00000000e+00, 3.09449290e+02, 1.07459871e+03, 2.05925305e+03,
            3.02877926e+03, 3.79724994e+03, 4.23486826e+03, 4.26190303e+03,
            3.88196432e+03, 3.15173930e+03, 2.20016363e+03, 1.20431941e+03,
            3.94499732e+02, 7.23850312e+00, 2.77854025e+02, 1.40711827e+03,
            3.50375319e+03, 6.57899330e+03, 1.05208190e+04, 1.50847596e+04,
            1.99269611e+04, 2.46164740e+04, 2.86972097e+04, 3.17430074e+04,
            3.34148312e+04, 3.35267510e+04, 3.20643922e+04, 2.91936786e+04,
            2.52325029e+04, 2.06007483e+04, 1.57531541e+04, 1.11268684e+04,
            7.08120679e+03, 3.87554760e+03, 1.63995547e+03, 3.80256396e+02,
            6.14367036e-01])
        abs_arkane_path = os.path.abspath(os.path.dirname(__file__))
        scanpath1 = os.path.join(abs_arkane_path, 'data', 'H2O2', 'scan.txt')
        scanlog1 = ScanLog(scanpath1)
        angles1, energies1 = scanlog1.load()
        self.assertTrue(np.allclose(angles, angles1, atol=1e-6))
        self.assertTrue(np.allclose(energies, energies1, atol=1e-6))

        scanpath2 = os.path.join(abs_arkane_path, 'data', 'H2O2', 'scan.yml')
        scanlog2 = ScanLog(scanpath2)
        angles2, energies2 = scanlog2.load()
        self.assertTrue(np.allclose(angles, angles2, atol=1e-6))
        print(energies, energies2)
        self.assertTrue(np.allclose(energies, energies2, atol=1e-6))

        scanpath3 = os.path.join(abs_arkane_path, 'data', 'H2O2', 'scan.csv')
        scanlog3 = ScanLog(scanpath3)
        angles3, energies3 = scanlog3.load()
        self.assertTrue(np.allclose(angles, angles3, atol=1e-6))
        self.assertTrue(np.allclose(energies, energies3, atol=1e-6))

    def test_hindered_rotor_from_scan_logs(self):
        """
        Test assigning hindered rotor 1D PES profile via ScanLog to HinderedRotor in statmech jobs.
        """
        angles = np.array([0.        , 0.17453293, 0.34906585, 0.52359878, 0.6981317 ,
            0.87266463, 1.04719755, 1.22173048, 1.3962634 , 1.57079633,
            1.74532925, 1.91986218, 2.0943951 , 2.26892803, 2.44346095,
            2.61799388, 2.7925268 , 2.96705973, 3.14159265, 3.31612558,
            3.4906585 , 3.66519143, 3.83972435, 4.01425728, 4.1887902 ,
            4.36332313, 4.53785606, 4.71238898, 4.88692191, 5.06145483,
            5.23598776, 5.41052068, 5.58505361, 5.75958653, 5.93411946,
            6.10865238, 6.28318531])
        energies = np.array([0.00000000e+00, 3.09449290e+02, 1.07459871e+03, 2.05925305e+03,
            3.02877926e+03, 3.79724994e+03, 4.23486826e+03, 4.26190303e+03,
            3.88196432e+03, 3.15173930e+03, 2.20016363e+03, 1.20431941e+03,
            3.94499732e+02, 7.23850312e+00, 2.77854025e+02, 1.40711827e+03,
            3.50375319e+03, 6.57899330e+03, 1.05208190e+04, 1.50847596e+04,
            1.99269611e+04, 2.46164740e+04, 2.86972097e+04, 3.17430074e+04,
            3.34148312e+04, 3.35267510e+04, 3.20643922e+04, 2.91936786e+04,
            2.52325029e+04, 2.06007483e+04, 1.57531541e+04, 1.11268684e+04,
            7.08120679e+03, 3.87554760e+03, 1.63995547e+03, 3.80256396e+02,
            6.14367036e-01])
        h2o2_input = """#!/usr/bin/env python
# -*- coding: utf-8 -*-

bonds = {{'H-O': 2, 'O-O': 1}}

externalSymmetry = 2

spinMultiplicity = 1

opticalIsomers = 1

energy = {{'b3lyp/6-311+g(3df,2p)': Log('{energy}')}}

geometry = Log('{freq}')

frequencies = Log('{freq}')

rotors = [HinderedRotor(scanLog=ScanLog('{scan}'), pivots=[1, 2], top=[1, 3], symmetry=1, fit='fourier')]

"""
        abs_arkane_path = os.path.abspath(os.path.dirname(__file__))  # this is the absolute path to `.../RMG-Py/arkane`
        energy_path = os.path.join(abs_arkane_path, 'data', 'H2O2', 'sp_a19032.out')
        freq_path = os.path.join(abs_arkane_path, 'data', 'H2O2', 'freq_a19031.out')
        h2o2_path = os.path.join(abs_arkane_path, 'data', 'H2O2', 'H2O2.py')
        h2o2 = Species(label='H2O2', smiles='OO')
        os.makedirs(os.path.dirname(h2o2_path), exist_ok=True)

        for file in ['scan.txt', 'scan.csv', 'scan.yml']:
            scan_path = os.path.join(abs_arkane_path, 'data', 'H2O2', file)
            h2o2_input_tmp = h2o2_input.format(energy=energy_path, freq=freq_path, scan=scan_path)
            with open(h2o2_path, 'w') as f:
                f.write(h2o2_input_tmp)
            statmech_job = StatMechJob(species=h2o2, path=h2o2_path)
            statmech_job.level_of_theory = LevelOfTheory('b3lyp', '6-311+g(3df,2p)')
            statmech_job.load(pdep=False, plot=False)
            self.assertEqual(len(statmech_job.raw_hindered_rotor_data), 1)
            self.assertTrue(np.allclose(statmech_job.raw_hindered_rotor_data[0][3], angles, atol=1e-6))
            self.assertTrue(np.allclose(statmech_job.raw_hindered_rotor_data[0][4], energies, atol=1e-6))
            os.remove(h2o2_path)



################################################################################


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
