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
This module contains unit tests of the :mod:`arkane.input` module.
"""

import os
import unittest

import rmgpy
from rmgpy.exceptions import InputError
from rmgpy.kinetics.tunneling import Eckart
from rmgpy.pdep.collision import SingleExponentialDown
from rmgpy.statmech.rotation import NonlinearRotor
from rmgpy.statmech.translation import IdealGasTranslation
from rmgpy.statmech.vibration import HarmonicOscillator
from rmgpy.thermo.nasa import NASAPolynomial, NASA
from rmgpy.transport import TransportData

from arkane.input import species, transitionState, reaction, SMILES, load_input_file, process_model_chemistry
from arkane.modelchem import LevelOfTheory, CompositeLevelOfTheory

################################################################################


class InputTest(unittest.TestCase):
    """
    Contains unit tests for the Arkane input module
    """

    def test_species(self):
        """
        Test loading a species from input file-like kew word arguments
        """
        label0 = 'CH2O'
        kwargs = {'E0': (28.69, 'kcal/mol'),
                  'structure': SMILES('C=O'),
                  'collisionModel': TransportData(sigma=(3.69e-10, 'm'), epsilon=(4.0, 'kJ/mol')),
                  'energyTransferModel': SingleExponentialDown(alpha0=(0.956, 'kJ/mol'), T0=(300, 'K'), n=0.95),
                  'spinMultiplicity': 1,
                  'opticalIsomers': 1,
                  'modes': [HarmonicOscillator(frequencies=([1180, 1261, 1529, 1764, 2931, 2999], 'cm^-1')),
                            NonlinearRotor(rotationalConstant=([1.15498821005263, 1.3156969584727, 9.45570474524524],
                                                               "cm^-1"), symmetry=2, quantum=False),
                            IdealGasTranslation(mass=(30.0106, "g/mol")),
                            ]}

        spc0 = species(label0, **kwargs)
        self.assertEqual(spc0.label, 'CH2O')
        self.assertEqual(spc0.smiles, 'C=O')
        self.assertAlmostEqual(spc0.conformer.E0.value_si, 120038.96)
        self.assertEqual(spc0.conformer.spin_multiplicity, 1)
        self.assertEqual(spc0.conformer.optical_isomers, 1)
        self.assertEqual(len(spc0.conformer.modes), 3)
        self.assertIsInstance(spc0.transport_data, TransportData)
        self.assertIsInstance(spc0.energy_transfer_model, SingleExponentialDown)

    def test_species_atomic_nasa_polynomial(self):
        """
        Test loading a atom with NASA polynomials
        """
        label0 = "H(1)"
        kwargs = {"structure": SMILES('[H]'),
                  "thermo": NASA(polynomials=[
                      NASAPolynomial(coeffs=[2.5, 0, 0, 0, 0, 25473.7, -0.446683], Tmin=(200, 'K'), Tmax=(1000, 'K')),
                      NASAPolynomial(coeffs=[2.5, 0, 0, 0, 0, 25473.7, -0.446683], Tmin=(1000, 'K'), Tmax=(6000, 'K'))],
                      Tmin=(200, 'K'), Tmax=(6000, 'K'), comment="""Thermo library: FFCM1(-)"""),
                  "energyTransferModel": SingleExponentialDown(alpha0=(3.5886, 'kJ/mol'), T0=(300, 'K'), n=0.85)}
        spc0 = species(label0, **kwargs)
        self.assertEqual(spc0.label, label0)
        self.assertEqual(spc0.smiles, '[H]')
        self.assertTrue(spc0.has_statmech())
        self.assertEqual(spc0.thermo, kwargs['thermo'])

    def test_species_polyatomic_nasa_polynomial(self):
        """
        Test loading a species with NASA polynomials
        """
        label0 = "benzyl"
        kwargs = {"structure": SMILES('[c]1ccccc1'),
                  "thermo": NASA(polynomials=[NASAPolynomial(
                      coeffs=[2.78632, 0.00784632, 7.97887e-05, -1.11617e-07, 4.39429e-11, 39695, 11.5114],
                      Tmin=(100, 'K'), Tmax=(943.73, 'K')),
                      NASAPolynomial(
                          coeffs=[13.2455, 0.0115667, -2.49996e-06, 4.66496e-10, -4.12376e-14,
                                  35581.1, -49.6793], Tmin=(943.73, 'K'), Tmax=(5000, 'K'))],
                      Tmin=(100, 'K'), Tmax=(5000, 'K'),
                      comment="""Thermo library: Fulvene_H + radical(CbJ)"""),
                  "energyTransferModel": SingleExponentialDown(alpha0=(3.5886, 'kJ/mol'), T0=(300, 'K'), n=0.85)}
        spc0 = species(label0, **kwargs)
        self.assertEqual(spc0.label, label0)
        self.assertTrue(spc0.has_statmech())
        self.assertEqual(spc0.thermo, kwargs['thermo'])

    def test_transition_state(self):
        """
        Test loading a transition state from input file-like kew word arguments
        """
        label0 = 'TS1'
        kwargs = {'E0': (39.95, 'kcal/mol'),
                  'spinMultiplicity': 2,
                  'opticalIsomers': 1,
                  'frequency': (-1934, 'cm^-1'),
                  'modes': [HarmonicOscillator(frequencies=([792, 987, 1136, 1142, 1482, 2441, 3096, 3183], 'cm^-1')),
                            NonlinearRotor(rotationalConstant=([0.928, 0.962, 5.807], "cm^-1"), symmetry=1,
                                           quantum=False),
                            IdealGasTranslation(mass=(31.01843, "g/mol"))]}

        ts0 = transitionState(label0, **kwargs)
        self.assertEqual(ts0.label, 'TS1')
        self.assertAlmostEqual(ts0.conformer.E0.value_si, 167150.8)
        self.assertEqual(ts0.conformer.spin_multiplicity, 2)
        self.assertEqual(ts0.conformer.optical_isomers, 1)
        self.assertEqual(ts0.frequency.value_si, -1934.0)
        self.assertEqual(len(ts0.conformer.modes), 3)

    def test_reaction(self):
        """
        Test loading a reaction from input file-like kew word arguments
        """

        species(
            label='methoxy',
            structure=SMILES('C[O]'),
            E0=(9.44, 'kcal/mol'),
            modes=[
                HarmonicOscillator(frequencies=([758, 960, 1106, 1393, 1403, 1518, 2940, 3019, 3065], 'cm^-1')),
                NonlinearRotor(rotationalConstant=([0.916, 0.921, 5.251], "cm^-1"), symmetry=3, quantum=False),
                IdealGasTranslation(mass=(31.01843, "g/mol"))],
            spinMultiplicity=2,
            opticalIsomers=1,
            molecularWeight=(31.01843, 'amu'),
            collisionModel=TransportData(sigma=(3.69e-10, 'm'), epsilon=(4.0, 'kJ/mol')),
            energyTransferModel=SingleExponentialDown(alpha0=(0.956, 'kJ/mol'), T0=(300, 'K'), n=0.95))

        species(
            label='formaldehyde',
            E0=(28.69, 'kcal/mol'),
            molecularWeight=(30.0106, "g/mol"),
            collisionModel=TransportData(sigma=(3.69e-10, 'm'), epsilon=(4.0, 'kJ/mol')),
            energyTransferModel=SingleExponentialDown(alpha0=(0.956, 'kJ/mol'), T0=(300, 'K'), n=0.95),
            spinMultiplicity=1,
            opticalIsomers=1,
            modes=[HarmonicOscillator(frequencies=([1180, 1261, 1529, 1764, 2931, 2999], 'cm^-1')),
                   NonlinearRotor(rotationalConstant=([1.15498821005263, 1.3156969584727, 9.45570474524524], "cm^-1"),
                                  symmetry=2, quantum=False),
                   IdealGasTranslation(mass=(30.0106, "g/mol"))])

        species(
            label='H',
            E0=(0.000, 'kcal/mol'),
            molecularWeight=(1.00783, "g/mol"),
            collisionModel=TransportData(sigma=(3.69e-10, 'm'), epsilon=(4.0, 'kJ/mol')),
            energyTransferModel=SingleExponentialDown(alpha0=(0.956, 'kJ/mol'), T0=(300, 'K'), n=0.95),
            modes=[IdealGasTranslation(mass=(1.00783, "g/mol"))],
            spinMultiplicity=2,
            opticalIsomers=1)

        transitionState(
            label='TS3',
            E0=(34.1, 'kcal/mol'),
            spinMultiplicity=2,
            opticalIsomers=1,
            frequency=(-967, 'cm^-1'),
            modes=[HarmonicOscillator(frequencies=([466, 581, 1169, 1242, 1499, 1659, 2933, 3000], 'cm^-1')),
                   NonlinearRotor(rotationalConstant=([0.970, 1.029, 3.717], "cm^-1"), symmetry=1, quantum=False),
                   IdealGasTranslation(mass=(31.01843, "g/mol"))])

        reactants = ['formaldehyde', 'H']
        products = ['methoxy']
        tunneling = 'Eckart'

        rxn = reaction('CH2O+H=Methoxy', reactants, products, 'TS3', tunneling=tunneling)
        self.assertEqual(rxn.label, 'CH2O+H=Methoxy')
        self.assertEqual(len(rxn.reactants), 2)
        self.assertEqual(len(rxn.products), 1)
        self.assertAlmostEqual(rxn.reactants[0].conformer.E0.value_si, 0)
        self.assertAlmostEqual(rxn.reactants[1].conformer.E0.value_si, 120038.96)
        self.assertAlmostEqual(rxn.products[0].conformer.E0.value_si, 39496.96)
        self.assertAlmostEqual(rxn.transition_state.conformer.E0.value_si, 142674.4)
        self.assertAlmostEqual(rxn.transition_state.frequency.value_si, -967.0)
        self.assertIsInstance(rxn.transition_state.tunneling, Eckart)

    def test_load_input_file(self):
        """Test loading an Arkane input file"""
        path = os.path.join(os.path.dirname(os.path.dirname(rmgpy.__file__)), 'examples', 'arkane', 'networks',
                            'acetyl+O2', 'input.py')
        job_list, reaction_dict, species_dict, transition_state_dict, network_dict, model_chemistry \
            = load_input_file(path)

        self.assertEqual(len(job_list), 1)

        self.assertEqual(len(reaction_dict), 5)
        self.assertTrue('entrance1' in reaction_dict)
        self.assertTrue('exit2' in reaction_dict)

        self.assertEqual(len(species_dict), 9)
        self.assertTrue('acetyl' in species_dict)
        self.assertTrue('hydroperoxyl' in species_dict)

        self.assertEqual(len(transition_state_dict), 5)
        self.assertTrue('entrance1' in transition_state_dict)
        self.assertTrue('isom1' in transition_state_dict)

        self.assertEqual(len(network_dict), 1)
        self.assertTrue('acetyl + O2' in network_dict)

        self.assertIsNone(model_chemistry)

    def test_process_model_chemistry(self):
        """
        Test processing the model chemistry to derive the sp and freq levels
        """
        mc = 'ccsd(t)-f12a/aug-cc-pvtz//b3lyp/6-311++g(3df,3pd)'
        lot = process_model_chemistry(mc)
        self.assertIsInstance(lot, CompositeLevelOfTheory)
        self.assertEqual(lot.energy, LevelOfTheory('ccsd(t)-f12a', 'aug-cc-pvtz'))
        self.assertEqual(lot.freq, LevelOfTheory('b3lyp', '6-311++g(3df,3pd)'))

        mc = 'b3lyp-d3/def2-tzvp'
        lot = process_model_chemistry(mc)
        self.assertIsInstance(lot, LevelOfTheory)
        self.assertEqual(lot, LevelOfTheory('b3lyp-d3', 'def2-tzvp'))

        mc = 'cbs-qb3'
        lot = process_model_chemistry(mc)
        self.assertIsInstance(lot, LevelOfTheory)
        self.assertEqual(lot, LevelOfTheory('cbs-qb3'))

        mc = LevelOfTheory('test')
        lot = process_model_chemistry(mc)
        self.assertIs(mc, lot)

        with self.assertRaises(InputError):
            process_model_chemistry('CCSD(T)-F12a/aug-cc-pVTZ//CCSD(T)-F12a/aug-cc-pVTZ//B3LYP/6-311++G(3df,3pd)')

################################################################################


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
