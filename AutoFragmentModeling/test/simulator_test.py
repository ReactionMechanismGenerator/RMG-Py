import os
import shutil
import unittest

import afm.utils
import afm.simulator
import afm.molecule

class TestSimulator(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        chemkin_path = os.path.join(os.path.dirname(__file__), 
                                    'data', 
                                    'mc_simulator_data',
                                    'chem.inp')

        dictionary_path = os.path.join(os.path.dirname(__file__), 
                                    'data', 
                                    'mc_simulator_data',
                                    'species_dictionary.txt')

        fragment_smiles_path = os.path.join(os.path.dirname(__file__), 
                                    'data', 
                                    'mc_simulator_data',
                                    'fragment_smiles.txt')

        self.simulator = afm.simulator.Simulator(chemkin_path, 
                                                dictionary_path,
                                                fragment_smiles_path)

    def test_fragment_chemistry(self):

        fragment_dict = self.simulator.fragment_dict
        fragment_reaction_list = self.simulator.fragment_reaction_list
        self.assertEqual(40, len(fragment_dict))

        # change from 5 to 4 due to no pseudo fragment rxn
        self.assertEqual(4, len(fragment_reaction_list))

class TestOdeSimulator(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        chemkin_path = os.path.join(os.path.dirname(__file__), 
                                    'data', 
                                    'ode_simulator_data',
                                    'chem.inp')

        dictionary_path = os.path.join(os.path.dirname(__file__), 
                                    'data', 
                                    'ode_simulator_data',
                                    'species_dictionary.txt')

        fragment_smiles_path = os.path.join(os.path.dirname(__file__), 
                                    'data', 
                                    'ode_simulator_data',
                                    'fragment_smiles.txt')

        temperature = 673.15 # unit: K
        pressure = 350*3.75 # unit: bar
        self.outputDirectory = 'temp'
        self.odes = afm.simulator.OdeSimulator(chemkin_path,
                                              dictionary_path,
                                              fragment_smiles_path,
                                              temperature,
                                              pressure,
                                              self.outputDirectory)
    @classmethod
    def tearDownClass(self):

        shutil.rmtree(self.outputDirectory)


    def test_simulate(self):

        initial_mol_fraction = {
                                "ArCCCCR":1.0,
                                "LCCCCR":1.75,
                                "LC":1.0
                                }

        termination_time = 3600*14 # unit: sec
        all_data = self.odes.simulate(initial_mol_fraction, termination_time)

        # only one condition
        self.assertEqual(len(all_data), 1)

        # there's three parts of data
        self.assertEqual(len(all_data[0]), 3)

        time, dataList, reactionSensitivityData = all_data[0]
        TData = dataList[0]
        PData = dataList[1]
        VData = dataList[2]
        total_moles = PData.data*VData.data/8.314/TData.data
        
        arccccr_mf = dataList[3].data
        arccccr_moles = arccccr_mf*total_moles
        arccccr_conv = (arccccr_moles[0]-arccccr_moles[-1])/arccccr_moles[0]
        self.assertAlmostEqual(arccccr_conv, 0.82, 2)

    def test_reattach_fragments(self):

        initial_mol_fraction = {
                                "ArCCCCR":1.0,
                                "LCCCCR":1.75,
                                "LC":1.0
                                }

        termination_time = 3600*14 # unit: sec
        all_data = self.odes.simulate(initial_mol_fraction, termination_time)

        _, dataList, _ = all_data[0]
        TData = dataList[0]
        PData = dataList[1]
        VData = dataList[2]
        total_moles = PData.data*VData.data/8.314/TData.data
        moles_dict = {}
        for data in dataList[3:]:
            spe_label = data.label
            if '*' in spe_label:
                continue
            r_count = spe_label.count('R')
            l_count = spe_label.count('L')
            label_count = r_count + l_count

            if label_count == 0:
                continue
            if abs(data.data[-1]*total_moles[-1]) <= 1e-6:
                continue
            moles_dict[spe_label] = max(data.data[-1]*total_moles[-1],0)

        r_moles, l_moles, r_l_moles, _, rr_ll_list = afm.simulator.categorize_fragments(moles_dict)
        matches = self.odes.reattach_fragments(r_moles, 
                                               l_moles, 
                                               r_l_moles,
                                               rr_ll_list)

        moles_dict_after_match = {}
        for match in matches:
            combo, val = match
            frags = afm.utils.flatten(combo)
            for frag in frags:
                if frag not in moles_dict_after_match:
                    moles_dict_after_match[frag] = val
                else:
                    moles_dict_after_match[frag] += val

        for frag, val in moles_dict.iteritems():
            val_after_match = moles_dict_after_match[frag]
            diff_pct = abs(val - val_after_match)/moles_dict[frag]
            self.assertAlmostEqual(diff_pct, 0.0, 6)

# 1 label (R label only)
class TestOdeSimulator_1_label(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        chemkin_path = os.path.join(os.path.dirname(__file__),
                                    'data',
                                    '2mo_simulator_data',
                                    'chem.inp')

        dictionary_path = os.path.join(os.path.dirname(__file__),
                                       'data',
                                       '2mo_simulator_data',
                                       'species_dictionary.txt')

        fragment_smiles_path = os.path.join(os.path.dirname(__file__),
                                            'data',
                                            '2mo_simulator_data',
                                            'fragment_smiles.txt')

        temperature = 673.15  # unit: K
        pressure = 350 * 3.75  # unit: bar
        self.outputDirectory = 'temp'
        self.odes = afm.simulator.OdeSimulator(chemkin_path,
                                                   dictionary_path,
                                                   fragment_smiles_path,
                                                   temperature,
                                                   pressure,
                                                   self.outputDirectory)

    @classmethod
    def tearDownClass(self):

        shutil.rmtree(self.outputDirectory)

    def test_simulate(self):

        initial_mol_fraction = {
            "ArCC(C)R": 1.0,
            "RCCCCR": 1.0,
            "RCC": 1.0
        }

        termination_time = 3600 * 14  # unit: sec
        all_data = self.odes.simulate(initial_mol_fraction, termination_time)

        # only one condition
        self.assertEqual(len(all_data), 1)

        # there's three parts of data
        self.assertEqual(len(all_data[0]), 3)

        time, dataList, reactionSensitivityData = all_data[0]
        TData = dataList[0]
        PData = dataList[1]
        VData = dataList[2]
        total_moles = PData.data * VData.data / 8.314 / TData.data

        arcc_c_r_mf = dataList[3].data
        arcc_c_r_moles = arcc_c_r_mf * total_moles
        arcc_c_r_conv = (arcc_c_r_moles[0] - arcc_c_r_moles[-1]) / arcc_c_r_moles[0]
        self.assertAlmostEqual(arcc_c_r_conv, 0.18, 2)

    def test_reattach_fragments_1_label(self):

        initial_mol_fraction = {
            "ArCC(C)R": 1.0,
            "RCCCCR": 1.0,
            "RCC": 1.0
        }

        termination_time = 3600 * 14  # unit: sec
        all_data = self.odes.simulate(initial_mol_fraction, termination_time)

        _, dataList, _ = all_data[0]
        TData = dataList[0]
        PData = dataList[1]
        VData = dataList[2]
        total_moles = PData.data * VData.data / 8.314 / TData.data
        moles_dict = {}

        for data in dataList[3:]:
            spe_label = data.label
            if '*' in spe_label:
                continue
            r_count = spe_label.count('R')

            label_count = r_count

            if label_count == 0:
                continue
            if abs(data.data[-1] * total_moles[-1]) <= 1e-6:
                continue
            moles_dict[spe_label] = max(data.data[-1] * total_moles[-1], 0)

        # Try simple example of fragment sequence
        moles_dict = {}
        spe_label = ['ArCC(C)R', 'RCC', 'RCCCCR', 'RCC__C', 'ArC__CR', 'RCC__CCR']
        values = [6, 2, 3, 6, 4, 3]
        i = 0

        for spe_name in spe_label:
            moles_dict[spe_name] = values[i]
            i += 1

        r_moles, _, rr_list = afm.simulator.categorize_fragments_1_label(moles_dict)
        matches = self.odes.reattach_fragments_1_label(r_moles,
                                                       rr_list)

        # estimated matches:

        # [((('ArCC(C)R', 'RCC__CCR', 'ArCC(C)R'), ('RCCCCR', 'RCCCCR')), 0.5),
        # ((('RCC__C', 'RCC__CCR', 'RCC__C'), ('RCCCCR', 'RCCCCR')), 0.5),
        # ((('ArC__CR', 'RCC__C'), ('RCC__CCR', 'RCC__CCR')), 0.5),
        # (('ArC__CR', 'RCC__C'), 0.5),
        # (('RCC__C', 'RCC__CCR', 'RCC__C'), 0.5),
        # (('RCC__C', 'ArCC(C)R'), 1),
        # (('RCC__C', 'ArCC(C)R'), 1),
        # (('ArCC(C)R', 'RCCCCR', 'ArCC(C)R'), 0.5),
        # (('RCC', 'RCC__C'), 1),
        # (('ArCC(C)R', 'RCC__CCR', 'ArCC(C)R'), 0.5),
        # (('ArCC(C)R', 'RCCCCR', 'ArCC(C)R'), 0.5),
        # (('ArC__CR', 'ArC__CR'), 1),
        # (('RCC', 'ArC__CR'), 1)]

        moles_dict_after_match = {}
        for match in matches:
            combo, val = match
            frags = afm.utils.flatten(combo)
            for frag in frags:
                if frag not in moles_dict_after_match:
                    moles_dict_after_match[frag] = val
                else:
                    moles_dict_after_match[frag] += val

        for frag, val in moles_dict.iteritems():
            val_after_match = moles_dict_after_match[frag]
            diff_pct = abs(val - val_after_match) / moles_dict[frag]
            self.assertAlmostEqual(diff_pct, 0.0, 6)

class TestMonteCarloSimulator(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        chemkin_path = os.path.join(os.path.dirname(__file__), 
                                    'data', 
                                    'mc_simulator_data',
                                    'chem.inp')

        dictionary_path = os.path.join(os.path.dirname(__file__), 
                                    'data', 
                                    'mc_simulator_data',
                                    'species_dictionary.txt')

        fragment_smiles_path = os.path.join(os.path.dirname(__file__), 
                                    'data', 
                                    'mc_simulator_data',
                                    'fragment_smiles.txt')

        composition1 = {'ArCCCCR': 500, 'ArC__C': 1}
        composition2 = {'ArCCCCR': 300, 'ArC__C': 9}
        composition3 = {'ArCCCCR': 200}
        mol1 = afm.molecule.FragmentMolecule(composition1)
        mol2 = afm.molecule.FragmentMolecule(composition2)
        mol3 = afm.molecule.FragmentMolecule(composition3)
        initial_molecules = [mol1, mol2, mol3]

        volume = 4.8e-25 # unit: m^3
        temperature = 700 # unit: K
        self.mcs = afm.simulator.MonteCarloSimulator(chemkin_path, 
                                                     dictionary_path,
                                                     fragment_smiles_path,
                                                     initial_molecules,
                                                     volume, 
                                                     temperature)

        self.mcs.update_reaction_fluxes()

    def test_initialize_fragment_counts(self):

        # test all thre fragment labels are in the 
        # fragment_count_dict
        for fragment_label in self.mcs.fragment_dict:
            self.assertIn(fragment_label, self.mcs.fragment_count_dict)

        # check ArCCCCR column, order follows the order of molecules
        # in the initial_molecules list
        self.assertEqual(1000, self.mcs.fragment_count_dict["ArCCCCR"])
        self.assertEqual(10, self.mcs.fragment_count_dict["ArC__C"])

    def test_initialize_molecule_fragment_df(self):

        # test all thre fragment labels are in the dataframe
        # columns
        for fragment_label in self.mcs.fragment_dict:
            self.assertIn(fragment_label, self.mcs.molecule_fragment_df.columns.values)

        # check ArCCCCR column, order follows the order of molecules
        # in the initial_molecules list
        self.assertEqual(500, self.mcs.molecule_fragment_df["ArCCCCR"].values[0])
        self.assertEqual(300, self.mcs.molecule_fragment_df["ArCCCCR"].values[1])
        self.assertEqual(200, self.mcs.molecule_fragment_df["ArCCCCR"].values[2])
        # check ArC__C column
        self.assertEqual(1, self.mcs.molecule_fragment_df["ArC__C"].values[0])
        self.assertEqual(9, self.mcs.molecule_fragment_df["ArC__C"].values[1])

    def test_calculate_unimolucular_rate(self):

        # select 157th reaction: ArCCCCR => ArC* + RCCC*
        uni_reaction = self.mcs.fragment_reaction_list[2]
        rate_u = self.mcs.calculate_unimolucular_rate(uni_reaction)

        expected_rate_u = 1.51e-4
        self.assertAlmostEqual(expected_rate_u, rate_u, 6)

    def test_calculate_bimolucular_rate(self):

        # select 182th reaction: ArC__C + ArCCCCR => ArC*CCCR + ArC*C
        bi_reaction = self.mcs.fragment_reaction_list[3]
        rate_b = self.mcs.calculate_bimolucular_rate(bi_reaction)

        expected_rate_b = 1.43e-3
        self.assertAlmostEqual(expected_rate_b, rate_b, 5)

    def test_time_step(self):

        time_step = self.mcs.time_step()

        expected_time_step = 633.0

        self.assertAlmostEqual(time_step, expected_time_step, 1)

    def test_random_reaction_selection(self):

        random_reaction_idx = self.mcs.random_reaction_selection(seed=4)

        expected_idx = 3
        self.assertEqual(expected_idx, random_reaction_idx)

    def test_update_fragment_counts(self):

        reaction_idx = 3
        self.mcs.update_fragment_counts(reaction_idx)
        self.assertEqual(self.mcs.fragment_count_dict['ArCCCCR'], 999)
        self.assertEqual(self.mcs.fragment_count_dict['ArC__C'], 9)
        self.assertEqual(self.mcs.fragment_count_dict['ArC*CCCR'], 1)
        self.assertEqual(self.mcs.fragment_count_dict['ArC*C'], 1)
