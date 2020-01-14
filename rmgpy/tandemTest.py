#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains unit tests of the tandem root module.
"""

import os
import shutil
import unittest

import tandem
from rmgpy import settings
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import NASAPolynomial, NASA, ThermoData, Wilhoit
from rmgpy.thermo.model import HeatCapacityModel

################################################################################


class TestTandem(unittest.TestCase):
    """
    Contains unit tests for the RMG-ARC Tandem tool
    """

    @classmethod
    def setUpClass(cls):
        """
        A method that is run ONCE before all unit tests in this class.
        """
        cls.maxDiff = None
        cls.base_path = os.path.join(settings['test_data.directory'], 'tandem')
        cls.arc_input_file_path = os.path.join(cls.base_path, 'tandem_1.yml')
        cls.rmg_input_file_path = os.path.join(cls.base_path, 'rmg', 'input.py')

        cls.spc1 = Species().from_smiles('CC')
        cls.spc1.thermo = HeatCapacityModel()
        cls.spc1.thermo.comment = 'Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH)' \
                                  ' + group(Cs-CsHHH) + radical(RCCJ)'

        cls.spc2 = Species().from_smiles('CC')
        cls.spc2.thermo = HeatCapacityModel()
        cls.spc2.thermo.comment = 'Thermo library: primaryThermoLibrary + radical(CH3)'

        cls.spc3 = Species().from_smiles('CCO')
        cls.spc3.thermo = HeatCapacityModel()
        cls.spc3.thermo.comment = 'Thermo library: primaryThermoLibrary'

        cls.species_to_calc0 = [Species(label='i-C3H7', smiles='C[CH]C'),
                                Species(label='n-C3H7', smiles='[CH2]CC')]
        cls.all_species0 = list()

        cls.species_to_calc1 = [Species(label='C4H9a', smiles='[CH2]CCC'),
                                Species(label='C4H9b', smiles='C[CH]CC'),
                                Species(label='NH3', smiles='N')]
        cls.all_species1 = [Species(label='C4H9_0', smiles='C[C](C)(C)'),
                            Species(label='H2O', smiles='O')]

    def test_run_rmg(self):
        """Test the ability to run RMG from The Tandem Tool"""
        run_rmg_path = os.path.join(self.base_path, 'rmg', 'run_rmg')
        if os.path.isdir(run_rmg_path):
            shutil.rmtree(run_rmg_path)
        os.makedirs(run_rmg_path)
        rmg_input_file_path = os.path.join(run_rmg_path, 'input.py')
        shutil.copyfile(src=self.rmg_input_file_path, dst=rmg_input_file_path)
        tandem.run_rmg(input_file=rmg_input_file_path,
                       output_directory=run_rmg_path,
                       kwargs={'restart': '',
                               'walltime': '00:00:00:00',
                               'maxproc': 1,
                               'kineticsdatastore': False},
                       arguments={'max RMG walltime': '00:00:01:00',
                                  'max RMG exceptions allowed': 0},
                       tolerance=0.01,
                       thermo_library=None,
                       verbose=False)
        with open(os.path.join(run_rmg_path, 'RMG.log'), 'r') as f:
            line = f.readline()
        self.assertIn('RMG execution initiated', line)
        shutil.rmtree(run_rmg_path)

    def test_run_arc(self):
        """Test the ability to run ARC from The Tandem Tool"""
        tandem.run_arc(input_dict={'level_of_theory': 'b3lyp/6-31g**//b3lyp/6-31g**',
                                   'calc_freq_factor': False},
                       run_directory=self.base_path,
                       species_to_calc=list(),
                       verbose=False)
        with open(os.path.join(self.base_path, 'ARC', 'arc.log'), 'r') as f:
            line = f.readline()
        self.assertIn('ARC execution initiated', line)
        shutil.rmtree(os.path.join(self.base_path, 'ARC'))

    def test_parse_arc_input_file(self):
        """Test parsing the arc input file"""
        arguments, input_dict = tandem.parse_arc_input_file(self.arc_input_file_path, has_sa=True, has_pdep=False)
        expected_arguments = {'SA observables': [{'label': 'OH', 'smiles': '[OH]'}],
                              'SA method': 'RMG',
                              'SA threshold': 0.001,
                              'SA species': 10,
                              'SA reactions': 10,
                              'SA pdep threshold': 0.1,
                              'collision violators': True,
                              'all core species': True,
                              'RMG tolerances': [0.1, 0.01],
                              'max tandem iterations': 10,
                              'max RMG exceptions allowed': 10,
                              'max RMG walltime': '01:00:00:00'}
        expected_input_dict = {'level_of_theory': 'b3lyp/6-31g**//b3lyp/6-31g**',
                               'job_types': {'rotors': False, 'conformers': True, 'fine': False, 'freq': True,
                                             'opt': True, 'sp': True, 'onedmin': False, 'orbitals': False},
                               'allow_nonisomorphic_2d': True}
        self.assertEqual(arguments, expected_arguments)
        self.assertEqual(input_dict, expected_input_dict)

    def test_set_legal_species_labels(self):
        """Test setting legal species labels"""
        # test two species with the same formula
        updated_species_to_calc = tandem.set_legal_species_labels(self.species_to_calc0, self.all_species0)
        updated_labels = [spc.label for spc in updated_species_to_calc]
        self.assertEqual(updated_labels, ['C3H7_0', 'C3H7_1'])

        # test having a species with the same formula in all_species
        updated_species_to_calc = tandem.set_legal_species_labels(self.species_to_calc1, self.all_species1)
        updated_labels = [spc.label for spc in updated_species_to_calc]
        self.assertEqual(updated_labels, ['C4H9_1', 'C4H9_2', 'H3N_0'])

    def test_get_species_label_by_structure(self):
        """Test getting the specie slabel attribute from a list by its structure"""
        adj = """1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}"""
        spc_list = [Species(label='propane').from_smiles('CCC'), Species(label='ethane').from_smiles('CC')]
        label = tandem.get_species_label_by_structure(adj=adj, species_list=spc_list)
        self.assertEqual(label, 'ethane')

    def test_species_not_in_list(self):
        """Test determining whether a species is NOT in a list of species"""
        # calling set_legal_species_labels() sets the global species_labels_dict
        updated_species_to_calc = tandem.set_legal_species_labels(self.species_to_calc1, all_species=list())
        species_not_in_list = tandem.species_not_in_list('C4H9_1', updated_species_to_calc)
        self.assertFalse(species_not_in_list)
        species_not_in_list = tandem.species_not_in_list('C4H9_5', updated_species_to_calc)
        self.assertTrue(species_not_in_list)

    def test_get_species_by_label(self):
        """Test getting a species from a list by its label"""
        spc = tandem.get_species_by_label(label='C4H9a', species_list=self.species_to_calc1)
        self.assertIsInstance(spc, Species)
        self.assertEqual(spc.label, 'C4H9a')

        spc = tandem.get_species_by_label(label='x', species_list=self.species_to_calc1)
        self.assertIsNone(spc)

    def test_get_reaction_by_index(self):
        """Test getting a reaction from a list by its index"""
        reaction_list = [Reaction(index=0,
                                  reactants=[Species().from_smiles('[CH2]CC')],
                                  products=[Species().from_smiles('C[CH]C')])]
        rxn = tandem.get_reaction_by_index(index=0, reaction_list=reaction_list)
        self.assertIsInstance(rxn, Reaction)
        self.assertEqual(rxn.index, 0)
        self.assertEqual(str(rxn), '[CH2]CC <=> C[CH]C')

        rxn = tandem.get_reaction_by_index(index=5, reaction_list=reaction_list)
        self.assertIsNone(rxn)

    def test_calc_based_on_thermo_comment(self):
        """Test which species are selected for calculation based on their thermo comment"""
        self.assertTrue(tandem.calc_based_on_thermo_comment(self.spc1))
        self.assertTrue(tandem.calc_based_on_thermo_comment(self.spc2))
        self.assertFalse(tandem.calc_based_on_thermo_comment(self.spc3))

    def test_has_high_uncertainty(self):
        """Test determining whether a species thermo should be calculated"""
        should_species_be_calculated = tandem.has_high_uncertainty(
            species=self.spc1, unconverged_species=list(), species_to_calc=dict())
        self.assertTrue(should_species_be_calculated)

        should_species_be_calculated = tandem.has_high_uncertainty(
            species=self.spc1, unconverged_species=[self.spc1.copy()], species_to_calc=dict())
        self.assertFalse(should_species_be_calculated)

        should_species_be_calculated = tandem.has_high_uncertainty(
            species=self.spc1, unconverged_species=list(), species_to_calc={'label': {'spc': self.spc2}})
        self.assertFalse(should_species_be_calculated)

        should_species_be_calculated = tandem.has_high_uncertainty(
            species=self.spc1, unconverged_species=list(), species_to_calc={'label': {'spc': self.spc3}})
        self.assertTrue(should_species_be_calculated)

    def test_load_species_and_reactions_from_chemkin_file(self):
        """Test loading species and reactions from a Chemkin file"""
        run_directory = os.path.join(self.base_path, 'iteration_0')
        rmg_species, rmg_reactions = tandem.load_species_and_reactions_from_chemkin_file(run_directory)
        self.assertEqual(len(rmg_species), 27)
        self.assertEqual(len(rmg_reactions), 227)
        self.assertIsInstance(rmg_species[0], Species)
        self.assertIsInstance(rmg_reactions[0], Reaction)

    # def test_determine_species_to_calculate(self):
    #     """Test that we correctly determine the species to be calculated from an RMG job"""
    #     # TODO: add an actual case with SA and PDep and coll violators
    #     pass

    def test_determine_species_based_on_sensitivity(self):
        """Test determining species to calculate based on sensitivity analysis"""
        run_directory = os.path.join(self.base_path, 'iteration_1')
        arguments = tandem.parse_arc_input_file(self.arc_input_file_path, has_sa=True, has_pdep=False)[0]
        rmg_species, rmg_reactions = tandem.load_species_and_reactions_from_chemkin_file(run_directory)
        unconverged_species = list()
        species_to_calc = tandem.determine_species_based_on_sensitivity(
            run_directory, arguments, rmg_species, rmg_reactions, unconverged_species,
            iteration=1, executed_networks=list(), verbose=False)[0]
        self.assertEqual(len(list(species_to_calc.values())), 7)
        species_to_calc_str = tandem.dict_to_str(species_to_calc)
        expected_species_to_calc_str = """ethane(1):
  spc: ethane(1)
  reason: observable
C2H4(17):
  spc: C=C(17)
  reason: (iteration 1) participates in the 2nd most sensitive reaction for ethane(1): C=C(17) + H(5) <=> C[CH2](4)
C2H5(4):
  spc: C[CH2](4)
  reason: (iteration 1) participates in the 2nd most sensitive reaction for ethane(1): C=C(17) + H(5) <=> C[CH2](4)
HO2(8):
  spc: [O]O(8)
  reason: (iteration 1) participates in the 5th most sensitive reaction for ethane(1): H(5) + O2(2) <=> [O]O(8)
C2H4O(41):
  spc: C=CO(41)
  reason: (iteration 1) participates in the 6th most sensitive reaction for ethane(1): C=C(17) + O(T)(14) <=> C=CO(41)
OO(34):
  spc: OO(34)
  reason: (iteration 1) participates in the 8th most sensitive reaction for ethane(1): OH(D)(33) + OH(D)(33) <=> OO(34)
CH3(3):
  spc: [CH3](3)
  reason: (iteration 1) the 5th most sensitive species thermo for ethane(1)
"""
        self.assertEqual(species_to_calc_str, expected_species_to_calc_str)

#     def test_determine_species_from_pdep_network(self):
#         """"""
    # use /home/alongd/Code/RMG-Py/rmgpy/test_data/tandem/sa_coefficients.yml with '+' in species names
#         ch2_adj = """1 C u0 p1 c0 {2,S} {3,S}
# 2 H u0 p0 c0 {1,S}
# 3 H u0 p0 c0 {1,S}"""
#         rmg_species = [Species(label='CO[O](9)').from_smiles('CO[O]'),
#                        Species(label='[CH2]OO(10)').from_smiles('[CH2]OO'),
#                        Species(label='O2(2)').from_smiles('[O][O]'),
#                        Species(label='CH3_0(5)').from_smiles('[CH3]'),
#                        Species(label='OH(D)(27)').from_smiles('[OH]'),
#                        Species(label='C=O(26)').from_smiles('C=O'),
#                        Species(label='[O]CO(28)').from_smiles('[O]CO'),
#                        Species(label='[O]O(8)').from_smiles('[O]O'),
#                        Species(label='CH2(S)(3)').from_adjacency_list(ch2_adj),
#                        Species(label='N2').from_smiles('N#N'),
#                        Species(label='Ar').from_smiles('[Ar]'),
#                        Species(label='He').from_smiles('[He]'),
#                        Species(label='Ne').from_smiles('[Ne]')]
#         pdep_rxn = PDepReaction(index=0,
#                                 reactants=[Species(label='[O]O(8)').from_smiles('[O]O'),
#                                            Species(label='CH2(S)(3)').from_adjacency_list(ch2_adj)],
#                                 products=[Species(label='CO[O](9)').from_smiles('CO[O]')],
#                                 network=PDepNetwork(index=27))
#         pdep_rxns_to_explore = [(pdep_rxn, 1, 'CH2(S)(3)')]
#         species_to_calc, executed_networks = \
#             tandem.determine_species_from_pdep_network(
#                 run_directory=os.path.join(self.base_path, 'iteration_0'),
#                 pdep_rxns_to_explore=pdep_rxns_to_explore,
#                 unconverged_species=list(),
#                 species_to_calc=dict(),
#                 iteration=1,
#                 threshold=0.1,
#                 executed_networks=list(),
#                 rmg_species=rmg_species,
#                 verbose=False)
#         # print(species_to_calc)
#         # print(executed_networks)
#         # raise

    def test_modify_pdep_network_file(self):
        """Test modifying an Arkane P-dep network input file"""
        pdep_sa_path = os.path.join(self.base_path, 'iteration_0', 'pdep_sa')
        if os.path.isdir(pdep_sa_path):
            shutil.rmtree(pdep_sa_path)

        input_file_path, output_file_path, isomer_labels = \
            tandem.modify_pdep_network_file(run_directory=os.path.join(self.base_path, 'iteration_0'),
                                            network_name='network27_2', method='CSE')
        self.assertIn('tandem/iteration_0/pdep_sa/network27_2/CSE/input.py', input_file_path)
        self.assertIn('tandem/iteration_0/pdep_sa/network27_2/CSE/sensitivity/sa_coefficients.yml', output_file_path)
        self.assertEqual(isomer_labels, ('CO[O](9)', '[CH2]OO(10)'))

        with open(input_file_path, 'r') as f:
            lines = f.readlines()
        sensitivity_conditions, cse = False, False
        for line in lines:
            if "    sensitivity_conditions = [[(300, 'K'), (0.01, 'bar')]," in line:
                sensitivity_conditions = True
            elif "    method = 'chemically-significant eigenvalues'," in line:
                cse = True
        self.assertTrue(sensitivity_conditions)
        self.assertTrue(cse)

        input_file_path, output_file_path, isomer_labels = \
            tandem.modify_pdep_network_file(run_directory=os.path.join(self.base_path, 'iteration_0'),
                                            network_name='network3_1', method='MSC')
        self.assertIn('tandem/iteration_0/pdep_sa/network3_1/MSC/input.py', input_file_path)
        self.assertIn('tandem/iteration_0/pdep_sa/network3_1/MSC/sensitivity/sa_coefficients.yml', output_file_path)
        self.assertEqual(isomer_labels, ('ethane(1)',))

        with open(input_file_path, 'r') as f:
            lines = f.readlines()
        sensitivity_conditions, msc, wrong_msc = False, False, False
        for line in lines:
            if "    sensitivity_conditions = [[(300, 'K'), (0.01, 'bar')]," in line:
                sensitivity_conditions = True
            elif "    method = 'modified strong collision'," in line:
                msc = True
            elif "modified strong collision" in line:
                # this string must not appear in any other line
                wrong_msc = True
        self.assertTrue(sensitivity_conditions)
        self.assertTrue(msc)
        self.assertFalse(wrong_msc)

        shutil.rmtree(pdep_sa_path)

    def test_determine_species_based_on_collision_violators(self):
        """Test determining species to calculate based on collision rate violating reactions"""
        run_directory = os.path.join(self.base_path, 'iteration_2')
        rmg_species, rmg_reactions = tandem.load_species_and_reactions_from_chemkin_file(run_directory)
        unconverged_species = list()
        species_to_calc = tandem.determine_species_based_on_collision_violators(
            run_directory, rmg_species, unconverged_species)
        self.assertEqual(len(list(species_to_calc.values())), 1)
        species_to_calc_str = tandem.dict_to_str(species_to_calc)
        expected_species_to_calc_str = """C3H6(19):
  spc: [CH2]C[CH2](19)
  reason: species participates in a collision rate violating reaction, C3H6(19)+H(4)=C3H7_0(15)
"""
        self.assertEqual(species_to_calc_str, expected_species_to_calc_str)

    def test_add_rmg_libraries(self):
        """Test adding an RMG library to thr RMG database repository"""
        rmg_thermo_lib_1_path = os.path.join(settings['database.directory'], 'thermo', 'libraries', 't3_thermo.py')
        rmg_thermo_lib_2_path = os.path.join(settings['database.directory'], 'thermo', 'libraries', 't3_thermo_0.py')
        libraries_path = os.path.join(self.base_path, 'iteration_0')
        local_context = {'ThermoData': ThermoData, 'Wilhoit': Wilhoit, 'NASAPolynomial': NASAPolynomial, 'NASA': NASA}

        if os.path.isfile(rmg_thermo_lib_1_path):
            os.remove(rmg_thermo_lib_1_path)
        if os.path.isfile(rmg_thermo_lib_2_path):
            os.remove(rmg_thermo_lib_2_path)

        # test adding a library for the fist time
        library_name = tandem.add_rmg_libraries(run_directory=libraries_path, library_name=None, verbose=False)
        self.assertEqual(library_name, 't3_thermo')
        thermo_lib = ThermoLibrary()
        thermo_lib.load(path=rmg_thermo_lib_1_path, local_context=local_context, global_context=dict())
        self.assertEqual(len(list(thermo_lib.entries.values())), 10)

        # test adding a library for the fist time with a different name ('t3_thermo' is occupied)
        library_name = tandem.add_rmg_libraries(run_directory=libraries_path, library_name=None, verbose=False)
        self.assertEqual(library_name, 't3_thermo_0')
        thermo_lib = ThermoLibrary()
        thermo_lib.load(path=rmg_thermo_lib_1_path, local_context=local_context, global_context=dict())
        self.assertEqual(len(list(thermo_lib.entries.values())), 10)

        # test appending entries to an existing library
        libraries_path = os.path.join(self.base_path, 'iteration_1')
        library_name = tandem.add_rmg_libraries(run_directory=libraries_path, library_name='t3_thermo', verbose=False)
        self.assertEqual(library_name, 't3_thermo')
        thermo_lib = ThermoLibrary()
        thermo_lib.load(path=rmg_thermo_lib_1_path, local_context=local_context, global_context=dict())
        self.assertEqual(len(list(thermo_lib.entries.values())), 11)  # extended with one additional entry

        os.remove(rmg_thermo_lib_1_path)
        os.remove(rmg_thermo_lib_2_path)

    def test_get_unconverged_species(self):
        """Test attaining a list of unconverged species from the ARC project info file"""
        labels = ['C3H8_0', 'C3H7_1', 'C3H6_0', 'C3H6_1', 'C3H5_1',
                  'C3H5_2', 'C3H6_2', 'C3H4_0', 'C3H3_0', 'C3H4_1', 'C3H4_2']
        all_species = [Species(label=label) for label in labels]
        unconverged_species = tandem.get_unconverged_species(run_directory=os.path.join(self.base_path, 'iteration_0'),
                                                             all_species=all_species,
                                                             log_species=False)
        self.assertEqual(len(unconverged_species), 1)
        self.assertEqual(unconverged_species[0].label, 'C3H4_1')

    def test_dict_to_str(self):
        """Test prettifying a dictionary"""
        dictionary = {'label1': {'spc': 'Species1', 'reason': 'Reason1'},
                      'label2': {'spc': 'Species2', 'reason': 'Reason2'}}
        output = tandem.dict_to_str(dictionary)
        expected_output = """label1:
  spc: Species1
  reason: Reason1
label2:
  spc: Species2
  reason: Reason2
"""
        self.assertEqual(output, expected_output)

    def test_combine_dicts(self):
        """Test combining two dictionaries"""
        dict1 = {1: 7, 3: 9}
        dict2 = {5: 2, 4: 8, 1: 7}
        combined = tandem.combine_dicts(dict1=dict1, dict2=None)
        self.assertEqual(combined, dict1)
        combined = tandem.combine_dicts(dict1=dict1, dict2=dict())
        self.assertEqual(combined, dict1)
        combined = tandem.combine_dicts(dict1=None, dict2=dict1)
        self.assertEqual(combined, dict1)
        combined = tandem.combine_dicts(dict1=dict1, dict2=dict2)
        self.assertEqual(combined, {1: 7, 3: 9, 5: 2, 4: 8})

################################################################################


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
