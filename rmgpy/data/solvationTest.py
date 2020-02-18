#!/usr/bin/env python3

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

import os
from unittest import TestCase, TestLoader, TextTestRunner

from rmgpy import settings
from rmgpy.data.solvation import DatabaseError, SoluteData, SolvationDatabase, SolventLibrary
from rmgpy.molecule import Molecule
from rmgpy.rmg.main import RMG
from rmgpy.rmg.main import Species
from rmgpy.exceptions import InputError

###################################################

class TestSoluteDatabase(TestCase):

    def setUp(self):
        self.database = SolvationDatabase()
        self.database.load(os.path.join(settings['database.directory'], 'solvation'))

    def tearDown(self):
        """
        Reset the database & liquid parameters for solution
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

    def test_solute_library(self):
        """Test we can obtain solute parameters from a library"""
        species = Species(molecule=[Molecule(smiles='COC=O')])  # methyl formate - we know this is in the solute library

        library_data = self.database.get_solute_data_from_library(species, self.database.libraries['solute'])
        self.assertEqual(len(library_data), 3)

        solute_data = self.database.get_solute_data(species)
        self.assertTrue(isinstance(solute_data, SoluteData))

        s = solute_data.S
        self.assertEqual(s, 0.68)
        self.assertTrue(solute_data.V is not None)

    def test_mcgowan(self):
        """Test we can calculate and set the McGowan volume for species containing H,C,O,N or S"""
        self.testCases = [
            ['CCCCCCCC', 1.2358],  # n-octane, in library
            ['C(CO)O', 0.5078],  # ethylene glycol
            ['CC#N', 0.4042],  # acetonitrile
            ['CCS', 0.5539]  # ethanethiol
        ]

        for smiles, volume in self.testCases:
            species = Species(molecule=[Molecule(smiles=smiles)])
            solute_data = self.database.get_solute_data(species)
            solute_data.set_mcgowan_volume(species)  # even if it was found in library, recalculate
            self.assertIsNotNone(solute_data.V)  # so if it wasn't found in library, we should have calculated it
            self.assertAlmostEqual(solute_data.V, volume)  # the volume is what we expect given the atoms and bonds

    def test_diffusivity(self):
        """Test that for a given solvent viscosity and temperature we can calculate a solute's diffusivity"""
        species = Species(molecule=[Molecule(smiles='O')])  # water
        solute_data = self.database.get_solute_data(species)
        temperature = 298.
        solvent_viscosity = 0.00089  # water is about 8.9e-4 Pa.s
        d = solute_data.get_stokes_diffusivity(temperature, solvent_viscosity)  # m2/s
        self.assertAlmostEqual((d * 1e9), 1.3, 1)
        # self-diffusivity of water is about 2e-9 m2/s

    def test_solvent_library(self):
        """Test we can obtain solvent parameters from a library"""
        solvent_data = self.database.get_solvent_data('water')
        self.assertIsNotNone(solvent_data)
        self.assertEqual(solvent_data.s_h, 2.836)
        self.assertRaises(DatabaseError, self.database.get_solvent_data, 'orange_juice')
        solvent_data = self.database.get_solvent_data('cyclohexane')
        self.assertEqual(solvent_data.name_in_coolprop, 'CycloHexane')

    def test_viscosity(self):
        """Test we can calculate the solvent viscosity given a temperature and its A-E correlation parameters"""
        solvent_data = self.database.get_solvent_data('water')
        self.assertAlmostEqual(solvent_data.get_solvent_viscosity(298), 0.0009155)

    def test_critical_temperature(self):
        """
        Test we can calculate the solvent critical temperature given the solvent's name_in_coolprop
        and we can raise DatabaseError when the solvent's name_in_coolprop is None.
        """
        solvent_data = self.database.get_solvent_data('water')
        self.assertAlmostEqual(solvent_data.get_solvent_critical_temperature(), 647.096)
        solvent_data = self.database.get_solvent_data('dibutylether')
        self.assertRaises(DatabaseError, solvent_data.get_solvent_critical_temperature)

    def test_solute_generation(self):
        """Test we can estimate Abraham solute parameters correctly using group contributions"""

        self.testCases = [
            ['1,2-ethanediol', 'C(CO)O', 0.823, 0.685, 0.327, 2.572, 0.693, None],
        ]

        for name, smiles, S, B, E, L, A, V in self.testCases:
            species = Species(smiles=smiles)
            solute_data = self.database.get_solute_data_from_groups(species)
            self.assertAlmostEqual(solute_data.S, S, places=2)
            self.assertAlmostEqual(solute_data.B, B, places=2)
            self.assertAlmostEqual(solute_data.E, E, places=2)
            self.assertAlmostEqual(solute_data.L, L, places=2)
            self.assertAlmostEqual(solute_data.A, A, places=2)

    def test_solute_with_resonance_structures(self):
        """
        Test we can estimate Abraham solute parameters correctly using group contributions
        for the solute species with resonance structures.
        """
        smiles = "CC1=CC=CC=C1N"
        species = Species(smiles=smiles)
        species.generate_resonance_structures()
        solute_data = self.database.get_solute_data(species)
        solvent_data = self.database.get_solvent_data('water')
        solvation_correction = self.database.get_solvation_correction(solute_data, solvent_data)
        dGsolv_spc = solvation_correction.gibbs / 1000
        for mol in species.molecule:
            spc = Species(molecule=[mol])
            solute_data = self.database.get_solute_data_from_groups(spc)
            solvation_correction = self.database.get_solvation_correction(solute_data, solvent_data)
            dGsolv_mol = solvation_correction.gibbs / 1000
            if mol == species.molecule[0]:
                self.assertEqual(dGsolv_spc, dGsolv_mol)
            else:
                self.assertNotAlmostEqual(dGsolv_spc, dGsolv_mol)

    def test_lone_pair_solute_generation(self):
        """Test we can obtain solute parameters via group additivity for a molecule with lone pairs"""
        molecule = Molecule().from_adjacency_list(
            """
            CH2_singlet
            multiplicity 1
            1 C u0 p1 c0 {2,S} {3,S}
            2 H u0 p0 c0 {1,S}
            3 H u0 p0 c0 {1,S}
            """)
        species = Species(molecule=[molecule])
        solute_data = self.database.get_solute_data_from_groups(species)
        self.assertIsNotNone(solute_data)

    def test_solute_data_generation_ammonia(self):
        """Test we can obtain solute parameters via group additivity for ammonia"""
        molecule = Molecule().from_adjacency_list(
            """
            1 N u0 p1 c0 {2,S} {3,S} {4,S}
            2 H u0 p0 c0 {1,S}
            3 H u0 p0 c0 {1,S}
            4 H u0 p0 c0 {1,S}
            """)
        species = Species(molecule=[molecule])
        solute_data = self.database.get_solute_data_from_groups(species)
        self.assertIsNotNone(solute_data)

    def test_solute_data_generation_amide(self):
        """Test that we can obtain solute parameters via group additivity for an amide"""
        molecule = Molecule().from_adjacency_list(
            """
            1 N u0 p1 {2,S} {3,S} {4,S}
            2 H u0 {1,S}
            3 C u0 {1,S} {6,S} {7,S} {8,S}
            4 C u0 {1,S} {5,D} {9,S}
            5 O u0 p2 {4,D}
            6 H u0 {3,S}
            7 H u0 {3,S}
            8 H u0 {3,S}
            9 H u0 {4,S}
            """)
        species = Species(molecule=[molecule])
        solute_data = self.database.get_solute_data_from_groups(species)
        self.assertIsNotNone(solute_data)

    def test_solute_data_generation_co(self):
        """Test that we can obtain solute parameters via group additivity for CO."""
        molecule = Molecule().from_adjacency_list(
            """
            1  C u0 p1 c-1 {2,T}
            2  O u0 p1 c+1 {1,T}
            """)
        species = Species(molecule=[molecule])
        solute_data = self.database.get_solute_data_from_groups(species)
        self.assertIsNotNone(solute_data)

    def test_radical_and_lone_pair_generation(self):
        """
        Test we can obtain solute parameters via group additivity for a molecule with both lone 
        pairs and a radical
        """
        molecule = Molecule().from_adjacency_list(
            """
            [C]OH
            multiplicity 2
            1 C u1 p1 c0 {2,S}
            2 O u0 p2 c0 {1,S} {3,S}
            3 H u0 p0 c0 {2,S}
            """)
        species = Species(molecule=[molecule])
        solute_data = self.database.get_solute_data_from_groups(species)
        self.assertIsNotNone(solute_data)

    def test_radical_solute_group(self):
        """Test that the existing radical group is found for the radical species when using group additivity"""
        # First check whether the radical group is found for the radical species
        rad_species = Species(smiles='[OH]')
        rad_solute_data = self.database.get_solute_data_from_groups(rad_species)
        self.assertTrue('radical' in rad_solute_data.comment)
        # Then check that the radical and its saturated species give different solvation free energies
        saturated_struct = rad_species.molecule[0].copy(deep=True)
        saturated_struct.saturate_radicals()
        sat_species = Species(molecule=[saturated_struct])
        sat_solute_data = self.database.get_solute_data_from_groups(sat_species)
        solvent_data = self.database.get_solvent_data('water')
        rad_solvation_correction = self.database.get_solvation_correction(rad_solute_data, solvent_data)
        sat_solvation_correction = self.database.get_solvation_correction(sat_solute_data, solvent_data)
        self.assertNotAlmostEqual(rad_solvation_correction.gibbs / 1000, sat_solvation_correction.gibbs / 1000)

    def test_correction_generation(self):
        """Test we can estimate solvation thermochemistry."""
        self.testCases = [
            # solventName, soluteName, soluteSMILES, Hsolv, Gsolv
            ['water', 'acetic acid', 'C(C)(=O)O', -56500, -6700 * 4.184],
            ['water', 'naphthalene', 'C1=CC=CC2=CC=CC=C12', -42800, -2390 * 4.184],
            ['1-octanol', 'octane', 'CCCCCCCC', -40080, -4180 * 4.184],
            ['1-octanol', 'tetrahydrofuran', 'C1CCOC1', -28320, -3930 * 4.184],
            ['benzene', 'toluene', 'C1(=CC=CC=C1)C', -37660, -5320 * 4.184],
            ['benzene', '1,4-dioxane', 'C1COCCO1', -39030, -5210 * 4.184]
        ]

        for solventName, soluteName, smiles, H, G in self.testCases:
            species = Species(molecule=[Molecule(smiles=smiles)])
            solute_data = self.database.get_solute_data(species)
            solvent_data = self.database.get_solvent_data(solventName)
            solvation_correction = self.database.get_solvation_correction(solute_data, solvent_data)
            self.assertAlmostEqual(solvation_correction.enthalpy / 10000., H / 10000., 0,  # 0 decimal place, in 10kJ.
                                   msg="Solvation enthalpy discrepancy ({2:.0f}!={3:.0f}) for {0} in {1}"
                                       "".format(soluteName, solventName, solvation_correction.enthalpy, H))
            self.assertAlmostEqual(solvation_correction.gibbs / 10000., G / 10000., 0,
                                   msg="Solvation Gibbs free energy discrepancy ({2:.0f}!={3:.0f}) for {0} in {1}"
                                       "".format(soluteName, solventName, solvation_correction.gibbs, G))

    def test_Kfactor_parameters(self):
        """Test we can calculate the parameters for K-factor relationships"""
        species = Species().from_smiles('CCC(C)=O') # 2-Butanone for a solute
        solute_data = self.database.get_solute_data(species)
        solvent_data = self.database.get_solvent_data('water')
        kfactor_parameters = self.database.get_Kfactor_parameters(solute_data, solvent_data)
        self.assertAlmostEqual(kfactor_parameters.lower_T[0], -16.303, 3) # check up to 3 decimal places
        self.assertAlmostEqual(kfactor_parameters.lower_T[1], -0.930, 3)
        self.assertAlmostEqual(kfactor_parameters.lower_T[2], 17.550, 3)
        self.assertAlmostEqual(kfactor_parameters.higher_T, 1.308, 3)
        self.assertAlmostEqual(kfactor_parameters.T_transition, 485.3, 1)
        # check that DatabaseError is raised when the solvent's name_in_coolprop is None
        solvent_data = self.database.get_solvent_data('chloroform')
        self.assertRaises(DatabaseError, self.database.get_Kfactor_parameters, solute_data, solvent_data)

    def test_Tdep_solvation_calculation(self):
        '''Test we can calculate the temperature dependent K-factor and solvation free energy'''
        species = Species().from_smiles('CCC1=CC=CC=C1')  # ethylbenzene
        species.generate_resonance_structures()
        solute_data = self.database.get_solute_data(species)
        solvent_data = self.database.get_solvent_data('benzene')
        T = 500 # in K
        Kfactor = self.database.get_Kfactor(solute_data, solvent_data, T)
        delG = self.database.get_T_dep_solvation_energy(solute_data, solvent_data, T) / 1000 # in kJ/mol
        self.assertAlmostEqual(Kfactor, 0.416, 3)
        self.assertAlmostEqual(delG, -13.463, 3)
        # For temperature greater than or equal to the critical temperature of the solvent,
        # it should raise InputError
        T = 1000
        self.assertRaises(InputError, self.database.get_T_dep_solvation_energy, solute_data, solvent_data, T)

    def test_initial_species(self):
        """Test we can check whether the solvent is listed as one of the initial species in various scenarios"""

        # Case 1. when SMILES for solvent is available, the molecular structures of the initial species and the solvent
        # are compared to check whether the solvent is in the initial species list

        # Case 1-1: the solvent water is not in the initialSpecies list, so it raises Exception
        rmg = RMG()
        rmg.initial_species = []
        solute = Species(label='n-octane', molecule=[Molecule().from_smiles('C(CCCCC)CC')])
        rmg.initial_species.append(solute)
        rmg.solvent = 'water'
        solvent_structure = Species().from_smiles('O')
        self.assertRaises(Exception, self.database.check_solvent_in_initial_species, rmg, solvent_structure)

        # Case 1-2: the solvent is now octane and it is listed as the initialSpecies. Although the string
        # names of the solute and the solvent are different, because the solvent SMILES is provided,
        # it can identify the 'n-octane' as the solvent
        rmg.solvent = 'octane'
        solvent_structure = Species().from_smiles('CCCCCCCC')
        self.database.check_solvent_in_initial_species(rmg, solvent_structure)
        self.assertTrue(rmg.initial_species[0].is_solvent)

        # Case 2: the solvent SMILES is not provided. In this case, it can identify the species as the
        # solvent by looking at the string name.

        # Case 2-1: Since 'n-octane and 'octane' are not equal, it raises Exception
        solvent_structure = None
        self.assertRaises(Exception, self.database.check_solvent_in_initial_species, rmg, solvent_structure)

        # Case 2-2: The label 'n-ocatne' is corrected to 'octane', so it is identified as the solvent
        rmg.initial_species[0].label = 'octane'
        self.database.check_solvent_in_initial_species(rmg, solvent_structure)
        self.assertTrue(rmg.initial_species[0].is_solvent)

    def test_solvent_molecule(self):
        """Test that we can assign a proper solvent molecular structure when different formats are given"""

        # solventlibrary.entries['solvent_label'].item should be the instance of Species with the solvent's molecular
        # structure if the solvent database contains the solvent SMILES or adjacency list. If not, then item is None

        # Case 1: When the solventDatabase does not contain the solvent SMILES, the item attribute is None
        solventlibrary = SolventLibrary()
        solventlibrary.load_entry(index=1, label='water', solvent=None)
        self.assertTrue(solventlibrary.entries['water'].item is None)

        # Case 2: When the solventDatabase contains the correct solvent SMILES, the item attribute is the instance of
        # Species with the correct solvent molecular structure
        solventlibrary.load_entry(index=2, label='octane', solvent=None, molecule='CCCCCCCC')
        solvent_species = Species().from_smiles('C(CCCCC)CC')
        self.assertTrue(solvent_species.is_isomorphic(solventlibrary.entries['octane'].item[0]))

        # Case 3: When the solventDatabase contains the correct solvent adjacency list, the item attribute
        # is the instance of the species with the correct solvent molecular structure.
        # This will display the SMILES Parse Error message from the external function, but ignore it.
        solventlibrary.load_entry(index=3, label='ethanol', solvent=None, molecule="""
        1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
        2 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
        3 O u0 p2 c0 {2,S} {9,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {1,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {2,S}
        9 H u0 p0 c0 {3,S}
        """)
        solvent_species = Species().from_smiles('CCO')
        self.assertTrue(solvent_species.is_isomorphic(solventlibrary.entries['ethanol'].item[0]))

        # Case 4: when the solventDatabase contains incorrect values for the molecule attribute, it raises Exception
        # This will display the SMILES Parse Error message from the external function, but ignore it.
        self.assertRaises(Exception, solventlibrary.load_entry, index=4, label='benzene', solvent=None, molecule='ring')

        # Case 5: when the solventDatabase contains data for co-solvents.
        solventlibrary.load_entry(index=5, label='methanol_50_water_50', solvent=None, molecule=['CO', 'O'])
        solvent_species_list = [Species().from_smiles('CO'), Species().from_smiles('O')]
        self.assertEqual(len(solventlibrary.entries['methanol_50_water_50'].item), 2)
        for spc1 in solventlibrary.entries['methanol_50_water_50'].item:
            self.assertTrue(any([spc1.is_isomorphic(spc2) for spc2 in solvent_species_list]))


#####################################################


if __name__ == '__main__':
    suite = TestLoader().loadTestsFromTestCase(TestSoluteDatabase)
    TextTestRunner(verbosity=2).run(suite)
