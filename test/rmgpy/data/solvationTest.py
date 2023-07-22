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

import os
from unittest import TestCase, TestLoader, TextTestRunner

from rmgpy import settings
from rmgpy.data.solvation import (
    DatabaseError,
    SoluteData,
    SolvationDatabase,
    SolventLibrary,
    get_critical_temperature,
    get_liquid_saturation_density,
    get_gas_saturation_density,
)
from rmgpy.molecule import Molecule
from rmgpy.rmg.main import RMG
from rmgpy.rmg.main import Species
from rmgpy.exceptions import InputError

###################################################


class TestSoluteDatabase(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.database = SolvationDatabase()
        cls.database.load(os.path.join(settings["database.directory"], "solvation"))

    @classmethod
    def tearDownClass(cls):
        """
        Reset the database & liquid parameters for solution
        """
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None

    def test_solute_library(self):
        """Test we can obtain solute parameters from a library"""
        species = Species(
            molecule=[Molecule(smiles="COC=O")]
        )  # methyl formate - we know this is in the solute library

        library_data = self.database.get_solute_data_from_library(
            species, self.database.libraries["solute"]
        )
        self.assertEqual(len(library_data), 3)

        solute_data = self.database.get_solute_data(species)
        self.assertTrue(isinstance(solute_data, SoluteData))

        s = solute_data.S
        self.assertEqual(s, 0.68)
        self.assertTrue(solute_data.V is not None)

    def test_mcgowan(self):
        """Test we can calculate and set the McGowan volume for species containing H,C,O,N or S"""
        self.testCases = [
            ["CCCCCCCC", 1.2358],  # n-octane, in library
            ["C(CO)O", 0.5078],  # ethylene glycol
            ["CC#N", 0.4042],  # acetonitrile
            ["CCS", 0.5539],  # ethanethiol
        ]

        for smiles, volume in self.testCases:
            species = Species(molecule=[Molecule(smiles=smiles)])
            solute_data = self.database.get_solute_data(species)
            solute_data.set_mcgowan_volume(
                species
            )  # even if it was found in library, recalculate
            self.assertIsNotNone(
                solute_data.V
            )  # so if it wasn't found in library, we should have calculated it
            self.assertAlmostEqual(
                solute_data.V, volume
            )  # the volume is what we expect given the atoms and bonds

    def test_diffusivity(self):
        """Test that for a given solvent viscosity and temperature we can calculate a solute's diffusivity"""
        species = Species(molecule=[Molecule(smiles="O")])  # water
        solute_data = self.database.get_solute_data(species)
        temperature = 298.0
        solvent_viscosity = 0.00089  # water is about 8.9e-4 Pa.s
        d = solute_data.get_stokes_diffusivity(temperature, solvent_viscosity)  # m2/s
        self.assertAlmostEqual((d * 1e9), 1.3, 1)
        # self-diffusivity of water is about 2e-9 m2/s

    def test_solvent_library(self):
        """Test we can obtain solvent parameters and data count from a library"""
        solvent_data = self.database.get_solvent_data("water")
        self.assertIsNotNone(solvent_data)
        self.assertEqual(solvent_data.s_h, -0.75922)
        self.assertRaises(DatabaseError, self.database.get_solvent_data, "orange_juice")
        solvent_data = self.database.get_solvent_data("cyclohexane")
        self.assertEqual(solvent_data.name_in_coolprop, "CycloHexane")
        solvent_data_count = self.database.get_solvent_data_count("dodecan-1-ol")
        self.assertEqual(solvent_data_count.dGsolvCount, 11)
        dHsolvMAE = (0.05, "kcal/mol")
        self.assertTrue(solvent_data_count.dHsolvMAE == dHsolvMAE)

    def test_viscosity(self):
        """Test we can calculate the solvent viscosity given a temperature and its A-E correlation parameters"""
        solvent_data = self.database.get_solvent_data("water")
        self.assertAlmostEqual(solvent_data.get_solvent_viscosity(298), 0.0009155)

    def test_critical_temperature(self):
        """
        Test we can calculate the solvent critical temperature given the solvent's name_in_coolprop
        and we can raise DatabaseError when the solvent's name_in_coolprop is None.
        """
        solvent_data = self.database.get_solvent_data("water")
        solvent_name = solvent_data.name_in_coolprop
        self.assertAlmostEqual(get_critical_temperature(solvent_name), 647.096)
        solvent_data = self.database.get_solvent_data("dibutylether")
        solvent_name = solvent_data.name_in_coolprop
        self.assertRaises(DatabaseError, get_critical_temperature, solvent_name)

    def test_saturation_density(self):
        """
        Test we can calculate the solvent's liquid-phase and gas-phase saturation densities given the compound name
        and temperature and we can raise DatabaseError when the compound is not available in CoolProp or
        the temperature is out of the calculable range.
        """
        compound_name = "Hexane"
        temp = 400  # in K
        self.assertAlmostEqual(
            get_liquid_saturation_density(compound_name, temp), 6383.22, places=2
        )
        self.assertAlmostEqual(
            get_gas_saturation_density(compound_name, temp), 162.99, places=2
        )
        # Unsupported compound name
        self.assertRaises(DatabaseError, get_gas_saturation_density, "Hexadecane", temp)
        # Out of the valid temperature range
        self.assertRaises(DatabaseError, get_gas_saturation_density, compound_name, 700)

    def test_find_solvent(self):
        """Test we can find solvents from the solvent library using SMILES"""
        # Case 1: one solvent is matched
        solvent_smiles = "NC=O"
        match_list = self.database.find_solvent_from_smiles(solvent_smiles)
        self.assertEqual(len(match_list), 1)
        self.assertTrue(match_list[0][0] == "formamide")
        # Case 2: two solvents are matched
        solvent_smiles = "ClC=CCl"
        match_list = self.database.find_solvent_from_smiles(solvent_smiles)
        self.assertEqual(len(match_list), 2)
        self.assertTrue(match_list[0][0] == "cis-1,2-dichloroethene")
        self.assertTrue(match_list[1][0] == "trans-1,2-dichloroethene")
        # Case 3: no solvent is matched
        solvent_smiles = "C(CCl)O"
        match_list = self.database.find_solvent_from_smiles(solvent_smiles)
        self.assertEqual(len(match_list), 0)

    def test_solute_groups(self):
        """Test we can correctly load the solute groups from the solvation group database"""
        solute_group = self.database.groups["group"].entries["Cds-N3dCbCb"]
        self.assertEqual(solute_group.data_count.S, 28)
        self.assertEqual(solute_group.data.B, 0.06652)
        solute_group = self.database.groups["ring"].entries["FourMember"]
        self.assertIsNone(solute_group.data_count)
        self.assertEqual(solute_group.data, "Cyclobutane")

    def test_solute_generation(self):
        """Test we can estimate Abraham solute parameters correctly using group contributions"""

        self.testCases = [
            ["1,2-ethanediol", "C(CO)O", 0.809, 0.740, 0.393, 2.482, 0.584, 0.508]
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
        solvent_data = self.database.get_solvent_data("water")
        solvation_correction = self.database.get_solvation_correction(
            solute_data, solvent_data
        )
        dGsolv_spc = solvation_correction.gibbs / 1000
        for mol in species.molecule:
            spc = Species(molecule=[mol])
            solute_data = self.database.get_solute_data_from_groups(spc)
            solvation_correction = self.database.get_solvation_correction(
                solute_data, solvent_data
            )
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
            """
        )
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
            """
        )
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
            """
        )
        species = Species(molecule=[molecule])
        solute_data = self.database.get_solute_data_from_groups(species)
        self.assertIsNotNone(solute_data)

    def test_solute_data_generation_co(self):
        """Test that we can obtain solute parameters via group additivity for CO."""
        molecule = Molecule().from_adjacency_list(
            """
            1  C u0 p1 c-1 {2,T}
            2  O u0 p1 c+1 {1,T}
            """
        )
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
            """
        )
        species = Species(molecule=[molecule])
        solute_data = self.database.get_solute_data_from_groups(species)
        self.assertIsNotNone(solute_data)

    def test_radical_solute_group(self):
        """Test that the existing radical group is found for the radical species when using group additivity"""
        # First check whether the radical group is found for the radical species
        rad_species = Species(smiles="[OH]")
        rad_solute_data = self.database.get_solute_data_from_groups(rad_species)
        self.assertTrue("radical" in rad_solute_data.comment)
        # Then check that the radical and its saturated species give different solvation free energies
        saturated_struct = rad_species.molecule[0].copy(deep=True)
        saturated_struct.saturate_radicals()
        sat_species = Species(molecule=[saturated_struct])
        sat_solute_data = self.database.get_solute_data_from_groups(sat_species)
        solvent_data = self.database.get_solvent_data("water")
        rad_solvation_correction = self.database.get_solvation_correction(
            rad_solute_data, solvent_data
        )
        sat_solvation_correction = self.database.get_solvation_correction(
            sat_solute_data, solvent_data
        )
        self.assertNotAlmostEqual(
            rad_solvation_correction.gibbs / 1000, sat_solvation_correction.gibbs / 1000
        )

    def test_halogen_solute_group(self):
        """Test that the correct halogen groups can be found for the halogenated species using get_solute_data method"""
        # Check the species whose halogen-replaced form can be found from solute library
        species = Species().from_smiles("CCCCCCl")
        solute_data = self.database.get_solute_data(species)
        self.assertTrue(
            "Solute library: n-pentane + halogen(Cl-(Cs-CsHH))" in solute_data.comment
        )
        # Check the species whose halogen-replaced form cannot be found from solute library
        species = Species().from_smiles("OCCCCCCC(Br)CCCCCO")
        solute_data = self.database.get_solute_data(species)
        self.assertTrue(
            "+ group(Cs-Cs(Os-H)HH) + halogen(Br-(Cs-CsCsH))" in solute_data.comment
        )

    def test_radical_halogen_solute_group(self):
        """Test that the correct halogen and radical groups can be found for the halogenated radical species
        using get_solute_data method"""
        # Check the species whose saturated and halogenated form can be found from solute library
        species = Species().from_smiles("[O]CCCCl")
        solute_data = self.database.get_solute_data(species)
        self.assertTrue(
            "Solute library: 3-Chloropropan-1-ol + radical(ROJ)" == solute_data.comment
        )
        # Check the species whose saturated and halogen-replaced form can be found from solute library
        species = Species().from_smiles("[O]CCCC(Br)(I)Cl")
        solute_data = self.database.get_solute_data(species)
        self.assertTrue(
            "Solute library: butan-1-ol + halogen(I-(Cs-CsHH)) + halogen(Br-(Cs-CsFCl)) + halogen(Cl-(Cs-CsFBr)) + radical(ROJ)"
            == solute_data.comment
        )
        # Check the species whose saturated and halogen-replaced form cannot be found from solute library
        species = Species().from_smiles("[NH]C(=O)CCCl")
        solute_data = self.database.get_solute_data(species)
        self.assertTrue(
            "group(Cds-Od(N3s-HH)Cs) + halogen(Cl-(Cs-CsHH)) + radical(N3_amide_pri)"
            in solute_data.comment
        )
        # Check the species whose radical site is bonded to halogen
        species = Species().from_smiles("F[N]C(=O)CCCl")
        solute_data = self.database.get_solute_data(species)
        self.assertTrue(
            "group(Cds-Od(N3s-HH)Cs) + halogen(Cl-(Cs-CsHH)) + halogen(F-N3s) + radical(N3_amide_sec)"
            in solute_data.comment
        )

    def test_correction_generation(self):
        """Test we can estimate solvation thermochemistry."""
        self.testCases = [
            # solventName, soluteName, soluteSMILES, Hsolv, Gsolv in kJ/mol
            ["water", "acetic acid", "C(C)(=O)O", -48.48, -28.12],
            ["water", "naphthalene", "C1=CC=CC2=CC=CC=C12", -37.15, -11.21],
            ["1-octanol", "octane", "CCCCCCCC", -39.44, -16.83],
            ["1-octanol", "tetrahydrofuran", "C1CCOC1", -32.27, -17.81],
            ["benzene", "toluene", "C1(=CC=CC=C1)C", -39.33, -23.81],
            ["benzene", "1,4-dioxane", "C1COCCO1", -39.15, -22.01],
        ]

        for solventName, soluteName, smiles, H, G in self.testCases:
            species = Species().from_smiles(smiles)
            species.generate_resonance_structures()
            solute_data = self.database.get_solute_data(species)
            solvent_data = self.database.get_solvent_data(solventName)
            solvation_correction = self.database.get_solvation_correction(
                solute_data, solvent_data
            )
            self.assertAlmostEqual(
                solvation_correction.enthalpy / 1000,
                H,
                2,  # 2 decimal places, in kJ.
                msg="Solvation enthalpy discrepancy ({2:.2f}!={3:.2f}) for {0} in {1}"
                "".format(
                    soluteName, solventName, solvation_correction.enthalpy / 1000, H
                ),
            )
            self.assertAlmostEqual(
                solvation_correction.gibbs / 1000,
                G,
                2,  # 2 decimal places, in kJ.
                msg="Solvation Gibbs free energy discrepancy ({2:.2f}!={3:.2f}) for {0} in {1}"
                "".format(
                    soluteName, solventName, solvation_correction.gibbs / 1000, G
                ),
            )

    def test_Kfactor_parameters(self):
        """Test we can calculate the parameters for K-factor relationships"""
        species = Species().from_smiles("CCC(C)=O")  # 2-Butanone for a solute
        solute_data = self.database.get_solute_data(species)
        solvent_data = self.database.get_solvent_data("water")
        correction = self.database.get_solvation_correction(solute_data, solvent_data)
        delG298 = correction.gibbs  # in J/mol
        delH298 = correction.enthalpy  # in J/mol
        delS298 = correction.entropy  # in J/mol/K
        solvent_name = solvent_data.name_in_coolprop
        kfactor_parameters = self.database.get_Kfactor_parameters(
            delG298, delH298, delS298, solvent_name
        )
        self.assertAlmostEqual(
            kfactor_parameters.lower_T[0], -9.780, 3
        )  # check up to 3 decimal places
        self.assertAlmostEqual(kfactor_parameters.lower_T[1], 0.492, 3)
        self.assertAlmostEqual(kfactor_parameters.lower_T[2], 10.485, 3)
        self.assertAlmostEqual(kfactor_parameters.higher_T, 1.147, 3)
        self.assertAlmostEqual(kfactor_parameters.T_transition, 485.3, 1)
        # check that DatabaseError is raised when the solvent's name_in_coolprop is None
        solvent_data = self.database.get_solvent_data("chloroform")
        solvent_name = solvent_data.name_in_coolprop
        self.assertRaises(
            DatabaseError,
            self.database.get_Kfactor_parameters,
            delG298,
            delH298,
            delS298,
            solvent_name,
        )

    def test_Tdep_solvation_calculation(self):
        """
        Test we can calculate the temperature dependent solvation free energy and K-factor
        using both `get_T_dep_solvation_energy_from_LSER_298` and `get_T_dep_solvation_energy_from_input_298` methods.
        """
        # First, test `get_T_dep_solvation_energy_from_LSER_298` method.
        species = Species().from_smiles("CCC1=CC=CC=C1")  # ethylbenzene
        species.generate_resonance_structures()
        solute_data = self.database.get_solute_data(species)
        solvent_data = self.database.get_solvent_data("benzene")
        T = 500  # in K
        delG, Kfactor, kH = self.database.get_T_dep_solvation_energy_from_LSER_298(
            solute_data, solvent_data, T
        )
        self.assertAlmostEqual(Kfactor, 0.403, 3)
        self.assertAlmostEqual(delG / 1000, -13.59, 2)  # delG is in J/mol
        # For temperature greater than or equal to the critical temperature of the solvent,
        # it should raise InputError
        T = 1000
        self.assertRaises(
            InputError,
            self.database.get_T_dep_solvation_energy_from_LSER_298,
            solute_data,
            solvent_data,
            T,
        )

        # Now test `get_T_dep_solvation_energy_from_input_298` method.
        delG298 = -23570  # in J/mol
        delH298 = -40612  # in J/mol
        delS298 = (delH298 - delG298) / 298  # in J/mol/K
        solvent_name = "benzene"
        T = 500  # in K
        delG, Kfactor, kH = self.database.get_T_dep_solvation_energy_from_input_298(
            delG298, delH298, delS298, solvent_name, T
        )
        self.assertAlmostEqual(Kfactor, 0.567, 3)
        self.assertAlmostEqual(delG / 1000, -12.18, 2)  # delG is in J/mol
        # test that it raises InputError for T above the critical temperature
        T = 1000
        self.assertRaises(
            InputError,
            self.database.get_T_dep_solvation_energy_from_input_298,
            delG298,
            delH298,
            delS298,
            solvent_name,
            T,
        )

    def test_initial_species(self):
        """Test we can check whether the solvent is listed as one of the initial species in various scenarios"""

        # Case 1. when SMILES for solvent is available, the molecular structures of the initial species and the solvent
        # are compared to check whether the solvent is in the initial species list

        # Case 1-1: the solvent water is not in the initialSpecies list, so it raises Exception
        rmg = RMG()
        rmg.initial_species = []
        solute = Species(
            label="n-octane", molecule=[Molecule().from_smiles("C(CCCCC)CC")]
        )
        rmg.initial_species.append(solute)
        rmg.solvent = "water"
        solvent_structure = Species().from_smiles("O")
        self.assertRaises(
            Exception,
            self.database.check_solvent_in_initial_species,
            rmg,
            solvent_structure,
        )

        # Case 1-2: the solvent is now octane and it is listed as the initialSpecies. Although the string
        # names of the solute and the solvent are different, because the solvent SMILES is provided,
        # it can identify the 'n-octane' as the solvent
        rmg.solvent = "octane"
        solvent_structure = Species().from_smiles("CCCCCCCC")
        self.database.check_solvent_in_initial_species(rmg, solvent_structure)
        self.assertTrue(rmg.initial_species[0].is_solvent)

        # Case 2: the solvent SMILES is not provided. In this case, it can identify the species as the
        # solvent by looking at the string name.

        # Case 2-1: Since 'n-octane and 'octane' are not equal, it raises Exception
        solvent_structure = None
        self.assertRaises(
            Exception,
            self.database.check_solvent_in_initial_species,
            rmg,
            solvent_structure,
        )

        # Case 2-2: The label 'n-ocatne' is corrected to 'octane', so it is identified as the solvent
        rmg.initial_species[0].label = "octane"
        self.database.check_solvent_in_initial_species(rmg, solvent_structure)
        self.assertTrue(rmg.initial_species[0].is_solvent)

    def test_solvent_molecule(self):
        """Test that we can assign a proper solvent molecular structure when different formats are given"""

        # solventlibrary.entries['solvent_label'].item should be the instance of Species with the solvent's molecular
        # structure if the solvent database contains the solvent SMILES or adjacency list. If not, then item is None

        # Case 1: When the solventDatabase does not contain the solvent SMILES, the item attribute is None
        solventlibrary = SolventLibrary()
        solventlibrary.load_entry(index=1, label="water", solvent=None)
        self.assertTrue(solventlibrary.entries["water"].item is None)

        # Case 2: When the solventDatabase contains the correct solvent SMILES, the item attribute is the instance of
        # Species with the correct solvent molecular structure
        solventlibrary.load_entry(
            index=2, label="octane", solvent=None, molecule="CCCCCCCC"
        )
        solvent_species = Species().from_smiles("C(CCCCC)CC")
        self.assertTrue(
            solvent_species.is_isomorphic(solventlibrary.entries["octane"].item[0])
        )

        # Case 3: When the solventDatabase contains the correct solvent adjacency list, the item attribute
        # is the instance of the species with the correct solvent molecular structure.
        # This will display the SMILES Parse Error message from the external function, but ignore it.
        solventlibrary.load_entry(
            index=3,
            label="ethanol",
            solvent=None,
            molecule="""
        1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
        2 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
        3 O u0 p2 c0 {2,S} {9,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 H u0 p0 c0 {1,S}
        7 H u0 p0 c0 {2,S}
        8 H u0 p0 c0 {2,S}
        9 H u0 p0 c0 {3,S}
        """,
        )
        solvent_species = Species().from_smiles("CCO")
        self.assertTrue(
            solvent_species.is_isomorphic(solventlibrary.entries["ethanol"].item[0])
        )

        # Case 4: when the solventDatabase contains incorrect values for the molecule attribute, it raises Exception
        # This will display the SMILES Parse Error message from the external function, but ignore it.
        self.assertRaises(
            Exception,
            solventlibrary.load_entry,
            index=4,
            label="benzene",
            solvent=None,
            molecule="ring",
        )

        # Case 5: when the solventDatabase contains data for co-solvents.
        solventlibrary.load_entry(
            index=5, label="methanol_50_water_50", solvent=None, molecule=["CO", "O"]
        )
        solvent_species_list = [Species().from_smiles("CO"), Species().from_smiles("O")]
        self.assertEqual(len(solventlibrary.entries["methanol_50_water_50"].item), 2)
        for spc1 in solventlibrary.entries["methanol_50_water_50"].item:
            self.assertTrue(
                any([spc1.is_isomorphic(spc2) for spc2 in solvent_species_list])
            )


#####################################################


if __name__ == "__main__":
    suite = TestLoader().loadTestsFromTestCase(TestSoluteDatabase)
    TextTestRunner(verbosity=2).run(suite)
