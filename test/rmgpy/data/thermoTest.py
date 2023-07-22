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

import math
import os
import shutil
import unittest

import rmgpy
import rmgpy.constants as constants
from external.wip import work_in_progress
from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.thermo import (
    ThermoDatabase,
    ThermoData,
    ThermoCentralDatabaseInterface,
    convert_ring_to_sub_molecule,
    bicyclic_decomposition_for_polyring,
    combine_cycles,
    combine_two_rings_into_sub_molecule,
    find_aromatic_bonds_from_sub_molecule,
    get_copy_for_one_ring,
    is_aromatic_ring,
    is_bicyclic,
    is_ring_partial_matched,
    saturate_ring_bonds,
    split_bicyclic_into_single_rings,
)
from rmgpy.exceptions import DatabaseError
from rmgpy.ml.estimator import MLEstimator
from rmgpy.molecule.molecule import Molecule
from rmgpy.quantity import Quantity
from rmgpy.species import Species


################################################################################


def setUpModule():
    """A function that is run ONCE before all unit tests in this module."""
    global database
    database = RMGDatabase()
    database.load_thermo(
        os.path.join(settings["database.directory"], "thermo"),
        thermo_libraries=["DFT_QCI_thermo", "SABIC_aromatics", "primaryThermoLibrary"],
        surface=True,
    )
    database.load_surface(os.path.join(settings["database.directory"], "surface"))


def tearDownModule():
    """A function that is run ONCE after all unit tests in this module."""
    from rmgpy.data import rmg

    rmg.database = None


class TestThermoDatabaseLoading(unittest.TestCase):
    def test_failing_loads_thermo_libraries(self):
        database = ThermoDatabase()
        libraries = [
            "primaryThermoLibrary",
            "GRI-Mech3.0",
            "I am a library not existing in official RMG",
        ]
        path = os.path.join(settings["database.directory"], "thermo")

        with self.assertRaises(Exception):
            database.load_libraries(os.path.join(path, "libraries"), libraries)

    def test_loading_external_thermo_library(self):
        """This tests loading a thermo library which is not in the RMG-database repo"""
        thermo_lib_in_db_path = os.path.join(
            settings["database.directory"], "thermo", "libraries", "primaryNS.py"
        )
        thermo_lib_in_test_dir_path = os.path.join(
            os.path.dirname(rmgpy.__file__), "test_data", "copied_thermo_lib.py"
        )
        shutil.copyfile(src=thermo_lib_in_db_path, dst=thermo_lib_in_test_dir_path)
        database = ThermoDatabase()
        database.load_libraries(path="", libraries=[thermo_lib_in_test_dir_path])
        self.assertEqual(list(database.libraries.keys()), ["copied_thermo_lib"])
        os.remove(thermo_lib_in_test_dir_path)


class TestThermoDatabase(unittest.TestCase):
    """
    Contains unit tests of the ThermoDatabase class.
    """

    @classmethod
    def setUpClass(cls):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        cls.database = database.thermo
        cls.database.set_binding_energies("Pt111")

        cls.databaseWithoutLibraries = ThermoDatabase()
        cls.databaseWithoutLibraries.load(
            os.path.join(settings["database.directory"], "thermo"),
            libraries=[],
            surface=True,
        )
        cls.databaseWithoutLibraries.set_binding_energies("Pt111")

        # Set up ML estimator
        models_path = os.path.join(
            settings["database.directory"], "thermo", "ml", "main"
        )
        hf298_path = os.path.join(models_path, "hf298")
        s298_cp_path = os.path.join(models_path, "s298_cp")
        cls.ml_estimator = MLEstimator(hf298_path, s298_cp_path)

    def test_pickle(self):
        """
        Test that a ThermoDatabase object can be successfully pickled and
        unpickled with no loss of information.
        """
        import pickle

        thermodb0 = pickle.loads(pickle.dumps(self.database))

        self.assertEqual(thermodb0.library_order, self.database.library_order)
        self.assertEqual(
            sorted(thermodb0.depository.keys()), sorted(self.database.depository.keys())
        )

        self.assertEqual(
            sorted(thermodb0.libraries.keys()), sorted(self.database.libraries.keys())
        )
        self.assertEqual(
            sorted(thermodb0.groups.keys()), sorted(self.database.groups.keys())
        )

        for key, depository0 in thermodb0.depository.items():
            depository = self.database.depository[key]
            self.assertTrue(type(depository0), type(depository))
            self.assertEqual(
                sorted(depository0.entries.keys()), sorted(depository.entries.keys())
            )

        for key, library0 in thermodb0.libraries.items():
            library = self.database.libraries[key]
            self.assertTrue(type(library0), type(library))
            self.assertEqual(
                sorted(library0.entries.keys()), sorted(library.entries.keys())
            )

        for key, group0 in thermodb0.groups.items():
            group = self.database.groups[key]
            self.assertTrue(type(group0), type(group))
            self.assertEqual(
                sorted(group0.entries.keys()), sorted(group.entries.keys())
            )

    def test_symmetry_added_by_get_thermo_data(self):
        """
        Test that `get_thermo_data` properly accounts for symmetry in thermo
        by comping with the method `estimate_thermo_via_group_additivity`
        """

        spc = Species(molecule=[Molecule().from_smiles("C[CH]C=CC")])

        thermo_with_sym = self.databaseWithoutLibraries.get_thermo_data(spc)
        thermo_without_sym = (
            self.databaseWithoutLibraries.estimate_thermo_via_group_additivity(
                spc.molecule[0]
            )
        )

        symmetry_number = spc.get_symmetry_number()
        self.assertNotEqual(
            symmetry_number,
            spc.molecule[0].get_symmetry_number(),
            "For this test to be robust, species symmetry ({}) and molecule symmetry ({}) "
            "must be different".format(
                symmetry_number, spc.molecule[0].get_symmetry_number()
            ),
        )

        symmetry_contribution_to_entropy = -constants.R * math.log(symmetry_number)

        self.assertAlmostEqual(
            thermo_with_sym.get_entropy(298.0),
            thermo_without_sym.get_entropy(298.0) + symmetry_contribution_to_entropy,
            msg="The symmetry contribution is wrong {:.3f} /= {:.3f} + {:.3f}".format(
                thermo_with_sym.get_entropy(298.0),
                thermo_without_sym.get_entropy(298.0),
                symmetry_contribution_to_entropy,
            ),
        )

    def test_symmetry_contribution_radicals(self):
        """
        Test that the symmetry contribution is correctly added for radicals
        estimated via the HBI method.

        This is done by testing thermo_data from a database and from group
        additivity and ensuring they give the correct value.
        """
        spc = Species(molecule=[Molecule().from_smiles("[CH3]")])

        thermo_data_lib = self.database.get_thermo_data(spc)

        thermo_data_ga = self.databaseWithoutLibraries.get_thermo_data(spc)

        self.assertAlmostEqual(
            thermo_data_lib.get_entropy(298.0), thermo_data_ga.get_entropy(298.0), 0
        )

    def test_parse_thermo_comments(self):
        """
        Test that the ThermoDatabase.extract_source_from_comments function works properly
        on various thermo comments.
        """
        from rmgpy.thermo import NASA, NASAPolynomial

        # Pure group additivity thermo.
        gav_species = Species(
            index=3,
            label="c1c(O)c(O)c(CC(C)CC)cc1",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            -1.18833,
                            0.11272,
                            -4.26393e-05,
                            -2.12017e-08,
                            1.441e-11,
                            -51642.9,
                            38.8904,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(1078.35, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            26.6057,
                            0.0538434,
                            -2.22538e-05,
                            4.22393e-09,
                            -3.00808e-13,
                            -60208.4,
                            -109.218,
                        ],
                        Tmin=(1078.35, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""
Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) +
group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-O2s) + group(Cb-O2s) + group(Cb-H) +
group(Cb-H) + group(Cb-H) + group(O2s-CbH) + group(O2s-CbH) + longDistanceInteraction_cyclic(o_OH_OH) +
longDistanceInteraction_cyclic(o_OH_OH) + ring(Benzene)
""",
            ),
            molecule=[Molecule(smiles="c1c(O)c(O)c(CC(C)CC)cc1")],
        )

        source = self.database.extract_source_from_comments(gav_species)
        self.assertTrue(
            "GAV" in source, "Should have found that the thermo source is GAV."
        )
        self.assertEqual(len(source["GAV"]["group"]), 8)
        self.assertEqual(len(source["GAV"]["longDistanceInteraction_noncyclic"]), 1)
        self.assertEqual(len(source["GAV"]["longDistanceInteraction_cyclic"]), 1)
        self.assertEqual(len(source["GAV"]["ring"]), 1)

        # Pure library thermo
        dipk = Species(
            index=1,
            label="DIPK",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            3.35002,
                            0.017618,
                            -2.46235e-05,
                            1.7684e-08,
                            -4.87962e-12,
                            35555.7,
                            5.75335,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(888.28, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            6.36001,
                            0.00406378,
                            -1.73509e-06,
                            5.05949e-10,
                            -4.49975e-14,
                            35021,
                            -8.41155,
                        ],
                        Tmin=(888.28, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""Thermo library: DIPK""",
            ),
            molecule=[Molecule(smiles="CC(C)C(=O)C(C)C")],
        )

        source = self.database.extract_source_from_comments(dipk)
        self.assertTrue("Library" in source)

        # Mixed library and HBI thermo
        dipk_rad = Species(
            index=4,
            label="R_tert",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            2.90061,
                            0.0298018,
                            -7.06268e-05,
                            6.9636e-08,
                            -2.42414e-11,
                            54431,
                            5.44492,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(882.19, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            6.70999,
                            0.000201027,
                            6.65617e-07,
                            -7.99543e-11,
                            4.08721e-15,
                            54238.6,
                            -9.73662,
                        ],
                        Tmin=(882.19, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""Thermo library: DIPK + radical(C2CJCHO)""",
            ),
            molecule=[
                Molecule(smiles="C[C](C)C(=O)C(C)C"),
                Molecule(smiles="CC(C)=C([O])C(C)C"),
            ],
        )

        source = self.database.extract_source_from_comments(dipk_rad)
        self.assertTrue("Library" in source)
        self.assertTrue("GAV" in source)
        self.assertEqual(len(source["GAV"]["radical"]), 1)

        # Pure QM thermo
        cineole = Species(
            index=6,
            label="Cineole",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            -0.324129,
                            0.0619667,
                            9.71008e-05,
                            -1.60598e-07,
                            6.28285e-11,
                            -38699.9,
                            29.3686,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(985.52, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            20.6043,
                            0.0586913,
                            -2.22152e-05,
                            4.19949e-09,
                            -3.06046e-13,
                            -46791,
                            -91.4152,
                        ],
                        Tmin=(985.52, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""QM MopacMolPM3 calculation attempt 1""",
            ),
            molecule=[Molecule(smiles="CC12CCC(CC1)C(C)(C)O2")],
        )

        source = self.database.extract_source_from_comments(cineole)
        self.assertTrue("QM" in source)

        # Mixed QM and HBI thermo
        cineole_rad = Species(
            index=7,
            label="CineoleRad",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            -0.2897,
                            0.0627717,
                            8.63299e-05,
                            -1.47868e-07,
                            5.81665e-11,
                            -14017.6,
                            31.0266,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(988.76, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            20.4836,
                            0.0562555,
                            -2.13903e-05,
                            4.05725e-09,
                            -2.96023e-13,
                            -21915,
                            -88.1205,
                        ],
                        Tmin=(988.76, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_P)""",
            ),
            molecule=[Molecule(smiles="[CH2]C12CCC(CC1)C(C)(C)O2")],
        )

        source = self.database.extract_source_from_comments(cineole_rad)
        self.assertTrue("QM" in source)
        self.assertTrue("GAV" in source)
        self.assertEqual(len(source["GAV"]["radical"]), 1)

        # No thermo comments
        other = Species(
            index=7,
            label="CineoleRad",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            -0.2897,
                            0.0627717,
                            8.63299e-05,
                            -1.47868e-07,
                            5.81665e-11,
                            -14017.6,
                            31.0266,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(988.76, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            20.4836,
                            0.0562555,
                            -2.13903e-05,
                            4.05725e-09,
                            -2.96023e-13,
                            -21915,
                            -88.1205,
                        ],
                        Tmin=(988.76, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
            ),
            molecule=[Molecule(smiles="[CH2]C12CCC(CC1)C(C)(C)O2")],
        )

        # Function should complain if there's no thermo comments
        with self.assertRaises(ValueError):
            self.database.extract_source_from_comments(other)

        # Check a dummy species that has plus and minus thermo group contributions
        polycyclic = Species(
            index=7,
            label="dummy",
            thermo=NASA(
                polynomials=[
                    NASAPolynomial(
                        coeffs=[
                            -0.2897,
                            0.0627717,
                            8.63299e-05,
                            -1.47868e-07,
                            5.81665e-11,
                            -14017.6,
                            31.0266,
                        ],
                        Tmin=(100, "K"),
                        Tmax=(988.76, "K"),
                    ),
                    NASAPolynomial(
                        coeffs=[
                            20.4836,
                            0.0562555,
                            -2.13903e-05,
                            4.05725e-09,
                            -2.96023e-13,
                            -21915,
                            -88.1205,
                        ],
                        Tmin=(988.76, "K"),
                        Tmax=(5000, "K"),
                    ),
                ],
                Tmin=(100, "K"),
                Tmax=(5000, "K"),
                comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) - ring(Benzene)""",
            ),
            molecule=[Molecule(smiles="[CH2]C12CCC(CC1)C(C)(C)O2")],
        )

        source = self.database.extract_source_from_comments(polycyclic)
        self.assertTrue("GAV" in source)
        self.assertEqual(
            source["GAV"]["ring"][0][1], -1
        )  # the weight of benzene contribution should be -1
        self.assertEqual(
            source["GAV"]["group"][0][1], 2
        )  # weight of the group(Cs-CsCsHH) conbtribution should be 2

    def test_species_thermo_generation_hbi_library(self):
        """Test thermo generation for species objects for HBI correction on library value.

        Ensure that molecule list is only reordered, and not changed after matching library value
        """
        spec = Species().from_smiles("C[CH]c1ccccc1")
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.get_thermo_data(spec)

        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertTrue(
            "library" in thermo.comment,
            "Thermo not found from library, test purpose not fulfilled.",
        )

    def test_species_thermo_generation_hbi_gav(self):
        """Test thermo generation for species objects for HBI correction on group additivity value.

        Ensure that molecule list is only reordered, and not changed after group additivity
        """
        spec = Species().from_smiles("C[CH]c1ccccc1")
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.databaseWithoutLibraries.get_thermo_data(spec)

        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertTrue(
            "group additivity" in thermo.comment,
            "Thermo not found from GAV, test purpose not fulfilled.",
        )

    def test_species_thermo_generation_library(self):
        """Test thermo generation for species objects for library value.

        Ensure that the matched molecule is placed at the beginning of the list."""
        spec = Species().from_smiles("c12ccccc1c(C=[CH])ccc2")
        arom = Molecule().from_adjacency_list(
            """
multiplicity 2
1  C u0 p0 c0 {2,B} {3,B} {5,B}
2  C u0 p0 c0 {1,B} {4,B} {7,B}
3  C u0 p0 c0 {1,B} {6,B} {11,S}
4  C u0 p0 c0 {2,B} {8,B} {13,S}
5  C u0 p0 c0 {1,B} {9,B} {16,S}
6  C u0 p0 c0 {3,B} {10,B} {17,S}
7  C u0 p0 c0 {2,B} {10,B} {19,S}
8  C u0 p0 c0 {4,B} {9,B} {14,S}
9  C u0 p0 c0 {5,B} {8,B} {15,S}
10 C u0 p0 c0 {6,B} {7,B} {18,S}
11 C u0 p0 c0 {3,S} {12,D} {20,S}
12 C u1 p0 c0 {11,D} {21,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {12,S}
"""
        )
        spec.generate_resonance_structures()

        self.assertTrue(
            arom.is_isomorphic(spec.molecule[0])
        )  # The aromatic structure should be the first one
        # Move the aromatic structure to the end for testing
        spec.molecule.append(spec.molecule.pop(0))

        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.get_thermo_data(spec)

        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertTrue(
            arom.is_isomorphic(spec.molecule[0])
        )  # The aromatic structure should now be the first one
        self.assertTrue(
            "library" in thermo.comment,
            "Thermo not found from library, test purpose not fulfilled.",
        )

    def test_species_thermo_generation_ml(self):
        """Test thermo generation for species objects based on ML estimation."""

        # ML settings
        ml_settings = dict(
            min_heavy_atoms=1,
            max_heavy_atoms=None,
            min_carbon_atoms=0,
            max_carbon_atoms=None,
            min_oxygen_atoms=0,
            max_oxygen_atoms=None,
            min_nitrogen_atoms=0,
            max_nitrogen_atoms=None,
            only_cyclics=False,
            only_heterocyclics=False,
            min_cycle_overlap=0,
        )

        # Make these large so they don't influence estimation
        ml_uncertainty_cutoffs = dict(
            H298=Quantity(1e8, "kcal/mol"),
            S298=Quantity(1e8, "cal/(mol*K)"),
            Cp=Quantity(1e8, "cal/(mol*K)"),
        )
        ml_settings["uncertainty_cutoffs"] = ml_uncertainty_cutoffs

        spec1 = Species().from_smiles("C[CH]c1ccccc1")
        spec1.generate_resonance_structures()
        spec2 = Species().from_smiles("NC=O")

        thermo1 = self.database.get_thermo_data_from_ml(
            spec1, self.ml_estimator, ml_settings
        )
        thermo2 = self.database.get_thermo_data_from_ml(
            spec2, self.ml_estimator, ml_settings
        )
        self.assertIsInstance(thermo1, ThermoData)
        self.assertIsInstance(thermo2, ThermoData)
        self.assertTrue(
            "ML Estimation" in thermo1.comment,
            "Thermo not from ML estimation, test purpose not fulfilled",
        )
        self.assertTrue(
            "ML Estimation" in thermo2.comment,
            "Thermo not from ML estimation, test purpose not fulfilled",
        )

        # Now make these negative to make sure we don't use ML
        ml_uncertainty_cutoffs = dict(
            H298=Quantity(-1.0, "kcal/mol"),
            S298=Quantity(-1.0, "cal/(mol*K)"),
            Cp=Quantity(-1.0, "cal/(mol*K)"),
        )
        ml_settings["uncertainty_cutoffs"] = ml_uncertainty_cutoffs

        thermo1 = self.database.get_thermo_data_from_ml(
            spec1, self.ml_estimator, ml_settings
        )
        thermo2 = self.database.get_thermo_data_from_ml(
            spec2, self.ml_estimator, ml_settings
        )
        self.assertIsNone(thermo1)
        self.assertIsNone(thermo2)

    def test_thermo_generation_ml_settings(self):
        """Test that thermo generation with ML correctly respects settings"""

        # ML settings
        ml_settings = dict(
            min_heavy_atoms=5,
            max_heavy_atoms=6,
            min_carbon_atoms=5,
            max_carbon_atoms=5,
            min_oxygen_atoms=0,
            max_oxygen_atoms=None,
            min_nitrogen_atoms=0,
            max_nitrogen_atoms=None,
            only_cyclics=False,
            only_heterocyclics=False,
            min_cycle_overlap=0,
            uncertainty_cutoffs=dict(
                H298=Quantity(1e8, "kcal/mol"),
                S298=Quantity(1e8, "cal/(mol*K)"),
                Cp=Quantity(1e8, "cal/(mol*K)"),
            ),
        )

        spec1 = Species().from_smiles("CCCC")
        spec2 = Species().from_smiles("CCCCC")
        spec3 = Species().from_smiles("C1CC12CC2")
        spec4 = Species().from_smiles("C1CC2CC1O2")

        # Test atom limits
        thermo = self.database.get_thermo_data_from_ml(
            spec1, self.ml_estimator, ml_settings
        )
        self.assertIsNone(thermo)
        thermo = self.database.get_thermo_data_from_ml(
            spec2, self.ml_estimator, ml_settings
        )
        self.assertIsInstance(thermo, ThermoData)
        self.assertTrue(
            "ML Estimation" in thermo.comment,
            "Thermo not from ML estimation, test purpose not fulfilled",
        )

        # Test cyclic species
        ml_settings["only_cyclics"] = True
        thermo = self.database.get_thermo_data_from_ml(
            spec2, self.ml_estimator, ml_settings
        )
        self.assertIsNone(thermo)

        # Test cyclic species whether it outputs None when_both onlycyclics and heterocyclics are both True.
        ml_settings["only_heterocyclics"] = True
        thermo = self.database.get_thermo_data_from_ml(
            spec3, self.ml_estimator, ml_settings
        )
        self.assertIsNone(thermo)

        # Test heterocyclic species when both onlycyclics and heterocyclics are both True.
        thermo = self.database.get_thermo_data_from_ml(
            spec4, self.ml_estimator, ml_settings
        )
        self.assertIsInstance(thermo, ThermoData)
        self.assertTrue(
            "ML Estimation" in thermo.comment,
            "Thermo not from ML estimation, test purpose not fulfilled",
        )

        # Test cyclic species whether it outputs None when_only heterocyclics is True
        ml_settings["only_cyclics"] = False
        thermo = self.database.get_thermo_data_from_ml(
            spec3, self.ml_estimator, ml_settings
        )
        self.assertIsNone(thermo)

        # Test heterocyclic species when_only heterocyclics is True
        thermo = self.database.get_thermo_data_from_ml(
            spec4, self.ml_estimator, ml_settings
        )
        self.assertIsInstance(thermo, ThermoData)
        self.assertTrue(
            "ML Estimation" in thermo.comment,
            "Thermo not from ML estimation, test purpose not fulfilled",
        )

        # Test spiro species
        ml_settings["only_heterocyclics"] = False
        ml_settings["min_cycle_overlap"] = 1
        thermo = self.database.get_thermo_data_from_ml(
            spec3, self.ml_estimator, ml_settings
        )
        self.assertIsInstance(thermo, ThermoData)
        self.assertTrue(
            "ML Estimation" in thermo.comment,
            "Thermo not from ML estimation, test purpose not fulfilled",
        )

        # Test bridged species
        ml_settings["min_cycle_overlap"] = 3
        thermo = self.database.get_thermo_data_from_ml(
            spec3, self.ml_estimator, ml_settings
        )
        self.assertIsNone(thermo)
        thermo = self.database.get_thermo_data_from_ml(
            spec4, self.ml_estimator, ml_settings
        )
        self.assertIsInstance(thermo, ThermoData)
        self.assertTrue(
            "ML Estimation" in thermo.comment,
            "Thermo not from ML estimation, test purpose not fulfilled",
        )

    def test_thermo_estimation_not_affect_database(self):
        poly_root = self.database.groups["polycyclic"].entries["PolycyclicRing"]
        previous_enthalpy = poly_root.data.get_enthalpy(298) / 4184.0
        smiles = "C1C2CC1C=CC=C2"
        spec = Species().from_smiles(smiles)
        spec.generate_resonance_structures()

        thermo_gav = self.database.get_thermo_data_from_groups(spec)
        polycyclic_groups = self.database.get_ring_groups_from_comments(thermo_gav)[1]

        polycyclic_group_labels = [
            polycyclicGroup.label for polycyclicGroup in polycyclic_groups
        ]

        self.assertIn("PolycyclicRing", polycyclic_group_labels)

        latter_enthalpy = poly_root.data.get_enthalpy(298) / 4184.0

        self.assertAlmostEqual(previous_enthalpy, latter_enthalpy, 2)

    def test_get_all_thermo_data_fails_quietly(self):
        """Test that get_all_thermo_data doesn't break when GAV fails."""
        spec = Species().from_smiles("[Ne]")

        # Check that GAV fails
        with self.assertRaises(DatabaseError):
            self.database.get_thermo_data_from_groups(spec)

        # Check that get_all_thermo_data doesn't break
        thermo = self.database.get_all_thermo_data(spec)
        self.assertEqual(len(thermo), 1)

    def test_lowest_h298_for_resonance_structures(self):
        """Test that the thermo entry with the lowest H298 is selected for a species with resonance structures"""

        smiles = "[C]#C[O]"
        # has H298 ~= 640 kJ/mol; has resonance structure `[C]=C=O` with H298 ~= 380 kJ/mol
        spec = Species().from_smiles(smiles)
        thermo_gav1 = self.database.get_thermo_data_from_groups(spec)
        spec.generate_resonance_structures()
        thermo_gav2 = self.database.get_thermo_data_from_groups(spec)
        self.assertTrue(
            thermo_gav2.get_enthalpy(298) < thermo_gav1.get_enthalpy(298),
            msg="Did not select the molecule with the lowest H298 "
            "as a the thermo entry for [C]#C[O] / [C]=C=O",
        )

        smiles = "C=C[CH][O]"
        # has H298 ~= 209 kJ/mol; has (a reactive) resonance structure `C=CC=O` with H298 ~= -67 kJ/mol
        spec = Species().from_smiles(smiles)
        thermo_gav1 = self.database.get_thermo_data_from_groups(spec)
        spec.generate_resonance_structures()
        thermo_gav2 = self.database.get_thermo_data_from_groups(spec)
        self.assertTrue(
            thermo_gav2.get_enthalpy(298) < thermo_gav1.get_enthalpy(298),
            msg="Did not select the molecule with the lowest H298 "
            "as a the thermo entry for C=C[CH][O] / C=CC=O",
        )

        smiles = "CS(=O)#S(=O)C"
        # when using group additivity should instead select CS(=O)S(=O)C which has a lower enthalpy
        spec = Species().from_smiles(smiles)
        thermo_gav1 = self.database.get_thermo_data_from_groups(spec)
        spec.generate_resonance_structures()
        thermo_gav2 = self.database.get_thermo_data_from_groups(spec)
        self.assertTrue(
            thermo_gav2.get_enthalpy(298) < thermo_gav1.get_enthalpy(298),
            msg="Did not select the molecule with the lowest H298 "
            "as a the thermo entry for CS(=O)S(=O)C / CS(=O)#S(=O)C",
        )

    def test_thermo_for_mixed_reactive_and_nonreactive_molecules(self):
        """Test that the thermo entry of nonreactive molecules isn't selected for a species, even if it's more stable"""

        smiles = "[C]=C=O"  # has H298 ~= 640 kJ/mol; has resonance structure `[C]=C=O` with H298 ~= 380 kJ/mol
        spec = Species().from_smiles(smiles)
        thermo_gav1 = self.database.get_thermo_data_from_groups(
            spec
        )  # thermo of the stable molecule
        spec.generate_resonance_structures()
        spec.molecule[
            0
        ].reactive = False  # set the more stable molecule to nonreactive for this check
        thermo_gav2 = self.database.get_thermo_data_from_groups(
            spec
        )  # thermo of the speciesless stable molecule
        self.assertTrue(
            thermo_gav2.get_enthalpy(298) > thermo_gav1.get_enthalpy(298),
            msg="Did not select the reactive molecule for thermo",
        )

    def test_thermo_for_aromatic_radicals(self):
        """Test that we use the most aromatic resonance structure for thermo estimation"""
        spec = Species(smiles="C=[C]c1ccc2ccccc2c1")  # vinylnaphthalene radical
        spec.generate_resonance_structures()
        thermo_gav = self.database.get_thermo_data_from_groups(spec)
        self.assertAlmostEqual(thermo_gav.H298.value_si / 4184, 107, delta=1)
        self.assertIn(
            "group additivity",
            thermo_gav.comment,
            "Thermo not found from GAV, test purpose not fulfilled.",
        )
        self.assertIn("polycyclic(s2_6_6_naphthalene)", thermo_gav.comment)

    def test_identifying_missing_group(self):
        """Test identifying a missing GAV group"""
        # this test should be updated once data is added to the missing group
        spc = Species(smiles="S[N+]#[C-]")
        spc.generate_resonance_structures()
        thermo_gav = self.database.get_thermo_data_from_groups(spc)
        self.assertIn("missing(N5tc-C2tcS2s)", thermo_gav.comment)

    def test_adsorbate_thermo_generation_gav(self):
        """Test thermo generation for adsorbate from Group Additivity value.

        Ensure that molecule list is only reordered, and not changed after matching library value
        """
        spec = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 X u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}"""
                )
            ]
        )
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.databaseWithoutLibraries.get_thermo_data(spec)
        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertTrue(
            "group additivity estimation" in thermo.comment,
            "Thermo not found from group additivity, test purpose not fulfilled.",
        )
        self.assertTrue(
            "radical" in thermo.comment,
            "Didn't apply radical correction during group estimation.",
        )
        self.assertTrue(
            "Adsorption correction" in thermo.comment,
            "Adsorption correction not added to thermo.",
        )

    def test_adsorbate_thermo_generation_library(self):
        """Test thermo generation for adsorbate from gas phase library value.

        Ensure that molecule list is only reordered, and not changed after matching library value
        """
        spec = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 X u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}"""
                )
            ]
        )
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.get_thermo_data(spec)

        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertTrue(
            "library" in thermo.comment,
            "Thermo not found from library, test purpose not fulfilled.",
        )
        self.assertFalse(
            "radical" in thermo.comment,
            "Applied radical correction instead of finding CH3 in the library.",
        )
        self.assertTrue(
            "Adsorption correction" in thermo.comment,
            "Adsorption correction not added to thermo.",
        )

    def test_adsorbate_thermo_generation_bidentate(self):
        """Test thermo generation for a bidentate adsorbate, CH2XCH2X

        CH2-CH2
        |   |
        X   X
        """
        spec = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """
1 X u0 p0 c0 {3,S}
2 X u0 p0 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {4,S}"""
                )
            ]
        )
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.get_thermo_data(spec)

        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        # self.assertTrue('library' in thermo.comment, 'Thermo not found from library, test purpose not fulfilled.')
        self.assertFalse(
            "radical" in thermo.comment,
            "Applied radical correction instead of finding CH2=CH2",
        )
        self.assertTrue(
            "Adsorption correction" in thermo.comment,
            "Adsorption correction not added to thermo.",
        )

    def test_adsorbate_thermo_generation_bidentate_double(self):
        """Test thermo generation for a bidentate adsorbate, CH=XCH=X

        CH-CH
        ‖  ‖
        X  X
        """
        spec = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """
1 X u0 p0 c0 {3,D}
2 X u0 p0 c0 {4,D}
3 C u0 p0 c0 {1,D} {4,S} {5,S}
4 C u0 p0 c0 {2,D} {3,S} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}"""
                )
            ]
        )
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.get_thermo_data(spec)

        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertFalse(
            "radical" in thermo.comment,
            "Applied radical correction instead of finding CH#CH",
        )
        self.assertTrue(
            "Adsorption correction" in thermo.comment,
            "Adsorption correction not added to thermo.",
        )

    def test_adsorbate_thermo_generation_bidentate_C2(self):
        """Test thermo generation for a bidentate adsorbate [C]#[C]

        C#C
        | |
        X X
        """
        # Start with X-C#C-X
        spec = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """
1 X u0 p0 c0 {3,S}
2 X u0 p0 c0 {4,S}
3 C u0 p0 c0 {1,S} {4,T}
4 C u0 p0 c0 {2,S} {3,T}"""
                )
            ]
        )
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.get_thermo_data(spec)
        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertFalse(
            "radical" in thermo.comment,
            "Applied radical correction instead of finding C2(T) directly in library",
        )
        self.assertEqual(
            thermo.label,
            "C2(T)XX",
            "Should have found triplet C2 in the gas phase library",
        )
        self.assertTrue(
            "Adsorption correction" in thermo.comment,
            "Adsorption correction not added to thermo.",
        )
        # Now see what happens for X=C=C=X
        spec = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """
1 X u0 p0 c0 {3,D}
2 X u0 p0 c0 {4,D}
3 C u0 p0 c0 {1,D} {4,D}
4 C u0 p0 c0 {2,D} {3,D}"""
                )
            ]
        )
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.get_thermo_data(spec)
        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertFalse(
            "radical" in thermo.comment,
            "Applied radical correction instead of finding C2(T) directly in library",
        )
        self.assertEqual(
            thermo.label,
            "C2(T)XX",
            "Should have found triplet C2 in the gas phase library",
        )
        self.assertTrue(
            "Adsorption correction" in thermo.comment,
            "Adsorption correction not added to thermo.",
        )
        # Now see what happens for X#C-C#X
        spec = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """
1 X u0 p0 c0 {3,T}
2 X u0 p0 c0 {4,T}
3 C u0 p0 c0 {1,T} {4,S}
4 C u0 p0 c0 {2,T} {3,S}"""
                )
            ]
        )
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.get_thermo_data(spec)
        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertFalse(
            "radical" in thermo.comment,
            "Applied radical correction instead of finding C2(T) directly in library",
        )
        self.assertEqual(
            thermo.label,
            "C2(T)XX",
            "Should have found triplet C2 in the gas phase library",
        )
        self.assertTrue(
            "Adsorption correction" in thermo.comment,
            "Adsorption correction not added to thermo.",
        )

    def test_adsorbate_thermo_generation_bidentate_asymmetric(self):
        """Test thermo generation for a bidentate adsorbate, CH=XCH2X

        CH-CH2
        ‖  |
        X  X
        """
        spec = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """
1 X u0 p0 c0 {3,S}
2 X u0 p0 c0 {4,D}
3 C u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C u0 p0 c0 {2,D} {3,S} {7,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}"""
                )
            ]
        )
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.get_thermo_data(spec)
        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertFalse(
            "radical" in thermo.comment,
            "Applied radical correction instead of finding C2H3 in the DFT library",
        )
        self.assertTrue(
            "Adsorption correction" in thermo.comment,
            "Adsorption correction not added to thermo.",
        )

    def test_adsorbate_thermo_generation_bidentate_asymmetric_NNOH(self):
        """Test thermo generation for a bidentate adsorbate, N(=X)N(X)OH

        N--N--O-H
        ‖  |
        X  X
        """
        spec = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """
1 O u0 p2 c0 {2,S} {4,S}
2 N u0 p1 c0 {1,S} {3,S} {5,S}
3 N u0 p1 c0 {2,S} {6,D}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {2,S}
6 X u0 p0 c0 {3,D}"""
                )
            ]
        )
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.get_thermo_data(spec)
        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertTrue(
            "Adsorption correction" in thermo.comment,
            "Adsorption correction not added to thermo.",
        )
        match_str = "Gas phase thermo for ON=[N]"
        self.assertEqual(
            thermo.comment[0 : len(match_str)],
            match_str,
            "Gas phase species in thermo.comment should be ON=[N]",
        )

    def test_adsorbate_thermo_generation_bidentate_OO(self):
        """Test thermo generation for a bidentate adsorbate, [X]OO[X]

        O--O
        |  |
        X  X
        """
        spec = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {4,S}
3 X u0 p0 c0 {1,S}
4 X u0 p0 c0 {2,S}"""
                )
            ]
        )
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.get_thermo_data(spec)
        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertTrue(
            "Adsorption correction" in thermo.comment,
            "Adsorption correction not added to thermo.",
        )
        self.assertEqual(thermo.label, "O2XX", "thermo.label should be O2XX")

    def test_adsorbate_thermo_generation_bidentate_CO(self):
        """Test thermo generation for a bidentate adsorbate, [X][C-]=[O+][X]

        C- = O+
        |    |
        X    X
        """
        spec = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """
[Pt][C-]=[O+][Pt]
1 O u0 p1 c+1 {2,D} {4,S}
2 C u0 p1 c-1 {1,D} {3,S}
3 X u0 p0 c0 {2,S}
4 X u0 p0 c0 {1,S}"""
                )
            ]
        )
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.get_thermo_data(spec)
        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertTrue(
            "Adsorption correction" in thermo.comment,
            "Adsorption correction not added to thermo.",
        )
        self.assertEqual(thermo.label, "COXX", "thermo.label should be COXX")

    def test_adsorbate_thermo_raises_error(self):
        """Test thermo generation group tree error handling."""
        spec = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """
1 N u0 p1 c0 {2,D} {4,S}
2 N u0 p0 c+1 {1,D} {3,D}
3 N u0 p2 c-1 {2,D}
4 X u0 p0 c0 {1,S}
"""
                )
            ]
        )
        groups = self.database.groups["adsorptionPt111"].entries
        # save a few things we're about to change
        _nstarparent = groups["N*"].parent
        _nstardata = groups["N*"].data
        _ostardata = groups["O*"].data
        # change the database to cause errors
        groups["N*"].data = None
        groups["N*"].parent = None
        with self.assertRaisesRegex(
            DatabaseError, "Could not find an adsorption correction"
        ):
            thermo = self.database.get_thermo_data(spec)
        groups["N*"].data = "O*"
        groups["O*"].data = "N*"
        with self.assertRaisesRegex(DatabaseError, "circular reference"):
            thermo = self.database.get_thermo_data(spec)
        groups["N*"].data = "O*"
        groups["O*"].data = "foobar"
        with self.assertRaisesRegex(DatabaseError, "non-existing group"):
            thermo = self.database.get_thermo_data(spec)
        # Now restore the database to working order
        groups["N*"].parent = _nstarparent
        groups["N*"].data = _nstardata
        groups["O*"].data = _ostardata

    def test_adsorbate_thermo_generation_bidentate_weird_CO(self):
        """Test thermo generation for a bidentate adsorbate weird resonance of CO

        C-:O:
        #  |
        X  X
        """
        spec = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """
1 O u0 p2 c0 {2,S} {4,S}
2 C u0 p0 c0 {1,S} {3,T}
3 X u0 p0 c0 {2,T}
4 X u0 p0 c0 {1,S}"""
                )
            ]
        )
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.get_thermo_data(spec)
        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertTrue(
            "Adsorption correction" in thermo.comment,
            "Adsorption correction not added to thermo.",
        )

    def test_adsorbate_thermo_generation_bidentate_nonadjacent(self):
        """Test thermo generation for a bidentate adsorbate, CH2X-CH2-CH2X

        CH2-CH2-CH2
        |       |
        X       X
        """
        spec = Species(
            molecule=[
                Molecule().from_adjacency_list(
                    """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
3 C u0 p0 c0 {1,S} {8,S} {9,S} {11,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
10 X u0 p0 c0 {2,S}
11 X u0 p0 c0 {3,S}"""
                )
            ]
        )
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.get_thermo_data(spec)

        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        # self.assertTrue('library' in thermo.comment, 'Thermo not found from library, test purpose not fulfilled.')
        self.assertTrue(
            "radical" in thermo.comment,
            "Expected to need radical correctinos to find CH2j-CH2-CH2j",
        )
        self.assertTrue(
            "Adsorption correction" in thermo.comment,
            "Adsorption correction not added to thermo.",
        )


class TestThermoAccuracy(unittest.TestCase):
    """
    Contains tests for accuracy of thermo estimates and symmetry calculations.
    """

    @classmethod
    def setUpClass(cls):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        cls.database = database.thermo

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.Tlist = [300, 400, 500, 600, 800, 1000, 1500]

        self.testCases = [
            # SMILES         symm    H298   S298  Cp300  Cp400  Cp500  Cp600  Cp800 Cp1000 Cp1500
            # 1,3-hexadiene decomposition products
            [
                "C=CC=CCC",
                3,
                13.45,
                86.37,
                29.49,
                37.67,
                44.54,
                50.12,
                58.66,
                64.95,
                74.71,
            ],
            [
                "[CH]=CC=CCC",
                3,
                72.55,
                87.76,
                29.30,
                36.92,
                43.18,
                48.20,
                55.84,
                61.46,
                70.18,
            ],
            [
                "C=[C]C=CCC",
                3,
                61.15,
                87.08,
                29.68,
                36.91,
                43.03,
                48.11,
                55.96,
                61.78,
                71.54,
            ],
            [
                "C=C[C]=CCC",
                3,
                61.15,
                87.08,
                29.68,
                36.91,
                43.03,
                48.11,
                55.96,
                61.78,
                71.54,
            ],
            [
                "C=CC=[C]CC",
                3,
                70.35,
                88.18,
                29.15,
                36.46,
                42.6,
                47.6,
                55.32,
                61.04,
                69.95,
            ],
            [
                "C=CC=C[CH]C",
                3,
                38.24,
                84.41,
                27.79,
                35.46,
                41.94,
                47.43,
                55.74,
                61.92,
                71.86,
            ],
            [
                "C=CC=CC[CH2]",
                2,
                62.45,
                89.78,
                28.72,
                36.31,
                42.63,
                47.72,
                55.50,
                61.21,
                70.05,
            ],
            ["[CH3]", 6, 34.81, 46.37, 9.14, 10.18, 10.81, 11.34, 12.57, 13.71, 15.2],
            [
                "C=CC=C[CH2]",
                2,
                46.11,
                75.82,
                22.54,
                28.95,
                34.24,
                38.64,
                45.14,
                49.97,
                57.85,
            ],
            ["[CH2]C", 6, 28.6, 59.87, 11.73, 14.47, 17.05, 19.34, 23.02, 25.91, 31.53],
            [
                "C=CC=[CH]",
                1,
                85.18,
                69.37,
                18.93,
                23.55,
                27.16,
                29.92,
                34.02,
                37.03,
                41.81,
            ],
            [
                "C=[CH]",
                1,
                71.62,
                56.61,
                10.01,
                11.97,
                13.66,
                15.08,
                17.32,
                19.05,
                21.85,
            ],
            [
                "[CH]=CCC",
                3,
                58.99,
                75.0,
                20.38,
                25.34,
                29.68,
                33.36,
                39.14,
                43.48,
                50.22,
            ],
            # Cyclic Structures
            [
                "C1CCCCC1",
                12,
                -29.45,
                69.71,
                27.20,
                37.60,
                46.60,
                54.80,
                67.50,
                76.20,
                88.50,
            ],
            ["C1CCC1", 8, 6.51, 63.35, 17.39, 23.91, 29.86, 34.76, 42.40, 47.98, 56.33],
            [
                "C1C=CC=C1",
                2,
                32.5,
                65.5,
                18.16,
                24.71,
                30.25,
                34.7,
                41.25,
                45.83,
                52.61,
            ],
        ]

    @work_in_progress
    def test_new_thermo_generation(self):
        """
        Test that the new ThermoDatabase generates appropriate thermo data.
        """
        for (
            smiles,
            symm,
            H298,
            S298,
            Cp300,
            Cp400,
            Cp500,
            Cp600,
            Cp800,
            Cp1000,
            Cp1500,
        ) in self.testCases:
            cp_list = [Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500]
            species = Species().from_smiles(smiles)
            species.generate_resonance_structures()
            thermo_data = self.database.get_thermo_data_from_groups(species)
            molecule = species.molecule[0]
            for mol in species.molecule[1:]:
                thermo_data0 = self.database.get_all_thermo_data(
                    Species(molecule=[mol])
                )[0][0]
                for data in self.database.get_all_thermo_data(Species(molecule=[mol]))[
                    1:
                ]:
                    if data[0].get_enthalpy(298) < thermo_data0.get_enthalpy(298):
                        thermo_data0 = data[0]
                if thermo_data0.get_enthalpy(298) < thermo_data.get_enthalpy(298):
                    thermo_data = thermo_data0
                    molecule = mol
            self.assertAlmostEqual(
                H298,
                thermo_data.get_enthalpy(298) / 4184,
                places=1,
                msg="H298 error for {0}. Expected {1}, but calculated {2}.".format(
                    smiles, H298, thermo_data.get_enthalpy(298) / 4184
                ),
            )
            self.assertAlmostEqual(
                S298,
                thermo_data.get_entropy(298) / 4.184,
                places=1,
                msg="S298 error for {0}. Expected {1}, but calculated {2}.".format(
                    smiles, S298, thermo_data.get_entropy(298) / 4.184
                ),
            )
            for T, Cp in zip(self.Tlist, cp_list):
                self.assertAlmostEqual(
                    Cp,
                    thermo_data.get_heat_capacity(T) / 4.184,
                    places=1,
                    msg="Cp{3} error for {0}. Expected {1} but calculated {2}.".format(
                        smiles, Cp, thermo_data.get_heat_capacity(T) / 4.184, T
                    ),
                )

    def test_symmetry_number_generation(self):
        """
        Test we generate symmetry numbers correctly.

        This uses the new thermo database to generate the H298, used
        to select the stablest resonance isomer.
        """
        for (
            smiles,
            symm,
            H298,
            S298,
            Cp300,
            Cp400,
            Cp500,
            Cp600,
            Cp800,
            Cp1000,
            Cp1500,
        ) in self.testCases:
            species = Species().from_smiles(smiles)
            calc_symm = species.get_symmetry_number()
            self.assertEqual(
                symm,
                calc_symm,
                msg="Symmetry number error for {0}. Expected {1} but calculated {2}.".format(
                    smiles, symm, calc_symm
                ),
            )


class TestThermoAccuracyAromatics(TestThermoAccuracy):
    """
    Contains tests for accuracy of thermo estimates and symmetry calculations for aromatics only.

    A copy of the above class, but with different test compounds.
    """

    def setUp(self):
        self.Tlist = [300, 400, 500, 600, 800, 1000, 1500]
        self.testCases = [
            # SMILES         symm    H298   S298  Cp300  Cp400  Cp500  Cp600  Cp800 Cp1000 Cp1500
            [
                "c1ccccc1",
                12,
                19.80,
                64.24,
                19.44,
                26.64,
                32.76,
                37.80,
                45.24,
                50.46,
                58.38,
            ],
            [
                "c1ccc2ccccc2c1",
                4,
                36.0,
                79.49,
                31.94,
                42.88,
                52.08,
                59.62,
                70.72,
                78.68,
                90.24,
            ],
        ]

    def __init__(self, *args, **kwargs):
        super(TestThermoAccuracyAromatics, self).__init__(*args, **kwargs)
        self._testMethodDoc = (
            self._testMethodDoc.strip().split("\n")[0] + " for Aromatics.\n"
        )

    def test_long_distance_interaction_in_aromatic_molecule(self):
        """
        Test long distance interaction is properly calculated for aromatic molecule.
        """
        spec = Species().from_smiles("c(O)1c(O)c(C=O)c(C=O)c(O)c(C=O)1")
        spec.generate_resonance_structures()
        thermo = self.database.get_thermo_data_from_groups(spec)

        self.assertIn("o_OH_OH", thermo.comment)
        self.assertIn("o_OH_CHO", thermo.comment)
        self.assertIn("o_CHO_CHO", thermo.comment)
        self.assertIn("m_CHO_CHO", thermo.comment)
        self.assertIn("p_OH_OH", thermo.comment)
        self.assertIn("p_OH_CHO", thermo.comment)
        self.assertIn("p_CHO_CHO", thermo.comment)

    def test_long_distance_interaction_in_aromatic_radical(self):
        """
        Test long distance interaction is properly calculated for aromatic radical.
        """
        spec = Species().from_smiles("c([O])1c(C=O)c(C=O)c(OC)cc1")
        spec.generate_resonance_structures()
        thermo = self.database.get_thermo_data_from_groups(spec)

        self.assertNotIn("o_OH_CHO", thermo.comment)
        self.assertNotIn("p_OH_MeO", thermo.comment)
        self.assertIn("o_Oj_CHO", thermo.comment)
        self.assertIn("m_Oj_CHO", thermo.comment)
        self.assertIn("p_Oj_OCH3", thermo.comment)
        self.assertIn("o_CHO_CHO", thermo.comment)
        self.assertIn("o_CHO_MeO", thermo.comment)

    def test_long_distance_interaction_in_aromatic_biradical(self):
        """
        Test long distance interaction is properly calculated for aromatic biradical.
        """
        spec = Species().from_smiles("c([O])1c([C]=O)cc(C=O)cc1")
        spec.generate_resonance_structures()
        thermo = self.database.get_thermo_data_from_groups(spec)

        self.assertNotIn("o_OH_CHO", thermo.comment)
        self.assertNotIn("m_CHO_CHO", thermo.comment)
        self.assertNotIn("p_OH_CHO", thermo.comment)
        self.assertNotIn("o_Oj_CHO", thermo.comment)
        self.assertIn("m_Cj=O_CHO", thermo.comment)


class TestCyclicThermo(unittest.TestCase):
    """
    Contains unit tests of the ThermoDatabase class.
    """

    @classmethod
    def setUpClass(cls):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        cls.database = database.thermo

    def test_compute_group_additivity_thermo_for_two_ring_molecule(self):
        """
        The molecule being tested has two rings, one is 13cyclohexadiene5methylene
        the other is benzene ring. This method is to test thermo estimation will
        give two different corrections accordingly.
        """
        spec = Species().from_smiles("CCCCCCCCCCCC(CC=C1C=CC=CC1)c1ccccc1")
        spec.generate_resonance_structures()
        thermo = self.database.get_thermo_data_from_groups(spec)

        ring_groups, polycyclic_groups = self.database.get_ring_groups_from_comments(
            thermo
        )
        self.assertEqual(len(ring_groups), 2)
        self.assertEqual(len(polycyclic_groups), 0)

        expected_matched_rings_labels = ["13cyclohexadiene5methylene", "Benzene"]
        expected_matched_rings = [
            self.database.groups["ring"].entries[label]
            for label in expected_matched_rings_labels
        ]

        self.assertEqual(set(ring_groups), set(expected_matched_rings))

    def test_thermo_for_monocyclic_and_polycyclic_same_molecule(self):
        """
        Test a molecule that has both a polycyclic and a monocyclic ring in the same molecule
        """
        spec = Species().from_smiles("C(CCC1C2CCC1CC2)CC1CCC1")
        spec.generate_resonance_structures()
        thermo = self.database.get_thermo_data_from_groups(spec)
        ring_groups, polycyclic_groups = self.database.get_ring_groups_from_comments(
            thermo
        )
        self.assertEqual(len(ring_groups), 1)
        self.assertEqual(len(polycyclic_groups), 1)

        expected_matched_rings_labels = ["Cyclobutane"]
        expected_matched_rings = [
            self.database.groups["ring"].entries[label]
            for label in expected_matched_rings_labels
        ]
        self.assertEqual(set(ring_groups), set(expected_matched_rings))

        expected_matched_polyrings_labels = ["s3_5_5_ane"]
        expected_matched_polyrings = [
            self.database.groups["polycyclic"].entries[label]
            for label in expected_matched_polyrings_labels
        ]

        self.assertEqual(set(polycyclic_groups), set(expected_matched_polyrings))

    def test_get_ring_groups_from_comments(self):
        """
        Test that get_ring_groups_from_comments method works for fused polycyclics.
        """
        from rmgpy.thermo.thermoengine import generate_thermo_data

        smi = "C12C(C3CCC2C3)C4CCC1C4"  # two norbornane rings fused together
        spc = Species().from_smiles(smi)

        spc.thermo = generate_thermo_data(spc)

        self.database.get_ring_groups_from_comments(spc.thermo)

    def test_remove_group(self):
        """
        Test that removing groups using nodes near the root of radical.py
        """
        # load up test data designed for this test
        database2 = ThermoDatabase()
        path = os.path.join(os.path.dirname(rmgpy.__file__), "data/test_data/")
        database2.load(os.path.join(path, "thermo"), depository=False)

        # load up the thermo radical database as a test
        rad_group = database2.groups["radical"]

        # use root as removed groups parent, which should be an LogicOr node
        root = rad_group.top[0]
        # use group to remove as
        group_to_remove = rad_group.entries["RJ"]
        children = group_to_remove.children

        # remove the group
        rad_group.remove_group(group_to_remove)

        # afterwards group_to_remove should not be in the database or root's children
        self.assertNotIn(group_to_remove, list(rad_group.entries.values()))
        self.assertNotIn(group_to_remove, root.children)

        for child in children:
            # group_to_remove children should all be in roots item.component and children attribuetes
            self.assertIn(child.label, root.item.components)
            self.assertIn(child, root.children)
            # the children should all have root a their parent now
            self.assertIs(child.parent, root)

        # Specific to ThermoDatabase, (above test apply to all base class Database)
        # we check that unicode entry.data pointers are correctly reassigned

        # if group_to_remove is a pointer and another node pointed to it, we copy
        # group_to_remove pointer
        self.assertIs(rad_group.entries["OJ"].data, group_to_remove.data)

        # Remove an entry with a ThermoData object
        group_to_remove2 = rad_group.entries["CsJ"]
        rad_group.remove_group(group_to_remove2)
        # If group_to_remove was a data object, we point toward parent instead
        self.assertEqual(
            rad_group.entries["RJ2_triplet"].data, group_to_remove2.parent.label
        )
        # If the parent pointed toward group_to_remove, we need should have copied data object
        Tlist = [300, 400, 500, 600, 800, 1000, 1500]
        self.assertNotIsInstance(group_to_remove2.parent.data, str)
        self.assertEqual(
            group_to_remove2.parent.data.get_enthalpy(298),
            group_to_remove2.data.get_enthalpy(298),
        )
        self.assertEqual(
            group_to_remove2.parent.data.get_entropy(298),
            group_to_remove2.data.get_entropy(298),
        )
        self.assertTrue(
            all(
                [
                    group_to_remove2.parent.data.get_heat_capacity(x)
                    == group_to_remove2.data.get_heat_capacity(x)
                    for x in Tlist
                ]
            )
        )

    def test_is_ring_partial_matched(self):
        # create testing molecule
        smiles = "C1CC2CCCC3CCCC(C1)C23"
        mol = Molecule().from_smiles(smiles)
        polyring = [atom for atom in mol.atoms if atom.is_non_hydrogen()]

        # create matched group
        matched_group = (
            self.database.groups["polycyclic"].entries["PolycyclicRing"].item
        )

        # test
        self.assertTrue(is_ring_partial_matched(polyring, matched_group))

    def test_add_ring_correction_thermo_data_from_tree_for_existing_tricyclic(self):
        # create testing molecule: C1CC2C3CCC(C3)C2C1
        # this tricyclic molecule is already in polycyclic database
        # so algorithm should give complete match: s2-3_5_5_5_ane
        smiles = "C1CC2C3CCC(C3)C2C1"
        mol = Molecule().from_smiles(smiles)
        polyring = mol.get_disparate_cycles()[1][0]

        poly_groups = self.database.groups["polycyclic"]
        _, matched_entry, _ = self.database._add_ring_correction_thermo_data_from_tree(
            None, poly_groups, mol, polyring
        )

        self.assertEqual(matched_entry.label, "s2-3_5_5_5_ane")

    def test_add_poly_ring_correction_thermo_data_from_heuristic_using_pyrene(self):
        # create testing molecule: Pyrene with two ring of aromatic version
        # the other two ring of kekulized version
        #
        # creating it seems not natural in RMG, that's because
        # RMG cannot parse the adjacencyList of that isomer correctly
        # so here we start with pyrene radical and get the two aromatic ring isomer
        # then saturate it.
        smiles = "C1C=C2C=CC=C3C=CC4=CC=CC=1C4=C23"
        spe = Species().from_smiles(smiles)
        spe.generate_resonance_structures()
        mols = []
        for mol in spe.molecule:
            sssr0 = mol.get_smallest_set_of_smallest_rings()
            aromatic_ring_num = 0
            for sr0 in sssr0:
                sr0mol = Molecule(atoms=sr0)
                if is_aromatic_ring(sr0mol):
                    aromatic_ring_num += 1
            if aromatic_ring_num == 2:
                mols.append(mol)

        ring_group_labels = []
        polycyclic_group_labels = []
        for mol in mols:
            polyring = mol.get_disparate_cycles()[1][0]

            thermo_data = ThermoData(
                Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
                Cpdata=([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "J/(mol*K)"),
                H298=(0.0, "kJ/mol"),
                S298=(0.0, "J/(mol*K)"),
            )

            self.database._add_poly_ring_correction_thermo_data_from_heuristic(
                thermo_data, polyring
            )

            (
                ring_groups,
                polycyclic_groups,
            ) = self.database.get_ring_groups_from_comments(thermo_data)

            ring_group_labels += [ringGroup.label for ringGroup in ring_groups]
            polycyclic_group_labels += [
                polycyclicGroup.label for polycyclicGroup in polycyclic_groups
            ]

        self.assertIn("Benzene", ring_group_labels)
        self.assertIn("Cyclohexene", ring_group_labels)
        self.assertIn("s2_6_6_ben_ene_1", polycyclic_group_labels)
        self.assertIn("s2_6_6_diene_2_7", polycyclic_group_labels)

    def test_add_poly_ring_correction_thermo_data_from_heuristic_using_aromatic_tricyclic(
        self,
    ):
        # create testing molecule
        #
        # creating it seems not natural in RMG, that's because
        # RMG cannot parse the adjacencyList of that isomer correctly
        # so here we start with kekulized version and generate_resonance_structures
        # and pick the one with two aromatic rings
        smiles = "C1=CC2C=CC=C3C=CC(=C1)C=23"
        spe = Species().from_smiles(smiles)
        spe.generate_resonance_structures()
        for mol in spe.molecule:
            sssr0 = mol.get_smallest_set_of_smallest_rings()
            aromatic_ring_num = 0
            for sr0 in sssr0:
                sr0mol = Molecule(atoms=sr0)
                if is_aromatic_ring(sr0mol):
                    aromatic_ring_num += 1
            if aromatic_ring_num == 2:
                break

        # extract polyring from the molecule
        polyring = mol.get_disparate_cycles()[1][0]

        thermo_data = ThermoData(
            Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
            Cpdata=([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "J/(mol*K)"),
            H298=(0.0, "kJ/mol"),
            S298=(0.0, "J/(mol*K)"),
        )

        self.database._add_poly_ring_correction_thermo_data_from_heuristic(
            thermo_data, polyring
        )

        ring_groups, polycyclic_groups = self.database.get_ring_groups_from_comments(
            thermo_data
        )

        ring_group_labels = [ringGroup.label for ringGroup in ring_groups]
        polycyclic_group_labels = [
            polycyclicGroup.label for polycyclicGroup in polycyclic_groups
        ]

        self.assertIn("Benzene", ring_group_labels)
        self.assertIn("Cyclopentene", ring_group_labels)
        self.assertIn("s2_5_6_indene", polycyclic_group_labels)
        self.assertIn("s2_6_6_naphthalene", polycyclic_group_labels)

    def test_add_poly_ring_correction_thermo_data_from_heuristic_using_alkane_tricyclic(
        self,
    ):
        # create testing molecule
        smiles = "C1CC2CCCC3C(C1)C23"
        mol = Molecule().from_smiles(smiles)

        # extract polyring from the molecule
        polyring = mol.get_disparate_cycles()[1][0]

        thermo_data = ThermoData(
            Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
            Cpdata=([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "J/(mol*K)"),
            H298=(0.0, "kJ/mol"),
            S298=(0.0, "J/(mol*K)"),
        )

        self.database._add_poly_ring_correction_thermo_data_from_heuristic(
            thermo_data, polyring
        )

        ring_groups, polycyclic_groups = self.database.get_ring_groups_from_comments(
            thermo_data
        )

        ring_group_labels = [ringGroup.label for ringGroup in ring_groups]
        polycyclic_group_labels = [
            polycyclicGroup.label for polycyclicGroup in polycyclic_groups
        ]

        self.assertIn("Cyclohexane", ring_group_labels)
        self.assertIn("Cyclopropane", ring_group_labels)
        self.assertIn("s2_6_6_ane", polycyclic_group_labels)
        self.assertIn("s2_3_6_ane", polycyclic_group_labels)

    def test_add_poly_ring_correction_thermo_data_from_heuristic_using_highly_unsaturated_polycyclics1(
        self,
    ):
        """
        Test proper thermo estimation for highly unsaturated polycyclic whose decomposed
        bicyclics are not stored in database. Those bicyclics thermo will be estimated through
        a heuristic formula.

        In the future, the test assertion may be updated if some of the decomposed bicyclics
        have been added to database.
        """
        # create testing molecule
        smiles = "[CH]=C1C2=C=C3C=CC1C=C32"
        mol = Molecule().from_smiles(smiles)

        # extract polyring from the molecule
        polyring = mol.get_disparate_cycles()[1][0]

        thermo_data = ThermoData(
            Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
            Cpdata=([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "J/(mol*K)"),
            H298=(0.0, "kJ/mol"),
            S298=(0.0, "J/(mol*K)"),
        )

        self.database._add_poly_ring_correction_thermo_data_from_heuristic(
            thermo_data, polyring
        )

        ring_groups, polycyclic_groups = self.database.get_ring_groups_from_comments(
            thermo_data
        )

        ring_group_labels = [ringGroup.label for ringGroup in ring_groups]
        polycyclic_group_labels = [
            polycyclicGroup.label for polycyclicGroup in polycyclic_groups
        ]

        self.assertIn("1,4-Cyclohexadiene", ring_group_labels)
        self.assertIn("Cyclopentene", ring_group_labels)
        self.assertIn("cyclobutadiene_13", ring_group_labels)
        self.assertIn("s3_5_6_ane", polycyclic_group_labels)
        self.assertIn("s2_4_6_ane", polycyclic_group_labels)
        self.assertIn("s2_4_5_ane", polycyclic_group_labels)

    def test_add_poly_ring_correction_thermo_data_from_heuristic_using_highly_unsaturated_polycyclics2(
        self,
    ):
        """
        Test proper thermo estimation for highly unsaturated polycyclic whose decomposed
        bicyclics are not stored in database. Those bicyclics thermo will be estimated through
        a heuristic formula.

        In the future, the test assertion may be updated if some of the decomposed bicyclics
        have been added to database.
        """
        # create testing molecule
        smiles = "C1=C2C#CC3C=CC1C=C23"
        mol = Molecule().from_smiles(smiles)

        # extract polyring from the molecule
        polyring = mol.get_disparate_cycles()[1][0]

        thermo_data = ThermoData(
            Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
            Cpdata=([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "J/(mol*K)"),
            H298=(0.0, "kJ/mol"),
            S298=(0.0, "J/(mol*K)"),
        )

        self.database._add_poly_ring_correction_thermo_data_from_heuristic(
            thermo_data, polyring
        )

        ring_groups, polycyclic_groups = self.database.get_ring_groups_from_comments(
            thermo_data
        )

        ring_group_labels = [ringGroup.label for ringGroup in ring_groups]
        polycyclic_group_labels = [
            polycyclicGroup.label for polycyclicGroup in polycyclic_groups
        ]

        self.assertIn("1,4-Cyclohexadiene", ring_group_labels)
        self.assertIn("Cyclopentyne", ring_group_labels)
        self.assertIn("Cyclopentadiene", ring_group_labels)
        self.assertIn("s3_5_6_ane", polycyclic_group_labels)
        self.assertIn("s2_5_6_ane", polycyclic_group_labels)
        self.assertIn("s2_5_5_ane", polycyclic_group_labels)

    def test_get_bicyclic_correction_thermo_data_from_heuristic1(self):
        """
        Test bicyclic correction estimated properly from heuristic formula
        The test molecule "C1=CCC2C1=C2" has a shared atom with Cd atomtype,
        but in the correction estimation we stil expect the five-member ring
        part to match Cyclopentene
        """
        smiles = "C1=CCC2C1=C2"
        mol = Molecule().from_smiles(smiles)

        # extract polyring from the molecule
        polyring = mol.get_disparate_cycles()[1][0]

        thermo_data = self.database.get_bicyclic_correction_thermo_data_from_heuristic(
            polyring
        )

        ring_groups, polycyclic_groups = self.database.get_ring_groups_from_comments(
            thermo_data
        )

        ring_group_labels = [ringGroup.label for ringGroup in ring_groups]
        polycyclic_group_labels = [
            polycyclicGroup.label for polycyclicGroup in polycyclic_groups
        ]

        self.assertIn("Cyclopentane", ring_group_labels)
        self.assertIn("Cyclopropane", ring_group_labels)
        self.assertIn("Cyclopentene", ring_group_labels)
        self.assertIn("Cyclopropene", ring_group_labels)
        self.assertIn("s2_3_5_ane", polycyclic_group_labels)

    def test_get_bicyclic_correction_thermo_data_from_heuristic2(self):
        """
        Test bicyclic correction estimated properly from heuristic formula
        The test molecule "C1=CCC2=C1C2" doesn't have controversial shared
        atomtypes in correction estimation, which is regarded as a simple case.
        """
        smiles = "C1=CCC2=C1C2"
        mol = Molecule().from_smiles(smiles)

        # extract polyring from the molecule
        polyring = mol.get_disparate_cycles()[1][0]

        thermo_data = self.database.get_bicyclic_correction_thermo_data_from_heuristic(
            polyring
        )

        ring_groups, polycyclic_groups = self.database.get_ring_groups_from_comments(
            thermo_data
        )

        ring_group_labels = [ringGroup.label for ringGroup in ring_groups]
        polycyclic_group_labels = [
            polycyclicGroup.label for polycyclicGroup in polycyclic_groups
        ]

        self.assertIn("Cyclopentane", ring_group_labels)
        self.assertIn("Cyclopropane", ring_group_labels)
        self.assertIn("Cyclopentadiene", ring_group_labels)
        self.assertIn("Cyclopropene", ring_group_labels)
        self.assertIn("s2_3_5_ane", polycyclic_group_labels)


class TestMolecularManipulationInvolvedInThermoEstimation(unittest.TestCase):
    """
    Contains unit tests for methods of molecular manipulations for thermo estimation
    """

    def test_convert_ring_to_sub_molecule(self):
        # list out testing moleculess
        smiles1 = "C1CCCCC1"
        smiles2 = "C1CCC2CCCCC2C1"
        smiles3 = "C1CC2CCCC3CCCC(C1)C23"
        mol1 = Molecule().from_smiles(smiles1)
        mol2 = Molecule().from_smiles(smiles2)
        mol3 = Molecule().from_smiles(smiles3)

        # get ring structure by only extracting non-hydrogens
        ring1 = [atom for atom in mol1.atoms if atom.is_non_hydrogen()]
        ring2 = [atom for atom in mol2.atoms if atom.is_non_hydrogen()]
        ring3 = [atom for atom in mol3.atoms if atom.is_non_hydrogen()]

        # convert to submolecules
        submol1, _ = convert_ring_to_sub_molecule(ring1)
        submol2, _ = convert_ring_to_sub_molecule(ring2)
        submol3, _ = convert_ring_to_sub_molecule(ring3)

        # test against expected submolecules
        self.assertEqual(len(submol1.atoms), 6)
        self.assertEqual(len(submol2.atoms), 10)
        self.assertEqual(len(submol3.atoms), 13)

        bonds1 = []
        for atom in submol1.atoms:
            for bond in atom.edges.values():
                if bond not in bonds1:
                    bonds1.append(bond)

        bonds2 = []
        for atom in submol2.atoms:
            for bond in atom.edges.values():
                if bond not in bonds2:
                    bonds2.append(bond)

        bonds3 = []
        for atom in submol3.atoms:
            for bond in atom.edges.values():
                if bond not in bonds3:
                    bonds3.append(bond)

        self.assertEqual(len(bonds1), 6)
        self.assertEqual(len(bonds2), 11)
        self.assertEqual(len(bonds3), 15)

    def test_get_copy_for_one_ring(self):
        """
        This method tests the get_copy_for_one_ring method, which returns
        an atom object list that contains deep copies of the atoms
        """

        test_atom_list = Molecule(smiles="C1CCCCC1").atoms
        copied_atom_list = get_copy_for_one_ring(test_atom_list)

        test_molecule = Molecule(atoms=test_atom_list)
        copied_molecule = Molecule(atoms=copied_atom_list)

        self.assertNotEqual(test_atom_list, copied_atom_list)
        self.assertEqual(len(test_atom_list), len(copied_atom_list))
        self.assertEqual(test_molecule, copied_molecule)

    def test_to_fail_combine_two_rings_into_sub_molecule(self):
        """
        Test that if two non-overlapped rings lead to AssertionError
        """

        smiles1 = "C1CCCCC1"
        smiles2 = "C1CCCCC1"
        mol1 = Molecule().from_smiles(smiles1)
        mol2 = Molecule().from_smiles(smiles2)

        ring1 = [atom for atom in mol1.atoms if atom.is_non_hydrogen()]
        ring2 = [atom for atom in mol2.atoms if atom.is_non_hydrogen()]
        with self.assertRaises(AssertionError):
            combine_two_rings_into_sub_molecule(ring1, ring2)

    def test_combine_two_rings_into_sub_molecule(self):
        # create testing molecule
        smiles1 = "C1CCC2CCCCC2C1"
        mol1 = Molecule().from_smiles(smiles1)

        # get two SSSRs
        sssr = mol1.get_smallest_set_of_smallest_rings()
        ring1 = sssr[0]
        ring2 = sssr[1]

        # combine two rings into submolecule
        submol, _ = combine_two_rings_into_sub_molecule(ring1, ring2)

        self.assertEqual(len(submol.atoms), 10)

        bonds = []
        for atom in submol.atoms:
            for bondAtom, bond in atom.edges.items():
                if bond not in bonds:
                    bonds.append(bond)

        self.assertEqual(len(bonds), 11)

    def test_is_aromatic_ring(self):
        # create testing rings
        smiles1 = "C1CCC1"
        smiles2 = "C1CCCCC1"
        adj3 = """1  C u0 p0 c0 {2,B} {6,B} {7,S}
2  C u0 p0 c0 {1,B} {3,B} {8,S}
3  C u0 p0 c0 {2,B} {4,B} {9,S}
4  C u0 p0 c0 {3,B} {5,B} {10,S}
5  C u0 p0 c0 {4,B} {6,B} {11,S}
6  C u0 p0 c0 {1,B} {5,B} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
        """
        mol1 = Molecule().from_smiles(smiles1)
        mol2 = Molecule().from_smiles(smiles2)
        mol3 = Molecule().from_adjacency_list(adj3)
        ring1mol = Molecule(
            atoms=[atom for atom in mol1.atoms if atom.is_non_hydrogen()]
        )
        ring2mol = Molecule(
            atoms=[atom for atom in mol2.atoms if atom.is_non_hydrogen()]
        )
        ring3mol = Molecule(
            atoms=[atom for atom in mol3.atoms if atom.is_non_hydrogen()]
        )

        # check with expected results
        self.assertEqual(is_aromatic_ring(ring1mol), False)
        self.assertEqual(is_aromatic_ring(ring2mol), False)
        self.assertEqual(is_aromatic_ring(ring3mol), True)

    def test_is_bicyclic1(self):
        """
        Test is_bicyclic identifies bicyclic correctly
        The test molecule is bicyclic, we expect is_bicyclic()
        returns True.
        """
        smiles = "C1=CCC2C1=C2"
        mol = Molecule().from_smiles(smiles)
        polyring = mol.get_disparate_cycles()[1][0]

        self.assertTrue(is_bicyclic(polyring))

    def test_is_bicyclic2(self):
        """
        Test is_bicyclic identifies bicyclic correctly
        The test molecule is tetracyclic, we expect
        is_bicyclic() returns False
        """
        smiles = "C1C=C2C=CC=C3C=CC4=CC=CC=1C4=C23"
        mol = Molecule().from_smiles(smiles)
        polyring = mol.get_disparate_cycles()[1][0]

        self.assertFalse(is_bicyclic(polyring))

    def test_find_aromatic_bonds_from_sub_molecule(self):
        smiles = "C1=CC=C2C=CC=CC2=C1"
        spe = Species().from_smiles(smiles)
        spe.generate_resonance_structures()
        mol = spe.molecule[0]

        # get two SSSRs
        sssr = mol.get_smallest_set_of_smallest_rings()
        ring1 = sssr[0]
        ring2 = sssr[1]

        # create two testing submols
        submol1 = Molecule(atoms=ring1)
        submol2 = Molecule(atoms=ring2)

        # check with expected results
        self.assertEqual(len(find_aromatic_bonds_from_sub_molecule(submol1)), 6)
        self.assertEqual(len(find_aromatic_bonds_from_sub_molecule(submol2)), 6)

    def test_bicyclic_decomposition_for_polyring_using_pyrene(self):
        # create testing molecule: Pyrene with two ring of aromatic version
        # the other two ring of kekulized version
        #
        # creating it seems not natural in RMG, that's because
        # RMG cannot parse the adjacencyList of that isomer correctly
        # so here we start with pyrene radical and get the two aromatic ring isomer
        # then saturate it.
        smiles = "C1C=C2C=CC=C3C=CC4=CC=CC=1C4=C23"
        spe = Species().from_smiles(smiles)
        spe.generate_resonance_structures()
        for mol in spe.molecule:
            sssr0 = mol.get_smallest_set_of_smallest_rings()
            aromatic_ring_num = 0
            for sr0 in sssr0:
                sr0mol = Molecule(atoms=sr0)
                if is_aromatic_ring(sr0mol):
                    aromatic_ring_num += 1
            if aromatic_ring_num == 2:
                break

        # extract polyring from the molecule
        polyring = mol.get_disparate_cycles()[1][0]

        bicyclic_list, ring_occurrences_dict = bicyclic_decomposition_for_polyring(
            polyring
        )
        for bicyclic in bicyclic_list:
            bicyclic.delete_hydrogens()

        # 1st test: number of cores
        self.assertEqual(len(bicyclic_list), 5)

        # 2nd test: ring_occurrences_dict
        ring_in_core_occurrences = sorted(ring_occurrences_dict.values())
        expected_ring_in_core_occurances = [2, 2, 3, 3]
        self.assertEqual(ring_in_core_occurrences, expected_ring_in_core_occurances)

        # 3rd test: size of each bicyclic core
        bicyclic_sizes = sorted([len(bicyclic.atoms) for bicyclic in bicyclic_list])
        expected_bicyclic_sizes = [10, 10, 10, 10, 10]
        self.assertEqual(bicyclic_sizes, expected_bicyclic_sizes)

        # 4th test: bond info for members of each core
        aromatic_bond_num_in_bicyclics = []
        for bicyclic in bicyclic_list:
            aromatic_bond_num = len(find_aromatic_bonds_from_sub_molecule(bicyclic))
            aromatic_bond_num_in_bicyclics.append(aromatic_bond_num)
        aromatic_bond_num_in_bicyclics = sorted(aromatic_bond_num_in_bicyclics)
        expected_aromatic_bond_num_in_bicyclics = [0, 6, 6, 6, 6]
        self.assertEqual(
            aromatic_bond_num_in_bicyclics, expected_aromatic_bond_num_in_bicyclics
        )

    def test_bicyclic_decomposition_for_polyring_using_aromatic_tricyclic(self):
        # create testing molecule
        #
        # creating it seems not natural in RMG, that's because
        # RMG cannot parse the adjacencyList of that isomer correctly
        # so here we start with kekulized version and generate_resonance_structures
        # and pick the one with two aromatic rings
        smiles = "C1=CC2C=CC=C3C=CC(=C1)C=23"
        spe = Species().from_smiles(smiles)
        spe.generate_resonance_structures()
        for mol in spe.molecule:
            sssr0 = mol.get_smallest_set_of_smallest_rings()
            aromatic_ring_num = 0
            for sr0 in sssr0:
                sr0mol = Molecule(atoms=sr0)
                if is_aromatic_ring(sr0mol):
                    aromatic_ring_num += 1
            if aromatic_ring_num == 2:
                break

        # extract polyring from the molecule
        polyring = mol.get_disparate_cycles()[1][0]

        bicyclic_list, ring_occurrences_dict = bicyclic_decomposition_for_polyring(
            polyring
        )
        for bicyclic in bicyclic_list:
            bicyclic.delete_hydrogens()

        # 1st test: number of cores
        self.assertEqual(len(bicyclic_list), 3)

        # 2nd test: ring_occurrences_dict
        ring_in_core_occurrences = sorted(ring_occurrences_dict.values())
        expected_ring_in_core_occurrences = [2, 2, 2]
        self.assertEqual(ring_in_core_occurrences, expected_ring_in_core_occurrences)

        # 3rd test: size of each bicyclic core
        bicyclic_sizes = sorted([len(bicyclic.atoms) for bicyclic in bicyclic_list])
        expected_bicyclic_sizes = [9, 9, 10]
        self.assertEqual(bicyclic_sizes, expected_bicyclic_sizes)

        # 4th test: bond info for members of each core
        aromatic_bond_num_in_bicyclics = []
        for bicyclic in bicyclic_list:
            aromatic_bond_num = len(find_aromatic_bonds_from_sub_molecule(bicyclic))
            aromatic_bond_num_in_bicyclics.append(aromatic_bond_num)
        aromatic_bond_num_in_bicyclics = sorted(aromatic_bond_num_in_bicyclics)
        expected_aromatic_bond_num_in_bicyclics = [6, 6, 11]
        self.assertEqual(
            aromatic_bond_num_in_bicyclics, expected_aromatic_bond_num_in_bicyclics
        )

    def test_bicyclic_decomposition_for_polyring_using_alkane_tricyclic(self):
        # create testing molecule
        smiles = "C1CC2CCCC3C(C1)C23"
        mol = Molecule().from_smiles(smiles)

        # extract polyring from the molecule
        polyring = mol.get_disparate_cycles()[1][0]

        bicyclic_list, ring_occurrences_dict = bicyclic_decomposition_for_polyring(
            polyring
        )
        for bicyclic in bicyclic_list:
            bicyclic.delete_hydrogens()

        # 1st test: number of cores
        self.assertEqual(len(bicyclic_list), 3)

        # 2nd test: ring_occurrences_dict
        ring_in_core_occurrences = sorted(ring_occurrences_dict.values())
        expected_ring_in_core_occurrences = [2, 2, 2]
        self.assertEqual(ring_in_core_occurrences, expected_ring_in_core_occurrences)

        # 3rd test: size of each bicyclic core
        bicyclic_sizes = sorted([len(bicyclic.atoms) for bicyclic in bicyclic_list])
        expected_bicyclic_sizes = [7, 7, 10]
        self.assertEqual(bicyclic_sizes, expected_bicyclic_sizes)

        # 4th test: bond info for members of each core
        aromatic_bond_num_in_bicyclics = []
        for bicyclic in bicyclic_list:
            aromatic_bond_num = len(find_aromatic_bonds_from_sub_molecule(bicyclic))
            aromatic_bond_num_in_bicyclics.append(aromatic_bond_num)
        aromatic_bond_num_in_bicyclics = sorted(aromatic_bond_num_in_bicyclics)
        expected_aromatic_bond_num_in_bicyclics = [0, 0, 0]
        self.assertEqual(
            aromatic_bond_num_in_bicyclics, expected_aromatic_bond_num_in_bicyclics
        )

    def test_combine_cycles(self):
        """
        This method tests the combine_cycles method, which simply joins two lists
        together without duplication.
        """
        main_cycle = Molecule(smiles="C1CCC2CCCCC2C1").atoms
        test_cycle1 = main_cycle[0:8]
        test_cycle2 = main_cycle[6:]
        joined_cycle = combine_cycles(test_cycle1, test_cycle2)
        self.assertEqual(set(main_cycle), set(joined_cycle))

    def test_split_bicyclic_into_single_rings1(self):
        """
        Test bicyclic molecule "C1=CCC2C1=C2" can be divided into
        individual rings properly
        """
        smiles = "C1=CCC2C1=C2"
        mol = Molecule().from_smiles(smiles)
        bicyclic = mol.get_disparate_cycles()[1][0]

        bicyclic_submol = convert_ring_to_sub_molecule(bicyclic)[0]
        single_ring_submols = split_bicyclic_into_single_rings(bicyclic_submol)
        self.assertEqual(len(single_ring_submols), 2)

        single_ring_submol_a, single_ring_submol_b = sorted(
            single_ring_submols, key=lambda submol: len(submol.atoms)
        )

        single_ring_submol_a.saturate_unfilled_valence()
        single_ring_submol_b.saturate_unfilled_valence()

        expected_submol_a = Molecule().from_smiles("C1=CC1")
        expected_submol_a.update_connectivity_values()

        expected_submol_b = Molecule().from_smiles("C1=CCCC1")
        expected_submol_b.update_connectivity_values()

        self.assertTrue(single_ring_submol_a.is_isomorphic(expected_submol_a))
        self.assertTrue(single_ring_submol_b.is_isomorphic(expected_submol_b))

    def test_split_bicyclic_into_single_rings2(self):
        """
        Test bicyclic molecule "C1=CCC2=C1C2" can be divided into
        individual rings properly
        """

        smiles = "C1=CCC2=C1C2"
        mol = Molecule().from_smiles(smiles)
        bicyclic = mol.get_disparate_cycles()[1][0]

        bicyclic_submol = convert_ring_to_sub_molecule(bicyclic)[0]
        single_ring_submols = split_bicyclic_into_single_rings(bicyclic_submol)
        self.assertEqual(len(single_ring_submols), 2)

        single_ring_submol_a, single_ring_submol_b = sorted(
            single_ring_submols, key=lambda submol: len(submol.atoms)
        )

        single_ring_submol_a.saturate_unfilled_valence()
        single_ring_submol_b.saturate_unfilled_valence()

        expected_submol_a = Molecule().from_smiles("C1=CC1")
        expected_submol_a.update_connectivity_values()

        expected_submol_b = Molecule().from_smiles("C1=CC=CC1")
        expected_submol_b.update_connectivity_values()

        self.assertTrue(single_ring_submol_a.is_isomorphic(expected_submol_a))
        self.assertTrue(single_ring_submol_b.is_isomorphic(expected_submol_b))

    def test_saturate_ring_bonds1(self):
        """
        Test unsaturated bonds of "C1=CCC2=C1C2" to be saturated properly
        """
        smiles = "C1=CCC2=C1C2"
        mol = Molecule().from_smiles(smiles)
        ring_submol = convert_ring_to_sub_molecule(mol.get_disparate_cycles()[1][0])[0]

        saturated_ring_submol, already_saturated = saturate_ring_bonds(ring_submol)

        expected_saturated_ring_submol = Molecule().from_smiles("C1CCC2C1C2")

        expected_saturated_ring_submol.update_connectivity_values()

        self.assertFalse(already_saturated)
        self.assertEqual(
            saturated_ring_submol.multiplicity,
            expected_saturated_ring_submol.multiplicity,
        )
        self.assertTrue(
            saturated_ring_submol.is_isomorphic(expected_saturated_ring_submol)
        )

    def test_saturate_ring_bonds2(self):
        """
        Test unsaturated bonds of "C1=CC=C2CCCCC2=C1" to be saturated properly
        """
        smiles = "C1=CC=C2CCCCC2=C1"
        spe = Species().from_smiles(smiles)
        spe.generate_resonance_structures()
        mol = spe.molecule[0]
        ring_submol = convert_ring_to_sub_molecule(mol.get_disparate_cycles()[1][0])[0]

        saturated_ring_submol, already_saturated = saturate_ring_bonds(ring_submol)

        expected_spe = Species().from_smiles("C1=CC=C2CCCCC2=C1")
        expected_spe.generate_resonance_structures()
        expected_saturated_ring_submol = expected_spe.molecule[0]

        expected_saturated_ring_submol.update_connectivity_values()

        self.assertTrue(already_saturated)
        self.assertEqual(
            saturated_ring_submol.multiplicity,
            expected_saturated_ring_submol.multiplicity,
        )
        self.assertTrue(
            saturated_ring_submol.is_isomorphic(expected_saturated_ring_submol)
        )

    def test_saturate_ring_bonds3(self):
        """
        Test unsaturated bonds of "C1=CC=C2CC=CCC2=C1" to be saturated properly
        """
        smiles = "C1=CC=C2CC=CCC2=C1"
        spe = Species().from_smiles(smiles)
        spe.generate_resonance_structures()
        mol = spe.molecule[0]
        ring_submol = convert_ring_to_sub_molecule(mol.get_disparate_cycles()[1][0])[0]

        saturated_ring_submol, already_saturated = saturate_ring_bonds(ring_submol)

        expected_spe = Species().from_smiles("C1=CC=C2CCCCC2=C1")
        expected_spe.generate_resonance_structures()
        expected_saturated_ring_submol = expected_spe.molecule[0]

        expected_saturated_ring_submol.update_connectivity_values()

        self.assertFalse(already_saturated)
        self.assertEqual(
            saturated_ring_submol.multiplicity,
            expected_saturated_ring_submol.multiplicity,
        )
        self.assertTrue(
            saturated_ring_submol.is_isomorphic(expected_saturated_ring_submol)
        )


def get_testing_tcd_authentication_info():
    try:
        host = os.environ["TCD_HOST"]
        port = int(os.environ["TCD_PORT"])
        username = os.environ["TCD_USER"]
        password = os.environ["TCD_PW"]
    except KeyError:
        print(
            "Thermo Central Database Authentication Environment Variables Not Completely Set!"
        )
        return None, 0, None, None

    return host, port, username, password


def is_tcd_available():
    """Check if TCD is available."""
    import platform
    import subprocess

    host = get_testing_tcd_authentication_info()[0]
    if host is not None:
        arg = "-n" if platform.system() == "Windows" else "-c"
        result = subprocess.call(["ping", arg, "1", host]) == 0
    else:
        result = False

    return result


@unittest.skip("Skipping TCD unit tests because database is offline.")
class TestThermoCentralDatabaseInterface(unittest.TestCase):
    """
    Contains unit tests for methods of ThermoCentralDatabaseInterface
    """

    @classmethod
    def setUpClass(cls):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        cls.database = database.thermo

    def connect_to_test_central_database(self):
        host, port, username, password = get_testing_tcd_authentication_info()
        application = "test"

        tcdi = ThermoCentralDatabaseInterface(
            host, port, username, password, application
        )
        return tcdi

    def test_connect_failure(self):
        host = "somehost"
        port = 27017
        username = "me"
        password = "pswd"
        application = "test"

        tcdi = ThermoCentralDatabaseInterface(
            host, port, username, password, application
        )

        self.assertTrue(tcdi.client is None)

    def test_connect_success(self):
        tcdi = self.connect_to_test_central_database()

        self.assertTrue(tcdi.client is not None)

    def test_satisfy_registration_requirements1(self):
        """
        the species is non-cyclic, currently regarded no need to
        register in thermo central database
        """
        tcdi = self.connect_to_test_central_database()

        species = Species().from_smiles("C[CH2]")

        thermo_data = self.database.get_thermo_data_from_groups(species)

        self.assertFalse(
            tcdi.satisfy_registration_requirements(species, thermo_data, self.database)
        )

    def test_satisfy_registration_requirements2(self):
        """
        the species is for non-cyclic, so no need to register in
        thermo central database
        """
        tcdi = self.connect_to_test_central_database()

        species = Species().from_smiles("CC")

        thermo_data = self.database.get_thermo_data_from_groups(species)

        self.assertFalse(
            tcdi.satisfy_registration_requirements(species, thermo_data, self.database)
        )

    def test_satisfy_registration_requirements3(self):
        """
        the thermo is exact match, so no need to register in
        thermo central database
        """
        tcdi = self.connect_to_test_central_database()

        species = Species().from_smiles("C1CC1")

        thermo_data = self.database.get_thermo_data_from_groups(species)

        self.assertFalse(
            tcdi.satisfy_registration_requirements(species, thermo_data, self.database)
        )

    def test_satisfy_registration_requirements4(self):
        """
        the thermo is from library, so no need to register in
        thermo central database
        """
        tcdi = self.connect_to_test_central_database()

        species = Species().from_smiles("[H][H]")

        thermo_data = self.database.get_thermo_data_from_libraries(species)

        self.assertFalse(
            tcdi.satisfy_registration_requirements(species, thermo_data, self.database)
        )

    def test_satisfy_registration_requirements5(self):
        """
        the thermo is matching generic node, so it needs to register in
        thermo central database

        In the future, if RMG-database includes corresponding exact match
        this test should be modified.
        """
        tcdi = self.connect_to_test_central_database()

        species = Species().from_smiles("C1C=CC2C=CC2=C1")

        thermo_data = self.database.get_thermo_data_from_groups(species)

        self.assertTrue(
            tcdi.satisfy_registration_requirements(species, thermo_data, self.database)
        )

    def test_satisfy_registration_requirements6(self):
        """
        the thermo is matching generic node, so it needs to register in
        thermo central database

        In the future, if RMG-database includes corresponding exact match
        this test should be modified.
        """
        tcdi = self.connect_to_test_central_database()

        species = Species().from_smiles("C1=C=C2CC23C=CC=1C=C3")

        thermo_data = self.database.get_thermo_data_from_groups(species)

        self.assertTrue(
            tcdi.satisfy_registration_requirements(species, thermo_data, self.database)
        )

    def test_register_in_central_thermo_db1(self):
        """
        Test situation where both registration_table and results_table have no
        species as the one going to be registered
        """
        # connect to thermo central database
        host, port, username, password = get_testing_tcd_authentication_info()
        application = "test"
        tcdi = ThermoCentralDatabaseInterface(
            host, port, username, password, application
        )

        # prepare species to register
        species = Species().from_smiles("C1=C=C2CC23C=CC=1C=C3")
        expected_aug_inchi = "InChI=1S/C10H6/c1-2-9-7-10(9)5-3-8(1)4-6-10/h3-6H,7H2"

        # select registration table
        # and clean previous data
        db = getattr(tcdi.client, "thermoCentralDB")
        registration_table = getattr(db, "registration_table")
        results_table = getattr(db, "results_table")
        registration_table.delete_many({"aug_inchi": expected_aug_inchi})
        results_table.delete_many({"aug_inchi": expected_aug_inchi})

        tcdi.register_in_central_thermo_db(species)
        registered_species_entries = list(
            registration_table.find({"aug_inchi": expected_aug_inchi})
        )

        # should expect only one registered such species
        self.assertEqual(len(registered_species_entries), 1)
        registered_species_entry = registered_species_entries[0]

        # check all the columns are expected
        registered_species = Species().from_smiles(
            str(registered_species_entry["SMILES_input"])
        )
        self.assertEqual(registered_species_entry["aug_inchi"], expected_aug_inchi)
        self.assertTrue(registered_species.is_isomorphic(species))
        self.assertIn(registered_species_entry["status"], ["pending", "submitted"])

        # clean up the table
        registration_table.delete_many({"aug_inchi": expected_aug_inchi})

    def test_register_in_central_thermo_db2(self):
        """
        Test situation where registration_table has species as the one going
        to be registered
        """

        # connect to thermo central database
        host, port, username, password = get_testing_tcd_authentication_info()
        application = "test"

        tcdi = ThermoCentralDatabaseInterface(
            host, port, username, password, application
        )

        # prepare species to register
        species = Species().from_smiles("C1=C=C2CC23C=CC=1C=C3")
        expected_aug_inchi = "InChI=1S/C10H6/c1-2-9-7-10(9)5-3-8(1)4-6-10/h3-6H,7H2"

        # select registration table
        # and clean previous data
        db = getattr(tcdi.client, "thermoCentralDB")
        registration_table = getattr(db, "registration_table")
        results_table = getattr(db, "results_table")
        registration_table.delete_many({"aug_inchi": expected_aug_inchi})
        registration_table.insert_one({"aug_inchi": expected_aug_inchi})
        results_table.delete_many({"aug_inchi": expected_aug_inchi})

        tcdi.register_in_central_thermo_db(species)
        registered_species_entries = list(
            registration_table.find({"aug_inchi": expected_aug_inchi})
        )

        # should expect only one registered such species
        self.assertEqual(len(registered_species_entries), 1)
        registered_species_entry = registered_species_entries[0]

        # check all the columns are expected
        self.assertEqual(registered_species_entry["aug_inchi"], expected_aug_inchi)
        self.assertTrue(len(registered_species_entry), 2)

        # clean up the table
        registration_table.delete_many({"aug_inchi": expected_aug_inchi})

    def test_register_in_central_thermo_db3(self):
        """
        Test situation where results_table has species as the one going
        to be registered
        """

        # connect to thermo central database
        host, port, username, password = get_testing_tcd_authentication_info()
        application = "test"

        tcdi = ThermoCentralDatabaseInterface(
            host, port, username, password, application
        )

        # prepare species to register
        species = Species().from_smiles("C1=C=C2CC23C=CC=1C=C3")
        expected_aug_inchi = "InChI=1S/C10H6/c1-2-9-7-10(9)5-3-8(1)4-6-10/h3-6H,7H2"

        # select registration table
        # and clean previous data
        db = getattr(tcdi.client, "thermoCentralDB")
        registration_table = getattr(db, "registration_table")
        results_table = getattr(db, "results_table")
        registration_table.delete_many({"aug_inchi": expected_aug_inchi})
        results_table.delete_many({"aug_inchi": expected_aug_inchi})
        results_table.insert_one({"aug_inchi": expected_aug_inchi})

        tcdi.register_in_central_thermo_db(species)
        registered_species_entries = list(
            registration_table.find({"aug_inchi": expected_aug_inchi})
        )

        # should expect only one registered such species
        self.assertEqual(len(registered_species_entries), 0)

        # clean up the table
        results_table.delete_many({"aug_inchi": expected_aug_inchi})


################################################################################
