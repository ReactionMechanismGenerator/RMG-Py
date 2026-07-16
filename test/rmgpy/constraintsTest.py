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
This script contains unit tests of the :mod:`rmgpy.constraints` module.
"""


from unittest import mock

from rmgpy.rmg.main import RMG
from rmgpy.constraints import (
    fails_species_constraints,
    is_polymer_constraint_member,
    _evaluate_constraints,
    reset_polymer_warning,
)
from rmgpy.species import Species
from rmgpy.molecule import Molecule
import rmgpy.rmg.input


class TestFailsSpeciesConstraints:
    """
    Contains unit tests of the fails_species_constraints function.
    """

    @classmethod
    def setup_class(cls):
        """
        A function run ONCE before all unit tests in this class.
        """
        cls.rmg = RMG()
        rmgpy.rmg.input.rmg = cls.rmg
        rmgpy.rmg.input.generated_species_constraints(
            maximumCarbonAtoms=2,
            maximumOxygenAtoms=1,
            maximumNitrogenAtoms=1,
            maximumSiliconAtoms=1,
            maximumSulfurAtoms=1,
            maximumSurfaceSites=2,
            maximumSurfaceBondOrder=3,
            maximumHeavyAtoms=3,
            maximumRadicalElectrons=2,
            maximumSingletCarbenes=1,
            maximumCarbeneRadicals=0,
        )

    @classmethod
    def teardown_class(cls):
        """
        A function run ONCE after all unit tests in this class.
        """
        rmgpy.rmg.input.rmg = None

    @mock.patch("rmgpy.constraints.logging")
    def test_constraints_not_loaded(self, mock_logging):
        """
        Test what happens when constraints are not loaded.
        """
        # Reset module level rmg variable in rmgpy.rmg.input
        rmgpy.rmg.input.rmg = None

        mol = Molecule(smiles="C")
        assert not fails_species_constraints(mol)

        mock_logging.debug.assert_called_with("Species constraints could not be found.")

        # Restore module level rmg variable in rmgpy.rmg.input
        rmgpy.rmg.input.rmg = self.rmg

    def test_species_input(self):
        """
        Test that fails_species_constraints can handle a Species object.
        """
        spc = Species().from_smiles("C")

        assert not fails_species_constraints(spc)

    def test_explicitly_allowed_molecules(self):
        """
        Test that we can explicitly allow molecules in species constraints.
        """
        mol = Molecule(smiles="CCCC")
        assert fails_species_constraints(mol)

        self.rmg.species_constraints["explicitlyAllowedMolecules"] = [Molecule(smiles="CCCC")]
        assert not fails_species_constraints(mol)

    def test_carbon_constraint(self):
        """
        Test that we can constrain the max number of carbon atoms.
        """
        mol1 = Molecule(smiles="CC")
        assert not fails_species_constraints(mol1)

        mol2 = Molecule(smiles="CCC")
        assert fails_species_constraints(mol2)

    def test_oxygen_constraint(self):
        """
        Test that we can constrain the max number of oxygen atoms.
        """
        mol1 = Molecule(smiles="C=O")
        assert not fails_species_constraints(mol1)

        mol2 = Molecule(smiles="OC=O")
        assert fails_species_constraints(mol2)

    def test_nitrogen_constraint(self):
        """
        Test that we can constrain the max number of nitrogen atoms.
        """
        mol1 = Molecule(smiles="CN")
        assert not fails_species_constraints(mol1)

        mol2 = Molecule(smiles="NCN")
        assert fails_species_constraints(mol2)

    def test_silicon_constraint(self):
        """
        Test that we can constrain the max number of silicon atoms.
        """
        mol1 = Molecule(smiles="[SiH4]")
        assert not fails_species_constraints(mol1)

        mol2 = Molecule(smiles="[SiH3][SiH3]")
        assert fails_species_constraints(mol2)

    def test_sulfur_constraint(self):
        """
        Test that we can constrain the max number of sulfur atoms.
        """
        mol1 = Molecule(smiles="CS")
        assert not fails_species_constraints(mol1)

        mol2 = Molecule(smiles="SCS")
        assert fails_species_constraints(mol2)

    def test_surface_site_constraint(self):
        """
        Test that we can constrain the max number of surface sites.
        """

        mol_1site = Molecule().from_adjacency_list(
            """
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,D}
3 X u0 p0 c0 {2,D}
"""
        )
        mol_2site = Molecule().from_adjacency_list(
            """
1 C u0 p0 c0 {2,D} {3,D}
2 C u0 p0 c0 {1,D} {4,D}
3 X u0 p0 c0 {1,D}
4 X u0 p0 c0 {2,D}
"""
        )

        mol_3site_vdW = Molecule().from_adjacency_list(
            """
1 C u0 p0 c0 {2,D} {3,D}
2 C u0 p0 c0 {1,D} {4,D}
3 X u0 p0 c0 {1,D}
4 X u0 p0 c0 {2,D}
6 X u0 p0 c0
"""
        )

        mol_3site = Molecule().from_adjacency_list(
            """
1 C u0 p0 c0 {4,S} {2,D} {7,S}
2 C u0 p0 c0 {1,D} {3,S} {8,S}
3 C u0 p0 c0 {2,S} {5,S} {6,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 X u0 p0 c0 {1,S}
8 X u0 p0 c0 {2,S}
9 X u0 p0 c0 {3,S}
"""
        )
        max_carbon = self.rmg.species_constraints["maximumCarbonAtoms"]
        max_heavy_atoms = self.rmg.species_constraints["maximumHeavyAtoms"]

        self.rmg.species_constraints["maximumCarbonAtoms"] = 3
        self.rmg.species_constraints["maximumHeavyAtoms"] = 6

        assert not fails_species_constraints(mol_1site)

        assert not fails_species_constraints(mol_2site)

        assert fails_species_constraints(mol_3site_vdW)

        assert fails_species_constraints(mol_3site)

        self.rmg.species_constraints["maximumCarbonAtoms"] = max_carbon
        self.rmg.species_constraints["maximumHeavyAtoms"] = max_heavy_atoms

    def test_surface_bond_order_constraint(self):
        """
        Test that we can constrain the max bond order of surface sites.
        """
        mol_1site = Molecule().from_adjacency_list(
            """
1 C u0 p0 c0 {2,Q}
2 X u0 p0 c0 {1,Q}
"""
        )
        assert fails_species_constraints(mol_1site)

    def test_heavy_constraint(self):
        """
        Test that we can constrain the max number of heavy atoms.
        """
        mol1 = Molecule(smiles="CCO")
        assert not fails_species_constraints(mol1)

        mol2 = Molecule(smiles="CCN=O")
        assert fails_species_constraints(mol2)

    def test_radical_constraint(self):
        """
        Test that we can constrain the max number of radical electrons.
        """
        mol1 = Molecule(smiles="[CH2][CH2]")
        assert not fails_species_constraints(mol1)

        mol2 = Molecule(smiles="[CH2][CH][CH2]")
        assert fails_species_constraints(mol2)

    def test_carbene_constraint(self):
        """
        Test that we can constrain the max number of singlet carbenes.
        """
        mol1 = Molecule().from_adjacency_list(
            """
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""
        )
        assert not fails_species_constraints(mol1)

        mol2 = Molecule().from_adjacency_list(
            """
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 C u0 p1 c0 {1,S} {4,S}
4 H u0 p0 c0 {3,S}
"""
        )
        assert fails_species_constraints(mol2)

    def test_carbene_radical_constraint(self):
        """
        Test that we can constrain the max number of radical electrons with a carbene.
        """
        mol1 = Molecule().from_adjacency_list(
            """
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""
        )
        assert not fails_species_constraints(mol1)

        mol2 = Molecule().from_adjacency_list(
            """
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 C u1 p0 c0 {1,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
"""
        )
        assert fails_species_constraints(mol2)


class TestPolymerConstraints:
    """Unit tests for size-based polymer/gas constraint routing."""

    @classmethod
    def setup_class(cls):
        cls.rmg = RMG()
        rmgpy.rmg.input.rmg = cls.rmg
        # Strict gas cap so polymer-routed species are discriminated from gas-routed ones.
        rmgpy.rmg.input.generated_species_constraints(maximumCarbonAtoms=2)

    @classmethod
    def teardown_class(cls):
        rmgpy.rmg.input.rmg = None

    def setup_method(self):
        self.rmg.polymer_constraints = None
        reset_polymer_warning()

    # --- membership predicate (size-based) ---

    def test_small_species_not_polymer_member(self):
        pc = {"maximumCarbonAtoms": 30}
        assert not is_polymer_constraint_member(Molecule(smiles="CC"), pc)          # 2 heavy
        assert not is_polymer_constraint_member(Species().from_smiles("CCCCCCCCCC"), pc)  # 10 heavy < 15

    def test_large_species_is_polymer_member(self):
        pc = {"maximumCarbonAtoms": 30}
        assert is_polymer_constraint_member(Molecule(smiles="C" * 20), pc)          # 20 heavy >= 15
        assert is_polymer_constraint_member(Species().from_smiles("C" * 20), pc)

    def test_threshold_override(self):
        pc = {"maximumCarbonAtoms": 30, "polymerSizeThreshold": 8}
        assert is_polymer_constraint_member(Molecule(smiles="CCCCCCCCCC"), pc)      # 10 heavy >= 8

    def test_threshold_boundary_is_inclusive(self):
        # The heavy >= threshold comparison is inclusive at the exact line, and the
        # 15-heavy default is pinned to CKMG's max_heavy_atoms_for_gas_species -- an
        # off-by-one here silently misaligns the polymer/pyrolysis handoff.
        default_pc = {"maximumCarbonAtoms": 30}
        assert is_polymer_constraint_member(Molecule(smiles="C" * 15), default_pc)       # 15 heavy == 15
        assert not is_polymer_constraint_member(Molecule(smiles="C" * 14), default_pc)   # 14 heavy < 15
        override_pc = {"maximumCarbonAtoms": 30, "polymerSizeThreshold": 10}
        assert is_polymer_constraint_member(Molecule(smiles="C" * 10), override_pc)      # 10 heavy == 10
        assert not is_polymer_constraint_member(Molecule(smiles="C" * 9), override_pc)   # 9 heavy < 10

    def test_is_polymer_flag_fast_path(self):
        # Real Polymer objects carry is_polymer=True and route to polymer regardless of size
        # and regardless of whether a polymer block is configured.
        assert is_polymer_constraint_member(mock.Mock(is_polymer=True), None)

    def test_proxy_flag_no_longer_routes(self):
        # Regression: the leaky is_polymer_proxy flag must NOT route a tiny species to polymer.
        pc = {"maximumCarbonAtoms": 30}
        mol = Molecule(smiles="CC")
        mol.is_polymer_proxy = True
        assert not is_polymer_constraint_member(mol, pc)                            # 2 heavy, flag ignored

    def test_no_block_only_real_polymer_is_member(self):
        # With no polymer block, size routing is inert: only real Polymer objects are members.
        assert not is_polymer_constraint_member(Molecule(smiles="C" * 20), None)

    # --- _evaluate_constraints purity (unchanged) ---

    def test_evaluate_constraints_is_pure(self):
        mol = Molecule(smiles="CCC")  # 3 carbons
        assert _evaluate_constraints(mol, {"maximumCarbonAtoms": 2})
        assert not _evaluate_constraints(mol, {"maximumCarbonAtoms": 10})

    # --- fails_species_constraints end-to-end routing ---

    def test_large_species_admitted_under_polymer_cap(self):
        # C20 (20 heavy) fails the gas cap (2) but passes a polymer cap of 30 -> took polymer path.
        self.rmg.polymer_constraints = {"maximumCarbonAtoms": 30}
        assert not fails_species_constraints(Molecule(smiles="C" * 20))

    def test_large_species_rejected_above_polymer_cap(self):
        self.rmg.polymer_constraints = {"maximumCarbonAtoms": 5}
        reason = fails_species_constraints(Molecule(smiles="C" * 20))
        assert reason and "maximumCarbonAtoms" in reason

    def test_species_level_large_routes_to_polymer_cap(self):
        self.rmg.polymer_constraints = {"maximumCarbonAtoms": 5}
        reason = fails_species_constraints(Species().from_smiles("C" * 20))
        assert reason and "maximumCarbonAtoms" in reason

    def test_species_level_large_admitted_under_polymer_cap(self):
        self.rmg.polymer_constraints = {"maximumCarbonAtoms": 30}
        assert not fails_species_constraints(Species().from_smiles("C" * 20))

    def test_small_species_uses_gas_cap_regardless_of_polymer_block(self):
        self.rmg.polymer_constraints = {"maximumCarbonAtoms": 30}
        assert not fails_species_constraints(Molecule(smiles="CC"))    # 2 C ok under gas cap 2
        assert fails_species_constraints(Molecule(smiles="CCC"))        # 3 C fails gas cap 2

    def test_threshold_override_reroutes_to_polymer(self):
        # Lower threshold to 8 -> C10 routes to polymer and escapes the strict gas cap.
        self.rmg.polymer_constraints = {"maximumCarbonAtoms": 30, "polymerSizeThreshold": 8}
        assert not fails_species_constraints(Molecule(smiles="CCCCCCCCCC"))  # 10 C, gas cap 2 would reject

    @mock.patch("rmgpy.constraints.logging")
    def test_absent_block_bypasses_only_real_polymer_and_warns_once(self, mock_logging):
        self.rmg.polymer_constraints = None
        poly = mock.Mock(is_polymer=True)
        assert not fails_species_constraints(poly)   # unbounded bypass
        assert not fails_species_constraints(poly)   # still bypass, no second warn
        mock_logging.warning.assert_called_once()
        reset_polymer_warning()
        assert not fails_species_constraints(poly)
        assert mock_logging.warning.call_count == 2
