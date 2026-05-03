#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
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
Tests for rmgpy.yaml_cantera1 module.
"""

import copy
import os
import pytest
import numpy as np
import yaml

from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.transport import TransportData
from rmgpy.kinetics import (
    Arrhenius,
    PDepArrhenius,
    Troe,
    ThirdBody,
)
from rmgpy.kinetics.surface import SurfaceArrhenius, StickingCoefficient
from rmgpy.yaml_cantera1 import (
    CanteraWriter1,
    species_to_dict,
    reaction_to_dicts,
)


def _make_nasa_thermo():
    coeffs = [1.0, 0.0, 0.0, 0.0, 0.0, -100.0, 1.0]
    return NASA(
        polynomials=[
            NASAPolynomial(coeffs=coeffs, Tmin=(200, "K"), Tmax=(1000, "K")),
            NASAPolynomial(coeffs=coeffs, Tmin=(1000, "K"), Tmax=(6000, "K")),
        ],
        Tmin=(200, "K"),
        Tmax=(6000, "K"),
    )


def _make_gas_species(label, smiles, index):
    sp = Species(label=label, index=index)
    sp.from_smiles(smiles)
    sp.thermo = _make_nasa_thermo()
    sp.transport_data = TransportData(
        shapeIndex=1,
        sigma=(3.0, "angstrom"),
        epsilon=(100.0, "K"),
        dipoleMoment=(0.0, "De"),
        polarizability=(0.0, "angstrom^3"),
        rotrelaxcollnum=1.0,
    )
    return sp


def _make_surface_species(label, adjlist, index):
    sp = Species(label=label, index=index)
    sp.from_adjacency_list(adjlist)
    sp.thermo = _make_nasa_thermo()
    return sp


class TestYamlCantera1Functions:
    """Unit tests for the individual helper functions in yaml_cantera1."""

    # ------------------------------------------------------------------
    # Shared fixtures
    # ------------------------------------------------------------------
    @pytest.fixture(autouse=True)
    def _build_species(self):
        self.h2 = _make_gas_species("H2", "[H][H]", index=1)
        self.h = _make_gas_species("H", "[H]", index=2)
        self.ar = _make_gas_species("Ar", "[Ar]", index=3)
        self.x = _make_surface_species("X", "1 X u0 p0", index=4)
        self.hx = _make_surface_species(
            "H_X", "1 H u0 p0 {2,S}\n2 X u0 p0 {1,S}", index=5
        )
        self.all_gas = [self.h2, self.h, self.ar]
        self.all_surface = [self.x, self.hx, self.h2]

    # ------------------------------------------------------------------
    # species_to_dict
    # ------------------------------------------------------------------
    def test_species_to_dict_gas_name_and_thermo(self):
        """Gas species: correct name, NASA7 thermo with two polynomial ranges."""
        d = species_to_dict(self.h2)
        assert d["name"] == "H2(1)"
        assert d["composition"] == {"H": 2.0}
        assert d["thermo"]["model"] == "NASA7"
        assert len(d["thermo"]["temperature-ranges"]) == 3  # low, mid, high
        assert len(d["thermo"]["data"]) == 2

    def test_species_to_dict_gas_transport(self):
        """Gas species transport data is present and contains geometry."""
        d = species_to_dict(self.h2)
        assert "transport" in d
        assert d["transport"]["model"] == "gas"
        assert d["transport"]["geometry"] == "linear"
        assert np.isclose(d["transport"]["diameter"], 3.0)
        assert np.isclose(d["transport"]["well-depth"], 100.0)

    def test_species_to_dict_surface_composition(self):
        """Surface species has X in composition and no transport block."""
        d = species_to_dict(self.hx)
        assert "X" in d["composition"]
        assert d["composition"]["X"] == 1.0
        assert "H" in d["composition"]
        assert "transport" not in d

    def test_species_to_dict_surface_thermo_model(self):
        """Surface species reports NASA7 thermo."""
        d = species_to_dict(self.x)
        assert d["thermo"]["model"] == "NASA7"

    # ------------------------------------------------------------------
    # reaction_to_dicts — gas-phase kinetics
    # Units declared in yaml_cantera1: activation-energy: J/kmol,
    # so all Ea values come through multiplied by 1000 relative to J/mol.
    # ------------------------------------------------------------------
    def test_reaction_to_dicts_arrhenius_equation_and_rate(self):
        """Arrhenius: equation string and rate-constant keys present with J/kmol Ea."""
        kin = Arrhenius(A=(1e13, "s^-1"), n=0.5, Ea=(10, "kJ/mol"), T0=(1, "K"))
        rxn = Reaction(reactants=[self.h2], products=[self.h, self.h], kinetics=kin)
        entries = reaction_to_dicts(rxn, self.all_gas)
        assert len(entries) == 1
        d = entries[0]
        assert d["equation"] == "H2(1) <=> 2 H(2)"
        assert "rate-constant" in d
        assert np.isclose(d["rate-constant"]["A"], 1e13)
        assert np.isclose(d["rate-constant"]["b"], 0.5)
        assert np.isclose(d["rate-constant"]["Ea"], 10e6)  # 10 kJ/mol → 1e7 J/kmol

    def test_reaction_to_dicts_thirdbody(self):
        """ThirdBody: equation uses M, efficiencies map present."""
        kin = ThirdBody(
            arrheniusLow=Arrhenius(
                A=(1e18, "cm^6/(mol^2*s)"), n=-1, Ea=(0, "J/mol"), T0=(1, "K")
            ),
            efficiencies={self.ar.molecule[0]: 0.7},
        )
        rxn = Reaction(
            reactants=[self.h, self.h], products=[self.h2], kinetics=kin
        )
        entries = reaction_to_dicts(rxn, self.all_gas)
        d = entries[0]
        assert "M" in d["equation"]
        assert "rate-constant" in d
        assert "efficiencies" in d
        assert np.isclose(d["efficiencies"]["Ar(3)"], 0.7)

    def test_reaction_to_dicts_pdep_arrhenius(self):
        """PDepArrhenius: type is pressure-dependent-Arrhenius, rate-constants list."""
        kin = PDepArrhenius(
            pressures=([0.1, 1.0], "atm"),
            arrhenius=[
                Arrhenius(A=(1e10, "s^-1"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K")),
                Arrhenius(A=(1e12, "s^-1"), n=0, Ea=(15, "kJ/mol"), T0=(1, "K")),
            ],
        )
        rxn = Reaction(reactants=[self.h2], products=[self.h, self.h], kinetics=kin)
        entries = reaction_to_dicts(rxn, self.all_gas)
        d = entries[0]
        assert d["type"] == "pressure-dependent-Arrhenius"
        rates = d["rate-constants"]
        assert len(rates) == 2
        assert np.isclose(rates[0]["P"], 0.1 * 101325.0)
        assert np.isclose(rates[0]["A"], 1e10)
        assert np.isclose(rates[0]["Ea"], 10e6)   # J/kmol
        assert np.isclose(rates[1]["P"], 1.0 * 101325.0)
        assert np.isclose(rates[1]["A"], 1e12)
        assert np.isclose(rates[1]["Ea"], 15e6)   # J/kmol

    def test_reaction_to_dicts_troe(self):
        """Troe: type falloff, Troe block present, high/low rate constants."""
        kin = Troe(
            arrheniusHigh=Arrhenius(
                A=(1e14, "s^-1"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K")
            ),
            arrheniusLow=Arrhenius(
                A=(1e20, "cm^3/(mol*s)"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K")
            ),
            alpha=0.5,
            T3=(100, "K"),
            T1=(200, "K"),
            T2=(300, "K"),
            efficiencies={self.ar.molecule[0]: 2.0},
        )
        rxn = Reaction(reactants=[self.h], products=[self.h], kinetics=kin)
        entries = reaction_to_dicts(rxn, self.all_gas)
        d = entries[0]
        assert d["type"] == "falloff"
        assert "Troe" in d
        assert np.isclose(d["Troe"]["A"], 0.5)
        assert np.isclose(d["Troe"]["T3"], 100.0)
        assert np.isclose(d["Troe"]["T1"], 200.0)
        assert np.isclose(d["Troe"]["T2"], 300.0)
        assert np.isclose(d["high-P-rate-constant"]["A"], 1e14)
        assert np.isclose(d["high-P-rate-constant"]["Ea"], 10e6)  # J/kmol
        assert "efficiencies" in d
        assert np.isclose(d["efficiencies"]["Ar(3)"], 2.0)

    # ------------------------------------------------------------------
    # reaction_to_dicts — surface kinetics
    # ------------------------------------------------------------------
    def test_reaction_to_dicts_surface_arrhenius(self):
        """SurfaceArrhenius: rate-constant present, A converted to /kmol."""
        kin = SurfaceArrhenius(
            A=(1e13, "m^2/(mol*s)"), n=0.5, Ea=(50, "kJ/mol"), T0=(1, "K")
        )
        rxn = Reaction(
            reactants=[self.h2, self.x],
            products=[self.hx, self.hx],
            kinetics=kin,
        )
        entries = reaction_to_dicts(rxn, self.all_surface)
        d = entries[0]
        assert "rate-constant" in d
        # A in m^2/(mol*s) → ×1000 → m^2/(kmol*s)
        assert np.isclose(d["rate-constant"]["A"], 1e16)
        assert np.isclose(d["rate-constant"]["b"], 0.5)
        assert np.isclose(d["rate-constant"]["Ea"], 50e6)  # J/kmol

    def test_reaction_to_dicts_sticking_coefficient(self):
        """StickingCoefficient: sticking-coefficient block present, A dimensionless."""
        kin = StickingCoefficient(
            A=(0.1, ""), n=0, Ea=(0, "kJ/mol"), T0=(1, "K")
        )
        rxn = Reaction(
            reactants=[self.h2, self.x, self.x],
            products=[self.hx, self.hx],
            kinetics=kin,
        )
        entries = reaction_to_dicts(rxn, self.all_surface)
        d = entries[0]
        assert "sticking-coefficient" in d
        assert np.isclose(d["sticking-coefficient"]["A"], 0.1)
        assert np.isclose(d["sticking-coefficient"]["Ea"], 0.0)

    def test_reaction_to_dicts_coverage_dependence(self):
        """Coverage-dependent kinetics: coverage-dependencies block present with correct units.

        yaml_cantera1 declares activation-energy: J/kmol, so E must be ×1000
        relative to RMG's J/mol value_si.
        """
        kin = StickingCoefficient(
            A=(0.1, ""),
            n=0,
            Ea=(0, "kJ/mol"),
            T0=(1, "K"),
            coverage_dependence={
                self.hx: {"a": 0.5, "m": -1.0, "E": (5.0, "kJ/mol")}
            },
        )
        rxn = Reaction(
            reactants=[self.h2, self.x, self.x],
            products=[self.hx, self.hx],
            kinetics=kin,
        )
        entries = reaction_to_dicts(rxn, self.all_surface)
        d = entries[0]
        assert "coverage-dependencies" in d
        cov = d["coverage-dependencies"]["H_X(5)"]
        assert np.isclose(cov["a"], 0.5)
        assert np.isclose(cov["m"], -1.0)
        # 5 kJ/mol = 5000 J/mol → ×1000 → 5 000 000 J/kmol
        assert np.isclose(cov["E"], 5e6)


class TestCanteraWriter1:
    """Tests for the CanteraWriter1 class."""

    def test_can_instantiate(self):
        """Test that CanteraWriter1 can be instantiated."""
        writer = CanteraWriter1()
        assert writer is not None

class CanteraYamlFileComparer:
    """
    For comparing two Cantera YAML files. 
    This class provides methods to compare species and reactions between the two files.

        Args:
            yaml_path_1: Path to the first YAML file, converted from Chemkin by ck2yaml.
            yaml_path_2: Path to the second YAML file, written directly by RMG.
    """
    yaml_path_1 = None
    yaml_path_2 = None

    @pytest.fixture(autouse=True, scope="class") # loaded once per Class
    def load_yaml_files(self, request):
        """Load the two YAML files to be compared."""
        with open(request.cls.yaml_path_1, 'r') as file:
            request.cls.yaml1 = yaml.safe_load(file)
        with open(request.cls.yaml_path_2, 'r') as file:
            request.cls.yaml2 = yaml.safe_load(file)

    @pytest.fixture(autouse=True)  # runs before each test method
    def copy_yaml_dicts(self):
        """Make deep copies so tests can modify without affecting other tests."""
        self.yaml1 = copy.deepcopy(self.__class__.yaml1)
        self.yaml2 = copy.deepcopy(self.__class__.yaml2)

    def testGeneratorsAsExpected(self):
        "Check the two yaml files were generated by the expected tools (ck2yaml vs RMG)."
        assert self.yaml1['generator'] == 'ck2yaml', "First YAML file should be generated by ck2yaml."
        assert self.yaml2['generator'] == 'RMG', "Second YAML file should be generated by RMG."

    def testKeysMatch(self):
        """Test that the top-level keys in both YAML files match, except those expected not to."""
        # Remove keys from ck2yaml output that are not present in RMG output
        self.yaml1.pop('input-files', None)
        self.yaml1.pop('cantera-version', None)
        for model in [self.yaml1, self.yaml2]:
            for phase in model['phases']:
                for reactions_block in phase.get('reactions', []): # for multi-phase mechanisms, reactions are under each phase
                    assert reactions_block in model, f"Expected reactions block '{reactions_block}' not found in YAML file."
                    model.pop(reactions_block, None)  # Remove reactions block to allow keys to match
        assert self.yaml1.keys() == self.yaml2.keys(), "YAML files have different top-level keys."
    
    def testPhasesMatch(self):
        """Test that the phase definitions in both YAML files match."""
        assert len(self.yaml1['phases']) == len(self.yaml2['phases']), "YAML files have different numbers of phases"

        for phase1, phase2 in zip(self.yaml1['phases'], self.yaml2['phases']):
            assert phase1['name'] == phase2['name'], f"Phase names do not match: {phase1['name']} vs {phase2['name']}."
            assert phase1['thermo'] == phase2['thermo'], f"Thermo definitions for phase {phase1['name']} do not match."
            assert phase1.get('transport', '') == phase2.get('transport', ''), f"Transport definitions for phase {phase1['name']} do not match."
            assert phase1.get('adjacent-phases', []) == phase2.get('adjacent-phases', []), f"Adjacent phases for phase {phase1['name']} do not match."
            assert phase1.get('species', []) == phase2.get('species', []), f"Species lists for phase {phase1['name']} do not match."
            assert phase1.get('reactions', []) == phase2.get('reactions', []), f"Reactions blocks for phase {phase1['name']} do not match."
            # the ck2yaml has all elements in Titlecase, while RMG lets some isotopes be CI and OI (not Ci and Oi).
            assert sorted(phase1.get('elements', [])) == sorted(e.title() for e in phase2.get('elements', [])), f"Element lists for phase {phase1['name']} do not match."
            assert phase1.get('state', {}) == phase2.get('state', {}), f"State definitions for phase {phase1['name']} do not match."

    def testElementsMatch(self):
        """Test that the element definitions in both YAML files match."""
        ck2yaml_elements = sorted(self.yaml1['elements'], key=lambda e: e['symbol'])
        # Put symbol into Titlecase to match ck2yaml's formatting
        rmg_elements = [{'symbol': e['symbol'].title(), 'atomic-weight': e['atomic-weight']} for e in self.yaml2['elements']]
        # Sort by the 'symbol' key.
        rmg_elements = sorted(rmg_elements, key=lambda e: e['symbol'])
        # Compare symbols exactly, and atomic weights approximately
        assert [e['symbol'] for e in ck2yaml_elements] == [e['symbol'] for e in rmg_elements], \
            "YAML files have different element symbols."
        assert [e['atomic-weight'] for e in ck2yaml_elements] == pytest.approx(
            [e['atomic-weight'] for e in rmg_elements], abs=1e-3
        ), "YAML files have different element atomic weights."

    def testSpeciesMatch(self):
        """Test that species definitions match between the two YAML files."""
        species1 = {s['name']: s for s in self.yaml1['species']}
        species2 = {s['name']: s for s in self.yaml2['species']}
        assert species1.keys() == species2.keys(), "Species names do not match."

        for name in species1:
            s1 = species1[name]
            s2 = species2[name]

            # Composition: ck2yaml uses int values, RMG uses float
            assert {k: int(v) for k, v in s2['composition'].items()} == s1['composition'], \
                f"Composition mismatch for {name}."

            # Thermo model
            assert s1['thermo']['model'] == s2['thermo']['model'], \
                f"Thermo model mismatch for {name}."

            # Temperature ranges and polynomial data
            # ck2yaml may collapse single-polynomial NASA7 (e.g. Ar) into one range
            # while RMG always writes two polynomials with a midpoint temperature.
            t_ranges1 = s1['thermo'].get('temperature-ranges', [])
            t_ranges2 = s2['thermo'].get('temperature-ranges', [])
            data1 = s1['thermo'].get('data', [])
            data2 = s2['thermo'].get('data', [])

            if len(t_ranges1) == 2 and len(t_ranges2) == 3:
                # ck2yaml collapsed to single polynomial; RMG has two identical ones
                assert t_ranges1[0] == pytest.approx(t_ranges2[0], rel=1e-4), \
                    f"Temperature range lower bound mismatch for {name}."
                assert t_ranges1[1] == pytest.approx(t_ranges2[2], rel=1e-4), \
                    f"Temperature range upper bound mismatch for {name}."
                assert len(data1) == 1 and len(data2) == 2, \
                    f"Expected 1 vs 2 polynomials for collapsed species {name}."
                assert data1[0] == pytest.approx(data2[0], rel=1e-4), \
                    f"Thermo polynomial mismatch for {name} (low range)."
                assert data1[0] == pytest.approx(data2[1], rel=1e-4), \
                    f"Thermo polynomial mismatch for {name} (high range should match low)."
            else:
                assert t_ranges1 == pytest.approx(t_ranges2, rel=1e-4), \
                    f"Temperature ranges mismatch for {name}."
                assert len(data1) == len(data2), \
                    f"Number of thermo polynomial ranges differs for {name}."
                for i, (poly1, poly2) in enumerate(zip(data1, data2)):
                    assert poly1 == pytest.approx(poly2, rel=1e-4), \
                        f"Thermo polynomial {i} mismatch for {name}."
            # Ideally thermo data would have notes.

            # RMG includes reference-pressure but ck2yaml does not (when it's non-default)
            # (no assertion needed, just noting the known difference)

            # Transport data
            assert ('transport' in s1) == ('transport' in s2), f"Transport data presence mismatch for {name}."
            if 'transport' in s1 and 'transport' in s2:
                t1 = s1['transport']
                t2 = s2['transport']
                assert t1['model'] == t2['model'], f"Transport model mismatch for {name}."
                assert t1['geometry'] == t2['geometry'], f"Transport geometry mismatch for {name}."
                assert t1.get('well-depth', 0) == pytest.approx(
                    t2.get('well-depth', 0), rel=1e-3
                ), f"Transport well-depth mismatch for {name}."
                assert t1.get('diameter', 0) == pytest.approx(
                    t2.get('diameter', 0), rel=1e-3
                ), f"Transport diameter mismatch for {name}."
                assert t1.get('polarizability', 0) == pytest.approx(
                    t2.get('polarizability', 0), rel=1e-3
                ), f"Transport polarizability mismatch for {name}."
                assert t1.get('dipole', 0) == pytest.approx(
                    t2.get('dipole', 0), rel=1e-3
                ), f"Transport dipole mismatch for {name}."
                assert t1.get('rotational-relaxation', 0) == pytest.approx(
                    t2.get('rotational-relaxation', 0), rel=1e-3
                ), f"Transport rotational-relaxation mismatch for {name}."
                assert t1.get('note', '') == t2.get('note', ''), \
                    f"Transport note mismatch for {name}."

@pytest.mark.skip(reason="These files are out of date and have been removed.")
class TestPreviouslyWrittenCanteraYamlGasOnly(CanteraYamlFileComparer):
    """Tests for comparing previously written Cantera YAML files, gas-only mechanism.

    These are stored in the testing data directory.
    """
    test_data_folder='test/rmgpy/test_data/yaml_writer_data/'
    # generated on the fly in recent functional test
    yaml_path_1 = os.path.join(test_data_folder, 'chemkin/chem37.yaml')
    yaml_path_2 = os.path.join(test_data_folder, 'cantera1/chem37.yaml')

class TestRecentlyGeneratedCanteraYamlGasOnly(CanteraYamlFileComparer):
    """Tests for comparing recently generated Cantera YAML files, gas-only mechanism.

    These are generated on the fly in the mainTest.py functional test and stored in the testing data directory.
    """
    test_data_folder='test/rmgpy/test_data/yaml_writer_data/'
    
    @pytest.fixture(autouse=True, scope="class")
    def find_recent_files(self, request):
        """Find the YAML files generated by mainTest."""
        cantera_dir = os.path.join(self.test_data_folder, 'cantera1')
        chemkin_dir = os.path.join(self.test_data_folder, 'chemkin')
        
        if not os.path.exists(cantera_dir) or not os.path.exists(chemkin_dir):
            pytest.skip("YAML test data directories not found. Run mainTest first.")
        
        # Look for specifically named files from mainTest
        cantera_file = os.path.join(cantera_dir, 'from_main_test.yaml')
        chemkin_file = os.path.join(chemkin_dir, 'from_main_test.yaml')
        
        if not os.path.exists(cantera_file) or not os.path.exists(chemkin_file):
            pytest.skip("from_main_test.yaml files not found. Run mainTest first.")
        
        request.cls.yaml_path_1 = chemkin_file
        request.cls.yaml_path_2 = cantera_file

@pytest.mark.skip(reason="These files are out of date and have been removed.")
class TestPreviouslyWrittenCanteraYamlWithSurface(CanteraYamlFileComparer):
    """Tests for comparing previously written Cantera YAML files, with surface mechanism.

    These are stored in the testing data directory.
    """
    test_data_folder='test/rmgpy/test_data/yaml_writer_data/'
    # saved by Prosper in earlier commit
    yaml_path_1 = os.path.join(test_data_folder, 'chemkin/chem0047-gas.yaml')
    yaml_path_2 = os.path.join(test_data_folder, 'cantera1/chem47.yaml')
