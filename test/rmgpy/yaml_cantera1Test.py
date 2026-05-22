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

import os
import pytest
import numpy as np

from cantera_yaml_comparer import CanteraYamlFileComparer
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.transport import TransportData
from rmgpy.kinetics import (
    Arrhenius,
    Lindemann,
    PDepArrhenius,
    Troe,
    ThirdBody,
)
from rmgpy.kinetics.surface import SurfaceArrhenius, StickingCoefficient
from rmgpy.yaml_cantera1 import (
    CanteraWriter1,
    species_to_dict,
    reaction_to_dicts,
    get_elements_block,
    get_phases_gas_only,
    get_phases_with_surface,
    get_mech_dict_nonsurface,
    get_mech_dict_surface,
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

    def test_species_to_dict_charged_species_uses_electron_pseudo_element(self):
        """Charged species composition uses Cantera's E pseudo-element."""
        lithium_ion = Species(label="Li+", index=10)
        lithium_ion.from_adjacency_list("1 Li u0 p0 c+1")
        lithium_ion.thermo = _make_nasa_thermo()

        d = species_to_dict(lithium_ion)

        assert d["composition"] == {"E": -1.0, "Li": 1.0}

    def test_species_to_dict_electron_uses_electron_pseudo_element(self):
        """The electron species exports composition as E, not internal e."""
        electron = Species(label="e", index=11)
        electron.from_adjacency_list("1 e u1 p0 c-1")
        electron.thermo = _make_nasa_thermo()

        d = species_to_dict(electron)

        assert d["composition"] == {"E": 1.0}

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
        assert d["equation"] == "H2(1) <=> H(2) + H(2)"
        assert "rate-constant" in d
        assert np.isclose(d["rate-constant"]["A"], 1e13)
        assert np.isclose(d["rate-constant"]["b"], 0.5)
        assert np.isclose(d["rate-constant"]["Ea"], 10e6)  # 10 kJ/mol → 1e7 J/kmol

    def test_reaction_to_dicts_negative_a_arrhenius(self):
        """Negative Arrhenius A factors are marked for Cantera."""
        kin = Arrhenius(A=(-1e13, "s^-1"), n=0.5, Ea=(10, "kJ/mol"), T0=(1, "K"))
        rxn = Reaction(reactants=[self.h2], products=[self.h, self.h], kinetics=kin)
        entries = reaction_to_dicts(rxn, self.all_gas)
        assert entries[0]["negative-A"] is True

    def test_reaction_to_dicts_thirdbody(self):
        """ThirdBody: type is three-body, equation uses M, efficiencies map present."""
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
        assert d["type"] == "three-body"
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

    def test_reaction_to_dicts_negative_a_falloff(self):
        """Falloff rates with negative high- and low-pressure A factors are marked for Cantera."""
        kin = Lindemann(
            arrheniusHigh=Arrhenius(
                A=(-1e14, "s^-1"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K")
            ),
            arrheniusLow=Arrhenius(
                A=(-1e20, "cm^3/(mol*s)"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K")
            ),
            efficiencies={self.ar.molecule[0]: 2.0},
        )
        rxn = Reaction(reactants=[self.h], products=[self.h], kinetics=kin)
        entries = reaction_to_dicts(rxn, self.all_gas)
        d = entries[0]
        assert d["negative-A"] is True
        assert d["low-P-rate-constant"]["A"] < 0
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

    def test_reaction_to_dicts_negative_a_sticking_coefficient(self):
        """Negative sticking-coefficient A factors are marked for Cantera."""
        kin = StickingCoefficient(
            A=(-0.1, ""), n=0, Ea=(0, "kJ/mol"), T0=(1, "K")
        )
        rxn = Reaction(
            reactants=[self.h2, self.x, self.x],
            products=[self.hx, self.hx],
            kinetics=kin,
        )
        entries = reaction_to_dicts(rxn, self.all_surface)
        assert entries[0]["negative-A"] is True

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

    def test_reaction_to_dicts_surface_spectator_species(self):
        """SurfaceArrhenius with a spectator species on both sides must not produce efficiencies.

        Cantera's Python API misidentifies a species with equal stoichiometry on
        both sides (a net spectator / surface catalyst) as a third-body collider,
        causing input_data to include a spurious 'efficiencies' entry and doubled
        stoichiometry in the equation string.  reaction_to_dicts must detect and
        fix this.
        """
        # Build a second surface species to act as spectator (adsorbed oxygen)
        ox = _make_surface_species(
            "O_X",
            "1 O u0 p2 c0 {2,D}\n2 X u0 p0 c0 {1,D}",
            index=6,
        )
        # H_X(5) + O_X(6) <=> X(4) + O_X(6)  — O_X is spectator on both sides
        kin = SurfaceArrhenius(
            A=(4.18e20, "m^2/(mol*s)"), n=0.0, Ea=(148.7, "kJ/mol"), T0=(1, "K")
        )
        rxn = Reaction(
            reactants=[self.hx, ox],
            products=[self.x, ox],
            kinetics=kin,
        )
        species_list = self.all_surface + [ox]
        entries = reaction_to_dicts(rxn, species_list)
        d = entries[0]
        assert "efficiencies" not in d, (
            "Spurious 'efficiencies' must not appear for a SurfaceArrhenius reaction "
            "even when a species appears on both sides."
        )
        eq = d["equation"]
        # O_X(6) should appear exactly once on each side, not doubled
        assert eq.count("O_X(6)") == 2, (
            f"Spectator O_X should appear once per side in: {eq}"
        )
        assert "2 O_X" not in eq, (
            f"Spectator stoichiometry should not be doubled in: {eq}"
        )

    def test_reaction_to_dicts_three_species_one_side_spectator(self):
        """Spectator on both sides with 3+ stoichiometric items on a side must still drop efficiencies.

        Cantera's API (v. 3.1 and 3.2) misidentifies a species with net-zero stoichiometry as
        a third-body collider whenever the reaction has three or more
        stoichiometric items on one side. Routing ct.Reaction through an
        equation string avoids this only for the 2-each-side case. For
        wider reactions the writer must strip the resulting spurious
        'efficiencies' from input_data so the YAML round-trips through
        ct.Solution.
        See https://github.com/Cantera/cantera/issues/2115#issuecomment-4465564540
        """
        ox = _make_surface_species(
            "O_X",
            "1 O u0 p2 c0 {2,D}\n2 X u0 p0 c0 {1,D}",
            index=6,
        )
        # Recombinative H2 desorption with an adjacent O_X site as a chemically
        # inert spectator:
        #   2 H_X + O_X <=> H2 + O_X + 2 X    (atom balance: 2 H + O + 3 X)
        # O_X has stoichiometry 1 on each side (net zero), so Cantera flags it
        # as a third-body collider even though the underlying kinetics is a
        # plain SurfaceArrhenius.
        kin = SurfaceArrhenius(
            A=(3.73e19, "m^2/(mol*s)"), n=0.475, Ea=(33.6, "kJ/mol"), T0=(1, "K")
        )
        rxn = Reaction(
            reactants=[self.hx, self.hx, ox],
            products=[self.h2, ox, self.x, self.x],
            kinetics=kin,
        )
        species_list = self.all_surface + [ox]
        entries = reaction_to_dicts(rxn, species_list)
        d = entries[0]
        assert "efficiencies" not in d, (
            "Spurious 'efficiencies' must not appear for a SurfaceArrhenius "
            f"reaction with a net-zero-stoichiometry spectator. Got: {d}"
        )

    def test_get_elements_block_isotopes_and_surface_site(self):
        """get_elements_block emits isotope and X definitions only when requested."""
        from rmgpy.molecule.element import H, C, D, T, X, e

        # With D, T, X in use, the block names them with the right masses
        elements_block, elements_line = get_elements_block({H, C, D, T, X})
        assert 'X' in elements_line
        assert 'D' in elements_block    # H-2
        assert 'T' in elements_block    # H-3
        assert 'X' in elements_block
        assert '195.083' in elements_block  # X atomic weight

        # Without isotopes / X in use, neither appears
        elements_block, elements_line = get_elements_block({H, C})
        assert 'X' not in elements_line
        assert ' D ' not in elements_line and '[D' not in elements_line and ', D' not in elements_line
        assert 'D' not in elements_block
        assert 'T' not in elements_block
        assert 'X' not in elements_block

        # The RMG electron singleton is lowercase e internally, but Cantera
        # expects the uppercase pseudo-element E in phase element lists.
        elements_block, elements_line = get_elements_block({H, e})
        assert elements_block == ''
        assert elements_line == 'elements: [E, H]'

    def test_get_phases_gas_only_has_state(self):
        """Gas-only phases block includes state with T and P."""
        phases_block = get_phases_gas_only([self.h2, self.h], 'elements: [H]')
        assert 'state:' in phases_block
        assert 'T: 300.0' in phases_block
        assert 'P: 1 atm' in phases_block

    def test_get_phases_gas_only_has_elements(self):
        """Gas-only phases block includes elements line."""
        phases_block = get_phases_gas_only([self.h2], 'elements: [H]')
        assert 'elements:' in phases_block

    def test_get_mech_dict_nonsurface_reactions_key(self):
        """Gas-only mech dict uses 'reactions' key."""
        rxn = Reaction(
            reactants=[self.h2], products=[self.h, self.h],
            kinetics=Arrhenius(A=(1e13, "s^-1"), n=0, Ea=(400, "kJ/mol"), T0=(1, "K"))
        )
        result = get_mech_dict_nonsurface([self.h2, self.h], [rxn])
        assert 'species' in result
        assert 'reactions' in result
        assert len(result['reactions']) == 1

    def test_get_mech_dict_surface_has_gas_and_surface_reactions(self):
        """Surface mech dict separates gas-reactions and site0-reactions."""
        kin = SurfaceArrhenius(
            A=(1e13, "m^2/(mol*s)"), n=0, Ea=(50, "kJ/mol"), T0=(1, "K")
        )
        rxn = Reaction(reactants=[self.h2, self.x], products=[self.hx, self.hx], kinetics=kin)
        result = get_mech_dict_surface(
            [self.h2, self.x, self.hx], [rxn]
        )
        assert 'species' in result
        assert 'gas-reactions' in result
        assert 'site0-reactions' in result
        assert len(result['site0-reactions']) == 1

    def test_get_phases_with_surface_has_state(self):
        """Surface phases block includes state for both gas and surface phases."""
        phases_block = get_phases_with_surface(
            [self.h2, self.x, self.hx],
            surface_site_density=2.5e-9,  # mol/cm^2-equivalent SI
            elements_line='elements: [H, X]',
        )
        # Should have two phase definitions each with state
        assert phases_block.count('state:') == 2
        assert 'T: 300.0' in phases_block
        assert 'P: 1 atm' in phases_block

    def test_species_to_dict_surface_sites_count(self):
        """species_to_dict reports correct 'sites' count for multi-site surface species."""
        # Bidentate glyoxal adsorbed via C and O (2 X atoms)
        glyoxal_xx = _make_surface_species(
            "glyoxalXX",
            "1 O u0 p2 c0 {3,S} {8,S}\n"
            "2 O u0 p2 c0 {4,D}\n"
            "3 C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}\n"
            "4 C u0 p0 c0 {2,D} {3,S} {6,S}\n"
            "5 H u0 p0 c0 {3,S}\n"
            "6 H u0 p0 c0 {4,S}\n"
            "7 X u0 p0 c0 {3,S}\n"
            "8 X u0 p0 c0 {1,S}",
            index=215,
        )
        d2 = species_to_dict(glyoxal_xx)
        assert 'sites' in d2
        assert d2['sites'] == 2

    def test_species_to_dict_gas_no_sites_field(self):
        """Gas species must not have a 'sites' field."""
        d = species_to_dict(self.h2)
        assert 'sites' not in d

    def test_species_to_dict_transport_note_always_present(self):
        """Transport 'note' is always written when transport_data.comment is set."""
        self.h2.transport_data.comment = "from GRI-Mech"
        d = species_to_dict(self.h2)
        assert d['transport']['note'] == "from GRI-Mech"


class TestCanteraWriter1:
    """Tests for the CanteraWriter1 class."""

    def test_can_instantiate(self):
        """Test that CanteraWriter1 can be instantiated."""
        writer = CanteraWriter1()
        assert writer is not None


class TestRecentlyGeneratedCanteraYaml1GasOnly(CanteraYamlFileComparer):
    """Tests for comparing recently generated Cantera YAML files, gas-only mechanism.

    These are generated on the fly in the mainTest.py functional test and stored in the testing data directory.
    """
    test_data_folder='test/rmgpy/test_data/yaml_writer_data/'
    
    @pytest.fixture(autouse=True, scope="class")
    def find_recent_files(self, request):
        """
        Find the YAML files generated by mainTest.
        """
        cantera_file = os.path.join(self.test_data_folder, 'cantera1', 'from_main_test.yaml')
        chemkin_file = os.path.join(self.test_data_folder, 'ck2yaml', 'from_main_test.yaml')

        if not (os.path.exists(cantera_file) and os.path.exists(chemkin_file)):
            # If mainTest's copy step was collected for this pytest session but the
            # files are still missing, treat that as a failure rather than a skip —
            # it means mainTest ran but didn't produce the expected output.
            # (if using pytest-randomly or pytest-ordering to reorder them this will need altering).
            main_test_collected = any(
                item.name == "test_cantera_input_files_match_chemkin_later"
                for item in request.session.items
            )
            if main_test_collected:
                pytest.fail(
                    "from_main_test.yaml files missing even though mainTest's "
                    "copy step was collected — it likely failed before copying."
                )
            # If mainTest wasn't collected, it's likely that we're running this test 
            # in isolation, so skip without failing.
            pytest.skip("from_main_test.yaml files not found. Run mainTest first.")

        request.cls.yaml_path_1 = chemkin_file
        request.cls.yaml_path_2 = cantera_file

# This is kept here as an example of how to write a test comparing previously generated YAML files,
# but the files themselves are now out of date and have been removed,
# so the test is skipped and the class is commented out to avoid confusion.
# @pytest.mark.skip(reason="These files are out of date and have been removed.")
# class TestPreviouslyWrittenCanteraYamlWithSurface(CanteraYamlFileComparer):
#     """Tests for comparing previously written Cantera YAML files.

#     These are stored in the testing data directory.
#     """
#     test_data_folder='test/rmgpy/test_data/yaml_writer_data/'
#     # saved by Prosper in earlier commit
#     yaml_path_1 = os.path.join(test_data_folder, 'chemkin/chem1.yaml')
#     yaml_path_2 = os.path.join(test_data_folder, 'cantera1/chem1.yaml')
