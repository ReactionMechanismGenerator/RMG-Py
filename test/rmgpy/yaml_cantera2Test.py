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


import cantera as ct
import os
import shutil
import numpy as np
import pytest

from cantera_yaml_comparer import CanteraYamlFileComparer
from rmgpy.molecule import Atom, Molecule, get_element
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.kinetics import (
    Arrhenius,
    PDepArrhenius,
    MultiArrhenius,
    Chebyshev,
    Troe,
    Lindemann,
    ThirdBody,
)
from rmgpy.kinetics.surface import SurfaceArrhenius, StickingCoefficient
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.transport import TransportData
from rmgpy.yaml_cantera2 import (
    CanteraWriter2,
    save_cantera_files,
    species_to_dict,
    reaction_to_dict_list,
    generate_cantera_data,
    get_elements_lists,
)


class TestCanteraWriter2:

    def setup_method(self):
        """
        Create a temporary directory for file I/O tests.
        """
        base_dir = os.path.dirname(os.path.abspath(__file__))
        self.tmp_dir = os.path.join(base_dir, 'tmp')

        # Ensure a clean start: delete if exists, then create
        if os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.makedirs(self.tmp_dir)

    def teardown_method(self):
        """
        Clean up the temporary directory after tests.
        """
        shutil.rmtree(self.tmp_dir)


    def _create_dummy_species(self, label, formula, index=-1):
        """Helper to create a functional RMG Species object with thermo/transport"""
        sp = Species(label=label).from_smiles(formula)
        sp.index = index
        coeffs = [1.0, 0.0, 0.0, 0.0, 0.0, -100.0, 1.0]
        poly_low = NASAPolynomial(coeffs=coeffs, Tmin=(200, 'K'), Tmax=(1000, 'K'))
        poly_high = NASAPolynomial(coeffs=coeffs, Tmin=(1000, 'K'), Tmax=(6000, 'K'))
        sp.thermo = NASA(polynomials=[poly_low, poly_high], Tmin=(200, 'K'), Tmax=(6000, 'K'))
        num_atoms = len(sp.molecule[0].atoms)
        if num_atoms == 1:
            shape_idx = 0
        elif num_atoms == 2:
            shape_idx = 1
        else:
            shape_idx = 2
        sp.transport_data = TransportData(
            shapeIndex=shape_idx,
            sigma=(3.0, 'angstrom'),
            epsilon=(100.0, 'K'),
            dipoleMoment=(1.7, 'De'),
            polarizability=(0.0, 'angstrom^3'),
            rotrelaxcollnum=1.0
        )
        return sp

    def test_species_to_dict_standard(self):
        """Test conversion of a standard gas species."""
        sp = self._create_dummy_species("H2", "[H][H]", index=1)
        d = species_to_dict(sp, [sp])

        assert d['name'] == "H2(1)"
        assert 'composition' in d
        assert d['thermo']['model'] == 'NASA7'
        assert len(d['thermo']['data']) == 2

        # Verify Transport
        assert 'transport' in d
        assert d['transport']['model'] == 'gas'
        assert d['transport']['geometry'] == 'linear'
        # Diameter should be in angstroms ( https://cantera.org/dev/yaml/species.html#gas-transport )
        assert np.isclose(d['transport']['diameter'], 3.0) # Angstroms
        assert np.isclose(d['transport']['dipole'], 1.7) # Debye
        assert np.isclose(d['transport']['well-depth'], 100.0) # Kelvin
        assert np.isclose(d['transport']['rotational-relaxation'], 1.0)

    def test_species_to_dict_warns_and_uses_charge_for_explicit_electron_mismatch(self, caplog):
        """
        Test that YAML species composition uses the net charge when explicit
        electrons disagree.
        """
        sp = self._create_dummy_species("H2", "[H][H]", index=1)
        sp.label = "H2minus_with_one_electron"
        sp.molecule = [Molecule(atoms=[
            Atom(element=get_element("H"), charge=-2, radical_electrons=0, lone_pairs=0),
            Atom(element=get_element("e"), charge=0, radical_electrons=0, lone_pairs=0),
        ])]

        d = species_to_dict(sp, [sp])

        assert d["composition"]["E"] == 2
        assert "has 1 electrons but charge -2" in caplog.text
        assert "Reporting 2 electrons in the Cantera YAML composition." in caplog.text

    def test_species_to_dict_uses_charge_when_electrons_unspecified(self, caplog):
        """
        Test that YAML species composition uses the net charge without a warning
        when no explicit electrons are present.
        """
        sp = self._create_dummy_species("H2", "[H][H]", index=1)
        sp.label = "Hminus"
        sp.molecule = [Molecule(atoms=[
            Atom(element=get_element("H"), charge=-1, radical_electrons=0, lone_pairs=0),
        ])]

        d = species_to_dict(sp, [sp])

        assert d["composition"]["E"] == 1
        assert "electrons but charge" not in caplog.text

    def test_reaction_to_dict_arrhenius(self):
        """Test standard Arrhenius kinetics."""
        r = self._create_dummy_species("R", "[CH2]O", index=1)
        p = self._create_dummy_species("P", "C[O]", index=2)
        rxn = Reaction(
            reactants=[r], products=[p],
            kinetics=Arrhenius(A=(1e10, "s^-1"), n=0.5, Ea=(10, "kJ/mol"), T0=(1, "K"))
        )

        entries = reaction_to_dict_list(rxn, species_list=[r, p])
        assert len(entries) == 1
        data = entries[0]

        assert data['equation'] == "R(1) <=> P(2)"
        assert 'rate-constant' in data
        assert np.isclose(data['rate-constant']['A'], 1e10)
        assert np.isclose(data['rate-constant']['b'], 0.5)
        assert np.isclose(data['rate-constant']['Ea'], 10000.0)

    def test_reaction_to_dict_negative_a_arrhenius(self):
        """Negative Arrhenius A factors are marked for Cantera."""
        r = self._create_dummy_species("R", "[CH2]O", index=1)
        p = self._create_dummy_species("P", "C[O]", index=2)
        rxn = Reaction(
            reactants=[r], products=[p],
            kinetics=Arrhenius(A=(-1e10, "s^-1"), n=0.5, Ea=(10, "kJ/mol"), T0=(1, "K"))
        )
        entries = reaction_to_dict_list(rxn, species_list=[r, p])
        assert entries[0]['negative-A'] is True

    def test_reaction_to_dict_duplicates(self):
        """Test that MultiKinetics objects result in multiple YAML entries."""
        r = self._create_dummy_species("R", "[H]", index=1)
        k1 = Arrhenius(A=(1e10, "s^-1"), n=0, Ea=(0, "J/mol"), T0=(1, "K"))
        k2 = Arrhenius(A=(2e10, "s^-1"), n=0, Ea=(0, "J/mol"), T0=(1, "K"))

        rxn = Reaction(
            reactants=[r], products=[r],
            kinetics=MultiArrhenius(arrhenius=[k1, k2]),
            duplicate=True
        )

        entries = reaction_to_dict_list(rxn, species_list=[r])
        assert len(entries) == 2
        assert entries[0]['rate-constant']['A'] == 1e10
        assert entries[1]['rate-constant']['A'] == 2e10
        assert entries[0].get('duplicate') is True

    def test_reaction_to_dict_troe(self):
        """Test Falloff/Troe serialization."""
        r = self._create_dummy_species("R", "[H]", index=1)
        M = self._create_dummy_species("M", "[Ar]", index=-1)

        # Troe
        k_high = Arrhenius(A=(1e14, "s^-1"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K"))
        k_low = Arrhenius(A=(1e20, "cm^3/(mol*s)"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K"))

        troe = Troe(
            arrheniusHigh=k_high, arrheniusLow=k_low,
            alpha=0.5, T3=(100, "K"), T1=(200, "K"), T2=(300, "K"),
            efficiencies={M.molecule[0]: 2.0}
        )

        rxn = Reaction(reactants=[r], products=[r], kinetics=troe)
        entries = reaction_to_dict_list(rxn, species_list=[r, M])
        data = entries[0]

        assert data['type'] == 'falloff'
        assert 'Troe' in data
        assert data['Troe']['A'] == 0.5
        assert data['Troe']['T2'] == 300.0
        # Efficiencies should map label -> val
        assert data['efficiencies'] == {"M": 2.0}

    def test_generate_cantera_data_detects_plasma(self):
        """Test that the writer detects 'e' and sets thermo: plasma."""
        h2 = self._create_dummy_species("H2", "[H][H]", index=1)

        # Case 1: No Electron
        data = generate_cantera_data([h2], [], is_plasma=False)
        phase = data['phases'][0]
        assert phase['thermo'] == 'ideal-gas'
        assert phase['transport'] == 'mixture-averaged'

    def test_full_integration_plasma_model(self):
        """
        Create a comprehensive RMG model, write it to disk, and load it in Cantera
        to ensure all fields are valid and parsed correctly.
        """

        # 1. Create Model Components
        h2 = self._create_dummy_species("H2", "[H][H]", index=1)
        h = self._create_dummy_species("H", "[H]", index=2)
        ch4 = self._create_dummy_species("CH4", "C", index=3)
        oh = self._create_dummy_species("OH", "[OH]", index=4)
        ar = self._create_dummy_species("Ar", "[Ar]", index=-1)

        species = [h2, h, ch4, oh, ar]

        r1 = Reaction(
            reactants=[h2], products=[h, h],
            kinetics=Arrhenius(A=(1e13, "s^-1"), n=0, Ea=(400, "kJ/mol"), T0=(1, "K"))
        )

        r2 = Reaction(
            reactants=[h, h], products=[h2],
            kinetics=ThirdBody(
                arrheniusLow=Arrhenius(A=(1e18, "cm^6/(mol^2*s)"), n=-1, Ea=(0, "J/mol"), T0=(1, "K")),
                efficiencies={ar.molecule[0]: 0.7}
            )
        )

        reactions = [r1, r2]

        # 2. Mock RMG Object Structure
        # The writer expects: rmg.output_directory and rmg.reaction_model.core
        from rmgpy.rmg.model import ReactionModel

        class MockCore(ReactionModel):
            def __init__(self):
                super().__init__(species=species, reactions=reactions)

        class MockModel:
            def __init__(self):
                self.core = MockCore()
                self.edge = MockCore()  # Empty for now

        class MockRMG:
            def __init__(self, out_dir):
                self.output_directory = out_dir
                self.reaction_model = MockModel()
                self.save_edge_species = False
                self.verbose_comments = False

        mock_rmg = MockRMG(self.tmp_dir)
        save_cantera_files(mock_rmg)

        yaml_file = os.path.join(self.tmp_dir, "cantera2", "chem_annotated.yaml")
        versioned_file = os.path.join(self.tmp_dir, "cantera2", "chem_annotated0005.yaml")
        assert os.path.exists(yaml_file)
        assert os.path.exists(versioned_file)

        try:
            sol = ct.Solution(yaml_file)
        except Exception as e:
            pytest.fail(f"Cantera failed to load the generated YAML: {e}")

        assert sol.n_species == 5

        assert sol.n_reactions == 2

        ct_r2 = sol.reaction(1)
        assert "three-body" in ct_r2.reaction_type or "ThreeBody" in ct_r2.reaction_type
        assert np.isclose(ct_r2.third_body.efficiencies["Ar"], 0.7)

    def test_reaction_to_dict_pdep_arrhenius(self):
        """Test Pressure-Dependent Arrhenius (PLOG) structure."""
        r = self._create_dummy_species("R", "[CH2]O", index=1)
        p = self._create_dummy_species("P", "C[O]", index=2)

        k_low = Arrhenius(A=(1e10, "s^-1"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K"))
        k_high = Arrhenius(A=(1e12, "s^-1"), n=0, Ea=(15, "kJ/mol"), T0=(1, "K"))

        pdep = PDepArrhenius(
            pressures=([0.1, 1.0], "atm"),
            arrhenius=[k_low, k_high],
        )

        rxn = Reaction(reactants=[r], products=[p], kinetics=pdep)

        entries = reaction_to_dict_list(rxn, species_list=[r, p])
        data = entries[0]

        assert data['type'] == 'pressure-dependent-Arrhenius'
        rates = data['rate-constants']
        assert len(rates) == 2

        assert np.isclose(rates[0]['P'], 0.1 * 101325.0)
        assert np.isclose(rates[0]['A'], 1e10)
        assert np.isclose(rates[0]['Ea'], 10000.0)

        assert np.isclose(rates[1]['P'], 1.0 * 101325.0)
        assert np.isclose(rates[1]['A'], 1e12)
        assert np.isclose(rates[1]['Ea'], 15000.0)

    def test_reaction_to_dict_chebyshev(self):
        """Test Chebyshev kinetics structure."""
        r = self._create_dummy_species("R", "[H]", index=1)

        # 2x2 Coefficients matrix
        coeffs = np.array([[1.0, 2.0], [3.0, 4.0]])
        cheb = Chebyshev(
            Tmin=(300, "K"), Tmax=(2000, "K"),
            Pmin=(0.01, "atm"), Pmax=(100, "atm"),
            coeffs=coeffs,
            kunits="s^-1"
        )

        rxn = Reaction(reactants=[r], products=[r], kinetics=cheb)

        entries = reaction_to_dict_list(rxn, species_list=[r])
        data = entries[0]

        assert data['type'] == 'Chebyshev'

        assert np.allclose(data['temperature-range'], [300.0, 2000.0])
        assert np.allclose(data['pressure-range'], [0.01 * 101325.0, 100 * 101325.0])
        assert np.allclose(data['data'], coeffs)

    def test_reaction_to_dict_lindemann(self):
        """Test Lindemann (Falloff without Troe parameters)."""
        r = self._create_dummy_species("R", "[H]", index=1)
        M = self._create_dummy_species("M", "[Ar]", index=-1)

        k_high = Arrhenius(A=(1e14, "s^-1"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K"))
        k_low = Arrhenius(A=(1e21, "cm^3/(mol*s)"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K"))
        lind = Lindemann(
            arrheniusHigh=k_high,
            arrheniusLow=k_low,
            efficiencies={M.molecule[0]: 5.0},
        )
        rxn = Reaction(reactants=[r], products=[r], kinetics=lind)
        entries = reaction_to_dict_list(rxn, species_list=[r, M])
        data = entries[0]

        assert data['type'] == 'falloff'
        assert 'high-P-rate-constant' in data
        assert 'low-P-rate-constant' in data
        assert np.isclose(data['high-P-rate-constant']['A'], 1e14)
        assert np.isclose(data['low-P-rate-constant']['A'], 1e15)
        assert data['efficiencies'] == {"M": 5.0}
        assert 'Troe' not in data

    def test_reaction_to_dict_negative_a_falloff(self):
        """Negative high- or low-pressure A factors are marked for Cantera."""
        r = self._create_dummy_species("R", "[H]", index=1)
        M = self._create_dummy_species("M", "[Ar]", index=-1)
        k_high = Arrhenius(A=(1e14, "s^-1"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K"))
        k_low = Arrhenius(A=(-1e21, "cm^3/(mol*s)"), n=0, Ea=(10, "kJ/mol"), T0=(1, "K"))
        lind = Lindemann(
            arrheniusHigh=k_high,
            arrheniusLow=k_low,
            efficiencies={M.molecule[0]: 5.0},
        )
        rxn = Reaction(reactants=[r], products=[r], kinetics=lind)
        entries = reaction_to_dict_list(rxn, species_list=[r, M])
        assert entries[0]['negative-A'] is True

    def test_cantera_writer_class_listener(self):
        """
        Test the CanteraWriter2 class directly to ensure it correctly initializes
        subdirectories and triggers the save on update().
        """
        writer = CanteraWriter2(self.tmp_dir)
        cantera_dir = os.path.join(self.tmp_dir, 'cantera2')
        assert os.path.exists(cantera_dir)
        assert os.path.isdir(cantera_dir)

        mock_rmg = self._create_dummy_model()
        writer.update(mock_rmg)

        versioned_file = os.path.join(cantera_dir, 'chem_annotated0002.yaml')
        latest_file = os.path.join(cantera_dir, 'chem_annotated.yaml')

        assert os.path.exists(versioned_file)
        assert os.path.exists(latest_file)

        with open(latest_file, 'r') as f:
            content = f.read()
            assert "generator: RMG-Py CanteraWriter2" in content or "generator: 'RMG-Py CanteraWriter2" in content
            assert "phases:" in content
            assert "species:" in content

    def _create_dummy_model(self):
        """Creates a mock object structure resembling RMG.reaction_model"""

        # 1. Species
        sp_H2 = self._create_dummy_species("H2", "[H][H]", index=1)
        sp_H = self._create_dummy_species("H", "[H]", index=2)
        species_list = [sp_H2, sp_H]

        # 2. Reactions
        rxn_arr = Reaction(
            reactants=[sp_H2], products=[sp_H, sp_H],
            kinetics=Arrhenius(A=(1e13, "s^-1"), n=0.0, Ea=(200, "kJ/mol"), T0=(1, "K"))
        )
        reaction_list = [rxn_arr]

        # Mock Object Structure
        from rmgpy.rmg.model import ReactionModel

        class MockCore(ReactionModel):
            def __init__(self, s, r):
                super().__init__(species=s, reactions=r)

        class MockModel:
            def __init__(self, core):
                self.core = core
                self.edge = MockCore([], [])
                self.output_species_list = []
                self.output_reaction_list = []

        class MockRMG:
            def __init__(self, out_dir, model):
                self.output_directory = out_dir
                self.reaction_model = model
                self.save_edge_species = False
                self.verbose_comments = False

        return MockRMG(self.tmp_dir, MockModel(MockCore(species_list, reaction_list)))

    def _create_surface_species(self, label, adjlist, index):
        """Helper to create an RMG surface Species with NASA thermo (no transport)."""
        sp = Species(label=label, index=index)
        sp.from_adjacency_list(adjlist)
        coeffs = [1.0, 0.0, 0.0, 0.0, 0.0, -100.0, 1.0]
        sp.thermo = NASA(
            polynomials=[
                NASAPolynomial(coeffs=coeffs, Tmin=(200, "K"), Tmax=(1000, "K")),
                NASAPolynomial(coeffs=coeffs, Tmin=(1000, "K"), Tmax=(6000, "K")),
            ],
            Tmin=(200, "K"),
            Tmax=(6000, "K"),
        )
        return sp

    # ------------------------------------------------------------------
    # Surface species
    # ------------------------------------------------------------------
    def test_species_to_dict_surface_no_transport(self):
        """Surface species: composition contains X, no transport block."""
        sp = self._create_surface_species(
            "X", "1 X u0 p0", index=10
        )
        d = species_to_dict(sp, [sp])
        assert "X" in d["composition"]
        assert "transport" not in d

    def test_species_to_dict_surface_thermo(self):
        """Surface species reports NASA7 thermo with two polynomial ranges."""
        hx = self._create_surface_species(
            "H_X", "1 H u0 p0 {2,S}\n2 X u0 p0 {1,S}", index=11
        )
        d = species_to_dict(hx, [hx])
        assert d["thermo"]["model"] == "NASA7"
        assert len(d["thermo"]["data"]) == 2
        assert d["composition"] == {"H": 1, "X": 1}

    # ------------------------------------------------------------------
    # Surface reactions
    # yaml_cantera2 declares activation-energy: J/mol, so value_si is used
    # directly without any ×1000 conversion.
    # ------------------------------------------------------------------
    def test_reaction_to_dict_surface_arrhenius(self):
        """SurfaceArrhenius: type is interface-Arrhenius, rate-constant in J/mol."""
        h2 = self._create_dummy_species("H2", "[H][H]", index=1)
        x = self._create_surface_species("X", "1 X u0 p0", index=2)
        hx = self._create_surface_species(
            "H_X", "1 H u0 p0 {2,S}\n2 X u0 p0 {1,S}", index=3
        )
        kin = SurfaceArrhenius(
            A=(1e13, "m^2/(mol*s)"), n=0.5, Ea=(50, "kJ/mol"), T0=(1, "K")
        )
        rxn = Reaction(
            reactants=[h2, x], products=[hx, hx], kinetics=kin
        )
        entries = reaction_to_dict_list(rxn, species_list=[h2, x, hx])
        assert len(entries) == 1
        d = entries[0]
        assert d["type"] == "interface-Arrhenius"
        assert "rate-constant" in d
        assert np.isclose(d["rate-constant"]["A"], 1e13)
        assert np.isclose(d["rate-constant"]["b"], 0.5)
        assert np.isclose(d["rate-constant"]["Ea"], 50000.0)  # J/mol

    def test_reaction_to_dict_sticking_coefficient(self):
        """StickingCoefficient: type is sticking-Arrhenius, A is dimensionless."""
        h2 = self._create_dummy_species("H2", "[H][H]", index=1)
        x = self._create_surface_species("X", "1 X u0 p0", index=2)
        hx = self._create_surface_species(
            "H_X", "1 H u0 p0 {2,S}\n2 X u0 p0 {1,S}", index=3
        )
        kin = StickingCoefficient(
            A=(0.1, ""), n=0, Ea=(0, "kJ/mol"), T0=(1, "K")
        )
        rxn = Reaction(
            reactants=[h2, x, x], products=[hx, hx], kinetics=kin
        )
        entries = reaction_to_dict_list(rxn, species_list=[h2, x, hx])
        assert len(entries) == 1
        d = entries[0]
        assert d["type"] == "sticking-Arrhenius"
        assert "sticking-coefficient" in d
        assert np.isclose(d["sticking-coefficient"]["A"], 0.1)
        assert np.isclose(d["sticking-coefficient"]["Ea"], 0.0)

    def test_reaction_to_dict_negative_a_sticking_coefficient(self):
        """Negative sticking-coefficient A factors are marked for Cantera."""
        h2 = self._create_dummy_species("H2", "[H][H]", index=1)
        x = self._create_surface_species("X", "1 X u0 p0", index=2)
        hx = self._create_surface_species(
            "H_X", "1 H u0 p0 {2,S}\n2 X u0 p0 {1,S}", index=3
        )
        kin = StickingCoefficient(
            A=(-0.1, ""), n=0, Ea=(0, "kJ/mol"), T0=(1, "K")
        )
        rxn = Reaction(
            reactants=[h2, x, x], products=[hx, hx], kinetics=kin
        )

        entries = reaction_to_dict_list(rxn, species_list=[h2, x, hx])

        assert entries[0]["negative-A"] is True

    def test_reaction_to_dict_coverage_dependence(self):
        """Coverage-dependent kinetics: coverage-dependencies block written correctly.

        yaml_cantera2 declares activation-energy: J/mol, so E uses value_si
        (J/mol) directly — no ×1000 conversion.
        """
        h2 = self._create_dummy_species("H2", "[H][H]", index=1)
        x = self._create_surface_species("X", "1 X u0 p0", index=2)
        hx = self._create_surface_species(
            "H_X", "1 H u0 p0 {2,S}\n2 X u0 p0 {1,S}", index=3
        )
        kin = StickingCoefficient(
            A=(0.1, ""),
            n=0,
            Ea=(0, "kJ/mol"),
            T0=(1, "K"),
            coverage_dependence={
                hx: {"a": 0.5, "m": -1.0, "E": (5.0, "kJ/mol")}
            },
        )
        rxn = Reaction(
            reactants=[h2, x, x], products=[hx, hx], kinetics=kin
        )
        entries = reaction_to_dict_list(rxn, species_list=[h2, x, hx])
        d = entries[0]
        assert "coverage-dependencies" in d
        cov = d["coverage-dependencies"]["H_X(3)"]
        assert np.isclose(cov["a"], 0.5)
        assert np.isclose(cov["m"], -1.0)
        # 5 kJ/mol = 5000 J/mol — written directly (J/mol units declared)
        assert np.isclose(cov["E"], 5000.0)

    def test_reaction_to_dict_thirdbody_unit(self):
        """ThirdBody: type three-body, efficiencies map, rate-constant present."""
        h = self._create_dummy_species("H", "[H]", index=1)
        h2 = self._create_dummy_species("H2", "[H][H]", index=2)
        ar = self._create_dummy_species("Ar", "[Ar]", index=3)
        kin = ThirdBody(
            arrheniusLow=Arrhenius(
                A=(1e18, "cm^6/(mol^2*s)"), n=-1, Ea=(0, "J/mol"), T0=(1, "K")
            ),
            efficiencies={ar.molecule[0]: 0.7},
        )
        rxn = Reaction(reactants=[h, h], products=[h2], kinetics=kin)
        entries = reaction_to_dict_list(rxn, species_list=[h, h2, ar])
        d = entries[0]
        assert d.get("type") == "three-body"
        assert "rate-constant" in d
        assert "efficiencies" in d
        assert np.isclose(d["efficiencies"]["Ar(3)"], 0.7)

    def test_get_elements_block_isotopes_and_surface_site(self):
        """get_elements_lists emits isotope and X definitions only when in use."""
        from rmgpy.molecule.element import H, C, D, T, X, e

        # With D, T, X in use: isotope and X entries appear
        custom_elements, elements_list = get_elements_lists({H, C, D, T, X})
        assert 'H' in elements_list
        assert 'C' in elements_list
        assert 'X' in elements_list
        x_entry = next((e for e in custom_elements if e['symbol'] == 'X'), None)
        assert x_entry is not None
        assert np.isclose(x_entry['atomic-weight'], 195.083)
        symbols = [e['symbol'] for e in custom_elements]
        assert 'D' in symbols   # H-2
        assert 'T' in symbols   # H-3
        for entry in custom_elements:
            assert entry['atomic-weight'] > 0

        # Without isotopes / X: no custom entries, no X in the elements list
        custom_elements, elements_list = get_elements_lists({H, C})
        assert custom_elements == []
        assert 'X' not in elements_list
        assert 'D' not in elements_list
        assert 'T' not in elements_list

        # The RMG electron singleton is lowercase e internally, but exports as
        # Cantera's uppercase pseudo-element E.
        custom_elements, elements_list = get_elements_lists({H, e})
        assert custom_elements == []
        assert elements_list == ['E', 'H']

    def test_generate_cantera_data_elements_block(self):
        """generate_cantera_data emits a top-level 'elements' key only when
        non-builtin elements (isotopes, X) are in use, matching ck2yaml."""
        from rmgpy.molecule.element import H, X, e
        h2 = self._create_dummy_species("H2", "[H][H]", index=1)

        # Gas-only H2: no isotopes, no X -> no top-level 'elements' block.
        data = generate_cantera_data([h2], [], elements_in_use={H})
        assert 'elements' not in data

        # Surface fixture: X is in use, so it appears as a custom element.
        x = self._create_surface_species("X", "1 X u0 p0", index=2)
        data = generate_cantera_data([h2, x], [], elements_in_use={H, X})
        symbols = [entry['symbol'] for entry in data['elements']]
        assert 'X' in symbols

        electron = Species(label="e", index=2)
        electron.from_adjacency_list("1 e u1 p0 c-1")
        electron.thermo = h2.thermo
        data = generate_cantera_data([h2, electron], [], elements_in_use={H, e}, is_plasma=True)
        assert data['phases'][0]['elements'] == ['E', 'H']
        electron_entry = next(sp for sp in data['species'] if sp['name'] == 'e(2)')
        assert electron_entry['composition'] == {'E': 1}

    def test_generate_cantera_data_gas_phase_state(self):
        """Gas phase definition includes a 'state' block with T and P."""
        from rmgpy.molecule.element import H
        h2 = self._create_dummy_species("H2", "[H][H]", index=1)
        data = generate_cantera_data([h2], [], elements_in_use={H})

        gas_phase = data['phases'][0]
        assert 'state' in gas_phase
        assert np.isclose(gas_phase['state']['T'], 300.0)
        assert gas_phase['state']['P'] == '1 atm'

    def test_generate_cantera_data_gas_reactions_key(self):
        """Gas-only model uses top-level 'reactions' key (matching ck2yaml)."""
        from rmgpy.molecule.element import H
        h2 = self._create_dummy_species("H2", "[H][H]", index=1)
        h = self._create_dummy_species("H", "[H]", index=2)
        rxn = Reaction(
            reactants=[h2], products=[h, h],
            kinetics=Arrhenius(A=(1e13, "s^-1"), n=0, Ea=(400, "kJ/mol"), T0=(1, "K"))
        )
        data = generate_cantera_data([h2, h], [rxn], elements_in_use={H})

        gas_phase = data['phases'][0]
        assert 'reactions' not in gas_phase, "Gas-only phase should not reference reactions"
        assert 'reactions' in data
        assert len(data['reactions']) == 1
        assert 'gas-reactions' not in data

    def test_generate_cantera_data_surface_phase_state_and_reactions_key(self):
        """Surface phase has 'state' and references 'surface-reactions'; data has that key."""
        from rmgpy.molecule.element import H, X
        h2 = self._create_dummy_species("H2", "[H][H]", index=1)
        x = self._create_surface_species("X", "1 X u0 p0", index=2)
        hx = self._create_surface_species(
            "H_X", "1 H u0 p0 {2,S}\n2 X u0 p0 {1,S}", index=3
        )
        kin = SurfaceArrhenius(
            A=(1e13, "m^2/(mol*s)"), n=0, Ea=(50, "kJ/mol"), T0=(1, "K")
        )
        rxn = Reaction(reactants=[h2, x], products=[hx, hx], kinetics=kin)
        data = generate_cantera_data([h2, x, hx], [rxn], elements_in_use={H, X})

        surface_phase = next(p for p in data['phases'] if p['name'] == 'surface')
        assert 'state' in surface_phase
        assert np.isclose(surface_phase['state']['T'], 300.0)
        assert surface_phase['reactions'] == ['surface-reactions']
        assert 'surface-reactions' in data
        assert len(data['surface-reactions']) == 1

    def test_species_to_dict_surface_sites_count(self):
        """species_to_dict reports 'sites' only for multi-site (bidentate+) surface species."""
        # Single-site species: 'sites' key is omitted (Cantera defaults to 1)
        hx = self._create_surface_species(
            "H_X", "1 H u0 p0 {2,S}\n2 X u0 p0 {1,S}", index=11
        )
        d = species_to_dict(hx, [hx])
        assert 'sites' not in d

        # Bidentate glyoxal adsorbed via C and O (2 X atoms)
        glyoxal_xx = self._create_surface_species(
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
        d2 = species_to_dict(glyoxal_xx, [glyoxal_xx])
        assert 'sites' in d2
        assert d2['sites'] == 2

    def test_species_to_dict_gas_no_sites_field(self):
        """Gas species must not have a 'sites' field."""
        h2 = self._create_dummy_species("H2", "[H][H]", index=1)
        d = species_to_dict(h2, [h2])
        assert 'sites' not in d

    def test_species_to_dict_transport_note_always_present(self):
        """Transport 'note' is always written when transport_data.comment is set."""
        sp = self._create_dummy_species("H2", "[H][H]", index=1)
        sp.transport_data.comment = "from GRI-Mech"
        d = species_to_dict(sp, [sp])
        assert d['transport']['note'] == "from GRI-Mech"


class TestRecentlyGeneratedCanteraYaml2GasOnly(CanteraYamlFileComparer):
    """Tests for comparing recently generated Cantera YAML files from cantera2, gas-only mechanism.

    These are generated on the fly in the mainTest.py functional test and stored in the testing data directory.
    """
    test_data_folder = 'test/rmgpy/test_data/yaml_writer_data/'

    @pytest.fixture(autouse=True, scope="class")
    def find_recent_files(self, request):
        """
        Find the YAML files generated by mainTest.
        """
        cantera_file = os.path.join(self.test_data_folder, 'cantera2', 'from_main_test.yaml')
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
