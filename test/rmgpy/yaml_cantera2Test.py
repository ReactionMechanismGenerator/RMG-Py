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
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.transport import TransportData
from rmgpy.yaml_cantera2 import (
    CanteraWriter2,
    save_cantera_files,
    species_to_dict,
    reaction_to_dict_list,
    generate_cantera_data
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
            dipoleMoment=(0.0, 'De'),
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
        # Diameter should be in meters (SI)
        assert np.isclose(d['transport']['diameter'], 3.0e-10)

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
        class MockCore:
            def __init__(self):
                self.species = species
                self.reactions = reactions

        class MockModel:
            def __init__(self):
                self.core = MockCore()
                self.edge = MockCore()  # Empty for now

        class MockRMG:
            def __init__(self, out_dir):
                self.output_directory = out_dir
                self.reaction_model = MockModel()
                self.save_edge_species = False

        mock_rmg = MockRMG(self.tmp_dir)
        save_cantera_files(mock_rmg)

        yaml_file = os.path.join(self.tmp_dir, "cantera", "chem.yaml")
        versioned_file = os.path.join(self.tmp_dir, "cantera", "chem0005.yaml")
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

    def test_cantera_writer_class_listener(self):
        """
        Test the CanteraWriter2 class directly to ensure it correctly initializes
        subdirectories and triggers the save on update().
        """
        writer = CanteraWriter2(self.tmp_dir)
        cantera_dir = os.path.join(self.tmp_dir, 'cantera')
        assert os.path.exists(cantera_dir)
        assert os.path.isdir(cantera_dir)

        mock_rmg = self._create_dummy_model()
        writer.update(mock_rmg)

        versioned_file = os.path.join(cantera_dir, 'chem0002.yaml')
        latest_file = os.path.join(cantera_dir, 'chem.yaml')

        assert os.path.exists(versioned_file)
        assert os.path.exists(latest_file)

        with open(latest_file, 'r') as f:
            content = f.read()
            assert "generator: RMG-Py CanteraWriter2" in content
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
        class MockCore:
            def __init__(self, s, r):
                self.species = s
                self.reactions = r

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

        return MockRMG(self.tmp_dir, MockModel(MockCore(species_list, reaction_list)))
