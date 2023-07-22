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

"""
This script contains unit tests for the :mod:`arkane.modelchem` module.
"""


from dataclasses import FrozenInstanceError

from arkane.modelchem import (
    LOT,
    LevelOfTheory,
    CompositeLevelOfTheory,
    model_chem_to_lot,
    str_to_lot,
    get_software_id,
)
import pytest

# Instances for use in tests
FREQ = LevelOfTheory(method="wB97X-D", basis="def2-TZVP", software="Gaussian 16", args="very-tight")
ENERGY = LevelOfTheory(method="DLPNO-CCSD(T)-F12", basis="def2-TZVP", software="Orca")
COMPOSITE = CompositeLevelOfTheory(freq=FREQ, energy=ENERGY)

# Representations corresponding to instances
FREQ_REPR = "LevelOfTheory(method='wb97xd',basis='def2tzvp',software='gaussian',args=('verytight',))"
ENERGY_REPR = "LevelOfTheory(method='dlpnoccsd(t)f12',basis='def2tzvp',software='orca')"
COMPOSITE_REPR = f"CompositeLevelOfTheory(freq={FREQ_REPR},energy={ENERGY_REPR})"

# Dictionaries corresponding to instances
FREQ_DICT = {
    "class": "LevelOfTheory",
    "method": "wb97xd",
    "basis": "def2tzvp",
    "software": "gaussian",
    "args": ["verytight"],  # This is a list instead of tuple because that's what YAML files expect
}
ENERGY_DICT = {
    "class": "LevelOfTheory",
    "method": "dlpnoccsd(t)f12",
    "basis": "def2tzvp",
    "software": "orca",
}
COMPOSITE_DICT = {
    "class": "CompositeLevelOfTheory",
    "freq": FREQ_DICT,
    "energy": ENERGY_DICT,
}

# Model chemistries corresponding to instances
FREQ_MODELCHEM = "wb97xd/def2tzvp"
ENERGY_MODELCHEM = "dlpnoccsd(t)f12/def2tzvp"
COMPOSITE_MODELCHEM = f"{ENERGY_MODELCHEM}//{FREQ_MODELCHEM}"


class TestLevelOfTheory:
    """
    A class for testing that the LevelOfTheory class functions properly.
    """

    def test_attrs(self):
        """
        Test that instance behaves correctly.
        """
        assert FREQ.method == "wb97xd"
        assert FREQ.basis == "def2tzvp"
        assert FREQ.software == "gaussian"
        assert FREQ.args == ("verytight",)
        with pytest.raises(FrozenInstanceError):
            FREQ.method = ""

        assert repr(FREQ) == FREQ_REPR
        assert repr(ENERGY) == ENERGY_REPR

        with pytest.raises(ValueError):
            _ = LevelOfTheory(method=FREQ.method)
        lot = LevelOfTheory(method=FREQ.method, software=FREQ.software)
        assert lot.basis is None
        assert lot.auxiliary_basis is None
        assert lot.cabs is None
        assert lot.software_version is None
        assert lot.solvent is None
        assert lot.solvation_method is None
        assert lot.args is None

        assert isinstance(FREQ, LOT)

    def test_comparison(self):
        """
        Test comparisons between instances.
        """
        assert isinstance(hash(FREQ), int)
        assert FREQ != ENERGY
        with pytest.raises(TypeError):
            _ = ENERGY > FREQ

        # Test args in different order
        lot1 = LevelOfTheory("method", args=("arg1", "arg2"))
        lot2 = LevelOfTheory("method", args=("arg2", "arg1"))
        assert lot1 == lot2

    def test_simple(self):
        """
        Test that simple level of theory can be obtained.
        """
        lot = FREQ.simple()
        assert lot is not FREQ
        assert lot.method == FREQ.method
        assert lot.basis == FREQ.basis
        assert lot.software == FREQ.software
        for attr, val in lot.__dict__.items():
            if attr not in {"method", "basis", "software"}:
                assert val is None

    def test_to_model_chem(self):
        """
        Test conversion to model chemistry.
        """
        assert FREQ.to_model_chem() == FREQ_MODELCHEM
        assert ENERGY.to_model_chem() == ENERGY_MODELCHEM

        lot = LevelOfTheory(method="CBS-QB3", software="g16")
        assert lot.to_model_chem() == "cbsqb3"

    def test_update(self):
        """
        Test updating attributes.
        """
        lot = FREQ.update(software="Q-Chem")
        assert lot is not FREQ
        assert lot.software == "qchem"
        with pytest.raises(TypeError):
            FREQ.update(test="test")

    def test_as_dict(self):
        """
        Test conversion to dictionary.
        """
        assert FREQ.as_dict() == FREQ_DICT
        assert ENERGY.as_dict() == ENERGY_DICT


class TestCompositeLevelOfTheory:
    """
    A class for testing that the CompositeLevelOfTheory class functions properly.
    """

    def test_attrs(self):
        """
        Test that instance behaves correctly.
        """
        assert COMPOSITE.freq is FREQ
        assert COMPOSITE.energy is ENERGY
        assert repr(COMPOSITE) == COMPOSITE_REPR
        with pytest.raises(FrozenInstanceError):
            COMPOSITE.energy = ""

        assert isinstance(COMPOSITE, LOT)

    def test_comparison(self):
        """
        Test comparisons between instances.
        """
        other = CompositeLevelOfTheory(freq=ENERGY, energy=FREQ)
        assert isinstance(hash(COMPOSITE), int)
        assert COMPOSITE != other
        with pytest.raises(TypeError):
            _ = COMPOSITE > other

    def test_simple(self):
        """
        Test that simple level of theory can be obtained.
        """
        lot = COMPOSITE.simple()
        assert lot is not COMPOSITE
        assert lot.freq.method == COMPOSITE.freq.method
        assert lot.freq.basis == COMPOSITE.freq.basis
        assert lot.freq.software == COMPOSITE.freq.software
        assert lot.energy.method == COMPOSITE.energy.method
        assert lot.energy.basis == COMPOSITE.energy.basis
        for attr, val in lot.freq.__dict__.items():
            if attr not in {"method", "basis", "software"}:
                assert val is None
        for attr, val in lot.energy.__dict__.items():
            if attr not in {"method", "basis"}:
                assert val is None

    def test_to_model_chem(self):
        """
        Test conversion to model chemistry.
        """
        assert COMPOSITE.to_model_chem() == COMPOSITE_MODELCHEM

    def test_as_dict(self):
        """
        Test conversion to dictionary.
        """
        assert COMPOSITE.as_dict() == COMPOSITE_DICT


class TestFuncs:
    """
    A class for testing that the functions in the modelchem module work.
    """

    def test_model_chem_to_lot(self):
        """
        Test model chemistry to quantum calculation settings conversion.
        """
        assert model_chem_to_lot(FREQ_MODELCHEM, software="gaussian", args="verytight") == FREQ
        assert (
            model_chem_to_lot(
                FREQ_MODELCHEM,
                freq_settings={"software": "gaussian", "args": "verytight"},
            )
            == FREQ
        )
        assert (
            model_chem_to_lot(
                FREQ_MODELCHEM,
                freq_settings={"software": "gaussian", "args": "verytight"},
                energy_settings={"unused setting": None},
            )
            == FREQ
        )
        assert model_chem_to_lot(ENERGY_MODELCHEM, energy_settings={"software": "orca"}) == ENERGY
        assert (
            model_chem_to_lot(
                COMPOSITE_MODELCHEM,
                freq_settings={"software": "gaussian", "args": "verytight"},
                energy_settings={"software": "orca"},
            )
            == COMPOSITE
        )

    def test_str_to_lot(self):
        """
        Test key to quantum calculation settings conversion.
        """
        assert str_to_lot(FREQ_REPR) == FREQ
        assert str_to_lot(ENERGY_REPR) == ENERGY
        assert str_to_lot(COMPOSITE_REPR) == COMPOSITE

    def test_get_software_id(self):
        """
        Test standardized software identifiers.
        """
        test_names = ["gaussian", "Gaussian 09", "g-16", "Gau  03"]
        for name in test_names:
            assert get_software_id(name) == "gaussian"

        test_names = ["qchem", "QChem", "Q-Chem"]
        for name in test_names:
            assert get_software_id(name) == "qchem"

        test_names = ["molpro", "Molpro", "MOLPRO"]
        for name in test_names:
            assert get_software_id(name) == "molpro"

        test_names = ["orca", "Orca", "ORCA"]
        for name in test_names:
            assert get_software_id(name) == "orca"

        test_names = ["terachem", "Terachem", "TeraChem", "Tera-Chem", "Tera Chem"]
        for name in test_names:
            assert get_software_id(name) == "terachem"

        with pytest.raises(ValueError):
            get_software_id("g")
