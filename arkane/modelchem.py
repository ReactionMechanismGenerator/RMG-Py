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
Model chemistry standardization module
"""

from __future__ import annotations

import os
import yaml
from dataclasses import dataclass, Field, fields, MISSING, replace
from typing import Iterable, Union

from rmgpy import settings
from rmgpy.rmgobject import RMGObject


def standardize_name(name: str) -> str:
    """Remove spaces, hyphens, and make lowercase"""
    if isinstance(name, str):
        return name.replace('-', '').replace(' ', '').lower()
    else:
        return name


with open(os.path.join(settings['database.directory'], 'quantum_corrections', 'lot_constraints.yml')) as f:
    METHODS_THAT_REQUIRE_SOFTWARE = yaml.safe_load(f)['METHODS_THAT_REQUIRE_SOFTWARE']
METHODS_THAT_REQUIRE_SOFTWARE = {standardize_name(method) for method in METHODS_THAT_REQUIRE_SOFTWARE}


@dataclass(frozen=True)
class LOT(RMGObject):
    """
    Base class for quantum chemistry settings.
    """
    def _is_default(self, f: Field) -> bool:
        """Check if field is set to its default value"""
        if f.default is not MISSING:
            if f.default is None:  # Check using 'is'
                return True if getattr(self, f.name) is None else False
            else:
                return True if getattr(self, f.name) == f.default else False
        else:
            return False

    def __repr__(self) -> str:
        """
        Don't include attributes that are set to their default values,
        and don't include spaces.
        """
        r = (self.__class__.__qualname__ + '('
             + ','.join(f"{f.name}={getattr(self, f.name)!r}" for f in fields(self)
                        if f.repr and f.init and not self._is_default(f))
             + ')')
        return r.replace(' ', '')

    def update(self, **kwargs) -> LOT:
        """
        Return a new instance with updated attributes.

        Args:
            kwargs: Keyword arguments specifying attributes and their new values.

        Returns:
            New instance of same type.
        """
        return replace(self, **kwargs)


@dataclass(repr=False, frozen=True)  # repr=False uses parent repr; frozen makes instances quasi-immutable and hashable
class LevelOfTheory(LOT):
    """
    Uniquely defines the settings used for a quantum calculation.

    Attributes:
        method: Quantum chemistry method.
        basis: Basis set.
        auxiliary_basis: Auxiliary basis set for correlated methods.
        cabs: Complementary auxiliary basis set for F12 calculations.
        software: Quantum chemistry software.
        software_version: Quantum chemistry software version.
        solvent: Solvent.
        solvation_method: Solvation method.
        args: Tuple of additional arguments provided to the software.
    """
    method: str
    basis: str = None
    auxiliary_basis: str = None
    cabs: str = None
    software: str = None
    software_version: Union[int, float, str] = None
    solvent: str = None
    solvation_method: str = None
    args: Union[str, Iterable[str]] = None

    def __post_init__(self):
        """
        Standardize attribute values and check if software is set.

        __post_init__ should allow mutating attributes of frozen
        instances, but it doesn't, so use object.__setattr__ to get
        around that issue.
        """
        for attr, val in self.__dict__.items():
            if val is not None:
                if attr == 'software':
                    std_fn = get_software_id
                elif attr == 'args':  # Standardize and sort args to make unique; convert to tuple
                    def std_fn(args):
                        args = (args,) if isinstance(args, str) else args
                        return tuple(sorted(standardize_name(a) for a in args))
                else:
                    std_fn = standardize_name
                object.__setattr__(self, attr, std_fn(val))

        if self.method in METHODS_THAT_REQUIRE_SOFTWARE and self.software is None:
            raise ValueError(f'Software must be set when using {self.method} method')

    def simple(self) -> LevelOfTheory:
        """
        Convert to level of theory containing only method and basis.
        If the method requires the software attribute to be set,
        include it in the simple representation.

        Returns:
            New instance with only method and basis attributes set.
        """
        if self.method in METHODS_THAT_REQUIRE_SOFTWARE:
            return LevelOfTheory(method=self.method, basis=self.basis, software=self.software)
        else:
            return LevelOfTheory(method=self.method, basis=self.basis)

    def to_model_chem(self) -> str:
        """
        Return model chemistry containing method and basis.

        Returns:
            Model chemistry string.
        """
        if self.basis is None:
            return self.method
        else:
            return f'{self.method}/{self.basis}'


@dataclass(repr=False, frozen=True)
class CompositeLevelOfTheory(LOT):
    """
    Uniquely defines the settings used for a combination of frequency
    and single-point energy quantum calculations.

    Notes:
        Assumes that the geometry optimization was performed using the
        same settings as the frequency calculation.
        CBS-QB3 and similar methods are technically composite methods,
        but they should use the LevelOfTheory interface instead because
        they are specified using a single method.

    Attributes:
        freq: Settings for the frequency calculation.
        energy: Settings for the single-point energy calculation.
    """
    freq: LevelOfTheory
    energy: LevelOfTheory

    def simple(self) -> CompositeLevelOfTheory:
        """
        Convert to composite level of theory containing only methods
        and bases.

        Returns:
            New instance with only method and basis attributes set.
        """
        return CompositeLevelOfTheory(freq=self.freq.simple(), energy=self.energy.simple())

    def to_model_chem(self) -> str:
        """
        Return model chemistry containing methods and bases.

        Returns:
            Model chemistry string.
        """
        return f'{self.energy.to_model_chem()}//{self.freq.to_model_chem()}'


def _to_lot_helper(model_chem: str, **kwargs) -> LevelOfTheory:
    """Convert non-composite level of theory"""
    try:
        method, basis = model_chem.split('/')
    except ValueError:  # no explicit basis
        return LevelOfTheory(method=model_chem, **kwargs)
    else:
        return LevelOfTheory(method=method, basis=basis, **kwargs)


def model_chem_to_lot(model_chem: str,
                      freq_settings: dict = None,
                      energy_settings: dict = None,
                      **kwargs) -> Union[LevelOfTheory, CompositeLevelOfTheory]:
    """
    Convert model chemistry to standardized level of theory object.

    Notes:
        If `model_chem` only contains one level of theory, either
        settings dictionary can be used with `freq_settings` being
        preferred if both are provided.
        Additional settings can also be provided with `kwargs` which
        are applied to all levels of theory.

    Args:
        model_chem: Model chemistry string.
        freq_settings: Additional parameters for frequency level of theory.
        energy_settings: Additional parameters for energy level of theory.
        kwargs: Additional parameters applied to all levels of theory.

    Returns:
        LevelOfTheory or CompositeLevelOfTheory instance.
    """
    if freq_settings is None:
        freq_settings = {}
    if energy_settings is None:
        energy_settings = {}

    try:
        energy, freq = model_chem.split('//')
    except ValueError:  # No energy level
        _settings = freq_settings or energy_settings
        return _to_lot_helper(model_chem, **_settings, **kwargs)
    else:
        return CompositeLevelOfTheory(
            freq=_to_lot_helper(freq, **freq_settings, **kwargs),
            energy=_to_lot_helper(energy, **energy_settings, **kwargs)
        )


def str_to_lot(s: str) -> Union[LevelOfTheory, CompositeLevelOfTheory]:
    """
    Convert string to standardized level of theory object.

    Args:
        s: String representation of object.

    Returns:
        LevelOfTheory or CompositeLevelOfTheory instance.
    """
    return eval(s, {'__builtins__': None},
                {'LevelOfTheory': LevelOfTheory, 'CompositeLevelOfTheory': CompositeLevelOfTheory})


# Map from multiple names for a software to a single standard identifier
# No need for hyphens, spaces, or capitalization in the possible names
_valid_software_names = {
    'gaussian': ('gaussian', 'gaussian16', 'gaussian09', 'gaussian03',
                 'gau', 'gau16', 'gau09', 'gau03',
                 'g16', 'g09', 'g03'),
    'qchem': ('qchem',),
    'molpro': ('molpro',),
    'orca': ('orca',),
    'terachem': ('terachem',),
    'mopac': ('mopac',),
    'psi4': ('psi4',)
}
_software_ids = {_name: _id for _id, _names in _valid_software_names.items() for _name in _names}


def get_software_id(name: str) -> str:
    """
    Given a name for a quantum chemistry software, return its standard
    identifier.

    Args:
        name: Commonly used name of software.

    Returns:
        Unique software identifier.
    """
    try:
        return _software_ids[standardize_name(name)]
    except KeyError:
        raise ValueError(f'"{name}" is an invalid quantum chemistry software')
