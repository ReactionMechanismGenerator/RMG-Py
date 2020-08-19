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
This module provides methods for applying energy, frequency scale factor, and bond additivity corrections.
"""

import logging
from typing import Dict, Iterable, Union

import numpy as np

import rmgpy.constants as constants

import arkane.encorr.data as data
from arkane.encorr.bac import BAC
from arkane.exceptions import AtomEnergyCorrectionError, BondAdditivityCorrectionError
from arkane.modelchem import LevelOfTheory, CompositeLevelOfTheory

################################################################################


def get_energy_correction(level_of_theory: Union[LevelOfTheory, CompositeLevelOfTheory],
                          atoms: Dict[str, int],
                          bonds: Dict[str, int],
                          coords: np.ndarray,
                          nums: Iterable[int],
                          multiplicity: int = 1,
                          atom_energies: Dict[str, float] = None,
                          apply_atom_corrections: bool = True,
                          apply_bac: bool = False,
                          bac_type: str = 'p') -> float:
    """
    Calculate a correction to the electronic energy obtained from a
    quantum chemistry calculation at a given level of theory such that
    it is consistent with the normal gas-phase reference states.
    Optionally, correct the energy using bond additivity corrections.

    Args:
        level_of_theory: The level of theory.
        atoms: A dictionary of element symbols with their associated counts.
        bonds: A dictionary of bond types (e.g., 'C=O') with their associated counts.
        coords: A Numpy array of Cartesian molecular coordinates.
        nums: A sequence of atomic numbers.
        multiplicity: The spin multiplicity of the molecule.
        atom_energies: A dictionary of element symbols with their associated atomic energies in Hartree.
        apply_atom_corrections: Include the atom correction to the electronic energy.
        apply_bac: Include the bond additivity correction to the electronic energy.
        bac_type: The type of bond additivity correction to use.

    Returns:
        The correction to the electronic energy in J/mol.
    """
    logging.warning('get_energy_correction has be deprecated, use get_atom_correction '
                    'and get_bac instead')

    corr = 0.0
    if apply_atom_corrections:
        corr += get_atom_correction(level_of_theory, atoms, atom_energies=atom_energies)
    if apply_bac:
        corr += get_bac(level_of_theory, bonds, coords, nums, bac_type=bac_type, multiplicity=multiplicity)

    return corr


def get_atom_correction(level_of_theory: Union[LevelOfTheory, CompositeLevelOfTheory],
                        atoms: Dict[str, int],
                        atom_energies: Dict[str, float] = None) -> float:
    """
    Calculate a correction to the electronic energy obtained from a
    quantum chemistry calculation at a given level of theory such that
    it is consistent with the normal gas-phase reference states.

    Args:
        level_of_theory: The level of theory.
        atoms: A dictionary of element symbols with their associated counts.
        atom_energies: A dictionary of element symbols with their associated atomic energies in Hartree.

    Returns:
        The atom correction to the electronic energy in J/mol.

    The assumption for the multiplicity of each atom is:
    H doublet, C triplet, N quartet, O triplet, F doublet, Si triplet,
    P quartet, S triplet, Cl doublet, Br doublet, I doublet.
    """
    corr = 0.0
    # Step 1: Reference all energies to a level of theory-independent
    # basis by subtracting out that level of theory's atomic energies
    if atom_energies is None:
        energy_level = getattr(level_of_theory, 'energy', level_of_theory)
        try:
            atom_energies = data.atom_energies[energy_level]
        except KeyError:
            try:
                atom_energies = data.atom_energies[energy_level.simple()]
            except KeyError:
                return corr
            else:
                logging.warning(f'No exact atom energy match found for {energy_level}.'
                                f' Using {energy_level.simple()} instead.')

    for symbol, count in atoms.items():
        if symbol in atom_energies:
            corr -= count * atom_energies[symbol] * 4.35974394e-18 * constants.Na  # Convert Hartree to J/mol
        else:
            raise AtomEnergyCorrectionError(
                f'An energy correction for element "{symbol}" is unavailable for {level_of_theory}.'
                ' Turn off atom corrections if only running a kinetics jobs'
                ' or supply a dictionary of atom energies'
                ' as `atomEnergies` in the input file.'
            )

    # Step 2: Atom energy corrections to reach gas-phase reference state
    atom_enthalpy_corrections = {symbol: data.atom_hf[symbol] - data.atom_thermal[symbol] for symbol in data.atom_hf}
    for symbol, count in atoms.items():
        if symbol in atom_enthalpy_corrections:
            corr += count * atom_enthalpy_corrections[symbol] * 4184.0  # Convert kcal/mol to J/mol
        else:
            raise AtomEnergyCorrectionError(
                f'Element "{symbol}" is not yet supported in Arkane.'
                ' To include it, add its experimental heat of formation in the atom_hf'
                ' and atom_thermal dictionaries in input/quantum_corrections/data.py in the RMG database.'
            )

    return corr


def get_bac(level_of_theory: Union[LevelOfTheory, CompositeLevelOfTheory],
            bonds: Dict[str, int],
            coords: np.ndarray,
            nums: Iterable[int],
            bac_type: str = 'p',
            multiplicity: int = 1) -> float:
    """
    Returns the bond additivity correction in J/mol.

    There are two bond additivity corrections currently supported. Peterson-type
    corrections can be specified by setting `bac_type` to 'p'. This will use the
    `bonds` attribute, which is a dictionary associating bond types with the number
    of that bond in the molecule.

    The Melius-type BAC is specified with 'm' and utilizes the atom xyz coordinates
    in `coords` and array of atomic numbers of atoms as well as the structure's multiplicity.

    Args:
        level_of_theory: The level of theory.
        bonds: A dictionary of bond types (e.g., 'C=O') with their associated counts.
        coords: A Numpy array of Cartesian molecular coordinates.
        nums: A sequence of atomic numbers.
        multiplicity: The spin multiplicity of the molecule.
        bac_type: The type of bond additivity correction to use.

    Returns:
        The bond correction to the electronic energy in J/mol.
    """
    def _get_bac(_lot):
        """Helper function to get BACs"""
        bac = BAC(_lot, bac_type=bac_type)
        return bac.get_correction(bonds=bonds, coords=coords, nums=nums, multiplicity=multiplicity).value_si

    # Try to match each of these levels of theory, but issue warning if full level of theory cannot be matched
    lots_to_attempt = [
        level_of_theory,  # Full level of theory
        level_of_theory.simple()  # Only method and basis
    ]
    if isinstance(level_of_theory, CompositeLevelOfTheory):
        lots_to_attempt.extend([
            level_of_theory.energy,  # Full energy level
            level_of_theory.energy.simple()  # Energy level with only method and basis
        ])
    for lot in lots_to_attempt:
        try:
            corr = _get_bac(lot)
        except BondAdditivityCorrectionError as e:
            if lot is not lots_to_attempt[-1]:
                continue
            else:
                if 'BAC parameters' in str(e):
                    bac_type_str = 'Melius' if bac_type == 'm' else 'Petersson'
                    raise BondAdditivityCorrectionError(
                        f'Missing {bac_type_str}-type BAC parameters for {level_of_theory}'
                    )
                else:
                    raise
        else:
            if lot is not lots_to_attempt[0]:
                logging.warning(f'No exact BAC match found for {level_of_theory}. Using {lot} instead.')
            return corr


def assign_frequency_scale_factor(level_of_theory: Union[LevelOfTheory, CompositeLevelOfTheory]) -> Union[int, float]:
    """
    Assign the frequency scaling factor according to the level of theory.
    Refer to https://comp.chem.umn.edu/freqscale/index.html for future updates of these factors

    Sources:
        [1] I.M. Alecu, J. Zheng, Y. Zhao, D.G. Truhlar, J. Chem. Theory Comput. 2010, 6, 2872, DOI: 10.1021/ct100326h
        [2] http://cccbdb.nist.gov/vibscalejust.asp
        [3] http://comp.chem.umn.edu/freqscale/190107_Database_of_Freq_Scale_Factors_v4.pdf
        [4] Calculated as described in 10.1021/ct100326h
        [5] J.A. Montgomery, M.J. Frisch, J. Chem. Phys. 1999, 110, 2822–2827, DOI: 10.1063/1.477924

    Args:
        level_of_theory (LevelOfTheory, CompositeLevelOfTheory): The level of theory.

    Returns:
        float: The frequency scaling factor (1 by default).
    """
    if level_of_theory is None:
        return 1
    freq_level = getattr(level_of_theory, 'freq', level_of_theory)
    try:
        scaling_factor = data.freq_dict[freq_level]
    except KeyError:
        try:
            scaling_factor = data.freq_dict[freq_level.simple()]
        except KeyError:
            scaling_factor = 1
        else:
            logging.warning(f'No exact frequency scaling factor match found for {freq_level}.'
                            f' Using {freq_level.simple()} instead.')
    if scaling_factor == 1:
        logging.warning(
            f'No frequency scaling factor found for {freq_level}. Assuming a value of unity.'
            ' This will affect the partition function and all quantities derived from it '
            ' (thermo quantities and rate coefficients).')
    else:
        logging.info(
            f'Assigned a frequency scale factor of {scaling_factor} for {level_of_theory}'
        )
    return scaling_factor
