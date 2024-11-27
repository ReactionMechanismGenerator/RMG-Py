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
This module provides classes for deriving and applying two types of
bond additivity corrections (BACs).

The first type, Petersson-type BACs, are described in:
Petersson et al., J. Chem. Phys. 1998, 109, 10570-10579

The second type, Melius-type BACs, are described in:
Anantharaman and Melius, J. Phys. Chem. A 2005, 109, 1734-1747
"""

import csv
import importlib
import json
import logging
import os
import re
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Sequence, Set, Tuple, Union

import numpy as np
import scipy.optimize as optimize
from scipy.stats import distributions
from sklearn.model_selection import KFold

try:
    import matplotlib.pyplot as plt
except ImportError as e:
    plt = None
    matplotlib_exception = e

from rmgpy.quantity import ScalarQuantity

import arkane.encorr.data as data
from arkane.encorr.data import Molecule, BACDatapoint, BACDataset, extract_dataset, geo_to_mol
from arkane.encorr.reference import ReferenceSpecies, ReferenceDatabase
from arkane.exceptions import BondAdditivityCorrectionError
from arkane.modelchem import LevelOfTheory, CompositeLevelOfTheory


class BACJob:
    """
    A representation of an Arkane BAC job. This job is used to fit and
    save bond additivity corrections.
    """

    def __init__(self,
                 level_of_theory: Union[LevelOfTheory, CompositeLevelOfTheory],
                 bac_type: str = 'p',
                 crossval_n_folds: int = 1,
                 write_to_database: bool = False,
                 overwrite: bool = False,
                 **kwargs):
        """
        Initialize a BACJob instance.

        Args:
            level_of_theory: The level of theory that will be used to get training data from the RMG database.
            bac_type: 'p' for Petersson-style BACs, 'm' for Melius-style BACs.
            crossval_n_folds: Performs k-fold cross-validation.
                              If k does not equal 1, the fitted BACs are not written to RMG database.
                              1 indicates to not do cross-validation.
                              -1 indicates to do leave-one-out cross-validation.
                              Any other positive integer defines the number of folds.
            write_to_database: Save the fitted BACs directly to the RMG database if `crossval_n_folds=1`
                               i.e. BACs are only saved if fit to all training data.
            overwrite: Overwrite BACs in the RMG database if they already exist.
            kwargs: Additional parameters passed to BAC.fit or CrossVal.fit.
        """
        self.level_of_theory = level_of_theory
        self.bac_type = bac_type
        self.crossval_n_folds = crossval_n_folds
        self.write_to_database = write_to_database
        self.overwrite = overwrite
        self.kwargs = kwargs
        if self.crossval_n_folds != 1:
            self.bac = CrossVal(level_of_theory, bac_type=bac_type, n_folds=crossval_n_folds)
        else:
            self.bac = BAC(level_of_theory, bac_type=bac_type)

    def execute(self, output_directory: str = None, plot: bool = False, jobnum: int = 1):
        """
        Execute the BAC job.

        Args:
            output_directory: Save the results in this directory.
            plot: Save plots of results.
            jobnum: Job number.
        """
        logging.info(f'Running BAC job {jobnum}')
        self.bac.fit(**self.kwargs)

        if output_directory is not None and self.crossval_n_folds == 1:
            os.makedirs(output_directory, exist_ok=True)
            self.write_output(output_directory, jobnum=jobnum)

            if plot:
                self.plot(output_directory, jobnum=jobnum)

        if self.write_to_database and self.crossval_n_folds == 1:
            try:
                self.bac.write_to_database(overwrite=self.overwrite)
            except IOError as e:
                logging.warning('Could not write BACs to database. Captured error:')
                logging.warning(str(e))

    def write_output(self, output_directory: str, jobnum: int = 1):
        """
        Save the BACs to the `output.py` file located in
        `output_directory` and save a CSV file of the results.

        Args:
            output_directory: Save the results in this directory.
            jobnum: Job number.
        """
        model_chemistry_formatted = self.level_of_theory.to_model_chem().replace('//', '__').replace('/', '_')
        output_file1 = os.path.join(output_directory, 'output.py')
        output_file2 = os.path.join(output_directory, f'{jobnum}_{model_chemistry_formatted}.csv')
        logging.info(f'Saving results for {self.level_of_theory}...')

        with open(output_file1, 'a') as f:
            stats_before = self.bac.dataset.calculate_stats()
            stats_after = self.bac.dataset.calculate_stats(for_bac_data=True)
            bac_type_str = f'{"Melius" if self.bac.bac_type == "m" else "Petersson"}-type BACs'
            f.write(f'# Job {jobnum}: {bac_type_str}:\n')
            f.write(f'# Training RMSE/MAE before fitting: {stats_before.rmse:.2f}/{stats_before.mae:.2f} kcal/mol\n')
            f.write(f'# Training RMSE/MAE after fitting: {stats_after.rmse:.2f}/{stats_after.mae:.2f} kcal/mol\n')
            f.writelines(self.bac.format_bacs(ci=True))
            f.write('\n')

        with open(output_file2, 'w') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Index',
                'Label',
                'Smiles',
                'InChI',
                'Formula',
                'Multiplicity',
                'Charge',
                'Reference Enthalpy',
                'Calculated Enthalpy',
                'Corrected Enthalpy',
                'Source'
            ])
            for d in self.bac.dataset:
                writer.writerow([
                    d.spc.index,
                    d.spc.label,
                    d.spc.smiles,
                    d.spc.inchi,
                    d.spc.formula,
                    d.spc.multiplicity,
                    d.spc.charge,
                    f'{d.ref_data:.3f}',
                    f'{d.calc_data:.3f}',
                    f'{d.bac_data:.3f}',
                    d.spc.get_preferred_source()
                ])

    def plot(self, output_directory: str, jobnum: int = 1):
        """
        Plot the distribution of errors before and after fitting BACs
        and plot the parameter correlation matrix.

        Args:
            output_directory: Save the plots in this directory.
            jobnum: Job number
        """
        if plt is None:
            raise matplotlib_exception

        model_chemistry_formatted = self.level_of_theory.to_model_chem().replace('//', '__').replace('/', '_')
        if self.crossval_n_folds == 1:
            correlation_path = os.path.join(output_directory, f'{jobnum}_{model_chemistry_formatted}_correlation.pdf')
            self.bac.save_correlation_mat(correlation_path)

        plt.rcParams.update({'font.size': 16})
        fig_path = os.path.join(output_directory, f'{jobnum}_{model_chemistry_formatted}_errors.pdf')

        fig = plt.figure(figsize=(10, 7))
        ax = fig.gca()

        error_before = self.bac.dataset.calc_data - self.bac.dataset.ref_data
        error_after = self.bac.dataset.bac_data - self.bac.dataset.ref_data
        _, _, patches = ax.hist(
            (error_before, error_after),
            bins=50,
            label=('before fitting', 'after fitting'),
            edgecolor='black',
            linewidth=0.5
        )
        ax.set_xlabel('Error (kcal/mol)')
        ax.set_ylabel('Count')

        hatches = ('////', '----')
        for patch_set, hatch in zip(patches, hatches):
            plt.setp(patch_set, hatch=hatch)
        ax.tick_params(bottom=False)
        ax.set_axisbelow(True)
        ax.grid()
        ax.legend()

        fig.savefig(fig_path, bbox_inches='tight', pad_inches=0)


class BAC:
    """
    A class for deriving and applying bond additivity corrections.
    """

    ref_databases = {}
    atom_spins = {
        'H': 0.5, 'C': 1.0, 'N': 1.5, 'O': 1.0, 'F': 0.5,
        'Si': 1.0, 'P': 1.5, 'S': 1.0, 'Cl': 0.5, 'Br': 0.5, 'I': 0.5, 'Li': 0.5,
    }
    exp_coeff = 3.0  # Melius-type parameter (Angstrom^-1)

    def __init__(self, level_of_theory: Union[LevelOfTheory, CompositeLevelOfTheory], bac_type: str = 'p'):
        """
        Initialize a BAC instance.

        There are two implemented BAC types:
            Petersson-type: Petersson et al., J. Chem. Phys. 1998, 109, 10570-10579
            Melius-type: Anantharaman and Melius, J. Phys. Chem. A 2005, 109, 1734-1747

        Args:
            level_of_theory: Level of theory to get preexisting BACs or data from reference database.
            bac_type: Type of BACs to get/fit ('p' for Petersson and 'm' for Melius).
        """
        self._level_of_theory = self._bac_type = None  # Set these first to avoid errors in setters
        self.level_of_theory = level_of_theory
        self.bac_type = bac_type

        # Attributes related to fitting BACs for a given model chemistry
        self.database_key = None  # Dictionary key to access reference database
        self.dataset = None  # Collection of BACDatapoints in BACDataset
        self.correlation = None  # Correlation matrix for BAC parameters
        self.confidence_intervals = None  # 95% confidence intervals for BAC parameters

        # Define attributes for memoization during fitting
        self._reset_memoization()

    @property
    def bac_type(self) -> str:
        return self._bac_type

    @bac_type.setter
    def bac_type(self, val: str):
        """Check validity and update BACs every time the BAC type is changed."""
        if val not in {'m', 'p'}:
            raise BondAdditivityCorrectionError(f'Invalid BAC type: {val}')
        self._bac_type = val
        self._update_bacs()

    @property
    def level_of_theory(self) -> Union[LevelOfTheory, CompositeLevelOfTheory]:
        return self._level_of_theory

    @level_of_theory.setter
    def level_of_theory(self, val: Union[LevelOfTheory, CompositeLevelOfTheory]):
        """Update BACs every time the level of theory is changed."""
        self._level_of_theory = val
        self._update_bacs()

    def _update_bacs(self):
        self.bacs = None
        try:
            if self.bac_type == 'm':
                self.bacs = data.mbac[self.level_of_theory]
            elif self.bac_type == 'p':
                self.bacs = data.pbac[self.level_of_theory]
        except KeyError:
            pass

    @classmethod
    def load_database(cls,
                      paths: Union[str, List[str]] = None,
                      names: Union[str, List[str]] = None,
                      reload: bool = False) -> str:
        """
        Load a reference database.

        Args:
            paths: Paths to database folders.
            names: Names of database folders in RMG database.
            reload: Force reload of database.

        Returns:
            Key to access just loaded database.
        """
        paths = ReferenceDatabase.get_database_paths(paths=paths, names=names)
        key = cls.get_database_key(paths)
        if key not in cls.ref_databases or reload:
            logging.info(f'Loading reference database from {paths}')
            cls.ref_databases[key] = ReferenceDatabase()
            cls.ref_databases[key].load(paths=paths)
        return key

    @staticmethod
    def get_database_key(paths: Union[str, List[str]]) -> Union[str, Tuple[str, ...]]:
        """Get a key to access a stored reference database based on the database paths."""
        if not (isinstance(paths, str) or (isinstance(paths, list) and all(isinstance(p, str) for p in paths))):
            raise ValueError(f'{paths} paths is invalid')
        return tuple(sorted(paths)) if isinstance(paths, list) else paths

    def _reset_memoization(self):
        self._alpha_coeffs = {}
        self._beta_coeffs = {}
        self._gamma_coeffs = {}
        self._k_coeffs = {}

    def get_correction(self,
                       bonds: Dict[str, int] = None,
                       coords: np.ndarray = None,
                       nums: Iterable[int] = None,
                       datapoint: BACDatapoint = None,
                       spc: ReferenceSpecies = None,
                       multiplicity: int = None) -> ScalarQuantity:
        """
        Returns the bond additivity correction.

        There are two bond additivity corrections currently supported.
        Peterson-type corrections can be specified by setting
        `self.bac_type` to 'p'. This will use the `bonds` variable,
        which is a dictionary associating bond types with the number of
        that bond in the molecule.

        The Melius-type BAC is specified with 'm' and utilizes the atom
        coordinates in `coords` and the structure's multiplicity.

        Args:
            bonds: A dictionary of bond types (e.g., 'C=O') with their associated counts.
            coords: A Numpy array of Cartesian molecular coordinates.
            nums: A sequence of atomic numbers.
            datapoint: If not using bonds, coords, nums, use BACDatapoint.
            spc: Alternatively, use ReferenceSpecies.
            multiplicity: The spin multiplicity of the molecule.

        Returns:
            The bond correction to the electronic energy.
        """
        if self.bacs is None:
            bac_type_str = 'Melius' if self.bac_type == 'm' else 'Petersson'
            raise BondAdditivityCorrectionError(
                f'Missing {bac_type_str}-type BAC parameters for {self.level_of_theory}'
            )

        if datapoint is None and spc is not None:
            datapoint = BACDatapoint(spc, level_of_theory=self.level_of_theory)

        if self.bac_type == 'm':
            return self._get_melius_correction(coords=coords, nums=nums, datapoint=datapoint, multiplicity=multiplicity)
        elif self.bac_type == 'p':
            return self._get_petersson_correction(bonds=bonds, datapoint=datapoint)

    def _get_petersson_correction(self, bonds: Dict[str, int] = None, datapoint: BACDatapoint = None) -> ScalarQuantity:
        """
        Given the level of theory and a dictionary of bonds, return the
        total BAC.

        Args:
            bonds: Dictionary of bonds with the following format:
                bonds = {
                    'C-H': C-H_bond_count,
                    'C-C': C-C_bond_count,
                    'C=C': C=C_bond_count,
                    ...
                }
            datapoint: BACDatapoint instead of bonds.

        Returns:
            Petersson-type bond additivity correction.
        """
        if datapoint is not None:
            if bonds is None:
                bonds = datapoint.bonds
            else:
                logging.warning(f'Species {datapoint.spc.label} will not be used because `bonds` was specified')

        # Sum up corrections for all bonds
        bac = 0.0
        for symbol, count in bonds.items():
            if symbol in self.bacs:
                bac += count * self.bacs[symbol]
            else:
                symbol_flipped = ''.join(re.findall(r'[a-zA-Z]+|[^a-zA-Z]+', symbol)[::-1])  # Check reversed symbol
                if symbol_flipped in self.bacs:
                    bac += count * self.bacs[symbol_flipped]
                else:
                    logging.warning(f'Bond correction not applied for unknown bond type {symbol}.')

        return ScalarQuantity(bac, 'kcal/mol')

    def _get_melius_correction(self,
                               coords: np.ndarray = None,
                               nums: Iterable[int] = None,
                               datapoint: BACDatapoint = None,
                               multiplicity: int = None,
                               params: Dict[str, Union[float, Dict[str, float]]] = None) -> ScalarQuantity:
        """
        Given the level of theory, molecular coordinates, atomic numbers,
        and dictionaries of BAC parameters, return the total BAC.

        Notes:
            A molecular correction term other than 0 destroys the size
            consistency of the quantum chemistry method. This correction
            also requires the multiplicity of the molecule.

            The negative of the total correction described in
            Anantharaman and Melius (JPCA 2005) is returned so that it
            can be added to the energy.

        Args:
            coords: Numpy array of Cartesian atomic coordinates.
            nums: Sequence of atomic numbers.
            datapoint: BACDatapoint instead of molecule.
            multiplicity: Multiplicity of the molecule (not necessary if using datapoint).
            params: Optionally provide parameters other than those stored in self.

        Returns:
            Melius-type bond additivity correction.
        """
        if params is None:
            params = self.bacs
        atom_corr = params['atom_corr']
        bond_corr_length = params['bond_corr_length']
        bond_corr_neighbor = params['bond_corr_neighbor']
        mol_corr = params.get('mol_corr', 0.0)

        # Get single-bonded RMG molecule
        mol = None
        if datapoint is not None:
            if nums is None or coords is None:
                mol = datapoint.to_mol(from_geo=True)
                multiplicity = datapoint.spc.multiplicity  # Use species multiplicity instead
            else:
                logging.warning(
                    f'Species {datapoint.spc.label} will not be used because `nums` and `coords` were specified'
                )
        if mol is None:
            mol = geo_to_mol(coords, nums=nums)

        # Molecular correction
        if mol_corr != 0 and multiplicity is None:
            raise BondAdditivityCorrectionError(f'Missing multiplicity for {mol}')
        bac_mol = mol_corr * self._get_mol_coeff(mol, multiplicity=multiplicity)

        # Atomic correction
        bac_atom = sum(count * atom_corr[symbol] for symbol, count in self._get_atom_counts(mol).items())

        # Bond correction
        bac_length = sum(
            coeff * (bond_corr_length[symbol[0]] * bond_corr_length[symbol[1]]) ** 0.5 if isinstance(symbol, tuple)
            else coeff * bond_corr_length[symbol]
            for symbol, coeff in self._get_length_coeffs(mol).items()
        )
        bac_neighbor = sum(count * bond_corr_neighbor[symbol] for
                           symbol, count in self._get_neighbor_coeffs(mol).items())
        bac_bond = bac_length + bac_neighbor

        # Note the minus sign
        return ScalarQuantity(-(bac_mol + bac_atom + bac_bond), 'kcal/mol')

    def _get_atom_counts(self, mol: Molecule) -> Counter:
        """
        Get a counter containing how many atoms of each type are
        present in the molecule.

        Args:
            mol: RMG-Py molecule.

        Returns:
            Counter containing atom counts.
        """
        if hasattr(mol, 'id') and mol.id is not None:
            if mol.id in self._alpha_coeffs:
                return self._alpha_coeffs[mol.id]

        atom_counts = Counter(atom.element.symbol for atom in mol.atoms)

        if hasattr(mol, 'id'):
            self._alpha_coeffs[mol.id] = atom_counts
        return atom_counts

    def _get_length_coeffs(self, mol: Molecule) -> defaultdict:
        """
        Get a dictionary containing the coefficients for the beta
        (bond_corr_length) variables. There is one coefficient per atom
        type and an additional coefficient for each combination of atom
        types.

        Example: If the atoms are H, C, and O, there are (at most)
        coefficients for H, C, O, (C, H), (H, O), and (C, O).

        Args:
            mol: RMG-Py molecule.

        Returns:
            Defaultdict containing beta coefficients.
        """
        if hasattr(mol, 'id') and mol.id is not None:
            if mol.id in self._beta_coeffs:
                return self._beta_coeffs[mol.id]

        coeffs = defaultdict(float)

        for bond in mol.get_all_edges():
            atom1 = bond.atom1
            atom2 = bond.atom2
            symbol1 = atom1.element.symbol
            symbol2 = atom2.element.symbol

            c = np.exp(-self.exp_coeff * np.linalg.norm(atom1.coords - atom2.coords))
            k = symbol1 if symbol1 == symbol2 else tuple(sorted([symbol1, symbol2]))
            coeffs[k] += c

        if hasattr(mol, 'id'):
            self._beta_coeffs[mol.id] = coeffs
        return coeffs

    def _get_neighbor_coeffs(self, mol: Molecule) -> Counter:
        """
        Get a counter containing the coefficients for the gamma
        (bond_corr_neighbor) variables.

        Args:
            mol: RMG-Py molecule.

        Returns:
            Counter containing gamma coefficients.
        """
        if hasattr(mol, 'id') and mol.id is not None:
            if mol.id in self._gamma_coeffs:
                return self._gamma_coeffs[mol.id]

        coeffs = Counter()

        for bond in mol.get_all_edges():
            atom1 = bond.atom1
            atom2 = bond.atom2

            # Atoms adjacent to atom1
            counts1 = Counter(a.element.symbol for a, b in atom1.bonds.items() if b is not bond)
            counts1[atom1.element.symbol] += max(0, len(atom1.bonds) - 1)

            # Atoms adjacent to atom2
            counts2 = Counter(a.element.symbol for a, b in atom2.bonds.items() if b is not bond)
            counts2[atom2.element.symbol] += max(0, len(atom2.bonds) - 1)

            coeffs += counts1 + counts2

        if hasattr(mol, 'id'):
            self._gamma_coeffs[mol.id] = coeffs
        return coeffs

    def _get_mol_coeff(self, mol: Molecule, multiplicity: int = 1) -> float:
        """
        Get the coefficient for the K (mol_corr) variable.

        Args:
            mol: RMG-Py molecule.
            multiplicity: Multiplicity of the molecule.

        Returns:
            K coefficient.
        """
        if hasattr(mol, 'id') and mol.id is not None:
            if mol.id in self._k_coeffs:
                return self._k_coeffs[mol.id]

        spin = 0.5 * (multiplicity - 1)
        coeff = spin - sum(self.atom_spins[atom.element.symbol] for atom in mol.atoms)

        if hasattr(mol, 'id'):
            self._k_coeffs[mol.id] = coeff
        return coeff

    def fit(self,
            weighted: bool = False,
            db_names: Union[str, List[str]] = 'main',
            idxs: Union[Sequence[int], Set[int], int] = None,
            exclude_idxs: Union[Sequence[int], Set[int], int] = None,
            exclude_elements: Union[Sequence[str], Set[str], str] = None,
            charge: Union[Sequence[Union[str, int]], Set[Union[str, int]], str, int] = 'all',
            multiplicity: Union[Sequence[int], Set[int], int, str] = 'all',
            **kwargs):
        """
        Fits bond additivity corrections using calculated and reference
        data available in the RMG database. The resulting BACs stored
        in self.bacs will be based on kcal/mol.

        Args:
            weighted: Perform weighted least squares by balancing training data.
            db_names: Optionally specify database names to train on (defaults to main).
            idxs: Only include reference species with these indices in the training data.
            exclude_idxs: Exclude reference species with these indices from the training data.
            exclude_elements: Molecules with any of the elements in this sequence are excluded from training data.
            charge: Allowable charges for molecules in training data.
            multiplicity: Allowable multiplicities for molecules in training data.
            kwargs: Keyword arguments for fitting Melius-type BACs (see self._fit_melius).
        """
        self._reset_memoization()
        self.database_key = self.load_database(names=db_names)

        self.dataset = extract_dataset(self.ref_databases[self.database_key], self.level_of_theory,
                                       idxs=idxs, exclude_idxs=exclude_idxs,
                                       exclude_elements=exclude_elements, charge=charge, multiplicity=multiplicity)
        if len(self.dataset) == 0:
            raise BondAdditivityCorrectionError(f'No species available for {self.level_of_theory}')

        if weighted:
            self.dataset.compute_weights()

        if self.bac_type == 'm':
            logging.info(f'Fitting Melius-type BACs for {self.level_of_theory}...')
            self._fit_melius(**kwargs)
        elif self.bac_type == 'p':
            logging.info(f'Fitting Petersson-type BACs for {self.level_of_theory}...')
            self._fit_petersson()

        stats_before = self.dataset.calculate_stats()
        stats_after = self.dataset.calculate_stats(for_bac_data=True)
        logging.info(f'Training RMSE/MAE before fitting: {stats_before.rmse:.2f}/{stats_before.mae:.2f} kcal/mol')
        logging.info(f'Training RMSE/MAE after fitting: {stats_after.rmse:.2f}/{stats_after.mae:.2f} kcal/mol')

    def test(self,
             species: List[ReferenceSpecies] = None,
             dataset: BACDataset = None,
             db_names: Union[str, List[str]] = None) -> BACDataset:
        """
        Test on data.

        Note:
            Only one of `species`, `dataset`, or `db_names` can be specified.

        Args:
            species: Species to test on.
            dataset: BACDataset to test on.
            db_names: Database names to test on..

        Returns:
            BACDataset containing the calculated BAC enthalpies in `bac_data`.
        """
        if sum(1 for arg in (species, dataset, db_names) if arg is not None) > 1:
            raise BondAdditivityCorrectionError('Cannot specify several data sources')

        if species is not None:
            dataset = BACDataset([BACDatapoint(spc, level_of_theory=self.level_of_theory) for spc in species])
        elif db_names is not None:
            database_key = self.load_database(names=db_names)
            dataset = extract_dataset(self.ref_databases[database_key], self.level_of_theory)

        if dataset is None or len(dataset) == 0:
            raise BondAdditivityCorrectionError('No data available for evaluation')

        corr = np.array([self.get_correction(datapoint=d).value_si / 4184 for d in dataset])
        dataset.bac_data = dataset.calc_data + corr
        return dataset

    def _fit_petersson(self):
        """
        Fit Petersson-type BACs.
        """
        features = self.dataset.bonds
        feature_keys = list({k for f in features for k in f})
        feature_keys.sort()

        def make_feature_mat(_features: List[Dict[str, int]]) -> np.ndarray:
            _x = np.zeros((len(_features), len(feature_keys)))
            for idx, f in enumerate(_features):
                flist = []
                for k in feature_keys:
                    try:
                        flist.append(f[k])
                    except KeyError:
                        flist.append(0.0)
                _x[idx] = np.array(flist)
            return _x

        # Assume that variance of observations is unity. This is clearly
        # not true because we know the uncertainties but we often care
        # more about less certain molecules.
        x = make_feature_mat(features)
        y = self.dataset.ref_data - self.dataset.calc_data
        weights = np.diag(self.dataset.weight)
        w = np.linalg.solve(x.T @ weights @ x, x.T @ weights @ y)
        ypred = x @ w

        ci, covariance = get_confidence_intervals(x, y, ypred, weights=weights)
        self.confidence_intervals = dict(zip(feature_keys, ci))  # Parameter estimates are w +/- ci
        self.correlation = _covariance_to_correlation(covariance)

        self.dataset.bac_data = self.dataset.calc_data + ypred
        self.bacs = {fk: wi for fk, wi in zip(feature_keys, w)}

    def _fit_melius(self,
                    fit_mol_corr: bool = True,
                    global_opt: bool = True,
                    global_opt_iter: int = 10,
                    lsq_max_nfev: int = 500):
        """
        Fit Melius-type BACs.

        Args:
            fit_mol_corr: Also fit molecular correction term.
            global_opt: Perform a global optimization.
            global_opt_iter: Number of global opt iterations.
            lsq_max_nfev: Maximum function evaluations in least squares optimizer.
        """
        mols = self.dataset.get_mols(from_geo=True)
        for i, mol in enumerate(mols):
            mol.id = i

        all_atom_symbols = list({atom.element.symbol for mol in mols for atom in mol.atoms})
        all_atom_symbols.sort()
        nelements = len(all_atom_symbols)

        # The order of parameters is
        #     atom_corr (alpha)
        #     bond_corr_length (beta)
        #     bond_corr_neighbor (gamma)
        #     optional: mol_corr (k)
        # where atom_corr are the atomic corrections, bond_corr_length are the bondwise corrections
        # due to bond lengths (bounded by 0 below), bond_corr_neighbor are the bondwise corrections
        # due to neighboring atoms, and mol_corr (optional) is a molecular correction.

        # Choose reasonable bounds depending on the parameter
        lim_alpha = (-5.0, 5.0)
        lim_beta = (0.0, 1e4)
        lim_gamma = (-1.0, 1.0)
        lim_k = (-5.0, 5.0)

        wmin = [lim_alpha[0]] * nelements + [lim_beta[0]] * nelements + [lim_gamma[0]] * nelements
        wmax = [lim_alpha[1]] * nelements + [lim_beta[1]] * nelements + [lim_gamma[1]] * nelements
        if fit_mol_corr:
            wmin.append(lim_k[0])
            wmax.append(lim_k[1])

        def get_params(_w: np.ndarray) -> Dict[str, Union[float, Dict[str, float]]]:
            _atom_corr = dict(zip(all_atom_symbols, _w[:nelements]))
            _bond_corr_length = dict(zip(all_atom_symbols, _w[nelements:2*nelements]))
            _bond_corr_neighbor = dict(zip(all_atom_symbols, _w[2*nelements:3*nelements]))
            _mol_corr = _w[3*nelements] if fit_mol_corr else 0.0
            return dict(
                atom_corr=_atom_corr,
                bond_corr_length=_bond_corr_length,
                bond_corr_neighbor=_bond_corr_neighbor,
                mol_corr=_mol_corr
            )

        def get_bac_data(_w: np.ndarray) -> np.ndarray:
            corr = np.array(
                [self._get_melius_correction(datapoint=d, params=get_params(_w)).value_si / 4184 for d in self.dataset]
            )
            return self.dataset.calc_data + corr

        # Construct weight matrix
        weights = np.diag(self.dataset.weight)

        def residuals(_w: np.ndarray) -> Union[float, np.ndarray]:
            """Calculate residuals"""
            bac_data = get_bac_data(_w)
            return np.sqrt(weights) @ (self.dataset.ref_data - bac_data)

        global_opt_iter = global_opt_iter if global_opt else 1
        results = []

        for it in range(global_opt_iter):
            if global_opt:
                logging.info(f'Global opt iteration {it}')

            # Get random initial guess
            w_alpha = np.random.uniform(*lim_alpha, nelements)
            w_beta = np.exp(np.random.uniform(-5, np.log(lim_beta[1]), nelements))
            w_gamma = np.random.uniform(*lim_gamma, nelements)
            w = np.concatenate((w_alpha, w_beta, w_gamma))
            if fit_mol_corr:
                w_k = np.random.uniform(*lim_k, 1)
                w = np.concatenate((w, w_k))

            res = optimize.least_squares(residuals, w, jac='3-point', bounds=(wmin, wmax),
                                         max_nfev=lsq_max_nfev, verbose=1)
            results.append(res)

        res = min(results, key=lambda r: r.cost)
        w = res.x

        self.dataset.bac_data = get_bac_data(w)
        self.bacs = get_params(w)

        # Estimate parameter covariance matrix using Jacobian
        ci, covariance = get_confidence_intervals(res.jac, self.dataset.ref_data, self.dataset.bac_data,
                                                  weights=weights)
        self.confidence_intervals = get_params(ci)
        self.correlation = _covariance_to_correlation(covariance)

    def write_to_database(self, overwrite: bool = False, alternate_path: str = None):
        """
        Write BACs to database.

        Args:
            overwrite: Overwrite existing BACs.
            alternate_path: Write BACs to this path instead.
        """
        if self.bacs is None:
            raise BondAdditivityCorrectionError('No BACs available for writing')

        data_path = data.quantum_corrections_path
        with open(data_path) as f:
            lines = f.readlines()

        bacs_formatted = self.format_bacs(indent=True)

        bac_dict = data.mbac if self.bac_type == 'm' else data.pbac
        keyword = 'mbac' if self.bac_type == 'm' else 'pbac'
        has_entries = bool(data.mbac) if self.bac_type == 'm' else bool(data.pbac)

        # Add new BACs to file without changing existing formatting
        # First: find the BACs dict in the file
        for i, line in enumerate(lines):
            if keyword in line:
                break
        else:
            # 'pbac' and 'mbac' should both be found at `data_path`
            raise RuntimeError(f'Keyword "{keyword} is not found in the data file. '
                               f'Please check the database file at {data_path} and '
                               f'make sure an up-to-date RMG-database branch is used.')

        # Second: Write the BACs block into the BACs dict
        # Does not overwrite comments
        if self.level_of_theory in bac_dict and overwrite:
            del_idx_start = del_idx_end = None
            lot_repr = repr(self.level_of_theory)
            for j, line2 in enumerate(lines[i:]):
                if lot_repr in line2 and 'Composite' not in lot_repr and 'Composite' not in line2:
                    del_idx_start = i + j
                elif lot_repr in line2 and 'Composite' in lot_repr:
                    del_idx_start = i + j

                if del_idx_start is not None and line2.rstrip() == '    },':  # Can't have comment after final brace
                    del_idx_end = i + j + 1
                    if (lines[del_idx_start - 1].lstrip().startswith('#')
                            or lines[del_idx_end + 1].lstrip().startswith('#')):
                        logging.warning('There may be left over comments from previous BACs')
                    lines[del_idx_start:del_idx_end] = bacs_formatted
                    break

            # Check if the entry is successfully inserted to the `lines`
            if del_idx_start is None or del_idx_end is None:
                raise RuntimeError(f'The script cannot identify the corresponding block for the given BACs. '
                                   f'It is possible that the database file at {data_path} is not correctly '
                                   f'formatted. Please check the file.')

        elif self.level_of_theory in bac_dict and not overwrite:
            raise IOError(
                    f'{self.level_of_theory} already exists. Set `overwrite` to True.'
                )
        else:
            # Either empty BACs dict or adding BACs for a new level of theory
            if not has_entries and '}' in lines[i]:  # Empty BACs dict
                lines[i] = f'{keyword} = {{\n'
                lines[(i+1):(i+1)] = ['\n}\n']
            lines[(i+1):(i+1)] = ['\n'] + bacs_formatted

        with open(data_path if alternate_path is None else alternate_path, 'w') as f:
            f.writelines(lines)

        # Reload data to update BAC dictionaries
        if alternate_path is None:
            importlib.reload(data)

    def format_bacs(self, indent: bool = False, ci: bool = False) -> List[str]:
        """
        Obtain a list of nicely formatted BACs suitable for writelines.

        Args:
            indent: Indent each line for printing in database.
            ci: Append confidence intervals.

        Returns:
            Formatted list of BACs.
        """
        bacs_formatted = json.dumps(self.bacs, indent=4).replace('"', "'").split('\n')
        bacs_formatted[0] = f'"{self.level_of_theory!r}": ' + bacs_formatted[0]
        bacs_formatted[-1] += ','
        bacs_formatted = [e + '\n' for e in bacs_formatted]
        if indent:
            bacs_formatted = ['    ' + e for e in bacs_formatted]

        if ci:
            ci_formatted = ['95% Confidence interval half-widths:']
            ci_formatted += json.dumps(self.confidence_intervals, indent=4).replace('"', "'").split('\n')
            ci_formatted = ['# ' + e + '\n' for e in ci_formatted]
            bacs_formatted.extend(ci_formatted)

        return bacs_formatted

    def save_correlation_mat(self, path: str, labels: List[str] = None):
        """
        Save a visual representation of the parameter correlation matrix.

        Args:
            path: Path to save figure to.
            labels: Parameter labels.
        """

        if plt is None:
            raise matplotlib_exception

        if self.correlation is None:
            raise BondAdditivityCorrectionError('Fit BACs before saving correlation matrix!')

        if labels is None:
            if self.bac_type == 'm':
                param_types = list(self.bacs.keys())
                atom_symbols = list(self.bacs[param_types[0]])
                labels = [r'$\alpha_{' + s + r'}$' for s in atom_symbols]      # atom_corr is alpha
                labels.extend(r'$\beta_{' + s + r'}$' for s in atom_symbols)   # bond_corr_length is beta
                labels.extend(r'$\gamma_{' + s + r'}$' for s in atom_symbols)  # bond_corr_neighbor is gamma
                if len(self.correlation) == 3 * len(atom_symbols) + 1:
                    labels.append('K')  # mol_corr is K
            elif self.bac_type == 'p':
                labels = list(self.bacs.keys())

        fig, ax = plt.subplots(figsize=(11, 11) if self.bac_type == 'm' else (18, 18))
        ax.matshow(self.correlation, cmap=plt.cm.PiYG)

        # Superimpose values as text
        for i in range(len(self.correlation)):
            for j in range(len(self.correlation)):
                c = self.correlation[j, i]
                ax.text(i, j, f'{c: .2f}', va='center', ha='center', fontsize=8)

        # Save lims because they get changed when modifying labels
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        ax.set_xticks(list(range(len(self.correlation))))
        ax.set_yticks(list(range(len(self.correlation))))
        ax.set_xticklabels(labels, fontsize=14, rotation='vertical' if self.bac_type == 'p' else None)
        ax.set_yticklabels(labels, fontsize=14)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.tick_params(bottom=False, top=False, left=False, right=False)

        fig.savefig(path, dpi=600, bbox_inches='tight', pad_inches=0)


class CrossVal:
    """
    A class for BAC fitting with cross-validation.
    """

    def __init__(self, level_of_theory: Union[LevelOfTheory, CompositeLevelOfTheory],
                 bac_type: str = 'p',
                 n_folds: int = -1):
        """
        Initialize a CrossVal instance.

        Args:
            level_of_theory: Level of theory for getting data from reference database.
            bac_type: Type of BACs to fit ('p' for Petersson and 'm' for Melius).
            n_folds: Number of folds to use during cross-validation.
                     Default value is -1 i.e. use leave-one-out cross-validation.
        """
        self.level_of_theory = level_of_theory
        self.bac_type = bac_type
        self.n_folds = n_folds

        self.dataset = None  # Complete dataset containing cross-validation estimates for each data point
        self.bacs = None  # List of BAC instances, one for each fold


    def fit(self,
            db_names: Union[str, List[str]] = 'main',
            idxs: Union[Sequence[int], Set[int], int] = None,
            exclude_idxs: Union[Sequence[int], Set[int], int] = None,
            exclude_elements: Union[Sequence[str], Set[str], str] = None,
            charge: Union[Sequence[Union[str, int]], Set[Union[str, int]], str, int] = 'all',
            multiplicity: Union[Sequence[int], Set[int], int, str] = 'all',
            **kwargs):
        """
        Run cross-validation.

        Args:
            db_names: Optionally specify database names to train on (defaults to main).
            idxs: Only include reference species with these indices in the training data.
            exclude_idxs: Exclude reference species with these indices from the training data.
            exclude_elements: Molecules with any of the elements in this sequence are excluded from training data.
            charge: Allowable charges for molecules in training data.
            multiplicity: Allowable multiplicities for molecules in training data.
            kwargs: Parameters passed to BAC.fit.
        """
        database_key = BAC.load_database(names=db_names)
        self.dataset = extract_dataset(BAC.ref_databases[database_key], self.level_of_theory,
                                       idxs=idxs, exclude_idxs=exclude_idxs,
                                       exclude_elements=exclude_elements, charge=charge, multiplicity=multiplicity)
        self.bacs = []
        test_data_results = []
        if self.n_folds == -1:
            logging.info(f'Starting leave-one-out cross-validation for {self.level_of_theory}')
            folds = KFold(n_splits=len(self.dataset)).split(self.dataset)
        else:
            logging.info(f'Starting {self.n_folds}-fold cross-validation for {self.level_of_theory}')
            folds = KFold(n_splits=self.n_folds).split(self.dataset)
        for i, (train, test) in enumerate(folds):
            logging.info(f'\nFold {i}')
            train_idxs = [self.dataset[i].spc.index for i in train]
            test_data = BACDataset([self.dataset[i] for i in test])
            logging.info(f'Testing on species {", ".join(str(d.spc.index) for d in test_data)}')
            bac = BAC(self.level_of_theory, self.bac_type)
            bac.fit(db_names=db_names, idxs=train_idxs, **kwargs)
            bac.test(dataset=test_data)  # Stores predictions in each BACDataset
            self.bacs.append(bac)
            test_data_results.append(test_data)

            stats_before = test_data.calculate_stats()
            stats_after = test_data.calculate_stats(for_bac_data=True)
            logging.info('Testing results:')
            logging.info(f'RMSE/MAE before fitting: {stats_before.rmse:.2f}/{stats_before.mae:.2f} kcal/mol')
            logging.info(f'RMSE/MAE after fitting: {stats_after.rmse:.2f}/{stats_after.mae:.2f} kcal/mol')

        num_test_data = sum(len(test_data) for test_data in test_data_results)
        rmse_before = np.sqrt(np.sum([test_data.calculate_stats().rmse**2 * len(test_data) for test_data in test_data_results]) / num_test_data)
        mae_before = np.sum([test_data.calculate_stats().mae * len(test_data) for test_data in test_data_results]) / num_test_data
        rmse_after = np.sqrt(np.sum([test_data.calculate_stats(for_bac_data=True).rmse**2 * len(test_data) for test_data in test_data_results]) / num_test_data)
        mae_after = np.sum([test_data.calculate_stats(for_bac_data=True).mae * len(test_data) for test_data in test_data_results]) / num_test_data

        logging.info('\nCross-validation results:')
        logging.info(f'Testing RMSE before fitting: '
                    f'{rmse_before:.2f} kcal/mol')
        logging.info(f'Testing MAE before fitting: '
                    f'{mae_before:.2f} kcal/mol')
        logging.info(f'Testing RMSE after fitting: '
                    f'{rmse_after:.2f} kcal/mol')
        logging.info(f'Testing MAE after fitting: '
                    f'{mae_after:.2f} kcal/mol')


def get_confidence_intervals(x: np.ndarray,
                             y: np.ndarray,
                             ypred: np.ndarray,
                             weights: np.ndarray = None,
                             alpha: float = 0.05):
    """
    Compute confidence intervals with two-sided t-test.

    Args:
        x: Feature matrix.
        y: Target vector.
        ypred: Vector of predictions.
        weights: Weight matrix.
        alpha: Significance level (e.g., alpha=0.05 are 95% confidence intervals).

    Returns:
        Vector of confidence interval half-widths
        and variance-covariance matrix.
    """
    n = len(y)  # Ndata
    p = len(x.T)  # Nparam
    if weights is None:
        weights = np.eye(n)

    e = y - ypred  # Residuals
    sigma2 = e.T @ weights @ e / (n - p)  # MSE
    cov = sigma2 * np.linalg.inv(x.T @ weights @ x)  # covariance matrix
    se = np.sqrt(np.diag(cov))  # standard error
    tdist = distributions.t.ppf(1 - alpha / 2, n - p)  # student-t
    ci = tdist * se  # confidence interval half-width
    return ci, cov


def _covariance_to_correlation(cov: np.ndarray) -> np.ndarray:
    """Convert (unscaled) covariance matrix to correlation matrix"""
    v = np.sqrt(np.diag(cov))
    corr = cov / np.outer(v, v)
    corr[cov == 0] = 0
    return corr
