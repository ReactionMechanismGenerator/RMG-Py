#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
This module provides a class for deriving and applying two types of
bond additivity corrections (BACs).

The first type, Petersson-type BACs, are described in:
Petersson et al., J. Chem. Phys. 1998, 109, 10570-10579

The second type, Melius-type BACs, are described in:
Anantharaman and Melius, J. Phys. Chem. A 2005, 109, 1734-1747
"""

import importlib
import json
import logging
import re
from typing import Dict, Iterable, List, Tuple, Union

import numpy as np
import scipy.optimize as optimize

import arkane.encorr.data as data
from arkane.encorr.data import BACDatapoint, extract_dataset, geo_to_mol
from arkane.encorr.reference import ReferenceSpecies, ReferenceDatabase
from arkane.exceptions import BondAdditivityCorrectionError


class BAC:
    """
    A class for deriving and applying bond additivity corrections.
    """

    ref_databases = {}
    atom_spins = {
        'H': 0.5, 'C': 1.0, 'N': 1.5, 'O': 1.0, 'F': 0.5,
        'Si': 1.0, 'P': 1.5, 'S': 1.0, 'Cl': 0.5, 'Br': 0.5, 'I': 0.5
    }
    exp_coeff = 3.0  # Melius-type parameter (Angstrom^-1)

    def __init__(self, model_chemistry: str, bac_type: str = 'p'):
        """
        Initialize a BAC instance.

        There are two implemented BAC types:
            Petersson-type: Petersson et al., J. Chem. Phys. 1998, 109, 10570-10579
            Melius-type: Anantharaman and Melius, J. Phys. Chem. A 2005, 109, 1734-1747

        Args:
            model_chemistry: Model chemistry to get preexisting BACs or data from reference database.
            bac_type: Type of BACs to get/fit ('p' for Petersson and 'm' for Melius).
        """
        self._model_chemistry = self._bac_type = None  # Set these first to avoid errors in setters
        self.model_chemistry = model_chemistry
        self.bac_type = bac_type

        # Attributes related to fitting BACs for a given model chemistry
        self.database_key = None  # Dictionary key to access reference database
        self.dataset = None  # Collection of BACDatapoints in BACDataset

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
    def model_chemistry(self) -> str:
        return self._model_chemistry

    @model_chemistry.setter
    def model_chemistry(self, val: str):
        """Update BACs every time the model chemistry is changed."""
        self._model_chemistry = val
        self._update_bacs()

    def _update_bacs(self):
        self.bacs = None
        try:
            if self.bac_type == 'm':
                self.bacs = data.mbac[self.model_chemistry]
            elif self.bac_type == 'p':
                self.bacs = data.pbac[self.model_chemistry]
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

    def get_correction(self,
                       bonds: Dict[str, int] = None,
                       coords: np.ndarray = None,
                       nums: Iterable[int] = None,
                       datapoint: BACDatapoint = None,
                       spc: ReferenceSpecies = None,
                       multiplicity: int = None) -> float:
        """
        Returns the bond additivity correction in J/mol.

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
            The bond correction to the electronic energy in J/mol.
        """
        if self.bacs is None:
            bac_type_str = 'Melius' if self.bac_type == 'm' else 'Petersson'
            raise BondAdditivityCorrectionError(
                f'Missing {bac_type_str}-type BAC parameters for model chemistry {self.model_chemistry}'
            )

        if datapoint is None and spc is not None:
            datapoint = BACDatapoint(spc, model_chemistry=self.model_chemistry)

        if self.bac_type == 'm':
            return self._get_melius_correction(coords=coords, nums=nums, datapoint=datapoint, multiplicity=multiplicity)
        elif self.bac_type == 'p':
            return self._get_petersson_correction(bonds=bonds, datapoint=datapoint)

    def _get_petersson_correction(self, bonds: Dict[str, int] = None, datapoint: BACDatapoint = None) -> float:
        """
        Given the model_chemistry and a dictionary of bonds, return the
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
            Petersson-type bond additivity correction in J/mol.
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
                symbol_flipped = ''.join(re.findall('[a-zA-Z]+|[^a-zA-Z]+', symbol)[::-1])  # Check reversed symbol
                if symbol_flipped in self.bacs:
                    bac += count * self.bacs[symbol_flipped]
                else:
                    logging.warning(f'Bond correction not applied for unknown bond type {symbol}.')

        return bac * 4184.0  # Convert kcal/mol to J/mol

    def _get_melius_correction(self,
                               coords: np.ndarray = None,
                               nums: Iterable[int] = None,
                               datapoint: BACDatapoint = None,
                               multiplicity: int = None,
                               params: Dict[str, Union[float, Dict[str, float]]] = None) -> float:
        """
        Given the model chemistry, molecular coordinates, atomic numbers,
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
            Melius-type bond additivity correction in J/mol.
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
        spin = 0.5 * (multiplicity - 1)
        bac_mol = mol_corr * (spin - sum(self.atom_spins[atom.element.symbol] for atom in mol.atoms))

        # Atomic correction
        bac_atom = sum(atom_corr[atom.element.symbol] for atom in mol.atoms)

        # Bond correction
        bac_bond = 0.0
        for bond in mol.get_all_edges():
            atom1 = bond.atom1
            atom2 = bond.atom2
            symbol1 = atom1.element.symbol
            symbol2 = atom2.element.symbol

            # Bond length correction
            length_corr = (bond_corr_length[symbol1] * bond_corr_length[symbol2]) ** 0.5
            length = np.linalg.norm(atom1.coords - atom2.coords)
            bac_bond += length_corr * np.exp(-self.exp_coeff * length)

            # Neighbor correction
            for other_atom, other_bond in mol.get_bonds(atom1).items():  # Atoms adjacent to atom1
                if other_bond is not bond:
                    other_symbol = other_atom.element.symbol
                    bac_bond += bond_corr_neighbor[symbol1] + bond_corr_neighbor[other_symbol]
            for other_atom, other_bond in mol.get_bonds(atom2).items():  # Atoms adjacent to atom2
                if other_bond is not bond:
                    other_symbol = other_atom.element.symbol
                    bac_bond += bond_corr_neighbor[symbol2] + bond_corr_neighbor[other_symbol]

        # Note the minus sign
        return -(bac_mol + bac_atom + bac_bond) * 4184.0  # Convert kcal/mol to J/mol

    def fit(self, db_names: Union[str, List[str]] = 'main', **kwargs):
        """
        Fits bond additivity corrections using calculated and reference
        data available in the RMG database. The resulting BACs stored
        in self.bacs will be based on kcal/mol.

        Args:
            db_names: Optionally specify database names to train on (defaults to main).
            kwargs: Keyword arguments for fitting Melius-type BACs (see self._fit_melius).
        """
        self.database_key = self.load_database(names=db_names)

        self.dataset = extract_dataset(self.ref_databases[self.database_key], self.model_chemistry)
        if len(self.dataset) == 0:
            raise BondAdditivityCorrectionError(f'No species available for {self.model_chemistry} model chemistry')

        if self.bac_type == 'm':
            logging.info(f'Fitting Melius-type BACs for {self.model_chemistry}...')
            self._fit_melius(**kwargs)
        elif self.bac_type == 'p':
            logging.info(f'Fitting Petersson-type BACs for {self.model_chemistry}...')
            self._fit_petersson()

        stats_before = self.dataset.calculate_stats()
        stats_after = self.dataset.calculate_stats(for_bac_data=True)
        logging.info(f'RMSE/MAE before fitting: {stats_before.rmse:.2f}/{stats_before.mae:.2f} kcal/mol')
        logging.info(f'RMSE/MAE after fitting: {stats_after.rmse:.2f}/{stats_after.mae:.2f} kcal/mol')

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

        x = make_feature_mat(features)
        y = self.dataset.ref_data - self.dataset.calc_data
        w = np.linalg.solve(x.T @ x, x.T @ y)
        ypred = x @ w

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
                [self._get_melius_correction(datapoint=d, params=get_params(_w)) / 4184 for d in self.dataset]
            )
            return self.dataset.calc_data + corr

        def residuals(_w: np.ndarray) -> Union[float, np.ndarray]:
            """Calculate residuals"""
            bac_data = get_bac_data(_w)
            return self.dataset.ref_data - bac_data

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
        for i, line in enumerate(lines):
            if keyword in line:
                if has_entries:
                    if self.model_chemistry in bac_dict:
                        if overwrite:
                            # Does not overwrite comments
                            del_idx_start = del_idx_end = None
                            for j, line2 in enumerate(lines[i:]):
                                if self.model_chemistry in line2:
                                    del_idx_start = i + j
                                    del_idx_end = None
                                elif line2.rstrip() == '    },':  # Can't have comment after final brace
                                    del_idx_end = i + j + 1
                                if del_idx_start is not None and del_idx_end is not None:
                                    if (lines[del_idx_start - 1].lstrip().startswith('#')
                                            or lines[del_idx_end + 1].lstrip().startswith('#')):
                                        logging.warning('There may be left over comments from previous BACs')
                                    lines[del_idx_start:del_idx_end] = bacs_formatted
                                    break
                        else:
                            raise IOError(
                                f'{self.model_chemistry} model chemistry already exists. Set `overwrite` to True.'
                            )
                    else:
                        lines[(i+1):(i+1)] = ['\n'] + bacs_formatted
                else:
                    lines[i] = f'{keyword} = {{\n'
                    lines[(i+1):(i+1)] = ['\n'] + bacs_formatted + ['\n}\n']
                break

        with open(data_path if alternate_path is None else alternate_path, 'w') as f:
            f.writelines(lines)

        # Reload data to update BAC dictionaries
        if alternate_path is None:
            importlib.reload(data)

    def format_bacs(self, indent=False):
        """
        Obtain a list of nicely formatted BACs suitable for writelines.

        Args:
            indent: Indent each line for printing in database.

        Returns:
            Formatted list of BACs.
        """
        bacs_formatted = json.dumps(self.bacs, indent=4).replace('"', "'").split('\n')
        bacs_formatted[0] = f"'{self.model_chemistry}': " + bacs_formatted[0]
        bacs_formatted[-1] += ','
        bacs_formatted = [e + '\n' for e in bacs_formatted]
        if indent:
            bacs_formatted = ['    ' + e for e in bacs_formatted]
        return bacs_formatted
