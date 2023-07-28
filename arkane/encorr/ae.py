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
This module provides classes for fitting atom energies based on a very
small, predetermined set of molecules.
"""

import importlib
import json
import logging
from collections import Counter
from typing import Dict, Hashable, List, Union

import numpy as np
from scipy.stats import distributions

from rmgpy import constants
from rmgpy.molecule import get_element, Molecule

import arkane.encorr.data as data
from arkane.encorr.reference import ReferenceDatabase
from arkane.modelchem import LevelOfTheory, CompositeLevelOfTheory

# List of species labels that will be used for fitting (labels should match reference database)
SPECIES_LABELS = [
    'Dihydrogen',
    'Dinitrogen',
    'Dioxygen',
    'Disulfur',
    'Difluorine',
    'Dichlorine',
    'Dibromine',
    'Hydrogen fluoride',
    'Hydrogen chloride',
    'Hydrogen bromide',
    'Hydrogen sulfide',
    'Water',
    'Methane',
    'Methyl',
    'Ammonia',
    'Chloromethane',
    'Lithium Hydride',
    'Lithium Fluoride'
]


class AEJob:
    """
    A job for fitting atom energies.
    """

    def __init__(self,
                 species_energies: Dict[str, float],
                 level_of_theory: Union[LevelOfTheory, CompositeLevelOfTheory] = None,
                 write_to_database: bool = False,
                 overwrite: bool = False):
        """
        Initialize an AEJob instance.

        Notes:
            The species energies should be provided as a dictionary
            containing the species labels as keys and their single-
            point electronic energies in Hartree as values. The
            energies should be calculated using the experimental
            geometry provided for the species in the reference
            database, and the zero-point energy should not be included
            in the electronic energy.

        Args:
            species_energies: Dictionary of species labels with single-point electronic energies (Hartree).
            level_of_theory: Dictionary key for saving atom energies to the database.
            write_to_database: Save the fitted atom energies directly to the RMG database.
            overwrite: Overwrite atom energies in the RMG database if they already exist.
        """
        self.spcs_energies = species_energies
        self.level_of_theory = level_of_theory
        self.write_to_database = write_to_database
        self.overwrite = overwrite
        self.ae = AE(species_energies)

    def execute(self, output_file: str = None):
        """
        Execute the atom energy job.

        Args:
            output_file: Write the fitted energies to this file.
        """
        if self.level_of_theory is None:
            logging.info('Fitting atom energies')
        else:
            logging.info(f'Fitting atom energies for {self.level_of_theory}')
        self.ae.fit()

        if output_file is not None:
            with open(output_file, 'a') as f:
                if self.level_of_theory is not None:
                    f.write(f'#  {self.level_of_theory}\n')
                for element, energy in self.ae.atom_energies.items():
                    f.write(f'# {element:2}: {energy:15.8f} +/- {self.ae.confidence_intervals[element]:.8f} Hartree\n')
                f.writelines(self.ae.format_atom_energies(
                    'atom_energies' if self.level_of_theory is None else self.level_of_theory))

        if self.write_to_database:
            if self.level_of_theory is None:
                raise Exception('Level of theory is required for writing to database')
            try:
                self.ae.write_to_database(self.level_of_theory, overwrite=self.overwrite)
            except ValueError as e:
                logging.warning('Could not write atom energies to database. Captured error:')
                logging.warning(str(e))


class AE:
    """
    A class for fitting atom energies.
    """

    ref_data_src = 'CCCBDB'  # Use CCCBDB data
    ref_data = None  # Dictionary of reference data entries

    def __init__(self, species_energies: Dict[str, float]):
        self.species_energies = species_energies  # Hartree
        self.atom_energies = None
        self.confidence_intervals = None

        for lbl in SPECIES_LABELS:
            if lbl not in self.species_energies:
                logging.warning(f'{lbl} missing from provided species energies!')

    @classmethod
    def _load_refdata(cls):
        if cls.ref_data is None:
            logging.info('Loading reference database')
            db = ReferenceDatabase()
            db.load()
            cls.ref_data = {lbl: spc for lbl, spc in zip(SPECIES_LABELS, db.get_species_from_label(SPECIES_LABELS))}

    def fit(self):
        """
        Fit atom energies using the provided species energies and
        corresponding atomization energies from the reference data.
        """
        self._load_refdata()

        mols = [
            Molecule().from_adjacency_list(
                self.ref_data[lbl].adjacency_list,
                raise_atomtype_exception=False,
                raise_charge_exception=False
            ) for lbl in self.species_energies
        ]
        atom_counts = [Counter(atom.element.symbol for atom in mol.atoms) for mol in mols]
        elements = sorted({element for ac in atom_counts for element in ac}, key=lambda s: get_element(s).number)
        x = np.array([[ac[element] for element in elements] for ac in atom_counts])  # Nmols x Nelements

        atomization_energies = np.array([
            self.ref_data[lbl].reference_data[self.ref_data_src].atomization_energy.value_si
            / constants.E_h / constants.Na for lbl in self.species_energies
        ])
        zpes = np.array([
            self.ref_data[lbl].reference_data[self.ref_data_src].zpe.value_si
            / constants.E_h / constants.Na for lbl in self.species_energies
        ])
        elec_energies = np.array(list(self.species_energies.values()))  # Should already be in Hartree
        y = atomization_energies + elec_energies + zpes

        w = np.linalg.solve(x.T @ x, x.T @ y)
        self.atom_energies = dict(zip(elements, w))

        # Get confidence intervals
        n = len(y)  # Ndata
        k = len(w)  # Nparam
        ypred = x @ w
        sigma2 = np.sum((y - ypred)**2) / (n - k)  # MSE
        cov = sigma2 * np.linalg.inv(x.T @ x)  # covariance matrix
        se = np.sqrt(np.diag(cov))  # standard error
        alpha = 0.05  # 95% confidence level
        tdist = distributions.t.ppf(1 - alpha/2, n - k)  # student-t
        ci = tdist * se  # confidence interval half-width
        self.confidence_intervals = dict(zip(elements, ci))  # Parameter estimates are w +/- ci

    def write_to_database(self, key: Hashable, overwrite: bool = False, alternate_path: str = None):
        """
        Write atom energies to database.

        Args:
            key: Dictionary key to use for atom energies in database.
            overwrite: Overwrite existing atom energies.
            alternate_path: Write atom energies and existing database to this path instead.
        """
        if self.atom_energies is None:
            raise ValueError('No atom energies available for writing')

        data_path = data.quantum_corrections_path
        with open(data_path) as f:
            lines = f.readlines()

        ae_formatted = self.format_atom_energies(key, indent=True)

        # Add new atom energies to file without changing existing formatting
        for i, line in enumerate(lines):
            if 'atom_energies' in line:
                if key in data.atom_energies:
                    if overwrite:
                        # Does not overwrite comments
                        del_idx_start = del_idx_end = None
                        for j, line2 in enumerate(lines[i:]):
                            if repr(key) in line2:
                                del_idx_start = i + j
                                del_idx_end = None
                            elif line2.rstrip() == '    },':  # Can't have a comment after final brace
                                del_idx_end = i + j + 1
                            if del_idx_start is not None and del_idx_end is not None:
                                if (lines[del_idx_start - 1].lstrip().startswith('#')
                                        or lines[del_idx_end + 1].lstrip().startswith('#')):
                                    logging.warning('There may be left over comments from previous atom energies')
                                lines[del_idx_start:del_idx_end] = ae_formatted
                                break
                    else:
                        raise ValueError(f'{key} already exists. Set `overwrite` to True.')
                else:
                    lines[(i+1):(i+1)] = ['\n'] + ae_formatted
                break

        with open(data_path if alternate_path is None else alternate_path, 'w') as f:
            f.writelines(lines)

        # Reload data to update atom energy dictionary
        if alternate_path is None:
            importlib.reload(data)

    def format_atom_energies(self, key: Hashable, indent: bool = False) -> List[str]:
        """
        Obtain a list of nicely formatted atom energies suitable for
        writelines.

        Args:
            key: Dictionary key to use for formatting dictionary.
            indent: Indent each line.

        Returns:
            Formatted list of atom energies.
        """
        ae_formatted = json.dumps(self.atom_energies, indent=4).replace('"', "'").split('\n')
        ae_formatted[0] = f'"{key}": ' + ae_formatted[0]
        ae_formatted[-1] += ','
        ae_formatted = [e + '\n' for e in ae_formatted]
        if indent:
            ae_formatted = ['    ' + e for e in ae_formatted]
        return ae_formatted
