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
A module for generating ESS adapters
"""

import os
from typing import Type

from arkane.ess.adapter import ESSAdapter

from rmgpy.exceptions import InputError

_registered_ess_adapters = {}


def register_ess_adapter(ess: str,
                         ess_class: Type[ESSAdapter],
                         ) -> None:
    """
    A register for the ESS adapters.

    Args:
        ess: A string representation for an ESS adapter.
        ess_class: The ESS adapter class.

    Raises:
        TypeError: If ``ess_class`` is not an ``ESSAdapter`` instance.
    """
    if not issubclass(ess_class, ESSAdapter):
        raise TypeError(f'{ess_class} is not an ESSAdapter')
    _registered_ess_adapters[ess] = ess_class


def ess_factory(fullpath: str,
                check_for_errors: bool = True,
                scratch_directory: str = None,
                ) -> Type[ESSAdapter]:
    """
    A factory generating the ESS adapter corresponding to ``ess_adapter``.
    Given a path to the log file of a QM software, determine whether it is
    Gaussian, Molpro, Orca, Psi4, QChem, or TeraChem.

    Args:
        fullpath (str): The disk location of the output file of interest.
        check_for_errors (bool): Boolean indicating whether to check the QM log for common errors
                                 before parsing relevant information.

    Returns:
        Type[ESSAdapter]: The requested ESSAdapter child, initialized with the respective arguments.
    """

    ess_name = None
    if os.path.splitext(fullpath)[-1] in ['.xyz', '.dat', '.geometry']:
        ess_name = 'TeraChemLog'
    else:
        with open(fullpath, 'r') as f:
            line = f.readline().lower()
            while ess_name is None and line != '':
                if 'gaussian' in line:
                    ess_name = 'GaussianLog'
                    break
                elif 'molpro' in line:
                    ess_name = 'MolproLog'
                    break
                elif 'o   r   c   a' in line or 'orca' in line:
                    ess_name = 'OrcaLog'
                    break
                elif 'psi4' in line or 'rob parrish' in line:
                    ess_name = 'Psi4Log'
                    break
                elif 'qchem' in line:
                    ess_name = 'QChemLog'
                    break
                elif 'terachem' in line:
                    ess_name = 'TeraChemLog'
                    break
                line = f.readline().lower()
    if ess_name is None:
        raise InputError(f'The file at {fullpath} could not be identified as a '
                         f'Gaussian, Molpro, Orca, Psi4, QChem, or TeraChem log file.')

    return _registered_ess_adapters[ess_name](path=fullpath, check_for_errors=check_for_errors, scratch_directory=scratch_directory)
