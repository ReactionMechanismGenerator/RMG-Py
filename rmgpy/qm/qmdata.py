#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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

import logging
import re

from rmgpy.quantity import Energy, Mass, Length, Frequency


class QMData(object):
    """
    General class for data extracted from a QM calculation
    """

    def __init__(self,
                 groundStateDegeneracy=-1,
                 numberOfAtoms=None,
                 stericEnergy=None,
                 molecularMass=None,
                 energy=0,
                 atomicNumbers=None,
                 rotationalConstants=None,
                 atomCoords=None,
                 frequencies=None,
                 source=None,
                 ):
        #: Electronic ground state degeneracy in RMG taken as number of radicals +1
        self.groundStateDegeneracy = groundStateDegeneracy
        self.numberOfAtoms = numberOfAtoms  #: Number of atoms.
        self.stericEnergy = Energy(stericEnergy)
        self.molecularMass = Mass(molecularMass)
        self.energy = Energy(energy)
        self.atomicNumbers = atomicNumbers
        self.rotationalConstants = Frequency(rotationalConstants)
        self.atomCoords = Length(atomCoords)
        self.frequencies = Frequency(frequencies)
        self.source = source

        self.test_valid()

    def test_valid(self):
        assert self.groundStateDegeneracy > 0

    def __repr__(self):
        things = []
        for attribute in ['groundStateDegeneracy',
                          'numberOfAtoms',
                          'stericEnergy',
                          'molecularMass',
                          'energy',
                          'atomicNumbers',
                          'rotationalConstants',
                          'atomCoords',
                          'frequencies',
                          'source',
                          ]:
            things.append("{0!s}={1!r}".format(attribute, getattr(self, attribute)))
        string = ', '.join(things)
        string = re.sub('\s+', ' ', string)
        return 'QMData({0!s})'.format(string)


def parse_cclib_data(cclib_data, ground_state_degeneracy):
    """
    Parses a CCLib data object and returns QMData object.
    """
    try:
        number_of_atoms = cclib_data.natom
        if hasattr(cclib_data, 'molmass'):
            molecular_mass = (cclib_data.molmass, 'amu')
        else:
            molecular_mass = None
        energy = (cclib_data.scfenergies[-1], 'eV/molecule')
        atomic_numbers = cclib_data.atomnos
        rotational_constants = ([i * 1e9 for i in cclib_data.rotcons[-1]], 'hertz')
        atom_coords = (cclib_data.atomcoords[-1], 'angstrom')
        frequencies = (cclib_data.vibfreqs, 'cm^-1')

    except AttributeError:
        logging.error("The passed in cclibData has these attributes: {0!r}".format(cclib_data._attrlist))
        raise

    if hasattr(cclib_data, 'stericenergy'):
        stericEnergy = (cclib_data.stericenergy, 'eV/molecule')
    else:
        stericEnergy = None

    return QMData(
        groundStateDegeneracy=ground_state_degeneracy,
        numberOfAtoms=number_of_atoms,
        stericEnergy=stericEnergy,
        molecularMass=molecular_mass,
        energy=energy,
        atomicNumbers=atomic_numbers,
        rotationalConstants=rotational_constants,
        atomCoords=atom_coords,
        frequencies=frequencies
    )
