# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Calculate properties of nuclei based on data parsed by cclib."""

import logging

import numpy

from .calculationmethod import Method


class Nuclear(Method):
    """A container for methods pertaining to atomic nuclei."""
    
    def __init__(self, data, progress=None, loglevel=logging.INFO, logname="Log"):

        super(Nuclear, self).__init__(data, progress, loglevel, logname)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Nuclear"

    def __repr__(self):
        """Return a representation of the object."""
        return "Nuclear"
    
    def repulsion_energy(self):
        """Return the nuclear repulsion energy."""

        nre = 0.0
        for i in range(self.data.natom):
            ri = self.data.atomcoords[0][i]
            zi = self.data.atomnos[i]
            for j in range(i+1, self.data.natom):
                rj = self.data.atomcoords[0][j]
                zj = self.data.atomnos[j]
                d = numpy.linalg.norm(ri-rj)
                nre += zi*zj/d
        return nre


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
