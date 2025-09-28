"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 737 $"

from PyQuante.Molecule import Molecule

def makepyquante(atomcoords, atomnos, charge=0, mult=1):
    """Create a PyQuante Molecule.

    >>> import numpy
    >>> from PyQuante.hartree_fock import hf
    >>> atomnos = numpy.array([1,8,1],"i")
    >>> a = numpy.array([[-1,1,0],[0,0,0],[1,1,0]],"f")
    >>> pyqmol = makepyquante(a,atomnos)
    >>> en,orbe,orbs = hf(pyqmol)
    >>> print int(en * 10) / 10. # Should be around -73.8
    -73.8
    """
    return Molecule("notitle", zip(atomnos, atomcoords), units="Angstrom",
                    charge=charge, multiplicity=mult)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
