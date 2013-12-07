"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 709 $"

from Bio.PDB.Atom import Atom
from cclib.parser.utils import PeriodicTable

def makebiopython(atomcoords, atomnos):
    """Create a list of BioPython Atoms.

    This creates a list of BioPython Atoms suitable for use
    by Bio.PDB.Superimposer, for example.

    >>> import numpy
    >>> from Bio.PDB.Superimposer import Superimposer
    >>> atomnos = numpy.array([1,8,1],"i")
    >>> a = numpy.array([[-1,1,0],[0,0,0],[1,1,0]],"f")
    >>> b = numpy.array([[1.1,2,0],[1,1,0],[2,1,0]],"f")
    >>> si = Superimposer()
    >>> si.set_atoms(makebiopython(a,atomnos),makebiopython(b,atomnos))
    >>> print si.rms
    0.29337859596
    """
    pt = PeriodicTable()
    bioatoms = []
    for coords, atomno in zip(atomcoords, atomnos):
        bioatoms.append(Atom(pt.element[atomno], coords, 0, 0, 0, 0, 0))
    return bioatoms

if __name__ == "__main__":
    import doctest
    doctest.testmod()
