"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 867 $"

import openbabel as ob

def makeopenbabel(atomcoords, atomnos, charge=0, mult=1):
    """Create an Open Babel molecule.

    >>> import numpy, openbabel
    >>> atomnos = numpy.array([1,8,1],"i")
    >>> coords = numpy.array([[-1.,1.,0.],[0.,0.,0.],[1.,1.,0.]])
    >>> obmol = makeopenbabel(coords, atomnos)
    >>> obconversion = openbabel.OBConversion()
    >>> formatok = obconversion.SetOutFormat("inchi")
    >>> print obconversion.WriteString(obmol).strip()
    InChI=1/H2O/h1H2
    """
    obmol = ob.OBMol()
    for i in range(len(atomnos)):
        # Note that list(atomcoords[i]) is not equivalent!!!
        coords = atomcoords[i].tolist()
        atomno = int(atomnos[i])
        obatom = ob.OBAtom()
        obatom.SetAtomicNum(atomno)
        obatom.SetVector(*coords)
        obmol.AddAtom(obatom)
    obmol.ConnectTheDots()
    obmol.PerceiveBondOrders()
    obmol.SetTotalSpinMultiplicity(mult)
    obmol.SetTotalCharge(charge)
    return obmol

if __name__ == "__main__":
    import doctest
    doctest.testmod()
