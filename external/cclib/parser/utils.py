"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 865 $"


def convertor(value, fromunits, tounits):
    """Convert from one set of units to another.

    >>> print "%.1f" % convertor(8, "eV", "cm-1")
    64524.8
    """
    
    _convertor = {"eV_to_cm-1": lambda x: x*8065.6,
                  "hartree_to_eV": lambda x: x*27.2113845,
                  "bohr_to_Angstrom": lambda x: x*0.529177,
                  "Angstrom_to_bohr": lambda x: x*1.889716,
                  "nm_to_cm-1": lambda x: 1e7/x,
                  "cm-1_to_nm": lambda x: 1e7/x,
                  "hartree_to_cm-1": lambda x: x*219474.6,
                  # Taken from GAMESS docs, "Further information",
                  # "Molecular Properties and Conversion Factors"
                  "Debye^2/amu-Angstrom^2_to_km/mol": lambda x: x*42.255}

    return _convertor["%s_to_%s" % (fromunits, tounits)] (value)


class PeriodicTable(object):
    """Allows conversion between element name and atomic no.

    >>> t = PeriodicTable()
    >>> t.element[6]
    'C'
    >>> t.number['C']
    6
    >>> t.element[44]
    'Ru'
    >>> t.number['Au']
    79
    """
    
    def __init__(self):
        self.element = [None,
            'H', 'He',
            'Li', 'Be',
            'B', 'C', 'N', 'O', 'F', 'Ne',
            'Na', 'Mg',
            'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
            'K', 'Ca',
            'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
            'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
            'Rb', 'Sr',
            'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
            'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
            'Cs', 'Ba',
            'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
            'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
            'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
            'Fr', 'Ra',
            'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
            'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Uub']
        self.number = {}
        for i in range(1, len(self.element)):
            self.number[self.element[i]] = i


if __name__ == "__main__":
    import doctest, utils
    doctest.testmod(utils, verbose=False)
