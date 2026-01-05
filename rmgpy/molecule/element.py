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
This module defined the chemical elements available in RMG. Information for
each element is stored as attributes of an object of the :class:`Element`
class.

Element objects for each chemical element (1-112) have also been declared as
module-level variables, using each element's symbol as its variable name. The
:meth:`get_element` method can also be used to retrieve the :class:`Element`
object associated with an atomic number or symbol. Generally applications will
want to use these objects, both to conserve memory and to make for easy
comparisons.
"""

import cython
from rdkit.Chem import GetPeriodicTable

from rmgpy.exceptions import ElementError
from rmgpy.quantity import Quantity

################################################################################

_rdkit_periodic_table = GetPeriodicTable()


class Element(object):
    """
    A chemical element. The attributes are:

    ============== =============== ================================================
    Attribute      Type            Description
    ============== =============== ================================================
    `number`       ``int``         The atomic number of the element
    `symbol`       ``str``         The symbol used for the element
    `name`         ``str``         The IUPAC name of the element
    `mass`         ``float``       The mass of the element in kg/mol
    `cov_radius`   ``float``       Covalent bond radius in Angstrom
    `isotope`      ``int``         The isotope integer of the element
    `chemkin_name` ``str``         The chemkin compatible representation of the element
    ============== =============== ================================================

    This class is specifically for properties that all atoms of the same element
    share. Ideally there is only one instance of this class for each element.
    """

    def __init__(self, number, symbol, name, mass, isotope=-1, chemkin_name=None):
        self.number = number
        self.symbol = symbol
        self.name = name
        self.mass = mass
        self.isotope = isotope
        self.chemkin_name = chemkin_name or self.name
        if symbol in {'X','L','R','e'}:
            self.cov_radius = 0
        else:
            try:
                self.cov_radius = _rdkit_periodic_table.GetRcovalent(symbol)
            except RuntimeError:
                import logging
                logging.error("RDkit doesn't know element {0} so covalent radius unknown".format(symbol))
                self.cov_radius = 0

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return self.symbol

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "Element(%s, '%s', '%s', %s, %s)" % (self.number, self.symbol, self.name, self.mass, self.isotope)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Element, (self.number, self.symbol, self.name, self.mass, self.isotope))


class PeriodicSystem(object):
    """
    Collects hard-coded information of elements in periodic table.

    Currently it has static attributes:
    `valences`: the number of bonds an element is able to form around it.
    `valence_electrons`: the number of electrons in  the outermost shell of the element that can
    participate in bond formation. for instance, `S` has 6 outermost electrons (valence_electrons)
    but 4 of them form lone pairs, the remaining 2 electrons can form bonds so the normal valence
    for `S` is 2.
    `lone_pairs`: the number of lone pairs an element has
    'electronegativity': The atom's electronegativity (how well can an atom attract electrons), taken from
    https://sciencenotes.org/list-of-electronegativity-values-of-the-elements/
    isotopes of the same element may have slight different electronegativities, which is not reflected below
    """
    valences = {'H+':0, 'e': 0, 'H': 1, 'He': 0, 'C': 4, 'N': 3, 'O': 2, 'F': 1, 'Ne': 0,
                'Si': 4, 'P': 3, 'S': 2, 'Cl': 1, 'Br': 1, 'Ar': 0, 'I': 1, 'X': 4, 'Li': 1}
    valence_electrons = {'H+':0, 'e': 1, 'H': 1, 'He': 2, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'Ne': 8,
                         'Si': 4, 'P': 5, 'S': 6, 'Cl': 7, 'Br': 7, 'Ar': 8, 'I': 7, 'X': 4, 'Li': 1}
    lone_pairs = {'H+':0, 'e': 0, 'H': 0, 'He': 1, 'C': 0, 'N': 1, 'O': 2, 'F': 3, 'Ne': 4,
                  'Si': 0, 'P': 1, 'S': 2, 'Cl': 3, 'Br': 3, 'Ar': 4, 'I': 3, 'X': 0, 'Li': 0}
    electronegativity = {'H': 2.20, 'D': 2.20, 'T': 2.20, 'C': 2.55, 'C13': 2.55, 'N': 3.04, 'O': 3.44, 'O18': 3.44,
                         'F': 3.98, 'Si': 1.90, 'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Br': 2.96, 'I': 2.66, 'X': 0.0,
                         'Li' : 0.98}


################################################################################


def get_element(value, isotope=-1):
    """
    Return the :class:`Element` object corresponding to the given parameter
    `value`. If an integer is provided, the value is treated as the atomic
    number. If a string is provided, the value is treated as the symbol. An
    :class:`ElementError` is raised if no matching element is found.
    """
    cython.declare(element=Element, number=cython.int, symbol=str)
    if isinstance(value, int):
        # The parameter is an integer; assume this is the atomic number
        number = value
        for element in element_list:
            if element.number == number and element.isotope == isotope:
                return element
        # If we reach this point that means we did not find an appropriate element,
        # so we raise an exception
        raise ElementError("No element found with atomic number %i, and isotope %i" % (number, isotope))
    elif isinstance(value, str):
        # The parameter is a string; assume this is the element symbol
        symbol = value
        for element in element_list:
            if element.symbol == symbol and element.isotope == isotope:
                return element
        # If we reach this point that means we did not find an appropriate element,
        # so we raise an exception
        raise ElementError("No element found with symbol %s, and isotope %i." % (symbol, isotope))
    else:
        raise ElementError('No element found based on parameter %s "%r", isotope: %i.' % (type(value), value, isotope))


################################################################################

# The variable names correspond to each element's symbol
# The elements are sorted by increasing atomic number and grouped by period
#Taken from https://iupac.qmul.ac.uk/AtWt/ (IUPAC 2023, accessed November 2025)
#The molar mass constant Mu = 1.0000000010531e-3 kg/mol was used to convert from dimensionaless atomic weights to kg/mol. 
# Recommended IUPAC nomenclature is used throughout (including 'aluminium' and
# 'caesium')

# electron
e = Element(-1,   'e', 'electron'      , 5.485799090441e-7)

# Surface site
X = Element(0,    'X', 'surface_site'   , 0.0)

# Period 1
H  = Element(1,   'H' , 'hydrogen'	, 0.00100802000106)
D  = Element(1,   'H' , 'deuterium'     , 0.00201410177818, 2, 'D')  
T  = Element(1,   'H' , 'tritium'       , 0.00301604928, 3, 'T')
He = Element(2,   'He', 'helium'	, 0.00400260220422)

# Period 2
Li = Element(3,   'Li', 'lithium'	, 0.00694600000731)
Be = Element(4,   'Be', 'beryllium'	, 0.00901218315949)
B  = Element(5,   'B' , 'boron'		, 0.01081200001139)
C  = Element(6,   'C' , 'carbon'	, 0.01201120001265)
C13= Element(6,   'C' , 'carbon-13'     , 0.0130033548352, 13, 'C13') 
N  = Element(7,   'N' , 'nitrogen'	, 0.01400710001475)
N15= Element(7,   'N' , 'nitrogen-15'	, 0.0150001088994, 15, 'N15')
O  = Element(8,   'O' , 'oxygen'	, 0.01599910001685)
O17= Element(8,   'O' , 'oxygen-17'     , 0.0169991317575, 17, 'O17')
O18= Element(8,   'O' , 'oxygen-18'     , 0.0179991596135, 18, 'O18')
F  = Element(9,   'F' , 'fluorine'	, 0.01899840318251)
Ne = Element(10,  'Ne', 'neon'		, 0.02017976002125)

# Period 3
Na = Element(11,  'Na', 'sodium'	, 0.02298976930621)
Mg = Element(12,  'Mg', 'magnesium'	, 0.0243052000256)
Al = Element(13,  'Al', 'aluminium'	, 0.02698153845841)
Si = Element(14,  'Si', 'silicon'	, 0.02808510002958)
P  = Element(15,  'P' , 'phosphorus'	, 0.03097376203112)
S  = Element(16,  'S' , 'sulfur'	, 0.03206200003376)
Cl = Element(17,  'Cl', 'chlorine'	, 0.03545100003733)
Ar = Element(18,  'Ar', 'argon'		, 0.03995160004207)

# Period 4
K  = Element(19,  'K' , 'potassium'	, 0.03909831004117)
Ca = Element(20,  'Ca', 'calcium'	, 0.04007840004221)
Sc = Element(21,  'Sc', 'scandium'	, 0.04495590744734)
Ti = Element(22,  'Ti', 'titanium'	, 0.04786710005041)
V  = Element(23,  'V' , 'vanadium'	, 0.05094151005365)
Cr = Element(24,  'Cr', 'chromium'	, 0.05199616005476)
Mn = Element(25,  'Mn', 'manganese'	, 0.05493804325786)
Fe = Element(26,  'Fe', 'iron'		, 0.05584520005881)
Co = Element(27,  'Co', 'cobalt'	, 0.05893319436206)
Ni = Element(28,  'Ni', 'nickel'	, 0.05869344006181)
Cu = Element(29,  'Cu', 'copper'	, 0.06354630006692)
Zn = Element(30,  'Zn', 'zinc'		, 0.06538200006885)
Ga = Element(31,  'Ga', 'gallium'	, 0.06972310007343)
Ge = Element(32,  'Ge', 'germanium'	, 0.07263080007649)
As = Element(33,  'As', 'arsenic'	, 0.0749215956789)
Se = Element(34,  'Se', 'selenium'	, 0.07897180008317)
Br = Element(35,  'Br', 'bromine'	, 0.07990430008415)
Kr = Element(36,  'Kr', 'krypton'	, 0.08379820008825)

# Period 5
Rb = Element(37,  'Rb', 'rubidium'	, 0.08546783009001)
Sr = Element(38,  'Sr', 'strontium'	, 0.08762100009227)
Y  = Element(39,  'Y' , 'yttrium'	, 0.08890583829363)
Zr = Element(40,  'Zr', 'zirconium'	, 0.09122230009607)
Nb = Element(41,  'Nb', 'niobium'	, 0.09290637109784)
Mo = Element(42,  'Mo', 'molybdenum'	, 0.09595100010105)
Ru = Element(44,  'Ru', 'ruthenium'	, 0.10107200010644)
Rh = Element(45,  'Rh', 'rhodium'	, 0.10290549210837)
Pd = Element(46,  'Pd', 'palladium'	, 0.10642100011207)
Ag = Element(47,  'Ag', 'silver'	, 0.1078682201136)
Cd = Element(48,  'Cd', 'cadmium'	, 0.11241440011838)
In = Element(49,  'In', 'indium'	, 0.11481810012091)
Sn = Element(50,  'Sn', 'tin'		, 0.11871070012501)
Sb = Element(51,  'Sb', 'antimony'	, 0.12176010012823)
Te = Element(52,  'Te', 'tellurium'	, 0.12760300013438)
I  = Element(53,  'I' , 'iodine'	, 0.12690447313364)
Xe = Element(54,  'Xe', 'xenon'		, 0.13129360013827)

# Period 6
Cs = Element(55,  'Cs', 'caesium'	, 0.13290545210596)
Ba = Element(56,  'Ba', 'barium'	, 0.13732770014462)
La = Element(57,  'La', 'lanthanum'	, 0.13890547714628)
Ce = Element(58,  'Ce', 'cerium'	, 0.14011610014756)
Pr = Element(59,  'Pr', 'praseodymium'	, 0.14090766114839)
Nd = Element(60,  'Nd', 'neodymium'	, 0.1442423001519)
Sm = Element(62,  'Sm', 'samarium'	, 0.15036200015835)
Eu = Element(63,  'Eu', 'europium'	, 0.15196410016003)
Gd = Element(64,  'Gd', 'gadolinium'	, 0.1572492001656)
Tb = Element(65,  'Tb', 'terbium'	, 0.15892535486736)
Dy = Element(66,  'Dy', 'dysprosium'	, 0.16250010017113)
Ho = Element(67,  'Ho', 'holmium'	, 0.16493032967369)
Er = Element(68,  'Er', 'erbium'	, 0.16725930017614)
Tm = Element(69,  'Tm', 'thulium'	, 0.1689342196779)
Yb = Element(70,  'Yb', 'ytterbium'	, 0.17304510018223)
Lu = Element(71,  'Lu', 'lutetium'	, 0.17496669518426)
Hf = Element(72,  'Hf', 'hafnium'	, 0.17848660018796)
Ta = Element(73,  'Ta', 'tantalum'	, 0.18094788219056)
W  = Element(74,  'W' , 'tungsten'	, 0.1838410001936)
Re = Element(75,  'Re', 'rhenium'	, 0.18620710019609)
Os = Element(76,  'Os', 'osmium'	, 0.19023300020033)
Ir = Element(77,  'Ir', 'iridium'	, 0.19221720020242)
Pt = Element(78,  'Pt', 'platinum'	, 0.19508490020544)
Au = Element(79,  'Au', 'gold'		, 0.19696657060743)
Hg = Element(80,  'Hg', 'mercury'	, 0.20059230021124)
Tl = Element(81,  'Tl', 'thallium'	, 0.20438100021523)
Pb = Element(82,  'Pb', 'lead'		, 0.2072000002182)
Bi = Element(83,  'Bi', 'bismuth'	, 0.20898040122008)

# Period 7
Th = Element(90,  'Th', 'thorium'	, 0.23203774024436)
Pa = Element(91,  'Pa', 'protactinium'	, 0.2310358812433)
U  = Element(92,  'U' , 'uranium'	, 0.23802891325067)

# A list of the elements, sorted by increasing atomic number
element_list = [
    e,
    X,
    H, D, T, He,
    Li, Be, B, C, C13, N, N15, O, O17, O18, F, Ne,
    Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
    Rb, Sr, Y, Zr, Nb, Mo, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe,
    Cs, Ba, La, Ce, Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi,
    Th, Pa, U
]

# Bond Dissociation Energies
# References:
# Huheey, pps. A-21 to A-34; T.L. Cottrell, "The Strengths of Chemical Bonds," 2nd ed., Butterworths, London, 1958
# B. deB. Darwent, "National Standard Reference Data Series," National Bureau of Standards, No. 31, Washington, DC, 1970
# S.W. Benson, J. Chem. Educ., 42, 502 (1965)
# (C,C,1.5) was taken from an unsourced table that had similar values to those used below, should be replaced
# if a sourced value becomes available
# (C,C,2.5) is C#C - (CbenzeneC - C-C)
# P=P value is from: https://www2.chemistry.msu.edu/faculty/reusch/OrgPage/bndenrgy.htm
# C#S is the value for [C+]#[S-] from 10.1002/chem.201002840 referenced relative to 0 K
# X-O and X-X (X=F,Cl,Br) taken from https://labs.chem.ucsb.edu/zakarian/armen/11---bonddissociationenergy.pdf
# Li-C and Li-S are taken referenced from 0K from from http://staff.ustc.edu.cn/~luo971/2010-91-CRC-BDEs-Tables.pdf
# Li-N is a G3 calculation taken from https://doi.org/10.1021/jp050857o
# The reference state is gaseous state at 298 K, but some of the values in the bde_dict might be coming from 0 K.
# The bond dissociation energy at 298 K is greater than the bond dissociation energy at 0 K by 0.6 to 0.9 kcal/mol
# (between RT and 3/2 RT), and this difference is usually much smaller than the uncertainty in the bond dissociation
# energy itself. Therefore, the discrepancy between 0 K and 298 K shouldn't matter too much.
# But for any new entries, try to use the consistent reference state of 298 K.
bde_elements = ['C', 'N', 'H', 'O', 'S', 'Cl', 'Si', 'P', 'F', 'Br', 'I', 'Li']  # elements supported by BDE
bde_dict = {('H', 'H', 1.0): (432.0, 'kJ/mol'), ('H', 'C', 1): (411.0, 'kJ/mol'),
            ('H', 'N', 1): (386.0, 'kJ/mol'), ('H', 'O', 1.0): (459.0, 'kJ/mol'),
            ('H', 'P', 1): (322.0, 'kJ/mol'), ('H', 'S', 1): (363.0, 'kJ/mol'),
            ('H', 'F', 1): (565.0, 'kJ/mol'), ('H', 'Cl', 1): (428.0, 'kJ/mol'),
            ('H', 'Br', 1): (362.0, 'kJ/mol'), ('H', 'I', 1): (295.0, 'kJ/mol'),
            ('C', 'C', 1): (346.0, 'kJ/mol'), ('C', 'C', 2): (602.0, 'kJ/mol'),
            ('C', 'C', 3): (835.0, 'kJ/mol'), ('C', 'C', 1.5): (518.0, 'kJ/mol'),
            ('C', 'C', 2.5): (663.0, 'kJ/mol'), ('C', 'Si', 1): (318.0, 'kJ/mol'),
            ('C', 'N', 1): (305.0, 'kJ/mol'), ('C', 'N', 2): (615.0, 'kJ/mol'),
            ('C', 'N', 3): (887.0, 'kJ/mol'), ('C', 'O', 1): (358.0, 'kJ/mol'),
            ('C', 'O', 2): (799.0, 'kJ/mol'), ('C', 'O', 3): (1072.0, 'kJ/mol'),
            ('C', 'P', 1): (264.0, 'kJ/mol'), ('C', 'S', 1): (272.0, 'kJ/mol'),
            ('C', 'S', 2): (573.0, 'kJ/mol'), ('C', 'S', 3): (709.6, 'kJ/mol'),
            ('C', 'F', 1): (485.0, 'kJ/mol'),
            ('C', 'Cl', 1): (327.0, 'kJ/mol'), ('C', 'Br', 1): (285.0, 'kJ/mol'),
            ('C', 'I', 1): (213.0, 'kJ/mol'),
            ('Si', 'Si', 1): (222.0, 'kJ/mol'), ('Si', 'N', 1): (355.0, 'kJ/mol'),
            ('Si', 'O', 1): (452.0, 'kJ/mol'), ('Si', 'S', 1): (293.0, 'kJ/mol'),
            ('Si', 'F', 1): (565.0, 'kJ/mol'), ('Si', 'Cl', 1): (381.0, 'kJ/mol'),
            ('Si', 'Br', 1): (310.0, 'kJ/mol'), ('Si', 'I', 1): (234.0, 'kJ/mol'),
            ('N', 'N', 1): (167.0, 'kJ/mol'), ('N', 'N', 2): (418.0, 'kJ/mol'),
            ('N', 'N', 3): (942.0, 'kJ/mol'), ('N', 'O', 1): (201.0, 'kJ/mol'),
            ('N', 'O', 2): (607.0, 'kJ/mol'), ('N', 'F', 1): (283.0, 'kJ/mol'),
            ('N', 'Cl', 1): (313.0, 'kJ/mol'), ('N', 'S', 1): (467.0, 'kJ/mol'),
            ('O', 'O', 1): (142.0, 'kJ/mol'), ('O', 'O', 2): (494.0, 'kJ/mol'),
            ('O', 'P', 1): (335.0, 'kJ/mol'), ('O', 'P', 2): (544.0, 'kJ/mol'),
            ('O', 'S', 1): (265.0, 'kJ/mol'), ('O', 'S', 2): (522.0, 'kJ/mol'),
            ('P', 'P', 1): (201.0, 'kJ/mol'), ('P', 'P', 2): (351.0, 'kJ/mol'),
            ('P', 'P', 3): (489.0, 'kJ/mol'), ('P', 'S', 2): (335.0, 'kJ/mol'),
            ('P', 'F', 1): (490.0, 'kJ/mol'), ('P', 'Cl', 1): (326.0, 'kJ/mol'),
            ('P', 'Br', 1): (264.0, 'kJ/mol'), ('P', 'I', 2): (184.0, 'kJ/mol'),
            ('S', 'S', 1): (226.0, 'kJ/mol'), ('S', 'S', 2): (425.0, 'kJ/mol'),
            ('S', 'Cl', 1): (255.0, 'kJ/mol'),
            ('F', 'F', 1): (155.0, 'kJ/mol'), ('Cl', 'Cl', 1): (240.0, 'kJ/mol'),
            ('Br', 'Br', 1): (190.0, 'kJ/mol'), ('I', 'I', 1): (148.0, 'kJ/mol'),
            ('F', 'O', 1): (222.0, 'kJ/mol'), ('Cl', 'O', 1): (272.0, 'kJ/mol'),
            ('Br', 'O', 1): (235.1, 'kJ/mol'), ('Cl', 'F', 1): (250.54, 'kJ/mol'),
            ('Br', 'F', 1): (233.8, 'kJ/mol'), ('Br', 'Cl', 1): (218.84, 'kJ/mol'),
            ('Li', 'Br', 1): (423.0, 'kJ/mol'), ('Li', 'F', 1): (577.0, 'kJ/mol'),
            ('Li', 'Cl', 1): (469.0, 'kJ/mol'), ('Li', 'H', 1): (247.0, 'kJ/mol'),
            ('Li', 'I', 1): (352.0, 'kJ/mol'), ('Li', 'O', 1): (341.6, 'kJ/mol'),
            ('Li', 'C', 1): (214.6, 'kJ/mol'), ('Li', 'S', 1): (312.5, 'kJ/mol'),
            ('Li', 'N', 1): (302.4, 'kJ/mol')}

bdes = {}
for key, value in bde_dict.items():
    q = Quantity(value).value_si
    bdes[(key[0], key[1], key[2])] = q
    bdes[(key[1], key[0], key[2])] = q
