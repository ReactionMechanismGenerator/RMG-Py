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
This creates atomTypes and assigns them with the correct bond/lone_pairs/charge
Used in isomorphismTest.py to create group_atomtypes
"""


class AbstractAtomType(object):
    def __init__(self, element=None, label=None, double=-1, triple=-1, quadruple=-1, benzene=-1, lp=-1, chrg=-1):
        self.element = element
        self.label = label
        self.double = double
        self.triple = triple
        self.quadruple = quadruple
        self.benzene = benzene
        self.lp = lp
        self.chrg = chrg


class Column4(AbstractAtomType):  # C
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.lp = 0


class Column5(AbstractAtomType):  # N
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.lp = 1


class Column6(AbstractAtomType):  # O, S
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.lp = 2


class Xs(AbstractAtomType):
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.double, self.triple, self.benzene, self.quadruple = 0, 0, 0, 0
        self.label = 's'


class Xd(AbstractAtomType):
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.double, self.triple, self.benzene, self.quadruple = 1, 0, 0, 0
        self.label = 'd'


class Xdd(AbstractAtomType):
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.double, self.triple, self.benzene, self.quadruple = 2, 0, 0, 0
        self.label = 'dd'


class Xt(AbstractAtomType):
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.double, self.triple, self.benzene, self.quadruple = 0, 1, 0, 0
        self.label = 't'


class Xq(AbstractAtomType):
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.double, self.triple, self.benzene, self.quadruple = 0, 0, 0, 1
        self.label = 'q'


class Xb(AbstractAtomType):
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.double, self.triple, self.benzene, self.quadruple = 0, 0, 2, 0
        self.label = 'b'


class Xbf(AbstractAtomType):
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.double, self.triple, self.benzene, self.quadruple = 0, 0, 3, 0
        self.label = 'bf'


def create_atom_types():
    atomtypes = []

    # tetravalent:

    # tetravalent = []
    # for type in [Xs, Xd, Xdd, Xt, Xb, Xbf]:
    #     #tetravalent.extend(create_types(type, ['C', 'Si']))
    #     tetravalent.extend(create_types(type, ['C']))

    # for at in tetravalent: at.lp = 0

    # atomtypes.extend(tetravalent)

    # bivalent:
    # bivalent = []
    # for type in [Xs, Xd]:
    #     #bivalent.extend(create_types(type, ['O', 'S']))
    #     bivalent.extend(create_types(type, ['O']))

    # for at in bivalent: at.lp = 2

    # atomtypes.extend(bivalent)

    # trivalent nitrogen:
    # trivalent_N = []
    # for type in [Xs, Xd, Xt, Xb]:
    #     trivalent_N.extend(create_types(type, ['N'], ['N3']))

    # for at in trivalent_N: at.lp = 1
    # atomtypes.extend(trivalent_N)

    # pentavalent nitrogen:
    # pentavalent_N = []
    # for type in [Xs, Xd, Xdd, Xt, Xb]:
    #     pentavalent_N.extend(create_types(type, ['N'], ['N5']))

    # for at in pentavalent_N: at.lp = 0
    # atomtypes.extend(pentavalent_N)

    return atomtypes


def create_types(Type, elements, labels=None):
    if labels is None:
        labels = elements
    atomtypes = []
    for el, label in zip(elements, labels):
        at = Type(element=el)
        at.label = label + at.label
        atomtypes.append(at)
    return atomtypes
