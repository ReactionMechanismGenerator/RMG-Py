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

"""
This module defines the atom types that are available for representing
molecular functional groups and substructure patterns. Each available atom type
is defined as an instance of the :class:`AtomType` class. The atom types 
themselves are available in the ``ATOMTYPES`` module-level variable, or as
the return value from the :meth:`get_atomtype()` method.

If you want to change which atom types are available in RMG and/or what they
represent, this should be the only module you need to change to do so.
"""

import cython
from rmgpy.exceptions import AtomTypeError

################################################################################


class AtomType:
    """
    A class for internal representation of atom types. Using unique objects
    rather than strings allows us to use fast pointer comparisons instead of
    slow string comparisons, as well as store extra metadata. In particular,
    we store metadata describing the atom type's hierarchy with regard to other
    atom types, and the atom types that can result when various actions
    involving this atom type are taken. The attributes are:

    ===================== =================== ====================================
    Attribute             Type                Description
    ===================== =================== ====================================
    `label`               ``str``             A unique label for the atom type
    `generic`             ``list``            The atom types that are more generic than this one
    `specific`            ``list``            The atom types that are more specific than this one
    `increment_bond`      ``list``            The atom type(s) that result when an adjacent bond's order is incremented
    `decrement_bond`      ``list``            The atom type(s) that result when an adjacent bond's order is decremented
    `form_bond`           ``list``            The atom type(s) that result when a new single bond is formed to this atom type
    `break_bond`          ``list``            The atom type(s) that result when an existing single bond to this atom type is broken
    `increment_radical`   ``list``            The atom type(s) that result when the number of radical electrons is incremented
    `decrement_radical`   ``list``            The atom type(s) that result when the number of radical electrons is decremented
    `increment_lone_pair` ``list``            The atom type(s) that result when the number of lone electron pairs is incremented
    `decrement_lone_pair` ``list``            The atom type(s) that result when the number of lone electron pairs is decremented

    The following features are what are required in a given atomtype. Any int in the list is acceptable. An empty list is a wildcard
    ------------------------------------------------------------------------------
    `single`              ``list``            The total number of single bonds on the atom
    `all_double`          ``list``            The total number of double bonds on the atom
    `r_double`            ``list``            The number of double bonds to any non-oxygen, nonsulfur
    `o_double`            ``list``            The number of double bonds to oxygen
    `s_double`            ``list``            The number of double bonds to sulfur
    `triple`              ``list``            The total number of triple bonds on the atom
    `quadruple`           ``list``            The total number of quadruple bonds on the atom
    `benzene`             ``list``            The total number of benzene bonds on the atom
    `lone_pairs`          ``list``            The number of lone pairs on the atom
    `charge`              ``list``            The partial charge of the atom
    ===================== =================== ====================================

    """

    def __init__(self, label='',
                 generic=None,
                 specific=None,
                 single=None,
                 all_double=None,
                 r_double=None,
                 o_double=None,
                 s_double=None,
                 triple=None,
                 quadruple=None,
                 benzene=None,
                 lone_pairs=None,
                 charge=None):
        self.label = label
        self.generic = generic or []
        self.specific = specific or []
        self.increment_bond = []
        self.decrement_bond = []
        self.form_bond = []
        self.break_bond = []
        self.increment_radical = []
        self.decrement_radical = []
        self.increment_lone_pair = []
        self.decrement_lone_pair = []
        self.single = single or []
        self.all_double = all_double or []
        self.r_double = r_double or []
        self.o_double = o_double or []
        self.s_double = s_double or []
        self.triple = triple or []
        self.quadruple = quadruple or []
        self.benzene = benzene or []
        self.lone_pairs = lone_pairs or []
        self.charge = charge or []

    def __repr__(self):
        return '<AtomType "%s">' % self.label

    def __reduce__(self):
        """
        A helper function used when pickling an AtomType object.
        """
        d = {
            'label': self.label,
            'generic': self.generic,
            'specific': self.specific,
            'increment_bond': self.increment_bond,
            'decrement_bond': self.decrement_bond,
            'form_bond': self.form_bond,
            'break_bond': self.break_bond,
            'increment_radical': self.increment_radical,
            'decrement_radical': self.decrement_radical,
            'increment_lone_pair': self.increment_lone_pair,
            'decrement_lone_pair': self.decrement_lone_pair,
            'single': self.single,
            'all_double': self.all_double,
            'r_double': self.r_double,
            'o_double': self.o_double,
            's_double': self.s_double,
            'triple': self.triple,
            'quadruple': self.quadruple,
            'benzene': self.benzene,
            'lone_pairs': self.lone_pairs,
            'charge': self.charge
        }
        return AtomType, (), d

    def __setstate__(self, d):
        """
        A helper function used when unpickling an AtomType object.
        """
        self.label = d['label']
        self.generic = d['generic']
        self.specific = d['specific']
        self.increment_bond = d['increment_bond']
        self.decrement_bond = d['decrement_bond']
        self.form_bond = d['form_bond']
        self.break_bond = d['break_bond']
        self.increment_radical = d['increment_radical']
        self.decrement_radical = d['decrement_radical']
        self.increment_lone_pair = d['increment_lone_pair']
        self.decrement_lone_pair = d['decrement_lone_pair']
        self.single = d['single']
        self.all_double = d['all_double']
        self.r_double = d['r_double']
        self.o_double = d['o_double']
        self.s_double = d['s_double']
        self.triple = d['triple']
        self.quadruple = d['quadruple']
        self.benzene = d['benzene']
        self.lone_pairs = d['lone_pairs']
        self.charge = d['charge']

    def set_actions(self, increment_bond, decrement_bond, form_bond, break_bond, increment_radical, decrement_radical,
                    increment_lone_pair, decrement_lone_pair):
        self.increment_bond = increment_bond
        self.decrement_bond = decrement_bond
        self.form_bond = form_bond
        self.break_bond = break_bond
        self.increment_radical = increment_radical
        self.decrement_radical = decrement_radical
        self.increment_lone_pair = increment_lone_pair
        self.decrement_lone_pair = decrement_lone_pair

    def equivalent(self, other):
        """
        Returns ``True`` if two atom types `atomType1` and `atomType2` are
        equivalent or ``False``  otherwise. This function respects wildcards,
        e.g. ``R!H`` is equivalent to ``C``.
        """
        return self is other or self in other.specific or other in self.specific

    def is_specific_case_of(self, other):
        """
        Returns ``True`` if atom type `atomType1` is a specific case of
        atom type `atomType2` or ``False``  otherwise.
        """
        return self is other or self in other.specific
    
    def get_features(self):
        """
        Returns a list of the features that are checked to determine atomtype
        """
        features = [self.single,
                    self.all_double,
                    self.r_double,
                    self.o_double,
                    self.s_double,
                    self.triple,
                    self.quadruple,
                    self.benzene,
                    self.lone_pairs,
                    self.charge]
        return features

################################################################################


"""
Note: function to read adjacency lists assumes that all atom types begin
with a capital letter [A-Z]

For making sample atoms, we use the first atomtype under specific,
therefore the first one in the list should always be an element.

The atomTypes naming convention is:
<element> <valence> <characteristic bonds> <charge(optional)>
For example:
- N3d is nitrogen with valence=3 (i.e., 3 electrons are able to form bonds or remain as radicals) with one double bond
- S2tc is a charged sulfur with valence=2 with a triple bonds
- Oa is atomic oxygen, i.e., a closed shell atom
Some charged atom types were merged together, and are marked as '*Composite atomtype'
"""

ATOMTYPES = {}

# Surface sites:
ATOMTYPES['X']   = AtomType(label='X', generic=[], specific=['Xv', 'Xo'])

# Vacant surface site:
ATOMTYPES['Xv']   = AtomType('Xv', generic=['X'], specific=[],
                             single=[0], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[0],
                             benzene=[0], lone_pairs=[0])
# Occupied surface site:
ATOMTYPES['Xo']   = AtomType('Xo', generic=['X'], specific=[],
                             single=[0, 1], all_double=[0, 1], r_double=[], o_double=[], s_double=[], triple=[0, 1],
                             quadruple=[0, 1], benzene=[0], lone_pairs=[0])

# Non-surface atomTypes, R being the most generic:
ATOMTYPES['R']    = AtomType(label='R', generic=[], specific=[
    'H',
    'R!H',
    'Val4','Val5','Val6','Val7',    
    'He','Ne','Ar',
    'C','Ca','Cs','Csc','Cd','CO','CS','Cdd','Cdc','Ct','Cb','Cbf','Cq','C2s','C2sc','C2d','C2dc','C2tc',
    'N','N0sc','N1s','N1sc','N1dc','N3s','N3sc','N3d','N3t','N3b','N5sc','N5dc','N5ddc','N5dddc','N5tc','N5b','N5bd',
    'O','Oa','O0sc','O2s','O2sc','O2d','O4sc','O4dc','O4tc','O4b',
    'Si','Sis','Sid','Sidd','Sit','SiO','Sib','Sibf','Siq',
    'P','P0sc','P1s','P1sc','P1dc','P3s','P3d','P3t','P3b','P5s','P5sc','P5d','P5dd','P5dc','P5ddc','P5t','P5td','P5tc','P5b','P5bd',
    'S','Sa','S0sc','S2s','S2sc','S2d','S2dc','S2tc','S4s','S4sc','S4d','S4dd','S4dc','S4b','S4t','S4tdc','S6s','S6sc','S6d','S6dd','S6ddd','S6dc','S6t','S6td','S6tt','S6tdc',
    'Cl','Cl1s',
    'Br','Br1s',
    'I','I1s',
    'F','F1s'])

ATOMTYPES['R!H']  = AtomType(label='R!H', generic=['R'], specific=[
    'Val4','Val5','Val6','Val7',
    'He','Ne','Ar',
    'C','Ca','Cs','Csc','Cd','CO','CS','Cdd','Cdc','Ct','Cb','Cbf','Cq','C2s','C2sc','C2d','C2dc','C2tc',
    'N','N0sc','N1s','N1sc','N1dc','N3s','N3sc','N3d','N3t','N3b','N5sc','N5dc','N5ddc','N5dddc','N5tc','N5b','N5bd',
    'O','Oa','O0sc','O2s','O2sc','O2d','O4sc','O4dc','O4tc','O4b',
    'Si','Sis','Sid','Sidd','Sit','SiO','Sib','Sibf','Siq',
    'P','P0sc','P1s','P1sc','P1dc','P3s','P3d','P3t','P3b','P5s','P5sc','P5d','P5dd','P5dc','P5ddc','P5t','P5td','P5tc','P5b','P5bd',
    'S','Sa','S0sc','S2s','S2sc','S2d','S2dc','S2tc','S4s','S4sc','S4d','S4dd','S4dc','S4b','S4t','S4tdc','S6s','S6sc','S6d','S6dd','S6ddd','S6dc','S6t','S6td','S6tt','S6tdc',
    'Cl','Cl1s',
    'Br','Br1s',
    'I','I1s',
    'F','F1s'])

ATOMTYPES['Val4'] = AtomType(label='Val4', generic=['R', 'R!H'], specific=[
    'C','Ca','Cs','Csc','Cd','CO','CS','Cq','Cdd','Cdc','Ct','Cb','Cbf','C2s','C2sc','C2d','C2dc','C2tc',
    'Si','Sis','Sid','Sidd','Sit','SiO','Sib','Sibf','Siq'])

ATOMTYPES['Val5'] = AtomType(label='Val5', generic=['R', 'R!H'], specific=[
    'N','N0sc','N1s','N1sc','N1dc','N3s','N3sc','N3d','N3t','N3b','N5sc','N5dc','N5ddc','N5dddc','N5tc','N5b','N5bd',
    'P','P0sc','P1s','P1sc','P1dc','P3s','P3d','P3t','P3b','P5s','P5sc','P5d','P5dd','P5dc','P5ddc','P5t','P5td','P5tc','P5b','P5bd'])

ATOMTYPES['Val6'] = AtomType(label='Val6', generic=['R', 'R!H'], specific=[
    'O','Oa','O0sc','O2s','O2sc','O2d','O4sc','O4dc','O4tc','O4b',
    'S','Sa','S0sc','S2s','S2sc','S2d','S2dc','S2tc','S4s','S4sc','S4d','S4dd','S4dc','S4b','S4t','S4tdc','S6s','S6sc','S6d','S6dd','S6ddd','S6dc','S6t','S6td','S6tt','S6tdc'])

ATOMTYPES['Val7'] = AtomType(label='Val7', generic=['R', 'R!H'], specific=[
    'Cl','Cl1s',
    'Br','Br1s',
    'I','I1s',
    'F','F1s'])

ATOMTYPES['H'] = AtomType('H', generic=['R'], specific=[])

ATOMTYPES['He'] = AtomType('He', generic=['R', 'R!H'], specific=[])
ATOMTYPES['Ne'] = AtomType('Ne', generic=['R', 'R!H'], specific=[])
ATOMTYPES['Ar'] = AtomType('Ar', generic=['R', 'R!H'], specific=[])

ATOMTYPES['C'] = AtomType('C', generic=['R', 'R!H', 'Val4'], specific=['Ca', 'Cs', 'Csc', 'Cd', 'CO', 'Cq', 'CS', 'Cdd', 'Cdc', 'Ct', 'Cb', 'Cbf', 'C2s', 'C2sc', 'C2d', 'C2dc', 'C2tc'],
                          single=[], all_double=[], r_double=[], o_double=[], s_double=[], triple=[], quadruple=[], benzene=[], lone_pairs=[], charge=[])  # todo: double check to see if quadruple should be blank or 0 for all of these as well as being 1 for quadruple
ATOMTYPES['Ca'] = AtomType('Ca', generic=['R', 'R!H', 'C', 'Val4'], specific=[],  # (shared electrons = 4)
                           single=[0], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[2], charge=[0])
# examples for Ca: atomic carbon (closed shell)
ATOMTYPES['Cs'] = AtomType('Cs', generic=['R', 'R!H', 'C', 'Val4'], specific=[],  # (shared electrons = 4-8)
                           single=[0,1,2,3,4], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[0])
# examples for Cs: C, CC,
ATOMTYPES['Csc'] = AtomType('Csc', generic=['R', 'R!H', 'C', 'Val4'], specific=[],  # (shared electrons = 3-6)
                            single=[0,1,2,3], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[+1])
# examples for Csc: C1=CCC([O-])[CH+]1, O[O+]=C[C+]C([O-])[O-]
ATOMTYPES['Cd'] = AtomType('Cd', generic=['R', 'R!H', 'C', 'Val4'], specific=[],  # (shared electrons = 6-8)
                           single=[0,1,2], all_double=[1], r_double=[1], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[0])
# examples for Cd: C=C, C=N
ATOMTYPES['Cdc'] = AtomType('Cdc', generic=['R', 'R!H', 'C', 'Val4'], specific=[],  # (shared electrons = 6)
                            single=[0,1], all_double=[1], r_double=[0, 1], o_double=[0, 1], s_double=[0, 1], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[+1])
# examples for Cdc: [CH+]=C=[CH-], [CH+]=N[O-] (one of the res structures of Fulminic acid)
ATOMTYPES['CO'] = AtomType('CO', generic=['R', 'R!H', 'C', 'Val4'], specific=[],  # (shared electrons = 6-8)
                           single=[0,1,2], all_double=[1], r_double=[0], o_double=[1], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[0])
# examples for CO: C=O
ATOMTYPES['CS'] = AtomType('CS', generic=['R', 'R!H', 'C', 'Val4'], specific=[],  # (shared electrons = 6-8)
                           single=[0,1,2], all_double=[1], r_double=[0], o_double=[0], s_double=[1], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[0])
# examples for CS: C=S
ATOMTYPES['Cdd'] = AtomType('Cdd', generic=['R', 'R!H', 'C', 'Val4'], specific=[],  # (shared electrons = 8)
                            single=[0], all_double=[2], r_double=[0, 1, 2], o_double=[0, 1, 2], s_double=[0, 1, 2], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[0])
# examples for Cdd: O=C=O, C=C=C
ATOMTYPES['Ct'] = AtomType('Ct', generic=['R', 'R!H', 'C', 'Val4'], specific=[],  # (shared electrons = 7-8)
                           single=[0,1], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[1], quadruple=[], benzene=[0], lone_pairs=[0], charge=[0])
# examples for Ct: C#C, C#N
ATOMTYPES['Cb'] = AtomType('Cb', generic=['R', 'R!H', 'C', 'Val4'], specific=[],
                           single=[0,1], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[1, 2], lone_pairs=[], charge=[])
# examples for Cb: benzene (C6H6)
ATOMTYPES['Cbf'] = AtomType('Cbf', generic=['R', 'R!H', 'C', 'Val4'], specific=[],
                            single=[0], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[3], lone_pairs=[], charge=[])
# examples for Cbf: Naphthalene
ATOMTYPES['Cq'] = AtomType('Cq', generic=['R', 'R!H', 'C', 'Val4'], specific=[],
                           single=[0], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[1], benzene=[0], lone_pairs=[], charge=[])
# examples for Cq: C2
ATOMTYPES['C2s'] = AtomType('C2s', generic=['R', 'R!H', 'C', 'Val4'], specific=[],  # (shared electrons = 4-6)
                            single=[0,1,2], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[0])
# examples for C2s: singlet[CH2]
ATOMTYPES['C2sc'] = AtomType('C2sc', generic=['R', 'R!H', 'C', 'Val4'], specific=[],  # (shared electrons = 5-8)
                             single=[0,1,2,3], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[-1])
# examples for C2sc: [CH2-][N+]#N
ATOMTYPES['C2d'] = AtomType('C2d', generic=['R', 'R!H', 'C', 'Val4'], specific=[],  # (shared electrons = 6)
                            single=[0], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[0])
# examples for C2d: singlet[C]=C
ATOMTYPES['C2dc'] = AtomType('C2dc', generic=['R', 'R!H', 'C', 'Val4'], specific=[],  # (shared electrons = 7-8)
                             single=[0,1], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[-1])
# examples for C2dc: C=[C-][N+]#N, [CH-]=[N+]=O, [CH+]=C=[CH-]
ATOMTYPES['C2tc'] = AtomType('C2tc', generic=['R', 'R!H', 'C', 'Val4'], specific=[],  # (shared electrons = 8)
                             single=[0], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[1], quadruple=[], benzene=[0], lone_pairs=[1], charge=[-1])
# examples for C2tc: [C-]#[O+], H[N+]#[C-]

ATOMTYPES['N'] = AtomType('N', generic=['R', 'R!H', 'Val5'], specific=['N0sc', 'N1s', 'N1sc', 'N1dc', 'N3s', 'N3sc', 'N3d', 'N3t', 'N3b', 'N5sc', 'N5dc', 'N5ddc', 'N5dddc', 'N5tc', 'N5b', 'N5bd'],
                          single=[], all_double=[], r_double=[], o_double=[], s_double=[], triple=[], quadruple=[], benzene=[], lone_pairs=[], charge=[])
ATOMTYPES['N0sc'] = AtomType('N0sc', generic=['R', 'R!H', 'N', 'Val5'], specific=[],  # (shared electrons = 7-8)
                             single=[0,1], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[3], charge=[-2])
# examples for N0sc: [NH+]#[N+][N-2] with adjList 1 N u0 p0 c+1 {2,S} {3,T}; 2 H u0 p0 c0 {1,S}; 3 N u0 p0 c+1 {1,T} {4,S}; 4 N u0 p3 c-2 {3,S}
ATOMTYPES['N1s'] = AtomType('N1s', generic=['R', 'R!H', 'N', 'Val5'], specific=[],  # (shared electrons = 5-6)
                            single=[0,1], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[2], charge=[0])
# examples for N1s: closed shell N-N, closed shell NH
ATOMTYPES['N1sc'] = AtomType('N1sc', generic=['R', 'R!H', 'N', 'Val5'], specific=[],  # (shared electrons = 6-8)
                             single=[0,1,2], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[2], charge=[-1])
# examples for N1sc: [NH-][S+]=C, [NH-][N+]#C
ATOMTYPES['N1dc'] = AtomType('N1dc', generic=['R', 'R!H', 'N', 'Val5'], specific=[],  # (shared electrons = 8)
                             single=[0], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[2], charge=[-1])
# examples for N1dc: [N-]=[N+]=N terminal nitrogen on azide (two lone pairs), [N-]=[NH+], [N-]=[SH+]
ATOMTYPES['N3s'] = AtomType('N3s', generic=['R', 'R!H', 'N', 'Val5'], specific=[],  # (shared electrons = 5-8)
                            single=[0,1,2,3], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[0])
# examples for N3s: NH3, NH2, NH, N, C[NH]...
ATOMTYPES['N3sc'] = AtomType('N3sc', generic=['R', 'R!H', 'N', 'Val5'], specific=[],  # (shared electrons = 4-6)
                             single=[0,1,2], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[+1])
# examples for N3sc: !! N3sc should eventually be deleted, see #1206
ATOMTYPES['N3d'] = AtomType('N3d', generic=['R', 'R!H', 'N', 'Val5'], specific=[],  # (shared electrons = 7-8)
                            single=[0,1], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[0])
# examples for N3d: N=O, N=N, C=N, [O]N=O, [N]=O, [N]=C
ATOMTYPES['N3t'] = AtomType('N3t', generic=['R', 'R!H', 'N', 'Val5'], specific=[],  # (shared electrons = 8)
                            single=[0], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[1], quadruple=[], benzene=[0], lone_pairs=[1], charge=[0])
# examples for N3t: N2, N#C, N#[C], N#CC
ATOMTYPES['N3b'] = AtomType('N3b', generic=['R', 'R!H', 'N', 'Val5'], specific=[],
                            single=[0], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[2], lone_pairs=[1], charge=[0])
# examples for N3b: Oxazole, Pyradine, Pyrazine, 1,3,5-Triazine, Benzimidazole, Purine
ATOMTYPES['N5sc'] = AtomType('N5sc', generic=['R', 'R!H', 'N', 'Val5'], specific=[],  # (shared electrons = 4-8)
                             single=[0,1,2,3,4], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], benzene=[0], lone_pairs=[0], charge=[+1, +2])
# examples for N5sc: [NH3+][O-]
ATOMTYPES['N5dc'] = AtomType('N5dc', generic=['R', 'R!H', 'N', 'Val5'], specific=[],  # (shared electrons = 6-8)
                             single=[0,1,2], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[+1])
# examples for N5dc: O[N+](=O)(O-) nitrate group, [N+](=O)(O)[O-], O=[N+][O-], [N+](=O)(O[N+](=O)[O-])[O-], C=[N+]=[SH-], [NH2+]=[SH-]
ATOMTYPES['N5ddc'] = AtomType('N5ddc', generic=['R', 'R!H', 'N', 'Val5'], specific=[],  # (shared electrons = 8)
                              single=[0], all_double=[2], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[+1])
# examples for N5ddc: N=[N+]=[N-] center nitrogen on azide, [N-]=[N+]=O, C=[N+]=[SH-]
ATOMTYPES['N5dddc'] = AtomType('N5dddc', generic=['R', 'R!H', 'N', 'Val5'], specific=[],  # (shared electrons = 6)
                               single=[0], all_double=[3], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[-1])
# examples for N5dddc: C=[N-](=C)=[NH2+]
ATOMTYPES['N5tc'] = AtomType('N5tc', generic=['R', 'R!H', 'N', 'Val5'], specific=[],  # (shared electrons = 7-8)
                             single=[0,1], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[1], quadruple=[], benzene=[0], lone_pairs=[0], charge=[+1])
# examples for N5tc: C[N+]#[C-] isocyano group, N#[N+][O-], [NH+]#[C-] (note that C- has p1 here), [N+]#[C-] (note that C- has p1 here), [O-][N+]#C (one of the res structures of Fulminic acid), C[N+]#[C-] (note that C- has p1 here)
ATOMTYPES['N5b'] = AtomType('N5b', generic=['R', 'R!H', 'N', 'Val5'], specific=[],
                            single=[0,1], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[2], lone_pairs=[0], charge=[0, +1])
# examples for N5b: Pyrrole, Indole, Benzimidazole, Purine; Note that this is the only N atomtype with valence 5 which isn't necessarily charged.
ATOMTYPES['N5bd'] = AtomType('N5bd', generic=['R', 'R!H', 'N', 'Val5'], specific=[],
                             single=[0], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[2], lone_pairs=[0], charge=[0])
# examples for N5bd: AdjList """1 N u0 p0 c0 {2,B} {6,B} {7,D} 2 C u0 p0 {1,B} {3,B} {8,S} 3 C u0 p0 {2,B} {4,B} {9,S} 4 C u0 p0 {3,B} {5,B} {10,S} 5 C u0 p0 {4,B} {6,B} {11,S} 6 N u0 p1 {1,B} {5,B} 7 O u0 p2 c0 {1,D} 8 H u0 p0 {2,S} 9 H u0 p0 {3,S} 10 H u0 p0 {4,S} 11 H u0 p0 {5,S}"""

ATOMTYPES['O'] = AtomType('O', generic=['R', 'R!H', 'Val6'], specific=['Oa', 'O0sc', 'O2s', 'O2sc', 'O2d', 'O4sc', 'O4dc', 'O4tc', 'O4b'],
                          single=[], all_double=[], r_double=[], o_double=[], s_double=[], triple=[], quadruple=[], benzene=[], lone_pairs=[], charge=[])
ATOMTYPES['Oa'] = AtomType('Oa', generic=['R', 'R!H', 'O', 'Val6'], specific=[],  # (shared electrons = 6)
                           single=[0], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[0], benzene=[0], lone_pairs=[3], charge=[0])
# examples for Oa: atomic oxygen (closed shell)
ATOMTYPES['O0sc'] = AtomType('O0sc', generic=['R', 'R!H', 'O', 'Val6'], specific=[],  # (shared electrons = 8)
                             single=[0,1], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[3], charge=[-1])
# examples for O0sc: Nitric acid O[N+](=O)([O-])
ATOMTYPES['O2s'] = AtomType('O2s', generic=['R', 'R!H', 'O', 'Val6'], specific=[],  # (shared electrons = 8)
                            single=[0,1,2], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[2], charge=[0])
# examples for O2s: H2O, OH, CH3OH
ATOMTYPES['O2sc'] = AtomType('O2sc', generic=['R', 'R!H', 'O', 'Val6'], specific=[],  # (shared electrons = 6)
                             single=[0,1], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[2], charge=[+1])
# examples for O2sc: C=[S-][O+]
ATOMTYPES['O2d'] = AtomType('O2d', generic=['R', 'R!H', 'O', 'Val6'], specific=[],  # (shared electrons = 8)
                            single=[0], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[2], charge=[0])
# examples for O2d: CO2, CH2O
ATOMTYPES['O4sc'] = AtomType('O4sc', generic=['R', 'R!H', 'O', 'Val6'], specific=[],  # (shared electrons = 5-8)
                             single=[0,1,2,3], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[+1])
# examples for O4sc: [O-][OH+]C
ATOMTYPES['O4dc'] = AtomType('O4dc', generic=['R', 'R!H', 'O', 'Val6'], specific=[],  # (shared electrons = 7-8)
                             single=[0,1], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[+1])
# examples for O4dc: the positively charged O in ozone [O-][O+]=O
ATOMTYPES['O4tc'] = AtomType('O4tc', generic=['R', 'R!H', 'O', 'Val6'], specific=[],  # (shared electrons = 8)
                             single=[0], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[1], quadruple=[], benzene=[0], lone_pairs=[1], charge=[+1])
# examples for O4tc: [C-]#[O+]
ATOMTYPES['O4b'] = AtomType('O4b', generic=['R', 'R!H', 'O', 'Val6'], specific=[],
                            single=[0], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[2], lone_pairs=[1], charge=[0])
# examples for O4b: Furane, Benzofurane, Oxazole...

ATOMTYPES['Ne'] = AtomType('Ne', generic=['R', 'R!H'], specific=[])
ATOMTYPES['Si'] = AtomType('Si', generic=['R', 'R!H', 'Val4'], specific=['Sis', 'Sid', 'Sidd', 'Sit', 'SiO', 'Sib', 'Sibf', 'Siq'],
                           single=[], all_double=[], r_double=[], o_double=[], s_double=[], triple=[], quadruple=[], benzene=[], lone_pairs=[], charge=[])
ATOMTYPES['Sis'] = AtomType('Sis', generic=['R', 'R!H', 'Si', 'Val4'], specific=[],
                            single=[], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[], charge=[])
ATOMTYPES['SiO'] = AtomType('SiO', generic=['R', 'R!H', 'Si', 'Val4'], specific=[],
                            single=[], all_double=[1], r_double=[], o_double=[1], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[], charge=[])
ATOMTYPES['Sid'] = AtomType('Sid', generic=['R', 'R!H', 'Si', 'Val4'], specific=[],
                            single=[], all_double=[1], r_double=[], o_double=[0], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[], charge=[])
ATOMTYPES['Sidd'] = AtomType('Sidd', generic=['R', 'R!H', 'Si', 'Val4'], specific=[],
                             single=[], all_double=[2], r_double=[0, 1, 2], o_double=[0, 1, 2], s_double=[0, 1, 2], triple=[0], quadruple=[], benzene=[0], lone_pairs=[], charge=[])
ATOMTYPES['Sit'] = AtomType('Sit', generic=['R', 'R!H', 'Si', 'Val4'], specific=[],
                            single=[], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[1], quadruple=[], benzene=[0], lone_pairs=[], charge=[])
ATOMTYPES['Sib'] = AtomType('Sib', generic=['R', 'R!H', 'Si', 'Val4'], specific=[],
                            single=[], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[2], lone_pairs=[], charge=[])
ATOMTYPES['Sibf'] = AtomType('Sibf', generic=['R', 'R!H', 'Si', 'Val4'], specific=[],
                             single=[], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[3], lone_pairs=[], charge=[])
ATOMTYPES['Siq'] = AtomType('Siq', generic=['R', 'R!H', 'Si', 'Val4'], specific=[],
                            single=[0], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[1], benzene=[0], lone_pairs=[], charge=[])

ATOMTYPES['P'] = AtomType('P', generic=['R', 'R!H', 'Val5'], specific=['P0sc', 'P1s', 'P1sc', 'P1dc', 'P3s', 'P3d', 'P3t', 'P3b', 'P5s', 'P5sc', 'P5d', 'P5dd', 'P5dc', 'P5ddc', 'P5t', 'P5td', 'P5tc', 'P5b', 'P5bd'],
                          single=[], all_double=[], r_double=[], o_double=[], s_double=[], triple=[], quadruple=[], benzene=[], lone_pairs=[], charge=[])
ATOMTYPES['P0sc'] = AtomType('P0sc', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                             single=[0,1], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[3], charge=[-2])
# examples for P0sc: [PH-2] (Phosphanediide), [P-2][P+]#[PH+] with adjList '''1 P u0 p3 c-2 {2,S}  2 P u0 p0 c+1 {1,S} {3,T} 3 P u0 p0 c+1 {2,T} {4,S}  4 H u0 p0 c0 {3,S}'''
ATOMTYPES['P1s'] = AtomType('P1s', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                            single=[0,1], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[2], charge=[0])
# examples for P1s: closed shell [PH] (Phosphinidene)
ATOMTYPES['P1sc'] = AtomType('P1sc', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                             single=[0,1,2], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[2], charge=[-1])
# examples for P1sc: C[PH-] (methylphosphanide)
ATOMTYPES['P1dc'] = AtomType('P1dc', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                             single=[0], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[2], charge=[-1])
# examples for P1dc: C=[P-] (methylidenephosphanide)
ATOMTYPES['P3s'] = AtomType('P3s', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                            single=[0,1,2,3], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[0])
# examples for P3s: PH3, PCl3
ATOMTYPES['P3d'] = AtomType('P3d', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                            single=[0,1], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[0])
# examples for P3d: O=[PH] with adjList '''1 O u0 p2 c0 {2,D} 2 P u0 p1 c0 {1,D} {3,S} 3 H u0 p0 c0 {2,S}'''
ATOMTYPES['P3t'] = AtomType('P3t', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                            single=[0], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[1], quadruple=[], benzene=[0], lone_pairs=[1], charge=[0])
# examples for P3t: P#P (diphosphorus)
ATOMTYPES['P3b'] = AtomType('P3b', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                            single=[0], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[2], lone_pairs=[1], charge=[0])
# examples for P3b: c1ccpcc1 (phosphorine) with InChI 'InChI=1S/C5H5P/c1-2-4-6-5-3-1/h1-5H'
ATOMTYPES['P5s'] = AtomType('P5s', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                            single=[0,1,2,3,4,5], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[0], benzene=[0], lone_pairs=[0], charge=[0])
# examples for P5s: P(Cl)(Cl)(Cl)(Cl)Cl (phosphorus pentachloride)
ATOMTYPES['P5sc'] = AtomType('P5sc', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                             single=[0,1,2,3,4,5,6], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], benzene=[0], lone_pairs=[0], charge=[-1, +1, +2])
# examples for P5sc: [O-][PH3+] (oxidophosphanium), F[P-](F)(F)(F)(F)F (Hexafluorophosphate)
ATOMTYPES['P5d'] = AtomType('P5d', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                            single=[0,1,2,3], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[0], benzene=[0], lone_pairs=[0], charge=[0])
# examples for P5d: OP(=O)(O)O (phosphoric acid)
ATOMTYPES['P5dd'] = AtomType('P5dd', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                            single=[0,1], all_double=[2], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[0], benzene=[0], lone_pairs=[0], charge=[0])
# examples for P5dd: CP(=O)=O (methylphosphinate)
ATOMTYPES['P5dc'] = AtomType('P5dc', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                             single=[0,1,2], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[+1])
# examples for P5dc: C=C[P+](=N)[O-] (ethenyl-imino-oxidophosphanium), C[P+](=C)C (methylenedimethylphosphorane)
ATOMTYPES['P5ddc'] = AtomType('P5ddc', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                              single=[0], all_double=[2], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[+1])
# examples for P5ddc: C=[P+]=N (imino(methylidene)phosphanium)
ATOMTYPES['P5t'] = AtomType('P5t', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                             single=[0,1,2], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[1], quadruple=[], benzene=[0], lone_pairs=[0], charge=[0])
# examples for P5t: N#P(Cl)Cl (phosphonitrile chloride)
ATOMTYPES['P5td'] = AtomType('P5td', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                             single=[0], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[1], quadruple=[], benzene=[0], lone_pairs=[0], charge=[0])
# examples for P5td: COC(=O)C#P=O (methyl phosphorylacetate)
ATOMTYPES['P5tc'] = AtomType('P5tc', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                             single=[0,1], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[1], quadruple=[], benzene=[0], lone_pairs=[0], charge=[+1])
# examples for P5tc: C[P+]#C (methyl(methylidyne)phosphanium), C#[P+]O (hydroxy(methylidyne)phosphanium)
ATOMTYPES['P5b'] = AtomType('P5b', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                            single=[0,1], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[2], lone_pairs=[0], charge=[0, +1])
# examples for P5b: C1=CC=[PH+]C=C1 (Phosphoniabenzene)
ATOMTYPES['P5bd'] = AtomType('P5bd', generic=['R', 'R!H', 'P', 'Val5'], specific=[],
                             single=[0], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[2], lone_pairs=[0], charge=[0])
# examples for P5bd: C1=CC=P(=S)C=C1 (Phosphorin 1-sulfide), C1=CC=P(=O)C=C1 (Phosphorin 1-oxide)

ATOMTYPES['S'] = AtomType('S', generic=['R', 'R!H', 'Val6'], specific=['Sa', 'S0sc', 'S2s', 'S2sc', 'S2d', 'S2dc', 'S2tc', 'S4s', 'S4sc', 'S4d', 'S4dd', 'S4dc', 'S4b', 'S4t', 'S4tdc', 'S6s', 'S6sc', 'S6d', 'S6dd', 'S6ddd', 'S6dc', 'S6t', 'S6td', 'S6tt', 'S6tdc'],
                          single=[], all_double=[], r_double=[], o_double=[], s_double=[], triple=[], quadruple=[], benzene=[], lone_pairs=[], charge=[])
ATOMTYPES['Sa'] = AtomType('Sa', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 6)
                           single=[0], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[3], charge=[0])
# examples for Sa: atomic sulfur (closed shell)
ATOMTYPES['S0sc'] = AtomType('S0sc', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 7-8)
                             single=[0,1], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[3], charge=[-1])
# examples for S0sc: [S-][S+]=S
ATOMTYPES['S2s'] = AtomType('S2s', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 6-8)
                            single=[0,1,2], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[2], charge=[0])
# examples for S2s: [S], [SH], S {H2S}, [S][S], SS {H2S2}, SSC, CSSC, SO {HSOH}...
ATOMTYPES['S2sc'] = AtomType('S2sc', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 7-10)
                             single=[0,1,2,3], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[2], charge=[-1, +1])
# examples for S2sc: N#[N+][S-](O)O
ATOMTYPES['S2d'] = AtomType('S2d', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 8)
                            single=[0], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[2], charge=[0])
# examples for S2d: S=S, C=S, S=O, S=N, S=C=S, S=C=O, S=C=S...
ATOMTYPES['S2dc'] = AtomType('S2dc', generic=['R', 'R!H', 'S', 'Val6'], specific=[],
                             single=[0,1], all_double=[1, 2], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[2], charge=[-1])
# *Composite atomtype; examples for S2dc: [SH-]=[N+]
ATOMTYPES['S2tc'] = AtomType('S2tc', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 10)
                             single=[0], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[1], quadruple=[], benzene=[0], lone_pairs=[2], charge=[-1])
# examples for S2tc: [S-]#[NH+]
ATOMTYPES['S4s'] = AtomType('S4s', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 6-10)
                            single=[0,1,2,3,4], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[0])
# examples for S4s: H4S, SH3CH3...
ATOMTYPES['S4sc'] = AtomType('S4sc', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 5-8)
                             single=[0,1,2,3,4,5], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[-1, +1])
# examples for S4sc: CS[S+]([O-])C, O[SH..-][N+]#N
ATOMTYPES['S4d'] = AtomType('S4d', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 8-10)
                            single=[0,1,2], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[0])
# examples for S4d: O=S(O)O {Sulfurous acid}
ATOMTYPES['S4dd'] = AtomType('S4dd', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 10)
                             single=[0], all_double=[2], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[0])
# examples for S4dd: O=S=O
ATOMTYPES['S4dc'] = AtomType('S4dc', generic=['R', 'R!H', 'S', 'Val6'], specific=[],
                             single=[0,1,2,3,4,5], all_double=[1, 2], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[1], charge=[-1, +1])
# *Composite atomtype; examples for S4dc: [CH2-][S+]=C {where the [CH2-] has a lone pair}, [O+][S-](=O)=O, [O-][S+]=C, [NH-][S+]=C {where the [NH-] has two lone pairs}, [O-][S+]=O
ATOMTYPES['S4b'] = AtomType('S4b', generic=['R', 'R!H', 'S', 'Val6'], specific=[],
                            single=[0], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[], benzene=[2], lone_pairs=[1], charge=[0])
# examples for S4b: Thiophene, Benzothiophene, Benzo[c]thiophene, Thiazole, Benzothiazole...
ATOMTYPES['S4t'] = AtomType('S4t', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 10)
                            single=[0,1], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[1], quadruple=[], benzene=[0], lone_pairs=[1], charge=[0])
# examples for S4t: C#S, C#SO, C#[S]
ATOMTYPES['S4tdc'] = AtomType('S4tdc', generic=['R', 'R!H', 'S', 'Val6'], specific=[],
                              single=[0,1,2], all_double=[0, 1, 2], r_double=[], o_double=[], s_double=[], triple=[1, 2], quadruple=[], benzene=[0], lone_pairs=[1], charge=[-1, +1])
# *Composite atomtype; examples for S4tdc: [C-]#[S+]
ATOMTYPES['S6s'] = AtomType('S6s', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 6-12)
                            single=[0,1,2,3,4,5,6], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[0], benzene=[0], lone_pairs=[0], charge=[0])
# examples for S6s: H6S, F6S
ATOMTYPES['S6sc'] = AtomType('S6sc', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 7-14)
                             single=[0,1,2,3,4,5,6,7], all_double=[0], r_double=[0], o_double=[0], s_double=[0], triple=[0], quadruple=[0], benzene=[0], lone_pairs=[0], charge=[-1, +1, +2])
# examples for S6sc: [O-][S+2](O)(O)[O-]CS(=O)
ATOMTYPES['S6d'] = AtomType('S6d', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 8-12)
                            single=[0,1,2,3,4], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[0])
# examples for S6d: [SH4]=O, SF4=O, [SH4]=C, C[SH3]=C...
ATOMTYPES['S6dd'] = AtomType('S6dd', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 10-12)
                             single=[0,1,2], all_double=[2], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[0])
# examples for S6dd: S(=O)(=O)(O)O {H2SO4, Sulfuric acid}, Perfluorooctanesulfonic acid, Pyrosulfuric acid, Thiosulfuric acid {middle S}, OS(=O)(=O)OOS(=O)(=O)O
ATOMTYPES['S6ddd'] = AtomType('S6ddd', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 12)
                              single=[0], all_double=[3], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[], benzene=[0], lone_pairs=[0], charge=[0])
# examples for S6ddd: O=S(=O)(=O)
ATOMTYPES['S6dc'] = AtomType('S6dc', generic=['R', 'R!H', 'S', 'Val6'], specific=[],
                             single=[0,1,2,3,4,5], all_double=[1, 2, 3], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[0], benzene=[0], lone_pairs=[0], charge=[-1, +1, +2])
# *Composite atomtype; examples for S6dc: O=[S+2]([O-])[O-], [CH-]=[SH3+], [CH-]=[SH2+]O, [CH-][SH2+], O=[S+](=O)[O-], [OH+]=[S-](=O)=O
ATOMTYPES['S6t'] = AtomType('S6t', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 9-12)
                            single=[0,1,2,3], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[1], quadruple=[], benzene=[0], lone_pairs=[0], charge=[0])
# examples for S6t: H3S#N
ATOMTYPES['S6td'] = AtomType('S6td', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 11-12)
                             single=[0,1], all_double=[1], r_double=[], o_double=[], s_double=[], triple=[1], quadruple=[], benzene=[0], lone_pairs=[0], charge=[0])
# examples for S6td: HS(=O)#N
ATOMTYPES['S6tt'] = AtomType('S6tt', generic=['R', 'R!H', 'S', 'Val6'], specific=[],  # (shared electrons = 12)
                             single=[0], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[2], quadruple=[], benzene=[0], lone_pairs=[0], charge=[0])
# examples for S6tt: N#S#N
ATOMTYPES['S6tdc'] = AtomType('S6tdc', generic=['R', 'R!H', 'S', 'Val6'], specific=[],
                              single=[0,1,2,3,4], all_double=[0, 1, 2], r_double=[], o_double=[], s_double=[], triple=[1, 2], quadruple=[], benzene=[0], lone_pairs=[0], charge=[-1, +1])
# *Composite atomtype; examples for S6tdc: [SH2+]#[C-], [N-]=[S+]#N

ATOMTYPES['Cl'] = AtomType('Cl', generic=['R', 'R!H', 'Val7'], specific=['Cl1s'])
ATOMTYPES['Cl1s'] = AtomType('Cl1s', generic=['R', 'R!H', 'Cl', 'Val7'], specific=[],
                             single=[0,1], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[0], benzene=[0], lone_pairs=[3], charge=[0])
# examples for Cl1s: HCl, [Cl]
ATOMTYPES['Br'] = AtomType('Br', generic=['R', 'R!H', 'Val7'], specific=['Br1s'])
ATOMTYPES['Br1s'] = AtomType('Br1s', generic=['R', 'R!H', 'Br', 'Val7'], specific=[],
                             single=[0,1], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[0], benzene=[0], lone_pairs=[3], charge=[0])
# examples for Br1s: HBr, [Br]

ATOMTYPES['I'] = AtomType('I', generic=['R', 'R!H', 'Val7'], specific=['I1s'])
ATOMTYPES['I1s'] = AtomType('I1s', generic=['R', 'R!H', 'I', 'Val7'], specific=[],
                            single=[0,1], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[0], benzene=[0], lone_pairs=[3], charge=[0])
# examples for I1s: HI, [I], IO, CH3I, I2

ATOMTYPES['F'] = AtomType('F', generic=['R', 'R!H', 'Val7'], specific=['F1s'])
ATOMTYPES['F1s'] = AtomType('F1s', generic=['R', 'R!H', 'F', 'Val7'], specific=[],
                            single=[0,1], all_double=[0], r_double=[], o_double=[], s_double=[], triple=[0], quadruple=[0], benzene=[0], lone_pairs=[3], charge=[0])
# examples for F1s: HF, [F], FO, CH3F, F2

ATOMTYPES['X'].set_actions(increment_bond=['X'], decrement_bond=['X'], form_bond=['X'], break_bond=['X'], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Xv'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['Xo'], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Xo'].set_actions(increment_bond=['Xo'], decrement_bond=['Xo'], form_bond=[], break_bond=['Xv'], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])

ATOMTYPES['R'].set_actions(increment_bond=['R'], decrement_bond=['R'], form_bond=['R'], break_bond=['R'], increment_radical=['R'], decrement_radical=['R'], increment_lone_pair=['R'], decrement_lone_pair=['R'])
ATOMTYPES['R!H'].set_actions(increment_bond=['R!H'], decrement_bond=['R!H'], form_bond=['R!H'], break_bond=['R!H'], increment_radical=['R!H'], decrement_radical=['R!H'], increment_lone_pair=['R!H'], decrement_lone_pair=['R!H'])
ATOMTYPES['Val4'].set_actions(increment_bond=['Val4'], decrement_bond=['Val4'], form_bond=['Val4'], break_bond=['Val4'], increment_radical=['Val4'], decrement_radical=['Val4'], increment_lone_pair=['Val4'], decrement_lone_pair=['Val4'])
ATOMTYPES['Val5'].set_actions(increment_bond=['Val5'], decrement_bond=['Val5'], form_bond=['Val5'], break_bond=['Val5'], increment_radical=['Val5'], decrement_radical=['Val5'], increment_lone_pair=['Val5'], decrement_lone_pair=['Val5'])
ATOMTYPES['Val6'].set_actions(increment_bond=['Val6'], decrement_bond=['Val6'], form_bond=['Val6'], break_bond=['Val6'], increment_radical=['Val6'], decrement_radical=['Val6'], increment_lone_pair=['Val6'], decrement_lone_pair=['Val6'])
ATOMTYPES['Val7'].set_actions(increment_bond=['Val7'], decrement_bond=['Val7'], form_bond=['Val7'], break_bond=['Val7'], increment_radical=['Val7'], decrement_radical=['Val7'], increment_lone_pair=['Val7'], decrement_lone_pair=['Val7'])

ATOMTYPES['H'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['H'], break_bond=['H'], increment_radical=['H'], decrement_radical=['H'], increment_lone_pair=[], decrement_lone_pair=[])

ATOMTYPES['He'].set_actions(increment_bond=[], decrement_bond=[], form_bond=[], break_bond=[], increment_radical=['He'], decrement_radical=['He'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Ne'].set_actions(increment_bond=[], decrement_bond=[], form_bond=[], break_bond=[], increment_radical=['Ne'], decrement_radical=['Ne'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Ar'].set_actions(increment_bond=[], decrement_bond=[], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])

ATOMTYPES['C'].set_actions(increment_bond=['C'], decrement_bond=['C'], form_bond=['C'], break_bond=['C'], increment_radical=['C'], decrement_radical=['C'], increment_lone_pair=['C'], decrement_lone_pair=['C'])
ATOMTYPES['Ca'].set_actions(increment_bond=[], decrement_bond=[], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=['C2s'])
ATOMTYPES['Cs'].set_actions(increment_bond=['Cd', 'CO', 'CS'], decrement_bond=[], form_bond=['Cs', 'Csc'], break_bond=['Cs'], increment_radical=['Cs'], decrement_radical=['Cs'], increment_lone_pair=['C2s'], decrement_lone_pair=['C2s'])
ATOMTYPES['Csc'].set_actions(increment_bond=['Cdc'], decrement_bond=[], form_bond=['Csc'], break_bond=['Csc', 'Cs'], increment_radical=['Csc'], decrement_radical=['Csc'], increment_lone_pair=['C2sc'], decrement_lone_pair=['C2sc'])
ATOMTYPES['Cd'].set_actions(increment_bond=['Cdd', 'Ct', 'C2tc'], decrement_bond=['Cs'], form_bond=['Cd', 'Cdc'], break_bond=['Cd'], increment_radical=['Cd'], decrement_radical=['Cd'], increment_lone_pair=['C2d'], decrement_lone_pair=[])
ATOMTYPES['Cdc'].set_actions(increment_bond=[], decrement_bond=['Csc'], form_bond=['Cdc'], break_bond=['Cdc', 'Cd', 'CO', 'CS'], increment_radical=['Cdc'], decrement_radical=['Cdc'], increment_lone_pair=['C2dc'], decrement_lone_pair=[])
ATOMTYPES['CO'].set_actions(increment_bond=['Cdd', 'C2tc'], decrement_bond=['Cs'], form_bond=['CO', 'Cdc'], break_bond=['CO'], increment_radical=['CO'], decrement_radical=['CO'], increment_lone_pair=['C2d'], decrement_lone_pair=[])
ATOMTYPES['CS'].set_actions(increment_bond=['Cdd', 'C2tc'], decrement_bond=['Cs'], form_bond=['CS', 'Cdc'], break_bond=['CS'], increment_radical=['CS'], decrement_radical=['CS'], increment_lone_pair=['C2d'], decrement_lone_pair=[])
ATOMTYPES['Cdd'].set_actions(increment_bond=[], decrement_bond=['Cd', 'CO', 'CS'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Ct'].set_actions(increment_bond=['Cq'], decrement_bond=['Cd', 'CO', 'CS'], form_bond=['Ct'], break_bond=['Ct'], increment_radical=['Ct'], decrement_radical=['Ct'], increment_lone_pair=['C2tc'], decrement_lone_pair=[])
ATOMTYPES['Cb'].set_actions(increment_bond=['Cbf'], decrement_bond=[], form_bond=['Cb'], break_bond=['Cb'], increment_radical=['Cb'], decrement_radical=['Cb'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Cbf'].set_actions(increment_bond=[], decrement_bond=['Cb'], form_bond=[], break_bond=['Cb'], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['C2s'].set_actions(increment_bond=['C2d'], decrement_bond=[], form_bond=['C2s'], break_bond=['C2s'], increment_radical=['C2s'], decrement_radical=['C2s'], increment_lone_pair=['Ca'], decrement_lone_pair=['Cs'])
ATOMTYPES['C2sc'].set_actions(increment_bond=['C2dc'], decrement_bond=[], form_bond=['C2sc'], break_bond=['C2sc'], increment_radical=['C2sc'], decrement_radical=['C2sc'], increment_lone_pair=[], decrement_lone_pair=['Cs'])
ATOMTYPES['C2d'].set_actions(increment_bond=['C2tc'], decrement_bond=['C2s'], form_bond=['C2dc'], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=['Cd', 'CO', 'CS'])
ATOMTYPES['C2dc'].set_actions(increment_bond=[], decrement_bond=['C2sc'], form_bond=[], break_bond=['C2d'], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=['Cdc'])
ATOMTYPES['C2tc'].set_actions(increment_bond=[], decrement_bond=['C2d'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=['Ct'])
ATOMTYPES['Cq'].set_actions(increment_bond=[], decrement_bond=['Ct'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])

ATOMTYPES['N'].set_actions(increment_bond=['N'], decrement_bond=['N'], form_bond=['N'], break_bond=['N'], increment_radical=['N'], decrement_radical=['N'], increment_lone_pair=['N'], decrement_lone_pair=['N'])
ATOMTYPES['N0sc'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['N0sc'], break_bond=['N0sc'], increment_radical=['N0sc'], decrement_radical=['N0sc'], increment_lone_pair=[], decrement_lone_pair=['N1s', 'N1sc'])
ATOMTYPES['N1s'].set_actions(increment_bond=['N1dc'], decrement_bond=[], form_bond=['N1s'], break_bond=['N1s'], increment_radical=['N1s'], decrement_radical=['N1s'], increment_lone_pair=['N0sc'], decrement_lone_pair=['N3s', 'N3sc'])
ATOMTYPES['N1sc'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['N1sc'], break_bond=['N1sc'], increment_radical=['N1sc'], decrement_radical=['N1sc'], increment_lone_pair=[], decrement_lone_pair=['N3s', 'N3sc'])
ATOMTYPES['N1dc'].set_actions(increment_bond=['N1dc'], decrement_bond=['N1s', 'N1dc'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=['N3d'])
ATOMTYPES['N3s'].set_actions(increment_bond=['N3d'], decrement_bond=[], form_bond=['N3s'], break_bond=['N3s'], increment_radical=['N3s'], decrement_radical=['N3s'], increment_lone_pair=['N1s', 'N1sc'], decrement_lone_pair=['N5sc'])
ATOMTYPES['N3sc'].set_actions(increment_bond=['N3d'], decrement_bond=[], form_bond=['N3sc'], break_bond=['N3sc'], increment_radical=['N3sc'], decrement_radical=['N3sc'], increment_lone_pair=['N1s', 'N1sc'], decrement_lone_pair=['N5sc'])
ATOMTYPES['N3d'].set_actions(increment_bond=['N3t'], decrement_bond=['N3s', 'N3sc'], form_bond=['N3d'], break_bond=['N3d'], increment_radical=['N3d'], decrement_radical=['N3d'], increment_lone_pair=['N1dc'], decrement_lone_pair=['N5dc'])
ATOMTYPES['N3t'].set_actions(increment_bond=[], decrement_bond=['N3d'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['N3b'].set_actions(increment_bond=[], decrement_bond=[], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['N5sc'].set_actions(increment_bond=['N5dc'], decrement_bond=[], form_bond=['N5sc'], break_bond=['N5sc'], increment_radical=['N5sc'], decrement_radical=['N5sc'], increment_lone_pair=['N3s', 'N3sc'], decrement_lone_pair=[])
ATOMTYPES['N5dc'].set_actions(increment_bond=['N5ddc', 'N5tc'], decrement_bond=['N5sc'], form_bond=['N5dc'], break_bond=['N5dc'], increment_radical=['N5dc'], decrement_radical=['N5dc'], increment_lone_pair=['N3d'], decrement_lone_pair=[])
ATOMTYPES['N5ddc'].set_actions(increment_bond=['N5dddc'], decrement_bond=['N5dc'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['N5dddc'].set_actions(increment_bond=[], decrement_bond=['N5ddc'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['N5tc'].set_actions(increment_bond=[], decrement_bond=['N5dc'], form_bond=['N5tc'], break_bond=['N5tc'], increment_radical=['N5tc'], decrement_radical=['N5tc'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['N5b'].set_actions(increment_bond=['N5bd'], decrement_bond=[], form_bond=['N5b'], break_bond=['N5b'], increment_radical=['N5b'], decrement_radical=['N5b'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['N5bd'].set_actions(increment_bond=[], decrement_bond=['N5b'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])

ATOMTYPES['O'].set_actions(increment_bond=['O'], decrement_bond=['O'], form_bond=['O'], break_bond=['O'], increment_radical=['O'], decrement_radical=['O'], increment_lone_pair=['O'], decrement_lone_pair=['O'])
ATOMTYPES['Oa'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['O0sc'], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=['O2s', 'O2sc'])
ATOMTYPES['O0sc'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['O0sc'], break_bond=['Oa', 'O0sc'], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=['O2s', 'O2sc'])
ATOMTYPES['O2s'].set_actions(increment_bond=['O2d'], decrement_bond=[], form_bond=['O2s', 'O2sc'], break_bond=['O2s'], increment_radical=['O2s'], decrement_radical=['O2s'], increment_lone_pair=['Oa', 'O0sc'], decrement_lone_pair=['O4sc'])
ATOMTYPES['O2sc'].set_actions(increment_bond=['O2d'], decrement_bond=[], form_bond=[], break_bond=['O2s'], increment_radical=['O2sc'], decrement_radical=['O2sc'], increment_lone_pair=[], decrement_lone_pair=['O4sc'])
ATOMTYPES['O2d'].set_actions(increment_bond=[], decrement_bond=['O2s', 'O2sc'], form_bond=[], break_bond=[], increment_radical=['O2d'], decrement_radical=['O2d'], increment_lone_pair=[], decrement_lone_pair=['O4dc', 'O4tc'])
ATOMTYPES['O4sc'].set_actions(increment_bond=['O4dc'], decrement_bond=[], form_bond=['O4sc'], break_bond=['O4sc'], increment_radical=['O4sc'], decrement_radical=['O4sc'], increment_lone_pair=['O2s', 'O2sc'], decrement_lone_pair=[])
ATOMTYPES['O4dc'].set_actions(increment_bond=['O4tc'], decrement_bond=['O4sc'], form_bond=['O4dc'], break_bond=['O4dc'], increment_radical=['O4dc'], decrement_radical=['O4dc'], increment_lone_pair=['O2d'], decrement_lone_pair=[])
ATOMTYPES['O4tc'].set_actions(increment_bond=[], decrement_bond=['O4dc'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['O4b'].set_actions(increment_bond=[], decrement_bond=[], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])

ATOMTYPES['Ne'].set_actions(increment_bond=[], decrement_bond=[], form_bond=[], break_bond=[], increment_radical=['Ne'], decrement_radical=['Ne'], increment_lone_pair=[], decrement_lone_pair=[])

ATOMTYPES['Si'].set_actions(increment_bond=['Si'], decrement_bond=['Si'], form_bond=['Si'], break_bond=['Si'], increment_radical=['Si'], decrement_radical=['Si'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Sis'].set_actions(increment_bond=['Sid', 'SiO'], decrement_bond=[], form_bond=['Sis'], break_bond=['Sis'], increment_radical=['Sis'], decrement_radical=['Sis'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Sid'].set_actions(increment_bond=['Sidd', 'Sit'], decrement_bond=['Sis'], form_bond=['Sid'], break_bond=['Sid'], increment_radical=['Sid'], decrement_radical=['Sid'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Sidd'].set_actions(increment_bond=[], decrement_bond=['Sid', 'SiO'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Sit'].set_actions(increment_bond=['Siq'], decrement_bond=['Sid'], form_bond=['Sit'], break_bond=['Sit'], increment_radical=['Sit'], decrement_radical=['Sit'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['SiO'].set_actions(increment_bond=['Sidd'], decrement_bond=['Sis'], form_bond=['SiO'], break_bond=['SiO'], increment_radical=['SiO'], decrement_radical=['SiO'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Sib'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['Sib'], break_bond=['Sib'], increment_radical=['Sib'], decrement_radical=['Sib'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Sibf'].set_actions(increment_bond=[], decrement_bond=[], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Siq'].set_actions(increment_bond=[], decrement_bond=['Sit'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])

ATOMTYPES['P'].set_actions(increment_bond=['P'], decrement_bond=['P'], form_bond=['P'], break_bond=['P'], increment_radical=['P'], decrement_radical=['P'], increment_lone_pair=['P'], decrement_lone_pair=['P'])
ATOMTYPES['P0sc'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['P0sc'], break_bond=['P0sc'], increment_radical=['P0sc'], decrement_radical=['P0sc'], increment_lone_pair=[], decrement_lone_pair=['P1s', 'P1sc'])
ATOMTYPES['P1s'].set_actions(increment_bond=['P1dc'], decrement_bond=[], form_bond=['P1s'], break_bond=['P1s'], increment_radical=['P1s'], decrement_radical=['P1s'], increment_lone_pair=['P0sc'], decrement_lone_pair=['P3s'])
ATOMTYPES['P1sc'].set_actions(increment_bond=['P1dc'], decrement_bond=[], form_bond=['P1sc'], break_bond=['P1sc'], increment_radical=['P1sc'], decrement_radical=['P1sc'], increment_lone_pair=['P0sc'], decrement_lone_pair=['P3s'])
ATOMTYPES['P1dc'].set_actions(increment_bond=[], decrement_bond=['P1s'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=['P3d'])
ATOMTYPES['P3s'].set_actions(increment_bond=['P3d'], decrement_bond=[], form_bond=['P3s'], break_bond=['P3s'], increment_radical=['P3s'], decrement_radical=['P3s'], increment_lone_pair=['P1s', 'P1sc'], decrement_lone_pair=['P5s', 'P5sc'])
ATOMTYPES['P3d'].set_actions(increment_bond=['P3t'], decrement_bond=['P3s'], form_bond=['P3d'], break_bond=['P3d'], increment_radical=['P3d'], decrement_radical=['P3d'], increment_lone_pair=['P1dc'], decrement_lone_pair=['P5d', 'P5dc'])
ATOMTYPES['P3t'].set_actions(increment_bond=[], decrement_bond=['P3d'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=['P5t', 'P5tc'])
ATOMTYPES['P3b'].set_actions(increment_bond=[], decrement_bond=[], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['P5s'].set_actions(increment_bond=['P5d', 'P5dc'], decrement_bond=[], form_bond=['P5s'], break_bond=['P5s'], increment_radical=['P5s'], decrement_radical=['P5s'], increment_lone_pair=['P3s'], decrement_lone_pair=[])
ATOMTYPES['P5sc'].set_actions(increment_bond=['P5dc'], decrement_bond=[], form_bond=['P5sc'], break_bond=['P5sc'], increment_radical=['P5sc'], decrement_radical=['P5sc'], increment_lone_pair=['P3s'], decrement_lone_pair=[])
ATOMTYPES['P5d'].set_actions(increment_bond=['P5dd', 'P5ddc', 'P5t', 'P5tc'], decrement_bond=['P5s'], form_bond=['P5d'], break_bond=['P5d'], increment_radical=['P5d'], decrement_radical=['P5d'], increment_lone_pair=['P3d'], decrement_lone_pair=[])
ATOMTYPES['P5dd'].set_actions(increment_bond=['P5td'], decrement_bond=['P5d', 'P5dc'], form_bond=['P5dd'], break_bond=['P5dd'], increment_radical=['P5dd'], decrement_radical=['P5dd'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['P5dc'].set_actions(increment_bond=['P5dd', 'P5ddc', 'P5tc'], decrement_bond=['P5sc'], form_bond=['P5dc'], break_bond=['P5dc'], increment_radical=['P5dc'], decrement_radical=['P5dc'], increment_lone_pair=['P3d'], decrement_lone_pair=[])
ATOMTYPES['P5ddc'].set_actions(increment_bond=[], decrement_bond=['P5dc'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['P5t'].set_actions(increment_bond=['P5td'], decrement_bond=['P5d'], form_bond=['P5t'], break_bond=['P5t'], increment_radical=['P5t'], decrement_radical=['P5t'], increment_lone_pair=['P3t'], decrement_lone_pair=[])
ATOMTYPES['P5td'].set_actions(increment_bond=[], decrement_bond=['P5t', 'P5dd'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['P5tc'].set_actions(increment_bond=[], decrement_bond=['P5dc'], form_bond=['P5tc'], break_bond=['P5tc'], increment_radical=['P5tc'], decrement_radical=['P5tc'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['P5b'].set_actions(increment_bond=['P5bd'], decrement_bond=[], form_bond=['P5b'], break_bond=['P5b'], increment_radical=['P5b'], decrement_radical=['P5b'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['P5bd'].set_actions(increment_bond=[], decrement_bond=['P5b'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])

ATOMTYPES['S'].set_actions(increment_bond=['S'], decrement_bond=['S'], form_bond=['S'], break_bond=['S'], increment_radical=['S'], decrement_radical=['S'], increment_lone_pair=['S'], decrement_lone_pair=['S'])
ATOMTYPES['S0sc'].set_actions(increment_bond=['S0sc'], decrement_bond=['S0sc'], form_bond=['S0sc'], break_bond=['Sa', 'S0sc'], increment_radical=['S0sc'], decrement_radical=['S0sc'], increment_lone_pair=[], decrement_lone_pair=['S2s', 'S2sc', 'S2dc', 'S2tc'])
ATOMTYPES['Sa'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['S0sc'], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=['S2s'])
ATOMTYPES['S2s'].set_actions(increment_bond=['S2d', 'S2dc'], decrement_bond=[], form_bond=['S2s', 'S2sc'], break_bond=['S2s'], increment_radical=['S2s'], decrement_radical=['S2s'], increment_lone_pair=['Sa', 'S0sc'], decrement_lone_pair=['S4s', 'S4sc'])
ATOMTYPES['S2sc'].set_actions(increment_bond=['S2dc'], decrement_bond=[], form_bond=['S2sc'], break_bond=['S2sc', 'S2s'], increment_radical=['S2sc'], decrement_radical=['S2sc'], increment_lone_pair=['S0sc'], decrement_lone_pair=['S4s', 'S4sc'])
ATOMTYPES['S2d'].set_actions(increment_bond=['S2tc'], decrement_bond=['S2s'], form_bond=['S2d'], break_bond=['S2d'], increment_radical=['S2d'], decrement_radical=['S2d'], increment_lone_pair=[], decrement_lone_pair=['S4dc', 'S4d'])
ATOMTYPES['S2dc'].set_actions(increment_bond=['S2tc', 'S2dc'], decrement_bond=['S2sc', 'S2s', 'S2dc'], form_bond=['S2dc'], break_bond=['S2dc'], increment_radical=['S2dc'], decrement_radical=['S2dc'], increment_lone_pair=['S0sc'], decrement_lone_pair=['S4d', 'S4dc'])
ATOMTYPES['S2tc'].set_actions(increment_bond=[], decrement_bond=['S2d', 'S2dc'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=['S0sc'], decrement_lone_pair=['S4t'])
ATOMTYPES['S4s'].set_actions(increment_bond=['S4d', 'S4dc'], decrement_bond=[], form_bond=['S4s'], break_bond=['S4s'], increment_radical=['S4s'], decrement_radical=['S4s'], increment_lone_pair=['S2s', 'S2sc'], decrement_lone_pair=['S6s'])
ATOMTYPES['S4sc'].set_actions(increment_bond=['S4d', 'S4dc'], decrement_bond=[], form_bond=['S4s', 'S4sc'], break_bond=['S4sc'], increment_radical=['S4sc'], decrement_radical=['S4sc'], increment_lone_pair=['S2s', 'S2sc'], decrement_lone_pair=['S6s'])
ATOMTYPES['S4d'].set_actions(increment_bond=['S4dd', 'S4dc', 'S4t', 'S4tdc'], decrement_bond=['S4s', 'S4sc'], form_bond=['S4dc', 'S4d'], break_bond=['S4d', 'S4dc'], increment_radical=['S4d'], decrement_radical=['S4d'], increment_lone_pair=['S2d', 'S2dc'], decrement_lone_pair=['S6d', 'S6dc'])
ATOMTYPES['S4dc'].set_actions(increment_bond=['S4dd', 'S4dc', 'S4tdc'], decrement_bond=['S4sc', 'S4dc'], form_bond=['S4d', 'S4dc'], break_bond=['S4d', 'S4dc'], increment_radical=['S4dc'], decrement_radical=['S4dc'], increment_lone_pair=['S2d', 'S2dc'], decrement_lone_pair=['S6d', 'S6dc'])
ATOMTYPES['S4b'].set_actions(increment_bond=[], decrement_bond=[], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['S4dd'].set_actions(increment_bond=['S4dc'], decrement_bond=['S4dc', 'S4d'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=['S6dd'])
ATOMTYPES['S4t'].set_actions(increment_bond=[], decrement_bond=['S4d'], form_bond=['S4t'], break_bond=['S4t'], increment_radical=['S4t'], decrement_radical=['S4t'], increment_lone_pair=['S2tc'], decrement_lone_pair=['S6t', 'S6tdc'])
ATOMTYPES['S4tdc'].set_actions(increment_bond=['S4tdc'], decrement_bond=['S4d', 'S4tdc'], form_bond=['S4tdc'], break_bond=['S4tdc'], increment_radical=['S4tdc'], decrement_radical=['S4tdc'], increment_lone_pair=['S6tdc'], decrement_lone_pair=['S6td', 'S6tdc'])
ATOMTYPES['S6s'].set_actions(increment_bond=['S6d', 'S6dc'], decrement_bond=[], form_bond=['S6s'], break_bond=['S6s'], increment_radical=['S6s'], decrement_radical=['S6s'], increment_lone_pair=['S4s', 'S4sc'], decrement_lone_pair=[])
ATOMTYPES['S6sc'].set_actions(increment_bond=['S6dc'], decrement_bond=[], form_bond=['S6sc'], break_bond=['S6sc'], increment_radical=['S6sc'], decrement_radical=['S6sc'], increment_lone_pair=['S4s', 'S4sc'], decrement_lone_pair=[])
ATOMTYPES['S6d'].set_actions(increment_bond=['S6dd', 'S6t', 'S6tdc'], decrement_bond=['S6s'], form_bond=['S6d', 'S6dc'], break_bond=['S6d', 'S6dc'], increment_radical=['S6d'], decrement_radical=['S6d'], increment_lone_pair=['S4d', 'S4dc'], decrement_lone_pair=[])
ATOMTYPES['S6dc'].set_actions(increment_bond=['S6dd', 'S6ddd', 'S6dc', 'S6t', 'S6td', 'S6tdc'], decrement_bond=['S6sc', 'S6dc'], form_bond=['S6d', 'S6dc'], break_bond=['S6d', 'S6dc'], increment_radical=['S6dc'], decrement_radical=['S6dc'], increment_lone_pair=['S4d', 'S4dc'], decrement_lone_pair=[])
ATOMTYPES['S6dd'].set_actions(increment_bond=['S6ddd', 'S6td'], decrement_bond=['S6d', 'S6dc'], form_bond=['S6dd', 'S6dc'], break_bond=['S6dd'], increment_radical=['S6dd'], decrement_radical=['S6dd'], increment_lone_pair=['S4dd'], decrement_lone_pair=[])
ATOMTYPES['S6ddd'].set_actions(increment_bond=[], decrement_bond=['S6dd', 'S6dc'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['S6t'].set_actions(increment_bond=['S6td'], decrement_bond=['S6d', 'S6dc'], form_bond=['S6t'], break_bond=['S6t'], increment_radical=['S6t'], decrement_radical=['S6t'], increment_lone_pair=['S4t'], decrement_lone_pair=[])
ATOMTYPES['S6td'].set_actions(increment_bond=['S6tt', 'S6tdc'], decrement_bond=['S6dc', 'S6t', 'S6dd', 'S6tdc'], form_bond=['S6td'], break_bond=['S6td'], increment_radical=['S6td'], decrement_radical=['S6td'], increment_lone_pair=['S4tdc'], decrement_lone_pair=[])
ATOMTYPES['S6tt'].set_actions(increment_bond=[], decrement_bond=['S6td', 'S6tdc'], form_bond=[], break_bond=[], increment_radical=[], decrement_radical=[], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['S6tdc'].set_actions(increment_bond=['S6td', 'S6tdc', 'S6tt'], decrement_bond=['S6dc', 'S6tdc'], form_bond=['S6tdc'], break_bond=['S6tdc'], increment_radical=['S6tdc'], decrement_radical=['S6tdc'], increment_lone_pair=['S4t', 'S4tdc'], decrement_lone_pair=[])

ATOMTYPES['Cl'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['Cl'], break_bond=['Cl'], increment_radical=['Cl'], decrement_radical=['Cl'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Cl1s'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['Cl1s'], break_bond=['Cl1s'], increment_radical=['Cl1s'], decrement_radical=['Cl1s'], increment_lone_pair=[], decrement_lone_pair=[])

ATOMTYPES['Br'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['Br'], break_bond=['Br'], increment_radical=['Br'], decrement_radical=['Br'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['Br1s'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['Br1s'], break_bond=['Br1s'], increment_radical=['Br1s'], decrement_radical=['Br1s'], increment_lone_pair=[], decrement_lone_pair=[])

ATOMTYPES['I'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['I'], break_bond=['I'], increment_radical=['I'], decrement_radical=['I'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['I1s'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['I1s'], break_bond=['I1s'], increment_radical=['I1s'], decrement_radical=['I1s'], increment_lone_pair=[], decrement_lone_pair=[])

ATOMTYPES['F'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['F'], break_bond=['F'], increment_radical=['F'], decrement_radical=['F'], increment_lone_pair=[], decrement_lone_pair=[])
ATOMTYPES['F1s'].set_actions(increment_bond=[], decrement_bond=[], form_bond=['F1s'], break_bond=['F1s'], increment_radical=['F1s'], decrement_radical=['F1s'], increment_lone_pair=[], decrement_lone_pair=[])

# these are ordered in priority of picking if a more general atomtype is encountered
allElements = ['H', 'C', 'O', 'N', 'S', 'P', 'Si', 'F', 'Cl', 'Br', 'I', 'Ne', 'Ar', 'He', 'X']
# list of elements that do not have more specific atomTypes
nonSpecifics = ['H', 'He', 'Ne', 'Ar',]

for atomtype in ATOMTYPES.values():
    for items in [atomtype.generic, atomtype.specific,
                  atomtype.increment_bond, atomtype.decrement_bond, atomtype.form_bond,
                  atomtype.break_bond, atomtype.increment_radical, atomtype.decrement_radical, atomtype.increment_lone_pair,
                  atomtype.decrement_lone_pair]:
        for index in range(len(items)):
            items[index] = ATOMTYPES[items[index]]


def get_features(atom, bonds):
    """
    Returns a list of features needed to determine atomtype for :class:'Atom'
    or :class:'GroupAtom' object 'atom and with local bond structure `bonds`,
    a ``dict`` containing atom-bond pairs.
    """
    cython.declare(single=cython.int, all_double=cython.int, r_double=cython.int,
                   s_double=cython.int, o_double=cython.int, triple=cython.int,
                   benzene=cython.int, quadruple=cython.int)
    cython.declare(features=cython.list)

    # Count numbers of each higher-order bond type
    single = r_double = o_double = s_double = triple = benzene = quadruple = 0
    for atom2, bond12 in bonds.items():
        if bond12.is_single():
            single += 1
        elif bond12.is_double():
            if atom2.is_oxygen():
                o_double += 1
            elif atom2.is_sulfur():
                s_double += 1
            else:
                # r_double is for double bonds NOT to oxygen or Sulfur
                r_double += 1
        elif bond12.is_triple():
            triple += 1
        elif bond12.is_benzene():
            benzene += 1
        elif bond12.is_quadruple():
            quadruple += 1

    # all_double is for all double bonds, to anything
    all_double = r_double + o_double + s_double
    # Warning: some parts of code assume this list matches the list returned by count_bonds()
    # possibly the two methods could be merged or one could call the other.
    features = [single, all_double, r_double, o_double, s_double, triple, quadruple, benzene, atom.lone_pairs, atom.charge]

    return features


def get_atomtype(atom, bonds):
    """
    Determine the appropriate atom type for an :class:`Atom` object `atom`
    with local bond structure `bonds`, a ``dict`` containing atom-bond pairs.
    """

    cython.declare(atom_symbol=str)
    cython.declare(mol_feature_list=cython.list, atomtype_feature_list=cython.list)

    # Use element and counts to determine proper atom type
    atom_symbol = atom.symbol
    # These elements do not do not have a more specific atomtype
    if atom_symbol in nonSpecifics:
        return ATOMTYPES[atom_symbol]

    mol_feature_list = get_features(atom, bonds)
    for specific_atom_type in ATOMTYPES[atom_symbol].specific:
        atomtype_feature_list = specific_atom_type.get_features()
        for mol_feature, atomtype_feature in zip(mol_feature_list, atomtype_feature_list):
            if atomtype_feature == []:
                continue
            elif mol_feature not in atomtype_feature:
                break
        else:
            return specific_atom_type
    else:
        single, all_double, r_double, o_double, s_double, triple, quadruple, benzene, lone_pairs, charge = mol_feature_list

        raise AtomTypeError(
            f'Unable to determine atom type for atom {atom}, which has {single:d} single bonds, '
            f'{all_double:d} double bonds ({o_double:d} to O, {s_double:d} to S, '
            f'{r_double:d} others), {triple:d} triple bonds, {quadruple:d} quadruple bonds, '
            f'{benzene:d} benzene bonds, {lone_pairs:d} lone pairs, and {charge:+d} charge.')
