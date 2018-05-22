#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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
themselves are available in the ``atomTypes`` module-level variable, or as
the return value from the :meth:`getAtomType()` method.

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

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `label`             ``str``             A unique label for the atom type
    `generic`           ``list``            The atom types that are more generic than this one
    `specific`          ``list``            The atom types that are more specific than this one
    `incrementBond`     ``list``            The atom type(s) that result when an adjacent bond's order is incremented
    `decrementBond`     ``list``            The atom type(s) that result when an adjacent bond's order is decremented
    `formBond`          ``list``            The atom type(s) that result when a new single bond is formed to this atom type
    `breakBond`         ``list``            The atom type(s) that result when an existing single bond to this atom type is broken
    `incrementRadical`  ``list``            The atom type(s) that result when the number of radical electrons is incremented
    `decrementRadical`  ``list``            The atom type(s) that result when the number of radical electrons is decremented
    `incrementLonePair` ``list``            The atom type(s) that result when the number of lone electron pairs is incremented
    `decrementLonePair` ``list``            The atom type(s) that result when the number of lone electron pairs is decremented

    The following features are what are required in a given atomtype. Any int in the list is acceptable. An empty list is a wildcard
    ----------------------------------------------------------------------------
    'single'            ''list''            The total number of single bonds on the atom
    'allDouble'         ''list''            The total number of double bonds on the atom
    'rDouble'           ''list''            The number of double bonds to any non-oxygen, nonsulfur
    'oDouble'           ''list''            The number of double bonds to oxygen
    'sDouble'           ''list''            The number of double bonds to sulfur
    'triple'            ''list''            The total number of triple bonds on the atom
    'benzene'           ''list''            The total number of benzene bonds on the atom
    'lonePairs'         ''list''            The number of lone pairs on the atom
    'charge'            ''list''            The partial charge of the atom
    =================== =================== ====================================

    """
    

    def __init__(self, label='', generic=None, specific=None,
                 single=None,
                 allDouble=None,
                 rDouble=None,
                 oDouble=None,
                 sDouble=None,
                 triple=None,
                 benzene=None,
                 lonePairs=None,
                 charge=None):
        self.label = label
        self.generic = generic or []
        self.specific = specific or []
        self.incrementBond = []
        self.decrementBond = []
        self.formBond = []
        self.breakBond = []
        self.incrementRadical = []
        self.decrementRadical = []
        self.incrementLonePair = []
        self.decrementLonePair = []
        self.single = single or []
        self.allDouble = allDouble or []
        self.rDouble = rDouble or []
        self.oDouble = oDouble or []
        self.sDouble = sDouble or []
        self.triple = triple or []
        self.benzene = benzene or []
        self.lonePairs = lonePairs or []
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
            'incrementBond': self.incrementBond,
            'decrementBond': self.decrementBond,
            'formBond': self.formBond,
            'breakBond': self.breakBond,
            'incrementRadical': self.incrementRadical,
            'decrementRadical': self.decrementRadical,
            'incrementLonePair': self.incrementLonePair,
            'decrementLonePair': self.decrementLonePair,
            'single': self.single,
            'allDouble': self.allDouble,
            'rDouble': self.rDouble,
            'oDouble': self.oDouble,
            'sDouble': self.sDouble,
            'triple': self.triple,
            'benzene': self.benzene,
            'lonePairs': self.lonePairs,
            'charge': self.charge
        }
        return (AtomType, (), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling an AtomType object.
        """
        self.label = d['label']
        self.generic = d['generic']
        self.specific = d['specific']
        self.incrementBond = d['incrementBond']
        self.decrementBond = d['decrementBond']
        self.formBond = d['formBond']
        self.breakBond = d['breakBond']
        self.incrementRadical = d['incrementRadical']
        self.decrementRadical = d['decrementRadical']
        self.incrementLonePair = d['incrementLonePair']
        self.decrementLonePair = d['decrementLonePair']
        self.single = d['single']
        self.allDouble = d['allDouble']
        self.rDouble = d['rDouble']
        self.oDouble = d['oDouble']
        self.sDouble = d['sDouble']
        self.triple = d['triple']
        self.benzene = d['benzene']
        self.lonePairs = d['lonePairs']
        self.charge = d['charge']

    def setActions(self, incrementBond, decrementBond, formBond, breakBond, incrementRadical, decrementRadical, incrementLonePair, decrementLonePair):
        self.incrementBond = incrementBond
        self.decrementBond = decrementBond
        self.formBond = formBond
        self.breakBond = breakBond
        self.incrementRadical = incrementRadical
        self.decrementRadical = decrementRadical
        self.incrementLonePair = incrementLonePair
        self.decrementLonePair = decrementLonePair

    def equivalent(self, other):
        """
        Returns ``True`` if two atom types `atomType1` and `atomType2` are
        equivalent or ``False``  otherwise. This function respects wildcards,
        e.g. ``R!H`` is equivalent to ``C``.
        """
        return self is other or self in other.specific or other in self.specific

    def isSpecificCaseOf(self, other):
        """
        Returns ``True`` if atom type `atomType1` is a specific case of
        atom type `atomType2` or ``False``  otherwise.
        """
        return self is other or self in other.specific
    
    def getFeatures(self):
        """
        Returns a list of the features that are checked to determine atomtype
        """
        features=[self.single,
                  self.allDouble,
                  self.rDouble,
                  self.oDouble,
                  self.sDouble,
                  self.triple,
                  self.benzene,
                  self.lonePairs,
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
- N3d is nitrogen with valence=3 (i.e., 3 electronce are able to form bonds or remain as radicals) with one double bond
- S2tc is a charged sulful with valence=2 with a triple bonds
- Oa is atomic oxygen, i.e., a closed shell atom
Some charged atom types were merged together, and are marked as '*Composite atomType'
"""

atomTypes = {}
atomTypes['R']    = AtomType(label='R', generic=[], specific=[
    'H',
    'R!H',
    'Val4','Val5','Val6','Val7',    
    'He','Ne','Ar',
    'C','Ca','Cs','Csc','Cd','CO','CS','Cdd','Cdc','Ct','Cb','Cbf','C2s','C2sc','C2d','C2dc','C2tc',
    'N','N0sc','N1s','N1sc','N1dc','N3s','N3sc','N3d','N3t','N3b','N5sc','N5dc','N5ddc','N5dddc','N5t','N5tc','N5b','N5bd',
    'O','Oa','O0sc','O2s','O2sc','O2d','O4sc','O4dc','O4tc','O4b',
    'Si','Sis','Sid','Sidd','Sit','SiO','Sib','Sibf',
    'S','Sa','S0sc','S2s','S2sc','S2d','S2dc','S2tc','S4s','S4sc','S4d','S4dd','S4dc','S4b','S4t','S4tdc','S6s','S6sc','S6d','S6dd','S6ddd','S6dc','S6t','S6td','S6tt','S6tdc',
    'Cl','Cl1s',
    'I','I1s'])

atomTypes['R!H']  = AtomType(label='R!H', generic=['R'], specific=[
    'Val4','Val5','Val6','Val7',
    'He','Ne','Ar',
    'C','Ca','Cs','Csc','Cd','CO','CS','Cdd','Cdc','Ct','Cb','Cbf','C2s','C2sc','C2d','C2dc','C2tc',
    'N','N0sc','N1s','N1sc','N1dc','N3s','N3sc','N3d','N3t','N3b','N5sc','N5dc','N5ddc','N5dddc','N5t','N5tc','N5b','N5bd',
    'O','Oa','O0sc','O2s','O2sc','O2d','O4sc','O4dc','O4tc','O4b',
    'Si','Sis','Sid','Sidd','Sit','SiO','Sib','Sibf',
    'S','Sa','S0sc','S2s','S2sc','S2d','S2dc','S2tc','S4s','S4sc','S4d','S4dd','S4dc','S4b','S4t','S4tdc','S6s','S6sc','S6d','S6dd','S6ddd','S6dc','S6t','S6td','S6tt','S6tdc',
    'Cl','Cl1s',
    'I','I1s'])

atomTypes['Val4'] = AtomType(label='Val4', generic=['R','R!H'], specific=[
    'C','Ca','Cs','Csc','Cd','CO','CS','Cdd','Cdc','Ct','Cb','Cbf','C2s','C2sc','C2d','C2dc','C2tc',
    'Si','Sis','Sid','Sidd','Sit','SiO','Sib','Sibf'])

atomTypes['Val5'] = AtomType(label='Val5', generic=['R','R!H'], specific=[
    'N','N0sc','N1s','N1sc','N1dc','N3s','N3sc','N3d','N3t','N3b','N5sc','N5dc','N5ddc','N5dddc','N5t','N5tc','N5b','N5bd'])

atomTypes['Val6'] = AtomType(label='Val6', generic=['R','R!H'], specific=[
    'O','Oa','O0sc','O2s','O2sc','O2d','O4sc','O4dc','O4tc','O4b',
    'S','Sa','S0sc','S2s','S2sc','S2d','S2dc','S2tc','S4s','S4sc','S4d','S4dd','S4dc','S4b','S4t','S4tdc','S6s','S6sc','S6d','S6dd','S6ddd','S6dc','S6t','S6td','S6tt','S6tdc'])

atomTypes['Val7'] = AtomType(label='Val7', generic=['R','R!H'], specific=[
    'Cl','Cl1s',
    'I','I1s'])

atomTypes['H'   ] = AtomType('H',    generic=['R'],            specific=[])

atomTypes['He'  ] = AtomType('He',   generic=['R','R!H'],      specific=[])
atomTypes['Ne'  ] = AtomType('Ne',   generic=['R','R!H'],      specific=[])
atomTypes['Ar'  ] = AtomType('Ar',   generic=['R','R!H'],      specific=[])

atomTypes['C'   ] = AtomType('C',    generic=['R','R!H','Val4'],      specific=['Ca','Cs','Csc','Cd','CO','CS','Cdd','Cdc','Ct','Cb','Cbf','C2s','C2sc','C2d','C2dc','C2tc'],
                             single=[], allDouble=[], rDouble=[], oDouble=[], sDouble=[], triple=[], benzene=[], lonePairs=[], charge=[])
atomTypes['Ca'  ] = AtomType('Ca',   generic=['R','R!H','C','Val4'],  specific=[],  # (shared electrons = 4)
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[0])
# examples for Ca: atomic carbon (closed shell)
atomTypes['Cs'  ] = AtomType('Cs',   generic=['R','R!H','C','Val4'],  specific=[],  # (shared electrons = 4-8)
                             single=[0,1,2,3,4], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for Cs: C, CC,
atomTypes['Csc' ] = AtomType('Csc',  generic=['R','R!H','C','Val4'],  specific=[],  # (shared electrons = 3-6)
                             single=[0,1,2,3], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[+1])
# examples for Csc: C1=CCC([O-])[CH+]1, O[O+]=C[C+]C([O-])[O-]
atomTypes['Cd'  ] = AtomType('Cd',   generic=['R','R!H','C','Val4'],  specific=[],  # (shared electrons = 6-8)
                             single=[0,1,2], allDouble=[1], rDouble=[1], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for Cd: C=C, C=N
atomTypes['Cdc' ] = AtomType('Cdc',  generic=['R','R!H','C','Val4'],  specific=[],  # (shared electrons = 6)
                             single=[0,1], allDouble=[1], rDouble=[0,1], oDouble=[0,1], sDouble=[0,1], triple=[0], benzene=[0], lonePairs=[0], charge=[+1])
# examples for Cdc: [CH+]=C=[CH-], [CH+]=N[O-] (one of the res structures of Fulminic acid)
atomTypes['CO'  ] = AtomType('CO',   generic=['R','R!H','C','Val4'],  specific=[],  # (shared electrons = 6-8)
                             single=[0,1,2], allDouble=[1], rDouble=[0], oDouble=[1], sDouble=[0], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for CO: C=O
atomTypes['CS'  ] = AtomType('CS',   generic=['R','R!H','C','Val4'],  specific=[],  # (shared electrons = 6-8)
                             single=[0,1,2], allDouble=[1], rDouble=[0], oDouble=[0], sDouble=[1], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for CS: C=S
atomTypes['Cdd' ] = AtomType('Cdd',  generic=['R','R!H','C','Val4'],  specific=[],  # (shared electrons = 8)
                             single=[0], allDouble=[2], rDouble=[0,1,2], oDouble=[0,1,2], sDouble=[0,1,2], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for Cdd: O=C=O, C=C=C
atomTypes['Ct'  ] = AtomType('Ct',   generic=['R','R!H','C','Val4'],  specific=[],  # (shared electrons = 7-8)
                             single=[0,1], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[0], charge=[0])
# examples for Ct: C#C, C#N
atomTypes['Cb'  ] = AtomType('Cb',   generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[1,2], lonePairs=[], charge=[])
# examples for Cb: benzene (C6H6)
atomTypes['Cbf' ] = AtomType('Cbf',  generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[3], lonePairs=[], charge=[])
# examples for Cbf: Naphthalene
atomTypes['C2s' ] = AtomType('C2s',  generic=['R','R!H','C','Val4'],  specific=[],  # (shared electrons = 4-6)
                             single=[0,1,2], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[0])
# examples for C2s: singlet[CH2]
atomTypes['C2sc'] = AtomType('C2sc', generic=['R','R!H','C','Val4'],  specific=[],  # (shared electrons = 5-8)
                             single=[0,1,2,3], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[-1])
# examples for C2sc: [CH2-][N+]#N
atomTypes['C2d' ] = AtomType('C2d',  generic=['R','R!H','C','Val4'],  specific=[],  # (shared electrons = 6)
                             single=[0], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[0])
# examples for C2d: singlet[C]=C
atomTypes['C2dc'] = AtomType('C2dc', generic=['R','R!H','C','Val4'],  specific=[],  # (shared electrons = 7-8)
                             single=[0,1], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[-1])
# examples for C2dc: C=[C-][N+]#N, [CH-]=[N+]=O, [CH+]=C=[CH-]
atomTypes['C2tc'] = AtomType('C2tc',  generic=['R','R!H','C','Val4'],  specific=[],  # (shared electrons = 8)
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[1], charge=[-1])
# examples for C2tc: [C-]#[O+]

atomTypes['N'   ] = AtomType('N',    generic=['R','R!H','Val5'],      specific=['N0sc','N1s','N1sc','N1dc','N3s','N3sc','N3d','N3t','N3b','N5sc','N5dc','N5ddc','N5dddc','N5t','N5tc','N5b','N5bd'],
                             single=[], allDouble=[], rDouble=[], oDouble=[], sDouble=[], triple=[], benzene=[], lonePairs=[], charge=[])
atomTypes['N0sc'] = AtomType('N0sc',  generic=['R','R!H','N','Val5'],  specific=[],  # (shared electrons = 7-8)
                             single=[0,1], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[3], charge=[-2])
# examples for N0sc: [NH+]#[N+][N-2] with adjList 1 N u0 p0 c+1 {2,S} {3,T}; 2 H u0 p0 c0 {1,S}; 3 N u0 p0 c+1 {1,T} {4,S}; 4 N u0 p3 c-2 {3,S}
atomTypes['N1s' ] = AtomType('N1s',  generic=['R','R!H','N','Val5'],  specific=[],  # (shared electrons = 5-6)
                             single=[0,1], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[2], charge=[0])
# examples for N1s: closed shell N-N, closed shell NH
atomTypes['N1sc'] = AtomType('N1sc',  generic=['R','R!H','N','Val5'],  specific=[],  # (shared electrons = 6-8)
                             single=[0,1,2], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[-1])
# examples for N1sc: [NH-][S+]=C, [NH-][N+]#C
atomTypes['N1dc'] = AtomType('N1dc',  generic=['R','R!H','N','Val5'],  specific=[],  # (shared electrons = 8)
                             single=[0], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[-1])
# examples for N1dc: [N-]=[N+]=N terminal nitrogen on azide (two lone pairs), [N-]=[NH+], [N-]=[SH+]
atomTypes['N3s' ] = AtomType('N3s',  generic=['R','R!H','N','Val5'],  specific=[],  # (shared electrons = 5-8)
                             single=[0,1,2,3], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[1], charge=[0])
# examples for N3s: NH3, NH2, NH, N, C[NH]...
atomTypes['N3sc'] = AtomType('N3sc', generic=['R','R!H','N','Val5'],  specific=[],  # (shared electrons = 4-6)
                             single=[0,1,2], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[1], charge=[+1])
# examples for N3sc: !! N3sc should eventually be deleted, see #1206
atomTypes['N3d' ] = AtomType('N3d',  generic=['R','R!H','N','Val5'],  specific=[],  # (shared electrons = 7-8)
                             single=[0,1], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[0])
# examples for N3d: N=O, N=N, C=N, [O]N=O, [N]=O, [N]=C
atomTypes['N3t' ] = AtomType('N3t',  generic=['R','R!H','N','Val5'],  specific=[],  # (shared electrons = 8)
                             single=[0], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[1], benzene=[0], lonePairs=[1], charge=[0])
# examples for N3t: N2, N#C, N#[C], N#CC
atomTypes['N3b' ] = AtomType('N3b',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[2], lonePairs=[1], charge=[0])
# examples for N3b: Oxazole, Pyradine, Pyrazine, 1,3,5-Triazine, Benzimidazole, Purine
atomTypes['N5sc'] = AtomType('N5sc',  generic=['R','R!H','N','Val5'],  specific=[],  # (shared electrons = 4-8)
                             single=[0,1,2,3,4], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[0], charge=[+1,+2])
# examples for N5sc: !! N5sc should eventually be deleted, see #1206   [NH4+], [NH3+][O-] {N has u1 p0}
atomTypes['N5dc'] = AtomType('N5dc',  generic=['R','R!H','N','Val5'],  specific=[],  # (shared electrons = 6-8)
                             single=[0,1,2], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[+1])
# examples for N5dc: O[N+](=O)(O-) nitrate group, [N+](=O)(O)[O-], O=[N+][O-], [N+](=O)(O[N+](=O)[O-])[O-], C=[N+]=[SH-], [NH2+]=[SH-]
atomTypes['N5ddc'] = AtomType('N5ddc', generic=['R','R!H','N','Val5'],  specific=[],  # (shared electrons = 8)
                             single=[0], allDouble=[2], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[+1])
# examples for N5ddc: N=[N+]=[N-] center nitrogen on azide, [N-]=[N+]=O, C=[N+]=[SH-]
atomTypes['N5dddc'] = AtomType('N5dddc', generic=['R','R!H','N','Val5'],  specific=[],  # (shared electrons = 6)
                             single=[0], allDouble=[3], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[-1])
# examples for N5dddc: C=[N-](=C)=[NH2+]
atomTypes['N5t' ] = AtomType('N5t',   generic=['R','R!H','N','Val5'],  specific=[],  # (shared electrons = 8-10)
                             single=[0,1,2], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[1], benzene=[0], lonePairs=[0], charge=[0])
# examples for N5t: C#[NH2]
atomTypes['N5tc'] = AtomType('N5tc',  generic=['R','R!H','N','Val5'],  specific=[],  # (shared electrons = 7-8)
                             single=[0,1], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[1], benzene=[0], lonePairs=[0], charge=[+1])
# examples for N5tc: C[N+]#[C-] isocyano group, N#[N+][O-], [NH+]#[C-] (note that C- has p1 here), [N+]#[C-] (note that C- has p1 here), [O-][N+]#C (one of the res structures of Fulminic acid), C[N+]#[C-] (note that C- has p1 here)
atomTypes['N5b' ] = AtomType('N5b',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[2], lonePairs=[0], charge=[0,+1])
# examples for N5b: Pyrrole, Indole, Benzimidazole, Purine; Note that this is the only N atomType with valence 5 which isn't necessarily charged.
atomTypes['N5bd'] = AtomType('N5bd',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[2], lonePairs=[0], charge=[0])
# examples for N5bd: AdjList """1 N u0 p0 c0 {2,B} {6,B} {7,D} 2 C u0 p0 {1,B} {3,B} {8,S} 3 C u0 p0 {2,B} {4,B} {9,S} 4 C u0 p0 {3,B} {5,B} {10,S} 5 C u0 p0 {4,B} {6,B} {11,S} 6 N u0 p1 {1,B} {5,B} 7 O u0 p2 c0 {1,D} 8 H u0 p0 {2,S} 9 H u0 p0 {3,S} 10 H u0 p0 {4,S} 11 H u0 p0 {5,S}"""

atomTypes['O'   ] = AtomType('O',    generic=['R','R!H','Val6'],      specific=['Oa','O0sc','O2s','O2sc','O2d','O4sc','O4dc','O4tc','O4b'],
                             single=[], allDouble=[], rDouble=[], oDouble=[], sDouble=[], triple=[], benzene=[], lonePairs=[], charge=[])
atomTypes['Oa'  ] = AtomType('Oa',   generic=['R','R!H','O','Val6'],  specific=[],  # (shared electrons = 6)
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[3], charge=[0])
# examples for Oa: atomic oxygen (closed shell)
atomTypes['O0sc'] = AtomType('O0sc', generic=['R','R!H','O','Val6'],  specific=[],  # (shared electrons = 8)
                             single=[0,1], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[3], charge=[-1])
# examples for O0sc: Nitric acid O[N+](=O)([O-])
atomTypes['O2s' ] = AtomType('O2s',  generic=['R','R!H','O','Val6'],  specific=[],  # (shared electrons = 8)
                             single=[0,1,2], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[0])
# examples for O2s: H2O, OH, CH3OH
atomTypes['O2sc'] = AtomType('O2sc', generic=['R','R!H','O','Val6'],  specific=[],  # (shared electrons = 6)
                             single=[0,1], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[+1])
# examples for O2sc: C=[S-][O+]
atomTypes['O2d' ] = AtomType('O2d',  generic=['R','R!H','O','Val6'],  specific=[],  # (shared electrons = 8)
                             single=[0], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[0])
# examples for O2d: CO2, CH2O
atomTypes['O4sc'] = AtomType('O4sc', generic=['R','R!H','O','Val6'],  specific=[],  # (shared electrons = 5-8)
                             single=[0,1,2,3], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[+1])
# examples for O4sc: [O-][OH+]C
atomTypes['O4dc'] = AtomType('O4dc', generic=['R','R!H','O','Val6'],  specific=[],  # (shared electrons = 7-8)
                             single=[0,1], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[+1])
# examples for O4dc: the positively charged O in ozone [O-][O+]=O
atomTypes['O4tc'] = AtomType('O4tc',  generic=['R','R!H','O','Val6'],  specific=[],  # (shared electrons = 8)
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[1], charge=[+1])
# examples for O4tc: [C-]#[O+]
atomTypes['O4b' ] = AtomType('O4b',  generic=['R','R!H','O','Val6'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[2], lonePairs=[1], charge=[0])
# examples for S4b: Furane, Benzofurane, Benzo[c]thiophene, Oxazole...

atomTypes['Si'  ] = AtomType('Si',   generic=['R','R!H','Val4'],      specific=['Sis','Sid','Sidd','Sit','SiO','Sib','Sibf'],
                             single=[], allDouble=[], rDouble=[], oDouble=[], sDouble=[], triple=[], benzene=[], lonePairs=[], charge=[])
atomTypes['Sis' ] = AtomType('Sis',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[], charge=[])
atomTypes['SiO' ] = AtomType('SiO',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[1], rDouble=[], oDouble=[1], sDouble=[], triple=[0], benzene=[0], lonePairs=[], charge=[])
atomTypes['Sid' ] = AtomType('Sid',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[1], rDouble=[], oDouble=[0], sDouble=[], triple=[0], benzene=[0], lonePairs=[], charge=[])
atomTypes['Sidd'] = AtomType('Sidd', generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[2], rDouble=[0,1,2], oDouble=[0,1,2], sDouble=[0,1,2], triple=[0], benzene=[0], lonePairs=[], charge=[])
atomTypes['Sit' ] = AtomType('Sit',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[], charge=[])
atomTypes['Sib' ] = AtomType('Sib',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[2], lonePairs=[], charge=[])
atomTypes['Sibf'] = AtomType('Sibf', generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[3], lonePairs=[], charge=[])

atomTypes['S'   ] = AtomType('S',    generic=['R','R!H','Val6'],      specific=['Sa','S0sc','S2s','S2sc','S2d','S2dc','S2tc','S4s','S4sc','S4d','S4dd','S4dc','S4b','S4t','S4tdc','S6s','S6sc','S6d','S6dd','S6ddd','S6dc','S6t','S6td','S6tt','S6tdc'],
                             single=[], allDouble=[], rDouble=[], oDouble=[], sDouble=[], triple=[], benzene=[], lonePairs=[], charge=[])
atomTypes['Sa'  ] = AtomType('Sa',   generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 6)
                            single=[0], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[3], charge=[0])
# examples for Sa: atomic sulfur (closed shell)
atomTypes['S0sc'] = AtomType('S0sc', generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 7-8)
                             single=[0,1], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[3], charge=[-1])
# examples for S0sc: [S-][S+]=S
atomTypes['S2s' ] = AtomType('S2s',  generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 6-8)
                             single=[0,1,2], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[2], charge=[0])
# examples for S2s: [S], [SH], S {H2S}, [S][S], SS {H2S2}, SSC, CSSC, SO {HSOH}...
atomTypes['S2sc'] = AtomType('S2sc', generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 7-10)
                            single=[0,1,2,3], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[2], charge=[-1,+1])
# examples for S2sc: [S-][S+], N#[N+][S-](O)O
atomTypes['S2d' ] = AtomType('S2d',  generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 8)
                             single=[0], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[0])
# examples for S2d: S=S, C=S, S=O, S=N, S=C=S, S=C=O, S=C=S...
atomTypes['S2dc'] = AtomType('S2dc', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1], allDouble=[1,2], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[-1])
# *Composite atomType; examples for S2dc: [SH-]=[N+]
atomTypes['S2tc'] = AtomType('S2tc', generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 10)
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[2], charge=[-1])
# examples for S2tc: [S-]#[NH+]
atomTypes['S4s' ] = AtomType('S4s',  generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 6-10)
                             single=[0,1,2,3,4], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[1], charge=[0])
# examples for S4s: H4S, SH3CH3...
atomTypes['S4sc'] = AtomType('S4sc', generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 5-8)
                             single=[0,1,2,3,4,5], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[1], charge=[-1,+1])
# examples for S4sc: CS[S+]([O-])C, O[SH..-][N+]#N
atomTypes['S4d' ] = AtomType('S4d',  generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 8-10)
                             single=[0,1,2], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[0])
# examples for S4d: O=S(O)O {Sulfurous acid}
atomTypes['S4dd'] = AtomType('S4dd', generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 10)
                             single=[0], allDouble=[2], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[0])
# examples for S4dd: O=S=O
atomTypes['S4dc'] = AtomType('S4dc', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2,3,4,5], allDouble=[1,2], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[-1,+1])
# *Composite atomType; examples for S4dc: [CH2-][S+]=C {where the [CH2-] has a lone pair}, [O+][S-](=O)=O, [O-][S+]=C, [NH-][S+]=C {where the [NH-] has two lone pairs}, [O-][S+]=O
atomTypes['S4b' ] = AtomType('S4b',  generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[2], lonePairs=[1], charge=[0])
# examples for S4b: Thiophene, Benzothiophene, Benzo[c]thiophene, Thiazole, Benzothiazole...
atomTypes['S4t' ] = AtomType('S4t',  generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 10)
                             single=[0,1], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[1], benzene=[0], lonePairs=[1], charge=[0])
# examples for S4t: C#S, C#SO, C#[S]
atomTypes['S4tdc'] = AtomType('S4tdc',generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2], allDouble=[0,1,2], rDouble=[], oDouble=[], sDouble=[], triple=[1,2], benzene=[0], lonePairs=[1], charge=[-1,+1])
# *Composite atomType; examples for S4tdc: [C-]#[S+]
atomTypes['S6s' ] = AtomType('S6s',  generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 6-12)
                             single=[0,1,2,3,4,5,6], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for S6s: H6S, F6S
atomTypes['S6sc'] = AtomType('S6sc', generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 7-14)
                             single=[0,1,2,3,4,5,6,7], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[0], charge=[-1,+1])
# examples for S6sc: O[SH3+][O-], [O-][S+3]([O-])[O-]
atomTypes['S6d' ] = AtomType('S6d',  generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 8-12)
                             single=[0,1,2,3,4], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for S6d: [SH4]=O, SF4=O, [SH4]=C, C[SH3]=C...
atomTypes['S6dd'] = AtomType('S6dd', generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 10-12)
                             single=[0,1,2], allDouble=[2], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for S6dd: S(=O)(=O)(O)O {H2SO4, Sulfuric acid}, Perfluorooctanesulfonic acid, Pyrosulfuric acid, Thiosulfuric acid {middle S}, OS(=O)(=O)OOS(=O)(=O)O
atomTypes['S6ddd'] = AtomType('S6ddd', generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 12)
                             single=[0], allDouble=[3], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for S6ddd: O=S(=O)(=O)
atomTypes['S6dc'] = AtomType('S6dc', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2,3,4,5], allDouble=[1,2,3], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[-1,+1])
# *Composite atomType; examples for S6dc: [CH-]=[SH3+], [CH-]=[SH2+]O, [CH-][SH2+], O=[S+](=O)[O-], [OH+]=[S-](=O)=O
atomTypes['S6t' ] = AtomType('S6t', generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 9-12)
                             single=[0,1,2,3], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[0], charge=[0])
# examples for S6t: H3S#N
atomTypes['S6td'] = AtomType('S6td', generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 11-12)
                             single=[0,1], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[0], charge=[0])
# examples for S6td: HS(=O)#N
atomTypes['S6tt'] = AtomType('S6tt', generic=['R','R!H','S','Val6'],  specific=[],  # (shared electrons = 12)
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[2], benzene=[0], lonePairs=[0], charge=[0])
# examples for S6tt: N#S#N
atomTypes['S6tdc'] = AtomType('S6tdc',generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2,3,4], allDouble=[0,1,2], rDouble=[], oDouble=[], sDouble=[], triple=[1,2], benzene=[0], lonePairs=[0], charge=[-1,+1])
# *Composite atomType; examples for S6tdc: [SH2+]#[C-], [N-]=[S+]#N

atomTypes['Cl'  ] = AtomType('Cl',   generic=['R','R!H','Val7'],      specific=['Cl1s'])
atomTypes['Cl1s'] = AtomType('Cl1s', generic=['R','R!H','Cl','Val7'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[3], charge=[0])
# examples for Cl1s: HCl, [Cl]

atomTypes['I'  ] = AtomType('I',   generic=['R','R!H','Val7'],      specific=['I1s'])
atomTypes['I1s'] = AtomType('I1s', generic=['R','R!H','I','Val7'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[3], charge=[0])
# examples for I1s: HI, [I], IO, CH3I, I2

atomTypes['R'   ].setActions(incrementBond=['R'],            decrementBond=['R'],            formBond=['R'],         breakBond=['R'],         incrementRadical=['R'],    decrementRadical=['R'],    incrementLonePair=['R'],   decrementLonePair=['R'])
atomTypes['R!H' ].setActions(incrementBond=['R!H'],          decrementBond=['R!H'],          formBond=['R!H'],       breakBond=['R!H'],       incrementRadical=['R!H'],  decrementRadical=['R!H'],  incrementLonePair=['R!H'], decrementLonePair=['R!H'])
atomTypes['Val4'].setActions(incrementBond=['Val4'],         decrementBond=['Val4'],         formBond=['Val4'],      breakBond=['Val4'],      incrementRadical=['Val4'], decrementRadical=['Val4'], incrementLonePair=['Val4'],decrementLonePair=['Val4'])
atomTypes['Val5'].setActions(incrementBond=['Val5'],         decrementBond=['Val5'],         formBond=['Val5'],      breakBond=['Val5'],      incrementRadical=['Val5'], decrementRadical=['Val5'], incrementLonePair=['Val5'],decrementLonePair=['Val5'])
atomTypes['Val6'].setActions(incrementBond=['Val6'],         decrementBond=['Val6'],         formBond=['Val6'],      breakBond=['Val6'],      incrementRadical=['Val6'], decrementRadical=['Val6'], incrementLonePair=['Val6'],decrementLonePair=['Val6'])
atomTypes['Val7'].setActions(incrementBond=['Val7'],         decrementBond=['Val7'],         formBond=['Val7'],      breakBond=['Val7'],      incrementRadical=['Val7'], decrementRadical=['Val7'], incrementLonePair=['Val7'],decrementLonePair=['Val7'])

atomTypes['H'   ].setActions(incrementBond=[],               decrementBond=[],               formBond=['H'],         breakBond=['H'],         incrementRadical=['H'],    decrementRadical=['H'],    incrementLonePair=[],      decrementLonePair=[])

atomTypes['He'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=['He'],   decrementRadical=['He'],   incrementLonePair=[],      decrementLonePair=[])
atomTypes['Ne'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=['Ne'],   decrementRadical=['Ne'],   incrementLonePair=[],      decrementLonePair=[])
atomTypes['Ar'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])

atomTypes['C'   ].setActions(incrementBond=['C'],            decrementBond=['C'],            formBond=['C'],         breakBond=['C'],         incrementRadical=['C'],    decrementRadical=['C'],    incrementLonePair=['C'],   decrementLonePair=['C'])
atomTypes['Ca'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['C2s'])
atomTypes['Cs'  ].setActions(incrementBond=['Cd','CO','CS'], decrementBond=[],               formBond=['Cs','Csc'],  breakBond=['Cs'],        incrementRadical=['Cs'],   decrementRadical=['Cs'],   incrementLonePair=['C2s'], decrementLonePair=['C2s'])
atomTypes['Csc' ].setActions(incrementBond=['Cdc'],          decrementBond=[],               formBond=['Csc'],       breakBond=['Csc','Cs'],  incrementRadical=['Csc'],  decrementRadical=['Csc'],  incrementLonePair=['C2sc'],decrementLonePair=['C2sc'])
atomTypes['Cd'  ].setActions(incrementBond=['Cdd','Ct','C2tc'],decrementBond=['Cs'],         formBond=['Cd','Cdc'],  breakBond=['Cd'],        incrementRadical=['Cd'],   decrementRadical=['Cd'],   incrementLonePair=['C2d'], decrementLonePair=[])
atomTypes['Cdc' ].setActions(incrementBond=[],               decrementBond=['Csc'],          formBond=['Cdc'],       breakBond=['Cdc','Cd','CO','CS'],incrementRadical=['Cdc'],decrementRadical=['Cdc'],incrementLonePair=['C2dc'],decrementLonePair=[])
atomTypes['CO'  ].setActions(incrementBond=['Cdd','C2tc'],   decrementBond=['Cs'],           formBond=['CO','Cdc'],  breakBond=['CO'],        incrementRadical=['CO'],   decrementRadical=['CO'],   incrementLonePair=['C2d'], decrementLonePair=[])
atomTypes['CS'  ].setActions(incrementBond=['Cdd','C2tc'],   decrementBond=['Cs'],           formBond=['CS','Cdc'],  breakBond=['CS'],        incrementRadical=['CS'],   decrementRadical=['CS'],   incrementLonePair=['C2d'], decrementLonePair=[])
atomTypes['Cdd' ].setActions(incrementBond=[],               decrementBond=['Cd','CO','CS'], formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['Ct'  ].setActions(incrementBond=[],               decrementBond=['Cd','CO','CS'], formBond=['Ct'],        breakBond=['Ct'],        incrementRadical=['Ct'],   decrementRadical=['Ct'],   incrementLonePair=['C2tc'],decrementLonePair=[])
atomTypes['Cb'  ].setActions(incrementBond=['Cbf'],          decrementBond=[],               formBond=['Cb'],        breakBond=['Cb'],        incrementRadical=['Cb'],   decrementRadical=['Cb'],   incrementLonePair=[],      decrementLonePair=[])
atomTypes['Cbf' ].setActions(incrementBond=[],               decrementBond=['Cb'],           formBond=[],            breakBond=['Cb'],        incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['C2s' ].setActions(incrementBond=['C2d'],          decrementBond=[],               formBond=['C2s'],       breakBond=['C2s'],       incrementRadical=['C2s'],  decrementRadical=['C2s'],  incrementLonePair=['Ca'],  decrementLonePair=['Cs'])
atomTypes['C2sc'].setActions(incrementBond=['C2dc'],         decrementBond=[],               formBond=['C2sc'],      breakBond=['C2sc'],      incrementRadical=['C2sc'], decrementRadical=['C2sc'], incrementLonePair=[],      decrementLonePair=['Cs'])
atomTypes['C2d' ].setActions(incrementBond=['C2tc'],         decrementBond=['C2s'],          formBond=['C2dc'],      breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['Cd','CO','CS'])
atomTypes['C2dc'].setActions(incrementBond=[],               decrementBond=['C2sc'],         formBond=[],            breakBond=['C2d'],       incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['Cdc'])
atomTypes['C2tc'].setActions(incrementBond=[],               decrementBond=['C2d'],          formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['Ct'])

atomTypes['N'   ].setActions(incrementBond=['N'],            decrementBond=['N'],            formBond=['N'],         breakBond=['N'],         incrementRadical=['N'],    decrementRadical=['N'],    incrementLonePair=['N'],   decrementLonePair=['N'])
atomTypes['N0sc'].setActions(incrementBond=[],               decrementBond=[],               formBond=['N0sc'],      breakBond=['N0sc'],      incrementRadical=['N0sc'], decrementRadical=['N0sc'], incrementLonePair=[],      decrementLonePair=['N1s','N1sc'])
atomTypes['N1s' ].setActions(incrementBond=['N1dc'],         decrementBond=[],               formBond=['N1s'],       breakBond=['N1s'],       incrementRadical=['N1s'],  decrementRadical=['N1s'],  incrementLonePair=['N0sc'],decrementLonePair=['N3s','N3sc'])
atomTypes['N1sc'].setActions(incrementBond=[],               decrementBond=[],               formBond=['N1sc'],      breakBond=['N1sc'],      incrementRadical=['N1sc'], decrementRadical=['N1sc'], incrementLonePair=[],      decrementLonePair=['N3s','N3sc'])
atomTypes['N1dc'].setActions(incrementBond=['N1dc'],         decrementBond=['N1s','N1dc'],   formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['N3d'])
atomTypes['N3s' ].setActions(incrementBond=['N3d'],          decrementBond=[],               formBond=['N3s'],       breakBond=['N3s'],       incrementRadical=['N3s'],  decrementRadical=['N3s'],  incrementLonePair=['N1s','N1sc'],decrementLonePair=['N5sc'])
atomTypes['N3sc'].setActions(incrementBond=['N3d'],          decrementBond=[],               formBond=['N3sc'],      breakBond=['N3sc'],      incrementRadical=['N3sc'], decrementRadical=['N3sc'], incrementLonePair=['N1s','N1sc'],decrementLonePair=['N5sc'])
atomTypes['N3d' ].setActions(incrementBond=['N3t'],          decrementBond=['N3s','N3sc'],   formBond=['N3d'],       breakBond=['N3d'],       incrementRadical=['N3d'],  decrementRadical=['N3d'],  incrementLonePair=['N1dc'],decrementLonePair=['N5dc'])
atomTypes['N3t' ].setActions(incrementBond=[],               decrementBond=['N3d'],          formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['N5t'])
atomTypes['N3b' ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['N5sc'].setActions(incrementBond=['N5dc'],         decrementBond=[],               formBond=['N5sc'],      breakBond=['N5sc'],      incrementRadical=['N5sc'], decrementRadical=['N5sc'], incrementLonePair=['N3s','N3sc'],decrementLonePair=[])
atomTypes['N5dc'].setActions(incrementBond=['N5ddc','N5tc','N5t'],decrementBond=['N5sc'],    formBond=['N5dc'],      breakBond=['N5dc'],      incrementRadical=['N5dc'], decrementRadical=['N5dc'], incrementLonePair=['N3d'], decrementLonePair=[])
atomTypes['N5ddc'].setActions(incrementBond=['N5dddc'],      decrementBond=['N5dc'],         formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['N5dddc'].setActions(incrementBond=[],             decrementBond=['N5ddc'],        formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['N5t' ].setActions(incrementBond=[],               decrementBond=['N5dc'],         formBond=['N5t'],       breakBond=['N5t'],       incrementRadical=['N5t'],  decrementRadical=['N5t'],  incrementLonePair=['N3t'], decrementLonePair=[])
atomTypes['N5tc'].setActions(incrementBond=[],               decrementBond=['N5dc'],         formBond=['N5tc'],      breakBond=['N5tc'],      incrementRadical=['N5tc'], decrementRadical=['N5tc'], incrementLonePair=[],      decrementLonePair=[])
atomTypes['N5b' ].setActions(incrementBond=['N5bd'],         decrementBond=[],               formBond=['N5b'],       breakBond=['N5b'],       incrementRadical=['N5b'],  decrementRadical=['N5b'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['N5bd'].setActions(incrementBond=[],               decrementBond=['N5b'],          formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])

atomTypes['O'   ].setActions(incrementBond=['O'],            decrementBond=['O'],            formBond=['O'],         breakBond=['O'],         incrementRadical=['O'],    decrementRadical=['O'],    incrementLonePair=['O'],   decrementLonePair=['O'])
atomTypes['Oa'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=['O0sc'],      breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['O2s','O2sc'])
atomTypes['O0sc'].setActions(incrementBond=[],               decrementBond=[],               formBond=['O0sc'],      breakBond=['Oa','O0sc'], incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['O2s','O2sc'])
atomTypes['O2s' ].setActions(incrementBond=['O2d'],          decrementBond=[],               formBond=['O2s','O2sc'],breakBond=['O2s'],       incrementRadical=['O2s'],  decrementRadical=['O2s'],  incrementLonePair=['Oa','O0sc'],decrementLonePair=['O4sc'])
atomTypes['O2sc'].setActions(incrementBond=['O2d'],          decrementBond=[],               formBond=[],            breakBond=['O2s'],       incrementRadical=['O2sc'], decrementRadical=['O2sc'], incrementLonePair=[],      decrementLonePair=['O4sc'])
atomTypes['O2d' ].setActions(incrementBond=[],               decrementBond=['O2s','O2sc'],   formBond=[],            breakBond=[],            incrementRadical=['O2d'],  decrementRadical=['O2d'],  incrementLonePair=[],      decrementLonePair=['O4dc','O4tc'])
atomTypes['O4sc'].setActions(incrementBond=['O4dc'],         decrementBond=[],               formBond=['O4sc'],      breakBond=['O4sc'],      incrementRadical=['O4sc'], decrementRadical=['O4sc'], incrementLonePair=['O2s','O2sc'],decrementLonePair=[])
atomTypes['O4dc'].setActions(incrementBond=['O4tc'],         decrementBond=['O4sc'],         formBond=['O4dc'],      breakBond=['O4dc'],      incrementRadical=['O4dc'], decrementRadical=['O4dc'], incrementLonePair=['O2d'], decrementLonePair=[])
atomTypes['O4tc'].setActions(incrementBond=[],               decrementBond=['O4dc'],         formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[], decrementLonePair=[])
atomTypes['O4b' ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])

atomTypes['Si'  ].setActions(incrementBond=['Si'],           decrementBond=['Si'],           formBond=['Si'],        breakBond=['Si'],        incrementRadical=['Si'],   decrementRadical=['Si'],   incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sis' ].setActions(incrementBond=['Sid','SiO'],    decrementBond=[],               formBond=['Sis'],       breakBond=['Sis'],       incrementRadical=['Sis'],  decrementRadical=['Sis'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sid' ].setActions(incrementBond=['Sidd','Sit'],   decrementBond=['Sis'],          formBond=['Sid'],       breakBond=['Sid'],       incrementRadical=['Sid'],  decrementRadical=['Sid'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sidd'].setActions(incrementBond=[],               decrementBond=['Sid','SiO'],    formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sit' ].setActions(incrementBond=[],               decrementBond=['Sid'],          formBond=['Sit'],       breakBond=['Sit'],       incrementRadical=['Sit'],  decrementRadical=['Sit'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['SiO' ].setActions(incrementBond=['Sidd'],         decrementBond=['Sis'],          formBond=['SiO'],       breakBond=['SiO'],       incrementRadical=['SiO'],  decrementRadical=['SiO'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sib' ].setActions(incrementBond=[],               decrementBond=[],               formBond=['Sib'],       breakBond=['Sib'],       incrementRadical=['Sib'],  decrementRadical=['Sib'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sibf'].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])

atomTypes['S'   ].setActions(incrementBond=['S'],            decrementBond=['S'],            formBond=['S'],         breakBond=['S'],         incrementRadical=['S'],    decrementRadical=['S'],    incrementLonePair=['S'],   decrementLonePair=['S'])
atomTypes['S0sc'].setActions(incrementBond=['S0sc'],         decrementBond=['S0sc'],         formBond=['S0sc'],      breakBond=['Sa','S0sc'], incrementRadical=['S0sc'], decrementRadical=['S0sc'], incrementLonePair=[],      decrementLonePair=['S2s','S2sc','S2dc','S2tc'])
atomTypes['Sa'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=['S0sc'],      breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['S2s'])
atomTypes['S2s' ].setActions(incrementBond=['S2d','S2dc'],   decrementBond=[],               formBond=['S2s','S2sc'],breakBond=['S2s'],       incrementRadical=['S2s'],  decrementRadical=['S2s'],  incrementLonePair=['Sa','S0sc'],decrementLonePair=['S4s','S4sc'])
atomTypes['S2sc'].setActions(incrementBond=['S2dc'],         decrementBond=[],               formBond=['S2sc'],      breakBond=['S2sc','S2s'],incrementRadical=['S2sc'], decrementRadical=['S2sc'], incrementLonePair=['S0sc'],decrementLonePair=['S4s','S4sc'])
atomTypes['S2d' ].setActions(incrementBond=['S2tc'],         decrementBond=['S2s'],          formBond=['S2d'],       breakBond=['S2d'],       incrementRadical=['S2d'],  decrementRadical=['S2d'],  incrementLonePair=[],      decrementLonePair=['S4dc','S4d'])
atomTypes['S2dc'].setActions(incrementBond=['S2tc','S2dc'],  decrementBond=['S2sc','S2s','S2dc'],formBond=['S2dc'],  breakBond=['S2dc'],      incrementRadical=['S2dc'], decrementRadical=['S2dc'], incrementLonePair=['S0sc'],decrementLonePair=['S4d','S4dc'])
atomTypes['S2tc'].setActions(incrementBond=[],               decrementBond=['S2d','S2dc'],   formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=['S0sc'],decrementLonePair=['S4t'])
atomTypes['S4s' ].setActions(incrementBond=['S4d','S4dc'],   decrementBond=[],               formBond=['S4s'],       breakBond=['S4s'],       incrementRadical=['S4s'],  decrementRadical=['S4s'],  incrementLonePair=['S2s','S2sc'],decrementLonePair=['S6s'])
atomTypes['S4sc'].setActions(incrementBond=['S4d','S4dc'],   decrementBond=[],               formBond=['S4s','S4sc'],breakBond=['S4sc'],      incrementRadical=['S4sc'], decrementRadical=['S4sc'], incrementLonePair=['S2s','S2sc'],decrementLonePair=['S6s'])
atomTypes['S4d' ].setActions(incrementBond=['S4dd','S4dc','S4t','S4tdc'],decrementBond=['S4s','S4sc'],formBond=['S4dc','S4d'],breakBond=['S4d','S4dc'],incrementRadical=['S4d'],decrementRadical=['S4d'],incrementLonePair=['S2d','S2dc'],decrementLonePair=['S6d','S6dc'])
atomTypes['S4dc'].setActions(incrementBond=['S4dd','S4dc','S4tdc'],decrementBond=['S4sc','S4dc'],formBond=['S4d','S4dc'],breakBond=['S4d','S4dc'],incrementRadical=['S4dc'],decrementRadical=['S4dc'],incrementLonePair=['S2d','S2dc'],decrementLonePair=['S6d','S6dc'])
atomTypes['S4b' ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['S4dd'].setActions(incrementBond=['S4dc'],         decrementBond=['S4dc','S4d'],   formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['S6dd'])
atomTypes['S4t' ].setActions(incrementBond=[],               decrementBond=['S4d'],          formBond=['S4t'],       breakBond=['S4t'],       incrementRadical=['S4t'],  decrementRadical=['S4t'],  incrementLonePair=['S2tc'],decrementLonePair=['S6t','S6tdc'])
atomTypes['S4tdc'].setActions(incrementBond=['S4tdc'],       decrementBond=['S4d','S4tdc'],  formBond=['S4tdc'],     breakBond=['S4tdc'],     incrementRadical=['S4tdc'],decrementRadical=['S4tdc'],incrementLonePair=['S6tdc'],decrementLonePair=['S6td','S6tdc'])
atomTypes['S6s' ].setActions(incrementBond=['S6d','S6dc'],   decrementBond=[],               formBond=['S6s'],       breakBond=['S6s'],       incrementRadical=['S6s'],  decrementRadical=['S6s'],  incrementLonePair=['S4s','S4sc'],decrementLonePair=[])
atomTypes['S6sc'].setActions(incrementBond=['S6dc'],         decrementBond=[],               formBond=['S6sc'],      breakBond=['S6sc'],      incrementRadical=['S6sc'], decrementRadical=['S6sc'], incrementLonePair=['S4s','S4sc'],decrementLonePair=[])
atomTypes['S6d' ].setActions(incrementBond=['S6dd','S6t','S6tdc'],decrementBond=['S6s'],     formBond=['S6d','S6dc'],breakBond=['S6d','S6dc'],incrementRadical=['S6d'],  decrementRadical=['S6d'],  incrementLonePair=['S4d','S4dc'], decrementLonePair=[])
atomTypes['S6dc'].setActions(incrementBond=['S6dd','S6ddd','S6dc','S6t','S6td','S6tdc'],decrementBond=['S6sc','S6dc'],formBond=['S6d','S6dc'],breakBond=['S6d','S6dc'],incrementRadical=['S6dc'],decrementRadical=['S6dc'],incrementLonePair=['S4d','S4dc'],decrementLonePair=[])
atomTypes['S6dd'].setActions(incrementBond=['S6ddd','S6td'], decrementBond=['S6d','S6dc'],   formBond=['S6dd','S6dc'],breakBond=['S6dd'],     incrementRadical=['S6dd'], decrementRadical=['S6dd'], incrementLonePair=['S4dd'],decrementLonePair=[])
atomTypes['S6ddd'].setActions(incrementBond=[],              decrementBond=['S6dd','S6dc'],  formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['S6t' ].setActions(incrementBond=['S6td'],         decrementBond=['S6d','S6dc'],   formBond=['S6t'],       breakBond=['S6t'],       incrementRadical=['S6t'],  decrementRadical=['S6t'],  incrementLonePair=['S4t'], decrementLonePair=[])
atomTypes['S6td'].setActions(incrementBond=['S6tt','S6tdc'], decrementBond=['S6dc','S6t','S6dd','S6tdc'],formBond=['S6td'],breakBond=['S6td'],incrementRadical=['S6td'], decrementRadical=['S6td'], incrementLonePair=['S4tdc'],decrementLonePair=[])
atomTypes['S6tt'].setActions(incrementBond=[],               decrementBond=['S6td','S6tdc'], formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['S6tdc'].setActions(incrementBond=['S6td','S6tdc','S6tt'],decrementBond=['S6dc','S6tdc'],formBond=['S6tdc'],breakBond=['S6tdc'],    incrementRadical=['S6tdc'],decrementRadical=['S6tdc'],incrementLonePair=['S4t','S4tdc'],decrementLonePair=[])

atomTypes['Cl'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=['Cl'],        breakBond=['Cl'],        incrementRadical=['Cl'],   decrementRadical=['Cl'],   incrementLonePair=[],      decrementLonePair=[])
atomTypes['Cl1s'].setActions(incrementBond=[],               decrementBond=[],               formBond=['Cl1s'],      breakBond=['Cl1s'],      incrementRadical=['Cl1s'], decrementRadical=['Cl1s'], incrementLonePair=[],      decrementLonePair=[])

atomTypes['I'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=['I'],        breakBond=['I'],        incrementRadical=['I'],   decrementRadical=['I'],   incrementLonePair=[],      decrementLonePair=[])
atomTypes['I1s'].setActions(incrementBond=[],               decrementBond=[],               formBond=['I1s'],      breakBond=['I1s'],      incrementRadical=['I1s'], decrementRadical=['I1s'], incrementLonePair=[],      decrementLonePair=[])

#list of elements that do not have more specific atomTypes
#these are ordered on priority of picking if we encounter a more general atomType for make
allElements=['H', 'C', 'O', 'N', 'S', 'Si', 'Cl', 'Ne', 'Ar', 'He',]
nonSpecifics=['H', 'He', 'Ne', 'Ar',]

for atomType in atomTypes.values():
    for items in [atomType.generic, atomType.specific,
      atomType.incrementBond, atomType.decrementBond, atomType.formBond,
      atomType.breakBond, atomType.incrementRadical, atomType.decrementRadical, atomType.incrementLonePair, atomType.decrementLonePair]:
        for index in range(len(items)):
            items[index] = atomTypes[items[index]]


def getFeatures(atom, bonds):
    """
    Returns a list of features needed to determine atomType for :class:'Atom'
    or :class:'GroupAtom' object 'atom and with local bond structure `bonds`,
    a ``dict`` containing atom-bond pairs.

    """
    cython.declare(single=cython.int, allDouble=cython.int, rDouble=cython.int,
                   sDouble=cython.int, oDouble=cython.int, triple=cython.int,
                   benzene=cython.int)
    cython.declare(features=cython.list)

    # Count numbers of each higher-order bond type
    single = 0; rDouble = 0; oDouble = 0; sDouble = 0; triple = 0; benzene = 0
    for atom2, bond12 in bonds.iteritems():
        if bond12.isSingle():
            single += 1
        elif bond12.isDouble():
            if atom2.isOxygen():
                oDouble += 1
            elif atom2.isSulfur():
                sDouble += 1
            else:
                # rDouble is for double bonds NOT to oxygen or Sulfur
                rDouble += 1
        elif bond12.isTriple(): triple += 1
        elif bond12.isBenzene(): benzene += 1

    # allDouble is for all double bonds, to anything
    allDouble = rDouble + oDouble + sDouble
    features = [single, allDouble, rDouble, oDouble, sDouble, triple, benzene, atom.lonePairs, atom.charge]

    return features

def getAtomType(atom, bonds):
    """
    Determine the appropriate atom type for an :class:`Atom` object `atom`
    with local bond structure `bonds`, a ``dict`` containing atom-bond pairs.
    """

    cython.declare(atomSymbol=str)
    cython.declare(molFeatureList=cython.list, atomTypeFeatureList=cython.list)

    # Use element and counts to determine proper atom type
    atomSymbol = atom.symbol
    #These elements do not do not have a more specific atomType
    if atomSymbol in nonSpecifics:
        return atomTypes[atomSymbol]

    molFeatureList = getFeatures(atom, bonds)
    for specificAtomType in atomTypes[atomSymbol].specific:
        atomtypeFeatureList = specificAtomType.getFeatures()
        for molFeature, atomtypeFeature in zip(molFeatureList, atomtypeFeatureList):
            if atomtypeFeature == []:
                continue
            elif molFeature not in atomtypeFeature:
                break
        else:
            return specificAtomType
    else:
        single = molFeatureList[0]
        rDouble = molFeatureList[2]
        oDouble = molFeatureList[3]
        sDouble = molFeatureList[4]
        triple = molFeatureList[5]
        benzene = molFeatureList[6]
        lonePairs = molFeatureList[7]
        charge = molFeatureList[8]

        raise AtomTypeError(
            'Unable to determine atom type for atom {0}, which has {1:d} single bonds, {2:d} double bonds to C, {3:d} double bonds to O, {4:d} double bonds to S, {5:d} triple bonds, {6:d} benzene bonds, {7:d} lone pairs, and {8:d} charge.'.format(
                atom, single, rDouble, oDouble, sDouble, triple, benzene, lonePairs, charge))

