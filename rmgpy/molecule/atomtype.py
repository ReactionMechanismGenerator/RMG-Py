#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

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

    The following features are what are required in a given atomtype. Any int in the list is acceptable.
    An empty list is a wildcard
    'single'            ''list''            The total number of single bonds on the atom
    'allDouble'         ''list''            The total number of double bonds on the atom
    'rDouble'           ''list''            The number of double bonds to any non-oxygen, nonsulfur
    'oDouble'           ''list''            The number of double bonds to oxygen
    'sDouble'           ''list''            The number of double bonds to sulfur
    'triple'            ''list''            The total number of triple bonds on the atom
    'benzene'           ''list''            The total number of benzene bonds on the atom
    'lonePairs'         ''list''            The number of lone pairs on the atom
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
                 lonePairs=None,):
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
            'lonePairs': self.lonePairs
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
                  self.lonePairs,]
        return features

################################################################################


"""
Note: function to read adjacency lists assumes that all atom types begin
with a capital letter [A-Z]

For making sample atoms, we use the first atomtype under specific,
therefore the first one in the list should always be an element.
"""

atomTypes = {}
atomTypes['R']    = AtomType(label='R', generic=[], specific=[
    'H',
    'R!H',
    'Val4','Val5','Val6','Val7',    
    'He',
    'C','Cs','Cd','Cdd','Ct','CO','Cb','Cbf','CS',
    'N','N1sc','N1s','N1d','N2s','N3s','N3d','N3t','N3b','N5s','N5d','N5dd','N5t','N5b',
    'O','Os','Od','Oa','Ot',
    'Ne',
    'Si','Sis','Sid','Sidd','Sit','SiO','Sib','Sibf',
    'S','Ss','Sd','Sa','St',
    'Cl','Ar']
)
atomTypes['R!H']  = AtomType(label='R!H', generic=['R'], specific=[
    'He',
    'Val4','Val5','Val6','Val7',
    'C', 'Cs','Cd','Cdd','Ct','CO','Cb','Cbf','CS',
    'N','N1sc','N1s','N1d','N2s','N3s','N3d','N3t','N3b','N5s','N5d','N5dd','N5t','N5b',
    'O','Os','Od','Oa','Ot',
    'Ne',
    'Si','Sis','Sid','Sidd','Sit','SiO','Sib','Sibf',
    'S','Ss','Sd','Sa','St',
    'Cl','Ar'])

atomTypes['Val4'] = AtomType(label='Val4', generic=['R','R!H'], specific=[
    'C','Cs','Cd','Cdd','Ct','CO','Cb','Cbf','CS',
    'Si','Sis','Sid','Sidd','Sit','SiO','Sib','Sibf'])

atomTypes['Val5'] = AtomType(label='Val5', generic=['R','R!H'], specific=[
    'N','N1sc','N1s','N1d','N2s','N3s','N3d','N3t','N3b','N5s','N5d','N5dd','N5t','N5b'])

atomTypes['Val6'] = AtomType(label='Val6', generic=['R','R!H'], specific=[
    'O','Os','Od','Oa','Ot',
    'S','Ss','Sd','Sa', 'St'])

atomTypes['Val7'] = AtomType(label='Val7', generic=['R','R!H'], specific=[
    'Cl'])

atomTypes['H'   ] = AtomType('H',    generic=['R'],            specific=[])

atomTypes['He'   ] = AtomType('He',  generic=['R','R!H'],      specific=[])

atomTypes['C'   ] = AtomType('C',    generic=['R','R!H','Val4'],      specific=['Cs','Cd','Cdd','Ct','CO','Cb','Cbf','CS'],
                             single=[], allDouble=[], rDouble=[], oDouble=[], sDouble=[], triple=[], benzene=[],)
atomTypes['Cs'  ] = AtomType('Cs',   generic=['R','R!H','C','Val4'],  specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0],)
atomTypes['Cd'  ] = AtomType('Cd',   generic=['R','R!H','C','Val4'],  specific=[],
                             single=[], allDouble=[1], rDouble=[1], oDouble=[0], sDouble=[0], triple=[0], benzene=[0])
atomTypes['Cdd' ] = AtomType('Cdd',  generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0], allDouble=[2], rDouble=[0,1,2], oDouble=[0,1,2], sDouble=[0,1,2], triple=[0], benzene=[0])
atomTypes['Ct'  ] = AtomType('Ct',   generic=['R','R!H','C','Val4'],  specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0])
atomTypes['CO'  ] = AtomType('CO',   generic=['R','R!H','C','Val4'],  specific=[],
                             single=[], allDouble=[1], rDouble=[0], oDouble=[1], sDouble=[0], triple=[0], benzene=[0])
atomTypes['Cb'  ] = AtomType('Cb',   generic=['R','R!H','C','Val4'],  specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[0], triple=[0], benzene=[1,2])
atomTypes['Cbf' ] = AtomType('Cbf',  generic=['R','R!H','C','Val4'],  specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[3])
atomTypes['CS'  ] = AtomType('CS',   generic=['R','R!H','C','Val4'],  specific=[],
                             single=[], allDouble=[1], rDouble=[0], oDouble=[0], sDouble=[1], triple=[0], benzene=[0])

atomTypes['N'   ] = AtomType('N',    generic=['R','R!H','Val5'],      specific=['N1sc','N1s','N1d','N2s','N3s','N3d','N3t','N3b','N5s','N5d','N5dd','N5t','N5b'],
                             single=[], allDouble=[], rDouble=[], oDouble=[], sDouble=[], triple=[], benzene=[])
atomTypes['N1sc'] = AtomType('N1sc',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[3])
    #examples for N1sc: [NH+]#[N+][N-2] with adjList 1 N u0 p0 c+1 {2,S} {3,T}; 2 H u0 p0 c0 {1,S}; 3 N u0 p0 c+1 {1,T} {4,S}; 4 N u0 p3 c-2 {3,S}
atomTypes['N1s' ] = AtomType('N1s',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[2])
    #examples for N1s: closed shell N-N, closed shell NH
atomTypes['N1d' ] = AtomType('N1d',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2])
    #examples for N1d: [N-]=[N+]=N terminal nitrogen on azide (two lone pairs), [N-]=[NH+]
atomTypes['N2s' ] = AtomType('N2s',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[1,2], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2])
    #examples for N2s: [NH-][S+]=C, [NH-][N+]#C
atomTypes['N3s' ] = AtomType('N3s',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1,2,3], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[1])
    #examples for N3s: NH3, NH2, NH, N, C[NH]...
atomTypes['N3d' ] = AtomType('N3d',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1])
    #examples for N3d: N=O, N=N, C=N, [O]N=O, [N]=O, [N]=C
atomTypes['N3t' ] = AtomType('N3t',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[1], benzene=[0], lonePairs=[1])
    #examples for N3t: N2, N#C, N#[C], N#CC
atomTypes['N3b' ] = AtomType('N3b',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[2], lonePairs=[1])
    #examples for N3b: Oxazole, Pyradine, Pyrazine, 1,3,5-Triazine, Benzimidazole, Purine
atomTypes['N5s' ] = AtomType('N5s',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1,2,3,4], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[0])
    #examples for N5s: [NH4+], [NH2+][O-] {N has u1 p0}, [NH3+][O-]
atomTypes['N5d' ] = AtomType('N5d',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1,2], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0])
    #examples for N5d: O[N+](=O)(O-) nitrate group, [N+](=O)(O)[O-], O=[N+][O-], [N+](=O)(O[N+](=O)[O-])[O-]
atomTypes['N5dd'] = AtomType('N5dd', generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0], allDouble=[2], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0])
    #examples for N5dd: N=[N+]=[N-] center nitrogen on azide, [N-]=[N+]=O
atomTypes['N5t' ] = AtomType('N5t',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[1], benzene=[0], lonePairs=[0])
    #examples for N5t: C[N+]#[C-] isocyano group, N#[N+][O-], [NH+]#[C-] (note that C- has p1 here), [N+]#[C-] (note that C- has p1 here), [O-][N+]#C, C[N+]#[C-] (note that C- has p1 here)
atomTypes['N5b' ] = AtomType('N5b',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[2], lonePairs=[0])
    #examples for N5b: Pyrrole, Indole, Benzimidazole, Purine (or any of these examples where the H is abstracted from N and leaves a radical?)

atomTypes['O'   ] = AtomType('O',    generic=['R','R!H','Val6'],      specific=['Os','Od','Oa','Ot'])
atomTypes['Os'  ] = AtomType('Os',   generic=['R','R!H','O','Val6'],  specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0])
atomTypes['Od'  ] = AtomType('Od',   generic=['R','R!H','O','Val6'],  specific=[],
                             single=[], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[], benzene=[])
atomTypes['Oa'  ] = AtomType('Oa',   generic=['R','R!H','O','Val6'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0])
atomTypes['Ot'  ] = AtomType('Ot',   generic=['R','R!H','O','Val6'],  specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0])

atomTypes['Ne'  ] = AtomType('Ne',   generic=['R','R!H'],      specific=[])
atomTypes['Si'  ] = AtomType('Si',   generic=['R','R!H','Val4'],      specific=['Sis','Sid','Sidd','Sit','SiO','Sib','Sibf'])
atomTypes['Sis' ] = AtomType('Sis',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0])
atomTypes['SiO' ] = AtomType('SiO',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[1], rDouble=[], oDouble=[1], sDouble=[], triple=[0], benzene=[0])
atomTypes['Sid' ] = AtomType('Sid',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[1], rDouble=[], oDouble=[0], sDouble=[], triple=[0], benzene=[0])
atomTypes['Sidd'] = AtomType('Sidd', generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[2], rDouble=[0,1,2], oDouble=[0,1,2], sDouble=[0,1,2], triple=[0], benzene=[0])
atomTypes['Sit' ] = AtomType('Sit',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0])
atomTypes['Sib' ] = AtomType('Sib',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[2])
atomTypes['Sibf'] = AtomType('Sibf', generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[3])

atomTypes['S'   ] = AtomType('S',    generic=['R','R!H','Val6'],      specific=['Ss','Sd','Sa','St'])
atomTypes['Ss'  ] = AtomType('Ss',   generic=['R','R!H','S','Val6'],  specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0])
atomTypes['Sd'  ] = AtomType('Sd',   generic=['R','R!H','S','Val6'],  specific=[],
                             single=[], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0])
atomTypes['Sa'  ] = AtomType('Sa',   generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0])
atomTypes['St'  ] = AtomType('St',   generic=['R','R!H','S','Val6'],  specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0])


atomTypes['Cl'  ] = AtomType('Cl',   generic=['R','R!H','Val7'],      specific=[])

atomTypes['Ar'  ] = AtomType('Ar',   generic=['R','R!H'],      specific=[])

atomTypes['R'   ].setActions(incrementBond=['R'],            decrementBond=['R'],            formBond=['R'],         breakBond=['R'],         incrementRadical=['R'],    decrementRadical=['R'],    incrementLonePair=['R'],   decrementLonePair=['R'])
atomTypes['R!H' ].setActions(incrementBond=['R!H'],          decrementBond=['R!H'],          formBond=['R!H'],       breakBond=['R!H'],       incrementRadical=['R!H'],  decrementRadical=['R!H'],  incrementLonePair=['R!H'], decrementLonePair=['R!H'])
atomTypes['Val4'].setActions(incrementBond=['Val4'],         decrementBond=['Val4'],         formBond=['Val4'],      breakBond=['Val4'],      incrementRadical=['Val4'], decrementRadical=['Val4'], incrementLonePair=['Val4'],decrementLonePair=['Val4'])
atomTypes['Val5'].setActions(incrementBond=['Val5'],         decrementBond=['Val5'],         formBond=['Val5'],      breakBond=['Val5'],      incrementRadical=['Val5'], decrementRadical=['Val5'], incrementLonePair=['Val5'],decrementLonePair=['Val5'])
atomTypes['Val6'].setActions(incrementBond=['Val6'],         decrementBond=['Val6'],         formBond=['Val6'],      breakBond=['Val6'],      incrementRadical=['Val6'], decrementRadical=['Val6'], incrementLonePair=['Val6'],decrementLonePair=['Val6'])
atomTypes['Val7'].setActions(incrementBond=['Val7'],         decrementBond=['Val7'],         formBond=['Val7'],      breakBond=['Val7'],      incrementRadical=['Val7'], decrementRadical=['Val7'], incrementLonePair=['Val7'],decrementLonePair=['Val7'])


atomTypes['H'   ].setActions(incrementBond=[],               decrementBond=[],               formBond=['H'],         breakBond=['H'],         incrementRadical=['H'],    decrementRadical=['H'],    incrementLonePair=[],      decrementLonePair=[])

atomTypes['He'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=['He'],   decrementRadical=['He'],   incrementLonePair=[],      decrementLonePair=[])

atomTypes['C'   ].setActions(incrementBond=['C'],            decrementBond=['C'],            formBond=['C'],         breakBond=['C'],         incrementRadical=['C'],    decrementRadical=['C'],    incrementLonePair=['C'],      decrementLonePair=['C'])
atomTypes['Cs'  ].setActions(incrementBond=['Cd','CO','CS'], decrementBond=[],               formBond=['Cs'],        breakBond=['Cs'],        incrementRadical=['Cs'],   decrementRadical=['Cs'],   incrementLonePair=['Cs'],      decrementLonePair=['Cs'])
atomTypes['Cd'  ].setActions(incrementBond=['Cdd','Ct'],     decrementBond=['Cs'],           formBond=['Cd'],        breakBond=['Cd'],        incrementRadical=['Cd'],   decrementRadical=['Cd'],   incrementLonePair=['Cd'],      decrementLonePair=['Cd'])
atomTypes['Cdd' ].setActions(incrementBond=[],               decrementBond=['Cd','CO','CS'], formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['Ct'  ].setActions(incrementBond=[],               decrementBond=['Cd','CO','CS'], formBond=['Ct'],        breakBond=['Ct'],        incrementRadical=['Ct'],   decrementRadical=['Ct'],   incrementLonePair=['Ct'],  decrementLonePair=['Ct'])
atomTypes['CO'  ].setActions(incrementBond=['Cdd'],          decrementBond=['Cs'],           formBond=['CO'],        breakBond=['CO'],        incrementRadical=['CO'],   decrementRadical=['CO'],   incrementLonePair=['CO'],  decrementLonePair=['CO'])
atomTypes['CS'  ].setActions(incrementBond=['Cdd'],          decrementBond=['Cs'],           formBond=['CS'],        breakBond=['CS'],        incrementRadical=['CS'],   decrementRadical=['CS'],   incrementLonePair=['CS'],  decrementLonePair=['CS'])
atomTypes['Cb'  ].setActions(incrementBond=['Cb'],           decrementBond=['Cb'],           formBond=['Cb'],        breakBond=['Cb'],        incrementRadical=['Cb'],   decrementRadical=['Cb'],   incrementLonePair=[],      decrementLonePair=[])
atomTypes['Cbf' ].setActions(incrementBond=['Cbf'],          decrementBond=['Cb'],           formBond=[],            breakBond=['Cb'],        incrementRadical=['Cbf'],  decrementRadical=['Cbf'],  incrementLonePair=[],      decrementLonePair=[])

atomTypes['N'   ].setActions(incrementBond=['N'],            decrementBond=['N'],            formBond=['N'],         breakBond=['N'],         incrementRadical=['N'],    decrementRadical=['N'],    incrementLonePair=['N'],   decrementLonePair=['N'])
atomTypes['N1sc'].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=['N1s'],       incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['N1s','N2s'])
atomTypes['N1s' ].setActions(incrementBond=['N1d'],          decrementBond=[],               formBond=['N1s'],       breakBond=['N1s'],       incrementRadical=['N1s'],  decrementRadical=['N1s'],  incrementLonePair=['N1sc'],decrementLonePair=['N3s'])
atomTypes['N1d' ].setActions(incrementBond=[],               decrementBond=['N1s'],          formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['N3d'])
atomTypes['N2s' ].setActions(incrementBond=[],               decrementBond=[],               formBond=['N2s'],       breakBond=['N2s'],       incrementRadical=['N2s'],  decrementRadical=['N2s'],  incrementLonePair=[],      decrementLonePair=['N3s'])
atomTypes['N3s' ].setActions(incrementBond=['N3d'],          decrementBond=[],               formBond=['N3s'],       breakBond=['N3s'],       incrementRadical=['N3s'],  decrementRadical=['N3s'],  incrementLonePair=['N1s','N2s'],decrementLonePair=['N5s'])
atomTypes['N3d' ].setActions(incrementBond=['N3t'],          decrementBond=['N3s'],          formBond=['N3d'],       breakBond=['N3d'],       incrementRadical=['N3d'],  decrementRadical=['N3d'],  incrementLonePair=['N1d'], decrementLonePair=['N5d'])
atomTypes['N3t' ].setActions(incrementBond=[],               decrementBond=['N3d'],          formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['N5t'])
atomTypes['N3b' ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['N5s' ].setActions(incrementBond=['N5d'],          decrementBond=[],               formBond=['N5s'],       breakBond=['N5s'],       incrementRadical=['N5s'],  decrementRadical=['N5s'],  incrementLonePair=['N3s'], decrementLonePair=[])
atomTypes['N5d' ].setActions(incrementBond=['N5dd','N5t'],   decrementBond=['N5s'],          formBond=['N5d'],       breakBond=['N5d'],       incrementRadical=['N5d'],  decrementRadical=['N5d'],  incrementLonePair=['N3d'], decrementLonePair=[])
atomTypes['N5dd'].setActions(incrementBond=[],               decrementBond=['N5d'],          formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['N5t' ].setActions(incrementBond=[],               decrementBond=['N5d'],          formBond=['N5t'],       breakBond=['N5t'],       incrementRadical=['N5t'],  decrementRadical=['N5t'],  incrementLonePair=['N3t'], decrementLonePair=[])
atomTypes['N5b' ].setActions(incrementBond=[],               decrementBond=[],               formBond=['N5b'],       breakBond=['N5b'],       incrementRadical=['N5b'],  decrementRadical=['N5b'],  incrementLonePair=[],      decrementLonePair=[])

atomTypes['O'   ].setActions(incrementBond=['O'],            decrementBond=['O'],            formBond=['O'],         breakBond=['O'],         incrementRadical=['O'],    decrementRadical=['O'],    incrementLonePair=['O'],  decrementLonePair=['O'])
atomTypes['Os'  ].setActions(incrementBond=['Od'],           decrementBond=[],               formBond=['Os'],        breakBond=['Os'],        incrementRadical=['Os'],   decrementRadical=['Os'],   incrementLonePair=['Os'],  decrementLonePair=['Os'])
atomTypes['Od'  ].setActions(incrementBond=[],               decrementBond=['Os'],           formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=['Od'],  decrementLonePair=['Od'])
atomTypes['Oa'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['Ot'  ].setActions(incrementBond=[],               decrementBond=['Od'],           formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=['Ot'],  decrementLonePair=['Ot'])

atomTypes['Ne'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=['Ne'],   decrementRadical=['Ne'],   incrementLonePair=[],    decrementLonePair=[])

atomTypes['Si'  ].setActions(incrementBond=['Si'],           decrementBond=['Si'],           formBond=['Si'],        breakBond=['Si'],        incrementRadical=['Si'],   decrementRadical=['Si'],   incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sis' ].setActions(incrementBond=['Sid','SiO'],    decrementBond=[],               formBond=['Sis'],       breakBond=['Sis'],       incrementRadical=['Sis'],  decrementRadical=['Sis'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sid' ].setActions(incrementBond=['Sidd','Sit'],   decrementBond=['Sis'],          formBond=['Sid'],       breakBond=['Sid'],       incrementRadical=['Sid'],  decrementRadical=['Sid'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sidd'].setActions(incrementBond=[],               decrementBond=['Sid','SiO'],    formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sit' ].setActions(incrementBond=[],               decrementBond=['Sid'],          formBond=['Sit'],       breakBond=['Sit'],       incrementRadical=['Sit'],  decrementRadical=['Sit'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['SiO' ].setActions(incrementBond=['Sidd'],         decrementBond=['Sis'],          formBond=['SiO'],       breakBond=['SiO'],       incrementRadical=['SiO'],  decrementRadical=['SiO'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sib' ].setActions(incrementBond=[],               decrementBond=[],               formBond=['Sib'],       breakBond=['Sib'],       incrementRadical=['Sib'],  decrementRadical=['Sib'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sibf'].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])

atomTypes['S'   ].setActions(incrementBond=['S'],            decrementBond=['S'],            formBond=['S'],         breakBond=['S'],         incrementRadical=['S'],    decrementRadical=['S'],    incrementLonePair=['S'],   decrementLonePair=['S'])
atomTypes['Ss'  ].setActions(incrementBond=['Sd'],           decrementBond=[],               formBond=['Ss'],        breakBond=['Ss'],        incrementRadical=['Ss'],   decrementRadical=['Ss'],   incrementLonePair=['Ss'],  decrementLonePair=['Ss'])
atomTypes['Sd'  ].setActions(incrementBond=[],               decrementBond=['Ss'],           formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=['Sd'],  decrementLonePair=['Sd'])
atomTypes['Sa'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['St'  ].setActions(incrementBond=[],               decrementBond=['Sd'],           formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=['St'],  decrementLonePair=['St'])


atomTypes['Cl'  ].setActions(incrementBond=[],               decrementBond=['Cl'],           formBond=['Cl'],        breakBond=['Cl'],        incrementRadical=['Cl'],   decrementRadical=['Cl'],   incrementLonePair=[],      decrementLonePair=[])

atomTypes['Ar'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])

#list of elements that do not have more specific atomTypes
#these are ordered on priority of picking if we encounter a more general atomType for make
allElements=['H', 'C', 'O', 'N', 'S', 'Si', 'Cl', 'Ne', 'Ar', 'He',]
nonSpecifics=['H', 'He', 'Ne', 'Cl', 'Ar',]

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
    features = [single, allDouble, rDouble, oDouble, sDouble, triple, benzene, atom.lonePairs]

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

        raise AtomTypeError(
            'Unable to determine atom type for atom {0}, which has {1:d} single bonds, {2:d} double bonds to C, {3:d} double bonds to O, {4:d} double bonds to S, {5:d} triple bonds, {6:d} benzene bonds, and {7:d} lone pairs.'.format(
                atom, single, rDouble, oDouble, sDouble, triple, benzene, lonePairs))

