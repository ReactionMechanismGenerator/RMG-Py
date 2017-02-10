#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
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
This module provides classes and methods for working with molecules and
molecular configurations. A molecule is represented internally using a graph
data type, where atoms correspond to vertices and bonds correspond to edges.
Both :class:`Atom` and :class:`Bond` objects store semantic information that
describe the corresponding atom or bond.
"""

import cython
import logging
import os
import re
import numpy
import urllib
from collections import OrderedDict
import itertools

import element as elements
try:
    import openbabel
except:
    pass
from rdkit import Chem
from .graph import Vertex, Edge, Graph, getVertexConnectivityValue
import rmgpy.molecule.group as gr
from .atomtype import AtomType, atomTypes, getAtomType, AtomTypeError
import rmgpy.constants as constants
import rmgpy.molecule.parser as parser
import rmgpy.molecule.generator as generator
import rmgpy.molecule.resonance as resonance
from .kekulize import kekulize

################################################################################

bond_orders = {'S': 1, 'D': 2, 'T': 3, 'B': 1.5}

globals().update({
    'bond_orders': bond_orders,
})

class Atom(Vertex):
    """
    An atom. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `atomType`          :class:`AtomType`   The :ref:`atom type <atom-types>`
    `element`           :class:`Element`    The chemical element the atom represents
    `radicalElectrons`  ``short``           The number of radical electrons
    `charge`            ``short``           The formal charge of the atom
    `label`             ``str``             A string label that can be used to tag individual atoms
    `coords`            ``numpy array``     The (x,y,z) coordinates in Angstrom
    `lonePairs`         ``short``           The number of lone electron pairs
    =================== =================== ====================================

    Additionally, the ``mass``, ``number``, and ``symbol`` attributes of the
    atom's element can be read (but not written) directly from the atom object,
    e.g. ``atom.symbol`` instead of ``atom.element.symbol``.
    """

    def __init__(self, element=None, radicalElectrons=0, charge=0, label='', lonePairs=-100, coords=numpy.array([])):
        Vertex.__init__(self)
        if isinstance(element, str):
            self.element = elements.__dict__[element]
        else:
            self.element = element
        self.radicalElectrons = radicalElectrons
        self.charge = charge
        self.label = label
        self.atomType = None
        self.lonePairs = lonePairs
        self.coords = coords

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return '{0}{1}{2}'.format(
            str(self.element),
            '.' * self.radicalElectrons,
            '+' * self.charge if self.charge > 0 else '-' * -self.charge,
        )

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "<Atom '{0}'>".format(str(self))

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        d = {
            'edges': self.edges,
            'connectivity1': self.connectivity1,
            'connectivity2': self.connectivity2,
            'connectivity3': self.connectivity3,
            'sortingLabel': self.sortingLabel,
            'atomType': self.atomType.label if self.atomType else None,
            'lonePairs': self.lonePairs,
        }
        return (Atom, (self.element.symbol, self.radicalElectrons, self.charge, self.label), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling an object.
        """
        self.edges = d['edges']
        self.connectivity1 = d['connectivity1']
        self.connectivity2 = d['connectivity2']
        self.connectivity3 = d['connectivity3']
        self.sortingLabel = d['sortingLabel']
        self.atomType = atomTypes[d['atomType']] if d['atomType'] else None
        self.lonePairs = d['lonePairs']
    
    @property
    def mass(self): return self.element.mass
    
    @property
    def number(self): return self.element.number

    @property
    def symbol(self): return self.element.symbol

    @property
    def bonds(self): return self.edges

    def equivalent(self, other):
        """
        Return ``True`` if `other` is indistinguishable from this atom, or
        ``False`` otherwise. If `other` is an :class:`Atom` object, then all
        attributes except `label` must match exactly. If `other` is an
        :class:`GroupAtom` object, then the atom must match any of the
        combinations in the atom pattern.
        """
        cython.declare(atom=Atom, ap=gr.GroupAtom)
        if isinstance(other, Atom):
            atom = other
            return (
                self.element                is atom.element and
                self.radicalElectrons       == atom.radicalElectrons   and
                self.lonePairs              == atom.lonePairs           and
                self.charge                 == atom.charge
                )
        elif isinstance(other, gr.GroupAtom):
            cython.declare(a=AtomType, radical=cython.short, lp=cython.short, charge=cython.short)
            ap = other
            for a in ap.atomType:
                if self.atomType.equivalent(a): break
            else:
                return False
            if ap.radicalElectrons:
                for radical in ap.radicalElectrons:
                    if self.radicalElectrons == radical: break
                else:
                    return False
            if ap.lonePairs:
                for lp in ap.lonePairs:
                    if self.lonePairs == lp: break
                else:
                    return False
            if ap.charge:
                for charge in ap.charge:
                    if self.charge == charge: break
                else:
                    return False
            return True
    
    def getDescriptor(self):
        return (self.getAtomConnectivityValue(), self.number)

    def getAtomConnectivityValue(self):
        return -1*self.connectivity

    def isSpecificCaseOf(self, other):
        """
        Return ``True`` if `self` is a specific case of `other`, or ``False``
        otherwise. If `other` is an :class:`Atom` object, then this is the same
        as the :meth:`equivalent()` method. If `other` is an
        :class:`GroupAtom` object, then the atom must match or be more
        specific than any of the combinations in the atom pattern.
        """
        if isinstance(other, Atom):
            return self.equivalent(other)
        elif isinstance(other, gr.GroupAtom):
            cython.declare(atom=gr.GroupAtom, a=AtomType, radical=cython.short, lp = cython.short, charge=cython.short)
            atom = other
            if self.atomType is None:
                return False
            for a in atom.atomType: 
                if self.atomType.isSpecificCaseOf(a): break
            else:
                return False
            if atom.radicalElectrons:
                for radical in atom.radicalElectrons:
                    if self.radicalElectrons == radical: break
                else:
                    return False
            if atom.lonePairs:
                for lp in atom.lonePairs:
                    if self.lonePairs == lp: break
                else:
                    return False
            if atom.charge:
                for charge in atom.charge:
                    if self.charge == charge: break
                else:
                    return False
            return True

    def copy(self):
        """
        Generate a deep copy of the current atom. Modifying the
        attributes of the copy will not affect the original.
        """
        cython.declare(a=Atom)
        #a = Atom(self.element, self.radicalElectrons, self.spinMultiplicity, self.charge, self.label)
        a = Atom.__new__(Atom)
        a.edges = {}
        a.resetConnectivityValues()
        a.element = self.element
        a.radicalElectrons = self.radicalElectrons
        a.charge = self.charge
        a.label = self.label
        a.atomType = self.atomType
        a.lonePairs = self.lonePairs
        a.coords = self.coords[:]
        return a

    def isHydrogen(self):
        """
        Return ``True`` if the atom represents a hydrogen atom or ``False`` if
        not.
        """
        return self.element.number == 1

    def isNonHydrogen(self):
        """
        Return ``True`` if the atom does not represent a hydrogen atom or
        ``False`` if not.
        """
        return self.element.number > 1

    def isCarbon(self):
        """
        Return ``True`` if the atom represents a carbon atom or ``False`` if
        not.
        """
        return self.element.number == 6
    
    def isNitrogen(self):
        """
        Return ``True`` if the atom represents a nitrogen atom or ``False`` if
        not.
        """
        return self.element.number == 7

    def isOxygen(self):
        """
        Return ``True`` if the atom represents an oxygen atom or ``False`` if
        not.
        """
        return self.element.number == 8

    def isSilicon(self):
        """
        Return ``True`` if the atom represents an silicon atom or ``False`` if
        not.
        """
        return self.element.number == 14

    def isSulfur(self):
        """
        Return ``True`` if the atom represents an sulfur atom or ``False`` if
        not.
        """
        return self.element.number == 16

    def incrementRadical(self):
        """
        Update the atom pattern as a result of applying a GAIN_RADICAL action,
        where `radical` specifies the number of radical electrons to add.
        """
        # Set the new radical electron count
        self.radicalElectrons += 1
        if self.radicalElectrons <= 0:
            raise gr.ActionError('Unable to update Atom due to GAIN_RADICAL action: Invalid radical electron set "{0}".'.format(self.radicalElectrons))

    def decrementRadical(self):
        """
        Update the atom pattern as a result of applying a LOSE_RADICAL action,
        where `radical` specifies the number of radical electrons to remove.
        """
        cython.declare(radicalElectrons=cython.short)
        # Set the new radical electron count
        radicalElectrons = self.radicalElectrons = self.radicalElectrons - 1
        if radicalElectrons  < 0:
            raise gr.ActionError('Unable to update Atom due to LOSE_RADICAL action: Invalid radical electron set "{0}".'.format(self.radicalElectrons))

    def setLonePairs(self,lonePairs):
        """
        Set the number of lone electron pairs.
        """
        # Set the number of electron pairs
        self.lonePairs = lonePairs
        if self.lonePairs < 0:
            raise gr.ActionError('Unable to update Atom due to setLonePairs : Invalid lone electron pairs set "{0}".'.format(self.setLonePairs))
        self.updateCharge()

    def incrementLonePairs(self):
        """
        Update the lone electron pairs pattern as a result of applying a GAIN_PAIR action.
        """
        # Set the new lone electron pairs count
        self.lonePairs += 1
        if self.lonePairs <= 0:
            raise gr.ActionError('Unable to update Atom due to GAIN_PAIR action: Invalid lone electron pairs set "{0}".'.format(self.lonePairs))
        self.updateCharge()

    def decrementLonePairs(self):
        """
        Update the lone electron pairs pattern as a result of applying a LOSE_PAIR action.
        """
        # Set the new lone electron pairs count
        self.lonePairs -= 1
        if self.lonePairs  < 0:
            raise gr.ActionError('Unable to update Atom due to LOSE_PAIR action: Invalid lone electron pairs set "{0}".'.format(self.lonePairs))
        self.updateCharge()

    def updateCharge(self):
        """
        Update self.charge, according to the valence, and the
        number and types of bonds, radicals, and lone pairs.
        """
        valence_electron = elements.PeriodicSystem.valence_electrons[self.symbol]
        order = self.getBondOrdersForAtom()
        self.charge = valence_electron - order - self.radicalElectrons - 2*self.lonePairs

    def applyAction(self, action):
        """
        Update the atom pattern as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        # Invalidate current atom type
        self.atomType = None
        act = action[0].upper()
        # Modify attributes if necessary
        if act in ['CHANGE_BOND', 'FORM_BOND', 'BREAK_BOND']:
            # Nothing else to do here
            pass
        elif act == 'GAIN_RADICAL':
            for i in range(action[2]): self.incrementRadical()
        elif act == 'LOSE_RADICAL':
            for i in range(abs(action[2])): self.decrementRadical()
        elif action[0].upper() == 'GAIN_PAIR':
            for i in range(action[2]): self.incrementLonePairs()
        elif action[0].upper() == 'LOSE_PAIR':
            for i in range(abs(action[2])): self.decrementLonePairs()
        else:
            raise gr.ActionError('Unable to update Atom: Invalid action {0}".'.format(action))

    def setSpinMultiplicity(self,spinMultiplicity):
        """
        Set the spin multiplicity.
        """
        raise NotImplementedError("I thought multiplicity was now a molecule attribute not atom?")
        # Set the spin multiplicity
        self.spinMultiplicity = spinMultiplicity
        if self.spinMultiplicity < 0:
            raise gr.ActionError('Unable to update Atom due to spin multiplicity : Invalid spin multiplicity set "{0}".'.format(self.spinMultiplicity))
        self.updateCharge()

    def getBondOrdersForAtom(self):
        """
        This helper function is to help calculate total bond orders for an
        input atom.

        Some special consideration for the order `B` bond. For atoms having 
        three `B` bonds, the order for each is 4/3.0, while for atoms having other
        than three `B` bonds, the order for  each is 3/2.0
        """
        num_B_bond = 0
        order = 0
        for _, bond in self.bonds.iteritems():
            if bond.isBenzene():
                num_B_bond += 1
            else:
                order += bond.order

        if num_B_bond == 3:
            order += num_B_bond * 4/3.0
        else:
            order += num_B_bond * 3/2.0

        return order

################################################################################

class Bond(Edge):
    """
    A chemical bond. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `order`             ``float``             The :ref:`bond type <bond-types>`
    =================== =================== ====================================

    """

    def __init__(self, atom1, atom2, order=1):
        Edge.__init__(self, atom1, atom2)
        if isinstance(order, str):
            self.setOrderStr(order)
        else:
            self.order = order

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return self.getOrderStr()

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return '<Bond "{0}">'.format(self.order)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Bond, (self.vertex1, self.vertex2, self.order))

    @property
    def atom1(self):
        return self.vertex1

    @property
    def atom2(self):
        return self.vertex2

    def equivalent(self, other):
        """
        Return ``True`` if `other` is indistinguishable from this bond, or
        ``False`` otherwise. `other` can be either a :class:`Bond` or a
        :class:`GroupBond` object.
        """
        cython.declare(bond=Bond, bp=gr.GroupBond)
        if isinstance(other, Bond):
            bond = other
            return (self.getOrderNum() == bond.getOrderNum())
        elif isinstance(other, gr.GroupBond):
            bp = other
            return (self.getOrderNum() in bp.getOrderNum())

    def isSpecificCaseOf(self, other):
        """
        Return ``True`` if `self` is a specific case of `other`, or ``False``
        otherwise. `other` can be either a :class:`Bond` or a
        :class:`GroupBond` object.
        """
        # There are no generic bond types, so isSpecificCaseOf is the same as equivalent
        return self.equivalent(other)

    def getOrderStr(self):
        """
        returns a string representing the bond order
        """
        if self.isSingle():
            return 'S'
        elif self.isBenzene():
            return 'B'
        elif self.isDouble():
            return 'D'
        elif self.isTriple():
            return 'T'
        else:
            raise ValueError("Bond order {} does not have string representation." +  \
            "".format(self.order))
        
    def setOrderStr(self, newOrder):
        """
        set the bond order using a valid bond-order character
        """
        if newOrder == 'S':
            self.order = 1
        elif newOrder == 'D':
            self.order = 2
        elif newOrder == 'T':
            self.order = 3
        elif newOrder == 'B':
            self.order = 1.5
        else:
            raise TypeError('Bond order {} is not hardcoded into this method'.format(newOrder))
            
            
    def getOrderNum(self):
        """
        returns the bond order as a number
        """
        
        return self.order
            
    def setOrderNum(self, newOrder):
        """
        change the bond order with a number
        """
        
        self.order = newOrder
        
    def copy(self):
        """
        Generate a deep copy of the current bond. Modifying the
        attributes of the copy will not affect the original.
        """
        #return Bond(self.vertex1, self.vertex2, self.order)
        cython.declare(b=Bond)
        b = Bond.__new__(Bond)
        b.vertex1 = self.vertex1
        b.vertex2 = self.vertex2
        b.order = self.order
        return b

    def isOrder(self, otherOrder):
        """
        Return ``True`` if the bond represents a single bond or ``False`` if
        not. This compares floats that takes into account floating point error
        
        NOTE: we can replace the absolute value relation with math.isclose when
        we swtich to python 3.5+
        """
        return abs(self.order - otherOrder) <= 1e-9

        
    def isSingle(self):
        """
        Return ``True`` if the bond represents a single bond or ``False`` if
        not.
        """
        return abs(self.order-1) <= 1e-9

    def isDouble(self):
        """
        Return ``True`` if the bond represents a double bond or ``False`` if
        not.
        """
        return abs(self.order-2) <= 1e-9

    def isTriple(self):
        """
        Return ``True`` if the bond represents a triple bond or ``False`` if
        not.
        """
        return abs(self.order-3) <= 1e-9

    def isBenzene(self):
        """
        Return ``True`` if the bond represents a benzene bond or ``False`` if
        not.
        """
        return abs(self.order-1.5) <= 1e-9

    def incrementOrder(self):
        """
        Update the bond as a result of applying a CHANGE_BOND action to
        increase the order by one.
        """
        if self.order <=2: 
            self.order += 1
        else:
            raise gr.ActionError('Unable to increment Bond due to CHANGE_BOND action: '+\
            'Bond order "{0}" is greater than 2.'.format(self.order))

    def decrementOrder(self):
        """
        Update the bond as a result of applying a CHANGE_BOND action to
        decrease the order by one.
        """
        if self.order >=1: 
            self.order -= 1
        else:
            raise gr.ActionError('Unable to decrease Bond due to CHANGE_BOND action: '+\
            'bond order "{0}" is less than 1.'.format(self.order))

    def __changeBond(self, order):
        """
        Update the bond as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order, and can be any real number.
        """
        self.order += order
        if self.order < 0 or self.order >3:
            raise gr.ActionError('Unable to update Bond due to CHANGE_BOND action: Invalid resulting order "{0}".'.format(self.order))

    def applyAction(self, action):
        """
        Update the bond as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        if action[0].upper() == 'CHANGE_BOND':
            if isinstance(action[2],str):
                self.setOrderStr(action[2])
            else:
                try: # try to see if addable
                   self.__changeBond(action[2])
                except TypeError:
                    raise gr.ActionError('Unable to update Bond due to CHANGE_BOND action: Invalid order "{0}".'.format(action[2]))
        else:
            raise gr.ActionError('Unable to update GroupBond: Invalid action {0}.'.format(action))
        

#################################################################################
    

class Molecule(Graph):
    """
    A representation of a molecular structure using a graph data type, extending
    the :class:`Graph` class. The `atoms` and `bonds` attributes are aliases
    for the `vertices` and `edges` attributes. Other attributes are:

    ======================= =========== ========================================
    Attribute               Type        Description
    ======================= =========== ========================================
    `symmetryNumber`        ``int``     The (estimated) external + internal symmetry number of the molecule
    `multiplicity`          ``int``     The multiplicity of this species, multiplicity = 2*total_spin+1
    ======================= =========== ========================================

    A new molecule object can be easily instantiated by passing the `SMILES` or
    `InChI` string representing the molecular structure.
    """

    def __init__(self, atoms=None, symmetry=-1, multiplicity=-187, props=None, SMILES=''):
        Graph.__init__(self, atoms)
        self.symmetryNumber = symmetry
        self.multiplicity = multiplicity
        self._fingerprint = None
        self.InChI = ''
        if SMILES != '': self.fromSMILES(SMILES)
        self.props = props or {}
        if multiplicity != -187:  # it was set explicitly, so re-set it (fromSMILES etc may have changed it)
            self.multiplicity = multiplicity
    
    
    def __hash__(self):
        return hash((self.getFingerprint()))
            
    def __richcmp__(x, y, op):
        if op == 2:#Py_EQ
            return x.is_equal(y)
        if op == 3:#Py_NE
            return not x.is_equal(y)
        else:
            raise NotImplementedError("Can only check equality of molecules, not > or <")
    
    def is_equal(self,other):
        """Method to test equality of two Molecule objects."""
        if not isinstance(other, Molecule): return False #different type
        elif self is other: return True #same reference in memory
        elif self.getFingerprint() != other.getFingerprint(): return False
        else:
            return self.isIsomorphic(other)   

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return '<Molecule "{0}">'.format(self.toSMILES())

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        cython.declare(multiplicity=cython.int)
        multiplicity = self.multiplicity
        if multiplicity != self.getRadicalCount() + 1:
            return 'Molecule(SMILES="{0}", multiplicity={1:d})'.format(self.toSMILES(), multiplicity)
        return 'Molecule(SMILES="{0}")'.format(self.toSMILES())

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Molecule, (self.vertices, self.symmetryNumber, self.multiplicity, self.props))

    def __getAtoms(self): return self.vertices
    def __setAtoms(self, atoms): self.vertices = atoms
    atoms = property(__getAtoms, __setAtoms)

    def addAtom(self, atom):
        """
        Add an `atom` to the graph. The atom is initialized with no bonds.
        """
        self._fingerprint = None
        return self.addVertex(atom)
    
    def addBond(self, bond):
        """
        Add a `bond` to the graph as an edge connecting the two atoms `atom1`
        and `atom2`.
        """
        self._fingerprint = None
        return self.addEdge(bond)

    def getBonds(self, atom):
        """
        Return a list of the bonds involving the specified `atom`.
        """
        return self.getEdges(atom)

    def getBond(self, atom1, atom2):
        """
        Returns the bond connecting atoms `atom1` and `atom2`.
        """
        return self.getEdge(atom1, atom2)

    def hasAtom(self, atom):
        """
        Returns ``True`` if `atom` is an atom in the graph, or ``False`` if
        not.
        """
        return self.hasVertex(atom)

    def hasBond(self, atom1, atom2):
        """
        Returns ``True`` if atoms `atom1` and `atom2` are connected
        by an bond, or ``False`` if not.
        """
        return self.hasEdge(atom1, atom2)

    def removeAtom(self, atom):
        """
        Remove `atom` and all bonds associated with it from the graph. Does
        not remove atoms that no longer have any bonds as a result of this
        removal.
        """
        self._fingerprint = None
        return self.removeVertex(atom)

    def removeBond(self, bond):
        """
        Remove the bond between atoms `atom1` and `atom2` from the graph.
        Does not remove atoms that no longer have any bonds as a result of
        this removal.
        """
        self._fingerprint = None
        return self.removeEdge(bond)

    def sortAtoms(self):
        """
        Sort the atoms in the graph. This can make certain operations, e.g.
        the isomorphism functions, much more efficient.
        """
        return self.sortVertices()

    def update(self):
        """
        Update connectivity values, atom types of atoms.
        Update multiplicity, and sort atoms using the new
        connectivity values.
        """
        self.updateAtomTypes()
        self.updateMultiplicity()
        self.sortVertices()

        for atom in self.atoms:
            atom.updateCharge()

    def getFormula(self):
        """
        Return the molecular formula for the molecule.
        """
        cython.declare(atom=Atom, symbol=str, elements=dict, keys=list, formula=str)
        cython.declare(hasCarbon=cython.bint, hasHydrogen=cython.bint)
        
        # Count the number of each element in the molecule
        elements = {}
        for atom in self.vertices:
            symbol = atom.element.symbol
            elements[symbol] = elements.get(symbol, 0) + 1
        
        # Use the Hill system to generate the formula
        formula = ''
        
        # Carbon and hydrogen always come first if carbon is present
        if 'C' in elements.keys():
            count = elements['C']
            formula += 'C{0:d}'.format(count) if count > 1 else 'C'
            del elements['C']
            if 'H' in elements.keys():
                count = elements['H']
                formula += 'H{0:d}'.format(count) if count > 1 else 'H'
                del elements['H']

        # Other atoms are in alphabetical order
        # (This includes hydrogen if carbon is not present)
        keys = elements.keys()
        keys.sort()
        for key in keys:
            count = elements[key]
            formula += '{0}{1:d}'.format(key, count) if count > 1 else key
        
        return formula

    def getMolecularWeight(self):
        """
        Return the molecular weight of the molecule in kg/mol.
        """
        cython.declare(atom=Atom, mass=cython.double)
        mass = 0
        for atom in self.vertices:
            mass += atom.element.mass
        return mass
    
    def getRadicalCount(self):
        """
        Return the total number of radical electrons on all atoms in the
        molecule. In this function, monoradical atoms count as one, biradicals
        count as two, etc.
        """
        cython.declare(atom=Atom, radicals=cython.short)
        radicals = 0
        for atom in self.vertices:
            radicals += atom.radicalElectrons
        return radicals

    def getNumAtoms(self, element = None):
        """
        Return the number of atoms in molecule.  If element is given, ie. "H" or "C",
        the number of atoms of that element is returned.
        """
        cython.declare(numAtoms=cython.int, atom=Atom)
        if element == None:
            return len(self.vertices)
        else:
            numAtoms = 0
            for atom in self.vertices:
                if atom.element.symbol == element:
                    numAtoms += 1
            return numAtoms

    def copy(self, deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """
        other = cython.declare(Molecule)
        g = Graph.copy(self, deep)
        other = Molecule(g.vertices)
        other.multiplicity = self.multiplicity
        return other

    def merge(self, other):
        """
        Merge two molecules so as to store them in a single :class:`Molecule`
        object. The merged :class:`Molecule` object is returned.
        """
        g = Graph.merge(self, other)
        molecule = Molecule(atoms=g.vertices)
        return molecule

    def split(self):
        """
        Convert a single :class:`Molecule` object containing two or more
        unconnected molecules into separate class:`Molecule` objects.
        """
        graphs = Graph.split(self)
        molecules = []
        for g in graphs:
            molecule = Molecule(atoms=g.vertices)
            molecules.append(molecule)
        return molecules

    def deleteHydrogens(self):
        """
        Irreversibly delete all non-labeled hydrogens without updating
        connectivity values. If there's nothing but hydrogens, it does nothing.
        It destroys information; be careful with it.
        """
        cython.declare(atom=Atom, hydrogens=list)
        # Check that the structure contains at least one heavy atom
        for atom in self.vertices:
            if not atom.isHydrogen():
                break
        else:
            # No heavy atoms, so leave explicit
            return
        hydrogens = []
        for atom in self.vertices:
            if atom.isHydrogen() and atom.label == '':
                hydrogens.append(atom)
        # Remove the hydrogen atoms from the structure
        for atom in hydrogens:
            self.removeAtom(atom)

    def connectTheDots(self):
        """
        Delete all bonds, and set them again based on the Atoms' coords.
        Does not detect bond type.
        """
        cython.declare(criticalDistance=float, i=int, atom1=Atom, atom2=Atom,
                       bond=Bond, atoms=list, zBoundary=float)
                       # groupBond=GroupBond, 
        self._fingerprint = None
        
        atoms = self.vertices
        
        # Ensure there are coordinates to work with
        for atom in atoms:
            assert len(atom.coords) != 0
        
        # If there are any bonds, remove them
        for atom1 in atoms:
            for bond in self.getBonds(atom1):
                self.removeEdge(bond)
        
        # Sort atoms by distance on the z-axis
        sortedAtoms = sorted(atoms, key=lambda x: x.coords[2])
        
        for i, atom1 in enumerate(sortedAtoms):
            for atom2 in sortedAtoms[i+1:]:
                # Set upper limit for bond distance
                criticalDistance = (atom1.element.covRadius + atom2.element.covRadius + 0.45)**2
                
                # First atom that is more than 4.0 Anstroms away in the z-axis, break the loop
                # Atoms are sorted along the z-axis, so all following atoms should be even further
                zBoundary = (atom1.coords[2] - atom2.coords[2])**2
                if zBoundary > 16.0:
                    break
                
                distanceSquared = sum((atom1.coords - atom2.coords)**2)
                
                if distanceSquared > criticalDistance or distanceSquared < 0.40:
                    continue
                else:
                    # groupBond = GroupBond(atom1, atom2, [1,2,3,1.5])
                    bond = Bond(atom1, atom2, 1)
                    self.addBond(bond)
        self.updateAtomTypes()
        
    def updateAtomTypes(self, logSpecies=True):
        """
        Iterate through the atoms in the structure, checking their atom types
        to ensure they are correct (i.e. accurately describe their local bond
        environment) and complete (i.e. are as detailed as possible).
        """
        for atom in self.vertices:
            try:
                atom.atomType = getAtomType(atom, atom.edges)
            except AtomTypeError:
                if logSpecies:
                    logging.error("Could not update atomtypes for {0}.\n{1}".format(self, self.toAdjacencyList()))
                raise
            
    def updateMultiplicity(self):
        """
        Update the multiplicity of a newly formed molecule.
        """
        # Assume this is always true
        # There are cases where 2 radicalElectrons is a singlet, but
        # the triplet is often more stable, 
        self.multiplicity = self.getRadicalCount() + 1

    def clearLabeledAtoms(self):
        """
        Remove the labels from all atoms in the molecule.
        """
        for atom in self.vertices:
            atom.label = ''

    def containsLabeledAtom(self, label):
        """
        Return :data:`True` if the molecule contains an atom with the label
        `label` and :data:`False` otherwise.
        """
        for atom in self.vertices:
            if atom.label == label: return True
        return False

    def getLabeledAtom(self, label):
        """
        Return the atoms in the molecule that are labeled.
        """
        for atom in self.vertices:
            if atom.label == label: return atom
        raise ValueError('No atom in the molecule has the label "{0}".'.format(label))

    def getLabeledAtoms(self):
        """
        Return the labeled atoms as a ``dict`` with the keys being the labels
        and the values the atoms themselves. If two or more atoms have the
        same label, the value is converted to a list of these atoms.
        """
        labeled = {}
        for atom in self.vertices:
            if atom.label != '':
                if atom.label in labeled:
                    if isinstance(labeled[atom.label],list):
                        labeled[atom.label].append(atom)
                    else:
                        labeled[atom.label] = [labeled[atom.label]]
                        labeled[atom.label].append(atom)
                else:
                    labeled[atom.label] = atom
        return labeled

    def getFingerprint(self):
        """
        Return a string containing the "fingerprint" used to accelerate graph
        isomorphism comparisons with other molecules. The fingerprint is a
        short string containing a summary of selected information about the 
        molecule. Two fingerprint strings matching is a necessary (but not
        sufficient) condition for the associated molecules to be isomorphic.
        """
        if self._fingerprint is None:
            self._fingerprint = self.getFormula()
        return self._fingerprint
    
    def isIsomorphic(self, other, initialMap=None):
        """
        Returns :data:`True` if two graphs are isomorphic and :data:`False`
        otherwise. The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Molecule` object, or a :class:`TypeError` is raised.
        Also ensures multiplicities are also equal.
        """
        # It only makes sense to compare a Molecule to a Molecule for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Molecule):
            raise TypeError('Got a {0} object for parameter "other", when a Molecule object is required.'.format(other.__class__))
        # Do the quick isomorphism comparison using the fingerprint
        # Two fingerprint strings matching is a necessary (but not
        # sufficient!) condition for the associated molecules to be isomorphic
        if self.getFingerprint() != other.getFingerprint():
            return False
        # check multiplicity
        if self.multiplicity != other.multiplicity:
            return False
        # Do the full isomorphism comparison
        result = Graph.isIsomorphic(self, other, initialMap)
        return result

    def findIsomorphism(self, other, initialMap=None):
        """
        Returns :data:`True` if `other` is isomorphic and :data:`False`
        otherwise, and the matching mapping. The `initialMap` attribute can be
        used to specify a required mapping from `self` to `other` (i.e. the
        atoms of `self` are the keys, while the atoms of `other` are the
        values). The returned mapping also uses the atoms of `self` for the keys
        and the atoms of `other` for the values. The `other` parameter must
        be a :class:`Molecule` object, or a :class:`TypeError` is raised.
        """
        # It only makes sense to compare a Molecule to a Molecule for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Molecule):
            raise TypeError('Got a {0} object for parameter "other", when a Molecule object is required.'.format(other.__class__))
        # Do the quick isomorphism comparison using the fingerprint
        # Two fingerprint strings matching is a necessary (but not
        # sufficient!) condition for the associated molecules to be isomorphic
        if self.getFingerprint() != other.getFingerprint():
            return []
        # check multiplicity
        if self.multiplicity != other.multiplicity:
            return []
            
        # Do the isomorphism comparison
        result = Graph.findIsomorphism(self, other, initialMap)
        return result

    def isSubgraphIsomorphic(self, other, initialMap=None):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise. The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """
        cython.declare(group=gr.Group, atom=Atom)
        cython.declare(carbonCount=cython.short, nitrogenCount=cython.short, oxygenCount=cython.short, sulfurCount=cython.short, radicalCount=cython.short)
        
        # It only makes sense to compare a Molecule to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, gr.Group):
            raise TypeError('Got a {0} object for parameter "other", when a Molecule object is required.'.format(other.__class__))
        group = other
        
        # Count the number of carbons, oxygens, and radicals in the molecule
        carbonCount = 0; nitrogenCount = 0; oxygenCount = 0; sulfurCount = 0; radicalCount = 0
        for atom in self.vertices:
            if atom.element.symbol == 'C':
                carbonCount += 1
            elif atom.element.symbol == 'N':
                nitrogenCount += 1
            elif atom.element.symbol == 'O':
                oxygenCount += 1
            elif atom.element.symbol == 'S':
                sulfurCount += 1
            radicalCount += atom.radicalElectrons
        
        
        if group.multiplicity:
            if self.multiplicity not in group.multiplicity: return False
        # If the molecule has fewer of any of these things than the functional
        # group does, then we know the subgraph isomorphism fails without
        # needing to perform the full isomorphism check
        if (radicalCount < group.radicalCount or
            carbonCount < group.carbonCount or
            nitrogenCount < group.nitrogenCount or
            oxygenCount < group.oxygenCount or
            sulfurCount < group.sulfurCount):
            return False

        # Do the isomorphism comparison
        result = Graph.isSubgraphIsomorphic(self, other, initialMap)
        return result

    def findSubgraphIsomorphisms(self, other, initialMap=None):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise. Also returns the lists all of valid mappings. The
        `initialMap` attribute can be used to specify a required mapping from
        `self` to `other` (i.e. the atoms of `self` are the keys, while the
        atoms of `other` are the values). The returned mappings also use the
        atoms of `self` for the keys and the atoms of `other` for the values.
        The `other` parameter must be a :class:`Group` object, or a
        :class:`TypeError` is raised.
        """
        cython.declare(group=gr.Group, atom=Atom)
        cython.declare(carbonCount=cython.short, nitrogenCount=cython.short, oxygenCount=cython.short, sulfurCount=cython.short, radicalCount=cython.short)

        # It only makes sense to compare a Molecule to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, gr.Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        group = other
                # Count the number of carbons, oxygens, and radicals in the molecule
        carbonCount = 0; nitrogenCount = 0; oxygenCount = 0; sulfurCount = 0; radicalCount = 0
        for atom in self.vertices:
            if atom.element.symbol == 'C':
                carbonCount += 1
            elif atom.element.symbol == 'N':
                nitrogenCount += 1
            elif atom.element.symbol == 'O':
                oxygenCount += 1
            elif atom.element.symbol == 'S':
                sulfurCount += 1
            radicalCount += atom.radicalElectrons
        
        
        if group.multiplicity:
            if self.multiplicity not in group.multiplicity: return []
        # If the molecule has fewer of any of these things than the functional
        # group does, then we know the subgraph isomorphism fails without
        # needing to perform the full isomorphism check
        if (radicalCount < group.radicalCount or
            carbonCount < group.carbonCount or
            nitrogenCount < group.nitrogenCount or
            oxygenCount < group.oxygenCount or
            sulfurCount < group.sulfurCount):
            return []
        # Do the isomorphism comparison
        result = Graph.findSubgraphIsomorphisms(self, other, initialMap)
        return result

    def isAtomInCycle(self, atom):
        """
        Return :data:`True` if `atom` is in one or more cycles in the structure,
        and :data:`False` if not.
        """
        return self.isVertexInCycle(atom)

    def isBondInCycle(self, bond):
        """
        Return :data:`True` if the bond between atoms `atom1` and `atom2`
        is in one or more cycles in the graph, or :data:`False` if not.
        """
        return self.isEdgeInCycle(bond)

    def draw(self, path):
        """
        Generate a pictorial representation of the chemical graph using the
        :mod:`draw` module. Use `path` to specify the file to save
        the generated image to; the image type is automatically determined by
        extension. Valid extensions are ``.png``, ``.svg``, ``.pdf``, and
        ``.ps``; of these, the first is a raster format and the remainder are
        vector formats.
        """
        from .draw import MoleculeDrawer
        format = os.path.splitext(path)[-1][1:].lower()
        MoleculeDrawer().draw(self, format, target=path)
    
    def _repr_png_(self):
        """
        Return a png picture of the molecule, useful for ipython-qtconsole.
        """
        from .draw import MoleculeDrawer
        tempFileName = 'temp_molecule.png'
        MoleculeDrawer().draw(self, 'png', tempFileName)
        png = open(tempFileName,'rb').read()
        os.unlink(tempFileName)
        return png

    def fromInChI(self, inchistr, backend='try-all'):
        """
        Convert an InChI string `inchistr` to a molecular structure.
        """
        parser.fromInChI(self, inchistr, backend)
        return self

    def fromAugmentedInChI(self, aug_inchi):
        """
        Convert an Augmented InChI string `aug_inchi` to a molecular structure.
        """
        parser.fromAugmentedInChI(self, aug_inchi)
        return self

    def fromSMILES(self, smilesstr, backend='try-all'):
        """
        Convert a SMILES string `smilesstr` to a molecular structure.
        """
        parser.fromSMILES(self, smilesstr, backend)
        return self
        
    def fromSMARTS(self, smartsstr):
        """
        Convert a SMARTS string `smartsstr` to a molecular structure. Uses
        `RDKit <http://rdkit.org/>`_ to perform the conversion.
        This Kekulizes everything, removing all aromatic atom types.
        """
        parser.fromSMARTS(self, smartsstr)
        return self

    def fromAdjacencyList(self, adjlist, saturateH=False):
        """
        Convert a string adjacency list `adjlist` to a molecular structure.
        Skips the first line (assuming it's a label) unless `withLabel` is
        ``False``.
        """
        from .adjlist import fromAdjacencyList
        
        self.vertices, self.multiplicity = fromAdjacencyList(adjlist, group=False, saturateH=saturateH)
        self.updateAtomTypes()
        
        # Check if multiplicity is possible
        n_rad = self.getRadicalCount() 
        multiplicity = self.multiplicity
        if not (n_rad + 1 == multiplicity or n_rad - 1 == multiplicity or n_rad - 3 == multiplicity or n_rad - 5 == multiplicity):
            raise ValueError('Impossible multiplicity for molecule\n{0}\n multiplicity = {1} and number of unpaired electrons = {2}'.format(self.toAdjacencyList(),multiplicity,n_rad))
        if self.getNetCharge() != 0:
            raise ValueError('Non-neutral molecule encountered. Currently, RMG does not support ion chemistry.\n {0}'.format(adjlist))
        return self
        
    def fromXYZ(self, atomicNums, coordinates):
        """
        Create an RMG molecule from a list of coordinates and a corresponding
        list of atomic numbers. These are typically received from CCLib and the
        molecule is sent to `ConnectTheDots` so will only contain single bonds.
        """
        
        _rdkit_periodic_table = elements.GetPeriodicTable()
        
        for i, atNum in enumerate(atomicNums):
            atom = Atom(_rdkit_periodic_table.GetElementSymbol(int(atNum)))
            atom.coords = coordinates[i]
            self.addAtom(atom)
        return self.connectTheDots()
    
    def toSingleBonds(self):
        """
        Returns a copy of the current molecule, consisting of only single bonds.
        
        This is useful for isomorphism comparison against something that was made
        via fromXYZ, which does not attempt to perceive bond orders
        """
        cython.declare(atom1=Atom, atom2=Atom, bond=Bond, newMol=Molecule, atoms=list, mapping=dict)
    
        newMol = Molecule()
        atoms = self.atoms
        mapping = {}
        for atom1 in atoms:
            atom2 = newMol.addAtom(Atom(atom1.element))
            mapping[atom1] = atom2
    
        for atom1 in atoms:
            for atom2 in atom1.bonds:
                bond = Bond(mapping[atom1], mapping[atom2], 1)
                newMol.addBond(bond)
        newMol.updateAtomTypes()
        return newMol

    def toInChI(self):
        """
        Convert a molecular structure to an InChI string. Uses
        `RDKit <http://rdkit.org/>`_ to perform the conversion.
        Perceives aromaticity.
        
        or
        
        Convert a molecular structure to an InChI string. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        return generator.toInChI(self)            
        
    def toAugmentedInChI(self):
        """
        Adds an extra layer to the InChI denoting the multiplicity
        of the molecule.
        
        Separate layer with a forward slash character.
        """
        return generator.toAugmentedInChI(self)
        
    
    def toInChIKey(self):
        """
        Convert a molecular structure to an InChI Key string. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        
        or 
        
        Convert a molecular structure to an InChI Key string. Uses
        `RDKit <http://rdkit.org/>`_ to perform the conversion.
        
        Removes check-sum dash (-) and character so that only 
        the 14 + 9 characters remain.
        """
        return generator.toInChIKey(self)
    
    def toAugmentedInChIKey(self):
        """
        Adds an extra layer to the InChIKey denoting the multiplicity
        of the molecule.

        Simply append the multiplicity string, do not separate by a
        character like forward slash.
        """
        return generator.toAugmentedInChIKey(self)
    

    def toSMARTS(self):
        """
        Convert a molecular structure to an SMARTS string. Uses
        `RDKit <http://rdkit.org/>`_ to perform the conversion.
        Perceives aromaticity and removes Hydrogen atoms.
        """
        return generator.toSMARTS(self)
    
    
    def toSMILES(self):
        """
        Convert a molecular structure to an SMILES string. 
        
        If there is a Nitrogen atom present it uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion,
        and the SMILES may or may not be canonical.
        
        Otherwise, it uses `RDKit <http://rdkit.org/>`_ to perform the 
        conversion, so it will be canonical SMILES.
        While converting to an RDMolecule it will perceive aromaticity
        and removes Hydrogen atoms.
        """
        
        return generator.toSMILES(self)

    def toRDKitMol(self, *args, **kwargs):
        """
        Convert a molecular structure to a RDKit rdmol object.
        """
        return generator.toRDKitMol(self, *args, **kwargs)

    def toAdjacencyList(self, label='', removeH=False, removeLonePairs=False, oldStyle=False):
        """
        Convert the molecular structure to a string adjacency list.
        """
        from .adjlist import toAdjacencyList
        result = toAdjacencyList(self.vertices, self.multiplicity,  label=label, group=False, removeH=removeH, removeLonePairs=removeLonePairs, oldStyle=oldStyle)
        return result

    def isLinear(self):
        """
        Return :data:`True` if the structure is linear and :data:`False`
        otherwise.
        """

        atomCount = len(self.vertices)

        # Monatomic molecules are definitely nonlinear
        if atomCount == 1:
            return False
        # Diatomic molecules are definitely linear
        elif atomCount == 2:
            return True
        # Cyclic molecules are definitely nonlinear
        elif self.isCyclic():
            return False

        # True if all bonds are double bonds (e.g. O=C=O)
        allDoubleBonds = True
        for atom1 in self.vertices:
            for bond in atom1.edges.values():
                if not bond.isDouble(): allDoubleBonds = False
        if allDoubleBonds: return True

        # True if alternating single-triple bonds (e.g. H-C#C-H)
        # This test requires explicit hydrogen atoms
        for atom in self.vertices:
            bonds = atom.edges.values()
            if len(bonds)==1:
                continue # ok, next atom
            if len(bonds)>2:
                break # fail!
            if bonds[0].isSingle() and bonds[1].isTriple():
                continue # ok, next atom
            if bonds[1].isSingle() and bonds[0].isTriple():
                continue # ok, next atom
            break # fail if we haven't continued
        else:
            # didn't fail
            return True
        
        # not returned yet? must be nonlinear
        return False
    
    def isAromatic(self):
        """ 
        Returns ``True`` if the molecule is aromatic, or ``False`` if not.  
        Iterates over the SSSR's and searches for rings that consist solely of Cb 
        atoms.  Assumes that aromatic rings always consist of 6 atoms. 
        In cases of naphthalene, where a 6 + 4 aromatic system exists,
        there will be at least one 6 membered aromatic ring so this algorithm
        will not fail for fused aromatic rings.
        """
        cython.declare(SSSR=list, vertices=list, polycyclicVertices=list)
        SSSR = self.getSmallestSetOfSmallestRings()
        if SSSR:
            for cycle in SSSR:
                if len(cycle) == 6:
                    for atom in cycle:
                        #print atom.atomType.label
                        if atom.atomType.label == 'Cb' or atom.atomType.label == 'Cbf':
                            continue                        
                        # Go onto next cycle if a non Cb atomtype was discovered in this cycle
                        break 
                    else:
                        # Molecule is aromatic when all 6 atoms are type 'Cb'
                        return True    
        return False

    def countInternalRotors(self):
        """
        Determine the number of internal rotors in the structure. Any single
        bond not in a cycle and between two atoms that also have other bonds
        are considered to be internal rotors.
        """
        count = 0
        for atom1 in self.vertices:
            for atom2, bond in atom1.edges.items():
                if self.vertices.index(atom1) < self.vertices.index(atom2) and bond.isSingle() and not self.isBondInCycle(bond):
                    if len(atom1.edges) > 1 and len(atom2.edges) > 1:
                        count += 1
        return count

    def calculateCp0(self):
        """
        Return the value of the heat capacity at zero temperature in J/mol*K.
        """
        if len(self.atoms) == 1:
            return 2.5 * constants.R
        else:
            return (3.5 if self.isLinear() else 4.0) * constants.R

    def calculateCpInf(self):
        """
        Return the value of the heat capacity at infinite temperature in J/mol*K.
        """
        cython.declare(Natoms=cython.int, Nvib=cython.int, Nrotors=cython.int)
        
        if len(self.vertices) == 1:
            return self.calculateCp0()
        else:
            
            Natoms = len(self.vertices)
            Nvib = 3 * Natoms - (5 if self.isLinear() else 6)
            Nrotors = self.countInternalRotors()
            Nvib -= Nrotors
            
            return self.calculateCp0() + (Nvib + 0.5 * Nrotors) * constants.R

    def getSymmetryNumber(self):
        """
        Returns the symmetry number of Molecule.
        First checks whether the value is stored as an attribute of Molecule.
        If not, it calls the calculateSymmetryNumber method. 
        """
        if self.symmetryNumber == -1:
            self.calculateSymmetryNumber()
        return self.symmetryNumber
        
        
    def calculateSymmetryNumber(self):
        """
        Return the symmetry number for the structure. The symmetry number
        includes both external and internal modes.
        """
        from rmgpy.molecule.symmetry import calculateSymmetryNumber
        self.symmetryNumber = calculateSymmetryNumber(self)
        return self.symmetryNumber
    
    def isRadical(self):
        """
        Return ``True`` if the molecule contains at least one radical electron,
        or ``False`` otherwise.
        """
        cython.declare(atom=Atom)
        for atom in self.vertices:
            if atom.radicalElectrons > 0:
                return True
        return False

    def isArylRadical(self, ASSSR=None):
        """
        Return ``True`` if the molecule only contains aryl radicals,
        ie. radical on an aromatic ring, or ``False`` otherwise.
        """
        cython.declare(atom=Atom, radList=list)
        if ASSSR is None:
            ASSSR = self.getAromaticSSSR()[0]

        total = self.getRadicalCount()
        aromaticAtoms = set([atom for atom in itertools.chain.from_iterable(ASSSR)])
        aryl = sum([atom.radicalElectrons for atom in aromaticAtoms])

        return total == aryl

    def generateResonanceIsomers(self):
        return resonance.generateResonanceStructures(self)

    def getURL(self):
        """
        Get a URL to the molecule's info page on the RMG website.
        """
        # eg. http://dev.rmg.mit.edu/database/kinetics/reaction/reactant1=1%20C%200%20%7B2,S%7D;2%20O%200%20%7B1,S%7D;__reactant2=1%20C%202T;__product1=1%20C%201;__product2=1%20C%200%20%7B2,S%7D;2%20O%201%20%7B1,S%7D;

        base_url = "http://rmg.mit.edu/database/molecule/"
        adjlist = self.toAdjacencyList(removeH=False)
        url = base_url + urllib.quote(adjlist)
        return url.strip('_')
                    
    def getRadicalAtoms(self):
        """
        Return the atoms in the molecule that have unpaired electrons.
        """
        radicalAtomsList = []
        for atom in self.vertices:
            if atom.radicalElectrons > 0:
                radicalAtomsList.append(atom)
        return radicalAtomsList
    
    def updateLonePairs(self):
        """
        Iterate through the atoms in the structure and calculate the
        number of lone electron pairs, assuming a neutral molecule.
        """
        cython.declare(atom1=Atom, atom2=Atom, bond12=Bond, order=float)
        for atom1 in self.vertices:
            if not atom1.isHydrogen():
                order = atom1.getBondOrdersForAtom()
                atom1.lonePairs = (elements.PeriodicSystem.valence_electrons[atom1.symbol] - atom1.radicalElectrons - atom1.charge - int(order)) / 2.0
                if atom1.lonePairs % 1 > 0 or atom1.lonePairs > 4:
                    logging.error("Unable to determine the number of lone pairs for element {0} in {1}".format(atom1,self))
            else:
                atom1.lonePairs = 0
                
    def getNetCharge(self):
        """
        Iterate through the atoms in the structure and calculate the net charge
        on the overall molecule.
        """
        charge = 0
        for atom in self.vertices:
            charge += atom.charge
        return charge

    def saturate(self):
        """
        Saturate the molecule by replacing all radicals with bonds to hydrogen atoms.  Changes self molecule object.  
        """
        cython.declare(added=dict, atom=Atom, i=int, H=Atom, bond=Bond)
        added = {}
        for atom in self.atoms:
            for i in range(atom.radicalElectrons):
                H = Atom('H', radicalElectrons=0, lonePairs=0, charge=0)
                bond = Bond(atom, H, 1)
                self.addAtom(H)
                self.addBond(bond)
                if atom not in added:
                    added[atom] = []
                added[atom].append([H, bond])
                atom.decrementRadical()
      
        # Update the atom types of the saturated structure (not sure why
        # this is necessary, because saturating with H shouldn't be
        # changing atom types, but it doesn't hurt anything and is not
        # very expensive, so will do it anyway)
        self.sortVertices()
        self.updateAtomTypes()
        self.updateLonePairs()
        self.multiplicity = 1

        return added

    def toGroup(self):
        """
        This method converts a list of atoms in a Molecule to a Group object.
        """
        
        # Create GroupAtom object for each atom in the molecule
        groupAtoms = OrderedDict()# preserver order of atoms in original container
        for atom in self.atoms:
            groupAtoms[atom] = gr.GroupAtom(atomType=[atom.atomType],
                                         radicalElectrons=[atom.radicalElectrons],
                                         charge=[atom.charge],
                                         lonePairs=[atom.lonePairs]
                                         )
                    
        group = gr.Group(atoms=groupAtoms.values(), multiplicity=[self.multiplicity])
        
        # Create GroupBond for each bond between atoms in the molecule
        for atom in self.atoms:
            for bondedAtom, bond in atom.edges.iteritems():
                group.addBond(gr.GroupBond(groupAtoms[atom],groupAtoms[bondedAtom], order=[bond.order]))
            
        group.update()
        
        return group

    def getAromaticSSSR(self, SSSR=None):
        """
        Returns the smallest set of smallest aromatic rings as a list of atoms and a list of bonds

        Identifies rings using `Graph.getSmallestSetOfSmallestRings()`, then uses RDKit to perceive aromaticity.
        RDKit uses an atom-based pi-electron counting algorithm to check aromaticity based on Huckel's Rule.
        Therefore, this method identifies "true" aromaticity, rather than simply the RMG bond type.

        The method currently restricts aromaticity to six-membered carbon-only rings. This is a limitation imposed
        by RMG, and not by RDKit.
        """
        cython.declare(rdAtomIndices=dict, aromaticRings=list, aromaticBonds=list)
        cython.declare(rings=list, ring0=list, i=cython.int, atom1=Atom, atom2=Atom)

        from rdkit.Chem.rdchem import BondType

        AROMATIC = BondType.AROMATIC

        if SSSR is None:
            SSSR = self.getSmallestSetOfSmallestRings()

        rings = [ring0 for ring0 in SSSR if len(ring0) == 6]
        if not rings:
            return [], []

        try:
            rdkitmol, rdAtomIndices = generator.toRDKitMol(self, removeHs=False, returnMapping=True)
        except ValueError:
            return [], []

        aromaticRings = []
        aromaticBonds = []
        for ring0 in rings:
            aromaticBondsInRing = []
            # Figure out which atoms and bonds are aromatic and reassign appropriately:
            for i, atom1 in enumerate(ring0):
                if not atom1.isCarbon():
                    # all atoms in the ring must be carbon in RMG for our definition of aromatic
                    break
                for atom2 in ring0[i + 1:]:
                    if self.hasBond(atom1, atom2):
                        if rdkitmol.GetBondBetweenAtoms(rdAtomIndices[atom1],
                                                        rdAtomIndices[atom2]).GetBondType() is AROMATIC:
                            aromaticBondsInRing.append(self.getBond(atom1, atom2))
            else:  # didn't break so all atoms are carbon
                if len(aromaticBondsInRing) == 6:
                    aromaticRings.append(ring0)
                    aromaticBonds.append(aromaticBondsInRing)

        return aromaticRings, aromaticBonds

    def getDeterministicSmallestSetOfSmallestRings(self):
        """
        Modified `Graph` method `getSmallestSetOfSmallestRings` by sorting calculated cycles
        by short lenth and then high atomic number instead of just short length (for cases where
        multiple cycles with same length are found, `getSmallestSetOfSmallestRings` outputs 
        non-determinstically ). 
        
        For instance, molecule with this SMILES: C1CC2C3CSC(CO3)C2C1, will have non-deterministic
        output from `getSmallestSetOfSmallestRings`, which leads to non-deterministic bycyclic decomposition
        Using this new method can effectively prevent this situation.
        """
        cython.declare(vertices=list, verticesToRemove=list, rootCandidates_tups=list, graphs=list)
        cython.declare(cycleList=list, cycleCandidate_tups=list, cycles=list, cycle0=list, originConnDict=dict)

        cython.declare(graph=Molecule, graph0=Molecule, vertex=Atom, rootVertex=Atom)
        
        # Make a copy of the graph so we don't modify the original
        graph = self.copy(deep=True)
        vertices = graph.vertices[:]
        
        # Step 1: Remove all terminal vertices
        done = False
        while not done:
            verticesToRemove = []
            for vertex in graph.vertices:
                if len(vertex.edges) == 1: verticesToRemove.append(vertex)
            done = len(verticesToRemove) == 0
            # Remove identified vertices from graph
            for vertex in verticesToRemove:
                graph.removeVertex(vertex)

        graph.updateConnectivityValues()
        # get original connectivity values
        originConnDict = {}
        for v in graph.vertices:
            originConnDict[v] = getVertexConnectivityValue(v)

        # Step 2: Remove all other vertices that are not part of cycles
        verticesToRemove = []
        for vertex in graph.vertices:
            found = graph.isVertexInCycle(vertex)
            if not found:
                verticesToRemove.append(vertex)
        # Remove identified vertices from graph
        for vertex in verticesToRemove:
            graph.removeVertex(vertex)

        # Step 3: Split graph into remaining subgraphs
        graphs = graph.split()

        # Step 4: Find ring sets in each subgraph
        cycleList = []
        for graph0 in graphs:

            while len(graph0.vertices) > 0:

                # Choose root vertex as vertex with smallest number of edges
                rootVertex = None
                graph0.updateConnectivityValues()

                rootCandidates_tups = []
                for vertex in graph0.vertices:
                    tup = (vertex, getVertexConnectivityValue(vertex), -originConnDict[vertex])
                    rootCandidates_tups.append(tup)

                rootVertex = sorted(rootCandidates_tups, key=lambda tup0: tup0[1:], reverse=True)[0][0]

                # Get all cycles involving the root vertex
                cycles = graph0.getAllCycles(rootVertex)
                if len(cycles) == 0:
                    # This vertex is no longer in a ring, so remove it
                    graph0.removeVertex(rootVertex)
                    continue

                # Keep the smallest of the cycles found above
                cycleCandidate_tups = []
                for cycle0 in cycles:
                    tup = (cycle0, len(cycle0), -sum([originConnDict[v] for v in cycle0]), -sum([v.element.number for v in cycle0]))
                    cycleCandidate_tups.append(tup)
                
                cycle = sorted(cycleCandidate_tups, key=lambda tup0: tup0[1:])[0][0]

                cycleList.append(cycle)

                # Remove the root vertex to create single edges, note this will not
                # function properly if there is no vertex with 2 edges (i.e. cubane)
                graph0.removeVertex(rootVertex)

                # Remove from the graph all vertices in the cycle that have only one edge
                loneCarbon = True
                while loneCarbon:
                    loneCarbon = False
                    verticesToRemove = []

                    for vertex in cycle:
                        if len(vertex.edges) == 1:
                            loneCarbon = True
                            verticesToRemove.append(vertex)
                    else:
                        for vertex in verticesToRemove:
                            graph0.removeVertex(vertex)

        # Map atoms in cycles back to atoms in original graph
        for i in range(len(cycleList)):
            cycleList[i] = [self.vertices[vertices.index(v)] for v in cycleList[i]]

        return cycleList

    def kekulize(self):
        """
        Kekulizes an aromatic molecule.
        """
        try:
            kekulize(self)
        except AtomTypeError:
            logging.error('Unable to kekulize molecule:/n{0}'.format(self.toAdjacencyList()))
