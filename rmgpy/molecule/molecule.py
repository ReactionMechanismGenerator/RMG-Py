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
import element as elements
try:
    import openbabel
except:
    pass
from rdkit import Chem
from .graph import Vertex, Edge, Graph
from .group import GroupAtom, GroupBond, Group, ActionError
from .atomtype import AtomType, atomTypes, getAtomType
import rmgpy.constants as constants

import numpy

#: This dictionary is used to shortcut lookups of a molecule's SMILES string from its chemical formula.
_known_smiles_molecules = {
                 'N2': 'N#N',
                 'CH4': 'C',
                 'H2O': 'O',
                 'C2H6': 'CC',
                 'H2': '[H][H]',
                 'H2O2': 'OO',
                 'C3H8': 'CCC',
                 'Ar': '[Ar]',
                 'He': '[He]',
                 'CH4O': 'CO',
                 'CO2': 'O=C=O',
                 'CO': '[C-]#[O+]',
                 'C2H4': 'C=C',
                 'O2': 'O=O'
             }

_known_smiles_radicals = {
                 'CH3': '[CH3]',
                 'HO': '[OH]',
                 'C2H5': 'C[CH2]',
                 'O': '[O]',
                 'HO2': '[O]O',
                 'CH': '[CH]',
                 'H': '[H]',
                 'C': '[C]',
                 #'CO2': it could be [O][C][O] or O=[C][O]
                 #'CO': '[C]=O', could also be [C][O]
                 #'C2H4': could  be [CH3][CH] or [CH2][CH2]
                 'O2': '[O][O]',
             }

################################################################################

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
        cython.declare(atom=Atom, ap=GroupAtom)
        if isinstance(other, Atom):
            atom = other
            return (
                self.element                is atom.element and
                self.radicalElectrons       == atom.radicalElectrons   and
                self.lonePairs              == atom.lonePairs           and
                self.charge                 == atom.charge
                )
        elif isinstance(other, GroupAtom):
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
        elif isinstance(other, GroupAtom):
            cython.declare(atom=GroupAtom, a=AtomType, radical=cython.short, lp = cython.short, charge=cython.short)
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

    def incrementRadical(self):
        """
        Update the atom pattern as a result of applying a GAIN_RADICAL action,
        where `radical` specifies the number of radical electrons to add.
        """
        # Set the new radical electron count
        self.radicalElectrons += 1
        if self.radicalElectrons <= 0:
            raise ActionError('Unable to update Atom due to GAIN_RADICAL action: Invalid radical electron set "{0}".'.format(self.radicalElectrons))

    def decrementRadical(self):
        """
        Update the atom pattern as a result of applying a LOSE_RADICAL action,
        where `radical` specifies the number of radical electrons to remove.
        """
        cython.declare(radicalElectrons=cython.short)
        # Set the new radical electron count
        radicalElectrons = self.radicalElectrons = self.radicalElectrons - 1
        if radicalElectrons  < 0:
            raise ActionError('Unable to update Atom due to LOSE_RADICAL action: Invalid radical electron set "{0}".'.format(self.radicalElectrons))

    def setLonePairs(self,lonePairs):
        """
        Set the number of lone electron pairs.
        """
        # Set the number of electron pairs
        self.lonePairs = lonePairs
        if self.lonePairs < 0:
            raise ActionError('Unable to update Atom due to setLonePairs : Invalid lone electron pairs set "{0}".'.format(self.setLonePairs))
        self.updateCharge()

    def incrementLonePairs(self):
        """
        Update the lone electron pairs pattern as a result of applying a GAIN_PAIR action.
        """
        # Set the new lone electron pairs count
        self.lonePairs += 1
        if self.lonePairs <= 0:
            raise ActionError('Unable to update Atom due to GAIN_PAIR action: Invalid lone electron pairs set "{0}".'.format(self.lonePairs))
        self.updateCharge()

    def decrementLonePairs(self):
        """
        Update the lone electron pairs pattern as a result of applying a LOSE_PAIR action.
        """
        # Set the new lone electron pairs count
        self.lonePairs -= 1
        if self.lonePairs  < 0:
            raise ActionError('Unable to update Atom due to LOSE_PAIR action: Invalid lone electron pairs set "{0}".'.format(self.lonePairs))
        self.updateCharge()
        
    def updateCharge(self):
        """
        Update self.charge, according to the valence, and the 
        number and types of bonds, radicals, and lone pairs.
        """
        valences = {'H': 1, 'C': 4, 'O': 2, 'N': 3, 'S': 2, 'Si': 4, 'He': 0, 'Ne': 0, 'Ar': 0, 'Cl': 1}
        orders = {'S': 1, 'D': 2, 'T': 3, 'B': 1.5}
        valence = valences[self.symbol]
        order = 0
        for atom2, bond in self.bonds.items():
            order += orders[bond.order]
        if self.symbol == 'H' or self.symbol == 'He':
            self.charge = 2 - valence - order - self.radicalElectrons - 2*self.lonePairs
        else:
            self.charge = 8 - valence - order - self.radicalElectrons - 2*self.lonePairs
        
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
            raise ActionError('Unable to update Atom: Invalid action {0}".'.format(action))
        
    def setSpinMultiplicity(self,spinMultiplicity):
        """
        Set the spin multiplicity.
        """
        raise NotImplementedError("I thought multiplicity was now a molecule attribute not atom?")
        # Set the spin multiplicity
        self.spinMultiplicity = spinMultiplicity
        if self.spinMultiplicity < 0:
            raise ActionError('Unable to update Atom due to spin multiplicity : Invalid spin multiplicity set "{0}".'.format(self.spinMultiplicity))
        self.updateCharge()
        

################################################################################

class Bond(Edge):
    """
    A chemical bond. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `order`             ``str``             The :ref:`bond type <bond-types>`
    =================== =================== ====================================

    """

    def __init__(self, atom1, atom2, order=1):
        Edge.__init__(self, atom1, atom2)
        self.order = order

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return self.order

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
        cython.declare(bond=Bond, bp=GroupBond)
        if isinstance(other, Bond):
            bond = other
            return (self.order == bond.order)
        elif isinstance(other, GroupBond):
            bp = other
            return (self.order in bp.order)

    def isSpecificCaseOf(self, other):
        """
        Return ``True`` if `self` is a specific case of `other`, or ``False``
        otherwise. `other` can be either a :class:`Bond` or a
        :class:`GroupBond` object.
        """
        # There are no generic bond types, so isSpecificCaseOf is the same as equivalent
        return self.equivalent(other)

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

    def isSingle(self):
        """
        Return ``True`` if the bond represents a single bond or ``False`` if
        not.
        """
        return self.order == 'S'

    def isDouble(self):
        """
        Return ``True`` if the bond represents a double bond or ``False`` if
        not.
        """
        return self.order == 'D'

    def isTriple(self):
        """
        Return ``True`` if the bond represents a triple bond or ``False`` if
        not.
        """
        return self.order == 'T'

    def isBenzene(self):
        """
        Return ``True`` if the bond represents a benzene bond or ``False`` if
        not.
        """
        return self.order == 'B'

    def incrementOrder(self):
        """
        Update the bond as a result of applying a CHANGE_BOND action to
        increase the order by one.
        """
        if self.order == 'S': self.order = 'D'
        elif self.order == 'D': self.order = 'T'
        else:
            raise ActionError('Unable to update Bond due to CHANGE_BOND action: Invalid bond order "{0}".'.format(self.order))
        
    def decrementOrder(self):
        """
        Update the bond as a result of applying a CHANGE_BOND action to
        decrease the order by one.
        """
        if self.order == 'D': self.order = 'S'
        elif self.order == 'T': self.order = 'D'
        else:
            raise ActionError('Unable to update Bond due to CHANGE_BOND action: Invalid bond order "{0}".'.format(self.order))
        
    def __changeBond(self, order):
        """
        Update the bond as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order, and should be 1 or -1.
        """
        if order == 1:
            if self.order == 'S': self.order = 'D'
            elif self.order == 'D': self.order = 'T'
            else:
                raise ActionError('Unable to update Bond due to CHANGE_BOND action: Invalid bond order "{0}".'.format(self.order))
        elif order == -1:
            if self.order == 'D': self.order = 'S'
            elif self.order == 'T': self.order = 'D'
            else:
                raise ActionError('Unable to update Bond due to CHANGE_BOND action: Invalid bond order "{0}".'.format(self.order))
        else:
            raise ActionError('Unable to update Bond due to CHANGE_BOND action: Invalid order "{0}".'.format(order))

    def applyAction(self, action):
        """
        Update the bond as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        if action[0].upper() == 'CHANGE_BOND':
            if action[2] == 1:
                self.incrementOrder()
            elif action[2] == -1:
                self.decrementOrder()
            elif action[2] == 'B':
                self.order = 'B'
            else:
                raise ActionError('Unable to update Bond due to CHANGE_BOND action: Invalid order "{0}".'.format(action[2]))
        else:
            raise ActionError('Unable to update GroupBond: Invalid action {0}.'.format(action))

#################################################################################
try:
    SMILEwriter = openbabel.OBConversion()
    SMILEwriter.SetOutFormat('smi')
    SMILEwriter.SetOptions("i",SMILEwriter.OUTOPTIONS) # turn off isomer and stereochemistry information (the @ signs!)
except:
    pass
    
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

    def __init__(self, atoms=None, symmetry=1, multiplicity=-187, SMILES='', InChI='', SMARTS=''):
        Graph.__init__(self, atoms)
        self.symmetryNumber = symmetry
        self.multiplicity = multiplicity
        self._fingerprint = None
        if SMILES != '': self.fromSMILES(SMILES)
        elif InChI != '': self.fromInChI(InChI)
        elif SMARTS != '': self.fromSMARTS(SMARTS)
        if multiplicity != -187:  # it was set explicitly, so re-set it (fromSMILES etc may have changed it)
            self.multiplicity = multiplicity
    
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
        return (Molecule, (self.vertices, self.symmetryNumber, self.multiplicity))

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

    def getFormula(self):
        """
        Return the molecular formula for the molecule.
        """
        cython.declare(atom=Atom, symbol=str, elements=dict, keys=list, formula=str)
        cython.declare(hasCarbon=cython.bint, hasHydrogen=cython.bint)
        
        # Count the number of each element in the molecule
        hasCarbon = False; hasHydrogen = False
        elements = {}
        for atom in self.vertices:
            symbol = atom.element.symbol
            elements[symbol] = elements.get(symbol, 0) + 1
        
        # Use the Hill system to generate the formula
        formula = ''
        
        # Carbon and hydrogen always come first if carbon is present
        if hasCarbon:
            count = elements['C']
            formula += 'C{0:d}'.format(count) if count > 1 else 'C'
            del elements['C']
            if hasHydrogen:
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
        Return the number of unpaired electrons.
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

    def getNumberOfRadicalElectrons(self):
        """
        Return the total number of radical electrons on all atoms in the
        molecule. In this function, monoradical atoms count as one, biradicals
        count as two, etc. 
        """
        cython.declare(numRadicals=cython.int, atom=Atom)
        numRadicals = 0
        for atom in self.vertices:
            numRadicals += atom.radicalElectrons
        return numRadicals

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
            assert atom.coords != None
        
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
                    # groupBond = GroupBond(atom1, atom2, ['S','D','T','B'])
                    bond = Bond(atom1, atom2, 'S')
                    self.addBond(bond)
        self.updateAtomTypes()
        
    def updateAtomTypes(self):
        """
        Iterate through the atoms in the structure, checking their atom types
        to ensure they are correct (i.e. accurately describe their local bond
        environment) and complete (i.e. are as detailed as possible).
        """
        for atom in self.vertices:
            atom.atomType = getAtomType(atom, atom.edges)
            
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
        cython.declare(group=Group, atom=Atom)
        cython.declare(carbonCount=cython.short, nitrogenCount=cython.short, oxygenCount=cython.short, sulfurCount=cython.short, radicalCount=cython.short)
        
        # It only makes sense to compare a Molecule to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
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
        # If the molecule has fewer of any of these things than the functional
        # group does, then we know the subgraph isomorphism fails without
        # needing to perform the full isomorphism check
        if (radicalCount < group.radicalCount or
            carbonCount < group.carbonCount or
            nitrogenCount < group.nitrogenCount or
            oxygenCount < group.oxygenCount or
            sulfurCount < group.sulfurCount or
            self.multiplicity not in group.multiplicity):
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
        # It only makes sense to compare a Molecule to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
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
        MoleculeDrawer().draw(self, format, path=path)
    
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

    def fromInChI(self, inchistr):
        """
        Convert an InChI string `inchistr` to a molecular structure. Uses
        `RDKit <http://rdkit.org/>`_ to perform the conversion.
        This Kekulizes everything, removing all aromatic atom types.
        """
        #RDkit was improperly handling the Hydrogen radical from InChI
        if inchistr == 'InChI=1/H' or inchistr == 'InChI=1S/H':
            self.fromSMILES('[H]')
            return self          
        elif inchistr == 'InChI=1/He' or inchistr == 'InChI=1S/He':
            self.fromSMILES('[He]')
            return self
        else:
            rdkitmol = Chem.inchi.MolFromInchi(inchistr)
            self.fromRDKitMol(rdkitmol)
            return self

    def fromSMILES(self, smilesstr):
        """
        Convert a SMILES string `smilesstr` to a molecular structure. Uses
        `RDKit <http://rdkit.org/>`_ to perform the conversion.
        This Kekulizes everything, removing all aromatic atom types.
        """
        # Special handling of helium
        if smilesstr == '[He]':
            # RDKit improperly handles helium and returns it in a triplet state
            self.fromAdjacencyList(
            """
            He
            multiplicity 1
            1 He u0 p1
            """)
            return self
        
        else:
            rdkitmol = Chem.MolFromSmiles(smilesstr)
            if rdkitmol is None:
                raise ValueError("Could not interpret the SMILES string {0!r}".format(smilesstr))
            self.fromRDKitMol(rdkitmol)
            return self
        
    def fromSMARTS(self, smartsstr):
        """
        Convert a SMARTS string `smartsstr` to a molecular structure. Uses
        `RDKit <http://rdkit.org/>`_ to perform the conversion.
        This Kekulizes everything, removing all aromatic atom types.
        """
        rdkitmol = Chem.MolFromSmarts(smartsstr)
        self.fromRDKitMol(rdkitmol)
        return self
        
    def fromRDKitMol(self, rdkitmol):
        """
        Convert a RDKit Mol object `rdkitmol` to a molecular structure. Uses
        `RDKit <http://rdkit.org/>`_ to perform the conversion.
        This Kekulizes everything, removing all aromatic atom types.
        """
        # Below are the declared variables for cythonizing the module
        cython.declare(i=cython.int)
        cython.declare(radicalElectrons=cython.int, charge=cython.int, lonePairs=cython.int)
        cython.declare(atom=Atom, atom1=Atom, atom2=Atom, bond=Bond)
        
        self.vertices = []
        
        # Add hydrogen atoms to complete molecule if needed
        rdkitmol = Chem.AddHs(rdkitmol)
        Chem.rdmolops.Kekulize(rdkitmol, clearAromaticFlags=True)
        
        # iterate through atoms in rdkitmol
        for i in range(rdkitmol.GetNumAtoms()):
            rdkitatom = rdkitmol.GetAtomWithIdx(i)
            
            # Use atomic number as key for element
            number = rdkitatom.GetAtomicNum()
            element = elements.getElement(number)
                
            # Process charge
            charge = rdkitatom.GetFormalCharge()
            radicalElectrons = rdkitatom.GetNumRadicalElectrons()
            
            atom = Atom(element, radicalElectrons, charge, '', 0)
            self.vertices.append(atom)
            
            # Add bonds by iterating again through atoms
            for j in range(0, i):
                rdkitatom2 = rdkitmol.GetAtomWithIdx(j + 1)
                rdkitbond = rdkitmol.GetBondBetweenAtoms(i, j)
                if rdkitbond is not None:
                    order = 0
        
                    # Process bond type
                    rdbondtype = rdkitbond.GetBondType()
                    if rdbondtype.name == 'SINGLE': order = 'S'
                    elif rdbondtype.name == 'DOUBLE': order = 'D'
                    elif rdbondtype.name == 'TRIPLE': order = 'T'
                    elif rdbondtype.name == 'AROMATIC': order = 'B'
        
                    bond = Bond(self.vertices[i], self.vertices[j], order)
                    self.addBond(bond)
        
        # Set atom types and connectivity values
        self.updateConnectivityValues()
        self.updateLonePairs()
        self.updateAtomTypes()
        
        # Assume this is always true
        # There are cases where 2 radicalElectrons is a singlet, but
        # the triplet is often more stable, 
        self.multiplicity = self.getRadicalCount() + 1
        
        return self

    def fromAdjacencyList(self, adjlist, saturateH=False):
        """
        Convert a string adjacency list `adjlist` to a molecular structure.
        Skips the first line (assuming it's a label) unless `withLabel` is
        ``False``.
        """
        from .adjlist import fromAdjacencyList
        
        self.vertices, self.multiplicity = fromAdjacencyList(adjlist, group=False, saturateH=saturateH)
        self.updateConnectivityValues()
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
                bond = Bond(mapping[atom1], mapping[atom2], 'S')
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
        try:
            if not Chem.inchi.INCHI_AVAILABLE:
                return "RDKitInstalledWithoutInChI"
            rdkitmol = self.toRDKitMol()
            return Chem.inchi.MolToInchi(rdkitmol, options='-SNon')
        except:
            pass

        obmol = self.toOBMol()
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat('inchi')
        obConversion.SetOptions('w', openbabel.OBConversion.OUTOPTIONS)
        return obConversion.WriteString(obmol).strip()

    def toAugmentedInChI(self):
        """
        Adds an extra layer to the InChI denoting the number of unpaired electrons in case
        more than 1 ( >= 2) unpaired electrons are present in the molecule.
        """
        inchi = self.toInChI()
        
        radicalNumber = sum([i.radicalElectrons for i in self.atoms])
        
        if radicalNumber >= 2:
            return inchi+'/mult'+str(radicalNumber+1)
        else:
            return inchi
    
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
        try:
            if not Chem.inchi.INCHI_AVAILABLE:
                return "RDKitInstalledWithoutInChI"
            inchi = self.toInChI()
            return Chem.inchi.InchiToInchiKey(inchi)[:-2]
        except:
            pass
        
        import openbabel

#        for atom in self.vertices:
 #           if atom.isNitrogen():
        obmol = self.toOBMol()
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat('inchi')
        obConversion.SetOptions('w', openbabel.OBConversion.OUTOPTIONS)
        obConversion.SetOptions('K', openbabel.OBConversion.OUTOPTIONS)
        return obConversion.WriteString(obmol).strip()[:-2]
    
    def toAugmentedInChIKey(self):
        """
        Adds an extra layer to the InChIKey denoting the number of unpaired electrons in case
        more than 1 ( >= 2) unpaired electrons are present in the molecule.
        """
        key = self.toInChIKey()
        
        radicalNumber = sum([i.radicalElectrons for i in self.atoms])
        
        if radicalNumber >= 2:
            return key+'mult'+str(radicalNumber+1)
        else:
            return key

    def toSMARTS(self):
        """
        Convert a molecular structure to an SMARTS string. Uses
        `RDKit <http://rdkit.org/>`_ to perform the conversion.
        Perceives aromaticity and removes Hydrogen atoms.
        """
        rdkitmol = self.toRDKitMol()
        
        return Chem.MolToSmarts(rdkitmol)
    
    
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
        
        # If we're going to have to check the formula anyway,
        # we may as well shortcut a few small known molecules.
        # Dictionary lookups are O(1) so this should be fast:
        # The dictionary is defined at the top of this file.
        try:
            if self.isRadical():
                return _known_smiles_radicals[self.getFormula()]
            else:
                return _known_smiles_molecules[self.getFormula()]
        except KeyError:
            # It wasn't in the above list.
            pass
        for atom in self.vertices:
            if atom.isNitrogen():
                mol = self.toOBMol()
                return SMILEwriter.WriteString(mol).strip()

        rdkitmol = self.toRDKitMol()
        return Chem.MolToSmiles(rdkitmol)

    def toOBMol(self):
        """
        Convert a molecular structure to an OpenBabel OBMol object. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        cython.declare(atom=Atom, atom1=Atom, bonds=dict, atom2=Atom, bond=Bond)
        cython.declare(index1=cython.int, index2=cython.int, order=cython.int)

        # Sort the atoms before converting to ensure output is consistent
        # between different runs
        self.sortAtoms()

        atoms = self.vertices

        obmol = openbabel.OBMol()
        for atom in atoms:
            a = obmol.NewAtom()
            a.SetAtomicNum(atom.number)
            a.SetFormalCharge(atom.charge)
        orders = {'S': 1, 'D': 2, 'T': 3, 'B': 5}
        for atom1 in self.vertices:
            for atom2, bond in atom1.edges.iteritems():
                index1 = atoms.index(atom1)
                index2 = atoms.index(atom2)
                if index1 < index2:
                    order = orders[bond.order]
                    obmol.AddBond(index1+1, index2+1, order)

        obmol.AssignSpinMultiplicity(True)

        return obmol
    
    def toRDKitMol(self, removeHs=True, returnMapping=False):
        """
        Convert a molecular structure to a RDKit rdmol object. Uses
        `RDKit <http://rdkit.org/>`_ to perform the conversion.
        Perceives aromaticity and, unless removeHs==False, removes Hydrogen atoms.
        
        If returnMapping==True then it also returns a dictionary mapping the 
        atoms to RDKit's atom indices.
        """
        # Sort the atoms before converting to ensure output is consistent
        # between different runs
        self.sortAtoms()
        atoms = self.vertices
        rdAtomIndices = {} # dictionary of RDKit atom indices
        rdkitmol = Chem.rdchem.EditableMol(Chem.rdchem.Mol())
        for index, atom in enumerate(self.vertices):
            rdAtom = Chem.rdchem.Atom(atom.element.symbol)
            rdAtom.SetNumRadicalElectrons(atom.radicalElectrons)
            rdkitmol.AddAtom(rdAtom)
            rdAtomIndices[atom] = index
        
        rdBonds = Chem.rdchem.BondType
        orders = {'S': rdBonds.SINGLE, 'D': rdBonds.DOUBLE, 'T': rdBonds.TRIPLE, 'B': rdBonds.AROMATIC}
        # Add the bonds
        for atom1 in self.vertices:
            for atom2, bond in atom1.edges.iteritems():
                index1 = atoms.index(atom1)
                index2 = atoms.index(atom2)
                if index1 < index2:
                    order = orders[bond.order]
                    rdkitmol.AddBond(index1, index2, order)
        
        # Make editable mol into a mol and rectify the molecule
        rdkitmol = rdkitmol.GetMol()
        Chem.SanitizeMol(rdkitmol)
        if removeHs:
            rdkitmol = Chem.RemoveHs(rdkitmol)
        
        if returnMapping:
            return rdkitmol, rdAtomIndices
        return rdkitmol

    def toAdjacencyList(self, label='', removeH=False, removeLonePairs=False):
        """
        Convert the molecular structure to a string adjacency list.
        """
        from .adjlist import toAdjacencyList
        result = toAdjacencyList(self.vertices, self.multiplicity,  label=label, group=False, removeH=removeH, removeLonePairs=removeLonePairs)
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
    
    def generateResonanceIsomers(self):
        """
        Generate and return all of the resonance isomers of this molecule.
        """
        cython.declare(isomers=list, newIsomers=list, index=cython.int, atom=Atom)
        cython.declare(isomer=Molecule, newIsomer=Molecule, isom=Molecule)
        
        isomers = [self]

        # Iterate over resonance isomers
        index = 0
        while index < len(isomers):
            isomer = isomers[index]
                
            newIsomers = isomer.getAdjacentResonanceIsomers()
            newIsomers += isomer.getLonePairRadicalResonanceIsomers()
            newIsomers += isomer.getN5dd_N5tsResonanceIsomers()
            for newIsomer in newIsomers:
                newIsomer.updateAtomTypes()
                # Append to isomer list if unique
                for isom in isomers:
                    if isom.isIsomorphic(newIsomer):
                        break
                else:
                    isomers.append(newIsomer)
                        
            # Move to next resonance isomer
            index += 1
        
        return isomers

    def getAdjacentResonanceIsomers(self):
        """
        Generate all of the resonance isomers formed by one allyl radical shift.
        """
        cython.declare(isomers=list, paths=list, index=cython.int, isomer=Molecule)
        cython.declare(atom=Atom, atom1=Atom, atom2=Atom, atom3=Atom, bond12=Bond, bond23=Bond)
        cython.declare(v1=Vertex, v2=Vertex)
        
        isomers = []

        # Radicals
        if self.isRadical():
            # Iterate over radicals in structure
            for atom in self.vertices:
                paths = self.findAllDelocalizationPaths(atom)
                for atom1, atom2, atom3, bond12, bond23 in paths:
                    # Adjust to (potentially) new resonance isomer
                    atom1.decrementRadical()
                    atom3.incrementRadical()
                    bond12.incrementOrder()
                    bond23.decrementOrder()
                    # Make a copy of isomer
                    isomer = self.copy(deep=True)
                    # Also copy the connectivity values, since they are the same
                    # for all resonance forms
                    for index in range(len(self.vertices)):
                        v1 = self.vertices[index]
                        v2 = isomer.vertices[index]
                        v2.connectivity1 = v1.connectivity1
                        v2.connectivity2 = v1.connectivity2
                        v2.connectivity3 = v1.connectivity3
                        v2.sortingLabel = v1.sortingLabel
                    # Restore current isomer
                    atom1.incrementRadical()
                    atom3.decrementRadical()
                    bond12.decrementOrder()
                    bond23.incrementOrder()
                    # Append to isomer list if unique
                    isomers.append(isomer)

        return isomers
    
    def getLonePairRadicalResonanceIsomers(self):
        """
        Generate all of the resonance isomers formed by lone electron pair - radical shifts.
        """
        cython.declare(isomers=list, paths=list, index=cython.int, isomer=Molecule)
        cython.declare(atom=Atom, atom1=Atom, atom2=Atom)
        cython.declare(v1=Vertex, v2=Vertex)
        
        isomers = []

        # Radicals
        if self.isRadical():
            # Iterate over radicals in structure
            for atom in self.vertices:
                paths = self.findAllDelocalizationPathsLonePairRadical(atom)
                for atom1, atom2 in paths:
                    # Adjust to (potentially) new resonance isomer
                    atom1.decrementRadical()
                    atom1.incrementLonePairs()
                    atom1.updateCharge()
                    atom2.incrementRadical()
                    atom2.decrementLonePairs()
                    atom2.updateCharge()
                    # Make a copy of isomer
                    isomer = self.copy(deep=True)
                    # Also copy the connectivity values, since they are the same
                    # for all resonance forms
                    for index in range(len(self.vertices)):
                        v1 = self.vertices[index]
                        v2 = isomer.vertices[index]
                        v2.connectivity1 = v1.connectivity1
                        v2.connectivity2 = v1.connectivity2
                        v2.connectivity3 = v1.connectivity3
                        v2.sortingLabel = v1.sortingLabel
                    # Restore current isomer
                    atom1.incrementRadical()
                    atom1.decrementLonePairs()
                    atom1.updateCharge()
                    atom2.decrementRadical()
                    atom2.incrementLonePairs()
                    atom2.updateCharge()
                    # Append to isomer list if unique
                    isomers.append(isomer)

        return isomers
    
    def getN5dd_N5tsResonanceIsomers(self):
        """
        Generate all of the resonance isomers formed by shifts between N5dd and N5ts.
        """
        cython.declare(isomers=list, paths=list, index=cython.int, isomer=Molecule)
        cython.declare(atom=Atom, atom1=Atom, atom2=Atom, atom3=Atom)
        cython.declare(bond12=Bond, bond13=Bond)
        cython.declare(v1=Vertex, v2=Vertex)
        
        isomers = []
        
        # Iterate over nitrogen atoms in structure
        for atom in self.vertices:
            paths = self.findAllDelocalizationPathsN5dd_N5ts(atom)
            for atom1, atom2, atom3, bond12, bond13, direction in paths:
                # from N5dd to N5ts
                if direction == 1:
                    # Adjust to (potentially) new resonance isomer
                    bond12.decrementOrder()
                    bond13.incrementOrder()
                    atom2.incrementLonePairs()
                    atom3.decrementLonePairs()
                    atom1.updateCharge()
                    atom2.updateCharge()
                    atom3.updateCharge()
                    # Make a copy of isomer
                    isomer = self.copy(deep=True)
                    # Also copy the connectivity values, since they are the same
                    # for all resonance forms
                    for index in range(len(self.vertices)):
                        v1 = self.vertices[index]
                        v2 = isomer.vertices[index]
                        v2.connectivity1 = v1.connectivity1
                        v2.connectivity2 = v1.connectivity2
                        v2.connectivity3 = v1.connectivity3
                        v2.sortingLabel = v1.sortingLabel
                    # Restore current isomer
                    bond12.incrementOrder()
                    bond13.decrementOrder()
                    atom2.decrementLonePairs()
                    atom3.incrementLonePairs()
                    atom1.updateCharge()
                    atom2.updateCharge()
                    atom3.updateCharge()
                    # Append to isomer list if unique
                    isomers.append(isomer)
                
                # from N5ts to N5dd
                if direction == 2:
                    # Adjust to (potentially) new resonance isomer
                    bond12.decrementOrder()
                    bond13.incrementOrder()
                    atom2.incrementLonePairs()
                    atom3.decrementLonePairs()
                    atom1.updateCharge()
                    atom2.updateCharge()
                    atom3.updateCharge()
                    # Make a copy of isomer
                    isomer = self.copy(deep=True)
                    # Also copy the connectivity values, since they are the same
                    # for all resonance forms
                    for index in range(len(self.vertices)):
                        v1 = self.vertices[index]
                        v2 = isomer.vertices[index]
                        v2.connectivity1 = v1.connectivity1
                        v2.connectivity2 = v1.connectivity2
                        v2.connectivity3 = v1.connectivity3
                        v2.sortingLabel = v1.sortingLabel
                    # Restore current isomer
                    bond12.incrementOrder()
                    bond13.decrementOrder()
                    atom2.decrementLonePairs()
                    atom3.incrementLonePairs()
                    atom1.updateCharge()
                    atom2.updateCharge()
                    atom3.updateCharge()
                    # Append to isomer list if unique
                    isomers.append(isomer)
                    
        return isomers

    def findAllDelocalizationPaths(self, atom1):
        """
        Find all the delocalization paths allyl to the radical center indicated
        by `atom1`. Used to generate resonance isomers.
        """
        cython.declare(paths=list)
        cython.declare(atom2=Atom, atom3=Atom, bond12=Bond, bond23=Bond)
        
        # No paths if atom1 is not a radical
        if atom1.radicalElectrons <= 0:
            return []

        # Find all delocalization paths
        paths = []
        for atom2, bond12 in atom1.edges.items():
            # Vinyl bond must be capable of gaining an order
            if (bond12.isSingle() or bond12.isDouble()) and atom1.radicalElectrons == 1:
                for atom3, bond23 in atom2.edges.items():
                    # Allyl bond must be capable of losing an order without breaking
                    if atom1 is not atom3 and (bond23.isDouble() or bond23.isTriple()):
                        paths.append([atom1, atom2, atom3, bond12, bond23])
        return paths
    
    def findAllDelocalizationPathsLonePairRadical(self, atom1):
        """
        Find all the delocalization paths of lone electron pairs next to the radical center indicated
        by `atom1`. Used to generate resonance isomers.
        """
        cython.declare(paths=list)
        cython.declare(atom2=Atom, bond12=Bond)
        
        # No paths if atom1 is not a radical
        if atom1.radicalElectrons <= 0:
            return []
        
        # In a first step we only consider nitrogen and oxygen atoms as possible radical centers
        if not ((atom1.lonePairs == 0 and atom1.isNitrogen()) or(atom1.lonePairs == 2 and atom1.isOxygen())):
            return []
        
        # Find all delocalization paths
        paths = []
        for atom2, bond12 in atom1.edges.items():
            # Only single bonds are considered
            if bond12.isSingle():
                # Neighboring atom must posses a lone electron pair to loose it
                if ((atom2.lonePairs == 1 and atom2.isNitrogen()) or (atom2.lonePairs == 3 and atom2.isOxygen())) and (atom2.radicalElectrons == 0):
                    paths.append([atom1, atom2])
                    
        return paths
    
    def findAllDelocalizationPathsN5dd_N5ts(self, atom1):
        """
        Find all the resonance structures of nitrogen atoms with two double bonds (N5dd)
        and nitrogen atoms with one triple and one single bond (N5ts)
        """
        cython.declare(paths=list)
        cython.declare(atom2=Atom, bond12=Bond)
        
        # No paths if atom1 is not nitrogen
        if not (atom1.isNitrogen()):
            return []
        
        # Find all delocalization paths
        paths = []
        index_atom_2 = 0
        index_atom_3 = 0
        
        for atom2, bond12 in atom1.edges.items():
            index_atom_2 = index_atom_2 + 1
            # Only double bonds are considered
            if bond12.isDouble():
                for atom3, bond13 in atom1.edges.items():
                    index_atom_3 = index_atom_3 + 1
                    # Only double bonds are considered, at the moment we only consider non-radical nitrogen and oxygen atoms
                    if (bond13.isDouble() and atom3.radicalElectrons == 0 and atom3.lonePairs > 0 and not atom3.isOxygen() and not atom3.isCarbon() and (index_atom_2 != index_atom_3)):
                        paths.append([atom1, atom2, atom3, bond12, bond13, 1])
        
        for atom2, bond12 in atom1.edges.items():
            # Only triple bonds are considered
            if bond12.isTriple():
                for atom3, bond13 in atom1.edges.items():
                    # Only single bonds are considered, at the moment we only consider negatively charged nitrogen and oxygen
                    if (bond13.isSingle() and ((atom3.isNitrogen() and atom3.lonePairs >= 2) or (atom3.isOxygen() and atom3.lonePairs >= 3))):
                        paths.append([atom1, atom2, atom3, bond12, bond13, 2])
        
        return paths

    def getURL(self):
        """
        Get a URL to the molecule's info page on the RMG website.
        """
        # eg. http://dev.rmg.mit.edu/database/kinetics/reaction/reactant1=1%20C%200%20%7B2,S%7D;2%20O%200%20%7B1,S%7D;__reactant2=1%20C%202T;__product1=1%20C%201;__product2=1%20C%200%20%7B2,S%7D;2%20O%201%20%7B1,S%7D;

        url = "http://rmg.mit.edu/database/molecule/"
        adjlist = self.toAdjacencyList(removeH=False)
        url += "{0}".format(re.sub('\s+', '%20', adjlist.replace('\n', ';')))
        return url.strip('_')
    
    def isBiradicalSinglet(self):
        """
        Return ``True`` if the molecule is a 1-centered biradical, and the molecule is in singlet state,
        or ``False`` otherwise.
        """
        cython.declare(atom=Atom)
        for atom in self.vertices:
            if atom.radicalElectrons == 2:
                if self.multiplicity == 1:
                    return True
        return False
    
    def isBiradicalTriplet(self):
        """
        Return ``True`` if the molecule is a 1-centered biradical, and the molecule is in triplet state,
        or ``False`` otherwise.
        """
        cython.declare(atom=Atom)
        for atom in self.vertices:
            if atom.radicalElectrons == 2:
                if self.multiplicity == 3:
                    return True
        return False
    
    def changeTripletSinglet(self):
        """
        If the molecule is a 1-centered biradical, and the molecule is in a triplet state,
        change it to singlet state.
        """
        cython.declare(atom=Atom)
        for atom in self.vertices:
            if atom.radicalElectrons == 2:
                if self.multiplicity == 3:
                    self.multiplicity = 1
                    
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
        Iterate through the atoms in the structure and calcualte the
        number of lone electron pairs, assumin a neutral molecule.
        """
        for atom1 in self.vertices:
            order = 0
            if not atom1.isHydrogen():
                for atom2, bond12 in atom1.edges.items():
                    if bond12.isSingle():
                        order = order + 1
                    if bond12.isDouble():
                        order = order + 2
                    if bond12.isTriple():
                        order = order + 3
                        
                atom1.lonePairs = 4 - atom1.radicalElectrons - order
        
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

