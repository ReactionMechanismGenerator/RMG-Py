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
import os
import re
import element as elements
from .graph import Vertex, Edge, Graph
from .group import GroupAtom, GroupBond, Group, ActionError, fromAdjacencyList, toAdjacencyList
from .atomtype import AtomType, atomTypes, getAtomType
        
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
    `spinMultiplicity`  ``short``           The spin multiplicity of the atom
    `charge`            ``short``           The formal charge of the atom
    `label`             ``str``             A string label that can be used to tag individual atoms
    =================== =================== ====================================

    Additionally, the ``mass``, ``number``, and ``symbol`` attributes of the
    atom's element can be read (but not written) directly from the atom object,
    e.g. ``atom.symbol`` instead of ``atom.element.symbol``.
    """

    def __init__(self, element=None, radicalElectrons=0, spinMultiplicity=1, charge=0, label=''):
        Vertex.__init__(self)
        if isinstance(element, str):
            self.element = elements.__dict__[element]
        else:
            self.element = element
        self.radicalElectrons = radicalElectrons
        self.spinMultiplicity = spinMultiplicity
        self.charge = charge
        self.label = label
        self.atomType = None

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return "<Atom '{0}'>".format(
            str(self.element) +
            '.' * self.radicalElectrons +
            '+' * self.charge if self.charge > 0 else '-' * -self.charge
        )

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return 'Atom(element="{0}", radicalElectrons={1}, spinMultiplicity={2}, charge={3}, label="{4}")'.format(self.element, self.radicalElectrons, self.spinMultiplicity, self.charge, self.label)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        d = {
            'connectivity1': self.connectivity1,
            'connectivity2': self.connectivity2,
            'connectivity3': self.connectivity3,
            'sortingLabel': self.sortingLabel,
            'atomType': self.atomType.label if self.atomType else None,
        }
        return (Atom, (self.element.symbol, self.radicalElectrons, self.spinMultiplicity, self.charge, self.label), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling an object.
        """
        self.connectivity1 = d['connectivity1']
        self.connectivity2 = d['connectivity2']
        self.connectivity3 = d['connectivity3']
        self.sortingLabel = d['sortingLabel']
        self.atomType = atomTypes[d['atomType']] if d['atomType'] else None
    
    @property
    def mass(self): return self.element.mass
    
    @property
    def number(self): return self.element.number

    @property
    def symbol(self): return self.element.symbol

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
            return (self.element is atom.element and
                self.radicalElectrons == atom.radicalElectrons and
                self.spinMultiplicity == atom.spinMultiplicity and
                self.charge == atom.charge)
        elif isinstance(other, GroupAtom):
            cython.declare(a=AtomType, radical=cython.short, spin=cython.short, charge=cython.short)
            ap = other
            for a in ap.atomType:
                if self.atomType.equivalent(a): break
            else:
                return False
            for radical, spin in zip(ap.radicalElectrons, ap.spinMultiplicity):
                if self.radicalElectrons == radical and self.spinMultiplicity == spin: break
            else:
                return False
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
            cython.declare(atom=GroupAtom, a=AtomType, radical=cython.short, spin=cython.short, charge=cython.short)
            atom = other
            for a in atom.atomType: 
                if self.atomType.isSpecificCaseOf(a): break
            else:
                return False
            for radical, spin in zip(atom.radicalElectrons, atom.spinMultiplicity):
                if self.radicalElectrons == radical and self.spinMultiplicity == spin: break
            else:
                return False
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
        a = Atom(self.element, self.radicalElectrons, self.spinMultiplicity, self.charge, self.label)
        a.atomType = self.atomType
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
        # Set the new radical electron counts and spin multiplicities
        self.radicalElectrons += 1
        if self.radicalElectrons <= 0:
            raise ActionError('Unable to update Atom due to GAIN_RADICAL action: Invalid radical electron set "{0}".'.format(self.radicalElectrons))
        elif self.radicalElectrons == 1:
            self.spinMultiplicity = 2
        elif self.radicalElectrons == 2:
            self.spinMultiplicity = 3   # Assume this always results in the triplet, as they tend to be more stable than the singlet (though there are exceptions!)
        else:
            self.spinMultiplicity = self.radicalElectrons + 1

    def decrementRadical(self):
        """
        Update the atom pattern as a result of applying a LOSE_RADICAL action,
        where `radical` specifies the number of radical electrons to remove.
        """
        # Set the new radical electron counts and spin multiplicities
        self.radicalElectrons -= 1
        if self.radicalElectrons  < 0:
            raise ActionError('Unable to update Atom due to LOSE_RADICAL action: Invalid radical electron set "{0}".'.format(self.radicalElectrons))
        elif self.radicalElectrons == 0:
            self.spinMultiplicity = 1
        elif self.radicalElectrons == 1:
            self.spinMultiplicity = 2
        elif self.radicalElectrons == 2:
            self.spinMultiplicity = 3   # Assume this always results in the triplet, as they tend to be more stable than the singlet (though there are exceptions!)
        else:
            self.spinMultiplicity = self.radicalElectrons + 1

    def applyAction(self, action):
        """
        Update the atom pattern as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        # Invalidate current atom type
        self.atomType = None
        # Modify attributes if necessary
        if action[0].upper() in ['CHANGE_BOND', 'FORM_BOND', 'BREAK_BOND']:
            # Nothing else to do here
            pass
        elif action[0].upper() == 'GAIN_RADICAL':
            for i in range(action[2]): self.incrementRadical()
        elif action[0].upper() == 'LOSE_RADICAL':
            for i in range(abs(action[2])): self.decrementRadical()
        else:
            raise ActionError('Unable to update Atom: Invalid action {0}".'.format(action))

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

    def __init__(self, order=1):
        Edge.__init__(self)
        self.order = order

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return '<Bond "{0}">'.format(self.order)

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return 'Bond(order="{0}")'.format(self.order)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Bond, (self.order,))

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
        return Bond(self.order)

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
            else:
                raise ActionError('Unable to update Bond due to CHANGE_BOND action: Invalid order "{0}".'.format(action[2]))
        else:
            raise ActionError('Unable to update GroupBond: Invalid action {0}.'.format(action))

################################################################################

class Molecule(Graph):
    """
    A representation of a molecular structure using a graph data type, extending
    the :class:`Graph` class. The `atoms` and `bonds` attributes are aliases
    for the `vertices` and `edges` attributes. Other attributes are:

    ======================= =========== ========================================
    Attribute               Type        Description
    ======================= =========== ========================================
    `symmetryNumber`        ``int``     The (estimated) external + internal symmetry number of the molecule
    ======================= =========== ========================================

    A new molecule object can be easily instantiated by passing the `SMILES` or
    `InChI` string representing the molecular structure.
    """

    def __init__(self, atoms=None, bonds=None, symmetry=1, SMILES='', InChI=''):
        Graph.__init__(self, atoms, bonds)
        self.symmetryNumber = symmetry
        if SMILES != '': self.fromSMILES(SMILES)
        elif InChI != '': self.fromInChI(InChI)
    
    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return '<Molecule "{0}">'.format(self.toSMILES())

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return 'Molecule(SMILES="{0}")'.format(self.toSMILES())

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Molecule, (self.vertices, self.edges, self.symmetryNumber))

    def __getAtoms(self): return self.vertices
    def __setAtoms(self, atoms): self.vertices = atoms
    atoms = property(__getAtoms, __setAtoms)

    def __getBonds(self): return self.edges
    def __setBonds(self, bonds): self.edges = bonds
    bonds = property(__getBonds, __setBonds)

    def addAtom(self, atom):
        """
        Add an `atom` to the graph. The atom is initialized with no bonds.
        """
        return self.addVertex(atom)
    
    def addBond(self, atom1, atom2, bond):
        """
        Add a `bond` to the graph as an edge connecting the two atoms `atom1`
        and `atom2`.
        """
        return self.addEdge(atom1, atom2, bond)

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
        return self.removeVertex(atom)

    def removeBond(self, atom1, atom2):
        """
        Remove the bond between atoms `atom1` and `atom2` from the graph.
        Does not remove atoms that no longer have any bonds as a result of
        this removal.
        """
        return self.removeEdge(atom1, atom2)

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
        import pybel
        mol = pybel.Molecule(self.toOBMol())
        return mol.formula

    def getMolecularWeight(self):
        """
        Return the molecular weight of the molecule in kg/mol.
        """
        mass = 0
        H = elements.getElement('H')
        for atom in self.vertices:
            mass += atom.element.mass
        return mass

    def copy(self, deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """
        other = cython.declare(Molecule)
        g = Graph.copy(self, deep)
        other = Molecule(g.vertices, g.edges)
        return other

    def merge(self, other):
        """
        Merge two molecules so as to store them in a single :class:`Molecule`
        object. The merged :class:`Molecule` object is returned.
        """
        g = Graph.merge(self, other)
        molecule = Molecule(atoms=g.vertices, bonds=g.edges)
        return molecule

    def split(self):
        """
        Convert a single :class:`Molecule` object containing two or more
        unconnected molecules into separate class:`Molecule` objects.
        """
        graphs = Graph.split(self)
        molecules = []
        for g in graphs:
            molecule = Molecule(atoms=g.vertices, bonds=g.edges)
            molecules.append(molecule)
        return molecules

    def deleteHydrogens(self):
        """
        Irreversibly delete all non-labeled hydrogens without updating
        connectivity values. If there's nothing but hydrogens, it does nothing.
        It destroys information; be careful with it.
        """
        cython.declare(atom=Atom, neighbor=Atom, hydrogens=list)
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
                neighbor = self.edges[atom].keys()[0]
                hydrogens.append(atom)
        # Remove the hydrogen atoms from the structure
        for atom in hydrogens:
            self.removeAtom(atom)

    def updateAtomTypes(self):
        """
        Iterate through the atoms in the structure, checking their atom types
        to ensure they are correct (i.e. accurately describe their local bond
        environment) and complete (i.e. are as detailed as possible).
        """
        for atom in self.vertices:
            atom.atomType = getAtomType(atom, self.edges[atom])

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

    def isIsomorphic(self, other, initialMap=None):
        """
        Returns :data:`True` if two graphs are isomorphic and :data:`False`
        otherwise. The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Molecule` object, or a :class:`TypeError` is raised.
        """
        # It only makes sense to compare a Molecule to a Molecule for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Molecule):
            raise TypeError('Got a {0} object for parameter "other", when a Molecule object is required.'.format(other.__class__))
        # Do the isomorphism comparison
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
        # It only makes sense to compare a Molecule to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Molecule object is required.'.format(other.__class__))
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
            raise TypeError('Got a {0} object for parameter "other", when a Molecule object is required.'.format(other.__class__))
        # Do the isomorphism comparison
        result = Graph.findSubgraphIsomorphisms(self, other, initialMap)
        return result

    def isAtomInCycle(self, atom):
        """
        Return :data:`True` if `atom` is in one or more cycles in the structure,
        and :data:`False` if not.
        """
        return self.isVertexInCycle(atom)

    def isBondInCycle(self, atom1, atom2):
        """
        Return :data:`True` if the bond between atoms `atom1` and `atom2`
        is in one or more cycles in the graph, or :data:`False` if not.
        """
        return self.isEdgeInCycle(atom1, atom2)

    def draw(self, path):
        """
        Generate a pictorial representation of the chemical graph using the
        :mod:`ext.molecule_draw` module. Use `path` to specify the file to save
        the generated image to; the image type is automatically determined by
        extension. Valid extensions are ``.png``, ``.svg``, ``.pdf``, and
        ``.ps``; of these, the first is a raster format and the remainder are
        vector formats.
        """
        from molecule_draw import drawMolecule
        drawMolecule(self, path=path)
    
    def _repr_png_(self):
        """
        Return a png picture of the molecule, useful for ipython-qtconsole.
        """
        from molecule_draw import drawMolecule
        tempFileName = 'temp_molecule.png'
        drawMolecule(self, path=tempFileName)
        png = open(tempFileName,'rb').read()
        os.unlink(tempFileName)
        return png
        

    def fromCML(self, cmlstr):
        """
        Convert a string of CML `cmlstr` to a molecular structure. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        import pybel
        cmlstr = cmlstr.replace('\t', '')
        mol = pybel.readstring('cml', cmlstr)
        self.fromOBMol(mol.OBMol)
        return self

    def fromInChI(self, inchistr):
        """
        Convert an InChI string `inchistr` to a molecular structure. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        if inchistr in ['InChI=1/H', 'InChI=1S/H']:
            # cope with a bug in OpenBabel
            return self.fromAdjacencyList('1 H 1')
        import pybel
        mol = pybel.readstring('inchi', inchistr)
        self.fromOBMol(mol.OBMol)
        return self

    def fromSMILES(self, smilesstr):
        """
        Convert a SMILES string `smilesstr` to a molecular structure. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        if smilesstr == '[H]':
            # cope with a bug in OpenBabel < 2.3
            return self.fromAdjacencyList('1 H 1')
        import pybel
        mol = pybel.readstring('smiles', smilesstr)
        self.fromOBMol(mol.OBMol)
        return self

    def fromOBMol(self, obmol):
        """
        Convert an OpenBabel OBMol object `obmol` to a molecular structure. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """

        cython.declare(i=cython.int)
        cython.declare(radicalElectrons=cython.int, spinMultiplicity=cython.int, charge=cython.int)
        cython.declare(atom=Atom, atom1=Atom, atom2=Atom, bond=Bond)

        self.vertices = []
        self.edges = {}

        # Add hydrogen atoms to complete molecule if needed
        obmol.AddHydrogens()

        # Iterate through atoms in obmol
        for i in range(0, obmol.NumAtoms()):
            obatom = obmol.GetAtom(i + 1)

            # Use atomic number as key for element
            number = obatom.GetAtomicNum()
            element = elements.getElement(number)
            
            # Process spin multiplicity
            radicalElectrons = 0
            spinMultiplicity = obatom.GetSpinMultiplicity()
            if spinMultiplicity == 0:
                radicalElectrons = 0; spinMultiplicity = 1
            elif spinMultiplicity == 1:
                radicalElectrons = 2; spinMultiplicity = 1
            elif spinMultiplicity == 2:
                radicalElectrons = 1; spinMultiplicity = 2
            elif spinMultiplicity == 3:
                radicalElectrons = 2; spinMultiplicity = 3
            elif spinMultiplicity == 4:
                radicalElectrons = 3; spinMultiplicity = 4
            elif spinMultiplicity == 5:
                radicalElectrons = 4; spinMultiplicity = 5
            
            # Process charge
            charge = obatom.GetFormalCharge()

            atom = Atom(element, radicalElectrons, spinMultiplicity, charge)
            self.vertices.append(atom)
            self.edges[atom] = {}
            
            # Add bonds by iterating again through atoms
            for j in range(0, i):
                obatom2 = obmol.GetAtom(j + 1)
                obbond = obatom.GetBond(obatom2)
                if obbond is not None:
                    order = 0

                    # Process bond type
                    if obbond.IsSingle(): order = 'S'
                    elif obbond.IsDouble(): order = 'D'
                    elif obbond.IsTriple(): order = 'T'
                    elif obbond.IsAromatic(): order = 'B'

                    bond = Bond(order)
                    atom1 = self.vertices[i]
                    atom2 = self.vertices[j]
                    self.edges[atom1][atom2] = bond
                    self.edges[atom2][atom1] = bond

        # Set atom types and connectivity values
        self.updateConnectivityValues()
        self.updateAtomTypes()

        return self

    def fromAdjacencyList(self, adjlist):
        """
        Convert a string adjacency list `adjlist` to a molecular structure.
        Skips the first line (assuming it's a label) unless `withLabel` is
        ``False``.
        """
        self.vertices, self.edges = fromAdjacencyList(adjlist, False, True)
        self.updateConnectivityValues()
        self.updateAtomTypes()
        return self

    def toCML(self):
        """
        Convert the molecular structure to CML. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        import pybel
        mol = pybel.Molecule(self.toOBMol())
        cml = mol.write('cml').strip()
        return '\n'.join([l for l in cml.split('\n') if l.strip()])

    def toInChI(self):
        """
        Convert a molecular structure to an InChI string. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        import openbabel
        # This version does not write a warning to stderr if stereochemistry is undefined
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
        
        Removes check-sum dash (-) and character so that only 
        the 14 + 9 characters remain.
        """
        import openbabel
        # This version does not write a warning to stderr if stereochemistry is undefined
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
        
    def toSMILES(self):
        """
        Convert a molecular structure to an SMILES string. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        import pybel
        mol = pybel.Molecule(self.toOBMol())
        return mol.write('smiles').strip()

    def toOBMol(self):
        """
        Convert a molecular structure to an OpenBabel OBMol object. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """

        import openbabel
        
        cython.declare(atom=Atom, atom1=Atom, bonds=dict, atom2=Atom, bond=Bond)
        cython.declare(index1=cython.int, index2=cython.int, order=cython.int)

        # Sort the atoms before converting to ensure output is consistent
        # between different runs
        self.sortAtoms()

        atoms = self.vertices
        bonds = self.edges

        obmol = openbabel.OBMol()
        for atom in atoms:
            a = obmol.NewAtom()
            a.SetAtomicNum(atom.number)
            a.SetFormalCharge(atom.charge)
        orders = {'S': 1, 'D': 2, 'T': 3, 'B': 5}
        for atom1, bonds in bonds.iteritems():
            for atom2, bond in bonds.iteritems():
                index1 = atoms.index(atom1)
                index2 = atoms.index(atom2)
                if index1 < index2:
                    order = orders[bond.order]
                    obmol.AddBond(index1+1, index2+1, order)

        obmol.AssignSpinMultiplicity(True)

        return obmol

    def toAdjacencyList(self, label='', removeH=False):
        """
        Convert the molecular structure to a string adjacency list.
        """
        result = toAdjacencyList(self, label=label, group=False, removeH=removeH)
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
        for atom1 in self.edges:
            for bond in self.edges[atom1].values():
                if not bond.isDouble(): allDoubleBonds = False
        if allDoubleBonds: return True

        # True if alternating single-triple bonds (e.g. H-C#C-H)
        # This test requires explicit hydrogen atoms
        for atom in self.vertices:
            bonds = self.edges[atom].values()
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

    def countInternalRotors(self):
        """
        Determine the number of internal rotors in the structure. Any single
        bond not in a cycle and between two atoms that also have other bonds
        are considered to be internal rotors.
        """
        count = 0
        for atom1 in self.edges:
            for atom2, bond in self.edges[atom1].iteritems():
                if self.vertices.index(atom1) < self.vertices.index(atom2) and bond.isSingle() and not self.isBondInCycle(atom1, atom2):
                    if len(self.edges[atom1]) > 1 and len(self.edges[atom2]) > 1:
                        count += 1
        return count

    def calculateSymmetryNumber(self):
        """
        Return the symmetry number for the structure. The symmetry number
        includes both external and internal modes.
        """
        from rmgpy.molecule.symmetry import calculateSymmetryNumber
        self.symmetryNumber = calculateSymmetryNumber(self)
        return self.symmetryNumber
    
    def generateResonanceIsomers(self):
        """
        Generate and return all of the resonance isomers of this molecule.
        """
        
        isomers = [self]

        # Radicals
        if sum([atom.radicalElectrons for atom in isomers[0].atoms]) > 0:
            # Iterate over resonance isomers
            index = 0
            while index < len(isomers):
                isomer = isomers[index]
                newIsomers = isomer.getAdjacentResonanceIsomers()
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

        isomers = []

        # Radicals
        if sum([atom.radicalElectrons for atom in self.vertices]) > 0:
            # Iterate over radicals in structure
            for atom in self.vertices:
                paths = self.findAllDelocalizationPaths(atom)
                for path in paths:
                    atom1, atom2, atom3, bond12, bond23 = path
                    # Adjust to (potentially) new resonance isomer
                    atom1.decrementRadical()
                    atom3.incrementRadical()
                    bond12.incrementOrder()
                    bond23.decrementOrder()
                    # Make a copy of isomer
                    isomer = self.copy(deep=True)
                    # Also copy the connectivity values, since they are the same
                    # for all resonance forms
                    for v1, v2 in zip(self.vertices, isomer.vertices):
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

    def findAllDelocalizationPaths(self, atom1):
        """
        Find all the delocalization paths allyl to the radical center indicated
        by `atom1`. Used to generate resonance isomers.
        """

        # No paths if atom1 is not a radical
        if atom1.radicalElectrons <= 0:
            return []

        # Find all delocalization paths
        paths = []
        for atom2, bond12 in self.edges[atom1].iteritems():
            # Vinyl bond must be capable of gaining an order
            if bond12.order in ['S', 'D']:
                for atom3, bond23 in self.getBonds(atom2).iteritems():
                    # Allyl bond must be capable of losing an order without breaking
                    if atom1 is not atom3 and bond23.order in ['D', 'T']:
                        paths.append([atom1, atom2, atom3, bond12, bond23])
        return paths

    def getURL(self):
        """
        Get a URL to the molecule's info page on the RMG website.
        """
        # eg. http://dev.rmg.mit.edu/database/kinetics/reaction/reactant1=1%20C%200%20%7B2,S%7D;2%20O%200%20%7B1,S%7D;__reactant2=1%20C%202T;__product1=1%20C%201;__product2=1%20C%200%20%7B2,S%7D;2%20O%201%20%7B1,S%7D;

        url = "http://rmg.mit.edu/database/molecule/"
        adjlist = self.toAdjacencyList(removeH=True)
        url += "{0}".format(re.sub('\s+', '%20', adjlist.replace('\n', ';')))
        return url.strip('_')
        