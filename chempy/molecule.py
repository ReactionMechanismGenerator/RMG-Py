#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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

import element as elements
from graph import Vertex, Edge, Graph
from exception import ChemPyError
from pattern import AtomPattern, BondPattern, MoleculePattern, AtomType
from pattern import getAtomType, fromAdjacencyList, toAdjacencyList

################################################################################

class Atom(Vertex):
    """
    An atom. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `element`           :class:`Element`    The chemical element the atom represents
    `radicalElectrons`  ``short``           The number of radical electrons
    `spinMultiplicity`  ``short``           The spin multiplicity of the atom
    `implicitHydrogens` ``short``           The number of implicit hydrogen atoms bonded to this atom
    `charge`            ``short``           The formal charge of the atom
    `label`             ``str``             A string label that can be used to tag individual atoms
    =================== =================== ====================================

    Additionally, the ``mass``, ``number``, and ``symbol`` attributes of the
    atom's element can be read (but not written) directly from the atom object,
    e.g. ``atom.symbol`` instead of ``atom.element.symbol``.
    """

    def __init__(self, element=None, radicalElectrons=0, spinMultiplicity=1, implicitHydrogens=0, charge=0, label=''):
        Vertex.__init__(self)
        if isinstance(element, str):
            self.element = elements.__dict__[element]
        else:
            self.element = element
        self.radicalElectrons = radicalElectrons
        self.spinMultiplicity = spinMultiplicity
        self.implicitHydrogens = implicitHydrogens
        self.charge = charge
        self.label = label
        self.atomType = None

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return "<Atom '%s'>" % (
            str(self.element) +
            ''.join(['.' for i in range(self.radicalElectrons)]) +
            ''.join(['+' for i in range(self.charge)]) +
            ''.join(['-' for i in range(-self.charge)])
        )

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "Atom(element='%s', radicalElectrons=%s, spinMultiplicity=%s, implicitHydrogens=%s, charge=%s, label='%s')" % (self.element, self.radicalElectrons, self.spinMultiplicity, self.implicitHydrogens, self.charge, self.label)

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
        :class:`AtomPattern` object, then the atom must match any of the
        combinations in the atom pattern.
        """
        cython.declare(atom=Atom, ap=AtomPattern)
        if isinstance(other, Atom):
            atom = other
            return (self.element is atom.element and
                self.radicalElectrons == atom.radicalElectrons and
                self.spinMultiplicity == atom.spinMultiplicity and
                self.implicitHydrogens == atom.implicitHydrogens and
                self.charge == atom.charge)
        elif isinstance(other, AtomPattern):
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
        :class:`AtomPattern` object, then the atom must match or be more
        specific than any of the combinations in the atom pattern.
        """
        if isinstance(other, Atom):
            return self.equivalent(other)
        elif isinstance(other, AtomPattern):
            cython.declare(atom=AtomPattern, a=AtomType, radical=cython.short, spin=cython.short, charge=cython.short)
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
        a = Atom(self.element, self.radicalElectrons, self.spinMultiplicity, self.implicitHydrogens, self.charge, self.label)
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
        self.spinMultiplicity += 1

    def decrementRadical(self):
        """
        Update the atom pattern as a result of applying a LOSE_RADICAL action,
        where `radical` specifies the number of radical electrons to remove.
        """
        # Set the new radical electron counts and spin multiplicities
        if self.radicalElectrons - 1 < 0:
            raise ChemPyError('Unable to update Atom due to LOSE_RADICAL action: Invalid radical electron set "%s".' % (self.radicalElectrons))
        self.radicalElectrons -= 1
        if self.spinMultiplicity - 1 < 0:
            self.spinMultiplicity -= 1 - 2
        else:
            self.spinMultiplicity -= 1

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
            raise ChemPyError('Unable to update Atom: Invalid action %s".' % (action))

################################################################################

class Bond(Edge):
    """
    A chemical bond. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `order`             ``str``             The bond order (``S`` = single, `D`` = double, ``T`` = triple, ``B`` = benzene)
    =================== =================== ====================================

    """

    def __init__(self, order=1):
        Edge.__init__(self)
        self.order = order

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return "<Bond '%s'>" % (self.order)

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "Bond(order='%s')" % (self.order)

    def equivalent(self, other):
        """
        Return ``True`` if `other` is indistinguishable from this bond, or
        ``False`` otherwise. `other` can be either a :class:`Bond` or a
        :class:`BondPattern` object.
        """
        cython.declare(bond=Bond, bp=BondPattern)
        if isinstance(other, Bond):
            bond = other
            return (self.order == bond.order)
        elif isinstance(other, BondPattern):
            bp = other
            return (self.order in bp.order)

    def isSpecificCaseOf(self, other):
        """
        Return ``True`` if `self` is a specific case of `other`, or ``False``
        otherwise. `other` can be either a :class:`Bond` or a
        :class:`BondPattern` object.
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
            raise ChemPyError('Unable to update Bond due to CHANGE_BOND action: Invalid bond order "%s".' % (self.order))
        
    def decrementOrder(self):
        """
        Update the bond as a result of applying a CHANGE_BOND action to
        decrease the order by one.
        """
        if self.order == 'D': self.order = 'S'
        elif self.order == 'T': self.order = 'D'
        else:
            raise ChemPyError('Unable to update Bond due to CHANGE_BOND action: Invalid bond order "%s".' % (self.order))
        
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
                raise ChemPyError('Unable to update Bond due to CHANGE_BOND action: Invalid bond order "%s".' % (self.order))
        elif order == -1:
            if self.order == 'D': self.order = 'S'
            elif self.order == 'T': self.order = 'D'
            else:
                raise ChemPyError('Unable to update Bond due to CHANGE_BOND action: Invalid bond order "%s".' % (self.order))
        else:
            raise ChemPyError('Unable to update Bond due to CHANGE_BOND action: Invalid order "%g".' % order)

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
                raise ChemPyError('Unable to update Bond due to CHANGE_BOND action: Invalid order "%g".' % action[2])
        else:
            raise ChemPyError('Unable to update BondPattern: Invalid action %s".' % (action))

################################################################################

class Molecule(Graph):
    """
    A representation of a molecular structure using a graph data type, extending
    the :class:`Graph` class. The `atoms` and `bonds` attributes are aliases
    for the `vertices` and `edges` attributes. Corresponding alias methods have
    also been provided.
    """

    def __init__(self, atoms=None, bonds=None, SMILES='', InChI=''):
        Graph.__init__(self, atoms, bonds)
        self.implicitHydrogens = False
        if SMILES != '': self.fromSMILES(SMILES)
        elif InChI != '': self.fromInChI(InChI)
    
    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return "<Molecule '%s'>" % (self.toSMILES())

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "Molecule(SMILES='%s')" % (self.toSMILES())

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
        return sum([atom.mass for atom in self.vertices])

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

    def makeHydrogensImplicit(self):
        """
        Convert all explicitly stored hydrogen atoms to be stored implicitly.
        An implicit hydrogen atom is stored on the heavy atom it is connected
        to as a single integer counter. This is done to save memory.
        """

        cython.declare(atom=Atom, neighbor=Atom, hydrogens=list)

        # Check that the structure contains at least one heavy atom
        for atom in self.vertices:
            if not atom.isHydrogen():
                break
        else:
            # No heavy atoms, so leave explicit
            return
        
        # Count the hydrogen atoms on each non-hydrogen atom and set the
        # `implicitHydrogens` attribute accordingly
        hydrogens = []
        for atom in self.vertices:
            if atom.isHydrogen():
                neighbor = self.edges[atom].keys()[0]
                neighbor.implicitHydrogens += 1
                hydrogens.append(atom)

        # Remove the hydrogen atoms from the structure
        for atom in hydrogens:
            self.removeAtom(atom)

        # Set implicitHydrogens flag to True
        self.implicitHydrogens = True

    def makeHydrogensExplicit(self):
        """
        Convert all implicitly stored hydrogen atoms to be stored explicitly.
        An explicit hydrogen atom is stored as its own atom in the graph, with
        a single bond to the heavy atom it is attached to. This consumes more
        memory, but may be required for certain tasks (e.g. subgraph matching).
        """

        cython.declare(atom=Atom, H=Atom, bond=Bond, hydrogens=list, numAtoms=cython.short)

        # Create new hydrogen atoms for each implicit hydrogen
        hydrogens = []
        for atom in self.vertices:
            while atom.implicitHydrogens > 0:
                H = Atom(element='H')
                bond = Bond(order='S')
                hydrogens.append((H, atom, bond))
                atom.implicitHydrogens -= 1

        # Add the hydrogens to the graph
        numAtoms = len(self.vertices)
        for H, atom, bond in hydrogens:
            self.addAtom(H)
            self.addBond(H, atom, bond)
            H.atomType = getAtomType(H, {atom:bond})
            # If known, set the connectivity information
            H.connectivity1 = 1
            H.connectivity2 = atom.connectivity1
            H.connectivity3 = atom.connectivity2
            H.sortingLabel = numAtoms
            numAtoms += 1

        # Set implicitHydrogens flag to False
        self.implicitHydrogens = False

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
        return None

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
            raise TypeError('Got a %s object for parameter "other", when a Molecule object is required.' % other.__class__)
        # Ensure that both self and other have the same implicit hydrogen status
        # If not, make them both explicit just to be safe
        implicitH = [self.implicitHydrogens, other.implicitHydrogens]
        if not all(implicitH):
            self.makeHydrogensExplicit()
            other.makeHydrogensExplicit()
        # Do the isomorphism comparison
        result = Graph.isIsomorphic(self, other, initialMap)
        # Restore implicit status if needed
        if implicitH[0]: self.makeHydrogensImplicit()
        if implicitH[1]: other.makeHydrogensImplicit()
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
            raise TypeError('Got a %s object for parameter "other", when a Molecule object is required.' % other.__class__)
        # Ensure that both self and other have the same implicit hydrogen status
        # If not, make them both explicit just to be safe
        implicitH = [self.implicitHydrogens, other.implicitHydrogens]
        if not all(implicitH):
            self.makeHydrogensExplicit()
            other.makeHydrogensExplicit()
        # Do the isomorphism comparison
        result = Graph.findIsomorphism(self, other, initialMap)
        # Restore implicit status if needed
        if implicitH[0]: self.makeHydrogensImplicit()
        if implicitH[1]: other.makeHydrogensImplicit()
        return result

    def isSubgraphIsomorphic(self, other, initialMap=None):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise. The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`MoleculePattern` object, or a :class:`TypeError` is raised.
        """
        # It only makes sense to compare a Molecule to a MoleculePattern for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, MoleculePattern):
            raise TypeError('Got a %s object for parameter "other", when a MoleculePattern object is required.' % other.__class__)
        # Ensure that self is explicit (assume other is explicit)
        implicitH = self.implicitHydrogens
        self.makeHydrogensExplicit()
        # Do the isomorphism comparison
        result = Graph.isSubgraphIsomorphic(self, other, initialMap)
        # Restore implicit status if needed
        if implicitH: self.makeHydrogensImplicit()
        return result

    def findSubgraphIsomorphisms(self, other, initialMap=None):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise. Also returns the lists all of valid mappings. The
        `initialMap` attribute can be used to specify a required mapping from
        `self` to `other` (i.e. the atoms of `self` are the keys, while the
        atoms of `other` are the values). The returned mappings also use the
        atoms of `self` for the keys and the atoms of `other` for the values.
        The `other` parameter must be a :class:`MoleculePattern` object, or a
        :class:`TypeError` is raised.
        """
        # It only makes sense to compare a Molecule to a MoleculePattern for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, MoleculePattern):
            raise TypeError('Got a %s object for parameter "other", when a MoleculePattern object is required.' % other.__class__)
        # Ensure that self is explicit (assume other is explicit)
        implicitH = self.implicitHydrogens
        self.makeHydrogensExplicit()
        # Do the isomorphism comparison
        result = Graph.findSubgraphIsomorphisms(self, other, initialMap)
        # Restore implicit status if needed
        if implicitH: self.makeHydrogensImplicit()
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
        from ext.molecule_draw import drawMolecule
        drawMolecule(self, path=path)

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
        import pybel
        mol = pybel.readstring('inchi', inchistr)
        self.fromOBMol(mol.OBMol)
        return self

    def fromSMILES(self, smilesstr):
        """
        Convert a SMILES string `smilesstr` to a molecular structure. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
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

        # Add hydrogen atoms to complete molecule if needed
        obmol.AddHydrogens()

        # Iterate through atoms in obmol
        for i in range(0, obmol.NumAtoms()):
            obatom = obmol.GetAtom(i + 1)

            # Use atomic number as key for element
            number = obatom.GetAtomicNum()
            element = elements.getElement(number=number)
            
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

            # Process charge
            charge = obatom.GetFormalCharge()

            atom = Atom(element, radicalElectrons, spinMultiplicity, 0, charge)
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

        # Make hydrogens implicit to conserve memory
        self.makeHydrogensImplicit()

        return self

    def fromAdjacencyList(self, adjlist, withLabel=True):
        """
        Convert a string adjacency list `adjlist` to a molecular structure.
        Skips the first line (assuming it's a label) unless `withLabel` is
        ``False``.
        """
        self.vertices, self.edges = fromAdjacencyList(adjlist, False, True, withLabel)
        self.updateConnectivityValues()
        self.updateAtomTypes()
        self.makeHydrogensImplicit()
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
        
        cython.declare(implicitH=cython.bint)
        cython.declare(atom=Atom, atom1=Atom, bonds=dict, atom2=Atom, bond=Bond)
        cython.declare(index1=cython.int, index2=cython.int, order=cython.int)

        # Make hydrogens explicit while we perform the conversion
        implicitH = self.implicitHydrogens
        self.makeHydrogensExplicit()

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

        # Restore implicit hydrogens if necessary
        if implicitH: self.makeHydrogensImplicit()

        return obmol

    def toAdjacencyList(self):
        """
        Convert the molecular structure to a string adjacency list.
        """
        return toAdjacencyList(self)

    def isLinear(self):
        """
        Return :data:`True` if the structure is linear and :data:`False`
        otherwise.
        """

        atomCount = len(self.vertices) + sum([atom.implicitHydrogens for atom in self.vertices])

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
            if atom1.implicitHydrogens > 0: allDoubleBonds = False
            for bond in self.edges[atom1].values():
                if not bond.isDouble(): allDoubleBonds = False
        if allDoubleBonds: return True

        # True if alternating single-triple bonds (e.g. H-C#C-H)
        # This test requires explicit hydrogen atoms
        implicitH = self.implicitHydrogens
        self.makeHydrogensExplicit()
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
            if implicitH: self.makeHydrogensImplicit()
            return True
        
        # not returned yet? must be nonlinear
        if implicitH: self.makeHydrogensImplicit()
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
                    if len(self.edges[atom1]) + atom1.implicitHydrogens > 1 and len(self.edges[atom2]) + atom2.implicitHydrogens > 1:
                        count += 1
        return count

    def calculateAtomSymmetryNumber(self, atom):
        """
        Return the symmetry number centered at `atom` in the structure. The
        `atom` of interest must not be in a cycle.
        """
        symmetryNumber = 1

        single = 0; double = 0; triple = 0; benzene = 0
        numNeighbors = 0
        for bond in self.edges[atom].values():
            if bond.isSingle(): single += 1
            elif bond.isDouble(): double += 1
            elif bond.isTriple(): triple += 1
            elif bond.isBenzene(): benzene += 1
            numNeighbors += 1
        
        # If atom has zero or one neighbors, the symmetry number is 1
        if numNeighbors < 2: return symmetryNumber

        # Create temporary structures for each functional group attached to atom
        molecule = self.copy()
        for atom2 in molecule.bonds[atom].keys(): molecule.removeBond(atom, atom2)
        molecule.removeAtom(atom)
        groups = molecule.split()

        # Determine equivalence of functional groups around atom
        groupIsomorphism = dict([(group, dict()) for group in groups])
        for group1 in groups:
            for group2 in groups:
                if group1 is not group2 and group2 not in groupIsomorphism[group1]:
                    groupIsomorphism[group1][group2] = group1.isIsomorphic(group2)
                    groupIsomorphism[group2][group1] = groupIsomorphism[group1][group2]
                elif group1 is group2:
                    groupIsomorphism[group1][group1] = True
        count = [sum([int(groupIsomorphism[group1][group2]) for group2 in groups]) for group1 in groups]
        for i in range(count.count(2) / 2):
            count.remove(2)
        for i in range(count.count(3) / 3):
            count.remove(3); count.remove(3)
        for i in range(count.count(4) / 4):
            count.remove(4); count.remove(4); count.remove(4)
        count.sort(); count.reverse()
        
        if atom.radicalElectrons == 0:
            if single == 4:
                # Four single bonds
                if count == [4]: symmetryNumber *= 12
                elif count == [3, 1]: symmetryNumber *= 3
                elif count == [2, 2]: symmetryNumber *= 2
                elif count == [2, 1, 1]: symmetryNumber *= 1
                elif count == [1, 1, 1, 1]: symmetryNumber *= 1
            elif single == 2:
                # Two single bonds
                if count == [2]: symmetryNumber *= 2
            elif double == 2:
                # Two double bonds
                if count == [2]: symmetryNumber *= 2
        elif atom.radicalElectrons == 1:
            if single == 3:
                # Three single bonds
                if count == [3]: symmetryNumber *= 6
                elif count == [2, 1]: symmetryNumber *= 2
                elif count == [1, 1, 1]: symmetryNumber *= 1
        elif atom.radicalElectrons == 2:
            if single == 2:
                # Two single bonds
                if count == [2]: symmetryNumber *= 2

        return symmetryNumber

    def calculateBondSymmetryNumber(self, atom1, atom2):
        """
        Return the symmetry number centered at `bond` in the structure.
        """
        bond = self.edges[atom1][atom2]
        symmetryNumber = 1
        if bond.isSingle() or bond.isDouble() or bond.isTriple():
            if atom1.equivalent(atom2):
                # An O-O bond is considered to be an "optical isomer" and so no
                # symmetry correction will be applied
                if atom1.atomType == atom2.atomType == 'Os' and \
                    atom1.radicalElectrons == atom2.radicalElectrons == 0:
                    pass
                # If the molecule is diatomic, then we don't have to check the
                # ligands on the two atoms in this bond (since we know there
                # aren't any)
                elif len(self.vertices) == 2:
                    symmetryNumber = 2
                else:
                    molecule = self.copy()
                    molecule.removeBond(atom1, atom2)
                    fragments = molecule.split()
                    if len(fragments) != 2: return symmetryNumber

                    fragment1, fragment2 = fragments
                    if atom1 in fragment1.atoms: fragment1.removeAtom(atom1)
                    if atom2 in fragment1.atoms: fragment1.removeAtom(atom2)
                    if atom1 in fragment2.atoms: fragment2.removeAtom(atom1)
                    if atom2 in fragment2.atoms: fragment2.removeAtom(atom2)
                    groups1 = fragment1.split()
                    groups2 = fragment2.split()

                    # Test functional groups for symmetry
                    if len(groups1) == len(groups2) == 1:
                        if groups1[0].isIsomorphic(groups2[0]): symmetryNumber *= 2
                    elif len(groups1) == len(groups2) == 2:
                        if groups1[0].isIsomorphic(groups2[0]) and groups1[1].isIsomorphic(groups2[1]): symmetryNumber *= 2
                        elif groups1[1].isIsomorphic(groups2[0]) and groups1[0].isIsomorphic(groups2[1]): symmetryNumber *= 2
                    elif len(groups1) == len(groups2) == 3:
                        if groups1[0].isIsomorphic(groups2[0]) and groups1[1].isIsomorphic(groups2[1]) and groups1[2].isIsomorphic(groups2[2]): symmetryNumber *= 2
                        elif groups1[0].isIsomorphic(groups2[0]) and groups1[1].isIsomorphic(groups2[2]) and groups1[2].isIsomorphic(groups2[1]): symmetryNumber *= 2
                        elif groups1[0].isIsomorphic(groups2[1]) and groups1[1].isIsomorphic(groups2[2]) and groups1[2].isIsomorphic(groups2[0]): symmetryNumber *= 2
                        elif groups1[0].isIsomorphic(groups2[1]) and groups1[1].isIsomorphic(groups2[0]) and groups1[2].isIsomorphic(groups2[2]): symmetryNumber *= 2
                        elif groups1[0].isIsomorphic(groups2[2]) and groups1[1].isIsomorphic(groups2[0]) and groups1[2].isIsomorphic(groups2[1]): symmetryNumber *= 2
                        elif groups1[0].isIsomorphic(groups2[2]) and groups1[1].isIsomorphic(groups2[1]) and groups1[2].isIsomorphic(groups2[0]): symmetryNumber *= 2

        return symmetryNumber

    def calculateAxisSymmetryNumber(self):
        """
        Get the axis symmetry number correction. The "axis" refers to a series
        of two or more cumulated double bonds (e.g. C=C=C, etc.). Corrections
        for single C=C bonds are handled in getBondSymmetryNumber().
        
        Each axis (C=C=C) has the potential to double the symmetry number.
        If an end has 0 or 1 groups (eg. =C=CJJ or =C=C-R) then it cannot 
        alter the axis symmetry and is disregarded::
        
            A=C=C=C..        A-C=C=C=C-A
            
              s=1                s=1
        
        If an end has 2 groups that are different then it breaks the symmetry 
        and the symmetry for that axis is 1, no matter what's at the other end::
        
            A\               A\         /A
              T=C=C=C=C-A      T=C=C=C=T
            B/               A/         \B
                  s=1             s=1
        
        If you have one or more ends with 2 groups, and neither end breaks the 
        symmetry, then you have an axis symmetry number of 2::
        
            A\         /B      A\         
              C=C=C=C=C          C=C=C=C-B
            A/         \B      A/         
                  s=2                s=2
        """

        symmetryNumber = 1

        # List all double bonds in the structure
        doubleBonds = []
        for atom1 in self.edges:
            for atom2 in self.edges[atom1]:
                if self.edges[atom1][atom2].isDouble() and self.vertices.index(atom1) < self.vertices.index(atom2):
                    doubleBonds.append((atom1, atom2))

        # Search for adjacent double bonds
        cumulatedBonds = []
        for i, bond1 in enumerate(doubleBonds):
            atom11, atom12 = bond1
            for bond2 in doubleBonds[i+1:]:
                atom21, atom22 = bond2
                if atom11 is atom21 or atom11 is atom22 or atom12 is atom21 or atom12 is atom22:
                    listToAddTo = None
                    for cumBonds in cumulatedBonds:
                        if (atom11, atom12) in cumBonds or (atom21, atom22) in cumBonds:
                            listToAddTo = cumBonds
                    if listToAddTo is not None:
                        if (atom11, atom12) not in listToAddTo: listToAddTo.append((atom11, atom12))
                        if (atom21, atom22) not in listToAddTo: listToAddTo.append((atom21, atom22))
                    else:
                        cumulatedBonds.append([(atom11, atom12), (atom21, atom22)])

        # For each set of adjacent double bonds, check for axis symmetry
        for bonds in cumulatedBonds:
            
            # Do nothing if less than two cumulated bonds
            if len(bonds) < 2: continue

            # Do nothing if axis is in cycle
            found = False
            for atom1, atom2 in bonds:
               if self.isBondInCycle(atom1, atom2): found = True
            if found: continue

            # Find terminal atoms in axis
            # Terminal atoms labelled T:  T=C=C=C=T
            axis = []
            for bond in bonds: axis.extend(bond)
            terminalAtoms = []
            for atom in axis:
                if axis.count(atom) == 1: terminalAtoms.append(atom)
            if len(terminalAtoms) != 2: continue
            
            # Remove axis from (copy of) structure
            structure = self.copy()
            for atom1, atom2 in bonds:
                structure.removeBond(atom1, atom2)
            atomsToRemove = []
            for atom in structure.atoms:
                if len(structure.bonds[atom]) == 0: # it's not bonded to anything
                    atomsToRemove.append(atom)
            for atom in atomsToRemove: structure.removeAtom(atom)

            # Split remaining fragments of structure
            end_fragments = structure.split()
            # you may have only one end fragment,
            # eg. if you started with H2C=C=C..
            
            # 
            # there can be two groups at each end     A\         /B
            #                                           T=C=C=C=T
            #                                         A/         \B
            
            # to start with nothing has broken symmetry about the axis
            symmetry_broken=False 
            for fragment in end_fragments: # a fragment is one end of the axis
                
                # remove the atom that was at the end of the axis and split what's left into groups
                for atom in terminalAtoms:
                    if atom in fragment.atoms: fragment.removeAtom(atom)
                groups = fragment.split()
                
                # If end has only one group then it can't contribute to (nor break) axial symmetry
                #   Eg. this has no axis symmetry:   A-T=C=C=C=T-A
                # so we remove this end from the list of interesting end fragments
                if len(groups)==1:
                    end_fragments.remove(fragment)
                    continue # next end fragment
                if len(groups)==2:
                    if not groups[0].isIsomorphic(groups[1]):
                        # this end has broken the symmetry of the axis
                        symmetry_broken = True
                        
            # If there are end fragments left that can contribute to symmetry,
            # and none of them broke it, then double the symmetry number
            # NB>> This assumes coordination number of 4 (eg. Carbon).
            #      And would be wrong if we had /B
            #                         =C=C=C=C=T-B
            #                                   \B
            #      (for some T with coordination number 5).
            if end_fragments and not symmetry_broken:
                symmetryNumber *= 2
                    
        return symmetryNumber

    def calculateCyclicSymmetryNumber(self):
        """
        Get the symmetry number correction for cyclic regions of a molecule.
        For complicated fused rings the smallest set of smallest rings is used.
        """

        symmetryNumber = 1

        # Get symmetry number for each ring in structure
        rings = self.getSmallestSetOfSmallestRings()
        for ring in rings:

            # Make copy of structure
            structure = self.copy()

            # Remove bonds of ring from structure
            for i, atom1 in enumerate(ring):
                for atom2 in ring[i+1:]:
                    if structure.hasBond(atom1, atom2):
                        structure.removeBond(atom1, atom2)

            structures = structure.split()
            groups = []
            for struct in structures:
                for atom in ring:
                    if atom in struct.atoms(): struct.removeAtom(atom)
                groups.append(struct.split())

            # Find equivalent functional groups on ring
            equivalentGroups = []
            for group in groups:
                found = False
                for eqGroup in equivalentGroups:
                    if not found:
                        if group.isIsomorphic(eqGroup[0]):
                            eqGroup.append(group)
                            found = True
                if not found:
                    equivalentGroups.append([group])

            # Find equivalent bonds on ring
            equivalentBonds = []
            for i, atom1 in enumerate(ring):
                for atom2 in ring[i+1:]:
                    if self.hasBond(atom1, atom2):
                        bond = self.getBond(atom1, atom2)
                        found = False
                        for eqBond in equivalentBonds:
                            if not found:
                                if bond.equivalent(eqBond[0]):
                                    eqBond.append(group)
                                    found = True
                        if not found:
                            equivalentBonds.append([bond])

            # Find maximum number of equivalent groups and bonds
            maxEquivalentGroups = 0
            for groups in equivalentGroups:
                if len(groups) > maxEquivalentGroups:
                    maxEquivalentGroups = len(groups)
            maxEquivalentBonds = 0
            for bonds in equivalentBonds:
                if len(bonds) > maxEquivalentBonds:
                    maxEquivalentBonds = len(bonds)

            if maxEquivalentGroups == maxEquivalentBonds == len(ring):
                symmetryNumber *= len(ring)
            else:
                symmetryNumber *= max(maxEquivalentGroups, maxEquivalentBonds)

            print len(ring), maxEquivalentGroups, maxEquivalentBonds, symmetryNumber


        return symmetryNumber

    def calculateSymmetryNumber(self):
        """
        Return the symmetry number for the structure. The symmetry number
        includes both external and internal modes.
        """
        symmetryNumber = 1

        implicitH = self.implicitHydrogens
        self.makeHydrogensExplicit()

        for atom in self.vertices:
            if not self.isAtomInCycle(atom):
                symmetryNumber *= self.calculateAtomSymmetryNumber(atom)

        for atom1 in self.edges:
            for atom2 in self.edges[atom1]:
                if self.vertices.index(atom1) < self.vertices.index(atom2) and not self.isBondInCycle(atom1, atom2):
                    symmetryNumber *= self.calculateBondSymmetryNumber(atom1, atom2)

        symmetryNumber *= self.calculateAxisSymmetryNumber()

        #if self.isCyclic():
        #   symmetryNumber *= self.calculateCyclicSymmetryNumber()

        self.symmetryNumber = symmetryNumber

        if implicitH: self.makeHydrogensImplicit()

        return symmetryNumber

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
