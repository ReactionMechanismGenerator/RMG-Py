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
molecular configurations.
"""

import cython

import element as elements
import graph
from exception import ChemPyError
from pattern import *

################################################################################

class Atom(graph.Vertex):
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
        graph.Vertex.__init__(self)
        if isinstance(element, str):
            self.element = elements.__dict__[element]
        else:
            self.element = element
        self.radicalElectrons = radicalElectrons
        self.spinMultiplicity = spinMultiplicity
        self.implicitHydrogens = implicitHydrogens
        self.charge = charge
        self.label = label
        self.atomType = ''

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
        ``False`` otherwise.
        """
        if isinstance(other, Atom):
            cython.declare(atom=Atom)
            atom = other
            return (self.element is atom.element and
                self.radicalElectrons == atom.radicalElectrons and
                self.spinMultiplicity == atom.spinMultiplicity and
                self.implicitHydrogens == atom.implicitHydrogens and
                self.charge == atom.charge)
        elif isinstance(other, AtomPattern):
            cython.declare(atom=AtomPattern)
            atom = other
            return (any([atomTypesEquivalent(self.atomType, a) for a in atom.atomType]) and
                [self.radicalElectrons, self.spinMultiplicity] in zip(atom.radicalElectrons, atom.spinMultiplicity) and
                self.charge in atom.charge)

    def isSpecificCaseOf(self, other):
        """
        Return ``True`` if `self` is a specific case of `other`, or ``False``
        otherwise.
        """
        if isinstance(other, Atom):
            return self.equivalent(other)
        elif isinstance(other, AtomPattern):
            cython.declare(atom=AtomPattern)
            atom = other
            return (any([atomTypesSpecificCaseOf(self.atomType, a) for a in atom.atomType]) and
                (self.radicalElectrons, self.spinMultiplicity) in zip(atom.radicalElectrons, atom.spinMultiplicity) and
                self.charge in atom.charge)

    def copy(self):
        """
        Generate a copy of the current atom.
        """
        return Atom(self.element, self.radicalElectrons, self.spinMultiplicity, self.implicitHydrogens, self.charge, self.label)

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

################################################################################

class Bond(graph.Edge):
    """
    A chemical bond. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `order`             ``str``             The bond order ('S' = single, 'D' = double, 'T' = triple, 'B' = benzene)
    =================== =================== ====================================

    """

    def __init__(self, order=1):
        graph.Edge.__init__(self)
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
        ``False`` otherwise.
        """
        if isinstance(other, Bond):
            cython.declare(bond=Bond)
            bond = other
            return (self.order == bond.order)
        elif isinstance(other, BondPattern):
            cython.declare(bond=BondPattern)
            bond = other
            return (self.order in bond.order)

    def isSpecificCaseOf(self, other):
        """
        Return ``True`` if `self` is a specific case of `other`, or ``False``
        otherwise.
        """
        # There are no generic bond types, so isSpecificCaseOf is the same as equivalent
        return self.equivalent(other)

    def copy(self):
        """
        Generate a copy of the current bond.
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

################################################################################

class Molecule(graph.Graph):
    """
    A representation of a molecular structure using a graph data type, extending
    the :class:`Graph` class. The `atoms` and `bonds` attributes are aliases
    for the `vertices` and `edges` attributes. Corresponding alias methods have
    also been provided.
    """

    def __init__(self, atoms=None, bonds=None, SMILES='', InChI=''):
        graph.Graph.__init__(self, atoms, bonds)
        if SMILES != '': self.fromSMILES(SMILES)
        elif InChI != '': self.fromInChI(InChI)
        self.implicitHydrogens = False

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return "<Molecule '%s'>" % (self.toSMILES())

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "Bond(SMILES='%s')" % (self.toSMILES())


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

    def copy(self, deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """
        other = cython.declare(Molecule)
        g = graph.Graph.copy(self, deep)
        other = Molecule(g.vertices, g.edges)
        return other

    def makeHydrogensImplicit(self):
        """
        Convert all explicitly stored hydrogen atoms to be stored implicitly.
        An implicit hydrogen atom is stored on the heavy atom it is connected
        to as a single integer counter. This is done to save memory.
        """

        cython.declare(atom=Atom, neighbor=Atom, hydrogens=list)

        # Check that the structure contains at least one heavy atom
        if all([atom.isHydrogen() for atom in self.atoms]):
            return
        
        # Count the hydrogen atoms on each non-hydrogen atom and set the
        # `implicitHydrogens` attribute accordingly
        hydrogens = []
        for atom in self.atoms:
            if atom.isHydrogen():
                neighbor = self.bonds[atom].keys()[0]
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

        cython.declare(atom=Atom, H=Atom, bond=Bond, hydrogens=list)

        # Create new hydrogen atoms for each implicit hydrogen
        hydrogens = []
        for atom in self.atoms:
            while atom.implicitHydrogens > 0:
                H = Atom(element='H')
                H.atomType = 'H'
                bond = Bond(order='S')
                hydrogens.append((H, atom, bond))
                atom.implicitHydrogens -= 1

        # Add the hydrogens to the graph
        for H, atom, bond in hydrogens:
            self.addAtom(H)
            self.addBond(H, atom, bond)

        # Set implicitHydrogens flag to False
        self.implicitHydrogens = False

    def updateAtomTypes(self):
        """
        Iterate through the atoms in the structure, checking their atom types
        to ensure they are correct (i.e. accurately describe their local bond
        environment) and complete (i.e. are as detailed as possible).
        """
        for atom in self.atoms:
            atom.atomType = getAtomType(atom, self.bonds[atom])

    def getLabeledAtoms(self):
        """
        Return the labeled atoms as a ``dict`` with the keys being the labels
        and the values the atoms themselves. If two or more atoms have the
        same label, the value is converted to a list of these atoms.
        """
        labeled = {}
        for atom in self.atoms:
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
        result = graph.Graph.isIsomorphic(self, other, initialMap)
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
        result = graph.Graph.findIsomorphism(self, other, initialMap)
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
        result = graph.Graph.isSubgraphIsomorphic(self, other, initialMap)
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
        result = graph.Graph.findSubgraphIsomorphisms(self, other, initialMap)
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
            self.atoms.append(atom)

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
                    atom1 = self.atoms[i]
                    atom2 = self.atoms[j]
                    if atom1 not in self.bonds: self.bonds[atom1] = {}
                    if atom2 not in self.bonds: self.bonds[atom2] = {}
                    self.bonds[atom1][atom2] = bond
                    self.bonds[atom2][atom1] = bond

        # Make hydrogens implicit to conserve memory
        self.makeHydrogensImplicit()

        # Set atom types and connectivity values
        self.updateConnectivityValues()
        self.updateAtomTypes()

        return self

    def fromAdjacencyList(self, adjlist, withLabel=True):
        """
        Convert a string adjacency list `adjlist` to a molecular structure.
        Skips the first line (assuming it's a label) unless `withLabel` is
        ``False``.
        """
        self.vertices, self.edges = fromAdjacencyList(adjlist, pattern=False, addH=True, withLabel=withLabel)
        self.makeHydrogensImplicit()
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

        atoms = self.atoms
        bonds = self.bonds

        obmol = openbabel.OBMol()
        for atom in atoms:
            a = obmol.NewAtom()
            a.SetAtomicNum(atom.number)
            a.SetFormalCharge(atom.charge)
        for atom1, bonds in bonds.iteritems():
            for atom2, bond in bonds.iteritems():
                index1 = atoms.index(atom1)
                index2 = atoms.index(atom2)
                if index1 < index2:
                    order = bond.order
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
