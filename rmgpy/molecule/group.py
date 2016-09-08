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
This module provides classes and methods for working with molecular substructure
groups. These enable molecules to be searched for common motifs (e.g.
reaction sites).
"""

import cython

from .graph import Vertex, Edge, Graph
from .atomtype import atomTypes, allElements, nonSpecifics, getFeatures

################################################################################

class ActionError(Exception):
    """
    An exception class for errors that occur while applying reaction recipe
    actions. Pass a string describing the circumstances that caused the
    exceptional behavior.
    """
    pass

################################################################################

class GroupAtom(Vertex):
    """
    An atom group. This class is based on the :class:`Atom` class, except that
    it uses :ref:`atom types <atom-types>` instead of elements, and all
    attributes are lists rather than individual values. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `atomType`          ``list``            The allowed atom types (as :class:`AtomType` objects)
    `radicalElectrons`  ``list``            The allowed numbers of radical electrons (as short integers)
    `charge`            ``list``            The allowed formal charges (as short integers)
    `label`             ``str``             A string label that can be used to tag individual atoms
    `lonePairs`         ``list``            The number of lone electron pairs
    =================== =================== ====================================

    Each list represents a logical OR construct, i.e. an atom will match the
    group if it matches *any* item in the list. However, the
    `radicalElectrons`, and `charge` attributes are linked
    such that an atom must match values from the same index in each of these in
    order to match.
    """

    def __init__(self, atomType=None, radicalElectrons=None, charge=None, label='', lonePairs=None):
        Vertex.__init__(self)
        self.atomType = atomType or []
        for index in range(len(self.atomType)):
            if isinstance(self.atomType[index], str):
                self.atomType[index] = atomTypes[self.atomType[index]]
        self.radicalElectrons = radicalElectrons or []
        self.charge = charge or []
        self.label = label
        self.lonePairs = lonePairs or []

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
        }
        atomType = self.atomType
        if atomType is not None:
            atomType = [a.label for a in atomType]
        return (GroupAtom, (atomType, self.radicalElectrons, self.charge, self.label, self.lonePairs), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling an object.
        """
        self.edges = d['edges']
        self.connectivity1 = d['connectivity1']
        self.connectivity2 = d['connectivity2']
        self.connectivity3 = d['connectivity3']
        self.sortingLabel = d['sortingLabel']

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return '[{0}]'.format(','.join([repr(a.label) for a in self.atomType]))

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "<GroupAtom {0!s}>".format(self)

    @property
    def bonds(self): return self.edges

    def copy(self):
        """
        Return a deep copy of the :class:`GroupAtom` object. Modifying the
        attributes of the copy will not affect the original.
        """
        return GroupAtom(self.atomType[:], self.radicalElectrons[:], self.charge[:], self.label)

    def __changeBond(self, order):
        """
        Update the atom group as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order, and should be 1 or -1.
        """
        atomType = []
        for atom in self.atomType:
            if order == 1:
                atomType.extend(atom.incrementBond)
            elif order == -1:
                atomType.extend(atom.decrementBond)
            else:
                raise ActionError('Unable to update GroupAtom due to CHANGE_BOND action: Invalid order "{0}".'.format(order))
        if len(atomType) == 0:
            raise ActionError('Unable to update GroupAtom due to CHANGE_BOND action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __formBond(self, order):
        """
        Update the atom group as a result of applying a FORM_BOND action,
        where `order` specifies the order of the forming bond, and should be
        'S' (since we only allow forming of single bonds).
        """
        if order != 'S':
            raise ActionError('Unable to update GroupAtom due to FORM_BOND action: Invalid order "{0}".'.format(order))
        atomType = []
        for atom in self.atomType:
            atomType.extend(atom.formBond)
        if len(atomType) == 0:
            raise ActionError('Unable to update GroupAtom due to FORM_BOND action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __breakBond(self, order):
        """
        Update the atom group as a result of applying a BREAK_BOND action,
        where `order` specifies the order of the breaking bond, and should be
        'S' (since we only allow breaking of single bonds).
        """
        if order != 'S':
            raise ActionError('Unable to update GroupAtom due to BREAK_BOND action: Invalid order "{0}".'.format(order))
        atomType = []
        for atom in self.atomType:
            atomType.extend(atom.breakBond)
        if len(atomType) == 0:
            raise ActionError('Unable to update GroupAtom due to BREAK_BOND action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __gainRadical(self, radical):
        """
        Update the atom group as a result of applying a GAIN_RADICAL action,
        where `radical` specifies the number of radical electrons to add.
        """
        radicalElectrons = []
        if any([len(atomType.incrementRadical) == 0 for atomType in self.atomType]):
            raise ActionError('Unable to update GroupAtom due to GAIN_RADICAL action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        for electron in self.radicalElectrons:
            radicalElectrons.append(electron + radical)
        # Set the new radical electron counts
        self.radicalElectrons = radicalElectrons

    def __loseRadical(self, radical):
        """
        Update the atom group as a result of applying a LOSE_RADICAL action,
        where `radical` specifies the number of radical electrons to remove.
        """
        radicalElectrons = []
        pairs = set()
        if any([len(atomType.decrementRadical) == 0 for atomType in self.atomType]):
            raise ActionError('Unable to update GroupAtom due to LOSE_RADICAL action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        for electron in self.radicalElectrons:
            electron = electron - radical
            if electron < 0:
                raise ActionError('Unable to update GroupAtom due to LOSE_RADICAL action: Invalid radical electron set "{0}".'.format(self.radicalElectrons))    
            radicalElectrons.append(electron)
            
        # Set the new radical electron counts
        self.radicalElectrons = radicalElectrons
        
    def __gainPair(self, pair):
        """
        Update the atom group as a result of applying a GAIN_PAIR action,
        where `pair` specifies the number of lone electron pairs to add.
        """
        lonePairs = []
        if any([len(atomType.incrementLonePair) == 0 for atomType in self.atomType]):
            raise ActionError('Unable to update GroupAtom due to GAIN_PAIR action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        for lonePairs in zip(self.lonePairs):
            lonePairs.append(lonePairs + pair)
        # Set the new lone electron pair count
        self.lonePairs = lonePairs
        
    def __losePair(self, pair):
        """
        Update the atom group as a result of applying a LOSE_PAIR action,
        where `pair` specifies the number of lone electron pairs to remove.
        """
        lonePairs = []
        if any([len(atomType.decrementLonePair) == 0 for atomType in self.atomType]):
            raise ActionError('Unable to update GroupAtom due to LOSE_PAIR action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        for lonePairs in zip(self.lonePairs):
            if lonePairs - pair < 0:
                raise ActionError('Unable to update GroupAtom due to LOSE_PAIR action: Invalid lone electron pairs set "{0}".'.format(self.lonePairs))
            lonePairs.append(lonePairs - pair)
        # Set the new lone electron pair count
        self.lonePairs = lonePairs

    def applyAction(self, action):
        """
        Update the atom group as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        act = action[0].upper()
        if act == 'CHANGE_BOND':
            self.__changeBond(action[2])
        elif act == 'FORM_BOND':
            self.__formBond(action[2])
        elif act == 'BREAK_BOND':
            self.__breakBond(action[2])
        elif act == 'GAIN_RADICAL':
            self.__gainRadical(action[2])
        elif act == 'LOSE_RADICAL':
            self.__loseRadical(action[2])
        elif action[0].upper() == 'GAIN_PAIR':
            self.__gainPair(action[2])
        elif action[0].upper() == 'LOSE_PAIR':
            self.__losePair(action[2])
        else:
            raise ActionError('Unable to update GroupAtom: Invalid action {0}".'.format(action))

    def equivalent(self, other):
        """
        Returns ``True`` if `other` is equivalent to `self` or ``False`` if not,
        where `other` can be either an :class:`Atom` or an :class:`GroupAtom`
        object. When comparing two :class:`GroupAtom` objects, this function
        respects wildcards, e.g. ``R!H`` is equivalent to ``C``.
        
        """
        cython.declare(group=GroupAtom)
        if not isinstance(other, GroupAtom):
            # Let the equivalent method of other handle it
            # We expect self to be an Atom object, but can't test for it here
            # because that would create an import cycle
            return other.equivalent(self)
        group=other
        
        cython.declare(atomType1=AtomType, atomtype2=AtomType, radical1=cython.short, radical2=cython.short,
                       lp1=cython.short, lp2=cython.short, charge1=cython.short, charge2=cython.short)
        # Compare two atom groups for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for atomType1 in self.atomType:
            for atomType2 in group.atomType:
                if atomType1.equivalent(atomType2): break
            else:
                return False
        for atomType1 in group.atomType:
            for atomType2 in self.atomType:
                if atomType1.equivalent(atomType2): break
            else:
                return False
        # Each free radical electron state in self must have an equivalent in other (and vice versa)
        for radical1 in self.radicalElectrons:
            if group.radicalElectrons:  # Only check if the list is non-empty.  An empty list indicates a wildcard.
                for radical2  in group.radicalElectrons:
                    if radical1 == radical2: break
                else:
                    return False
        for radical1 in group.radicalElectrons:
            if self.radicalElectrons:
                for radical2 in self.radicalElectrons:
                    if radical1 == radical2: break
                else:
                    return False
        for lp1 in self.lonePairs:
            if group.lonePairs:
                for lp2 in group.lonePairs:
                    if lp1 == lp2: break
                else:
                    return False
        #Each charge in self must have an equivalent in other (and vice versa)
        for charge1 in self.charge:
            if group.charge:
                for charge2 in group.charge:
                    if charge1 == charge2: break
                else:
                    return False
        for charge1 in group.charge:
            if self.charge:
                for charge2 in self.charge:
                    if charge1 == charge2: break
                else:
                    return False
        # Otherwise the two atom groups are equivalent
        return True

    def isSpecificCaseOf(self, other):
        """
        Returns ``True`` if `other` is the same as `self` or is a more
        specific case of `self`. Returns ``False`` if some of `self` is not
        included in `other` or they are mutually exclusive. 
        """
        cython.declare(group=GroupAtom)
        if not isinstance(other, GroupAtom):
            # Let the isSpecificCaseOf method of other handle it
            # We expect self to be an Atom object, but can't test for it here
            # because that would create an import cycle
            return other.isSpecificCaseOf(self)
        group=other
        
        cython.declare(atomType1=AtomType, atomtype2=AtomType, radical1=cython.short, radical2=cython.short, 
                       lp1=cython.short, lp2=cython.short, charge1=cython.short, charge2=cython.short)
        # Compare two atom groups for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for atomType1 in self.atomType: # all these must match
            for atomType2 in group.atomType: # can match any of these
                if atomType1.isSpecificCaseOf(atomType2): break
            else:
                return False
        # Each free radical electron state in self must have an equivalent in other (and vice versa)
        if self.radicalElectrons:
            for radical1 in self.radicalElectrons:
                if group.radicalElectrons:
                    for radical2 in group.radicalElectrons:
                        if radical1 == radical2: break
                    else:
                        return False
        else:
            if group.radicalElectrons: return False
        if self.lonePairs:
            for lp1 in self.lonePairs:
                if group.lonePairs:
                    for lp2 in group.lonePairs:
                        if lp1 == lp2: break
                    else:
                        return False
        else:
            if group.lonePairs: return False
        #Each charge in self must have an equivalent in other
        if self.charge:
            for charge1 in self.charge:
                if group.charge:
                    for charge2 in group.charge:
                        if charge1 == charge2: break
                    else:
                        return False
        else:
            if group.charge: return False
        # Otherwise self is in fact a specific case of other
        return True

    def isOxygen(self):
        """
        Return ``True`` if the atom represents an oxygen atom or ``False`` if
        not.
        """
        allOxygens = [atomTypes['O']] + atomTypes['O'].specific
        checkList=[x in allOxygens for x in self.atomType]

        return all(checkList)

    def isSulfur(self):
        """
        Return ``True`` if the atom represents an sulfur atom or ``False`` if
        not.
        """
        allSulfur = [atomTypes['S']] + atomTypes['S'].specific
        checkList=[x in allSulfur for x in self.atomType]

        return all(checkList)
################################################################################

class GroupBond(Edge):
    """
    A bond group. This class is based on the :class:`Bond` class, except that
    all attributes are lists rather than individual values. The allowed bond
    types are given :ref:`here <bond-types>`. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `order`             ``list``            The allowed bond orders (as character strings)
    =================== =================== ====================================

    Each list represents a logical OR construct, i.e. a bond will match the
    group if it matches *any* item in the list.
    """

    def __init__(self, atom1, atom2, order=None):
        Edge.__init__(self, atom1, atom2)
        self.order = order or []

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return str(self.order)

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "<GroupBond {0!r}>".format(self.order)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (GroupBond, (self.vertex1, self.vertex2, self.order))

    def copy(self):
        """
        Return a deep copy of the :class:`GroupBond` object. Modifying the
        attributes of the copy will not affect the original.
        """
        return GroupBond(self.vertex1, self.vertex2, self.order[:])

    def isSingle(self):
        """
        Return ``True`` if the bond represents a single bond or ``False`` if
        not. Bonds with any wildcards will return  ``False``.
        """
        return self.order[0] == 'S' and len(self.order) == 1

    def isDouble(self):
        """
        Return ``True`` if the bond represents a double bond or ``False`` if
        not. Bonds with any wildcards will return  ``False``.
        """
        return self.order[0] == 'D' and len(self.order) == 1

    def isTriple(self):
        """
        Return ``True`` if the bond represents a triple bond or ``False`` if
        not. Bonds with any wildcards will return  ``False``.
        """
        return self.order[0] == 'T' and len(self.order) == 1

    def isBenzene(self):
        """
        Return ``True`` if the bond represents a benzene bond or ``False`` if
        not. Bonds with any wildcards will return  ``False``.
        """
        return self.order[0] == 'B' and len(self.order) == 1

    def __changeBond(self, order):
        """
        Update the bond group as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order, and should be 1 or -1.
        """
        newOrder = []
        for bond in self.order:
            if order == 1:
                if bond == 'S':         newOrder.append('D')
                elif bond == 'D':       newOrder.append('T')
                else:
                    raise ActionError('Unable to update GroupBond due to CHANGE_BOND action: Invalid bond order "{0}" in set {1}".'.format(bond, self.order))
            elif order == -1:
                if bond == 'D':         newOrder.append('S')
                elif bond == 'T':       newOrder.append('D')
                else:
                    raise ActionError('Unable to update GroupBond due to CHANGE_BOND action: Invalid bond order "{0}" in set {1}".'.format(bond, self.order))
            else:
                raise ActionError('Unable to update GroupBond due to CHANGE_BOND action: Invalid order "{0}".'.format(order))
        # Set the new bond orders, removing any duplicates
        self.order = list(set(newOrder))

    def applyAction(self, action):
        """
        Update the bond group as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        if action[0].upper() == 'CHANGE_BOND':
            self.__changeBond(action[2])
        else:
            raise ActionError('Unable to update GroupBond: Invalid action {0}".'.format(action))

    def equivalent(self, other):
        """
        Returns ``True`` if `other` is equivalent to `self` or ``False`` if not,
        where `other` can be either an :class:`Bond` or an :class:`GroupBond`
        object.
        """
        cython.declare(gb=GroupBond)
        if not isinstance(other, GroupBond):
            # Let the equivalent method of other handle it
            # We expect self to be a Bond object, but can't test for it here
            # because that would create an import cycle
            return other.equivalent(self)
        gb = other
        
        cython.declare(order1=str, order2=str)
        # Compare two bond groups for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for order1 in self.order:
            for order2 in gb.order:
                if order1 == order2: break
            else:
                return False
        for order1 in gb.order:
            for order2 in self.order:
                if order1 == order2: break
            else:
                return False
        # Otherwise the two bond groups are equivalent
        return True

    def isSpecificCaseOf(self, other):
        """
        Returns ``True`` if `other` is the same as `self` or is a more
        specific case of `self`. Returns ``False`` if some of `self` is not
        included in `other` or they are mutually exclusive.
        """
        cython.declare(gb=GroupBond)
        if not isinstance(other, GroupBond):
            # Let the isSpecificCaseOf method of other handle it
            # We expect self to be a Bond object, but can't test for it here
            # because that would create an import cycle
            return other.isSpecificCaseOf(self)
        gb = other
        
        cython.declare(order1=str, order2=str)
        # Compare two bond groups for equivalence
        # Each atom type in self must have an equivalent in other
        for order1 in self.order: # all these must match
            for order2 in gb.order: # can match any of these
                if order1 == order2: break
            else:
                return False
        # Otherwise self is in fact a specific case of other
        return True

################################################################################

class Group(Graph):
    """
    A representation of a molecular substructure group using a graph data
    type, extending the :class:`Graph` class. The `atoms` and `bonds` attributes
    are aliases for the `vertices` and `edges` attributes, and store 
    :class:`GroupAtom` and :class:`GroupBond` objects, respectively.
    Corresponding alias methods have also been provided.
    """

    def __init__(self, atoms=None, multiplicity=None):
        Graph.__init__(self, atoms)
        self.multiplicity = multiplicity if multiplicity else []
        self.update()

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Group, (self.vertices,))

    def __getAtoms(self): return self.vertices
    def __setAtoms(self, atoms): self.vertices = atoms
    atoms = property(__getAtoms, __setAtoms)

    def addAtom(self, atom):
        """
        Add an `atom` to the graph. The atom is initialized with no bonds.
        """
        return self.addVertex(atom)

    def addBond(self, bond):
        """
        Add a `bond` to the graph as an edge connecting the two atoms `atom1`
        and `atom2`.
        """
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
        return self.removeVertex(atom)

    def removeBond(self, bond):
        """
        Remove the bond between atoms `atom1` and `atom2` from the graph.
        Does not remove atoms that no longer have any bonds as a result of
        this removal.
        """
        return self.removeEdge(bond)

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
        other = cython.declare(Group)
        g = Graph.copy(self, deep)
        other = Group(g.vertices)
        return other

    def update(self):

        self.updateConnectivityValues()
        self.updateFingerprint()


    def merge(self, other):
        """
        Merge two groups so as to store them in a single
        :class:`Group` object. The merged :class:`Group`
        object is returned.
        """
        g = Graph.merge(self, other)
        molecule = Group(atoms=g.vertices)
        return molecule

    def split(self):
        """
        Convert a single :class:`Group` object containing two or more
        unconnected groups into separate class:`Group` objects.
        """
        graphs = Graph.split(self)
        molecules = []
        for g in graphs:
            molecule = Group(atoms=g.vertices)
            molecules.append(molecule)
        return molecules

    def clearLabeledAtoms(self):
        """
        Remove the labels from all atoms in the molecular group.
        """
        for atom in self.vertices:
            atom.label = ''

    def containsLabeledAtom(self, label):
        """
        Return ``True`` if the group contains an atom with the label
        `label` and ``False`` otherwise.
        """
        for atom in self.vertices:
            if atom.label == label: return True
        return False

    def getLabeledAtom(self, label):
        """
        Return the atom in the group that is labeled with the given `label`.
        Raises :class:`ValueError` if no atom in the group has that label.
        """
        for atom in self.vertices:
            if atom.label == label: return atom
        raise ValueError('No atom in the functional group has the label "{0}".'.format(label))

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

    def fromAdjacencyList(self, adjlist):
        """
        Convert a string adjacency list `adjlist` to a molecular structure.
        Skips the first line (assuming it's a label) unless `withLabel` is
        ``False``.
        """
        from .adjlist import fromAdjacencyList
        self.vertices, multiplicity = fromAdjacencyList(adjlist, group=True)
        if multiplicity is not None:
            self.multiplicity = multiplicity
        self.update()
        return self

    def toAdjacencyList(self, label=''):
        """
        Convert the molecular structure to a string adjacency list.
        """
        from .adjlist import toAdjacencyList
        return toAdjacencyList(self.vertices, multiplicity=self.multiplicity, label='', group=True)


    def updateFingerprint(self):
        """
        Update the molecular fingerprint used to accelerate the subgraph
        isomorphism checks.
        """
        cython.declare(atom=GroupAtom, atomType=AtomType)
        cython.declare(carbon=AtomType, nitrogen=AtomType, oxygen=AtomType, sulfur=AtomType)
        cython.declare(isCarbon=cython.bint, isNitrogen=cython.bint, isOxygen=cython.bint, isSulfur=cython.bint, radical=cython.int)
        
        carbon   = atomTypes['C']
        nitrogen = atomTypes['N']
        oxygen   = atomTypes['O']
        sulfur   = atomTypes['S']
        
        self.carbonCount   = 0
        self.nitrogenCount = 0
        self.oxygenCount   = 0
        self.sulfurCount   = 0
        self.radicalCount  = 0
        for atom in self.vertices:
            if len(atom.atomType) == 1:
                atomType   = atom.atomType[0]
                isCarbon   = atomType.equivalent(carbon)
                isNitrogen = atomType.equivalent(nitrogen)
                isOxygen   = atomType.equivalent(oxygen)
                isSulfur   = atomType.equivalent(sulfur)
                if isCarbon and not isNitrogen and not isOxygen and not isSulfur:
                    self.carbonCount += 1
                elif isNitrogen and not isCarbon and not isOxygen and not isSulfur:
                    self.nitrogenCount += 1
                elif isOxygen and not isCarbon and not isNitrogen and not isSulfur:
                    self.oxygenCount += 1
                elif isSulfur and not isCarbon and not isNitrogen and not isOxygen:
                    self.sulfurCount += 1
            if len(atom.radicalElectrons) == 1:
                radical = atom.radicalElectrons[0]
                self.radicalCount += radical

    def isIsomorphic(self, other, initialMap=None):
        """
        Returns ``True`` if two graphs are isomorphic and ``False``
        otherwise. The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """
        # It only makes sense to compare a Group to a Group for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        # Do the isomorphism comparison
        return Graph.isIsomorphic(self, other, initialMap)

    def findIsomorphism(self, other, initialMap=None):
        """
        Returns ``True`` if `other` is isomorphic and ``False``
        otherwise, and the matching mapping. The `initialMap` attribute can be
        used to specify a required mapping from `self` to `other` (i.e. the
        atoms of `self` are the keys, while the atoms of `other` are the
        values). The returned mapping also uses the atoms of `self` for the keys
        and the atoms of `other` for the values. The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """
        # It only makes sense to compare a Group to a Group for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        # Do the isomorphism comparison
        return Graph.findIsomorphism(self, other, initialMap)

    def isSubgraphIsomorphic(self, other, initialMap=None):
        """
        Returns ``True`` if `other` is subgraph isomorphic and ``False``
        otherwise. In other words, return ``True`` if self is more specific than other.
        The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """        
        cython.declare(group=Group)
        cython.declare(mult1=cython.short, mult2=cython.short)
        # It only makes sense to compare a Group to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        group = other
        
        if self.multiplicity:
            for mult1 in self.multiplicity:
                if group.multiplicity:
                    for mult2 in group.multiplicity:
                        if mult1 == mult2: break
                    else:
                        return False
        else:
            if group.multiplicity: return False
        # Do the isomorphism comparison
        return Graph.isSubgraphIsomorphic(self, other, initialMap)

    def findSubgraphIsomorphisms(self, other, initialMap=None):
        """
        Returns ``True`` if `other` is subgraph isomorphic and ``False``
        otherwise. In other words, return ``True`` is self is more specific than other.
        Also returns the lists all of valid mappings. The
        `initialMap` attribute can be used to specify a required mapping from
        `self` to `other` (i.e. the atoms of `self` are the keys, while the
        atoms of `other` are the values). The returned mappings also use the
        atoms of `self` for the keys and the atoms of `other` for the values.
        The `other` parameter must be a :class:`Group` object, or a
        :class:`TypeError` is raised.
        """
        cython.declare(group=Group)
        cython.declare(mult1=cython.short, mult2=cython.short)

        # It only makes sense to compare a Group to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        group = other
        
        if self.multiplicity:
            for mult1 in self.multiplicity:
                if group.multiplicity:
                    for mult2 in group.multiplicity:
                        if mult1 == mult2: break
                    else:
                        return []
        else:
            if group.multiplicity: return []
                
        # Do the isomorphism comparison
        return Graph.findSubgraphIsomorphisms(self, other, initialMap)
    
    def isIdentical(self, other):
        """
        Returns ``True`` if `other` is identical and ``False`` otherwise.
        The function `isIsomorphic` respects wildcards, while this function
        does not, make it more useful for checking groups to groups (as
        opposed to molecules to groups)
        """
        # It only makes sense to compare a Group to a Group for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        # An identical group is always a child of itself and 
        # is the only case where that is true. Therefore
        # if we do both directions of isSubgraphIsmorphic, we need
        # to get True twice for it to be identical
        if not self.isSubgraphIsomorphic(other):
            return False
        elif not other.isSubgraphIsomorphic(self):
            return False
        else:
            return True

    def standardizeAtomType(self):
        """
        This function changes the atomTypes in a group if the atom must
        be a specific atomType based on its bonds and valency.

        Currently only standardizes oxygen, carbon and sulfur atomTypes

        We also only check when there is exactly one atomType,
        one bondType, one radical setting.
        For any group where there are wildcards or multiple attributes,
        we cannot apply this check.

        In the case where the atomType is ambigious based on bonds
        and valency, this function will not change the type.

        Returns a 'True' if the group was modified otherwise returns 'False'
        """

        modified = False

        #dictionary of element to expected valency
        valency = {atomTypes['C'] : 4,
                   atomTypes['O'] : 2,
                   atomTypes['S']: 2,
                   atomTypes['Si']:4
                   }

        #list of :class:AtomType which are elements with more sub-divided atomtypes beneath them
        specifics= [elementLabel for elementLabel in allElements if elementLabel not in nonSpecifics]
        for index, atom in enumerate(self.atoms):
            claimedAtomType = atom.atomType[0]
            newAtomType = None
            element = None
            #Ignore elements that do not have more than one atomtype
            if claimedAtomType.label in nonSpecifics: continue
            for elementLabel in specifics:
                if claimedAtomType.label == elementLabel or atomTypes[claimedAtomType.label] in atomTypes[elementLabel].specific:
                    element = atomTypes[elementLabel]
                    break

            #claimedAtomType is not in one of the specified elements
            if not element: continue
            #Don't standardize atomtypes for nitrogen for now. My feeling is that
            # the work on the nitrogen atomtypes is still incomplete
            elif element is atomTypes['N']: continue

            groupFeatures = getFeatures(atom, atom.bonds)

            single = groupFeatures[0]
            allDouble = groupFeatures[1]
            triple = groupFeatures[5]
            benzene = groupFeatures[6]
            if benzene == 3:
                bondValency = single + 2 * allDouble + 3 * triple + 4.0/3.0 * benzene
            else:
                bondValency =  single + 2 * allDouble + 3 * triple + 3.0/2.0 * benzene
            filledValency =  atom.radicalElectrons[0] + bondValency

            #For an atomtype to be known for certain, the valency must be filled
            #within 1 of the total valency available
            if filledValency >= valency[element] - 1:
                for specificAtomType in element.specific:
                    atomtypeFeatureList = specificAtomType.getFeatures()
                    for molFeature, atomtypeFeature in zip(groupFeatures, atomtypeFeatureList):
                        if atomtypeFeature == []:
                            continue
                        elif molFeature not in atomtypeFeature:
                            break
                    else:
                        if specificAtomType is atomTypes['Oa'] or specificAtomType is atomTypes['Sa']:
                            if atom.lonePairs == 3 or atom.radicalElectrons == 2:
                                newAtomType = specificAtomType
                                break
                        else:
                            newAtomType = specificAtomType
                            break

            #set the new atom type if the algorithm found one
            if newAtomType and newAtomType is not claimedAtomType:
                atom.atomType[0] = newAtomType
                modified = True

        return modified

    def addExplicitLigands(self):
        """
        This function Od/Sd ligand to CO or CS atomtypes if they are not already there.

        Returns a 'True' if the group was modified otherwise returns 'False'
        """

        modified = False

        atomsToAddTo=[]

        for index, atom in enumerate(self.atoms):
            claimedAtomType = atom.atomType[0]
            if claimedAtomType is atomTypes['CO'] or claimedAtomType is atomTypes['CS']:
                for atom2, bond12 in atom.bonds.iteritems():
                    if bond12.isDouble():
                        break
                else: atomsToAddTo.append(index)

        for atomIndex in atomsToAddTo:
            modified = True
            if self.atoms[atomIndex].atomType[0] is atomTypes['CO']:
                newAtom = GroupAtom(atomType=[atomTypes['Od']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
            elif self.atoms[atomIndex].atomType[0] is atomTypes['CS']:
                newAtom = GroupAtom(atomType=[atomTypes['Sd']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
            self.addAtom(newAtom)
            newBond = GroupBond(self.atoms[atomIndex], newAtom, order=['D'])
            self.addBond(newBond)

        return modified

    def standardizeGroup(self):
        """
        This function modifies groups to make them have a standard AdjList form.

        Currently it makes atomtypes as specific as possible and makes CO/CS atomtypes
        have explicit Od/Sd ligands. Other functions can be added as necessary

        We also only check when there is exactly one atomType, one bondType, one
        radical setting. For any group where there are wildcards or multiple
        attributes, we do not apply this check.

        Returns a 'True' if the group was modified otherwise returns 'False'
        """

        modified = False

        #see if this is a group we can check, must not have any OR groups in
        #its bonds or atomtypes
        viableToCheck = True
        for atom in self.atoms:
            if len(atom.atomType) > 1:
                viableToCheck = False
                break
            elif len(atom.radicalElectrons) > 1:
                viableToCheck = False
                break
            elif len(atom.lonePairs) > 1:
                viableToCheck = False
                break
            for bond in atom.bonds.values():
                if len(bond.order) > 1:
                    viableToCheck = False
                    break
        if not viableToCheck: return modified

        #If viable then we apply current conventions:
        checkList=[]
        checkList.append(self.standardizeAtomType())
        checkList.append(self.addExplicitLigands())

        return any(checkList)
