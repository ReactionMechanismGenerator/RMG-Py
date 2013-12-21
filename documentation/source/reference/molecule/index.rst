*************************************************
Molecular representations (:mod:`rmgpy.molecule`)
*************************************************

.. module:: rmgpy.molecule

The :mod:`rmgpy.molecule` subpackage contains classes and functions for working
with molecular representations, particularly using chemical graph theory.



Graphs
======

.. currentmodule:: rmgpy.molecule.graph

======================= ========================================================
Class                   Description
======================= ========================================================
:class:`Vertex`         A generic vertex (node) in a graph
:class:`Edge`           A generic edge (arc) in a graph
:class:`Graph`          A generic graph data type
======================= ========================================================



Graph isomorphism
=================

.. currentmodule:: rmgpy.molecule.vf2

======================= ========================================================
Class                   Description
======================= ========================================================
:class:`VF2`            Graph isomorphism using the VF2 algorithm
======================= ========================================================



Elements and atom types
=======================

.. currentmodule:: rmgpy.molecule

======================= ========================================================
Class/Function          Description
======================= ========================================================
:class:`Element`        A model of a chemical element
:func:`getElement`      Return the :class:`Element` object for a given atomic number or symbol
:class:`AtomType`       A model of an atom type: an element and local bond structure
:func:`getAtomType`     Return the :class:`AtomType` object for a given atom in a molecule
======================= ========================================================



Molecules
=========

.. currentmodule:: rmgpy.molecule

======================= ========================================================
Class                   Description
======================= ========================================================
:class:`Atom`           An atom in a molecule
:class:`Bond`           A bond in a molecule
:class:`Molecule`       A molecular structure represented using a chemical graph
======================= ========================================================



Functional groups
=================

.. currentmodule:: rmgpy.molecule

======================= ========================================================
Class                   Description
======================= ========================================================
:class:`GroupAtom`      An atom in a functional group
:class:`GroupBond`      A bond in a functional group
:class:`Group`          A functional group structure represented using a chemical graph
======================= ========================================================


Adjacency lists
===============

.. currentmodule:: rmgpy.molecule.adjlist

=========================== ====================================================
Function                    Description
=========================== ====================================================
:func:`fromAdjacencyList`   Convert an adjacency list to a set of atoms and bonds
:func:`toAdjacencyList`     Convert a set of atoms and bonds to an adjacency list
=========================== ====================================================



Symmetry numbers
================

.. currentmodule:: rmgpy.molecule.symmetry

======================================= ========================================
Class                                   Description
======================================= ========================================
:func:`calculateAtomSymmetryNumber`     Calculate the atom-centered symmetry number for an atom in a molecule
:func:`calculateBondSymmetryNumber`     Calculate the bond-centered symmetry number for a bond in a molecule
:func:`calculateAxisSymmetryNumber`     Calculate the axis-centered symmetry number for a double bond axis in a molecule
:func:`calculateCyclicSymmetryNumber`   Calculate the ring-centered symmetry number for a ring in a molecule
:func:`calculateSymmetryNumber`         Calculate the total internal + external symmetry number for a molecule
======================================= ========================================



Molecule and reaction drawing
=============================

.. currentmodule:: rmgpy.molecule.draw

======================== =======================================================
Class                    Description
======================== =======================================================
:class:`MoleculeDrawer`  Draw the skeletal formula of a molecule
:class:`ReactionDrawer`  Draw a chemical reaction
======================== =======================================================


Exceptions
==========

.. currentmodule:: rmgpy.molecule

=================================== ============================================
Exception                           Description
=================================== ============================================
:exc:`ElementError`                 Raised when an error occurs while working with chemical elements
:exc:`AtomTypeError`                Raised when an error occurs while working with atom types
:exc:`InvalidAdjacencyListError`    Raised when an invalid adjacency list is encountered
:exc:`ActionError`                  Raised when an error occurs while working with a reaction recipe action
=================================== ============================================


.. toctree::
    :hidden:
    
    vertex
    edge
    graph
    vf2
    element
    atomtype
    recipe
    atom
    bond
    molecule
    groupatom
    groupbond
    group
    adjlist
    symmetry
    moleculedrawer
    reactiondrawer
