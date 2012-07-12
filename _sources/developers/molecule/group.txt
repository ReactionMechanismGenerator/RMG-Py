*************************************************
:mod:`rmgpy.molecule.group` --- Functional Groups
*************************************************

.. automodule:: rmgpy.molecule.group

Bond Types
==========

.. _bond-types:

The bond type simply indicates the order of a chemical bond. We define the 
following bond types:

=============== ============================================================
Bond type       Description
=============== ============================================================
``S``           a single bond
``D``           a double bond
``T``           a triple bond
``B``           a benzene bond
=============== ============================================================

Reaction Recipe Actions
=======================

.. _reaction-recipe-actions:

A reaction recipe is a procedure for applying a reaction to a set of chemical
species. Each reaction recipe is made up of a set of actions that, when applied 
sequentially, a set of chemical reactants to chemical products via that
reaction's characteristic chemical process. Each action requires a small set of
parameters in order to be fully defined.

We define the following reaction recipe actions:

============= ============================= ================================
Action name   Arguments                     Action
============= ============================= ================================
CHANGE_BOND   `center1`, `order`, `center2` change the bond order of the bond between `center1` and `center2` by `order`; do not break or form bonds
FORM_BOND     `center1`, `order`, `center2` form a new bond between `center1` and `center2` of type `order`
BREAK_BOND    `center1`, `order`, `center2` break the bond between `center1` and `center2`, which should be of type `order`
GAIN_RADICAL  `center`, `radical`           increase the number of free electrons on `center` by `radical`
LOSE_RADICAL  `center`, `radical`           decrease the number of free electrons on `center` by `radical`
============= ============================= ================================

.. autoclass:: rmgpy.molecule.group.ActionError
    :members:

Adjacency Lists
===============

An adjacency list is the most general way of specifying a chemical molecule or
molecular pattern in RMG. It is based on the adjacency list representation of
the graph data type --  the underlying data type for molecules and patterns in 
RMG -- but extended to allow for specification of extra semantic information.

The first line of most adjacency lists is a unique identifier for the molecule
or pattern the adjacency list represents. This is not strictly required, but
is recommended in most cases. Generally the identifier should only use
alphanumeric characters and the underscore, as if an identifer in many popular
programming languages. However, strictly speaking any non-space ASCII character
is allowed.

After the identifier line, each subsequent line describes a single atom and its
local bond structure. The format of these lines is a whitespace-delimited list
with tokens ::

    <number> [<label>] <element> <radicals> <bondlist>

The first item is the number used to identify that atom. Any number may be used,
though it is recommended to number the atoms sequentially starting from one.
Next is an optional label used to tag that atom; this should be an
asterisk followed by a unique number for the label, e.g. ``*1``. After that is
the atom's element, indicated by its atomic symbol, followed by the number of
radical electrons on the atom. The last set of tokens is the list of bonds.
To indicate a bond, place the number of the atom at the other end of the bond
and the bond type within curly braces and separated by a comma, e.g. ``{2,S}``.
Multiple bonds to the same atom should be separated by whitespace.

.. note::
    You must take care to make sure each bond is listed on the lines of *both*
    atoms in the bond, and that these entries have the same bond type. RMG will
    raise an exception if it encounters such an invalid adjacency list.

When writing a molecular substructure pattern, you may specify multiple 
elements, radical counts, and bond types as a comma-separated list inside curly
braces. For example, to specify any carbon or oxygen atom, use the syntax 
``{C,O}``. Atom types may also be used as a shorthand. (Atom types can also be
used in full molecules, but this use is discouraged, as RMG can compute them
automatically for full molecules.)

Below is an example adjacency list, for 1,3-hexadiene, with the weakest bond in
the molecule labeled with ``*1`` and ``*2``. Note that hydrogen atoms
can be omitted if desired, as their presence is inferred::

    HXD13
    1    C 0       {2,D}
    2    C 0 {1,D} {3,S}
    3    C 0 {2,S} {4,D}
    4    C 0 {3,D} {5,S}
    5 *1 C 0 {4,S} {6,S}
    6 *2 C 0 {5,S}

.. autofunction:: rmgpy.molecule.adjlist.fromAdjacencyList

.. autofunction:: rmgpy.molecule.adjlist.toAdjacencyList


GroupAtom Objects
=================

.. autoclass:: rmgpy.molecule.group.GroupAtom
    :members:

GroupBond Objects
=================

.. autoclass:: rmgpy.molecule.group.GroupBond
    :members:

Group Objects
=============

.. autoclass:: rmgpy.molecule.group.Group
    :members:

