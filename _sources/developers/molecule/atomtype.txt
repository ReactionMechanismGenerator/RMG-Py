*********************************************
:mod:`rmgpy.molecule.atomtype` --- Atom Types
*********************************************

.. automodule:: rmgpy.molecule.atomtype

Atom Types
==========

.. _atom-types:

The atom type of an atom describes the atom itself and (often) something about
the local bond structure around that atom. This is a useful semantic tool for
accelerating graph isomorphism queries, and a useful shorthand when specifying
molecular substructure patterns via an RMG-style adjacency list.

We define the following basic atom types:

=============== ============================================================
Atom type       Description
=============== ============================================================
*General atom types*
----------------------------------------------------------------------------
``R``           any atom with any local bond structure
``R!H``         any non-hydrogen atom with any local bond structure
*Carbon atom types*
----------------------------------------------------------------------------
``C``           carbon atom with any local bond structure
``Cs``          carbon atom with four single bonds
``Cd``          carbon atom with one double bond (to carbon) and two single bonds
``Cdd``         carbon atom with two double bonds
``Ct``          carbon atom with one triple bond and one single bond
``CO``          carbon atom with one double bond (to oxygen) and two single bonds
``Cb``          carbon atom with two benzene bonds and one single bond
``Cbf``         carbon atom with three benzene bonds
*Hydrogen atom types*
----------------------------------------------------------------------------
``H``           hydrogen atom with one single bond
*Oxygen atom types*
----------------------------------------------------------------------------
``O``           oxygen atom with any local bond structure
``Os``          oxygen atom with two single bonds
``Od``          oxygen atom with one double bond
``Oa``          oxygen atom with no bonds
*Silicon atom types*
----------------------------------------------------------------------------
``Si``          silicon atom with any local bond structure
``Sis``         silicon atom with four single bonds
``Sid``         silicon atom with one double bond (to carbon) and two single bonds
``Sidd``        silicon atom with two double bonds
``Sit``         silicon atom with one triple bond and one single bond
``SiO``         silicon atom with one double bond (to oxygen) and two single bonds
``Sib``         silicon atom with two benzene bonds and one single bond
``Sibf``        silicon atom with three benzene bonds
*Sulfur atom types*
----------------------------------------------------------------------------
``S``           sulfur atom with any local bond structure
``Ss``          sulfur atom with two single bonds
``Sd``          sulfur atom with one double bond
``Sa``          sulfur atom with no bonds
=============== ============================================================

.. autoclass:: rmgpy.molecule.atomtype.AtomTypeError
    :members:

.. autoclass:: rmgpy.molecule.atomtype.AtomType
    :members:

.. autofunction:: rmgpy.molecule.atomtype.getAtomType

