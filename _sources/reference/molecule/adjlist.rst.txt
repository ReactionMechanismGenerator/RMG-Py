.. _rmgpy.molecule.adjlist:


***************
Adjacency Lists
***************

.. module:: rmgpy.molecule.adjlist


.. note::
    The adjacency list syntax changed in July 2014.
    The minimal requirement for most translations is to prefix the number
    of unpaired electrons with the letter `u`.
    The new syntax, however, allows much
    greater flexibility, including definition of lone pairs, partial charges, 
    wildcards, and molecule multiplicities.

.. note::
    To quickly visualize any adjacency list, or to generate an adjacency list from
    other types of molecular representations such as SMILES, InChI, or even common
    species names, use the Molecule Search tool found here: http://rmg.mit.edu/molecule_search


An adjacency list is the most general way of specifying a chemical molecule or
molecular pattern in RMG. It is based on the adjacency list representation of
the graph data type --  the underlying data type for molecules and patterns in 
RMG -- but extended to allow for specification of extra semantic information.

The first line of most adjacency lists is a unique identifier for the molecule
or pattern the adjacency list represents. This is not strictly required, but
is recommended in most cases. Generally the identifier should only use
alphanumeric characters and the underscore, as if an identifier in many popular
programming languages. However, strictly speaking any non-space ASCII character
is allowed.

The subsequent lines may contain keyword-value pairs. Currently there is only
one keyword, ``multiplicity``.

For species or molecule declarations, the value after ``multiplicity`` defines 
the spin multiplicity of the molecule. E.g. ``multiplicity 1`` for most ground state 
closed shell species, ``multiplicity 2`` for most radical species, 
and ``multiplicity 3`` for a triplet biradical.
If the ``multiplicity`` line is not present then a value of 
(1 + number of unpaired electrons) is assumed. 
Thus, it can usually be omitted, but if present can be used to distinguish,
for example, singlet CH2 from triplet CH2.

If defining a Functional :class:`~rmgpy.molecule.Group`, then the value must be a list,
which defines the multiplicities that will be matched by the group, eg. 
``multiplicity [1,2,3]`` or, for a single value, ``multiplicity [1]``. 
If a wildcard is desired, the line ``'multiplicity x`` can be used instead to accept
all multiplicities.  If the multiplicity line is omitted altogether, then a wildcard 
is assumed.

e.g. the following two group adjlists represent identical groups. ::

    group1
    multiplicity x
    1    R!H u0

::

    group2
    1    R!H u0

After the identifier line and keyword-value lines,
each subsequent line describes a single atom and its
local bond structure. The format of these lines is a whitespace-delimited list
with tokens ::

    <number> [<label>] <element> u<unpaired> [p<pairs>] [c<charge>] <bondlist>

The first item is the number used to identify that atom. Any number may be used,
though it is recommended to number the atoms sequentially starting from one.
Next is an optional label used to tag that atom; this should be an
asterisk followed by a unique number for the label, e.g. ``*1``.
In some cases (e.g. thermodynamics groups) there is only one labeled atom, and the label 
is just an asterisk with no number: ``*``.

After that is
the atom's element or atom type, indicated by its atomic symbol, followed by 
a sequence of tokens describing the electronic state of the atom:

* ``u0`` number of **unpaired** electrons (eg. radicals)
* ``p0`` number of lone **pairs** of electrons, common on oxygen and nitrogen.
* ``c0`` formal **charge** on the atom, e.g. ``c-1`` (negatively charged),
  ``c0``, ``c+1`` (positively charged)

For :class:`~rmgpy.molecule.Molecule` definitions:
The value must be a single integer (and for charge must have a + or - sign if not equal to 0)
The number of unpaired electrons (i.e. radical electrons) is required, even if zero.
The number of lone pairs and the formal charge are assumed to be zero if omitted.

For :class:`~rmgpy.molecule.Group` definitions:
The value can be an integer or a list of integers (with signs, for charges),
eg. ``u[0,1,2]`` or ``c[0,+1,+2,+3,+4]``, or may be a wildcard ``x`` 
which matches any valid value, 
eg. ``px`` is the same as ``p[0,1,2,3,4, ...]`` and ``cx`` is the same as
``c[...,-4,-3,-2,-1,0,+1,+2,+3,+4,...]``. Lists must be enclosed is square brackets,
and separated by commas, without spaces.
If lone pairs or formal charges are omitted from a group definition,
the wildcard is assumed.


The last set of tokens is the list of bonds.
To indicate a bond, place the number of the atom at the other end of the bond
and the bond type within curly braces and separated by a comma, e.g. ``{2,S}``.
Multiple bonds from the same atom should be separated by whitespace.

.. note::
    You must take care to make sure each bond is listed on the lines of *both*
    atoms in the bond, and that these entries have the same bond type. RMG will
    raise an exception if it encounters such an invalid adjacency list.


When writing a molecular substructure pattern, you may specify multiple 
elements, radical counts, and bond types as a comma-separated list inside square
brackets. For example, to specify any carbon or oxygen atom, use the syntax 
``[C,O]``. For a single or double bond to atom 2, write ``{2,[S,D]}``.

Atom types such as ``R!H`` or ``Cdd`` may also be used as a shorthand. (Atom types 
like ``Cdd`` can also be
used in full molecules, but this use is discouraged, as RMG can compute them
automatically for full molecules.)

Below is an example adjacency list, for 1,3-hexadiene, with the weakest bond in
the molecule labeled with ``*1`` and ``*2``. Note that hydrogen atoms
can be omitted if desired, as their presence is inferred, provided that unpaired
electrons, lone pairs, and charges are all correctly defined::

    HXD13
    multiplicity 1
    1    C u0       {2,D}
    2    C u0 {1,D} {3,S}
    3    C u0 {2,S} {4,D}
    4    C u0 {3,D} {5,S}
    5 *1 C u0 {4,S} {6,S}
    6 *2 C u0 {5,S}
    
The allowed element types, radicals, and bonds are listed in the following table:

 +----------------------+----------+---------------------+
 |                      | Notation | Explanation         |
 +======================+==========+=====================+
 | Chemical Element     | C        | Carbon atom         |
 |                      +----------+---------------------+
 |                      | O        | Oxygen atom         | 
 |                      +----------+---------------------+
 |                      | H        | Hydrogen atom       |
 |                      +----------+---------------------+
 |                      | S        | Sulfur atom         |
 |                      +----------+---------------------+
 |                      | N        | Nitrogen atom       |
 +----------------------+----------+---------------------+
 | Nonreactive Elements | Si       | Silicon atom        |
 |                      +----------+---------------------+
 |                      | Cl       | Chlorine atom       |
 |                      +----------+---------------------+
 |                      | He       | Helium atom         |
 |                      +----------+---------------------+
 |                      | Ar       | Argon atom          |
 |                      +----------+---------------------+
 +----------------------+----------+---------------------+
 | Chemical Bond        | S        | Single Bond         |
 |                      +----------+---------------------+
 |                      | D        | Double Bond         |
 |                      +----------+---------------------+
 |                      | T        | Triple bond         |
 |                      +----------+---------------------+
 |                      | B        | Benzene bond        |
 +----------------------+----------+---------------------+


.. autofunction:: rmgpy.molecule.adjlist.fromAdjacencyList

.. autofunction:: rmgpy.molecule.adjlist.toAdjacencyList
