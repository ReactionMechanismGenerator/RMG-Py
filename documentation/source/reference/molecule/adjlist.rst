.. _rmgpy.molecule.adjlist:


***************
Adjacency Lists
***************

.. module:: rmgpy.molecule.adjlist

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
 +----------------------+----------+---------------------+
 | Nonreactive Elements | N        | Nitrogen atom       |
 |                      +----------+---------------------+
 |                      | Si       | Silicon atom        |
 |                      +----------+---------------------+
 |                      | Cl       | Chlorine atom       |
 |                      +----------+---------------------+
 |                      | He       | Helium atom         |
 |                      +----------+---------------------+
 |                      | Ar       | Argon atom          |
 |                      +----------+---------------------+
 +----------------------+----------+---------------------+
 | Free Electrons       | 0        | Non-radical         |
 |                      +----------+---------------------+
 |                      | 1        | Mono-radical        |
 |                      +----------+---------------------+
 |                      | 2        | Bi-radical          |
 |                      +----------+---------------------+
 |                      | 2T       | Triplet             |
 |                      +----------+---------------------+
 |                      | 2S       | Singlet             |
 |                      +----------+---------------------+
 |                      | 3        | Tri-radical         |
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
