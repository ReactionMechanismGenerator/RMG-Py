.. _species:

**********************
Species Representation
**********************

The base class for chemical structures in RMG is ``Graph`` (see :ref:`rmgpy.molecule.Graph`), 
which is a basic implementation of a 2D mathematical graph. 
A graph is comprised of a set of vertices connected by a set of edges. 
In RMG, the Graph class does not store any chemical information on its own, 
but it is the parent class of Group and Molecule.

The ``Group`` class is used to represent a molecular fragment, whereas the ``Molecule`` class 
is used to represent a specific molecular structure. 
A ``Group`` object can have a list of allowed specifications for each atom or bond property, 
while a ``Molecule`` object  can only have one. Additionally, a ``Group`` object does not 
have to have a complete molecule while a Molecule object must be a full chemical structure.
See the :ref:`rmgpy.molecule.Group` and :ref:`rmgpy.molecule.Molecule` class documentation.

Finally, there is the ``Species`` class. Although colloquially we use molecules and species 
interchangeably, these terms have precise meanings in RMG. 
A ``Species`` object contains a list of molecule objects which are different representations of 
the same chemical compound (i.e., resonance structures). 
It also contains a descriptive label for the species, its thermochemical properties, transport data, 
and molecular weight, among other information. This distinction is important because many chemical 
compounds have resonance structures. Therefore, a ``Molecule`` object is a graph that denotes the structure 
of a specific resonance isomer, while a ``Species`` incorporates all resonance structures found into one object. 
This prevents duplicate reactions from being generated for two molecule objects that are really the same chemical compound.
See the :ref:`rmgpy.species.Species` class documentation for more information.

RMG considers each species to be unique, and comprised of a set of molecular structural 
isomers, including resonance isomers.  RMG uses the list of resonance isomers to 
compare whether two species are the same. Each molecular structure is stored in RMG using
graph representations, using vertices and edges, where the vertices are the atoms and the 
edges are the bonds.  This form of representation is known as an adjacency list.  For
more information on adjacency lists, see the :ref:`rmgpy.molecule.adjlist` page.  

Species objects in the input file can also be constructed using other common representations
such as SMILES, SMARTS, and InChIs.  The following can all be used to represent the methane 
species: ::

    species(
        label='CH4',
        reactive=True,
        structure=SMILES("C"),
    )

Replacing the structure with any of the following representations will also produce
the same species: ::

    structure=adjacencyList("
    1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
    2 H u0 p0 c0 {1,S}
    3 H u0 p0 c0 {1,S}
    4 H u0 p0 c0 {1,S}
    5 H u0 p0 c0 {1,S}
    "),
   
    structure=SMARTS("[CH4]"),
   
    structure=SMILES("C"),
   
    structure=InChI("InChI=1S/CH4/h1H4"),

Be careful using the SMILES shorthand with lowercase letters for aromatics,
radicals, or double bonds, because these can be ambiguous and the resulting
molecule may depend on the version of OpenBabel, RDKit, and RMG in use.
To quickly generate any adjacency list, or to generate an adjacency list from
other types of molecular representations such as SMILES, InChI, or even common
species names, use the Molecule Search tool found here: https://rmg.mit.edu/molecule_search

.. _representing_oxygen:

Representing Oxygen
===================

Special care should be taken when constructing a mechanism that involves 
molecular oxygen. The ground electronic state of molecular oxygen,
:math:`^3\Sigma^-_g`, does *not* contain a double bond, but instead a single
bond and two lone electrons. In RMG's adjaceny list notation the ground state
of oxygen is represented as ::

   1 O u1 p2 {2,S}
   2 O u1 p2 {1,S}

You should use the above adjacency list to represent molecular oxygen in
your condition files, seed mechanisms, etc. The triplet form is 22 kcal/mol
more stable than the first singlet excited state, :math:`^1\Delta_g`, which 
does contain a double bond. The adjacency list for singlet oxygen is ::

   1 O u0 p2 {2,D}
   2 O u0 p2 {1,D}

Selecting the correct structure for oxygen is important, as the reactions
generated from a double bond are significantly different than those generated
from a radical or diradical. For example, the reaction

.. math:: \mathrm{CH_4} + \mathrm{O_2} \rightarrow \mathrm{CH_3} + \mathrm{HO_2}

would occur for both triplet and singlet oxygen, but in entirely different
families. For triplet oxygen the above represents a hydrogen abstraction, while
for singlet oxygen it represents the reverse of a disproportionation reaction.

The RMG databases have been modified to make all of the
oxygen-related chemistry that was present in RMG databases consistent with the
single-bonded biradical representation.

Conversion between triplet and singlet forms is possible through the primary
reaction library ``OxygenSingTrip``; the reactions involved are very slow, however,
and are likely to be absent from any mechanisms generated. At this point, no other
reactions of singlet oxygen have been included in RMG.

In order to allow the singlet form of O2 to be used in RMG, please allow it explicitly by
setting ``allowSingletO2`` to ``True`` in the ``generateSpeciesConstraints`` section of the
RMG input file. ::

    generatedSpeciesConstraints(
        allowSingletO2 = True,
    )
