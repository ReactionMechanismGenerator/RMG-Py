.. _species:

**********************
Species Representation
**********************

Species objects in RMG contain a variety of attributes, including user given names, 
thermochemistry, as well as structural isomers.  See the :ref:`rmgpy.species.Species` class
documentation for more information.

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

    structure=adjacencyList("1 C 0"),
   
    structure=SMARTS("[CH4]"),
   
    structure=SMILES("C"),
   
    structure=InChI("InChI=1S/CH4/h1H4"),


To quickly generate any adjacency list, or to generate an adjacency list from
other types of molecular representations such as SMILES, InChI, or even common
species names, use the Molecule Search tool found here: http://rmg.mit.edu/molecule_search


Representing Oxygen
===================

Special care should be taken when constructing a mechanism that involves 
molecular oxygen. The ground electronic state of molecular oxygen,
:math:`^3\Sigma^-_g`, does *not* contain a double bond, but instead a single
bond and two lone electrons. In RMG's adjaceny list notation the ground state
of oxygen is represented as ::

   1 O 1 {2,S}
   2 O 1 {1,S}

You should use the above adjacency list to represent molecular oxygen in
your condition files, seed mechanisms, etc. The triplet form is 22 kcal/mol
more stable than the first singlet excited state, :math:`^1\Delta_g`, which 
does contain a double bond. The adjacency list for singlet oxygen is ::

   1 O 0 {2,D}
   2 O 0 {1,D}

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