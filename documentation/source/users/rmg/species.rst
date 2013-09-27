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
