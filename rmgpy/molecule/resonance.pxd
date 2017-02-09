from .graph cimport Vertex, Edge, Graph
from .molecule cimport Atom, Bond, Molecule

cpdef tuple populate_resonance_generation_algorithm()

cpdef list populateResonanceAlgorithms(dict features=?)

cpdef dict analyzeMolecule(Molecule mol)

cpdef list generateResonanceIsomers(Molecule mol)

cpdef list generateAdjacentResonanceIsomers(Molecule mol)

cpdef list generateLonePairRadicalResonanceIsomers(Molecule mol)

cpdef list generateN5dd_N5tsResonanceIsomers(Molecule mol)

cpdef list generate_isomorphic_isomers(Molecule mol)

cpdef list generateKekulizedResonanceIsomers(Molecule mol)

cpdef list generateOppositeKekuleStructure(Molecule mol)

cpdef list generateClarStructures(Molecule mol)

cpdef list clarOptimization(Molecule mol, list constraints=?, maxNum=?)

cpdef list clarTransformation(Molecule mol, list ring)
