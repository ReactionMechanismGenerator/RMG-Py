from .graph cimport Vertex, Edge, Graph
from .molecule cimport Atom, Bond, Molecule

cpdef list generateResonanceIsomers(Molecule mol)

cpdef list generateAdjacentResonanceIsomers(Molecule mol)

cpdef list generateLonePairRadicalResonanceIsomers(Molecule mol)

cpdef list generateN5dd_N5tsResonanceIsomers(Molecule mol)

cpdef list generate_isomorphic_isomers(Molecule mol)

cpdef list generateKekulizedResonanceIsomers(Molecule mol)

cpdef tuple populate_resonance_generation_algorithm()
