from .graph cimport Vertex, Edge, Graph
from .molecule cimport Atom, Bond, Molecule

cpdef list generateResonanceIsomers(Molecule mol)

cpdef list getAdjacentResonanceIsomers(Molecule mol)

cpdef list getLonePairRadicalResonanceIsomers(Molecule mol)

cpdef list getN5dd_N5tsResonanceIsomers(Molecule mol)

cpdef list findAllDelocalizationPaths(Molecule mol, Atom atom1)

cpdef list findAllDelocalizationPathsLonePairRadical(Molecule mol, Atom atom1)

cpdef list findAllDelocalizationPathsN5dd_N5ts(Molecule mol, Atom atom1)

cpdef list generate_isomorphic_isomers(Molecule mol)

cpdef list getKekulizedResonanceIsomers(Molecule mol)
