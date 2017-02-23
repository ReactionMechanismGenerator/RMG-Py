from .graph cimport Vertex, Edge, Graph
from .molecule cimport Atom, Bond, Molecule

cpdef list populateResonanceAlgorithms(dict features=?)

cpdef dict analyzeMolecule(Molecule mol)

cpdef list generateResonanceStructures(Molecule mol)

cpdef list _generateResonanceStructures(list molList, list methodList, bint copy=?)

cpdef list generateAdjacentResonanceStructures(Molecule mol)

cpdef list generateLonePairRadicalResonanceStructures(Molecule mol)

cpdef list generateN5dd_N5tsResonanceStructures(Molecule mol)

cpdef list generateIsomorphicResonanceStructures(Molecule mol)

cpdef list generateAromaticResonanceStructures(Molecule mol, dict features=?)

cpdef list generateKekuleStructure(Molecule mol)

cpdef list generateOppositeKekuleStructure(Molecule mol)

cpdef list generateClarStructures(Molecule mol)

cpdef list _clarOptimization(Molecule mol, list constraints=?, maxNum=?)

cpdef list _clarTransformation(Molecule mol, list ring)
