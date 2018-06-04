from .graph cimport Vertex, Edge, Graph
from .molecule cimport Atom, Bond, Molecule

cpdef list populate_resonance_algorithms(dict features=?)

cpdef dict analyze_molecule(Molecule mol)

cpdef list generate_resonance_structures(Molecule mol, bint clarStructures=?, bint keepIsomorphic=?)

cpdef list _generate_resonance_structures(list molList, list methodList, bint keepIsomorphic=?, bint copy=?)

cpdef list generate_adjacent_resonance_structures(Molecule mol)

cpdef list generate_lone_pair_radical_resonance_structures(Molecule mol)

cpdef list generate_N5dd_N5ts_resonance_structures(Molecule mol)

cpdef list generate_isomorphic_resonance_structures(Molecule mol)

cpdef list generate_aromatic_resonance_structures(Molecule mol, dict features=?)

cpdef list generate_kekule_structure(Molecule mol)

cpdef list generate_opposite_kekule_structure(Molecule mol)

cpdef list generate_clar_structures(Molecule mol)

cpdef list _clar_optimization(Molecule mol, list constraints=?, maxNum=?)

cpdef list _clar_transformation(Molecule mol, list ring)
