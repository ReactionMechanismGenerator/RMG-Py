from .molecule cimport Atom, Bond, Molecule
from .graph cimport Vertex, Edge

cpdef list find_butadiene(Atom start, Atom end)

cpdef list find_butadiene_end_with_charge(Atom start)

cpdef list find_allyl_end_with_charge(Atom start)

cpdef list find_shortest_path(Vertex start, Vertex end, list path=*)

cpdef list add_unsaturated_bonds(list path)

cpdef list add_allyls(list path)

cpdef list add_inverse_allyls(list path)

cpdef dict compute_atom_distance(list atom_indices, Molecule mol)

cpdef list findAllDelocalizationPaths(Atom atom1)

cpdef list findAllDelocalizationPathsLonePairRadicalCharge(Atom atom1)

cpdef list findAllDelocalizationPathsLonePairRadicalDoubleBond(Atom atom1)

cpdef list findAllDelocalizationPathsN5dd_N5ts(Atom atom1)
