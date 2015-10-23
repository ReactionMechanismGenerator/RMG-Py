from .molecule cimport Atom

cpdef find_butadiene(Atom start, Atom end)

cpdef find_butadiene_end_with_charge(Atom start)

cpdef find_allyl_end_with_charge(Atom start)

cpdef list add_unsaturated_bonds(list path)

cpdef list add_allyls(list path)

cpdef list add_inverse_allyls(list path)