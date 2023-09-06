#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import itertools
import logging

import rmgpy.molecule.group as gr
import rmgpy.molecule.element as elements
from rmgpy.molecule.element import Element
import rmgpy.molecule.converter as converter
from rmgpy.molecule.element import get_element
from rmgpy.molecule.graph import Graph, Vertex
from rmgpy.molecule.molecule import Atom, Bond, Molecule
from rmgpy.molecule.atomtype import get_atomtype, AtomTypeError, ATOMTYPES, AtomType
from rdkit import Chem

# this variable is used to name atom IDs so that there are as few conflicts by
# using the entire space of integer objects
ATOM_ID_COUNTER = -(2**15)


class CuttingLabel(Atom):
    def __init__(self, name="", label="", id=-1):
        super().__init__(
            element=Element(0, name, "cutting label", 0.0, -1),
            radical_electrons=0,
            charge=0,
            label=label,
            lone_pairs=0,
            site="",
            coords=None,
            id=id,
            props=None,
        )
        # define a new custom atom type
        self.atomtype = AtomType(label=label, generic=[], specific=[])
        # accessed by other methods in Fragment
        self.name = name
        self.isotope = -1
        
    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return str(self.name)

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "<CuttingLabel '{0}'>".format(str(self))

    @property
    def symbol(self):
        return self.name

    def is_specific_case_of(self, other):
        """
        Return ``True`` if `self` is a specific case of `other`, or ``False``
        otherwise. At this moment, this is the same as the :math:`equivalent()`.
        """
        return self.equivalent(other)

    def equivalent(self, other, strict=True):
        """
        Return ``True`` if `other` is indistinguishable from this CuttingLabel, or
        ``False`` otherwise. If `other` is an :class:`CuttingLabel` object, then all
        attributes must match exactly.
        """
        if isinstance(other, CuttingLabel):
            return self.name == other.name
        else:
            return False

    def copy(self):
        """
        Generate a deep copy of the current CuttingLabel. Modifying the
        attributes of the copy will not affect the original.
        """

        c = CuttingLabel.__new__(CuttingLabel)
        c.edges = {}
        c.reset_connectivity_values()
        c.name = self.name
        c.label = self.label
        c.charge = self.charge
        c.radical_electrons = self.radical_electrons
        c.lone_pairs = self.lone_pairs
        c.isotope = self.isotope
        c.id = self.id
        c.element = self.element
        c.site = ""
        c.morphology = ""
        return c


class Fragment(Molecule):
    def __init__(
        self,
        label="",
        species_repr=None,
        vertices=None,
        symmetry=-1,
        multiplicity=-187,
        reactive=True,
        props=None,
        inchi="",
        smiles="",
    ):
        if inchi and smiles:
            logging.warning(
                "Both InChI and SMILES provided for Fragment instantiation, "
                "using InChI and ignoring SMILES."
            )
            smiles = ""

        super().__init__(
            atoms=vertices,
            symmetry=symmetry,
            multiplicity=multiplicity,
            reactive=reactive,
            props=props,
            inchi=inchi,
            smiles=smiles,
            metal="",
            facet="",
        )
        self.index = -1
        self.label = label
        self.species_repr = species_repr

    def __str__(self):
        """
        Return a string representation, in the form 'label(id)'.
        """
        if self.index == -1:
            return self.label
        else:
            return "{0}({1:d})".format(self.label, self.index)

    # override methods
    def copy(self, deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """
        g = Graph.copy(self, deep)
        other = Fragment(vertices=g.vertices)
        # Copy connectivity values and sorting labels
        for i in range(len(self.vertices)):
            v1 = self.vertices[i]
            v2 = other.vertices[i]
            v2.connectivity1 = v1.connectivity1
            v2.connectivity2 = v1.connectivity2
            v2.connectivity3 = v1.connectivity3
            v2.sorting_label = v1.sorting_label
        other.multiplicity = self.multiplicity
        other.reactive = self.reactive
        return other

    def contains_surface_site(self):
        """
        Returns ``True`` iff the molecule contains an 'X' surface site.
        """
        for atom in self.vertices:
            if atom.symbol == "X":
                return True
        return False

    def is_surface_site(self):
        "Returns ``True`` iff the molecule is nothing but a surface site 'X'."
        return len(self.vertices) == 1 and self.vertices[0].is_surface_site()
    
    def merge(self, other):
        """
        Merge two fragments so as to store them in a single :class:`Fragment`
        object. The merged :class:`Fragment` object is returned.
        """
        g = Graph.merge(self, other)
        fragment = Fragment(vertices=g.vertices)
        return fragment

    def split(self):
        """
        Convert a single :class:`Fragment` object containing two or more
        unconnected fragment into separate class:`Fragment` objects.
        """
        graphs = Graph.split(self)
        fragments = []
        for g in graphs:
            fragment = Fragment(vertices=g.vertices)
            fragments.append(fragment)
        return fragments

    def get_singlet_carbene_count(self):
        """
        Return the total number of singlet carbenes (lone pair on a carbon atom)
        in the fragment. Counts the number of carbon atoms with a lone pair.
        In the case of [C] with two lone pairs, this method will return 1.
        """
        carbenes = 0
        for v in self.vertices:
            if isinstance(v, Atom) and v.is_carbon() and v.lone_pairs > 0:
                carbenes += 1
        return carbenes

    def get_radical_count(self):
        """
        Return the total number of radical electrons on all atoms in the
        molecule. In this function, monoradical atoms count as one, biradicals
        count as two, etc.
        """
        radicals = 0
        for v in self.vertices:
            radicals += v.radical_electrons
        return radicals

    def is_subgraph_isomorphic(
        self, other, initial_map=None, generate_initial_map=False, save_order=False
    ):
        """
        Fragment's subgraph isomorphism check is done by first creating
        a representative molecule of fragment, and then following same procedure
        of subgraph isomorphism check of `Molecule` object aganist `Group` object
        """
        if not isinstance(other, gr.Group):
            raise TypeError(
                'Got a {0} object for parameter "other", when a Molecule object is required.'.format(
                    other.__class__
                )
            )
        group = other

        mapping = self.assign_representative_molecule()

        # Check multiplicity
        if group.multiplicity:
            if self.mol_repr.multiplicity not in group.multiplicity:
                return False

        # Compare radical counts
        if self.mol_repr.get_radical_count() < group.radicalCount:
            return False

        # Compare element counts
        element_count = self.mol_repr.get_element_count()
        for element, count in group.elementCount.items():
            if element not in element_count:
                return False
            elif element_count[element] < count:
                return False

        if generate_initial_map:
            keys = []
            atms = []
            initial_map = dict()
            for atom in self.atoms:
                if atom.label and atom.label != "":
                    atom_label_map = [a for a in other.atoms if a.label == atom.label]
                    if not atom_label_map:
                        return False
                    elif len(atom_label_map) == 1:
                        initial_map[atom] = atom_label_map[0]
                    else:
                        keys.append(atom)
                        atms.append(atom_label_map)
            if atms:
                for atmlist in itertools.product(*atms):
                    # skip entries that map multiple graph atoms to the same subgraph atom
                    if len(set(atmlist)) != len(atmlist):
                        continue
                    for i, key in enumerate(keys):
                        initial_map[key] = atmlist[i]
                    if self.is_mapping_valid(
                        other, initial_map, equivalent=False
                    ) and Graph.is_subgraph_isomorphic(
                        self, other, initial_map, save_order=save_order
                    ):
                        return True
                else:
                    return False
            else:
                if not self.is_mapping_valid(other, initial_map, equivalent=False):
                    return False

        # Do the isomorphism comparison
        new_initial_map = None
        if initial_map:
            new_initial_map = {}
            for fragment_vertex in initial_map:
                repr_mol_vertex = mapping[fragment_vertex]
                new_initial_map[repr_mol_vertex] = initial_map[fragment_vertex]

        result = Graph.is_subgraph_isomorphic(self.mol_repr, other, new_initial_map)
        return result

    def get_molecular_weight(self):
        """
        Return the fragmental weight of the fragment in kg/mol.
        """
        mass = 0
        for vertex in self.vertices:
            mass += vertex.element.mass
        return mass

    def calculate_cp0(self):
        """
        Return the value of the heat capacity at zero temperature in J/mol*K.
        """
        self.assign_representative_molecule()
        # currently linear fragment will be treated as non-linear molecule
        return self.mol_repr.calculate_cp0()

    def calculate_cpinf(self):
        """
        Return the value of the heat capacity at infinite temperature in J/mol*K.
        """
        self.assign_representative_molecule()
        # currently linear fragment will be treated as non-linear molecule
        return self.mol_repr.calculate_cpinf()

    def calculate_symmetry_number(self):
        """
        Return the symmetry number for the structure. The symmetry number
        includes both external and internal modes. First replace Cuttinglabel
        with different elements and then calculate symmetry number
        """
        from rmgpy.molecule.symmetry import calculate_symmetry_number

        smiles = self.to_smiles()

        _, cutting_label_list = self.detect_cutting_label(smiles)

        metal_list = [
            "[Cl]",
            "[I]",
            "[Si]",
            "[F]",
            "[Si+]",
            "[Si-]",
            "[Br]",
            "[He+]",
            "[Ne+]",
            "[Ar+]",
            "[He-]",
            "[Ne-]",
            "[Ar-]",
            "[P]",
            "[P+]",
            "[P-]",
        ]

        for index, element in enumerate(cutting_label_list):
            smiles = smiles.replace(element, metal_list[index], 1)

        frag_sym = Molecule().from_smiles(smiles)

        frag_sym.update_connectivity_values()  # for consistent results
        self.symmetry_number = calculate_symmetry_number(frag_sym)
        return self.symmetry_number

    def is_radical(self):
        """
        Return ``True`` if the fragment contains at least one radical electron,
        or ``False`` otherwise.
        """
        for vertex in self.vertices:
            if isinstance(vertex, Atom) and vertex.radical_electrons > 0:
                return True
        return False

    def update(self, sort_atoms=True):
        # currently sort_atoms does not work for fragments
        for v in self.vertices:
            if not isinstance(v, CuttingLabel):
                v.update_charge()

        self.update_atomtypes()
        self.update_multiplicity()
        self.sort_vertices()

    def update_atomtypes(self, log_species=True, raise_exception=True):
        """
        Iterate through the atoms in the structure, checking their atom types
        to ensure they are correct (i.e. accurately describe their local bond
        environment) and complete (i.e. are as detailed as possible).

        If `raise_exception` is `False`, then the generic atomType 'R' will
        be prescribed to any atom when get_atomtype fails. Currently used for
        resonance hybrid atom types.
        """
        # Because we use lonepairs to match atomtypes and default is -100 when unspecified,
        # we should update before getting the atomtype.
        self.update_lone_pairs()

        for v in self.vertices:
            if isinstance(v, CuttingLabel):
                continue
            try:
                v.atomtype = get_atomtype(v, v.edges)
            except AtomTypeError:
                if log_species:
                    logging.error(
                        "Could not update atomtypes for this fragment:\n{0}".format(
                            self.to_adjacency_list()
                        )
                    )
                if raise_exception:
                    raise
                v.atomtype = ATOMTYPES["R"]

    def update_lone_pairs(self):
        """
        Iterate through the atoms in the structure and calculate the
        number of lone electron pairs, assuming a neutral molecule.
        """
        for v in self.vertices:
            if isinstance(v, CuttingLabel):
                continue
            if not v.is_hydrogen():
                order = v.get_total_bond_order()
                v.lone_pairs = (
                    elements.PeriodicSystem.valence_electrons[v.symbol]
                    - v.radical_electrons
                    - v.charge
                    - int(order)
                ) / 2.0
                if v.lone_pairs % 1 > 0 or v.lone_pairs > 4:
                    logging.error(
                        "Unable to determine the number of lone pairs for element {0} in {1}".format(
                            v, self
                        )
                    )
            else:
                v.lone_pairs = 0

    def get_formula(self):
        """
        Return the molecular formula for the fragment.
        """

        # Count the number of each element in the molecule
        elements = {}
        cuttinglabels = {}
        for atom in self.vertices:
            symbol = atom.symbol
            elements[symbol] = elements.get(symbol, 0) + 1

        # Use the Hill system to generate the formula
        formula = ""

        # Carbon and hydrogen always come first if carbon is present
        if "C" in elements.keys():
            count = elements["C"]
            formula += "C{0:d}".format(count) if count > 1 else "C"
            del elements["C"]
            if "H" in elements.keys():
                count = elements["H"]
                formula += "H{0:d}".format(count) if count > 1 else "H"
                del elements["H"]

        # Other atoms are in alphabetical order
        # (This includes hydrogen if carbon is not present)
        keys = sorted(elements.keys())
        for key in keys:
            count = elements[key]
            formula += "{0}{1:d}".format(key, count) if count > 1 else key

        return formula

    def to_rdkit_mol(self, remove_h=False, return_mapping=True, save_order=False):
        """
        Convert a molecular structure to a RDKit rdmol object.
        """
        if remove_h:
            # because we're replacing
            # cutting labels with hydrogens
            # so do not allow removeHs to be True
            raise "Currently fragment to_rdkit_mol only allows keeping all the hydrogens."

        mol0, mapping = self.get_representative_molecule("minimal", update=False)

        rdmol, rdAtomIdx_mol0 = converter.to_rdkit_mol(
            mol0,
            remove_h=remove_h,
            return_mapping=return_mapping,
            sanitize=True,
            save_order=save_order,
        )

        rdAtomIdx_frag = {}
        for frag_atom, mol0_atom in mapping.items():
            rd_idx = rdAtomIdx_mol0[mol0_atom]
            rdAtomIdx_frag[frag_atom] = rd_idx

        # sync the order of fragment vertices with the order
        # of mol0.atoms since mol0.atoms is changed/sorted in
        # converter.to_rdkit_mol().
        # Since the rdmol's atoms order is same as the order of mol0's atoms,
        # the synchronization between fragment.atoms order and mol0.atoms order
        # is necessary to make sure the order of fragment vertices
        # reflects the order of rdmol's atoms
        vertices_order = []
        for v in self.vertices:
            a = mapping[v]
            idx = mol0.atoms.index(a)
            vertices_order.append((v, idx))

        adapted_vertices = [
            tup[0] for tup in sorted(vertices_order, key=lambda tup: tup[1])
        ]

        self.vertices = adapted_vertices

        return rdmol, rdAtomIdx_frag

    def from_adjacency_list(
        self,
        adjlist,
        saturate_h=False,
        raise_atomtype_exception=True,
        raise_charge_exception=True,
    ):
        """
        Convert a string adjacency list `adjlist` to a fragment structure.
        Skips the first line (assuming it's a label) unless `withLabel` is
        ``False``.
        """
        from rmgpy.molecule.adjlist import from_adjacency_list

        (
            self.vertices,
            self.multiplicity,
            self.site,
            self.morphology,
        ) = from_adjacency_list(adjlist, group=False, saturate_h=saturate_h)
        self.update_atomtypes(raise_exception=raise_atomtype_exception)

        # Check if multiplicity is possible
        n_rad = self.get_radical_count()
        multiplicity = self.multiplicity
        if not (
            n_rad + 1 == multiplicity
            or n_rad - 1 == multiplicity
            or n_rad - 3 == multiplicity
            or n_rad - 5 == multiplicity
        ):
            raise ValueError(
                "Impossible multiplicity for molecule\n{0}\n multiplicity = {1} and number of "
                "unpaired electrons = {2}".format(
                    self.to_adjacency_list(), multiplicity, n_rad
                )
            )
        if raise_charge_exception:
            if self.get_net_charge() != 0:
                raise ValueError(
                    "Non-neutral molecule encountered. "
                    "Currently, AFM does not support ion chemistry.\n {0}".format(
                        adjlist
                    )
                )
        return self

    def get_aromatic_rings(self, rings=None, save_order=False):
        """
        Returns all aromatic rings as a list of atoms and a list of bonds.

        Identifies rings using `Graph.get_smallest_set_of_smallest_rings()`, then uses RDKit to perceive aromaticity.
        RDKit uses an atom-based pi-electron counting algorithm to check aromaticity based on Huckel's Rule.
        Therefore, this method identifies "true" aromaticity, rather than simply the RMG bond type.

        The method currently restricts aromaticity to six-membered carbon-only rings. This is a limitation imposed
        by RMG, and not by RDKit.
        """

        from rdkit.Chem.rdchem import BondType

        AROMATIC = BondType.AROMATIC

        if rings is None:
            rings = self.get_relevant_cycles()
            rings = [ring for ring in rings if len(ring) == 6]
        if not rings:
            return [], []

        try:
            rdkitmol, rdAtomIndices = self.to_rdkit_mol(
                remove_h=False, return_mapping=True, save_order=save_order
            )
        except ValueError:
            logging.warning("Unable to check aromaticity by converting to RDKit Mol.")
        else:
            aromatic_rings = []
            aromatic_bonds = []
            for ring0 in rings:
                aromatic_bonds_in_ring = []
                # Figure out which atoms and bonds are aromatic and reassign appropriately:
                for i, atom1 in enumerate(ring0):
                    if not atom1.is_carbon():
                        # all atoms in the ring must be carbon in RMG for our definition of aromatic
                        break
                    for atom2 in ring0[i + 1 :]:
                        if self.has_bond(atom1, atom2):
                            if (
                                rdkitmol.GetBondBetweenAtoms(
                                    rdAtomIndices[atom1], rdAtomIndices[atom2]
                                ).GetBondType()
                                is AROMATIC
                            ):
                                aromatic_bonds_in_ring.append(
                                    self.get_bond(atom1, atom2)
                                )
                else:  # didn't break so all atoms are carbon
                    if len(aromatic_bonds_in_ring) == 6:
                        aromatic_rings.append(ring0)
                        aromatic_bonds.append(aromatic_bonds_in_ring)

            return aromatic_rings, aromatic_bonds

    def is_aromatic(self):
        """
        Returns ``True`` if the fragment is aromatic, or ``False`` if not.
        """
        mol0, _ = self.get_representative_molecule("minimal")
        return mol0.is_aromatic()

    def to_smiles(self):
        cutting_label_list = []
        for vertex in self.vertices:
            if isinstance(vertex, CuttingLabel):
                cutting_label_list.append(vertex.symbol)

        smiles_before = self.copy(deep=True)
        final_vertices = []
        for ind, atom in enumerate(smiles_before.atoms):
            element_symbol = atom.symbol
            if isinstance(atom, CuttingLabel):
                substi = Atom(
                    element=get_element("Si"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=3,
                )
                substi.label = element_symbol

                for bonded_atom, bond in atom.edges.items():
                    new_bond = Bond(bonded_atom, substi, order=bond.order)

                    bonded_atom.edges[substi] = new_bond
                    del bonded_atom.edges[atom]

                    substi.edges[bonded_atom] = new_bond

                final_vertices.append(substi)
            else:
                final_vertices.append(atom)

        smiles_before.vertices = final_vertices
        mol_repr = Molecule()
        mol_repr.atoms = smiles_before.vertices
        mol_repr.update()
        smiles_after = mol_repr.to_smiles()
        import re

        smiles = re.sub(r"\[Si-3\]", "", smiles_after)

        return smiles

    def get_element_count(self):
        """
        Returns the element count for the fragment as a dictionary.
        """
        element_count = {}
        for atom in self.vertices:
            if isinstance(atom, CuttingLabel):
                continue
            symbol = atom.element.symbol
            isotope = atom.element.isotope
            key = symbol  # if isotope == -1 else (symbol, isotope)
            if key in element_count:
                element_count[key] += 1
            else:
                element_count[key] = 1

        return element_count

    def saturate_radicals(self):
        """
        Saturate the fragment by replacing all radicals with bonds to hydrogen atoms.  Changes self molecule object.
        """
        added = {}
        for atom in self.vertices:
            for i in range(atom.radical_electrons):
                H = Atom("H", radical_electrons=0, lone_pairs=0, charge=0)
                bond = Bond(atom, H, 1)
                self.add_atom(H)
                self.add_bond(bond)
                if atom not in added:
                    added[atom] = []
                added[atom].append([H, bond])
                atom.decrement_radical()

        # Update the atom types of the saturated structure (not sure why
        # this is necessary, because saturating with H shouldn't be
        # changing atom types, but it doesn't hurt anything and is not
        # very expensive, so will do it anyway)
        self.sort_vertices()
        self.update_atomtypes()
        self.multiplicity = 1

        return added

    def is_aryl_radical(self, aromatic_rings=None, save_order=False):
        """
        Return ``True`` if the molecule only contains aryl radicals,
        i.e., radical on an aromatic ring, or ``False`` otherwise.
        If no ``aromatic_rings`` provided, aromatic rings will be searched in-place,
        and this process may involve atom order change by default. Set ``save_order`` to
        ``True`` to force the atom order unchanged.
        """
        if aromatic_rings is None:
            aromatic_rings = self.get_aromatic_rings(save_order=save_order)[0]

        total = self.get_radical_count()
        aromatic_atoms = set(
            [atom for atom in itertools.chain.from_iterable(aromatic_rings)]
        )
        aryl = sum([atom.radical_electrons for atom in aromatic_atoms])

        return total == aryl

    def get_num_atoms(self, element=None):
        """
        Return the number of atoms in molecule.  If element is given, ie. "H" or "C",
        the number of atoms of that element is returned.
        """
        num_atoms = 0
        if element is None:
            for vertex in self.vertices:
                if not isinstance(vertex, CuttingLabel):
                    num_atoms += 1
        else:
            for vertex in self.vertices:
                if not isinstance(vertex, CuttingLabel):
                    if vertex.element.symbol == element:
                        num_atoms += 1
        return num_atoms

    # extension methods
    def from_rdkit_mol(self, rdkitmol, atom_replace_dict=None):
        """
        Convert a RDKit Mol object `rdkitmol` to a molecular structure. Uses
        `RDKit <http://rdkit.org/>`_ to perform the conversion.
        This Kekulizes everything, removing all aromatic atom types.
        """

        from rdkit import Chem

        self.vertices = []

        # Add hydrogen atoms to complete molecule if needed
        rdkitmol.UpdatePropertyCache(strict=False)
        rdkitmol = Chem.AddHs(rdkitmol)
        Chem.rdmolops.Kekulize(rdkitmol, clearAromaticFlags=True)

        # iterate through atoms in rdkitmol
        for i in range(rdkitmol.GetNumAtoms()):
            rdkitatom = rdkitmol.GetAtomWithIdx(i)

            # Use atomic number as key for element
            number = rdkitatom.GetAtomicNum()
            element = get_element(number)

            # Process charge
            charge = rdkitatom.GetFormalCharge()
            radical_electrons = rdkitatom.GetNumRadicalElectrons()

            ELE = element.symbol
            if "[" + ELE + "]" in atom_replace_dict:
                cutting_label_name = atom_replace_dict["[" + ELE + "]"]
                cutting_label = CuttingLabel(name=cutting_label_name)
                self.vertices.append(cutting_label)
            else:
                atom = Atom(element, radical_electrons, charge, "", 0)
                self.vertices.append(atom)

            # Add bonds by iterating again through atoms
            for j in range(0, i):
                rdkitbond = rdkitmol.GetBondBetweenAtoms(i, j)
                if rdkitbond is not None:
                    order = 0

                    # Process bond type
                    rdbondtype = rdkitbond.GetBondType()
                    if rdbondtype.name == "SINGLE":
                        order = 1
                    elif rdbondtype.name == "DOUBLE":
                        order = 2
                    elif rdbondtype.name == "TRIPLE":
                        order = 3
                    elif rdbondtype.name == "AROMATIC":
                        order = 1.5

                    bond = Bond(self.vertices[i], self.vertices[j], order)
                    self.add_bond(bond)

        # We need to update lone pairs first because the charge was set by RDKit
        self.update_lone_pairs()
        # Set atom types and connectivity values
        self.update()

        # Assume this is always true
        # There are cases where 2 radicalElectrons is a singlet, but
        # the triplet is often more stable,
        self.update_multiplicity()
        # mol.update_atomtypes()

        return self

    def cut_molecule(self, output_smiles=False, cut_through=True, size_threshold=None):
        """
        For given input, output a list of cut fragments (either string or Fragment).
        if output_smiles = True, the output list of fragments will be smiles.
        cut_through is whether to have multiple cuts (True) or just cut once (False).
        For initial input species cutting, it's better to have multiple cutting so that each
        fragment lies in appropriate size range. As for auto-cutting during model generation,
        the number of products should be less (or equal) than 3, due to the 3 body collision
        rxn format limit.
        size_threshold is the minimum size for each aliphatic fragment size. Default value is 5.
        """
        mol = self.generate_resonance_structures()[0]

        # slice mol
        frag_smiles_list = []
        if cut_through:
            arom_cut_frag = self.sliceitup_arom(
                mol.to_smiles(), size_threshold=size_threshold
            )
            for frag in arom_cut_frag:
                aliph_cut_frag = self.sliceitup_aliph(
                    frag, size_threshold=size_threshold
                )
                for ele in aliph_cut_frag:
                    frag_smiles_list.append(ele)
        else:
            # if aromatic, only perform sliceitup_arom, if aliphatic, only sliceitup_aliph
            if mol.is_aromatic():
                # try aromatic cut first, if no cut found, try aliphatic cut then
                frag_smiles_list = self.sliceitup_arom(
                    mol.to_smiles(), size_threshold=size_threshold
                )
                if len(frag_smiles_list) == 1:
                    # try aliphatic cut then
                    frag_smiles_list = self.sliceitup_aliph(
                        mol.to_smiles(), size_threshold=size_threshold
                    )
            else:
                frag_smiles_list = self.sliceitup_aliph(
                    mol.to_smiles(), size_threshold=size_threshold
                )

        if output_smiles:
            return frag_smiles_list
        else:
            frag_list = []
            for frag in frag_smiles_list:
                frag = Fragment().from_smiles_like_string(frag)
                res_frag = frag.generate_resonance_structures()[0]
                frag_list.append(res_frag)
            return frag_list

    def sliceitup_arom(self, molecule, size_threshold=None):
        """
        Several specified aromatic patterns
        """
        # set min size for each aliphatic fragment size
        if size_threshold:
            size_threshold = size_threshold
        else:
            size_threshold = 5
        # if input is smiles string, output smiles
        if isinstance(molecule, str):
            molecule_smiles = molecule
        # if input is fragment, output frafgment
        elif isinstance(molecule, Fragment):
            mol = molecule.generate_resonance_structures()[0]
            molecule_smiles = mol.to_smiles()

        # if input has Cuttinglabel, need to transform to later
        _, cutting_label_list = self.detect_cutting_label(molecule_smiles)
        # transfer to rdkit molecule for substruct matching
        f = self.from_smiles_like_string(molecule_smiles)
        molecule_to_cut, rdAtomIdx_frag = f.to_rdkit_mol()

        # replace CuttingLabel to special Atom (metal) in rdkit
        for atom, idx in rdAtomIdx_frag.items():
            if isinstance(atom, CuttingLabel):
                cuttinglabel_atom = molecule_to_cut.GetAtomWithIdx(idx)
                if atom.symbol == "R":
                    cuttinglabel_atom.SetAtomicNum(
                        11
                    )  # [Na], will replace back to CuttingLabel later
                else:
                    cuttinglabel_atom.SetAtomicNum(19)  # [K]

        # substructure matching
        pattern_list = ["pattern_1", "pattern_2", "pattern_3", "pattern_4"]
        frag_list = []
        for pattern in pattern_list:
            emol, atom_map_index = self.pattern_call("Arom", pattern)
            # start pattern matching
            atom_map = molecule_to_cut.GetSubstructMatches(emol)
            if atom_map:
                # go through all matches and see if it agree with size threshold
                for matched_atom_map in atom_map:
                    # find the correct bond to break
                    bonds_to_cut = []
                    for ind in atom_map_index:
                        b1 = matched_atom_map[ind]
                        b2 = matched_atom_map[ind + 1]
                        bond = molecule_to_cut.GetBondBetweenAtoms(b1, b2)
                        bonds_to_cut.append(bond)
                    # Break bonds
                    newmol = Chem.RWMol(molecule_to_cut)
                    # fragmentize
                    new_mol = Chem.FragmentOnBonds(
                        newmol,
                        [bond.GetIdx() for bond in bonds_to_cut],
                        dummyLabels=[(0, 0)] * len(bonds_to_cut),
                    )
                    # mol_set contains new set of fragments
                    mol_set = Chem.GetMolFrags(new_mol, asMols=True)
                    # check all fragments' size
                    if all(
                        sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
                        >= size_threshold
                        for mol in mol_set
                    ):
                        # replace * at cutting position with cutting label
                        for ind, rdmol in enumerate(mol_set):
                            frag = Chem.MolToSmiles(rdmol)
                            if len(mol_set) > 2:  # means it cut into 3 fragments
                                if frag.count("*") > 1:
                                    # replace both with R
                                    frag_smi = frag.replace("*", "R")
                                else:
                                    frag_smi = frag.replace("*", "L")
                            else:  # means it only cut once, generate 2 fragments
                                if ind == 0:
                                    frag_smi = frag.replace("*", "R")
                                else:
                                    frag_smi = frag.replace("*", "L")
                            frag_list.append(frag_smi)
                        break
                    else:
                        # turn to next matched_atom_map
                        continue
            else:
                # no match for this pattern
                # cannot cut this molecule with this pattern
                continue
            # if there is appropariate cutting, stop further pattern matching
            if frag_list != []:
                break

        if frag_list == []:
            # this means no appropriate match for all patterns
            # cannot cut this molecule
            frag_list.append(molecule_smiles)
        else:
            if cutting_label_list != []:
                frag_list_replaced = []
                # replace metal atom back to Cuttinglabel
                for metal_frag in frag_list:
                    n_frag_smiles = metal_frag.replace("[Na]", "R")
                    nn_frag_smiles = n_frag_smiles.replace("[K]", "L")
                    frag_list_replaced.append(nn_frag_smiles)
                frag_list = frag_list_replaced

        if isinstance(molecule, str):
            frag_list_new = []
            for frag in frag_list:
                n_frag_smiles = frag.replace("[Na]", "R")
                new_frag_smiles = n_frag_smiles.replace("[K]", "L")
                frag_list_new.append(new_frag_smiles)
            return frag_list_new
        elif isinstance(molecule, Fragment):
            frag_list_new = []
            for frag in frag_list:
                n_frag_smiles = frag.replace("[Na]", "R")
                new_frag_smiles = n_frag_smiles.replace("[K]", "L")

                frag = Fragment().from_smiles_like_string(new_frag_smiles)
                res_frag = frag.generate_resonance_structures()[0]
                frag_list_new.append(res_frag)
            return frag_list_new

    def sliceitup_aliph(self, molecule, size_threshold=None):
        """
        Several specified aliphatic patterns
        """
        # set min size for each aliphatic fragment size
        if size_threshold:
            size_threshold = size_threshold
        else:
            size_threshold = 5
        # if input is smiles string, output smiles
        if isinstance(molecule, str):
            molecule_smiles = molecule
        # if input is fragment, output frafgment
        elif isinstance(molecule, Fragment):
            mol = molecule.generate_resonance_structures()[0]
            molecule_smiles = mol.to_smiles()

        # if input has Cuttinglabel, need to transform to later
        _, cutting_label_list = self.detect_cutting_label(molecule_smiles)
        # transfer to rdkit molecule for substruct matching
        f = self.from_smiles_like_string(molecule_smiles)
        molecule_to_cut, rdAtomIdx_frag = f.to_rdkit_mol()

        # replace CuttingLabel to special Atom (metal) in rdkit
        for atom, idx in rdAtomIdx_frag.items():
            if isinstance(atom, CuttingLabel):
                cuttinglabel_atom = molecule_to_cut.GetAtomWithIdx(idx)
                if atom.symbol == "R":
                    cuttinglabel_atom.SetAtomicNum(
                        11
                    )  # [Na], will replace back to CuttingLabel later
                else:
                    cuttinglabel_atom.SetAtomicNum(19)  # [K]

        # substructure matching
        pattern_list = ["pattern_1", "pattern_2", "pattern_3"]
        frag_list = []
        for pattern in pattern_list:
            emol, atom_map_index = self.pattern_call("Aliph", pattern)
            # start pattern matching
            atom_map = molecule_to_cut.GetSubstructMatches(emol)
            if atom_map:
                # go through all matches and see if it agree with size threshold or in ring
                for matched_atom_map in atom_map:
                    if self.check_in_ring(molecule_to_cut, matched_atom_map):
                        # do not cut ring
                        continue
                    # find the correct bond to break
                    bonds_to_cut = []
                    for ind in atom_map_index:
                        b1 = matched_atom_map[ind]
                        b2 = matched_atom_map[ind + 1]
                        bond = molecule_to_cut.GetBondBetweenAtoms(b1, b2)
                        bonds_to_cut.append(bond)
                    # Break bonds
                    newmol = Chem.RWMol(molecule_to_cut)
                    # fragmentize
                    new_mol = Chem.FragmentOnBonds(
                        newmol,
                        [bond.GetIdx() for bond in bonds_to_cut],
                        dummyLabels=[(0, 0)] * len(bonds_to_cut),
                    )
                    # mol_set contains new set of fragments
                    mol_set = Chem.GetMolFrags(new_mol, asMols=True)
                    # check all fragments' size
                    if all(
                        sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
                        >= size_threshold
                        for mol in mol_set
                    ):
                        # replace * at cutting position with cutting label
                        for ind, rdmol in enumerate(mol_set):
                            frag = Chem.MolToSmiles(rdmol)
                            if len(mol_set) > 2:  # means it cut into 3 fragments
                                if frag.count("*") > 1:
                                    # replace both with R
                                    frag_smi = frag.replace("*", "R")
                                else:
                                    frag_smi = frag.replace("*", "L")
                            else:  # means it only cut once, generate 2 fragments
                                if ind == 0:
                                    frag_smi = frag.replace("*", "R")
                                else:
                                    frag_smi = frag.replace("*", "L")
                            frag_list.append(frag_smi)
                        break
                    else:
                        # turn to next matched_atom_map
                        continue
            else:
                # no match for this pattern
                # cannot cut this molecule with this pattern
                continue
            # if there is appropariate cutting, stop further pattern matching
            if frag_list != []:
                break

        if frag_list == []:
            # this means no appropriate match for all patterns
            # cannot cut this molecule
            frag_list.append(molecule_smiles)
        else:
            if cutting_label_list != []:
                frag_list_replaced = []
                # replace metal atom back to Cuttinglabel
                for metal_frag in frag_list:
                    n_frag_smiles = metal_frag.replace("[Na]", "R")
                    nn_frag_smiles = n_frag_smiles.replace("[K]", "L")
                    frag_list_replaced.append(nn_frag_smiles)
                frag_list = frag_list_replaced

        if isinstance(molecule, str):
            frag_list_new = []
            for frag in frag_list:
                n_frag_smiles = frag.replace("[Na]", "R")
                new_frag_smiles = n_frag_smiles.replace("[K]", "L")
                frag_list_new.append(new_frag_smiles)
            return frag_list_new
        elif isinstance(molecule, Fragment):
            frag_list_new = []
            for frag in frag_list:
                n_frag_smiles = frag.replace("[Na]", "R")
                new_frag_smiles = n_frag_smiles.replace("[K]", "L")

                frag = Fragment().from_smiles_like_string(new_frag_smiles)
                res_frag = frag.generate_resonance_structures()[0]
                frag_list_new.append(res_frag)
            return frag_list_new

    def pattern_call(self, pattern_type, pattern):
        """
        pattern_type currently only supports 'Aliph' and 'Arom', specifying which pattern will
        give output as pattern in rdkit mol and atom_map_index for substructure matching
        """
        if pattern_type == "Arom":
            if pattern == "pattern_1":
                mol = Chem.MolFromSmiles("c1ccccc1C(C(c2ccccc2)CCCCC)CCCCC")
                emol = Chem.RWMol(mol)
                # add some H at terminal aliphatic C to avoid cutting at potential allylic C
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(16, 24, Chem.rdchem.BondType.SINGLE)
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(16, 25, Chem.rdchem.BondType.SINGLE)
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(21, 26, Chem.rdchem.BondType.SINGLE)
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(21, 27, Chem.rdchem.BondType.SINGLE)
                atom_map_index = [16, 21]

            if pattern == "pattern_2":
                mol = Chem.MolFromSmiles("c1ccccc1C(CCCCC)CCCCC")
                emol = Chem.RWMol(mol)
                # add some H at terminal aliphatic C to avoid cutting at potential allylic C
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(14, 17, Chem.rdchem.BondType.SINGLE)
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(14, 18, Chem.rdchem.BondType.SINGLE)
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(9, 19, Chem.rdchem.BondType.SINGLE)
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(9, 20, Chem.rdchem.BondType.SINGLE)
                atom_map_index = [9, 14]

            if pattern == "pattern_3":
                mol = Chem.MolFromSmiles("c1ccccc1C(=C)CCCCCC")
                emol = Chem.RWMol(mol)
                # add some H at terminal aliphatic C to avoid cutting at potential allylic C
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(11, 14, Chem.rdchem.BondType.SINGLE)
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(11, 15, Chem.rdchem.BondType.SINGLE)
                atom_map_index = [11]

            if pattern == "pattern_4":
                mol = Chem.MolFromSmiles("c1ccccc1CCCCCC")
                emol = Chem.RWMol(mol)
                # add some H at terminal aliphatic C to avoid cutting at potential allylic C
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(9, 12, Chem.rdchem.BondType.SINGLE)
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(9, 13, Chem.rdchem.BondType.SINGLE)
                atom_map_index = [9]

        elif pattern_type == "Aliph":
            if pattern == "pattern_1":
                mol = Chem.MolFromSmiles("CCC=CCCCCC")
                emol = Chem.RWMol(mol)
                # add some H at terminal aliphatic C to avoid cutting at potential allylic C
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(6, 9, Chem.rdchem.BondType.SINGLE)
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(6, 10, Chem.rdchem.BondType.SINGLE)
                atom_map_index = [6]

            if pattern == "pattern_2":
                mol = Chem.MolFromSmiles("C=CCCCCC")
                emol = Chem.RWMol(mol)
                # add some H at terminal aliphatic C to avoid cutting at potential allylic C
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(4, 7, Chem.rdchem.BondType.SINGLE)
                emol.AddAtom(Chem.rdchem.Atom("H"))
                emol.AddBond(4, 8, Chem.rdchem.BondType.SINGLE)
                atom_map_index = [4]

            if pattern == "pattern_3":
                mol = Chem.MolFromSmiles("CCCCCC")
                mol = Chem.AddHs(mol, onlyOnAtoms=[2, 3])
                emol = Chem.RWMol(mol)
                atom_map_index = [2]
        else:
            raise NameError("Currently only Arom or Aliph type patterns are included")

        return emol, atom_map_index

    def check_in_ring(self, rd_mol, mapped_atom_idx):
        """
        Check if the mapped structure is in ring
        """
        for idx in mapped_atom_idx:
            atom = rd_mol.GetAtomWithIdx(idx)
            if atom.IsInRing():
                return True
        return False

    def from_smiles_like_string(self, smiles_like_string):
        smiles = smiles_like_string

        # input: smiles
        # output: ind_ranger & cutting_label_list
        ind_ranger, cutting_label_list = self.detect_cutting_label(smiles)

        smiles_replace_dict = {}
        metal_list = [
            "[Na]",
            "[K]",
            "[Cs]",
            "[Fr]",
            "[Be]",
            "[Mg]",
            "[Ca]",
            "[Sr]",
            "[Ba]",
            "[Hf]",
            "[Nb]",
            "[Ta]",
            "[Db]",
            "[Mo]",
        ]
        for index, label_str in enumerate(cutting_label_list):
            smiles_replace_dict[label_str] = metal_list[index]

        atom_replace_dict = {}
        for key, value in smiles_replace_dict.items():
            atom_replace_dict[value] = key

        # replace cutting labels with elements in smiles
        # to generate rdkit compatible SMILES
        new_smi = self.replace_cutting_label(
            smiles, ind_ranger, cutting_label_list, smiles_replace_dict
        )

        from rdkit import Chem

        rdkitmol = Chem.MolFromSmiles(new_smi)

        self.from_rdkit_mol(rdkitmol, atom_replace_dict)

        return self

    def detect_cutting_label(self, smiles):
        import re
        from rmgpy.molecule.element import element_list

        # store elements' symbol
        all_element_list = []
        for element in element_list[1:]:
            all_element_list.append(element.symbol)

        # store the tuple of matched indexes, however,
        # the index might contain redundant elements such as C, Ra, (), Li, ...
        index_indicator = [x.span() for x in re.finditer(r"(\w?[LR][^:()]?)", smiles)]
        possible_cutting_label_list = re.findall(r"(\w?[LR][^:()]?)", smiles)

        cutting_label_list = []
        ind_ranger = []

        # check if the matched items are cutting labels indeed
        for i, strs in enumerate(possible_cutting_label_list):
            # initialize "add" for every possible cutting label
            add = False
            if len(strs) == 1:
                # it should be a cutting label either R or L
                # add it to cutting label list
                add = True
                # add the index span
                new_index_ranger = index_indicator[i]
            elif len(strs) == 2:
                # it's possible to be L+digit, L+C, C+L, R+a
                # check if it is a metal, if yes then don't add to cutting_label_list
                if strs in all_element_list:
                    # do not add it to cutting label list
                    add = False
                else:
                    add = True
                    # keep digit and remove the other non-metalic elements such as C
                    if strs[0] in ["L", "R"] and strs[1].isdigit():
                        # keep strs as it is
                        strs = strs
                        new_index_ranger = index_indicator[i]
                    elif strs[0] in ["L", "R"] and not strs[1].isdigit():
                        strs = strs[0]
                        # keep the first index but subtract 1 for the end index
                        ind_tup = index_indicator[i]
                        int_ind = ind_tup[0]
                        end_ind = ind_tup[1] - 1
                        new_index_ranger = (int_ind, end_ind)
                    else:
                        strs = strs[1]
                        # add 1 for the start index and keep the end index
                        ind_tup = index_indicator[i]
                        int_ind = ind_tup[0] + 1
                        end_ind = ind_tup[1]
                        new_index_ranger = (int_ind, end_ind)
            elif len(strs) == 3:
                # it's possible to be C+R+digit, C+L+i(metal), C+R+a(metal)
                # only C+R+digit has cutting label
                if strs[2].isdigit():
                    add = True
                    strs = strs.replace(strs[0], "")
                    # add 1 for the start index and keep the end index
                    ind_tup = index_indicator[i]
                    int_ind = ind_tup[0] + 1
                    end_ind = ind_tup[1]
                    new_index_ranger = (int_ind, end_ind)
                else:
                    # do not add this element to cutting_label_list
                    add = False
            if add:
                cutting_label_list.append(strs)
                ind_ranger.append(new_index_ranger)
        return ind_ranger, cutting_label_list

    def replace_cutting_label(
        self, smiles, ind_ranger, cutting_label_list, smiles_replace_dict
    ):
        last_end_ind = 0
        new_smi = ""

        for ind, label_str in enumerate(cutting_label_list):
            tup = ind_ranger[ind]
            int_ind = tup[0]
            end_ind = tup[1]

            element = smiles_replace_dict[label_str]

            if ind == len(cutting_label_list) - 1:
                new_smi = (
                    new_smi + smiles[last_end_ind:int_ind] + element + smiles[end_ind:]
                )
            else:
                new_smi = new_smi + smiles[last_end_ind:int_ind] + element

            last_end_ind = end_ind
        # if the smiles does not include cutting label
        if new_smi == "":
            return smiles
        return new_smi

    def assign_representative_molecule(self):
        # create a molecule from fragment.vertices.copy
        mapping = self.copy_and_map()

        # replace CuttingLabel with C14 structure to avoid fragment symmetry number issue
        atoms = []
        additional_atoms = []
        additional_bonds = []
        for vertex in self.vertices:
            mapped_vertex = mapping[vertex]
            if isinstance(mapped_vertex, CuttingLabel):
                # replace cutting label with atom C
                atom_C1 = Atom(
                    element=get_element("C"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                for bondedAtom, bond in mapped_vertex.edges.items():
                    new_bond = Bond(bondedAtom, atom_C1, order=bond.order)

                    bondedAtom.edges[atom_C1] = new_bond
                    del bondedAtom.edges[mapped_vertex]

                    atom_C1.edges[bondedAtom] = new_bond

                # add hydrogens and carbon to make it CC
                atom_H1 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H2 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_C2 = Atom(
                    element=get_element("C"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H3 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H4 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_C3 = Atom(
                    element=get_element("C"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H5 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H6 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_C4 = Atom(
                    element=get_element("C"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_C5 = Atom(
                    element=get_element("C"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H7 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H8 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H9 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_C6 = Atom(
                    element=get_element("C"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H10 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_C7 = Atom(
                    element=get_element("C"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H11 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H12 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_C8 = Atom(
                    element=get_element("C"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H13 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_C9 = Atom(
                    element=get_element("C"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H14 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H15 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H16 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_C10 = Atom(
                    element=get_element("C"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H17 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_C11 = Atom(
                    element=get_element("C"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H18 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H19 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H20 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_C12 = Atom(
                    element=get_element("C"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H21 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_C13 = Atom(
                    element=get_element("C"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H22 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_C14 = Atom(
                    element=get_element("C"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H23 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H24 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atom_H25 = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                atoms.append(atom_C1)

                additional_atoms.extend(
                    [
                        atom_H1,
                        atom_H2,
                        atom_H3,
                        atom_H4,
                        atom_H5,
                        atom_H6,
                        atom_H7,
                        atom_H8,
                        atom_H9,
                        atom_H10,
                        atom_H11,
                        atom_H12,
                        atom_H13,
                        atom_H14,
                        atom_H15,
                        atom_H16,
                        atom_H17,
                        atom_H18,
                        atom_H19,
                        atom_H20,
                        atom_H21,
                        atom_H22,
                        atom_H23,
                        atom_H24,
                        atom_H25,
                        atom_C2,
                        atom_C3,
                        atom_C4,
                        atom_C5,
                        atom_C6,
                        atom_C7,
                        atom_C8,
                        atom_C9,
                        atom_C10,
                        atom_C11,
                        atom_C12,
                        atom_C13,
                        atom_C14,
                    ]
                )

                additional_bonds.extend(
                    [
                        Bond(atom_C1, atom_H1, 1),
                        Bond(atom_C1, atom_H2, 1),
                        Bond(atom_C2, atom_H3, 1),
                        Bond(atom_C2, atom_H4, 1),
                        Bond(atom_C1, atom_C2, 1),
                        Bond(atom_C2, atom_C3, 1),
                        Bond(atom_C3, atom_H5, 1),
                        Bond(atom_C3, atom_H6, 1),
                        Bond(atom_C3, atom_C4, 1),
                        Bond(atom_C4, atom_C5, 1),
                        Bond(atom_C5, atom_H7, 1),
                        Bond(atom_C5, atom_H8, 1),
                        Bond(atom_C5, atom_H9, 1),
                        Bond(atom_C4, atom_C6, 1),
                        Bond(atom_C6, atom_C7, 2),
                        Bond(atom_C6, atom_H10, 1),
                        Bond(atom_C7, atom_H11, 1),
                        Bond(atom_C7, atom_H12, 1),
                        Bond(atom_C4, atom_C8, 1),
                        Bond(atom_C8, atom_H13, 1),
                        Bond(atom_C8, atom_C9, 1),
                        Bond(atom_C9, atom_H14, 1),
                        Bond(atom_C9, atom_H15, 1),
                        Bond(atom_C9, atom_H16, 1),
                        Bond(atom_C8, atom_C10, 1),
                        Bond(atom_C10, atom_H17, 1),
                        Bond(atom_C10, atom_C11, 1),
                        Bond(atom_C11, atom_H18, 1),
                        Bond(atom_C11, atom_H19, 1),
                        Bond(atom_C11, atom_H20, 1),
                        Bond(atom_C10, atom_C12, 1),
                        Bond(atom_C12, atom_H21, 1),
                        Bond(atom_C12, atom_C13, 2),
                        Bond(atom_C13, atom_H22, 1),
                        Bond(atom_C13, atom_C14, 1),
                        Bond(atom_C14, atom_H23, 1),
                        Bond(atom_C14, atom_H24, 1),
                        Bond(atom_C14, atom_H25, 1),
                    ]
                )

            else:
                atoms.append(mapped_vertex)

        mol_repr = Molecule()
        mol_repr.atoms = atoms
        for atom in additional_atoms:
            mol_repr.add_atom(atom)
        for bond in additional_bonds:
            mol_repr.add_bond(bond)
        # update connectivity
        mol_repr.update()

        # create a species object from molecule
        self.mol_repr = mol_repr

        return mapping

    def assign_representative_species(self):
        from rmgpy.species import Species

        self.assign_representative_molecule()
        self.species_repr = Species(molecule=[self.mol_repr])
        self.symmetry_number = self.get_symmetry_number()
        self.species_repr.symmetry_number = self.symmetry_number

    def get_representative_molecule(self, mode="minimal", update=True):
        if mode != "minimal":
            raise RuntimeError(
                'Fragment.get_representative_molecule onyl supports mode="minimal"'
            )
        # create a molecule from fragment.vertices.copy
        mapping = self.copy_and_map()

        # replace CuttingLabel with H
        atoms = []
        for vertex in self.vertices:
            mapped_vertex = mapping[vertex]
            if isinstance(mapped_vertex, CuttingLabel):
                # replace cutting label with atom H
                atom_H = Atom(
                    element=get_element("H"),
                    radical_electrons=0,
                    charge=0,
                    lone_pairs=0,
                )

                for bondedAtom, bond in mapped_vertex.edges.items():
                    new_bond = Bond(bondedAtom, atom_H, order=bond.order)

                    bondedAtom.edges[atom_H] = new_bond
                    del bondedAtom.edges[mapped_vertex]

                    atom_H.edges[bondedAtom] = new_bond

                mapping[vertex] = atom_H
                atoms.append(atom_H)

            else:
                atoms.append(mapped_vertex)

        # Note: mapping is a dict with
        # key: self.vertex and value: mol_repr.atom
        mol_repr = Molecule()
        mol_repr.atoms = atoms
        if update:
            mol_repr.update()

        return mol_repr, mapping
