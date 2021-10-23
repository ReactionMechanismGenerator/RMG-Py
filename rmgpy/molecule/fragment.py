import os
from urllib.parse import quote
import itertools

import rmgpy.molecule.group as gr
import rmgpy.molecule.element as elements
import rmgpy.molecule.converter as converter
import rmgpy.molecule.resonance as resonance
from rmgpy.molecule.element import get_element
from rmgpy.molecule.graph import Graph, Vertex
from rmgpy.molecule.molecule import Atom, Bond, Molecule
from rmgpy.molecule.atomtype import get_atomtype, AtomTypeError
from rmgpy.molecule.kekulize import kekulize
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

class CuttingLabel(Vertex):

    def __init__(self, name='', label='', id=-1):
        Vertex.__init__(self)
        self.name = name # equivalent to Atom element symbol
        self.label = label # equivalent to Atom label attribute
        self.charge = 0
        self.radical_electrons = 0
        self.lone_pairs = 0
        self.isotope = -1
        self.id = id
        self.mass = 0

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
    def symbol(self): return self.name

    @property
    def bonds(self): return self.edges

    def is_specific_case_of(self, other):
        """
        Return ``True`` if `self` is a specific case of `other`, or ``False``
        otherwise. At this moment, this is the same as the :meth:`equivalent()`.
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
        return c

    def is_carbon(self):
        return False

    def is_nitrogen(self):
        return False

    def is_oxygen(self):
        return False

    def is_fluorine(self):
        return False

    def is_surface_site(self):
        return False

    def is_silicon(self):
        return False

    def is_sulfur(self):
        return False

    def is_chlorine(self):
        return False

    def is_iodine(self):
        return False

    def is_nos(self):
        """
        Return ``True`` if the atom represent either nitrogen, sulfur, or oxygen
        ``False`` if it does not.
        """
        return False

    def is_non_hydrogen(self):
        """
        Return ``True`` if the atom does not represent a hydrogen atom or
        ``False`` if it does.
        """
        return True

class Fragment(Graph):

    def __init__(self,
                label='',
                species_repr=None,
                vertices=None,
                symmetry=-1,
                multiplicity=-187,
                reactive=True,
                props=None,
                inchi='',
                smiles=''):
        self.index = -1
        self.label = label
        self.species_repr = species_repr
        Graph.__init__(self, vertices)
        self.symmetry_number = symmetry
        self._fingerprint = None
        self._inchi = None
        self._smiles = None
        self.props = props or {}
        self.multiplicity = multiplicity
        self.reactive = reactive

        if inchi and smiles:
            logging.warning('Both InChI and SMILES provided for Fragment instantiation, '
                            'using InChI and ignoring SMILES.')
        if inchi:
            self.from_inchi(inchi)
            self._inchi = inchi
        elif smiles:
            self.from_smiles_like_string(smiles)
            self._smiles = smiles

    def __deepcopy__(self, memo):
        return self.copy(deep=True)

    def __str__(self):
        """
        Return a string representation, in the form 'label(id)'.
        """
        if self.index == -1: return self.label
        else: return '{0}({1:d})'.format(self.label, self.index)

    def __getAtoms(self): return self.vertices
    def __setAtoms(self, atoms): self.vertices = atoms
    atoms = property(__getAtoms, __setAtoms)

    def draw(self, path):
        """
        Generate a pictorial representation of the chemical graph using the
        :mod:`draw` module. Use `path` to specify the file to save
        the generated image to; the image type is automatically determined by
        extension. Valid extensions are ``.png``, ``.svg``, ``.pdf``, and
        ``.ps``; of these, the first is a raster format and the remainder are
        vector formats.
        """
        from rmgpy.molecule.draw import MoleculeDrawer
        format = os.path.splitext(path)[-1][1:].lower()
        MoleculeDrawer().draw(self, format, target=path)

    def _repr_png_(self):
        """
        Return a png picture of the molecule, useful for ipython-qtconsole.
        """
        from rmgpy.molecule.draw import MoleculeDrawer
        tempFileName = 'temp_molecule.png'
        MoleculeDrawer().draw(self, 'png', tempFileName)
        png = open(tempFileName,'rb').read()
        os.unlink(tempFileName)
        return png

    def copy(self, deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """
        g = Graph.copy(self, deep)
        other = Fragment(vertices=g.vertices)
        return other

    def clear_labeled_atoms(self):
        """
        Remove the labels from all atoms in the fragment.
        """
        for v in self.vertices:
            v.label = ''

    def contains_labeled_atom(self, label):
        """
        Return :data:`True` if the fragment contains an atom with the label
        `label` and :data:`False` otherwise.
        """
        for v in self.vertices:
            if v.label == label: return True
        return False

    def get_labeled_atoms(self, label):
        """
        Return the atoms in the fragment that are labeled.
        """
        alist = [v for v in self.vertices if v.label == label]
        if alist == []:
            raise ValueError(
                'No vertex in the fragment \n{1}\n has the label "{0}".'.format(label, self.to_adjacency_list()))
        return alist

    def get_all_labeled_atoms(self):
        """
        Return the labeled atoms as a ``dict`` with the keys being the labels
        and the values the atoms themselves. If two or more atoms have the
        same label, the value is converted to a list of these atoms.
        """
        labeled = {}
        for v in self.vertices:
            if v.label != '':
                if v.label in labeled:
                    if isinstance(labeled[v.label],list):
                        labeled[v.label].append(v)
                    else:
                        labeled[v.label] = [labeled[v.label]]
                        labeled[v.label].append(v)
                else:
                    labeled[v.label] = v
        return labeled

    @property
    def fingerprint(self):
        """
        Fingerprint used to accelerate graph isomorphism comparisons with
        other molecules. The fingerprint is a short string containing a
        summary of selected information about the molecule. Two fingerprint
        strings matching is a necessary (but not sufficient) condition for
        the associated molecules to be isomorphic.

        Use an expanded molecular formula to also enable sorting.
        """
        if self._fingerprint is None:
            # Include these elements in this order at minimum
            element_dict = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0}
            all_elements = sorted(self.get_element_count().items(), key=lambda x: x[0])  # Sort alphabetically
            element_dict.update(all_elements)
            self._fingerprint = ''.join([f'{symbol}{num:0>2}' for symbol, num in element_dict.items()])
        return self._fingerprint

    @fingerprint.setter
    def fingerprint(self, fingerprint):
        self._fingerprint = fingerprint

    @property
    def inchi(self):
        """InChI string for this molecule. Read-only."""
        if self._inchi is None:
            self._inchi = self.to_inchi()
        return self._inchi

    @property
    def smiles(self):
        """SMILES string for this molecule. Read-only."""
        if self._smiles is None:
            self._smiles = self.to_smiles()
        return self._smiles

    def add_atom(self, atom):
        """
        Add an `atom` to the graph. The atom is initialized with no bonds.
        """
        self._fingerprint = self._inchi = self._smiles = None
        return self.add_vertex(atom)

    def remove_atom(self, atom):
        """
        Remove `atom` and all bonds associated with it from the graph. Does
        not remove atoms that no longer have any bonds as a result of this
        removal.
        """
        self._fingerprint = self._inchi = self._smiles = None
        return self.remove_vertex(atom)

    def contains_surface_site(self):
        """
        Returns ``True`` iff the molecule contains an 'X' surface site.
        """
        for atom in self.vertices:
            if atom.symbol == 'X':
                return True
        return False

    def is_surface_site(self):
        "Returns ``True`` iff the molecule is nothing but a surface site 'X'."
        return len(self.vertices) == 1 and self.vertices[0].is_surface_site()

    def has_bond(self, atom1, atom2):
        """
        Returns ``True`` if atoms `atom1` and `atom2` are connected
        by an bond, or ``False`` if not.
        """
        return self.has_edge(atom1, atom2)

    def get_bond(self, atom1, atom2):
        """
        Returns the bond connecting atoms `atom1` and `atom2`.
        """
        return self.get_edge(atom1, atom2)

    def add_bond(self, bond):
        """
        Add a `bond` to the graph as an edge connecting the two atoms `atom1`
        and `atom2`.
        """
        self._fingerprint = self._inchi = self._smiles = None
        return self.add_edge(bond)

    def remove_bond(self, bond):
        """
        Remove the bond between atoms `atom1` and `atom2` from the graph.
        Does not remove atoms that no longer have any bonds as a result of
        this removal.
        """
        self._fingerprint = self._inchi = self._smiles = None
        return self.remove_edge(bond)

    def get_net_charge(self):
        """
        Iterate through the atoms in the structure and calculate the net charge
        on the overall fragment.
        """
        charge = 0
        for v in self.vertices:
            charge += v.charge
        return charge

    def get_charge_span(self):
        """
        Iterate through the atoms in the structure and calculate the charge span
        on the overall molecule.
        The charge span is a measure of the number of charge separations in a molecule.
        """
        abs_net_charge = abs(self.get_net_charge())
        sum_of_abs_charges = sum([abs(atom.charge) for atom in self.vertices])
        return (sum_of_abs_charges - abs_net_charge) / 2

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

    def from_inchi(self, inchistr, backend='try-all'):
        """
        Convert an InChI string `inchistr` to a molecular structure.
        """
        translator.from_inchi(self, inchistr, backend)
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
        index_indicator = [x.span() for x in re.finditer(r'(\w?[LR][^:()]?)', smiles)]
        possible_cutting_label_list = re.findall(r'(\w?[LR][^:()]?)', smiles)

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
            elif len(strs)==2:
                # it's possible to be L+digit, L+C, C+L, R+a
                # check if it is a metal, if yes then don't add to cutting_label_list
                if strs in all_element_list:
                    # do not add it to cutting label list
                    add = False
                else:
                    add = True
                    # keep digit and remove the other non-metalic elements such as C
                    if strs[0] in ['L','R'] and strs[1].isdigit() == True:
                        # keep strs as it is
                        strs = strs
                        new_index_ranger = index_indicator[i]
                    elif strs[0] in ['L','R'] and strs[1].isdigit() != True:
                        strs = strs[0]
                        # keep the first index but subtract 1 for the end index
                        ind_tup = index_indicator[i]
                        int_ind = ind_tup[0]
                        end_ind = ind_tup[1]-1
                        new_index_ranger = (int_ind,end_ind)
                    else:
                        strs = strs[1]
                        # add 1 for the start index and keep the end index
                        ind_tup = index_indicator[i]
                        int_ind = ind_tup[0]+1
                        end_ind = ind_tup[1]
                        new_index_ranger = (int_ind,end_ind)
            elif len(strs)==3:
                # it's possible to be C+R+digit, C+L+i(metal), C+R+a(metal)
                # only C+R+digit has cutting label
                if strs[2].isdigit() == True:
                    add = True
                    strs = strs.replace(strs[0],"")
                    # add 1 for the start index and keep the end index
                    ind_tup = index_indicator[i]
                    int_ind = ind_tup[0]+1
                    end_ind = ind_tup[1]
                    new_index_ranger = (int_ind,end_ind)
                else:
                    # do not add this element to cutting_label_list
                    add = False
            if add == True:
                cutting_label_list.append(strs)
                ind_ranger.append(new_index_ranger)
        return ind_ranger, cutting_label_list

    def replace_cutting_label(self, smiles, ind_ranger, cutting_label_list, smiles_replace_dict):

        last_end_ind = 0
        new_smi = ""

        for ind, label_str in enumerate(cutting_label_list):
            tup = ind_ranger[ind]
            int_ind = tup[0]
            end_ind = tup[1]

            element = smiles_replace_dict[label_str]

            if ind == len(cutting_label_list)-1:
                new_smi = new_smi + smiles[last_end_ind:int_ind] + element + smiles[end_ind:]
            else:
                new_smi = new_smi + smiles[last_end_ind:int_ind] + element

            last_end_ind = end_ind
        # if the smiles does not include cutting label
        if new_smi == "":
            return smiles
        return new_smi

    def is_isomorphic(self, other, initial_map=None, generate_initial_map=False, save_order=False, strict=True):
        """
        Returns :data:`True` if two graphs are isomorphic and :data:`False`
        otherwise. The `initial_map` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Graph` object, or a :class:`TypeError` is raised.
        Also ensures multiplicities are also equal.

        Args:
            initial_map (dict, optional):         initial atom mapping to use
            generate_initial_map (bool, optional): if ``True``, initialize map by pairing atoms with same labels
            save_order (bool, optional):          if ``True``, reset atom order after performing atom isomorphism
            strict (bool, optional):             if ``False``, perform isomorphism ignoring electrons
        """
        # It only makes sense to compare a Molecule to a Molecule for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Graph):
            raise TypeError('Got a {0} object for parameter "other", when a Molecule object is required.'.format(other.__class__))
        # Do the quick isomorphism comparison using the fingerprint
        # Two fingerprint strings matching is a necessary (but not
        # sufficient!) condition for the associated molecules to be isomorphic
        if self.fingerprint != other.fingerprint:
            return False
        # check multiplicity
        if self.multiplicity != other.multiplicity:
            return False

        if generate_initial_map:
            initial_map = dict()
            for atom in self.vertices:
                if atom.label and atom.label != '':
                    for a in other.vertices:
                        if a.label == atom.label:
                            initial_map[atom] = a
                            break
                    else:
                        return False
            if not self.is_mapping_valid(other, initial_map, equivalent=True):
                return False

        # Do the full isomorphism comparison
        result = Graph.is_isomorphic(self, other, initial_map, save_order=save_order, strict=strict)
        return result

    def is_subgraph_isomorphic(self, other, initial_map=None, generate_initial_map=False, save_order=False):
        """
        Fragment's subgraph isomorphism check is done by first creating 
        a representative molecule of fragment, and then following same procedure
        of subgraph isomorphism check of `Molecule` object aganist `Group` object
        """
        if not isinstance(other, gr.Group):
            raise TypeError('Got a {0} object for parameter "other", when a Molecule object is required.'.format(other.__class__))
        group = other

        mapping = self.assign_representative_molecule()

        # Check multiplicity
        if group.multiplicity:
            if self.mol_repr.multiplicity not in group.multiplicity: return False

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
                if atom.label and atom.label != '':
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
                    for i,key in enumerate(keys):
                        initial_map[key] = atmlist[i]
                    if self.is_mapping_valid(other, initial_map, equivalent=False) and \
                            Graph.is_subgraph_isomorphic(self, other, initial_map, save_order=save_order):
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

    def assign_representative_molecule(self):

        # create a molecule from fragment.vertices.copy
        mapping = self.copy_and_map()

        # replace CuttingLabel with CC
        atoms = []
        additional_atoms = []
        additional_bonds = []
        for vertex in self.vertices:

            mapped_vertex = mapping[vertex]
            if isinstance(mapped_vertex, CuttingLabel):

                # replace cutting label with atom C
                atom_C1 = Atom(element=get_element('C'), 
                            radical_electrons=0, 
                            charge=0, 
                            lone_pairs=0)

                for bondedAtom, bond in mapped_vertex.edges.items():
                    new_bond = Bond(bondedAtom, atom_C1, order=bond.order)
                    
                    bondedAtom.edges[atom_C1] = new_bond
                    del bondedAtom.edges[mapped_vertex]

                    atom_C1.edges[bondedAtom] = new_bond

                # add hydrogens and carbon to make it CC
                atom_H1 = Atom(element=get_element('H'), 
                            radical_electrons=0, 
                            charge=0, 
                            lone_pairs=0)

                atom_H2 = Atom(element=get_element('H'), 
                            radical_electrons=0, 
                            charge=0, 
                            lone_pairs=0)

                atom_C2 = Atom(element=get_element('C'), 
                            radical_electrons=0, 
                            charge=0, 
                            lone_pairs=0)

                atom_H3 = Atom(element=get_element('H'), 
                            radical_electrons=0, 
                            charge=0, 
                            lone_pairs=0)

                atom_H4 = Atom(element=get_element('H'), 
                            radical_electrons=0, 
                            charge=0, 
                            lone_pairs=0)

                atom_C3 = Atom(element=get_element('C'),
                            radical_electrons=0,
                            charge=0,
                            lone_pairs=0)

                atom_H5 = Atom(element=get_element('H'),
                            radical_electrons=0,
                            charge=0, 
                            lone_pairs=0)

                atom_H6 = Atom(element=get_element('H'),
                            radical_electrons=0,
                            charge=0,
                            lone_pairs=0)

                atom_C4 = Atom(element=get_element('C'),
                            radical_electrons=0,
                            charge=0,
                            lone_pairs=0)

                atom_C5 = Atom(element=get_element('C'),
                            radical_electrons=0,
                            charge=0,
                            lone_pairs=0)

                atom_H7 = Atom(element=get_element('H'),
                               radical_electrons=0,
                               charge=0,
                               lone_pairs=0)

                atom_H8 = Atom(element=get_element('H'),
                               radical_electrons=0,
                               charge=0,
                               lone_pairs=0)

                atom_H9 = Atom(element=get_element('H'),
                               radical_electrons=0,
                               charge=0,
                               lone_pairs=0)

                atom_C6 = Atom(element=get_element('C'),
                            radical_electrons=0,
                            charge=0,
                            lone_pairs=0)

                atom_H10 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_C7 = Atom(element=get_element('C'),
                            radical_electrons=0,
                            charge=0,
                            lone_pairs=0)

                atom_H11 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_H12 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_C8 = Atom(element=get_element('C'),
                               radical_electrons=0,
                               charge=0,
                               lone_pairs=0)

                atom_H13 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_C9 = Atom(element=get_element('C'),
                               radical_electrons=0,
                               charge=0,
                               lone_pairs=0)

                atom_H14 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_H15 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_H16 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_C10 = Atom(element=get_element('C'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_H17 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_C11 = Atom(element=get_element('C'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_H18 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_H19 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_H20 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_C12 = Atom(element=get_element('C'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_H21 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_C13 = Atom(element=get_element('C'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_H22 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_C14 = Atom(element=get_element('C'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_H23 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_H24 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atom_H25 = Atom(element=get_element('H'),
                                radical_electrons=0,
                                charge=0,
                                lone_pairs=0)

                atoms.append(atom_C1)

                additional_atoms.extend([atom_H1,
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
                                    ])

                additional_bonds.extend([Bond(atom_C1, atom_H1, 1),
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
                                    ])

            else:
                atoms.append(mapped_vertex)

        mol_repr = Molecule()
        mol_repr.atoms = atoms
        for atom in additional_atoms: mol_repr.add_atom(atom)
        for bond in additional_bonds: mol_repr.add_bond(bond)
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

    def get_molecular_weight(self):
        """
        Return the fragmental weight of the fragment in kg/mol.
        """
        mass = 0
        for vertex in self.vertices:
            if isinstance(vertex, Atom):
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
            if isinstance(v, Atom):
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
        #Because we use lonepairs to match atomtypes and default is -100 when unspecified,
        #we should update before getting the atomtype.
        self.update_lone_pairs()

        for v in self.vertices:
            if not isinstance(v, Atom): continue
            try:
                v.atomtype = get_atomtype(v, v.edges)
            except AtomTypeError:
                if log_species:
                    logging.error("Could not update atomtypes for this fragment:\n{0}".format(self.to_adjacency_list()))
                if raise_exception:
                    raise
                v.atomtype = ATOMTYPES['R']

    def update_multiplicity(self):
        """
        Update the multiplicity of a newly formed molecule.
        """
        # Assume this is always true
        # There are cases where 2 radicalElectrons is a singlet, but
        # the triplet is often more stable, 
        self.multiplicity = self.get_radical_count() + 1


    def update_lone_pairs(self):
        """
        Iterate through the atoms in the structure and calculate the
        number of lone electron pairs, assuming a neutral molecule.
        """
        for v in self.vertices:
            if not isinstance(v, Atom): continue
            if not v.is_hydrogen():
                order = v.get_total_bond_order()
                v.lone_pairs = (elements.PeriodicSystem.valence_electrons[v.symbol] - v.radical_electrons - v.charge - int(order)) / 2.0
                if v.lone_pairs % 1 > 0 or v.lone_pairs > 4:
                    logging.error("Unable to determine the number of lone pairs for element {0} in {1}".format(v,self))
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
        formula = ''
        
        # Carbon and hydrogen always come first if carbon is present
        if 'C' in elements.keys():
            count = elements['C']
            formula += 'C{0:d}'.format(count) if count > 1 else 'C'
            del elements['C']
            if 'H' in elements.keys():
                count = elements['H']
                formula += 'H{0:d}'.format(count) if count > 1 else 'H'
                del elements['H']

        # Other atoms are in alphabetical order
        # (This includes hydrogen if carbon is not present)
        keys = sorted(elements.keys())
        for key in keys:
            count = elements[key]
            formula += '{0}{1:d}'.format(key, count) if count > 1 else key
        
        return formula

    def get_representative_molecule(self, mode='minimal', update=True):

        if mode == 'minimal':
            # create a molecule from fragment.vertices.copy
            mapping = self.copy_and_map()

            # replace CuttingLabel with H
            atoms = []
            for vertex in self.vertices:

                mapped_vertex = mapping[vertex]
                if isinstance(mapped_vertex, CuttingLabel):

                    # replace cutting label with atom H
                    atom_H = Atom(element=get_element('H'), 
                                radical_electrons=0, 
                                charge=0, 
                                lone_pairs=0)

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

    def get_aromatic_rings(self, rings=None):
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
            rdkitmol, rdAtomIndices = self.to_rdkit_mol(remove_h=False, return_mapping=True)
        except ValueError:
            logging.warning('Unable to check aromaticity by converting to RDKit Mol.')
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
                    for atom2 in ring0[i + 1:]:
                        if self.has_bond(atom1, atom2):
                            if rdkitmol.GetBondBetweenAtoms(rdAtomIndices[atom1],
                                                            rdAtomIndices[atom2]).GetBondType() is AROMATIC:
                                aromatic_bonds_in_ring.append(self.get_bond(atom1, atom2))
                else:  # didn't break so all atoms are carbon
                    if len(aromatic_bonds_in_ring) == 6:
                        aromatic_rings.append(ring0)
                        aromatic_bonds.append(aromatic_bonds_in_ring)

            return aromatic_rings, aromatic_bonds

    def is_aromatic(self):
        """ 
        Returns ``True`` if the fragment is aromatic, or ``False`` if not.  
        """
        mol0, _ = self.get_representative_molecule('minimal')   
        return mol0.is_aromatic()

    def atom_ids_valid(self):
        """
        Checks to see if the atom IDs are valid in this structure
        """
        num_atoms = len(self.atoms)
        num_IDs = len(set([atom.id for atom in self.atoms]))

        if num_atoms == num_IDs:
            # all are unique
            return True
        return False

    def kekulize(self):
        """
        Kekulizes an aromatic molecule.
        """
        kekulize(self)

    def assign_atom_ids(self):
        """
        Assigns an index to every atom in the fragment for tracking purposes.
        Uses entire range of cython's integer values to reduce chance of duplicates
        """

        global atom_id_counter

        for atom in self.atoms:
            atom.id = atom_id_counter
            atom_id_counter += 1
            if atom_id_counter == 2**15:
                atom_id_counter = -2**15

    def generate_resonance_structures(self, keep_isomorphic=False, filter_structures=True, save_order=False):
        """Returns a list of resonance structures of the fragment."""
        return resonance.generate_resonance_structures(self, keep_isomorphic=keep_isomorphic,
                                                        filter_structures = filter_structures,
                                                        save_order=save_order,
                                                        )

    def is_identical(self, other, strict=True):
        """
        Performs isomorphism checking, with the added constraint that atom IDs must match.

        Primary use case is tracking atoms in reactions for reaction degeneracy determination.

        Returns :data:`True` if two graphs are identical and :data:`False` otherwise.
        """

        if not isinstance(other, (Fragment, Molecule)):
            raise TypeError('Got a {0} object for parameter "other", when a Fragment object is required.'.format(other.__class__))

        # Get a set of atom indices for each molecule
        atom_ids = set([atom.id for atom in self.atoms])
        other_ids = set([atom.id for atom in other.atoms])

        if atom_ids == other_ids:
            # If the two molecules have the same indices, then they might be identical
            # Sort the atoms by ID
            atom_list = sorted(self.atoms, key=lambda x: x.id)
            other_list = sorted(other.atoms, key=lambda x: x.id)

            # If matching atom indices gives a valid mapping, then the molecules are fully identical
            mapping = {}
            for atom1, atom2 in zip(atom_list, other_list):
                mapping[atom1] = atom2

            return self.is_mapping_valid(other, mapping, equivalent=True, strict=strict)
        else:
            # The molecules don't have the same set of indices, so they are not identical
            return False

    def get_element_count(self):
        """
        Returns the element count for the fragment as a dictionary.
        """
        element_count = {}
        for atom in self.vertices:
            if not isinstance(atom, Atom): continue
            symbol = atom.element.symbol
            isotope = atom.element.isotope
            key = symbol # if isotope == -1 else (symbol, isotope)
            if key in element_count:
                element_count[key] += 1
            else:
                element_count[key] = 1

        return element_count

    def get_url(self):
        """
        Get a URL to the fragment's info page on the RMG website.
        """

        base_url = "http://rmg.mit.edu/database/molecule/"
        adjlist = self.to_adjacency_list(remove_h=False)
        url = base_url + quote(adjlist)
        return url.strip('_')

    def is_linear(self):
        """
        Return :data:`True` if the structure is linear and :data:`False`
        otherwise.
        """

        atom_count = len(self.vertices)

        # Monatomic molecules are definitely nonlinear
        if atom_count == 1:
            return False
        # Diatomic molecules are definitely linear
        elif atom_count == 2:
            return True
        # Cyclic molecules are definitely nonlinear
        elif self.is_cyclic():
            return False

        # True if all bonds are double bonds (e.g. O=C=O)
        all_double_bonds = True
        for atom1 in self.vertices:
            for bond in atom1.edges.values():
                if not bond.is_double(): all_double_bonds = False
        if all_double_bonds: return True

        # True if alternating single-triple bonds (e.g. H-C#C-H)
        # This test requires explicit hydrogen atoms
        for atom in self.vertices:
            bonds = list(atom.edges.values())
            if len(bonds)==1:
                continue # ok, next atom
            if len(bonds)>2:
                break # fail!
            if bonds[0].is_single() and bonds[1].is_triple():
                continue # ok, next atom
            if bonds[1].is_single() and bonds[0].is_triple():
                continue # ok, next atom
            break # fail if we haven't continued
        else:
            # didn't fail
            return True
        
        # not returned yet? must be nonlinear
        return False

    def saturate_radicals(self):
        """
        Saturate the fragment by replacing all radicals with bonds to hydrogen atoms.  Changes self molecule object.  
        """
        added = {}
        for atom in self.vertices:
            for i in range(atom.radical_electrons):
                H = Atom('H', radical_electrons=0, lone_pairs=0, charge=0)
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

    def is_aryl_radical(self, aromatic_rings=None):
        """
        Return ``True`` if the fragment only contains aryl radicals,
        ie. radical on an aromatic ring, or ``False`` otherwise.
        """
        if aromatic_rings is None:
            aromatic_rings = self.get_aromatic_rings()[0]

        total = self.get_radical_count()
        aromatic_atoms = set([atom for atom in itertools.chain.from_iterable(aromatic_rings)])
        aryl = sum([atom.radical_electrons for atom in aromatic_atoms])

        return total == aryl

    def get_num_atoms(self, element = None):
        """
        Return the number of atoms in molecule.  If element is given, ie. "H" or "C",
        the number of atoms of that element is returned.
        """
        num_atoms = 0
        if element == None:
            for vertex in self.vertices:
                if isinstance(vertex, Atom):
                    num_atoms += 1
        else:
            for vertex in self.vertices:
                if isinstance(vertex, Atom):
                    if vertex.element.symbol == element:
                        num_atoms += 1
        return num_atoms

# this variable is used to name atom IDs so that there are as few conflicts by 
# using the entire space of integer objects
atom_id_counter = -2**15
