import os
import urllib
import itertools

from rmgpy.species import Species
import rmgpy.molecule.group as gr
import rmgpy.molecule.element as elements
import rmgpy.molecule.converter as converter
import rmgpy.molecule.resonance as resonance
from rmgpy.molecule.element import getElement
from rmgpy.molecule.graph import Graph, Vertex
from rmgpy.molecule.molecule import Atom, Bond, Molecule
from rmgpy.molecule.atomtype import getAtomType, AtomTypeError
from rmgpy.molecule.kekulize import kekulize

class CuttingLabel(Vertex):

    def __init__(self, name='', label='', id=-1):
        Vertex.__init__(self)
        self.name = name # equivalent to Atom element symbol
        self.label = label # equivalent to Atom label attribute
        self.charge = 0
        self.radicalElectrons = 0
        self.lonePairs = 0
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

    def isSpecificCaseOf(self, other):
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
        c.resetConnectivityValues()
        c.name = self.name
        c.label = self.label
        c.charge = self.charge
        c.radicalElectrons = self.radicalElectrons
        c.lonePairs = self.lonePairs
        c.isotope = self.isotope
        c.id = self.id
        return c

    def isCarbon(self):
        return False

    def isNitrogen(self):
        return False

    def isOxygen(self):
        return False

    def isFluorine(self):
        return False

    def isSurfaceSite(self):
        return False

    def isSilicon(self):
        return False

    def isSulfur(self):
        return False

    def isChlorine(self):
        return False

    def isIodine(self):
        return False

    def isNOS(self):
        """
        Return ``True`` if the atom represent either nitrogen, sulfur, or oxygen
        ``False`` if it does not.
        """
        return False

class Fragment(Graph):

    def __init__(self,
                label='',
                species_repr=None,
                vertices=None,
                symmetry=-1,
                multiplicity=-187,
                reactive=True,
                props=None,
                InChI='',
                SMILES=''):
        self.index = -1
        self.label = label
        self.species_repr = species_repr
        Graph.__init__(self, vertices)
        self.symmetryNumber = symmetry
        self.fingerprint = None
        self.inchi = None
        self.smiles = None
        self.props = props or {}
        self.multiplicity = multiplicity
        self.reactive = reactive

        if InChI and SMILES:
            logging.warning('Both InChI and SMILES provided for Fragment instantiation, using InChI and ignoring SMILES.')
        if InChI:
            self.fromInChI(InChI)
            self.inchi = InChI
        elif SMILES:
            self.from_SMILES_like_string(SMILES)
            self.smiles = SMILES

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

    def clearLabeledAtoms(self):
        """
        Remove the labels from all atoms in the fragment.
        """
        for v in self.vertices:
            v.label = ''

    def containsLabeledAtom(self, label):
        """
        Return :data:`True` if the fragment contains an atom with the label
        `label` and :data:`False` otherwise.
        """
        for v in self.vertices:
            if v.label == label: return True
        return False

    def getLabeledAtom(self, label):
        """
        Return the atoms in the fragment that are labeled.
        """
        alist = [v for v in self.vertices if v.label == label]
        if alist == []:
            raise ValueError(
                'No vertex in the fragment \n{1}\n has the label "{0}".'.format(label, self.toAdjacencyList()))
        return alist

    def getLabeledAtoms(self):
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
    def InChI(self):
        """InChI string for this molecule. Read-only."""
        if self.inchi is None:
            self.inchi = self.toInChI()
        return self.inchi

    @property
    def SMILES(self):
        """SMILES string for this molecule. Read-only."""
        if self.smiles is None:
            self.smiles = self.toSMILES()
        return self.smiles

    def addAtom(self, atom):
        """
        Add an `atom` to the graph. The atom is initialized with no bonds.
        """
        self.fingerprint = None
        return self.addVertex(atom)

    def removeAtom(self, atom):
        """
        Remove `atom` and all bonds associated with it from the graph. Does
        not remove atoms that no longer have any bonds as a result of this
        removal.
        """
        self.fingerprint = None
        return self.removeVertex(atom)

    def containsSurfaceSite(self):
        """
        Returns ``True`` iff the molecule contains an 'X' surface site.
        """
        for atom in self.vertices:
            if atom.symbol == 'X':
                return True
        return False

    def isSurfaceSite(self):
        "Returns ``True`` iff the molecule is nothing but a surface site 'X'."
        return len(self.vertices) == 1 and self.vertices[0].isSurfaceSite()

    def hasBond(self, atom1, atom2):
        """
        Returns ``True`` if atoms `atom1` and `atom2` are connected
        by an bond, or ``False`` if not.
        """
        return self.hasEdge(atom1, atom2)

    def getBond(self, atom1, atom2):
        """
        Returns the bond connecting atoms `atom1` and `atom2`.
        """
        return self.getEdge(atom1, atom2)

    def addBond(self, bond):
        """
        Add a `bond` to the graph as an edge connecting the two atoms `atom1`
        and `atom2`.
        """
        self.fingerprint = None
        return self.addEdge(bond)

    def removeBond(self, bond):
        """
        Remove the bond between atoms `atom1` and `atom2` from the graph.
        Does not remove atoms that no longer have any bonds as a result of
        this removal.
        """
        self.fingerprint = None
        return self.removeEdge(bond)

    def getNetCharge(self):
        """
        Iterate through the atoms in the structure and calculate the net charge
        on the overall fragment.
        """
        charge = 0
        for v in self.vertices:
            charge += v.charge
        return charge

    def getChargeSpan(self):
        """
        Iterate through the atoms in the structure and calculate the charge span
        on the overall molecule.
        The charge span is a measure of the number of charge separations in a molecule.
        """
        abs_net_charge = abs(self.getNetCharge())
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

    def getSingletCarbeneCount(self):
        """
        Return the total number of singlet carbenes (lone pair on a carbon atom)
        in the fragment. Counts the number of carbon atoms with a lone pair.
        In the case of [C] with two lone pairs, this method will return 1.
        """
        carbenes = 0
        for v in self.vertices:
            if isinstance(v, Atom) and v.isCarbon() and v.lonePairs > 0:
                carbenes += 1
        return carbenes
    
    def getRadicalCount(self):
        """
        Return the total number of radical electrons on all atoms in the
        molecule. In this function, monoradical atoms count as one, biradicals
        count as two, etc.
        """
        radicals = 0
        for v in self.vertices:
            radicals += v.radicalElectrons
        return radicals

    def fromInChI(self, inchistr, backend='try-all'):
        """
        Convert an InChI string `inchistr` to a molecular structure.
        """
        translator.fromInChI(self, inchistr, backend)
        return self

    def from_SMILES_like_string(self, SMILES_like_string):

        import re
        smiles = SMILES_like_string

        cutting_label_list = re.findall(r'([LR]\d?)', smiles)

        smiles_replace_dict = {}
        metal_list = ['[Na]', '[K]', '[Cs]', '[Fr]', '[Be]', '[Mg]', '[Ca]', '[Sr]', '[Ba]']
        for index, label_str in enumerate(cutting_label_list):
            smiles_replace_dict[label_str] = metal_list[index]

        atom_replace_dict = {}
        for key, value in smiles_replace_dict.iteritems():
            atom_replace_dict[value] = key

        for label_str, element in smiles_replace_dict.iteritems():
            smiles = smiles.replace(label_str, element)

        from rdkit import Chem
        rdkitmol = Chem.MolFromSmiles(smiles)

        self.fromRDKitMol(rdkitmol, atom_replace_dict)

        return self

    def isIsomorphic(self, other, initialMap=None, generateInitialMap=False, saveOrder=False, strict=True):
        """
        Returns :data:`True` if two graphs are isomorphic and :data:`False`
        otherwise. The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Graph` object, or a :class:`TypeError` is raised.
        Also ensures multiplicities are also equal.

        Args:
            initialMap (dict, optional):         initial atom mapping to use
            generateInitialMap (bool, optional): if ``True``, initialize map by pairing atoms with same labels
            saveOrder (bool, optional):          if ``True``, reset atom order after performing atom isomorphism
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

        if generateInitialMap:
            initialMap = dict()
            for atom in self.vertices:
                if atom.label and atom.label != '':
                    for a in other.vertices:
                        if a.label == atom.label:
                            initialMap[atom] = a
                            break
                    else:
                        return False
            if not self.isMappingValid(other, initialMap, equivalent=True):
                return False

        # Do the full isomorphism comparison
        result = Graph.isIsomorphic(self, other, initialMap, saveOrder=saveOrder, strict=strict)
        return result

    def isSubgraphIsomorphic(self, other, initialMap=None, generateInitialMap=False, saveOrder=False):
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
        if self.mol_repr.getRadicalCount() < group.radicalCount:
            return False

        # Compare element counts
        element_count = self.mol_repr.get_element_count()
        for element, count in group.elementCount.iteritems():
            if element not in element_count:
                return False
            elif element_count[element] < count:
                return False

        if generateInitialMap:
            keys = []
            atms = []
            initialMap = dict()
            for atom in self.atoms:
                if atom.label and atom.label != '':
                    atom_label_map = [a for a in other.atoms if a.label == atom.label]
                    if not atom_label_map:
                        return False
                    elif len(atom_label_map) == 1:
                        initialMap[atom] = atom_label_map[0]
                    else:
                        keys.append(atom)
                        atms.append(atom_label_map)
            if atms:
                for atmlist in itertools.product(*atms):
                    # skip entries that map multiple graph atoms to the same subgraph atom
                    if len(set(atmlist)) != len(atmlist):
                        continue
                    for i,key in enumerate(keys):
                        initialMap[key] = atmlist[i]
                    if self.isMappingValid(other, initialMap, equivalent=False) and \
                            Graph.isSubgraphIsomorphic(self, other, initialMap, saveOrder=saveOrder):
                        return True
                else:
                    return False
            else:
                if not self.isMappingValid(other, initialMap, equivalent=False):
                    return False

        # Do the isomorphism comparison
        new_initial_map = None
        if initialMap:
            new_initial_map = {}
            for fragment_vertex in initialMap:
                repr_mol_vertex = mapping[fragment_vertex]
                new_initial_map[repr_mol_vertex] = initialMap[fragment_vertex]

        result = Graph.isSubgraphIsomorphic(self.mol_repr, other, new_initial_map)
        return result

    def assign_representative_molecule(self):

        # create a molecule from fragment.vertices.copy
        mapping = self.copyAndMap()

        # replace CuttingLabel with CC
        atoms = []
        additional_atoms = []
        additional_bonds = []
        for vertex in self.vertices:

            mapped_vertex = mapping[vertex]
            if isinstance(mapped_vertex, CuttingLabel):

                # replace cutting label with atom C
                atom_C1 = Atom(element=getElement('C'), 
                            radicalElectrons=0, 
                            charge=0, 
                            lonePairs=0)

                for bondedAtom, bond in mapped_vertex.edges.iteritems():
                    new_bond = Bond(bondedAtom, atom_C1, order=bond.order)
                    
                    bondedAtom.edges[atom_C1] = new_bond
                    del bondedAtom.edges[mapped_vertex]

                    atom_C1.edges[bondedAtom] = new_bond

                # add hydrogens and carbon to make it CC
                atom_H1 = Atom(element=getElement('H'), 
                            radicalElectrons=0, 
                            charge=0, 
                            lonePairs=0)

                atom_H2 = Atom(element=getElement('H'), 
                            radicalElectrons=0, 
                            charge=0, 
                            lonePairs=0)

                atom_C2 = Atom(element=getElement('C'), 
                            radicalElectrons=0, 
                            charge=0, 
                            lonePairs=0)

                atom_H3 = Atom(element=getElement('H'), 
                            radicalElectrons=0, 
                            charge=0, 
                            lonePairs=0)

                atom_H4 = Atom(element=getElement('H'), 
                            radicalElectrons=0, 
                            charge=0, 
                            lonePairs=0)

                atom_H5 = Atom(element=getElement('H'), 
                            radicalElectrons=0, 
                            charge=0, 
                            lonePairs=0)

                atoms.append(atom_C1)

                additional_atoms.extend([atom_H1,
                                    atom_H2,
                                    atom_H3,
                                    atom_H4,
                                    atom_H5,
                                    atom_C2])

                additional_bonds.extend([Bond(atom_C1, atom_H1, 1),
                                    Bond(atom_C1, atom_H2, 1),
                                    Bond(atom_C2, atom_H3, 1),
                                    Bond(atom_C2, atom_H4, 1),
                                    Bond(atom_C2, atom_H5, 1),
                                    Bond(atom_C1, atom_C2, 1)])

            else:
                atoms.append(mapped_vertex)

        mol_repr = Molecule()
        mol_repr.atoms = atoms
        for atom in additional_atoms: mol_repr.addAtom(atom)
        for bond in additional_bonds: mol_repr.addBond(bond)
        # update connectivity
        mol_repr.update()

        # create a species object from molecule
        self.mol_repr = mol_repr

        return mapping

    def assign_representative_species(self):

        self.assign_representative_molecule()
        self.species_repr = Species(molecule=[self.mol_repr])

    def getMolecularWeight(self):
        """
        Return the fragmental weight of the fragment in kg/mol.
        """
        mass = 0
        for vertex in self.vertices:
            if isinstance(vertex, Atom):
                mass += vertex.element.mass
        return mass

    def isRadical(self):
        """
        Return ``True`` if the fragment contains at least one radical electron,
        or ``False`` otherwise.
        """
        for vertex in self.vertices:
            if isinstance(vertex, Atom) and vertex.radicalElectrons > 0:
                return True
        return False

    def update(self):

        for v in self.vertices:
            if isinstance(v, Atom):
                v.updateCharge()

        self.updateAtomTypes()
        self.updateMultiplicity()
        self.sortVertices()

    def updateAtomTypes(self, logSpecies=True, raiseException=True):
        """
        Iterate through the atoms in the structure, checking their atom types
        to ensure they are correct (i.e. accurately describe their local bond
        environment) and complete (i.e. are as detailed as possible).
        
        If `raiseException` is `False`, then the generic atomType 'R' will
        be prescribed to any atom when getAtomType fails. Currently used for
        resonance hybrid atom types.
        """
        #Because we use lonepairs to match atomtypes and default is -100 when unspecified,
        #we should update before getting the atomtype.
        self.updateLonePairs()

        for v in self.vertices:
            if not isinstance(v, Atom): continue
            try:
                v.atomType = getAtomType(v, v.edges)
            except AtomTypeError:
                if logSpecies:
                    logging.error("Could not update atomtypes for this fragment:\n{0}".format(self.toAdjacencyList()))
                if raiseException:
                    raise
                v.atomType = atomTypes['R']

    def updateMultiplicity(self):
        """
        Update the multiplicity of a newly formed molecule.
        """
        # Assume this is always true
        # There are cases where 2 radicalElectrons is a singlet, but
        # the triplet is often more stable, 
        self.multiplicity = self.getRadicalCount() + 1


    def updateLonePairs(self):
        """
        Iterate through the atoms in the structure and calculate the
        number of lone electron pairs, assuming a neutral molecule.
        """
        for v in self.vertices:
            if not isinstance(v, Atom): continue
            if not v.isHydrogen():
                order = v.getBondOrdersForAtom()
                v.lonePairs = (elements.PeriodicSystem.valence_electrons[v.symbol] - v.radicalElectrons - v.charge - int(order)) / 2.0
                if v.lonePairs % 1 > 0 or v.lonePairs > 4:
                    logging.error("Unable to determine the number of lone pairs for element {0} in {1}".format(v,self))
            else:
                v.lonePairs = 0

    def getFormula(self):
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
        keys = elements.keys()
        keys.sort()
        for key in keys:
            count = elements[key]
            formula += '{0}{1:d}'.format(key, count) if count > 1 else key
        
        return formula

    def get_representative_molecule(self, mode='minimal', update=True):

        if mode == 'minimal':
            # create a molecule from fragment.vertices.copy
            mapping = self.copyAndMap()

            # replace CuttingLabel with H
            atoms = []
            for vertex in self.vertices:

                mapped_vertex = mapping[vertex]
                if isinstance(mapped_vertex, CuttingLabel):

                    # replace cutting label with atom H
                    atom_H = Atom(element=getElement('H'), 
                                radicalElectrons=0, 
                                charge=0, 
                                lonePairs=0)

                    for bondedAtom, bond in mapped_vertex.edges.iteritems():
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


    def toRDKitMol(self, removeHs=False, returnMapping=True):
        """
        Convert a molecular structure to a RDKit rdmol object.
        """
        if removeHs: 
            # because we're replacing
            # cutting labels with hydrogens
            # so do not allow removeHs to be True
            raise "Currently fragment toRDKitMol only allows keeping all the hydrogens."

        mol0, mapping = self.get_representative_molecule('minimal', update=False)

        rdmol, rdAtomIdx_mol0 = converter.toRDKitMol(mol0, removeHs=removeHs, 
                                                     returnMapping=returnMapping, 
                                                     sanitize=True)

        rdAtomIdx_frag = {}
        for frag_atom, mol0_atom in mapping.iteritems():
            rd_idx = rdAtomIdx_mol0[mol0_atom]
            rdAtomIdx_frag[frag_atom] = rd_idx

        # sync the order of fragment vertices with the order
        # of mol0.atoms since mol0.atoms is changed/sorted in 
        # converter.toRDKitMol().
        # Since the rdmol's atoms order is same as the order of mol0's atoms,
        # the synchronization between fragment.atoms order and mol0.atoms order
        # is necessary to make sure the order of fragment vertices
        # reflects the order of rdmol's atoms
        vertices_order = []
        for v in self.vertices:
            a = mapping[v]
            idx = mol0.atoms.index(a)
            vertices_order.append((v, idx))

        adapted_vertices = [tup[0] for tup in sorted(vertices_order, key=lambda tup: tup[1])]

        self.vertices = adapted_vertices

        return rdmol, rdAtomIdx_frag

    def toAdjacencyList(self, 
                        label='', 
                        removeH=False, 
                        removeLonePairs=False, 
                        oldStyle=False):
        """
        Convert the molecular structure to a string adjacency list.
        """
        from rmgpy.molecule.adjlist import toAdjacencyList
        result = toAdjacencyList(self.vertices, 
                                 self.multiplicity,  
                                 label=label, 
                                 group=False, 
                                 removeH=removeH, 
                                 removeLonePairs=removeLonePairs, 
                                 oldStyle=oldStyle)
        return result

    def fromAdjacencyList(self, adjlist, saturateH=False):
        """
        Convert a string adjacency list `adjlist` to a fragment structure.
        Skips the first line (assuming it's a label) unless `withLabel` is
        ``False``.
        """
        from rmgpy.molecule.adjlist import fromAdjacencyList
        
        self.vertices, self.multiplicity = fromAdjacencyList(adjlist, 
                                                             group=False, 
                                                             saturateH=saturateH)
        self.updateAtomTypes()
        
        # Check if multiplicity is possible
        n_rad = self.getRadicalCount() 
        multiplicity = self.multiplicity
        if not (n_rad + 1 == multiplicity or n_rad - 1 == multiplicity or n_rad - 3 == multiplicity or n_rad - 5 == multiplicity):
            raise ValueError('Impossible multiplicity for molecule\n{0}\n multiplicity = {1} and number of unpaired electrons = {2}'.format(self.toAdjacencyList(),multiplicity,n_rad))
        if self.getNetCharge() != 0:
            raise ValueError('Non-neutral molecule encountered. Currently, AFM does not support ion chemistry.\n {0}'.format(adjlist))
        return self

    def getAromaticRings(self, rings=None):
        """
        Returns all aromatic rings as a list of atoms and a list of bonds.

        Identifies rings using `Graph.getSmallestSetOfSmallestRings()`, then uses RDKit to perceive aromaticity.
        RDKit uses an atom-based pi-electron counting algorithm to check aromaticity based on Huckel's Rule.
        Therefore, this method identifies "true" aromaticity, rather than simply the RMG bond type.

        The method currently restricts aromaticity to six-membered carbon-only rings. This is a limitation imposed
        by RMG, and not by RDKit.
        """

        from rdkit.Chem.rdchem import BondType
        AROMATIC = BondType.AROMATIC

        if rings is None:
            rings = self.getRelevantCycles()
            rings = [ring for ring in rings if len(ring) == 6]
        if not rings:
            return [], []

        try:
            rdkitmol, rdAtomIndices = self.toRDKitMol(removeHs=False, returnMapping=True)
        except ValueError:
            logging.warning('Unable to check aromaticity by converting to RDKit Mol.')
        else:
            aromaticRings = []
            aromaticBonds = []
            for ring0 in rings:
                aromaticBondsInRing = []
                # Figure out which atoms and bonds are aromatic and reassign appropriately:
                for i, atom1 in enumerate(ring0):
                    if not atom1.isCarbon():
                        # all atoms in the ring must be carbon in RMG for our definition of aromatic
                        break
                    for atom2 in ring0[i + 1:]:
                        if self.hasBond(atom1, atom2):
                            if rdkitmol.GetBondBetweenAtoms(rdAtomIndices[atom1],
                                                            rdAtomIndices[atom2]).GetBondType() is AROMATIC:
                                aromaticBondsInRing.append(self.getBond(atom1, atom2))
                else:  # didn't break so all atoms are carbon
                    if len(aromaticBondsInRing) == 6:
                        aromaticRings.append(ring0)
                        aromaticBonds.append(aromaticBondsInRing)

            return aromaticRings, aromaticBonds

    def isAromatic(self):
        """ 
        Returns ``True`` if the fragment is aromatic, or ``False`` if not.  
        """
        mol0, _ = self.get_representative_molecule('minimal')   
        return mol0.isAromatic()

    def atomIDValid(self):
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

    def assignAtomIDs(self):
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

    def generate_resonance_structures(self, keep_isomorphic=False, filter_structures=True):
        """Returns a list of resonance structures of the fragment."""
        return resonance.generate_resonance_structures(self, keep_isomorphic=keep_isomorphic,
                                                        filter_structures = filter_structures)

    def isIdentical(self, other, strict=True):
        """
        Performs isomorphism checking, with the added constraint that atom IDs must match.

        Primary use case is tracking atoms in reactions for reaction degeneracy determination.

        Returns :data:`True` if two graphs are identical and :data:`False` otherwise.
        """

        if not isinstance(other, (Fragment, Molecule)):
            raise TypeError('Got a {0} object for parameter "other", when a Fragment object is required.'.format(other.__class__))

        # Get a set of atom indices for each molecule
        atomIDs = set([atom.id for atom in self.atoms])
        otherIDs = set([atom.id for atom in other.atoms])

        if atomIDs == otherIDs:
            # If the two molecules have the same indices, then they might be identical
            # Sort the atoms by ID
            atomList = sorted(self.atoms, key=lambda x: x.id)
            otherList = sorted(other.atoms, key=lambda x: x.id)

            # If matching atom indices gives a valid mapping, then the molecules are fully identical
            mapping = {}
            for atom1, atom2 in itertools.izip(atomList, otherList):
                mapping[atom1] = atom2

            return self.isMappingValid(other, mapping, equivalent=True, strict=strict)
        else:
            # The molecules don't have the same set of indices, so they are not identical
            return False

    def toSMILES(self):

        cutting_label_list = []
        for vertex in self.vertices:
            if isinstance(vertex, CuttingLabel):
                cutting_label_list.append(vertex.symbol)

        SMILES_before = self.copy(deep=True)
        final_vertices = []
        for ind, atom in enumerate(SMILES_before.atoms):
            element_symbol = atom.symbol
            if isinstance(atom, CuttingLabel):
                substi_name = 'Si'
                substi = Atom(element=substi_name)
                substi.label = element_symbol

                for bondedAtom, bond in atom.edges.iteritems():
                    new_bond = Bond(bondedAtom, substi, order=bond.order)

                    bondedAtom.edges[substi] = new_bond
                    del bondedAtom.edges[atom]

                    substi.edges[bondedAtom] = new_bond

                substi.radicalElectrons = 3

                final_vertices.append(substi)
            else:
                final_vertices.append(atom)

        SMILES_before.vertices = final_vertices
        mol_repr = Molecule()
        mol_repr.atoms = SMILES_before.vertices
        SMILES_after = mol_repr.toSMILES()
        import re
        smiles = re.sub('\[Si\]', '', SMILES_after)

        return smiles

    def get_element_count(self):
        """
        Returns the element count for the fragment as a dictionary.
        """
        element_count = {}
        for atom in self.vertices:
            if not isinstance(atom, Atom): continue
            symbol = atom.element.symbol
            isotope = atom.element.isotope
            key = symbol if isotope == -1 else (symbol, isotope)
            if key in element_count:
                element_count[key] += 1
            else:
                element_count[key] = 1

        return element_count

    def getURL(self):
        """
        Get a URL to the fragment's info page on the RMG website.
        """

        base_url = "http://rmg.mit.edu/database/molecule/"
        adjlist = self.toAdjacencyList(removeH=False)
        url = base_url + urllib.quote(adjlist)
        return url.strip('_')

    def isLinear(self):
        """
        Return :data:`True` if the structure is linear and :data:`False`
        otherwise.
        """

        atomCount = len(self.vertices)

        # Monatomic molecules are definitely nonlinear
        if atomCount == 1:
            return False
        # Diatomic molecules are definitely linear
        elif atomCount == 2:
            return True
        # Cyclic molecules are definitely nonlinear
        elif self.isCyclic():
            return False

        # True if all bonds are double bonds (e.g. O=C=O)
        allDoubleBonds = True
        for atom1 in self.vertices:
            for bond in atom1.edges.values():
                if not bond.isDouble(): allDoubleBonds = False
        if allDoubleBonds: return True

        # True if alternating single-triple bonds (e.g. H-C#C-H)
        # This test requires explicit hydrogen atoms
        for atom in self.vertices:
            bonds = atom.edges.values()
            if len(bonds)==1:
                continue # ok, next atom
            if len(bonds)>2:
                break # fail!
            if bonds[0].isSingle() and bonds[1].isTriple():
                continue # ok, next atom
            if bonds[1].isSingle() and bonds[0].isTriple():
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
            for i in range(atom.radicalElectrons):
                H = Atom('H', radicalElectrons=0, lonePairs=0, charge=0)
                bond = Bond(atom, H, 1)
                self.addAtom(H)
                self.addBond(bond)
                if atom not in added:
                    added[atom] = []
                added[atom].append([H, bond])
                atom.decrementRadical()
      
        # Update the atom types of the saturated structure (not sure why
        # this is necessary, because saturating with H shouldn't be
        # changing atom types, but it doesn't hurt anything and is not
        # very expensive, so will do it anyway)
        self.sortVertices()
        self.updateAtomTypes()
        self.multiplicity = 1

        return added

    def isArylRadical(self, aromaticRings=None):
        """
        Return ``True`` if the fragment only contains aryl radicals,
        ie. radical on an aromatic ring, or ``False`` otherwise.
        """
        if aromaticRings is None:
            aromaticRings = self.getAromaticRings()[0]

        total = self.getRadicalCount()
        aromaticAtoms = set([atom for atom in itertools.chain.from_iterable(aromaticRings)])
        aryl = sum([atom.radicalElectrons for atom in aromaticAtoms])

        return total == aryl

    def getNumAtoms(self, element = None):
        """
        Return the number of atoms in molecule.  If element is given, ie. "H" or "C",
        the number of atoms of that element is returned.
        """
        numAtoms = 0
        if element == None:
            for vertex in self.vertices:
                if isinstance(vertex, Atom):
                    numAtoms += 1
        else:
            for vertex in self.vertices:
                if isinstance(vertex, Atom):
                    if vertex.element.symbol == element:
                        numAtoms += 1
        return numAtoms

    def fromRDKitMol(self, rdkitmol, atom_replace_dict = None):
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
        for i in xrange(rdkitmol.GetNumAtoms()):
            rdkitatom = rdkitmol.GetAtomWithIdx(i)

            # Use atomic number as key for element
            number = rdkitatom.GetAtomicNum()
            element = getElement(number)

            # Process charge
            charge = rdkitatom.GetFormalCharge()
            radicalElectrons = rdkitatom.GetNumRadicalElectrons()

            ELE = element.symbol
            if atom_replace_dict.has_key('[' + ELE + ']'):
                cutting_label_name = atom_replace_dict['[' + ELE + ']']
                cutting_label = CuttingLabel(name=cutting_label_name)
                self.vertices.append(cutting_label)
            else:
                atom = Atom(element, radicalElectrons, charge, '', 0)
                self.vertices.append(atom)

            # Add bonds by iterating again through atoms
            for j in xrange(0, i):
                rdkitbond = rdkitmol.GetBondBetweenAtoms(i, j)
                if rdkitbond is not None:
                    order = 0

                    # Process bond type
                    rdbondtype = rdkitbond.GetBondType()
                    if rdbondtype.name == 'SINGLE':
                        order = 1
                    elif rdbondtype.name == 'DOUBLE':
                        order = 2
                    elif rdbondtype.name == 'TRIPLE':
                        order = 3
                    elif rdbondtype.name == 'AROMATIC':
                        order = 1.5

                    bond = Bond(self.vertices[i], self.vertices[j], order)
                    self.addBond(bond)

        # We need to update lone pairs first because the charge was set by RDKit
        self.updateLonePairs()
        # Set atom types and connectivity values
        self.update()

        # Assume this is always true
        # There are cases where 2 radicalElectrons is a singlet, but
        # the triplet is often more stable,
        self.updateMultiplicity()
        # mol.updateAtomTypes()

        return self

# this variable is used to name atom IDs so that there are as few conflicts by 
# using the entire space of integer objects
atom_id_counter = -2**15
