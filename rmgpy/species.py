#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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

"""
This module contains classes and functions for working with chemical species.

From the `IUPAC Compendium of Chemical Terminology 
<https://doi.org/10.1351/goldbook>`_, a chemical species is "an 
ensemble of chemically identical molecular entities that can explore the same 
set of molecular energy levels on the time scale of the experiment". This
definition is purposefully vague to allow the user flexibility in application.

In RMG Py, a chemical species -- a local minimum on a potential energy surface
-- is represented in memory as a :class:`Species` object. This module also
contains the :class:`TransitionState` class for representing chemical reaction
transition states (first-order saddle points on a potential energy surface).
"""

import logging
from copy import deepcopy
from operator import itemgetter

import cython
import numpy as np

import rmgpy.quantity as quantity
from rmgpy.exceptions import SpeciesError, StatmechError
from rmgpy.molecule.graph import Graph
from rmgpy.molecule.molecule import Atom, Bond, Molecule
from rmgpy.molecule.fragment import CuttingLabel, Fragment
from rmgpy.pdep import SingleExponentialDown
from rmgpy.statmech.conformer import Conformer
from rmgpy.thermo import Wilhoit, NASA, ThermoData
from rmgpy.data.vaporLiquidMassTransfer import vapor_liquid_mass_transfer

#: This dictionary is used to add multiplicity to species label
_multiplicity_labels = {1: 'S', 2: 'D', 3: 'T', 4: 'Q', 5: 'V', }


################################################################################

class Species(object):
    """
    A chemical species, representing a local minimum on a potential energy
    surface. The attributes are:

    ======================= ====================================================
    Attribute               Description
    ======================= ====================================================
    `index`                 A unique nonnegative integer index
    `label`                 A descriptive string label
    `thermo`                The heat capacity model for the species
    `conformer`             The molecular conformer for the species
    `molecule`              A list of the :class:`Molecule` objects describing the molecular structure
    `transport_data`        A set of transport collision parameters
    `molecular_weight`      The molecular weight of the species
    `energy_transfer_model` The collisional energy transfer model to use
    `reactive`              ``True`` if the species participates in reaction families, ``False`` if not
                                Reaction libraries and seed mechanisms that include the species are
                                always considered regardless of this variable
    `props`                 A generic 'properties' dictionary to store user-defined flags
    `aug_inchi`             Unique augmented inchi
    `symmetry_number`       Estimated symmetry number of the species, using the resonance hybrid
    `creation_iteration`    Iteration which the species is created within the reaction mechanism generation algorithm
    `explicitly_allowed`    Flag to exempt species from forbidden structure checks
    ======================= ====================================================

    """

    def __init__(self, index=-1, label='', thermo=None, conformer=None, molecule=None, transport_data=None,
                 molecular_weight=None, energy_transfer_model=None, reactive=True, props=None, smiles='', inchi='',
                 aug_inchi=None, symmetry_number=-1, creation_iteration=0, explicitly_allowed=False,
                 liquid_volumetric_mass_transfer_coefficient_data=None,henry_law_constant_data=None):
        self.index = index
        self.label = label
        self.thermo = thermo
        self.conformer = conformer
        self.molecule = molecule or []
        self.transport_data = transport_data
        self.reactive = reactive
        self.molecular_weight = molecular_weight
        self.energy_transfer_model = energy_transfer_model
        self.props = props or {}
        self.aug_inchi = aug_inchi
        self.symmetry_number = symmetry_number
        self.is_solvent = False
        self.creation_iteration = creation_iteration
        self.explicitly_allowed = explicitly_allowed
        self._fingerprint = None
        self._inchi = None
        self._smiles = None
        self.liquid_volumetric_mass_transfer_coefficient_data = liquid_volumetric_mass_transfer_coefficient_data
        self.henry_law_constant_data = henry_law_constant_data

        if inchi and smiles:
            logging.warning('Both InChI and SMILES provided for Species instantiation, '
                            'using InChI and ignoring SMILES.')
        if inchi:
            self.molecule = [Molecule(inchi=inchi)]
            self._inchi = inchi
        elif smiles:
            # check it is fragment or molecule
            _ , cutting_label_list = Fragment().detect_cutting_label(smiles)
            if cutting_label_list != []: # Fragment
                self.molecule = [Fragment(smiles=smiles)]
            else: # Molecule
                self.molecule = [Molecule(smiles=smiles)]
            self._smiles = smiles

        # Check multiplicity of each molecule is the same
        if molecule is not None and len(molecule) > 1:
            mult = molecule[0].multiplicity
            for m in molecule[1:]:
                if mult != m.multiplicity:
                    raise SpeciesError('Multiplicities of molecules in species {species} '
                                       'do not match.'.format(species=label))

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'Species('
        if self.index != -1:
            string += 'index={0:d}, '.format(self.index)
        if self.label != -1:
            string += 'label="{0}", '.format(self.label)
        if self.thermo is not None:
            string += 'thermo={0!r}, '.format(self.thermo)
        if self.conformer is not None:
            string += 'conformer={0!r}, '.format(self.conformer)
        if len(self.molecule) > 0:
            string += 'molecule={0!r}, '.format(self.molecule)
        if self.transport_data is not None:
            string += 'transport_data={0!r}, '.format(self.transport_data)
        if self.liquid_volumetric_mass_transfer_coefficient_data is not None:
            string += f'liquid_volumetric_mass_transfer_coefficient_data={self.liquid_volumetric_mass_transfer_coefficient_data}'
        if self.henry_law_constant_data is not None:
            string += f'henry_law_constant_data={self.henry_law_constant_data}'
        if not self.reactive:
            string += 'reactive={0}, '.format(self.reactive)
        if self.molecular_weight is not None:
            string += 'molecular_weight={0!r}, '.format(self.molecular_weight)
        if self.energy_transfer_model is not None:
            string += 'energy_transfer_model={0!r}, '.format(self.energy_transfer_model)
        string = string[:-2] + ')'
        return string

    def _repr_png_(self):
        if len(self.molecule) > 0:
            return self.molecule[0]._repr_png_()
        else:
            return None

    def __str__(self):
        """
        Return a string representation of the species, in the form 'label(id)'.
        """
        if not self.label:
            self.label = self.molecule[0].to_smiles()
        if self.index == -1:
            return self.label
        else:
            return '{0}({1:d})'.format(self.label, self.index)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Species, (self.index, self.label, self.thermo, self.conformer, self.molecule, self.transport_data,
                          self.molecular_weight, self.energy_transfer_model, self.reactive, self.props))

    def __hash__(self):
        """
        Define a custom hash method to allow Species objects to be used in dictionaries and sets.

        Use the fingerprint property, which is taken from the first molecule entry.
        This is currently defined as the molecular formula, which is not an ideal hash, since there will be significant
        hash collisions, leading to inefficient lookups.
        """
        return hash(('Species', self.fingerprint))

    def __eq__(self, other):
        """Define equality comparison. Define as a reference comparison"""
        return self is other

    def __lt__(self, other):
        """Define less than comparison. For comparing against other Species objects (e.g. when sorting)."""
        if isinstance(other, Species):
            return self.sorting_key < other.sorting_key
        else:
            raise NotImplementedError('Cannot perform less than comparison between Species and '
                                      '{0}.'.format(type(other).__name__))

    def __gt__(self, other):
        """Define greater than comparison. For comparing against other Species objects (e.g. when sorting)."""
        if isinstance(other, Species):
            return self.sorting_key > other.sorting_key
        else:
            raise NotImplementedError('Cannot perform greater than comparison between Species and '
                                      '{0}.'.format(type(other).__name__))

    @property
    def sorting_key(self):
        """Returns a sorting key for comparing Species objects. Read-only"""
        return self.fingerprint, self.label, self.index

    @property
    def fingerprint(self):
        """Fingerprint of this species, taken from molecule attribute. Read-only."""
        if self._fingerprint is None:
            if self.molecule:
                self._fingerprint = self.molecule[0].fingerprint
        return self._fingerprint

    @property
    def inchi(self):
        """InChI string representation of this species. Read-only."""
        if self._inchi is None:
            if self.molecule:
                self._inchi = self.molecule[0].inchi
        return self._inchi

    @property
    def smiles(self):
        """
        SMILES string representation of this species. Read-only.

        Note that SMILES representations for different resonance structures of the same species may be different.
        """
        if self._smiles is None:
            if self.molecule:
                self._smiles = self.molecule[0].smiles
        return self._smiles

    @property
    def multiplicity(self):
        """Fingerprint of this species, taken from molecule attribute. Read-only."""
        if self.molecule:
            return self.molecule[0].multiplicity
        else:
            return None

    @property
    def molecular_weight(self):
        """The molecular weight of the species. (Note: value_si is in kg/molecule not kg/mol)"""
        if self._molecular_weight is None and self.molecule is not None and len(self.molecule) > 0:
            self._molecular_weight = quantity.Mass(self.molecule[0].get_molecular_weight(), 'kg/mol')
        return self._molecular_weight

    @molecular_weight.setter
    def molecular_weight(self, value):
        self._molecular_weight = quantity.Mass(value)

    def generate_resonance_structures(self, keep_isomorphic=True, filter_structures=True, save_order=False):
        """
        Generate all of the resonance structures of this species. The isomers are
        stored as a list in the `molecule` attribute. If the length of
        `molecule` is already greater than one, it is assumed that all of the
        resonance structures have already been generated.
        If ``save_order`` is ``True`` the atom order is reset after performing atom isomorphism.
        """
        if len(self.molecule) == 1 or not self.molecule[0].atom_ids_valid():
            if not self.molecule[0].atom_ids_valid():
                self.molecule[0].assign_atom_ids()
            self.molecule = self.molecule[0].generate_resonance_structures(keep_isomorphic=keep_isomorphic,
                                                                           filter_structures=filter_structures,
                                                                           save_order=save_order
                                                                           )

    def is_isomorphic(self, other, generate_initial_map=False, save_order=False, strict=True):
        """
        Return ``True`` if the species is isomorphic to `other`, which can be
        either a :class:`Molecule` object or a :class:`Species` object.

        Args:
            generate_initial_map (bool, optional): If ``True``, make initial map by matching labeled atoms
            save_order (bool, optional):           if ``True``, reset atom order after performing atom isomorphism
            strict (bool, optional):               If ``False``, perform isomorphism ignoring electrons.
        """
        if isinstance(other, Molecule) or isinstance(other, Fragment):
            for molecule in self.molecule:
                if molecule.is_isomorphic(other, generate_initial_map=generate_initial_map,
                                          save_order=save_order, strict=strict):
                    return True
                elif not strict:
                    return False
        elif isinstance(other, Species):
            for molecule1 in self.molecule:
                for molecule2 in other.molecule:
                    if molecule1.is_isomorphic(molecule2, generate_initial_map=generate_initial_map,
                                               save_order=save_order, strict=strict):
                        return True
                    elif not strict:
                        return False
        else:
            raise ValueError('Unexpected value "{0!r}" for other parameter;'
                             ' should be a Molecule or Species object.'.format(other))
        return False

    def is_identical(self, other, strict=True):
        """
        Return ``True`` if at least one molecule of the species is identical to `other`,
        which can be either a :class:`Molecule` object or a :class:`Species` object.

        If ``strict=False``, performs the check ignoring electrons and resonance structures.
        """
        if isinstance(other, Molecule) or isinstance(other, Fragment):
            for molecule in self.molecule:
                if molecule.is_identical(other, strict=strict):
                    return True
        elif isinstance(other, Species):
            for molecule1 in self.molecule:
                for molecule2 in other.molecule:
                    if molecule1.is_identical(molecule2, strict=strict):
                        return True
        else:
            raise ValueError('Unexpected value "{0!r}" for other parameter;'
                             ' should be a Molecule or Species object.'.format(other))
        return False

    def is_structure_in_list(self, species_list):
        """
        Return ``True`` if at least one Molecule in self is isomorphic with at least one other Molecule in at least
        one Species in species list.
        """
        for species in species_list:
            if isinstance(species, Species):
                return self.is_isomorphic(species)
            else:
                raise TypeError('Unexpected value "{0!r}" for species_list parameter;'
                                ' should be a List of Species objects.'.format(species))
        return False

    def from_adjacency_list(self, adjlist, raise_atomtype_exception=True, raise_charge_exception=False):
        """
        Load the structure of a species as a :class:`Molecule` object from the
        given adjacency list `adjlist` and store it as the first entry of a 
        list in the `molecule` attribute. Does not generate resonance isomers
        of the loaded molecule.
        """
        # detect if it contains cutting label
        _ , cutting_label_list = Fragment().detect_cutting_label(adjlist)
        if cutting_label_list == []:
            self.molecule = [Molecule().from_adjacency_list(adjlist, saturate_h=False,
                                                            raise_atomtype_exception=raise_atomtype_exception,
                                                            raise_charge_exception=raise_charge_exception)]
        else:
            self.molecule = [Fragment().from_adjacency_list(adjlist, saturate_h=False,
                                                            raise_atomtype_exception=raise_atomtype_exception,
                                                            raise_charge_exception=raise_charge_exception)]
        # If the first line is a label, then save it to the label attribute
        for label in adjlist.splitlines():
            if label.strip():
                break
        else:
            label = ''
        if len(label.split()) > 0 and not label.split()[0].isdigit() and 'multiplicity' not in label:
            self.label = label.strip()
        # Return a reference to itself so we can use e.g. Species().from_adjacency_list()
        return self

    def from_smiles(self, smiles):
        """
        Load the structure of a species as a :class:`Molecule` object from the
        given SMILES string `smiles` and store it as the first entry of a 
        list in the `molecule` attribute. Does not generate resonance isomers
        of the loaded molecule.
        """
        self.molecule = [Molecule().from_smiles(smiles)]
        # Return a reference to itself so we can use e.g. Species().from_adjacency_list()
        return self

    def to_adjacency_list(self):
        """
        Return a string containing each of the molecules' adjacency lists.
        """
        output = '\n\n'.join([m.to_adjacency_list(label=self.label, remove_h=False) for m in self.molecule])
        return output

    def to_chemkin(self):
        """
        Return the chemkin-formatted string for this species.
        """
        from rmgpy.chemkin import get_species_identifier
        return get_species_identifier(self)

    def to_cantera(self, use_chemkin_identifier=False):
        """
        Converts the RMG Species object to a Cantera Species object
        with the appropriate thermo data.

        If use_chemkin_identifier is set to False, the species label is used
        instead. Be sure that species' labels are unique when setting it False.
        """
        import cantera as ct

        # Determine the number of each type of element in the molecule
        element_dict = {}  # element_counts = [0,0,0,0]
        for vertex in self.molecule[0].vertices:
            # The atom itself
            if not isinstance(vertex, CuttingLabel):
                symbol = vertex.element.symbol
            else: # that means this vertex is CuttingLabel
                continue
            if symbol not in element_dict:
                element_dict[symbol] = 1
            else:
                element_dict[symbol] += 1
        if use_chemkin_identifier:
            ct_species = ct.Species(self.to_chemkin(), element_dict)
        else:
            ct_species = ct.Species(self.label, element_dict)
        if self.thermo:
            try:
                ct_species.thermo = self.thermo.to_cantera()
            except Exception:
                logging.error('Could not convert thermo to create Cantera Species object. '
                              'Check that thermo is a NASA polynomial.')
                raise

        if self.transport_data:
            ct_species.transport = self.transport_data.to_cantera()

        return ct_species

    def has_statmech(self):
        """
        Return ``True`` if the species has statistical mechanical parameters,
        or ``False`` otherwise.
        """
        if (len(self.molecule) > 0 and len(self.molecule[0].atoms) == 1):
            # atomic molecules have no modes, check only E0
            return self.conformer is not None and self.conformer.E0 is not None
        else:
            # polyatomic molecules should have modes and E0, so check both
            return self.conformer is not None and len(self.conformer.modes) > 0 and self.conformer.E0 is not None

    def has_thermo(self):
        """
        Return ``True`` if the species has thermodynamic parameters, or 
        ``False`` otherwise.
        """
        return self.thermo is not None

    def contains_surface_site(self):
        """
        Return ``True`` if the species is adsorbed on a surface (or is itself a site), else ``False``.
        """
        return self.molecule[0].contains_surface_site()

    def is_surface_site(self):
        """Return ``True`` if the species is a vacant surface site."""
        return self.molecule[0].is_surface_site()

    def is_electron(self):
        """Return ``True`` if the species is an electron"""
        
        if len(self.molecule) == 0:
            return False
        else:
            return self.molecule[0].is_electron()

    def is_proton(self):
        """Return ``True`` if the species is a proton"""
        
        if len(self.molecule) == 0:
            return False
        else:
            return self.molecule[0].is_proton()

    def get_partition_function(self, T):
        """
        Return the partition function for the species at the specified
        temperature `T` in K.
        """
        cython.declare(Q=cython.double)
        if self.has_statmech():
            Q = self.conformer.get_partition_function(T)
        else:
            raise Exception('Unable to calculate partition function for species {0!r}: '
                            'no statmech data available.'.format(self.label))
        return Q

    def get_heat_capacity(self, T):
        """
        Return the heat capacity in J/mol*K for the species at the specified
        temperature `T` in K.
        """
        cython.declare(Cp=cython.double)
        Cp = 0.0
        if self.has_thermo():
            Cp = self.get_thermo_data().get_heat_capacity(T)
        elif self.has_statmech():
            Cp = self.conformer.get_heat_capacity(T)
        else:
            raise Exception('Unable to calculate heat capacity for species {0!r}: '
                            'no thermo or statmech data available.'.format(self.label))
        return Cp

    def get_enthalpy(self, T):
        """
        Return the enthalpy in J/mol for the species at the specified
        temperature `T` in K.
        """
        cython.declare(H=cython.double)
        H = 0.0
        if self.has_thermo():
            H = self.get_thermo_data().get_enthalpy(T)
        elif self.has_statmech():
            H = self.conformer.get_enthalpy(T) + self.conformer.E0.value_si
        else:
            raise Exception('Unable to calculate enthalpy for species {0!r}: '
                            'no thermo or statmech data available.'.format(self.label))
        return H

    def get_entropy(self, T):
        """
        Return the entropy in J/mol*K for the species at the specified
        temperature `T` in K.
        """
        cython.declare(S=cython.double)
        S = 0.0
        if self.has_thermo():
            S = self.get_thermo_data().get_entropy(T)
        elif self.has_statmech():
            S = self.conformer.get_entropy(T)
        else:
            raise Exception('Unable to calculate entropy for species {0!r}: '
                            'no thermo or statmech data available.'.format(self.label))
        return S

    def get_free_energy(self, T):
        """
        Return the Gibbs free energy in J/mol for the species at the specified
        temperature `T` in K.
        """
        cython.declare(G=cython.double)
        G = 0.0
        if self.has_thermo():
            G = self.get_thermo_data().get_free_energy(T)
        elif self.has_statmech():
            G = self.conformer.get_free_energy(T) + self.conformer.E0.value_si
        else:
            raise Exception('Unable to calculate free energy for species {0!r}: '
                            'no thermo or statmech data available.'.format(self.label))
        return G

    def get_sum_of_states(self, e_list):
        """
        Return the sum of states :math:`N(E)` at the specified energies `e_list`
        in J/mol.
        """
        if self.has_statmech():
            return self.conformer.get_sum_of_states(e_list)
        else:
            raise Exception('Unable to calculate sum of states for species {0!r}: '
                            'no statmech data available.'.format(self.label))

    def get_density_of_states(self, e_list):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `e_list` in J/mol above the ground state.
        """
        if self.has_statmech():
            try:
                return self.conformer.get_density_of_states(e_list)
            except StatmechError:
                logging.error('StatmechError raised for species {0}'.format(self.label))
                raise
        else:
            raise Exception('Unable to calculate density of states for species {0!r}: '
                            'no statmech data available.'.format(self.label))

    def get_symmetry_number(self):
        """
        Get the symmetry number for the species, which is the highest symmetry number amongst
        its resonance isomers and the resonance hybrid.  
        This function is currently used for website purposes and testing only as it
        requires additional calculate_symmetry_number calls.
        """
        if self.symmetry_number < 1:
            if isinstance(self.molecule[0], Molecule):
                cython.declare(resonanceHybrid=Molecule, maxSymmetryNum=cython.short)
                resonance_hybrid = self.get_resonance_hybrid()
                try:
                    self.symmetry_number = resonance_hybrid.get_symmetry_number()
                except KeyError:
                    logging.error('Wrong bond order generated by resonance hybrid.')
                    logging.error('Resonance Hybrid: {}'.format(resonance_hybrid.to_adjacency_list()))
                    for index, mol in enumerate(self.molecule):
                        logging.error("Resonance Structure {}: {}".format(index, mol.to_adjacency_list()))
                    raise
            else:
                self.symmetry_number = self.molecule[0].get_symmetry_number()
        return self.symmetry_number

    def get_resonance_hybrid(self):
        """
        Returns a molecule object with bond orders that are the average 
        of all the resonance structures.
        """
        # get labeled resonance isomers
        self.generate_resonance_structures(keep_isomorphic=True)

        # only consider reactive molecules as representative structures
        molecules = [mol for mol in self.molecule if mol.reactive]

        # return if no resonance
        if len(molecules) == 1:
            return molecules[0]

        # create a sorted list of atom objects for each resonance structure
        cython.declare(atomsFromStructures=list, oldAtoms=list, newAtoms=list,
                       numResonanceStructures=cython.short, structureNum=cython.short,
                       oldBondOrder=cython.float,
                       index1=cython.short, index2=cython.short,
                       newMol=Molecule, oldMol=Molecule,
                       atom1=Atom, atom2=Atom,
                       bond=Bond,
                       atoms=list, )

        atoms_from_structures = []
        for new_mol in molecules:
            new_mol.atoms.sort(key=lambda atom: atom.id)
            atoms_from_structures.append(new_mol.atoms)

        num_resonance_structures = len(molecules)

        # make original structure with no bonds
        new_mol = Molecule()
        original_atoms = atoms_from_structures[0]
        for atom1 in original_atoms:
            atom = new_mol.add_atom(Atom(atom1.element))
            atom.id = atom1.id

        new_atoms = new_mol.atoms

        # initialize bonds to zero order
        for index1, atom1 in enumerate(original_atoms):
            for atom2 in atom1.bonds:
                index2 = original_atoms.index(atom2)
                bond = Bond(new_atoms[index1], new_atoms[index2], 0)
                new_mol.add_bond(bond)

        # set bonds to the proper value
        for structureNum, oldMol in enumerate(molecules):
            old_atoms = atoms_from_structures[structureNum]

            for index1, atom1 in enumerate(old_atoms):
                # make bond orders average of resonance structures
                for atom2 in atom1.bonds:
                    index2 = old_atoms.index(atom2)

                    new_bond = new_mol.get_bond(new_atoms[index1], new_atoms[index2])
                    old_bond_order = oldMol.get_bond(old_atoms[index1], old_atoms[index2]).get_order_num()
                    new_bond.apply_action(('CHANGE_BOND', None, old_bond_order / num_resonance_structures / 2))
                # set radicals in resonance hybrid to maximum of all structures
                if atom1.radical_electrons > 0:
                    new_atoms[index1].radical_electrons = max(atom1.radical_electrons,
                                                              new_atoms[index1].radical_electrons)
        new_mol.update_atomtypes(log_species=False, raise_exception=False)
        return new_mol

    def calculate_cp0(self):
        """
        Return the value of the heat capacity at zero temperature in J/mol*K.
        """
        return self.molecule[0].calculate_cp0()

    def calculate_cpinf(self):
        """
        Return the value of the heat capacity at infinite temperature in J/mol*K.
        """
        return self.molecule[0].calculate_cpinf()

    def has_reactive_molecule(self):
        """
        `True` if the species has at least one reactive molecule, `False` otherwise
        """
        cython.declare(molecule=Graph)
        return any([molecule.reactive for molecule in self.molecule])

    def copy(self, deep=False):
        """
        Create a copy of the current species. If the 
        kw argument 'deep' is True, then a deep copy will be made of the 
        Molecule objects in self.molecule.

        For other complex attributes, a deep copy will always be made.
        """
        cython.declare(other=Species)

        other = Species.__new__(Species)

        other.index = self.index

        other.label = self.label

        other.thermo = deepcopy(self.thermo)

        other.molecule = []
        for mol in self.molecule:
            other.molecule.append(mol.copy(deep=deep))

        other.conformer = deepcopy(self.conformer)

        other.transport_data = deepcopy(self.transport_data)

        other.molecular_weight = deepcopy(self.molecular_weight)
        other.energy_transfer_model = deepcopy(self.energy_transfer_model)
        other.reactive = self.reactive
        other.props = deepcopy(self.props)

        return other

    def get_augmented_inchi(self):
        if self.aug_inchi is None:
            self.aug_inchi = self.generate_aug_inchi()
        return self.aug_inchi

    def generate_aug_inchi(self):
        candidates = []
        self.generate_resonance_structures()
        for mol in self.molecule:
            try:
                cand = [mol.to_augmented_inchi(), mol]
            except ValueError:
                pass  # not all resonance structures can be parsed into InChI (e.g. if containing a hypervalance atom)
            else:
                candidates.append(cand)
        candidates = sorted(candidates, key=itemgetter(0))
        for cand in candidates:
            if all(atom.charge == 0 for atom in cand[1].vertices):
                return cand[0]
        return candidates[0][0]

    def get_thermo_data(self, solvent_name=''):
        """
        Returns a `thermoData` object of the current Species object.

        If the thermo object already exists, it is either of the (Wilhoit, ThermoData)
        type, or it is a Future.

        If the type of the thermo attribute is Wilhoit, or ThermoData,
        then it is converted into a NASA format.

        If it is a Future, then a blocking call is made to retrieve the NASA object.
        If the thermo object did not exist yet, the thermo object is generated.        
        """

        from rmgpy.thermo.thermoengine import submit

        if self.thermo:
            if not isinstance(self.thermo, (NASA, Wilhoit, ThermoData)):
                self.thermo = self.thermo.result()
        else:
            submit(self, solvent_name)
            if not isinstance(self.thermo, (NASA, Wilhoit, ThermoData)):
                self.thermo = self.thermo.result()

        return self.thermo

    def generate_transport_data(self):
        """
        Generate the transport_data parameters for the species.
        """
        from rmgpy.data.rmg import get_db
        try:
            transport_db = get_db('transport')
            if not transport_db: raise Exception
        except Exception:
            logging.debug('Could not obtain the transport database. Not generating transport...')
            raise

        # count = sum([1 for atom in self.molecule[0].vertices if atom.is_non_hydrogen()])
        if isinstance(self.molecule[0], Molecule):
            self.transport_data = transport_db.get_transport_properties(self)[0]
        else:
            # assume it's a species for Fragment
            self.molecule[0].assign_representative_species()
            self.transport_data = transport_db.get_transport_properties(self.molecule[0].species_repr)[0]

    def get_transport_data(self):
        """
        Returns the transport data associated with this species, and
        calculates it if it is not yet available.
        """

        if not self.transport_data:
            self.generate_transport_data()

        return self.transport_data

    def generate_statmech(self):
        """
        Generate molecular degree of freedom data for the species. You must
        have already provided a thermodynamics model using e.g.
        :meth:`generate_thermo_data()`.
        """
        logging.debug("Generating statmech for species {}".format(self.label))
        from rmgpy.data.rmg import get_db
        try:
            statmech_db = get_db('statmech')
            if not statmech_db: raise Exception
        except Exception:
            logging.debug('Could not obtain the stat. mech database. Not generating stat. mech...')
            raise

        molecule = self.molecule[0]
        conformer = statmech_db.get_statmech_data(molecule, self.get_thermo_data())

        if self.conformer is None:
            self.conformer = Conformer()

        if self.conformer.E0 is None:
            self.set_e0_with_thermo()

        self.conformer.modes = conformer.modes
        self.conformer.spin_multiplicity = conformer.spin_multiplicity
        if self.conformer.E0 is None or not self.has_statmech():
            logging.error('The conformer in question is {}'.format(self.conformer))
            raise StatmechError('Species {0} does not have stat mech after generate_statmech called'.format(self.label))

    def set_e0_with_thermo(self):
        """
        Helper method that sets species' E0 using the species' thermo data
        """
        if self.get_thermo_data().E0 is not None:
            self.conformer.E0 = self.get_thermo_data().E0
        else:
            if not self.thermo.Cp0 or not self.thermo.CpInf:
                # set Cp0 and CpInf
                from rmgpy.data.thermo import find_cp0_and_cpinf
                find_cp0_and_cpinf(self, self.thermo)
            self.conformer.E0 = self.get_thermo_data().to_wilhoit().E0

    def generate_energy_transfer_model(self):
        """
        Generate the collisional energy transfer model parameters for the
        species. This "algorithm" is *very* much in need of improvement.
        """
        self.energy_transfer_model = SingleExponentialDown(
            alpha0=(300 * 0.011962, "kJ/mol"),
            T0=(300, "K"),
            n=0.85,
        )

    def set_structure(self, structure):
        """
        Set self.molecule from `structure` which could be either a SMILES string or an adjacency list multi-line string
        """
        if not self.molecule:
            try:
                self.molecule = [Molecule(smiles=structure)]
            except ValueError:
                try:
                    self.molecule = [Molecule().from_adjacency_list(structure)]
                except ValueError:
                    logging.error("Cannot understand the given structure '{0}' of species {1}. Could not "
                                  "interpret it as SMILES nor as adjacency list".format(structure, self.label))
                    raise
            self.generate_resonance_structures()

    def get_henry_law_constant_data(self, Ts=[]):

        if (not Ts) and self.henry_law_constant_data:
            return self.henry_law_constant_data

        if vapor_liquid_mass_transfer.enabled:
            self.henry_law_constant_data = vapor_liquid_mass_transfer.get_henry_law_constant_data(self, Ts=Ts)
            return self.henry_law_constant_data
        else:
            raise Exception('Unable to calculate Henry\'s law coefficients when the vapor liquid mass transfer is not enabled '
                            'or liquid volumetric mass transfer coefficient power law is not provided.')
        

    def get_liquid_volumetric_mass_transfer_coefficient_data(self, Ts=[]):

        if (not Ts) and self.liquid_volumetric_mass_transfer_coefficient_data:
            return self.liquid_volumetric_mass_transfer_coefficient_data

        if vapor_liquid_mass_transfer.enabled:
            self.liquid_volumetric_mass_transfer_coefficient_data = vapor_liquid_mass_transfer.get_liquid_volumetric_mass_transfer_coefficient_data(self, Ts=Ts)
            return self.liquid_volumetric_mass_transfer_coefficient_data
        else:
            raise Exception('Unable to calculate liquid volumetric mass transfer coefficient when the diffusion limiter is not enabled '
                            'or liquid volumetric mass transfer coefficient power law is not provided.')


################################################################################


class TransitionState(object):
    """
    A chemical transition state, representing a first-order saddle point on a
    potential energy surface. The attributes are:

    =============== ============================================================
    Attribute       TDescription
    =============== ============================================================
    `label`         A descriptive string label
    `conformer`     The molecular degrees of freedom model for the species
    `frequency`     The negative frequency of the first-order saddle point
    `tunneling`     The type of tunneling model to use for tunneling through the reaction barrier
    `degeneracy`    The reaction path degeneracy
    =============== ============================================================
    """

    def __init__(self, label='', conformer=None, frequency=None, tunneling=None, degeneracy=1):
        self.label = label
        self.conformer = conformer
        self.frequency = frequency
        self.tunneling = tunneling
        self.degeneracy = degeneracy

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'TransitionState('
        if self.label != '': string += 'label="{0}", '.format(self.label)
        if self.conformer is not None: string += 'conformer={0!r}, '.format(self.conformer)
        if self.frequency is not None: string += 'frequency={0!r}, '.format(self.frequency)
        if self.tunneling is not None: string += 'tunneling={0!r}, '.format(self.tunneling)
        if self.degeneracy != 1: string += 'degeneracy={0}, '.format(self.degeneracy)
        string = string[:-2] + ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (TransitionState, (self.label, self.conformer, self.frequency, self.tunneling, self.degeneracy))

    @property
    def frequency(self):
        """The negative frequency of the first-order saddle point."""
        return self._frequency

    @frequency.setter
    def frequency(self, value):
        self._frequency = quantity.Frequency(value)

    def get_partition_function(self, T):
        """
        Return the partition function for the transition state at the
        specified temperature `T` in K.
        """
        cython.declare(Q=cython.double)
        if self.conformer is not None and len(self.conformer.modes) > 0:
            Q = self.conformer.get_partition_function(T)
        else:
            raise SpeciesError('Unable to calculate partition function for transition state {0!r}: '
                               'no statmech data available.'.format(self.label))
        return Q

    def get_heat_capacity(self, T):
        """
        Return the heat capacity in J/mol*K for the transition state at the
        specified temperature `T` in K.
        """
        cython.declare(Cp=cython.double)
        Cp = 0.0

        if self.get_thermo_data() is not None:
            Cp = self.get_thermo_data().get_heat_capacity(T)
        elif self.conformer is not None and len(self.conformer.modes) > 0:
            Cp = self.conformer.get_heat_capacity(T)
        else:
            raise Exception('Unable to calculate heat capacity for transition state {0!r}: '
                            'no thermo or statmech data available.'.format(self.label))
        return Cp

    def get_enthalpy(self, T):
        """
        Return the enthalpy in J/mol for the transition state at the
        specified temperature `T` in K.
        """
        cython.declare(H=cython.double)
        H = 0.0

        if self.get_thermo_data() is not None:
            H = self.get_thermo_data().get_enthalpy(T)
        elif self.conformer is not None and len(self.conformer.modes) > 0:
            H = self.conformer.get_enthalpy(T)
        else:
            raise Exception('Unable to calculate enthalpy for transition state {0!r}: '
                            'no thermo or statmech data available.'.format(self.label))
        return H

    def get_entropy(self, T):
        """
        Return the entropy in J/mol*K for the transition state at the
        specified temperature `T` in K.
        """
        cython.declare(S=cython.double)
        S = 0.0

        if self.get_thermo_data() is not None:
            S = self.get_thermo_data().get_entropy(T)
        elif self.conformer is not None and len(self.conformer.modes) > 0:
            S = self.conformer.get_entropy(T)
        else:
            raise Exception('Unable to calculate entropy for transition state {0!r}: '
                            'no thermo or statmech data available.'.format(self.label))
        return S

    def get_free_energy(self, T):
        """
        Return the Gibbs free energy in J/mol for the transition state at the
        specified temperature `T` in K.
        """
        cython.declare(G=cython.double)
        G = 0.0

        if self.get_thermo_data() is not None:
            G = self.get_thermo_data().get_free_energy(T)
        elif self.conformer is not None and len(self.conformer.modes) > 0:
            G = self.conformer.get_free_energy(T)
        else:
            raise Exception('Unable to calculate free energy for transition state {0!r}: '
                            'no thermo or statmech data available.'.format(self.label))
        return G

    def get_sum_of_states(self, e_list):
        """
        Return the sum of states :math:`N(E)` at the specified energies `e_list`
        in J/mol.
        """
        if self.conformer is not None and len(self.conformer.modes) > 0:
            return self.conformer.get_sum_of_states(e_list)
        else:
            raise Exception('Unable to calculate sum of states for transition state {0!r}: '
                            'no statmech data available.'.format(self.label))

    def get_density_of_states(self, e_list):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `e_list` in J/mol above the ground state.
        """
        if self.conformer is not None and len(self.conformer.modes) > 0:
            return self.conformer.get_density_of_states(e_list)
        else:
            raise Exception('Unable to calculate density of states for transition state {0!r}: '
                            'no statmech data available.'.format(self.label))

    def calculate_tunneling_factor(self, T):
        """
        Calculate and return the value of the canonical tunneling correction 
        factor for the reaction at the given temperature `T` in K.
        """
        if self.tunneling is not None:
            return self.tunneling.calculate_tunneling_factor(T)
        else:
            # Return unity
            return 1.0

    def calculate_tunneling_function(self, e_list):
        """
        Calculate and return the value of the microcanonical tunneling 
        correction for the reaction at the given energies `e_list` in J/mol.
        """
        if self.tunneling is not None:
            return self.tunneling.calculate_tunneling_function(e_list)
        else:
            # Return step function
            kappa = np.ones_like(e_list)
            E0 = float(self.conformer.E0.value_si)
            for r in range(e_list.shape[0]):
                if e_list[r] >= E0:
                    break
                kappa[r] = 0.0
            return kappa
