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

# this variable is used to name atom IDs so that there are as few conflicts by 
# using the entire space of integer objects
atom_id_counter = -2**15
