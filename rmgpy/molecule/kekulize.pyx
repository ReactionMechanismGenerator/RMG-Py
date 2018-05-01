###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains functions for kekulization of a aromatic molecule.
The only function that should be used outside of this module is the main
`kekulize()` function. The remaining functions and classes are designed only
to support the kekulization algorithm, and should not be used on their own.

The basic algorithm is as follows:
1. Identify all aromatic rings in the molecule, based on bond types.
2. For each ring, identify endocyclic and exocyclic bonds.
3. Determine if any bonds in the ring are already defined (not benzene bonds).
4. For the remaining bonds, determine whether or not they can be double bonds.
5. If a clear determination cannot be made, make heuristic based assumption.
6. Continue until all bonds in the ring are determined.
7. Continue until all rings in the molecule are determined.

Here, `endo` refers to bonds that comprise a given ring, while `exo` refers to
bonds that are connected to atoms in the ring, but not part of the ring itself.

A key part of the algorithm is use of degree of freedom (DOF) analysis in order
to determine the optimal order to solve the system. Rings and bonds with fewer
DOFs have fewer ways to be to be kekulized, and are generally easier to solve.
Each ring or bond that is fixed reduces the DOF of adjacent rings and bonds,
and the process continues until the entire molecule can be solved.
"""

import logging

from .graph cimport Graph
from .molecule cimport Atom, Bond, Molecule
from .element import PeriodicSystem
from rmgpy.exceptions import KekulizationError, AtomTypeError

cpdef kekulize(Graph mol):
    """
    Kekulize an aromatic molecule in place. If the molecule cannot be kekulized,
    a KekulizationError will be raised. However, the molecule will be left in
    a semi-kekulized state. Therefore, if the original molecule needs to be kept,
    it is advisable to create a copy before kekulizing.

    Args: :class:`Molecule` object to be kekulized
    """
    cdef list ring, rings, aromatic_rings, resolved_rings
    cdef set endo_bonds, exo_bonds
    cdef Atom atom1, atom2, atom
    cdef Bond bond
    cdef bint aromatic, successful, bridged
    cdef int itercount, maxiter
    cdef AromaticRing aromatic_ring

    # Get all potentially aromatic rings
    rings = mol.getAllCyclesOfSize(6)

    # Identify aromatic rings and categorize endocyclic and exocyclic bonds for each ring
    aromatic_rings = []
    for ring in rings:
        endo_bonds = set()
        exo_bonds = set()
        aromatic = True
        for atom1 in ring:
            # Check if this is a bridged ring
            bridged = sum([1 if atom in ring else 0 for atom in atom1.bonds.iterkeys()]) > 2
            for atom2, bond in atom1.bonds.iteritems():
                if bridged and sum([1 if atom in ring else 0 for atom in atom2.bonds.iterkeys()]) > 2:
                    # This atom2 is the other end of the bridging bond, so don't consider it as a part of the ring
                    exo_bonds.add(bond)
                    continue
                elif atom2 in ring:
                    if abs(round(bond.order) - bond.order) < 1e-9:
                        aromatic = False
                        break
                    endo_bonds.add(bond)
                else:
                    exo_bonds.add(bond)
            if not aromatic:
                break
        if aromatic:
            # Use an AromaticRing object to store the info about this ring
            aromatic_rings.append(AromaticRing(atoms=ring, endo_bonds=endo_bonds, exo_bonds=exo_bonds))

    resolved_rings = []
    itercount = 0
    maxiter = 2 * len(aromatic_rings)
    while aromatic_rings and itercount < maxiter:
        # Update and sort the remaining rings
        prioritize_rings(aromatic_rings)
        # Take the next ring off the stack
        aromatic_ring = aromatic_rings.pop()
        # Try to kekulize this ring
        successful = aromatic_ring.kekulize()
        if successful:
            resolved_rings.append(aromatic_ring)
        else:
            # Put it back in the list, which will get resorted by DOF in the next iteration
            aromatic_rings.append(aromatic_ring)
        itercount += 1

    if aromatic_rings:
        raise KekulizationError('Unable to kekulize molecule, reached maximum attempts:/n{0}'.format(mol.toAdjacencyList()))

    try:
        mol.updateAtomTypes(logSpecies=False)
    except AtomTypeError:
        logging.debug('Unable to kekulize molecule, final result was invalid:/n{0}'.format(mol.toAdjacencyList()))
        raise KekulizationError('Unable to kekulize molecule, final result was invalid.')

cdef list prioritize_rings(list item_list):
    """Update list of AromaticRing objects, then sort by DOF."""
    cdef AromaticRing item, x
    for item in item_list:
        item.update()
    return item_list.sort(key=lambda x: (x.endo_dof, x.exo_dof), reverse=True)

cdef list prioritize_bonds(list item_list):
    """Update list of Aromatic Bond objects, then sort by DOF."""
    cdef AromaticBond item, x
    for item in item_list:
        item.update()
    return item_list.sort(key=lambda x: (x.double_possible, not x.double_required, x.endo_dof, x.exo_dof), reverse=True)

cdef class AromaticRing(object):
    """
    Helper class containing information about a single aromatic ring in a molecule.

    DO NOT use outside of this module. This class does not do any aromaticity perception.
    """
    cdef public list atoms, resolved, unresolved
    cdef set endo_bonds, exo_bonds
    cdef public int endo_dof, exo_dof

    def __init__(self, atoms=None, endo_bonds=None, exo_bonds=None, endo_dof=-1, exo_dof=-1):
        self.atoms = atoms
        self.endo_bonds = endo_bonds
        self.exo_bonds = exo_bonds
        self.endo_dof = endo_dof
        self.exo_dof = exo_dof
        self.resolved = []
        self.unresolved = []

    cpdef update(self):
        """
        Update the degree of freedom information for this aromatic ring.

        `endo_dof` refers to the number of bonds in the ring without fixed bond orders.
        `exo_dof`  refers to the number of bonds outside the ring without fixed bond orders.
        """
        cdef int endo_dof, exo_dof
        cdef Bond bond

        endo_dof = 0
        for bond in self.endo_bonds:
            if bond.isBenzene():
                # Add one dof for each aromatic bond
                endo_dof += 1
        exo_dof = 0
        for bond in self.exo_bonds:
            if bond.isBenzene():
                # Add one dof for each aromatic bond
                exo_dof += 1
        self.endo_dof = endo_dof
        self.exo_dof = exo_dof

        self.process_bonds()

    cpdef tuple process_bonds(self):
        """Create AromaticBond objects for each endocyclic bond."""
        cdef Bond bond0
        cdef int i

        if not self.unresolved and not self.resolved:
            # We just started on this ring
            for bond0 in self.endo_bonds:
                self.unresolved.append(AromaticBond(bond=bond0, ring_bonds=self.endo_bonds))

        i = 0
        while i < len(self.unresolved):
            bond0 = self.unresolved[i].bond
            if bond0.isOrder(round(bond0.order)):
                # Bond has already been assigned, so mark as resolved
                self.resolved.append(self.unresolved.pop(i))
            elif bond0.isOrder(2.5):
                # Bond was incremented, so it must be a double bond
                bond0.order = 2
                self.resolved.append(self.unresolved.pop(i))
            elif bond0.isOrder(0.5):
                # Bond was decremented, so it must be a single bond
                bond0.order = 1
                self.resolved.append(self.unresolved.pop(i))
            else:
                i += 1

        assert len(self.resolved) + len(self.unresolved) == len(self.endo_bonds)

    cpdef bint kekulize(self) except -2:
        """
        Attempts to kekulize a single aromatic ring in a molecule.

        Returns True if successful, and False otherwise.
        """
        cdef list resolved, unresolved
        cdef int itercount, maxiter
        cdef AromaticBond bond

        # Check status
        if len(self.unresolved) == 0:
            return True

        itercount = 0
        maxiter = 2 * len(self.unresolved)
        while self.unresolved and itercount < maxiter:
            # Update and sort the unresolved bonds
            prioritize_bonds(self.unresolved)
            # Take the next bond off the stack
            bond = self.unresolved.pop()

            if bond.double_possible and bond.double_required:
                # This bond must be a double bond to satisfy atom valence
                bond.bond.order = 2
                self.resolved.append(bond)
                self.endo_dof -= 1
            elif bond.double_possible and not bond.double_required:
                # This could be a double bond, but we don't know for sure
                # There are a few cases where it's safe to assume that it is a double bond:
                #   - All exo bonds are defined, and no endo bonds have been defined
                #   - Exo bonds adjacent to the current bond are defined, and no endo bonds have been defined
                #   - Exo bonds adjacent to the current bond are not defined, but one adjacent endo bond is defined
                #   - This is the last undefined endo bond
                if ((self.endo_dof == 6 and self.exo_dof == 0)
                        or (self.endo_dof == 6 and bond.exo_dof == 0)
                        or (bond.endo_dof == 1 and (bond.exo_dof == 1 or bond.exo_dof == 2))
                        or self.endo_dof == 1):
                    # Go ahead an assume this bond is double
                    bond.bond.order = 2
                    self.resolved.append(bond)
                    self.endo_dof -= 1
                else:
                    # Come back to this bond later
                    self.unresolved.append(bond)
            else:
                # Double bond is not possible, so it must be a single bond
                bond.bond.order = 1
                self.resolved.append(bond)
                self.endo_dof -= 1

            itercount += 1

        if self.unresolved:
            # We've hit the iteration limit, but could not solve the ring
            return False

        return True

cdef class AromaticBond(object):
    """
    Helper class containing information about a single aromatic bond in a molecule.

    DO NOT use outside of this module. This class does not do any aromaticity perception.
    """
    cdef public Bond bond
    cdef public set ring_bonds
    cdef public int endo_dof, exo_dof
    cdef public bint double_possible, double_required

    def __init__(self, bond=None, ring_bonds=None, endo_dof=-1, exo_dof=-1, double_possible=True, double_required=False):
        self.bond = bond
        self.ring_bonds = ring_bonds
        self.endo_dof = endo_dof
        self.exo_dof = exo_dof
        self.double_possible = double_possible
        self.double_required = double_required

    cpdef update(self):
        """
        Update the local degree of freedom information for this aromatic bond.
        The DOF counts do not include the bond itself, only its adjacent bonds.

        `endo_dof` refers to the number of adjacent bonds in the ring without fixed bond orders.
        `exo_dof`  refers to the number of adjacent bonds outside the ring without fixed bond orders.
        """
        cdef dict valences
        cdef Atom atom
        cdef Bond bond
        cdef int endo_dof, exo_dof, occupied, uncertain, available

        valences = PeriodicSystem.valences

        endo_dof = 0
        exo_dof = 0
        for atom in [self.bond.atom1, self.bond.atom2]:
            occupied = 0
            uncertain = 0
            # Count electrons in bonds
            for bond in atom.bonds.itervalues():
                if abs(round(bond.order) - bond.order) < 1e-9:
                    # This is a fixed bond, either single or double
                    occupied += int(round(bond.order))
                elif bond.isBenzene():
                    # The atom has a benzene bond, so at least one electron is occupied, but there is a second uncertain electron
                    occupied += 1
                    uncertain += 1
                    if bond is not self.bond:
                        if bond in self.ring_bonds:
                            endo_dof += 1
                        else:
                            exo_dof += 1
                else:
                    raise KekulizationError('Unexpected bond order {0}.'.format(bond.order))
            # Count radicals and lone pairs
            occupied += atom.radicalElectrons
            occupied += 2 * atom.lonePairs
            # Valence calculation to determine available electrons
            available = valences[atom.element.symbol] - occupied
            if available < 0:
                raise KekulizationError('Atom {0} cannot have negative available valence.'.format(atom))
            elif available == 0:
                # There are no extra electrons available, so this bond cannot be a double bond
                self.double_possible = False
            elif available == 1 and uncertain == 1:
                # There is an extra electron available, but the current bond is the only uncertain one,
                # so it must be a double bond
                self.double_required = True

        self.endo_dof = endo_dof
        self.exo_dof = exo_dof
