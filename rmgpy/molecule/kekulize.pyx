# encoding: utf-8
# cython: linetrace=True

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains functions for kekulization of a aromatic molecule.

The basic algorithm is as follows:
1. Identify all aromatic rings in the molecule, based on bond types.
2. For each ring, identify endocyclic and exocyclic bonds.
3. Determine if any bonds in the ring are already defined (not benzene bonds).
4. For the remaining bonds, determine whether or not they can be double bonds.
5. If a clear determination cannot be made, make heuristic based assumption.
6. Continue until all bonds in the ring are determined.
7. Continue until all rings in the molecule are determined.
"""

from .molecule cimport Atom, Bond, Molecule
from .element import PeriodicSystem


cpdef kekulize(Molecule mol):
    """Kekulize an aromatic molecule."""
    cdef list ring, rings, aromaticRings, resolvedRings
    cdef set endoBonds, exoBonds
    cdef Atom atom1, atom2
    cdef Bond bond
    cdef bint aromatic, successful
    cdef int itercount, maxiter
    cdef AromaticRing aromaticRing

    # Get all potentially aromatic rings
    rings = mol.getAllCyclesOfSize(6)

    # Identify aromatic rings and categorize endocyclic and exocyclic bonds for each ring
    aromaticRings = []
    for ring in rings:
        endoBonds = set()
        exoBonds = set()
        aromatic = True
        for atom1 in ring:
            for atom2, bond in atom1.bonds.iteritems():
                if atom2 in ring:
                    if abs(round(bond.order) - bond.order) < 1e-9:
                        aromatic = False
                        break
                    endoBonds.add(bond)
                else:
                    exoBonds.add(bond)
            if not aromatic:
                break
        if aromatic:
            aromaticRings.append(AromaticRing(atoms=ring, endoBonds=endoBonds, exoBonds=exoBonds))

    resolvedRings = []
    itercount = 0
    maxiter = 2 * len(aromaticRings)
    while aromaticRings and itercount < maxiter:
        # Update and sort the remaining rings
        prioritizeRings(aromaticRings)
        # Take the next ring off the stack
        aromaticRing = aromaticRings.pop()
        # Try to kekulize this ring
        successful = aromaticRing.kekulize()
        if successful:
            resolvedRings.append(aromaticRing)
        else:
            aromaticRings.append(aromaticRing)
        itercount += 1

    mol.updateAtomTypes(logSpecies=False)

cdef list prioritizeRings(list aromaticList):
    """Update list of AromaticRing objects, then sort by DOF."""
    cdef AromaticRing item, x
    for item in aromaticList:
        item.update()
    return aromaticList.sort(key=lambda x: (x.endoDOF, x.exoDOF), reverse=True)

cdef list prioritizeBonds(list aromaticList):
    """Update list of Aromatic Bond objects, then sort by DOF."""
    cdef AromaticBond item, x
    for item in aromaticList:
        item.update()
    return aromaticList.sort(key=lambda x: (x.endoDOF, x.exoDOF), reverse=True)

cdef class AromaticRing:
    """Helper class containing information about a single aromatic ring in a molecule."""
    cdef list atoms
    cdef set endoBonds, exoBonds
    cdef public int endoDOF, exoDOF

    def __init__(self, atoms=None, endoBonds=None, exoBonds=None, endoDOF=-1, exoDOF=-1):
        self.atoms = atoms
        self.endoBonds = endoBonds
        self.exoBonds = exoBonds
        self.endoDOF = endoDOF
        self.exoDOF = exoDOF

    cdef update(self):
        """Update the degree of freedom information for this aromatic ring."""
        cdef int endoDOF, exoDOF
        cdef Bond bond

        endoDOF = 0
        for bond in self.endoBonds:
            if bond.isBenzene():
                # Add one dof for each aromatic bond
                endoDOF += 1
        exoDOF = 0
        for bond in self.exoBonds:
            if bond.isBenzene():
                # Add one dof for each aromatic bond
                exoDOF += 1
        self.endoDOF = endoDOF
        self.exoDOF = exoDOF

    cdef tuple processBonds(self):
        """Create AromaticBond objects for each endocyclic bond."""
        cdef list resolved, unresolved
        cdef Bond bond0
        cdef AromaticBond aromaticBond

        resolved = []
        unresolved = []
        for bond0 in self.endoBonds:
            aromaticBond = AromaticBond(bond=bond0, ringBonds=self.endoBonds)

            if abs(round(bond0.order) - bond0.order) < 1e-9:
                # Bond has already been assigned, so mark as resolved
                resolved.append(aromaticBond)
            elif bond0.isOrder(2.5):
                # Bond was incremented, so it must be a double bond
                bond0.order = 2
                resolved.append(aromaticBond)
            elif bond0.isOrder(0.5):
                # Bond was decremented, so it must be a single bond
                bond0.order = 1
                resolved.append(aromaticBond)
            else:
                unresolved.append(aromaticBond)

        assert len(resolved) + len(unresolved) == len(self.endoBonds)

        return resolved, unresolved

    cdef bint kekulize(self) except -2:
        """
        Attempts to kekulize a single aromatic ring in a molecule.

        Returns True if successful, and False otherwise.
        """
        cdef list resolved, unresolved
        cdef int itercount, maxiter
        cdef AromaticBond bond

        resolved, unresolved = self.processBonds()

        # Check status
        if len(unresolved) == 0:
            return True

        itercount = 0
        maxiter = 2 * len(unresolved)
        while unresolved and itercount < maxiter:
            # Update and sort the unresolved bonds
            prioritizeBonds(unresolved)
            # Take the next bond off the stack
            bond = unresolved.pop()

            if bond.doublePossible and bond.doubleRequired:
                bond.bond.order = 2
                resolved.append(bond)
                self.endoDOF -= 1
            elif bond.doublePossible and not bond.doubleRequired:
                # This could be a double bond, but we don't know for sure
                # There are a few cases where it's safe to assume that it is a double bond:
                #   - All exo bonds are defined, and no endo bonds have been defined
                #   - Exo bonds adjacent to the current bond are defined, and no endo bonds have been defined
                #   - Exo bonds adjacent to the current bond are not defined, but one adjacent endo bond is defined
                #   - This is the last undefined endo bond
                if ((self.endoDOF == 6 and self.exoDOF == 0)
                        or (self.endoDOF == 6 and bond.exoDOF == 0)
                        or (bond.endoDOF == 1 and bond.exoDOF == 2)
                        or self.endoDOF == 1):
                    # Go ahead an assume this bond is double
                    bond.bond.order = 2
                    resolved.append(bond)
                    self.endoDOF -= 1
                else:
                    # Come back to this bond later
                    unresolved.append(bond)
            else:
                # Double bond is not possible, so it must be a single bond
                bond.bond.order = 1
                resolved.append(bond)
                self.endoDOF -= 1

            itercount += 1

        if unresolved:
            # We've hit the iteration limit, but could not solve the ring
            return False

        return True

cdef class AromaticBond:
    """Helper class containing information about a single aromatic bond in a molecule."""
    cdef Bond bond
    cdef set ringBonds
    cdef public int endoDOF, exoDOF
    cdef public bint doublePossible, doubleRequired

    def __init__(self, bond=None, ringBonds=None, endoDOF=-1, exoDOF=-1, doublePossible=True, doubleRequired=False):
        self.bond = bond
        self.ringBonds = ringBonds
        self.endoDOF = endoDOF
        self.exoDOF = exoDOF
        self.doublePossible = doublePossible
        self.doubleRequired = doubleRequired

    cdef update(self):
        """Update the local degree of freedom information for this aromatic bond."""
        cdef dict bondOrders, valences
        cdef Atom atom
        cdef Bond bond
        cdef int endoDOF, exoDOF, occupied, uncertain, available

        valences = PeriodicSystem.valences

        endoDOF = 0
        exoDOF = 0
        for atom in [self.bond.atom1, self.bond.atom2]:
            occupied = 0
            uncertain = 0
            for bond in atom.bonds.itervalues():
                if abs(round(bond.order) - bond.order) < 1e-9:
                    occupied += int(round(bond.order))
                elif bond.isBenzene():
                    # The atom has a benzene bond, so at least one electron is occupied, but there is a second uncertain electron
                    occupied += 1
                    uncertain += 1
                    if bond is not self.bond:
                        if bond in self.ringBonds:
                            endoDOF += 1
                        else:
                            exoDOF += 1
                else:
                    ValueError('Unexpected bond order {0}.'.format(bond.order))
            occupied += atom.radicalElectrons
            occupied += 2 * atom.lonePairs
            available = valences[atom.element.symbol] - occupied
            if available < 0:
                ValueError('Atom {0} cannot have negative available valence.'.format(atom))
            elif available == 0:
                # There are no extra electrons available, so this bond cannot be a double bond
                self.doublePossible = False
            elif available == 1 and uncertain == 1:
                # There is an extra electron available, but the current bond is the only uncertain one,
                # so it must be a double bond
                self.doubleRequired = True

        self.endoDOF = endoDOF
        self.exoDOF = exoDOF

