# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
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

from .molecule cimport Atom, Bond, Molecule
from .element import PeriodicSystem
from rmgpy.exceptions import KekulizationError, AtomTypeError

cpdef kekulize(Molecule mol):
    """
    Kekulize an aromatic molecule in place. If the molecule cannot be kekulized,
    an AtomTypeError will be raised. However, the molecule will be left in
    a semi-kekulized state. Therefore, if the original molecule needs to be kept,
    it is advisable to create a copy before kekulizing.

    Args: :class:`Molecule` object to be kekulized
    """
    cdef list ring, rings, aromaticRings, resolvedRings
    cdef set endoBonds, exoBonds
    cdef Atom atom1, atom2, atom
    cdef Bond bond
    cdef bint aromatic, successful, bridged
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
            # Check if this is a bridged ring
            bridged = sum([1 if atom in ring else 0 for atom in atom1.bonds.iterkeys()]) > 2
            for atom2, bond in atom1.bonds.iteritems():
                if bridged and sum([1 if atom in ring else 0 for atom in atom2.bonds.iterkeys()]) > 2:
                    # This atom2 is the other end of the bridging bond, so don't consider it as a part of the ring
                    exoBonds.add(bond)
                    continue
                elif atom2 in ring:
                    if abs(round(bond.order) - bond.order) < 1e-9:
                        aromatic = False
                        break
                    endoBonds.add(bond)
                else:
                    exoBonds.add(bond)
            if not aromatic:
                break
        if aromatic:
            # Use an AromaticRing object to store the info about this ring
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
            # Put it back in the list, which will get resorted by DOF in the next iteration
            aromaticRings.append(aromaticRing)
        itercount += 1

    try:
        mol.updateAtomTypes(logSpecies=False)
    except AtomTypeError:
        logging.debug('Unable to kekulize molecule, final result was invalid:/n{0}'.format(mol.toAdjacencyList()))
        raise KekulizationError('Unable to kekulize molecule, final result was invalid.')

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
    return aromaticList.sort(key=lambda x: (x.doublePossible, not x.doubleRequired, x.endoDOF, x.exoDOF), reverse=True)

cdef class AromaticRing(object):
    """
    Helper class containing information about a single aromatic ring in a molecule.

    DO NOT use outside of this module. This class does not do any aromaticity perception.
    """
    cdef list atoms
    cdef set endoBonds, exoBonds
    cdef public int endoDOF, exoDOF

    def __init__(self, atoms=None, endoBonds=None, exoBonds=None, endoDOF=-1, exoDOF=-1):
        self.atoms = atoms
        self.endoBonds = endoBonds
        self.exoBonds = exoBonds
        self.endoDOF = endoDOF
        self.exoDOF = exoDOF

    cpdef update(self):
        """
        Update the degree of freedom information for this aromatic ring.

        `endoDOF` refers to the number of bonds in the ring without fixed bond orders.
        `exoDOF`  refers to the number of bonds outside the ring without fixed bond orders.
        """
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

    cpdef tuple processBonds(self):
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

    cpdef bint kekulize(self) except -2:
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
                # This bond must be a double bond to satisfy atom valence
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
                        or (bond.endoDOF == 1 and (bond.exoDOF == 1 or bond.exoDOF == 2))
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

cdef class AromaticBond(object):
    """
    Helper class containing information about a single aromatic bond in a molecule.

    DO NOT use outside of this module. This class does not do any aromaticity perception.
    """
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

    cpdef update(self):
        """
        Update the local degree of freedom information for this aromatic bond.
        The DOF counts do not include the bond itself, only its adjacent bonds.

        `endoDOF` refers to the number of adjacent bonds in the ring without fixed bond orders.
        `exoDOF`  refers to the number of adjacent bonds outside the ring without fixed bond orders.
        """
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
                        if bond in self.ringBonds:
                            endoDOF += 1
                        else:
                            exoDOF += 1
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
                self.doublePossible = False
            elif available == 1 and uncertain == 1:
                # There is an extra electron available, but the current bond is the only uncertain one,
                # so it must be a double bond
                self.doubleRequired = True

        self.endoDOF = endoDOF
        self.exoDOF = exoDOF
