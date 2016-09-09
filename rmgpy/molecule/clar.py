#!/usr/bin/env python
# encoding: utf-8

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
Clar structure generation
"""

import numpy as np
from scipy.optimize import linprog


def generateClarStructures(molecule):

    output = clarOptimization(molecule)

    newmol = output[0]
    SSSR = output[1]
    bonds = output[2]
    result = output[3]

    # The solution includes a part corresponding to rings, y, and a part corresponding to bonds, x, using nomenclature
    # from the paper. In y, 1 means the ring as a sextet, 0 means it does not. In x, 1 corresponds to a double bond,
    # 0 either means a single bond or the bond is part of a sextet.
    y = result.x[0:len(SSSR)]
    x = result.x[len(SSSR):]

    # Apply results to molecule - double bond locations first
    for index, bond in enumerate(bonds):
        if x[index] == 0:
            bond.order = 'S'
        elif x[index] == 1:
            bond.order = 'D'
        else:
            raise ValueError('Unaccepted bond value obtained from optimization.')

    # Then apply locations of aromatic sextets by converting to benzene bonds
    for index, ring in enumerate(SSSR):
        if y[index] == 1:
            clarTransformation(newmol, ring)

    return newmol


def clarOptimization(molecule, constraints=None):
    """
    Generates Clar structures for a given molecule using linear programming. This algorithm maximizes the number
    of Clar sextets within the constraints of molecular geometry and atom valency.

    Method from:
        Hansen, P.; Zheng, M. The Clar Number of a Benzenoid Hydrocarbon and Linear Programming.
            J. Math. Chem. 1994, 15 (1), 93–107.
    """

    # Make a copy of the molecule so we don't destroy the original
    mol = molecule.copy(deep=True)

    SSSR = getAromaticSSSR(mol)

    # Get list of atoms that are in rings
    atoms = set()
    for ring in SSSR:
        atoms.update(ring)
    atoms = list(atoms)

    # Get list of bonds involving the ring atoms, ignoring bonds to hydrogen
    bonds = set()
    for atom in atoms:
        bonds.update([atom.bonds[key] for key in atom.bonds.keys() if key.isNonHydrogen()])
    bonds = list(bonds)

    # Connectivity matrix which indicates which rings and bonds each atom is in
    # Part of equality constraint Ax=b
    a = []
    for atom in atoms:
        inRing = [1 if atom in ring else 0 for ring in SSSR]
        inBond = [1 if atom in [bond.atom1, bond.atom2] else 0 for bond in bonds]
        a.append(inRing + inBond)
    a = np.array(a)
    # Each atom can only be a part of one double bond or one sextet
    # Other part of equality constraint
    b = np.ones(len(atoms))
    # Weighting vector for optimization: sextets have a weight of -1, double bonds have a weight of 0
    c = np.concatenate((-np.ones(len(SSSR)), np.zeros(len(bonds))))
    # Bounds on variables
    bounds = ((0, 1))

    # Solve linear program.
    result = linprog(c, A_eq=a, b_eq=b, bounds=bounds)

    # Check that optimization was successful
    if result.status != 0:
        raise Exception('Optimization exited with message: {0}'.format(result.message))
    if np.any(np.logical_and(np.not_equal(result.x, 1), np.not_equal(result.x, 0))):
        raise ValueError('Non-integer solution obtained from optimization.')

    return (mol, SSSR, bonds, result)


def getAromaticSSSR(molecule):
    """
    Returns the smallest set of smallest aromatic rings
    """
    import rmgpy.molecule.generator as generator
    from rdkit.Chem.rdchem import BondType
    AROMATIC = BondType.AROMATIC

    rings = [ring0 for ring0 in molecule.getSmallestSetOfSmallestRings() if len(ring0) == 6]
    if not rings:
        return []

    try:
        rdkitmol, rdAtomIndices = generator.toRDKitMol(molecule, removeHs=False, returnMapping=True)
    except ValueError:
        return []

    aromaticRings = []
    for ring0 in rings:
        aromaticBonds = []
        # Figure out which atoms and bonds are aromatic and reassign appropriately:
        for i, atom1 in enumerate(ring0):
            if not atom1.isCarbon():
                # all atoms in the ring must be carbon in RMG for our definition of aromatic
                break
            for atom2 in ring0[i + 1:]:
                if molecule.hasBond(atom1, atom2):
                    if rdkitmol.GetBondBetweenAtoms(rdAtomIndices[atom1],
                                                    rdAtomIndices[atom2]).GetBondType() is AROMATIC:
                        aromaticBonds.append(molecule.getBond(atom1, atom2))
        else:  # didn't break so all atoms are carbon
            if len(aromaticBonds) == 6:
                aromaticRings.append(ring0)

    return aromaticRings


def clarTransformation(molecule, ring):
    """
    Performs Clar transformation for given ring in a molecule, ie. conversion to aromatic sextet.
    """
    indexList = zip(range(len(ring)), range(1, len(ring)) + [0])

    bondList = []
    for index1, index2 in indexList:
        try:
            bondList.append(molecule.getBond(ring[index1], ring[index2]))
        except ValueError:
            raise Exception('Atoms in ring not in connected order.')

    for bond in bondList:
        bond.order = 'B'

    try:
        molecule.updateAtomTypes()
    except:
        return []

    return molecule
