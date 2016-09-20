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

import glpk

def generateClarStructures(molecule):

    output = clarOptimization(molecule)

    molList = []

    for o in output:
        newmol = o[0]
        SSSR = o[1]
        bonds = o[2]
        solution = o[3]

        # The solution includes a part corresponding to rings, y, and a part corresponding to bonds, x, using
        # nomenclature from the paper. In y, 1 means the ring as a sextet, 0 means it does not.
        # In x, 1 corresponds to a double bond, 0 either means a single bond or the bond is part of a sextet.
        y = solution[0:len(SSSR)]
        x = solution[len(SSSR):]

        # Apply results to molecule - double bond locations first
        for index, bond in enumerate(bonds):
            if x[index] == 0:
                bond.order = 'S'
            elif x[index] == 1:
                bond.order = 'D'
            else:
                raise ValueError('Unaccepted bond value {0} obtained from optimization.'.format(x[index]))

        # Then apply locations of aromatic sextets by converting to benzene bonds
        for index, ring in enumerate(SSSR):
            if y[index] == 1:
                clarTransformation(newmol, ring)

        molList.append(newmol)

    return molList


def clarOptimization(molecule, constraints=None, maxNum=None):
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
        a += (inRing + inBond) # Intentionally created as 1D list
    # Weighting vector for optimization: sextets have a weight of 1, double bonds have a weight of 0
    objective = [1] * len(SSSR) + [0] * len(bonds)

    # Initialize LP problem using pyglpk
    lp = glpk.LPX()
    lp.obj.maximize = True

    # Each row is a constraint, so we have one constraint for each atom
    lp.rows.add(len(atoms))
    for r in lp.rows:
        r.bounds = 1  # Constraint must be 1 since each atom can only be part of one sextet or double bond

    # Each column is a variable, corresponding to each ring and bond (that is part of a ring) in the molecule
    lp.cols.add(len(SSSR) + len(bonds))
    for c in lp.cols:
        c.kind = bool

    # Add constraints to problem if provided
    if constraints:
        first = lp.rows.add(len(constraints))
        for index, constraint in enumerate(constraints):
            a += constraint[0]
            lp.rows[first + index].bounds = constraint[1]

    lp.obj[:] = objective  # Set objective coefficients
    lp.matrix = a  # Set constraint coefficients, coefficients ordered left to right, top to bottom
    msg = lp.intopt()  # Solve the linear program

    # Check that optimization was successful
    if lp.status != 'opt':
        raise RuntimeError('Optimization exited with status {0} and message {1}'.format(lp.status, msg))

    if maxNum is None:
        maxNum = lp.obj.value  # This is the first solution, without constraints, so the result should be an upper limit
    elif lp.obj.value < maxNum:
        raise ValueError('Sub-optimal solution obtained.')  # We don't want solutions with fewer sextets

    solution = [c.value for c in lp.cols]

    if any([x != 1 and x != 0 for x in solution]):
        raise ValueError('Non-integer solution obtained from optimization.')

    # Generate constraints based on the solution obtained
    y = solution[0:len(SSSR)]
    new_a = y + [0] * len(bonds)
    new_b = (0, sum(y) - 1)
    if constraints:
        constraints.append((new_a, new_b))
    else:
        constraints = [(new_a, new_b)]

    # Run optimization with additional constraints
    try:
        innerSolutions = clarOptimization(molecule, constraints=constraints, maxNum=maxNum)
    except (ValueError, RuntimeError):
        innerSolutions = []

    return innerSolutions + [(mol, SSSR, bonds, solution)]


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
