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
This module contains methods for generation of resonance isomers of molecules.
"""

import cython

import rmgpy.molecule.generator as generator
import rmgpy.molecule.parser as parser

from .graph import Vertex, Edge, Graph, getVertexConnectivityValue
from .molecule import Atom, Bond, Molecule
from .atomtype import AtomTypeError
import rmgpy.molecule.pathfinder as pathfinder

def populate_resonance_generation_algorithm():
    """
    A list with the current set of resonance generation algorithms.
    """
    algorithms = (
        generateAdjacentResonanceIsomers,
        generateLonePairRadicalResonanceIsomers,
        generateN5dd_N5tsResonanceIsomers,
        generateKekulizedResonanceIsomers,
        generateAromaticResonanceIsomers,
    )

    return algorithms

def generateResonanceIsomers(mol):
    """
    Generate and return all of the resonance isomers of this molecule.
    """
    cython.declare(isomers=list, newIsomers=list, index=cython.int, atom=Atom)
    cython.declare(isomer=Molecule, newIsomer=Molecule, isom=Molecule)
        
    isomers = [mol]

    # Iterate over resonance isomers
    index = 0
    while index < len(isomers):
        isomer = isomers[index]
        
        newIsomers = []
        for algo in populate_resonance_generation_algorithm():
            newIsomers.extend(algo(isomer))

        for newIsomer in newIsomers:
            # Append to isomer list if unique
            for isom in isomers:
                if isom.isIsomorphic(newIsomer):
                    break
            else:
                isomers.append(newIsomer)
    
        # Move to next resonance isomer
        index += 1
    
    return isomers

def generateAdjacentResonanceIsomers(mol):
    """
    Generate all of the resonance isomers formed by one allyl radical shift.
    """
    cython.declare(isomers=list, paths=list, index=cython.int, isomer=Molecule)
    cython.declare(atom=Atom, atom1=Atom, atom2=Atom, atom3=Atom, bond12=Bond, bond23=Bond)
    cython.declare(v1=Vertex, v2=Vertex)
    
    isomers = []

    # Radicals
    if mol.isRadical():
        # Iterate over radicals in structure
        for atom in mol.vertices:
            paths = pathfinder.findAllDelocalizationPaths(atom)
            for atom1, atom2, atom3, bond12, bond23 in paths:
                # Adjust to (potentially) new resonance isomer
                atom1.decrementRadical()
                atom3.incrementRadical()
                bond12.incrementOrder()
                bond23.decrementOrder()
                # Make a copy of isomer
                isomer = mol.copy(deep=True)
                # Also copy the connectivity values, since they are the same
                # for all resonance forms
                for index in range(len(mol.vertices)):
                    v1 = mol.vertices[index]
                    v2 = isomer.vertices[index]
                    v2.connectivity1 = v1.connectivity1
                    v2.connectivity2 = v1.connectivity2
                    v2.connectivity3 = v1.connectivity3
                    v2.sortingLabel = v1.sortingLabel
                # Restore current isomer
                atom1.incrementRadical()
                atom3.decrementRadical()
                bond12.decrementOrder()
                bond23.incrementOrder()
                # Append to isomer list if unique
                isomer.updateAtomTypes(logSpecies=False)
                isomers.append(isomer)

    return isomers

def generateLonePairRadicalResonanceIsomers(mol):
    """
    Generate all of the resonance isomers formed by lone electron pair - radical shifts.
    """
    cython.declare(isomers=list, paths=list, index=cython.int, isomer=Molecule)
    cython.declare(atom=Atom, atom1=Atom, atom2=Atom)
    cython.declare(v1=Vertex, v2=Vertex)
    
    isomers = []

    # Radicals
    if mol.isRadical():
        # Iterate over radicals in structure
        for atom in mol.vertices:
            paths = pathfinder.findAllDelocalizationPathsLonePairRadical(atom)
            for atom1, atom2 in paths:
                # Adjust to (potentially) new resonance isomer
                atom1.decrementRadical()
                atom1.incrementLonePairs()
                atom1.updateCharge()
                atom2.incrementRadical()
                atom2.decrementLonePairs()
                atom2.updateCharge()
                # Make a copy of isomer
                isomer = mol.copy(deep=True)
                # Also copy the connectivity values, since they are the same
                # for all resonance forms
                for index in range(len(mol.vertices)):
                    v1 = mol.vertices[index]
                    v2 = isomer.vertices[index]
                    v2.connectivity1 = v1.connectivity1
                    v2.connectivity2 = v1.connectivity2
                    v2.connectivity3 = v1.connectivity3
                    v2.sortingLabel = v1.sortingLabel
                # Restore current isomer
                atom1.incrementRadical()
                atom1.decrementLonePairs()
                atom1.updateCharge()
                atom2.decrementRadical()
                atom2.incrementLonePairs()
                atom2.updateCharge()
                # Append to isomer list if unique
                isomer.updateAtomTypes(logSpecies=False)
                isomers.append(isomer)

    return isomers

def generateN5dd_N5tsResonanceIsomers(mol):
    """
    Generate all of the resonance isomers formed by shifts between N5dd and N5ts.
    """
    cython.declare(isomers=list, paths=list, index=cython.int, isomer=Molecule)
    cython.declare(atom=Atom, atom1=Atom, atom2=Atom, atom3=Atom)
    cython.declare(bond12=Bond, bond13=Bond)
    cython.declare(v1=Vertex, v2=Vertex)
    
    isomers = []
    
    # Iterate over nitrogen atoms in structure
    for atom in mol.vertices:
        paths = pathfinder.findAllDelocalizationPathsN5dd_N5ts(atom)
        for atom1, atom2, atom3, bond12, bond13, direction in paths:
            # from N5dd to N5ts
            if direction == 1:
                # Adjust to (potentially) new resonance isomer
                bond12.decrementOrder()
                bond13.incrementOrder()
                atom2.incrementLonePairs()
                atom3.decrementLonePairs()
                atom1.updateCharge()
                atom2.updateCharge()
                atom3.updateCharge()
                # Make a copy of isomer
                isomer = mol.copy(deep=True)
                # Also copy the connectivity values, since they are the same
                # for all resonance forms
                for index in range(len(mol.vertices)):
                    v1 = mol.vertices[index]
                    v2 = isomer.vertices[index]
                    v2.connectivity1 = v1.connectivity1
                    v2.connectivity2 = v1.connectivity2
                    v2.connectivity3 = v1.connectivity3
                    v2.sortingLabel = v1.sortingLabel
                # Restore current isomer
                bond12.incrementOrder()
                bond13.decrementOrder()
                atom2.decrementLonePairs()
                atom3.incrementLonePairs()
                atom1.updateCharge()
                atom2.updateCharge()
                atom3.updateCharge()
                # Append to isomer list if unique
                isomer.updateAtomTypes(logSpecies=False)
                isomers.append(isomer)
            
            # from N5ts to N5dd
            if direction == 2:
                # Adjust to (potentially) new resonance isomer
                bond12.decrementOrder()
                bond13.incrementOrder()
                atom2.incrementLonePairs()
                atom3.decrementLonePairs()
                atom1.updateCharge()
                atom2.updateCharge()
                atom3.updateCharge()
                # Make a copy of isomer
                isomer = mol.copy(deep=True)
                # Also copy the connectivity values, since they are the same
                # for all resonance forms
                for index in range(len(mol.vertices)):
                    v1 = mol.vertices[index]
                    v2 = isomer.vertices[index]
                    v2.connectivity1 = v1.connectivity1
                    v2.connectivity2 = v1.connectivity2
                    v2.connectivity3 = v1.connectivity3
                    v2.sortingLabel = v1.sortingLabel
                # Restore current isomer
                bond12.incrementOrder()
                bond13.decrementOrder()
                atom2.decrementLonePairs()
                atom3.incrementLonePairs()
                atom1.updateCharge()
                atom2.updateCharge()
                atom3.updateCharge()
                # Append to isomer list if unique
                isomer.updateAtomTypes(logSpecies=False)
                isomers.append(isomer)
                
    return isomers

def generateAromaticResonanceIsomers(mol):
    """
    Generate the aromatic form of the molecule.
    
    Returns it as a single element of a list.
    If there's an error (eg. in RDKit) it just returns an empty list.
    """
    cython.declare(molecule=Molecule, rdAtomIndices=dict, aromatic=cython.bint, aromaticBonds=list)
    cython.declare(rings=list, ring0=list, i=cython.int, atom1=Atom, atom2=Atom, bond=Bond)
    from rdkit.Chem.rdchem import BondType
    AROMATIC = BondType.AROMATIC

    # Radicals
    if not mol.isCyclic():
        return []

    molecule = mol.copy(deep=True)

    # In RMG, only 6-member rings can be considered aromatic, so ignore all other rings
    rings = [ring0 for ring0 in molecule.getSmallestSetOfSmallestRings() if len(ring0) == 6]
    if not rings:
        return []

    try:
        rdkitmol, rdAtomIndices = generator.toRDKitMol(molecule, removeHs=False, returnMapping=True)
    except ValueError:
        return []
    aromatic = False
    for ring0 in rings:
        aromaticBonds = []
        # Figure out which atoms and bonds are aromatic and reassign appropriately:
        for i, atom1 in enumerate(ring0):
            if not atom1.isCarbon():
                # all atoms in the ring must be carbon in RMG for our definition of aromatic
                break
            for atom2 in ring0[i + 1:]:
                if molecule.hasBond(atom1, atom2):
                    if rdkitmol.GetBondBetweenAtoms(rdAtomIndices[atom1], rdAtomIndices[atom2]).GetBondType() is AROMATIC:
                        aromaticBonds.append(molecule.getBond(atom1, atom2))
        else:  # didn't break so all atoms are carbon
            if len(aromaticBonds) == 6:
                aromatic = True
                # Only change bonds if there are all 6 are aromatic.  Otherwise don't do anything
                for bond in aromaticBonds:
                    bond.order = 1.5

    if aromatic:
        try:
            molecule.updateAtomTypes(logSpecies=False)
        except:
            # Something incorrect has happened, ie. 2 double bonds on a Cb atomtype
            # Do not add the new isomer since it is malformed
            return []
        else:
            # nothing bad happened
            return [molecule]
    else:
        return []

def generateKekulizedResonanceIsomers(mol):
    """
    Generate a kekulized (single-double bond) form of the molecule.
    
    Returns a single Kekule form, as an element of a list of length 1.
    If there's an error (eg. in RDKit) then it just returns an empty list.
    """
    cython.declare(atom=Atom)
    for atom in mol.atoms:
        if atom.atomType.label == 'Cb' or atom.atomType.label == 'Cbf':
            break
    else:
        return []
   
    try:
        rdkitmol = generator.toRDKitMol(mol)  # This perceives aromaticity
        isomer = parser.fromRDKitMol(Molecule(), rdkitmol)  # This step Kekulizes the molecule
    except ValueError:
        return []
    isomer.updateAtomTypes(logSpecies=False)
    return [isomer]

def generate_isomorphic_isomers(mol):
    """
    Select the resonance isomer that is isomorphic to the parameter isomer, with the lowest unpaired
    electrons descriptor.

    We generate over all resonance isomers (non-isomorphic as well as isomorphic) and retain isomorphic
    isomers.

    WIP: do not generate aromatic resonance isomers.
    """

    cython.declare(isomorphic_isomers=list,\
                   isomers=list,
                    )

    cython.declare(isomer=Molecule,\
                   newIsomer=Molecule,\
                   isom=Molecule
                   )

    cython.declare(index=int)

    isomorphic_isomers = [mol]# resonance isomers that are isomorphic to the parameter isomer.

    isomers = [mol]

    # Iterate over resonance isomers
    index = 0
    while index < len(isomers):
        isomer = isomers[index]
        
        newIsomers = []
        for algo in populate_resonance_generation_algorithm():
            newIsomers.extend(algo(isomer))
        
        for newIsomer in newIsomers:
            # Append to isomer list if unique
            for isom in isomers:
                if isom.copy(deep=True).isIsomorphic(newIsomer.copy(deep=True)):
                    isomorphic_isomers.append(newIsomer)
                    break
            else:
                isomers.append(newIsomer)        
                    
        # Move to next resonance isomer
        index += 1

    return isomorphic_isomers


def generateClarStructures(mol):
    """
    Generate Clar structures for a given molecule.

    Returns a list of :class:`Molecule` objects corresponding to the Clar structures.
    """
    cython.declare(output=list, molList=list, newmol=Molecule, asssr=list, bonds=list, solution=list,
                   y=list, x=list, index=cython.int, bond=Bond, ring=list)

    if not mol.isCyclic():
        return []

    output = clarOptimization(mol)

    molList = []

    for newmol, asssr, bonds, solution in output:

        # The solution includes a part corresponding to rings, y, and a part corresponding to bonds, x, using
        # nomenclature from the paper. In y, 1 means the ring as a sextet, 0 means it does not.
        # In x, 1 corresponds to a double bond, 0 either means a single bond or the bond is part of a sextet.
        y = solution[0:len(asssr)]
        x = solution[len(asssr):]

        # Apply results to molecule - double bond locations first
        for index, bond in enumerate(bonds):
            if x[index] == 0:
                bond.order = 1 # single
            elif x[index] == 1:
                bond.order = 2 # double
            else:
                raise ValueError('Unaccepted bond value {0} obtained from optimization.'.format(x[index]))

        # Then apply locations of aromatic sextets by converting to benzene bonds
        for index, ring in enumerate(asssr):
            if y[index] == 1:
                clarTransformation(newmol, ring)

        try:
            newmol.updateAtomTypes()
        except AtomTypeError:
            pass
        else:
            molList.append(newmol)

    return molList


def clarOptimization(mol, constraints=None, maxNum=None):
    """
    Implements linear programming algorithm for finding Clar structures. This algorithm maximizes the number
    of Clar sextets within the constraints of molecular geometry and atom valency.

    Returns a list of valid Clar solutions in the form of a tuple, with the following entries:
        [0] Molecule object
        [1] List of aromatic rings
        [2] List of bonds
        [3] Optimization solution

    The optimization solution is a list of boolean values with sextet assignments followed by double bond assignments,
    with indices corresponding to the list of aromatic rings and list of bonds, respectively.

    Method from:
        Hansen, P.; Zheng, M. The Clar Number of a Benzenoid Hydrocarbon and Linear Programming.
            J. Math. Chem. 1994, 15 (1), 93–107.
    """
    cython.declare(molecule=Molecule, asssr=list, exo=list, l=cython.int, m=cython.int, n=cython.int,
                   a=list, objective=list, status=cython.int, solution=list, innerSolutions=list)

    from lpsolve55 import lpsolve

    # Make a copy of the molecule so we don't destroy the original
    molecule = mol.copy(deep=True)

    asssr = molecule.getAromaticSSSR()

    if not asssr:
        return []

    # Get list of atoms that are in rings
    atoms = set()
    for ring in asssr:
        atoms.update(ring)
    atoms = list(atoms)

    # Get list of bonds involving the ring atoms, ignoring bonds to hydrogen
    bonds = set()
    for atom in atoms:
        bonds.update([atom.bonds[key] for key in atom.bonds.keys() if key.isNonHydrogen()])
    bonds = list(bonds)

    # Identify exocyclic bonds, and save their bond orders
    exo = []
    for bond in bonds:
        if bond.atom1 not in atoms or bond.atom2 not in atoms:
            if bond.isDouble():
                exo.append(1)
            else:
                exo.append(0)
        else:
            exo.append(None)

    # Dimensions
    l = len(asssr)
    m = len(atoms)
    n = l + len(bonds)

    # Connectivity matrix which indicates which rings and bonds each atom is in
    # Part of equality constraint Ax=b
    a = []
    for atom in atoms:
        inRing = [1 if atom in ring else 0 for ring in asssr]
        inBond = [1 if atom in [bond.atom1, bond.atom2] else 0 for bond in bonds]
        a.append(inRing + inBond)

    # Objective vector for optimization: sextets have a weight of 1, double bonds have a weight of 0
    objective = [1] * l + [0] * len(bonds)

    # Solve LP problem using lpsolve
    lp = lpsolve('make_lp', m, n)               # initialize lp with constraint matrix with m rows and n columns
    lpsolve('set_verbose', lp, 2)               # reduce messages from lpsolve
    lpsolve('set_obj_fn', lp, objective)        # set objective function
    lpsolve('set_maxim', lp)                    # set solver to maximize objective
    lpsolve('set_mat', lp, a)                   # set left hand side to constraint matrix
    lpsolve('set_rh_vec', lp, [1] * m)          # set right hand side to 1 for all constraints
    lpsolve('set_constr_type', lp, ['='] * m)   # set all constraints as equality constraints
    lpsolve('set_binary', lp, [True] * n)       # set all variables to be binary

    # Constrain values of exocyclic bonds, since we don't want to modify them
    for i in range(l, n):
        if exo[i - l] is not None:
            # NOTE: lpsolve indexes from 1, so the variable we're changing should be i + 1
            lpsolve('set_bounds', lp, i + 1, exo[i - l], exo[i - l])

    # Add constraints to problem if provided
    if constraints is not None:
        for constraint in constraints:
            lpsolve('add_constraint', lp, constraint[0], '<=', constraint[1])

    status = lpsolve('solve', lp)
    objVal, solution = lpsolve('get_solution', lp)[0:2]
    lpsolve('delete_lp', lp)  # Delete the LP problem to clear up memory

    # Check that optimization was successful
    if status != 0:
        raise ILPSolutionError('Optimization could not find a valid solution.')

    # Check that we the result contains at least one aromatic sextet
    if objVal == 0:
        return []

    # Check that the solution contains the maximum number of sextets possible
    if maxNum is None:
        maxNum = objVal  # This is the first solution, so the result should be an upper limit
    elif objVal < maxNum:
        raise ILPSolutionError('Optimization obtained a sub-optimal solution.')

    if any([x != 1 and x != 0 for x in solution]):
        raise ILPSolutionError('Optimization obtained a non-integer solution.')

    # Generate constraints based on the solution obtained
    y = solution[0:l]
    new_a = y + [0] * len(bonds)
    new_b = sum(y) - 1
    if constraints is not None:
        constraints.append((new_a, new_b))
    else:
        constraints = [(new_a, new_b)]

    # Run optimization with additional constraints
    try:
        innerSolutions = clarOptimization(mol, constraints=constraints, maxNum=maxNum)
    except ILPSolutionError:
        innerSolutions = []

    return innerSolutions + [(molecule, asssr, bonds, solution)]


def clarTransformation(mol, aromaticRing):
    """
    Performs Clar transformation for given ring in a molecule, ie. conversion to aromatic sextet.

    Args:
        mol             a :class:`Molecule` object
        aromaticRing    a list of :class:`Atom` objects corresponding to an aromatic ring in mol

    This function directly modifies the input molecule and does not return anything.
    """
    cython.declare(bondList=list, i=cython.int, atom1=Atom, atom2=Atom, bond=Bond)

    bondList = []

    for i, atom1 in enumerate(aromaticRing):
        for atom2 in aromaticRing[i + 1:]:
            if mol.hasBond(atom1, atom2):
                bondList.append(mol.getBond(atom1, atom2))

    for bond in bondList:
        bond.order = 1.5


class ILPSolutionError(Exception):
    """
    An exception to be raised when solving an integer linear programming problem if a solution
    could not be found or the solution is not valid. Can pass a string to indicate the reason
    that the solution is invalid.
    """
    pass
