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
                    bond.order = 'B'

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

    Returns all Clar structures as a list.
    """
    cython.declare(output=list, molList=list, o=tuple, y=list, x=list, index=cython.int, bond=Bond, ring=list)

    output = clarOptimization(mol)

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


def clarOptimization(mol, constraints=None, maxNum=None):
    """
    Implements linear programming algorithm for finding Clar structures. This algorithm maximizes the number
    of Clar sextets within the constraints of molecular geometry and atom valency.

    Method from:
        Hansen, P.; Zheng, M. The Clar Number of a Benzenoid Hydrocarbon and Linear Programming.
            J. Math. Chem. 1994, 15 (1), 93–107.
    """
    cython.declare(molecule=Molecule, SSSR=list, a=list, objective=list, solution=list, innerSolutions=list)

    import os, sys
    import glpk

    # Make a copy of the molecule so we don't destroy the original
    molecule = mol.copy(deep=True)

    SSSR = getAromaticSSSR(molecule)

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

    # Hack to suppress glpk output, as an alternative to modifying glpk
    with open(os.devnull, 'w') as dn:
        try:
            fd = sys.stdout.fileno()
        except IOError:  # Most likely stdout does not have a file descriptor (eg. running IPython notebook)
            msg = lp.intopt()  # Solve the linear program
        else:
            original = os.dup(fd)  # Save a copy of original file descriptor
            os.dup2(dn.fileno(), fd)  # Set stdout to devnull
            msg = lp.intopt()  # Solve the linear program
            os.dup2(original, fd)  # Restore stdout
            os.close(original)

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
        innerSolutions = clarOptimization(mol, constraints=constraints, maxNum=maxNum)
    except (ValueError, RuntimeError):
        innerSolutions = []

    return innerSolutions + [(molecule, SSSR, bonds, solution)]


def getAromaticSSSR(mol):
    """
    Returns the smallest set of smallest aromatic rings
    """
    cython.declare(rdAtomIndices=dict, aromaticRings=list, aromaticBonds=list)
    cython.declare(rings=list, ring0=list, i=cython.int, atom1=Atom, atom2=Atom)

    from rdkit.Chem.rdchem import BondType

    AROMATIC = BondType.AROMATIC

    rings = [ring0 for ring0 in mol.getSmallestSetOfSmallestRings() if len(ring0) == 6]
    if not rings:
        return []

    try:
        rdkitmol, rdAtomIndices = generator.toRDKitMol(mol, removeHs=False, returnMapping=True)
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
                if mol.hasBond(atom1, atom2):
                    if rdkitmol.GetBondBetweenAtoms(rdAtomIndices[atom1],
                                                    rdAtomIndices[atom2]).GetBondType() is AROMATIC:
                        aromaticBonds.append(mol.getBond(atom1, atom2))
        else:  # didn't break so all atoms are carbon
            if len(aromaticBonds) == 6:
                aromaticRings.append(ring0)

    return aromaticRings


def clarTransformation(mol, ring):
    """
    Performs Clar transformation for given ring in a molecule, ie. conversion to aromatic sextet.
    """
    cython.declare(indexList=list, bondList=list, index1=cython.int, index2=cython.int, bond=Bond)

    indexList = zip(range(len(ring)), range(1, len(ring)) + [0])

    bondList = []
    for index1, index2 in indexList:
        try:
            bondList.append(mol.getBond(ring[index1], ring[index2]))
        except ValueError:
            raise Exception('Atoms in ring not in connected order.')

    for bond in bondList:
        bond.order = 'B'

    try:
        mol.updateAtomTypes()
    except:
        return []

    return [mol]
