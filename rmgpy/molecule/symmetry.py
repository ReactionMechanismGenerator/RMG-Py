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
import rmgpy.molecule
"""
This module provides functionality for estimating the symmetry number of a
molecule from its chemical graph representation.
"""

def calculateAtomSymmetryNumber(molecule, atom):
    """
    Return the symmetry number centered at `atom` in the structure. The
    `atom` of interest must not be in a cycle.
    """
    symmetryNumber = 1

    single = 0; double = 0; triple = 0; benzene = 0
    numNeighbors = 0
    for bond in atom.edges.values():
        if bond.isSingle(): single += 1
        elif bond.isDouble(): double += 1
        elif bond.isTriple(): triple += 1
        elif bond.isBenzene(): benzene += 1
        numNeighbors += 1
    
    # If atom has zero or one neighbors, the symmetry number is 1
    if numNeighbors < 2: return symmetryNumber

    # Create temporary structures for each functional group attached to atom
    molecule0 = molecule
    molecule = molecule0.copy(True)
    atom = molecule.vertices[molecule0.vertices.index(atom)]
    molecule.removeAtom(atom)
    groups = molecule.split()

    # Determine equivalence of functional groups around atom
    groupIsomorphism = dict([(group, dict()) for group in groups])
    for group1 in groups:
        for group2 in groups:
            if group1 is not group2 and group2 not in groupIsomorphism[group1]:
                groupIsomorphism[group1][group2] = group1.isIsomorphic(group2)
                groupIsomorphism[group2][group1] = groupIsomorphism[group1][group2]
            elif group1 is group2:
                groupIsomorphism[group1][group1] = True
    count = [sum([int(groupIsomorphism[group1][group2]) for group2 in groups]) for group1 in groups]
    for i in range(count.count(2) / 2):
        count.remove(2)
    for i in range(count.count(3) / 3):
        count.remove(3); count.remove(3)
    for i in range(count.count(4) / 4):
        count.remove(4); count.remove(4); count.remove(4)
    count.sort(); count.reverse()
    
    if atom.radicalElectrons == 0:
        if single == 4:
            # Four single bonds
            if count == [4]: symmetryNumber *= 12
            elif count == [3, 1]: symmetryNumber *= 3
            elif count == [2, 2]: symmetryNumber *= 2
            elif count == [2, 1, 1]: symmetryNumber *= 1
            elif count == [1, 1, 1, 1]: symmetryNumber *= 1
        elif single == 3:
            # Three single bonds
            if count == [3]: symmetryNumber *= 3
            elif count == [2, 1]: symmetryNumber *= 1
            elif count == [1, 1, 1]: symmetryNumber *= 1

        elif single == 2:
            # Two single bonds
            if count == [2]: symmetryNumber *= 2
        elif double == 2:
            # Two double bonds
            if count == [2]: symmetryNumber *= 2
    elif atom.radicalElectrons == 1:
        if single == 3:
            # Three single bonds
            if count == [3]: symmetryNumber *= 6
            elif count == [2, 1]: symmetryNumber *= 2
            elif count == [1, 1, 1]: symmetryNumber *= 1
    elif atom.radicalElectrons == 2:
        if single == 2:
            # Two single bonds
            if count == [2]:
                symmetryNumber *= 2
    
    if atom.isNitrogen():
        for groupN in groups:
            if groupN.toSMILES() == "[N+](=O)[O-]":
                symmetryNumber *= 2
    
    return symmetryNumber

################################################################################

def calculateBondSymmetryNumber(molecule, atom1, atom2):
    """
    Return the symmetry number centered at `bond` in the structure.
    """
    bond = atom1.edges[atom2]
    symmetryNumber = 1
    if bond.isSingle() or bond.isDouble() or bond.isTriple():
        if atom1.equivalent(atom2):
            # An O-O bond is considered to be an "optical isomer" and so no
            # symmetry correction will be applied
            if atom1.atomType.label == 'Os' and atom2.atomType.label == 'Os' and atom1.radicalElectrons == atom2.radicalElectrons == 0:
                return symmetryNumber
            # If the molecule is diatomic, then we don't have to check the
            # ligands on the two atoms in this bond (since we know there
            # aren't any)
            elif len(molecule.vertices) == 2:
                symmetryNumber = 2
            else:
                molecule.removeBond(bond)
                structure = molecule.copy(True)
                molecule.addBond(bond)

                atom1 = structure.atoms[molecule.atoms.index(atom1)]
                atom2 = structure.atoms[molecule.atoms.index(atom2)]
                fragments = structure.split()
                
                if len(fragments) != 2: return symmetryNumber

                fragment1, fragment2 = fragments
                if atom1 in fragment1.atoms: fragment1.removeAtom(atom1)
                if atom2 in fragment1.atoms: fragment1.removeAtom(atom2)
                if atom1 in fragment2.atoms: fragment2.removeAtom(atom1)
                if atom2 in fragment2.atoms: fragment2.removeAtom(atom2)
                groups1 = fragment1.split()
                groups2 = fragment2.split()

                # Test functional groups for symmetry
                if len(groups1) == len(groups2) == 1:
                    if groups1[0].isIsomorphic(groups2[0]): symmetryNumber *= 2
                elif len(groups1) == len(groups2) == 2:
                    if groups1[0].isIsomorphic(groups2[0]) and groups1[1].isIsomorphic(groups2[1]): symmetryNumber *= 2
                    elif groups1[1].isIsomorphic(groups2[0]) and groups1[0].isIsomorphic(groups2[1]): symmetryNumber *= 2
                elif len(groups1) == len(groups2) == 3:
                    if groups1[0].isIsomorphic(groups2[0]) and groups1[1].isIsomorphic(groups2[1]) and groups1[2].isIsomorphic(groups2[2]): symmetryNumber *= 2
                    elif groups1[0].isIsomorphic(groups2[0]) and groups1[1].isIsomorphic(groups2[2]) and groups1[2].isIsomorphic(groups2[1]): symmetryNumber *= 2
                    elif groups1[0].isIsomorphic(groups2[1]) and groups1[1].isIsomorphic(groups2[2]) and groups1[2].isIsomorphic(groups2[0]): symmetryNumber *= 2
                    elif groups1[0].isIsomorphic(groups2[1]) and groups1[1].isIsomorphic(groups2[0]) and groups1[2].isIsomorphic(groups2[2]): symmetryNumber *= 2
                    elif groups1[0].isIsomorphic(groups2[2]) and groups1[1].isIsomorphic(groups2[0]) and groups1[2].isIsomorphic(groups2[1]): symmetryNumber *= 2
                    elif groups1[0].isIsomorphic(groups2[2]) and groups1[1].isIsomorphic(groups2[1]) and groups1[2].isIsomorphic(groups2[0]): symmetryNumber *= 2
                
                
    return symmetryNumber

################################################################################

def calculateAxisSymmetryNumber(molecule):
    """
    Get the axis symmetry number correction. The "axis" refers to a series
    of two or more cumulated double bonds (e.g. C=C=C, etc.). Corrections
    for single C=C bonds are handled in getBondSymmetryNumber().
    
    Each axis (C=C=C) has the potential to double the symmetry number.
    If an end has 0 or 1 groups (eg. =C=CJJ or =C=C-R) then it cannot 
    alter the axis symmetry and is disregarded::
    
        A=C=C=C..        A-C=C=C=C-A
        
          s=1                s=1
    
    If an end has 2 groups that are different then it breaks the symmetry 
    and the symmetry for that axis is 1, no matter what's at the other end::
    
        A\               A\         /A
          T=C=C=C=C-A      T=C=C=C=T
        B/               A/         \B
              s=1             s=1
    
    If you have one or more ends with 2 groups, and neither end breaks the 
    symmetry, then you have an axis symmetry number of 2::
    
        A\         /B      A\         
          C=C=C=C=C          C=C=C=C-B
        A/         \B      A/         
              s=2                s=2
    """

    symmetryNumber = 1

    # List all double bonds in the structure
    doubleBonds = []
    for atom1 in molecule.vertices:
        for atom2 in atom1.edges:
            if atom1.edges[atom2].isDouble() and molecule.vertices.index(atom1) < molecule.vertices.index(atom2):
                doubleBonds.append((atom1, atom2))

    # Search for adjacent double bonds
    cumulatedBonds = []
    for i, bond1 in enumerate(doubleBonds):
        atom11, atom12 = bond1
        for bond2 in doubleBonds[i+1:]:
            atom21, atom22 = bond2
            if atom11 is atom21 or atom11 is atom22 or atom12 is atom21 or atom12 is atom22:
                listToAddTo = None
                for cumBonds in cumulatedBonds:
                    if (atom11, atom12) in cumBonds or (atom21, atom22) in cumBonds:
                        listToAddTo = cumBonds
                if listToAddTo is not None:
                    if (atom11, atom12) not in listToAddTo: listToAddTo.append((atom11, atom12))
                    if (atom21, atom22) not in listToAddTo: listToAddTo.append((atom21, atom22))
                else:
                    cumulatedBonds.append([(atom11, atom12), (atom21, atom22)])
    
    # Also keep isolated double bonds
    for bond1 in doubleBonds:
        for bonds in cumulatedBonds:
            if bond1 in bonds:
                break
        else:
            cumulatedBonds.append([bond1])
           
    # For each set of adjacent double bonds, check for axis symmetry
    for bonds in cumulatedBonds:
        
        # Do nothing if less than two cumulated bonds
        if len(bonds) < 1: continue

        # Do nothing if axis is in cycle
        found = False
        for atom1, atom2 in bonds:
           if molecule.isBondInCycle(atom1.edges[atom2]): found = True
        if found: continue

        # Find terminal atoms in axis
        # Terminal atoms labelled T:  T=C=C=C=T
        axis = []
        for bond in bonds: axis.extend(bond)
        terminalAtoms = []
        for atom in axis:
            if axis.count(atom) == 1: terminalAtoms.append(atom)
        if len(terminalAtoms) != 2: continue
        
        # Remove axis from (copy of) structure
        bondlist = []
        for atom1, atom2 in bonds:
            bond = atom1.edges[atom2]
            bondlist.append(bond)
            molecule.removeBond(bond)
        structure = molecule.copy(True)
        terminalAtoms = [structure.vertices[molecule.vertices.index(atom)] for atom in terminalAtoms]
        for bond in bondlist:
            molecule.addBond(bond)
        
        atomsToRemove = []
        for atom in structure.vertices:
            if len(atom.edges) == 0 and atom not in terminalAtoms: # it's not bonded to anything
                atomsToRemove.append(atom)
        for atom in atomsToRemove: structure.removeAtom(atom)

        # Split remaining fragments of structure
        end_fragments = structure.split()
        
        # 
        # there can be two groups at each end     A\         /B
        #                                           T=C=C=C=T
        #                                         A/         \B
        
        # to start with nothing has broken symmetry about the axis
        symmetry_broken=False
        end_fragments_to_remove = []
        for fragment in end_fragments: # a fragment is one end of the axis
            # remove the atom that was at the end of the axis and split what's left into groups
            terminalAtom = None
            for atom in terminalAtoms:
                if atom in fragment.atoms: 
                    terminalAtom = atom
                    fragment.removeAtom(atom)
                    break
            else:
                continue
            
            groups = []
            if len(fragment.atoms) > 0:
                groups = fragment.split()
            
            # If end has only one group then it can't contribute to (nor break) axial symmetry
            #   Eg. this has no axis symmetry:   A-T=C=C=C=T-A
            # so we remove this end from the list of interesting end fragments
            if len(groups) == 0:
                end_fragments_to_remove.append(fragment)
                continue # next end fragment
            elif len(groups)==1 and terminalAtom.radicalElectrons == 0:
                if terminalAtom.atomType.label == 'N3d':
                    symmetry_broken = True
                else:
                    end_fragments_to_remove.append(fragment)
                    continue # next end fragment
            elif len(groups)==1 and terminalAtom.radicalElectrons != 0:
                symmetry_broken = True
            elif len(groups)==2:
                if not groups[0].isIsomorphic(groups[1]):
                    # this end has broken the symmetry of the axis
                    symmetry_broken = True
        
        for fragment in end_fragments_to_remove:
            end_fragments.remove(fragment)
                            
        # If there are end fragments left that can contribute to symmetry,
        # and none of them broke it, then double the symmetry number
        # NB>> This assumes coordination number of 4 (eg. Carbon).
        #      And would be wrong if we had /B
        #                         =C=C=C=C=T-B
        #                                   \B
        #      (for some T with coordination number 5).
        if end_fragments and not symmetry_broken:
            symmetryNumber *= 2
                
    return symmetryNumber

################################################################################

def calculateCyclicSymmetryNumber(molecule):
    """
    Get the symmetry number correction for cyclic regions of a molecule.
    For complicated fused rings the smallest set of smallest rings is used.
    """
    from rdkit.Chem.rdmolops import SanitizeMol
    from rdkit.Chem.rdchem import Mol 

    symmetryNumber = 1
    
    molecule = molecule.copy(True)
    rings = molecule.getSmallestSetOfSmallestRings()
    try: 
        mcopy = molecule.toRDKitMol(removeHs=True, returnMapping=False)
        SanitizeMol(mcopy)
    except:
        # RDKit failed
        mcopy = None
        pass
    
    # Get symmetry number for each ring in structure
    for ring0 in rings:

        # Make another copy structure
        structure = molecule.copy(True)
        ring = [structure.atoms[molecule.atoms.index(atom)] for atom in ring0]
        
        # Figure out which atoms and bonds are aromatic and reassign appropriately:
        if mcopy is not None:
            for i, atom1 in enumerate(ring0):
                for atom2 in ring0[i+1:]:
                    if molecule.hasBond(atom1, atom2):
                        if mcopy.GetBondBetweenAtoms(i,i+1) is not None:
                            if str(mcopy.GetBondBetweenAtoms(i,i+1).GetBondType()) == 'AROMATIC':
                                bond = molecule.getBond(atom1, atom2)
                                bond.applyAction(['CHANGE_BOND', atom1, 'B', atom2])
                                atom1.atomType = atom2.atomType = rmgpy.molecule.atomTypes['Cb']
        
        # Remove bonds of ring from structure
        for i, atom1 in enumerate(ring):
            for atom2 in ring[i+1:]:
                if structure.hasBond(atom1, atom2):
                    structure.removeBond(atom1.edges[atom2])

        structures = structure.split()
        groups = []
        for struct in structures:
            for atom in ring:
                if struct.hasAtom(atom): struct.removeAtom(atom)
            groups.append(struct.split())
        # Find equivalent functional groups on ring
        equivalentGroups = []; equivalentGroupCount = []
        for group in groups:
            found = False
            for i, eqGroup in enumerate(equivalentGroups):
                if not found and len(group) == len(eqGroup):
                    for g, eg in zip(group, eqGroup):
                        if not g.isIsomorphic(eg):
                            # The groups do not match
                            break
                    else:
                        # The groups match
                        found = True
                if found:
                    # We've found a matching group, so increment its count
                    equivalentGroupCount[i] += 1        
                    break
            else:
                # No matching group found, so add it as a new group
                equivalentGroups.append(group)
                equivalentGroupCount.append(1)

        # Find equivalent bonds on ring
        equivalentBonds = []
        for i, atom1 in enumerate(ring0):
            for atom2 in ring0[i+1:]:
                if molecule.hasBond(atom1, atom2):
                    bond = molecule.getBond(atom1, atom2)
                    found = False
                    for eqBond in equivalentBonds:
                        if not found:
                            if bond.equivalent(eqBond[0]):
                                eqBond.append(group)
                                found = True
                    if not found:
                        equivalentBonds.append([bond])

        # Find maximum number of equivalent groups and bonds
        minEquivalentGroups = min(equivalentGroupCount)
        maxEquivalentGroups = max(equivalentGroupCount)
        minEquivalentBonds = None
        maxEquivalentBonds = 0
        for bonds in equivalentBonds:
            N = len(bonds)
            if minEquivalentBonds is None or N < minEquivalentBonds:
                minEquivalentBonds = N
            if N > maxEquivalentBonds:
                maxEquivalentBonds = N

        if maxEquivalentGroups == maxEquivalentBonds == len(ring):
            symmetryNumber *= len(ring) * 2
        else:
            symmetryNumber *= min(minEquivalentGroups, minEquivalentBonds)


    return symmetryNumber

################################################################################

def calculateSymmetryNumber(molecule):
    """
    Return the symmetry number for the structure. The symmetry number
    includes both external and internal modes.
    """
    symmetryNumber = 1

    for atom in molecule.vertices:
        if not molecule.isAtomInCycle(atom):
            symmetryNumber *= calculateAtomSymmetryNumber(molecule, atom)

    for atom1 in molecule.vertices:
        for atom2 in atom1.edges:
            if molecule.vertices.index(atom1) < molecule.vertices.index(atom2) and not molecule.isBondInCycle(atom1.edges[atom2]):
                symmetryNumber *= calculateBondSymmetryNumber(molecule, atom1, atom2)

    symmetryNumber *= calculateAxisSymmetryNumber(molecule)

    if molecule.isCyclic():
       symmetryNumber *= calculateCyclicSymmetryNumber(molecule)

    return symmetryNumber
