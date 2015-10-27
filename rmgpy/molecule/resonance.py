import cython

import rmgpy.molecule.generator as generator
import rmgpy.molecule.parser as parser

from .graph import Vertex, Edge, Graph, getVertexConnectivityValue
from .molecule import Atom, Bond, Molecule

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
            
        newIsomers = getAdjacentResonanceIsomers(isomer)
        newIsomers += getLonePairRadicalResonanceIsomers(isomer)
        newIsomers += getN5dd_N5tsResonanceIsomers(isomer)
        newIsomers += getKekulizedResonanceIsomers(isomer)
        
        for newIsomer in newIsomers:
            newIsomer.updateAtomTypes()
            # Append to isomer list if unique
            for isom in isomers:
                if isom.isIsomorphic(newIsomer):
                    break
            else:
                isomers.append(newIsomer)
        
        newIsomers = getAromaticResonanceIsomers(isomer)
        # Perform extra check for aromatic isomers when updating atomtypes
        for newIsomer in newIsomers:
            try:
                newIsomer.updateAtomTypes()
            except:
                # Something incorrect has happened, ie. 2 double bonds on a Cb atomtype
                # Do not add the new isomer since it is malformed
                continue 
            # Append to isomer list if unique
            for isom in isomers:
                if isom.isIsomorphic(newIsomer):
                    break
            else:
                isomers.append(newIsomer)
        
                    
        # Move to next resonance isomer
        index += 1
    
    return isomers

def getAdjacentResonanceIsomers(mol):
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
            paths = findAllDelocalizationPaths(mol, atom)
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
                    v2.connectivity = v1.connectivity
                    v2.sortingLabel = v1.sortingLabel
                # Restore current isomer
                atom1.incrementRadical()
                atom3.decrementRadical()
                bond12.decrementOrder()
                bond23.incrementOrder()
                # Append to isomer list if unique
                isomers.append(isomer)

    return isomers

def getLonePairRadicalResonanceIsomers(mol):
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
            paths = findAllDelocalizationPathsLonePairRadical(mol, atom)
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
                    v2.connectivity = v1.connectivity
                    v2.sortingLabel = v1.sortingLabel
                # Restore current isomer
                atom1.incrementRadical()
                atom1.decrementLonePairs()
                atom1.updateCharge()
                atom2.decrementRadical()
                atom2.incrementLonePairs()
                atom2.updateCharge()
                # Append to isomer list if unique
                isomers.append(isomer)

    return isomers

def getN5dd_N5tsResonanceIsomers(mol):
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
        paths = findAllDelocalizationPathsN5dd_N5ts(mol, atom)
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
                    v2.connectivity = v1.connectivity
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
                    v2.connectivity = v1.connectivity
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
                isomers.append(isomer)
                
    return isomers

def getAromaticResonanceIsomers(mol):
    """
    Generate the aromatic form of the molecule.
    """
    cython.declare(isomers=list, molecule=Molecule, rdAtomIndices=dict, aromatic=cython.bint, aromaticBonds=list)
    cython.declare(rings=list, ring0=list, i=cython.int, atom1=Atom, atom2=Atom, bond=Bond)
    
    isomers = []

    # Radicals
    if mol.isCyclic():
        molecule = mol.copy(deep=True)
        try:
            rdkitmol, rdAtomIndices = generator.toRDKitMol(molecule, removeHs=False, returnMapping=True)
        except:
            return []
        aromatic = False
        rings = molecule.getSmallestSetOfSmallestRings()            
        for ring0 in rings:
            # In RMG, only 6-member rings can be considered aromatic, so ignore all other rings                
            aromaticBonds = []
            if len(ring0) == 6:
                # Figure out which atoms and bonds are aromatic and reassign appropriately:
                for i, atom1 in enumerate(ring0):
                    if not atom1.isCarbon():
                        # all atoms in the ring must be carbon in RMG for our definition of aromatic
                        break
                    for atom2 in ring0[i+1:]:
                        if molecule.hasBond(atom1, atom2):
                            if str(rdkitmol.GetBondBetweenAtoms(rdAtomIndices[atom1],rdAtomIndices[atom2]).GetBondType()) == 'AROMATIC':
                                aromaticBonds.append(molecule.getBond(atom1, atom2))
            if len(aromaticBonds) == 6:
                aromatic = True
                # Only change bonds if there are all 6 are aromatic.  Otherwise don't do anything
                for bond in aromaticBonds:
                    bond.order = 'B'
                    
        if aromatic:              
            isomers.append(molecule)

    return isomers

def getKekulizedResonanceIsomers(mol):
    """
    Generate the kekulized (single-double bond) form of the molecule.
    """
    cython.declare(isomers=list, atom=Atom)
    isomers = []
    for atom in mol.vertices:
        if atom.atomType.label == 'Cb' or atom.atomType.label == 'Cbf':
            break
    else:
        return isomers
    
    rdkitmol = generator.toRDKitMol(mol)  # This perceives aromaticity
    mol = Molecule()
    isomers.append(parser.fromRDKitMol(mol, rdkitmol))  # This step Kekulizes the molecule
    return isomers

def findAllDelocalizationPaths(mol, atom1):
    """
    Find all the delocalization paths allyl to the radical center indicated
    by `atom1`. Used to generate resonance isomers.
    """
    cython.declare(paths=list)
    cython.declare(atom2=Atom, atom3=Atom, bond12=Bond, bond23=Bond)
    
    # No paths if atom1 is not a radical
    if atom1.radicalElectrons <= 0:
        return []

    # Find all delocalization paths
    paths = []
    for atom2, bond12 in atom1.edges.items():
        # Vinyl bond must be capable of gaining an order
        if (bond12.isSingle() or bond12.isDouble()) and atom1.radicalElectrons == 1:
            for atom3, bond23 in atom2.edges.items():
                # Allyl bond must be capable of losing an order without breaking
                if atom1 is not atom3 and (bond23.isDouble() or bond23.isTriple()):
                    paths.append([atom1, atom2, atom3, bond12, bond23])
    return paths

def findAllDelocalizationPathsLonePairRadical(mol, atom1):
    """
    Find all the delocalization paths of lone electron pairs next to the radical center indicated
    by `atom1`. Used to generate resonance isomers.
    """
    cython.declare(paths=list)
    cython.declare(atom2=Atom, bond12=Bond)
    
    # No paths if atom1 is not a radical
    if atom1.radicalElectrons <= 0:
        return []
    
    # In a first step we only consider nitrogen and oxygen atoms as possible radical centers
    if not ((atom1.lonePairs == 0 and atom1.isNitrogen()) or(atom1.lonePairs == 2 and atom1.isOxygen())):
        return []
    
    # Find all delocalization paths
    paths = []
    for atom2, bond12 in atom1.edges.items():
        # Only single bonds are considered
        if bond12.isSingle():
            # Neighboring atom must posses a lone electron pair to loose it
            if ((atom2.lonePairs == 1 and atom2.isNitrogen()) or (atom2.lonePairs == 3 and atom2.isOxygen())) and (atom2.radicalElectrons == 0):
                paths.append([atom1, atom2])
                
    return paths

def findAllDelocalizationPathsN5dd_N5ts(mol, atom1):
    """
    Find all the resonance structures of nitrogen atoms with two double bonds (N5dd)
    and nitrogen atoms with one triple and one single bond (N5ts)
    """
    cython.declare(paths=list)
    cython.declare(atom2=Atom, bond12=Bond)
    
    # No paths if atom1 is not nitrogen
    if not (atom1.isNitrogen()):
        return []
    
    # Find all delocalization paths
    paths = []
    index_atom_2 = 0
    index_atom_3 = 0
    
    for atom2, bond12 in atom1.edges.items():
        index_atom_2 = index_atom_2 + 1
        # Only double bonds are considered
        if bond12.isDouble():
            for atom3, bond13 in atom1.edges.items():
                index_atom_3 = index_atom_3 + 1
                # Only double bonds are considered, at the moment we only consider non-radical nitrogen and oxygen atoms
                if (bond13.isDouble() and atom3.radicalElectrons == 0 and atom3.lonePairs > 0 and not atom3.isOxygen() and not atom3.isCarbon() and (index_atom_2 != index_atom_3)):
                    paths.append([atom1, atom2, atom3, bond12, bond13, 1])
    
    for atom2, bond12 in atom1.edges.items():
        # Only triple bonds are considered
        if bond12.isTriple():
            for atom3, bond13 in atom1.edges.items():
                # Only single bonds are considered, at the moment we only consider negatively charged nitrogen and oxygen
                if (bond13.isSingle() and ((atom3.isNitrogen() and atom3.lonePairs >= 2) or (atom3.isOxygen() and atom3.lonePairs >= 3))):
                    paths.append([atom1, atom2, atom3, bond12, bond13, 2])
    
    return paths
