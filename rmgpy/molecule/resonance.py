import cython

import rmgpy.molecule.generator as generator
import rmgpy.molecule.parser as parser

from .graph import Vertex, Edge, Graph, getVertexConnectivityValue
from .molecule import Atom, Bond, Molecule
import rmgpy.molecule.pathfinder as pathfinder

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
                isomer.updateAtomTypes()
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
                isomer.updateAtomTypes()
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
                isomer.updateAtomTypes()
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
                isomer.updateAtomTypes()
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
            molecule.updateAtomTypes()
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
    Generate the kekulized (single-double bond) form of the molecule.
    """
    cython.declare(isomers=list, atom=Atom)
    isomers = []
    for atom in mol.atoms:
        if atom.atomType.label == 'Cb' or atom.atomType.label == 'Cbf':
            break
    else:
        return isomers
   

    rdkitmol = generator.toRDKitMol(mol)  # This perceives aromaticit
    isomer = parser.fromRDKitMol(Molecule(), rdkitmol)# This step Kekulizes the molecule
    isomer.updateAtomTypes()
    isomers.append(isomer)  
    return isomers

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