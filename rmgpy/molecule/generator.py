# global imports

import logging
import itertools
from collections import Counter

# local imports
try:
    import openbabel
except:
    pass
from rdkit import Chem

from .inchi import compose_aug_inchi_key, compose_aug_inchi, parse_E_layer, parse_N_layer, INCHI_PREFIX, MULT_PREFIX, U_LAYER_PREFIX
from .molecule import Atom, Bond, Molecule
from .pathfinder import compute_atom_distance
from .util import partition, agglomerate, generate_combo

import rmgpy.molecule.resonance as resonance
# global variables:

#: This dictionary is used to shortcut lookups of a molecule's SMILES string from its chemical formula.
_known_smiles_molecules = {
                 'N2': 'N#N',
                 'CH4': 'C',
                 'H2O': 'O',
                 'C2H6': 'CC',
                 'H2': '[H][H]',
                 'H2O2': 'OO',
                 'C3H8': 'CCC',
                 'Ar': '[Ar]',
                 'He': '[He]',
                 'CH4O': 'CO',
                 'CO2': 'O=C=O',
                 'CO': '[C-]#[O+]',
                 'C2H4': 'C=C',
                 'O2': 'O=O'
             }

_known_smiles_radicals = {
                 'CH3': '[CH3]',
                 'HO': '[OH]',
                 'C2H5': 'C[CH2]',
                 'O': '[O]',
                 'HO2': '[O]O',
                 'CH': '[CH]',
                 'H': '[H]',
                 'C': '[C]',
                 #'CO2': it could be [O][C][O] or O=[C][O]
                 #'CO': '[C]=O', could also be [C][O]
                 #'C2H4': could  be [CH3][CH] or [CH2][CH2]
                 'O2': '[O][O]',
             }

def toInChI(mol):
    """
    Convert a molecular structure to an InChI string. Uses
    `RDKit <http://rdkit.org/>`_ to perform the conversion.
    Perceives aromaticity.
    
    or
    
    Convert a molecular structure to an InChI string. Uses
    `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
    """
    try:
        if not Chem.inchi.INCHI_AVAILABLE:
            return "RDKitInstalledWithoutInChI"
        rdkitmol = toRDKitMol(mol)
        return Chem.inchi.MolToInchi(rdkitmol, options='-SNon')
    except:
        pass

    obmol = toOBMol(mol)
    obConversion = openbabel.OBConversion()
    obConversion.SetOutFormat('inchi')
    obConversion.SetOptions('w', openbabel.OBConversion.OUTOPTIONS)
    return obConversion.WriteString(obmol).strip()

def create_U_layer(mol):
    """
    Creates a string with the positions of the atoms that bear unpaired electrons. The string
    can be used to complement the InChI with an additional layer that allows for the differentiation
    between structures with multiple unpaired electrons.

    The string is composed of a prefix ('u') followed by the positions of each of the unpaired electrons,
    sorted in numerical order.

    The indices in the string refer to the atom indices in the molecule, according to the atom order
    obtained by sorting the atoms using the InChI canonicalization algorithm.

    Example:
    - methyl radical ([CH3]) : u1
    - triplet methylene biradical ([CH2]) : u1,1
    - ethane-1,2-diyl biradical ([CH2][CH2]): u1,2
    
    First a deep copy is created of the original molecule and hydrogen atoms are removed from the molecule.
    Next, the molecule is converted into an InChI string, and the auxiliary information of the inchification 
    procedure is retrieved. 

    The N-layer is parsed and used to sort the atoms of the original order according
    to the order in the InChI. In case, the molecule contains atoms that cannot be distinguished
    with the InChI algorithm ('equivalent atoms'), the position of the unpaired electrons is changed
    as to ensure the atoms with the lowest indices are used to compose the string.

    When the molecule does not bear any unpaired electrons, None is returned.

    """

    if mol.getRadicalCount() == 0:
        return None
    elif mol.getFormula() == 'H':
        return U_LAYER_PREFIX + '1'

    molcopy = mol.copy(deep=True)

    hydrogens = [at for at in molcopy.atoms if at.number == 1]
    [molcopy.removeAtom(h) for h in hydrogens]

    rdkitmol = toRDKitMol(molcopy)
    _, auxinfo = Chem.MolToInchiAndAuxInfo(rdkitmol, options='-SNon')# suppress stereo warnings
    
    # extract the atom numbers from N-layer of auxiliary info:
    atom_indices = parse_N_layer(auxinfo)    
    atom_indices = [atom_indices.index(i + 1) for i, atom in enumerate(molcopy.atoms)]

    # sort the atoms based on the order of the atom indices
    molcopy.atoms = [x for (y,x) in sorted(zip(atom_indices, molcopy.atoms), key=lambda pair: pair[0])]

    # find the resonance isomer with the lowest u index:
    molcopy = normalize(molcopy)
    
    # create preliminary u-layer:
    u_layer = []
    for i, at in enumerate(molcopy.atoms):
        u_layer.extend([i+1] * at.radicalElectrons)
    
    # extract equivalent atom pairs from E-layer of auxiliary info:
    equivalent_atoms = parse_E_layer(auxinfo)
    if equivalent_atoms:
        # select lowest u-layer:
        u_layer = find_lowest_u_layer(molcopy, u_layer, equivalent_atoms)

    return (U_LAYER_PREFIX + ','.join(map(str, u_layer)))


def toAugmentedInChI(mol):
    """
    This function generates the augmented InChI canonical identifier, and that allows for the differentiation
    between structures with spin states and multiple unpaired electrons.

    Two additional layers are added to the InChI:
    - multiplicity layer: the total multiplicity of the molecule
    - unpaired electrons layer: the position of the unpaired electrons in the molecule

    """

    inchi = toInChI(mol)

    mult = createMultiplicityLayer(mol.multiplicity)    

    ulayer = create_U_layer(mol)

    aug_inchi = compose_aug_inchi(inchi, mult, ulayer)

    return aug_inchi

def toInChIKey(mol):
    """
    Convert a molecular structure to an InChI Key string. Uses
    `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
    
    or 
    
    Convert a molecular structure to an InChI Key string. Uses
    `RDKit <http://rdkit.org/>`_ to perform the conversion.
    
    Removes check-sum dash (-) and character so that only 
    the 14 + 9 characters remain.
    """
    try:
        if not Chem.inchi.INCHI_AVAILABLE:
            return "RDKitInstalledWithoutInChI"
        inchi = toInChI(mol)
        return Chem.inchi.InchiToInchiKey(inchi)[:-2]
    except:
        pass
    

#        for atom in mol.vertices:
#           if atom.isNitrogen():
    obmol = toOBMol(mol)
    obConversion = openbabel.OBConversion()
    obConversion.SetOutFormat('inchi')
    obConversion.SetOptions('w', openbabel.OBConversion.OUTOPTIONS)
    obConversion.SetOptions('K', openbabel.OBConversion.OUTOPTIONS)
    return obConversion.WriteString(obmol).strip()[:-2]

def toAugmentedInChIKey(mol):
    """
    Adds an extra layer to the InChIKey denoting the multiplicity
    of the molecule.

    Simply append the multiplicity string, do not separate by a
    character like forward slash.
    """
    key = toInChIKey(mol)
    
    mult_layer = '-mult'+str(mol.multiplicity) 
    ulayer = [str(i+1) for i, at in enumerate(mol.atoms) if at.radicalElectrons > 0]
    ulayer = '-u' + ','.join(ulayer) if mol.getNumberOfRadicalElectrons() > 1 else None

    return compose_aug_inchi_key(key, mult_layer, ulayer)

def createMultiplicityLayer(multiplicity):
    """
    Creates the string with the multiplicity information
    that is appended to the InChI to create an augmented InChI.
    
    """
    
    return MULT_PREFIX + str(multiplicity)


def toSMARTS(mol):
    """
    Convert a molecular structure to an SMARTS string. Uses
    `RDKit <http://rdkit.org/>`_ to perform the conversion.
    Perceives aromaticity and removes Hydrogen atoms.
    """
    rdkitmol = toRDKitMol(mol)
    
    return Chem.MolToSmarts(rdkitmol)


def toSMILES(mol):
    """
    Convert a molecular structure to an SMILES string. 
    
    If there is a Nitrogen atom present it uses
    `OpenBabel <http://openbabel.org/>`_ to perform the conversion,
    and the SMILES may or may not be canonical.
    
    Otherwise, it uses `RDKit <http://rdkit.org/>`_ to perform the 
    conversion, so it will be canonical SMILES.
    While converting to an RDMolecule it will perceive aromaticity
    and removes Hydrogen atoms.
    """
    
    # If we're going to have to check the formula anyway,
    # we may as well shortcut a few small known molecules.
    # Dictionary lookups are O(1) so this should be fast:
    # The dictionary is defined at the top of this file.
    try:
        if mol.isRadical():
            return _known_smiles_radicals[mol.getFormula()]
        else:
            return _known_smiles_molecules[mol.getFormula()]
    except KeyError:
        # It wasn't in the above list.
        pass
    for atom in mol.vertices:
        if atom.isNitrogen():
            obmol = toOBMol(mol)
            try:
                SMILEwriter = openbabel.OBConversion()
                SMILEwriter.SetOutFormat('smi')
                SMILEwriter.SetOptions("i",SMILEwriter.OUTOPTIONS) # turn off isomer and stereochemistry information (the @ signs!)
            except:
                pass
            return SMILEwriter.WriteString(obmol).strip()

    rdkitmol = toRDKitMol(mol, sanitize=False)
    if not mol.isAromatic():
        return Chem.MolToSmiles(rdkitmol, kekuleSmiles=True)
    return Chem.MolToSmiles(rdkitmol)

def toOBMol(mol):
    """
    Convert a molecular structure to an OpenBabel OBMol object. Uses
    `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
    """

    atoms = mol.vertices

    obmol = openbabel.OBMol()
    for atom in atoms:
        a = obmol.NewAtom()
        a.SetAtomicNum(atom.number)
        a.SetFormalCharge(atom.charge)
    orders = {'S': 1, 'D': 2, 'T': 3, 'B': 5}
    for atom1 in mol.vertices:
        for atom2, bond in atom1.edges.iteritems():
            index1 = atoms.index(atom1)
            index2 = atoms.index(atom2)
            if index1 < index2:
                order = orders[bond.order]
                obmol.AddBond(index1+1, index2+1, order)

    obmol.AssignSpinMultiplicity(True)

    return obmol

def toRDKitMol(mol, removeHs=True, returnMapping=False, sanitize=True):
    """
    Convert a molecular structure to a RDKit rdmol object. Uses
    `RDKit <http://rdkit.org/>`_ to perform the conversion.
    Perceives aromaticity and, unless removeHs==False, removes Hydrogen atoms.
    
    If returnMapping==True then it also returns a dictionary mapping the 
    atoms to RDKit's atom indices.
    """
    
    # Sort the atoms before converting to ensure output is consistent
    # between different runs
    mol.sortAtoms()           
    atoms = mol.vertices
    rdAtomIndices = {} # dictionary of RDKit atom indices
    rdkitmol = Chem.rdchem.EditableMol(Chem.rdchem.Mol())
    for index, atom in enumerate(mol.vertices):
        rdAtom = Chem.rdchem.Atom(atom.element.symbol)
        rdAtom.SetNumRadicalElectrons(atom.radicalElectrons)
        if atom.element.symbol == 'C' and atom.lonePairs == 1 and mol.multiplicity == 1: rdAtom.SetNumRadicalElectrons(2)
        rdkitmol.AddAtom(rdAtom)
        if removeHs and atom.symbol == 'H':
            pass
        else:
            rdAtomIndices[atom] = index
    
    rdBonds = Chem.rdchem.BondType
    orders = {'S': rdBonds.SINGLE, 'D': rdBonds.DOUBLE, 'T': rdBonds.TRIPLE, 'B': rdBonds.AROMATIC}
    # Add the bonds
    for atom1 in mol.vertices:
        for atom2, bond in atom1.edges.iteritems():
            index1 = atoms.index(atom1)
            index2 = atoms.index(atom2)
            if index1 < index2:
                order = orders[bond.order]
                rdkitmol.AddBond(index1, index2, order)
    
    # Make editable mol into a mol and rectify the molecule
    rdkitmol = rdkitmol.GetMol()
    if sanitize:
        Chem.SanitizeMol(rdkitmol)
    if removeHs:
        rdkitmol = Chem.RemoveHs(rdkitmol, sanitize=sanitize)
    
    if returnMapping:
        return rdkitmol, rdAtomIndices
    return rdkitmol

def is_valid_combo(combo, mol, distances):
    """
    Check if the combination of atom indices refers to
    atoms that are adjacent in the molecule.
    """

    # compute shortest path between atoms 
    agglomerates = agglomerate(combo)
    new_distances = compute_agglomerate_distance(agglomerates, mol)

    # combo is valid if the distance is equal to the parameter distance

    if len(distances) != len(new_distances): return False

    for orig_dist, new_dist in zip(distances, new_distances):
        # only compare the values of the dictionaries:
        if sorted(orig_dist.values()) != sorted(new_dist.values()):
            return False

    return True

def find_lowest_u_layer(mol, u_layer, equivalent_atoms):
    """..."""

    if not equivalent_atoms:
        return u_layer

    new_u_layer = []

    grouped_electrons, corresponding_E_layers = partition(u_layer, equivalent_atoms)

    # don't process atoms that do not belong to an equivalence layer
    grouped_electrons_copy = grouped_electrons[:]
    corresponding_E_layers_copy = corresponding_E_layers[:]
    for group, e_layer in zip(grouped_electrons_copy, corresponding_E_layers_copy):
        if not e_layer:
            new_u_layer.extend(group)
            grouped_electrons.remove(group)
            corresponding_E_layers.remove(e_layer)


    combos = generate_combo(grouped_electrons, corresponding_E_layers)
    # compute original distance:
    orig_agglomerates = agglomerate(grouped_electrons)
    orig_distances = compute_agglomerate_distance(orig_agglomerates, mol)

    # deflate the list of lists to be able to numerically compare them
    selected_group = sorted(itertools.chain.from_iterable(grouped_electrons))

    # see if any of the combos is valid and results in a lower numerical combination than the original 
    for combo in combos:    
        if is_valid_combo(combo, mol, orig_distances):
            combo = sorted(itertools.chain.from_iterable(combo))
            if combo < selected_group:
                selected_group = combo

    # add the minimized unpaired electron positions to the u-layer:
    new_u_layer.extend(selected_group)

    return sorted(new_u_layer)

def normalize(mol):
    """
    Select the resonance isomer that is isomorphic to the parameter isomer, with the lowest unpaired
    electrons descriptor.

    We generate over all resonance isomers (non-isomorphic as well as isomorphic) and add isomorphic
    isomers to a list (candidates).

    Next, we search through this list and return the candidate with the lowest unpaired electrons descriptor.

    """
    candidates = resonance.generate_isomorphic_isomers(mol)
    
    current_minimum = candidates[0]#isomer with the smallest u-layer descriptor.
    unpaired_electrons_current_minimum = get_unpaired_electrons(current_minimum)
    for cand in candidates[1:]:
       unpaired_electrons_candidate = get_unpaired_electrons(cand)
       if unpaired_electrons_candidate < unpaired_electrons_current_minimum:
            current_minimum = cand
            unpaired_electrons_current_minimum = unpaired_electrons_candidate


    return current_minimum


def get_unpaired_electrons(mol):
    """
    returns a sorted list of the indices of the atoms that bear one or more 
    unpaired electrons.
    """
    locations = []
    for index, at in enumerate(mol.atoms):
        if at.radicalElectrons >= 1:
            locations.append(index)

    return sorted(locations)  

def compute_agglomerate_distance(agglomerates, mol):
    distances = []
    for agglomerate in agglomerates:
        dist = compute_atom_distance(agglomerate, mol)
        distances.append(dist)

    return distances

