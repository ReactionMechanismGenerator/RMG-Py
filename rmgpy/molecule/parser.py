# global imports

import cython
import logging
import itertools

# local imports

# Assume that OB is not installed by default
INSTALLED_BACKENDS = {
    'OB': False,
}

try:
    import openbabel
    INSTALLED_BACKENDS['OB'] = True
except :
    pass

from rdkit import Chem

from rmgpy.molecule import element as elements

from rmgpy.molecule.util import retrieveElementCount, VALENCES, ORDERS
from rmgpy.molecule.inchi import AugmentedInChI, compose_aug_inchi_key, compose_aug_inchi, INCHI_PREFIX, MULT_PREFIX, U_LAYER_PREFIX

from .molecule import Atom, Bond, Molecule

# constants

BACKENDS = [
            'rdkit',
            ]

if INSTALLED_BACKENDS['OB']:
    BACKENDS.insert(0,'openbabel')

INCHI_LOOKUPS = {
            'H': '[H]',#RDkit was improperly handling the Hydrogen radical from InChI
            'He': '[He]',
        }
SMILES_LOOKUPS = {
            '[He]':# RDKit improperly handles helium and returns it in a triplet state
            """
            He
            multiplicity 1
            1 He u0 p1
            """
} 

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

# global variables:

try:
    SMILEwriter = openbabel.OBConversion()
    SMILEwriter.SetOutFormat('smi')
    SMILEwriter.SetOptions("i",SMILEwriter.OUTOPTIONS) # turn off isomer and stereochemistry information (the @ signs!)
except:
    pass

def reset_lone_pairs_to_default(at):
    """Resets the atom's lone pair count to its default value."""

    bondorder = 0
    bonds = at.edges.values()
    for bond in bonds:
        bondorder += ORDERS[bond.order]
    
    at.lonePairs = (VALENCES[at.element.symbol] - bondorder - at.radicalElectrons - at.charge) / 2

def convert_unsaturated_bond_to_biradical(mol, u_indices):
    """
    Convert an unsaturated bond (double, triple) into a bond
    with a lower bond order (single, double), and give an unpaired electron
    to each of the neighboring atoms, with indices referring to the 1-based
    index in the InChI string.
    """
    cython.declare(u1=cython.int, u2=cython.int)
    cython.declare(atom1=Atom, atom2=Atom)
    cython.declare(b=Bond)

    combos = itertools.combinations(u_indices, 2)

    for u1, u2 in combos:
        atom1 = mol.atoms[u1 - 1 ] # convert to 0-based index for atoms in molecule
        atom2 = mol.atoms[u2 - 1] # convert to 0-based index for atoms in molecule
        if mol.hasBond(atom1, atom2):
            b = mol.getBond(atom1, atom2)
            if not b.isSingle():
                atom1.radicalElectrons += 1
                atom2.radicalElectrons += 1
                
                if b.isTriple():
                    b.order = 'D'
                elif b.isDouble():
                    b.order = 'S'

                u_indices.remove(u1)
                u_indices.remove(u2)
                return mol        
    raise Exception('The indices {} did not refer to atoms that are connected in the molecule {}.'.format(u_indices, mol))    

def isUnsaturated(mol):
    """Does the molecule have a bond that's not single?
    
    (eg. a bond that is double or triple or beneze)"""
    cython.declare(atom1=Atom,
                   atom2=Atom,
                   bonds=dict,
                   bond=Bond)
    for atom1 in mol.atoms:
        bonds = mol.getBonds(atom1)
        for atom2, bond in bonds.iteritems():
            if not bond.isSingle():
                return True

    return False

def check_number_unpaired_electrons(mol):
    """Check if the number of unpaired electrons equals (m - 1)"""
    return mol.getNumberOfRadicalElectrons() == (mol.multiplicity - 1)        


def __fromSMILES(mol, smilesstr, backend):
    """Replace the Molecule `mol` with that given by the SMILES `smilesstr`
       using the backend `backend`"""
    if backend.lower() == 'rdkit':
        rdkitmol = Chem.MolFromSmiles(smilesstr)
        if rdkitmol is None:
            raise ValueError("Could not interpret the SMILES string {0!r}".format(smilesstr))
        fromRDKitMol(mol, rdkitmol)
        return mol
    elif backend.lower() == 'openbabel':
        parse_openbabel(mol, smilesstr, 'smi')
        return mol
    else:
        raise NotImplementedError('Unrecognized backend for SMILES parsing: {0}'.format(backend))

def __fromInChI(mol, inchistr, backend):
    """Replace the Molecule `mol` with that given by the InChI `inchistr`
       using the backend `backend`"""
    if backend.lower() == 'rdkit':
        rdkitmol = Chem.inchi.MolFromInchi(inchistr, removeHs=False)
        mol = fromRDKitMol(mol, rdkitmol)
        return mol 
    elif backend.lower() == 'openbabel':
        return parse_openbabel(mol, inchistr, 'inchi')
    else:
        raise NotImplementedError('Unrecognized backend for InChI parsing: {0}'.format(backend))


def __parse(mol, identifier, type_identifier, backend):
    """
    Parses the identifier based on the type of identifier (inchi/smi)
    and the backend used.
    
    First, look up the identifier in a dictionary to see if it can be processed
    this way.

    If not in the dictionary, parse it through the specified backed, 
    or try all backends.

    """

    if __lookup(mol, identifier, type_identifier) is not None:
        if isCorrectlyParsed(mol, identifier):
            return mol 

    for _backend in (BACKENDS if backend=='try-all' else [backend]):
        if type_identifier == 'smi':
            __fromSMILES(mol, identifier, _backend)
        elif type_identifier == 'inchi':
            __fromInChI(mol, identifier, _backend)
        else:
            raise NotImplementedError("Unknown identifier type {0}".format(type_identifier))

        if isCorrectlyParsed(mol, identifier):
            return mol
        else:
            logging.debug('Backend %s is not able to parse identifier %s', _backend, identifier)

    logging.error("Unable to correctly parse %s with backend %s", identifier, backend)
    raise Exception("Couldn't parse {0}".format(identifier))

def parse_openbabel(mol, identifier, type_identifier):
    """Converts the identifier to a Molecule using Openbabel."""
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(type_identifier, "smi")#SetInFormat(identifier) does not exist.
    obmol = openbabel.OBMol()
    obConversion.ReadString(obmol, identifier)
    obmol.AddHydrogens()
    obmol.AssignSpinMultiplicity(True)
    fromOBMol(mol, obmol)
    return mol


def isCorrectlyParsed(mol, identifier):
    """Check if molecule object has been correctly parsed."""
    conditions = []

    if mol.atoms:
        conditions.append(True)
    else:
        conditions.append(False)

    if 'InChI' in identifier:
        inchi_elementcount = retrieveElementCount(identifier)
        mol_elementcount = retrieveElementCount(mol)
        conditions.append(inchi_elementcount == mol_elementcount)

    return all(conditions)

def __lookup(mol, identifier, type_identifier):
    """
    Looks up the identifier and parses it the way we think is best.

    For troublesome inchis, we look up the smiles, and parse smiles.
    For troublesome smiles, we look up the adj list, and parse the adj list.

    """
    if type_identifier.lower() == 'inchi':
        try:
            smi = INCHI_LOOKUPS[identifier.split('/', 1)[1]]
            return mol.fromSMILES(smi)
        except KeyError:
            return None
    elif type_identifier.lower() == 'smi':
        try:
            adjList = SMILES_LOOKUPS[identifier]
            return mol.fromAdjacencyList(adjList)
        except KeyError:
            return None


def isZwitterIon(mol):
    """
    A zwitter ion is a molecule with total net charge 0, 
    but with some atoms bearing a formal charge.
    """

    total_charge = sum([at.charge for at in mol.atoms])
    abs_charge = sum([abs(at.charge) for at in mol.atoms])

    return (total_charge == 0) and (abs_charge != 0)
  
def check(mol, aug_inchi) :
    """Check if molecule corresponds to the aug. inchi"""
    cython.declare(conditions=list,
                   inchi=str,
                   multi=cython.int,
                   )
    conditions = []

    inchi, mult, u_indices = aug_inchi.inchi, aug_inchi.mult, aug_inchi.u_indices
    conditions.append(mult == mol.getRadicalCount() + 1)

    condition_electrons = True
    
    for at in mol.atoms:
        order = 0
        bonds = at.edges.values()
        for bond in bonds:
            order += ORDERS[bond.order]

        if (order + at.radicalElectrons + 2*at.lonePairs + at.charge) != VALENCES[at.symbol]:
            condition_electrons = False
            break

    conditions.append(condition_electrons)

    if not all(conditions):
        raise Exception('Molecule \n {0} does not correspond to aug. inchi {1}'.format(mol.toAdjacencyList(), aug_inchi))


def correct_O_triple_bond(mol):
    """
    Searches for a radical oxygen atom connected to 
    a closed-shell carbon via a triple bond.

    Converts the triple bond into a double bond,
    transfers the unpaired electron from O to C
    increases the lone pair count of O to 2.

    Only do this once per molecule.
    """

    for at in mol.atoms:
        if at.isOxygen() and at.radicalElectrons == 1 and at.lonePairs == 1:
            bonds = mol.getBonds(at)
            for atom2, bond in bonds.iteritems():
                if bond.isTriple():
                    bond.order = 'D'
                    at.radicalElectrons = 0
                    atom2.radicalElectrons = 1
                    at.lonePairs = 2
                    return

def fromInChI(mol, inchistr, backend='try-all'):
    """
    Convert an InChI string `inchistr` to a molecular structure. Uses 
    a user-specified backend for conversion, currently supporting
    rdkit (default) and openbabel.
    """

    mol.InChI = inchistr

    if INCHI_PREFIX in inchistr:
        return __parse(mol, inchistr, 'inchi', backend)
    else:
        return __parse(mol, INCHI_PREFIX + '/' + inchistr, 'inchi', backend)



def fromAugmentedInChI(mol, aug_inchi):
    """
    Creates a Molecule object from the augmented inchi.

    First, split off the multiplicity.
    Next, prepend the version layer to the inchi.
    Next, convert the inchi into a Molecule and
    set the multiplicity.

    Correct singlet one-center biradicals by replacing
    (u2, p0) by (u0, p1).

    Correct triplet two-center biradicals perceived as a
    zwitter ion by replacing (u0, c+/-1) by (u1, c0)

    Correct two-center triplet biradicals perceived as 
    a singlet double bond by replacing (u0)=(u0) by (u1)-(u1). 
    
    returns Molecule
    """

    if not isinstance(aug_inchi, AugmentedInChI):
        aug_inchi = AugmentedInChI(aug_inchi)

    mol = fromInChI(mol, aug_inchi.inchi)

    # multiplicity not specified in augmented InChI. Setting 
    if aug_inchi.mult == -1:
        logging.debug('Multiplicity not specified in augmented InChI.')
        logging.debug('Setting the multiplicity equal to the number of unpaired electrons + 1 of the parsed InChI.')
        mol.multiplicity = mol.getNumberOfRadicalElectrons() + 1
        return mol        

    mol.multiplicity = aug_inchi.mult


    #triplet to singlet conversion
    if mol.multiplicity == 1 and mol.getNumberOfRadicalElectrons() == 2:
        for at in mol.atoms:
            if at.radicalElectrons == 2:
                at.lonePairs = 1
                at.radicalElectrons = 0


    # zwitterion to triplet conversion
    if mol.multiplicity >= 3 and not check_number_unpaired_electrons(mol) and isZwitterIon(mol):
        for at in mol.atoms:
            if at.charge != 0 and at.radicalElectrons == 0:
                at.charge = 0
                at.radicalElectrons += 1
                if at.element.symbol == 'O':
                    bonds = mol.getBonds(at)
                    for atom2, bond in bonds.iteritems():
                        if bond.isDouble():
                            bond.order = 'S'
                            at.radicalElectrons = 0
                            atom2.radicalElectrons = 1
                            break

    # reset lone pairs                                
    for at in mol.atoms:
        reset_lone_pairs_to_default(at)

    # correct .O#C to O=C.
    correct_O_triple_bond(mol)


    # unsaturated bond to triplet conversion
    # sort the atoms so that they adopt the order of RMG's algorithm, which was used for the u indices:
    mol.sortAtoms()
    correct = check_number_unpaired_electrons(mol)

    unsaturated = isUnsaturated(mol)

    indices = aug_inchi.u_indices[:] if aug_inchi.u_indices is not None else None    
    
    if not correct and not indices:
        raise Exception('Cannot correct {} based on {} by converting unsaturated bonds into unpaired electrons...'\
            .format(mol.toAdjacencyList(), aug_inchi))

    while not correct and unsaturated and len(indices) > 1:
        mol = convert_unsaturated_bond_to_biradical(mol, indices)
        correct = check_number_unpaired_electrons(mol)
        unsaturated = isUnsaturated(mol)

    check(mol, aug_inchi)
    mol.sortAtoms()
    return mol

def fromSMILES(mol, smilesstr, backend='try-all'):
    """
    Convert a SMILES string `smilesstr` to a molecular structure. Uses 
    a user-specified backend for conversion, currently supporting
    rdkit (default) and openbabel.
    """
    return __parse(mol, smilesstr, 'smi', backend)


def fromSMARTS(mol, smartsstr):
    """
    Convert a SMARTS string `smartsstr` to a molecular structure. Uses
    `RDKit <http://rdkit.org/>`_ to perform the conversion.
    This Kekulizes everything, removing all aromatic atom types.
    """
    rdkitmol = Chem.MolFromSmarts(smartsstr)
    fromRDKitMol(mol, rdkitmol)
    return mol


def fromRDKitMol(mol, rdkitmol):
    """
    Convert a RDKit Mol object `rdkitmol` to a molecular structure. Uses
    `RDKit <http://rdkit.org/>`_ to perform the conversion.
    This Kekulizes everything, removing all aromatic atom types.
    """
    cython.declare(i=cython.int,
                   radicalElectrons=cython.int,
                   charge=cython.int,
                   lonePairs=cython.int,
                   number=cython.int,
                   order=cython.str,
                   atom=Atom,
                   atom1=Atom,
                   atom2=Atom,
                   bond=Bond)
    
    mol.vertices = []
    
    # Add hydrogen atoms to complete molecule if needed
    rdkitmol = Chem.AddHs(rdkitmol)
    Chem.rdmolops.Kekulize(rdkitmol, clearAromaticFlags=True)
    
    # iterate through atoms in rdkitmol
    for i in xrange(rdkitmol.GetNumAtoms()):
        rdkitatom = rdkitmol.GetAtomWithIdx(i)
        
        # Use atomic number as key for element
        number = rdkitatom.GetAtomicNum()
        element = elements.getElement(number)
            
        # Process charge
        charge = rdkitatom.GetFormalCharge()
        radicalElectrons = rdkitatom.GetNumRadicalElectrons()
        
        atom = Atom(element, radicalElectrons, charge, '', 0)
        mol.vertices.append(atom)
        
        # Add bonds by iterating again through atoms
        for j in xrange(0, i):
            rdkitatom2 = rdkitmol.GetAtomWithIdx(j + 1)
            rdkitbond = rdkitmol.GetBondBetweenAtoms(i, j)
            if rdkitbond is not None:
                order = ''
    
                # Process bond type
                rdbondtype = rdkitbond.GetBondType()
                if rdbondtype.name == 'SINGLE': order = 'S'
                elif rdbondtype.name == 'DOUBLE': order = 'D'
                elif rdbondtype.name == 'TRIPLE': order = 'T'
                elif rdbondtype.name == 'AROMATIC': order = 'B'
    
                bond = Bond(mol.vertices[i], mol.vertices[j], order)
                mol.addBond(bond)
    
    # Set atom types and connectivity values
    mol.update()
    mol.updateLonePairs()
    
    # Assume this is always true
    # There are cases where 2 radicalElectrons is a singlet, but
    # the triplet is often more stable, 
    mol.multiplicity = mol.getRadicalCount() + 1
    
    return mol

def fromOBMol(mol, obmol):
    """
    Convert a OpenBabel Mol object `obmol` to a molecular structure. Uses
    `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
    """
    # Below are the declared variables for cythonizing the module
    # cython.declare(i=cython.int)
    # cython.declare(radicalElectrons=cython.int, charge=cython.int, lonePairs=cython.int)
    # cython.declare(atom=Atom, atom1=Atom, atom2=Atom, bond=Bond)
    
    mol.vertices = []
    
    # Add hydrogen atoms to complete molecule if needed
    obmol.AddHydrogens()
    # TODO Chem.rdmolops.Kekulize(obmol, clearAromaticFlags=True)
    
    # iterate through atoms in obmol
    for obatom in openbabel.OBMolAtomIter(obmol):
        idx = obatom.GetIdx()#openbabel idx starts at 1!
        
        # Use atomic number as key for element
        number = obatom.GetAtomicNum()
        element = elements.getElement(number)
        # Process charge
        charge = obatom.GetFormalCharge()
        obatom_multiplicity = obatom.GetSpinMultiplicity()
        radicalElectrons =  obatom_multiplicity - 1 if obatom_multiplicity != 0 else 0
        
        atom = Atom(element, radicalElectrons, charge, '', 0)
        mol.vertices.append(atom)
    
    # iterate through bonds in obmol
    for obbond in openbabel.OBMolBondIter(obmol):
        order = 0
        # Process bond type
        oborder = obbond.GetBondOrder()
        if oborder == 1: order = 'S'
        elif oborder == 2: order = 'D'
        elif oborder == 3: order = 'T'
        elif obbond.IsAromatic() : order = 'B'

        bond = Bond(mol.vertices[obbond.GetBeginAtomIdx() - 1], mol.vertices[obbond.GetEndAtomIdx() - 1], order)#python array indices start at 0
        mol.addBond(bond)

    
    # Set atom types and connectivity values
    mol.update()
    mol.updateLonePairs()
    
    # Assume this is always true
    # There are cases where 2 radicalElectrons is a singlet, but
    # the triplet is often more stable, 
    mol.multiplicity = mol.getRadicalCount() + 1
    
    return mol


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
    # cython.declare(atom=Atom, atom1=Atom, bonds=dict, atom2=Atom, bond=Bond)
    # cython.declare(index1=cython.int, index2=cython.int, order=cython.int)

    # Sort the atoms before converting to ensure output is consistent
    # between different runs
    mol.sortAtoms()

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
    cython.declare(atoms=list,
                   rdAtomIndices=dict,
                   index=cython.int,
                   atom=Atom,
                   orders=dict,
                   atom1=Atom,
                   atom2=Atom,
                   index1=cython.int,
                   index2=cython.int,
                   )
                   
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

def createULayer(mol):
    """
    Creates a layer with the positions of the atoms that bear unpaired electrons.

    E.g.: u1,2 means that atoms with indices 1 and 2 bear one or more unpaired electrons

    Returns None if the molecule bears less than 2 unpaired electrons
    """
    ulayer = [str(i+1) for i, at in enumerate(mol.atoms) if at.radicalElectrons > 0]
    if ulayer:
        return (U_LAYER_PREFIX + ','.join(ulayer))
    else:
        return None


def toAugmentedInChI(mol):
    """
    Adds an extra layer to the InChI denoting the multiplicity
    of the molecule.
    
    Separate layer with a forward slash character.
    """
    inchi = toInChI(mol)

    mult = createMultiplicityLayer(mol.multiplicity)    
    ulayer = createULayer(mol)

    return compose_aug_inchi(inchi, mult, ulayer)

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
