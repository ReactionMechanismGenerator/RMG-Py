# global imports

import cython
import logging
import itertools
from collections import Counter
import re

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
from rmgpy.molecule.inchi import AugmentedInChI, compose_aug_inchi_key, compose_aug_inchi, parse_H_layer, INCHI_PREFIX, MULT_PREFIX, U_LAYER_PREFIX

from rmgpy.molecule.pathfinder import \
 find_butadiene,\
 find_butadiene_end_with_charge,\
 find_allyl_end_with_charge

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

def reset_lone_pairs_to_default(at):
    """Resets the atom's lone pair count to its default value."""

    bondorder = 0
    bonds = at.edges.values()
    for bond in bonds:
        bondorder += ORDERS[bond.order]
    
    at.lonePairs = (VALENCES[at.element.symbol] - bondorder - at.radicalElectrons - at.charge) / 2

def convert_unsaturated_bond_to_biradical(mol, inchi, u_indices):
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
        atom1 = mol.atoms[u1 - 1] # convert to 0-based index for atoms in molecule
        atom2 = mol.atoms[u2 - 1] # convert to 0-based index for atoms in molecule
        if mol.hasBond(atom1, atom2):
            b = mol.getBond(atom1, atom2)
            if not b.isSingle():
                atom1.radicalElectrons += 1
                atom2.radicalElectrons += 1
                b.decrementOrder()

                u_indices.remove(u1)
                u_indices.remove(u2)

                return mol
            else:#maybe it's a mobile H-layer problem
                all_mobile_h_atoms_couples = parse_H_layer(inchi)
                """
                assume that O2=C1-O3H is the keto-enol system
                    1) find its partner (O2)
                    2) transfer H atom to partner (O2)
                    3) change bond order between partner and central carbon
                    4) add unpaired electrons to central carbon and original O.

                """
                if all_mobile_h_atoms_couples:
                    #find central atom:
                    central, original_atom, new_partner = find_mobile_h_system(mol, 
                        all_mobile_h_atoms_couples, [u1, u2])

                    # search Hydrogen atom and bond
                    hydrogen = None
                    for at, bond in original_atom.bonds.iteritems():
                        if at.number == 1:
                            hydrogen = at
                            mol.removeBond(bond)
                            break

                    new_h_bond = Bond(new_partner, hydrogen, order='S')
                    mol.addBond(new_h_bond)
                    
                    mol.getBond(central, new_partner).decrementOrder()

                    central.radicalElectrons += 1
                    original_atom.radicalElectrons += 1

                    u_indices.remove(u1)
                    u_indices.remove(u2)
                    return mol
        else:
            path = find_butadiene(atom1, atom2)
            if path is not None:
                atom1.radicalElectrons += 1
                atom2.radicalElectrons += 1
                # filter bonds from path and convert bond orders:
                bonds = path[1::2]#odd elements
                for bond in bonds[::2]:# even bonds
                    assert isinstance(bond, Bond)
                    bond.decrementOrder()
                for bond in bonds[1::2]:# odd bonds
                    assert isinstance(bond, Bond)
                    bond.incrementOrder()    

                u_indices.remove(u1)
                u_indices.remove(u2)

                return mol

    raise Exception('The indices {} did not refer to atoms that are connected in the molecule {}.'
        .format(u_indices, mol.toAdjacencyList()))    

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


def contains_charge(mol):
    abs_charge = sum([abs(at.charge) for at in mol.atoms])

    return abs_charge != 0
  
def check(mol, aug_inchi) :
    """Check if molecule corresponds to the aug. inchi"""
    cython.declare(conditions=list,
                   inchi=str,
                   multi=cython.int,
                   )

    _, mult, __ = aug_inchi.inchi, aug_inchi.mult, aug_inchi.u_indices
    assert mult == mol.getRadicalCount() + 1,\
     'Multiplicity of molecule \n {0} does not correspond to aug. inchi {1}'.format(mol.toAdjacencyList(), aug_inchi)
    
    for at in mol.atoms:
        order = 0
        bonds = at.edges.values()
        for bond in bonds:
            order += ORDERS[bond.order]

        assert (order + at.radicalElectrons + 2*at.lonePairs + at.charge) == VALENCES[at.symbol],\
            'Valency for an atom of molecule \n {0} does not correspond to aug. inchi {1}'.format(mol.toAdjacencyList(), aug_inchi)

def correct_O_unsaturated_bond(mol, u_indices):
    """
    Searches for a radical or a charged oxygen atom connected to 
    a closed-shell carbon via an unsatured bond.

    Decrements the unsatured bond,
    transfers the unpaired electron from O to C or
    converts the charge from O to an unpaired electron on C, 
    increases the lone pair count of O to 2.

    Only do this once per molecule.
    """

    for at in mol.atoms:
        if at.isOxygen() and at.radicalElectrons == 1 and at.lonePairs == 1:
            bonds = mol.getBonds(at)
            oxygen = at
            for atom2, bond in bonds.iteritems():
                if bond.isTriple():
                    bond.decrementOrder()
                    oxygen.radicalElectrons -= 1
                    atom2.radicalElectrons += 1
                    oxygen.lonePairs += 1
                    return
        elif at.isOxygen() and at.charge == 1 and at.lonePairs == 1:
            bonds = mol.getBonds(at)
            oxygen = at

            start = oxygen
            # search for 3-atom-2-bond [X=X-X] paths
            paths = find_allyl_end_with_charge(start)
            for path in paths:    
                end = path[-1]
                start.charge += 1 if start.charge < 0 else -1
                end.charge += 1 if end.charge < 0 else -1
                start.lonePairs += 1
                # filter bonds from path and convert bond orders:
                bonds = path[1::2]#odd elements
                for bond in bonds[::2]:# even bonds
                    assert isinstance(bond, Bond)
                    bond.decrementOrder()
                for bond in bonds[1::2]:# odd bonds
                    assert isinstance(bond, Bond)
                    bond.incrementOrder()  
                return
            else:
                for atom2, bond in bonds.iteritems():
                    if not bond.isSingle() and atom2.charge == 0:
                        oxygen.charge -= 1
                        if (mol.atoms.index(atom2) + 1) in u_indices:
                            bond.decrementOrder()
                            atom2.radicalElectrons += 1
                            u_indices.remove(mol.atoms.index(atom2) + 1)
                        oxygen.lonePairs += 1
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

    indices = aug_inchi.u_indices[:] if aug_inchi.u_indices is not None else []
    c = Counter(indices)
    for k,v in c.iteritems():
        atom = mol.atoms[k - 1]
        [indices.remove(k) for _ in range(atom.radicalElectrons)]


    if mol.multiplicity >= 3 and not check_number_unpaired_electrons(mol) and contains_charge(mol):
        fixCharge(mol, indices)

    # reset lone pairs                                
    for at in mol.atoms:
        reset_lone_pairs_to_default(at)

    # correct .O#C to O=C.
    correct_O_unsaturated_bond(mol, indices)


    # unsaturated bond to triplet conversion
    correct = check_number_unpaired_electrons(mol)

    unsaturated = isUnsaturated(mol)
    
    if not correct and not indices:
        raise Exception( 'Cannot correct {} based on {} by converting unsaturated bonds into unpaired electrons...'\
            .format(mol.toAdjacencyList(), aug_inchi))

    while not correct and unsaturated and len(indices) > 1:
        mol = convert_unsaturated_bond_to_biradical(mol, aug_inchi.inchi, indices)
        correct = check_number_unpaired_electrons(mol)
        unsaturated = isUnsaturated(mol)

    check(mol, aug_inchi)
    mol.updateAtomTypes()
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
    mol.updateConnectivityValues()
    mol.updateAtomTypes()
    mol.updateMultiplicity()
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
            try:
                obmol = toOBMol(mol)
                SMILEwriter = openbabel.OBConversion()
                SMILEwriter.SetOutFormat('smi')
                SMILEwriter.SetOptions("i",SMILEwriter.OUTOPTIONS) # turn off isomer and stereochemistry information (the @ signs!)
                return SMILEwriter.WriteString(obmol).strip()
            except:
                break# break from for loop

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

    if INSTALLED_BACKENDS['OB']: 
        obmol = toOBMol(mol)
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat('inchi')
        obConversion.SetOptions('w', openbabel.OBConversion.OUTOPTIONS)
        return obConversion.WriteString(obmol).strip()
    else:
        raise Exception('Could not generate InChI, because Openbabel installation was not found. ')

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
    

    if INSTALLED_BACKENDS['OB']:         
        obmol = toOBMol(mol)
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat('inchi')
        obConversion.SetOptions('w', openbabel.OBConversion.OUTOPTIONS)
        obConversion.SetOptions('K', openbabel.OBConversion.OUTOPTIONS)
        return obConversion.WriteString(obmol).strip()[:-2]
    else:
        raise Exception('Could not generate InChI Key, because Openbabel installation was not found. ')
        
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

def fixCharge(mol, u_indices):
    """
    Fix molecules perceived as zwitterions that in reality are structures
    with multiple unpaired electrons.

    The simplest case converts atoms with a charge to atoms with one more
    unpaired electron.

    """
    if not u_indices:
        return

    # converting charges to unpaired electrons for atoms in the u-layer
    for at in mol.atoms:
        if at.charge != 0 and (mol.atoms.index(at) + 1) in u_indices:
            at.charge += 1 if at.charge < 0 else -1
            at.radicalElectrons += 1
            u_indices.remove(mol.atoms.index(at) + 1)

    # convert neighboring atoms (or delocalized paths) to unpaired electrons
    u_indices_copy = u_indices[:]
    for index in u_indices_copy:
        start = mol.atoms[index -1]

        # search for 4-atom-3-bond [X=X-X=X] paths
        path = find_butadiene_end_with_charge(start)
        if path is not None:    
            # we have found the atom we are looking for
            start.radicalElectrons += 1
            end = path[-1]
            end.charge += 1 if end.charge < 0 else -1
            end.lonePairs += 1
            # filter bonds from path and convert bond orders:
            bonds = path[1::2]#odd elements
            for bond in bonds[::2]:# even bonds
                assert isinstance(bond, Bond)
                bond.decrementOrder()
            for bond in bonds[1::2]:# odd bonds
                assert isinstance(bond, Bond)
                bond.incrementOrder()  
            u_indices.remove(mol.atoms.index(start) + 1)
            continue

        # search for 3-atom-2-bond [X=X-X] paths
        paths = find_allyl_end_with_charge(start)
        from rmgpy.data.kinetics.family import ReactionRecipe

        for path in paths:
            # label atoms so that we can use the labels in the actions of the recipe
            for i, at in enumerate(path[::2]):
                assert isinstance(at, Atom)
                at.label = str(i)
            # we have found the atom we are looking for
            fix_charge_recipe = ReactionRecipe()
            fix_charge_recipe.addAction(['GAIN_RADICAL', start.label, 1])

            end = path[-1]
            end_original_charge = end.charge
          
            # filter bonds from path and convert bond orders:
            bonds = path[1::2]#odd elements
            for bond in bonds[::2]:# even bonds
                assert isinstance(bond, Bond)
                fix_charge_recipe.addAction(['CHANGE_BOND', bond.atom1.label, -1, bond.atom2.label])
            for bond in bonds[1::2]:# odd bonds
                assert isinstance(bond, Bond)
                fix_charge_recipe.addAction(['CHANGE_BOND', bond.atom1.label, 1, bond.atom2.label])

            end.charge += 1 if end.charge < 0 else -1
            fix_charge_recipe.applyForward(mol, update=False)

            if check_bond_order_oxygen(mol):
                u_indices.remove(mol.atoms.index(start) + 1)
                # unlabel atoms so that they never cause trouble downstream
                for i, at in enumerate(path[::2]):
                    assert isinstance(at, Atom)
                    at.label = ''
                break
            else:
                fix_charge_recipe.applyReverse(mol, update=False)
                end.charge = end_original_charge

                # unlabel atoms so that they never cause trouble downstream
                for i, at in enumerate(path[::2]):
                    assert isinstance(at, Atom)
                    at.label = ''

                continue # to next path

            
        continue # to next index in u-layer

    # fix adjacent charges
    for at in mol.atoms:
        if at.charge != 0:
            for neigh, bond in at.bonds.iteritems():
                if neigh.charge != 0:
                    bond.incrementOrder()
                    at.charge += 1 if at.charge < 0 else -1
                    neigh.charge += 1 if neigh.charge < 0 else -1


def moveHs(mol):
    """
    Sorts the list of atoms so that heavy atoms always come before H atoms.
    """

    atoms = mol.atoms
    heavy_atom_count = sum([1 for at in atoms if at.number != 1])
    i = 0
    while True:# continue until the atoms are ordered
        # search for a (H,nonH) pair and switch:

        # search for H-atom in 0 < i < heavy atom count:
        while i < heavy_atom_count:
            if atoms[i].number != 1:
                i += 1 # increase and go to next
                continue # ignore heavy atoms
            h_index = i

            # search for non-H atom in i < j < atom count:
            j = i + 1
            while j < len(atoms):
                if atoms[j].number == 1:
                    h_index = j
                    j += 1
                else:# we have found a (H, nonH) pair
                    non_h_index = j
                    atoms[h_index], atoms[non_h_index] = atoms[non_h_index], atoms[h_index]
                    break

            break

        if i == heavy_atom_count:#all atoms are ordered.
            break

    mol.atoms = atoms

def updateAtomConnectivityValues(mol):
    """
    Update the connectivity values for each atom in the molecule.
    """
    cython.declare(atom=Atom, value=int, connectivityValues=list)

    connectivityValues = mol.update_recurse([atom.number * len(atom.bonds) for atom in mol.atoms], 0)
    for atom, value in zip(mol.atoms, connectivityValues):
        atom.connectivity = value

def normalize(mol):
    """
    Select the resonance isomer that is isomorphic to the parameter isomer, with the lowest unpaired
    electrons descriptor.

    We generate over all resonance isomers (non-isomorphic as well as isomorphic) and add isomorphic
    isomers to a list (candidates).

    Next, we search through this list and return the candidate with the lowest unpaired electrons descriptor.

    """
    candidates = [mol]# resonance isomers that are isomorphic to the parameter isomer.

    isomers = [mol]

    # Iterate over resonance isomers
    index = 0
    while index < len(isomers):
        isomer = isomers[index]
            
        newIsomers = isomer.getAdjacentResonanceIsomers()
        newIsomers += isomer.getLonePairRadicalResonanceIsomers()
        newIsomers += isomer.getN5dd_N5tsResonanceIsomers()
        newIsomers += isomer.getKekulizedResonanceIsomers()

        for newIsomer in newIsomers:
            newIsomer.updateAtomTypes()
            # Append to isomer list if unique
            for isom in isomers:
                isom_copy = isom.copy(deep=True)
                newIsomer_copy = newIsomer.copy(deep=True)
                if isom_copy.isIsomorphic(newIsomer_copy):
                    candidates.append(newIsomer)
                    break
            else:
                isomers.append(newIsomer)        
                    
        # Move to next resonance isomer
        index += 1
    
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

def find_mobile_h_system(mol, all_mobile_h_atoms_couples, test_indices):
    dummy = test_indices[:]

    for mobile_h_atom_couple in all_mobile_h_atoms_couples:
        for test_index in test_indices:
            if test_index in mobile_h_atom_couple:
                original_atom = test_index
                dummy.remove(test_index)
                mobile_h_atom_couple.remove(test_index)
                new_partner = mobile_h_atom_couple[0]
                central = dummy[0]
                return mol.atoms[central - 1], mol.atoms[original_atom - 1], mol.atoms[new_partner - 1]

    raise Exception('We should always have found the mobile-H system. All mobile H couples: {}, test indices: {}'
        .format(all_mobile_h_atoms_couples, test_indices))


def check_bond_order_oxygen(mol):
    """Check if total bond order of oxygen atoms is smaller than 4."""
    from rmgpy.molecule.util import ORDERS

    for at in mol.atoms:
        if at.number == 8:
            order = sum([ORDERS[b.order] for _, b in at.bonds.iteritems()])
            not_correct = order >= 4
            if not_correct:
                return False

    return True

def parse_N_layer(auxinfo):
    """
    Parses the layer with atom ordering information (N-layer) 
    and returns a list of atom indices that reflect how the atoms of the original
    molecule should be ordered according to the InChI algorithm.


    Example:
    Auxiliary info of SMILES OCCC (InChI=1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3):
    AuxInfo=1/0/N:4,3,2,1/rA:4OCCC/rB:s1;s2;s3;/rC:;;;;

    N-layer: 
    /N:4,3,2,1

    The original number of an atom with identification number n is given as the
    n-th member of this list for a component; the lists are separated with “;”. 

    Raises an exception when the N-layer could not be found.
    """

    pieces = auxinfo.split('/')
    atom_numbers = None
    for piece in pieces:
        if piece.startswith('N'):
            atom_numbers = piece[2:]#cut off N:
            break
    else:
        raise Exception('Could not find the N-layer in the auxiliary info: {}'.format(auxinfo))

    indices = map(int, atom_numbers.split(','))

    return indices


def parse_E_layer(auxinfo):
    """
    Converts the layer with equivalence information (E-layer) 
    on atoms into a list of lists of equivalent atom indices.

    Example:
    Auxiliary info of InChI=1S/C8H14/c1-5-7(3)8(4)6-2/h5-8H,1-2H2,3-4H3:
    AuxInfo=1/0/N:1,8,4,6,2,7,3,5/E:(1,2)(3,4)(5,6)(7,8)/rA:8C.2C.2CCCCCC/rB:s1;s2;s3;s3;s5;s5;d7;/rC:;;;;;;;;
    E-layer: 

    /E:(1,2)(3,4)(5,6)(7,8)/

    denotes that atoms (1,2), (3,4), (5,6), (7,8) are equivalent and cannot be distinguished based on the 
    implemented canonicalization algorithm.

    Returned object:

    [[1,2],[3,4],[5,6],[7,8]]

    Returns an empty list of the E-layer could not be found.

    """

    pieces = auxinfo.split('/')
    e_layer = None
    for piece in pieces:
        if piece.startswith('E'):
            e_layer = piece[2:]#cut off /E:
            break
    else:
        return []

    # search for (*) pattern
    pattern = re.compile( r'\((.[^\(\)]*)\)')

    equivalent_atoms = []
    for atomtuple in re.findall(pattern, e_layer):
        indices = list(map(int, atomtuple.split(',')))
        equivalent_atoms.append(indices)

    return equivalent_atoms

def group_adjacent_unpaired_electrons(mol, u_layer, equivalent_atoms):
    """
    Iterates through the list of unpaired electrons and forms couples of 
    adjancent unpaired electrons.

    If an adjancent unpaired electron cannot be found, the list is 
    expanded with a the single element, instead of the couple.

    Every time two adjacent electrons have been detected, they are removed
    from the u-layer. The loop is continued until all electrons have been either
    assigned to a partner or until all single electrons have been removed.
    """
    adjacent_electrons = []

    u_layer_copy = u_layer[:]

    while u_layer_copy:
        i = u_layer_copy.pop()
        if not any(i in x for x in equivalent_atoms):
            # add the atom by itself:
            adjacent_electrons.append([i])
            continue

        # iterate over neighbors and check if neighbor is in u_layer
        for at, bond in mol.atoms[i-1].bonds.iteritems():
            at_index = mol.atoms.index(at)+1
            
            if any(i in x for x in equivalent_atoms):
                if at_index in u_layer_copy:
                    pair = sorted([i, at_index])
                    adjacent_electrons.append(pair)
                    u_layer_copy.remove(at_index)
                    break
        else:
            adjacent_electrons.append([i])
            
    return adjacent_electrons


def generate_combos(unpaired_electrons, equivalent_atoms):
    """
    Generate all possible combinations of groups of unpaired electrons
    based on the information of the equivalent atoms.

    Example:

    unpaired electrons: [1,2]
    equivalent atoms: [[1,3], [2,4]]

    what is returned: the cross product of both lists:
    [[1,2], [1,4], [3,2], [3,4]]

    if the unpaired electrons appear in the same list of equivalent atoms, 
    [[1,2,3,4]]

    the combinations of two from that list is returned instead:
    [[1,2], [1,4], [3,2], [3,4]]

    Returns an empty list of one of the unpaired electrons could not be found
    in the list of equivalent atoms.
    """

    # assign the right list with equivalent atoms for each electron:
    corr_eq_atoms = []
    for electron in unpaired_electrons:
        for eq in equivalent_atoms:
            if electron in eq:
                corr_eq_atoms.append(eq)
                break
        else: return []

    if corr_eq_atoms[1:] == corr_eq_atoms[:-1]:#check if all elements in a list are identical
        combos = [list(tup) for tup in itertools.combinations(corr_eq_atoms[0], len(unpaired_electrons))]
    else:#cross product of the elements in the list
        combos = [list(tup) for tup in itertools.product(*corr_eq_atoms)]
    
    return combos

def valid_combo(combo, mol, u_layer):
    """
    Check if the combination of atom indices refers to
    atoms that are adjacent in the molecule.
    """
    assert len(combo) < 3

    if len(combo) == 1 and len(u_layer) > 1: return False

    atoms = [mol.atoms[index-1] for index in combo]

    conditions = []
    conditions.append(any([at.radicalElectrons == 0 for at in atoms]))
    if len(atoms) == 2:
        at1, at2 = atoms
        conditions.append(mol.hasBond(at1, at2))
    
    return all(conditions)

def find_lowest_u_layer(mol, u_layer, equivalent_atoms):
    """..."""

    if not equivalent_atoms:
        return u_layer

    new_u_layer = []
    groups = group_adjacent_unpaired_electrons(mol, u_layer, equivalent_atoms)

    for group in groups:
        selected_group = group
        combos = generate_combos(group, equivalent_atoms)
        combos = [combo for combo in combos if combo != group]
        for combo in combos:
            combo = sorted(combo)
            if valid_combo(combo, mol, u_layer):
                if combo < selected_group:
                    selected_group = combo

        new_u_layer.extend(selected_group)

    return sorted(new_u_layer)
    
