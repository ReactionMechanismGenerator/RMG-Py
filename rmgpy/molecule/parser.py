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
from .molecule import Atom, Bond, Molecule

import rmgpy.molecule.inchi as inchiutil
import rmgpy.molecule.util as util
import rmgpy.molecule.pathfinder as pathfinder

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

def reset_lone_pairs_to_default(mol):
    """Resets the atom's lone pair count to its default value."""

    for at in mol.atoms:
        order = sum([ORDERS[b.order] for _,b in mol.getBonds(at).iteritems()])
        at.lonePairs = (VALENCES[at.element.symbol] - order - at.radicalElectrons - at.charge) / 2

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

    isFixed = False
    for u1, u2 in combos:
        atom1 = mol.atoms[u1 - 1] # convert to 0-based index for atoms in molecule
        atom2 = mol.atoms[u2 - 1] # convert to 0-based index for atoms in molecule
        if mol.hasBond(atom1, atom2):
            b = mol.getBond(atom1, atom2)
            isFixed = fix_unsaturated_bond(b)
            if isFixed:
                break
                
            else:
                isFixed = fix_mobile_h(mol, inchi, u1, u2)
                if isFixed:
                    break
        else:
            isFixed = convert_butadiene_path(atom1, atom2)
            if isFixed:
                break

    if isFixed:
        u_indices.remove(u1)
        u_indices.remove(u2)
        return mol                
    else:
        raise Exception(
            'Could not convert an unsaturated bond into a biradical for the \
            indices {} provided in the molecule: {}.'
            .format(u_indices, mol.toAdjacencyList())
            )    

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
        inchi_elementcount = util.retrieveElementCount(identifier)
        mol_elementcount = util.retrieveElementCount(mol)
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

def check(mol, aug_inchi) :
    """
    Check if the molecular structure is correct.

    Checks whether the multiplicity contained in the augmented inchi, 
    corresponds to the number of unpaired electrons + 1 found in the molecule.

    Checks whether the valence of each atom is compatible with the bond order,
    number of unpaired electrons, lone pairs and charge.

    """
    cython.declare(inchi=str,
                   multi=cython.int,
                   at=Atom
                   )

    assert mol.multiplicity == mol.getRadicalCount() + 1,\
     'Multiplicity of molecule \n {0} does not correspond to aug. inchi {1}'.format(mol.toAdjacencyList(), aug_inchi)
    
    for at in mol.atoms:
        order = sum([util.ORDERS[b.order] for _,b in mol.getBonds(at).iteritems()])
        assert (order + at.radicalElectrons + 2*at.lonePairs + at.charge) == util.VALENCES[at.symbol],\
            'Valency for an atom of molecule \n {0} does not correspond to aug. inchi {1}'.format(mol.toAdjacencyList(), aug_inchi)

def fix_oxygen_unsaturated_bond(mol, u_indices):
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
            paths = pathfinder.find_allyl_end_with_charge(start)
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

    if inchiutil.INCHI_PREFIX in inchistr:
        return __parse(mol, inchistr, 'inchi', backend)
    else:
        return __parse(mol, inchiutil.INCHI_PREFIX + '/' + inchistr, 'inchi', backend)



def fromAugmentedInChI(mol, aug_inchi):
    """
    Creates a Molecule object from the augmented inchi.

    First, the inchi is converted into a Molecule using
    the backend parsers.

    Next, the multiplicity and unpaired electron information
    is used to fix a number of parsing errors made by the backends.

    Finally, the atom types of the corrected molecule are perceived.

    Returns a Molecule object
    """

    if not isinstance(aug_inchi, inchiutil.AugmentedInChI):
        aug_inchi = inchiutil.AugmentedInChI(aug_inchi)

    mol = fromInChI(mol, aug_inchi.inchi)

    mol.multiplicity = len(aug_inchi.u_indices) + 1 if aug_inchi.u_indices else 1

    fix(mol, aug_inchi)

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

def fixCharge(mol, u_indices):
    """
    Tries to fix a number of structural features in the molecule related to charge, 
    based on the information from the parameter list of atom indices with unpaired electrons.
    """

    if not u_indices:
        return

    is_charged = sum([abs(at.charge) for at in mol.atoms]) != 0
    is_correct = mol.getNumberOfRadicalElectrons() == (mol.multiplicity - 1)
    if mol.multiplicity < 3 or is_correct or not is_charged:
        return

    # converting charges to unpaired electrons for atoms in the u-layer
    convert_charge_to_unpaired_electron(mol, u_indices)

    # convert neighboring atoms (or delocalized paths) to unpaired electrons
    convert_delocalized_charge_to_unpaired_electron(mol, u_indices)

    fix_adjacent_charges(mol)

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

def find_mobile_h_system(mol, all_mobile_h_atoms_couples, test_indices):
    """
    
    """
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
    
def fix_adjacent_charges(mol):
    """
    Searches for pairs of charged atoms.
    Neutralizes one unit of charge on each atom,
    and increments the bond order of the bond in between
    the atoms.
    """
    for at in mol.atoms:
        if at.charge != 0:
            for neigh, bond in at.bonds.iteritems():
                if neigh.charge != 0:
                    bond.incrementOrder()
                    at.charge += 1 if at.charge < 0 else -1
                    neigh.charge += 1 if neigh.charge < 0 else -1

def convert_charge_to_unpaired_electron(mol, u_indices):
    """
    Iterates over the atoms foundin the parameter list and
    converts a unit of charge on atoms into an unpaired electron.

    Removes treated atoms from the parameter list.
    """
    for at in mol.atoms:
        at_index = mol.atoms.index(at) + 1
        if at.charge != 0 and at_index in u_indices:
            at.charge += 1 if at.charge < 0 else -1
            at.radicalElectrons += 1
            u_indices.remove(at_index)                    

def convert_delocalized_charge_to_unpaired_electron(mol, u_indices):
    """
    Iterates over the atom indices of the parameter list and searches 
    a charged atom that is connected to that atom via some kind of
    delocalization path.

    """
    u_indices_copy = u_indices[:]
    for index in u_indices_copy:
        start = mol.atoms[index -1]

        found = convert_4_atom_3_bond_path(start)
        if found: 
            u_indices.remove(index)
            continue

        found = convert_3_atom_2_bond_path(start, mol)
        if found:
            u_indices.remove(index)
            continue

def convert_4_atom_3_bond_path(start):
    """
    Searches for 4-atom-3-bond [X=X-X=X+] paths starting from the parameter atom.
    If a path is found, the starting atom receives an unpaired electron while
    the bonds in the delocalization path are "inverted". A unit of charge on the 
    end atom is neutralized and a lone pair is added.
    """
    path = pathfinder.find_butadiene_end_with_charge(start)

    if path is not None:    
        start.radicalElectrons += 1
        end = path[-1]
        end.charge += 1 if end.charge < 0 else -1
        end.lonePairs += 1

        # filter bonds from path and convert bond orders:
        bonds = path[1::2]#odd
        for bond in bonds[::2]:# even
            assert isinstance(bond, Bond)
            bond.decrementOrder()
        for bond in bonds[1::2]:# odd bonds
            assert isinstance(bond, Bond)
            bond.incrementOrder()  

        return True

    return False

def convert_3_atom_2_bond_path(start, mol):
    """
    Searches for 3-atom-2-bond [X=X-X+] paths paths starting from the parameter atom.
    If a correct path is found, the starting atom receives an unpaired electron while
    the bonds in the delocalization path are "inverted". A unit of charge on the 
    end atom is neutralized and a lone pair is added.

    If it turns out the path was invalid, the actions are reverted, and another path
    is tried instead.

    To facilitate reverting the changes, we use a reaction recipe and populate it
    with a number of actions that reflect the changes in bond orders and unpaired
    electrons that the molecule should undergo.
    """
    from rmgpy.data.kinetics.family import ReactionRecipe

    def is_valid(mol):
        """Check if total bond order of oxygen atoms is smaller than 4."""

        for at in mol.atoms:
            if at.number == 8:
                order = sum([util.ORDERS[b.order] for _, b in at.bonds.iteritems()])
                not_correct = order >= 4
                if not_correct:
                    return False

        return True

    index = mol.atoms.index(start) + 1

    paths = pathfinder.find_allyl_end_with_charge(start)

    for path in paths:
        # label atoms so that we can use the labels in the actions of the recipe
        for i, at in enumerate(path[::2]):
            at.label = str(i)
        # we have found the atom we are looking for
        recipe = ReactionRecipe()
        recipe.addAction(['GAIN_RADICAL', start.label, 1])

        end = path[-1]
        end_original_charge = end.charge
      
        # filter bonds from path and convert bond orders:
        bonds = path[1::2]#odd elements
        for bond in bonds[::2]:# even
            recipe.addAction(['CHANGE_BOND', bond.atom1.label, -1, bond.atom2.label])
        for bond in bonds[1::2]:# odd
            recipe.addAction(['CHANGE_BOND', bond.atom1.label, 1, bond.atom2.label])

        end.charge += 1 if end.charge < 0 else -1
        recipe.applyForward(mol, update=False)

        if is_valid(mol):
            # unlabel atoms so that they never cause trouble downstream
            for i, at in enumerate(path[::2]):
                at.label = ''
            return True
        else:
            recipe.applyReverse(mol, update=False)
            end.charge = end_original_charge

            # unlabel atoms so that they never cause trouble downstream
            for i, at in enumerate(path[::2]):
                assert isinstance(at, Atom)
                at.label = ''

    return False

def fix(mol, aug_inchi):
    """
    Fixes a number of structural features of the erroneous Molecule
    parsed by the backends, based on multiplicity and unpaired electron information
    stored in the augmented inchi.
    """   

    u_indices = aug_inchi.u_indices[:] if aug_inchi.u_indices else []
    p_indices = aug_inchi.p_indices[:] if aug_inchi.p_indices else []

    # ignore atoms that bear already unpaired electrons:
    for i in set(u_indices[:]):
        atom = mol.atoms[i - 1]
        [u_indices.remove(i) for _ in range(atom.radicalElectrons)]        

    # ignore atoms that bear already lone pairs:
    for i in set(p_indices[:]):
        atom = mol.atoms[i - 1]
        [p_indices.remove(i) for _ in range(atom.lonePairs)]   


    fix_triplet_to_singlet(mol, p_indices)
    
    fixCharge(mol, u_indices)
                                
    reset_lone_pairs(mol, p_indices)

    fix_oxygen_unsaturated_bond(mol, u_indices)

    fix_unsaturated_bond(mol, u_indices, aug_inchi)

    check(mol, aug_inchi)    


def fix_triplet_to_singlet(mol, p_indices):
    """
    Iterates over the atoms and checks whether atoms bearing two unpaired electrons are
    also present in the p_indices list.

    If so, convert to the two unpaired electrons into a lone pair, and remove that atom
    index from the p_indices list.
    """

    for at in mol.atoms:
        index = mol.atoms.index(at) + 1
        if mol.getNumberOfRadicalElectrons() == 2 and index in p_indices:
            at.lonePairs += 1
            at.radicalElectrons -= 2
            p_indices.remove(index)


def fix_butadiene_path(start, end):
    """
    Searches for a 1,3-butadiene path between the start and end atom.
    Adds an unpaired electron to start and end atom, and "inverts" the bonds
    in between them. 
    """
    path = pathfinder.find_butadiene(start, end)
    if path is not None:
        start.radicalElectrons += 1
        end.radicalElectrons += 1
        # filter bonds from path and convert bond orders:
        bonds = path[1::2]#odd elements
        for bond in bonds[::2]:# even bonds
            assert isinstance(bond, Bond)
            bond.decrementOrder()
        for bond in bonds[1::2]:# odd bonds
            assert isinstance(bond, Bond)
            bond.incrementOrder()    

        return True

    return False

def fix_mobile_h(mol, inchi, u1, u2):
    """

    Identifies a system of atoms bearing unpaired electrons and mobile hydrogens
    at the same time.

    The system will consist of a central atom that does not bear any mobile hydrogens,
    but that is bound to an atom that does bear a mobile hydrogen, called the "original atom".

    The algorithm identifies the "new partner" atom that is part of the mobile hydrogen 
    system.

    Next, the mobile hydrogen is transferred from the original atom, to the new partner,
    and a bond is removed and added respectively.

    Finally, the central atom and the original atom will each receive an unpaired electron,
    and the bond between them will decrease in order.
    """

    mobile_hydrogens = inchiutil.parse_H_layer(inchi)

    if mobile_hydrogens:
        # WIP: only consider the first system of mobile hydrogens:
        mobile_hydrogens = mobile_hydrogens[0]

        #find central atom:
        central, original, new_partner = util.swap(mobile_hydrogens, [u1, u2])

        central, original, new_partner = \
        mol.atoms[central - 1], mol.atoms[original - 1], mol.atoms[new_partner - 1]

        # search hydrogen atom and bond
        hydrogen = None
        for at, bond in original.bonds.iteritems():
            if at.number == 1:
                hydrogen = at
                mol.removeBond(bond)
                break

        new_h_bond = Bond(new_partner, hydrogen, order='S')
        mol.addBond(new_h_bond)
        
        mol.getBond(central, new_partner).decrementOrder()

        central.radicalElectrons += 1
        original.radicalElectrons += 1
        return True

    return False

def convert_unsaturated_bond_to_triplet(bond):
    """
    Decrements the bond if it is unsatured, and adds an unpaired
    electron to each of the atoms connected by the bond.
    """
    if not bond.isSingle():
        for at in (bond.atom1, bond.atom2):
            at.radicalElectrons += 1
        bond.decrementOrder()    
        return True
    return False

def reset_lone_pairs(mol, p_indices):
    """
    Iterates over the atoms of the molecule and
    resets the atom's lone pair count to the value stored in the p_indices list,
    or to the default value.

    """

    for at in mol.atoms:
        index = mol.atoms.index(at) + 1 #1-based index
        count = p_indices.count(index)
        if count != 0:
            at.lonePairs = count
        else:    
            order = sum([util.ORDERS[b.order] for _,b in mol.getBonds(at).iteritems()])
            at.lonePairs = (util.VALENCES[at.symbol] - order - at.radicalElectrons - at.charge) / 2

def fix_unsaturated_bond_to_biradical(mol, inchi, u_indices):
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

    isFixed = False
    for u1, u2 in combos:
        atom1 = mol.atoms[u1 - 1] # convert to 0-based index for atoms in molecule
        atom2 = mol.atoms[u2 - 1] # convert to 0-based index for atoms in molecule
        if mol.hasBond(atom1, atom2):
            b = mol.getBond(atom1, atom2)
            isFixed = convert_unsaturated_bond_to_triplet(b)
            if isFixed:
                break
                
            else:
                isFixed = fix_mobile_h(mol, inchi, u1, u2)
                if isFixed:
                    break
        else:
            isFixed = fix_butadiene_path(atom1, atom2)
            if isFixed:
                break

    if isFixed:
        u_indices.remove(u1)
        u_indices.remove(u2)
        return mol                
    else:
        raise Exception(
            'Could not convert an unsaturated bond into a biradical for the \
            indices {} provided in the molecule: {}.'
            .format(u_indices, mol.toAdjacencyList())
            )    

def isUnsaturated(mol):
    """
    Does the molecule have a bond that's not single?
    Eg. a bond that is double or triple or benzene
    """
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


def fix_unsaturated_bond(mol, indices, aug_inchi):
    """
    Adds unpaired electrons to the molecule by converting unsaturated bonds into triplets.

    It does so by converting an unsaturated bond into a triplet, and verifying whether
    the total number of unpaired electrons matches the multiplicity.

    Finishes when all unsaturated bonds have been tried, or when there are no pairs
    of atoms that should be unpaired electrons left. 
    """

    correct = mol.getNumberOfRadicalElectrons() == (mol.multiplicity - 1)
    
    if not correct and not indices:
        raise Exception( 'Cannot correct {} based on {} by converting unsaturated bonds into unpaired electrons...'\
            .format(mol.toAdjacencyList(), aug_inchi))

    unsaturated = isUnsaturated(mol)

    while not correct and unsaturated and len(indices) > 1:
        mol = fix_unsaturated_bond_to_biradical(mol, aug_inchi.inchi, indices)
        correct = mol.getNumberOfRadicalElectrons() == (mol.multiplicity - 1)
        unsaturated = isUnsaturated(mol)
