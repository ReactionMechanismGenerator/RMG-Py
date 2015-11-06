import cython
import re

# search for (*) PARENTHESES
PARENTHESES = re.compile( r'\((.[^\(\)]*)\)')

INCHI_PREFIX = 'InChI=1'


"""
The prefix with the information on the distribution of unpaired electrons across the atoms.

For example, a triplet biradical with unpaired electrons on atom 1 and atom 3 
will have the following unpaired electron layer in the augmented InChI:

InChI=1/.../u1,3

The indices refer to the 1-based indices in the InChI string (NOT the 0-based
    indices of the Molecule container!)

"""
U_LAYER_PREFIX = '/u'

"""The separator that separates the indices of the atoms that bear unpaired electrons."""
U_LAYER_SEPARATOR = ','


"""
The prefix with the information on the distribution of the atoms 
with an unexpected number of lone pairs.

For example, a singlet methylene with a lone pair on atom 1
will have the following lone pair layer in the augmented InChI:

InChI=1/.../p1

The indices refer to the 1-based indices in the InChI string (NOT the 0-based
    indices of the Molecule container!)

"""
P_LAYER_PREFIX = '/p'

"""The separator that separates the indices of the atoms that bear unpaired electrons."""
P_LAYER_SEPARATOR = ','

ulayer_pattern = re.compile(r'/u(.*)')
player_pattern = re.compile(r'/p(.*)')

def decompose(string):
    """
    Converts an augmented inchi into 
    - an inchi, 
    - indices array for the atoms bearing unpaired electrons.
    - indices array for the atoms bearing (unexpected) lone pairs.


    """
    cython.declare(
            inchi=str,
            u_indices=list,
            p_indices=list,
        )

    if U_LAYER_PREFIX in string:
        inchi = string.split(U_LAYER_PREFIX)[0]
    elif P_LAYER_PREFIX in string:
        inchi = string.split(P_LAYER_PREFIX)[0]
    else:
        inchi = string

    u_indices, p_indices = [], []
    matches = re.findall(ulayer_pattern, string)
    if matches:
        u_indices = map(int, matches.pop().split('/')[0].split(U_LAYER_SEPARATOR))

    matches = re.findall(player_pattern, string)
    if matches:
        p_indices = map(int, matches.pop().split('/')[0].split(P_LAYER_SEPARATOR))

    return inchi, u_indices, p_indices

def ignore_prefix(string):
    """
    Splits off the 'InChI=1S' or 'InChI=1' layer of an InChI
    and returns the last part.
    """

    if not INCHI_PREFIX in string:
        raise InchiException('Not a valid InChI: {}'.format(string))

    return re.split(r"(InChI=1+)(S*)/", string)[-1]

def compose_aug_inchi(inchi, ulayer=None, player=None):
    """
    Composes an augmented InChI by concatenating the different pieces
    as follows:

    InChI=1S/XXXX.../c.../h.../ux,x,...
    """
    cython.declare(
            temp=str,
        )

    aug_inchi = INCHI_PREFIX + '/' if not INCHI_PREFIX in inchi else ''
    aug_inchi += inchi
    
    for layer in filter(None, [ulayer, player]):
        aug_inchi += layer

    return aug_inchi

def compose_aug_inchi_key(inchi_key, ulayer=None, player=None):
    """
    Composes an augmented InChI Key by concatenating the different pieces
    as follows:

    XXXXXXXXXXXXXX-XXXXXXXXXX-ux,x,xxx

    Uses hyphens rather than forward slashes to avoid messing up file paths.
    """

    aug_inchi_key = inchi_key

    for layer in filter(None, [ulayer, player]):
        aug_inchi_key += '-' + layer[1:]#cut off the '/'

    return aug_inchi_key 

def parse_H_layer(inchi):
    """
    Converts the Mobile H layer of an inchi string into a 
    list of atom indices couples that carry a mobile hydrogen.

    Example:
    The hydrogen of the hydroxyl group can migrate to the carbonyl 
    oxygen.

    O=C-O
    InChI=1S/CH2O2/c2-1-3/h1H,(H,2,3)

    The atoms bearing a mobile hydrogen will be found in
    the hydrogen layer, within the parentheses.

    An empty list will be returned when there are no mobile hydrogens.

    """

    cython.declare(
            pieces=list,
            h_layer=str,
            piece=str,
            couples=list,
            match=str,
            mobile_h_atoms=list,
        )


    pieces = inchi.split('/')
    h_layer = None
    for piece in pieces:
        if piece.startswith('h'):
            h_layer = piece
            break
    else: 
        raise Exception('Could not find the hydrogen layer in the inchi: {}'.format(inchi))

    couples = []
    for match in re.findall(PARENTHESES, h_layer):
        mobile_h_atoms = map(int, match[2:].split(','))
        couples.append(mobile_h_atoms)

    return couples

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

    cython.declare(
            pieces=list,
            e_layer=str,
            piece=str,
            equivalent_atoms=list,
            atomtuple=str,
            indices=list,
        )

    pieces = auxinfo.split('/')
    e_layer = None
    for piece in pieces:
        if piece.startswith('E'):
            e_layer = piece[2:]#cut off /E:
            break
    else:
        return []

    equivalent_atoms = []
    for atomtuple in re.findall(PARENTHESES, e_layer):
        indices = list(map(int, atomtuple.split(',')))
        equivalent_atoms.append(indices)

    return equivalent_atoms

 
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


    cython.declare(
            pieces=list,
            atom_numbers=str,
            piece=str,
            indices=list,
        )

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
        
class InchiException(Exception):
    pass

class InChI(str):
    """InChI is a type of string in which the InChI=1 prefix is ignored."""
    def __new__(self, inchi):

      if not INCHI_PREFIX in inchi:
        raise InchiException('Not a valid InChI: {}'.format(inchi))

      return str.__new__(self, ignore_prefix(inchi))

class AugmentedInChI(InChI):
    """AugmentedInChI is an InChI with inchi, and unpaired electron attributes."""
    def __init__(self, aug_inchi):
        super(AugmentedInChI, self).__init__()
        inchi, u_indices, p_indices = decompose(str(self))

        self.inchi = str(inchi)

        # default to None
        self.u_indices = u_indices or None
        self.p_indices = p_indices or None
        
