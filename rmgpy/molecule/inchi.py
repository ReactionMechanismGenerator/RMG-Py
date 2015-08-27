import re

INCHI_PREFIX = 'InChI=1'

"""

The prefix with the multiplicity information.

For example, a triplet biradical will bear the 
following multiplicity layer in the augmented InChI:

InChI=1/.../mult3

The singlet equivalent of aforementioned biradical has the following
augmented InChI:

InChI=1/.../mult1

"""
MULT_PREFIX = '/mult'



"""
The prefix with the information on the distribution of unpaired electrons across the atoms.

For example, a triplet biradical with unpaired electrons on atom 1 and atom 3 
will have the following unpaired electron layer in the augmented InChI:

InChI=1/.../mult3/u1,3

The indices refer to the 1-based indices in the InChI string (NOT the 0-based
    indices of the Molecule container!)

"""
U_LAYER_PREFIX = '/u'

"""The separator that separates the indices of the atoms that bear unpaired electrons."""
U_LAYER_SEPARATOR = ','

def decompose(string):
    """
    Converts an augmented inchi into an inchi, multiplicity and indices array for the atoms
    bearing unpaired electrons.

    returns: 
    - multiplicity (int)
    - indices array for the atoms bearing unpaired electrons.

    """

    if not MULT_PREFIX in string:
        return string, None, []

    inchi, mult_ulayer = string.split(MULT_PREFIX)
    
    if not U_LAYER_PREFIX in mult_ulayer:
        mult = int(mult_ulayer)
        u_indices = []    
        return inchi, mult, u_indices

    mult, ulayer = mult_ulayer.split(U_LAYER_PREFIX)
    mult = int(mult)
    u_indices = [int(i) for i in ulayer.split(U_LAYER_SEPARATOR)]
    return inchi, mult, u_indices

def ignore_prefix(string):
    """
    Splits off the 'InChI=1S' or 'InChI=1' layer of an InChI
    and returns the last part.
    """

    if not INCHI_PREFIX in string:
        raise InchiException('Not a valid InChI: {}'.format(string))

    return re.split(r"(InChI=1+)(S*)/", string)[-1]

def compose_aug_inchi(inchi, mult_layer, ulayer=None):
    prefix = INCHI_PREFIX + '/' if not INCHI_PREFIX in inchi else ''
    if ulayer is not None:
        return prefix + inchi + mult_layer + ulayer
    else:
        return prefix + inchi + mult_layer 


def compose_aug_inchi_key(inchi_key, mult_layer, ulayer=None):
    if ulayer is not None:
        return inchi_key + mult_layer + ulayer
    else:
        return inchi_key + mult_layer 

class InchiException(Exception):
    pass

class InChI(str):
    """InChI is a type of string in which the InChI=1 prefix is ignored."""
    def __new__(self, inchi):

      if not INCHI_PREFIX in inchi:
        raise InchiException('Not a valid InChI: {}'.format(inchi))

      return str.__new__(self, ignore_prefix(inchi))

class AugmentedInChI(InChI):
    """AugmentedInChI is an InChI with inchi, multiplicity and unpaired electron attributes."""
    def __init__(self, aug_inchi):
        super(AugmentedInChI, self).__init__()
        inchi, mult, u_indices = decompose(self)

        self.inchi = str(inchi)

        # default to multiplicity -1
        self.mult = mult or -1

        # default to None
        self.u_indices = u_indices or None
        