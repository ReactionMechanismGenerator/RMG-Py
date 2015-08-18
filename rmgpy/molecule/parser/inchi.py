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

For closed shell molecules (multiplicity 1) and monoradicals (multiplicity 2),
the multiplicity layer is optional.

"""
MULT_PREFIX = '/mult'



"""
The prefix with the information on the distribution of unpaired electrons across the atoms.

For example, a triplet biradical with unpaired electrons on atom 1 and atom 3 
will have the following unpaired electron layer in the augmented InChI:

InChI=1/.../mult3/u1,3

For monoradicals (multiplicity 2),
the multiplicity layer is optional.

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

    assert INCHI_PREFIX in string, 'Not a valid InChI: {}'.format(string)
    return re.split(r"(InChI=1+)(S*)/", string)[-1]

class InChI(str):
    """InChI is a type of string in which the InChI=1 prefix is ignored."""
    def __new__(self, inchi):
      assert INCHI_PREFIX in inchi, 'Not a valid InChI: {}'.format(inchi)
      return str.__new__(self, ignore_prefix(inchi))

class AugmentedInChI(InChI):
    """AugmentedInChI is an InChI with inchi, multiplicity and unpaired electron attributes."""
    def __init__(self, aug_inchi):
        super(AugmentedInChI, self).__init__()
        inchi, mult, u_indices = decompose(self)

        self.inchi = inchi

        # default to multiplicity 1
        self.mult = mult or 1

        # default to None
        self.u_indices = u_indices or None
        