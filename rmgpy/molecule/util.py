import re
import itertools

from rmgpy.molecule.molecule import Molecule


"""
Number of valence electrons of the elements
"""
VALENCES = {'H': 1, 'He': 2, 'C': 4, 'N': 5, 'O': 6,  'Ne': 8, 
            'Si': 4, 'S': 6, 'Cl': 7, 'Ar': 8
            }

"""
Bond orders of single, double, triple bonds.
"""            
ORDERS = {'S': 1, 'D': 2, 'T': 3}

def retrieveElementCount(obj):
    """Converts an (augmented) inchi or Molecule into a dictionary element -> count"""
    element_count = {}

    if isinstance(obj, str):
        assert 'InChI=1' in obj
        mf = obj.split('/')[1]
        pieces = re.findall('[A-Z][^A-Z]*', mf)#split on capital letters

        for piece in pieces:
            match = re.match(r"([a-z]+)([0-9]*)", piece, re.I)
            if match:
                element, count = match.groups()
                if count is '':
                    count = 1
                element_count[element] = int(count)
        return element_count
    
    elif isinstance(obj, Molecule):
        for atom in obj.atoms:
            element = atom.element.symbol
            if element in element_count:
                updated_count = element_count[element] + 1
                element_count[element] = updated_count 
            else:
                element_count[element] = 1
        return element_count

    else:
        raise Exception

def partition(sample, list_of_samples):
    """
    Group indices from the parameter sample 
    that belong to the same list.

    Returns a list of lists with partitioned indices, and 
    a list of lists with the corresponding sample list they were found in.

    E.g.:
    sample : [1,3]
    list_of_samples : [[1,2,3,4],[5,6]]
    returns: [[1,3]], [[1,2,3,4]]

    Indices not part of any of the lists should be in singleton list and 
    have corresponding empty list:

    E.g.:
    sample : [7]
    list_of_samples : [[1,2,3,4],[5,6]]
    returns: [[7]], [[]]

    """

    partitions, sample_lists  = [], []

    for s in sample:
        for one_sample_list in list_of_samples:
            if s in one_sample_list:
                try:
                    index = sample_lists.index(one_sample_list)
                    partitions[index].append(s)
                except ValueError:
                    partitions.append([s])
                    sample_lists.append(one_sample_list)                    
                
                break
        else:# s does not belong to any list of samples
            partitions.append([s])
            sample_lists.append([])

    return partitions, sample_lists          


def agglomerate(groups):
    """
    Iterates over the parameter list of lists, and identifies all singletons.
    A new list of lists is created in which all singletons are combined together,
    while the other lists consisting of more than 1 element are simply copied.
    The newly created collapsed list singletons is appended at the end.

    Example:

    [[1,2,3], [4], [5,6], [7]]

    Returns:
    [[1,2,3], [5,6], [4,7]]
    """

    result = filter(lambda x: len(x) > 1, groups)
    singletons = filter(lambda x: len(x) == 1, groups)
    
    # collapse singletons:
    singletons = list(itertools.chain.from_iterable(singletons))

    result.append(singletons)
    return result

def generate_combo(samples, sample_spaces):
    """
    First, generate combinations of i samples from the corresponding sample_spaces.
    Next, generate the cross product between the elements in the previously generated list.

    Finally, filter out the original combination from the result.
    """
    combos = []
    for sample, sample_space in zip(samples, sample_spaces):
        combos_one_sample = [list(tup) for tup in itertools.combinations(sample_space, len(sample))]
        combos.append(combos_one_sample)

    # cross product of the elements in the list:
    combos = [list(tup) for tup in itertools.product(*combos)]

    # don't add the original samples
    combos = filter(lambda x: x != samples, combos)

    return combos    

def swap(to_be_swapped, sample):
    """
    Identifies which index of the list samples  is present in
    the list to be swapped. 

    E.g.:
    to be swapped: [2,3]
    sample: [1,3]

    Returns: 
    1, 3, 2
    
    """

    to_be_swapped = set(to_be_swapped)
    sample = set(sample)

    original = (sample.intersection(to_be_swapped)).pop()
    central = (sample - to_be_swapped).pop()
    new_partner = (to_be_swapped - sample).pop()
    
    return central, original, new_partner