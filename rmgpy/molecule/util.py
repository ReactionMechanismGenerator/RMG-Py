import re

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

    for index in sample:
        for one_sample_list in list_of_samples:
            if index in one_sample_list:
                try:
                    index = sample_lists.index(one_sample_list)
                    partitions[index].append(index)
                except ValueError:
                    partitions.append([index])
                    sample_lists.append(one_sample_list)                    
                
                break
        else:# index does not belong to any list of samples
            partitions.append([index])
            sample_lists.append([])

    return partitions, sample_lists          