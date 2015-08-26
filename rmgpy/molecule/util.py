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