#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
Contains classes for automatically breaking down unstable or non-existent species
"""
import logging
from rmgpy.molecule.group import Group
from rmgpy.species import Species
from rmgpy.data.kinetics.family import ReactionRecipe

decay_group_recipe_pairs = [(Group().from_adjacency_list("""
    1  *3 O u0 p2 c0 {2,S} {4,S}
    2  *2 O u0 p2 c0 {1,S} {3,S}
    3  *1 R!H u1 px c0 {2,S}
    4  H u0 p0  c0 {1,S}
    """),
    ReactionRecipe(actions=[
    ['BREAK_BOND', '*3', 1, '*2'],
    ['CHANGE_BOND', '*2', 1, '*1'],
    ['LOSE_RADICAL', '*1','1'],
    ['GAIN_RADICAL','*3','1']
    ])),
                            ]

def decay_species(spc, check_deterministic=True):
    """
    recursively decays a species object as long as there are valid decay recipes
    raises a warning if the decays could've been done in a different order
    """
    decay = None
    for mol in spc.molecule:
        d = decay_molecule(mol,check_deterministic=check_deterministic)
        if d != mol and decay is None:
            decay = d
            if not check_deterministic:
                break
        elif d != mol:
            sm = spc.molecule[0].to_smiles()
            logging.warning(f"found more than one possible decay for {sm}, may not be able to deterministically determine decay")
            break
    
    if type(decay) == list:
        dlist = []
        for dmol in decay:
            dspc = Species(molecule=[dmol])
            dspc.generate_resonance_structures()
            decay_spc = decay_species(dspc)
            dlist.extend(decay_spc)
        return dlist
    else:
        return [spc]

def decay_molecule(mol, check_deterministic=True):
    """
    breaks down a molecule object down one step using the first valid recipe and atom mappings
    raises a warning if there were multiple valid recipes/atom mappings
    """
    ind = None
    for i,(grp,_) in enumerate(decay_group_recipe_pairs):
        if mol.is_subgraph_isomorphic(grp):  #check is True at most twice
            if ind is None:
                ind = i
                if not check_deterministic:
                    break
            else:
                sm = mol.to_smiles()
                logging.warning(f"found more than one possible decay for {sm}, may not be able to deterministically determine decay")
                break
    
    if ind is None:
        return mol
    else:
        grp,recipe = decay_group_recipe_pairs[ind]
        mol2 = mol.copy(deep=True)
        subs = mol2.find_subgraph_isomorphisms(grp)
        
        if len(subs)>1:
            sm = mol.to_smiles()
            logging.warning(f"found more than one possible decay for {sm}, may not be able to deterministically determine decay")
        sub = subs[0]
        for key,item in sub.items():
            if item.label:
                key.label = item.label
        recipe.apply_forward(mol2)
        mol2.clear_labeled_atoms()
        mol2.update()
        return mol2.split()
        