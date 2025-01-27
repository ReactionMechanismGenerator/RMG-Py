from rmgpy.molecule.fragment import Fragment,CuttingLabel
from rmgpy.molecule.molecule import Bond
from rdkit import Chem
from numpy.random import randint
from rmgpy.tools.canteramodel import Cantera
from rmgpy.chemkin import load_chemkin_file
import re
import os
import numpy as np
import matplotlib.pyplot as plt


def match_sequences(seq1, seq2, rtol=1e-6):
    '''
    Given two lists (each item is int or float):
    seq1 and seq2 with same sum, the method returns
    matched indices and values.
    Example:
    seq1 = [1, 3, 1]
    seq2 = [2, 1, 2]
    return: [[(0,0),1],
                    [(1,0),1],
                    [(1,1),1],
                    [(1,2),1],
                    [(2,2),1]]
    '''
    sum_diff = sum(seq2) - sum(seq1)
    assert (
        np.isclose(sum(seq1), sum(seq2), rtol=rtol)
    ), "seq1 has different sum (diff={0}) than seq2.".format(sum_diff)

    # force the sum to be same if the difference
    # is small enough
    if sum_diff >= 0:
        seq1[-1] = seq1[-1] + sum_diff
    else:
        seq2[-1] = seq2[-1] - sum_diff

    # make cumulative sequences
    cum_seq1 = np.cumsum(seq1)
    cum_seq2 = np.cumsum(seq2)

    # add index tags two both cumulative seqs
    pin1 = 0
    pin2 = 0
    matched_indices = []
    matched_cum_values = []
    while pin1 < len(cum_seq1) and pin2 < len(cum_seq2):
        matched_indices.append((pin1, pin2))

        if cum_seq1[pin1] > cum_seq2[pin2]:
            matched_cum_values.append(cum_seq2[pin2])
            pin2 += 1
        elif cum_seq1[pin1] < cum_seq2[pin2]:
            matched_cum_values.append(cum_seq1[pin1])
            pin1 += 1
        else:
            matched_cum_values.append(cum_seq2[pin2])
            pin1 += 1
            pin2 += 1

    # get matches
    matches = []
    for i in range(len(matched_indices)):
        matched_index_tup = matched_indices[i]
        matched_cum_value = matched_cum_values[i]
        if i == 0:
            previous_cum_value = 0
        else:
            previous_cum_value = matched_cum_values[i - 1]

        matches.append(
            [matched_index_tup, matched_cum_value - previous_cum_value])

    return matches


def match_concentrations_with_same_sums(conc1, conc2, rtol=1e-6):
    '''match_concentrations_with_same_sums
    Given two lists with each item to be a tuple
    (species label, concentration)
    conc1 and conc2 with same total concentrations,
    the method returns matched species labels and
    concentrations.
    Example:
    conc1 = [('a', 1),
                    ('b', 3),
                    ('c', 1)]
    conc2 = [('x', 2),
                    ('y', 1),
                    ('z', 2)]
    return: [(('a','x'),1),
                    (('b','x'),1),
                    (('b','y'),1),
                    (('b','z'),1),
                    (('c','z'),1)]
    '''
    labels1 = [tup[0] for tup in conc1]
    labels2 = [tup[0] for tup in conc2]

    seq1 = [tup[1] for tup in conc1]
    seq2 = [tup[1] for tup in conc2]

    matches_seq = FragList.match_sequences(seq1, seq2, rtol)

    matches_conc = []
    for match_seq in matches_seq:
        matched_label_index1 = match_seq[0][0]
        matched_label_index2 = match_seq[0][1]
        matched_value = match_seq[1]

        matched_label1 = labels1[matched_label_index1]
        matched_label2 = labels2[matched_label_index2]
        match_conc = ((matched_label1, matched_label2), matched_value)
        matches_conc.append(match_conc)
    return matches_conc


def match_concentrations_with_different_sums(conc1, conc2):
    """
    Given two lists with each item to be a tuple
    (species label, concentration)
    conc1 and conc2 with different total concentrations,
    the method returns matched species labels and
    concentrations.
    Example:
    conc1 = [('a', 1),
                    ('b', 3),
                    ('c', 1)]
    conc2 = [('x', 2),
                    ('y', 1),
                    ('z', 10)]
    return: [(('a','x', 'z', 'z'),1),
                    (('b','x', 'z', 'z'),1),
                    (('b','y', 'z', 'z'),1),
                    (('b','z', 'z'),1),
                    (('c','z', 'z'),1)]
    """
    labels1 = [tup[0] for tup in conc1]
    labels2 = [tup[0] for tup in conc2]

    seq1 = [tup[1] for tup in conc1]
    seq2 = [tup[1] for tup in conc2]

    matches_conc = []
    pin1 = 0
    pin2 = 0
    val1 = seq1[pin1]
    val2 = seq2[pin2]

    while True:
        if val1 > val2:
            match = ((labels1[pin1], labels2[pin2]), val2)
            matches_conc.append(match)
            val1 = val1 - val2
            pin2 += 1
            if pin2 == len(seq2):
                break
            val2 = seq2[pin2]
        elif val1 < val2:
            match = ((labels1[pin1], labels2[pin2]), val1)
            matches_conc.append(match)
            val2 = val2 - val1
            pin1 += 1
            if pin1 == len(seq1):
                break
            val1 = seq1[pin1]
        else:
            match = ((labels1[pin1], labels2[pin2]), val1)
            matches_conc.append(match)
            pin1 += 1
            pin2 += 1
            if pin1 == len(seq1):
                break
            val1 = seq1[pin1]
            if pin2 == len(seq2):
                break
            val2 = seq2[pin2]

    # if pin2 first reaches the end
    # append all the remaining seq1 to matches_conc
    if pin2 == len(seq2) and pin1 < len(seq1):
        remain_conc1 = [(labels1[pin1], val1)] + conc1[(pin1 + 1):]
        matches_conc.extend(remain_conc1)

    # if pin1 first reaches the end
    # let matches_conc match with remaining seq2
    elif pin1 == len(seq1) and pin2 < len(seq2):
        remain_conc2 = [(labels2[pin2], val2)] + conc2[(pin2 + 1):]
        matches_conc = FragList.match_concentrations_with_different_sums(
            matches_conc, remain_conc2
        )

    # if pin1 and pin2 reach the ends at same time
    # matches_conc is ready to return
    return matches_conc


def shuffle(conc, seed=None):
    """
    Randomly shuffle a list of fragments
    """
    idx_arr = np.arange(len(conc))

    if seed is not None:
        np.random.seed(seed)
    np.random.shuffle(idx_arr)

    return [conc[idx] for idx in idx_arr]


def flatten(combo):
    """
    Given a combo nested `tuple`, e.g.,
    ((('LY', 'XR'), ('LWL', 'RUR'))
    return a list of labels contained in
    the combo ['LY', 'XR', 'LWL', 'RUR']
    """
    return_list = []
    for i in combo:
        if isinstance(i, tuple):
            return_list.extend(FragList.flatten(i))
        else:
            return_list.append(i)
    return return_list


# label should match the desired merging l/'abel on frag2
def merge_frag_to_frag(frag1, frag2, label):
    from rmgpy.molecule import Bond
    from rmgpy.molecule.fragment import Fragment, CuttingLabel

    frag_spe1 = Fragment().from_smiles_like_string(frag1)
    frag_spe2 = Fragment().from_smiles_like_string(frag2)
    # find position of desired CuttingLabel
    # need to find CuttingLabel on frag2 first
    for vertex in frag_spe2.vertices:
        if isinstance(vertex, CuttingLabel):
            if vertex.symbol == label:
                cut2 = vertex

                atom2 = list(cut2.edges.keys())[0]
                frag_spe2.remove_atom(cut2)
                break

    if cut2.symbol[0] == 'L':
        Ctl = cut2.symbol.replace('L', 'R')
    else:  # that means this CuttingLabel is R something

        Ctl = cut2.symbol.replace('R', 'L')

    # merge to frag_spe1
    for vertex in frag_spe1.vertices:
        if isinstance(vertex, CuttingLabel):
            if vertex.symbol == Ctl:
                cut1 = vertex
                atom1 = list(cut1.edges.keys())[0]
                frag_spe1.remove_atom(cut1)
                break

    # new merged fragment
    new_frag = frag_spe1.merge(frag_spe2)
    new_frag.add_bond(Bond(atom1=atom1, atom2=atom2, order=1))
    new_frag = new_frag.copy(deep=True)
    new_frag.update()
    return new_frag  # return Fragment obtl


def merge_frag_list(to_be_merged):
    import os
    # merges fragments in list from right to left
    species_list = []
    ethylene = []
    newlist = []
    warnings = []

    while len(to_be_merged) > 1:

        # second to last fragmentin list
        frag1 = to_be_merged[-2].smiles
        frag2 = to_be_merged[-1].smiles  # last fragment in list

        if 'R' in frag1 and 'L' in frag2:
            newfrag = FragList.merge_frag_to_frag(frag1, frag2, 'L')

        elif 'L' in frag1 and 'R' in frag2:
            newfrag = FragList.merge_frag_to_frag(frag1, frag2, 'R')

        # warn user if last two fragments in list cannot be merged (no R/L
        # combo to be made)
        else:
            print('Warning! Could not merge fragments {} and {}'.format(
                frag1, frag2))

            if 'L' in frag1 and 'L' in frag2:
                newfrag = FragList.merge_frag_to_frag(
                    frag1.replace('L', 'R'), frag2, 'L')
        if len(to_be_merged) > 2:
            cut = len(to_be_merged) - 2
            newfraglist = to_be_merged[:cut]

            newfraglist.append(newfrag)
        elif len(to_be_merged) == 2:

            newfraglist = [newfrag]

        to_be_merged = newfraglist

    to_be_merged = newfraglist
    # newlist.append(newfraglist) # if done merging list, write final
    # structure to list of smiles structures

    # print('{}% of fragments fully merged...'.format(np.round(100*(i+1)/len(flattened_matches_random)),1))
# print(newfraglist)
    return newfraglist

import re
def check_if_radical_near_cutting_label(species_smiles):
    species_fragment = Fragment().from_smiles_like_string(species_smiles)
    for i, atom in enumerate(species_fragment.atoms):
        species_fragment.atoms[i].id = i
    for atom in species_fragment.atoms:
        if atom.radical_electrons ==1:
            bonded_atoms = list(atom.bonds.keys())
            bonded_atom_types = [type(x) for x in bonded_atoms]
            if CuttingLabel in bonded_atom_types:
                radical_near_cutting_label = atom
                nearest_cutting_label = bonded_atoms[bonded_atom_types.index(CuttingLabel)]
                return species_fragment, nearest_cutting_label
            for atom2 in bonded_atoms:
                bonded_atoms_2 = list(atom2.bonds.keys())
                bonded_atom_types_2 = [type(x) for x in bonded_atoms_2]
                if CuttingLabel in bonded_atom_types_2:
                    radical_near_cutting_label = atom2
                    nearest_cutting_label = bonded_atoms_2[bonded_atom_types_2.index(CuttingLabel)]
                    return species_fragment, nearest_cutting_label
    else:
        return False

def check_if_pibond_near_cutting_label(species_smiles):
    species_fragment = Fragment().from_smiles_like_string(species_smiles)
    for i, atom in enumerate(species_fragment.atoms):
        species_fragment.atoms[i].id = i
    for atom in species_fragment.atoms:
        neighboring_bond_orders = [x.order for x in atom.bonds.values()]
        if 2 in neighboring_bond_orders:
            bonded_atoms = list(atom.bonds.keys())
            bonded_atom_types = [type(x) for x in bonded_atoms]
            if CuttingLabel in bonded_atom_types:
                nearest_cutting_label = bonded_atoms[bonded_atom_types.index(CuttingLabel)]
                return species_fragment, nearest_cutting_label
            for atom2 in bonded_atoms:
                bonded_atoms_2 = list(atom2.bonds.keys())
                bonded_atom_types_2 = [type(x) for x in bonded_atoms_2]
                if CuttingLabel in bonded_atom_types_2:
                    nearest_cutting_label = bonded_atoms_2[bonded_atom_types_2.index(CuttingLabel)]
                    return species_fragment, nearest_cutting_label
    else:
        return False

def merge_fragment_a_to_cutting_label_on_b(smiles_a,smiles_b,cuttinglabel):
    a_frag = Fragment().from_smiles_like_string(smiles_a)
    b_frag = Fragment().from_smiles_like_string(smiles_b)
    for vertex in b_frag.vertices:
        if isinstance(vertex, CuttingLabel):

            if vertex.symbol == cuttinglabel.symbol and str(vertex.bonds) == str(cuttinglabel.bonds):
                cutb = vertex

                atom2 = list(cutb.edges.keys())[0]
                b_frag.remove_atom(cutb)
                break
    if cutb.symbol[0] == 'L':
        Ctl = cutb.symbol.replace('L', 'R')
    else:  # that means this CuttingLabel is R something
        Ctl = cutb.symbol.replace('R', 'L')

    # merge to frag_spe1
    for vertex in a_frag.vertices:
        if isinstance(vertex, CuttingLabel):
            if vertex.symbol == Ctl:
                cuta = vertex
                atom1 = list(cuta.edges.keys())[0]
                a_frag.remove_atom(cuta)
                break

    # new merged fragment
    new_frag = a_frag.merge(b_frag)
    new_frag.add_bond(Bond(atom1=atom1, atom2=atom2, order=1))
    new_frag = new_frag.copy(deep=True)
    new_frag.update()
    for i, atom in enumerate(new_frag.atoms):
        new_frag.atoms[i].id = i
    return new_frag  # return Fragment obtl

def get_single_cc_bonds(frag):
    single_cc_bonds = [tuple(sorted([x.atom1.id, x.atom2.id])) for x in frag.get_all_edges() if x.atom1.symbol =="C" and x.atom2.symbol=="C" and x.order ==1]
    return frag, single_cc_bonds

def get_atoms_neighboring_bond(frag, bond_idx):
    atom1id,atom2id = bond_idx
    atom1 = frag.atoms[atom1id]
    atom2 = frag.atoms[atom2id]
    neighboring_atoms = [x for x in list(set(list(atom1.bonds.keys())+list(atom2.bonds.keys()))) if x!=atom1 and x!=atom2]
    return frag, neighboring_atoms

def flatten_list_of_lists(list_of_lists):
    return [item for lst in list_of_lists for item in lst]
    
def get_atoms_nearby_bond(frag, bond_idx):
    atom1id,atom2id = bond_idx
    atom1 = frag.atoms[atom1id]
    atom2 = frag.atoms[atom2id]
    frag, neighboring_atoms = get_atoms_neighboring_bond(frag, bond_idx)
    neighboring_atoms = [x for x in neighboring_atoms]
    neighboring_atoms_ = []
    for atom in neighboring_atoms:
        neighboring_atoms_ += [x for x in atom.bonds.keys() ]
    
    nearby_atoms = list(set(neighboring_atoms + neighboring_atoms_))
    return frag, nearby_atoms

def cut_specific_bond(frag, cut_bond_idx):

    rdfrag, atom_mapping = frag.to_rdkit_mol(remove_h=False, return_mapping=True,  save_order=True)

    atom_mapping_flipped = {v:k for k,v in atom_mapping.items()}


    rdfrag,atom_mapping_flipped = replace_cutting_labels_with_metals(rdfrag, atom_mapping_flipped)
    atom1 = atom_mapping[frag.atoms[cut_bond_idx[0]]]
    atom2 = atom_mapping[frag.atoms[cut_bond_idx[1]]]

    bond_to_cut = rdfrag.GetBondBetweenAtoms(atom1,atom2)

    frag_list = cut_and_place_cutting_labels(rdfrag, bond_to_cut)
    return frag_list
    
def replace_cutting_labels_with_metals(rdfrag,atom_mapping_flipped):
    for idx,atom in atom_mapping_flipped.items():
        if isinstance(atom, CuttingLabel):
            if atom.symbol == "R":
                rdfrag.GetAtomWithIdx(idx).SetAtomicNum(11) # [Na]
            else:
                rdfrag.GetAtomWithIdx(idx).SetAtomicNum(19)  # [K]

    return rdfrag, atom_mapping_flipped

def cut_and_place_cutting_labels(rdfrag, bond_to_cut):

    newmol = Chem.RWMol(rdfrag)
    new_mol = Chem.FragmentOnBonds(newmol, [bond_to_cut.GetIdx()], dummyLabels=[(0,0)])
    mol_set = Chem.GetMolFrags(new_mol, asMols=True)
    
    if len(mol_set) == 2:
        frag1 = Chem.MolToSmiles(mol_set[0])
        frag2 = Chem.MolToSmiles(mol_set[1])

        frag1_R = frag1.count("Na")
        frag1_L = frag1.count("K")
        frag2_R = frag2.count("Na")
        frag2_L = frag2.count("K")

        if frag1_R > frag2_R and frag1_L <= frag2_L:
            frag1_smi = frag1.replace("*", "L")
            frag2_smi = frag2.replace("*", "R")
        elif frag1_L > frag2_L and frag1_R <= frag2_R:
            frag1_smi = frag1.replace("*", "R")
            frag2_smi = frag2.replace("*", "L")
        elif frag2_L > frag1_L and frag2_R <= frag1_R:
            frag1_smi = frag1.replace("*", "L")
            frag2_smi = frag2.replace("*", "R")
        elif frag2_R > frag1_R and frag2_L <= frag1_L:
            frag1_smi = frag1.replace("*", "R")
            frag2_smi = frag2.replace("*", "L")
        elif randint(0,1)==1:
            frag1_smi = frag1.replace("*", "L")
            frag2_smi = frag2.replace("*", "R")
        else:
            frag1_smi = frag1.replace("*", "R")
            frag2_smi = frag2.replace("*", "L")
        frag1_smi = frag1_smi.replace("[Na]", "R")
        frag1_smi = frag1_smi.replace("[K]", "L")
        frag2_smi = frag2_smi.replace("[Na]","R")
        frag2_smi = frag2_smi.replace("[K]","L")
        frag_list = [frag1_smi,frag2_smi]
        frag_list = [Fragment().from_smiles_like_string(x).smiles for x in frag_list]

        return frag_list

def add_starting_fragment_cut_if_needed(species_fragment, nearest_cutting_label, starting_fragment_smiles,species_cutting_threshold = 20):
    merged_frag = merge_fragment_a_to_cutting_label_on_b(starting_fragment_smiles,species_fragment.smiles, nearest_cutting_label)
    merged_frag_smiles = merged_frag.smiles
    if merged_frag_smiles.count("C") + merged_frag_smiles.count("c") > species_cutting_threshold:
        cuttable_bonds = return_cuttable_bonds(merged_frag)
        options = []
        for cuttable_bond_idx_tuple in cuttable_bonds:
            frag1smi, frag2smi = cut_specific_bond(merged_frag, cuttable_bond_idx_tuple)
            options.append([frag1smi,frag2smi])
        
        return options
    else:
        return [[merged_frag_smiles]]
    
def check_if_bond_is_cuttable(frag, bond_idx_tuple):
    frag, nearby_atoms = get_atoms_nearby_bond(frag, bond_idx_tuple)
    
    nearby_bond_orders = flatten_list_of_lists([[x.order for x in atom.bonds.values()] for atom in nearby_atoms])
    nearby_cl_bool = (sum([1 for x in nearby_atoms if type(x)==CuttingLabel])>0)
    nearby_radical_bool = (sum([x.radical_electrons for x in nearby_atoms]) > 0)
    nearby_pi_bond_bool = any(x>1 for x in nearby_bond_orders)
    if nearby_radical_bool ==False and nearby_pi_bond_bool == False and nearby_cl_bool ==False:
        return True
    else:
        return False

def return_cuttable_bonds(frag):
    frag, single_cc_bonds = get_single_cc_bonds(frag)
    
    cuttable_bonds = []
    for bond_idx_tuple in single_cc_bonds:
        cuttable_bool = check_if_bond_is_cuttable(frag, bond_idx_tuple)
        if cuttable_bool:
            cuttable_bonds.append(bond_idx_tuple)
    cuttable_bonds = list(set(cuttable_bonds))
    return cuttable_bonds

def process_new_fragment(species_smiles, starting_fragment_smiles,species_cutting_threshold = 20):

    radical_result = check_if_radical_near_cutting_label(species_smiles)
    

    if radical_result != False:
        species_fragment, nearest_cutting_label = radical_result
        options = add_starting_fragment_cut_if_needed(species_fragment, nearest_cutting_label, starting_fragment_smiles,species_cutting_threshold = species_cutting_threshold)
        return options
    else:
        pibond_result = check_if_pibond_near_cutting_label(species_smiles)
        if pibond_result!=False:
            species_fragment, nearest_cutting_label = pibond_result
            options = add_starting_fragment_cut_if_needed(species_fragment, nearest_cutting_label, starting_fragment_smiles,species_cutting_threshold = species_cutting_threshold)
            return options
        else:
            return species_smiles

def make_new_reaction_string(species_smiles, starting_fragment_smiles, frag_list):
    frag_list_str = " + ".join(frag_list)
    return f"{species_smiles} + {starting_fragment_smiles} => {frag_list_str}"
