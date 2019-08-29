import os

from rmgpy.kinetics import Arrhenius
from rmgpy.chemkin import loadChemkinFile

from afm.fragment import Fragment
from afm.reaction import FragmentReaction

def load_fragment_reactions_from_chemkin(chemkin_path,
                                        dictionary_path, 
                                        fragment_smiles_path):
    """
    This method loads chemkin mechanism and 
    generate fragment reactions in irreversible
    format.
    """
    speciesList, reactionList = loadChemkinFile(chemkin_path, dictionary_path)
    
    # load fragment from smiles-like string
    fragments = []
    with open(fragment_smiles_path) as f_in:
        for line in f_in:
            if line.strip() and not line.startswith('#') and ':' in line:
                label, smiles = [token.strip() for token in line.split(":")]
                frag = Fragment(label=label).from_SMILES_like_string(smiles)
                frag.assign_representative_species()
                frag.species_repr.label = label
                for prev_frag in fragments:
                    if frag.isIsomorphic(prev_frag):
                        raise Exception('Isomorphic duplicate found: {0} and {1}'.format(label, prev_frag.label))
                fragments.append(frag)

    # construct label-key fragment dictionary
    fragments_dict = {}
    for frag0 in fragments:
        if frag0.label not in fragments_dict:
            fragments_dict[frag0.label] = frag0
        else:
            raise Exception('Fragment with duplicated labels found: {0}'.format(frag0.label))

    orig_fragrxns = []
    for rxn0 in reactionList:
        fragrxts = [fragments_dict[spec.label] for spec in rxn0.reactants]
        
        fragprds = [fragments_dict[spec.label] for spec in rxn0.products]
        
        fragpairs = [(fragments_dict[rxt.label], fragments_dict[prod.label]) for rxt, prod in rxn0.pairs]
        fragrxn = FragmentReaction(index=-1,
                                    reactants=fragrxts,
                                    products=fragprds,
                                    kinetics=rxn0.kinetics,
                                   reversible=False,
                                    pairs=fragpairs,
                                    family=rxn0.family)
        orig_fragrxns.append(fragrxn)


    revs_fragrxns = []
    for rxn0 in reactionList:
        fragrxts = [fragments_dict[spec.label] for spec in rxn0.products]
        
        fragprds = [fragments_dict[spec.label] for spec in rxn0.reactants]
        
        fragpairs = [(fragments_dict[prod.label], fragments_dict[rxt.label]) for rxt, prod in rxn0.pairs]
        
        revs_kinetics = rxn0.generateReverseRateCoefficient()
        fragrxn = FragmentReaction(index=-1,
                                    reactants=fragrxts,
                                    products=fragprds,
                                    kinetics=revs_kinetics,
                                   reversible=False,
                                    pairs=fragpairs,
                                    family=rxn0.family)
        revs_fragrxns.append(fragrxn)

    return fragments_dict, orig_fragrxns + revs_fragrxns

def load_pseudo_fragment_reactions(fragments_dict):
    """
    Currently only returns a pseudo reaction. It can be
    extended to generate multiple ractions in the future.
    """

    pseudo_fragrxts = ['RC*C__C', 'RCCCCR']
    pseudo_fragprds = ['RCCCCC__CC*']
    pseudo_frag_pairs = [('RC*C__C', 'RCCCCC__CC*'), ('RCCCCR', 'RCCCCC__CC*')]

    fragrxts = [fragments_dict[label] for label in pseudo_fragrxts]
        
    fragprds = [fragments_dict[label] for label in pseudo_fragprds]

    fragpairs = [(fragments_dict[rxt_label], fragments_dict[prod_label]) for rxt_label, prod_label in pseudo_frag_pairs]

    pseudo_kinetics = Arrhenius(A=(2.000e+05, 'cm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kcal/mol'), T0=(1, 'K'))

    pseudo_fragrxn = FragmentReaction(index=-1,
                                reactants=fragrxts,
                                products=fragprds,
                                kinetics=pseudo_kinetics,
                                reversible=False,
                                pairs=fragpairs,
                                family='pseudo_rxn')

    return [pseudo_fragrxn]