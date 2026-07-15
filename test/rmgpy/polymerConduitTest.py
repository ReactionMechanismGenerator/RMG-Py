"""Tests for the M18.2 census-only refusal classifier
(rmgpy/polymer_conduit.py, moment_credit_conduit/1 vocabulary).

Ported verbatim from the M18.2 working-set test_classifier.py (41 tests: 14
real poly_102_conduit3 census rows embedded as fixtures + synthetic edge
regimes), plus the round-36 landing-P1 additions at the bottom of this file
(candidate-key stability, both-census overlap precedence, append-only
annotation shape, label/isomorphism divergence path).

Real-row fixtures are embedded verbatim from the poly_102_conduit3 refusal
census (extract.py output over /home/alon/runs/RMG/poly_102_conduit3/RMG.log,
read-only forensic artifact). Synthetic fixtures are clearly marked and cover
gates the real data never exercises (reverse orientation, irreversible source
opposing the admitted direction, near-threshold gas products, unresolved
species) -- the real census has NO rows in those regimes (see census_report.md).
"""

import copy

import pytest

from rmgpy.polymer_conduit import (
    CHAIN_SCALE_HEAVY,
    CHAIN_SCALE_MW,
    GAS_MW_THRESHOLD,
    MONOMER_MW,
    REASON_DIRECTION,
    REASON_FEATURE_RADICAL,
    REASON_GAS_MW_FAIL,
    REASON_SHAPE_DEFERRED,
    REASON_UNKNOWN_SHAPE,
    REASON_UNRESOLVED,
    bucket_counts,
    classify_all,
    classify_record,
    gas_mw_ok,
    is_chain_scale,
    is_pool,
    is_pool_state_resolvable,
    near_threshold,
    species_role,
)

# ---------------------------------------------------------------------------
# REAL fixtures (poly_102_conduit3 census rows, embedded verbatim)
# ---------------------------------------------------------------------------

A_vinylcresol = {'census': 'r93_general',
 'log_reason': 'refused conduit-deferred: chain-scale-discrete/pool coupling (general). A '
               'chain-scale proxy-derived discrete not resolvable to pool state co-occurs with '
               'a polymer pool participant; the whole-row flux is zeroed (stamp-but-keep) '
               'pending the moment-credit conduit.',
 'products': [{'formula': 'C9H10O1',
               'heavy_atoms': 10,
               'index': 6914,
               'label': 'C=Cc1ccc(C)cc1O',
               'link_marker': False,
               'mw': 134.178,
               'source': 'seed_dictionary',
               'token': 'C=Cc1ccc(C)cc1O(6914)'},
              {'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'}],
 'reactants': [{'formula': 'C36H42O4',
                'heavy_atoms': 40,
                'index': 15102,
                'label': 'C=CC1=C(O)C(C2(C)[CH]C=C(CCc3c(C)ccc(C)c3O)C(O)=C2CCc2ccc(C)c(C)c2O)C(C)=C[CH]1',
                'link_marker': False,
                'mw': 538.728,
                'source': 'rdkit_smiles',
                'token': 'C=CC1=C(O)C(C2(C)[CH]C=C(CCc3c(C)ccc(C)c3O)C(O)=C2CCc2ccc(C)c(C)c2O)C(C)=C[CH]1(15102)'}],
 'reaction': 'C=CC1=C(O)C(C2(C)[CH]C=C(CCc3c(C)ccc(C)c3O)C(O)=C2CCc2ccc(C)c(C)c2O)C(C)=C[CH]1(15102) '
             '<=> C=Cc1ccc(C)cc1O(6914) + phenol_formaldehyde(2)',
 'reversible': True}

A_CO = {'census': 'r93_general',
 'log_reason': 'refused conduit-deferred: chain-scale-discrete/pool coupling (general). A '
               'chain-scale proxy-derived discrete not resolvable to pool state co-occurs with '
               'a polymer pool participant; the whole-row flux is zeroed (stamp-but-keep) '
               'pending the moment-credit conduit.',
 'products': [{'formula': 'C1O1',
               'heavy_atoms': 2,
               'index': 17,
               'label': 'CO',
               'link_marker': False,
               'mw': 28.01,
               'source': 'core_dictionary',
               'token': 'CO(17)'},
              {'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'}],
 'reactants': [{'formula': 'C28H32O4',
                'heavy_atoms': 32,
                'index': 1972,
                'label': '[CH-]=[O+]C1(CCc2c(C)ccc(CCc3c(C)ccc(C)c3O)c2O)C=CC(C)=C(C)C1=O',
                'link_marker': False,
                'mw': 432.56,
                'source': 'rdkit_smiles',
                'token': '[CH-]=[O+]C1(CCc2c(C)ccc(CCc3c(C)ccc(C)c3O)c2O)C=CC(C)=C(C)C1=O(1972)'}],
 'reaction': '[CH-]=[O+]C1(CCc2c(C)ccc(CCc3c(C)ccc(C)c3O)c2O)C=CC(C)=C(C)C1=O(1972) <=> CO(17) '
             '+ phenol_formaldehyde(2)',
 'reversible': True}

A_vinylxylenol = {'census': 'r93_general',
 'log_reason': 'refused conduit-deferred: chain-scale-discrete/pool coupling (general). A '
               'chain-scale proxy-derived discrete not resolvable to pool state co-occurs with '
               'a polymer pool participant; the whole-row flux is zeroed (stamp-but-keep) '
               'pending the moment-credit conduit.',
 'products': [{'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'},
              {'formula': 'C10H12O1',
               'heavy_atoms': 11,
               'index': 24,
               'label': 'vinylxylenol_a',
               'link_marker': False,
               'mw': 148.205,
               'source': 'seed_dictionary',
               'token': 'vinylxylenol_a(24)'}],
 'reactants': [{'formula': 'C37H44O4',
                'heavy_atoms': 41,
                'index': 2044,
                'label': '[CH2]C=C1C=CC(C)(C2(C)[CH]C=C(CCc3c(C)ccc(C)c3O)C(O)=C2CCc2ccc(C)c(C)c2O)C(C)=C1O',
                'link_marker': False,
                'mw': 552.755,
                'source': 'rdkit_smiles',
                'token': '[CH2]C=C1C=CC(C)(C2(C)[CH]C=C(CCc3c(C)ccc(C)c3O)C(O)=C2CCc2ccc(C)c(C)c2O)C(C)=C1O(2044)'}],
 'reaction': '[CH2]C=C1C=CC(C)(C2(C)[CH]C=C(CCc3c(C)ccc(C)c3O)C(O)=C2CCc2ccc(C)c(C)c2O)C(C)=C1O(2044) '
             '<=> phenol_formaldehyde(2) + vinylxylenol_a(24)',
 'reversible': True}

B_admissible = {'census': 'r93_general',
 'log_reason': 'refused conduit-deferred: chain-scale-discrete/pool coupling (general). A '
               'chain-scale proxy-derived discrete not resolvable to pool state co-occurs with '
               'a polymer pool participant; the whole-row flux is zeroed (stamp-but-keep) '
               'pending the moment-credit conduit.',
 'products': [{'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'},
              {'formula': 'H2',
               'heavy_atoms': 0,
               'index': 6,
               'label': 'H2',
               'link_marker': False,
               'mw': 2.016,
               'source': 'core_dictionary',
               'token': 'H2(6)'}],
 'reactants': [{'formula': 'H1',
                'heavy_atoms': 0,
                'index': 1,
                'label': 'H',
                'link_marker': False,
                'mw': 1.008,
                'source': 'core_dictionary',
                'token': 'H(1)'},
               {'formula': 'C27H33O3',
                'heavy_atoms': 30,
                'index': 78,
                'label': 'C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O',
                'link_marker': False,
                'mw': 405.558,
                'source': 'rdkit_smiles',
                'token': 'C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O(78)'}],
 'reaction': 'H(1) + C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O(78) <=> '
             'phenol_formaldehyde(2) + H2(6)',
 'reversible': True}

A_char_defer = {'census': 'r93_general',
 'log_reason': 'refused conduit-deferred: chain-scale-discrete/pool coupling (general). A '
               'chain-scale proxy-derived discrete not resolvable to pool state co-occurs with '
               'a polymer pool participant; the whole-row flux is zeroed (stamp-but-keep) '
               'pending the moment-credit conduit.',
 'products': [{'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'},
              {'formula': 'C17H18O2',
               'heavy_atoms': 19,
               'index': 31,
               'label': 'char_C17',
               'link_marker': False,
               'mw': 254.329,
               'source': 'core_dictionary',
               'token': 'char_C17(31)'}],
 'reactants': [{'formula': 'C44H50O5',
                'heavy_atoms': 49,
                'index': 11351,
                'label': 'CC1=C[CH]C(C)=C(O)C1(C=Cc1ccc(C)cc1O)C1(C)[CH]C=C(CCc2c(C)ccc(C)c2O)C(O)=C1CCc1ccc(C)c(C)c1O',
                'link_marker': False,
                'mw': 658.879,
                'source': 'rdkit_smiles',
                'token': 'CC1=C[CH]C(C)=C(O)C1(C=Cc1ccc(C)cc1O)C1(C)[CH]C=C(CCc2c(C)ccc(C)c2O)C(O)=C1CCc1ccc(C)c(C)c1O(11351)'}],
 'reaction': 'CC1=C[CH]C(C)=C(O)C1(C=Cc1ccc(C)cc1O)C1(C)[CH]C=C(CCc2c(C)ccc(C)c2O)C(O)=C1CCc1ccc(C)c(C)c1O(11351) '
             '<=> phenol_formaldehyde(2) + char_C17(31)',
 'reversible': True}

A_char224 = {'census': 'r93_general',
 'log_reason': 'refused conduit-deferred: chain-scale-discrete/pool coupling (general). A '
               'chain-scale proxy-derived discrete not resolvable to pool state co-occurs with '
               'a polymer pool participant; the whole-row flux is zeroed (stamp-but-keep) '
               'pending the moment-credit conduit.',
 'products': [{'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'},
              {'formula': 'C16H16O1',
               'heavy_atoms': 17,
               'index': 38,
               'label': 'char_mature17',
               'link_marker': False,
               'mw': 224.303,
               'source': 'seed_dictionary',
               'token': 'char_mature17(38)'}],
 'reactants': [{'formula': 'C43H48O4',
                'heavy_atoms': 47,
                'index': 7545,
                'label': 'CC1=C[CH]C2(C3(C)[CH]C=C(CCc4c(C)ccc(C)c4O)C(O)=C3CCc3ccc(C)c(C)c3O)C(=C1)Cc1cc(C)c(O)c(C)c12',
                'link_marker': False,
                'mw': 628.853,
                'source': 'rdkit_smiles',
                'token': 'CC1=C[CH]C2(C3(C)[CH]C=C(CCc4c(C)ccc(C)c4O)C(O)=C3CCc3ccc(C)c(C)c3O)C(=C1)Cc1cc(C)c(O)c(C)c12(7545)'}],
 'reaction': 'CC1=C[CH]C2(C3(C)[CH]C=C(CCc4c(C)ccc(C)c4O)C(O)=C3CCc3ccc(C)c(C)c3O)C(=C1)Cc1cc(C)c(O)c(C)c12(7545) '
             '<=> phenol_formaldehyde(2) + char_mature17(38)',
 'reversible': True}

B_heavy_defer = {'census': 'r93_general',
 'log_reason': 'refused conduit-deferred: chain-scale-discrete/pool coupling (general). A '
               'chain-scale proxy-derived discrete not resolvable to pool state co-occurs with '
               'a polymer pool participant; the whole-row flux is zeroed (stamp-but-keep) '
               'pending the moment-credit conduit.',
 'products': [{'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'},
              {'formula': 'C17H20O2',
               'heavy_atoms': 19,
               'index': 34,
               'label': 'C17arylA_H',
               'link_marker': False,
               'mw': 256.345,
               'source': 'core_dictionary',
               'token': 'C17arylA_H(34)'}],
 'reactants': [{'formula': 'C17H19O2',
                'heavy_atoms': 19,
                'index': 1675,
                'label': 'Cc1ccc(CCc2c(C)ccc(C)c2[O])c(O)c1',
                'link_marker': False,
                'mw': 255.337,
                'source': 'rdkit_smiles',
                'token': 'Cc1ccc(CCc2c(C)ccc(C)c2[O])c(O)c1(1675)'},
               {'formula': 'C27H33O3',
                'heavy_atoms': 30,
                'index': 78,
                'label': 'C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O',
                'link_marker': False,
                'mw': 405.558,
                'source': 'rdkit_smiles',
                'token': 'C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O(78)'}],
 'reaction': 'Cc1ccc(CCc2c(C)ccc(C)c2[O])c(O)c1(1675) + '
             'C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O(78) <=> '
             'phenol_formaldehyde(2) + C17arylA_H(34)',
 'reversible': True}

C_row = {'census': 'r93_general',
 'log_reason': 'refused conduit-deferred: chain-scale-discrete/pool coupling (general). A '
               'chain-scale proxy-derived discrete not resolvable to pool state co-occurs with '
               'a polymer pool participant; the whole-row flux is zeroed (stamp-but-keep) '
               'pending the moment-credit conduit.',
 'products': [{'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'},
              {'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'}],
 'reactants': [{'formula': 'C54H64O6',
                'heavy_atoms': 60,
                'index': 96,
                'label': 'Cc1ccc(CCC2=C(O)C(CCc3c(C)ccc(C)c3O)=C[CH]C2(C)C2(C)[CH]C=C(CCc3c(C)ccc(C)c3O)C(O)=C2CCc2ccc(C)c(C)c2O)c(O)c1C',
                'link_marker': False,
                'mw': 809.1,
                'source': 'rdkit_smiles',
                'token': 'Cc1ccc(CCC2=C(O)C(CCc3c(C)ccc(C)c3O)=C[CH]C2(C)C2(C)[CH]C=C(CCc3c(C)ccc(C)c3O)C(O)=C2CCc2ccc(C)c(C)c2O)c(O)c1C(96)'}],
 'reaction': 'Cc1ccc(CCC2=C(O)C(CCc3c(C)ccc(C)c3O)=C[CH]C2(C)C2(C)[CH]C=C(CCc3c(C)ccc(C)c3O)C(O)=C2CCc2ccc(C)c(C)c2O)c(O)c1C(96) '
             '<=> phenol_formaldehyde(2) + phenol_formaldehyde(2)',
 'reversible': True}

D_row = {'census': 'r93_general',
 'log_reason': 'refused conduit-deferred: chain-scale-discrete/pool coupling (general). A '
               'chain-scale proxy-derived discrete not resolvable to pool state co-occurs with '
               'a polymer pool participant; the whole-row flux is zeroed (stamp-but-keep) '
               'pending the moment-credit conduit.',
 'products': [{'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'},
              {'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'}],
 'reactants': [{'formula': 'C27H31O3',
                'heavy_atoms': 30,
                'index': 42,
                'label': 'Cc1ccc(CCc2c(C)ccc(CCc3c(C)ccc(C)c3O)c2O)c([O])c1C',
                'link_marker': False,
                'mw': 403.542,
                'source': 'rdkit_smiles',
                'token': 'Cc1ccc(CCc2c(C)ccc(CCc3c(C)ccc(C)c3O)c2O)c([O])c1C(42)'},
               {'formula': 'C27H33O3',
                'heavy_atoms': 30,
                'index': 78,
                'label': 'C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O',
                'link_marker': False,
                'mw': 405.558,
                'source': 'rdkit_smiles',
                'token': 'C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O(78)'}],
 'reaction': 'Cc1ccc(CCc2c(C)ccc(CCc3c(C)ccc(C)c3O)c2O)c([O])c1C(42) + '
             'C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O(78) <=> '
             'phenol_formaldehyde(2) + phenol_formaldehyde(2)',
 'reversible': True}

E_row = {'census': 'r93_general',
 'log_reason': 'refused conduit-deferred: chain-scale-discrete/pool coupling (general). A '
               'chain-scale proxy-derived discrete not resolvable to pool state co-occurs with '
               'a polymer pool participant; the whole-row flux is zeroed (stamp-but-keep) '
               'pending the moment-credit conduit.',
 'products': [{'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'},
              {'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'}],
 'reactants': [{'formula': 'C27H31O3',
                'heavy_atoms': 30,
                'index': 22,
                'label': 'trimer_rad33',
                'link_marker': False,
                'mw': 403.542,
                'source': 'core_dictionary',
                'token': 'trimer_rad33(22)'},
               {'formula': 'C27H33O3',
                'heavy_atoms': 30,
                'index': 78,
                'label': 'C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O',
                'link_marker': False,
                'mw': 405.558,
                'source': 'rdkit_smiles',
                'token': 'C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O(78)'}],
 'reaction': 'trimer_rad33(22) + C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O(78) <=> '
             'phenol_formaldehyde(2) + phenol_formaldehyde(2)',
 'reversible': True}

F_row = {'census': 'r93_general',
 'log_reason': 'refused conduit-deferred: chain-scale-discrete/pool coupling (general). A '
               'chain-scale proxy-derived discrete not resolvable to pool state co-occurs with '
               'a polymer pool participant; the whole-row flux is zeroed (stamp-but-keep) '
               'pending the moment-credit conduit.',
 'products': [{'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'}],
 'reactants': [{'formula': 'C27H31O2',
                'heavy_atoms': 29,
                'index': 40,
                'label': 'Cc1[c]c(CCc2c(C)ccc(CCc3c(C)ccc(C)c3O)c2O)ccc1C',
                'link_marker': False,
                'mw': 387.543,
                'source': 'rdkit_smiles',
                'token': 'Cc1[c]c(CCc2c(C)ccc(CCc3c(C)ccc(C)c3O)c2O)ccc1C(40)'},
               {'formula': 'HO',
                'heavy_atoms': 1,
                'index': 41,
                'label': '[OH]',
                'link_marker': False,
                'mw': 17.007,
                'source': 'rdkit_smiles',
                'token': '[OH](41)'}],
 'reaction': 'Cc1[c]c(CCc2c(C)ccc(CCc3c(C)ccc(C)c3O)c2O)ccc1C(40) + [OH](41) <=> '
             'phenol_formaldehyde(2)',
 'reversible': True}

FR_elim = {'census': 'feature_radical',
 'log_reason': 'eliminating radical refused (conduit-deferred); no flux applied '
               "(stamp-but-keep). Deferred to item 20's conduit.",
 'products': [{'formula': 'C27H32O3',
               'heavy_atoms': 30,
               'index': 2,
               'label': 'phenol_formaldehyde',
               'link_marker': False,
               'mw': 404.55,
               'source': 'seed_dictionary',
               'token': 'phenol_formaldehyde(2)'}],
 'reactants': [{'formula': 'H1',
                'heavy_atoms': 0,
                'index': 1,
                'label': 'H',
                'link_marker': False,
                'mw': 1.008,
                'source': 'core_dictionary',
                'token': 'H(1)'},
               {'formula': 'C27H31O3',
                'heavy_atoms': 30,
                'index': 22,
                'label': 'trimer_rad33',
                'link_marker': False,
                'mw': 403.542,
                'source': 'core_dictionary',
                'token': 'trimer_rad33(22)'}],
 'reaction': 'H(1) + trimer_rad33(22) <=> phenol_formaldehyde(2)',
 'reversible': True}

FR_qssa_irrev = {'census': 'feature_radical',
 'log_reason': 'accumulating radical refused (qssa-unassessable); no flux applied '
               "(stamp-but-keep). Deferred to item 20's conduit.",
 'products': [{'formula': 'C9H12O',
               'heavy_atoms': 10,
               'index': 6262,
               'label': 'Cc1ccc(CL)c(O)c1C',
               'link_marker': True,
               'mw': 136.194,
               'source': 'rdkit_smiles_linkmapped',
               'token': 'Cc1ccc(CL)c(O)c1C(6262)'},
              {'formula': 'C25H29O3',
               'heavy_atoms': 28,
               'index': 6261,
               'label': 'CC1=C(CR)[C](O)C(CCc2c(C)ccc(C)c2O)(c2ccc(C)cc2O)C=C1',
               'link_marker': True,
               'mw': 377.504,
               'source': 'rdkit_smiles_linkmapped',
               'token': 'CC1=C(CR)[C](O)C(CCc2c(C)ccc(C)c2O)(c2ccc(C)cc2O)C=C1(6261)'}],
 'reactants': [{'formula': 'C7H7O1',
                'heavy_atoms': 8,
                'index': 36,
                'label': 'C7H7O_cresyl',
                'link_marker': False,
                'mw': 107.132,
                'source': 'core_dictionary',
                'token': 'C7H7O_cresyl(36)'},
               {'formula': 'C27H32O3',
                'heavy_atoms': 30,
                'index': 2,
                'label': 'phenol_formaldehyde',
                'link_marker': False,
                'mw': 404.55,
                'source': 'seed_dictionary',
                'token': 'phenol_formaldehyde(2)'}],
 'reaction': 'C7H7O_cresyl(36) + phenol_formaldehyde(2) => Cc1ccc(CL)c(O)c1C(6262) + '
             'CC1=C(CR)[C](O)C(CCc2c(C)ccc(C)c2O)(c2ccc(C)cc2O)C=C1(6261)',
 'reversible': False}

FR_qssa_inval = {'census': 'feature_radical',
 'log_reason': 'accumulating radical refused (qssa-invalid); no flux applied (stamp-but-keep). '
               "Deferred to item 20's conduit.",
 'products': [{'formula': 'C27H33O3',
               'heavy_atoms': 30,
               'index': 78,
               'label': 'C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O',
               'link_marker': False,
               'mw': 405.558,
               'source': 'rdkit_smiles',
               'token': 'C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O(78)'}],
 'reactants': [{'formula': 'H1',
                'heavy_atoms': 0,
                'index': 1,
                'label': 'H',
                'link_marker': False,
                'mw': 1.008,
                'source': 'core_dictionary',
                'token': 'H(1)'},
               {'formula': 'C27H32O3',
                'heavy_atoms': 30,
                'index': 2,
                'label': 'phenol_formaldehyde',
                'link_marker': False,
                'mw': 404.55,
                'source': 'seed_dictionary',
                'token': 'phenol_formaldehyde(2)'}],
 'reaction': 'H(1) + phenol_formaldehyde(2) <=> '
             'C[C]1C=CC(CCc2c(C)ccc(C)c2O)=C(O)C1CCc1ccc(C)c(C)c1O(78)',
 'reversible': True}

REAL_FIXTURES = {
    "A_vinylcresol": A_vinylcresol,
    "A_CO": A_CO,
    "A_vinylxylenol": A_vinylxylenol,
    "B_admissible": B_admissible,
    "A_char_defer": A_char_defer,
    "A_char224": A_char224,
    "B_heavy_defer": B_heavy_defer,
    "C_row": C_row,
    "D_row": D_row,
    "E_row": E_row,
    "F_row": F_row,
    "FR_elim": FR_elim,
    "FR_qssa_irrev": FR_qssa_irrev,
    "FR_qssa_inval": FR_qssa_inval,
}

EXPECTED_BUCKETS = {
    "A_vinylcresol": "ADMISSIBLE_A",
    "A_CO": "ADMISSIBLE_A",
    "A_vinylxylenol": "ADMISSIBLE_A",
    "B_admissible": "ADMISSIBLE_B",
    "A_char_defer": "DEFERRED_A",
    "A_char224": "DEFERRED_A",
    "B_heavy_defer": "DEFERRED_B",
    "C_row": "DEFERRED_C",
    "D_row": "DEFERRED_D",
    "E_row": "DEFERRED_E",
    "F_row": "DEFERRED_F",
    "FR_elim": "FEATURE_RADICAL",
    "FR_qssa_irrev": "FEATURE_RADICAL",
    "FR_qssa_inval": "FEATURE_RADICAL",
}


# ---------------------------------------------------------------------------
# SYNTHETIC fixtures (marked; regimes absent from the real census)
# ---------------------------------------------------------------------------

def _species(token, mw, heavy, label=None, link_marker=False):
    return {
        "token": token,
        "label": label if label is not None else token.rsplit("(", 1)[0],
        "index": None,
        "formula": None,
        "mw": mw,
        "heavy_atoms": heavy,
        "source": "synthetic",
        "link_marker": link_marker,
    }


SYN_CHAIN = _species("SYN_chain(900001)", 405.5, 30)
SYN_POOL = _species("phenol_formaldehyde(2)", 404.55, 30, label="phenol_formaldehyde")
SYN_GAS_OK = _species("SYN_gas134(900002)", 134.178, 10)
# near-threshold band is +-10% of 201.267 -> [181.14, 221.39]
SYN_GAS_NEAR_UNDER = _species("SYN_gas195(900003)", 195.0, 14)   # pass + near flag
SYN_GAS_NEAR_OVER = _species("SYN_gas205(900004)", 205.0, 15)    # fail + near flag
SYN_GAS_UNRESOLVED = _species("SYN_unknown(900005)", None, None)


def _syn_record(reactants, products, reversible=True, census="r93_general"):
    return {
        "census": census,
        "reaction": " + ".join(s["token"] for s in reactants)
        + (" <=> " if reversible else " => ")
        + " + ".join(s["token"] for s in products),
        "reversible": reversible,
        "reactants": [copy.deepcopy(s) for s in reactants],
        "products": [copy.deepcopy(s) for s in products],
        "log_reason": "synthetic fixture",
    }


# ---------------------------------------------------------------------------
# Species-level tests
# ---------------------------------------------------------------------------

class TestSpeciesRoles:
    def test_monomer_constants(self):
        assert MONOMER_MW == pytest.approx(134.178, abs=1e-3)
        assert GAS_MW_THRESHOLD == pytest.approx(201.267, abs=1e-3)
        assert CHAIN_SCALE_MW == pytest.approx(335.445, abs=1e-3)
        assert CHAIN_SCALE_HEAVY == pytest.approx(25.0)

    def test_pool_by_label(self):
        assert is_pool(SYN_POOL)
        assert species_role(SYN_POOL) == "POOL"
        mod = _species("phenol_formaldehyde_mod_3(75)", 404.55, 30,
                       label="phenol_formaldehyde_mod_3")
        assert is_pool(mod)

    def test_chip_resolves_to_pool_state(self):
        # trimer_rad* chips are chain-SCALE by size but take the POOL role
        # (this is the empirical D/E split in the census).
        chip = E_row["reactants"][0]
        assert chip["label"].startswith("trimer_rad")
        assert is_pool_state_resolvable(chip)
        assert is_chain_scale(chip)  # size alone would say CHAIN...
        assert species_role(chip) == "POOL"  # ...but resolvability wins

    def test_chain_scale_requires_both_axes(self):
        # dimer-scale: heavy mass but < 25 heavy atoms -> DISC
        dimer = _species("SYN_dimer(900006)", 256.3, 19)
        assert not is_chain_scale(dimer)
        assert species_role(dimer) == "DISC"
        # boundary: exactly 2.5 monomer-equivalents on both axes -> CHAIN
        edge = _species("SYN_edge(900007)", CHAIN_SCALE_MW, 25)
        assert is_chain_scale(edge)

    def test_gas_mw_threshold(self):
        assert gas_mw_ok(SYN_GAS_OK)
        assert gas_mw_ok(_species("x(1)", GAS_MW_THRESHOLD, 10))  # <= passes
        assert not gas_mw_ok(_species("x(1)", GAS_MW_THRESHOLD + 0.01, 10))
        assert not gas_mw_ok(SYN_GAS_UNRESOLVED)

    def test_near_threshold_band(self):
        assert near_threshold(SYN_GAS_NEAR_UNDER)
        assert near_threshold(SYN_GAS_NEAR_OVER)
        assert not near_threshold(SYN_GAS_OK)          # 134.2, far below
        assert not near_threshold(_species("x(1)", 224.303, 17))  # char_mature17


# ---------------------------------------------------------------------------
# Real-row bucketing
# ---------------------------------------------------------------------------

class TestRealRowBuckets:
    @pytest.mark.parametrize("name", sorted(REAL_FIXTURES))
    def test_bucket(self, name):
        result = classify_record(REAL_FIXTURES[name])
        assert result["bucket"] == EXPECTED_BUCKETS[name], result

    def test_admissible_a_fields(self):
        r = classify_record(A_vinylcresol)
        assert r["admissible"] is True
        assert r["shape"] == "A"
        assert r["admitted_direction"] == "forward"
        assert r["irreversible_only"] is True
        assert r["destination_pool"] == "phenol_formaldehyde"
        assert r["gas_products"] == ["C=Cc1ccc(C)cc1O(6914)"]  # C9H10O, 134.178
        assert r["gas_products_mw_ok"] is True
        assert r["gas_products_near_threshold"] == []
        assert r["refusal_reason"] is None

    def test_admissible_a_small_gas_co(self):
        r = classify_record(A_CO)
        assert r["bucket"] == "ADMISSIBLE_A"
        assert r["gas_products"] == ["CO(17)"]

    def test_admissible_b_fields(self):
        r = classify_record(B_admissible)
        assert r["admissible"] is True
        assert r["shape"] == "B"
        assert r["destination_pool"] == "phenol_formaldehyde"
        assert r["gas_products"] == ["H2(6)"]
        assert r["gas_reactants_over_threshold"] == []

    def test_deferred_a_char_reason(self):
        r = classify_record(A_char_defer)
        assert r["admissible"] is False
        assert r["refusal_reason"] == REASON_GAS_MW_FAIL
        assert "char_C17(31)" in r["gas_products"]
        assert r["gas_products_mw_ok"] is False
        # char_C17 (254.3) is over threshold but NOT in the +-10% band
        assert r["gas_products_near_threshold"] == []

    def test_deferred_b_heavy_reactant_flag(self):
        # every DEFERRED_B row in the census consumes a >threshold discrete
        # AND produces a >threshold char; the reactant side is a flag only,
        # the product side is the gate.
        r = classify_record(B_heavy_defer)
        assert r["bucket"] == "DEFERRED_B"
        assert r["refusal_reason"] == REASON_GAS_MW_FAIL
        assert r["gas_reactants_over_threshold"]  # informational

    def test_deferred_shapes_reason(self):
        for fx in (C_row, D_row, E_row, F_row):
            r = classify_record(fx)
            assert r["refusal_reason"] == REASON_SHAPE_DEFERRED
            assert r["admissible"] is False

    def test_two_pool_rows_have_no_destination(self):
        for fx in (C_row, D_row, E_row):
            r = classify_record(fx)
            assert r["destination_pool"] is None  # ambiguous, part of deferral

    def test_f_row_gas_association_has_no_gas_products(self):
        r = classify_record(F_row)
        assert r["shape"] == "F"
        assert r["gas_products"] == []

    def test_feature_radical_always_deferred(self):
        for fx in (FR_elim, FR_qssa_irrev, FR_qssa_inval):
            r = classify_record(fx)
            assert r["bucket"] == "FEATURE_RADICAL"
            assert r["refusal_reason"] == REASON_FEATURE_RADICAL
            assert r["admissible"] is False

    def test_feature_radical_link_marker_flag(self):
        # the 8 irreversible qssa rows carry CL/CR proxy link pseudo-atoms
        r = classify_record(FR_qssa_irrev)
        assert "link-marker-species" in r["flags"]
        assert FR_qssa_irrev["reversible"] is False


# ---------------------------------------------------------------------------
# Synthetic edge cases (regimes absent from the real census)
# ---------------------------------------------------------------------------

class TestSyntheticEdgeCases:
    def test_reverse_oriented_row_is_normalized(self):
        # SYNTHETIC: pool + gas on the LEFT, chain on the RIGHT. All 18,827
        # real rows are written chain->pool; the classifier must still admit
        # the flipped orientation, one-way, in the reverse direction.
        rec = _syn_record([SYN_POOL, SYN_GAS_OK], [SYN_CHAIN], reversible=True)
        r = classify_record(rec)
        assert r["shape"] == "A"
        assert r["bucket"] == "ADMISSIBLE_A"
        assert r["admitted_direction"] == "reverse"

    def test_irreversible_source_opposing_admitted_direction(self):
        # SYNTHETIC: an irreversible row written pool -> chain can never run
        # in the chain-consuming direction.
        rec = _syn_record([SYN_POOL, SYN_GAS_OK], [SYN_CHAIN], reversible=False)
        r = classify_record(rec)
        assert r["bucket"] == "DEFERRED_A"
        assert r["refusal_reason"] == REASON_DIRECTION

    def test_irreversible_source_along_admitted_direction_is_fine(self):
        rec = _syn_record([SYN_CHAIN], [SYN_POOL, SYN_GAS_OK], reversible=False)
        r = classify_record(rec)
        assert r["bucket"] == "ADMISSIBLE_A"
        assert r["admitted_direction"] == "forward"

    def test_near_threshold_under_admits_with_flag(self):
        # SYNTHETIC: real census has ZERO gas products in the +-10% band
        # (gap between 168.24 and 224.30).
        rec = _syn_record([SYN_CHAIN], [SYN_POOL, SYN_GAS_NEAR_UNDER])
        r = classify_record(rec)
        assert r["bucket"] == "ADMISSIBLE_A"
        assert r["gas_products_near_threshold"] == [SYN_GAS_NEAR_UNDER["token"]]

    def test_near_threshold_over_defers_with_flag(self):
        rec = _syn_record([SYN_CHAIN], [SYN_POOL, SYN_GAS_NEAR_OVER])
        r = classify_record(rec)
        assert r["bucket"] == "DEFERRED_A"
        assert r["refusal_reason"] == REASON_GAS_MW_FAIL
        assert r["gas_products_near_threshold"] == [SYN_GAS_NEAR_OVER["token"]]

    def test_multiple_gas_products_fail_closed(self):
        # SYNTHETIC: no real r93 row has two gas products (round-32 shape A is
        # "pool + 1 discrete" by definition; A/B rows carry exactly one gas
        # product). A multi-gas row is OUTSIDE the adjudicated vocabulary and
        # must fail CLOSED (UNCLASSIFIED for review), never be silently
        # admitted.
        rec = _syn_record([SYN_CHAIN],
                          [SYN_POOL, SYN_GAS_OK, SYN_GAS_NEAR_OVER])
        r = classify_record(rec)
        assert r["admissible"] is False
        assert r["bucket"] == "UNCLASSIFIED"
        assert r["refusal_reason"] == REASON_UNKNOWN_SHAPE

    def test_unresolved_species_unclassified(self):
        rec = _syn_record([SYN_CHAIN], [SYN_POOL, SYN_GAS_UNRESOLVED])
        r = classify_record(rec)
        assert r["bucket"] == "UNCLASSIFIED"
        assert r["refusal_reason"] == REASON_UNRESOLVED

    def test_gas_association_shape_f(self):
        # SYNTHETIC mirror of the 23 real F rows: chain + small radical -> pool
        rad = _species("SYN_OH(900008)", 17.007, 1)
        rec = _syn_record([SYN_CHAIN, rad], [SYN_POOL])
        r = classify_record(rec)
        assert r["bucket"] == "DEFERRED_F"


# ---------------------------------------------------------------------------
# Purity / determinism
# ---------------------------------------------------------------------------

class TestDeterminism:
    def test_input_not_mutated(self):
        rec = copy.deepcopy(A_vinylcresol)
        classify_record(rec)
        assert rec == A_vinylcresol

    def test_repeat_classification_identical(self):
        fixtures = list(REAL_FIXTURES.values())
        first = classify_all(fixtures)
        second = classify_all(fixtures)
        assert first == second

    def test_bucket_counts_helper(self):
        counts = bucket_counts(classify_all(list(REAL_FIXTURES.values())))
        assert counts["ADMISSIBLE_A"] == 3
        assert counts["ADMISSIBLE_B"] == 1
        assert counts["FEATURE_RADICAL"] == 3


# ---------------------------------------------------------------------------
# Round-36 landing P1s (M18.2 in-repo additions)
# ---------------------------------------------------------------------------

from rmgpy.polymer_conduit import (  # noqa: E402
    candidate_key_from_label,
    census_summary,
    census_suffix,
    conduit_candidate_key,
    gas_mw_threshold_for_pools,
    lookup_candidate,
    register_candidate,
    reset_census_registry,
)


@pytest.fixture(autouse=False)
def clean_registry():
    reset_census_registry()
    yield
    reset_census_registry()


class TestCandidateKeyStability:
    """Round-36 P1(a): the candidate key is a deterministic, orientation-
    and census-independent identity, so M18.3 can detect already-classified
    candidates and never double-book."""

    def test_orientation_independent(self, clean_registry):
        k1 = candidate_key_from_label("A(1) + B(2) <=> C(3) + D(4)")
        k2 = candidate_key_from_label("D(4) + C(3) <=> B(2) + A(1)")
        assert k1 == k2

    def test_arrow_and_census_independent(self, clean_registry):
        # reversibility is a property of the row, not of its identity: the
        # overlap ledger must join r93/feature-radical sightings of one row
        assert (candidate_key_from_label("A(1) + B(2) => C(3)")
                == candidate_key_from_label("A(1) + B(2) <=> C(3)"))

    def test_distinct_rows_distinct_keys(self, clean_registry):
        assert (candidate_key_from_label("A(1) + B(2) <=> C(3)")
                != candidate_key_from_label("A(1) + B(2) <=> C(4)"))

    def test_record_and_label_paths_agree(self, clean_registry):
        # the string-only feature-radical hook and the full r93 record path
        # must produce IDENTICAL keys for the same row
        rec = dict(A_vinylcresol)
        assert conduit_candidate_key(rec) == candidate_key_from_label(
            rec["reaction"])

    def test_key_is_stable_across_reruns(self, clean_registry):
        label = A_vinylcresol["reaction"]
        assert (candidate_key_from_label(label)
                == candidate_key_from_label(label)
                == candidate_key_from_label(" " + label + " "))


class TestOverlapPrecedence:
    """Round-36 P1(a): 8,561 rows sit in BOTH censuses; FEATURE-RADICAL
    refusal WINS as the upstream blocker until admission exists."""

    KEY = "SYNKEY_chain(1)<>SYNKEY_pool(2)+SYNKEY_gas(3)"

    def test_feature_radical_wins_when_seen_second(self, clean_registry):
        register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A")
        entry = register_candidate(self.KEY, "feature_radical",
                                   "FEATURE_RADICAL")
        assert entry["effective_bucket"] == "FEATURE_RADICAL"
        assert entry["shadow_bucket"] == "ADMISSIBLE_A"
        assert entry["precedence"] == "feature_radical"
        assert entry["censuses"] == {"r93_general", "feature_radical"}

    def test_feature_radical_wins_when_seen_first(self, clean_registry):
        register_candidate(self.KEY, "feature_radical", "FEATURE_RADICAL")
        entry = register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A")
        assert entry["effective_bucket"] == "FEATURE_RADICAL"
        assert entry["shadow_bucket"] == "ADMISSIBLE_A"
        assert entry["precedence"] == "feature_radical"

    def test_single_census_key_has_no_precedence(self, clean_registry):
        entry = register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A")
        assert entry["effective_bucket"] == "ADMISSIBLE_A"
        assert entry["shadow_bucket"] is None
        assert entry["precedence"] is None

    def test_lookup_is_the_m183_no_double_book_gate(self, clean_registry):
        assert lookup_candidate(self.KEY) is None
        register_candidate(self.KEY, "feature_radical", "FEATURE_RADICAL")
        register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A")
        entry = lookup_candidate(self.KEY)
        # ONE ledger entry per key; M18.3 books effective_bucket only.
        assert entry["effective_bucket"] == "FEATURE_RADICAL"
        assert entry["bucket_by_census"] == {
            "feature_radical": "FEATURE_RADICAL",
            "r93_general": "ADMISSIBLE_A"}

    def test_summary_always_names_unclassified(self, clean_registry):
        # LOUD census (round-36 P1(c)): UNCLASSIFIED surfaced even at zero.
        register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A")
        s = census_summary()
        assert "UNCLASSIFIED=0" in s
        assert "candidates=1" in s and "overlap=0" in s


class TestAppendOnlyAnnotation:
    """Round-36 P1(b): warning edits are APPEND-ONLY -- the historical
    header/reason text is byte-identical and the structured annotation is a
    suffix. Integration-level: drive the REAL r93 refusal branch in
    rmgpy/polymer.py with live RMG objects and inspect the census line."""

    _R93_PREFIX = "POLYMER REFUSAL CENSUS (r93 general branch): "
    _R93_BODY = (
        "-- refused conduit-deferred: chain-scale-discrete/pool coupling "
        "(general). A chain-scale proxy-derived discrete not resolvable to "
        "pool state co-occurs with a polymer pool participant; the "
        "whole-row flux is zeroed (stamp-but-keep) pending the "
        "moment-credit conduit.")

    def _refused_row(self):
        """Minimal live fixture mirroring polymerTest's r93 run-5 shape-D
        cohort: two PP-scale pools + a gas-vetoed chain-scale discrete."""
        from rmgpy.molecule import Molecule
        from rmgpy.polymer import Polymer, set_polymer_gas_veto
        from rmgpy.reaction import Reaction
        from rmgpy.species import Species

        fr1 = Polymer(label="FR1", monomer="[CH2][CH]C",
                      Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        sidegrp = Polymer(label="FR1_sidegrp", monomer="[CH2][CH]C",
                          Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        discrete = Species(molecule=[Molecule().from_smiles(
            "[CH2]C(CBr)CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)Br")])
        discrete.label = "(5)"
        set_polymer_gas_veto(discrete)
        return Reaction(reactants=[fr1, sidegrp],
                        products=[discrete, fr1], reversible=True)

    def test_r93_census_line_gains_suffix_append_only(self, caplog,
                                                      clean_registry):
        import logging as _logging

        from rmgpy.polymer import (_general_chain_scale_pool_warned,
                                   stamp_gas_association_refusal)

        rxn = self._refused_row()
        _general_chain_scale_pool_warned.clear()
        with caplog.at_level(_logging.WARNING):
            stamp_gas_association_refusal(rxn)
        assert rxn.polymer_refused is True          # behavior unchanged
        assert rxn.polymer_refused_accumulating is False
        lines = [r.getMessage() for r in caplog.records
                 if r.getMessage().startswith(self._R93_PREFIX)]
        assert len(lines) == 1
        msg = lines[0]
        # append-only: the historical body text is intact, verbatim...
        assert self._R93_BODY in msg
        # ... and the conduit annotation is a SUFFIX after it.
        body_end = msg.index(self._R93_BODY) + len(self._R93_BODY)
        suffix = msg[body_end:]
        assert suffix.startswith(" [conduit-census/1 key=")
        assert "bucket=" in suffix and "censuses=r93_general" in suffix
        assert "unclassified=" in suffix          # loud even when zero
        assert msg.endswith("]")

    def test_feature_radical_census_line_gains_suffix(self, caplog,
                                                      clean_registry):
        import logging as _logging

        from rmgpy.polymer import _refused_census_warned, _warn_once_refused

        entry = {"reaction": "SYN_R(1) + SYN_P(2) <=> SYN_Q(3)",
                 "radical_class": "eliminating",
                 "reason": "conduit-deferred"}
        _refused_census_warned.discard(
            (entry["reaction"], entry["reason"]))
        with caplog.at_level(_logging.WARNING):
            _warn_once_refused(entry)
        msgs = [r.getMessage() for r in caplog.records
                if r.getMessage().startswith("FEATURE-RADICAL REFUSED "
                                             "CENSUS: ")]
        assert len(msgs) == 1
        msg = msgs[0]
        # historical text intact, suffix appended after it
        assert ("no flux applied (stamp-but-keep). Deferred to item 20's "
                "conduit.") in msg
        assert " [conduit-census/1 key=" in msg
        assert "bucket=FEATURE_RADICAL" in msg
        # and the ledger recorded the upstream blocker
        key = candidate_key_from_label(entry["reaction"])
        assert lookup_candidate(key)["effective_bucket"] == "FEATURE_RADICAL"

    def test_both_census_row_gets_precedence_fields(self, caplog,
                                                    clean_registry):
        """A row censused feature-radical FIRST (upstream) and then hitting
        the r93 branch must carry overlap + precedence in its r93 suffix
        and keep the r93 verdict only as shadow."""
        import logging as _logging

        from rmgpy.polymer import (_general_chain_scale_pool_warned,
                                   _reaction_census_label,
                                   _warn_once_refused,
                                   stamp_gas_association_refusal)

        rxn = self._refused_row()
        label = _reaction_census_label(rxn)
        from rmgpy.polymer import _refused_census_warned
        _refused_census_warned.discard((label, "conduit-deferred"))
        _warn_once_refused({"reaction": label,
                            "radical_class": "eliminating",
                            "reason": "conduit-deferred"})
        _general_chain_scale_pool_warned.clear()
        with caplog.at_level(_logging.WARNING):
            stamp_gas_association_refusal(rxn)
        r93_lines = [r.getMessage() for r in caplog.records
                     if r.getMessage().startswith(self._R93_PREFIX)]
        assert len(r93_lines) == 1
        assert "bucket=FEATURE_RADICAL" in r93_lines[0]
        assert "precedence=feature_radical" in r93_lines[0]
        assert "censuses=feature_radical+r93_general" in r93_lines[0]
        key = candidate_key_from_label(label)
        assert lookup_candidate(key)["effective_bucket"] == "FEATURE_RADICAL"


class TestAdapterDivergenceAndThreshold:
    """Round-36 P1(c) divergence path + the INTEGRATION.md threshold rule
    (gas MW threshold from the row's own pool participants)."""

    def test_gas_threshold_from_row_pools(self, clean_registry):
        from rmgpy.polymer import Polymer
        pp = Polymer(label="PP", monomer="[CH2][CH]C",
                     Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        thr = gas_mw_threshold_for_pools([pp])
        assert thr == pytest.approx(1.5 * pp.monomer_mw_g_mol, rel=1e-12)
        # multi-pool: the tightest (smallest-monomer) pool gates
        big = Polymer(label="BIG", monomer="[CH2][CH]c1ccccc1",
                      Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        assert gas_mw_threshold_for_pools([big, pp]) == pytest.approx(
            1.5 * min(pp.monomer_mw_g_mol, big.monomer_mw_g_mol), rel=1e-12)
        # no usable pool: falls back to the module default
        assert gas_mw_threshold_for_pools([]) == GAS_MW_THRESHOLD

    def test_label_isomorphism_divergence_is_logged_not_overridden(
            self, caplog, clean_registry):
        """A species whose LABEL claims pool-state resolvability
        (trimer_rad33) but whose STRUCTURE does not resolve to the row's
        pool must trigger the CONDUIT CLASSIFIER DIVERGENCE census line,
        and the classification must follow the ISOMORPHISM verdict."""
        import logging as _logging

        from rmgpy.molecule import Molecule
        from rmgpy.polymer import Polymer
        from rmgpy.polymer_conduit import record_from_reaction
        from rmgpy.reaction import Reaction
        from rmgpy.species import Species

        pp = Polymer(label="phenol_formaldehyde", monomer="[CH2][CH]C",
                     Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        impostor = Species(molecule=[Molecule().from_smiles("CCO")])
        impostor.label = "trimer_rad33"      # label says pool-resolvable
        rxn = Reaction(reactants=[impostor], products=[pp],
                       reversible=True)
        with caplog.at_level(_logging.WARNING):
            record = record_from_reaction(rxn, [pp])
        assert record["label_isomorphism_divergence"] is True
        assert any("CONDUIT CLASSIFIER DIVERGENCE" in r.getMessage()
                   for r in caplog.records)
        # isomorphism verdict wins in the role assignment: the impostor is
        # NOT pool-state resolvable, so its role falls back to MW-based
        from rmgpy.polymer_conduit import _apply_iso_overrides, species_role
        rec = _apply_iso_overrides(record)
        assert species_role(rec["reactants"][0]) == "DISC"   # 46 g/mol

    def test_agreeing_species_produces_no_divergence_line(self, caplog,
                                                          clean_registry):
        import logging as _logging

        from rmgpy.molecule import Molecule
        from rmgpy.polymer import Polymer
        from rmgpy.polymer_conduit import record_from_reaction
        from rmgpy.reaction import Reaction
        from rmgpy.species import Species

        pp = Polymer(label="phenol_formaldehyde", monomer="[CH2][CH]C",
                     Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        gas = Species(molecule=[Molecule().from_smiles("C=C")])
        gas.label = "ethylene"
        rxn = Reaction(reactants=[gas], products=[pp], reversible=True)
        with caplog.at_level(_logging.WARNING):
            record = record_from_reaction(rxn, [pp])
        assert record["label_isomorphism_divergence"] is False
        assert not any("CONDUIT CLASSIFIER DIVERGENCE" in r.getMessage()
                       for r in caplog.records)
