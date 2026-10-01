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

        # I-016: a GENUINE feature-radical is an accumulating/QSSA-invalid
        # row (see is_qssa_eliminating_radical); only it buckets FEATURE_RADICAL.
        entry = {"reaction": "SYN_R(1) + SYN_P(2) <=> SYN_Q(3)",
                 "radical_class": "accumulating",
                 "reason": "qssa-invalid"}
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
        # I-016: the feature_radical+r93 overlap precedence only arises for a
        # GENUINE feature-radical (accumulating/QSSA-invalid) upstream sighting.
        _refused_census_warned.discard((label, "qssa-invalid"))
        _warn_once_refused({"reaction": label,
                            "radical_class": "accumulating",
                            "reason": "qssa-invalid"})
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


# ---------------------------------------------------------------------------
# M18.3: admission gates G0-G7 (T3), ledger epochs + run-boundary reset (T8)
# ---------------------------------------------------------------------------

def _admissible_fixture(reversible=True, aligned=True, extra_gas=0,
                        pool_mw=True, kinetics=True, label="CHAIN",
                        gas_veto=True):
    """Live shape-A admission fixture: CHAIN => phenol_formaldehyde + CH2O,
    element-balanced against the pool's real representative molecule (the
    chain is the proxy molecule with one H replaced by CH2OH, so
    chain = proxy + CH2O exactly)."""
    from rmgpy.kinetics import Arrhenius
    from rmgpy.molecule import Molecule
    from rmgpy.polymer import Polymer, set_polymer_gas_veto
    from rmgpy.reaction import Reaction
    from rmgpy.species import Species

    pf = Polymer(label="phenol_formaldehyde",
                 monomer="[CH2]c1ccc(cc1)C([CH2])O",
                 Mn=5000.0, Mw=8000.0, initial_mass=1.0)
    if not pool_mw:
        pf.monomer_mw_g_mol = 0.0
    proxy_smiles = pf.molecule[0].to_smiles()
    assert proxy_smiles.startswith("C")
    chain = Species(molecule=[Molecule().from_smiles("OCC" + proxy_smiles[1:])])
    chain.label = label
    if gas_veto:
        set_polymer_gas_veto(chain)
    gas = Species(molecule=[Molecule().from_smiles("C=O")])
    gas.label = "CH2O"
    products = [pf, gas] + [gas] * extra_gas
    kin = (Arrhenius(A=(1.0, "s^-1"), n=0.0, Ea=(0.0, "J/mol"))
           if kinetics else None)
    if aligned:
        rxn = Reaction(reactants=[chain], products=products,
                       reversible=reversible, kinetics=kin)
    else:
        rxn = Reaction(reactants=products, products=[chain],
                       reversible=reversible, kinetics=kin)
    return rxn, pf


class TestAdmissionGates:
    """T3 -- evaluate_conduit_admission unit tests (fail-closed dead code:
    CONDUIT_ADMISSION_ENABLED stays False; these tests drive the evaluator
    directly)."""

    def test_flag_is_false(self):
        from rmgpy.polymer_conduit import CONDUIT_ADMISSION_ENABLED
        assert CONDUIT_ADMISSION_ENABLED is False

    def test_irreversible_aligned_admits_no_rewrite(self, clean_registry):
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        rxn, pf = _admissible_fixture(reversible=False)
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is True
        assert v.deny_reason is None
        assert v.needs_irreversible_rewrite is False
        assert v.dst_pool == "phenol_formaldehyde"
        # landing-cone numbers: u = (MW(chain) - MW(gas))/M, a = MW(gas)/M
        chain_mw = rxn.reactants[0].molecule[0].get_molecular_weight() * 1e3
        gas_mw = 30.026
        m = pf.monomer_mw_g_mol
        assert v.gas_units == pytest.approx(gas_mw / m, rel=1e-3)
        assert v.chain_units == pytest.approx((chain_mw - gas_mw) / m,
                                              rel=1e-3)
        assert v.chain_units >= 1.0
        assert v.gas_product[0] == "CH2O"

    def test_reversible_aligned_admits_with_rewrite(self, clean_registry):
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        rxn, pf = _admissible_fixture(reversible=True)
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is True
        assert v.needs_irreversible_rewrite is True

    def test_apply_rewrite_flips_and_is_idempotent(self, caplog,
                                                   clean_registry):
        import logging as _logging
        from rmgpy.polymer import _apply_conduit_irreversible_rewrite
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        rxn, pf = _admissible_fixture(reversible=True)
        v = evaluate_conduit_admission(rxn, [pf])
        with caplog.at_level(_logging.WARNING):
            _apply_conduit_irreversible_rewrite(rxn, v)
        assert rxn.reversible is False
        killed = [r.getMessage() for r in caplog.records
                  if "reverse killed" in r.getMessage()]
        assert len(killed) == 1
        # idempotent: second call neither flips back nor re-logs
        with caplog.at_level(_logging.WARNING):
            _apply_conduit_irreversible_rewrite(rxn, v)
        assert rxn.reversible is False
        killed = [r.getMessage() for r in caplog.records
                  if "reverse killed" in r.getMessage()]
        assert len(killed) == 1

    def test_rewrite_noop_without_flag(self, clean_registry):
        from rmgpy.polymer import _apply_conduit_irreversible_rewrite
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        rxn, pf = _admissible_fixture(reversible=False)
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.needs_irreversible_rewrite is False
        rxn.reversible = True  # simulate: rewrite must NOT touch it
        _apply_conduit_irreversible_rewrite(rxn, v)
        assert rxn.reversible is True

    def test_reversible_anti_aligned_denies_flip_rewrite(self,
                                                         clean_registry):
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        rxn, pf = _admissible_fixture(reversible=True, aligned=False)
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        assert v.deny_reason == "direction-requires-flip-rewrite"

    def test_irreversible_anti_aligned_denies_direction(self,
                                                        clean_registry):
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        rxn, pf = _admissible_fixture(reversible=False, aligned=False)
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        assert v.deny_reason == "direction-inadmissible"

    def test_multi_gas_denies_count(self, clean_registry):
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        rxn, pf = _admissible_fixture(extra_gas=1)  # stoich-2 gas product
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        assert v.deny_reason == "gas-product-count"

    def test_missing_pool_mw_denies_unresolvable_no_fallback(
            self, clean_registry):
        """No usable pool monomer MW must DENY -- the census-time module
        default (gas_mw_threshold_for_pools fallback) is forbidden on the
        admission path (uncertainty never admits)."""
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        rxn, pf = _admissible_fixture(pool_mw=False)
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        assert v.deny_reason == "gas-mw-threshold-unresolvable"

    def test_gas_over_threshold_denies(self, clean_registry):
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        rxn, pf = _admissible_fixture()
        pf.monomer_mw_g_mol = 15.0   # threshold 22.5 < CH2O's 30.03
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        assert v.deny_reason == "gas-mw-over-threshold"

    def test_fr_overlapped_key_denies(self, clean_registry):
        from rmgpy.polymer import _reaction_census_label
        from rmgpy.polymer_conduit import (candidate_key_from_label,
                                           evaluate_conduit_admission,
                                           register_candidate)
        rxn, pf = _admissible_fixture()
        key = candidate_key_from_label(_reaction_census_label(rxn))
        register_candidate(key, "feature_radical", "FEATURE_RADICAL")
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        assert v.deny_reason == "feature-radical-overlap"

    def test_classifier_divergence_denies(self, caplog, clean_registry):
        """A label/isomorphism divergence is a finding, never an admission
        basis: relabel the gas product as a pool-state-resolvable chip."""
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        rxn, pf = _admissible_fixture()
        gas = rxn.products[1]
        gas.label = "trimer_rad33"   # label claims pool-state; structure no
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        assert v.deny_reason == "classifier-divergence"

    def test_unbalanced_denies(self, clean_registry):
        from rmgpy.molecule import Molecule
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        from rmgpy.species import Species
        rxn, pf = _admissible_fixture()
        heavy_gas = Species(molecule=[Molecule().from_smiles("CC=O")])
        heavy_gas.label = "CH3CHO"
        rxn.products[1] = heavy_gas   # breaks element balance only
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        assert v.deny_reason == "not-balanced"

    def test_non_arrhenius_kinetics_denies(self, clean_registry):
        from rmgpy.kinetics import Chebyshev
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        rxn, pf = _admissible_fixture()
        rxn.kinetics = Chebyshev()
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        assert v.deny_reason == "kinetics-not-exportable"

    def test_internal_exception_denies_evaluation_error(self, caplog,
                                                        clean_registry):
        import logging as _logging
        from rmgpy.polymer_conduit import evaluate_conduit_admission

        class _Boom:
            def __getattr__(self, name):
                raise RuntimeError("boom")

        with caplog.at_level(_logging.WARNING):
            v = evaluate_conduit_admission(_Boom(), [])
        assert v.admitted is False
        assert v.deny_reason.startswith("admission-evaluation-error:")
        assert any("DENYING fail-closed" in r.getMessage()
                   for r in caplog.records)

    def test_non_chain_scale_chain_denies_classifier(self, clean_registry):
        """A chain that is NOT chain-scale relative to its OWN pool is
        demoted CHAIN->DISC by the pool-DERIVED yardstick (the fixed
        classifier derives the chain-scale bar from the row's live monomer,
        not the hardcoded PF constant). With no CHAIN role on either side the
        G2 direction/orientation gate denies at classification, BEFORE the G5
        phase gate is ever reached -- the row is denied fail-closed either
        way; only the deny label reflects the consistent pool-derived bar."""
        from rmgpy.polymer_conduit import (CHAIN_SCALE_FACTOR,
                                           evaluate_conduit_admission)
        rxn, pf = _admissible_fixture()
        # Chain proxy is 434.57 g/mol / 32 heavy atoms. A large pool monomer
        # MW lifts the pool-relative chain-scale mass bar
        # (CHAIN_SCALE_FACTOR*400 = 1000 g/mol) above the chain, so it is no
        # longer chain-scale against its own pool.
        pf.monomer_mw_g_mol = 400.0
        assert 434.57 < CHAIN_SCALE_FACTOR * pf.monomer_mw_g_mol
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        # Demoted to DISC -> no CHAIN role -> G2 direction gate denies at
        # classification, before G5's chain-not-condensed can ever fire.
        assert v.deny_reason == "classifier-not-admissible"

    def test_chain_scale_non_condensed_denies_g5(self, clean_registry):
        """G5 phase-gate coverage: a chain that IS chain-scale relative to
        its pool (so it passes the G2 direction gate as CHAIN) but carries NO
        melt/proxy-derivation evidence (no is_polymer_proxy tag, no gas veto)
        is NOT melt-classified, so the G5 phase gate denies it
        ``chain-not-condensed`` -- the consumer's step-2 phase gate
        (V_rxn = V_poly) must never be fed a non-condensed event."""
        from rmgpy.polymer import _discrete_is_chain_scale_proxy_derived
        from rmgpy.polymer_conduit import (CHAIN,
                                           chain_scale_thresholds_for_pools,
                                           evaluate_conduit_admission,
                                           record_from_reaction, species_role)
        # Default PF monomer MW keeps the 434.57 g/mol / 32-heavy chain
        # chain-scale (>= 2.5*134.178 g/mol AND >= 2.5*10 heavy), so it passes
        # G2 as CHAIN; gas_veto=False strips the only melt-evidence the
        # fixture stamps, so the G5 phase gate refuses it.
        rxn, pf = _admissible_fixture(reversible=False, gas_veto=False)
        cs_mw, cs_heavy = chain_scale_thresholds_for_pools([pf])
        record = record_from_reaction(rxn, [pf], census="r93_general")
        # precondition: still classifies CHAIN (reaches G5)...
        assert CHAIN in [species_role(s, cs_mw, cs_heavy)
                         for s in record["reactants"]]
        # ...but is NOT proxy-derived / melt-classified (fails the G5 gate).
        assert not _discrete_is_chain_scale_proxy_derived(rxn.reactants[0],
                                                          [pf])
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        assert v.deny_reason == "chain-not-condensed"


class TestAdmissionEndToEndDeadFlag:
    """T3 END-TO-END: with CONDUIT_ADMISSION_ENABLED False, the r93 stamp
    site refuses a WOULD-ADMIT row exactly as before (attrs, census
    reason), and the census suffix carries the verdict."""

    def test_would_admit_row_still_refuses(self, caplog, clean_registry):
        import logging as _logging
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   _general_chain_scale_pool_warned,
                                   stamp_gas_association_refusal)
        rxn, pf = _admissible_fixture(reversible=True)
        _general_chain_scale_pool_warned.clear()
        with caplog.at_level(_logging.WARNING):
            stamp_gas_association_refusal(rxn)
        # byte-for-byte today's refusal state
        assert rxn.polymer_refused is True
        assert rxn.polymer_refused_accumulating is False
        assert rxn.reversible is True                      # NO rewrite
        assert int(getattr(rxn, "polymer_flux_archetype", 0)) != int(
            PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT)
        assert getattr(rxn, "polymer_conduit_params", None) is None  # NO stamp
        msgs = [r.getMessage() for r in caplog.records
                if "[conduit-admission/1" in r.getMessage()]
        assert len(msgs) == 1
        assert "would_admit=1 deny=None stage=final rewrite=True" in msgs[0]
        # append-only: the M18.2 census suffix is intact before it
        assert "[conduit-census/1 key=" in msgs[0]
        assert msgs[0].index("[conduit-census/1") < msgs[0].index(
            "[conduit-admission/1")

    def test_denied_row_suffix_carries_reason(self, caplog, clean_registry):
        import logging as _logging
        from rmgpy.polymer import (_general_chain_scale_pool_warned,
                                   stamp_gas_association_refusal)
        rxn, pf = _admissible_fixture(reversible=True, extra_gas=1)
        _general_chain_scale_pool_warned.clear()
        with caplog.at_level(_logging.WARNING):
            stamp_gas_association_refusal(rxn)
        assert rxn.polymer_refused is True
        msgs = [r.getMessage() for r in caplog.records
                if "[conduit-admission/1" in r.getMessage()]
        assert len(msgs) == 1
        # (round-49) a plain one-shot deny at stamp time IS the final word
        # for its key and says so: stage=final rides adjacent to deny=.
        assert "would_admit=0 deny=gas-product-count stage=final" in msgs[0]


class TestG6ReadjudicationProvisionalKinetics:
    """G6 re-adjudication defect fix (adjudicated): the r93 stamp site
    (stamp_gas_association_refusal <- make_new_reaction) runs BEFORE
    make_new_reaction assigns kinetics, so family-generated rows always
    reached G6 with kinetics=None and were denied kinetics-not-exportable
    -- a reason that said nothing about their eventual kinetics (conduit4
    forensics: exactly 18 ADMISSIBLE_B rows denied solely on this
    artifact). Now: kinetics=None denies PROVISIONALLY
    (kinetics-not-yet-assigned, NEVER would_admit=1) and
    readjudicate_conduit_admission resolves the FINAL verdict once the
    kinetics exists -- on the new-reaction path after the kinetics
    conversion/barrier block, and across the check_existing canonical-dedup
    merge (the trap: a provisional stamp merged onto an existing canonical
    must not escape re-adjudication)."""

    def test_vocabulary_gains_exactly_the_provisional_token(self):
        from rmgpy.polymer_conduit import (ADMISSION_DENY_REASONS,
                                           PROVISIONAL_DENY_REASONS)
        assert "kinetics-not-yet-assigned" in ADMISSION_DENY_REASONS
        assert PROVISIONAL_DENY_REASONS == {"kinetics-not-yet-assigned"}
        # provisional reasons stay inside the closed deny vocabulary
        assert PROVISIONAL_DENY_REASONS <= ADMISSION_DENY_REASONS

    def test_stage_token_partitions_the_closed_deny_vocabulary(self):
        """(round-49) every census line carries exactly one stage= token,
        adjacent to the would_admit=/deny= tokens: stage=provisional iff
        the verdict is a PROVISIONAL_DENY_REASONS member (subject to G6
        re-adjudication), stage=final for every other deny AND for the
        admit line -- so grep tallies filtered on stage=final never
        double-count a re-adjudicated candidate."""
        from rmgpy.polymer_conduit import (ADMISSION_DENY_REASONS,
                                           PROVISIONAL_DENY_REASONS,
                                           AdmissionVerdict,
                                           admission_census_suffix, _deny)
        for reason in sorted(ADMISSION_DENY_REASONS):
            suffix = admission_census_suffix(_deny("K", reason))
            expected = ("provisional" if reason in PROVISIONAL_DENY_REASONS
                        else "final")
            assert suffix.endswith(f"deny={reason} stage={expected}]")
            assert suffix.count("stage=") == 1
        admit = admission_census_suffix(AdmissionVerdict(
            admitted=True, deny_reason=None, candidate_key="K"))
        assert ("would_admit=1 deny=None stage=final rewrite=False"
                in admit)
        assert admit.count("stage=") == 1

    def test_kinetics_none_denies_provisionally(self, clean_registry):
        """(a) kinetics=None -> provisional kinetics-not-yet-assigned;
        a provisional row NEVER has would_admit=1."""
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        rxn, pf = _admissible_fixture(kinetics=False)
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        assert v.deny_reason == "kinetics-not-yet-assigned"

    def test_stamp_site_marks_pending_and_censuses_provisional(
            self, caplog, clean_registry):
        import logging as _logging
        from rmgpy.polymer import (_general_chain_scale_pool_warned,
                                   stamp_gas_association_refusal)
        rxn, pf = _admissible_fixture(kinetics=False)
        _general_chain_scale_pool_warned.clear()
        with caplog.at_level(_logging.WARNING):
            stamp_gas_association_refusal(rxn)
        # refusal state byte-identical to any other r93 refusal
        assert rxn.polymer_refused is True
        assert rxn.polymer_refused_accumulating is False
        assert rxn.polymer_conduit_admission_pending is True
        msgs = [r.getMessage() for r in caplog.records
                if "[conduit-admission/1" in r.getMessage()]
        assert len(msgs) == 1
        # (round-49) the provisional line is machine-distinguishable from
        # the later FINAL re-adjudication line: stage=provisional.
        assert ("would_admit=0 deny=kinetics-not-yet-assigned "
                "stage=provisional") in msgs[0]
        assert "would_admit=1" not in msgs[0]
        assert "stage=final" not in msgs[0]

    def test_readjudication_with_final_arrhenius_would_admit(
            self, caplog, clean_registry):
        """(b) provisional row + final plain Arrhenius -> would_admit=1
        (census-only: admission stays disabled, refusal untouched)."""
        import logging as _logging
        from rmgpy.kinetics import Arrhenius
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   _general_chain_scale_pool_warned,
                                   readjudicate_conduit_admission,
                                   stamp_gas_association_refusal)
        rxn, pf = _admissible_fixture(kinetics=False)
        _general_chain_scale_pool_warned.clear()
        stamp_gas_association_refusal(rxn)
        rxn.kinetics = Arrhenius(A=(1.0, "s^-1"), n=0.0, Ea=(0.0, "J/mol"))
        with caplog.at_level(_logging.WARNING):
            readjudicate_conduit_admission(rxn)
        assert rxn.polymer_conduit_admission_pending is False
        assert rxn.polymer_conduit_admission_readjudicated is True
        msgs = [r.getMessage() for r in caplog.records
                if "G6 re-adjudication" in r.getMessage()
                and "[conduit-admission/1" in r.getMessage()]
        assert len(msgs) == 1
        assert "would_admit=1 deny=None stage=final rewrite=True" in msgs[0]
        # census-only while CONDUIT_ADMISSION_ENABLED is False: the refusal
        # state is untouched, nothing is stamped, nothing is rewritten.
        assert rxn.polymer_refused is True
        assert rxn.reversible is True
        assert int(getattr(rxn, "polymer_flux_archetype", 0)) != int(
            PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT)
        assert getattr(rxn, "polymer_conduit_params", None) is None

    @pytest.mark.parametrize("kin_name", ["Chebyshev", "MultiArrhenius"])
    def test_readjudication_non_arrhenius_final_denies_exportable(
            self, caplog, clean_registry, kin_name):
        """(c) provisional row + genuinely non-Arrhenius final kinetics ->
        FINAL kinetics-not-exportable (the strict type check stays the
        fail-closed last word)."""
        import logging as _logging
        import rmgpy.kinetics as _kinetics
        from rmgpy.polymer import (_general_chain_scale_pool_warned,
                                   readjudicate_conduit_admission,
                                   stamp_gas_association_refusal)
        rxn, pf = _admissible_fixture(kinetics=False)
        _general_chain_scale_pool_warned.clear()
        stamp_gas_association_refusal(rxn)
        rxn.kinetics = getattr(_kinetics, kin_name)()
        with caplog.at_level(_logging.WARNING):
            readjudicate_conduit_admission(rxn)
        assert rxn.polymer_conduit_admission_pending is False
        assert rxn.polymer_conduit_admission_readjudicated is True
        msgs = [r.getMessage() for r in caplog.records
                if "G6 re-adjudication" in r.getMessage()
                and "[conduit-admission/1" in r.getMessage()]
        assert len(msgs) == 1
        # (round-49) a re-adjudicated FINAL deny is the last word: stage=final.
        assert ("would_admit=0 deny=kinetics-not-exportable "
                "stage=final") in msgs[0]

    def test_readjudication_waits_while_kinetics_still_none(
            self, caplog, clean_registry):
        """generate_kinetics-disabled callers: the hook keeps the pending
        marker and the provisional deny stands (uncertainty never
        admits)."""
        import logging as _logging
        from rmgpy.polymer import (_general_chain_scale_pool_warned,
                                   readjudicate_conduit_admission,
                                   stamp_gas_association_refusal)
        rxn, pf = _admissible_fixture(kinetics=False)
        _general_chain_scale_pool_warned.clear()
        stamp_gas_association_refusal(rxn)
        with caplog.at_level(_logging.WARNING):
            readjudicate_conduit_admission(rxn)
        assert rxn.polymer_conduit_admission_pending is True
        assert getattr(rxn, "polymer_conduit_admission_readjudicated",
                       False) is False
        assert not any("G6 re-adjudication" in r.getMessage()
                       for r in caplog.records)

    def test_dedup_merge_transfers_pending_and_resolves_on_canonical(
            self, caplog, clean_registry):
        """(d) THE TRAP: check_for_existing_reaction's early return
        discards the freshly (provisionally) stamped candidate. Mirrors
        make_new_reaction's dedup arm exactly -- merge, then re-adjudicate
        the CANONICAL object against its own final kinetics."""
        import logging as _logging
        from rmgpy.polymer import (_general_chain_scale_pool_warned,
                                   merge_polymer_adjudication_stamps,
                                   readjudicate_conduit_admission,
                                   stamp_gas_association_refusal)
        canonical, pf = _admissible_fixture()  # final plain Arrhenius
        dup, _ = _admissible_fixture(kinetics=False)  # pre-kinetics twin
        _general_chain_scale_pool_warned.clear()
        stamp_gas_association_refusal(dup)
        assert dup.polymer_conduit_admission_pending is True
        # make_new_reaction dedup arm: merge stamps, then re-adjudicate rxn
        merge_polymer_adjudication_stamps(dup, canonical)
        assert canonical.polymer_conduit_admission_pending is True
        assert canonical.polymer_refused is True  # refusal merged too
        with caplog.at_level(_logging.WARNING):
            readjudicate_conduit_admission(canonical)
        assert canonical.polymer_conduit_admission_pending is False
        assert canonical.polymer_conduit_admission_readjudicated is True
        msgs = [r.getMessage() for r in caplog.records
                if "G6 re-adjudication" in r.getMessage()
                and "[conduit-admission/1" in r.getMessage()]
        assert len(msgs) == 1
        assert "would_admit=1 deny=None stage=final rewrite=True" in msgs[0]

    def test_merge_never_demotes_a_final_verdict(self, caplog,
                                                 clean_registry):
        """(d) corollary: a canonical whose verdict is already FINAL
        (re-adjudicated) is never demoted back to provisional by a later
        duplicate's merge -- and re-running the hook emits nothing new."""
        import logging as _logging
        from rmgpy.kinetics import Arrhenius
        from rmgpy.polymer import (_general_chain_scale_pool_warned,
                                   merge_polymer_adjudication_stamps,
                                   readjudicate_conduit_admission,
                                   stamp_gas_association_refusal)
        canonical, pf = _admissible_fixture(kinetics=False)
        _general_chain_scale_pool_warned.clear()
        stamp_gas_association_refusal(canonical)
        canonical.kinetics = Arrhenius(A=(1.0, "s^-1"), n=0.0,
                                       Ea=(0.0, "J/mol"))
        readjudicate_conduit_admission(canonical)
        assert canonical.polymer_conduit_admission_readjudicated is True
        dup, _ = _admissible_fixture(kinetics=False)
        stamp_gas_association_refusal(dup)
        assert dup.polymer_conduit_admission_pending is True
        merge_polymer_adjudication_stamps(dup, canonical)
        assert canonical.polymer_conduit_admission_pending is False
        caplog.clear()  # drop the first (legitimate) re-adjudication line
        with caplog.at_level(_logging.WARNING):
            readjudicate_conduit_admission(canonical)
        assert not any("G6 re-adjudication" in r.getMessage()
                       for r in caplog.records)

    def test_census_lines_deterministic_final_verdict_last(
            self, caplog, clean_registry):
        """(e) census-line determinism: the stamp-time line carries the
        provisional token; the re-adjudication line carries the FINAL
        verdict, appears LAST, and keeps the append-only census+admission
        suffix shape."""
        import logging as _logging
        from rmgpy.kinetics import Chebyshev
        from rmgpy.polymer import (_general_chain_scale_pool_warned,
                                   readjudicate_conduit_admission,
                                   stamp_gas_association_refusal)
        rxn, pf = _admissible_fixture(kinetics=False)
        _general_chain_scale_pool_warned.clear()
        with caplog.at_level(_logging.WARNING):
            stamp_gas_association_refusal(rxn)
            rxn.kinetics = Chebyshev()
            readjudicate_conduit_admission(rxn)
        msgs = [r.getMessage() for r in caplog.records
                if "[conduit-admission/1" in r.getMessage()]
        assert len(msgs) == 2
        assert ("would_admit=0 deny=kinetics-not-yet-assigned "
                "stage=provisional") in msgs[0]
        assert ("would_admit=0 deny=kinetics-not-exportable "
                "stage=final") in msgs[1]
        # (round-49) the stage= token is the double-count discriminator:
        # a grep tally filtered on stage=final sees this candidate exactly
        # once even though it emitted two census lines.
        assert sum("stage=final" in m for m in msgs) == 1
        assert sum("stage=provisional" in m for m in msgs) == 1
        # the FINAL line keeps the append-only suffix shape (census tokens
        # first, admission tokens after -- round-36 P1(b) ordering)
        assert "[conduit-census/1 key=" in msgs[1]
        assert msgs[1].index("[conduit-census/1") < msgs[1].index(
            "[conduit-admission/1")


class TestR42AdmissionCensusDeterminism:
    """r42 P1-4 -- (a) EVERY census line carries the admission verdict
    tokens (the FR line lacked them); (b) the would_admit evaluation must
    be deterministic w.r.t. ledger state that includes the current row
    (an FR sighting landing between evaluation and the census line / the
    stamping decision must win)."""

    def test_feature_radical_line_carries_admission_suffix(
            self, caplog, clean_registry):
        """r42 P1-4(a): the FEATURE-RADICAL census line carries the
        [conduit-admission/1 ...] tokens like every other census line."""
        import logging as _logging

        from rmgpy.polymer import _refused_census_warned, _warn_once_refused

        # I-016: only a genuine (accumulating/QSSA-invalid) feature-radical
        # carries the feature-radical-overlap admission token.
        entry = {"reaction": "SYN_R42(1) + SYN_P42(2) <=> SYN_Q42(3)",
                 "radical_class": "accumulating",
                 "reason": "qssa-invalid"}
        _refused_census_warned.discard((entry["reaction"], entry["reason"]))
        with caplog.at_level(_logging.WARNING):
            _warn_once_refused(entry)
        msgs = [r.getMessage() for r in caplog.records
                if r.getMessage().startswith("FEATURE-RADICAL REFUSED "
                                             "CENSUS: ")]
        assert len(msgs) == 1
        assert ("[conduit-admission/1 would_admit=0 "
                "deny=feature-radical-overlap stage=final]") in msgs[0]
        # append-only: the M18.2 census suffix rides before it
        assert msgs[0].index("[conduit-census/1") < msgs[0].index(
            "[conduit-admission/1")

    def test_i016_only_accumulating_gets_sticky_feature_radical_tag(
            self, clean_registry):
        """I-016 tag-narrowing lock: an accumulating/QSSA-invalid row
        registers the sticky ``feature_radical`` census (the admission
        gate-G1 blocker); an eliminating/conduit-deferred row registers the
        non-blocking ``conduit_deferred`` token instead. Because G1 keys on
        ``"feature_radical" in censuses``, the eliminating row is no longer
        permanently blocked -- the self-referential trap is gone. Admission
        is default-off, so this is ZERO generation change (the un-blocked
        row still re-denies downstream at G2/G3)."""
        from rmgpy.polymer import _refused_census_warned, _warn_once_refused
        from rmgpy.polymer_conduit import (lookup_candidate,
                                           candidate_key_from_label)

        acc_label = "ACC_R(1) <=> ACC_Q(2)"
        elim_label = "ELM_R(1) + ELM_P(2) <=> ELM_Q(3)"
        _refused_census_warned.discard((acc_label, "qssa-invalid"))
        _refused_census_warned.discard((elim_label, "conduit-deferred"))

        _warn_once_refused({"reaction": acc_label,
                            "radical_class": "accumulating",
                            "reason": "qssa-invalid"})
        _warn_once_refused({"reaction": elim_label,
                            "radical_class": "eliminating",
                            "reason": "conduit-deferred"})

        acc = lookup_candidate(candidate_key_from_label(acc_label))
        elim = lookup_candidate(candidate_key_from_label(elim_label))
        assert acc is not None and elim is not None
        # genuine feature-radical keeps the sticky G1-blocking tag ...
        assert acc["effective_bucket"] == "FEATURE_RADICAL"
        assert "feature_radical" in acc["censuses"]
        # ... conduit-deferred registers non-blocking; NOT feature_radical,
        # so G1 (which tests membership of "feature_radical") cannot fire.
        assert elim["effective_bucket"] == "CONDUIT_DEFERRED"
        assert "feature_radical" not in elim["censuses"]

    def test_fr_registered_after_evaluation_downgrades_would_admit(
            self, clean_registry):
        """r42 P1-4(b), census-line half: a verdict evaluated BEFORE the
        FR sighting landed must NOT print would_admit=1 -- the suffix is
        re-checked against the POST-registration ledger entry."""
        from rmgpy.polymer_conduit import (annotate_refused_row,
                                           evaluate_conduit_admission,
                                           register_candidate)
        rxn, pf = _admissible_fixture(reversible=True)
        verdict = evaluate_conduit_admission(rxn, [pf])
        assert verdict.admitted is True          # ledger has no FR yet
        # the FR census sights the same key AFTER evaluation:
        register_candidate(verdict.candidate_key, "feature_radical",
                           "FEATURE_RADICAL")
        suffix = annotate_refused_row(rxn, [pf], census="r93_general",
                                      verdict=verdict)
        assert "would_admit=1" not in suffix
        assert ("[conduit-admission/1 would_admit=0 "
                "deny=feature-radical-overlap stage=final]") in suffix

    def test_stamp_site_rechecks_ledger_before_admitting(
            self, caplog, clean_registry, monkeypatch):
        """r42 P1-4(b), decision-point half: even under a (test-scoped)
        flipped admission flag, an FR sighting that lands between
        evaluation and stamping deterministically blocks the admit arm --
        the row refuses exactly as today (fail-closed narrowing; the
        production flag stays False, pinned by test_flag_is_false)."""
        import logging as _logging

        import rmgpy.polymer_conduit as pc
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   _general_chain_scale_pool_warned,
                                   stamp_gas_association_refusal)

        real_evaluate = pc.evaluate_conduit_admission

        def evaluate_then_fr_sighting(forward, row_pools):
            v = real_evaluate(forward, row_pools)
            # the FR census fires for the same key right after evaluation
            pc.register_candidate(v.candidate_key, "feature_radical",
                                  "FEATURE_RADICAL")
            return v

        monkeypatch.setattr(pc, "evaluate_conduit_admission",
                            evaluate_then_fr_sighting)
        monkeypatch.setattr(pc, "CONDUIT_ADMISSION_ENABLED", True)
        rxn, pf = _admissible_fixture(reversible=True)
        _general_chain_scale_pool_warned.clear()
        with caplog.at_level(_logging.WARNING):
            stamp_gas_association_refusal(rxn)
        # the admit arm did NOT run: refusal state byte-identical to today
        assert rxn.polymer_refused is True
        assert rxn.reversible is True                       # NO rewrite
        assert int(getattr(rxn, "polymer_flux_archetype", 0)) != int(
            PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT)
        assert getattr(rxn, "polymer_conduit_params", None) is None
        msgs = [r.getMessage() for r in caplog.records
                if "[conduit-admission/1" in r.getMessage()]
        assert len(msgs) == 1
        assert "would_admit=0 deny=feature-radical-overlap" in msgs[0]


class TestM18LiveAdmissionRelocation:
    """OQ-2 relocation (M18.4 wiring): the admit-stamp + [r39-P1]
    irreversible export rewrite has moved OUT of the pre-dedup
    ``stamp_gas_association_refusal`` arm and INTO
    ``readjudicate_conduit_admission``, so it always lands on the
    CANONICAL, post-dedup, post-kinetics object -- never on a candidate
    dedup can silently discard.

    Two defects this closes:

    * CLOBBER: :func:`merge_polymer_adjudication_stamps` never carries
      ``.reversible``, so a live rewrite applied to a PRE-DEDUP candidate
      that dedup then discards (keeping a pre-existing REVERSIBLE
      canonical instead) left the surviving canonical stamped
      MOMENT_CREDIT_CONDUIT while still exporting a spurious ``<=>``.
    * FAMILY UNDER-WIRING: family-generated TemplateReactions reach the
      stamp site with ``kinetics is None`` (provisional deny) and never
      got the rewrite even once their final kinetics resolved.

    All three tests here monkeypatch ``CONDUIT_ADMISSION_ENABLED`` True
    for the duration only (the production flag itself stays False,
    pinned by ``TestAdmissionGates.test_flag_is_false``)."""

    def test_dedup_clobber_relocated_admit_rewrites_canonical(
            self, caplog, clean_registry, monkeypatch):
        """THE clobber-closing test: a pre-existing REVERSIBLE canonical
        duplicate is kept by dedup over a freshly (provisionally) stamped
        admissible twin. Before the relocation, the stamp site would have
        mutated ONLY the discarded twin (or, prior to that, never have
        been reached at all for the canonical); merge_polymer_adjudication_stamps
        does not carry .reversible, so the canonical would survive
        MOMENT_CREDIT_CONDUIT-stamped yet still reversible=True. After the
        relocation, readjudicate_conduit_admission on the canonical performs
        the live admit + rewrite, so it must end up reversible=False."""
        import logging as _logging

        import rmgpy.polymer_conduit as pc
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   _general_chain_scale_pool_warned,
                                   merge_polymer_adjudication_stamps,
                                   readjudicate_conduit_admission,
                                   stamp_gas_association_refusal)

        monkeypatch.setattr(pc, "CONDUIT_ADMISSION_ENABLED", True)
        # pre-existing REVERSIBLE canonical, final Arrhenius kinetics
        canonical, pf = _admissible_fixture(reversible=True)
        # a fresh admissible twin, still pre-kinetics (family-shaped)
        dup, _ = _admissible_fixture(reversible=True, kinetics=False)
        _general_chain_scale_pool_warned.clear()
        stamp_gas_association_refusal(dup)
        assert dup.polymer_conduit_admission_pending is True
        assert dup.reversible is True  # NO mutation ever lands on dup

        # make_new_reaction's dedup early-return: dup is discarded, its
        # stamps merge onto the surviving canonical, which is then
        # re-adjudicated against its OWN final kinetics.
        merge_polymer_adjudication_stamps(dup, canonical)
        assert canonical.polymer_conduit_admission_pending is True
        assert canonical.reversible is True  # not yet rewritten

        with caplog.at_level(_logging.WARNING):
            readjudicate_conduit_admission(canonical)

        assert canonical.reversible is False
        assert int(canonical.polymer_flux_archetype) == int(
            PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT)
        assert canonical.polymer_refused is False
        assert canonical.polymer_conduit_params is not None
        assert any("MOMENT-CREDIT CONDUIT ADMITTED" in r.getMessage()
                   for r in caplog.records)
        # the discarded dup itself was NEVER mutated -- the clobber this
        # relocation closes is exactly "mutate the object dedup discards".
        assert dup.reversible is True
        assert int(getattr(dup, "polymer_flux_archetype", 0)) != int(
            PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT)

    def test_whole_model_backstop_every_conduit_stamp_is_irreversible(
            self, caplog, clean_registry, monkeypatch):
        """Whole-model invariant: ANY reaction carrying the
        MOMENT_CREDIT_CONDUIT archetype must be irreversible -- across a
        small mixed population of an admitted row and a denied row, only
        the admitted one is stamped, and it alone is irreversible."""
        import rmgpy.polymer_conduit as pc
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   _general_chain_scale_pool_warned,
                                   readjudicate_conduit_admission,
                                   stamp_gas_association_refusal)

        monkeypatch.setattr(pc, "CONDUIT_ADMISSION_ENABLED", True)
        admitted_rxn, _ = _admissible_fixture(reversible=True, label="CHAIN1")
        denied_rxn, _ = _admissible_fixture(reversible=True, extra_gas=1,
                                            label="CHAIN2")  # gas-product-count deny
        reaction_list = [admitted_rxn, denied_rxn]
        for rxn in reaction_list:
            _general_chain_scale_pool_warned.clear()
            stamp_gas_association_refusal(rxn)
            readjudicate_conduit_admission(rxn)

        conduit_archetype = int(PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT)
        conduit_rows = [r for r in reaction_list
                        if int(getattr(r, "polymer_flux_archetype", 0))
                        == conduit_archetype]
        assert len(conduit_rows) == 1
        assert conduit_rows[0] is admitted_rxn
        for r in conduit_rows:
            assert r.reversible is False

    def test_family_row_kinetics_none_then_final_gets_admitted_and_rewritten(
            self, caplog, clean_registry, monkeypatch):
        """FAMILY under-wiring-closing test: a family-generated
        TemplateReaction reaches the stamp site with kinetics=None (denies
        provisionally, kinetics-not-yet-assigned) and only becomes
        admissible once its final kinetics is assigned (simulating
        make_new_reaction's kinetics conversion step); readjudicate must
        then perform the live admit + rewrite."""
        import logging as _logging

        import rmgpy.polymer_conduit as pc
        from rmgpy.kinetics import Arrhenius
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   _general_chain_scale_pool_warned,
                                   readjudicate_conduit_admission,
                                   stamp_gas_association_refusal)

        monkeypatch.setattr(pc, "CONDUIT_ADMISSION_ENABLED", True)
        rxn, pf = _admissible_fixture(reversible=True, kinetics=False)
        _general_chain_scale_pool_warned.clear()
        with caplog.at_level(_logging.WARNING):
            stamp_gas_association_refusal(rxn)
        assert rxn.polymer_conduit_admission_pending is True
        assert rxn.reversible is True  # no mutation while kinetics unknown
        assert rxn.polymer_conduit_params is None

        rxn.kinetics = Arrhenius(A=(1.0, "s^-1"), n=0.0, Ea=(0.0, "J/mol"))
        with caplog.at_level(_logging.WARNING):
            readjudicate_conduit_admission(rxn)

        assert rxn.polymer_conduit_admission_pending is False
        assert rxn.polymer_conduit_admission_readjudicated is True
        assert rxn.polymer_refused is False
        assert rxn.reversible is False
        assert int(rxn.polymer_flux_archetype) == int(
            PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT)
        assert rxn.polymer_conduit_params is not None
        assert any("MOMENT-CREDIT CONDUIT ADMITTED" in r.getMessage()
                   for r in caplog.records)


class TestM18OptInOverride:
    """M18.4 opt-in wiring: ``readjudicate_conduit_admission`` now takes a
    per-run ``admission_enabled`` override (threaded from
    ``CoreEdgeReactionModel.conduit_admission_enabled``, itself set from the
    ``options(polymerConduitAdmission=...)`` input flag). The precedence
    contract these tests LOCK:

    * ``admission_enabled=True`` ADMITS even while the module constant
      ``CONDUIT_ADMISSION_ENABLED`` stays False -- the whole point of the
      opt-in: a deck enables live admission WITHOUT the global flag being
      flipped (so no OTHER deck's generation behavior moves).
    * ``admission_enabled=False`` REFUSES even if the module constant were
      True -- an explicit per-deck opt-OUT wins over any global enable.
    * ``admission_enabled=None`` (the default, and what every non-model
      caller passes) INHERITS the module constant -- the default-off
      fallback, and the reason the existing relocation tests still pass.

    The production module constant is pinned False by
    ``TestAdmissionGates.test_flag_is_false``; these tests never flip it to
    prove the True case -- that is exactly what the per-call override buys."""

    def test_optin_true_admits_while_module_constant_false(
            self, caplog, clean_registry):
        """Core opt-in proof: admission_enabled=True admits with the module
        constant still at its production False -- no global flip needed."""
        import logging as _logging

        import rmgpy.polymer_conduit as pc
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   _general_chain_scale_pool_warned,
                                   readjudicate_conduit_admission,
                                   stamp_gas_association_refusal)

        assert pc.CONDUIT_ADMISSION_ENABLED is False  # not touched
        rxn, _ = _admissible_fixture(reversible=True)
        _general_chain_scale_pool_warned.clear()
        stamp_gas_association_refusal(rxn)
        assert rxn.polymer_conduit_admission_pending is True

        with caplog.at_level(_logging.WARNING):
            readjudicate_conduit_admission(rxn, admission_enabled=True)

        assert int(rxn.polymer_flux_archetype) == int(
            PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT)
        assert rxn.polymer_refused is False
        assert rxn.reversible is False
        assert rxn.polymer_conduit_params is not None
        assert any("MOMENT-CREDIT CONDUIT ADMITTED" in r.getMessage()
                   for r in caplog.records)

    def test_optin_false_refuses_even_when_module_constant_true(
            self, caplog, clean_registry, monkeypatch):
        """Explicit opt-OUT wins: admission_enabled=False refuses even with
        the module constant monkeypatched True (a hypothetical global
        enable) -- the row stays census-only and reversible."""
        import logging as _logging

        import rmgpy.polymer_conduit as pc
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   _general_chain_scale_pool_warned,
                                   readjudicate_conduit_admission,
                                   stamp_gas_association_refusal)

        monkeypatch.setattr(pc, "CONDUIT_ADMISSION_ENABLED", True)
        rxn, _ = _admissible_fixture(reversible=True)
        _general_chain_scale_pool_warned.clear()
        stamp_gas_association_refusal(rxn)
        assert rxn.polymer_conduit_admission_pending is True

        with caplog.at_level(_logging.WARNING):
            readjudicate_conduit_admission(rxn, admission_enabled=False)

        # re-adjudicated but NOT admitted: no conduit stamp, still reversible
        assert rxn.polymer_conduit_admission_readjudicated is True
        assert int(getattr(rxn, "polymer_flux_archetype", 0)) != int(
            PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT)
        assert rxn.reversible is True
        assert not any("MOMENT-CREDIT CONDUIT ADMITTED" in r.getMessage()
                       for r in caplog.records)

    def test_optin_none_inherits_module_constant_both_directions(
            self, caplog, clean_registry, monkeypatch):
        """admission_enabled=None inherits the module constant: refuses
        while it is False, admits once it is True -- the default-off
        fallback that keeps the existing (no-param) call sites correct."""
        import logging as _logging

        import rmgpy.polymer_conduit as pc
        from rmgpy.polymer import (PolymerFluxArchetype,
                                   _general_chain_scale_pool_warned,
                                   readjudicate_conduit_admission,
                                   stamp_gas_association_refusal)

        conduit_archetype = int(PolymerFluxArchetype.MOMENT_CREDIT_CONDUIT)

        # constant False -> None inherits -> refuse
        assert pc.CONDUIT_ADMISSION_ENABLED is False
        off_rxn, _ = _admissible_fixture(reversible=True, label="OFFCHAIN")
        _general_chain_scale_pool_warned.clear()
        stamp_gas_association_refusal(off_rxn)
        readjudicate_conduit_admission(off_rxn, admission_enabled=None)
        assert int(getattr(off_rxn, "polymer_flux_archetype", 0)) != \
            conduit_archetype
        assert off_rxn.reversible is True

        # constant True -> None inherits -> admit
        monkeypatch.setattr(pc, "CONDUIT_ADMISSION_ENABLED", True)
        on_rxn, _ = _admissible_fixture(reversible=True, label="ONCHAIN")
        _general_chain_scale_pool_warned.clear()
        stamp_gas_association_refusal(on_rxn)
        with caplog.at_level(_logging.WARNING):
            readjudicate_conduit_admission(on_rxn, admission_enabled=None)
        assert int(on_rxn.polymer_flux_archetype) == conduit_archetype
        assert on_rxn.reversible is False


class TestLandingConeEqualityBoundary:
    """T6 (admission-side leg) -- the ratified equality-boundary fixture:
    u_raw EXACTLY a + 1 (credit u == 1.0, point mass ON the cone). The
    closed `>=` semantics do NOT freeze until this fixture reports (OQ-4).

    Constructed at the exact class boundary: chain = 2.5 monomer-equiv of
    mass, gas = 1.5 monomer-equiv (proxy = 1.0), compositionally
    proportional molecules (C25H50 -> C10H20 + C15H30) so every ratio is
    exact up to float rounding."""

    def _boundary_fixture(self):
        from rmgpy.kinetics import Arrhenius
        from rmgpy.molecule import Molecule
        from rmgpy.polymer import Polymer, set_polymer_gas_veto
        from rmgpy.reaction import Reaction
        from rmgpy.species import Species

        pf = Polymer(label="phenol_formaldehyde", monomer="[CH2][CH2]",
                     Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        proxy = Molecule().from_smiles("C=C" + "C" * 8)      # C10H20
        pf.molecule = [proxy]
        gas = Species(molecule=[Molecule().from_smiles("C=C" + "C" * 13)])
        gas.label = "C15H30"
        chain = Species(molecule=[Molecule().from_smiles("C=C" + "C" * 23)])
        chain.label = "CHAIN25"
        set_polymer_gas_veto(chain)
        gas_mw = gas.molecule[0].get_molecular_weight() * 1000.0
        # M := gas_mw / 1.5 == proxy MW (compositional 2:3 proportionality)
        pf.monomer_mw_g_mol = gas_mw / 1.5
        rxn = Reaction(reactants=[chain], products=[pf, gas],
                       reversible=False,
                       kinetics=Arrhenius(A=(1.0, "s^-1"), n=0.0,
                                          Ea=(0.0, "J/mol")))
        return rxn, pf

    def test_boundary_row_report(self, clean_registry):
        """REPORTING test (rider): drive the exact-boundary row through the
        full gate chain and pin whatever the closed `>=` semantics produce.
        As built, u sits AT 1.0 up to float rounding; the assertion pins
        the ADMIT outcome and u == 1.0 -- if float rounding ever flips this
        to a landing-cone denial, that is the boundary-stiffness signal
        that moves the guard to `1 + eps` (OQ-4), not a regression."""
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        rxn, pf = self._boundary_fixture()
        v = evaluate_conduit_admission(rxn, [pf])
        # chain is EXACTLY 2.5 monomer-equiv (closed >=): still chain-scale
        # gas is EXACTLY 1.5 monomer-equiv (closed <=): still under
        assert v.admitted is True, (
            f"equality-boundary row denied ({v.deny_reason}): if this is "
            f"float rounding at the closed boundary, record it for OQ-4")
        assert v.chain_units == pytest.approx(1.0, abs=1e-12)
        assert v.gas_units == pytest.approx(1.5, abs=1e-12)


class TestLedgerEpochsAndReset:
    """T8 -- ledger epochs, sticky FR across epochs, run-boundary reset
    (both-or-neither with the warn-once census sets)."""

    KEY = "EPOCH_chain(1)<>EPOCH_gas(2)+EPOCH_pool(3)"

    def test_epoch_recorded_per_sighting(self, clean_registry):
        from rmgpy.polymer_conduit import (lookup_candidate,
                                           register_candidate)
        register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A",
                           epoch="sig-1")
        register_candidate(self.KEY, "feature_radical", "FEATURE_RADICAL",
                           epoch="sig-3")
        entry = lookup_candidate(self.KEY)
        assert entry["epochs"] == {"sig-1", "sig-3"}

    def test_epoch_none_allowed(self, clean_registry):
        from rmgpy.polymer_conduit import (lookup_candidate,
                                           register_candidate)
        register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A")
        assert lookup_candidate(self.KEY)["epochs"] == set()

    def test_fr_sticky_across_epochs(self, clean_registry):
        """An FR sighting in epoch 1 still blocks admission when the row is
        re-evaluated in epoch 3 (sticky, fail-closed: it can only defer)."""
        from rmgpy.polymer import _reaction_census_label
        from rmgpy.polymer_conduit import (candidate_key_from_label,
                                           evaluate_conduit_admission,
                                           register_candidate)
        rxn, pf = _admissible_fixture()
        key = candidate_key_from_label(_reaction_census_label(rxn))
        register_candidate(key, "feature_radical", "FEATURE_RADICAL",
                           epoch="epoch-1")
        register_candidate(key, "r93_general", "ADMISSIBLE_A",
                           epoch="epoch-3")
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        assert v.deny_reason == "feature-radical-overlap"

    def test_reset_clears_registry_and_warn_once_together(self):
        """DESIGN §3.3 'reset both or neither': reset_conduit_state clears
        the ledger AND _refused_census_warned, so a post-reset FR
        re-sighting re-registers AND re-warns."""
        import logging as _logging
        from rmgpy.polymer import (_refused_census_warned,
                                   _warn_once_refused)
        from rmgpy.polymer_conduit import (lookup_candidate,
                                           register_candidate,
                                           reset_conduit_state,
                                           candidate_key_from_label)
        label = "RESET_R(1) + RESET_P(2) <=> RESET_Q(3)"
        entry = {"reaction": label, "radical_class": "eliminating",
                 "reason": "conduit-deferred"}
        _warn_once_refused(entry)
        key = candidate_key_from_label(label)
        assert lookup_candidate(key) is not None
        assert (label, "conduit-deferred") in _refused_census_warned
        reset_conduit_state()
        assert lookup_candidate(key) is None
        assert (label, "conduit-deferred") not in _refused_census_warned
        # post-reset re-sighting re-registers (would be starved if the
        # warn-once set had survived the ledger reset). I-016: this entry is
        # radical_class="eliminating" (conduit-deferred), so it re-registers
        # under the non-blocking CONDUIT_DEFERRED bucket, not FEATURE_RADICAL.
        _warn_once_refused(entry)
        assert lookup_candidate(key)["effective_bucket"] == "CONDUIT_DEFERRED"
        reset_conduit_state()

    def test_reset_census_registry_is_alias(self):
        from rmgpy import polymer_conduit as pc
        from rmgpy.polymer import _refused_census_warned, _warn_once_refused
        label = "ALIAS_R(1) <=> ALIAS_Q(2)"
        _warn_once_refused({"reaction": label,
                            "radical_class": "eliminating",
                            "reason": "conduit-deferred"})
        pc.reset_census_registry()
        assert pc.lookup_candidate(
            pc.candidate_key_from_label(label)) is None
        assert (label, "conduit-deferred") not in _refused_census_warned

    def test_reset_clears_flux_totals(self):
        from rmgpy import polymer_conduit as pc
        pc._CONDUIT_FLUX_TOTALS["k"] = {"grams": 1.0, "revoked": False}
        pc.reset_conduit_state()
        assert pc.get_conduit_flux_totals() == {}

    def test_polymer_reexports_reset(self):
        import rmgpy.polymer as rp
        from rmgpy import polymer_conduit as pc
        assert rp.reset_conduit_state is pc.reset_conduit_state


# ---------------------------------------------------------------------------
# Latent classifier-census defects (Codex review, 2026-07): novolac pools
# mis-classified as non-POOL, and near-threshold non-PF chains mis-bucketed
# because chain-scale thresholds were hardcoded to the PF monomer.
# ---------------------------------------------------------------------------

import sys
import types

from rmgpy.polymer_conduit import (
    CHAIN,
    CHAIN_SCALE_HEAVY,
    CHAIN_SCALE_MW,
    DISC,
    POOL,
    chain_scale_thresholds_for_pools,
)


def _cs_species(label, mw, heavy, index=1, formula=None, token=None):
    """Minimal pure-core species record (no live RMG object)."""
    return {
        "token": token or f"{label}({index})",
        "label": label,
        "index": index,
        "formula": formula,
        "mw": mw,
        "heavy_atoms": heavy,
        "link_marker": False,
    }


def _cs_shape_a_record(chain_sp, gas_sp, pool_sp):
    """A shape-A row: [CHAIN] -> [DISC gas, POOL] (reversible source)."""
    return {
        "census": "r93_general",
        "reaction": f"{chain_sp['token']} <=> {gas_sp['token']} + {pool_sp['token']}",
        "reversible": True,
        "reactants": [chain_sp],
        "products": [gas_sp, pool_sp],
        "log_reason": "",
    }


class TestNovolacPoolClassification:
    """BUG 1: a live pool labeled ``novolac`` (and its ``novolac_mod*``
    daughters) must classify as POOL. Before the fix POOL_LABEL_PREFIXES only
    listed ``phenol_formaldehyde``, so ``is_pool`` returned False for every
    novolac pool -- the 329x classifier-not-admissible census-denial defect."""

    def test_novolac_and_daughters_are_pool(self):
        for label in ("novolac", "novolac_mod", "novolac_mod_2",
                      "novolac_mod_3", "novolac_mod_4", "novolac_mu1"):
            sp = _cs_species(label, mw=404.55, heavy=30)
            assert is_pool(sp), f"{label} should be recognized as a pool"
            # A pool participant is POOL regardless of its (chain-sized) MW.
            assert species_role(sp) == POOL, f"{label} role should be POOL"

    def test_phenol_formaldehyde_still_pool(self):
        # Regression guard: the PF prefix is untouched.
        assert is_pool(_cs_species("phenol_formaldehyde", 404.55, 30))
        assert is_pool(_cs_species("phenol_formaldehyde_mod_2", 404.55, 30))

    def test_novolac_shape_a_row_is_admissible_not_unclassified(self):
        # Chain reactant is chain-scale; product side is a small gas + the
        # novolac pool. With novolac recognized as POOL the shape is A and the
        # row is ADMISSIBLE_A. Before the fix the "novolac" product was scored
        # by size (chain-scale -> CHAIN), yielding a (CHAIN, DISC) product side
        # that is outside the shape vocabulary -> UNCLASSIFIED.
        chain = _cs_species("bigchain", mw=600.0, heavy=44, index=100)
        gas = _cs_species("CO", mw=28.01, heavy=2, index=8, formula="CO1")
        pool = _cs_species("novolac", mw=404.55, heavy=30, index=2)
        result = classify_record(_cs_shape_a_record(chain, gas, pool))
        assert result["bucket"] == "ADMISSIBLE_A"
        assert result["destination_pool"] == "novolac"


class TestNovolacChainScaleThresholds:
    """BUG 2: chain-scale role must be decided against the ROW'S OWN monomer,
    not the hardcoded PF monomer (MW 134.178 / 10 heavy atoms). The novolac
    monomer is MW 120.148 / 9 heavy atoms, so a chain of MW ~315 / 23 heavy is
    chain-scale for novolac (>= 2.5*120.148 and >= 2.5*9) but sub-threshold for
    PF (< 2.5*134.178). Before the fix such a near-threshold novolac chain was
    silently mis-bucketed as DISC."""

    # Live novolac monomer metadata (per the smoke-run polymer_pools.json).
    NOVOLAC_MW = 120.148
    NOVOLAC_HEAVY = 9
    CS_MW = 2.5 * 120.148    # 300.37
    CS_HEAVY = 2.5 * 9       # 22.5

    def _near_threshold_chain(self):
        # Between the novolac (300.37 / 22.5) and PF (335.445 / 25) bars.
        # A discrete chain fragment (NOT a pool-labeled species).
        return _cs_species("C=Cc1ccc(CC)cc1O_dimer", mw=315.0, heavy=23, index=200)

    def test_pf_default_thresholds_misbucket_novolac_chain_as_disc(self):
        # Documents the defect: with PF defaults the novolac chain is DISC.
        sp = self._near_threshold_chain()
        assert not is_chain_scale(sp)  # PF default: 315 < 335.445
        assert species_role(sp) == DISC

    def test_live_novolac_thresholds_make_it_chain(self):
        sp = self._near_threshold_chain()
        assert is_chain_scale(sp, self.CS_MW, self.CS_HEAVY)
        assert species_role(sp, self.CS_MW, self.CS_HEAVY) == CHAIN

    def test_classify_record_uses_live_thresholds(self):
        # Same shape-A row classified with PF defaults vs novolac thresholds.
        chain = self._near_threshold_chain()
        gas = _cs_species("CO", mw=28.01, heavy=2, index=8, formula="CO1")
        pool = _cs_species("phenol_formaldehyde", mw=404.55, heavy=30, index=2)
        record = _cs_shape_a_record(chain, gas, pool)

        # PF defaults: chain reads DISC -> (DISC, POOL) product side with a
        # DISC reactant side is outside the shape vocabulary -> UNCLASSIFIED.
        pf = classify_record(record)
        assert pf["bucket"] == "UNCLASSIFIED"

        # Live novolac thresholds: chain reads CHAIN -> shape A -> ADMISSIBLE_A.
        nv = classify_record(record, chain_scale_mw=self.CS_MW,
                             chain_scale_heavy=self.CS_HEAVY)
        assert nv["shape"] == "A"
        assert nv["bucket"] == "ADMISSIBLE_A"

    def test_pf_row_thresholds_unchanged_by_derivation(self):
        # Regression: deriving from a PF monomer reproduces the old constants.
        assert self.CS_MW != CHAIN_SCALE_MW  # novolac differs from PF
        pf_mw = 2.5 * MONOMER_MW
        assert pf_mw == CHAIN_SCALE_MW == 335.445


class TestChainScaleThresholdsForPools:
    """The adapter helper that supplies the live per-row thresholds. Uses a
    stub ``rmgpy.polymer._heavy_atom_count`` so it runs without a Cython
    build."""

    def _install_fake_polymer(self, monkeypatch):
        fake = types.ModuleType("rmgpy.polymer")
        fake._heavy_atom_count = lambda mol: (
            0 if mol is None else int(getattr(mol, "heavy", 0)))
        monkeypatch.setitem(sys.modules, "rmgpy.polymer", fake)

    def _pool(self, mw, heavy):
        return types.SimpleNamespace(
            monomer_mw_g_mol=mw,
            monomer=types.SimpleNamespace(heavy=heavy))

    def test_derives_from_live_monomer(self, monkeypatch):
        self._install_fake_polymer(monkeypatch)
        cs_mw, cs_heavy = chain_scale_thresholds_for_pools(
            [self._pool(120.148, 9)])
        assert cs_mw == pytest.approx(300.37)
        assert cs_heavy == pytest.approx(22.5)

    def test_min_over_multiple_pools(self, monkeypatch):
        self._install_fake_polymer(monkeypatch)
        cs_mw, cs_heavy = chain_scale_thresholds_for_pools(
            [self._pool(134.178, 10), self._pool(120.148, 9)])
        # min() over pools -> loosest (smallest) bar.
        assert cs_mw == pytest.approx(300.37)
        assert cs_heavy == pytest.approx(22.5)

    def test_falls_back_to_pf_defaults(self, monkeypatch):
        self._install_fake_polymer(monkeypatch)
        cs_mw, cs_heavy = chain_scale_thresholds_for_pools([])
        assert cs_mw == CHAIN_SCALE_MW
        assert cs_heavy == CHAIN_SCALE_HEAVY

    def test_pf_pool_reproduces_old_constants(self, monkeypatch):
        self._install_fake_polymer(monkeypatch)
        cs_mw, cs_heavy = chain_scale_thresholds_for_pools(
            [self._pool(134.178, 10)])
        assert cs_mw == pytest.approx(CHAIN_SCALE_MW)
        assert cs_heavy == pytest.approx(CHAIN_SCALE_HEAVY)
