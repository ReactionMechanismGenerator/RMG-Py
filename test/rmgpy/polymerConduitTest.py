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
    "E_row": "DEFERRED_D",
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

    def test_chip_resolves_to_pool_state(self, clean_registry):
        # round-53 REWRITE (was: "trimer_rad* chips are chain-SCALE by size
        # but take the POOL role -- the empirical D/E split"). That pinned
        # the OLD bug: a DECLARED label alone (trimer_rad33 in
        # POOL_STATE_RESOLVABLE_LABELS) was trusted as ground truth with no
        # isomorphism validation ever run. A real conduit4 run showed this
        # exact label is 108/108 a core discrete radical, never isomorphic
        # to any pool -- so a DECLARED label that was never validated
        # (clean_registry-reset process state, no adapter call touched it)
        # must NOT resolve to POOL: is_pool_state_resolvable now consults
        # the VALIDATED ACTIVE set, which starts empty. Chain-scale size
        # alone therefore wins, correctly reclassifying this chip as CHAIN
        # (see TestRealRowBuckets::test_bucket[E_row], now DEFERRED_D).
        chip = E_row["reactants"][0]
        assert chip["label"].startswith("trimer_rad")
        assert not is_pool_state_resolvable(chip)
        assert is_chain_scale(chip)
        assert species_role(chip) == "CHAIN"

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
    def test_bucket(self, name, clean_registry):
        # clean_registry (round-53): classify_record()'s pure-core role
        # test for a DECLARED-but-unvalidated label (trimer_rad33, in
        # E_row) now reads the process-wide ACTIVE label set, so this
        # parametrized test needs pristine oracle state regardless of
        # execution order -- otherwise an earlier adapter-touching test
        # that validated trimer_rad33 ACTIVE could flip E_row's shape.
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

    def test_conduit_echo_record_classifies_to_conduit_echo_bucket(self):
        # P2a fix: classify_record() must honor the documented
        # record["census"] == "conduit_echo" contract (module docstring)
        # the same way it already honors "feature_radical" -- a pure-core
        # caller must never see this fall through to shape/UNCLASSIFIED.
        from rmgpy.polymer_conduit import REASON_CONDUIT_ECHO

        echo_record = dict(FR_elim, census="conduit_echo")
        r = classify_record(echo_record)
        assert r["bucket"] == "CONDUIT_ECHO"
        assert r["refusal_reason"] == REASON_CONDUIT_ECHO
        assert r["admissible"] is False


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

    def test_conduit_echo_ranks_lowest_r93_then_echo(self, clean_registry):
        """round-50 determinism finding: 'conduit_echo' is the LOWEST
        precedence rank -- r93_general's bucket must survive a LATER
        conduit_echo sighting of the same key (echo never overwrites a
        genuine census's effective_bucket)."""
        register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A")
        entry = register_candidate(self.KEY, "conduit_echo", "CONDUIT_ECHO")
        assert entry["effective_bucket"] == "ADMISSIBLE_A"
        assert entry["censuses"] == {"r93_general", "conduit_echo"}

    def test_conduit_echo_ranks_lowest_echo_then_r93(self, clean_registry):
        """Same pair, OPPOSITE registration order: conduit_echo sighted
        FIRST must not stick -- a later r93_general sighting must still
        win for effective_bucket. Together with the r93-then-echo test
        above, this pins order-independence: the SAME sighting set must
        resolve to the SAME ledger state regardless of write order (the
        determinism bug this finding fixes)."""
        register_candidate(self.KEY, "conduit_echo", "CONDUIT_ECHO")
        entry = register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A")
        assert entry["effective_bucket"] == "ADMISSIBLE_A"
        assert entry["censuses"] == {"r93_general", "conduit_echo"}

    def test_conduit_echo_order_independence_end_state_identical(self):
        """Direct order-independence pin: registering the SAME two
        sightings in either order must leave the ledger in an IDENTICAL
        final state (effective_bucket, censuses, bucket_by_census) --
        not merely the same effective_bucket value in isolation."""
        from rmgpy.polymer_conduit import reset_conduit_state

        reset_conduit_state()
        register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A")
        register_candidate(self.KEY, "conduit_echo", "CONDUIT_ECHO")
        forward_order = lookup_candidate(self.KEY)
        reset_conduit_state()

        register_candidate(self.KEY, "conduit_echo", "CONDUIT_ECHO")
        register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A")
        reverse_order = lookup_candidate(self.KEY)
        reset_conduit_state()

        assert forward_order == reverse_order
        assert forward_order["effective_bucket"] == "ADMISSIBLE_A"

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


class TestSameCensusResightResolution:
    """round-50 P1 (main determinism finding fix): a (candidate, census)
    pair sighted with more than one DISTINCT bucket must resolve to a
    single deterministic bucket -- via BUCKET_DECLARATION_ORDER's
    most-conservative-wins rule -- independent of sighting order or which
    rebuild epoch each sighting landed in. Previously this was
    last-write-wins, so the SAME multiset of sightings in different
    orders produced DIFFERENT effective_bucket / summary output."""

    KEY = "RESIGHT_chain(1)<>RESIGHT_gas(2)+RESIGHT_pool(3)"

    def _register_permutation(self, sightings):
        """Register a (census, bucket, epoch) sighting sequence into a
        freshly reset registry; return (final_entry, summary_text)."""
        reset_census_registry()
        entry = None
        for census, bucket, epoch in sightings:
            entry = register_candidate(self.KEY, census, bucket, epoch=epoch)
        return entry, census_summary()

    def test_same_census_divergent_resight_is_order_independent(self):
        # Same multiset of sightings -- r93_general sighted TWICE with
        # DIFFERENT buckets (ADMISSIBLE_A then DEFERRED_A, or the reverse),
        # plus an unrelated conduit_echo sighting -- applied in 3 distinct
        # orders. One sighting lands under "epoch-1", the other under
        # "epoch-2" (the module's existing rebuild-epoch bookkeeping
        # mechanism), so every permutation exercises sightings straddling
        # an epoch boundary regardless of registration order.
        sightings_orders = [
            [("r93_general", "ADMISSIBLE_A", "epoch-1"),
             ("conduit_echo", "CONDUIT_ECHO", "epoch-1"),
             ("r93_general", "DEFERRED_A", "epoch-2")],
            [("r93_general", "DEFERRED_A", "epoch-2"),
             ("conduit_echo", "CONDUIT_ECHO", "epoch-1"),
             ("r93_general", "ADMISSIBLE_A", "epoch-1")],
            [("conduit_echo", "CONDUIT_ECHO", "epoch-1"),
             ("r93_general", "DEFERRED_A", "epoch-2"),
             ("r93_general", "ADMISSIBLE_A", "epoch-1")],
        ]
        results = [self._register_permutation(order)
                  for order in sightings_orders]
        entries = [e for e, _ in results]
        summaries = [s for _, s in results]

        # DEFERRED_A is strictly more conservative than ADMISSIBLE_A in
        # BUCKET_DECLARATION_ORDER, so it must win regardless of order.
        for entry in entries:
            assert entry["effective_bucket"] == "DEFERRED_A"
            assert entry["bucket_by_census"]["r93_general"] == "DEFERRED_A"
            assert entry["bucket_sightings_by_census"]["r93_general"] == {
                "ADMISSIBLE_A", "DEFERRED_A"}
            assert entry["epochs"] == {"epoch-1", "epoch-2"}

        # Full entries (not merely effective_bucket in isolation) are
        # byte-identical across all permutations.
        for entry in entries[1:]:
            assert entry == entries[0]

        # Summary text -- including the new resight_divergence token -- is
        # byte-identical across all permutations too.
        for s in summaries[1:]:
            assert s == summaries[0]
        assert "resight_divergence=1" in summaries[0]

        reset_census_registry()

    def test_no_divergence_control_same_bucket_repeated(self, clean_registry):
        # Re-sighting the SAME census with the SAME bucket repeatedly must
        # never count as a divergence.
        register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A")
        register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A")
        entry = register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A")
        assert entry["bucket_sightings_by_census"]["r93_general"] == {
            "ADMISSIBLE_A"}
        assert entry["effective_bucket"] == "ADMISSIBLE_A"
        assert "resight_divergence=0" in census_summary()

    def test_unclassified_total_recomputed_from_final_state(self,
                                                             clean_registry):
        # round-50 P1 fix #5: a candidate first sighted UNCLASSIFIED under
        # r93_general, then re-sighted (same census) as ADMISSIBLE_A, must
        # NOT leave a stale UNCLASSIFIED contribution -- unclassified_total
        # is computed fresh from final registry state, so it tracks the
        # RESOLVED (most-conservative) bucket, not the first sighting.
        from rmgpy.polymer_conduit import CENSUS_REGISTRY
        register_candidate(self.KEY, "r93_general", "UNCLASSIFIED")
        assert CENSUS_REGISTRY.unclassified_total == 1
        register_candidate(self.KEY, "r93_general", "ADMISSIBLE_A")
        # UNCLASSIFIED is the most conservative bucket in
        # BUCKET_DECLARATION_ORDER, so it still wins the per-census
        # resolution and the candidate's effective_bucket stays
        # UNCLASSIFIED -- this pins the resolution rule, not staleness.
        assert CENSUS_REGISTRY.unclassified_total == 1
        assert "unclassified=1" in census_suffix(
            {"candidate_key": self.KEY, "shape": None, "admissible": False},
            lookup_candidate(self.KEY))


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

        # round-50 FR-census scoping: only a GENUINE reason (qssa-invalid /
        # qssa-unassessable) registers into the real 'feature_radical'
        # census -- see TestFrCensusScoping for the conduit-deferred
        # (echo, non-genuine) counterpart of this test.
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

        # round-50 FR-census scoping: a genuine reason (qssa-invalid) is
        # what makes this an upstream FEATURE-RADICAL sighting; a
        # conduit-deferred echo of the SAME row must NOT do this (see
        # TestFrCensusScoping.test_r93_self_echo_never_creates_fr_membership).
        rxn = self._refused_row()
        label = _reaction_census_label(rxn)
        from rmgpy.polymer import _refused_census_warned
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
        """round-53 REWRITE: with the DECLARED/ACTIVE label-oracle split, a
        label's FIRST sighting is validated against THAT SAME row's own
        structural verdict (:meth:`_LabelOracleState.note_sighting`), so a
        same-row DECLARED-but-never-validated label can no longer produce a
        spurious divergence -- that was exactly the eliminated false
        positive (the OLD form of this test pinned it: an unvalidated
        DECLARED label alone was enough to "diverge" from a structurally
        unrelated species). A genuine divergence is now MID-RUN DRIFT:
        validate trimer_rad33 ACTIVE against a pool it genuinely resolves
        to, then show a LATER row where a species carries that same
        (now-ACTIVE) label but is NOT that pool -- the ACTIVE-set label
        test and this row's own isomorphism verdict disagree, and the
        CONDUIT CLASSIFIER DIVERGENCE/1 line fires with the isomorphism
        verdict winning the role assignment."""
        import logging as _logging

        from rmgpy.molecule import Molecule
        from rmgpy.polymer import Polymer
        from rmgpy.polymer_conduit import (active_pool_state_labels,
                                           record_from_reaction)
        from rmgpy.reaction import Reaction
        from rmgpy.species import Species

        pp = Polymer(label="phenol_formaldehyde", monomer="[CH2][CH]C",
                     Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        # First sighting: a species genuinely isomorphic to pp -> validates
        # trimer_rad33 ACTIVE.
        proxy = Species(molecule=[Molecule().from_smiles(
            pp.molecule[0].to_smiles())])
        proxy.label = "trimer_rad33"
        record_from_reaction(
            Reaction(reactants=[proxy], products=[pp], reversible=True),
            [pp])
        assert "trimer_rad33" in active_pool_state_labels()

        # Later sighting: same ACTIVE label, but THIS species is not pp --
        # the ACTIVE-set label test (True) and this row's own isomorphism
        # verdict (False) disagree.
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
                        pool_mw=True, kinetics=True, label="CHAIN"):
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
        """round-53 REWRITE: a label/isomorphism divergence is a finding,
        never an admission basis. With the DECLARED/ACTIVE label-oracle
        split, trimer_rad33 must first be validated ACTIVE (a same-row
        DECLARED-but-unvalidated label can no longer disagree with its own
        row's isomorphism verdict -- see round-53 note on
        test_label_isomorphism_divergence_is_logged_not_overridden). Seed
        ACTIVE via a genuinely-isomorphic proxy, then relabel the gas
        product of a SEPARATE admissible reaction as the now-ACTIVE
        trimer_rad33 -- the gas is structurally CH2O, not the pf pool, so
        the ACTIVE-set label test and this row's isomorphism verdict
        disagree: genuine mid-run drift."""
        from rmgpy.molecule import Molecule
        from rmgpy.polymer_conduit import (active_pool_state_labels,
                                           evaluate_conduit_admission,
                                           record_from_reaction)
        from rmgpy.reaction import Reaction
        from rmgpy.species import Species

        rxn, pf = _admissible_fixture()

        proxy = Species(molecule=[Molecule().from_smiles(
            pf.molecule[0].to_smiles())])
        proxy.label = "trimer_rad33"
        record_from_reaction(
            Reaction(reactants=[proxy], products=[pf], reversible=True),
            [pf])
        assert "trimer_rad33" in active_pool_state_labels()

        gas = rxn.products[1]
        gas.label = "trimer_rad33"   # ACTIVE label; structure says CH2O
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

    def test_non_chain_scale_chain_denies_condensed(self, clean_registry):
        """The consumed species must be melt-classified (chain-scale
        proxy-derived) so the consumer's phase gate passes the event."""
        from rmgpy.polymer_conduit import evaluate_conduit_admission
        rxn, pf = _admissible_fixture()
        # A huge pool monomer MW pushes the pool-relative chain-scale bar
        # above the chain while the record-level CHAIN role (module
        # constants) still holds.
        pf.monomer_mw_g_mol = 400.0
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        # chain 464.6 g/mol < 2.5*400: fails the pool-relative dual-axis
        # chain-scale conjunct -> chain-not-condensed
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

    def test_stage_token_partitions_the_closed_deny_vocabulary(
            self, clean_registry):
        """(round-49) every census line carries exactly one stage= token,
        adjacent to the would_admit=/deny= tokens: stage=provisional iff
        the verdict is a PROVISIONAL_DENY_REASONS member (subject to G6
        re-adjudication), stage=final for every other deny AND for the
        admit line -- so grep tallies filtered on stage=final never
        double-count a re-adjudicated candidate."""
        from rmgpy.polymer_conduit import (ADMISSION_DENY_REASONS,
                                           PROVISIONAL_DENY_REASONS,
                                           AdmissionVerdict,
                                           admission_census_suffix, _deny,
                                           register_candidate)
        # round-55 P2(a): echo-not-evaluated now REQUIRES a recorded echo
        # sighting for the key (unregistered keys raise), so give 'K' the
        # echo-only ledger entry the token documents before the sweep.
        register_candidate("K", "conduit_echo", "CONDUIT_ECHO")
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

    def test_echo_not_evaluated_unregistered_key_anomaly_fail_closed(
            self, caplog, clean_registry):
        # round-56 F3 (flipped from the round-55 raise pin, itself never
        # pushed): an unregistered-key 'echo-not-evaluated' token is a
        # caller bug, but admission_census_suffix is a census-only
        # chokepoint that must NEVER raise (a raise here was swallowed by
        # the production wrappers into a generic 'annotation-failed' line,
        # losing the invariant info). It now emits a loud VERSIONED anomaly
        # line and falls back to a fail-closed deny that does NOT propagate
        # the reserved token.
        import logging as _logging
        from rmgpy.polymer_conduit import admission_census_suffix, _deny
        with caplog.at_level(_logging.WARNING):
            suffix = admission_census_suffix(
                _deny("NEVER_SEEN(1)<>NOWHERE(2)", "echo-not-evaluated"))
        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 1
        assert "key=NEVER_SEEN(1)<>NOWHERE(2)" in anomalies[0]
        assert "token=echo-not-evaluated" in anomalies[0]
        assert "reason=unregistered-key" in anomalies[0]
        # fail-closed deny, reserved token NOT propagated into the suffix
        assert "would_admit=0" in suffix
        assert "deny=admission-token-anomaly stage=final" in suffix
        assert "echo-not-evaluated" not in suffix

    def test_echo_not_evaluated_mixed_membership_anomaly_fail_closed(
            self, caplog, clean_registry):
        # round-56 F3: 'echo-not-evaluated' is reserved for candidates whose
        # ONLY census membership is conduit_echo. A candidate that ALSO
        # carries a non-echo census membership (e.g. r93_general) WAS
        # evaluated under a real census, so emitting 'not-evaluated' would
        # be a lie. Previously this raised ValueError; it is now kept
        # EQUALLY loud/versioned (not downgraded) as an anomaly line, with a
        # fail-closed deny that does not propagate the reserved token.
        import logging as _logging
        from rmgpy.polymer_conduit import (register_candidate,
                                           admission_census_suffix, _deny)
        key = "MIXED_chain(1)<>MIXED_pool(2)+MIXED_gas(3)"
        register_candidate(key, "r93_general", "ADMISSIBLE_A")
        register_candidate(key, "conduit_echo", "CONDUIT_ECHO")
        with caplog.at_level(_logging.WARNING):
            suffix = admission_census_suffix(_deny(key, "echo-not-evaluated"))
        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 1
        # the anomaly line names the key and its actual non-echo membership
        assert f"key={key}" in anomalies[0]
        assert "reason=mixed-census-membership" in anomalies[0]
        assert "r93_general" in anomalies[0]
        # fail-closed deny; reserved token never propagated
        assert "deny=admission-token-anomaly stage=final" in suffix
        assert "echo-not-evaluated" not in suffix

    def test_echo_not_evaluated_passes_for_echo_only_candidate(
            self, clean_registry):
        # Positive direction: a candidate registered ONLY under
        # conduit_echo is exactly the case the reserved deny reason
        # documents, so admission_census_suffix must emit the token
        # normally (never raise) and must preserve stage=final.
        from rmgpy.polymer_conduit import (register_candidate,
                                           admission_census_suffix, _deny)
        key = "PUREECHO_chain(1)<>PUREECHO_pool(2)+PUREECHO_gas(3)"
        register_candidate(key, "conduit_echo", "CONDUIT_ECHO")
        suffix = admission_census_suffix(_deny(key, "echo-not-evaluated"))
        assert "deny=echo-not-evaluated stage=final" in suffix

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

        # round-50 FR-census scoping: a genuine reason is required for the
        # FR line to self-deny -- see TestFrCensusScoping for the
        # conduit-deferred (echo) counterpart, which must NOT self-deny.
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

    def test_fr_sticky_across_epochs_scoped_by_reason(self, clean_registry):
        """round-50 FR-census scoping: the warn-once hook's SCOPED
        stickiness -- a genuine reason (qssa-invalid / qssa-unassessable)
        sticks as a real 'feature_radical' sighting and blocks admission
        across epochs exactly like test_fr_sticky_across_epochs above; a
        conduit-deferred sighting of a DIFFERENT row echoes into
        CONDUIT_ECHO instead and G1 ignores it -- it must never block
        admission of its own (or any other) candidate key."""
        from rmgpy.polymer import _reaction_census_label, _warn_once_refused
        from rmgpy.polymer_conduit import (candidate_key_from_label,
                                           evaluate_conduit_admission,
                                           lookup_candidate)

        genuine_rxn, genuine_pf = _admissible_fixture()
        genuine_label = _reaction_census_label(genuine_rxn)
        _warn_once_refused({"reaction": genuine_label,
                            "radical_class": "accumulating",
                            "reason": "qssa-invalid"})
        genuine_key = candidate_key_from_label(genuine_label)
        assert lookup_candidate(genuine_key)["effective_bucket"] == (
            "FEATURE_RADICAL")
        v_genuine = evaluate_conduit_admission(genuine_rxn, [genuine_pf])
        assert v_genuine.admitted is False
        assert v_genuine.deny_reason == "feature-radical-overlap"

        # Distinct species pair (own chain label) so echo_key != genuine_key
        # -- otherwise the echo sighting would land on the genuine row's
        # key, which already carries a genuine feature_radical sighting.
        echo_rxn, echo_pf = _admissible_fixture(reversible=True,
                                                 label="ECHO_CHAIN")
        echo_label = _reaction_census_label(echo_rxn)
        _warn_once_refused({"reaction": echo_label,
                            "radical_class": "eliminating",
                            "reason": "conduit-deferred"})
        echo_key = candidate_key_from_label(echo_label)
        assert lookup_candidate(echo_key)["effective_bucket"] == (
            "CONDUIT_ECHO")
        v_echo = evaluate_conduit_admission(echo_rxn, [echo_pf])
        assert v_echo.admitted is True
        assert v_echo.deny_reason is None

    def test_fr_and_echo_censuses_on_same_key_stay_feature_radical(
            self, clean_registry):
        """Mixed-membership case: a candidate key that accumulates BOTH a
        genuine feature_radical sighting AND a conduit_echo sighting must
        stay FEATURE_RADICAL with precedence='feature_radical' -- the
        upstream blocker wins regardless of sighting order. This is the
        exact scenario the old blanket (unscoped) warn-once bug would have
        gotten wrong by downgrading/echoing over a genuine FR row."""
        from rmgpy.polymer import _reaction_census_label, _warn_once_refused
        from rmgpy.polymer_conduit import (candidate_key_from_label,
                                           lookup_candidate)

        rxn, _pf = _admissible_fixture()
        label = _reaction_census_label(rxn)
        key = candidate_key_from_label(label)

        _warn_once_refused({"reaction": label,
                            "radical_class": "accumulating",
                            "reason": "qssa-invalid"})
        entry = lookup_candidate(key)
        assert entry["censuses"] == {"feature_radical"}
        assert entry["effective_bucket"] == "FEATURE_RADICAL"
        assert entry["precedence"] is None

        _warn_once_refused({"reaction": label,
                            "radical_class": "eliminating",
                            "reason": "conduit-deferred"})
        entry = lookup_candidate(key)
        assert entry["censuses"] == {"feature_radical", "conduit_echo"}
        assert entry["effective_bucket"] == "FEATURE_RADICAL"
        assert entry["precedence"] == "feature_radical"

    def test_fr_and_echo_censuses_on_same_key_stay_feature_radical_reverse_order(
            self, clean_registry):
        """Same mixed-membership scenario as
        :meth:`test_fr_and_echo_censuses_on_same_key_stay_feature_radical`,
        with sighting order REVERSED (echo first, genuine FR second) --
        the ledger must land in the IDENTICAL final state either way: the
        upstream blocker wins regardless of which sighting arrives first."""
        from rmgpy.polymer import _reaction_census_label, _warn_once_refused
        from rmgpy.polymer_conduit import (candidate_key_from_label,
                                           lookup_candidate)

        rxn, _pf = _admissible_fixture()
        label = _reaction_census_label(rxn)
        key = candidate_key_from_label(label)

        _warn_once_refused({"reaction": label,
                            "radical_class": "eliminating",
                            "reason": "conduit-deferred"})
        entry = lookup_candidate(key)
        assert entry["censuses"] == {"conduit_echo"}
        assert entry["effective_bucket"] == "CONDUIT_ECHO"
        assert entry["precedence"] is None

        _warn_once_refused({"reaction": label,
                            "radical_class": "accumulating",
                            "reason": "qssa-invalid"})
        entry = lookup_candidate(key)
        assert entry["censuses"] == {"feature_radical", "conduit_echo"}
        assert entry["effective_bucket"] == "FEATURE_RADICAL"
        assert entry["precedence"] == "feature_radical"

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
        # round-50 FR-census scoping: 'conduit-deferred' is a non-genuine
        # (echo) reason, so this row lands in CONDUIT_ECHO, not
        # FEATURE_RADICAL -- reset must clear that census entry too.
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
        # warn-once set had survived the ledger reset)
        _warn_once_refused(entry)
        assert lookup_candidate(key)["effective_bucket"] == "CONDUIT_ECHO"
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


class TestFrCensusScoping:
    """round-50 (ratified, adversarial review round 50) FR-census scoping
    fix: a refusal echoing through the warn-once hook (:func:`rmgpy.polymer.
    _warn_once_refused` / :func:`rmgpy.polymer_conduit.annotate_feature_radical`)
    for a NON-genuine reason (chiefly 'conduit-deferred', the accumulating=
    False branch that M18.4 will eventually try to admit) must never join
    the real 'feature_radical' ledger census -- only qssa-invalid /
    qssa-unassessable rows do. CONDUIT_ADMISSION_ENABLED stays False
    throughout: this is a pure census/ledger bookkeeping fix, never a flux
    change."""

    def test_r93_self_echo_never_creates_fr_membership(self, clean_registry):
        """A row refused with 'conduit-deferred' (echoed through the same
        warn-once hook that genuine FR rows go through) must register into
        CONDUIT_ECHO, not FEATURE_RADICAL -- the downgrade trap this split
        fixes: previously a candidate's own refusal echo poisoned its own
        key (and any co-keyed row) as an upstream feature-radical blocker."""
        from rmgpy.polymer import _reaction_census_label, _warn_once_refused
        from rmgpy.polymer_conduit import (candidate_key_from_label,
                                           lookup_candidate)

        rxn, _pf = _admissible_fixture()
        label = _reaction_census_label(rxn)
        _warn_once_refused({"reaction": label,
                            "radical_class": "eliminating",
                            "reason": "conduit-deferred"})
        entry = lookup_candidate(candidate_key_from_label(label))
        assert entry is not None
        assert entry["censuses"] == {"conduit_echo"}
        assert "feature_radical" not in entry["censuses"]
        assert entry["effective_bucket"] == "CONDUIT_ECHO"

    def test_pure_echo_line_carries_admission_verdict_token(
            self, caplog, clean_registry):
        """Contract fix: EVERY census line -- including a pure conduit_echo
        line -- carries the ``[conduit-admission/1 ...]`` verdict token
        (the module docstring/BUILD_SPEC W1.6 promise). A pure echo
        sighting never sees a live reaction object, so it cannot fabricate
        a real G0/G2-G7 verdict; it must stay fail-closed (would_admit=0)
        with the HONEST 'echo-not-evaluated' reason, never the
        'feature-radical-overlap' reason (that would misrepresent a key
        G1 does not actually block as upstream-blocked)."""
        import logging as _logging

        from rmgpy.polymer import _reaction_census_label, _warn_once_refused

        rxn, _pf = _admissible_fixture()
        label = _reaction_census_label(rxn)
        with caplog.at_level(_logging.WARNING):
            _warn_once_refused({"reaction": label,
                                "radical_class": "eliminating",
                                "reason": "conduit-deferred"})
        msgs = [r.getMessage() for r in caplog.records
                if r.getMessage().startswith("FEATURE-RADICAL REFUSED "
                                             "CENSUS: ")]
        assert len(msgs) == 1
        assert ("[conduit-admission/1 would_admit=0 "
                "deny=echo-not-evaluated stage=final]") in msgs[0]
        assert "deny=feature-radical-overlap" not in msgs[0]
        # append-only: the M18.2 census suffix rides before the verdict.
        assert msgs[0].index("[conduit-census/1") < msgs[0].index(
            "[conduit-admission/1")

    def test_mixed_membership_echo_line_still_self_denies_fr_overlap(
            self, caplog, clean_registry):
        """When the SAME key already carries genuine 'feature_radical'
        membership (from an earlier sighting under a genuine reason), a
        LATER pure-echo-reason sighting of that same key must still print
        the real G1 verdict -- deny=feature-radical-overlap -- because G1
        genuinely blocks this key; 'echo-not-evaluated' is reserved for
        keys G1 does NOT block."""
        import logging as _logging

        from rmgpy.polymer import _reaction_census_label, _warn_once_refused

        rxn, _pf = _admissible_fixture()
        label = _reaction_census_label(rxn)
        with caplog.at_level(_logging.WARNING):
            _warn_once_refused({"reaction": label,
                                "radical_class": "accumulating",
                                "reason": "qssa-invalid"})
        caplog.clear()  # drop the genuine-reason line; only the later
                        # echo-reason line is under test here.
        with caplog.at_level(_logging.WARNING):
            _warn_once_refused({"reaction": label,
                                "radical_class": "eliminating",
                                "reason": "conduit-deferred"})
        msgs = [r.getMessage() for r in caplog.records
                if r.getMessage().startswith("FEATURE-RADICAL REFUSED "
                                             "CENSUS: ")]
        assert len(msgs) == 1
        assert ("[conduit-admission/1 would_admit=0 "
                "deny=feature-radical-overlap stage=final]") in msgs[0]
        assert "echo-not-evaluated" not in msgs[0]

    def test_admitted_row_not_downgraded_by_echo_census_entry(
            self, clean_registry):
        """An echo-census sighting (conduit-deferred) of a candidate key
        must never downgrade a WOULD-ADMIT verdict for that same key --
        only a genuine 'feature_radical' sighting may. Exercised at both
        consumer sites: G1 (:func:`evaluate_conduit_admission`) and the r42
        post-registration re-check (:func:`annotate_refused_row`)."""
        from rmgpy.polymer import _reaction_census_label, _warn_once_refused
        from rmgpy.polymer_conduit import (annotate_refused_row,
                                           candidate_key_from_label,
                                           evaluate_conduit_admission)

        rxn, pf = _admissible_fixture()
        label = _reaction_census_label(rxn)
        key = candidate_key_from_label(label)

        # An echo sighting of the SAME key lands first (e.g. from a prior
        # rebuild epoch's accumulating=False branch).
        _warn_once_refused({"reaction": label,
                            "radical_class": "eliminating",
                            "reason": "conduit-deferred"})

        # G1 still admits: the echo census is invisible to it.
        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is True
        assert v.deny_reason is None

        # The r42 post-registration re-check (annotate_refused_row) must
        # likewise refuse to downgrade a would-admit verdict on account of
        # the echo-only sighting.
        suffix = annotate_refused_row(rxn, [pf], census="r93_general",
                                      verdict=v)
        assert "would_admit=1" in suffix
        assert "deny=feature-radical-overlap" not in suffix

    def test_reason_threading_deterministic_across_rebuild_epochs(
            self, clean_registry):
        """Re-running the SAME input rows (same reaction labels, same
        reasons) across successive rebuild epochs must produce IDENTICAL
        census contents -- the warn-once dedup keys on (reaction, reason),
        so re-sighting the identical row in a later epoch is a no-op for
        the ledger, and a fresh key seen with the same reason in a later
        epoch resolves to the same census/bucket every time."""
        from rmgpy.polymer import _reaction_census_label, _warn_once_refused
        from rmgpy.polymer_conduit import (candidate_key_from_label,
                                           lookup_candidate,
                                           reset_conduit_state)

        genuine_rxn, _genuine_pf = _admissible_fixture()
        genuine_label = _reaction_census_label(genuine_rxn)
        # Distinct species pair (own chain label) so echo_key != genuine_key
        # -- otherwise the echo sighting would land on the genuine row's
        # key and the genuine feature_radical sighting would poison it.
        echo_rxn, _echo_pf = _admissible_fixture(reversible=True,
                                                  label="ECHO_CHAIN")
        echo_label = _reaction_census_label(echo_rxn)

        def run_one_rebuild_epoch():
            _warn_once_refused({"reaction": genuine_label,
                                "radical_class": "accumulating",
                                "reason": "qssa-invalid"})
            _warn_once_refused({"reaction": echo_label,
                                "radical_class": "eliminating",
                                "reason": "conduit-deferred"})

        run_one_rebuild_epoch()
        genuine_key = candidate_key_from_label(genuine_label)
        echo_key = candidate_key_from_label(echo_label)
        first_genuine = lookup_candidate(genuine_key)
        first_echo = lookup_candidate(echo_key)
        assert first_genuine["effective_bucket"] == "FEATURE_RADICAL"
        assert first_echo["effective_bucket"] == "CONDUIT_ECHO"

        # A later rebuild epoch re-sighting the SAME rows (warn-once
        # dedup makes this a no-op) must reproduce the identical census
        # state -- not merely an equivalent one.
        run_one_rebuild_epoch()
        second_genuine = lookup_candidate(genuine_key)
        second_echo = lookup_candidate(echo_key)
        assert second_genuine == first_genuine
        assert second_echo == first_echo

        # A hard reset + fresh rebuild epoch (as happens at a new RMG run
        # boundary) must reproduce the SAME census contents from the same
        # input rows -- reason threading is deterministic, not order- or
        # run-dependent.
        reset_conduit_state()
        run_one_rebuild_epoch()
        assert lookup_candidate(genuine_key) == first_genuine
        assert lookup_candidate(echo_key) == first_echo


class TestLabelOracleAcceptanceBattery:
    """round-53 ratified acceptance battery for the DECLARED-vs-ACTIVE
    pool-state label oracle (:class:`rmgpy.polymer_conduit._LabelOracleState`).
    The production side landed with no tests; this class is the missing
    coverage. It pins: the one-shot ``CONDUIT CLASSIFIER ORACLE HEALTH/1``
    line, the per-label ``CONDUIT CLASSIFIER LABEL VALIDATION/1`` line
    (all four reasons), the reshaped versioned runtime
    ``CONDUIT CLASSIFIER DIVERGENCE/1`` line, the boolean
    ``pool_role_override`` crash fix (the old
    ``POOL_STATE_RESOLVABLE_LABELS[0]`` surrogate raised IndexError on an
    empty resolvable set), and finalize() determinism across sighting
    order and rebuild-epoch boundaries."""

    #: Three SMILES that are structurally unrelated to the phenol-
    #: formaldehyde proxy used throughout this module's fixtures -- never
    #: isomorphic to any pool built from it.
    NONISO_SMILES = ["CCO", "CCN", "CC=O"]

    @staticmethod
    def _pf():
        from rmgpy.polymer import Polymer
        return Polymer(label="phenol_formaldehyde", monomer="[CH2][CH]C",
                       Mn=5000.0, Mw=8000.0, initial_mass=1.0)

    @staticmethod
    def _iso_proxy_species(pool, label):
        """A discrete Species whose molecule is exactly the pool's own
        reactive proxy image -- genuinely isomorphic to that pool."""
        from rmgpy.molecule import Molecule
        from rmgpy.species import Species
        sp = Species(molecule=[Molecule().from_smiles(
            pool.molecule[0].to_smiles())])
        sp.label = label
        return sp

    @staticmethod
    def _noniso_species(smiles, label):
        from rmgpy.molecule import Molecule
        from rmgpy.species import Species
        sp = Species(molecule=[Molecule().from_smiles(smiles)])
        sp.label = label
        return sp

    @staticmethod
    def _lines(caplog, prefix):
        return [r.getMessage() for r in caplog.records
                if r.getMessage().startswith(prefix)]

    def _label_lines(self, caplog):
        return self._lines(caplog, "CONDUIT CLASSIFIER LABEL VALIDATION/1")

    def _health_lines(self, caplog):
        return self._lines(caplog, "CONDUIT CLASSIFIER ORACLE HEALTH/1")

    def _divergence_lines(self, caplog):
        return self._lines(caplog, "CONDUIT CLASSIFIER DIVERGENCE/1")

    def test_a_conduit4_shaped_none_resolvable(self, caplog, clean_registry):
        """The real conduit4 shape (round-53's own motivating regression):
        all three declared labels sighted on genuinely non-isomorphic
        discretes -> every one drops as iso-fail (the species WAS seen,
        it just isn't isomorphic to the row's pool -- 'missing-species'
        would be the wrong reason here), health line
        configured=3 active=0 dropped=3, and ZERO runtime divergence: the
        label test and the isomorphism test agree (both say False) once
        the label has been validated, so no deny=classifier-divergence is
        ever produced."""
        import logging as _logging

        from rmgpy.polymer_conduit import (POOL_STATE_RESOLVABLE_LABELS,
                                           evaluate_conduit_admission,
                                           finalize_label_oracle,
                                           record_from_reaction)
        from rmgpy.reaction import Reaction

        pf = self._pf()
        assert len(POOL_STATE_RESOLVABLE_LABELS) == 3
        rows = []
        with caplog.at_level(_logging.WARNING):
            for label, smiles in zip(POOL_STATE_RESOLVABLE_LABELS,
                                     self.NONISO_SMILES):
                sp = self._noniso_species(smiles, label)
                rxn = Reaction(reactants=[sp], products=[pf], reversible=True)
                record = record_from_reaction(rxn, [pf])
                assert record["label_isomorphism_divergence"] is False
                rows.append(rxn)
            # round-55 P1-3: census_summary() is read-only now; the
            # end-of-lifecycle finalization is explicit.
            finalize_label_oracle()

        label_lines = self._label_lines(caplog)
        assert len(label_lines) == 3
        for label, line in zip(POOL_STATE_RESOLVABLE_LABELS, label_lines):
            assert f"label={label}" in line
            assert "status=dropped" in line
            assert "reason=iso-fail" in line

        health_lines = self._health_lines(caplog)
        assert len(health_lines) == 1
        expected_health = (
            "CONDUIT CLASSIFIER ORACLE HEALTH/1 configured=3 active=0 "
            "dropped=3 source=deck validation=lazy-first-sight")
        assert health_lines[0] == expected_health

        assert not self._divergence_lines(caplog)

        for rxn in rows:
            v = evaluate_conduit_admission(rxn, [pf])
            assert v.deny_reason != "classifier-divergence"

    def test_b_positive_control_iso_resolvable(self, caplog, clean_registry):
        """A declared label that IS genuinely iso-resolvable -> LABEL
        VALIDATION status=active reason=iso-pass, health active=1, no
        runtime divergence for that species."""
        import logging as _logging

        from rmgpy.polymer_conduit import (active_pool_state_labels,
                                           finalize_label_oracle,
                                           record_from_reaction)
        from rmgpy.reaction import Reaction

        pf = self._pf()
        proxy = self._iso_proxy_species(pf, "trimer_rad33")
        rxn = Reaction(reactants=[proxy], products=[pf], reversible=True)

        with caplog.at_level(_logging.WARNING):
            record = record_from_reaction(rxn, [pf])
            finalize_label_oracle()

        assert record["label_isomorphism_divergence"] is False
        assert "trimer_rad33" in active_pool_state_labels()

        label_lines = self._label_lines(caplog)
        assert any("label=trimer_rad33" in l and "status=active" in l
                   and "reason=iso-pass" in l for l in label_lines)

        health_lines = self._health_lines(caplog)
        assert len(health_lines) == 1
        assert "configured=3" in health_lines[0]
        assert "active=1" in health_lines[0]

    def test_c_mutation_undeclared_label_diverges_and_denies(
            self, caplog, clean_registry):
        """Mutation check 1: same species shape as (b) (genuinely
        isomorphic) but the label is NOT declared -> the ACTIVE-set label
        test says False while the isomorphism test says True: a genuine
        runtime DIVERGENCE/1 with label_says=0 iso_says=1 action=use-iso.
        The isomorphism verdict still wins the role assignment (the
        species is classified POOL), and admission denies with
        deny=classifier-divergence on the affected row."""
        import logging as _logging

        from rmgpy.polymer_conduit import (POOL_STATE_RESOLVABLE_LABELS,
                                           _apply_iso_overrides,
                                           annotate_refused_row,
                                           evaluate_conduit_admission,
                                           record_from_reaction,
                                           species_role)
        from rmgpy.reaction import Reaction

        pf = self._pf()
        undeclared_label = "trimer_rad_not_declared"
        assert undeclared_label not in POOL_STATE_RESOLVABLE_LABELS
        proxy = self._iso_proxy_species(pf, undeclared_label)
        rxn = Reaction(reactants=[proxy], products=[pf], reversible=True)

        with caplog.at_level(_logging.WARNING):
            record = record_from_reaction(rxn, [pf])
        assert record["label_isomorphism_divergence"] is True
        div_lines = self._divergence_lines(caplog)
        assert len(div_lines) == 1
        assert f"label={undeclared_label}" in div_lines[0]
        assert "label_says=0" in div_lines[0]
        assert "iso_says=1" in div_lines[0]
        assert "active_label=0" in div_lines[0]
        assert "action=use-iso" in div_lines[0]

        # Isomorphism verdict still applied: the species is classified
        # POOL-role via the boolean override, never silently dropped.
        applied = _apply_iso_overrides(copy.deepcopy(record))
        assert applied["reactants"][0]["pool_role_override"] is True
        assert species_role(applied["reactants"][0]) == "POOL"

        v = evaluate_conduit_admission(rxn, [pf])
        assert v.admitted is False
        assert v.deny_reason == "classifier-divergence"

        suffix = annotate_refused_row(rxn, [pf], verdict=v)
        assert "deny=classifier-divergence" in suffix

    def test_d_mutation_bogus_declared_label_dropped_others_unaffected(
            self, caplog, clean_registry):
        """Mutation check 2: alongside a (b)-shaped valid sighting, sight
        ANOTHER declared label on a non-isomorphic species -> that label
        drops (iso-fail: it was seen, just not isomorphic), health
        dropped count increases, and runtime divergence stays 0 for BOTH
        the valid and the bogus parts (each row's label test and
        isomorphism test agree once validated)."""
        import logging as _logging

        from rmgpy.polymer_conduit import (POOL_STATE_RESOLVABLE_LABELS,
                                           active_pool_state_labels,
                                           finalize_label_oracle,
                                           record_from_reaction)
        from rmgpy.reaction import Reaction

        pf = self._pf()
        valid_label = POOL_STATE_RESOLVABLE_LABELS[0]
        bogus_label = POOL_STATE_RESOLVABLE_LABELS[1]
        valid = self._iso_proxy_species(pf, valid_label)
        bogus = self._noniso_species("CCO", bogus_label)

        with caplog.at_level(_logging.WARNING):
            r1 = record_from_reaction(
                Reaction(reactants=[valid], products=[pf], reversible=True),
                [pf])
            r2 = record_from_reaction(
                Reaction(reactants=[bogus], products=[pf], reversible=True),
                [pf])
            finalize_label_oracle()

        assert r1["label_isomorphism_divergence"] is False
        assert r2["label_isomorphism_divergence"] is False
        assert valid_label in active_pool_state_labels()
        assert bogus_label not in active_pool_state_labels()

        label_lines = self._label_lines(caplog)
        assert any(f"label={bogus_label}" in l and "status=dropped" in l
                   and "reason=iso-fail" in l for l in label_lines)

        health_lines = self._health_lines(caplog)
        assert len(health_lines) == 1
        assert "active=1" in health_lines[0]
        # bogus_label (iso-fail) + the third, never-sighted declared
        # label (missing-species) = 2 dropped.
        assert "dropped=2" in health_lines[0]
        assert not self._divergence_lines(caplog)

    def test_e_missing_species_reason(self, caplog, clean_registry):
        """A declared label never sighted at all during the run resolves
        as missing-species at finalize -- there is no species to validate
        it against, so it can never earn iso-pass."""
        import logging as _logging

        from rmgpy.polymer_conduit import (POOL_STATE_RESOLVABLE_LABELS,
                                           finalize_label_oracle)

        with caplog.at_level(_logging.WARNING):
            finalize_label_oracle()

        label_lines = self._label_lines(caplog)
        assert len(label_lines) == len(POOL_STATE_RESOLVABLE_LABELS)
        for label in POOL_STATE_RESOLVABLE_LABELS:
            assert any(
                f"label={label} status=dropped reason=missing-species"
                in l for l in label_lines)

        health_lines = self._health_lines(caplog)
        assert len(health_lines) == 1
        assert "configured=3 active=0 dropped=3" in health_lines[0]

    def test_f_crash_regression_empty_active_set_iso_true(
            self, caplog, clean_registry):
        """Old surrogate-label bug: ``label_for_roles =
        POOL_STATE_RESOLVABLE_LABELS[0]`` raised IndexError whenever the
        resolvable-label state was empty and a species genuinely resolved
        via isomorphism under an undeclared/unvalidated label -- exactly
        the state a fresh process (nothing validated yet) starts in.
        round-53 replaces the surrogate with a boolean
        ``pool_role_override``: this pins that no exception is raised,
        POOL role is applied via the override, and the DIVERGENCE/1 +
        deny=classifier-divergence tokens still fire."""
        import logging as _logging

        from rmgpy.polymer_conduit import (_apply_iso_overrides,
                                           active_pool_state_labels,
                                           annotate_refused_row,
                                           evaluate_conduit_admission,
                                           record_from_reaction,
                                           species_role)
        from rmgpy.reaction import Reaction

        # Fresh clean_registry: the ACTIVE set is genuinely empty -- this
        # is the exact state that used to crash.
        assert active_pool_state_labels() == frozenset()

        pf = self._pf()
        proxy = self._iso_proxy_species(pf, "wholly_undeclared_label")
        rxn = Reaction(reactants=[proxy], products=[pf], reversible=True)

        with caplog.at_level(_logging.WARNING):
            record = record_from_reaction(rxn, [pf])  # must not raise

        assert record["label_isomorphism_divergence"] is True
        applied = _apply_iso_overrides(record)
        assert applied["reactants"][0]["pool_role_override"] is True
        assert species_role(applied["reactants"][0]) == "POOL"  # no IndexError

        v = evaluate_conduit_admission(rxn, [pf])
        assert v.deny_reason == "classifier-divergence"
        suffix = annotate_refused_row(rxn, [pf], verdict=v)
        assert "deny=classifier-divergence" in suffix
        assert self._divergence_lines(caplog)

    def test_ambiguous_label_reason_when_implemented(
            self, caplog, clean_registry):
        """The production code implements an 'ambiguous' validation
        reason (:data:`rmgpy.polymer_conduit.LABEL_VALIDATION_AMBIGUOUS`):
        a declared label whose species is isomorphic to more than one
        distinct pool in the row is dropped as ambiguous rather than
        accepted as active."""
        import logging as _logging

        from rmgpy.polymer import Polymer
        from rmgpy.polymer_conduit import (LABEL_VALIDATION_AMBIGUOUS,
                                           record_from_reaction)
        from rmgpy.reaction import Reaction

        pool1 = Polymer(label="phenol_formaldehyde", monomer="[CH2][CH]C",
                        Mn=5000.0, Mw=8000.0, initial_mass=1.0)
        pool2 = Polymer(label="phenol_formaldehyde_dup",
                        monomer="[CH2][CH]C", Mn=5000.0, Mw=8000.0,
                        initial_mass=1.0)
        proxy = self._iso_proxy_species(pool1, "trimer_rad44")
        rxn = Reaction(reactants=[proxy], products=[pool1], reversible=True)

        with caplog.at_level(_logging.WARNING):
            record_from_reaction(rxn, [pool1, pool2])

        label_lines = self._label_lines(caplog)
        assert any("label=trimer_rad44" in l and "status=dropped" in l
                   and f"reason={LABEL_VALIDATION_AMBIGUOUS}" in l
                   for l in label_lines)

    def test_g_determinism_across_orderings_and_epoch_boundary(
            self, caplog, clean_registry):
        """Identical sighting multiset delivered in different orders, and
        again across a hard reset_conduit_state() (rebuild-epoch)
        boundary, must produce byte-identical final summary text
        including the health line."""
        import logging as _logging

        from rmgpy.polymer_conduit import (POOL_STATE_RESOLVABLE_LABELS,
                                           census_summary,
                                           finalize_label_oracle,
                                           record_from_reaction,
                                           reset_conduit_state)
        from rmgpy.reaction import Reaction

        pf = self._pf()
        valid_label, bogus_label = POOL_STATE_RESOLVABLE_LABELS[:2]

        def run_epoch(order):
            reset_conduit_state()
            sightings = {
                "valid": self._iso_proxy_species(pf, valid_label),
                "bogus": self._noniso_species("CCO", bogus_label),
            }
            caplog.clear()
            with caplog.at_level(_logging.WARNING):
                for key in order:
                    sp = sightings[key]
                    record_from_reaction(
                        Reaction(reactants=[sp], products=[pf],
                                reversible=True),
                        [pf])
                summary = census_summary()
                # round-55 P1-3: the summary is read-only; the health
                # line comes from the explicit lifecycle finalization.
                finalize_label_oracle()
            health = self._health_lines(caplog)
            assert len(health) == 1
            return summary + "\n" + health[0]

        forward = run_epoch(["valid", "bogus"])
        reversed_order = run_epoch(["bogus", "valid"])
        assert forward == reversed_order

        epoch_1 = run_epoch(["valid", "bogus"])
        epoch_2 = run_epoch(["valid", "bogus"])
        assert epoch_1 == epoch_2


# ---------------------------------------------------------------------------
# round-55 M3: REAL conduit4 artifact fixtures (poly_102_conduit4, read-only
# forensic artifact). The three DECLARED pool-state labels' species, extracted
# VERBATIM from /home/alon/runs/RMG/poly_102_conduit4/chemkin/
# species_dictionary.txt (adjacency lists, C27H31O3 trimer radicals), and the
# run's own phenol_formaldehyde pool definition, copied from the deck
# (/home/alon/runs/RMG/poly_102_conduit4/input.py polymer() block). Embedded
# here as literal fixture data -- the artifact itself is never written.
# ---------------------------------------------------------------------------

#: Deck pool parameters (poly_102_conduit4 input.py, polymer() block).
CONDUIT4_POOL_MONOMER_SMILES = "[CH2]c1c(O)c([CH2])c(C)cc1"
CONDUIT4_POOL_KWARGS = dict(
    label="phenol_formaldehyde",
    end_groups=["[H]", "[H]"],
    cutoff=3,
    Mn=1000.0,
    Mw=3000.0,
    initial_mass=1.0,
)

#: chemkin/species_dictionary.txt indices of the three trimer species.
CONDUIT4_TRIMER_INDICES = {"trimer_rad33": 22, "trimer_rad38": 25,
                           "trimer_rad44": 28}

#: Verbatim adjacency-list blocks (header line included, exactly as written
#: in the artifact's species_dictionary.txt).
CONDUIT4_TRIMER_ADJLISTS = {
    "trimer_rad33": """\
trimer_rad33(22)
multiplicity 2
1  O u0 p2 c0 {23,S} {59,S}
2  O u0 p2 c0 {21,S} {60,S}
3  O u0 p2 c0 {22,S} {61,S}
4  C u0 p0 c0 {5,S} {15,S} {31,S} {32,S}
5  C u0 p0 c0 {4,S} {12,S} {33,S} {34,S}
6  C u0 p0 c0 {13,S} {24,S} {35,S} {36,S}
7  C u0 p0 c0 {14,S} {37,S} {38,S} {39,S}
8  C u0 p0 c0 {16,S} {40,S} {41,S} {42,S}
9  C u0 p0 c0 {17,S} {43,S} {44,S} {45,S}
10 C u0 p0 c0 {19,S} {46,S} {47,S} {51,S}
11 C u0 p0 c0 {18,S} {48,S} {49,S} {50,S}
12 C u0 p0 c0 {5,S} {18,B} {22,B}
13 C u0 p0 c0 {6,S} {17,B} {21,B}
14 C u0 p0 c0 {7,S} {16,B} {23,B}
15 C u0 p0 c0 {4,S} {21,B} {27,B}
16 C u0 p0 c0 {8,S} {14,B} {25,B}
17 C u0 p0 c0 {9,S} {13,B} {26,B}
18 C u0 p0 c0 {11,S} {12,B} {28,B}
19 C u0 p0 c0 {10,S} {22,B} {29,B}
20 C u0 p0 c0 {23,B} {24,S} {30,B}
21 C u0 p0 c0 {2,S} {13,B} {15,B}
22 C u0 p0 c0 {3,S} {12,B} {19,B}
23 C u0 p0 c0 {1,S} {14,B} {20,B}
24 C u1 p0 c0 {6,S} {20,S} {52,S}
25 C u0 p0 c0 {16,B} {30,B} {53,S}
26 C u0 p0 c0 {17,B} {27,B} {55,S}
27 C u0 p0 c0 {15,B} {26,B} {56,S}
28 C u0 p0 c0 {18,B} {29,B} {57,S}
29 C u0 p0 c0 {19,B} {28,B} {58,S}
30 C u0 p0 c0 {20,B} {25,B} {54,S}
31 H u0 p0 c0 {4,S}
32 H u0 p0 c0 {4,S}
33 H u0 p0 c0 {5,S}
34 H u0 p0 c0 {5,S}
35 H u0 p0 c0 {6,S}
36 H u0 p0 c0 {6,S}
37 H u0 p0 c0 {7,S}
38 H u0 p0 c0 {7,S}
39 H u0 p0 c0 {7,S}
40 H u0 p0 c0 {8,S}
41 H u0 p0 c0 {8,S}
42 H u0 p0 c0 {8,S}
43 H u0 p0 c0 {9,S}
44 H u0 p0 c0 {9,S}
45 H u0 p0 c0 {9,S}
46 H u0 p0 c0 {10,S}
47 H u0 p0 c0 {10,S}
48 H u0 p0 c0 {11,S}
49 H u0 p0 c0 {11,S}
50 H u0 p0 c0 {11,S}
51 H u0 p0 c0 {10,S}
52 H u0 p0 c0 {24,S}
53 H u0 p0 c0 {25,S}
54 H u0 p0 c0 {30,S}
55 H u0 p0 c0 {26,S}
56 H u0 p0 c0 {27,S}
57 H u0 p0 c0 {28,S}
58 H u0 p0 c0 {29,S}
59 H u0 p0 c0 {1,S}
60 H u0 p0 c0 {2,S}
61 H u0 p0 c0 {3,S}
""",
    "trimer_rad38": """\
trimer_rad38(25)
multiplicity 2
1  O u0 p2 c0 {21,S} {59,S}
2  O u0 p2 c0 {23,S} {60,S}
3  O u0 p2 c0 {22,S} {61,S}
4  C u0 p0 c0 {5,S} {14,S} {31,S} {32,S}
5  C u0 p0 c0 {4,S} {12,S} {33,S} {34,S}
6  C u0 p0 c0 {13,S} {24,S} {35,S} {36,S}
7  C u0 p0 c0 {15,S} {37,S} {38,S} {39,S}
8  C u0 p0 c0 {16,S} {40,S} {41,S} {42,S}
9  C u0 p0 c0 {17,S} {43,S} {44,S} {45,S}
10 C u0 p0 c0 {19,S} {46,S} {47,S} {51,S}
11 C u0 p0 c0 {18,S} {48,S} {49,S} {50,S}
12 C u0 p0 c0 {5,S} {17,B} {23,B}
13 C u0 p0 c0 {6,S} {18,B} {22,B}
14 C u0 p0 c0 {4,S} {21,B} {26,B}
15 C u0 p0 c0 {7,S} {16,B} {21,B}
16 C u0 p0 c0 {8,S} {15,B} {25,B}
17 C u0 p0 c0 {9,S} {12,B} {27,B}
18 C u0 p0 c0 {11,S} {13,B} {28,B}
19 C u0 p0 c0 {10,S} {22,B} {29,B}
20 C u0 p0 c0 {23,B} {24,S} {30,B}
21 C u0 p0 c0 {1,S} {14,B} {15,B}
22 C u0 p0 c0 {3,S} {13,B} {19,B}
23 C u0 p0 c0 {2,S} {12,B} {20,B}
24 C u1 p0 c0 {6,S} {20,S} {52,S}
25 C u0 p0 c0 {16,B} {26,B} {53,S}
26 C u0 p0 c0 {14,B} {25,B} {54,S}
27 C u0 p0 c0 {17,B} {30,B} {55,S}
28 C u0 p0 c0 {18,B} {29,B} {57,S}
29 C u0 p0 c0 {19,B} {28,B} {58,S}
30 C u0 p0 c0 {20,B} {27,B} {56,S}
31 H u0 p0 c0 {4,S}
32 H u0 p0 c0 {4,S}
33 H u0 p0 c0 {5,S}
34 H u0 p0 c0 {5,S}
35 H u0 p0 c0 {6,S}
36 H u0 p0 c0 {6,S}
37 H u0 p0 c0 {7,S}
38 H u0 p0 c0 {7,S}
39 H u0 p0 c0 {7,S}
40 H u0 p0 c0 {8,S}
41 H u0 p0 c0 {8,S}
42 H u0 p0 c0 {8,S}
43 H u0 p0 c0 {9,S}
44 H u0 p0 c0 {9,S}
45 H u0 p0 c0 {9,S}
46 H u0 p0 c0 {10,S}
47 H u0 p0 c0 {10,S}
48 H u0 p0 c0 {11,S}
49 H u0 p0 c0 {11,S}
50 H u0 p0 c0 {11,S}
51 H u0 p0 c0 {10,S}
52 H u0 p0 c0 {24,S}
53 H u0 p0 c0 {25,S}
54 H u0 p0 c0 {26,S}
55 H u0 p0 c0 {27,S}
56 H u0 p0 c0 {30,S}
57 H u0 p0 c0 {28,S}
58 H u0 p0 c0 {29,S}
59 H u0 p0 c0 {1,S}
60 H u0 p0 c0 {2,S}
61 H u0 p0 c0 {3,S}
""",
    "trimer_rad44": """\
trimer_rad44(28)
multiplicity 2
1  O u0 p2 c0 {21,S} {59,S}
2  O u0 p2 c0 {22,S} {60,S}
3  O u0 p2 c0 {23,S} {61,S}
4  C u0 p0 c0 {5,S} {13,S} {31,S} {32,S}
5  C u0 p0 c0 {4,S} {12,S} {33,S} {34,S}
6  C u0 p0 c0 {15,S} {24,S} {35,S} {36,S}
7  C u0 p0 c0 {14,S} {37,S} {38,S} {39,S}
8  C u0 p0 c0 {16,S} {40,S} {41,S} {42,S}
9  C u0 p0 c0 {17,S} {43,S} {44,S} {45,S}
10 C u0 p0 c0 {19,S} {46,S} {47,S} {51,S}
11 C u0 p0 c0 {18,S} {48,S} {49,S} {50,S}
12 C u0 p0 c0 {5,S} {17,B} {22,B}
13 C u0 p0 c0 {4,S} {21,B} {26,B}
14 C u0 p0 c0 {7,S} {16,B} {21,B}
15 C u0 p0 c0 {6,S} {22,B} {28,B}
16 C u0 p0 c0 {8,S} {14,B} {25,B}
17 C u0 p0 c0 {9,S} {12,B} {27,B}
18 C u0 p0 c0 {11,S} {20,B} {29,B}
19 C u0 p0 c0 {10,S} {23,B} {30,B}
20 C u0 p0 c0 {18,B} {23,B} {24,S}
21 C u0 p0 c0 {1,S} {13,B} {14,B}
22 C u0 p0 c0 {2,S} {12,B} {15,B}
23 C u0 p0 c0 {3,S} {19,B} {20,B}
24 C u1 p0 c0 {6,S} {20,S} {52,S}
25 C u0 p0 c0 {16,B} {26,B} {53,S}
26 C u0 p0 c0 {13,B} {25,B} {54,S}
27 C u0 p0 c0 {17,B} {28,B} {55,S}
28 C u0 p0 c0 {15,B} {27,B} {56,S}
29 C u0 p0 c0 {18,B} {30,B} {57,S}
30 C u0 p0 c0 {19,B} {29,B} {58,S}
31 H u0 p0 c0 {4,S}
32 H u0 p0 c0 {4,S}
33 H u0 p0 c0 {5,S}
34 H u0 p0 c0 {5,S}
35 H u0 p0 c0 {6,S}
36 H u0 p0 c0 {6,S}
37 H u0 p0 c0 {7,S}
38 H u0 p0 c0 {7,S}
39 H u0 p0 c0 {7,S}
40 H u0 p0 c0 {8,S}
41 H u0 p0 c0 {8,S}
42 H u0 p0 c0 {8,S}
43 H u0 p0 c0 {9,S}
44 H u0 p0 c0 {9,S}
45 H u0 p0 c0 {9,S}
46 H u0 p0 c0 {10,S}
47 H u0 p0 c0 {10,S}
48 H u0 p0 c0 {11,S}
49 H u0 p0 c0 {11,S}
50 H u0 p0 c0 {11,S}
51 H u0 p0 c0 {10,S}
52 H u0 p0 c0 {24,S}
53 H u0 p0 c0 {25,S}
54 H u0 p0 c0 {26,S}
55 H u0 p0 c0 {27,S}
56 H u0 p0 c0 {28,S}
57 H u0 p0 c0 {29,S}
58 H u0 p0 c0 {30,S}
59 H u0 p0 c0 {1,S}
60 H u0 p0 c0 {2,S}
61 H u0 p0 c0 {3,S}
""",
}


# ---------------------------------------------------------------------------
# round-55 adversarial-review closures (P1-1..P1-4, M1-M4, P2a)
# ---------------------------------------------------------------------------

def _r55_lines(caplog, prefix):
    return [r.getMessage() for r in caplog.records
            if r.getMessage().startswith(prefix)]


class TestR55RegistryEntryIsolation:
    """round-55 P1-1: entries returned by register_candidate()/
    lookup_candidate() are FULLY detached defensive copies -- no nested
    mutable container (censuses/epochs sets, bucket_by_census dict --
    previously returned LIVE -- or the per-census sighting sets) reaches
    internal registry state."""

    def test_returned_entries_fully_detached(self, clean_registry):
        from rmgpy.polymer_conduit import (CENSUS_REGISTRY,
                                           lookup_candidate,
                                           register_candidate)

        returned = register_candidate("KD(1)<>KP(2)", "r93_general",
                                      "DEFERRED_A", epoch=7)
        register_candidate("KD(1)<>KP(2)", "conduit_echo", "CONDUIT_ECHO",
                           epoch=8)
        # Independent pre-mutation snapshot (deepcopy detaches even if the
        # lookup were returning live internals).
        expected = copy.deepcopy(lookup_candidate("KD(1)<>KP(2)"))
        victim = lookup_candidate("KD(1)<>KP(2)")

        # Mutate EVERY mutable field, top-level and nested, on BOTH the
        # register()-returned and the lookup()-returned views.
        for entry in (returned, victim):
            entry["censuses"].add("evil_census")
            entry["epochs"].add(999)
            entry["bucket_by_census"]["evil_census"] = "EVIL"
            entry["bucket_by_census"]["r93_general"] = "EVIL"
            for sightings in entry["bucket_sightings_by_census"].values():
                sightings.add("EVIL")
            entry["bucket_sightings_by_census"]["evil_census"] = {"EVIL"}
            entry["effective_bucket"] = "EVIL"
            entry["shadow_bucket"] = "EVIL"
            entry["precedence"] = "EVIL"

        refetched = lookup_candidate("KD(1)<>KP(2)")
        assert refetched == expected
        assert "evil_census" not in refetched["censuses"]
        assert "EVIL" not in refetched["bucket_by_census"].values()
        assert all("EVIL" not in s for s in
                   refetched["bucket_sightings_by_census"].values())
        assert refetched["epochs"] == {7, 8}
        # Derived registry aggregates are untouched too.
        eff, overlap, total, resight = CENSUS_REGISTRY.counts()
        assert total == 1
        assert eff["DEFERRED_A"] == 1
        assert "EVIL" not in eff


class TestR55HealthLineLifecycle:
    """round-55 P1-2/P1-4: the ORACLE HEALTH/1 line fires from the
    PRODUCTION lifecycle boundary (reset_conduit_state, the exact hook
    rmgpy/rmg/main.py RMG.initialize calls; the process-exit atexit hook
    closes the final lifecycle), exactly once per epoch-reset cycle, even
    for a zero-row lifecycle -- and carries the corrected
    ``validation=lazy-first-sight`` token."""

    def test_exactly_one_health_line_per_production_lifecycle(
            self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import (annotate_feature_radical,
                                           census_summary,
                                           reset_conduit_state)

        with caplog.at_level(_logging.WARNING):
            # Production run boundary (RMG.initialize): opens a lifecycle
            # (closing whatever preceded it -- drop those lines).
            reset_conduit_state()
            caplog.clear()
            # Production-shaped census traffic: the FR warn-once hook
            # (rmgpy.polymer._warn_once_refused) drives exactly this call.
            annotate_feature_radical("FRP(1) => Q(2)", reason="qssa-invalid")
            # P1-3: summaries are READ-ONLY -- no health line, twice over.
            census_summary()
            census_summary()
            assert not _r55_lines(caplog,
                                  "CONDUIT CLASSIFIER ORACLE HEALTH/1")
            # Next epoch boundary closes the lifecycle: one health line.
            reset_conduit_state()
            health = _r55_lines(caplog, "CONDUIT CLASSIFIER ORACLE HEALTH/1")
            assert len(health) == 1
            assert "validation=lazy-first-sight" in health[0]
            # A second summary inside the NEW lifecycle adds nothing.
            census_summary()
            assert len(_r55_lines(
                caplog, "CONDUIT CLASSIFIER ORACLE HEALTH/1")) == 1

    def test_zero_row_lifecycle_still_emits_exactly_one_health_line(
            self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import reset_conduit_state

        with caplog.at_level(_logging.WARNING):
            reset_conduit_state()   # open
            caplog.clear()
            reset_conduit_state()   # close, zero census rows in between
            health = _r55_lines(caplog, "CONDUIT CLASSIFIER ORACLE HEALTH/1")
            assert health == [
                "CONDUIT CLASSIFIER ORACLE HEALTH/1 configured=3 active=0 "
                "dropped=3 source=deck validation=lazy-first-sight"]
            # All three never-sighted declared labels resolved at the
            # boundary as missing-species.
            label_lines = _r55_lines(
                caplog, "CONDUIT CLASSIFIER LABEL VALIDATION/1")
            assert len(label_lines) == 3
            assert all("reason=missing-species" in l for l in label_lines)


class TestR55ReadOnlySummaryAndLateSightings:
    """round-55 P1-3: census_summary() no longer finalizes (an early
    diagnostic call must not poison later legitimate sightings), and a
    sighting after TRUE finalization of a dropped label surfaces the loud
    versioned ORACLE ANOMALY/1 line while keeping the finalized verdict."""

    def test_early_summary_does_not_poison_later_sighting(
            self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import (active_pool_state_labels,
                                           census_summary,
                                           record_from_reaction)
        from rmgpy.reaction import Reaction

        # EARLY diagnostic summary, before any sighting: with the old
        # mutating finalizer this permanently resolved every declared
        # label as missing-species.
        census_summary()
        assert active_pool_state_labels() == frozenset()

        pf = TestLabelOracleAcceptanceBattery._pf()
        proxy = TestLabelOracleAcceptanceBattery._iso_proxy_species(
            pf, "trimer_rad33")
        with caplog.at_level(_logging.WARNING):
            record = record_from_reaction(
                Reaction(reactants=[proxy], products=[pf], reversible=True),
                [pf])

        # The later legitimate first sighting still validates correctly.
        assert record["label_isomorphism_divergence"] is False
        assert "trimer_rad33" in active_pool_state_labels()
        label_lines = _r55_lines(
            caplog, "CONDUIT CLASSIFIER LABEL VALIDATION/1")
        assert any("label=trimer_rad33" in l and "status=active" in l
                   and "reason=iso-pass" in l for l in label_lines)
        assert not any("label=trimer_rad33" in l
                       and "reason=missing-species" in l
                       for l in label_lines)
        assert not _r55_lines(caplog, "CONDUIT CLASSIFIER ORACLE ANOMALY/1")

    def test_post_finalization_sighting_emits_versioned_anomaly(
            self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import (active_pool_state_labels,
                                           finalize_label_oracle,
                                           record_from_reaction)
        from rmgpy.reaction import Reaction

        # TRUE end-of-lifecycle finalization with no sightings: all three
        # declared labels drop as missing-species.
        finalize_label_oracle()

        pf = TestLabelOracleAcceptanceBattery._pf()
        with caplog.at_level(_logging.WARNING):
            for _ in range(2):  # repeat sighting: anomaly is once-per-label
                proxy = TestLabelOracleAcceptanceBattery._iso_proxy_species(
                    pf, "trimer_rad33")
                record_from_reaction(
                    Reaction(reactants=[proxy], products=[pf],
                             reversible=True),
                    [pf])

        anomalies = _r55_lines(caplog, "CONDUIT CLASSIFIER ORACLE ANOMALY/1")
        assert len(anomalies) == 1
        assert "label=trimer_rad33" in anomalies[0]
        assert "event=post-finalization-sighting" in anomalies[0]
        assert "reason=missing-species" in anomalies[0]
        assert "action=keep-finalized-verdict" in anomalies[0]
        # The finalized dropped verdict is KEPT -- no post-hoc
        # re-validation contradicting the already-emitted health line.
        assert "trimer_rad33" not in active_pool_state_labels()


class TestR55UnclassifiedTotalDrift:
    """round-55 M1: unclassified_total is recomputed from final registry
    state, pinned against BOTH drift directions -- a key drifting INTO
    UNCLASSIFIED on a same-census resight (a stale first-sighting
    incremental counter would under-count) and a key drifting OUT of it
    via feature-radical precedence (an incremental counter could never
    decrement)."""

    def test_recompute_tracks_drift_in_both_directions(self,
                                                       clean_registry):
        from rmgpy.polymer_conduit import (CENSUS_REGISTRY, census_summary,
                                           lookup_candidate,
                                           register_candidate)

        # Direction 1: classified -> UNCLASSIFIED (same-census resight;
        # UNCLASSIFIED is the most conservative bucket, so it wins the
        # sighting-set resolution).
        register_candidate("KA(1)<>KB(2)", "r93_general", "DEFERRED_A")
        assert CENSUS_REGISTRY.unclassified_total == 0
        register_candidate("KA(1)<>KB(2)", "r93_general", "UNCLASSIFIED")
        assert lookup_candidate(
            "KA(1)<>KB(2)")["effective_bucket"] == "UNCLASSIFIED"
        assert CENSUS_REGISTRY.unclassified_total == 1  # stale counter: 0

        # Direction 2: UNCLASSIFIED -> classified (FR census lands on the
        # key; feature_radical precedence takes effective_bucket over).
        register_candidate("KC(3)<>KD(4)", "r93_general", "UNCLASSIFIED")
        assert CENSUS_REGISTRY.unclassified_total == 2
        register_candidate("KC(3)<>KD(4)", "feature_radical",
                           "FEATURE_RADICAL")
        assert lookup_candidate(
            "KC(3)<>KD(4)")["effective_bucket"] == "FEATURE_RADICAL"
        assert CENSUS_REGISTRY.unclassified_total == 1  # stale counter: 2

        # The loud census reports the recomputed value.
        assert "UNCLASSIFIED=1" in census_summary()


class TestR55TieBreakDeterminism:
    """round-55 M2: _most_conservative_bucket() is a pure function of the
    sighting SET. Known buckets have pairwise-distinct ranks (no ties by
    construction); UNKNOWN bucket strings (reachable through the public
    register_candidate API) share the defensive after-everything rank and
    used to resolve by set-iteration order, which str-hash randomization
    varies across processes -- rank ties now break lexicographically."""

    def test_known_buckets_latest_declaration_wins_any_order(self):
        import itertools

        from rmgpy.polymer_conduit import _most_conservative_bucket

        for perm in itertools.permutations(
                ("ADMISSIBLE_A", "DEFERRED_C", "UNCLASSIFIED")):
            assert _most_conservative_bucket(list(perm)) == "UNCLASSIFIED"
            assert _most_conservative_bucket(set(perm)) == "UNCLASSIFIED"
        for perm in itertools.permutations(
                ("ADMISSIBLE_B", "DEFERRED_F")):
            assert _most_conservative_bucket(list(perm)) == "DEFERRED_F"

    def test_unknown_bucket_ties_break_deterministically(self):
        import itertools

        from rmgpy.polymer_conduit import (BUCKET_DECLARATION_ORDER,
                                           _most_conservative_bucket)

        # Two unknown names share the defensive rank: the lexicographic
        # max wins, in EVERY input ordering and container shape.
        for perm in itertools.permutations(
                ("ZZ_NOT_A_BUCKET", "AA_NOT_A_BUCKET", "MM_NOT_A_BUCKET")):
            assert _most_conservative_bucket(list(perm)) == "ZZ_NOT_A_BUCKET"
            assert _most_conservative_bucket(set(perm)) == "ZZ_NOT_A_BUCKET"
        # The defensive rank still never loses to any KNOWN bucket.
        for known in BUCKET_DECLARATION_ORDER:
            assert _most_conservative_bucket(
                {known, "AA_NOT_A_BUCKET"}) == "AA_NOT_A_BUCKET"

    def test_registry_end_state_order_independent_with_unknown_buckets(
            self, clean_registry):
        # round-56 F4 (flipped from the round-55 pin, itself never pushed):
        # register_candidate now ENFORCES the closed BUCKET_DECLARATION_ORDER
        # vocabulary at registration -- an unknown bucket string is coerced
        # to UNCLASSIFIED (so it can never dominate _most_conservative_bucket
        # by sorting after every known bucket). Determinism is preserved:
        # registering the unknown names in either order leaves an IDENTICAL
        # end state, now resolving to UNCLASSIFIED rather than to the (typo'd)
        # unknown string.
        from rmgpy.polymer_conduit import (lookup_candidate,
                                           register_candidate,
                                           reset_conduit_state)

        def end_state(order):
            reset_conduit_state()
            for bucket in order:
                register_candidate("KT(1)<>KU(2)", "r93_general", bucket)
            return lookup_candidate("KT(1)<>KU(2)")

        forward = end_state(["B_WEIRD_BUCKET", "A_WEIRD_BUCKET"])
        backward = end_state(["A_WEIRD_BUCKET", "B_WEIRD_BUCKET"])
        assert forward == backward
        assert forward["effective_bucket"] == "UNCLASSIFIED"


class TestR56UnknownBucketCoercion:
    """round-56 F4: register_candidate validates bucket names against the
    ENFORCED closed BUCKET_DECLARATION_ORDER vocabulary. An unknown name is
    coerced to UNCLASSIFIED with a loud versioned anomaly line, so a typo'd
    bucket can never reach _most_conservative_bucket and dominate the
    conservative classification."""

    KEY = "COERCE_chain(1)<>COERCE_pool(2)+COERCE_gas(3)"

    def test_unknown_bucket_coerced_and_anomaly_emitted(
            self, caplog, clean_registry):
        import logging as _logging
        from rmgpy.polymer_conduit import (lookup_candidate,
                                           register_candidate)
        with caplog.at_level(_logging.WARNING):
            entry = register_candidate(self.KEY, "r93_general",
                                       "DEFRRED_A")  # typo of DEFERRED_A
        anomalies = _r55_lines(caplog, "CONDUIT CENSUS BUCKET ANOMALY/1")
        assert len(anomalies) == 1
        assert "census=r93_general" in anomalies[0]
        assert "bucket=DEFRRED_A" in anomalies[0]
        assert "action=coerced-unclassified" in anomalies[0]
        # the census entry resolves exactly as if UNCLASSIFIED was registered
        assert entry["effective_bucket"] == "UNCLASSIFIED"
        assert entry["bucket_by_census"]["r93_general"] == "UNCLASSIFIED"
        assert lookup_candidate(self.KEY)["effective_bucket"] == "UNCLASSIFIED"

    def test_unknown_bucket_matches_explicit_unclassified_registration(
            self, clean_registry):
        # An unknown bucket registered under a key must leave the registry
        # in the SAME state as registering UNCLASSIFIED directly under an
        # otherwise-identical key -- and be deterministic across repeats.
        from rmgpy.polymer_conduit import (lookup_candidate,
                                           register_candidate,
                                           reset_conduit_state)

        def run(bucket):
            reset_conduit_state()
            register_candidate("RK(1)<>RP(2)", "r93_general", bucket)
            return lookup_candidate("RK(1)<>RP(2)")

        coerced_a = run("TOTALLY_BOGUS_BUCKET")
        coerced_b = run("TOTALLY_BOGUS_BUCKET")  # determinism: repeat run
        explicit = run("UNCLASSIFIED")
        assert coerced_a == coerced_b            # deterministic
        assert coerced_a["effective_bucket"] == explicit["effective_bucket"]
        assert coerced_a["bucket_by_census"] == explicit["bucket_by_census"]

    def test_known_buckets_never_emit_the_coercion_anomaly(
            self, caplog, clean_registry):
        # sanity: every declared bucket registers WITHOUT a coercion line.
        import logging as _logging
        from rmgpy.polymer_conduit import (BUCKET_DECLARATION_ORDER,
                                           register_candidate)
        with caplog.at_level(_logging.WARNING):
            for i, bucket in enumerate(BUCKET_DECLARATION_ORDER):
                register_candidate(f"OK{i}(1)<>OK{i}(2)", "r93_general",
                                   bucket)
        assert not _r55_lines(caplog, "CONDUIT CENSUS BUCKET ANOMALY/1")


class TestR55Conduit4ArtifactReplay:
    """round-55 M3: end-to-end replay of the conduit4 label-oracle finding
    with the REAL structures -- the three declared trimer_rad* species
    (verbatim adjacency lists from poly_102_conduit4's
    species_dictionary.txt: C27H31O3 chain-scale RADICALS, not pool
    species) sighted against the run's own phenol_formaldehyde pool
    (deck monomer SMILES). Every declared label validates and DROPS as
    iso-fail; health line configured=3 active=0 dropped=3 exact; ZERO
    DIVERGENCE/1 lines and ZERO deny=classifier-divergence tokens (the
    validated label test and the isomorphism test agree on every
    sighting). Row shape is a minimal synthetic refused row; the species
    and pool structures are the artifact's own."""

    def _real_pool(self):
        from rmgpy.polymer import Polymer
        return Polymer(monomer=CONDUIT4_POOL_MONOMER_SMILES,
                       **CONDUIT4_POOL_KWARGS)

    def _real_trimer(self, label):
        import rmgpy.molecule  # M3 precondition: importable in rmg_env
        from rmgpy.molecule import Molecule
        from rmgpy.species import Species
        adj = CONDUIT4_TRIMER_ADJLISTS[label]
        body = "\n".join(adj.splitlines()[1:])  # drop the name(idx) header
        sp = Species(molecule=[Molecule().from_adjacency_list(body)])
        sp.label = label
        sp.index = CONDUIT4_TRIMER_INDICES[label]
        return sp

    def test_conduit4_replay_real_structures(self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import (CHAIN_SCALE_MW,
                                           POOL_STATE_RESOLVABLE_LABELS,
                                           annotate_refused_row,
                                           evaluate_conduit_admission,
                                           finalize_label_oracle)
        from rmgpy.reaction import Reaction

        pf = self._real_pool()
        assert getattr(pf, "monomer_mw_g_mol", 0.0) == pytest.approx(
            134.17, abs=0.05)

        suffixes = []
        with caplog.at_level(_logging.WARNING):
            for label in POOL_STATE_RESOLVABLE_LABELS:
                sp = self._real_trimer(label)
                mol = sp.molecule[0]
                # Real-structure pins: the artifact species are C27H31O3
                # trimer radicals, genuinely chain-scale.
                assert mol.get_formula() == "C27H31O3"
                assert mol.get_molecular_weight() * 1000.0 == pytest.approx(
                    403.53, abs=0.05)
                assert mol.get_molecular_weight() * 1000.0 >= CHAIN_SCALE_MW

                rxn = Reaction(reactants=[sp], products=[pf],
                               reversible=True)
                verdict = evaluate_conduit_admission(rxn, [pf])
                assert verdict.admitted is False
                assert verdict.deny_reason != "classifier-divergence"
                suffix = annotate_refused_row(rxn, [pf], verdict=verdict)
                assert suffix  # annotation genuinely produced
                assert "deny=classifier-divergence" not in suffix
                suffixes.append(suffix)
            finalize_label_oracle()

        # configured=3 active=0 dropped=3, exact token match.
        health = _r55_lines(caplog, "CONDUIT CLASSIFIER ORACLE HEALTH/1")
        assert health == [
            "CONDUIT CLASSIFIER ORACLE HEALTH/1 configured=3 active=0 "
            "dropped=3 source=deck validation=lazy-first-sight"]
        # Every declared label was sighted with a REAL species and dropped
        # as iso-fail -- never missing-species.
        label_lines = _r55_lines(
            caplog, "CONDUIT CLASSIFIER LABEL VALIDATION/1")
        assert len(label_lines) == 3
        assert all("status=dropped" in l and "reason=iso-fail" in l
                   for l in label_lines)
        # Zero divergence, in the log and in every emitted suffix.
        assert not _r55_lines(caplog, "CONDUIT CLASSIFIER DIVERGENCE/1")
        assert not any("deny=classifier-divergence" in s for s in suffixes)


class TestR55FrCountPreservation:
    """round-55 M4: the FR census count (genuine qssa-* upstream blockers)
    is bookkept independently of the label-oracle/divergence feature --
    divergence-dropped labels and echo sightings never leak into (or out
    of) the FEATURE_RADICAL effective-bucket population."""

    def test_fr_count_unaffected_by_label_oracle_divergence(
            self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import (CENSUS_REGISTRY,
                                           annotate_feature_radical,
                                           annotate_refused_row,
                                           candidate_key_from_label,
                                           evaluate_conduit_admission,
                                           lookup_candidate)
        from rmgpy.reaction import Reaction

        # Genuine qssa-invalid/unassessable FR rows (the real FR census).
        annotate_feature_radical("FRA(1) + HX(2) => YA(3)",
                                 reason="qssa-invalid")
        annotate_feature_radical("FRB(4) => YB(5)",
                                 reason="qssa-unassessable")
        eff0, _overlap0, total0, _r0 = CENSUS_REGISTRY.counts()
        assert eff0["FEATURE_RADICAL"] == 2

        pf = TestLabelOracleAcceptanceBattery._pf()
        with caplog.at_level(_logging.WARNING):
            # (i) a divergence-dropped sighting: undeclared label on a
            # genuinely isomorphic chip -> DIVERGENCE/1 +
            # deny=classifier-divergence on its row.
            chip = TestLabelOracleAcceptanceBattery._iso_proxy_species(
                pf, "undeclared_chip")
            rxn_div = Reaction(reactants=[chip], products=[pf],
                               reversible=True)
            v = evaluate_conduit_admission(rxn_div, [pf])
            assert v.deny_reason == "classifier-divergence"
            annotate_refused_row(rxn_div, [pf], verdict=v)
            # (ii) a declared label dropping as iso-fail.
            bogus = TestLabelOracleAcceptanceBattery._noniso_species(
                "CCO", "trimer_rad33")
            rxn_drop = Reaction(reactants=[bogus], products=[pf],
                                reversible=True)
            annotate_refused_row(rxn_drop, [pf],
                                 verdict=evaluate_conduit_admission(
                                     rxn_drop, [pf]))
            # (iii) a conduit-deferred ECHO through the FR warn-once hook:
            # registers as conduit_echo, never as feature_radical.
            annotate_feature_radical("FRC(6) => YC(7)",
                                     reason="conduit-deferred")

        assert _r55_lines(caplog, "CONDUIT CLASSIFIER DIVERGENCE/1")

        eff1, _overlap1, total1, _r1 = CENSUS_REGISTRY.counts()
        # The FR census count is EXACTLY what the qssa rows put there.
        assert eff1["FEATURE_RADICAL"] == eff0["FEATURE_RADICAL"] == 2
        assert eff1["CONDUIT_ECHO"] == 1
        assert total1 == total0 + 3
        # The genuine FR keys' ledger entries are untouched by (i)-(iii).
        for fr_label in ("FRA(1) + HX(2) => YA(3)", "FRB(4) => YB(5)"):
            entry = lookup_candidate(candidate_key_from_label(fr_label))
            assert entry["censuses"] == {"feature_radical"}
            assert entry["effective_bucket"] == "FEATURE_RADICAL"
        # And neither divergence row acquired FR membership.
        for rxn in (rxn_div, rxn_drop):
            from rmgpy.polymer_conduit import conduit_candidate_key, \
                record_from_reaction
            key = conduit_candidate_key(
                record_from_reaction(rxn, [pf]))
            assert "feature_radical" not in lookup_candidate(key)["censuses"]


# ---------------------------------------------------------------------------
# round-56 additions (F1 lifecycle close, F2 race, optional artifact guard)
# ---------------------------------------------------------------------------

class TestR56ProductionFinishClose:
    """round-56 F1: the label-oracle lifecycle is closed at the PRODUCTION
    end-of-run point (close_conduit_lifecycle, wired at rmgpy/rmg/main.py
    RMG.finish), so each run's ORACLE HEALTH/1 line lands at its OWN finish
    -- not deferred to the next run's initialize, and not left to the
    fragile atexit path. close_conduit_lifecycle() is idempotent per
    lifecycle, so a double close (or a following initialize) never
    re-emits."""

    def test_single_run_finish_emits_exactly_one_health_line(
            self, caplog, clean_registry):
        import logging as _logging
        from rmgpy.polymer_conduit import (close_conduit_lifecycle,
                                           reset_conduit_state)
        with caplog.at_level(_logging.WARNING):
            reset_conduit_state()          # RMG.initialize: opens lifecycle
            caplog.clear()
            close_conduit_lifecycle()      # RMG.finish: closes + emits once
            health = _r55_lines(caplog, "CONDUIT CLASSIFIER ORACLE HEALTH/1")
            assert len(health) == 1
            assert "validation=lazy-first-sight" in health[0]

    def test_back_to_back_runs_each_emit_at_own_finish(
            self, caplog, clean_registry):
        import logging as _logging
        from rmgpy.polymer_conduit import (close_conduit_lifecycle,
                                           reset_conduit_state)
        with caplog.at_level(_logging.WARNING):
            # Run 1: initialize then finish.
            reset_conduit_state()
            caplog.clear()
            close_conduit_lifecycle()
            assert len(_r55_lines(
                caplog, "CONDUIT CLASSIFIER ORACLE HEALTH/1")) == 1
            caplog.clear()
            # Run 2 initialize must NOT re-emit run 1's (already-closed) line.
            reset_conduit_state()
            assert not _r55_lines(
                caplog, "CONDUIT CLASSIFIER ORACLE HEALTH/1")
            # Run 2 finish emits ITS OWN line.
            close_conduit_lifecycle()
            assert len(_r55_lines(
                caplog, "CONDUIT CLASSIFIER ORACLE HEALTH/1")) == 1

    def test_double_close_emits_exactly_one_health_line(
            self, caplog, clean_registry):
        import logging as _logging
        from rmgpy.polymer_conduit import (close_conduit_lifecycle,
                                           reset_conduit_state)
        with caplog.at_level(_logging.WARNING):
            reset_conduit_state()
            caplog.clear()
            close_conduit_lifecycle()
            close_conduit_lifecycle()      # second close: no-op
            assert len(_r55_lines(
                caplog, "CONDUIT CLASSIFIER ORACLE HEALTH/1")) == 1


class TestR56NoteSightingFinalizeRace:
    """round-56 F2: note_sighting() validates the declared label OUTSIDE the
    oracle lock; if finalize() interleaves in that gap and marks the label
    dropped/missing-species, the sighting must NOT silently return the
    dropped verdict -- it routes through the loud post-finalization
    ANOMALY/1 path. Pinned DETERMINISTICALLY (no threads, no sleeps): the
    validation callable is monkeypatched to drive finalize() itself
    mid-call, exactly simulating the interleave."""

    def test_finalize_interleave_during_validation_fires_anomaly(
            self, caplog, clean_registry, monkeypatch):
        import logging as _logging
        import rmgpy.polymer_conduit as pc

        label = pc.POOL_STATE_RESOLVABLE_LABELS[0]
        oracle = pc._LABEL_ORACLE
        oracle.reset()  # fresh, open, un-finalized lifecycle

        def racing_validate(species, row_pools):
            # finalize() wins the race in the gap between the cached-check
            # and the setdefault: it marks every never-sighted declared
            # label missing-species, THEN this validator returns normally
            # with a would-be ACTIVE verdict.
            oracle.finalize()
            return ("active", pc.LABEL_VALIDATION_ISO_PASS, "somepool")

        monkeypatch.setattr(pc, "_validate_declared_label", racing_validate)
        with caplog.at_level(_logging.WARNING):
            result = oracle.note_sighting(label, object(), [])

        # The finalized dropped verdict is KEPT -- not the active verdict
        # the validator returned.
        assert result[0] == "dropped"
        assert result[1] == pc.LABEL_VALIDATION_MISSING_SPECIES
        # The loud anomaly fired instead of a silent drop.
        anomalies = _r55_lines(caplog, "CONDUIT CLASSIFIER ORACLE ANOMALY/1")
        assert len(anomalies) == 1
        assert f"label={label}" in anomalies[0]
        assert "event=post-finalization-sighting" in anomalies[0]
        assert "action=keep-finalized-verdict" in anomalies[0]
        # The label was NOT promoted to active; the finalized verdict stands.
        assert label not in oracle.active_labels()


class TestR56Conduit4ArtifactDriftGuard:
    """round-56 (optional): guard the embedded CONDUIT4_TRIMER_ADJLISTS
    fixtures against drift from the on-disk conduit4 artifact. Hermetic:
    skips when the READ-ONLY artifact is absent (CI / other machines)."""

    ARTIFACT = ("/home/alon/runs/RMG/poly_102_conduit4/chemkin/"
                "species_dictionary.txt")

    def test_embedded_trimer_blocks_match_on_disk_artifact(self):
        import os
        if not os.path.isfile(self.ARTIFACT):
            pytest.skip("conduit4 artifact not present on this machine")
        with open(self.ARTIFACT) as f:
            lines = f.read().splitlines()
        for label, idx in CONDUIT4_TRIMER_INDICES.items():
            header = "{0}({1})".format(label, idx)
            assert header in lines, (
                "{0} block missing from {1}".format(header, self.ARTIFACT))
            start = lines.index(header)
            block = [lines[start]]
            for ln in lines[start + 1:]:
                if ln.strip() == "":
                    break
                block.append(ln)
            on_disk = "\n".join(block) + "\n"
            assert on_disk == CONDUIT4_TRIMER_ADJLISTS[label], (
                "{0} block in {1} has drifted from the embedded fixture"
                .format(label, self.ARTIFACT))


# ---------------------------------------------------------------------------
# Round-58: NO-GO review closure (P1-A execute() lifecycle-close
# determinism, P1-B same-label concurrent first-sighting nondeterminism,
# P2-A close-failure exc_info, P2-B anomaly rate-limiting, P2-C deny-reason
# emission validation).
# ---------------------------------------------------------------------------

class TestR58ExecuteLifecycleCloseDeterminism:
    """P1-A: RMG.execute() previously called self.finish() (which closes
    the conduit lifecycle) only on its SUCCESS tail; the two early-return
    termination paths (walltime exhaustion, max-iterations) and any
    exception skip finish() entirely, so a caller that catches an
    exception from execute() defeats the atexit guard too.

    Full end-to-end drives of RMG.execute() require the real generation
    loop/database and are out of scope here; the exception path (the most
    safety-critical -- "exceptions also skip the close") is pinned
    end-to-end against a minimal RMG stub, idempotency of the close helper
    is pinned directly, and the early-return shape is pinned structurally:
    both early-return statements are lexically inside the try body whose
    finally unconditionally calls the guarded close helper, so Python's
    own try/finally semantics guarantee exactly one close on those paths
    too."""

    def test_exception_propagates_and_close_runs_exactly_once(
            self, monkeypatch):
        from rmgpy.rmg.main import RMG

        rmg = RMG.__new__(RMG)
        closes = []
        monkeypatch.setattr(RMG, "_conduit_lifecycle_close",
                            lambda self: closes.append(1))

        def boom(self, requires_rms=False):
            raise RuntimeError("boom-p1a")

        monkeypatch.setattr(RMG, "register_listeners", boom)
        with pytest.raises(RuntimeError, match="boom-p1a"):
            rmg.execute(initialize=False)
        # the finally ran exactly once; the exception was not swallowed.
        assert closes == [1]

    def test_close_helper_is_idempotent_under_a_double_call(
            self, monkeypatch):
        from rmgpy.rmg.main import RMG

        calls = []

        def fake_close_conduit_lifecycle():
            calls.append(1)

        monkeypatch.setattr("rmgpy.polymer_conduit.close_conduit_lifecycle",
                            fake_close_conduit_lifecycle)
        rmg = RMG.__new__(RMG)
        rmg._conduit_lifecycle_close()
        rmg._conduit_lifecycle_close()
        # RMG._conduit_lifecycle_close is a thin, unconditional guarded
        # wrapper -- it is close_conduit_lifecycle() itself
        # (round-56 F1(b)) that is idempotent per lifecycle; this pins that
        # the wrapper calls through on every invocation without adding its
        # own (redundant) guard state.
        assert calls == [1, 1]

    def test_close_failure_logs_with_exc_info(self, caplog, monkeypatch):
        """P2-A: the guarded close's except clause must log with
        exc_info so a real close defect is diagnosable, not just a bare
        one-line warning."""
        import logging as _logging

        from rmgpy.rmg.main import RMG

        def boom_close():
            raise RuntimeError("close-boom-p2a")

        monkeypatch.setattr("rmgpy.polymer_conduit.close_conduit_lifecycle",
                            boom_close)
        rmg = RMG.__new__(RMG)
        with caplog.at_level(_logging.WARNING):
            rmg._conduit_lifecycle_close()  # must not raise
        records = [r for r in caplog.records
                  if "conduit lifecycle close failed" in r.getMessage()]
        assert len(records) == 1
        assert records[0].exc_info is not None
        assert records[0].exc_info[1].args[0] == "close-boom-p2a"

    def test_close_wrapper_try_finally_encloses_initialize_and_all_returns(
            self):
        """P1-2 (round-60): robust structural pin, replacing a prior weak
        pin that grabbed whichever ``ast.Try`` node ``ast.walk()`` visited
        FIRST and never verified it was the INTENDED close wrapper -- a
        pin that would have stayed green even if it had locked onto the
        wrong try block, and that went stale the moment P1-1 moved
        ``self.initialize(**kwargs)`` inside the try (a pin keyed on
        position, not identity, cannot see that kind of restructuring).

        This version instead: (1) finds the try/finally block(s) inside
        RMG.execute() by MATCHING a call to the lifecycle-close helper
        appearing in the finally body (never by position), (2) asserts
        EXACTLY ONE such wrapper exists, and (3) asserts BOTH the
        ``self.initialize(...)`` call AND every ``return`` statement in
        execute() are lexically inside that try's BODY (by AST descendant
        identity -- never inside its finally or orelse, never outside the
        try). That directly pins the P1-1 shape: initialize() opens a
        conduit lifecycle near its own start, so it must be covered by the
        same try/finally that closes it, same as every other exit path."""
        import ast
        import inspect
        import textwrap

        from rmgpy.rmg.main import RMG

        source = textwrap.dedent(inspect.getsource(RMG.execute))
        func = ast.parse(source).body[0]
        assert isinstance(func, ast.FunctionDef)

        def call_name(node):
            """Best-effort dotted-tail name of a Call's func, e.g. both
            'self._conduit_lifecycle_close()' and a bare
            '_conduit_lifecycle_close()' resolve to
            '_conduit_lifecycle_close'."""
            target = node.func
            if isinstance(target, ast.Name):
                return target.id
            if isinstance(target, ast.Attribute):
                return target.attr
            return None

        def contains_close_call(stmts):
            """round-65 P3 tightening (Finding 5): the qualifying call must
            be the FIRST top-level statement in ``stmts`` (the finally
            body's own statement list) -- ``stmts and stmts[0]`` is an
            ``ast.Expr`` whose ``.value`` is an ``ast.Call`` to the
            lifecycle-close helper. round-64's version accepted the call
            ANYWHERE among top-level finally statements, which would still
            pass even if an earlier top-level statement in the same finally
            body could terminate control flow first (e.g. a ``raise``
            before the close call would mean close() sometimes never
            runs). The real production finally body
            (rmgpy/rmg/main.py execute()) is confirmed to be exactly
            ``self._conduit_lifecycle_close()`` as its sole/first
            statement, so this simpler, stricter, purely positional rule
            matches production as-is -- no main.py reordering needed
            (option 2 of the Finding-5 decision tree). Deliberately does
            NOT descend into nested compound statements (If, While, With,
            Try, FunctionDef, AsyncFunctionDef, Lambda): a call buried in a
            dead ``if False:`` branch or inside a nested def lexically
            inside finally -- code that would never actually run on the
            real exit path -- must not satisfy the pin either."""
            if not stmts:
                return False
            stmt = stmts[0]
            if isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Call):
                name = call_name(stmt.value)
                if name and name.endswith("_conduit_lifecycle_close"):
                    return True
            return False

        close_wrappers = [n for n in ast.walk(func)
                          if isinstance(n, ast.Try)
                          and contains_close_call(n.finalbody)]
        assert len(close_wrappers) == 1, (
            "expected exactly one try/finally in execute() whose finally "
            "calls the conduit lifecycle-close helper, found %d"
            % len(close_wrappers))
        try_node = close_wrappers[0]

        # Every node lexically inside the try's BODY (not finally/orelse),
        # identified by AST object identity so membership below is exact.
        try_body_nodes = [n for stmt in try_node.body
                          for n in ast.walk(stmt)]

        initialize_calls = [
            n for n in ast.walk(func)
            if isinstance(n, ast.Call) and call_name(n) == "initialize"]
        assert initialize_calls, "no call to initialize(...) found in execute()"
        for call in initialize_calls:
            assert call in try_body_nodes, (
                "self.initialize(...) must be inside the try body whose "
                "finally closes the conduit lifecycle (P1-1)")

        all_returns = [n for n in ast.walk(func) if isinstance(n, ast.Return)]
        assert len(all_returns) >= 2  # at least the two early-return paths
        for ret in all_returns:
            assert ret in try_body_nodes, (
                "every return in execute() must be inside the try body "
                "whose finally closes the conduit lifecycle")
        # ... and never inside the finally itself (would mask exceptions).
        finally_returns = [n for stmt in try_node.finalbody
                           for n in ast.walk(stmt)
                           if isinstance(n, ast.Return)]
        assert not finally_returns

    def test_matcher_rejects_close_call_nested_in_dead_branch(self):
        """round-64 P3 matcher self-test: pins the tightened
        ``contains_close_call`` helper above against the exact adversarial
        shape the finding calls out -- a ``_conduit_lifecycle_close(...)``
        call that is lexically present inside the finally body's AST
        subtree, but only inside a dead ``if False:`` branch that would
        never actually execute. The (deliberately re-derived, matching the
        production test's round-65 positional matcher) must REJECT this
        snippet: ``finalbody[0]`` is an ``If``, not the qualifying call."""
        import ast
        import textwrap

        source = textwrap.dedent("""
            def f():
                try:
                    pass
                finally:
                    if False:
                        _conduit_lifecycle_close(x)
            """)
        func = ast.parse(source).body[0]
        try_node = [n for n in ast.walk(func) if isinstance(n, ast.Try)][0]

        def call_name(node):
            target = node.func
            if isinstance(target, ast.Name):
                return target.id
            if isinstance(target, ast.Attribute):
                return target.attr
            return None

        def contains_close_call(stmts):
            # round-65 positional matcher (Finding 5): only finalbody[0]
            # qualifies -- see the production test's contains_close_call
            # docstring above for the full rationale.
            if not stmts:
                return False
            stmt = stmts[0]
            if isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Call):
                name = call_name(stmt.value)
                if name and name.endswith("_conduit_lifecycle_close"):
                    return True
            return False

        assert not contains_close_call(try_node.finalbody), (
            "matcher must reject a lifecycle-close call buried in a dead "
            "if-False branch nested inside finally -- it is not a "
            "top-level finalbody statement")

    def test_matcher_rejects_close_call_preceded_by_raise(self):
        """round-65 P3 matcher self-test (Finding 5): a ``finally:`` body
        where a top-level ``raise`` precedes the lifecycle-close call must
        be REJECTED -- the raise can terminate control flow before close()
        ever runs, so 'the close call appears somewhere in finalbody' is
        not a binding pin. The tightened positional matcher (finalbody[0]
        must itself be the qualifying call) rejects this by construction,
        since the raise sits at finalbody[0] and the close call is only
        finalbody[1]."""
        import ast
        import textwrap

        source = textwrap.dedent("""
            def f():
                try:
                    pass
                finally:
                    raise RuntimeError("boom")
                    _conduit_lifecycle_close()
            """)
        func = ast.parse(source).body[0]
        try_node = [n for n in ast.walk(func) if isinstance(n, ast.Try)][0]

        def call_name(node):
            target = node.func
            if isinstance(target, ast.Name):
                return target.id
            if isinstance(target, ast.Attribute):
                return target.attr
            return None

        def contains_close_call(stmts):
            if not stmts:
                return False
            stmt = stmts[0]
            if isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Call):
                name = call_name(stmt.value)
                if name and name.endswith("_conduit_lifecycle_close"):
                    return True
            return False

        assert not contains_close_call(try_node.finalbody), (
            "matcher must reject a finally body where a raise precedes "
            "the lifecycle-close call -- close() could be skipped")

    def test_initialize_exception_after_reset_conduit_state_closes_lifecycle(
            self, monkeypatch):
        """P1-1 direct regression pin (round-60): initialize() opens a
        fresh conduit lifecycle via reset_conduit_state() near its OWN
        start, well before load_input/load_database/walltime-parsing/
        constraints-setup run. Before P1-1, execute()'s try/finally began
        only AFTER self.initialize(**kwargs) returned, so an exception
        raised by any LATER initialize() step -- after
        reset_conduit_state() had already opened a lifecycle -- left that
        lifecycle unclosed until the next run's reset or the last-resort
        atexit guard.

        This monkeypatches self.initialize with a stand-in that
        reproduces exactly that shape: it calls the REAL
        reset_conduit_state() (genuinely opening a lifecycle) and then
        raises, mimicking a failure in one of initialize()'s later steps.
        The exception must still propagate (never swallowed), and the
        lifecycle-close helper must have run exactly once."""
        from rmgpy.polymer_conduit import reset_conduit_state
        from rmgpy.rmg.main import RMG

        rmg = RMG.__new__(RMG)
        closes = []
        monkeypatch.setattr(RMG, "_conduit_lifecycle_close",
                            lambda self: closes.append(1))

        def fake_initialize(self, **kwargs):
            # Reproduces the real initialize()'s shape: open a lifecycle
            # first (this IS the production reset_conduit_state() call),
            # then fail as if a later initialize() step raised.
            reset_conduit_state()
            raise RuntimeError("boom-after-reset-p1-1")

        monkeypatch.setattr(RMG, "initialize", fake_initialize)
        with pytest.raises(RuntimeError, match="boom-after-reset-p1-1"):
            rmg.execute(initialize=True)
        # the finally ran exactly once; the exception was not swallowed.
        assert closes == [1]


class TestR58ConcurrentFirstSightingDisagreement:
    """P1-B: two concurrent FIRST sightings of the SAME declared label both
    validate outside the oracle lock; whichever setdefault() wins silently
    defines the lifecycle verdict. Previously, if the two validations
    DISAGREED, the losing call's result simply vanished with no anomaly
    logged. Pinned deterministically (no threads): the validation callable
    is monkeypatched to plant a DIFFERING "winning" verdict directly into
    the oracle's ledger (simulating the other concurrent first-sighting
    winning the race) before returning this call's own (different)
    computed verdict."""

    def test_disagreeing_concurrent_first_sighting_logs_anomaly_and_keeps_first_stored(
            self, caplog, monkeypatch, clean_registry):
        import logging as _logging

        import rmgpy.polymer_conduit as pc

        label = pc.POOL_STATE_RESOLVABLE_LABELS[0]
        oracle = pc._LABEL_ORACLE
        oracle.reset()  # fresh, open, un-finalized lifecycle

        winning_verdict = ("active", pc.LABEL_VALIDATION_ISO_PASS,
                           "winner_pool")

        def racing_validate(species, row_pools):
            # Simulate a concurrent first sighting of the SAME label
            # winning the setdefault race, in the gap between this call's
            # own validation and its setdefault, with a DIFFERENT verdict.
            with oracle._lock:
                oracle._validated.setdefault(label, winning_verdict)
            return ("dropped", pc.LABEL_VALIDATION_ISO_FAIL, None)

        monkeypatch.setattr(pc, "_validate_declared_label", racing_validate)
        with caplog.at_level(_logging.WARNING):
            result = oracle.note_sighting(label, object(), [])

        # first-write-wins: the cached (first-stored) verdict is kept, NOT
        # this (losing) call's freshly-computed one.
        assert result == winning_verdict

        anomalies = [r.getMessage() for r in caplog.records
                    if r.getMessage().startswith(
                        "CONDUIT CLASSIFIER ORACLE ANOMALY/1")]
        assert len(anomalies) == 1
        assert f"label={label}" in anomalies[0]
        assert "reason=concurrent-validation-disagreement" in anomalies[0]
        assert "kept=" in anomalies[0] and "winner_pool" in anomalies[0]
        assert "discarded=" in anomalies[0]
        # the label was genuinely promoted active via the KEPT verdict.
        assert label in oracle.active_labels()

    def test_agreeing_concurrent_first_sighting_stays_silent(
            self, caplog, monkeypatch, clean_registry):
        """Sanity: when the raced-in verdict and this call's computed
        verdict AGREE, no disagreement anomaly fires (only the normal
        LABEL VALIDATION/1 line, from whichever call actually won)."""
        import logging as _logging

        import rmgpy.polymer_conduit as pc

        label = pc.POOL_STATE_RESOLVABLE_LABELS[1]
        oracle = pc._LABEL_ORACLE
        oracle.reset()

        agreed_verdict = ("dropped", pc.LABEL_VALIDATION_ISO_FAIL, None)

        def racing_validate(species, row_pools):
            with oracle._lock:
                oracle._validated.setdefault(label, agreed_verdict)
            return agreed_verdict

        monkeypatch.setattr(pc, "_validate_declared_label", racing_validate)
        with caplog.at_level(_logging.WARNING):
            result = oracle.note_sighting(label, object(), [])

        assert result == agreed_verdict
        anomalies = [r.getMessage() for r in caplog.records
                    if "concurrent-validation-disagreement" in r.getMessage()]
        assert not anomalies


class TestR58BucketAnomalyRateLimiting:
    """P2-B: the unknown-bucket coercion anomaly in register_candidate's
    coercion path previously logged once PER REGISTRATION -- a hot
    mistyped-bucket producer spams the log. The LOG LINE is now
    rate-limited to once per (census, bucket-name) per lifecycle; the
    underlying COERCION must still apply on every registration."""

    def test_repeated_same_typo_logs_once_but_coerces_every_call(
            self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import lookup_candidate, register_candidate

        keys = [f"RATE{i}_chain(1)<>RATE{i}_pool(2)" for i in range(5)]
        with caplog.at_level(_logging.WARNING):
            for key in keys:
                register_candidate(key, "r93_general", "TYPO_BUCKET_RL")

        anomalies = _r55_lines(caplog, "CONDUIT CENSUS BUCKET ANOMALY/1")
        assert len(anomalies) == 1
        # the coercion side-effect happened for EVERY registration, even
        # though only the first emitted a log line.
        for key in keys:
            assert lookup_candidate(key)["effective_bucket"] == "UNCLASSIFIED"

    def test_lifecycle_reset_rearms_the_line(self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import register_candidate, reset_conduit_state

        key = "RATE_RESET_chain(1)<>RATE_RESET_pool(2)"
        with caplog.at_level(_logging.WARNING):
            register_candidate(key, "r93_general", "TYPO_BUCKET_RESET")
            reset_conduit_state()
            register_candidate(key, "r93_general", "TYPO_BUCKET_RESET")

        anomalies = _r55_lines(caplog, "CONDUIT CENSUS BUCKET ANOMALY/1")
        assert len(anomalies) == 2

    def test_different_census_same_bucket_name_logs_separately(
            self, caplog, clean_registry):
        # the dedup key is (census, bucket) -- a typo repeated under a
        # DIFFERENT census must still get its own line.
        import logging as _logging

        from rmgpy.polymer_conduit import register_candidate

        with caplog.at_level(_logging.WARNING):
            register_candidate("SPLIT1(1)<>SPLIT1(2)", "r93_general",
                               "TYPO_BUCKET_SPLIT")
            register_candidate("SPLIT2(1)<>SPLIT2(2)", "feature_radical",
                               "TYPO_BUCKET_SPLIT")

        anomalies = _r55_lines(caplog, "CONDUIT CENSUS BUCKET ANOMALY/1")
        assert len(anomalies) == 2


class TestR58AdmissionTokenAnomalyRateLimiting:
    """P2-B: the F3 token-anomaly lines in admission_census_suffix
    (unregistered-key and mixed-census-membership echo-token misuse)
    previously logged once PER CALL; a hot producer of the same misuse
    spams the log. Now rate-limited to once per (candidate_key, reason)
    per lifecycle -- the fail-closed substitution still applies on every
    call."""

    def test_unregistered_key_misuse_logs_once_but_denies_every_call(
            self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import _deny, admission_census_suffix

        key = "UNREG_RL(1)<>UNREG_RL(2)"
        with caplog.at_level(_logging.WARNING):
            for _ in range(4):
                suffix = admission_census_suffix(_deny(key, "echo-not-evaluated"))
                assert "deny=admission-token-anomaly stage=final" in suffix

        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 1

    def test_lifecycle_reset_rearms_the_line(self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import (_deny, admission_census_suffix,
                                           reset_conduit_state)

        key = "UNREG_RL_RESET(1)<>UNREG_RL_RESET(2)"
        with caplog.at_level(_logging.WARNING):
            admission_census_suffix(_deny(key, "echo-not-evaluated"))
            reset_conduit_state()
            admission_census_suffix(_deny(key, "echo-not-evaluated"))

        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 2

    def test_mixed_membership_misuse_logs_once_but_denies_every_call(
            self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import (_deny, admission_census_suffix,
                                           register_candidate)

        key = "MIXED_RL_chain(1)<>MIXED_RL_pool(2)+MIXED_RL_gas(3)"
        register_candidate(key, "r93_general", "ADMISSIBLE_A")
        register_candidate(key, "conduit_echo", "CONDUIT_ECHO")
        with caplog.at_level(_logging.WARNING):
            for _ in range(3):
                suffix = admission_census_suffix(_deny(key, "echo-not-evaluated"))
                assert "deny=admission-token-anomaly stage=final" in suffix

        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 1


class TestR58DenyReasonEmissionValidation:
    """P2-C: admission_census_suffix() previously emitted
    verdict.deny_reason without validating membership in the closed
    ADMISSION_DENY_REASONS vocabulary, so a future typo'd reason would
    leak as a structured token into the census output. Now validated at
    emission time: a known reason passes through untouched; an unknown one
    is substituted with the SAME conservative in-vocabulary fail-closed
    token ('admission-token-anomaly') this function already uses for the
    echo-token misuses, with a loud anomaly line."""

    def test_known_deny_reason_passes_through_untouched(self, clean_registry):
        from rmgpy.polymer_conduit import _deny, admission_census_suffix

        suffix = admission_census_suffix(
            _deny("KNOWN(1)<>KNOWN(2)", "gas-mw-over-threshold"))
        assert "deny=gas-mw-over-threshold stage=final" in suffix

    def test_g7_dynamically_suffixed_reason_passes_through_untouched(
            self, clean_registry):
        # the admission-evaluation-error:* G7 family is a documented
        # exception to exact-membership: it is suffixed with the caught
        # exception's type name and must never be flagged as unknown.
        from rmgpy.polymer_conduit import _deny, admission_census_suffix

        suffix = admission_census_suffix(
            _deny("G7KEY(1)<>G7KEY(2)",
                 "admission-evaluation-error:RuntimeError"))
        assert ("deny=admission-evaluation-error:RuntimeError stage=final"
                in suffix)

    def test_unknown_deny_reason_substituted_with_conservative_fallback(
            self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import _deny, admission_census_suffix

        with caplog.at_level(_logging.WARNING):
            suffix = admission_census_suffix(
                _deny("TYPO_KEY(1)<>TYPO_KEY(2)", "gas-mw-over-threshhold"))

        assert "deny=admission-token-anomaly stage=final" in suffix
        assert "gas-mw-over-threshhold" not in suffix

        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 1
        assert "key=TYPO_KEY(1)<>TYPO_KEY(2)" in anomalies[0]
        assert "reason=unknown-deny-reason" in anomalies[0]
        assert "value=gas-mw-over-threshhold" in anomalies[0]
        assert "action=fail-closed-deny" in anomalies[0]

    def test_every_declared_reason_is_accepted_without_anomaly(
            self, caplog, clean_registry):
        # sanity sweep: no member of the closed vocabulary ever trips the
        # unknown-deny-reason anomaly.
        import logging as _logging

        from rmgpy.polymer_conduit import (ADMISSION_DENY_REASONS, _deny,
                                           admission_census_suffix,
                                           register_candidate)

        # give the echo-not-evaluated member the registered echo sighting
        # it requires so this sweep exercises every OTHER reason cleanly.
        register_candidate("SWEEP_KEY", "conduit_echo", "CONDUIT_ECHO")
        with caplog.at_level(_logging.WARNING):
            for reason in sorted(ADMISSION_DENY_REASONS):
                admission_census_suffix(_deny("SWEEP_KEY", reason))
        assert not _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")


# ---------------------------------------------------------------------------
# round-60 NO-GO closures (P1-1/P1-2, P2-1..P2-3, P3-1)
# ---------------------------------------------------------------------------

class TestP21AdmissionEvaluationErrorTokenSanitization:
    """P2-1: the dynamic ``admission-evaluation-error:*`` deny-reason family
    (G7's fail-closed catch-all, suffixed with the failing exception's TYPE
    NAME) bypassed the closed-vocabulary check in
    admission_census_suffix() by design -- but that also meant it bypassed
    ANY charset validation, so a crafted exception type name could inject
    embedded whitespace/'='/brackets that forge extra structured key=value
    tokens into the census line. The chokepoint now sanitizes that dynamic
    suffix to a strict [A-Za-z0-9_.-]+ charset before serialization."""

    def test_forged_reason_is_sanitized_and_cannot_inject_extra_tokens(
            self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import _deny, admission_census_suffix

        forged = "admission-evaluation-error: x deny=would_admit=1"
        with caplog.at_level(_logging.WARNING):
            suffix = admission_census_suffix(
                _deny("INJECT_KEY(1)<>INJECT_KEY(2)", forged))

        # The injected sequence must never survive as a parsable extra
        # 'deny=' or 'would_admit=' token: a naive parser splitting the
        # suffix on whitespace and '=' must not find either as a
        # standalone key.
        tokens = suffix.split()
        assert not any(t.startswith("deny=would_admit") for t in tokens)
        assert not any(t == "deny=" for t in tokens)
        assert "would_admit=1" not in suffix.replace(
            "would_admit=1 deny=None", "")  # not the legitimate admit token
        # No embedded space or '=' survives inside the sanitized suffix
        # itself (only the surrounding, legitimate token structure may
        # contain them).
        assert " x " not in suffix
        assert "deny=would_admit=1" not in suffix

        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 1
        assert "reason=unsanitized-error-token" in anomalies[0]
        assert "key=INJECT_KEY(1)<>INJECT_KEY(2)" in anomalies[0]
        # (iii) round-66: charset-dirty-but-short raw keeps the original
        # round-65 `action=sanitized` label unchanged.
        assert "action=sanitized" in anomalies[0]

    def test_clean_g7_reason_still_passes_through_untouched(
            self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import _deny, admission_census_suffix

        with caplog.at_level(_logging.WARNING):
            suffix = admission_census_suffix(
                _deny("CLEAN_KEY(1)<>CLEAN_KEY(2)",
                      "admission-evaluation-error:RuntimeError"))
        assert ("deny=admission-evaluation-error:RuntimeError stage=final"
                in suffix)
        assert not _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")

    def test_sanitization_anomaly_is_rate_limited(
            self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import _deny, admission_census_suffix

        key = "INJECT_RL_KEY(1)<>INJECT_RL_KEY(2)"
        with caplog.at_level(_logging.WARNING):
            for _ in range(4):
                admission_census_suffix(
                    _deny(key, "admission-evaluation-error: bad=token here"))
        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 1


# ---------------------------------------------------------------------------
# round-64 P2 closure: sanitized deny reasons must not silently collapse
# distinct raw causes
# ---------------------------------------------------------------------------

class TestP64SanitizedErrorTokenAuditTrail:
    """P2 (round-64): the sanitized-suffix TOKEN emitted into the census
    (``admission-evaluation-error:<sanitized>``) deliberately stays a
    closed-charset value -- that vocabulary must not change -- but two
    DISTINCT raw pre-sanitization strings can collide onto the SAME
    sanitized token (e.g. "bad token" and "bad=token" both become
    "bad_token"), and prior to this fix the raw value was never recorded
    anywhere auditable, and the anomaly line's rate-limit key was a fixed
    literal that would suppress the SECOND distinct raw's anomaly line for
    the same candidate. The anomaly line now carries a ``raw_sha=`` short
    digest of the raw string, keyed into the warn-once rate limiter, so
    each distinct raw gets its own audit-trail line while exact repeats of
    the same raw are still deduped."""

    def test_colliding_raws_share_token_but_get_distinct_audited_anomalies(
            self, caplog, clean_registry):
        import hashlib
        import logging as _logging
        import re

        from rmgpy.polymer_conduit import _deny, admission_census_suffix

        raw_a = "bad token"
        raw_b = "bad=token"
        # Pin the premise: both raws really do collide under the module's
        # own sanitizer charset before relying on that collision below.
        sanitize = lambda s: re.sub(r"[^A-Za-z0-9_.-]", "_", s) or "unsanitizable"
        assert sanitize(raw_a) == sanitize(raw_b) == "bad_token"

        key = "COLLIDE_KEY(1)<>COLLIDE_KEY(2)"
        with caplog.at_level(_logging.WARNING):
            suffix_a = admission_census_suffix(
                _deny(key, f"admission-evaluation-error:{raw_a}"))
            suffix_b = admission_census_suffix(
                _deny(key, f"admission-evaluation-error:{raw_b}"))
            # Repeating the exact same raw a second time must not add a
            # third anomaly line -- the rate limiter still suppresses
            # exact repeats of the same raw.
            admission_census_suffix(
                _deny(key, f"admission-evaluation-error:{raw_a}"))

        # (i) both raws produce the SAME sanitized census deny token --
        # the token vocabulary/cardinality must not change.
        assert "deny=admission-evaluation-error:bad_token" in suffix_a
        assert "deny=admission-evaluation-error:bad_token" in suffix_b

        # (ii) processing both raws produces TWO distinct anomaly lines
        # with DIFFERENT raw_sha values.
        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 2

        sha_a = hashlib.sha256(raw_a.encode("utf-8")).hexdigest()[:12]
        sha_b = hashlib.sha256(raw_b.encode("utf-8")).hexdigest()[:12]
        assert sha_a != sha_b
        assert any(f"raw_sha={sha_a}" in line for line in anomalies)
        assert any(f"raw_sha={sha_b}" in line for line in anomalies)


# ---------------------------------------------------------------------------
# round-66 P2 closure (Finding 1): the sanitized suffix itself was still
# unbounded -- cap it at 64 chars, apply the SAME cap to the census token
# and the `sanitized=` log field, and surface the truncation in `action=`
# even when the raw was already charset-clean.
# ---------------------------------------------------------------------------

class TestR66SanitizedSuffixLengthCap:
    """round-66 Finding 1: admission_census_suffix() capped raw_b64 to the
    first 96 encoded bytes (round-64/65), but never capped the SANITIZED
    suffix that is actually emitted as the census deny token and as the
    anomaly line's `sanitized=` field -- a huge raw_suffix still produced
    an equally huge sanitized string in both places, so the log/census-line
    DoS the earlier round claimed to close was never actually closed for
    the sanitized value itself. Now: the sanitized form is capped at 64
    chars, the identical capped string is used for BOTH the census token
    and `sanitized=`, and the anomaly line's `action=` field records
    whether the raw was charset-dirty (`sanitized`), merely over-length
    (`truncated`), or both (`sanitized-truncated`)."""

    def test_long_clean_raw_is_truncated_to_64_chars_with_truncated_action(
            self, caplog, clean_registry):
        """(i) a >64-char charset-CLEAN raw: the census token is exactly
        the first 64 chars, the anomaly line fires (this is new -- prior
        to round-66 a clean-but-long raw sailed through with NO anomaly at
        all), `action=truncated` (over-length only, not charset-dirty),
        and `raw_len=` reports the FULL original length, not the capped
        one."""
        import logging as _logging
        import re

        from rmgpy.polymer_conduit import _deny, admission_census_suffix

        raw = "y" * 100  # charset-clean (matches [A-Za-z0-9_.-]+ already)
        key = "LONGCLEAN_KEY(1)<>LONGCLEAN_KEY(2)"
        with caplog.at_level(_logging.WARNING):
            suffix = admission_census_suffix(
                _deny(key, f"admission-evaluation-error:{raw}"))

        expected_token = "y" * 64
        assert re.fullmatch(r"[A-Za-z0-9_.-]+", expected_token)
        assert f"deny=admission-evaluation-error:{expected_token}" in suffix
        # the un-truncated 100-char token must NOT appear anywhere.
        assert f"admission-evaluation-error:{raw}" not in suffix

        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 1
        assert "action=truncated" in anomalies[0]
        assert "raw_len=100" in anomalies[0]
        assert f"sanitized={expected_token}" in anomalies[0]
        assert f"key={key}" in anomalies[0]

    def test_dirty_and_long_raw_gets_single_combined_action(
            self, caplog, clean_registry):
        """(ii) a raw that is BOTH charset-dirty and over-length (after
        sanitization) produces exactly one anomaly line with the combined
        `action=sanitized-truncated`, and the census token is the 64-char
        prefix of the SANITIZED (not raw) string."""
        import logging as _logging
        import re as _re

        from rmgpy.polymer_conduit import _deny, admission_census_suffix

        raw = ("bad token " * 10).strip()  # dirty (spaces) and >64 chars
        sanitized_full = _re.sub(r"[^A-Za-z0-9_.-]", "_", raw)
        assert sanitized_full != raw  # premise: really is charset-dirty
        assert len(sanitized_full) > 64  # premise: really is over-length
        expected_token = sanitized_full[:64]

        key = "DIRTYLONG_KEY(1)<>DIRTYLONG_KEY(2)"
        with caplog.at_level(_logging.WARNING):
            suffix = admission_census_suffix(
                _deny(key, f"admission-evaluation-error:{raw}"))

        assert f"deny=admission-evaluation-error:{expected_token}" in suffix
        assert _re.fullmatch(r"[A-Za-z0-9_.-]+", expected_token)

        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 1
        assert "action=sanitized-truncated" in anomalies[0]
        assert f"raw_len={len(raw)}" in anomalies[0]
        assert f"sanitized={expected_token}" in anomalies[0]

    def test_short_clean_raw_unchanged_no_anomaly_no_truncation(
            self, caplog, clean_registry):
        """(iv) regression guard: a <=64-char charset-clean raw still
        produces NO anomaly line at all and passes through untouched,
        exactly as before round-66."""
        import logging as _logging

        from rmgpy.polymer_conduit import _deny, admission_census_suffix

        raw = "RuntimeError"  # short and clean
        assert len(raw) <= 64
        key = "SHORTCLEAN_KEY(1)<>SHORTCLEAN_KEY(2)"
        with caplog.at_level(_logging.WARNING):
            suffix = admission_census_suffix(
                _deny(key, f"admission-evaluation-error:{raw}"))

        assert f"deny=admission-evaluation-error:{raw}" in suffix
        assert not _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")


# ---------------------------------------------------------------------------
# round-66 P3 closure (Finding 3): raw_b64 switched from padded standard
# base64 to unpadded base64url, so the field can never contain a literal
# '=' that could confuse a naive split-on-first-'=' key=value parser.
# ---------------------------------------------------------------------------

class TestR66RawB64IsUnpaddedBase64Url:
    """round-66 Finding 3: raw_b64 previously used standard
    `base64.b64encode`, which can emit trailing `=` padding characters.
    The log line's key=value contract splits on the FIRST `=`, so this was
    never actually a parse-breaking bug for THIS field (it's the last
    field on the line) -- but it is an unnecessary ambiguity given the
    line otherwise promises whitespace/'='-clean values. Switched to
    unpadded base64url (`base64.urlsafe_b64encode(...).rstrip(b"=")`)."""

    def test_raw_b64_field_matches_unpadded_base64url_recipe(
            self, caplog, clean_registry):
        import base64
        import hashlib
        import logging as _logging

        from rmgpy.polymer_conduit import _deny, admission_census_suffix

        # A raw crafted so its standard-base64 encoding would need padding
        # ('=') -- exercises the exact case this fix removes ambiguity for.
        # ("bad token" is 9 bytes -- a multiple of 3 -- and pads with
        # nothing; one extra char makes the byte count 10, which needs 2
        # '=' padding chars in standard base64.)
        raw = "bad tokens"
        key = "B64_KEY(1)<>B64_KEY(2)"
        with caplog.at_level(_logging.WARNING):
            admission_census_suffix(
                _deny(key, f"admission-evaluation-error:{raw}"))

        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 1

        raw_bytes = raw.encode("utf-8", errors="backslashreplace")
        expected_b64 = base64.urlsafe_b64encode(
            raw_bytes[:96]).rstrip(b"=").decode("ascii")
        # Premise check: the padded standard-base64 form of this raw really
        # does contain '=', so this test actually exercises the fix.
        assert base64.b64encode(raw_bytes[:96]).decode("ascii").endswith("=")
        assert "=" not in expected_b64
        assert f"raw_b64={expected_b64}" in anomalies[0]

        # Also confirm the reset_conduit_state / rate-limit dedup key
        # (raw_sha) is unaffected by the base64 encoding change.
        raw_sha = hashlib.sha256(raw_bytes).hexdigest()[:12]
        assert f"raw_sha={raw_sha}" in anomalies[0]


class TestP22ConcurrentValidationDisagreementRateLimiting:
    """P2-2: the concurrent-validation-disagreement anomaly in
    _LabelOracleState.note_sighting() previously fired on every LOSING
    thread of a race, producing log spam. It is now deduped once per
    (label, lifecycle) via the same warn-once-set + lock pattern used
    elsewhere in this module/class."""

    def test_repeated_disagreement_on_same_label_logs_once(
            self, caplog, monkeypatch, clean_registry):
        """Drives the REAL note_sighting() code path (not a reimplemented
        stand-in): the validation callable is monkeypatched exactly as in
        TestR58ConcurrentFirstSightingDisagreement to plant a differing
        "winning" verdict and force a losing race. Between iterations the
        cached verdict is popped (WITHOUT touching the per-lifecycle
        _disagreement_warned set) to simulate several losing "first
        sighting" threads landing across the same lifecycle -- exactly the
        log-spam shape P2-2 fixes -- and asserts only the first logs."""
        import logging as _logging

        import rmgpy.polymer_conduit as pc

        label = pc.POOL_STATE_RESOLVABLE_LABELS[2]
        oracle = pc._LABEL_ORACLE
        oracle.reset()

        winning_verdict = ("active", pc.LABEL_VALIDATION_ISO_PASS,
                           "winner_pool")

        def racing_validate(species, row_pools):
            with oracle._lock:
                oracle._validated.setdefault(label, winning_verdict)
            return ("dropped", pc.LABEL_VALIDATION_ISO_FAIL, None)

        monkeypatch.setattr(pc, "_validate_declared_label", racing_validate)
        with caplog.at_level(_logging.WARNING):
            for _ in range(3):
                with oracle._lock:
                    oracle._validated.pop(label, None)
                result = oracle.note_sighting(label, object(), [])
                assert result == winning_verdict

        anomalies = [r.getMessage() for r in caplog.records
                    if "concurrent-validation-disagreement" in r.getMessage()]
        assert len(anomalies) == 1

    def test_reset_rearms_the_disagreement_line(self, clean_registry):
        from rmgpy.polymer_conduit import _LABEL_ORACLE, reset_conduit_state

        label = "trimer_rad38"
        _LABEL_ORACLE._disagreement_warned.add(label)
        assert label in _LABEL_ORACLE._disagreement_warned
        reset_conduit_state()
        assert label not in _LABEL_ORACLE._disagreement_warned


class TestP31UnknownDenyReasonRateLimiting:
    """P3-1: the unknown-deny-reason TOKEN ANOMALY line in
    admission_census_suffix() previously fired unconditionally on every
    call. It is now rate-limited once per (candidate_key, reason) per
    lifecycle, matching every other TOKEN ANOMALY line in this module."""

    def test_repeated_same_unknown_reason_logs_once_but_denies_every_call(
            self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import _deny, admission_census_suffix

        key = "UNKNOWN_REASON_RL(1)<>UNKNOWN_REASON_RL(2)"
        with caplog.at_level(_logging.WARNING):
            for _ in range(4):
                suffix = admission_census_suffix(
                    _deny(key, "totally-made-up-reason"))
                assert "deny=admission-token-anomaly" in suffix
        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 1

    def test_different_reason_same_key_logs_separately(
            self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import _deny, admission_census_suffix

        key = "UNKNOWN_REASON_SPLIT(1)<>UNKNOWN_REASON_SPLIT(2)"
        with caplog.at_level(_logging.WARNING):
            admission_census_suffix(_deny(key, "made-up-reason-a"))
            admission_census_suffix(_deny(key, "made-up-reason-b"))
        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 2

    def test_lifecycle_reset_rearms_the_line(self, caplog, clean_registry):
        import logging as _logging

        from rmgpy.polymer_conduit import (_deny, admission_census_suffix,
                                           reset_conduit_state)

        key = "UNKNOWN_REASON_RESET(1)<>UNKNOWN_REASON_RESET(2)"
        with caplog.at_level(_logging.WARNING):
            admission_census_suffix(_deny(key, "made-up-reason"))
            reset_conduit_state()
            admission_census_suffix(_deny(key, "made-up-reason"))
        anomalies = _r55_lines(caplog, "CONDUIT ADMISSION TOKEN ANOMALY/1")
        assert len(anomalies) == 2


class TestP23WarnOnceLockDiscipline:
    """P2-3: every check-then-add sequence against the module-level
    warn-once sets (_bucket_anomaly_warned, _admission_token_anomaly_warned)
    -- and the corresponding per-instance set on the label oracle
    (_disagreement_warned) -- is now guarded by a lock, so concurrent
    threads can never both pass the membership check and double-emit.
    Lifecycle reset clears these sets under the SAME lock discipline."""

    def test_warn_once_helper_wins_exactly_once_under_concurrency(self):
        import threading

        from rmgpy.polymer_conduit import _warn_once

        warned = set()
        wins = []
        wins_lock = threading.Lock()

        def attempt():
            if _warn_once(warned, "shared-key", "unused %s", "x"):
                with wins_lock:
                    wins.append(1)

        threads = [threading.Thread(target=attempt) for _ in range(32)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        assert len(wins) == 1

    def test_register_candidate_bucket_anomaly_warns_once_under_barrier_race(
            self, caplog, clean_registry):
        """round-64 P3: the test above pins `_warn_once` directly, launched
        via a plain `Thread.start()` loop with no start barrier -- it can
        pass under fully serial scheduling and never really forces a
        concurrent race, and it never touches a real production caller.
        This test strengthens that pin: it drives the actual PUBLIC
        production entry point (`register_candidate`, the documented
        round-56 F4 chokepoint that funnels an out-of-vocabulary bucket
        through `_warn_once` to emit the ``CONDUIT CENSUS BUCKET
        ANOMALY/1`` line) from N threads held on a `threading.Barrier` so
        they are released SIMULTANEOUSLY -- a real start-line race, not
        launch-order interleaving -- and asserts the anomaly line is
        emitted exactly once despite that."""
        import logging as _logging
        import threading

        from rmgpy.polymer_conduit import register_candidate

        n = 32
        barrier = threading.Barrier(n)
        census = "r93_general"
        bad_bucket = "NOT_A_REAL_BUCKET"

        def attempt(i):
            barrier.wait(timeout=5)
            register_candidate(f"BARRIER_RACE({i})<>BARRIER_RACE_OTHER",
                                census, bad_bucket)

        threads = [threading.Thread(target=attempt, args=(i,))
                  for i in range(n)]
        with caplog.at_level(_logging.WARNING):
            for t in threads:
                t.start()
            for t in threads:
                t.join(timeout=5)
        assert all(not t.is_alive() for t in threads), (
            "a thread failed to join within the deterministic timeout")

        anomalies = _r55_lines(caplog, "CONDUIT CENSUS BUCKET ANOMALY/1")
        assert len(anomalies) == 1

    def test_register_candidate_bucket_anomaly_warns_once_under_forced_interleave(
            self, caplog, clean_registry, monkeypatch):
        """round-65 P3 (Finding 4): the barrier test above starts N=32
        threads on a `threading.Barrier` so they are RELEASED
        simultaneously, but it does not force the interleave BETWEEN the
        `key not in warned_set` membership check and `warned_set.add(key)`
        inside `_warn_once` -- under CPython GIL scheduling, an
        unsynchronized ("if key not in set: add; log") version of that
        code could still, by luck, serialize its check-then-add steps and
        pass the test above. It is not a binding pin on the LOCK itself.

        This test forces that exact interleave deterministically: the
        module-level `_bucket_anomaly_warned` seen-set (a plain
        module-level `set()` -- confirmed by reading
        rmgpy/polymer_conduit.py, `register_candidate` -> `_warn_once`) is
        monkeypatched with a subclass whose `__contains__` -- the
        production membership-check hook `_warn_once` calls via
        `key in warned_set` -- blocks on a `threading.Barrier(N)` (N=4)
        before delegating to the real set behavior.

        Why this is a BINDING pin: the membership snapshot is taken
        (`super().__contains__`) BEFORE the barrier wait, not after. With
        the production `_WARN_ONCE_LOCK` held across the check+add in
        `_warn_once`, threads enter `__contains__` ONE AT A TIME
        (serialized by the lock) -- so the barrier of size N never fills;
        each thread's `barrier.wait(0.1)` times out via
        `BrokenBarrierError` (swallowed), and that thread proceeds alone
        to add+emit, then the lock is released for the next thread. Net
        effect: exactly one emission, same as today. Against a
        mentally-reverted UNSYNCHRONIZED version (no lock around the
        check+add), all N threads instead call `__contains__`
        CONCURRENTLY and EACH snapshots "not present" before any of them
        can reach `warned_set.add(key)` -- so the size-N barrier actually
        FILLS (all N arrive before any times out) with every thread
        already holding an "absent" snapshot, all N then proceed to
        add+emit, and the "exactly one emission" assertion below fails.
        (An earlier draft of this pin waited on the barrier BEFORE
        delegating to `super().__contains__` -- but by the time the
        barrier released, one unsynchronized thread could already have
        added the key, so the remaining threads' real membership checks
        would see "present" and skip emitting, letting the unsynchronized
        code pass the test by scheduling luck. Snapshotting first and
        returning the snapshot after the wait closes that hole.) That is
        the failure this test is designed to catch if the lock discipline
        in `_warn_once` ever regresses."""
        import logging as _logging
        import threading

        import rmgpy.polymer_conduit as pc

        n = 4
        interleave_barrier = threading.Barrier(n)

        class ForcedInterleaveSet(set):
            def __contains__(self, item):
                # Snapshot the real membership result FIRST, then wait --
                # so every thread's decision is locked in before any
                # thread can have added the key (see docstring above).
                result = super().__contains__(item)
                try:
                    interleave_barrier.wait(timeout=0.1)
                except threading.BrokenBarrierError:
                    pass
                return result

        forced_set = ForcedInterleaveSet()
        monkeypatch.setattr(pc, "_bucket_anomaly_warned", forced_set)

        census = "r93_general"
        bad_bucket = "NOT_A_REAL_BUCKET_FORCED"
        start_barrier = threading.Barrier(n)

        def attempt(i):
            start_barrier.wait(timeout=5)
            pc.register_candidate(
                f"FORCED_INTERLEAVE({i})<>FORCED_INTERLEAVE_OTHER",
                census, bad_bucket)

        threads = [threading.Thread(target=attempt, args=(i,))
                  for i in range(n)]
        with caplog.at_level(_logging.WARNING):
            for t in threads:
                t.start()
            for t in threads:
                t.join(timeout=5)
        assert all(not t.is_alive() for t in threads), (
            "a thread failed to join within the deterministic timeout")

        anomalies = _r55_lines(caplog, "CONDUIT CENSUS BUCKET ANOMALY/1")
        assert len(anomalies) == 1

    def test_reset_clears_bucket_and_admission_warn_once_sets_atomically(
            self, clean_registry):
        from rmgpy.polymer_conduit import (_admission_token_anomaly_warned,
                                           _bucket_anomaly_warned,
                                           reset_conduit_state)

        _bucket_anomaly_warned.add(("some_census", "some_bucket"))
        _admission_token_anomaly_warned.add(("some_key", "some_reason"))
        reset_conduit_state()
        assert not _bucket_anomaly_warned
        assert not _admission_token_anomaly_warned
