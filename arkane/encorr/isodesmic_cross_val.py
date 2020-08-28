import multiprocessing as mp
import pickle
import sys

from arkane.encorr.reference import ReferenceDatabase
from arkane.encorr.isodesmic import ErrorCancelingScheme
from arkane.modelchem import CompositeLevelOfTheory, LevelOfTheory

JOB_ID = int(sys.argv[1])
N_TASKS = 8
CLASSES = ['rc3', 'rc4']
RC_CLASS = CLASSES[JOB_ID // 18]
LABELS = ['A', 'B', 'C', 'D', 'E', 'F']
LOT_LABEL = LABELS[JOB_ID % 6]
CONSTRAINT_LIST = [(False, False), (True, False), (True, True)]
CONTRAINT = CONSTRAINT_LIST[(JOB_ID // 6) % 3]

database = ReferenceDatabase()
database.load()

LOT_DICT = {'A': LevelOfTheory(method='wb97mv', basis='def2tzvpd', software='qchem'),
            'B': LevelOfTheory(method='b3lypd3bj', basis='def2tzvp', software='gaussian'),
            'C': LevelOfTheory(method='wb97xd', basis='def2tzvp', software='gaussian'),
            'D': CompositeLevelOfTheory(freq=LevelOfTheory(method='wb97xd', basis='def2tzvp', software='gaussian'),
                                        energy=LevelOfTheory(method='dlpnoccsd(t)', basis='def2tzvp',
                                                             auxiliary_basis='def2tzvp/c', software='orca',
                                                             args=('normalpno',))),
            'E': CompositeLevelOfTheory(freq=LevelOfTheory(method='wb97mv', basis='def2tzvpd', software='qchem'),
                                        energy=LevelOfTheory(method='dlpnoccsd(t)f12', basis='ccpvdzf12',
                                                             auxiliary_basis='augccpvdz/c', cabs='ccpvdzf12cabs',
                                                             software='orca', args=('tightpno',))),
            'F': CompositeLevelOfTheory(freq=LevelOfTheory(method='wb97mv', basis='def2tzvpd', software='qchem'),
                                        energy=LevelOfTheory(method='dlpnoccsd(t)f12', basis='ccpvtzf12',
                                                             auxiliary_basis='augccpvtz/c', cabs='ccpvtzf12cabs',
                                                             software='orca', args=('tightpno',)))
            }

lot = LOT_DICT[LOT_LABEL]

ref_set = database.reference_sets['main']
ref_mapping = {}
for spcs in ref_set:
    try:
        ref_mapping[spcs.index] = spcs.to_error_canceling_spcs(lot)
    except KeyError:
        pass


def process(tup):
    indx, spc = tup
    error_spcs = [value for key, value in ref_mapping.items() if key != indx]
    scheme = ErrorCancelingScheme(target=spc,
                                  reference_set=error_spcs,
                                  isodesmic_class=RC_CLASS,
                                  conserve_ring_size=False,
                                  limit_charges=CONTRAINT[0],
                                  limit_scope=CONTRAINT[1],
                                  )

    h298, rxns = scheme.calculate_target_enthalpy(n_reactions_max=100)
    lit_h298 = spc.high_level_hf298
    if h298 is not None:
        data = {'Index': indx,
                'CalcHf298': h298.value_si,
                'LitHf298': lit_h298.value_si,
                'Error': (h298.value_si - lit_h298.value_si)/4184.0,
                'Reactions': rxns,
                }
    else:
        data = {'Index': indx,
                'CalcHf298': h298,
                'LitHf298': lit_h298.value_si,
                'Error': None,
                'Reactions': rxns,
                }
    full_data[indx] = data


with mp.Manager() as manager:
    full_data = manager.dict()

    pool = mp.Pool(N_TASKS)
    results = pool.map(process, list(ref_mapping.items())[0])
    pool.close()
    pool.join()

    full_data = dict(full_data)

if any(CONTRAINT):
    constraint_str = ''
    if CONTRAINT[0]:
        constraint_str += 'limit_charges'
    if CONTRAINT[1]:
        constraint_str += '_and_scope'
else:
    constraint_str = 'full'

with open(f'./full_data_pickles/{LOT_LABEL}_{RC_CLASS}_{constraint_str}.pkl', 'wb') as f:
    pickle.dump(full_data, f)
