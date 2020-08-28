from arkane.encorr.reference import ReferenceDatabase
from arkane.encorr.isodesmic import ErrorCancelingScheme, IsodesmicRingScheme
from arkane.modelchem import CompositeLevelOfTheory, LevelOfTheory

database = ReferenceDatabase()
database.load()

lot = CompositeLevelOfTheory(freq=LevelOfTheory(method='wb97mv', basis='def2tzvpd', software='qchem'),
                             energy=LevelOfTheory(method='dlpnoccsd(t)f12', basis='ccpvtzf12',
                                                  auxiliary_basis='augccpvtz/c', cabs='ccpvtzf12cabs', software='orca',
                                                  args=('tightpno',)))

ref_set = database.reference_sets['main']
ref_mapping = {}
for spcs in ref_set:
    try:
        ref_mapping[spcs.index] = spcs.to_error_canceling_spcs(lot)
    except KeyError:
        pass

full_data = {}


def process(tup):
    indx, spc = tup
    error_spcs = [value for key, value in ref_mapping.items() if key != indx]
    scheme = ErrorCancelingScheme(target=spc,
                                  reference_set=error_spcs,
                                  isodesmic_class='rc4',
                                  conserve_ring_size=False,
                                  limit_charges=False,
                                  limit_scope=False,
                                  )

    h298, rxns = scheme.calculate_target_enthalpy(n_reactions_max=100)
    full_data[indx] = (h298, rxns)
    print(f'Done with species {indx}')


spcs = list(ref_mapping.items())
process(spcs[318])
for i in range(20):
    process(spcs[i])
print(full_data)
