title = 'Bond additivity correction fitting for wB97M-V/def2-TZVPD level of theory'

lot = LevelOfTheory(
    method='wB97M-V',
    basis='def2-TZVPD',
    software='Q-Chem'
)

# Petersson-type
bac(
    level_of_theory=lot,
    bac_type='p',  # Petersson
    weighted=True,
    write_to_database=False,
    overwrite=False
)

# Melius-type
bac(
    level_of_theory=lot,
    bac_type='m',  # Melius
    write_to_database=False,
    overwrite=False,
    fit_mol_corr=True,
    global_opt=True,
    global_opt_iter=2,  # Recommended: >=10
)
