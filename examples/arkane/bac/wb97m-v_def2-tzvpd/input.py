title = 'Bond additivity correction fitting for wB97M-V/def2-TZVPD level of theory'

# Petersson-type
bac(
    model_chemistry='wb97m-v/def2-tzvpd',
    bac_type='p',  # Petersson
    weighted=True,
    write_to_database=False,
    overwrite=False
)

# Melius-type
bac(
    model_chemistry='wb97m-v/def2-tzvpd',
    bac_type='m',  # Melius
    write_to_database=False,
    overwrite=False,
    fit_mol_corr=True,
    global_opt=True,
    global_opt_iter=2,  # Recommended: >=10
)
