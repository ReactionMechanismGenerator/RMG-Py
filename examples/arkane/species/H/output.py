# Coordinates for H in Input Orientation (angstroms):
#   H    0.0000    0.0000    0.0000
conformer(
    label = 'H',
    E0 = (211.794, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(1.00783, 'amu')),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
)

# Thermodynamics for H:
#   Enthalpy of formation (298 K)   =    52.100 kcal/mol
#   Entropy of formation (298 K)    =    27.387 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300       4.968      52.110      27.420      43.884
#            400       4.968      52.606      28.849      41.067
#            500       4.968      53.103      29.958      38.124
#            600       4.968      53.600      30.863      35.082
#            800       4.968      54.594      32.293      28.760
#           1000       4.968      55.587      33.401      22.186
#           1500       4.968      58.071      35.416       4.948
#           2000       4.968      60.555      36.845     -13.134
#           2400       4.968      62.542      37.751     -28.059
#    =========== =========== =========== =========== ===========
thermo(
    label = 'H',
    thermo = Wilhoit(
        Cp0 = (20.7862, 'J/(mol*K)'),
        CpInf = (20.7862, 'J/(mol*K)'),
        a0 = 0,
        a1 = 0,
        a2 = 0,
        a3 = 0,
        H0 = (211791, 'J/mol'),
        S0 = (
            -3.83529,
            'J/(mol*K)',
        ),
        B = (500, 'K'),
    ),
)

