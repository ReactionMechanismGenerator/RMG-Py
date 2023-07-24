# Coordinates for C2H5 in Input Orientation (angstroms):
#   C   -0.0001    0.0011   -0.0019
#   C   -0.0001    0.0003    1.4856
#   H    0.9284   -0.0035   -0.5595
#   H   -0.9205   -0.1214   -0.5595
#   H   -0.9148    0.4413    1.8935
#   H    0.0610   -1.0218    1.8949
#   H    0.8552    0.5475    1.8936
conformer(
    label = 'C2H5',
    E0 = (109.427, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(29.0391, 'amu')),
        NonlinearRotor(
            inertia = ([4.8709, 22.2353, 23.9925], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([492.752, 815.246, 982.075, 1065.92, 1194.8, 1404.37, 1469.55, 1485.46, 1487.13, 2955.2, 3046.47, 3090.29, 3150.75, 3251.16], 'cm^-1'),
        ),
        HinderedRotor(
            inertia = (1.11737, 'amu*angstrom^2'),
            symmetry = 6,
            barrier = (0.244029, 'kJ/mol'),
            quantum = None,
            semiclassical = None,
        ),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
)

# Thermodynamics for C2H5:
#   Enthalpy of formation (298 K)   =    29.076 kcal/mol
#   Entropy of formation (298 K)    =    59.128 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      12.353      29.101      59.211      11.338
#            400      14.656      30.451      63.079       5.220
#            500      16.915      32.031      66.594      -1.266
#            600      19.023      33.829      69.867      -8.091
#            800      22.639      38.009      75.855     -22.675
#           1000      25.467      42.832      81.224     -38.392
#           1500      29.998      56.821      92.506     -81.938
#           2000      32.417      72.484     101.502    -130.519
#           2400      33.580      85.700     107.522    -172.353
#    =========== =========== =========== =========== ===========
thermo(
    label = 'C2H5',
    thermo = Wilhoit(
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (153.818, 'J/(mol*K)'),
        a0 = -2.17764,
        a1 = -8.467,
        a2 = 15.934,
        a3 = -7.49338,
        H0 = (
            -91.9983,
            'kJ/mol',
        ),
        S0 = (
            -803.481,
            'J/(mol*K)',
        ),
        B = (1161.93, 'K'),
    ),
)

