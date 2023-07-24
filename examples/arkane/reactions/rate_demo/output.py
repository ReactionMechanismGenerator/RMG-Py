# Coordinates for CH4 in Input Orientation (angstroms):
#   C    0.0000    0.0000    0.0000
#   H   -0.6478   -0.7710   -0.4108
#   H   -0.3478    0.9775   -0.3263
#   H   -0.0222   -0.0487    1.0863
#   H    1.0179   -0.1578   -0.3492
conformer(
    label = 'CH4',
    E0 = (
        -82.1866,
        'kJ/mol',
    ),
    modes = [
        IdealGasTranslation(mass=(16.0313, 'amu')),
        NonlinearRotor(
            inertia = ([3.17924, 3.17924, 3.17924], 'amu*angstrom^2'),
            symmetry = 12,
        ),
        HarmonicOscillator(
            frequencies = ([1352.14, 1352.14, 1352.14, 1570.75, 1570.75, 3047.91, 3160.77, 3160.77, 3160.77], 'cm^-1'),
        ),
    ],
    spin_multiplicity = 1,
    optical_isomers = 1,
)

# Coordinates for NH2 in Input Orientation (angstroms):
#   N    0.0002    0.4237    0.0000
#   H   -0.8049   -0.2114    0.0000
#   H    0.8047   -0.2123    0.0000
conformer(
    label = 'NH2',
    E0 = (176.409, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(16.0187, 'amu')),
        NonlinearRotor(
            inertia = ([0.711671, 1.30561, 2.01728], 'amu*angstrom^2'),
            symmetry = 2,
        ),
        HarmonicOscillator(frequencies=([1540.1, 3370.28, 3466.5], 'cm^-1')),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
)

# Coordinates for CH3 in Input Orientation (angstroms):
#   C    0.0000    0.0000    0.0000
#   H    1.0615   -0.1743    0.0539
#   H   -0.6818   -0.8333   -0.0279
#   H   -0.3796    1.0076   -0.0259
conformer(
    label = 'CH3',
    E0 = (136.807, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(15.0235, 'amu')),
        NonlinearRotor(
            inertia = ([1.75359, 1.75359, 3.50717], 'amu*angstrom^2'),
            symmetry = 6,
        ),
        HarmonicOscillator(
            frequencies = ([522.824, 1424.32, 1424.33, 3135.76, 3316.64, 3316.65], 'cm^-1'),
        ),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
)

# Coordinates for NH3 in Input Orientation (angstroms):
#   N    0.0006   -0.0009    0.2796
#   H   -0.4185    0.8433   -0.0894
#   H   -0.5211   -0.7834   -0.0947
#   H    0.9390   -0.0590   -0.0955
conformer(
    label = 'NH3',
    E0 = (
        -53.1186,
        'kJ/mol',
    ),
    modes = [
        IdealGasTranslation(mass=(17.0266, 'amu')),
        NonlinearRotor(
            inertia = ([1.68433, 1.68433, 2.67769], 'amu*angstrom^2'),
            symmetry = 3,
        ),
        HarmonicOscillator(
            frequencies = ([1040.01, 1680.33, 1680.33, 3477.84, 3606.68, 3606.68], 'cm^-1'),
        ),
    ],
    spin_multiplicity = 1,
    optical_isomers = 1,
)

# Coordinates for C2H6 in Input Orientation (angstroms):
#   C    0.7621    0.0250    0.0051
#   C   -0.7621   -0.0250   -0.0051
#   H    1.1915   -0.9227   -0.3213
#   H    1.1371    0.8034   -0.6603
#   H    1.1447    0.2334    1.0048
#   H   -1.1447   -0.2334   -1.0048
#   H   -1.1915    0.9227    0.3213
#   H   -1.1371   -0.8034    0.6603
conformer(
    label = 'C2H6',
    E0 = (
        -93.9171,
        'kJ/mol',
    ),
    modes = [
        IdealGasTranslation(mass=(30.0469, 'amu')),
        NonlinearRotor(
            inertia = ([6.24311, 25.1916, 25.1916], 'amu*angstrom^2'),
            symmetry = 6,
        ),
        HarmonicOscillator(
            frequencies = ([829.963, 829.989, 1009.6, 1230.58, 1230.61, 1419.13, 1432.7, 1515.11, 1515.12, 1517.36, 1517.39, 3048.05, 3048.78, 3101.4, 3101.42, 3125.63, 3125.65], 'cm^-1'),
        ),
        HinderedRotor(
            inertia = (1.56078, 'amu*angstrom^2'),
            symmetry = 3,
            fourier = (
                [
                    [0.00483331, 0.00590149, -5.91353, 0.00405931, 0.00424175],
                    [0.00077504, -0.002541, 0.000169903, -0.000615555, 0.00090592],
                ],
                'kJ/mol',
            ),
            quantum = None,
            semiclassical = None,
        ),
    ],
    spin_multiplicity = 1,
    optical_isomers = 1,
)

# Coordinates for C2H5 in Input Orientation (angstroms):
#   C    0.8578    0.0568    0.0848
#   C   -0.6200   -0.0359   -0.0180
#   H    1.4627   -0.8311    0.1906
#   H    1.3674    0.9931   -0.0853
#   H   -1.0153   -0.8622    0.5744
#   H   -0.9458   -0.2095   -1.0529
#   H   -1.1065    0.8837    0.3103
conformer(
    label = 'C2H5',
    E0 = (109.592, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(29.0391, 'amu')),
        NonlinearRotor(
            inertia = ([4.84707, 22.1121, 23.8553], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([483.363, 816.754, 989.436, 1072.24, 1205.56, 1411.67, 1480.75, 1493.46, 1495.6, 2979.52, 3066.51, 3111.44, 3169.91, 3271.77], 'cm^-1'),
        ),
        FreeRotor(inertia=(1.11199, 'amu*angstrom^2'), symmetry=6),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
)

# Coordinates for TS1 in Input Orientation (angstroms):
#   N    0.0276    1.1818    0.3450
#   H   -0.8177    1.4794   -0.1476
#   H   -0.2910    0.9920    1.2980
#   C    0.0447   -1.2556   -0.4830
#   H    0.1444   -0.0064   -0.0982
#   H    0.8284   -1.7591    0.0736
#   H   -0.9509   -1.6095   -0.2397
#   H    0.2371   -1.2117   -1.5501
conformer(
    label = 'TS1',
    E0 = (151.345, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(32.05, 'amu')),
        NonlinearRotor(
            inertia = ([5.29874, 56.9697, 57.6076], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([41.0363, 378.039, 382.742, 572.586, 793.867, 868.802, 1180.8, 1351.25, 1424.51, 1451.82, 1494.58, 1568.77, 3078.82, 3207.86, 3212.19, 3396.69, 3492.79], 'cm^-1'),
        ),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
    frequency = (
        -1787.6,
        'cm^-1',
    ),
)

# Coordinates for TS2 in Input Orientation (angstroms):
#   N   -0.3259    1.6474    0.7949
#   H   -0.7018    1.4120    1.7167
#   H    0.6175    1.9905    0.9920
#   C   -0.3351   -0.9778   -1.3170
#   C    0.3035   -0.7461    0.0301
#   H   -1.4198   -0.8879   -1.2602
#   H    0.0212   -0.2560   -2.0520
#   H   -0.1047   -1.9778   -1.6959
#   H   -0.0359    0.4468    0.3643
#   H    1.3915   -0.7365    0.0224
#   H   -0.0761   -1.3800    0.8288
conformer(
    label = 'TS2',
    E0 = (126.602, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(46.0657, 'amu')),
        NonlinearRotor(
            inertia = ([20.0004, 104.897, 115.906], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([18.1651, 144.892, 400.771, 496.354, 758.259, 767.847, 884.821, 943.302, 1051.07, 1210.24, 1239.66, 1314.1, 1417.01, 1474.12, 1476.44, 1499.95, 1503.87, 1560.36, 3030.89, 3083.93, 3097.99, 3119.74, 3159.81, 3391.6, 3487.83], 'cm^-1'),
        ),
        HinderedRotor(
            inertia = (2.52632, 'amu*angstrom^2'),
            symmetry = 3,
            fourier = (
                [
                    [0.0369385, 0.0138552, -4.2297, 0.00812861, -0.00397751],
                    [0.036598, 0.000687098, -0.0240438, -0.0179961, 0.0269291],
                ],
                'kJ/mol',
            ),
            quantum = None,
            semiclassical = None,
        ),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
    frequency = (
        -1708.28,
        'cm^-1',
    ),
)

#   ======= =========== =========== =========== ===============
#   Temp.   k (TST)     Tunneling   k (TST+T)   Units
#   ======= =========== =========== =========== ===============
#    2500 K   2.118e+13      1.0596   2.244e+13 cm^3/(mol*s)
#    1914.89 K   4.320e+12     1.09715   4.740e+12 cm^3/(mol*s)
#    1551.72 K   1.059e+12     1.14529   1.213e+12 cm^3/(mol*s)
#    1304.35 K   2.935e+11      1.2055   3.539e+11 cm^3/(mol*s)
#    1125 K   8.885e+10     1.27973   1.137e+11 cm^3/(mol*s)
#    989.011 K   2.874e+10     1.37058   3.939e+10 cm^3/(mol*s)
#    882.353 K   9.787e+09     1.48146   1.450e+10 cm^3/(mol*s)
#    796.46 K   3.473e+09      1.6169   5.615e+09 cm^3/(mol*s)
#    725.806 K   1.274e+09      1.7829   2.271e+09 cm^3/(mol*s)
#    666.667 K   4.799e+08      1.9875   9.539e+08 cm^3/(mol*s)
#    616.438 K   1.850e+08     2.24163   4.146e+08 cm^3/(mol*s)
#    573.248 K   7.264e+07     2.56028   1.860e+08 cm^3/(mol*s)
#    535.714 K   2.898e+07      2.9643   8.591e+07 cm^3/(mol*s)
#    502.793 K   1.172e+07     3.48313   4.081e+07 cm^3/(mol*s)
#    473.684 K   4.792e+06     4.15902   1.993e+07 cm^3/(mol*s)
#    447.761 K   1.979e+06     5.05351   9.999e+06 cm^3/(mol*s)
#    424.528 K   8.239e+05     6.25781   5.156e+06 cm^3/(mol*s)
#    403.587 K   3.455e+05     7.90949   2.733e+06 cm^3/(mol*s)
#    384.615 K   1.458e+05     10.2196   1.490e+06 cm^3/(mol*s)
#    367.347 K   6.186e+04     13.5177   8.362e+05 cm^3/(mol*s)
#    351.562 K   2.637e+04     18.3276   4.833e+05 cm^3/(mol*s)
#    337.079 K   1.129e+04     25.4967   2.877e+05 cm^3/(mol*s)
#    323.741 K   4.847e+03     36.4189   1.765e+05 cm^3/(mol*s)
#    311.419 K   2.088e+03      53.427   1.116e+05 cm^3/(mol*s)
#     300 K   9.021e+02     80.4886   7.261e+04 cm^3/(mol*s)
#   ======= =========== =========== =========== ===============


#   ======= ============ =========== ============ ============= =========
#   Temp.    Kc (eq)        Units     k_rev (TST) k_rev (TST+T)   Units
#   ======= ============ =========== ============ ============= =========
#    2500 K   5.387e+00              3.931e+12     4.165e+12      cm^3/(mol*s)
#    1914.89 K   6.084e+00              7.101e+11     7.791e+11      cm^3/(mol*s)
#    1551.72 K   6.830e+00              1.551e+11     1.776e+11      cm^3/(mol*s)
#    1304.35 K   7.639e+00              3.842e+10     4.632e+10      cm^3/(mol*s)
#    1125 K   8.532e+00              1.041e+10     1.333e+10      cm^3/(mol*s)
#    989.011 K   9.528e+00              3.016e+09     4.134e+09      cm^3/(mol*s)
#    882.353 K   1.065e+01              9.186e+08     1.361e+09      cm^3/(mol*s)
#    796.46 K   1.194e+01              2.909e+08     4.704e+08      cm^3/(mol*s)
#    725.806 K   1.340e+01              9.503e+07     1.694e+08      cm^3/(mol*s)
#    666.667 K   1.509e+01              3.181e+07     6.322e+07      cm^3/(mol*s)
#    616.438 K   1.703e+01              1.086e+07     2.435e+07      cm^3/(mol*s)
#    573.248 K   1.927e+01              3.769e+06     9.649e+06      cm^3/(mol*s)
#    535.714 K   2.187e+01              1.325e+06     3.928e+06      cm^3/(mol*s)
#    502.793 K   2.487e+01              4.711e+05     1.641e+06      cm^3/(mol*s)
#    473.684 K   2.835e+01              1.690e+05     7.028e+05      cm^3/(mol*s)
#    447.761 K   3.239e+01              6.108e+04     3.087e+05      cm^3/(mol*s)
#    424.528 K   3.708e+01              2.222e+04     1.391e+05      cm^3/(mol*s)
#    403.587 K   4.252e+01              8.127e+03     6.428e+04      cm^3/(mol*s)
#    384.615 K   4.883e+01              2.986e+03     3.051e+04      cm^3/(mol*s)
#    367.347 K   5.618e+01              1.101e+03     1.489e+04      cm^3/(mol*s)
#    351.562 K   6.471e+01              4.075e+02     7.468e+03      cm^3/(mol*s)
#    337.079 K   7.464e+01              1.512e+02     3.855e+03      cm^3/(mol*s)
#    323.741 K   8.618e+01              5.624e+01     2.048e+03      cm^3/(mol*s)
#    311.419 K   9.962e+01              2.096e+01     1.120e+03      cm^3/(mol*s)
#     300 K   1.153e+02              7.827e+00     6.300e+02      cm^3/(mol*s)
#   ======= ============ =========== ============ ============= =========


# k_rev (TST) = Arrhenius(A=(2890.53,'cm^3/(mol*s)'), n=3.03416, Ea=(57.7448,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Fitted to 25 data points; dA = *|/ 1.43776, dn = +|- 0.0475538, dEa = +|- 0.263072 kJ/mol""") 
# k_rev (TST+T) = Arrhenius(A=(0.0001138,'cm^3/(mol*s)'), n=5.12165, Ea=(34.9471,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Fitted to 25 data points; dA = *|/ 4.21476, dn = +|- 0.188413, dEa = +|- 1.04232 kJ/mol""") 

# kinetics fitted using the modified three-parameter Arrhenius equation k = A * (T/T0)^n * exp(-Ea/RT) 
kinetics(
    label = 'CH4 + NH2 <=> CH3 + NH3',
    kinetics = Arrhenius(
        A = (5.37365e-05, 'cm^3/(mol*s)'),
        n = 5.37728,
        Ea = (24.9432, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
        comment = 'Fitted to 25 data points; dA = *|/ 5.06424, dn = +|- 0.21246, dEa = +|- 1.17535 kJ/mol',
    ),
)

#   ======= =========== =========== =========== ===============
#   Temp.   k (TST)     Tunneling   k (TST+T)   Units
#   ======= =========== =========== =========== ===============
#    2500 K   3.310e+13     1.05541   3.494e+13 cm^3/(mol*s)
#    1914.89 K   8.094e+12     1.08972   8.820e+12 cm^3/(mol*s)
#    1551.72 K   2.375e+12     1.13339   2.691e+12 cm^3/(mol*s)
#    1304.35 K   7.865e+11     1.18757   9.341e+11 cm^3/(mol*s)
#    1125 K   2.841e+11     1.25379   3.563e+11 cm^3/(mol*s)
#    989.011 K   1.095e+11     1.33404   1.461e+11 cm^3/(mol*s)
#    882.353 K   4.441e+10     1.43087   6.355e+10 cm^3/(mol*s)
#    796.46 K   1.874e+10     1.54762   2.900e+10 cm^3/(mol*s)
#    725.806 K   8.165e+09     1.68859   1.379e+10 cm^3/(mol*s)
#    666.667 K   3.652e+09     1.85942   6.791e+09 cm^3/(mol*s)
#    616.438 K   1.669e+09      2.0675   3.452e+09 cm^3/(mol*s)
#    573.248 K   7.771e+08     2.32264   1.805e+09 cm^3/(mol*s)
#    535.714 K   3.673e+08     2.63799   9.689e+08 cm^3/(mol*s)
#    502.793 K   1.759e+08     3.03132   5.331e+08 cm^3/(mol*s)
#    473.684 K   8.513e+07     3.52696   3.003e+08 cm^3/(mol*s)
#    447.761 K   4.161e+07     4.15858   1.730e+08 cm^3/(mol*s)
#    424.528 K   2.050e+07     4.97335   1.020e+08 cm^3/(mol*s)
#    403.587 K   1.017e+07     6.03813   6.144e+07 cm^3/(mol*s)
#    384.615 K   5.081e+06     7.44886   3.785e+07 cm^3/(mol*s)
#    367.347 K   2.551e+06     9.34495   2.384e+07 cm^3/(mol*s)
#    351.562 K   1.287e+06     11.9313   1.536e+07 cm^3/(mol*s)
#    337.079 K   6.521e+05     15.5131   1.012e+07 cm^3/(mol*s)
#    323.741 K   3.316e+05     20.5495   6.814e+06 cm^3/(mol*s)
#    311.419 K   1.692e+05     27.7404   4.694e+06 cm^3/(mol*s)
#     300 K   8.659e+04     38.1635   3.305e+06 cm^3/(mol*s)
#   ======= =========== =========== =========== ===============


#   ======= ============ =========== ============ ============= =========
#   Temp.    Kc (eq)        Units     k_rev (TST) k_rev (TST+T)   Units
#   ======= ============ =========== ============ ============= =========
#    2500 K   2.255e+01              1.468e+12     1.549e+12      cm^3/(mol*s)
#    1914.89 K   3.381e+01              2.394e+11     2.609e+11      cm^3/(mol*s)
#    1551.72 K   5.023e+01              4.728e+10     5.358e+10      cm^3/(mol*s)
#    1304.35 K   7.412e+01              1.061e+10     1.260e+10      cm^3/(mol*s)
#    1125 K   1.088e+02              2.611e+09     3.273e+09      cm^3/(mol*s)
#    989.011 K   1.593e+02              6.877e+08     9.174e+08      cm^3/(mol*s)
#    882.353 K   2.326e+02              1.910e+08     2.732e+08      cm^3/(mol*s)
#    796.46 K   3.391e+02              5.526e+07     8.553e+07      cm^3/(mol*s)
#    725.806 K   4.939e+02              1.653e+07     2.791e+07      cm^3/(mol*s)
#    666.667 K   7.190e+02              5.080e+06     9.445e+06      cm^3/(mol*s)
#    616.438 K   1.046e+03              1.596e+06     3.299e+06      cm^3/(mol*s)
#    573.248 K   1.522e+03              5.105e+05     1.186e+06      cm^3/(mol*s)
#    535.714 K   2.215e+03              1.658e+05     4.374e+05      cm^3/(mol*s)
#    502.793 K   3.222e+03              5.457e+04     1.654e+05      cm^3/(mol*s)
#    473.684 K   4.689e+03              1.816e+04     6.404e+04      cm^3/(mol*s)
#    447.761 K   6.822e+03              6.099e+03     2.536e+04      cm^3/(mol*s)
#    424.528 K   9.928e+03              2.065e+03     1.027e+04      cm^3/(mol*s)
#    403.587 K   1.445e+04              7.043e+02     4.252e+03      cm^3/(mol*s)
#    384.615 K   2.103e+04              2.416e+02     1.800e+03      cm^3/(mol*s)
#    367.347 K   3.060e+04              8.336e+01     7.790e+02      cm^3/(mol*s)
#    351.562 K   4.455e+04              2.889e+01     3.447e+02      cm^3/(mol*s)
#    337.079 K   6.485e+04              1.006e+01     1.560e+02      cm^3/(mol*s)
#    323.741 K   9.440e+04              3.513e+00     7.219e+01      cm^3/(mol*s)
#    311.419 K   1.374e+05              1.231e+00     3.415e+01      cm^3/(mol*s)
#     300 K   2.001e+05              4.328e-01     1.652e+01      cm^3/(mol*s)
#   ======= ============ =========== ============ ============= =========


# k_rev (TST) = Arrhenius(A=(158.285,'cm^3/(mol*s)'), n=3.30773, Ea=(61.6715,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Fitted to 25 data points; dA = *|/ 1.25153, dn = +|- 0.0293849, dEa = +|- 0.16256 kJ/mol""") 
# k_rev (TST+T) = Arrhenius(A=(0.00025397,'cm^3/(mol*s)'), n=4.93712, Ea=(43.1885,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Fitted to 25 data points; dA = *|/ 2.9743, dn = +|- 0.142759, dEa = +|- 0.789755 kJ/mol""") 

# kinetics fitted using the modified three-parameter Arrhenius equation k = A * (T/T0)^n * exp(-Ea/RT) 
kinetics(
    label = 'C2H6 + NH2 <=> C2H5 + NH3',
    kinetics = Arrhenius(
        A = (0.00358032, 'cm^3/(mol*s)'),
        n = 4.84415,
        Ea = (18.0398, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
        comment = 'Fitted to 25 data points; dA = *|/ 3.13894, dn = +|- 0.149815, dEa = +|- 0.828791 kJ/mol',
    ),
)

