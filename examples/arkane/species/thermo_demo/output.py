# Coordinates for methoxy in Input Orientation (angstroms):
#   C   -0.0396    0.0295    0.0290
#   O    1.3091   -0.0148   -0.1926
#   H   -0.3091   -0.1705    1.0721
#   H   -0.5035    0.9527   -0.3359
#   H   -0.4570   -0.7969   -0.5726
conformer(
    label = 'methoxy',
    E0 = (9.20815, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(31.0184, 'amu')),
        NonlinearRotor(
            inertia = ([3.1945, 18.0652, 18.1781], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([742.343, 961.897, 1113.04, 1386.19, 1388.45, 1526.85, 2934.41, 3012.1, 3056.43], 'cm^-1'),
        ),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
)

# Coordinates for aziridine in Input Orientation (angstroms):
#   C   -0.6796   -0.3917    0.0399
#   N   -0.1846    0.8620   -0.5488
#   C    0.7735   -0.1309   -0.0393
#   H   -1.1427   -1.0682   -0.6652
#   H   -1.1487   -0.3557    1.0148
#   H   -0.2840    1.6186    0.1191
#   H    1.3008    0.0840    0.8814
#   H    1.3654   -0.6180   -0.8019
conformer(
    label = 'aziridine',
    E0 = (112.79, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(43.0422, 'amu')),
        NonlinearRotor(
            inertia = ([22.0165, 23.6557, 37.4124], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([782.498, 856.802, 875.995, 921.658, 1013.05, 1121.6, 1125.28, 1162.59, 1243.07, 1271.83, 1301.54, 1507.86, 1535.22, 3130.15, 3135.33, 3212.57, 3225.83, 3507.07], 'cm^-1'),
        ),
    ],
    spin_multiplicity = 1,
    optical_isomers = 1,
)

# Coordinates for hydrazino in Input Orientation (angstroms):
#   N   -0.4249    0.0019    0.1697
#   N    0.7889   -0.5269   -0.0829
#   H   -1.1837   -0.6451    0.0297
#   H   -0.6312    0.9585   -0.0852
#   H    1.4357    0.2625   -0.1118
conformer(
    label = 'hydrazino',
    E0 = (213.92, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(31.0296, 'amu')),
        NonlinearRotor(
            inertia = ([2.46079, 16.402, 18.5898], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([589.511, 1143.9, 1244.3, 1487.9, 1671.65, 3426.46, 3494.86, 3640.51], 'cm^-1'),
        ),
        HinderedRotor(
            inertia = (0.566181, 'amu*angstrom^2'),
            symmetry = 1,
            fourier = (
                [
                    [-3.9529, -45.0639, -0.234001, 8.50694, -0.0475748],
                    [-4.64412, 9.41636, -1.80644, -3.49548, 0.790767],
                ],
                'kJ/mol',
            ),
            quantum = None,
            semiclassical = None,
        ),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
)

# Coordinates for 1-propene-12-diol in Input Orientation (angstroms):
#   C   -1.4099   -0.0270   -0.3127
#   C    0.0213    0.2194    0.0012
#   O    0.2915    1.5162    0.3309
#   C    0.9831   -0.7041   -0.0387
#   O    2.2978   -0.3054    0.1859
#   H   -1.7247    0.5805   -1.1630
#   H   -2.0361    0.2471    0.5379
#   H   -1.5777   -1.0757   -0.5506
#   H    1.2436    1.5703    0.5014
#   H    0.8045   -1.7310   -0.3211
#   H    2.7198   -0.9282    0.7862
conformer(
    label = '1-propene-12-diol',
    E0 = (
        -355.219,
        'kJ/mol',
    ),
    modes = [
        IdealGasTranslation(mass=(74.0368, 'amu')),
        NonlinearRotor(
            inertia = ([50.6864, 131.492, 178.232], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([257.363, 271.836, 387.972, 574.571, 618.642, 740.341, 907.827, 1010.25, 1072.14, 1116.77, 1196.51, 1222.47, 1370.86, 1386.86, 1441.2, 1479.81, 1502.15, 1763.14, 3047.83, 3104.81, 3141.82, 3216.33, 3730.7, 3809.33], 'cm^-1'),
        ),
        HinderedRotor(
            inertia = (2.98048, 'amu*angstrom^2'),
            symmetry = 3,
            fourier = (
                [
                    [-0.00547243, 0.0543069, -5.80519, -0.0494455, 0.0116301],
                    [0.0101382, 0.0833978, -0.125129, 0.0756802, -0.00581562],
                ],
                'kJ/mol',
            ),
            quantum = None,
            semiclassical = None,
        ),
        HinderedRotor(
            inertia = (0.883325, 'amu*angstrom^2'),
            symmetry = 1,
            fourier = (
                [
                    [-9.98378, -9.60876, -2.16058, 0.132521, -0.0586424],
                    [-1.01986, -0.0695887, 0.206067, 0.0366542, 0.0501131],
                ],
                'kJ/mol',
            ),
            quantum = None,
            semiclassical = None,
        ),
        HinderedRotor(
            inertia = (0.505793, 'amu*angstrom^2'),
            symmetry = 1,
            fourier = (
                [
                    [-7.9719, 0.0374592, -0.364815, -0.373962, -0.000601827],
                    [-7.91287, 3.80599, 0.330654, -0.0396571, -0.0554915],
                ],
                'kJ/mol',
            ),
            quantum = None,
            semiclassical = None,
        ),
    ],
    spin_multiplicity = 1,
    optical_isomers = 2,
)

# Coordinates for hydroxyiminomethyl in Input Orientation (angstroms):
#   C   -0.9590   -0.5263   -0.0663
#   N    0.1675   -0.0986   -0.1887
#   O    0.8429    0.5123    0.9460
#   H   -1.7451   -0.6344    0.6680
#   H    1.6936    0.7470    0.5525
conformer(
    label = 'hydroxyiminomethyl',
    E0 = (245.459, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(44.0136, 'amu')),
        NonlinearRotor(
            inertia = ([5.87762, 42.9064, 48.784], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([374.292, 609.788, 728.187, 868.399, 1316.21, 1765.44, 3212.2, 3775.9], 'cm^-1'),
        ),
        HinderedRotor(
            inertia = (0.599431, 'amu*angstrom^2'),
            symmetry = 1,
            fourier = (
                [
                    [-11.8145, -7.86023, 0.613799, 0.0710175, -0.0452699],
                    [-0.00182629, -0.000878852, 0.00328663, 0.00221045, -0.00239991],
                ],
                'kJ/mol',
            ),
            quantum = None,
            semiclassical = None,
        ),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
)

# Coordinates for 2-methyl-2-propanamine in Input Orientation (angstroms):
#   C    0.6206    1.2785   -0.0064
#   C   -0.0793   -0.0750    0.0928
#   C   -1.2411   -0.1340   -0.9070
#   C    0.9115   -1.1998   -0.1981
#   N   -0.5296   -0.2355    1.4849
#   H   -0.0728    2.0900    0.2265
#   H    1.4500    1.3279    0.6984
#   H    1.0039    1.4456   -1.0138
#   H   -0.8928   -0.0137   -1.9350
#   H   -1.7594   -1.0927   -0.8371
#   H   -1.9649    0.6575   -0.7017
#   H    1.7425   -1.1631    0.5058
#   H    1.3050   -1.1193   -1.2121
#   H    0.4279   -2.1748   -0.1033
#   H   -1.0118   -1.1219    1.5930
#   H   -1.2010    0.4893    1.7176
conformer(
    label = '2-methyl-2-propanamine',
    E0 = (
        -146.64,
        'kJ/mol',
    ),
    modes = [
        IdealGasTranslation(mass=(73.0892, 'amu')),
        NonlinearRotor(
            inertia = ([109.351, 110.074, 112.084], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([337.037, 345.115, 420.33, 449.026, 449.994, 751.2, 872.41, 921.597, 940.985, 961.496, 964.415, 1017.84, 1062.76, 1134.77, 1250.06, 1282.22, 1356.69, 1403.39, 1408.7, 1432.12, 1486.59, 1490.48, 1496.13, 1507.93, 1515.79, 1524.22, 1658.98, 3028.2, 3037.08, 3042.68, 3099.58, 3100, 3105.88, 3109.5, 3125.21, 3130.41, 3466.75, 3549.45], 'cm^-1'),
        ),
        HinderedRotor(
            inertia = (3.0607, 'amu*angstrom^2'),
            symmetry = 3,
            fourier = (
                [
                    [-0.0308892, 0.0229415, -7.90065, -0.00141094, 0.0601841],
                    [0.0493114, -0.041063, 0.0199175, -0.0515469, 0.0421529],
                ],
                'kJ/mol',
            ),
            quantum = None,
            semiclassical = None,
        ),
        HinderedRotor(
            inertia = (3.05594, 'amu*angstrom^2'),
            symmetry = 3,
            fourier = (
                [
                    [-0.00235268, -0.0359168, -8.51865, 0.0569597, 0.01784],
                    [-0.000590203, -0.00271914, 0.0522713, -0.0103954, -0.0115878],
                ],
                'kJ/mol',
            ),
            quantum = None,
            semiclassical = None,
        ),
        HinderedRotor(
            inertia = (3.0607, 'amu*angstrom^2'),
            symmetry = 3,
            fourier = (
                [
                    [0.0715403, 0.0368228, -7.89987, -0.0229025, -0.0587407],
                    [-0.00769607, -0.0282783, 0.0720633, -0.0280901, -0.0124928],
                ],
                'kJ/mol',
            ),
            quantum = None,
            semiclassical = None,
        ),
        HinderedRotor(
            inertia = (1.71757, 'amu*angstrom^2'),
            symmetry = 3,
            fourier = (
                [
                    [-0.0135542, 0.0563181, -5.10933, -0.0394663, 0.0176944],
                    [-0.0242171, -0.0583021, -0.126264, -0.0201532, 0.0490503],
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

# Coordinates for ethynol in Input Orientation (angstroms):
#   C    1.1129   -0.0536   -0.0883
#   C   -0.0738    0.1373   -0.1167
#   O   -1.3543    0.4133   -0.1045
#   H    2.1581   -0.2279   -0.0671
#   H   -1.8429   -0.2691   -0.5829
conformer(
    label = 'ethynol',
    E0 = (73.1755, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(42.0106, 'amu')),
        NonlinearRotor(
            inertia = ([0.751417, 51.899, 52.6504], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([384.802, 421.302, 538.057, 628.682, 1077.52, 1253.76, 2249.36, 3496.25, 3766.8], 'cm^-1'),
        ),
    ],
    spin_multiplicity = 1,
    optical_isomers = 1,
)

# Coordinates for 2-iminoethyl in Input Orientation (angstroms):
#   C   -1.0429   -0.1733   -0.0142
#   C    0.2986    0.2182    0.0006
#   N    0.6236    1.4830   -0.0199
#   H   -1.3227   -1.2156    0.0024
#   H   -1.8148    0.5821   -0.0429
#   H    1.0518   -0.5723    0.0297
#   H    1.6375    1.5868   -0.0047
conformer(
    label = '2-iminoethyl',
    E0 = (194.056, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(42.0344, 'amu')),
        NonlinearRotor(
            inertia = ([7.97796, 45.7172, 53.6952], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([493.458, 805.229, 1005.99, 1008.08, 1097.93, 1251.74, 1368.87, 1473.77, 1513.71, 3069.56, 3172.94, 3283.63, 3457.39], 'cm^-1'),
        ),
        HinderedRotor(
            inertia = (1.34374, 'amu*angstrom^2'),
            symmetry = 2,
            fourier = (
                [
                    [0.0109211, -31.9907, 0.0124177, -1.91606, 0.0122431],
                    [-0.000117336, 0.000370217, 7.64567e-06, -0.000668938, 0.000135694],
                ],
                'kJ/mol',
            ),
            quantum = None,
            semiclassical = None,
        ),
        HinderedRotor(
            inertia = (0.511557, 'amu*angstrom^2'),
            symmetry = 1,
            fourier = (
                [
                    [-0.503883, -34.9578, -0.374585, -0.260442, -0.360165],
                    [-0.00176791, 0.00436538, 0.000137125, -0.000837771, -0.000177683],
                ],
                'kJ/mol',
            ),
            quantum = None,
            semiclassical = None,
        ),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
)

# Thermodynamics for methoxy:
#   Enthalpy of formation (298 K)   =     4.684 kcal/mol
#   Entropy of formation (298 K)    =    56.586 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300       9.795       4.704      56.651     -12.292
#            400      11.379       5.761      59.682     -18.112
#            500      13.010       6.981      62.398     -24.218
#            600      14.528       8.359      64.906     -30.584
#            800      16.992      11.524      69.441     -44.029
#           1000      18.847      15.116      73.442     -58.325
#           1500      21.743      25.348      81.699     -97.201
#           2000      23.214      36.625      88.177    -139.730
#           2400      23.935      46.063      92.477    -175.883
#    =========== =========== =========== =========== ===========
thermo(
    label = 'methoxy',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [4.07167, -0.00476973, 3.68981e-05, -4.32415e-08, 1.66496e-11, 1107.57, 5.40999],
                Tmin = (10, 'K'),
                Tmax = (774.398, 'K'),
            ),
            NASAPolynomial(
                coeffs = [1.26204, 0.01394, -7.47207e-06, 1.95474e-09, -2.00626e-13, 1416.87, 17.434],
                Tmin = (774.398, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (9.225, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (108.088, 'J/(mol*K)'),
    ),
)

# Thermodynamics for aziridine:
#   Enthalpy of formation (298 K)   =    29.578 kcal/mol
#   Entropy of formation (298 K)    =    59.644 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      12.267      29.603      59.726      11.685
#            400      15.922      31.011      63.751       5.510
#            500      19.499      32.784      67.694      -1.063
#            600      22.633      34.895      71.534      -8.025
#            800      27.378      39.924      78.736     -23.064
#           1000      30.848      45.764      85.237     -39.473
#           1500      36.050      62.652      98.860     -85.639
#           2000      38.669      81.393     109.625    -137.858
#           2400      40.072      97.153     116.805    -183.179
#    =========== =========== =========== =========== ===========
thermo(
    label = 'aziridine',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [4.18439, -0.0129536, 9.86204e-05, -1.27273e-07, 5.37393e-11, 13568.2, 6.67328],
                Tmin = (10, 'K'),
                Tmax = (732.475, 'K'),
            ),
            NASAPolynomial(
                coeffs = [-0.839572, 0.0289451, -1.68003e-05, 4.73523e-09, -5.17013e-13, 13916.2, 26.6984],
                Tmin = (732.475, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (112.846, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (182.918, 'J/(mol*K)'),
    ),
)

# Thermodynamics for hydrazino:
#   Enthalpy of formation (298 K)   =    53.828 kcal/mol
#   Entropy of formation (298 K)    =    57.271 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      11.341      53.851      57.347      36.647
#            400      12.793      55.059      60.811      30.734
#            500      14.078      56.404      63.808      24.500
#            600      15.206      57.869      66.476      17.984
#            800      17.095      61.107      71.120       4.211
#           1000      18.583      64.680      75.101     -10.421
#           1500      21.117      74.662      83.161     -50.079
#           2000      22.696      85.641      89.467     -93.293
#           2400      23.570      94.905      93.687    -129.944
#    =========== =========== =========== =========== ===========
thermo(
    label = 'hydrazino',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [3.98757, 0.000752, 2.9792e-05, -5.31982e-08, 3.07588e-11, 25693.3, 4.96394],
                Tmin = (10, 'K'),
                Tmax = (512.722, 'K'),
            ),
            NASAPolynomial(
                coeffs = [3.00718, 0.0106617, -5.81457e-06, 1.70064e-09, -2.03514e-13, 25764.1, 8.74901],
                Tmin = (512.722, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (213.622, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (103.931, 'J/(mol*K)'),
    ),
)

# Thermodynamics for 1-propene-12-diol:
#   Enthalpy of formation (298 K)   =   -80.486 kcal/mol
#   Entropy of formation (298 K)    =    77.438 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      23.551     -80.439      77.595    -103.717
#            400      28.853     -77.809      85.125    -111.859
#            500      33.091     -74.704      92.034    -120.721
#            600      36.662     -71.212      98.392    -130.247
#            800      42.222     -63.291     109.746    -151.088
#           1000      46.107     -54.435     119.612    -174.046
#           1500      51.286     -29.887     139.448    -239.059
#           2000      53.676      -3.597     154.558    -312.712
#           2400      55.102      18.168     164.474    -376.570
#    =========== =========== =========== =========== ===========
thermo(
    label = '1-propene-12-diol',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [3.89618, 0.0068893, 0.000140192, -3.14252e-07, 2.16762e-10, -42687.8, 10.8383],
                Tmin = (10, 'K'),
                Tmax = (467.501, 'K'),
            ),
            NASAPolynomial(
                coeffs = [2.14227, 0.0403926, -2.66526e-05, 8.30273e-09, -9.83016e-13, -42726, 15.8046],
                Tmin = (467.501, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (
            -354.945,
            'kJ/mol',
        ),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (245.277, 'J/(mol*K)'),
    ),
)

# Thermodynamics for hydroxyiminomethyl:
#   Enthalpy of formation (298 K)   =    61.651 kcal/mol
#   Entropy of formation (298 K)    =    62.540 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      13.454      61.678      62.630      42.889
#            400      15.601      63.134      66.802      36.413
#            500      17.327      64.784      70.478      29.545
#            600      18.621      66.585      73.757      22.331
#            800      20.466      70.507      79.385       6.998
#           1000      21.639      74.726      84.089      -9.362
#           1500      22.916      85.932      93.158     -53.806
#           2000      23.449      97.527      99.826    -102.125
#           2400      23.905     106.997     104.141    -142.942
#    =========== =========== =========== =========== ===========
thermo(
    label = 'hydroxyiminomethyl',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [3.9626, 0.00210306, 4.51447e-05, -8.36646e-08, 4.60304e-11, 29495, 6.91197],
                Tmin = (10, 'K'),
                Tmax = (603.513, 'K'),
            ),
            NASAPolynomial(
                coeffs = [3.21138, 0.0160647, -1.18823e-05, 3.9922e-09, -4.96815e-13, 29422.1, 8.80152],
                Tmin = (603.513, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (245.218, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (103.931, 'J/(mol*K)'),
    ),
)

# Thermodynamics for 2-methyl-2-propanamine:
#   Enthalpy of formation (298 K)   =   -30.004 kcal/mol
#   Entropy of formation (298 K)    =    78.193 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      29.134     -29.946      78.387     -53.462
#            400      36.329     -26.660      87.790     -61.776
#            500      42.326     -22.720      96.556     -70.998
#            600      47.557     -18.220     104.746     -81.068
#            800      56.070      -7.817     119.651    -103.538
#           1000      62.466       4.067     132.884    -128.817
#           1500      72.312      38.047     160.305    -202.411
#           2000      77.607      75.638     181.897    -288.157
#           2400      80.458     107.282     196.313    -363.869
#    =========== =========== =========== =========== ===========
thermo(
    label = '2-methyl-2-propanamine',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [3.87476, 0.00856898, 0.000199706, -4.62508e-07, 3.36983e-10, -17641.8, 9.26793],
                Tmin = (10, 'K'),
                Tmax = (437.548, 'K'),
            ),
            NASAPolynomial(
                coeffs = [1.54289, 0.0528722, -3.09731e-05, 9.02657e-09, -1.03462e-12, -17657.8, 16.0761],
                Tmin = (437.548, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (
            -146.693,
            'kJ/mol',
        ),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (365.837, 'J/(mol*K)'),
    ),
)

# Thermodynamics for ethynol:
#   Enthalpy of formation (298 K)   =    20.483 kcal/mol
#   Entropy of formation (298 K)    =    59.284 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      13.339      20.509      59.373       2.698
#            400      15.126      21.937      63.466      -3.450
#            500      16.443      23.519      66.991      -9.977
#            600      17.449      25.215      70.080     -16.833
#            800      19.067      28.874      75.333     -31.392
#           1000      20.282      32.814      79.724     -46.910
#           1500      22.226      43.490      88.355     -89.042
#           2000      23.438      54.923      94.925    -134.926
#           2400      24.160      64.451      99.265    -173.786
#    =========== =========== =========== =========== ===========
thermo(
    label = 'ethynol',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [3.9512, 0.00295995, 4.42163e-05, -9.64785e-08, 6.1581e-11, 8769.68, 5.2066],
                Tmin = (10, 'K'),
                Tmax = (546.862, 'K'),
            ),
            NASAPolynomial(
                coeffs = [4.42481, 0.0104758, -6.5163e-06, 2.0836e-09, -2.61588e-13, 8553.7, 1.70644],
                Tmin = (546.862, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (72.8979, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (108.088, 'J/(mol*K)'),
    ),
)

# Thermodynamics for 2-iminoethyl:
#   Enthalpy of formation (298 K)   =    49.336 kcal/mol
#   Entropy of formation (298 K)    =    63.042 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      13.506      49.363      63.132      30.423
#            400      16.463      50.862      67.423      23.892
#            500      19.244      52.649      71.401      16.949
#            600      21.669      54.698      75.130       9.620
#            800      25.508      59.435      81.917      -6.099
#           1000      28.353      64.836      87.931     -23.095
#           1500      32.531      80.195     100.329     -70.299
#           2000      34.475      97.001     109.986    -122.970
#           2400      35.403     110.988     116.358    -168.272
#    =========== =========== =========== =========== ===========
thermo(
    label = '2-iminoethyl',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [4.02742, -0.00164503, 5.72236e-05, -8.08628e-08, 3.64964e-11, 23337.2, 7.37022],
                Tmin = (10, 'K'),
                Tmax = (665.027, 'K'),
            ),
            NASAPolynomial(
                coeffs = [0.776356, 0.0238763, -1.37995e-05, 3.8271e-09, -4.12476e-13, 23637.6, 20.7365],
                Tmin = (665.027, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (194.04, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (149.66, 'J/(mol*K)'),
    ),
)

