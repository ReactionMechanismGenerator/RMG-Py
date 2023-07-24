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

# Coordinates for 1,2-butadiene in Input Orientation (angstroms):
#   C   -1.5700   -0.3138   -0.1093
#   C   -0.3985    0.6217    0.0128
#   C   -0.1016    1.5528   -0.8513
#   C    0.1888    2.4820   -1.7204
#   H   -1.2333   -1.3520   -0.1359
#   H   -2.2368   -0.2109    0.7490
#   H   -2.1391   -0.1122   -1.0147
#   H    0.2407    0.5085    0.8837
#   H   -0.2132    3.4837   -1.6285
#   H    0.8433    2.2823   -2.5600
conformer(
    label = '1,2-butadiene',
    E0 = (136.794, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(54.047, 'amu')),
        NonlinearRotor(
            inertia = ([14.6894, 119.932, 128.051], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([209.132, 335.04, 538.985, 571.018, 876.154, 884.908, 898.45, 1027.74, 1067.39, 1098.09, 1151.73, 1367.87, 1416.67, 1480.7, 1494.79, 1513.21, 2042.37, 3041.59, 3096.1, 3131.81, 3137.71, 3147.39, 3214.07], 'cm^-1'),
        ),
        HinderedRotor(
            inertia = (2.49882, 'amu*angstrom^2'),
            symmetry = 3,
            quantum = None,
            semiclassical = None,
        ),
    ],
    spin_multiplicity = 1,
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

# Coordinates for nitrosodioxaziridine in Input Orientation (angstroms):
#   O    1.9457   -0.1411    0.1302
#   N    0.8735   -0.2562    0.5458
#   N   -0.2399    0.2806   -0.5930
#   O   -1.2506    0.7770    0.2591
#   O   -1.3246   -0.5771   -0.3068
conformer(
    label = 'nitrosodioxaziridine',
    E0 = (392.61, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(75.9909, 'amu')),
        NonlinearRotor(
            inertia = ([26.3301, 137.076, 145.284], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([301.865, 311.665, 380.893, 729.846, 730.204, 856.613, 1167.13, 1769.65], 'cm^-1'),
        ),
        HinderedRotor(
            inertia = (3.01662, 'amu*angstrom^2'),
            symmetry = 1,
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

