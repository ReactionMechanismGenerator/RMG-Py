# Coordinates for NH in Input Orientation (angstroms):
#   N    0.0000    0.0000    0.1297
#   H    0.0000    0.0000   -0.9077
conformer(
    label = 'NH',
    E0 = (350.042, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(15.0109, 'amu')),
        LinearRotor(inertia=(1.01168, 'amu*angstrom^2'), symmetry=1),
        HarmonicOscillator(frequencies=([3226.97], 'cm^-1')),
    ],
    spin_multiplicity = 3,
    optical_isomers = 1,
)

# Coordinates for NH2 in Input Orientation (angstroms):
#   N    0.0000    0.0000    0.1410
#   H   -0.8043    0.0000   -0.4935
#   H    0.8043    0.0000   -0.4935
conformer(
    label = 'NH2',
    E0 = (177.216, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(16.0187, 'amu')),
        NonlinearRotor(
            inertia = ([0.709304, 1.30383, 2.01313], 'amu*angstrom^2'),
            symmetry = 2,
        ),
        HarmonicOscillator(frequencies=([1474.3, 3315.45, 3405.08], 'cm^-1')),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
)

# Coordinates for N2H3 in Input Orientation (angstroms):
#   N    0.5911    0.0244   -0.0683
#   H    1.0135    0.9034    0.1931
#   H    1.1322   -0.7941    0.1537
#   N   -0.7336   -0.1512    0.0220
#   H   -1.1480    0.7787   -0.0224
conformer(
    label = 'N2H3',
    E0 = (217.602, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(31.0296, 'amu')),
        NonlinearRotor(
            inertia = ([2.44132, 16.2349, 18.4358], 'amu*angstrom^2'),
            symmetry = 2,
        ),
        HarmonicOscillator(
            frequencies = ([503.478, 689.495, 1119.83, 1255.02, 1451.08, 1618.69, 3380.58, 3442.68, 3577.12], 'cm^-1'),
        ),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
)

# Coordinates for N2H4 in Input Orientation (angstroms):
#   N    0.7034    0.0975   -0.0730
#   N   -0.7034   -0.0975   -0.0730
#   H    1.0540    0.3872    0.8316
#   H   -1.0540   -0.3872    0.8316
#   H    1.1435   -0.7766   -0.3202
#   H   -1.1435    0.7766   -0.3202
conformer(
    label = 'N2H4',
    E0 = (89.6694, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(32.0375, 'amu')),
        NonlinearRotor(
            inertia = ([3.4483, 20.5012, 20.5139], 'amu*angstrom^2'),
            symmetry = 2,
        ),
        HarmonicOscillator(
            frequencies = ([804.703, 944.371, 1119.47, 1272.46, 1302.32, 1630.68, 1642.24, 3409.86, 3418.21, 3510.51, 3514.93], 'cm^-1'),
        ),
        HinderedRotor(
            inertia = (0.86462, 'amu*angstrom^2'),
            symmetry = 1,
            fourier = (
                [
                    [0.162229, -12.2027, -0.568401, -0.159396, 0.480391],
                    [-7.93829, -1.62775, 3.25806, 0.627004, -0.250363],
                ],
                'kJ/mol',
            ),
        ),
    ],
    spin_multiplicity = 1,
    optical_isomers = 1,
)

# Coordinates for TS1 in Input Orientation (angstroms):
#   N   -0.4465    0.6830   -0.0933
#   H   -0.4574    1.1483    0.8105
#   H    0.6774    0.3821   -0.2197
#   N   -1.2239   -0.4696   -0.0070
#   H   -1.8039   -0.5112    0.8167
#   H   -1.7837   -0.5686   -0.8405
#   N    1.9039   -0.1568   -0.0766
#   H    1.7333   -0.8469    0.6712
conformer(
    label = 'TS1',
    E0 = (482.108, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(47.0484, 'amu')),
        NonlinearRotor(
            inertia = ([15.4654, 87.3694, 97.4227], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([24.5119, 158.789, 236.806, 449.624, 640.278, 726.2, 819.01, 1128.63, 1140.18, 1322.51, 1399.92, 1511, 1617.4, 3279.34, 3414.81, 3454.27, 3552.09], 'cm^-1'),
        ),
    ],
    spin_multiplicity = 3,
    optical_isomers = 2,
    frequency = (
        -1843.03,
        'cm^-1',
    ),
)

#   ======= =========== =========== =========== ===============
#   Temp.   k (TST)     Tunneling   k (TST+T)   Units
#   ======= =========== =========== =========== ===============
#    3000 K   6.429e+12     1.04515   6.719e+12 cm^3/(mol*s)
#    2558.82 K   2.217e+12     1.05982   2.349e+12 cm^3/(mol*s)
#    2230.77 K   8.195e+11      1.0767   8.824e+11 cm^3/(mol*s)
#    1977.27 K   3.197e+11     1.09588   3.503e+11 cm^3/(mol*s)
#    1775.51 K   1.301e+11     1.11749   1.454e+11 cm^3/(mol*s)
#    1611.11 K   5.482e+10     1.14166   6.259e+10 cm^3/(mol*s)
#    1474.58 K   2.378e+10     1.16857   2.778e+10 cm^3/(mol*s)
#    1359.37 K   1.057e+10     1.19839   1.266e+10 cm^3/(mol*s)
#    1260.87 K   4.794e+09     1.23136   5.903e+09 cm^3/(mol*s)
#    1175.68 K   2.214e+09      1.2677   2.807e+09 cm^3/(mol*s)
#    1101.27 K   1.039e+09     1.30771   1.359e+09 cm^3/(mol*s)
#    1035.71 K   4.943e+08      1.3517   6.682e+08 cm^3/(mol*s)
#    977.528 K   2.380e+08     1.40003   3.333e+08 cm^3/(mol*s)
#    925.532 K   1.159e+08     1.45312   1.684e+08 cm^3/(mol*s)
#    878.788 K   5.695e+07     1.51143   8.607e+07 cm^3/(mol*s)
#    836.538 K   2.823e+07     1.57549   4.448e+07 cm^3/(mol*s)
#    798.165 K   1.410e+07     1.64592   2.321e+07 cm^3/(mol*s)
#    763.158 K   7.095e+06      1.7234   1.223e+07 cm^3/(mol*s)
#    731.092 K   3.592e+06     1.80872   6.497e+06 cm^3/(mol*s)
#    701.613 K   1.829e+06      1.9028   3.480e+06 cm^3/(mol*s)
#    674.419 K   9.362e+05     2.00668   1.879e+06 cm^3/(mol*s)
#    649.254 K   4.815e+05     2.12154   1.021e+06 cm^3/(mol*s)
#    625.899 K   2.487e+05     2.24877   5.593e+05 cm^3/(mol*s)
#    604.167 K   1.290e+05     2.38996   3.082e+05 cm^3/(mol*s)
#    583.893 K   6.713e+04     2.54696   1.710e+05 cm^3/(mol*s)
#    564.935 K   3.506e+04      2.7219   9.543e+04 cm^3/(mol*s)
#    547.17 K   1.837e+04     2.91727   5.358e+04 cm^3/(mol*s)
#    530.488 K   9.650e+03     3.13598   3.026e+04 cm^3/(mol*s)
#    514.793 K   5.084e+03     3.38142   1.719e+04 cm^3/(mol*s)
#     500 K   2.685e+03     3.65755   9.819e+03 cm^3/(mol*s)
#   ======= =========== =========== =========== ===============


#   ======= ============ =========== ============ ============= =========
#   Temp.    Kc (eq)        Units     k_rev (TST) k_rev (TST+T)   Units
#   ======= ============ =========== ============ ============= =========
#    3000 K   8.544e-02              7.525e+13     7.864e+13      cm^3/(mol*s)
#    2558.82 K   6.128e-02              3.617e+13     3.834e+13      cm^3/(mol*s)
#    2230.77 K   4.407e-02              1.860e+13     2.002e+13      cm^3/(mol*s)
#    1977.27 K   3.177e-02              1.006e+13     1.103e+13      cm^3/(mol*s)
#    1775.51 K   2.295e-02              5.669e+12     6.335e+12      cm^3/(mol*s)
#    1611.11 K   1.661e-02              3.300e+12     3.767e+12      cm^3/(mol*s)
#    1474.58 K   1.205e-02              1.974e+12     2.306e+12      cm^3/(mol*s)
#    1359.37 K   8.749e-03              1.208e+12     1.447e+12      cm^3/(mol*s)
#    1260.87 K   6.363e-03              7.534e+11     9.276e+11      cm^3/(mol*s)
#    1175.68 K   4.634e-03              4.779e+11     6.058e+11      cm^3/(mol*s)
#    1101.27 K   3.379e-03              3.076e+11     4.022e+11      cm^3/(mol*s)
#    1035.71 K   2.466e-03              2.005e+11     2.710e+11      cm^3/(mol*s)
#    977.528 K   1.802e-03              1.321e+11     1.850e+11      cm^3/(mol*s)
#    925.532 K   1.317e-03              8.796e+10     1.278e+11      cm^3/(mol*s)
#    878.788 K   9.640e-04              5.908e+10     8.929e+10      cm^3/(mol*s)
#    836.538 K   7.059e-04              3.999e+10     6.300e+10      cm^3/(mol*s)
#    798.165 K   5.173e-04              2.726e+10     4.487e+10      cm^3/(mol*s)
#    763.158 K   3.793e-04              1.870e+10     3.224e+10      cm^3/(mol*s)
#    731.092 K   2.783e-04              1.291e+10     2.335e+10      cm^3/(mol*s)
#    701.613 K   2.043e-04              8.954e+09     1.704e+10      cm^3/(mol*s)
#    674.419 K   1.500e-04              6.241e+09     1.252e+10      cm^3/(mol*s)
#    649.254 K   1.102e-04              4.369e+09     9.269e+09      cm^3/(mol*s)
#    625.899 K   8.099e-05              3.071e+09     6.905e+09      cm^3/(mol*s)
#    604.167 K   5.955e-05              2.166e+09     5.177e+09      cm^3/(mol*s)
#    583.893 K   4.379e-05              1.533e+09     3.904e+09      cm^3/(mol*s)
#    564.935 K   3.222e-05              1.088e+09     2.962e+09      cm^3/(mol*s)
#    547.17 K   2.371e-05              7.748e+08     2.260e+09      cm^3/(mol*s)
#    530.488 K   1.745e-05              5.530e+08     1.734e+09      cm^3/(mol*s)
#    514.793 K   1.285e-05              3.957e+08     1.338e+09      cm^3/(mol*s)
#     500 K   9.460e-06              2.838e+08     1.038e+09      cm^3/(mol*s)
#   ======= ============ =========== ============ ============= =========


# k_rev (TST) = Arrhenius(A=(5034.8,'cm^3/(mol*s)'), n=3.09853, Ea=(34.5325,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(3000,'K'), comment="""Fitted to 30 data points; dA = *|/ 1.05473, dn = +|- 0.00664664, dEa = +|- 0.0557567 kJ/mol""") 
# k_rev (TST+T) = Arrhenius(A=(24.7454,'cm^3/(mol*s)'), n=3.71865, Ea=(23.3052,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(3000,'K'), comment="""Fitted to 30 data points; dA = *|/ 1.35013, dn = +|- 0.0374473, dEa = +|- 0.314135 kJ/mol""") 

# kinetics fitted using the modified three-parameter Arrhenius equation k = A * (T/T0)^n * exp(-Ea/RT) 
kinetics(
    label = 'NH2 + N2H3 = NH + N2H4',
    kinetics = Arrhenius(
        A = (1.7382, 'cm^3/(mol*s)'),
        n = 3.9609,
        Ea = (66.6259, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (500, 'K'),
        Tmax = (3000, 'K'),
        comment = 'Fitted to 30 data points; dA = *|/ 1.38681, dn = +|- 0.0407906, dEa = +|- 0.34218 kJ/mol',
    ),
)

