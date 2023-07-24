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

# Coordinates for C2H4 in Input Orientation (angstroms):
#   C    0.0017    0.0000    0.0011
#   H   -0.0002    0.0000    1.0862
#   H    0.9750    0.0000   -0.4788
#   C   -1.1246    0.0000   -0.7008
#   H   -1.1210    0.0000   -1.7858
#   H   -2.0970    0.0000   -0.2195
conformer(
    label = 'C2H4',
    E0 = (43.5843, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(28.0313, 'amu')),
        NonlinearRotor(
            inertia = ([3.41526, 16.6498, 20.065], 'amu*angstrom^2'),
            symmetry = 4,
        ),
        HarmonicOscillator(
            frequencies = ([839.994, 984.241, 990.905, 1067.67, 1250.82, 1386.7, 1485.6, 1695.66, 3141.83, 3155.27, 3210.11, 3238.25], 'cm^-1'),
        ),
    ],
    spin_multiplicity = 1,
    optical_isomers = 1,
)

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
            frequencies = ([101.951, 492.752, 816.307, 982.078, 1065.92, 1194.81, 1404.38, 1469.56, 1486.34, 1487.15, 2955.21, 3046.49, 3090.43, 3150.77, 3251.19], 'cm^-1'),
        ),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
)

# Coordinates for TS in Input Orientation (angstroms):
#   C   -0.1043    0.7446    0.0000
#   C   -0.1043   -0.6010    0.0000
#   H   -0.0837    1.3123    0.9233
#   H   -0.0837    1.3123   -0.9233
#   H   -0.1986   -1.1655   -0.9205
#   H   -0.1986   -1.1655    0.9205
#   H    1.8162   -1.1553    0.0000
conformer(
    label = 'TS',
    E0 = (267.403, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(29.0391, 'amu')),
        NonlinearRotor(
            inertia = ([6.78512, 22.1437, 22.2114], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([418.528, 421.019, 832.996, 937.382, 996.472, 1038.5, 1241.35, 1344.93, 1475.43, 1622.75, 3144.88, 3154.1, 3219.79, 3246.71], 'cm^-1'),
        ),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
    frequency = (
        -760.735,
        'cm^-1',
    ),
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
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [2.5, 4.66012e-16, -1.89077e-18, 2.26862e-21, -8.121e-25, 25472.6, -0.461279],
                Tmin = (10, 'K'),
                Tmax = (1259.07, 'K'),
            ),
            NASAPolynomial(
                coeffs = [2.5, 1.42187e-14, -1.19788e-17, 4.27646e-21, -5.48837e-25, 25472.6, -0.461279],
                Tmin = (1259.07, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (211.791, 'kJ/mol'),
        Cp0 = (20.7862, 'J/(mol*K)'),
        CpInf = (20.7862, 'J/(mol*K)'),
    ),
)

# Thermodynamics for C2H4:
#   Enthalpy of formation (298 K)   =    12.912 kcal/mol
#   Entropy of formation (298 K)    =    52.274 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      10.228      12.933      52.342      -2.770
#            400      12.296      14.057      55.563      -8.168
#            500      14.431      15.394      58.537     -13.875
#            600      16.421      16.938      61.347     -19.870
#            800      19.649      20.562      66.537     -32.668
#           1000      22.093      24.747      71.196     -46.449
#           1500      25.971      36.871      80.975     -84.592
#           2000      28.007      50.415      88.754    -127.093
#           2400      29.035      61.835      93.956    -163.660
#    =========== =========== =========== =========== ===========
thermo(
    label = 'C2H4',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [4.10211, -0.0068058, 4.97293e-05, -5.81716e-08, 2.24308e-11, 5242.89, 3.22404],
                Tmin = (10, 'K'),
                Tmax = (775.817, 'K'),
            ),
            NASAPolynomial(
                coeffs = [0.436833, 0.0179601, -9.50042e-06, 2.47486e-09, -2.53793e-13, 5635.01, 18.8383],
                Tmin = (775.817, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (43.6147, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (133.032, 'J/(mol*K)'),
    ),
)

# Thermodynamics for C2H5:
#   Enthalpy of formation (298 K)   =    29.239 kcal/mol
#   Entropy of formation (298 K)    =    61.371 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      13.179      29.265      61.459      10.827
#            400      15.466      30.695      65.557       4.472
#            500      17.853      32.362      69.266      -2.271
#            600      20.027      34.258      72.717      -9.372
#            800      23.658      38.641      78.997     -24.557
#           1000      26.483      43.668      84.593     -40.926
#           1500      31.011      58.166      96.289     -86.267
#           2000      33.384      74.324     105.569    -136.814
#           2400      34.551      87.925     111.765    -180.312
#    =========== =========== =========== =========== ===========
thermo(
    label = 'C2H5',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [3.90884, 0.00976281, -1.53091e-05, 6.02494e-08, -5.61617e-11, 13157.8, 5.96352],
                Tmin = (10, 'K'),
                Tmax = (484.856, 'K'),
            ),
            NASAPolynomial(
                coeffs = [1.21225, 0.0200755, -1.02939e-05, 2.5897e-09, -2.5663e-13, 13559.5, 18.4675],
                Tmin = (484.856, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (109.403, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (157.975, 'J/(mol*K)'),
    ),
)

#   ======= =========== =========== =========== ===============
#   Temp.   k (TST)     Tunneling   k (TST+T)   Units
#   ======= =========== =========== =========== ===============
#     400 K   1.067e+12     1.34032   1.429e+12 cm^3/(mol*s)
#     500 K   2.448e+12     1.20737   2.956e+12 cm^3/(mol*s)
#     700 K   7.086e+12     1.10202   7.809e+12 cm^3/(mol*s)
#     900 K   1.405e+13     1.06105   1.491e+13 cm^3/(mol*s)
#    1100 K   2.309e+13     1.04074   2.403e+13 cm^3/(mol*s)
#    1200 K   2.832e+13     1.03423   2.929e+13 cm^3/(mol*s)
#   ======= =========== =========== =========== ===============


#   ======= ============ =========== ============ ============= =========
#   Temp.    Kc (eq)        Units     k_rev (TST) k_rev (TST+T)   Units
#   ======= ============ =========== ============ ============= =========
#     400 K   1.130e+20   cm^3/mol    9.437e-09     1.265e-08      s^-1
#     500 K   1.625e+16   cm^3/mol    1.506e-04     1.818e-04      s^-1
#     700 K   6.738e+11   cm^3/mol    1.052e+01     1.159e+01      s^-1
#     900 K   2.542e+09   cm^3/mol    5.528e+03     5.866e+03      s^-1
#    1100 K   7.478e+07   cm^3/mol    3.088e+05     3.214e+05      s^-1
#    1200 K   2.011e+07   cm^3/mol    1.408e+06     1.456e+06      s^-1
#   ======= ============ =========== ============ ============= =========


# k_rev (TST) = Arrhenius(A=(7.00354e+09,'s^-1'), n=1.0226, Ea=(157.234,'kJ/mol'), T0=(1,'K'), Tmin=(400,'K'), Tmax=(1200,'K'), comment="""Fitted to 6 data points; dA = *|/ 1.23596, dn = +|- 0.0279805, dEa = +|- 0.154308 kJ/mol""") 
# k_rev (TST+T) = Arrhenius(A=(1.31064e+09,'s^-1'), n=1.22993, Ea=(154.824,'kJ/mol'), T0=(1,'K'), Tmin=(400,'K'), Tmax=(1200,'K'), comment="""Fitted to 6 data points; dA = *|/ 1.46667, dn = +|- 0.0505844, dEa = +|- 0.278965 kJ/mol""") 

# kinetics fitted using the modified three-parameter Arrhenius equation k = A * (T/T0)^n * exp(-Ea/RT) 
kinetics(
    label = 'H + C2H4 <=> C2H5',
    kinetics = Arrhenius(
        A = (4.00615e+08, 'cm^3/(mol*s)'),
        n = 1.66398,
        Ea = (5.95613, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (400, 'K'),
        Tmax = (1200, 'K'),
        comment = 'Fitted to 6 data points; dA = *|/ 1.0977, dn = +|- 0.0123114, dEa = +|- 0.0678956 kJ/mol',
    ),
)

