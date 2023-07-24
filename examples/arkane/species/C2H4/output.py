# Coordinates for C2H4 in Input Orientation (angstroms):
#   C    0.0017    0.0000    0.0011
#   H   -0.0002    0.0000    1.0862
#   H    0.9750    0.0000   -0.4788
#   C   -1.1246    0.0000   -0.7008
#   H   -1.1210    0.0000   -1.7858
#   H   -2.0970    0.0000   -0.2195
conformer(
    label = 'C2H4',
    E0 = (45.76, 'kJ/mol'),
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

# Thermodynamics for C2H4:
#   Enthalpy of formation (298 K)   =    13.432 kcal/mol
#   Entropy of formation (298 K)    =    52.274 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      10.228      13.453      52.342      -2.250
#            400      12.296      14.577      55.563      -7.648
#            500      14.431      15.914      58.537     -13.355
#            600      16.421      17.458      61.347     -19.350
#            800      19.649      21.082      66.537     -32.148
#           1000      22.093      25.267      71.196     -45.929
#           1500      25.971      37.391      80.975     -84.072
#           2000      28.007      50.935      88.754    -126.573
#           2400      29.035      62.355      93.956    -163.140
#    =========== =========== =========== =========== ===========
thermo(
    label = 'C2H4',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [4.10211, -0.0068058, 4.97293e-05, -5.81716e-08, 2.24308e-11, 5504.57, 3.22404],
                Tmin = (10, 'K'),
                Tmax = (775.817, 'K'),
            ),
            NASAPolynomial(
                coeffs = [0.436833, 0.0179601, -9.50042e-06, 2.47486e-09, -2.53793e-13, 5896.68, 18.8383],
                Tmin = (775.817, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (45.7904, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (133.032, 'J/(mol*K)'),
    ),
)

