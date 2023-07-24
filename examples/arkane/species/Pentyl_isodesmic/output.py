# Coordinates for Pentyl in Input Orientation (angstroms):
#   C    2.4514   -0.3262   -0.0025
#   C    1.1905    0.5233    0.0067
#   C   -0.0695   -0.3426   -0.0097
#   C   -1.3320    0.5173    0.0037
#   C   -2.5633   -0.3230    0.0029
#   H    3.3393    0.3135    0.0113
#   H    2.4968   -0.9513   -0.9000
#   H    2.4904   -0.9802    0.8745
#   H    1.1976    1.1872   -0.8655
#   H    1.1930    1.1598    0.8993
#   H   -0.0657   -1.0110    0.8609
#   H   -0.0643   -0.9809   -0.9025
#   H   -1.3503    1.1759   -0.8721
#   H   -1.3429    1.1607    0.8907
#   H   -3.0125   -0.6325   -0.9336
#   H   -2.9437   -0.7339    0.9309
conformer(
    label = 'Pentyl',
    E0 = (43.1689, 'kJ/mol'),
    modes = [
        IdealGasTranslation(mass=(71.0861, 'amu')),
        NonlinearRotor(
            inertia = ([28.3304, 249.55, 262.024], 'amu*angstrom^2'),
            symmetry = 1,
        ),
        HarmonicOscillator(
            frequencies = ([38.7974, 117.414, 122.004, 178.326, 244.328, 368.816, 401.221, 556.942, 733.945, 761.931, 856.81, 908.644, 983.453, 1002.12, 1075.79, 1088.61, 1121, 1196.01, 1258.82, 1265.29, 1321.8, 1326.8, 1342.56, 1392.52, 1411.35, 1462.8, 1483.42, 1484.68, 1492.15, 1494.3, 1506.87, 3031.71, 3037.79, 3046.35, 3054.11, 3058.45, 3069.67, 3091.2, 3117.83, 3125.41, 3172.13, 3271.83], 'cm^-1'),
        ),
    ],
    spin_multiplicity = 2,
    optical_isomers = 1,
)


#Isodesmic Reactions Used:
#------------------------
#Reaction 1:    17.739 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH]
#	Products:
#		1.0*[CH2]C
#		1.0*C[C](C)C
#
#Reaction 2:    17.458 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH2]
#	Products:
#		1.0*[CH2]C
#		1.0*CC(C)C
#
#Reaction 3:    15.323 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH2]
#	Products:
#		1.0*CCCCC
#		1.0*[CH]
#
#Reaction 4:    17.315 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH]
#	Products:
#		1.0*[CH2]C
#		1.0*[CH2]C(C)C
#
#Reaction 5:    15.476 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH3]
#	Products:
#		1.0*CC(C)C
#		1.0*[CH]C
#
#Reaction 6:    15.667 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH2]
#	Products:
#		1.0*[CH]C
#		1.0*C[C](C)C
#
#Reaction 7:    17.899 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH2]
#	Products:
#		1.0*[CH2]C
#		1.0*CCCC
#
#Reaction 8:    15.437 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH3]
#	Products:
#		1.0*C[C](C)C
#		1.0*C[CH2]
#
#Reaction 9:    14.763 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH2]
#	Products:
#		1.0*[CH]
#		1.0*CC(C)(C)C
#
#Reaction 10:    15.667 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH2]
#	Products:
#		1.0*[CH]C
#		1.0*C[C](C)C
#
#Reaction 11:    17.458 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH2]
#	Products:
#		1.0*[CH2]C
#		1.0*CC(C)C
#
#Reaction 12:    15.323 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH2]
#	Products:
#		1.0*CCCCC
#		1.0*[CH]
#
#Reaction 13:    17.899 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH]
#	Products:
#		1.0*[CH2]C
#		1.0*[CH2]CCC
#
#Reaction 14:    19.625 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH-]
#	Products:
#		1.0*CC(C)C
#		1.0*[CH]C
#
#Reaction 15:    15.916 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH3]
#	Products:
#		1.0*CCCC
#		1.0*[CH]C
#
#Reaction 16:    15.437 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH3]
#	Products:
#		1.0*C[C](C)C
#		1.0*C[CH2]
#
#Reaction 17:    15.476 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH3]
#	Products:
#		1.0*CC(C)C
#		1.0*[CH]C
#
#Reaction 18:    15.412 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH3]
#	Products:
#		1.0*CCCCC
#		1.0*[CH2]
#
#Reaction 19:    15.244 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH2]
#	Products:
#		1.0*[CH]C
#		1.0*[CH2]C(C)C
#
#Reaction 20:    15.916 kcal/mol
#	Reactants:
#		1.0*[CH2]CCCC
#		1.0*[CH3]
#	Products:
#		1.0*CCCC
#		1.0*[CH]C
#
#
# Thermodynamics for Pentyl:
#   Enthalpy of formation (298 K)   =    14.769 kcal/mol
#   Entropy of formation (298 K)    =    85.704 cal/(mol*K)
#    =========== =========== =========== =========== ===========
#    Temperature Heat cap.   Enthalpy    Entropy     Free energy
#    (K)         (cal/mol*K) (kcal/mol)  (cal/mol*K) (kcal/mol)
#    =========== =========== =========== =========== ===========
#            300      27.699      14.825      85.889     -10.942
#            400      34.674      17.940      94.802     -19.981
#            500      41.430      21.752     103.280     -29.888
#            600      47.371      26.198     111.371     -40.624
#            800      57.086      36.689     126.393     -64.426
#           1000      64.411      48.874     139.958     -91.084
#           1500      75.554      84.205     168.457    -168.480
#           2000      81.125     123.515     191.036    -258.558
#           2400      83.942     156.559     206.090    -338.057
#    =========== =========== =========== =========== ===========
thermo(
    label = 'Pentyl',
    thermo = NASA(
        polynomials = [
            NASAPolynomial(
                coeffs = [3.51823, 0.0511227, -0.000134204, 3.60074e-07, -3.16044e-10, 4736.33, 11.2552],
                Tmin = (10, 'K'),
                Tmax = (433.56, 'K'),
            ),
            NASAPolynomial(
                coeffs = [-1.27148, 0.0585866, -3.29658e-05, 9.02833e-09, -9.65105e-13, 5496.83, 34.3407],
                Tmin = (433.56, 'K'),
                Tmax = (3000, 'K'),
            ),
        ],
        Tmin = (10, 'K'),
        Tmax = (3000, 'K'),
        E0 = (39.3687, 'kJ/mol'),
        Cp0 = (33.2579, 'J/(mol*K)'),
        CpInf = (382.466, 'J/(mol*K)'),
    ),
)

