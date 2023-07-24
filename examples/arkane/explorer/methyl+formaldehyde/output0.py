#Thermo used:  
#CH2O  SMILES: C=O
#Thermo group additivity estimation: group(Cds-OdHH)
#CH3  SMILES: [CH3]
#Thermo library: primaryThermoLibrary + radical(CH3)
#[CH2]OC  SMILES: [CH2]OC
#Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-OsHHH) + group(Cs-OsHHH) + radical(CsJOCH3)
#[CH2]  SMILES: [CH2]
#Thermo library: primaryThermoLibrary
#[CH2]O  SMILES: [CH2]O
#Thermo group additivity estimation: group(O2s-CsH) + group(Cs-OsHHH) + radical(CsJOH)
#[CH2]OC  SMILES: [CH2]OC
#Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-OsHHH) + group(Cs-OsHHH) + radical(CsJOCH3)

#Path Reactions used:  
#CH2O(2) + CH3(1) <=> [CH2]OC(5)
#Estimated using average of templates [R_R;CsJ-HHH] + [Od_CO-HH;YJ] for rate rule [Od_CO-HH;CsJ-HHH]
#Euclidian distance = 3.0
#family: R_Addition_MultipleBond
#[CH2](7) + [CH2]O(8) <=> [CH2]OC(5)
#Estimated using an average for rate rule [carbene;R_H]
#Euclidian distance = 0
#family: 1,2_Insertion_carbene
#Ea raised from -5.1 to 0 kJ/mol.

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   2.187e+08   3.138e+08   4.490e+08   6.403e+08   9.096e+08   1.287e+09   1.812e+09
#           800   1.153e+08   1.591e+08   2.180e+08   2.967e+08   4.006e+08   5.365e+08   7.119e+08
#           600   5.512e+07   7.240e+07   9.417e+07   1.212e+08   1.542e+08   1.938e+08   2.405e+08
#           480   1.961e+07   2.445e+07   3.010e+07   3.656e+07   4.379e+07   5.171e+07   6.015e+07
#           400   5.380e+06   6.381e+06   7.462e+06   8.601e+06   9.769e+06   1.094e+07   1.207e+07
#       342.857   1.237e+06   1.405e+06   1.573e+06   1.736e+06   1.890e+06   2.032e+06   2.158e+06
#           300   2.552e+05   2.795e+05   3.022e+05   3.228e+05   3.409e+05   3.563e+05   3.691e+05
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['CH2O', 'CH3'],
    products = ['[CH2]OC'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (2.33725e+22, 'cm^3/(mol*s)'),
                n = -3.9962,
                Ea = (40.3861, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 31.6293, dn = +|- 0.467407, dEa = +|- 2.12628 kJ/mol',
            ),
            Arrhenius(
                A = (3.39592e+22, 'cm^3/(mol*s)'),
                n = -3.98507,
                Ea = (41.2636, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 26.9718, dn = +|- 0.445852, dEa = +|- 2.02822 kJ/mol',
            ),
            Arrhenius(
                A = (4.30086e+22, 'cm^3/(mol*s)'),
                n = -3.95554,
                Ea = (42.0949, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 22.4077, dn = +|- 0.420765, dEa = +|- 1.9141 kJ/mol',
            ),
            Arrhenius(
                A = (4.62141e+22, 'cm^3/(mol*s)'),
                n = -3.90412,
                Ea = (42.8606, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 18.1253, dn = +|- 0.392065, dEa = +|- 1.78354 kJ/mol',
            ),
            Arrhenius(
                A = (4.09781e+22, 'cm^3/(mol*s)'),
                n = -3.82725,
                Ea = (43.5397, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 14.2898, dn = +|- 0.35989, dEa = +|- 1.63718 kJ/mol',
            ),
            Arrhenius(
                A = (2.91782e+22, 'cm^3/(mol*s)'),
                n = -3.72145,
                Ea = (44.1107, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 11.0189, dn = +|- 0.324715, dEa = +|- 1.47716 kJ/mol',
            ),
            Arrhenius(
                A = (1.62716e+22, 'cm^3/(mol*s)'),
                n = -3.5836,
                Ea = (44.5525, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 8.3686, dn = +|- 0.287486, dEa = +|- 1.3078 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   6.177e+08   8.863e+08   1.268e+09   1.808e+09   2.569e+09   3.635e+09   5.117e+09
#           800   9.427e+07   1.300e+08   1.782e+08   2.425e+08   3.274e+08   4.385e+08   5.818e+08
#           600   1.053e+07   1.383e+07   1.799e+07   2.316e+07   2.946e+07   3.703e+07   4.595e+07
#           480   8.088e+05   1.008e+06   1.241e+06   1.508e+06   1.806e+06   2.132e+06   2.481e+06
#           400   4.636e+04   5.498e+04   6.430e+04   7.411e+04   8.418e+04   9.423e+04   1.040e+05
#       342.857   2.201e+03   2.500e+03   2.798e+03   3.088e+03   3.363e+03   3.615e+03   3.840e+03
#           300   9.359e+01   1.025e+02   1.108e+02   1.184e+02   1.250e+02   1.307e+02   1.354e+02
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#   pdepreaction(
#       reactants = ['[CH2]OC'],
#       products = ['CH2O', 'CH3'],
#       kinetics = PDepArrhenius(
#           pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
#           arrhenius = [
#               Arrhenius(
#                   A = (1.03804e+28, 's^-1'),
#                   n = -5.18132,
#                   Ea = (75.7516, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 13.958, dn = +|- 0.356711, dEa = +|- 1.62271 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (1.50823e+28, 's^-1'),
#                   n = -5.17019,
#                   Ea = (76.629, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 12.0083, dn = +|- 0.336352, dEa = +|- 1.5301 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (1.91013e+28, 's^-1'),
#                   n = -5.14066,
#                   Ea = (77.4604, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 10.1027, dn = +|- 0.312968, dEa = +|- 1.42372 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (2.0525e+28, 's^-1'),
#                   n = -5.08924,
#                   Ea = (78.2261, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 8.32385, dn = +|- 0.28676, dEa = +|- 1.3045 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (1.81996e+28, 's^-1'),
#                   n = -5.01236,
#                   Ea = (78.9052, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 6.74741, dn = +|- 0.258348, dEa = +|- 1.17525 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (1.29589e+28, 's^-1'),
#                   n = -4.90657,
#                   Ea = (79.4762, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 5.43325, dn = +|- 0.229035, dEa = +|- 1.0419 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (7.22667e+27, 's^-1'),
#                   n = -4.76872,
#                   Ea = (79.918, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 4.42322, dn = +|- 0.201203, dEa = +|- 0.915293 kJ/mol',
#               ),
#           ],
#       ),
#   )
#   
#   #   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   1.050e-07   1.542e-07   2.263e-07   3.322e-07   4.875e-07   7.156e-07   1.050e-06
#           800   2.243e-16   3.292e-16   4.833e-16   7.093e-16   1.041e-15   1.528e-15   2.243e-15
#           600   6.495e-25   9.534e-25   1.399e-24   2.054e-24   3.015e-24   4.425e-24   6.494e-24
#           480   1.882e-33   2.763e-33   4.055e-33   5.951e-33   8.735e-33   1.282e-32   1.882e-32
#           400   4.980e-42   7.310e-42   1.073e-41   1.575e-41   2.311e-41   3.392e-41   4.978e-41
#       342.857   1.214e-50   1.781e-50   2.614e-50   3.837e-50   5.632e-50   8.265e-50   1.213e-49
#           300   2.800e-59   4.110e-59   6.032e-59   8.852e-59   1.299e-58   1.906e-58   2.797e-58
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['[CH2]OC'],
    products = ['[CH2]', '[CH2]O'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (2.96662e+10, 's^-1'),
                n = -0.0959095,
                Ea = (394.873, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 29.3918, dn = +|- 0.457479, dEa = +|- 2.08112 kJ/mol',
            ),
            Arrhenius(
                A = (4.35621e+10, 's^-1'),
                n = -0.0959632,
                Ea = (394.873, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 29.3931, dn = +|- 0.457485, dEa = +|- 2.08114 kJ/mol',
            ),
            Arrhenius(
                A = (6.39793e+10, 's^-1'),
                n = -0.0960421,
                Ea = (394.874, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 29.395, dn = +|- 0.457494, dEa = +|- 2.08118 kJ/mol',
            ),
            Arrhenius(
                A = (9.39926e+10, 's^-1'),
                n = -0.0961578,
                Ea = (394.875, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 29.3978, dn = +|- 0.457507, dEa = +|- 2.08124 kJ/mol',
            ),
            Arrhenius(
                A = (1.38143e+11, 's^-1'),
                n = -0.0963277,
                Ea = (394.876, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 29.4019, dn = +|- 0.457526, dEa = +|- 2.08133 kJ/mol',
            ),
            Arrhenius(
                A = (2.03156e+11, 's^-1'),
                n = -0.0965769,
                Ea = (394.878, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 29.408, dn = +|- 0.457553, dEa = +|- 2.08145 kJ/mol',
            ),
            Arrhenius(
                A = (2.99034e+11, 's^-1'),
                n = -0.0969426,
                Ea = (394.88, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 29.4168, dn = +|- 0.457594, dEa = +|- 2.08164 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   3.222e-03   3.222e-03   3.222e-03   3.222e-03   3.222e-03   3.222e-03   3.222e-03
#           800   2.157e-11   2.157e-11   2.157e-11   2.157e-11   2.157e-11   2.157e-11   2.157e-11
#           600   1.660e-19   1.660e-19   1.660e-19   1.660e-19   1.660e-19   1.660e-19   1.660e-19
#           480   1.328e-27   1.328e-27   1.328e-27   1.328e-27   1.328e-27   1.328e-27   1.327e-27
#           400   1.067e-35   1.067e-35   1.067e-35   1.067e-35   1.067e-35   1.067e-35   1.066e-35
#       342.857   8.544e-44   8.544e-44   8.543e-44   8.543e-44   8.542e-44   8.540e-44   8.538e-44
#           300   6.828e-52   6.827e-52   6.827e-52   6.826e-52   6.825e-52   6.823e-52   6.821e-52
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['CH2O', 'CH3'],
    products = ['[CH2]', '[CH2]O'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (7.15526e+11, 'cm^3/(mol*s)'),
                n = 0.56871,
                Ea = (369.996, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.19375, dn = +|- 0.106309, dEa = +|- 0.483609 kJ/mol',
            ),
            Arrhenius(
                A = (7.15822e+11, 'cm^3/(mol*s)'),
                n = 0.568656,
                Ea = (369.996, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.19384, dn = +|- 0.106315, dEa = +|- 0.483637 kJ/mol',
            ),
            Arrhenius(
                A = (7.16256e+11, 'cm^3/(mol*s)'),
                n = 0.568578,
                Ea = (369.997, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.19399, dn = +|- 0.106324, dEa = +|- 0.483677 kJ/mol',
            ),
            Arrhenius(
                A = (7.16894e+11, 'cm^3/(mol*s)'),
                n = 0.568462,
                Ea = (369.998, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.1942, dn = +|- 0.106337, dEa = +|- 0.483737 kJ/mol',
            ),
            Arrhenius(
                A = (7.17831e+11, 'cm^3/(mol*s)'),
                n = 0.568293,
                Ea = (369.999, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.19451, dn = +|- 0.106356, dEa = +|- 0.483824 kJ/mol',
            ),
            Arrhenius(
                A = (7.19208e+11, 'cm^3/(mol*s)'),
                n = 0.568044,
                Ea = (370.001, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.19497, dn = +|- 0.106384, dEa = +|- 0.483952 kJ/mol',
            ),
            Arrhenius(
                A = (7.21233e+11, 'cm^3/(mol*s)'),
                n = 0.56768,
                Ea = (370.003, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.19563, dn = +|- 0.106425, dEa = +|- 0.484139 kJ/mol',
            ),
        ],
    ),
)

