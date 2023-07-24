#Thermo used:  
#CH2O  SMILES: C=O
#Thermo group additivity estimation: group(Cds-OdHH)
#CH3  SMILES: [CH3]
#Thermo library: primaryThermoLibrary + radical(CH3)
#CC[O]  SMILES: CC[O]
#Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ)
#[CH2]  SMILES: [CH2]
#Thermo library: primaryThermoLibrary
#C[O]  SMILES: C[O]
#Thermo group additivity estimation: group(O2s-CsH) + group(Cs-OsHHH) + radical(H3COJ)
#CC[O]  SMILES: CC[O]
#Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ)
#[H]  SMILES: [H]
#Thermo library: primaryThermoLibrary
#CC=O  SMILES: CC=O
#Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsH)
#CC[O]  SMILES: CC[O]
#Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ)
#C[CH]O  SMILES: C[CH]O
#Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH)
#CC[O]  SMILES: CC[O]
#Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ)
#CC[O]  SMILES: CC[O]
#Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ)
#[CH2]CO  SMILES: [CH2]CO
#Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CJCO)
#[H]  SMILES: [H]
#Thermo library: primaryThermoLibrary
#CC=O  SMILES: CC=O
#Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsH)
#C[CH]O  SMILES: C[CH]O
#Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH)
#[H]  SMILES: [H]
#Thermo library: primaryThermoLibrary
#C=CO  SMILES: C=CO
#Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH)
#C[CH]O  SMILES: C[CH]O
#Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH)
#C[CH]O  SMILES: C[CH]O
#Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH)
#[CH2]CO  SMILES: [CH2]CO
#Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CJCO)

#Path Reactions used:  
#CH2O(2) + CH3(1) <=> CC[O](6)
#Matched reaction 2817 CH2O + CH3 <=> C2H5O-2 in R_Addition_MultipleBond/training
#This reaction matched rate rule [CO-HH_O;CsJ-HHH]
#family: R_Addition_MultipleBond
#[CH2](7) + C[O](9) <=> CC[O](6)
#Estimated using an average for rate rule [carbene;R_H]
#Euclidian distance = 0
#Multiplied by reaction path degeneracy 3.0
#family: 1,2_Insertion_carbene
#Ea raised from -5.1 to 0 kJ/mol.
#[H](11) + CC=O(10) <=> CC[O](6)
#Matched reaction 2782 H + C2H4O <=> C2H5O-3 in R_Addition_MultipleBond/training
#This reaction matched rate rule [CO-CsH_O;HJ]
#family: R_Addition_MultipleBond
#C[CH]O(12) <=> CC[O](6)
#Matched reaction 322 C2H5O <=> C2H5O-2 in intra_H_migration/training
#This reaction matched rate rule [R2H_S;O_rad_out;Cs_H_out_H/NonDeC]
#family: intra_H_migration
#CC[O](6) <=> [CH2]CO(13)
#Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
#Euclidian distance = 1.0
#Multiplied by reaction path degeneracy 3.0
#family: intra_H_migration
#[H](11) + CC=O(10) <=> C[CH]O(12)
#Matched reaction 2818 H + C2H4O-5 <=> C2H5O-5 in R_Addition_MultipleBond/training
#This reaction matched rate rule [Od_CO-CsH;HJ]
#family: R_Addition_MultipleBond
#[H](11) + C=CO(14) <=> C[CH]O(12)
#Matched reaction 2816 H + C2H4O-2 <=> C2H5O-4 in R_Addition_MultipleBond/training
#This reaction matched rate rule [Cds-HH_Cds-OsH;HJ]
#family: R_Addition_MultipleBond
#C[CH]O(12) <=> [CH2]CO(13)
#Matched reaction 346 C2H5O-3 <=> C2H5O-4 in intra_H_migration/training
#This reaction matched rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_2H]
#family: intra_H_migration

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   1.858e+03   2.995e+03   4.771e+03   7.507e+03   1.166e+04   1.786e+04   2.694e+04
#           800   8.755e+00   1.293e+01   1.884e+01   2.707e+01   3.830e+01   5.330e+01   7.286e+01
#           600   2.635e-02   3.680e-02   5.052e-02   6.805e-02   8.976e-02   1.158e-01   1.459e-01
#           480   5.750e-05   7.630e-05   9.896e-05   1.253e-04   1.546e-04   1.860e-04   2.181e-04
#           400   9.999e-08   1.264e-07   1.556e-07   1.865e-07   2.176e-07   2.477e-07   2.754e-07
#       342.857   1.519e-10   1.839e-10   2.165e-10   2.483e-10   2.777e-10   3.038e-10   3.259e-10
#           300   2.135e-13   2.492e-13   2.831e-13   3.139e-13   3.405e-13   3.626e-13   3.803e-13
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['C[CH]O'],
    products = ['CC[O]'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (7.44879e+23, 's^-1'),
                n = -4.66518,
                Ea = (143.491, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.27525, dn = +|- 0.111245, dEa = +|- 0.506066 kJ/mol',
            ),
            Arrhenius(
                A = (4.11414e+23, 's^-1'),
                n = -4.50946,
                Ea = (143.84, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.51155, dn = +|- 0.124616, dEa = +|- 0.56689 kJ/mol',
            ),
            Arrhenius(
                A = (2.09797e+23, 's^-1'),
                n = -4.34372,
                Ea = (144.205, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.65578, dn = +|- 0.132173, dEa = +|- 0.601266 kJ/mol',
            ),
            Arrhenius(
                A = (9.19129e+22, 's^-1'),
                n = -4.1586,
                Ea = (144.533, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.69396, dn = +|- 0.134104, dEa = +|- 0.610051 kJ/mol',
            ),
            Arrhenius(
                A = (3.23801e+22, 's^-1'),
                n = -3.9456,
                Ea = (144.773, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.65759, dn = +|- 0.132265, dEa = +|- 0.601685 kJ/mol',
            ),
            Arrhenius(
                A = (8.68752e+21, 's^-1'),
                n = -3.69785,
                Ea = (144.877, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.61834, dn = +|- 0.130251, dEa = +|- 0.592524 kJ/mol',
            ),
            Arrhenius(
                A = (1.70949e+21, 's^-1'),
                n = -3.41074,
                Ea = (144.809, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.66522, dn = +|- 0.132653, dEa = +|- 0.603449 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   9.126e+08   1.259e+09   1.726e+09   2.347e+09   3.165e+09   4.226e+09   5.584e+09
#           800   7.156e+08   9.187e+08   1.166e+09   1.462e+09   1.810e+09   2.211e+09   2.662e+09
#           600   4.385e+08   5.226e+08   6.141e+08   7.115e+08   8.124e+08   9.143e+08   1.014e+09
#           480   1.923e+08   2.153e+08   2.379e+08   2.594e+08   2.792e+08   2.970e+08   3.124e+08
#           400   6.970e+07   7.468e+07   7.911e+07   8.294e+07   8.616e+07   8.880e+07   9.089e+07
#       342.857   2.336e+07   2.432e+07   2.511e+07   2.575e+07   2.625e+07   2.664e+07   2.692e+07
#           300   7.678e+06   7.857e+06   7.997e+06   8.103e+06   8.183e+06   8.241e+06   8.283e+06
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['CH2O', 'CH3'],
    products = ['CC[O]'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (4.095e+21, 'cm^3/(mol*s)'),
                n = -3.64985,
                Ea = (32.6838, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 3.94653, dn = +|- 0.185772, dEa = +|- 0.845097 kJ/mol',
            ),
            Arrhenius(
                A = (1.08572e+21, 'cm^3/(mol*s)'),
                n = -3.41651,
                Ea = (32.657, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 3.56842, dn = +|- 0.172144, dEa = +|- 0.783099 kJ/mol',
            ),
            Arrhenius(
                A = (2.23865e+20, 'cm^3/(mol*s)'),
                n = -3.15082,
                Ea = (32.4763, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 3.48131, dn = +|- 0.168799, dEa = +|- 0.767884 kJ/mol',
            ),
            Arrhenius(
                A = (3.59746e+19, 'cm^3/(mol*s)'),
                n = -2.85331,
                Ea = (32.1363, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 3.67873, dn = +|- 0.176264, dEa = +|- 0.80184 kJ/mol',
            ),
            Arrhenius(
                A = (4.55147e+18, 'cm^3/(mol*s)'),
                n = -2.52556,
                Ea = (31.6366, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 4.13472, dn = +|- 0.192076, dEa = +|- 0.873773 kJ/mol',
            ),
            Arrhenius(
                A = (4.61967e+17, 'cm^3/(mol*s)'),
                n = -2.17033,
                Ea = (30.9823, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 4.79932, dn = +|- 0.212246, dEa = +|- 0.965528 kJ/mol',
            ),
            Arrhenius(
                A = (3.8667e+16, 'cm^3/(mol*s)'),
                n = -1.79152,
                Ea = (30.184, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 5.59129, dn = +|- 0.232915, dEa = +|- 1.05955 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   3.116e+05   5.022e+05   8.001e+05   1.259e+06   1.955e+06   2.995e+06   4.518e+06
#           800   1.020e+04   1.506e+04   2.195e+04   3.154e+04   4.463e+04   6.211e+04   8.490e+04
#           600   2.248e+02   3.140e+02   4.311e+02   5.806e+02   7.659e+02   9.880e+02   1.245e+03
#           480   3.849e+00   5.107e+00   6.625e+00   8.385e+00   1.035e+01   1.245e+01   1.460e+01
#           400   5.543e-02   7.008e-02   8.627e-02   1.034e-01   1.207e-01   1.373e-01   1.527e-01
#       342.857   7.254e-04   8.782e-04   1.034e-03   1.186e-03   1.326e-03   1.451e-03   1.556e-03
#           300   9.038e-06   1.055e-05   1.199e-05   1.329e-05   1.441e-05   1.535e-05   1.610e-05
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#   pdepreaction(
#       reactants = ['CC[O]'],
#       products = ['C[CH]O'],
#       kinetics = PDepArrhenius(
#           pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
#           arrhenius = [
#               Arrhenius(
#                   A = (1.26743e+21, 's^-1'),
#                   n = -3.68868,
#                   Ea = (97.7096, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 1.47095, dn = +|- 0.0522215, dEa = +|- 0.237561 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (7.00032e+20, 's^-1'),
#                   n = -3.53296,
#                   Ea = (98.0591, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 1.61585, dn = +|- 0.0649346, dEa = +|- 0.295394 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (3.56975e+20, 's^-1'),
#                   n = -3.36722,
#                   Ea = (98.4241, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 1.7655, dn = +|- 0.0769203, dEa = +|- 0.349918 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (1.56392e+20, 's^-1'),
#                   n = -3.1821,
#                   Ea = (98.7523, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 1.91985, dn = +|- 0.0882623, dEa = +|- 0.401514 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (5.50956e+19, 's^-1'),
#                   n = -2.9691,
#                   Ea = (98.9919, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 2.11641, dn = +|- 0.101452, dEa = +|- 0.461517 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (1.47821e+19, 's^-1'),
#                   n = -2.72135,
#                   Ea = (99.0961, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 2.41788, dn = +|- 0.119473, dEa = +|- 0.543494 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (2.90875e+18, 's^-1'),
#                   n = -2.43424,
#                   Ea = (99.0275, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 2.89376, dn = +|- 0.143785, dEa = +|- 0.654092 kJ/mol',
#               ),
#           ],
#       ),
#   )
#   
#   #   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   9.970e+07   1.221e+08   1.477e+08   1.761e+08   2.070e+08   2.396e+08   2.728e+08
#           800   1.529e+07   1.655e+07   1.766e+07   1.858e+07   1.924e+07   1.962e+07   1.968e+07
#           600   1.450e+06   1.456e+06   1.441e+06   1.405e+06   1.347e+06   1.268e+06   1.169e+06
#           480   1.345e+05   1.281e+05   1.199e+05   1.100e+05   9.868e+04   8.647e+04   7.390e+04
#           400   1.307e+04   1.188e+04   1.056e+04   9.148e+03   7.725e+03   6.351e+03   5.085e+03
#       342.857   1.337e+03   1.162e+03   9.835e+02   8.099e+02   6.490e+02   5.067e+02   3.860e+02
#           300   1.413e+02   1.178e+02   9.555e+01   7.540e+01   5.799e+01   4.356e+01   3.205e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['CH2O', 'CH3'],
    products = ['C[CH]O'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (8.23447e+12, 'cm^3/(mol*s)'),
                n = -0.886674,
                Ea = (49.4228, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 37.2938, dn = +|- 0.4897, dEa = +|- 2.22769 kJ/mol',
            ),
            Arrhenius(
                A = (1.31598e+12, 'cm^3/(mol*s)'),
                n = -0.600979,
                Ea = (49.3538, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 30.8182, dn = +|- 0.463892, dEa = +|- 2.11029 kJ/mol',
            ),
            Arrhenius(
                A = (2.28443e+11, 'cm^3/(mol*s)'),
                n = -0.327219,
                Ea = (49.3943, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 25.6937, dn = +|- 0.439283, dEa = +|- 1.99834 kJ/mol',
            ),
            Arrhenius(
                A = (4.10171e+10, 'cm^3/(mol*s)'),
                n = -0.0589651,
                Ea = (49.5133, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 22.1415, dn = +|- 0.419148, dEa = +|- 1.90674 kJ/mol',
            ),
            Arrhenius(
                A = (7.16371e+09, 'cm^3/(mol*s)'),
                n = 0.211805,
                Ea = (49.669, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 20.2041, dn = +|- 0.406757, dEa = +|- 1.85038 kJ/mol',
            ),
            Arrhenius(
                A = (1.13653e+09, 'cm^3/(mol*s)'),
                n = 0.493978,
                Ea = (49.8117, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 19.9357, dn = +|- 0.404947, dEa = +|- 1.84214 kJ/mol',
            ),
            Arrhenius(
                A = (1.52883e+08, 'cm^3/(mol*s)'),
                n = 0.796434,
                Ea = (49.889, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 21.5558, dn = +|- 0.41552, dEa = +|- 1.89024 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   2.562e+08   3.536e+08   4.846e+08   6.591e+08   8.886e+08   1.187e+09   1.568e+09
#           800   1.917e+07   2.461e+07   3.124e+07   3.918e+07   4.850e+07   5.923e+07   7.130e+07
#           600   9.144e+05   1.090e+06   1.281e+06   1.484e+06   1.694e+06   1.906e+06   2.115e+06
#           480   2.945e+04   3.298e+04   3.644e+04   3.973e+04   4.276e+04   4.548e+04   4.785e+04
#           400   7.719e+02   8.270e+02   8.761e+02   9.185e+02   9.542e+02   9.834e+02   1.007e+03
#       342.857   1.870e+01   1.947e+01   2.010e+01   2.061e+01   2.101e+01   2.132e+01   2.155e+01
#           300   4.469e-01   4.574e-01   4.655e-01   4.717e-01   4.763e-01   4.797e-01   4.821e-01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#   pdepreaction(
#       reactants = ['CC[O]'],
#       products = ['CH2O', 'CH3'],
#       kinetics = PDepArrhenius(
#           pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
#           arrhenius = [
#               Arrhenius(
#                   A = (1.85123e+26, 's^-1'),
#                   n = -4.55284,
#                   Ea = (88.196, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 3.76831, dn = +|- 0.179519, dEa = +|- 0.81665 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (4.90823e+25, 's^-1'),
#                   n = -4.3195,
#                   Ea = (88.1692, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 4.44178, dn = +|- 0.20177, dEa = +|- 0.91787 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (1.01203e+25, 's^-1'),
#                   n = -4.05381,
#                   Ea = (87.9884, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 5.49947, dn = +|- 0.230674, dEa = +|- 1.04936 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (1.62631e+24, 's^-1'),
#                   n = -3.7563,
#                   Ea = (87.6485, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 7.00045, dn = +|- 0.26333, dEa = +|- 1.19791 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (2.05759e+23, 's^-1'),
#                   n = -3.42855,
#                   Ea = (87.1488, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 8.9795, dn = +|- 0.29702, dEa = +|- 1.35117 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (2.08842e+22, 's^-1'),
#                   n = -3.07332,
#                   Ea = (86.4945, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 11.4029, dn = +|- 0.329352, dEa = +|- 1.49825 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (1.74802e+21, 's^-1'),
#                   n = -2.69451,
#                   Ea = (85.6962, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 14.1192, dn = +|- 0.358265, dEa = +|- 1.62978 kJ/mol',
#               ),
#           ],
#       ),
#   )
#   
#   #   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   1.669e+05   2.045e+05   2.472e+05   2.949e+05   3.466e+05   4.012e+05   4.568e+05
#           800   3.515e+02   3.805e+02   4.061e+02   4.270e+02   4.424e+02   4.511e+02   4.523e+02
#           600   3.543e-01   3.558e-01   3.522e-01   3.434e-01   3.292e-01   3.098e-01   2.857e-01
#           480   3.079e-04   2.932e-04   2.743e-04   2.516e-04   2.258e-04   1.979e-04   1.691e-04
#           400   2.611e-07   2.373e-07   2.109e-07   1.828e-07   1.543e-07   1.269e-07   1.016e-07
#       342.857   2.241e-10   1.948e-10   1.648e-10   1.357e-10   1.088e-10   8.492e-11   6.469e-11
#           300   1.943e-13   1.620e-13   1.314e-13   1.037e-13   7.974e-14   5.990e-14   4.407e-14
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#   pdepreaction(
#       reactants = ['C[CH]O'],
#       products = ['CH2O', 'CH3'],
#       kinetics = PDepArrhenius(
#           pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
#           arrhenius = [
#               Arrhenius(
#                   A = (2.18778e+20, 's^-1'),
#                   n = -2.76616,
#                   Ea = (150.716, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 58.6884, dn = +|- 0.551057, dEa = +|- 2.50681 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (3.49638e+19, 's^-1'),
#                   n = -2.48047,
#                   Ea = (150.647, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 48.4387, dn = +|- 0.525083, dEa = +|- 2.38865 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (6.06941e+18, 's^-1'),
#                   n = -2.20671,
#                   Ea = (150.688, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 40.3958, dn = +|- 0.500512, dEa = +|- 2.27688 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (1.08977e+18, 's^-1'),
#                   n = -1.93845,
#                   Ea = (150.807, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 34.8832, dn = +|- 0.480658, dEa = +|- 2.18656 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (1.9033e+17, 's^-1'),
#                   n = -1.66768,
#                   Ea = (150.962, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 31.946, dn = +|- 0.468755, dEa = +|- 2.13241 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (3.01959e+16, 's^-1'),
#                   n = -1.38551,
#                   Ea = (151.105, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 31.6567, dn = +|- 0.467524, dEa = +|- 2.12681 kJ/mol',
#               ),
#               Arrhenius(
#                   A = (4.06189e+15, 's^-1'),
#                   n = -1.08305,
#                   Ea = (151.182, 'kJ/mol'),
#                   T0 = (1, 'K'),
#                   Tmin = (300, 'K'),
#                   Tmax = (1200, 'K'),
#                   comment = 'Fitted to 7 data points; dA = *|/ 34.3613, dn = +|- 0.478618, dEa = +|- 2.17728 kJ/mol',
#               ),
#           ],
#       ),
#   )
#   
#   #   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   3.291e-09   4.831e-09   7.091e-09   1.041e-08   1.527e-08   2.242e-08   3.290e-08
#           800   5.383e-19   7.901e-19   1.160e-18   1.702e-18   2.498e-18   3.666e-18   5.379e-18
#           600   1.107e-28   1.624e-28   2.384e-28   3.499e-28   5.135e-28   7.535e-28   1.106e-27
#           480   2.058e-38   3.021e-38   4.433e-38   6.506e-38   9.548e-38   1.401e-37   2.055e-37
#           400   3.312e-48   4.861e-48   7.134e-48   1.047e-47   1.536e-47   2.254e-47   3.306e-47
#       342.857   4.972e-58   7.297e-58   1.071e-57   1.571e-57   2.305e-57   3.382e-57   4.959e-57
#           300   6.859e-68   1.007e-67   1.477e-67   2.168e-67   3.180e-67   4.664e-67   6.838e-67
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['CC[O]'],
    products = ['[CH2]', 'C[O]'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (3.62851e+14, 's^-1'),
                n = -1.0784,
                Ea = (453.842, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 28.8358, dn = +|- 0.454895, dEa = +|- 2.06936 kJ/mol',
            ),
            Arrhenius(
                A = (5.3294e+14, 's^-1'),
                n = -1.07848,
                Ea = (453.843, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 28.8338, dn = +|- 0.454885, dEa = +|- 2.06932 kJ/mol',
            ),
            Arrhenius(
                A = (7.82998e+14, 's^-1'),
                n = -1.0786,
                Ea = (453.844, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 28.8308, dn = +|- 0.454871, dEa = +|- 2.06925 kJ/mol',
            ),
            Arrhenius(
                A = (1.1509e+15, 's^-1'),
                n = -1.07878,
                Ea = (453.846, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 28.8264, dn = +|- 0.45485, dEa = +|- 2.06916 kJ/mol',
            ),
            Arrhenius(
                A = (1.69276e+15, 's^-1'),
                n = -1.07905,
                Ea = (453.849, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 28.8199, dn = +|- 0.45482, dEa = +|- 2.06902 kJ/mol',
            ),
            Arrhenius(
                A = (2.49212e+15, 's^-1'),
                n = -1.07944,
                Ea = (453.852, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 28.8104, dn = +|- 0.454775, dEa = +|- 2.06882 kJ/mol',
            ),
            Arrhenius(
                A = (3.67409e+15, 's^-1'),
                n = -1.08,
                Ea = (453.858, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 28.7964, dn = +|- 0.454709, dEa = +|- 2.06852 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   8.820e-12   1.294e-11   1.900e-11   2.788e-11   4.091e-11   6.002e-11   8.804e-11
#           800   4.125e-22   6.053e-22   8.882e-22   1.303e-21   1.912e-21   2.804e-21   4.110e-21
#           600   1.281e-32   1.879e-32   2.757e-32   4.045e-32   5.931e-32   8.694e-32   1.274e-31
#           480   3.275e-43   4.805e-43   7.048e-43   1.033e-42   1.515e-42   2.219e-42   3.247e-42
#           400   6.554e-54   9.614e-54   1.410e-53   2.067e-53   3.028e-53   4.433e-53   6.480e-53
#       342.857   1.106e-64   1.622e-64   2.378e-64   3.485e-64   5.105e-64   7.469e-64   1.091e-63
#           300   1.791e-75   2.627e-75   3.851e-75   5.643e-75   8.261e-75   1.208e-74   1.763e-74
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['C[CH]O'],
    products = ['[CH2]', 'C[O]'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (4.51005e+23, 's^-1'),
                n = -4.11482,
                Ea = (506.565, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 4.01388, dn = +|- 0.188062, dEa = +|- 0.855513 kJ/mol',
            ),
            Arrhenius(
                A = (6.62546e+23, 's^-1'),
                n = -4.11491,
                Ea = (506.568, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 4.01153, dn = +|- 0.187983, dEa = +|- 0.855153 kJ/mol',
            ),
            Arrhenius(
                A = (9.73687e+23, 's^-1'),
                n = -4.11505,
                Ea = (506.572, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 4.00809, dn = +|- 0.187867, dEa = +|- 0.854624 kJ/mol',
            ),
            Arrhenius(
                A = (1.43176e+24, 's^-1'),
                n = -4.11524,
                Ea = (506.578, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 4.00304, dn = +|- 0.187696, dEa = +|- 0.853848 kJ/mol',
            ),
            Arrhenius(
                A = (2.10704e+24, 's^-1'),
                n = -4.11552,
                Ea = (506.588, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 3.99564, dn = +|- 0.187446, dEa = +|- 0.85271 kJ/mol',
            ),
            Arrhenius(
                A = (3.10445e+24, 's^-1'),
                n = -4.11592,
                Ea = (506.601, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 3.98481, dn = +|- 0.187079, dEa = +|- 0.851039 kJ/mol',
            ),
            Arrhenius(
                A = (4.58155e+24, 's^-1'),
                n = -4.11649,
                Ea = (506.62, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 3.96898, dn = +|- 0.18654, dEa = +|- 0.848589 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   1.149e-04   1.149e-04   1.149e-04   1.149e-04   1.149e-04   1.149e-04   1.149e-04
#           800   1.762e-13   1.762e-13   1.762e-13   1.762e-13   1.761e-13   1.761e-13   1.761e-13
#           600   2.894e-22   2.894e-22   2.894e-22   2.893e-22   2.893e-22   2.892e-22   2.891e-22
#           480   4.768e-31   4.767e-31   4.767e-31   4.766e-31   4.765e-31   4.763e-31   4.761e-31
#           400   7.536e-40   7.535e-40   7.534e-40   7.533e-40   7.530e-40   7.526e-40   7.521e-40
#       342.857   1.182e-48   1.182e-48   1.182e-48   1.182e-48   1.181e-48   1.180e-48   1.179e-48
#           300   1.762e-57   1.762e-57   1.761e-57   1.761e-57   1.760e-57   1.759e-57   1.757e-57
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['CH2O', 'CH3'],
    products = ['[CH2]', 'C[O]'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (1.34184e+14, 'cm^3/(mol*s)'),
                n = -0.146268,
                Ea = (404.947, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.88032, dn = +|- 0.143155, dEa = +|- 0.651227 kJ/mol',
            ),
            Arrhenius(
                A = (1.34271e+14, 'cm^3/(mol*s)'),
                n = -0.146351,
                Ea = (404.948, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.88011, dn = +|- 0.143145, dEa = +|- 0.651182 kJ/mol',
            ),
            Arrhenius(
                A = (1.34398e+14, 'cm^3/(mol*s)'),
                n = -0.146473,
                Ea = (404.95, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.8798, dn = +|- 0.143131, dEa = +|- 0.651116 kJ/mol',
            ),
            Arrhenius(
                A = (1.34585e+14, 'cm^3/(mol*s)'),
                n = -0.146651,
                Ea = (404.951, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.87934, dn = +|- 0.143109, dEa = +|- 0.651018 kJ/mol',
            ),
            Arrhenius(
                A = (1.3486e+14, 'cm^3/(mol*s)'),
                n = -0.146913,
                Ea = (404.954, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.87867, dn = +|- 0.143078, dEa = +|- 0.650875 kJ/mol',
            ),
            Arrhenius(
                A = (1.35263e+14, 'cm^3/(mol*s)'),
                n = -0.147297,
                Ea = (404.958, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.87769, dn = +|- 0.143032, dEa = +|- 0.650665 kJ/mol',
            ),
            Arrhenius(
                A = (1.35856e+14, 'cm^3/(mol*s)'),
                n = -0.147858,
                Ea = (404.963, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.87625, dn = +|- 0.142964, dEa = +|- 0.650356 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   1.627e+07   2.325e+07   3.298e+07   4.640e+07   6.468e+07   8.922e+07   1.216e+08
#           800   8.499e+05   1.172e+06   1.592e+06   2.129e+06   2.796e+06   3.607e+06   4.565e+06
#           600   3.448e+04   4.499e+04   5.743e+04   7.167e+04   8.744e+04   1.043e+05   1.218e+05
#           480   7.127e+02   8.842e+02   1.069e+03   1.261e+03   1.452e+03   1.635e+03   1.803e+03
#           400   1.144e+01   1.361e+01   1.576e+01   1.780e+01   1.965e+01   2.128e+01   2.264e+01
#       342.857   1.931e-01   2.207e-01   2.460e-01   2.683e-01   2.871e-01   3.024e-01   3.145e-01
#           300   4.439e-03   4.891e-03   5.277e-03   5.594e-03   5.845e-03   6.037e-03   6.181e-03
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['CC[O]'],
    products = ['[H]', 'CC=O'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (1.56115e+25, 's^-1'),
                n = -4.4985,
                Ea = (94.6625, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 312.76, dn = +|- 0.777474, dEa = +|- 3.5368 kJ/mol',
            ),
            Arrhenius(
                A = (1.95507e+25, 's^-1'),
                n = -4.46855,
                Ea = (95.4182, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 328.307, dn = +|- 0.784038, dEa = +|- 3.56666 kJ/mol',
            ),
            Arrhenius(
                A = (1.71099e+25, 's^-1'),
                n = -4.39122,
                Ea = (96.0099, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 356.492, dn = +|- 0.795183, dEa = +|- 3.61736 kJ/mol',
            ),
            Arrhenius(
                A = (1.01674e+25, 's^-1'),
                n = -4.2631,
                Ea = (96.4052, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 398.782, dn = +|- 0.810353, dEa = +|- 3.68637 kJ/mol',
            ),
            Arrhenius(
                A = (4.06104e+24, 's^-1'),
                n = -4.08326,
                Ea = (96.5826, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 456.133, dn = +|- 0.828536, dEa = +|- 3.76909 kJ/mol',
            ),
            Arrhenius(
                A = (1.09864e+24, 's^-1'),
                n = -3.85311,
                Ea = (96.5324, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 527.84, dn = +|- 0.848294, dEa = +|- 3.85897 kJ/mol',
            ),
            Arrhenius(
                A = (2.06237e+23, 's^-1'),
                n = -3.57627,
                Ea = (96.2558, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 610.017, dn = +|- 0.867874, dEa = +|- 3.94804 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   2.490e+06   3.126e+06   3.869e+06   4.719e+06   5.670e+06   6.709e+06   7.817e+06
#           800   5.172e+03   5.777e+03   6.358e+03   6.900e+03   7.393e+03   7.828e+03   8.201e+03
#           600   4.740e+00   4.958e+00   5.137e+00   5.280e+00   5.392e+00   5.477e+00   5.540e+00
#           480   2.003e-03   2.037e-03   2.061e-03   2.078e-03   2.090e-03   2.097e-03   2.101e-03
#           400   9.680e-07   9.735e-07   9.767e-07   9.781e-07   9.783e-07   9.777e-07   9.767e-07
#       342.857   8.253e-10   8.261e-10   8.259e-10   8.249e-10   8.237e-10   8.223e-10   8.209e-10
#           300   4.357e-13   4.334e-13   4.309e-13   4.284e-13   4.261e-13   4.241e-13   4.225e-13
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['C[CH]O'],
    products = ['[H]', 'CC=O'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (2.30686e+25, 's^-1'),
                n = -3.8398,
                Ea = (162.445, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 19049.2, dn = +|- 1.33355, dEa = +|- 6.06645 kJ/mol',
            ),
            Arrhenius(
                A = (1.17633e+24, 's^-1'),
                n = -3.40598,
                Ea = (161.192, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 16742.6, dn = +|- 1.31608, dEa = +|- 5.98699 kJ/mol',
            ),
            Arrhenius(
                A = (6.13771e+22, 's^-1'),
                n = -2.978,
                Ea = (159.911, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 14023.4, dn = +|- 1.2921, dEa = +|- 5.87789 kJ/mol',
            ),
            Arrhenius(
                A = (3.40414e+21, 's^-1'),
                n = -2.56097,
                Ea = (158.625, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 11280.6, dn = +|- 1.26265, dEa = +|- 5.74392 kJ/mol',
            ),
            Arrhenius(
                A = (2.08022e+20, 's^-1'),
                n = -2.15966,
                Ea = (157.354, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 8797.77, dn = +|- 1.22901, dEa = +|- 5.59089 kJ/mol',
            ),
            Arrhenius(
                A = (1.44791e+19, 's^-1'),
                n = -1.77847,
                Ea = (156.119, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 6722.26, dn = +|- 1.1926, dEa = +|- 5.42525 kJ/mol',
            ),
            Arrhenius(
                A = (1.18275e+18, 's^-1'),
                n = -1.42134,
                Ea = (154.938, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 5086.79, dn = +|- 1.15488, dEa = +|- 5.25364 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   6.170e+09   6.134e+09   6.087e+09   6.025e+09   5.946e+09   5.844e+09   5.715e+09
#           800   6.352e+08   6.232e+08   6.079e+08   5.886e+08   5.650e+08   5.365e+08   5.031e+08
#           600   8.961e+07   8.491e+07   7.936e+07   7.303e+07   6.605e+07   5.860e+07   5.094e+07
#           480   1.086e+07   9.813e+06   8.683e+06   7.516e+06   6.358e+06   5.252e+06   4.237e+06
#           400   1.315e+06   1.131e+06   9.481e+05   7.753e+05   6.185e+05   4.817e+05   3.668e+05
#       342.857   1.804e+05   1.478e+05   1.181e+05   9.194e+04   6.996e+04   5.212e+04   3.812e+04
#           300   3.326e+04   2.605e+04   1.991e+04   1.489e+04   1.092e+04   7.875e+03   5.602e+03
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['CH2O', 'CH3'],
    products = ['[H]', 'CC=O'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (5.17077e+07, 'cm^3/(mol*s)'),
                n = 1.17568,
                Ea = (35.3404, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 35.0989, dn = +|- 0.481492, dEa = +|- 2.19035 kJ/mol',
            ),
            Arrhenius(
                A = (1.36962e+08, 'cm^3/(mol*s)'),
                n = 1.05716,
                Ea = (36.7046, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 38.8293, dn = +|- 0.49516, dEa = +|- 2.25253 kJ/mol',
            ),
            Arrhenius(
                A = (3.00734e+08, 'cm^3/(mol*s)'),
                n = 0.964641,
                Ea = (38.0378, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 45.3418, dn = +|- 0.516142, dEa = +|- 2.34798 kJ/mol',
            ),
            Arrhenius(
                A = (5.07511e+08, 'cm^3/(mol*s)'),
                n = 0.907964,
                Ea = (39.2852, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 56.1785, dn = +|- 0.545142, dEa = +|- 2.4799 kJ/mol',
            ),
            Arrhenius(
                A = (6.13729e+08, 'cm^3/(mol*s)'),
                n = 0.896142,
                Ea = (40.3931, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 74.0817, dn = +|- 0.582576, dEa = +|- 2.65019 kJ/mol',
            ),
            Arrhenius(
                A = (5.00836e+08, 'cm^3/(mol*s)'),
                n = 0.936805,
                Ea = (41.3121, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 103.823, dn = +|- 0.628249, dEa = +|- 2.85797 kJ/mol',
            ),
            Arrhenius(
                A = (2.63475e+08, 'cm^3/(mol*s)'),
                n = 1.03562,
                Ea = (42.0002, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 153.426, dn = +|- 0.681096, dEa = +|- 3.09837 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   5.536e+03   8.045e+03   1.165e+04   1.682e+04   2.417e+04   3.457e+04   4.919e+04
#           800   4.351e+00   6.184e+00   8.737e+00   1.228e+01   1.718e+01   2.396e+01   3.335e+01
#           600   4.415e-03   6.162e-03   8.568e-03   1.189e-02   1.651e-02   2.294e-02   3.188e-02
#           480   4.470e-06   6.142e-06   8.450e-06   1.166e-05   1.614e-05   2.237e-05   3.096e-05
#           400   4.140e-09   5.674e-09   7.807e-09   1.079e-08   1.493e-08   2.063e-08   2.832e-08
#       342.857   3.677e-12   5.062e-12   6.995e-12   9.683e-12   1.338e-11   1.838e-11   2.494e-11
#           300   3.214e-15   4.433e-15   6.132e-15   8.480e-15   1.167e-14   1.590e-14   2.130e-14
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['CC[O]'],
    products = ['[CH2]CO'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (5.83388e+09, 's^-1'),
                n = 0.00483944,
                Ea = (139.27, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 12.3051, dn = +|- 0.339655, dEa = +|- 1.54512 kJ/mol',
            ),
            Arrhenius(
                A = (3.58874e+09, 's^-1'),
                n = 0.121491,
                Ea = (138.919, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 10.8862, dn = +|- 0.323076, dEa = +|- 1.4697 kJ/mol',
            ),
            Arrhenius(
                A = (2.17299e+09, 's^-1'),
                n = 0.239198,
                Ea = (138.532, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 10.3771, dn = +|- 0.316596, dEa = +|- 1.44022 kJ/mol',
            ),
            Arrhenius(
                A = (1.43773e+09, 's^-1'),
                n = 0.34412,
                Ea = (138.18, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 10.7528, dn = +|- 0.321408, dEa = +|- 1.46211 kJ/mol',
            ),
            Arrhenius(
                A = (1.14509e+09, 's^-1'),
                n = 0.423585,
                Ea = (137.937, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 12.011, dn = +|- 0.336382, dEa = +|- 1.53023 kJ/mol',
            ),
            Arrhenius(
                A = (1.18165e+09, 's^-1'),
                n = 0.468136,
                Ea = (137.868, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 14.1195, dn = +|- 0.358268, dEa = +|- 1.62979 kJ/mol',
            ),
            Arrhenius(
                A = (1.64328e+09, 's^-1'),
                n = 0.472939,
                Ea = (138.02, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 16.8813, dn = +|- 0.382442, dEa = +|- 1.73977 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   9.576e+03   1.345e+04   1.866e+04   2.550e+04   3.426e+04   4.519e+04   5.843e+04
#           800   3.653e+00   4.799e+00   6.155e+00   7.695e+00   9.374e+00   1.113e+01   1.288e+01
#           600   6.223e-04   7.650e-04   9.135e-04   1.061e-03   1.199e-03   1.324e-03   1.432e-03
#           480   7.853e-08   9.046e-08   1.014e-07   1.109e-07   1.189e-07   1.253e-07   1.302e-07
#           400   7.699e-12   8.537e-12   9.245e-12   9.818e-12   1.026e-11   1.060e-11   1.085e-11
#       342.857   6.639e-16   7.244e-16   7.736e-16   8.121e-16   8.410e-16   8.623e-16   8.777e-16
#           300   5.724e-20   6.133e-20   6.453e-20   6.694e-20   6.871e-20   6.998e-20   7.088e-20
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['C[CH]O'],
    products = ['[CH2]CO'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (5.54968e+29, 's^-1'),
                n = -5.48625,
                Ea = (203.436, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 5.6389, dn = +|- 0.234062, dEa = +|- 1.06477 kJ/mol',
            ),
            Arrhenius(
                A = (1.45037e+29, 's^-1'),
                n = -5.24985,
                Ea = (203.298, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 7.70574, dn = +|- 0.276319, dEa = +|- 1.257 kJ/mol',
            ),
            Arrhenius(
                A = (2.45699e+28, 's^-1'),
                n = -4.95763,
                Ea = (202.916, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 10.4545, dn = +|- 0.317601, dEa = +|- 1.4448 kJ/mol',
            ),
            Arrhenius(
                A = (2.82074e+27, 's^-1'),
                n = -4.61617,
                Ea = (202.295, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 13.5775, dn = +|- 0.352971, dEa = +|- 1.6057 kJ/mol',
            ),
            Arrhenius(
                A = (2.35537e+26, 's^-1'),
                n = -4.23543,
                Ea = (201.46, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 16.4471, dn = +|- 0.378917, dEa = +|- 1.72373 kJ/mol',
            ),
            Arrhenius(
                A = (1.56166e+25, 's^-1'),
                n = -3.82752,
                Ea = (200.45, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 18.2932, dn = +|- 0.393312, dEa = +|- 1.78921 kJ/mol',
            ),
            Arrhenius(
                A = (9.0443e+23, 's^-1'),
                n = -3.40545,
                Ea = (199.312, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 18.5692, dn = +|- 0.395338, dEa = +|- 1.79843 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   1.794e+07   1.785e+07   1.772e+07   1.756e+07   1.734e+07   1.707e+07   1.673e+07
#           800   6.850e+04   6.691e+04   6.501e+04   6.282e+04   6.038e+04   5.777e+04   5.506e+04
#           600   3.839e+02   3.673e+02   3.496e+02   3.317e+02   3.141e+02   2.972e+02   2.812e+02
#           480   2.741e+00   2.570e+00   2.408e+00   2.260e+00   2.125e+00   2.003e+00   1.887e+00
#           400   2.133e-02   1.988e-02   1.859e-02   1.745e-02   1.642e-02   1.544e-02   1.445e-02
#       342.857   1.766e-04   1.652e-04   1.551e-04   1.460e-04   1.373e-04   1.285e-04   1.189e-04
#           300   1.534e-06   1.438e-06   1.352e-06   1.272e-06   1.192e-06   1.107e-06   1.012e-06
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['CH2O', 'CH3'],
    products = ['[CH2]CO'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (54.6296, 'cm^3/(mol*s)'),
                n = 3.00648,
                Ea = (86.0924, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.24801, dn = +|- 0.109615, dEa = +|- 0.498651 kJ/mol',
            ),
            Arrhenius(
                A = (24.7546, 'cm^3/(mol*s)'),
                n = 3.11413,
                Ea = (85.8197, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 1.87706, dn = +|- 0.0852117, dEa = +|- 0.387636 kJ/mol',
            ),
            Arrhenius(
                A = (10.2647, 'cm^3/(mol*s)'),
                n = 3.23276,
                Ea = (85.4689, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 1.70731, dn = +|- 0.0723851, dEa = +|- 0.329287 kJ/mol',
            ),
            Arrhenius(
                A = (4.29156, 'cm^3/(mol*s)'),
                n = 3.34944,
                Ea = (85.1045, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 1.70806, dn = +|- 0.0724445, dEa = +|- 0.329557 kJ/mol',
            ),
            Arrhenius(
                A = (2.00379, 'cm^3/(mol*s)'),
                n = 3.45072,
                Ea = (84.8007, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 1.83513, dn = +|- 0.0821553, dEa = +|- 0.373732 kJ/mol',
            ),
            Arrhenius(
                A = (1.1439, 'cm^3/(mol*s)'),
                n = 3.52478,
                Ea = (84.6307, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.06397, dn = +|- 0.0980571, dEa = +|- 0.446071 kJ/mol',
            ),
            Arrhenius(
                A = (0.852027, 'cm^3/(mol*s)'),
                n = 3.5633,
                Ea = (84.6558, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 2.38512, dn = +|- 0.117627, dEa = +|- 0.535097 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   1.426e+05   1.852e+05   2.373e+05   3.001e+05   3.740e+05   4.588e+05   5.535e+05
#           800   7.494e+02   8.501e+02   9.490e+02   1.042e+03   1.126e+03   1.196e+03   1.248e+03
#           600   3.669e+00   3.805e+00   3.889e+00   3.914e+00   3.876e+00   3.773e+00   3.603e+00
#           480   1.857e-02   1.819e-02   1.754e-02   1.662e-02   1.545e-02   1.404e-02   1.246e-02
#           400   1.106e-04   1.037e-04   9.530e-05   8.561e-05   7.502e-05   6.403e-05   5.317e-05
#       342.857   7.267e-07   6.537e-07   5.734e-07   4.896e-07   4.066e-07   3.284e-07   2.581e-07
#           300   5.387e-09   4.650e-09   3.902e-09   3.183e-09   2.524e-09   1.949e-09   1.469e-09
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['CC[O]'],
    products = ['[H]', 'C=CO'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (76492.4, 's^-1'),
                n = 1.46228,
                Ea = (96.5578, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 24.0114, dn = +|- 0.430119, dEa = +|- 1.95665 kJ/mol',
            ),
            Arrhenius(
                A = (7542.57, 's^-1'),
                n = 1.82084,
                Ea = (96.2334, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 20.4088, dn = +|- 0.408121, dEa = +|- 1.85658 kJ/mol',
            ),
            Arrhenius(
                A = (823.582, 's^-1'),
                n = 2.16473,
                Ea = (96.0245, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 16.7847, dn = +|- 0.381666, dEa = +|- 1.73624 kJ/mol',
            ),
            Arrhenius(
                A = (97.6593, 's^-1'),
                n = 2.4965,
                Ea = (95.9218, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 13.6863, dn = +|- 0.354051, dEa = +|- 1.61061 kJ/mol',
            ),
            Arrhenius(
                A = (12.0918, 's^-1'),
                n = 2.82128,
                Ea = (95.901, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 11.3793, dn = +|- 0.329071, dEa = +|- 1.49697 kJ/mol',
            ),
            Arrhenius(
                A = (1.47864, 's^-1'),
                n = 3.14633,
                Ea = (95.9243, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 9.9369, dn = +|- 0.31073, dEa = +|- 1.41354 kJ/mol',
            ),
            Arrhenius(
                A = (0.167061, 's^-1'),
                n = 3.48028,
                Ea = (95.9439, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 9.36397, dn = +|- 0.302694, dEa = +|- 1.37698 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   2.541e+06   3.000e+06   3.504e+06   4.047e+06   4.623e+06   5.220e+06   5.826e+06
#           800   8.525e+03   9.038e+03   9.501e+03   9.910e+03   1.026e+04   1.056e+04   1.080e+04
#           600   1.317e+01   1.338e+01   1.354e+01   1.367e+01   1.376e+01   1.383e+01   1.389e+01
#           480   1.631e-02   1.636e-02   1.641e-02   1.644e-02   1.646e-02   1.647e-02   1.648e-02
#           400   2.284e-05   2.286e-05   2.287e-05   2.288e-05   2.289e-05   2.290e-05   2.290e-05
#       342.857   3.358e-08   3.359e-08   3.360e-08   3.360e-08   3.361e-08   3.361e-08   3.361e-08
#           300   5.641e-11   5.642e-11   5.643e-11   5.643e-11   5.643e-11   5.643e-11   5.643e-11
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['C[CH]O'],
    products = ['[H]', 'C=CO'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (6.55338e+18, 's^-1'),
                n = -2.05968,
                Ea = (138.055, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 637.932, dn = +|- 0.873929, dEa = +|- 3.97559 kJ/mol',
            ),
            Arrhenius(
                A = (4.50447e+17, 's^-1'),
                n = -1.67723,
                Ea = (136.797, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 480.083, dn = +|- 0.835461, dEa = +|- 3.80059 kJ/mol',
            ),
            Arrhenius(
                A = (3.45053e+16, 's^-1'),
                n = -1.31118,
                Ea = (135.576, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 350.654, dn = +|- 0.792949, dEa = +|- 3.6072 kJ/mol',
            ),
            Arrhenius(
                A = (3.00776e+15, 's^-1'),
                n = -0.964306,
                Ea = (134.403, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 250.515, dn = +|- 0.747444, dEa = +|- 3.40019 kJ/mol',
            ),
            Arrhenius(
                A = (3.03915e+14, 's^-1'),
                n = -0.639066,
                Ea = (133.29, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 176.578, dn = +|- 0.700114, dEa = +|- 3.18489 kJ/mol',
            ),
            Arrhenius(
                A = (3.61613e+13, 's^-1'),
                n = -0.337533,
                Ea = (132.249, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 123.924, dn = +|- 0.652199, dEa = +|- 2.96691 kJ/mol',
            ),
            Arrhenius(
                A = (5.13056e+12, 's^-1'),
                n = -0.06134,
                Ea = (131.286, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 87.3943, dn = +|- 0.604939, dEa = +|- 2.75193 kJ/mol',
            ),
        ],
    ),
)

#   =========== =========== =========== =========== =========== =========== =========== =========== 
#         T / P   1.013e+00   1.487e+00   2.183e+00   3.204e+00   4.703e+00   6.903e+00   1.013e+01
#   =========== =========== =========== =========== =========== =========== =========== =========== 
#          1200   1.693e+08   1.599e+08   1.496e+08   1.383e+08   1.262e+08   1.135e+08   1.005e+08
#           800   3.837e+06   3.206e+06   2.630e+06   2.118e+06   1.673e+06   1.295e+06   9.825e+05
#           600   8.450e+04   6.336e+04   4.666e+04   3.376e+04   2.399e+04   1.673e+04   1.146e+04
#           480   2.279e+03   1.586e+03   1.086e+03   7.320e+02   4.848e+02   3.149e+02   2.003e+02
#           400   8.259e+01   5.485e+01   3.582e+01   2.294e+01   1.438e+01   8.805e+00   5.251e+00
#       342.857   3.797e+00   2.427e+00   1.518e+00   9.272e-01   5.517e-01   3.193e-01   1.795e-01
#           300   2.164e-01   1.331e-01   7.973e-02   4.648e-02   2.633e-02   1.449e-02   7.758e-03
#   =========== =========== =========== =========== =========== =========== =========== =========== 
pdepreaction(
    reactants = ['CH2O', 'CH3'],
    products = ['[H]', 'C=CO'],
    kinetics = PDepArrhenius(
        pressures = ([1.01325, 1.48725, 2.18298, 3.20418, 4.70309, 6.90319, 10.1325], 'bar'),
        arrhenius = [
            Arrhenius(
                A = (0.047569, 'cm^3/(mol*s)'),
                n = 3.84394,
                Ea = (51.2806, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 282.969, dn = +|- 0.763928, dEa = +|- 3.47518 kJ/mol',
            ),
            Arrhenius(
                A = (0.00286999, 'cm^3/(mol*s)'),
                n = 4.22803,
                Ea = (50.9548, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 330.929, dn = +|- 0.785114, dEa = +|- 3.57156 kJ/mol',
            ),
            Arrhenius(
                A = (0.000171349, 'cm^3/(mol*s)'),
                n = 4.61236,
                Ea = (50.6648, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 353.278, dn = +|- 0.793958, dEa = +|- 3.61179 kJ/mol',
            ),
            Arrhenius(
                A = (1.04471e-05, 'cm^3/(mol*s)'),
                n = 4.99267,
                Ea = (50.435, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 347.919, dn = +|- 0.79189, dEa = +|- 3.60238 kJ/mol',
            ),
            Arrhenius(
                A = (6.58614e-07, 'cm^3/(mol*s)'),
                n = 5.3672,
                Ea = (50.2763, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 321.94, dn = +|- 0.781388, dEa = +|- 3.55461 kJ/mol',
            ),
            Arrhenius(
                A = (4.25427e-08, 'cm^3/(mol*s)'),
                n = 5.73705,
                Ea = (50.1845, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 286.998, dn = +|- 0.765841, dEa = +|- 3.48388 kJ/mol',
            ),
            Arrhenius(
                A = (2.73084e-09, 'cm^3/(mol*s)'),
                n = 6.10611,
                Ea = (50.1387, 'kJ/mol'),
                T0 = (1, 'K'),
                Tmin = (300, 'K'),
                Tmax = (1200, 'K'),
                comment = 'Fitted to 7 data points; dA = *|/ 253.98, dn = +|- 0.749302, dEa = +|- 3.40865 kJ/mol',
            ),
        ],
    ),
)

