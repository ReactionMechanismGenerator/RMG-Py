#!/usr/bin/env python
# encoding: utf-8

name = "ethane-oxidation"
shortDesc = u""
longDesc = u"""

"""
duplicatesChecked=False
entry(
    index = 0,
    label = "[CH3] + [CH3] <=> ethane",
    degeneracy = 1.0,
    kinetics = Chebyshev(
        coeffs = [
            [12.8791, 0.560009, -0.115177, 0.0105093],
            [-1.04561, 0.857712, -0.116582, -0.0134747],
            [-0.599769, 0.449793, 0.0096248, -0.0271013],
            [-0.314421, 0.177204, 0.0458462, -0.0107188],
            [-0.155503, 0.0517415, 0.0323201, 0.00343069],
            [-0.0727158, 0.00622718, 0.0138549, 0.00630263],
        ],
        kunits = 'cm^3/(mol*s)',
        Tmin = (300, 'K'),
        Tmax = (3000, 'K'),
        Pmin = (0.001, 'bar'),
        Pmax = (100, 'bar'),
    ),
)

entry(
    index = 1,
    label = "O2 + [CH3] <=> CO[O]",
    degeneracy = 1.0,
    kinetics = Chebyshev(
        coeffs = [
            [9.90191, 1.69355, -0.276112, -0.014107],
            [-0.201801, 0.693442, 0.183308, -0.0255194],
            [-0.446451, 0.0857639, 0.0676062, 0.0247624],
            [-0.137012, -0.010582, 0.00986909, 0.00873805],
            [-0.0347555, 0.00594742, -0.00642878, -0.00257428],
            [0.0288809, -0.0156872, -0.00508067, 0.000142012],
        ],
        kunits = 'cm^3/(mol*s)',
        Tmin = (300, 'K'),
        Tmax = (3000, 'K'),
        Pmin = (0.001, 'bar'),
        Pmax = (100, 'bar'),
    ),
)

entry(
    index = 2,
    label = "ethane + [CH3] <=> C[CH2] + CH4",
    degeneracy = 6.0,
    kinetics = Arrhenius(
        A = (2.94005e-11, 'm^3/(mol*s)'),
        n = 5.135,
        Ea = (33.0118, 'kJ/mol'),
        T0 = (1, 'K'),
        comment = 'Estimated using average of templates [C/H3/Cs;C_methyl] + [C/H3/Cs\\H3;Cs_rad] for rate rule [C/H3/Cs\\H3;C_methyl]\nEuclidian distance = 1.0\nMultiplied by reaction path degeneracy 6.0\nfamily: H_Abstraction',
    ),
    longDesc = 
u"""
Estimated using average of templates [C/H3/Cs;C_methyl] + [C/H3/Cs\H3;Cs_rad] for rate rule [C/H3/Cs\H3;C_methyl]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: H_Abstraction
""",
)

entry(
    index = 3,
    label = "[CH3] + C[CH2] <=> C=C + CH4",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (6.57e+14, 'cm^3/(mol*s)', '*|/', 1.1),
        n = -0.68,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
        comment = 'Exact match found for rate rule [C_methyl;Cmethyl_Csrad]\nEuclidian distance = 0\nMultiplied by reaction path degeneracy 3.0\nfamily: Disproportionation',
    ),
    longDesc = 
u"""
Exact match found for rate rule [C_methyl;Cmethyl_Csrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: Disproportionation
""",
)

entry(
    index = 4,
    label = "C[CH2] + C[CH2] <=> ethane + C=C",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (6.9e+13, 'cm^3/(mol*s)', '*|/', 1.1),
        n = -0.35,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
        comment = 'Exact match found for rate rule [C_rad/H2/Cs;Cmethyl_Csrad]\nEuclidian distance = 0\nMultiplied by reaction path degeneracy 3.0\nfamily: Disproportionation',
    ),
    longDesc = 
u"""
Exact match found for rate rule [C_rad/H2/Cs;Cmethyl_Csrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: Disproportionation
""",
)

entry(
    index = 5,
    label = "C[CH2] + H <=> ethane",
    degeneracy = 1.0,
    kinetics = Chebyshev(
        coeffs = [
            [13.5315, 0.449416, -0.0774212, 0.00616325],
            [-0.82609, 0.750471, -0.0982394, -0.00339778],
            [-0.593837, 0.463725, -0.0136024, -0.0157443],
            [-0.368451, 0.22034, 0.0283534, -0.0104839],
            [-0.206783, 0.0792869, 0.0290969, -0.000446421],
            [-0.108198, 0.0181184, 0.0155923, 0.00424989],
        ],
        kunits = 'cm^3/(mol*s)',
        Tmin = (300, 'K'),
        Tmax = (3000, 'K'),
        Pmin = (0.001, 'bar'),
        Pmax = (100, 'bar'),
    ),
)

entry(
    index = 6,
    label = "C=C + H <=> C[CH2]",
    degeneracy = 1.0,
    kinetics = Chebyshev(
        coeffs = [
            [12.0332, 0.873833, -0.162407, 0.0168639],
            [-0.164957, 1.00913, -0.0324404, -0.0377042],
            [-0.453805, 0.37093, 0.0671946, -0.0148132],
            [-0.204108, 0.0773128, 0.043048, 0.00631233],
            [-0.0559888, -0.0120163, 0.0102436, 0.00671826],
            [0.00490435, -0.021391, -0.00366922, 0.00178308],
        ],
        kunits = 'cm^3/(mol*s)',
        Tmin = (300, 'K'),
        Tmax = (3000, 'K'),
        Pmin = (0.001, 'bar'),
        Pmax = (100, 'bar'),
    ),
)

entry(
    index = 7,
    label = "[CH3] + H <=> CH4",
    degeneracy = 1.0,
    kinetics = Chebyshev(
        coeffs = [
            [12.8642, 1.35989, -0.256864, -0.00362492],
            [-0.9798, 0.658274, 0.0852943, -0.0312454],
            [-0.472679, 0.155862, 0.0540246, 0.00851828],
            [-0.223854, 0.0394752, 0.0179535, 0.00554265],
            [-0.106076, 0.00938481, 0.00466234, 0.00169057],
            [-0.0493405, 0.00170846, 0.000954967, 0.000351117],
        ],
        kunits = 'cm^3/(mol*s)',
        Tmin = (300, 'K'),
        Tmax = (3000, 'K'),
        Pmin = (0.001, 'bar'),
        Pmax = (100, 'bar'),
    ),
)

entry(
    index = 8,
    label = "ethane + H <=> H2 + C[CH2]",
    degeneracy = 6.0,
    kinetics = Arrhenius(
        A = (2.17494e-06, 'm^3/(mol*s)'),
        n = 4.07,
        Ea = (25.4387, 'kJ/mol'),
        T0 = (1, 'K'),
        comment = 'Estimated using average of templates [C/H3/Cs;H_rad] + [C/H3/Cs\\H3;Y_rad] for rate rule [C/H3/Cs\\H3;H_rad]\nEuclidian distance = 1.0\nMultiplied by reaction path degeneracy 6.0\nfamily: H_Abstraction',
    ),
    longDesc = 
u"""
Estimated using average of templates [C/H3/Cs;H_rad] + [C/H3/Cs\H3;Y_rad] for rate rule [C/H3/Cs\H3;H_rad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: H_Abstraction
""",
)

entry(
    index = 9,
    label = "C[CH2] + H <=> C=C + H2",
    degeneracy = 3.0,
    kinetics = Arrhenius(
        A = (1.083e+13, 'cm^3/(mol*s)', '*|/', 2),
        n = 0,
        Ea = (0, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (2500, 'K'),
        comment = 'Exact match found for rate rule [H_rad;Cmethyl_Csrad]\nEuclidian distance = 0\nMultiplied by reaction path degeneracy 3.0\nfamily: Disproportionation',
    ),
    longDesc = 
u"""
Exact match found for rate rule [H_rad;Cmethyl_Csrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: Disproportionation
""",
)

entry(
    index = 10,
    label = "H + CH4 <=> H2 + [CH3]",
    degeneracy = 4.0,
    kinetics = Arrhenius(
        A = (0.876, 'cm^3/(mol*s)'),
        n = 4.34,
        Ea = (34.3088, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (1500, 'K'),
        comment = 'Exact match found for rate rule [C_methane;H_rad]\nEuclidian distance = 0\nMultiplied by reaction path degeneracy 4.0\nfamily: H_Abstraction',
    ),
    longDesc = 
u"""
Exact match found for rate rule [C_methane;H_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: H_Abstraction
""",
)

entry(
    index = 11,
    label = "H + H <=> H2",
    degeneracy = 1.0,
    kinetics = Chebyshev(
        coeffs = [
            [8.48588, 1.72816, -0.232832, -0.0184726],
            [0.0439119, 0.262399, 0.0740112, -0.0228139],
            [-0.107787, 0.0499097, 0.00176484, -0.00151681],
            [-0.0405237, 0.0129564, 0.00214129, 0.00242554],
            [-0.0159905, 0.00413936, 0.0016664, 0.00147354],
            [-0.00614586, 0.00151107, 0.000719274, 0.000598049],
        ],
        kunits = 'cm^3/(mol*s)',
        Tmin = (300, 'K'),
        Tmax = (3000, 'K'),
        Pmin = (0.001, 'bar'),
        Pmax = (100, 'bar'),
    ),
)

entry(
    index = 12,
    label = "C[CH2] + [O]O <=> ethane + O2",
    degeneracy = 1.0,
    kinetics = Arrhenius(
        A = (7.74112e-12, 'm^3/(mol*s)'),
        n = 4.92,
        Ea = (19.623, 'kJ/mol'),
        T0 = (1, 'K'),
        comment = 'Estimated using template [X_H;C_rad/H2/Cs\\H3] for rate rule [Orad_O_H;C_rad/H2/Cs\\H3]\nEuclidian distance = 1.0\nfamily: H_Abstraction',
    ),
    longDesc = 
u"""
Estimated using template [X_H;C_rad/H2/Cs\H3] for rate rule [Orad_O_H;C_rad/H2/Cs\H3]
Euclidian distance = 1.0
family: H_Abstraction
""",
)

entry(
    index = 13,
    label = "O2 + C[CH2] <=> C=C + [O]O",
    degeneracy = 6.0,
    kinetics = Arrhenius(
        A = (4.338e+13, 'cm^3/(mol*s)', '*|/', 3),
        n = 0,
        Ea = (66.9022, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (700, 'K'),
        Tmax = (2500, 'K'),
        comment = 'Exact match found for rate rule [O2b;Cmethyl_Csrad]\nEuclidian distance = 0\nMultiplied by reaction path degeneracy 6.0\nfamily: Disproportionation',
    ),
    longDesc = 
u"""
Exact match found for rate rule [O2b;Cmethyl_Csrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: Disproportionation
""",
)

entry(
    index = 14,
    label = "O2 + CH4 <=> [CH3] + [O]O",
    degeneracy = 8.0,
    kinetics = Arrhenius(
        A = (7.94e+13, 'cm^3/(mol*s)', '*|/', 10),
        n = 0,
        Ea = (237.777, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (500, 'K'),
        Tmax = (2000, 'K'),
        comment = 'Exact match found for rate rule [C_methane;O2b]\nEuclidian distance = 0\nMultiplied by reaction path degeneracy 8.0\nfamily: H_Abstraction',
    ),
    longDesc = 
u"""
Exact match found for rate rule [C_methane;O2b]
Euclidian distance = 0
Multiplied by reaction path degeneracy 8.0
family: H_Abstraction
""",
)

entry(
    index = 15,
    label = "O2 + C[CH2] <=> C=C + [O]O",
    degeneracy = 1.0,
    kinetics = Chebyshev(
        coeffs = [
            [9.44348, -0.969233, -0.242864, -0.000655854],
            [1.68707, 0.910295, 0.146246, -0.0476217],
            [0.118536, 0.124642, 0.0773244, 0.0178854],
            [-0.0263916, -0.0713451, 0.00709469, 0.014242],
            [-0.0315143, -0.0469276, -0.0141423, 0.00204769],
            [-0.0174421, -0.00689436, -0.00820692, -0.00214934],
        ],
        kunits = 'cm^3/(mol*s)',
        Tmin = (300, 'K'),
        Tmax = (3000, 'K'),
        Pmin = (0.001, 'bar'),
        Pmax = (100, 'bar'),
    ),
)

entry(
    index = 16,
    label = "O2 + H <=> [O]O",
    degeneracy = 1.0,
    kinetics = Chebyshev(
        coeffs = [
            [10.9076, 1.98606, -0.233099, -0.0481274],
            [-0.571345, 0.31254, 0.126352, 0.0124235],
            [-0.307278, 0.0783359, 0.0239632, 0.00365003],
            [-0.0857745, -0.0227595, 0.00823973, 0.00871158],
            [-0.0494686, 0.00279831, 0.00361649, 0.000666411],
            [-0.0263897, 0.0102559, -0.00268875, -0.00159869],
        ],
        kunits = 'cm^3/(mol*s)',
        Tmin = (300, 'K'),
        Tmax = (3000, 'K'),
        Pmin = (0.001, 'bar'),
        Pmax = (100, 'bar'),
    ),
)

entry(
    index = 17,
    label = "O2 + H2 <=> H + [O]O",
    degeneracy = 4.0,
    kinetics = Arrhenius(
        A = (2.9e+14, 'cm^3/(mol*s)', '*|/', 5),
        n = 0,
        Ea = (236.982, 'kJ/mol'),
        T0 = (1, 'K'),
        Tmin = (300, 'K'),
        Tmax = (800, 'K'),
        comment = 'Exact match found for rate rule [H2;O2b]\nEuclidian distance = 0\nMultiplied by reaction path degeneracy 4.0\nfamily: H_Abstraction',
    ),
    longDesc = 
u"""
Exact match found for rate rule [H2;O2b]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: H_Abstraction
""",
)

