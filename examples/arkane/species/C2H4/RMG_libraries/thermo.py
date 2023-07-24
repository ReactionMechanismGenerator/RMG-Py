#!/usr/bin/env python
# encoding: utf-8

name = "thermo"
shortDesc = ""
longDesc = """
Calculated using Arkane v3.1.0 using LevelOfTheory(method='cbsqb3').
"""
entry(
    index = 0,
    label = "C2H4",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.10211,-0.0068058,4.97293e-05,-5.81716e-08,2.24308e-11,5504.57,3.22404], Tmin=(10,'K'), Tmax=(775.817,'K')),
            NASAPolynomial(coeffs=[0.436833,0.0179601,-9.50042e-06,2.47486e-09,-2.53793e-13,5896.68,18.8383], Tmin=(775.817,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (45.7904,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (133.032,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Spin multiplicity: 1
External symmetry: -1.0
Optical isomers: 1

Geometry:
C       0.00172300    0.00000000    0.00107000
H      -0.00018100    0.00000000    1.08623200
H       0.97503900    0.00000000   -0.47876200
C      -1.12459600    0.00000000   -0.70077800
H      -1.12099200    0.00000000   -1.78578100
H      -2.09702200    0.00000000   -0.21949100
""",
)

