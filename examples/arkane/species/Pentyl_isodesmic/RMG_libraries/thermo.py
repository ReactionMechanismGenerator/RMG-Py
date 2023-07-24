#!/usr/bin/env python
# encoding: utf-8

name = "thermo"
shortDesc = ""
longDesc = """
Calculated using Arkane v3.1.0 using LevelOfTheory(method='wb97mv',basis='def2tzvpd',software='qchem').
"""
entry(
    index = 0,
    label = "Pentyl",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {3,S} {15,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.51823,0.0511227,-0.000134204,3.60074e-07,-3.16044e-10,4736.33,11.2552], Tmin=(10,'K'), Tmax=(433.56,'K')),
            NASAPolynomial(coeffs=[-1.27148,0.0585866,-3.29658e-05,9.02833e-09,-9.65105e-13,5496.83,34.3407], Tmin=(433.56,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (39.3687,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (382.466,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Spin multiplicity: 2
External symmetry: -1.0
Optical isomers: 1

Geometry:
C       2.45139259   -0.32620661   -0.00254336
C       1.19045885    0.52326516    0.00668734
C      -0.06949942   -0.34262836   -0.00966205
C      -1.33201159    0.51731764    0.00366158
C      -2.56330119   -0.32295505    0.00289586
H       3.33934460    0.31354145    0.01128396
H       2.49682513   -0.95126819   -0.90003707
H       2.49035800   -0.98020259    0.87447088
H       1.19760942    1.18716466   -0.86554746
H       1.19297733    1.15975846    0.89926367
H      -0.06570761   -1.01096590    0.86088115
H      -0.06426154   -0.98094919   -0.90245127
H      -1.35026429    1.17589679   -0.87209010
H      -1.34288089    1.16065340    0.89068816
H      -3.01248882   -0.63245458   -0.93359075
H      -2.94374681   -0.73393098    0.93089260
""",
)

