#!/usr/bin/env python
# encoding: utf-8

name = "/home/alongd/Code/RMG-Py/importer/pentaneNO"
shortDesc = u"/home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat"
longDesc = u"""
Unknown source
"""
entry(
    index = 1,
    label = "H2O2",
    molecule = 
"""
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.31515,-0.000847391,1.76404e-05,-2.26763e-08,9.0895e-12,-17706.7,3.27373], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.57977,0.00405326,-1.29845e-06,1.98211e-10,-1.13969e-14,-18007.2,0.664971], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 8/03""",
    longDesc = 
u"""
T 8/03.
OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 2,
    label = "OH",
    molecule = 
"""
multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.99198,-0.00240107,4.61664e-06,-3.87916e-09,1.3632e-12,3368.9,-0.103998], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.83853,0.00110741,-2.94e-07,4.20699e-11,-2.4229e-15,3697.81,5.84495], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""IU3/03""",
    longDesc = 
u"""
IU3/03.
[OH]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 3,
    label = "H2",
    molecule = 
"""
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.34433,0.00798052,-1.94782e-05,2.01572e-08,-7.37612e-12,-917.935,0.68301], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.93287,0.000826608,-1.46402e-07,1.541e-11,-6.88805e-16,-813.066,-1.02433], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""TPIS78""",
    longDesc = 
u"""
TPIS78.
[H][H]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 4,
    label = "H",
    molecule = 
"""
multiplicity 2
1 H u1 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5,0,0,0,0,25473.7,-0.446683], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.5,0,0,0,0,25473.7,-0.446683], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 6/94""",
    longDesc = 
u"""
L 6/94.
[H]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 5,
    label = "O",
    molecule = 
"""
multiplicity 3
1 O u2 p2 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.16827,-0.00327932,6.64306e-06,-6.12807e-09,2.11266e-12,29122.3,2.05193], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.54364,-2.73162e-05,-4.1903e-09,4.95482e-12,-4.79554e-16,29226,4.92229], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 1/90""",
    longDesc = 
u"""
L 1/90.
[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 6,
    label = "H2O",
    molecule = 
"""
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.19864,-0.0020364,6.52034e-06,-5.48793e-09,1.77197e-12,-30293.7,-0.849009], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.67704,0.00297318,-7.73769e-07,9.44335e-11,-4.269e-15,-29885.9,6.88255], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 5/89""",
    longDesc = 
u"""
L 5/89.
O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 7,
    label = "HO2",
    molecule = 
"""
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.3018,-0.00474912,2.11583e-05,-2.42764e-08,9.29225e-12,264.018,3.71666], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.17229,0.00188118,-3.46277e-07,1.94658e-11,1.76257e-16,31.0207,2.95768], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""T 1/09""",
    longDesc = 
u"""
T 1/09.
[O]O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 8,
    label = "O2",
    molecule = 
"""
multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.78246,-0.00299673,9.8473e-06,-9.6813e-09,3.24373e-12,-1063.94,3.65768], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.66096,0.000656366,-1.4115e-07,2.05798e-11,-1.29913e-15,-1215.98,3.41536], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""RUS 89""",
    longDesc = 
u"""
RUS 89.
[O][O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 9,
    label = "AR",
    molecule = 
"""
1 Ar u0 p4 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""G 5/97""",
    longDesc = 
u"""
G 5/97.
[Ar]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 10,
    label = "N2",
    molecule = 
"""
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""G 8/02""",
    longDesc = 
u"""
G 8/02.
N#N
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 11,
    label = "HE",
    molecule = 
"""
1 He u0 p1 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""G 5/97""",
    longDesc = 
u"""
G 5/97.
[He]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 12,
    label = "OHV",
    molecule = 
"""
multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.63727,0.000185091,-1.67616e-06,2.3872e-09,-8.43144e-13,50021.3,1.35886], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.88273,0.00101397,-2.27688e-07,2.17468e-11,-5.12631e-16,50265,5.59571], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""121286""",
    longDesc = 
u"""
121286
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
Duplicate of species OH (i.e. same molecular structure according to RMG)
[OH]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 13,
    label = "CO",
    molecule = 
"""
multiplicity 3
1 C u2 p0 c0 {2,D}
2 O u0 p2 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.57953,-0.000610354,1.01681e-06,9.07006e-10,-9.04424e-13,-14344.1,3.50841], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.04849,0.00135173,-4.85794e-07,7.88536e-11,-4.69807e-15,-14266.1,6.0171], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""RUS 79""",
    longDesc = 
u"""
RUS 79.
[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 14,
    label = "CO2",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,D}
2 O u0 p2 c0 {1,D}
3 O u0 p2 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.35681,0.00898413,-7.12206e-06,2.4573e-09,-1.42885e-13,-48372,9.9009], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.63651,0.00274146,-9.95898e-07,1.60387e-10,-9.16199e-15,-49024.9,-1.9349], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 7/88""",
    longDesc = 
u"""
L 7/88.
O=C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 15,
    label = "CH4",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.14911,-0.0136622,4.91454e-05,-4.84247e-08,1.66603e-11,-10246.6,-4.63849], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[1.65326,0.0100263,-3.31661e-06,5.36483e-10,-3.14697e-14,-10009.6,9.90506], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""G 8/99""",
    longDesc = 
u"""
G 8/99.
C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 16,
    label = "CH3",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""IU0702""",
    longDesc = 
u"""
IU0702.
[CH3]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 17,
    label = "C",
    molecule = 
"""
1 C u0 p2 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.55424,-0.000321538,7.33792e-07,-7.32235e-10,2.66521e-13,85442.7,4.53131], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.60558,-0.000195934,1.06737e-07,-1.64239e-11,8.18706e-16,85411.7,4.19239], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 7/88""",
    longDesc = 
u"""
L 7/88.
[C]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 18,
    label = "CH",
    molecule = 
"""
multiplicity 4
1 C u3 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.48976,0.000324322,-1.68998e-06,3.16284e-09,-1.40618e-12,70612.6,2.08428], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.52094,0.00176536,-4.61477e-07,5.92897e-11,-3.34745e-15,70946.8,7.40518], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""IU3/03""",
    longDesc = 
u"""
IU3/03.
[CH]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 19,
    label = "CHV",
    molecule = 
"""
multiplicity 4
1 C u3 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.2002,0.00207288,-5.13443e-06,5.73389e-09,-1.95553e-12,103937,3.33159], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.19622,0.00234038,-7.0582e-07,9.00758e-11,-3.85504e-15,104196,9.17837], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""073003""",
    longDesc = 
u"""
073003
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
Duplicate of species CH (i.e. same molecular structure according to RMG)
[CH]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 20,
    label = "CH3OH",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.65851,-0.0162983,6.91938e-05,-7.58373e-08,2.80428e-11,-25612,-0.897331], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.52727,0.0103179,-3.62893e-06,5.77448e-10,-3.42183e-14,-26002.9,5.16759], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T06/02""",
    longDesc = 
u"""
T06/02.
CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 21,
    label = "C2H6",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.29143,-0.00550155,5.99438e-05,-7.08466e-08,2.68686e-11,-11522.2,2.66679], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.04666,0.0153539,-5.47039e-06,8.77827e-10,-5.23168e-14,-12447.3,-0.968698], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""G 8/88""",
    longDesc = 
u"""
G 8/88.
CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 22,
    label = "C2H5",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u1 p0 c0 {1,S} {6,S} {7,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.3273,0.0176657,-6.14927e-06,-3.01143e-10,4.38618e-13,13428.4,17.1789], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[5.88784,0.0103077,-3.46844e-06,5.32499e-10,-3.06513e-14,11506.5,-8.49652], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/ 4/ 4 THERM""",
    longDesc = 
u"""
8/ 4/ 4 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 23,
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
            NASAPolynomial(coeffs=[0.481118,0.0183778,-9.99634e-06,2.73211e-09,-3.01837e-13,5443.87,18.5867], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[5.07061,0.00911141,-3.10507e-06,4.80734e-10,-2.78321e-14,3663.91,-6.64501], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 24,
    label = "C2H",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,T} {3,S}
2 C u1 p0 c0 {1,T}
3 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.89868,0.0132988,-2.80733e-05,2.89485e-08,-1.07502e-11,67061.6,6.18548], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.6627,0.00382492,-1.36633e-06,2.13455e-10,-1.23217e-14,67168.4,3.92206], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 5/10""",
    longDesc = 
u"""
T 5/10.
[C]#C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 25,
    label = "C3H8",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.240878,0.0339549,-1.60931e-05,2.83481e-09,2.78195e-14,-14036.3,21.6501], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[9.15541,0.0172574,-5.85615e-06,9.0419e-10,-5.22524e-14,-17576.2,-27.7419], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 26,
    label = "C2H2",
    molecule = 
"""
1 C u0 p0 c0 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.80868,0.0233616,-3.55172e-05,2.80153e-08,-8.50075e-12,26429,13.9397], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.65878,0.00488397,-1.60829e-06,2.46975e-10,-1.38606e-14,25759.4,-3.99838], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""G 1/91""",
    longDesc = 
u"""
G 1/91.
C#C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 27,
    label = "CH3OCH3",
    molecule = 
"""
1 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3 O u0 p2 c0 {1,S} {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.05597,0.0207019,-5.00382e-06,-1.6228e-09,6.8433e-13,-23549.4,14.503], Tmin=(298,'K'), Tmax=(1999,'K')),
            NASAPolynomial(coeffs=[6.03233,0.0156155,-5.50761e-06,8.75666e-10,-5.17181e-14,-25269,-8.25885], Tmin=(1999,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/11/14 THERM""",
    longDesc = 
u"""
2/11/14 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
COC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 28,
    label = "O3",
    molecule = 
"""
1 O u0 p1 c+1 {2,S} {3,D}
2 O u0 p3 c-1 {1,S}
3 O u0 p2 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.40738,0.00205379,1.38486e-05,-2.23312e-08,9.76073e-12,15864.5,8.28248], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[12.3303,-0.0119325,7.98741e-06,-1.77195e-09,1.26076e-13,12675.6,-40.8823], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 5/90""",
    longDesc = 
u"""
L 5/90.
O=OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 29,
    label = "N",
    molecule = 
"""
multiplicity 4
1 N u3 p1 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5,0,0,0,0,56104.6,4.19391], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.41594,0.000174891,-1.19024e-07,3.02262e-11,-2.0361e-15,56133.8,4.64961], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 6/88""",
    longDesc = 
u"""
L 6/88.
[N]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 30,
    label = "NO2",
    molecule = 
"""
multiplicity 2
1 N u0 p1 c0 {2,D} {3,S}
2 O u0 p2 c0 {1,D}
3 O u1 p2 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.6706,0.0078385,-8.06386e-06,6.16171e-09,-2.32015e-12,2896.3,11.612], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.68286,0.00246243,-1.04226e-06,1.9769e-10,-1.392e-14,2261.3,0.989], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
N(=O)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 31,
    label = "CH3O",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 O u1 p2 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.71181,-0.00280463,3.76551e-05,-4.73072e-08,1.86588e-11,1295.7,6.57241], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.75779,0.00744142,-2.69705e-06,4.38091e-10,-2.63537e-14,378.112,-1.9668], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""IU1/03""",
    longDesc = 
u"""
IU1/03.
C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 32,
    label = "CH2OH",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 O u0 p2 c0 {1,S} {5,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.47834,-0.0013507,2.78485e-05,-3.64869e-08,1.47907e-11,-3500.73,3.30913], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.09314,0.00594761,-2.06497e-06,3.23008e-10,-1.88126e-14,-4034.1,-1.84691], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""IU2/03""",
    longDesc = 
u"""
IU2/03.
[CH2]O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 33,
    label = "CH2O",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 O u0 p2 c0 {1,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14379.2,0.602798], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.16953,0.00619321,-2.25056e-06,3.65976e-10,-2.20149e-14,-14548.7,6.04208], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 5/11""",
    longDesc = 
u"""
T 5/11.
C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 34,
    label = "NC3H7",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u1 p0 c0 {1,S} {9,S} {10,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.20121,0.0529642,-7.23641e-05,6.36997e-08,-2.29333e-11,11513.1,34.3669], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.48614,0.0165769,-5.74876e-06,9.04104e-10,-5.30867e-14,8937.1,-14.2595], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15.
[CH2]CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 35,
    label = "HCO",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,S} {3,D}
2 H u0 p0 c0 {1,S}
3 O u0 p2 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.23755,-0.00332075,1.4003e-05,-1.3424e-08,4.37416e-12,3872.41,3.30835], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.92002,0.00252279,-6.71004e-07,1.05616e-10,-7.43798e-15,3653.43,3.58077], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 5/03""",
    longDesc = 
u"""
T 5/03.
[CH]=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 36,
    label = "HONO",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {3,D}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p2 c0 {1,D}
4 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.29041,0.0140992,-1.36787e-05,7.49878e-09,-1.87691e-12,-10431.9,13.281], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.48689,0.00421806,-1.64914e-06,2.9719e-10,-2.021e-14,-11268.6,-2.997], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
N(=O)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 37,
    label = "C2H3",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u1 p0 c0 {1,D} {5,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.25545,0.0157482,-1.12218e-05,4.50916e-09,-7.74862e-13,34743.6,16.9664], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[4.99675,0.00655838,-2.20922e-06,3.393e-10,-1.95317e-14,33460.4,-3.01451], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 38,
    label = "CH2",
    molecule = 
"""
multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.71758,0.00127391,2.17347e-06,-3.48858e-09,1.65209e-12,45872.4,1.75298], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.14632,0.00303671,-9.96474e-07,1.50484e-10,-8.57336e-15,46041.3,4.72342], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""IU3/03""",
    longDesc = 
u"""
IU3/03.
[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 39,
    label = "CH3OCH2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2 C u1 p0 c0 {3,S} {7,S} {8,S}
3 O u0 p2 c0 {1,S} {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.58875,0.0224414,-1.19435e-05,3.3716e-09,-4.15077e-13,-1372.08,18.7549], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[6.62622,0.0122219,-4.12417e-06,6.34128e-10,-3.65317e-14,-3339.66,-8.95306], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/11/14 THERM""",
    longDesc = 
u"""
2/11/14 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]OC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 40,
    label = "C3H6",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,D} {7,S}
3 C u0 p0 c0 {2,D} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.54607,0.0436553,-5.61392e-05,4.98422e-08,-1.84799e-11,2070.56,29.9232], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.59032,0.0152593,-5.30369e-06,8.35511e-10,-4.91216e-14,-247.481,-11.5748], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15.
C=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 41,
    label = "IC3H7",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {2,S} {10,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.897467,0.0415744,-4.94778e-05,4.56494e-08,-1.79085e-11,9939.5,29.2642], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.70776,0.0174048,-6.07616e-06,9.60084e-10,-5.65656e-14,7553.78,-10.3687], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15.
C[CH]C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 42,
    label = "C3H5-A",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {4,S}
2 C u1 p0 c0 {1,S} {5,S} {6,S}
3 C u0 p0 c0 {1,D} {7,S} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.32899,0.0538423,-7.65501e-05,6.35512e-08,-2.14283e-11,20342.1,36.8038], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.37604,0.012345,-4.26464e-06,6.69046e-10,-3.92203e-14,17733.3,-16.1758], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15.
[CH2]C=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 43,
    label = "C3H5-T",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {3,D} {7,S} {8,S}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.29257,0.0198528,-6.42636e-06,-5.90016e-10,5.05491e-13,28577.3,13.9407], Tmin=(298,'K'), Tmax=(1376,'K')),
            NASAPolynomial(coeffs=[7.69949,0.0117804,-4.07792e-06,6.38119e-10,-3.7223e-14,26174.7,-16.8306], Tmin=(1376,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=[C]C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 44,
    label = "C3H5-S",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,D} {7,S}
3 C u1 p0 c0 {2,D} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.61793,0.0244804,-1.41857e-05,4.16402e-09,-4.90905e-13,30429.1,16.6341], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[7.95954,0.0111163,-3.75198e-06,5.77246e-10,-3.32769e-14,28056.8,-17.98], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 45,
    label = "C3H4-A",
    molecule = 
"""
1 C u0 p0 c0 {3,D} {4,S} {5,S}
2 C u0 p0 c0 {3,D} {6,S} {7,S}
3 C u0 p0 c0 {1,D} {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.61304,0.0121226,1.85399e-05,-3.45251e-08,1.53351e-11,21541.6,10.2261], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.31687,0.0111337,-3.96294e-06,6.35642e-10,-3.78755e-14,20117.5,-10.9958], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 8/89""",
    longDesc = 
u"""
L 8/89.
C=C=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 46,
    label = "C3H4-P",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.68039,0.0157997,2.50706e-06,-1.36576e-08,6.61543e-12,20802.4,9.87694], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.02524,0.0113365,-4.02234e-06,6.43761e-10,-3.82996e-14,19620.9,-8.60438], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 2/90""",
    longDesc = 
u"""
T 2/90.
C#CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 47,
    label = "NO",
    molecule = 
"""
multiplicity 2
1 N u1 p1 c0 {2,D}
2 O u0 p2 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.37654,0.00125306,-3.30275e-06,5.21781e-09,-2.44626e-12,9818,5.83], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.24544,0.00126914,-5.0159e-07,9.169e-11,-6.28e-15,9800.8,6.417], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""------------------------------------------------------------------------------------------------------""",
    longDesc = 
u"""
------------------------------------------------------------------------------------------------------
-------------------------------------- AJOUT DAGAUT --------------------------------------------------
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[N]=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 48,
    label = "HNO",
    molecule = 
"""
1 N u0 p1 c0 {2,D} {3,S}
2 O u0 p2 c0 {1,D}
3 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.7844,0.00660965,-9.30022e-06,9.43798e-09,-3.75315e-12,10918.8,9.036], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.61514,0.00321249,-1.26034e-06,2.2673e-10,-1.536e-14,10661.9,4.81], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
N=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 49,
    label = "NO3",
    molecule = 
"""
multiplicity 2
1 N u0 p0 c+1 {2,D} {3,S} {4,S}
2 O u0 p2 c0 {1,D}
3 O u0 p3 c-1 {1,S}
4 O u1 p2 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.17359,0.0104903,1.10473e-05,-2.81562e-08,1.36584e-11,7812.91,14.6022], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.48348,0.00257772,-1.00946e-06,1.72314e-10,-1.07154e-14,6129.9,-14.1618], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT/A""",
    longDesc = 
u"""
ATcT/A.
[N+](=O)([O-])[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 50,
    label = "C6H101-5",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,D} {11,S}
4  C u0 p0 c0 {2,S} {6,D} {12,S}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  C u0 p0 c0 {4,D} {15,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.01375,0.0638243,-4.40654e-05,1.58295e-08,-2.30831e-12,7940.34,32.5056], Tmin=(298,'K'), Tmax=(1413,'K')),
            NASAPolynomial(coeffs=[16.0456,0.0234774,-7.85798e-06,1.20201e-09,-6.901e-14,2118.99,-58.8452], Tmin=(1413,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""4/12/13 THERM""",
    longDesc = 
u"""
4/12/13 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCC=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 51,
    label = "C4H7O1-4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  O u1 p2 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.44895,0.0369424,-1.77511e-05,2.64538e-09,2.43666e-13,6163.86,17.3175], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[13.3251,0.0165058,-5.65236e-06,8.7832e-10,-5.09915e-14,1912.81,-42.8181], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 52,
    label = "HNO3",
    molecule = 
"""
1 N u0 p0 c+1 {2,S} {3,D} {4,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u0 p2 c0 {1,D}
4 O u0 p3 c-1 {1,S}
5 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.69329,0.0190168,-8.25177e-06,-6.06114e-09,4.65237e-12,-17419.9,17.184], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.03099,0.00446959,-1.72459e-06,2.91556e-10,-1.80103e-14,-19313.8,-16.2617], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 8/03""",
    longDesc = 
u"""
T 8/03.
[N+](=O)(O)[O-]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 53,
    label = "CH3O2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 O u1 p2 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.97339,0.0153542,-6.37315e-06,3.19931e-10,2.82194e-13,254.279,16.9194], Tmin=(298,'K'), Tmax=(1374,'K')),
            NASAPolynomial(coeffs=[6.4797,0.00744401,-2.52349e-06,3.89577e-10,-2.25182e-14,-1562.85,-8.19477], Tmin=(1374,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 54,
    label = "IC3H7O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  O u0 p2 c0 {1,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u1 p2 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.58518,0.0416107,-2.92194e-05,1.08615e-08,-1.66312e-12,-9670.13,14.4731], Tmin=(298,'K'), Tmax=(1407,'K')),
            NASAPolynomial(coeffs=[13.5268,0.0154307,-5.17464e-06,7.92549e-10,-4.55415e-14,-13394.6,-44.0461], Tmin=(1407,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 55,
    label = "C3H5O",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,D} {7,S}
3 C u0 p0 c0 {2,D} {8,S} {9,S}
4 O u1 p2 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.824069,0.034675,-2.51787e-05,9.56782e-09,-1.48085e-12,10420.4,22.8283], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[10.2638,0.011761,-3.89838e-06,5.92651e-10,-3.38867e-14,7259.38,-27.5109], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""KPS12""",
    longDesc = 
u"""
KPS12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 56,
    label = "NC3H7O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  O u1 p2 c0 {3,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.57487,0.0307101,-1.20049e-05,3.40807e-12,7.25275e-13,-6209.13,14.5966], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[11.5279,0.0153776,-5.23946e-06,8.11383e-10,-4.69928e-14,-9851,-35.4233], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 57,
    label = "C3H6OOH2-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {10,S} {11,S}
4  O u0 p2 c0 {1,S} {5,S}
5  O u0 p2 c0 {4,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.38466,0.0442929,-3.50977e-05,1.53695e-08,-2.81168e-12,-1809.8,16.9923], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[13.6645,0.015433,-5.29286e-06,8.23001e-10,-4.77931e-14,-5582.96,-42.8758], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 58,
    label = "C2H5CHO",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,D} {10,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  O u0 p2 c0 {3,D}
10 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.18896,0.025829,-6.0417e-06,-3.70703e-09,1.57131e-12,-24267.1,16.1496], Tmin=(298,'K'), Tmax=(1449,'K')),
            NASAPolynomial(coeffs=[10.6224,0.0135569,-4.60755e-06,7.12755e-10,-4.12632e-14,-27869.2,-31.6629], Tmin=(1449,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 59,
    label = "IC3H5OCH2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,D} {5,S}
3  C u0 p0 c0 {2,D} {9,S} {10,S}
4  C u1 p0 c0 {5,S} {11,S} {12,S}
5  O u0 p2 c0 {2,S} {4,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.14799,0.0700226,-8.21595e-05,5.2259e-08,-1.34177e-11,3264.75,35.5146], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.64732,0.0308191,-1.73209e-05,5.011e-09,-6.00089e-13,1706.8,-5.91215], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = u"""6/ 2/14 CZHOU""",
    longDesc = 
u"""
6/ 2/14 CZHOU
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
[CH2]OC(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 60,
    label = "C2H5CO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 C u1 p0 c0 {1,S} {9,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 O u0 p2 c0 {3,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[6.25722,-0.00917612,7.6119e-05,-9.05515e-08,3.46198e-11,-5916.16,2.23331], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.52325,0.0154212,-5.50898e-06,8.8589e-10,-5.28846e-14,-7196.32,-5.19862], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""A10/04""",
    longDesc = 
u"""
A10/04.
CC[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 61,
    label = "CH3O2H",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 O u0 p2 c0 {1,S} {3,S}
3 O u0 p2 c0 {2,S} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.90541,0.0174995,5.28244e-06,-2.52827e-08,1.34368e-11,-16889.5,11.3742], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.76538,0.008615,-2.98007e-06,4.68638e-10,-2.75339e-14,-18298,-14.3993], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""A 7/05""",
    longDesc = 
u"""
A 7/05.
COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 62,
    label = "IC3H7O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  O u1 p2 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.36108,0.034565,-1.9458e-05,4.71537e-09,-2.64705e-13,-8287.91,13.3112], Tmin=(298,'K'), Tmax=(1527,'K')),
            NASAPolynomial(coeffs=[11.9648,0.0142944,-4.71413e-06,7.14027e-10,-4.07161e-14,-11751.9,-38.8861], Tmin=(1527,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 63,
    label = "CH2CH2CHO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u1 p0 c0 {1,S} {6,S} {7,S}
3 C u0 p0 c0 {1,S} {8,D} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 O u0 p2 c0 {3,D}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.55799,0.0223392,-4.89741e-06,-3.58874e-09,1.47175e-12,453.128,16.7016], Tmin=(298,'K'), Tmax=(1437,'K')),
            NASAPolynomial(coeffs=[10.0673,0.0114971,-3.90138e-06,6.03029e-10,-3.48958e-14,-2750.81,-25.8818], Tmin=(1437,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 64,
    label = "NC3H7O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  O u0 p2 c0 {2,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u1 p2 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.13312,0.0396692,-2.3757e-05,6.9602e-09,-7.82577e-13,-7466.87,19.2445], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[13.2753,0.0161303,-5.52348e-06,8.58197e-10,-4.98173e-14,-11603.3,-41.5091], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 65,
    label = "C3H6OH1-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {2,S} {10,S}
4  O u0 p2 c0 {1,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.505208,0.036387,-2.15531e-05,6.45585e-09,-7.71267e-13,-9269.81,27.9804], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[10.0338,0.0160227,-5.41658e-06,8.34191e-10,-4.81216e-14,-12791.2,-23.9034], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 66,
    label = "C2H3OO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 O u0 p2 c0 {1,S} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 O u1 p2 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.09785,0.0295333,-2.27744e-05,7.20559e-09,-3.07929e-13,11399.6,21.3564], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.04484,0.0145511,-7.50975e-06,1.83488e-09,-1.6669e-13,10169.9,-3.71145], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C=CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 67,
    label = "CH3COCH3",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {10,D}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 O u0 p2 c0 {3,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.20008,0.027402,-1.31342e-05,2.5715e-09,-6.21509e-14,-27993.4,15.5884], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[8.87619,0.01457,-4.84823e-06,7.38615e-10,-4.22831e-14,-30604.6,-21.273], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 68,
    label = "C2H5O",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3 O u1 p2 c0 {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.90354,0.0177257,-2.69625e-06,-3.45831e-09,1.25225e-12,-3289.3,11.3546], Tmin=(298,'K'), Tmax=(1467,'K')),
            NASAPolynomial(coeffs=[8.19121,0.0110392,-3.75271e-06,5.80276e-10,-3.35735e-14,-5668.47,-19.0131], Tmin=(1467,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 69,
    label = "C3H6CHO-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {2,S} {10,S}
4  C u0 p0 c0 {1,S} {11,D} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 O u0 p2 c0 {4,D}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.01373,0.0233174,5.72464e-06,-1.27183e-08,3.69913e-12,-4167.66,13.0546], Tmin=(298,'K'), Tmax=(1439,'K')),
            NASAPolynomial(coeffs=[12.7129,0.0168633,-5.8378e-06,9.1406e-10,-5.33576e-14,-8501.66,-38.462], Tmin=(1439,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 70,
    label = "PC4H8CHO-4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u1 p0 c0 {2,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {14,D} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 O u0 p2 c0 {5,D}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.89963,0.0483599,-2.37381e-05,3.58762e-09,3.40202e-13,-5387.02,23.184], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[16.6333,0.020362,-6.95845e-06,1.08002e-09,-6.26592e-14,-11082.5,-58.1167], Tmin=(1383,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 71,
    label = "C5H91-5",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {5,D} {10,S}
4  C u1 p0 c0 {2,S} {11,S} {12,S}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.207677,0.0496049,-3.00622e-05,9.19962e-09,-1.13246e-12,20125.2,28.5856], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[13.7335,0.0207003,-7.06963e-06,1.09639e-09,-6.35591e-14,15117.6,-45.0795], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 72,
    label = "CH3NO2",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 N u0 p0 c+1 {1,S} {6,D} {7,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 O u0 p2 c0 {2,D}
7 O u0 p3 c-1 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.54054,0.0018656,4.44947e-05,-5.87057e-08,2.30684e-11,-11138.6,10.6885], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.73035,0.0109601,-4.05358e-06,6.67102e-10,-4.04687e-14,-12914.3,-10.1801], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T01/00""",
    longDesc = 
u"""
T01/00.
C[N+](=O)[O-]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 73,
    label = "C4H8-1",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.1626,0.0401053,-2.18039e-05,5.47071e-09,-4.54073e-13,-1654.03,24.8169], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[11.0189,0.0182714,-6.21802e-06,9.62039e-10,-5.56791e-14,-5809.99,-34.7942], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 74,
    label = "C5H10-1",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {5,D} {13,S}
5  C u0 p0 c0 {4,D} {14,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.165024,0.0530727,-3.10862e-05,8.92413e-09,-9.8162e-13,-4573.63,28.057], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[14.3625,0.0226076,-7.70501e-06,1.1933e-09,-6.91126e-14,-9999.16,-51.2512], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 75,
    label = "NC3H7O2H",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  O u0 p2 c0 {2,S} {5,S}
5  O u0 p2 c0 {4,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.35816,0.0456684,-2.91646e-05,9.41701e-09,-1.22337e-12,-24152.8,22.3323], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[14.2246,0.0174341,-5.97064e-06,9.27754e-10,-5.38585e-14,-28816,-47.4358], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 76,
    label = "C5H91-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u1 p0 c0 {1,S} {4,S} {11,S}
4  C u0 p0 c0 {3,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.46816,0.0513473,-3.05649e-05,8.8281e-09,-9.61458e-13,12378.1,28.3588], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[13.9594,0.020918,-7.14307e-06,1.10781e-09,-6.42268e-14,7026.51,-50.3104], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[CH]CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 77,
    label = "C4H71-4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,D} {7,S}
3  C u1 p0 c0 {1,S} {8,S} {9,S}
4  C u0 p0 c0 {2,D} {10,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.536903,0.0366356,-2.07815e-05,5.74895e-09,-6.05743e-13,23064.5,25.337], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[10.3875,0.0163677,-5.58416e-06,8.65389e-10,-5.01415e-14,19328.3,-28.6081], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 78,
    label = "C4H71-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u1 p0 c0 {1,S} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {10,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.94035,0.035683,-1.74385e-05,2.78965e-09,1.78069e-13,14930.3,21.1349], Tmin=(298,'K'), Tmax=(1367,'K')),
            NASAPolynomial(coeffs=[11.6978,0.0153405,-5.16929e-06,7.95431e-10,-4.58914e-14,10739.5,-38.2993], Tmin=(1367,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/13/16""",
    longDesc = 
u"""
1/13/16
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[CH]C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 79,
    label = "C5H91-4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u1 p0 c0 {1,S} {2,S} {11,S}
4  C u0 p0 c0 {1,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.58798,0.0401577,-1.50063e-05,-3.94608e-10,1.02064e-12,18468.4,23.4273], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[12.923,0.0210932,-7.14227e-06,1.10138e-09,-6.35984e-14,13819.2,-40.0293], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC[CH]C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 80,
    label = "C3H3",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {6,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.35111,0.0327411,-4.73827e-05,3.7631e-08,-1.18541e-11,40768,15.2059], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.14222,0.00761902,-2.6746e-06,4.24915e-10,-2.51475e-14,39571,-12.5849], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 7/11""",
    longDesc = 
u"""
T 7/11.
C#C[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 81,
    label = "C4H8-2",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {11,S}
4  C u0 p0 c0 {2,S} {3,D} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.30796,0.0353137,-1.51866e-05,1.64112e-09,3.44258e-13,-3197.68,18.1595], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[10.8652,0.0184123,-6.26887e-06,9.70206e-10,-5.61639e-14,-7096.26,-35.1547], Tmin=(1383,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 82,
    label = "C5H10-2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {15,S}
5  C u0 p0 c0 {3,S} {4,D} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.590529,0.0485275,-2.43232e-05,4.86096e-09,-1.13099e-13,-5997.78,24.1816], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[13.9426,0.0228735,-7.778e-06,1.20285e-09,-6.95972e-14,-11216.6,-49.538], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CCC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 83,
    label = "C8H142-6",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {7,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {8,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {7,D} {20,S}
6  C u0 p0 c0 {2,S} {8,D} {21,S}
7  C u0 p0 c0 {3,S} {5,D} {19,S}
8  C u0 p0 c0 {4,S} {6,D} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.343525,0.0772359,-4.34861e-05,1.12949e-08,-9.77995e-13,-837.798,29.8241], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[22.0961,0.0327887,-1.117e-05,1.7297e-09,-1.00179e-13,-9049.52,-89.2943], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/25/15""",
    longDesc = 
u"""
8/25/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CCCC=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 84,
    label = "IC3H7O2H",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  O u0 p2 c0 {1,S} {5,S}
5  O u0 p2 c0 {4,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.77385,0.0475813,-3.43745e-05,1.31405e-08,-2.06923e-12,-26345.9,17.767], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[14.4896,0.0168268,-5.67601e-06,8.72851e-10,-5.02994e-14,-30647.8,-50.1352], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 85,
    label = "C4H71-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {10,S}
4  C u1 p0 c0 {3,D} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.897231,0.0377004,-2.33195e-05,7.38468e-09,-9.50028e-13,27649.8,21.9835], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[11.0532,0.0155669,-5.25853e-06,8.09627e-10,-4.67015e-14,23945.6,-33.1548], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=CCC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 86,
    label = "C6H6",
    molecule = 
"""
1  C u0 p0 c0 {2,D} {6,S} {7,S}
2  C u0 p0 c0 {1,D} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {1,S} {5,D} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.504819,0.0185021,7.38346e-05,-1.18136e-07,5.0721e-11,8552.48,21.6413], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[11.081,0.0207177,-7.52146e-06,1.22321e-09,-7.36091e-14,4306.41,-40.0413], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""G 6/01""",
    longDesc = 
u"""
G 6/01.
C1=CC=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 87,
    label = "FULVENE",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,D}
2  C u0 p0 c0 {1,S} {4,D} {7,S}
3  C u0 p0 c0 {1,S} {5,D} {10,S}
4  C u0 p0 c0 {2,D} {5,S} {8,S}
5  C u0 p0 c0 {3,D} {4,S} {9,S}
6  C u0 p0 c0 {1,D} {11,S} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.718132,0.0379343,1.13988e-05,-4.13335e-08,1.80559e-11,24223.8,27.8557], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[11.1035,0.0206007,-7.53022e-06,1.23887e-09,-7.5416e-14,20361.8,-36.6652], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""0""",
    longDesc = 
u"""
0.
C=C1C=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 88,
    label = "C6H5",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,D} {3,S} {8,S}
2  C u0 p0 c0 {1,D} {4,S} {7,S}
3  C u0 p0 c0 {1,S} {5,D} {9,S}
4  C u0 p0 c0 {2,S} {6,D} {10,S}
5  C u0 p0 c0 {3,D} {6,S} {11,S}
6  C u1 p0 c0 {4,D} {5,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.210307,0.0204746,5.89743e-05,-1.01534e-07,4.47106e-11,39546.9,25.291], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[10.8445,0.0173212,-6.29233e-06,1.0237e-09,-6.16217e-14,35559.8,-35.3735], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T04/02""",
    longDesc = 
u"""
T04/02.
[C]1=CC=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 89,
    label = "C2H3CHO",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {6,S} {7,S}
3 C u0 p0 c0 {1,S} {5,D} {8,S}
4 H u0 p0 c0 {1,S}
5 O u0 p2 c0 {3,D}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.733844,0.0317483,-2.29599e-05,8.42104e-09,-1.23613e-12,-9384.74,21.0309], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[9.99155,0.00982348,-3.31203e-06,5.09524e-10,-2.93822e-14,-12530.4,-28.5169], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""KPS12""",
    longDesc = 
u"""
KPS12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 90,
    label = "C4H7O2-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {11,S}
4  C u0 p0 c0 {2,S} {3,D} {12,S}
5  O u1 p2 c0 {2,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.16435,0.0253428,-1.39669e-06,-6.00714e-09,1.76697e-12,5032.11,4.37441], Tmin=(298,'K'), Tmax=(1684,'K')),
            NASAPolynomial(coeffs=[11.2686,0.0201725,-7.41926e-06,1.21276e-09,-7.30115e-14,2138.76,-31.5208], Tmin=(1684,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 91,
    label = "C4H612",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {4,D} {9,S} {10,S}
4  C u0 p0 c0 {2,D} {3,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.945516,0.0346162,-1.98591e-05,5.02139e-09,-3.67977e-13,18143.9,20.2191], Tmin=(298,'K'), Tmax=(1374,'K')),
            NASAPolynomial(coeffs=[11.406,0.013149,-4.43542e-06,6.83029e-10,-3.94289e-14,14242.7,-36.9674], Tmin=(1374,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""A 8/83""",
    longDesc = 
u"""
A 8/83
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 92,
    label = "C4H6",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,D} {5,S}
2  C u0 p0 c0 {1,S} {4,D} {6,S}
3  C u0 p0 c0 {1,D} {7,S} {8,S}
4  C u0 p0 c0 {2,D} {9,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.01356,0.0335723,-1.96279e-05,5.74804e-09,-6.75029e-13,13395.7,19.9957], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[10.1065,0.0146248,-5.01374e-06,7.79511e-10,-4.52676e-14,9961.34,-29.7311], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 93,
    label = "C8H132-6,SA",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
2  C u0 p0 c0 {7,S} {11,S} {12,S} {13,S}
3  C u0 p0 c0 {6,S} {14,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {6,D} {19,S}
5  C u1 p0 c0 {1,S} {8,S} {17,S}
6  C u0 p0 c0 {3,S} {4,D} {20,S}
7  C u0 p0 c0 {2,S} {8,D} {18,S}
8  C u0 p0 c0 {5,S} {7,D} {21,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.220018,0.0748621,-4.2124e-05,1.07414e-08,-8.6656e-13,16102.8,29.9517], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[21.6167,0.0312261,-1.0664e-05,1.6541e-09,-9.59112e-14,8033.85,-87.2194], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C[CH]CC=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 94,
    label = "C4H5-N",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {5,S}
2 C u0 p0 c0 {1,S} {4,D} {6,S}
3 C u0 p0 c0 {1,D} {7,S} {8,S}
4 C u1 p0 c0 {2,D} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.163053,0.0398301,-3.40001e-05,1.51472e-08,-2.46658e-12,41429.8,23.5362], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[9.8502,0.010779,-1.36721e-06,-7.72005e-10,1.83663e-13,38840.3,-26.0018], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""H6W/94""",
    longDesc = 
u"""
H6W/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=CC=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 95,
    label = "C4H5-I",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {4,S} {5,S}
2 C u0 p0 c0 {1,D} {6,S} {7,S}
3 C u0 p0 c0 {4,D} {8,S} {9,S}
4 C u1 p0 c0 {1,S} {3,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.0199329,0.0380057,-2.75594e-05,7.78356e-09,4.02094e-13,37496.2,24.3942], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[10.2291,0.00948501,-9.04064e-08,-1.25961e-09,2.47815e-13,34642.8,-28.5645], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""H6W/94""",
    longDesc = 
u"""
H6W/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=[C]C=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 96,
    label = "C8H141-5,3-4",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {7,D} {17,S}
6  C u0 p0 c0 {2,S} {8,D} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  C u0 p0 c0 {6,D} {21,S} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.52713,0.0928838,-6.67177e-05,2.53684e-08,-4.01213e-12,-45.5966,42.3181], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[23.0691,0.0323154,-1.10796e-05,1.72276e-09,-1.00053e-13,-8912.4,-94.9656], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)C(C)C=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 97,
    label = "C8H141-5,3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {7,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {8,D} {18,S}
6  C u0 p0 c0 {2,S} {7,D} {19,S}
7  C u0 p0 c0 {4,S} {6,D} {20,S}
8  C u0 p0 c0 {5,D} {21,S} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.02072,0.0847809,-5.47568e-05,1.81792e-08,-2.47336e-12,-453.249,36.4247], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[22.5174,0.0326749,-1.11808e-05,1.7363e-09,-1.00752e-13,-8960.63,-91.0858], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/25/15""",
    longDesc = 
u"""
8/25/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)CC=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 98,
    label = "AC3H5OCH2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u1 p0 c0 {1,S} {11,S} {12,S}
5  O u0 p2 c0 {2,S} {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.20089,0.0424761,-2.94115e-05,1.11825e-08,-1.8131e-12,7267.92,24.3629], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[11.988,0.0170264,-5.78617e-06,8.9415e-10,-5.16995e-14,3483.67,-33.5861], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C1COC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 99,
    label = "C3H3O2H",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2 C u0 p0 c0 {1,S} {4,T}
3 O u0 p2 c0 {1,S} {5,S}
4 C u0 p0 c0 {2,T} {8,S}
5 O u0 p2 c0 {3,S} {9,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.09787,0.0422718,-3.83969e-05,1.77405e-08,-3.27674e-12,10359.2,23.0652], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[13.8152,0.00862175,-3.0671e-06,4.88874e-10,-2.88888e-14,6291.83,-43.9151], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/31/13""",
    longDesc = 
u"""
1/31/13
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 100,
    label = "C4H72-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {11,S}
4  C u1 p0 c0 {2,S} {3,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.465,0.0294957,-1.08905e-05,9.17747e-11,5.46418e-13,24290.4,14.4728], Tmin=(298,'K'), Tmax=(1378,'K')),
            NASAPolynomial(coeffs=[10.536,0.0165536,-5.71669e-06,8.93155e-10,-5.20433e-14,20816.1,-31.1046], Tmin=(1378,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C]=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 101,
    label = "C4H71-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {4,D} {10,S} {11,S}
4  C u1 p0 c0 {1,S} {3,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.56406,0.0332162,-1.59178e-05,2.92638e-09,-3.02645e-14,25796.6,19.3052], Tmin=(298,'K'), Tmax=(1381,'K')),
            NASAPolynomial(coeffs=[10.7106,0.0163539,-5.63688e-06,8.79592e-10,-5.12099e-14,22101.1,-31.53], Tmin=(1381,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=[C]CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 102,
    label = "SC4H9",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {3,S} {13,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.942662,0.0377415,-1.58912e-05,1.75489e-09,2.89726e-13,6205.43,24.2127], Tmin=(298,'K'), Tmax=(1682,'K')),
            NASAPolynomial(coeffs=[9.25139,0.0224301,-7.82649e-06,1.23559e-09,-7.2625e-14,3111.49,-21.608], Tmin=(1682,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 103,
    label = "C4H4",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,S} {5,S}
2 C u0 p0 c0 {1,D} {6,S} {7,S}
3 C u0 p0 c0 {1,S} {4,T}
4 C u0 p0 c0 {3,T} {8,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.91525,0.0527509,-7.16559e-05,5.50724e-08,-1.72862e-11,32978.5,31.42], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.65071,0.0161294,-7.19389e-06,1.49818e-09,-1.18641e-13,31196,-9.79521], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""H6W/94""",
    longDesc = 
u"""
H6W/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#CC=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 104,
    label = "C2H3CO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 C u1 p0 c0 {1,S} {7,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 O u0 p2 c0 {3,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.65335,0.0257403,-1.8901e-05,7.29175e-09,-1.16083e-12,10202.1,17.8706], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[8.86033,0.00848985,-2.9035e-06,4.50764e-10,-2.61524e-14,7734.89,-20.6979], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""KPS12""",
    longDesc = 
u"""
KPS12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 105,
    label = "C4H6-2",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,T}
4  C u0 p0 c0 {2,S} {3,T}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.97152,0.0276791,-1.13397e-05,1.02971e-09,2.75291e-13,15728.4,14.2147], Tmin=(298,'K'), Tmax=(1377,'K')),
            NASAPolynomial(coeffs=[9.60306,0.0148972,-5.16751e-06,8.09757e-10,-4.72818e-14,12483.1,-28.713], Tmin=(1377,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""A 8/83""",
    longDesc = 
u"""
A 8/83
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC#CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 106,
    label = "C5H92-5",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {11,S}
4  C u0 p0 c0 {2,S} {3,D} {12,S}
5  C u1 p0 c0 {1,S} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.996364,0.0449149,-2.30951e-05,5.01699e-09,-2.38646e-13,18716.2,24.5616], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[13.2989,0.0209846,-7.14999e-06,1.10717e-09,-6.41184e-14,13927,-43.2759], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 107,
    label = "C5H92-4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
3  C u1 p0 c0 {1,S} {5,S} {12,S}
4  C u0 p0 c0 {2,S} {5,D} {13,S}
5  C u0 p0 c0 {3,S} {4,D} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.255028,0.0469861,-2.40448e-05,4.89007e-09,-1.15558e-13,10977.5,24.6227], Tmin=(298,'K'), Tmax=(1384,'K')),
            NASAPolynomial(coeffs=[13.5638,0.0211621,-7.20826e-06,1.11611e-09,-6.4637e-14,5822.59,-48.7324], Tmin=(1384,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 108,
    label = "C8H132-6,PA",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {6,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {6,D} {17,S}
5  C u0 p0 c0 {1,S} {7,D} {16,S}
6  C u0 p0 c0 {3,S} {4,D} {18,S}
7  C u0 p0 c0 {5,D} {8,S} {19,S}
8  C u1 p0 c0 {7,S} {20,S} {21,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.10512,0.077015,-4.48891e-05,1.19517e-08,-1.04378e-12,17271.3,32.2059], Tmin=(298,'K'), Tmax=(1378,'K')),
            NASAPolynomial(coeffs=[22.9294,0.0297187,-1.00712e-05,1.55505e-09,-8.9913e-14,8779.76,-92.4517], Tmin=(1378,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CCCC=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 109,
    label = "C6H9-A",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u1 p0 c0 {1,S} {4,S} {9,S}
3  C u0 p0 c0 {1,S} {5,D} {10,S}
4  C u0 p0 c0 {2,S} {6,D} {11,S}
5  C u0 p0 c0 {3,D} {14,S} {15,S}
6  C u0 p0 c0 {4,D} {12,S} {13,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.66715,0.0726196,-6.05324e-05,2.66001e-08,-4.74613e-12,26441.5,40.222], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[17.0843,0.0208843,-7.14529e-06,1.10944e-09,-6.43677e-14,20104,-63.9326], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""12/ 5/12 THERM""",
    longDesc = 
u"""
12/ 5/12 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[CH]CC=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 110,
    label = "C8H131-5,3-4,TA",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
3  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
4  C u1 p0 c0 {1,S} {3,S} {6,S}
5  C u0 p0 c0 {1,S} {7,D} {16,S}
6  C u0 p0 c0 {4,S} {8,D} {17,S}
7  C u0 p0 c0 {5,D} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {18,S} {19,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.23997,0.0913495,-6.60037e-05,2.51706e-08,-3.99031e-12,16434.5,44.9499], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[22.2823,0.0310652,-1.06941e-05,1.6674e-09,-9.70243e-14,7576.81,-91.9859], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[C](C)C(C)C=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 111,
    label = "C8H131-5,3,PA",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {7,D} {16,S}
5  C u0 p0 c0 {2,S} {6,D} {15,S}
6  C u0 p0 c0 {5,D} {8,S} {17,S}
7  C u0 p0 c0 {4,D} {20,S} {21,S}
8  C u1 p0 c0 {6,S} {18,S} {19,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.57063,0.0860081,-5.83884e-05,2.02048e-08,-2.83116e-12,17701.5,39.5476], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[23.2086,0.0297529,-1.01334e-05,1.56926e-09,-9.09001e-14,8953.39,-94.0778], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CCC(C)C=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 112,
    label = "C8H131-5,3,SA",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
3  C u0 p0 c0 {6,S} {13,S} {14,S} {15,S}
4  C u1 p0 c0 {1,S} {7,S} {16,S}
5  C u0 p0 c0 {1,S} {8,D} {17,S}
6  C u0 p0 c0 {3,S} {7,D} {18,S}
7  C u0 p0 c0 {4,S} {6,D} {19,S}
8  C u0 p0 c0 {5,D} {20,S} {21,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.32186,0.0830474,-5.42213e-05,1.80737e-08,-2.45105e-12,16518.3,36.717], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[22.1142,0.0309871,-1.06198e-05,1.651e-09,-9.58788e-14,8085.1,-90.1448], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)[CH]C=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 113,
    label = "C8H131-5,3,TA",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
2  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
3  C u0 p0 c0 {6,S} {14,S} {15,S} {16,S}
4  C u1 p0 c0 {1,S} {2,S} {7,S}
5  C u0 p0 c0 {1,S} {6,D} {17,S}
6  C u0 p0 c0 {3,S} {5,D} {18,S}
7  C u0 p0 c0 {4,S} {8,D} {19,S}
8  C u0 p0 c0 {7,D} {20,S} {21,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.503201,0.0778854,-4.63156e-05,1.32834e-08,-1.42623e-12,15838.9,32.6514], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[21.6162,0.0314327,-1.07786e-05,1.67633e-09,-9.73755e-14,7608.97,-88.0384], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[C](C)CC=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 114,
    label = "C-C6H4",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,D} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {1,D} {6,S} {9,S}
4  C u0 p0 c0 {2,D} {5,S} {10,S}
5  C u0 p0 c0 {4,S} {6,T}
6  C u0 p0 c0 {3,S} {5,T}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.09913,0.0540306,-4.0839e-05,1.07388e-08,9.80785e-13,52205.7,37.4152], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[13.8492,0.00788079,1.82438e-06,-2.11692e-09,3.746e-13,47446.3,-50.405], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""H6W/94""",
    longDesc = 
u"""
H6W/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C1#CC=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 115,
    label = "C4H10",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.0920862,0.0469704,-2.54762e-05,6.35895e-09,-5.16006e-13,-16955.7,24.9102], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[12.4924,0.0215952,-7.34278e-06,1.1353e-09,-6.5673e-14,-21759.9,-44.1547], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 116,
    label = "C5H81-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,D} {9,S}
3  C u0 p0 c0 {2,D} {4,S} {11,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {4,D} {12,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.54882,0.0415043,-2.1436e-05,4.71146e-09,-2.42143e-13,9056.36,19.5666], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[12.9945,0.0192678,-6.58967e-06,1.02296e-09,-5.93441e-14,4590.4,-43.569], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 117,
    label = "C4H3-N",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {5,S}
2 C u0 p0 c0 {1,S} {4,T}
3 C u1 p0 c0 {1,D} {6,S}
4 C u0 p0 c0 {2,T} {7,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.316841,0.0469121,-6.80938e-05,5.31799e-08,-1.6523e-11,62476.2,24.6226], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.43283,0.016861,-9.43131e-06,2.57039e-09,-2.74563e-13,61600.7,-1.5674], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""H6W/94""",
    longDesc = 
u"""
H6W/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=CC#C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 118,
    label = "C2H5OH",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 O u0 p2 c0 {1,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.215806,0.0295228,-1.68271e-05,4.49485e-09,-4.02452e-13,-29485.2,24.5725], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[8.14484,0.0128314,-4.29053e-06,6.55972e-10,-3.76507e-14,-32400.6,-18.6241], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 119,
    label = "C4H72-1OOH",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,D} {13,S}
4  C u0 p0 c0 {2,S} {3,D} {12,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.29755,0.0559252,-4.0889e-05,1.54881e-08,-2.42412e-12,-11604.7,24.3621], Tmin=(298,'K'), Tmax=(1381,'K')),
            NASAPolynomial(coeffs=[18.0123,0.0170341,-5.89884e-06,9.23962e-10,-5.3954e-14,-17458.5,-65.521], Tmin=(1381,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 120,
    label = "IC3H6CHO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u1 p0 c0 {1,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {11,D} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 O u0 p2 c0 {4,D}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.521482,0.0443114,-2.86617e-05,9.3032e-09,-1.20762e-12,-2996.77,26.8182], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[13.3102,0.0162098,-5.57576e-06,8.69004e-10,-5.05554e-14,-7621.78,-42.5051], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/22/96 THERM""",
    longDesc = 
u"""
2/22/96 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 121,
    label = "SC3H5CHO",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,D} {8,S}
3  C u0 p0 c0 {2,D} {4,S} {9,S}
4  C u0 p0 c0 {3,S} {10,D} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 O u0 p2 c0 {4,D}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.09373,0.0443315,-3.41918e-05,1.3937e-08,-2.33791e-12,-15674.6,19.4458], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[13.3892,0.0139115,-4.75821e-06,7.38737e-10,-4.28607e-14,-19791.7,-46.0146], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 122,
    label = "C5H9O1-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {5,D} {13,S}
5  C u0 p0 c0 {4,D} {14,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  O u1 p2 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.49666,0.0520326,-3.09828e-05,9.02567e-09,-1.04311e-12,1547.59,17.0868], Tmin=(298,'K'), Tmax=(1379,'K')),
            NASAPolynomial(coeffs=[17.886,0.0204837,-7.16613e-06,1.12949e-09,-6.62221e-14,-4331.4,-67.2851], Tmin=(1379,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC([O])CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 123,
    label = "C8H132-6,SAO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
3  C u0 p0 c0 {7,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {8,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {7,D} {20,S}
6  C u0 p0 c0 {2,S} {8,D} {21,S}
7  C u0 p0 c0 {3,S} {5,D} {19,S}
8  C u0 p0 c0 {4,S} {6,D} {22,S}
9  O u1 p2 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.34169,0.0749918,-4.18383e-05,1.06108e-08,-8.96488e-13,5244.58,18.6128], Tmin=(298,'K'), Tmax=(1381,'K')),
            NASAPolynomial(coeffs=[25.4268,0.0310422,-1.08046e-05,1.69717e-09,-9.92702e-14,-3286.42,-102.873], Tmin=(1381,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CCC([O])C=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 124,
    label = "C8H131-5,3,SAO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {7,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {8,D} {18,S}
6  C u0 p0 c0 {2,S} {7,D} {19,S}
7  C u0 p0 c0 {4,S} {6,D} {20,S}
8  C u0 p0 c0 {5,D} {21,S} {22,S}
9  H u0 p0 c0 {1,S}
10 O u1 p2 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.60954,0.0839739,-5.49469e-05,1.84255e-08,-2.5604e-12,5337.13,26.2654], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[26.1024,0.0305037,-1.06268e-05,1.67027e-09,-9.77389e-14,-3647.44,-106.781], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)C([O])C=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 125,
    label = "CH2CHCHCHO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,D} {3,S} {5,S}
2  C u0 p0 c0 {1,D} {4,S} {6,S}
3  C u1 p0 c0 {1,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {9,D} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  O u0 p2 c0 {4,D}
10 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.451553,0.043753,-3.37274e-05,1.31111e-08,-2.04732e-12,2149.43,24.0836], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[13.93,0.0113229,-3.87394e-06,6.02379e-10,-3.5012e-14,-2395.9,-47.9027], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 126,
    label = "PC4H9",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {2,S} {12,S} {13,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.409645,0.0429511,-2.36583e-05,6.15745e-09,-5.64301e-13,7743.19,25.5313], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[11.8548,0.0196962,-6.71054e-06,1.03891e-09,-6.01514e-14,3381.82,-37.2343], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 127,
    label = "CVCCVCCJ",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,D} {3,S} {7,S}
2  C u0 p0 c0 {1,D} {4,S} {6,S}
3  C u0 p0 c0 {1,S} {5,D} {8,S}
4  C u1 p0 c0 {2,S} {9,S} {10,S}
5  C u0 p0 c0 {3,D} {11,S} {12,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.60087,0.0538765,-3.96302e-05,1.49599e-08,-2.31995e-12,23120,33.5493], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[14.7303,0.0159031,-5.5773e-06,8.80605e-10,-5.16964e-14,17405.1,-54.2671], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""Z&B""",
    longDesc = 
u"""
Z&B
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CC=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 128,
    label = "PC2H4OH",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u1 p0 c0 {1,S} {6,S} {7,S}
3 O u0 p2 c0 {1,S} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5948,0.0227101,-1.39474e-05,4.70096e-09,-6.90044e-13,-4914.87,14.3241], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[8.0675,0.0106144,-3.57999e-06,5.50364e-10,-3.17052e-14,-6927.48,-15.3833], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 129,
    label = "C5H11-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {2,S} {4,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.817689,0.0492654,-2.13786e-05,1.85532e-09,7.06259e-13,3262.21,26.5876], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[14.7177,0.0241668,-8.18648e-06,1.26272e-09,-7.29268e-14,-2260.28,-50.6071], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 130,
    label = "AC3H5OOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,D} {8,S}
3  C u0 p0 c0 {2,D} {9,S} {10,S}
4  O u0 p2 c0 {1,S} {5,S}
5  O u0 p2 c0 {4,S} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.18125,0.0435233,-5.16277e-05,4.32011e-08,-1.57715e-11,-7635.22,12.1726], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[12.0839,0.0147947,-5.13213e-06,8.07505e-10,-4.74395e-14,-10218.4,-33.6435], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""GOLDSMITH""",
    longDesc = 
u"""
GOLDSMITH.
C=CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 131,
    label = "C5H11O2-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {12,S} {13,S} {14,S}
6  O u0 p2 c0 {1,S} {18,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.95945,0.0673013,-4.74262e-05,1.76453e-08,-2.71599e-12,-15510.2,21.5326], Tmin=(298,'K'), Tmax=(1401,'K')),
            NASAPolynomial(coeffs=[20.1379,0.0243049,-8.24745e-06,1.27345e-09,-7.35954e-14,-21777.8,-75.9194], Tmin=(1401,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC(C)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 132,
    label = "CH3COCH2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,S} {7,D}
3 C u1 p0 c0 {2,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 O u0 p2 c0 {2,D}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.13382,0.0325095,-2.10425e-05,6.64421e-09,-8.12619e-13,-6048.68,21.7159], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[10.9524,0.0111459,-3.86263e-06,6.05089e-10,-3.53293e-14,-9608.34,-31.5623], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/14/13 THERM""",
    longDesc = 
u"""
2/14/13 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 133,
    label = "CH3CHO",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,D} {7,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 O u0 p2 c0 {2,D}
7 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.72946,-0.00319329,4.75349e-05,-5.74586e-08,2.19311e-11,-21572.9,4.10302], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.40411,0.0117231,-4.22631e-06,6.83725e-10,-4.09849e-14,-22593.1,-3.48079], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 8/88""",
    longDesc = 
u"""
L 8/88.
CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 134,
    label = "SC4H9O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
5  O u0 p2 c0 {1,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 O u1 p2 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.27751,0.0544335,-3.82988e-05,1.42448e-08,-2.18865e-12,-12590.9,18.3237], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[16.8209,0.0198834,-6.71756e-06,1.03412e-09,-5.96371e-14,-17581.5,-59.5731], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 135,
    label = "C4H2",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,T}
2 C u0 p0 c0 {1,S} {4,T}
3 C u0 p0 c0 {1,T} {5,S}
4 C u0 p0 c0 {2,T} {6,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.0544,0.041627,-6.58718e-05,5.32571e-08,-1.66832e-11,54185.2,14.8666], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[9.15763,0.00554305,-1.35916e-06,1.87801e-11,2.31895e-14,52588,-23.7115], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""D11/99""",
    longDesc = 
u"""
D11/99
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#CC#C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 136,
    label = "C4H71-3OOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {11,S}
4  C u0 p0 c0 {3,D} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.50977,0.0685369,-5.75194e-05,2.43179e-08,-4.09788e-12,-11801.9,35.042], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[19.2985,0.0154534,-5.2546e-06,8.13772e-10,-4.7169e-14,-18500.3,-74.9927], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 137,
    label = "H2CNO2",
    molecule = 
"""
multiplicity 2
1 N u0 p0 c+1 {2,S} {3,D} {6,S}
2 C u1 p0 c0 {1,S} {4,S} {5,S}
3 O u0 p2 c0 {1,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 O u0 p3 c-1 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.42742,0.0160496,2.84728e-06,-1.82218e-08,9.35384e-12,14012.1,16.1086], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.57505,0.00701471,-2.51481e-06,4.05671e-10,-2.42797e-14,12388,-11.5986], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""RADICAL   T08/07""",
    longDesc = 
u"""
RADICAL   T08/07.
[N+](=O)([CH2])[O-]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 138,
    label = "C2H3OH",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 O u0 p2 c0 {1,S} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.127972,0.0338506,-3.30645e-05,1.64859e-08,-3.19935e-12,-15991.5,23.0439], Tmin=(298,'K'), Tmax=(1410,'K')),
            NASAPolynomial(coeffs=[8.32598,0.00803387,-2.63928e-06,3.98411e-10,-2.26551e-14,-18322.1,-20.208], Tmin=(1410,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/ 3/ 9 THERM""",
    longDesc = 
u"""
2/ 3/ 9 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 139,
    label = "C4H71-4O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  O u0 p2 c0 {2,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 O u1 p2 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.09043,0.045713,-2.94816e-05,9.69605e-09,-1.30061e-12,4889.33,21.5456], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[15.024,0.0172861,-5.94297e-06,9.25861e-10,-5.38462e-14,192.461,-48.5941], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/25/15""",
    longDesc = 
u"""
9/25/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 140,
    label = "C5H10OOH2-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
2  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  C u1 p0 c0 {1,S} {2,S} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.03845,0.0635523,-4.29455e-05,1.50666e-08,-2.16139e-12,-9066.91,18.1118], Tmin=(298,'K'), Tmax=(1410,'K')),
            NASAPolynomial(coeffs=[20.36,0.0236444,-8.00997e-06,1.23544e-09,-7.13457e-14,-15126.4,-75.0913], Tmin=(1410,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 141,
    label = "C5H10OOH2-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
3  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {2,S} {16,S} {17,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.97483,0.0685027,-4.88793e-05,1.80505e-08,-2.71516e-12,-7589.36,23.171], Tmin=(298,'K'), Tmax=(1401,'K')),
            NASAPolynomial(coeffs=[21.1684,0.0230996,-7.85422e-06,1.21454e-09,-7.02697e-14,-14165.3,-79.6524], Tmin=(1401,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(CCC)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 142,
    label = "C4H8OOH2-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.32205,0.0555032,-3.97238e-05,1.47354e-08,-2.22482e-12,-4673.32,19.8332], Tmin=(298,'K'), Tmax=(1404,'K')),
            NASAPolynomial(coeffs=[17.7648,0.0187777,-6.3641e-06,9.81958e-10,-5.67252e-14,-9938.52,-62.817], Tmin=(1404,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(CC)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 143,
    label = "CH2CH2CH2COCH3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {3,S} {13,D}
5  C u1 p0 c0 {2,S} {14,S} {15,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.03929,0.0476716,-2.60374e-05,6.44111e-09,-5.11323e-13,-9142.43,22.5011], Tmin=(298,'K'), Tmax=(2010,'K')),
            NASAPolynomial(coeffs=[13.6122,0.0229083,-7.89569e-06,1.23743e-09,-7.24019e-14,-13133.9,-40.0995], Tmin=(2010,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/6/15""",
    longDesc = 
u"""
10/6/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 144,
    label = "C2H5O2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 O u0 p2 c0 {1,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 O u1 p2 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.90352,0.0222599,-1.0161e-05,1.7171e-09,1.88167e-14,-5096.54,8.98723], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[9.50283,0.012043,-4.09492e-06,6.33049e-10,-3.66134e-14,-7370.69,-22.1717], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 145,
    label = "CC3H6",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
3 C u0 p0 c0 {1,S} {2,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.83279,-0.00521029,9.29583e-05,-1.22753e-07,4.99191e-11,5195.2,10.8306], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.21663,0.0165394,-5.90076e-06,9.48095e-10,-5.65662e-14,2959.37,-13.6041], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
C1CC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 146,
    label = "C4H71-4OOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {11,S}
4  C u0 p0 c0 {3,D} {12,S} {13,S}
5  O u0 p2 c0 {2,S} {6,S}
6  O u0 p2 c0 {5,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.31653,0.0516546,-3.47311e-05,1.20406e-08,-1.71651e-12,-11775.5,24.6385], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[15.9871,0.0186028,-6.40003e-06,9.97532e-10,-5.80334e-14,-17015.2,-54.6227], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 147,
    label = "CH2CHO",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 C u0 p0 c0 {1,S} {5,D} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 O u0 p2 c0 {2,D}
6 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.79503,0.0101099,1.61751e-05,-3.10303e-08,1.39436e-11,162.945,12.3647], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.53928,0.00780239,-2.76414e-06,4.42099e-10,-2.62954e-14,-1188.59,-8.72091], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T03/10""",
    longDesc = 
u"""
T03/10.
[CH2]C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 148,
    label = "CH3CO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u1 p0 c0 {1,S} {6,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 O u0 p2 c0 {2,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.03587,0.000877295,3.071e-05,-3.92476e-08,1.52969e-11,-2682.07,7.86177], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.31372,0.00917378,-3.32204e-06,5.39475e-10,-3.24524e-14,-3645.04,-1.67576], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""IU2/03""",
    longDesc = 
u"""
IU2/03.
C[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 149,
    label = "SC4H9O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  O u1 p2 c0 {2,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.01773,0.0470084,-2.74727e-05,7.3629e-09,-6.30414e-13,-11189.2,17.4372], Tmin=(298,'K'), Tmax=(1411,'K')),
            NASAPolynomial(coeffs=[15.213,0.019003,-6.39005e-06,9.80774e-10,-5.64494e-14,-15988.9,-54.3195], Tmin=(1411,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 150,
    label = "C4H71-O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  O u1 p2 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.60619,0.0558563,-4.35596e-05,1.70589e-08,-2.65635e-12,4850.9,34.7113], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[15.3138,0.0143427,-4.81626e-06,7.39575e-10,-4.26141e-14,-729.343,-55.2938], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""4/ 3/ 0 THERM""",
    longDesc = 
u"""
4/ 3/ 0 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 151,
    label = "CH2CO",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 O u0 p2 c0 {2,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.81423,0.0199009,-2.21416e-05,1.45029e-08,-3.98877e-12,-7053.95,13.6079], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.35869,0.00695642,-2.64803e-06,4.65068e-10,-3.08642e-14,-7902.94,-3.98526], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 152,
    label = "SC4H9O2H",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.44011,0.06053,-4.36263e-05,1.66226e-08,-2.61557e-12,-29263.1,21.7354], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[17.8076,0.0212546,-7.2096e-06,1.11291e-09,-6.4306e-14,-34845.6,-65.8011], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 153,
    label = "CH2(S)",
    molecule = 
"""
multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.19331,-0.00233105,8.15676e-06,-6.62986e-09,1.93233e-12,50366.2,-0.746734], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.13502,0.00289594,-8.16668e-07,1.13573e-10,-6.36263e-15,50504.1,4.06031], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""IU6/03""",
    longDesc = 
u"""
IU6/03.
Duplicate of species CH2 (i.e. same molecular structure according to RMG)
[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 154,
    label = "HOCHO",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,D} {4,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u0 p2 c0 {1,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.89836,-0.00355878,3.55205e-05,-4.385e-08,1.71078e-11,-46770.6,7.34954], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.61383,0.00644964,-2.29083e-06,3.6716e-10,-2.18737e-14,-47514.8,0.847884], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 8/88""",
    longDesc = 
u"""
L 8/88.
O=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 155,
    label = "AC5H9-D",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {5,D}
4  C u1 p0 c0 {1,S} {11,S} {12,S}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.167398,0.0509877,-3.22616e-05,1.05912e-08,-1.43705e-12,18179.3,29.5711], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[13.5608,0.0207515,-7.06609e-06,1.09362e-09,-6.33083e-14,13189.7,-44.8719], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 156,
    label = "NC5H12",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.299552,0.0594963,-3.41764e-05,9.47896e-09,-9.73675e-13,-19896,27.5742], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[15.8289,0.0259345,-8.83016e-06,1.36655e-09,-7.91029e-14,-25939.7,-60.5558], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 157,
    label = "NC3H7COCH3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {4,S} {16,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.22466,0.0482792,-2.2102e-05,2.6403e-09,5.2372e-13,-33937.9,19.4495], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[15.9633,0.0226395,-7.5643e-06,1.15632e-09,-6.63757e-14,-39282.2,-56.5033], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/6/15 THERM""",
    longDesc = 
u"""
10/6/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 158,
    label = "DC5H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {2,S} {15,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
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
            NASAPolynomial(coeffs=[-1.22778,0.0632926,-4.52076e-05,1.75515e-08,-2.87257e-12,4325.62,33.1703], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[15.2292,0.0238483,-8.09849e-06,1.25097e-09,-7.23134e-14,-1358.28,-54.9727], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 159,
    label = "IC4H10",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.07414,0.0524618,-3.42408e-05,1.18818e-08,-1.73238e-12,-17921.9,27.4852], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[12.6423,0.0214134,-7.26712e-06,1.12207e-09,-6.48434e-14,-22829.4,-46.606], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 160,
    label = "IC4H7",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u1 p0 c0 {2,S} {8,S} {9,S}
4  C u0 p0 c0 {2,D} {10,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.229579,0.0417843,-2.66886e-05,8.42206e-09,-1.03175e-12,14394.7,25.4798], Tmin=(298,'K'), Tmax=(1384,'K')),
            NASAPolynomial(coeffs=[11.8999,0.015157,-5.09995e-06,7.83722e-10,-4.5166e-14,10036.4,-40.2287], Tmin=(1384,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 161,
    label = "CC5H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {1,S} {4,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.197886,0.0535219,-2.95962e-05,7.61177e-09,-6.45329e-13,2665.43,27.8232], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[14.3629,0.024346,-8.21904e-06,1.26458e-09,-7.28991e-14,-2645.94,-49.6386], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C(C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 162,
    label = "DC5H11O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {3,S} {18,S}
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
17 H u0 p0 c0 {5,S}
18 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.354375,0.0721759,-5.35556e-05,2.12668e-08,-3.51578e-12,-13832,29.255], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[19.8701,0.0247772,-8.46414e-06,1.31284e-09,-7.61129e-14,-20483.2,-75.0051], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)CCO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 163,
    label = "PC4H9CHO",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,D} {16,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.55168,0.0517581,-2.47522e-05,3.39586e-09,4.54215e-13,-30110.5,22.534], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[17.1928,0.0224205,-7.66525e-06,1.19e-09,-6.90481e-14,-36204.8,-63.9265], Tmin=(1383,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 164,
    label = "C5H10OOH2-1O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {8,S}
7  O u0 p2 c0 {4,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.92124,0.0767602,-5.10112e-05,1.59643e-08,-1.8298e-12,-25548.7,22.7306], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[27.1106,0.0234082,-8.1172e-06,1.27231e-09,-7.43222e-14,-34140.6,-108.093], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC(CO[O])OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 165,
    label = "O2CHO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {4,S}
2 O u0 p2 c0 {1,S} {5,S}
3 O u0 p2 c0 {1,D}
4 H u0 p0 c0 {1,S}
5 O u1 p2 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.96059,0.0106002,-5.25713e-06,1.01717e-09,-2.87488e-14,-17359.9,11.7807], Tmin=(298,'K'), Tmax=(1368,'K')),
            NASAPolynomial(coeffs=[7.24075,0.00463313,-1.63694e-06,2.59707e-10,-1.52965e-14,-18702.8,-6.49547], Tmin=(1368,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""6/26/95 THERM""",
    longDesc = 
u"""
6/26/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 166,
    label = "IC5H12",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.61701,0.0668963,-4.65395e-05,1.75637e-08,-2.80513e-12,-20392.2,32.7017], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[15.8024,0.0258962,-8.80211e-06,1.36051e-09,-7.86805e-14,-26487.7,-60.8691], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 167,
    label = "HCCO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,T} {4,S}
2 C u0 p0 c0 {1,T} {3,S}
3 O u1 p2 c0 {2,S}
4 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.87608,0.0221205,-3.58869e-05,3.05403e-08,-1.01281e-11,20163.4,13.6968], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.91479,0.00371409,-1.30137e-06,2.06473e-10,-1.21477e-14,19359.6,-5.50567], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 4/09""",
    longDesc = 
u"""
T 4/09.
C#C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 168,
    label = "SC2H4OH",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {3,S} {7,S}
3 O u0 p2 c0 {2,S} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.46281,0.0239194,-1.30667e-05,3.10615e-09,-1.85896e-13,-8007.9,19.2547], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[8.15007,0.0102549,-3.40138e-06,5.1751e-10,-2.96129e-14,-10501.4,-17.3135], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 169,
    label = "C5H10OOH2-3O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {1,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.94959,0.0794697,-6.123e-05,2.45361e-08,-3.99244e-12,-27626.3,21.0066], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[25.1554,0.0242477,-8.22109e-06,1.26886e-09,-7.33158e-14,-34952.2,-96.9776], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(O[O])C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 170,
    label = "HCCOH",
    molecule = 
"""
1 C u0 p0 c0 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {4,S}
3 O u0 p2 c0 {1,S} {5,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.05541,0.0252003,-3.80822e-05,3.09891e-08,-9.898e-12,9768.72,12.2272], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.3751,0.00549429,-1.88137e-06,2.93804e-10,-1.71772e-14,8932.78,-8.24498], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T12/09""",
    longDesc = 
u"""
T12/09.
C#CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 171,
    label = "DC5H11O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {9,S} {16,S} {17,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  O u1 p2 c0 {5,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.843925,0.0629724,-4.12079e-05,1.37396e-08,-1.84098e-12,-12581.9,24.3859], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[18.2339,0.0237924,-8.07056e-06,1.24593e-09,-7.19996e-14,-18745,-69.4853], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)CC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 172,
    label = "IC3H5OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,D} {4,S}
3  C u0 p0 c0 {2,D} {8,S} {9,S}
4  O u0 p2 c0 {2,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.58376,0.0316215,-1.73665e-05,4.18928e-09,-2.799e-13,-21264.3,18.8314], Tmin=(298,'K'), Tmax=(1374,'K')),
            NASAPolynomial(coeffs=[10.7381,0.0131698,-4.4153e-06,6.7701e-10,-3.89609e-14,-24729.8,-31.3634], Tmin=(1374,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/ 1/95 THERM""",
    longDesc = 
u"""
8/ 1/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 173,
    label = "C2H5O2H",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  O u0 p2 c0 {1,S} {4,S}
4  O u0 p2 c0 {3,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.83755,0.0338054,-2.37548e-05,9.31975e-09,-1.58003e-12,-21581.4,18.0978], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[10.4824,0.013478,-4.62179e-06,7.18619e-10,-4.17307e-14,-24657.8,-28.4294], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 174,
    label = "CC5H11O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {2,S} {18,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.715982,0.0744249,-5.92759e-05,2.52116e-08,-4.38573e-12,-16019.8,25.6129], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[20.2076,0.0239554,-8.06287e-06,1.23798e-09,-7.12608e-14,-22307.4,-77.3358], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)C(C)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 175,
    label = "DC5H11O2H",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {19,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.476375,0.0783789,-5.91713e-05,2.38264e-08,-3.97732e-12,-30508.8,32.6083], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[20.8513,0.0260811,-8.91736e-06,1.38397e-09,-8.02717e-14,-37719.4,-81.1417], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 176,
    label = "CC5H11O2H",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {19,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.00630157,0.0801496,-6.42109e-05,2.73711e-08,-4.76269e-12,-32712.9,28.4639], Tmin=(298,'K'), Tmax=(1406,'K')),
            NASAPolynomial(coeffs=[21.1504,0.0253099,-8.53704e-06,1.31271e-09,-7.56412e-14,-39527.5,-83.2516], Tmin=(1406,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 177,
    label = "O2C2H4OH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  O u0 p2 c0 {1,S} {9,S}
4  O u0 p2 c0 {2,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  O u1 p2 c0 {3,S}
10 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[7.0401,0.0159564,2.21097e-06,-7.05197e-09,2.08266e-12,-22452.4,-1.75362], Tmin=(298,'K'), Tmax=(1506,'K')),
            NASAPolynomial(coeffs=[12.7504,0.0111514,-3.83474e-06,5.98156e-10,-3.48372e-14,-25277.1,-35.4318], Tmin=(1506,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12 THERM""",
    longDesc = 
u"""
9/ 1/12 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OCCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 178,
    label = "CH3OCHO",
    molecule = 
"""
1 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {3,S} {7,D} {8,S}
3 O u0 p2 c0 {1,S} {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 O u0 p2 c0 {2,D}
8 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.96757,-0.00938085,7.07648e-05,-8.29932e-08,3.13523e-11,-45571.3,0.750341], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.33361,0.0134851,-4.84306e-06,7.81719e-10,-4.67917e-14,-46831.7,-6.91543], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 6/08""",
    longDesc = 
u"""
T 6/08.
COC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 179,
    label = "C5H11O2H-2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {19,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.06841,0.0736472,-5.31234e-05,2.02381e-08,-3.18672e-12,-32174.7,25.1902], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[21.1348,0.0256494,-8.7265e-06,1.34983e-09,-7.81087e-14,-39041.6,-82.1964], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 180,
    label = "CVCCJCVC",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,D} {3,S} {7,S}
2  C u0 p0 c0 {1,D} {4,S} {6,S}
3  C u0 p0 c0 {1,S} {5,D} {8,S}
4  C u1 p0 c0 {2,S} {9,S} {10,S}
5  C u0 p0 c0 {3,D} {11,S} {12,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.94596,0.0568784,-4.31336e-05,1.6817e-08,-2.67926e-12,23515.7,39.8189], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[14.0879,0.0162399,-5.64769e-06,8.86858e-10,-5.18699e-14,17679.9,-51.3735], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""3/1/95  Z&B""",
    longDesc = 
u"""
3/1/95  Z&B
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
Duplicate of species CVCCVCCJ (i.e. same molecular structure according to RMG)
[CH2]C=CC=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 181,
    label = "NH",
    molecule = 
"""
multiplicity 3
1 N u2 p1 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.49295,0.000311796,-1.48907e-06,2.48167e-09,-1.03571e-12,42106,1.84835], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.78373,0.00132986,-4.24786e-07,7.83494e-11,-5.50451e-15,42346.2,5.74085], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT/A""",
    longDesc = 
u"""
ATcT/A.
[NH]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 182,
    label = "CJVCCVO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {4,S}
2 C u0 p0 c0 {1,S} {5,D} {6,S}
3 C u1 p0 c0 {1,D} {7,S}
4 H u0 p0 c0 {1,S}
5 O u0 p2 c0 {2,D}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.46654,0.032339,-3.05588e-05,1.44082e-08,-2.65601e-12,17885,18.085], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[10.7483,0.00619823,-2.06131e-06,3.14419e-10,-1.8031e-14,15141,-30.1266], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""4/ 8/94 THERM""",
    longDesc = 
u"""
4/ 8/94 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 183,
    label = "HO2CHO",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {4,D} {5,S}
2 O u0 p2 c0 {1,S} {3,S}
3 O u0 p2 c0 {2,S} {6,S}
4 O u0 p2 c0 {1,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.42465,0.0219706,-1.68706e-05,6.25612e-09,-9.11646e-13,-35482.8,17.5028], Tmin=(298,'K'), Tmax=(1378,'K')),
            NASAPolynomial(coeffs=[9.87504,0.00464664,-1.67231e-06,2.68624e-10,-1.59595e-14,-38050.2,-22.4939], Tmin=(1378,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""6/26/95 THERM""",
    longDesc = 
u"""
6/26/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 184,
    label = "C5H11-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {1,S} {2,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.817689,0.0492654,-2.13786e-05,1.85532e-09,7.06259e-13,3262.21,25.8981], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[14.7177,0.0241668,-8.18648e-06,1.26272e-09,-7.29268e-14,-2260.28,-51.2966], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 185,
    label = "C5H11-1",
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
            NASAPolynomial(coeffs=[0.098319,0.0558654,-3.28856e-05,9.58367e-09,-1.08641e-12,4820.66,28.6921], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[15.1918,0.0240339,-8.19718e-06,1.27003e-09,-7.35728e-14,-802.149,-53.6479], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 186,
    label = "CH2OCHO",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {3,S} {4,S} {5,S}
2 C u0 p0 c0 {3,S} {6,D} {7,S}
3 O u0 p2 c0 {1,S} {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 O u0 p2 c0 {2,D}
7 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.31032,0.0180474,-2.7152e-06,-4.60919e-09,1.70037e-12,-20291.1,17.155], Tmin=(298,'K'), Tmax=(1442,'K')),
            NASAPolynomial(coeffs=[10.096,0.00719887,-2.59813e-06,4.18111e-10,-2.48723e-14,-23638.9,-27.1144], Tmin=(1442,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""4/15/ 8 THERM""",
    longDesc = 
u"""
4/15/ 8 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]OC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 187,
    label = "CH3OCO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 O u0 p2 c0 {1,S} {3,S}
3 C u1 p0 c0 {2,S} {7,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 O u0 p2 c0 {3,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.16215,0.0138038,-3.08486e-07,-4.56431e-09,1.4691e-12,-21013,8.64301], Tmin=(298,'K'), Tmax=(1601,'K')),
            NASAPolynomial(coeffs=[9.7366,0.00742433,-2.65642e-06,4.25031e-10,-2.51825e-14,-23601.6,-23.6353], Tmin=(1601,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""5/ 8/ 3 THERM""",
    longDesc = 
u"""
5/ 8/ 3 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CO[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 188,
    label = "NC5KET24O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {4,S} {16,D}
6  O u1 p2 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.40946,0.0475587,-2.19238e-05,2.4325e-09,5.88505e-13,-27653.6,11.7394], Tmin=(298,'K'), Tmax=(1382,'K')),
            NASAPolynomial(coeffs=[18.7901,0.021001,-7.16188e-06,1.11001e-09,-6.43339e-14,-33289.9,-67.8863], Tmin=(1382,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)CC(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 189,
    label = "BC5H10",
    molecule = 
"""
1  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {12,S} {13,S} {14,S}
3  C u0 p0 c0 {5,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {3,S} {4,D} {15,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.604883,0.0496635,-2.66571e-05,6.47312e-09,-4.84018e-13,-8063.12,22.3818], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[14.0426,0.0227915,-7.74903e-06,1.19815e-09,-6.93127e-14,-13216,-51.447], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C(C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 190,
    label = "IC3H5OOCH2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,D} {5,S}
3  C u0 p0 c0 {2,D} {10,S} {11,S}
4  C u1 p0 c0 {6,S} {12,S} {13,S}
5  O u0 p2 c0 {2,S} {6,S}
6  O u0 p2 c0 {4,S} {5,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.26548,-0.00361338,7.31804e-05,-5.49846e-08,1.20447e-11,13310.9,53.3429], Tmin=(298,'K'), Tmax=(1410,'K')),
            NASAPolynomial(coeffs=[5.05336,0.0362904,-1.54013e-05,2.74342e-09,-1.74719e-13,4083.26,-2.5539], Tmin=(1410,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""L 2/00""",
    longDesc = 
u"""
L 2/00
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]OOC(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 191,
    label = "IC5KETDBO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {15,D} {16,S}
6  O u1 p2 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.21028,0.0444018,-1.51327e-05,-2.59537e-09,1.79446e-12,-26464,5.32662], Tmin=(298,'K'), Tmax=(1519,'K')),
            NASAPolynomial(coeffs=[20.0669,0.0202402,-6.97501e-06,1.08902e-09,-6.34525e-14,-32585.9,-78.0518], Tmin=(1519,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)([O])CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 192,
    label = "AC5H10",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {14,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.539429,0.054449,-3.32708e-05,1.03048e-08,-1.28363e-12,-6539.67,29.035], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[14.1931,0.0226551,-7.70009e-06,1.19031e-09,-6.8849e-14,-11949.1,-51.0689], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 193,
    label = "NC5KET13O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,D} {16,S}
6  O u1 p2 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.92095,0.0494039,-1.996e-05,-1.05667e-09,1.70331e-12,-24364.9,13.3503], Tmin=(298,'K'), Tmax=(1323,'K')),
            NASAPolynomial(coeffs=[21.1225,0.019083,-6.52403e-06,1.01368e-09,-5.88793e-14,-31115.2,-82.1812], Tmin=(1323,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC([O])CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 194,
    label = "OCHO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {4,S}
2 O u1 p2 c0 {1,S}
3 O u0 p2 c0 {1,D}
4 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.68826,-0.00414872,2.55066e-05,-2.84474e-08,1.04423e-11,-16986.7,4.28426], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.14394,0.00559739,-1.99794e-06,3.16179e-10,-1.85614e-14,-17246,5.07779], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATCT/A""",
    longDesc = 
u"""
ATCT/A.
[O]C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 195,
    label = "C4H61-3OOH4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
2  C u1 p0 c0 {1,S} {3,S} {9,S}
3  C u0 p0 c0 {2,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {13,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.0032,0.0499894,-3.43101e-05,1.20104e-08,-1.71133e-12,5177.44,24.9842], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[15.5818,0.0169176,-5.84009e-06,9.12401e-10,-5.31699e-14,12.0753,-53.6684], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[CH]COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 196,
    label = "CC5H11O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
6  H u0 p0 c0 {1,S}
7  O u1 p2 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.691614,0.0660697,-4.75286e-05,1.81072e-08,-2.845e-12,-14654.7,23.6234], Tmin=(298,'K'), Tmax=(1409,'K')),
            NASAPolynomial(coeffs=[18.3214,0.0235429,-7.9447e-06,1.22205e-09,-7.04357e-14,-20634,-70.5656], Tmin=(1409,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)C(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 197,
    label = "NC3H7COCH2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {5,S} {13,D}
5  C u1 p0 c0 {4,S} {14,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.14965,0.0509716,-2.97745e-05,7.64275e-09,-5.47754e-13,-11983.9,19.1969], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[17.2611,0.0191829,-6.4168e-06,9.82027e-10,-5.6427e-14,-17483,-63.0398], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/6/15 THERM""",
    longDesc = 
u"""
10/6/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=O)CCC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 198,
    label = "C5H11O-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {12,S} {16,S} {17,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u1 p2 c0 {5,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.61443,0.0584449,-3.26028e-05,7.76525e-09,-4.4874e-13,-12021,22.3621], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[18.7118,0.0236055,-8.05611e-06,1.24901e-09,-7.23989e-14,-18413.3,-71.1769], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 199,
    label = "C5H9O2-4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {15,S}
5  C u0 p0 c0 {3,S} {4,D} {14,S}
6  O u1 p2 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.389,0.0471705,-2.4039e-05,4.99167e-09,-2.03243e-13,115.233,12.5184], Tmin=(298,'K'), Tmax=(1376,'K')),
            NASAPolynomial(coeffs=[17.3707,0.0209045,-7.30724e-06,1.15107e-09,-6.74601e-14,-5480.2,-65.0109], Tmin=(1376,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CC(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 200,
    label = "PC4H9O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {9,S} {13,S} {14,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  O u1 p2 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.84659,0.0461054,-2.44857e-05,5.11293e-09,-1.07538e-13,-9092.07,19.5237], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[15.3372,0.019279,-6.56857e-06,1.01724e-09,-5.89183e-14,-14195.9,-54.5072], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 201,
    label = "C5H9O24-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  C u1 p0 c0 {3,S} {14,S} {15,S}
6  O u0 p2 c0 {1,S} {3,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-6.45049,0.0825685,-7.30636e-05,3.23768e-08,-5.7248e-12,4474.17,57.4525], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[14.0158,0.0251554,-1.11149e-05,2.14885e-09,-1.4456e-13,-1547.69,-48.8521], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C1CC(C)O1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 202,
    label = "IC5KETDAO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {9,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {15,D} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  O u1 p2 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.77855,0.0549043,-3.12957e-05,7.3334e-09,-3.357e-13,-22810.3,18.95], Tmin=(298,'K'), Tmax=(1485,'K')),
            NASAPolynomial(coeffs=[19.5836,0.0202789,-6.90313e-06,1.0687e-09,-6.18931e-14,-28993.5,-72.7515], Tmin=(1485,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C[O])CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 203,
    label = "C4H6-1",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,T}
4  C u0 p0 c0 {3,T} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.42819,0.0249822,6.27371e-06,-2.61748e-08,1.26585e-11,18024.9,13.6684], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.81179,0.0179734,-6.61044e-06,1.05501e-09,-6.19297e-14,16177,-15.9658], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L10/93""",
    longDesc = 
u"""
L10/93.
C#CCC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 204,
    label = "NC3H7CHO",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,D} {13,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u0 p2 c0 {4,D}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.24209,0.0421278,-2.13832e-05,4.22615e-09,-1.03711e-13,-27104.9,22.1567], Tmin=(298,'K'), Tmax=(1679,'K')),
            NASAPolynomial(coeffs=[11.9789,0.0204894,-7.24832e-06,1.15562e-09,-6.8412e-14,-30927.2,-36.393], Tmin=(1679,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 205,
    label = "C4H72-1,3OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {5,S} {9,S}
2  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {2,S} {15,S}
5  O u0 p2 c0 {1,S} {8,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {16,S}
8  O u0 p2 c0 {5,S} {17,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.85542,0.0605524,-4.18122e-05,1.4668e-08,-2.07882e-12,-16035.9,18.7735], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[21.4626,0.0202946,-6.97889e-06,1.08751e-09,-6.32605e-14,-22219.7,-76.0672], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC([CH]COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 206,
    label = "OCH2OCHO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {3,S} {7,D} {8,S}
3 O u0 p2 c0 {1,S} {2,S}
4 O u1 p2 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 O u0 p2 c0 {2,D}
8 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.8954,0.0274119,-1.36476e-05,1.26326e-09,5.1797e-13,-42787.9,20.2333], Tmin=(298,'K'), Tmax=(1523,'K')),
            NASAPolynomial(coeffs=[12.4013,0.00783738,-2.82993e-06,4.55559e-10,-2.71061e-14,-46845.3,-37.8085], Tmin=(1523,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""5/29/14 THERM""",
    longDesc = 
u"""
5/29/14 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]COC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 207,
    label = "O2CH2CHO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,D} {7,S}
3 O u0 p2 c0 {1,S} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 O u0 p2 c0 {2,D}
7 H u0 p0 c0 {2,S}
8 O u1 p2 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.29466,0.0444936,-4.26577e-05,2.07392e-08,-3.96829e-12,-11827.6,36.0779], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[11.1808,0.00914479,-3.1509e-06,4.91944e-10,-2.86639e-14,-15579,-28.7893], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""BOZ_03""",
    longDesc = 
u"""
BOZ_03
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 208,
    label = "NC5KET14O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,D} {16,S}
6  O u1 p2 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.92095,0.0494039,-1.996e-05,-1.05667e-09,1.70331e-12,-24364.9,13.3503], Tmin=(298,'K'), Tmax=(1323,'K')),
            NASAPolynomial(coeffs=[21.1225,0.019083,-6.52403e-06,1.01368e-09,-5.88793e-14,-31115.2,-82.1812], Tmin=(1323,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC([O])CCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 209,
    label = "H15DE25DM",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {6,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {3,S} {7,D}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  C u0 p0 c0 {6,D} {21,S} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.71853,0.0882614,-6.03141e-05,2.15862e-08,-3.19691e-12,-1952.56,38.6049], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[22.5356,0.0323956,-1.10271e-05,1.70641e-09,-9.87758e-14,-10480.9,-91.9785], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)CCC(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 210,
    label = "C4H72-1OOCH3",
    molecule = 
"""
1  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
2  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
3  C u0 p0 c0 {7,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {5,D} {17,S}
5  C u0 p0 c0 {2,S} {4,D} {16,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {3,S} {6,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.44675,0.062143,-4.39326e-05,1.72567e-08,-2.94951e-12,-11143.2,19.4029], Tmin=(298,'K'), Tmax=(1382,'K')),
            NASAPolynomial(coeffs=[19.2915,0.0234475,-8.16714e-06,1.28338e-09,-7.50838e-14,-17268.9,-71.651], Tmin=(1382,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CCOOC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 211,
    label = "IC5KETADO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {9,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,D} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  O u1 p2 c0 {4,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[6.20341,0.0378489,-8.21632e-06,-4.45999e-09,1.69104e-12,-23370.7,3.88282], Tmin=(298,'K'), Tmax=(1670,'K')),
            NASAPolynomial(coeffs=[15.4491,0.025883,-9.38195e-06,1.5192e-09,-9.08764e-14,-27312.9,-48.9349], Tmin=(1670,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C=O)CC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 212,
    label = "B13DE2M",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u0 p0 c0 {2,S} {5,D} {9,S}
4  C u0 p0 c0 {2,D} {12,S} {13,S}
5  C u0 p0 c0 {3,D} {10,S} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.47217,0.0559413,-4.50459e-05,1.96181e-08,-3.5201e-12,7149.6,25.6277], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[14.0758,0.0182376,-6.20748e-06,9.60351e-10,-5.55742e-14,2400.34,-51.2988], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 213,
    label = "CH3OCH2O",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2 C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3 O u0 p2 c0 {1,S} {2,S}
4 O u1 p2 c0 {2,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.63414,0.0089283,1.37226e-05,-1.40497e-08,3.54626e-12,-22282.5,1.93589], Tmin=(298,'K'), Tmax=(1523,'K')),
            NASAPolynomial(coeffs=[9.81289,0.0121313,-4.30286e-06,6.84443e-10,-4.03863e-14,-25076.1,-25.1866], Tmin=(1523,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""5/15/14 THERM""",
    longDesc = 
u"""
5/15/14 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
COC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 214,
    label = "BC5H11O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
6  O u0 p2 c0 {1,S} {18,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.53424,0.0694505,-5.16918e-05,2.0869e-08,-3.52553e-12,-18280.7,20.372], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[19.7546,0.0247325,-8.41612e-06,1.30187e-09,-7.53317e-14,-24467.2,-76.8414], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)(C)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 215,
    label = "AC5H9-C",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u1 p0 c0 {2,S} {3,S} {12,S}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.841373,0.0527157,-3.27335e-05,1.01974e-08,-1.26085e-12,10432.1,29.3317], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[13.7918,0.0209645,-7.13793e-06,1.10481e-09,-6.39628e-14,5095.56,-50.1391], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)[CH]C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 216,
    label = "AC5H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {1,S} {15,S} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
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
            NASAPolynomial(coeffs=[-1.22778,0.0632926,-4.52076e-05,1.75515e-08,-2.87257e-12,4325.62,33.8598], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[15.2292,0.0238483,-8.09849e-06,1.25097e-09,-7.23134e-14,-1358.28,-54.2832], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 217,
    label = "BC5H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {1,S} {3,S} {4,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.65843,0.0460868,-1.82957e-05,8.82319e-10,7.56538e-13,1290.55,21.2275], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[14.1465,0.024834,-8.45503e-06,1.30849e-09,-7.57429e-14,-3880.02,-48.7193], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[C](C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 218,
    label = "IC4H9",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {1,S} {12,S} {13,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.120802,0.0473187,-3.1644e-05,1.1423e-08,-1.74785e-12,6840.33,24.4291], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[12.3262,0.0192058,-6.52064e-06,1.00704e-09,-5.82039e-14,2509.96,-41.3479], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 219,
    label = "TC4H9",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.05842,0.0341134,-9.03157e-06,-2.95313e-09,1.41437e-12,4226.99,22.3965], Tmin=(298,'K'), Tmax=(1380,'K')),
            NASAPolynomial(coeffs=[10.2683,0.0209965,-7.14946e-06,1.10648e-09,-6.40498e-14,157.543,-30.0961], Tmin=(1380,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C](C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 220,
    label = "C4H3-I",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {5,S} {6,S}
2 C u1 p0 c0 {1,D} {3,S}
3 C u0 p0 c0 {2,S} {4,T}
4 C u0 p0 c0 {3,T} {7,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.08304,0.0408343,-6.21597e-05,5.16794e-08,-1.70292e-11,58005.1,13.6175], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[9.09782,0.00922071,-3.38784e-06,4.91605e-10,-1.45298e-14,56600.6,-19.8026], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""AB1/93""",
    longDesc = 
u"""
AB1/93
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#C[C]=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 221,
    label = "C5H11O-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  O u1 p2 c0 {3,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.80798,0.0595687,-3.62937e-05,1.05786e-08,-1.11087e-12,-14129.5,20.1085], Tmin=(298,'K'), Tmax=(1407,'K')),
            NASAPolynomial(coeffs=[18.5178,0.023356,-7.87897e-06,1.21191e-09,-6.98573e-14,-20148.9,-70.5224], Tmin=(1407,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 222,
    label = "CC5H9-B",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
3  C u1 p0 c0 {1,S} {2,S} {4,S}
4  C u0 p0 c0 {3,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.691502,0.0509805,-2.96448e-05,8.28273e-09,-8.57575e-13,10747.8,27.6712], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[13.525,0.0214364,-7.35304e-06,1.14369e-09,-6.64364e-14,5413.82,-50.0302], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[C](C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 223,
    label = "BC5H11O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
6  O u1 p2 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.67366,0.0517139,-2.49458e-05,4.54085e-09,-3.31185e-14,-16322.8,9.40461], Tmin=(298,'K'), Tmax=(1378,'K')),
            NASAPolynomial(coeffs=[18.5549,0.0244787,-8.51927e-06,1.33803e-09,-7.82547e-14,-22366.4,-73.3946], Tmin=(1378,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 224,
    label = "CH2OCH2CHO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,D} {10,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  O u1 p2 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  O u0 p2 c0 {3,D}
10 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.09797,0.0250971,-1.0438e-05,6.95156e-10,3.93229e-13,-15120.3,11.6942], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[11.1193,0.0127913,-4.35464e-06,6.73922e-10,-3.90119e-14,-17969.4,-27.4688], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]CCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 225,
    label = "C4H8OOH1-3O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {3,S} {7,S}
6  O u0 p2 c0 {1,S} {16,S}
7  O u0 p2 c0 {5,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u1 p2 c0 {6,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.02241,0.0653863,-4.89646e-05,1.90438e-08,-3.02309e-12,-22580.6,20.5729], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[21.5735,0.0204529,-6.99498e-06,1.08598e-09,-6.30071e-14,-28842.8,-78.4561], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CCOO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 226,
    label = "PC4H9O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  O u0 p2 c0 {3,S} {15,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u1 p2 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.80541,0.0526493,-3.3087e-05,1.04593e-08,-1.32306e-12,-10386.7,22.483], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[16.612,0.0204752,-7.01415e-06,1.09011e-09,-6.32914e-14,-15790.2,-57.9272], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 227,
    label = "HO2CH2CO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 O u0 p2 c0 {1,S} {4,S}
3 C u1 p0 c0 {1,S} {7,D}
4 O u0 p2 c0 {2,S} {8,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 O u0 p2 c0 {3,D}
8 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.22682,0.0356781,-3.26402e-05,1.47652e-08,-2.64794e-12,-11873.5,19.1581], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[10.4146,0.011268,-5.17495e-06,1.00333e-09,-6.68166e-14,-14095.6,-22.7894], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""BOZ_03""",
    longDesc = 
u"""
BOZ_03
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=[C]COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 228,
    label = "HOCH2OCO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p2 c0 {1,S} {7,S}
4 C u1 p0 c0 {2,S} {8,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
8 O u0 p2 c0 {4,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.95255,0.00842196,1.36742e-05,-1.46786e-08,3.84144e-12,-44447,2.85657], Tmin=(298,'K'), Tmax=(1443,'K')),
            NASAPolynomial(coeffs=[11.1498,0.00934737,-3.35542e-06,5.38037e-10,-3.1926e-14,-47501.2,-29.5984], Tmin=(1443,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""5/29/14 THERM""",
    longDesc = 
u"""
5/29/14 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=[C]OCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 229,
    label = "BC5H11O2H",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {19,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.748348,0.0754437,-5.69961e-05,2.32342e-08,-3.94363e-12,-34963.9,23.5201], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[20.7291,0.0260359,-8.86783e-06,1.37263e-09,-7.94629e-14,-41699.2,-82.9349], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 230,
    label = "C4H8OOH2-4O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {3,S} {16,S}
7  O u0 p2 c0 {5,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u1 p2 c0 {6,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.02241,0.0653863,-4.89646e-05,1.90438e-08,-3.02309e-12,-22580.6,20.5729], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[21.5735,0.0204529,-6.99498e-06,1.08598e-09,-6.30071e-14,-28842.8,-78.4561], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CCO[O])OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 231,
    label = "AC5H9O-C",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {14,S} {15,S}
6  O u1 p2 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.04951,0.0536585,-3.35087e-05,1.06174e-08,-1.39224e-12,-385.113,19.1077], Tmin=(298,'K'), Tmax=(1382,'K')),
            NASAPolynomial(coeffs=[17.679,0.0205915,-7.18764e-06,1.13116e-09,-6.62505e-14,-6251.47,-66.2132], Tmin=(1382,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)C(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 232,
    label = "IC4H9O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
5  H u0 p0 c0 {1,S}
6  O u1 p2 c0 {4,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.80297,0.0156874,6.81105e-05,-9.83347e-08,3.95262e-11,-10083.2,9.78963], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[11.631,0.0247982,-9.01551e-06,1.46715e-09,-8.83215e-14,-13785.5,-38.1956], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""A08/04""",
    longDesc = 
u"""
A08/04.
CC(C)C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 233,
    label = "C5H11O2-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {4,S} {18,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.4785,0.0656223,-4.24002e-05,1.39456e-08,-1.86036e-12,-13306.5,25.7185], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[19.9509,0.0248182,-8.50415e-06,1.32191e-09,-7.67597e-14,-19978.3,-74.3587], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCCO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 234,
    label = "CC5H9O-B",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {5,D} {13,S}
5  C u0 p0 c0 {4,D} {14,S} {15,S}
6  O u1 p2 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.75615,0.0520514,-3.10967e-05,8.25895e-09,-6.59917e-13,-1351.63,14.099], Tmin=(298,'K'), Tmax=(1377,'K')),
            NASAPolynomial(coeffs=[18.6974,0.0184543,-6.18683e-06,9.48851e-10,-5.46204e-14,-7158.98,-72.654], Tmin=(1377,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 235,
    label = "TC4H9O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  O u1 p2 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.77057,0.0268033,4.12718e-05,-7.22055e-08,3.02642e-11,-12707.9,12.1533], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[12.7372,0.0233707,-8.50517e-06,1.3852e-09,-8.34398e-14,-16694,-45.3156], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T08/04""",
    longDesc = 
u"""
T08/04.
CC(C)(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 236,
    label = "C5H6",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,D} {8,S}
3  C u0 p0 c0 {1,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {5,S} {10,S}
5  C u0 p0 c0 {2,D} {4,S} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.86109,0.014804,7.21089e-05,-1.13381e-07,4.869e-11,14801.8,21.3535], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[9.97578,0.0189055,-6.84115e-06,1.10993e-09,-6.66802e-14,11081.7,-32.2095], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 1/90""",
    longDesc = 
u"""
T 1/90.
C1C=CCC=1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 237,
    label = "C4H71-3OOCH3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {7,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {5,D} {15,S}
5  C u0 p0 c0 {4,D} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {3,S} {6,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.542405,0.0751511,-6.00685e-05,2.49701e-08,-4.2232e-12,-11306.1,30.2982], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[21.2299,0.0207004,-6.99464e-06,1.07837e-09,-6.23e-14,-18500.6,-85.35], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)OOC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 238,
    label = "AC5H11O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {9,S} {16,S} {17,S}
5  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  O u1 p2 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.843925,0.0629724,-4.12079e-05,1.37396e-08,-1.84098e-12,-12581.9,25.0754], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[18.2339,0.0237924,-8.07056e-06,1.24593e-09,-7.19996e-14,-18745,-68.7959], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 239,
    label = "C2H3OOH",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,S} {5,S}
2 C u0 p0 c0 {1,D} {6,S} {7,S}
3 O u0 p2 c0 {1,S} {4,S}
4 O u0 p2 c0 {3,S} {8,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.35644,0.0337002,-2.75989e-05,1.14223e-08,-1.89489e-12,-5499.97,19.8354], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[11.575,0.00809909,-2.81809e-06,4.42698e-10,-2.58998e-14,-8848.53,-34.3859], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""4/18/ 8 THERM""",
    longDesc = 
u"""
4/18/ 8 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 240,
    label = "C4H7O1-2OOH-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {3,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.944964,0.0731027,-6.71228e-05,3.12252e-08,-5.67207e-12,-27451.1,32.5427], Tmin=(298,'K'), Tmax=(1435,'K')),
            NASAPolynomial(coeffs=[18.3476,0.0172628,-5.5344e-06,8.21396e-10,-4.6148e-14,-32922.4,-66.992], Tmin=(1435,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(OO)C1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 241,
    label = "CH3COCH2O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {10,D}
4  O u1 p2 c0 {2,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 O u0 p2 c0 {3,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.8596,0.0178955,7.41506e-07,-5.40033e-09,1.47393e-12,-19071.5,2.70988], Tmin=(298,'K'), Tmax=(2002,'K')),
            NASAPolynomial(coeffs=[9.84062,0.0159181,-5.85165e-06,9.5616e-10,-5.75477e-14,-21121.5,-21.2331], Tmin=(2002,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/ 8/13 THERM""",
    longDesc = 
u"""
2/ 8/13 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 242,
    label = "C5H11O-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  O u1 p2 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.80798,0.0595687,-3.62937e-05,1.05786e-08,-1.11087e-12,-14129.5,18.7296], Tmin=(298,'K'), Tmax=(1407,'K')),
            NASAPolynomial(coeffs=[18.5178,0.023356,-7.87897e-06,1.21191e-09,-6.98573e-14,-20148.9,-71.9014], Tmin=(1407,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC([O])CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 243,
    label = "C5H91-2,3OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
3  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {2,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {2,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.91226,0.0795651,-6.1323e-05,2.44461e-08,-3.95034e-12,-19651.7,23.2188], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[25.5056,0.023699,-8.07033e-06,1.24935e-09,-7.23427e-14,-27135.3,-96.9349], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(OO)C(CC)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 244,
    label = "C5H91-4,5OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
5  C u1 p0 c0 {3,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.50133,0.0793141,-5.51609e-05,1.83964e-08,-2.32881e-12,-17513,26.433], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[27.5716,0.0226282,-7.85652e-06,1.23252e-09,-7.20427e-14,-26270.3,-108.648], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC(COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 245,
    label = "SC2H2OH",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {4,S} {5,S}
2 C u1 p0 c0 {1,D} {3,S}
3 O u0 p2 c0 {2,S} {6,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.63792,0.0264969,-2.74821e-05,1.43111e-08,-2.85795e-12,11146.7,15.8715], Tmin=(298,'K'), Tmax=(1410,'K')),
            NASAPolynomial(coeffs=[7.99235,0.00583109,-1.89243e-06,2.83129e-10,-1.59933e-14,9512.37,-16.2058], Tmin=(1410,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/25/15""",
    longDesc = 
u"""
9/25/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=[C]O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 246,
    label = "CH3CHCH2COCH3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {2,S} {14,S}
5  C u0 p0 c0 {1,S} {3,S} {15,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.82775,0.0359153,-6.78571e-06,-6.2482e-09,2.43917e-12,-10849.6,15.5399], Tmin=(298,'K'), Tmax=(1313,'K')),
            NASAPolynomial(coeffs=[14.5648,0.0210277,-6.95711e-06,1.05654e-09,-6.03738e-14,-15461.8,-45.5118], Tmin=(1313,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/6/15""",
    longDesc = 
u"""
10/6/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CC(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 247,
    label = "C5H9O1-3OOH-2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {4,S}
7  O u0 p2 c0 {2,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.67709,0.0925192,-7.56246e-05,3.08656e-08,-4.96975e-12,-30517.9,51.9808], Tmin=(298,'K'), Tmax=(1421,'K')),
            NASAPolynomial(coeffs=[22.9239,0.0222245,-7.5093e-06,1.15689e-09,-6.67794e-14,-39314.9,-93.8427], Tmin=(1421,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1OCC1OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 248,
    label = "IC3H5CHO",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,D} {4,S}
3  C u0 p0 c0 {2,D} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {8,D} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  O u0 p2 c0 {4,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.09373,0.0443315,-3.41918e-05,1.3937e-08,-2.33791e-12,-15674.6,19.4458], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[13.3892,0.0139115,-4.75821e-06,7.38737e-10,-4.28607e-14,-19791.7,-46.0146], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 249,
    label = "C6H5OO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,D} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {1,D} {6,S} {12,S}
4  C u0 p0 c0 {2,D} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {3,S} {5,D} {11,S}
7  O u0 p2 c0 {1,S} {13,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {3,S}
13 O u1 p2 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.99165,0.0703857,-6.34401e-05,2.91549e-08,-5.30707e-12,14132,42.0143], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[16.7078,0.0162326,-5.4797e-06,8.4351e-10,-4.86562e-14,8142.43,-60.8347], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""3/26/ 9 THERM""",
    longDesc = 
u"""
3/26/ 9 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OC1C=CC=CC=1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 250,
    label = "CC5H10OOH-B",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {1,S} {3,S} {4,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.95783,0.0596256,-3.9221e-05,1.37409e-08,-2.01546e-12,-11335.8,13.3114], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[19.382,0.0243097,-8.20157e-06,1.26129e-09,-7.26837e-14,-16796.9,-69.8331], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C](C)C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 251,
    label = "BC5H10OOH-A",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {1,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.67959,0.0703877,-5.38977e-05,2.22635e-08,-3.82053e-12,-10382,22.0634], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[20.3178,0.0239006,-8.14494e-06,1.26119e-09,-7.30299e-14,-16627.1,-77.1009], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)(CC)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 252,
    label = "C5H9O1-5OOH-2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {6,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
6  O u0 p2 c0 {4,S} {5,S}
7  O u0 p2 c0 {1,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-9.63082,0.10654,-8.95954e-05,3.77712e-08,-6.32249e-12,-39749.9,71.8166], Tmin=(298,'K'), Tmax=(1407,'K')),
            NASAPolynomial(coeffs=[22.4619,0.0242252,-8.36381e-06,1.30713e-09,-7.6201e-14,-50010.2,-97.6836], Tmin=(1407,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOC1CCCOC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 253,
    label = "H15DE2M-T",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
3  C u0 p0 c0 {4,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {15,S} {16,S}
6  C u0 p0 c0 {7,D} {17,S} {18,S}
7  C u1 p0 c0 {2,S} {6,D}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.388856,0.0670216,-4.29671e-05,1.4249e-08,-1.96106e-12,30362.7,29.5444], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[19.0145,0.0261932,-9.01184e-06,1.40456e-09,-8.17075e-14,23555.1,-71.5569], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=[C]CCC(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 254,
    label = "C4H5-2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2 C u1 p0 c0 {4,S} {8,S} {9,S}
3 C u0 p0 c0 {1,S} {4,T}
4 C u0 p0 c0 {2,S} {3,T}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.31011,0.0283747,-1.63837e-05,4.46252e-09,-4.30511e-13,35584.3,14.9106], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[10.3231,0.0117626,-4.00005e-06,6.18728e-10,-3.58084e-14,32586.1,-28.8794], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""H6W/94""",
    longDesc = 
u"""
H6W/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C#CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 255,
    label = "CH2CCH2OH",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {3,D} {7,S} {8,S}
3 C u1 p0 c0 {1,S} {2,D}
4 O u0 p2 c0 {1,S} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.88423,0.0242428,-1.14152e-05,1.71775e-09,1.42177e-13,11793.6,15.2102], Tmin=(298,'K'), Tmax=(1372,'K')),
            NASAPolynomial(coeffs=[9.70702,0.0113973,-3.77994e-06,5.75209e-10,-3.29229e-14,9132.13,-22.5013], Tmin=(1372,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/95 THERM""",
    longDesc = 
u"""
9/ 8/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=[C]CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 256,
    label = "C3H6O1-2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
4  O u0 p2 c0 {1,S} {2,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.42807,0.00625177,6.13196e-05,-8.60387e-08,3.51371e-11,-12844.7,10.4245], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.01491,0.017392,-6.26028e-06,1.01188e-09,-6.06239e-14,-15198.1,-18.828], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""A01/05""",
    longDesc = 
u"""
A01/05.
CC1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 257,
    label = "C5H11O2H-1",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {4,S} {7,S}
7  O u0 p2 c0 {6,S} {19,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.756967,0.0713994,-4.75006e-05,1.62337e-08,-2.26834e-12,-30001,28.5559], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[20.8827,0.0261587,-8.9686e-06,1.39463e-09,-8.10033e-14,-37186.4,-80.1917], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 258,
    label = "CVCCVCCOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,D} {9,S}
3  C u0 p0 c0 {2,D} {4,S} {11,S}
4  C u0 p0 c0 {3,S} {5,D} {10,S}
5  C u0 p0 c0 {4,D} {12,S} {13,S}
6  O u0 p2 c0 {1,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.531488,0.0606984,-4.815e-05,2.00308e-08,-3.38987e-12,-10330.1,30.7961], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[16.308,0.0179958,-6.03116e-06,9.23992e-10,-5.31254e-14,-15820.5,-58.4137], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/23/ 9 WKM""",
    longDesc = 
u"""
1/23/ 9 WKM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC=CCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 259,
    label = "C5H10OOH1-5",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u1 p0 c0 {3,S} {16,S} {17,S}
6  O u0 p2 c0 {4,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.960306,0.0685338,-4.72581e-05,1.69913e-08,-2.53056e-12,-5252.27,29.9006], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[20.2627,0.0243041,-8.36491e-06,1.30409e-09,-7.58782e-14,-12079.2,-74.1333], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 260,
    label = "HOCH2O",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u1 p2 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.11183,0.00753851,3.77337e-06,-5.38746e-09,1.45616e-12,-22802.3,7.46807], Tmin=(298,'K'), Tmax=(1452,'K')),
            NASAPolynomial(coeffs=[6.39522,0.00743673,-2.50422e-06,3.8488e-10,-2.21779e-14,-24110.9,-6.63866], Tmin=(1452,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/16/99 THERM""",
    longDesc = 
u"""
2/16/99 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 261,
    label = "C3H6OH1-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {4,S} {10,S}
4  O u0 p2 c0 {3,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.13615,0.0368851,-2.24074e-05,6.62399e-09,-7.32206e-13,-10927.3,22.4919], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[11.4795,0.0145881,-4.88359e-06,7.47607e-10,-4.29608e-14,-14683.2,-33.6826], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 2/15""",
    longDesc = 
u"""
10/ 2/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 262,
    label = "SC3H5OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,D} {8,S}
3  C u0 p0 c0 {2,D} {4,S} {9,S}
4  O u0 p2 c0 {3,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.0353977,0.0434969,-3.7448e-05,1.70906e-08,-3.13775e-12,-20250.3,24.1528], Tmin=(298,'K'), Tmax=(1404,'K')),
            NASAPolynomial(coeffs=[11.1222,0.0127745,-4.25316e-06,6.48216e-10,-3.71191e-14,-23669.1,-34.1335], Tmin=(1404,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/ 3/ 9""",
    longDesc = 
u"""
2/ 3/ 9
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 263,
    label = "C5H91-3,4OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {3,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {8,S}
7  O u0 p2 c0 {2,S} {9,S}
8  O u0 p2 c0 {6,S} {19,S}
9  O u0 p2 c0 {7,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.48775,0.0815272,-6.44013e-05,2.63444e-08,-4.35804e-12,-19563.5,25.0627], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[25.4681,0.0236307,-8.02474e-06,1.23995e-09,-7.17032e-14,-27063.6,-96.7606], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(OO)C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 264,
    label = "C2H4O2H",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {7,S} {8,S}
3 O u0 p2 c0 {1,S} {4,S}
4 O u0 p2 c0 {3,S} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.75788,0.0288272,-2.08302e-05,8.47401e-09,-1.48618e-12,3001.54,15.9922], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[10.0591,0.0113379,-3.89403e-06,6.06091e-10,-3.52212e-14,424.049,-23.2087], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 265,
    label = "NH2",
    molecule = 
"""
multiplicity 2
1 N u1 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.19198,-0.00204603,6.67756e-06,-5.24907e-09,1.5559e-12,21499.1,-0.0904785], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.59263,0.00347684,-1.08272e-06,1.49343e-10,-5.75241e-15,21886.5,7.90565], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""AMIDOGEN RAD IU3/03""",
    longDesc = 
u"""
AMIDOGEN RAD IU3/03.
[NH2]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 266,
    label = "C4H72-1,4OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
4  C u1 p0 c0 {1,S} {3,S} {15,S}
5  O u0 p2 c0 {2,S} {8,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {16,S}
8  O u0 p2 c0 {5,S} {17,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.41715,0.0589886,-3.76843e-05,1.1845e-08,-1.46347e-12,-13843.8,22.7095], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[21.0229,0.0210128,-7.3035e-06,1.14623e-09,-6.70076e-14,-20290.8,-72.9927], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOC[CH]CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 267,
    label = "IC4H8",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.0572478,0.0417769,-2.49096e-05,7.54294e-09,-9.23202e-13,-3721.66,23.5699], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[11.1444,0.0181609,-6.17791e-06,9.55482e-10,-5.52826e-14,-7840.25,-36.8509], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 268,
    label = "C5H5",
    molecule = 
"""
multiplicity 2
1  C u1 p0 c0 {2,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,D} {7,S}
3  C u0 p0 c0 {2,D} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {5,D} {9,S}
5  C u0 p0 c0 {1,S} {4,D} {10,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.97555,0.0741371,-0.000111803,9.04629e-08,-2.81e-11,30176.9,36.7154], Tmin=(298,'K'), Tmax=(969.35,'K')),
            NASAPolynomial(coeffs=[1.33676,0.0324794,-1.67588e-05,4.03514e-09,-3.70739e-13,30073.1,16.0316], Tmin=(969.35,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""TAK0505""",
    longDesc = 
u"""
TAK0505
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
[CH]1C=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 269,
    label = "IC4H7O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  O u1 p2 c0 {2,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.74701,0.0407783,-2.4475e-05,7.06503e-09,-7.51571e-13,4869.79,19.4536], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[13.3458,0.0161219,-5.44376e-06,8.38199e-10,-4.83608e-14,611.444,-43.6819], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""4/ 3/ 0 THERM""",
    longDesc = 
u"""
4/ 3/ 0 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 270,
    label = "IC3H7CHO",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,D} {13,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u0 p2 c0 {4,D}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.273021,0.0489696,-3.1277e-05,1.00053e-08,-1.27512e-12,-27605.5,28.3451], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[13.7502,0.0183127,-6.28573e-06,9.78251e-10,-5.68539e-14,-32693.7,-47.7271], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/22/96 THERM""",
    longDesc = 
u"""
2/22/96 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 271,
    label = "CC5H10",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {5,D} {13,S}
5  C u0 p0 c0 {4,D} {14,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.964397,0.0568483,-3.66324e-05,1.23009e-08,-1.71393e-12,-5935.02,30.2998], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[14.4243,0.0226455,-7.73646e-06,1.2e-09,-6.95711e-14,-11498.1,-53.0391], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 272,
    label = "PC4H9CO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {3,S} {15,D}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.32674,0.0469568,-2.17003e-05,2.33167e-09,6.15109e-13,-11582.4,19.9866], Tmin=(298,'K'), Tmax=(1380,'K')),
            NASAPolynomial(coeffs=[16.7761,0.0202898,-6.94418e-06,1.07888e-09,-6.2636e-14,-17242.7,-60.0185], Tmin=(1380,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 273,
    label = "B13DE2MJ",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,D}
2  C u0 p0 c0 {1,S} {5,D} {6,S}
3  C u1 p0 c0 {1,S} {7,S} {8,S}
4  C u0 p0 c0 {1,D} {9,S} {10,S}
5  C u0 p0 c0 {2,D} {11,S} {12,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.898037,0.0565174,-4.76108e-05,2.09672e-08,-3.73076e-12,25288.2,27.4967], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[14.8413,0.0152426,-5.1356e-06,7.89787e-10,-4.55353e-14,20266.2,-55.4475], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=C)C=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 274,
    label = "H15DE25DM-A",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {3,S} {6,D}
5  C u0 p0 c0 {2,S} {7,S} {8,D}
6  C u0 p0 c0 {4,D} {20,S} {21,S}
7  C u1 p0 c0 {5,S} {16,S} {17,S}
8  C u0 p0 c0 {5,D} {18,S} {19,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.21863,0.0891653,-6.33361e-05,2.32007e-08,-3.4633e-12,16197.1,41.5117], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[23.3283,0.0293694,-9.94259e-06,1.53368e-09,-8.8603e-14,7373.19,-95.5902], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=C)CCC(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 275,
    label = "C5H7",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  C u1 p0 c0 {2,S} {5,S} {10,S}
4  C u0 p0 c0 {1,S} {5,D} {11,S}
5  C u0 p0 c0 {3,S} {4,D} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-6.75118,0.0606462,-4.0126e-05,1.22052e-08,-1.3346e-12,20136.5,56.2695], Tmin=(298,'K'), Tmax=(1377,'K')),
            NASAPolynomial(coeffs=[13.663,0.0168061,-5.98747e-06,9.55341e-10,-5.64952e-14,12723.9,-54.6331], Tmin=(1377,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/22/ 9 WKM""",
    longDesc = 
u"""
1/22/ 9 WKM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]1C=CCC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 276,
    label = "PC4H9O2H",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {3,S} {6,S}
6  O u0 p2 c0 {5,S} {16,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.04178,0.0585997,-3.84212e-05,1.28728e-08,-1.75491e-12,-27074.4,25.518], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[17.5611,0.0217833,-7.46366e-06,1.16012e-09,-6.7363e-14,-33003.6,-63.8547], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 277,
    label = "AC5H9-A2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u1 p0 c0 {3,S} {11,S} {12,S}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.954048,0.0550153,-3.58308e-05,1.16444e-08,-1.49087e-12,11596,30.8477], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[15.002,0.0195965,-6.60206e-06,1.01536e-09,-5.85455e-14,5900.81,-55.4529], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=C)CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 278,
    label = "CC5H9-A",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,D} {10,S}
4  C u1 p0 c0 {1,S} {11,S} {12,S}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.67891,0.0537496,-3.61117e-05,1.28563e-08,-1.92156e-12,18797.4,31.2403], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[13.8636,0.0206363,-7.05717e-06,1.0954e-09,-6.35382e-14,13610.4,-47.2507], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)C=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 279,
    label = "IC3H4CHO-A",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,D} {4,S}
2  C u1 p0 c0 {1,S} {5,S} {6,S}
3  C u0 p0 c0 {1,D} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {7,D} {10,S}
5  H u0 p0 c0 {2,S}
6  H u0 p0 c0 {2,S}
7  O u0 p2 c0 {4,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.764345,0.0445242,-3.61034e-05,1.48295e-08,-2.43809e-12,2447.33,20.8542], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[14.1737,0.0109162,-3.69021e-06,5.69228e-10,-3.29023e-14,-1928.68,-50.2664], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=C)C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 280,
    label = "TC3H6CHO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u1 p0 c0 {1,S} {2,S} {4,S}
4  C u0 p0 c0 {3,S} {11,D} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 O u0 p2 c0 {4,D}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.87053,0.041487,-2.66816e-05,9.01532e-09,-1.27871e-12,-8977.31,16.6174], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[13.1013,0.0166392,-5.68458e-06,8.81808e-10,-5.1129e-14,-13063.9,-44.2706], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/22/96 THERM""",
    longDesc = 
u"""
2/22/96 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C](C)C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 281,
    label = "IC4H7OH",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  O u0 p2 c0 {1,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.04124,0.0414387,-2.55229e-05,8.28133e-09,-1.10654e-12,-22271,18.8699], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[12.3304,0.0183885,-6.06722e-06,9.19055e-10,-5.24036e-14,-25945.2,-36.7418], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 282,
    label = "IC4H9O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  O u0 p2 c0 {2,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u1 p2 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.7722,0.0534033,-3.31042e-05,9.24466e-09,-8.01707e-13,-11776.9,20.9581], Tmin=(298,'K'), Tmax=(1432,'K')),
            NASAPolynomial(coeffs=[17.8794,0.0182475,-6.01252e-06,9.11107e-10,-5.20019e-14,-17457,-66.1553], Tmin=(1432,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 283,
    label = "TC4H9O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {15,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u1 p2 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.63892,0.0544717,-3.75505e-05,1.40479e-08,-2.27969e-12,-14759.9,14.0326], Tmin=(298,'K'), Tmax=(1380,'K')),
            NASAPolynomial(coeffs=[18.0863,0.0199283,-6.98287e-06,1.10172e-09,-6.46381e-14,-20442.1,-69.7533], Tmin=(1380,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(C)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 284,
    label = "C#CCVCCJ",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,D} {3,S} {6,S}
2  C u0 p0 c0 {1,D} {4,S} {7,S}
3  C u1 p0 c0 {1,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {5,T}
5  C u0 p0 c0 {4,T} {10,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.616144,0.0506467,-4.48562e-05,2.02459e-08,-3.64542e-12,47153.2,27.1623], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[14.1231,0.0114233,-3.95851e-06,6.20129e-10,-3.62098e-14,42515.8,-50.2943], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""GLAR""",
    longDesc = 
u"""
GLAR
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#CC=C[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 285,
    label = "C5H93-1,2OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {10,S}
2  C u0 p0 c0 {4,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {1,S} {2,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.42923,0.0733174,-5.05533e-05,1.77709e-08,-2.53421e-12,-18921.3,22.6144], Tmin=(298,'K'), Tmax=(1413,'K')),
            NASAPolynomial(coeffs=[24.6778,0.0247897,-8.53168e-06,1.33018e-09,-7.7405e-14,-26400.9,-91.8794], Tmin=(1413,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]C(COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 286,
    label = "C5H10OOH3-2O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
6  O u0 p2 c0 {1,S} {8,S}
7  O u0 p2 c0 {2,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.94959,0.0794697,-6.123e-05,2.45361e-08,-3.99244e-12,-27626.3,21.0066], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[25.1554,0.0242477,-8.22109e-06,1.26886e-09,-7.33158e-14,-34952.2,-96.9776], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(OO)C(C)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 287,
    label = "C5H91-2,5OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {7,S} {15,S} {16,S}
5  C u1 p0 c0 {2,S} {17,S} {18,S}
6  O u0 p2 c0 {2,S} {9,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.56086,0.0780409,-5.74863e-05,2.18761e-08,-3.39924e-12,-17485.1,26.6589], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[25.0741,0.024438,-8.40705e-06,1.31038e-09,-7.62382e-14,-25187.1,-93.8668], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(CCCOO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 288,
    label = "NEOC5KETOX",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,D} {16,S}
6  O u1 p2 c0 {4,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.76959,0.0484037,-2.22216e-05,2.98251e-09,2.98065e-13,-25299.4,12.7919], Tmin=(298,'K'), Tmax=(1677,'K')),
            NASAPolynomial(coeffs=[16.2834,0.0247668,-8.89982e-06,1.43325e-09,-8.54233e-14,-29910,-56.0198], Tmin=(1677,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(C=O)C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 289,
    label = "IC3H5CO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,D} {4,S}
3  C u0 p0 c0 {2,D} {8,S} {9,S}
4  C u1 p0 c0 {2,S} {10,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 O u0 p2 c0 {4,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.87307,0.0395189,-3.11404e-05,1.28844e-08,-2.18165e-12,2852.71,16.8774], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[12.9634,0.0117955,-4.04361e-06,6.28772e-10,-3.6521e-14,-826.519,-42.0563], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 290,
    label = "C5H10OOH1-2O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {4,S} {8,S}
7  O u0 p2 c0 {1,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.92124,0.0767602,-5.10112e-05,1.59643e-08,-1.8298e-12,-25548.7,22.7306], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[27.1106,0.0234082,-8.1172e-06,1.27231e-09,-7.43222e-14,-34140.6,-108.093], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC(COO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 291,
    label = "IC4H7O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  O u0 p2 c0 {1,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 O u1 p2 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.43532,0.0489027,-3.63601e-05,1.41906e-08,-2.24558e-12,3092.85,24.132], Tmin=(298,'K'), Tmax=(1404,'K')),
            NASAPolynomial(coeffs=[14.5792,0.0162136,-5.26957e-06,7.90454e-10,-4.47755e-14,-1208.48,-45.6459], Tmin=(1404,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""L 2/00""",
    longDesc = 
u"""
L 2/00
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 292,
    label = "C2H5COC2H5",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {2,S} {16,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.04336,0.0424812,-1.12088e-05,-5.05052e-09,2.35367e-12,-34061.6,14.8305], Tmin=(298,'K'), Tmax=(1421,'K')),
            NASAPolynomial(coeffs=[16.5156,0.0218966,-7.26184e-06,1.10538e-09,-6.3291e-14,-39672.2,-61.1244], Tmin=(1421,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/6/15 THERM""",
    longDesc = 
u"""
10/6/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(=O)CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 293,
    label = "IC4H7-I1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u1 p0 c0 {3,D} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.912465,0.0388654,-2.57576e-05,9.0776e-09,-1.33947e-12,25563.6,20.8635], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[11.1159,0.0155127,-5.23769e-06,8.05998e-10,-4.64703e-14,21948.8,-34.144], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""5/13/15""",
    longDesc = 
u"""
5/13/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=C(C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 294,
    label = "C3H6OOH1-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u1 p0 c0 {1,S} {2,S} {11,S}
4  O u0 p2 c0 {1,S} {5,S}
5  O u0 p2 c0 {4,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.83631,0.038823,-2.47944e-05,7.85645e-09,-9.58634e-13,-1260.03,17.255], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[13.8089,0.0143846,-4.74441e-06,7.19308e-10,-4.10654e-14,-5143.53,-42.0211], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 295,
    label = "C4H8OOH1-4O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
5  O u0 p2 c0 {3,S} {7,S}
6  O u0 p2 c0 {4,S} {16,S}
7  O u0 p2 c0 {5,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 O u1 p2 c0 {6,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.91974,0.0634948,-4.29699e-05,1.42283e-08,-1.84506e-12,-20486.7,22.6279], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[22.6393,0.0198017,-6.9235e-06,1.091e-09,-6.39618e-14,-27544.2,-84.0748], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OCCCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 296,
    label = "C5H11O2-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {18,S}
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
17 H u0 p0 c0 {5,S}
18 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.95945,0.0673013,-4.74262e-05,1.76453e-08,-2.71599e-12,-15510.2,20.1536], Tmin=(298,'K'), Tmax=(1401,'K')),
            NASAPolynomial(coeffs=[20.1379,0.0243049,-8.24745e-06,1.27345e-09,-7.35954e-14,-21777.8,-77.2984], Tmin=(1401,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(CC)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 297,
    label = "H15DE25DM-AO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
4  C u0 p0 c0 {6,S} {13,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {3,S} {7,D}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  C u0 p0 c0 {6,D} {21,S} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 O u1 p2 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.03692,0.075064,-4.1778e-05,1.07916e-08,-9.97734e-13,6116.5,20.4746], Tmin=(298,'K'), Tmax=(1378,'K')),
            NASAPolynomial(coeffs=[25.0986,0.0315086,-1.10077e-05,1.73332e-09,-1.01557e-13,-2494.65,-101.098], Tmin=(1378,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)CCC(=C)C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 298,
    label = "C5H91-3OOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {14,S}
5  C u0 p0 c0 {4,D} {15,S} {16,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.8959,0.0818484,-6.75161e-05,2.83332e-08,-4.76827e-12,-14714.1,38.5399], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[22.466,0.0200408,-6.8455e-06,1.06266e-09,-6.16742e-14,-22614.9,-90.447], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(CC)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 299,
    label = "C5H10OOH1-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {2,S} {3,S} {17,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.65531,0.0613777,-3.72527e-05,1.10683e-08,-1.26576e-12,-6875.69,21.8611], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[20.0796,0.0243589,-8.36203e-06,1.3014e-09,-7.5633e-14,-13316.6,-73.0605], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[CH]COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 300,
    label = "C2H3COCH3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,D}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {10,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  O u0 p2 c0 {2,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.655911,0.0419405,-3.0027e-05,1.16578e-08,-1.92101e-12,-17221.7,23.8412], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[11.999,0.0152616,-5.26019e-06,8.20774e-10,-4.77835e-14,-21217.1,-37.1384], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 301,
    label = "C4H71-2,3OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {2,S} {14,S} {15,S}
5  O u0 p2 c0 {1,S} {8,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {16,S}
8  O u0 p2 c0 {5,S} {17,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.30791,0.0666513,-5.21745e-05,2.10434e-08,-3.42495e-12,-16753.7,19.5807], Tmin=(298,'K'), Tmax=(1406,'K')),
            NASAPolynomial(coeffs=[22.2829,0.0192254,-6.52861e-06,1.00882e-09,-5.83408e-14,-22972,-81.1232], Tmin=(1406,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(OO)C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 302,
    label = "TQJC4H8OH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {15,S}
6  O u0 p2 c0 {2,S} {16,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u1 p2 c0 {5,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.64342,0.0849132,-8.17211e-05,3.9098e-08,-7.27093e-12,-34237.6,28.4394], Tmin=(298,'K'), Tmax=(1415,'K')),
            NASAPolynomial(coeffs=[22.9682,0.0165163,-5.50247e-06,8.39335e-10,-4.81031e-14,-41005.1,-93.4898], Tmin=(1415,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(CO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 303,
    label = "C4H5OH-13",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,D} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {6,S}
3  C u0 p0 c0 {1,D} {5,S} {8,S}
4  C u0 p0 c0 {2,D} {9,S} {10,S}
5  O u0 p2 c0 {3,S} {11,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.60023,0.0612386,-6.18934e-05,3.14927e-08,-6.20397e-12,-9774.35,30.8328], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[14.0975,0.0129827,-4.36357e-06,6.6926e-10,-3.8492e-14,-14109.9,-49.4341], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/24/15""",
    longDesc = 
u"""
9/24/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 304,
    label = "C5H92-1OOH",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {5,D} {15,S}
5  C u0 p0 c0 {2,S} {4,D} {16,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.57788,0.0691529,-5.00924e-05,1.881e-08,-2.91765e-12,-14404.3,29.7057], Tmin=(298,'K'), Tmax=(1382,'K')),
            NASAPolynomial(coeffs=[21.0624,0.0215772,-7.44877e-06,1.16424e-09,-6.7882e-14,-21581.7,-80.4706], Tmin=(1382,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC=CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 305,
    label = "C5H10OOH2-5O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {8,S}
7  O u0 p2 c0 {4,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.92124,0.0767602,-5.10112e-05,1.59643e-08,-1.8298e-12,-25548.7,22.7306], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[27.1106,0.0234082,-8.1172e-06,1.27231e-09,-7.43222e-14,-34140.6,-108.093], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CCCO[O])OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 306,
    label = "C4H71-2,4OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {14,S} {15,S}
5  O u0 p2 c0 {1,S} {8,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {16,S}
8  O u0 p2 c0 {5,S} {17,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.99387,0.0652914,-4.87959e-05,1.88108e-08,-2.95257e-12,-14602.1,22.7764], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[21.863,0.0199359,-6.85104e-06,1.06712e-09,-6.20562e-14,-21004.3,-78.0707], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(CCOO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 307,
    label = "IC3H5COCH3",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u0 p0 c0 {2,S} {3,S} {12,D}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 O u0 p2 c0 {4,D}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.54292,0.0505432,-3.42434e-05,1.25252e-08,-1.95302e-12,-21270.3,19.9856], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[14.9936,0.0197788,-6.78313e-06,1.05481e-09,-6.12615e-14,-26087.8,-52.6202], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)C(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 308,
    label = "AC5H11O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
6  O u0 p2 c0 {3,S} {18,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.354375,0.0721759,-5.35556e-05,2.12668e-08,-3.51578e-12,-13832,29.9445], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[19.8701,0.0247772,-8.46414e-06,1.31284e-09,-7.61129e-14,-20483.2,-74.3157], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 309,
    label = "C3H5OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,D} {7,S}
3  C u0 p0 c0 {2,D} {8,S} {9,S}
4  O u0 p2 c0 {1,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.15012,0.0128538,4.28438e-05,-6.67819e-08,2.80408e-11,-16641.4,13.5066], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.72477,0.0163943,-5.90853e-06,9.53262e-10,-5.70318e-14,-19049.7,-19.7199], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T06/10""",
    longDesc = 
u"""
T06/10.
C=CCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 310,
    label = "CCYCCOOC-T1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {2,S} {5,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.29768,0.0665082,-4.67054e-05,1.5103e-08,-1.74322e-12,3979.36,51.6412], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[16.827,0.0164472,-5.63767e-06,8.78002e-10,-5.1096e-14,-3657.11,-67.4097], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C]1COOC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 311,
    label = "C5H9O1-4OOH-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {2,S} {4,S}
7  O u0 p2 c0 {1,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-9.18085,0.108786,-9.86756e-05,4.49208e-08,-8.03088e-12,-39801.3,72.6761], Tmin=(298,'K'), Tmax=(1417,'K')),
            NASAPolynomial(coeffs=[21.726,0.0234146,-7.86188e-06,1.20572e-09,-6.93654e-14,-49050.9,-88.3591], Tmin=(1417,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1OCCC1OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 312,
    label = "C5H5OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {1,S} {5,D} {9,S}
4  C u0 p0 c0 {2,D} {5,S} {10,S}
5  C u0 p0 c0 {3,D} {4,S} {11,S}
6  O u0 p2 c0 {1,S} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.26822,0.0662447,-5.68494e-05,2.46859e-08,-4.26821e-12,-5755.81,44.7963], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[15.3433,0.0150754,-5.13554e-06,7.95808e-10,-4.61312e-14,-11964.5,-58.5204], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""5/ 2/91 THE.M""",
    longDesc = 
u"""
5/ 2/91 THE.M
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OC1C=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 313,
    label = "QC4H7OHP",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {2,S} {15,S}
7  O u0 p2 c0 {5,S} {16,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.27864,0.0894493,-8.78565e-05,4.22111e-08,-7.83451e-12,-25897.5,35.3964], Tmin=(298,'K'), Tmax=(1416,'K')),
            NASAPolynomial(coeffs=[24.3481,0.0150316,-5.01788e-06,7.66774e-10,-4.40093e-14,-33192.2,-96.8211], Tmin=(1416,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)(CO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 314,
    label = "IC4H7OOCH3",
    molecule = 
"""
1  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
3  C u0 p0 c0 {7,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {4,D} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {3,S} {6,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.26784,0.0682442,-5.30808e-05,2.27496e-08,-4.11475e-12,-11676.8,24.4901], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[19.5897,0.0231057,-8.02911e-06,1.2597e-09,-7.36169e-14,-18006.9,-73.4192], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)COOC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 315,
    label = "DC5H10OOH-C",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
5  C u1 p0 c0 {1,S} {4,S} {17,S}
6  O u0 p2 c0 {4,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.46802,0.0681102,-4.85321e-05,1.84083e-08,-2.9172e-12,-7389.26,25.7095], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[20.0395,0.0242833,-8.31028e-06,1.29054e-09,-7.48846e-14,-13846,-73.9591], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)[CH]COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 316,
    label = "C5H92-3,4OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {2,S} {4,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {2,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.76788,0.0747668,-5.4178e-05,2.01841e-08,-3.04868e-12,-21082.3,19.2562], Tmin=(298,'K'), Tmax=(1418,'K')),
            NASAPolynomial(coeffs=[25.107,0.0240686,-8.20467e-06,1.27103e-09,-7.36339e-14,-28359.1,-94.9638], Tmin=(1418,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C(OO)C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 317,
    label = "C5H10OH14",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {2,S} {4,S} {16,S}
6  O u0 p2 c0 {3,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.541168,0.0566862,-2.94678e-05,5.84905e-09,-3.31866e-14,-15094.4,31.8827], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[16.8143,0.0244876,-8.29373e-06,1.27921e-09,-7.38805e-14,-21256.8,-57.4573], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 7/15 THERM""",
    longDesc = 
u"""
10/ 7/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 318,
    label = "C5H9O2-3OOH-4",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {9,S}
4  C u0 p0 c0 {3,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {2,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.57026,0.105196,-0.000105064,5.15547e-08,-9.70822e-12,-32044.2,49.0446], Tmin=(298,'K'), Tmax=(1430,'K')),
            NASAPolynomial(coeffs=[23.6692,0.019256,-6.07857e-06,8.92372e-10,-4.97449e-14,-39629.6,-95.1544], Tmin=(1430,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(OO)C1OC1C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 319,
    label = "C4H7O1-3OOH-4",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
5  O u0 p2 c0 {1,S} {3,S}
6  O u0 p2 c0 {4,S} {7,S}
7  O u0 p2 c0 {6,S} {15,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.4095,0.0805103,-6.64646e-05,2.73712e-08,-4.44857e-12,-25300.5,55.6327], Tmin=(298,'K'), Tmax=(1418,'K')),
            NASAPolynomial(coeffs=[18.725,0.0188463,-6.40166e-06,9.89755e-10,-5.72728e-14,-32982.6,-71.8225], Tmin=(1418,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOCC1CCO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 320,
    label = "C6H5O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,D}
2  C u0 p0 c0 {1,S} {3,D} {9,S}
3  C u0 p0 c0 {2,D} {4,S} {10,S}
4  C u1 p0 c0 {3,S} {5,S} {8,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {1,S} {5,D} {12,S}
7  O u0 p2 c0 {1,D}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.466204,0.0413444,1.32413e-05,-5.72873e-08,2.89764e-11,4778.58,27.699], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[13.7222,0.0174689,-6.35505e-06,1.03492e-09,-6.23411e-14,287.275,-48.8182], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T05/02""",
    longDesc = 
u"""
T05/02.
O=C1C=C[CH]C=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 321,
    label = "IC3H7CO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {1,S} {12,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u0 p2 c0 {4,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.503453,0.0441608,-2.82139e-05,8.93549e-09,-1.11327e-12,-9077.55,26.1991], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[13.3306,0.0161874,-5.56711e-06,8.67576e-10,-5.04697e-14,-13730.7,-43.3959], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/22/96 THERM""",
    longDesc = 
u"""
2/22/96 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 322,
    label = "C4H7O1-4OOH-2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
5  O u0 p2 c0 {3,S} {4,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {15,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.17427,0.0789226,-6.62654e-05,2.77699e-08,-4.54378e-12,-35378.1,53.5508], Tmin=(298,'K'), Tmax=(1470,'K')),
            NASAPolynomial(coeffs=[17.4906,0.0187795,-6.05492e-06,9.03186e-10,-5.09572e-14,-42271.7,-65.186], Tmin=(1470,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOC1CCOC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 323,
    label = "AC5H10OOH-B",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
2  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {1,S} {2,S} {4,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.65134,0.0570948,-3.30635e-05,9.5243e-09,-1.08604e-12,-9155.25,17.3958], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[19.0362,0.0251438,-8.60836e-06,1.33716e-09,-7.76006e-14,-14972.1,-66.7702], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[C](C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 324,
    label = "C5H92-4,5OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
4  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {2,S} {4,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.06621,0.0720274,-4.95956e-05,1.70909e-08,-2.31584e-12,-19252.3,19.9317], Tmin=(298,'K'), Tmax=(1413,'K')),
            NASAPolynomial(coeffs=[24.7515,0.024146,-8.18239e-06,1.26251e-09,-7.29352e-14,-26378.3,-91.1444], Tmin=(1413,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CC(COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 325,
    label = "C5H10OOH1-4O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {4,S} {8,S}
7  O u0 p2 c0 {1,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.92124,0.0767602,-5.10112e-05,1.59643e-08,-1.8298e-12,-25548.7,22.7306], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[27.1106,0.0234082,-8.1172e-06,1.27231e-09,-7.43222e-14,-34140.6,-108.093], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CCCOO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 326,
    label = "IC3H5O2HCHO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,D} {13,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u0 p2 c0 {4,D}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.05985,0.0582332,-4.37672e-05,1.6325e-08,-2.43462e-12,-16349.6,21.3688], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[20.6289,0.0148626,-5.25305e-06,8.33773e-10,-4.91277e-14,-22758.9,-78.2963], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/2/95 THERM""",
    longDesc = 
u"""
8/2/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)(C=O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 327,
    label = "C5H10O1-2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
6  O u0 p2 c0 {1,S} {4,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.22385,0.0795281,-6.27918e-05,2.51187e-08,-3.99106e-12,-17615.2,51.5328], Tmin=(298,'K'), Tmax=(1424,'K')),
            NASAPolynomial(coeffs=[17.6453,0.0217647,-7.28831e-06,1.11572e-09,-6.41056e-14,-24962.8,-69.4781], Tmin=(1424,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 328,
    label = "C5H11O2H-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {19,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.06841,0.0736472,-5.31234e-05,2.02381e-08,-3.18672e-12,-32174.7,23.8112], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[21.1348,0.0256494,-8.7265e-06,1.34983e-09,-7.81087e-14,-39041.6,-83.5754], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(CC)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 329,
    label = "CCY(CCO)COH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {2,S}
6  O u0 p2 c0 {3,S} {14,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.10049,0.0953372,-9.78702e-05,4.90006e-08,-9.41686e-12,-39898.7,56.0925], Tmin=(298,'K'), Tmax=(1412,'K')),
            NASAPolynomial(coeffs=[19.1885,0.0156256,-5.2257e-06,7.99171e-10,-4.58832e-14,-47112,-78.4579], Tmin=(1412,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(CO)CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 330,
    label = "C5H10OOH2-5",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {3,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.46239,0.0700774,-5.18935e-05,2.03024e-08,-3.27212e-12,-7458.96,25.6256], Tmin=(298,'K'), Tmax=(1401,'K')),
            NASAPolynomial(coeffs=[20.5595,0.0235924,-8.01737e-06,1.23919e-09,-7.16674e-14,-13905.9,-76.2752], Tmin=(1401,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 331,
    label = "CH3CH2CHCOCH3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {15,S}
5  C u0 p0 c0 {3,S} {4,D} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 O u1 p2 c0 {5,S}
15 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.6218,0.0378869,-9.52099e-06,-4.70141e-09,2.11457e-12,-14664.3,16.3247], Tmin=(298,'K'), Tmax=(1428,'K')),
            NASAPolynomial(coeffs=[15.304,0.0203924,-6.74389e-06,1.02444e-09,-5.85693e-14,-19576.5,-49.6835], Tmin=(1428,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/6/15""",
    longDesc = 
u"""
10/6/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC=C(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 332,
    label = "NC5KET14",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,D} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.00527,0.0657859,-4.36349e-05,1.39803e-08,-1.68002e-12,-42409.5,19.0078], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[22.5902,0.0217812,-7.38064e-06,1.13911e-09,-6.58296e-14,-49267.9,-86.5971], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CCC=O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 333,
    label = "IC4H8O2H-T",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {4,S} {12,S} {13,S} {14,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.84375,0.0436801,-2.076e-05,2.51709e-09,5.41307e-13,-6507.66,13.4245], Tmin=(298,'K'), Tmax=(1413,'K')),
            NASAPolynomial(coeffs=[16.9754,0.0185198,-6.09075e-06,9.21674e-10,-5.25503e-14,-11481.3,-58.8259], Tmin=(1413,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C](C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 334,
    label = "C5H10OOH3-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {1,S} {4,S} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.03845,0.0635523,-4.29455e-05,1.50666e-08,-2.16139e-12,-9066.91,17.4223], Tmin=(298,'K'), Tmax=(1410,'K')),
            NASAPolynomial(coeffs=[20.36,0.0236444,-8.00997e-06,1.23544e-09,-7.13457e-14,-15126.4,-75.7807], Tmin=(1410,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C(CC)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 335,
    label = "SC4H7OH-I",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {5,S} {12,S}
5  O u0 p2 c0 {4,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.70103,0.041795,-2.67861e-05,9.38191e-09,-1.41171e-12,-27356.1,13.5316], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[13.0299,0.0183782,-6.1853e-06,9.49578e-10,-5.46526e-14,-31072.3,-42.2892], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""L 2/00""",
    longDesc = 
u"""
L 2/00
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 336,
    label = "C4H8OOH1-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {3,S} {14,S}
5  O u0 p2 c0 {2,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.46293,0.0470131,-2.42145e-05,4.65408e-09,1.62199e-14,-3957.88,22.457], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[16.1248,0.0202421,-6.88631e-06,1.06535e-09,-6.166e-14,-9168.02,-52.6575], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 337,
    label = "C4H8OOH1-4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u1 p0 c0 {2,S} {13,S} {14,S}
5  O u0 p2 c0 {3,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.35984,0.0552812,-3.75641e-05,1.32624e-08,-1.93623e-12,-2344.56,26.3198], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[16.9218,0.0199306,-6.85728e-06,1.06879e-09,-6.21775e-14,-7879.21,-57.663], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 338,
    label = "C4H8OOH1-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {1,S} {2,S} {14,S}
5  O u0 p2 c0 {2,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.91878,0.0486144,-2.81342e-05,7.64456e-09,-7.32889e-13,-3944.65,18.9326], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[16.781,0.0199678,-6.85257e-06,1.06629e-09,-6.19616e-14,-9147.63,-56.8668], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 339,
    label = "C3H6O1-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  O u0 p2 c0 {2,S} {3,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.15284,-0.0186402,0.000129981,-1.5863e-07,6.20669e-11,-11324.4,4.73561], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.80717,0.0188825,-6.79082e-06,1.09714e-09,-6.57155e-14,-13654.8,-13.5382], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""A11/04""",
    longDesc = 
u"""
A11/04.
C1COC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 340,
    label = "C5H9O13-5",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u1 p0 c0 {3,S} {14,S} {15,S}
6  O u0 p2 c0 {1,S} {4,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-6.25394,0.0792185,-6.70994e-05,2.87608e-08,-4.83752e-12,6546.69,58.0477], Tmin=(298,'K'), Tmax=(1458,'K')),
            NASAPolynomial(coeffs=[15.599,0.0200371,-6.41155e-06,9.50808e-10,-5.34042e-14,3.37602,-56.0383], Tmin=(1458,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC1CCO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 341,
    label = "C4H8O1-4",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
5  O u0 p2 c0 {3,S} {4,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.78118,0.0698405,-5.45316e-05,2.1303e-08,-3.24667e-12,-22899.6,61.7621], Tmin=(298,'K'), Tmax=(1484,'K')),
            NASAPolynomial(coeffs=[12.2763,0.0189106,-6.09114e-06,9.08066e-10,-5.1215e-14,-29226.1,-44.1236], Tmin=(1484,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C1CCOC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 342,
    label = "C3KET21",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {11,D}
4  O u0 p2 c0 {1,S} {5,S}
5  O u0 p2 c0 {4,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 O u0 p2 c0 {3,D}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.874353,0.0612501,-5.51475e-05,2.48491e-08,-4.42613e-12,-35806.1,35.9306], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[17.5768,0.0120312,-4.11634e-06,6.40149e-10,-3.72128e-14,-41550.2,-60.9097], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/14/13 THERM""",
    longDesc = 
u"""
2/14/13 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 343,
    label = "C3H3O",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {7,S}
4 O u1 p2 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.875024,0.0351184,-3.89901e-05,2.40256e-08,-6.10884e-12,32042.8,20.4717], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.19356,0.0195625,-1.22336e-05,3.90615e-09,-5.08539e-13,31493.2,5.03216], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = u"""2/17/14 CZHOU""",
    longDesc = 
u"""
2/17/14 CZHOU
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C#CC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 344,
    label = "CC3H4",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {3,D} {6,S}
3 C u0 p0 c0 {1,S} {2,D} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.024621,0.0231972,-1.84744e-06,-1.59276e-08,8.68462e-12,32334.1,22.7298], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.69999,0.0103574,-3.45512e-06,5.06529e-10,-2.66823e-14,30199.1,-13.3788], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""T12/81""",
    longDesc = 
u"""
T12/81
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C1=CC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 345,
    label = "C2H2OH",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u1 p0 c0 {1,D} {5,S}
3 O u0 p2 c0 {1,S} {6,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.641643,0.0261904,-2.30385e-05,1.02805e-08,-1.81971e-12,14827.7,20.6751], Tmin=(298,'K'), Tmax=(1401,'K')),
            NASAPolynomial(coeffs=[8.20268,0.00592989,-1.99194e-06,3.05794e-10,-1.76115e-14,12488.1,-18.967], Tmin=(1401,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 346,
    label = "C2H5COOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {10,D}
4  O u0 p2 c0 {3,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 O u0 p2 c0 {3,D}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.51515,-1.02654e-05,7.54897e-05,-9.81745e-08,3.8883e-11,-56361.8,6.11744], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.61037,0.0187895,-6.70711e-06,1.07428e-09,-6.3887e-14,-58480.8,-16.311], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T11/07""",
    longDesc = 
u"""
T11/07.
CCC(=O)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 347,
    label = "CH3COCH2O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {10,D}
4  O u0 p2 c0 {1,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 O u0 p2 c0 {3,D}
11 O u1 p2 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.19378,0.0498027,-4.18e-05,1.74528e-08,-2.88199e-12,-19324.4,26.7877], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[16.5756,0.0106465,-3.61369e-06,5.59054e-10,-3.23832e-14,-24254.1,-54.5305], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/14/13 THERM""",
    longDesc = 
u"""
2/14/13 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 348,
    label = "C2CYCOOC-I1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {1,S} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {2,S} {5,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.04718,0.0726608,-6.8296e-05,3.27506e-08,-6.23802e-12,-371.476,19.2056], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[18.9745,0.0150114,-5.30728e-06,8.42455e-10,-4.96393e-14,-6931.11,-90.8581], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/14""",
    longDesc = 
u"""
7/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C1(C)COO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 349,
    label = "C5H92-3,5OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {7,S} {13,S} {14,S}
4  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {1,S} {4,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.42923,0.0733174,-5.05533e-05,1.77709e-08,-2.53421e-12,-18921.3,22.6144], Tmin=(298,'K'), Tmax=(1413,'K')),
            NASAPolynomial(coeffs=[24.6778,0.0247897,-8.53168e-06,1.33018e-09,-7.7405e-14,-26400.9,-91.8794], Tmin=(1413,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C(CCOO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 350,
    label = "C4H8OOH2-3O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
3  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {2,S} {16,S}
7  O u0 p2 c0 {5,S} {17,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u1 p2 c0 {6,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.34684,0.0667468,-5.25089e-05,2.14288e-08,-3.53382e-12,-24734.1,17.3189], Tmin=(298,'K'), Tmax=(1408,'K')),
            NASAPolynomial(coeffs=[21.9463,0.0197308,-6.65765e-06,1.02425e-09,-5.90486e-14,-30777.2,-81.2054], Tmin=(1408,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(O[O])C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 351,
    label = "C6H5OOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,D} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {9,S}
3  C u0 p0 c0 {1,D} {6,S} {13,S}
4  C u0 p0 c0 {2,D} {5,S} {10,S}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {3,S} {5,D} {12,S}
7  O u0 p2 c0 {1,S} {8,S}
8  O u0 p2 c0 {7,S} {14,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.03106,0.0796102,-7.21655e-05,3.27611e-08,-5.85584e-12,-3109.73,45.4325], Tmin=(298,'K'), Tmax=(1404,'K')),
            NASAPolynomial(coeffs=[19.2317,0.0163155,-5.53449e-06,8.5506e-10,-4.94584e-14,-10197.1,-76.1674], Tmin=(1404,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""3/26/ 9 THERM""",
    longDesc = 
u"""
3/26/ 9 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOC1C=CC=CC=1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 352,
    label = "C6H5OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,D} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {1,D} {6,S} {12,S}
4  C u0 p0 c0 {2,D} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {3,S} {5,D} {11,S}
7  O u0 p2 c0 {1,S} {13,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.69565,0.0522713,-7.20241e-06,-3.58596e-08,2.04491e-11,-13284.1,32.5422], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[14.9121,0.0183781,-6.19831e-06,9.19832e-10,-4.92096e-14,-18375.2,-55.9241], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""L 4/84""",
    longDesc = 
u"""
L 4/84
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OC1C=CC=CC=1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 353,
    label = "C4H71-3,4OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u1 p0 c0 {2,S} {14,S} {15,S}
5  O u0 p2 c0 {1,S} {8,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {16,S}
8  O u0 p2 c0 {5,S} {17,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5689,0.0671962,-5.17e-05,2.0577e-08,-3.32914e-12,-14512.6,24.6333], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[21.8395,0.0198803,-6.81505e-06,1.05975e-09,-6.1556e-14,-20946.9,-77.9985], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 354,
    label = "C6H4OH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,D} {7,S}
2  C u0 p0 c0 {1,S} {3,D} {8,S}
3  C u0 p0 c0 {2,D} {5,S} {9,S}
4  C u0 p0 c0 {1,D} {6,S} {11,S}
5  C u0 p0 c0 {3,S} {6,D} {10,S}
6  C u1 p0 c0 {4,S} {5,D}
7  O u0 p2 c0 {1,S} {12,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.99875,0.0859063,-9.12526e-05,4.72276e-08,-9.35577e-12,17862.2,49.9931], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[17.3188,0.0136367,-4.68316e-06,7.29071e-10,-4.23805e-14,11499,-68.9987], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""4/ 9/ 9 THERM""",
    longDesc = 
u"""
4/ 9/ 9 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OC1C=CC=[C]C=1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 355,
    label = "H15DE25DM-SO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {6,S} {9,S} {12,S}
3  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {6,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {3,S} {7,D}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  C u0 p0 c0 {6,D} {21,S} {22,S}
9  O u1 p2 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.930464,0.087167,-5.99583e-05,2.15026e-08,-3.21688e-12,4193.1,29.0925], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[26.0152,0.0304296,-1.05677e-05,1.65747e-09,-9.6847e-14,-4795.69,-106.442], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)CC([O])C(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 356,
    label = "H15DE25DM-S",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
2  C u0 p0 c0 {4,S} {14,S} {15,S} {16,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {2,S} {7,D}
5  C u0 p0 c0 {3,S} {6,D} {8,S}
6  C u0 p0 c0 {1,S} {5,D} {17,S}
7  C u0 p0 c0 {4,D} {20,S} {21,S}
8  C u1 p0 c0 {5,S} {18,S} {19,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.04236,0.0866262,-5.99111e-05,2.15552e-08,-3.18977e-12,15022.4,39.6918], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[22.1423,0.0306966,-1.04617e-05,1.62038e-09,-9.38579e-14,6560.46,-90.406], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)=CCC(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 357,
    label = "TC3H6O2CHO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,D} {13,S}
5  O u0 p2 c0 {1,S} {14,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u0 p2 c0 {4,D}
13 H u0 p0 c0 {4,S}
14 O u1 p2 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.17883,0.0541596,-3.83436e-05,1.38308e-08,-2.0419e-12,-22739.4,20.0751], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[18.5534,0.0168774,-5.90753e-06,9.31518e-10,-5.46345e-14,-28544.7,-68.2487], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/ 2/95 THERM""",
    longDesc = 
u"""
8/ 2/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(C=O)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 358,
    label = "C5H6-L",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,D} {9,S}
3  C u0 p0 c0 {2,D} {4,S} {10,S}
4  C u0 p0 c0 {3,S} {5,T}
5  C u0 p0 c0 {4,T} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.58448,0.032446,-1.70151e-05,4.22716e-09,-4.18453e-13,27651.5,9.60644], Tmin=(298,'K'), Tmax=(1372,'K')),
            NASAPolynomial(coeffs=[12.9601,0.0148954,-5.23623e-06,8.27916e-10,-4.86465e-14,23818.1,-42.5312], Tmin=(1372,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/ 5/ 9 THERM""",
    longDesc = 
u"""
2/ 5/ 9 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#CC=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 359,
    label = "NEOC5KEJOL",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {1,S} {15,D}
6  O u0 p2 c0 {2,S} {16,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.30299,0.0584133,-3.81961e-05,1.24622e-08,-1.59437e-12,-32897.7,25.9155], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[17.8858,0.0212423,-7.21948e-06,1.11609e-09,-6.45626e-14,-38768.1,-63.6234], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)([C]=O)CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 360,
    label = "IC4H8OH-IT",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  O u0 p2 c0 {1,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.05276,0.0393926,-1.90686e-05,3.86408e-09,-1.48005e-13,-14226.4,16.0841], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[12.9137,0.0206583,-6.98446e-06,1.07563e-09,-6.20444e-14,-18139.5,-38.4972], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C](C)CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 361,
    label = "IC4H7OOH",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.131852,0.0619561,-4.99344e-05,2.09628e-08,-3.59718e-12,-12140,29.3906], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[18.2897,0.0167816,-5.80668e-06,9.08949e-10,-5.30513e-14,-18204.7,-67.2111], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""4/15/15""",
    longDesc = 
u"""
4/15/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 362,
    label = "C5H10OOH1-3O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {4,S} {8,S}
7  O u0 p2 c0 {1,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.92124,0.0767602,-5.10112e-05,1.59643e-08,-1.8298e-12,-25548.7,22.7306], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[27.1106,0.0234082,-8.1172e-06,1.27231e-09,-7.43222e-14,-34140.6,-108.093], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(CCOO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 363,
    label = "C4H8OOH1-2O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {3,S} {7,S}
6  O u0 p2 c0 {1,S} {16,S}
7  O u0 p2 c0 {5,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u1 p2 c0 {6,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.02241,0.0653863,-4.89646e-05,1.90438e-08,-3.02309e-12,-22580.6,20.5729], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[21.5735,0.0204529,-6.99498e-06,1.08598e-09,-6.30071e-14,-28842.8,-78.4561], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(COO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 364,
    label = "C5H10OOH3-1O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {8,S}
7  O u0 p2 c0 {4,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.92124,0.0767602,-5.10112e-05,1.59643e-08,-1.8298e-12,-25548.7,22.7306], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[27.1106,0.0234082,-8.1172e-06,1.27231e-09,-7.43222e-14,-34140.6,-108.093], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(CCO[O])OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 365,
    label = "C4H7O1-2OOH-4",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
5  O u0 p2 c0 {1,S} {3,S}
6  O u0 p2 c0 {4,S} {7,S}
7  O u0 p2 c0 {6,S} {15,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.06057,0.0752777,-6.16382e-05,2.51521e-08,-4.05109e-12,-24931,44.8871], Tmin=(298,'K'), Tmax=(1417,'K')),
            NASAPolynomial(coeffs=[19.6187,0.0177382,-6.03564e-06,9.34268e-10,-5.41073e-14,-32191.8,-75.0296], Tmin=(1417,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOCCC1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 366,
    label = "C5H91-5OOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {5,D} {14,S}
5  C u0 p0 c0 {4,D} {15,S} {16,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.993374,0.0646025,-4.40192e-05,1.55375e-08,-2.26236e-12,-14715.7,27.8588], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[19.3103,0.0229913,-7.9127e-06,1.23359e-09,-7.17773e-14,-21224.7,-70.9826], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 367,
    label = "C5H92-1,4OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
5  C u1 p0 c0 {2,S} {4,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.42923,0.0733174,-5.05533e-05,1.77709e-08,-2.53421e-12,-18921.3,22.6144], Tmin=(298,'K'), Tmax=(1413,'K')),
            NASAPolynomial(coeffs=[24.6778,0.0247897,-8.53168e-06,1.33018e-09,-7.7405e-14,-26400.9,-91.8794], Tmin=(1413,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C[CH]COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 368,
    label = "IC4H6OH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,D} {4,S}
3  C u0 p0 c0 {2,D} {5,S} {9,S}
4  C u1 p0 c0 {2,S} {10,S} {11,S}
5  O u0 p2 c0 {3,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.46664,0.0603352,-5.43113e-05,2.493e-08,-4.52282e-12,-6950.12,32.0768], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[15.3491,0.0138857,-4.56428e-06,6.90419e-10,-3.9354e-14,-12016.5,-55.5976], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 369,
    label = "C8H131-5,3,TAO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {7,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {8,D} {18,S}
6  C u0 p0 c0 {2,S} {7,D} {19,S}
7  C u0 p0 c0 {4,S} {6,D} {20,S}
8  C u0 p0 c0 {5,D} {21,S} {22,S}
9  O u1 p2 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.46772,0.0767801,-5.6319e-05,2.2269e-08,-3.64814e-12,1174.76,14.562], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[21.6639,0.0277253,-9.69376e-06,1.67582e-09,-1.08059e-13,-5703.01,-93.3238], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)([O])CC=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 370,
    label = "C8H131-5,3-4,TAO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {7,D} {17,S}
6  C u0 p0 c0 {2,S} {8,D} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  C u0 p0 c0 {6,D} {21,S} {22,S}
9  O u1 p2 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.0068,0.0890191,-6.27964e-05,2.24162e-08,-3.20587e-12,4564.23,27.6547], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[27.0579,0.0284727,-9.66435e-06,1.49337e-09,-8.63788e-14,-4431.11,-112.228], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)C(C)([O])C=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 371,
    label = "TC4H9O2H",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {16,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.445574,0.0666154,-5.20932e-05,2.22302e-08,-4.0019e-12,-31226.1,23.7278], Tmin=(298,'K'), Tmax=(1382,'K')),
            NASAPolynomial(coeffs=[19.0927,0.0212698,-7.46253e-06,1.17841e-09,-6.91795e-14,-37727.8,-76.1321], Tmin=(1382,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 372,
    label = "SC3H5CO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,D} {8,S}
3  C u0 p0 c0 {2,D} {4,S} {9,S}
4  C u1 p0 c0 {3,S} {10,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 O u0 p2 c0 {4,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.776401,0.0426828,-3.37881e-05,1.39128e-08,-2.33332e-12,1299.03,19.8102], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[12.9926,0.0122141,-4.19305e-06,6.52698e-10,-3.79407e-14,-2747.82,-45.1092], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 373,
    label = "IC4H9O2H",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {2,S} {6,S}
6  O u0 p2 c0 {5,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.547969,0.0657429,-4.70123e-05,1.65991e-08,-2.27331e-12,-28218,31.2961], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[19.1652,0.0193679,-6.41984e-06,9.76948e-10,-5.593e-14,-34879.1,-74.2064], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 374,
    label = "NC3H7CO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {2,S} {12,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u0 p2 c0 {4,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.63538,0.0340369,-1.24119e-05,-1.17887e-09,1.16488e-12,-8659.2,16.8408], Tmin=(298,'K'), Tmax=(1496,'K')),
            NASAPolynomial(coeffs=[13.487,0.0158627,-5.41699e-06,8.40509e-10,-4.8757e-14,-13072.5,-43.8634], Tmin=(1496,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 375,
    label = "C5H4OH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,D} {6,S}
2  C u1 p0 c0 {1,S} {4,S} {7,S}
3  C u0 p0 c0 {1,D} {5,S} {10,S}
4  C u0 p0 c0 {2,S} {5,D} {8,S}
5  C u0 p0 c0 {3,S} {4,D} {9,S}
6  O u0 p2 c0 {1,S} {11,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.28398,0.0490299,-1.35844e-05,-2.92984e-08,1.90821e-11,6373.65,30.8074], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[13.3741,0.0151996,-5.45685e-06,8.80945e-10,-5.27493e-14,2203.58,-45.9569], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 8/99""",
    longDesc = 
u"""
T 8/99.
OC1[CH]C=CC=1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 376,
    label = "IC3H5COCH2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u0 p0 c0 {2,S} {5,S} {9,D}
4  C u0 p0 c0 {2,D} {12,S} {13,S}
5  C u1 p0 c0 {3,S} {10,S} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  O u0 p2 c0 {3,D}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.42098,0.0534337,-4.21904e-05,1.76938e-08,-3.06075e-12,691.083,19.9522], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[16.2815,0.0163493,-5.64821e-06,8.82756e-10,-5.14522e-14,-4286.51,-59.1076], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=O)C(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 377,
    label = "AC3H4COCH3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {3,S} {4,S} {5,D}
3  C u0 p0 c0 {1,S} {2,S} {9,D}
4  C u1 p0 c0 {2,S} {10,S} {11,S}
5  C u0 p0 c0 {2,D} {12,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  O u0 p2 c0 {3,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.948806,0.0517869,-3.74868e-05,1.40295e-08,-2.14699e-12,-3105.54,22.6479], Tmin=(298,'K'), Tmax=(1384,'K')),
            NASAPolynomial(coeffs=[16.0406,0.0162939,-5.49275e-06,8.45442e-10,-4.87866e-14,-8309.73,-58.3003], Tmin=(1384,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=C)C(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 378,
    label = "C3H6CHO-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u1 p0 c0 {1,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,D} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 O u0 p2 c0 {4,D}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.21483,0.0354113,-1.44083e-05,5.27606e-11,8.95003e-13,-2464.83,20.0076], Tmin=(298,'K'), Tmax=(1538,'K')),
            NASAPolynomial(coeffs=[13.3449,0.0159347,-5.43144e-06,8.41706e-10,-4.87847e-14,-6912.82,-41.9662], Tmin=(1538,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 379,
    label = "CH2CHOCH2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3 C u1 p0 c0 {1,S} {4,S} {9,S}
4 O u0 p2 c0 {2,S} {3,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.1535,0.0351254,-2.50072e-05,9.00716e-09,-1.32377e-12,10830.1,21.0607], Tmin=(298,'K'), Tmax=(1384,'K')),
            NASAPolynomial(coeffs=[12.0077,0.0105055,-3.69921e-06,5.8563e-10,-3.44432e-14,6973.12,-37.519], Tmin=(1384,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/ 8/15""",
    longDesc = 
u"""
8/ 8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]1CCO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 380,
    label = "C5H9O2-4OOH-1",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {6,S} {10,S}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {2,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.54952,0.0955182,-8.0162e-05,3.35588e-08,-5.52327e-12,-30361.1,56.4175], Tmin=(298,'K'), Tmax=(1429,'K')),
            NASAPolynomial(coeffs=[22.6535,0.0220609,-7.36975e-06,1.12674e-09,-6.46941e-14,-39161.6,-91.9654], Tmin=(1429,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1CC(COO)O1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 381,
    label = "C5H92-5OOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {5,D} {16,S}
5  C u0 p0 c0 {3,S} {4,D} {15,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.76945,0.0600153,-3.72031e-05,1.14551e-08,-1.39226e-12,-16124.4,23.8819], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[18.8657,0.0233222,-8.01684e-06,1.24883e-09,-7.26245e-14,-22413,-69.1317], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 382,
    label = "OC6H4OH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,D} {7,S}
2  C u1 p0 c0 {1,S} {6,S} {9,S}
3  C u0 p0 c0 {1,D} {4,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {8,D}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {2,S} {5,D} {11,S}
7  O u0 p2 c0 {1,S} {13,S}
8  O u0 p2 c0 {4,D}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.02206,0.109403,-0.000123489,6.56287e-08,-1.31528e-11,-15594.9,57.2175], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[22.2718,0.0121039,-4.1843e-06,6.54475e-10,-3.81747e-14,-23482.8,-96.1035], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""4/ 9/ 9 THERM""",
    longDesc = 
u"""
4/ 9/ 9 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=C1C=C[CH]C(O)=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 383,
    label = "C4H6O23",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {5,S} {11,S}
5  O u0 p2 c0 {2,S} {4,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.67053,0.00492586,8.86967e-05,-1.26219e-07,5.23991e-11,-10278.8,14.5722], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.60658,0.020831,-8.42229e-06,1.56718e-09,-1.09391e-13,-13239.3,-23.2465], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""T 3/97""",
    longDesc = 
u"""
T 3/97.
C1=COCC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 384,
    label = "C5H92-1,3OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {1,S} {3,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.42923,0.0733174,-5.05533e-05,1.77709e-08,-2.53421e-12,-18921.3,22.6144], Tmin=(298,'K'), Tmax=(1413,'K')),
            NASAPolynomial(coeffs=[24.6778,0.0247897,-8.53168e-06,1.33018e-09,-7.7405e-14,-26400.9,-91.8794], Tmin=(1413,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC([CH]COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 385,
    label = "CCYC2OCO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {2,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  O u1 p2 c0 {4,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.97381,0.0675615,-5.79902e-05,2.55213e-08,-4.49992e-12,-9936.07,42.5128], Tmin=(298,'K'), Tmax=(1401,'K')),
            NASAPolynomial(coeffs=[16.4353,0.0165843,-5.71791e-06,8.92587e-10,-5.19865e-14,-16093.8,-59.6946], Tmin=(1401,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""L 2/00""",
    longDesc = 
u"""
L 2/00
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(C[O])CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 386,
    label = "NC5KET15O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,D} {16,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u1 p2 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.41522,0.0472024,-1.66453e-05,-2.39549e-09,1.83815e-12,-22114.8,17.7035], Tmin=(298,'K'), Tmax=(1489,'K')),
            NASAPolynomial(coeffs=[19.2273,0.0211349,-7.32089e-06,1.14667e-09,-6.69508e-14,-28586.7,-70.8965], Tmin=(1489,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]CCCCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 387,
    label = "C5H92-4OOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {5,D} {16,S}
5  C u0 p0 c0 {3,S} {4,D} {15,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.10664,0.0770652,-6.01782e-05,2.37846e-08,-3.77299e-12,-16121.7,34.5219], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[22.2715,0.0199891,-6.7876e-06,1.05005e-09,-6.0813e-14,-23909.8,-90.0214], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CC(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 388,
    label = "C3H6OOH1-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {10,S} {11,S}
4  O u0 p2 c0 {2,S} {5,S}
5  O u0 p2 c0 {4,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.74271,0.0453734,-3.5758e-05,1.4854e-08,-2.49982e-12,232.581,22.0973], Tmin=(298,'K'), Tmax=(1401,'K')),
            NASAPolynomial(coeffs=[13.9131,0.0140218,-4.55921e-06,6.84182e-10,-3.87696e-14,-3656.51,-42.1533], Tmin=(1401,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 389,
    label = "NC5KET31",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {2,S} {17,D}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 O u0 p2 c0 {5,D}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.06394,0.0546275,-2.50666e-05,2.02313e-09,1.0015e-12,-44162.5,15.9655], Tmin=(298,'K'), Tmax=(1431,'K')),
            NASAPolynomial(coeffs=[21.4963,0.0220659,-7.3549e-06,1.12359e-09,-6.45028e-14,-50847.1,-80.243], Tmin=(1431,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(=O)CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 390,
    label = "NC5KET24",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {4,S} {17,D}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 O u0 p2 c0 {5,D}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.70378,0.062214,-4.08904e-05,1.31907e-08,-1.60742e-12,-46241.1,15.8017], Tmin=(298,'K'), Tmax=(1417,'K')),
            NASAPolynomial(coeffs=[21.3558,0.0220059,-7.28286e-06,1.10611e-09,-6.32057e-14,-52344.6,-79.1469], Tmin=(1417,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)CC(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 391,
    label = "CH2CH2COCH3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {10,D}
4  C u1 p0 c0 {1,S} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 O u0 p2 c0 {3,D}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.36192,0.0350401,-1.70487e-05,3.0907e-09,2.80365e-14,-6027.54,19.5185], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[11.7916,0.0170296,-5.75444e-06,8.86035e-10,-5.11071e-14,-9706.83,-32.5533], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 392,
    label = "IC4H7OOIC4H7",
    molecule = 
"""
1  C u0 p0 c0 {5,S} {9,S} {11,S} {12,S}
2  C u0 p0 c0 {6,S} {10,S} {13,S} {14,S}
3  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
4  C u0 p0 c0 {6,S} {18,S} {19,S} {20,S}
5  C u0 p0 c0 {1,S} {3,S} {7,D}
6  C u0 p0 c0 {2,S} {4,S} {8,D}
7  C u0 p0 c0 {5,D} {21,S} {22,S}
8  C u0 p0 c0 {6,D} {23,S} {24,S}
9  O u0 p2 c0 {1,S} {10,S}
10 O u0 p2 c0 {2,S} {9,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.207882,0.104076,-8.33973e-05,3.60898e-08,-6.47907e-12,-7790.63,35.0447], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[28.0447,0.0327505,-1.1322e-05,1.77028e-09,-1.03212e-13,-17268.3,-115.135], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)COOCC(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 393,
    label = "C5H5O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {1,S} {5,D} {9,S}
4  C u0 p0 c0 {2,D} {5,S} {10,S}
5  C u0 p0 c0 {3,D} {4,S} {11,S}
6  O u1 p2 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.83113,0.0567277,-4.44757e-05,1.74924e-08,-2.76005e-12,20499.2,36.9634], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[14.8323,0.0140483,-4.92302e-06,7.77041e-10,-4.56104e-14,14552.4,-57.3228], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""5/16/90 THERM""",
    longDesc = 
u"""
5/16/90 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]C1C=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 394,
    label = "CH3COOH",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,S} {7,D}
3 O u0 p2 c0 {2,S} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 O u0 p2 c0 {2,D}
8 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.7895,0.00999942,3.42572e-05,-5.09031e-08,2.06222e-11,-53475.2,14.1053], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.67085,0.0135153,-5.25874e-06,8.93184e-10,-5.53181e-14,-55756.1,-15.4677], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""G 6/00""",
    longDesc = 
u"""
G 6/00.
CC(=O)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 395,
    label = "NC5CYCPER24",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {2,S} {6,S}
8  O u0 p2 c0 {1,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.2587,0.0926924,-7.79864e-05,3.34971e-08,-5.74391e-12,-47690.2,46.2164], Tmin=(298,'K'), Tmax=(1408,'K')),
            NASAPolynomial(coeffs=[22.2821,0.023331,-7.91636e-06,1.22269e-09,-7.06907e-14,-56069.8,-93.5419], Tmin=(1408,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""11/15 THERM""",
    longDesc = 
u"""
11/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1CC(C)(O)OO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 396,
    label = "C5H4O",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {6,D}
2  C u0 p0 c0 {1,S} {3,D} {7,S}
3  C u0 p0 c0 {2,D} {4,S} {8,S}
4  C u0 p0 c0 {3,S} {5,D} {9,S}
5  C u0 p0 c0 {1,S} {4,D} {10,S}
6  O u0 p2 c0 {1,D}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.264576,0.0334874,1.67738e-06,-2.96207e-08,1.54431e-11,5111.59,23.541], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[10.0807,0.0161143,-5.83315e-06,9.46759e-10,-5.68972e-14,1943.65,-29.4522], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 8/99""",
    longDesc = 
u"""
T 8/99.
O=C1C=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 397,
    label = "C5H9O1-2OOH-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {2,S} {4,S}
7  O u0 p2 c0 {1,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.63459,0.0837429,-6.46449e-05,2.45435e-08,-3.63638e-12,-30244.4,38.101], Tmin=(298,'K'), Tmax=(1422,'K')),
            NASAPolynomial(coeffs=[24.0742,0.021005,-7.12966e-06,1.10196e-09,-6.37577e-14,-38690.6,-98.6615], Tmin=(1422,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(OO)C1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 398,
    label = "C5H9O2-3OOH-1",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {2,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.95212,0.0890845,-7.32339e-05,2.98425e-08,-4.76173e-12,-30028.1,44.5335], Tmin=(298,'K'), Tmax=(1436,'K')),
            NASAPolynomial(coeffs=[23.731,0.0206752,-6.88444e-06,1.05045e-09,-6.02362e-14,-38436.3,-96.1946], Tmin=(1436,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1OC1COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 399,
    label = "C4H8OOH2-1O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {3,S} {16,S}
7  O u0 p2 c0 {5,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u1 p2 c0 {6,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.02241,0.0653863,-4.89646e-05,1.90438e-08,-3.02309e-12,-22580.6,20.5729], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[21.5735,0.0204529,-6.99498e-06,1.08598e-09,-6.30071e-14,-28842.8,-78.4561], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(CO[O])OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 400,
    label = "AC5H11O2H",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {19,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.476375,0.0783789,-5.91713e-05,2.38264e-08,-3.97732e-12,-30508.8,33.2978], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[20.8513,0.0260811,-8.91736e-06,1.38397e-09,-8.02717e-14,-37719.4,-80.4522], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 401,
    label = "C5H91-4OOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {5,D} {14,S}
5  C u0 p0 c0 {4,D} {15,S} {16,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.38243,0.0667661,-4.96344e-05,1.94278e-08,-3.12176e-12,-16908.7,24.0816], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[19.6307,0.0222,-7.52489e-06,1.16109e-09,-6.70722e-14,-23040,-73.2135], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCC(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 402,
    label = "C5H93-1,4OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {6,S} {10,S}
2  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {7,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {1,S} {2,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.42923,0.0733174,-5.05533e-05,1.77709e-08,-2.53421e-12,-18921.3,22.6144], Tmin=(298,'K'), Tmax=(1413,'K')),
            NASAPolynomial(coeffs=[24.6778,0.0247897,-8.53168e-06,1.33018e-09,-7.7405e-14,-26400.9,-91.8794], Tmin=(1413,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC([CH]CCOO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 403,
    label = "C4H72-3,4OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {4,S} {12,S} {13,S} {14,S}
4  C u1 p0 c0 {1,S} {3,S} {15,S}
5  O u0 p2 c0 {1,S} {8,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {16,S}
8  O u0 p2 c0 {5,S} {17,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.85542,0.0605524,-4.18122e-05,1.4668e-08,-2.07882e-12,-16035.9,18.7735], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[21.4626,0.0202946,-6.97889e-06,1.08751e-09,-6.32605e-14,-22219.7,-76.0672], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C(COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 404,
    label = "C5H10OOH1-4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {2,S} {4,S} {17,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.21642,0.0596427,-3.30151e-05,7.82802e-09,-4.53627e-13,-6890.38,25.3166], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[19.4584,0.0246032,-8.38535e-06,1.29884e-09,-7.52367e-14,-13354.9,-69.0616], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 405,
    label = "A-BC5H10O",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
6  O u0 p2 c0 {1,S} {3,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.81578,0.070415,-5.23594e-05,1.99929e-08,-3.0557e-12,-19908,38.8311], Tmin=(298,'K'), Tmax=(1424,'K')),
            NASAPolynomial(coeffs=[16.6888,0.0223673,-7.44829e-06,1.13558e-09,-6.50517e-14,-26312.7,-64.8349], Tmin=(1424,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1(C)CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 406,
    label = "C3H6OH2-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u1 p0 c0 {1,S} {9,S} {10,S}
4  O u0 p2 c0 {1,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.0967,0.0380728,-2.75022e-05,1.07477e-08,-1.74896e-12,-14076.4,22.2476], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[11.2222,0.0136444,-4.51407e-06,7.10523e-10,-4.2269e-14,-17535,-31.8912], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/ 9/ 4 THERM""",
    longDesc = 
u"""
8/ 9/ 4 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 407,
    label = "C4H8O1-2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.29657,0.0675907,-5.89614e-05,2.59401e-08,-4.45927e-12,-14880,44.7756], Tmin=(298,'K'), Tmax=(1463,'K')),
            NASAPolynomial(coeffs=[14.1886,0.0163163,-5.16581e-06,7.60174e-10,-4.24548e-14,-20283.9,-51.2818], Tmin=(1463,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 408,
    label = "C5H91-3,5OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {7,S} {15,S} {16,S}
5  C u1 p0 c0 {3,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.50133,0.0793141,-5.51609e-05,1.83964e-08,-2.32881e-12,-17513,26.433], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[27.5716,0.0226282,-7.85652e-06,1.23252e-09,-7.20427e-14,-26270.3,-108.648], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(CCOO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 409,
    label = "C5H10O2-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {11,S} {12,S} {13,S}
6  O u0 p2 c0 {1,S} {2,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-6.83688,0.0900167,-8.09486e-05,3.71523e-08,-6.73312e-12,-19490.7,57.3459], Tmin=(298,'K'), Tmax=(1415,'K')),
            NASAPolynomial(coeffs=[17.6695,0.0216741,-7.23761e-06,1.10546e-09,-6.34044e-14,-26800,-70.1762], Tmin=(1415,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1OC1C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 410,
    label = "C5H10O1-4",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
6  O u0 p2 c0 {1,S} {4,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.99511,0.0848028,-6.84821e-05,2.79647e-08,-4.48602e-12,-27934.7,63.6835], Tmin=(298,'K'), Tmax=(1470,'K')),
            NASAPolynomial(coeffs=[15.6838,0.0228597,-7.36869e-06,1.09856e-09,-6.19426e-14,-35236.5,-60.6987], Tmin=(1470,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1CCCO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 411,
    label = "NC5KET13",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,D} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.00527,0.0657859,-4.36349e-05,1.39803e-08,-1.68002e-12,-42409.5,19.0078], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[22.5902,0.0217812,-7.38064e-06,1.13911e-09,-6.58296e-14,-49267.9,-86.5971], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(CC=O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 412,
    label = "TQC4H7OHI",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {6,S} {14,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {4,S} {15,S}
7  O u0 p2 c0 {5,S} {16,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.55844,0.0637086,-4.9717e-05,1.99225e-08,-3.21373e-12,-27752.7,19.5001], Tmin=(298,'K'), Tmax=(1404,'K')),
            NASAPolynomial(coeffs=[20.8281,0.0181675,-6.12943e-06,9.43194e-10,-5.43937e-14,-33738.7,-77.4824], Tmin=(1404,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""L 2/00""",
    longDesc = 
u"""
L 2/00
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)([CH]O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 413,
    label = "C4H8OOH2-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {3,S} {14,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.3029,0.0508232,-3.39615e-05,1.17622e-08,-1.66013e-12,-6136.71,15.1726], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[17.0354,0.019273,-6.50685e-06,1.00126e-09,-5.77271e-14,-10943.7,-58.7403], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 414,
    label = "CH2CH2COC2H5",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {2,S} {13,D}
5  C u1 p0 c0 {2,S} {14,S} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.36988,0.0393334,-1.08692e-05,-4.28049e-09,2.08636e-12,-9337.95,16.2476], Tmin=(298,'K'), Tmax=(1418,'K')),
            NASAPolynomial(coeffs=[15.803,0.0200451,-6.64102e-06,1.01015e-09,-5.78092e-14,-14484,-53.7408], Tmin=(1418,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/6/15""",
    longDesc = 
u"""
10/6/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(=O)CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 415,
    label = "IC4H7CHO",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {5,D}
4  C u0 p0 c0 {1,S} {11,D} {12,S}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 O u0 p2 c0 {4,D}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.20983,0.0592603,-4.2696e-05,1.60356e-08,-2.49348e-12,-15620.5,35.3137], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[15.9172,0.0193357,-6.70858e-06,1.05155e-09,-6.14176e-14,-21614.1,-56.7617], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 416,
    label = "C5H9O1-4OOH-2",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {6,S} {9,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {4,S}
7  O u0 p2 c0 {2,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-9.18085,0.108786,-9.86756e-05,4.49208e-08,-8.03088e-12,-39801.3,72.6761], Tmin=(298,'K'), Tmax=(1417,'K')),
            NASAPolynomial(coeffs=[21.726,0.0234146,-7.86188e-06,1.20572e-09,-6.93654e-14,-49050.9,-88.3591], Tmin=(1417,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1CC(CO1)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 417,
    label = "CC5H10OOH-D",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {2,S} {16,S} {17,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.913422,0.0751514,-6.12063e-05,2.64606e-08,-4.65309e-12,-8129.47,26.3698], Tmin=(298,'K'), Tmax=(1407,'K')),
            NASAPolynomial(coeffs=[20.7407,0.0231707,-7.81226e-06,1.20092e-09,-6.91855e-14,-14455.5,-78.1145], Tmin=(1407,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(OO)C(C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 418,
    label = "C5H9O1-2OOH-4",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {6,S} {10,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {4,S}
7  O u0 p2 c0 {2,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.63459,0.0837429,-6.46449e-05,2.45435e-08,-3.63638e-12,-30244.4,38.101], Tmin=(298,'K'), Tmax=(1422,'K')),
            NASAPolynomial(coeffs=[24.0742,0.021005,-7.12966e-06,1.10196e-09,-6.37577e-14,-38690.6,-98.6615], Tmin=(1422,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CC1CO1)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 419,
    label = "B2E2M1OJ",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
2  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {5,S} {7,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {3,S} {4,D} {15,S}
6  O u1 p2 c0 {2,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.64763,0.034361,-5.6065e-06,-5.19111e-09,1.73155e-12,-6.65738,3.13825], Tmin=(298,'K'), Tmax=(2003,'K')),
            NASAPolynomial(coeffs=[13.5666,0.0256548,-9.37226e-06,1.52518e-09,-9.15379e-14,-3587.67,-42.7273], Tmin=(2003,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C(C)C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 420,
    label = "C3H51-2,3OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u1 p0 c0 {1,S} {11,S} {12,S}
4  O u0 p2 c0 {1,S} {7,S}
5  O u0 p2 c0 {2,S} {6,S}
6  O u0 p2 c0 {5,S} {13,S}
7  O u0 p2 c0 {4,S} {14,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5562,0.0613504,-5.23205e-05,2.28208e-08,-4.02232e-12,-13135.3,22.1044], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[21.2378,0.013952,-4.94539e-06,7.86381e-10,-4.63926e-14,-19286.5,-76.9637], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/26/3 THRM""",
    longDesc = 
u"""
8/26/3 THRM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 421,
    label = "CHOCHO",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,D} {5,S}
2 C u0 p0 c0 {1,S} {4,D} {6,S}
3 O u0 p2 c0 {1,D}
4 O u0 p2 c0 {2,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.88105,0.0236386,-1.83443e-05,6.84843e-09,-9.92734e-13,-26928,15.9155], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[9.75439,0.00497646,-1.7441e-06,2.75587e-10,-1.6197e-14,-29583.3,-26.1878], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 422,
    label = "C6H101-3,3",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {4,D} {5,S}
4  C u0 p0 c0 {2,S} {3,D} {13,S}
5  C u0 p0 c0 {3,S} {6,D} {14,S}
6  C u0 p0 c0 {5,D} {15,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.030131,0.0643105,-4.75083e-05,1.89878e-08,-3.17681e-12,2823.67,24.9255], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[16.9678,0.0228868,-7.78758e-06,1.20465e-09,-6.97081e-14,-2972.77,-65.8633], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 423,
    label = "C6H10D24",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,D} {13,S}
4  C u0 p0 c0 {2,S} {6,D} {14,S}
5  C u0 p0 c0 {3,D} {6,S} {15,S}
6  C u0 p0 c0 {4,D} {5,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.363783,0.0636719,-4.5309e-05,1.73368e-08,-2.79133e-12,3040.87,27.5469], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[16.7663,0.023211,-7.93355e-06,1.23102e-09,-7.13893e-14,-2933.25,-64.4104], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CC=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 424,
    label = "B12DE3M",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {2,S} {5,D}
4  C u0 p0 c0 {5,D} {12,S} {13,S}
5  C u0 p0 c0 {3,D} {4,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.2786,0.0454917,-2.82334e-05,8.98449e-09,-1.1722e-12,13163.7,18.7039], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[13.7093,0.0185727,-6.33116e-06,9.80719e-10,-5.681e-14,8595.19,-48.875], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""11/12/12 THERM""",
    longDesc = 
u"""
11/12/12 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C=C(C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 425,
    label = "NC5H11OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {4,S} {18,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.740083,0.068388,-4.47177e-05,1.49115e-08,-1.99982e-12,-38250.1,34.1516], Tmin=(298,'K'), Tmax=(1399,'K')),
            NASAPolynomial(coeffs=[18.1326,0.0258734,-8.76148e-06,1.35108e-09,-7.8017e-14,-44939.4,-67.7259], Tmin=(1399,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 7/15 THERM""",
    longDesc = 
u"""
10/ 7/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 426,
    label = "HOCHCHO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 O u0 p2 c0 {1,S} {7,S}
4 H u0 p0 c0 {1,S}
5 O u1 p2 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.542856,0.0295296,-2.08315e-05,6.72777e-09,-8.06116e-13,-18937.8,23.1681], Tmin=(298,'K'), Tmax=(1680,'K')),
            NASAPolynomial(coeffs=[9.99252,0.00777561,-2.94239e-06,4.90136e-10,-2.9898e-14,-22044.1,-27.3957], Tmin=(1680,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/25/15""",
    longDesc = 
u"""
9/25/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]C=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 427,
    label = "C5H92-1,5OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {6,S} {14,S} {15,S}
4  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
5  C u1 p0 c0 {2,S} {4,S} {18,S}
6  O u0 p2 c0 {3,S} {9,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.43107,0.0725134,-4.53395e-05,1.29499e-08,-1.26535e-12,-16900.1,23.7844], Tmin=(298,'K'), Tmax=(1399,'K')),
            NASAPolynomial(coeffs=[26.3836,0.0238373,-8.31679e-06,1.30878e-09,-7.66591e-14,-25293.1,-101.119], Tmin=(1399,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOC[CH]CCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 428,
    label = "C5H9O1-4OOH-5",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {6,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
6  O u0 p2 c0 {1,S} {4,S}
7  O u0 p2 c0 {5,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-9.29661,0.104959,-9.00893e-05,3.89844e-08,-6.6882e-12,-37638.3,75.2774], Tmin=(298,'K'), Tmax=(1411,'K')),
            NASAPolynomial(coeffs=[21.3868,0.0242467,-8.26745e-06,1.28123e-09,-7.42528e-14,-47241.1,-86.0552], Tmin=(1411,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOCC1CCCO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 429,
    label = "C4H7O1-3OOH-2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {3,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {15,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.68244,0.080342,-6.67265e-05,2.74081e-08,-4.41624e-12,-27530.7,50.4193], Tmin=(298,'K'), Tmax=(1425,'K')),
            NASAPolynomial(coeffs=[19.711,0.017906,-6.05725e-06,9.3411e-10,-5.39628e-14,-35253.3,-78.3129], Tmin=(1425,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1OCC1OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 430,
    label = "C3H6OOH2-1O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  O u0 p2 c0 {1,S} {6,S}
5  O u0 p2 c0 {2,S} {13,S}
6  O u0 p2 c0 {4,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u1 p2 c0 {5,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.99085,0.0531865,-4.28598e-05,1.77187e-08,-2.92769e-12,-20214.4,13.4151], Tmin=(298,'K'), Tmax=(1404,'K')),
            NASAPolynomial(coeffs=[19.1045,0.0144076,-4.72128e-06,7.12632e-10,-4.05578e-14,-25027.1,-66.3748], Tmin=(1404,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CO[O])OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 431,
    label = "C5H93-2,4OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {7,S} {10,S}
2  C u0 p0 c0 {4,S} {5,S} {6,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {1,S} {2,S} {18,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {1,S} {9,S}
8  O u0 p2 c0 {6,S} {19,S}
9  O u0 p2 c0 {7,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.76788,0.0747668,-5.4178e-05,2.01841e-08,-3.04868e-12,-21082.3,18.5667], Tmin=(298,'K'), Tmax=(1418,'K')),
            NASAPolynomial(coeffs=[25.107,0.0240686,-8.20467e-06,1.27103e-09,-7.36339e-14,-28359.1,-95.6533], Tmin=(1418,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC([CH]C(C)OO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 432,
    label = "C5H3O",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {4,S} {6,S}
2 C u0 p0 c0 {1,D} {3,S} {7,S}
3 C u0 p0 c0 {2,S} {5,S} {8,D}
4 C u0 p0 c0 {1,S} {5,D} {9,S}
5 C u1 p0 c0 {3,S} {4,D}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 O u0 p2 c0 {3,D}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.03243,0.0543937,-4.95018e-05,2.25524e-08,-4.10728e-12,33564.4,37.8375], Tmin=(298,'K'), Tmax=(1500,'K')),
            NASAPolynomial(coeffs=[11.9962,0.0134287,-5.90045e-06,1.22554e-09,-9.86115e-14,28959.2,-40.7548], Tmin=(1500,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""TAK0905""",
    longDesc = 
u"""
TAK0905
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=C1[C]=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 433,
    label = "C4H7O2-3OOH-1",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {8,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {2,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {15,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.04171,0.0775255,-6.55777e-05,2.74944e-08,-4.52477e-12,-27034.3,43.3137], Tmin=(298,'K'), Tmax=(1424,'K')),
            NASAPolynomial(coeffs=[20.3028,0.0169332,-5.71129e-06,8.79006e-10,-5.07088e-14,-34345.1,-79.5919], Tmin=(1424,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1OC1COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 434,
    label = "C5H10OH15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u1 p0 c0 {3,S} {15,S} {16,S}
6  O u0 p2 c0 {4,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.394307,0.0650012,-4.38043e-05,1.5254e-08,-2.16516e-12,-13525.9,34.8138], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[17.4869,0.0239843,-8.13298e-06,1.25529e-09,-7.25301e-14,-19796.9,-61.4612], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 7/15 THERM""",
    longDesc = 
u"""
10/ 7/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 435,
    label = "NC5KET25O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {4,S} {16,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 O u1 p2 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.40785,0.0434772,-1.40285e-05,-3.0231e-09,1.86769e-12,-26026.5,12.8766], Tmin=(298,'K'), Tmax=(1433,'K')),
            NASAPolynomial(coeffs=[18.1304,0.0213611,-7.24506e-06,1.11897e-09,-6.47006e-14,-31668.3,-64.1505], Tmin=(1433,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)CCC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 436,
    label = "C5H9O1-2OOH-5",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {6,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {7,S} {14,S} {15,S}
6  O u0 p2 c0 {1,S} {4,S}
7  O u0 p2 c0 {5,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.91503,0.086372,-6.80859e-05,2.66151e-08,-4.08265e-12,-27929.8,45.8804], Tmin=(298,'K'), Tmax=(1425,'K')),
            NASAPolynomial(coeffs=[23.1522,0.0215717,-7.27788e-06,1.12022e-09,-6.46249e-14,-36375.4,-92.3673], Tmin=(1425,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOCCCC1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 437,
    label = "C5H9O1-3OOH-5",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {6,S} {16,S} {17,S}
5  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
6  O u0 p2 c0 {1,S} {4,S}
7  O u0 p2 c0 {5,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.55071,0.0929081,-7.49705e-05,3.02211e-08,-4.80853e-12,-28256.1,57.9511], Tmin=(298,'K'), Tmax=(1421,'K')),
            NASAPolynomial(coeffs=[22.1703,0.022876,-7.73631e-06,1.19254e-09,-6.8863e-14,-37150.3,-88.7114], Tmin=(1421,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOCCC1CCO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 438,
    label = "C3H6OOH1-2O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  O u0 p2 c0 {2,S} {6,S}
5  O u0 p2 c0 {1,S} {13,S}
6  O u0 p2 c0 {4,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u1 p2 c0 {5,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.99085,0.0531865,-4.28598e-05,1.77187e-08,-2.92769e-12,-20214.4,13.4151], Tmin=(298,'K'), Tmax=(1404,'K')),
            NASAPolynomial(coeffs=[19.1045,0.0144076,-4.72128e-06,7.12632e-10,-4.05578e-14,-25027.1,-66.3748], Tmin=(1404,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(COO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 439,
    label = "C3H5O1-3OOH-2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
4  O u0 p2 c0 {2,S} {3,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.43959,0.0716533,-7.15032e-05,3.51738e-08,-6.63683e-12,-22725,45.6394], Tmin=(298,'K'), Tmax=(1434,'K')),
            NASAPolynomial(coeffs=[14.4493,0.0136373,-4.25837e-06,6.20006e-10,-3.43452e-14,-27737.2,-50.6099], Tmin=(1434,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/13 THER""",
    longDesc = 
u"""
10/13 THER
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOC1COC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 440,
    label = "TC3H6O2HCO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  O u0 p2 c0 {1,S} {6,S}
5  C u1 p0 c0 {1,S} {13,D}
6  O u0 p2 c0 {4,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {5,D}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.03864,0.0580421,-4.32124e-05,1.58792e-08,-2.3221e-12,-22428.5,20.3681], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[20.6473,0.0148527,-5.25105e-06,8.33619e-10,-4.91256e-14,-28872,-79.5951], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/ 2/95 THERM""",
    longDesc = 
u"""
8/ 2/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)([C]=O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 441,
    label = "C-DC5H10O",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
3  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
6  O u0 p2 c0 {2,S} {3,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.94335,0.093635,-8.55714e-05,3.97013e-08,-7.24457e-12,-17911.7,62.3627], Tmin=(298,'K'), Tmax=(1415,'K')),
            NASAPolynomial(coeffs=[17.5722,0.0216806,-7.22379e-06,1.1017e-09,-6.31222e-14,-25433.7,-70.1132], Tmin=(1415,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)C1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 442,
    label = "NC5KET25",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {4,S} {17,D}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 O u0 p2 c0 {5,D}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.30376,0.0601763,-3.55337e-05,9.4996e-09,-7.96339e-13,-44048.6,19.6187], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[20.9052,0.0229956,-7.75343e-06,1.19263e-09,-6.87613e-14,-50467.8,-76.1624], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)CCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 443,
    label = "C5H9O1-5OOH-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {6,S} {16,S} {17,S}
6  O u0 p2 c0 {4,S} {5,S}
7  O u0 p2 c0 {1,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-9.63082,0.10654,-8.95954e-05,3.77712e-08,-6.32249e-12,-39749.9,71.1221], Tmin=(298,'K'), Tmax=(1407,'K')),
            NASAPolynomial(coeffs=[22.4619,0.0242252,-8.36381e-06,1.30713e-09,-7.6201e-14,-50010.2,-98.3781], Tmin=(1407,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOC1CCOCC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 444,
    label = "TC3H6OCHO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,D} {13,S}
5  O u1 p2 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u0 p2 c0 {4,D}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.37083,0.0538476,-3.82478e-05,1.32882e-08,-1.79229e-12,-21839.1,25.8142], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[17.0371,0.0154401,-5.28333e-06,8.21085e-10,-4.76898e-14,-27587.2,-63.7271], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/25/95 THERM""",
    longDesc = 
u"""
8/25/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)([O])C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 445,
    label = "IIC4H7Q2-T",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
2  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
4  C u1 p0 c0 {1,S} {2,S} {3,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {5,S} {16,S}
8  O u0 p2 c0 {6,S} {17,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[8.16274,0.0434463,-1.76972e-05,4.88791e-10,9.03915e-13,-19650.2,0.262067], Tmin=(298,'K'), Tmax=(1377,'K')),
            NASAPolynomial(coeffs=[21.507,0.020536,-7.12383e-06,1.11655e-09,-6.52112e-14,-25111.8,-74.338], Tmin=(1377,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/15/96 THERM""",
    longDesc = 
u"""
7/15/96 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C](COO)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 446,
    label = "C5H9O2-4OOH-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {2,S} {3,S}
7  O u0 p2 c0 {1,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.69459,0.0949457,-7.98763e-05,3.34049e-08,-5.48568e-12,-32615.9,49.8789], Tmin=(298,'K'), Tmax=(1427,'K')),
            NASAPolynomial(coeffs=[23.6121,0.0213802,-7.16433e-06,1.09773e-09,-6.3128e-14,-41463.5,-99.1051], Tmin=(1427,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1OC(C)C1OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 447,
    label = "AC5H10OOH-C",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {1,S} {4,S} {17,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.24167,0.0655076,-4.30153e-05,1.43717e-08,-1.92813e-12,-7438.58,28.852], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[19.427,0.0245167,-8.32902e-06,1.28717e-09,-7.4437e-14,-13880.9,-69.3043], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C(C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 448,
    label = "AC5H10OOH-A",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {1,S} {16,S} {17,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.0893968,0.0748083,-5.79088e-05,2.38676e-08,-4.05836e-12,-5791.21,33.772], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[20.2759,0.0240402,-8.21726e-06,1.27509e-09,-7.3947e-14,-12589,-74.5433], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(CC)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 449,
    label = "SC4H8OH-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {3,S} {13,S}
5  O u0 p2 c0 {1,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.93186,0.0479506,-3.23719e-05,1.16233e-08,-1.72499e-12,-13943.1,20.3995], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[14.2454,0.0189491,-6.27404e-06,9.52693e-10,-5.44158e-14,-18180.4,-45.6206], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C(C)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 450,
    label = "C5H10O1-5",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
6  O u0 p2 c0 {4,S} {5,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.45067,0.0836042,-6.24375e-05,2.32498e-08,-3.38258e-12,-27906,61.9634], Tmin=(298,'K'), Tmax=(1447,'K')),
            NASAPolynomial(coeffs=[15.8723,0.0241485,-8.02831e-06,1.22306e-09,-7.00385e-14,-35866.7,-67.3549], Tmin=(1447,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C1CCOCC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 451,
    label = "IC5KETCDO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {4,S} {16,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u1 p2 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[7.66114,0.0331213,-2.76253e-06,-8.32131e-09,2.78276e-12,-26882.3,-4.11947], Tmin=(298,'K'), Tmax=(1489,'K')),
            NASAPolynomial(coeffs=[17.8986,0.0213374,-7.19313e-06,1.10668e-09,-6.38267e-14,-31606.4,-63.332], Tmin=(1489,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)C(=O)C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 452,
    label = "TC4H8O2H-I",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.54378,0.0525201,-3.69898e-05,1.44635e-08,-2.47536e-12,-6981.83,12.7624], Tmin=(298,'K'), Tmax=(1379,'K')),
            NASAPolynomial(coeffs=[18.1415,0.0194699,-6.8275e-06,1.07773e-09,-6.32519e-14,-12357.1,-66.3492], Tmin=(1379,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 453,
    label = "IQC4H8OT",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {2,S} {6,S}
6  O u0 p2 c0 {5,S} {16,S}
7  O u1 p2 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.72212,0.0642865,-5.52809e-05,2.50037e-08,-4.53472e-12,-24311.1,12.1981], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[20.4824,0.0182967,-6.04413e-06,9.16381e-10,-5.22867e-14,-29428.7,-75.3563], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)([O])COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 454,
    label = "NC5KET12",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {12,S}
3  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
4  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {16,D} {17,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.29175,0.0864139,-7.02456e-05,2.86657e-08,-4.68322e-12,-40171.9,40.0248], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[25.6959,0.0197582,-6.83009e-06,1.06938e-09,-6.24522e-14,-49102.8,-103.505], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC(C=O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 455,
    label = "CC5H10OOH-A",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {1,S} {16,S} {17,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.328494,0.0768086,-6.328e-05,2.76118e-08,-4.8869e-12,-7987.41,29.8696], Tmin=(298,'K'), Tmax=(1407,'K')),
            NASAPolynomial(coeffs=[20.5911,0.0232494,-7.8291e-06,1.20251e-09,-6.92367e-14,-14404.7,-76.7477], Tmin=(1407,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 456,
    label = "NC5KET21O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {4,S} {16,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u1 p2 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[7.57725,0.0309672,2.94906e-06,-1.1097e-08,3.01675e-12,-25634.1,-2.1672], Tmin=(298,'K'), Tmax=(1669,'K')),
            NASAPolynomial(coeffs=[15.5374,0.0265967,-9.8252e-06,1.61087e-09,-9.71864e-14,-29622,-49.7983], Tmin=(1669,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC(=O)C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 457,
    label = "NC53ONEO2-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {2,S} {16,D}
6  O u0 p2 c0 {3,S} {17,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.51723,0.0463562,-1.26144e-05,-5.62684e-09,2.63012e-12,-27332.2,15.1182], Tmin=(298,'K'), Tmax=(1472,'K')),
            NASAPolynomial(coeffs=[20.8298,0.0214529,-7.34745e-06,1.14331e-09,-6.64861e-14,-34182.9,-76.9527], Tmin=(1472,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/8/15""",
    longDesc = 
u"""
10/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(=O)CCO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 458,
    label = "IC5KETCD",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {4,S} {17,D}
6  O u0 p2 c0 {4,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 O u0 p2 c0 {5,D}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.283544,0.0806299,-6.46192e-05,2.62468e-08,-4.26141e-12,-43048.6,35.1464], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[23.7899,0.0204032,-6.84168e-06,1.04955e-09,-6.04314e-14,-50895.7,-92.5391], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)C(=O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 459,
    label = "C2H4O1-2",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
3 O u0 p2 c0 {1,S} {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.75905,-0.00944122,8.03097e-05,-1.00808e-07,4.00399e-11,-7560.81,7.84975], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.48876,0.0120462,-4.33369e-06,7.00283e-10,-4.19491e-14,-9180.43,-7.07996], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 8/88""",
    longDesc = 
u"""
L 8/88.
C1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 460,
    label = "NC52ONEO2-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {4,S} {16,D}
6  O u0 p2 c0 {1,S} {17,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.367213,0.0761783,-6.13905e-05,2.51923e-08,-4.15674e-12,-27349.1,32.8428], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[23.2385,0.0191384,-6.48796e-06,1.00265e-09,-5.80286e-14,-34864.4,-88.6015], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/8/15""",
    longDesc = 
u"""
10/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(O[O])C(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 461,
    label = "IC3H6CO",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {11,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 O u0 p2 c0 {4,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.28039,0.0417017,-3.2509e-05,1.37243e-08,-2.40573e-12,-16394,13.8188], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[13.2548,0.0140143,-4.7891e-06,7.42924e-10,-4.30738e-14,-20053,-44.481], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""03/03/95 THERM""",
    longDesc = 
u"""
03/03/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)=C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 462,
    label = "DC5H10OOH-A",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {1,S} {16,S} {17,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.0893968,0.0748083,-5.79088e-05,2.38676e-08,-4.05836e-12,-5791.21,33.772], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[20.2759,0.0240402,-8.21726e-06,1.27509e-09,-7.3947e-14,-12589,-74.5433], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 463,
    label = "C4H8O1-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.37284,0.0662224,-5.35318e-05,2.19452e-08,-3.5448e-12,-15414.5,49.8345], Tmin=(298,'K'), Tmax=(1447,'K')),
            NASAPolynomial(coeffs=[13.2077,0.0177468,-5.69934e-06,8.47771e-10,-4.77346e-14,-21171.8,-47.8386], Tmin=(1447,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1CCO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 464,
    label = "C2CY(COC)OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {2,S}
6  O u0 p2 c0 {2,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.592324,0.0552429,-4.02419e-05,1.57152e-08,-2.57388e-12,-35824.1,22.3378], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[15.683,0.0192911,-6.63718e-06,1.03441e-09,-6.01715e-14,-41059.8,-58.5686], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(C)OC1O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 465,
    label = "C4H8O2-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {7,S}
3  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {2,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.48183,0.0689313,-6.15372e-05,2.73744e-08,-4.85891e-12,-16892.4,43.8883], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[10.6342,0.0241442,-1.13124e-05,2.25481e-09,-1.54043e-13,-21038.3,-33.6764], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1OC1C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 466,
    label = "PC4H8OH-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {2,S} {5,S} {13,S}
5  O u0 p2 c0 {4,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.772567,0.0500429,-3.20164e-05,1.02913e-08,-1.30709e-12,-13842.2,25.8928], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[14.7813,0.0189693,-6.38671e-06,9.81341e-10,-5.65332e-14,-18848.7,-49.8892], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[CH]O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 467,
    label = "C3KET12",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,D} {11,S}
4  O u0 p2 c0 {1,S} {5,S}
5  O u0 p2 c0 {4,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 O u0 p2 c0 {3,D}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.03883,0.053418,-4.47684e-05,1.94652e-08,-3.45055e-12,-37030.9,25.6511], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[17.0188,0.0132097,-4.67055e-06,7.41412e-10,-4.3687e-14,-42357.3,-59.2616], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/17/12""",
    longDesc = 
u"""
10/17/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C=O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 468,
    label = "A-DC5H10O",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
6  O u0 p2 c0 {3,S} {4,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.99709,0.0881388,-7.28826e-05,3.06897e-08,-5.10686e-12,-26779.2,68.8674], Tmin=(298,'K'), Tmax=(1453,'K')),
            NASAPolynomial(coeffs=[15.3605,0.02351,-7.66793e-06,1.15246e-09,-6.53534e-14,-34248.8,-58.8396], Tmin=(1453,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1CCOC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 469,
    label = "TIC4H7Q2-I",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {14,S} {15,S}
5  O u0 p2 c0 {1,S} {8,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {16,S}
8  O u0 p2 c0 {5,S} {17,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.48426,0.0661225,-5.27349e-05,2.18216e-08,-3.66789e-12,-19890.7,12.672], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[23.3849,0.018707,-6.44022e-06,1.00428e-09,-5.84468e-14,-26118.1,-87.661], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""5/ 6/96 THERM""",
    longDesc = 
u"""
5/ 6/96 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)(COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 470,
    label = "C4H6O2-1OOH4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,D} {13,S}
4  C u0 p0 c0 {2,S} {3,D} {12,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  O u1 p2 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.38345,0.044709,-2.4363e-05,5.07835e-09,-1.20806e-13,-3407.11,8.84731], Tmin=(298,'K'), Tmax=(1364,'K')),
            NASAPolynomial(coeffs=[20.5165,0.0155728,-5.54197e-06,8.83593e-10,-5.22242e-14,-9300.31,-74.6631], Tmin=(1364,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]CC=CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 471,
    label = "NEOC5H11O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {6,S} {16,S} {17,S}
6  O u1 p2 c0 {5,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[11.14,0.00995551,4.37547e-05,-3.60664e-08,8.13723e-12,-15774,-25.4644], Tmin=(298,'K'), Tmax=(1374,'K')),
            NASAPolynomial(coeffs=[14.792,0.036024,-1.44549e-05,2.47285e-09,-1.52829e-13,-21288.1,-58.6606], Tmin=(1374,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(C)C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 472,
    label = "C5H10OOH2-4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {2,S} {4,S} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.69411,0.0613273,-3.77803e-05,1.1181e-08,-1.20044e-12,-9095,21.1389], Tmin=(298,'K'), Tmax=(1410,'K')),
            NASAPolynomial(coeffs=[19.8844,0.0237471,-7.98274e-06,1.22497e-09,-7.04947e-14,-15241.3,-71.9594], Tmin=(1410,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CC(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 473,
    label = "IC5KETCBO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {4,S} {16,D}
6  O u1 p2 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.60922,0.0463105,-1.85547e-05,-2.86159e-10,1.27532e-12,-28789.6,7.87702], Tmin=(298,'K'), Tmax=(1510,'K')),
            NASAPolynomial(coeffs=[19.3742,0.0206425,-7.06984e-06,1.09903e-09,-6.38342e-14,-34700.4,-74.3848], Tmin=(1510,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)C(C)(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 474,
    label = "C4H64,2-1OH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u1 p0 c0 {1,S} {3,S} {8,S}
3  C u0 p0 c0 {2,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {10,S} {11,S}
5  O u0 p2 c0 {1,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.50805,0.0353393,-1.41418e-05,5.65502e-10,5.94838e-13,-3085.74,18.1062], Tmin=(298,'K'), Tmax=(1366,'K')),
            NASAPolynomial(coeffs=[13.6259,0.0170233,-6.0046e-06,9.51604e-10,-5.60056e-14,-7814.99,-44.4948], Tmin=(1366,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[CH]CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 475,
    label = "IC5KETDCO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {15,D} {16,S}
6  H u0 p0 c0 {1,S}
7  O u1 p2 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.74668,0.0615642,-4.06489e-05,1.31251e-08,-1.65755e-12,-24025.3,22.6218], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[20.675,0.0201791,-7.04935e-06,1.11015e-09,-6.50552e-14,-30873.8,-80.0273], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)C([O])C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 476,
    label = "IC4H8O2H-I",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {13,S} {14,S}
5  O u0 p2 c0 {2,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.186433,0.062643,-4.83691e-05,1.88657e-08,-2.91189e-12,-3590.87,29.7635], Tmin=(298,'K'), Tmax=(1414,'K')),
            NASAPolynomial(coeffs=[18.3915,0.0173043,-5.66841e-06,8.55414e-10,-4.86782e-14,-9485.7,-66.7673], Tmin=(1414,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 477,
    label = "NC52ONEO2-5",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {4,S} {16,D}
6  O u0 p2 c0 {3,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.08985,0.0540453,-2.97971e-05,6.76826e-09,-2.84375e-13,-27362.5,16.4944], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[20.0577,0.0215456,-7.245e-06,1.11256e-09,-6.4075e-14,-33300.7,-70.8203], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/8/15""",
    longDesc = 
u"""
10/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)CCCO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 478,
    label = "NC52ONEO2-4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {4,S} {16,D}
6  O u0 p2 c0 {1,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.46745,0.0562678,-3.55389e-05,1.0713e-08,-1.1504e-12,-29553.6,12.765], Tmin=(298,'K'), Tmax=(1419,'K')),
            NASAPolynomial(coeffs=[20.5186,0.0204864,-6.73858e-06,1.01932e-09,-5.80865e-14,-35167.2,-73.8258], Tmin=(1419,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/8/15""",
    longDesc = 
u"""
10/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)CC(C)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 479,
    label = "CH3CHCOC2H5",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,S} {15,D}
5  C u1 p0 c0 {3,S} {4,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 O u0 p2 c0 {4,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.2615,0.0332446,-1.19633e-06,-1.01813e-08,3.34233e-12,-14767.4,13.1662], Tmin=(298,'K'), Tmax=(1426,'K')),
            NASAPolynomial(coeffs=[15.2591,0.0206874,-6.90185e-06,1.05491e-09,-6.05776e-14,-19770.1,-50.3404], Tmin=(1426,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/6/15""",
    longDesc = 
u"""
10/6/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C(=O)CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 480,
    label = "NC5KET23O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {4,S} {16,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  O u1 p2 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.5646,0.0518213,-2.70331e-05,5.19426e-09,1.695e-14,-27344.9,15.8031], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[19.5246,0.0208859,-7.23257e-06,1.13238e-09,-6.60902e-14,-33492.6,-72.1145], Tmin=(1383,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC([O])C(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 481,
    label = "CCYCCOOC-I2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {6,S} {13,S}
5  O u0 p2 c0 {2,S} {6,S}
6  O u0 p2 c0 {4,S} {5,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.55773,0.0710221,-6.05228e-05,2.61774e-08,-4.51391e-12,13292.2,46.0829], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[16.3389,0.0165879,-5.61829e-06,8.67331e-10,-5.01477e-14,6669.47,-64.0322], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""L 2/00""",
    longDesc = 
u"""
L 2/00
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1[CH]OOC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 482,
    label = "CVCCVCCJVO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,D} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {6,S}
3  C u0 p0 c0 {1,D} {5,S} {8,S}
4  C u0 p0 c0 {2,D} {9,S} {10,S}
5  C u1 p0 c0 {3,S} {11,D}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.218492,0.05921,-5.89241e-05,2.97412e-08,-5.85245e-12,12060.1,25.5969], Tmin=(298,'K'), Tmax=(1399,'K')),
            NASAPolynomial(coeffs=[15.3178,0.0127353,-4.35883e-06,6.76913e-10,-3.92771e-14,7605.83,-54.36], Tmin=(1399,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/ 5/ 9 THERM""",
    longDesc = 
u"""
2/ 5/ 9 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC=C[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 483,
    label = "NC5KET32",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
2  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {2,S} {17,D}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 O u0 p2 c0 {5,D}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.380342,0.0767073,-5.67635e-05,2.07068e-08,-2.96498e-12,-44157.9,32.1098], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[24.3326,0.0202927,-6.88031e-06,1.06352e-09,-6.15648e-14,-52289,-96.1326], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(=O)C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 484,
    label = "NC4KET23",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {3,S} {14,D}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 O u0 p2 c0 {4,D}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.56927,0.0520286,-3.64134e-05,1.32008e-08,-1.93577e-12,-43065.1,14.8669], Tmin=(298,'K'), Tmax=(1411,'K')),
            NASAPolynomial(coeffs=[17.6878,0.018482,-6.18385e-06,9.45792e-10,-5.42982e-14,-47859.4,-60.6682], Tmin=(1411,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 485,
    label = "C2HCHO",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {4,D} {5,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {6,S}
4 O u0 p2 c0 {1,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.20777,0.0134383,-5.15442e-06,-2.24571e-11,2.74111e-13,10211.7,5.43872], Tmin=(298,'K'), Tmax=(2012,'K')),
            NASAPolynomial(coeffs=[7.99952,0.00707825,-2.63087e-06,4.33073e-10,-2.62003e-14,8718.63,-15.7226], Tmin=(2012,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/31/13""",
    longDesc = 
u"""
1/31/13
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 486,
    label = "NC4KET21",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {2,S} {14,D}
5  O u0 p2 c0 {2,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 O u0 p2 c0 {4,D}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.388449,0.0692735,-5.58731e-05,2.25791e-08,-3.639e-12,-38687.6,35.2948], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[21.0786,0.0161788,-5.53076e-06,8.59641e-10,-4.99535e-14,-45756.4,-78.8015], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(=O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 487,
    label = "NC4KET24",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {3,S} {14,D}
5  O u0 p2 c0 {2,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 O u0 p2 c0 {4,D}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.12063,0.0501344,-3.10195e-05,9.36512e-09,-1.07549e-12,-40864.4,19.6072], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[17.4146,0.0192744,-6.57971e-06,1.02024e-09,-5.91418e-14,-46066.3,-58.0321], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 488,
    label = "HO2CH2CHO",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {7,D} {8,S}
3 O u0 p2 c0 {1,S} {4,S}
4 O u0 p2 c0 {3,S} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 O u0 p2 c0 {2,D}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.32769,0.0521619,-4.97328e-05,2.31272e-08,-4.20788e-12,-29060.9,36.186], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[15.1555,0.0075724,-2.72693e-06,4.38217e-10,-2.60434e-14,-34142,-50.1255], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 489,
    label = "IC4KETIT",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,D} {14,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.14244,0.0633841,-4.73085e-05,1.77145e-08,-2.67265e-12,-40936.7,23.4845], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[20.937,0.0171091,-6.01892e-06,9.52354e-10,-5.59926e-14,-47782,-82.7718], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/19/ 0 THERM""",
    longDesc = 
u"""
7/19/ 0 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(C=O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 490,
    label = "C3H52-1,3OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
3  C u1 p0 c0 {1,S} {2,S} {12,S}
4  O u0 p2 c0 {1,S} {6,S}
5  O u0 p2 c0 {2,S} {7,S}
6  O u0 p2 c0 {4,S} {13,S}
7  O u0 p2 c0 {5,S} {14,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.12254,0.0519554,-3.83734e-05,1.45852e-08,-2.29821e-12,-12275.9,14.8367], Tmin=(298,'K'), Tmax=(1379,'K')),
            NASAPolynomial(coeffs=[20.2818,0.0148155,-5.25503e-06,8.35963e-10,-4.93309e-14,-18008.5,-72.2688], Tmin=(1379,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/26/3 THRM""",
    longDesc = 
u"""
8/26/3 THRM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOC[CH]COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 491,
    label = "C8H132-6,PAO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {7,S} {14,S} {15,S} {16,S}
4  C u0 p0 c0 {8,S} {13,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {7,D} {20,S}
6  C u0 p0 c0 {2,S} {8,D} {21,S}
7  C u0 p0 c0 {3,S} {5,D} {19,S}
8  C u0 p0 c0 {4,S} {6,D} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 O u1 p2 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[6.30808,0.0580374,-1.75638e-05,-2.70888e-09,1.61702e-12,7078.39,6.29345], Tmin=(298,'K'), Tmax=(2014,'K')),
            NASAPolynomial(coeffs=[19.9392,0.0379082,-1.37893e-05,2.23771e-09,-1.34044e-13,1429.81,-70.7923], Tmin=(2014,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CCCC=CC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 492,
    label = "NC4KET14",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,D} {14,S}
5  O u0 p2 c0 {3,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.92379,0.0507578,-2.88361e-05,6.6465e-09,-2.83907e-13,-37294.6,19.6328], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[18.9232,0.018227,-6.27434e-06,9.78729e-10,-5.69844e-14,-43250.9,-67.8717], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=CCCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 493,
    label = "NC4KET12",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,D} {14,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.724232,0.0726649,-6.04779e-05,2.54349e-08,-4.30153e-12,-37293.7,35.6277], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[21.7577,0.0164473,-5.79962e-06,9.1915e-10,-5.41037e-14,-44711.5,-83.7725], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C=O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 494,
    label = "NC4KET13",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,D} {14,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.31776,0.0528482,-3.43212e-05,1.04563e-08,-1.12797e-12,-39486.8,15.8443], Tmin=(298,'K'), Tmax=(1411,'K')),
            NASAPolynomial(coeffs=[19.3085,0.0173455,-5.85047e-06,9.00298e-10,-5.19275e-14,-45102.4,-70.487], Tmin=(1411,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CC=O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 495,
    label = "CH3OCH2O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  O u0 p2 c0 {1,S} {2,S}
4  O u0 p2 c0 {1,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 O u1 p2 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.39931,0.030946,-1.92548e-05,5.76034e-09,-6.16082e-13,-20443.3,13.943], Tmin=(298,'K'), Tmax=(1441,'K')),
            NASAPolynomial(coeffs=[11.9179,0.0119413,-3.93526e-06,5.95756e-10,-3.39598e-14,-23423.2,-32.0097], Tmin=(1441,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/12/14 THERM""",
    longDesc = 
u"""
2/12/14 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
COCO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 496,
    label = "DC5H10OOH-B",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {1,S} {3,S} {4,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.74458,0.0578236,-3.12391e-05,7.31571e-09,-4.5089e-13,-8818.62,22.0757], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[19.2228,0.0250236,-8.57782e-06,1.33372e-09,-7.74605e-14,-15128.3,-68.4715], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C](C)CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 497,
    label = "C2H5COCH2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {10,D}
4  C u1 p0 c0 {3,S} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 O u0 p2 c0 {3,D}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.96643,0.0410271,-2.56194e-05,7.86244e-09,-9.26826e-13,-8801.49,18.4804], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[13.5979,0.0157188,-5.35201e-06,8.28428e-10,-4.79646e-14,-13011.2,-44.6216], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=O)CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 498,
    label = "CH3CHCOCH3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {12,S}
4  C u0 p0 c0 {2,S} {3,D} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 O u1 p2 c0 {4,S}
12 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.31082,0.0283712,-5.64604e-06,-4.62626e-09,1.82954e-12,-11460.2,16.2217], Tmin=(298,'K'), Tmax=(1438,'K')),
            NASAPolynomial(coeffs=[11.9651,0.0163142,-5.3988e-06,8.20455e-10,-4.69198e-14,-15199.1,-33.0143], Tmin=(1438,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 499,
    label = "NEOC5H11O2H",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {19,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.355226,0.0691131,-3.61243e-05,6.97537e-09,-9.54856e-14,-32338.9,26.1109], Tmin=(298,'K'), Tmax=(2019,'K')),
            NASAPolynomial(coeffs=[20.1975,0.029738,-1.09931e-05,1.80338e-09,-1.08852e-13,-39516.6,-82.3754], Tmin=(2019,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 500,
    label = "OCH2CHO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,D} {7,S}
3 O u1 p2 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 O u0 p2 c0 {2,D}
7 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[6.12785,0.00667998,1.00125e-05,-1.06888e-08,2.76653e-12,-12821.7,-0.809763], Tmin=(298,'K'), Tmax=(1468,'K')),
            NASAPolynomial(coeffs=[9.83473,0.00777896,-2.78175e-06,4.44912e-10,-2.6353e-14,-15074,-24.1645], Tmin=(1468,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 501,
    label = "C5H9O1-2O-5",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
6  O u0 p2 c0 {1,S} {4,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 O u1 p2 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.646244,0.0653342,-4.13295e-05,1.15535e-08,-9.90013e-13,-10103.7,33.5842], Tmin=(298,'K'), Tmax=(1421,'K')),
            NASAPolynomial(coeffs=[20.122,0.0203729,-6.96642e-06,1.08202e-09,-6.28152e-14,-17488.6,-78.9087], Tmin=(1421,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]CCCC1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 502,
    label = "CH3NO",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 N u0 p1 c0 {1,S} {6,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 O u0 p2 c0 {2,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.18535,-0.00634086,4.57171e-05,-5.30422e-08,1.99502e-11,6937.72,2.18493], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.04712,0.00921544,-3.29035e-06,5.2894e-10,-3.1569e-14,6237.18,-0.774396], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""NitrosomethyT12/09""",
    longDesc = 
u"""
NitrosomethyT12/09.
CN=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 503,
    label = "CH2CH2OCH2CH2CHO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
6  O u0 p2 c0 {4,S} {5,S}
7  O u1 p2 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.618509,0.0643128,-4.6037e-05,1.65229e-08,-2.31383e-12,-20210.3,37.5876], Tmin=(298,'K'), Tmax=(1447,'K')),
            NASAPolynomial(coeffs=[17.5698,0.020546,-6.77594e-06,1.02645e-09,-5.85405e-14,-26236.7,-59.3673], Tmin=(1447,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""4""",
    longDesc = 
u"""
4
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]C1CCOCC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 504,
    label = "DC5H10OOH-CO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {10,S}
3  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {3,S} {8,S}
7  O u0 p2 c0 {2,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.83496,0.0860948,-7.22651e-05,3.17312e-08,-5.6201e-12,-26135.5,25.6008], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[25.2421,0.0241694,-8.19141e-06,1.26382e-09,-7.30015e-14,-33530.5,-97.5222], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)C(COO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 505,
    label = "CH2O2H",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,S} {4,S} {5,S}
2 O u0 p2 c0 {1,S} {3,S}
3 O u0 p2 c0 {2,S} {6,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.88976,0.0209466,-1.75191e-05,7.2782e-09,-1.18912e-12,6123.91,12.3802], Tmin=(298,'K'), Tmax=(1410,'K')),
            NASAPolynomial(coeffs=[9.24698,0.00460846,-1.53501e-06,2.34435e-10,-1.34573e-14,4115.3,-21.1503], Tmin=(1410,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 506,
    label = "HOCO",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,S} {3,D}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p2 c0 {1,D}
4 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.92208,0.00762454,3.29884e-06,-1.07135e-08,5.11587e-12,-23028.2,11.2926], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.39206,0.00411221,-1.48195e-06,2.39875e-10,-1.43903e-14,-23860.7,-2.23529], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T05/06""",
    longDesc = 
u"""
T05/06.
O=[C]O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 507,
    label = "CC5H10OOH-AO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
3  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {3,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.83496,0.0860948,-7.22651e-05,3.17312e-08,-5.6201e-12,-26135.5,26.2903], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[25.2421,0.0241694,-8.19141e-06,1.26382e-09,-7.30015e-14,-33530.5,-96.8327], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CO[O])C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 508,
    label = "DC5H10OOH-AO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {4,S} {8,S}
7  O u0 p2 c0 {3,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.47959,0.0837669,-6.63706e-05,2.7648e-08,-4.71461e-12,-23947.4,29.9188], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[24.9285,0.0249502,-8.57455e-06,1.33548e-09,-7.76538e-14,-31716.8,-94.6333], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CCOO)CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 509,
    label = "C5H10OH-5O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
6  O u0 p2 c0 {4,S} {18,S}
7  O u0 p2 c0 {5,S} {19,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 O u1 p2 c0 {6,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.13789,0.0741287,-5.25165e-05,1.92381e-08,-2.87963e-12,-31676.8,31.1293], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[22.145,0.0249622,-8.52929e-06,1.32337e-09,-7.67473e-14,-38946.8,-81.6311], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 8/15""",
    longDesc = 
u"""
10/ 8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OCCCCCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 510,
    label = "HOCH2O2H",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p2 c0 {1,S} {7,S}
4 O u0 p2 c0 {2,S} {8,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.53519,0.0373267,-3.153e-05,1.30353e-08,-2.11473e-12,-38660.9,27.1776], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[12.4532,0.00718221,-2.4703e-06,3.85612e-10,-2.24774e-14,-42486.3,-35.8745], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 511,
    label = "AC5H10OOH-AO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {3,S} {8,S}
7  O u0 p2 c0 {4,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.47959,0.0837669,-6.63706e-05,2.7648e-08,-4.71461e-12,-23947.4,29.9188], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[24.9285,0.0249502,-8.57455e-06,1.33548e-09,-7.76538e-14,-31716.8,-94.6333], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(CO[O])COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 512,
    label = "CCYCCO-T1",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 O u0 p2 c0 {1,S} {3,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.39311,0.0381789,-2.62316e-05,8.74877e-09,-1.11297e-12,11253.5,34.1877], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[10.3395,0.011518,-3.87497e-06,5.95744e-10,-3.4354e-14,7176.59,-28.9688], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""L 2/00""",
    longDesc = 
u"""
L 2/00
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C]1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 513,
    label = "AC5H10OOH-CO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {3,S} {8,S}
7  O u0 p2 c0 {2,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.83496,0.0860948,-7.22651e-05,3.17312e-08,-5.6201e-12,-26135.5,26.2903], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[25.2421,0.0241694,-8.19141e-06,1.26382e-09,-7.30015e-14,-33530.5,-96.8327], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(COO)C(C)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 514,
    label = "NEOC5H12",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.860482,0.0463444,2.26558e-06,-1.9174e-08,6.07292e-12,-22402.1,17.8236], Tmin=(298,'K'), Tmax=(1459,'K')),
            NASAPolynomial(coeffs=[21.341,0.02325,-8.36388e-06,1.34307e-09,-7.97753e-14,-31805.6,-100.65], Tmin=(1459,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 515,
    label = "C5H9OH1-4",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {5,D} {13,S}
5  C u0 p0 c0 {4,D} {14,S} {15,S}
6  O u0 p2 c0 {3,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.542248,0.0616234,-4.1047e-05,1.39788e-08,-1.92326e-12,-22956,33.6582], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[16.679,0.0225497,-7.64069e-06,1.17888e-09,-6.81025e-14,-29031.7,-59.2057], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 8/15""",
    longDesc = 
u"""
10/ 8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 516,
    label = "TQC4H8OI",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {16,S}
7  O u1 p2 c0 {4,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.0745748,0.07465,-6.42255e-05,2.80909e-08,-4.87692e-12,-24718.3,29.4512], Tmin=(298,'K'), Tmax=(1411,'K')),
            NASAPolynomial(coeffs=[21.3201,0.018049,-6.06124e-06,9.29741e-10,-5.34977e-14,-31296.7,-82.0047], Tmin=(1411,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(C[O])OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 517,
    label = "CC5H10OOH-BO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {1,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.92758,0.0837498,-7.09477e-05,3.16521e-08,-5.69638e-12,-30570.8,17.1231], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[25.1614,0.0240667,-8.11774e-06,1.2483e-09,-7.19347e-14,-37527.4,-99.5534], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(OO)C(C)(C)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 518,
    label = "DC5H10OOH-BO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {3,S} {8,S}
7  O u0 p2 c0 {1,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.66225,0.0810932,-6.46324e-05,2.73539e-08,-4.75157e-12,-28398.3,20.3147], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[24.8033,0.0249329,-8.54e-06,1.32701e-09,-7.70334e-14,-35696.4,-97.1095], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(CCOO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 519,
    label = "CH2COHCHO",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {6,S} {7,S}
3 C u0 p0 c0 {1,S} {5,D} {8,S}
4 O u0 p2 c0 {1,S} {9,S}
5 O u0 p2 c0 {3,D}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.057664,0.051659,-5.49774e-05,2.83487e-08,-5.57464e-12,-33959.5,22.9757], Tmin=(298,'K'), Tmax=(1406,'K')),
            NASAPolynomial(coeffs=[14.02,0.00795093,-2.63245e-06,3.99896e-10,-2.28538e-14,-37680.8,-48.0458], Tmin=(1406,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/24/15""",
    longDesc = 
u"""
9/24/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(O)C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 520,
    label = "NC53ONE-1",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {11,D}
4  C u0 p0 c0 {3,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 O u0 p2 c0 {3,D}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.06846,0.049705,-2.97718e-05,8.51066e-09,-8.95411e-13,-20258,23.473], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[15.4699,0.0194419,-6.68487e-06,1.04158e-09,-6.05842e-14,-25603.8,-55.0816], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 8/15 THERM""",
    longDesc = 
u"""
10/ 8/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(=O)CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 521,
    label = "IC4H8OH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {1,S} {12,S} {13,S}
5  O u0 p2 c0 {2,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.29613,0.034765,-1.02506e-05,-2.04642e-09,1.18879e-12,-14562.7,15.8606], Tmin=(298,'K'), Tmax=(1376,'K')),
            NASAPolynomial(coeffs=[12.5606,0.0210637,-7.1502e-06,1.10439e-09,-6.38429e-14,-18620.3,-36.7889], Tmin=(1376,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/14/95 THERM""",
    longDesc = 
u"""
2/14/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 522,
    label = "BC5H10OOH-AO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
6  O u0 p2 c0 {1,S} {8,S}
7  O u0 p2 c0 {3,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.66225,0.0810932,-6.46324e-05,2.73539e-08,-4.75157e-12,-28398.3,21.0042], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[24.8033,0.0249329,-8.54e-06,1.32701e-09,-7.70334e-14,-35696.4,-96.4201], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)(CO[O])OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 523,
    label = "IIC4H7Q2-I",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {14,S} {15,S}
5  O u0 p2 c0 {2,S} {7,S}
6  O u0 p2 c0 {3,S} {8,S}
7  O u0 p2 c0 {5,S} {16,S}
8  O u0 p2 c0 {6,S} {17,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.93056,0.0605819,-4.23666e-05,1.49122e-08,-2.10979e-12,-16841.5,13.6228], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[23.05,0.0192149,-6.66623e-06,1.04496e-09,-6.10371e-14,-23208.7,-83.995], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/15/96 THERM""",
    longDesc = 
u"""
7/15/96 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(COO)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 524,
    label = "C4H71-1O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {11,S}
4  C u0 p0 c0 {3,D} {5,S} {12,S}
5  O u0 p2 c0 {4,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 O u1 p2 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.02223,0.0485578,-3.30294e-05,1.13407e-08,-1.57078e-12,4041.44,20.6107], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[16.3739,0.0162685,-5.62192e-06,8.78983e-10,-5.12511e-14,-1063.18,-56.9006], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/29/15""",
    longDesc = 
u"""
9/29/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC=CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 525,
    label = "CH3CO3",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,S} {7,D}
3 O u0 p2 c0 {2,S} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 O u0 p2 c0 {2,D}
8 O u1 p2 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.60373,0.027008,-2.08293e-05,8.50541e-09,-1.43846e-12,-23420.5,11.2015], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[11.2522,0.00833653,-2.89015e-06,4.52782e-10,-2.64354e-14,-26023.9,-29.637], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""4/ 3/ 0 THERM""",
    longDesc = 
u"""
4/ 3/ 0 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 526,
    label = "CH3CO2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 O u1 p2 c0 {2,S}
7 O u0 p2 c0 {2,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.37441,0.0249116,-1.74309e-05,6.248e-09,-9.09517e-13,-27233,18.1405], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[8.5406,0.00832951,-2.84722e-06,4.41927e-10,-2.56373e-14,-29729.1,-20.3884], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/14/95 THERM""",
    longDesc = 
u"""
2/14/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC([O])=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 527,
    label = "IC4H8OOH-IO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {2,S} {7,S}
6  O u0 p2 c0 {3,S} {16,S}
7  O u0 p2 c0 {5,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u1 p2 c0 {6,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.23355,0.0563089,-3.15673e-05,7.79537e-09,-6.21665e-13,-22278.3,15.2623], Tmin=(298,'K'), Tmax=(1367,'K')),
            NASAPolynomial(coeffs=[22.4665,0.0209351,-7.44324e-06,1.18589e-09,-7.00547e-14,-29449.5,-85.4241], Tmin=(1367,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CO[O])COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 528,
    label = "NEOC5H11O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {2,S} {18,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.08145,0.0633535,-3.11939e-05,4.88884e-09,2.48584e-13,-15645.5,23.2475], Tmin=(298,'K'), Tmax=(1682,'K')),
            NASAPolynomial(coeffs=[19.3713,0.0282613,-1.0467e-05,1.71917e-09,-1.03854e-13,-22354.4,-77.1662], Tmin=(1682,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(C)CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 529,
    label = "C5H9A-COOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {4,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {15,S} {16,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.23219,0.0831069,-6.96695e-05,2.98091e-08,-5.11291e-12,-16666.6,39.3337], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[22.1763,0.0202998,-6.93528e-06,1.07654e-09,-6.24713e-14,-24506.1,-89.6059], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 530,
    label = "C5H9A-DOOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {4,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {4,D} {15,S} {16,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.584051,0.0661435,-4.64794e-05,1.71321e-08,-2.62129e-12,-16656.5,29.0021], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[19.1192,0.0231084,-7.94263e-06,1.23714e-09,-7.19377e-14,-23156.3,-70.6977], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 531,
    label = "C3H5O1-2OOH-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  O u0 p2 c0 {1,S} {2,S}
5  O u0 p2 c0 {3,S} {6,S}
6  O u0 p2 c0 {5,S} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.25001,0.0665787,-6.1886e-05,2.84639e-08,-5.08512e-12,-22237.1,43.0381], Tmin=(298,'K'), Tmax=(1432,'K')),
            NASAPolynomial(coeffs=[15.7042,0.0130256,-4.23544e-06,6.35556e-10,-3.6011e-14,-27726.9,-55.1895], Tmin=(1432,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/13 THER""",
    longDesc = 
u"""
10/13 THER
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOCC1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 532,
    label = "C4H71-2O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {5,S}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  O u0 p2 c0 {3,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 O u1 p2 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.4246,0.0520548,-3.89441e-05,1.51726e-08,-2.42612e-12,1688.54,22.3231], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[16.3566,0.0161768,-5.5661e-06,8.67701e-10,-5.04886e-14,-3402.45,-57.5261], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/29/15""",
    longDesc = 
u"""
9/29/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(CC)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 533,
    label = "C5H9A-B,DOOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {1,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.87833,0.0817085,-6.64694e-05,2.85545e-08,-5.00711e-12,-20509.3,21.6843], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[25.3148,0.024105,-8.26152e-06,1.28429e-09,-7.45769e-14,-27825.1,-97.0436], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)(CCOO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 534,
    label = "C5H9A-B,COOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {1,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {8,S}
7  O u0 p2 c0 {2,S} {9,S}
8  O u0 p2 c0 {6,S} {19,S}
9  O u0 p2 c0 {7,S} {20,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.12195,0.0844868,-7.2887e-05,3.29013e-08,-5.96275e-12,-22680,18.5842], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[25.6998,0.0232725,-7.86284e-06,1.21048e-09,-6.98108e-14,-29677.2,-99.6712], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)(OO)C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 535,
    label = "C5H9OA-DOOH-B",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {3,S} {4,S}
7  O u0 p2 c0 {1,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.94627,0.105725,-9.40624e-05,4.23682e-08,-7.53895e-12,-40979,70.6043], Tmin=(298,'K'), Tmax=(1424,'K')),
            NASAPolynomial(coeffs=[20.9161,0.0242618,-8.18692e-06,1.25968e-09,-7.26328e-14,-50050,-85.4044], Tmin=(1424,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(CCOC1)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 536,
    label = "CH2CHOOHCOCH3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {5,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {11,D}
4  C u1 p0 c0 {1,S} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 O u0 p2 c0 {3,D}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.45962,0.04672,-3.15879e-05,1.06246e-08,-1.38857e-12,-18472.1,12.9656], Tmin=(298,'K'), Tmax=(1413,'K')),
            NASAPolynomial(coeffs=[17.8031,0.0160286,-5.38152e-06,8.25242e-10,-4.74719e-14,-23076.8,-58.737], Tmin=(1413,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(OO)C(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 537,
    label = "C2H5COCH3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {3,S} {13,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.57048,0.0351447,-1.2385e-05,-1.21281e-09,1.16164e-12,-31019.4,16.1395], Tmin=(298,'K'), Tmax=(1454,'K')),
            NASAPolynomial(coeffs=[12.8183,0.0179874,-5.94195e-06,9.01636e-10,-5.14994e-14,-35171.2,-41.1609], Tmin=(1454,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/12/15""",
    longDesc = 
u"""
8/12/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 538,
    label = "CJVCCVCCVO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,D} {3,S} {6,S}
2  C u0 p0 c0 {1,D} {4,S} {7,S}
3  C u0 p0 c0 {1,S} {5,D} {8,S}
4  C u0 p0 c0 {2,S} {9,D} {10,S}
5  C u1 p0 c0 {3,D} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  O u0 p2 c0 {4,D}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.506629,0.0604672,-5.97397e-05,2.96804e-08,-5.7624e-12,24276.6,28.2994], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[16.2361,0.0118297,-4.11454e-06,6.46027e-10,-3.77768e-14,19350,-58.3499], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/ 5/ 9 THERM""",
    longDesc = 
u"""
2/ 5/ 9 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=CC=CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 539,
    label = "CH3CO3H",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2 C u0 p0 c0 {1,S} {3,S} {8,D}
3 O u0 p2 c0 {2,S} {4,S}
4 O u0 p2 c0 {3,S} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
8 O u0 p2 c0 {2,D}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.24136,0.0337964,-2.53887e-05,9.67584e-09,-1.49266e-12,-42467.8,17.0668], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[12.506,0.0094779,-3.30402e-06,5.19631e-10,-3.04234e-14,-45985.7,-37.9196], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""6/26/95 THERM""",
    longDesc = 
u"""
6/26/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 540,
    label = "NC5CYCPER13",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {7,S} {8,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {4,S} {6,S}
8  O u0 p2 c0 {4,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.99454,0.0945988,-7.58465e-05,3.02153e-08,-4.76417e-12,-45988.1,51.3024], Tmin=(298,'K'), Tmax=(1407,'K')),
            NASAPolynomial(coeffs=[24.1561,0.0224096,-7.74874e-06,1.21234e-09,-7.07323e-14,-55534.5,-103.522], Tmin=(1407,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""11/15 THERM""",
    longDesc = 
u"""
11/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1CC(O)OO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 541,
    label = "CH2COHCO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 O u0 p2 c0 {1,S} {8,S}
4 C u1 p0 c0 {1,S} {7,D}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 O u0 p2 c0 {4,D}
8 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.24341,0.0471485,-5.442e-05,2.9597e-08,-6.01447e-12,-16555.5,16.8301], Tmin=(298,'K'), Tmax=(1411,'K')),
            NASAPolynomial(coeffs=[13.1285,0.00617949,-1.93388e-06,2.81893e-10,-1.56244e-14,-19350.8,-42.2958], Tmin=(1411,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/25/15""",
    longDesc = 
u"""
9/25/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(O)[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 542,
    label = "NEOC5H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {1,S} {15,S} {16,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
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
            NASAPolynomial(coeffs=[9.68072,0.00707718,4.53302e-05,-3.61952e-08,8.06654e-12,1051.17,-19.4358], Tmin=(298,'K'), Tmax=(1375,'K')),
            NASAPolynomial(coeffs=[11.7983,0.035958,-1.44109e-05,2.46308e-09,-1.52118e-13,-3861.01,-44.1297], Tmin=(1375,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)(C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 543,
    label = "NC5CYCPER31",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {4,S} {6,S}
8  O u0 p2 c0 {1,S} {18,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.01648,0.0890856,-7.09557e-05,2.89247e-08,-4.75384e-12,-45627.4,46.5939], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[22.0861,0.0239837,-8.24799e-06,1.28547e-09,-7.47886e-14,-54216.2,-92.0233], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""11/15 THERM""",
    longDesc = 
u"""
11/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1(O)CCOO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 544,
    label = "C5H91-2,4OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {5,S} {7,S} {11,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {3,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {8,S}
7  O u0 p2 c0 {3,S} {9,S}
8  O u0 p2 c0 {6,S} {19,S}
9  O u0 p2 c0 {7,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.91226,0.0795651,-6.1323e-05,2.44461e-08,-3.95034e-12,-19651.7,23.2188], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[25.5056,0.023699,-8.07033e-06,1.24935e-09,-7.23427e-14,-27135.3,-96.9349], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(CC(C)OO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 545,
    label = "C5H9D-A,BOOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {2,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.27985,0.0834372,-6.86098e-05,2.97416e-08,-5.24931e-12,-20365.8,24.5532], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[25.1744,0.0242094,-8.29427e-06,1.28907e-09,-7.48414e-14,-27783.7,-96.4353], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(C)(COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 546,
    label = "NC5KET31O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {10,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {2,S} {16,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 O u1 p2 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.37525,0.0349659,3.6792e-06,-1.56751e-08,4.84743e-12,-26124.8,8.64756], Tmin=(298,'K'), Tmax=(1426,'K')),
            NASAPolynomial(coeffs=[19.2416,0.0207625,-7.12773e-06,1.11093e-09,-6.46776e-14,-32571.1,-71.9638], Tmin=(1426,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(=O)CC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 547,
    label = "CH3CHOOCOCH3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {3,S} {13,D}
5  O u0 p2 c0 {1,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
14 O u1 p2 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.30172,0.0462274,-3.12564e-05,1.08628e-08,-1.51478e-12,-26373.2,11.9724], Tmin=(298,'K'), Tmax=(1411,'K')),
            NASAPolynomial(coeffs=[16.8056,0.0170791,-5.69439e-06,8.68879e-10,-4.98008e-14,-30671.9,-55.1179], Tmin=(1411,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)C(C)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 548,
    label = "IC5KETCAO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {4,S} {16,D}
6  H u0 p0 c0 {1,S}
7  O u1 p2 c0 {3,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.75737,0.0452906,-2.09521e-05,2.65372e-09,4.60703e-13,-27333.3,10.2325], Tmin=(298,'K'), Tmax=(1514,'K')),
            NASAPolynomial(coeffs=[17.2486,0.0216327,-7.22891e-06,1.10468e-09,-6.33813e-14,-32137.1,-58.6643], Tmin=(1514,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)C(C)C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 549,
    label = "PC4H8CHO-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {2,S} {5,S} {13,S}
5  C u0 p0 c0 {4,S} {14,D} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 O u0 p2 c0 {5,D}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.88447,0.0526233,-2.83374e-05,5.66304e-09,-3.63459e-15,-9892.41,25.025], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[17.3204,0.0200051,-6.88463e-06,1.07357e-09,-6.24869e-14,-16087.9,-65.1443], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[CH]C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 550,
    label = "C5H10OOH1-5O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {17,S} {18,S}
5  C u0 p0 c0 {2,S} {7,S} {15,S} {16,S}
6  O u0 p2 c0 {4,S} {8,S}
7  O u0 p2 c0 {5,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.61115,0.0758125,-5.0394e-05,1.63241e-08,-2.04839e-12,-23397,25.8821], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[25.8212,0.0246141,-8.55566e-06,1.34293e-09,-7.85178e-14,-31715.8,-99.7698], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OCCCCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 551,
    label = "AC5H10OOH-D",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {2,S} {16,S} {17,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.0893968,0.0748083,-5.79088e-05,2.38676e-08,-4.05836e-12,-5791.21,33.772], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[20.2759,0.0240402,-8.21726e-06,1.27509e-09,-7.3947e-14,-12589,-74.5433], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 552,
    label = "PC3H4OH-3",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {5,S}
2 C u1 p0 c0 {1,S} {4,S} {6,S}
3 C u0 p0 c0 {1,D} {7,S} {8,S}
4 O u0 p2 c0 {2,S} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.406368,0.0370231,-2.72203e-05,1.05783e-08,-1.73633e-12,527.161,23.4782], Tmin=(298,'K'), Tmax=(1378,'K')),
            NASAPolynomial(coeffs=[11.5298,0.0114381,-4.04415e-06,6.41949e-10,-3.78242e-14,-3454.86,-36.5385], Tmin=(1378,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/25/15""",
    longDesc = 
u"""
9/25/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[CH]O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 553,
    label = "C5H9B-A,COOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {10,S}
2  C u0 p0 c0 {5,S} {7,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {5,S} {16,S} {17,S} {18,S}
5  C u1 p0 c0 {1,S} {2,S} {4,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {2,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.06775,0.071334,-5.22649e-05,2.02919e-08,-3.25625e-12,-21450.1,14.0311], Tmin=(298,'K'), Tmax=(1675,'K')),
            NASAPolynomial(coeffs=[24.4196,0.0245183,-8.32769e-06,1.2867e-09,-7.43971e-14,-28021.1,-89.3467], Tmin=(1675,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C](COO)C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 554,
    label = "C4H8OOH2-4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {2,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.82673,0.0569912,-4.24389e-05,1.67119e-08,-2.70636e-12,-4546.3,22.2052], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[17.2354,0.0191946,-6.49947e-06,1.00211e-09,-5.7856e-14,-9711.58,-59.8986], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 555,
    label = "C5H93-1,5OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {6,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {7,S} {16,S} {17,S}
5  C u1 p0 c0 {1,S} {2,S} {18,S}
6  O u0 p2 c0 {3,S} {8,S}
7  O u0 p2 c0 {4,S} {9,S}
8  O u0 p2 c0 {6,S} {19,S}
9  O u0 p2 c0 {7,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.67512,0.0698557,-4.39194e-05,1.3129e-08,-1.43359e-12,-17058.8,23.0263], Tmin=(298,'K'), Tmax=(1399,'K')),
            NASAPolynomial(coeffs=[24.4185,0.0249502,-8.57441e-06,1.33559e-09,-7.76698e-14,-24553,-89.5288], Tmin=(1399,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOCC[CH]CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 556,
    label = "C5H9O1-3OOH-4",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {4,S}
7  O u0 p2 c0 {2,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.67709,0.0925192,-7.56246e-05,3.08656e-08,-4.96975e-12,-30517.9,51.9808], Tmin=(298,'K'), Tmax=(1421,'K')),
            NASAPolynomial(coeffs=[22.9239,0.0222245,-7.5093e-06,1.15689e-09,-6.67794e-14,-39314.9,-93.8427], Tmin=(1421,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(OO)C1CCO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 557,
    label = "IC4H6OOH-I",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u1 p0 c0 {2,S} {9,S} {10,S}
4  C u0 p0 c0 {2,D} {11,S} {12,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {13,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[6.06669,0.039495,-2.52722e-05,7.84486e-09,-8.98877e-13,2874.48,-0.378377], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[17.3601,0.0142046,-4.65348e-06,7.0221e-10,-3.99553e-14,-1075.09,-61.2822], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""L 2/00""",
    longDesc = 
u"""
L 2/00
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 558,
    label = "IC4H8OOH-TO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {2,S} {7,S}
6  O u0 p2 c0 {1,S} {16,S}
7  O u0 p2 c0 {5,S} {17,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u1 p2 c0 {6,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.36414,0.0693743,-5.70416e-05,2.4604e-08,-4.32849e-12,-25113.8,16.5767], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[23.2465,0.0188385,-6.40938e-06,9.92649e-10,-5.75276e-14,-31653.3,-88.8302], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(COO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 559,
    label = "PC4H8CHO-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {2,S} {13,S}
5  C u0 p0 c0 {2,S} {14,D} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 O u0 p2 c0 {5,D}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.19129,0.0391487,-8.88574e-06,-5.8634e-09,2.44961e-12,-7026.34,18.4683], Tmin=(298,'K'), Tmax=(1426,'K')),
            NASAPolynomial(coeffs=[15.8456,0.0208431,-7.08563e-06,1.09605e-09,-6.34446e-14,-12430,-53.2908], Tmin=(1426,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 560,
    label = "C2H3O1-2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u1 p0 c0 {1,S} {3,S} {6,S}
3 O u0 p2 c0 {1,S} {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.58349,-0.00602276,6.32427e-05,-8.18541e-08,3.30445e-11,18568.1,9.59726], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.60158,0.00917614,-3.28029e-06,5.27904e-10,-3.15362e-14,17144.6,-5.47229], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""A 1/05""",
    longDesc = 
u"""
A 1/05.
[CH]1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 561,
    label = "NC53ONEO2-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {2,S} {16,D}
6  O u0 p2 c0 {1,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 O u0 p2 c0 {5,D}
17 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.18594,0.0706411,-5.13728e-05,1.83401e-08,-2.55856e-12,-27478.1,28.8674], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[23.3147,0.01911,-6.48745e-06,1.00363e-09,-5.81318e-14,-35052.8,-89.8266], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/8/15""",
    longDesc = 
u"""
10/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(=O)C(C)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 562,
    label = "C5H10OH12",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {6,S} {14,S} {15,S}
5  C u1 p0 c0 {2,S} {4,S} {16,S}
6  O u0 p2 c0 {4,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.84754,0.0573087,-3.47271e-05,1.06516e-08,-1.30778e-12,-14799.6,24.015], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[17.1781,0.0242601,-8.22987e-06,1.2705e-09,-7.34156e-14,-20429,-59.3472], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 7/15 THERM""",
    longDesc = 
u"""
10/ 7/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[CH]CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 563,
    label = "C5H10OH13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
2  C u0 p0 c0 {4,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {1,S} {2,S} {16,S}
6  O u0 p2 c0 {3,S} {17,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.541168,0.0566862,-2.94678e-05,5.84905e-09,-3.31866e-14,-15094.4,31.8827], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[16.8143,0.0244876,-8.29373e-06,1.27921e-09,-7.38805e-14,-21256.8,-57.4573], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 7/15 THERM""",
    longDesc = 
u"""
10/ 7/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]CCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 564,
    label = "C5H9B-A,DOOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
4  C u0 p0 c0 {5,S} {16,S} {17,S} {18,S}
5  C u1 p0 c0 {1,S} {3,S} {4,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {3,S} {9,S}
8  O u0 p2 c0 {6,S} {19,S}
9  O u0 p2 c0 {7,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.2835,0.0683511,-4.42293e-05,1.46615e-08,-1.99883e-12,-19121,20.2163], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[23.7033,0.0257554,-8.89409e-06,1.3897e-09,-8.09867e-14,-26201.5,-85.162], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C](CCOO)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 565,
    label = "C5H9OA-COOH-B",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
6  O u0 p2 c0 {2,S} {3,S}
7  O u0 p2 c0 {1,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-6.82292,0.103424,-9.5179e-05,4.41171e-08,-8.01855e-12,-33450.8,59.2929], Tmin=(298,'K'), Tmax=(1427,'K')),
            NASAPolynomial(coeffs=[21.8581,0.0226583,-7.55183e-06,1.15215e-09,-6.60356e-14,-41894,-89.6205], Tmin=(1427,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1OCC1(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 566,
    label = "C5H9OB-COOH-A",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
6  O u0 p2 c0 {1,S} {2,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.03627,0.0844316,-6.63448e-05,2.5862e-08,-3.95367e-12,-32045,39.1517], Tmin=(298,'K'), Tmax=(1427,'K')),
            NASAPolynomial(coeffs=[23.366,0.0212968,-7.16312e-06,1.10028e-09,-6.33839e-14,-40271.1,-95.565], Tmin=(1427,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1OC1(C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 567,
    label = "CH3COCHOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,D} {8,S}
3  C u0 p0 c0 {2,D} {4,S} {9,S}
4  O u0 p2 c0 {3,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  O u1 p2 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.11184,0.0332521,-1.96834e-05,5.40737e-09,-5.32451e-13,-25854.5,17.7643], Tmin=(298,'K'), Tmax=(1378,'K')),
            NASAPolynomial(coeffs=[12.3885,0.01231,-4.32264e-06,6.83044e-10,-4.01193e-14,-29760.9,-38.5683], Tmin=(1378,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/25/15""",
    longDesc = 
u"""
9/25/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC([O])=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 568,
    label = "A-CC5H10O",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
3  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
6  O u0 p2 c0 {2,S} {3,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.16061,0.087014,-7.59661e-05,3.37198e-08,-5.87393e-12,-19202.8,58.9165], Tmin=(298,'K'), Tmax=(1429,'K')),
            NASAPolynomial(coeffs=[16.2176,0.0218776,-7.00341e-06,1.03851e-09,-5.83155e-14,-26045.3,-62.5263], Tmin=(1429,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1COC1C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 569,
    label = "C4H72-2O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,D} {5,S}
4  C u0 p0 c0 {2,S} {3,D} {12,S}
5  O u0 p2 c0 {3,S} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {4,S}
13 O u1 p2 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.49537,0.0475681,-3.26453e-05,1.1407e-08,-1.61495e-12,176.384,16.7027], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[16.3301,0.0160699,-5.50238e-06,8.55099e-10,-4.96517e-14,-4700.4,-57.8736], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/29/15""",
    longDesc = 
u"""
9/29/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C(C)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 570,
    label = "NEO-C5H10O",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
6  O u0 p2 c0 {2,S} {3,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.80125,0.0719986,-4.12061e-05,8.52948e-09,2.15532e-14,-18591.8,51.3449], Tmin=(298,'K'), Tmax=(1465,'K')),
            NASAPolynomial(coeffs=[18.8696,0.0218446,-7.48657e-06,1.16574e-09,-6.78305e-14,-27653.8,-83.3848], Tmin=(1465,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(C)COC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 571,
    label = "C5H9B-C,DOOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
3  C u0 p0 c0 {5,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {5,S} {16,S} {17,S} {18,S}
5  C u1 p0 c0 {1,S} {3,S} {4,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {2,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.06775,0.071334,-5.22649e-05,2.02919e-08,-3.25625e-12,-21450.1,12.6521], Tmin=(298,'K'), Tmax=(1675,'K')),
            NASAPolynomial(coeffs=[24.4196,0.0245183,-8.32769e-06,1.2867e-09,-7.43971e-14,-28021.1,-90.7256], Tmin=(1675,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C](C)C(COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 572,
    label = "C5H9A-C,DOOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {11,S}
3  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {1,S} {17,S} {18,S}
6  O u0 p2 c0 {2,S} {9,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.46543,0.0883867,-7.6201e-05,3.41129e-08,-6.11972e-12,-18105,29.7797], Tmin=(298,'K'), Tmax=(1406,'K')),
            NASAPolynomial(coeffs=[25.5933,0.0234677,-7.95419e-06,1.2273e-09,-7.0896e-14,-25609.6,-96.7335], Tmin=(1406,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)C(COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 573,
    label = "C5H10OH-4O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {18,S}
7  O u0 p2 c0 {4,S} {19,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 O u1 p2 c0 {6,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.54218,0.0757844,-5.63028e-05,2.16375e-08,-3.37143e-12,-33865.1,27.3478], Tmin=(298,'K'), Tmax=(1406,'K')),
            NASAPolynomial(coeffs=[22.8162,0.024075,-8.15387e-06,1.2576e-09,-7.26293e-14,-40989.4,-86.0866], Tmin=(1406,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 8/15""",
    longDesc = 
u"""
10/ 8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CCCO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 574,
    label = "C4H72-1OH",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {12,S}
4  C u0 p0 c0 {2,S} {3,D} {11,S}
5  O u0 p2 c0 {1,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.94638,0.0345853,-1.10025e-05,-1.33281e-09,9.54964e-13,-21224.2,15.5037], Tmin=(298,'K'), Tmax=(1367,'K')),
            NASAPolynomial(coeffs=[13.2893,0.0193844,-6.80733e-06,1.0756e-09,-6.31699e-14,-25862.8,-43.4782], Tmin=(1367,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 575,
    label = "C7H111-5,1P",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {5,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {5,D} {15,S}
5  C u0 p0 c0 {3,S} {4,D} {16,S}
6  C u0 p0 c0 {2,S} {7,D} {17,S}
7  C u1 p0 c0 {6,D} {18,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.816832,0.0656745,-4.14491e-05,1.32235e-08,-1.6929e-12,32764.6,27.5443], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[19.1916,0.0255013,-8.65882e-06,1.33797e-09,-7.73785e-14,26094.6,-72.1388], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=CCCC=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 576,
    label = "C5H9C-AOOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {14,S}
5  C u0 p0 c0 {4,D} {15,S} {16,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.0133188,0.0692568,-5.06647e-05,1.94154e-08,-3.06921e-12,-16025.4,31.76], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[19.5116,0.0227723,-7.8263e-06,1.21895e-09,-7.08784e-14,-22737.7,-72.8225], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 577,
    label = "C5H9B-A,AOOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
3  C u0 p0 c0 {5,S} {7,S} {17,S} {18,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {1,S} {2,S} {3,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {3,S} {9,S}
8  O u0 p2 c0 {6,S} {19,S}
9  O u0 p2 c0 {7,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.2835,0.0683511,-4.42293e-05,1.46615e-08,-1.99883e-12,-19121,19.5268], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[23.7033,0.0257554,-8.89409e-06,1.3897e-09,-8.09867e-14,-26201.5,-85.8515], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[C](COO)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 578,
    label = "C5H10OH11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {3,S} {6,S} {16,S}
6  O u0 p2 c0 {5,S} {17,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.42605,0.0631258,-4.15123e-05,1.38896e-08,-1.86732e-12,-16759.7,29.2144], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[18.0953,0.0233361,-7.88464e-06,1.21428e-09,-7.00616e-14,-23020.7,-66.1676], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 7/15 THERM""",
    longDesc = 
u"""
10/ 7/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC[CH]O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 579,
    label = "O-C6H4O2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {6,S} {7,D}
2  C u0 p0 c0 {1,S} {3,D} {9,S}
3  C u0 p0 c0 {2,D} {4,S} {10,S}
4  C u0 p0 c0 {3,S} {5,D} {11,S}
5  C u0 p0 c0 {4,D} {6,S} {12,S}
6  C u0 p0 c0 {1,S} {5,S} {8,D}
7  O u0 p2 c0 {1,D}
8  O u0 p2 c0 {6,D}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.3618,0.0686058,-6.3913e-05,3.06903e-08,-5.97358e-12,-12670.4,35.3724], Tmin=(270,'K'), Tmax=(1370,'K')),
            NASAPolynomial(coeffs=[12.3614,0.0240491,-1.16529e-05,2.71333e-09,-2.47593e-13,-16708,-40.0311], Tmin=(1370,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (270,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""AK0405""",
    longDesc = 
u"""
AK0405.
O=C1C=CC=CC1=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 580,
    label = "P-C6H4O2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {6,S} {7,D}
2  C u0 p0 c0 {1,S} {3,D} {9,S}
3  C u0 p0 c0 {2,D} {4,S} {10,S}
4  C u0 p0 c0 {3,S} {5,S} {8,D}
5  C u0 p0 c0 {4,S} {6,D} {11,S}
6  C u0 p0 c0 {1,S} {5,D} {12,S}
7  O u0 p2 c0 {1,D}
8  O u0 p2 c0 {4,D}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.4317,0.0687938,-6.41383e-05,3.08127e-08,-5.99832e-12,-16569.7,34.8309], Tmin=(270,'K'), Tmax=(1370,'K')),
            NASAPolynomial(coeffs=[12.3424,0.0240613,-1.16565e-05,2.71394e-09,-2.47643e-13,-20618.5,-40.8244], Tmin=(1370,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (270,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""AK0405""",
    longDesc = 
u"""
AK0405.
O=C1C=CC(=O)C=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 581,
    label = "C3H6OOH1-3O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  O u0 p2 c0 {2,S} {6,S}
5  O u0 p2 c0 {3,S} {13,S}
6  O u0 p2 c0 {4,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 O u1 p2 c0 {5,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.56933,0.0468523,-3.58918e-05,1.43315e-08,-2.29776e-12,-18606.6,7.18655], Tmin=(298,'K'), Tmax=(1416,'K')),
            NASAPolynomial(coeffs=[18.1662,0.0147645,-4.74843e-06,7.06972e-10,-3.98306e-14,-22625.6,-59.3719], Tmin=(1416,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OCCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 582,
    label = "B2E3M1OJ",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {7,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
3  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {3,S} {4,D} {15,S}
6  O u1 p2 c0 {3,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.64763,0.034361,-5.6065e-06,-5.19111e-09,1.73155e-12,-6.65738,3.13825], Tmin=(298,'K'), Tmax=(2003,'K')),
            NASAPolynomial(coeffs=[13.5666,0.0256548,-9.37226e-06,1.52518e-09,-9.15379e-14,-3587.67,-42.7273], Tmin=(2003,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)=CC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 583,
    label = "C5H9D-A,AOOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
5  C u1 p0 c0 {2,S} {17,S} {18,S}
6  O u0 p2 c0 {3,S} {8,S}
7  O u0 p2 c0 {4,S} {9,S}
8  O u0 p2 c0 {6,S} {19,S}
9  O u0 p2 c0 {7,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.03378,0.0863931,-7.07032e-05,3.02282e-08,-5.2511e-12,-15905.8,33.0645], Tmin=(298,'K'), Tmax=(1399,'K')),
            NASAPolynomial(coeffs=[25.3362,0.0242043,-8.32306e-06,1.29685e-09,-7.54307e-14,-23822.8,-95.5642], Tmin=(1399,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(COO)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 584,
    label = "CVCYCCOC",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {10,S} {11,S}
5  O u0 p2 c0 {1,S} {2,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.45567,0.0488377,-3.46186e-05,1.26591e-08,-1.88738e-12,-643.461,35.8434], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[11.4364,0.0163251,-5.57146e-06,8.63796e-10,-5.00702e-14,-5446.53,-38.7179], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""L 2/00""",
    longDesc = 
u"""
L 2/00
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C1COC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 585,
    label = "C4H4O",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {4,D} {6,S}
2 C u0 p0 c0 {1,S} {3,D} {7,S}
3 C u0 p0 c0 {2,D} {5,S} {8,S}
4 C u0 p0 c0 {1,D} {5,S} {9,S}
5 O u0 p2 c0 {3,S} {4,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.847469,0.0131774,5.99736e-05,-9.71563e-08,4.22734e-11,-5367.85,21.4945], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[9.38935,0.0140291,-5.07755e-06,8.24137e-10,-4.9532e-14,-8682.42,-27.9163], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T03/97""",
    longDesc = 
u"""
T03/97.
C1C=COC=1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 586,
    label = "AC3H5CHO",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,D} {7,S}
3  C u0 p0 c0 {1,S} {8,D} {9,S}
4  C u0 p0 c0 {2,D} {10,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  O u0 p2 c0 {3,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.135628,0.0432103,-2.86019e-05,9.48373e-09,-1.26554e-12,-10810.1,29.3765], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[12.7609,0.0148417,-5.16797e-06,8.12053e-10,-4.75112e-14,-15480.2,-40.5363], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 587,
    label = "C2H4OCHO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,D} {10,S}
4  O u1 p2 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  O u0 p2 c0 {3,D}
10 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.11985,0.0319514,-1.67509e-05,3.25915e-09,-4.96819e-14,-17623.9,14.2845], Tmin=(298,'K'), Tmax=(1684,'K')),
            NASAPolynomial(coeffs=[12.2682,0.0137337,-5.06733e-06,8.30266e-10,-5.00728e-14,-20912.6,-35.6985], Tmin=(1684,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC([O])C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 588,
    label = "OCH2O2H",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 O u0 p2 c0 {1,S} {3,S}
3 O u0 p2 c0 {2,S} {7,S}
4 O u1 p2 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.446349,0.036305,-3.26131e-05,1.37051e-08,-2.20873e-12,-14197.3,27.296], Tmin=(298,'K'), Tmax=(1418,'K')),
            NASAPolynomial(coeffs=[12.9622,0.00421949,-1.54275e-06,2.50413e-10,-1.49856e-14,-18132.6,-38.7016], Tmin=(1418,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/21/14 THERM""",
    longDesc = 
u"""
7/21/14 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 589,
    label = "C3H2C",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {3,S} {5,S}
3 C u2 p0 c0 {1,S} {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.12959,0.0172874,-1.13668e-05,3.45693e-09,-3.6616e-13,58419.1,17.3314], Tmin=(200,'K'), Tmax=(1500,'K')),
            NASAPolynomial(coeffs=[6.56327,0.00523633,-1.75448e-06,2.68661e-10,-1.54285e-14,56514.6,-12.0006], Tmin=(1500,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""0""",
    longDesc = 
u"""
0.
[C]1C=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 590,
    label = "CH3ONO",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 O u0 p2 c0 {1,S} {3,S}
3 N u0 p1 c0 {2,S} {7,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 O u0 p2 c0 {3,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[6.15261,-0.00291937,4.14527e-05,-4.93955e-08,1.85608e-11,-9852.6,0.804057], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.93605,0.00997319,-3.60643e-06,5.83462e-10,-3.50059e-14,-10838.2,-6.98145], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""A 5/05""",
    longDesc = 
u"""
A 5/05.
CON=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 591,
    label = "PC4H8CHO-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {3,S} {13,S}
5  C u0 p0 c0 {2,S} {14,D} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 O u0 p2 c0 {5,D}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.19129,0.0391487,-8.88574e-06,-5.8634e-09,2.44961e-12,-7026.34,18.4683], Tmin=(298,'K'), Tmax=(1426,'K')),
            NASAPolynomial(coeffs=[15.8456,0.0208431,-7.08563e-06,1.09605e-09,-6.34446e-14,-12430,-53.2908], Tmin=(1426,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 592,
    label = "CC5H10OOH-DO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {10,S}
3  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {3,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.83496,0.0860948,-7.22651e-05,3.17312e-08,-5.6201e-12,-26135.5,25.6008], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[25.2421,0.0241694,-8.19141e-06,1.26382e-09,-7.30015e-14,-33530.5,-97.5222], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)C(CO[O])OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 593,
    label = "AC5H10OOH-DO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {3,S} {8,S}
7  O u0 p2 c0 {4,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.47959,0.0837669,-6.63706e-05,2.7648e-08,-4.71461e-12,-23947.4,29.9188], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[24.9285,0.0249502,-8.57455e-06,1.33548e-09,-7.76538e-14,-31716.8,-94.6333], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CCO[O])COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 594,
    label = "C5H9D-A,COOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
3  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {2,S} {17,S} {18,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {3,S} {9,S}
8  O u0 p2 c0 {6,S} {19,S}
9  O u0 p2 c0 {7,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.13973,0.0863108,-7.36972e-05,3.27652e-08,-5.85074e-12,-18258.5,26.5675], Tmin=(298,'K'), Tmax=(1406,'K')),
            NASAPolynomial(coeffs=[25.6538,0.0234057,-7.93033e-06,1.22331e-09,-7.06529e-14,-25610.5,-96.8592], Tmin=(1406,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(OO)C(C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 595,
    label = "C5H9B-DOOH",
    molecule = 
"""
1  C u0 p0 c0 {5,S} {6,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
3  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {3,S} {5,D}
5  C u0 p0 c0 {1,S} {4,D} {16,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.549656,0.070499,-5.27826e-05,2.06471e-08,-3.33664e-12,-16463.6,28.7892], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[21.133,0.0215087,-7.42113e-06,1.15935e-09,-6.75687e-14,-23560.8,-81.5041], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)=CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 596,
    label = "C5H9OA-BOOH-D",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {3,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.90416,0.0828867,-6.64999e-05,2.72222e-08,-4.44711e-12,-29988.7,39.8117], Tmin=(298,'K'), Tmax=(1417,'K')),
            NASAPolynomial(coeffs=[21.7023,0.022612,-7.58982e-06,1.16366e-09,-6.69308e-14,-37536.5,-84.9179], Tmin=(1417,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(CCOO)CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 597,
    label = "C5H9OB-DOOH-A",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {4,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.52561,0.0853243,-6.713e-05,2.68114e-08,-4.26214e-12,-30885.6,47.0812], Tmin=(298,'K'), Tmax=(1420,'K')),
            NASAPolynomial(coeffs=[20.999,0.0235749,-7.91024e-06,1.21256e-09,-6.9736e-14,-38796.5,-82.7788], Tmin=(1420,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(CCO1)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 598,
    label = "C4H63,1-1OH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,D} {9,S}
3  C u0 p0 c0 {2,D} {4,S} {10,S}
4  C u1 p0 c0 {3,S} {5,S} {11,S}
5  O u0 p2 c0 {4,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.596641,0.0531885,-4.34792e-05,1.86514e-08,-3.23781e-12,-6155.44,27.4734], Tmin=(298,'K'), Tmax=(1401,'K')),
            NASAPolynomial(coeffs=[13.8675,0.0154796,-5.1772e-06,7.91828e-10,-4.54657e-14,-10765.5,-48.7821], Tmin=(1401,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C[CH]O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 599,
    label = "C5H9C-A,DOOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
5  C u1 p0 c0 {1,S} {4,S} {18,S}
6  O u0 p2 c0 {2,S} {9,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.40737,0.0731845,-4.34551e-05,1.06844e-08,-6.06073e-13,-17482.7,29.2965], Tmin=(298,'K'), Tmax=(1681,'K')),
            NASAPolynomial(coeffs=[26.3011,0.0238232,-8.29527e-06,1.30384e-09,-7.63134e-14,-26280.2,-101.075], Tmin=(1681,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC([CH]COO)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 600,
    label = "C5H9A-A,DOOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u1 p0 c0 {1,S} {17,S} {18,S}
6  O u0 p2 c0 {4,S} {8,S}
7  O u0 p2 c0 {3,S} {9,S}
8  O u0 p2 c0 {6,S} {19,S}
9  O u0 p2 c0 {7,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.03378,0.0863931,-7.07032e-05,3.02282e-08,-5.2511e-12,-15905.8,33.754], Tmin=(298,'K'), Tmax=(1399,'K')),
            NASAPolynomial(coeffs=[25.3362,0.0242043,-8.32306e-06,1.29685e-09,-7.54307e-14,-23822.8,-94.8747], Tmin=(1399,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(CCOO)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 601,
    label = "C4H71-1OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {11,S}
4  C u0 p0 c0 {3,D} {5,S} {12,S}
5  O u0 p2 c0 {4,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.465981,0.0556573,-4.50578e-05,1.93726e-08,-3.38968e-12,-23100.8,27.9745], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[14.2587,0.0171933,-5.74956e-06,8.7909e-10,-5.04586e-14,-27805.1,-49.6582], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 602,
    label = "CH3CHCO",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,D} {7,S}
3 C u0 p0 c0 {2,D} {8,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 O u0 p2 c0 {3,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.4838,0.0322203,-2.7025e-05,1.20499e-08,-2.18366e-12,-11527.7,17.1552], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[10.0219,0.00956966,-3.26222e-06,5.05232e-10,-2.92593e-14,-14248.3,-27.783], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""03/03/95 THERM""",
    longDesc = 
u"""
03/03/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 603,
    label = "O-OC6H5OJ",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,D} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {9,D}
4  C u0 p0 c0 {2,D} {6,S} {11,S}
5  C u0 p0 c0 {3,S} {6,D} {13,S}
6  C u0 p0 c0 {4,S} {5,D} {12,S}
7  O u1 p2 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  O u0 p2 c0 {3,D}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.65459,0.0717179,-6.31552e-05,2.81133e-08,-4.97463e-12,6452.83,38.1123], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[18.4626,0.0157607,-5.44671e-06,8.51766e-10,-4.9676e-14,-172.77,-72.8742], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""WKM""",
    longDesc = 
u"""
WKM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]C1C=CC=CC1=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 604,
    label = "C5H9OB-COOH-D",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
3  C u0 p0 c0 {2,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {2,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.03627,0.0844316,-6.63448e-05,2.5862e-08,-3.95367e-12,-32045,38.4622], Tmin=(298,'K'), Tmax=(1427,'K')),
            NASAPolynomial(coeffs=[23.366,0.0212968,-7.16312e-06,1.10028e-09,-6.33839e-14,-40271.1,-96.2545], Tmin=(1427,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(C)OC1COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 605,
    label = "C5H9OA-AOOH-D",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {6,S} {16,S} {17,S}
5  C u0 p0 c0 {2,S} {7,S} {12,S} {13,S}
6  O u0 p2 c0 {3,S} {4,S}
7  O u0 p2 c0 {5,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.39507,0.105103,-9.37686e-05,4.22174e-08,-7.50607e-12,-26432.5,70.169], Tmin=(298,'K'), Tmax=(1411,'K')),
            NASAPolynomial(coeffs=[21.6523,0.0234456,-7.95666e-06,1.22905e-09,-7.10635e-14,-35596.3,-86.9264], Tmin=(1411,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOCCC1COC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 606,
    label = "C5H9A-A,COOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {11,S}
3  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {1,S} {17,S} {18,S}
6  O u0 p2 c0 {2,S} {9,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.46543,0.0883867,-7.6201e-05,3.41129e-08,-6.11972e-12,-18105,29.7797], Tmin=(298,'K'), Tmax=(1406,'K')),
            NASAPolynomial(coeffs=[25.5933,0.0234677,-7.95419e-06,1.2273e-09,-7.0896e-14,-25609.6,-96.7335], Tmin=(1406,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(COO)C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 607,
    label = "C5H9B-AOOH",
    molecule = 
"""
1  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {13,S} {14,S} {15,S}
3  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {3,S} {4,D} {16,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.549656,0.070499,-5.27826e-05,2.06471e-08,-3.33664e-12,-16463.6,28.7892], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[21.133,0.0215087,-7.42113e-06,1.15935e-09,-6.75687e-14,-23560.8,-81.5041], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C(C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 608,
    label = "C5H9A-AOOH",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {2,S} {5,D}
5  C u0 p0 c0 {4,D} {15,S} {16,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.566437,0.0751593,-5.92588e-05,2.44465e-08,-4.14012e-12,-14944.1,34.6236], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[21.2339,0.0214696,-7.41711e-06,1.15964e-09,-6.76214e-14,-22282.5,-81.5535], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(CC)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 609,
    label = "IC5KETCB",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {4,S} {17,D}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 O u0 p2 c0 {5,D}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.499341,0.0779929,-6.02226e-05,2.34361e-08,-3.65044e-12,-46424.4,29.3526], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[24.3081,0.0203117,-6.88557e-06,1.06414e-09,-6.15901e-14,-54389,-97.6241], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)C(C)(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 610,
    label = "PC4H8OH-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
4  C u1 p0 c0 {1,S} {2,S} {13,S}
5  O u0 p2 c0 {2,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.89477,0.0487803,-3.26703e-05,1.1712e-08,-1.771e-12,-10745.5,20.6226], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[14.7199,0.0193692,-6.59357e-06,1.02019e-09,-5.90411e-14,-15302.3,-48.5296], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 611,
    label = "NC5KET15",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,D} {17,S}
6  O u0 p2 c0 {4,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.63079,0.0636042,-3.8005e-05,1.00772e-08,-8.15256e-13,-40220.1,22.7069], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[22.2257,0.0226347,-7.79348e-06,1.21574e-09,-7.07808e-14,-47426.6,-84.1025], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=CCCCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 612,
    label = "BC5H10OOH-CO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
3  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
6  O u0 p2 c0 {1,S} {8,S}
7  O u0 p2 c0 {2,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.92758,0.0837498,-7.09477e-05,3.16521e-08,-5.69638e-12,-30570.8,17.1231], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[25.1614,0.0240667,-8.11774e-06,1.2483e-09,-7.19347e-14,-37527.4,-99.5534], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(O[O])C(C)(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 613,
    label = "IC5KETDC",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {16,D} {17,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.47014,0.0930872,-8.13203e-05,3.57768e-08,-6.26401e-12,-40686.8,43.8329], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[25.746,0.0195157,-6.7026e-06,1.04493e-09,-6.08454e-14,-49655.6,-104.886], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)C(C=O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 614,
    label = "IC5KETDA",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {16,D} {17,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.3189,0.0708408,-4.98675e-05,1.7632e-08,-2.48535e-12,-40714.2,27.8334], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[22.3253,0.022301,-7.62391e-06,1.18365e-09,-6.86863e-14,-47995.4,-85.0606], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CC=O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 615,
    label = "IC5KETAD",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,D} {17,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.55453,0.0648725,-4.2145e-05,1.35396e-08,-1.69772e-12,-40993.6,22.551], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[21.4062,0.0231815,-7.94585e-06,1.2355e-09,-7.17602e-14,-47741.6,-79.4699], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C=O)CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 616,
    label = "BC5H10OOH-DO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {8,S}
7  O u0 p2 c0 {3,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.66225,0.0810932,-6.46324e-05,2.73539e-08,-4.75157e-12,-28398.3,20.3147], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[24.8033,0.0249329,-8.54e-06,1.32701e-09,-7.70334e-14,-35696.4,-97.1095], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(CCO[O])OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 617,
    label = "TC4H8OOH-IO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {2,S} {16,S}
7  O u0 p2 c0 {5,S} {17,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u1 p2 c0 {6,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.36414,0.0693743,-5.70416e-05,2.4604e-08,-4.32849e-12,-25113.8,16.5767], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[23.2465,0.0188385,-6.40938e-06,9.92649e-10,-5.75276e-14,-31653.3,-88.8302], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(CO[O])OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 618,
    label = "L-C6H4",
    molecule = 
"""
1  C u0 p0 c0 {2,D} {3,S} {7,S}
2  C u0 p0 c0 {1,D} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,T}
4  C u0 p0 c0 {3,T} {5,S}
5  C u0 p0 c0 {4,S} {6,T}
6  C u0 p0 c0 {5,T} {10,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.295902,0.0580533,-6.77668e-05,4.33768e-08,-1.14189e-11,60001.4,22.319], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[12.7152,0.0138397,-4.37654e-06,3.15416e-10,4.6619e-14,57031.1,-39.4646], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""H6W/94""",
    longDesc = 
u"""
H6W/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#CC#CC=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 619,
    label = "CH3OCH2O2H",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  O u0 p2 c0 {1,S} {2,S}
4  O u0 p2 c0 {1,S} {5,S}
5  O u0 p2 c0 {4,S} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.05787,0.0436787,-3.46384e-05,1.44809e-08,-2.46101e-12,-36885.1,24.3392], Tmin=(298,'K'), Tmax=(1404,'K')),
            NASAPolynomial(coeffs=[12.8159,0.0134818,-4.50398e-06,6.88229e-10,-3.94884e-14,-40674.6,-37.8048], Tmin=(1404,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/12/14 THERM""",
    longDesc = 
u"""
2/12/14 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
COCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 620,
    label = "IC5KETCA",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {4,S} {17,D}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 O u0 p2 c0 {5,D}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.30505,0.0612212,-3.95733e-05,1.30418e-08,-1.72022e-12,-45238.9,19.0775], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[19.9169,0.0237976,-8.01511e-06,1.2314e-09,-7.09188e-14,-51116.9,-70.578], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)C(C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 621,
    label = "IC5KETDB",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {16,D} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.45036,0.0683214,-4.82823e-05,1.7384e-08,-2.52419e-12,-45155.4,18.4815], Tmin=(298,'K'), Tmax=(1381,'K')),
            NASAPolynomial(coeffs=[22.2545,0.0221858,-7.5454e-06,1.16735e-09,-6.75735e-14,-51997,-87.8432], Tmin=(1381,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(CC=O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 622,
    label = "NC5KET21",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {3,S} {17,D}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 O u0 p2 c0 {5,D}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.197529,0.0792953,-6.01376e-05,2.26193e-08,-3.36573e-12,-41875.2,35.9472], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[24.5823,0.0202064,-6.87782e-06,1.06593e-09,-6.18172e-14,-50228.8,-96.4925], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC(=O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 623,
    label = "NC5KET23",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {4,S} {17,D}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 O u0 p2 c0 {5,D}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.459063,0.0823364,-6.692e-05,2.76548e-08,-4.58693e-12,-44025.7,36.181], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[24.2608,0.0203282,-6.88576e-06,1.06357e-09,-6.15325e-14,-52106.4,-94.9417], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(OO)C(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 624,
    label = "C4H71-4OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  O u0 p2 c0 {2,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.212276,0.0486358,-3.17242e-05,1.04987e-08,-1.39007e-12,-20036.5,30.4092], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[13.3437,0.0182064,-6.15151e-06,9.47319e-10,-5.46543e-14,-24847.9,-42.7999], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 625,
    label = "P-C6H3O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,D} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {7,D}
3  C u0 p0 c0 {1,D} {4,S} {9,S}
4  C u0 p0 c0 {3,S} {6,S} {10,D}
5  C u0 p0 c0 {2,S} {6,D} {11,S}
6  C u1 p0 c0 {4,S} {5,D}
7  O u0 p2 c0 {2,D}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 O u0 p2 c0 {4,D}
11 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.57852,0.0655376,-6.50309e-05,3.32027e-08,-6.86666e-12,15175,33.1519], Tmin=(270,'K'), Tmax=(1290,'K')),
            NASAPolynomial(coeffs=[12.2964,0.0215055,-1.07516e-05,2.57528e-09,-2.41024e-13,11542.9,-37.2584], Tmin=(1290,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (270,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""AK0505""",
    longDesc = 
u"""
AK0505.
O=C1[C]=CC(=O)C=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 626,
    label = "BC5H10OOH-C",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {1,S} {4,S} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.66397,0.0653351,-4.66109e-05,1.79783e-08,-2.91993e-12,-11840.9,16.7474], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[19.9231,0.0242353,-8.26025e-06,1.27917e-09,-7.40757e-14,-17827.3,-75.7848], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C(C)(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 627,
    label = "HOCH2O2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 O u0 p2 c0 {1,S} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 O u1 p2 c0 {2,S}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.82069,0.0247857,-1.66186e-05,4.79633e-09,-4.28088e-13,-22207.7,17.06], Tmin=(298,'K'), Tmax=(1377,'K')),
            NASAPolynomial(coeffs=[11.6406,0.00572826,-2.05362e-06,3.29071e-10,-1.95188e-14,-25350.6,-30.7332], Tmin=(1377,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 628,
    label = "CH3CHCHO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,D} {7,S}
3 C u0 p0 c0 {2,D} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 O u1 p2 c0 {3,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.47167,0.0269252,-1.00248e-05,-1.13421e-09,1.03417e-12,-4041.42,18.8722], Tmin=(298,'K'), Tmax=(1424,'K')),
            NASAPolynomial(coeffs=[10.6781,0.0112806,-3.89011e-06,6.07617e-10,-3.54121e-14,-7732.34,-32.4971], Tmin=(1424,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 629,
    label = "C5H9OB-DOOH-C",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {9,S}
3  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {3,S}
7  O u0 p2 c0 {2,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.01431,0.0855622,-6.67458e-05,2.61563e-08,-4.06383e-12,-33073.9,42.2585], Tmin=(298,'K'), Tmax=(1416,'K')),
            NASAPolynomial(coeffs=[22.3011,0.0228353,-7.73563e-06,1.19366e-09,-6.89722e-14,-41340.7,-92.1483], Tmin=(1416,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(C)OCC1OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 630,
    label = "C5H9D-B,COOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {2,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {2,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.12195,0.0844868,-7.2887e-05,3.29013e-08,-5.96275e-12,-22680,17.8948], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[25.6998,0.0232725,-7.86284e-06,1.21048e-09,-6.98108e-14,-29677.2,-100.361], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(OO)C(C)(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 631,
    label = "AC5H10OOH-BO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
6  O u0 p2 c0 {3,S} {8,S}
7  O u0 p2 c0 {1,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.66225,0.0810932,-6.46324e-05,2.73539e-08,-4.75157e-12,-28398.3,21.0042], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[24.8033,0.0249329,-8.54e-06,1.32701e-09,-7.70334e-14,-35696.4,-96.4201], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)(COO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 632,
    label = "C5H9C-B,DOOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
5  C u1 p0 c0 {1,S} {4,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.83903,0.0767535,-5.9302e-05,2.4335e-08,-4.12011e-12,-21964.5,16.4789], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[24.9289,0.0244398,-8.37848e-06,1.3027e-09,-7.56544e-14,-29031.2,-95.7845], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)([CH]COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 633,
    label = "C5H9C-BOOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {14,S}
5  C u0 p0 c0 {4,D} {15,S} {16,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {17,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.04862,0.08267,-7.17942e-05,3.18245e-08,-5.62913e-12,-19038.3,34.1117], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[23.0089,0.0191356,-6.44358e-06,9.91225e-10,-5.71817e-14,-26610.3,-92.4108], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 634,
    label = "C5H9OC-DOOH-A",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
3  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {2,S} {4,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.71788,0.0963723,-8.44994e-05,3.7085e-08,-6.39906e-12,-28365.6,53.1581], Tmin=(298,'K'), Tmax=(1424,'K')),
            NASAPolynomial(coeffs=[23.2181,0.021238,-7.10088e-06,1.0861e-09,-6.23742e-14,-36881.2,-93.0382], Tmin=(1424,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(COO)C1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 635,
    label = "C5H9OA-DOOH-C",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {3,S} {4,S}
7  O u0 p2 c0 {2,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-10.0859,0.111154,-0.000100718,4.57406e-08,-8.16325e-12,-38648.9,77.5075], Tmin=(298,'K'), Tmax=(1418,'K')),
            NASAPolynomial(coeffs=[21.7932,0.0235823,-7.96857e-06,1.22739e-09,-7.08291e-14,-48248.3,-88.7831], Tmin=(1418,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1COCC1OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 636,
    label = "C5H9OA-COOH-D",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {2,S} {3,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.23567,0.106958,-9.80175e-05,4.50653e-08,-8.12668e-12,-28961.1,69.3053], Tmin=(298,'K'), Tmax=(1416,'K')),
            NASAPolynomial(coeffs=[21.9787,0.0227762,-7.64127e-06,1.17119e-09,-6.73502e-14,-37940,-87.8801], Tmin=(1416,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1COC1COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 637,
    label = "C5H9O2-3OOH-5",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {2,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.95212,0.0890845,-7.32339e-05,2.98425e-08,-4.76173e-12,-30028.1,44.5335], Tmin=(298,'K'), Tmax=(1436,'K')),
            NASAPolynomial(coeffs=[23.731,0.0206752,-6.88444e-06,1.05045e-09,-6.02362e-14,-38436.3,-96.1946], Tmin=(1436,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1OC1CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 638,
    label = "C5H9OA-DOOH-A",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {16,S} {17,S}
4  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
6  O u0 p2 c0 {3,S} {5,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-10.5446,0.109795,-9.75108e-05,4.3619e-08,-7.70229e-12,-36452.9,81.5338], Tmin=(298,'K'), Tmax=(1412,'K')),
            NASAPolynomial(coeffs=[21.0306,0.0243652,-8.26753e-06,1.27699e-09,-7.38333e-14,-46105.5,-83.6589], Tmin=(1412,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOCC1CCOC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 639,
    label = "C5H9A-A,BOOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {1,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {3,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.87833,0.0817085,-6.64694e-05,2.85545e-08,-5.00711e-12,-20509.3,21.6843], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[25.3148,0.024105,-8.26152e-06,1.28429e-09,-7.45769e-14,-27825.1,-97.0436], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(CC)(COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 640,
    label = "CH3COCHO",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,S} {7,D}
3 C u0 p0 c0 {2,S} {8,D} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 O u0 p2 c0 {2,D}
8 O u0 p2 c0 {3,D}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.08731,0.0309032,-1.98794e-05,6.26175e-09,-7.69946e-13,-34398.9,18.284], Tmin=(298,'K'), Tmax=(1381,'K')),
            NASAPolynomial(coeffs=[11.4371,0.0106774,-3.68968e-06,5.77007e-10,-3.36532e-14,-37807.9,-32.5054], Tmin=(1381,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 641,
    label = "C5H9OA-AOOH-C",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {3,S} {4,S}
7  O u0 p2 c0 {2,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-6.41151,0.0994318,-8.58616e-05,3.70273e-08,-6.27783e-12,-28855.2,59.8124], Tmin=(298,'K'), Tmax=(1426,'K')),
            NASAPolynomial(coeffs=[22.7244,0.0220779,-7.391e-06,1.13153e-09,-6.50287e-14,-37821.7,-92.9948], Tmin=(1426,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(OO)C1COC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 642,
    label = "C5H9OA-AOOH-B",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {3,S} {4,S}
7  O u0 p2 c0 {1,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.49709,0.0947289,-8.00296e-05,3.42125e-08,-5.80023e-12,-31144.2,54.0079], Tmin=(298,'K'), Tmax=(1421,'K')),
            NASAPolynomial(coeffs=[21.7927,0.0230721,-7.77386e-06,1.19508e-09,-6.887e-14,-39671.2,-89.4814], Tmin=(1421,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1(COC1)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 643,
    label = "B-CC5H10O",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
6  O u0 p2 c0 {1,S} {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.02734,0.0775766,-6.31068e-05,2.66345e-08,-4.52085e-12,-22213.3,42.3801], Tmin=(298,'K'), Tmax=(1413,'K')),
            NASAPolynomial(coeffs=[17.2822,0.0221812,-7.45086e-06,1.14267e-09,-6.57264e-14,-28975.7,-69.9426], Tmin=(1413,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1OC1(C)C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 644,
    label = "C5H9OA-COOH-A",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {2,S} {3,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.23567,0.106958,-9.80175e-05,4.50653e-08,-8.12668e-12,-28961.1,69.3053], Tmin=(298,'K'), Tmax=(1416,'K')),
            NASAPolynomial(coeffs=[21.9787,0.0227762,-7.64127e-06,1.17119e-09,-6.73502e-14,-37940,-87.8801], Tmin=(1416,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1OCC1COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 645,
    label = "P-OC6H5OJ",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,D} {9,S}
3  C u0 p0 c0 {1,S} {5,D} {10,S}
4  C u0 p0 c0 {2,D} {6,S} {12,S}
5  C u0 p0 c0 {3,D} {6,S} {13,S}
6  C u0 p0 c0 {4,S} {5,S} {11,D}
7  O u1 p2 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 O u0 p2 c0 {6,D}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.29683,0.0727366,-6.36158e-05,2.80684e-08,-4.92279e-12,6734.02,40.935], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[18.28,0.0159281,-5.50765e-06,8.6165e-10,-5.02678e-14,-62.5908,-72.5809], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""WKM""",
    longDesc = 
u"""
WKM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]C1C=CC(=O)C=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 646,
    label = "C5H10OOH1-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {1,S} {2,S} {17,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.21642,0.0596427,-3.30151e-05,7.82802e-09,-4.53627e-13,-6890.38,25.3166], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[19.4584,0.0246032,-8.38535e-06,1.29884e-09,-7.52367e-14,-13354.9,-69.0616], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 647,
    label = "PC4H8OH-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {3,S} {13,S}
5  O u0 p2 c0 {2,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.21685,0.047134,-3.04255e-05,1.03165e-08,-1.45044e-12,-10820.9,18.3498], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[14.5945,0.0192223,-6.48869e-06,9.98309e-10,-5.75489e-14,-15239.3,-48.5128], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 648,
    label = "NC5KET12O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,D} {16,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 O u1 p2 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.82328,0.0556014,-3.02956e-05,6.22735e-09,-8.2538e-14,-23505.6,19.212], Tmin=(298,'K'), Tmax=(1377,'K')),
            NASAPolynomial(coeffs=[20.9368,0.0200874,-7.0483e-06,1.1133e-09,-6.53771e-14,-30431.4,-80.4082], Tmin=(1377,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC([O])C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 649,
    label = "A-AC5H10O",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
6  O u0 p2 c0 {3,S} {4,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.03246,0.0838427,-6.94901e-05,2.94315e-08,-4.93197e-12,-16718.4,59.1468], Tmin=(298,'K'), Tmax=(1452,'K')),
            NASAPolynomial(coeffs=[15.8864,0.0227073,-7.40384e-06,1.11239e-09,-6.30612e-14,-23720.2,-60.9114], Tmin=(1452,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1COC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 650,
    label = "NC5KET32O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {5,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {2,S} {16,D}
6  O u1 p2 c0 {2,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.55783,0.0455546,-1.73485e-05,-9.24152e-10,1.38889e-12,-27491.9,11.0597], Tmin=(298,'K'), Tmax=(1539,'K')),
            NASAPolynomial(coeffs=[19.0206,0.0210023,-7.20629e-06,1.12155e-09,-6.51937e-14,-33355.3,-69.7524], Tmin=(1539,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(=O)C(C)[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 651,
    label = "TC3H6OH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u1 p0 c0 {1,S} {2,S} {4,S}
4  O u0 p2 c0 {3,S} {11,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.0967,0.0380728,-2.75022e-05,1.07477e-08,-1.74896e-12,-14076.4,22.2476], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[11.2222,0.0136444,-4.51407e-06,7.10523e-10,-4.2269e-14,-17535,-31.8912], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/ 9/ 4 THERM""",
    longDesc = 
u"""
8/ 9/ 4 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C](C)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 652,
    label = "CH3ONO2",
    molecule = 
"""
1 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2 N u0 p0 c+1 {3,S} {7,D} {8,S}
3 O u0 p2 c0 {1,S} {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 O u0 p2 c0 {2,D}
8 O u0 p3 c-1 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.91364,0.0152138,1.73479e-05,-3.37074e-08,1.44322e-11,-16610.3,9.44208], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[9.77845,0.011007,-4.25929e-06,7.18198e-10,-4.42042e-14,-18880.4,-23.9163], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T05/98""",
    longDesc = 
u"""
T05/98.
CO[N+](=O)[O-]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 653,
    label = "NEOC5KET",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,D} {17,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.993,0.06638,-4.41434e-05,1.47656e-08,-1.97164e-12,-43172.8,23.0091], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[21.1457,0.0234263,-8.03515e-06,1.24991e-09,-7.26173e-14,-49981.3,-80.4502], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(C=O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 654,
    label = "C5H10OOH3-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {3,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.46239,0.0700774,-5.18935e-05,2.03024e-08,-3.27212e-12,-7458.96,24.9361], Tmin=(298,'K'), Tmax=(1401,'K')),
            NASAPolynomial(coeffs=[20.5595,0.0235924,-8.01737e-06,1.23919e-09,-7.16674e-14,-13905.9,-76.9647], Tmin=(1401,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(CC)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 655,
    label = "C5H10O2-4",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
3  C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
6  O u0 p2 c0 {1,S} {2,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-6.91381,0.0865695,-7.51445e-05,3.29307e-08,-5.77631e-12,-20233.7,56.6267], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[15.8873,0.0248574,-1.03963e-05,1.95193e-09,-1.29465e-13,-27222.4,-62.6943], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1CC(C)O1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 656,
    label = "SC4H8OH-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {1,S} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.78283,0.0518872,-3.89601e-05,1.57819e-08,-2.64134e-12,-12830.5,20.0308], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[15.0517,0.0184749,-6.15328e-06,9.38017e-10,-5.37214e-14,-17197.4,-50.3734], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(O)CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 657,
    label = "NEOC5H10OOH-O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {3,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.12491,0.0751548,-4.3935e-05,1.09212e-08,-8.06704e-13,-25745.9,24.7165], Tmin=(298,'K'), Tmax=(2038,'K')),
            NASAPolynomial(coeffs=[24.3461,0.0284841,-1.05961e-05,1.74545e-09,-1.05655e-13,-33534.8,-95.8566], Tmin=(2038,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""914""",
    longDesc = 
u"""
914
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(CO[O])COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 658,
    label = "B-DC5H10O",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
6  O u0 p2 c0 {1,S} {3,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.57208,0.0733363,-5.36914e-05,2.00361e-08,-2.9759e-12,-20780.8,46.0627], Tmin=(298,'K'), Tmax=(1426,'K')),
            NASAPolynomial(coeffs=[15.9398,0.0234036,-7.80111e-06,1.19023e-09,-6.82188e-14,-27566,-63.1519], Tmin=(1426,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(C)CCO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 659,
    label = "IQJC4H8OH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {16,S}
6  O u0 p2 c0 {2,S} {15,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u1 p2 c0 {6,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.81449,0.0747453,-7.10895e-05,3.44974e-08,-6.54647e-12,-34402.4,17.738], Tmin=(298,'K'), Tmax=(1410,'K')),
            NASAPolynomial(coeffs=[21.1752,0.0175144,-5.73227e-06,8.63387e-10,-4.90282e-14,-39888.2,-81.9187], Tmin=(1410,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""L 2/00""",
    longDesc = 
u"""
L 2/00
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(O)CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 660,
    label = "CVCCJCVCOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,D} {8,S}
2  C u0 p0 c0 {1,S} {4,D} {9,S}
3  C u0 p0 c0 {1,D} {5,S} {7,S}
4  C u0 p0 c0 {2,D} {6,S} {10,S}
5  C u1 p0 c0 {3,S} {11,S} {12,S}
6  O u0 p2 c0 {4,S} {13,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.91175,0.0669362,-5.71603e-05,2.48754e-08,-4.33244e-12,1964.42,41.7454], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[16.7466,0.0158357,-5.44955e-06,8.49881e-10,-4.94743e-14,-4309.73,-61.9379], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/6/95 Z&B""",
    longDesc = 
u"""
10/6/95 Z&B
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CC=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 661,
    label = "CC4H8O",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  O u0 p2 c0 {2,S} {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-6.56747,0.0787299,-7.33065e-05,3.40603e-08,-6.15675e-12,-27158.3,38.0876], Tmin=(298,'K'), Tmax=(1431,'K')),
            NASAPolynomial(coeffs=[15.1842,0.0164657,-5.33483e-06,7.9815e-10,-4.5116e-14,-33392.3,-74.3747], Tmin=(1431,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1COC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 662,
    label = "HO2CH2OCHO",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,D} {9,S}
3  O u0 p2 c0 {1,S} {2,S}
4  O u0 p2 c0 {1,S} {5,S}
5  O u0 p2 c0 {4,S} {10,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  O u0 p2 c0 {2,D}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.2191,0.0428858,-3.17634e-05,1.11543e-08,-1.49753e-12,-57928.8,24.9759], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[15.7136,0.0096443,-3.44136e-06,5.49722e-10,-3.2536e-14,-62940.9,-52.9505], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/12/14 THERM""",
    longDesc = 
u"""
2/12/14 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=COCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 663,
    label = "CH2OCH2O2H",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u1 p0 c0 {3,S} {8,S} {9,S}
3  O u0 p2 c0 {1,S} {2,S}
4  O u0 p2 c0 {1,S} {5,S}
5  O u0 p2 c0 {4,S} {10,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.162245,0.0476101,-4.52047e-05,2.18379e-08,-4.11296e-12,-14649.8,29.8253], Tmin=(298,'K'), Tmax=(1418,'K')),
            NASAPolynomial(coeffs=[12.3893,0.0111759,-3.59249e-06,5.34196e-10,-3.00537e-14,-18055.2,-32.9577], Tmin=(1418,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/12/14 THERM""",
    longDesc = 
u"""
2/12/14 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]OCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 664,
    label = "PQC4H8OS",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {3,S} {6,S}
6  O u0 p2 c0 {5,S} {16,S}
7  O u1 p2 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.6898,0.0610664,-4.46911e-05,1.68939e-08,-2.59692e-12,-21246.1,19.4272], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[19.9994,0.0196391,-6.70756e-06,1.04036e-09,-6.03188e-14,-27125.8,-73.1285], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC([O])COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 665,
    label = "IC5KETAA",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,D} {17,S}
6  O u0 p2 c0 {3,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.55453,0.0648725,-4.2145e-05,1.35396e-08,-1.69772e-12,-40993.6,22.551], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[21.4062,0.0231815,-7.94585e-06,1.2355e-09,-7.17602e-14,-47741.6,-79.4699], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C=O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 666,
    label = "IC5KETAB",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {16,D} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 O u0 p2 c0 {5,D}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.276645,0.0819422,-6.36166e-05,2.46646e-08,-3.82375e-12,-42178.9,33.5973], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[25.4564,0.0201871,-7.02057e-06,1.10298e-09,-6.45488e-14,-50866.6,-103.885], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)(C=O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 667,
    label = "IC5KETAC",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,D} {17,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {5,D}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.38365,0.0608703,-3.15373e-05,4.94037e-09,4.82023e-13,-43163.9,17.4831], Tmin=(298,'K'), Tmax=(1530,'K')),
            NASAPolynomial(coeffs=[23.6202,0.0216575,-7.50494e-06,1.17597e-09,-6.86856e-14,-50839.9,-93.7965], Tmin=(1530,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C=O)C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 668,
    label = "C4H",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,T}
2 C u0 p0 c0 {1,S} {4,T}
3 C u0 p0 c0 {1,T} {5,S}
4 C u1 p0 c0 {2,T}
5 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.02325,0.00709238,-6.07376e-09,-2.27575e-09,8.08699e-13,76238.1,-0.0694259], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.24288,0.00619368,-2.08593e-06,3.0822e-10,-1.63648e-14,75680.2,-7.21081], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""121686""",
    longDesc = 
u"""
121686
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[C]#CC#C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 669,
    label = "PC4H8CHO-2O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,D} {16,S}
6  O u0 p2 c0 {1,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
17 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.75699,0.059922,-3.84711e-05,1.16653e-08,-1.26665e-12,-25720.7,16.0218], Tmin=(298,'K'), Tmax=(1404,'K')),
            NASAPolynomial(coeffs=[21.6725,0.0203884,-6.89025e-06,1.06159e-09,-6.12789e-14,-32058.2,-80.8233], Tmin=(1404,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/8/15""",
    longDesc = 
u"""
10/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(CC=O)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 670,
    label = "PC4H8CHO-3O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,D} {16,S}
6  O u0 p2 c0 {1,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
17 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.75699,0.059922,-3.84711e-05,1.16653e-08,-1.26665e-12,-25720.7,16.0218], Tmin=(298,'K'), Tmax=(1404,'K')),
            NASAPolynomial(coeffs=[21.6725,0.0203884,-6.89025e-06,1.06159e-09,-6.12789e-14,-32058.2,-80.8233], Tmin=(1404,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/8/15""",
    longDesc = 
u"""
10/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CCC=O)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 671,
    label = "PC4H8CHO-4O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,D} {16,S}
6  O u0 p2 c0 {4,S} {17,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
17 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.37851,0.0577659,-3.28975e-05,7.80996e-09,-4.14572e-13,-23530.9,19.738], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[21.2856,0.0212771,-7.31815e-06,1.14082e-09,-6.639e-14,-30207.9,-78.203], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/8/15""",
    longDesc = 
u"""
10/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OCCCCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 672,
    label = "IC4KETII",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,D} {14,S}
5  O u0 p2 c0 {2,S} {6,S}
6  O u0 p2 c0 {5,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.15502,0.0610622,-4.49711e-05,1.70515e-08,-2.65949e-12,-38274.8,26.9612], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[19.5143,0.0182377,-6.38909e-06,1.00802e-09,-5.9144e-14,-44688.5,-71.7168], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/19/ 0 THERM""",
    longDesc = 
u"""
7/19/ 0 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C=O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 673,
    label = "C5H10OH-3O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {18,S}
7  O u0 p2 c0 {4,S} {19,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 O u1 p2 c0 {6,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.54218,0.0757844,-5.63028e-05,2.16375e-08,-3.37143e-12,-33865.1,27.3478], Tmin=(298,'K'), Tmax=(1406,'K')),
            NASAPolynomial(coeffs=[22.8162,0.024075,-8.15387e-06,1.2576e-09,-7.26293e-14,-40989.4,-86.0866], Tmin=(1406,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 8/15""",
    longDesc = 
u"""
10/ 8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(CCO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 674,
    label = "IC5KETACO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {15,D} {16,S}
6  H u0 p0 c0 {1,S}
7  O u1 p2 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.17102,0.0552253,-3.48046e-05,1.07471e-08,-1.26265e-12,-25463.5,16.0764], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[18.776,0.0209364,-7.11851e-06,1.10073e-09,-6.36813e-14,-31056.1,-68.4257], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC([O])C(C)C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 675,
    label = "BC5H10OOH-D",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {2,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.1161,0.0719622,-5.58678e-05,2.33618e-08,-4.0445e-12,-10243.5,24.0826], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[20.1548,0.023999,-8.17022e-06,1.26425e-09,-7.31719e-14,-16570.5,-77.0357], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(C)(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 676,
    label = "NEOC5H10OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u1 p0 c0 {1,S} {16,S} {17,S}
6  O u0 p2 c0 {2,S} {7,S}
7  O u0 p2 c0 {6,S} {18,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.726389,0.0654764,-3.46434e-05,6.8222e-09,-1.22217e-13,-7615.22,27.7841], Tmin=(298,'K'), Tmax=(1686,'K')),
            NASAPolynomial(coeffs=[19.8217,0.0274351,-1.01764e-05,1.67311e-09,-1.01142e-13,-14484.7,-76.55], Tmin=(1686,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/8/14""",
    longDesc = 
u"""
9/8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)(C)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 677,
    label = "C5H10OOH2-4O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {6,S} {10,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {8,S}
7  O u0 p2 c0 {2,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 O u1 p2 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.94959,0.0794697,-6.123e-05,2.45361e-08,-3.99244e-12,-27626.3,21.0066], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[25.1554,0.0242477,-8.22109e-06,1.26886e-09,-7.33158e-14,-34952.2,-96.9776], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CC(C)OO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 678,
    label = "C5H10O1-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
6  O u0 p2 c0 {1,S} {4,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.3877,0.0862948,-7.34566e-05,3.18221e-08,-5.43058e-12,-18064.9,60.9884], Tmin=(298,'K'), Tmax=(1445,'K')),
            NASAPolynomial(coeffs=[16.1537,0.0222699,-7.20968e-06,1.07781e-09,-6.08828e-14,-25121.7,-61.8744], Tmin=(1445,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1CCO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 679,
    label = "C4H7CHO-3",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {12,S}
4  C u0 p0 c0 {2,S} {3,D} {11,S}
5  C u0 p0 c0 {1,S} {13,D} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {5,D}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.41419,0.051215,-3.06682e-05,8.72169e-09,-9.15725e-13,-15155.3,28.1575], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[15.5931,0.0196345,-6.81722e-06,1.06912e-09,-6.24665e-14,-20840.8,-54.783], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/8/15 THERM""",
    longDesc = 
u"""
10/8/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 680,
    label = "NC52ONE-3",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,D} {13,S}
4  C u0 p0 c0 {2,S} {5,S} {12,D}
5  C u0 p0 c0 {3,D} {4,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 O u0 p2 c0 {4,D}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.818407,0.0517196,-3.47852e-05,1.26179e-08,-1.96289e-12,-21509.6,24.4001], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[14.8313,0.0201259,-6.9497e-06,1.08574e-09,-6.32626e-14,-26593.9,-51.4362], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 8/15 THERM""",
    longDesc = 
u"""
10/ 8/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CC(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 681,
    label = "NC52ONE-4",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {11,D}
4  C u0 p0 c0 {1,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 O u0 p2 c0 {3,D}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.246588,0.0527914,-3.58755e-05,1.29117e-08,-1.94285e-12,-17567.9,29.322], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[14.5047,0.0199955,-6.81519e-06,1.05547e-09,-6.11274e-14,-22607.2,-47.4948], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 8/15 THERM""",
    longDesc = 
u"""
10/ 8/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCC(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 682,
    label = "C4H7CHO-4",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,D} {10,S}
4  C u0 p0 c0 {2,S} {11,D} {12,S}
5  C u0 p0 c0 {3,D} {13,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 O u0 p2 c0 {4,D}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.15814,0.0481747,-2.67363e-05,6.50397e-09,-4.78417e-13,-14739,24.6792], Tmin=(298,'K'), Tmax=(2020,'K')),
            NASAPolynomial(coeffs=[13.7156,0.0216425,-7.67112e-06,1.22483e-09,-7.25907e-14,-19111.1,-43.3807], Tmin=(2020,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/8/15 THERM""",
    longDesc = 
u"""
10/8/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 683,
    label = "C4H7CHO-2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {11,S}
4  C u0 p0 c0 {3,D} {5,S} {12,S}
5  C u0 p0 c0 {4,S} {13,D} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 O u0 p2 c0 {5,D}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.000369689,0.0571383,-4.15288e-05,1.5829e-08,-2.49867e-12,-18765.6,27.587], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[16.2589,0.0187788,-6.45805e-06,1.00641e-09,-5.85466e-14,-24402.3,-59.645], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/8/15 THERM""",
    longDesc = 
u"""
10/8/15 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC=CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 684,
    label = "TQC4H7OHIO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {1,S} {8,S}
6  O u0 p2 c0 {2,S} {17,S}
7  O u0 p2 c0 {2,S} {16,S}
8  O u0 p2 c0 {5,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u1 p2 c0 {7,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.17336,0.0794006,-6.51166e-05,2.62036e-08,-4.13406e-12,-48494.3,17.7867], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[28.2565,0.016697,-5.67315e-06,8.7835e-10,-5.0909e-14,-56601.7,-115.148], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(OO)C(O)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 685,
    label = "SC4H8OH-3O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {16,S}
6  O u0 p2 c0 {2,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u1 p2 c0 {6,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.78918,0.0662988,-5.52281e-05,2.40248e-08,-4.18762e-12,-33042.8,22.8886], Tmin=(298,'K'), Tmax=(1417,'K')),
            NASAPolynomial(coeffs=[19.3083,0.0192157,-6.34251e-06,9.60913e-10,-5.4794e-14,-38437.1,-68.8888], Tmin=(1417,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(O)C(C)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 686,
    label = "SC4H8OH-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {3,S} {5,S}
5  O u0 p2 c0 {4,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.96361,0.0477674,-3.11652e-05,1.07499e-08,-1.54834e-12,-16500,19.7287], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[14.5515,0.0192722,-6.50955e-06,1.002e-09,-5.77831e-14,-20991,-48.2454], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[C](C)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 687,
    label = "CH3COCH2OCH2CH2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
6  O u0 p2 c0 {3,S} {4,S}
7  O u1 p2 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.46265,0.0531606,-3.44159e-05,1.17107e-08,-1.63587e-12,-23514.9,23.2545], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[16.1916,0.021726,-7.18256e-06,1.08968e-09,-6.22059e-14,-28331.2,-50.6779], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""14""",
    longDesc = 
u"""
14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1([O])CCOC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 688,
    label = "C5H10OH-1O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {6,S} {7,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {4,S} {19,S}
7  O u0 p2 c0 {4,S} {18,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 O u1 p2 c0 {7,S}
19 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.87267,0.0735121,-5.1977e-05,1.89528e-08,-2.82611e-12,-36160.6,23.9948], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[22.9459,0.024518,-8.42765e-06,1.31283e-09,-7.63478e-14,-43498.2,-89.257], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 8/15""",
    longDesc = 
u"""
10/ 8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC(O)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 689,
    label = "C5H10OH-2O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {18,S}
7  O u0 p2 c0 {4,S} {19,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 O u1 p2 c0 {6,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.54218,0.0757844,-5.63028e-05,2.16375e-08,-3.37143e-12,-33865.1,27.3478], Tmin=(298,'K'), Tmax=(1406,'K')),
            NASAPolynomial(coeffs=[22.8162,0.024075,-8.15387e-06,1.2576e-09,-7.26293e-14,-40989.4,-86.0866], Tmin=(1406,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 8/15""",
    longDesc = 
u"""
10/ 8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC(CO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 690,
    label = "IC5KETAAO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {15,D} {16,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  O u1 p2 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[6.20341,0.0378489,-8.21632e-06,-4.45999e-09,1.69104e-12,-23370.7,3.88282], Tmin=(298,'K'), Tmax=(1670,'K')),
            NASAPolynomial(coeffs=[15.4491,0.025883,-9.38195e-06,1.5192e-09,-9.08764e-14,-27312.9,-48.9349], Tmin=(1670,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C=O)C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 691,
    label = "C2H3CHOCH2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,D} {9,S}
4  C u0 p0 c0 {3,D} {10,S} {11,S}
5  O u0 p2 c0 {1,S} {2,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.59388,0.0579063,-4.97163e-05,2.15819e-08,-3.69199e-12,858.853,42.8475], Tmin=(298,'K'), Tmax=(1431,'K')),
            NASAPolynomial(coeffs=[12.6763,0.014082,-4.63474e-06,7.01091e-10,-3.99438e-14,-4080.65,-42.2516], Tmin=(1431,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 692,
    label = "C5H9C-A,BOOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {1,S} {4,S} {18,S}
6  O u0 p2 c0 {1,S} {9,S}
7  O u0 p2 c0 {2,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  O u0 p2 c0 {6,S} {20,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.83903,0.0767535,-5.9302e-05,2.4335e-08,-4.12011e-12,-21964.5,17.1684], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[24.9289,0.0244398,-8.37848e-06,1.3027e-09,-7.56544e-14,-29031.2,-95.095], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C(C)(COO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 693,
    label = "IC4H8O",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
3  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {2,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.02574,0.0751341,-6.88669e-05,3.12223e-08,-5.60129e-12,-30748.1,29.6284], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[14.0434,0.0205734,-9.09519e-06,1.73417e-09,-1.14909e-13,-36227.5,-69.001], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(C)CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 694,
    label = "CHOIC3H6O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,D} {13,S}
5  H u0 p0 c0 {1,S}
6  O u1 p2 c0 {3,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u0 p2 c0 {4,D}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.255559,0.0522086,-3.72284e-05,1.36714e-08,-2.05639e-12,-19728.5,29.921], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[15.5512,0.016736,-5.73573e-06,8.92095e-10,-5.18352e-14,-25067.4,-52.3216], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""L 2/00""",
    longDesc = 
u"""
L 2/00
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C=O)C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 695,
    label = "TC4H8CHO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {1,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {14,D} {15,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 O u0 p2 c0 {5,D}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.958078,0.0642003,-4.70777e-05,1.75738e-08,-2.64896e-12,-6865.83,33.3781], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[17.9664,0.0194207,-6.67409e-06,1.03969e-09,-6.04703e-14,-13336.9,-67.9819], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 7/95 THERM""",
    longDesc = 
u"""
9/ 7/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)(C)C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 696,
    label = "C5H9C-A,AOOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
4  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {1,S} {4,S} {18,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {3,S} {9,S}
8  O u0 p2 c0 {6,S} {19,S}
9  O u0 p2 c0 {7,S} {20,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.40737,0.0731845,-4.34551e-05,1.06844e-08,-6.06073e-13,-17482.7,28.607], Tmin=(298,'K'), Tmax=(1681,'K')),
            NASAPolynomial(coeffs=[26.3011,0.0238232,-8.29527e-06,1.30384e-09,-7.63134e-14,-26280.2,-101.764], Tmin=(1681,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C(COO)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 697,
    label = "PC4H8OH-2O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {15,S}
6  O u0 p2 c0 {3,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u1 p2 c0 {5,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.97781,0.0625692,-4.7673e-05,1.90562e-08,-3.10597e-12,-30962.5,23.588], Tmin=(298,'K'), Tmax=(1407,'K')),
            NASAPolynomial(coeffs=[19.0514,0.0200624,-6.7697e-06,1.04137e-09,-6.00272e-14,-36596.7,-67.1202], Tmin=(1407,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(CO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 698,
    label = "PC4H8CHO-1O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,D} {16,S}
6  O u0 p2 c0 {3,S} {17,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
17 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.463291,0.0802505,-6.47237e-05,2.62037e-08,-4.25128e-12,-23495.5,36.6761], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[24.6719,0.0185504,-6.4221e-06,1.00647e-09,-5.88171e-14,-31855.5,-97.1441], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/8/15""",
    longDesc = 
u"""
10/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC(C=O)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 699,
    label = "C5H9OA-BOOH-C",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
6  O u0 p2 c0 {1,S} {3,S}
7  O u0 p2 c0 {2,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.448273,0.0790873,-5.9708e-05,2.24184e-08,-3.31049e-12,-32328.1,31.2503], Tmin=(298,'K'), Tmax=(1417,'K')),
            NASAPolynomial(coeffs=[23.2823,0.021634,-7.33435e-06,1.13243e-09,-6.54655e-14,-40198.4,-95.1976], Tmin=(1417,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(OO)C1(C)CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 700,
    label = "C5H9OC-DOOH-B",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {9,S}
3  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {2,S} {3,S}
7  O u0 p2 c0 {1,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.31968,0.0928313,-8.14722e-05,3.59633e-08,-6.24713e-12,-32852,42.5342], Tmin=(298,'K'), Tmax=(1424,'K')),
            NASAPolynomial(coeffs=[23.1592,0.0211056,-7.01495e-06,1.0685e-09,-6.11803e-14,-40875.8,-95.8609], Tmin=(1424,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(OO)C1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 701,
    label = "SC4H8OH-1O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {16,S}
6  O u0 p2 c0 {3,S} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u1 p2 c0 {6,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.47038,0.0638295,-4.92334e-05,1.99705e-08,-3.30499e-12,-30860.2,26.3463], Tmin=(298,'K'), Tmax=(1405,'K')),
            NASAPolynomial(coeffs=[18.8547,0.0202218,-6.82329e-06,1.04962e-09,-6.05033e-14,-36572.1,-65.9105], Tmin=(1405,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(O)CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 702,
    label = "SQC4H8OP",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 O u1 p2 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.92966,0.0567741,-3.88078e-05,1.31554e-08,-1.72823e-12,-21464.9,13.2313], Tmin=(298,'K'), Tmax=(1421,'K')),
            NASAPolynomial(coeffs=[20.2482,0.0190403,-6.41819e-06,9.86792e-10,-5.68666e-14,-27068.6,-74.3735], Tmin=(1421,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C[O])OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 703,
    label = "SQC4H8OS",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {5,S} {16,S}
7  H u0 p0 c0 {1,S}
8  O u1 p2 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.45579,0.061042,-4.64107e-05,1.80687e-08,-2.81124e-12,-23483.6,13.9992], Tmin=(298,'K'), Tmax=(1424,'K')),
            NASAPolynomial(coeffs=[20.5364,0.0185262,-6.18258e-06,9.44045e-10,-5.4139e-14,-29045.4,-76.6224], Tmin=(1424,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC([O])C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 704,
    label = "TQC4H7OHIQ-P",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {14,S} {15,S}
5  O u0 p2 c0 {1,S} {8,S}
6  O u0 p2 c0 {2,S} {9,S}
7  O u0 p2 c0 {2,S} {16,S}
8  O u0 p2 c0 {5,S} {17,S}
9  O u0 p2 c0 {6,S} {18,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.36948,0.0796855,-6.74366e-05,2.82233e-08,-4.65271e-12,-40556.8,19.2726], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[28.1439,0.0163525,-5.54998e-06,8.58604e-10,-4.97359e-14,-48450.3,-111.573], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)(OO)C(O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 705,
    label = "IQC4H7OHT",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {1,S} {13,S} {14,S}
5  O u0 p2 c0 {2,S} {7,S}
6  O u0 p2 c0 {1,S} {15,S}
7  O u0 p2 c0 {5,S} {16,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.589,0.0725591,-7.1408e-05,3.55251e-08,-6.84992e-12,-25724.1,12.3767], Tmin=(298,'K'), Tmax=(1413,'K')),
            NASAPolynomial(coeffs=[21.9946,0.0162011,-5.23758e-06,7.81898e-10,-4.41126e-14,-30738.4,-81.6614], Tmin=(1413,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)(O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 706,
    label = "C6H3",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,D} {7,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {5,S}
4 C u1 p0 c0 {1,D} {8,S}
5 C u0 p0 c0 {3,S} {6,T}
6 C u0 p0 c0 {5,T} {9,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.17906,0.0555474,-7.30762e-05,5.20767e-08,-1.5047e-11,85647.3,19.1792], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.81883,0.0279334,-1.78254e-05,5.37025e-09,-6.17076e-13,85188.2,-0.921478], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""H6W/94""",
    longDesc = 
u"""
H6W/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=CC#CC#C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 707,
    label = "QCYC(CCOC)OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  O u0 p2 c0 {2,S} {3,S}
6  O u0 p2 c0 {1,S} {8,S}
7  O u0 p2 c0 {2,S} {15,S}
8  O u0 p2 c0 {6,S} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.091,0.0802657,-6.69756e-05,2.79282e-08,-4.60982e-12,-51252.9,41.1897], Tmin=(298,'K'), Tmax=(1411,'K')),
            NASAPolynomial(coeffs=[22.086,0.0183744,-6.29479e-06,9.78788e-10,-5.68621e-14,-58965.2,-86.4964], Tmin=(1411,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(COC1O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 708,
    label = "C3H6CHO-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 O u1 p2 c0 {4,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.529001,0.0430707,-2.49474e-05,6.40935e-09,-5.2777e-13,-6876.67,24.8856], Tmin=(298,'K'), Tmax=(1678,'K')),
            NASAPolynomial(coeffs=[12.173,0.0180057,-6.43783e-06,1.03362e-09,-6.1485e-14,-10864.2,-38.0322], Tmin=(1678,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC=C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 709,
    label = "IC3H6OHCHO",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,D} {13,S}
5  O u0 p2 c0 {1,S} {14,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u0 p2 c0 {4,D}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.84081,0.0529601,-3.94262e-05,1.59063e-08,-2.69565e-12,-50143.7,17.5483], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[16.0254,0.0185402,-6.36974e-06,9.91733e-10,-5.76473e-14,-55019.9,-58.3075], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(O)C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 710,
    label = "C5H9OA-BOOH-A",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {1,S} {3,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.90416,0.0828867,-6.64999e-05,2.72222e-08,-4.44711e-12,-29988.7,39.8117], Tmin=(298,'K'), Tmax=(1417,'K')),
            NASAPolynomial(coeffs=[21.7023,0.022612,-7.58982e-06,1.16366e-09,-6.69308e-14,-37536.5,-84.9179], Tmin=(1417,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/15""",
    longDesc = 
u"""
8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1(COO)CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 711,
    label = "C4H71-2OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {5,S}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  O u0 p2 c0 {3,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.28561,0.0635758,-5.85328e-05,2.80453e-08,-5.32187e-12,-25153.9,29.3803], Tmin=(298,'K'), Tmax=(1404,'K')),
            NASAPolynomial(coeffs=[15.0658,0.0165276,-5.52198e-06,8.43592e-10,-4.83871e-14,-29973.8,-55.3405], Tmin=(1404,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(O)CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 712,
    label = "COHQCYC(COC)",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {3,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {2,S} {15,S}
8  O u0 p2 c0 {6,S} {16,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.40688,0.059518,-3.16914e-05,4.23695e-09,8.42033e-13,-51969.5,18.8941], Tmin=(298,'K'), Tmax=(1319,'K')),
            NASAPolynomial(coeffs=[24.4599,0.0164188,-5.74024e-06,9.0571e-10,-5.31829e-14,-60181.1,-102.05], Tmin=(1319,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(CO1)C(O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 713,
    label = "C2H5CHCO",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {11,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 O u0 p2 c0 {4,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-22.8307,0.170978,-0.000353394,2.78222e-07,-6.77325e-11,-10412.5,131.233], Tmin=(298,'K'), Tmax=(1550,'K')),
            NASAPolynomial(coeffs=[-204.041,0.293467,-0.000115885,1.95254e-08,-1.19031e-12,82738,1212.33], Tmin=(1550,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC=C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 714,
    label = "PC3H4OH-2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {3,S} {5,S} {6,S} {7,S}
2 C u0 p0 c0 {3,D} {4,S} {8,S}
3 C u1 p0 c0 {1,S} {2,D}
4 O u0 p2 c0 {2,S} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.42757,0.0364826,-3.18007e-05,1.46915e-08,-2.72331e-12,7803.43,18.589], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[10.7164,0.0106066,-3.51374e-06,5.33714e-10,-3.04902e-14,4984.87,-29.8329], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""4/ 2/13 THERM""",
    longDesc = 
u"""
4/ 2/13 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C]=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 715,
    label = "IC3H5COHQ",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {6,S} {8,S}
2  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {1,S} {14,S}
7  O u0 p2 c0 {5,S} {15,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.64361,0.0572485,-3.95908e-05,1.27776e-08,-1.47241e-12,-38010,17.7322], Tmin=(298,'K'), Tmax=(1504,'K')),
            NASAPolynomial(coeffs=[20.7388,0.0158361,-5.27614e-06,8.05932e-10,-4.62682e-14,-44182.1,-79.424], Tmin=(1504,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)C(O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 716,
    label = "CCY(CCOC)OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  O u0 p2 c0 {2,S} {3,S}
6  O u0 p2 c0 {1,S} {14,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.84915,0.0623737,-4.70327e-05,1.84908e-08,-2.94976e-12,-37425.7,38.3164], Tmin=(298,'K'), Tmax=(1404,'K')),
            NASAPolynomial(coeffs=[14.3405,0.0198312,-6.60658e-06,1.0078e-09,-5.77614e-14,-43095.9,-53.0535], Tmin=(1404,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""L 2/00""",
    longDesc = 
u"""
L 2/00
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(O)COC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 717,
    label = "TQC4H7OHIQ-I",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {5,S} {17,S}
8  O u0 p2 c0 {6,S} {18,S}
9  O u1 p2 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[6.09882,0.0693745,-5.12049e-05,1.81276e-08,-2.4633e-12,-39139.4,2.93863], Tmin=(298,'K'), Tmax=(1384,'K')),
            NASAPolynomial(coeffs=[28.8467,0.016629,-5.74302e-06,8.98791e-10,-5.24795e-14,-46924.9,-119.118], Tmin=(1384,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(OO)C([O])OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 718,
    label = "SQC4H7OHP-4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u1 p0 c0 {2,S} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {3,S} {15,S}
7  O u0 p2 c0 {5,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.44518,0.0652944,-5.1038e-05,2.05239e-08,-3.32051e-12,-22906,27.8546], Tmin=(298,'K'), Tmax=(1410,'K')),
            NASAPolynomial(coeffs=[19.9207,0.0189423,-6.39488e-06,9.84231e-10,-5.67606e-14,-28925.1,-70.1061], Tmin=(1410,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(CO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 719,
    label = "SQC4H7OHS-4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {5,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u1 p0 c0 {2,S} {13,S} {14,S}
5  O u0 p2 c0 {2,S} {7,S}
6  O u0 p2 c0 {1,S} {15,S}
7  O u0 p2 c0 {5,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.37953,0.0650949,-5.5911e-05,2.49644e-08,-4.44897e-12,-25150.2,15.8597], Tmin=(298,'K'), Tmax=(1414,'K')),
            NASAPolynomial(coeffs=[20.6664,0.0178772,-5.93194e-06,9.02006e-10,-5.15692e-14,-30413.2,-74.4595], Tmin=(1414,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(OO)C(C)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 720,
    label = "C5H9OH1-2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {14,S}
5  C u0 p0 c0 {3,S} {4,D} {15,S}
6  O u0 p2 c0 {3,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.71305,0.0505618,-2.46539e-05,4.87082e-09,-1.89221e-13,-23958.9,23.1278], Tmin=(298,'K'), Tmax=(1375,'K')),
            NASAPolynomial(coeffs=[16.2865,0.0240392,-8.40886e-06,1.32516e-09,-7.76838e-14,-29943.5,-58.0986], Tmin=(1375,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 8/15""",
    longDesc = 
u"""
10/ 8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC=CCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 721,
    label = "PC3H4OH-1",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2 C u0 p0 c0 {1,S} {3,D} {8,S}
3 C u1 p0 c0 {2,D} {4,S}
4 O u0 p2 c0 {3,S} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.0715,0.0350017,-3.01922e-05,1.38468e-08,-2.5538e-12,6819.14,15.1917], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[10.9469,0.0104015,-3.44082e-06,5.22106e-10,-2.98049e-14,4115.31,-31.1162], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/25/15""",
    longDesc = 
u"""
9/25/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=[C]O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 722,
    label = "O2CH2OCH2O2H",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
3  O u0 p2 c0 {1,S} {2,S}
4  O u0 p2 c0 {1,S} {6,S}
5  O u0 p2 c0 {2,S} {11,S}
6  O u0 p2 c0 {4,S} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 O u1 p2 c0 {5,S}
12 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.39978,0.0539882,-4.8797e-05,2.19792e-08,-3.86107e-12,-33782.5,23.0683], Tmin=(298,'K'), Tmax=(1433,'K')),
            NASAPolynomial(coeffs=[17.7378,0.011359,-3.67383e-06,5.49256e-10,-3.10406e-14,-38290.3,-56.661], Tmin=(1433,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""2/12/14 ERM""",
    longDesc = 
u"""
2/12/14 ERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OCOCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 723,
    label = "NC4KET21OH",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {2,S} {13,D}
5  O u0 p2 c0 {2,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.46677,0.030112,-2.59549e-06,-7.61701e-09,2.55573e-12,-48990.9,6.33308], Tmin=(298,'K'), Tmax=(1507,'K')),
            NASAPolynomial(coeffs=[15.2253,0.0186615,-6.30659e-06,9.72356e-10,-5.6177e-14,-53476,-50.0397], Tmin=(1507,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(=O)CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 724,
    label = "C2H4COCH2OH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {1,S} {12,D}
5  O u0 p2 c0 {2,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u0 p2 c0 {4,D}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[7.98874,0.0141152,1.863e-05,-1.99223e-08,5.12984e-12,-29870.6,-1.87192], Tmin=(298,'K'), Tmax=(1462,'K')),
            NASAPolynomial(coeffs=[14.0828,0.0176833,-6.09759e-06,9.52795e-10,-5.55572e-14,-33751.8,-40.8944], Tmin=(1462,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC([C]=O)CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 725,
    label = "C4H63,1-2OH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {3,D} {4,S} {5,S}
3  C u0 p0 c0 {1,S} {2,D} {9,S}
4  C u1 p0 c0 {2,S} {10,S} {11,S}
5  O u0 p2 c0 {2,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.58904,0.0617723,-5.78264e-05,2.7809e-08,-5.2678e-12,-8179.77,29.6999], Tmin=(298,'K'), Tmax=(1404,'K')),
            NASAPolynomial(coeffs=[14.69,0.0148008,-4.94484e-06,7.55543e-10,-4.33462e-14,-12944.8,-54.5653], Tmin=(1404,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(O)=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 726,
    label = "C5H9O23-5",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {5,S} {12,S} {13,S} {14,S}
5  C u1 p0 c0 {1,S} {4,S} {15,S}
6  O u0 p2 c0 {1,S} {3,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.67636,0.0774809,-6.42244e-05,2.6844e-08,-4.43868e-12,4979.27,49.7878], Tmin=(298,'K'), Tmax=(1424,'K')),
            NASAPolynomial(coeffs=[17.7402,0.0192298,-6.43473e-06,9.84654e-10,-5.65625e-14,-2054.48,-68.2409], Tmin=(1424,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C1CCO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 727,
    label = "PQC4H7OHS-3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {3,S} {14,S}
5  O u0 p2 c0 {2,S} {7,S}
6  O u0 p2 c0 {1,S} {15,S}
7  O u0 p2 c0 {5,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.14059,0.0594343,-4.55415e-05,1.85507e-08,-3.10109e-12,-24076.5,19.9412], Tmin=(298,'K'), Tmax=(1404,'K')),
            NASAPolynomial(coeffs=[18.9704,0.0196512,-6.61523e-06,1.01591e-09,-5.8489e-14,-29287.6,-64.0683], Tmin=(1404,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C(O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 728,
    label = "CCY(COCC)OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {6,S} {8,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {3,S}
6  O u0 p2 c0 {2,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-6.13261,0.0775376,-6.73066e-05,2.82842e-08,-4.6547e-12,-35782.5,55.3089], Tmin=(298,'K'), Tmax=(1328,'K')),
            NASAPolynomial(coeffs=[10.3817,0.0253053,-1.00512e-05,1.68205e-09,-9.58093e-14,-39434.9,-27.3523], Tmin=(1328,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1OCC1O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 729,
    label = "OC5H7O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u1 p0 c0 {1,S} {5,S} {10,S}
4  C u0 p0 c0 {2,S} {11,D} {12,S}
5  C u0 p0 c0 {3,S} {13,D} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 O u0 p2 c0 {4,D}
12 H u0 p0 c0 {4,S}
13 O u0 p2 c0 {5,D}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.88395,0.0403401,-1.97774e-05,3.68904e-09,-3.40202e-14,-23529.6,9.9707], Tmin=(298,'K'), Tmax=(1375,'K')),
            NASAPolynomial(coeffs=[16.5417,0.0186678,-6.44836e-06,1.00788e-09,-5.87522e-14,-28201.7,-54.7258], Tmin=(1375,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/22/ 9 WKM""",
    longDesc = 
u"""
1/22/ 9 WKM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=C[CH]CCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 730,
    label = "CY(CCCO)COH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
5  O u0 p2 c0 {1,S} {3,S}
6  O u0 p2 c0 {4,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-6.9759,0.0798393,-7.11132e-05,3.16612e-08,-5.61324e-12,-33577.8,61.2093], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[13.4636,0.0228746,-9.9362e-06,1.89605e-09,-1.26467e-13,-39639.7,-45.1064], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OCC1CCO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 731,
    label = "C3KET13",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,D} {11,S}
4  O u0 p2 c0 {2,S} {5,S}
5  O u0 p2 c0 {4,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 O u0 p2 c0 {3,D}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.74957,0.0314081,-6.83838e-06,-5.67124e-09,2.27687e-12,-35192.5,9.83754], Tmin=(298,'K'), Tmax=(1508,'K')),
            NASAPolynomial(coeffs=[17.3613,0.0132331,-4.75332e-06,7.62529e-10,-4.52614e-14,-40624.8,-61.7768], Tmin=(1508,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/17/12""",
    longDesc = 
u"""
10/17/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=CCCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 732,
    label = "NEOC5H9Q2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {1,S} {17,S} {18,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {3,S} {9,S}
8  O u0 p2 c0 {6,S} {19,S}
9  O u0 p2 c0 {7,S} {20,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.90419,0.0768024,-4.67667e-05,1.25175e-08,-1.11212e-12,-17739.3,27.5076], Tmin=(298,'K'), Tmax=(2050,'K')),
            NASAPolynomial(coeffs=[24.8113,0.0276462,-1.03018e-05,1.69885e-09,-1.02913e-13,-25666.3,-96.4119], Tmin=(2050,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)(COO)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 733,
    label = "IC4H8OH-TI",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u1 p0 c0 {1,S} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {14,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.33169,0.0513017,-4.02699e-05,1.7515e-08,-3.16002e-12,-14831.9,15.5368], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[14.6324,0.0188896,-6.30561e-06,9.62474e-10,-5.5164e-14,-18797.6,-49.3219], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)(C)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 734,
    label = "NEOC5H9O-OOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  O u0 p2 c0 {2,S} {3,S}
7  O u0 p2 c0 {4,S} {8,S}
8  O u0 p2 c0 {7,S} {18,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.33879,0.0823095,-5.24635e-05,1.40859e-08,-1.00626e-12,-28763.1,50.4066], Tmin=(298,'K'), Tmax=(1429,'K')),
            NASAPolynomial(coeffs=[23.4301,0.0226474,-7.84731e-06,1.22995e-09,-7.18634e-14,-38637.8,-100.093], Tmin=(1429,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/14""",
    longDesc = 
u"""
9/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1(COO)COC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 735,
    label = "IC5KETABO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {1,S} {15,D} {16,S}
6  O u1 p2 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.44493,0.0539988,-2.76353e-05,5.14987e-09,-2.28549e-14,-24253.9,19.1348], Tmin=(298,'K'), Tmax=(1680,'K')),
            NASAPolynomial(coeffs=[17.4004,0.0243634,-8.89218e-06,1.44673e-09,-8.68336e-14,-29635.9,-62.6151], Tmin=(1680,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""17/8/15""",
    longDesc = 
u"""
17/8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(C)([O])C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 736,
    label = "PC4H8OH-4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u1 p0 c0 {2,S} {12,S} {13,S}
5  O u0 p2 c0 {3,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.0695635,0.0520263,-3.44849e-05,1.17704e-08,-1.63047e-12,-10605.4,31.591], Tmin=(298,'K'), Tmax=(1401,'K')),
            NASAPolynomial(coeffs=[14.152,0.019644,-6.64538e-06,1.02403e-09,-5.91002e-14,-15614.1,-45.0604], Tmin=(1401,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 737,
    label = "CH2CQCOHQ",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {3,D} {5,S}
3  C u0 p0 c0 {2,D} {10,S} {11,S}
4  O u0 p2 c0 {1,S} {7,S}
5  O u0 p2 c0 {2,S} {8,S}
6  O u0 p2 c0 {1,S} {12,S}
7  O u0 p2 c0 {4,S} {14,S}
8  O u0 p2 c0 {5,S} {13,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-16.1172,0.178866,-0.000216751,1.1545e-07,-2.27163e-11,-52116.5,114.441], Tmin=(298,'K'), Tmax=(1418,'K')),
            NASAPolynomial(coeffs=[38.6574,0.000483815,-2.10414e-07,3.81491e-11,-2.46188e-15,-65639.2,-161.084], Tmin=(1418,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/ 1/14""",
    longDesc = 
u"""
7/ 1/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(OO)C(O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 738,
    label = "HOC3H6O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  O u0 p2 c0 {1,S} {12,S}
5  O u0 p2 c0 {2,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u1 p2 c0 {4,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.8496,0.0477245,-3.60393e-05,1.4348e-08,-2.33508e-12,-28210.6,17.6479], Tmin=(298,'K'), Tmax=(1407,'K')),
            NASAPolynomial(coeffs=[15.6948,0.0157704,-5.30502e-06,8.14308e-10,-4.68666e-14,-32454.1,-50.6084], Tmin=(1407,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 1/12""",
    longDesc = 
u"""
9/ 1/12
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(CO)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 739,
    label = "C4H71-3OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,D} {10,S}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  O u0 p2 c0 {1,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.0716708,0.0504773,-3.50357e-05,1.30539e-08,-2.07783e-12,-21698.4,28.2073], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[14.2402,0.0183039,-6.37213e-06,1.001e-09,-5.85511e-14,-26818.4,-48.395], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 740,
    label = "IC4H6Q2-II",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
2  C u0 p0 c0 {3,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {2,S} {4,D}
4  C u0 p0 c0 {3,D} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {2,S} {8,S}
7  O u0 p2 c0 {5,S} {15,S}
8  O u0 p2 c0 {6,S} {16,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.216408,0.0825572,-7.6454e-05,3.58212e-08,-6.67719e-12,-20569.5,33.6943], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[25.0361,0.016023,-5.70705e-06,9.10412e-10,-5.38295e-14,-28419.6,-96.7186], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8/14""",
    longDesc = 
u"""
9/ 8/14
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(COO)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 741,
    label = "C8H131-5,3,PAO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {7,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {8,D} {18,S}
6  C u0 p0 c0 {2,S} {7,D} {19,S}
7  C u0 p0 c0 {4,S} {6,D} {20,S}
8  C u0 p0 c0 {5,D} {21,S} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 O u1 p2 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.91532,0.071032,-3.55624e-05,6.98595e-09,-1.83466e-13,7229.12,16.7146], Tmin=(298,'K'), Tmax=(1374,'K')),
            NASAPolynomial(coeffs=[25.1903,0.031563,-1.10568e-05,1.74424e-09,-1.02327e-13,-1344.02,-101.448], Tmin=(1374,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(C)CC=CC[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 742,
    label = "HOCH2CHO",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,D} {7,S}
3 O u0 p2 c0 {1,S} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 O u0 p2 c0 {2,D}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.35084,0.0150413,-5.95083e-07,-3.75736e-09,1.09354e-12,-39165.5,8.0961], Tmin=(298,'K'), Tmax=(1680,'K')),
            NASAPolynomial(coeffs=[7.99599,0.0120664,-4.43693e-06,7.25167e-10,-4.36535e-14,-40902.1,-13.3754], Tmin=(1680,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/24/15""",
    longDesc = 
u"""
9/24/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=CCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 743,
    label = "SQC4H7OHP-4O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {8,S}
6  O u0 p2 c0 {4,S} {16,S}
7  O u0 p2 c0 {3,S} {17,S}
8  O u0 p2 c0 {5,S} {18,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 O u1 p2 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.62618,0.0726813,-5.64957e-05,2.25892e-08,-3.65247e-12,-40900,26.5615], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[23.5936,0.0208588,-7.13602e-06,1.10821e-09,-6.43133e-14,-47845.2,-84.9506], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OCCC(CO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 744,
    label = "C4H72-2OH",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,D} {5,S}
4  C u0 p0 c0 {2,S} {3,D} {12,S}
5  O u0 p2 c0 {3,S} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.184034,0.0590105,-5.22507e-05,2.44136e-08,-4.56519e-12,-26671.8,23.6075], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[14.9412,0.016629,-5.55628e-06,8.48899e-10,-4.86949e-14,-31249.7,-55.1731], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C(C)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 745,
    label = "IQC4H7OHTO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {2,S} {8,S}
6  O u0 p2 c0 {1,S} {17,S}
7  O u0 p2 c0 {3,S} {16,S}
8  O u0 p2 c0 {5,S} {18,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u1 p2 c0 {7,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.83287,0.0690926,-5.15031e-05,1.99236e-08,-3.17457e-12,-43444.3,19.4649], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[24.172,0.0209397,-7.29061e-06,1.14554e-09,-6.70218e-14,-50477.2,-89.5997], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(O)(CO[O])COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 746,
    label = "PQC4H7OHS-3O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {3,S} {8,S}
6  O u0 p2 c0 {2,S} {16,S}
7  O u0 p2 c0 {1,S} {17,S}
8  O u0 p2 c0 {5,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u1 p2 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.90987,0.0776596,-6.68514e-05,2.93628e-08,-5.12008e-12,-43155.2,22.9131], Tmin=(298,'K'), Tmax=(1414,'K')),
            NASAPolynomial(coeffs=[24.6608,0.0192948,-6.45434e-06,9.87246e-10,-5.6689e-14,-49846.7,-91.0425], Tmin=(1414,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(O[O])C(O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 747,
    label = "SQC4H7OHS-4O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {1,S} {8,S}
6  O u0 p2 c0 {2,S} {17,S}
7  O u0 p2 c0 {3,S} {16,S}
8  O u0 p2 c0 {5,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u1 p2 c0 {7,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.90987,0.0776596,-6.68514e-05,2.93628e-08,-5.12008e-12,-43155.2,22.9131], Tmin=(298,'K'), Tmax=(1414,'K')),
            NASAPolynomial(coeffs=[24.6608,0.0192948,-6.45434e-06,9.87246e-10,-5.6689e-14,-49846.7,-91.0425], Tmin=(1414,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(O)C(CO[O])OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 748,
    label = "C6H2",
    molecule = 
"""
1 C u0 p0 c0 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {4,S}
3 C u0 p0 c0 {1,S} {5,T}
4 C u0 p0 c0 {2,S} {6,T}
5 C u0 p0 c0 {3,T} {7,S}
6 C u0 p0 c0 {4,T} {8,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.59326,0.0805301,-0.000148006,1.33e-07,-4.53323e-11,83273.2,27.9809], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[13.2263,0.00739043,-2.27154e-06,2.58752e-10,-5.53567e-15,80565.3,-41.2012], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""P 1/93""",
    longDesc = 
u"""
P 1/93
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#CC#CC#C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 749,
    label = "HOCH2CO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u1 p0 c0 {1,S} {6,D}
3 O u0 p2 c0 {1,S} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 O u0 p2 c0 {2,D}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.12917,0.00907173,6.49228e-06,-8.56592e-09,2.31071e-12,-20615.2,5.73008], Tmin=(298,'K'), Tmax=(1487,'K')),
            NASAPolynomial(coeffs=[9.43497,0.00768897,-2.74959e-06,4.39772e-10,-2.60488e-14,-22970,-20.4619], Tmin=(1487,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/25/15""",
    longDesc = 
u"""
9/25/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=[C]CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 750,
    label = "NC4KET13OH-2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {6,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,D} {14,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {2,S} {15,S}
7  O u0 p2 c0 {5,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.13332,0.0669841,-5.19672e-05,2.03864e-08,-3.22075e-12,-59096.3,24.0965], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[22.5434,0.0176377,-6.14926e-06,9.67201e-10,-5.66322e-14,-65972.4,-84.8523], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(OO)C(O)C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 751,
    label = "NC4KET24OH-3",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {4,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {2,S} {14,D}
5  O u0 p2 c0 {2,S} {7,S}
6  O u0 p2 c0 {1,S} {15,S}
7  O u0 p2 c0 {5,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 O u0 p2 c0 {4,D}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.6098,0.061153,-4.50409e-05,1.73313e-08,-2.75574e-12,-60770.8,24.4892], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[20.2237,0.0194541,-6.73928e-06,1.05529e-09,-6.15931e-14,-66865.1,-69.9677], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(O)C(=O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 752,
    label = "IQC4H7OHTQ-P",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {14,S} {15,S}
5  O u0 p2 c0 {2,S} {8,S}
6  O u0 p2 c0 {3,S} {9,S}
7  O u0 p2 c0 {1,S} {16,S}
8  O u0 p2 c0 {5,S} {17,S}
9  O u0 p2 c0 {6,S} {18,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.45243,0.0714541,-5.54617e-05,2.22755e-08,-3.66067e-12,-35413.3,22.3024], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[24.577,0.0201889,-7.03585e-06,1.10623e-09,-6.47523e-14,-42578.3,-90.5053], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(O)(COO)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 753,
    label = "NC4KET24OH-1",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {3,S} {14,D}
5  O u0 p2 c0 {2,S} {7,S}
6  O u0 p2 c0 {3,S} {15,S}
7  O u0 p2 c0 {5,S} {16,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 O u0 p2 c0 {4,D}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[7.17383,0.0372831,-7.3e-06,-5.42424e-09,1.94299e-12,-59147,3.98041], Tmin=(298,'K'), Tmax=(1672,'K')),
            NASAPolynomial(coeffs=[17.664,0.0238922,-8.87612e-06,1.46065e-09,-8.83477e-14,-63656.1,-56.0405], Tmin=(1672,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=C(CO)CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 754,
    label = "C7H111-5,3,6P",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {6,D} {14,S}
5  C u0 p0 c0 {2,S} {7,D} {15,S}
6  C u0 p0 c0 {4,D} {16,S} {17,S}
7  C u1 p0 c0 {5,D} {18,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.676849,0.0738083,-5.35945e-05,2.06294e-08,-3.29701e-12,32816,34.0518], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[19.5787,0.0254214,-8.68127e-06,1.3463e-09,-7.80465e-14,25853.8,-74.4112], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=CCC(C)C=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 755,
    label = "CHOC(CH3)OHCH2Q",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,D} {14,S}
5  O u0 p2 c0 {2,S} {7,S}
6  O u0 p2 c0 {1,S} {15,S}
7  O u0 p2 c0 {5,S} {16,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.40934,0.0613326,-4.03875e-05,1.22738e-08,-1.32997e-12,-58305,22.9871], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[22.4481,0.0178754,-6.26905e-06,9.90084e-10,-5.81417e-14,-65509.3,-85.6747], Tmin=(1383,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(O)(C=O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 756,
    label = "C4H6OHOOH2-2-1",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {4,D} {6,S}
4  C u0 p0 c0 {2,S} {3,D} {13,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {3,S} {14,S}
7  O u0 p2 c0 {5,S} {15,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.0472577,0.0785787,-7.68006e-05,3.79206e-08,-7.33204e-12,-35114,28.0065], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[21.6656,0.0159965,-5.52215e-06,8.6274e-10,-5.02768e-14,-41487.5,-83.9378], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C(O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 757,
    label = "C4H6OHOOH1-3-4",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {11,S}
4  C u0 p0 c0 {3,D} {12,S} {13,S}
5  O u0 p2 c0 {2,S} {7,S}
6  O u0 p2 c0 {1,S} {14,S}
7  O u0 p2 c0 {5,S} {15,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.944278,0.0634331,-4.99456e-05,2.05721e-08,-3.48144e-12,-31782.7,29.2924], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[19.3445,0.0182186,-6.34492e-06,9.97098e-10,-5.8342e-14,-37995.2,-68.8145], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 758,
    label = "C4H6OHOOH1-4-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {11,S}
4  C u0 p0 c0 {3,D} {12,S} {13,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {2,S} {14,S}
7  O u0 p2 c0 {5,S} {15,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.90016,0.0772727,-6.79668e-05,2.98158e-08,-5.16343e-12,-30165,40.681], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[21.5164,0.0155812,-5.2757e-06,8.14811e-10,-4.71425e-14,-37482.6,-82.4101], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(CO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 759,
    label = "C4H6OHOOH1-2-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,D} {5,S}
4  C u0 p0 c0 {3,D} {12,S} {13,S}
5  O u0 p2 c0 {3,S} {7,S}
6  O u0 p2 c0 {1,S} {14,S}
7  O u0 p2 c0 {5,S} {15,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.88476,0.0916099,-9.37152e-05,4.66407e-08,-8.92736e-12,-35290,39.2852], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[23.2793,0.0138281,-4.61331e-06,7.0532e-10,-4.05175e-14,-42634.8,-95.198], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(OO)C(C)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 760,
    label = "C4H6O1-3OOH4",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {11,S}
4  C u0 p0 c0 {3,D} {12,S} {13,S}
5  O u0 p2 c0 {2,S} {6,S}
6  O u0 p2 c0 {5,S} {14,S}
7  O u1 p2 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.83046,0.0507388,-3.42226e-05,1.13486e-08,-1.48529e-12,-5617.86,15.151], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[19.5456,0.0160359,-5.61769e-06,8.86342e-10,-5.20073e-14,-11267.9,-69.949], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC([O])COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 761,
    label = "IQC4H8OTQ-I",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  O u0 p2 c0 {2,S} {7,S}
6  O u0 p2 c0 {3,S} {8,S}
7  O u0 p2 c0 {5,S} {17,S}
8  O u0 p2 c0 {6,S} {18,S}
9  O u1 p2 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.9022,0.0626952,-4.17352e-05,1.35422e-08,-1.7117e-12,-33962,7.8929], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[25.3249,0.0200635,-7.0125e-06,1.10475e-09,-6.47562e-14,-40961.1,-97.3593], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC([O])(COO)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 762,
    label = "IC3H5Q",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,D} {4,S}
3  C u0 p0 c0 {2,D} {9,S} {10,S}
4  O u0 p2 c0 {2,S} {5,S}
5  O u0 p2 c0 {4,S} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.32903,0.0449171,-3.51235e-05,1.41982e-08,-2.33335e-12,-12189.8,20.2697], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[14.3424,0.0128054,-4.40585e-06,6.86848e-10,-3.99675e-14,-16526.1,-48.9935], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 763,
    label = "C5H9OH1-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {15,S}
5  C u0 p0 c0 {3,S} {4,D} {14,S}
6  O u0 p2 c0 {2,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.188733,0.0571778,-3.44215e-05,9.99837e-09,-1.07316e-12,-24356.1,29.9036], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[16.2749,0.0227984,-7.70752e-06,1.18748e-09,-6.85336e-14,-30238,-57.5828], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 8/15""",
    longDesc = 
u"""
10/ 8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CCCO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 764,
    label = "C5H9OH1-1",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {5,D} {14,S}
5  C u0 p0 c0 {4,D} {6,S} {15,S}
6  O u0 p2 c0 {5,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.752968,0.0684304,-5.40649e-05,2.26641e-08,-3.88275e-12,-26026,31.0301], Tmin=(298,'K'), Tmax=(1401,'K')),
            NASAPolynomial(coeffs=[17.58,0.0215511,-7.24423e-06,1.11154e-09,-6.39603e-14,-31984.3,-65.9854], Tmin=(1401,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/ 8/15""",
    longDesc = 
u"""
10/ 8/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 765,
    label = "C4H63,1-3OH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,D} {5,S}
3  C u0 p0 c0 {2,D} {4,S} {9,S}
4  C u1 p0 c0 {3,S} {10,S} {11,S}
5  O u0 p2 c0 {2,S} {12,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.596155,0.0435893,-2.65985e-05,8.1425e-09,-1.02316e-12,-5056.31,24.3924], Tmin=(298,'K'), Tmax=(1380,'K')),
            NASAPolynomial(coeffs=[13.3633,0.0170223,-5.95245e-06,9.37896e-10,-5.49767e-14,-9901.13,-45.4753], Tmin=(1380,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=C(C)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 766,
    label = "CO(CH2OOH)2",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {2,S} {12,D}
4  O u0 p2 c0 {1,S} {6,S}
5  O u0 p2 c0 {2,S} {7,S}
6  O u0 p2 c0 {4,S} {13,S}
7  O u0 p2 c0 {5,S} {14,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 O u0 p2 c0 {3,D}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.47627,0.0893737,-9.25891e-05,4.63168e-08,-8.933e-12,-43892.4,48.4479], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[24.3376,0.0114074,-4.08932e-06,6.55183e-10,-3.88571e-14,-51686.3,-90.1518], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=C(COO)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 767,
    label = "C4H6O25",
    molecule = 
"""
1  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {4,D} {10,S}
4  C u0 p0 c0 {1,S} {3,D} {11,S}
5  O u0 p2 c0 {1,S} {2,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.67053,0.00492586,8.86967e-05,-1.26219e-07,5.23991e-11,-14657.2,14.5722], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.60658,0.020831,-8.42229e-06,1.56718e-09,-1.09391e-13,-17617.7,-23.2465], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""T 3/97""",
    longDesc = 
u"""
T 3/97.
C1=CCOC1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 768,
    label = "C4H7O2-1,3OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {3,S} {8,S}
7  O u0 p2 c0 {5,S} {17,S}
8  O u0 p2 c0 {6,S} {18,S}
9  H u0 p0 c0 {1,S}
10 O u1 p2 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.88339,0.0718219,-5.86792e-05,2.42824e-08,-3.98539e-12,-33657.6,12.4696], Tmin=(298,'K'), Tmax=(1425,'K')),
            NASAPolynomial(coeffs=[25.4691,0.0186441,-6.23826e-06,9.54423e-10,-5.48154e-14,-40153.2,-96.039], Tmin=(1425,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(OO)C([O])COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 769,
    label = "AC4H7OOH",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {6,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,D} {12,S}
4  C u0 p0 c0 {3,D} {13,S} {14,S}
5  O u0 p2 c0 {1,S} {6,S}
6  O u0 p2 c0 {2,S} {5,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.33471,0.0527831,-3.58861e-05,1.32495e-08,-2.06619e-12,-8878.92,24.3857], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[14.7661,0.0212235,-7.09403e-06,1.08424e-09,-6.22146e-14,-13561.7,-47.7449], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""6/17/13 THERM""",
    longDesc = 
u"""
6/17/13 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCOOC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 770,
    label = "CH2COHCH2OOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,D} {5,S}
3  C u0 p0 c0 {2,D} {9,S} {10,S}
4  O u0 p2 c0 {1,S} {6,S}
5  O u0 p2 c0 {2,S} {11,S}
6  O u0 p2 c0 {4,S} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.389823,0.0701531,-7.42037e-05,3.84181e-08,-7.63556e-12,-30788,28.6874], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[18.7971,0.0112783,-3.90789e-06,6.12065e-10,-3.57305e-14,-36115.5,-69.4914], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 771,
    label = "SC3H4OH",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {4,S}
2 C u1 p0 c0 {1,S} {5,S} {6,S}
3 C u0 p0 c0 {1,D} {7,S} {8,S}
4 O u0 p2 c0 {1,S} {9,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.72871,0.0441016,-4.72014e-05,2.52074e-08,-5.13376e-12,2227.21,14.3928], Tmin=(298,'K'), Tmax=(1407,'K')),
            NASAPolynomial(coeffs=[12.0968,0.00943977,-3.10774e-06,4.69609e-10,-2.67166e-14,-385.855,-37.6796], Tmin=(1407,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""3/28/13""",
    longDesc = 
u"""
3/28/13
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=C)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 772,
    label = "O2CCHOOJ",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,D}
2 C u0 p0 c0 {1,S} {5,S} {6,D}
3 O u0 p2 c0 {1,S} {7,S}
4 O u0 p2 c0 {1,D}
5 O u1 p2 c0 {2,S}
6 O u0 p2 c0 {2,D}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[8.91498,0.00860572,5.24417e-07,-2.79301e-09,7.62963e-13,-34086.8,-8.72978], Tmin=(298,'K'), Tmax=(1682,'K')),
            NASAPolynomial(coeffs=[10.9911,0.00746986,-2.75568e-06,4.51353e-10,-2.72109e-14,-35133.5,-21.1652], Tmin=(1682,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""Z&B""",
    longDesc = 
u"""
Z&B
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]C(=O)C(=O)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 773,
    label = "O2C4H8CHO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,D} {16,S}
6  O u0 p2 c0 {2,S} {17,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 O u0 p2 c0 {5,D}
16 H u0 p0 c0 {5,S}
17 O u1 p2 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.91848,0.0667246,-4.80871e-05,1.78589e-08,-2.71164e-12,-24983.8,23.8578], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[21.263,0.0214072,-7.38343e-06,1.15282e-09,-6.71508e-14,-31685.5,-79.9829], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 7/95 THERM""",
    longDesc = 
u"""
9/ 7/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(C=O)CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 774,
    label = "HOCH2COCH2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,D}
3  C u1 p0 c0 {2,S} {8,S} {9,S}
4  O u0 p2 c0 {1,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  O u0 p2 c0 {2,D}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.98635,0.0251767,-1.10893e-05,1.12011e-09,2.89491e-13,-24004.8,6.38259], Tmin=(298,'K'), Tmax=(1363,'K')),
            NASAPolynomial(coeffs=[12.7107,0.0115514,-3.95952e-06,6.16237e-10,-3.58328e-14,-27153.8,-36.7129], Tmin=(1363,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=O)CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 775,
    label = "TQC4H7OHTO2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {2,S} {8,S}
6  O u0 p2 c0 {1,S} {16,S}
7  O u0 p2 c0 {2,S} {17,S}
8  O u0 p2 c0 {5,S} {18,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u1 p2 c0 {6,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.17336,0.0794006,-6.51166e-05,2.62036e-08,-4.13406e-12,-48494.3,17.7867], Tmin=(298,'K'), Tmax=(1402,'K')),
            NASAPolynomial(coeffs=[28.2565,0.016697,-5.67315e-06,8.7835e-10,-5.0909e-14,-56601.7,-115.148], Tmin=(1402,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(O[O])C(O)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 776,
    label = "O2HC4H8CO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  O u0 p2 c0 {2,S} {7,S}
6  C u1 p0 c0 {1,S} {16,D}
7  O u0 p2 c0 {5,S} {17,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 O u0 p2 c0 {6,D}
17 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.82607,0.0693466,-4.93125e-05,1.69848e-08,-2.26118e-12,-24657.8,24.1168], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[23.822,0.0191411,-6.67919e-06,1.05127e-09,-6.15877e-14,-32309.4,-94.2581], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 7/95 THERM""",
    longDesc = 
u"""
9/ 7/95 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)([C]=O)COO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 777,
    label = "HOCOCQ(CH3)2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {6,S} {14,D}
5  O u0 p2 c0 {1,S} {7,S}
6  O u0 p2 c0 {4,S} {15,S}
7  O u0 p2 c0 {5,S} {16,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 O u0 p2 c0 {4,D}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.56326,0.0677166,-5.146e-05,1.99185e-08,-3.16007e-12,-73094.3,23.7255], Tmin=(298,'K'), Tmax=(1380,'K')),
            NASAPolynomial(coeffs=[22.8935,0.0178628,-6.34628e-06,1.01069e-09,-5.96886e-14,-80540.1,-90.8954], Tmin=(1380,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)(OO)C(=O)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 778,
    label = "NC4KET12OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {12,D} {13,S}
5  O u0 p2 c0 {1,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 O u0 p2 c0 {4,D}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.52027,0.0540317,-3.57074e-05,1.18022e-08,-1.57486e-12,-46778.1,28.9999], Tmin=(298,'K'), Tmax=(1384,'K')),
            NASAPolynomial(coeffs=[16.8651,0.0183319,-6.4139e-06,1.01105e-09,-5.92859e-14,-52735.3,-59.7188], Tmin=(1384,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(O)C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 779,
    label = "NC4KET23OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  C u0 p0 c0 {4,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {3,S} {13,D}
5  O u0 p2 c0 {1,S} {14,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
14 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.64761,0.0474122,-2.72509e-05,7.25756e-09,-6.68683e-13,-50494.6,23.78], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[15.4636,0.019212,-6.64681e-06,1.03989e-09,-6.06553e-14,-55729.3,-51.9147], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/27/15""",
    longDesc = 
u"""
7/27/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)C(C)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 780,
    label = "C5H9O14-5",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {11,S}
3  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
4  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
5  C u1 p0 c0 {2,S} {14,S} {15,S}
6  O u0 p2 c0 {2,S} {4,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.49219,0.0807086,-6.60337e-05,2.72424e-08,-4.4175e-12,-3237.5,63.602], Tmin=(298,'K'), Tmax=(1447,'K')),
            NASAPolynomial(coeffs=[15.4611,0.020681,-6.66142e-06,9.93099e-10,-5.60118e-14,-10332.7,-57.0023], Tmin=(1447,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C1CCCO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 781,
    label = "C5H9O12-5",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
5  C u1 p0 c0 {3,S} {14,S} {15,S}
6  O u0 p2 c0 {1,S} {4,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.74571,0.0755875,-6.08535e-05,2.46895e-08,-3.96233e-12,7086.54,51.5669], Tmin=(298,'K'), Tmax=(1428,'K')),
            NASAPolynomial(coeffs=[17.1693,0.019644,-6.56271e-06,1.00314e-09,-5.75798e-14,120.373,-64.1594], Tmin=(1428,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC1CO1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 782,
    label = "C5H9O23-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {6,S} {8,S}
3  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {11,S} {12,S} {13,S}
5  C u1 p0 c0 {2,S} {14,S} {15,S}
6  O u0 p2 c0 {1,S} {2,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.67636,0.0774809,-6.42244e-05,2.6844e-08,-4.43868e-12,4979.27,49.7878], Tmin=(298,'K'), Tmax=(1424,'K')),
            NASAPolynomial(coeffs=[17.7402,0.0192298,-6.43473e-06,9.84654e-10,-5.65625e-14,-2054.48,-68.2409], Tmin=(1424,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C1OC1CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 783,
    label = "H2NO",
    molecule = 
"""
multiplicity 2
1 N u1 p0 c+1 {2,S} {3,S} {4,S}
2 O u0 p3 c-1 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.53059,0.00859604,-5.47103e-06,2.27625e-09,-4.64807e-13,6868.03,11.2665], Tmin=(298,'K'), Tmax=(1500,'K')),
            NASAPolynomial(coeffs=[5.67335,0.00229884,-1.77445e-07,-1.10348e-10,1.85976e-14,5569.32,-6.15354], Tmin=(1500,'K'), Tmax=(4000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (4000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[NH2+][O-]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 784,
    label = "C2H5CHOHCO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,D} {5,S}
4  C u0 p0 c0 {3,D} {11,S} {12,S}
5  O u0 p2 c0 {3,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 O u1 p2 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.44749,0.0646308,-5.10458e-05,1.99851e-08,-3.12214e-12,-26667,39.5567], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[18.2917,0.0148007,-5.24975e-06,8.3522e-10,-4.92942e-14,-33670.4,-71.2364], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(O)=C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 785,
    label = "CH3COCOHCH3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {9,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {4,D} {5,S}
4  C u0 p0 c0 {2,S} {3,D} {12,S}
5  O u0 p2 c0 {3,S} {13,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 O u1 p2 c0 {4,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.65608,0.060841,-4.85215e-05,1.97794e-08,-3.26305e-12,-30515.8,35.8912], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[16.5625,0.0159037,-5.54554e-06,8.72275e-10,-5.10735e-14,-36588.5,-61.0834], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC([O])=C(C)O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 786,
    label = "CH3COHCO",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2 C u0 p0 c0 {1,S} {3,S} {4,D}
3 O u0 p2 c0 {2,S} {9,S}
4 C u0 p0 c0 {2,D} {8,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
8 O u0 p2 c0 {4,D}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.825855,0.0522217,-5.66992e-05,2.94816e-08,-5.81715e-12,-35015.7,17.8489], Tmin=(298,'K'), Tmax=(1410,'K')),
            NASAPolynomial(coeffs=[15.0111,0.00726697,-2.42872e-06,3.71247e-10,-2.13069e-14,-38735.5,-54.1006], Tmin=(1410,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/24/15""",
    longDesc = 
u"""
9/24/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(O)=C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 787,
    label = "C4H5C2H",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,D}
2  C u0 p0 c0 {1,S} {4,D} {7,S}
3  C u0 p0 c0 {1,S} {5,D} {10,S}
4  C u0 p0 c0 {2,D} {5,S} {8,S}
5  C u0 p0 c0 {3,D} {4,S} {9,S}
6  C u0 p0 c0 {1,D} {11,S} {12,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.718132,0.0379343,1.13988e-05,-4.13335e-08,1.80559e-11,24223.8,27.8557], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[11.1035,0.0206007,-7.53022e-06,1.23887e-09,-7.5416e-14,20361.8,-36.6652], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""0""",
    longDesc = 
u"""
0.
Duplicate of species FULVENE (i.e. same molecular structure according to RMG)
C=C1C=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 788,
    label = "C3H2",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 C u2 p0 c0 {2,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.43417,0.0173013,-1.18294e-05,1.02756e-09,1.62626e-12,76907.5,12.1012], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.67325,0.00557729,-1.9918e-06,3.20289e-10,-1.91216e-14,75757.1,-9.72894], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T12/00""",
    longDesc = 
u"""
T12/00.
[C]=C=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 789,
    label = "H2CCC(S)",
    molecule = 
"""
multiplicity 3
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 C u0 p0 c0 {1,S} {5,T}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 C u1 p0 c0 {2,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.72297,0.00925899,-2.30062e-06,-1.02008e-09,4.53744e-13,64877.3,5.68659], Tmin=(200,'K'), Tmax=(1500,'K')),
            NASAPolynomial(coeffs=[6.48888,0.00531128,-1.78095e-06,2.72526e-10,-1.56196e-14,63661.9,-10.0643], Tmin=(1500,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""0""",
    longDesc = 
u"""
0.
[C]#C[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 790,
    label = "C3H2(S)",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 C u2 p0 c0 {2,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.29765,0.0169875,-2.42665e-05,1.86537e-08,-5.5763e-12,67240.5,-3.754], Tmin=(200,'K'), Tmax=(900,'K')),
            NASAPolynomial(coeffs=[7.76426,0.00471128,-1.61706e-06,2.54724e-10,-1.50386e-14,66849.7,-15.0985], Tmin=(900,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""0""",
    longDesc = 
u"""
0.
Duplicate of species C3H2 (i.e. same molecular structure according to RMG)
[C]=C=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 791,
    label = "H2CC",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u2 p0 c0 {1,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.28155,0.00697648,-2.38552e-06,-1.21044e-09,9.81895e-13,48621.8,5.92039], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.27803,0.00475628,-1.6301e-06,2.54628e-10,-1.48864e-14,48316.7,0.640237], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L12/89""",
    longDesc = 
u"""
L12/89.
[C]=C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 792,
    label = "C2O",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,T} {3,S}
2 C u1 p0 c0 {1,T}
3 O u1 p2 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.36885,0.0082418,-8.76514e-06,5.56926e-09,-1.54001e-12,33170.8,6.71331], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.84981,0.00294758,-1.09073e-06,1.79256e-10,-1.11576e-14,32820.6,-0.645323], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""121286""",
    longDesc = 
u"""
121286
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[C]#C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 793,
    label = "SC3H5OCH2-1",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,D} {9,S}
3  C u0 p0 c0 {2,D} {5,S} {10,S}
4  C u1 p0 c0 {5,S} {11,S} {12,S}
5  O u0 p2 c0 {3,S} {4,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.235694,0.0467368,-3.04881e-05,9.75216e-09,-1.23281e-12,7116.68,28.1379], Tmin=(298,'K'), Tmax=(1382,'K')),
            NASAPolynomial(coeffs=[14.7022,0.0155342,-5.45701e-06,8.62545e-10,-5.06734e-14,1812.95,-50.512], Tmin=(1382,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]OC=CC
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 794,
    label = "NH3",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.30178,-0.00477127,2.19342e-05,-2.29856e-08,8.28992e-12,-6748.06,-0.690644], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.7171,0.00556856,-1.76886e-06,2.67417e-10,-1.52731e-14,-6584.52,6.0929], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Anharmonic    RUS 89""",
    longDesc = 
u"""
Anharmonic    RUS 89.
N
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 795,
    label = "N2H4",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 N u0 p1 c0 {1,S} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.83472,-0.00064913,3.76848e-05,-5.00709e-08,2.03362e-11,10089.4,5.75272], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.93957,0.00875017,-2.99399e-06,4.67278e-10,-2.73069e-14,9282.66,-2.6944], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""HYDRAZINE    L 5/90""",
    longDesc = 
u"""
HYDRAZINE    L 5/90.
NN
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 796,
    label = "NNH",
    molecule = 
"""
multiplicity 2
1 N u0 p1 c0 {2,D} {3,S}
2 N u1 p1 c0 {1,D}
3 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.25475,-0.00345098,1.37789e-05,-1.33264e-08,4.41023e-12,28793.2,3.28552], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.42744,0.00323295,-1.17296e-06,1.90508e-10,-1.14492e-14,28767.6,6.39209], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 8/11""",
    longDesc = 
u"""
T 8/11.
N=[N]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 797,
    label = "N2H2",
    molecule = 
"""
1 N u0 p1 c0 {2,D} {3,S}
2 N u0 p1 c0 {1,D} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.91066,-0.0107792,3.86516e-05,-3.86502e-08,1.34852e-11,22824.2,0.0910273], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[1.31115,0.00900187,-3.14912e-06,4.8145e-10,-2.71898e-14,23386.3,16.4091], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""equil & transT 9/11""",
    longDesc = 
u"""
equil & transT 9/11.
N=N
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 798,
    label = "N2O",
    molecule = 
"""
1 N u0 p0 c+1 {2,D} {3,D}
2 N u0 p2 c-1 {1,D}
3 O u0 p2 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.54306,0.00949219,-9.79278e-06,6.26384e-09,-1.90183e-12,8765.1,9.511], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.71898,0.00287371,-1.1975e-06,2.2506e-10,-1.575e-14,8165.8,-1.657], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[N+](=[N-])=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 799,
    label = "OC4H6O",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,D} {11,S}
4  C u0 p0 c0 {2,S} {10,D} {12,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  O u0 p2 c0 {3,D}
10 O u0 p2 c0 {4,D}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.21629,0.0357423,-2.04226e-05,5.63821e-09,-5.88889e-13,-37205.6,10.2815], Tmin=(298,'K'), Tmax=(1382,'K')),
            NASAPolynomial(coeffs=[14.1895,0.0153346,-5.24595e-06,8.14655e-10,-4.72759e-14,-41000.2,-44.3772], Tmin=(1382,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/23/ 9 WKM""",
    longDesc = 
u"""
1/23/ 9 WKM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=CCCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 800,
    label = "OC4H5O",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {9,D} {10,S}
4  C u1 p0 c0 {2,S} {11,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  O u0 p2 c0 {3,D}
10 H u0 p0 c0 {3,S}
11 O u0 p2 c0 {4,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.60551,0.0330499,-2.13102e-05,7.37021e-09,-1.08289e-12,-18546.1,10.1599], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[13.2139,0.0137339,-4.6264e-06,7.10941e-10,-4.09538e-14,-21653.5,-36.4185], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/23/ 9 WKM""",
    longDesc = 
u"""
1/23/ 9 WKM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=[C]CCC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 801,
    label = "H2C4O",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {5,S} {6,S}
2 C u0 p0 c0 {1,D} {3,D}
3 C u0 p0 c0 {2,D} {4,D}
4 C u0 p0 c0 {3,D} {7,D}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 O u0 p2 c0 {4,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.81097,0.01314,9.86507e-07,-6.12072e-09,1.64e-12,25458,2.11342], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[10.2689,0.00489616,-4.88508e-07,-2.70857e-10,5.10701e-14,23469,-28.1598], Tmin=(1000,'K'), Tmax=(4000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (4000,'K'),
    ),
    shortDesc = u"""120189""",
    longDesc = 
u"""
120189
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C=C=C=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 802,
    label = "HOCVCCJVO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {5,S}
2 C u0 p0 c0 {1,D} {4,S} {6,S}
3 C u1 p0 c0 {1,S} {7,D}
4 O u0 p2 c0 {2,S} {8,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 O u0 p2 c0 {3,D}
8 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.60727,0.0496011,-5.32301e-05,2.68393e-08,-5.13095e-12,-15881.5,19.4817], Tmin=(298,'K'), Tmax=(1414,'K')),
            NASAPolynomial(coeffs=[15.2721,0.00502586,-1.68409e-06,2.58391e-10,-1.48849e-14,-19850.7,-55.4642], Tmin=(1414,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/26/ 9 WKM""",
    longDesc = 
u"""
1/26/ 9 WKM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=[C]C=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 803,
    label = "HOCVCCVO",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,S} {5,S}
2 C u0 p0 c0 {1,D} {4,S} {6,S}
3 C u0 p0 c0 {1,S} {7,D} {8,S}
4 O u0 p2 c0 {2,S} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 O u0 p2 c0 {3,D}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.01837,0.062654,-6.73359e-05,3.3943e-08,-6.48918e-12,-33136.8,31.8163], Tmin=(298,'K'), Tmax=(1413,'K')),
            NASAPolynomial(coeffs=[16.6505,0.00611745,-2.09081e-06,3.24986e-10,-1.88875e-14,-38218,-63.6795], Tmin=(1413,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/26/ 9 WKM""",
    longDesc = 
u"""
1/26/ 9 WKM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=CC=CO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 804,
    label = "HCOH",
    molecule = 
"""
multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.82157,0.0357332,-3.80862e-05,1.86206e-08,-3.45958e-12,11295.7,34.8488], Tmin=(298,'K'), Tmax=(1398,'K')),
            NASAPolynomial(coeffs=[9.18749,0.00152011,-6.27604e-07,1.09728e-10,-6.89655e-15,7813.65,-27.3434], Tmin=(1398,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""MAR94""",
    longDesc = 
u"""
MAR94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 805,
    label = "CHCHO",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u1 p0 c0 {1,D} {5,S}
3 O u1 p2 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.33257,0.0162953,-9.72052e-06,5.15124e-10,1.03837e-12,29658.5,13.9905], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.92633,0.00971712,-5.54856e-06,1.53069e-09,-1.64742e-13,28949.9,0.527875], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
[CH]=C[O]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 806,
    label = "N2H3",
    molecule = 
"""
multiplicity 2
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 N u1 p1 c0 {1,S} {5,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.42126,0.00134902,2.23459e-05,-2.99728e-08,1.20979e-11,25819.9,7.83176], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.04484,0.0073113,-2.47626e-06,3.83733e-10,-2.23108e-14,25324.1,2.88423], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Rad.       T 7/11""",
    longDesc = 
u"""
Rad.       T 7/11.
N[NH]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 807,
    label = "C5H10OOH1-3OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {4,S} {8,S}
7  O u0 p2 c0 {1,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.35467,0.0827825,-6.38328e-05,2.5859e-08,-4.28053e-12,-50462,32.7492], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[23.1193,0.0259987,-8.8286e-06,1.36396e-09,-7.88607e-14,-57989.2,-88.1999], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/14/15""",
    longDesc = 
u"""
10/14/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(O)CCOO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 808,
    label = "C5H10OOH2-4OH",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {6,S} {10,S}
2  C u0 p0 c0 {3,S} {4,S} {7,S} {9,S}
3  C u0 p0 c0 {1,S} {2,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {8,S}
7  O u0 p2 c0 {2,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.769421,0.0847003,-6.89331e-05,2.93939e-08,-5.06044e-12,-52655.6,28.874], Tmin=(298,'K'), Tmax=(1411,'K')),
            NASAPolynomial(coeffs=[23.5163,0.0251436,-8.42079e-06,1.28868e-09,-7.40108e-14,-59855.8,-90.9143], Tmin=(1411,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/14/15""",
    longDesc = 
u"""
10/14/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(O)CC(C)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 809,
    label = "C5H10OOH3-1OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  O u0 p2 c0 {1,S} {8,S}
7  O u0 p2 c0 {4,S} {19,S}
8  O u0 p2 c0 {6,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.753539,0.0816717,-6.11918e-05,2.3649e-08,-3.69941e-12,-50545.8,30.5276], Tmin=(298,'K'), Tmax=(1406,'K')),
            NASAPolynomial(coeffs=[23.8717,0.0253397,-8.59895e-06,1.32803e-09,-7.6771e-14,-58270.4,-92.6846], Tmin=(1406,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/14/15""",
    longDesc = 
u"""
10/14/15
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(CCO)OO
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 810,
    label = "NC5DIONE13",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {2,S} {13,D}
5  C u0 p0 c0 {2,S} {14,D} {15,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {4,D}
14 O u0 p2 c0 {5,D}
15 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.16338,0.0532068,-2.90298e-05,5.99808e-09,-7.60468e-14,-44671.4,25.3181], Tmin=(298,'K'), Tmax=(1380,'K')),
            NASAPolynomial(coeffs=[18.1672,0.0195478,-6.78557e-06,1.06428e-09,-6.21973e-14,-51113.5,-68.0388], Tmin=(1380,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/15  THERM""",
    longDesc = 
u"""
10/15  THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC(=O)CC=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 811,
    label = "NC5DIONE24",
    molecule = 
"""
1  C u0 p0 c0 {4,S} {5,S} {6,S} {7,S}
2  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {2,S} {14,D}
5  C u0 p0 c0 {1,S} {3,S} {15,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 O u0 p2 c0 {4,D}
15 O u0 p2 c0 {5,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.23574,0.0544885,-3.57895e-05,1.21944e-08,-1.7181e-12,-48405.8,24.4807], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[16.326,0.0207678,-7.12494e-06,1.10834e-09,-6.43883e-14,-53836.1,-57.1668], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""10/15  THERM""",
    longDesc = 
u"""
10/15  THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(=O)CC(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 812,
    label = "HCCCHCCH",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {6,S}
2 C u0 p0 c0 {1,D} {4,D}
3 C u0 p0 c0 {1,S} {5,T}
4 C u1 p0 c0 {2,D} {7,S}
5 C u0 p0 c0 {3,T} {8,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {4,S}
8 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.49152,0.0303631,-2.5254e-05,1.08322e-08,-1.86418e-12,65975.1,7.59522], Tmin=(298,'K'), Tmax=(1399,'K')),
            NASAPolynomial(coeffs=[12.0901,0.00797136,-2.65858e-06,4.0607e-10,-2.33012e-14,63246.9,-37.7215], Tmin=(1399,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""160387""",
    longDesc = 
u"""
160387
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=C=CC#C
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 813,
    label = "C5H2",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,D} {3,D}
2 C u0 p0 c0 {1,D} {4,D}
3 C u0 p0 c0 {1,D} {5,D}
4 C u1 p0 c0 {2,D} {6,S}
5 C u1 p0 c0 {3,D} {7,S}
6 H u0 p0 c0 {4,S}
7 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.2934,0.0341865,-3.10647e-05,1.37482e-08,-2.37886e-12,81220.8,10.0728], Tmin=(298,'K'), Tmax=(1403,'K')),
            NASAPolynomial(coeffs=[13.1045,0.0055835,-1.96405e-06,3.10896e-10,-1.82888e-14,77864.7,-46.6964], Tmin=(1403,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""160387""",
    longDesc = 
u"""
160387
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=C=C=C=[CH]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

entry(
    index = 814,
    label = "H2CCCCCH",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,S} {6,S} {7,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {4,S}
4 C u0 p0 c0 {3,S} {5,T}
5 C u0 p0 c0 {4,T} {8,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.49152,0.0303631,-2.5254e-05,1.08322e-08,-1.86418e-12,65975.1,7.59522], Tmin=(298,'K'), Tmax=(1399,'K')),
            NASAPolynomial(coeffs=[12.0901,0.00797136,-2.65858e-06,4.0607e-10,-2.33012e-14,63246.9,-37.7215], Tmin=(1399,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""160387""",
    longDesc = 
u"""
160387
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#CC#C[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/pentaneNO/therm.dat.
""",
)

