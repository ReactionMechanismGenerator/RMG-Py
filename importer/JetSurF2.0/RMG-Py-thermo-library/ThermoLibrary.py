#!/usr/bin/env python
# encoding: utf-8

name = "/home/alongd/Code/RMG-Py/importer/JetSurF2.0"
shortDesc = u"/home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt"
longDesc = u"""
Unknown source
"""
entry(
    index = 1,
    label = "C12H24",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {7,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {17,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
5  C u0 p0 c0 {4,S} {6,S} {23,S} {24,S}
6  C u0 p0 c0 {5,S} {8,S} {25,S} {26,S}
7  C u0 p0 c0 {1,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {6,S} {10,S} {27,S} {28,S}
9  C u0 p0 c0 {7,S} {11,S} {29,S} {30,S}
10 C u0 p0 c0 {8,S} {31,S} {32,S} {33,S}
11 C u0 p0 c0 {9,S} {12,D} {34,S}
12 C u0 p0 c0 {11,D} {35,S} {36,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {4,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {5,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {6,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {12,S}
36 H u0 p0 c0 {12,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.96343,0.143992,-9.61384e-05,3.30174e-08,-4.62398e-12,-24634.5,52.9159], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[37.4002,0.0526231,-1.78624e-05,2.7595e-09,-1.59562e-13,-38940.6,-164.893], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/22/ 7 THERM""",
    longDesc = 
u"""
1/22/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCCCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 2,
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
            NASAPolynomial(coeffs=[1.02347,0.0349592,-2.2009e-05,6.94227e-09,-7.87919e-13,18118,19.7507], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[17.8156,-0.0042575,1.05118e-05,-4.47384e-09,5.84814e-13,12673.4,-69.8266], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""A 8/83""",
    longDesc = 
u"""
A 8/83
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C=CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 3,
    label = "PXC12H25",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {8,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {17,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
5  C u0 p0 c0 {4,S} {6,S} {23,S} {24,S}
6  C u0 p0 c0 {5,S} {7,S} {25,S} {26,S}
7  C u0 p0 c0 {6,S} {9,S} {27,S} {28,S}
8  C u0 p0 c0 {1,S} {10,S} {13,S} {14,S}
9  C u0 p0 c0 {7,S} {11,S} {29,S} {30,S}
10 C u0 p0 c0 {8,S} {12,S} {31,S} {32,S}
11 C u0 p0 c0 {9,S} {33,S} {34,S} {35,S}
12 C u1 p0 c0 {10,S} {36,S} {37,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {4,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {5,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {6,S}
27 H u0 p0 c0 {7,S}
28 H u0 p0 c0 {7,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {11,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {12,S}
37 H u0 p0 c0 {12,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.85029,0.142671,-9.18917e-05,3.00883e-08,-3.97454e-12,-15453,49.3702], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[38.0922,0.0542108,-1.84206e-05,2.84762e-09,-1.64732e-13,-29819.4,-166.883], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 4,
    label = "SXC5H11",
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
            NASAPolynomial(coeffs=[-0.479207,0.0574354,-4.84093e-05,2.5387e-08,-6.16554e-12,18841.4,30.043], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[13.6419,0.018558,-5.11222e-06,7.02796e-10,-3.91545e-14,14672.4,-43.577], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERG""",
    longDesc = 
u"""
1/ 2/ 7 THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 5,
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
            NASAPolynomial(coeffs=[5.7154,-0.0152309,6.52441e-05,-7.10807e-08,2.61353e-11,-25642.8,-1.5041], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[1.78971,0.0140938,-6.36501e-06,1.38171e-09,-1.1706e-13,-25374.9,14.5024], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L 8/88""",
    longDesc = 
u"""
L 8/88.
CO
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 6,
    label = "PXC12H23",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {17,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
5  C u0 p0 c0 {4,S} {7,S} {23,S} {24,S}
6  C u0 p0 c0 {1,S} {9,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {8,S} {25,S} {26,S}
8  C u0 p0 c0 {7,S} {10,S} {29,S} {30,S}
9  C u0 p0 c0 {6,S} {11,S} {27,S} {28,S}
10 C u0 p0 c0 {8,S} {12,D} {31,S}
11 C u1 p0 c0 {9,S} {32,S} {33,S}
12 C u0 p0 c0 {10,D} {34,S} {35,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {4,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {5,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {8,S}
30 H u0 p0 c0 {8,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {11,S}
33 H u0 p0 c0 {11,S}
34 H u0 p0 c0 {12,S}
35 H u0 p0 c0 {12,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.04898,0.141285,-9.41858e-05,3.21335e-08,-4.4665e-12,-7644.66,52.2139], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[37.1516,0.0509575,-1.7432e-05,2.70699e-09,-1.57088e-13,-21983.3,-164.992], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/22/ 7 THERM""",
    longDesc = 
u"""
1/22/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCCCCCCC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 7,
    label = "SXC12H25",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {7,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {17,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
5  C u0 p0 c0 {4,S} {6,S} {23,S} {24,S}
6  C u0 p0 c0 {5,S} {8,S} {25,S} {26,S}
7  C u0 p0 c0 {1,S} {9,S} {13,S} {14,S}
8  C u0 p0 c0 {6,S} {10,S} {27,S} {28,S}
9  C u0 p0 c0 {7,S} {12,S} {29,S} {30,S}
10 C u0 p0 c0 {8,S} {31,S} {32,S} {33,S}
11 C u0 p0 c0 {12,S} {34,S} {35,S} {36,S}
12 C u1 p0 c0 {9,S} {11,S} {37,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {4,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {5,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {6,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {11,S}
37 H u0 p0 c0 {12,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.36787,0.137355,-8.24076e-05,2.36422e-08,-2.47436e-12,-16766.1,48.3522], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[37.9688,0.0538719,-1.82171e-05,2.80775e-09,-1.62108e-13,-31214.5,-165.806], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCCCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 8,
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
            NASAPolynomial(coeffs=[4.16343,-0.000232616,3.42678e-05,-4.41052e-08,1.72756e-11,-2657.45,7.34683], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.94477,0.00786672,-2.88659e-06,4.72709e-10,-2.85999e-14,-3787.31,-5.01368], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 9/92""",
    longDesc = 
u"""
T 9/92.
C[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 9,
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
            NASAPolynomial(coeffs=[3.28254,0.00148309,-7.57967e-07,2.09471e-10,-2.16718e-14,-1088.46,5.45323], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""TPIS89""",
    longDesc = 
u"""
TPIS89.
[O][O]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 10,
    label = "C6H4CH3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u1 p0 c0 {2,S} {5,S} {11,S}
4  C u0 p0 c0 {2,D} {6,S} {12,S}
5  C u0 p0 c0 {3,S} {7,D} {13,S}
6  C u0 p0 c0 {4,S} {7,D} {14,S}
7  C u0 p0 c0 {5,D} {6,D}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.14159,0.0567231,-8.68851e-06,-3.42496e-08,1.92669e-11,35738.5,39.7428], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[11.6155,0.0274318,-1.08993e-05,1.86418e-09,-1.01916e-13,31209.3,-38.9946], Tmin=(1000,'K'), Tmax=(2500,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2500,'K'),
    ),
    shortDesc = u"""P 1/93""",
    longDesc = 
u"""
P 1/93
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1[CH]C=C=CC=1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 11,
    label = "S4XC12H25",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {19,S} {20,S}
2  C u0 p0 c0 {1,S} {3,S} {21,S} {22,S}
3  C u0 p0 c0 {2,S} {7,S} {23,S} {24,S}
4  C u0 p0 c0 {6,S} {8,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {9,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {11,S} {25,S} {26,S}
8  C u0 p0 c0 {4,S} {12,S} {27,S} {28,S}
9  C u0 p0 c0 {5,S} {12,S} {29,S} {30,S}
10 C u0 p0 c0 {6,S} {31,S} {32,S} {33,S}
11 C u0 p0 c0 {7,S} {34,S} {35,S} {36,S}
12 C u1 p0 c0 {8,S} {9,S} {37,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {1,S}
21 H u0 p0 c0 {2,S}
22 H u0 p0 c0 {2,S}
23 H u0 p0 c0 {3,S}
24 H u0 p0 c0 {3,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {11,S}
37 H u0 p0 c0 {12,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.36787,0.137355,-8.24076e-05,2.36422e-08,-2.47436e-12,-16766.1,48.3522], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[37.9688,0.0538719,-1.82171e-05,2.80775e-09,-1.62108e-13,-31214.5,-165.806], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC[CH]CCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 12,
    label = "SXC8H17",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {8,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {19,S} {20,S} {21,S}
7  C u0 p0 c0 {8,S} {22,S} {23,S} {24,S}
8  C u1 p0 c0 {5,S} {7,S} {25,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.304233,0.0880077,-4.90743e-05,1.21858e-08,-8.87773e-13,-5237.93,36.6583], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[24.9044,0.0366394,-1.2385e-05,1.90835e-09,-1.10161e-13,-14713.5,-101.345], Tmin=(1383,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 13,
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
            NASAPolynomial(coeffs=[1.24237,0.0310722,-5.08669e-05,4.31371e-08,-1.40146e-11,8031.61,13.8743], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.92383,0.00679236,-2.56586e-06,4.49878e-10,-2.99401e-14,7264.63,-7.60177], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""SRI91""",
    longDesc = 
u"""
SRI91
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#CO
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 14,
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
            NASAPolynomial(coeffs=[-0.736766,0.0607201,-3.57593e-05,1.04907e-08,-1.21487e-12,-19893.5,29.5358], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[15.7258,0.0261086,-8.90971e-06,1.38102e-09,-8.00297e-14,-26052,-60.3365], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 15,
    label = "PXC8H17",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {3,S} {13,S} {14,S}
3  C u0 p0 c0 {2,S} {5,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
5  C u0 p0 c0 {3,S} {7,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {8,S} {19,S} {20,S}
7  C u0 p0 c0 {5,S} {21,S} {22,S} {23,S}
8  C u1 p0 c0 {6,S} {24,S} {25,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.772759,0.093255,-5.84447e-05,1.8557e-08,-2.37127e-12,-3926.9,37.6131], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[25.051,0.036948,-1.25765e-05,1.94628e-09,-1.12669e-13,-13330.1,-102.557], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 16,
    label = "S3XC10H21",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {4,S} {17,S} {18,S}
3  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {9,S} {19,S} {20,S}
5  C u0 p0 c0 {6,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {5,S} {10,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {10,S} {23,S} {24,S}
8  C u0 p0 c0 {5,S} {25,S} {26,S} {27,S}
9  C u0 p0 c0 {4,S} {28,S} {29,S} {30,S}
10 C u1 p0 c0 {6,S} {7,S} {31,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.930537,0.113138,-6.64034e-05,1.83221e-08,-1.77128e-12,-10989,42.9335], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[31.4448,0.0452779,-1.53146e-05,2.36072e-09,-1.36312e-13,-22970.3,-133.634], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[CH]CCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 17,
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
            NASAPolynomial(coeffs=[-5.04302,0.0712535,-7.09182e-05,3.86802e-08,-8.78883e-12,-6416.78,48.6171], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.4894,0.0380526,-2.16545e-05,5.92386e-09,-6.27635e-13,-8213.1,7.12481], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""HWZD99""",
    longDesc = 
u"""
HWZD99
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OC1C=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 18,
    label = "NC7H16",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {6,S} {8,S} {9,S}
5  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {18,S} {19,S} {20,S}
7  C u0 p0 c0 {5,S} {21,S} {22,S} {23,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.26836,0.0854356,-5.25347e-05,1.62946e-08,-2.02395e-12,-25658.7,35.3733], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[22.2149,0.0347676,-1.18407e-05,1.83298e-09,-1.0613e-13,-34276,-92.304], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 19,
    label = "SXC7H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {7,S} {19,S} {20,S} {21,S}
7  C u1 p0 c0 {4,S} {6,S} {22,S}
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
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.0379156,0.0756727,-4.07474e-05,9.32679e-09,-4.92361e-13,-2356.05,33.7322], Tmin=(298,'K'), Tmax=(1382,'K')),
            NASAPolynomial(coeffs=[21.6369,0.0323325,-1.09274e-05,1.68357e-09,-9.71774e-14,-10587.4,-85.221], Tmin=(1382,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 20,
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
            NASAPolynomial(coeffs=[3.76268,0.000968872,2.7949e-06,-3.85091e-09,1.68742e-12,46004,1.56253], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.8741,0.00365639,-1.40895e-06,2.6018e-10,-1.87728e-14,46263.6,6.17119], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L S/93""",
    longDesc = 
u"""
L S/93.
[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 21,
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
            NASAPolynomial(coeffs=[3.67359,0.00201095,5.73022e-06,-6.87117e-09,2.54386e-12,16445,1.60456], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.28572,0.0072399,-2.98714e-06,5.95685e-10,-4.67154e-14,16775.6,8.48007], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L11/89""",
    longDesc = 
u"""
L11/89.
[CH3]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 22,
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
            NASAPolynomial(coeffs=[5.14988,-0.013671,4.91801e-05,-4.84743e-08,1.66694e-11,-10246.6,-4.6413], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[0.0748515,0.0133909,-5.73286e-06,1.22293e-09,-1.01815e-13,-9468.34,18.4373], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L 8/88""",
    longDesc = 
u"""
L 8/88.
C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 23,
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
            NASAPolynomial(coeffs=[2.88966,0.01341,-2.8477e-05,2.94791e-08,-1.09332e-11,66839.4,6.22296], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.16781,0.00475222,-1.83787e-06,3.0419e-10,-1.77233e-14,67121.1,6.63589], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L 1/91""",
    longDesc = 
u"""
L 1/91.
[C]#C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 24,
    label = "C2O",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,D} {3,D}
2 C u2 p0 c0 {1,D}
3 O u0 p2 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.86486,0.0119902,-1.83624e-05,1.57697e-08,-5.38975e-12,33749.9,8.88678], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.15127,0.00237267,-7.6136e-07,1.17064e-10,-7.02578e-15,33241.9,-2.21831], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""RUS 79""",
    longDesc = 
u"""
RUS 79.
[C]=C=O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 25,
    label = "CH3CHCH",
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
            NASAPolynomial(coeffs=[0.913729,0.0264323,-1.1759e-05,-2.30357e-09,2.77155e-12,30916.9,19.9893], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.37253,0.0157805,-5.99229e-06,9.30897e-10,-3.6551e-14,29614.8,-3.41865], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""PD5/98""",
    longDesc = 
u"""
PD5/98
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 26,
    label = "S2XC10H21",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {5,S} {17,S} {18,S}
4  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {8,S} {19,S} {20,S}
6  C u0 p0 c0 {4,S} {10,S} {23,S} {24,S}
7  C u0 p0 c0 {9,S} {10,S} {21,S} {22,S}
8  C u0 p0 c0 {5,S} {28,S} {29,S} {30,S}
9  C u0 p0 c0 {7,S} {25,S} {26,S} {27,S}
10 C u1 p0 c0 {6,S} {7,S} {31,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {8,S}
30 H u0 p0 c0 {8,S}
31 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.930537,0.113138,-6.64034e-05,1.83221e-08,-1.77128e-12,-10989,42.9335], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[31.4448,0.0452779,-1.53146e-05,2.36072e-09,-1.36312e-13,-22970.3,-133.634], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]CCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 27,
    label = "C7H14",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {7,D} {19,S}
7  C u0 p0 c0 {6,D} {20,S} {21,S}
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
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.67721,0.0824612,-5.46504e-05,1.87862e-08,-2.65738e-12,-10216.9,38.5068], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[21.0898,0.0310608,-1.05645e-05,1.63406e-09,-9.45598e-14,-18326,-84.4391], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 28,
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
            NASAPolynomial(coeffs=[0.797985,0.0344034,-1.24599e-05,-5.18063e-18,1.9936e-21,-648.928,21.8897], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[-4.72093,0.0391414,-6.52873e-06,-7.68209e-09,2.51473e-12,1753.52,51.719], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""A 8/83""",
    longDesc = 
u"""
A 8/83
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC1CO1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 29,
    label = "H",
    molecule = 
"""
multiplicity 2
1 H u1 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5,7.05333e-13,-1.99592e-15,2.30082e-18,-9.27732e-22,25473.7,-0.446683], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.5,-2.30843e-11,1.61562e-14,-4.73515e-18,4.98197e-22,25473.7,-0.446683], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L 7/88""",
    longDesc = 
u"""
L 7/88.
[H]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 30,
    label = "NC9H20",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {14,S} {15,S}
3  C u0 p0 c0 {2,S} {4,S} {16,S} {17,S}
4  C u0 p0 c0 {3,S} {5,S} {18,S} {19,S}
5  C u0 p0 c0 {4,S} {7,S} {20,S} {21,S}
6  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
7  C u0 p0 c0 {5,S} {9,S} {22,S} {23,S}
8  C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
9  C u0 p0 c0 {7,S} {27,S} {28,S} {29,S}
10 H u0 p0 c0 {6,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.8139,0.110177,-6.93124e-05,2.20958e-08,-2.83356e-12,-31420.8,41.2827], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[28.729,0.0434075,-1.47661e-05,2.28421e-09,-1.32195e-13,-42517.4,-124.429], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 31,
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
            NASAPolynomial(coeffs=[1.56854,0.0346523,6.81681e-06,-2.79951e-08,1.23077e-11,-17130,17.908], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[10.5268,0.0235907,-7.85225e-06,1.14484e-09,-5.98277e-14,-20479.2,-32.1986], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""P11/94""",
    longDesc = 
u"""
P11/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 32,
    label = "CH3CHOCH2",
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
            NASAPolynomial(coeffs=[0.487338,0.0285197,3.00962e-06,-2.26526e-08,1.07067e-11,-12556.4,22.6053], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.69006,0.016021,-5.39718e-06,7.99415e-10,-4.26564e-14,-15420.7,-22.485], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""T 6/92""",
    longDesc = 
u"""
T 6/92
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
CC1CO1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 33,
    label = "S3XC8H17",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
3  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {3,S} {8,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {8,S} {17,S} {18,S}
6  C u0 p0 c0 {3,S} {19,S} {20,S} {21,S}
7  C u0 p0 c0 {2,S} {22,S} {23,S} {24,S}
8  C u1 p0 c0 {4,S} {5,S} {25,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.304233,0.0880077,-4.90743e-05,1.21858e-08,-8.87773e-13,-5237.93,36.6583], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[24.9044,0.0366394,-1.2385e-05,1.90835e-09,-1.10161e-13,-14713.5,-101.345], Tmin=(1383,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[CH]CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 34,
    label = "C4H81",
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
            NASAPolynomial(coeffs=[1.18114,0.0308534,5.08652e-06,-2.46549e-08,1.11102e-11,-1790.4,21.0625], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.05358,0.0343505,-1.58832e-05,3.30897e-09,-2.5361e-13,-2139.72,15.5432], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""T 6/83""",
    longDesc = 
u"""
T 6/83
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 35,
    label = "OH",
    molecule = 
"""
multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.12531,-0.00322545,6.52765e-06,-5.79854e-09,2.06237e-12,3381.54,-0.690433], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.86473,0.0010565,-2.59083e-07,3.05219e-11,-1.33196e-15,3718.86,5.70164], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""S 9/01""",
    longDesc = 
u"""
S 9/01.
[OH]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 36,
    label = "NC10H22",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {7,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {4,S} {17,S} {18,S}
4  C u0 p0 c0 {3,S} {5,S} {19,S} {20,S}
5  C u0 p0 c0 {4,S} {6,S} {21,S} {22,S}
6  C u0 p0 c0 {5,S} {8,S} {23,S} {24,S}
7  C u0 p0 c0 {1,S} {9,S} {11,S} {12,S}
8  C u0 p0 c0 {6,S} {10,S} {25,S} {26,S}
9  C u0 p0 c0 {7,S} {27,S} {28,S} {29,S}
10 C u0 p0 c0 {8,S} {30,S} {31,S} {32,S}
11 H u0 p0 c0 {7,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.08417,0.122535,-7.76816e-05,2.49835e-08,-3.23548e-12,-34302.2,44.226], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[31.9882,0.0477245,-1.62276e-05,2.50963e-09,-1.45216e-13,-46639.3,-140.504], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 37,
    label = "S4XC9H19",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {12,S} {13,S}
2  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {8,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {9,S} {18,S} {19,S}
6  C u0 p0 c0 {2,S} {9,S} {20,S} {21,S}
7  C u0 p0 c0 {3,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {4,S} {25,S} {26,S} {27,S}
9  C u1 p0 c0 {5,S} {6,S} {28,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.576046,0.100512,-5.77755e-05,1.53277e-08,-1.35057e-12,-8122.9,39.5796], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[28.0393,0.041144,-1.3926e-05,2.14746e-09,-1.24022e-13,-18772.8,-116.697], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC[CH]CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 38,
    label = "H2",
    molecule = 
"""
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.34433,0.00798052,-1.94782e-05,2.01572e-08,-7.37612e-12,-917.935,0.68301], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.33728,-4.94025e-05,4.99457e-07,-1.79566e-10,2.00255e-14,-950.159,-3.20502], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""TPIS78""",
    longDesc = 
u"""
TPIS78.
[H][H]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 39,
    label = "pC3H4",
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 40,
    label = "S4XC11H23",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
2  C u0 p0 c0 {1,S} {6,S} {20,S} {21,S}
3  C u0 p0 c0 {5,S} {7,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {8,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {10,S} {22,S} {23,S}
7  C u0 p0 c0 {3,S} {11,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {11,S} {26,S} {27,S}
9  C u0 p0 c0 {5,S} {28,S} {29,S} {30,S}
10 C u0 p0 c0 {6,S} {31,S} {32,S} {33,S}
11 C u1 p0 c0 {7,S} {8,S} {34,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {1,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {2,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.1025,0.125022,-7.40802e-05,2.07819e-08,-2.07856e-12,-13884,45.431], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[34.7028,0.0495634,-1.67589e-05,2.58285e-09,-1.49119e-13,-27089.2,-149.691], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/22/ 7 THERM""",
    longDesc = 
u"""
1/22/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC[CH]CCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 41,
    label = "SAXC5H9",
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
            NASAPolynomial(coeffs=[-1.08167,0.0572223,-4.01417e-05,1.46881e-08,-2.19956e-12,10788.2,31.1232], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[15.1441,0.0193269,-6.34197e-06,9.62601e-10,-5.51208e-14,5169.43,-55.9956], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/29/ 8 G""",
    longDesc = 
u"""
8/29/ 8 G
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[CH]CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 42,
    label = "C6H12",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {6,D} {16,S}
6  C u0 p0 c0 {5,D} {17,S} {18,S}
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
18 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.555477,0.0646057,-3.42004e-05,5.04218e-09,1.14774e-12,-7486.02,31.4621], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[15.0917,0.0290769,-9.5016e-06,1.4765e-09,-8.94808e-14,-12462.9,-51.9476], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/01/7THERG""",
    longDesc = 
u"""
9/01/7THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 43,
    label = "CH3CH2CHO",
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
            NASAPolynomial(coeffs=[2.72557,0.023236,2.97407e-06,-1.66134e-08,7.42501e-12,-24556.7,14.1663], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.26374,0.0199763,-7.61951e-06,1.16871e-09,-4.196e-14,-25886,-5.77865], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""USC/07""",
    longDesc = 
u"""
USC/07
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC=O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 44,
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
            NASAPolynomial(coeffs=[1.49331,0.0209252,4.48679e-06,-1.66891e-08,7.15815e-12,1074.83,16.1453], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.73226,0.0149083,-4.9499e-06,7.21202e-10,-3.7662e-14,-923.57,-13.3133], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""120186""",
    longDesc = 
u"""
120186
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 45,
    label = "S2XC9H19",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {7,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {9,S} {20,S} {21,S}
6  C u0 p0 c0 {8,S} {9,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {25,S} {26,S} {27,S}
8  C u0 p0 c0 {6,S} {22,S} {23,S} {24,S}
9  C u1 p0 c0 {5,S} {6,S} {28,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {7,S}
28 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.576046,0.100512,-5.77755e-05,1.53277e-08,-1.35057e-12,-8122.9,39.5796], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[28.0393,0.041144,-1.3926e-05,2.14746e-09,-1.24022e-13,-18772.8,-116.697], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]CCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 46,
    label = "l-C6H4",
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 47,
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
            NASAPolynomial(coeffs=[1.35111,0.0327411,-4.73827e-05,3.7631e-08,-1.18541e-11,40105.8,15.2059], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.14222,0.00761902,-2.6746e-06,4.24915e-10,-2.51475e-14,38908.7,-12.5848], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 5/97""",
    longDesc = 
u"""
T 5/97.
C#C[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 48,
    label = "S2XC11H23",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
3  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
4  C u0 p0 c0 {3,S} {6,S} {20,S} {21,S}
5  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {9,S} {22,S} {23,S}
7  C u0 p0 c0 {5,S} {11,S} {26,S} {27,S}
8  C u0 p0 c0 {10,S} {11,S} {24,S} {25,S}
9  C u0 p0 c0 {6,S} {31,S} {32,S} {33,S}
10 C u0 p0 c0 {8,S} {28,S} {29,S} {30,S}
11 C u1 p0 c0 {7,S} {8,S} {34,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {7,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {9,S}
32 H u0 p0 c0 {9,S}
33 H u0 p0 c0 {9,S}
34 H u0 p0 c0 {11,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.1025,0.125022,-7.40802e-05,2.07819e-08,-2.07856e-12,-13884,45.431], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[34.7028,0.0495634,-1.67589e-05,2.58285e-09,-1.49119e-13,-27089.2,-149.691], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/22/ 7 THERM""",
    longDesc = 
u"""
1/22/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]CCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 49,
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
            NASAPolynomial(coeffs=[0.928511,0.0264606,6.03324e-06,-2.1915e-08,9.49615e-12,-14057.9,19.2255], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.52442,0.0188983,-6.2921e-06,9.21615e-10,-4.86845e-14,-16564.4,-17.8384], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""P11/94""",
    longDesc = 
u"""
P11/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 50,
    label = "CH2OCH2",
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
            NASAPolynomial(coeffs=[3.75905,-0.00944122,8.03097e-05,-1.00808e-07,4.00399e-11,-7560.81,7.84975], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.48876,0.0120462,-4.33369e-06,7.00283e-10,-4.19491e-14,-9180.43,-7.07996], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""T 6/92""",
    longDesc = 
u"""
T 6/92
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C1CO1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 51,
    label = "C6H10-13",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,D} {12,S}
4  C u0 p0 c0 {3,D} {5,S} {14,S}
5  C u0 p0 c0 {4,S} {6,D} {13,S}
6  C u0 p0 c0 {5,D} {15,S} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.7417,0.0744459,-6.04686e-05,2.60488e-08,-4.55731e-12,4787.62,39.3239], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[17.5264,0.0218448,-7.07365e-06,1.06375e-09,-6.05245e-14,-1718.56,-67.654], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/29/ 8 G""",
    longDesc = 
u"""
8/29/ 8 G
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC=CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 52,
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
            NASAPolynomial(coeffs=[4.79372,-0.00990833,3.7322e-05,-3.79285e-08,1.31773e-11,-14309,0.602813], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[1.76069,0.0092,-4.42259e-06,1.00641e-09,-8.83856e-14,-13995.8,13.6563], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L 8/88""",
    longDesc = 
u"""
L 8/88.
C=O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 53,
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 54,
    label = "C11H22",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
3  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
4  C u0 p0 c0 {3,S} {5,S} {20,S} {21,S}
5  C u0 p0 c0 {4,S} {7,S} {22,S} {23,S}
6  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {24,S} {25,S}
8  C u0 p0 c0 {6,S} {10,S} {26,S} {27,S}
9  C u0 p0 c0 {7,S} {28,S} {29,S} {30,S}
10 C u0 p0 c0 {8,S} {11,D} {31,S}
11 C u0 p0 c0 {10,D} {32,S} {33,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {11,S}
33 H u0 p0 c0 {11,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.69654,0.13165,-8.77957e-05,3.01468e-08,-4.22584e-12,-21752.6,49.988], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[34.1377,0.0483102,-1.64025e-05,2.53434e-09,-1.46557e-13,-34817.1,-148.798], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 55,
    label = "aC3H4",
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 56,
    label = "aC3H5",
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
            NASAPolynomial(coeffs=[1.36318,0.0198138,1.24971e-05,-3.33556e-08,1.58466e-11,19245.6,17.1732], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.50079,0.0143247,-5.67816e-06,1.10808e-09,-9.03639e-14,17482.4,-11.2431], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""PD5/98""",
    longDesc = 
u"""
PD5/98
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 57,
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
            NASAPolynomial(coeffs=[4.27611,-0.000542822,1.67336e-05,-2.15771e-08,8.62454e-12,-17702.6,3.43505], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.165,0.00490832,-1.90139e-06,3.71186e-10,-2.87908e-14,-17861.8,2.91616], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L 7/88""",
    longDesc = 
u"""
L 7/88.
OO
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 58,
    label = "C6H5CH2OH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u0 p0 c0 {2,S} {5,D} {11,S}
4  C u0 p0 c0 {2,D} {7,S} {15,S}
5  C u0 p0 c0 {3,D} {6,S} {12,S}
6  C u0 p0 c0 {5,S} {7,D} {13,S}
7  C u0 p0 c0 {4,S} {6,D} {14,S}
8  O u0 p2 c0 {1,S} {16,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.0642,0.0227751,9.59721e-05,-1.50851e-07,6.41758e-11,-14285,18.1483], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[15.2812,0.0272085,-9.85847e-06,1.60122e-09,-9.62781e-14,-19700.5,-59.4187], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 7/87""",
    longDesc = 
u"""
L 7/87.
OCC1C=CC=CC=1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 59,
    label = "CH3CHCHCO",
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
            NASAPolynomial(coeffs=[5.30535,0.0157494,2.16239e-05,-3.66078e-08,1.49325e-11,5758.86,4.20435], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.76082,0.0200318,-8.0631e-06,1.33614e-09,-6.23084e-14,4570.83,-11.0956], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""USC/07""",
    longDesc = 
u"""
USC/07
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 60,
    label = "pC4H9",
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
            NASAPolynomial(coeffs=[1.2087,0.0382975,-7.26605e-06,-1.54285e-08,8.68594e-12,7322.1,22.1693], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.68224,0.0236911,-7.59489e-06,6.64271e-10,5.48451e-14,4964.41,-17.8917], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""USC/07""",
    longDesc = 
u"""
USC/07
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 61,
    label = "SXC11H23",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
3  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
4  C u0 p0 c0 {3,S} {5,S} {20,S} {21,S}
5  C u0 p0 c0 {4,S} {7,S} {22,S} {23,S}
6  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
7  C u0 p0 c0 {5,S} {9,S} {24,S} {25,S}
8  C u0 p0 c0 {6,S} {11,S} {26,S} {27,S}
9  C u0 p0 c0 {7,S} {28,S} {29,S} {30,S}
10 C u0 p0 c0 {11,S} {31,S} {32,S} {33,S}
11 C u1 p0 c0 {8,S} {10,S} {34,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.1025,0.125022,-7.40802e-05,2.07819e-08,-2.07856e-12,-13884,45.431], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[34.7028,0.0495634,-1.67589e-05,2.58285e-09,-1.49119e-13,-27089.2,-149.691], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/22/ 7 THERM""",
    longDesc = 
u"""
1/22/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 62,
    label = "C4H82",
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
            NASAPolynomial(coeffs=[1.25943,0.0278084,8.70139e-06,-2.44022e-08,9.89777e-12,-2964.77,20.5011], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[0.827977,0.0358645,-1.66345e-05,3.47328e-09,-2.66574e-13,-3052.1,21.3425], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""T 6/83""",
    longDesc = 
u"""
T 6/83
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 63,
    label = "SXC6H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {6,S} {16,S} {17,S} {18,S}
6  C u1 p0 c0 {3,S} {5,S} {19,S}
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
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.22956,0.0633327,-3.24135e-05,6.46388e-09,-9.6142e-14,525.639,30.8006], Tmin=(298,'K'), Tmax=(1380,'K')),
            NASAPolynomial(coeffs=[18.3687,0.0280268,-9.47032e-06,1.45889e-09,-8.42002e-14,-6460.94,-69.0934], Tmin=(1380,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 64,
    label = "C5H10",
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
            NASAPolynomial(coeffs=[-0.500287,0.0537811,-2.94505e-05,5.61977e-09,3.87026e-13,-4574.86,29.4517], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[12.2563,0.0244603,-8.0484e-06,1.25615e-09,-7.63205e-14,-8621.53,-38.4776], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERG""",
    longDesc = 
u"""
1/ 2/ 7 THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 65,
    label = "SXC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
3  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {3,S} {14,S}
5  C u0 p0 c0 {2,S} {6,D} {15,S}
6  C u0 p0 c0 {5,D} {16,S} {17,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.737326,0.0696677,-5.65392e-05,2.81807e-08,-6.59685e-12,15951.9,32.92], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[16.5055,0.0231158,-6.52726e-06,9.14692e-10,-5.16951e-14,10821.8,-57.1912], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERGAS G""",
    longDesc = 
u"""
THERGAS G
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCC[CH]C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 66,
    label = "cC3H4",
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 67,
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
            NASAPolynomial(coeffs=[2.25172,0.017655,-2.37291e-05,1.72758e-08,-5.06648e-12,20059.4,12.4904], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.62821,0.00408534,-1.59345e-06,2.86261e-10,-1.94078e-14,19327.2,-3.93026], Tmin=(1000,'K'), Tmax=(4000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (4000,'K'),
    ),
    shortDesc = u"""SRIC91""",
    longDesc = 
u"""
SRIC91
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#C[O]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 68,
    label = "PXC5H11",
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
            NASAPolynomial(coeffs=[0.0524384,0.0560797,-3.31546e-05,9.77534e-09,-1.1401e-12,4716.11,28.7239], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[15.2977,0.0239735,-8.18393e-06,1.26883e-09,-7.35409e-14,-980.712,-54.4829], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 69,
    label = "C",
    molecule = 
"""
1 C u0 p2 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.55424,-0.000321538,7.33792e-07,-7.32235e-10,2.66521e-13,85443.9,4.53131], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.49267,4.79889e-05,-7.24335e-08,3.74291e-11,-4.87278e-15,85451.3,4.8015], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L11/88""",
    longDesc = 
u"""
L11/88.
[C]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 70,
    label = "CH2OCH",
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
            NASAPolynomial(coeffs=[-0.383961,0.023879,-1.24676e-05,-1.76864e-09,2.81424e-12,18836.2,25.7417], Tmin=(298,'K'), Tmax=(500,'K')),
            NASAPolynomial(coeffs=[4.49941,0.0115526,-4.81441e-06,8.92349e-10,-5.68706e-14,17474,0.339255], Tmin=(500,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""A12/04""",
    longDesc = 
u"""
A12/04
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
[CH]1CO1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 71,
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
            NASAPolynomial(coeffs=[1.21087,0.0352059,-1.09391e-05,-1.17206e-08,7.61749e-12,2266.57,20.6135], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.30106,0.0199453,-8.29038e-06,1.51008e-09,-9.15812e-14,157.884,-16.9106], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""USC/07""",
    longDesc = 
u"""
USC/07
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CC=O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 72,
    label = "S3XC7H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {8,S} {9,S}
2  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {2,S} {19,S} {20,S} {21,S}
7  C u1 p0 c0 {3,S} {4,S} {22,S}
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
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.0379156,0.0756727,-4.07474e-05,9.32679e-09,-4.92361e-13,-2356.05,33.7322], Tmin=(298,'K'), Tmax=(1382,'K')),
            NASAPolynomial(coeffs=[21.6369,0.0323325,-1.09274e-05,1.68357e-09,-9.71774e-14,-10587.4,-85.221], Tmin=(1382,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[CH]CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 73,
    label = "O",
    molecule = 
"""
multiplicity 3
1 O u2 p2 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.16827,-0.00327932,6.64306e-06,-6.12807e-09,2.11266e-12,29122.3,2.05193], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.56942,-8.59741e-05,4.19485e-08,-1.00178e-11,1.22834e-15,29217.6,4.78434], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L 1/90""",
    longDesc = 
u"""
L 1/90.
[O]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 74,
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
            NASAPolynomial(coeffs=[0.58857,0.0365467,-3.4107e-05,1.66526e-08,-3.00646e-12,33359.5,20.6579], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.25396,0.0139141,-5.29322e-06,8.34805e-10,-3.51979e-14,31766,-12.6295], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""USC/07""",
    longDesc = 
u"""
USC/07
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#CC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 75,
    label = "PXC11H21",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
3  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
4  C u0 p0 c0 {3,S} {6,S} {20,S} {21,S}
5  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
6  C u0 p0 c0 {4,S} {7,S} {22,S} {23,S}
7  C u0 p0 c0 {6,S} {9,S} {26,S} {27,S}
8  C u0 p0 c0 {5,S} {10,S} {24,S} {25,S}
9  C u0 p0 c0 {7,S} {11,D} {28,S}
10 C u1 p0 c0 {8,S} {29,S} {30,S}
11 C u0 p0 c0 {9,D} {31,S} {32,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {7,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {11,S}
32 H u0 p0 c0 {11,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.77277,0.128898,-8.57717e-05,2.92169e-08,-4.05817e-12,-4764.1,49.2434], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[33.8968,0.0466347,-1.59682e-05,2.4812e-09,-1.44046e-13,-17863.8,-148.943], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCCCCCC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 76,
    label = "PXC11H23",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
3  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
4  C u0 p0 c0 {3,S} {5,S} {20,S} {21,S}
5  C u0 p0 c0 {4,S} {6,S} {22,S} {23,S}
6  C u0 p0 c0 {5,S} {8,S} {24,S} {25,S}
7  C u0 p0 c0 {1,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {6,S} {10,S} {26,S} {27,S}
9  C u0 p0 c0 {7,S} {11,S} {28,S} {29,S}
10 C u0 p0 c0 {8,S} {30,S} {31,S} {32,S}
11 C u1 p0 c0 {9,S} {33,S} {34,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {11,S}
34 H u0 p0 c0 {11,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.5823,0.130324,-8.35408e-05,2.72126e-08,-3.5753e-12,-12571.3,46.4373], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[34.8306,0.0498968,-1.69602e-05,2.62239e-09,-1.51722e-13,-25696.4,-150.794], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/22/ 7 THERM""",
    longDesc = 
u"""
1/22/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 77,
    label = "C8H16",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {7,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {19,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {8,D} {22,S}
8  C u0 p0 c0 {7,D} {23,S} {24,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.89227,0.0946066,-6.27386e-05,2.15158e-08,-3.02719e-12,-13107.5,41.1879], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[24.354,0.0353666,-1.20208e-05,1.85855e-09,-1.07522e-13,-22448.6,-100.538], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 78,
    label = "iC3H7",
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
            NASAPolynomial(coeffs=[1.44492,0.0209991,7.70362e-06,-1.84763e-08,7.1283e-12,9422.37,20.1163], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.51927,0.0172201,-5.73642e-06,8.41307e-10,-4.45659e-14,7322.72,-9.08302], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""P11/94""",
    longDesc = 
u"""
P11/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 79,
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
            NASAPolynomial(coeffs=[4.29142,-0.00550154,5.99438e-05,-7.08466e-08,2.68686e-11,-11522.2,2.66682], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[1.07188,0.0216853,-1.00256e-05,2.21412e-09,-1.90003e-13,-11426.4,15.1156], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L 8/88""",
    longDesc = 
u"""
L 8/88.
CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 80,
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 81,
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
            NASAPolynomial(coeffs=[3.9592,-0.00757052,5.7099e-05,-6.91589e-08,2.69884e-11,5089.78,4.09733], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.03611,0.0146454,-6.71078e-06,1.47223e-09,-1.25706e-13,4939.89,10.3054], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L 1/91""",
    longDesc = 
u"""
L 1/91.
C=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 82,
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 83,
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 84,
    label = "C6H4O2",
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
            NASAPolynomial(coeffs=[-0.95193,0.0578424,-3.82144e-05,4.63127e-09,3.62967e-12,-17611,29.2395], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[11.7308,0.023615,-1.02346e-05,1.95322e-09,-1.2746e-13,-21085.8,-36.3005], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""PUML96""",
    longDesc = 
u"""
PUML96
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=C1C=CC(=O)C=C1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 85,
    label = "nC3H7",
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
            NASAPolynomial(coeffs=[1.04912,0.026009,2.35425e-06,-1.95951e-08,9.37202e-12,10312.3,21.136], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.70975,0.0160315,-5.27202e-06,7.58884e-10,-3.88627e-14,7976.22,-15.5153], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""P11/94""",
    longDesc = 
u"""
P11/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 86,
    label = "SAXC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {2,S} {5,S} {14,S}
5  C u0 p0 c0 {4,S} {6,D} {15,S}
6  C u0 p0 c0 {5,D} {16,S} {17,S}
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
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.102852,0.0623649,-3.81471e-05,1.15531e-08,-1.44216e-12,8354.81,29.0126], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[16.0125,0.0243247,-7.05642e-06,1.00566e-09,-5.73513e-14,3255.88,-56.644], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[CH]CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 87,
    label = "S3XC11H23",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {16,S} {17,S}
2  C u0 p0 c0 {1,S} {3,S} {18,S} {19,S}
3  C u0 p0 c0 {2,S} {5,S} {20,S} {21,S}
4  C u0 p0 c0 {1,S} {8,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {10,S} {22,S} {23,S}
6  C u0 p0 c0 {7,S} {9,S} {12,S} {13,S}
7  C u0 p0 c0 {6,S} {11,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {11,S} {26,S} {27,S}
9  C u0 p0 c0 {6,S} {28,S} {29,S} {30,S}
10 C u0 p0 c0 {5,S} {31,S} {32,S} {33,S}
11 C u1 p0 c0 {7,S} {8,S} {34,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.1025,0.125022,-7.40802e-05,2.07819e-08,-2.07856e-12,-13884,45.431], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[34.7028,0.0495634,-1.67589e-05,2.58285e-09,-1.49119e-13,-27089.2,-149.691], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/22/ 7 THERM""",
    longDesc = 
u"""
1/22/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[CH]CCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 88,
    label = "NC8H18",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {3,S} {13,S} {14,S}
3  C u0 p0 c0 {2,S} {4,S} {15,S} {16,S}
4  C u0 p0 c0 {3,S} {6,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {7,S} {9,S} {10,S}
6  C u0 p0 c0 {4,S} {8,S} {19,S} {20,S}
7  C u0 p0 c0 {5,S} {21,S} {22,S} {23,S}
8  C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.54218,0.0978112,-6.09318e-05,1.92006e-08,-2.42996e-12,-28539.6,38.3328], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[25.471,0.0390887,-1.33039e-05,2.05868e-09,-1.19167e-13,-38396.3,-108.361], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 89,
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
            NASAPolynomial(coeffs=[0.451,0.0674752,-0.000118099,1.03676e-07,-3.4851e-11,82173.1,17.7041], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[12.8939,0.00791451,-2.40272e-06,2.43401e-10,3.13832e-15,79832.4,-40.772], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""D11/99""",
    longDesc = 
u"""
D11/99
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#CC#CC#C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 90,
    label = "S2XC7H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
4  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {19,S} {20,S} {21,S}
6  C u0 p0 c0 {4,S} {16,S} {17,S} {18,S}
7  C u1 p0 c0 {3,S} {4,S} {22,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.0379156,0.0756727,-4.07474e-05,9.32679e-09,-4.92361e-13,-2356.05,33.7322], Tmin=(298,'K'), Tmax=(1382,'K')),
            NASAPolynomial(coeffs=[21.6369,0.0323325,-1.09274e-05,1.68357e-09,-9.71774e-14,-10587.4,-85.221], Tmin=(1382,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 91,
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
            NASAPolynomial(coeffs=[-4.84377,0.0584276,-2.94859e-05,-6.93904e-09,8.21253e-12,9181.78,43.8898], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[9.13812,0.0238544,-8.81277e-06,1.2099e-09,-1.82215e-14,5204.35,-29.1157], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""D11/99""",
    longDesc = 
u"""
D11/99
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C1=CC=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 92,
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
            NASAPolynomial(coeffs=[-3.69315,0.052179,-2.55584e-05,-7.06611e-09,7.5834e-12,39779.6,41.3325], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.59731,0.0222416,-8.72e-06,1.37888e-09,-5.31461e-14,36261,-22.9546], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""D11/99""",
    longDesc = 
u"""
D11/99
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[C]1=CC=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 93,
    label = "S3XC9H19",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {8,S} {16,S} {17,S}
4  C u0 p0 c0 {5,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {9,S} {18,S} {19,S}
6  C u0 p0 c0 {2,S} {9,S} {20,S} {21,S}
7  C u0 p0 c0 {4,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {3,S} {25,S} {26,S} {27,S}
9  C u1 p0 c0 {5,S} {6,S} {28,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.576046,0.100512,-5.77755e-05,1.53277e-08,-1.35057e-12,-8122.9,39.5796], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[28.0393,0.041144,-1.3926e-05,2.14746e-09,-1.24022e-13,-18772.8,-116.697], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[CH]CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 94,
    label = "c-C4H5",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {4,S} {7,S}
3 C u0 p0 c0 {1,S} {4,D} {8,S}
4 C u0 p0 c0 {2,S} {3,D} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.63976,0.0415492,-2.1921e-05,-4.6559e-09,6.13489e-12,35373.8,35.7018], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.74672,0.017283,-6.51686e-06,9.89176e-10,-3.46049e-14,32808.4,-12.9129], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""PUPM3""",
    longDesc = 
u"""
PUPM3
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]1C=CC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 95,
    label = "PXC8H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {6,S} {17,S} {18,S}
5  C u0 p0 c0 {2,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {8,D} {19,S}
7  C u1 p0 c0 {5,S} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {22,S} {23,S}
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
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.90099,0.0915068,-6.01588e-05,2.02338e-08,-2.78235e-12,3872.14,40.1405], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[24.1485,0.0336466,-1.15693e-05,1.80262e-09,-1.04847e-13,-5516.31,-100.895], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCCC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 96,
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
            NASAPolynomial(coeffs=[3.18119,0.0298408,-3.28324e-05,2.06318e-08,-5.42006e-12,24125.6,9.42101], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.42922,0.0105027,-4.20668e-06,7.11849e-10,-3.57966e-14,22907.8,-16.512], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""USC/07""",
    longDesc = 
u"""
USC/07
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C=C=C=O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 97,
    label = "PXC7H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,D} {16,S}
6  C u1 p0 c0 {4,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.66946,0.0793202,-5.20232e-05,1.74804e-08,-2.40913e-12,6759.27,37.376], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[20.9278,0.0292841,-1.009e-05,1.57426e-09,-9.16506e-14,-1412.96,-85.0447], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 98,
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 99,
    label = "PXC7H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {18,S} {19,S} {20,S}
7  C u1 p0 c0 {5,S} {21,S} {22,S}
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
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.49957,0.0808826,-5.00533e-05,1.56549e-08,-1.96616e-12,-1045.9,34.6564], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[21.7941,0.032628,-1.11138e-05,1.72067e-09,-9.96367e-14,-9209.38,-86.4954], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 100,
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
            NASAPolynomial(coeffs=[2.13733,0.0264862,-9.05687e-06,-5.53864e-19,2.12819e-22,15710.9,13.5294], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[9.03381,0.00821245,7.1754e-06,-5.88343e-09,1.03439e-12,14335.1,-20.9858], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""A 8/83""",
    longDesc = 
u"""
A 8/83
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC#CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 101,
    label = "S3XC12H25",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {17,S} {18,S}
2  C u0 p0 c0 {1,S} {3,S} {19,S} {20,S}
3  C u0 p0 c0 {2,S} {4,S} {21,S} {22,S}
4  C u0 p0 c0 {3,S} {6,S} {23,S} {24,S}
5  C u0 p0 c0 {1,S} {9,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {11,S} {25,S} {26,S}
7  C u0 p0 c0 {8,S} {10,S} {13,S} {14,S}
8  C u0 p0 c0 {7,S} {12,S} {27,S} {28,S}
9  C u0 p0 c0 {5,S} {12,S} {29,S} {30,S}
10 C u0 p0 c0 {7,S} {31,S} {32,S} {33,S}
11 C u0 p0 c0 {6,S} {34,S} {35,S} {36,S}
12 C u1 p0 c0 {8,S} {9,S} {37,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {3,S}
23 H u0 p0 c0 {4,S}
24 H u0 p0 c0 {4,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {6,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {11,S}
37 H u0 p0 c0 {12,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.36787,0.137355,-8.24076e-05,2.36422e-08,-2.47436e-12,-16766.1,48.3522], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[37.9688,0.0538719,-1.82171e-05,2.80775e-09,-1.62108e-13,-31214.5,-165.806], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[CH]CCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 102,
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 103,
    label = "C5H5O(2,4)",
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
            NASAPolynomial(coeffs=[-3.07776,0.0525817,-2.88565e-05,-3.38855e-09,6.33614e-12,25510.5,39.5915], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.54053,0.0229895,-9.54376e-06,1.70616e-09,-9.74594e-14,22263.7,-20.8188], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""D 9/97""",
    longDesc = 
u"""
D 9/97
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]C1C=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 104,
    label = "SXC5H9",
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
            NASAPolynomial(coeffs=[-0.479207,0.0574354,-4.84093e-05,2.5387e-08,-6.16554e-12,18841.4,30.043], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[13.6419,0.018558,-5.11222e-06,7.02796e-10,-3.91545e-14,14672.4,-43.577], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERGAS G""",
    longDesc = 
u"""
THERGAS G
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC[CH]C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 105,
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 106,
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 107,
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
            NASAPolynomial(coeffs=[3.21247,0.00151479,2.59209e-05,-3.57658e-08,1.47151e-11,34859.8,8.51054], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.01672,0.0103302,-4.68082e-06,1.01763e-09,-8.62607e-14,34612.9,7.78732], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L 2/92""",
    longDesc = 
u"""
L 2/92.
[CH]=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 108,
    label = "SXC10H21",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {4,S} {17,S} {18,S}
4  C u0 p0 c0 {3,S} {6,S} {19,S} {20,S}
5  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {21,S} {22,S}
7  C u0 p0 c0 {5,S} {10,S} {23,S} {24,S}
8  C u0 p0 c0 {6,S} {25,S} {26,S} {27,S}
9  C u0 p0 c0 {10,S} {28,S} {29,S} {30,S}
10 C u1 p0 c0 {7,S} {9,S} {31,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.930537,0.113138,-6.64034e-05,1.83221e-08,-1.77128e-12,-10989,42.9335], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[31.4448,0.0452779,-1.53146e-05,2.36072e-09,-1.36312e-13,-22970.3,-133.634], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 109,
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
            NASAPolynomial(coeffs=[4.22119,-0.00324393,1.37799e-05,-1.33144e-08,4.33769e-12,3839.56,3.39437], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.77217,0.00495696,-2.48446e-06,5.89162e-10,-5.33509e-14,4011.92,9.79834], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L12/89""",
    longDesc = 
u"""
L12/89.
[CH]=O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 110,
    label = "NC11H24",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {8,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
3  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
4  C u0 p0 c0 {3,S} {5,S} {20,S} {21,S}
5  C u0 p0 c0 {4,S} {6,S} {22,S} {23,S}
6  C u0 p0 c0 {5,S} {7,S} {24,S} {25,S}
7  C u0 p0 c0 {6,S} {9,S} {26,S} {27,S}
8  C u0 p0 c0 {1,S} {10,S} {12,S} {13,S}
9  C u0 p0 c0 {7,S} {11,S} {28,S} {29,S}
10 C u0 p0 c0 {8,S} {30,S} {31,S} {32,S}
11 C u0 p0 c0 {9,S} {33,S} {34,S} {35,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {7,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {11,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.35338,0.134888,-8.60424e-05,2.78658e-08,-3.6362e-12,-37183.8,47.1645], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[35.2485,0.0520402,-1.76887e-05,2.73497e-09,-1.58232e-13,-50761.6,-156.585], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/22/ 7 THERM""",
    longDesc = 
u"""
1/22/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 111,
    label = "PXC5H9",
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
            NASAPolynomial(coeffs=[-1.22721,0.0549419,-3.56284e-05,1.18245e-08,-1.61716e-12,12539.2,31.9626], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[14.5547,0.0204332,-7.07391e-06,1.1072e-09,-6.46021e-14,6768.13,-53.7244], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 112,
    label = "CH2CHCO",
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
            NASAPolynomial(coeffs=[3.21169,0.0118422,1.67463e-05,-3.06947e-08,1.33049e-11,7128.16,10.0882], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.95842,0.0107193,-3.85218e-06,6.22009e-10,-3.72402e-14,5648.26,-11.4746], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T05/99""",
    longDesc = 
u"""
T05/99.
C=C[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 113,
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
            NASAPolynomial(coeffs=[5.55579,-0.00283654,7.05689e-05,-8.78105e-08,3.40283e-11,-28113.3,2.32266], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.2976,0.0175662,-6.31705e-06,1.02031e-09,-6.1094e-14,-29817.7,-12.757], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 5/92""",
    longDesc = 
u"""
T 5/92.
CC(C)=O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 114,
    label = "CH",
    molecule = 
"""
multiplicity 4
1 C u3 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.48982,0.000323836,-1.68899e-06,3.16217e-09,-1.40609e-12,70797.3,2.08401], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.87846,0.000970914,1.44446e-07,-1.30688e-10,1.76079e-14,71012.4,5.48498], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""TPIS79""",
    longDesc = 
u"""
TPIS79.
[CH]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 115,
    label = "CO",
    molecule = 
"""
1 C u0 p1 c-1 {2,T}
2 O u0 p1 c+1 {1,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.57953,-0.000610354,1.01681e-06,9.07006e-10,-9.04424e-13,-14344.1,3.50841], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.71519,0.00206253,-9.98826e-07,2.30053e-10,-2.03648e-14,-14151.9,7.81869], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""TPIS79""",
    longDesc = 
u"""
TPIS79.
[C-]#[O+]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 116,
    label = "cC5H8",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {5,D} {12,S}
5  C u0 p0 c0 {2,S} {4,D} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.68981,0.00209545,0.000113037,-1.54081e-07,6.27637e-11,2313.97,15.2941], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.72448,0.0283223,-1.15452e-05,2.15408e-09,-1.50542e-13,-782.616,-19.7697], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T03/97""",
    longDesc = 
u"""
T03/97.
C1=CCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 117,
    label = "CH3CCH2",
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
            NASAPolynomial(coeffs=[1.73292,0.0223946,-5.14906e-06,-6.75965e-09,3.82532e-12,29040.5,16.5689], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.42555,0.0155111,-5.66784e-06,7.92244e-10,-1.6878e-14,27843,-3.35272], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""PD5/98""",
    longDesc = 
u"""
PD5/98
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=[C]C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 118,
    label = "o-C6H4",
    molecule = 
"""
multiplicity 3
1  C u0 p0 c0 {2,S} {3,D} {8,S}
2  C u1 p0 c0 {1,S} {4,S} {7,S}
3  C u0 p0 c0 {1,D} {6,S} {10,S}
4  C u1 p0 c0 {2,S} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,T}
6  C u0 p0 c0 {3,S} {5,T}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.84542,0.0583916,-4.86448e-05,1.67703e-08,-7.85807e-13,52592.5,40.5871], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.8433,0.0203015,-8.86743e-06,1.72643e-09,-1.1786e-13,49317.1,-24.0143], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""D11/99""",
    longDesc = 
u"""
D11/99
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C1#CC=C[CH][CH]1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 119,
    label = "C9H18",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {14,S} {15,S}
3  C u0 p0 c0 {2,S} {4,S} {16,S} {17,S}
4  C u0 p0 c0 {3,S} {6,S} {18,S} {19,S}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {8,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {6,S} {9,D} {25,S}
9  C u0 p0 c0 {8,D} {26,S} {27,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.16108,0.106958,-7.10973e-05,2.43971e-08,-3.42772e-12,-15989.1,44.1245], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[27.6142,0.0396825,-1.34819e-05,2.0839e-09,-1.20539e-13,-26570.9,-116.619], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 120,
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
            NASAPolynomial(coeffs=[1.27135,0.0262311,-9.29123e-06,-4.78373e-09,3.34805e-12,-9335.73,19.4981], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.81119,0.0171143,-7.48342e-06,1.42522e-09,-9.17468e-14,-10784.1,-4.8588], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""USC/07""",
    longDesc = 
u"""
USC/07
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC=O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 121,
    label = "iC4H10",
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
            NASAPolynomial(coeffs=[0.541095,0.0378603,5.54598e-06,-3.05001e-08,1.40334e-11,-17977.6,21.1509], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[10.8462,0.0233384,-7.7834e-06,1.13938e-09,-5.99183e-14,-21669.9,-35.8706], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""P11/94""",
    longDesc = 
u"""
P11/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC(C)C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 122,
    label = "N2",
    molecule = 
"""
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.29868,0.00140824,-3.96322e-06,5.64152e-09,-2.44485e-12,-1020.9,3.95037], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.92664,0.00148798,-5.68476e-07,1.0097e-10,-6.75335e-15,-922.798,5.98053], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""121286""",
    longDesc = 
u"""
121286
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
N#N
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 123,
    label = "iC4H3",
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
            NASAPolynomial(coeffs=[3.72215,0.0259575,-2.63563e-05,1.55089e-08,-3.80406e-12,58837.1,7.56372], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.65385,0.0112041,-4.64013e-06,8.67866e-10,-5.74306e-14,57954.4,-11.7565], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""USC/07""",
    longDesc = 
u"""
USC/07
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C#C[C]=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 124,
    label = "CO2",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,D}
2 O u0 p2 c0 {1,D}
3 O u0 p2 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.35677,0.0089846,-7.12356e-06,2.45919e-09,-1.437e-13,-48372,9.90105], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.85746,0.00441437,-2.21481e-06,5.2349e-10,-4.72084e-14,-48759.2,2.27164], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L 7/88""",
    longDesc = 
u"""
L 7/88.
O=C=O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 125,
    label = "iC4H5",
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
            NASAPolynomial(coeffs=[0.113081,0.0409506,-3.54136e-05,1.5531e-08,-2.33551e-12,36383.4,23.6925], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.9646,0.0182743,-7.81337e-06,1.52922e-09,-1.09205e-13,34725.1,-10.6493], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""USC/07""",
    longDesc = 
u"""
USC/07
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=[C]C=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 126,
    label = "iC4H7",
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
            NASAPolynomial(coeffs=[-1.03759,0.0455667,-3.04762e-05,7.11026e-09,9.96857e-13,14896.5,29.8637], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.14859,0.0221897,-8.44002e-06,1.31334e-09,-5.16179e-14,12712.3,-12.1312], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""USC/07""",
    longDesc = 
u"""
USC/07
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 127,
    label = "iC4H9",
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
            NASAPolynomial(coeffs=[0.975279,0.0416138,-1.44673e-05,-9.38524e-09,6.87974e-12,6668.83,21.2776], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.49817,0.0246895,-8.64876e-06,1.07793e-09,-6.43406e-16,4428.82,-18.4414], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""USC/07""",
    longDesc = 
u"""
USC/07
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 128,
    label = "iC4H8",
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
            NASAPolynomial(coeffs=[2.64714,0.025903,8.19854e-06,-2.21933e-08,8.89586e-12,-4037.31,12.6764], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.46095,0.0296115,-1.30771e-05,2.65719e-09,-2.01347e-13,-5006.68,1.06715], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""T 6/83""",
    longDesc = 
u"""
T 6/83
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 129,
    label = "S4XC10H21",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {17,S} {18,S}
2  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {9,S} {19,S} {20,S}
6  C u0 p0 c0 {2,S} {10,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {10,S} {23,S} {24,S}
8  C u0 p0 c0 {4,S} {25,S} {26,S} {27,S}
9  C u0 p0 c0 {5,S} {28,S} {29,S} {30,S}
10 C u1 p0 c0 {6,S} {7,S} {31,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.930537,0.113138,-6.64034e-05,1.83221e-08,-1.77128e-12,-10989,42.9335], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[31.4448,0.0452779,-1.53146e-05,2.36072e-09,-1.36312e-13,-22970.3,-133.634], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC[CH]CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 130,
    label = "H2O",
    molecule = 
"""
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.19864,-0.00203643,6.5204e-06,-5.48797e-09,1.77198e-12,-30293.7,-0.849032], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.03399,0.00217692,-1.64073e-07,-9.7042e-11,1.68201e-14,-30004.3,4.96677], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L 8/89""",
    longDesc = 
u"""
L 8/89.
O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 131,
    label = "PXC9H17",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
5  C u0 p0 c0 {4,S} {7,S} {20,S} {21,S}
6  C u0 p0 c0 {3,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {5,S} {9,D} {22,S}
8  C u1 p0 c0 {6,S} {23,S} {24,S}
9  C u0 p0 c0 {7,D} {25,S} {26,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.20061,0.103997,-6.87317e-05,2.32482e-08,-3.21153e-12,995.159,42.3691], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[27.3846,0.0379931,-1.30426e-05,2.02999e-09,-1.17985e-13,-9626.53,-117.686], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 3/ 7 THERM""",
    longDesc = 
u"""
1/ 3/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCCCC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 132,
    label = "cC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {6,S} {15,S} {16,S}
6  C u1 p0 c0 {4,S} {5,S} {17,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
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
            NASAPolynomial(coeffs=[-7.87845,0.0829126,-5.41053e-05,1.87956e-08,-3.25634e-12,8102.21,60.6987], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[12.2446,0.0344416,-1.19965e-05,1.93711e-09,-1.19843e-13,1765.76,-46.0541], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 133,
    label = "cC6H10",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,D} {15,S}
6  C u0 p0 c0 {3,S} {5,D} {16,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.76213,0.0722665,-3.4207e-05,-6.55031e-09,7.77014e-12,-1870.35,50.2213], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[10.0457,0.0342604,-1.28037e-05,2.17489e-09,-1.39539e-13,-6408.85,-32.6015], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERGA""",
    longDesc = 
u"""
THERGA
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C1=CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 134,
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
            NASAPolynomial(coeffs=[4.30647,-0.00418659,4.97143e-05,-5.99127e-08,2.30509e-11,12841.6,4.70721], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[1.95466,0.0173973,-7.98207e-06,1.75218e-09,-1.49642e-13,12857.5,13.4624], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L12/92""",
    longDesc = 
u"""
L12/92.
C[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 135,
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 136,
    label = "S2XC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u1 p0 c0 {1,S} {2,S} {14,S}
5  C u0 p0 c0 {2,S} {6,D} {15,S}
6  C u0 p0 c0 {5,D} {16,S} {17,S}
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
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.737326,0.0696677,-5.65392e-05,2.81807e-08,-6.59685e-12,15951.9,32.92], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[16.5055,0.0231158,-6.52726e-06,9.14692e-10,-5.16951e-14,10821.8,-57.1912], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERGAS G""",
    longDesc = 
u"""
THERGAS G
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC[CH]CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 137,
    label = "C5H5O(1,3)",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,D}
3  C u1 p0 c0 {1,S} {5,S} {9,S}
4  C u0 p0 c0 {2,S} {5,D} {10,S}
5  C u0 p0 c0 {3,S} {4,D} {11,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  O u0 p2 c0 {2,D}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.9567,0.0558519,-3.72416e-05,4.16244e-09,3.9272e-12,4857.32,38.6767], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[9.24314,0.0222013,-9.31059e-06,1.71552e-09,-1.0614e-13,1590.84,-24.0877], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""DU0997""",
    longDesc = 
u"""
DU0997
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=C1C=C[CH]C1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 138,
    label = "S2XC6H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
5  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
6  C u1 p0 c0 {2,S} {3,S} {19,S}
7  H u0 p0 c0 {1,S}
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
19 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.22956,0.0633327,-3.24135e-05,6.46388e-09,-9.6142e-14,525.639,30.8006], Tmin=(298,'K'), Tmax=(1380,'K')),
            NASAPolynomial(coeffs=[18.3687,0.0280268,-9.47032e-06,1.45889e-09,-8.42002e-14,-6460.94,-69.0934], Tmin=(1380,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 139,
    label = "HOC6H4CH3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {5,D}
3  C u0 p0 c0 {4,D} {6,S} {8,S}
4  C u0 p0 c0 {2,S} {3,D} {15,S}
5  C u0 p0 c0 {2,D} {7,S} {12,S}
6  C u0 p0 c0 {3,S} {7,D} {14,S}
7  C u0 p0 c0 {5,S} {6,D} {13,S}
8  O u0 p2 c0 {3,S} {16,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.422583,0.0455516,3.20125e-05,-8.1122e-08,3.76657e-11,-18202.6,26.0329], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[15.933,0.0270112,-9.94487e-06,1.62967e-09,-9.85133e-14,-23592.1,-59.7328], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""AVG CRESOL6/87""",
    longDesc = 
u"""
AVG CRESOL6/87.
CC1=CC=CC(O)=C1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 140,
    label = "PXC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {6,D} {13,S}
5  C u1 p0 c0 {3,S} {14,S} {15,S}
6  C u0 p0 c0 {4,D} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.810229,0.0688802,-4.93843e-05,1.71881e-08,-1.64486e-12,17071.9,33.2599], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[13.6624,0.0278801,-8.49489e-06,1.24384e-09,-7.17067e-14,13141,-41.2369], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/22/ 7 THERG""",
    longDesc = 
u"""
1/22/ 7 THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 141,
    label = "PXC6H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {7,S} {8,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  C u1 p0 c0 {4,S} {18,S} {19,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.204871,0.0683801,-4.14448e-05,1.26156e-08,-1.5312e-12,1832.8,31.6075], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[18.5385,0.0283108,-9.65307e-06,1.49548e-09,-8.66336e-14,-5092.99,-70.4491], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 142,
    label = "C6H5CH2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,D} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {1,D} {6,S} {12,S}
4  C u0 p0 c0 {2,D} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {3,S} {5,D} {11,S}
7  C u1 p0 c0 {1,S} {13,S} {14,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.481115,0.0385128,3.28615e-05,-7.69727e-08,3.54231e-11,23307,23.5488], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[14.044,0.0234939,-8.53754e-06,1.38908e-09,-8.36144e-14,18564.2,-51.6656], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T08/90""",
    longDesc = 
u"""
T08/90.
[CH2]C1C=CC=CC=1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 143,
    label = "C6H5CH3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u0 p0 c0 {2,S} {5,D} {11,S}
4  C u0 p0 c0 {2,D} {7,S} {15,S}
5  C u0 p0 c0 {3,D} {6,S} {12,S}
6  C u0 p0 c0 {5,S} {7,D} {13,S}
7  C u0 p0 c0 {4,S} {6,D} {14,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.61527,0.0210994,8.5366e-05,-1.32611e-07,5.59566e-11,4075.63,20.2822], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[12.94,0.0266913,-9.68385e-06,1.57386e-09,-9.46636e-14,-697.649,-46.7288], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 6/87""",
    longDesc = 
u"""
L 6/87.
CC1C=CC=CC=1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 144,
    label = "lC5H7",
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
            NASAPolynomial(coeffs=[-4.09743,0.061832,-4.87708e-05,1.66964e-08,-7.53349e-13,23683.6,45.1481], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.22465,0.0396013,-2.23456e-05,6.06497e-09,-6.384e-13,22303.4,14.01], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""HWZD99""",
    longDesc = 
u"""
HWZD99
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 145,
    label = "S2XC5H11",
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
            NASAPolynomial(coeffs=[0.498944,0.050985,-2.40687e-05,3.59465e-09,3.01383e-13,3407.02,27.8601], Tmin=(298,'K'), Tmax=(1377,'K')),
            NASAPolynomial(coeffs=[15.0998,0.0237225,-8.01389e-06,1.23431e-09,-7.123e-14,-2334.2,-52.9614], Tmin=(1377,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 146,
    label = "sC4H9",
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
            NASAPolynomial(coeffs=[0.694284,0.0331133,6.29426e-06,-2.70253e-08,1.19893e-11,6417.57,26.2798], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[9.42638,0.021919,-7.28684e-06,1.06303e-09,-5.56495e-14,3196.59,-22.4061], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""P11/94""",
    longDesc = 
u"""
P11/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 147,
    label = "OC12OOH",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {8,S} {17,S} {18,S}
2  C u0 p0 c0 {1,S} {3,S} {19,S} {20,S}
3  C u0 p0 c0 {2,S} {4,S} {21,S} {22,S}
4  C u0 p0 c0 {3,S} {5,S} {23,S} {24,S}
5  C u0 p0 c0 {4,S} {6,S} {25,S} {26,S}
6  C u0 p0 c0 {5,S} {7,S} {27,S} {28,S}
7  C u0 p0 c0 {6,S} {9,S} {29,S} {30,S}
8  C u0 p0 c0 {1,S} {11,S} {15,S} {16,S}
9  C u0 p0 c0 {7,S} {12,S} {31,S} {32,S}
10 C u0 p0 c0 {12,S} {13,S} {36,S} {37,S}
11 C u0 p0 c0 {8,S} {33,S} {34,S} {35,S}
12 C u0 p0 c0 {9,S} {10,S} {38,D}
13 O u0 p2 c0 {10,S} {14,S}
14 O u0 p2 c0 {13,S} {39,S}
15 H u0 p0 c0 {8,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {3,S}
23 H u0 p0 c0 {4,S}
24 H u0 p0 c0 {4,S}
25 H u0 p0 c0 {5,S}
26 H u0 p0 c0 {5,S}
27 H u0 p0 c0 {6,S}
28 H u0 p0 c0 {6,S}
29 H u0 p0 c0 {7,S}
30 H u0 p0 c0 {7,S}
31 H u0 p0 c0 {9,S}
32 H u0 p0 c0 {9,S}
33 H u0 p0 c0 {11,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {10,S}
37 H u0 p0 c0 {10,S}
38 O u0 p2 c0 {12,D}
39 H u0 p0 c0 {14,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.815247,0.161762,-0.000115472,4.29871e-08,-6.63427e-12,-62108,49.0272], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[45.827,0.0536022,-1.85842e-05,2.91165e-09,-1.70003e-13,-78495.2,-201.938], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/23/ 7 THERM""",
    longDesc = 
u"""
7/23/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCCCCCCCC(=O)COO
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 148,
    label = "PXC10H19",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {5,S} {17,S} {18,S}
4  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {6,S} {19,S} {20,S}
6  C u0 p0 c0 {5,S} {8,S} {23,S} {24,S}
7  C u0 p0 c0 {4,S} {9,S} {21,S} {22,S}
8  C u0 p0 c0 {6,S} {10,D} {25,S}
9  C u1 p0 c0 {7,S} {26,S} {27,S}
10 C u0 p0 c0 {8,D} {28,S} {29,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.4934,0.116496,-7.73332e-05,2.62845e-08,-3.64631e-12,-1884.01,46.2585], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[30.6443,0.0423089,-1.45032e-05,2.25521e-09,-1.30991e-13,-13745.3,-132.907], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/22/ 7 THERM""",
    longDesc = 
u"""
1/22/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCCCCC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 149,
    label = "PC12H25O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {10,S} {16,S} {17,S}
2  C u0 p0 c0 {1,S} {3,S} {18,S} {19,S}
3  C u0 p0 c0 {2,S} {4,S} {20,S} {21,S}
4  C u0 p0 c0 {3,S} {5,S} {22,S} {23,S}
5  C u0 p0 c0 {4,S} {6,S} {24,S} {25,S}
6  C u0 p0 c0 {5,S} {7,S} {26,S} {27,S}
7  C u0 p0 c0 {6,S} {8,S} {28,S} {29,S}
8  C u0 p0 c0 {7,S} {9,S} {30,S} {31,S}
9  C u0 p0 c0 {8,S} {11,S} {32,S} {33,S}
10 C u0 p0 c0 {1,S} {12,S} {14,S} {15,S}
11 C u0 p0 c0 {9,S} {13,S} {34,S} {35,S}
12 C u0 p0 c0 {10,S} {36,S} {37,S} {38,S}
13 O u0 p2 c0 {11,S} {39,S}
14 H u0 p0 c0 {10,S}
15 H u0 p0 c0 {10,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {4,S}
23 H u0 p0 c0 {4,S}
24 H u0 p0 c0 {5,S}
25 H u0 p0 c0 {5,S}
26 H u0 p0 c0 {6,S}
27 H u0 p0 c0 {6,S}
28 H u0 p0 c0 {7,S}
29 H u0 p0 c0 {7,S}
30 H u0 p0 c0 {8,S}
31 H u0 p0 c0 {8,S}
32 H u0 p0 c0 {9,S}
33 H u0 p0 c0 {9,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {12,S}
37 H u0 p0 c0 {12,S}
38 H u0 p0 c0 {12,S}
39 O u1 p2 c0 {13,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.413302,0.147893,-9.69755e-05,3.35724e-08,-4.91455e-12,-33684.4,43.2259], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[41.0721,0.0576224,-1.98682e-05,3.10098e-09,-1.80567e-13,-48496.1,-177.177], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/23/ 7 THERM""",
    longDesc = 
u"""
7/23/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCCCCCCCCCO[O]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 150,
    label = "C10H20",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {4,S} {17,S} {18,S}
4  C u0 p0 c0 {3,S} {6,S} {19,S} {20,S}
5  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
6  C u0 p0 c0 {4,S} {8,S} {21,S} {22,S}
7  C u0 p0 c0 {5,S} {9,S} {23,S} {24,S}
8  C u0 p0 c0 {6,S} {25,S} {26,S} {27,S}
9  C u0 p0 c0 {7,S} {10,D} {28,S}
10 C u0 p0 c0 {9,D} {29,S} {30,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
30 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.42902,0.119306,-7.94489e-05,2.72737e-08,-3.82718e-12,-18870.8,47.0571], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[30.8754,0.0439972,-1.49426e-05,2.30918e-09,-1.33551e-13,-30693.7,-132.705], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/22/ 7 THERM""",
    longDesc = 
u"""
1/22/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 151,
    label = "CH2CHO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {4,S} {5,S}
2 C u0 p0 c0 {1,D} {3,S} {6,S}
3 O u1 p2 c0 {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.40906,0.0107386,1.89149e-06,-7.15858e-09,2.86739e-12,62,9.57145], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.97567,0.00813059,-2.74362e-06,4.0703e-10,-2.17602e-14,-969.5,-5.03209], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""D05/83""",
    longDesc = 
u"""
D05/83
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[O]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 152,
    label = "tC4H9",
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
            NASAPolynomial(coeffs=[0.961676,0.0257359,1.5609e-05,-2.66565e-08,8.9418e-12,4656.44,24.8054], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.66073,0.0238794,-8.08904e-06,1.20575e-09,-6.50098e-14,1620.76,-14.8003], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""P11/94""",
    longDesc = 
u"""
P11/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C](C)C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 153,
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
            NASAPolynomial(coeffs=[2.96963,0.0244422,-9.12514e-06,-4.24669e-18,1.63047e-21,35503.3,12.0361], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[14.5382,-0.00856771,2.35595e-05,-1.36764e-08,2.44369e-12,33259.1,-45.3695], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""H6W/94""",
    longDesc = 
u"""
H6W/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C#CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 154,
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 155,
    label = "cC5H9",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
5  C u1 p0 c0 {3,S} {4,S} {14,S}
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
            NASAPolynomial(coeffs=[0.294271,0.0138234,9.08477e-05,-1.30087e-07,5.30518e-11,12565.7,27.3898], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[11.4068,0.022564,-7.02356e-06,1.1322e-09,-7.34382e-14,7526.88,-39.6363], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T03/97""",
    longDesc = 
u"""
T03/97.
[CH]1CCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 156,
    label = "NC12H26",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {9,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {17,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
5  C u0 p0 c0 {4,S} {6,S} {23,S} {24,S}
6  C u0 p0 c0 {5,S} {7,S} {25,S} {26,S}
7  C u0 p0 c0 {6,S} {8,S} {27,S} {28,S}
8  C u0 p0 c0 {7,S} {10,S} {29,S} {30,S}
9  C u0 p0 c0 {1,S} {11,S} {13,S} {14,S}
10 C u0 p0 c0 {8,S} {12,S} {31,S} {32,S}
11 C u0 p0 c0 {9,S} {33,S} {34,S} {35,S}
12 C u0 p0 c0 {10,S} {36,S} {37,S} {38,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {4,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {5,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {6,S}
27 H u0 p0 c0 {7,S}
28 H u0 p0 c0 {7,S}
29 H u0 p0 c0 {8,S}
30 H u0 p0 c0 {8,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {11,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {12,S}
37 H u0 p0 c0 {12,S}
38 H u0 p0 c0 {12,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.62182,0.147238,-9.4397e-05,3.07441e-08,-4.03602e-12,-40065.4,50.0995], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[38.5095,0.056355,-1.91493e-05,2.96025e-09,-1.71244e-13,-54884.3,-172.671], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCCCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 157,
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
            NASAPolynomial(coeffs=[4.3018,-0.00474912,2.11583e-05,-2.42764e-08,9.29225e-12,294.808,3.71666], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.01721,0.00223982,-6.33658e-07,1.14246e-10,-1.07909e-14,111.857,3.7851], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L 5/89""",
    longDesc = 
u"""
L 5/89.
[O]O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 158,
    label = "OC6H4CH3",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {4,D}
3  C u1 p0 c0 {2,S} {5,S} {11,S}
4  C u0 p0 c0 {2,D} {7,S} {15,S}
5  C u0 p0 c0 {3,S} {6,D} {13,S}
6  C u0 p0 c0 {5,D} {7,S} {14,S}
7  C u0 p0 c0 {4,S} {6,S} {12,D}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 O u0 p2 c0 {7,D}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.288558,0.0480035,1.8033e-05,-6.17415e-08,2.88526e-11,-689.456,26.7201], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[22.6094,0.00756462,6.59609e-06,-4.71509e-09,8.04091e-13,-8202.52,-97.2925], Tmin=(1000,'K'), Tmax=(2500,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2500,'K'),
    ),
    shortDesc = u"""EST/BUR P 1/93""",
    longDesc = 
u"""
EST/BUR P 1/93
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1[CH]C=CC(=O)C=1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 159,
    label = "nC4H3",
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
            NASAPolynomial(coeffs=[0.816677,0.0387162,-4.80457e-05,3.20668e-08,-8.56282e-12,64455.8,19.7405], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.80457,0.0107124,-4.19391e-06,7.04463e-10,-3.62713e-14,62987.8,-14.1297], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""USC/07""",
    longDesc = 
u"""
USC/07
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=CC#C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 160,
    label = "nC4H5",
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
            NASAPolynomial(coeffs=[0.226113,0.0367424,-2.21205e-05,1.43901e-09,2.64358e-12,42428.4,24.0664], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.40873,0.0177527,-7.56015e-06,1.42038e-09,-9.11002e-14,40438.8,-13.15], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""USC/07""",
    longDesc = 
u"""
USC/07
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]=CC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 161,
    label = "S5XC12H25",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
2  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
3  C u0 p0 c0 {2,S} {7,S} {23,S} {24,S}
4  C u0 p0 c0 {1,S} {8,S} {17,S} {18,S}
5  C u0 p0 c0 {2,S} {9,S} {19,S} {20,S}
6  C u0 p0 c0 {1,S} {10,S} {13,S} {14,S}
7  C u0 p0 c0 {3,S} {11,S} {25,S} {26,S}
8  C u0 p0 c0 {4,S} {12,S} {27,S} {28,S}
9  C u0 p0 c0 {5,S} {12,S} {29,S} {30,S}
10 C u0 p0 c0 {6,S} {31,S} {32,S} {33,S}
11 C u0 p0 c0 {7,S} {34,S} {35,S} {36,S}
12 C u1 p0 c0 {8,S} {9,S} {37,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {2,S}
22 H u0 p0 c0 {2,S}
23 H u0 p0 c0 {3,S}
24 H u0 p0 c0 {3,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {11,S}
37 H u0 p0 c0 {12,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.36787,0.137355,-8.24076e-05,2.36422e-08,-2.47436e-12,-16766.1,48.3522], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[37.9688,0.0538719,-1.82171e-05,2.80775e-09,-1.62108e-13,-31214.5,-165.806], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC[CH]CCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 162,
    label = "P12OOHX2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {7,S} {17,S} {18,S}
2  C u0 p0 c0 {1,S} {3,S} {19,S} {20,S}
3  C u0 p0 c0 {2,S} {4,S} {21,S} {22,S}
4  C u0 p0 c0 {3,S} {5,S} {23,S} {24,S}
5  C u0 p0 c0 {4,S} {6,S} {25,S} {26,S}
6  C u0 p0 c0 {5,S} {8,S} {27,S} {28,S}
7  C u0 p0 c0 {1,S} {11,S} {15,S} {16,S}
8  C u0 p0 c0 {6,S} {12,S} {29,S} {30,S}
9  C u0 p0 c0 {10,S} {12,S} {31,S} {32,S}
10 C u0 p0 c0 {9,S} {13,S} {33,S} {34,S}
11 C u0 p0 c0 {7,S} {35,S} {36,S} {37,S}
12 C u1 p0 c0 {8,S} {9,S} {38,S}
13 O u0 p2 c0 {10,S} {14,S}
14 O u0 p2 c0 {13,S} {39,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {3,S}
23 H u0 p0 c0 {4,S}
24 H u0 p0 c0 {4,S}
25 H u0 p0 c0 {5,S}
26 H u0 p0 c0 {5,S}
27 H u0 p0 c0 {6,S}
28 H u0 p0 c0 {6,S}
29 H u0 p0 c0 {8,S}
30 H u0 p0 c0 {8,S}
31 H u0 p0 c0 {9,S}
32 H u0 p0 c0 {9,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {10,S}
35 H u0 p0 c0 {11,S}
36 H u0 p0 c0 {11,S}
37 H u0 p0 c0 {11,S}
38 H u0 p0 c0 {12,S}
39 H u0 p0 c0 {14,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.160721,0.149326,-9.67267e-05,3.17154e-08,-4.20449e-12,-27556.9,46.3656], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[42.9995,0.0550834,-1.88956e-05,2.93978e-09,-1.70824e-13,-43065.8,-185.857], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/23/ 7 THERM""",
    longDesc = 
u"""
7/23/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCCCCCC[CH]CCOO
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 163,
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
            NASAPolynomial(coeffs=[0.983498,0.0336515,-1.10542e-07,-3.67434e-08,2.31412e-11,29626,16.5855], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.47439,0.0160127,-6.48231e-09,-3.58197e-09,9.23651e-13,28086,-16.133], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = u"""T12/89""",
    longDesc = 
u"""
T12/89
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]1C=CC=C1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 164,
    label = "NC6H14",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {18,S} {19,S} {20,S}
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
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.969606,0.0729086,-4.38854e-05,1.32313e-08,-1.58437e-12,-22780.4,32.307], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[18.9634,0.030448,-1.03795e-05,1.60775e-09,-9.3127e-14,-30162.9,-76.2839], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 165,
    label = "S5XC11H23",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
2  C u0 p0 c0 {4,S} {6,S} {20,S} {21,S}
3  C u0 p0 c0 {1,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {2,S} {8,S} {18,S} {19,S}
5  C u0 p0 c0 {1,S} {9,S} {12,S} {13,S}
6  C u0 p0 c0 {2,S} {10,S} {22,S} {23,S}
7  C u0 p0 c0 {3,S} {11,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {11,S} {26,S} {27,S}
9  C u0 p0 c0 {5,S} {28,S} {29,S} {30,S}
10 C u0 p0 c0 {6,S} {31,S} {32,S} {33,S}
11 C u1 p0 c0 {7,S} {8,S} {34,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {2,S}
21 H u0 p0 c0 {2,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {9,S}
31 H u0 p0 c0 {10,S}
32 H u0 p0 c0 {10,S}
33 H u0 p0 c0 {10,S}
34 H u0 p0 c0 {11,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.1025,0.125022,-7.40802e-05,2.07819e-08,-2.07856e-12,-13884,45.431], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[34.7028,0.0495634,-1.67589e-05,2.58285e-09,-1.49119e-13,-27089.2,-149.691], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/22/ 7 THERM""",
    longDesc = 
u"""
1/22/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC[CH]CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 166,
    label = "C5H8-13",
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
            NASAPolynomial(coeffs=[3.44317,0.033499,-7.9124e-06,-6.55475e-10,-3.17004e-13,6892.23,9.11999], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[11.5705,0.0240576,-8.98579e-06,1.53219e-09,-9.89135e-14,3252.9,-37.8272], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERGAS G""",
    longDesc = 
u"""
THERGAS G
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC=CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 167,
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 168,
    label = "C5H8-14",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {1,S} {5,D} {9,S}
4  C u0 p0 c0 {2,D} {10,S} {11,S}
5  C u0 p0 c0 {3,D} {12,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.08677,0.0451422,-2.39204e-05,3.98173e-09,4.99385e-13,10626.9,21.4715], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[11.9449,0.0209191,-7.13315e-06,1.1413e-09,-7.05198e-14,7108.63,-36.6149], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/29/ 8 G""",
    longDesc = 
u"""
8/29/ 8 G
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 169,
    label = "C4H7",
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
            NASAPolynomial(coeffs=[0.744494,0.0396789,-2.28981e-05,2.1353e-09,2.30964e-12,22653.3,23.4379], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.01348,0.0226346,-9.25455e-06,1.68079e-09,-1.04086e-13,20955,-8.88931], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""USC/07""",
    longDesc = 
u"""
USC/07
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 170,
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
            NASAPolynomial(coeffs=[0.112845,0.034369,-1.11074e-05,-9.21067e-09,6.20652e-12,11802.3,23.09], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.86731,0.0149187,-3.15487e-06,-4.18413e-10,1.57613e-13,9133.85,-23.3282], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""H6W/94""",
    longDesc = 
u"""
H6W/94
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 171,
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 172,
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
            NASAPolynomial(coeffs=[0.808681,0.0233616,-3.55172e-05,2.80152e-08,-8.50073e-12,26429,13.9397], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.14757,0.00596167,-2.37295e-06,4.67412e-10,-3.61235e-14,25936,-1.23028], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L 1/91""",
    longDesc = 
u"""
L 1/91.
C#C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 173,
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 174,
    label = "SXC9H19",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {14,S} {15,S}
3  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
4  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {9,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {9,S} {25,S} {26,S} {27,S}
9  C u1 p0 c0 {6,S} {8,S} {28,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.576046,0.100512,-5.77755e-05,1.53277e-08,-1.35057e-12,-8122.9,39.5796], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[28.0393,0.041144,-1.3926e-05,2.14746e-09,-1.24022e-13,-18772.8,-116.697], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 175,
    label = "PXC10H21",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {4,S} {17,S} {18,S}
4  C u0 p0 c0 {3,S} {5,S} {19,S} {20,S}
5  C u0 p0 c0 {4,S} {7,S} {21,S} {22,S}
6  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {5,S} {9,S} {23,S} {24,S}
8  C u0 p0 c0 {6,S} {10,S} {25,S} {26,S}
9  C u0 p0 c0 {7,S} {27,S} {28,S} {29,S}
10 C u1 p0 c0 {8,S} {30,S} {31,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.31358,0.117973,-7.51843e-05,2.43331e-08,-3.17523e-12,-9689.68,43.501], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[31.5697,0.0455818,-1.54995e-05,2.39711e-09,-1.3871e-13,-21573.8,-134.709], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 176,
    label = "CH3CHCHCHO",
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
            NASAPolynomial(coeffs=[-1.55578,0.0409641,-1.69869e-05,-6.00928e-18,2.31369e-21,-14139.5,37.4708], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[19.8795,-0.0209131,4.45361e-05,-2.60375e-08,4.86836e-12,-19527.9,-68.72], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""T 5/92""",
    longDesc = 
u"""
T 5/92
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
CC=CC=O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 177,
    label = "S2XC8H17",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {8,S} {17,S} {18,S}
5  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {22,S} {23,S} {24,S}
7  C u0 p0 c0 {5,S} {19,S} {20,S} {21,S}
8  C u1 p0 c0 {4,S} {5,S} {25,S}
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
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.304233,0.0880077,-4.90743e-05,1.21858e-08,-8.87773e-13,-5237.93,36.6583], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[24.9044,0.0366394,-1.2385e-05,1.90835e-09,-1.10161e-13,-14713.5,-101.345], Tmin=(1383,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 178,
    label = "PXC9H19",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {3,S} {14,S} {15,S}
3  C u0 p0 c0 {2,S} {4,S} {16,S} {17,S}
4  C u0 p0 c0 {3,S} {6,S} {18,S} {19,S}
5  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
6  C u0 p0 c0 {4,S} {8,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {9,S} {22,S} {23,S}
8  C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
9  C u1 p0 c0 {7,S} {27,S} {28,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.04387,0.105617,-6.682e-05,2.14486e-08,-2.77404e-12,-6808.19,42.3519], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[28.3098,0.0412657,-1.40383e-05,2.17175e-09,-1.25692e-13,-17451.6,-116.838], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 179,
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
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 180,
    label = "C6H5CHO",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,D} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {1,D} {6,S} {12,S}
4  C u0 p0 c0 {2,D} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {3,S} {5,D} {11,S}
7  C u0 p0 c0 {1,S} {13,D} {14,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {7,D}
14 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.16273,0.0663692,-3.48164e-05,-6.29994e-09,8.58071e-12,-6116.93,40.2317], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[13.6507,0.0256804,-1.04667e-05,1.94134e-09,-1.34838e-13,-11019.7,-47.9658], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""L 3/86""",
    longDesc = 
u"""
L 3/86
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
O=CC1C=CC=CC=1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 181,
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
            NASAPolynomial(coeffs=[2.13584,0.0181189,-1.73947e-05,9.34398e-09,-2.01458e-12,-7270,12.2156], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.5113,0.0090036,-4.1694e-06,9.23346e-10,-7.94838e-14,-7778.5,0.632247], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""D05/90""",
    longDesc = 
u"""
D05/90.
C=C=O
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 182,
    label = "S2XC12H25",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {17,S} {18,S}
3  C u0 p0 c0 {2,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {21,S} {22,S}
5  C u0 p0 c0 {4,S} {7,S} {23,S} {24,S}
6  C u0 p0 c0 {1,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {5,S} {10,S} {25,S} {26,S}
8  C u0 p0 c0 {6,S} {12,S} {29,S} {30,S}
9  C u0 p0 c0 {11,S} {12,S} {27,S} {28,S}
10 C u0 p0 c0 {7,S} {34,S} {35,S} {36,S}
11 C u0 p0 c0 {9,S} {31,S} {32,S} {33,S}
12 C u1 p0 c0 {8,S} {9,S} {37,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {4,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {5,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {8,S}
30 H u0 p0 c0 {8,S}
31 H u0 p0 c0 {11,S}
32 H u0 p0 c0 {11,S}
33 H u0 p0 c0 {11,S}
34 H u0 p0 c0 {10,S}
35 H u0 p0 c0 {10,S}
36 H u0 p0 c0 {10,S}
37 H u0 p0 c0 {12,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.36787,0.137355,-8.24076e-05,2.36422e-08,-2.47436e-12,-16766.1,48.3522], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[37.9688,0.0538719,-1.82171e-05,2.80775e-09,-1.62108e-13,-31214.5,-165.806], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""1/ 2/ 7 THERM""",
    longDesc = 
u"""
1/ 2/ 7 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]CCCCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 183,
    label = "C6H5CO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,D} {7,S}
2  C u0 p0 c0 {1,S} {4,D} {8,S}
3  C u0 p0 c0 {1,D} {6,S} {12,S}
4  C u0 p0 c0 {2,D} {5,S} {9,S}
5  C u0 p0 c0 {4,S} {6,D} {10,S}
6  C u0 p0 c0 {3,S} {5,D} {11,S}
7  C u1 p0 c0 {1,S} {13,D}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {3,S}
13 O u0 p2 c0 {7,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.02512,0.0615125,-3.16037e-05,-6.97246e-09,7.98351e-12,11255.8,35.7782], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[13.3744,0.0239993,-1.04657e-05,2.16691e-09,-1.8007e-13,6914.78,-44.6592], Tmin=(1000,'K'), Tmax=(2500,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2500,'K'),
    ),
    shortDesc = u"""EST/BUR P 1/93""",
    longDesc = 
u"""
EST/BUR P 1/93
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
O=[C]C1C=CC=CC=1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 184,
    label = "AR",
    molecule = 
"""
1 Ar u0 p4 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.366], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.366], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""120186""",
    longDesc = 
u"""
120186
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[Ar]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 185,
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
    shortDesc = u"""L10/90""",
    longDesc = 
u"""
L10/90.
[He]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 186,
    label = "OH*",
    molecule = 
"""
multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.46084,0.000501872,-2.00254e-06,3.18902e-09,-1.35452e-12,50734.9,1.73976], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.75583,0.00139849,-4.19428e-07,6.33453e-11,-3.56042e-15,50975.2,5.62581], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT A""",
    longDesc = 
u"""
ATcT A.
Duplicate of species OH (i.e. same molecular structure according to RMG)
[OH]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 187,
    label = "CH*",
    molecule = 
"""
multiplicity 4
1 C u3 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.4725,0.000426444,-1.95182e-06,3.51755e-09,-1.60436e-12,104335,1.448], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.78221,0.00147247,-4.63436e-07,7.32736e-11,-4.19705e-15,104547,5.17421], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""TPIS79""",
    longDesc = 
u"""
TPIS79.
Duplicate of species CH (i.e. same molecular structure according to RMG)
[CH]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 188,
    label = "SAXC4H7",
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
            NASAPolynomial(coeffs=[0.477014,0.0361447,-1.63111e-05,-1.57565e-10,1.72942e-12,14199.3,23.2006], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[10.0858,0.0156299,-4.4259e-06,6.21461e-10,-3.51689e-14,11044.4,-28.4205], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[CH]C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 189,
    label = "CH3cC6H11",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {17,S} {18,S}
4  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {5,S} {15,S} {16,S}
7  C u0 p0 c0 {1,S} {19,S} {20,S} {21,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.90849,0.0969227,-5.76086e-05,1.48744e-08,-1.11358e-12,-19264.3,65.7805], Tmin=(298,'K'), Tmax=(1381,'K')),
            NASAPolynomial(coeffs=[22.0212,0.0332077,-1.15858e-05,1.82325e-09,-1.06797e-13,-30769.4,-103.212], Tmin=(1381,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/28/ 9 THERM""",
    longDesc = 
u"""
7/28/ 9 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 190,
    label = "PXCH2cC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {17,S} {18,S}
4  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {5,S} {15,S} {16,S}
7  C u1 p0 c0 {1,S} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.11266,0.0922538,-5.49739e-05,1.41483e-08,-1.03789e-12,5344.31,63.8409], Tmin=(298,'K'), Tmax=(1380,'K')),
            NASAPolynomial(coeffs=[21.5892,0.031077,-1.08618e-05,1.71136e-09,-1.00327e-13,-5697.06,-98.4332], Tmin=(1380,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 191,
    label = "cC6H12",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {5,S} {17,S} {18,S}
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
18 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.09894,0.0821969,-4.45302e-05,8.07342e-09,5.67231e-13,-15694.9,59.3353], Tmin=(298,'K'), Tmax=(1375,'K')),
            NASAPolynomial(coeffs=[10.2209,0.0414192,-1.6003e-05,2.79457e-09,-1.83461e-13,-21665,-38.752], Tmin=(1375,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 192,
    label = "PXCH2cC5H9",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {4,S} {12,S} {13,S}
6  C u1 p0 c0 {1,S} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.01089,0.0859363,-7.12283e-05,3.50314e-08,-8.11113e-12,10160.6,59.1252], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[15.1875,0.0262978,-7.97774e-06,1.17628e-09,-6.87891e-14,3553.22,-56.9284], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C1CCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 193,
    label = "cC6H10O2H-2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {6,S} {7,S} {9,S}
2  C u0 p0 c0 {3,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {2,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {5,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {6,S} {16,S} {17,S}
6  C u1 p0 c0 {1,S} {5,S} {18,S}
7  O u0 p2 c0 {1,S} {8,S}
8  O u0 p2 c0 {7,S} {19,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.91892,0.0890964,-5.70299e-05,1.65057e-08,-1.68121e-12,-5689.08,57.9514], Tmin=(298,'K'), Tmax=(1674,'K')),
            NASAPolynomial(coeffs=[20.8647,0.0296342,-1.07902e-05,1.75532e-09,-1.05424e-13,-14698.9,-86.1399], Tmin=(1674,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""therm""",
    longDesc = 
u"""
therm
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
OOC1[CH]CCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 194,
    label = "cC6H11O2",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,S} {17,S} {18,S}
4  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
5  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
6  C u0 p0 c0 {3,S} {5,S} {15,S} {16,S}
7  O u0 p2 c0 {1,S} {19,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 O u1 p2 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.91664,0.0836157,-4.98339e-05,1.27484e-08,-9.18574e-13,-12318.2,50.7735], Tmin=(298,'K'), Tmax=(1379,'K')),
            NASAPolynomial(coeffs=[22.3724,0.0276799,-9.74217e-06,1.54195e-09,-9.06768e-14,-22504.3,-98.4391], Tmin=(1379,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""therm""",
    longDesc = 
u"""
therm
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OC1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 195,
    label = "C7H12-16",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,D} {14,S}
5  C u0 p0 c0 {3,S} {7,D} {15,S}
6  C u0 p0 c0 {4,D} {16,S} {17,S}
7  C u0 p0 c0 {5,D} {18,S} {19,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.97047,0.0790572,-5.6207e-05,2.09665e-08,-3.22821e-12,5205.22,39.6986], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[19.9812,0.0273257,-9.2759e-06,1.43298e-09,-8.28566e-14,-2379.88,-78.0467], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCCC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 196,
    label = "C7H14-2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {7,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {7,D} {21,S}
7  C u0 p0 c0 {5,S} {6,D} {20,S}
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
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.16533,0.079044,-4.96102e-05,1.58569e-08,-2.05346e-12,-11736.2,38.1864], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[20.6192,0.0314853,-1.07162e-05,1.65828e-09,-9.59912e-14,-19671.3,-80.0526], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 197,
    label = "C8H16-3",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {8,S} {15,S} {16,S}
4  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {20,S} {21,S} {22,S}
6  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {8,D} {23,S}
8  C u0 p0 c0 {3,S} {7,D} {24,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.21774,0.0882852,-4.86946e-05,9.3159e-09,6.45244e-13,-14627.3,37.6708], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[20.5157,0.038473,-1.24391e-05,1.91727e-09,-1.15482e-13,-21538,-78.1175], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC=CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 198,
    label = "SOOcC6O2H",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {7,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {6,S} {18,S} {19,S}
5  C u0 p0 c0 {3,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {4,S} {5,S} {16,S} {17,S}
7  O u0 p2 c0 {1,S} {9,S}
8  O u0 p2 c0 {2,S} {20,S}
9  O u0 p2 c0 {7,S} {21,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 O u1 p2 c0 {8,S}
21 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.28956,0.0963507,-6.55057e-05,2.09993e-08,-2.51429e-12,-25652.1,46.803], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[28.1644,0.0270788,-9.55221e-06,1.51424e-09,-8.91458e-14,-36875,-123.417], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""therm""",
    longDesc = 
u"""
therm
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OC1CCCCC1OO
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 199,
    label = "CH3-4-SAXC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {6,D} {17,S}
6  C u0 p0 c0 {5,D} {7,S} {18,S}
7  C u1 p0 c0 {6,S} {19,S} {20,S}
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
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.11469,0.0810101,-5.40865e-05,1.86778e-08,-2.67405e-12,6108.19,40.336], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[21.0132,0.0294229,-1.01841e-05,1.59373e-09,-9.29771e-14,-2240.32,-84.8628], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CC(C)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 200,
    label = "PXCH2-5-1C6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {7,D} {16,S}
6  C u1 p0 c0 {1,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
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
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.37028,0.0798095,-5.47529e-05,1.97898e-08,-2.98366e-12,14146.8,39.704], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[20.6259,0.0292698,-1.00239e-05,1.55747e-09,-9.04079e-14,6367.1,-78.8223], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)CCC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 201,
    label = "PXCH2-5-1C7H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {8,D} {19,S}
7  C u1 p0 c0 {1,S} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {22,S} {23,S}
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
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.58551,0.0919236,-6.27829e-05,2.24798e-08,-3.34343e-12,10854.6,42.3931], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[23.849,0.033629,-1.15018e-05,1.7856e-09,-1.03592e-13,1859.76,-94.6861], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(CC)CCC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 202,
    label = "PXCH2-5-1C8H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {19,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {9,D} {22,S}
8  C u1 p0 c0 {1,S} {23,S} {24,S}
9  C u0 p0 c0 {7,D} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.86646,0.104331,-7.12367e-05,2.54237e-08,-3.75769e-12,7572.17,45.386], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[27.0848,0.0379747,-1.29747e-05,2.0129e-09,-1.16725e-13,-2653.19,-110.625], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(CCC)CCC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 203,
    label = "PXCH2-5-1C9H17",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
5  C u0 p0 c0 {4,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {8,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {6,S} {10,D} {25,S}
9  C u1 p0 c0 {1,S} {26,S} {27,S}
10 C u0 p0 c0 {8,D} {28,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.18648,0.11691,-7.99349e-05,2.85126e-08,-4.20287e-12,5100.92,48.5598], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[30.3291,0.0423114,-1.44442e-05,2.23965e-09,-1.29826e-13,-6364.82,-126.614], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(CCC=C)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 204,
    label = "C3H5cC6H11",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {19,S} {20,S}
4  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {5,S} {17,S} {18,S}
7  C u0 p0 c0 {1,S} {8,S} {21,S} {22,S}
8  C u0 p0 c0 {7,S} {9,D} {23,S}
9  C u0 p0 c0 {8,D} {24,S} {25,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {9,S}
25 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-9.73004,0.118168,-7.58416e-05,2.28048e-08,-2.48241e-12,-10008.9,73.0984], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[27.3375,0.0382152,-1.32622e-05,2.07988e-09,-1.21544e-13,-23424,-128.131], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCC1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 205,
    label = "PX6-2C6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {14,S}
5  C u0 p0 c0 {3,S} {4,D} {15,S}
6  C u1 p0 c0 {2,S} {16,S} {17,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.106007,0.0620508,-3.86461e-05,1.22632e-08,-1.57993e-12,15754.6,33.3335], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[16.9431,0.0250225,-8.52515e-06,1.32009e-09,-7.64493e-14,9515.57,-59.2876], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC=CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 206,
    label = "PXCH2-4-1C5H9",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,D} {13,S}
5  C u1 p0 c0 {1,S} {14,S} {15,S}
6  C u0 p0 c0 {4,D} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.12566,0.07235,-5.28656e-05,2.01124e-08,-3.1193e-12,16518.2,40.0257], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[18.7248,0.0226408,-7.40593e-06,1.12566e-09,-6.45784e-14,9411.05,-71.5448], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8 G""",
    longDesc = 
u"""
9/ 8 G
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)CC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 207,
    label = "PX1-4C8H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
3  C u0 p0 c0 {4,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {3,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {1,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {7,D} {20,S}
7  C u0 p0 c0 {2,S} {6,D} {21,S}
8  C u1 p0 c0 {4,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.29095,0.0907951,-6.12886e-05,2.14028e-08,-3.04869e-12,11530.4,41.2455], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[23.8416,0.0331279,-1.1134e-05,1.70873e-09,-9.83549e-14,2724.59,-94.0576], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC=CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 208,
    label = "PXCH2-3-1C5H9",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {3,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {10,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {6,D} {13,S}
5  C u1 p0 c0 {1,S} {14,S} {15,S}
6  C u0 p0 c0 {4,D} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.18772,0.067972,-4.72509e-05,1.74233e-08,-2.68819e-12,16761.1,37.1393], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[17.3862,0.0248639,-8.51628e-06,1.3233e-09,-7.68158e-14,10227.7,-62.8048], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C=C)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 209,
    label = "PXCH2-3-1C4H7",
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
            NASAPolynomial(coeffs=[-0.4713,0.0637137,-5.73101e-05,3.08682e-08,-7.51184e-12,18554,27.9406], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[14.6956,0.0195986,-5.79712e-06,8.38458e-10,-4.83385e-14,14288.4,-50.3257], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""9/ 8 G""",
    longDesc = 
u"""
9/ 8 G
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)C=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 210,
    label = "CH3-2-C4H5-13",
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
            NASAPolynomial(coeffs=[-1.29615,0.0600088,-5.08076e-05,2.29593e-08,-4.20161e-12,7576.17,30.5152], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[14.1556,0.0181481,-6.08113e-06,9.30934e-10,-5.34767e-14,2751.51,-50.4831], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 211,
    label = "C2H3cC6H11",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {18,S} {19,S}
4  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {5,S} {16,S} {17,S}
7  C u0 p0 c0 {1,S} {8,D} {20,S}
8  C u0 p0 c0 {7,D} {21,S} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-9.56052,0.106426,-6.8516e-05,2.05745e-08,-2.22221e-12,-6991.03,70.5863], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[24.0767,0.0338603,-1.1779e-05,1.85017e-09,-1.08237e-13,-19154.9,-112.003], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 212,
    label = "C7H12-13",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {5,D} {15,S}
5  C u0 p0 c0 {4,D} {6,S} {17,S}
6  C u0 p0 c0 {5,S} {7,D} {16,S}
7  C u0 p0 c0 {6,D} {18,S} {19,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.71195,0.0856732,-6.7431e-05,2.8366e-08,-4.92153e-12,1830.55,41.9698], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[20.4739,0.0272296,-9.30926e-06,1.44474e-09,-8.37944e-14,-5865.38,-81.1539], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC=CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 213,
    label = "C3H7S2XcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {4,S} {19,S} {20,S}
4  C u0 p0 c0 {3,S} {5,S} {17,S} {18,S}
5  C u0 p0 c0 {4,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
7  C u0 p0 c0 {5,S} {9,S} {21,S} {22,S}
8  C u0 p0 c0 {6,S} {23,S} {24,S} {25,S}
9  C u1 p0 c0 {1,S} {7,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.16097,0.111676,-6.2366e-05,1.36019e-08,-3.71583e-13,-2137,69.7547], Tmin=(298,'K'), Tmax=(1374,'K')),
            NASAPolynomial(coeffs=[27.8484,0.0395366,-1.3645e-05,2.13266e-09,-1.24355e-13,-15668.3,-127.604], Tmin=(1374,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC1[CH]CCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 214,
    label = "C4H9S2XcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {16,S} {17,S}
3  C u0 p0 c0 {1,S} {5,S} {22,S} {23,S}
4  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {6,S} {20,S} {21,S}
6  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {6,S} {10,S} {24,S} {25,S}
9  C u0 p0 c0 {7,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {1,S} {8,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {3,S}
23 H u0 p0 c0 {3,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.43326,0.12404,-7.07409e-05,1.64921e-08,-7.74026e-13,-5018.04,72.7079], Tmin=(298,'K'), Tmax=(1376,'K')),
            NASAPolynomial(coeffs=[31.1137,0.0438458,-1.51034e-05,2.35758e-09,-1.37346e-13,-19793.2,-143.715], Tmin=(1376,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC1[CH]CCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 215,
    label = "PX1-3C7H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {6,S} {7,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {6,D} {18,S}
6  C u0 p0 c0 {3,S} {5,D} {17,S}
7  C u1 p0 c0 {3,S} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.28817,0.078204,-5.22552e-05,1.82428e-08,-2.64174e-12,13102.5,40.7391], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[20.2595,0.0294019,-1.00313e-05,1.55475e-09,-9.0097e-14,5423.36,-75.6006], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC=CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 216,
    label = "PXC2H4-4-1C6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {3,S} {8,D} {19,S}
7  C u1 p0 c0 {4,S} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {22,S} {23,S}
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
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.47916,0.092823,-6.50277e-05,2.40441e-08,-3.68237e-12,9838.52,42.2096], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[23.7873,0.0334277,-1.12875e-05,1.73745e-09,-1.00205e-13,1088.69,-93.3753], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(CC)CC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 217,
    label = "SAX6-4C10H19",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
3  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {8,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {9,S} {19,S} {20,S}
6  C u0 p0 c0 {3,S} {21,S} {22,S} {23,S}
7  C u0 p0 c0 {2,S} {24,S} {25,S} {26,S}
8  C u0 p0 c0 {4,S} {10,D} {28,S}
9  C u1 p0 c0 {5,S} {10,S} {27,S}
10 C u0 p0 c0 {8,D} {9,S} {29,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.86658,0.11671,-7.72536e-05,2.62645e-08,-3.6614e-12,-3176.94,50.2855], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[30.2417,0.0427836,-1.46951e-05,2.288e-09,-1.33013e-13,-15069.6,-128.827], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC=C[CH]CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 218,
    label = "SAX4-2C8H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {7,S} {18,S} {19,S} {20,S}
6  C u1 p0 c0 {3,S} {8,S} {21,S}
7  C u0 p0 c0 {5,S} {8,D} {22,S}
8  C u0 p0 c0 {6,S} {7,D} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.35404,0.0921118,-6.0654e-05,2.05583e-08,-2.87162e-12,2590.03,44.5285], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[23.7926,0.0340643,-1.17383e-05,1.83155e-09,-1.06634e-13,-6859.31,-97.0816], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C[CH]CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 219,
    label = "PX8-2C8H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {8,S} {13,S} {14,S}
5  C u0 p0 c0 {7,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {7,D} {20,S}
7  C u0 p0 c0 {5,S} {6,D} {21,S}
8  C u1 p0 c0 {4,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.29095,0.0907951,-6.12886e-05,2.14028e-08,-3.04869e-12,11530.4,41.2455], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[23.8416,0.0331279,-1.1134e-05,1.70873e-09,-9.83549e-14,2724.59,-94.0576], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCC=CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 220,
    label = "CH3S3XcC5H8",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {7,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
6  C u1 p0 c0 {3,S} {4,S} {17,S}
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
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.37989,0.0804871,-5.54172e-05,1.94407e-08,-2.75597e-12,10695.9,60.6844], Tmin=(298,'K'), Tmax=(1396,'K')),
            NASAPolynomial(coeffs=[15.7571,0.0274323,-9.38653e-06,1.45782e-09,-8.46056e-14,2592.4,-63.8757], Tmin=(1396,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1C[CH]CC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 221,
    label = "C4H9S4XcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
2  C u0 p0 c0 {1,S} {5,S} {16,S} {17,S}
3  C u0 p0 c0 {1,S} {7,S} {18,S} {19,S}
4  C u0 p0 c0 {1,S} {8,S} {20,S} {21,S}
5  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {5,S} {9,S} {12,S} {13,S}
7  C u0 p0 c0 {3,S} {10,S} {22,S} {23,S}
8  C u0 p0 c0 {4,S} {10,S} {24,S} {25,S}
9  C u0 p0 c0 {6,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {7,S} {8,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.43326,0.12404,-7.07409e-05,1.64921e-08,-7.74026e-13,-5018.04,72.7079], Tmin=(298,'K'), Tmax=(1376,'K')),
            NASAPolynomial(coeffs=[31.1137,0.0438458,-1.51034e-05,2.35758e-09,-1.37346e-13,-19793.2,-143.715], Tmin=(1376,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC1CC[CH]CC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 222,
    label = "C2H5cC6H11",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {20,S} {21,S}
4  C u0 p0 c0 {1,S} {8,S} {10,S} {11,S}
5  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {5,S} {7,S} {16,S} {17,S}
7  C u0 p0 c0 {3,S} {6,S} {18,S} {19,S}
8  C u0 p0 c0 {4,S} {22,S} {23,S} {24,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.36504,0.104031,-5.5678e-05,8.2959e-09,1.62359e-12,-22309,58.8226], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[18.9267,0.0460284,-1.65981e-05,2.76141e-09,-1.75252e-13,-30882.6,-82.0369], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 223,
    label = "CH3-5-SAX1C7H12",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {7,D} {20,S}
7  C u0 p0 c0 {6,D} {8,S} {21,S}
8  C u1 p0 c0 {7,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.67782,0.0948369,-6.48578e-05,2.30902e-08,-3.39537e-12,3267.18,44.604], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[23.8752,0.0340728,-1.1665e-05,1.81218e-09,-1.05186e-13,-6115.6,-98.5086], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CCC(C)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 224,
    label = "C3H7S4XcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {6,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {7,S} {17,S} {18,S}
5  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {3,S} {9,S} {19,S} {20,S}
7  C u0 p0 c0 {4,S} {9,S} {21,S} {22,S}
8  C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
9  C u1 p0 c0 {6,S} {7,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.16097,0.111676,-6.2366e-05,1.36019e-08,-3.71583e-13,-2137,69.7547], Tmin=(298,'K'), Tmax=(1374,'K')),
            NASAPolynomial(coeffs=[27.8484,0.0395366,-1.3645e-05,2.13266e-09,-1.24355e-13,-15668.3,-127.604], Tmin=(1374,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC1CC[CH]CC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 225,
    label = "C2H5-4-SAX1C6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {1,S} {7,D} {20,S}
7  C u0 p0 c0 {6,D} {8,S} {21,S}
8  C u1 p0 c0 {7,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.40917,0.0935024,-6.2723e-05,2.17703e-08,-3.12711e-12,2827.4,43.3875], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[24.2081,0.0338319,-1.16844e-05,1.82585e-09,-1.06411e-13,-6737.56,-100.572], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CC(CC)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 226,
    label = "C4H9cC6H11",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
2  C u0 p0 c0 {1,S} {5,S} {16,S} {17,S}
3  C u0 p0 c0 {1,S} {6,S} {18,S} {19,S}
4  C u0 p0 c0 {1,S} {8,S} {26,S} {27,S}
5  C u0 p0 c0 {2,S} {9,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {7,S} {20,S} {21,S}
7  C u0 p0 c0 {6,S} {8,S} {22,S} {23,S}
8  C u0 p0 c0 {4,S} {7,S} {24,S} {25,S}
9  C u0 p0 c0 {5,S} {10,S} {12,S} {13,S}
10 C u0 p0 c0 {9,S} {28,S} {29,S} {30,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {9,S}
13 H u0 p0 c0 {9,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {4,S}
27 H u0 p0 c0 {4,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
30 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-9.65592,0.133757,-8.24457e-05,2.34054e-08,-2.29382e-12,-28321.5,75.0048], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[31.7165,0.04625,-1.6005e-05,2.50511e-09,-1.4619e-13,-43495,-150.255], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 227,
    label = "C3H7cC6H11",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
2  C u0 p0 c0 {1,S} {8,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {5,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {7,S} {23,S} {24,S}
5  C u0 p0 c0 {3,S} {6,S} {17,S} {18,S}
6  C u0 p0 c0 {5,S} {7,S} {19,S} {20,S}
7  C u0 p0 c0 {4,S} {6,S} {21,S} {22,S}
8  C u0 p0 c0 {2,S} {9,S} {11,S} {12,S}
9  C u0 p0 c0 {8,S} {25,S} {26,S} {27,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {8,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {4,S}
24 H u0 p0 c0 {4,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.56398,0.115852,-6.28198e-05,1.01053e-08,1.53959e-12,-25204.7,61.4419], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[21.8211,0.0505114,-1.79565e-05,2.95858e-09,-1.86558e-13,-34742.1,-95.8203], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 228,
    label = "CH3-5-SAXC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {6,D} {17,S}
6  C u0 p0 c0 {5,D} {7,S} {18,S}
7  C u1 p0 c0 {6,S} {19,S} {20,S}
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
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.11469,0.0810101,-5.40865e-05,1.86778e-08,-2.67405e-12,6108.19,40.336], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[21.0132,0.0294229,-1.01841e-05,1.59373e-09,-9.29771e-14,-2240.32,-84.8628], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CCC(C)C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 229,
    label = "C9H18-4",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {8,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {9,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {2,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {9,D} {26,S}
9  C u0 p0 c0 {5,S} {8,D} {27,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.10105,0.0979164,-5.05751e-05,5.87235e-09,2.42065e-12,-17556.8,38.9535], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[23.3474,0.0430964,-1.38965e-05,2.13847e-09,-1.28697e-13,-25378.2,-91.568], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC=CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 230,
    label = "C6H10-15",
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
            NASAPolynomial(coeffs=[-0.242008,0.0584661,-3.14826e-05,4.43728e-09,1.21856e-12,7741.97,29.5274], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[14.5037,0.0248335,-8.10295e-06,1.25523e-09,-7.57703e-14,3073.7,-49.0055], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 231,
    label = "C10H20-5",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
3  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {8,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {9,S} {19,S} {20,S}
6  C u0 p0 c0 {2,S} {10,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {26,S} {27,S} {28,S}
9  C u0 p0 c0 {5,S} {10,D} {29,S}
10 C u0 p0 c0 {6,S} {9,D} {30,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.85367,0.119797,-7.98206e-05,2.7542e-08,-3.90606e-12,-20157.3,51.3122], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[30.4543,0.0444949,-1.51432e-05,2.34335e-09,-1.35652e-13,-32005.2,-128.511], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/31/ 9 THERM""",
    longDesc = 
u"""
7/31/ 9 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC=CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 232,
    label = "PAXCH2-2-C4H5",
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
            NASAPolynomial(coeffs=[-2.02197,0.063629,-5.92706e-05,2.8561e-08,-5.45166e-12,25864.1,32.6121], Tmin=(298,'K'), Tmax=(1400,'K')),
            NASAPolynomial(coeffs=[14.8523,0.0155764,-5.3169e-06,8.24224e-10,-4.77657e-14,20807.7,-55.0503], Tmin=(1400,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=C)C=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 233,
    label = "S2XC8H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
6  C u1 p0 c0 {2,S} {4,S} {20,S}
7  C u0 p0 c0 {3,S} {8,D} {21,S}
8  C u0 p0 c0 {7,D} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.659887,0.0848568,-5.09517e-05,1.45505e-08,-1.49836e-12,10194,39.9395], Tmin=(298,'K'), Tmax=(1384,'K')),
            NASAPolynomial(coeffs=[23.8138,0.0329112,-1.11043e-05,1.70913e-09,-9.85919e-14,1217.67,-93.2788], Tmin=(1384,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCC[CH]CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 234,
    label = "S3XC9H17",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {2,S} {7,S} {18,S} {19,S}
5  C u0 p0 c0 {1,S} {8,S} {14,S} {15,S}
6  C u0 p0 c0 {2,S} {20,S} {21,S} {22,S}
7  C u1 p0 c0 {3,S} {4,S} {23,S}
8  C u0 p0 c0 {5,S} {9,D} {24,S}
9  C u0 p0 c0 {8,D} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.953199,0.097325,-5.94746e-05,1.75306e-08,-1.92049e-12,7315.75,42.9874], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[27.0832,0.0372236,-1.25652e-05,1.93456e-09,-1.11618e-13,-2909.67,-109.416], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCC[CH]CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 235,
    label = "S4XC10H19",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {8,S} {19,S} {20,S}
5  C u0 p0 c0 {1,S} {8,S} {21,S} {22,S}
6  C u0 p0 c0 {2,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
8  C u1 p0 c0 {4,S} {5,S} {26,S}
9  C u0 p0 c0 {6,S} {10,D} {27,S}
10 C u0 p0 c0 {9,D} {28,S} {29,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.14353,0.109291,-6.72748e-05,2.00666e-08,-2.24431e-12,4423.59,45.57], Tmin=(298,'K'), Tmax=(1385,'K')),
            NASAPolynomial(coeffs=[30.3374,0.0415135,-1.40112e-05,2.15696e-09,-1.2444e-13,-7026.53,-125.45], Tmin=(1385,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCC[CH]CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 236,
    label = "C2H5-2-1C6H11",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {6,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {2,S} {20,S} {21,S} {22,S}
6  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
7  C u0 p0 c0 {3,S} {4,S} {8,D}
8  C u0 p0 c0 {7,D} {23,S} {24,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.77723,0.0955634,-5.59149e-05,9.8885e-09,1.71932e-12,-14686.4,44.5293], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[20.4993,0.0387739,-1.26042e-05,1.94964e-09,-1.17693e-13,-21653,-78.02], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(CC)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 237,
    label = "C4H9-2-1C6H11",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {5,S} {11,S} {12,S}
2  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
3  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {8,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {9,S} {19,S} {20,S}
6  C u0 p0 c0 {2,S} {9,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {26,S} {27,S} {28,S}
9  C u0 p0 c0 {5,S} {6,S} {10,D}
10 C u0 p0 c0 {9,D} {29,S} {30,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {10,S}
30 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.55941,0.121061,-8.27599e-05,2.96444e-08,-4.40185e-12,-20571.1,48.1589], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[30.8677,0.0442388,-1.50739e-05,2.33436e-09,-1.35197e-13,-32354.1,-131.894], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(CCCC)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 238,
    label = "PXC3H6-3-1C7H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
5  C u0 p0 c0 {4,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {9,S} {20,S} {21,S}
7  C u0 p0 c0 {5,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {1,S} {10,D} {25,S}
9  C u1 p0 c0 {6,S} {26,S} {27,S}
10 C u0 p0 c0 {8,D} {28,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.26798,0.117503,-8.10012e-05,2.92025e-08,-4.35371e-12,4830.2,48.889], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[30.3178,0.0422783,-1.44231e-05,2.23533e-09,-1.29532e-13,-6610.76,-126.492], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC(C=C)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 239,
    label = "C3H7-2-1C6H11",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
3  C u0 p0 c0 {4,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {3,S} {8,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {8,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {2,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {5,S} {9,D}
9  C u0 p0 c0 {8,D} {26,S} {27,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.00968,0.107624,-6.36488e-05,1.22993e-08,1.42032e-12,-17578.7,47.2893], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[23.3221,0.0434132,-1.40707e-05,2.17304e-09,-1.3109e-13,-25489,-91.4189], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(CCC)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 240,
    label = "C2H5S3XcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {5,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {8,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {5,S} {6,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.88541,0.0992968,-5.39692e-05,1.06981e-08,3.38822e-14,743.536,66.7862], Tmin=(298,'K'), Tmax=(1371,'K')),
            NASAPolynomial(coeffs=[24.5833,0.0352272,-1.21865e-05,1.90774e-09,-1.11364e-13,-11543.4,-111.494], Tmin=(1371,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1C[CH]CCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 241,
    label = "C2H5S4XcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {6,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {2,S} {8,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {5,S} {6,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.88541,0.0992968,-5.39692e-05,1.06981e-08,3.38822e-14,743.536,66.7862], Tmin=(298,'K'), Tmax=(1371,'K')),
            NASAPolynomial(coeffs=[24.5833,0.0352272,-1.21865e-05,1.90774e-09,-1.11364e-13,-11543.4,-111.494], Tmin=(1371,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1CC[CH]CC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 242,
    label = "S2XC4H8cC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {6,S} {20,S} {21,S}
4  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {6,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {5,S} {18,S} {19,S}
7  C u0 p0 c0 {1,S} {10,S} {24,S} {25,S}
8  C u0 p0 c0 {9,S} {10,S} {22,S} {23,S}
9  C u0 p0 c0 {8,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {7,S} {8,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.43326,0.12404,-7.07409e-05,1.64921e-08,-7.74026e-13,-5018.04,72.7079], Tmin=(298,'K'), Tmax=(1376,'K')),
            NASAPolynomial(coeffs=[31.1137,0.0438458,-1.51034e-05,2.35758e-09,-1.37346e-13,-19793.2,-143.715], Tmin=(1376,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]CC1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 243,
    label = "C4H9S3XcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {16,S} {17,S}
3  C u0 p0 c0 {1,S} {5,S} {20,S} {21,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {8,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {9,S} {12,S} {13,S}
7  C u0 p0 c0 {1,S} {10,S} {22,S} {23,S}
8  C u0 p0 c0 {5,S} {10,S} {24,S} {25,S}
9  C u0 p0 c0 {6,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {7,S} {8,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.43326,0.12404,-7.07409e-05,1.64921e-08,-7.74026e-13,-5018.04,72.7079], Tmin=(298,'K'), Tmax=(1376,'K')),
            NASAPolynomial(coeffs=[31.1137,0.0438458,-1.51034e-05,2.35758e-09,-1.37346e-13,-19793.2,-143.715], Tmin=(1376,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC1C[CH]CCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 244,
    label = "CH3S4XcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {17,S} {18,S} {19,S}
7  C u1 p0 c0 {4,S} {5,S} {20,S}
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
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.5397,0.0919152,-5.4639e-05,1.44889e-08,-1.27842e-12,4149.49,67.2785], Tmin=(298,'K'), Tmax=(2030,'K')),
            NASAPolynomial(coeffs=[17.3197,0.0361145,-1.30439e-05,2.10911e-09,-1.26099e-13,-4744.32,-72.4615], Tmin=(2030,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1CC[CH]CC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 245,
    label = "SXC2H4cC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {18,S} {19,S}
4  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {5,S} {16,S} {17,S}
7  C u0 p0 c0 {8,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {1,S} {7,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.88541,0.0992968,-5.39692e-05,1.06981e-08,3.38822e-14,743.536,66.7862], Tmin=(298,'K'), Tmax=(1371,'K')),
            NASAPolynomial(coeffs=[24.5833,0.0352272,-1.21865e-05,1.90774e-09,-1.11364e-13,-11543.4,-111.494], Tmin=(1371,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]C1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 246,
    label = "C3H7S3XcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {4,S} {17,S} {18,S}
4  C u0 p0 c0 {3,S} {7,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {8,S} {11,S} {12,S}
6  C u0 p0 c0 {1,S} {9,S} {19,S} {20,S}
7  C u0 p0 c0 {4,S} {9,S} {21,S} {22,S}
8  C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
9  C u1 p0 c0 {6,S} {7,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.16097,0.111676,-6.2366e-05,1.36019e-08,-3.71583e-13,-2137,69.7547], Tmin=(298,'K'), Tmax=(1374,'K')),
            NASAPolynomial(coeffs=[27.8484,0.0395366,-1.3645e-05,2.13266e-09,-1.24355e-13,-15668.3,-127.604], Tmin=(1374,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC1C[CH]CCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 247,
    label = "PXC3H6cC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {7,S} {21,S} {22,S}
4  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {5,S} {7,S} {17,S} {18,S}
7  C u0 p0 c0 {3,S} {6,S} {19,S} {20,S}
8  C u0 p0 c0 {4,S} {9,S} {23,S} {24,S}
9  C u1 p0 c0 {8,S} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {3,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.61025,0.116809,-7.15343e-05,1.98377e-08,-1.82506e-12,-828.177,70.2183], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[28.0463,0.0397796,-1.38109e-05,2.16643e-09,-1.26619e-13,-14311.6,-129.526], Tmin=(1383,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 248,
    label = "SXC3H6cC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {19,S} {20,S}
4  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {5,S} {17,S} {18,S}
7  C u0 p0 c0 {1,S} {9,S} {21,S} {22,S}
8  C u0 p0 c0 {9,S} {23,S} {24,S} {25,S}
9  C u1 p0 c0 {7,S} {8,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.16097,0.111676,-6.2366e-05,1.36019e-08,-3.71583e-13,-2137,69.7547], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[27.8484,0.0395366,-1.3645e-05,2.13266e-09,-1.24355e-13,-15668.3,-127.604], Tmin=(1383,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CC1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 249,
    label = "SXC4H8cC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
2  C u0 p0 c0 {1,S} {5,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {7,S} {22,S} {23,S}
4  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {6,S} {16,S} {17,S}
6  C u0 p0 c0 {5,S} {7,S} {18,S} {19,S}
7  C u0 p0 c0 {3,S} {6,S} {20,S} {21,S}
8  C u0 p0 c0 {4,S} {10,S} {24,S} {25,S}
9  C u0 p0 c0 {10,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {8,S} {9,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {3,S}
23 H u0 p0 c0 {3,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.43326,0.12404,-7.07409e-05,1.64921e-08,-7.74026e-13,-5018.04,72.7079], Tmin=(298,'K'), Tmax=(1376,'K')),
            NASAPolynomial(coeffs=[31.1137,0.0438458,-1.51034e-05,2.35758e-09,-1.37346e-13,-19793.2,-143.715], Tmin=(1376,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCC1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 250,
    label = "S3XC4H8cC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {6,S} {22,S} {23,S}
4  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
5  C u0 p0 c0 {4,S} {6,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {5,S} {20,S} {21,S}
7  C u0 p0 c0 {8,S} {9,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {10,S} {24,S} {25,S}
9  C u0 p0 c0 {7,S} {26,S} {27,S} {28,S}
10 C u1 p0 c0 {1,S} {8,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {3,S}
23 H u0 p0 c0 {3,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.43326,0.12404,-7.07409e-05,1.64921e-08,-7.74026e-13,-5018.04,72.7079], Tmin=(298,'K'), Tmax=(1376,'K')),
            NASAPolynomial(coeffs=[31.1137,0.0438458,-1.51034e-05,2.35758e-09,-1.37346e-13,-19793.2,-143.715], Tmin=(1376,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[CH]C1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 251,
    label = "PXC4H8cC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {11,S}
2  C u0 p0 c0 {1,S} {8,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {5,S} {16,S} {17,S}
4  C u0 p0 c0 {1,S} {7,S} {24,S} {25,S}
5  C u0 p0 c0 {3,S} {6,S} {18,S} {19,S}
6  C u0 p0 c0 {5,S} {7,S} {20,S} {21,S}
7  C u0 p0 c0 {4,S} {6,S} {22,S} {23,S}
8  C u0 p0 c0 {2,S} {9,S} {12,S} {13,S}
9  C u0 p0 c0 {8,S} {10,S} {26,S} {27,S}
10 C u1 p0 c0 {9,S} {28,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {8,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {4,S}
25 H u0 p0 c0 {4,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.88803,0.12921,-7.99739e-05,2.27722e-08,-2.23752e-12,-3708.69,73.1945], Tmin=(298,'K'), Tmax=(1384,'K')),
            NASAPolynomial(coeffs=[31.2925,0.0441161,-1.52806e-05,2.39325e-09,-1.39725e-13,-18427.1,-145.526], Tmin=(1384,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCC1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 252,
    label = "PX10-4C10H19",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
3  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {8,S} {19,S} {20,S}
5  C u0 p0 c0 {3,S} {9,S} {21,S} {22,S}
6  C u0 p0 c0 {1,S} {10,S} {17,S} {18,S}
7  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {9,D} {26,S}
9  C u0 p0 c0 {5,S} {8,D} {27,S}
10 C u1 p0 c0 {6,S} {28,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.08366,0.115238,-7.7328e-05,2.68948e-08,-3.84652e-12,4455.24,49.4927], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[30.0352,0.042353,-1.44153e-05,2.23087e-09,-1.29149e-13,-6939.47,-123.81], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCC=CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 253,
    label = "PXC2H4cC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {18,S} {19,S}
4  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {5,S} {16,S} {17,S}
7  C u0 p0 c0 {1,S} {8,S} {20,S} {21,S}
8  C u1 p0 c0 {7,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.34951,0.104479,-6.3186e-05,1.69547e-08,-1.4235e-12,2054.98,67.3216], Tmin=(298,'K'), Tmax=(1382,'K')),
            NASAPolynomial(coeffs=[24.8102,0.0354349,-1.23386e-05,1.93923e-09,-1.13492e-13,-10201.8,-113.587], Tmin=(1382,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 254,
    label = "CH3-5-SAX1C9H16",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {4,S} {14,S} {15,S}
4  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {8,S} {18,S} {19,S}
6  C u0 p0 c0 {1,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {4,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {5,S} {9,D} {26,S}
9  C u0 p0 c0 {8,D} {10,S} {27,S}
10 C u1 p0 c0 {9,S} {28,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.98842,0.118472,-7.995e-05,2.7918e-08,-4.02451e-12,-2528.89,49.4378], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[30.6358,0.0426181,-1.46746e-05,2.28852e-09,-1.33192e-13,-14542.2,-132.213], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CCC(C)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 255,
    label = "CH3-5-SAX1C8H14",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {4,S} {8,D} {23,S}
8  C u0 p0 c0 {7,D} {9,S} {24,S}
9  C u1 p0 c0 {8,S} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.72628,0.106086,-7.14846e-05,2.49338e-08,-3.59431e-12,355.531,46.5458], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[27.3935,0.0382539,-1.31903e-05,2.05892e-09,-1.19904e-13,-10425.2,-116.227], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CCC(C)CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 256,
    label = "C4H9TXcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {17,S} {18,S}
2  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {8,S} {19,S} {20,S}
5  C u0 p0 c0 {2,S} {9,S} {11,S} {12,S}
6  C u0 p0 c0 {2,S} {10,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {10,S} {23,S} {24,S}
8  C u0 p0 c0 {4,S} {10,S} {25,S} {26,S}
9  C u0 p0 c0 {5,S} {27,S} {28,S} {29,S}
10 C u1 p0 c0 {6,S} {7,S} {8,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {1,S}
18 H u0 p0 c0 {1,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.57636,0.111756,-5.48518e-05,8.07611e-09,8.35132e-13,-6429.69,59.8607], Tmin=(298,'K'), Tmax=(1375,'K')),
            NASAPolynomial(coeffs=[29.7531,0.0454307,-1.57353e-05,2.46433e-09,-1.43867e-13,-20228.8,-135.481], Tmin=(1375,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC[C]1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 257,
    label = "C2H5S2XcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {16,S} {17,S}
3  C u0 p0 c0 {2,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
5  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
6  C u0 p0 c0 {5,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
8  C u1 p0 c0 {1,S} {6,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.88541,0.0992968,-5.39692e-05,1.06981e-08,3.38822e-14,743.536,66.7862], Tmin=(298,'K'), Tmax=(1371,'K')),
            NASAPolynomial(coeffs=[24.5833,0.0352272,-1.21865e-05,1.90774e-09,-1.11364e-13,-11543.4,-111.494], Tmin=(1371,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1[CH]CCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 258,
    label = "S2XC3H6cC6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {19,S} {20,S}
4  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
5  C u0 p0 c0 {4,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {5,S} {17,S} {18,S}
7  C u0 p0 c0 {8,S} {9,S} {21,S} {22,S}
8  C u0 p0 c0 {7,S} {23,S} {24,S} {25,S}
9  C u1 p0 c0 {1,S} {7,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.16097,0.111676,-6.2366e-05,1.36019e-08,-3.71583e-13,-2137,69.7547], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[27.8484,0.0395366,-1.3645e-05,2.13266e-09,-1.24355e-13,-15668.3,-127.604], Tmin=(1383,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[CH]C1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 259,
    label = "C3H7TXcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {14,S} {15,S}
2  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {5,S} {8,S} {10,S} {11,S}
5  C u0 p0 c0 {4,S} {9,S} {18,S} {19,S}
6  C u0 p0 c0 {2,S} {9,S} {20,S} {21,S}
7  C u0 p0 c0 {3,S} {9,S} {22,S} {23,S}
8  C u0 p0 c0 {4,S} {24,S} {25,S} {26,S}
9  C u1 p0 c0 {5,S} {6,S} {7,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.30898,0.0993922,-4.64546e-05,5.16496e-09,1.24251e-12,-3547.4,56.9343], Tmin=(298,'K'), Tmax=(1372,'K')),
            NASAPolynomial(coeffs=[26.5101,0.041096,-1.42672e-05,2.23786e-09,-1.30785e-13,-16116.4,-119.504], Tmin=(1372,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[C]1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 260,
    label = "PX9-3C9H17",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {1,S} {9,S} {14,S} {15,S}
5  C u0 p0 c0 {6,S} {8,S} {18,S} {19,S}
6  C u0 p0 c0 {5,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {8,D} {23,S}
8  C u0 p0 c0 {5,S} {7,D} {24,S}
9  C u1 p0 c0 {4,S} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.81323,0.102878,-6.89568e-05,2.40055e-08,-3.44422e-12,7336.64,46.5487], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[26.7765,0.0380351,-1.29535e-05,2.0054e-09,-1.16124e-13,-2817.92,-107.738], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCCC=CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 261,
    label = "SAXC6H9-15",
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
            NASAPolynomial(coeffs=[0.14774,0.056651,-3.64346e-05,1.1938e-08,-1.71764e-12,23589.7,28.0327], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[15.2564,0.0203886,-5.84222e-06,8.27671e-10,-4.71388e-14,18858.3,-52.0823], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/28/  THERGA""",
    longDesc = 
u"""
7/28/  THERGA
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[CH]CC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 262,
    label = "CH2cC6H10",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {6,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {18,S} {19,S}
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
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.78092,0.0929268,-5.9179e-05,1.73251e-08,-1.77259e-12,-4903.05,64.9507], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[20.9535,0.0293619,-1.02414e-05,1.61157e-09,-9.43983e-14,-15713.1,-96.6591], Tmin=(1383,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 263,
    label = "C2H5TXcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {8,S} {17,S} {18,S}
5  C u0 p0 c0 {3,S} {8,S} {19,S} {20,S}
6  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
7  C u0 p0 c0 {6,S} {21,S} {22,S} {23,S}
8  C u1 p0 c0 {4,S} {5,S} {6,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.76244,0.109462,-7.54295e-05,2.6916e-08,-4.20357e-12,1196.23,61.5953], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[19.1764,0.0418358,-1.41184e-05,2.24297e-09,-1.3794e-13,-7011.23,-80.3126], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC[C]1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 264,
    label = "C2H5-3-PXC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {6,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {1,S} {8,D} {19,S}
7  C u1 p0 c0 {4,S} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {22,S} {23,S}
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
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.47916,0.092823,-6.50277e-05,2.40441e-08,-3.68237e-12,9838.52,42.2096], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[23.7873,0.0334277,-1.12875e-05,1.73745e-09,-1.00205e-13,1088.69,-93.3753], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC(C=C)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 265,
    label = "PAXCH2-2-1C6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {6,S} {7,D}
6  C u1 p0 c0 {5,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.6261,0.0840492,-6.12142e-05,2.37212e-08,-3.84407e-12,6142.7,37.9624], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[21.7766,0.0286182,-9.87112e-06,1.54112e-09,-8.97595e-14,-1997.18,-87.6034], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=C)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 266,
    label = "CH3-2-PXC4H6",
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
            NASAPolynomial(coeffs=[0.516043,0.0506942,-3.28148e-05,1.12296e-08,-1.62047e-12,18616.3,27.6661], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[14.1325,0.0204049,-6.97672e-06,1.08277e-09,-6.28003e-14,13673,-46.1058], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 267,
    label = "CH3-3-TAXC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {2,S} {4,S} {6,S}
6  C u0 p0 c0 {5,S} {7,D} {18,S}
7  C u0 p0 c0 {6,D} {19,S} {20,S}
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
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.11204,0.0699305,-3.27762e-05,4.29403e-09,5.4659e-13,3210.07,42.2783], Tmin=(298,'K'), Tmax=(1364,'K')),
            NASAPolynomial(coeffs=[22.964,0.0288693,-1.02474e-05,1.6309e-09,-9.62704e-14,-6214.13,-86.4295], Tmin=(1364,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[C](C)CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 268,
    label = "C2H5-4-SAX1C8H14",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {14,S} {15,S}
3  C u0 p0 c0 {2,S} {5,S} {16,S} {17,S}
4  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {3,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {1,S} {9,D} {26,S}
9  C u0 p0 c0 {8,D} {10,S} {27,S}
10 C u1 p0 c0 {9,S} {28,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.98842,0.118472,-7.995e-05,2.7918e-08,-4.02451e-12,-2931.51,49.4378], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[30.6358,0.0426181,-1.46746e-05,2.28852e-09,-1.33192e-13,-14944.8,-132.213], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CC(CC)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 269,
    label = "C2H5-4-SAX1C7H12",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {1,S} {8,D} {23,S}
8  C u0 p0 c0 {7,D} {9,S} {24,S}
9  C u1 p0 c0 {8,S} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.72628,0.106086,-7.14846e-05,2.49338e-08,-3.59431e-12,-47.0865,46.5458], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[27.3935,0.0382539,-1.31903e-05,2.05892e-09,-1.19904e-13,-10827.8,-116.227], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C=CC(CC)CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 270,
    label = "PAXCH2-2-1C5H9",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {5,S} {6,D}
5  C u1 p0 c0 {4,S} {14,S} {15,S}
6  C u0 p0 c0 {4,D} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.29812,0.0713005,-5.21275e-05,2.03174e-08,-3.31747e-12,9018.87,34.7791], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[18.5996,0.0241598,-8.34657e-06,1.30451e-09,-7.60357e-14,2088.35,-71.9954], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=C)CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 271,
    label = "SAX5-3C8H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
3  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {8,D} {22,S}
7  C u1 p0 c0 {2,S} {8,S} {21,S}
8  C u0 p0 c0 {6,D} {7,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.35404,0.0921118,-6.0654e-05,2.05583e-08,-2.87162e-12,2590.03,44.5285], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[23.7926,0.0340643,-1.17383e-05,1.83155e-09,-1.06634e-13,-6859.31,-97.0816], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC=C[CH]CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 272,
    label = "SAX4-2C7H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {6,S} {15,S} {16,S} {17,S}
5  C u1 p0 c0 {2,S} {7,S} {18,S}
6  C u0 p0 c0 {4,S} {7,D} {19,S}
7  C u0 p0 c0 {5,S} {6,D} {20,S}
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
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.13711,0.0758043,-4.6837e-05,1.44624e-08,-1.786e-12,5237.06,36.9626], Tmin=(298,'K'), Tmax=(1387,'K')),
            NASAPolynomial(coeffs=[20.4542,0.0297102,-1.02421e-05,1.59852e-09,-9.30844e-14,-2757.27,-80.6405], Tmin=(1387,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC=C[CH]CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 273,
    label = "PAXCH2-2-1C4H7",
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
            NASAPolynomial(coeffs=[-0.945638,0.0585063,-4.30091e-05,1.69152e-08,-2.79442e-12,12292.1,31.4683], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[15.3933,0.0197642,-6.85172e-06,1.07333e-09,-6.26601e-14,6588.87,-56.2239], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=C)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 274,
    label = "PXCH2-3-1C7H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {9,S}
2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  C u0 p0 c0 {2,S} {4,S} {12,S} {13,S}
4  C u0 p0 c0 {3,S} {5,S} {14,S} {15,S}
5  C u0 p0 c0 {4,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {1,S} {8,D} {19,S}
7  C u1 p0 c0 {1,S} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {22,S} {23,S}
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
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.58551,0.0919236,-6.27829e-05,2.24798e-08,-3.34343e-12,10854.6,42.3931], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[23.849,0.033629,-1.15018e-05,1.7856e-09,-1.03592e-13,1859.76,-94.6861], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C=C)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 275,
    label = "PX1-3C6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {5,S} {6,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {5,D} {15,S}
5  C u0 p0 c0 {2,S} {4,D} {14,S}
6  C u1 p0 c0 {2,S} {16,S} {17,S}
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
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.977624,0.0656665,-4.36212e-05,1.51942e-08,-2.20544e-12,15977.8,37.6094], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[17.0133,0.0250666,-8.56224e-06,1.32807e-09,-7.69999e-14,9538.6,-59.6004], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC=CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 276,
    label = "PXCH2-3-1C6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {7,D} {16,S}
6  C u1 p0 c0 {1,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
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
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.0510707,0.0729893,-4.0495e-05,1.00132e-08,-8.10665e-13,11990.9,43.9029], Tmin=(298,'K'), Tmax=(1372,'K')),
            NASAPolynomial(coeffs=[22.7491,0.0287965,-1.01622e-05,1.61105e-09,-9.48402e-14,3077.27,-81.3859], Tmin=(1372,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C=C)CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 277,
    label = "S4XC8H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {6,S} {7,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
6  C u1 p0 c0 {3,S} {4,S} {20,S}
7  C u0 p0 c0 {4,S} {8,D} {21,S}
8  C u0 p0 c0 {7,D} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.659887,0.0848568,-5.09517e-05,1.45505e-08,-1.49836e-12,10194,39.9395], Tmin=(298,'K'), Tmax=(1384,'K')),
            NASAPolynomial(coeffs=[23.8138,0.0329112,-1.11043e-05,1.70913e-09,-9.85919e-14,1217.67,-93.2788], Tmin=(1384,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC[CH]CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 278,
    label = "PXCH2-2-C6H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {6,S} {13,S} {14,S}
5  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6  C u0 p0 c0 {4,S} {18,S} {19,S} {20,S}
7  C u1 p0 c0 {1,S} {21,S} {22,S}
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
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.11073,0.0833825,-5.34421e-05,1.7763e-08,-2.44822e-12,-1270.89,39.0662], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[21.7823,0.0329533,-1.12932e-05,1.75545e-09,-1.0193e-13,-9602.9,-85.093], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 279,
    label = "PXCH2-2-C4H9",
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
            NASAPolynomial(coeffs=[-0.592409,0.0586735,-3.65981e-05,1.18413e-08,-1.59846e-12,4094.94,33.3003], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[15.3819,0.0241141,-8.27848e-06,1.28832e-09,-7.48649e-14,-1807.73,-53.602], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 280,
    label = "PX10-5C10H19",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {8,S} {19,S} {20,S}
5  C u0 p0 c0 {1,S} {9,S} {21,S} {22,S}
6  C u0 p0 c0 {2,S} {10,S} {17,S} {18,S}
7  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {9,D} {26,S}
9  C u0 p0 c0 {5,S} {8,D} {27,S}
10 C u1 p0 c0 {6,S} {28,S} {29,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.78077,0.113827,-7.52241e-05,2.56837e-08,-3.61636e-12,4411.35,48.1043], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[30.2162,0.042342,-1.45334e-05,2.26171e-09,-1.31437e-13,-7095.16,-125.011], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCC=CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 281,
    label = "PXC4H8-2-1C6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
4  C u0 p0 c0 {2,S} {8,S} {19,S} {20,S}
5  C u0 p0 c0 {1,S} {8,S} {21,S} {22,S}
6  C u0 p0 c0 {2,S} {9,S} {17,S} {18,S}
7  C u0 p0 c0 {3,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {4,S} {5,S} {10,D}
9  C u1 p0 c0 {6,S} {26,S} {27,S}
10 C u0 p0 c0 {8,D} {28,S} {29,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.79074,0.116509,-8.02793e-05,2.90054e-08,-4.34419e-12,4041.66,47.0399], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[30.4469,0.0420992,-1.4347e-05,2.22204e-09,-1.28704e-13,-7287.44,-126.488], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCC(=C)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 282,
    label = "C2H5-2-PXC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {5,S} {6,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {2,S} {3,S} {8,D}
7  C u1 p0 c0 {4,S} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.47916,0.092823,-6.50277e-05,2.40441e-08,-3.68237e-12,9838.52,42.2096], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[23.7873,0.0334277,-1.12875e-05,1.73745e-09,-1.00205e-13,1088.69,-93.3753], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCC(=C)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 283,
    label = "CH2*",
    molecule = 
"""
multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.1986,-0.00236661,8.23296e-06,-6.68816e-09,1.94315e-12,50496.8,-0.769119], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.29204,0.00465589,-2.01192e-06,4.17906e-10,-3.39716e-14,50926,8.6265], Tmin=(1000,'K'), Tmax=(3500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3500,'K'),
    ),
    shortDesc = u"""L S/93""",
    longDesc = 
u"""
L S/93.
Duplicate of species CH2 (i.e. same molecular structure according to RMG)
[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 284,
    label = "PXC2H4-4-1C8H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {11,S}
2  C u0 p0 c0 {1,S} {3,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {4,S} {14,S} {15,S}
4  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {8,S} {20,S} {21,S}
6  C u0 p0 c0 {1,S} {9,S} {18,S} {19,S}
7  C u0 p0 c0 {4,S} {22,S} {23,S} {24,S}
8  C u0 p0 c0 {5,S} {10,D} {25,S}
9  C u1 p0 c0 {6,S} {26,S} {27,S}
10 C u0 p0 c0 {8,D} {28,S} {29,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.18648,0.11691,-7.99349e-05,2.85126e-08,-4.20287e-12,4698.3,48.5598], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[30.3291,0.0423114,-1.44442e-05,2.23965e-09,-1.29826e-13,-6767.44,-126.614], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(CC=C)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 285,
    label = "PXC2H4-4-1C7H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {10,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {7,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {8,S} {15,S} {16,S}
6  C u0 p0 c0 {3,S} {19,S} {20,S} {21,S}
7  C u0 p0 c0 {4,S} {9,D} {22,S}
8  C u1 p0 c0 {5,S} {23,S} {24,S}
9  C u0 p0 c0 {7,D} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.86646,0.104331,-7.12367e-05,2.54237e-08,-3.75769e-12,7572.17,45.386], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[27.0848,0.0379747,-1.29747e-05,2.0129e-09,-1.16725e-13,-2653.19,-110.625], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(CC=C)CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 286,
    label = "CH3-4-PXC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {4,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {7,D} {16,S}
6  C u1 p0 c0 {3,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.94472,0.0871524,-7.24735e-05,3.62361e-08,-8.46811e-12,13746,38.9027], Tmin=(298,'K'), Tmax=(1372,'K')),
            NASAPolynomial(coeffs=[19.8708,0.0274743,-7.90331e-06,1.12063e-09,-6.37471e-14,7344.05,-74.7968], Tmin=(1372,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(C)CC=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 287,
    label = "CH3S3XcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {11,S} {12,S}
3  C u0 p0 c0 {2,S} {5,S} {9,S} {10,S}
4  C u0 p0 c0 {1,S} {7,S} {13,S} {14,S}
5  C u0 p0 c0 {3,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {17,S} {18,S} {19,S}
7  C u1 p0 c0 {4,S} {5,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.5397,0.0919152,-5.4639e-05,1.44889e-08,-1.27842e-12,4149.49,67.973], Tmin=(298,'K'), Tmax=(2030,'K')),
            NASAPolynomial(coeffs=[17.3197,0.0361145,-1.30439e-05,2.10911e-09,-1.26099e-13,-4744.32,-71.767], Tmin=(2030,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1C[CH]CCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 288,
    label = "PX7-2C7H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
4  C u0 p0 c0 {6,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {6,D} {17,S}
6  C u0 p0 c0 {4,S} {5,D} {18,S}
7  C u1 p0 c0 {3,S} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.65701,0.080941,-6.47455e-05,3.26716e-08,-7.96074e-12,12720.2,33.6426], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[19.6992,0.0273057,-7.79661e-06,1.1004e-09,-6.24402e-14,6511.8,-73.2632], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCC=CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 289,
    label = "CH3S2XcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {13,S} {14,S}
3  C u0 p0 c0 {2,S} {4,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {7,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {17,S} {18,S} {19,S}
7  C u1 p0 c0 {1,S} {5,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-8.5397,0.0919152,-5.4639e-05,1.44889e-08,-1.27842e-12,4149.49,67.973], Tmin=(298,'K'), Tmax=(2030,'K')),
            NASAPolynomial(coeffs=[17.3197,0.0361145,-1.30439e-05,2.10911e-09,-1.26099e-13,-4744.32,-71.767], Tmin=(2030,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1[CH]CCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 290,
    label = "CH3-3-PXC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {7,D} {16,S}
6  C u1 p0 c0 {3,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
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
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.925345,0.087751,-7.25544e-05,3.54518e-08,-8.02218e-12,13392.3,33.4286], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[20.4343,0.0287095,-8.63114e-06,1.26195e-09,-7.32388e-14,7203.41,-77.6291], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC(C)C=C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 291,
    label = "CH3TXcC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {4,S} {8,S} {9,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {3,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {7,S} {18,S} {19,S} {20,S}
7  C u1 p0 c0 {4,S} {5,S} {6,S}
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
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-7.21964,0.0957038,-6.71531e-05,2.66733e-08,-5.32466e-12,4159.66,56.9924], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[18.2642,0.0339358,-1.11432e-05,1.73561e-09,-1.05251e-13,-3919.62,-78.273], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[C]1CCCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 292,
    label = "SAXcC6H9",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
4  C u1 p0 c0 {3,S} {6,S} {13,S}
5  C u0 p0 c0 {2,S} {6,D} {14,S}
6  C u0 p0 c0 {4,S} {5,D} {15,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.86858,0.0738153,-5.33135e-05,1.95943e-08,-2.95965e-12,13547.7,51.457], Tmin=(298,'K'), Tmax=(1381,'K')),
            NASAPolynomial(coeffs=[17.1128,0.0215146,-7.60952e-06,1.20824e-09,-7.12064e-14,5375.52,-72.5512], Tmin=(1381,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/28/  CBSQB3""",
    longDesc = 
u"""
7/28/  CBSQB3
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]1C=CCCC1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 293,
    label = "SXC7H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {2,S} {4,S} {17,S}
6  C u0 p0 c0 {3,S} {7,D} {18,S}
7  C u0 p0 c0 {6,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.337871,0.0722549,-4.22499e-05,1.14646e-08,-1.05307e-12,13068.2,36.7608], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[20.5321,0.0285989,-9.64126e-06,1.48314e-09,-8.55246e-14,5352.6,-77.0636], Tmin=(1383,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CCCC[CH]C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 294,
    label = "CH3-2-SXC6H12",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {5,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u0 p0 c0 {2,S} {7,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {16,S} {17,S} {18,S}
6  C u0 p0 c0 {7,S} {19,S} {20,S} {21,S}
7  C u1 p0 c0 {3,S} {6,S} {22,S}
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
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.656082,0.0782944,-4.44575e-05,1.16721e-08,-1.02741e-12,-2581.97,38.5655], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[21.4813,0.0328091,-1.11611e-05,1.72676e-09,-9.99478e-14,-10899.9,-82.5476], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C[CH]CCC(C)C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 295,
    label = "PXC3H6-3-1C6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {8,S} {17,S} {18,S}
6  C u0 p0 c0 {4,S} {19,S} {20,S} {21,S}
7  C u0 p0 c0 {1,S} {9,D} {22,S}
8  C u1 p0 c0 {5,S} {23,S} {24,S}
9  C u0 p0 c0 {7,D} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.06986,0.101314,-6.36018e-05,2.02116e-08,-2.66583e-12,5886.16,52.2665], Tmin=(298,'K'), Tmax=(1379,'K')),
            NASAPolynomial(coeffs=[29.055,0.0378577,-1.3283e-05,2.09769e-09,-1.23156e-13,-5464.61,-112.312], Tmin=(1379,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCC(C=C)CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 296,
    label = "PX9-4C9H17",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {2,S} {8,S} {18,S} {19,S}
5  C u0 p0 c0 {1,S} {9,S} {14,S} {15,S}
6  C u0 p0 c0 {2,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {8,D} {23,S}
8  C u0 p0 c0 {4,S} {7,D} {24,S}
9  C u1 p0 c0 {5,S} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.81323,0.102878,-6.89568e-05,2.40055e-08,-3.44422e-12,7336.64,46.5487], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[26.7765,0.0380351,-1.29535e-05,2.0054e-09,-1.16124e-13,-2817.92,-107.738], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCC=CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 297,
    label = "C8H14-13",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {6,D} {18,S}
6  C u0 p0 c0 {5,D} {7,S} {20,S}
7  C u0 p0 c0 {6,S} {8,D} {19,S}
8  C u0 p0 c0 {7,D} {21,S} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {8,S}
22 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.00851,0.098151,-7.59831e-05,3.1368e-08,-5.34847e-12,-1046.93,45.0345], Tmin=(298,'K'), Tmax=(1395,'K')),
            NASAPolynomial(coeffs=[23.7241,0.0315576,-1.0775e-05,1.67082e-09,-9.68539e-14,-9982.5,-97.1766], Tmin=(1395,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC=CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 298,
    label = "CH3-2-PXC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {4,S} {7,D}
6  C u1 p0 c0 {3,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.114306,0.0747239,-4.84643e-05,1.63867e-08,-2.31984e-12,12431.6,32.9237], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[20.922,0.0289166,-9.97311e-06,1.55692e-09,-9.06736e-14,4830.04,-79.9799], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCC(=C)C
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 299,
    label = "SAXcC6H7",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,D} {9,S}
3  C u0 p0 c0 {1,S} {5,D} {10,S}
4  C u0 p0 c0 {2,D} {6,S} {12,S}
5  C u0 p0 c0 {3,D} {6,S} {13,S}
6  C u1 p0 c0 {4,S} {5,S} {11,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.08823,0.0728412,-6.56192e-05,3.02959e-08,-5.58316e-12,26539.9,33.5253], Tmin=(298,'K'), Tmax=(1394,'K')),
            NASAPolynomial(coeffs=[17.8037,0.016554,-5.77071e-06,9.07479e-10,-5.31246e-14,19987.5,-76.1002], Tmin=(1394,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""7/31/ 9 THERM""",
    longDesc = 
u"""
7/31/ 9 THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH]1C=CCC=C1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 300,
    label = "S3XC7H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {2,S} {3,S} {17,S}
6  C u0 p0 c0 {3,S} {7,D} {18,S}
7  C u0 p0 c0 {6,D} {19,S} {20,S}
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
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.337871,0.0722549,-4.22499e-05,1.14646e-08,-1.05307e-12,13068.2,36.7608], Tmin=(298,'K'), Tmax=(1383,'K')),
            NASAPolynomial(coeffs=[20.5321,0.0285989,-9.64126e-06,1.48314e-09,-8.55246e-14,5352.6,-77.0636], Tmin=(1383,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC[CH]CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 301,
    label = "SAXC7H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u1 p0 c0 {3,S} {6,S} {17,S}
6  C u0 p0 c0 {5,S} {7,D} {18,S}
7  C u0 p0 c0 {6,D} {19,S} {20,S}
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
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.66946,0.0793202,-5.20232e-05,1.74804e-08,-2.40913e-12,6759.27,38.4731], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[20.9278,0.0292841,-1.009e-05,1.57426e-09,-9.16506e-14,-1412.96,-83.9476], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[CH]CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 302,
    label = "C3H7-2-PXC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {2,S} {7,S} {18,S} {19,S}
5  C u0 p0 c0 {1,S} {8,S} {14,S} {15,S}
6  C u0 p0 c0 {2,S} {20,S} {21,S} {22,S}
7  C u0 p0 c0 {3,S} {4,S} {9,D}
8  C u1 p0 c0 {5,S} {23,S} {24,S}
9  C u0 p0 c0 {7,D} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.51445,0.104122,-7.18643e-05,2.6088e-08,-3.93565e-12,6922.2,44.0691], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[27.192,0.0377765,-1.28832e-05,1.99625e-09,-1.15661e-13,-3167.81,-110.438], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCC(=C)CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 303,
    label = "PX1-3C8H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {7,S} {8,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {7,D} {21,S}
7  C u0 p0 c0 {4,S} {6,D} {20,S}
8  C u1 p0 c0 {4,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.29095,0.0907951,-6.12886e-05,2.14028e-08,-3.04869e-12,11530.4,41.2455], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[23.8416,0.0331279,-1.1134e-05,1.70873e-09,-9.83549e-14,2724.59,-94.0576], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC=CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 304,
    label = "PXC2H4-2-1C6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {6,S} {7,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {4,S} {8,D}
7  C u1 p0 c0 {4,S} {20,S} {21,S}
8  C u0 p0 c0 {6,D} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.24549,0.0917391,-6.34284e-05,2.31491e-08,-3.52127e-12,9804.59,41.1379], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[23.9347,0.0334623,-1.14237e-05,1.77127e-09,-1.02671e-13,950.665,-94.3797], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(=C)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 305,
    label = "PXC2H4-2-1C4H7",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
2  C u0 p0 c0 {4,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {2,S} {6,D}
5  C u1 p0 c0 {2,S} {14,S} {15,S}
6  C u0 p0 c0 {4,D} {16,S} {17,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.698169,0.0669899,-4.66161e-05,1.73167e-08,-2.70423e-12,15566.3,35.2195], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[17.4696,0.0247603,-8.47337e-06,1.31586e-09,-7.63531e-14,9169.05,-62.5388], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(=C)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 306,
    label = "C3H7-3-TAX1C7H13",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
3  C u0 p0 c0 {4,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {3,S} {8,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {8,S} {19,S} {20,S}
6  C u0 p0 c0 {3,S} {21,S} {22,S} {23,S}
7  C u0 p0 c0 {2,S} {24,S} {25,S} {26,S}
8  C u1 p0 c0 {4,S} {5,S} {9,S}
9  C u0 p0 c0 {8,S} {10,D} {27,S}
10 C u0 p0 c0 {9,D} {28,S} {29,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.25794,0.113906,-7.30951e-05,2.37879e-08,-3.13888e-12,-3980.52,46.4279], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[30.1,0.0430491,-1.48175e-05,2.3102e-09,-1.34427e-13,-15751.6,-129.134], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[C](CCC)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 307,
    label = "SAX5-3C9H17",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {10,S} {11,S}
2  C u0 p0 c0 {1,S} {5,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {8,S} {16,S} {17,S}
4  C u0 p0 c0 {6,S} {7,S} {14,S} {15,S}
5  C u0 p0 c0 {2,S} {21,S} {22,S} {23,S}
6  C u0 p0 c0 {4,S} {18,S} {19,S} {20,S}
7  C u0 p0 c0 {4,S} {9,D} {25,S}
8  C u1 p0 c0 {3,S} {9,S} {24,S}
9  C u0 p0 c0 {7,D} {8,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.57755,0.104264,-6.8742e-05,2.32847e-08,-3.23934e-12,-298.328,46.5613], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[27.0099,0.0384317,-1.32197e-05,2.06026e-09,-1.19852e-13,-10961.2,-113.606], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC=C[CH]CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 308,
    label = "SAX6-4C9H17",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {8,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
6  C u0 p0 c0 {2,S} {21,S} {22,S} {23,S}
7  C u1 p0 c0 {3,S} {9,S} {24,S}
8  C u0 p0 c0 {4,S} {9,D} {25,S}
9  C u0 p0 c0 {7,S} {8,D} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.57755,0.104264,-6.8742e-05,2.32847e-08,-3.23934e-12,-298.328,46.5613], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[27.0099,0.0384317,-1.32197e-05,2.06026e-09,-1.19852e-13,-10961.2,-113.606], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[CH]C=CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 309,
    label = "C3H7-3-TAX1C6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {4,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {1,S} {7,S} {14,S} {15,S}
4  C u0 p0 c0 {2,S} {7,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
6  C u0 p0 c0 {2,S} {21,S} {22,S} {23,S}
7  C u1 p0 c0 {3,S} {4,S} {8,S}
8  C u0 p0 c0 {7,S} {9,D} {24,S}
9  C u0 p0 c0 {8,D} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.06986,0.101314,-6.36018e-05,2.02116e-08,-2.66583e-12,5886.16,52.2665], Tmin=(298,'K'), Tmax=(1379,'K')),
            NASAPolynomial(coeffs=[29.055,0.0378577,-1.3283e-05,2.09769e-09,-1.23156e-13,-5464.61,-112.312], Tmin=(1379,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[C](CCC)CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 310,
    label = "PX8-3C8H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {8,S} {11,S} {12,S}
4  C u0 p0 c0 {5,S} {7,S} {15,S} {16,S}
5  C u0 p0 c0 {4,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {2,S} {7,D} {20,S}
7  C u0 p0 c0 {4,S} {6,D} {21,S}
8  C u1 p0 c0 {3,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.29095,0.0907951,-6.12886e-05,2.14028e-08,-3.04869e-12,11530.4,41.2455], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[23.8416,0.0331279,-1.1134e-05,1.70873e-09,-9.83549e-14,2724.59,-94.0576], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CCCC=CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 311,
    label = "cC6H8-13",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {6,D} {11,S}
4  C u0 p0 c0 {2,S} {5,D} {12,S}
5  C u0 p0 c0 {4,D} {6,S} {13,S}
6  C u0 p0 c0 {3,D} {5,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-9.28596,0.0790151,-5.35608e-05,1.65155e-08,-1.81087e-12,10377.5,64.4512], Tmin=(298,'K'), Tmax=(1374,'K')),
            NASAPolynomial(coeffs=[18.0508,0.0197853,-7.09477e-06,1.1369e-09,-6.74331e-14,544.318,-83.8058], Tmin=(1374,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""CBSQB3""",
    longDesc = 
u"""
CBSQB3
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C1C=CCCC=1
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 312,
    label = "C2H5-2-SAX1C6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {7,S} {11,S} {12,S}
3  C u0 p0 c0 {5,S} {6,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {18,S} {19,S} {20,S}
6  C u0 p0 c0 {3,S} {7,D} {8,S}
7  C u0 p0 c0 {2,S} {6,D} {21,S}
8  C u1 p0 c0 {6,S} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.77567,0.091789,-6.08038e-05,2.0792e-08,-2.95671e-12,2140.94,40.8158], Tmin=(298,'K'), Tmax=(1386,'K')),
            NASAPolynomial(coeffs=[24.7273,0.0332979,-1.15745e-05,1.8165e-09,-1.06184e-13,-7512.48,-102.916], Tmin=(1386,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=CCCC)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 313,
    label = "CH3-2-SAXC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {6,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {5,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {4,S} {6,D} {7,S}
6  C u0 p0 c0 {2,S} {5,D} {18,S}
7  C u1 p0 c0 {5,S} {19,S} {20,S}
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
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.985041,0.0776868,-5.06249e-05,1.70539e-08,-2.38473e-12,4844.78,35.1632], Tmin=(298,'K'), Tmax=(1388,'K')),
            NASAPolynomial(coeffs=[20.9478,0.0293621,-1.01371e-05,1.58364e-09,-9.22776e-14,-3144.16,-83.7983], Tmin=(1388,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(C)=CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 314,
    label = "C6H10-12",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {6,D} {14,S}
5  C u0 p0 c0 {6,D} {15,S} {16,S}
6  C u0 p0 c0 {4,D} {5,D}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.464396,0.0580709,-3.27949e-05,6.37514e-09,4.80422e-13,12183.2,25.8374], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[14.6061,0.0250445,-8.25216e-06,1.28846e-09,-7.82653e-14,7757.19,-49.2595], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/29/ 8 G""",
    longDesc = 
u"""
8/29/ 8 G
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C=CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 315,
    label = "CH2-3-1C6H10",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {5,S} {6,D}
5  C u0 p0 c0 {4,S} {7,D} {15,S}
6  C u0 p0 c0 {4,D} {18,S} {19,S}
7  C u0 p0 c0 {5,D} {16,S} {17,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {7,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.50244,0.0873984,-7.11094e-05,3.08712e-08,-5.48967e-12,1603.47,39.7299], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[20.9347,0.026808,-9.15645e-06,1.42006e-09,-8.23222e-14,-6027.82,-84.197], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(=C)CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 316,
    label = "CH2-3-1C7H12",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {11,S} {12,S}
3  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {6,S} {7,D}
6  C u0 p0 c0 {5,S} {8,D} {18,S}
7  C u0 p0 c0 {5,D} {21,S} {22,S}
8  C u0 p0 c0 {6,D} {19,S} {20,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {8,S}
20 H u0 p0 c0 {8,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.75518,0.0997015,-7.94484e-05,3.37635e-08,-5.8954e-12,-1280.99,42.5884], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[24.1478,0.03118,-1.06396e-05,1.64908e-09,-9.55598e-14,-10127.1,-100.001], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(=C)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 317,
    label = "C2H3-2-1C4H7",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u0 p0 c0 {3,S} {6,D} {12,S}
5  C u0 p0 c0 {3,D} {15,S} {16,S}
6  C u0 p0 c0 {4,D} {13,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.36165,0.0756316,-6.35843e-05,2.84791e-08,-5.19264e-12,4503.72,37.3812], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[17.7401,0.0224088,-7.66144e-06,1.18895e-09,-6.89534e-14,-1932.99,-68.4923], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CC(=C)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 318,
    label = "C2H5-3-TAX1C6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {9,S} {10,S}
2  C u0 p0 c0 {1,S} {6,S} {13,S} {14,S}
3  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
5  C u0 p0 c0 {3,S} {15,S} {16,S} {17,S}
6  C u1 p0 c0 {2,S} {3,S} {7,S}
7  C u0 p0 c0 {6,S} {8,D} {21,S}
8  C u0 p0 c0 {7,D} {22,S} {23,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.148128,0.0931725,-7.24799e-05,3.56574e-08,-8.69997e-12,786.101,30.3662], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[21.59,0.0352142,-1.08204e-05,1.60817e-09,-9.44551e-14,-5734.79,-83.4683], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[C](CC)CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 319,
    label = "CH3-2-1C6H11",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {2,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {6,S} {17,S} {18,S} {19,S}
6  C u0 p0 c0 {3,S} {5,S} {7,D}
7  C u0 p0 c0 {6,D} {20,S} {21,S}
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
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.812686,0.0800476,-5.21693e-05,1.77322e-08,-2.50123e-12,-12158.8,35.4576], Tmin=(298,'K'), Tmax=(1392,'K')),
            NASAPolynomial(coeffs=[21.0362,0.0312252,-1.06455e-05,1.64909e-09,-9.55267e-14,-20023.2,-82.7626], Tmin=(1392,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 320,
    label = "PXC2H4-2-1C5H9",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {8,S} {9,S}
2  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
3  C u0 p0 c0 {5,S} {6,S} {12,S} {13,S}
4  C u0 p0 c0 {1,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {3,S} {7,D}
6  C u1 p0 c0 {3,S} {17,S} {18,S}
7  C u0 p0 c0 {5,D} {19,S} {20,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.936837,0.0792075,-5.47955e-05,2.00972e-08,-3.08369e-12,12680.2,38.0172], Tmin=(298,'K'), Tmax=(1393,'K')),
            NASAPolynomial(coeffs=[20.6949,0.0291192,-9.95159e-06,1.54407e-09,-8.95422e-14,5062.97,-78.4167], Tmin=(1393,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]CC(=C)CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 321,
    label = "CH2-3-1C5H8",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u0 p0 c0 {3,S} {6,D} {12,S}
5  C u0 p0 c0 {3,D} {15,S} {16,S}
6  C u0 p0 c0 {4,D} {13,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[5.69172,0.0348843,9.72853e-06,-1.891e-08,5.1071e-12,-8914.98,4.26763], Tmin=(298,'K'), Tmax=(1368,'K')),
            NASAPolynomial(coeffs=[20.8165,0.0261803,-1.05808e-05,1.81995e-09,-1.12948e-13,-17278.2,-87.159], Tmin=(1368,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
Duplicate of species C2H3-2-1C4H7 (i.e. same molecular structure according to RMG)
C=CC(=C)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 322,
    label = "SAX4-5C10H19",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
3  C u0 p0 c0 {1,S} {8,S} {17,S} {18,S}
4  C u0 p0 c0 {6,S} {7,S} {11,S} {12,S}
5  C u0 p0 c0 {2,S} {9,S} {19,S} {20,S}
6  C u0 p0 c0 {4,S} {10,S} {21,S} {22,S}
7  C u0 p0 c0 {4,S} {23,S} {24,S} {25,S}
8  C u0 p0 c0 {3,S} {26,S} {27,S} {28,S}
9  C u0 p0 c0 {5,S} {10,D} {29,S}
10 C u1 p0 c0 {6,S} {9,D}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {8,S}
28 H u0 p0 c0 {8,S}
29 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.86658,0.11671,-7.72536e-05,2.62645e-08,-3.6614e-12,-3176.94,50.2855], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[30.2417,0.0427836,-1.46951e-05,2.288e-09,-1.33013e-13,-15069.6,-128.827], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC[C]=CCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 323,
    label = "C5H8-12",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {6,S} {7,S}
2  C u0 p0 c0 {1,S} {8,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,D} {11,S}
4  C u0 p0 c0 {5,D} {12,S} {13,S}
5  C u0 p0 c0 {3,D} {4,D}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.732228,0.0457648,-2.44733e-05,3.37843e-09,9.86638e-13,15071.7,22.9274], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[11.7029,0.0205508,-6.87242e-06,1.08557e-09,-6.65381e-14,11624.9,-35.4139], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""8/29/ 8 G""",
    longDesc = 
u"""
8/29/ 8 G
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C=CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 324,
    label = "SAXC4H8-2-1C6H11",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
3  C u0 p0 c0 {5,S} {6,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {8,S} {19,S} {20,S}
5  C u0 p0 c0 {3,S} {9,S} {17,S} {18,S}
6  C u0 p0 c0 {3,S} {21,S} {22,S} {23,S}
7  C u0 p0 c0 {2,S} {24,S} {25,S} {26,S}
8  C u0 p0 c0 {4,S} {9,D} {10,S}
9  C u0 p0 c0 {5,S} {8,D} {27,S}
10 C u1 p0 c0 {8,S} {28,S} {29,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {7,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.5334,0.117787,-7.98831e-05,2.81643e-08,-4.11206e-12,-3596.39,47.6485], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[30.6964,0.0424789,-1.4607e-05,2.27594e-09,-1.32376e-13,-15439.4,-131.758], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=CCCC)CCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 325,
    label = "C3H7-2-SAXC6H10",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {3,S} {6,S} {12,S} {13,S}
3  C u0 p0 c0 {2,S} {7,S} {16,S} {17,S}
4  C u0 p0 c0 {1,S} {8,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
6  C u0 p0 c0 {2,S} {21,S} {22,S} {23,S}
7  C u0 p0 c0 {3,S} {8,D} {9,S}
8  C u0 p0 c0 {4,S} {7,D} {24,S}
9  C u1 p0 c0 {7,S} {25,S} {26,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.63471,0.107773,-7.60328e-05,2.78298e-08,-4.26644e-12,-2064.68,45.2305], Tmin=(298,'K'), Tmax=(1379,'K')),
            NASAPolynomial(coeffs=[23.9718,0.0387116,-1.16059e-05,1.69698e-09,-9.86775e-14,-9916.41,-94.0437], Tmin=(1379,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERG""",
    longDesc = 
u"""
THERG
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=CCCC)CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 326,
    label = "SAXC10H19",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {13,S} {14,S}
2  C u0 p0 c0 {1,S} {3,S} {15,S} {16,S}
3  C u0 p0 c0 {2,S} {5,S} {17,S} {18,S}
4  C u0 p0 c0 {1,S} {6,S} {11,S} {12,S}
5  C u0 p0 c0 {3,S} {7,S} {19,S} {20,S}
6  C u0 p0 c0 {4,S} {8,S} {21,S} {22,S}
7  C u0 p0 c0 {5,S} {23,S} {24,S} {25,S}
8  C u1 p0 c0 {6,S} {9,S} {26,S}
9  C u0 p0 c0 {8,S} {10,D} {27,S}
10 C u0 p0 c0 {9,D} {28,S} {29,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {2,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.4934,0.116496,-7.73332e-05,2.62845e-08,-3.64631e-12,-1884.01,47.3556], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[30.6443,0.0423089,-1.45032e-05,2.25521e-09,-1.30991e-13,-13745.3,-131.81], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[CH]CCCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 327,
    label = "SAXC9H17",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {12,S} {13,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {1,S} {5,S} {10,S} {11,S}
4  C u0 p0 c0 {2,S} {6,S} {16,S} {17,S}
5  C u0 p0 c0 {3,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {4,S} {20,S} {21,S} {22,S}
7  C u1 p0 c0 {5,S} {8,S} {23,S}
8  C u0 p0 c0 {7,S} {9,D} {24,S}
9  C u0 p0 c0 {8,D} {25,S} {26,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.20061,0.103997,-6.87317e-05,2.32482e-08,-3.21153e-12,995.159,44.3168], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[27.3846,0.0379931,-1.30426e-05,2.02999e-09,-1.17985e-13,-9626.53,-115.738], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[CH]CCCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 328,
    label = "SAXC8H15",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {11,S} {12,S}
2  C u0 p0 c0 {1,S} {4,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {13,S} {14,S}
4  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
5  C u0 p0 c0 {3,S} {17,S} {18,S} {19,S}
6  C u1 p0 c0 {4,S} {7,S} {20,S}
7  C u0 p0 c0 {6,S} {8,D} {21,S}
8  C u0 p0 c0 {7,D} {22,S} {23,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.90099,0.0915068,-6.01588e-05,2.02338e-08,-2.78235e-12,3872.14,41.2376], Tmin=(298,'K'), Tmax=(1390,'K')),
            NASAPolynomial(coeffs=[24.1485,0.0336466,-1.15693e-05,1.80262e-09,-1.04847e-13,-5516.31,-99.7977], Tmin=(1390,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C[CH]CCCCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 329,
    label = "C2H5-2-SAX1C5H9",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {4,S} {5,S} {10,S} {11,S}
2  C u0 p0 c0 {3,S} {6,S} {8,S} {9,S}
3  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
4  C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
5  C u0 p0 c0 {1,S} {6,D} {7,S}
6  C u0 p0 c0 {2,S} {5,D} {18,S}
7  C u1 p0 c0 {5,S} {19,S} {20,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.74716,0.0806992,-5.45801e-05,1.92864e-08,-2.84516e-12,5054.42,38.9555], Tmin=(298,'K'), Tmax=(1389,'K')),
            NASAPolynomial(coeffs=[21.0643,0.0293146,-1.01321e-05,1.58405e-09,-9.23493e-14,-3139.65,-84.3685], Tmin=(1389,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[CH2]C(=CCC)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 330,
    label = "SOOcC6O",
    molecule = 
"""
multiplicity 3
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {4,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {6,S} {9,S} {18,S}
4  C u0 p0 c0 {2,S} {5,S} {12,S} {13,S}
5  C u0 p0 c0 {4,S} {6,S} {14,S} {15,S}
6  C u0 p0 c0 {3,S} {5,S} {16,S} {17,S}
7  O u0 p2 c0 {1,S} {19,S}
8  H u0 p0 c0 {1,S}
9  O u1 p2 c0 {3,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {3,S}
19 O u1 p2 c0 {7,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-6.33046,0.0986108,-7.17945e-05,2.51724e-08,-3.43846e-12,-41677.8,58.259], Tmin=(298,'K'), Tmax=(1381,'K')),
            NASAPolynomial(coeffs=[26.5003,0.0241046,-8.59209e-06,1.37146e-09,-8.11269e-14,-53187,-118.659], Tmin=(1381,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""therm""",
    longDesc = 
u"""
therm
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OC1CCCCC1[O]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 331,
    label = "C2H5cC6H10OO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {7,S} {10,S}
2  C u0 p0 c0 {5,S} {6,S} {9,S} {11,S}
3  C u0 p0 c0 {1,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {6,S} {20,S} {21,S}
5  C u0 p0 c0 {2,S} {3,S} {16,S} {17,S}
6  C u0 p0 c0 {2,S} {4,S} {18,S} {19,S}
7  C u0 p0 c0 {1,S} {8,S} {12,S} {13,S}
8  C u0 p0 c0 {7,S} {22,S} {23,S} {24,S}
9  O u0 p2 c0 {2,S} {25,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {7,S}
13 H u0 p0 c0 {7,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 O u1 p2 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.13112,0.109046,-6.14824e-05,1.36791e-08,-1.98043e-12,-17893.1,52.4648], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[48.3903,0.00637459,-6.61953e-07,3.05349e-11,-2.74783e-17,-37336.3,-240.933], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1CCC(CC1)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 332,
    label = "CH3cC6H10OO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {7,S} {9,S}
2  C u0 p0 c0 {5,S} {6,S} {8,S} {10,S}
3  C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
4  C u0 p0 c0 {1,S} {6,S} {17,S} {18,S}
5  C u0 p0 c0 {2,S} {3,S} {13,S} {14,S}
6  C u0 p0 c0 {2,S} {4,S} {15,S} {16,S}
7  C u0 p0 c0 {1,S} {19,S} {20,S} {21,S}
8  O u0 p2 c0 {2,S} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {6,S}
16 H u0 p0 c0 {6,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 O u1 p2 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.72202,0.0957705,-5.08545e-05,8.39858e-09,-6.70857e-13,-15291.6,48.9529], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[44.5644,0.00365816,-3.204e-07,2.21912e-11,-1.02989e-15,-33406.5,-222.036], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1CCC(CC1)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 333,
    label = "C4H9cC6H10OO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {5,S} {12,S}
2  C u0 p0 c0 {6,S} {7,S} {11,S} {13,S}
3  C u0 p0 c0 {1,S} {8,S} {18,S} {19,S}
4  C u0 p0 c0 {1,S} {6,S} {20,S} {21,S}
5  C u0 p0 c0 {1,S} {7,S} {26,S} {27,S}
6  C u0 p0 c0 {2,S} {4,S} {22,S} {23,S}
7  C u0 p0 c0 {2,S} {5,S} {24,S} {25,S}
8  C u0 p0 c0 {3,S} {9,S} {16,S} {17,S}
9  C u0 p0 c0 {8,S} {10,S} {14,S} {15,S}
10 C u0 p0 c0 {9,S} {28,S} {29,S} {30,S}
11 O u0 p2 c0 {2,S} {31,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {9,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {3,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {4,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {6,S}
24 H u0 p0 c0 {7,S}
25 H u0 p0 c0 {7,S}
26 H u0 p0 c0 {5,S}
27 H u0 p0 c0 {5,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
30 H u0 p0 c0 {10,S}
31 O u1 p2 c0 {11,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-5.38391,0.131707,-7.34607e-05,1.50356e-08,-1.35753e-12,-23700.6,57.0886], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[55.1728,0.0133965,-2.23711e-06,2.1634e-10,-9.80726e-15,-45408.1,-273.921], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC1CCC(CC1)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 334,
    label = "C2H5cC6H9OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {4,S} {6,S} {11,S}
2  C u0 p0 c0 {1,S} {5,S} {17,S} {18,S}
3  C u0 p0 c0 {5,S} {8,S} {9,S} {14,S}
4  C u0 p0 c0 {1,S} {7,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {3,S} {15,S} {16,S}
6  C u0 p0 c0 {1,S} {8,S} {19,S} {20,S}
7  C u0 p0 c0 {4,S} {21,S} {22,S} {23,S}
8  C u1 p0 c0 {3,S} {6,S} {24,S}
9  O u0 p2 c0 {3,S} {10,S}
10 O u0 p2 c0 {9,S} {25,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {2,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.5224,0.105375,-5.88346e-05,1.29048e-08,-1.93808e-12,-12041.8,47.1613], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[48.59,0.00602046,-6.77608e-07,5.12377e-11,-2.18745e-15,-31036.5,-238.742], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1C[CH]C(CC1)OO
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 335,
    label = "C2H5-2-C4H513",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
3  C u0 p0 c0 {1,S} {4,S} {5,D}
4  C u0 p0 c0 {3,S} {6,D} {12,S}
5  C u0 p0 c0 {3,D} {15,S} {16,S}
6  C u0 p0 c0 {4,D} {13,S} {14,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {2,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.36165,0.0756316,-6.35843e-05,2.84791e-08,-5.19264e-12,4503.72,37.3812], Tmin=(298,'K'), Tmax=(1397,'K')),
            NASAPolynomial(coeffs=[17.7401,0.0224088,-7.66144e-06,1.18895e-09,-6.89534e-14,-1932.99,-68.4923], Tmin=(1397,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
Duplicate of species C2H3-2-1C4H7 (i.e. same molecular structure according to RMG)
C=CC(=C)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 336,
    label = "C4H9cC6H9OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {3,S} {8,S} {13,S}
2  C u0 p0 c0 {1,S} {5,S} {18,S} {19,S}
3  C u0 p0 c0 {1,S} {6,S} {23,S} {24,S}
4  C u0 p0 c0 {6,S} {10,S} {11,S} {20,S}
5  C u0 p0 c0 {2,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {3,S} {4,S} {21,S} {22,S}
7  C u0 p0 c0 {5,S} {9,S} {14,S} {15,S}
8  C u0 p0 c0 {1,S} {10,S} {25,S} {26,S}
9  C u0 p0 c0 {7,S} {27,S} {28,S} {29,S}
10 C u1 p0 c0 {4,S} {8,S} {30,S}
11 O u0 p2 c0 {4,S} {12,S}
12 O u0 p2 c0 {11,S} {31,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {7,S}
15 H u0 p0 c0 {7,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {2,S}
19 H u0 p0 c0 {2,S}
20 H u0 p0 c0 {4,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {3,S}
24 H u0 p0 c0 {3,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {9,S}
29 H u0 p0 c0 {9,S}
30 H u0 p0 c0 {10,S}
31 H u0 p0 c0 {12,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.63307,0.127038,-6.83934e-05,1.18328e-08,-4.53671e-13,-17864.3,51.1852], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[55.6605,0.0125194,-1.95018e-06,1.72035e-10,-6.98824e-15,-39220.6,-273.328], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC1C[CH]C(CC1)OO
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 337,
    label = "CH3cC6H9OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {10,S}
2  C u0 p0 c0 {1,S} {4,S} {14,S} {15,S}
3  C u0 p0 c0 {4,S} {7,S} {8,S} {11,S}
4  C u0 p0 c0 {2,S} {3,S} {12,S} {13,S}
5  C u0 p0 c0 {1,S} {7,S} {16,S} {17,S}
6  C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
7  C u1 p0 c0 {3,S} {5,S} {21,S}
8  O u0 p2 c0 {3,S} {9,S}
9  O u0 p2 c0 {8,S} {22,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {2,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {5,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.17721,0.092551,-4.9312e-05,8.74579e-09,-1.03103e-12,-9433.61,43.9187], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[44.8411,0.00316995,-2.50622e-07,1.73628e-11,-8.67259e-16,-27139.1,-220.278], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1C[CH]C(CC1)OO
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 338,
    label = "C3H7cC6H10OO",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
2  C u0 p0 c0 {6,S} {7,S} {10,S} {12,S}
3  C u0 p0 c0 {1,S} {8,S} {15,S} {16,S}
4  C u0 p0 c0 {1,S} {6,S} {17,S} {18,S}
5  C u0 p0 c0 {1,S} {7,S} {23,S} {24,S}
6  C u0 p0 c0 {2,S} {4,S} {19,S} {20,S}
7  C u0 p0 c0 {2,S} {5,S} {21,S} {22,S}
8  C u0 p0 c0 {3,S} {9,S} {13,S} {14,S}
9  C u0 p0 c0 {8,S} {25,S} {26,S} {27,S}
10 O u0 p2 c0 {2,S} {28,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 H u0 p0 c0 {8,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {6,S}
21 H u0 p0 c0 {7,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {5,S}
24 H u0 p0 c0 {5,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 O u1 p2 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-4.97381,0.118408,-6.27424e-05,9.63862e-09,-6.0619e-17,-20827.2,53.5755], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[51.8071,0.00983107,-1.42033e-06,1.21343e-10,-5.01572e-15,-41380.1,-257.564], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC1CCC(CC1)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 339,
    label = "CH3-2-1C5H9",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {5,S} {9,S} {10,S}
3  C u0 p0 c0 {1,S} {11,S} {12,S} {13,S}
4  C u0 p0 c0 {5,S} {14,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {4,S} {6,D}
6  C u0 p0 c0 {5,D} {17,S} {18,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
11 H u0 p0 c0 {3,S}
12 H u0 p0 c0 {3,S}
13 H u0 p0 c0 {3,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {6,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.511775,0.0675516,-4.35931e-05,1.47175e-08,-2.07209e-12,-9281.96,32.3726], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[17.7919,0.0268883,-9.17598e-06,1.42233e-09,-8.24251e-14,-15908.9,-66.7741], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)CCC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 340,
    label = "CH3-2-1C4H7",
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
            NASAPolynomial(coeffs=[-0.276631,0.0553539,-3.545e-05,1.19627e-08,-1.69878e-12,-5992.86,29.5905], Tmin=(298,'K'), Tmax=(1391,'K')),
            NASAPolynomial(coeffs=[14.5614,0.0225371,-7.701e-06,1.19468e-09,-6.92705e-14,-11397.8,-50.8659], Tmin=(1391,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""THERM""",
    longDesc = 
u"""
THERM
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=C(C)CC
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 341,
    label = "C3H7cC6H9OOH",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {3,S} {4,S} {9,S} {12,S}
2  C u0 p0 c0 {5,S} {7,S} {10,S} {13,S}
3  C u0 p0 c0 {1,S} {5,S} {14,S} {15,S}
4  C u0 p0 c0 {1,S} {6,S} {16,S} {17,S}
5  C u0 p0 c0 {2,S} {3,S} {20,S} {21,S}
6  C u0 p0 c0 {4,S} {8,S} {18,S} {19,S}
7  C u0 p0 c0 {2,S} {9,S} {22,S} {23,S}
8  C u0 p0 c0 {6,S} {24,S} {25,S} {26,S}
9  C u1 p0 c0 {1,S} {7,S} {27,S}
10 O u0 p2 c0 {2,S} {11,S}
11 O u0 p2 c0 {10,S} {28,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 H u0 p0 c0 {3,S}
15 H u0 p0 c0 {3,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {4,S}
18 H u0 p0 c0 {6,S}
19 H u0 p0 c0 {6,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {7,S}
23 H u0 p0 c0 {7,S}
24 H u0 p0 c0 {8,S}
25 H u0 p0 c0 {8,S}
26 H u0 p0 c0 {8,S}
27 H u0 p0 c0 {9,S}
28 H u0 p0 c0 {11,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-3.80959,0.117823,-6.75149e-05,1.62749e-08,-2.58101e-12,-14928.4,50.154], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[52.1916,0.00914445,-1.24443e-06,1.00587e-10,-3.96989e-15,-35152.2,-256.398], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC1[CH]CC(CC1)OO
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 342,
    label = "CH3cC6H9O3",
    molecule = 
"""
multiplicity 3
1  C u0 p0 c0 {3,S} {4,S} {7,S} {9,S}
2  C u0 p0 c0 {3,S} {5,S} {8,S} {10,S}
3  C u0 p0 c0 {1,S} {2,S} {17,S} {18,S}
4  C u0 p0 c0 {1,S} {6,S} {12,S} {13,S}
5  C u0 p0 c0 {2,S} {6,S} {15,S} {16,S}
6  C u0 p0 c0 {4,S} {5,S} {11,S} {14,S}
7  C u0 p0 c0 {1,S} {19,S} {20,S} {21,S}
8  O u0 p2 c0 {2,S} {22,S}
9  H u0 p0 c0 {1,S}
10 H u0 p0 c0 {2,S}
11 O u1 p2 c0 {6,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {5,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {3,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {7,S}
22 O u1 p2 c0 {8,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.83831,0.0971288,-7.32222e-05,2.15381e-08,6.69395e-13,-46976.3,21.6509], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[-2.99562,0.0777889,-3.36105e-05,6.34831e-09,-4.40556e-13,-42656.9,53.0096], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CC1CC([O])CC(C1)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 343,
    label = "C2H5cC6H9O3",
    molecule = 
"""
multiplicity 3
1  C u0 p0 c0 {3,S} {4,S} {6,S} {10,S}
2  C u0 p0 c0 {3,S} {5,S} {9,S} {11,S}
3  C u0 p0 c0 {1,S} {2,S} {20,S} {21,S}
4  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
5  C u0 p0 c0 {2,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {1,S} {8,S} {13,S} {14,S}
7  C u0 p0 c0 {4,S} {5,S} {12,S} {17,S}
8  C u0 p0 c0 {6,S} {22,S} {23,S} {24,S}
9  O u0 p2 c0 {2,S} {25,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
12 O u1 p2 c0 {7,S}
13 H u0 p0 c0 {6,S}
14 H u0 p0 c0 {6,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {3,S}
21 H u0 p0 c0 {3,S}
22 H u0 p0 c0 {8,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {8,S}
25 O u1 p2 c0 {9,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.654059,0.108765,-7.96518e-05,2.2302e-08,1.09137e-12,-49601,24.2182], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[-1.28729,0.0840545,-3.56692e-05,6.65288e-09,-4.5759e-13,-45743.2,45.9514], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCC1CC([O])CC(C1)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 344,
    label = "C3H7cC6H9O3",
    molecule = 
"""
multiplicity 3
1  C u0 p0 c0 {3,S} {4,S} {5,S} {11,S}
2  C u0 p0 c0 {3,S} {6,S} {10,S} {12,S}
3  C u0 p0 c0 {1,S} {2,S} {23,S} {24,S}
4  C u0 p0 c0 {1,S} {8,S} {14,S} {15,S}
5  C u0 p0 c0 {1,S} {7,S} {18,S} {19,S}
6  C u0 p0 c0 {2,S} {7,S} {21,S} {22,S}
7  C u0 p0 c0 {5,S} {6,S} {13,S} {20,S}
8  C u0 p0 c0 {4,S} {9,S} {16,S} {17,S}
9  C u0 p0 c0 {8,S} {25,S} {26,S} {27,S}
10 O u0 p2 c0 {2,S} {28,S}
11 H u0 p0 c0 {1,S}
12 H u0 p0 c0 {2,S}
13 O u1 p2 c0 {7,S}
14 H u0 p0 c0 {4,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {8,S}
17 H u0 p0 c0 {8,S}
18 H u0 p0 c0 {5,S}
19 H u0 p0 c0 {5,S}
20 H u0 p0 c0 {7,S}
21 H u0 p0 c0 {6,S}
22 H u0 p0 c0 {6,S}
23 H u0 p0 c0 {3,S}
24 H u0 p0 c0 {3,S}
25 H u0 p0 c0 {9,S}
26 H u0 p0 c0 {9,S}
27 H u0 p0 c0 {9,S}
28 O u1 p2 c0 {10,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.31432,0.121469,-8.8618e-05,2.55715e-08,6.35157e-13,-52480.6,27.4455], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[0.392137,0.090407,-3.78021e-05,6.97832e-09,-4.76484e-13,-49095.6,39.0374], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCC1CC([O])CC(C1)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

entry(
    index = 345,
    label = "C4H9cC6H9O3",
    molecule = 
"""
multiplicity 3
1  C u0 p0 c0 {3,S} {4,S} {5,S} {12,S}
2  C u0 p0 c0 {3,S} {6,S} {11,S} {13,S}
3  C u0 p0 c0 {1,S} {2,S} {26,S} {27,S}
4  C u0 p0 c0 {1,S} {7,S} {15,S} {16,S}
5  C u0 p0 c0 {1,S} {8,S} {21,S} {22,S}
6  C u0 p0 c0 {2,S} {8,S} {24,S} {25,S}
7  C u0 p0 c0 {4,S} {9,S} {17,S} {18,S}
8  C u0 p0 c0 {5,S} {6,S} {14,S} {23,S}
9  C u0 p0 c0 {7,S} {10,S} {19,S} {20,S}
10 C u0 p0 c0 {9,S} {28,S} {29,S} {30,S}
11 O u0 p2 c0 {2,S} {31,S}
12 H u0 p0 c0 {1,S}
13 H u0 p0 c0 {2,S}
14 O u1 p2 c0 {8,S}
15 H u0 p0 c0 {4,S}
16 H u0 p0 c0 {4,S}
17 H u0 p0 c0 {7,S}
18 H u0 p0 c0 {7,S}
19 H u0 p0 c0 {9,S}
20 H u0 p0 c0 {9,S}
21 H u0 p0 c0 {5,S}
22 H u0 p0 c0 {5,S}
23 H u0 p0 c0 {8,S}
24 H u0 p0 c0 {6,S}
25 H u0 p0 c0 {6,S}
26 H u0 p0 c0 {3,S}
27 H u0 p0 c0 {3,S}
28 H u0 p0 c0 {10,S}
29 H u0 p0 c0 {10,S}
30 H u0 p0 c0 {10,S}
31 O u1 p2 c0 {11,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.119856,0.133167,-9.51817e-05,2.64605e-08,1.01394e-12,-55375.8,30.0574], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.21875,0.0965298,-3.98312e-05,7.28456e-09,-4.94096e-13,-52510.6,31.2898], Tmin=(1000,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CCCCC1CC([O])CC(C1)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/JetSurF2.0/Thermdat.txt.
""",
)

