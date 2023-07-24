#!/usr/bin/env python
# encoding: utf-8

name = "thermo"
shortDesc = ""
longDesc = """
Calculated using Arkane v3.1.0 using LevelOfTheory(method='ccsd(t)f12',basis='ccpvtzf12',software='molpro').
"""
entry(
    index = 0,
    label = "methoxy",
    molecule = 
"""
multiplicity 2
1 O u1 p2 c0 {2,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.07167,-0.00476973,3.68981e-05,-4.32415e-08,1.66496e-11,1107.57,5.40999], Tmin=(10,'K'), Tmax=(774.398,'K')),
            NASAPolynomial(coeffs=[1.26204,0.01394,-7.47207e-06,1.95474e-09,-2.00626e-13,1416.87,17.434], Tmin=(774.398,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (9.225,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (108.088,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Spin multiplicity: 2
External symmetry: -1.0
Optical isomers: 1

Geometry:
C      -0.03960000    0.02954200    0.02901200
O       1.30913600   -0.01479900   -0.19258600
H      -0.30905400   -0.17052800    1.07206200
H      -0.50351900    0.95270000   -0.33586400
H      -0.45696300   -0.79691700   -0.57262300
""",
)

entry(
    index = 1,
    label = "aziridine",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {3,S} {8,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C u0 p0 c0 {1,S} {2,S} {6,S} {7,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.18439,-0.0129536,9.86204e-05,-1.27273e-07,5.37393e-11,13568.2,6.67328], Tmin=(10,'K'), Tmax=(732.475,'K')),
            NASAPolynomial(coeffs=[-0.839572,0.0289451,-1.68003e-05,4.73523e-09,-5.17013e-13,13916.2,26.6984], Tmin=(732.475,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (112.846,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (182.918,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Spin multiplicity: 1
External symmetry: -1.0
Optical isomers: 1

Geometry:
C      -0.67962600   -0.39174100    0.03988900
N      -0.18463800    0.86200800   -0.54883200
C       0.77347400   -0.13091200   -0.03929800
H      -1.14266500   -1.06822400   -0.66520400
H      -1.14865200   -0.35569200    1.01484000
H      -0.28404700    1.61862000    0.11912700
H       1.30075900    0.08397300    0.88135900
H       1.36539600   -0.61803200   -0.80188100
""",
)

entry(
    index = 2,
    label = "hydrazino",
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
            NASAPolynomial(coeffs=[3.98757,0.000752,2.9792e-05,-5.31982e-08,3.07588e-11,25693.3,4.96394], Tmin=(10,'K'), Tmax=(512.722,'K')),
            NASAPolynomial(coeffs=[3.00718,0.0106617,-5.81457e-06,1.70064e-09,-2.03514e-13,25764.1,8.74901], Tmin=(512.722,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (213.622,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (103.931,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Spin multiplicity: 2
External symmetry: -1.0
Optical isomers: 1

Geometry:
N      -0.42485700    0.00192800    0.16965000
N       0.78887900   -0.52686600   -0.08286700
H      -1.18373200   -0.64514000    0.02966600
H      -0.63118400    0.95854400   -0.08521800
H       1.43572400    0.26253400   -0.11178700
""",
)

entry(
    index = 3,
    label = "1-propene-12-diol",
    molecule = 
"""
1  O u0 p2 c0 {4,S} {10,S}
2  O u0 p2 c0 {5,S} {11,S}
3  C u0 p0 c0 {4,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {1,S} {3,S} {5,D}
5  C u0 p0 c0 {2,S} {4,D} {9,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {1,S}
11 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.89618,0.0068893,0.000140192,-3.14252e-07,2.16762e-10,-42687.8,10.8383], Tmin=(10,'K'), Tmax=(467.501,'K')),
            NASAPolynomial(coeffs=[2.14227,0.0403926,-2.66526e-05,8.30273e-09,-9.83016e-13,-42726,15.8046], Tmin=(467.501,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (-354.945,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (245.277,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Spin multiplicity: 1
External symmetry: -1.0
Optical isomers: 2

Geometry:
C      -1.40993900   -0.02702000   -0.31272400
C       0.02134800    0.21936400    0.00119300
O       0.29152500    1.51620900    0.33085500
C       0.98308600   -0.70412300   -0.03871600
O       2.29775300   -0.30540400    0.18585800
H      -1.72468500    0.58052200   -1.16300500
H      -2.03611000    0.24714400    0.53787200
H      -1.57772700   -1.07568100   -0.55056000
H       1.24358400    1.57033400    0.50142600
H       0.80450100   -1.73104200   -0.32113800
H       2.71984000   -0.92815700    0.78618100
""",
)

entry(
    index = 4,
    label = "hydroxyiminomethyl",
    molecule = 
"""
multiplicity 2
1 O u0 p2 c0 {2,S} {5,S}
2 N u0 p1 c0 {1,S} {3,D}
3 C u1 p0 c0 {2,D} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.9626,0.00210306,4.51447e-05,-8.36646e-08,4.60304e-11,29495,6.91197], Tmin=(10,'K'), Tmax=(603.513,'K')),
            NASAPolynomial(coeffs=[3.21138,0.0160647,-1.18823e-05,3.9922e-09,-4.96815e-13,29422.1,8.80152], Tmin=(603.513,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (245.218,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (103.931,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Spin multiplicity: 2
External symmetry: -1.0
Optical isomers: 1

Geometry:
C      -0.95896200   -0.52629100   -0.06633600
N       0.16750300   -0.09855900   -0.18865700
O       0.84289500    0.51227800    0.94604900
H      -1.74506400   -0.63440100    0.66801500
H       1.69362900    0.74697400    0.55253800
""",
)

entry(
    index = 5,
    label = "2-methyl-2-propanamine",
    molecule = 
"""
1  N u0 p1 c0 {2,S} {15,S} {16,S}
2  C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
4  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
5  C u0 p0 c0 {2,S} {12,S} {13,S} {14,S}
6  H u0 p0 c0 {3,S}
7  H u0 p0 c0 {3,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
14 H u0 p0 c0 {5,S}
15 H u0 p0 c0 {1,S}
16 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.87476,0.00856898,0.000199706,-4.62508e-07,3.36983e-10,-17641.8,9.26793], Tmin=(10,'K'), Tmax=(437.548,'K')),
            NASAPolynomial(coeffs=[1.54289,0.0528722,-3.09731e-05,9.02657e-09,-1.03462e-12,-17657.8,16.0761], Tmin=(437.548,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (-146.693,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (365.837,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Spin multiplicity: 1
External symmetry: -1.0
Optical isomers: 1

Geometry:
C       0.62055100    1.27852200   -0.00644000
C      -0.07932400   -0.07495300    0.09283100
C      -1.24106700   -0.13404100   -0.90703500
C       0.91151800   -1.19976100   -0.19805600
N      -0.52963200   -0.23545700    1.48493800
H      -0.07280500    2.09002900    0.22645900
H       1.45000300    1.32791000    0.69842700
H       1.00386900    1.44560500   -1.01382200
H      -0.89275700   -0.01366600   -1.93501000
H      -1.75942400   -1.09266700   -0.83705900
H      -1.96490000    0.65746100   -0.70174200
H       1.74246700   -1.16312900    0.50582400
H       1.30501100   -1.11934500   -1.21214000
H       0.42791100   -2.17478500   -0.10328900
H      -1.01179300   -1.12194100    1.59298700
H      -1.20095900    0.48926000    1.71756200
""",
)

entry(
    index = 6,
    label = "ethynol",
    molecule = 
"""
1 O u0 p2 c0 {2,S} {5,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.9512,0.00295995,4.42163e-05,-9.64785e-08,6.1581e-11,8769.68,5.2066], Tmin=(10,'K'), Tmax=(546.862,'K')),
            NASAPolynomial(coeffs=[4.42481,0.0104758,-6.5163e-06,2.0836e-09,-2.61588e-13,8553.7,1.70644], Tmin=(546.862,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (72.8979,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (108.088,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Spin multiplicity: 1
External symmetry: -1.0
Optical isomers: 1

Geometry:
C       1.11285800   -0.05358000   -0.08834700
C      -0.07375200    0.13731600   -0.11674400
O      -1.35426300    0.41331900   -0.10446300
H       2.15805200   -0.22790800   -0.06712500
H      -1.84289700   -0.26914900   -0.58285200
""",
)

entry(
    index = 7,
    label = "2-iminoethyl",
    molecule = 
"""
multiplicity 2
1 N u0 p1 c0 {2,D} {7,S}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 C u1 p0 c0 {2,S} {5,S} {6,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.02742,-0.00164503,5.72236e-05,-8.08628e-08,3.64964e-11,23337.2,7.37022], Tmin=(10,'K'), Tmax=(665.027,'K')),
            NASAPolynomial(coeffs=[0.776356,0.0238763,-1.37995e-05,3.8271e-09,-4.12476e-13,23637.6,20.7365], Tmin=(665.027,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (194.04,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (149.66,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Spin multiplicity: 2
External symmetry: -1.0
Optical isomers: 1

Geometry:
C      -1.04285900   -0.17328000   -0.01422600
C       0.29857400    0.21824400    0.00060400
N       0.62355300    1.48299300   -0.01986500
H      -1.32272300   -1.21564600    0.00243900
H      -1.81477900    0.58205500   -0.04292800
H       1.05175300   -0.57230500    0.02971000
H       1.63747200    1.58682000   -0.00470500
""",
)

