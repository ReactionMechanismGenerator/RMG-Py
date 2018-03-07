#!/usr/bin/env python
# encoding: utf-8

name = "/home/alongd/Code/RMG-Py/importer/NOx2018"
shortDesc = u"/home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt"
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
            NASAPolynomial(coeffs=[4.23854,-0.000249611,1.59858e-05,-2.0692e-08,8.29766e-12,-17648.6,3.5885], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.54017,0.00415971,-1.30877e-06,2.00824e-10,-1.15509e-14,-17951.4,0.855882], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

H2O2 <g> ATcT ver. 1.122, DHf298 = -135.457 ? 0.064 kJ/mol - fit JAN17.
OO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
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
            NASAPolynomial(coeffs=[3.97585,-0.00228555,4.33443e-06,-3.59927e-09,1.26707e-12,3393.41,-0.0355397], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.84582,0.00109724,-2.89121e-07,4.091e-11,-2.31382e-15,3717.07,5.8034], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E
O <g, triplet> ATcT ver. 1.122, DHf298 = 249.229 ? 0.002 kJ/mol - fit JAN17
O <g, singlet> ATcT ver. 1.122, DHf298 = 438.523 ? 0.002 kJ/mol - fit JAN17

OH <g> ATcT ver. 1.122, DHf298 = 37.490 ? 0.027 kJ/mol - fit JAN17.
[OH]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
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
            NASAPolynomial(coeffs=[2.37694,0.00773917,-1.88735e-05,1.95517e-08,-7.17096e-12,-921.173,0.547185], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.90208,0.000868993,-1.65864e-07,1.90852e-11,-9.31122e-16,-797.949,-0.845591], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E


********************************************************************************
********************************************************************************
***** H2 O2 O3 H O OH HO2 H2O H2O2
********************************************************************************
********************************************************************************

H2 <g> ATcT ver. 1.122, DHf298 = 0.000 ? 0.000 kJ/mol - fit JAN17.
[H][H]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
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
            NASAPolynomial(coeffs=[2.49976,6.73824e-07,1.11807e-09,-3.70192e-12,2.14234e-15,25473.8,-0.445574], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.49985,2.34583e-07,-1.16172e-10,2.25708e-14,-1.52992e-18,25473.8,-0.445865], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

H <g> ATcT ver. 1.122, DHf298 = 217.998 ? 0.000 kJ/mol - fit JAN17.
[H]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
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
            NASAPolynomial(coeffs=[3.15907,-0.0032151,6.49256e-06,-5.98755e-09,2.06876e-12,29129.8,2.09078], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.5516,-3.83085e-05,8.43197e-10,4.01267e-12,-4.17477e-16,29228.8,4.87617], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

O <g> ATcT ver. 1.122, DHf298 = 249.229 ? 0.002 kJ/mol - fit JAN17.
[O]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
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
            NASAPolynomial(coeffs=[4.20148,-0.00205584,6.56547e-06,-5.52907e-09,1.78283e-12,-30295,-0.860611], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.73118,0.00295137,-8.3536e-07,1.26089e-10,-8.40532e-15,-29916.9,6.55183], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

H2O <g> ATcT ver. 1.122, DHf298 = -241.833 ? 0.027 kJ/mol - fit JAN17.
O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
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
            NASAPolynomial(coeffs=[4.26251,-0.00445642,2.05165e-05,-2.35794e-08,9.05614e-12,262.442,3.88224], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.10564,0.00204047,-3.65878e-07,1.85973e-11,4.98818e-16,43.2899,3.30808], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

HO2 <g> ATcT ver. 1.122, DHf298 = 12.26 ? 0.16 kJ/mol - fit JAN17.
[O]O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
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
            NASAPolynomial(coeffs=[3.78498,-0.00302002,9.92029e-06,-9.7784e-09,3.28878e-12,-1064.14,3.64781], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.6598,0.000659877,-1.44158e-07,2.14656e-11,-1.36504e-15,-1216.03,3.42074], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

O2 <g> ATcT ver. 1.122, DHf298 = 0.000 ? 0.000 kJ/mol - fit JAN17.
[O][O]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
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
            NASAPolynomial(coeffs=[2.49989,2.13038e-07,8.97321e-10,-2.31396e-12,1.30201e-15,-745.354,4.38024], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.49989,1.56135e-07,-7.76109e-11,1.52928e-14,-1.05304e-18,-745.328,4.3803], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

Ar <g> ATcT ver. 1.122, DHf298 = 0.000 ? 0.000 kJ/mol - fit JAN17.
[Ar]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
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
            NASAPolynomial(coeffs=[3.53604,-0.000158271,-4.26984e-07,2.37543e-09,-1.39708e-12,-1047.5,2.94604], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.93803,0.00141838,-5.03281e-07,8.07555e-11,-4.76064e-15,-917.181,5.95522], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E


********************************************************************************
********************************************************************************
***** N2 AR HE
********************************************************************************
********************************************************************************

N2 <g> ATcT ver. 1.122, DHf298 = 0.000 ? 0.000 kJ/mol - fit JAN17.
N#N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
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
            NASAPolynomial(coeffs=[2.49976,1.01013e-06,-8.24578e-10,-6.85983e-13,7.24752e-16,-745.341,0.9298], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.49986,2.19365e-07,-1.07525e-10,2.07198e-14,-1.39359e-18,-745.309,0.929535], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

He <g> ATcT ver. 1.122, DHf298 = 0.000 ? 0.000 kJ/mol - fit JAN17.
[He]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 12,
    label = "CO",
    molecule = 
"""
1 C u0 p1 c-1 {2,T}
2 O u0 p1 c+1 {1,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.59508,-0.000721197,1.28238e-06,6.52429e-10,-8.21715e-13,-14344.9,3.44356], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.03397,0.00137328,-4.96445e-07,8.10281e-11,-4.85332e-15,-14258.6,6.10076], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E


********************************************************************************
********************************************************************************
***** CO CO2
********************************************************************************
********************************************************************************

CO <g> ATcT ver. 1.122, DHf298 = -110.523 ? 0.026 kJ/mol - fit JAN17.
[C-]#[O+]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 13,
    label = "CO2",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,D}
2 O u0 p2 c0 {1,D}
3 O u0 p2 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.20664,0.010097,-9.96339e-06,5.47156e-09,-1.27734e-12,-48353,10.5262], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.63537,0.00274559,-9.98282e-07,1.61014e-10,-9.22019e-15,-49020.4,-1.92888], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

CO2 <g> ATcT ver. 1.122, DHf298 = -393.475 ? 0.015 kJ/mol - fit JAN17.
O=C=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 14,
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
            NASAPolynomial(coeffs=[5.00702,-0.0126484,4.66821e-05,-4.59211e-08,1.57634e-11,-10222.4,-4.04227], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[1.68377,0.010013,-3.31268e-06,5.30234e-10,-3.13372e-14,-10018.8,9.71477], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E


********************************************************************************
********************************************************************************
***** CH4 CH3 CH2 CH2(S) CH C
********************************************************************************
********************************************************************************

CH4 <g> ATcT ver. 1.122, DHf298 = -74.519 ? 0.057 kJ/mol - fit JAN17.
C
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 15,
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
            NASAPolynomial(coeffs=[3.61264,0.00309209,9.25475e-07,-1.65777e-09,6.07244e-13,16385,1.79995], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.92198,0.00537479,-1.99748e-06,2.97585e-10,-1.7186e-14,16544.7,5.25397], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

CH3 <g> ATcT ver. 1.122, DHf298 = 146.374 ? 0.080 kJ/mol - fit JAN17.
[CH3]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 16,
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
            NASAPolynomial(coeffs=[5.52796,-0.0138922,6.21716e-05,-6.82045e-08,2.52019e-11,-25596.8,-0.715601], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.71776,0.0100962,-3.53785e-06,5.61843e-10,-3.32714e-14,-26022.2,4.04183], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E
C <g, triplet> ATcT ver. 1.122, DHf298 = 716.886 ? 0.050 kJ/mol - fit JAN17
C <g, singlet> ATcT ver. 1.122, DHf298 = 838.478 ? 0.050 kJ/mol - fit JAN17


********************************************************************************
********************************************************************************
***** CH3OH CH3O CH2OH CH2O HCOH HCO HOCO
********************************************************************************
********************************************************************************

CH3OH <g> ATcT ver. 1.122, DHf298 = -200.71 ? 0.16 kJ/mol - fit JAN17.
CO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 17,
    label = "CH",
    molecule = 
"""
multiplicity 4
1 C u3 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.49602,0.000277424,-1.57434e-06,3.07039e-09,-1.40468e-12,70650,2.05802], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.74989,0.00136248,-2.71161e-07,2.47828e-11,-1.18068e-15,70899.9,6.13979], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

CH <g> ATcT ver. 1.122, DHf298 = 596.12 ? 0.11 kJ/mol - fit JAN17.
[CH]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 18,
    label = "C",
    molecule = 
"""
1 C u0 p2 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.55138,-0.000299926,6.79035e-07,-6.75577e-10,2.46001e-13,85469.3,4.54305], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.60607,-0.000196384,1.0681e-07,-1.6408e-11,8.15424e-16,85437.8,4.18935], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E
CH <g, doublet> ATcT ver. 1.122, DHf298 = 596.12 ? 0.11 kJ/mol - fit JAN17
CH <g, quartet> ATcT ver. 1.122, DHf298 = 667.88 ? 0.59 kJ/mol - fit JAN17

C <g> ATcT ver. 1.122, DHf298 = 716.886 ? 0.050 kJ/mol - fit JAN17.
[C]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 19,
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
            NASAPolynomial(coeffs=[4.14633,-0.00442384,5.70676e-05,-6.76419e-08,2.57224e-11,-11514.5,3.26636], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.97869,0.0154587,-5.37527e-06,8.59587e-10,-5.11298e-14,-12444.6,-0.62717], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E


********************************************************************************
********************************************************************************
***** C2H6 C2H5 C2H4 C2H3 C2H2 H2CC C2H C2
********************************************************************************
********************************************************************************

C2H6 <g> ATcT ver. 1.122, DHf298 = -83.91 ? 0.14 kJ/mol - fit JAN17.
CC
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 20,
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
            NASAPolynomial(coeffs=[4.11105,-0.00263844,4.60475e-05,-5.64273e-08,2.19028e-11,13001.4,4.99418], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.41717,0.0122704,-4.34558e-06,6.94963e-10,-4.13531e-14,12146.8,-0.387435], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

CH3CH2 <g> ATcT ver. 1.122, DHf298 = 119.86 ? 0.28 kJ/mol - fit JAN17.
C[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 21,
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
            NASAPolynomial(coeffs=[3.65151,-0.00535067,5.16486e-05,-6.36869e-08,2.50743e-11,5114.51,5.38561], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.14446,0.0102648,-3.61247e-06,5.74009e-10,-3.39296e-14,4190.59,-1.14778], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

C2H4 <g> ATcT ver. 1.122, DHf298 = 52.45 ? 0.13 kJ/mol - fit JAN17.
C=C
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 22,
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
    shortDesc = u"""Ethynyl Rad   T 5/10""",
    longDesc = 
u"""
Ethynyl Rad   T 5/10.
[C]#C
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 23,
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
            NASAPolynomial(coeffs=[3.71623,-0.00283501,3.77227e-05,-4.73648e-08,1.86739e-11,1358.74,6.55373], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.75205,0.00744944,-2.70076e-06,4.38793e-10,-2.64006e-14,444.371,-1.93362], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

CH3O <g> ATcT ver. 1.122, DHf298 = 21.53 ? 0.34 kJ/mol - fit JAN17.
C[O]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 24,
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
            NASAPolynomial(coeffs=[4.40537,-0.000758358,2.61926e-05,-3.46319e-08,1.4082e-11,-3442.28,3.60691], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.20504,0.00243666,1.94377e-06,-1.61013e-09,3.14355e-13,-4258.44,-7.46071], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

CH2OH <g> ATcT ver. 1.122, DHf298 = -16.57 ? 0.33 kJ/mol - fit JAN17.
[CH2]O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 25,
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
            NASAPolynomial(coeffs=[4.77187,-0.00976266,3.70122e-05,-3.76922e-08,1.31327e-11,-14379.8,0.696586], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.91333,0.0067004,-2.55521e-06,4.27795e-10,-2.44073e-14,-14462.2,7.43823], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

H2CO <g> ATcT ver. 1.122, DHf298 = -109.188 ? 0.099 kJ/mol - fit JAN17.
C=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 26,
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
            NASAPolynomial(coeffs=[3.12502,0.00235137,2.36803e-05,-3.35092e-08,1.39444e-11,34524.3,8.81538], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.37211,0.00746869,-2.64716e-06,4.22753e-10,-2.44958e-14,33805.2,0.428772], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

C2H3 <g> ATcT ver. 1.122, DHf298 = 296.91 ? 0.33 kJ/mol - fit JAN17.
[CH]=C
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 27,
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
            NASAPolynomial(coeffs=[3.97075,-0.00149122,9.54042e-06,-8.8272e-09,2.67645e-12,3842.03,4.4466], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.85781,0.00264114,-7.44177e-07,1.23313e-10,-8.88959e-15,3616.43,3.92451], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E
REF=ATcT C 2011}  Max Lst Sq Error Cp @ 6000 K 0.49%.
from Xiong et al., Combust. Flame 161 (2014) 885?897
3-HCOH                  C   1H   2O   1    0G   298.150  6000.000 1000.00      1
0.62184153E+01 0.18701090E-02-0.17117529E-06-0.44242689E-10 0.63999362E-14    2
0.23486822E+05-0.97996712E+01 0.36513441E+01 0.46834047E-02 0.11223500E-05    3
-0.80289814E-09-0.77469460E-12 0.24561953E+05 0.49211311E+01                   4
1-HCOH                  C   1H   2O   1    0G   298.150  6000.000 1000.00      1
0.57648468E+01 0.21295433E-02-0.19747641E-06-0.51074096E-10 0.74327861E-14    2
0.10621689E+05-0.77241540E+01 0.32130404E+01 0.32250737E-02 0.24648079E-05    3
0.15449875E-08-0.27946393E-11 0.11899702E+05 0.76449280E+01                   4
CHOH trans Hydroxymethylene trans  SIGMA=1  STATWT=1  IA=0.2970  IB=2.3197
IC=2.6167  Nu=3634,2827,1536,1330,1224,1099  REF=Burcat G3B3  HF298=108.16+/-.43
kJ  REF=Ruscic ATcT D 2013  {HF298=106.47+/-8 kJ  REF=Burcat G3B3}  Max Lst Sq
Error Cp @ 6000 K 0.45%.
CHOH trans  Hydr  T12/14C  1.H  2.O  1.   0.G   200.000  6000.000  B  30.02598 1
3.63500546E+00 5.45479239E-03-1.91285077E-06 3.03882849E-10-1.79904440E-14    2
5.28968001E+04 4.34543184E+00 4.72862471E+00-1.02082828E-02 4.21639528E-05    3
-4.64774522E-08 1.72559970E-11 5.31829869E+04 1.69093263E+00 5.44279146E+04    4

HCO <g> ATcT ver. 1.122, DHf298 = 41.803 ? 0.099 kJ/mol - fit JAN17.
[CH]=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 28,
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
            NASAPolynomial(coeffs=[0.626264,0.0247264,-3.89556e-05,3.16485e-08,-9.87351e-12,26455.6,14.6968], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.66166,0.00488862,-1.61257e-06,2.48284e-10,-1.39737e-14,25769.5,-4.01062], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

HCCH <g> ATcT ver. 1.122, DHf298 = 228.27 ? 0.13 kJ/mol - fit JAN17.
C#C
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 29,
    label = "CH2CHOO",
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
            NASAPolynomial(coeffs=[1.41461,0.0297546,-2.58394e-05,1.12221e-08,-1.86781e-12,11371.2,19.0226], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[1.41603,0.0299337,-2.63852e-05,1.17707e-08,-2.05116e-12,11361.7,18.9695], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = u"""0""",
    longDesc = 
u"""
0
J Gimenez CL Rasmussen MU Alzueta P Marshall P Glarborg Proc. Combust. Inst. 32 (2009) 367-375
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 30,
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
            NASAPolynomial(coeffs=[4.1849,-0.00216397,7.68524e-06,-5.98243e-09,1.6791e-12,50392.3,-0.717772], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.86931,0.00339757,-1.006e-06,1.50563e-10,-8.61238e-15,50628.5,5.53164], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E
CH2 <g, singlet> ATcT ver. 1.122, DHf298 = 429.03 ? 0.13 kJ/mol - fit JAN17.
[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 31,
    label = "CH2CHOOH",
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
            NASAPolynomial(coeffs=[2.79821,0.0274377,-2.03468e-05,7.62127e-09,-1.12671e-12,-6315.53,9.11829], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.82519,0.0274167,-2.04456e-05,7.77399e-09,-1.18661e-12,-6325.27,8.96641], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = u"""0""",
    longDesc = 
u"""
0
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
CH2CH2OOH               H  5 C  2 O  2      G     298.0    6000.0    1000.0    1
9.94719358E+00 1.14730988E-02-3.96870805E-06 6.23121554E-10-3.65488844E-14    2
2.33352925E+03-2.17163906E+01 6.49571381E-01 4.51634030E-02-6.02115954E-05    3
4.92121132E-08-1.67753354E-11 4.73413790E+03 2.49288138E+01                   4
! Goldsmith et al., J. Phys. Chem. A 2012, 116, 3325?3346

JIM/GLA08
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=COO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 32,
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
            NASAPolynomial(coeffs=[3.74843,0.00107656,3.62462e-06,-4.46166e-09,1.95225e-12,45899,1.63509], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.81413,0.00379513,-7.38425e-07,7.19949e-11,-2.63189e-15,46186,6.52946], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""(t)           ATcT3E""",
    longDesc = 
u"""
(t)           ATcT3E

CH2 <g, triplet> ATcT ver. 1.122, DHf298 = 391.52 ? 0.12 kJ/mol - fit JAN17.
Duplicate of species CH2(S) (i.e. same molecular structure according to RMG)
[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 33,
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
            NASAPolynomial(coeffs=[2.66874,0.0096233,1.60617e-05,-2.87682e-08,1.2503e-11,219.438,12.5694], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.91637,0.0088465,-3.14955e-06,5.05413e-10,-3.01305e-14,-1047.8,-6.1065], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Vinyl-  T04/06""",
    longDesc = 
u"""
Vinyl-  T04/06
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016.
[CH2]C=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 34,
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
            NASAPolynomial(coeffs=[2.82191,0.00966218,-2.7856e-06,-4.12692e-09,2.61472e-12,-23546.5,11.4285], Tmin=(200,'K'), Tmax=(998.4,'K')),
            NASAPolynomial(coeffs=[4.63989,0.00566363,-2.67855e-06,6.17049e-10,-5.60954e-14,-24052.7,1.90175], Tmin=(998.4,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""FAB/JAN05""",
    longDesc = 
u"""
FAB/JAN05

H298 =-44.33 kcal/mol [FAB/JAN05]
S298 = 60.07 cal/mol/K [FAB/JAN05]
Cp [FAB/JAN05] (polyfit RAS/GLA08a).
O=[C]O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 35,
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
L 8/88
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016.
CC=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 36,
    label = "CH2CHOH",
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
            NASAPolynomial(coeffs=[2.28758,0.0197013,1.96383e-06,-1.9439e-08,1.02617e-11,-16537.3,14.1333], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.49818,0.0103957,-3.66891e-06,5.85206e-10,-3.47374e-14,-18164.3,-13.8388], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T03/10""",
    longDesc = 
u"""
T03/10
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016.
C=CO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 37,
    label = "CH3CH2O",
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
            NASAPolynomial(coeffs=[3.26906,0.00933563,2.96317e-05,-4.53411e-08,1.88796e-11,-2950.23,10.4201], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.55054,0.0132526,-4.74726e-06,7.64699e-10,-4.57008e-14,-4471.92,-9.61231], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T06/11""",
    longDesc = 
u"""
T06/11
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016.
CC[O]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 38,
    label = "CH3CHOH",
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
            NASAPolynomial(coeffs=[4.22283,0.00512175,3.48387e-05,-4.91944e-08,2.01184e-11,-8356.22,8.01676], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.35842,0.0124356,-4.33097e-06,6.8453e-10,-4.03713e-14,-9530.19,-6.05106], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T06/11""",
    longDesc = 
u"""
T06/11
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016.
C[CH]O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 39,
    label = "CHCHOH",
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
            NASAPolynomial(coeffs=[-1.99181,0.0473233,-6.66059e-05,4.68997e-08,-1.30686e-11,15200.1,31.4259], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.78246,0.00524797,-1.71857e-06,2.59722e-10,-1.48227e-14,12883.6,-21.0851], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical""",
    longDesc = 
u"""
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016.
[CH]=CO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 40,
    label = "CH3CH2OH",
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
            NASAPolynomial(coeffs=[4.8587,-0.00374017,6.95554e-05,-8.86548e-08,3.51688e-11,-29996.1,4.80185], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.56244,0.0152042,-5.38968e-06,8.6225e-10,-5.12898e-14,-31525.6,-9.47302], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 8/88""",
    longDesc = 
u"""
L 8/88
C2 <g, singlet> ATcT ver. 1.122, DHf298 = 826.78 ? 0.27 kJ/mol - fit JAN17
C2 <g, triplet> ATcT ver. 1.122, DHf298 = 834.08 ? 0.27 kJ/mol - fit JAN17


********************************************************************************
********************************************************************************
***** CH3CH2OH CH3CHOH CH3CH2O CH2CH2OH CH3CHO CH2CHOH
***** cC2H4O CHCHOH cC2H3O CH2CHO CH3CO CH2CO CHCHO HCCOH HCCO C2O
********************************************************************************
********************************************************************************.
CCO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 41,
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
    shortDesc = u"""Ethynol   T12/09""",
    longDesc = 
u"""
Ethynol   T12/09
Lopez et al., Experimental and Kinetic Modeling Study of C2H2 Oxidation at High Pressure, Int. J. Chem. Kin., 2016.
C#CO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 42,
    label = "CH2CH2OH",
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
            NASAPolynomial(coeffs=[4.20954,0.00912965,2.47462e-05,-3.92946e-08,1.66541e-11,-4915.11,8.30445], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.01349,0.0120204,-4.21992e-06,6.70676e-10,-3.97135e-14,-6161.62,-8.62052], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T05/11""",
    longDesc = 
u"""
T05/11
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016.
[CH2]CO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 43,
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
    shortDesc = u"""RADICAL    IU2/03""",
    longDesc = 
u"""
RADICAL    IU2/03
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016.
C[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 44,
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
            NASAPolynomial(coeffs=[3.89836,-0.00355878,3.55205e-05,-4.385e-08,1.71078e-11,-46778.6,7.34954], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.61383,0.00644964,-2.29083e-06,3.6716e-10,-2.18737e-14,-45330.3,0.847884], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""FORMIC ACID A 5/14""",
    longDesc = 
u"""
FORMIC ACID A 5/14
Fabian WMF Janoschek R J Mol Struct THEOCHEM 2005, 713, 227?234
CL Rasmussen J Hansen P Marshall P Glarborg Int J Chem Kinet 40 (2008) 454-480


********************************************************************************
********************************************************************************
***** HOCHO OCHO
********************************************************************************
********************************************************************************.
O=CO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 45,
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
            NASAPolynomial(coeffs=[2.13241,0.0181319,-1.74093e-05,9.35336e-09,-2.01725e-12,-7148.09,13.3808], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.75871,0.00635124,-2.25955e-06,3.62322e-10,-2.15856e-14,-8085.33,-4.9649], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""g 4/02""",
    longDesc = 
u"""
g 4/02
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016.
C=C=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 46,
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
            NASAPolynomial(coeffs=[3.6286,0.00812496,-1.41561e-06,-3.27952e-09,1.61554e-12,-16747.8,7.8317], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.41052,0.00750888,-4.2589e-06,1.12761e-09,-1.14144e-13,-17029.8,3.43148], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1104""",
    longDesc = 
u"""
1104
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016

fitted to litt data
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
[O]C=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 47,
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
T 4/09
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016.
C#C[O]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 48,
    label = "CH3CH2OO",
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
            NASAPolynomial(coeffs=[4.50099,0.00687965,4.74144e-05,-6.92287e-08,2.87395e-11,-5395.48,7.9149], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.88872,0.0135833,-4.91117e-06,7.92343e-10,-4.73526e-14,-7441.07,-19.079], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical""",
    longDesc = 
u"""
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016.
CCO[O]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 49,
    label = "HOCH2CH2OO",
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
            NASAPolynomial(coeffs=[0.556418,0.0412519,-3.10647e-05,1.2293e-08,-1.99691e-12,-22047.4,24.9983], Tmin=(298,'K'), Tmax=(800,'K')),
            NASAPolynomial(coeffs=[0.556418,0.0412519,-3.10647e-05,1.2293e-08,-1.99691e-12,-22047.4,24.9983], Tmin=(800,'K'), Tmax=(2500,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2500,'K'),
    ),
    shortDesc = u"""J Gimenez CL Rasmussen MU Alzueta P Marshall P Glarborg Proc. Combust. Inst. 32 (2009) 367-375""",
    longDesc = 
u"""
J Gimenez CL Rasmussen MU Alzueta P Marshall P Glarborg Proc. Combust. Inst. 32 (2009) 367-375
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OCCO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 50,
    label = "CH3CH2OOH",
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
            NASAPolynomial(coeffs=[4.14672,0.00978668,4.91492e-05,-7.42532e-08,3.11169e-11,-21467.1,9.84025], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[9.58691,0.0148604,-5.29788e-06,8.47317e-10,-5.03436e-14,-23836.8,-22.8311], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T02/10""",
    longDesc = 
u"""
T02/10
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016


********************************************************************************
********************************************************************************
***** CH3CH2OOH CH3CH2OO CH3CHOOH CH2CH2OOH CH2CHOOH CH2CHOO
***** HOCH2CH2OO
********************************************************************************
********************************************************************************.
CCOO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 51,
    label = "OCHCHO",
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
            NASAPolynomial(coeffs=[4.68412,0.000478013,4.26391e-05,-5.79018e-08,2.31669e-11,-27198.5,4.51187], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.72507,0.00633097,-2.35575e-06,3.89783e-10,-2.37487e-14,-29102.4,-20.3904], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""S Olivella A Sole J Phys Chem A 108 (2004) 11651?11663""",
    longDesc = 
u"""
S Olivella A Sole J Phys Chem A 108 (2004) 11651?11663
J Gimenez CL Rasmussen MU Alzueta P Marshall P Glarborg Proc. Combust. Inst. 32 (2009) 367-375


********************************************************************************
********************************************************************************
***** OCHCHO OCHCO
********************************************************************************
********************************************************************************.
O=CC=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 52,
    label = "CH3CHOOH",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2 C u1 p0 c0 {1,S} {3,S} {8,S}
3 O u0 p2 c0 {2,S} {4,S}
4 O u0 p2 c0 {3,S} {9,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[6.26854,0.000111351,6.15944e-05,-7.95035e-08,3.12333e-11,962.571,-0.180971], Tmin=(200,'K'), Tmax=(999.99,'K')),
            NASAPolynomial(coeffs=[9.10837,0.0152964,-5.54864e-06,9.02286e-10,-5.4412e-14,-932.703,-20.3914], Tmin=(999.99,'K'), Tmax=(2500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (2500,'K'),
    ),
    shortDesc = u"""E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical""",
    longDesc = 
u"""
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016


H298 =  6.43 kcal/mol [JAN/ROS04]
S298 = 74.84 cal/mol/K [JAN/ROS04]
Cp(T) scaled Cp[CH3CH2OO](T) to Cp298 = 19.70 cal/mol/K [JAN/ROS04].
C[CH]OO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 53,
    label = "OCHCO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,D}
3 O u1 p2 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 O u0 p2 c0 {2,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.39406,0.0143629,-8.84138e-06,1.60961e-09,3.20389e-13,-8684.14,10.592], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.94994,0.010163,-5.5772e-06,1.4572e-09,-1.47435e-13,-9096.51,2.57988], Tmin=(1000,'K'), Tmax=(2900,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2900,'K'),
    ),
    shortDesc = u"""0""",
    longDesc = 
u"""
0
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
Low T polynomial Tmin changed from 350.0 to 298.0 K when importing to RMG.
[O]C=C=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 54,
    label = "CH3OO",
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
            NASAPolynomial(coeffs=[2.93065,0.00868504,8.80315e-06,-1.39561e-08,5.0294e-12,227.483,12.8755], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.88425,0.0140068,-6.88364e-06,1.6379e-09,-1.53129e-13,-20.0433,11.8153], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""H. Hashemi, et al., High-Pressure Oxidation of Natural Gas: Methane, 2016""",
    longDesc = 
u"""
H. Hashemi, et al., High-Pressure Oxidation of Natural Gas: Methane, 2016

H. Hashemi, et al., High-Pressure Oxidation of Natural Gas: Methane, 2016
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 55,
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
            NASAPolynomial(coeffs=[-1.97309,0.0426372,-5.84444e-05,4.12219e-08,-1.13585e-11,-21382.8,32.869], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.64142,0.0146739,-8.24143e-06,2.2479e-09,-2.38671e-13,-22230.4,7.15858], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1311""",
    longDesc = 
u"""
1311

Marshall and Glarborg, Proc. Combust. Inst. 35 (2015) 153?160.
[O]CO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 56,
    label = "CH3OOH",
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
            NASAPolynomial(coeffs=[2.85837,0.0159701,-2.52104e-06,-5.736e-09,2.80535e-12,-16883,12.897], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.45931,0.0141051,-6.53177e-06,1.47663e-09,-1.32512e-13,-17430.1,4.03867], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""Marshall and Glarborg, Proc. Combust. Inst. 35 (2015) 153?160""",
    longDesc = 
u"""
Marshall and Glarborg, Proc. Combust. Inst. 35 (2015) 153?160


********************************************************************************
********************************************************************************
***** CH3OOH CH3OO CH2OOH HOCH2O
********************************************************************************
********************************************************************************
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
COO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 57,
    label = "cC2H3O",
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
A 1/05
Goldsmith et al., J. Phys. Chem. A 2012, 116, 3325?3346.
[CH]1CO1
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 58,
    label = "CH2OOH",
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
            NASAPolynomial(coeffs=[5.83127,-0.00351771,4.54551e-05,-5.66903e-08,2.21633e-11,6061.87,-0.579143], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.98746,0.00900484,-3.24367e-06,5.24325e-10,-3.13587e-14,5012.58,-10.2619], Tmin=(1000,'K'), Tmax=(2500,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (2500,'K'),
    ),
    shortDesc = u"""Aranda, V., et al., Int. J. Chemical Kinet. 45.5 (2013) 283-294.""",
    longDesc = 
u"""
Aranda, V., et al., Int. J. Chemical Kinet. 45.5 (2013) 283-294.
H298 = 15.79 kcal/mol [JAN/ROS04]
S298 = 65.89 cal/mol/K [JAN/ROS04]
Cp(T) scaled Cp[CH3OO](T) to Cp298 = 14.89 cal/mol/K [JAN/ROS04].
[CH2]OO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 59,
    label = "cC2H4O",
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
    shortDesc = u"""OXYRANE    L 8/88""",
    longDesc = 
u"""
OXYRANE    L 8/88
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016.
C1CO1
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 60,
    label = "CH2CH2OOH",
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
            NASAPolynomial(coeffs=[5.07456,0.013239,2.53759e-05,-4.3809e-08,1.89062e-11,3710.63,2.7229], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[9.4138,0.0133542,-4.68068e-06,7.43231e-10,-4.39824e-14,1984.6,-22.4517], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 4/15""",
    longDesc = 
u"""
T 4/15
Janoschek R Rossi MJ Int J Chem Kinet 2004 36 661?686
CL Rasmussen JG Jacobsen P Glarborg Int J Chem Kinet 40 (2008) 778-807.
[CH2]COO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 61,
    label = "C2",
    molecule = 
"""
multiplicity 3
1 C u1 p0 c0 {2,T}
2 C u1 p0 c0 {1,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.33655,0.0359693,-0.00011185,1.31762e-07,-5.29904e-11,98422.1,9.54587], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.5587,0.000786412,-1.25558e-07,8.97793e-12,-1.32453e-16,98890.4,4.18695], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016

C2 <g> ATcT ver. 1.122, DHf298 = 828.67 ? 0.26 kJ/mol - fit JAN17.
[C]#[C]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 62,
    label = "CH3C(O)OO",
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
            NASAPolynomial(coeffs=[2.12838,0.0294919,-2.00221e-05,5.32494e-09,3.03452e-14,-21203.4,18.0452], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.76604,0.0188598,-9.95181e-06,2.52994e-09,-2.50541e-13,-22126.8,-0.48408], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1508""",
    longDesc = 
u"""
1508
Bozzelli 2015; PM.
CC(=O)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 63,
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
            NASAPolynomial(coeffs=[2.86278,0.0119701,-1.80851e-05,1.52778e-08,-5.20063e-12,44312.6,8.89759], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.42468,0.00185394,-5.17933e-07,6.77646e-11,-3.53315e-15,43716.1,-3.69608], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 8/11""",
    longDesc = 
u"""
T 8/11
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016.
[C]=C=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 64,
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
            NASAPolynomial(coeffs=[4.65733,-0.00953742,4.04679e-05,-4.45318e-08,1.64762e-11,13861.5,1.97861], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.65237,0.00555807,-1.97617e-06,3.16823e-10,-1.88748e-14,13553.6,4.22141], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""CH**-OH cis T 9/09""",
    longDesc = 
u"""
CH**-OH cis T 9/09.
[CH]O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 65,
    label = "NO",
    molecule = 
"""
multiplicity 2
1 N u1 p1 c0 {2,D}
2 O u0 p2 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.08518,-0.00364693,8.49608e-06,-6.62406e-09,1.77647e-12,9840.61,2.83578], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.25487,0.0011987,-4.33029e-07,7.02943e-11,-4.09789e-15,9907,6.40395], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

NO <g> ATcT ver. 1.122, DHf298 = 91.121 ? 0.065 kJ/mol - fit JAN17.
[N]=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 66,
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
            NASAPolynomial(coeffs=[4.14028,-0.00358489,1.89476e-05,-1.98834e-08,7.15268e-12,-6685.45,-0.0166755], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.36074,0.0063185,-2.28967e-06,4.11767e-10,-2.90837e-14,-6415.96,8.02154], Tmin=(1000,'K'), Tmax=(4000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (4000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E
in Marshall and Glarborg, Proc. Combust. Inst. 35 (2015) 153?160


********************************************************************************
********************************************************************************
***** NH3 NH2 NH N N2H4 N2H3 N2H2 H2NN NNH
********************************************************************************
********************************************************************************

NH3 <g> ATcT ver. 1.122, DHf298 = -45.554 ? 0.030 kJ/mol - fit JAN17.
N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 67,
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
            NASAPolynomial(coeffs=[4.06463,-0.00110021,4.25849e-06,-2.68224e-09,5.89267e-13,21176.9,0.439851], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.62499,0.00339841,-1.01631e-06,1.25511e-10,-2.66501e-15,21541.9,7.73537], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

NH2 <g> ATcT ver. 1.122, DHf298 = 186.02 ? 0.12 kJ/mol - fit JAN17.
[NH2]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 68,
    label = "NH",
    molecule = 
"""
multiplicity 3
1 N u2 p1 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.45887,0.000493904,-1.87863e-06,2.85542e-09,-1.16865e-12,42108.8,2.00373], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.79499,0.0012926,-3.85559e-07,6.26028e-11,-3.70422e-15,42340.9,5.68414], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

NH <g> ATcT ver. 1.122, DHf298 = 358.77 ? 0.17 kJ/mol - fit JAN17.
[NH]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 69,
    label = "N",
    molecule = 
"""
multiplicity 4
1 N u3 p1 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.49977,5.0215e-07,1.93091e-09,-4.94633e-12,2.7409e-15,56076.1,4.19499], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.41604,0.000174664,-1.18865e-07,3.0185e-11,-2.0326e-15,56105.2,4.64906], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

N <g> ATcT ver. 1.122, DHf298 = 472.440 ? 0.024 kJ/mol - fit JAN17.
[N]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 70,
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
            NASAPolynomial(coeffs=[4.25475,-0.00345098,1.37789e-05,-1.33264e-08,4.41023e-12,28832.4,3.28552], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.42744,0.00323295,-1.17296e-06,1.90508e-10,-1.14492e-14,28806.8,6.39209], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 1/06""",
    longDesc = 
u"""
T 1/06
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
N=[N]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 71,
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
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 72,
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
Rad.       T 7/11
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
N[NH]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 73,
    label = "HNO",
    molecule = 
"""
1 N u0 p1 c0 {2,D} {3,S}
2 O u0 p2 c0 {1,D}
3 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.55326,-0.00584532,1.88854e-05,-1.7604e-08,5.7289e-12,11631.6,1.66851], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.24129,0.00272377,-1.60633e-07,-9.79135e-11,1.17104e-14,11774.6,7.27914], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.
HNOH cis  ATcT C  T11/11H  2.N  1.O  1.   0.G   200.000  6000.000 1000.        1
4.11664692E+00 4.81707273E-03-1.63507639E-06 2.53797646E-10-1.47744717E-14    2
1.25020921E+04 3.12195287E+00 3.80983976E+00 4.35965662E-04 1.51571801E-05    3
-1.96181113E-08 7.75279218E-12 1.28164979E+04 5.90835846E+00 1.40705826E+04    4
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.

HNO <g> ATcT ver. 1.122, DHf298 = 106.96 ? 0.11 kJ/mol - fit JAN17.
N=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 74,
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
            NASAPolynomial(coeffs=[3.78713,-0.000429577,1.37384e-05,-1.74264e-08,6.7125e-12,2895,6.96592], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.90482,0.00214474,-8.12654e-07,1.55512e-10,-1.04114e-14,2289.59,-0.233567], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

ONO <g> ATcT ver. 1.122, DHf298 = 34.049 ? 0.065 kJ/mol - fit JAN17.
N(=O)[O]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 75,
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
equil & transT 9/11
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
N=N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 76,
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
            NASAPolynomial(coeffs=[3.16416,0.00850518,5.48562e-07,-8.27656e-09,4.39957e-12,-10774.4,10.0232], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.79145,0.00364631,-1.29113e-06,2.06498e-10,-1.22139e-14,-11597.4,-4.07145], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

HONO <g> ATcT ver. 1.122, DHf298 = -78.675 ? 0.079 kJ/mol - fit JAN17.
N(=O)O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 77,
    label = "HNO2",
    molecule = 
"""
1 N u0 p0 c+1 {2,S} {3,D} {4,S}
2 H u0 p0 c0 {1,S}
3 O u0 p2 c0 {1,D}
4 O u0 p3 c-1 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.03779,-0.00446123,3.19441e-05,-3.79359e-08,1.44571e-11,-6530.88,5.9062], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.66359,0.00489854,-1.79694e-06,2.9442e-10,-1.78236e-14,-7252.16,-0.0306054], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E
HONO <g, trans> ATcT ver. 1.122, DHf298 = -79.161 ? 0.079 kJ/mol - fit JAN17
HONO <g, cis> ATcT ver. 1.122, DHf298 = -77.72 ? 0.29 kJ/mol - fit JAN17

HN(O)O <g> ATcT ver. 1.122, DHf298 = -44.2 ? 1.5 kJ/mol - fit JAN17.
[NH+](=O)[O-]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 78,
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
            NASAPolynomial(coeffs=[3.35587,0.0106545,-2.8669e-06,-5.14712e-09,3.08532e-12,7475.35,8.94787], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.37569,0.00221733,-5.75696e-07,6.69775e-11,-2.58935e-15,6224.46,-12.4945], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

NO3 <g> ATcT ver. 1.122, DHf298 = 74.13 ? 0.19 kJ/mol - fit JAN17.
[N+](=O)([O-])[O]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 79,
    label = "C2H5NO2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  N u0 p0 c+1 {1,S} {9,D} {10,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {2,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  O u0 p2 c0 {3,D}
10 O u0 p3 c-1 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.4372,0.00825502,5.01014e-05,-7.12335e-08,2.88208e-11,-14991.3,9.09921], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[9.29092,0.0161242,-5.95123e-06,9.76015e-10,-5.90148e-14,-17371.9,-21.1218], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""NitroEth  T06/10""",
    longDesc = 
u"""
NitroEth  T06/10
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
C(C)[N+](=O)[O-]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 80,
    label = "CH2CH2NO2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 N u0 p0 c+1 {1,S} {6,D} {7,S}
3 C u1 p0 c0 {1,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 O u0 p2 c0 {2,D}
7 O u0 p3 c-1 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.19815,0.0383712,-3.03909e-05,1.17681e-08,-1.57118e-12,11332.9,22.1405], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.21841,0.022024,-1.1471e-05,2.8888e-09,-2.84836e-13,10142.1,-3.01283], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""0""",
    longDesc = 
u"""
0
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.

pm 1002
Low T polynomial Tmin changed from 350.0 to 298.0 K when importing to RMG.
C([N+](=O)[O-])[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 81,
    label = "CH3CH2ONO",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  O u0 p2 c0 {1,S} {4,S}
4  N u0 p1 c0 {3,S} {10,D}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 O u0 p2 c0 {4,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.37138,0.0137914,3.84688e-05,-6.02381e-08,2.49655e-11,-14019.9,15.8651], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[9.21849,0.0162002,-5.9816e-06,9.81277e-10,-5.93456e-14,-16554.4,-18.8592], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T04/98""",
    longDesc = 
u"""
T04/98
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.


H298 =-24.18 kcal/mol
S298 = 80.30 cal/mol/K
Cp(C2H5ONO) = Cp(C2H5NO2)[BURCAT].
C(C)ON=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 82,
    label = "HONO2",
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
            NASAPolynomial(coeffs=[1.55975,0.0201502,-1.15217e-05,-2.31891e-09,3.17581e-12,-17395.6,17.7295], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.03061,0.00446368,-1.72273e-06,2.91612e-10,-1.80487e-14,-19303.4,-16.2543], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

HNO3 <g> ATcT ver. 1.122, DHf298 = -134.19 ? 0.18 kJ/mol - fit JAN17.
[N+](=O)(O)[O-]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 83,
    label = "CH3CHNO2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {3,S} {7,S}
3 N u0 p0 c+1 {2,S} {8,D} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 O u0 p2 c0 {3,D}
9 O u0 p3 c-1 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.948311,0.037343,-2.65343e-05,7.9225e-09,-3.22058e-13,5946.38,23.8211], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.75403,0.0229073,-1.20613e-05,3.06099e-09,-3.0346e-13,4745.88,-0.560534], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""0""",
    longDesc = 
u"""
0

pm 1002
Low T polynomial Tmin changed from 350.0 to 298.0 K when importing to RMG.
C[CH][N+](=O)[O-]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 84,
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
T05/98
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
CO[N+](=O)[O-]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 85,
    label = "CH3CH2ONO2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  N u0 p0 c+1 {4,S} {10,D} {11,S}
4  O u0 p2 c0 {1,S} {3,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {2,S}
10 O u0 p2 c0 {3,D}
11 O u0 p3 c-1 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.75722,0.0193623,3.87534e-05,-6.6409e-08,2.82506e-11,-20844.4,11.1813], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[12.1361,0.0170091,-6.4374e-06,1.0722e-09,-6.54951e-14,-24190.2,-37.1641], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T05/98""",
    longDesc = 
u"""
T05/98

BURCAT
H298 =-37.04 kcal/mol
S298 = 78.59 cal/mol/K.
C(C)O[N+](=O)[O-]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 86,
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
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 87,
    label = "CH2NO2",
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
    shortDesc = u"""RADICAL  T08/07""",
    longDesc = 
u"""
RADICAL  T08/07
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
[N+](=O)([CH2])[O-]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 88,
    label = "HNOH",
    molecule = 
"""
multiplicity 2
1 N u1 p1 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.95608,-0.00302611,2.56874e-05,-3.15645e-08,1.24085e-11,10920,5.55951], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.98322,0.00488846,-1.65087e-06,2.55371e-10,-1.48309e-14,10578,3.62583], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""trans & Equ  T11/11""",
    longDesc = 
u"""
trans & Equ  T11/11
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
[NH]O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 89,
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
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 90,
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
NitrosomethyT12/09
Burcat 2006 access


********************************************************************************
********************************************************************************
***** CH3NO CH2NO C2H5NO CH3CHNO CH2CHNO CHCHNO
********************************************************************************
********************************************************************************.
CN=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 91,
    label = "C2H5NO",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 N u0 p1 c0 {1,S} {9,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 O u0 p2 c0 {3,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.19503,0.0284612,-1.06193e-05,-3.12526e-09,2.4706e-12,3551.62,20.6492], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.22313,0.0251564,-1.28736e-05,3.18527e-09,-3.0895e-13,2905.61,9.66286], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""0""",
    longDesc = 
u"""
0
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.

pm 0910
Low T polynomial Tmin changed from 350.0 to 298.0 K when importing to RMG.
C(C)N=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 92,
    label = "CH2NO",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 N u0 p1 c0 {1,S} {5,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 O u0 p2 c0 {2,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.87819,-0.00665309,5.39476e-05,-6.81768e-08,2.71817e-11,25716.9,7.46188], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.40282,0.0069057,-2.5163e-06,4.10141e-10,-2.47183e-14,24528.7,-4.45743], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""H2C*N=O     T 9/96""",
    longDesc = 
u"""
H2C*N=O     T 9/96
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
[CH2]N=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 93,
    label = "CH3CHNO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,D} {7,S}
3 N u0 p1 c0 {2,D} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 O u1 p2 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.92228,0.0244337,-9.49024e-06,-2.40398e-09,2.05477e-12,12313.6,16.9273], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.74755,0.0212797,-1.09796e-05,2.73592e-09,-2.66954e-13,11741.2,7.08475], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""0""",
    longDesc = 
u"""
0

pm 1002
Low T polynomial Tmin changed from 350.0 to 298.0 K when importing to RMG.
CC=N[O]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 94,
    label = "HCNO",
    molecule = 
"""
1 N u0 p0 c+1 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {4,S}
3 O u0 p3 c-1 {1,S}
4 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.607949,0.0282182,-4.60452e-05,3.82559e-08,-1.23227e-11,19071.4,16.9199], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.9198,0.00400115,-1.42063e-06,2.2757e-10,-1.35505e-14,18038.6,-8.26935], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Fulminic AcidA 5/05""",
    longDesc = 
u"""
Fulminic AcidA 5/05
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
[N+](#C)[O-]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 95,
    label = "HNCO",
    molecule = 
"""
1 N u0 p1 c0 {2,D} {4,S}
2 C u0 p0 c0 {1,D} {3,D}
3 O u0 p2 c0 {2,D}
4 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.24009,0.01456,-1.54352e-05,8.55535e-09,-1.79632e-12,-15459,12.1664], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.30045,0.00402251,-1.40962e-06,2.23855e-10,-1.325e-14,-16199.5,-3.11771], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Isocyanic AciA 5/05""",
    longDesc = 
u"""
Isocyanic AciA 5/05.
N=C=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 96,
    label = "NCO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {3,D}
2 N u1 p1 c0 {1,D}
3 O u0 p2 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.77405,0.00924523,-9.91774e-06,6.68461e-09,-2.09521e-12,14237,9.75459], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.08064,0.00237444,-9.07099e-07,1.52287e-10,-9.31009e-15,13578.1,-2.15734], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""(NCO)        A 5/05""",
    longDesc = 
u"""
(NCO)        A 5/05
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
C(=[N])=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 97,
    label = "H2NN",
    molecule = 
"""
multiplicity 3
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 N u2 p1 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.53204,-0.00732419,3.00804e-05,-3.04001e-08,1.04701e-11,34958,1.51074], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.05904,0.00618382,-2.22171e-06,3.58539e-10,-2.14533e-14,34853,6.69894], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Isodiazene   T 9/11""",
    longDesc = 
u"""
Isodiazene   T 9/11
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
N[N]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 98,
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
            NASAPolynomial(coeffs=[3.93201,-0.000164028,1.39161e-05,-1.62748e-08,6.00353e-12,6711.79,4.58837], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.75556,0.00516219,-1.76387e-06,2.75053e-10,-1.60643e-14,6518.26,4.30933], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""RADICAL     T09/09""",
    longDesc = 
u"""
RADICAL     T09/09
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
[NH2+][O-]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 99,
    label = "HON",
    molecule = 
"""
multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 N u2 p1 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.15272,-0.00387826,2.05476e-05,-2.49049e-08,9.87365e-12,24603.7,4.56636], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.12045,0.00228738,-7.14685e-07,1.03332e-10,-5.70484e-15,24364.4,3.38858], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

NOH <g> ATcT ver. 1.122, DHf298 = 214.57 ? 0.87 kJ/mol - fit JAN17.
O[N]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 100,
    label = "N2O",
    molecule = 
"""
1 N u0 p0 c+1 {2,D} {3,D}
2 N u0 p2 c-1 {1,D}
3 O u0 p2 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.13942,0.0121801,-1.59189e-05,1.2092e-08,-3.85126e-12,8870.09,11.2478], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.80641,0.00265307,-9.70797e-07,1.6259e-10,-9.96738e-15,8197.98,-2.10608], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

NNO <g> ATcT ver. 1.122, DHf298 = 82.569 ? 0.097 kJ/mol - fit JAN17.
[N+](=[N-])=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 101,
    label = "CN",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,T}
2 N u0 p1 c0 {1,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.6292,-0.00107173,2.42398e-06,-5.85516e-10,-3.75364e-13,51866.8,3.91231], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.73484,5.83996e-05,2.90444e-07,-6.74404e-11,4.33382e-15,51701.9,2.85163], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

CN <g> ATcT ver. 1.122, DHf298 = 440.01 ? 0.15 kJ/mol - fit JAN17.
[C]#N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 102,
    label = "HCN",
    molecule = 
"""
1 C u0 p0 c0 {2,T} {3,S}
2 N u0 p1 c0 {1,T}
3 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.10572,0.0111841,-1.62402e-05,1.31663e-08,-4.17254e-12,14542.7,9.55501], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.79525,0.0031633,-1.07317e-06,1.67852e-10,-9.83948e-15,14225,1.61357], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E


********************************************************************************
********************************************************************************
***** HCN HNC CN HNCO HOCN HCNO NCO
********************************************************************************
********************************************************************************

HCN <g> ATcT ver. 1.122, DHf298 = 129.275 ? 0.091 kJ/mol - fit JAN17.
C#N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 103,
    label = "HOCN",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,T}
2 O u0 p2 c0 {1,S} {4,S}
3 N u0 p1 c0 {1,T}
4 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.88944,0.0116487,-1.08005e-05,5.44139e-09,-1.06857e-12,-3152.97,9.51296], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.28768,0.00401747,-1.40407e-06,2.22563e-10,-1.31562e-14,-3774.1,-2.64471], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Cyanic Acid  A 5/05""",
    longDesc = 
u"""
Cyanic Acid  A 5/05
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
C(#N)O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 104,
    label = "HCNH",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,D} {3,S}
2 N u0 p1 c0 {1,D} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.97115,-0.00388876,2.92919e-05,-3.57482e-08,1.40304e-11,31578.9,5.06389], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.04015,0.00516592,-1.82277e-06,2.90299e-10,-1.71615e-14,31154,2.58894], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""H*C=NH TransT11/11""",
    longDesc = 
u"""
H*C=NH TransT11/11
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.

# trans or cis?.
[CH]=N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 105,
    label = "H2CN",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 N u1 p1 c0 {1,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.978,-0.00343276,2.59134e-05,-3.04692e-08,1.16273e-11,27485.4,4.43067], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.80316,0.00547197,-1.95315e-06,3.13362e-10,-1.86249e-14,27130.3,3.31759], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""H2C=N*      T 11/1""",
    longDesc = 
u"""
H2C=N*      T 11/1
DB00

E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
C=[N]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 106,
    label = "CH2NH",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 N u0 p1 c0 {1,D} {5,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.8185,0.00511983,6.38887e-06,-6.61375e-09,1.65532e-12,9884.43,10.3391], Tmin=(298,'K'), Tmax=(1577,'K')),
            NASAPolynomial(coeffs=[4.54738,0.00717721,-2.47935e-06,3.87692e-10,-2.26113e-14,8640.57,-1.16687], Tmin=(1577,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""MELIUS 88""",
    longDesc = 
u"""
MELIUS 88
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.

DB00
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 107,
    label = "CH3CHN",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,D} {7,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 N u1 p1 c0 {2,D}
7 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.81409,0.0132951,3.533e-06,-9.87941e-09,3.74397e-12,22555.6,12.3157], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.67758,0.0176422,-8.68917e-06,2.06973e-09,-1.93574e-13,22392.9,12.024], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1106""",
    longDesc = 
u"""
1106
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
CC=[N]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 108,
    label = "CH3CH2NH",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 N u1 p1 c0 {1,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.92366,0.0238729,-7.86535e-06,-2.71982e-09,1.89902e-12,16912.9,16.7828], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.09795,0.0224758,-1.07197e-05,2.48314e-09,-2.26728e-13,16513,10.2925], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1104""",
    longDesc = 
u"""
1104
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C(C)[NH]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 109,
    label = "CH3CHNH",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {1,S} {3,D} {7,S}
3 N u0 p1 c0 {2,D} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.24231,0.0146071,8.25387e-06,-1.56259e-08,5.70542e-12,4797.01,14.5454], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[1.99419,0.0212775,-1.02687e-05,2.40039e-09,-2.20648e-13,4562.73,14.3229], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1104""",
    longDesc = 
u"""
1104
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
CC=N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 110,
    label = "CH3CN",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,T}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 N u0 p1 c0 {2,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.82393,0.00408202,2.1621e-05,-2.89808e-08,1.12963e-11,7444.3,5.52656], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.09922,0.00969586,-3.48052e-06,5.6142e-10,-3.35836e-14,6609.67,-3.36087], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Methyl-Cya  T01/03""",
    longDesc = 
u"""
Methyl-Cya  T01/03


********************************************************************************
********************************************************************************
***** CH3CN CHCNH2 CH2CNH CH2CHN CH2CHN(S) c-C2H3N CH2CN CHCNH H2NCHO H2NCO
********************************************************************************
********************************************************************************.
CC#N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 111,
    label = "CH2CHNH",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {4,S}
2 C u1 p0 c0 {1,S} {5,S} {6,S}
3 N u0 p1 c0 {1,D} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.168729,0.0291658,-2.2416e-05,8.51021e-09,-1.09138e-12,23563.8,24.4877], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.67016,0.0167406,-8.17363e-06,1.94564e-09,-1.82823e-13,22649.5,5.23472], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1104""",
    longDesc = 
u"""
1104
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C(=N)[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 112,
    label = "CH3CNH",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {3,D}
3 N u0 p1 c0 {2,D} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.23656,0.0158045,1.22159e-06,-9.86661e-09,4.14655e-12,25585.4,15.0702], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.2967,0.0169751,-8.65109e-06,2.12636e-09,-2.04493e-13,25102.8,8.60277], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1106""",
    longDesc = 
u"""
1106
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C[C]=N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 113,
    label = "CH2CNH",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {4,S} {5,S}
2 C u0 p0 c0 {1,D} {3,D}
3 N u0 p1 c0 {2,D} {6,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.674453,0.0266535,-2.72826e-05,1.57426e-08,-3.7599e-12,21320.3,19.3645], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.08696,0.0127931,-6.17629e-06,1.46136e-09,-1.3698e-13,20648.3,2.95368], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1107""",
    longDesc = 
u"""
1107
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C=C=N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 114,
    label = "CH2CN",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 C u0 p0 c0 {1,S} {5,T}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 N u0 p1 c0 {2,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.63064,0.0173644,-1.70284e-05,9.86551e-09,-2.46034e-12,29579.2,11.2776], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.14874,0.006066,-2.17175e-06,3.4975e-10,-2.09004e-14,28649.1,-6.59236], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Radical     T01/03""",
    longDesc = 
u"""
Radical     T01/03.
[CH2]C#N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 115,
    label = "CHCNH",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {5,S}
3 N u1 p1 c0 {1,S} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.87301,0.0198561,-2.4152e-05,1.60559e-08,-4.26088e-12,45853.8,9.34467], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.54663,0.00765158,-3.57998e-06,8.31109e-10,-7.71667e-14,45394.6,-3.17646], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1107""",
    longDesc = 
u"""
1107
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.

E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C(#C)[NH]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 116,
    label = "CH3NH",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 N u1 p1 c0 {1,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.70974,-0.00731947,4.95106e-05,-5.7279e-08,2.16524e-11,20061.9,1.8546], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.02244,0.0103512,-3.6456e-06,5.80492e-10,-3.44104e-14,19505.1,1.64484], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""radical     T03/10""",
    longDesc = 
u"""
radical     T03/10
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
C[NH]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 117,
    label = "CH3NCH",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 N u0 p1 c0 {1,S} {3,D}
3 C u1 p0 c0 {2,D} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.23656,0.0158045,1.22159e-06,-9.86661e-09,4.14655e-12,25585.4,15.0702], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.2967,0.0169751,-8.65109e-06,2.12636e-09,-2.04493e-13,25102.8,8.60277], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1106""",
    longDesc = 
u"""
1106
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
CN=[CH]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 118,
    label = "CH2CHN",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,S} {3,D} {6,S}
2 C u1 p0 c0 {1,S} {4,S} {5,S}
3 N u1 p1 c0 {1,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.385115,0.0263156,-2.31525e-05,1.07917e-08,-1.99498e-12,41885.6,22.3684], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.14454,0.0131732,-6.28186e-06,1.4399e-09,-1.30842e-13,41038.9,3.75755], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1110""",
    longDesc = 
u"""
1110
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C(=[N])[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 119,
    label = "CHCNH2",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {3,T}
3 C u0 p0 c0 {2,T} {6,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.3481,0.0260893,-3.14512e-05,2.11896e-08,-5.70653e-12,26621.5,14.8874], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.67245,0.010499,-4.6263e-06,1.01346e-09,-8.92617e-14,26071.3,-0.577431], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1107""",
    longDesc = 
u"""
1107
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.

E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
NC#C
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 120,
    label = "CH2CHNH2",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,D} {4,S}
2 N u0 p1 c0 {1,S} {5,S} {6,S}
3 C u0 p0 c0 {1,D} {7,S} {8,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.710306,0.0374245,-3.66191e-05,1.98714e-08,-4.41925e-12,3580.17,25.7562], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.59217,0.0172596,-7.93921e-06,1.79645e-09,-1.61749e-13,2467.43,-0.08638], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1104""",
    longDesc = 
u"""
1104

E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C(=C)N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 121,
    label = "CH2CHN(S)",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,S} {3,D} {6,S}
2 C u1 p0 c0 {1,S} {4,S} {5,S}
3 N u1 p1 c0 {1,D}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.385115,0.0263156,-2.31525e-05,1.07917e-08,-1.99498e-12,49434,22.3684], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.14454,0.0131732,-6.28186e-06,1.4399e-09,-1.30842e-13,48587.4,3.75755], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1110""",
    longDesc = 
u"""
1110
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
Duplicate of species CH2CHN (i.e. same molecular structure according to RMG)
C(=[N])[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 122,
    label = "CH3NCH2",
    molecule = 
"""
1 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {3,D} {7,S} {8,S}
3 N u0 p1 c0 {1,S} {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.12431,0.0161642,4.57771e-06,-1.2116e-08,4.53506e-12,7335.42,14.5092], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[1.74096,0.0220935,-1.09102e-05,2.60524e-09,-2.44225e-13,7192.29,15.2596], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1104""",
    longDesc = 
u"""
1104
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
CN=C
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 123,
    label = "c-C2H3N",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {3,D} {6,S}
3 N u0 p1 c0 {1,S} {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.0828524,0.0220081,-1.21834e-05,1.00059e-09,1.0153e-12,31399.6,23.0893], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.80979,0.0150202,-7.58142e-06,1.85247e-09,-1.77642e-13,30658.2,8.95343], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1110""",
    longDesc = 
u"""
1110
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C1C=N1
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 124,
    label = "NCNOH",
    molecule = 
"""
multiplicity 2
1 N u1 p1 c0 {2,S} {3,S}
2 C u0 p0 c0 {1,S} {4,T}
3 O u0 p2 c0 {1,S} {5,S}
4 N u0 p1 c0 {2,T}
5 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.56418,-0.0166876,4.46094e-05,-4.38687e-08,1.49943e-11,25410.3,-7.23938], Tmin=(100,'K'), Tmax=(942.36,'K')),
            NASAPolynomial(coeffs=[3.4173,-0.000492619,8.01276e-07,-1.21784e-10,4.77703e-15,25123.5,-4.44293], Tmin=(942.36,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (100,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""Thermo group additivity estimation: group(N3s-CsHH) + radical(NHJ_C)""",
    longDesc = 
u"""
Thermo group additivity estimation: group(N3s-CsHH) + radical(NHJ_C).
[N](C#N)O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 125,
    label = "H2NCO",
    molecule = 
"""
multiplicity 2
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 C u1 p0 c0 {1,S} {5,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 O u0 p2 c0 {2,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.5364,0.00973407,-3.87293e-07,-5.90128e-09,3.01182e-12,-3096.24,8.47952], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.69169,0.00608718,-2.09434e-06,3.28449e-10,-1.92704e-14,-3810.29,-3.2271], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""H2N-C*=O  T09/09""",
    longDesc = 
u"""
H2N-C*=O  T09/09
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
N[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 126,
    label = "CH2NCH2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {3,D} {6,S} {7,S}
2 C u1 p0 c0 {3,S} {4,S} {5,S}
3 N u0 p1 c0 {1,D} {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.869061,0.0252858,-1.51436e-05,3.01037e-09,3.70027e-13,25778.2,19.1044], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.32005,0.0185257,-9.56913e-06,2.33792e-09,-2.22817e-13,25135.8,6.51881], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1104""",
    longDesc = 
u"""
1104
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C=N[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 127,
    label = "HNC",
    molecule = 
"""
1 N u0 p0 c+1 {2,S} {3,T}
2 H u0 p0 c0 {1,S}
3 C u0 p1 c-1 {1,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.36665,0.0149943,-3.03606e-05,2.99846e-08,-1.09095e-11,21980.8,7.87288], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.17272,0.0026657,-8.923e-07,1.37331e-10,-7.95546e-15,21797.7,0.215196], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

HNC <g> ATcT ver. 1.122, DHf298 = 192.39 ? 0.38 kJ/mol - fit JAN17.
[NH+]#[C-]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 128,
    label = "CH2NH2",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 N u0 p1 c0 {1,S} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.85538,0.00727364,1.65713e-05,-2.70977e-08,1.16328e-11,16616.6,10.2445], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.5533,0.00942003,-3.20781e-06,4.98962e-10,-2.90872e-14,15871.7,-0.0245945], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T 8/11""",
    longDesc = 
u"""
T 8/11
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013..
[CH2]N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 129,
    label = "CH3NHCH2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 N u0 p1 c0 {1,S} {3,S} {7,S}
3 C u1 p0 c0 {2,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.71402,0.0275573,-1.6131e-05,4.10455e-09,-8.07562e-14,16625.5,16.377], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.95751,0.0209452,-9.75571e-06,2.21626e-09,-1.99156e-13,16058.8,4.96297], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1104""",
    longDesc = 
u"""
1104
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
CN[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 130,
    label = "NCNO",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {3,D}
2 C u0 p0 c0 {1,S} {4,T}
3 O u0 p2 c0 {1,D}
4 N u0 p1 c0 {2,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.27607,0.0174736,-2.73909e-05,2.38161e-08,-8.27602e-12,37083.6,9.634], Tmin=(100,'K'), Tmax=(798.02,'K')),
            NASAPolynomial(coeffs=[5.02882,0.00670212,-3.41119e-06,6.64942e-10,-4.6398e-14,36867.1,1.96968], Tmin=(798.02,'K'), Tmax=(5000,'K')),
        ],
        Tmin = (100,'K'),
        Tmax = (5000,'K'),
    ),
    shortDesc = u"""thermo_DFT_CCSDTF12_BAC""",
    longDesc = 
u"""
thermo_DFT_CCSDTF12_BAC.
N(=O)C#N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 131,
    label = "CH3NCH3",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3 N u1 p1 c0 {1,S} {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.61165,0.00906653,2.27492e-05,-2.78392e-08,9.27267e-12,17306.4,11.2643], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[1.54968,0.0247836,-1.20303e-05,2.81641e-09,-2.5854e-13,17345.3,19.3448], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1104""",
    longDesc = 
u"""
1104
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C[N]C
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 132,
    label = "NCN",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,T} {3,S}
2 N u0 p1 c0 {1,T}
3 N u2 p1 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.00703,0.00838874,-6.90467e-06,2.84869e-09,-5.74273e-13,53004.5,7.78171], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.65712,0.00164211,-6.37923e-07,1.11209e-10,-6.99328e-15,52209.7,-6.1405], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed December 2016.

NCN <g> ATcT ver. 1.122, DHf298 = 450.80 ? 0.68 kJ/mol - fit JAN17.
C(#N)[N]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 133,
    label = "HNCN",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,T}
2 N u1 p1 c0 {1,S} {4,S}
3 N u0 p1 c0 {1,T}
4 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.06754,0.010679,-7.96224e-06,2.59883e-09,-1.27058e-13,36662.4,9.41075], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.53846,0.00389054,-1.38105e-06,2.21295e-10,-1.31827e-14,35963.5,-3.39587], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Cyanamide    T03/10""",
    longDesc = 
u"""
Cyanamide    T03/10.
C(#N)[NH]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 134,
    label = "HNCNH",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,D}
2 N u0 p1 c0 {1,D} {4,S}
3 N u0 p1 c0 {1,D} {5,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.00198,0.0246402,-2.75977e-05,1.53247e-08,-3.26828e-12,16793.7,16.9432], Tmin=(298,'K'), Tmax=(1500,'K')),
            NASAPolynomial(coeffs=[8.37414,0.00236614,-3.50232e-07,-4.3911e-11,1.09686e-14,14610.9,-21.0739], Tmin=(1500,'K'), Tmax=(4000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (4000,'K'),
    ),
    shortDesc = u"""62790""",
    longDesc = 
u"""
62790
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.


********************************************************************************
********************************************************************************
***** HNCN NCN NCCN
********************************************************************************
********************************************************************************
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C(=N)=N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 135,
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
            NASAPolynomial(coeffs=[3.78898,0.00591795,-1.08086e-06,-1.85638e-09,1.06887e-12,48196,3.71766], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.47848,0.00484708,-1.76213e-06,2.92788e-10,-1.76554e-14,47949,-0.0784806], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT3E""",
    longDesc = 
u"""
ATcT3E

CCH2 <g> ATcT ver. 1.122, DHf298 = 412.20 ? 0.33 kJ/mol - fit JAN17.
[C]=C
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 136,
    label = "CHCHO",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,S} {3,D} {4,S}
2 C u2 p0 c0 {1,S} {5,S}
3 O u0 p2 c0 {1,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.06864,0.0187233,-1.21319e-05,-3.33727e-10,2.32882e-12,29739.4,14.7866], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.96288,0.00799899,-4.30606e-06,1.11076e-09,-1.11415e-13,28725.6,-5.17392], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical""",
    longDesc = 
u"""
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016.
[CH]C=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 137,
    label = "NH2OH",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 O u0 p2 c0 {1,S} {5,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.21016,0.00619672,1.10595e-05,-1.96668e-08,8.82517e-12,-6581.48,7.93294], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.88113,0.00815708,-2.82616e-06,4.37931e-10,-2.52725e-14,-6860.18,3.79156], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ATcT/A""",
    longDesc = 
u"""
ATcT/A
Klippenstein et al, 2011 paper

********************************************************************************
********************************************************************************
***** NH2OH H2NO HNOH HNO HON NO HONO HNO2 NO2 HONO2 NO3 N2O
********************************************************************************
********************************************************************************.
NO
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 138,
    label = "NCCN",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,S} {3,T}
2 C u1 p0 c0 {1,S} {4,D}
3 N u0 p1 c0 {1,T}
4 N u1 p1 c0 {2,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.32928,0.0261541,-4.9001e-05,4.61923e-08,-1.64326e-11,35690.1,9.86348], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.7055,0.00364271,-1.3094e-06,2.16421e-10,-1.31194e-14,34882.4,-10.4803], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Dicyanogen   ATcT/A""",
    longDesc = 
u"""
Dicyanogen   ATcT/A.
C(#N)[C]=[N]
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 139,
    label = "CH2CHNO",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,S} {6,S}
3 N u0 p1 c0 {1,S} {7,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 O u0 p2 c0 {3,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.782134,0.0286251,-2.14511e-05,7.40944e-09,-7.1095e-13,18547.8,21.0235], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.48741,0.0170182,-8.86195e-06,2.23083e-09,-2.19836e-13,17646,2.34417], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""0""",
    longDesc = 
u"""
0

pm 0910
Low T polynomial Tmin changed from 350.0 to 298.0 K when importing to RMG.
C(=C)N=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 140,
    label = "CHCHNO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u1 p0 c0 {1,D} {5,S}
3 N u0 p1 c0 {1,S} {6,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 O u0 p2 c0 {3,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.66767,0.0279986,-2.86416e-05,1.53555e-08,-3.31674e-12,48996,17.1916], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.91529,0.012105,-6.44647e-06,1.65557e-09,-1.65983e-13,48091.6,-3.57487], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""0""",
    longDesc = 
u"""
0

pm 1002
Low T polynomial Tmin changed from 350.0 to 298.0 K when importing to RMG.
C(=[CH])N=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 141,
    label = "CH3NH2",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 N u0 p1 c0 {1,S} {6,S} {7,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.65718,-0.00737953,5.45297e-05,-6.24994e-08,2.33762e-11,-3760.73,1.63875], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.7861,0.012697,-4.46779e-06,7.11022e-10,-4.21333e-14,-4381,1.86269], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""T09/09""",
    longDesc = 
u"""
T09/09


********************************************************************************
********************************************************************************
***** CH3NH2 CH2NH2 CH3NH CH2NH H2CN HCNH
********************************************************************************
********************************************************************************.
CN
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 142,
    label = "CH3CH2NH2",
    molecule = 
"""
1  C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2  C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  N u0 p1 c0 {1,S} {9,S} {10,S}
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
            NASAPolynomial(coeffs=[0.64408,0.0314867,-1.71052e-05,3.57316e-09,1.95649e-13,-7471.2,21.3696], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.82805,0.025375,-1.1874e-05,2.70991e-09,-2.44591e-13,-8039.21,10.1772], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1104""",
    longDesc = 
u"""
1104
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.

E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.



********************************************************************************
********************************************************************************
***** CH3CH2NH2 CH2CH2NH2 CH3CHNH2 CH3CH2NH CH2CHNH2 CH3CHNH CH2CHNH CH2CNH2
***** CH3CNH CH3CHN CHCHNH2
********************************************************************************
********************************************************************************
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C(C)N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 143,
    label = "CH2CH2NH2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u1 p0 c0 {1,S} {6,S} {7,S}
3 N u0 p1 c0 {1,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.25929,0.0301867,-2.18269e-05,8.94597e-09,-1.52898e-12,17667,19.6477], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.10303,0.0203085,-9.25456e-06,2.06112e-09,-1.8195e-13,17023.4,5.55411], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1104""",
    longDesc = 
u"""
1104

E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C([CH2])N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 144,
    label = "CH3CHNH2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 C u1 p0 c0 {1,S} {3,S} {7,S}
3 N u0 p1 c0 {2,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.64075,0.0285827,-1.89469e-05,6.7295e-09,-9.08174e-13,12515.5,17.4881], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.08789,0.0205228,-9.45004e-06,2.12695e-09,-1.89706e-13,11939.7,5.24988], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1104""",
    longDesc = 
u"""
1104
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C[CH]N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 145,
    label = "CH2CNH2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {3,D} {4,S} {5,S}
2 N u0 p1 c0 {3,S} {6,S} {7,S}
3 C u1 p0 c0 {1,D} {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.01233,0.0256068,-2.42451e-05,1.39182e-08,-3.46136e-12,30777.7,15.07], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.21874,0.0156981,-7.7573e-06,1.84326e-09,-1.71838e-13,30390.6,4.69612], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1107""",
    longDesc = 
u"""
1107
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C=[C]N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 146,
    label = "CHCHNH2",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {4,S}
2 N u0 p1 c0 {1,S} {5,S} {6,S}
3 C u1 p0 c0 {1,D} {7,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.583298,0.0310085,-3.08423e-05,1.72291e-08,-3.96983e-12,36294.5,20.6175], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.77535,0.0146042,-6.78153e-06,1.55237e-09,-1.41604e-13,35437.9,0.30232], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1104""",
    longDesc = 
u"""
1104
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
C(=[CH])N
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 147,
    label = "H2NCHO",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {3,D} {6,S}
3 O u0 p2 c0 {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.17952,0.0019459,3.4333e-05,-4.64826e-08,1.87268e-11,-24058,9.94755], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.93695,0.00960834,-3.32935e-06,5.16598e-10,-2.99203e-14,-25091,-2.00078], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""O=CH-NH2   T12/09""",
    longDesc = 
u"""
O=CH-NH2   T12/09.
NC=O
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

entry(
    index = 148,
    label = "CH3NHCH3",
    molecule = 
"""
1  C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3  N u0 p1 c0 {1,S} {2,S} {10,S}
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
            NASAPolynomial(coeffs=[1.22882,0.0260431,-4.95078e-06,-6.74954e-09,3.27739e-12,-3768.13,18.4832], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.13913,0.0268643,-1.28763e-05,2.99663e-09,-2.74772e-13,-4173.32,12.9759], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1104""",
    longDesc = 
u"""
1104
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
tables (ftp://ftp.technion.ac.il/pub/supported/aetdd/thermodynamics; mirrored at
http://garfield.chem.elte.hu/burcat/burcat.html. Accessed July 2013.


********************************************************************************
********************************************************************************
***** CH3NHCH3 CH3NHCH2 CH3NCH3 CH3NCH2 CH2NCH2 CH3NCH
********************************************************************************
********************************************************************************
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
CNC
Imported from /home/alongd/Code/RMG-Py/importer/NOx2018/thermo.txt.
""",
)

