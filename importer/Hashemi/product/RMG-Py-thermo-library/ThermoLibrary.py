#!/usr/bin/env python
# encoding: utf-8

name = "/home/alongd/Code/RMG-Py/importer/Hashemi"
shortDesc = u"/home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT"
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
    shortDesc = u"""DOROFEEVA e  T 8/03""",
    longDesc = 
u"""
DOROFEEVA e  T 8/03
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
7722-84-1
H2O2 Hydrogen Peroxide  SIGMA=2  IAIBIC=2.976026 E-117 Ir=0.071093 NU=3599,
1388,875,3611,1266  X11=90.9 x12=11 x13=11 x15=167.6 x16=3 x22=10 x23=7
x25=11.5 x26=4  x33=10 x35=11.1 x36=2 x55=90.2 x56=3 x66=3  ALFAA1=0.234
ALFAA2=-0.104  ALFAA3=-0.034  ALFAA5=0.234 ALFAA6=-0.174  ALFAB1=0.003
ALFAB2=0.004  ALFAB3=.007 ALFAB5=.003 ALFAB6=.003 ALFAC1=-.001 ALFAC2=.004
ALFAC3=.013 ALFAC5=-.001 ALFAC6=.007  HF298= -135.88+/-0.2 kJ HF0=-129.9 kJ
REF=Dorofeeva et al JPCRD 32 (2003),879  {HF298=-135.77+/-0.07  kJ REF=ATcT A;
HF298=-135.44+/-0.07 kJ  REF=ATcT B}  Max Lst Sq Error CP @ 6000 K 0.29%
Calculated directly from Original tables.
OO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
    shortDesc = u"""HYDROXYL RADI  IU3/03""",
    longDesc = 
u"""
HYDROXYL RADI  IU3/03
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
3352-57-6
OH  HYDROXYL RADICAL T0=0 STATWT=2 WE=3738   WEXE=-84.88  WEYE=0.54 B0=18.550
T0=139.7 STATWT=2  T0=32403 STATWT=2 WE=3184.3 WEXE=97.845 BE=17.355  REF=Poly-
nomials were calculated directly from Gurvich's Tables wich are direct summation,
HF298=37.3+/-0.3 kJ  HF0=37.1+/-0.3 kJ   REF=Ruscic et al JPC 106,(2002),2727.
Ruscic et al JPCRD 2005  {HF298=37.5+/-0.03 kJ  REF=ATcT B}   Max Lst Sq Error
Cp @ 1200 K 0.28%.
[OH]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
    shortDesc = u"""REF ELEMENT    tpis78""",
    longDesc = 
u"""
REF ELEMENT    tpis78
1333-74-0
H2   Calc from Gurvic's table  HF298= 0.0  REF=Gurvich 78  Max Lst Sq Error Cp
@ 400 & 1300 K 0.33%.
[H][H]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
    shortDesc = u"""""",
    longDesc = 
u"""
[H]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
L 1/90
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
17778-80-2
O  HF298=249.175+/-0.1 KJ  REF=C.E. Moore "Selected Tables of Atomic Spectra"
NSRDS-NBS Sec 7 1976 p. A8 I {HF298=249.229+/-0.002 kJ  REF=ATcT B}.
[O]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
L 5/89
7732-18-5
H2O REF=Wooley J. RES. NBS 92 (1987), 35. HF298=-241.826+/-0.04 KJ based on
HF298(L) from Cox, Wagman & Medvedev CODATA KEY VAL. for THERMO, Hemisphere 1989
p.21 and heat of vap. from Haar, Gallagher & Kell NBS/NRC Tables, Hemisphere
1984.   {HF298=-241.822+/-0.03 kJ  REF=ATcT A}..
O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
T 1/09
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
3170-83-0
HO2 RADICAL  GROUND STATE STATWT=2  A0=20.356 B0=1.118 C0=1.056 NU=3436.2,
1391.75,1097.63 EXCITED STATE T0=7029.684 A0=20.486  B0=1.021 C0=0.968
NU=3268.5,1285,929.068 STATWT=2 REF=JACOX JPCRD 17, (1988) P.269   HF298=12.296
+/-0.25  REF=ATcT A   {HF298=12.27+/-0.16  kJ  REF=ATcT B;  HF298=12.55 kJ
REF=Hills & Howard J. CHEM. PHYS 81,(1984),4458;  HF298=-12.41 +/-0.58 kJ
REF=Karton, Parthiban, Martin JPC A 2009}  Max lst Sq Error Cp @ 700 K 0.27%.
[O]O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
    shortDesc = u"""REF ELEMENT    TPIS89""",
    longDesc = 
u"""
REF ELEMENT    TPIS89
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
7782-44-7
O2  CALCULATED FROM ORIGINAL VALUES  HF298=0 KJ  REF=Gurvich 1989. Corrected by
B.McBride NASA TP-2002-211556 Max Lst Sq Error Cp @ 1200 K 0.31%.
[O][O]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
    shortDesc = u"""REF ELEMENT    g 5/97""",
    longDesc = 
u"""
REF ELEMENT    g 5/97
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
7440-37-1
Ar  HF298=0.  REF=C.E. Moore "Atomic Energy Levels" NSRDS-NBS 35 (1971) p.211
Max Lst Sq Error 0.%.
[Ar]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
    shortDesc = u"""REF ELEMENT   G 8/02""",
    longDesc = 
u"""
REF ELEMENT   G 8/02
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
7727-37-9
N2  REFERENCE ELEMENT HF=0. from TSIV Tables  Max Lst Sq Error Cp @ 6000 K 0.29%.
N#N
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
    shortDesc = u"""REF ELEMENT    g 5/97""",
    longDesc = 
u"""
REF ELEMENT    g 5/97
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
7440-59-7
He  HF298=0.0 KJ  Atomic Energy levels Moore  NSRDS-NBS 35 1971  Max Lst Sq
Error Cp 0.0%.
[He]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
            NASAPolynomial(coeffs=[3.57953,-0.000610354,1.01681e-06,9.07006e-10,-9.04424e-13,-14344.1,3.50841], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.04849,0.00135173,-4.85794e-07,7.88536e-11,-4.69807e-15,-14266.1,6.0171], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""RUS 79""",
    longDesc = 
u"""
RUS 79
74-82-8
CH4 METHANE  STATWT=1. SIGMA=12.  IA=IB=IC=0.52410356  NU=2916.7,1533.295(2),
3019.491(3),1310.756(3)  X11=-26 X12=-3 X13=-75 X14=-4 X22=-.4,X23=-9 X24=-20
X33=-17 X34=-17 X44=-11  ALFA1=.01 ALFA2=-.09 ALFA3=.04 ALFA4=.07 D0=1.10864E-4
HF298=-74.6+/-0.3 KJ HF0=66.63 kJ REF=TSIV 91 {HF298=-74.554+/-0.60 kJ  REF=ATcT
C}   MAX LST SQ ERROR CP @ 1300 K 0.54%.
CH4   ANHARMONIC  g 8/99C  1.H  4.   0.   0.G   200.000  6000.000  1000.000,    1
1.65326226E+00 1.00263099E-02-3.31661238E-06 5.36483138E-10-3.14696758E-14    2
-1.00095936E+04 9.90506283E+00 5.14911468E+00-1.36622009E-02 4.91453921E-05    3
-4.84246767E-08 1.66603441E-11-1.02465983E+04-4.63848842E+00-8.97226656E+03    4
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
630-08-0
CO  CARBON-MONOXIDE  CALCULATED FROM TSIV TABLE. REF=TSIV 79  HF298=-110.53+/-
0.17 kJ {HF298=-110.53+/-0.026  REF=ATcT C} Max Lst Sq Error Cp @ 1300 K 0.12%..
[C-]#[O+]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
            NASAPolynomial(coeffs=[2.35681,0.00898413,-7.12206e-06,2.4573e-09,-1.42885e-13,-48372,9.9009], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.63651,0.00274146,-9.95898e-07,1.60387e-10,-9.16199e-15,-49024.9,-1.9349], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""L 7/88""",
    longDesc = 
u"""
L 7/88
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
124-38-9
CO2  CARBON-DIOXIDE   SIGMA=2   B0=0.39027   NU=1333.5,667(2),2351   X11=-3.014
X12=-5.058   X12=-19.048   X22=1.521   X23=-12.616   X33=-12.597   G22=-1.422
Y111=.0184   Y112=-.0667   Y113=-.0944    Y122=-.0657  Y123=.0880   Y133=.0268
Y222=.0105   Y223=-.0168   Y233=.0320   Y333=.0115   W0=51.834   ALPHA1=.00115
ALPHA2=-.000715   ALPHA3=.00311   D000=.129E-6   T0=30000 STATWT=3;   T0=33000
STATWT=6  T0=36000  STATWT=3;  T0=45000  STATWT=2;  REF=Gurvich Vol 2 1991 p.27
HF298=-393.51 KJ  {HF298=-393.472+/-0.014 kJ  REF=ATcT A}   Max Lst Sq Error Cp
@ 1400 K 0.4%.
O=C=O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
            NASAPolynomial(coeffs=[5.14826,-0.0137002,4.93749e-05,-4.91952e-08,1.70097e-11,-10245.3,-4.63323], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[1.91179,0.00960268,-3.38388e-06,5.38797e-10,-3.19307e-14,-10099.2,8.48242], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""RRHO        g 8/99""",
    longDesc = 
u"""
RRHO        g 8/99
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
74-82-8
CH4  METHANE  Same as the Anharmonic but calculated Using the RRHO method rather
than the NRRAO2.  Max Lst Sq Error Cp @ 6000. K 0.62%..
C
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
            NASAPolynomial(coeffs=[3.65718,0.0021266,5.45839e-06,-6.6181e-09,2.46571e-12,16422.7,1.67354], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.97812,0.00579785,-1.97558e-06,3.07298e-10,-1.79174e-14,16509.5,4.72248], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""IU0702""",
    longDesc = 
u"""
IU0702
DHf(0K) =   102.62 [kcal/mol], taken from SJK v0.9.
Q(T) from B3LYP/6-311++G(d,p) by CFG on  06Sep2013
2229-07-4
CH3 METHYL-RAD STATWT=1. SIGMA=6. IA=IB=.2923  IC=.5846  NU=3004,606.4,3161(2),
1396(2)  HF298=146.7 +/-0.3 KJ HF0=150.0+/-0.3 kJ REF= Ruscic et al JPCRD 2003.
{HF298=146.5+/-0.08 kJ  REF=ATcT C}   Max Lst Sq Error Cp @ 6000 K 0.44%.
METHYL RADICAL    IU0702C  1.H  3.   0.   0.G   200.000  6000.000  B  15.03452 1.
[CH3]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
            NASAPolynomial(coeffs=[5.65851,-0.0162983,6.91938e-05,-7.58373e-08,2.80428e-11,-25612,-0.897331], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.52727,0.0103179,-3.62893e-06,5.77448e-10,-3.42183e-14,-26002.9,5.16759], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Methyl alc  T06/02""",
    longDesc = 
u"""
Methyl alc  T06/02
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
67-56-1
CH4O  METHANOL (CH3OH) STATWT=1.  SIGMA=1. IA=.6578  IB=3.4004  IC=3.5306
Brot=28.182 ROSYM=3  V3=373.21  V6=-0.521 cm-1  NU=3681,3000,2844,1477,1455,
1345,1060,1033,2960,1477,1165    HF298=-201. KJ  REF=CHEN WILHOIT &  ZWOLINSKI
JPCRD 6,(1977),105  {HF298=-200.70+/-0.17 kJ REF=ATcT C}  MAX LST SQ ERROR Cp
@ 1300 K  0.82%..
CO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
            NASAPolynomial(coeffs=[2.49975,5.02139e-07,1.93093e-09,-4.94633e-12,2.74089e-15,85485.8,3.66241], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.49891,2.41252e-06,-1.88942e-09,6.24309e-13,-7.45464e-17,85486.2,3.6671], Tmin=(1000,'K'), Tmax=(3000,'K')),
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
Accessed April 2016
<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C1 species:
<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
7440-44-0
C  Amorphous Carbon, Acetylene black, Lamp black   HF298=716.68+/-0.45 kJ
REF=C.E. Moore "Selected Tables of Atomic Spectra"  NSRDS-NBS Sec 3 (1970)
p A6 I. {HF298=716.87+/-0.06 kJ  REF=ATcT B; HF0=711.53 kJ ATcT rev value as
quoted by Karton Tarnopolsky and Martin Mol Phys 2009}
C                 L 7/88C   1    0    0    0G   200.000  6000.000  1000.000,    1
0.26055830E+01-0.19593434E-03 0.10673722E-06-0.16423940E-10 0.81870580E-15    2
0.85411742E+05 0.41923868E+01 0.25542395E+01-0.32153772E-03 0.73379223E-06    3
-0.73223487E-09 0.26652144E-12 0.85442681E+05 0.45313085E+01 0.86195097E+05    4
! E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
! http://burcat.technion.ac.il/dir/
! mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
! Accessed April 2016.
[C]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
            NASAPolynomial(coeffs=[3.47159,0.000404745,-1.89388e-06,3.3321e-09,-1.51271e-12,70701.8,2.14542], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.66195,0.00170245,-6.83824e-07,1.3098e-10,-9.71843e-15,70958.9,6.52673], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""DHf(0K) =   170.13 [kcal/mol], taken from SJK v0.9.""",
    longDesc = 
u"""
DHf(0K) =   170.13 [kcal/mol], taken from SJK v0.9.
Q(T) from B3LYP/6-311++G(d,p) by CFG on  06Sep2013
3315-37-5
CH  METHYLIDYNE RADICAL CALCULATED FROM Gurvich's Tables IB=01973   We=2732.46
HF298=595.8+/-0.6 kJ HF0=592.5+/-0.6 kJ  REF= Ruscic et al JPCRD 2005
{HF298=596.17+/-0.12 kJ  REF-ATcT A}
CH               IU3/03 C  1.H  1.   0.   0.G   200.000  6000.000  1000.000,    1
0.25209369E+01 0.17653639E-02-0.46147660E-06 0.59289675E-10-0.33474501E-14    2
0.70946769E+05 0.74051829E+01 0.34897583E+01 0.32432160E-03-0.16899751E-05    3
0.31628420E-08-0.14061803E-11 0.70612646E+05 0.20842841E+01 0.71658188E+05    4
! E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
! http://burcat.technion.ac.il/dir/
! mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
! Accessed April 2016.
[CH]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 19,
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
Ethynyl Rad   T 5/10

12070-15-4
C2 triplet  T0=716.2  SIGMA=2  STATWT=3  Be=1.6324  We=1641.35  WeXe=11.67
ALPHAE=0.01661  De=6.44E-6  REF=Hubert & Herzberg Webbook and Gurvich 91
HF0=827.26+/-8. kJ  REF=Karton Tarnopolsky Martin Mol Phys 107,(2009),977
Max Lst Sq Error Cp @ 1300 K 0.30%
C2 triplet        T05/09C  2.   0.   0.   0.G   200.000  6000.000  B  24.02140 1
3.43350371E+00 1.07185010E-03-3.97897382E-07 6.67457391E-11-4.10152154E-15    2
1.00178987E+05 4.10588356E+00 3.76163273E+00-2.72143299E-03 8.69879462E-06    3
-8.19304667E-09 2.62415296E-12 1.00254566E+05 3.18038623E+00 1.01317039E+05    4
2122-48-7
C2H ETHYNYL RADICAL SIGMA=1  STATWT=2    B0=1.457 NU=3328,372(2),1841  T0=4000
STATWT=4  B0=1.457 NU=3460,560(2),1850  REF=Kiefer, Sidhu, Kern, Xie,Chen,
Harding 1992  HF298=568.056+/-0.3 kJ  REF=ATcT A {HF298=567.9+/-0.2 kJ
REF=ATcT C 2011;  HF298=568.522+/-4 kJ REF= NIST Webbook 1999. HF298=567.4+/-1.5
kJ  REF=Szalay Tajti & Stanton Mol Phys 103,(2005),xxx}   MAX LST SQ ERROR Cp @
400 K 0.34 % ..
[C]#C
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 20,
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
            NASAPolynomial(coeffs=[3.9592,-0.00757051,5.7099e-05,-6.91588e-08,2.69884e-11,5089.78,4.0973], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.99183,0.0104834,-3.71721e-06,5.94628e-10,-3.5363e-14,4268.66,-0.269082], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""g 1/00""",
    longDesc = 
u"""
g 1/00
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
74-85-1
C2H4  ETHYLENE STATWT=1.  SIGMA=4.  A0=4.86596  B0=1.001329  C0=0.828424
NU=3021,1625,1344,1026,3083,1222,949,940,3105,826,2989,1444   REF=CHAO &
ZWOLINSKY, JPCRD 4,(1975),251  HF298=52.5 kJ HF0=61.025 kJ   REF=TRC 4/1988
{HF298=52.55 +/-0.14 kJ  REF=ATcT C}  MAX LST SQ ERROR Cp 20K 0.80 ..
C=C
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 21,
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
            NASAPolynomial(coeffs=[4.24186,-0.00356905,4.82667e-05,-5.85401e-08,2.25805e-11,12969,4.44704], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.32196,0.0123931,-4.39681e-06,7.0352e-10,-4.18435e-14,12175.9,0.171104], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""ethyl radic  IU1/07""",
    longDesc = 
u"""
ethyl radic  IU1/07
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
2025-56-1
C2H5  ETHYL RAD.  STATWT=2  SIGMA=1  IA=0.2917  IB=3.6874  IC=3.9786  IR=0.1812
ROSYM=6.   V6=16.6 cm-1  NU=3128.7,3037,2987,2920,2842,1440(2),1425,1366,1175,
1138,958,780,528.1  HF298=119.7+/-0.7 kJ  REF=Ruscic IUPAC Task Group
{HF298=119.9+/-0.4 kJ  REF=ATcT C 2011;  HF298=28.36 Kcal. REF= Chen, Rauk &
Tschuikow-Roux (1990)}     Max Lst Sq Error Cp & 6000 K 0.58%.
C[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 22,
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
    shortDesc = u"""g 8/88""",
    longDesc = 
u"""
g 8/88
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
74-84-0
C2H6  ETHANE  STATWT=1.  SIGMA=6.  IA=1.0481  IB=IC=4.22486 Ir=.26203 ROSYM=3
V0=2.96 kcal   NU=2954,1388,995,2896,1379,2969(2),1468(2),1190(2),2985(2),
1469(2),822(2) HF298=-83.863 kJ  REF=CHAO WILHOIT & ZWOLINSKI JPCRD 2,(1973),
427   {HF298=-83.77 +/-0.16 kJ  REF=ATcT C}  MAX LST SQ ERROR Cp @ 6000K 0.63%..
CC
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
            NASAPolynomial(coeffs=[3.71181,-0.00280463,3.76551e-05,-4.73072e-08,1.86588e-11,1295.7,6.57241], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.75779,0.00744142,-2.69705e-06,4.38091e-10,-2.63537e-14,378.112,-1.9668], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""METHOXY RA  IU1/03""",
    longDesc = 
u"""
METHOXY RA  IU1/03
2143-68-2
CH3O  METHOXY RADICAL  SYMNO=3.  A0=5.2  B0=C0=0.93  T0 STATWT=3 NU=2840,1417,
1047,2774(2),1465,1210,914,653 T0=61.97 STATWT=1  Specific calculations perfo-
rmed and the polynomials were calculated from tabular values obtained.
HF298=21.0+/-2.1 kJ   HF0=28.4 +/-2.1 kJ  REF=Ruscic et al IUPAC Group JPCRD
2003  {HF298=21.85+/-0.4 kJ  REF=ATcT C}   Max Lst Sq Error Cp @ 200 K 0.93%..
C[O]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
            NASAPolynomial(coeffs=[4.47834,-0.0013507,2.78485e-05,-3.64869e-08,1.47907e-11,-3500.73,3.30913], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.09314,0.00594761,-2.06497e-06,3.23008e-10,-1.88126e-14,-4034.1,-1.84691], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""RADICAL     IU2/03""",
    longDesc = 
u"""
RADICAL     IU2/03
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
2597-43-5
H3CO  HYDROXYMETHYLENE RAD (CH2OH)  STATWT=2.  SIGMA=1.  IA=.4274  IB=2.789
IC=3.2164    NU=3650,3169,3071,1459,1334,1176,1048,420,234   HF298=-17.0+/-0.7
kJ.  HF0=-10.7+/-0.7  Polynomials calculated from original Tables of Johnson &
Hudgens JPC 100 (1996),19874  extrapolated to 6000 K    REF=Ruscic et al JPCRD
2003 IUPAC Group   {HF298=-16.1+/-0.4 kJ  REF=ATcT C}   Max Lst Sq Error Cp @
200 & 6000 K 0.38%.
[CH2]O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 25,
    label = "CH2(S)",
    molecule = 
"""
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.19195,-0.00230793,8.0509e-06,-6.60123e-09,1.95638e-12,50484.3,-0.754589], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.28556,0.00460255,-1.97412e-06,4.09548e-10,-3.34695e-14,50922.4,8.67684], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""DHf(0K) =    93.49 [kcal/mol], taken from SJK v0.9.""",
    longDesc = 
u"""
DHf(0K) =    93.49 [kcal/mol], taken from SJK v0.9.
Q(T) from B3LYP/6-311++G(d,p) by CFG on  06Sep2013
2465-56-7
CH2 METHYLENE RADICAL SINGLET   SIGMA=2  STATWT=1   T0=0  IA=0.1391  IB=0.2498
IC=0.3960   NU=2806,1353,2865.   T0(b 1B1)=8350.  NU=3000,570,3000  A0=73.8
B0=8.59  C0=7.2   HF298=428.8+/-1.6 kJ  HF0=428.3+/-1.6 kJ  REF=Ruscic et al
JPCRD IUPAC Task Group 2003. {HF298=429.04+/-0.14 kJ REF=ATcT A}  Max Lst Sq
Error Cp @ 1300 K 0.40%
CH2(S) SINGLET    IU6/03C  1.H  2.   0.   0.G   200.000  6000.000  1000.000,    1
3.13501686E+00 2.89593926E-03-8.16668090E-07 1.13572697E-10-6.36262835E-15    2
5.05040504E+04 4.06030621E+00 4.19331325E+00-2.33105184E-03 8.15676451E-06    3
-6.62985981E-09 1.93233199E-12 5.03662246E+04-7.46734310E-01 5.15727280E+04    4
! E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
! http://burcat.technion.ac.il/dir/
! mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
! Accessed April 2016.
[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
            NASAPolynomial(coeffs=[3.36378,0.000265766,2.79621e-05,-3.72987e-08,1.5159e-11,34475,7.9151], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.15027,0.00754021,-2.62998e-06,4.15974e-10,-2.45408e-14,33856.6,1.72812], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Vinyl Radi  ATcT/A""",
    longDesc = 
u"""
Vinyl Radi  ATcT/A
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
2669-89-8
C2H3  VINYL-RAD  STATWT=2.  SIGMA=1.  A0=7.49   B0=1.07   C0=0.93  Nu=3265,3190,
3115,1670,1445,1185,920,825,785  REF=Ervin JACS 112 (1990),5750}  HF298=296.58
+/-0.92 kJ  HF0=300.867 kJ  REF=ATcT A   {HF298=297.272+/-0.5  kJ  REF=ATcT C;
HF298=299.74+/-5 kJ  REF=Ervin JACS 112,(1990),5750; also Kromkin Chimicheskaya
Fizika 22,(2002),30;   HF298=295.4+/-1.7 kJ  REF=Russell & Gutman JPC 93,(1989),
5184 also  Kaiser & Wallington JPC 100,(1996),4111 also Parthiban & Martin JCP
114,(2001),6014; HF298=299.6+/-3 kJ  REF=Tsang Energetics of Organic Free Rad
1996;  HF298=297.1+/-4.2  REF=De Moore et al JPL 97-4 1997}  Max Lst Sq Error
Cp @ 400 K 0.54%..
[CH]=C
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 27,
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
    shortDesc = u"""Formaldehyde T 5/11""",
    longDesc = 
u"""
Formaldehyde T 5/11
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
50-00-0
CH2O  FORMALDEHYDE  SIGMA=2 A0=9.40546 B0=1.295407 C0=1.134216 NU=2782.4,1746.1,
1500.1,1167.2,2843.2,1249.1   HF298=-109.164+/-0.1 kJ     REF=ATcT C 2011.
{HF298=-108.58 kJ  REF=Gurvich 89;  HF298=-112.5+/-8. kJ  REF=Burcat G3B3}  Max
Lst Sq Error Cp @ 1300 K 0.61%.
HCHO Formaldehyde T 5/11H  2.C  1.O  1.   0.G   200.000  6000.000  B  30.02598 1.
C=O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
            NASAPolynomial(coeffs=[0.80868,0.0233616,-3.55172e-05,2.80153e-08,-8.50075e-12,26429,13.9397], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.65878,0.00488397,-1.60829e-06,2.46975e-10,-1.38606e-14,25759.4,-3.99838], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""g 1/91""",
    longDesc = 
u"""
g 1/91
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
74-86-2
C2H2 ACETYLENE SIGMA=2 B0=1.1766 NU=3372.83,1973.8,3283.83,612.88(2),730.29(2)
X11=-18.57,X12=-13.09 X13=-102.39 X14=-16.54 X15=-10.85 X22=-7.92 X23=-2.83
X24=-12.70 X25=-1.38 X33=-30.95 X34=-8.22 X35=-8.68 X44=3.3 X45=-5.24 X55=-2.27
G44=-1.36 G55=3.45 ALPHA1=6.83E-3 ALPHA2=6.3E-3 ALPHA3=5.6E-3 ALPHA4=-1.3E-3
ALPHA5=-2.2E-3 D0=1.598E-6 T0=25000(3),35000(6),42198(1),50000(3),54116(1)
REF=TSIV HF298=228.2+/-0.8 kJ  HF0=228.769 {HF298=228.32+/-0.15 kJ  REF=ATcT C}
Max Lst Sq Error Cp @ 6000 K 0.24%
C2H2,acetylene    g 1/91C  2.H  2.   0.   0.G   200.000  6000.000  A  26.03728 1.
C#C
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 29,
    label = "HCO",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,D} {3,S}
2 O u0 p2 c0 {1,D}
3 H u0 p0 c0 {1,S}
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
T 5/03
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
2597-44-6
CHO SIGMA=1 STATWT=2 A0=24.562 B0=1.498 C0=1.403 NU=2435,1878,1087 REF=Marenich
& Boggs JPC 107 (2003),2343-2350 Direct summation using CCSD(T) method. Calc.
from their tables  HF298=42.3+/-2.0 kJ   HF0=41.9 kJ  {HF298=41.83+/-0.1 kJ
REF=ATcT C;  HF298=38.8+/-3.3 kJ  G4  REF=Marochkin & Dorofeeva Comp Theor Chem
991,(2012},182}   Max Lst Sq Error Cp @ 1500 K 0.63%.
CHO               T 5/03C  1.H  1.O  1.   0.G   200.000  6000.000  A  29.01804 1.
[CH]=O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 30,
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
            NASAPolynomial(coeffs=[3.8328,0.000224446,4.68033e-06,-6.04743e-09,2.59009e-12,45920.8,1.40666], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.16229,0.00281798,-7.56235e-07,5.05446e-11,5.65236e-15,46099.1,4.77656], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""DHf(0K) =   141.76 [kcal/mol], taken from SJK v0.9.""",
    longDesc = 
u"""
DHf(0K) =   141.76 [kcal/mol], taken from SJK v0.9.
Q(T) from B3LYP/6-311++G(d,p) by CFG on  06Sep2013
2465-56-7
CH2 METHYLENE RAD TRIPLET This is for applications where triplet methylene is
not equilibrated with single methylene. Only HF0 is identical to the one given
for the equilibrium Singlet and Triplet. HF298 is o,oo5 kJ lower than given for
the equilibrium HF1000=1.1 kJ lower and HF3000 is 6.8 kJ lower than the
equilibrium    STATWT=3 SIGMA=2   A0=73.811  B0=8.450  C0=7.184   NU=3031,963,
3190  T0=3500  SIGMA=2   HF298=391.2+/-1.6   HF0=390.7 +/-1.6 kJ  REF=Ruscic et
al JPCRD 2003 IUPAC Task Group {HF298=391.46+/-0.13  REF=ATcT C}  Max Lst Sq
Error Cp @ 6000 K 0.27%
CH2 TRIPLET RAD   IU3/03C  1.H  2.   0.   0.G   200.000  6000.000  1000.000,    1
3.14631886E+00 3.03671259E-03-9.96474439E-07 1.50483580E-10-8.57335515E-15    2
4.60412605E+04 4.72341711E+00 3.71757846E+00 1.27391260E-03 2.17347251E-06    3
-3.48858500E-09 1.65208866E-12 4.58723866E+04 1.75297945E+00 4.70504920E+04    4
! E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
! http://burcat.technion.ac.il/dir/
! mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
! Accessed April 2016.
[CH2]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 31,
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
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
2564-86-5
COOH  HydroxyFormil Rad Equil   HOCO   SIGMA=1 STATWT=2 IAIBIC=35. NU=3316,1797,
1261,1088,620,615  REF=TSIV   HF298=-183.971+/-0.533 kJ equil  REF=ATcT C 2011;
{HF298=-181.32+/-2.3 kJ  REF=ATcT A;  HF298=-213.+/-13 KJ  REF=Gurvich 91;
HF298=-185.7 kJ  REF=Burcat G3B3}  Max Lst Sq Error Cp @ 6000 K 0.39%.  The
polynomials were adjusted by B. Ruscic.
COOH   equilib    A 5/14C  1.O  2.H  1.   0.G   200.000  6000.000  B  45.01744 1
HOCO   equilib    A 5/14C  1.O  2.H  1.   0.G   200.000  6000.000 1000.00      1
5.39206152E+00 4.11221455E-03-1.48194900E-06 2.39875460E-10-1.43903104E-14    2
-2.36480760E+04-2.23529091E+00 2.92207919E+00 7.62453859E-03 3.29884437E-06    3
-1.07135205E-08 5.11587057E-12-2.33512327E+04 1.12925886E+01-2.21265065E+04    4
! E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
! http://burcat.technion.ac.il/dir/
! mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
! Accessed April 2016
H298 =-44.33 kcal/mol [FAB/JAN05]
S298 = 60.07 cal/mol/K [FAB/JAN05]
Cp [FAB/JAN05] (polyfit RAS/GLA08a).
O=[C]O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 32,
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
    shortDesc = u"""E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical""",
    longDesc = 
u"""
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
2143-58-0
CH3O2  METHYLPEROXIDE RAD (CH3OO)  SIGMA=3  STATWT=2  IA= 1.6120 IB= 7.4226
IC=8.4941    NU=3380,3167,3075,1511,1500,1462,1228,1169,1140,929,493,137.8
HF298=2.854+/-2. kcal  REF=Burcat G3B3  {HF298=12.055+/-0.9 kJ  REF=ATcT C 2011;
HF298-9.0+/-5.1 kJ   REF=Janoscheck IUPAC Sheets JPC 102,(1998) 1770.}  MAX LST
SQ ERROR Cp @ 6000 K 0.50 %
CH3O2 Methyl Per  T04/10C  1.H  3.O  2.   0.G   200.000  6000.000  B  47.03332 1
5.55530486E+00 9.12236137E-03-3.23851661E-06 5.18713798E-10-3.08834151E-14    2
-1.03569402E+03-3.99158547E+00 4.97169544E+00-5.29356557E-03 4.77334149E-05    3
-5.77065617E-08 2.22219969E-11-1.29022161E+02 2.81501182E+00 1.43618036E+03    4
H. Hashemi, et al., High-Pressure Oxidation of  Methane, 2016
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 33,
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
Marshall and Glarborg, Proc. Combust. Inst. 35 (2015) 153?160
64-18-6
CH2O2  METHANOIC (FORMIC) ACID  HCOOH MONOMER  STATWT=1 SIGMA=1  IA=1.0953
IB=6.9125  IC=8.0078  BROT=24.96 ROSYM=1 V(1)=2011.  V(2)=3123. V(3)=192.
NU=3570,2943,1770,1387,1229,1105,625,1033  NEL=100  HF298=-378.5+/-0.25 kJ for
equil. mix.  REF=ATcT C   {HF298=-378.6 kJ  REF=CHAO & ZWOLINSKI JPCRD 7.(1978),
363   HF298=-363.9 kJ  REF=Burcat G3B3 calc}  Max Lst Sq Error Cp @6000 K 0.47%
The polynomials were adjusted by B Ruscic..
O=CO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 34,
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
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
COO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 35,
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
Fabian WMF Janoschek R J Mol Struct THEOCHEM 2005, 713, 227?234
CL Rasmussen J Hansen P Marshall P Glarborg Int J Chem Kinet 40 (2008) 454-480
fitted to litt data
Low T polynomial Tmin changed from 298.15 to 298.0 K when importing to RMG.
[O]C=O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 36,
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
Accessed April 2016
CH2CHOH           BOZ06 C   2H   4O   1    0G   300.00   1500.00  1500.        1 !SIL/BOZ06
-0.09611513E+01 0.39459205E-01-0.45157340E-04 0.26804884E-07-0.62063875E-11    2 !H298 =-30.000kcal/mol
-0.16213850E+05 0.26561595E+02-0.09611513E+01 0.39459205E-01-0.45157340E-04    3 !S298 = 61.73 cal/mol/K
0.26804884E-07-0.62063875E-11-0.16213850E+05 0.26561595E+02                   4 !
! G da Silva CH Kim J Bozzelli W J Phys Chem A 110 (2006) 7925?7934
2154-50-9
C2H5O  ETHYL-OXIDE RAD (CH3CH2O) SIGMA=1 STATWT=2  IA=2.1281  IB=8.8117
IC=9.9060  Ir=0.4303  ROSYM=3  V(3)=737.5 cm-1     NU=3015,3004,2937,2824,
2790,1468,1458,1378,1360,1321,1206,1064,1046,872,856,475,406 T0=355 IA=2.3996
IB=8.1338  IC=9.4591  IR=0.4375  V(3)=1029.5 cm-1 ROSYM=3 NU=3040,3028,2951,
2866,2850,1514,1471,1445,1356,1268,1216,1107,934,912,874,577,369,(249 torrsion)
HF298=-11.47+/-0.5 kJ  REF=ATcT C 2011  {HF298=-13.6+/-4.0 kJ HF0=-0.2+/-4.0 kJ
REF=DeTuri & Ervin JPC 103 (1999),6911 for HF298 and G3MP2B3 calculations for
the vibrations and moments of inertia.  Ruscic et al JPCRD 34,(2005),573}  MAX
LST SQ ERROR @ 6000 K 0.56 %.
C2H5O  CH3CH2O*   T06/11C  2.H  5.O  1.   0.G   200.000  6000.000  B  45.06050 1.
CC[O]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 37,
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
    shortDesc = u"""H298 = 15.79 kcal/mol [JAN/ROS04]""",
    longDesc = 
u"""
H298 = 15.79 kcal/mol [JAN/ROS04]
S298 = 65.89 cal/mol/K [JAN/ROS04]
Cp(T) scaled Cp[CH3OO](T) to Cp298 = 14.89 cal/mol/K [JAN/ROS04].
[CH2]OO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 38,
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
Accessed April 2016
75-07-0
C2OH4  ACETALDEHYDE (CH3CHO)  STATWT=1.  SIGMA=1.  IA=2.76748  IB=6.9781
IC=9.03498  Ir=0.44 ROSYM=3 V(3)=412.03 cm-1  Nu=3005,2967,2917,2822,1743,1441,
1420,1400,1352,1113,919,867,763,509   HF298=-166.19 kJ  REF=CHAO, HALL,MARSH &
WILHOIT JCPRD 15, (1986) p.1369   {HF298=-165.364+/-0.3 kJ REF=ATcT C;
HF298=-170.7+/-1.5 kJ  REF=Wiberg & Croker JACS 113,(1991),3447}  Max Lst
Sq Error Cp @ 6000 K 0.59%..
CC=O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 39,
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
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
16520-04-0
C2H5O  CH2-O-CH3 RAD   SIGMA=1  STATWT=2.   IA=1.7787  IB=7.8857  IC=9.0727
IR(CH2)=0.30289  V(3)=700 cm-1   ROSYM=2   IR(CH3)=0.47197  V(3)=951 cm-1
ROSYM=3  NU=3262,3155,3112,3079,3020,1530,1521,1515,1479,1301,1264,1183,1151,
976,678,431  HF298=0.98 kJ  HF0=14.08 kJ  REF=Janoshcek Rossi IJCK 36,(2004),
661  {HF298=-2 kcal REF=Benson;  HF298=-2.8+/-1.2 kcal REF=MacMillen Golden 1982
HF298=-1.2 kcal REF=NIST 94} MAX LST SQ ERROR Cp @ 6000 K 0.52 %.
C2H5O  CH3-O-CH2  A10/04C  2.H  5.O  1.   0.G   200.000  6000.000  B  45.06050 1
5.94067593E+00 1.29906358E-02-4.56921036E-06 7.26888932E-10-4.30599587E-14    2
-2.58503562E+03-4.52841964E+00 4.53195381E+00 7.81884271E-03 1.94968539E-05    3
-2.74538336E-08 1.06521135E-11-1.70629244E+03 5.06122980E+00 1.15460803E+02    4

CH3OCH2    7/20/98 THERMC   2H   5O   1    0G   300.000  5000.000 1376.000    21
8.17137842E+00 1.10086181E-02-3.82352277E-06 5.99637202E-10-3.50317513E-14    2
-3.41941605E+03-1.78650856E+01 2.91327415E+00 2.03364659E-02-9.59712342E-06    3
2.07478525E-09-1.71343362E-13-1.18844240E+03 1.16066817E+01                   4
Dooley et al., Int. J. Chem. Kinet. 42 (2010) 527?549 :
Goldsmith et al., J. Phys. Chem. A 2012, 116, 3325?3346
64-17-5
C2H6O  ETHANOL (C2H5OH)    STATWT=1.  SIGMA=3 ROSYM(CH3)=3.  ROSYM(OH)=1
This is an equilibrium mixture of one trans and two gauche isomers, therefore
sigma was set artificialy to 3. The two gauche isomers are equal. The trans
values are: IAIBIC=218.459  Brot(CH3)=6.4144 cm-1 V(3)CH3=1166. Brot(OH)=21.07
cm-1  V(OH) (1)=57 (2)=8.025 (3)=395.  NU=3659,2985,2939,2900,1460,1430,1395,
1320,1245,1055,1026,883,422,2887(2),1460,1270,1117,801 The cis values are:
IAIBIC=233.455E-117 Brot(CH3)=6.416 cm-1 Brot(OH)=20.94 cm-1  V(3) CH3=1331
V(OH) as for trans NU=3675,2985,2939,2900,1460(2),1430,1395,1320,1245,1055,1026,
887,596,2887(2)1270,1070,801  HF298=-234.95 kJ  REF=CHAO, HALL, MARSH & WILHOIT
JPCRD 15 (1986),1369.  {HF298=-234.56+/-0.2 kJ  REF=ATcT A}.
CCO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 40,
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
Accessed April 2016
6912-06-7
C2H3O Vinyl oxy radical (CH2=CHO*)  This radical is in resonance with the next
radical CH2=CHO*<==>*CH2-CH=O which is probably preferred, therefore both
radicals have the same polynomial   SIGMA=1  STATWT=2  Ia=1.2527  Ib=7.3645
Ic=8.6172  Nu=3284,3174,2965,1572,1495,1416,1168,985,976,751,505,459
HF298=12.75 kJ  HF0=20.189 kJ  REF=Burcat G3B3 calc  Max Lst Sq Error Cp @
6000 K 0.49%
CH2=CHO*  Vinyl-  T04/06C  2.H  3.O  1.   0.G   200.000  6000.000  B  43.04462 1.
[CH2]C=O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 41,
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
Accessed April 2016
4422-54-2
C2H5O  (CH2CH2OH) RADICAL  SIGMA=1  STATWT=2  IA=2.1001  IB=8.6720  IC=10.0014
REF=Chem3D  Ir(CH2)=0.79  Ir(OH)=.1363  ROSYM(OH)=1  V(3)=70.3 cm-1 ROSYN(CH2)=2
V(3)=1049  cm-1  REF= Burcat, Miller & Gardiner  TAE Report 504 1983  NU=3705,
3093,2985,2855,2811,1501,1458,1409,1254,1223,1102,1042,951,853,433,376
REF=Yamada, Bozzelli, Lay JPC A 103 (1999),7646  Vib=scaled x 0.9;  HF298=-25.82
+/-0.6 kJ  REF=ATcT C 2011 {HF298=5.70+/-0.85 kcal  REF=Bozzelli JCP 105,(2001),
9543     MAX LST SQ ERROR Cp @ 6000 K 0.48 %.
[CH2]CO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 42,
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
Accessed April 2016
2348-46-1
C2H5O  (CH3CHOH)    RADICAL  SIGMA =1  STATWT=2  IA=1.8971    IB=8.9667
IC=10.2405    IR(CH3)=.47087 IR(OH)=.14477  ROSYM(CH3)=3  V(3)=1158. cm-1
ROSYM(OH)=1  V(3)=70.3 cm-1  NU=3734,3203,3164,3027,2956,1519,1500,1459,1425,
1327,1213,1072,1037,923,612,407.  REF=Janoschek & Rossi Int.J. Chem. Kinet 36
(2004),661  HF298=-55.29+/-0.6 kJ  REF=ATcT C 2011   {HF298=-54.03+/-4.0 kJ
REF=Janoschek & Rossi Int.J. Chem. Kinet 36 (2004),661;  HF298=-13.34+/-.85 kcal
REF=Bozzelli et al JCP 105,(2001),9543; HF298=-5.0 KCAL  REF= Benson.}  Max  Lst
Sq Error Cp @ 6000 K  0.48%.
C[CH]O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 43,
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
Accessed April 2016
557-75-5
C2H4O Vinyl Alcohol IA=1.363243 IB=7.9930197 IC=9.3562625  NU=412,470.5,693.4,
922,926.3,1029,1054.3,1293,1315.3,1434.7,1645.4,2964.6,3022,3050,3461
BROT=6.414  ROSYM=1  V(3)=1067. V(6)=-1.49 Ref= Ab-Initio Calc Karni,Oref &
Burcat TAE Report 643 1989. HF298=-29.8 kcal REF=Holm & Losing JACS 104
(1982) 2648. {HF298=-28.85+/-2. kcal  REF=Burcat G3B3}     Max Lst Sq Error Cp @
6000 K 0.45%
C2H4O vinyl alco  T03/10C  2.H  4.O  1.   0.G   200.000  6000.000  B  44.05256 1.
C=CO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 44,
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
Accessed April 2016
3170-69-2
C2H3O  Acetyl Radical  CH3*CO  SIGMA=1  STATWT=2  A=2.9436  B=0.334  C=0.3186
BROT=10.51589   ROSYM=3  V(3)=92 1/cm   NU=2904,2903,2826,1886,1405,1402,1325,
1025,925,817,454  REF=NIMLOS SODERQUIST & ELLISON JACS 111,(1989),7675
HF298=-10.3+/-1.8 kJ  REF=Niiaranen, Gutman & Krasnoperov J. Phys. Chem. 96
(1992) 5881.; Ruscic et al JPCRD 2003  {HF298=-9.68+/-0.4 kJ  REF=ATcT C 2011;
HF298=-13.3 kJ G4  REF=Marochkin & Dorofeeva Comp. Theor. Chem 991,(2012},182}
Max Lst Sq Error Cp @ 6000 K 0.62%.
C[C]=O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 45,
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
    shortDesc = u"""Lopez et al., Experimental and Kinetic Modeling Study of C2H2 Oxidation at High Pressure, Int. J. Chem. Kin., 2016""",
    longDesc = 
u"""
Lopez et al., Experimental and Kinetic Modeling Study of C2H2 Oxidation at High Pressure, Int. J. Chem. Kin., 2016.
[CH]=CO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 46,
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
Accessed April 2016
463-51-4
C2H2O  KETENE SIGMA=2  A0=.940922  B0=.343370  C0=.330736  NU=3070.4,2152.6,
1387.5,1116,3158.7,977.8,439,587.3,528.4  REF=Duncan et al J Mol Spec 125,(1987)
196  HF298=-48.57+/-0.14 kJ  REF=ATcT A   {HF298=-11.4+/-0.4 kcal REF= Vogt,
Williamson & Beauchamp JACS 100 (1978),3478;  HF298=-49.576 kJ  REF=B McBride
NASA TP-2002-211556}  MAX Lst Sq ERROR Cp @ 6000 K 0.43%.
CH2CO,ketene      g 4/02C  2.H  2.O  1.   0.G   200.000  6000.000  B  42.03668 1.
C=C=O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 47,
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
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
32038-79-2
C2H2O  ETHYNOL  HCC-OH   SIGMA=1  IA=0.121318  IB=8.4583675  IC=8.5796856
NU=346,383,523,600,1072,1232,2198,3339,3501  REF= M. JACOX JPCRD 19,(1990),1469
HF298=22.3+/-4.37 KCAL REF=Allendorf & C. Melius BAC/MP4 Sandia Database 2002.
{HF298=92.73+/-1.33 kJ  REF=ATcT C 2011; HF298=21.536+/-2 kcal  REF=Burcat G3B3}
Max Lst Sq Error Cp @ 6000 K 0.31%
HCC-OH  Ethynol   T12/09C  2.H  2.O  1.   0.G   200.000  6000.000  B  42.03668 1.
C#CO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 48,
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
in Marshall and Glarborg, Proc. Combust. Inst. 35 (2015) 153?160
CH3OCO     4/15/ 8 THERMC   2H   3O   2    0G   300.000  5000.000 1423.000    21
1.02053021E+01 7.12269668E-03-2.57594255E-06 4.15104358E-10-2.47168227E-14    2
-2.46870635E+04-2.86407262E+01 2.83313145E+00 1.53447505E-02 1.89583962E-06    3
-7.70200413E-09 2.41564410E-12-2.13431832E+04 1.39524183E+01                   4
Dooley et al., Int. J. Chem. Kinet. 42 (2010) 527?549
Goldsmith et al., J. Phys. Chem. A 2012, 116, 3325?3346
Goldsmith et al., J. Phys. Chem. A 2012, 116, 3325?3346
Dooley et al., Int. J. Chem. Kinet. 42 (2010) 527?549 ==> FROM NUIG
JIM/GLA08
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=COO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 49,
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
Accessed April 2016
51095-15-9
C2HO  KETENYL RAD  T0=0  SIGMA=1  STATWT=2  A0=34.99  B0=0.3640  C0=0.3602
NU=3371,2099,1249,568,506,511  T0=643  STATWT=2  B0=0.3597 Nu=3485,2108,1302,
398(2),522(2)  T0=33423.92  STATWT=4  SIGMA=1    B0=0.324  Nu=3485,2338,1037,
416(2),365(2)  REF=Osborn, Mordaunt, Choi, Bise, Neumark, JCP 106,(1997),10087
Molecular data for the ground states and nu1 and nu2 for the excited
states were taken from ab initio calculations of Szalay, Fogarasi, Nemes, Chem.
Phys. Lett. 263,(1996) 91. For the lowest non-linear level the CCSD(T)/PVTZ
basis set data were used except for the missing nu5 frequency which was taken
from EOMIP-CCSD/TZ2P (table 1). The Renner-Teller split yields a linear
conformation. For that conformation CCSD(T)/TZ2P data were used (table2). For the
excited state, nu1 and nu2 were taken from the EOM-IP/TZ2P data (table 3).
The remaining data for the excited state were taken from: Jacox JPCRD 27,(1998),
115-393. HF298=178.242+/-0.68  REF=ATcT A  {HF298=177.97+/-0.59 kJ  REF=ATcT C
2011;  HF298=178.3+/-1.5 kJ  REF=Szalay, Tajti & Stanton  Mol Phys. 103,(2005),
xxx;  HF298=42.4+/- 2.1 kcal  REF=Oakes, Jones, Blerbaum & Ellison  JPC 87,
(1983),4810}   MAX LST SQ ERROR Cp @ 1300 K 0.31%.
C#C[O]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 50,
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
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
C=CO[O]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 51,
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
    shortDesc = u"""J Gimenez CL Rasmussen MU Alzueta P Marshall P Glarborg Proc. Combust. Inst. 32 (2009) 367-375""",
    longDesc = 
u"""
J Gimenez CL Rasmussen MU Alzueta P Marshall P Glarborg Proc. Combust. Inst. 32 (2009) 367-375
107-31-3
C2H4O2  MethylFormate HCOOCH3  SIGMA=1  STATWT=1  IA=4.1878  IB=12.2832
IC=15.9404  Ir(CH3)=0.51926  ROSYM=3  V(3)=272. cm-1  Nu=3188,3154,3079,3071,
1832,1524,1514,1488,1417,12461197,1185,1037,947,760,348,300  HF298=-357.796+/-
0.6  kJ  REF=ATcT C 2011  {HF298=-363.6+/-8. kJ;  REF=Burcat G3B3 calc;
HF298=-362. kJ REF=Hine & Klueppel JACS 96,(1974),2924;  HF298=-355.5 kJ
REF=Hall & Baldt JACS 93,(1971),140; HF298=-358.15+/-7. kJ  REF=Allendorf-Melius,
Sandia Database BAC/MP4 2008; HF298=-371.5+/-15. kJ/mol  REF=Catoire et al IJCK
39,(2007),481}  Max Lst Sq Error Cp @ 6000 K 0.58 %.
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
CH3OCHO    4/15/ 8 THERMC   2H   4O   2    0G   300.000  5000.000 1452.000    21
1.06569185E+01 9.38112723E-03-3.37677206E-06 5.42455596E-10-3.22295330E-14    2
-4.84651758E+04-3.27579443E+01 2.43526200E+00 1.81948325E-02 1.91236590E-06    3
-8.44248945E-09 2.61758910E-12-4.46514732E+04 1.49522265E+01                   4
! Dooley et al., Int. J. Chem. Kinet. 42 (2010) 527?549
3170-61-4
C2H5OO PEROXYETHYL RADICAL  STATWT=2  IA=2.4505  IB=18.5705  IC=19.984
Ir(CH3)=0.49487  ROSYM=3  V(3)=1144. cm-1  Ir(C2H5)=2.859  ROSYM=1  V(3)=1479.
cm-1  NU=2955,2936,2934,2901,2874,1493,1467,1454,1410,1371,1259,1152,1145,1129,
1006,860,786,491,300,231,91.9  REF=Melius MP4 A40 1988   HF298=-6.86 kcal
REF=Atkinson et. al, JPCRD 28 (1999),191 {HF298= -2.32 kcal  REF=Melius 1988;
HF298=-4. kcal REF=NIST 1994 estimate} Max Lst sq Error Cp @ 200 & 6000 K 0.52%
C2H5O2  C2H5OO    T10/10C  2.H  5.O  2.   0.G   200.000  6000.000  B  61.05990 1.
CCO[O]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 52,
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
9  H u0 p0 c0 {3,S}
10 O u1 p2 c0 {4,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.556418,0.0412519,-3.10647e-05,1.2293e-08,-1.99691e-12,-22047.4,24.9983], Tmin=(298,'K'), Tmax=(800,'K')),
            NASAPolynomial(coeffs=[0.556418,0.0412519,-3.10647e-05,1.2293e-08,-1.99691e-12,-22047.4,24.9983], Tmin=(800,'K'), Tmax=(2500,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2500,'K'),
    ),
    shortDesc = u"""Bozzelli 2015; PM""",
    longDesc = 
u"""
Bozzelli 2015; PM
40135-01-1
H298 =-40.30 kcal/mol
S298 = 78.1 cal/mol/K
Low T polynomial Tmin changed from 300.0 to 298.0 K when importing to RMG.
[O]OCCO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 53,
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
CH2CH2OOH CBS-Q//B3LYP  C   2H   5O   2    0G   300.000  5000.000 1405.000    31
1.02538840E+01 1.07067444E-02-3.57347410E-06 5.45355058E-10-3.12522454E-14    2
1.31664933E+03-2.19734246E+01 3.26495428E+00 2.91609353E-02-2.29001252E-05    3
9.93668113E-09-1.78518267E-12 3.56067276E+03 1.48547673E+01                   4
! W-C Ing CY Sheng JW Bozzelli Fuel Proc Tech 83 (2003) 111?145
CH2CH2OOH               H  5 C  2 O  2      G     298.0    6000.0    1000.0    1
9.94719358E+00 1.14730988E-02-3.96870805E-06 6.23121554E-10-3.65488844E-14    2
2.33352925E+03-2.17163906E+01 6.49571381E-01 4.51634030E-02-6.02115954E-05    3
4.92121132E-08-1.67753354E-11 4.73413790E+03 2.49288138E+01                   4
! Goldsmith et al., J. Phys. Chem. A 2012, 116, 3325?3346
Dooley et al., Int. J. Chem. Kinet. 42 (2010) 527?549 ==> FROM NUIG
3031-74-1
C2H6O2 PEROXYETHANE C2H5O-OH  SIGMA=1 STATWT=1   IA=21.607  IB=20.135  IC=2.7509
IR(CH3)=0.49484  ROSYM=3  V(3)(CH3)=1143.7 cm-1 IR(C2H5)=2.859  ROSYM=1
V(3)(C2H5)=1479.46 cm-1 IR(OH)=0.1428  ROSYM=1 V(3)(OH)=2427.3 cm-1  REF=Lay et
al JPC 100,(1996),8240   NU=3697,3141,3135,3068,3063,3021,1561,1532,1514,1441,
1407,1387,1286,1201,1171,1051,960,836,478,296   HF298=-38.738+/-2 kcal
REF=Burcat G3B3  {HF298=-41.5 Kcal REF=Lay et al JPC 100,(1996),8240}   Max Lst
Sq Error Cp @ 200 & 6000 K 0.50%.
CCOO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 54,
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
CL Rasmussen JG Jacobsen P Glarborg Int J Chem Kinet 40 (2008) 778-807
126898-98-4
C2H4OOH HydroperoxylEthyl Radical  SIGMA=1  STATWT=2  IA=2.5032  IB=19.3176
IC=20.7280  Nu=3699,3295,3187,3029,2972,1529,1477,1398,1380,1239,1160,1070,1004,
932,823,510,461,299,206,130,121  HF298=11.841+/-2. kcal  REF=Burcat G3B3
{HF298=12.2+/-0.9 kcal  REF=Green et al JPC A 116,(2012),9033}  Max Lst Sq Error
Cp @ 6000 K 0.44%.
!C2H4OOH hydroper  T 4/15C  2.H  5.O  2.   0.G   200.000  6000.000  B  61.05990 1.
[CH2]COO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 55,
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
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 56,
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
Accessed April 2016
75-21-8
C2H4O  OXYRANE (ETHYLENE OXIDE)    SIGMA=2  IA=3.2793  IB=3.8059  IC=5.9511
NU=3006,1498,1271,1120,877,3063,1300,860,3006,1472,1151,892,3065,1142,822
REF=SHIMANOUCHI  HF298=-52.635 kJ FROM JANAF 1985. HF0=-40.082 kJ {HF298=-53.668
kJ REF=Burcat G3B3 calc 1/2005; HF298=-52.681+/-1.32  REF=ATcT C 2011)   Max Lst
Sq Error Cp @ 200 K ***1.17%*** @ 6000 K  0.59%..
C1CO1
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
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
Goldsmith et al., J. Phys. Chem. A 2012, 116, 3325?3346
31586-84-2
C2H3O  OXIRANE (ETHYLENE OXIDE) RADICAL SIGMA=1  STATWT=2  IA=2.8160  IB=3.5503
IC=5.6365    Nu=3204,3144,3114,1551,1366,1195,1133,1089,1049,949,817,793
HF298=164.473 kJ  HF0=172.90 kJ  REF=Burcat G3B3 calc {HF298 = 139.83 KJ est of
THERM}.  Max Lst Sq Error Cp @ 200 K *** 1.0% ***  @ 6000 K 0.51%.
C2H3O Oxyrane Rad A 1/05C  2.H  3.O  1.   0.G   200.000  6000.000  B  43.04462 1.
[CH]1CO1
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 58,
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
    shortDesc = u"""Marshall and Glarborg, Proc. Combust. Inst. 35 (2015) 153?160""",
    longDesc = 
u"""
Marshall and Glarborg, Proc. Combust. Inst. 35 (2015) 153?160
107-22-2
C2H2O2 (CHO-CHO) Trans-Cis-GLYOXAL   SIGMA=2   T0=0 (trans)  STATWT=1
IAIBIC=504.42  ROSYM=1  Brot1=4.213  Brot2=-1.117  Brot3=0.421  Brot4==0.126
Brot5=0.040  Brot6=-0.015 ROSYM=1  V(1)=1588.   V2=1140. V(3)=-59.0  V(4)-110.9
V(5)=40. V(6)=0  NEL=150   REF=Dorofeeva  JPCRD 30,(2001),475   NU=2843,1744,
1353,1066,551,801,1048,2835,1732,1312,339   REF= SCUSERIA  & SCHAEFER JACS 111,
(1989),7761
T0=1555. (Cis) SIGMA=2  STATWT=1  IAIBIC=710.17 (No internal rotation for the
cis exited state B. McBride and Zeleznik) Nu=2841,1746,1369,827,284.5,1050,750,
2810,1761,1360,825,10**10 (for the missing frequency or rotation!)
HF298=-212.082+/-0.8 kJ REF=Dorofeeva JPCRD 30,(2001),475 & ATcT A  HF0=-213.38
kJ  {HF298=-212.0+/-0.79 KJ  REF=Fletcher & Pilcher  Trans Faraday Soc 66(1970),
794}  HF298=-193.249+/-0.8 for Cis only  REF=ATcT A   Max Lst Sq Error Cp @
1300 K 0.48%
!O(CH)2O Glyoxal   g 3/02C  2.H  2.O  2.   0.G   200.000  6000.000  B  58.03608 1.
O=CC=O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 59,
    label = "HOCH2O",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 O u1 p2 c0 {1,S}
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
Aranda, V., et al., Int. J. Chemical Kinet. 45.5 (2013): 283-294.
PM
Marshall and Glarborg, Proc. Combust. Inst. 35 (2015) 153?160.
[O]CO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 60,
    label = "OCHCO",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {5,D}
3 H u0 p0 c0 {1,S}
4 O u1 p2 c0 {1,S}
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
115-10-6
C2H6O  DIMETHYL-ETHER SIGMA=2  IAIBIC=170.493 IR=0.4291   ROSYM=3
V(3)=903.4 cal  Nu=2999(2),2935,2920,2820(2),1485,1467,1463,1459,1449,1432,1250,
1179,1178,1148,1104,931,424  HF298=-183.935+/-0.46 kJ REF=ATcT A  {HF298=-184.02
+/-0.43  kJ  REF=ATcT C 2011;  HF298=-184.05 kJ  REF=CHAO, HALL, MARSH & WILHOIT
JPCRD 15, (1986),1369}   Max Lst Sq Error Cp @ 1300 K 0.64%
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
CH3OCH3    7/20/98 THERMC   2H   6O   1    0G   300.000  5000.000 1368.000    21
8.27745836E+00 1.32135539E-02-4.53264362E-06 7.05316507E-10-4.09933283E-14    2
-2.61982700E+04-2.15190894E+01 1.50763450E+00 2.39914228E-02-8.68910500E-06    3
-9.66835762E-11 4.89319361E-13-2.32810894E+04 1.67317297E+01                   4
Dooley et al., Int. J. Chem. Kinet. 42 (2010) 527?549
Low T polynomial Tmin changed from 350.0 to 298.0 K when importing to RMG.
[O]C=C=O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 61,
    label = "CH3C(O)O",
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
            NASAPolynomial(coeffs=[1.30016,0.0237699,-1.23328e-05,5.87355e-11,1.44737e-12,-24433.8,20.8502], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.15457,0.0168675,-8.75191e-06,2.18661e-09,-2.13394e-13,-25230.4,5.95053], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1508""",
    longDesc = 
u"""
1508
J Gimenez CL Rasmussen MU Alzueta P Marshall P Glarborg Proc. Combust. Inst. 32 (2009) 367-375
W-C Ing CY Sheng JW Bozzelli Fuel Proc Tech 83 (2003) 111?145
Bozzelli 2015; PM.
CC([O])=O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 62,
    label = "C2",
    molecule = 
"""
multiplicity 3
1 C u1 p0 c0 {2,T}
2 C u1 p0 c0 {1,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.69386,-0.00184767,5.23713e-06,-3.83965e-09,8.61136e-13,98382.2,2.23677], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[3.25289,0.0012319,-4.50354e-07,7.49357e-11,-4.57925e-15,98373.7,3.9586], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""singlet        T05/09""",
    longDesc = 
u"""
singlet        T05/09
H. Hashemi, et al., High-Pressure Oxidation of Methane, 2016
Dooley, et al., Int. J. Chem. Kinet. 42 (2010) 527?549.
<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C2 species:
<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
12070-15-4
C2  singlet T0=0  SIGMA=2  STATWT=1 Be=1.7930  We=1854.71  WeXe=13.34
WeYe=-0.17  ALPHAE=0.0421  De=6.92E-6  REF=Hubert & Herzberg Webbook 2009 and
Gurvich 91  HF298=826.8+/-8. kJ  REF=Karton & Martin Mol Phys 107,(2009),977
{HF298=828.374+/-0.3 kJ  REF=ATcT C 2011;  HF298=830.457+/-10 kJ
REF=Gurvich 91; HF298=815.9 kJ REF=Van-Orden  Saykally Chem. Rev 98,(1998),2313}
Max Lst Sq Error Cp @ 1300 K 0.32%..
[C]#[C]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 63,
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
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
Wijaya, Sumathi, Green, J. Phys. Chem. A  107 (2003) 4908-4920
Dooley et al., Int. J. Chem. Kinet. 42 (2010) 527?549
Wijaya, Sumathi, Green, J. Phys. Chem. A  107 (2003) 4908-4920
Dooley et al., Int. J. Chem. Kinet. 42 (2010) 527?549
INGBOZ03 HC*OCOO.
INGBOZ03 C.*OCOOH
Dooley et al., Int. J. Chem. Kinet. 42 (2010) 527?549 ==> FROM NUIG
Dooley et al., Int. J. Chem. Kinet. 42 (2010) 527?549 ==> FROM NUIG
Wijaya, Sumathi, Green, J. Phys. Chem. A  107 (2003) 4908-4920
Dooley et al., Int. J. Chem. Kinet. 42 (2010) 527?549.
CC(=O)O[O]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 64,
    label = "CH3C(O)OOH",
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
            NASAPolynomial(coeffs=[1.43748,0.0339384,-2.25467e-05,5.77205e-09,1.28452e-13,-44070.1,19.8194], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.64132,0.0217231,-1.11238e-05,2.75682e-09,-2.67754e-13,-45140.9,-1.61172], Tmin=(1000,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (3000,'K'),
    ),
    shortDesc = u"""1508""",
    longDesc = 
u"""
1508
Bozzelli 2015; PM.
CC(=O)OO
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 65,
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
Accessed April 2016
12071-23-7
C2O        SIGMA=1  STATWT=3  B0=0.385  NU=1971,379.53(2),1063
T0=5310.   SIGMA=1  STATWT=2  B0=0.385  NU=1950,379.53(2),1063
T0=8190.   SIGMA=1  STATWT=1  B0=0.385  NU=2010,379.53(2),1063
T0=11651.  SIGMA=1  STATWT=6  B0=0.407  Nu=2046,594.75(2),1284  REF=Jacox 98
HF298=378.86+/-1.2 kJ  REF=ATcT C 2011  {HF298=291.04+/-12 kJ HF0=287.0 kJ
REF=Gurvich 91; HF298=286.6 kJ REF=JANAF; both values were declared erroneous by
Newmark JCP 108,(1998),4070 and Williams & Fleming Proc. Comb. Inst 31,(2007),
1109}  Max Lst Sq Error Cp @ 1200 K 0.23%..
[C]#C[O]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 66,
    label = "HCOH",
    molecule = 
"""
multiplicity 3
1 O u0 p2 c0 {2,S} {4,S}
2 C u2 p0 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {1,S}
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
CH**-OH cis T 9/09
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
from Xiong et al., Combust. Flame 161 (2014) 885?897
3-HCOH                  C   1H   2O   1    0G   298.150  6000.000 1000.00      1
0.62184153E+01 0.18701090E-02-0.17117529E-06-0.44242689E-10 0.63999362E-14    2
0.23486822E+05-0.97996712E+01 0.36513441E+01 0.46834047E-02 0.11223500E-05    3
-0.80289814E-09-0.77469460E-12 0.24561953E+05 0.49211311E+01                   4
1-HCOH                  C   1H   2O   1    0G   298.150  6000.000 1000.00      1
0.57648468E+01 0.21295433E-02-0.19747641E-06-0.51074096E-10 0.74327861E-14    2
0.10621689E+05-0.77241540E+01 0.32130404E+01 0.32250737E-02 0.24648079E-05    3
0.15449875E-08-0.27946393E-11 0.11899702E+05 0.76449280E+01                   4
19710-56-6
CHOH trans Hydroxymethylene trans  SIGMA=1  STATWT=1  IA=0.2970  IB=2.3197
IC=2.6167  Nu=3634,2827,1536,1330,1224,1099  REF=Burcat G3B3  HF298=108.16+/-.43
kJ  REF=Ruscic ATcT D 2013  {HF298=106.47+/-8 kJ  REF=Burcat G3B3}  Max Lst Sq
Error Cp @ 6000 K 0.45%.
CHOH trans  Hydr  T12/14C  1.H  2.O  1.   0.G   200.000  6000.000  B  30.02598 1
3.63500546E+00 5.45479239E-03-1.91285077E-06 3.03882849E-10-1.79904440E-14    2
5.28968001E+04 4.34543184E+00 4.72862471E+00-1.02082828E-02 4.21639528E-05    3
-4.64774522E-08 1.72559970E-11 5.31829869E+04 1.69093263E+00 5.44279146E+04    4

19710-56-6
CH2O cis  Hydroxymethylene  CH**OH cis  SIGMA=1 STATWT=1  IA=0.3036  IB=2.3235
IC=2.6260  Nu=3410,2727,1492,1346,1284.5,1019  HF298=30.0139+/-2 kcal
REF=Burcat G3B3  {HF298(trans)=108.19+/-0.42 kJ  HF298(cis)=126.66+/-0.42  kJ
REF=ATcT C 2011}  Max Lst Sq Error Cp @ 6000 K 0.49%..
[CH]O
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 67,
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
[CH]=C[O]
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

entry(
    index = 68,
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
            NASAPolynomial(coeffs=[3.28155,0.00697643,-2.38528e-06,-1.21078e-09,9.82042e-13,48319.2,5.92036], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.27807,0.00475623,-1.63007e-06,2.54623e-10,-1.4886e-14,48014,0.639979], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Vinylidene   T 7/11""",
    longDesc = 
u"""
Vinylidene   T 7/11
E Goos A Burcat B Ruscic Ideal gas thermochemical database with updates from active thermochemical
http://burcat.technion.ac.il/dir/
mirrored at http://garfield.chem.elte.hu/burcat/burcat.html.
Accessed April 2016
H2CC                    H  2 C  2 O  0      G     298.0    6000.0    1000.0    1
4.64241230E+00 4.45101029E-03-1.51900528E-06 2.36126946E-10-1.37457668E-14    2
4.83798914E+04-5.81740342E-01 2.62922318E+00 1.33909598E-02-1.98275249E-05    3
1.80813401E-08-6.47719980E-12 4.88573332E+04 9.20664687E+00                   4
! Goldsmith et al., J. Phys. Chem. A 2012, 116, 3325?3346
2143-69-3
CH2C   VINYLIDENE RADICAL  SIGMA=2  STATWT=1  IA=0.29432  IB=2.17453  IC=2.46885
NU=3344,3239,1710,1288,787,444 REF=OSAMURA, SCHAFER, GRAY & MILER J.A.C.S. 103
(1981) 1904.  HF298=412.272+/-0.61 kJ  REF=ATcT C 2011  {HF0=414.489 kJ
REF=Chen, Jonas,Kinsey & Field J Chem Phys 91,(1989),3976;  HF298=413.36+/-1.8
kJ REF=ATcT A}   Max Lst Sq Error Cp @ 6000 K 0.35%..
[C]=C
Imported from /home/alongd/Code/RMG-Py/importer/Hashemi/therm.DAT.
""",
)

