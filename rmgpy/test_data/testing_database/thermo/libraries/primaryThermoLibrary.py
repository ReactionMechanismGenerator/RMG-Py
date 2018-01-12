#!/usr/bin/env python
# encoding: utf-8

name = "primaryThermoLibrary"
shortDesc = u""
longDesc = u"""

"""
entry(
    index = 1,
    label = "H2",
    molecule = 
"""
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.895,6.975,6.994,7.009,7.081,7.219,7.72],'cal/(mol*K)'),
        H298 = (0,'kcal/mol'),
        S298 = (31.233,'cal/(mol*K)','+|-',0.0007),
    ),
    shortDesc = u"""library value for H2""",
    longDesc = 
u"""

""",
)

entry(
    index = 2,
    label = "H",
    molecule = 
"""
multiplicity 2
1 H u1 p0 c0
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.968,4.968,4.968,4.968,4.968,4.968,4.968],'cal/(mol*K)'),
        H298 = (52.103,'kcal/mol','+|-',0.001),
        S298 = (27.419,'cal/(mol*K)','+|-',0.0005),
    ),
    shortDesc = u"""library value for H radical""",
    longDesc = 
u"""

""",
)

entry(
    index = 3,
    label = "O2",
    molecule = 
"""
multiplicity 3
1 O u1 p2 c0 {2,S}
2 O u1 p2 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.0233,7.1986,7.4285,7.6673,8.0656,8.3363,8.7407],'cal/(mol*K)'),
        H298 = (-0.0010244,'kcal/mol'),
        S298 = (49.0236,'cal/(mol*K)'),
    ),
    shortDesc = u"""from GRI-Mech 3.0""",
    longDesc = 
u"""

""",
)

entry(
    index = 5,
    label = "CO3s1",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,D}
2 O u0 p2 c0 {1,S} {3,S}
3 O u0 p2 c0 {1,S} {2,S}
4 O u0 p2 c0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.82,13.53,14.8,15.76,17.05,17.85,18.84],'cal/(mol*K)'),
        H298 = (-37.9,'kcal/mol'),
        S298 = (61.28,'cal/(mol*K)'),
    ),
    shortDesc = u"""Mebel et al (2004) http:""",
    longDesc = 
u"""

""",
)

entry(
    index = 6,
    label = "CO3t1",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 O u0 p2 c0 {1,D}
3 O u1 p2 c0 {1,S}
4 O u1 p2 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([12.43,14.02,15.23,16.16,17.4,18.14,19.01],'cal/(mol*K)'),
        H298 = (-11.5,'kcal/mol'),
        S298 = (64.53,'cal/(mol*K)'),
    ),
    shortDesc = u"""Mebel et al (2004) http:""",
    longDesc = 
u"""

""",
)

entry(
    index = 7,
    label = "CO3t2",
    molecule = 
"""
multiplicity 3
1 C u1 p0 c0 {2,S} {3,D}
2 O u0 p2 c0 {1,S} {4,S}
3 O u0 p2 c0 {1,D}
4 O u1 p2 c0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.48,14.83,15.83,16.59,17.62,18.26,19.05],'cal/(mol*K)'),
        H298 = (28.7,'kcal/mol'),
        S298 = (67.06,'cal/(mol*K)'),
    ),
    shortDesc = u"""Mebel et al (2004) http:""",
    longDesc = 
u"""

""",
)

entry(
    index = 8,
    label = "cyclopropene12diyl",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u1 p0 c0 {1,S} {3,D}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([11.32,13.42,15.16,16.58,18.74,20.29,22.63],'cal/(mol*K)'),
        H298 = (188.13,'kcal/mol'),
        S298 = (59.26,'cal/(mol*K)'),
    ),
    shortDesc = u"""doi:10.1016/j.chemphys.2008.01.057 and doi:10.1021/jp003224c""",
    longDesc = 
u"""

""",
)

entry(
    index = 9,
    label = "cyclopropynylidyne",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.43,11.93,13.18,14.18,15.61,16.59,17.99],'cal/(mol*K)'),
        H298 = (169.8,'kcal/mol'),
        S298 = (57.47,'cal/(mol*K)'),
    ),
    shortDesc = u"""doi:10.1016/j.chemphys.2008.01.057 and doi:10.1021/jp003224c""",
    longDesc = 
u"""

""",
)

entry(
    index = 10,
    label = "OCCO(S)",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,D}
2 C u0 p0 c0 {1,D} {4,D}
3 O u0 p2 c0 {1,D}
4 O u0 p2 c0 {2,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([13.66,14.78,15.66,16.41,17.57,18.39,19.53],'cal/(mol*K)'),
        H298 = (18.31,'kcal/mol'),
        S298 = (90.63,'cal/(mol*K)'),
    ),
    shortDesc = u"""QCI""",
    longDesc = 
u"""

""",
)

entry(
    index = 11,
    label = "OCCO",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,T} {3,S}
2 C u0 p0 c0 {1,T} {4,S}
3 O u1 p2 c0 {1,S}
4 O u1 p2 c0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([14.21,15.2,15.98,16.65,17.72,18.49,19.58],'cal/(mol*K)'),
        H298 = (5.6,'kcal/mol'),
        S298 = (61.18,'cal/(mol*K)'),
    ),
    shortDesc = u"""QCI""",
    longDesc = 
u"""

""",
)

entry(
    index = 12,
    label = "C3H2",
    molecule = 
"""
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 {1,D} {3,S} {5,S}
3 C u0 p1 c0 {1,S} {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.61,12.79,14.81,16.48,18.71,20.33,22.59],'cal/(mol*K)'),
        H298 = (113.99,'kcal/mol'),
        S298 = (56.44,'cal/(mol*K)'),
    ),
    shortDesc = u"""Burcat's recommended value for 16165-40-5: C3H2 CYCLOPROPENYLIDENE BI-RADICAL SINGLET on 3/25/2011 (MRH converted from NASA-7 to RMG format)""",
    longDesc = 
u"""

""",
)

entry(
    index = 13,
    label = "S2",
    molecule = 
"""
1 S u0 p2 c0 {2,D}
2 S u0 p2 c0 {1,D}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.81,8.17,8.4,8.54,8.7,8.78,8.87],'cal/(mol*K)'),
        H298 = (47.12,'kcal/mol'),
        S298 = (53.78,'cal/(mol*K)'),
    ),
    shortDesc = u"""CBS-QB3 value A.G. Vandeputte""",
    longDesc = 
u"""

""",
)

entry(
    index = 14,
    label = "S2JJ",
    molecule = 
"""
multiplicity 3
1 S u1 p2 c0 {2,S}
2 S u1 p2 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.79,8.14,8.35,8.51,8.75,8.94,9.31],'cal/(mol*K)'),
        H298 = (30.74,'kcal/mol'),
        S298 = (54.54,'cal/(mol*K)'),
    ),
    shortDesc = u"""from Chase thermo database""",
    longDesc = 
u"""

""",
)

entry(
    index = 15,
    label = "HCS",
    molecule = 
"""
multiplicity 2
1 C u1 p0 c0 {2,D} {3,S}
2 S u0 p2 c0 {1,D}
3 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.77,9.46,10.04,10.51,11.22,11.76,12.62],'cal/(mol*K)'),
        H298 = (71.7,'kcal/mol'),
        S298 = (56.37,'cal/(mol*K)'),
    ),
    shortDesc = u"""NIST value + B3LYP/cbsb7 entropy and heat cap A.G. Vandeputte""",
    longDesc = 
u"""

""",
)

entry(
    index = 16,
    label = "Ar",
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
    shortDesc = u"""Burcat Thermo Data""",
    longDesc = 
u"""
Ar HF298=0.  REF=C.E. Moore 'Atomic Energy Levels' NSRDS-NBS 35 (1971) p.211
""",
)

entry(
    index = 17,
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
    shortDesc = u"""Burcat Thermo Data""",
    longDesc = 
u"""
N2  HF298= 0.0 KJ  REF=TSIV  Max Lst Sq Error Cp @ 6000 K 0.29%
""",
)

entry(
    index = 18,
    label = "He",
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
    shortDesc = u"""Burcat Thermo Data""",
    longDesc = 
u"""
McBride, Heimel, Ehlers & Gordon "Thermodynamic Properties to 6000 K", 1963.
""",
)

entry(
    index = 19,
    label = "C(S)",
    molecule = 
"""
1 C u0 p2 c0
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.968,4.968,4.968,4.968,4.968,4.968,4.968],'cal/(mol*K)'),
        H298 = (200.397,'kcal/mol'),
        S298 = (33.393,'cal/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
H298: ATcT version 1.110
level of theory energy: CCSD(T)F12A/cc-pVQZ-F12//CCSD(T)/cc-pVQZ
level of theory frequency: B3LYP/6-311++g(d,p)//B3LYP/6-311++g(d,p)
""",
)

entry(
    index = 20,
    label = "C(T)",
    molecule = 
"""
multiplicity 3
1 C u2 p1 c0
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.968,4.968,4.968,4.968,4.968,4.968,4.968],'cal/(mol*K)'),
        H298 = (171.336,'kcal/mol'),
        S298 = (35.576,'cal/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
H298: ATcT version 1.110
level of theory energy: CCSD(T)F12A/cc-pVQZ-F12//CCSD(T)/cc-pVQZ
level of theory frequency: B3LYP/6-311++g(d,p)//B3LYP/6-311++g(d,p)
""",
)

entry(
    index = 21,
    label = "CH2(S)",
    molecule = 
"""
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.089,8.316,8.622,8.99,9.787,10.502,11.832],'cal/(mol*K)'),
        H298 = (102.541,'kcal/mol'),
        S298 = (45.197,'cal/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
H298: ATcT version 1.110
level of theory energy: CCSD(T)F12A/cc-pVQZ-F12//CCSD(T)/cc-pVQZ
level of theory frequency: B3LYP/6-311++g(d,p)//B3LYP/6-311++g(d,p)
""",
)

entry(
    index = 22,
    label = "CH2(T)",
    molecule = 
"""
multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.342,8.61,8.914,9.238,9.893,10.5,11.68],'cal/(mol*K)'),
        H298 = (93.559,'kcal/mol'),
        S298 = (46.636,'cal/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
H298: ATcT version 1.110
level of theory energy: CCSD(T)F12A/cc-pVQZ-F12//CCSD(T)/cc-pVQZ
level of theory frequency: B3LYP/6-311++g(d,p)//B3LYP/6-311++g(d,p)
""",
)

entry(
    index = 23,
    label = "CH4",
    molecule = 
"""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.621,9.648,10.941,12.344,14.887,16.967,20.535],'cal/(mol*K)'),
        H298 = (-17.814,'kcal/mol'),
        S298 = (44.473,'cal/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
H298: ATcT version 1.110
level of theory energy: CCSD(T)F12A/cc-pVQZ-F12//CCSD(T)/cc-pVQZ
level of theory frequency: B3LYP/6-311++g(d,p)//B3LYP/6-311++g(d,p)
""",
)

entry(
    index = 24,
    label = "NH(T)",
    molecule = 
"""
multiplicity 3
1 N u2 p1 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.95,6.945,6.962,7.006,7.173,7.387,7.891],'cal/(mol*K)'),
        H298 = (85.753,'kcal/mol'),
        S298 = (43.265,'cal/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
H298: ATcT version 1.110
level of theory energy: CCSD(T)F12A/cc-pVQZ-F12//CCSD(T)/cc-pVQZ
level of theory frequency: B3LYP/6-311++g(d,p)//B3LYP/6-311++g(d,p)
""",
)

entry(
    index = 25,
    label = "NH2(D)",
    molecule = 
"""
multiplicity 2
1 N u1 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.055,8.227,8.467,8.762,9.435,10.083,11.402],'cal/(mol*K)'),
        H298 = (44.467,'kcal/mol'),
        S298 = (46.516,'cal/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
H298: ATcT version 1.110
level of theory energy: CCSD(T)F12A/cc-pVQZ-F12//CCSD(T)/cc-pVQZ
level of theory frequency: B3LYP/6-311++g(d,p)//B3LYP/6-311++g(d,p)
""",
)

entry(
    index = 26,
    label = "NH3",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.486,9.051,9.737,10.484,11.921,13.163,15.503],'cal/(mol*K)'),
        H298 = (-10.889,'kcal/mol'),
        S298 = (45.986,'cal/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
H298: ATcT version 1.110
level of theory energy: CCSD(T)F12A/cc-pVQZ-F12//CCSD(T)/cc-pVQZ
level of theory frequency: B3LYP/6-311++g(d,p)//B3LYP/6-311++g(d,p)
""",
)

entry(
    index = 27,
    label = "O(S)",
    molecule = 
"""
1 O u0 p3 c0
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.968,4.968,4.968,4.968,4.968,4.968,4.968],'cal/(mol*K)'),
        H298 = (104.81,'kcal/mol'),
        S298 = (34.25,'cal/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
H298: ATcT version 1.110
level of theory energy: CCSD(T)F12A/cc-pVQZ-F12//CCSD(T)/cc-pVQZ
level of theory frequency: B3LYP/6-311++g(d,p)//B3LYP/6-311++g(d,p)
""",
)

entry(
    index = 28,
    label = "O(T)",
    molecule = 
"""
multiplicity 3
1 O u2 p2 c0
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.968,4.968,4.968,4.968,4.968,4.968,4.968],'cal/(mol*K)'),
        H298 = (59.567,'kcal/mol'),
        S298 = (36.433,'cal/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
H298: ATcT version 1.110
level of theory energy: CCSD(T)F12A/cc-pVQZ-F12//CCSD(T)/cc-pVQZ
level of theory frequency: B3LYP/6-311++g(d,p)//B3LYP/6-311++g(d,p)
""",
)

entry(
    index = 29,
    label = "OH(D)",
    molecule = 
"""
multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.954,6.946,6.951,6.973,7.08,7.251,7.719],'cal/(mol*K)'),
        H298 = (8.863,'kcal/mol'),
        S298 = (43.958,'cal/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
H298: ATcT version 1.110
level of theory energy: CCSD(T)F12A/cc-pVQZ-F12//CCSD(T)/cc-pVQZ
level of theory frequency: B3LYP/6-311++g(d,p)//B3LYP/6-311++g(d,p)
""",
)

entry(
    index = 30,
    label = "H2O",
    molecule = 
"""
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.038,8.18,8.379,8.624,9.195,9.766,11.019],'cal/(mol*K)'),
        H298 = (-57.797,'kcal/mol'),
        S298 = (45.084,'cal/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""
H298: ATcT version 1.110
level of theory energy: CCSD(T)F12A/cc-pVQZ-F12//CCSD(T)/cc-pVQZ
level of theory frequency: B3LYP/6-311++g(d,p)//B3LYP/6-311++g(d,p)
""",
)

entry(
    index = 31,
    label = "Cl2",
    molecule = 
"""
1 Cl u0 p3 c0 {2,S}
2 Cl u0 p3 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.73638,0.00783526,-1.45105e-05,1.25731e-08,-4.13247e-12,-1058.8,9.44557], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.74728,-0.000488582,2.68445e-07,-2.43476e-11,-1.03683e-15,-1511.02,-0.344539], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Burcat Thermo Data""",
    longDesc = 
u"""
Burcat Thermo Data

REFERENCE ELEMENT REF=Gurvich 1989 V1 py.1 p.177 HF298=0.00 kcal Max Lst 
Sq Error Cp @ 6000 **1.26%** (Cp @ 700 K 0.08%)
""",
)

entry(
    index = 32,
    label = "Cl",
    molecule = 
"""
multiplicity 2
1 Cl u1 p3 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.26062,0.00154154,-6.80284e-07,-1.59973e-09,1.15417e-12,13855.3,6.57021], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.94658,-0.000385985,1.36139e-07,-2.17033e-11,1.28751e-15,13697,3.1133], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Burcat Thermo Data""",
    longDesc = 
u"""
Burcat Thermo Data

HF298=121.302+/-0.008 kJ HF0=119.633+/- 0.008 kJ  REF=JANAF  {HF298=121.302
+/-0.002 kJ  REF=ATcT A}
""",
)

entry(
    index = 33,
    label = "HCl",
    molecule = 
"""
1 Cl u0 p3 c0 {2,S}
2 H  u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.46376,0.000476484,-2.00301e-06,3.31714e-09,-1.44958e-12,-12144.4,2.66428], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.75758,0.00145387,-4.79647e-07,7.77909e-11,-4.79574e-15,-11913.8,6.52197], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Burcat Thermo Data""",
    longDesc = 
u"""
Burcat Thermo Data

HYDROCHLORIC ACID CALCULATED FROM ORIGINAL TABLES  REF=Gurvich 1989
HF298=-92.31 kJ {HF298=-92.17+/-0.006 kJ   REF=ATcT C}  Max Lst Sq Error Cp @ 
6000 K 0.17%
""",
)

entry(
    index = 34,
    label = "Ne",
    molecule = 
"""
1 Ne u0 p4 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    shortDesc = u"""Burcat Thermo Data""",
    longDesc = 
u"""
McBride, Heimel, Ehlers & Gordon, "Thermodynamic Properties to 6000 K", 1963.
""",
)

