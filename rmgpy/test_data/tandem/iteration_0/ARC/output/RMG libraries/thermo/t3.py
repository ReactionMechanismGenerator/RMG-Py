#!/usr/bin/env python
# encoding: utf-8

name = "rmg_arc"
shortDesc = ""
longDesc = """
ARC v1.1.0
ARC project rmg_arc

Levels of theory used:

Conformers:       apfd/def2svp
TS guesses:       apfd/def2svp
Optimization:     b3lyp/6-31g(d,p) (NOT using a fine grid)
Frequencies:      b3lyp/6-31g(d,p)
Single point:     b3lyp/6-31g(d,p)
Rotor scans:      
Using bond additivity corrections for thermo

Using the following ESS settings: {'gaussian': ['rmg', 'pharos', 'c3ddb'], 'molpro': ['rmg', 'pharos'], 'qchem': ['pharos'], 'orca': ['c3ddb'], 'onedmin': ['pharos'], 'gromacs': ['rmg']}

Considered the following species and TSs:
Species C3H8_0 (run time: 0:01:16)
Species C3H7_1 (run time: 0:04:01)
Species C3H6_0 (run time: 0:01:46)
Species C3H6_1 (run time: 0:03:14)
Species C3H5_1 (run time: 0:02:11)
Species C3H5_2 (run time: 0:03:52)
Species C3H6_2 (run time: 0:04:27)
Species C3H4_0 (run time: 0:00:43)
Species C3H3_0 (run time: 0:00:38)
Species C3H4_1 (Failed!) (run time: None)
Species C3H4_2 (run time: 0:01:00)

Overall time since project initiation: 00:28:38
"""
entry(
    index = 0,
    label = "C3H8_0",
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
            NASAPolynomial(coeffs=[3.97109,0.00137612,7.92112e-05,-1.1372e-07,5.27266e-11,-15373.6,6.78537], Tmin=(10,'K'), Tmax=(558.573,'K')),
            NASAPolynomial(coeffs=[-1.21311,0.0385009,-2.04846e-05,5.26979e-09,-5.29621e-13,-14794.5,28.777], Tmin=(558.573,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (-127.841,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (257.749,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Bond corrections: {'C-C': 2, 'C-H': 8}

External symmetry: 2, optical isomers: 1

Geometry:
{'symbols': ('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'), 'isotopes': (12, 12, 12, 1, 1, 1, 1, 1, 1, 1, 1), 'coords': ((-0.023532, 0.574557, -0.222108), (-1.249429, -0.336352, -0.100521), (1.266957, -0.091617, 0.265961), (-0.196257, 1.496157, 0.348442), (0.098804, 0.88327, -1.268261), (-1.41858, -0.63466, 0.940516), (-1.121233, -1.252295, -0.688714), (-2.157653, 0.162057, -0.455224), (1.48548, -0.998776, -0.309076), (2.12731, 0.578797, 0.168831), (1.188133, -0.38114, 1.320154))}
""",
)

entry(
    index = 1,
    label = "C3H7_1",
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
            NASAPolynomial(coeffs=[3.7715,0.0248421,-6.52002e-05,1.93999e-07,-1.75965e-10,7029.17,7.99842], Tmin=(10,'K'), Tmax=(437.752,'K')),
            NASAPolynomial(coeffs=[0.29535,0.0339487,-1.87679e-05,5.05002e-09,-5.31144e-13,7550.6,24.3767], Tmin=(437.752,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (58.4395,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (232.805,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Bond corrections: {'C-C': 2, 'C-H': 7}

External symmetry: 1, optical isomers: 1

Geometry:
{'symbols': ('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H'), 'isotopes': (12, 12, 12, 1, 1, 1, 1, 1, 1, 1), 'coords': ((0.995037, -0.622692, 0.893694), (0.494701, -0.462675, -1.652624), (0.037854, -0.682138, -0.249305), (1.304502, 0.412374, 1.125671), (1.919629, -1.173398, 0.674056), (0.564152, -1.036272, 1.811261), (1.420557, -1.013785, -1.865828), (-0.259938, -0.772711, -2.382717), (0.718321, 0.599847, -1.857531), (-1.027566, -0.636849, -0.03711))}
""",
)

entry(
    index = 2,
    label = "C3H6_0",
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
            NASAPolynomial(coeffs=[4.01913,-0.00168835,7.97088e-05,-1.22632e-07,6.20359e-11,2159.96,6.97913], Tmin=(10,'K'), Tmax=(509.608,'K')),
            NASAPolynomial(coeffs=[-0.200681,0.0314332,-1.77811e-05,4.90209e-09,-5.27987e-13,2590.06,24.4927], Tmin=(509.608,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (17.9556,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (207.862,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Bond corrections: {'C-C': 1, 'C-H': 6, 'C=C': 1}

External symmetry: 1, optical isomers: 1

Geometry:
{'symbols': ('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'), 'isotopes': (12, 12, 12, 1, 1, 1, 1, 1, 1), 'coords': ((-1.091212, -0.306822, -0.147969), (0.142987, 0.309266, 0.443127), (0.615535, 0.064863, 1.665168), (-0.857297, -0.861767, -1.065333), (-1.572133, -0.995935, 0.552631), (-1.82346, 0.461384, -0.42709), (0.679907, 1.010116, -0.197058), (1.516552, 0.543951, 2.035914), (0.11667, -0.625591, 2.341407))}
""",
)

entry(
    index = 3,
    label = "C3H6_1",
    molecule = 
"""
multiplicity 3
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u1 p0 c0 {1,S} {6,S} {7,S}
3 C u1 p0 c0 {1,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.83352,0.014184,4.6734e-05,-1.2723e-07,1.04135e-10,33791.3,8.35983], Tmin=(10,'K'), Tmax=(314.027,'K')),
            NASAPolynomial(coeffs=[2.81659,0.0271374,-1.51401e-05,4.12627e-09,-4.39346e-13,33855.1,12.088], Tmin=(314.027,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (280.942,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (207.862,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Bond corrections: {'C-C': 2, 'C-H': 6}

External symmetry: 2, optical isomers: 2

Geometry:
{'symbols': ('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'), 'isotopes': (12, 12, 12, 1, 1, 1, 1, 1, 1), 'coords': ((0.070108, -0.421727, -0.150623), (0.961618, 0.725147, -0.497548), (-1.072629, -0.053636, 0.737708), (-0.297366, -0.889336, -1.086486), (0.650356, -1.237209, 0.326701), (1.065572, 1.565562, 0.180896), (1.643324, 0.662222, -1.339261), (-1.482524, 0.950859, 0.718581), (-1.617452, -0.816145, 1.284435))}
""",
)

entry(
    index = 4,
    label = "C3H5_1",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {3,D} {4,S}
2 C u1 p0 c0 {1,S} {7,S} {8,S}
3 C u0 p0 c0 {1,D} {5,S} {6,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {3,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.11736,-0.0114127,0.000142317,-2.61881e-07,1.59394e-10,19197.7,6.75381], Tmin=(10,'K'), Tmax=(507.949,'K')),
            NASAPolynomial(coeffs=[1.09606,0.0273951,-1.66267e-05,4.92428e-09,-5.64706e-13,19310.9,17.3765], Tmin=(507.949,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (159.611,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (182.918,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Bond corrections: {'C-C': 1, 'C-H': 5, 'C=C': 1}

External symmetry: 2, optical isomers: 1

Geometry:
{'symbols': ('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H'), 'isotopes': (12, 12, 12, 1, 1, 1, 1, 1), 'coords': ((0.093707, -0.438417, 0.300451), (-1.232035, -0.037106, 0.280018), (1.170619, 0.351266, -0.067967), (0.30604, -1.454536, 0.632204), (2.187902, -0.020744, -0.030178), (1.025951, 1.372505, -0.407513), (-1.509543, 0.962649, -0.040285), (-2.032357, -0.703066, 0.581109))}
""",
)

entry(
    index = 5,
    label = "C3H5_2",
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
            NASAPolynomial(coeffs=[3.92566,0.00538718,4.62523e-05,-7.73624e-08,4.05623e-11,29256.4,7.45522], Tmin=(10,'K'), Tmax=(495.277,'K')),
            NASAPolynomial(coeffs=[1.46093,0.0252931,-1.40349e-05,3.78709e-09,-3.99374e-13,29500.6,17.6143], Tmin=(495.277,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (243.233,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (182.918,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Bond corrections: {'C-C': 1, 'C-H': 5, 'C=C': 1}

External symmetry: 1, optical isomers: 1

Geometry:
{'symbols': ('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H'), 'isotopes': (12, 12, 12, 1, 1, 1, 1, 1), 'coords': ((-1.122861, -0.070107, 0.242567), (1.411558, 0.091788, -0.404162), (0.10938, 0.026711, -0.555895), (-0.898668, -0.098124, 1.321326), (-1.68812, -0.976443, -0.004283), (-1.784002, 0.785445, 0.061406), (1.877923, 0.080079, 0.587759), (2.094783, 0.160488, -1.24877))}
""",
)

entry(
    index = 6,
    label = "C3H6_2",
    molecule = 
"""
1 C u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
2 C u0 p0 c0 {3,S} {7,S} {8,S} {9,S}
3 C u0 p1 c0 {1,S} {2,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
9 H u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.80645,0.0216817,-5.90914e-05,1.89649e-07,-1.82158e-10,35079.5,7.37548], Tmin=(10,'K'), Tmax=(417.329,'K')),
            NASAPolynomial(coeffs=[0.710028,0.0307386,-1.75244e-05,4.84307e-09,-5.21043e-13,35517.6,21.7595], Tmin=(417.329,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (291.661,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (207.862,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Bond corrections: {'C-C': 2, 'C-H': 6}

External symmetry: 1, optical isomers: 1

Geometry:
{'symbols': ('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'), 'isotopes': (12, 12, 12, 1, 1, 1, 1, 1, 1), 'coords': ((1.213909, 0.091937, 0.012913), (-1.228068, -0.098459, -0.260403), (0.01353, -0.750521, 0.212475), (1.085631, 1.180325, -0.128724), (1.996076, -0.113661, 0.752779), (1.624248, -0.321356, -0.928821), (-1.240169, 0.999101, -0.388129), (-2.113162, -0.434812, 0.291777), (-1.351996, -0.552555, -1.262679))}
""",
)

entry(
    index = 7,
    label = "C3H4_0",
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
            NASAPolynomial(coeffs=[3.98165,0.00126166,6.72715e-05,-1.39195e-07,9.57736e-11,24046.8,4.88762], Tmin=(10,'K'), Tmax=(370.984,'K')),
            NASAPolynomial(coeffs=[2.16018,0.0209018,-1.21426e-05,3.51913e-09,-4.02965e-13,24181.9,11.8689], Tmin=(370.984,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (199.936,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (157.975,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Bond corrections: {'C#C': 1, 'C-C': 1, 'C-H': 4}

External symmetry: 3, optical isomers: 1

Geometry:
{'symbols': ('C', 'C', 'C', 'H', 'H', 'H', 'H'), 'isotopes': (12, 12, 12, 1, 1, 1, 1), 'coords': ((-0.953096, 0.029754, 0.003462), (0.505776, -0.015789, -0.001837), (1.712662, -0.053466, -0.006221), (-1.378383, -0.90276, -0.38299), (-1.339537, 0.178629, 1.017511), (-1.324516, 0.850327, -0.619837), (2.777249, -0.0867, -0.010089))}
""",
)

entry(
    index = 8,
    label = "C3H3_0",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,D} {4,S} {5,S}
2 C u0 p0 c0 {1,D} {3,D}
3 C u1 p0 c0 {2,D} {6,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {3,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.93842,0.00386959,6.10992e-05,-1.37628e-07,9.19255e-11,42279.9,5.50457], Tmin=(10,'K'), Tmax=(515.253,'K')),
            NASAPolynomial(coeffs=[4.24546,0.014702,-8.91006e-06,2.73485e-09,-3.31295e-13,42072.8,2.5245], Tmin=(515.253,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (351.517,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (133.032,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Bond corrections: {'C#C': 1, 'C-C': 1, 'C-H': 3}

External symmetry: 2, optical isomers: 1

Geometry:
{'symbols': ('C', 'C', 'C', 'H', 'H', 'H'), 'isotopes': (12, 12, 12, 1, 1, 1), 'coords': ((-1.077564, -0.018332, -0.125497), (0.282058, 0.004799, 0.03285), (1.502272, 0.025558, 0.174961), (-1.61217, 0.848381, -0.500242), (-1.654569, -0.903957, 0.119785), (2.559973, 0.043552, 0.298145))}
""",
)

entry(
    index = 9,
    label = "C3H4_2",
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
            NASAPolynomial(coeffs=[4.06307,-0.00656032,0.000105529,-2.03192e-07,1.30748e-10,23832.6,4.89256], Tmin=(10,'K'), Tmax=(463.005,'K')),
            NASAPolynomial(coeffs=[1.59406,0.022459,-1.33953e-05,3.91037e-09,-4.43737e-13,23978.8,14.0129], Tmin=(463.005,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (198.152,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (157.975,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Bond corrections: {'C-H': 4, 'C=C': 2}

External symmetry: 4, optical isomers: 1

Geometry:
{'symbols': ('C', 'C', 'C', 'H', 'H', 'H', 'H'), 'isotopes': (12, 12, 12, 1, 1, 1, 1), 'coords': ((1.306454, 0.008864, 0.002277), (-1.306453, -0.008841, -0.002265), (1e-06, 2e-06, -2e-06), (1.867663, 0.699263, 0.627218), (1.879153, -0.673861, -0.620668), (-1.876443, 0.611344, -0.689754), (-1.870375, -0.636771, 0.683194))}
""",
)

