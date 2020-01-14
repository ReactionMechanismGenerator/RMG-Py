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
Species C3H6_3 (run time: 0:00:54)

Overall time since project initiation: 00:01:59
"""
entry(
    index = 0,
    label = "C3H6_3",
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
            NASAPolynomial(coeffs=[4.2074,-0.0154876,0.000128181,-1.78091e-07,8.11606e-11,4872.65,4.96395], Tmin=(10,'K'), Tmax=(678.21,'K')),
            NASAPolynomial(coeffs=[-1.23035,0.0340837,-2.01609e-05,5.77304e-09,-6.39364e-13,5207.76,26.1193], Tmin=(678.21,'K'), Tmax=(3000,'K')),
        ],
        Tmin = (10,'K'),
        Tmax = (3000,'K'),
        E0 = (40.539,'kJ/mol'),
        Cp0 = (33.2579,'J/(mol*K)'),
        CpInf = (207.862,'J/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
Bond corrections: {'C-C': 3, 'C-H': 6}

External symmetry: 6, optical isomers: 1

Geometry:
{'symbols': ('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'), 'isotopes': (12, 12, 12, 1, 1, 1, 1, 1, 1), 'coords': ((0.789858, -0.365811, -0.016399), (-0.078234, 0.863594, 0.07776), (-0.711624, -0.497784, -0.06136), (1.345841, -0.530678, -0.934365), (1.306202, -0.697575, 0.879302), (-0.15116, 1.366362, 1.037377), (-0.111521, 1.533259, -0.77629), (-1.174861, -0.752235, -1.009846), (-1.214501, -0.919133, 0.803821))}
""",
)

