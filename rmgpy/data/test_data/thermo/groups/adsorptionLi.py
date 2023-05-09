#!/usr/bin/env python
# encoding: utf-8

name = "Surface Adsorption Corrections Lithium"
shortDesc = "Li"
longDesc = """
Changes due to adsorbing on a Lithium electrode.
Note: "-h" means "horizontal"

Li enthalpy of formation 159.30 kJ/mol https://webbook.nist.gov/cgi/cbook.cgi?ID=C7439932&Units=SI&Mask=1#ref-2

https://labs.chem.ucsb.edu/zakarian/armen/11---bonddissociationenergy.pdf
LiF = 577
LiH = 247
LiO = 341
LiOH = 427

https://pubs.acs.org/doi/abs/10.1021/om950966x
Li-CH3 = 190.37

LiF
4.329154975994253
LiH
0.9089463523933348
LiO
1.8831875966917782
LiOH
2.7745146925392903
LiC
0.3220178240463045
"""

entry(
    index = 0,
    label = "R*",
    group = 
"""
1   R ux
2 * X ux
""",
    thermo = None,
    shortDesc = """Anything adsorbed anyhow.""",
    longDesc = 
"""
R
   X
***********
This node should be empty, ensuring that one of the nodes below is used.


The group could well be defined as:

    1 R ux
    2 * Xux

but then it is identical with the R*vdW node, and the database tests
do not like that. It should be OK, because things would check the
tree in order, and if there *was* a bond it would match either
R*bidentate or R*single_chemisorbed and thus not R*vdW.
""",
)

entry(
    index = 1,
    label = "R*single_chemisorbed",
    group = 
"""
1 * X u0 {2,S}
2   R ux {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([-0.07,1.05,1.77,2.43,2.8,3.08,3.39],'cal/(mol*K)'),
        H298 = (-2,'eV/molecule'),
        S298 = (-38.17,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 2,
    label = "C*",
    group = 
"""
1 * X u0 {2,S}
2   C ux {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([-0.45,0.61,1.42,2.02,2.81,3.26,3.73],'cal/(mol*K)'),
        H298 = (-0.5,'eV/molecule'),
        S298 = (-32.73,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
   CR3
   |
***********
""",
)

entry(
    index = 3,
    label = "O*",
    group = 
"""
1 * X u0 {2,S}
2   O ux {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([1.09,1.82,2.2,2.42,2.65,2.75,2.86],'cal/(mol*K)'),
        H298 = (-1.88,'eV/molecule'),
        S298 = (-33.89,'cal/(mol*K)'),
    ),
    shortDesc = """Came from OH single-bonded on Pt(111)""",
    longDesc = 
"""
   R
   |
   O
   |
***********
""",
)

entry(
    index = 4,
    label = "F*",
    group = 
"""
1 * X u0 {2,S}
2   F ux {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([1.09,1.82,2.2,2.42,2.65,2.75,2.86],'cal/(mol*K)'),
        H298 = (-4.33,'eV/molecule'),
        S298 = (-26,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
   F
   |
***********
""",

)

entry(
    index = 5,
    label = "H*",
    group = 
"""
1 * X u0 {2,S}
2   H ux {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([1.09,1.82,2.2,2.42,2.65,2.75,2.86],'cal/(mol*K)'),
        H298 = (-0.9,'eV/molecule'),
        S298 = (-26,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
   H
   |
***********
""",
)

entry(
    index = 6,
    label = "OH*",
    group = 
"""
1 * X u0 {2,S}
2   O u0 {1,S} {3,S}
3   H u0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([1.09,1.82,2.2,2.42,2.65,2.75,2.86],'cal/(mol*K)'),
        H298 = (-2.77,'eV/molecule'),
        S298 = (-34,'cal/(mol*K)'),
    ),
    shortDesc = """Came from OH single-bonded on Pt(111)""",
    longDesc = 
"""
   H
   |
   O
   |
***********
""",
)

entry(
    index = 7,
    label = "R*vdW",
    group = 
"""
1 * X u0
2   R u0
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([1.23,1.71,2,2.19,2.39,2.5,2.61],'cal/(mol*K)'),
        H298 = (-0.5,'eV/molecule'),
        S298 = (-20.48,'cal/(mol*K)'),
    ),
    shortDesc = """Average of (CR4)*, (NR3)* and (OR2)* thermo.""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "111",
)

entry(
    index = 8,
    label = "O*vdW",
    group = 
"""
1 * X u0 p0
2   O u0 p2
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0.71,1.22,1.49,1.65,1.81,1.9,1.98],'cal/(mol*K)'),
        H298 = (-0.8,'eV/molecule'),
        S298 = (-22.53,'cal/(mol*K)'),
    ),
    shortDesc = """""",
    longDesc = 
"""
""",
)

tree(
"""
L1: R*
    L2: R*single_chemisorbed
        L3: C*
        L3: O*
            L4: OH*
        L3: F*
        L3: H*
    L2: R*vdW
        L3: O*vdW
"""
)

