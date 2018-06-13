#!/usr/bin/env python
# encoding: utf-8

name = "Singlet_Carbene_Intra_Disproportionation/groups"
shortDesc = u"Convert a singlet carbene to a closed-shell molecule through a concerted 1,2-H shift + 1,2-bond formation"
longDesc = u"""
Reaction site *1 should always be a singlet in this family.
"""

template(reactants=["singletcarbene_CH"], products=["CH_C_unsaturated"], ownReverse=False)

reverse = "SingletCarbenefromMultipleBond"

recipe(actions=[
    ['LOSE_PAIR', '*1', '1'],
    ['FORM_BOND', '*1', 1, '*3'],
    ['BREAK_BOND', '*2', 1, '*3'],
    ['CHANGE_BOND', '*1', 1, '*2'],
])

boundaryAtoms = ["*1", "*2"]

entry(
    index = 1,
    label = "singletcarbene_CH",
    group=
    """
    1 *1 C u0 p1 c0 {2,[S,D]}
    2 *2 C u0 {1,[S,D]} {3,S}
    3 *3 H u0 {2,S}
    """,
    kinetics = None,
)

entry(
    index = 2,
    label = "singletcarbene",
    group =
"""
1 *1 C u0 p1 c0
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 3,
    label = "CH",
    group =
"""
1 *2 C u0 {2,S}
2 *3 H u0 {1,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 4,
    label = "fulvene_backbone",
    group = 
"""
1 *2 C u0 {2,S} {6,S} {7,S}
2 C u0 {1,S} {3,S} {5,D}
3 C u0 {2,S} {4,D}
4 C u0 {3,D} {6,S}
5 C u0 {2,D}
6 *1 C u0 p1 c0 {1,S} {4,S}
7 *3 H u0 {1,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 5,
    label = "benzene_backbone",
    group =
"""
1  *2 C u0 {2,S} {6,S} {7,S}
2  C u0 {1,S} {3,D}
3  C u0 {2,D} {4,S}
4  C u0 {3,S} {5,D}
5  C u0 {4,D} {6,S}
6  *1 C u0 p1 c0 {1,S} {5,S}
7 *3 H u0 {1,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 6,
    label = "CsJ2-C",
    group =
"""
1 *1 C u0 p1 c0 {2,S}
2 *2 C u0 {1,S} {3,S}
3 *3 H u0 {2,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 7,
    label = "CdJ2=C",
    group =
"""
1 *1 C u0 p1 c0 {2,D}
2 *2 C u0 {1,D} {3,S}
3 *3 H u0 {2,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 8,
    label = "CdJ2",
    group =
"""
1 *1 Cd u0 p1 c0
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 9,
    label = "CsJ2H",
    group =
"""
1 *1 Cs u0 p1 c0 {2,S}
2 H u0 {1,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 10,
    label = "CsJ2C",
    group =
"""
1 *1 Cs u0 p1 c0 {2,S}
2 C u0 {1,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 11,
    label = "CsJ2(CsC)",
    group =
"""
1 *1 Cs u0 p1 c0 {2,S}
2 Cs u0 {1,S} {3,S}
3 C u0 {2,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 12,
    label = "CsJ2(C=C)",
    group =
"""
1 *1 Cs u0 p1 c0 {2,S}
2 Cd u0 {1,S} {3,D}
3 C u0 {2,D}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 13,
    label = "CdH2",
    group =
"""
1 *2 Cd u0 {2,S} {3,S}
2 *3 H u0 {1,S}
3 H u0 {1,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 14,
    label = "CdHC",
    group =
"""
1 *2 Cd u0 {2,S} {3,S}
2 *3 H u0 {1,S}
3 C u0 {1,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 15,
    label = "CH3",
    group =
"""
1 *2 Cs u0 {2,S} {3,S} {4,S}
2 *3 H u0 {1,S}
3 H u0 {1,S}
4 H u0 {1,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 16,
    label = "CH2(C)",
    group =
"""
1 *2 Cs u0 {2,S} {3,S} {4,S}
2 *3 H u0 {1,S}
3 H u0 {1,S}
4 C u0 {1,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 17,
    label = "CH2(C=C)",
    group =
"""
1 *2 Cs u0 {2,S} {3,S} {4,S}
2 *3 H u0 {1,S}
3 H u0 {1,S}
4 Cd u0 {1,S} {5,D}
5 C u0 {4,D}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 18,
    label = "CH(C)C",
    group =
"""
1 *2 Cs u0 {2,S} {3,S} {4,S}
2 *3 H u0 {1,S}
3 C u0 {1,S}
4 C u0 {1,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 19,
    label = "CH=C",
    group =
"""
1 *2 Cd u0 {2,S} {3,D}
2 *3 H u0 {1,S}
3 C u0 {1,D}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

tree(
"""
L1: singletcarbene_CH
    L2: fulvene_backbone
    L2: benzene_backbone
    L2: CsJ2-C
    L2: CdJ2=C
L1: singletcarbene
    L2: CdJ2
    L2: CsJ2H
    L2: CsJ2C
        L3: CsJ2(CsC)
        L3: CsJ2(C=C)
L1: CH
    L2: CdH2
    L2: CdHC
    L2: CH3
    L2: CH2(C)
        L3: CH2(C=C)
    L2: CH(C)C
    L2: CH=C
"""
)
