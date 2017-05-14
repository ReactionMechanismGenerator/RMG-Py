#!/usr/bin/env python
# encoding: utf-8

name = "Radical Groups"
shortDesc = u"Radical corrections to A"
longDesc = u"""
H-bonding parameter A should be modified for when we saturate 
radical molecules with hydrogens and look up the saturated
structure.
"""
entry(
    index = 0,
    label = "R_rad",
    group = 
"""
1 * R u1
""",
    solute = None,
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = 1,
    label = "O_rad",
    group = 
"""
1 * O u1 p2
""",
    solute = None,
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = 2,
    label = "ROJ",
    group = 
"""
1 * O u1 p2 c0 {2,S}
2   R u0 {1,S}
""",
    solute = SoluteData(
        S = 0.0,
        B = 0.0,
        E = 0.0,
        L = 0.0,
        A = -0.345,
    ),
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = 5,
    label = "N3s_rad",
    group = 
"""
1 * N3s u1 p1 
""",
    solute = SoluteData(
        S = 0.0,
        B = 0.0,
        E = 0.0,
        L = 0.0,
        A = -0.087,
    ),
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = 14,
    label = "N3d_guanidine",
    group = 
"""
1   Cd  u0 {2,D} {3,S} {4,S}
2 * N3d u1 {1,D}
3   N3s u0 {1,S} {5,S} {6,S}
4   N3s u0 {1,S} {7,S} {8,S}
5   H   u0 {3,S}
6   H   u0 {3,S}
7   H   u0 {4,S}
8   H   u0 {4,S}
""",
    solute = SoluteData(
        S = 0.0,
        B = 0.0,
        E = 0.0,
        L = 0.0,
        A = -0.17
    ),
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)
tree(
"""
L1: R_rad
    L2: O_rad
        L3: ROJ
    L2: N3s_rad
    L2: N3d_guanidine
"""
)

