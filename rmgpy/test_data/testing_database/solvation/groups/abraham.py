#!/usr/bin/env python
# encoding: utf-8

name = "Abraham Solute Descriptors"
shortDesc = u""
longDesc = u"""

"""
entry(
    index = -3,
    label = "R",
    group = 
"""
1 * R u0
""",
    solute = None,
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = -2,
    label = "C",
    group = 
"""
1 * C u0
""",
    solute = None,
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = 8,
    label = "Cbf",
    group = 
"""
1 * Cbf u0
""",
    solute = SoluteData(
        S = 0.121,
        B = 0.019,
        E = 0.3,
        L = 0.744,
        A = 0,
    ),
    shortDesc = u"""Platts' fragment 8 fused aromatic""",
    longDesc = 
u"""

""",
)

entry(
    index = -4,
    label = "O",
    group = 
"""
1 * [Os,Od,Ot] u0
""",
    solute = None,
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = -10,
    label = "Oss",
    group = 
"""
1 * Os u0
""",
    solute = None,
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = 26,
    label = "OssH",
    group = 
"""
1 * Os u0 {2,S}
2   H  u0 {1,S}
""",
    solute = SoluteData(
        S = 0.247,
        B = 0.307,
        E = 0.061,
        L = 0.672,
        A = 0,
    ),
    shortDesc = u"""Platts fragment 26 -OH""",
    longDesc = 
u"""

""",
)

entry(
    index = -8,
    label = "N",
    group = 
"""
1 * N u0
""",
    solute = None,
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = -7,
    label = "N3s",
    group = 
"""
1 * N3s u0
""",
    solute = None,
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = 23,
    label = "NO2",
    group = 
"""
1 * N3s u0 {2,S} {3,S} {4,S}
2   R   u0 {1,S}
3   Os  u0 {1,S}
4   Os  u0 {1,S}
""",
    solute = SoluteData(
        S = 0.0,
        B = -0.476,
        E = 0.2,
        L = 0.278,
        A = 0.0,
    ),
    shortDesc = u"""Platts fragment 23 -NO2""",
    longDesc = 
u"""

""",
)

entry(
    index = 36,
    label = "S",
    group = 
"""
1 * [S,Ss,Sd,Sa] u0
""",
    solute = SoluteData(
        S = 0.643,
        B = 0.0,
        E = 0.465,
        L = 0.554,
        A = 0,
    ),
    shortDesc = u"""Platts fragment 36 (any other sulfur)""",
    longDesc = 
u"""

""",
)


tree(
"""
L1: R
    L2: C
        L3: Cbf
    L2: O
        L3: Oss
            L4: OssH
    L2: N
        L3: N3s
            L4: NO2
    L2: S
"""
)

