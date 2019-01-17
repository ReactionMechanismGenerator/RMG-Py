#!/usr/bin/env python
# encoding: utf-8

name = "Non Atom Centered Platts Groups"
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
    label = "CO",
    group = 
"""
1 * CO u0
""",
    solute = None,
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = 1,
    label = "Oss(CdsOd)",
    group = 
"""
1 * CO                         u0 {2,S} {3,S} {4,D}
2   Os                         u0 {1,S} {5,S}
3   [Cs,Cd,Cdd,Ct,Cb,Cbf,CO,H] u0 {1,S}
4   Od                         u0 {1,D}
5   R!H                        u0 {2,S}
""",
    solute = SoluteData(
        S = -0.225,
        B = -0.206,
        E = -0.113,
        L = -0.39,
        A = 0,
    ),
    shortDesc = u"""Platts fragment 43 non-cyclic ester""",
    longDesc = 
u"""

""",
)

entry(
    index = 4,
    label = "Cs(OssH)Cs(OssH)",
    group = 
"""
1  * Cs u0 {2,S} {3,S} {4,S} {5,S}
2    R  u0 {1,S}
3    R  u0 {1,S}
4    Os u0 {1,S} {6,S}
5    Cs u0 {1,S} {7,S} {8,S} {9,S}
6    H  u0 {4,S}
7    R  u0 {5,S}
8    R  u0 {5,S}
9    Os u0 {5,S} {10,S}
10   H  u0 {9,S}
""",
    solute = SoluteData(
        S = 0.026,
        B = 0,
        E = -0.0215,
        L = 0.05,
        A = 0,
    ),
    shortDesc = u"""Platts fragment 68 1,2 diol""",
    longDesc = 
u"""

""",
)

entry(
    index = 6,
    label = "CbCsOssH",
    group = 
"""
1 * Cb u0 {2,S}
2   Cs u0 {1,S} {3,S}
3   Os u0 {2,S} {4,S}
4   H  u0 {3,S}
""",
    solute = SoluteData(
        S = 0,
        B = 0.131,
        E = 0,
        L = -0.145,
        A = 0,
    ),
    shortDesc = u"""Platts fragment 79 benzyl alcohol""",
    longDesc = 
u"""

""",
)

entry(
    index = 8,
    label = "OssH",
    group = 
"""
1 * Os                 u0 {2,S} {3,S}
2   H                  u0 {1,S}
3   [Cs,Cd,Ct,CO,Os,H] u0 {1,S}
""",
    solute = SoluteData(
        S = 0,
        B = 0,
        E = 0,
        L = 0,
        A = 0.345,
    ),
    shortDesc = u"""-OH (connected to aliphatic) correction for A""",
    longDesc = 
u"""

""",
)

entry(
    index = 7,
    label = "phenol",
    group = 
"""
1 * Os       u0 {2,S} {3,S}
2   H        u0 {1,S}
3   [Cb,Cbf] u0 {1,S}
""",
    solute = SoluteData(
        S = 0,
        B = 0,
        E = 0,
        L = 0,
        A = 0.543,
    ),
    shortDesc = u"""phenol correction for A""",
    longDesc = 
u"""

""",
)

entry(
    index = 9,
    label = "OxRing",
    group = "OR{OxR3, OxR4, OxR5, OxR6, OxR7}",
    solute = SoluteData(
        S = 0,
        B = 0.12,
        E = -0.001,
        L = -0.001,
        A = 0,
    ),
    shortDesc = u"""Correction for an Oss group in a ring""",
    longDesc = 
u"""

""",
)

entry(
    index = 30,
    label = "N3sH2-benz",
    group = 
"""
1 * N3s u0 {2,S} {3,S} {4,S}
2   H   u0 {1,S}
3   H   u0 {1,S}
4   Cb  u0 {1,S} {5,B} {9,B}
5   Cb  u0 {4,B} {6,B}
6   Cb  u0 {5,B} {7,B}
7   Cb  u0 {6,B} {8,B}
8   Cb  u0 {7,B} {9,B}
9   Cb  u0 {4,B} {8,B}
""",
    solute = SoluteData(
        S = 0.0,
        B = 0.0,
        E = 0.0,
        L = 0.0,
        A = 0.247,
    ),
    shortDesc = u"""aniline correction for A (fragment 4)""",
    longDesc = 
u"""

""",
)

entry(
    index = 31,
    label = "Cd(Od)NH2",
    group = 
"""
1 * N3s u0 {2,S} {3,S} {4,S}
2   H   u0 {1,S}
3   H   u0 {1,S}
4   CO  u0 {1,S} {5,D}
5   Od  u0 {4,D}
""",
    solute = SoluteData(
        S = 0.0,
        B = 0.0,
        E = 0.0,
        L = 0.0,
        A = 0.275,
    ),
    shortDesc = u"""primary amide correction for A (fragment 10)""",
    longDesc = 
u"""

""",
)

entry(
    index = 32,
    label = "Cd(Od)NHR",
    group = 
"""
1 * N3s u0 {2,S} {3,S} {4,S}
2   H   u0 {1,S}
3   R!H u0 {1,S}
4   CO  u0 {1,S} {5,D}
5   Od  u0 {4,D}
""",
    solute = SoluteData(
        S = 0.0,
        B = 0.0,
        E = 0.0,
        L = 0.0,
        A = 0.281,
    ),
    shortDesc = u"""secondary amide correction for A (fragment 11)""",
    longDesc = 
u"""

""",
)

entry(
    index = 33,
    label = "Cd(Od)NH-arom",
    group = 
"""
1 * N3s      u0 {2,S} {3,S} {4,S}
2   H        u0 {1,S}
3   [Cb,N3b] u0 {1,S}
4   CO       u0 {1,S} {5,D}
5   Od       u0 {4,D}
""",
    solute = SoluteData(
        S = 0.0,
        B = 0.0,
        E = 0.0,
        L = 0.0,
        A = -0.091,
    ),
    shortDesc = u"""aromatic amide correction for A (fragment 12)""",
    longDesc = 
u"""

""",
)

entry(
    index = 34,
    label = "N3sHCd(Od)N3sH",
    group = 
"""
1 * N3s u0 {2,S} {3,S}
2   H   u0 {1,S}
3   CO  u0 {1,S} {4,D} {5,S}
4   Od  u0 {3,D}
5   N3s u0 {3,S} {6,S}
6   H   u0 {5,S}
""",
    solute = SoluteData(
        S = 0.0,
        B = 0.0,
        E = 0.0,
        L = 0.0,
        A = -0.0825,
    ),
    shortDesc = u"""urea correction for A (fragment 14)""",
    longDesc = 
u"""

""",
)

entry(
    index = 35,
    label = "N3sCd(Od)N3sH",
    group = 
"""
1 * N3s u0 {2,S} {3,S} {7,S}
2   R!H u0 {1,S}
3   CO  u0 {1,S} {4,D} {5,S}
4   Od  u0 {3,D}
5   N3s u0 {3,S} {6,S}
6   H   u0 {5,S}
7   R!H u0 {1,S}
""",
    solute = SoluteData(
        S = 0.0,
        B = 0.0,
        E = 0.0,
        L = 0.0,
        A = -0.119,
    ),
    shortDesc = u"""urea correction for A (fragment 15)""",
    longDesc = 
u"""

""",
)

entry(
    index = 36,
    label = "CdsNdNsNs",
    group = 
"""
1 * Cd  u0 {2,D} {3,S} {4,S}
2   N3d u0 {1,D} {5,S}
3   N3s u0 {1,S} {6,S} {7,S}
4   N3s u0 {1,S} {8,S} {9,S}
5   H   u0 {2,S}
6   H   u0 {3,S}
7   H   u0 {3,S}
8   H   u0 {4,S}
9   H   u0 {4,S}
""",
    solute = SoluteData(
        S = 0.0,
        B = 0.0,
        E = 0.0,
        L = 0.0,
        A = 0.17,
    ),
    shortDesc = u"""guanidine correction for A (fragment 17)""",
    longDesc = 
u"""

""",
)

tree(
"""
L1: R
    L2: CO
        L3: Oss(CdsOd)
    L2: Cs(OssH)Cs(OssH)
    L2: CbCsOssH
    L2: OssH
    L2: phenol
    L2: OxRing

    L2: N3sH2-benz
    L2: N3sHCd(Od)N3sH
    L2: Cd(Od)NH2
    L2: Cd(Od)NHR
        L3: Cd(Od)NH-arom
    L2: N3sCd(Od)N3sH
    L2: CdsNdNsNs
"""
)

