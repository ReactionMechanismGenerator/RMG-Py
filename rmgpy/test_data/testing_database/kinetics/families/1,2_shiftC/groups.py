#!/usr/bin/env python
# encoding: utf-8

name = "1,2_shiftC/groups"
shortDesc = u""
longDesc = u"""
1,2-methyl shift from a Carbon to another Carbon.
Could possibly be generalized to include  1,2-methyl shift from an Oxygen to Carbon and vice versa.
1,2-methyl shift from Sulfur to Carbon has separate family (1,2_shiftS).
"""

template(reactants=["CH3CCJ"], products=["CH3CCJ"], ownReverse=True)

recipe(actions=[
    ['BREAK_BOND', '*1', 1, '*2'],
    ['FORM_BOND', '*1', 1, '*3'],
    ['GAIN_RADICAL', '*2', '1'],
    ['LOSE_RADICAL', '*3', '1'],
])

boundaryAtoms = ["*1", "*3"]

entry(
    index = 1,
    label = "CH3CCJ",
    group = 
"""
1 *1 Cs  u0 {2,S} {4,S} {5,S} {6,S}
2 *2 C   u0 {1,S} {3,S}
3 *3 C   u1 {2,S}
4    H   u0 {1,S}
5    H   u0 {1,S}
6    H   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 2,
    label = "CJ",
    group = 
"""
1 *3 C   u1
""",
    kinetics = None,
)

entry(
    index = 22,
    label = "CH3CsCJ",
    group =
"""
1 *1 Cs  u0 {2,S} {4,S} {5,S} {6,S}
2 *2 Cs   u0 {1,S} {3,S}
3 *3 C   u1 {2,S}
4    H   u0 {1,S}
5    H   u0 {1,S}
6    H   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 23,
    label = "CH3Cs(-HH)CJ",
    group =
"""
1 *1 Cs  u0 {2,S} {4,S} {5,S} {6,S}
2 *2 Cs   u0 {1,S} {3,S} {7,S} {8,S}
3 *3 C   u1 {2,S}
4    H   u0 {1,S}
5    H   u0 {1,S}
6    H   u0 {1,S}
7    H   u0 {2,S}
8    H   u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 24,
    label = "CH3Cs(-HR!H)CJ",
    group =
"""
1 *1 Cs  u0 {2,S} {4,S} {5,S} {6,S}
2 *2 Cs   u0 {1,S} {3,S} {7,S} {8,S}
3 *3 C   u1 {2,S}
4    H   u0 {1,S}
5    H   u0 {1,S}
6    H   u0 {1,S}
7    H   u0 {2,S}
8    R!H   u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 25,
    label = "CH3Cs(-R!HR!H)CJ",
    group =
"""
1 *1 Cs  u0 {2,S} {4,S} {5,S} {6,S}
2 *2 Cs   u0 {1,S} {3,S} {7,S} {8,S}
3 *3 C   u1 {2,S}
4    H   u0 {1,S}
5    H   u0 {1,S}
6    H   u0 {1,S}
7    R!H   u0 {2,S}
8    R!H   u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 26,
    label = "CH3CdCJ",
    group =
"""
1 *1 Cs  u0 {2,S} {4,S} {5,S} {6,S}
2 *2 Cd   u0 {1,S} {3,S}
3 *3 C   u1 {2,S}
4    H   u0 {1,S}
5    H   u0 {1,S}
6    H   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 5,
    label = "CdsJ",
    group = 
"""
1 *3 Cd u1 {2,D}
2    C  u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 6,
    label = "CsJ",
    group = 
"""
1 *3 Cs u1 {2,S} {3,S}
2    R  u0 {1,S}
3    R  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 7,
    label = "CsJ-HH",
    group = 
"""
1 *3 Cs u1 {2,S} {3,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 8,
    label = "CsJ-CsH",
    group = 
"""
1 *3 Cs u1 {2,S} {3,S}
2    Cs u0 {1,S}
3    H  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 9,
    label = "CsJ-CsCs",
    group = 
"""
1 *3 Cs u1 {2,S} {3,S}
2    Cs u0 {1,S}
3    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 10,
    label = "CsJ-SsH",
    group = 
"""
1 *3 Cs u1 {2,S} {3,S}
2    Ss u0 {1,S}
3    H  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 11,
    label = "CsJ-SsSs",
    group = 
"""
1 *3 Cs u1 {2,S} {3,S}
2    Ss u0 {1,S}
3    Ss u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 12,
    label = "CsJ-CsSs",
    group = 
"""
1 *3 Cs u1 {2,S} {3,S}
2    Cs u0 {1,S}
3    Ss u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 13,
    label = "CsJ-OneDe",
    group = 
"""
1 *3 Cs               u1 {2,S} {3,S}
2    R                u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 14,
    label = "CsJ-OneDeH",
    group = 
"""
1 *3 Cs               u1 {2,S} {3,S}
2    H                u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 15,
    label = "CsJ-CdH",
    group = 
"""
1 *3 Cs u1 {2,S} {3,S}
2    H  u0 {1,S}
3    Cd u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 16,
    label = "CsJ-OneDeCs",
    group = 
"""
1 *3 Cs            u1 {2,S} {3,S}
2    Cs            u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 17,
    label = "CsJ-CdCs",
    group = 
"""
1 *3 Cs u1 {2,S} {3,S}
2    Cs u0 {1,S}
3    Cd u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 18,
    label = "CsJ-OneDeSs",
    group = 
"""
1 *3 Cs            u1 {2,S} {3,S}
2    Ss            u0 {1,S}
3    [Cd,Ct,Cb,CO] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 19,
    label = "CsJ-CdSs",
    group = 
"""
1 *3 Cs u1 {2,S} {3,S}
2    Ss u0 {1,S}
3    Cd u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 20,
    label = "CsJ-TwoDe",
    group =
"""
1 *3 Cs               u1 {2,S} {3,S}
2    [Cd,Ct,Cb,CO,CS] u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 21,
    label = "CsJ-CdCd",
    group =
"""
1 *3 Cs u1 {2,S} {3,S}
2    Cd u0 {1,S}
3    Cd u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 22,
    label = "CsJ-(CdCdCd)H",
    group =
"""
1 *3 Cs u1 {2,S} {3,S}
2    H  u0 {1,S}
3    Cd u0 {1,S} {4,D}
4    Cd u0 {5,S} {3,D}
5    Cd u0 {4,S}
""",
    kinetics = None,
)

tree(
"""
L1: CH3CCJ
    L2: CH3CsCJ
        L3: CH3Cs(-HH)CJ
        L3: CH3Cs(-HR!H)CJ
        L3: CH3Cs(-R!HR!H)CJ
    L2: CH3CdCJ
L1: CJ
    L2: CdsJ
    L2: CsJ
        L3: CsJ-HH
        L3: CsJ-CsH
        L3: CsJ-CsCs
        L3: CsJ-SsH
        L3: CsJ-SsSs
        L3: CsJ-CsSs
        L3: CsJ-TwoDe
            L4: CsJ-CdCd
        L3: CsJ-OneDe
            L4: CsJ-OneDeH
                L5: CsJ-(CdCdCd)H
                L5: CsJ-CdH
            L4: CsJ-OneDeCs
                L5: CsJ-CdCs
            L4: CsJ-OneDeSs
                L5: CsJ-CdSs
"""
)

