#!/usr/bin/env python
# encoding: utf-8

name = "intra_H_migration/groups"
shortDesc = u""
longDesc = u"""

"""

template(reactants=["RnH"], products=["RnH"], ownReverse=True)

recipe(actions=[
    ['BREAK_BOND', '*2', 1, '*3'],
    ['FORM_BOND', '*1', 1, '*3'],
    ['GAIN_RADICAL', '*2', '1'],
    ['LOSE_RADICAL', '*1', '1'],
])

boundaryAtoms = ["*1", "*2"]

entry(
    index = 0,
    label = "RnH",
    group = "OR{R5Hall, R6Hall}",
    kinetics = None,
)

entry(
    index = 1,
    label = "Y_rad_out",
    group = 
"""
1 *1 R!H u1
""",
    kinetics = None,
)

entry(
    index = 2,
    label = "XH_out",
    group = 
"""
1 *2 R!H u0 {2,S}
2 *3 H   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 3,
    label = "R5Hall",
    group = 
"""
1 *1 R!H u1 {2,[S,D,T,B]}
2 *4 R!H ux {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *6 R!H ux {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *5 R!H ux {3,[S,D,T,B]} {5,[S,D,T,B]}
5 *2 R!H u0 {4,[S,D,T,B]} {6,S}
6 *3 H   u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 4,
    label = "R5HJ_1",
    group = 
"""
1 *1 R!H u1 {2,[S,D,T,B]}
2 *4 R!H u1 {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *6 R!H u0 {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *5 R!H u0 {3,[S,D,T,B]} {5,[S,D,T,B]}
5 *2 R!H u0 {4,[S,D,T,B]} {6,S}
6 *3 H   u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 5,
    label = "R5HJ_2",
    group = 
"""
1 *1 R!H u1 {2,[S,D,T,B]}
2 *4 R!H u0 {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *6 R!H u1 {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *5 R!H u0 {3,[S,D,T,B]} {5,[S,D,T,B]}
5 *2 R!H u0 {4,[S,D,T,B]} {6,S}
6 *3 H   u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 6,
    label = "R5HJ_3",
    group = 
"""
1 *1 R!H u1 {2,[S,D,T,B]}
2 *4 R!H u0 {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *6 R!H u0 {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *5 R!H u1 {3,[S,D,T,B]} {5,[S,D,T,B]}
5 *2 R!H u0 {4,[S,D,T,B]} {6,S}
6 *3 H   u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 7,
    label = "R5H",
    group = 
"""
1 *1 R!H u1 {2,[S,D,T,B]}
2 *4 R!H u0 {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *6 R!H u0 {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *5 R!H u0 {3,[S,D,T,B]} {5,[S,D,T,B]}
5 *2 R!H u0 {4,[S,D,T,B]} {6,S}
6 *3 H   u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 8,
    label = "R6Hall",
    group = 
"""
1 *1 R!H u1 {2,[S,D,T,B]}
2 *4 R!H ux {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *6 R!H ux {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *7 R!H ux {3,[S,D,T,B]} {5,[S,D,T,B]}
5 *5 R!H ux {4,[S,D,T,B]} {6,[S,D,T,B]}
6 *2 R!H u0 {5,[S,D,T,B]} {7,S}
7 *3 H   u0 {6,S}
""",
    kinetics = None,
)

entry(
    index = 9,
    label = "R6HJ_1",
    group = 
"""
1 *1 R!H u1 {2,[S,D,T,B]}
2 *4 R!H u1 {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *6 R!H u0 {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *7 R!H u0 {3,[S,D,T,B]} {5,[S,D,T,B]}
5 *5 R!H u0 {4,[S,D,T,B]} {6,[S,D,T,B]}
6 *2 R!H u0 {5,[S,D,T,B]} {7,S}
7 *3 H   u0 {6,S}
""",
    kinetics = None,
)

entry(
    index = 10,
    label = "R6HJ_2",
    group = 
"""
1 *1 R!H u1 {2,[S,D,T,B]}
2 *4 R!H u0 {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *6 R!H u1 {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *7 R!H u0 {3,[S,D,T,B]} {5,[S,D,T,B]}
5 *5 R!H u0 {4,[S,D,T,B]} {6,[S,D,T,B]}
6 *2 R!H u0 {5,[S,D,T,B]} {7,S}
7 *3 H   u0 {6,S}
""",
    kinetics = None,
)

entry(
    index = 11,
    label = "R6HJ_3",
    group = 
"""
1 *1 R!H u1 {2,[S,D,T,B]}
2 *4 R!H u0 {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *6 R!H u0 {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *7 R!H u1 {3,[S,D,T,B]} {5,[S,D,T,B]}
5 *5 R!H u0 {4,[S,D,T,B]} {6,[S,D,T,B]}
6 *2 R!H u0 {5,[S,D,T,B]} {7,S}
7 *3 H   u0 {6,S}
""",
    kinetics = None,
)

entry(
    index = 12,
    label = "R6HJ_4",
    group = 
"""
1 *1 R!H u1 {2,[S,D,T,B]}
2 *4 R!H u0 {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *6 R!H u0 {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *7 R!H u0 {3,[S,D,T,B]} {5,[S,D,T,B]}
5 *5 R!H u1 {4,[S,D,T,B]} {6,[S,D,T,B]}
6 *2 R!H u0 {5,[S,D,T,B]} {7,S}
7 *3 H   u0 {6,S}
""",
    kinetics = None,
)

entry(
    index = 13,
    label = "R6H",
    group = 
"""
1 *1 R!H u1 {2,[S,D,T,B]}
2 *4 R!H u0 {1,[S,D,T,B]} {3,[S,D,T,B]}
3 *6 R!H u0 {2,[S,D,T,B]} {4,[S,D,T,B]}
4 *7 R!H u0 {3,[S,D,T,B]} {5,[S,D,T,B]}
5 *5 R!H u0 {4,[S,D,T,B]} {6,[S,D,T,B]}
6 *2 R!H u0 {5,[S,D,T,B]} {7,S}
7 *3 H   u0 {6,S}
""",
    kinetics = None,
)

entry(
    index = 14,
    label = "O_rad_out",
    group = 
"""
1 *1 O u1
""",
    kinetics = None,
)

entry(
    index = 15,
    label = "Cd_rad_out",
    group = 
"""
1 *1 Cd u1
""",
    kinetics = None,
)

entry(
    index = 16,
    label = "C_rad_out_single",
    group = 
"""
1 *1 C u1 {2,S} {3,S}
2    R u0 {1,S}
3    R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 17,
    label = "O_H_out",
    group = 
"""
1 *2 O u0 {2,S}
2 *3 H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 18,
    label = "Cd_H_out_double",
    group = 
"""
1 *2 Cd         u0 {2,S} {3,D}
2 *3 H          u0 {1,S}
3    [Cd,Cdd,O] u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 19,
    label = "Cd_H_out_single",
    group = 
"""
1 *2 Cd u0 {2,S} {3,S}
2 *3 H  u0 {1,S}
3    R  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 20,
    label = "Cs_H_out",
    group = 
"""
1 *2 Cs u0 {2,S} {3,S} {4,S}
2 *3 H  u0 {1,S}
3    R  u0 {1,S}
4    R  u0 {1,S}
""",
    kinetics = None,
)

tree(
"""
L1: RnH
    L2: R5Hall
        L3: R5HJ_1
        L3: R5HJ_2
        L3: R5HJ_3
        L3: R5H
    L2: R6Hall
        L3: R6HJ_1
        L3: R6HJ_2
        L3: R6HJ_3
        L3: R6HJ_4
        L3: R6H
L1: Y_rad_out
    L2: O_rad_out
    L2: Cd_rad_out
    L2: C_rad_out_single
L1: XH_out
    L2: O_H_out
    L2: Cd_H_out_double
    L2: Cd_H_out_single
    L2: Cs_H_out
"""
)

forbidden(
    label = "[CH2]C1=CC(C)CC=C1_1",
    group = 
"""
1 *5 C u0 {2,S} {3,S} {8,S}
2 *2 C u0 {1,S} {9,S}
3    C u0 {1,S} {4,S}
4    C u0 {3,S} {5,D}
5    C u0 {4,D} {6,S}
6 *4 C u0 {5,S} {7,S} {8,D}
7 *1 C u1 {6,S}
8 *6 C u0 {1,S} {6,D}
9 *3 H u0 {2,S}
""",
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

forbidden(
    label = "[CH2]C1=CC(C)CC=C1_2",
    group = 
"""
1 *5 C u0 {2,S} {3,S} {8,S}
2    C u0 {1,S}
3 *2 C u0 {1,S} {4,S} {9,S}
4    C u0 {3,S} {5,D}
5    C u0 {4,D} {6,S}
6 *4 C u0 {5,S} {7,S} {8,D}
7 *1 C u1 {6,S}
8 *6 C u0 {1,S} {6,D}
9 *3 H u0 {3,S}
""",
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

forbidden(
    label = "[CH2]C1=CC(C)CC=C1_3",
    group = 
"""
1    C u0 {2,S} {3,S} {8,S}
2    C u0 {1,S}
3 *2 C u0 {1,S} {4,S} {9,S}
4 *5 C u0 {3,S} {5,D}
5 *6 C u0 {4,D} {6,S}
6 *4 C u0 {5,S} {7,S} {8,D}
7 *1 C u1 {6,S}
8    C u0 {1,S} {6,D}
9 *3 H u0 {3,S}
""",
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

forbidden(
    label = "bridged56_1243",
    group = 
"""
1 *1 C u1 {2,S} {6,S}
2 *4 C u0 {1,S} {3,S} {7,S}
3 *2 C u0 {2,S} {4,S} {8,S}
4 *5 C u0 {3,S} {5,S}
5    C u0 {4,S} {6,S} {7,S}
6    C u0 {1,S} {5,S}
7    C u0 {2,S} {5,S}
8 *3 H u0 {3,S}
""",
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

forbidden(
    label = "bridged56_1254",
    group = 
"""
1 *1 C u1 {2,S} {6,S}
2 *4 C u0 {1,S} {3,S} {7,S}
3    C u0 {2,S} {4,S}
4 *2 C u0 {3,S} {5,S} {8,S}
5 *5 C u0 {4,S} {6,S} {7,S}
6    C u0 {1,S} {5,S}
7    C u0 {2,S} {5,S}
8 *3 H u0 {4,S}
""",
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

