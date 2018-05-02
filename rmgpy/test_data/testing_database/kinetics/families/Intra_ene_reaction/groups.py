#!/usr/bin/env python
# encoding: utf-8

name = "Intra_ene_reaction/groups"
shortDesc = u"6-membered Intramolecular H-shift from an allylic to an unsaturated endgroup (like in cyclopentadiene)"
longDesc = u"""

"""

template(reactants=["1_3_unsaturated_pentane_backbone"], products=["1_3_unsaturated_pentane_backbone"], ownReverse=True)

recipe(actions=[
    ['BREAK_BOND', '*1', 1, '*6'],
    ['FORM_BOND', '*2', 1, '*6'],
    ['CHANGE_BOND', '*2', -1, '*3'],
    ['CHANGE_BOND', '*4', -1, '*5'],
    ['CHANGE_BOND', '*1', 1, '*5'],
    ['CHANGE_BOND', '*4', 1, '*3'],
])

boundaryAtoms = ["*1", "*2"]

entry(
    index = 1,
    label = "1_3_unsaturated_pentane_backbone",
    group=
    """
    1 *1 C u0 {5,[S,D]} {6,S}
    2 *2 C u0 {3,[D,T]}
    3 *3 C u0 {2,[D,T]} {4,[S,D]}
    4 *4 C u0 {3,[S,D]} {5,[D,T]}
    5 *5 C u0 {1,[S,D]} {4,[D,T]}
    6 *6 H u0 {1,S}
    """,
    kinetics = None,
)

entry(
    index = 2,
    label = "CH_end",
    group =
"""
1 *1 C u0 {2,S}
2 *6 H u0 {1,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 3,
    label = "unsaturated_end",
    group =
"""
1 *2 C u0
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 4,
    label = "cyclopentadiene",
    group =
"""
1 *1 C u0 {2,S} {5,S} {6,S}
2 *2 C u0 {1,S} {3,D}
3 *3 Cd u0 {2,D} {4,S}
4 *4 Cd u0 {3,S} {5,D}
5 *5 Cd u0 {1,S} {4,D}
6 *6 H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 5,
    label = "1_3_4_pentatriene",
    group =
"""
1 *1 C u0 {5,D} {6,S}
2 *2 C u0 {3,D}
3 *3 Cd u0 {2,D} {4,S} {7,S}
4 *4 Cd u0 {3,S} {5,D} {8,S}
5 *5 Cdd u0 {1,D} {4,D}
6 *6 H u0 {1,S}
7 R u0 {3,S}
8 R u0 {4,S}
""",
    kinetics = None,
)

entry(
    index = 6,
    label = "1_3_pentadiene",
    group =
"""
1 *1 C u0 {5,S} {6,S}
2 *2 C u0 {3,D}
3 *3 Cd u0 {2,D} {4,S} {7,S}
4 *4 Cd u0 {3,S} {5,D} {8,S}
5 *5 Cd u0 {1,S} {4,D} {9,S}
6 *6 H u0 {1,S}
7 R u0 {3,S}
8 R u0 {4,S}
9 R u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 7,
    label = "1_pentyn_3_ene",
    group =
"""
1 *1 C u0 {5,S} {6,S}
2 *2 C u0 {3,T}
3 *3 Ct u0 {2,T} {4,S}
4 *4 Cd u0 {3,S} {5,D} {8,S}
5 *5 Cd u0 {1,S} {4,D} {9,S}
6 *6 H u0 {1,S}
8 R u0 {4,S}
9 R u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 8,
    label = "CdH2_1",
    group =
"""
1 *1 Cd u0 {2,S} {3,S}
2 *6 H u0 {1,S}
3 H u0 {1,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 9,
    label = "CdHC_1",
    group =
"""
1 *1 Cd u0 {2,S} {3,S}
2 *6 H u0 {1,S}
3 C u0 {1,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 10,
    label = "CH3_1",
    group =
"""
1 *1 Cs u0 {2,S} {3,S} {4,S}
2 *6 H u0 {1,S}
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
    index = 11,
    label = "CH2(C)_1",
    group =
"""
1 *1 Cs u0 {2,S} {3,S} {4,S}
2 *6 H u0 {1,S}
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
    index = 12,
    label = "CH(C)C_1",
    group =
"""
1 *1 Cs u0 {2,S} {3,S} {4,S}
2 *6 H u0 {1,S}
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
    index = 13,
    label = "CH=C_1",
    group =
"""
1 *1 Cd u0 {2,S} {3,D}
2 *6 H u0 {1,S}
3 C u0 {1,D}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 14,
    label = "CdH2_2",
    group =
"""
1  *2 Cd u0 {2,S} {3,S}
2  H u0 {1,S}
3  H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 15,
    label = "CdHC_2",
    group =
"""
1  *2 Cd u0 {2,S} {3,S}
2  H u0 {1,S}
3  C u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 16,
    label = "Cd(C)C_2",
    group =
"""
1  *2 Cd u0 {2,S} {3,S}
2  C u0 {1,S}
3  C u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 17,
    label = "CddC_2",
    group =
"""
1  *2 Cdd u0 {2,D}
2  C u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 18,
    label = "CtH_2",
    group =
"""
1  *2 Ct u0 {2,S}
2  H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 19,
    label = "CtC_2",
    group =
"""
1  *2 Ct u0 {2,S}
2  C u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 20,
    label = "CdCJ_2",
    group =
"""
1  *2 Cd u0 {2,S}
2  C u1 {1,S}
""",
    kinetics = None,
)

entry(
    index = 21,
    label = "CH(CJ)_1",
    group =
"""
1 *1 Cs u0 {2,S} {3,S}
2 *6 H u0 {1,S}
3 C u1 {1,S}
""",
    kinetics = None,
    shortDesc = "",
    longDesc =
u"""

""",
)

entry(
    index = 22,
    label = "indene",
    group =
"""
1 *1 C u0 {2,S} {5,S} {6,S} {7,S}
2 *2 C u0 {1,S} {3,D}
3 *3 Cd u0 {2,D} {4,S}
4 *4 Cd u0 {3,S} {5,D}
5 *5 Cd u0 {1,S} {4,D} {10,S}
6 *6 H  u0 {1,S}
7    Cd  u0 {1,S} {8,D}
8    Cd  u0 {7,D} {9,S}
9    Cd  u0 {8,S} {10,D}
10   Cd  u0 {9,D} {5,S}
""",
    kinetics = None,
)

entry(
    index = 23,
    label = "cyclopentadiene_cyc6",
    group =
"""
1 *1 C u0 {2,S} {5,S} {6,S}
2 *2 C u0 {1,S} {3,D}
3 *3 Cd u0 {2,D} {4,S} {7,S}
4 *4 Cd u0 {3,S} {5,D} {10,S}
5 *5 Cd u0 {1,S} {4,D}
6 *6 H  u0 {1,S}
7    Cd u0 {3,S} {8,D}
8    Cd u0 {7,D} {9,S}
9    Cd u0 {8,S} {10,D}
10   Cd u0 {9,D} {4,S}
""",
    kinetics = None,
)

tree(
"""
L1: 1_3_unsaturated_pentane_backbone
    L2: cyclopentadiene
        L3: indene
	L3: cyclopentadiene_cyc6
    L2: 1_3_4_pentatriene
    L2: 1_3_pentadiene
    L2: 1_pentyn_3_ene
L1: CH_end
    L2: CdH2_1
    L2: CdHC_1
    L2: CH3_1
    L2: CH2(C)_1
    L2: CH(C)C_1
    L2: CH(CJ)_1
    L2: CH=C_1
L1: unsaturated_end
    L2: CdH2_2
    L2: CdHC_2
    L2: Cd(C)C_2
    L2: CdCJ_2
    L2: CddC_2
    L2: CtH_2
    L2: CtC_2
"""
)

forbidden(
    label = "fulvene_H_shift_ring_edge_to_tail",
    group =
"""
1 *2 C u0 {2,D}
2 *3 C u0 {1,D} {3,S} {4,S}
3 C ux {2,S} {5,S}
4 *4 C u0 {2,S} {6,D}
5 *1 C u0 {3,S} {6,S} {7,S}
6 *5 C u0 {4,D} {5,S}
7 *6 H u0 {5,S}
""",
    shortDesc = u"""""",
    longDesc =
u"""
Prevents an H on the far edge of a fulvene-like ring molecule from shifting to the tail
""",
)

forbidden(
    label = "fulvene_H_shift_tail_to_ring_edge",
    group =
"""
1 *1 C u0 {2,S} {7,S}
2 *5 C u0 {1,S} {3,S} {4,D}
3 C ux {2,S} {5,S}
4 *4 C u0 {2,D} {6,S}
5 *2 C u0 {3,S} {6,D}
6 *3 C u0 {4,S} {5,D}
7 *6 H u0 {1,S}
""",
    shortDesc = u"""""",
    longDesc =
u"""
Prevents an H on the tail of a fulvene-like molecule from shifting to the far edge of the ring
""",
)

forbidden(
    label = "H_shift_to_single_resonant_radical_CPD",
    group =
"""
1 *1 C u0 {5,[S,D]} {6,S} {2,S}
2 *2 C u0 {3,[D,T]} {1,S}
3 *3 C u0 {2,[D,T]} {4,[S,D]} {7,S}
4 *4 C u0 {3,[S,D]} {5,[D,T]}
5 *5 C u0 {1,[S,D]} {4,[D,T]}
6 *6 H u0 {1,S}
7    R!H u1 {3,S}
""",
    shortDesc = u"""""",
    longDesc =
u"""
Forbid an H from shifting to a resonant radical site on a CPD ring,
in order to avoid redundancy with Intra_H_migration family, the radical atom can be anything element
""",
)

forbidden(
    label = "H_shift_to_double_resonant_radical_CPD",
    group =
"""
1 *1 C u0 {5,[S,D]} {6,S} {2,S}
2 *2 C u0 {3,[D,T]} {1,S}
3 *3 C u0 {2,[D,T]} {4,[S,D]}
4 *4 C u0 {3,[S,D]} {5,[D,T]}
5 *5 C u0 {1,[S,D]} {4,[D,T]} {7,S}
6 *6 H u0 {1,S}
7    R!H u1 {5,S}
""",
    shortDesc = u"""""",
    longDesc =
u"""
Forbid an H from shifting to a doubly resonant radical site on a CPD ring,
in order to avoid redundancy with Intra_H_migration family, the radical atom can be anything element
""",
)


forbidden(
    label = "H_shift_to_single_resonant_radical_linear",
    group =
"""
1 *1 C u0 {5,[S,D]} {6,S}
2 *2 C u0 {3,[D,T]}
3 *3 C u0 {2,[D,T]} {4,[S,D]} {7,S}
4 *4 C u0 {3,[S,D]} {5,[D,T]}
5 *5 C u0 {1,[S,D]} {4,[D,T]}
6 *6 H u0 {1,S}
7    R!H u1 {3,S}
""",
    shortDesc = u"""""",
    longDesc =
u"""
Forbid an H from shifting to a resonant radical site on a linear 1,3-unsaturated hexane chain,
in order to avoid redundancy with Intra_H_migration family, the radical atom can be anything element
""",
)

forbidden(
    label = "H_shift_to_double_resonant_radical_linear",
    group =
"""
1 *1 C u0 {5,[S,D]} {6,S}
2 *2 C u0 {3,[D,T]}
3 *3 C u0 {2,[D,T]} {4,[S,D]}
4 *4 C u0 {3,[S,D]} {5,[D,T]}
5 *5 C u0 {1,[S,D]} {4,[D,T]} {7,S}
6 *6 H u0 {1,S}
7    R!H u1 {5,S}
""",
    shortDesc = u"""""",
    longDesc =
u"""
Forbid an H from shifting to a doubly resonant radical site on a linear 1,3-unsaturated hexane chain,
in order to avoid redundancy with Intra_H_migration family, the radical atom can be anything element
""",
)

forbidden(
    label = "para_to_benzylic_shift",
    group =
"""
1  *1 C u0 {2,S} {8,S}
2  *5 C u0 {1,S} {3,D} {7,[S,D]}
3  *4 C u0 {2,D} {4,S}
4  *3 C u0 {3,S} {5,D}
5  *2 C u0 {4,D} {6,[S,D]}
6  R!H ux {5,[S,D]} {7,[D,T]}
7  R!H ux {2,[S,D]} {6,[D,T]}
8  *6 H u0 {1,S}
""",
    shortDesc = u"""""",
    longDesc =
u"""
Prevents an H on para position of a benzene ring from shifting to the benzylic position of a tail
""",
)

forbidden(
    label = "benzylic_to_para_shift",
    group =
"""
1  *2 C u0 {2,D}
2  *3 C u0 {1,D} {3,S} {7,[S,D]}
3  *4 C u0 {2,S} {4,D}
4  *5 C u0 {3,D} {5,S}
5  *1 C u0 {4,S} {6,[S,D]} {8,S}
6  R!H ux {5,[S,D]} {7,[D,T]}
7  R!H ux {2,[S,D]} {6,[D,T]}
8  *6 H u0 {5,S}
""",
    shortDesc = u"""""",
    longDesc =
u"""
Prevents an H on the benzylic position of a tail from shifting to the para position of the benzene ring
""",
)
