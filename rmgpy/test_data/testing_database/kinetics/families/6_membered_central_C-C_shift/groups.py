#!/usr/bin/env python
# encoding: utf-8

name = "6_membered_central_C-C_shift/groups"
shortDesc = u"Concerted shift of the central C-C bond in an 1,5-unsaturated hexane to between the end atoms"
longDesc = u"""
Taken from:

Miller, J. A.; Klippenstein, S. J., The Recombination of Propargyl Radicals and Other Reactions on a C6H6 Potential.
J. Phys. Chem. A 2003, 107, 7783-7799.
"""

template(reactants=["1_5_unsaturated_hexane"], products=["1_5_unsaturated_hexane"], ownReverse=True)

recipe(actions=[
    ['BREAK_BOND', '*3', 1, '*4'],
    ['CHANGE_BOND', '*1', -1, '*2'],
    ['CHANGE_BOND', '*5', -1, '*6'],
    ['CHANGE_BOND', '*2', 1, '*3'],
    ['CHANGE_BOND', '*4', 1, '*5'],
    ['FORM_BOND', '*1', 1, '*6'],
])

boundaryAtoms = ["*1", "*6"]

entry(
    index = 1,
    label = "1_5_unsaturated_hexane",
    group=
    """
    1 *3 C u0 {2,S} {3,[S,D]}
    2 *4 C u0 {1,S} {4,[S,D]}
    3 *2 C u0 {1,[S,D]} {5,[D,T]}
    4 *5 C u0 {2,[S,D]} {6,[D,T]}
    5 *1 C u0 {3,[D,T]}
    6 *6 C u0 {4,[D,T]}
    """,
    kinetics = None,
)

entry(
    index = 2,
    label = "1_5_hexadiyne",
    group=
"""
1 *3 C u0 {2,S} {3,S}
2 *4 C u0 {1,S} {4,S}
3 *2 Ct u0 {1,S} {5,T}
4 *5 Ct u0 {2,S} {6,T}
5 *1 Ct u0 {3,T}
6 *6 Ct u0 {4,T}
""",
    kinetics = None,
)

entry(
    index = 2,
    label = "1_2_4_5_hexatetraene",
    group=
"""
1 *3 Cd u0 {2,S} {3,D}
2 *4 Cd u0 {1,S} {4,D}
3 *2 Cdd u0 {1,D} {5,D}
4 *5 Cdd u0 {2,D} {6,D}
5 *1 C u0 {3,D}
6 *6 C u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 2,
    label = "1_5_hexadiene",
    group=
"""
1 *3 C u0 {2,S} {3,S}
2 *4 C u0 {1,S} {4,S}
3 *2 Cd u0 {1,S} {5,D}
4 *5 Cd u0 {2,S} {6,D}
5 *1 C u0 {3,D}
6 *6 C u0 {4,D}
""",
    kinetics = None,
)

tree(
"""
L1: 1_5_unsaturated_hexane
    L2: 1_5_hexadiyne
    L2: 1_2_4_5_hexatetraene
    L2: 1_5_hexadiene
"""
)

