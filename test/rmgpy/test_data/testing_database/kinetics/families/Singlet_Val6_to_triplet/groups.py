#!/usr/bin/env python
# encoding: utf-8

name = "Singlet_Val6_to_triplet/groups"
shortDesc = ""
longDesc = """
This family makes the transition from the excited O2 singlet state to the ground triplet state.
It also covers S2 and SO which are isoelectronic with O2 and behave similarly.
This family consists of these three transitions only, and cannot be generalized for structures other than O2, S2, or SO.
"""

template(reactants=["singlet"], products=["triplet"], ownReverse=False)

reverse = None
reversible = False

recipe(
    actions=[
        ["CHANGE_BOND", "*1", -1, "*2"],
        ["GAIN_RADICAL", "*1", "1"],
        ["GAIN_RADICAL", "*2", "1"],
    ]
)

entry(
    index=1,
    label="singlet",
    group="""
1 *1 [O2d,S2d] u0 p2 c0 {2,D}
2 *2 [O2d,S2d] u0 p2 c0 {1,D}
""",
    kinetics=None,
)

entry(
    index=2,
    label="O2",
    group="""
1 *1 O2d u0 p2 c0 {2,D}
2 *2 O2d u0 p2 c0 {1,D}
""",
    kinetics=None,
)

entry(
    index=2,
    label="S2",
    group="""
1 *1 S2d u0 p2 c0 {2,D}
2 *2 S2d u0 p2 c0 {1,D}
""",
    kinetics=None,
)

entry(
    index=2,
    label="SO",
    group="""
1 *1 S2d u0 p2 c0 {2,D}
2 *2 O2d u0 p2 c0 {1,D}
""",
    kinetics=None,
)

tree(
    """
L1: singlet
    L2: O2
    L2: S2
    L2: SO
"""
)
