#!/usr/bin/env python
# encoding: utf-8

name = "R_Recombination/groups"
shortDesc = u""
longDesc = u"""

"""

template(reactants=["Y_rad", "Y_rad"], products=["Y_Y"], ownReverse=False)

reverse = "Bond_Dissociation"

recipe(actions=[
    ['FORM_BOND', '*1', 1, '*2'],
    ['LOSE_RADICAL', '*1', '1'],
    ['LOSE_RADICAL', '*2', '1'],
])

entry(
    index = 0,
    label = "Y_rad",
    group = 
"""
1 * R u1
""",
    kinetics = None,
)

tree(
"""
L1: Y_rad
"""
)

forbidden(
    label = "Cl",
    group = 
"""
1 *1 Cl u1
""",
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

forbidden(
    label = "Cl_2",
    group = 
"""
1 *2 Cl u1
""",
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

forbidden(
    label = "O4",
    group = 
"""
1    O u1 {2,S}
2 *1 O u0 {1,S} {3,S}
3 *2 O u0 {2,S} {4,S}
4    O u1 {3,S}
""",
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

