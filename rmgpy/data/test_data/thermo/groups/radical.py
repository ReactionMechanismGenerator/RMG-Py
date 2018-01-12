#!/usr/bin/env python
# encoding: utf-8

name = "Radical Corrections"
shortDesc = u""
longDesc = u"""

"""
entry(
    index = 0,
    label = "Radical",
    group = "OR{RJ, RJ2_singlet}",
    thermo = u'RJ',
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = 1,
    label = "RJ",
    group = 
"""
1 * R u1
""",
    thermo = u'CJ',
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = 2,
    label = "CJ",
    group = 
"""
1 * C u1
""",
    thermo = u'CsJ',
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = 3,
    label = "CsJ",
    group = 
"""
1 * Cs u1
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([0.71,0.34,-0.33,-1.07,-2.43,-3.54,-5.43],'cal/(mol*K)'),
        H298 = (104.81,'kcal/mol','+|-',0.1),
        S298 = (0.52,'cal/(mol*K)'),
    ),
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = 94,
    label = "OJ",
    group = 
"""
1 * O u1
""",
    thermo = u'RJ',
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)


entry(
    index = 106,
    label = "RJ2_triplet",
    group = 
"""
1 * R u2
""",
    thermo = u'CsJ',
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)
tree(
"""
L1: Radical
    L2: RJ
        L3: CJ
            L4: CsJ
        L3: OJ
    L2: RJ2_triplet

"""
)

