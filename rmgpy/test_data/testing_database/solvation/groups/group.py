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
    index = -9,
    label = "Cb",
    group = 
"""
1 * Cb u0
""",
    solute = None,
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

entry(
    index = 38,
    label = "Cb-H",
    group = 
"""
1 * Cb  u0 {2,B} {3,B} {4,S}
2   R!H u0 {1,B}
3   R!H u0 {1,B}
4   H   u0 {1,S}
""",
    solute = SoluteData(
        S = 0.05,
        B = 0.011,
        E = 0.068,
        L = 0.469,
        A = 0,
    ),
    shortDesc = u"""same as Platts group 6""",
    longDesc = 
u"""

""",
)


entry(
	index = 736,
	label = "O",
	group = 
"""
1 * O u0
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 747,
	label = "O2s",
	group = 
"""
1 * O2s u0
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 768,
	label = "O2s-OsH",
	group = 
"""
1 * O2s u0 {2,S} {3,S}
2   O2s u0 {1,S}
3   H  u0 {1,S}
""",
	solute = SoluteData(
		S = 0.15824,
		B = 0.05902,
		E = 0.12734,
		L = 0.46756,
		A = 0.23700,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 3
B: 3
E: 4
L: 3
A: 6
""",
)

entry(
	index = 807,
	label = "S",
	group = 
"""
1 * S ux
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 812,
	label = "S2s",
	group = 
"""
1 * S2s u0
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 813,
	label = "S2s-CH",
	group = 
"""
1 * S2s u0 {2,S} {3,S}
2   C  u0 {1,S}
3   H  u0 {1,S}
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 814,
	label = "S2s-CsH",
	group = 
"""
1 * S2s u0 {2,S} {3,S}
2   Cs u0 {1,S}
3   H  u0 {1,S}
""",
	solute = SoluteData(
		S = 0.20433,
		B = 0.14840,
		E = 0.22845,
		L = 0.91094,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 37
B: 37
E: 38
L: 37
A: 38
""",
)

entry(
	index = 870,
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
	index = 873,
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
	index = 874,
	label = "N3s-CHH",
	group = 
"""
1 * N3s u0 {2,S} {3,S} {4,S}
2   C   u0 {1,S}
3   H   u0 {1,S}
4   H   u0 {1,S}
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 877,
	label = "N3s-CbHH",
	group = 
"""
1 * N3s u0 {2,S} {3,S} {4,S}
2   Cb  u0 {1,S}
3   H   u0 {1,S}
4   H   u0 {1,S}
""",
	solute = SoluteData(
		S = 0.40793,
		B = 0.16116,
		E = 0.22723,
		L = 0.66930,
		A = 0.23902,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 247
B: 229
E: 257
L: 204
A: 248
""",
)

entry(
	index = 1125,
	label = "P",
	group = 
"""
1 * P u0
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 1126,
	label = "P3s",
	group = 
"""
1 * P3s u0
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 1127,
	label = "P3s-O2sO2sO2s",
	group = 
"""
1 * P3s u0 {2,S} {3,S} {4,S}
2   O2s u0 {1,S}
3   O2s u0 {1,S}
4   O2s u0 {1,S}
""",
	solute = SoluteData(
		S = 0.03939,
		B = 0.06846,
		E = 0.00644,
		L = 0.17857,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 2
B: 2
E: 3
L: 2
A: 3
""",
)


tree(
"""
L1: R
	L2: C
		L3: Cb
			L4: Cb-H
	L2: O
		L3: O2s
			L4: O2s-OsH
	L2: S
		L3: S2s
			L4: S2s-CH
				L5: S2s-CsH
	L2: N
		L3: N3s
			L4: N3s-CHH
				L5: N3s-CbHH
	L2: P
		L3: P3s
			L4: P3s-O2sO2sO2s
"""
)
