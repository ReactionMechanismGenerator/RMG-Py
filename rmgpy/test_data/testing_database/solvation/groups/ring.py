#!/usr/bin/env python
# encoding: utf-8

name = "ring"
shortDesc = u""
longDesc = u""" 

"""

entry(
	index = 1,
	label = "Ring",
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
	index = 2,
	label = "Aromatic",
	group = 
"""
1 * Cb u0
""",
	solute = u'Benzene',
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 3,
	label = "Benzene",
	group = 
"""
1 * Cb u0 {2,B} {6,B}
2   Cb u0 {1,B} {3,B}
3   Cb u0 {2,B} {4,B}
4   Cb u0 {3,B} {5,B}
5   Cb u0 {4,B} {6,B}
6   Cb u0 {1,B} {5,B}
""",
	solute = SoluteData(
		S = -0.08211,
		B = 0.03455,
		E = -0.09106,
		L = -0.06839,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 3037
B: 2925
E: 3195
L: 2736
A: 3192
""",
)

entry(
	index = 4,
	label = "ThreeMember",
	group = 
"""
1 * R!H u0 {2,[S,D,T]} {3,[S,D]}
2   R!H u0 {1,[S,D,T]} {3,[S,D]}
3   R!H u0 {1,[S,D]} {2,[S,D]}
""",
	solute = u'Cyclopropane',
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 5,
	label = "Cyclopropane",
	group = 
"""
1 * Cs u0 {2,S} {3,S}
2   Cs u0 {1,S} {3,S}
3   Cs u0 {1,S} {2,S}
""",
	solute = SoluteData(
		S = 0.08944,
		B = 0.00000,
		E = 0.12286,
		L = -0.02938,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 46
B: 45
E: 45
L: 36
A: 46
""",
)

entry(
	index = 6,
	label = "Ethylene_oxide",
	group = 
"""
1 * O2s    u0 {2,S} {3,S}
2   [Cs,N] u0 {1,S} {3,S}
3   [Cs,N] u0 {1,S} {2,S}
""",
	solute = SoluteData(
		S = 0.26498,
		B = -0.00442,
		E = 0.07866,
		L = 0.07046,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 20
B: 20
E: 22
L: 20
A: 21
""",
)

entry(
	index = 7,
	label = "Ethyleneimine",
	group = 
"""
1 * N       u0 {2,S} {3,S}
2   [Cs,N,S]  u0 {1,S} {3,S}
3   [Cs,N,S]  u0 {1,S} {2,S}
""",
	solute = SoluteData(
		S = 0.07140,
		B = -0.02734,
		E = 0.04383,
		L = 0.07610,
		A = 0.00704,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 1
B: 1
E: 1
L: 1
A: 1
""",
)

tree(
"""
L1: Ring
	L2: Aromatic
		L3: Benzene
	L2: ThreeMember
		L3: Cyclopropane
		L3: Ethylene_oxide
		L3: Ethyleneimine
"""
)
