#!/usr/bin/env python
# encoding: utf-8

name = "polycyclic"
shortDesc = u""
longDesc = u""" 

"""

entry(
	index = 1,
	label = "PolycyclicRing",
	group = 
"""
1 * R u0
""",
	solute = SoluteData(
		S = 0.00086,
		B = 0.04137,
		E = -0.08617,
		L = 0.16153,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 17
B: 17
E: 17
L: 15
A: 17
""",
)

entry(
	index = 2,
	label = "polycyclic_7fused",
	group =  "OR{Strychnine_general}",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 3,
	label = "Strychnine_general",
	group = 
"""
1  * R!H u0 {9,[S,D,B]} {16,[S,D,B]}
2  R!H u0 {7,[S,D,B]} {13,[S,D,B]} {15,[S,D,B]}
3  R!H u0 {6,[S,D,B]} {14,[S,D,B]} {19,[S,D,B]}
4  R!H u0 {6,[S,D,B]} {7,[S,D,B]} {10,[S,D,B]} {18,[S,D,B]}
5  R!H u0 {6,[S,D,B]} {8,[S,D,B]} {9,[S,D,B]}
6  R!H u0 {3,[S,D,B]} {4,[S,D,B]} {5,[S,D,B]}
7  R!H u0 {2,[S,D,B]} {4,[S,D,B]} {11,[S,D,B]}
8  R!H u0 {5,[S,D,B]} {11,[S,D,B]} {17,[S,D,B]}
9  R!H u0 {1,[S,D,B]} {5,[S,D,B]} {12,[S,D,B]}
10 R!H u0 {4,[S,D,B]} {13,[S,D,B]}
11 R!H u0 {7,[S,D,B]} {8,[S,D,B]}
12 R!H u0 {9,[S,D,B]} {14,[S,D,B]}
13 R!H u0 {2,[S,D,B]} {10,[S,D,B]}
14 R!H u0 {3,[S,D,B]} {12,[S,D,B]}
15 R!H u0 {2,[S,D,B]} {17,[S,D,B]}
16 R!H u0 {1,[S,D,B]} {20,[S,D,B]}
17 R!H u0 {8,[S,D,B]} {15,[S,D,B]} {20,[S,D,B]}
18 R!H u0 {4,[S,D,B]} {19,[S,D,B]} {21,[S,D,B]}
19 R!H u0 {3,[S,D,B]} {18,[S,D,B]} {22,[S,D,B]}
20 R!H u0 {16,[S,D,B]} {17,[S,D,B]}
21 R!H u0 {18,[S,D,B]} {23,[S,D,B]}
22 R!H u0 {19,[S,D,B]} {24,[S,D,B]}
23 R!H u0 {21,[S,D,B]} {24,[S,D,B]}
24 R!H u0 {22,[S,D,B]} {23,[S,D,B]}
""",
	solute = u'Strychnine',
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 4,
	label = "Strychnine",
	group = 
"""
1  * O u0 {9,S} {16,S}
2  N u0 {7,S} {13,S} {15,S}
3  N u0 {6,S} {14,S} {19,S}
4  C u0 {6,S} {7,S} {10,S} {18,S}
5  C u0 {6,S} {8,S} {9,S}
6  C u0 {3,S} {4,S} {5,S}
7  C u0 {2,S} {4,S} {11,S}
8  C u0 {5,S} {11,S} {17,S}
9  C u0 {1,S} {5,S} {12,S}
10 C u0 {4,S} {13,S}
11 C u0 {7,S} {8,S}
12 C u0 {9,S} {14,S}
13 C u0 {2,S} {10,S}
14 CO u0 {3,S} {12,S}
15 C u0 {2,S} {17,S}
16 C u0 {1,S} {20,S}
17 C u0 {8,S} {15,S} {20,D}
18 C u0 {4,S} {19,B} {21,B}
19 C u0 {3,S} {18,B} {22,B}
20 C u0 {16,S} {17,D}
21 C u0 {18,B} {23,B}
22 C u0 {19,B} {24,B}
23 C u0 {21,B} {24,B}
24 C u0 {22,B} {23,B}
""",
	solute = SoluteData(
		S = -0.03901,
		B = 0.04494,
		E = 0.24179,
		L = 0.19867,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 2
B: 2
E: 2
L: 2
A: 2
""",
)

entry(
	index = 138,
	label = "s1_3_6",
	group = 
"""
1   R!H u0 {2,[S,D,T,B]} {3,[S,D,T,B]} {4,[S,D,T,B]} {5,[S,D,T,B]}
2   R!H u0 {1,[S,D,T,B]} {3,[S,D,T,B]}
3   R!H u0 {1,[S,D,T,B]} {2,[S,D,T,B]}
4   R!H u0 {1,[S,D,T,B]} {6,[S,D,T,B]}
5   R!H u0 {1,[S,D,T,B]} {7,[S,D,T,B]}
6   R!H u0 {4,[S,D,T,B]} {8,[S,D,T,B]}
7   R!H u0 {5,[S,D,T,B]} {8,[S,D,T,B]}
8 * R!H u0 {6,[S,D,T,B]} {7,[S,D,T,B]}
""",
	solute = u's1_3_6_ane',
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 139,
	label = "s1_3_6_ane",
	group = 
"""
1   R!H u0 {2,S} {3,S} {4,S} {5,S}
2   R!H u0 {1,S} {3,S}
3   R!H u0 {1,S} {2,S}
4   R!H u0 {1,S} {6,S}
5   R!H u0 {1,S} {7,S}
6   R!H u0 {4,S} {8,S}
7   R!H u0 {5,S} {8,S}
8 * R!H u0 {6,S} {7,S}
""",
	solute = SoluteData(
		S = -0.00921,
		B = 0.02340,
		E = -0.03777,
		L = -0.03540,
		A = 0.00000,
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
L1: PolycyclicRing
	L2: polycyclic_7fused
		L3: Strychnine_general
			L4: Strychnine
	L2: s1_3_6
		L3: s1_3_6_ane
"""
)
