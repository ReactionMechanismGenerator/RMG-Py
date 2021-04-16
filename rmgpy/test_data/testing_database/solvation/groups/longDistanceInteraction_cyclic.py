#!/usr/bin/env python
# encoding: utf-8

name = "longDistanceInteraction_cyclic"
shortDesc = u""
longDesc = u""" 

"""

entry(
	index = -1,
	label = "R",
	group = 
"""
1 *1 R u0
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 1,
	label = "aromatic-ortho",
	group = 
"""
1 *1 Cb u0 {2,B}
2 *2 Cb u0 {1,B}
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 2,
	label = "o_OH",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B}
3    O u0 {1,S} {4,S}
4    H u0 {3,S}
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 3,
	label = "o_OH_OH",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B} {5,S}
3    O u0 {1,S} {4,S}
4    H u0 {3,S}
5    O u0 {2,S} {6,S}
6    H u0 {5,S}
""",
	solute = SoluteData(
		S = -0.01250,
		B = -0.01125,
		E = -0.00857,
		L = -0.04727,
		A = -0.04812,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 45
B: 45
E: 46
L: 45
A: 45
""",
)

entry(
	index = 4,
	label = "o_OH_MeO",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B} {5,S}
3    O u0 {1,S} {4,S}
4    H u0 {3,S}
5    O u0 {2,S} {6,S}
6    C u0 {5,S} {7,S} {8,S} {9,S}
7    H u0 {6,S}
8    H u0 {6,S}
9    H u0 {6,S}
""",
	solute = SoluteData(
		S = -0.04689,
		B = 0.03429,
		E = -0.04782,
		L = -0.08123,
		A = -0.15375,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 27
B: 26
E: 30
L: 30
A: 27
""",
)

entry(
	index = 5,
	label = "o_OH_CHO",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B} {5,S}
3    O u0 {1,S} {4,S}
4    H u0 {3,S}
5    C u0 {2,S} {6,S} {7,D}
6    H u0 {5,S}
7    O u0 {5,D}
""",
	solute = SoluteData(
		S = -0.00219,
		B = -0.08964,
		E = -0.02234,
		L = -0.12099,
		A = -0.12705,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 3
B: 3
E: 3
L: 3
A: 3
""",
)

entry(
	index = 6,
	label = "o_CHO",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B}
3    C u0 {1,S} {4,S} {5,D}
4    H u0 {3,S}
5    O u0 {3,D}
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 7,
	label = "o_CHO_CHO",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B} {5,S}
3    C u0 {1,S} {4,S} {8,D}
4    H u0 {3,S}
5    C u0 {2,S} {6,S} {7,D}
6    H u0 {5,S}
7    O u0 {5,D}
8    O u0 {3,D}
""",
	solute = SoluteData(
		S = -0.02599,
		B = 0.02166,
		E = -0.02119,
		L = 0.05049,
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

entry(
	index = 8,
	label = "o_CHO_CH3",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B} {6,S}
3    C u0 {1,S} {4,S} {5,D}
4    H u0 {3,S}
5    O u0 {3,D}
6    C u0 {2,S} {7,S} {8,S} {9,S}
7    H u0 {6,S}
8    H u0 {6,S}
9    H u0 {6,S}
""",
	solute = SoluteData(
		S = -0.02032,
		B = 0.01155,
		E = 0.02170,
		L = 0.05993,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 4
B: 3
E: 4
L: 4
A: 4
""",
)

entry(
	index = 9,
	label = "o_CHO_MeO",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B} {6,S}
3    C u0 {1,S} {4,S} {5,D}
4    H u0 {3,S}
5    O u0 {3,D}
6    O u0 {2,S} {7,S}
7    C u0 {6,S} {8,S} {9,S} {10,S}
8    H u0 {7,S}
9    H u0 {7,S}
10   H u0 {7,S}
""",
	solute = SoluteData(
		S = -0.00691,
		B = -0.02058,
		E = 0.03712,
		L = 0.03491,
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
	index = 10,
	label = "o_vinyl",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B}
3    C u0 {1,S} {4,S} {5,D}
4    H u0 {3,S}
5    C u0 {3,D} {6,S} {7,S}
6   H u0 {5,S}
7   H u0 {5,S}
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 11,
	label = "o_vinyl_CH3",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B} {8,S}
3    C u0 {1,S} {4,S} {5,D}
4    H u0 {3,S}
5    C u0 {3,D} {6,S} {7,S}
6    H u0 {5,S}
7    H u0 {5,S}
8    C u0 {2,S} {9,S} {10,S} {11,S}
9    H u0 {8,S}
10   H u0 {8,S}
11   H u0 {8,S}
""",
	solute = SoluteData(
		S = 0.00322,
		B = -0.01187,
		E = 0.01917,
		L = -0.01980,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 1
B: 1
E: 2
L: 2
A: 2
""",
)

entry(
	index = 12,
	label = "o_MeO",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B}
3    O u0 {1,S} {4,S}
4    C u0 {3,S} {5,S} {6,S} {7,S}
5    H u0 {4,S}
6    H u0 {4,S}
7    H u0 {4,S}
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 13,
	label = "o_MeO_MeO",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B} {8,S}
3    O u0 {1,S} {4,S}
4    C u0 {3,S} {5,S} {6,S} {7,S}
5    H u0 {4,S}
6    H u0 {4,S}
7    H u0 {4,S}
8    O u0 {2,S} {9,S}
9    C u0 {8,S} {10,S} {11,S} {12,S}
10   H u0 {9,S}
11   H u0 {9,S}
12   H u0 {9,S}
""",
	solute = SoluteData(
		S = 0.00232,
		B = 0.03528,
		E = -0.00533,
		L = 0.02199,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 37
B: 27
E: 38
L: 35
A: 37
""",
)

entry(
	index = 14,
	label = "o_CH3",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B}
3    C u0 {1,S} {4,S} {5,S} {6,S}
4    H u0 {3,S}
5    H u0 {3,S}
6    H u0 {3,S}
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 15,
	label = "o_CH3_CH3",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B} {7,S}
3    C u0 {1,S} {4,S} {5,S} {6,S}
4    H u0 {3,S}
5    H u0 {3,S}
6    H u0 {3,S}
7    C u0 {2,S} {8,S} {9,S} {10,S}
8    H u0 {7,S}
9    H u0 {7,S}
10    H u0 {7,S}
""",
	solute = SoluteData(
		S = 0.00252,
		B = 0.00305,
		E = 0.01519,
		L = 0.05823,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 78
B: 72
E: 83
L: 73
A: 83
""",
)

entry(
	index = 16,
	label = "o_CH3_C2H5",
	group = 
"""
1 *1 Cb u0 {2,B} {3,S}
2 *2 Cb u0 {1,B} {7,S}
3    C u0 {1,S} {4,S} {5,S} {6,S}
4    H u0 {3,S}
5    H u0 {3,S}
6    H u0 {3,S}
7    C u0 {2,S} {8,S} {9,S} {10,S}
8    H u0 {7,S}
9    H u0 {7,S}
10   C u0 {7,S} {11,S} {12,S} {13,S}
11   H u0 {10,S}
12   H u0 {10,S}
13   H u0 {10,S}
""",
	solute = SoluteData(
		S = 0.02407,
		B = -0.01130,
		E = 0.04752,
		L = 0.08493,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 13
B: 12
E: 14
L: 13
A: 14
""",
)

entry(
	index = 17,
	label = "aromatic-meta",
	group = 
"""
1 *1 Cb u0 {2,B}
2    Cb u0 {1,B} {3,B}
3 *2 Cb u0 {2,B}
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 18,
	label = "m_CHO_CHO",
	group = 
"""
1 *1 Cb u0 {2,B} {4,S}
2    Cb u0 {1,B} {3,B}
3 *2 Cb u0 {2,B} {7,S}
4    C  u0 {1,S} {5,D} {6,S}
5    O  u0 {4,D}
6    H  u0 {4,S}
7    C  u0 {3,S} {8,D} {9,S}
8    O  u0 {7,D}
9    H  u0 {7,S}
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 19,
	label = "aromatic-para",
	group = 
"""
1 *1 Cb u0 {2,B}
2    Cb u0 {1,B} {3,B}
3    Cb u0 {2,B} {4,B}
4 *2 Cb u0 {3,B}
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 20,
	label = "p_OH",
	group = 
"""
1 *1 Cb u0 {2,B} {5,S}
2    Cb u0 {1,B} {3,B}
3    Cb u0 {2,B} {4,B}
4 *2 Cb u0 {3,B}
5    O u0 {1,S} {6,S}
6    H u0 {5,S}
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 21,
	label = "p_OH_OH",
	group = 
"""
1 *1 Cb u0 {2,B} {5,S}
2    Cb u0 {1,B} {3,B}
3    Cb u0 {2,B} {4,B}
4 *2 Cb u0 {3,B} {7,S}
5    O u0 {1,S} {6,S}
6    H u0 {5,S}
7    O u0 {4,S} {8,S}
8    H u0 {7,S}
""",
	solute = SoluteData(
		S = -0.01768,
		B = -0.01463,
		E = 0.00979,
		L = -0.07210,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 20
B: 20
E: 20
L: 20
A: 20
""",
)

entry(
	index = 22,
	label = "p_OH_MeO",
	group = 
"""
1 *1 Cb u0 {2,B} {5,S}
2    Cb u0 {1,B} {3,B}
3    Cb u0 {2,B} {4,B}
4 *2 Cb u0 {3,B} {7,S}
5    O u0 {1,S} {6,S}
6    H u0 {5,S}
7    O u0 {4,S} {8,S}
8    C u0 {7,S} {9,S} {10,S} {11,S}
9    H u0 {8,S}
10   H u0 {8,S}
11   H u0 {8,S}
""",
	solute = SoluteData(
		S = 0.05626,
		B = -0.02969,
		E = -0.02441,
		L = -0.06934,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 4
B: 4
E: 6
L: 5
A: 4
""",
)

entry(
	index = 23,
	label = "p_OH_CHO",
	group = 
"""
1 *1 Cb u0 {2,B} {5,S}
2    Cb u0 {1,B} {3,B}
3    Cb u0 {2,B} {4,B}
4 *2 Cb u0 {3,B} {7,S}
5    O u0 {1,S} {6,S}
6    H u0 {5,S}
7    C u0 {4,S} {8,D} {9,S}
8    O u0 {7,D}
9    H u0 {7,S}
""",
	solute = SoluteData(
		S = -0.00084,
		B = -0.02017,
		E = -0.00036,
		L = 0.17408,
		A = 0.03002,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 5
B: 5
E: 14
L: 14
A: 5
""",
)

entry(
	index = 24,
	label = "p_MeO",
	group = 
"""
1 *1 Cb u0 {2,B} {5,S}
2    Cb u0 {1,B} {3,B}
3    Cb u0 {2,B} {4,B}
4 *2 Cb u0 {3,B}
5    O u0 {1,S} {6,S}
6    C u0 {5,S} {7,S} {8,S} {9,S}
7    H u0 {6,S}
8    H u0 {6,S}
9    H u0 {6,S}
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 25,
	label = "p_MeO_MeO",
	group = 
"""
1 *1 Cb u0 {2,B} {5,S}
2    Cb u0 {1,B} {3,B}
3    Cb u0 {2,B} {4,B}
4 *2 Cb u0 {3,B} {7,S}
5    O u0 {1,S} {6,S}
6    C u0 {5,S} {12,S} {13,S} {14,S}
7    O u0 {4,S} {8,S}
8    C u0 {7,S} {9,S} {10,S} {11,S}
9    H u0 {8,S}
10   H u0 {8,S}
11   H u0 {8,S}
12   H u0 {6,S}
13   H u0 {6,S}
14   H u0 {6,S}
""",
	solute = SoluteData(
		S = 0.01741,
		B = -0.02377,
		E = -0.01924,
		L = -0.06033,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 5
B: 5
E: 6
L: 6
A: 6
""",
)

entry(
	index = 26,
	label = "p_MeO_CHO",
	group = 
"""
1 *1 Cb u0 {2,B} {5,S}
2    Cb u0 {1,B} {3,B}
3    Cb u0 {2,B} {4,B}
4 *2 Cb u0 {3,B} {7,S}
5    O u0 {1,S} {6,S}
6    C u0 {5,S} {10,S} {11,S} {12,S}
7    C u0 {4,S} {8,D} {9,S}
8    O u0 {7,D}
9    H u0 {7,S}
10   H u0 {6,S}
11   H u0 {6,S}
12   H u0 {6,S}
""",
	solute = SoluteData(
		S = 0.03761,
		B = -0.04243,
		E = 0.01332,
		L = -0.01369,
		A = 0.00000,
	),
	shortDesc = u"""""",
	longDesc = 
u"""
Number of data used to fit each solute parameter:
S: 4
B: 4
E: 4
L: 4
A: 4
""",
)

entry(
	index = 27,
	label = "p_CHO",
	group = 
"""
1 *1 Cb u0 {2,B} {5,S}
2    Cb u0 {1,B} {3,B}
3    Cb u0 {2,B} {4,B}
4 *2 Cb u0 {3,B}
5    C u0 {1,S} {6,S} {7,D}
6    H u0 {5,S}
7    O u0 {5,D}
""",
	solute = None,
	shortDesc = u"""""",
	longDesc = 
u"""
""",
)

entry(
	index = 28,
	label = "p_CHO_CHO",
	group = 
"""
1 *1 Cb u0 {2,B} {5,S}
2    Cb u0 {1,B} {3,B}
3    Cb u0 {2,B} {4,B}
4 *2 Cb u0 {3,B} {7,S}
5    C u0 {1,S} {6,S} {10,D}
6    H u0 {5,S}
7    C u0 {4,S} {8,D} {9,S}
8    O u0 {7,D}
9    H u0 {7,S}
10   O u0 {5,D}
""",
	solute = SoluteData(
		S = -0.03508,
		B = -0.00120,
		E = -0.00785,
		L = -0.05185,
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
L1: R
	L2: aromatic-ortho
		L3: o_OH
			L4: o_OH_OH
			L4: o_OH_MeO
			L4: o_OH_CHO
		L3: o_CHO
			L4: o_CHO_CHO
			L4: o_CHO_CH3
			L4: o_CHO_MeO
		L3: o_vinyl
			L4: o_vinyl_CH3
		L3: o_MeO
			L4: o_MeO_MeO
		L3: o_CH3
			L4: o_CH3_CH3
			L4: o_CH3_C2H5
	L2: aromatic-meta
		L3: m_CHO_CHO
	L2: aromatic-para
		L3: p_OH
			L4: p_OH_OH
			L4: p_OH_MeO
			L4: p_OH_CHO
		L3: p_MeO
			L4: p_MeO_MeO
			L4: p_MeO_CHO
		L3: p_CHO
			L4: p_CHO_CHO
"""
)