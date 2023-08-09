#!/usr/bin/env python
# encoding: utf-8

name = "halogen"
shortDesc = ""
longDesc = """ 

"""

entry(
    index=1,
    label="X",
    group="""
1 * [F1s,Cl1s,Br1s,I1s] ux
""",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=2,
    label="F",
    group="""
1 * F1s u0
""",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=3,
    label="F-C",
    group="""
1 * F1s u0 {2,S}
2   C   u0 {1,S}
""",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=4,
    label="F-Cb",
    group="""
1 * F1s u0 {2,S}
2   Cb  u0 {1,S}
""",
    solute=SoluteData(
        S=-0.02870,
        B=-0.01608,
        E=-0.10588,
        L=-0.12204,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 135
B: 131
E: 146
L: 104
A: 150
""",
)

entry(
    index=5,
    label="F-Phenol",
    group="OR{F-Phenol(ortho), F-Phenol(meta), F-Phenol(para)}",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=6,
    label="F-Phenol(ortho)",
    group="""
1 * F1s u0 {2,S}
2   Cb  u0 {1,S} {3,B}
3   Cb  u0 {2,B} {4,S}
4   O2s u0 {3,S} {5,S}
5   H   u0 {4,S}
""",
    solute=SoluteData(
        S=-0.01903,
        B=-0.03461,
        E=-0.10427,
        L=-0.12259,
        A=0.09756,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 5
B: 5
E: 8
L: 5
A: 7
""",
)

entry(
    index=7,
    label="F-Phenol(meta)",
    group="""
1 * F1s u0 {2,S}
2   Cb  u0 {1,S} {3,B}
3   Cb  u0 {2,B} {4,B}
4   Cb  u0 {3,B} {5,S}
5   O2s u0 {4,S} {6,S}
6   H   u0 {5,S}
""",
    solute=SoluteData(
        S=0.02174,
        B=-0.03645,
        E=-0.09013,
        L=-0.02170,
        A=0.07464,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 3
B: 3
E: 5
L: 3
A: 5
""",
)

entry(
    index=8,
    label="F-Phenol(para)",
    group="""
1 * F1s u0 {2,S}
2   Cb  u0 {1,S} {3,B}
3   Cb  u0 {2,B} {4,B}
4   Cb  u0 {3,B} {5,B}
5   Cb  u0 {4,B} {6,S}
6   O2s u0 {5,S} {7,S}
7   H   u0 {6,S}
""",
    solute=SoluteData(
        S=-0.00778,
        B=-0.02958,
        E=-0.09518,
        L=-0.13622,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 4
B: 4
E: 6
L: 4
A: 5
""",
)

entry(
    index=9,
    label="F-BenzoicAcid",
    group="OR{F-BenzoicAcid(ortho), F-BenzoicAcid(meta), F-BenzoicAcid(para)}",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=10,
    label="F-BenzoicAcid(ortho)",
    group="""
1 * F1s u0 {2,S}
2   Cb  u0 {1,S} {3,B}
3   Cb  u0 {2,B} {4,S}
4   CO  u0 {3,S} {5,S}
5   O2s u0 {4,S} {6,S}
6   H   u0 {5,S}
""",
    solute=SoluteData(
        S=-0.01126,
        B=0.00000,
        E=-0.07247,
        L=-0.04312,
        A=0.07852,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 6
B: 6
E: 6
L: 6
A: 6
""",
)

entry(
    index=11,
    label="F-BenzoicAcid(meta)",
    group="""
1 * F1s u0 {2,S}
2   Cb  u0 {1,S} {3,B}
3   Cb  u0 {2,B} {4,B}
4   Cb  u0 {3,B} {5,S}
5   CO  u0 {4,S} {6,S}
6   O2s u0 {5,S} {7,S}
7   H   u0 {6,S}
""",
    solute=SoluteData(
        S=-0.05665,
        B=0.00000,
        E=-0.07247,
        L=-0.11925,
        A=0.05052,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 6
B: 6
E: 6
L: 6
A: 6
""",
)

entry(
    index=12,
    label="F-BenzoicAcid(para)",
    group="""
1 * F1s u0 {2,S}
2   Cb  u0 {1,S} {3,B}
3   Cb  u0 {2,B} {4,B}
4   Cb  u0 {3,B} {5,B}
5   Cb  u0 {4,B} {6,S}
6   CO  u0 {5,S} {7,S}
7   O2s u0 {6,S} {8,S}
8   H   u0 {7,S}
""",
    solute=SoluteData(
        S=-0.04977,
        B=0.01219,
        E=-0.04977,
        L=-0.11742,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 5
B: 5
E: 5
L: 5
A: 5
""",
)

entry(
    index=13,
    label="F-Aniline",
    group="OR{F-Aniline(ortho), F-Aniline(meta), F-Aniline(para)}",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=14,
    label="F-Aniline(ortho)",
    group="""
1 * F1s u0 {2,S}
2   Cb  u0 {1,S} {3,B}
3   Cb  u0 {2,B} {4,S}
4   N3s u0 {3,S} {5,S} {6,S}
5   H   u0 {4,S}
6   H   u0 {4,S}
""",
    solute=SoluteData(
        S=-0.00377,
        B=-0.00166,
        E=-0.18386,
        L=-0.20968,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 4
B: 2
E: 6
L: 4
A: 4
""",
)

entry(
    index=15,
    label="F-Aniline(meta)",
    group="""
1 * F1s u0 {2,S}
2   Cb  u0 {1,S} {3,B}
3   Cb  u0 {2,B} {4,B}
4   Cb  u0 {3,B} {5,S}
5   N3s u0 {4,S} {6,S} {7,S}
6   H   u0 {5,S}
7   H   u0 {5,S}
""",
    solute=SoluteData(
        S=0.03842,
        B=-0.00651,
        E=-0.17201,
        L=-0.14335,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 3
B: 1
E: 4
L: 3
A: 3
""",
)

entry(
    index=16,
    label="F-Aniline(para)",
    group="""
1 * F1s u0 {2,S}
2   Cb  u0 {1,S} {3,B}
3   Cb  u0 {2,B} {4,B}
4   Cb  u0 {3,B} {5,B}
5   Cb  u0 {4,B} {6,S}
6   N3s u0 {5,S} {7,S} {8,S}
7   H   u0 {6,S}
8   H   u0 {6,S}
""",
    solute=SoluteData(
        S=0.02190,
        B=-0.00166,
        E=-0.14880,
        L=-0.13381,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 4
B: 2
E: 6
L: 4
A: 4
""",
)

tree(
    """
L1: X
	L2: F
		L3: F-C
			L4: F-Cb
				L5: F-Phenol
					L6: F-Phenol(ortho)
					L6: F-Phenol(meta)
					L6: F-Phenol(para)
				L5: F-BenzoicAcid
					L6: F-BenzoicAcid(ortho)
					L6: F-BenzoicAcid(meta)
					L6: F-BenzoicAcid(para)
				L5: F-Aniline
					L6: F-Aniline(ortho)
					L6: F-Aniline(meta)
					L6: F-Aniline(para)
"""
)
