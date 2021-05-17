#!/usr/bin/env python
# encoding: utf-8

name = "Surface_Abstraction_vdW/groups"
shortDesc = u""
longDesc = u"""
A vdW species splitting, adsorbing to the surface, and transferring a functional group to a double, triple, or
quadruple bonded surface species.

 *2-*3  *4             *2     *4-*3
  :     ||     ---->    |      |
~*1~ + ~*5~~          ~*1~ + ~*5~~

The rate, which should be in mol/m2/s,
will be given by k * (mol/m2) * (mol/m2)
so k should be in (m2/mol/s).

"""

template(reactants=["AdsorbateVdW","Adsorbate1"], products=["Adsorbate2","Adsorbate3"], ownReverse=False)

reverse = "Surface_Reverse_Abstraction_vdW"  # maybe needs a better name?

reactantNum=2
productNum=2

recipe(actions=[
    ['BREAK_BOND', '*2', 1, '*3'],
    ['CHANGE_BOND', '*1', 1, '*2'],
    ['CHANGE_BOND', '*4', -1, '*5'],
    ['FORM_BOND', '*3', 1, '*4'],
])

entry(
    index = 1,
    label = "AdsorbateVdW",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 R  ux px cx {3,S}
3 *3 R  ux px cx {2,S}
""",
    kinetics = None,
)

entry(
    index = 2,
    label = "Adsorbate1",
    group =
"""
1 *5 X   u0 p0 c0 {2,[D,T,Q]}
2 *4 R!H ux px cx {1,[D,T,Q]}
""",
    kinetics = None,
)

entry(
    index = 3,
    label = "O-R",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 O  ux p2 cx {3,S}
3 *3 R  ux px cx {2,S}
""",
    kinetics = None,
)

entry(
    index = 4,
    label = "C-R",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 C  ux px cx {3,S}
3 *3 R  ux px cx {2,S}
""",
    kinetics = None,
)

entry(
    index = 5,
    label = "*=O",
    group =
"""
1 *5 X u0 p0 c0 {2,D}
2 *4 O u0 p2 c0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 6,
    label = "*C",
    group =
"""
1 *5 X u0 p0 c0 {2,[D,T,Q]}
2 *4 C ux px cx {1,[D,T,Q]}
""",
    kinetics = None,
)

entry(
    index = 7,
    label = "H-H",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 H  u0 p0 c0 {3,S}
3 *3 H  u0 p0 c0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 8,
    label = "N-R",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 N  ux px cx {3,S}
3 *3 R  ux px cx {2,S}
""",
    kinetics = None,
)

entry(
    index = 9,
    label = "*N",
    group =
"""
1 *5 X u0 p0 c0 {2,[D,T]}
2 *4 N ux px cx {1,[D,T]}
""",
    kinetics = None,
)

entry(
    index = 10,
    label = "*=C",
    group =
"""
1 *5 X u0 p0 c0 {2,D}
2 *4 C ux px cx {1,D}
""",
    kinetics = None,
)

entry(
    index = 11,
    label = "*#C",
    group =
"""
1 *5 X u0 p0 c0 {2,T}
2 *4 C u0 p0 c0 {1,T}
""",
    kinetics = None,
)

entry(
    index = 12,
    label = "*$C",
    group =
"""
1 *5 X u0 p0 c0 {2,Q}
2 *4 C u0 p0 c0 {1,Q}
""",
    kinetics = None,
)

entry(
    index = 13,
    label = "*#CH",
    group =
"""
1 *5 X u0 p0 c0 {2,T}
2 *4 C u0 p0 c0 {1,T} {3,S}
3    H u0 p0 c0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 14,
    label = "*=C=R",
    group =
"""
1 *5 X   u0 p0 c0 {2,D}
2 *4 C   u0 p0 c0 {1,D} {3,D}
3    R!H ux px cx {2,D}
""",
    kinetics = None,
)

entry(
    index = 15,
    label = "*=C-2R",
    group =
"""
1 *5 X u0 p0 c0 {2,D}
2 *4 C u0 p0 c0 {1,D} {3,S} {4,S}
3    R ux px cx {2,S}
4    R ux px cx {2,S}
""",
    kinetics = None,
)

entry(
    index = 16,
    label = "*=N",
    group =
"""
1 *5 X u0 p0 c0 {2,D}
2 *4 N u0 p1 c0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 17,
    label = "*#N",
    group =
"""
1 *5 X u0 p0 c0 {2,T}
2 *4 N u0 p1 c0 {1,T}
""",
    kinetics = None,
)

entry(
    index = 18,
    label = "*=NH",
    group =
"""
1 *5 X u0 p0 c0 {2,D}
2 *4 N u0 p1 c0 {1,D} {3,S}
3    H u0 p0 c0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 19,
    label = "*=NO",
    group =
"""
1 *5 X u0 p0 c0 {2,D}
2 *4 N u0 p1 c0 {1,D} {3,S}
3    O u0 p2 c0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 20,
    label = "*=NOH",
    group =
"""
1 *5 X u0 p0 c0 {2,D}
2 *4 N u0 p1 c0 {1,D} {3,S}
3    O u0 p2 c0 {2,S} {4,S}
4    H u0 p0 c0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 21,
    label = "O-O",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 O  u0 p2 c0 {3,S}
3 *3 O  u0 p2 c0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 22,
    label = "HO-OH",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 O  u0 p2 c0 {3,S} {5,S}
3 *3 O  u0 p2 c0 {2,S} {4,S}
4    H  u0 p0 c0 {3,S}
5    H  u0 p0 c0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 23,
    label = "O-N",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 O  u0 p2 c0 {3,S}
3 *3 N  u0 p1 c0 {2,S} {4,[S,D]}
4    R  ux px cx {3,[S,D]}
""",
    kinetics = None,
)

entry(
    index = 24,
    label = "O-C",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 O  u0 p2 c0 {3,S}
3 *3 C  u0 p0 c0 {2,S} {4,[S,D,T]}
4    R  ux px cx {3,[S,D,T]}
""",
    kinetics = None,
)

entry(
    index = 25,
    label = "O-C-3R",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 O  u0 p2 c0 {3,S}
3 *3 C  u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
4    R  ux px cx {3,S}
5    R  ux px cx {3,S}
6    R  ux px cx {3,S}
""",
    kinetics = None,
)

entry(
    index = 26,
    label = "O-C=R",
    group =
"""
multiplicity [1]
1 *1 Xv  u0 p0 c0
2 *2 O   u0 p2 c0 {3,S}
3 *3 C   u0 p0 c0 {2,S} {4,D} {5,S}
4    R!H ux px cx {3,D}
5    R   ux px cx {3,S}
""",
    kinetics = None,
)

entry(
    index = 27,
    label = "O-C#R",
    group =
"""
multiplicity [1]
1 *1 Xv  u0 p0 c0
2 *2 O   u0 p2 c0 {3,S}
3 *3 C   u0 p0 c0 {2,S} {4,T}
4    R!H ux px cx {3,T}
""",
    kinetics = None,
)

entry(
    index = 28,
    label = "O-N-2R",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 O  u0 p2 c0 {3,S}
3 *3 N  u0 p1 c0 {2,S} {4,S} {5,S}
4    R  ux px cx {3,S}
5    R  ux px cx {3,S}
""",
    kinetics = None,
)

entry(
    index = 29,
    label = "O-NHH",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 O  u0 p2 c0 {3,S}
3 *3 N  u0 p1 c0 {2,S} {4,S} {5,S}
4    H  u0 p0 c0 {3,S}
5    H  u0 p0 c0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 30,
    label = "O-N=R",
    group =
"""
multiplicity [1]
1 *1 Xv  u0 p0 c0
2 *2 O   u0 p2 c0 {3,S}
3 *3 N   u0 p1 c0 {2,S} {4,D}
4    R!H ux px cx {3,D}
""",
    kinetics = None,
)

entry(
    index = 31,
    label = "C-C",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 C  ux px cx {3,S}
3 *3 C  ux px cx {2,S}
""",
    kinetics = None,
)

entry(
    index = 32,
    label = "C-O",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 C  ux px cx {3,S}
3 *3 O  ux p2 cx {2,S}
""",
    kinetics = None,
)

entry(
    index = 33,
    label = "C-N",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 C  ux px cx {3,S}
3 *3 N  u0 p1 c0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 34,
    label = "C-OH",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 C  ux px cx {3,S}
3 *3 O  u0 p2 c0 {2,S} {4,S}
4    H  u0 p0 c0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 35,
    label = "N-N",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 N  ux px cx {3,S}
3 *3 N  ux px cx {2,S}
""",
    kinetics = None,
)

entry(
    index = 36,
    label = "N-O",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 N  ux px cx {3,S}
3 *3 O  u0 p2 c0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 37,
    label = "N-C",
    group =
"""
multiplicity [1]
1 *1 Xv u0 p0 c0
2 *2 N  ux px cx {3,S}
3 *3 C  ux px cx {2,S}
""",
    kinetics = None,
)

tree(
"""
L1: AdsorbateVdW
    L2: H-H
    L2: O-R
        L3: O-O
            L4: HO-OH
        L3: O-N
            L4: O-N-2R
                L5: O-NHH
            L4: O-N=R
        L3: O-C
            L4: O-C-3R
            L4: O-C=R
            L4: O-C#R
    L2: C-R
        L3: C-C
        L3: C-O
            L4: C-OH
        L3: C-N
    L2: N-R
        L3: N-N
        L3: N-O
        L3: N-C

L1: Adsorbate1
    L2: *=O
    L2: *C
        L3: *=C
            L4: *=C=R
            L4: *=C-2R
        L3: *#C
            L4:*#CH
        L3: *$C
    L2: *N
        L3: *=N
            L4: *=NH
            L4: *=NO
                L5: *=NOH
        L3: *#N
"""
)
