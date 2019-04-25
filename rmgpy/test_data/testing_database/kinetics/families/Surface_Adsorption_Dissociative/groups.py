#!/usr/bin/env python
# encoding: utf-8

name = "Surface_Adsorption_Dissociative/groups"
shortDesc = u""
longDesc = u"""
Dissociative adsorption of a gas-phase species onto the surface. The single-bond in the gas-phase species is split; the resulting fragments each are singled bonded to the surface.

 *1-*2               *1      *2
             ---->    |       |
~*3~ + ~*4~         ~*3~~ + ~*4~~

The rate, which should be in mol/m2/s,
will be given by k * (mol/m2) * (mol/m2) * (mol/m3)
so k should be in (m5/mol2/s). We will use sticking coefficients.
"""

template(reactants=["Adsorbate", "VacantSite1", "VacantSite2"], products=["Adsorbed1", "Adsorbed2"], ownReverse=False)

reverse = "Surface_Desorption_Associative"

recipe(actions=[
    ['BREAK_BOND', '*1', 1, '*2'],
    ['FORM_BOND', '*1', 1, '*3'],
    ['FORM_BOND', '*2', 1, '*4']
])

entry(
    index = 1,
    label = "Adsorbate",
    group =
"""
1 *1 R u0 {2,S}
2 *2 R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 2,
    label="VacantSite1",
    group =
"""
1 *3 Xv u0
""",
    kinetics = None,
)

entry(
    index = 3,
    label="VacantSite2",
    group =
"""
1 *4 Xv u0
""",
    kinetics = None,
)

tree(
"""
L1: Adsorbate

L1: VacantSite1

L1: VacantSite2
"""
)


forbidden(
    label = "adjacentradical1",
    group =
"""
1 *1 R u0 {2,[S,D,T]}
2    R u1 {1,[S,D,T]}
""",
    shortDesc = u"""""",
    longDesc =
u"""
The adsorbing atom should not be adjacent to a radical.
e.g. this is not allowed:

CH2.-CH3    -->   CH2.-CH2   +   H
                       |         |
     X X               X         X
""",
)

forbidden(
    label = "adjacentradical2",
    group =
"""
1 *2 R u0 {2,[S,D,T]}
2    R u1 {1,[S,D,T]}
""",
    shortDesc = u"""""",
    longDesc =
u"""
Neither adsorbing atom should be adjacent to a radical
e.g. this is not allowed:

CH2.-CH3    -->   CH2.-CH2   +   H
                       |         |
     X X               X         X
""",
)


forbidden(
    label = "disigma1",
    group =
"""
1 *1 R u0 {2,[S,D,T]}
2    R u0 {1,[S,D,T]} {3,[S,D,T]}
3    X u0 {2,[S,D,T]}
""",
    shortDesc = u"""""",
    longDesc =
u"""
The adsorbing atom should not be adjacent to an atom that is already adsorbed.
e.g. this is not allowed:

H  O=C-O   <-->  O=CH-O
|    | |              |
X    X X         X X  X
""",
)

forbidden(
    label = "disigma2",
    group =
"""
1 *2 R u0 {2,[S,D,T]}
2    R u0 {1,[S,D,T]} {3,[S,D,T]}
3    X u0 {2,[S,D,T]}
""",
    shortDesc = u"""""",
    longDesc =
u"""
The adsorbing atom should not be adjacent to an atom that is already adsorbed.
e.g. this is not allowed:

H  O=C-O   <-->  O=CH-O
|    | |              |
X    X X         X X  X
""",
)

forbidden(
    label = "disigma3",
    group =
"""
1 *1 R u0 {2,[S,D,T]}
2    R u0 {1,[S,D,T]} {3,[S,D,T]}
3    R u0 {2,[S,D,T]} {4,[S,D,T]}
4    X u0 {3,[S,D,T]}
""",
    shortDesc = u"""""",
    longDesc =
u"""
The adsorbing atom should not be next-nearest neighbor to an atom that is already adsorbed.
""",
)


forbidden(
    label = "disigma4",
    group =
"""
1 *2 R u0 {2,[S,D,T]}
2    R u0 {1,[S,D,T]} {3,[S,D,T]}
3    R u0 {2,[S,D,T]} {4,[S,D,T]}
4    X u0 {3,[S,D,T]}
""",
    shortDesc = u"""""",
    longDesc =
u"""
The adsorbing atom should not be next-nearest neighbor to an atom that is already adsorbed.
""",
)
