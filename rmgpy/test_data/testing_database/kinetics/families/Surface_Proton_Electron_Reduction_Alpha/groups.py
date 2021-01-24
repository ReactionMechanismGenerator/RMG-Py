#!/usr/bin/env python
# encoding: utf-8

name = "Surface_Proton_Electron_Reduction_Alpha/groups"
shortDesc = u""
longDesc = u"""

   *1                        *1-*3H
   ||  + *3H+ + *e-  ---->   |
  ~*2~                      ~*2~~

The rate, which should be in mol/m2/s,
will be given by k * (mol/m2) * (mol/m3) * 1
so k should be in (m3/mol/s).
"""

template(reactants=["Adsorbate", "Proton", "Electron"], products=["Reduced"], ownReverse=False)

reverse = "Surface_Proton_Electron_Oxidation_Alpha"

reactantNum = 3
productNum = 1

recipe(actions=[
    ['LOSE_CHARGE', '*3', 1],
    ['CHANGE_BOND', '*1', -1, '*2'],
    ['FORM_BOND', '*1', 1, '*3'],
])

entry(
    index = 1,
    label = "Adsorbate",
    group =
"""
1 *1 R u0 {2,[S,D,T,Q]}
2 *2 X u0 {1,[S,D,T,Q]}
""",
    kinetics = None,
)

entry(
    index = 2,
    label = "Proton",
    group =
"""
1 *3 H+ u0 p0 c+1
""",
    kinetics = None,
)

entry(
    index = 3,
    label = "Electron",
    group =
"""
1 * e u0 p0 c-1
""",
    kinetics = None,
)

entry(
    index = 4,
    label = "CX",
    group =
"""
1 *1 C u0 {2,[S,D,T,Q]}
2 *2 X u0 {1,[S,D,T,Q]}
""",
    kinetics = None,
)

entry(
    index = 5,
    label = "OX",
    group =
"""
1 *1 O u0 {2,[S,D]}
2 *2 X u0 {1,[S,D]}
""",
    kinetics = None,
)


entry(
    index = 6,
    label = "NX",
    group =
"""
1 *1 N u0 {2,[S,D,T]}
2 *2 X u0 {1,[S,D,T]}
""",
    kinetics = None,
)

entry(
    index = 7,
    label = "O=C=X",
    group =
"""
1 *1 C u0 {2,D} {3,D}
2 *2 X u0 {1,D}
3    O2d u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 8,
    label = "N-N=X",
    group =
"""
1 *1 N u0 {2,D} {3,S}
2 *2 X u0 {1,D}
3    N u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 9,
    label = "N=N-X",
    group =
"""
1 *1 N u0 {2,S} {3,D}
2 *2 X u0 {1,S}
3    N u0 {1,D}
""",
    kinetics = None,
)


tree(
"""
L1: Adsorbate
    L2: CX
        L3: O=C=X
    L2: OX
    L2: NX
        L3: N=N-X
        L3: N-N=X

L1: Proton

L1: Electron
"""
)
