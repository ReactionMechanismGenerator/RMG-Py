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

template(reactants=["Adsorbate", "Proton"], products=["Reduced"], ownReverse=False)

reverse = "Surface_Proton_Electron_Oxidation_Alpha"

reactantNum = 2
productNum = 1
allowChargedSpecies = True
electrons = -1

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
1 *1 R!H u0 {2,[D,T,Q]}
2 *2 X u0 {1,[D,T,Q]}
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
    index = 4,
    label = "CX",
    group =
"""
1 *1 C u0 {2,[D,T,Q]}
2 *2 X u0 {1,[D,T,Q]}
""",
    kinetics = None,
)

entry(
    index = 5,
    label = "CTX",
    group =
"""
1 *1 C u0 {2,T}
2 *2 X u0 {1,T}
""",
    kinetics = None,
)

entry(
    index = 6,
    label = "HCX",
    group =
"""
1 *1 C u0 {2,T} {3,S}
2 *2 X u0 {1,T}
3    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 7,
    label = "C=X",
    group =
"""
1 *1 C u0 {2,D}
2 *2 X u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 8,
    label = "H2C=X",
    group =
"""
1 *1 C u0 {2,D} {3,S} {4,S}
2 *2 X u0 {1,D}
3    H u0 {1,S}
4    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 9,
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
    index = 10,
    label = "OX",
    group =
"""
1 *1 O u0 {2,D}
2 *2 X u0 {1,D}
""",
    kinetics = None,
)


entry(
    index = 11,
    label = "NX",
    group =
"""
1 *1 N u0 {2,[D,T]}
2 *2 X u0 {1,[D,T]}
""",
    kinetics = None,
)

entry(
    index = 12,
    label = "NTX",
    group =
"""
1 *1 N u0 {2,T}
2 *2 X u0 {1,T}
""",
    kinetics = None,
)

entry(
    index = 13,
    label = "N=X",
    group =
"""
1 *1 N u0 {2,D}
2 *2 X u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 14,
    label = "HN=X",
    group =
"""
1 *1 N u0 {2,D} {3,S}
2 *2 X u0 {1,D}
3    N u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 15,
    label = "N-N=X",
    group =
"""
1 *1 N u0 {2,D} {3,S}
2 *2 X u0 {1,D}
3    N u0 {1,S}
""",
    kinetics = None,
)

tree(
"""
L1: Adsorbate
    L2: CX
        L3: CTX
            L4: HCX
        L3: C=X
            L4: O=C=X
            L4: H2C=X
    L2: OX
    L2: NX
        L3: NTX
        L3: N=X
            L4: HN=X
            L4: N-N=X

L1: Proton
"""
)
