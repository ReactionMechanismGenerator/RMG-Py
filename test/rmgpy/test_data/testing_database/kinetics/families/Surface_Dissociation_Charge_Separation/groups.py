#!/usr/bin/env python
# encoding: utf-8

name = "Surface_Dissociation_Charge_Separation/groups"
shortDesc = u""
longDesc = u"""
Surface bond fission of one species into two distinct adsorbates.
Atom *1 is bonded to the surface (*3). The image below shows a single bond,
but single and double are possible. What matters is that the bond
between *1 and *2 must be single and have charge separation across it.

   *1[+]--*2[-]          *1      *2
    |            ---->    |      ||
  ~*3~ + ~*4~~           ~*3~ + ~*4~~

The rate, which should be in mol/m2/s,
will be given by k * (mol/m2) * (mol/m2)
so k should be in (m2/mol/s)
"""

template(reactants=["Combined", "VacantSite"], products=["Adsorbate1", "Adsorbate2"], ownReverse=False)

reverse = "Surface_Association_Charge_Separation"

reactantNum=2
productNum=2

recipe(actions=[
    ['BREAK_BOND', '*1', 1, '*2'],
    ['FORM_BOND', '*2', 1, '*4'],
    ['CHANGE_BOND', '*2', 1, '*4'],
    ['LOSE_PAIR','*2','1'],
    ['GAIN_PAIR', '*1', '1'],
    ['LOSE_CHARGE','*1','1'],
    ['GAIN_CHARGE', '*2', '1'],
])

entry(
    index = 1,
    label = "Combined",
    group =
"""
1 *1 Val5 u0 p0 c+1 {2,S} {3,[S,D]}
2 *2 R!H u0 p[1,2,3] c-1 {1,S}
3 *3 Xo u0 p0 c0 {1,[S,D]}
""",
    kinetics = None,
)

entry(
    index = 2,
    label="VacantSite",
    group =
"""
1 *4 Xv u0 p0 c0
""",
    kinetics = None,
)

tree(
"""
L1: Combined
L1: VacantSite
"""
)

forbidden(
    label = "Surf",
    group =
"""
1 *1 Val5 u0 p0 c+1 {2,S} {3,[S,D]}
2 *2 R!H u0 p[1,2,3] c-1 {1,S} {4,[S,D,T]}
3 *3 Xo u0 p0 c0 {1,[S,D]}
4 Xo u0 c0 {2,[S,D,T]}
""",
)
