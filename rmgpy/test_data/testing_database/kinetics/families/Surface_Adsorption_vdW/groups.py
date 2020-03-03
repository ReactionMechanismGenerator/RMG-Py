#!/usr/bin/env python
# encoding: utf-8

name = "Surface_Adsorption_vdW/groups"
shortDesc = u""
longDesc = u"""
Physisorption of a gas-phase species onto the surface.

 *1         *1
     ---->   :
~*2~       ~*2~~

The rate, which should be in mol/m2/s,
will be given by k * (mol/m2) * (mol/m3)
so k should be in (m3/mol/s). We will use sticking coefficients.
"""

template(reactants=["Adsorbate", "VacantSite"], products=["Adsorbed"], ownReverse=False)

reverse = "Surface_Desorption_vdW"

reactantNum=2
productNum=1

recipe(actions=[
    ['FORM_BOND', '*1', 0, '*2']
])

entry(
    index = 1,
    label = "Adsorbate",
    group =
"""
multiplicity [1,3]
1 *1 R u0 px c0
""",
    kinetics = None,
)

entry(
    index = 2,
    label="VacantSite",
    group =
"""
1 *2 Xv u0 p0 c0
""",
    kinetics = None,
)



tree(
"""
L1: Adsorbate

L1: VacantSite
"""
)

forbidden(
    label = "radical1",
    group =
"""
1 R u[1,2,3]
""",
    shortDesc = u"""Radicals not allowed""",
    longDesc =
u"""
""",
)


forbidden(
    label = "charge1",
    group =
"""
1 R ux c+1
""",
    shortDesc = u"""Charges not allowed""",
    longDesc =
u"""
""",
)

forbidden(
    label = "charge2",
    group =
"""
1 R ux c-1
""",
    shortDesc = u"""Charges not allowed""",
    longDesc =
u"""
""",
)
