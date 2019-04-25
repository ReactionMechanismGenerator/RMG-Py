#!/usr/bin/env python
# encoding: utf-8

name = "Surface_Dissociation_vdW/groups"
shortDesc = u""
longDesc = u"""
Surface bond fission of one vdW species into two distinct adsorbates. Atom *3 is vdW, but we leave that out of the graph.

 *1--*2                 *1     *2
  :            ---->     |      |
~*3~ + ~*4~~           ~*3~ + ~*4~~


The rate, which should be in mol/m2/s,
will be given by k * (mol/m2) * (mol/m2)
so k should be in (m2/mol/s)
"""

template(reactants=["Combined", "VacantSite"], products=["Adsorbate1", "Adsorbate2"], ownReverse=False)

reverse = "Surface_Association_vdW"

recipe(actions=[
    ['FORM_BOND', '*1', 1, '*3'],
    ['FORM_BOND', '*2', 1, '*4'],
    ['BREAK_BOND', '*1', 1, '*2'],
])

entry(
    index = 1,
    label = "Combined",
    group = 
"""
1 *1 R  u0 {2,S} 
2 *2 R  u0 {1,S}
3 *3 X  u0 
""",
    kinetics = None,
)

entry(
    index = 2,
    label="VacantSite",
    group = 
"""
1 *4 Xv u0
""",
    kinetics = None,
)


tree(
"""
L1: Combined

L1: VacantSite
"""
)

