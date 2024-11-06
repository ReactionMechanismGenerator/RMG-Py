name = "Surface Adsorption Corrections"
shortDesc = ""
longDesc = """
Changes due to adsorbing on a surface.
Here, Pt(111)
Note: "-h" means "horizontal".
"""

entry(
    index = 1,
    label = "R*",
    group=
"""
1 R u0
2 X u0
""",
    thermo=None,
    shortDesc="""Anything adsorbed anyhow.""",
    longDesc="""
   R
   x
***********
This node should be empty, ensuring that one of the nodes below is used.
""",
    metal = "Pt",
    facet = "111",
)

entry(
    index = 1,
    label = "R-*",
    group =
"""
1 X  u0 p0 c0 {2,S}
2 R  u0 p0 c0 {1,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-3.01, -1.78, -0.96, -0.41, 0.23, 0.56, 0.91], 'cal/(mol*K)'),
        H298=(-86.29, 'kcal/mol'),
        S298=(-26.39, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from H single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   R
   |
***********
""",
    metal = "Pt",
    facet = "111",
)
