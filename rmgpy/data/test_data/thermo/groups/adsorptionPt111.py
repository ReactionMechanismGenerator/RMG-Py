#!/usr/bin/env python
# encoding: utf-8

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
    metal="Pt",
    facet="111",
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
    metal="Pt",
    facet="111",
)

entry(
    index = 2,
    label = "(R2)*",
    group =
"""
1 X  u0 p0 c0
2 R  u0 p0 c0 {3,S}
3 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.13, 1.17, 1.19, 1.2, 1.21, 1.21, 1.22], 'cal/(mol*K)'),
        H298=(-1.22, 'kcal/mol'),
        S298=(-7.73, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from H2 vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

  R-R
   :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 3,
    label = "(OR2)*",
    group =
"""
1 X  u0 p0 c0
2 O  u0 p2 c0 {3,S} {4,S}
3 R  u0 p0 c0 {2,S}
4 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.88, 1.49, 1.82, 2.02, 2.22, 2.33, 2.43], 'cal/(mol*K)'),
        H298=(-4.85, 'kcal/mol'),
        S298=(-22.53, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from H2O vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 RO-R
   :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 4,
    label = "O-*R",
    group =
"""
1 X  u0 p0 c0 {2,S}
2 O  u0 p2 c0 {1,S} {3,S}
3 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.34, 2.23, 2.7, 2.97, 3.25, 3.38, 3.5], 'cal/(mol*K)'),
        H298=(-46.18, 'kcal/mol'),
        S298=(-33.89, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from OH single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   R
   |
   O
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 5,
    label = "(OROR)*",
    group =
"""
1 X  u0 p0 c0
2 O  u0 p2 c0 {3,S} {4,S}
3 O  u0 p2 c0 {2,S} {5,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.85, 2.13, 2.27, 2.35, 2.45, 2.51, 2.57], 'cal/(mol*K)'),
        H298=(-6.72, 'kcal/mol'),
        S298=(-26.31, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HO-OH vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 RO-OR
   :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 6,
    label = "O-*O-*",
    group =
"""
1 X  u0 p0 c0 {2,S} {3,S}
2 X  u0 p0 c0 {1,S} {4,S}
3 O  u0 p2 c0 {1,S} {4,S}
4 O  u0 p2 c0 {2,S} {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([2.43, 3.46, 3.9, 4.05, 4.07, 4.0, 3.85], 'cal/(mol*K)'),
        H298=(-8.59, 'kcal/mol'),
        S298=(-40.49, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from O2 bidentate, twice single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   O--O
   |  |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 7,
    label = "O-*OR",
    group =
"""
1 X  u0 p0 c0 {2,S}
2 O  u0 p2 c0 {1,S} {3,S}
3 O  u0 p2 c0 {2,S} {4,S}
4 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([2.69, 3.13, 3.18, 3.12, 2.92, 2.76, 2.56], 'cal/(mol*K)'),
        H298=(-17.47, 'kcal/mol'),
        S298=(-31.56, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from OOH single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   OR
   |
   O
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 8,
    label = "O=*",
    group =
"""
1 X  u0 p0 c0 {2,D}
2 O  u0 p2 c0 {1,D}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-0.39, 0.26, 0.58, 0.77, 0.96, 1.05, 1.14], 'cal/(mol*K)'),
        H298=(-99.97, 'kcal/mol'),
        S298=(-30.95, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from O double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   O
   ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 9,
    label = "O-*NR2",
    group =
"""
1 X  u0 p0 c0 {3,S}
2 N  u0 p1 c0 {3,S} {4,S} {5,S}
3 O  u0 p2 c0 {1,S} {2,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.33, 1.18, 1.65, 1.93, 2.2, 2.32, 2.41], 'cal/(mol*K)'),
        H298=(-16.75, 'kcal/mol'),
        S298=(-33.37, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from O-NH2 single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   NR2
   |
   O
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 10,
    label = "O-*CR3",
    group =
"""
1 X  u0 p0 c0 {3,S}
2 C  u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 O  u0 p2 c0 {1,S} {2,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.28, 0.58, 0.85, 1.08, 1.4, 1.61, 1.93], 'cal/(mol*K)'),
        H298=(-32.28, 'kcal/mol'),
        S298=(-34.6, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from O-CH3 single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   CR3
   |
   O
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 11,
    label = "(NR3)*",
    group =
"""
1 X  u0 p0 c0
2 N  u0 p1 c0 {3,S} {4,S} {5,S}
3 R  u0 p0 c0 {2,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.23, 2.36, 3.08, 3.56, 4.11, 4.4, 4.69], 'cal/(mol*K)'),
        H298=(-16.11, 'kcal/mol'),
        S298=(-32.0, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from NH3 vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R2N-R
    :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 12,
    label = "N-*R2",
    group =
"""
1 X  u0 p0 c0 {2,S}
2 N  u0 p1 c0 {1,S} {3,S} {4,S}
3 R  u0 p0 c0 {2,S}
4 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-1.05, 0.88, 2.07, 2.81, 3.6, 3.99, 4.4], 'cal/(mol*K)'),
        H298=(-48.33, 'kcal/mol'),
        S298=(-47.88, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from NH2 single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   NR2
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 13,
    label = "N=*R",
    group =
"""
1 X  u0 p0 c0 {2,D}
2 N  u0 p1 c0 {1,D} {3,S}
3 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-2.14, -0.29, 0.86, 1.58, 2.37, 2.76, 3.18], 'cal/(mol*K)'),
        H298=(-80.92, 'kcal/mol'),
        S298=(-40.72, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from NH double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 14,
    label = "N#*",
    group =
"""
1 X  u0 p0 c0 {2,T}
2 N  u0 p1 c0 {1,T}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-1.14, -0.25, 0.24, 0.52, 0.81, 0.96, 1.1], 'cal/(mol*K)'),
        H298=(-147.51, 'kcal/mol'),
        S298=(-32.92, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from N triple-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

    N
   |||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 15,
    label = "(NR2OR)*",
    group =
"""
1 X  u0 p0 c0
2 N  u0 p1 c0 {3,S} {4,S} {5,S}
3 O  u0 p2 c0 {2,S} {6,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-0.44, 0.2, 0.73, 1.14, 1.69, 2.01, 2.35], 'cal/(mol*K)'),
        H298=(-15.69, 'kcal/mol'),
        S298=(-32.2, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from H2N-OH vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R2N-OR
    :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 16,
    label = "(NRO)*",
    group =
"""
1 X  u0 p0 c0
2 N  u0 p1 c0 {3,D} {4,S}
3 O  u0 p2 c0 {2,D}
4 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-0.27, 0.8, 1.4, 1.72, 1.98, 2.06, 2.17], 'cal/(mol*K)'),
        H298=(-30.08, 'kcal/mol'),
        S298=(-32.78, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HN-O vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

  RN=O
    :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 17,
    label = "N-*ROR",
    group =
"""
1 X  u0 p0 c0 {2,S}
2 N  u0 p1 c0 {1,S} {3,S} {4,S}
3 O  u0 p2 c0 {2,S} {5,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([2.24, 3.32, 3.89, 4.22, 4.56, 4.73, 4.88], 'cal/(mol*K)'),
        H298=(-32.32, 'kcal/mol'),
        S298=(-45.51, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HN-OH single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R-N-OR
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 18,
    label = "N-*O",
    group =
"""
1 X  u0 p0 c0 {2,S}
2 N  u0 p1 c0 {1,S} {3,D}
3 O  u0 p2 c0 {2,D}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.81, 2.7, 3.19, 3.47, 3.71, 3.77, 3.75], 'cal/(mol*K)'),
        H298=(-37.18, 'kcal/mol'),
        S298=(-40.63, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from NO single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   O
   ||
   N
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 19,
    label = "N=*O-*",
    group =
"""
1 X  u0 p0 c0 {2,S} {3,D}
2 X  u0 p0 c0 {1,S} {4,S}
3 N  u0 p1 c0 {1,D} {4,S}
4 O  u0 p2 c0 {2,S} {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.05, 0.57, 0.86, 1.04, 1.2, 1.25, 1.26], 'cal/(mol*K)'),
        H298=(-32.66, 'kcal/mol'),
        S298=(-29.32, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from NO-h bidentate, double- and single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   N--O
  ||  |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 20,
    label = "N=*OR",
    group =
"""
1 X  u0 p0 c0 {2,D}
2 N  u0 p1 c0 {1,D} {3,S}
3 O  u0 p2 c0 {2,S} {4,S}
4 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([2.65, 3.78, 4.29, 4.49, 4.55, 4.5, 4.48], 'cal/(mol*K)'),
        H298=(-75.72, 'kcal/mol'),
        S298=(-44.7, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from NOH double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   OR
   |
   N
  ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 21,
    label = "(NR2NR2)*",
    group =
"""
1 X  u0 p0 c0
2 N  u0 p1 c0 {3,S} {4,S} {5,S}
3 N  u0 p1 c0 {2,S} {6,S} {7,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {3,S}
7 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.13, 0.74, 1.16, 1.46, 1.84, 2.06, 2.3], 'cal/(mol*K)'),
        H298=(-23.19, 'kcal/mol'),
        S298=(-31.95, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from H2N-NH2 vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R2N-NR2
    :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 22,
    label = "(NRNR)*",
    group =
"""
1 X  u0 p0 c0
2 N  u0 p1 c0 {3,D} {4,S}
3 N  u0 p1 c0 {2,D} {5,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([3.21, 4.62, 5.24, 5.46, 5.43, 5.28, 5.02], 'cal/(mol*K)'),
        H298=(-20.58, 'kcal/mol'),
        S298=(-42.07, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HN-NH vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 RN=NR
   :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 23,
    label = "N-*N-*",
    group =
"""
1 X  u0 p0 c0 {2,S} {3,S}
2 X  u0 p0 c0 {1,S} {4,S}
3 N  u0 p1 c0 {1,S} {4,D}
4 N  u0 p1 c0 {2,S} {3,D}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.2, 1.21, 1.21, 1.21, 1.22, 1.22, 1.22], 'cal/(mol*K)'),
        H298=(-2.39, 'kcal/mol'),
        S298=(-13.89, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from NN bidentate, twice single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

  N==N
  |  |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 24,
    label = "N-*RNR2",
    group =
"""
1 X  u0 p0 c0 {2,S}
2 N  u0 p1 c0 {1,S} {3,S} {4,S}
3 N  u0 p1 c0 {2,S} {5,S} {6,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {3,S}
6 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.92, 2.91, 3.52, 3.91, 4.35, 4.57, 4.79], 'cal/(mol*K)'),
        H298=(-29.97, 'kcal/mol'),
        S298=(-45.43, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HN-NH2 single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R-N-NR2
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 25,
    label = "N-*NR",
    group =
"""
1 X  u0 p0 c0 {2,S}
2 N  u0 p1 c0 {1,S} {3,D}
3 N  u0 p1 c0 {2,D} {4,S}
4 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.74, 2.91, 3.56, 3.93, 4.25, 4.37, 4.52], 'cal/(mol*K)'),
        H298=(-25.14, 'kcal/mol'),
        S298=(-43.45, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from N-NH single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   NR
  ||
   N
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 26,
    label = "N=*NR2",
    group =
"""
1 X  u0 p0 c0 {2,D}
2 N  u0 p1 c0 {1,D} {3,S}
3 N  u0 p1 c0 {2,S} {4,S} {5,S}
4 R  u0 p0 c0 {3,S}
5 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([3.32, 4.56, 5.07, 5.2, 5.03, 4.79, 4.55], 'cal/(mol*K)'),
        H298=(-47.66, 'kcal/mol'),
        S298=(-43.17, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from N-NH2 double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   NR2
   |
   N
  ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 27,
    label = "N-*RN-*R",
    group =
"""
1 X  u0 p0 c0 {2,S} {3,S}
2 X  u0 p0 c0 {1,S} {4,S}
3 N  u0 p1 c0 {1,S} {4,S} {5,S}
4 N  u0 p1 c0 {2,S} {3,S} {6,S}
5 R  u0 p0 c0 {3,S}
6 R  u0 p0 c0 {4,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([2.53, 4.03, 4.78, 5.11, 5.24, 5.17, 5.0], 'cal/(mol*K)'),
        H298=(-23.37, 'kcal/mol'),
        S298=(-43.91, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HN-NH-h bidentate, twice single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 RN--NR
  |  |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 28,
    label = "N-*RCR3",
    group =
"""
1 X  u0 p0 c0 {3,S}
2 C  u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 N  u0 p1 c0 {1,S} {2,S} {7,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {2,S}
7 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.18, 2.22, 2.89, 3.33, 3.85, 4.12, 4.45], 'cal/(mol*K)'),
        H298=(-43.5, 'kcal/mol'),
        S298=(-46.63, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HN-CH3 single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R-N-CR3
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 29,
    label = "N-*CR2",
    group =
"""
1 X  u0 p0 c0 {3,S}
2 C  u0 p0 c0 {3,D} {4,S} {5,S}
3 N  u0 p1 c0 {1,S} {2,D}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.99, 2.96, 3.5, 3.83, 4.17, 4.33, 4.54], 'cal/(mol*K)'),
        H298=(-39.07, 'kcal/mol'),
        S298=(-44.16, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from N-CH2 single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   CR2
  ||
   N
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 30,
    label = "N=*CR3",
    group =
"""
1 X  u0 p0 c0 {3,D}
2 C  u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 N  u0 p1 c0 {1,D} {2,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.12, 1.89, 2.43, 2.81, 3.29, 3.59, 4.07], 'cal/(mol*K)'),
        H298=(-71.1, 'kcal/mol'),
        S298=(-47.17, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from N-CH3 double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   CR3
   |
   N
  ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 31,
    label = "N-*O2",
    group =
"""
1 X  u0  p0 c0 {2,S}
2 N  u0  p0 c+1 {1,S} {3,S} {4,D}
3 O  u0  p2 c-1 {2,S}
4 O  u0  p2 c0 {2,D}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([2.35, 2.6, 2.67, 2.66, 2.61, 2.57, 2.5], 'cal/(mol*K)'),
        H298=(-16.1, 'kcal/mol'),
        S298=(-33.93, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from ON-O single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 O-N=O
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 32,
    label = "Cq*",
    group =
"""
1 X  u0 p0 c0 {2,Q}
2 C  u0 p0 c0 {1,Q}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-2.0, -0.88, -0.22, 0.18, 0.61, 0.82, 1.04], 'cal/(mol*K)'),
        H298=(-156.9, 'kcal/mol'),
        S298=(-31.82, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from C quadruple-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   C
 ||||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 33,
    label = "C-*C-*",
    group =
"""
1 X  u0  p0 c0 {2,S} {3,D}
2 X  u0  p0 c0 {1,S} {4,D}
3 C  u0  p0 c0 {1,D} {4,D}
4 C  u0  p0 c0 {2,D} {3,D}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.98, 2.17, 2.79, 3.13, 3.44, 3.55, 3.63], 'cal/(mol*K)'),
        H298=(-137.31, 'kcal/mol'),
        S298=(-41.99, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from C-C bidentate, twice double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

  C--C
  |  |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 34,
    label = "C=*CR2",
    group =
"""
1 X  u0  p0 c0 {2,D}
2 C  u0  p0 c0 {1,D} {3,D}
3 C  u0  p0 c0 {2,D} {4,S} {5,S}
4 R  u0  p0 c0 {3,S}
5 R  u0  p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-1.2, 0.67, 1.9, 2.71, 3.62, 4.07, 4.52], 'cal/(mol*K)'),
        H298=(-93.15, 'kcal/mol'),
        S298=(-48.06, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from C-CH2 double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   CR2
  ||
   C
  ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 35,
    label = "C#*CR3",
    group =
"""
1 X  u0 p0 c0 {3,T}
2 C  u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C  u0 p0 c0 {1,T} {2,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.98, 1.94, 2.58, 3.04, 3.6, 3.92, 4.33], 'cal/(mol*K)'),
        H298=(-129.74, 'kcal/mol'),
        S298=(-45.92, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from C-CH3 triple-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   CR3
   |
   C
  |||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 36,
    label = "C#*R",
    group =
"""
1 X  u0 p0 c0 {2,T}
2 C  u0 p0 c0 {1,T} {3,S}
3 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-2.35, -0.42, 0.76, 1.49, 2.29, 2.68, 3.14], 'cal/(mol*K)'),
        H298=(-145.5, 'kcal/mol'),
        S298=(-40.0, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from CH triple-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   R
   |
   C
  |||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 37,
    label = "C=*RC=*R",
    group =
"""
1 X  u0 p0 c0 {2,S} {3,D}
2 X  u0 p0 c0 {1,S} {4,D}
3 C  u0 p0 c0 {1,D} {4,S} {5,S}
4 C  u0 p0 c0 {2,D} {3,S} {6,S}
5 R  u0 p0 c0 {3,S}
6 R  u0 p0 c0 {4,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.05, 1.59, 2.7, 3.47, 4.37, 4.8, 5.11], 'cal/(mol*K)'),
        H298=(-47.33, 'kcal/mol'),
        S298=(-31.36, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from CH-CH bidentate, twice double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R-C--C-R
  ||  ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 38,
    label = "C=*R2",
    group =
"""
1 X  u0 p0 c0 {2,D}
2 C  u0 p0 c0 {1,D} {3,S} {4,S}
3 R  u0 p0 c0 {2,S}
4 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-1.54, 0.48, 1.74, 2.53, 3.38, 3.8, 4.29], 'cal/(mol*K)'),
        H298=(-85.5, 'kcal/mol'),
        S298=(-42.7, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from CH2 double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R-C-R
  ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 39,
    label = "C-*R2C-*R2",
    group =
"""
1 X  u0 p0 c0 {2,S} {3,S}
2 X  u0 p0 c0 {1,S} {4,S}
3 C  u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C  u0 p0 c0 {2,S} {3,S} {7,S} {8,S}
5 R  u0 p0 c0 {3,S}
6 R  u0 p0 c0 {3,S}
7 R  u0 p0 c0 {4,S}
8 R  u0 p0 c0 {4,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.94, 3.35, 4.13, 4.56, 4.94, 5.08, 5.11], 'cal/(mol*K)'),
        H298=(-22.63, 'kcal/mol'),
        S298=(-41.46, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from CH2-CH2 bidentate, twice single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R2C--CR2
   |  |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 40,
    label = "C-*R3",
    group =
"""
1 X  u0 p0 c0 {2,S}
2 C  u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 R  u0 p0 c0 {2,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-0.56, 0.74, 1.74, 2.48, 3.45, 4.0, 4.58], 'cal/(mol*K)'),
        H298=(-41.63, 'kcal/mol'),
        S298=(-32.73, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from CH3 single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   CR3
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 41,
    label = "(CR3CR3)*",
    group =
"""
1 X  u0 p0 c0
2 C  u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C  u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {2,S}
7 R  u0 p0 c0 {3,S}
8 R  u0 p0 c0 {3,S}
9 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([2.46, 2.5, 2.52, 2.53, 2.53, 2.53, 2.51], 'cal/(mol*K)'),
        H298=(-4.64, 'kcal/mol'),
        S298=(-15.11, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from CH3-CH3 vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R3C-CR3
    :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 42,
    label = "(CR4)*",
    group =
"""
1 X  u0 p0 c0
2 C  u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 R  u0 p0 c0 {2,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([2.42, 2.45, 2.46, 2.47, 2.48, 2.48, 2.47], 'cal/(mol*K)'),
        H298=(-2.4, 'kcal/mol'),
        S298=(-6.92, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from CH4 vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

  R3C-R
     :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 43,
    label = "C=*N-*",
    group =
"""
1 X  u0  p0 c0 {2,S} {3,D}
2 X  u0  p0 c0 {1,S} {4,S}
3 C  u0  p0 c0 {1,D} {4,D}
4 N  u0  p1 c0 {2,S} {3,D}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([2.99, 3.32, 3.51, 3.63, 3.74, 3.76, 3.74], 'cal/(mol*K)'),
        H298=(-77.01, 'kcal/mol'),
        S298=(-34.98, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from CN bidentate, double- and single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

  C==N
 ||  |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 44,
    label = "C=*NR",
    group =
"""
1 X  u0  p0 c0 {2,D}
2 C  u0  p0 c0 {1,D} {3,D}
3 N  u0  p1 c0 {2,D} {4,S}
4 R  u0  p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([2.64, 3.53, 4.09, 4.44, 4.81, 4.96, 5.04], 'cal/(mol*K)'),
        H298=(-40.56, 'kcal/mol'),
        S298=(-30.68, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from CNH double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

    NR
   ||
    C
   ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 45,
    label = "C#*NR2",
    group =
"""
1 X  u0 p0 c0 {2,T}
2 C  u0 p0 c0 {1,T} {3,S}
3 N  u0 p1 c0 {2,S} {4,S} {5,S}
4 R  u0 p0 c0 {3,S}
5 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([3.38, 4.13, 4.45, 4.59, 4.65, 4.62, 4.6], 'cal/(mol*K)'),
        H298=(-94.24, 'kcal/mol'),
        S298=(-49.82, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from CNH2 triple-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   NR2
   |
   C
  |||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 46,
    label = "C=*O",
    group =
"""
1 X  u0  p0 c0 {2,D}
2 C  u0  p0 c0 {1,D} {3,D}
3 O  u0  p2 c0 {2,D}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([2.21, 2.9, 3.29, 3.53, 3.74, 3.8, 3.78], 'cal/(mol*K)'),
        H298=(-34.7, 'kcal/mol'),
        S298=(-38.09, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from CO-f double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   O
  ||
   C
  ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 47,
    label = "C#*OR",
    group =
"""
1 X  u0 p0 c0 {2,T}
2 C  u0 p0 c0 {1,T} {3,S}
3 O  u0 p2 c0 {2,S} {4,S}
4 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([2.23, 3.29, 3.83, 4.12, 4.34, 4.41, 4.5], 'cal/(mol*K)'),
        H298=(-99.0, 'kcal/mol'),
        S298=(-43.75, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from COH triple-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   OR
   |
   C
  |||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 48,
    label = "C-*R2C=*R",
    group =
"""
1 X  u0 p0 c0 {2,S} {3,S}
2 X  u0 p0 c0 {1,S} {4,D}
3 C  u0 p0 c0 {1,S} {4,S} {5,S} {6,S}
4 C  u0 p0 c0 {2,D} {3,S} {7,S}
5 R  u0 p0 c0 {3,S}
6 R  u0 p0 c0 {3,S}
7 R  u0 p0 c0 {4,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-1.55, 0.65, 2.05, 2.94, 3.91, 4.38, 4.78], 'cal/(mol*K)'),
        H298=(-65.6, 'kcal/mol'),
        S298=(-53.04, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from H2C-CH bidentate, single- and double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R2C--CR
   |  ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 49,
    label = "C-*R2CR3",
    group =
"""
1 X  u0 p0 c0 {2,S}
2 C  u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 C  u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {3,S}
7 R  u0 p0 c0 {3,S}
8 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-1.67, -0.64, 0.12, 0.68, 1.38, 1.78, 2.19], 'cal/(mol*K)'),
        H298=(-41.42, 'kcal/mol'),
        S298=(-38.35, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from H2C-CH3 single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   R
   |
 R-C-CR3
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 50,
    label = "(CR2NR)*",
    group =
"""
1 X  u0 p0 c0
2 C  u0 p0 c0 {3,D} {4,S} {5,S}
3 N  u0 p1 c0 {2,D} {6,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.62, 1.68, 2.22, 2.47, 2.62, 2.62, 2.54], 'cal/(mol*K)'),
        H298=(-5.98, 'kcal/mol'),
        S298=(-33.14, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from H2C-NH vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R2C=NR
    :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 51,
    label = "C-*R2NR2",
    group =
"""
1 X  u0 p0 c0 {2,S}
2 C  u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 N  u0 p1 c0 {2,S} {6,S} {7,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {3,S}
7 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-1.06, -0.27, 0.37, 0.85, 1.47, 1.83, 2.21], 'cal/(mol*K)'),
        H298=(-46.51, 'kcal/mol'),
        S298=(-35.43, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from H2C-NH2 single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   R
   |
 R-C-NR2
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 52,
    label = "(CR2O)*",
    group =
"""
1 X  u0 p0 c0
2 C  u0 p0 c0 {3,D} {4,S} {5,S}
3 O  u0 p2 c0 {2,D}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.12, 1.56, 2.29, 2.63, 2.79, 2.75, 2.6], 'cal/(mol*K)'),
        H298=(-5.18, 'kcal/mol'),
        S298=(-34.64, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from H2C-O vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R2C=O
    :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 53,
    label = "C-*R2OR",
    group =
"""
1 X  u0 p0 c0 {2,S}
2 C  u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 O  u0 p2 c0 {2,S} {6,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-1.37, -0.44, 0.27, 0.8, 1.48, 1.87, 2.27], 'cal/(mol*K)'),
        H298=(-44.42, 'kcal/mol'),
        S298=(-35.7, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from H2C-OH single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

    R
    |
  R-C-OR
    |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 54,
    label = "(CR3NR2)*",
    group =
"""
1 X  u0 p0 c0
2 C  u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 N  u0 p1 c0 {2,S} {7,S} {8,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {2,S}
7 R  u0 p0 c0 {3,S}
8 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.09, 0.79, 1.24, 1.54, 1.88, 2.06, 2.26], 'cal/(mol*K)'),
        H298=(-20.93, 'kcal/mol'),
        S298=(-33.73, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from H3C-NH2 vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R3C-NR2
    :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 55,
    label = "(CR3OR)*",
    group =
"""
1 X  u0 p0 c0
2 C  u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 O  u0 p2 c0 {2,S} {7,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {2,S}
7 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.69, 2.05, 2.23, 2.33, 2.4, 2.43, 2.45], 'cal/(mol*K)'),
        H298=(-7.47, 'kcal/mol'),
        S298=(-28.83, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from H3C-OH vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 R3C-OR
    :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 56,
    label = "C-*RC=*",
    group =
"""
1 X  u0  p0 c0 {2,S} {3,S}
2 X  u0  p0 c0 {1,S} {4,D}
3 C  u0  p0 c0 {1,S} {4,D} {5,S}
4 C  u0  p0 c0 {2,D} {3,D}
5 R  u0  p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.85, 2.25, 3.11, 3.66, 4.25, 4.53, 4.78], 'cal/(mol*K)'),
        H298=(-95.45, 'kcal/mol'),
        S298=(-42.29, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HC-C bidentate, single- and double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 RC--C
  |  ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 57,
    label = "C-*RCR2",
    group =
"""
1 X  u0  p0 c0 {2,S}
2 C  u0  p0 c0 {1,S} {3,D} {4,S}
3 C  u0  p0 c0 {2,D} {5,S} {6,S}
4 R  u0  p0 c0 {2,S}
5 R  u0  p0 c0 {3,S}
6 R  u0  p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.22, 1.96, 3.02, 3.67, 4.35, 4.65, 4.89], 'cal/(mol*K)'),
        H298=(-65.44, 'kcal/mol'),
        S298=(-48.91, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HC-CH2 single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   CR2
  ||
   C-R
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 58,
    label = "C=*RCR3",
    group =
"""
1 X  u0 p0 c0 {3,D}
2 C  u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
3 C  u0 p0 c0 {1,D} {2,S} {7,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
6 R  u0 p0 c0 {2,S}
7 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.11, 2.09, 2.75, 3.2, 3.72, 4.0, 4.37], 'cal/(mol*K)'),
        H298=(-83.24, 'kcal/mol'),
        S298=(-44.11, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HC-CH3 double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   CR3
   |
   C-R
  ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 59,
    label = "(CRN)*",
    group =
"""
1 X  u0  p0 c0
2 C  u0  p0 c0 {3,T} {4,S}
3 N  u0  p1 c0 {2,T}
4 R  u0  p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-0.03, 0.84, 1.41, 1.79, 2.22, 2.4, 2.53], 'cal/(mol*K)'),
        H298=(-0.94, 'kcal/mol'),
        S298=(-22.92, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HCN vdW-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 RC#N
   :
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 60,
    label = "C=*RN=*",
    group =
"""
1 X  u0 p0 c0 {2,S} {3,D}
2 X  u0 p0 c0 {1,S} {4,D}
3 C  u0 p0 c0 {1,D} {4,S} {5,S}
4 N  u0 p1 c0 {2,D} {3,S}
5 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.72, 2.17, 3.13, 3.78, 4.5, 4.82, 5.03], 'cal/(mol*K)'),
        H298=(-15.96, 'kcal/mol'),
        S298=(-35.76, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HCN-h bidentate, twice double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

  R
  |
  C--N
 ||  ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 61,
    label = "C-*RNR",
    group =
"""
1 X  u0  p0 c0 {2,S}
2 C  u0  p0 c0 {1,S} {3,D} {4,S}
3 N  u0  p1 c0 {2,D} {5,S}
4 R  u0  p0 c0 {2,S}
5 R  u0  p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-0.28, 0.62, 1.17, 1.52, 1.89, 2.07, 2.25], 'cal/(mol*K)'),
        H298=(-51.9, 'kcal/mol'),
        S298=(-33.8, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HCNH single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   NR
  ||
   C-R
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 62,
    label = "C=*RN-*R",
    group =
"""
1 X  u0 p0 c0 {2,S} {3,D}
2 X  u0 p0 c0 {1,S} {4,S}
3 C  u0 p0 c0 {1,D} {4,S} {5,S}
4 N  u0 p1 c0 {2,S} {3,S} {6,S}
5 R  u0 p0 c0 {3,S}
6 R  u0 p0 c0 {4,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.09, 2.47, 3.27, 3.76, 4.25, 4.47, 4.67], 'cal/(mol*K)'),
        H298=(-58.44, 'kcal/mol'),
        S298=(-46.17, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HCNH-h bidentate, double- and single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

 RC--NR
 ||  |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 63,
    label = "C=*RNR2",
    group =
"""
1 X  u0 p0 c0 {2,D}
2 C  u0 p0 c0 {1,D} {3,S} {4,S}
3 N  u0 p1 c0 {2,S} {5,S} {6,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {3,S}
6 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([0.47, 1.4, 1.85, 2.07, 2.21, 2.24, 2.27], 'cal/(mol*K)'),
        H298=(-62.35, 'kcal/mol'),
        S298=(-33.34, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HCNH2 double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   NR2
   |
   C-R
  ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 64,
    label = "C-*RO",
    group =
"""
1 X  u0  p0 c0 {2,S}
2 C  u0  p0 c0 {1,S} {3,D} {4,S}
3 O  u0  p2 c0 {2,D}
4 R  u0  p0 c0 {2,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-0.71, 0.27, 0.89, 1.28, 1.69, 1.9, 2.14], 'cal/(mol*K)'),
        H298=(-51.82, 'kcal/mol'),
        S298=(-33.46, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HCO single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   R
   |
   C=O
   |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 65,
    label = "C=*RO-*",
    group =
"""
1 X  u0 p0 c0 {2,S} {3,D}
2 X  u0 p0 c0 {1,S} {4,S}
3 C  u0 p0 c0 {1,D} {4,S} {5,S}
4 O  u0 p2 c0 {2,S} {3,S}
5 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.21, 2.72, 3.6, 4.09, 4.5, 4.63, 4.73], 'cal/(mol*K)'),
        H298=(-44.71, 'kcal/mol'),
        S298=(-45.92, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HCO-h bidentate, double- and single-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

  R
  |
  C--O
 ||  |
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 66,
    label = "C=*ROR",
    group =
"""
1 X  u0 p0 c0 {2,D}
2 C  u0 p0 c0 {1,D} {3,S} {4,S}
3 O  u0 p2 c0 {2,S} {5,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {3,S}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-0.51, 0.59, 1.21, 1.58, 1.97, 2.17, 2.38], 'cal/(mol*K)'),
        H298=(-69.06, 'kcal/mol'),
        S298=(-33.81, 'cal/(mol*K)'),
    ),
    shortDesc="""Came from HCOH double-bonded on Pt(111)""",
    longDesc="""Calculated by Katrin Blondal at Brown University using statistical mechanics (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). Based on DFT calculations by Jelena Jelic at KIT.

   OR
   |
   C-R
  ||
***********
""",
    metal="Pt",
    facet="111",
)

entry(
    index = 67,
    label = "C*",
    group =
"""
1 X  u0 {2,[S,D,T,Q]}
2 C  u0 {1,[S,D,T,Q]}
""",
    thermo='C-*R3',
)

entry(
    index = 68,
    label = "N*",
    group =
"""
1 X  u0 {2,[S,D,T,Q]}
2 N  u0 {1,[S,D,T,Q]}
""",
    thermo='N-*R2',
)

entry(
    index = 69,
    label = "O*",
    group =
"""
1 X  u0 {2,[S,D,T,Q]}
2 O  u0 {1,[S,D,T,Q]}
""",
    thermo='O-*R',
)

entry(
    index = 70,
    label = "R*single_chemisorbed",
    group =
"""
1 X  u0  {2,[S,D,T,Q]}
2 R  u0  {1,[S,D,T,Q]}
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([-0.09, 1.28, 2.17, 2.75, 3.43, 3.79, 4.16], 'cal/(mol*K)'),
        H298=(-45.38, 'kcal/mol'),
        S298=(-38.17, 'cal/(mol*K)'),
    ),
    shortDesc="""Average of C-*R3, N-*R2 and O-*R thermo. """,
    metal="Pt",
    facet="111",
)

entry(
    index = 71,
    label = "C*C*",
    group =
"""
1 X  u0 
2 X  u0 
3 C  u0
4 C  u0
""",
    thermo='C-*R2C-*R2',
)

entry(
    index = 72,
    label = "C*N*",
    group =
"""
1 X  u0 
2 X  u0 
3 C  u0
4 N  u0
""",
    thermo='C=*RN-*R',
)

entry(
    index = 73,
    label = "C*O*",
    group =
"""
1 X  u0 
2 X  u0 
3 C  u0
4 O  u0
""",
    thermo='C=*RO-*R',
)

entry(
    index = 74,
    label = "N*N*",
    group =
"""
1 X  u0 
2 X  u0 
3 N  u0
4 N  u0
""",
    thermo='N-*RN-*R',
)

entry(
    index = 75,
    label = "R*bidentate",
    group =
"""
1 X  u0 
2 X  u0
3 R  u0
4 R  u0
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.69, 3.14, 3.95, 4.38, 4.73, 4.84, 4,88], 'cal/(mol*K)'),
        H298=(-37.29, 'kcal/mol'),
        S298=(-44.37, 'cal/(mol*K)'),
    ),
    shortDesc="""Average of C-*R2C-*R2, C=*RN-*R, C=*RO-* and N-*RN-*R thermo. """,
    metal="Pt",
    facet="111",
)

entry(
    index = 76,
    label = "R*vdW",
    group =
"""
1 X  u0 
2 R  u0
""",
    thermo=ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
        Cpdata=([1.51, 2.1, 2.45, 2.68, 2.94, 3.07, 3.2], 'cal/(mol*K)'),
        H298=(-7.79, 'kcal/mol'),
        S298=(-20.48, 'cal/(mol*K)'),
    ),
    shortDesc="""Average of (CR4)*, (NR3)* and (OR2)* thermo. """,
    metal="Pt",
    facet="111",
)

entry(
    index = 77,
    label = "N*O*",
    group =
"""
1 X  u0 p0 c0 {2,S} {3,[S,D]}
2 X  u0 p0 c0 {1,S} {4,[S,D]}
3 N  u0 p1 c0 {1,[S,D]} {4,[S,D]}
4 O  u0 p2 c0 {2,[S,D]} {3,[S,D]}
""",
    thermo='N=*O-*',
    longDesc="""Is there really any way to do N*O* besides N=*O-* ?""",
    metal="Pt",
    facet="111",
)

entry(
    index = 78,
    label = "O*O*",
    group =
"""
1 X  u0 p0 c0 {2,S} {3,S}
2 X  u0 p0 c0 {1,S} {4,S}
3 O  u0 p2 c0 {1,S} {4,S}
4 O  u0 p2 c0 {2,S} {3,S}
""",
    thermo='O-*O-*',
    longDesc="""Is there really any way to do O*O* besides O-*O-* ?""",
    metal="Pt",
    facet="111",
)

entry(
    index = 79,
    label = "N#*R",
    group =
"""
1 X  u0 {2,T}
2 N  u0 {1,T} {3,[S,D]}
3 R  u0 {2,[S,D]}
""",
    thermo='N*'
)

entry(
    index = 80,
    label = "(CR3)*",
    group =
"""
1 X  u0 
2 C  u0 {3,D} {4,S} {5,S}
3 R  u0 {2,D}
4 R  u0 {2,S}
5 R  u0 {2,S}
""",
    thermo='(CR2NR)*',
    longDesc="""Perhaps should be an average?""",
    metal="Pt",
    facet="111",
)

entry(
    index = 81,
    label = "(CR2)*",
    group =
"""
1 X  u0
2 C  u0 {3,[S,D,T]} {4,[S,D,T]}
3 R  u0 {2,[S,D,T]}
4 R  u0 {2,[S,D,T]}
""",
    thermo='(CRN)*'
)

entry(
    index = 82,
    label = "(NR2CR3)*",
    group =
"""
1 X  u0 p0 c0
2 N  u0 p1 c0 {3,S} {4,S} {5,S}
3 Cs u0 p0 c0 {2,S}
4 R  u0 p0 c0 {2,S}
5 R  u0 p0 c0 {2,S}
""",
    thermo='(NR3)*',
    longDesc="""Do we have data for this?""",
    metal="Pt",
    facet="111",
)

entry(
    index = 83,
    label = "(NR2)*",
    group =
"""
1 X  u0 
2 N  u0 {3,D} {4,S}
3 R  u0 {2,D}
4 R  u0 {2,S}
""",
    thermo='(NRO)*',
    longDesc="""Parent of (RN=O)* and (RN=NR)*. Should it be an average?""",
    metal="Pt",
    facet="111",
)

tree(
"""
L1: R*
    L2: R*bidentate
        L3: C*C*
            L4: C-*C-*
            L4: C=*RC=*R
            L4: C-*R2C-*R2
            L4: C-*R2C=*R
            L4: C-*RC=*
        L3: C*N*
            L4: C=*N-*
            L4: C=*RN=*
            L4: C=*RN-*R
        L3: C*O*
            L4: C=*RO-*
        L3: N*N*
            L4: N-*N-*
            L4: N-*RN-*R
        L3: N*O*
            L4: N=*O-*
        L3: O*O*
            L4: O-*O-*
    L2: R*single_chemisorbed
        L3: C*
            L4: Cq*
            L4: C#*R
                L5: C#*CR3
                L5: C#*NR2
                L5: C#*OR
            L4: C=*R2
                L5: C=*RCR3
                L5: C=*RNR2
                L5: C=*ROR
                L5: C=*CR2
                L5: C=*NR
            L4: C-*R3
                L5: C-*R2CR3
                L5: C-*R2NR2
                L5: C-*R2OR
                L5: C-*RCR2
                L5: C-*RNR
                L5: C-*RO
        L3: N*
            L4: N#*R
            L4: N=*R
                L5: N=*CR3
                L5: N=*NR2
                L5: N=*OR
            L4: N-*R2
                L5: N-*RCR3
                L5: N-*RNR2
                L5: N-*ROR
                L5: N-*CR2
                L5: N-*NR
        L3: O*
            L4: O=*
            L4: O-*R
                L5: O-*CR3
                L5: O-*NR2
                L5: O-*OR
    L2: R*vdW
        L3: (CR4)*
            L4: (CR3CR3)*
            L4: (CR3NR2)*
            L4: (CR3OR)*
        L3: (CR3)*
            L4: (CR2NR)*
            L4: (CR2O)*
        L3: (CR2)*
            L4: (CRN)*
        L3: (NR3)*
            L4: (NR2CR3)*
            L4: (NR2NR2)*
            L4: (NR2OR)*
        L3: (NR2)*
            L4: (NRO)*
            L4: (NRNR)*
        L3: (OR2)*
            L4: (OROR)*

"""
)
