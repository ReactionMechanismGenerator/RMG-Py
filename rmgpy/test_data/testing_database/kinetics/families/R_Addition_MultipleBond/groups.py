#!/usr/bin/env python
# encoding: utf-8

name = "R_Addition_MultipleBond/groups"
shortDesc = u""
longDesc = u"""
The reaction site *3 should be a triplet, otherwise it will react via the 1+2_Cycloaddition family instead.
"""

template(reactants=["R_R", "YJ"], products=["RJ_R_Y"], ownReverse=False)

reverse = "Beta_Scission"

recipe(actions=[
    ['CHANGE_BOND', '*1', -1, '*2'],
    ['FORM_BOND', '*1', 1, '*3'],
    ['GAIN_RADICAL', '*2', '1'],
    ['LOSE_RADICAL', '*3', '1'],
])

entry(
    index = 0,
    label = "R_R",
    group = "OR{Cd_R, Ct_R, Od_R, Sd_R, Nd_R, Nt_R}",
    kinetics = None,
)

entry(
    index = 1,
    label = "YJ",
    group = "OR{HJ, Y_1centerquadrad, Y_1centertrirad, Y_1centerbirad, CJ, OJ, SJ, NJ}",
    kinetics = None,
)

entry(
    index = 2,
    label = "Cd_R",
    group = 
"""
1 *1 C   u0 {2,D}
2 *2 R!H u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 3,
    label = "Ck_O",
    group = 
"""
1 *1 Cdd u0 {2,D} {3,D}
2 *2 Od  u0 {1,D}
3    C   u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 127,
    label = "Ck_Ca",
    group = 
"""
1 *1 Cdd u0 {2,D} {3,D}
2 *2 Cdd u0 {1,D} {4,D}
3    Od  u0 {1,D}
4    C   u0 {2,D}
""",
    kinetics = None,
)

entry(
    index = 838,
    label = "Ct_R",
    group = 
"""
1 *1 Ct  u0 {2,T}
2 *2 R!H u0 {1,T}
""",
    kinetics = None,
)

entry(
    index = 4,
    label = "Od_R",
    group = 
"""
1 *1 Od  u0 {2,D}
2 *2 R!H u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 5,
    label = "Nd_R",
    group = "OR{N1d_R, N3d_R}",
    kinetics = None,
)

entry(
    index = 6,
    label = "N1d_R",
    group = 
"""
1 *1 N1d u0 p2 {2,D}
2 *2 R!H u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 7,
    label = "N3d_R",
    group = 
"""
1 *1 N3d u0 {2,D}
2 *2 R!H u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 8,
    label = "Nt_R",
    group = "OR{N3t_R, N5t_R}",
    kinetics = None,
)

entry(
    index = 9,
    label = "N3t_R",
    group = 
"""
1 *1 N3t u0 {2,T}
2 *2 R!H u0 {1,T}
""",
    kinetics = None,
)

entry(
    index = 10,
    label = "N5t_R",
    group = 
"""
1 *1 N5t u0 {2,T}
2 *2 R!H u0 {1,T}
""",
    kinetics = None,
)

entry(
    index = 11,
    label = "Sd_R",
    group = 
"""
1 *1 Sd  u0 {2,D}
2 *2 R!H u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 12,
    label = "HJ",
    group = 
"""
1 *3 H u1
""",
    kinetics = None,
)

entry(
    index = 13,
    label = "Y_1centerquadrad",
    group = "OR{C_quintet, C_triplet}",
    kinetics = None,
)

entry(
    index = 14,
    label = "C_quintet",
    group = 
"""
1 *3 C u4 p0
""",
    kinetics = None,
)

entry(
    index = 15,
    label = "C_triplet",
    group = 
"""
1 *3 C u2 p1
""",
    kinetics = None,
)

entry(
    index = 16,
    label = "Y_1centertrirad",
    group = "OR{N_atom_quartet, N_atom_doublet, CH_quartet, CH_doublet}",
    kinetics = None,
)

entry(
    index = 17,
    label = "N_atom_quartet",
    group = 
"""
1 *3 N u3 p1
""",
    kinetics = None,
)

entry(
    index = 18,
    label = "N_atom_doublet",
    group = 
"""
1 *3 N u1 p2
""",
    kinetics = None,
)

entry(
    index = 19,
    label = "CH_quartet",
    group = 
"""
1 *3 Cs u3 p0 {2,S}
2    H  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 20,
    label = "CH_doublet",
    group = 
"""
1 *3 C u1 p1 {2,S}
2    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 21,
    label = "Y_1centerbirad",
    group = 
"""
1 *3 R!H u2
""",
    kinetics = None,
)

entry(
    index = 22,
    label = "CJ",
    group = 
"""
1 *3 C u1 p0
""",
    kinetics = None,
)

entry(
    index = 23,
    label = "CbJ",
    group = 
"""
1 *3 Cb u1 p0
""",
    kinetics = None,
)

entry(
    index = 24,
    label = "CtJ",
    group = 
"""
1 *3 Ct  u1 p0 {2,T}
2    R!H u0 {1,T}
""",
    kinetics = None,
)

entry(
    index = 25,
    label = "C2b",
    group = 
"""
1 *3 C u1 p0 {2,T}
2    C u1 {1,T}
""",
    kinetics = None,
)

entry(
    index = 26,
    label = "C=SJ",
    group = 
"""
1 *3 CS u1 p0 {2,S}
2    R  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 27,
    label = "CO_rad",
    group = 
"""
1 *3 C u1 p0 {2,D} {3,S}
2    O u0 {1,D}
3    R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 28,
    label = "CsJ",
    group = 
"""
1 *3 C u1 p0 {2,S} {3,S} {4,S}
2    R u0 {1,S}
3    R u0 {1,S}
4    R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 29,
    label = "OJ",
    group = "OR{OJ_pri, OJ_sec, O2b}",
    kinetics = None,
)

entry(
    index = 30,
    label = "OJ_pri",
    group = 
"""
1 *3 O u1 {2,S}
2    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 31,
    label = "OJ_sec",
    group = 
"""
1 *3 O   u1 {2,S}
2    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 32,
    label = "O2b",
    group = 
"""
1 *3 O u1 {2,S}
2    O u1 {1,S}
""",
    kinetics = None,
)

entry(
    index = 33,
    label = "SJ",
    group = 
"""
1 *3 S u1
""",
    kinetics = None,
)

entry(
    index = 34,
    label = "SsJ",
    group = 
"""
1 *3 Ss u1 {2,S}
2    R  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 35,
    label = "NJ",
    group = "OR{N3J}",
    kinetics = None,
)

entry(
    index = 36,
    label = "N3J",
    group = 
"""
1 *3 [N3s,N3d] u1
""",
    kinetics = None,
)

entry(
    index = 37,
    label = "N3sJ",
    group = 
"""
1 *3 N3s u1
""",
    kinetics = None,
)

entry(
    index = 38,
    label = "N3dJ",
    group = 
"""
1 *3 N3d u1
""",
    kinetics = None,
)

tree(
"""
L1: R_R
    L2: Cd_R
        L3: Ck_O
        L3: Ck_Ca
    L2: Ct_R
    L2: Od_R
    L2: Nd_R
        L3: N1d_R
        L3: N3d_R
    L2: Nt_R
        L3: N3t_R
        L3: N5t_R
    L2: Sd_R
L1: YJ
    L2: HJ
    L2: Y_1centerquadrad
        L3: C_quintet
        L3: C_triplet
    L2: Y_1centertrirad
        L3: N_atom_quartet
        L3: N_atom_doublet
        L3: CH_quartet
        L3: CH_doublet
    L2: Y_1centerbirad
    L2: CJ
        L3: CbJ
        L3: CtJ
        L3: C2b
        L3: C=SJ
        L3: CO_rad
        L3: CsJ
    L2: OJ
        L3: OJ_pri
        L3: OJ_sec
        L3: O2b
    L2: SJ
        L3: SsJ
    L2: NJ
        L3: N3J
            L4: N3sJ
            L4: N3dJ
"""
)

forbidden(
    label = "O2d",
    group = 
"""
1 *1 O u0 {2,D}
2 *2 O u0 {1,D}
""",
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

