#!/usr/bin/env python
# encoding: utf-8

name = "H_Abstraction/groups"
shortDesc = u""
longDesc = u"""
The reaction site *3 needs a lone pair in order to react. It cannot be 2S or 4S.
"""

template(reactants=["X_H_or_Xrad_H_Xbirad_H_Xtrirad_H", "Y_rad_birad_trirad_quadrad"], products=["X_H_or_Xrad_H_Xbirad_H_Xtrirad_H", "Y_rad_birad_trirad_quadrad"], ownReverse=True)

recipe(actions=[
    ['BREAK_BOND', '*1', 1, '*2'],
    ['FORM_BOND', '*2', 1, '*3'],
    ['GAIN_RADICAL', '*1', '1'],
    ['LOSE_RADICAL', '*3', '1'],
])

entry(
    index = 1,
    label = "X_H_or_Xrad_H_Xbirad_H_Xtrirad_H",
    group = "OR{Xtrirad_H, Xbirad_H, Xrad_H, X_H}",
    kinetics = None,
)

entry(
    index = 2,
    label = "Y_rad_birad_trirad_quadrad",
    group = "OR{Y_rad, Y_1centerbirad, Y_1centertrirad, Y_1centerquadrad}",
    kinetics = None,
)

entry(
    index = 3,
    label = "X_H",
    group = 
"""
1 *1 R u0 {2,S}
2 *2 H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 4,
    label = "H2",
    group = 
"""
1 *1 H u0 {2,S}
2 *2 H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 5,
    label = "Ct_H",
    group = 
"""
1 *1 Ct    u0 {2,S} {3,T}
2 *2 H     u0 {1,S}
3    [C,N] u0 {1,T}
""",
    kinetics = None,
)

entry(
    index = 457,
    label = "Ct/H/NonDeC",
    group = 
"""
1 *1 Ct u0 {2,S} {3,T}
2 *2 H  u0 {1,S}
3    Ct u0 {1,T}
""",
    kinetics = None,
)

entry(
    index = 458,
    label = "Ct/H/NonDeN",
    group = 
"""
1 *1 Ct  u0 {2,S} {3,T}
2 *2 H   u0 {1,S}
3    N3t u0 {1,T}
""",
    kinetics = None,
)

entry(
    index = 6,
    label = "O_H",
    group = 
"""
1 *1 O u0 {2,S} {3,S}
2 *2 H u0 {1,S}
3    R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 7,
    label = "O_pri",
    group = 
"""
1 *1 O u0 {2,S} {3,S}
2 *2 H u0 {1,S}
3    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 8,
    label = "O_sec",
    group = 
"""
1 *1 O   u0 {2,S} {3,S}
2 *2 H   u0 {1,S}
3    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 9,
    label = "O/H/NonDeC",
    group = 
"""
1 *1 O  u0 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 10,
    label = "O/H/NonDeO",
    group = 
"""
1 *1 O u0 {2,S} {3,S}
2 *2 H u0 {1,S}
3    O u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 11,
    label = "H2O2",
    group = 
"""
1 *1 O u0 {2,S} {3,S}
2    O u0 {1,S} {4,S}
3 *2 H u0 {1,S}
4    H u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 12,
    label = "ROOH_pri",
    group = 
"""
1 *1 O  u0 {2,S} {7,S}
2    O  u0 {1,S} {3,S}
3    C  u0 {2,S} {4,S} {5,S} {6,S}
4    Cs u0 {3,S}
5    H  u0 {3,S}
6    H  u0 {3,S}
7 *2 H  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 13,
    label = "ROOH_sec",
    group = 
"""
1 *1 O  u0 {2,S} {7,S}
2    O  u0 {1,S} {3,S}
3    C  u0 {2,S} {4,S} {5,S} {6,S}
4    Cs u0 {3,S}
5    Cs u0 {3,S}
6    H  u0 {3,S}
7 *2 H  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 14,
    label = "ROOH_ter",
    group = 
"""
1 *1 O  u0 {2,S} {7,S}
2    O  u0 {1,S} {3,S}
3    C  u0 {2,S} {4,S} {5,S} {6,S}
4    Cs u0 {3,S}
5    Cs u0 {3,S}
6    Cs u0 {3,S}
7 *2 H  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 471,
    label = "O/H/NonDeN",
    group = 
"""
1 *1 O   u0 {2,S} {3,S}
2 *2 H   u0 {1,S}
3    N3s u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 15,
    label = "O/H/OneDe",
    group = 
"""
1 *1 O                        u0 {2,S} {3,S}
2 *2 H                        u0 {1,S}
3    [Cd,Ct,Cb,CO,CS,N3d,N5d] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 454,
    label = "O/H/OneDeC",
    group = 
"""
1 *1 O                u0 {2,S} {3,S}
2 *2 H                u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 455,
    label = "O/H/OneDeN",
    group = 
"""
1 *1 O         u0 {2,S} {3,S}
2 *2 H         u0 {1,S}
3    [N3d,N5d] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 16,
    label = "Orad_O_H",
    group = 
"""
1 *1 O u0 {2,S} {3,S}
2 *2 H u0 {1,S}
3    O u1 {1,S}
""",
    kinetics = None,
)

entry(
    index = 17,
    label = "S_H",
    group = 
"""
1 *1 S u0 {2,S} {3,S}
2 *2 H u0 {1,S}
3    R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 18,
    label = "S_pri",
    group = 
"""
1 *1 S u0 {2,S} {3,S}
2 *2 H u0 {1,S}
3    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 19,
    label = "S_sec",
    group = 
"""
1 *1 S   u0 {2,S} {3,S}
2 *2 H   u0 {1,S}
3    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 20,
    label = "S/H/NonDeC",
    group = 
"""
1 *1 S  u0 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 21,
    label = "S/H/NonDeS",
    group = 
"""
1 *1 S u0 {2,S} {3,S}
2 *2 H u0 {1,S}
3    S u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 22,
    label = "S/H/OneDe",
    group = 
"""
1 *1 S                u0 {2,S} {3,S}
2 *2 H                u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 23,
    label = "S/H/Ct",
    group = 
"""
1 *1 S  u0 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 24,
    label = "S/H/Cb",
    group = 
"""
1 *1 S  u0 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Cb u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 25,
    label = "S/H/CO",
    group = 
"""
1 *1 S  u0 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    CO u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 26,
    label = "S/H/Cd",
    group = 
"""
1 *1 S  u0 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    Cd u0 {1,S} {4,D}
4    C  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 27,
    label = "S/H/CS",
    group = 
"""
1 *1 S  u0 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    CS u0 {1,S} {4,D}
4    S  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 28,
    label = "Cd_H",
    group = 
"""
1 *1 C     u0 {2,D} {3,S} {4,S}
2    [C,N] u0 {1,D}
3 *2 H     u0 {1,S}
4    R     u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 29,
    label = "Cd_pri",
    group = 
"""
1 *1 C      u0 {2,D} {3,S} {4,S}
2    [Cd,N] u0 {1,D} {5,S}
3 *2 H      u0 {1,S}
4    H      u0 {1,S}
5    R      u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 459,
    label = "Cd/H2/NonDeC",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    Cd u0 {1,D} {5,S}
3 *2 H  u0 {1,S}
4    H  u0 {1,S}
5    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 460,
    label = "Cd/H2/NonDeN",
    group = 
"""
1 *1 C   u0 {2,D} {3,S} {4,S}
2    N3d u0 {1,D} {5,S}
3 *2 H   u0 {1,S}
4    H   u0 {1,S}
5    R   u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 30,
    label = "Cd_sec",
    group = 
"""
1 *1 C      u0 {2,D} {3,S} {4,S}
2    [Cd,N] u0 {1,D} {5,S}
3 *2 H      u0 {1,S}
4    R!H    u0 {1,S}
5    R      u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 31,
    label = "Cd/H/NonDeC",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    Cd u0 {1,D} {5,S}
3 *2 H  u0 {1,S}
4    Cs u0 {1,S}
5    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 32,
    label = "Cd/H/NonDeO",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    Cd u0 {1,D} {5,S}
3 *2 H  u0 {1,S}
4    O  u0 {1,S}
5    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 33,
    label = "Cd/H/NonDeS",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    Cd u0 {1,D} {5,S}
3 *2 H  u0 {1,S}
4    S  u0 {1,S}
5    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 404,
    label = "Cd/H/NonDeN",
    group = 
"""
1 *1 C         u0 {2,D} {3,S} {4,S}
2    Cd        u0 {1,D} {5,S}
3 *2 H         u0 {1,S}
4    [N3s,N5s] u0 {1,S}
5    R         u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 34,
    label = "Cd/H/OneDe",
    group = 
"""
1 *1 C                                             u0 {2,D} {3,S} {4,S}
2    Cd                                            u0 {1,D} {5,S}
3 *2 H                                             u0 {1,S}
4    [Cd,Ct,Cb,CO,CS,N3d,N3t,N3b,N5d,N5dd,N5t,N5b] u0 {1,S}
5    R                                             u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 35,
    label = "Cd/H/Ct",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    Cd u0 {1,D} {5,S}
3 *2 H  u0 {1,S}
4    Ct u0 {1,S}
5    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 36,
    label = "Cd/H/Cb",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    Cd u0 {1,D} {5,S}
3 *2 H  u0 {1,S}
4    Cb u0 {1,S}
5    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 37,
    label = "Cd/H/CO",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    Cd u0 {1,D} {5,S}
3 *2 H  u0 {1,S}
4    CO u0 {1,S}
5    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 38,
    label = "Cd/H/Cd",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    Cd u0 {1,D} {6,S}
3 *2 H  u0 {1,S}
4    Cd u0 {1,S} {5,D}
5    C  u0 {4,D}
6    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 39,
    label = "Cd/H/CS",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    Cd u0 {1,D} {5,S}
3 *2 H  u0 {1,S}
4    CS u0 {1,S}
5    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 405,
    label = "Cd/H/DeN",
    group = 
"""
1 *1 C                              u0 {2,D} {3,S} {4,S}
2    Cd                             u0 {1,D} {5,S}
3 *2 H                              u0 {1,S}
4    [N3d,N3t,N3b,N5d,N5dd,N5t,N5b] u0 {1,S}
5    R                              u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 40,
    label = "Cd_allenic",
    group = 
"""
1 *1 C   u0 {2,D} {3,S} {4,S}
2    Cdd u0 {1,D}
3 *2 H   u0 {1,S}
4    R   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 41,
    label = "Cd_Cdd/H2",
    group = 
"""
1 *1 C   u0 {2,D} {3,S} {4,S}
2    Cdd u0 {1,D}
3 *2 H   u0 {1,S}
4    H   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 40,
    label = "Cb_H",
    group = 
"""
1 *1 Cb       u0 {2,B} {3,B} {4,S}
2    [Cb,Cbf] u0 {1,B}
3    [Cb,Cbf] u0 {1,B}
4 *2 H        u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 41,
    label = "CO_H",
    group = 
"""
1 *1 C u0 {2,D} {3,S} {4,S}
2    O u0 {1,D}
3 *2 H u0 {1,S}
4    R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 42,
    label = "CO_pri",
    group = 
"""
1 *1 C u0 {2,D} {3,S} {4,S}
2    O u0 {1,D}
3 *2 H u0 {1,S}
4    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 43,
    label = "CO_sec",
    group = 
"""
1 *1 C   u0 {2,D} {3,S} {4,S}
2    O   u0 {1,D}
3 *2 H   u0 {1,S}
4    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 44,
    label = "CO/H/NonDe",
    group = 
"""
1 *1 C        u0 {2,D} {3,S} {4,S}
2    O        u0 {1,D}
3 *2 H        u0 {1,S}
4    [Cs,O,S] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 45,
    label = "CO/H/Cs",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    O  u0 {1,D}
3 *2 H  u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 46,
    label = "CO/H/Cs\Cs|Cs",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    O  u0 {1,D}
3 *2 H  u0 {1,S}
4    Cs u0 {1,S} {5,S}
5    Cs u0 {4,S} {6,S}
6    Cs u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 47,
    label = "CO/H/OneDe",
    group = 
"""
1 *1 C             u0 {2,D} {3,S} {4,S}
2    O             u0 {1,D}
3 *2 H             u0 {1,S}
4    [Cd,Ct,Cb,CO] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 48,
    label = "CS_H",
    group = 
"""
1 *1 C u0 {2,D} {3,S} {4,S}
2    S u0 {1,D}
3 *2 H u0 {1,S}
4    R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 49,
    label = "CS_pri",
    group = 
"""
1 *1 C u0 {2,D} {3,S} {4,S}
2    S u0 {1,D}
3 *2 H u0 {1,S}
4    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 50,
    label = "CS_sec",
    group = 
"""
1 *1 C   u0 {2,D} {3,S} {4,S}
2    S   u0 {1,D}
3 *2 H   u0 {1,S}
4    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 51,
    label = "CS/H/NonDeC",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    S  u0 {1,D}
3 *2 H  u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 52,
    label = "CS/H/NonDeO",
    group = 
"""
1 *1 C u0 {2,D} {3,S} {4,S}
2    S u0 {1,D}
3 *2 H u0 {1,S}
4    O u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 53,
    label = "CS/H/NonDeS",
    group = 
"""
1 *1 C u0 {2,D} {3,S} {4,S}
2    S u0 {1,D}
3 *2 H u0 {1,S}
4    S u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 54,
    label = "CS/H/OneDe",
    group = 
"""
1 *1 C                u0 {2,D} {3,S} {4,S}
2    S                u0 {1,D}
3 *2 H                u0 {1,S}
4    [Cd,Ct,Cb,CO,CS] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 55,
    label = "CS/H/Ct",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    S  u0 {1,D}
3 *2 H  u0 {1,S}
4    Ct u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 56,
    label = "CS/H/Cb",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    S  u0 {1,D}
3 *2 H  u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 57,
    label = "CS/H/CO",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    S  u0 {1,D}
3 *2 H  u0 {1,S}
4    CO u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 58,
    label = "CS/H/Cd",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    S  u0 {1,D}
3 *2 H  u0 {1,S}
4    Cd u0 {1,S} {5,D}
5    C  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 59,
    label = "CS/H/CS",
    group = 
"""
1 *1 C  u0 {2,D} {3,S} {4,S}
2    S  u0 {1,D}
3 *2 H  u0 {1,S}
4    CS u0 {1,S} {5,D}
5    S  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 60,
    label = "Cs_H",
    group = 
"""
1 *1 C u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H u0 {1,S}
3    R u0 {1,S}
4    R u0 {1,S}
5    R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 61,
    label = "C_methane",
    group = 
"""
1 *1 C u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H u0 {1,S}
3    H u0 {1,S}
4    H u0 {1,S}
5    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 62,
    label = "C_pri",
    group = 
"""
1 *1 C   u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H   u0 {1,S}
3    H   u0 {1,S}
4    H   u0 {1,S}
5    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 63,
    label = "C/H3/Cs",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 64,
    label = "C/H3/Cs\H3",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2    Cs u0 {1,S} {6,S} {7,S} {8,S}
3 *2 H  u0 {1,S}
4    H  u0 {1,S}
5    H  u0 {1,S}
6    H  u0 {2,S}
7    H  u0 {2,S}
8    H  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 65,
    label = "C/H3/Cs\OneNonDe",
    group = 
"""
1 *1 C        u0 {2,S} {3,S} {4,S} {5,S}
2    Cs       u0 {1,S} {6,S} {7,S} {8,S}
3 *2 H        u0 {1,S}
4    H        u0 {1,S}
5    H        u0 {1,S}
6    [Cs,O,S] u0 {2,S}
7    H        u0 {2,S}
8    H        u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 66,
    label = "C/H3/Cs\H2\Cs",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2    Cs u0 {1,S} {6,S} {7,S} {8,S}
3 *2 H  u0 {1,S}
4    H  u0 {1,S}
5    H  u0 {1,S}
6    Cs u0 {2,S}
7    H  u0 {2,S}
8    H  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 67,
    label = "C/H3/Cs\H2\Cs|O",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2    Cs u0 {1,S} {6,S} {7,S} {8,S}
3 *2 H  u0 {1,S}
4    H  u0 {1,S}
5    H  u0 {1,S}
6    Cs u0 {2,S} {9,S}
7    H  u0 {2,S}
8    H  u0 {2,S}
9    O  u0 {6,S}
""",
    kinetics = None,
)

entry(
    index = 68,
    label = "C/H3/Cs\H2\O",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2    Cs u0 {1,S} {6,S} {7,S} {8,S}
3 *2 H  u0 {1,S}
4    H  u0 {1,S}
5    H  u0 {1,S}
6    O  u0 {2,S}
7    H  u0 {2,S}
8    H  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 69,
    label = "C/H3/Cs\TwoNonDe",
    group = 
"""
1 *1 C        u0 {2,S} {3,S} {4,S} {5,S}
2    Cs       u0 {1,S} {6,S} {7,S} {8,S}
3 *2 H        u0 {1,S}
4    H        u0 {1,S}
5    H        u0 {1,S}
6    [Cs,O,S] u0 {2,S}
7    [Cs,O,S] u0 {2,S}
8    H        u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 70,
    label = "C/H3/Cs\H\Cs\O",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2    Cs u0 {1,S} {6,S} {7,S} {8,S}
3 *2 H  u0 {1,S}
4    H  u0 {1,S}
5    H  u0 {1,S}
6    Cs u0 {2,S}
7    O  u0 {2,S}
8    H  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 71,
    label = "C/H3/Cs\H\Cs\Cs|O",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2    Cs u0 {1,S} {6,S} {7,S} {8,S}
3 *2 H  u0 {1,S}
4    H  u0 {1,S}
5    H  u0 {1,S}
6    Cs u0 {2,S} {9,S}
7    Cs u0 {2,S}
8    H  u0 {2,S}
9    O  u0 {6,S}
""",
    kinetics = None,
)

entry(
    index = 72,
    label = "C/H3/O",
    group = 
"""
1 *1 C u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H u0 {1,S}
3    H u0 {1,S}
4    H u0 {1,S}
5    O u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 73,
    label = "C/H3/S",
    group = 
"""
1 *1 C u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H u0 {1,S}
3    H u0 {1,S}
4    H u0 {1,S}
5    S u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 74,
    label = "C/H3/OneDe",
    group = 
"""
1 *1 C                u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H                u0 {1,S}
3    H                u0 {1,S}
4    H                u0 {1,S}
5    [Cd,Ct,Cb,CO,CS] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 75,
    label = "C/H3/Ct",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    Ct u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 76,
    label = "C/H3/Cb",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    Cb u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 77,
    label = "C/H3/CO",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    CO u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 78,
    label = "C/H3/Cd",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    Cd u0 {1,S} {6,D}
6    C  u0 {5,D}
""",
    kinetics = None,
)

entry(
    index = 79,
    label = "C/H3/Cd\H_Cd\H2",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    Cd u0 {1,S} {6,D} {7,S}
6    Cd u0 {5,D} {8,S} {9,S}
7    H  u0 {5,S}
8    H  u0 {6,S}
9    H  u0 {6,S}
""",
    kinetics = None,
)

entry(
    index = 80,
    label = "C/H3/Cd\H_Cd\H\Cs",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    Cd u0 {1,S} {6,D} {7,S}
6    Cd u0 {5,D} {8,S} {9,S}
7    H  u0 {5,S}
8    Cs u0 {6,S}
9    H  u0 {6,S}
""",
    kinetics = None,
)

entry(
    index = 81,
    label = "C/H3/Cd\Cs_Cd\H2",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    Cd u0 {1,S} {6,D} {7,S}
6    Cd u0 {5,D} {8,S} {9,S}
7    Cs u0 {5,S}
8    H  u0 {6,S}
9    H  u0 {6,S}
""",
    kinetics = None,
)

entry(
    index = 82,
    label = "C/H3/CS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    CS u0 {1,S} {6,D}
6    S  u0 {5,D}
""",
    kinetics = None,
)

entry(
    index = 461,
    label = "Cs/H3/NonDeN",
    group = 
"""
1 *1 C   u0 {2,S} {3,S} {4,S} {5,S}
2    N3s u0 {1,S}
3 *2 H   u0 {1,S}
4    H   u0 {1,S}
5    H   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 471,
    label = "Cs/H3/OneDeN",
    group = 
"""
1 *1 C         u0 {2,S} {3,S} {4,S} {5,S}
2    [N3d,N5d] u0 {1,S}
3 *2 H         u0 {1,S}
4    H         u0 {1,S}
5    H         u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 83,
    label = "C_sec",
    group = 
"""
1 *1 C   u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H   u0 {1,S}
3    H   u0 {1,S}
4    R!H u0 {1,S}
5    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 84,
    label = "C/H2/NonDeC",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 85,
    label = "C/H2/Cs/Cs\O",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S} {6,S}
6    O  u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 86,
    label = "C/H2/Cs/Cs\Cs|O",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S} {6,S}
6    Cs u0 {5,S} {7,S}
7    O  u0 {6,S}
""",
    kinetics = None,
)

entry(
    index = 87,
    label = "C/H2/NonDeC_5ring",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cs u0 {1,S} {6,S}
5    Cs u0 {1,S} {7,S}
6    Cs u0 {4,S} {7,S}
7    Cs u0 {5,S} {6,S}
""",
    kinetics = None,
)

entry(
    index = 88,
    label = "C/H2/NonDeC_5ring_fused6_1",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cs u0 {1,S} {6,S} {8,S}
5    Cs u0 {1,S} {7,S}
6    Cs u0 {4,S} {7,S}
7    Cs u0 {5,S} {6,S} {9,S}
8    Cs u0 {4,S} {9,S}
9    Cs u0 {7,S} {8,S}
""",
    kinetics = None,
)

entry(
    index = 89,
    label = "C/H2/NonDeC_5ring_fused6_2",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cs u0 {1,S} {6,S} {8,S}
5    Cs u0 {1,S} {7,S} {9,S}
6    Cs u0 {4,S} {7,S}
7    Cs u0 {5,S} {6,S}
8    Cs u0 {4,S} {9,S}
9    Cs u0 {5,S} {8,S}
""",
    kinetics = None,
)

entry(
    index = 90,
    label = "C/H2/NonDeC_5ring_alpha6ring",
    group = 
"""
1  *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2  *2 H  u0 {1,S}
3     H  u0 {1,S}
4     Cs u0 {1,S} {6,S} {8,S}
5     Cs u0 {1,S} {7,S}
6     Cs u0 {4,S} {7,S} {11,S}
7     Cs u0 {5,S} {6,S}
8     C  u0 {4,S} {9,S}
9     C  u0 {8,S} {10,S}
10    C  u0 {9,S} {11,S}
11    C  u0 {6,S} {10,S}
""",
    kinetics = None,
)

entry(
    index = 91,
    label = "C/H2/NonDeC_5ring_beta6ring",
    group = 
"""
1  *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2  *2 H  u0 {1,S}
3     H  u0 {1,S}
4     Cs u0 {1,S} {6,S}
5     Cs u0 {1,S} {7,S}
6     Cs u0 {4,S} {7,S} {8,S}
7     Cs u0 {5,S} {6,S} {11,S}
8     C  u0 {6,S} {9,S}
9     C  u0 {8,S} {10,S}
10    C  u0 {9,S} {11,S}
11    C  u0 {7,S} {10,S}
""",
    kinetics = None,
)

entry(
    index = 92,
    label = "C/H2/Cs\H3/Cs\H3",
    group = 
"""
1     Cs u0 {2,S} {4,S} {5,S} {6,S}
2  *1 C  u0 {1,S} {3,S} {7,S} {8,S}
3     Cs u0 {2,S} {9,S} {10,S} {11,S}
4     H  u0 {1,S}
5     H  u0 {1,S}
6     H  u0 {1,S}
7  *2 H  u0 {2,S}
8     H  u0 {2,S}
9     H  u0 {3,S}
10    H  u0 {3,S}
11    H  u0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 93,
    label = "C/H2/NonDeO",
    group = 
"""
1 *1 C        u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H        u0 {1,S}
3    H        u0 {1,S}
4    O        u0 {1,S}
5    [Cs,O,S] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 94,
    label = "C/H2/CsO",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    O  u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 95,
    label = "C/H2/Cs\Cs2/O",
    group = 
"""
1     C  u0 {2,S} {6,S} {7,S} {8,S}
2     Cs u0 {1,S} {3,S} {5,S} {9,S}
3  *1 C  u0 {2,S} {4,S} {10,S} {11,S}
4     O  u0 {3,S} {12,S}
5     C  u0 {2,S} {13,S} {14,S} {15,S}
6     H  u0 {1,S}
7     H  u0 {1,S}
8     H  u0 {1,S}
9     H  u0 {2,S}
10 *2 H  u0 {3,S}
11    H  u0 {3,S}
12    H  u0 {4,S}
13    H  u0 {5,S}
14    H  u0 {5,S}
15    H  u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 96,
    label = "C/H2/O2",
    group = 
"""
1 *1 C u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H u0 {1,S}
3    H u0 {1,S}
4    O u0 {1,S}
5    O u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 97,
    label = "C/H2/NonDeS",
    group = 
"""
1 *1 C        u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H        u0 {1,S}
3    H        u0 {1,S}
4    S        u0 {1,S}
5    [Cs,O,S] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 98,
    label = "C/H2/CsS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    S  u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 407,
    label = "C/H2/NonDeN",
    group = 
"""
1 *1 C          u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H          u0 {1,S}
3    H          u0 {1,S}
4    N          u0 {1,S}
5    [Cs,O,S,N] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 99,
    label = "C/H2/OneDe",
    group = 
"""
1 *1 C                u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H                u0 {1,S}
3    H                u0 {1,S}
4    [Cd,Ct,CO,Cb,CS] u0 {1,S}
5    [Cs,O,S]         u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 100,
    label = "C/H2/OneDeC",
    group = 
"""
1 *1 C                u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H                u0 {1,S}
3    H                u0 {1,S}
4    [Cd,Ct,CO,Cb,CS] u0 {1,S}
5    Cs               u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 101,
    label = "C/H2/CtCs",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Ct u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 102,
    label = "C/H2/CbCs",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cb u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 103,
    label = "C/H2/COCs",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    CO u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 104,
    label = "C/H2/CO\H/Cs\H3",
    group = 
"""
1     Cs u0 {2,S} {5,S} {6,S} {7,S}
2  *1 C  u0 {1,S} {3,S} {8,S} {9,S}
3     CO u0 {2,S} {4,D} {10,S}
4     O  u0 {3,D}
5     H  u0 {1,S}
6     H  u0 {1,S}
7     H  u0 {1,S}
8  *2 H  u0 {2,S}
9     H  u0 {2,S}
10    H  u0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 105,
    label = "C/H2/CdCs",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cd u0 {1,S} {6,D}
5    Cs u0 {1,S}
6    Cd u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 106,
    label = "C/H2/Cd\H_Cd\H2/Cs\H3",
    group = 
"""
1     Cs u0 {2,S} {5,S} {6,S} {7,S}
2  *1 C  u0 {1,S} {3,S} {8,S} {9,S}
3     Cd u0 {2,S} {4,D} {10,S}
4     Cd u0 {3,D} {11,S} {12,S}
5     H  u0 {1,S}
6     H  u0 {1,S}
7     H  u0 {1,S}
8  *2 H  u0 {2,S}
9     H  u0 {2,S}
10    H  u0 {3,S}
11    H  u0 {4,S}
12    H  u0 {4,S}
""",
    kinetics = None,
)

entry(
    index = 107,
    label = "C/H2/CSCs",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    CS u0 {1,S} {6,D}
5    Cs u0 {1,S}
6    S  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 108,
    label = "C/H2/OneDeO",
    group = 
"""
1 *1 C             u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H             u0 {1,S}
3    H             u0 {1,S}
4    [Cd,Ct,CO,Cb] u0 {1,S}
5    O             u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 109,
    label = "C/H2/OneDeS",
    group = 
"""
1 *1 C                u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H                u0 {1,S}
3    H                u0 {1,S}
4    [Cd,Ct,CO,Cb,CS] u0 {1,S}
5    S                u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 110,
    label = "C/H2/CbS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cb u0 {1,S}
5    S  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 111,
    label = "C/H2/CtS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Ct u0 {1,S}
5    S  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 112,
    label = "C/H2/CdS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cd u0 {1,S} {6,D}
5    S  u0 {1,S}
6    C  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 113,
    label = "C/H2/CSS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    CS u0 {1,S} {6,D}
5    S  u0 {1,S}
6    S  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 114,
    label = "C/H2/TwoDe",
    group = 
"""
1 *1 C                u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H                u0 {1,S}
3    H                u0 {1,S}
4    [Cd,Ct,CO,Cb,CS] u0 {1,S}
5    [Cd,Ct,CO,Cb,CS] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 115,
    label = "C/H2/CtCt",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Ct u0 {1,S}
5    Ct u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 116,
    label = "C/H2/CtCb",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Ct u0 {1,S}
5    Cb u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 117,
    label = "C/H2/CtCO",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Ct u0 {1,S}
5    CO u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 118,
    label = "C/H2/CbCb",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cb u0 {1,S}
5    Cb u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 119,
    label = "C/H2/CbCO",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cb u0 {1,S}
5    CO u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 120,
    label = "C/H2/COCO",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    CO u0 {1,S}
5    CO u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 121,
    label = "C/H2/CdCt",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cd u0 {1,S} {6,D}
5    Ct u0 {1,S}
6    C  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 122,
    label = "C/H2/CtCS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Ct u0 {1,S}
5    CS u0 {1,S} {6,D}
6    S  u0 {5,D}
""",
    kinetics = None,
)

entry(
    index = 123,
    label = "C/H2/CdCb",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cd u0 {1,S} {6,D}
5    Cb u0 {1,S}
6    C  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 124,
    label = "C/H2/CbCS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cb u0 {1,S}
5    CS u0 {1,S} {6,D}
6    S  u0 {5,D}
""",
    kinetics = None,
)

entry(
    index = 125,
    label = "C/H2/CdCO",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cd u0 {1,S} {6,D}
5    CO u0 {1,S}
6    C  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 126,
    label = "C/H2/COCS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    CO u0 {1,S}
5    CS u0 {1,S} {6,D}
6    S  u0 {5,D}
""",
    kinetics = None,
)

entry(
    index = 127,
    label = "C/H2/CdCd",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cd u0 {1,S} {6,D}
5    Cd u0 {1,S} {7,D}
6    C  u0 {4,D}
7    C  u0 {5,D}
""",
    kinetics = None,
)

entry(
    index = 128,
    label = "C/H2/CdCS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    Cd u0 {1,S} {6,D}
5    CS u0 {1,S} {7,D}
6    C  u0 {4,D}
7    S  u0 {5,D}
""",
    kinetics = None,
)

entry(
    index = 129,
    label = "C/H2/CSCS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    CS u0 {1,S} {6,D}
5    CS u0 {1,S} {7,D}
6    S  u0 {4,D}
7    S  u0 {5,D}
""",
    kinetics = None,
)

entry(
    index = 131,
    label = "C_ter",
    group = 
"""
1 *1 C   u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H   u0 {1,S}
3    R!H u0 {1,S}
4    R!H u0 {1,S}
5    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 132,
    label = "C/H/NonDe",
    group = 
"""
1 *1 C          u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H          u0 {1,S}
3    [Cs,O,S,N] u0 {1,S}
4    [Cs,O,S,N] u0 {1,S}
5    [Cs,O,S,N] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 133,
    label = "C/H/Cs3",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs u0 {1,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 137,
    label = "C/H/Cs2/Cs\O",
    group = 
"""
1     Cs u0 {2,S} {6,S} {7,S} {8,S}
2  *1 C  u0 {1,S} {3,S} {5,S} {9,S}
3     Cs u0 {2,S} {4,S} {10,S} {11,S}
4     O  u0 {3,S} {12,S}
5     Cs u0 {2,S} {13,S} {14,S} {15,S}
6     H  u0 {1,S}
7     H  u0 {1,S}
8     H  u0 {1,S}
9  *2 H  u0 {2,S}
10    H  u0 {3,S}
11    H  u0 {3,S}
12    H  u0 {4,S}
13    H  u0 {5,S}
14    H  u0 {5,S}
15    H  u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 138,
    label = "C/H/Cs2/Cs\Cs|O",
    group = 
"""
1     Cs u0 {2,S} {6,S} {7,S} {8,S}
2  *1 Cs u0 {1,S} {3,S} {5,S} {9,S}
3     Cs u0 {2,S} {4,S} {10,S} {11,S}
4     Cs u0 {3,S} {12,S} {16,S} {17,S}
5     Cs u0 {2,S} {13,S} {14,S} {15,S}
6     H  u0 {1,S}
7     H  u0 {1,S}
8     H  u0 {1,S}
9  *2 H  u0 {2,S}
10    H  u0 {3,S}
11    H  u0 {3,S}
12    H  u0 {4,S}
13    H  u0 {5,S}
14    H  u0 {5,S}
15    H  u0 {5,S}
16    H  u0 {4,S}
17    O  u0 {4,S}
""",
    kinetics = None,
)

entry(
    index = 134,
    label = "C/H/Cs3_5ring",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs u0 {1,S} {6,S}
4    Cs u0 {1,S} {7,S}
5    Cs u0 {1,S}
6    Cs u0 {3,S} {7,S}
7    Cs u0 {4,S} {6,S}
""",
    kinetics = None,
)

entry(
    index = 135,
    label = "C/H/Cs3_5ring_fused6",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs u0 {1,S} {6,S}
4    Cs u0 {1,S} {7,S}
5    Cs u0 {1,S} {8,S}
6    Cs u0 {3,S} {7,S}
7    Cs u0 {4,S} {6,S} {8,S}
8    Cs u0 {5,S} {7,S}
""",
    kinetics = None,
)

entry(
    index = 136,
    label = "C/H/Cs3_5ring_adj5",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cs u0 {1,S} {6,S}
4    Cs u0 {1,S} {7,S} {9,S}
5    Cs u0 {1,S} {8,S}
6    Cs u0 {3,S} {7,S}
7    Cs u0 {4,S} {6,S}
8    Cs u0 {5,S} {9,S}
9    Cs u0 {4,S} {8,S}
""",
    kinetics = None,
)

entry(
    index = 408,
    label = "C/H/Cs2N",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    N  u0 {1,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 138,
    label = "C/H/NDMustO",
    group = 
"""
1 *1 C      u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H      u0 {1,S}
3    O      u0 {1,S}
4    [Cs,O] u0 {1,S}
5    [Cs,O] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 139,
    label = "C/H/Cs2O",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    O  u0 {1,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 140,
    label = "C/H/CsO2",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    O  u0 {1,S}
4    O  u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 141,
    label = "C/H/O3",
    group = 
"""
1 *1 C u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H u0 {1,S}
3    O u0 {1,S}
4    O u0 {1,S}
5    O u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 142,
    label = "C/H/NDMustS",
    group = 
"""
1 *1 C      u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H      u0 {1,S}
3    S      u0 {1,S}
4    [Cs,S] u0 {1,S}
5    [Cs,S] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 143,
    label = "C/H/Cs2S",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    S  u0 {1,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 144,
    label = "C/H/CsS2",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    S  u0 {1,S}
4    S  u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 145,
    label = "C/H/S3",
    group = 
"""
1 *1 C u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H u0 {1,S}
3    S u0 {1,S}
4    S u0 {1,S}
5    S u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 146,
    label = "C/H/NDMustOS",
    group = 
"""
1 *1 C        u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H        u0 {1,S}
3    O        u0 {1,S}
4    S        u0 {1,S}
5    [Cs,O,S] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 147,
    label = "C/H/CsOS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    O  u0 {1,S}
4    S  u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 148,
    label = "C/H/OneDe",
    group = 
"""
1 *1 C                u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H                u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
4    [Cs,O,S]         u0 {1,S}
5    [Cs,O,S]         u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 149,
    label = "C/H/Cs2",
    group = 
"""
1 *1 C                u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H                u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
4    Cs               u0 {1,S}
5    Cs               u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 150,
    label = "C/H/Cs2Ct",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 151,
    label = "C/H/Cs2Cb",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cb u0 {1,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 152,
    label = "C/H/Cs2CO",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    CO u0 {1,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 153,
    label = "C/H/Cs2Cd",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd u0 {1,S} {6,D}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
6    C  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 154,
    label = "C/H/Cs2CS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    CS u0 {1,S} {6,D}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
6    S  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 155,
    label = "C/H/CsO",
    group = 
"""
1 *1 C             u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H             u0 {1,S}
3    [Cd,Ct,Cb,CO] u0 {1,S}
4    O             u0 {1,S}
5    Cs            u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 156,
    label = "C/H/CsS",
    group = 
"""
1 *1 C                u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H                u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
4    S                u0 {1,S}
5    Cs               u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 157,
    label = "C/H/CbCsS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cb u0 {1,S}
4    S  u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 158,
    label = "C/H/CtCsS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    S  u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 159,
    label = "C/H/CdCsS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd u0 {1,S} {6,D}
4    S  u0 {1,S}
5    Cs u0 {1,S}
6    C  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 160,
    label = "C/H/CSCsS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    CS u0 {1,S} {6,D}
4    S  u0 {1,S}
5    Cs u0 {1,S}
6    S  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 161,
    label = "C/H/OO",
    group = 
"""
1 *1 C             u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H             u0 {1,S}
3    [Cd,Ct,Cb,CO] u0 {1,S}
4    O             u0 {1,S}
5    O             u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 162,
    label = "C/H/OS",
    group = 
"""
1 *1 C             u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H             u0 {1,S}
3    [Cd,Ct,Cb,CO] u0 {1,S}
4    O             u0 {1,S}
5    S             u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 163,
    label = "C/H/SS",
    group = 
"""
1 *1 C             u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H             u0 {1,S}
3    [Cd,Ct,Cb,CO] u0 {1,S}
4    S             u0 {1,S}
5    S             u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 164,
    label = "C/H/TwoDe",
    group = 
"""
1 *1 C                u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H                u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
4    [Cd,Ct,Cb,CO,CS] u0 {1,S}
5    [Cs,O,S]         u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 165,
    label = "C/H/Cs",
    group = 
"""
1 *1 C                u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H                u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
4    [Cd,Ct,Cb,CO,CS] u0 {1,S}
5    Cs               u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 166,
    label = "C/H/CtCt",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    Ct u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 167,
    label = "C/H/CtCb",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    Cb u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 168,
    label = "C/H/CtCO",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    CO u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 169,
    label = "C/H/CbCb",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cb u0 {1,S}
4    Cb u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 170,
    label = "C/H/CbCO",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cb u0 {1,S}
4    CO u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 171,
    label = "C/H/COCO",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    CO u0 {1,S}
4    CO u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 172,
    label = "C/H/CdCt",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd u0 {1,S} {6,D}
4    Ct u0 {1,S}
5    Cs u0 {1,S}
6    C  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 173,
    label = "C/H/CtCS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Ct u0 {1,S}
4    CS u0 {1,S} {6,D}
5    Cs u0 {1,S}
6    S  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 174,
    label = "C/H/CdCb",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd u0 {1,S} {6,D}
4    Cb u0 {1,S}
5    Cs u0 {1,S}
6    C  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 175,
    label = "C/H/CbCS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cb u0 {1,S}
4    CS u0 {1,S}
5    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 176,
    label = "C/H/CdCO",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd u0 {1,S} {6,D}
4    CO u0 {1,S}
5    Cs u0 {1,S}
6    C  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 177,
    label = "C/H/COCS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    CO u0 {1,S}
4    CS u0 {1,S} {6,D}
5    Cs u0 {1,S}
6    S  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 178,
    label = "C/H/CdCd",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd u0 {1,S} {6,D}
4    Cd u0 {1,S} {7,D}
5    Cs u0 {1,S}
6    C  u0 {3,D}
7    C  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 179,
    label = "C/H/CdCS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    Cd u0 {1,S} {6,D}
4    CS u0 {1,S} {7,D}
5    Cs u0 {1,S}
6    C  u0 {3,D}
7    S  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 180,
    label = "C/H/CSCS",
    group = 
"""
1 *1 C  u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H  u0 {1,S}
3    CS u0 {1,S} {6,D}
4    CS u0 {1,S} {7,D}
5    Cs u0 {1,S}
6    S  u0 {3,D}
7    S  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 181,
    label = "C/H/TDMustO",
    group = 
"""
1 *1 C             u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H             u0 {1,S}
3    [Cd,Ct,Cb,CO] u0 {1,S}
4    [Cd,Ct,Cb,CO] u0 {1,S}
5    O             u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 182,
    label = "C/H/TDMustS",
    group = 
"""
1 *1 C             u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H             u0 {1,S}
3    [Cd,Ct,Cb,CO] u0 {1,S}
4    [Cd,Ct,Cb,CO] u0 {1,S}
5    S             u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 183,
    label = "C/H/ThreeDe",
    group = 
"""
1 *1 C             u0 {2,S} {3,S} {4,S} {5,S}
2 *2 H             u0 {1,S}
3    [Cd,Ct,Cb,CO] u0 {1,S}
4    [Cd,Ct,Cb,CO] u0 {1,S}
5    [Cd,Ct,Cb,CO] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 409,
    label = "N3_H",
    group = 
"""
1 *1 [N3s,N3d] u0 {2,S}
2 *2 H         u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 402,
    label = "N3s_H",
    group = 
"""
1 *1 N3s u0 {2,S} {3,S} {4,S}
2 *2 H   u0 {1,S}
3    R   u0 {1,S}
4    R   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 410,
    label = "NH3",
    group = 
"""
1 *1 N3s u0 {2,S} {3,S} {4,S}
2 *2 H   u0 {1,S}
3    H   u0 {1,S}
4    H   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 411,
    label = "N3s_pri_H",
    group = 
"""
1 *1 N3s u0 {2,S} {3,S} {4,S}
2 *2 H   u0 {1,S}
3    H   u0 {1,S}
4    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 450,
    label = "N3s/H2/NonDe",
    group = 
"""
1 *1 N3s         u0 {2,S} {3,S} {4,S}
2 *2 H           u0 {1,S}
3    H           u0 {1,S}
4    [N3s,Cs,Os] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 451,
    label = "N3s/H2/NonDeC",
    group = 
"""
1 *1 N3s u0 {2,S} {3,S} {4,S}
2 *2 H   u0 {1,S}
3    H   u0 {1,S}
4    Cs  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 452,
    label = "N3s/H2/NonDeO",
    group = 
"""
1 *1 N3s u0 {2,S} {3,S} {4,S}
2 *2 H   u0 {1,S}
3    H   u0 {1,S}
4    Os  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 453,
    label = "N3s/H2/NonDeN",
    group = 
"""
1 *1 N3s u0 {2,S} {3,S} {4,S}
2 *2 H   u0 {1,S}
3    H   u0 {1,S}
4    N3s u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 472,
    label = "N3s/H2/OneDe",
    group = 
"""
1 *1 N3s                       u0 {2,S} {3,S} {4,S}
2 *2 H                         u0 {1,S}
3    H                         u0 {1,S}
4    [Cd,Cdd,Ct,CO,CS,N3d,N5d] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 473,
    label = "N3s/H2/OneDeN",
    group = 
"""
1 *1 N3s       u0 {2,S} {3,S} {4,S}
2 *2 H         u0 {1,S}
3    H         u0 {1,S}
4    [N3d,N5d] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 412,
    label = "N3s_sec_H",
    group = 
"""
1 *1 N3s u0 {2,S} {3,S} {4,S}
2 *2 H   u0 {1,S}
3    R!H u0 {1,S}
4    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 413,
    label = "N3d_H",
    group = 
"""
1 *1 N3d u0 {2,S} {3,D}
2 *2 H   u0 {1,S}
3    R!H u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 465,
    label = "N3d/H/NonDe",
    group = 
"""
1 *1 N3d         u0 {2,S} {3,D}
2 *2 H           u0 {1,S}
3    [N3d,Od,Cd] u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 446,
    label = "N3d/H/NonDeC",
    group = 
"""
1 *1 N3d u0 {2,S} {3,D}
2 *2 H   u0 {1,S}
3    Cd  u0 {1,D} {4,S} {5,S}
4    R   u0 {3,S}
5    R   u0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 447,
    label = "N3d/H/NonDeO",
    group = 
"""
1 *1 N3d u0 {2,S} {3,D}
2 *2 H   u0 {1,S}
3    Od  u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 448,
    label = "N3d/H/NonDeN",
    group = 
"""
1 *1 N3d u0 {2,S} {3,D}
2 *2 H   u0 {1,S}
3    N3d u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 466,
    label = "N3d/H/OneDe",
    group = 
"""
1 *1 N3d  u0 {2,S} {3,D}
2 *2 H    u0 {1,S}
3    Cdd  u0 {1,D} {4,D}
4    R!H  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 467,
    label = "N3d/H/CddO",
    group = 
"""
1 *1 N3d  u0 {2,S} {3,D}
2 *2 H    u0 {1,S}
3    Cdd  u0 {1,D} {4,D}
4    Od   u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 414,
    label = "N5_H",
    group = 
"""
1 *1 [N5s,N5d,N5dd,N5t,N5b] u0 {2,S}
2 *2 H                      u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 415,
    label = "N5d_H",
    group = 
"""
1 *1 N5d u0 {2,S}
2 *2 H   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 456,
    label = "N5d/H/NonDeOO",
    group = 
"""
1 *1 N5d u0 {2,S} {3,S} {4,D}
2 *2 H   u0 {1,S}
3    Os  u0 {1,S}
4    Od  u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 185,
    label = "Xrad_H",
    group = 
"""
1 *1 R!H u1 {2,S}
2 *2 H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 469,
    label = "C_rad_H",
    group = 
"""
1 *1 C u1 {2,S}
2 *2 H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 442,
    label = "CH3_rad_H",
    group = 
"""
1 *1 Cs u1 {2,S} {3,S} {4,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 470,
    label = "Cs/H2/OneDeN",
    group = 
"""
1 *1 C         u1 {2,S} {3,S} {4,S}
2 *2 H         u0 {1,S}
3    H         u0 {1,S}
4    [N3d,N5d] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 443,
    label = "OH_rad_H",
    group = 
"""
1 *1 O u1 {2,S}
2 *2 H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 186,
    label = "Srad_H",
    group = 
"""
1 *1 S u1 {2,S}
2 *2 H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 416,
    label = "N3s_rad_H",
    group = 
"""
1 *1 N3s u1 {2,S}
2 *2 H   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 417,
    label = "NH2_rad_H",
    group = 
"""
1 *1 N3s u1 {2,S} {3,S}
2 *2 H   u0 {1,S}
3    H   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 418,
    label = "N3s_rad_H_pri",
    group = 
"""
1 *1 N3s     u1 {2,S} {3,S}
2 *2 H       u0 {1,S}
3    [C,N,O] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 449,
    label = "N3s_rad_H/H/NonDeN",
    group = 
"""
1 *1 N3s u1 {2,S} {3,S}
2 *2 H   u0 {1,S}
3    N3s u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 444,
    label = "Xbirad_H",
    group = "OR{CH2_triplet_H, CH2_singlet_H, NH_triplet_H, NH_singlet_H}",
    kinetics = None,
)

entry(
    index = 476,
    label = "CH2_triplet_H",
    group = 
"""
1 *1 Cs u2 {2,S} {3,S}
2 *2 H  u0 {1,S}
3    H  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 477,
    label = "CH2_singlet_H",
    group = 
"""
1 *1 C u0 p1 {2,S} {3,S}
2 *2 H u0 {1,S}
3    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 478,
    label = "NH_triplet_H",
    group = 
"""
1 *1 N u2 p1 {2,S}
2 *2 H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 479,
    label = "NH_singlet_H",
    group = 
"""
1 *1 N u0 p2 {2,S}
2 *2 H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 474,
    label = "Xtrirad_H",
    group = "OR{C_quartet_H, C_doublet_H}",
    kinetics = None,
)

entry(
    index = 480,
    label = "C_quartet_H",
    group = 
"""
1 *1 C u3 p0 {2,S}
2 *2 H u0 p0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 481,
    label = "C_doublet_H",
    group = 
"""
1 *1 C u1 p1 {2,S}
2 *2 H u0 p0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 475,
    label = "Y_1centerquadrad",
    group = "OR{C_quintet, C_triplet}",
    kinetics = None,
)

entry(
    index = 482,
    label = "C_quintet",
    group = 
"""
1 *3 C u4 p0
""",
    kinetics = None,
)

entry(
    index = 483,
    label = "C_triplet",
    group = 
"""
1 *3 C u2 p1
""",
    kinetics = None,
)

entry(
    index = 419,
    label = "Y_1centertrirad",
    group = "OR{N_atom_quartet, N_atom_doublet, CH_quartet, CH_doublet}",
    kinetics = None,
)

entry(
    index = 485,
    label = "N_atom_quartet",
    group = 
"""
1 *3 N u3 p1
""",
    kinetics = None,
)

entry(
    index = 486,
    label = "N_atom_doublet",
    group = 
"""
1 *3 N u1 p2
""",
    kinetics = None,
)

entry(
    index = 487,
    label = "CH_quartet",
    group = 
"""
1 *3 C u3 p0 {2,S}
2    H u0 p0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 487,
    label = "CH_doublet",
    group = 
"""
1 *3 C u1 p1 {2,S}
2    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 188,
    label = "Y_1centerbirad",
    group = 
"""
1 *3 [Cs,Cd,CO,CS,O,S,N] u2
""",
    kinetics = None,
)

entry(
    index = 189,
    label = "O_atom_triplet",
    group = 
"""
1 *3 O u2
""",
    kinetics = None,
)
 
entry(
    index = 600,
    label = "S_atom_triplet",
    group = 
"""
1 *3 S u2
""",
    kinetics = None,
)

entry(
    index = 190,
    label = "CH2_triplet",
    group = 
"""
1 *3 Cs u2 {2,S} {3,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 422,
    label = "NH_triplet",
    group = 
"""
1 *3 N3s u2 {2,S}
2    H   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 191,
    label = "Y_rad",
    group = 
"""
1 *3 R u1
""",
    kinetics = None,
)

entry(
    index = 192,
    label = "H_rad",
    group = 
"""
1 *3 H u1
""",
    kinetics = None,
)

entry(
    index = 193,
    label = "Y_2centeradjbirad",
    group = 
"""
1 *3 [Ct,Os,Ss] u1 {2,[S,T]}
2    [Ct,Os,Ss] u1 {1,[S,T]}
""",
    kinetics = None,
)

entry(
    index = 194,
    label = "O2b",
    group = 
"""
1 *3 Os u1 {2,S}
2    Os u1 {1,S}
""",
    kinetics = None,
)

entry(
    index = 195,
    label = "C2b",
    group = 
"""
1 *3 Ct u1 {2,T}
2    Ct u1 {1,T}
""",
    kinetics = None,
)

entry(
    index = 196,
    label = "Ct_rad",
    group = 
"""
1 *3 C     u1 {2,T}
2    [C,N] u0 {1,T}
""",
    kinetics = None,
)

entry(
    index = 423,
    label = "Ct_rad/Ct",
    group = 
"""
1 *3 Ct u1 {2,T}
2    Ct u0 {1,T}
""",
    kinetics = None,
)

entry(
    index = 424,
    label = "Ct_rad/N",
    group = 
"""
1 *3 Ct        u1 {2,T}
2    [N3t,N5t] u0 {1,T}
""",
    kinetics = None,
)

entry(
    index = 197,
    label = "O_rad",
    group = 
"""
1 *3 O u1 {2,S}
2    R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 198,
    label = "O_pri_rad",
    group = 
"""
1 *3 O u1 {2,S}
2    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 199,
    label = "O_sec_rad",
    group = 
"""
1 *3 O   u1 {2,S}
2    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 200,
    label = "O_rad/NonDeC",
    group = 
"""
1 *3 O  u1 {2,S}
2    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 201,
    label = "O_rad/Cs\H2\Cs|H|Cs2",
    group = 
"""
1     C  u0 {2,S} {6,S} {7,S} {8,S}
2     C  u0 {1,S} {3,S} {5,S} {9,S}
3     Cs u0 {2,S} {4,S} {10,S} {11,S}
4  *3 O  u1 {3,S}
5     C  u0 {2,S} {12,S} {13,S} {14,S}
6     H  u0 {1,S}
7     H  u0 {1,S}
8     H  u0 {1,S}
9     H  u0 {2,S}
10    H  u0 {3,S}
11    H  u0 {3,S}
12    H  u0 {5,S}
13    H  u0 {5,S}
14    H  u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 202,
    label = "O_rad/NonDeO",
    group = 
"""
1 *3 O u1 {2,S}
2    O u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 203,
    label = "OOC",
    group = 
"""
1 *3 O u1 {2,S}
2    O u0 {1,S} {3,S}
3    C u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 425,
    label = "O_rad/NonDeN",
    group = 
"""
1 *3 O   u1 {2,S}
2    N3s u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 204,
    label = "O_rad/OneDe",
    group = 
"""
1 *3 O                        u1 {2,S}
2    [Cd,Ct,Cb,CO,CS,N3d,N5d] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 491,
    label = "O_rad/OneDeC",
    group = 
"""
1 *3 O             u1 {2,S}
2    [Cd,Ct,Cb,CO] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 205,
    label = "O_rad/Cd",
    group = 
"""
1 *3 O        u1 {2,S}
2    Cd       u0 {1,S} {3,D}
3    [Cd,Cdd] u0 {2,D}
""",
    kinetics = None,
)

entry(
    index = 206,
    label = "O_rad/Cd\H_Cd\H2",
    group = 
"""
1 *3 O  u1 {2,S}
2    Cd u0 {1,S} {3,D} {4,S}
3    Cd u0 {2,D} {5,S} {6,S}
4    H  u0 {2,S}
5    H  u0 {3,S}
6    H  u0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 207,
    label = "O_rad/Cd\H_Cd\H\Cs",
    group = 
"""
1 *3 O  u1 {2,S}
2    Cd u0 {1,S} {3,D} {4,S}
3    Cd u0 {2,D} {5,S} {6,S}
4    H  u0 {2,S}
5    Cs u0 {3,S}
6    H  u0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 208,
    label = "O_rad/Cd\H_Cd\Cs2",
    group = 
"""
1 *3 O  u1 {2,S}
2    Cd u0 {1,S} {3,D} {4,S}
3    Cd u0 {2,D} {5,S} {6,S}
4    H  u0 {2,S}
5    Cs u0 {3,S}
6    Cs u0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 209,
    label = "O_rad/Cd\Cs_Cd\H2",
    group = 
"""
1 *3 O  u1 {2,S}
2    Cd u0 {1,S} {3,D} {4,S}
3    Cd u0 {2,D} {5,S} {6,S}
4    Cs u0 {2,S}
5    H  u0 {3,S}
6    H  u0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 210,
    label = "O_rad/Cd\Cs_Cd\H\Cs",
    group = 
"""
1 *3 O  u1 {2,S}
2    Cd u0 {1,S} {3,D} {4,S}
3    Cd u0 {2,D} {5,S} {6,S}
4    Cs u0 {2,S}
5    Cs u0 {3,S}
6    H  u0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 211,
    label = "O_rad/Cd\Cs_Cd\Cs2",
    group = 
"""
1 *3 O  u1 {2,S}
2    Cd u0 {1,S} {3,D} {4,S}
3    Cd u0 {2,D} {5,S} {6,S}
4    Cs u0 {2,S}
5    Cs u0 {3,S}
6    Cs u0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 426,
    label = "O_rad/OneDeN",
    group = 
"""
1 *3 O         u1 {2,S}
2    [N3d,N5d] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 403,
    label = "InChI=1S/NO3/c2-1(3)4",
    group = 
"""
1 *3 Os  u1 {2,S}
2    N5d u0 {1,S} {3,D} {4,S}
3    Od  u0 {2,D}
4    Os  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 212,
    label = "S_rad",
    group = 
"""
1 *3 S u1 {2,S}
2    R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 213,
    label = "S_pri_rad",
    group = 
"""
1 *3 S u1 {2,S}
2    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 214,
    label = "S_sec_rad",
    group = 
"""
1 *3 S   u1 {2,S}
2    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 215,
    label = "S_rad/NonDeC",
    group = 
"""
1 *3 S  u1 {2,S}
2    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 216,
    label = "S_rad/NonDeS",
    group = 
"""
1 *3 S u1 {2,S}
2    S u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 217,
    label = "S_rad/OneDe",
    group = 
"""
1 *3 S                u1 {2,S}
2    [Cd,Ct,Cb,CO,CS] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 218,
    label = "S_rad/Ct",
    group = 
"""
1 *3 S  u1 {2,S}
2    Ct u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 219,
    label = "S_rad/Cb",
    group = 
"""
1 *3 S  u1 {2,S}
2    Cb u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 220,
    label = "S_rad/CO",
    group = 
"""
1 *3 S  u1 {2,S}
2    CO u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 221,
    label = "S_rad/Cd",
    group = 
"""
1 *3 S  u1 {2,S}
2    Cd u0 {1,S} {3,D}
3    C  u0 {2,D}
""",
    kinetics = None,
)

entry(
    index = 222,
    label = "S_rad/CS",
    group = 
"""
1 *3 S  u1 {2,S}
2    CS u0 {1,S} {3,D}
3    S  u0 {2,D}
""",
    kinetics = None,
)

entry(
    index = 223,
    label = "Cd_rad",
    group = 
"""
1 *3 C u1 {2,D} {3,S}
2    C u0 {1,D}
3    R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 224,
    label = "Cd_pri_rad",
    group = 
"""
1 *3 C  u1 {2,D} {3,S}
2    Cd u0 {1,D} {4,S}
3    H  u0 {1,S}
4    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 225,
    label = "Cd_Cd\H2_pri_rad",
    group = 
"""
1 *3 C  u1 {2,D} {3,S}
2    Cd u0 {1,D} {4,S} {5,S}
3    H  u0 {1,S}
4    H  u0 {2,S}
5    H  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 226,
    label = "Cd_Cd\H\Cs_pri_rad",
    group = 
"""
1 *3 C  u1 {2,D} {3,S}
2    Cd u0 {1,D} {4,S} {5,S}
3    H  u0 {1,S}
4    Cs u0 {2,S}
5    H  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 227,
    label = "Cd_Cd\H\Cs|H2|Cs_pri_rad",
    group = 
"""
1    C  u0 {2,S}
2    Cs u0 {1,S} {3,S} {5,S} {6,S}
3    Cd u0 {2,S} {4,D} {7,S}
4 *3 Cd u1 {3,D} {8,S}
5    H  u0 {2,S}
6    H  u0 {2,S}
7    H  u0 {3,S}
8    H  u0 {4,S}
""",
    kinetics = None,
)

entry(
    index = 228,
    label = "Cd_Cd\Cs2_pri_rad",
    group = 
"""
1 *3 Cd u1 {2,D} {3,S}
2    Cd u0 {1,D} {4,S} {5,S}
3    H  u0 {1,S}
4    Cs u0 {2,S}
5    Cs u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 229,
    label = "Cd_sec_rad",
    group = 
"""
1 *3 Cd  u1 {2,D} {3,S}
2    Cd  u0 {1,D} {4,S}
3    R!H u0 {1,S}
4    R   u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 230,
    label = "Cd_rad/NonDeC",
    group = 
"""
1 *3 Cd u1 {2,D} {3,S}
2    Cd u0 {1,D} {4,S}
3    Cs u0 {1,S}
4    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 231,
    label = "Cd_Cd\H2_rad/Cs",
    group = 
"""
1    Cd u0 {2,D} {4,S} {5,S}
2 *3 Cd u1 {1,D} {3,S}
3    Cs u0 {2,S}
4    H  u0 {1,S}
5    H  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 232,
    label = "Cd_Cd\H\Cs_rad/Cs",
    group = 
"""
1    Cs u0 {2,S} {5,S} {6,S} {7,S}
2 *3 Cd u1 {1,S} {3,D}
3    Cd u0 {2,D} {4,S} {8,S}
4    Cs u0 {3,S}
5    H  u0 {1,S}
6    H  u0 {1,S}
7    H  u0 {1,S}
8    H  u0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 233,
    label = "Cd_rad/NonDeO",
    group = 
"""
1 *3 Cd u1 {2,D} {3,S}
2    Cd u0 {1,D} {4,S}
3    O  u0 {1,S}
4    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 234,
    label = "Cd_rad/NonDeS",
    group = 
"""
1 *3 Cd u1 {2,D} {3,S}
2    Cd u0 {1,D} {4,S}
3    S  u0 {1,S}
4    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 427,
    label = "Cd_rad/NonDeN",
    group = 
"""
1 *3 Cd u1 {2,D} {3,S}
2    Cd u0 {1,D} {4,S}
3    N  u0 {1,S}
4    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 235,
    label = "Cd_rad/OneDe",
    group = 
"""
1 *3 Cd               u1 {2,D} {3,S}
2    Cd               u0 {1,D} {4,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
4    R                u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 236,
    label = "Cd_rad/Ct",
    group = 
"""
1 *3 Cd u1 {2,D} {3,S}
2    Cd u0 {1,D} {4,S}
3    Ct u0 {1,S}
4    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 237,
    label = "Cd_rad/Cb",
    group = 
"""
1 *3 Cd u1 {2,D} {3,S}
2    Cd u0 {1,D} {4,S}
3    Cb u0 {1,S}
4    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 238,
    label = "Cd_rad/CO",
    group = 
"""
1 *3 Cd u1 {2,D} {3,S}
2    Cd u0 {1,D} {4,S}
3    CO u0 {1,S}
4    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 239,
    label = "Cd_rad/Cd",
    group = 
"""
1 *3 Cd u1 {2,D} {3,S}
2    Cd u0 {1,D} {5,S}
3    Cd u0 {1,S} {4,D}
4    C  u0 {3,D}
5    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 240,
    label = "Cd_rad/CS",
    group = 
"""
1 *3 Cd u1 {2,D} {3,S}
2    Cd u0 {1,D} {5,S}
3    CS u0 {1,S} {4,D}
4    S  u0 {3,D}
5    R  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 240,
    label = "Cd_allenic_rad",
    group = 
"""
1 *3 Cd  u1 {2,D} {3,S}
2    Cdd u0 {1,D}
3    R   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 240,
    label = "Cd_Cdd_rad/H",
    group = 
"""
1 *3 Cd  u1 {2,D} {3,S}
2    Cdd u0 {1,D}
3    H   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 241,
    label = "Cb_rad",
    group = 
"""
1 *3 Cb       u1 {2,B} {3,B}
2    [Cb,Cbf] u0 {1,B}
3    [Cb,Cbf] u0 {1,B}
""",
    kinetics = None,
)

entry(
    index = 242,
    label = "CO_rad",
    group = 
"""
1 *3 C u1 {2,D} {3,S}
2    O u0 {1,D}
3    R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 243,
    label = "CO_pri_rad",
    group = 
"""
1 *3 C u1 {2,D} {3,S}
2    O u0 {1,D}
3    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 244,
    label = "CO_sec_rad",
    group = 
"""
1 *3 C   u1 {2,D} {3,S}
2    O   u0 {1,D}
3    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 245,
    label = "CO_rad/NonDe",
    group = 
"""
1 *3 C        u1 {2,D} {3,S}
2    O        u0 {1,D}
3    [Cs,O,S] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 246,
    label = "CO_rad/Cs",
    group = 
"""
1 *3 C  u1 {2,D} {3,S}
2    O  u0 {1,D}
3    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 246,
    label = "CO_rad/OneDe",
    group = 
"""
1 *3 C             u1 {2,D} {3,S}
2    O             u0 {1,D}
3    [Cd,Ct,Cb,CO] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 247,
    label = "CS_rad",
    group = 
"""
1 *3 C u1 {2,D} {3,S}
2    S u0 {1,D}
3    R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 248,
    label = "CS_pri_rad",
    group = 
"""
1 *3 C u1 {2,D} {3,S}
2    S u0 {1,D}
3    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 249,
    label = "CS_sec_rad",
    group = 
"""
1 *3 C   u1 {2,D} {3,S}
2    S   u0 {1,D}
3    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 250,
    label = "CS_rad/NonDe",
    group = 
"""
1 *3 C        u1 {2,D} {3,S}
2    S        u0 {1,D}
3    [Cs,O,S] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 251,
    label = "CS_rad/Cs",
    group = 
"""
1 *3 C  u1 {2,D} {3,S}
2    S  u0 {1,D}
3    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 252,
    label = "CS_rad/O",
    group = 
"""
1 *3 C u1 {2,D} {3,S}
2    S u0 {1,D}
3    O u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 253,
    label = "CS_rad/S",
    group = 
"""
1 *3 C u1 {2,D} {3,S}
2    S u0 {1,D}
3    S u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 254,
    label = "CS_rad/OneDe",
    group = 
"""
1 *3 C                u1 {2,D} {3,S}
2    S                u0 {1,D}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 255,
    label = "CS_rad/Ct",
    group = 
"""
1 *3 C  u1 {2,D} {3,S}
2    S  u0 {1,D}
3    Ct u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 256,
    label = "CS_rad/Cb",
    group = 
"""
1 *3 C  u1 {2,D} {3,S}
2    S  u0 {1,D}
3    Cb u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 257,
    label = "CS_rad/CO",
    group = 
"""
1 *3 C  u1 {2,D} {3,S}
2    S  u0 {1,D}
3    CO u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 258,
    label = "CS_rad/Cd",
    group = 
"""
1 *3 C  u1 {2,D} {3,S}
2    S  u0 {1,D}
3    Cd u0 {1,S} {4,D}
4    C  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 259,
    label = "CS_rad/CS",
    group = 
"""
1 *3 C  u1 {2,D} {3,S}
2    S  u0 {1,D}
3    CS u0 {1,S} {4,D}
4    S  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 260,
    label = "Cs_rad",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    R u0 {1,S}
3    R u0 {1,S}
4    R u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 261,
    label = "C_methyl",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    H u0 {1,S}
4    H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 262,
    label = "C_pri_rad",
    group = 
"""
1 *3 C   u1 {2,S} {3,S} {4,S}
2    H   u0 {1,S}
3    H   u0 {1,S}
4    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 263,
    label = "C_rad/H2/Cs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 264,
    label = "C_rad/H2/Cs\H3",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cs u0 {1,S} {5,S} {6,S} {7,S}
3    H  u0 {1,S}
4    H  u0 {1,S}
5    H  u0 {2,S}
6    H  u0 {2,S}
7    H  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 265,
    label = "C_rad/H2/Cs\Cs2\O",
    group = 
"""
1 *3 C  u1 {2,S} {6,S} {7,S}
2    Cs u0 {1,S} {3,S} {4,S} {5,S}
3    C  u0 {2,S}
4    O  u0 {2,S}
5    C  u0 {2,S}
6    H  u0 {1,S}
7    H  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 266,
    label = "C_rad/H2/Cs\H\Cs\Cs|O",
    group = 
"""
1 *3 C  u1 {2,S} {6,S} {7,S}
2    Cs u0 {1,S} {3,S} {5,S} {8,S}
3    C  u0 {2,S} {4,S}
4    O  u0 {3,S}
5    C  u0 {2,S}
6    H  u0 {1,S}
7    H  u0 {1,S}
8    H  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 267,
    label = "C_rad/H2/Cs\H\Cs|Cs\O",
    group = 
"""
1 *3 C  u1 {2,S} {6,S} {7,S}
2    Cs u0 {1,S} {3,S} {4,S} {8,S}
3    C  u0 {2,S} {5,S}
4    O  u0 {2,S}
5    C  u0 {3,S}
6    H  u0 {1,S}
7    H  u0 {1,S}
8    H  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 268,
    label = "C_rad/H2/Cs\H2\Cs|Cs|O",
    group = 
"""
1 *3 C  u1 {2,S} {6,S} {7,S}
2    Cs u0 {1,S} {3,S} {8,S} {9,S}
3    C  u0 {2,S} {4,S} {5,S}
4    C  u0 {3,S}
5    O  u0 {3,S}
6    H  u0 {1,S}
7    H  u0 {1,S}
8    H  u0 {2,S}
9    H  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 269,
    label = "C_rad/H2/Cs\H2\Cs|Cs#O",
    group = 
"""
1 *3 C  u1 {2,S} {6,S} {7,S}
2    Cs u0 {1,S} {3,S} {8,S} {9,S}
3    C  u0 {2,S} {4,S}
4    C  u0 {3,S} {5,S}
5    O  u0 {4,S}
6    H  u0 {1,S}
7    H  u0 {1,S}
8    H  u0 {2,S}
9    H  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 270,
    label = "C_rad/H2/Ct",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
4    Ct u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 271,
    label = "C_rad/H2/Cb",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 272,
    label = "C_rad/H2/CO",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    H  u0 {1,S}
4    CO u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 273,
    label = "C_rad/H2/O",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    H u0 {1,S}
4    O u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 274,
    label = "C_rad/H2/S",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    H u0 {1,S}
4    S u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 275,
    label = "C_rad/H2/Cd",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    H u0 {1,S}
4    C u0 {1,S} {5,D}
5    C u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 276,
    label = "C_rad/H2/Cd\H_Cd\H2",
    group = 
"""
1 *3 C u1 {2,S} {4,S} {5,S}
2    C u0 {1,S} {3,D} {6,S}
3    C u0 {2,D}
4    H u0 {1,S}
5    H u0 {1,S}
6    H u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 277,
    label = "C_rad/H2/Cd\Cs_Cd\H2",
    group = 
"""
1     C u0 {2,D} {5,S} {6,S}
2     C u0 {1,D} {3,S} {4,S}
3  *3 C u1 {2,S} {7,S} {8,S}
4     C u0 {2,S} {9,S} {10,S} {11,S}
5     H u0 {1,S}
6     H u0 {1,S}
7     H u0 {3,S}
8     H u0 {3,S}
9     H u0 {4,S}
10    H u0 {4,S}
11    H u0 {4,S}
""",
    kinetics = None,
)

entry(
    index = 278,
    label = "C_rad/H2/CS",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    H u0 {1,S}
4    C u0 {1,S} {5,D}
5    S u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 428,
    label = "C_rad/H2/N",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    H u0 {1,S}
4    N u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 279,
    label = "C_sec_rad",
    group = 
"""
1 *3 C   u1 {2,S} {3,S} {4,S}
2    H   u0 {1,S}
3    R!H u0 {1,S}
4    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 280,
    label = "C_rad/H/NonDeC",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cs u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 281,
    label = "C_rad/H/NonDeC_5ring_fused6_1",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cs u0 {1,S} {5,S} {7,S}
4    Cs u0 {1,S} {6,S}
5    Cs u0 {3,S} {6,S}
6    Cs u0 {4,S} {5,S} {8,S}
7    Cs u0 {3,S} {8,S}
8    Cs u0 {6,S} {7,S}
""",
    kinetics = None,
)

entry(
    index = 282,
    label = "C_rad/H/NonDeC_5ring_fused6_2",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cs u0 {1,S} {5,S} {7,S}
4    Cs u0 {1,S} {6,S} {8,S}
5    Cs u0 {3,S} {6,S}
6    Cs u0 {4,S} {5,S}
7    Cs u0 {3,S} {8,S}
8    Cs u0 {4,S} {7,S}
""",
    kinetics = None,
)

entry(
    index = 283,
    label = "C_rad/H/Cs\H3/Cs\H3",
    group = 
"""
1     Cs u0 {2,S} {4,S} {5,S} {6,S}
2  *3 C  u1 {1,S} {3,S} {7,S}
3     Cs u0 {2,S} {8,S} {9,S} {10,S}
4     H  u0 {1,S}
5     H  u0 {1,S}
6     H  u0 {1,S}
7     H  u0 {2,S}
8     H  u0 {3,S}
9     H  u0 {3,S}
10    H  u0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 284,
    label = "C_rad/H/NonDeC_5ring_alpha6ring",
    group = 
"""
1  *3 C  u1 {2,S} {3,S} {4,S}
2     H  u0 {1,S}
3     Cs u0 {1,S} {5,S} {7,S}
4     Cs u0 {1,S} {6,S}
5     Cs u0 {3,S} {6,S} {10,S}
6     Cs u0 {4,S} {5,S}
7     C  u0 {3,S} {8,S}
8     C  u0 {7,S} {9,S}
9     C  u0 {8,S} {10,S}
10    C  u0 {5,S} {9,S}
""",
    kinetics = None,
)

entry(
    index = 285,
    label = "C_rad/H/NonDeC_5ring_beta6ring",
    group = 
"""
1  *3 C  u1 {2,S} {3,S} {4,S}
2     H  u0 {1,S}
3     Cs u0 {1,S} {5,S}
4     Cs u0 {1,S} {6,S}
5     Cs u0 {3,S} {6,S} {7,S}
6     Cs u0 {4,S} {5,S} {10,S}
7     C  u0 {5,S} {8,S}
8     C  u0 {7,S} {9,S}
9     C  u0 {8,S} {10,S}
10    C  u0 {6,S} {9,S}
""",
    kinetics = None,
)

entry(
    index = 288,
    label = "C_rad/H/Cs\H2\CO/Cs",
    group = 
"""
1    Cs u0 {2,S}
2 *3 C  u1 {1,S} {3,S} {5,S}
3    Cs u0 {2,S} {4,S} {6,S} {7,S}
4    CO u0 {3,S}
5    H  u0 {2,S}
6    H  u0 {3,S}
7    H  u0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 286,
    label = "C_rad/H/Cs\H2\Cs/Cs\H2\O",
    group = 
"""
1     C  u0 {2,S} {6,S} {7,S} {8,S}
2     Cs u0 {1,S} {3,S} {9,S} {10,S}
3  *3 C  u1 {2,S} {4,S} {11,S}
4     Cs u0 {3,S} {5,S} {12,S} {13,S}
5     O  u0 {4,S} {14,S}
6     H  u0 {1,S}
7     H  u0 {1,S}
8     H  u0 {1,S}
9     H  u0 {2,S}
10    H  u0 {2,S}
11    H  u0 {3,S}
12    H  u0 {4,S}
13    H  u0 {4,S}
14    H  u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 287,
    label = "C_rad/H/Cs\H\Cs\O/Cs",
    group = 
"""
1     Cs u0 {3,S} {6,S} {7,S} {8,S}
2     Cs u0 {4,S} {9,S} {10,S} {11,S}
3  *3 C  u1 {1,S} {4,S} {12,S}
4     Cs u0 {2,S} {3,S} {5,S} {13,S}
5     Os u0 {4,S} {14,S}
6     H  u0 {1,S}
7     H  u0 {1,S}
8     H  u0 {1,S}
9     H  u0 {2,S}
10    H  u0 {2,S}
11    H  u0 {2,S}
12    H  u0 {3,S}
13    H  u0 {4,S}
14    H  u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 288,
    label = "C_rad/H/Cs\H2\Cs|O/Cs",
    group = 
"""
1     Cs u0 {2,S} {6,S} {7,S} {8,S}
2  *3 C  u1 {1,S} {3,S} {9,S}
3     Cs u0 {2,S} {4,S} {10,S} {11,S}
4     C  u0 {3,S} {5,S} {12,S} {13,S}
5     O  u0 {4,S} {14,S}
6     H  u0 {1,S}
7     H  u0 {1,S}
8     H  u0 {1,S}
9     H  u0 {2,S}
10    H  u0 {3,S}
11    H  u0 {3,S}
12    H  u0 {4,S}
13    H  u0 {4,S}
14    H  u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 289,
    label = "C_rad/H/NonDeO",
    group = 
"""
1 *3 C      u1 {2,S} {3,S} {4,S}
2    H      u0 {1,S}
3    O      u0 {1,S}
4    [Cs,O] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 290,
    label = "C_rad/H/CsO",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cs u0 {1,S}
4    O  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 291,
    label = "C_rad/H/Cs\H2\Cs/O",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cs u0 {1,S} {5,S} {6,S} {7,S}
4    O  u0 {1,S}
5    H  u0 {3,S}
6    H  u0 {3,S}
7    Cs u0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 292,
    label = "C_rad/H/Cs\H2\Cs|H2|Cs/O",
    group = 
"""
1     C  u0 {2,S} {6,S} {7,S} {8,S}
2     Cs u0 {1,S} {3,S} {9,S} {10,S}
3     Cs u0 {2,S} {4,S} {11,S} {12,S}
4  *3 C  u1 {3,S} {5,S} {13,S}
5     O  u0 {4,S} {14,S}
6     H  u0 {1,S}
7     H  u0 {1,S}
8     H  u0 {1,S}
9     H  u0 {2,S}
10    H  u0 {2,S}
11    H  u0 {3,S}
12    H  u0 {3,S}
13    H  u0 {4,S}
14    H  u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 293,
    label = "C_rad/H/Cs\H\Cs2/O",
    group = 
"""
1    C  u0 {2,S}
2    Cs u0 {1,S} {3,S} {5,S} {6,S}
3 *3 C  u1 {2,S} {4,S} {7,S}
4    O  u0 {3,S} {8,S}
5    C  u0 {2,S}
6    H  u0 {2,S}
7    H  u0 {3,S}
8    H  u0 {4,S}
""",
    kinetics = None,
)

entry(
    index = 294,
    label = "C_rad/H/O2",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    O u0 {1,S}
4    O u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 295,
    label = "C_rad/H/NonDeS",
    group = 
"""
1 *3 C      u1 {2,S} {3,S} {4,S}
2    H      u0 {1,S}
3    S      u0 {1,S}
4    [Cs,S] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 296,
    label = "C_rad/H/CsS",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    S  u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 297,
    label = "C_rad/H/S2",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    S u0 {1,S}
4    S u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 429,
    label = "C_rad/H/NonDeCN",
    group = 
"""
1 *3 Cs u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cs u0 {1,S}
4    N  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 431,
    label = "C_rad/H/NonDeON",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    O u0 {1,S}
4    N u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 432,
    label = "C_rad/H/NonDeNN",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    N u0 {1,S}
4    N u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 298,
    label = "C_rad/H/OneDe",
    group = 
"""
1 *3 C                u1 {2,S} {3,S} {4,S}
2    H                u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
4    [Cs,O,S,N]       u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 299,
    label = "C_rad/H/OneDeC",
    group = 
"""
1 *3 C                u1 {2,S} {3,S} {4,S}
2    H                u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
4    Cs               u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 300,
    label = "C_rad/H/CtCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Ct u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 301,
    label = "C_rad/H/CbCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cb u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 302,
    label = "C_rad/H/CO/Cs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cs u0 {1,S}
4    CO u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 303,
    label = "C_rad/H/CO\H/Cs\H3",
    group = 
"""
1    Cs u0 {2,S} {5,S} {6,S} {7,S}
2 *3 C  u1 {1,S} {3,S} {8,S}
3    CO u0 {2,S} {4,D} {9,S}
4    O  u0 {3,D}
5    H  u0 {1,S}
6    H  u0 {1,S}
7    H  u0 {1,S}
8    H  u0 {2,S}
9    H  u0 {3,S}
""",
    kinetics = None,
)

entry(
    index = 304,
    label = "C_rad/H/CdCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cd u0 {1,S} {5,D}
4    Cs u0 {1,S}
5    C  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 305,
    label = "C_rad/H/CSCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    CS u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 306,
    label = "C_rad/H/OneDeO",
    group = 
"""
1 *3 C             u1 {2,S} {3,S} {4,S}
2    H             u0 {1,S}
3    [Cd,Ct,Cb,CO] u0 {1,S}
4    O             u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 307,
    label = "C_rad/H/OneDeS",
    group = 
"""
1 *3 C                u1 {2,S} {3,S} {4,S}
2    H                u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
4    S                u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 308,
    label = "C_rad/H/CtS",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Ct u0 {1,S}
4    S  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 309,
    label = "C_rad/H/CbS",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cb u0 {1,S}
4    S  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 310,
    label = "C_rad/H/CdS",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cd u0 {1,S} {5,D}
4    S  u0 {1,S}
5    C  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 311,
    label = "C_rad/H/CSS",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    CS u0 {1,S} {5,D}
4    S  u0 {1,S}
5    S  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 433,
    label = "C_rad/H/OneDeN",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cd u0 {1,S}
4    N  u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 312,
    label = "C_rad/H/TwoDe",
    group = 
"""
1 *3 C                u1 {2,S} {3,S} {4,S}
2    H                u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
4    [Cd,Ct,Cb,CO,CS] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 313,
    label = "C_rad/H/CtCt",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Ct u0 {1,S}
4    Ct u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 314,
    label = "C_rad/H/CtCb",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Ct u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 315,
    label = "C_rad/H/CtCO",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Ct u0 {1,S}
4    CO u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 316,
    label = "C_rad/H/CbCb",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cb u0 {1,S}
4    Cb u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 317,
    label = "C_rad/H/CbCO",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cb u0 {1,S}
4    CO u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 318,
    label = "C_rad/H/COCO",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    CO u0 {1,S}
4    CO u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 319,
    label = "C_rad/H/CdCt",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cd u0 {1,S} {5,D}
4    Ct u0 {1,S}
5    C  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 320,
    label = "C_rad/H/CtCS",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Ct u0 {1,S}
4    CS u0 {1,S} {5,D}
5    S  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 321,
    label = "C_rad/H/CdCb",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cd u0 {1,S} {5,D}
4    Cb u0 {1,S}
5    C  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 322,
    label = "C_rad/H/CbCS",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cb u0 {1,S}
4    CS u0 {1,S} {5,D}
5    S  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 323,
    label = "C_rad/H/CdCO",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cd u0 {1,S} {5,D}
4    CO u0 {1,S}
5    C  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 324,
    label = "C_rad/H/COCS",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    CO u0 {1,S}
4    CS u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 325,
    label = "C_rad/H/CdCd",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cd u0 {1,S} {5,D}
4    Cd u0 {1,S} {6,D}
5    C  u0 {3,D}
6    C  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 326,
    label = "C_rad/H/CdCS",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    Cd u0 {1,S} {5,D}
4    CS u0 {1,S} {6,D}
5    C  u0 {3,D}
6    S  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 327,
    label = "C_rad/H/CSCS",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    H  u0 {1,S}
3    CS u0 {1,S} {5,D}
4    CS u0 {1,S} {6,D}
5    S  u0 {3,D}
6    S  u0 {4,D}
""",
    kinetics = None,
)

entry(
    index = 328,
    label = "C_ter_rad",
    group = 
"""
1 *3 C   u1 {2,S} {3,S} {4,S}
2    R!H u0 {1,S}
3    R!H u0 {1,S}
4    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 329,
    label = "C_rad/NonDe",
    group = 
"""
1 *3 C        u1 {2,S} {3,S} {4,S}
2    [Cs,O,S] u0 {1,S}
3    [Cs,O,S] u0 {1,S}
4    [Cs,O,S] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 330,
    label = "C_rad/Cs3",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cs u0 {1,S}
3    Cs u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 331,
    label = "C_rad/Cs2/Cs\O",
    group = 
"""
1    Cs u0 {2,S}
2 *3 C  u1 {1,S} {3,S} {5,S}
3    Cs u0 {2,S} {4,S}
4    O  u0 {3,S}
5    Cs u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 332,
    label = "C_rad/Cs3_5ring_fused6",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cs u0 {1,S} {5,S}
3    Cs u0 {1,S} {6,S}
4    Cs u0 {1,S} {7,S}
5    Cs u0 {2,S} {6,S}
6    Cs u0 {3,S} {5,S} {7,S}
7    Cs u0 {4,S} {6,S}
""",
    kinetics = None,
)

entry(
    index = 333,
    label = "C_rad/Cs3_5ring_adj5",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cs u0 {1,S} {5,S}
3    Cs u0 {1,S} {6,S} {8,S}
4    Cs u0 {1,S} {7,S}
5    Cs u0 {2,S} {6,S}
6    Cs u0 {3,S} {5,S}
7    Cs u0 {4,S} {8,S}
8    Cs u0 {3,S} {7,S}
""",
    kinetics = None,
)

entry(
    index = 334,
    label = "C_rad/NDMustO",
    group = 
"""
1 *3 C      u1 {2,S} {3,S} {4,S}
2    O      u0 {1,S}
3    [Cs,O] u0 {1,S}
4    [Cs,O] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 335,
    label = "C_rad/Cs2O",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    O  u0 {1,S}
3    Cs u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 336,
    label = "C_rad/OOH/Cs/Cs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    O  u0 {1,S} {5,S}
3    Cs u0 {1,S}
4    Cs u0 {1,S}
5    O  u0 {2,S}
""",
    kinetics = None,
)

entry(
    index = 337,
    label = "C_rad/O/Cs/Cs\Cs",
    group = 
"""
1     C  u0 {3,S} {6,S} {7,S} {8,S}
2     Cs u0 {4,S} {9,S} {10,S} {11,S}
3     Cs u0 {1,S} {4,S} {12,S} {13,S}
4  *3 C  u1 {2,S} {3,S} {5,S}
5     O  u0 {4,S} {14,S}
6     H  u0 {1,S}
7     H  u0 {1,S}
8     H  u0 {1,S}
9     H  u0 {2,S}
10    H  u0 {2,S}
11    H  u0 {2,S}
12    H  u0 {3,S}
13    H  u0 {3,S}
14    H  u0 {5,S}
""",
    kinetics = None,
)

entry(
    index = 338,
    label = "C_rad/CsO2",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    O  u0 {1,S}
3    O  u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 339,
    label = "C_rad/O3",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    O u0 {1,S}
3    O u0 {1,S}
4    O u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 340,
    label = "C_rad/NDMustS",
    group = 
"""
1 *3 C      u1 {2,S} {3,S} {4,S}
2    S      u0 {1,S}
3    [Cs,S] u0 {1,S}
4    [Cs,S] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 341,
    label = "C_rad/Cs2S",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    S  u0 {1,S}
3    Cs u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 342,
    label = "C_rad/CsS2",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    S  u0 {1,S}
3    S  u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 343,
    label = "C_rad/S3",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    S u0 {1,S}
3    S u0 {1,S}
4    S u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 344,
    label = "C_rad/OneDe",
    group = 
"""
1 *3 C                u1 {2,S} {3,S} {4,S}
2    [Cd,Ct,Cb,CO,CS] u0 {1,S}
3    [Cs,O,S]         u0 {1,S}
4    [Cs,O,S]         u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 345,
    label = "C_rad/Cs2",
    group = 
"""
1 *3 C                u1 {2,S} {3,S} {4,S}
2    [Cd,Ct,Cb,CO,CS] u0 {1,S}
3    Cs               u0 {1,S}
4    Cs               u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 346,
    label = "C_rad/CtCs2",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Ct u0 {1,S}
3    Cs u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 347,
    label = "C_rad/CbCs2",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cb u0 {1,S}
3    Cs u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 348,
    label = "C_rad/COCs2",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    CO u0 {1,S}
3    Cs u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 349,
    label = "C_rad/CdCs2",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cd u0 {1,S} {5,D}
3    Cs u0 {1,S}
4    Cs u0 {1,S}
5    C  u0 {2,D}
""",
    kinetics = None,
)

entry(
    index = 350,
    label = "C_rad/CSCs2",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    CS u0 {1,S}
3    Cs u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 351,
    label = "C_rad/CsO",
    group = 
"""
1 *3 C             u1 {2,S} {3,S} {4,S}
2    [Cd,Ct,Cb,CO] u0 {1,S}
3    O             u0 {1,S}
4    Cs            u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 352,
    label = "C_rad/CsS",
    group = 
"""
1 *3 C                u1 {2,S} {3,S} {4,S}
2    [Cd,Ct,Cb,CO,CS] u0 {1,S}
3    S                u0 {1,S}
4    Cs               u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 353,
    label = "C_rad/CtCsS",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Ct u0 {1,S}
3    S  u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 354,
    label = "C_rad/CbCsS",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cb u0 {1,S}
3    S  u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 355,
    label = "C_rad/CdCsS",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cd u0 {1,S} {5,D}
3    S  u0 {1,S}
4    Cs u0 {1,S}
5    C  u0 {2,D}
""",
    kinetics = None,
)

entry(
    index = 356,
    label = "C_rad/CSCsS",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    CS u0 {1,S}
3    S  u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 357,
    label = "C_rad/O2",
    group = 
"""
1 *3 C             u1 {2,S} {3,S} {4,S}
2    [Cd,Ct,Cb,CO] u0 {1,S}
3    O             u0 {1,S}
4    O             u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 358,
    label = "C_rad/OS",
    group = 
"""
1 *3 C             u1 {2,S} {3,S} {4,S}
2    [Cd,Ct,Cb,CO] u0 {1,S}
3    S             u0 {1,S}
4    O             u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 359,
    label = "C_rad/S2",
    group = 
"""
1 *3 C             u1 {2,S} {3,S} {4,S}
2    [Cd,Ct,Cb,CO] u0 {1,S}
3    S             u0 {1,S}
4    S             u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 360,
    label = "C_rad/TwoDe",
    group = 
"""
1 *3 C                u1 {2,S} {3,S} {4,S}
2    [Cd,Ct,Cb,CO,CS] u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
4    [Cs,O,S]         u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 361,
    label = "C_rad/Cs",
    group = 
"""
1 *3 C                u1 {2,S} {3,S} {4,S}
2    [Cd,Ct,Cb,CO,CS] u0 {1,S}
3    [Cd,Ct,Cb,CO,CS] u0 {1,S}
4    Cs               u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 362,
    label = "C_rad/CtCtCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Ct u0 {1,S}
3    Ct u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 363,
    label = "C_rad/CtCbCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Ct u0 {1,S}
3    Cb u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 364,
    label = "C_rad/CtCOCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Ct u0 {1,S}
3    CO u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 365,
    label = "C_rad/CbCbCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cb u0 {1,S}
3    Cb u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 366,
    label = "C_rad/CbCOCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cb u0 {1,S}
3    CO u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 367,
    label = "C_rad/COCOCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    CO u0 {1,S}
3    CO u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 368,
    label = "C_rad/CdCtCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cd u0 {1,S} {5,D}
3    Ct u0 {1,S}
4    Cs u0 {1,S}
5    C  u0 {2,D}
""",
    kinetics = None,
)

entry(
    index = 369,
    label = "C_rad/CtCSCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Ct u0 {1,S}
3    CS u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 370,
    label = "C_rad/CdCbCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cd u0 {1,S} {5,D}
3    Cb u0 {1,S}
4    Cs u0 {1,S}
5    C  u0 {2,D}
""",
    kinetics = None,
)

entry(
    index = 371,
    label = "C_rad/CbCSCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cb u0 {1,S}
3    CS u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 372,
    label = "C_rad/CdCOCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cd u0 {1,S} {5,D}
3    CO u0 {1,S}
4    Cs u0 {1,S}
5    C  u0 {2,D}
""",
    kinetics = None,
)

entry(
    index = 373,
    label = "C_rad/COCSCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    CO u0 {1,S}
3    CS u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 374,
    label = "C_rad/CdCdCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cd u0 {1,S} {5,D}
3    Cd u0 {1,S} {6,D}
4    Cs u0 {1,S}
5    C  u0 {2,D}
6    C  u0 {3,D}
""",
    kinetics = None,
)

entry(
    index = 375,
    label = "C_rad/CdCSCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    Cd u0 {1,S} {5,D}
3    CS u0 {1,S}
4    Cs u0 {1,S}
5    C  u0 {2,D}
""",
    kinetics = None,
)

entry(
    index = 376,
    label = "C_rad/CSCSCs",
    group = 
"""
1 *3 C  u1 {2,S} {3,S} {4,S}
2    CS u0 {1,S}
3    CS u0 {1,S}
4    Cs u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 377,
    label = "C_rad/TDMustO",
    group = 
"""
1 *3 C             u1 {2,S} {3,S} {4,S}
2    [Cd,Ct,Cb,CO] u0 {1,S}
3    [Cd,Ct,Cb,CO] u0 {1,S}
4    O             u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 378,
    label = "C_rad/TDMustS",
    group = 
"""
1 *3 C             u1 {2,S} {3,S} {4,S}
2    [Cd,Ct,Cb,CO] u0 {1,S}
3    [Cd,Ct,Cb,CO] u0 {1,S}
4    S             u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 379,
    label = "C_rad/ThreeDe",
    group = 
"""
1 *3 C             u1 {2,S} {3,S} {4,S}
2    [Cd,Ct,Cb,CO] u0 {1,S}
3    [Cd,Ct,Cb,CO] u0 {1,S}
4    [Cd,Ct,Cb,CO] u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 434,
    label = "N3_rad",
    group = 
"""
1 *3 [N3s,N3d] u1
""",
    kinetics = None,
)

entry(
    index = 435,
    label = "N3s_rad",
    group = 
"""
1 *3 N3s u1 {2,S} {3,S}
2    R   u0 {1,S}
3    R   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 436,
    label = "NH2_rad",
    group = 
"""
1 *3 N3s u1 {2,S} {3,S}
2    H   u0 {1,S}
3    H   u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 437,
    label = "N3s_rad_pri",
    group = 
"""
1 *3 N3s u1 {2,S} {3,S}
2    H   u0 {1,S}
3    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 438,
    label = "N3s_rad_sec",
    group = 
"""
1 *3 N3s u1 {2,S} {3,S}
2    R!H u0 {1,S}
3    R!H u0 {1,S}
""",
    kinetics = None,
)

entry(
    index = 439,
    label = "N3d_rad",
    group = 
"""
1 *3 N3d u1 {2,D}
2    R!H u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 462,
    label = "N3d_rad/OneDe",
    group = 
"""
1 *3 N3d      u1 {2,D}
2    [Cd,Cdd] u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 463,
    label = "N3d_rad/OneDeC",
    group = 
"""
1 *3 N3d u1 {2,D}
2    Cdd u0 {1,D}
""",
    kinetics = None,
)

entry(
    index = 464,
    label = "N3d_rad/OneDeCdd_O",
    group = 
"""
1 *3 N3d u1 {2,D}
2    Cdd u0 {1,D} {3,D}
3    Od  u0 {2,D}
""",
    kinetics = None,
)

entry(
    index = 440,
    label = "N5_rad",
    group = 
"""
1 *3 [N5s,N5d,N5t] u1
""",
    kinetics = None,
)

entry(
    index = 441,
    label = "N5d_rad",
    group = 
"""
1 *3 N5d u1
""",
    kinetics = None,
)

tree(
"""
L1: X_H_or_Xrad_H_Xbirad_H_Xtrirad_H
    L2: Xtrirad_H
        L3: C_quartet_H
        L3: C_doublet_H
    L2: Xbirad_H
        L3: CH2_triplet_H
        L3: CH2_singlet_H
        L3: NH_triplet_H
        L3: NH_singlet_H
    L2: Xrad_H
        L3: C_rad_H
            L4: CH3_rad_H
            L4: Cs/H2/OneDeN
        L3: OH_rad_H
        L3: Srad_H
        L3: N3s_rad_H
            L4: NH2_rad_H
            L4: N3s_rad_H_pri
                L5: N3s_rad_H/H/NonDeN
    L2: X_H
        L3: H2
        L3: Ct_H
            L4: Ct/H/NonDeC
            L4: Ct/H/NonDeN
        L3: O_H
            L4: O_pri
            L4: O_sec
                L5: O/H/NonDeC
                L5: O/H/NonDeO
                    L6: H2O2
                    L6: ROOH_pri
                    L6: ROOH_sec
                    L6: ROOH_ter
                L5: O/H/NonDeN
                L5: O/H/OneDe
                    L6: O/H/OneDeC
                    L6: O/H/OneDeN
        L3: Orad_O_H
        L3: S_H
            L4: S_pri
            L4: S_sec
                L5: S/H/NonDeC
                L5: S/H/NonDeS
                L5: S/H/OneDe
                    L6: S/H/Ct
                    L6: S/H/Cb
                    L6: S/H/CO
                    L6: S/H/Cd
                    L6: S/H/CS
        L3: Cd_H
            L4: Cd_pri
                L5: Cd/H2/NonDeC
                L5: Cd/H2/NonDeN
            L4: Cd_sec
                L5: Cd/H/NonDeC
                L5: Cd/H/NonDeO
                L5: Cd/H/NonDeS
                L5: Cd/H/NonDeN
                L5: Cd/H/OneDe
                    L6: Cd/H/Ct
                    L6: Cd/H/Cb
                    L6: Cd/H/CO
                    L6: Cd/H/Cd
                    L6: Cd/H/CS
                    L6: Cd/H/DeN
            L4: Cd_allenic
                L5: Cd_Cdd/H2
        L3: Cb_H
        L3: CO_H
            L4: CO_pri
            L4: CO_sec
                L5: CO/H/NonDe
                    L6: CO/H/Cs
                        L7: CO/H/Cs\Cs|Cs
                L5: CO/H/OneDe
        L3: CS_H
            L4: CS_pri
            L4: CS_sec
                L5: CS/H/NonDeC
                L5: CS/H/NonDeO
                L5: CS/H/NonDeS
                L5: CS/H/OneDe
                    L6: CS/H/Ct
                    L6: CS/H/Cb
                    L6: CS/H/CO
                    L6: CS/H/Cd
                    L6: CS/H/CS
        L3: Cs_H
            L4: C_methane
            L4: C_pri
                L5: C/H3/Cs
                    L6: C/H3/Cs\H3
                    L6: C/H3/Cs\OneNonDe
                        L7: C/H3/Cs\H2\Cs
                            L8: C/H3/Cs\H2\Cs|O
                        L7: C/H3/Cs\H2\O
                    L6: C/H3/Cs\TwoNonDe
                        L7: C/H3/Cs\H\Cs\O
                        L7: C/H3/Cs\H\Cs\Cs|O
                L5: C/H3/O
                L5: C/H3/S
                L5: C/H3/OneDe
                    L6: C/H3/Ct
                    L6: C/H3/Cb
                    L6: C/H3/CO
                    L6: C/H3/Cd
                        L7: C/H3/Cd\H_Cd\H2
                        L7: C/H3/Cd\H_Cd\H\Cs
                        L7: C/H3/Cd\Cs_Cd\H2
                    L6: C/H3/CS
                L5: Cs/H3/NonDeN
                L5: Cs/H3/OneDeN
            L4: C_sec
                L5: C/H2/NonDeC
                    L6: C/H2/Cs/Cs\O
                    L6: C/H2/Cs/Cs\Cs|O
                    L6: C/H2/NonDeC_5ring
                        L7: C/H2/NonDeC_5ring_fused6_1
                        L7: C/H2/NonDeC_5ring_fused6_2
                        L7: C/H2/NonDeC_5ring_alpha6ring
                        L7: C/H2/NonDeC_5ring_beta6ring
                    L6: C/H2/Cs\H3/Cs\H3
                L5: C/H2/NonDeO
                    L6: C/H2/CsO
                        L7: C/H2/Cs\Cs2/O
                    L6: C/H2/O2
                L5: C/H2/NonDeS
                    L6: C/H2/CsS
                L5: C/H2/NonDeN
                L5: C/H2/OneDe
                    L6: C/H2/OneDeC
                        L7: C/H2/CtCs
                        L7: C/H2/CbCs
                        L7: C/H2/COCs
                            L8: C/H2/CO\H/Cs\H3
                        L7: C/H2/CdCs
                            L8: C/H2/Cd\H_Cd\H2/Cs\H3
                        L7: C/H2/CSCs
                    L6: C/H2/OneDeO
                    L6: C/H2/OneDeS
                        L7: C/H2/CbS
                        L7: C/H2/CtS
                        L7: C/H2/CdS
                        L7: C/H2/CSS
                L5: C/H2/TwoDe
                    L6: C/H2/CtCt
                    L6: C/H2/CtCb
                    L6: C/H2/CtCO
                    L6: C/H2/CbCb
                    L6: C/H2/CbCO
                    L6: C/H2/COCO
                    L6: C/H2/CdCt
                    L6: C/H2/CtCS
                    L6: C/H2/CdCb
                    L6: C/H2/CbCS
                    L6: C/H2/CdCO
                    L6: C/H2/COCS
                    L6: C/H2/CdCd
                    L6: C/H2/CdCS
                    L6: C/H2/CSCS
            L4: C_ter
                L5: C/H/NonDe
                    L6: C/H/Cs3
                        L7: C/H/Cs2/Cs\O
                        L7: C/H/Cs2/Cs\Cs|O
                        L7: C/H/Cs3_5ring
                            L8: C/H/Cs3_5ring_fused6
                            L8: C/H/Cs3_5ring_adj5
                    L6: C/H/Cs2N
                    L6: C/H/NDMustO
                        L7: C/H/Cs2O
                        L7: C/H/CsO2
                        L7: C/H/O3
                    L6: C/H/NDMustS
                        L7: C/H/Cs2S
                        L7: C/H/CsS2
                        L7: C/H/S3
                    L6: C/H/NDMustOS
                        L7: C/H/CsOS
                L5: C/H/OneDe
                    L6: C/H/Cs2
                        L7: C/H/Cs2Ct
                        L7: C/H/Cs2Cb
                        L7: C/H/Cs2CO
                        L7: C/H/Cs2Cd
                        L7: C/H/Cs2CS
                    L6: C/H/CsO
                    L6: C/H/CsS
                        L7: C/H/CbCsS
                        L7: C/H/CtCsS
                        L7: C/H/CdCsS
                        L7: C/H/CSCsS
                    L6: C/H/OO
                    L6: C/H/OS
                    L6: C/H/SS
                L5: C/H/TwoDe
                    L6: C/H/Cs
                        L7: C/H/CtCt
                        L7: C/H/CtCb
                        L7: C/H/CtCO
                        L7: C/H/CbCb
                        L7: C/H/CbCO
                        L7: C/H/COCO
                        L7: C/H/CdCt
                        L7: C/H/CtCS
                        L7: C/H/CdCb
                        L7: C/H/CbCS
                        L7: C/H/CdCO
                        L7: C/H/COCS
                        L7: C/H/CdCd
                        L7: C/H/CdCS
                        L7: C/H/CSCS
                    L6: C/H/TDMustO
                    L6: C/H/TDMustS
                L5: C/H/ThreeDe
        L3: N3_H
            L4: N3s_H
                L5: NH3
                L5: N3s_pri_H
                    L6: N3s/H2/NonDe
                        L7: N3s/H2/NonDeC
                        L7: N3s/H2/NonDeO
                        L7: N3s/H2/NonDeN
                    L6: N3s/H2/OneDe
                        L7: N3s/H2/OneDeN
                L5: N3s_sec_H
            L4: N3d_H
                L5: N3d/H/NonDe
                    L6: N3d/H/NonDeC
                    L6: N3d/H/NonDeO
                    L6: N3d/H/NonDeN
                L5: N3d/H/OneDe
                    L6: N3d/H/CddO
        L3: N5_H
            L4: N5d_H
                L5: N5d/H/NonDeOO
L1: Y_rad_birad_trirad_quadrad
    L2: Y_1centerquadrad
        L3: C_quintet
        L3: C_triplet
    L2: Y_1centertrirad
        L3: N_atom_quartet
        L3: N_atom_doublet
        L3: CH_quartet
        L3: CH_doublet
    L2: Y_1centerbirad
        L3: O_atom_triplet
        L3: S_atom_triplet
        L3: CH2_triplet
        L3: NH_triplet
    L2: Y_rad
        L3: H_rad
        L3: Y_2centeradjbirad
            L4: O2b
            L4: C2b
        L3: Ct_rad
            L4: Ct_rad/Ct
            L4: Ct_rad/N
        L3: O_rad
            L4: O_pri_rad
            L4: O_sec_rad
                L5: O_rad/NonDeC
                    L6: O_rad/Cs\H2\Cs|H|Cs2
                L5: O_rad/NonDeO
                    L6: OOC
                L5: O_rad/NonDeN
                L5: O_rad/OneDe
                    L6: O_rad/OneDeC
                        L7: O_rad/Cd
                            L8: O_rad/Cd\H_Cd\H2
                            L8: O_rad/Cd\H_Cd\H\Cs
                            L8: O_rad/Cd\H_Cd\Cs2
                            L8: O_rad/Cd\Cs_Cd\H2
                            L8: O_rad/Cd\Cs_Cd\H\Cs
                            L8: O_rad/Cd\Cs_Cd\Cs2
                    L6: O_rad/OneDeN
                        L7: InChI=1S/NO3/c2-1(3)4
        L3: S_rad
            L4: S_pri_rad
            L4: S_sec_rad
                L5: S_rad/NonDeC
                L5: S_rad/NonDeS
                L5: S_rad/OneDe
                    L6: S_rad/Ct
                    L6: S_rad/Cb
                    L6: S_rad/CO
                    L6: S_rad/Cd
                    L6: S_rad/CS
        L3: Cd_rad
            L4: Cd_pri_rad
                L5: Cd_Cd\H2_pri_rad
                L5: Cd_Cd\H\Cs_pri_rad
                    L6: Cd_Cd\H\Cs|H2|Cs_pri_rad
                L5: Cd_Cd\Cs2_pri_rad
            L4: Cd_sec_rad
                L5: Cd_rad/NonDeC
                    L6: Cd_Cd\H2_rad/Cs
                    L6: Cd_Cd\H\Cs_rad/Cs
                L5: Cd_rad/NonDeO
                L5: Cd_rad/NonDeS
                L5: Cd_rad/NonDeN
                L5: Cd_rad/OneDe
                    L6: Cd_rad/Ct
                    L6: Cd_rad/Cb
                    L6: Cd_rad/CO
                    L6: Cd_rad/Cd
                    L6: Cd_rad/CS
            L4: Cd_allenic_rad
                L5: Cd_Cdd_rad/H
        L3: Cb_rad
        L3: CO_rad
            L4: CO_pri_rad
            L4: CO_sec_rad
                L5: CO_rad/NonDe
                    L6: CO_rad/Cs
                L5: CO_rad/OneDe
        L3: CS_rad
            L4: CS_pri_rad
            L4: CS_sec_rad
                L5: CS_rad/NonDe
                    L6: CS_rad/Cs
                    L6: CS_rad/O
                    L6: CS_rad/S
                L5: CS_rad/OneDe
                    L6: CS_rad/Ct
                    L6: CS_rad/Cb
                    L6: CS_rad/CO
                    L6: CS_rad/Cd
                    L6: CS_rad/CS
        L3: Cs_rad
            L4: C_methyl
            L4: C_pri_rad
                L5: C_rad/H2/Cs
                    L6: C_rad/H2/Cs\H3
                    L6: C_rad/H2/Cs\Cs2\O
                    L6: C_rad/H2/Cs\H\Cs\Cs|O
                    L6: C_rad/H2/Cs\H\Cs|Cs\O
                    L6: C_rad/H2/Cs\H2\Cs|Cs|O
                    L6: C_rad/H2/Cs\H2\Cs|Cs#O
                L5: C_rad/H2/Ct
                L5: C_rad/H2/Cb
                L5: C_rad/H2/CO
                L5: C_rad/H2/O
                L5: C_rad/H2/S
                L5: C_rad/H2/Cd
                    L6: C_rad/H2/Cd\H_Cd\H2
                    L6: C_rad/H2/Cd\Cs_Cd\H2
                L5: C_rad/H2/CS
                L5: C_rad/H2/N
            L4: C_sec_rad
                L5: C_rad/H/NonDeC
                    L6: C_rad/H/NonDeC_5ring_fused6_1
                    L6: C_rad/H/NonDeC_5ring_fused6_2
                    L6: C_rad/H/Cs\H3/Cs\H3
                    L6: C_rad/H/NonDeC_5ring_alpha6ring
                    L6: C_rad/H/NonDeC_5ring_beta6ring
                    L6: C_rad/H/Cs\H2\CO/Cs
                    L6: C_rad/H/Cs\H2\Cs/Cs\H2\O
                    L6: C_rad/H/Cs\H\Cs\O/Cs
                    L6: C_rad/H/Cs\H2\Cs|O/Cs
                L5: C_rad/H/NonDeO
                    L6: C_rad/H/CsO
                        L7: C_rad/H/Cs\H2\Cs/O
                            L8: C_rad/H/Cs\H2\Cs|H2|Cs/O
                        L7: C_rad/H/Cs\H\Cs2/O
                    L6: C_rad/H/O2
                L5: C_rad/H/NonDeS
                    L6: C_rad/H/CsS
                    L6: C_rad/H/S2
                L5: C_rad/H/NonDeCN
                L5: C_rad/H/NonDeON
                L5: C_rad/H/NonDeNN
                L5: C_rad/H/OneDe
                    L6: C_rad/H/OneDeC
                        L7: C_rad/H/CtCs
                        L7: C_rad/H/CbCs
                        L7: C_rad/H/CO/Cs
                            L8: C_rad/H/CO\H/Cs\H3
                        L7: C_rad/H/CdCs
                        L7: C_rad/H/CSCs
                    L6: C_rad/H/OneDeO
                    L6: C_rad/H/OneDeS
                        L7: C_rad/H/CtS
                        L7: C_rad/H/CbS
                        L7: C_rad/H/CdS
                        L7: C_rad/H/CSS
                    L6: C_rad/H/OneDeN
                L5: C_rad/H/TwoDe
                    L6: C_rad/H/CtCt
                    L6: C_rad/H/CtCb
                    L6: C_rad/H/CtCO
                    L6: C_rad/H/CbCb
                    L6: C_rad/H/CbCO
                    L6: C_rad/H/COCO
                    L6: C_rad/H/CdCt
                    L6: C_rad/H/CtCS
                    L6: C_rad/H/CdCb
                    L6: C_rad/H/CbCS
                    L6: C_rad/H/CdCO
                    L6: C_rad/H/COCS
                    L6: C_rad/H/CdCd
                    L6: C_rad/H/CdCS
                    L6: C_rad/H/CSCS
            L4: C_ter_rad
                L5: C_rad/NonDe
                    L6: C_rad/Cs3
                        L7: C_rad/Cs2/Cs\O
                        L7: C_rad/Cs3_5ring_fused6
                        L7: C_rad/Cs3_5ring_adj5
                    L6: C_rad/NDMustO
                        L7: C_rad/Cs2O
                            L8: C_rad/OOH/Cs/Cs
                            L8: C_rad/O/Cs/Cs\Cs
                        L7: C_rad/CsO2
                        L7: C_rad/O3
                    L6: C_rad/NDMustS
                        L7: C_rad/Cs2S
                        L7: C_rad/CsS2
                        L7: C_rad/S3
                L5: C_rad/OneDe
                    L6: C_rad/Cs2
                        L7: C_rad/CtCs2
                        L7: C_rad/CbCs2
                        L7: C_rad/COCs2
                        L7: C_rad/CdCs2
                        L7: C_rad/CSCs2
                    L6: C_rad/CsO
                    L6: C_rad/CsS
                        L7: C_rad/CtCsS
                        L7: C_rad/CbCsS
                        L7: C_rad/CdCsS
                        L7: C_rad/CSCsS
                    L6: C_rad/O2
                    L6: C_rad/OS
                    L6: C_rad/S2
                L5: C_rad/TwoDe
                    L6: C_rad/Cs
                        L7: C_rad/CtCtCs
                        L7: C_rad/CtCbCs
                        L7: C_rad/CtCOCs
                        L7: C_rad/CbCbCs
                        L7: C_rad/CbCOCs
                        L7: C_rad/COCOCs
                        L7: C_rad/CdCtCs
                        L7: C_rad/CtCSCs
                        L7: C_rad/CdCbCs
                        L7: C_rad/CbCSCs
                        L7: C_rad/CdCOCs
                        L7: C_rad/COCSCs
                        L7: C_rad/CdCdCs
                        L7: C_rad/CdCSCs
                        L7: C_rad/CSCSCs
                    L6: C_rad/TDMustO
                    L6: C_rad/TDMustS
                L5: C_rad/ThreeDe
        L3: N3_rad
            L4: N3s_rad
                L5: NH2_rad
                L5: N3s_rad_pri
                L5: N3s_rad_sec
            L4: N3d_rad
                L5: N3d_rad/OneDe
                    L6: N3d_rad/OneDeC
                        L7: N3d_rad/OneDeCdd_O
        L3: N5_rad
            L4: N5d_rad
"""
)

forbidden(
    label = "disprop1",
    group = 
"""
1 *1 R u0 {2,S} {3,S}
2    C u1 {1,S}
3 *2 H u0 {1,S}
""",
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

forbidden(
    label = "disprop2",
    group = 
"""
1 *1 R u0 {2,S} {3,S}
2    R u0 {1,S} {4,D}
3 *2 H u0 {1,S}
4    R u0 {2,D} {5,S}
5    R u1 {4,S}
""",
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

forbidden(
    label = "disprop3",
    group = 
"""
1 *1 R u0 {2,S} {3,S}
2    R u0 {1,S} {4,T}
3 *2 H u0 {1,S}
4    R u0 {2,T} {5,S}
5    R u1 {4,S}
""",
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

forbidden(
    label = "disprop4",
    group = 
"""
1 *1 R u0 {2,S} {3,S}
2    R u0 {1,S} {4,D}
3 *2 H u0 {1,S}
4    R u0 {2,D} {5,S}
5    R u0 {4,S} {6,D}
6    R u0 {5,D} {7,S}
7    R u1 {6,S}
""",
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

forbidden(
    label = "disprop5",
    group = 
"""
1 *1 R u0 {2,S} {3,S}
2    R u0 {1,S} {4,D}
3 *2 H u0 {1,S}
4    R u0 {2,D} {5,S}
5    R u0 {4,S} {6,D}
6    R u0 {5,D} {7,S}
7    R u1 {6,S}
""",
    shortDesc = u"""""",
    longDesc = 
u"""

""",
)

