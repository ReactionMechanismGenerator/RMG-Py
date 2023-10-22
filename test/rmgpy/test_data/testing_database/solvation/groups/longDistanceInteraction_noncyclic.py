#!/usr/bin/env python
# encoding: utf-8

name = "longDistanceInteraction_noncyclic"
shortDesc = ""
longDesc = """ 

"""

entry(
    index=-1,
    label="R",
    group="""
1 *1 R u0
""",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=1,
    label="int14_gauche",
    group="""
1 *1 [Cs,O2s,Cd,S2s] u0 {2,S}
2 *2 Cs u0 {1,S}
""",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=2,
    label="CsCs",
    group="""
1 *1 Cs u0 {2,S}
2 *2 Cs u0 {1,S}
""",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=3,
    label="CsCs-P",
    group="""
1 *1 Cs                         u0 {2,S} {3,S} {4,S} {5,S}
2 *2 Cs                         u0 {1,S}
3   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
4   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
5   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 3049
B: 2937
E: 3215
L: 2767
A: 3240
""",
)

entry(
    index=4,
    label="CsCs-S",
    group="""
1 *1 Cs                         u0 {2,S} {3,S} {4,S} {5,S}
2 *2 Cs                         u0 {1,S}
3   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
4   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
5   Cs u0 {1,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 2073
B: 2023
E: 2170
L: 1887
A: 2169
""",
)

entry(
    index=5,
    label="CsCs-SS",
    group="""
1 *1 Cs                         u0 {2,S} {3,S} {4,S} {5,S}
2 *2 Cs                         u0 {1,S} {6,S} {7,S} {8,S}
3   Cs                         u0 {1,S}
4   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
5   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
6   Cs                         u0 {2,S}
7   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {2,S}
8   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {2,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 1351
B: 1338
E: 1404
L: 1226
A: 1410
""",
)

entry(
    index=6,
    label="CsCs-ST",
    group="""
1 *1 Cs                         u0 {2,S} {3,S} {4,S} {5,S}
2 *2 Cs                         u0 {1,S} {6,S} {7,S} {8,S}
3   Cs                         u0 {1,S}
4   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
5   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
6   Cs                         u0 {2,S}
7   Cs                         u0 {2,S}
8   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {2,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 262
B: 260
E: 263
L: 245
A: 268
""",
)

entry(
    index=7,
    label="CsCs-SQ",
    group="""
1 *1 Cs                         u0 {2,S} {3,S} {4,S} {5,S}
2 *2 Cs                         u0 {1,S} {6,S} {7,S} {8,S}
3   Cs                         u0 {1,S}
4   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
5   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
6   Cs                         u0 {2,S}
7   Cs                         u0 {2,S}
8   Cs                         u0 {2,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 51
B: 51
E: 51
L: 45
A: 51
""",
)

entry(
    index=8,
    label="CsCs-T",
    group="""
1 *1 Cs                         u0 {2,S} {3,S} {4,S} {5,S}
2 *2 Cs                         u0 {1,S}
3   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 526
B: 514
E: 552
L: 495
A: 567
""",
)

entry(
    index=9,
    label="CsCs-TT",
    group="""
1 *1 Cs                         u0 {2,S} {3,S} {4,S} {5,S}
2 *2 Cs                         u0 {1,S} {6,S} {7,S} {8,S}
3   Cs                         u0 {1,S}
4   Cs                         u0 {1,S}
5   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
6   Cs                         u0 {2,S}
7   Cs                         u0 {2,S}
8   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {2,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 27
B: 27
E: 25
L: 24
A: 27
""",
)

entry(
    index=10,
    label="CsCs-T(TTP)",
    group="""
1 *1 Cs                         u0 {2,S} {3,S} {4,S} {5,S}
2 *2 Cs                         u0 {1,S} {6,S} {7,S} {8,S}
3   Cs                         u0 {1,S} {9,S} {10,S} {11,S}
4   Cs                         u0 {1,S} {12,S} {13,S} {14,S}
5   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
6   Cs                         u0 {2,S}
7   Cs                         u0 {2,S}
8   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {2,S}
9   Cs                         u0 {3,S}
10   Cs                         u0 {3,S}
11   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {3,S}
12   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {4,S}
13   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {4,S}
14   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {4,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 4
B: 4
E: 4
L: 4
A: 4
""",
)

entry(
    index=11,
    label="CsCs-T(TTS)",
    group="""
1 *1 Cs                         u0 {2,S} {3,S} {4,S} {5,S}
2 *2 Cs                         u0 {1,S} {6,S} {7,S} {8,S}
3   Cs                         u0 {1,S} {9,S} {10,S} {11,S}
4   Cs                         u0 {1,S} {12,S} {13,S} {14,S}
5   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
6   Cs                         u0 {2,S}
7   Cs                         u0 {2,S}
8   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {2,S}
9   Cs                         u0 {3,S}
10   Cs                         u0 {3,S}
11   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {3,S}
12   Cs u0 {4,S}
13   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {4,S}
14   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {4,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 2
B: 2
E: 1
L: 1
A: 2
""",
)

entry(
    index=12,
    label="CsCs-TQ",
    group="""
1 *1 Cs                         u0 {2,S} {3,S} {4,S} {5,S}
2 *2 Cs                         u0 {1,S} {6,S} {7,S} {8,S}
3   Cs                         u0 {1,S}
4   Cs                         u0 {1,S}
5   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
6   Cs                         u0 {2,S}
7   Cs                         u0 {2,S}
8   Cs                         u0 {2,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 13
B: 13
E: 12
L: 12
A: 13
""",
)

entry(
    index=13,
    label="CsCs-Q",
    group="""
1 *1 Cs                         u0 {2,S} {3,S} {4,S} {5,S}
2 *2 Cs                         u0 {1,S}
3   Cs u0 {1,S}
4   Cs u0 {1,S}
5   Cs u0 {1,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 75
B: 75
E: 73
L: 64
A: 75
""",
)

entry(
    index=14,
    label="CsCs-QQ",
    group="""
1 *1 Cs                         u0 {2,S} {3,S} {4,S} {5,S}
2 *2 Cs                         u0 {1,S} {6,S} {7,S} {8,S}
3   Cs                         u0 {1,S}
4   Cs                         u0 {1,S}
5   Cs                         u0 {1,S}
6   Cs                         u0 {2,S}
7   Cs                         u0 {2,S}
8   Cs                         u0 {2,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 2
B: 2
E: 2
L: 2
A: 2
""",
)

entry(
    index=15,
    label="OsCs",
    group="""
1 *1 O2s u0 {2,S}
2 *2 Cs u0 {1,S}
""",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=16,
    label="OsCs-P",
    group="""
1 *1 O2s                         u0 {2,S} {3,S}
2 *2 Cs                         u0 {1,S}
3   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 1592
B: 1556
E: 1645
L: 1454
A: 1622
""",
)

entry(
    index=17,
    label="OsCs-S",
    group="""
1 *1 O2s                         u0 {2,S} {3,S}
2 *2 Cs                         u0 {1,S}
3    Cs                         u0 {1,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 10
B: 10
E: 10
L: 10
A: 10
""",
)

entry(
    index=18,
    label="OsCs-SP",
    group="""
1 *1 O2s                         u0 {2,S} {3,S}
2 *2 Cs                         u0 {1,S} {4,S} {5,S} {6,S}
3   Cs                         u0 {1,S}
4   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {2,S}
5   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {2,S}
6   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {2,S}
""",
    solute=SoluteData(
        S=0.05053,
        B=-0.03293,
        E=0.00000,
        L=-0.20990,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 184
B: 169
E: 199
L: 174
A: 204
""",
)

entry(
    index=19,
    label="OsCs-SS",
    group="""
1 *1 O2s                         u0 {2,S} {3,S}
2 *2 Cs                         u0 {1,S} {4,S} {5,S} {6,S}
3   Cs                         u0 {1,S}
4   Cs                         u0 {2,S}
5   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {2,S}
6   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {2,S}
""",
    solute=SoluteData(
        S=0.01821,
        B=-0.05133,
        E=0.00000,
        L=-0.26050,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 231
B: 212
E: 248
L: 221
A: 251
""",
)

entry(
    index=20,
    label="OsCs-ST",
    group="""
1 *1 O2s                         u0 {2,S} {3,S}
2 *2 Cs                         u0 {1,S} {4,S} {5,S} {6,S}
3   Cs                         u0 {1,S}
4   Cs                         u0 {2,S}
5   Cs                         u0 {2,S}
6   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {2,S}
""",
    solute=SoluteData(
        S=0.03351,
        B=0.02734,
        E=0.00000,
        L=-0.21073,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 34
B: 23
E: 40
L: 33
A: 39
""",
)

entry(
    index=21,
    label="OsCs-SQ",
    group="""
1 *1 O2s                         u0 {2,S} {3,S}
2 *2 Cs                         u0 {1,S} {4,S} {5,S} {6,S}
3   Cs                         u0 {1,S}
4   Cs                         u0 {2,S}
5   Cs                         u0 {2,S}
6   Cs                         u0 {2,S}
""",
    solute=SoluteData(
        S=0.03210,
        B=0.04690,
        E=0.00000,
        L=-0.09839,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 10
B: 10
E: 14
L: 9
A: 16
""",
)

entry(
    index=22,
    label="CdCs",
    group="""
1 *1 Cd u0 {2,D} {3,S}
2   Cd u0 {1,D}
3 *2 Cs u0 {1,S}
""",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=23,
    label="CdCs-P",
    group="""
1 *1 Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
4 *2 Cs u0 {1,S}
""",
    solute=SoluteData(
        S=0.03660,
        B=0.03194,
        E=0.00896,
        L=0.25534,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 496
B: 489
E: 510
L: 461
A: 513
""",
)

entry(
    index=24,
    label="CdCs-S",
    group="""
1 *1 Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Cs u0 {1,S}
4 *2 Cs u0 {1,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 0
B: 0
E: 0
L: 0
A: 0
""",
)

entry(
    index=25,
    label="CdCs-SP",
    group="""
1 *1 Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Cs u0 {1,S}
4 *2 Cs u0 {1,S} {5,S} {6,S} {7,S}
5   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {4,S}
6   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {4,S}
7   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {4,S}
""",
    solute=SoluteData(
        S=-0.06923,
        B=0.00629,
        E=0.03339,
        L=0.20427,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 151
B: 145
E: 155
L: 142
A: 155
""",
)

entry(
    index=26,
    label="CdCs-SS",
    group="""
1 *1 Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Cs u0 {1,S}
4 *2 Cs u0 {1,S} {5,S} {6,S} {7,S}
5   Cs                         u0 {4,S}
6   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {4,S}
7   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {4,S}
""",
    solute=SoluteData(
        S=-0.02028,
        B=0.00012,
        E=0.00999,
        L=0.22238,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 46
B: 46
E: 47
L: 45
A: 46
""",
)

entry(
    index=27,
    label="CdCs-ST",
    group="""
1 *1 Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Cs u0 {1,S}
4 *2 Cs u0 {1,S} {5,S} {6,S} {7,S}
5   Cs                         u0 {4,S}
6   Cs                         u0 {4,S}
7   [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {4,S}
""",
    solute=SoluteData(
        S=-0.02097,
        B=-0.00294,
        E=0.01043,
        L=0.08786,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 7
B: 7
E: 7
L: 7
A: 7
""",
)

entry(
    index=28,
    label="CdCs-SQ",
    group="""
1 *1 Cd u0 {2,D} {3,S} {4,S}
2   Cd u0 {1,D}
3   Cs u0 {1,S}
4 *2 Cs u0 {1,S} {5,S} {6,S} {7,S}
5   Cs                         u0 {4,S}
6   Cs                         u0 {4,S}
7   Cs                         u0 {4,S}
""",
    solute=SoluteData(
        S=-0.00597,
        B=0.00942,
        E=0.01083,
        L=0.05317,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 1
B: 1
E: 1
L: 1
A: 1
""",
)

entry(
    index=29,
    label="int15",
    group="""
1 *1 Cs u0 {2,S} {4,S} {5,S}
2    [Cs,O2s,S2s] u0 {1,S} {3,S}
3 *2 Cs u0 {2,S} {6,S} {7,S} {8,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
6    Cs u0 {3,S}
7    Cs u0 {3,S}
8    Cs u0 {3,S}
""",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=30,
    label="CsCsCs",
    group="""
1 *1 Cs u0 {2,S} {4,S} {5,S}
2    Cs u0 {1,S} {3,S}
3 *2 Cs u0 {2,S} {6,S} {7,S} {8,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
6    Cs u0 {3,S}
7    Cs u0 {3,S}
8    Cs u0 {3,S}
""",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=31,
    label="CsCsCs-TQ",
    group="""
1 *1 Cs u0 {2,S} {4,S} {5,S} {9,S}
2    Cs u0 {1,S} {3,S}
3 *2 Cs u0 {2,S} {6,S} {7,S} {8,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
6    Cs u0 {3,S}
7    Cs u0 {3,S}
8    Cs u0 {3,S}
9    [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 17
B: 17
E: 16
L: 13
A: 17
""",
)

entry(
    index=32,
    label="CsCsCs-QQ",
    group="""
1 *1 Cs u0 {2,S} {4,S} {5,S} {6,S}
2    Cs u0 {1,S} {3,S}
3 *2 Cs u0 {2,S} {7,S} {8,S} {9,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
6    Cs u0 {1,S}
7    Cs u0 {3,S}
8    Cs u0 {3,S}
9    Cs u0 {3,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 4
B: 4
E: 2
L: 1
A: 4
""",
)

entry(
    index=33,
    label="CsOsCs",
    group="""
1 *1 Cs u0 {2,S} {4,S} {5,S}
2    O2s u0 {1,S} {3,S}
3 *2 Cs u0 {2,S} {6,S} {7,S} {8,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
6    Cs u0 {3,S}
7    Cs u0 {3,S}
8    Cs u0 {3,S}
""",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=34,
    label="CsOsCs-TQ",
    group="""
1 *1 Cs u0 {2,S} {4,S} {5,S} {9,S}
2    O2s u0 {1,S} {3,S}
3 *2 Cs u0 {2,S} {6,S} {7,S} {8,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
6    Cs u0 {3,S}
7    Cs u0 {3,S}
8    Cs u0 {3,S}
9    [Cd,Cdd,Ct,Cb,Cbf,O2s,CO,H] u0 {1,S}
""",
    solute=SoluteData(
        S=-0.00592,
        B=-0.00434,
        E=-0.04000,
        L=-0.05596,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 1
B: 1
E: 1
L: 1
A: 1
""",
)

entry(
    index=35,
    label="CsOsCs-QQ",
    group="""
1 *1 Cs u0 {2,S} {4,S} {5,S} {6,S}
2    O2s u0 {1,S} {3,S}
3 *2 Cs u0 {2,S} {7,S} {8,S} {9,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
6    Cs u0 {1,S}
7    Cs u0 {3,S}
8    Cs u0 {3,S}
9    Cs u0 {3,S}
""",
    solute=SoluteData(
        S=0.00000,
        B=0.00000,
        E=0.00000,
        L=0.00000,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 1
B: 1
E: 1
L: 0
A: 2
""",
)

entry(
    index=36,
    label="CsSsCs",
    group="""
1 *1 Cs u0 {2,S} {4,S} {5,S}
2    S2s u0 {1,S} {3,S}
3 *2 Cs u0 {2,S} {6,S} {7,S} {8,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
6    Cs u0 {3,S}
7    Cs u0 {3,S}
8    Cs u0 {3,S}
""",
    solute=None,
    shortDesc="""""",
    longDesc="""
""",
)

entry(
    index=37,
    label="CsSsCs-QQ",
    group="""
1 *1 Cs u0 {2,S} {4,S} {5,S} {6,S}
2    S2s u0 {1,S} {3,S}
3 *2 Cs u0 {2,S} {7,S} {8,S} {9,S}
4    Cs u0 {1,S}
5    Cs u0 {1,S}
6    Cs u0 {1,S}
7    Cs u0 {3,S}
8    Cs u0 {3,S}
9    Cs u0 {3,S}
""",
    solute=SoluteData(
        S=0.01584,
        B=0.01615,
        E=0.00654,
        L=-0.01076,
        A=0.00000,
    ),
    shortDesc="""""",
    longDesc="""
Number of data used to fit each solute parameter:
S: 1
B: 1
E: 1
L: 1
A: 1
""",
)

tree(
    """
L1: R
	L2: int14_gauche
		L3: CsCs
			L4: CsCs-P
			L4: CsCs-S
				L5: CsCs-SS
				L5: CsCs-ST
				L5: CsCs-SQ
			L4: CsCs-T
				L5: CsCs-TT
					L6: CsCs-T(TTP)
					L6: CsCs-T(TTS)
				L5: CsCs-TQ
			L4: CsCs-Q
				L5: CsCs-QQ
		L3: OsCs
			L4: OsCs-P
			L4: OsCs-S
				L5: OsCs-SP
				L5: OsCs-SS
				L5: OsCs-ST
				L5: OsCs-SQ
		L3: CdCs
			L4: CdCs-P
			L4: CdCs-S
				L5: CdCs-SP
				L5: CdCs-SS
				L5: CdCs-ST
				L5: CdCs-SQ
	L2: int15
		L3: CsCsCs
			L4: CsCsCs-TQ
			L4: CsCsCs-QQ
		L3: CsOsCs
			L4: CsOsCs-TQ
			L4: CsOsCs-QQ
		L3: CsSsCs
			L4: CsSsCs-QQ
"""
)
