tree(
"""
L1: Radical
    L2: RJ
        L3: CJ
            L4: CsJ
                L5: CH3
                L5: Cs_P
                    L6: CsCsJ
                        L7: CJCOOH
                        L7: CCJ
                        L7: RCCJ
                        L7: Isobutyl
                        L7: Neopentyl
                    L6: Benzyl_P
                    L6: Allyl_P
                        L7: C=CC=CCJ
                        L7: CTCC=CCJ
                        L7: C2JC=O
                    L6: Propargyl
                L5: Cs_S
                    L6: (Cs)2CsJ
                        L7: CCJCOOH
                        L7: CCJC
                        L7: RCCJC
                        L7: RCCJCC
                    L6: Benzyl_S
                    L6: Allyl_S
                        L7: CCJCHO
                    L6: C=CCJC=C
                    L6: Sec_Propargyl
                L5: Cs_T
                    L6: Tertalkyl
                        L7: C2CJCOOH
                    L6: Benzyl_T
                    L6: Allyl_T
                        L7: C2CJCHO
                    L6: Tert_Propargyl
                L5: CsJO
                    L6: CsJOH
                    L6: CsJOC
                        L7: CsJOCs
                            L8: CsJOCH3
                            L8: CsJOCC
                            L8: CsJOCC2
                            L8: CsJOCC3
                        L7: CsJOCds
                            L8: CsJOC(O)
                                L9: CsJOC(O)H
                                L9: CsJOC(O)C
                    L6: CsJOO
                        L7: CsJOOH
                        L7: CsJOOC
                L5: CCsJO
                    L6: CCsJOH
                    L6: CCsJOC
                        L7: CCsJOCs
                        L7: CCsJOCds
                            L8: CCsJOC(O)
                                L9: CCsJOC(O)H
                                L9: CCsJOC(O)C
                    L6: CCsJOO
                        L7: CCsJOOH
                        L7: CCsJOOC
                L5: C2CsJO
                    L6: C2CsJOH
                    L6: C2CsJOC
                        L7: C2CsJOCs
                        L7: C2CsJOCds
                            L8: C2CsJOC(O)
                                L9: C2CsJOC(O)H
                                L9: C2CsJOC(O)C
                    L6: C2CsJOO
                        L7: C2CsJOOH
                        L7: C2CsJOOC
            L4: CdsJ
                L5: CdsJO
                    L6: HCdsJO
                    L6: CCJ=O
                        L7: CsCJ=O
                        L7: C=CCJ=O
                    L6: (O)CJO
                        L7: (O)CJOH
                        L7: (O)CJOC
                            L8: (O)CJOCH3
                            L8: (O)CJOCC
                            L8: (O)CJOCC2
                            L8: (O)CJOCC3
                L5: Cds_P
                    L6: C=C=CJ
                L5: Cds_S
                    L6: C=CJC=C
            L4: CtJ
                L5: Acetyl
            L4: CbJ
        L3: OJ
            L4: HOJ
            L4: COJ
                L5: CsOJ
                    L6: H3COJ
                L5: CdsOJ
                    L6: RC=COJ
                    L6: OJC=O
            L4: OOJ
                L5: ROOJ
                    L6: C(=O)OOJ
                    L6: C3COOJ
                L5: HOOJ
        L3: SiJ
        L3: SJ
    L2: RJ2
        L3: CJ2
            L4: CsJ2
                L5: CH2
                    L6: CH2_t
                    L6: CH2_s
                L5: CsJ2_P
                    L6: CsCsJ2
                        L7: CCJ2
                            L8: CCJ2_t
                            L8: CCJ2_s
                    L6: PhCH
                        L7: PhCH_t
                        L7: PhCH_s
                    L6: AllylJ2
                        L7: AllylJ2_t
                        L7: AllylJ2_s
                L5: CsJ2_S
            L4: CdJ2
                L5: CCdJ2
                    L6: CCdJ2_t
                    L6: CCdJ2_s
                L5: CO
        L3: Oa
            L4: Oa_t
            L4: Oa_s
        L3: SiJ2
        L3: Sa
    L2: RJ3
        L3: CJ3
        L3: SiJ3
"""
)

thermo(
    label="Radical",
    group="OR{RJ, RJ2, RJ3}",
    node="RJ",
    index=0,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="RJ",
    group=
        """
        1  *  R 1
        """,
    node="CJ",
    index=1,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CJ",
    group=
        """
        1  *  C 1
        """,
    node="CsJ",
    index=2,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJ",
    group=
        """
        1  *  Cs 1
        """,
    node="Cs_P",
    index=3,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CH3",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     H 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.71,0.34,-0.33,-1.07,-2.43,-3.54,-5.43],"cal/(mol*K)"),
        H298=(104.81,"kcal/mol"),
        S298=(0.52,"cal/(mol*K)"),
    ),
    index=4,
    short_comment="Calculated in relation to methane from NIST values",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs_P",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     C 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.77,-1.36,-1.91,-2.4,-3.16,-3.74,-4.66],"cal/(mol*K)"),
        H298=(101.1,"kcal/mol"),
        S298=(2.61,"cal/(mol*K)"),
    ),
    index=5,
    short_comment="Generic primary radical. (CHEN & BOZZELLI) #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsCsJ",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cs 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cs_P",
    index=6,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CJCOOH",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cs 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     O 0 {2,S} {6,S}
        6     Os 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.25,-0.76,-1.34,-1.91,-2.87,-3.6,-4.69],"cal/(mol*K)"),
        H298=(103.26,"kcal/mol"),
        S298=(3.54,"cal/(mol*K)"),
    ),
    index=11,
    short_comment="WIJAYA et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCJ",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     C 0 {1,S} {5,S} {7,S} {6,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {2,S}
        6     H 0 {2,S}
        7     H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.65,-1.21,-1.75,-2.24,-3.02,-3.63,-3.63],"cal/(mol*K)"),
        H298=(101.1,"kcal/mol"),
        S298=(2.61,"cal/(mol*K)"),
    ),
    index=7,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="RCCJ",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     C 0 {1,S} {5,S} {6,S} {7,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S}
        6     H 0 {2,S}
        7     H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.77,-1.36,-1.91,-2.4,-3.16,-3.74,-4.66],"cal/(mol*K)"),
        H298=(101.1,"kcal/mol"),
        S298=(2.61,"cal/(mol*K)"),
    ),
    index=8,
    short_comment="LAY et al. CHEN & BOZZELLI #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Isobutyl",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     C 0 {1,S} {5,S} {6,S} {7,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S}
        6     C 0 {2,S}
        7     H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.54,-1.26,-1.92,-2.46,-3.27,-3.84,-3.84],"cal/(mol*K)"),
        H298=(101.1,"kcal/mol"),
        S298=(2.91,"cal/(mol*K)"),
    ),
    index=9,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Neopentyl",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     C 0 {1,S} {5,S} {7,S} {6,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S}
        6     C 0 {2,S}
        7     C 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.59,-1.32,-2.05,-2.65,-3.5,-4.06,-4.87],"cal/(mol*K)"),
        H298=(101.1,"kcal/mol"),
        S298=(3.03,"cal/(mol*K)"),
    ),
    index=10,
    short_comment="LAY et al. CHEN & BOZZELLI #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Benzyl_P",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cb 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.75,0.6,0.13,-0.42,-1.41,-2.18,-2.18],"cal/(mol*K)"),
        H298=(88.5,"kcal/mol"),
        S298=(-4.74,"cal/(mol*K)"),
    ),
    index=12,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Allyl_P",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cd 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.62,-0.56,-0.78,-1.12,-1.84,-2.46,-3.49],"cal/(mol*K)"),
        H298=(88.2,"kcal/mol"),
        S298=(-2.56,"cal/(mol*K)"),
    ),
    index=13,
    short_comment="LAY et al. CHEN & BOZZELLI #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C=CC=CCJ",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cd 0 {1,S} {5,D}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,D} {6,S}
        6     Cd 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.83,-1.86,-1.98,-1.99,-2.3,-2.5,-2.5],"cal/(mol*K)"),
        H298=(80,"kcal/mol"),
        S298=(-1.55,"cal/(mol*K)"),
    ),
    index=14,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CTCC=CCJ",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cd 0 {1,S} {5,D}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,D} {6,S}
        6     Ct 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.09,-1.62,-2.01,-2.63,-3.07,-3.48,-3.48],"cal/(mol*K)"),
        H298=(81,"kcal/mol"),
        S298=(-3.55,"cal/(mol*K)"),
    ),
    index=15,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C2JC=O",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cd 0 {1,S} {5,D} {6,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     O 0 {2,D}
        6     C 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.32,0.19,-0.15,-0.57,-1.43,-2.22,-3.67],"cal/(mol*K)"),
        H298=(94.4,"kcal/mol"),
        S298=(-1.16,"cal/(mol*K)"),
    ),
    index=16,
    short_comment="CHEN & BOZZELLI",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Propargyl",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Ct 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.84,-1.17,-1.56,-1.95,-2.7,-3.31,-5.31],"cal/(mol*K)"),
        H298=(89.4,"kcal/mol"),
        S298=(-0.51,"cal/(mol*K)"),
    ),
    index=17,
    short_comment="LAY et al. CHEN & BOZZELLI #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs_S",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     C 0 {1,S}
        3     C 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.5,-2.33,-3.1,-3.39,-3.75,-4.45,-5.2],"cal/(mol*K)"),
        H298=(98.45,"kcal/mol"),
        S298=(4.44,"cal/(mol*K)"),
    ),
    index=18,
    short_comment="Generic secondary radical. (CHEN & BOZZELLI) #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="(Cs)2CsJ",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cs 0 {1,S}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cs_S",
    index=19,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCJCOOH",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cs 0 {1,S} {5,S}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        5     O 0 {2,S} {6,S}
        6     O 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.65,-1.4,-2,-2.5,-3.27,-3.84,-4.73],"cal/(mol*K)"),
        H298=(99.98,"kcal/mol"),
        S298=(4.79,"cal/(mol*K)"),
    ),
    index=23,
    short_comment="WIJAYA et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCJC",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     C 0 {1,S} {5,S} {6,S} {7,S}
        3     C 0 {1,S} {8,S} {10,S} {9,S}
        4     H 0 {1,S}
        5     H 0 {2,S}
        6     H 0 {2,S}
        7     H 0 {2,S}
        8     H 0 {3,S}
        9     H 0 {3,S}
        10    H 0 {3,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.3,-2.36,-3.02,-3.44,-3.98,-4.36,-4.99],"cal/(mol*K)"),
        H298=(98.45,"kcal/mol"),
        S298=(4.51,"cal/(mol*K)"),
    ),
    index=20,
    short_comment="LAY et al. CHEN & BOZZELLI #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="RCCJC",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     C 0 {1,S} {5,S} {7,S} {6,S}
        3     C 0 {1,S} {8,S} {9,S} {10,S}
        4     H 0 {1,S}
        5     C 0 {2,S}
        6     H 0 {2,S}
        7     H 0 {2,S}
        8     H 0 {3,S}
        9     H 0 {3,S}
        10    H 0 {3,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.54,-2.77,-3.49,-3.9,-4.35,-4.64,-4.64],"cal/(mol*K)"),
        H298=(98.45,"kcal/mol"),
        S298=(5.13,"cal/(mol*K)"),
    ),
    index=21,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="RCCJCC",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     C 0 {1,S} {5,S} {7,S} {6,S}
        3     C 0 {1,S} {8,S} {9,S} {10,S}
        4     H 0 {1,S}
        5     C 0 {2,S}
        6     H 0 {2,S}
        7     H 0 {2,S}
        8     C 0 {3,S}
        9     H 0 {3,S}
        10    H 0 {3,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.71,-3.14,-3.92,-4.33,-4.71,-4.92,-4.92],"cal/(mol*K)"),
        H298=(98.45,"kcal/mol"),
        S298=(4.9,"cal/(mol*K)"),
    ),
    index=22,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Benzyl_S",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cb 0 {1,S}
        3     C 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.87,0.09,-0.63,-1.21,-2.07,-2.69,-2.69],"cal/(mol*K)"),
        H298=(85.9,"kcal/mol"),
        S298=(-5.04,"cal/(mol*K)"),
    ),
    index=24,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Allyl_S",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cd 0 {1,S}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.54,-1.82,-2.08,-2.32,-2.75,-3.14,-3.85],"cal/(mol*K)"),
        H298=(85.6,"kcal/mol"),
        S298=(-3.81,"cal/(mol*K)"),
    ),
    index=25,
    short_comment="LAY et al. CHEN & BOZZELLI #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCJCHO",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cd 0 {1,S} {5,D} {6,S}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        5     O 0 {2,D}
        6     H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.36,-1.57,-1.73,-1.89,-2.66,-3.11,-3.5],"cal/(mol*K)"),
        H298=(91.9,"kcal/mol"),
        S298=(-2.37,"cal/(mol*K)"),
    ),
    index=26,
    short_comment="CHEN & BOZZELLI #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C=CCJC=C",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-2.13,-1.96,-1.88,-1.89,-2.2,-2.6,-2.6],"cal/(mol*K)"),
        H298=(76,"kcal/mol"),
        S298=(-4.05,"cal/(mol*K)"),
    ),
    index=27,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Sec_Propargyl",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Ct 0 {1,S}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.59,-1.2,-1.75,-2.19,-2.91,-3.49,-3.49],"cal/(mol*K)"),
        H298=(87,"kcal/mol"),
        S298=(-0.45,"cal/(mol*K)"),
    ),
    index=28,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs_T",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     C 0 {1,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        """,
    node="Tertalkyl",
    index=29,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Tertalkyl",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cs 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.78,-2.48,-3.55,-4.15,-4.75,-5.02,-5.39],"cal/(mol*K)"),
        H298=(96.5,"kcal/mol"),
        S298=(5.24,"cal/(mol*K)"),
    ),
    index=30,
    short_comment="LAY et al. CHEN & BOZZELLI #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C2CJCOOH",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cs 0 {1,S} {5,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     O 0 {2,S} {6,S}
        6     O 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-3.54,-4.16,-4.44,-4.58,-4.74,-4.88,-5.23],"cal/(mol*K)"),
        H298=(97.2,"kcal/mol"),
        S298=(7.31,"cal/(mol*K)"),
    ),
    index=31,
    short_comment="WIJAYA et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Benzyl_T",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cb 0 {1,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.27,-0.78,-1.54,-2.06,-2.74,-3.19,-3.19],"cal/(mol*K)"),
        H298=(83.8,"kcal/mol"),
        S298=(-5.34,"cal/(mol*K)"),
    ),
    index=32,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Allyl_T",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cd 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.79,-2.38,-2.74,-2.97,-3.28,-3.55,-3.55],"cal/(mol*K)"),
        H298=(83.4,"kcal/mol"),
        S298=(-3.69,"cal/(mol*K)"),
    ),
    index=33,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C2CJCHO",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Cd 0 {1,S} {5,D} {6,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     O 0 {2,D}
        6     H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.62,-0.2,-1.23,-1.82,-2.87,-3.47,-3.47],"cal/(mol*K)"),
        H298=(89.8,"kcal/mol"),
        S298=(-1.71,"cal/(mol*K)"),
    ),
    index=34,
    short_comment="CHEN & BOZZELLI #. Value for Cp1500 taken as equal to Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Tert_Propargyl",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     Ct 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.04,-1.01,-1.74,-2.41,-3.19,-3.65,-3.65],"cal/(mol*K)"),
        H298=(84.5,"kcal/mol"),
        S298=(1.48,"cal/(mol*K)"),
    ),
    index=35,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJO",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        """,
    node="CsJOH",
    index=36,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJOH",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.25,0.18,-0.26,-0.83,-1.95,-2.85,-4.22],"cal/(mol*K)"),
        H298=(96.51,"kcal/mol"),
        S298=(0.09,"cal/(mol*K)"),
    ),
    index=37,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJOC",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S}
        """,
    node="CsJOCs",
    index=38,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJOCs",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     Cs 0 {2,S}
        """,
    node="CsJOCH3",
    index=39,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJOCH3",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {6,S} {7,S} {8,S}
        6     H 0 {5,S}
        7     H 0 {5,S}
        8     H 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.16,-0.4,-0.82,-1.33,-2.32,-3.13,-4.37],"cal/(mol*K)"),
        H298=(97,"kcal/mol"),
        S298=(0.78,"cal/(mol*K)"),
    ),
    index=40,
    short_comment="SUMATHI & GREEN #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJOCC",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {6,S} {7,S} {8,S}
        6     C 0 {5,S}
        7     H 0 {5,S}
        8     H 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.01,-1.22,-1.4,-1.71,-3.5,-3.24,-4.42],"cal/(mol*K)"),
        H298=(96.83,"kcal/mol"),
        S298=(1.41,"cal/(mol*K)"),
    ),
    index=41,
    short_comment="Calculated from data in SUMATHI & GREEN. Values might have large error bars.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJOCC2",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {6,S} {8,S} {7,S}
        6     C 0 {5,S}
        7     C 0 {5,S}
        8     H 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.95,0.75,0.23,-0.43,-1.71,-2.72,-4.19],"cal/(mol*K)"),
        H298=(96.16,"kcal/mol"),
        S298=(-0.59,"cal/(mol*K)"),
    ),
    index=42,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJOCC3",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {6,S} {8,S} {7,S}
        6     C 0 {5,S}
        7     C 0 {5,S}
        8     C 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.08,-0.09,-0.52,-1.06,-2.11,-2.96,-4.27],"cal/(mol*K)"),
        H298=(95.75,"kcal/mol"),
        S298=(0.27,"cal/(mol*K)"),
    ),
    index=43,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJOCds",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     {Cd,CO} 0 {2,S}
        """,
    node="CsJOC(O)",
    index=44,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJOC(O)",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {6,D}
        6     O 0 {5,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.91,0.89,0.42,-0.21,-1.5,-2.62,-4.43],"cal/(mol*K)"),
        H298=(100.7,"kcal/mol"),
        S298=(-0.18,"cal/(mol*K)"),
    ),
    index=45,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJOC(O)H",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {6,D} {7,S}
        6     O 0 {5,D}
        7     H 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.95,0.97,0.53,-0.12,-1.54,-2.76,-4.53],"cal/(mol*K)"),
        H298=(100.88,"kcal/mol"),
        S298=(-0.18,"cal/(mol*K)"),
    ),
    index=46,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJOC(O)C",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {6,D} {7,S}
        6     O 0 {5,D}
        7     C 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.88,0.81,0.31,-0.3,-1.45,-2.47,-4.33],"cal/(mol*K)"),
        H298=(100.48,"kcal/mol"),
        S298=(-0.17,"cal/(mol*K)"),
    ),
    index=47,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJOO",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     O 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.18,-0.42,-0.79,-1.2,-1.99,-2.63,-3.65],"cal/(mol*K)"),
        H298=(98.5,"kcal/mol"),
        S298=(-1.57,"cal/(mol*K)"),
    ),
    index=48,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJOOH",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     O 0 {2,S} {6,S}
        6     H 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.06,-0.35,-0.76,-1.19,-1.99,-2.64,-3.68],"cal/(mol*K)"),
        H298=(98.91,"kcal/mol"),
        S298=(-1.52,"cal/(mol*K)"),
    ),
    index=49,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJOOC",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     O 0 {2,S} {6,S}
        6     C 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.31,-0.48,-0.82,-1.22,-1.99,-2.62,-3.63],"cal/(mol*K)"),
        H298=(98.34,"kcal/mol"),
        S298=(-1.62,"cal/(mol*K)"),
    ),
    index=50,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCsJO",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S}
        3     C 0 {1,S}
        4     H 0 {1,S}
        """,
    node="CCsJOC",
    index=51,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCsJOH",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.65,-0.01,-0.75,-1.43,-2.52,-3.31,-4.47],"cal/(mol*K)"),
        H298=(95.39,"kcal/mol"),
        S298=(0.92,"cal/(mol*K)"),
    ),
    index=52,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCsJOC",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S}
        """,
    node="CCsJOCs",
    index=53,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCsJOCs",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     H 0 {1,S}
        5     Cs 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.82,0.53,-0.11,-0.86,-2.2,-3.18,-4.51],"cal/(mol*K)"),
        H298=(95.41,"kcal/mol"),
        S298=(0.33,"cal/(mol*K)"),
    ),
    index=54,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCsJOCds",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     H 0 {1,S}
        5     {CO,Cd} 0 {2,S}
        """,
    node="CCsJOC(O)",
    index=55,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCsJOC(O)",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {6,D}
        6     O 0 {5,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([1.16,0.78,0.05,-0.73,-2.13,-3.24,-4.9],"cal/(mol*K)"),
        H298=(98.7,"kcal/mol"),
        S298=(0.98,"cal/(mol*K)"),
    ),
    index=56,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCsJOC(O)H",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {6,D} {7,S}
        6     O 0 {5,D}
        7     H 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([1.2,0.88,0.16,-0.67,-2.22,-3.43,-5],"cal/(mol*K)"),
        H298=(98.87,"kcal/mol"),
        S298=(0.98,"cal/(mol*K)"),
    ),
    index=57,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCsJOC(O)C",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {6,D} {7,S}
        6     O 0 {5,D}
        7     C 0 {5,S}
        """,
    node="CCsJOC(O)",
    index=0,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCsJOO",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     H 0 {1,S}
        5     O 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.48,-1.15,-1.68,-2.11,-2.77,-3.26,-4.02],"cal/(mol*K)"),
        H298=(96.9,"kcal/mol"),
        S298=(0.76,"cal/(mol*K)"),
    ),
    index=58,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCsJOOH",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     H 0 {1,S}
        5     O 0 {2,S} {6,S}
        6     H 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.39,-1.08,-1.64,-2.08,-2.75,-3.26,-4.03],"cal/(mol*K)"),
        H298=(97.19,"kcal/mol"),
        S298=(0.77,"cal/(mol*K)"),
    ),
    index=59,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCsJOOC",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     H 0 {1,S}
        5     O 0 {2,S} {6,S}
        6     C 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.58,-1.21,-1.73,-2.15,-2.8,-3.27,-4.01],"cal/(mol*K)"),
        H298=(96.64,"kcal/mol"),
        S298=(0.74,"cal/(mol*K)"),
    ),
    index=60,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C2CsJO",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        """,
    node="C2CsJOC",
    index=61,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C2CsJOH",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        5     H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.31,-0.66,-1.54,-2.23,-3.17,-3.8,-4.72],"cal/(mol*K)"),
        H298=(94.5,"kcal/mol"),
        S298=(2.17,"cal/(mol*K)"),
    ),
    index=62,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C2CsJOC",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        5     C 0 {2,S}
        """,
    node="C2CsJOCs",
    index=63,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C2CsJOCs",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        5     Cs 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.09,-1.37,-2.49,-3.26,-4.15,-4.63,-5.23],"cal/(mol*K)"),
        H298=(95.5,"kcal/mol"),
        S298=(3.71,"cal/(mol*K)"),
    ),
    index=64,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C2CsJOCds",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        5     {Cd,CO} 0 {2,S}
        """,
    node="C2CsJOC(O)",
    index=65,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C2CsJOC(O)",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        5     C 0 {2,S} {6,D}
        6     O 0 {5,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.04,-1.34,-2.3,-2.99,-3.99,-4.77,-5.98],"cal/(mol*K)"),
        H298=(100.1,"kcal/mol"),
        S298=(4.77,"cal/(mol*K)"),
    ),
    index=66,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C2CsJOC(O)H",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        5     C 0 {2,S} {6,D} {7,S}
        6     O 0 {5,D}
        7     H 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.03,-1.28,-2.28,-3.1,-4.35,-5.19,-6.06],"cal/(mol*K)"),
        H298=(99.97,"kcal/mol"),
        S298=(4.88,"cal/(mol*K)"),
    ),
    index=67,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C2CsJOC(O)C",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        5     C 0 {2,S} {6,D} {7,S}
        6     O 0 {5,D}
        7     C 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.04,-1.4,-2.32,-2.89,-3.62,-4.36,-5.9],"cal/(mol*K)"),
        H298=(100.25,"kcal/mol"),
        S298=(4.66,"cal/(mol*K)"),
    ),
    index=68,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C2CsJOO",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        5     O 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.89,-2.09,-2.81,-3.24,-3.69,-3.97,-4.43],"cal/(mol*K)"),
        H298=(96.7,"kcal/mol"),
        S298=(2.22,"cal/(mol*K)"),
    ),
    index=69,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C2CsJOOH",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        5     O 0 {2,S} {6,S}
        6     H 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.01,-2.17,-2.87,-3.3,-3.77,-4.05,-4.49],"cal/(mol*K)"),
        H298=(96.74,"kcal/mol"),
        S298=(2.37,"cal/(mol*K)"),
    ),
    index=70,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C2CsJOOC",
    group=
        """
        1  *  C 1 {2,S} {3,S} {4,S}
        2     O 0 {1,S} {5,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        5     O 0 {2,S} {6,S}
        6     C 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.78,-2.02,-2.75,-3.18,-3.62,-3.88,-4.37],"cal/(mol*K)"),
        H298=(96.58,"kcal/mol"),
        S298=(2.08,"cal/(mol*K)"),
    ),
    index=71,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CdsJ",
    group=
        """
        1  *  {Cd,CO} 1
        """,
    node="Cds_P",
    index=72,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CdsJO",
    group=
        """
        1  *  C 1 {2,D}
        2     O 0 {1,D}
        """,
    node="CCJ=O",
    index=79,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="HCdsJO",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     O 0 {1,D}
        3     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.19,-0.65,-1.19,-1.73,-2.63,-3.32,-4.42],"cal/(mol*K)"),
        H298=(88.45,"kcal/mol"),
        S298=(-0.01,"cal/(mol*K)"),
    ),
    index=80,
    short_comment="Calculated in relation to formaldehyde from NIST values",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCJ=O",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     O 0 {1,D}
        3     C 0 {1,S}
        """,
    node="CsCJ=O",
    index=81,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsCJ=O",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     O 0 {1,D}
        3     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.83,-1.43,-1.96,-2.42,-3.16,-3.73,-4.64],"cal/(mol*K)"),
        H298=(89,"kcal/mol"),
        S298=(1.12,"cal/(mol*K)"),
    ),
    index=82,
    short_comment="CHEN & BOZZELLI #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C=CCJ=O",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     O 0 {1,D}
        3     Cd 0 {1,S} {4,D}
        4     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.19,-0.85,-1.59,-2.21,-3.21,-3.89,-4.61],"cal/(mol*K)"),
        H298=(83,"kcal/mol"),
        S298=(-1.39,"cal/(mol*K)"),
    ),
    index=83,
    short_comment="CHEN & BOZZELLI #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="(O)CJO",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     O 0 {1,D}
        3     O 0 {1,S}
        """,
    node="(O)CJOC",
    index=84,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="(O)CJOH",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     O 0 {1,D}
        3     O 0 {1,S} {4,S}
        4     H 0 {3,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.02,-0.66,-1.4,-2.12,-3.41,-4.44,-5.79],"cal/(mol*K)"),
        H298=(100.75,"kcal/mol"),
        S298=(0.78,"cal/(mol*K)"),
    ),
    index=85,
    short_comment="SUMATHI & GREEN #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="(O)CJOC",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     O 0 {1,D}
        3     O 0 {1,S} {4,S}
        4     C 0 {3,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.45,-0.27,-1.19,-2.1,-3.63,-4.69,-5.8],"cal/(mol*K)"),
        H298=(98.99,"kcal/mol"),
        S298=(0.72,"cal/(mol*K)"),
    ),
    index=86,
    short_comment="SUMATHI & GREEN (Hf assigned value of (O)CJOCH(CH3)2)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="(O)CJOCH3",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     O 0 {1,D}
        3     O 0 {1,S} {4,S}
        4     C 0 {3,S} {5,S} {7,S} {6,S}
        5     H 0 {4,S}
        6     H 0 {4,S}
        7     H 0 {4,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.51,-0.11,-0.94,-1.8,-3.34,-4.48,-5.79],"cal/(mol*K)"),
        H298=(100.1,"kcal/mol"),
        S298=(0.72,"cal/(mol*K)"),
    ),
    index=87,
    short_comment="SUMATHI & GREEN",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="(O)CJOCC",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     O 0 {1,D}
        3     O 0 {1,S} {4,S}
        4     C 0 {3,S} {5,S} {6,S} {7,S}
        5     C 0 {4,S}
        6     H 0 {4,S}
        7     H 0 {4,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.45,-0.13,-0.98,-1.86,-3.43,-4.56,-5.79],"cal/(mol*K)"),
        H298=(99.49,"kcal/mol"),
        S298=(0.55,"cal/(mol*K)"),
    ),
    index=88,
    short_comment="SUMATHI & GREEN (values from (O)CJOCH2CH3)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="(O)CJOCC2",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     O 0 {1,D}
        3     O 0 {1,S} {4,S}
        4     C 0 {3,S} {5,S} {7,S} {6,S}
        5     C 0 {4,S}
        6     C 0 {4,S}
        7     H 0 {4,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.74,-0.06,-1.04,-2.01,-3.6,-4.66,-5.77],"cal/(mol*K)"),
        H298=(98.99,"kcal/mol"),
        S298=(0.82,"cal/(mol*K)"),
    ),
    index=89,
    short_comment="SUMATHI & GREEN (values from (O)CJOCH(CH3)2)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="(O)CJOCC3",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     O 0 {1,D}
        3     O 0 {1,S} {4,S}
        4     C 0 {3,S} {5,S} {7,S} {6,S}
        5     C 0 {4,S}
        6     C 0 {4,S}
        7     C 0 {4,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.11,-0.79,-1.8,-2.73,-4.17,-5.06,-5.87],"cal/(mol*K)"),
        H298=(97.98,"kcal/mol"),
        S298=(0.76,"cal/(mol*K)"),
    ),
    index=90,
    short_comment="SUMATHI & GREEN (values from (O)CJOC(CH3)3)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds_P",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     C 0 {1,D}
        3     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.19,-0.75,-1.36,-1.92,-2.82,-3.49,-4.53],"cal/(mol*K)"),
        H298=(111.2,"kcal/mol"),
        S298=(1.39,"cal/(mol*K)"),
    ),
    index=74,
    short_comment="LAY et al. CHEN & BOZZELLI #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C=C=CJ",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     C 0 {1,D} {4,D}
        3     H 0 {1,S}
        4     C 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.45,-1.05,-1.64,-2.15,-2.98,-3.6,-3.6],"cal/(mol*K)"),
        H298=(89,"kcal/mol"),
        S298=(1.29,"cal/(mol*K)"),
    ),
    index=75,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds_S",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     C 0 {1,D}
        3     C 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.34,-1.21,-1.94,-2.52,-3.34,-3.91,-4.76],"cal/(mol*K)"),
        H298=(109,"kcal/mol"),
        S298=(1.81,"cal/(mol*K)"),
    ),
    index=77,
    short_comment="LAY et al. CHEN & BOZZELLI #",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C=CJC=C",
    group=
        """
        1  *  C 1 {2,D} {3,S}
        2     Cd 0 {1,D}
        3     {Cd,CO} 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0.19,-0.76,-1.51,-2.01,-2.7,-3.17,-3.17],"cal/(mol*K)"),
        H298=(99.8,"kcal/mol"),
        S298=(0.71,"cal/(mol*K)"),
    ),
    index=78,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CtJ",
    group=
        """
        1  *  C 1 {2,T}
        2     C 0 {1,T}
        """,
    node="Acetyl",
    index=91,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Acetyl",
    group=
        """
        1  *  C 1 {2,T}
        2     C 0 {1,T} {3,S}
        3     H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.51,-1.56,-2.27,-2.78,-3.47,-3.97,-3.97],"cal/(mol*K)"),
        H298=(132.7,"kcal/mol"),
        S298=(2.11,"cal/(mol*K)"),
    ),
    index=92,
    short_comment="LAY et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CbJ",
    group=
        """
        1  *  C 1 {2,B} {3,B}
        2     C 0 {1,B}
        3     C 0 {1,B}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.41,-1.18,-1.93,-2.69,-3.75,-4.48,-5.24],"cal/(mol*K)"),
        H298=(113,"kcal/mol"),
        S298=(1.48,"cal/(mol*K)"),
    ),
    index=93,
    short_comment="BDE from TSANG, S and Cp from THERM",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="OJ",
    group=
        """
        1  *  O 1
        """,
    node="COJ",
    index=94,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="HOJ",
    group=
        """
        1  *  O 1 {2,S}
        2     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.87,-1.1,-1.36,-1.62,-2.11,-2.53,-3.38],"cal/(mol*K)"),
        H298=(119.22,"kcal/mol"),
        S298=(-2.6,"cal/(mol*K)"),
    ),
    index=95,
    short_comment="Calculated from NIST values for H2O, OH and H",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="COJ",
    group=
        """
        1  *  O 1 {2,S}
        2     C 0 {1,S}
        """,
    node="CsOJ",
    index=135,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsOJ",
    group=
        """
        1  *  O 1 {2,S}
        2     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.98,-1.3,-1.61,-1.89,-2.38,-2.8,-3.59],"cal/(mol*K)"),
        H298=(104.06,"kcal/mol"),
        S298=(-1.46,"cal/(mol*K)"),
    ),
    index=96,
    short_comment="CHEN & BOZZELLI(ROJ)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="H3COJ",
    group=
        """
        1  *  O 1 {2,S}
        2     C 0 {1,S} {3,S} {5,S} {4,S}
        3     H 0 {2,S}
        4     H 0 {2,S}
        5     H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.11,-1.29,-1.62,-1.97,-2.59,-3.07,-3.84],"cal/(mol*K)"),
        H298=(104.27,"kcal/mol"),
        S298=(0.51,"cal/(mol*K)"),
    ),
    index=97,
    short_comment="Enthalpy HBI calculated from NIST values, entropy and Cp from B3LYP/6-31G* for CH3OH, CH3O and H",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CdsOJ",
    group=
        """
        1  *  O 1 {2,S}
        2     {Cd,CO} 0 {1,S}
        """,
    node="RC=COJ",
    index=98,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="RC=COJ",
    group=
        """
        1  *  O 1 {2,S}
        2     Cd 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.34,-1.99,-2.48,-2.79,-3.13,-3.33,-3.79],"cal/(mol*K)"),
        H298=(88,"kcal/mol"),
        S298=(-1.11,"cal/(mol*K)"),
    ),
    index=99,
    short_comment="CHEN & BOZZELLI",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="OJC=O",
    group=
        """
        1  *  O 1 {2,S}
        2     CO 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.31,-1.87,-2.32,-2.69,-3.28,-3.74,-4.56],"cal/(mol*K)"),
        H298=(104,"kcal/mol"),
        S298=(0.79,"cal/(mol*K)"),
    ),
    index=100,
    short_comment="CHEN & BOZZELLI",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="OOJ",
    group=
        """
        1  *  O 1 {2,S}
        2     O 0 {1,S}
        """,
    node="ROOJ",
    index=101,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="ROOJ",
    group=
        """
        1  *  O 1 {2,S}
        2     O 0 {1,S} {3,S}
        3     R!H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-2.05,-2.84,-3.55,-4.09,-4.72,-4.97,-5.08],"cal/(mol*K)"),
        H298=(88.2,"kcal/mol"),
        S298=(0.22,"cal/(mol*K)"),
    ),
    index=102,
    short_comment="CHEN & BOZZELLI",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C(=O)OOJ",
    group=
        """
        1  *  O 1 {2,S}
        2     O 0 {1,S} {3,S}
        3     C 0 {2,S} {4,D}
        4     O 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-2.05,-2.84,-3.55,-4.09,-4.72,-4.97,-5.08],"cal/(mol*K)"),
        H298=(98.33,"kcal/mol"),
        S298=(0.22,"cal/(mol*K)"),
    ),
    index=104,
    short_comment="HBI for enthalpy from CHEN & BOZZELLI. Cp and S values taken from ROOJ",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C3COOJ",
    group=
        """
        1  *  O 1 {2,S}
        2     O 0 {1,S} {3,S}
        3     C 0 {2,S} {4,S} {6,S} {5,S}
        4     C 0 {3,S}
        5     C 0 {3,S}
        6     C 0 {3,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-2.05,-2.84,-3.55,-4.09,-4.72,-4.97,-5.08],"cal/(mol*K)"),
        H298=(85.3,"kcal/mol"),
        S298=(0.22,"cal/(mol*K)"),
    ),
    index=103,
    short_comment="CHEN & BOZZELLI",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="HOOJ",
    group=
        """
        1  *  O 1 {2,S}
        2     O 0 {1,S} {3,S}
        3     H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.99,-2.68,-3.07,-3.3,-3.55,-3.66,-3.9],"cal/(mol*K)"),
        H298=(85.13,"kcal/mol"),
        S298=(-0.92,"cal/(mol*K)"),
    ),
    index=105,
    short_comment="Calculated from NIST values for H2O2, O2H and H",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="SiJ",
    group=
        """
        1  *  Si 1
        """,
    node="CJ",
    index=134,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="SJ",
    group=
        """
        1  *  S 1
        """,
    node="OJ",
    index=137,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="RJ2",
    group=
        """
        1  *  R {2S,2T}
        """,
    node="CJ2",
    index=106,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CJ2",
    group=
        """
        1  *  C {2S,2T}
        """,
    node="CsJ2",
    index=107,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJ2",
    group=
        """
        1  *  Cs {2S,2T}
        """,
    node="CH2",
    index=108,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CH2",
    group=
        """
        1  *  C {2S,2T} {2,S} {3,S}
        2     H 0 {1,S}
        3     H 0 {1,S}
        """,
    node="CH2_t",
    index=109,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CH2_t",
    group=
        """
        1  *  C 2T {2,S} {3,S}
        2     H 0 {1,S}
        3     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.27,-1.08,-2.14,-3.23,-5.18,-6.74,-9.47],"cal/(mol*K)"),
        H298=(214.44,"kcal/mol"),
        S298=(-1.73,"cal/(mol*K)"),
    ),
    index=110,
    short_comment="Calculated for methylene in relation to methane from NIST values",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CH2_s",
    group=
        """
        1  *  C 2S {2,S} {3,S}
        2     H 0 {1,S}
        3     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.27,-1.08,-2.14,-3.23,-5.18,-6.74,-9.47],"cal/(mol*K)"),
        H298=(223.7,"kcal/mol"),
        S298=(-1.73,"cal/(mol*K)"),
    ),
    index=111,
    short_comment="BDE JANOSCHEK & ROSSI. S and Cp from CH2_t.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJ2_P",
    group=
        """
        1  *  C {2S,2T} {2,S} {3,S}
        2     C 0 {1,S}
        3     H 0 {1,S}
        """,
    node="CsCsJ2",
    index=112,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsCsJ2",
    group=
        """
        1  *  C {2S,2T} {2,S} {3,S}
        2     Cs 0 {1,S}
        3     H 0 {1,S}
        """,
    node="CCJ2",
    index=113,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCJ2",
    group=
        """
        1  *  C {2S,2T} {2,S} {3,S}
        2     Cs 0 {1,S} {4,S} {6,S} {5,S}
        3     H 0 {1,S}
        4     H 0 {2,S}
        5     H 0 {2,S}
        6     H 0 {2,S}
        """,
    node="CCJ2_t",
    index=114,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCJ2_t",
    group=
        """
        1  *  C 2T {2,S} {3,S}
        2     Cs 0 {1,S} {4,S} {5,S} {6,S}
        3     H 0 {1,S}
        4     H 0 {2,S}
        5     H 0 {2,S}
        6     H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-0.81,-1.74,-2.69,-3.61,-5.18,-6.42,-8.36],"cal/(mol*K)"),
        H298=(211.3,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=115,
    short_comment="BDE and Cp calculated from data in KIM et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCJ2_s",
    group=
        """
        1  *  C 2S {2,S} {3,S}
        2     Cs 0 {1,S} {4,S} {5,S} {6,S}
        3     H 0 {1,S}
        4     H 0 {2,S}
        5     H 0 {2,S}
        6     H 0 {2,S}
        """,
    node="CCJ2_t",
    index=116,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="PhCH",
    group=
        """
        1  *  C {2S,2T} {2,S} {3,S}
        2     Cb 0 {1,S}
        3     H 0 {1,S}
        """,
    node="PhCH_t",
    index=117,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="PhCH_t",
    group=
        """
        1  *  C 2T {2,S} {3,S}
        2     Cb 0 {1,S}
        3     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(195,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=118,
    short_comment="BDE from PUTSMA et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="PhCH_s",
    group=
        """
        1  *  C 2S {2,S} {3,S}
        2     Cb 0 {1,S}
        3     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(205.8,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=119,
    short_comment="BDE from NGUYEN et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="AllylJ2",
    group=
        """
        1  *  C {2S,2T} {2,S} {3,S}
        2     Cd 0 {1,S}
        3     H 0 {1,S}
        """,
    node="AllylJ2_t",
    index=120,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="AllylJ2_t",
    group=
        """
        1  *  C 2T {2,S} {3,S}
        2     Cd 0 {1,S}
        3     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(192.8,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=121,
    short_comment="BDE from PUTSMA et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="AllylJ2_s",
    group=
        """
        1  *  C 2S {2,S} {3,S}
        2     Cd 0 {1,S}
        3     H 0 {1,S}
        """,
    node="AllylJ2_t",
    index=122,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CsJ2_S",
    group=
        """
        1  *  C {2S,2T} {2,S} {3,S}
        2     C 0 {1,S}
        3     C 0 {1,S}
        """,
    node="CsJ2_P",
    index=123,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CdJ2",
    group=
        """
        1  *  {Cd,CO} {2S,2T}
        """,
    node="CCdJ2",
    index=124,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCdJ2",
    group=
        """
        1  *  C {2S,2T} {2,D}
        2     C 0 {1,D}
        """,
    node="CCdJ2_s",
    index=125,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCdJ2_t",
    group=
        """
        1  *  C 2T {2,D}
        2     C 0 {1,D}
        """,
    node="CCdJ2_s",
    index=126,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CCdJ2_s",
    group=
        """
        1  *  C 2S {2,D}
        2     C 0 {1,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(190.7,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=127,
    short_comment="BDE from ERWIN et al.",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CO",
    group=
        """
        1  *  C {2S,2T} {2,D}
        2     O 0 {1,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.5,-2.38,-3.32,-4.24,-5.75,-6.88,-8.59],"cal/(mol*K)"),
        H298=(103.73,"kcal/mol"),
        S298=(-6.47,"cal/(mol*K)"),
    ),
    index=128,
    short_comment="Value for carbon monoxide calculated in relation to formaldehyde from NIST values",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Oa",
    group=
        """
        1  *  O {2S,2T}
        """,
    node="Oa_t",
    index=129,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Oa_t",
    group=
        """
        1  *  O 2T
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-2.8,-3.05,-3.33,-3.62,-4.24,-4.86,-6.28],"cal/(mol*K)"),
        H298=(221.55,"kcal/mol"),
        S298=(-8.02,"cal/(mol*K)"),
    ),
    index=130,
    short_comment="Calculated for atomic oxygen in relation to water from NIST values",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Oa_s",
    group=
        """
        1  *  O 2S
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-2.8,-3.05,-3.33,-3.62,-4.24,-4.86,-6.28],"cal/(mol*K)"),
        H298=(266.9,"kcal/mol"),
        S298=(-8.02,"cal/(mol*K)"),
    ),
    index=131,
    short_comment="BDE from SCHALLEY et al. S and Cp values taken from Oa_t",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="SiJ2",
    group=
        """
        1  *  Si {2S,2T}
        """,
    node="CJ2",
    index=135,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Sa",
    group=
        """
        1  *  S {2S,2T}
        """,
    node="Oa",
    index=138,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="RJ3",
    group=
        """
        1  *  R 3
        """,
    node="CJ3",
    index=132,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="CJ3",
    group=
        """
        1  *  C 3
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.57,-2.73,-4.11,-5.5,-7.92,-9.85,-12.95],"cal/(mol*K)"),
        H298=(316.19,"kcal/mol"),
        S298=(-5.7,"cal/(mol*K)"),
    ),
    index=133,
    short_comment="Calculated for methylidyene in relation to methane from NIST values",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="SiJ3",
    group=
        """
        1  *  Si 3
        """,
    node="CJ3",
    index=136,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

