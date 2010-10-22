tree(
"""
L1: R
    L2: ketene
        L3: ketene_2C-C
        L3: ketene_1C-C_1C-H
        L3: biketene
        L3: ketene_2C-H
    L2: cis
        L3: 2-ene_cis
            L4: 2-butene_cis
            L4: t-butyl_cis_2-ene
        L3: higher-ene_cis
            L4: t-butyl_cis
                L5: t-butyl_cis_t-butyl
    L2: double_cis
    L2: ortho
"""
)

thermo(
    label="R",
    group="OR{ketene,cis,double_cis,ortho}",
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(0,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=0,
    short_comment="dummy root",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="ketene",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     C 0 {1,D} {5,D}
        3     R 0 {1,S}
        4     R 0 {1,S}
        5     O 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(0,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=10,
    short_comment="All the corrections from this family are from Sumathi & Green, J. Phys. Chem. A, 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="ketene_2C-C",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     C 0 {1,D} {5,D}
        3     {Cs,Cd} 0 {1,S} {6,S}
        4     {Cs,Cd} 0 {1,S} {7,S}
        5     O 0 {2,D}
        6     C 0 {3,S}
        7     C 0 {4,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(-1.6,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=13,
    short_comment="This is correction NN2 from Sumathi & Green, J. Phys. Chem. A, 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="ketene_1C-C_1C-H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     C 0 {1,D} {5,D}
        3     {Cs,Cd} 0 {1,S} {6,S}
        4     C 0 {1,S} {7,S} {8,S} {9,S}
        5     O 0 {2,D}
        6     C 0 {3,S}
        7     H 0 {4,S}
        8     H 0 {4,S}
        9     H 0 {4,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(-0.5,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=11,
    short_comment="This is correction NN1 from Sumathi & Green, J. Phys. Chem. A, 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="biketene",
    group=
        """
        1     C 0 {2,S} {3,S} {4,S} {5,S}
        2     C 0 {1,S} {6,D}
        3  *  C 0 {1,S} {7,D} {10,S}
        4     R!H 0 {1,S}
        5     R!H 0 {1,S}
        6     C 0 {2,D} {8,D}
        7     C 0 {3,D} {9,D}
        8     O 0 {6,D}
        9     O 0 {7,D}
        10    R 0 {3,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(-0.9,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=14,
    short_comment="This is correction NN3 from Sumathi & Green, J. Phys. Chem. A, 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="ketene_2C-H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     C 0 {1,D} {5,D}
        3     C 0 {1,S} {6,S} {7,S} {8,S}
        4     C 0 {1,S} {9,S} {11,S} {10,S}
        5     O 0 {2,D}
        6     H 0 {3,S}
        7     H 0 {3,S}
        8     H 0 {3,S}
        9     H 0 {4,S}
        10    H 0 {4,S}
        11    H 0 {4,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(0,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=12,
    short_comment="This is correction NN0 from Sumathi & Green, J. Phys. Chem. A, 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="cis",
    group=
        """
        1  *  C 0 {2,Dcis} {3,S} {4,S}
        2  *  C 0 {1,Dcis} {5,S} {6,S}
        3     R!H 0 {1,S}
        4     H 0 {1,S}
        5     R!H 0 {2,S}
        6     H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.34,-1.09,-0.81,-0.61,-0.39,-0.26,0],"cal/(mol*K)"),
        H298=(1,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=1,
    short_comment="Cis double bond interaction BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="2-ene_cis",
    group=
        """
        1  *  C 0 {2,Dcis} {3,S} {4,S}
        2  *  C 0 {1,Dcis} {5,S} {6,S}
        3     C 0 {1,S} {7,S} {8,S} {9,S}
        4     H 0 {1,S}
        5     R!H 0 {2,S}
        6     H 0 {2,S}
        7     H 0 {3,S}
        8     H 0 {3,S}
        9     H 0 {3,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.34,-1.09,-0.81,-0.61,-0.39,-0.26,0],"cal/(mol*K)"),
        H298=(1,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=2,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="2-butene_cis",
    group=
        """
        1  *  C 0 {2,Dcis} {3,S} {4,S}
        2  *  C 0 {1,Dcis} {5,S} {6,S}
        3     C 0 {1,S} {7,S} {8,S} {9,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {10,S} {12,S} {11,S}
        6     H 0 {2,S}
        7     H 0 {3,S}
        8     H 0 {3,S}
        9     H 0 {3,S}
        10    H 0 {5,S}
        11    H 0 {5,S}
        12    H 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.34,-1.09,-0.81,-0.61,-0.39,-0.26,0],"cal/(mol*K)"),
        H298=(1,"kcal/mol"),
        S298=(1.2,"cal/(mol*K)"),
    ),
    index=3,
    short_comment="The entropy correction for 2-cis-butene is not zero BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="t-butyl_cis_2-ene",
    group=
        """
        1  *  C 0 {2,Dcis} {3,S} {4,S}
        2  *  C 0 {1,Dcis} {5,S} {6,S}
        3     C 0 {1,S} {7,S} {8,S} {9,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {10,S} {11,S} {12,S}
        6     H 0 {2,S}
        7     H 0 {3,S}
        8     H 0 {3,S}
        9     H 0 {3,S}
        10    R!H 0 {5,S}
        11    R!H 0 {5,S}
        12    R!H 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.34,-1.09,-0.81,-0.61,-0.39,-0.26,0],"cal/(mol*K)"),
        H298=(4,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=4,
    short_comment="Cis double bond interaction involving tertiary butyl group BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="higher-ene_cis",
    group=
        """
        1  *  C 0 {2,Dcis} {3,S} {4,S}
        2  *  C 0 {1,Dcis} {5,S} {6,S}
        3     C 0 {1,S} {7,S} {9,S} {8,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {10,S} {11,S} {12,S}
        6     H 0 {2,S}
        7     R!H 0 {3,S}
        8     R 0 {3,S}
        9     R 0 {3,S}
        10    R!H 0 {5,S}
        11    R 0 {5,S}
        12    R 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.34,-1.09,-0.81,-0.61,-0.39,-0.26,0],"cal/(mol*K)"),
        H298=(1,"kcal/mol"),
        S298=(-0.6,"cal/(mol*K)"),
    ),
    index=5,
    short_comment="The entropy correction for 2-cis-butene is not zero BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="t-butyl_cis",
    group=
        """
        1  *  C 0 {2,Dcis} {3,S} {4,S}
        2  *  C 0 {1,Dcis} {5,S} {6,S}
        3     C 0 {1,S} {7,S} {8,S} {9,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {10,S} {12,S} {11,S}
        6     H 0 {2,S}
        7     R!H 0 {3,S}
        8     R!H 0 {3,S}
        9     R!H 0 {3,S}
        10    R!H 0 {5,S}
        11    R 0 {5,S}
        12    R 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.34,-1.09,-0.81,-0.61,-0.39,-0.26,0],"cal/(mol*K)"),
        H298=(4,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=7,
    short_comment="Cis double bond interaction involving tertiary butyl group BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="t-butyl_cis_t-butyl",
    group=
        """
        1  *  C 0 {2,Dcis} {3,S} {4,S}
        2  *  C 0 {1,Dcis} {5,S} {6,S}
        3     C 0 {1,S} {7,S} {8,S} {9,S}
        4     H 0 {1,S}
        5     C 0 {2,S} {10,S} {11,S} {12,S}
        6     H 0 {2,S}
        7     R!H 0 {3,S}
        8     R!H 0 {3,S}
        9     R!H 0 {3,S}
        10    R!H 0 {5,S}
        11    R!H 0 {5,S}
        12    R!H 0 {5,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.34,-1.09,-0.81,-0.61,-0.39,-0.26,0],"cal/(mol*K)"),
        H298=(10,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=6,
    short_comment="Cis double bond interaction invloving two tertiary butyl groups BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="double_cis",
    group=
        """
        1  *  C 0 {2,Dcis} {3,S} {4,S}
        2  *  C 0 {1,Dcis} {5,S} {6,S}
        3     R!H 0 {1,S}
        4     R!H 0 {1,S}
        5     R!H 0 {2,S}
        6     R!H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([-1.34,-1.09,-0.81,-0.61,-0.39,-0.26,0],"cal/(mol*K)"),
        H298=(3,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=8,
    short_comment="2 Cis interactions around a double bond BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="ortho",
    group=
        """
        1  *  C 0 {2,B} {3,B} {4,S}
        2  *  C 0 {1,B} {5,B} {6,S}
        3     C 0 {1,B}
        4     R!H 0 {1,S}
        5     C 0 {2,B}
        6     R!H 0 {2,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([1.12,1.35,1.3,1.17,0.88,0.66,-0.05],"cal/(mol*K)"),
        H298=(0.57,"kcal/mol"),
        S298=(-1.61,"cal/(mol*K)"),
    ),
    index=9,
    short_comment="Ortho correction from BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

