tree(
"""
L1: CsOs
    L2: Cs(Cs(CsCsCs)Cs(CsCsR)RR)
        L3: Cs(Cs(CsCsCs)Cs(CsCsCs)RR)
    L2: Os(Cs(CsCsCs)Cs(CsCsR))
        L3: Os(Cs(CsCsCs)Cs(CsCsCs))
"""
)

thermo(
    label="CsOs",
    group=
        """
        1  *  {Cs,Os} 0 {2,S} {3,S}
        2     Cs 0 {1,S} {4,S} {6,S} {5,S}
        3     Cs 0 {1,S} {7,S} {8,S}
        4     Cs 0 {2,S}
        5     Cs 0 {2,S}
        6     Cs 0 {2,S}
        7     Cs 0 {3,S}
        8     Cs 0 {3,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(0,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=0,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs(Cs(CsCsCs)Cs(CsCsR)RR)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S}
        2     Cs 0 {1,S} {4,S} {6,S} {5,S}
        3     Cs 0 {1,S} {7,S} {8,S}
        4     Cs 0 {2,S}
        5     Cs 0 {2,S}
        6     Cs 0 {2,S}
        7     Cs 0 {3,S}
        8     Cs 0 {3,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(1.5,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=1,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs(Cs(CsCsCs)Cs(CsCsCs)RR)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S}
        2     Cs 0 {1,S} {4,S} {6,S} {5,S}
        3     Cs 0 {1,S} {7,S} {8,S} {9,S}
        4     Cs 0 {2,S}
        5     Cs 0 {2,S}
        6     Cs 0 {2,S}
        7     Cs 0 {3,S}
        8     Cs 0 {3,S}
        9     Cs 0 {3,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(3,"kcal/mol"),
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
    label="Os(Cs(CsCsCs)Cs(CsCsR))",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Cs 0 {1,S} {4,S} {5,S} {6,S}
        3     Cs 0 {1,S} {7,S} {8,S}
        4     Cs 0 {2,S}
        5     Cs 0 {2,S}
        6     Cs 0 {2,S}
        7     Cs 0 {3,S}
        8     Cs 0 {3,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(3.5,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=3,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os(Cs(CsCsCs)Cs(CsCsCs))",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Cs 0 {1,S} {4,S} {5,S} {6,S}
        3     Cs 0 {1,S} {7,S} {8,S} {9,S}
        4     Cs 0 {2,S}
        5     Cs 0 {2,S}
        6     Cs 0 {2,S}
        7     Cs 0 {3,S}
        8     Cs 0 {3,S}
        9     Cs 0 {3,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(7,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=4,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

