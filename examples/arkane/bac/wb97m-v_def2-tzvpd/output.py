# Job 1: Petersson-type BACs:
# Training RMSE/MAE before fitting: 2.53/1.83 kcal/mol
# Training RMSE/MAE after fitting: 1.20/0.84 kcal/mol
"LevelOfTheory(method='wb97mv',basis='def2tzvpd',software='qchem')": {
    'Br-Br': -0.29900827104945105,
    'Br-C': 0.14519985383404127,
    'Br-Cl': 0.1564729358531758,
    'Br-F': -0.02198546881436414,
    'Br-H': -0.1877677960662112,
    'Br-O': -1.1705645162069829,
    'C#C': -1.4854470064856449,
    'C#N': 1.027764758058689,
    'C#O': 0.49218815065369625,
    'C-C': 0.5185443380878181,
    'C-Cl': 0.6293437113748294,
    'C-F': 0.5526648416071944,
    'C-H': -0.18369129390336966,
    'C-N': 1.5626088446555433,
    'C-O': 0.5539551755111399,
    'C-S': 0.7175710667178974,
    'C=C': 0.46679214572357686,
    'C=N': 1.17558540581856,
    'C=O': 0.4255962723777596,
    'C=S': 0.45344449204819365,
    'Cl-Cl': 0.7496140200477119,
    'Cl-F': 0.5089573574917736,
    'Cl-H': -0.258183629712196,
    'Cl-N': 1.4624406438788022,
    'Cl-O': 0.5419935093794896,
    'Cl-S': 0.24441683834083586,
    'F-F': -1.9508363888276188,
    'F-H': -2.5355256095650915,
    'F-O': -0.18223254017394275,
    'F-S': -0.06348672229322953,
    'H-H': -1.3532226208469873,
    'H-N': 0.20089687198226713,
    'H-O': -0.3553760324106493,
    'H-S': 1.0734603159144247,
    'N#N': 2.6437230556931492,
    'N-N': 3.508725977748223,
    'N-O': 1.5796389736297316,
    'N=N': 2.4303888779776615,
    'N=O': -0.8913561543694676,
    'O-O': -0.6278622651632636,
    'O-S': -0.6041291940722772,
    'O=O': -8.096244697146293,
    'O=S': -1.137823941569575,
    'S-S': 1.581099653809384,
    'S=S': -0.9106389050072522
},
# 95% Confidence interval half-widths:
# {
#     'Br-Br': 1.698244283517215,
#     'Br-C': 0.7432019263958128,
#     'Br-Cl': 1.698244283517215,
#     'Br-F': 1.698244283517215,
#     'Br-H': 1.698244283517215,
#     'Br-O': 1.1870300889782919,
#     'C#C': 0.7843168533379397,
#     'C#N': 0.6013681485762876,
#     'C#O': 1.698244283517215,
#     'C-C': 0.18148111644019385,
#     'C-Cl': 0.27633905592032315,
#     'C-F': 0.26627785118495734,
#     'C-H': 0.09288689324199018,
#     'C-N': 0.30910902586421124,
#     'C-O': 0.20412613303064972,
#     'C-S': 0.2592999409572966,
#     'C=C': 0.1927717926836692,
#     'C=N': 0.449800268575733,
#     'C=O': 0.3767727021101349,
#     'C=S': 0.6920909091163048,
#     'Cl-Cl': 1.698244283517215,
#     'Cl-F': 1.698244283517215,
#     'Cl-H': 1.698244283517215,
#     'Cl-N': 1.1939357137033106,
#     'Cl-O': 0.7556221678266412,
#     'Cl-S': 1.058024099815034,
#     'F-F': 1.698244283517215,
#     'F-H': 1.698244283517215,
#     'F-O': 0.6737783494519384,
#     'F-S': 0.34712459931271133,
#     'H-H': 1.698244283517215,
#     'H-N': 0.20361144006445098,
#     'H-O': 0.333538242333418,
#     'H-S': 0.8116191386638535,
#     'N#N': 1.180086389613317,
#     'N-N': 0.6910727280056931,
#     'N-O': 0.4516955485855899,
#     'N=N': 0.5910477502990007,
#     'N=O': 0.5410350156815463,
#     'O-O': 0.5715629648534439,
#     'O-S': 0.8461358448938238,
#     'O=O': 1.5018865435943265,
#     'O=S': 0.5680539924016892,
#     'S-S': 0.21023442820896468,
#     'S=S': 1.6321528063440098
# }

# Job 2: Melius-type BACs:
# Training RMSE/MAE before fitting: 2.53/1.83 kcal/mol
# Training RMSE/MAE after fitting: 1.34/0.88 kcal/mol
"LevelOfTheory(method='wb97mv',basis='def2tzvpd',software='qchem')": {
    'atom_corr': {
        'Br': -1.0069345550691,
        'C': -1.5084054119118733,
        'Cl': -1.1958375172653044,
        'F': -1.3797266986385623,
        'H': -0.35935274487563706,
        'N': -3.0907554240618076,
        'O': -1.7769694500519555,
        'S': -2.894625842421483
    },
    'bond_corr_length': {
        'Br': 1086.522209037005,
        'C': 0.9411405881823928,
        'Cl': 13.865120686557972,
        'F': 53.69588288801864,
        'H': 0.10792704037665543,
        'N': 6.415529611567397,
        'O': 116.36915086380722,
        'S': 280.7716960911814
    },
    'bond_corr_neighbor': {
        'Br': 0.06608177997464289,
        'C': -0.048677749345191326,
        'Cl': 0.012053484346398953,
        'F': 0.05223732919606618,
        'H': 0.009833437139509741,
        'N': -0.22755394557192776,
        'O': -0.0976858947906553,
        'S': -0.07751080835227722
    },
    'mol_corr': -1.416916143712511
},
# 95% Confidence interval half-widths:
# {
#     'atom_corr': {
#         'Br': 1.8387713358030398,
#         'C': 0.6924596750093204,
#         'Cl': 0.8163744026248948,
#         'F': 0.6942325003578922,
#         'H': 0.46252211009786415,
#         'N': 1.0159737675679787,
#         'O': 1.0650022173657137,
#         'S': 0.7565035051922301
#     },
#     'bond_corr_length': {
#         'Br': 3059.334299611757,
#         'C': 6.242092073505208,
#         'Cl': 115.37303257608322,
#         'F': 76.9573317309656,
#         'H': 0.7450599210452848,
#         'N': 16.118707894644015,
#         'O': 71.9790654412699,
#         'S': 159.4179964314026
#     },
#     'bond_corr_neighbor': {
#         'Br': 0.6989663235292405,
#         'C': 0.022209396236866328,
#         'Cl': 0.28261759032945116,
#         'F': 0.17990689437278468,
#         'H': 0.12318705019325947,
#         'N': 0.05865804713231371,
#         'O': 0.07746518953759682,
#         'S': 0.10222022287246929
#     },
#     'mol_corr': 0.6590311381791739
# }
