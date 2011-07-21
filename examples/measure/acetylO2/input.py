################################################################################
#
#   MEASURE input file for acetyl + O2 reaction network
#
################################################################################

title = 'acetyl + oxygen'

description = \
"""
The chemically-activated reaction of acetyl with oxygen. This system is of
interest in atmospheric chemistry as a step in the conversion of acetaldehyde
to the secondary pollutant peroxyacetylnitrate (PAN); it is also potentially
important in the ignition chemistry of ethanol.
"""

species(
    label='acetylperoxy',
    SMILES='CC(=O)O[O]',
    E0=(-34.6,'kcal/mol'),
    states=States(
        rotations=RigidRotor(
            linear=False,
            inertia=([54.2978, 104.8364, 156.0495],"amu*angstrom^2"),
            symmetry=1,
        ),
        vibrations=HarmonicOscillator(
            frequencies=([321.607, 503.468, 539.885, 547.148, 731.506, 979.187, 1043.981, 1126.416, 1188.619, 1399.432, 1458.200, 1463.423, 1881.701, 3055.285, 3115.447, 3155.144], 'cm^-1'),
        ),
        torsions=[
            HinderedRotor(inertia=(7.38359,"amu*angstrom^2"), barrier=(6.11665,"kcal/mol"), symmetry=1),
            HinderedRotor(inertia=(2.94725,"amu*angstrom^2"), barrier=(1.22157,"kcal/mol"), symmetry=3),
        ],
        frequencyScaleFactor=0.99,
        spinMultiplicity=2,
    ),
    lennardJones=LennardJones(sigma=(5.09,'angstrom'), epsilon=(473,'K')),
)

species(
    label='hydroperoxylvinoxy',
    SMILES='[CH2]C(=O)OO',
    E0=(-32.4,'kcal/mol'),
    states=States(
        rotations=RigidRotor(
            linear=False,
            inertia=([44.8035, 110.2250, 155.0285],"amu*angstrom^2"),
            symmetry=1,
        ),
        vibrations=HarmonicOscillator(
            frequencies=([320.665, 423.425, 670.208, 680.006, 757.327, 869.833, 1004.445, 1025.639, 1243.606, 1446.509, 1494.622, 1697.996, 3164.246, 3282.391, 3454.878], 'cm^-1'),
        ),
        torsions=[
            HinderedRotor(inertia=(1.68465,"amu*angstrom^2"), barrier=(8.06130,"kcal/mol"), symmetry=2),
            HinderedRotor(inertia=(8.50433,"amu*angstrom^2"), barrier=(14.7920,"kcal/mol"), symmetry=1),
            HinderedRotor(inertia=(0.803313,"amu*angstrom^2"), barrier=(4.78136,"kcal/mol"), symmetry=1),
        ],
        frequencyScaleFactor=0.99,
        spinMultiplicity=2,
    ),
    lennardJones=LennardJones(sigma=(5.09,'angstrom'), epsilon=(473,'K')),
)

species(
    label='acetyl',
    SMILES='C[C]=O',
    E0=(0.0,'kcal/mol'),
    states=States(
        rotations=RigidRotor(
            linear=False,
            inertia=([5.94519, 50.8166, 53.6436],"amu*angstrom^2"),
            symmetry=1,
        ),
        vibrations=HarmonicOscillator(
            frequencies=([467.090, 850.182, 1016.581, 1044.643, 1351.572, 1443.265, 1450.874, 1917.585, 3003.317, 3094.973, 3097.873], 'cm^-1'),
        ),
        torsions=[
            HinderedRotor(inertia=(1.61753,"amu*angstrom^2"), barrier=(0.510278,"kcal/mol"), symmetry=3),
        ],
        frequencyScaleFactor=0.99,
        spinMultiplicity=2,
    ),
)

species(
    label='oxygen',
    SMILES='[O][O]',
    E0=(0.0,'kcal/mol'),
    states=States(
        rotations=RigidRotor(
            linear=True,
            inertia=(11.6056,"amu*angstrom^2"),
            symmetry=2,
        ),
        vibrations=HarmonicOscillator(
            frequencies=([1631.232],"cm^-1"),
        ),
        frequencyScaleFactor=0.99,
        spinMultiplicity=3,
    ),
)

species(
    label='ketene',
    SMILES='C=C=O',
    E0=(-6.6,'kcal/mol'),
)

species(
    label='lactone',
    SMILES='C1OC1(=O)',
    E0=(-30.8,'kcal/mol'),
)

species(
    label='hydroxyl',
    SMILES='[OH]',
    E0=(0.0,'kcal/mol'),
)

species(
    label='hydroperoxyl',
    SMILES='O[O]',
    E0=(0.0,'kcal/mol'),
)

species(
    label='nitrogen',
    SMILES='N#N',
    lennardJones=LennardJones(sigma=(3.70,'angstrom'), epsilon=(94.9,'K')),
    collisionModel = SingleExponentialDown(
        alpha0 = (0.5718,'kcal/mol'),
        T0 = (300,'K'),
        n = 0.85,
    ),
)

################################################################################

isomer('acetylperoxy')

isomer('hydroperoxylvinoxy')

reactants('acetyl', 'oxygen')

################################################################################

reaction(
    reactants=['acetyl', 'oxygen'],
    products=['acetylperoxy'],
    kinetics=Arrhenius(
        A=(2.65e6,'m^3/(mol*s)'), 
        n=0.0, 
        Ea=(0.0,'kcal/mol')
    ),
    transitionState=TransitionState(
        E0=(0.0,'kcal/mol'),
    )
)

reaction(
    reactants=['acetylperoxy'],
    products=['hydroperoxylvinoxy'],
    kinetics=Arrhenius(
        A=(2.31e9,'s^-1'),
        n=0.75,
        Ea=(23.21,'kcal/mol')
    ),
    transitionState=TransitionState(
        E0=(-5.8,'kcal/mol'),
        #states=States(
            #rotations=RigidRotor(
                #linear=False,
                #inertia=([49.3418, 103.6975, 149.6820],"amu*angstrom^2"),
                #symmetry=1,
            #),
            #vibrations=HarmonicOscillator(
                #frequencies=([150.516,  309.455,  487.595,  540.400,  602.963,  679.579,  837.549,  923.888, 1028.647, 1037.741, 1107.523, 1136.811, 1409.904, 1711.388, 1855.194, 3096.997, 3181.924], 'cm^-1'),
            #),
            #frequencyScaleFactor=0.99,
            #spinMultiplicity=2,
        #),
        frequency=(-1672.207,'cm^-1'),
    )
)

reaction(
    reactants=['acetylperoxy'],
    products=['ketene', 'hydroperoxyl'],
    kinetics=Arrhenius(
        A=(2.62e9,'s^-1'),
        n=1.24,
        Ea=(34.06,'kcal/mol')
    ),
    transitionState=TransitionState(
        E0=(0.6,'kcal/mol'),
        #states=States(
            #rotations=RigidRotor(
                #linear=False,
                #inertia=([55.4256, 136.1886, 188.2442],"amu*angstrom^2"),
                #symmetry=1,
            #),
            #vibrations=HarmonicOscillator(
                #frequencies=([59.306,  205.421,  354.483,  468.861,  482.875,  545.574,  657.825,  891.898, 1023.947, 1085.617, 1257.494, 1316.937, 1378.552, 1688.566, 2175.346, 3079.822, 3154.325], 'cm^-1'),
            #),
            #frequencyScaleFactor=0.99,
            #spinMultiplicity=2,
        #),
        frequency=(-1048.9950,'cm^-1'),
    )
)

reaction(
    reactants=['hydroperoxylvinoxy'],
    products=['ketene', 'hydroperoxyl'],
    kinetics=Arrhenius(
        A=(5.33e16,'s^-1'),
        n=-1.02,
        Ea=(29.51,'kcal/mol')
    ),
    transitionState=TransitionState(
        E0=(-4.6,'kcal/mol'),
        #states=States(
            #rotations=RigidRotor(
                #linear=False,
                #inertia=([51.7432, 119.3733, 171.1165],"amu*angstrom^2"),
                #symmetry=1,
            #),
            #vibrations=HarmonicOscillator(
                #frequencies=([   251.808,  386.126,  547.638,  582.452,  598.884,  709.641,  970.482, 1109.848, 1153.775, 1424.448, 1492.202, 1995.657, 3146.702, 3162.637, 3274.761], 'cm^-1'),
            #),
            #frequencyScaleFactor=0.99,
            #spinMultiplicity=2,
        #),
        frequency=(-402.5870,'cm^-1'),
    )
)

reaction(
    reactants=['hydroperoxylvinoxy'],
    products=['lactone', 'hydroxyl'],
    kinetics=Arrhenius(
        A=(1.92e17,'s^-1'),
        n=-1.07,
        Ea=(27.19,'kcal/mol')
    ),
    transitionState=TransitionState(
        E0=(-7.2,'kcal/mol'),
        #states=States(
            #rotations=RigidRotor(
                #linear=False,
                #inertia=([53.2821, 120.4050, 170.1570],"amu*angstrom^2"),
                #symmetry=1,
            #),
            #vibrations=HarmonicOscillator(
                #frequencies=([   153.065,  182.997,  313.611,  350.732,  612.127,  628.112,  810.164,  954.550, 1001.209, 1002.946, 1176.094, 1421.050, 1845.809, 3143.113, 3264.608, 3656.183], 'cm^-1'),
            #),
            #frequencyScaleFactor=0.99,
            #spinMultiplicity=2,
        #),
        frequency=(-615.7468,'cm^-1'),
    )
)

#################################################################################

bathGas = {
    'nitrogen': 1.0,
}

temperatures(Tmin=(300.0,'K'), Tmax=(2000.0,'K'), count=7)
pressures(Pmin=(0.01,'bar'), Pmax=(100.0,'bar'), count=5)
energies(dE=(0.5,'kcal/mol'), count=250)
method('modified strong collision')
#method('reservoir state')
interpolationModel('chebyshev', 6, 4)
