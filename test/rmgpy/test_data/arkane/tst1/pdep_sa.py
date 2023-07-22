title = "acetyl + oxygen PDep sensitivity analysis test"

modelChemistry = "CCSD(T)-F12/cc-pVTZ-F12"

species(
    label="acetylperoxy",
    structure=SMILES("CC(=O)O[O]"),
    E0=(-34.6, "kcal/mol"),
    modes=[
        IdealGasTranslation(mass=(75.04, "g/mol")),
        NonlinearRotor(
            inertia=([54.2977, 104.836, 156.05], "amu*angstrom^2"), symmetry=1
        ),
        HarmonicOscillator(
            frequencies=(
                [
                    319.695,
                    500.474,
                    536.674,
                    543.894,
                    727.156,
                    973.365,
                    1037.77,
                    1119.72,
                    1181.55,
                    1391.11,
                    1449.53,
                    1454.72,
                    1870.51,
                    3037.12,
                    3096.93,
                    3136.39,
                ],
                "cm^-1",
            )
        ),
        HinderedRotor(
            inertia=(7.38359, "amu*angstrom^2"),
            symmetry=1,
            fourier=(
                [
                    [-1.95191, -11.8215, 0.740041, -0.049118, -0.464522],
                    [0.000227764, 0.00410782, -0.000805364, -0.000548218, -0.000266277],
                ],
                "kJ/mol",
            ),
        ),
        HinderedRotor(
            inertia=(2.94723, "amu*angstrom^2"),
            symmetry=3,
            fourier=(
                [
                    [0.130647, 0.0401507, -2.54582, -0.0436065, -0.120982],
                    [-0.000701659, -0.000989654, 0.00783349, -0.00140978, -0.00145843],
                ],
                "kJ/mol",
            ),
        ),
    ],
    spinMultiplicity=2,
    opticalIsomers=1,
    molecularWeight=(75.04, "g/mol"),
    collisionModel=TransportData(sigma=(5.09, "angstrom"), epsilon=(473, "K")),
    energyTransferModel=SingleExponentialDown(
        alpha0=(0.5718, "kcal/mol"),
        T0=(300, "K"),
        n=0.85,
    ),
)

species(
    label="hydroperoxylvinoxy",
    structure=SMILES("[CH2]C(=O)OO"),
    E0=(-32.4, "kcal/mol"),
    modes=[
        IdealGasTranslation(mass=(75.04, "g/mol")),
        NonlinearRotor(
            inertia=([44.8034, 110.225, 155.029], "u*angstrom**2"), symmetry=1
        ),
        HarmonicOscillator(
            frequencies=(
                [
                    318.758,
                    420.907,
                    666.223,
                    675.962,
                    752.824,
                    864.66,
                    998.471,
                    1019.54,
                    1236.21,
                    1437.91,
                    1485.74,
                    1687.9,
                    3145.44,
                    3262.88,
                    3434.34,
                ],
                "cm^-1",
            )
        ),
        HinderedRotor(
            inertia=(1.68464, "u*angstrom**2"),
            symmetry=2,
            fourier=(
                [
                    [0.359649, -16.1155, -0.593311, 1.72918, 0.256314],
                    [-7.42981e-06, -0.000238057, 3.29276e-05, -6.62608e-05, 8.8443e-05],
                ],
                "kJ/mol",
            ),
        ),
        HinderedRotor(
            inertia=(8.50433, "u*angstrom**2"),
            symmetry=1,
            fourier=(
                [
                    [-7.53504, -23.4471, -3.3974, -0.940593, -0.313674],
                    [-4.58248, -2.47177, 0.550012, 1.03771, 0.844226],
                ],
                "kJ/mol",
            ),
        ),
        HinderedRotor(
            inertia=(0.803309, "u*angstrom**2"),
            symmetry=1,
            fourier=(
                [
                    [-8.65946, -3.97888, -1.13469, -0.402197, -0.145101],
                    [4.41884e-05, 4.83249e-05, 1.30275e-05, -1.31353e-05, -6.66878e-06],
                ],
                "kJ/mol",
            ),
        ),
    ],
    spinMultiplicity=2,
    opticalIsomers=1,
    molecularWeight=(75.04, "g/mol"),
    collisionModel=TransportData(sigma=(5.09, "angstrom"), epsilon=(473, "K")),
    energyTransferModel=SingleExponentialDown(
        alpha0=(0.5718, "kcal/mol"),
        T0=(300, "K"),
        n=0.85,
    ),
)

species(
    label="nitrogen",
    structure=SMILES("N#N"),
    molecularWeight=(28.04, "g/mol"),
    collisionModel=TransportData(sigma=(3.70, "angstrom"), epsilon=(94.9, "K")),
)

transitionState(
    label="ts1",
    E0=(-5.8, "kcal/mol"),
    modes=[
        IdealGasTranslation(mass=(75.04, "g/mol")),
        NonlinearRotor(
            inertia=([49.3418, 103.697, 149.682], "u*angstrom**2"),
            symmetry=1,
            quantum=False,
        ),
        HarmonicOscillator(
            frequencies=(
                [
                    148.551,
                    306.791,
                    484.573,
                    536.709,
                    599.366,
                    675.538,
                    832.594,
                    918.413,
                    1022.28,
                    1031.45,
                    1101.01,
                    1130.05,
                    1401.51,
                    1701.26,
                    1844.17,
                    3078.6,
                    3163.07,
                ],
                "cm^-1",
            ),
            quantum=True,
        ),
    ],
    spinMultiplicity=2,
    opticalIsomers=1,
    frequency=(-1679.04, "cm^-1"),
)

reaction(
    label="acetylperoxy <=> hydroperoxylvinoxy",
    reactants=["acetylperoxy"],
    products=["hydroperoxylvinoxy"],
    transitionState="ts1",
    tunneling="Eckart",
)

network(
    label="network1",
    isomers=[
        "acetylperoxy",
        "hydroperoxylvinoxy",
    ],
    bathGas={
        "nitrogen": 1.0,
    },
)

pressureDependence(
    "network1",
    Tmin=(300.0, "K"),
    Tmax=(2000.0, "K"),
    Tcount=8,
    Pmin=(0.1, "bar"),
    Pmax=(10, "bar"),
    Pcount=5,
    maximumGrainSize=(5.0, "kcal/mol"),
    minimumGrainCount=100,
    method="modified strong collision",
    interpolationModel=("chebyshev", 6, 4),
    activeJRotor=True,
    sensitivity_conditions=[(1000, "K"), (1, "bar")],
)
