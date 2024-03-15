# Data sources
database(
    thermoLibraries=['LithiumSurface','electrocatLiThermo','primaryThermoLibrary', 'LithiumPrimaryThermo', 'LithiumAdditionalThermo', 'thermo_DFT_CCSDTF12_BAC','DFT_QCI_thermo'], # 'surfaceThermoPt' is the default. Thermo data is derived using bindingEnergies for other metals
    reactionLibraries = ['LithiumPrimaryKinetics',"LithiumSurface"], # when Ni is used change the library to Surface/Deutschmann_Ni
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['surface',
			'1+2_Cycloaddition',
    'Surface_Carbonate_Deposition',
    'Surface_Carbonate_F_CO_Decomposition',
    'Surface_Carbonate_2F_Decomposition',
    'Surface_Carbonate_CO_Decomposition',
    '1,2_Elimination_LiR',
    '1,2_Intra_Elimination_LiR',
    'Li_Addition_MultipleBond',
    'Li_NO_Substitution',
    'Li_NO_Ring_Opening',
    'Li_Abstraction',
    'R_Addition_MultipleBond_Disprop',
    'Cation_R_Recombination',
    'Cation_Addition_MultipleBond',
    '1,2-Birad_to_alkene',
    '1,2_Insertion_CO',
    '1,2_Insertion_carbene',
    '1,2_shiftS',
    '1,3_Insertion_CO2',
    '1,3_Insertion_ROR',
    '1,3_Insertion_RSR',
    '1,4_Cyclic_birad_scission',
    '1,4_Linear_birad_scission',
    '2+2_cycloaddition',
    'Birad_recombination',
    'CO_Disproportionation',
    'Birad_R_Recombination',
    'Cyclic_Ether_Formation',
    'Cyclic_Thioether_Formation',
    'Diels_alder_addition',
    'Diels_alder_addition_Aromatic',
    #'Disproportionation',
    'HO2_Elimination_from_PeroxyRadical',
    'H_Abstraction',
    'Intra_Retro_Diels_alder_bicyclic',
    'Intra_Disproportionation',
    'Intra_R_Add_Endocyclic',
    'Intra_R_Add_Exocyclic',
    'R_Addition_COm',
    'R_Addition_MultipleBond',
    'R_Recombination',
    'intra_H_migration',
    'intra_NO2_ONO_conversion',
    'intra_OH_migration',
    'intra_substitutionCS_cyclization',
    'intra_substitutionCS_isomerization',
    'intra_substitutionS_cyclization',
    'intra_substitutionS_isomerization',
    #'ketoenol',
    'Singlet_Carbene_Intra_Disproportionation',
    'Singlet_Val6_to_triplet',
    'Intra_5_membered_conjugated_C=C_C=C_addition',
    'Intra_Diels_alder_monocyclic',
    'Concerted_Intra_Diels_alder_monocyclic_1,2_shiftH',
    'Intra_2+2_cycloaddition_Cd',
    'Intra_ene_reaction',
    'Cyclopentadiene_scission',
    '6_membered_central_C-C_shift',
    'Intra_R_Add_Exo_scission',
    '1,2_shiftC',
    '1,2_NH3_elimination',
    '1,3_NH3_elimination',
    'Retroene',],
    kineticsEstimator = 'rate rules',
    adsorptionGroups='adsorptionLi',
)

catalystProperties(
    metal = 'Li110',
)

# List of species
species(
    label="Lip",
    reactive=True,
    structure=SMILES("[Li+]"),
)

species(
    label='ACN',
    reactive=True,
    structure=SMILES("CC#N"),
)

species(
    label='vacantX',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)

liquidSurfaceReactor(
    temperature=(298.15,'K'),
    distance=(10.0e-10,"m"),
    viscosity=(5e7,"Pa*s"),
    liqPotential=(0.3,'V'),
    surfPotential=(0.0,'V'),
    initialConcentrations={
        "ACN": (0.019146,'mol/cm^3'),
        "Lip": (15.0,'mol/m^3'),
    },
	initialSurfaceCoverages={
       "vacantX": 1.0, 
    },
    surfaceVolumeRatio=(1.0e-5, 'm^-1'),
    terminationTime=(1e3,'sec'),
	constantSpecies=["ACN","Lip"],
)

liquidSurfaceReactor(
    temperature=(298.15,'K'),
    distance=(0.0,"m"),
    liqPotential=(0.0,'V'),
    surfPotential=(0.0,'V'),
    initialConcentrations={
        "ACN": (0.019146,'mol/cm^3'),
        "Lip": (15.0,'mol/m^3'),
    },
	initialSurfaceCoverages={
        "vacantX": 1.0, 
    },
    surfaceVolumeRatio=(1.0e5, 'm^-1'),
    terminationTime=(1e3,'sec'),
	constantSpecies=["ACN","Lip"],
)

solvation(
	solvent='acetonitrile'
)

simulator(
    atol=1e-16,
    rtol=1e-6,
)

model(
    toleranceKeepInEdge=1E-20,
    toleranceMoveToCore=0.1,
    toleranceRadMoveToCore=0.1,
    toleranceInterruptSimulation=1e10,
    maximumEdgeSpecies=100000,
    filterReactions=False,
    maxNumObjsPerIter=1,
    terminateAtMaxObjects=True,
    toleranceBranchReactionToCore=0.001,
    branchingIndex=0.5,
    branchingRatioMax=1.0,
)

options(
    units='si',
    saveEdgeSpecies=False,
)

forbidden(
    label='vacancies',
    structure=adjacencyListGroup("""
1 Xv u0 p0 c0
"""),
)

forbidden(
    label='Li2',
    structure=adjacencyList("""
1 Li u0 p0 c0 {2,S}
2 Li u0 p0 c0 {1,S}"""),
)

generatedSpeciesConstraints(
    allowed=['input species','reaction libraries'],
    maximumSurfaceSites=1,
    maximumCarbonAtoms=7,
    maximumOxygenAtoms=4,
    maximumRadicalElectrons=1,
)
