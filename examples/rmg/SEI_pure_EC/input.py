# Data sources
database(
    thermoLibraries=['electrocatLiThermo','primaryThermoLibrary', 'LithiumPrimaryThermo', 'LithiumAdditionalThermo', 'thermo_DFT_CCSDTF12_BAC','DFT_QCI_thermo'], # 'surfaceThermoPt' is the default. Thermo data is derived using bindingEnergies for other metals
    reactionLibraries = ['LithiumPrimaryKinetics','LithiumAnalogyKinetics'], # when Ni is used change the library to Surface/Deutschmann_Ni
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['surface','electrochem',
			'1+2_Cycloaddition',
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
    metal = 'Li110'
)

# List of species
species(
    label="Lip",
    reactive=True,
    structure=SMILES("[Li+]"),
)

species(
    label='ethylene-carbonate',
    reactive=True,
    structure=SMILES("C1COC(=O)O1"),
)

species(
    label='vacantX',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)

species(
     label="Li",
     reactive=True,
     structure=SMILES("[Li]"),
)

species(
      label='[Li]O[C]1OCCO1',
      reactive=True,
      structure=SMILES("[Li]O[C]1OCCO1"),
)

species(
      label='[Li]OC(=O)OC[CH2]',
      reactive=True,
      structure=SMILES("[Li]OC(=O)OC[CH2]"),
)

#species(
#      label='[Li]OC(=O)O[Li]',
#      reactive=True,
#      structure=SMILES("[Li]OC(=O)O[Li]"),
# )

#species(
#     label='[Li]OC(=O)OCCOC(=O)OC[CH2]',
#     reactive=True,
#     structure=SMILES("[Li]OC(=O)OCCOC(=O)OC[CH2]"),
# )

#species(
#     label='[Li]OC(=O)OCCOC(=O)O[Li]',
#     reactive=True,
#     structure=SMILES("[Li]OC(=O)OCCOC(=O)O[Li]"),
# )

#species(
#     label='[Li]OC(=O)OCCCCOC(=O)O[Li]',
#     reactive=True,
#     structure=SMILES("[Li]OC(=O)OCCCCOC(=O)O[Li]"),)

#species(
#     label='[Li]OCCOC(=O)CCOC(=O)O[Li]',
#     reactive=True,
#     structure=SMILES("[Li]OCCOC(=O)CCOC(=O)O[Li]"),
# )

#species(
#     label='[Li]OCCOC(=O)OC(=O)O[Li]',
#     reactive=True,
#     structure=SMILES("[Li]OCCOC(=O)OC(=O)O[Li]"),
#)

#species(
#      label='C2H4',
#      reactive=True,
#      structure=SMILES("C=C"),
#)

#species(
#       label='O=[C]OCCO[Li]',
#       reactive=True,
#       structure=SMILES("O=[C]OCCO[Li]"),
#)

#species(
#     label='CO2',
#     reactive=True,
#     structure=SMILES("O=C=O"),
#)

#species(
#      label='[Li]OC[CH2]',
#      reactive=True,
#      structure=SMILES("[Li]OC[CH2]"),
#)

#species(
#     label='O1CCO[C]1OC2(O[Li])OCCO2',
#     reactive=True,
#     structure=SMILES("O1CCO[C]1OC2(O[Li])OCCO2"),
#)

#species(
#     label='O1CCOC1(O[Li])OC(=O)OC[CH2]',
#     reactive=True,
#     structure=SMILES("O1CCOC1(O[Li])OC(=O)OC[CH2]"),
#)

#species(
#     label='O1CCOC1(O[Li])OC(=O)OCCOC(=O)O[Li]',
#     reactive=True,
#     structure=SMILES("O1CCOC1(O[Li])OC(=O)OCCOC(=O)O[Li]"),
#)

species(
    label='CO3X2',
    reactive=True,
    structure=adjacencyList("""1     O u0 p2 {2,D}
2     C u0 p0 {1,D} {3,S} {4,S}
3     O u0 p2 {2,S} {5,S}
4     O u0 p2 {2,S} {6,S}
5     X  u0 p0 c0 {3,S}
6     X  u0 p0 c0 {4,S}
"""),
)


#species(
# label="CO",
# reactive=True,
# structure=SMILES("[C-]#[O+]"),
#)

#species(
#     label='[Li]OC(=O)OCCX',
#     reactive=True,
#     structure=adjacencyList("""1  O u0 p2 c0 {2,S} {7,S}
# 2  C u0 p0 c0 {1,S} {3,D} {4,S}
# 3  O u0 p2 c0 {2,D}
# 4  O u0 p2 c0 {2,S} {5,S}
# 5  C u0 p0 c0 {4,S} {6,S} {8,S} {9,S}
# 6  C u0 p0 c0 {5,S} {10,S} {11,S} {12,S}
# 7  Li u0 p0 c0 {1,S}
# 8  H u0 p0 c0 {5,S}
# 9  H u0 p0 c0 {5,S}
# 10 H u0 p0 c0 {6,S}
# 11 H u0 p0 c0 {6,S}
# 12 X u0 p0 c0 {6,S}
# """),
#)
#species(
#     label='O=C(X)OCCO[Li]',
#     reactive=True,
#     structure=adjacencyList("""1  O u0 p2 c0 {2,D}
# 2  C u0 p0 c0 {1,D} {3,S} {12,S}
# 3  O u0 p2 c0 {2,S} {4,S}
# 4  C u0 p0 c0 {3,S} {5,S} {7,S} {8,S}
# 5  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
# 6  O u0 p2 c0 {5,S} {11,S}
# 7  H u0 p0 c0 {4,S}
# 8  H u0 p0 c0 {4,S}
# 9  H u0 p0 c0 {5,S}
# 10 H u0 p0 c0 {5,S}
# 11 Li u0 p0 c0 {6,S}
# 12 X u0 p0 c0 {2,S}
# """),
#)

liquidSurfaceReactor(
    temperature=(298.15,'K'),
    liqPotential=(-1.0,'V'),
    surfPotential=(0.0,'V'),
    initialConcentrations={
        "ethylene-carbonate": (7.585e-3*2.0,'mol/cm^3'),
        "Lip": (15.0,'mol/m^3'),
    },
	initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.0e5, 'm^-1'),
    terminationTime=(1.0e3,'sec'),
	constantSpecies=["ethylene-carbonate","Lip"],
)

solvation(
	solvent='ethylene carbonate'
)

simulator(
    atol=1e-16,
    rtol=1e-6,
)

model(
    toleranceKeepInEdge=1E-20,
    toleranceMoveToCore=0.000001,
    toleranceRadMoveToCore=0.00000000001,
    toleranceInterruptSimulation=1e10,
    maximumEdgeSpecies=100000,
    filterReactions=False,
    maxNumObjsPerIter=1,
    terminateAtMaxObjects=True,
    toleranceBranchReactionToCore=0.000001,
    branchingIndex=0.3,
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

generatedSpeciesConstraints(
    allowed=['input species','reaction libraries'],
    maximumSurfaceSites=1,
    maximumCarbonAtoms=8,
    maximumOxygenAtoms=8,
    maximumRadicalElectrons=1,
)
