database(
    thermoLibraries = ['DFT_QCI_thermo', 'primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = 'default',
    kineticsFamilies = ['H_Abstraction', 'R_Addition_MultipleBond', 'R_Recombination', 'Disproportionation'],
    kineticsEstimator = 'rate rules',
)

generatedSpeciesConstraints(
    maximumRadicalElectrons = 1,
    maximumCarbonAtoms = 18,
)

species(
    label = "PDD",
    structure = SMILES("CCCCCCCCCCCCc1ccccc1"))

species(
    label = "TOLUENE",
    structure = SMILES("Cc1ccccc1"))

species(
    label = "STYRENE",
    structure = SMILES("C=Cc1ccccc1"))

species(
    label = "ETHBENZ",
    structure = SMILES("CCc1ccccc1"))

species(
    label = "BENZ3",
    structure = SMILES("CCCc1ccccc1"))

species(
    label = "BENZ4",
    structure = SMILES("CCCCc1ccccc1"))

species(
    label = "BENZ5",
    structure = SMILES("CCCCCc1ccccc1"))

species(
    label = "BENZ6",
    structure = SMILES("CCCCCCc1ccccc1"))

species(
    label = "BENZ7",
    structure = SMILES("CCCCCCCc1ccccc1"))

species(
    label = "BENZ8",
    structure = SMILES("CCCCCCCCc1ccccc1"))

species(
    label = "BENZ9",
    structure = SMILES("CCCCCCCCCc1ccccc1"))

species(
    label = "BENZ10",
    structure = SMILES("CCCCCCCCCCc1ccccc1"))

species(
    label = "BENZ11",
    structure = SMILES("CCCCCCCCCCCc1ccccc1"))

species(
    label = "RAD1",
    structure = SMILES("CCCCCCCCCCC[CH]c1ccccc1"))

species(
    label = "RAD2",
    structure = SMILES("CCCCCCCCCC[CH]Cc1ccccc1"))

species(
    label = "RAD3",
    structure = SMILES("CCCCCCCCC[CH]CCc1ccccc1"))

species(
    label = "RAD4",
    structure = SMILES("CCCCCCCC[CH]CCCc1ccccc1"))

species(
    label = "RAD5",
    structure = SMILES("CCCCCCC[CH]CCCCc1ccccc1"))

species(
    label = "RAD6",
    structure = SMILES("CCCCCC[CH]CCCCCc1ccccc1"))

species(
    label = "RAD7",
    structure = SMILES("CCCCC[CH]CCCCCCc1ccccc1"))

species(
    label = "RAD8",
    structure = SMILES("CCCC[CH]CCCCCCCc1ccccc1"))

species(
    label = "RAD9",
    structure = SMILES("CCC[CH]CCCCCCCCc1ccccc1"))

species(
    label = "RAD10",
    structure = SMILES("CC[CH]CCCCCCCCCc1ccccc1"))

species(
    label = "RAD11",
    structure = SMILES("C[CH]CCCCCCCCCCc1ccccc1"))

species(
    label = "RAD12",
    structure = SMILES("[CH2]CCCCCCCCCCCc1ccccc1"))

species(
    label = "C1",
    structure = SMILES("C"))

species(
    label = "C2",
    structure = SMILES("CC"))

species(
    label = "C3",
    structure = SMILES("CCC"))

species(
    label = "C4",
    structure = SMILES("CCCC"))

species(
    label = "C5",
    structure = SMILES("CCCCC"))

species(
    label = "C6",
    structure = SMILES("CCCCCC"))

species(
    label = "C7",
    structure = SMILES("CCCCCCC"))

species(
    label = "C8",
    structure = SMILES("CCCCCCCC"))

species(
    label = "C9",
    structure = SMILES("CCCCCCCCC"))

species(
    label = "C10",
    structure = SMILES("CCCCCCCCCC"))

species(
    label = "C11",
    structure = SMILES("CCCCCCCCCCC"))

species(
    label = "C2ene",
    structure = SMILES("C=C"))

species(
    label = "C3ene",
    structure = SMILES("CC=C"))

species(
    label = "C4ene",
    structure = SMILES("CCC=C"))

species(
    label = "C5ene",
    structure = SMILES("CCCC=C"))

species(
    label = "C6ene",
    structure = SMILES("CCCCC=C"))

species(
    label = "C7ene",
    structure = SMILES("CCCCCC=C"))

species(
    label = "C8ene",
    structure = SMILES("CCCCCCC=C"))

species(
    label = "C9ene",
    structure = SMILES("CCCCCCCC=C"))

species(
    label = "C10ene",
    structure = SMILES("CCCCCCCCC=C"))

species(
    label = "C11ene",
    structure = SMILES("CCCCCCCCCC=C"))

species(
    label = "METHYL",
    structure = SMILES("[CH3]"))

species(
    label = "ETHYL",
    structure = SMILES("C[CH2]"))

species(
    label = "PROPYL",
    structure = SMILES("CC[CH2]"))

species(
    label = "BUTYL",
    structure = SMILES("CCC[CH2]"))

species(
    label = "PENTYL",
    structure = SMILES("CCCC[CH2]"))

species(
    label = "HEXYL",
    structure = SMILES("CCCCC[CH2]"))

species(
    label = "HEPTYL",
    structure = SMILES("CCCCCC[CH2]"))

species(
    label = "OCTYL",
    structure = SMILES("CCCCCCC[CH2]"))

species(
    label = "NONYL",
    structure = SMILES("CCCCCCCC[CH2]"))

species(
    label = "DECYL",
    structure = SMILES("CCCCCCCCC[CH2]"))

species(
    label = "UDECYL",
    structure = SMILES("CCCCCCCCCC[CH2]"))

species(
    label = "BENZYL",
    structure = SMILES("[CH2]c1ccccc1"))

species(
    label = "EBZYL",
    structure = SMILES("[CH2]Cc1ccccc1"))

species(
    label = "EBZYL2",
    structure = SMILES("C[CH]c1ccccc1"))

species(
    label = "A3yl",
    structure = SMILES("[CH2]CCc1ccccc1"))

species(
    label = "A4yl",
    structure = SMILES("[CH2]CCCc1ccccc1"))

species(
    label = "A5yl",
    structure = SMILES("[CH2]CCCCc1ccccc1"))

species(
    label = "A6yl",
    structure = SMILES("[CH2]CCCCCc1ccccc1"))

species(
    label = "A7yl",
    structure = SMILES("[CH2]CCCCCCc1ccccc1"))

species(
    label = "A8yl",
    structure = SMILES("[CH2]CCCCCCCc1ccccc1"))

species(
    label = "A9yl",
    structure = SMILES("[CH2]CCCCCCCCc1ccccc1"))

species(
    label = "A10yl",
    structure = SMILES("[CH2]CCCCCCCCCc1ccccc1"))

species(
    label = "A3ene",
    structure = SMILES("C=CCc1ccccc1"))

species(
    label = "A4ene",
    structure = SMILES("C=CCCc1ccccc1"))

species(
    label = "A5ene",
    structure = SMILES("C=CCCCc1ccccc1"))

species(
    label = "A6ene",
    structure = SMILES("C=CCCCCc1ccccc1"))

species(
    label = "A7ene",
    structure = SMILES("C=CCCCCCc1ccccc1"))

species(
    label = "A8ene",
    structure = SMILES("C=CCCCCCCc1ccccc1"))

species(
    label = "A9ene",
    structure = SMILES("C=CCCCCCCCc1ccccc1"))

species(
    label = "A10ene",
    structure = SMILES("C=CCCCCCCCCc1ccccc1"))


simpleReactor(
    temperature=(623, "K"),  # 350 C
    pressure=(350, "bar"),
    initialMoleFractions={
        "PDD": 1
    },
    terminationTime=(72, "h"),
    terminationConversion={
        "PDD": 0.2,
    },
    sensitivity=['PDD', 'C11ene'],
)

uncertainty(
    globalAnalysis=True,
    pceRunTime=3600,
)

simulator(
    atol=1e-16,
    rtol=1e-08,
    sens_atol=1e-06,
    sens_rtol=0.0001,
)

model(
    toleranceMoveToCore=0.5,
    toleranceInterruptSimulation=0.5,
    filterReactions=True,
)

options(
    units="si",
    saveRestartPeriod=None,
    generateOutputHTML=False,
    generatePlots=False,
    saveSimulationProfiles=True,
    saveEdgeSpecies=False,
    verboseComments=False,
)

