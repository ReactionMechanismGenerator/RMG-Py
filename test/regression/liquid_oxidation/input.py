# Test liquid-phase simulation with reaction filtering and flux pruning
# Also test adding multiple species in each iteration and imposing species constraints.

database(
    thermoLibraries=['primaryThermoLibrary'],
    reactionLibraries=[],
    seedMechanisms=[],
    kineticsDepositories='default',
    kineticsFamilies=['default', 'liquid_peroxide'],
    kineticsEstimator='rate rules',
)

generatedSpeciesConstraints(
    allowed=['input species'],
    maximumRadicalElectrons=1,
    maximumCarbonAtoms=5,
    maximumOxygenAtoms=4,
)

species(
    label='oxygen',
    reactive=True,
    structure=SMILES('[O][O]')
)

species(
    label='pentane',
    reactive=True,
    structure=SMILES('CCCCC'),
)

liquidReactor(
    temperature=(600, 'K'),
    initialConcentrations={
        'pentane': (3.0e-3, 'mol/cm^3'),
        'oxygen': (1.0e-4, 'mol/cm^3'),
    },
    constantSpecies=['oxygen'],
    terminationConversion={'pentane': 0.3},
    terminationTime=(1000, 's'),
)

solvation(
    solvent='pentane'
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceMoveToCore=0.01,
    toleranceKeepInEdge=0.001,
    toleranceInterruptSimulation=1e8,
    maximumEdgeSpecies=10000,
    minCoreSizeForPrune=10,
    minSpeciesExistIterationsForPrune=2,
    maxNumObjsPerIter=3,
    filterReactions=True,
    maxNumSpecies=35,
)

options(
    saveEdgeSpecies=True,
)
