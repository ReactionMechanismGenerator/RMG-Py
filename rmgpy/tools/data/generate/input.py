# Data sources for kinetics
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],  
    seedMechanisms = [],
    kineticsDepositories = 'default', 
    #this section lists possible reaction families to find reactioons with
    kineticsFamilies = ['R_Recombination'],
    kineticsEstimator = 'rate rules',
)

# List all species you want reactions between
species(
    label='Propyl',
    reactive=True,
    structure=SMILES("CC[CH3]"),
)

species(
    label='H',
    reactive=True,
    structure=SMILES("[H]"),
)


# you must list reactor conditions (though this may not effect the output)
simpleReactor(
    temperature=(650,'K'),
    pressure=(10.0,'bar'),
    initialMoleFractions={
        "Propyl": 1,
    },
    terminationConversion={
        'Propyl': .99,
    },
    terminationTime=(40,'s'),
)
