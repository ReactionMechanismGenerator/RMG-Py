# Data sources for kinetics
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],  
    seedMechanisms = [],
    kineticsDepositories = 'default', 
    #this section lists possible reaction families to find reactioons with
    kineticsFamilies = ['!Intra_Disproportionation','!Substitution_O'],
    kineticsEstimator = 'rate rules',
)

# List all species you want reactions between
species(
    label='ethane',
    reactive=True,
    structure=SMILES("CC"),
)

species(
    label='H',
    reactive=True,
    structure=SMILES("[H]"),
)

species(
    label='butane',
    reactive=True,
    structure=SMILES("CCCC"),
)


# you must list reactor conditions (though this may not effect the output)
simpleReactor(
    temperature=(650,'K'),
    pressure=(10.0,'bar'),
    initialMoleFractions={
        "ethane": 1,
    },
    terminationConversion={
        'butane': .99,
    },
    terminationTime=(40,'s'),
)