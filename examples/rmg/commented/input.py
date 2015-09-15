# Data sources
database(
	#overrides RMG thermo calculation of RMG with these values.
	#libraries found at http://rmg.mit.edu/database/thermo/libraries/
	#if species exist in multiple libraries, the earlier libraries overwrite the 
	#previous values
    thermoLibraries = ['KlippensteinH2O2','primaryThermoLibrary','DFT_QCI_thermo','CBS_QB3_1dHR'],
	#overrides RMG kinetics estimation if needed in the core of RMG. 
	#list of libraries found at http://rmg.mit.edu/database/kinetics/libraries/
	#input each library as a ('library_name',True/False) where a True means that all 
	#unused reactions will be automatically added to the chemkin file
    reactionLibraries = [],
	#seed mechanisms are reactionLibraries that are forced into the initial mechanism 
	#in addition to species listed in this input file.  
	#This is helpful for reducing run time for species you know will appear in 
	#the mechanism.  
    seedMechanisms = ['KlippensteinH2O2','ERC-FoundationFuelv0.9'],
	#this is normally not changed in general RMG runs.  Usually used for testing with 
	#outside kinetics databases
    kineticsDepositories = 'default', 
	#lists specific families used to generate the model. 'default' uses a list of 
	#families from RMG-Database/input/families/recommended.py
	#a visual list of families is available in PDF form at RMG-database/families
    kineticsFamilies = 'default',
	#specifies how RMG calculates rates.  currently, the only option is 'rate rules'
    kineticsEstimator = 'rate rules',
)

# List of species
#list initial and expected species below to automatically put them into the core mechanism.  
#'structure' can utilize method of SMILES("put_SMILES_here"), 
#adjacencyList("""put_adj_list_here"""), or InChI("put_InChI_here")
#for molecular oxygen, use the smiles string [O][O] so the triplet form is used
species(
    label='butane',
    reactive=True,		#this parameter is optional if true
    structure=SMILES("CCCC"),
)
species(
    label='O2',
    structure=SMILES("[O][O]"),
)
species(
    label='N2',
    reactive=False,
    structure=adjacencyList("""
    1 N u0 p1 c0 {2,T}
	2 N u0 p1 c0 {1,T}
	"""),
)
# You can list species not initially in reactor to make sure RMG includes them in the mechanism
species(
    label='QOOH',
    reactive=True,
    structure=SMILES("OOCC[CH]C")
)
species(
    label='CO2',
    reactive=True,
    structure=SMILES("O=C=O")
)

#Reaction systems
#currently RMG models only constant temperature and pressure as homogeneous batch reactors.  
#two options are: simpleReactor for gas phase or liquidReactor for liquid phase
#use can use multiple reactors in an input file for each condition you want to test.  
simpleReactor(
	#specifies reaction temperature with units
    temperature=(700,'K'),
	#specifies reaction pressure with units
    pressure=(10.0,'bar'),
	#list initial mole fractions of compounds using the label from the 'species' label.  
	#RMG will normalize if sum/=1
    initialMoleFractions={
        "N2": 4,
        "O2": 1,
        "butane": 1./6.5,
    },
	#the following two values specify when to determine the final output model
	#only one must be specified
	#the first condition to be satisfied will terminate the process
    terminationConversion={
        'butane': .99,
    },
    terminationTime=(40,'s'),
	#the next two optional values specify how RMG computes sensitivities of 
	#rate coefficients with respect to species concentrations.  
	#sensitivity contains a list of species' labels to conduct sensitivity analysis on.
	#sensitvityThreshold is the required sensitiviy to be recorded in the csv output file
#    sensitivity=['CH4'],
#    sensitivityThreshold=0.0001,
)

# liquidReactor(
# 	temperature=(500,'K'),
#     initialConcentrations={
# 		"N2": 4,
# 		"O2": 1,
# 		"CO": 1,
#     },
#     terminationConversion=None,
#     terminationTime=(3600,'s'),
#     sensitivity=None,
#     sensitivityThreshold=1e-3
# )
	#liquid reactors also have solvents, you can specify one solvent
	#list of solvents available at : http://rmg.mit.edu/database/solvation/libraries/solvent/
# solvation('water')

#determines absolute and relative tolerances for ODE solver and sensitivities.
#normally this doesn't cause many issues and is modified after other issues are
#ruled out
simulator(
    atol=1e-16,
    rtol=1e-8,
#    sens_atol=1e-6,
#    sens_rtol=1e-4,
)

#used to add species to the model and to reduce memory usage by removing unimportant additional species.
#all relative values are normalized by a characteristic flux at that time point
model(
	#determines the relative flux to put a species into the core.  
	#A higher value will result in a larger, more complex model
	#when running a new model, it is recommended to start with higher values and then decrease to converge on the model
    toleranceMoveToCore=0.1,
    #comment out the next three terms to disable pruning
	   #determines the relative flux needed to not remove species from the model.  
	   #Lower values will keep more species and utilize more memory
    toleranceKeepInEdge=0.01,
	   #determines when to stop a ODE run to add a species.  
	   #Lower values will improve speed. 
	   #if it is too low, may never get to the end simulation to prune species.  
    toleranceInterruptSimulation=1,
	   #number of edge species needed to accumulate before pruning occurs
	   #larger values require more memory and will prune less often
    maximumEdgeSpecies=100000
)

options(
	#only option is 'si'
    units='si',
	#how often you want to save restart files.  
	#takes significant amount of time. comment out if you don't want to save 
    saveRestartPeriod=None,
	#Draws images of species and reactions and saves the model output to HTML.  
	#May consume extra memory when running large models.
    drawMolecules=True,
	#generates plots of the RMG's performance statistics. Not helpful if you just want a model.
    generatePlots=False,
	#saves mole fraction of species in 'solver/' to help you create plots
    saveSimulationProfiles=False,
	#gets RMG to output comments on where kinetics were obtained in the chemkin file.  
	#useful for debugging kinetics but increases memory usage of the chemkin output file
    verboseComments=False,
	#gets RMG to generate edge species chemkin files. Uses lots of memory in output.
	#Helpful for seeing why some reaction are not appearing in core model.  
    saveEdgeSpecies=False,
)

# optional module allows for correction to unimolecular reaction rates at low pressures and/or temperatures.
pressureDependence(
 	#two methods available: 'modified strong collision' is faster and less accurate than 'reservoir state'
 	method='modified strong collision',
 	#these two categories determine how fine energy is descretized.
 	#more grains increases accuracy but takes longer
 	maximumGrainSize=(0.5,'kcal/mol'),
 	minimumNumberOfGrains=250,
 	#the conditions for the rate to be output over
 	#parameter order is: low_value, high_value, units, internal points
 	temperatures=(300,2200,'K',2),
 	pressures=(0.01,100,'bar',3),
 	#The two options for interpolation are 'PDepArrhenius' (no extra arguments) and 
 	#'Chebyshev' which is followed by the number of basis sets in 
 	#Temperature and Pressure. These values must be less than the number of 
 	#internal points specified above
 	interpolation=('Chebyshev', 6, 4),
 	#turns off pressure dependence for molecules with number of atoms greater than the number specified below
 	#this is due to faster internal rate of energy transfer for larger molecules
 	maximumAtoms=15,
 )

#optional block adds constraints on what RMG can output.  
#This is helpful for improving the efficiency of RMG, but wrong inputs can lead to many errors.
generatedSpeciesConstraints(
	#allows exceptions to the following restrictions
    allowed=['input species','seed mechanisms','reaction libraries'],
	#maximum number of each atom in a molecule
    maximumCarbonAtoms=4,
    maximumHydrogenAtoms=10,
    maximumOxygenAtoms=7,
    maximumNitrogenAtoms=0,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
	#max number of non-hydrogen atoms
    #maximumHeavyAtoms=20,
	#maximum radicals on a molecule
    maximumRadicalElectrons=1,
    #If this is false or missing, RMG will throw an error if the more less-stable form of O2 is entered 
    #which doesn't react in the RMG system. normally input O2 as triplet with SMILES [O][O]
    #allowSingletO2=False,
)

#optional block allows thermo to be estimated through quantum calculations
# quantumMechanics(
# 	#the software package for calculations...can use 'mopac' or 'gaussian' if installed
# 	software='mopac',
# 	#methods available for calculations. 'pm2' 'pm3' or 'pm7' (last for mopac only)        
# 	method='pm3',
# 	#where to store calculations
# 	fileStore='QMfiles',
# 	#where to store temporary run files
# 	scratchDirectory = None,
# 	#onlyCyclics allows linear molecules to be calculated using bensen group addivity....need to verify
# 	onlyCyclics = True,
# 	#how many radicals should be utilized in the calculation.  
# 	#If the amount of radicals is more than this, RMG will use hydrogen bond incrementation method
# 	maxRadicalNumber = 0,
# )
