#This syngas example is for testing Model Resurrection handling (recovery from solver errors)
#The parameters chosen for this example were chosen to generate solver errors
#(so this is not a good example to base an input file for a real job on)
database(
	#overrides RMG thermo calculation of RMG with these values.
	#libraries found at http://rmg.mit.edu/database/thermo/libraries/
	#if species exist in multiple libraries, the earlier libraries overwrite the 
	#previous values
    thermoLibraries = ['primaryThermoLibrary'],
	#overrides RMG kinetics estimation if needed in the core of RMG. 
	#list of libraries found at http://rmg.mit.edu/database/kinetics/libraries/
	#libraries can be input as either a string or tuple of form ('library_name',True/False) 
     #where a `True` indicates that all unused reactions will be automatically added 
	#to the chemkin file at the end of the simulation. Placing just string values
     #defaults the tuple to `False`. The string input is sufficient in almost
     #all situations
    reactionLibraries = [],
	#seed mechanisms are reactionLibraries that are forced into the initial mechanism 
	#in addition to species listed in this input file.  
	#This is helpful for reducing run time for species you know will appear in 
	#the mechanism.  
    seedMechanisms = [],
	#this is normally not changed in general RMG runs.  Usually used for testing with 
	#outside kinetics databases
    kineticsDepositories = 'default', 
	#lists specific families used to generate the model. 'default' uses a list of 
	#families from RMG-Database/input/kinetics/families/recommended.py
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
   label='N2',
   reactive=False,
  structure=SMILES("N#N")
)
species(
    label='CO',
    reactive=True,
    structure=SMILES("[C-]#[O+]")
)
    
species(
    label='H2',
    reactive=True,
    structure=SMILES("[H][H]"),
)
species(
        label='O2',
        reactive=True,
        structure=SMILES("[O][O]")
        )

#Reaction systems
#currently RMG models only constant temperature and pressure as homogeneous batch reactors.  
#two options are: simpleReactor for gas phase or liquidReactor for liquid phase
#use can use multiple reactors in an input file for each condition you want to test.  

simpleReactor(
	#specifies reaction temperature with units
    temperature=(1000,'K'),
	#specifies reaction pressure with units
    pressure=(10.0,'bar'),
	#list initial mole fractions of compounds using the label from the 'species' label.  
	#RMG will normalize if sum/=1
    initialMoleFractions={
        "CO": .6,
        "H2": .4,
        "O2": .5,
 	"N2": 0,
    },
	#the following two values specify when to determine the final output model
	#only one must be specified
	#the first condition to be satisfied will terminate the process

    terminationTime=(12,'s'),
)



simpleReactor(
	#specifies reaction temperature with units
    temperature=(1000,'K'),
	#specifies reaction pressure with units
    pressure=(100.0,'bar'),
	#list initial mole fractions of compounds using the label from the 'species' label.  
	#RMG will normalize if sum/=1
    initialMoleFractions={
        "CO": .6,
        "H2": .4,
        "O2": .5,
	"N2": 0,
    },
	#the following two values specify when to determine the final output model
	#only one must be specified
	#the first condition to be satisfied will terminate the process

    terminationTime=(12,'s'),
)





simpleReactor(
              #specifies reaction temperature with units
              temperature=(2000,'K'),
              #specifies reaction pressure with units
              pressure=(100.0,'bar'),
              #list initial mole fractions of compounds using the label from the 'species' label.
              #RMG will normalize if sum/=1
              initialMoleFractions={
              "CO": .6,
              "H2": .4,
              "O2": .5,
	      "N2": 0,
              },
              #the following two values specify when to determine the final output model
              #only one must be specified
              #the first condition to be satisfied will terminate the process
              
              terminationTime=(12,'s'),
              )
simpleReactor(
              #specifies reaction temperature with units
              temperature=(2000,'K'),
              #specifies reaction pressure with units
              pressure=(10.0,'bar'),
              #list initial mole fractions of compounds using the label from the 'species' label.
              #RMG will normalize if sum/=1
              initialMoleFractions={
              "CO": .6,
              "H2": .4,
              "O2": .5,
	      "N2": 0,
              },
              #the following two values specify when to determine the final output model
              #only one must be specified
              #the first condition to be satisfied will terminate the process
              
              terminationTime=(12,'s'),
              )




#1500
simpleReactor(
              #specifies reaction temperature with units
              temperature=(1500,'K'),
              #specifies reaction pressure with units
              pressure=(100.0,'bar'),
              #list initial mole fractions of compounds using the label from the 'species' label.
              #RMG will normalize if sum/=1
              initialMoleFractions={
              "CO": .6,
              "H2": .4,
              "O2": .5,
	      "N2": 0,
              },
              #the following two values specify when to determine the final output model
              #only one must be specified
              #the first condition to be satisfied will terminate the process
              
              terminationTime=(12,'s'),
              )
simpleReactor(
              #specifies reaction temperature with units
              temperature=(1500,'K'),
              #specifies reaction pressure with units
              pressure=(10.0,'bar'),
              #list initial mole fractions of compounds using the label from the 'species' label.
              #RMG will normalize if sum/=1
              initialMoleFractions={
              "CO": .6,
              "H2": .4,
              "O2": .5,
              "N2": 0,
              },
              #the following two values specify when to determine the final output model
              #only one must be specified
              #the first condition to be satisfied will terminate the process
              
              terminationTime=(12,'s'),
              )

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
	#A smaller value will result in a larger, more complex model
	#when running a new model, it is recommended to start with higher values and then decrease to converge on the model
    toleranceMoveToCore=0.01,
    #comment out the next three terms to disable pruning
	   #determines the relative flux needed to not remove species from the model.  
	   #Lower values will keep more species and utilize more memory
    toleranceKeepInEdge=0.0,
	   #determines when to stop a ODE run to add a species.  
	   #Lower values will improve speed. 
	   #if it is too low, may never get to the end simulation to prune species.  
    toleranceInterruptSimulation=0.01,
	   #number of edge species needed to accumulate before pruning occurs
	   #larger values require more memory and will prune less often
    maximumEdgeSpecies=100000,
        #minimum number of core species needed before pruning occurs.
        #this prevents pruning when kinetic model is far away from completeness
    minCoreSizeForPrune=50,
        #make sure that the pruned edge species have existed for a set number of RMG iterations.  
        #the user can specify to increase it from the default value of 2
    minSpeciesExistIterationsForPrune=2,
        #filter the reactions during the enlarge step to omit species from reacting if their
        #concentration are deemed to be too low
    filterReactions=True,
    maxNumSpecies=24,
)

options(
    	#provides a name for the seed mechanism produced at the end of an rmg run default is 'Seed'
    name='Syngas',	
	#if True every iteration it saves the current model as libraries/seeds
	#(and deletes the old one)
	#Unlike HTML this is inexpensive time-wise
	#note a seed mechanism will be generated at the end of a completed run and some incomplete 
	#runs even if this is set as False
    generateSeedEachIteration=True,
	#If True the mechanism will also be saved directly as kinetics and thermo libraries in the database
    saveSeedToDatabase=False,
	#only option is 'si'
    units='si',
	#Draws images of species and reactions and saves the model output to HTML.  
	#May consume extra memory when running large models.
    generateOutputHTML=True,
	#generates plots of the RMG's performance statistics. Not helpful if you just want a model.
    generatePlots=False,
	#saves mole fraction of species in 'solver/' to help you create plots
    saveSimulationProfiles=True,
	#gets RMG to output comments on where kinetics were obtained in the chemkin file.  
	#useful for debugging kinetics but increases memory usage of the chemkin output file
    verboseComments=False,
	#gets RMG to generate edge species chemkin files. Uses lots of memory in output.
	#Helpful for seeing why some reaction are not appearing in core model.  
    saveEdgeSpecies=False,
    #Sets a time limit in the form DD:HH:MM:SS after which the RMG job will stop. Useful for profiling on jobs that
    #do not converge.
    #wallTime = '00:00:00', 
    keepIrreversible=False,
    #Forces RMG to import library reactions as reversible (default). Otherwise, if set to True, RMG will import library
    #reactions while keeping the reversibility as as.
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
 	pressures=(0.01,100.01,'bar',3),
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
    maximumOxygenAtoms=6,
    maximumNitrogenAtoms=0,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
	#max number of non-hydrogen atoms
    #maximumHeavyAtoms=20,
	#maximum radicals on a molecule
    maximumRadicalElectrons=2,
    #If this is false or missing, RMG will throw an error if the more less-stable form of O2 is entered 
    #which doesn't react in the RMG system. normally input O2 as triplet with SMILES [O][O]
    #allowSingletO2=False,
)

