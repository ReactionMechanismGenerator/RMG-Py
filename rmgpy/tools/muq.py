
    
from time import time 
import copy   
import numpy
from rmgpy.tools.canteraModel import Cantera
try:
    #import muq.libmuqModelling
    from libmuqModelling import ModPiece
    #import muq.libmuqUtilities as libmuqUtilities
    #import muq.libmuqApproximation as libmuqApproximation
    from libmuqUtilities import LegendrePolynomials1DRecursive, GaussPattersonQuadrature1D, VariableCollection
    from libmuqApproximation import SmolyakPCEFactory
except:
    ModPiece = object  
    print 'Could not import MUQ.  Please check that it is installed correctly before using the global uncertainty modules.'

# You must install the MUQ library before using this.  Add the folder containing
# libmuqUtilities.so, libmuqApproximation.so, etc to your $PYTHONPATH
# For linux users, you can install via 'conda install -c rmg muq' to your environment
# and add the ~/anaconda/envs/your_env/lib folder to your $PYTHONPATH

class ReactorModPiece(ModPiece):
    def __init__(self, cantera, outputSpeciesList, kReactions, kUncertainty, gSpecies, gUncertainty):
        """
        ======================= ====================================================
        Attribute               Description
        ======================= ====================================================
        `cantera`               A Cantera() object containing CanteraConditions and initialized species and reactions
        `outputSpeciesList`     A list of Species() objects corresponding to the desired observables for uncertainty analysis
        `kReactions`            A list of Reaction() objects corresponding to the uncertain input rate coefficients
        `kUncertainty`          A list of uncertainties dlnk corresponding to the reactions in kReactions
        `gSpecies`              A list of Species() objects corresponding to the uncertain input free energies of individual species
        `gUncertainty`          A list of uncertainties dG corresponding to the species in gSpecies in units of kcal/mol
        ============================================================================
        """
        self.cantera = cantera
        self.outputSpeciesList = outputSpeciesList
        self.outputSpeciesIndices = [cantera.speciesList.index(outputSpecies) for outputSpecies in outputSpeciesList]
        
        
        self.kReactions = kReactions
        kUncertaintyFactors = [val*numpy.sqrt(3)/numpy.log(10) for val in kUncertainty]
        self.kUncertaintyFactors = {}
        for i, rxnIndex in enumerate(kReactions):
            self.kUncertaintyFactors[rxnIndex] = kUncertaintyFactors[i]
            
        
        self.gSpecies = gSpecies
        gUncertaintyFactors = [val*numpy.sqrt(3) for val in gUncertainty]
        self.gUncertaintyFactors = {}
        for i, spcIndex in enumerate(gSpecies):
            self.gUncertaintyFactors[spcIndex] = gUncertaintyFactors[i]
        
        # The size of the uncertain inputs: [reaction rates k, species free energy G]         
        self.inputSize = [len(kReactions) + len(gSpecies)]
        
        # The size of the vector corresponding to the outputs to be analyzed for uncertainty analysis
        # is equal to the number of cantera conditions involved multiplied by the number of desired observables
        outputSize = len(cantera.conditions)*len(outputSpeciesList)
        
        self.numOutputSpecies = len(outputSpeciesList)
        self.numConditions = len(cantera.conditions)
        
        # Initialize the ModPiece with some input and output size specifications
        ModPiece.__init__(self,
                                         self.inputSize, 
                                         outputSize, 
                                         False, # No GradientImpl
                                         False, # No JacobianImpl
                                         False, # No HessianImpl
                                         False, # Not random)
                                         False, # Not random)
                                        )
        
    def EvaluateImpl(self, ins):
        """
        Evaluate the desired output mole fractions based on a set of inputs ins = [[k_rv], [G_rv]] which contains the 
        random normal variables attributed to the uncertain kinetics and free energy parameters, respectively.
        
        The output returned contains [Condition1_outputMoleFraction1, Condition1_outputMoleFraction2, Condition2_output.... ConditionN_output...]
        """
        assert len(ins[0]) == self.inputSize[0], "Number of inputs matches number of uncertain parameters"
        
        k_rv = ins[0][0:len(self.kReactions)]
        G_rv = ins[0][len(self.kReactions):]
        
        ## Check that the number of inputs is correct
        #assert len(k_rv) == len(self.kReactions), "Number of inputs matches number of kReactions"
        #assert len(G_rv) == len(self.gSpecies), "Number of inputs matches number of gSpecies"
        
        # Make deepcopies of the thermo and kinetics so as to not modify the originals in the speciesList and reactionList
        originalThermo = [copy.deepcopy(self.cantera.speciesList[index].thermo) for index in self.gSpecies]
        originalKinetics = [copy.deepcopy(self.cantera.reactionList[index].kinetics) for index in self.kReactions]
    
#         print ''
#         print 'Kinetics before'
#         ctReactions = self.cantera.model.reactions()
#         print ctReactions[0].rate
#         print ''
#         print 'Thermo before'
#         ctSpecies = self.cantera.model.species()
#         print ctSpecies[5].thermo.h(298)
        
        # Scale the thermo and kinetics of the current objects        
        for i, rv in enumerate(k_rv):
            self.scaleToKinetics(rv,self.kReactions[i])
        for i, rv in enumerate(G_rv):
            self.scaleToThermo(rv, self.gSpecies[i])
        
        # The model must be refreshed when there are any thermo changes
        # kinetics can be refreshed automatically so we don't need to recreate the Solution() object.
        if G_rv:
            self.cantera.refreshModel()
        
        # Run the cantera simulation
        allData = self.cantera.simulate()
        
        # Create a vector to hold the ModPiece output, which will be the mole fraction of the output species of interest
        output = numpy.zeros(self.outputSize)
        
        # Extract the final time point for each of the mole fractions within the outputSpeciesList
        
        for i in range(self.numConditions):
            for j in range(self.numOutputSpecies):
                speciesIndex = self.outputSpeciesIndices[j]                
                speciesGenericData = allData[i][1][2:]                
                output[i*self.numOutputSpecies+j]=speciesGenericData[speciesIndex].data[-1]
                
#         print ''
#         print 'Kinetics after'
#         ctReactions = self.cantera.model.reactions()
#         print ctReactions[0].rate
#         print ''
#         print 'Thermo after'
#         ctSpecies = self.cantera.model.species()
#         print ctSpecies[5].thermo.h(298)
        
        # Now reset the cantera object's speciesList and reactionList back to original thermo and kinetics 
        for i, thermo in enumerate(originalThermo):
            index = self.gSpecies[i]
            self.cantera.speciesList[index].thermo = thermo
            
        for i, kinetics in enumerate(originalKinetics):
            index = self.kReactions[i]
            self.cantera.reactionList[index].kinetics = kinetics
            
        return list(output)
    
            
    def scaleToKinetics(self, randomInput, reactionIndex):
        """
        This function takes a random uniform input X = Unif(-1,1) and scales the kinetics within a reaction to that value, given
        that the kinetics has a loguniform distribution where ln(k) = Unif[ln(k_min), ln(k_max)]
        
        k_sampled = 10^(randomInput * UncertaintyFactor) * k0
        
        The kinetics is permanently altered in the cantera model and must be reset to its original value after the evaluation is finished.
        """
        
        rxn = self.cantera.reactionList[reactionIndex]

        uncertaintyFactor = self.kUncertaintyFactors[reactionIndex]
        factor = randomInput*uncertaintyFactor
        
        
        # The rate is loguniform in k
        rxn.kinetics.changeRate(10**factor)
        self.cantera.modifyReactionKinetics(reactionIndex, rxn)
        

    def scaleToThermo(self,randomInput, speciesIndex):
        """
        This function takes a random normal input X = Unif(-1,1) and scales the thermodynamics for a species to that value,
        given that the thermo has a uniform distribution G = Unif(-Gmin,Gmax)
        
        G_sampled = randomInput * UncertaintyFactor + G0
        
        The thermo is permanently altered in the cantera model and must be reset to its original value after the evaluation is finished.
        """

        species = self.cantera.speciesList[speciesIndex]
        uncertaintyFactor = self.gUncertaintyFactors[speciesIndex]
        deltaH = randomInput*uncertaintyFactor*4184.0   # Convert kcal/mol to J/mol
        
        species.thermo.changeBaseEnthalpy(deltaH)
        self.cantera.modifySpeciesThermo(speciesIndex, species)
        
        
class ReactorPCEFactory:
    """
    This class uses MUQ to generate adaptive Polynomial Chaos Expansions for global uncertainty analysis in chemical reaction systems. 
    It uses RMG, Cantera, and MUQ dependencies.
    
    Methodology
        1. Set up reactor conditions and load chemical kinetic mechanism. Select desired outputs
        2. Run local uncertainty analysis
        3. Extract top N most uncertain parameters
        4. Input these set of parameters into a MUQ model class (which is a child of the ModPiece class)
        5. Create EvaluateImpl function within model class that runs simulation based on reactor conditions through Cantera
        6. Perform PCE analysis of desired outputs
    """
    
    def __init__(self, cantera, outputSpeciesList, kReactions, kUncertainty, gSpecies, gUncertainty):
        
        
        self.reactorMod = ReactorModPiece(cantera=cantera,
                            outputSpeciesList=outputSpeciesList,
                            kReactions=kReactions,
                            kUncertainty = kUncertainty,
                            gSpecies = gSpecies,
                            gUncertainty = gUncertainty,
                            )
        
        
        # Define the polynomials and quadrature rules in each dimension using a VariableCollection object. 
        # We can do this directly using the classes from libmuqUtilities. Build the PCE factory for the ReactorModPiece this way.
        # Uniform random variables used chemical kinetics uncertainty propagation uses Legendre polynomials
        # We select the Gauss-Patterson quadrature as it is recommended as the fastest in the Patrick Conrad, Youssef Marzouk paper
        
        # Select the polynomial and quadrature families 
        polyFamily = LegendrePolynomials1DRecursive()
        quadFamily = GaussPattersonQuadrature1D()
        
        # Create a random variable collection for each of the uncertain variables
        varCollection = VariableCollection()
        for rxnIndex in kReactions:
            varCollection.PushVariable("k{0}".format(rxnIndex+1), polyFamily, quadFamily)
        for speciesIndex in gSpecies:
            varCollection.PushVariable("G{0}".format(speciesIndex+1), polyFamily, quadFamily)
        
        # Initialize the PCE Factory
        self.factory = SmolyakPCEFactory(varCollection, self.reactorMod) 
        
        self.pce = None
        
    def generatePCE(self, runTime=None, startOrder=2, tolerance=None, fixedTerms=False):
        """
        Generate the PCEs adaptively. There are three methods for doing so. 
        `runTime` should be given in seconds
        Option 1: Adaptive for a pre-specified amount of time
        Option 2: Adaptively construct PCE to error tolerance
        Option 3: Used a fixed order, and (optionally) adapt later.  
        """
        
        # Also monitor the amount of time it takes
        start_time = time()
        if runTime:
            # Option 1: Adaptive for a pre-specified amount of time
            self.pce = self.factory.StartAdaptiveTimed(startOrder,runTime)
        elif tolerance:
            # Option 2: adaptively construct PCE to error tolerance
            self.pce = self.factory.StartAdaptiveToTolerance(startOrder,tolerance)
        elif fixedTerms:
            # Option 3: Used a fixed order, and (optionally) adapt later
            self.pce = self.factory.StartFixedTerms(startOrder)
        #     # Optionally adapt to tolerance later:
        #     pce = self.AdaptToTolerance(tolerance)
        else:
            raise Exception('Must have at least one chosen method')
    
        end_time = time()
        time_taken = end_time - start_time
        print 'Polynomial Chaos Expansion construction took {0:2f} seconds.'.format(time_taken)
        
    def compareOutput(self, testPoint):
        """
        Evaluate the PCEs against what the real output might give for a test point.
        testPoint is an array of all the values in terms of factor of f
        
        Returns a tuple containing the 
        (true output mole fractions, pce output mole fractions) evaluated at the test point.
        """
        
        trueOutput = self.reactorMod.Evaluate([testPoint])
        pceOutput = self.pce.Evaluate(testPoint)
        
        reactorMod = self.reactorMod
        
        for i in range(reactorMod.numConditions):
            print 'Condition {}'.format(i+1)
            print '======================================================='
            print str(reactorMod.cantera.conditions[i])
            print ''
            print 'Condition {} Mole Fractions Evaluated at Test Point'.format(i+1)
            print '========================================'
            print 'Species     True output     PCE output'
            print '========================================'
            for j, outputSpecies in enumerate(reactorMod.outputSpeciesList):
                outputIndex = i*reactorMod.numOutputSpecies+j
                print '{0:10} {1:11.2f} {2:14.2f}'.format(outputSpecies.toChemkin(),trueOutput[outputIndex],pceOutput[outputIndex])
            print ''
            
        return trueOutput, pceOutput
    
    def analyzeResults(self):
        """
        Obtain the results: the prediction mean and variance, as well as the global sensitivity indices
        Returns a tuple containing the following statistics
        
        (mean species mole fractions, variance, covariance, main sensitivity indices, total sensitivity indices)
        """
        reactorMod = self.reactorMod
        pce = self.pce
        # Compute the mean and variance for each of the uncertain parameters
        mean = numpy.array(pce.ComputeMean())
        
        var = numpy.array(pce.ComputeVariance())
        stddev = numpy.sqrt(var)
        stddev_percent = stddev/mean*100.0
        
        
        
        cov = pce.ComputeCovariance()
        # print "Covariance = ", cov
        
        
        # Extract the global sensitivity indices
        mainSens = numpy.array(pce.ComputeAllMainSensitivityIndices())
        totalSens = numpy.array(pce.ComputeAllSobolTotalSensitivityIndices())
        
        
        for i in range(reactorMod.numConditions):
            
            print 'Condition {}'.format(i+1)
            print '======================================================='
            print str(reactorMod.cantera.conditions[i])
            
            print ''
            print 'Condition {} Mole Fractions'.format(i+1)
            print '=============================================='
            print 'Species     Mean       Stddev     Stddev (%)'
            print '=============================================='
            for j, outputSpecies in enumerate(reactorMod.outputSpeciesList):
                outputIndex = i*reactorMod.numOutputSpecies+j
                print '{0:10} {1:10.3e} {2:10.3e} {3:10.3f}'.format(outputSpecies.toChemkin(), 
                                                          mean[outputIndex], 
                                                          stddev[outputIndex],
                                                          stddev_percent[outputIndex])
            print ''
        
            if reactorMod.kReactions:
                print ''
                print 'Condition {} Reaction Sensitivities'.format(i+1)
                print '==============================================================================='
                print 'Description                                             sens_main  sens_total'
                print '==============================================================================='
                for j, outputSpecies in enumerate(reactorMod.outputSpeciesList):
                    outputIndex = i*reactorMod.numOutputSpecies+j
                    for k, rxnIndex in enumerate(reactorMod.kReactions):
                        parameterIndex=k
                        description = 'dln[{0}]/dln[{1}]'.format(outputSpecies.toChemkin(),
                                                                 reactorMod.cantera.reactionList[rxnIndex].toChemkin(kinetics=False),
                                                                 )
                        print '{0:55} {1:10.3f} {2:10.3f}'.format(description,
                                                                    mainSens[outputIndex][parameterIndex],
                                                                 totalSens[outputIndex][parameterIndex],
                                                                 )
            if reactorMod.gSpecies:
                print ''
                print 'Condition {} Thermo Sensitivities'.format(i+1)
                print '==========================================================='
                print 'Description                         sens_main  sens_total'
                print '==========================================================='
                for j, outputSpecies in enumerate(reactorMod.outputSpeciesList):
                    outputIndex = i*reactorMod.numOutputSpecies+j
                    
                    for g, speciesIndex in enumerate(reactorMod.gSpecies):
                        parameterIndex = len(reactorMod.kReactions)+g
                        description = 'dln[{0}]/dlnG[{1}]'.format(outputSpecies.toChemkin(),
                                                             reactorMod.cantera.speciesList[speciesIndex].toChemkin(),)
                                                                 
                        print '{0:35} {1:10.3f} {2:10.3f}'.format(description,
                                                                 mainSens[outputIndex][parameterIndex],
                                                                 totalSens[outputIndex][parameterIndex],
                                                                 )
            print ''
            
        return mean, var, cov, mainSens, totalSens
        
