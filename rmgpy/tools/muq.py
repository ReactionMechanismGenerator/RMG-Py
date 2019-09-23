#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

from __future__ import division

import copy
import logging
from time import time

import numpy as np

try:
    from libmuqModelling import ModPiece
    from libmuqUtilities import LegendrePolynomials1DRecursive, GaussPattersonQuadrature1D, VariableCollection
    from libmuqApproximation import SmolyakPCEFactory
except:
    ModPiece = object
    print('Could not import MUQ. Please check that it is installed correctly before using the global uncertainty modules.')

# You must install the MUQ library before using this.  Add the folder containing
# libmuqUtilities.so, libmuqApproximation.so, etc to your $PYTHONPATH
# For linux users, you can install via 'conda install -c rmg muq' to your environment
# and add the ~/anaconda/envs/your_env/lib folder to your $PYTHONPATH


class ReactorModPiece(ModPiece):
    def __init__(self, cantera, outputSpeciesList, kParams, kUncertainty, gParams, gUncertainty, correlated=False, logx=True):
        """
        ======================= ====================================================
        Attribute               Description
        ======================= ====================================================
        `cantera`               A Cantera() object containing CanteraConditions and initialized species and reactions
        `outputSpeciesList`     A list of Species() objects corresponding to the desired observables for uncertainty analysis
        `kParams`               Uncorrelated: A list of indices of the Reaction() objects in cantera.reactionList corresponding to the uncertain input rate coefficients
                                Correlated: this is a list of strings of the uncertain kinetic parameter sources to be propagated in the model. i.e. 'H_Abstraction CHO/Oa'
        `kUncertainty`          Uncorrelated: A list of uncertainties dlnk corresponding to the reactions in kReactions
                                Correlated: A list of dictionaries corresponding to each reaction's partial uncertainties with respect to various kinetic sources.
                                            This is the output object from the Uncertainty.kineticInputUncertainties 
        `gParams `              Uncorrelated: A list of indices for the Species() objects in cantera.speciesList corresponding to the uncertain input free energies of individual species
                                Correlated: A list of strings corresponding to the uncertain thermo sources to be propagated.  i.e. 'Group(group) C=O'
        `gUncertainty`          Uncorrelated: A list of uncertainties dG corresponding to the species in gSpecies in units of kcal/mol
                                Correlated: A list of dictionaries corresponding to each specie's partial uncertainties with respect to various thermo sources.
                                            This is the output object from the Uncertainty.thermoInputUncertainties 
        `correlated`            A flag set to either `True` or `False` depending on whether the uncertainties are correlated or not
                                If correlated, the correlated uncertainties will propagate and affect multiple species thermo and reaction kinetics parameters
        
        `affectedReactions`      Used only in the correlated case, this is a list of the indices of the affected rxns for kParams
        `affectedSpecies`       Used only in the correlated case, this is a list of the indices of the affected species for gParams
        `logx`                  A flag to set whether to use mole fraction or ln(mole fraction) as the output variable
        ============================================================================
        """
        self.cantera = cantera
        self.outputSpeciesList = outputSpeciesList
        self.outputSpeciesIndices = [cantera.speciesList.index(outputSpecies) for outputSpecies in outputSpeciesList]
        self.correlated = correlated
        self.logx = logx

        # The size of the uncertain inputs: [parameters affecting k, parameters affecting free energy G]         
        self.inputSize = [len(kParams) + len(gParams)]

        if not self.correlated:
            # Convert input indices (RMG indices) into species list indices
            new_kParams = []
            for ind in kParams:
                for i, rxn in enumerate(cantera.reactionList):
                    if rxn.index == ind:
                        new_kParams.append(i)
                        logging.debug('Replacing reaction index {0} with {1} for reaction {2!s}'.format(ind, i, rxn))
                        break
                else:
                    raise ValueError('Could not find requested index {0} in reaction list.'.format(ind))
            new_gParams = []
            for ind in gParams:
                for i, spc in enumerate(cantera.speciesList):
                    if spc.index == ind:
                        new_gParams.append(i)
                        logging.debug('Replacing species index {0} with {1} for species {2!s}'.format(ind, i, spc))
                        break
                else:
                    raise ValueError('Could not find requested index {0} in species list.'.format(ind))
            kParams = new_kParams
            gParams = new_gParams

        # for uncorrelated case, these are indices of reactions
        # for correlated case, this is the list of labels for the uncertain parameters to be perturbed, i.e. 'H_Abstraction CHO/Oa'
        self.kParams = kParams
        # for uncorrelated case, these are indices of the species
        # for correlated case, this is the list of labels for the uncertain thermo parameters to be perturbed, i.e. 'Group(ring) cyclohexane'
        self.gParams = gParams

        if not self.correlated:
            kUncertaintyFactors = [val * np.sqrt(3) / np.log(10) for val in kUncertainty]
            self.kUncertaintyFactors = {}
            for i, rxnIndex in enumerate(kParams):
                self.kUncertaintyFactors[rxnIndex] = kUncertaintyFactors[rxnIndex]
                logging.debug('For {0}, set uncertainty factor to {1}'.format(cantera.reactionList[rxnIndex], kUncertaintyFactors[rxnIndex]))
            
            gUncertaintyFactors = [val * np.sqrt(3) for val in gUncertainty]
            self.gUncertaintyFactors = {}
            for i, spcIndex in enumerate(gParams):
                self.gUncertaintyFactors[spcIndex] = gUncertaintyFactors[spcIndex]
                logging.debug('For {0}, set uncertainty factor to {1}'.format(cantera.speciesList[spcIndex], gUncertaintyFactors[spcIndex]))
            
        else:
            # In the correlated case, keep track of which reactions and species each 
            # uncertain parameter affects 
            self.kUncertaintyFactors = {}
            affectedReactions = set()
            for kParam in kParams:
                rxnPartialUncertainty = []
                for rxnIndex, rxnInputDict in enumerate(kUncertainty):
                    if kParam in rxnInputDict:
                        affectedReactions.add(rxnIndex)
                        uncertaintyFactor = rxnInputDict[kParam] * np.sqrt(3) / np.log(10)
                        # If this parameter string contributes to the reaction, add this reaction index
                        rxnPartialUncertainty.append((rxnIndex, uncertaintyFactor))
                self.kUncertaintyFactors[kParam] = rxnPartialUncertainty  # list of indices of the reactions affected
            self.affectedReactions = list(affectedReactions)

            self.gUncertaintyFactors = {}
            affectedSpecies = set()
            for gParam in gParams:
                spcPartialUncertainty = []
                for spcIndex, spcInputDict in enumerate(gUncertainty):
                    if gParam in spcInputDict:
                        affectedSpecies.add(spcIndex)
                        uncertaintyFactor = spcInputDict[gParam] * np.sqrt(3)
                        # If this parameter contributes to the species thermo, add this species index
                        spcPartialUncertainty.append((spcIndex, uncertaintyFactor))
                self.gUncertaintyFactors[gParam] = spcPartialUncertainty
            self.affectedSpecies = list(affectedSpecies)

        # The size of the vector corresponding to the outputs to be analyzed for uncertainty analysis
        # is equal to the number of cantera conditions involved multiplied by the number of desired observables
        outputSize = len(cantera.conditions) * len(outputSpeciesList)

        self.numOutputSpecies = len(outputSpeciesList)
        self.numConditions = len(cantera.conditions)

        # Initialize the ModPiece with some input and output size specifications
        ModPiece.__init__(self,
                          self.inputSize,
                          outputSize,
                          False,  # No GradientImpl
                          False,  # No JacobianImpl
                          False,  # No HessianImpl
                          False,  # Not random)
                          False,  # Not random)
                          )

    def EvaluateImpl(self, ins):
        """
        Evaluate the desired output mole fractions based on a set of inputs ins = [[k_rv], [G_rv]] which contains the 
        random normal variables attributed to the uncertain kinetics and free energy parameters, respectively.
        
        The output returned contains [Condition1_outputMoleFraction1, Condition1_outputMoleFraction2, Condition2_output.... ConditionN_output...]
        """
        assert len(ins[0]) == self.inputSize[0], "Number of inputs matches number of uncertain parameters"

        k_rv = ins[0][0:len(self.kParams)]
        G_rv = ins[0][len(self.kParams):]

        ## Check that the number of inputs is correct
        # assert len(k_rv) == len(self.kParams), "Number of inputs matches number of kParams"
        # assert len(G_rv) == len(self.gParams), "Number of inputs matches number of gParams"

        if not self.correlated:
            # Make deepcopies of the thermo and kinetics so as to not modify the originals in the speciesList and reactionList
            originalThermo = [copy.deepcopy(self.cantera.speciesList[index].thermo) for index in self.gParams]
            originalKinetics = [copy.deepcopy(self.cantera.reactionList[index].kinetics) for index in self.kParams]

            # print ''
            # print 'Kinetics before'
            # ctReactions = self.cantera.model.reactions()
            # print ctReactions[0].rate
            # print ''
            # print 'Thermo before'
            # ctSpecies = self.cantera.model.species()
            # print ctSpecies[5].thermo.h(298)

            # Scale the thermo and kinetics of the current objects        
            for i, rv in enumerate(k_rv):
                rxnIndex = self.kParams[i]
                self.scaleToKinetics(rv, self.kUncertaintyFactors[rxnIndex], rxnIndex)
            for i, rv in enumerate(G_rv):
                spcIndex = self.gParams[i]
                self.scaleToThermo(rv, self.gUncertaintyFactors[spcIndex], spcIndex)

        else:
            # Make deepcopies of the thermo and kinetics so as to not modify the originals in the speciesList and reactionList
            originalKinetics = [copy.deepcopy(self.cantera.reactionList[index].kinetics) for index in
                                self.affectedReactions]
            originalThermo = [copy.deepcopy(self.cantera.speciesList[index].thermo) for index in self.affectedSpecies]

            mappedReactionScaling = {index: 0.0 for index in self.affectedReactions}
            mappedSpeciesScaling = {index: 0.0 for index in self.affectedSpecies}

            for i, kParam in enumerate(self.kParams):
                for rxnIndex, uncertaintyFactor in self.kUncertaintyFactors[kParam]:
                    mappedReactionScaling[rxnIndex] += k_rv[i] * uncertaintyFactor

            for i, gParam in enumerate(self.gParams):
                for spcIndex, uncertaintyFactor in self.gUncertaintyFactors[gParam]:
                    mappedSpeciesScaling[spcIndex] += G_rv[i] * uncertaintyFactor

            for rxnIndex, uncertaintyFactor in mappedReactionScaling.items():
                self.scaleToKinetics(1.0, uncertaintyFactor, rxnIndex)
            for spcIndex, uncertaintyFactor in mappedSpeciesScaling.items():
                self.scaleToThermo(1.0, uncertaintyFactor, spcIndex)

        # The model must be refreshed when there are any thermo changes
        # kinetics can be refreshed automatically so we don't need to recreate the Solution() object.
        if G_rv:
            self.cantera.refreshModel()

        # Run the cantera simulation
        allData = self.cantera.simulate()

        # Create a vector to hold the ModPiece output, which will be the mole fraction of the output species of interest
        output = np.zeros(self.outputSize)

        # Extract the final time point for each of the mole fractions within the outputSpeciesList

        for i in range(self.numConditions):
            for j in range(self.numOutputSpecies):
                speciesIndex = self.outputSpeciesIndices[j]
                speciesGenericData = allData[i][1][2:]
                if self.logx:
                    output[i * self.numOutputSpecies + j] = np.log(speciesGenericData[speciesIndex].data[-1])
                else:
                    output[i * self.numOutputSpecies + j] = speciesGenericData[speciesIndex].data[-1]

        # print ''
        # print 'Kinetics after'
        # ctReactions = self.cantera.model.reactions()
        # print ctReactions[0].rate
        # print ''
        # print 'Thermo after'
        # ctSpecies = self.cantera.model.species()
        # print ctSpecies[5].thermo.h(298)

        if not self.correlated:
            # Now reset the cantera object's speciesList and reactionList back to original thermo and kinetics 
            for i, thermo in enumerate(originalThermo):
                index = self.gParams[i]
                self.cantera.speciesList[index].thermo = thermo

            for i, kinetics in enumerate(originalKinetics):
                index = self.kParams[i]
                self.cantera.reactionList[index].kinetics = kinetics
        else:
            for i, thermo in enumerate(originalThermo):
                index = self.affectedSpecies[i]
                self.cantera.speciesList[index].thermo = thermo
            for i, kinetics in enumerate(originalKinetics):
                index = self.affectedReactions[i]
                self.cantera.reactionList[index].kinetics = kinetics

        return list(output)

    def scaleToKinetics(self, randomInput, uncertaintyFactor, reactionIndex):
        """
        This function takes a random uniform input X = Unif(-1,1) and scales the kinetics within a reaction to that value, given
        that the kinetics has a loguniform distribution where ln(k) = Unif[ln(k_min), ln(k_max)]
        
        k_sampled = 10^(randomInput * UncertaintyFactor) * k0
        
        The kinetics is permanently altered in the cantera model and must be reset to its original value after the evaluation is finished.
        """

        rxn = self.cantera.reactionList[reactionIndex]
        factor = randomInput * uncertaintyFactor

        # The rate is loguniform in k
        rxn.kinetics.change_rate(10 ** factor)
        self.cantera.modifyReactionKinetics(reactionIndex, rxn)

    def scaleToThermo(self, randomInput, uncertaintyFactor, speciesIndex):
        """
        This function takes a random normal input X = Unif(-1,1) and scales the thermodynamics for a species to that value,
        given that the thermo has a uniform distribution G = Unif(-Gmin,Gmax)
        
        G_sampled = randomInput * UncertaintyFactor + G0
        
        The thermo is permanently altered in the cantera model and must be reset to its original value after the evaluation is finished.
        """

        species = self.cantera.speciesList[speciesIndex]
        deltaH = randomInput * uncertaintyFactor * 4184.0  # Convert kcal/mol to J/mol

        species.thermo.changeBaseEnthalpy(deltaH)
        self.cantera.modifySpeciesThermo(speciesIndex, species, useChemkinIdentifier=True)


class ReactorPCEFactory(object):
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

    def __init__(self, cantera, outputSpeciesList, kParams, kUncertainty, gParams, gUncertainty, correlated=False, logx=True):

        self.reactorMod = ReactorModPiece(cantera=cantera,
                                          outputSpeciesList=outputSpeciesList,
                                          kParams=kParams,
                                          kUncertainty=kUncertainty,
                                          gParams=gParams,
                                          gUncertainty=gUncertainty,
                                          correlated=correlated,
                                          logx=logx,
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
        for i, rxnIndex in enumerate(kParams):
            varCollection.PushVariable("k{0}".format(i + 1), polyFamily, quadFamily)
        for i, speciesIndex in enumerate(gParams):
            varCollection.PushVariable("G{0}".format(i + 1), polyFamily, quadFamily)

        # Initialize the PCE Factory
        self.factory = SmolyakPCEFactory(varCollection, self.reactorMod)

        self.pce = None
        self.logx = logx

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
            self.pce = self.factory.StartAdaptiveTimed(startOrder, runTime)
        elif tolerance:
            # Option 2: adaptively construct PCE to error tolerance
            self.pce = self.factory.StartAdaptiveToTolerance(startOrder, tolerance)
        elif fixedTerms:
            # Option 3: Used a fixed order, and (optionally) adapt later
            self.pce = self.factory.StartFixedTerms(startOrder)
            # # Optionally adapt to tolerance later:
            # pce = self.AdaptToTolerance(tolerance)
        else:
            raise Exception('Must have at least one chosen method')

        end_time = time()
        time_taken = end_time - start_time
        logging.info('Polynomial Chaos Expansion construction took {0:2f} seconds.'.format(time_taken))

    def compareOutput(self, testPoint, log=True):
        """
        Evaluate the PCEs against what the real output might give for a test point.
        testPoint is an array of all the values in terms of factor of f
        
        Returns a tuple containing the 
        (true output mole fractions, pce output mole fractions) evaluated at the test point.
        """

        trueOutput = self.reactorMod.Evaluate([testPoint])
        pceOutput = self.pce.Evaluate(testPoint)

        reactorMod = self.reactorMod

        output = ''
        for i in range(reactorMod.numConditions):
            output += """============================================================
Condition {0}
------------------------------------------------------------
{1!s}
============================================================
Condition {0} {2}Mole Fractions Evaluated at Test Point
------------------------------------------------------------
Species                      True Output          PCE Output
------------------------------------------------------------
""".format(i + 1, reactorMod.cantera.conditions[i], 'Log ' if self.logx else '')

            for j, outputSpecies in enumerate(reactorMod.outputSpeciesList):
                outputIndex = i * reactorMod.numOutputSpecies + j
                output += '{0:<20}{1:>20.3f}{2:>20.3f}\n'.format(outputSpecies.to_chemkin(), trueOutput[outputIndex], pceOutput[outputIndex])
            output += '============================================================\n'

        if log:
            logging.info(output)
        else:
            print(output)

        return trueOutput, pceOutput

    def analyzeResults(self, log=True):
        """
        Obtain the results: the prediction mean and variance, as well as the global sensitivity indices
        Returns a tuple containing the following statistics
        
        (mean species mole fractions, variance, covariance, main sensitivity indices, total sensitivity indices)
        """
        reactorMod = self.reactorMod
        pce = self.pce
        # Compute the mean and variance for each of the uncertain parameters
        mean = np.array(pce.ComputeMean())

        var = np.array(pce.ComputeVariance())
        stddev = np.sqrt(var)
        stddev_percent = stddev / mean * 100.0

        cov = pce.ComputeCovariance()

        # Extract the global sensitivity indices
        mainSens = np.array(pce.ComputeAllMainSensitivityIndices())
        totalSens = np.array(pce.ComputeAllSobolTotalSensitivityIndices())

        output = ''
        for i in range(reactorMod.numConditions):
            output += """============================================================
Condition {0}
------------------------------------------------------------
{1!s}
============================================================
Condition {0} {2}Mole Fractions
------------------------------------------------------------
Species                   Mean         Stddev     Stddev (%)
------------------------------------------------------------
""".format(i + 1, reactorMod.cantera.conditions[i], 'Log ' if self.logx else '')

            for j, outputSpecies in enumerate(reactorMod.outputSpeciesList):
                outputIndex = i * reactorMod.numOutputSpecies + j
                output += '{0:<15}{1:>15.3e}{2:>15.3e}{3:>15.3f}\n'.format(outputSpecies.to_chemkin(),
                                                                           mean[outputIndex],
                                                                           stddev[outputIndex],
                                                                           stddev_percent[outputIndex])
            output += '============================================================\n\n'

            if reactorMod.kParams:
                output += """====================================================================================================
Condition {0} Reaction Rate Sensitivity Indices
----------------------------------------------------------------------------------------------------
Description                                                                 sens_main     sens_total
""".format(i + 1)

                for j, outputSpecies in enumerate(reactorMod.outputSpeciesList):
                    output += '----------------------------------------------------------------------------------------------------\n'
                    outputIndex = i * reactorMod.numOutputSpecies + j
                    for k, descriptor in enumerate(reactorMod.kParams):
                        parameterIndex = k
                        if not reactorMod.correlated:
                            description = 'dln[{0}]/dln[{1}]'.format(outputSpecies.to_chemkin(),
                                                                     reactorMod.cantera.reactionList[descriptor].to_chemkin(kinetics=False))
                        else:
                            description = 'dln[{0}]/dln[{1}]'.format(outputSpecies.to_chemkin(), descriptor)

                        output += '{0:<70}{1:>14.3f}%{2:>14.3f}%\n'.format(description,
                                                                           100 * mainSens[outputIndex][parameterIndex],
                                                                           100 * totalSens[outputIndex][parameterIndex])
                output += '====================================================================================================\n\n'

            if reactorMod.gParams:
                output += """====================================================================================================
Condition {0} Thermochemistry Sensitivity Indices
----------------------------------------------------------------------------------------------------
Description                                                                 sens_main     sens_total
""".format(i + 1)
                for j, outputSpecies in enumerate(reactorMod.outputSpeciesList):
                    output += '----------------------------------------------------------------------------------------------------\n'
                    outputIndex = i * reactorMod.numOutputSpecies + j
                    for g, descriptor in enumerate(reactorMod.gParams):
                        parameterIndex = len(reactorMod.kParams) + g
                        if not reactorMod.correlated:
                            description = 'dln[{0}]/dG[{1}]'.format(outputSpecies.to_chemkin(),
                                                                    reactorMod.cantera.speciesList[descriptor].to_chemkin())
                        else:
                            description = 'dln[{0}]/dG[{1}]'.format(outputSpecies.to_chemkin(), descriptor)

                        output += '{0:<70}{1:>14.3f}%{2:>14.3f}%\n'.format(description,
                                                                           100 * mainSens[outputIndex][parameterIndex],
                                                                           100 * totalSens[outputIndex][parameterIndex])
                output += '====================================================================================================\n\n'

        if log:
            logging.info(output)
        else:
            print(output)

        return mean, var, cov, mainSens, totalSens
