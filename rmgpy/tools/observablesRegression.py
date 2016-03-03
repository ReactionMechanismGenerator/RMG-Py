import sys
import os.path
import numpy 
from rmgpy.chemkin import loadChemkinFile
from rmgpy.tools.plot import *
from rmgpy.tools.canteraModel import *

def getNearestTime(timepoint, timeArray):
    """
    `timePoint`: the desired time point
    `timeArray`: the array of times in which to search for the nearest time to the one selected

    Returns a tuple containing (index, value of the time point closest to the timepoint desired)
    within the time array.
    """
    index = findNearest(timeArray, timepoint)
    return (index, timeArray[index])


def curvesSimilar(t1, y1, t2, y2, tol):
    """
    This function returns True if the two given curves are similar enough within tol. Otherwise returns False.

    t1: time/domain of standard curve we assume to be correct
    y1: values of standard curve, usually either temperature in (K) or log of a mol fraction
    t2: time/domain of test curve
    y2: values of test curve, usually either temperature in (K) or log of a mol fraction

    The test curve is first synchronized to the standard curve using geatNearestTime function. We then calculate the value of
    (y1-y2')^2/y1^2, giving us a normalized difference for every point. If the average value of these differences is less
    than tol, we say the curves are similar.

    We choose this criteria because it is compatible with step functions we expect to see in ignition systems.
    """
    # Make synchornized version of t2,y2 called t2sync,y2sync.
    t2sync=numpy.zeros_like(t1)
    y2sync=numpy.zeros_like(t1)
    for i, timepoint1 in enumerate(t1):
        (index, timepoint2)=getNearestTime(timepoint1, t2)
        t2sync[i]=timepoint2
        y2sync[i]=y2[index]

    # Get R^2 value equivalent:
    normalizedError=(y1-y2sync)**2/y1**2
    normalizedError=sum(normalizedError)/len(y1)

    if normalizedError > tol:
        return False
    else: 
        return True


class ObservablesTestCase:
    """
    We use this class to run regressive tests

    ======================= ==============================================================================================
    Attribute               Description
    ======================= ==============================================================================================
    `title`                 A string describing the test. For regressive tests, should be same as example's name.
    `oldDir`                A directory path containing the chem_annotated.inp and species_dictionary.txt of the old model
    `newDir`                A directory path containing the chem_annotated.inp and species_dictionary.txt of the new model
    `conditions`            A list of the :class: 'Condition' objects describing reaction conditions
    `observables`           A dictionary of observables
                            key: 'species', value: a list of smiles that correspond with species mole fraction observables
                            key: 'variable', value: a list of state variable observables, i.e. ['Temperature'] or ['Temperature','Pressure']
                            key: 'ignitionDelay', value: a tuple containing (ignition metric, yVar)
                                                         for example: ('maxDerivative','P')
                                                                      ('maxHalfConcentration', '[OH]')
                                                                      ('maxSpeciesConcentrations',['[CH2]','[O]'])
                                                        see findIgnitionDelay function for more details
    'exptData'              An array of GenericData objects
    ======================= ==============================================================================================


    """
    def __init__(self, title='', oldDir='', newDir='', observables = {}, exptData = []):
        self.title=title
        self.newDir=newDir
        self.oldDir=oldDir
        self.conditions=None
        self.exptData=exptData
        self.observables=observables

        # load the species and reactions from each model
        oldSpeciesList, oldReactionList = loadChemkinFile(os.path.join(oldDir,'chem_annotated.inp'),
                                                          os.path.join(oldDir,'species_dictionary.txt'))

        newSpeciesList, newReactionList = loadChemkinFile(os.path.join(newDir,'chem_annotated.inp'),
                                                          os.path.join(newDir,'species_dictionary.txt'))
        self.oldSim = Cantera(speciesList = oldSpeciesList,
                              reactionList = oldReactionList,
                              outputDirectory = oldDir)
        self.newSim = Cantera(speciesList = newSpeciesList,
                              reactionList = newReactionList,
                              outputDirectory = newDir)
        
        # load each chemkin file into the cantera model
        self.oldSim.loadChemkinModel(os.path.join(oldDir,'chem_annotated.inp'))
        self.newSim.loadChemkinModel(os.path.join(newDir,'chem_annotated.inp'))

    def __str__(self):
        """
        Return a string representation of this test case, using its title'.
        """
        return 'Observables Test Case: {0}'.format(self.title)

    def generateConditions(self, reactorType, reactionTime, molFracList, Tlist=None, Plist=None, Vlist=None):
        """
        Creates a list of conditions from from the lists provided. 
        
        `reactorType`: a string indicating the Cantera reactor type
        `reactionTime`: ScalarQuantity object for time
        `molFracList`: a list of dictionaries containing species smiles and their mole fraction values
        `Tlist`: ArrayQuantity object of temperatures
        `Plist`: ArrayQuantity object of pressures
        `Vlist`: ArrayQuantity object of volumes
        
        This saves all the reaction conditions into both the old and new cantera jobs.
        """
        # Store the conditions in the observables test case, for bookkeeping
        self.conditions = generateCanteraConditions(reactorType, reactionTime, molFracList, Tlist=Tlist, Plist=Plist, Vlist=Vlist)

        # Map the mole fractions dictionaries to species objects from the old and new models
        oldMolFracList = []
        newMolFracList = []

        for molFracCondition in molFracList:
            oldCondition = {}
            newCondition = {} 
            oldSpeciesDict = getRMGSpeciesFromSMILES(molFracCondition.keys(), self.oldSim.speciesList)
            newSpeciesDict = getRMGSpeciesFromSMILES(molFracCondition.keys(), self.newSim.speciesList)
            for smiles, molfrac in molFracCondition.iteritems():
                if oldSpeciesDict[smiles] is None:
                    raise Exception('SMILES {0} was not found in the old model!'.format(smiles))
                if newSpeciesDict[smiles] is None:
                    raise Exception('SMILES {0} was not found in the new model!'.format(smiles))

                oldCondition[oldSpeciesDict[smiles]] = molfrac
                newCondition[newSpeciesDict[smiles]] = molfrac
            oldMolFracList.append(oldCondition)
            newMolFracList.append(newCondition)
        
        # Generate the conditions in each simulation
        self.oldSim.generateConditions(reactorType, reactionTime, oldMolFracList, Tlist=Tlist, Plist=Plist, Vlist=Vlist)
        self.newSim.generateConditions(reactorType, reactionTime, newMolFracList, Tlist=Tlist, Plist=Plist, Vlist=Vlist)

    def compare(self, plot=False):
        """
        Compare the old and new model
        `plot`: if set to True, it will comparison plots of the two models comparing their species.

        Returns a list of variables failed in a list of tuples in the format:
        
        (CanteraCondition, variable label, variableOld, variableNew)

        """
        # Ignore Inerts
        inertList = ['[Ar]','[He]','[N#N]','[Ne]']

        oldConditionData, newConditionData = self.runSimulations()

        conditionsBroken=[]
        variablesFailed=[]
        
        print ''
        print '{0} Comparison'.format(self)
        print '================'
        # Check the species profile observables
        if 'species' in self.observables:
            oldSpeciesDict = getRMGSpeciesFromSMILES(self.observables['species'], self.oldSim.speciesList)
            newSpeciesDict = getRMGSpeciesFromSMILES(self.observables['species'], self.newSim.speciesList)
        
        # Check state variable observables 
        implementedVariables = ['temperature','pressure']
        if 'variable' in self.observables:
            for item in self.observables['variable']:
                if item.lower() not in implementedVariables:
                    print 'Observable variable {0} not yet implemented'.format(item)
                    
        print ''
        print 'The following observables did not match:'
        print ''
        for i in range(len(oldConditionData)):
            timeOld, dataListOld = oldConditionData[i]
            timeNew, dataListNew = newConditionData[i]

            # Compare species observables
            if 'species' in self.observables:
                for smiles in self.observables['species']:
                    
                    fail = False
                    oldRmgSpecies = oldSpeciesDict[smiles]
                    newRmgSpecies = newSpeciesDict[smiles]
                    
                    if oldRmgSpecies:
                        variableOld = next((data for data in dataListOld if data.species == oldRmgSpecies), None)
                    else:
                        print 'No RMG species found for observable species {0} in old model.'.format(smiles)
                        fail = True
                    if newRmgSpecies:
                        variableNew = next((data for data in dataListNew if data.species == newRmgSpecies), None)
                    else:
                        print 'No RMG species found for observable species {0} in new model.'.format(smiles)
                        fail = True
                    
                    if fail is False:
                        if not curvesSimilar(timeOld.data, variableOld.data, timeNew.data, variableNew.data, 0.05):
                            fail = True
                            
                        # Try plotting only when species are found in both models
                        if plot:
                            oldSpeciesPlot = SimulationPlot(xVar=timeOld, yVar=variableOld)
                            newSpeciesPlot = SimulationPlot(xVar=timeNew, yVar=variableNew)
                            oldSpeciesPlot.comparePlot(newSpeciesPlot,
                                                       title='Observable Species {0} Comparison'.format(smiles),
                                                       ylabel='Mole Fraction',
                                                       filename='condition_{0}_species_{1}.png'.format(i+1,smiles))
                    
                    # Append to failed variables or conditions if this test failed
                    if fail:
                        if i not in conditionsBroken: conditionsBroken.append(i)
                        print "Observable species {0} does not match between old model {1} and \
new model {2} in condition {3:d}.".format(smiles,
                                           variableOld.label, 
                                           variableNew.label,
                                           i+1)
                        variablesFailed.append((self.conditions[i], smiles, variableOld, variableNew))
                    
            
            # Compare state variable observables
            if 'variable' in self.observables:
                for varName in self.observables['variable']:
                    variableOld = next((data for data in dataListOld if data.label == varName), None)
                    variableNew = next((data for data in dataListNew if data.label == varName), None)
                    if not curvesSimilar(timeOld.data, variableOld.data, timeNew.data, variableNew.data, 0.05):
                        if i not in conditionsBroken: conditionsBroken.append(i)
                        print "Observable variable {0} does not match between old model and \
new model in condition {1:d}.".format(variableOld.label, i+1)
                        variablesFailed.append((self.conditions[i], varName, variableOld, variableNew))
                    
                    if plot:
                        oldVarPlot = GenericPlot(xVar=timeOld, yVar=variableOld)
                        newVarPlot = GenericPlot(xVar=timeNew, yVar=variableNew)
                        oldVarPlot.comparePlot(newSpeciesPlot,
                                                   title='Observable Variable {0} Comparison'.format(varName),
                                                   filename='condition_{0}_variable_{1}.png'.format(i+1, varName))
                        
            # Compare ignition delay observables
            if 'ignitionDelay' in self.observables:
                print 'Ignition delay observable comparison not implemented yet.'
                
                
        
        print ''
        print 'The following reaction conditions were broken:'
        print ''
        for index in conditionsBroken:
            print "Condition {0:d}:".format(index+1)
            print str(self.conditions[index])
            print ''

        return variablesFailed

    def runSimulations(self):
        """
        Run a selection of conditions in Cantera and return
        generic data objects containing the time, pressure, temperature,
        and mole fractions from the simulations.

        Returns (oldConditionData, newConditionData)
        where conditionData is a list of of tuples: (time, dataList) for each condition in the same order as conditions
        time is a GenericData object which gives the time domain for each profile
        dataList is a list of GenericData objects for the temperature, profile, and mole fraction of major species
        """
        oldConditionData = self.oldSim.simulate()
        newConditionData = self.newSim.simulate()
        return (oldConditionData, newConditionData)
