#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

import sys
import os.path
import numpy
import csv
from rmgpy.chemkin import loadChemkinFile
from rmgpy.tools.plot import GenericData, GenericPlot, SimulationPlot, findNearest
from rmgpy.tools.canteraModel import Cantera, generateCanteraConditions, getRMGSpeciesFromUserSpecies


def curvesSimilar(t1, y1, t2, y2, tol, method='maxDiff', ignore1=0.0, ignore2=0.0):
    """
    This function returns True if the two given curves are similar enough within tol. Otherwise returns False.

    t1: List for time/domain of standard curve we assume to be correct (curve 1)
    y1: List for values of standard curve, usually either temperature in (K) or log of a mol fraction
    t2: List for time/domain of test curve (curve 2)
    y2: List for values of test curve, usually either temperature in (K) or log of a mol fraction
    tol: maximum percent difference that curve 2 can be from curve 1 before we return False
    method: can be 'maxDiff' or 'avgDiff', describing method for determining similarity (see below)
    ignore1: First domain point defining interval to ignore (see below)
    ignore2: Second domain point defining interval to ignore (see below)

    For the maxDiff method:
    The test curve is first synchronized to the standard curve using geatNearestTime function. We then calculate the
    value of abs((y1-y2')/y1), giving us a normalized difference for every point. If any point has a greater percent
    difference than tol, then we say the curves are not similar

    For the avgDiff method:
    The test curve is first synchronized to the standard curve using geatNearestTime function. We then calculate the
    value of abs((y1-y2')/y1), giving us a normalized difference for every point. If the average value of these
    differences is greater than tol, we say the curves are not similar.

    Additionally, we can add a domain interval where we do not test the difference between the two curves. This could
    be useful in the case of ignition, where we expect a large difference in y-values between and after ignition. If
    the ignition event did not occur at exactly the same time, we would could potentially fail the above criteria when
    the curves are still quite similar. In this case you would set ignore1 to the ignition delay from curve 1 and
    ignore2 to the ignition delay of curve 2. (And rely on the compareIgnitionDelay function to test if that is similar)
    """
    # Make synchornized version of t2,y2 called t2sync,y2sync.
    t2sync=numpy.zeros_like(t1)
    y2sync=numpy.zeros_like(t1)
    for i, timepoint1 in enumerate(t1):
        time_index = findNearest(t2, timepoint1)
        t2sync[i]=t2[time_index]
        y2sync[i]=y2[time_index]

    #Determine which of the ignore intervals is larger
    if ignore2 < ignore1:
        sign=1.0
    #If ignore 2 is smaller than ignore1, then flip the sign on both so we can do the same operation to see if
    #ignore1<sign*t<ignore2 for both cases
    else:
        ignore1=-ignore1
        ignore2=-ignore2
        sign=-1.0

    if method=='maxDiff':
        for i, timepoint1 in enumerate(t1):
            if ignore1<sign*timepoint1 and ignore2>sign*timepoint1: pass
            else:
                if abs((y1[i]-y2sync[i])/y1[i])>tol: return False
        else: return True

    elif method=='avgDiff':
        #check to see that we are not in ignore interval
        totalError=0.0
        ignoredTimePoints=0.0
        for i, timepoint1 in enumerate(t1):
            if ignore1<sign*timepoint1 and ignore2>sign*timepoint1:
                ignoredTimePoints+=1
            else:totalError+=abs((y1[i]-y2sync[i])/y1[i])

        normalizedError=totalError/(len(y1)-ignoredTimePoints)

        if normalizedError > tol:
            return False
        else:
            return True
    else: raise Exception("Method {0} for comparing curves is not supported".format(method))

class ObservablesTestCase:
    """
    We use this class to run regressive tests

    ======================= ==============================================================================================
    Attribute               Description
    ======================= ==============================================================================================
    `title`                 A string describing the test. For regressive tests, should be same as example's name.
    `oldDir`                A directory path containing the chem_annotated.inp and species_dictionary.txt of the old model
    `newDir`                A directory path containing the chem_annotated.inp and species_dictionary.txt of the new model
    `conditions`            A list of the :class: 'CanteraCondition' objects describing reaction conditions
    `observables`           A dictionary of observables
                            key: 'species', value: a list of the "class" 'Species' that correspond with species mole fraction observables
                            key: 'variable', value: a list of state variable observables, i.e. ['Temperature'] or ['Temperature','Pressure']
                            key: 'ignitionDelay', value: a tuple containing (ignition metric, yVar)
                                                         for example: ('maxDerivative','P')
                                                                      ('maxHalfConcentration', '[OH]')
                                                                      ('maxSpeciesConcentrations',['[CH2]','[O]'])
                                                        see findIgnitionDelay function for more details
    'exptData'              An array of GenericData objects
    'ck2tci'                Indicates whether to convert chemkin to cti mechanism.  If set to False, RMG will convert the species
                            and reaction objects to Cantera objects internally
    ======================= ==============================================================================================


    """
    def __init__(self, title='', oldDir='', newDir='', observables = None, exptData = None, ck2cti=True):
        self.title=title
        self.newDir=newDir
        self.oldDir=oldDir
        self.conditions=None
        self.exptData=exptData if exptData else []
        self.observables=observables if observables else {}

        # Detect if the transport file exists
        oldTransportPath = None
        if os.path.exists(os.path.join(oldDir,'tran.dat')):
            oldTransportPath = os.path.join(oldDir,'tran.dat')
        newTransportPath = None
        if os.path.exists(os.path.join(newDir,'tran.dat')):
            newTransportPath = os.path.join(newDir,'tran.dat')

        # load the species and reactions from each model
        oldSpeciesList, oldReactionList = loadChemkinFile(os.path.join(oldDir,'chem_annotated.inp'),
                                                          os.path.join(oldDir,'species_dictionary.txt'),
                                                          oldTransportPath)

        newSpeciesList, newReactionList = loadChemkinFile(os.path.join(newDir,'chem_annotated.inp'),
                                                          os.path.join(newDir,'species_dictionary.txt'),
                                                          newTransportPath)

        self.oldSim = Cantera(speciesList = oldSpeciesList,
                              reactionList = oldReactionList,
                              outputDirectory = oldDir)
        self.newSim = Cantera(speciesList = newSpeciesList,
                              reactionList = newReactionList,
                              outputDirectory = newDir)
        
        # load each chemkin file into the cantera model
        if not ck2cti:
            self.oldSim.loadModel()
            self.newSim.loadModel()
        else:
            self.oldSim.loadChemkinModel(os.path.join(oldDir,'chem_annotated.inp'), transportFile=oldTransportPath, quiet=True)
            self.newSim.loadChemkinModel(os.path.join(newDir,'chem_annotated.inp'), transportFile=newTransportPath, quiet=True)

    def __str__(self):
        """
        Return a string representation of this test case, using its title'.
        """
        return 'Observables Test Case: {0}'.format(self.title)

    def generateConditions(self, reactorTypeList, reactionTimeList, molFracList, Tlist=None, Plist=None, Vlist=None):
        """
        Creates a list of conditions from from the lists provided. 
        
        ======================= ====================================================
        Argument                Description
        ======================= ====================================================
        `reactorTypeList`        A list of strings of the cantera reactor type. List of supported types below:
            IdealGasReactor: A constant volume, zero-dimensional reactor for ideal gas mixtures
            IdealGasConstPressureReactor: A homogeneous, constant pressure, zero-dimensional reactor for ideal gas mixtures

        `reactionTimeList`      A tuple object giving the ([list of reaction times], units)
        `molFracList`           A list of molfrac dictionaries with species object keys
                               and mole fraction values
        To specify the system for an ideal gas, you must define 2 of the following 3 parameters:
        `T0List`                A tuple giving the ([list of initial temperatures], units)
        'P0List'                A tuple giving the ([list of initial pressures], units)
        'V0List'                A tuple giving the ([list of initial specific volumes], units)
        
        This saves all the reaction conditions into both the old and new cantera jobs.
        """
        # Store the conditions in the observables test case, for bookkeeping
        self.conditions = generateCanteraConditions(reactorTypeList, reactionTimeList, molFracList, Tlist=Tlist, Plist=Plist, Vlist=Vlist)

        # Map the mole fractions dictionaries to species objects from the old and new models
        oldMolFracList = []
        newMolFracList = []

        for molFracCondition in molFracList:
            oldCondition = {}
            newCondition = {} 
            oldSpeciesDict = getRMGSpeciesFromUserSpecies(molFracCondition.keys(), self.oldSim.speciesList)
            newSpeciesDict = getRMGSpeciesFromUserSpecies(molFracCondition.keys(), self.newSim.speciesList)
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
        self.oldSim.generateConditions(reactorTypeList, reactionTimeList, oldMolFracList, Tlist=Tlist, Plist=Plist, Vlist=Vlist)
        self.newSim.generateConditions(reactorTypeList, reactionTimeList, newMolFracList, Tlist=Tlist, Plist=Plist, Vlist=Vlist)

    def compare(self, tol, method="maxDiff", plot=False):
        """
        Compare the old and new model
        'tol':  average error acceptable between old and new model for variables
        `plot`: if set to True, it will comparison plots of the two models comparing their species.

        Returns a list of variables failed in a list of tuples in the format:
        
        (CanteraCondition, variable label, variableOld, variableNew)

        """
        # Ignore Inerts
        inertList = ['[Ar]','[He]','[N#N]','[Ne]']

        #Dictionary that gives string to be printed depending on method of comparison
        methodString={"maxDiff": "at least once",
                      "avgDiff": "on average",
        }

        oldConditionData, newConditionData = self.runSimulations()

        conditionsBroken=[]
        variablesFailed=[]
        
        print ''
        print '{0} Comparison'.format(self)
        print '================'
        # Check the species profile observables
        if 'species' in self.observables:
            oldSpeciesDict = getRMGSpeciesFromUserSpecies(self.observables['species'], self.oldSim.speciesList)
            newSpeciesDict = getRMGSpeciesFromUserSpecies(self.observables['species'], self.newSim.speciesList)
        
        # Check state variable observables 
        implementedVariables = ['temperature','pressure']
        if 'variable' in self.observables:
            for item in self.observables['variable']:
                if item.lower() not in implementedVariables:
                    print 'Observable variable {0} not yet implemented'.format(item)
                    
        failHeader='\nThe following observables did not match:\n'
        failHeaderPrinted=False
        for i in range(len(oldConditionData)):
            timeOld, dataListOld, reactionSensitivityDataOld = oldConditionData[i]
            timeNew, dataListNew, reactionSensitivityDataOld = newConditionData[i]

            # Compare species observables
            if 'species' in self.observables:
                smilesList=[] #This is to make sure we don't have species with duplicate smiles
                multiplicityList=['','(S)','(D)','(T)','(Q)'] #list ot add multiplcity
                for species in self.observables['species']:

                    smiles=species.molecule[0].toSMILES() #For purpose of naming the plot only
                    if smiles in smilesList: smiles=smiles+multiplicityList[species.molecule[0].multiplicity]
                    smilesList.append(smiles)
                    
                    fail = False
                    oldRmgSpecies = oldSpeciesDict[species]
                    newRmgSpecies = newSpeciesDict[species]
                    
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
                        if not curvesSimilar(timeOld.data, variableOld.data, timeNew.data, variableNew.data, tol, method):
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
                        if not failHeaderPrinted:
                            print failHeader
                            failHeaderPrinted=True
                        if i not in conditionsBroken: conditionsBroken.append(i)
                        print "Observable species {0} varied by more than {1:.3f} {2} between old model {3} and " \
                              "new model {4} in condition {5:d}.".format(smiles,
                                                                         tol,
                                                                         methodString[method],
                                                                         variableOld.label,
                                                                         variableNew.label,
                                                                         i+1)
                        variablesFailed.append((self.conditions[i], smiles, variableOld, variableNew))
                    
            
            # Compare state variable observables
            if 'variable' in self.observables:
                for varName in self.observables['variable']:
                    variableOld = next((data for data in dataListOld if data.label == varName), None)
                    variableNew = next((data for data in dataListNew if data.label == varName), None)
                    if not curvesSimilar(timeOld.data, variableOld.data, timeNew.data, variableNew.data, tol, method):
                        if i not in conditionsBroken: conditionsBroken.append(i)
                        if not failHeaderPrinted:
                            failHeaderPrinted=True
                            print failHeader

                        print "Observable variable {0} varied by more than {1:.3f} {2} between old model and \
new model in condition {3:d}.".format(variableOld.label, methodString[method], i+1)
                        variablesFailed.append((self.conditions[i], tol, varName, variableOld, variableNew))
                    
                    if plot:
                        oldVarPlot = GenericPlot(xVar=timeOld, yVar=variableOld)
                        newVarPlot = GenericPlot(xVar=timeNew, yVar=variableNew)
                        oldVarPlot.comparePlot(newSpeciesPlot,
                                                   title='Observable Variable {0} Comparison'.format(varName),
                                                   filename='condition_{0}_variable_{1}.png'.format(i+1, varName))
                        
            # Compare ignition delay observables
            if 'ignitionDelay' in self.observables:
                print 'Ignition delay observable comparison not implemented yet.'
                
                
        if failHeaderPrinted:
            print ''
            print 'The following reaction conditions were had some discrepancies:'
            print ''
            for index in conditionsBroken:
                print "Condition {0:d}:".format(index+1)
                print str(self.conditions[index])
                print ''

            return variablesFailed
        else:
            print ''
            if method=='avgDiff':
                print 'All Observables varied by less than {0:.3f} on average between old model and ' \
                      'new model in all conditions!'.format(tol)
            elif method=='maxDiff':
                print 'All Observables always varied by less than {0:.3f} between old model and ' \
                      'new model in all conditions!'.format(tol)
            print ''

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
