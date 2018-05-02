################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import os.path
import numpy as np
import cantera as ct
from rmgpy.chemkin import getSpeciesIdentifier
from rmgpy.tools.data import GenericData
from rmgpy.tools.plot import GenericPlot, SimulationPlot, ReactionSensitivityPlot
from rmgpy.quantity import Quantity


class CanteraCondition:
    """
    This class organizes the inputs needed for a cantera simulation

    ======================= ====================================================
    Attribute               Description
    ======================= ====================================================
    `reactorType`           A string of the cantera reactor type. List of supported types below:
        IdealGasReactor: A constant volume, zero-dimensional reactor for ideal gas mixtures
        IdealGasConstPressureReactor: A homogeneous, constant pressure, zero-dimensional reactor for ideal gas mixtures
        IdealGasConstPressureTemperatureReactor: A homogenous, constant pressure and constant temperature, zero-dimensional reactor 
                            for ideal gas mixtures (the same as RMG's SimpleReactor)

    `reactionTime`          A tuple object giving the (reaction time, units)
    `molFrac`               A dictionary giving the initial mol Fractions. Keys are species objects and the values are floats

    To specify the system for an ideal gas, you must define 2 of the following 3 parameters:
    `T0`                    A tuple giving the (initial temperature, units) which reconstructs a Quantity object
    'P0'                    A tuple giving the (initial pressure, units) which reconstructs a Quantity object
    'V0'                    A tuple giving the (initial specific volume, units) which reconstructs a Quantity object
    ======================= ====================================================


    """
    def __init__(self, reactorType, reactionTime, molFrac, T0=None, P0=None, V0=None):
        self.reactorType=reactorType
        self.reactionTime=Quantity(reactionTime)
        
        # Normalize initialMolFrac if not already done:
        if sum(molFrac.values())!=1.00:
            total=sum(molFrac.values())
            for species, value in molFrac.iteritems():
                molFrac[species]= value / total

        self.molFrac=molFrac
        
        # Check to see that one of the three attributes T0, P0, and V0 is less unspecified
        props=[T0,P0,V0]
        total=0
        for prop in props:
            if prop is None: total+=1

        if not total==1:
            raise Exception("Cantera conditions must leave one of T0, P0, and V0 state variables unspecified")


        self.T0=Quantity(T0) if T0 else None
        self.P0=Quantity(P0) if P0 else None
        self.V0=Quantity(V0) if V0 else None


    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string="CanteraCondition("
        string += 'reactorType="{0}", '.format(self.reactorType)
        string += 'reactionTime={}, '.format(self.reactionTime.__repr__())
        string += 'molFrac={0}, '.format(self.molFrac.__repr__())
        if self.T0: string += 'T0={}, '.format(self.T0.__repr__())
        if self.P0: string += 'P0={}, '.format(self.P0.__repr__())
        if self.V0: string += 'V0={}, '.format(self.V0__repr__())
        string = string[:-2] + ')'
        return string

    def __str__(self):
        """
        Return a string representation of the condition.
        """
        string=""
        string += 'Reactor Type: {0}\n'.format(self.reactorType)
        string += 'Reaction Time: {}\n'.format(self.reactionTime)
        if self.T0: string += 'T0: {}\n'.format(self.T0)
        if self.P0: string += 'P0: {}\n'.format(self.P0)
        if self.V0: string += 'V0: {}\n'.format(self.V0)
        #ConvertMolFrac to SMILES for keys for display
        prettyMolFrac={}
        for key, value in self.molFrac.iteritems():
            prettyMolFrac[key.molecule[0].toSMILES()]=value
        string += 'Initial Mole Fractions: {0}'.format(prettyMolFrac.__repr__())
        return string


def generateCanteraConditions(reactorTypeList, reactionTimeList, molFracList, Tlist=None, Plist=None, Vlist=None):
        """
        Creates a list of cantera conditions from from the arguments provided. 
        
        ======================= ====================================================
        Argument                Description
        ======================= ====================================================
        `reactorTypeList`        A list of strings of the cantera reactor type. List of supported types below:
            IdealGasReactor: A constant volume, zero-dimensional reactor for ideal gas mixtures
            IdealGasConstPressureReactor: A homogeneous, constant pressure, zero-dimensional reactor for ideal gas mixtures
            IdealGasConstPressureTemperatureReactor: A homogenous, constant pressure and constant temperature, zero-dimensional reactor 
                                for ideal gas mixtures (the same as RMG's SimpleReactor)

        `reactionTimeList`      A tuple object giving the ([list of reaction times], units)
        `molFracList`           A list of molfrac dictionaries with species object keys
                               and mole fraction values
        To specify the system for an ideal gas, you must define 2 of the following 3 parameters:
        `T0List`                A tuple giving the ([list of initial temperatures], units)
        'P0List'                A tuple giving the ([list of initial pressures], units)
        'V0List'                A tuple giving the ([list of initial specific volumes], units)
    
        
        This saves all the reaction conditions into the Cantera class.
        """
        
        # Create individual ScalarQuantity objects for Tlist, Plist, Vlist, and reactionTimeList
        if Tlist:
            Tlist = Quantity(Tlist) # Be able to create a Quantity object from it first
            Tlist = [(Tlist.value[i],Tlist.units) for i in range(len(Tlist.value))]
        if Plist:
            Plist = Quantity(Plist)
            Plist = [(Plist.value[i],Plist.units) for i in range(len(Plist.value))]
        if Vlist:
            Vlist = Quantity(Vlist)
            Vlist = [(Vlist.value[i],Vlist.units) for i in range(len(Vlist.value))]
        if reactionTimeList:
            reactionTimeList = Quantity(reactionTimeList)
            reactionTimeList = [(reactionTimeList.value[i],reactionTimeList.units) for i in range(len(reactionTimeList.value))]
        
        conditions=[]
        
    
        if Tlist is None:
            for reactorType in reactorTypeList:
                for reactionTime in reactionTimeList:
                    for molFrac in molFracList:
                        for P in Plist:
                            for V in Vlist:
                                conditions.append(CanteraCondition(reactorType, reactionTime, molFrac, P0=P, V0=V))
    
        elif Plist is None:
            for reactorType in reactorTypeList:
                for reactionTime in reactionTimeList:
                    for molFrac in molFracList:
                        for T in Tlist:
                            for V in Vlist:
                                conditions.append(CanteraCondition(reactorType, reactionTime, molFrac, T0=T, V0=V))
    
        elif Vlist is None:
            for reactorType in reactorTypeList:
                for reactionTime in reactionTimeList:
                    for molFrac in molFracList:
                        for T in Tlist:
                            for P in Plist:
                                conditions.append(CanteraCondition(reactorType, reactionTime, molFrac, T0=T, P0=P))
    
        else: raise Exception("Cantera conditions must leave one of T0, P0, and V0 state variables unspecified")
        return conditions

class Cantera:
    """
    This class contains functions associated with an entire Cantera job
    """
    
    def __init__(self, speciesList=None, reactionList=None, canteraFile='', outputDirectory='', conditions=None, sensitiveSpecies = None):
        """
        `speciesList`: list of RMG species objects
        `reactionList`: list of RMG reaction objects
        `reactionMap`: dict mapping the RMG reaction index within the `reactionList` to cantera model reaction(s) indices
        `canteraFile` path of the chem.cti file associated with this job
        `conditions`: a list of `CanteraCondition` objects
        `sensitiveSpecies`: a list of RMG species objects for conductng sensitivity analysis on
        """
        self.speciesList = speciesList 
        self.reactionList = reactionList 
        self.reactionMap = {}
        self.model = ct.Solution(canteraFile) if canteraFile else None
        self.outputDirectory = outputDirectory if outputDirectory else os.getcwd()
        self.conditions = conditions if conditions else []
        self.sensitiveSpecies = sensitiveSpecies if sensitiveSpecies else []

        # Make output directory if it does not yet exist:
        if not os.path.exists(self.outputDirectory):
            try:
                os.makedirs(self.outputDirectory)
            except:
                raise Exception('Cantera output directory could not be created.')

    def generateConditions(self, reactorTypeList, reactionTimeList, molFracList, Tlist=None, Plist=None, Vlist=None):
        """
        This saves all the reaction conditions into the Cantera class.
        ======================= ====================================================
        Argument                Description
        ======================= ====================================================
        `reactorTypeList`        A list of strings of the cantera reactor type. List of supported types below:
            IdealGasReactor: A constant volume, zero-dimensional reactor for ideal gas mixtures
            IdealGasConstPressureReactor: A homogeneous, constant pressure, zero-dimensional reactor for ideal gas mixtures
            IdealGasConstPressureTemperatureReactor: A homogenous, constant pressure and constant temperature, zero-dimensional reactor 
                                for ideal gas mixtures (the same as RMG's SimpleReactor)

        `reactionTimeList`      A tuple object giving the ([list of reaction times], units)
        `molFracList`           A list of molfrac dictionaries with species object keys
                               and mole fraction values
        To specify the system for an ideal gas, you must define 2 of the following 3 parameters:
        `T0List`                A tuple giving the ([list of initial temperatures], units)
        'P0List'                A tuple giving the ([list of initial pressures], units)
        'V0List'                A tuple giving the ([list of initial specific volumes], units)
        """

        self.conditions = generateCanteraConditions(reactorTypeList, reactionTimeList, molFracList, Tlist, Plist)

    def loadModel(self):
        """
        Load a cantera Solution model from the job's own speciesList and reactionList attributes
        """

        ctSpecies =[spec.toCantera(useChemkinIdentifier = True) for spec in self.speciesList]

        self.reactionMap = {}
        ctReactions = []
        for rxn in self.reactionList:
            index = len(ctReactions)

            convertedReactions = rxn.toCantera(self.speciesList, useChemkinIdentifier = True)

            if isinstance(convertedReactions, list):
                indices = range(index, index+len(convertedReactions))
                ctReactions.extend(convertedReactions)
            else:
                indices = [index]
                ctReactions.append(convertedReactions)

            self.reactionMap[self.reactionList.index(rxn)] = indices

        self.model = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                          species=ctSpecies, reactions=ctReactions)

    def refreshModel(self):
        """
        Modification to thermo requires that the cantera model be refreshed to 
        recalculate reverse rate coefficients and equilibrium constants... 
        As soon as cantera has its own Kinetics().modify_thermo function in place,
        this function may be deprecated.
        """
        ctReactions = self.model.reactions()
        ctSpecies = self.model.species()

        self.model = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                          species=ctSpecies, reactions=ctReactions)

    def loadChemkinModel(self, chemkinFile, transportFile=None, **kwargs):
        """
        Convert a chemkin mechanism chem.inp file to a cantera mechanism file chem.cti 
        and save it in the outputDirectory
        Then load it into self.model
        """
        from cantera import ck2cti
        
        base = os.path.basename(chemkinFile)
        baseName = os.path.splitext(base)[0]
        outName = os.path.join(self.outputDirectory, baseName + ".cti")
        if os.path.exists(outName):
            os.remove(outName)
        parser = ck2cti.Parser()
        parser.convertMech(chemkinFile, transportFile=transportFile, outName=outName, **kwargs)
        self.model = ct.Solution(outName)

    def modifyReactionKinetics(self, rmgReactionIndex, rmgReaction):
        """
        Modify the corresponding cantera reaction's kinetics to match 
        the reaction kinetics of an `rmgReaction`, using the `rmgReactionIndex` to
        map to the corresponding reaction in the cantera model. Note that
        this method only works if there is a reactionMap available (therefore only when the cantera model
        is generated directly from rmg objects and not from a chemkin file)
        """
        indices = self.reactionMap[rmgReactionIndex]
        modified_ctReactions = rmgReaction.toCantera(self.speciesList, useChemkinIdentifier = True)
        if not isinstance(modified_ctReactions, list):
            modified_ctReactions = [modified_ctReactions]

        for i in range(len(indices)):
            self.model.modify_reaction(indices[i], modified_ctReactions[i])

    def modifySpeciesThermo(self, rmgSpeciesIndex, rmgSpecies, useChemkinIdentifier = False):
        """
        Modify the corresponding cantera species thermo to match that of a
        `rmgSpecies` object, given the `rmgSpeciesIndex` which indicates the
        index at which this species appears in the `speciesList`
        """
        modified_ctSpecies = rmgSpecies.toCantera(useChemkinIdentifier = useChemkinIdentifier)
        ctSpecies = self.model.species(rmgSpeciesIndex)
        ctSpecies.thermo = modified_ctSpecies.thermo

    def plot(self, data, topSpecies=10, topSensitiveReactions=10):
        """
        Plots data from the simulations from this cantera job.
        Takes data in the format of a list of tuples containing (time, [list of temperature, pressure, and species data]) 
        
        3 plots will be created for each condition:
        - T vs. time
        - P vs. time
        - Maximum species mole fractions (the number of species plotted is based on the `topSpecies` argument)
        
        Reaction sensitivity plots will also be plotted automatically if there were sensitivities evaluated.
        The number of reactions to be plotted is defined by the `topSensitiveReactions` argument.
        
        """
        numCtReactions = len(self.model.reactions())
        for i, conditionData in enumerate(data):
            time, dataList, reactionSensitivityData = conditionData
            # In RMG, any species with an index of -1 is an inert and should not be plotted
            inertList = [species for species in self.speciesList if species.index == -1 ]
            
            TData = dataList[0]
            PData = dataList[1]
            speciesData = [data for data in dataList if data.species not in inertList]
            
            # plot
            GenericPlot(xVar=time, yVar=TData).plot(os.path.join(self.outputDirectory,'{0}_temperature.png'.format(i+1)))
            GenericPlot(xVar=time, yVar=PData).plot(os.path.join(self.outputDirectory,'{0}_pressure.png'.format(i+1)))
            SimulationPlot(xVar=time, yVar=speciesData, numSpecies=topSpecies, ylabel='Mole Fraction').plot(os.path.join(self.outputDirectory,'{0}_mole_fractions.png'.format(i+1)))
            
            for j, species in enumerate(self.sensitiveSpecies):
                ReactionSensitivityPlot(xVar=time, yVar=reactionSensitivityData[j*numCtReactions:(j+1)*numCtReactions], numReactions=topSensitiveReactions).barplot(os.path.join(self.outputDirectory,'{0}_{1}_sensitivity.png'.format(i+1,species.toChemkin())))
            
    def simulate(self):
        """
        Run all the conditions as a cantera simulation.
        Returns the data as a list of tuples containing: (time, [list of temperature, pressure, and species data]) 
            for each reactor condition
        """
        # Get all the cantera names for the species
        speciesNamesList = [getSpeciesIdentifier(species) for species in self.speciesList]
        inertIndexList = [self.speciesList.index(species) for species in self.speciesList if species.index == -1]
        
        allData = []
        for condition in self.conditions:

            # First translate the molFrac from species objects to species names
            newMolFrac = {}
            for key, value in condition.molFrac.iteritems():
                newkey = getSpeciesIdentifier(key)
                newMolFrac[newkey] = value

            # Set Cantera simulation conditions
            if condition.V0 is None:
                self.model.TPX = condition.T0.value_si, condition.P0.value_si, newMolFrac
            elif condition.P0 is None:
                self.model.TDX = condition.T0.value_si, 1.0/condition.V0.value_si, newMolFrac
            else:
                raise Exception("Cantera conditions in which T0 and P0 or T0 and V0 are not the specified state variables are not yet implemented.")


            # Choose reactor
            if condition.reactorType == 'IdealGasReactor':
                canteraReactor=ct.IdealGasReactor(self.model)
            elif condition.reactorType == 'IdealGasConstPressureReactor':
                canteraReactor=ct.IdealGasConstPressureReactor(contents=self.model)
            elif condition.reactorType == 'IdealGasConstPressureTemperatureReactor':
                canteraReactor=ct.IdealGasConstPressureReactor(contents=self.model, energy='off')
            else:
                raise Exception('Other types of reactor conditions are currently not supported')
            
            # Run this individual condition as a simulation
            canteraSimulation=ct.ReactorNet([canteraReactor])
            
            numCtReactions = len(self.model.reactions())
            if self.sensitiveSpecies:
                if ct.__version__ == '2.2.1':
                    print 'Warning: Cantera version 2.2.1 may not support sensitivity analysis unless SUNDIALS was used during compilation.'
                    print 'Warning: Upgrade to newer of Cantera in anaconda using the command "conda update -c rmg cantera"'
                # Add all the reactions as part of the analysis
                for i in range(numCtReactions):
                    canteraReactor.add_sensitivity_reaction(i)
                # Set the tolerances for the sensitivity coefficients
                canteraSimulation.rtol_sensitivity = 1e-4
                canteraSimulation.atol_sensitivity = 1e-6
                
            # Initialize the variables to be saved
            times=[]
            temperature=[]
            pressure=[]
            speciesData=[]
            sensitivityData = []
            
            # Begin integration
            time = 0.0
            # Run the simulation over 100 time points
            while canteraSimulation.time<condition.reactionTime.value_si:

                # Advance the state of the reactor network in time from the current time to time t [s], taking as many integrator timesteps as necessary.
                canteraSimulation.step(condition.reactionTime.value_si)
                times.append(canteraSimulation.time)
                temperature.append(canteraReactor.T)
                pressure.append(canteraReactor.thermo.P)
                speciesData.append(canteraReactor.thermo[speciesNamesList].X)
                
                
                if self.sensitiveSpecies:
                    # Cantera returns mass-based sensitivities rather than molar concentration or mole fraction based sensitivities.
                    # The equation for converting between them is:
                    # 
                    # d ln xi = d ln wi - sum_(species i) (dln wi) (xi)
                    # 
                    # where xi is the mole fraction of species i and wi is the mass fraction of species i
                    
                    massFracSensitivityArray = canteraSimulation.sensitivities()
                    if condition.reactorType =='IdealGasReactor':
                        # Row 0: mass, Row 1: volume, Row 2: internal energy or temperature, Row 3+: mass fractions of species
                        massFracSensitivityArray = massFracSensitivityArray[3:,:]
                    elif condition.reactorType == 'IdealGasConstPressureReactor' or condition.reactorType == 'IdealGasConstPressureTemperatureReactor':
                        # Row 0: mass, Row 1: enthalpy or temperature, Row 2+: mass fractions of the species
                        massFracSensitivityArray = massFracSensitivityArray[2:,:]
                    else:
                        raise Exception('Other types of reactor conditions are currently not supported')
                    
                    for i in range(len(massFracSensitivityArray)):
                        massFracSensitivityArray[i] *= speciesData[-1][i]
                        
                    sensitivityArray= np.zeros(len(self.sensitiveSpecies)*len(self.model.reactions()))
                    for index, species in enumerate(self.sensitiveSpecies):
                        for j in range(numCtReactions):
                            sensitivityArray[numCtReactions*index+j] = canteraSimulation.sensitivity(species.toChemkin(),j)

                            for i in range(len(massFracSensitivityArray)):
                                if i not in inertIndexList:
                                    # massFracSensitivity for inerts are returned as nan in Cantera, so we must not include them here
                                    sensitivityArray[numCtReactions*index+j] -= massFracSensitivityArray[i][j]
                    sensitivityData.append(sensitivityArray)
                
            # Convert speciesData and sensitivityData to a numpy array
            speciesData=np.array(speciesData)
            sensitivityData = np.array(sensitivityData)

            # Resave data into generic data objects
            time = GenericData(label = 'Time', 
                               data = times,
                               units = 's')
            temperature = GenericData(label='Temperature',
                                      data = temperature,
                                      units = 'K')
            pressure = GenericData(label='Pressure',
                                      data = pressure,
                                      units = 'Pa')
            conditionData = []
            conditionData.append(temperature)
            conditionData.append(pressure)
            
            for index, species in enumerate(self.speciesList):
                # Create generic data object that saves the species object into the species object.  To allow easier manipulate later.
                speciesGenericData = GenericData(label=speciesNamesList[index],
                                          species = species,
                                          data = speciesData[:,index],
                                          index = species.index
                                          )
                conditionData.append(speciesGenericData)
            
            reactionSensitivityData = []
            for index, species in enumerate(self.sensitiveSpecies):
                for j in range(numCtReactions):
                    reactionSensitivityGenericData = GenericData(label = 'dln[{0}]/dln[k{1}]: {2}'.format(species.toChemkin(),j+1, self.model.reactions()[j]),
                                  species = species,
                                  reaction = self.model.reactions()[j],
                                  data = sensitivityData[:,numCtReactions*index+j],
                                  index = j+1,
                                  )
                    reactionSensitivityData.append(reactionSensitivityGenericData)
            
            allData.append((time,conditionData,reactionSensitivityData))
            
        return allData


def getRMGSpeciesFromUserSpecies(userList, RMGList):
    """
    Args:
        userList: list of generic Class Species Objects created by user
        speciesList: a list of RMG species objects

    This function takes a list of generic species objects and returns the species object generated from a loaded RMG
    dictionary, thereby gaining the correct label for a given mechanism.

    Returns: A dict containing the Species Object from userList and RMG Species objects as their values
    If the species is not found, the value will be returned as None
    """
    mapping = {}
    for userSpecies in userList:
        userSpecies.generate_resonance_structures()

        for rmgSpecies in RMGList:
            if userSpecies.isIsomorphic(rmgSpecies):
                if userSpecies in mapping:
                    raise KeyError("The Species with SMIlES {0} has appeared twice in the species list!".format(userSpecies.molecule[0].toSMILES()))
                mapping[userSpecies] = rmgSpecies
                break
        else: 
            mapping[userSpecies] = None

    return mapping

def findIgnitionDelay(time, yVar=None, metric='maxDerivative'):
    """
    Identify the ignition delay point based on the following parameters:

    `time`: an array containing different times
    `yVars`: either a single y array or a list of arrays, typically containing only a single array such as pressure,
             but can contain multiple arrays such as species
    `metric`: can be set to 
        'maxDerivative': This is selected by default for y(t_ign) = max(dY/dt), and is typically used for a yVar containing T or P data
        'maxHalfConcentration': This is selected for the case where a metric like [OH](t_ign) = [OH]_max/2 is desired
        'maxSpeciesConcentrations': This is selected for the case where the metric for ignition
            is y1*y2*...*yn(t_ign) = max(y1*y2*...*yn) such as when the time desired if for max([CH][O]).  This is
            the only metric that requires a list of arrays

    Note that numpy array must be used.
    """

    if not isinstance(yVar, list):
        yVar = [yVar]

    for y in yVar:
        if len(y) != len(time):
            raise Exception('Mismatch of array length for time and y variable.')

    if metric == 'maxDerivative':
        if len(yVar) != 1:
            raise Exception('Maximum derivative metric for ignition delay must be used with a single y variable.')

        y = yVar[0]
        dydt = (y[1:] - y[:-1]) / (time[1:] - time[:-1])
        index = next(i for i,d in enumerate(dydt) if d==max(dydt))
        
        return 0.5 * (time[index] + time[index+1])
    elif metric == 'maxHalfConcentration':
        if len(yVar) != 1:
            raise Exception('Max([OH]/2) metric for ignition delay must be used with a single y variable.')

        y = yVar[0]
        maxIndex = y.argmax()
        OHmetric = max(y)/2 
        mindata = OHmetric - y[0:maxIndex]
        index = mindata.argmin()
        return time[index]

    elif metric == 'maxSpeciesConcentrations':
        multdata = np.ones(len(yVar[0]))
        for spec in yVar:  
            multdata *= spec
        index = multdata.argmax()
        return time[index]


def checkNearlyEqual(value1, value2, dE = 1e-5):
    """
    Check that two values are nearly equivalent by abs(val1-val2) < abs(dE*val1)
    """
    
    if abs(value1-value2) <= abs(dE*value1) or abs(value1-value2) <= abs(dE*value2) or abs(value1-value2) <= dE:
        return True
    else:
        return False
    
    
def checkEquivalentCanteraSpecies(ctSpec1, ctSpec2, dE=1e-5):
    """
    Checks that the two cantera species are nearly equivalent
    """
    try:
        assert ctSpec1.name == ctSpec2.name, "Identical name"
        assert ctSpec1.composition == ctSpec2.composition, "Identical composition"
        assert ctSpec1.size == ctSpec2.size, "Identical species size"
        assert ctSpec1.charge == ctSpec2.charge, "Identical charge"

        if ctSpec1.transport or ctSpec2.transport:
            trans1 = ctSpec1.transport
            trans2 = ctSpec2.transport

            assert checkNearlyEqual(trans1.acentric_factor, trans2.acentric_factor, dE), "Identical acentric factor"
            assert checkNearlyEqual(trans1.diameter, trans2.diameter, dE), "Identical diameter"
            assert checkNearlyEqual(trans1.dipole, trans2.dipole, dE), "Identical dipole moment"
            assert trans1.geometry == trans2.geometry, "Identical geometry"
            assert checkNearlyEqual(trans1.polarizability, trans2.polarizability, dE), "Identical polarizibility"
            assert checkNearlyEqual(trans1.rotational_relaxation, trans2.rotational_relaxation, dE), "Identical rotational relaxation number"
            assert checkNearlyEqual(trans1.well_depth, trans2.well_depth, dE), "Identical well depth"

        if ctSpec1.thermo or ctSpec2.thermo:
            thermo1 = ctSpec1.thermo
            thermo2 = ctSpec2.thermo

            Tlist = [300,500,1000,1500,2000]
            for T in Tlist:
                assert checkNearlyEqual(thermo1.cp(T), thermo2.cp(T), dE),  "Similar heat capacity"
                assert checkNearlyEqual(thermo1.h(T), thermo2.h(T), dE), "Similar enthalpy"
                assert checkNearlyEqual(thermo1.s(T), thermo2.s(T), dE),  "Similar entropy"
    except Exception as e:
        print "Cantera species {0} failed equivalency check on: {1}".format(ctSpec1,e)
        return False

    return True
    
def checkEquivalentCanteraReaction(ctRxn1, ctRxn2, checkID=False, dE=1e-5):
    """
    Checks that the two cantera species are nearly equivalent
    if checkID is True, then ID's for the reactions will also be checked
    """
    def checkEquivalentArrhenius(arr1, arr2):
        assert checkNearlyEqual(arr1.activation_energy, arr2.activation_energy, dE), "Similar Arrhenius Ea"
        assert checkNearlyEqual(arr1.pre_exponential_factor, arr2.pre_exponential_factor, dE), "Similar Arrhenius A-factor"
        assert checkNearlyEqual(arr1.temperature_exponent, arr2.temperature_exponent, dE), "Similar Arrhenius temperature exponent"
    
    def checkEquivalentFalloff(fall1, fall2):
        assert len(fall1.parameters) == len(fall2.parameters), "Same number of falloff parameters"
        for i in range(len(fall1.parameters)):
            assert checkNearlyEqual(fall1.parameters[i], fall2.parameters[i], dE), "Similar falloff parameters"
        assert fall1.type == fall2.type, "Same falloff parameterization type"
    
    try:
        assert type(ctRxn1) == type(ctRxn2), "Same Cantera reaction type"

        if isinstance(ctRxn1, list):
            assert len(ctRxn1) == len(ctRxn2), "Same number of reactions"
            for i in range(len(ctRxn1)):
                checkEquivalentCanteraReaction(ctRxn1[i], ctRxn2[i], checkID=checkID)



        if checkID:
            assert ctRxn1.ID == ctRxn2.ID, "Same reaction ID"

        assert ctRxn1.duplicate == ctRxn2.duplicate, "Same duplicate attribute"
        assert ctRxn1.reversible == ctRxn2.reversible, "Same reversible attribute"
        assert ctRxn1.orders == ctRxn2.orders, "Same orders attribute"
        assert ctRxn1.allow_negative_orders == ctRxn2.allow_negative_orders, "Same allow_negative_orders attribute"
        assert ctRxn1.allow_nonreactant_orders == ctRxn2.allow_nonreactant_orders, "Same allow_nonreactant_orders attribute"
        assert ctRxn1.reactants == ctRxn2.reactants, "Same reactants"
        assert ctRxn1.products == ctRxn2.products, "Same products"


        if isinstance(ctRxn1, ct.ElementaryReaction):
            assert ctRxn1.allow_negative_pre_exponential_factor == ctRxn2.allow_negative_pre_exponential_factor, \
                "Same allow_negative_pre_exponential_factor attribute"
            if ctRxn1.rate or ctRxn2.rate:
                checkEquivalentArrhenius(ctRxn1.rate,ctRxn2.rate)

        elif isinstance(ctRxn1, ct.PlogReaction):
            if ctRxn1.rates or ctRxn2.rates:
                assert len(ctRxn1.rates) == len(ctRxn2.rates), "Same number of rates in PLOG reaction"

                for i in range(len(ctRxn1.rates)):
                    P1, arr1 = ctRxn1.rates[i]
                    P2, arr2 = ctRxn2.rates[i]
                    assert checkNearlyEqual(P1, P2, dE), "Similar pressures for PLOG rates"
                    checkEquivalentArrhenius(arr1, arr2)

        elif isinstance(ctRxn1, ct.ChebyshevReaction):
            assert ctRxn1.Pmax == ctRxn2.Pmax, "Same Pmax for Chebyshev reaction" 
            assert ctRxn1.Pmin == ctRxn2.Pmin, "Same Pmin for Chebyshev reaction" 
            assert ctRxn1.Tmax == ctRxn2.Tmax, "Same Tmax for Chebyshev reaction" 
            assert ctRxn1.Tmin == ctRxn2.Tmin, "Same Tmin for Chebyshev reaction" 
            assert ctRxn1.nPressure == ctRxn2.nPressure, "Same number of pressure interpolations"
            assert ctRxn1.nTemperature == ctRxn2.nTemperature, "Same number of temperature interpolations"
            for i in range(ctRxn1.coeffs.shape[0]):
                for j in range(ctRxn1.coeffs.shape[1]):
                    assert checkNearlyEqual(ctRxn1.coeffs[i,j], ctRxn2.coeffs[i,j], dE), \
                    "Similar Chebyshev coefficients" 

        elif isinstance(ctRxn1, ct.ThreeBodyReaction):
            assert ctRxn1.default_efficiency == ctRxn2.default_efficiency, "Same default efficiency" 
            assert ctRxn1.efficiencies == ctRxn2.efficiencies, "Same efficienciess" 

        elif isinstance(ctRxn1, ct.FalloffReaction):
            assert ctRxn1.default_efficiency == ctRxn2.default_efficiency, "Same default efficiency" 
            assert ctRxn1.efficiencies == ctRxn2.efficiencies, "Same efficienciess" 
            if ctRxn1.falloff or ctRxn2.falloff:
                checkEquivalentFalloff(ctRxn1.falloff,ctRxn2.falloff)
            if ctRxn1.high_rate or ctRxn2.high_rate:
                checkEquivalentArrhenius(ctRxn1.high_rate, ctRxn2.high_rate)
            if ctRxn1.low_rate or ctRxn2.low_rate:
                checkEquivalentArrhenius(ctRxn1.low_rate, ctRxn2.low_rate)
                
    except Exception as e:
        print "Cantera reaction {0} failed equivalency check on: {1}".format(ctRxn1, e)
        return False
        
    return True
