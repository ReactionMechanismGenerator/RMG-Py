import sys
import os.path
import numpy as np
import cantera as ct
from rmgpy.chemkin import loadSpeciesDictionary, getSpeciesIdentifier
from rmgpy.species import Species
from rmgpy.tools.plot import *

class CanteraCondition:
    """
    This class organizes the inputs needed for a cantera simulation

    ======================= ====================================================
    Attribute               Description
    ======================= ====================================================
    `reactorType`           A string of the cantera reactor type. List of supported types below:
        IdealGasReactor: A constant volume, zero-dimensional reactor for ideal gas mixtures
        IdealGasConstPressureReactor: A homogeneous, constant pressure, zero-dimensional reactor for ideal gas mixtures

    `reactionTime`          A float giving the reaction time in seconds
    `molFrac`               A dictionary giving the initial mol Fractions. Keys are species cantera names and the values are floats

    To specifiy the system for an ideal gas, you must define 2 of the following 3 parameters:
    `T0`                    A float giving the initial temperature in K
    'P0'                    A float giving the initial pressure in Pa
    'V0'                    A float giving the initial reactor volume in m^3
    ======================= ====================================================


    """
    def __init__(self, reactorType, reactionTime, molFrac, T0=None, P0=None, V0=None):
        self.reactorType=reactorType
        self.reactionTime=float(reactionTime)
        
        #Normalize initialMolFrac if not already done:
        if sum(molFrac.values())!=1.00:
            total=sum(molFrac.values())
            for species, value in molFrac.iteritems():
                molFrac[species]= value / total
        self.molFrac=molFrac
        self.T0=float(T0) if T0 else None
        self.P0=float(P0) if P0 else None
        self.V0=float(V0) if V0 else None

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string="Condition("
        string += 'reactorType="{0}", '.format(self.reactorType)
        string += 'reactionTime={:0.10f}, '.format(self.reactionTime)
        string += 'molFrac={0}, '.format(self.molFrac.__repr__())
        if self.T0: string += 'T0={:0.10f}, '.format(self.T0)
        if self.P0: string += 'P0={:0.10f}, '.format(self.P0)
        if self.V0: string += 'V0={:0.10f}, '.format(self.V0)
        string = string[:-2] + ')'
        return string

    def __str__(self):
        """
        Return a string representation of the condition.
        """
        string=""
        string += 'Reactor Type: {0}, '.format(self.reactorType)
        string += 'Reaction Time: {:0.10f}, '.format(self.reactionTime)
        if self.T0: string += 'T0: {:0.10f}, '.format(self.T0)
        if self.P0: string += 'P0: {:0.10f}, '.format(self.P0)
        if self.V0: string += 'V0: {:0.10f}, '.format(self.V0)
        string += 'Initial Mole Fractions: {0}, '.format(self.molFrac.__repr__())
        return string
    
class Cantera:
    """
    This class contains functions associated with an entire Cantera job
    """
    
    def __init__(self, speciesList=None, reactionList=None, canteraFile='', outputDirectory='', conditions=[]):
        """
        `speciesList`: list of RMG species objects
        `reactionList`: list of RMG reaction objects
        `canteraFile` path of the chem.cti file associated with this job
        `conditions`: a list of `CanteraCondition` objects
        """
        self.speciesList = speciesList 
        self.reactionList = reactionList 
        self.model = ct.Solution(canteraFile) if canteraFile else None
        self.outputDirectory = outputDirectory if outputDirectory else os.getcwd()
        self.conditions = conditions
        
    def loadChemkinModel(self, chemkinFile, **kwargs):
        """
        Convert a chemkin mechanism chem.inp file to a cantera mechanism file chem.cti 
        and save it in the outputDirectory
        Then load it into self.model
        """
        from cantera import ck2cti
        print 'Converting chem.inp to chem.cti...'
        outName=os.path.join(self.outputDirectory, "chem.cti")
        if os.path.exists(outName):
            os.remove(outName)
        parser = ck2cti.Parser()
        parser.convertMech(chemkinFile, outName=outName, **kwargs)
        print 'Saving into Cantera Model...'
        self.model = ct.Solution(outName)
    
    def generateConditions(self, reactorType, reactionTime, molFracList, Tlist=None, Plist=None, Vlist=None):
        """
        Creates a list of conditions from from the lists provided. 
        
        `reactorType`: a string indicating the Cantera reactor type
        `reactionTime`: ScalarQuantity object for time
        `molFracList`: a list of molfrac dictionaries wisth species keys and mole fraction values
        `Tlist`: ArrayQuantity object of temperatures
        `Plist`: ArrayQuantity object of pressures
        `Vlist`: ArrayQuantity object of volumes
        
        This saves all the reaction conditions into the Cantera class.
        """
        # First translate the molFracList from species objects to species names:
        newMolFracList = []
        for molFrac in molFracList:
            newMolFrac = {}
            for key, value in molFrac.iteritems():
                newkey = getSpeciesIdentifier(key)
                newMolFrac[newkey] = value
            newMolFracList.append(newMolFrac)
            
        molFracList = newMolFracList
        
        # Extract the si values from the arrays
        reactionTime = reactionTime.value_si
        Tlist = Tlist.value_si if Tlist is not None else None
        Plist = Plist.value_si if Plist is not None else None
        Vlist = Vlist.value_si if Vlist is not None else None
        
        
        conditions=[]
        
    
        if Tlist is None:
            for molFrac in molFracList:
                for P in Plist:
                    for V in Vlist:
                        conditions.append(CanteraCondition(reactorType, reactionTime, molFrac, P0=P, V0=V))
    
        elif Plist is None:
            for molFrac in molFracList:
                for T in Tlist:
                    for V in Vlist:
                        conditions.append(CanteraCondition(reactorType, reactionTime, molFrac, T0=T, V0=V))
    
        elif Vlist is None:
            for molFrac in molFracList:
                for T in Tlist:
                    for P in Plist:
                        conditions.append(CanteraCondition(reactorType, reactionTime, molFrac, T0=T, P0=P))
    
        else:
            for molFrac in molFracList:
                for T in Tlist:
                    for P in Plist:
                        for V in Vlist:
                            conditions.append(CanteraCondition(reactorType, reactionTime, molFrac, T0=T, P0=P, V0=V))
                            
        self.conditions = conditions

    def plot(self, data):
        """
        Plots data from the simulations from this cantera job.
        Takes data in the format of a list of tuples containing (time, [list of temperature, pressure, and species data]) 
        
        3 plots will be created for each condition:
        - T vs. time
        - P vs. time
        - Maximum 10 species mole fractions vs time
        """
        for i, conditionData in enumerate(data):
            time, dataList = conditionData
            # In RMG, any species with an index of -1 is an inert and should not be plotted
            inertList = [species for species in self.speciesList if species.index == -1 ]
            
            TData = dataList[0]
            PData = dataList[1]
            speciesData = [data for data in dataList if data.species not in inertList]
            
            # plot
            GenericPlot(xVar=time, yVar=TData).plot('{0}_temperature.png'.format(i+1))
            GenericPlot(xVar=time, yVar=PData).plot('{0}_pressure.png'.format(i+1))
            SimulationPlot(xVar=time, yVar=speciesData, ylabel='Mole Fraction').plot('{0}_mole_fractions.png'.format(i+1))
            
    def simulate(self):
        """
        Run all the conditions as a cantera simulation.
        Returns the data as a list of tuples containing: (time, [list of temperature, pressure, and species data]) 
            for each reactor condition
        """
        # Get all the cantera names for the species
        speciesNamesList = [getSpeciesIdentifier(species) for species in self.speciesList]
        
        allData = []
        for condition in self.conditions:
            
            # Set Cantera simulation conditions
            self.model.TPX = condition.T0, condition.P0, condition.molFrac
            
            # Choose reactor
            if condition.reactorType == 'IdealGasReactor':
                canteraReactor=ct.IdealGasReactor(self.model)
            else:
                raise Exception('Other types of reactor conditions are currently not supported')
            
            # Run this individual condition as a simulation
            canteraSimulation=ct.ReactorNet([canteraReactor])

            
            # Initialize the variables to be saved
            times = np.zeros(100)
            temperature = np.zeros(100, dtype=np.float64)
            pressure = np.zeros(100, dtype=np.float64)
            speciesData = np.zeros((100,len(speciesNamesList)),dtype=np.float64)

            
            # Begin integration
            time = 0.0
            # Run the simulation over 100 time points
            for n in range(100):
                time += condition.reactionTime/100
                
                # Advance the state of the reactor network in time from the current time to time t [s], taking as many integrator timesteps as necessary.
                canteraSimulation.advance(time)
                times[n] = time * 1e3  # time in ms
                temperature[n] = canteraReactor.T
                pressure[n] = canteraReactor.thermo.P
                speciesData[n,:] = canteraReactor.thermo[speciesNamesList].X
                
            # Resave data into generic data objects
            time = GenericData(label = 'Time', 
                               data = times,
                               units = 'ms')
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
            
            allData.append((time,conditionData))
            
        return allData