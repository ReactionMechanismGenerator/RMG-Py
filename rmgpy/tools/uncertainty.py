import os
import numpy
import rmgpy.util as util 
from rmgpy.tools.plot import *

def retrieveSaturatedSpeciesFromList(species, speciesList):
    """
    Given a radical `species`, this function retrieves the saturated species objects from a list of species objects.
    """
    molecule = species.molecule[0]
    assert molecule.isRadical(), "Method only valid for radicals."
    saturatedStruct = molecule.copy(deep=True)
    saturatedStruct.saturate()
    for otherSpecies in speciesList:
        if otherSpecies.isIsomorphic(saturatedStruct):
            return otherSpecies
    else:
        raise Exception('Could not retrieve saturated species form of {0} from the species list'.format(species))



class Uncertainty:
    """
    This class contains functions associated with running uncertainty analyses
    for a single RMG-generated mechanism.
    """

    def __init__(self, speciesList=None, reactionList=None, outputDirectory=''):
        """
        `speciesList`: list of RMG species objects
        `reactionList`: list of RMG reaction objects
        `outputDirectoy`: directory path for saving output files from the analyses
        """
        self.database = None
        self.speciesList = speciesList 
        self.reactionList = reactionList 
        self.speciesSourcesDict = None
        self.reactionSourcesDict = None
        self.allThermoSources = None
        self.allKineticSources = None
        self.thermoInputUncertainties = None
        self.kineticInputUncertainties = None
        self.outputDirectory = outputDirectory if outputDirectory else os.getcwd()
        
        # Make output directory if it does not yet exist:
        if not os.path.exists(self.outputDirectory):
            try:
                os.makedirs(self.outputDirectory)
            except:
                raise Exception('Uncertainty output directory could not be created.')
            
    def loadDatabase(self, kineticsFamilies='all',kineticsDepositories=None,thermoLibraries=None, reactionLibraries=None):
        """
        This function loads a single copy of the RMGDatabase with full verbose averaging
        of the rate rule to trace kinetics sources.  
        
        By default, this function loads all the kinetics families, only the training kinetics depository,
        the primaryThermoLibrary, and no reaction libraries.  
        """
        from rmgpy.data.rmg import RMGDatabase 
        from rmgpy import settings
        
        if not kineticsDepositories:
            kineticsDepositories = ['training']
        if not thermoLibraries:
            thermoLibraries = ['primaryThermoLibrary']
        if not reactionLibraries:
            reactionLibraries = []
            
            
        self.database = RMGDatabase()
        self.database.load(settings['database.directory'], 
                          kineticsFamilies=kineticsFamilies, 
                          kineticsDepositories=kineticsDepositories,
                          thermoLibraries=thermoLibraries,
                          reactionLibraries=reactionLibraries,
                          )
        
    
        # Prepare the database by loading training reactions but not averaging the rate rules
        for familyLabel, family in self.database.kinetics.families.iteritems():
            family.addKineticsRulesFromTrainingSet(thermoDatabase=self.database.thermo)
        
            family.fillKineticsRulesByAveragingUp(verbose=True)

    def loadModel(self, chemkinPath, dictionaryPath, transportPath=None):
        """
        Load a RMG-generated model into the Uncertainty class
        `chemkinPath`: path to the chem_annotated.inp CHEMKIN mechanism 
        `dictionaryPath`: path to the species_dictionary.txt file 
        `transportPath`: path to the tran.dat file (optional)

        Then create dictionaries stored in self.thermoGroups and self.rateRules
        containing information about the source of the thermodynamic and kinetic
        parameters
        """
        from rmgpy.chemkin import loadChemkinFile

        self.speciesList, self.reactionList = loadChemkinFile(chemkinPath,
                                                              dictionaryPath=dictionaryPath,
                                                              transportPath=transportPath)


    def extractSourcesFromModel(self):
        """
        Extract the source data from the model using its comments.
        Must be done after loading model and database to work.
        """
        self.speciesSourcesDict = {}
        for species in self.speciesList:
            source = self.database.thermo.extractSourceFromComments(species)
            
            # Now prep the source data
            # Do not alter the GAV information, but reassign QM and Library sources to the species indices that they came from
            if len(source.keys()) == 1:
                # The thermo came from a single source, so we know it comes from a value describing the exact species
                if 'Library' in source:
                    source['Library'] = self.speciesList.index(species)   # Use just the species index in self.speciesList, for better shorter printouts when debugging
                if 'QM' in source:
                    source['QM'] = self.speciesList.index(species)
                    
            elif len(source.keys()) == 2:
                # The thermo has two sources, which indicates it's an HBI correction on top of a library or QM value.  We must retrieve the original
                # saturated molecule's thermo instead of using the radical species as the source of thermo
                saturatedSpecies = retrieveSaturatedSpeciesFromList(species,self.speciesList)
                if 'Library' in source:
                    source['Library'] = self.speciesList.index(saturatedSpecies)
                if 'QM' in source:
                    source['QM'] = self.speciesList.index(saturatedSpecies)
            else:
                raise Exception('Source of thermo should not use more than two sources out of QM, Library, or GAV.')
            
            self.speciesSourcesDict[species] = source
        
        self.reactionSourcesDict = {}
        for reaction in self.reactionList:
            source = self.database.kinetics.extractSourceFromComments(reaction)
            # Prep the source data 
            # Consider any library or PDep reaction to be an independent parameter for now and assign the source to the index of the
            # reaction within self.reactionList
            if 'Library' in source:
                source['Library'] = self.reactionList.index(reaction)
            elif 'PDep' in source:
                source['PDep'] = self.reactionList.index(reaction)
            elif 'Training' in source:
                # Do nothing here because training source already saves the entry from the training reaction
                pass
            elif 'Rate Rules' in source:
                # Do nothing
                pass
            else:
                raise Exception('Source of kinetics must be either Library, PDep, Training, or Rate Rules')
            self.reactionSourcesDict[reaction] = source
        
        # Account for all the thermo sources
        allThermoSources = {'GAV':{}, 'Library':set(), 'QM':set()}
        for source in self.speciesSourcesDict.values():
            if 'GAV' in source:
                for groupType in source['GAV'].keys():
                    groupEntries = [groupTuple[0] for groupTuple in source['GAV'][groupType]]
                    if not groupType in allThermoSources['GAV']:
                        allThermoSources['GAV'][groupType] = set(groupEntries)
                    else:
                        allThermoSources['GAV'][groupType].update(groupEntries)
            if 'Library' in source:
                allThermoSources['Library'].add(source['Library'])  
            if 'QM' in source:
                allThermoSources['QM'].add(source['QM'])   

                
        # Convert to lists
        self.allThermoSources = {}
        self.allThermoSources['Library'] = list(allThermoSources['Library'])
        self.allThermoSources['QM'] = list(allThermoSources['QM'])
        self.allThermoSources['GAV'] = {}
        for groupType in allThermoSources['GAV'].keys():
            self.allThermoSources['GAV'][groupType] = list(allThermoSources['GAV'][groupType])
                
        # Account for all the kinetics sources
        allKineticSources = {'Rate Rules':{}, 'Training':{}, 'Library':[], 'PDep':[]}
        for source in self.reactionSourcesDict.values():
            if 'Training' in source:
                familyLabel = source['Training'][0]
                trainingEntry = source['Training'][1]
                if not familyLabel in allKineticSources['Training']:
                    allKineticSources['Training'][familyLabel] = set([trainingEntry])
                else:
                    allKineticSources['Training'][familyLabel].add(trainingEntry)
            elif 'Library' in source:
                allKineticSources['Library'].append(source['Library'])
            elif 'PDep' in source:
                allKineticSources['PDep'].append(source['PDep'])
            elif 'Rate Rules' in source:
                familyLabel = source['Rate Rules'][0]
                sourceDict = source['Rate Rules'][1]
                rules = sourceDict['rules']
                training = sourceDict['training']
                if rules:
                    ruleEntries = [ruleTuple[0] for ruleTuple in rules]
                    if not familyLabel in allKineticSources['Rate Rules']:
                        allKineticSources['Rate Rules'][familyLabel] = set(ruleEntries)
                    else:
                        allKineticSources['Rate Rules'][familyLabel].update(ruleEntries)
                if training:
                    # Even though they are from training reactions, we consider the rate rules derived from the training
                    # reactions to be noncorrelated, due to the fact that some may be reversed.
                    trainingRules = [trainingTuple[0] for trainingTuple in training]  # Pick the rate rule entries
                    if not familyLabel in allKineticSources['Rate Rules']:
                        allKineticSources['Rate Rules'][familyLabel] = set(trainingRules)
                    else:
                        allKineticSources['Rate Rules'][familyLabel].update(trainingRules)
        
        self.allKineticSources = {}
        self.allKineticSources['Library'] = allKineticSources['Library']
        self.allKineticSources['PDep'] = allKineticSources['PDep']
        # Convert to lists
        self.allKineticSources['Rate Rules'] = {}
        for familyLabel in allKineticSources['Rate Rules'].keys():
            self.allKineticSources['Rate Rules'][familyLabel] = list(allKineticSources['Rate Rules'][familyLabel])
            
        self.allKineticSources['Training'] = {}
        for familyLabel in allKineticSources['Training'].keys():
            self.allKineticSources['Training'][familyLabel] = list(allKineticSources['Training'][familyLabel])
            
 