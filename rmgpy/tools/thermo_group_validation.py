#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script compares the group additivity to species in a thermo library. It is also used to compare
differences in group additivity between different commits of RMG-database.
"""
from contextlib import contextmanager
import os
import csv
import subprocess
from numpy import std
import argparse
from rmgpy import settings
from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.species import Species
import matplotlib.pyplot as plt

#initalize input variables
global speciesThermo, inputDict
speciesThermo = {}
inputDict = {}


@contextmanager
def cd(newdir):
    """
    Context to change directory to newdir before reverting back to the original working directory
    """
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

@contextmanager
def previousCommit(commitStr, gitDirectory):
    """
    Context to go to a previous commit before reverting back to the the original commit
    """
    if commitStr == 'HEAD': yield
    else:
        try:
            #Check that the commit str is not a file, so that you don't accidently check out unsaved changes
            if os.path.exists(commitStr):
                raise Exception("Argument commitStr must be a sha1, not a file")
            else:
                subprocess.call(['git', 'checkout', commitStr])
                actualCommit = subprocess.check_output(['git', 'log','--format=%H%n%cd', '-1'], cwd=gitDirectory).splitlines()[0]
                assert commitStr in actualCommit, "Checkout to another commit failed, please make sure your changes are stashed or committed."
            yield
        finally:
            subprocess.call(['git', 'checkout', '-'])

class LibraryEntry:
    """
    Class to keep track of the where the original entry came from
    """
    def __init__(self, entryName, thermoData, libraryName):
        self.entryName = entryName
        self.thermoData = thermoData
        self.libraryName = libraryName

def groupValidationInputs(thermoLibraries= None,
                          commitStrs = None,
                          commitNames= None,
                          variable = 'H',
                          temperature = 298,):
    #Save the data from libraries into the speciesThermo dictionary
    if isinstance(thermoLibraries, str): libraryList = [thermoLibraries]
    groupsDatabase= ThermoDatabase()
    groupsDatabase.load(os.path.join(settings['database.directory'], 'thermo'))

    for libraryName in thermoLibraries:
        library = ThermoLibrary()
        library.load(os.path.join(settings['database.directory'], "thermo/libraries/" + libraryName + ".py"),
                     groupsDatabase.local_context,
                     groupsDatabase.global_context)
        for name, entry in library.entries.iteritems():
            if not entry.item.toSMILES() in speciesThermo:
                speciesThermo[entry.item.toSMILES()] = [LibraryEntry(entryName = name, thermoData = entry.data, libraryName = libraryName)]

    #Save other variables into inputDict
    if isinstance(commitStrs, str): commitStrs = [commitStrs]
    if isinstance(commitNames, str): commitNames = [commitNames]
    inputDict['commitStrs'] = commitStrs
    inputDict['commitNames'] = commitNames
    inputDict['variable'] = variable
    inputDict['temperature'] = float(temperature)


def loadInput(path):

    full_path = os.path.abspath(os.path.expandvars(path))
    f = open(full_path)

    global_context = { '__builtins__': None }
    local_context = {
        '__builtins__': None,
        'True': True,
        'False': False,
        'groupValidationInputs': groupValidationInputs,
    }

    try:
        exec f in global_context, local_context
    except (NameError, TypeError, SyntaxError), e:
        raise Exception(e)
    finally:
        f.close()

def parityFunction(variable, thermoData, T0):
    if variable == 'H':
        return thermoData.getEnthalpy(T0)
    elif variable =='S':
        return thermoData.getEntropy(T0)
    elif variable == 'G':
        return thermoData.getFreeEnergy(T0)
    elif variable == 'Cp':
        return thermoData.getHeatCapacity(T0)
    else:
        raise Exception("Variable {0} is not a supported thermo attribute for comparison".format(variable))

def execute(path):

    #key is RMG-generated SMILES, value is list with first entry as LibraryEntry
    #and the rest ThermoData from group estimation

    loadInput(path)
    variable = inputDict['variable']
    heading = ''
    units =''
    if variable == 'H':
        heading = 'Enthalpy'
        units = 'kcal/mol'
    elif variable == 'S':
        heading = 'Entropy'
        units = 'cal/mol'
    elif variable == 'Cp':
        heading = 'Heat_Capacity'
        units = 'cal/mol'
    elif variable == 'G':
        heading = 'Free_Energy'
        units = 'kcal/mol'
    T0 = inputDict['temperature']
    #load database with no libraries

    #change to new directory for RMG-database
    with cd(settings['database.directory']):
        for commitStr in inputDict['commitStrs']:
            with previousCommit(commitStr, settings['database.directory']):
                newGroups= ThermoDatabase()
                newGroups.load(os.path.join(settings['database.directory'], 'thermo'))
                for smiles in speciesThermo:
                    newSpecies = Species().fromSMILES(smiles)
                    speciesThermo[smiles].append(newGroups.getThermoDataFromGroups(newSpecies))

    parityData={}

    #get initial dictionary value: [data from library, data from groups of initial commit]
    for smiles, thermoList in speciesThermo.iteritems():
        parityData[smiles] = []
        for thermoData in thermoList:
            if isinstance(thermoData, LibraryEntry):
                parityData[smiles].append(parityFunction(variable, thermoData.thermoData, T0))
            else: parityData[smiles].append(parityFunction(variable, thermoData, T0))

    #make list of highst stdDev of groups (first element is name, second is library value)
    csvList = [[name]+data for name, data in parityData.iteritems()]

    #Put into kcal/cal
    conversionCsvList=[]
    for point in csvList:
        newPoint=[point[0]]
        for value in point[1:]:
            if variable == 'H' or variable == 'G':
                newPoint.append(value/1000.0/4.182)
            else:
                newPoint.append(value/4.182)
        conversionCsvList.append(newPoint)

    csvList= conversionCsvList

    #Sort based on largest difference in group estimation
    csvList = sorted(csvList, key = lambda x: std(x[2:]))
    csvList.reverse()

    with open(os.path.join(os.getcwd(), heading+'.csv'), 'wb') as output:
        spamwriter=csv.writer(output)
        #add Header
        spamwriter.writerow(['SMILES', 'Library Value at '+str(T0)+'K']+ inputDict['commitNames'])
        for line in csvList:
            spamwriter.writerow(line)

    #Remove points where group estimates are the same for plotting
    prunedCsvList=[]
    #don't include header in pruned csv
    for item in csvList:
        if not std(item[2:]) == 0:
            prunedCsvList.append(item)

    xvar = [x+1 for x in range(len(prunedCsvList))]
    for index, x in enumerate(range(len(prunedCsvList[0])-2)):
        yvar = [abs(x[1]-x[index+2]) for x in prunedCsvList]
        plt.semilogy(xvar, yvar, 's', label = inputDict['commitNames'][index])

    plt.xlabel("Species index")
    plt.ylabel("Abs Error from Library for $\Delta${0} ({1}".format(variable, units))
    plt.legend()
    plt.savefig(os.path.join(os.getcwd(),'comparison.png'))

    #Print some statistics:
    print ""
    print "Entries in library(s):", len(csvList)
    print "Entries with group value changed between commits:", len(prunedCsvList)
    print ""
    avgErrorAll = None #Average error across all library entries per commit
    avgErrorChanged = None #Average error for library entries where there was a difference in the group additivity
    for index, commitName in enumerate(inputDict['commitNames']):
        avgErrorAll = round(sum([abs(x[1]-x[index+2]) for x in csvList])/len(csvList), 2)
        avgErrorChanged = round(sum([abs(x[1]-x[index+2]) for x in prunedCsvList])/len(prunedCsvList), 2)
        print commitName, ": has an average error of", str(avgErrorAll), units, "for all entries, and", str(avgErrorChanged), units, "for entries with different group estimation."
        print ""

def main():
    """
    Driver function that parses command line arguments and passes them to the execute function.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('input', metavar='INPUT', type=str, nargs=1,
        help='Thermo Group Validation input file')
    args = parser.parse_args()

    path = os.path.abspath(args.input[0])

    execute(path)