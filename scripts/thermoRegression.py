#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script compares the group additivity to species in a thermo library
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
def previousCommit(commitStr):
    """
    Context to go to a previous commit before reverting back to the the original commit
    """
    try:
        #Check that the commit str is not a file, so that you don't accidently check out unsaved changes
        if os.path.exists(commitStr):
            raise Exception("Argument commitStr must be a sha1, not a file")
        subprocess.call(['git', 'checkout', commitStr])
        yield
    finally:
        subprocess.call(['git', 'checkout', '-'])


################################################################################

if __name__ == '__main__':
    libraryName="burcat"
    # libraryName = "SulfurLibrary"
    #load database with no libraries
    groupsDatabase= ThermoDatabase()
    groupsDatabase.load(os.path.join(settings['database.directory'], 'thermo'))
    # commitStrs=['eb15e4']
    commitStrs=['HEAD~1']

    library = ThermoLibrary()
    library.load(os.path.join(settings['database.directory'], "thermo/libraries/" + libraryName + ".py"),
                 groupsDatabase.local_context,
                 groupsDatabase.global_context)


    parityData={}
    commitNames=['HEAD'] + commitStrs

    #get initial dictionary value: [data from library, data from groups of initial commit]
    for name, entry in library.entries.iteritems():
        test = Species().fromSMILES(entry.item.toSMILES())
        test.thermo = groupsDatabase.getThermoDataFromGroups(test)
        parityData[name] = [entry.data.getFreeEnergy(1000), test.getFreeEnergy(1000)]

    #change to new directory for RMG-database
    with cd(settings['database.directory']):
        for commitStr in commitStrs:
            with previousCommit(commitStr):
                newGroups= ThermoDatabase()
                newGroups.load(os.path.join(settings['database.directory'], 'thermo'))
                for name, entry in library.entries.iteritems():
                    test = Species().fromSMILES(entry.item.toSMILES())
                    test.thermo = newGroups.getThermoDataFromGroups(test)
                    parityData[name].append(test.getFreeEnergy(1000))

    # subprocess.call(["git", 'checkout', 'HEAD~1'])
    # subprocess.call(["git", 'checkout', 'HEAD'])
    # subprocess.call(["ls"])

    #make list of highst stdDev of groups (first element is name, second is library value)
    csvList = [[name]+data for name, data in parityData.iteritems()]
    #Put into kJ

    newCsvList=[]
    for point in csvList:
        newPoint=[point[0]]
        for value in point[1:]:
            newPoint.append(value/1000.0/4.182)
        newCsvList.append(newPoint)

    csvList= newCsvList

    csvList = sorted(csvList, key = lambda x: std(x[2:]))
    csvList.reverse()
    for item in csvList: print item, abs(item[2]-item[3])

    with open(os.path.join(os.getcwd(),'output.csv'), 'wb') as output:
        spamwriter=csv.writer(output)
        for line in csvList:
            spamwriter.writerow(line)


    #remove same points:
    prunedCsvList=[]
    for item in csvList:
        if not std(item[2:]) == 0:
            prunedCsvList.append(item)
    print "Entries in library:", len(csvList), "Entries with group value changed between commits", len(prunedCsvList)
    #This is if I don't want to prune:
    # prunedCsvList = csvList
    #make plot

    xvar = [x+1 for x in range(len(prunedCsvList))]
    for index, x in enumerate(range(len(prunedCsvList[0])-2)):
        yvar = [abs(x[1]-x[index+2]) for x in prunedCsvList]
        # plt.plot(xvar, yvar, 's')
        plt.semilogy(xvar, yvar, 's', label = commitNames[index])

    plt.xlabel("Species index")
    plt.ylabel("$\Delta\Delta$G$^f$$_{1000}$ (kcal/mol)")
    plt.legend()
    plt.savefig(os.path.join(os.getcwd(),'parity.png'))