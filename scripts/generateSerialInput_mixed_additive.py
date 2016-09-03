#!/usr/bin/env python
# encoding: utf-8

"""
This script generates an input.py using species concentrations from a previous RMG job.
It copies the previous input file, but adds the old mechanism as an appended library.

Then it takes the output from the LAST simulation done by RMG and uses this as a base for the new mole fractions to
put into the next RMG job.

It will then also add the additives denoted by additive_i_SMILES with mole fractions (before mixing) of additive_i_mole_frac 
to the initial mole fraction of the new RMG job

The ratio of the additive to reactor 1 stream mass flow rates is stream2Mass.

    $ python generateSerialInput.py name /path/to/oldWorkDir /path/to/newWorkDir stream2Mass --additiveSMILES 'additive_1_SMILES' 'additive_2_SMILES' ... 'additive_N_SMILES' 
	--additivemolefrac additive_1_mole_frac additive_2_mole_frac ... additive_N_mole_frac --Reactor_2_P_dep on or off --Temperatures Temperature1 Temperature2 ... TemperatureN
	--Reactor_2_Aromatic_libraries on or off --Reactor_2_seed_with_aromatics on or off

"""
from rmgpy.tools.plot import parseCSVData
from rmgpy.chemkin import getSpeciesIdentifier
from rmgpy.chemkin import loadSpeciesDictionary
from rmgpy.species import Species
import os.path
import argparse
import re
import copy as cp

"""
Still need to set up submission script that runs first RMG job, this script, importChemkinLibrary script, and then
second RMG job in sequence.

"""
################################################################################
if __name__ == '__main__':
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('name', type=str, help='the name of the burner library')
    parser.add_argument('oldWorkDir', type=str, help='Path to the old working directory')
    parser.add_argument('newWorkDir', type=str, help='Path to the new working directory')
    parser.add_argument('stream2Mass', type=float, help='Mass in second stream normalized over mass of first stream')
    #eg, if they have the same mass flow rate, use 1 for last argument, if second stream is twice as massive, use 2.

    parser.add_argument('--additiveSMILES', nargs='+', help='List of SMILES for the additives being put into the second reactor')
    parser.add_argument('--additivemolefrac', nargs='+',type=float, help='List of mole fractions (before mixing) of the additives being put into the second reactor')

    parser.add_argument('--Reactor_2_P_dep', type=str, help='Specify whether P_dep in Reactor 2 should be on or off', default='on')

    parser.add_argument('--Temperatures', nargs='+',type=float, help='List of temperatures to simulate')

    parser.add_argument('--Reactor_2_Aromatic_libraries', type=str, help='Specify whether to use Aromatic Libraries in Reactor 2 (Cracker)', default='off')

    parser.add_argument('--Reactor_2_seed_with_aromatics', type=str, help='Specify whether to seed aromatics with 0.0 initial mole fraction in Reactor 2 (Cracker)', default='off')
    
    parser.add_argument('--Reactor_2_seed_with_bicpd', type=str, help='Specify whether to seed bicpd aromatics with 0.0 initial mole fraction in Reactor 2 (Cracker)', default='off')

    parser.add_argument('--Reactor_2_tolmovetocore', type=float, help='Option to specify tolmovetocore', default='0.1')
    
    args = parser.parse_args()

    #Base of reactor block

    newReactorBlockTemp=[]

    for Temperature in args.Temperatures:
    	newReactorBlockTemp.append("simpleReactor(\n")
	newReactorBlockTemp.append("    temperature= (" + str(Temperature) + ",'K'),\n")
	newReactorBlockTemp.append("    pressure=(1.5,'atm'),\n")
	newReactorBlockTemp.append("    initialMoleFractions={\n")
	newReactorBlockTemp.append("    },\n")
	newReactorBlockTemp.append("    terminationTime=(.005,'s'),\n")
	newReactorBlockTemp.append(")\n\n")


###########Find correct simulation file#################################################################################
    #First find last reactor solver:
    lastSimulation=0
    solverPath=os.path.join(args.oldWorkDir, 'solver')
    for filename in os.listdir(solverPath):
        if filename.endswith(".csv"):
            number=re.sub('simulation\_', '', filename)
            number =int(re.sub('\_.*', '', number))
            if number > lastSimulation:
                lastSimulation=cp.copy(number)


    #Now find the correct file, it should be the solver with the highest number at the end
    lastIteration=0
    lastFile=''
    for filename in os.listdir(solverPath):
        if filename.endswith(".csv"):
            if re.search('simulation\_'+str(lastSimulation)+'\_', filename):
                number = re.sub('.*\_', '', filename)
                number =int(re.sub('\.csv', '', number))
                if number > lastIteration:
                    lastIteration=cp.copy(number)
                    lastFile=cp.copy(filename)

    #Extract data from file
    time, dataList=parseCSVData(os.path.join(solverPath, lastFile))




########################Write out correct reaction stream###############################################################
    #Import species dictionary
    speciesDict=loadSpeciesDictionary(os.path.join(args.oldWorkDir, "chemkin/species_dictionary.txt"))
    bathGases={"Ar": Species().fromSMILES("[Ar]"),
               "He": Species().fromSMILES("[He]"),
               "Ne": Species().fromSMILES("[Ne]"),
               "N2": Species().fromSMILES("N#N")}
    speciesDict.update(bathGases)

    #mol fraction of stream1
    stream1={} #key is a Species Object, value is the final mol fraction from the simulation
    tol = 1e-6
    for item in dataList:
        #Check that it is a species and not volume or a sensitivity
        if item.species is not None or item.label in bathGases:
            if item.data[-1]>tol:
	        if item.label not in speciesDict:
                    continue
		elif speciesDict[item.label].molecule[0].toSMILES() == "O=O":
                    continue
                else:
                    stream1[speciesDict[item.label]]=item.data[-1]

    #mol fraction of stream2 (additive stream)
    stream2={} #key is a Species Object, value is the mol fraction of additive stream before mixing
    additivePresent={}
    counter=0
    for SMILES in args.additiveSMILES:
	additiveSpecies = Species().fromSMILES(SMILES)
        stream2[additiveSpecies]=args.additivemolefrac[counter]
	additivePresent[additiveSpecies] = False
	counter+=1

    counter2=0
    if re.search('on', args.Reactor_2_seed_with_aromatics):
	seeded_aromatic_SMILES = ['C1=CC=CC=C1', 'C1C=CCC=1', 'C1=CC=C2CC=CC2=C1', 'C1=CC=C2C=CC=CC2=C1','C1C=C[C](C=1)C1C=CC=C1','C1=CCC(=C1)[C]1C=CC=C1','C1C=C[C](C=1)C1C=CCC=1']
	seeded_aromatic_labels = ['A1', 'CPD', 'INDENE', 'A2','2a','2b','2c']
	for SMILES in seeded_aromatic_SMILES:
	    additiveSpecies = Species().fromSMILES(SMILES)
	    additiveSpecies.label = seeded_aromatic_labels[counter2]
            stream2[additiveSpecies]=0.0
	    additivePresent[additiveSpecies] = False
	    counter+=1
	    counter2+=1
    
    counter2=0
    if re.search('on', args.Reactor_2_seed_with_bicpd):
	seeded_aromatic_bicpd_SMILES = ['C1=CCC(=C1)C1=CC=CC1', 'C1=CCC(=C1)C1C=CCC=1', 'C1=CC(=CC1)C1C=CCC=1', 'C1C=CC(C=1)C1C=CC=C1', 'C1=CCC(=C1)C1C=CC=C1', 'C1C=CC(C=1)C1C=CCC=1']
	seeded_aromatic_bicpd_labels = ['1c', '1d', '1e', '1a', '1b', '1f']
	for SMILES in seeded_aromatic_bicpd_SMILES:
	    additiveSpecies = Species().fromSMILES(SMILES)
            additiveSpecies.label = seeded_aromatic_bicpd_labels[counter2]
	    stream2[additiveSpecies]=0.0
	    additivePresent[additiveSpecies] = False
	    counter+=1
	    counter2+=1

    newStream={} #For the final stream that goes into input of second RMG JOB

    #Calculate the effective molecular weight of stream 1
    eff_MolecularWeight_1=0
    for component in stream1:
        eff_MolecularWeight_1+=stream1[component]*component.molecule[0].getMolecularWeight()

    #Calculate the effective molecular weight of stream 2
    eff_MolecularWeight_2=0
    for component in stream2:
        eff_MolecularWeight_2+=stream2[component]*component.molecule[0].getMolecularWeight()

    #convert everything to number of mols in stream1 and stream 2
    stream1Mols={}
    for component in stream1:
        stream1Mols[component]=stream1[component]/eff_MolecularWeight_1

    stream2Mols={}
    for component in stream2:
        stream2Mols[component]=stream2[component]*args.stream2Mass/eff_MolecularWeight_2

    #get total number of mols for new stream:
    total=0
    for component, value in stream1Mols.iteritems():
        total+=value
    for component, value in stream2Mols.iteritems():
        total+=value

    #update newStream with additive
    additivematch = False
    for species in stream1:
	for additiveSpecies in stream2:
	    if species.isIsomorphic(additiveSpecies):
                newStream[species]=(stream1Mols[species]+stream2Mols[additiveSpecies])/total
                additivePresent[additiveSpecies]=True
		additivematch = True
		break
	    else: additivematch = False
        if not additivematch: 
	    newStream[species]=stream1Mols[species]/total

    additiveName={}
    for additiveSpecies in stream2:
    	if not additivePresent[additiveSpecies]:
            newStream[additiveSpecies]=stream2Mols[additiveSpecies]/total
	    if not additiveSpecies.label:
		additiveName[additiveSpecies]=additiveSpecies.molecule[0].toSMILES() +"add" #getSpeciesIdentifier(additiveSpecies)+"add"
            else:
		additiveName[additiveSpecies]=additiveSpecies.label+"add"
	else: additiveName[additiveSpecies]=getSpeciesIdentifier(additiveSpecies)

    #Add species into reactor block
    additivematch = False
    newReactorBlock=[]
    for line in newReactorBlockTemp:
        newReactorBlock.append(line)
        if re.search("initialMoleFractions",line):
            for species, value in newStream.iteritems():
		for additiveSpecies in stream2:
		    if species==additiveSpecies and not additivePresent[additiveSpecies]:
                        newReactorBlock.append("        '"+ additiveName[additiveSpecies] + "': "+ str(value)+ ",\n")
			additivematch = True
			break
		    else: additivematch = False
                if not additivematch:
                    newReactorBlock.append("        '"+ getSpeciesIdentifier(species) + "': "+ str(value)+ ",\n")

#Construct adjList dictionary with correct indentation
    adjListDict={}
    for component in newStream:
        newAdjListList=re.split("\n", component.molecule[0].toAdjacencyList())
        newAdjList=""
        for line in newAdjListList:
            newAdjList+="        "+line+"\n"
        adjListDict[component]=newAdjList

    for bathGas in bathGases:
        newAdjListList=re.split("\n", bathGases[bathGas].molecule[0].toAdjacencyList())
        newAdjList=""
        for line in newAdjListList:
            newAdjList+="        "+line+"\n"
        adjListDict[bathGases[bathGas]]=newAdjList


####################Copy in, edit, write out new input file############################################################
    #Copy in the old inputFile:
    inputList=[]
    bath_gas_needed = False
    with open(os.path.join(args.oldWorkDir, 'input.py'), 'rb') as inputFile:
        for line in inputFile:
            inputList.append(line)
	    if re.search('pressureDependence', line):
		bath_gas_needed = True

    newInputList=[]
    newReactorsInputted=False
    dontCopy=False
    species_input_yet=False

    #Start editting the input file
    for line in inputList:
        #change the thermo libraries
        if re.search('thermoLibraries', line):
	    #Add aromatic libraries in addition to burner libraries if specified by user
	    if re.search('on', args.Reactor_2_Aromatic_libraries):
            	newLine=re.sub('\]', ",'" +args.name+ "', 'C10H11', 'C3', 'Fulvene_H', 'naphthalene_H', 'vinylCPD_H', 'biCPD_H_shift', 'Swamy_AR2']", line)
            	newInputList.append(newLine)
	    else:
		newLine=re.sub('\]', ",'" +args.name+ "']", line)
            	newInputList.append(newLine)
        
	#change the reaction libraries
        elif re.search('reactionLibraries', line):
	    #Add aromatic libraries in addition to burner libraries if specified by user
	    if re.search('on', args.Reactor_2_Aromatic_libraries):
            	newLine=re.sub('\]', "('" +args.name+ "',True), ('C10H11',False), ('C3',False), ('Fulvene_H',False), ('naphthalene_H',False), ('vinylCPD_H',False), ('biCPD_H_shift',False), ('Swamy_AR2', False), ('fascella',False), ('kislovB',False)]", line)
            	newInputList.append(newLine)
	    else:
		newLine=re.sub('\]', "('" +args.name+ "',True)]", line)
            	newInputList.append(newLine)
        
        #Include an exception for 'input species' if it doesn't already exist
        elif re.search('allowed', line):
            if not re.search('input species', line):
                newLine=re.sub('\]', ",'input species']", line)
                newInputList.append(newLine)
	    else:
		newInputList.append(line)

	#Allow singlet O2
        elif re.search('allowSingletO2', line):
            if not re.search('True', line):
                newLine=re.sub('False', "True", line)
                newInputList.append(newLine)
	    else:
		newInputList.append(line)

        #Change the species
        elif re.search('species\(', line):
            dontCopy=True
	    if not species_input_yet:

                species_input_yet=True
		additivematch = False
                for component in newStream:
                    newInputList.append("species(\n")
		    for additiveSpecies in stream2:
                        if component==additiveSpecies and not additivePresent[additiveSpecies]:
                            newInputList.append("    label='"+additiveName[additiveSpecies]+"',\n")
			    additivematch = True
			    break
		        else: additivematch = False
                    if not additivematch:
                        newInputList.append("    label='"+getSpeciesIdentifier(component)+"',\n")
                    if component in bathGases.values():
                        newInputList.append("    reactive=False,\n")
			bath_gas_needed = False
                    elif component.molecule[0].toSMILES() in seeded_aromatic_SMILES:
                        newInputList.append("    reactive=False,\n")
                    else:
                        newInputList.append("    reactive=True,\n")
                    newInputList.append('    structure=adjacencyList(\n        """\n')
                    newInputList.append(adjListDict[component])
                    newInputList.append('        """),\n')
                    newInputList.append(")\n\n")
		#Insert N2 as an inert bath gas if a bath gas is required for the P-dep and there is nothing already there
		if bath_gas_needed:
		    newInputList.append("species(\n")
		    newInputList.append("    label='N2',\n")
		    newInputList.append("    reactive=False,\n")
		    newInputList.append('    structure=adjacencyList(\n        """\n')
                    newInputList.append(adjListDict[bathGases['N2']])
                    newInputList.append('        """),\n')
                    newInputList.append(")\n\n")

		    #Add N2 with 0.0 initial concentration to new reactor block
		    newReactorBlock.insert(newReactorBlock.index("    initialMoleFractions={\n") + 1,"        'N2': 0.0,\n")

	    else:
	        continue

        #Otherwise just copy the same list
        elif re.search("simpleReactor\(", line):
            dontCopy=True
            if not newReactorsInputted:
                newReactorsInputted=True
                for line3 in newReactorBlock:
                    newInputList.append(line3)
	#Don't copy P-dep settings to Reactor 2 if user specifies P_dep to be off for that reactor
	elif re.search("pressureDependence\(", line) and re.search('off', args.Reactor_2_P_dep):    
	    dontCopy=True
	elif re.search("toleranceMoveToCore", line):
            newInputList.append("    toleranceMoveToCore=" + str(args.Reactor_2_tolmovetocore) + ",\n")
	else:
            if not dontCopy: newInputList.append(line)
        #Start copying again after reaching end of any input block
        if re.search("\)\n", line): dontCopy=False

    #Make a new directory and save new input.py
    if not os.path.exists(args.newWorkDir):
        os.makedirs(args.newWorkDir)
    with open(os.path.join(args.newWorkDir, 'input.py'), 'wb') as outFile:
        for line in newInputList:
            outFile.write(line)
