#!/usr/bin/env python
# encoding: utf-8

#Import a variety of things
import os
import sys
import logging
import re

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.kinetics import PDepArrhenius, PDepKineticsModel

from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.kinetics import KineticsDepository, KineticsRules
from rmgpy.qm.main import QMCalculator


# from model_kinetics import info

#####

# This wasn't really specified, but I think this is something required for Discovery
"""
if len(sys.argv)>1:
	i = int(sys.argv[-1])
elif os.getenv('LSB_JOBINDEX'):
	i = int(os.getenv('LSB_JOBINDEX'))
else:
	raise Exception("Specify a TS number!")
"""
#####

rxnFamilies = ['H_abstraction'] # Only looking at H_abstraction via OOH

print 'Loading RMG Database ...'
rmgDatabase = RMGDatabase()
rmgDatabase.load(os.path.abspath(os.path.join(os.getenv('RMGpy'), '..', 'RMG-database', 'input')), kineticsFamilies=rxnFamilies) # unsure if this is the right environment
print 'RMG Database Loaded'

def calculate(reaction):
	rxnFamily = reaction.family
	tsDatabase = rmgDatabase.kinetics.families[rxnFamily].transitionStates
	reaction = qmCalc.getKineticData(reaction, tsDatabase)
	for files in os.listdir('./'): # TBH, I don't know what this part does.
		if files.startswith('core'):
			os.remove(files)
	return reaction

#######################
"""
This section of code is designed to go through each file listed as a SMILES.txt
and create a dictionary that has the location of the smiles file an identifying
name to go along with it.
"""

SMILES_files = []
for root, dirs, files in os.walk("../RMG-models", topdown=False):
	for name in files:
		if str(os.path.join(root,name)).endswith('SMILES.txt'):
			SMILES_files.append(os.path.join(root,name))

smiles_dict = {}
for name in SMILES_files:
	try:
		e,d,c,b,a = name.split('/')
		smiles_entry = c + "_" + b
	except ValueError:
		d,c,b,a = name.split('/')
		smiles_entry = b
	smiles_dict[name] = smiles_entry

# print smiles_dict

############

for entry in SMILES_files:
	smiles_file = str(entry)


	known_smiles = {}
	known_names = []
	identified_labels = []
	line = None

	################################################

	with open(smiles_file) as f:
	    for line in f:
	        if not line.strip():
	            continue
	        try:
	            user = None
	            if '!' in line:
	                line, comments = line.split("!",1)
	                if comments:
	                    usermatch = re.match("\s+Confirmed by (.*)",comments)
	                    if usermatch:
	                        user = usermatch.group(1)
	            tokens = line.split()
	            assert len(tokens) == 2, "Not two tokens on line (was expecting NAME    SMILES)"
	            name, smiles = tokens
	            if name in known_smiles:
	            	assert smiles == known_smiles[name], "{0} defined twice, as {1} and {2}".format(name, known_smiles[name], smiles)
	            known_smiles[name] = smiles
	            if name not in known_names:
	                known_names.append(name)
	        except Exception as e:
	            logging.warning("Error reading line '{0}'".format(line))
	            raise e
	break
	#####
	allRxns = []
	for chemkinRxn in info:
	    for rFam in chemkinRxn['possibleReactionFamilies']:
	        if rFam in rxnFamilies:
	            if len(chemkinRxn['possibleReactionFamilies']) > 1:
	                chemkinRxn['possibleReactionFamilies'] = [rFam]
		    chemkin_rxn_name = chemkinRxn['chemkinKinetics'].split()[0]
		    new_rxn = True
		    for chem_rxn in allRxns:
		        if chem_rxn['chemkinKinetics'].strip().startswith(chemkin_rxn_name):
					new_rxn = False
					break
		    if new_rxn:
		        allRxns.append(chemkinRxn)


	chemkinRxn = allRxns[i-1]
	reaction = chemkinRxn['reaction']

	chemkinRxn = allRxns[i-1]
	reaction = chemkinRxn['reaction']
	try:
	    rcts, prdts = reaction.split('<=>')
	except ValueError:
	    rcts, prdts = reaction.split('=>')
	react = rcts.split('+')
	products = prdts.split('+')
	reactants = []

	### This is one of the specific things that I changed from Pierre's script. ###
	for r in react:
		if "OOH" in r:
			reactants.append(r) # Although, I think reactants is a dictionary? Or it's a list of lists? If it's a list of lists, I think this should work.
	################################################################################

	reactant_molecules = [Molecule().fromSMILES(known_smiles[key.strip()]) for key in reactants]
	product_molecules = [Molecule().fromSMILES(known_smiles[key.strip()]) for key in products]

	reactants_species = [Species(molecule=mol.generateResonanceIsomers()) for mol in reactant_molecules]
	product_species = [Species(molecule=mol.generateResonanceIsomers()) for mol in product_molecules]

	testReaction = Reaction(reactants = reactant_species, products = product_species, reversible = True)

	checkRxn = rmgDatabase.kinetics.generateReactionsFromFamilies(reactant_molecules, product_molecules, only_families=chemkinRxn['possibleReactionFamilies'])
	if checkRxn==[]:
		reactIsoms = [mol.generateResonanceIsomers() for mol in reactant_molecules]
		if len(reactIsoms)==2:
			r1, r2 = reactIsoms
			while checkRxn==[]:
				for moleculeA in r1:
					for moleculeB in r2:
						reactant_molecules = [moleculeA, moleculeB]
						checkRxn = rmgDatabase.kinetics.generateReactionsFromFamilies(reactant_molecules, product_molecules, only_families=chemkinRxn['possibleReactionFamilies'])

	reaction = checkRxn[0]

	assert testReaction.isIsomorphic(reaction)

	atLblsR = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])
	atLblsP = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])

	gotOne = False
	for reactant in reaction.reactants:
	    reactant = reactant.molecule[0]
	    reactant.clearLabeledAtoms()
	    for atom in reactant.atoms:
	        for atomLabel in reaction.labeledAtoms:
	            if atom==atomLabel[1]:
	                atom.label = atomLabel[0]
	                atLblsR[atomLabel[0]] = True

	for product in reaction.products:
	    product = product.molecule[0]
	    product.clearLabeledAtoms()
	    for atom in product.atoms:
	        for atomLabel in reaction.labeledAtoms:
	            if atom==atomLabel[1]:
	                atom.label = atomLabel[0]
	                atLblsP[atomLabel[0]] = True

	if all( atLblsR.values() ) and all( atLblsP.values() ):
		gotOne=True
	rxnFamily = reaction.family.label
	assert gotOne
	qmCalc = QMCalculator(
	                        software='gaussian',
	                        method='m062x',
	                        fileStore='/scratch/bhoorasingh.p/QMfiles',
	                        scratchDirectory='/scratch/bhoorasingh.p/QMscratch',
	                        )
	reaction = calculate(reaction)

	if reaction.kinetics:
		"""
		REWRITE SARATHYKIN STUFF FOR VARYING IMPORTED MECHANISMS # Fixed?
		"""
		"""
		Return the rate coefficient in the appropriate combination of cm^3,
		mol, and s at temperature `T` in K and pressure `P` in Pa.
		"""
		importKin = chemkinRxn['rmgPyKinetics']
		Temp = 1000 # Kelvin
		idx = str(i)
		row = [idx, rxnFamily]
		row.extend([mol.label for mol in reaction.reactants])
		row.extend([mol.label for mol in reaction.products])

		rateCal = reaction.calculateTSTRateCoefficient(Temp)
		row.extend(['AutoTST_fwd', str(rateCal)])

		if rxnFamily == 'R_Addition_MultipleBond':
			if reverseRxn:
				revRate = reaction.generateReverseRateCoefficient()
				rev_rateCal = revRate.getRateCoefficient(Temp)
				row.extend(['AutoTST_rev', str(rev_rateCal)])
			else:
				row.extend(['AutoTST_rev', 'not reverse'])

		if isinstance(importKin, PDepArrhenius) or isinstance(importKin, PDepKineticsModel):
			row.extend([smiles_dict[entry], str(importKin.getRateCoefficient(Temp, 1000000))]) #1000000 Pa = 10 bar
		else:
			row.extend([smiles_dict[entry], str(importKin.getRateCoefficient(Temp))])

		famDatabase = rmgDatabase.kinetics.families[rxnFamily]
		famDatabase.addKineticsRulesFromTrainingSet(thermoDatabase=rmgDatabase.thermo)
		famDatabase.fillKineticsRulesByAveragingUp()

		# Using the `testReaction` will keep the reaction direction the same as the model
		rxnTemplate = famDatabase.getReactionTemplate(testReaction)
		kList = famDatabase.getKinetics(testReaction, rxnTemplate)

		allRates = []
		kComments = []
		for rate in kList:
			k = rate[0]
			if k:
				label = rate[1]
				if isinstance(label, KineticsDepository):
					label = label.label
				elif isinstance(label, KineticsRules):
					label = label.label
				if label.lower()=='rate rules':
					kComments = k.comment.split('\n')
				allRates.append((label, k.getRateCoefficient(Temp)))

		kDict = {'Rate Rules': []}
		for kTuple in allRates:
			if kTuple[0]=='AutoTST':
				kDict['AutoTST'] = kDict['AutoTST'] + [kTuple[1]]
			elif kTuple[0]=='rate rules':
				kDict['Rate Rules'] = kDict['Rate Rules'] + [kTuple[1]]
			elif kTuple[0].endswith('/NIST'):
				kDict['NIST'] = kDict['NIST'] + [kTuple[1]]
			elif kTuple[0].endswith('/rules'):
				kDict['KineitcsRules'] = kDict['KineitcsRules'] + [kTuple[1]]

		# Print the values in order
		kineticsTypes = ['Rate Rules']#['AutoTST', 'Rate Rules', 'NIST', 'KineitcsRules']
		for kinType in kineticsTypes:
			row.append(kinType)
			if kDict[kinType]==[]:
				row.append('no Val')
			else:
				for val in kDict[kinType]:
					row.append(str(val))

		# Store line containing reaction from corresponding chemkin file
		row.append(smiles_dict[entry])
		row.append(chemkinRxn['chemkinKinetics'].strip())

		# Store rmgpy kinetics group comments
		row = row + kComments

		folderPath = os.path.join('KinTxtFiles', idx)

		if not os.path.exists(folderPath):
			os.makedirs(folderPath)

		input_string = ','.join(row)
		with open(os.path.join(folderPath, smiles_dict[entry] + '_kinetics.txt'), 'w') as kinTxt:
	            kinTxt.write(input_string)


"""
if checkRxn == []:
  reactIsoms = [mol.generateResonanceIsomers() for mol in reactant_molecules]
  if len(reactIsoms) ==  2:
    r1, r2 = reactIsoms
    while checkRxn == []:
      for moleculeA in r1:
        for moleculeB in r2:
          reactant_molecules = [moleculeA, moleculeB]
          checkRxn = rmgDatabase.kinetics.generateReactionsFromFamilies(reactant_molecules, product_molecules, only_families = chemkinRxn)

reaction = checkRxn[0]

assert testRaction.isIsomorphic(reaction) # not really sure what this does, ask Pierre

# This adds labeled atoms to the molecule and checks that all the labels were added.
atLblsR = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])
atLblsP = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])

gotOne = False
for reactant in reaction.reactants:
  reactant = reactant.molecule[0]
  reactant.clearLabeledAtoms()
  for atoms in reactant.atoms:
    for atomLabel in reaction.labeledAtoms:
      if atom == atomLabel[1]
      atom.label = atomLabel[0]
      atLblsR[atomLabel[0]] = True

for product in reaction.products:
  product = product.molecule[0]
  product.clearLabeledAtoms()
  for atoms in product.atoms:
    for atomLabel in reaction.labeledAtoms:
      if atom == atomLabel[1]
      atom.label = atomLabel[0]
      atLblsP[atomLabel[0]] = True

if all( atLblsR.values() ) and all( atLblsP.values() ):
  gotOne = True # Still confused about what this thing is

rxnFamily = reaction.family
assert gotOne

qmCalc = QMCalculator( software = 'gaussian', method = 'm062x', fileStore = '/scratch/harms.n/QMfiles', scratchDirectory = '/scratch/harms.n',)

reaction = calculate(reaction)

if reaction.kinetics:
  importerKin = chemkinRxn['rmgPyKinetics'] # May need to adjust this depending on what kinetics we load. Sarathy is from Pierre's last kinetic model thing, I need to pick a new variable here. Importer was selected.
  Temp = 1000 # Kelvin
  idx = str(i)
  row = [idx, rxnFamily]
  row.extend([mol.label for mol in reaction.reactants])
  row.extend([mol.label for mol in reaction.products])

  rateCal = reaction.calculateTSTRateCoefficient(Temp)
  row.extend(['AutoTST_fwd',str(rateCal)])

  if isinstance(importerKin, PDepArrhenius) or isinstance(importerKin, PDepKineticsModel):
    row.extend(['Importer', str(importerKin.getRateCoefficient(Temp, 1000000))])
  else:
    row.extend(['Importer', str(importerKin.getRateCoefficient(Temp))])

  famDatabase = rmgDatabase.kinetics.families[rxnFamily]
  famDatabase.addKineticsRulesFromTrainingSet(thermoDatabase=rmgDatabase.thermo)
	famDatabase.fillKineticsRulesByAveragingUp()

  rxnTemplate = famDatabase.getReactionTemplate(testReaction)
  kList = famDatabase.getKinetics(testReaction, rxnTemplate)

  # Sort the calculated rates bases on how they were determined.
  allRates = []
  kComments = []
  for rate in kList:
    k = rate [0]
    if k:
      label = rate[1]
      if isinstance(label, KineticsDepository):
        label = label.label
      elif is instance(label, KineticsRules):
        label = label.label
      if label.lower() == 'rate rules':
        kComments = k.comment.split('\n')
      allRates.append(label, k.getRateCoefficient(Temp))

  kDict = {'Rate Rules': []}
  for kTuple in allRates:
    if kTuple[0] == 'AutoTST':
      kDict['AutoTST'] = kDict['AutoTST'] + [kTuple[1]]
    if kTuple[0] == 'rate rules':
      kDict['Rate Rules'] = kDict['Rate Rules'] + [kTuple[1]]
    if kTuple[0].endswith('/NIST'):
      kDict['NIST'] = kDict['NIST'] + [kTuple[1]]
    if kTuple[0].endswith('/rules'):
      kDict['KineticsRules'] = kDict['KineticsRules'] + [kTuple[1]]

  # Print the values in order
  kineticsTypes = ['Rate Rules']#['AutoTST', 'Rate Rules', 'NIST', 'KineitcsRules']
  for kinType in kineticsTypes:
    row.append(kinType)
    if kDict[kinType] == []:
      row.append('no Val')
    else:
      for val in kDict[kinType]:
        row.append(str(val))

  # Store line containing reaction from source of intrest's chemkin file.
  row.append('ImporterSource')
  row.append(chemkinRxn['chemkinKinetics'].strip())

  # Store rmgpy kinetics group comments
  row = row + kComments

  folderPath = os.path.join('KinTxtFiles', idx)

  if not os.path.exists(folderPath):
    os.makedirs(folderPath)

  input_string = ','.join(row)
  with open(os.path.join(folderPath, 'kinetic.txt'), 'w') as kinTxt:
    kinTxt.write(input_string)
"""
