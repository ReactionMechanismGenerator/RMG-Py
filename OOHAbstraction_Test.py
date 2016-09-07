#!/usr/bin/env python
# encoding: utf-8

#Import a variety of things
import os
import sys
import logging
import re
import imp


from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.kinetics import PDepArrhenius, PDepKineticsModel

from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.kinetics import KineticsDepository, KineticsRules
from rmgpy.qm.main import QMCalculator


#####
# /gss_gpfs_scratch/harms.n


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

rxnFamilies = ['H_Abstraction'] # Only looking at H_abstraction via OOH



print 'Loading RMG Database ...'
rmgDatabase = RMGDatabase()
rmgDatabase.load(os.path.abspath(os.path.join(os.getenv('RMGpy'), '..', 'RMG-database', 'input')), kineticsFamilies=rxnFamilies) # unsure if this is the right environment
print 'RMG Database Loaded'

qmCalc = QMCalculator(
						software='gaussian',
						method='m062x',
						fileStore='/scratch/harms.n/QMfiles',
						scratchDirectory='/scratch/harms.n/QMscratch',
						)

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


importer_files = []
for root, dirs, files in os.walk("../RMG-models", topdown=False):
	smiles_file = None
	kinetics_file = None
	for name in files:
		if str(os.path.join(root,name)).endswith('SMILES.txt'):
			smiles_file = os.path.join(root,name)
		if str(os.path.join(root,name)).endswith('model_kinetics.py'):
			kinetics_file = os.path.join(root,name)
	if kinetics_file and smiles_file:
		importer_files.append((smiles_file, kinetics_file))
	else:
		"""
		print root
		if smiles_file is None:
			print "Missing SMILES file"
		if kinetics_file is None:
			print "Missing kinetics file"
		"""

############

for item in importer_files:
	smiles_file, kinetics_file = item
	try:
		a,b,c,d,e = smiles_file.split('/')
		txt_file = c + '_' + d
	except ValueError:
		a,b,c,d = smiles_file.split('/')
		txt_file = c
	print txt_file

	break

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
	#break #do I need this break here? This might be why the script stops so abruptly...?
	#####
	allRxns = []

	#####
	# This section obtains the info from the appropriate kinetics file. Currently labled "model_kinetics.py".
	load = imp.load_source('info', kinetics_file)
	info = load.info
	#####

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

	# allRxns seems to be a dictionary... I'm a little drunk at the moment, so idk if I'll be able to decipher it rn


	# Now this portion is throwing errors ####!!!!
	reactions = [d["reaction"] for d in allRxns]

	reactant = []
	product = []

	for reaction in reactions:
		try:
		    rcts, prdts = reaction.split('<=>')
		except ValueError:
		    rcts, prdts = reaction.split('=>')
		reactants = rcts.split('+') # A list of reactants for a given rxn
		products = prdts.split('+') # A list of products for a given rxn

		# Obtains a list of SMILES structures for a reaction.
		reactant_molecules = [(known_smiles[key.strip()]) for key in reactants]
		product_molecules = [(known_smiles[key.strip()]) for key in products]



		### This section is designed to determine if OOH is in the reactants and if OO is in the products
		# i.e. if this is OOH abstraction

		reactant_mol = []
		product_mol = []
		if "[O]O" in reactant_molecules:
			if "OO" in product_molecules:
				reactant_mol = [Molecule().fromSMILES(key) for key in reactant_molecules]
				product_mol = [Molecule().fromSMILES(key) for key in product_molecules]

		if reactant_mol != [] and product_mol != []:


			reactant_species = [Species(molecule=mol.generateResonanceIsomers()) for mol in reactant_mol]
			product_species = [Species(molecule=mol.generateResonanceIsomers()) for mol in product_mol]

			testReaction = Reaction(reactants = reactant_species, products = product_species, reversible = True)
			checkRxn = rmgDatabase.kinetics.generateReactionsFromFamilies(reactant_mol, product_mol, only_families=rxnFamilies)
			#print [testReaction]
			#print checkRxn, type(checkRxn), len(checkRxn)


			# This part was in Pierre's script, but I have no use for it (I think)

			if checkRxn==[]:
				reactIsoms = [mol.generateResonanceIsomers() for mol in reactant_mol]
				if len(reactIsoms)==2:
					r1, r2 = reactIsoms
					while checkRxn==[]:
						for moleculeA in r1:
							for moleculeB in r2:
								reactant_mol = [moleculeA, moleculeB]
								checkRxn = rmgDatabase.kinetics.generateReactionsFromFamilies(reactant_mol, product_mol, only_families=chemkinRxn['possibleReactionFamilies'])

			reaction = checkRxn[0]
			print reaction

			assert testReaction.isIsomorphic(reaction)

			atLblsR = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])
			atLblsP = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])

			gotOne = False
			for reactant in reaction.reactants:
			    #reactant = reactant.molecule[0]
			    reactant.clearLabeledAtoms()
			    for atom in reactant.atoms:
			        for atomLabel in reaction.labeledAtoms:
			            if atom==atomLabel[1]:
			                atom.label = atomLabel[0]
			                atLblsR[atomLabel[0]] = True

			for product in reaction.products:
			    #product = product.molecule[0]
			    product.clearLabeledAtoms()
			    for atom in product.atoms:
			        for atomLabel in reaction.labeledAtoms:
			            if atom==atomLabel[1]:
			                atom.label = atomLabel[0]
			                atLblsP[atomLabel[0]] = True

			if all( atLblsR.values() ) and all( atLblsP.values() ):
				gotOne=True
			rxnFamily = reaction.family
			assert gotOne

			######### NO GUARENTEES IF ANYTHING AFTER THIS LINE WORKS

			reaction = calculate(reaction)

			if reaction.kinetics:
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

				"""
				# No reaction will be of this reaction family b/c H_Abstraction
				if rxnFamily == 'R_Addition_MultipleBond':
					if reverseRxn:
						revRate = reaction.generateReverseRateCoefficient()
						rev_rateCal = revRate.getRateCoefficient(Temp)
						row.extend(['AutoTST_rev', str(rev_rateCal)])
					else:
						row.extend(['AutoTST_rev', 'not reverse'])
				"""
				if isinstance(importKin, PDepArrhenius) or isinstance(importKin, PDepKineticsModel):
					row.extend([txt_file, str(importKin.getRateCoefficient(Temp, 1000000))]) #1000000 Pa = 10 bar
				else:
					row.extend([txt_file, str(importKin.getRateCoefficient(Temp))])

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
				row.append(txt_file)
				row.append(chemkinRxn['chemkinKinetics'].strip())

				# Store rmgpy kinetics group comments
				row = row + kComments

				folderPath = os.path.join('KinTxtFiles', idx)

				if not os.path.exists(folderPath):
					os.makedirs(folderPath)

				input_string = ','.join(row)
				with open(os.path.join(folderPath, txt_file + '_kinetics.txt'), 'a') as kinTxt:
			            kinTxt.write(input_string)
