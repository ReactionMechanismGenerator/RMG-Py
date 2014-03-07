import os
import sys
from collections import defaultdict, Counter

import rmgpy
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.qm.main import QMCalculator
from rmgpy.data.kinetics import KineticsFamily, ReactionRecipe
from rmgpy.qm.gaussian import GaussianTSB3LYP
from rmgpy.qm.mopac import MopacTSPM7

# script to prep ts structures
actions = [
			['BREAK_BOND', '*1', 'S', '*2'],
			['FORM_BOND', '*2', 'S', '*3'],
			['GAIN_RADICAL', '*1', '1'],
			['LOSE_RADICAL', '*3', '1']
			]

family = 'H_Abstraction'

reactRecipe = ReactionRecipe(actions)

template = KineticsFamily(forwardRecipe=reactRecipe)

trusted = open(os.path.join(os.getenv('HOME'),'Code/RMG-database/input/kinetics/families/H_Abstraction/NIST.py'))

lines = trusted.readlines()
k = 0
idx = 0
reactants1 = defaultdict(list)
reactants2 = defaultdict(list)
for line in lines:
	k += 1
	if line.startswith('    reactant1 ='):
		idx += 1
		for num in range(k, k+100):
			if lines[num].find('{') != -1 or lines[num].find('*') != -1:
				reactants1[idx].append(lines[num])
			elif lines[num].find(',') != -1:
				break
	elif line.startswith('    reactant2 ='):
		for num in range(k, k+100):
			if lines[num].find('{') != -1 or lines[num].find('*') != -1:
				reactants2[idx].append(lines[num])
			elif lines[num].find(',') != -1:
				break

prevReactions = list()
tsStructures = list()
for idx in range(1, len(reactants1) + 1):
	r1 = ''
	r2 = ''
	for line in reactants1[idx]:
		r1 = r1 + line
	for line in reactants2[idx]:
		r2 = r2 + line
	r1 = Molecule().fromAdjacencyList(r1, saturateH=True)
	r2 = Molecule().fromAdjacencyList(r2, saturateH=True)
	rStruct = [r1, r2]
	pStruct, tsStruct = template.applyRecipe(rStruct, getTS=True)
	rxnInChI = [rStruct[0].toInChI(), rStruct[1].toInChI(), pStruct[0].toInChI(), pStruct[1].toInChI()]
	doubleChk = 0
	for pair in prevReactions:
		if Counter(pair) == Counter(rxnInChI):
			doubleChk = 1
	if doubleChk == 0:
		prevReactions.append(rxnInChI)
		tsStructures.append(tsStruct)

quantumMechanics = QMCalculator()
quantumMechanics.settings.software = 'mopac'
quantumMechanics.settings.fileStore = 'QMfiles'
quantumMechanics.settings.scratchDirectory = 'scratch'
quantumMechanics.settings.onlyCyclics = False
quantumMechanics.settings.maxRadicalNumber = 0

########################################################################################    
def fixSortLabel(molecule):
	"""
	This may not be required anymore. Was needed as when molecules were created, the
	rmg sorting labels would be set after where we tried to generate the TS.
	"""
	sortLbl = 0
	for vertex in molecule.vertices:
		vertex.sortingLabel = sortLbl
		sortLbl += 1
	return molecule

def calculate(TS):

	reactant = fixSortLabel(TS[0])
	product = fixSortLabel(TS[1])

	TS = [reactant, product]

	reaction = Reaction(label='H_Abstraction', reactants=reactant.split(), products=product.split(), reversible=True)

	qmReaction = MopacTSPM7(reaction, quantumMechanics.settings)
	qmReaction.generateTSGeometryDoubleEnded(doubleEnd=TS)
	# mopac, fromDbl, labels, notes = qmReaction.generateTSGeometryDoubleEnded(doubleEnd=TS)

	# if mopac:
	# 	import shutil
	# 	quantumMechanics.settings.software = 'gaussian'
	# 	gaussian = GaussianTSB3LYP(reaction, quantumMechanics.settings)
	# 	gaussian.geometry = fromDbl
	# 	gaussian.writeInputFile(1, fromSddl=fromDbl.molecule.getRadicalCount())
	# 	converged, internalCoord = gaussian.run()
	# 	shutil.copy(gaussian.outputFilePath, gaussian.outputFilePath+'.TS1.log')
	# 	
	# 	if internalCoord and not converged:
	# 		print "Internal coordinate error, trying in cartesian"
	# 		gaussian.writeInputFile(2, fromQST2=True)
	# 		converged, internalCoord = gaussian.run()
	# 	
	# 	if converged:
	# 		notes = notes + 'Transition state converged\n'
	# 		if not os.path.exists(gaussian.ircOutputFilePath):
	# 			gaussian.writeIRCFile()
	# 			rightTS = gaussian.runIRC()
	# 		else:
	# 			rightTS = gaussian.verifyIRCOutputFile()
	# 		if rightTS:
	# 			gaussian.writeRxnOutputFile(labels)
	# 		else:
	# 			notes = notes + 'IRC failed\n'
	# 	else:
	# 		notes = notes + 'Transition state not converged\n'
	# 		
	# with open(os.path.join(quantumMechanics.settings.fileStore, qmReaction.uniqueID + '.error'), 'w') as errorFile:
	# 	errorFile.write(notes)

########################################################################################

# if len(sys.argv)>1:
# 	i = int(sys.argv[1])
# elif os.getenv('LSB_JOBINDEX'):
# 	i = int(os.getenv('LSB_JOBINDEX'))
# else:
# 	raise Exception("Specify a TS number!")
# 
# quantumMechanics.settings.fileStore = os.path.join('QMfiles',str(i))
# if not os.path.exists(quantumMechanics.settings.fileStore):
# 	os.makedirs(quantumMechanics.settings.fileStore)
# 
# TS = tsStructures[i-1]
# calculate(TS)

for i, TS in enumerate(tsStructures):
	print "*"*70
	print "Reaction number {0} of {1}".format(i, len(tsStructures))
	quantumMechanics.settings.fileStore = os.path.join('QMfiles',str(i))
	if not os.path.exists(quantumMechanics.settings.fileStore):
		os.makedirs(quantumMechanics.settings.fileStore)
	calculate(TS)
