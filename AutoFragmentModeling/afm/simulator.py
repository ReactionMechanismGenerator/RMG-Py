
# one iteration of the Monte Carlo simulation
# consists of the following steps
# s1: decide the time step to choose,
# need k, count of the related fragments
# s2: cast random number to decide which
# reaction rule to fire, need to create a
# dict which has <reaction, rate> pairs
# s3: realize the reaction rule on fragment level
# update fragment_count_dict
# s4: realize the reaction rule on molecule level
# need dict which has <fragment, [mol1, mol2, ...]> pairs
# update molecule list: some molecules to remove,
# some molecules to append.
import random
import numpy as np
import pandas as pd

import rmgpy.constants
from rmgpy.chemkin import loadChemkinFile

import afm.loader
import afm.utils
from afm.canteraModel import Cantera, CanteraCondition

class Simulator(object):

	def __init__(self, chemkin_path, dictionary_path, fragment_smiles_path):

		self.load_fragment_chemistry(chemkin_path, dictionary_path, fragment_smiles_path)

	def load_fragment_chemistry(self, chemkin_path, dictionary_path, fragment_smiles_path):

		fragment_dict, fragment_rxns = afm.loader.load_fragment_reactions_from_chemkin(chemkin_path,
                                        											   dictionary_path,
                                        											   fragment_smiles_path)

	#	pseudo_fragrxns = afm.loader.load_pseudo_fragment_reactions(fragment_dict)

	#	self.fragment_reaction_list = fragment_rxns + pseudo_fragrxns
		self.fragment_reaction_list = fragment_rxns
		self.fragment_dict = fragment_dict

class OdeSimulator(Simulator):

	########################
	# Construction section #
	########################
	def __init__(self, 
				chemkin_path, 
				dictionary_path,
				fragment_smiles_path,
				temperature,
				pressure,
				outputDirectory='temp'):
		super(OdeSimulator, self).__init__(chemkin_path, 
										   dictionary_path,
										   fragment_smiles_path)

		speciesList, reactionList = loadChemkinFile(chemkin_path, dictionary_path)
		self.speciesList = speciesList
		self.reactionList = reactionList
		self.temperature = temperature # unit: K
		self.pressure = pressure # unit: bar
		self.outputDirectory = outputDirectory

	def simulate(self, initial_mol_fraction, termination_time):

		cantera_job = Cantera(speciesList=self.speciesList, 
							  reactionList=self.reactionList, 
							  outputDirectory=self.outputDirectory)

		cantera_job.loadModel()

		# Generate the conditions
		reactorTypeList = ['IdealGasConstPressureTemperatureReactor']

		speciesDict = {}
		for spe in self.speciesList:
			speciesDict[spe.label] = spe

		# prepare mole fraction for simulation condition
		molFracDict = {}
		total_mol_frac = 0
		for spe_label, init_mol_frac in initial_mol_fraction.iteritems():
			spe = speciesDict[spe_label]
			molFracDict[spe] = init_mol_frac
			total_mol_frac += init_mol_frac

		# normalize initial mole fractions
		for spe in molFracDict:
			molFracDict[spe] = molFracDict[spe]/total_mol_frac

		Tlist = ([self.temperature],'K')
		Plist = ([self.pressure],'bar')
		reactionTimeList = ([termination_time], 's')
		cantera_job.generateConditions(reactorTypeList, reactionTimeList, [molFracDict], Tlist, Plist)
		
		# Simulate
		alldata = cantera_job.simulate()

		return alldata

	def reattach_fragments(self, 
						   r_moles, 
						   l_moles, 
						   r_l_moles, 
						   rr_ll_list, 
						   grind_size=10, 
						   shuffle_seed=0):

		"""
		Given simulation data, alldata, which contains concentrations
		of all the fragments, return a list of merged molecules 
		after re-attachment with their concentrations.
		"""
		# cut large moles into smaller pieces
		grinded_r_moles = afm.utils.grind(r_moles, grind_size)
		grinded_l_moles = afm.utils.grind(l_moles, grind_size)

		# random shuffle
		r_moles_shuffle = afm.utils.shuffle(grinded_r_moles, shuffle_seed)
		l_moles_shuffle = afm.utils.shuffle(grinded_l_moles, shuffle_seed)

		# match concentrations for single-labeled fragments
		# including RCCCCR
		matches0 = afm.utils.match_concentrations_with_same_sums(l_moles_shuffle, 
													   			 r_moles_shuffle, 
													   			 diff_tol=1e-3)

		matches1, new_r_l_moles = afm.utils.matches_resolve(matches0, rr_ll_list)

		r_l_moles.extend(new_r_l_moles)
		# insert double-labeled fragments into matches
		# e.g., LCCCCR
		matches = afm.utils.match_concentrations_with_different_sums(matches1, r_l_moles)

		return matches

	def reattach_fragments_1_label(self,
						   r_moles,
						   rr_list,
						   grind_size=1,
						   shuffle_seed=0):

		"""
		Given simulation data, alldata, which contains concentrations
		of all the fragments, return a list of merged molecules
		after re-attachment with their concentrations.
		"""
		# cut large moles into smaller pieces
		grinded_r_moles = afm.utils.grind(r_moles, grind_size)

		# random shuffle
		r_moles_shuffle = afm.utils.shuffle(grinded_r_moles, shuffle_seed)

		# match concentrations for single-labeled fragments
		# including RCCCCR
		if len(r_moles_shuffle) % 2 != 0:
			half_length = (len(r_moles_shuffle) - 1) / 2
		else:
			half_length = len(r_moles_shuffle) / 2

		half_r_moles_shuffle_1 = r_moles_shuffle[0:half_length]
		half_r_moles_shuffle_2 = r_moles_shuffle[half_length:len(r_moles_shuffle)]

		matches0 = afm.utils.match_concentrations_with_same_sums(half_r_moles_shuffle_1,
																		 half_r_moles_shuffle_2,
																		 diff_tol=1e-3)

		matches1, new_r_l_moles = afm.utils.matches_resolve_1_label(matches0, rr_list)

		r_r_moles = []
		r_r_moles.extend(new_r_l_moles)
		# insert double-labeled fragments into matches
		# e.g., LCCCCR

		matches_1_label = afm.utils.match_concentrations_with_different_sums(matches1, r_r_moles)

		return matches_1_label

	def get_molecular_weight_distribution(self, 
										  alldata, 
										  grind_size=10, 
										  shuffle_seed=0):

		# prepare moles data
		_, dataList, _ = alldata[0]
		TData = dataList[0]
		PData = dataList[1]
		VData = dataList[2]
		total_moles = PData.data*VData.data/8.314/TData.data

		moles_dict = {}
		for data in dataList[3:]:
			spe_label = data.label
			moles_dict[spe_label] = max(data.data[-1]*total_moles[-1],0)

		# prepare moles data for re-attachment
		r_moles, l_moles, r_l_moles, remain_moles, rr_ll_list = categorize_fragments(moles_dict)

		matches = self.reattach_fragments(r_moles, 
										  l_moles, 
										  r_l_moles, 
										  rr_ll_list, 
										  grind_size, 
										  shuffle_seed)

		flattened_matches = [(tuple(afm.utils.flatten(m[0])), m[1]) for m in matches]

		final_frags_moles = []
		for remain in remain_moles:
			label, val = remain
			final_frags_moles.append(((label, ), val))

		final_frags_moles.extend(flattened_matches)

		# calculate fragmental weight distribution
		fragmental_weight_distri = []
		for final_frag_mole in final_frags_moles:
			sub_frag_labels, mole = final_frag_mole
			total_frag_weight = 0
			for sub_frag_label in sub_frag_labels:
				sub_frag = self.fragment_dict[sub_frag_label]
				total_frag_weight += sub_frag.getMolecularWeight()
			
			fragmental_weight_distri.append((total_frag_weight, mole))

		return fragmental_weight_distri

# 1 label (R label only)
	def get_molecular_weight_distribution_1_label(self,
										         alldata,
										         grind_size=1,
								        	     shuffle_seed=0):

		# prepare moles data
		_, dataList, _ = alldata[0]
		TData = dataList[0]
		PData = dataList[1]
		VData = dataList[2]
		total_moles = PData.data * VData.data / 8.314 / TData.data

		moles_dict = {}
		for data in dataList[3:]:
			spe_label = data.label
			moles_dict[spe_label] = max(data.data[-1] * total_moles[-1], 0)

		# prepare moles data for re-attachment
		r_moles, remain_moles, rr_list = categorize_fragments_1_label(moles_dict)

		matches = self.reattach_fragments_1_label(r_moles,
										  rr_list,
										  grind_size,
										  shuffle_seed)

		flattened_matches = [(tuple(afm.utils.flatten(m[0])), m[1]) for m in matches]

		final_frags_moles = []
		for remain in remain_moles:
			label, val = remain
			final_frags_moles.append(((label,), val))

		final_frags_moles.extend(flattened_matches)

		# calculate fragmental weight distribution
		fragmental_weight_distri_1_label = []
		for final_frag_mole in final_frags_moles:
			sub_frag_labels, mole = final_frag_mole
			total_frag_weight = 0
			for sub_frag_label in sub_frag_labels:
				sub_frag = self.fragment_dict[sub_frag_label]
				total_frag_weight += sub_frag.getMolecularWeight()

			fragmental_weight_distri_1_label.append((total_frag_weight, mole))

		return fragmental_weight_distri_1_label


def categorize_fragments(moles_dict):

	r_moles = []
	l_moles = []
	r_l_moles = []
	rr_ll_list = []
	remain_moles = []
	for spe_label in moles_dict:
		if '*' in spe_label:
			remain_moles.append((spe_label, moles_dict[spe_label]))
			continue
		if abs(moles_dict[spe_label]) <= 1e-6:
			remain_moles.append((spe_label, moles_dict[spe_label]))
			continue
		
		r_count = spe_label.count('R')
		l_count = spe_label.count('L')
		label_count = r_count + l_count
		
		if label_count == 0:
			remain_moles.append((spe_label, moles_dict[spe_label]))
			continue

		if label_count == 1:
			if r_count == 1:
				r_moles.append((spe_label, moles_dict[spe_label]))
			elif l_count == 1:
				l_moles.append((spe_label, moles_dict[spe_label]))
		elif label_count == 2 and l_count == r_count:
			r_l_moles.append((spe_label, moles_dict[spe_label]))
		elif label_count == 2 and l_count != r_count:
			rr_ll_list.append(spe_label)
			if r_count == 2:
				r_moles.append((spe_label, 2*moles_dict[spe_label]))
			elif l_count == 2:
				l_moles.append((spe_label, 2*moles_dict[spe_label]))
		else:
			print "{0} has more than two cutting labels, which is not supported.".format(spe_label)

	return r_moles, l_moles, r_l_moles, remain_moles, rr_ll_list

# 1 label (R label only)
def categorize_fragments_1_label(moles_dict):
	r_moles_1_label = []
	rr_list_1_label = []
	remain_moles_1_label = []
	for spe_label in moles_dict:
		if '*' in spe_label:
			remain_moles_1_label.append((spe_label, moles_dict[spe_label]))
			continue
		if abs(moles_dict[spe_label]) <= 1e-6:
			remain_moles_1_label.append((spe_label, moles_dict[spe_label]))
			continue

		r_count = spe_label.count('R')

		label_count = r_count

		if label_count == 0:
			remain_moles_1_label.append((spe_label, moles_dict[spe_label]))
			continue

		if label_count == 1:
			if r_count == 1:
				r_moles_1_label.append((spe_label, moles_dict[spe_label]))

		elif label_count == 2:
			rr_list_1_label.append(spe_label)
			if r_count == 2:
				r_moles_1_label.append((spe_label, 2 * moles_dict[spe_label]))

		else:
			print "{0} has more than two cutting labels, which is not supported.".format(spe_label)

	return r_moles_1_label, remain_moles_1_label, rr_list_1_label


class MonteCarloSimulator(Simulator):

	########################
	# Construction section #
	########################
	def __init__(self, 
				chemkin_path, 
				dictionary_path,
				fragment_smiles_path,
				initial_molecules,
				volume, 
				temperature):
		super(MonteCarloSimulator, self).__init__(chemkin_path, 
												  dictionary_path,
												  fragment_smiles_path)

		self.initialize_fragment_counts(initial_molecules)

		self.initialize_molecule_fragment_df(initial_molecules)
		self.volume = volume # unit: m^3
		self.temperature = temperature
		self.reaction_flux_array = np.zeros(len(self.fragment_reaction_list))

	def initialize_fragment_counts(self, initial_molecules):

		# initialize the fragment_count
		self.fragment_count_dict = dict(zip(self.fragment_dict.keys(), 
											[0]*len(self.fragment_dict)))

		for mol in initial_molecules:
			for frag in mol.composition:
				self.fragment_count_dict[frag] += mol.composition[frag]

	def initialize_molecule_fragment_df(self, initial_molecules):

		self.molecule_fragment_df = pd.DataFrame(columns=self.fragment_count_dict.keys())

		for mol in initial_molecules:
			insert_row = {}
			for fragment_label  in self.fragment_count_dict:
				if fragment_label in mol.composition:
					insert_row[fragment_label] = mol.composition[fragment_label]
				else:
					insert_row[fragment_label] = 0

			if not self.molecule_fragment_df.index.empty:
				insert_index = 1 + max(self.molecule_fragment_df.index)
			else:
				insert_index = 0
			self.molecule_fragment_df.loc[insert_index] = insert_row


	######################
	# Simulation section #
	######################
	def calculate_unimolucular_rate(self, reaction):

		k_u = reaction.kinetics.getRateCoefficient(self.temperature)
		frag_label = reaction.reactants[0].label
		frag_count = self.fragment_count_dict[frag_label]
		rate_u = k_u * frag_count # unit: 1/s

		return rate_u

	def calculate_bimolucular_rate(self, reaction):

		Na = rmgpy.constants.Na
		k_b = reaction.kinetics.getRateCoefficient(self.temperature)/Na # unit: m^3/s
		frag_label1 = reaction.reactants[0].label
		frag_label2 = reaction.reactants[1].label

		frag_count1 = self.fragment_count_dict[frag_label1]
		frag_count2 = self.fragment_count_dict[frag_label2]
		rate_b = k_b * frag_count1 * frag_count2 / self.volume # unit: 1/s

		return rate_b

	def update_reaction_fluxes(self):

		for idx, frag_rxn in enumerate(self.fragment_reaction_list):
			# unimolecular reactions
			if len(frag_rxn.reactants) == 1:
				rate = self.calculate_unimolucular_rate(frag_rxn)

			# bimolecular reactions
			elif len(frag_rxn.reactants) == 2:
				rate = self.calculate_bimolucular_rate(frag_rxn)

			self.reaction_flux_array[idx] = rate

	def time_step(self):

		total_rate = np.sum(self.reaction_flux_array)

		return 1.0/total_rate # unit: s


	def random_reaction_selection(self, seed=None):

		# cast a random number
		random.seed(seed)
		rand_num = random.uniform(0, np.sum(self.reaction_flux_array))

		# pick the one fragment reaction
		# indicated by the random number
		current_fluxsum = 0
		for idx, flux in enumerate(self.reaction_flux_array):

			current_fluxsum += flux
			print current_fluxsum
			if current_fluxsum >= rand_num:
				return idx

	def update_fragment_counts(self, reaction_idx):

		reaction = self.fragment_reaction_list[reaction_idx]
		for reactant in reaction.reactants:
			self.fragment_count_dict[reactant.label] -= 1

		for product in reaction.products:
			self.fragment_count_dict[product.label] += 1

