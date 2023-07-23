#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
This script is to help users transition from RMG-Py 2.4.x to RMG-Py 3.x by
automatically detecting variables, functions, attributes, methods, and arguments
which have been renamed. It works by using basic regex searches, so it is not
syntax aware and should therefore be used with caution. However, one benefit to
this approach is that this script can be used on any text file, including
IPython notebooks.

There are a number of available options. The only required arguments are the
files to analyze and the stage of changes to apply. The name changes are roughly
categorized into stage 1 with commonly used functionality and stage 2 with less
commonly used or private names. In addition, there is a list of dangerous names
which have increased risk of resulting in undesired replacements. They are
not replaced by default, but can be included using the `-x` argument.

By default, names are matched as individual words, i.e. surrounded by non-word
characters. Additionally, argument replacements check for an equal sign after
the name, while attribute and method replacements check for a period before
the name. All of these additional requirements can be disabled if it is desired
to do a naive search for all matches, although it is not recommended.

Also by default, only a diff of the potential changes are displayed to stdout.
The changes can be written to file by providing the `-w` argument. A backup of
the file is saved by default if writing to file, which can be disabled by
providing the `-n` argument.

A suggested option set for transitioning code which depends on RMG is

    python rmg2to3.py -0wn filename

which is equivalent to

    python rmg2to3.py --both-stages --write --no-backups filename

The script also accepts multiple filenames by manual specification, e.g.

    python rmg2to3.py -0wn file1 file2 file3

or using glob patterns, e.g.

    python rmg2to3.py -0wn *.py

Always be sure to double-check the result!
"""

import argparse
import difflib
import os
import re
import shutil
import sys

from tqdm import tqdm


# Module names
MODULES = {
    "canteraModel": "canteramodel",
    "canteraTest": "canteramodelTest",
    "extractInfoFromckcsv": "ckcsvparser",
    "diff_models": "diffmodels",
    "diff_modelsTest": "diffmodelsTest",
    "fluxtest": "fluxdiagramTest",
    "generate_reactions": "generatereactions",
    "testGenerateReactions": "generatereactionsTest",
    "merge_models": "mergemodels",
    "merge_modelsTest": "mergemodelsTest",
    "observablesRegression": "observablesregression",
}

# Global variables and functions
GLOBALS1 = {
    # Arkane
    "jobList": "job_list",
    "transitionStateDict": "transition_state_dict",
    # rmgpy.data.base
    "makeLogicNode": "make_logic_node",
    "getAllCombinations": "get_all_combinations",
    # rmgpy.data.rmg
    "getDB": "get_db",
    # rmgpy.data.solvation
    "saveEntry": "save_entry",
    # rmgpy.data.thermo
    "findCp0andCpInf": "find_cp0_and_cpinf",
    # rmgpy.data.kinetics.family
    "informationGain": "information_gain",
    "getObjectiveFunction": "get_objective_function",
    # rmgpy.data.kinetics.rules
    "removeIdenticalKinetics": "remove_identical_kinetics",
    "getTemplateLabel": "get_template_label",
    # rmgpy.kinetics.model
    "getRateCoefficientUnitsFromReactionOrder": "get_rate_coefficient_units_from_reaction_order",
    "getReactionOrderFromRateCoefficientUnits": "get_reaction_order_from_rate_coefficient_units",
    # rmgpy.kinetics.arrhenius
    "getw0": "get_w0",
    "getw0s": "get_w0s",
    # rmgpy.kinetics.diffusionLimited
    "diffusionLimiter": "diffusion_limiter",
    # rmgpy.molecule.adjlist
    "fromAdjacencyList": "from_adjacency_list",
    "toAdjacencyList": "to_adjacency_list",
    # rmgpy.molecule.atomtype
    "getFeatures": "get_features",
    "getAtomType": "get_atomtype",
    # rmgpy.molecule.converter
    "toRDKitMol": "to_rdkit_mol",
    "fromRDKitMol": "from_rdkit_mol",
    "toOBMol": "to_ob_mol",
    "fromOBMol": "from_ob_mol",
    # rmgpy.molecule.element
    "getElement": "get_element",
    "BDE_elements": "bde_elements",
    "BDEDict": "bde_dict",
    "BDEs": "bdes",
    # rmgpy.molecule.graph
    "getVertexConnectivityValue": "get_vertex_connectivity_value",
    "getVertexSortingLabel": "get_vertex_sorting_label",
    # rmgpy.molecule.symmetry
    "calculateAtomSymmetryNumber": "calculate_atom_symmetry_number",
    "calculateBondSymmetryNumber": "calculate_bond_symmetry_number",
    "calculateAxisSymmetryNumber": "calculate_axis_symmetry_number",
    "calculateCyclicSymmetryNumber": "calculate_cyclic_symmetry_number",
    # rmgpy.molecule.util
    "retrieveElementCount": "get_element_count",
    # rmgpy.pdep.cse
    "applyChemicallySignificantEigenvaluesMethod": "apply_chemically_significant_eigenvalues_method",
    # rmgpy.pdep.me
    "generateFullMEMatrix": "generate_full_me_matrix",
    # rmgpy.pdep.msc
    "applyModifiedStrongCollisionMethod": "apply_modified_strong_collision_method",
    # rmgpy.pdep.re
    "applyReservoirStateMethod": "apply_reservoir_state_method",
    # rmgpy.pdep.reaction
    "calculateMicrocanonicalRateCoefficients": "calculate_microcanonical_rate_coefficient",
    "applyRRKMTheory": "apply_rrkm_theory",
    "applyInverseLaplaceTransformMethod": "apply_inverse_laplace_transform_method",
    "fitInterpolationModel": "fit_interpolation_model",
    # rmgpy.rmg.main
    "initializeLog": "initialize_log",
    "get_condaPackage": "get_conda_package",
    "processProfileStats": "process_profile_stats",
    "makeProfileGraph": "make_profile_graph",
    # rmgpy.rmg.input
    "setGlobalRMG": "set_global_rmg",
    "readInputFile": "read_input_file",
    "readThermoInputFile": "read_thermo_input_file",
    "saveInputFile": "save_input_file",
    "getInput": "get_input",
    # rmgpy.thermo.thermoengine
    "processThermoData": "process_thermo_data",
    "generateThermoData": "generate_thermo_data",
    # rmgpy.thermo.wilhoit
    "Wilhoit_to_NASA": "wilhoit_to_nasa",
    "Wilhoit_to_NASA_TintOpt": "wilhoit_to_nasa_t_int_opt",
    "Wilhoit_to_NASA_TintOpt_objFun": "wilhoit_to_nasa_t_int_opt_obj_fun",
    "Wilhoit_to_NASA_TintOpt_objFun_NW": "wilhoit_to_nasa_t_int_opt_obj_fun_nw",
    "Wilhoit_to_NASA_TintOpt_objFun_W": "wilhoit_to_nasa_t_int_opt_obj_fun_w",
    # rmgpy.quantity
    "conversionFactorsFromSItoCmMolS": "conversion_factors_from_si_to_cm_mol_s",
    # rmgpy.util
    "makeOutputSubdirectory": "make_output_subdirectory",
    # rmgpy.chemkin
    "readThermoEntry": "read_thermo_entry",
    "readKineticsEntry": "read_kinetics_entry",
    "readReactionComments": "read_reaction_comments",
    "loadSpeciesDictionary": "load_species_dictionary",
    "loadTransportFile": "load_transport_file",
    "loadChemkinFile": "load_chemkin_file",
    "readSpeciesBlock": "read_species_block",
    "readThermoBlock": "read_thermo_block",
    "readReactionsBlock": "read_reactions_block",
    "getSpeciesIdentifier": "get_species_identifier",
    "writeThermoEntry": "write_thermo_entry",
    "writeReactionString": "write_reaction_string",
    "writeKineticsEntry": "write_kinetics_entry",
    "markDuplicateReaction": "mark_duplicate_reaction",
    "markDuplicateReactions": "mark_duplicate_reactions",
    "saveSpeciesDictionary": "save_species_dictionary",
    "saveTransportFile": "save_transport_file",
    "saveChemkinFile": "save_chemkin_file",
    "saveChemkinSurfaceFile": "save_chemkin_surface_file",
    "saveChemkin": "save_chemkin",
    "saveChemkinFiles": "save_chemkin_files",
    "writeElementsSection": "write_elements_sections",
    # rmgpy.constraints
    "failsSpeciesConstraints": "fails_species_constraints",
    # rmgpy.tools.canteraModel
    "generateCanteraConditions": "generate_cantera_conditions",
    "getRMGSpeciesFromUserSpecies": "get_rmg_species_from_user_species",
    "findIgnitionDelay": "find_ignition_delay",
    "checkNearlyEqual": "check_nearly_equal",
    "checkEquivalentCanteraSpecies": "check_equivalent_cantera_species",
    "checkEquivalentCanteraReaction": "check_equivalent_cantera_reaction",
    # rmgpy.tools.diff_models
    "compareModelKinetics": "compare_model_kinetics",
    "compareModelSpecies": "compare_model_species",
    "compareModelReactions": "compare_model_reactions",
    "saveCompareHTML": "save_compare_html",
    "enthalpyDiff": "enthalpy_diff",
    "kineticsDiff": "kinetics_diff",
    "identicalThermo": "identical_thermo",
    "identicalKinetics": "identical_kinetics",
    "parseCommandLineArguments": "parse_command_line_arguments",
    # rmgpy.tools.fluxdiagram
    "maximumNodeCount": "max_node_count",
    "maximumEdgeCount": "max_edge_count",
    "concentrationTolerance": "concentration_tol",
    "speciesRateTolerance": "species_rate_tol",
    "maximumNodePenWidth": "max_node_pen_width",
    "maximumEdgePenWidth": "max_edge_pen_width",
    "centralReactionCount": "central_reaction_count",
    "initialTime": "initial_time",
    "timeStep": "time_step",
    "absoluteTolerance": "abs_tol",
    "relativeTolerance": "rel_tol",
    "framesPerSecond": "video_fps",
    "initialPadding": "initial_padding",
    "finalPadding": "final_padding",
    "generateFluxDiagram": "generate_flux_diagram",
    "addAdjacentNodes": "add_adjacent_nodes",
    "loadChemkinOutput": "load_chemkin_output",
    "createFluxDiagram": "create_flux_diagram",
    # rmgpy.tools.isotopes
    "generate_RMG_model": "generate_rmg_model",
    # rmgpy.tools.loader
    "loadRMGJob": "load_rmg_job",
    "loadRMGPyJob": "load_rmg_py_job",
    # rmgpy.tools.plot
    "parseCSVData": "parse_csv_data",
    "findNearest": "find_nearest",
    "linearlyInterpolatePoint": "linearly_interpolate_point",
}

GLOBALS2 = {
    # rmgpy.data.base
    "removeCommentFromLine": "remove_comment_from_line",
    "splitLineAndComment": "split_line_and_comment",
    # rmgpy.data.rmg
    "getDB": "get_db",
    # rmgpy.data.solvation
    "generateOldLibraryEntry": "generate_old_library_entry",
    "processOldLibraryEntry": "process_old_library_entry",
    # rmgpy.data.statmechfit
    "hoFreqLowerBound": "ho_freq_lower_bound",
    "hoFreqUpperBound": "ho_freq_upper_bound",
    "hrFreqLowerBound": "hr_freq_lower_bound",
    "hrFreqUpperBound": "hr_freq_upper_bound",
    "hrBarrLowerBound": "hr_barr_lower_bound",
    "hrBarrUpperBound": "hr_barr_upper_bound",
    "maxIter": "max_iter",
    "fitStatmechToHeatCapacity": "fit_statmech_to_heat_capacity",
    "fitStatmechDirect": "fit_statmech_direct",
    "fitStatmechPseudoRotors": "fit_statmech_pseudo_rotors",
    "fitStatmechPseudo": "fit_statmech_pseudo",
    "harmonicOscillator_heatCapacity": "harmonic_oscillator_heat_capacity",
    "harmonicOscillator_d_heatCapacity_d_freq": "harmonic_oscillator_d_heat_capacity_d_freq",
    "hinderedRotor_heatCapacity": "hindered_rotor_heat_capacity",
    "hinderedRotor_d_heatCapacity_d_freq": "hindered_rotor_d_heat_capacity_d_freq",
    "hinderedRotor_d_heatCapacity_d_barr": "hindered_rotor_d_heat_capacity_d_barr",
    # rmgpy.data.thermo
    "addThermoData": "add_thermo_data",
    "removeThermoData": "remove_thermo_data",
    "averageThermoData": "average_thermo_data",
    "commonAtoms": "common_atoms",
    "combineCycles": "combine_cycles",
    "isAromaticRing": "is_aromatic_ring",
    "isBicyclic": "is_bicyclic",
    "findAromaticBondsFromSubMolecule": "find_aromatic_bonds_from_sub_molecule",
    "convertRingToSubMolecule": "convert_ring_to_sub_molecule",
    "combineTwoRingsIntoSubMolecule": "combine_two_rings_into_sub_molecule",
    "getCopyForOneRing": "get_copy_for_one_ring",
    "getCopyFromTwoRingsWithCommonAtoms": "get_copy_from_two_rings_with_common_atoms",
    "isRingPartialMatched": "is_ring_partial_matched",
    "bicyclicDecompositionForPolyring": "bicyclic_decomposition_for_polyring",
    "splitBicyclicIntoSingleRings": "split_bicyclic_into_single_rings",
    "saturateRingBonds": "saturate_ring_bonds",
    # rmgpy.data.kinetics.family
    "_makeRule": "_make_rule",
    "_spawnTreeProcess": "_spawn_tree_process",
    "_childMakeTreeNodes": "_child_make_tree_nodes",
    # rmgpy.molecule.adjlist
    "fromOldAdjacencyList": "from_old_adjacency_list",
    "getOldElectronState": "get_old_electron_state",
    "toOldAdjacencyList": "to_old_adjacency_list",
    # rmgpy.molecule.atomtype
    "atomTypes": "ATOMTYPES",
    # rmgpy.molecule.converter
    "debugRDKitMol": "debug_rdkit_mol",
    # rmgpy.molecule.graph
    "_getEdgeVertex1": "_get_edge_vertex1",
    "_getEdgeVertex2": "_get_edge_vertex2",
    # rmgpy.molecule.inchi
    "_parse_H_layer": "_parse_h_layer",
    "_parse_E_layer": "_parse_e_layer",
    "_parse_N_layer": "_parse_n_layer",
    "_create_U_layer": "_create_u_layer",
    "_create_P_layer": "_create_p_layer",
    # rmgpy.qm.main
    "_write_QMfiles_star": "_write_qm_files_star",
    "_write_QMfiles": "_write_qm_files",
    # rmgpy.qm.molecule
    "loadThermoDataFile": "load_thermo_data_file",
    # rmgpy.qm.qmdata
    "parseCCLibData": "parse_cclib_data",
    # rmgpy.qm.symmetry'
    "makePointGroupDictionary": "make_point_group_dictionary",
    "pointGroupDictionary": "point_group_dictionary",
    # rmgpy.rmg.model
    "generateReactionKey": "generate_reaction_key",
    "generateReactionId": "generate_reaction_id",
    "getFamilyLibraryObject": "get_family_library_object",
    "getKey": "get_key",
    "areIdenticalSpeciesReferences": "are_identical_species_references",
    # rmgpy.rmg.output
    "saveOutputHTML": "save_output_html",
    "saveDiffHTML": "save_diff_html",
    "saveOutput": "save_output",
    # rmgpy.statmech.conformer
    "getDensityOfStatesForst": "get_density_of_states_forst",
    # rmgpy.statmech.schrodinger
    "unitDegeneracy": "unit_degeneracy",
    "convolveBS": "convolve_bs",
    "convolveBSSR": "convolve_bssr",
    # rmgpy.util
    "makeOutputSubdirectory": "make_output_subdirectory",
    # rmgpy.yml
    "convertChemkin2yml": "convert_chemkin_to_yml",
    "writeyml": "write_yml",
    "getMechDict": "get_mech_dict",
    "getRadicals": "get_radicals",
    "obj2dict": "obj_to_dict",
    # rmgpy.chemkin
    "__chemkin_reaction_count": "_chemkin_reaction_count",
    "Ffloat": "fortran_float",
    "_readKineticsReaction": "_read_kinetics_reaction",
    "_readKineticsLine": "_read_kinetics_line",
    "_removeLineBreaks": "_remove_line_breaks",
    "saveJavaKineticsLibrary": "save_java_kinetics_library",
    # rmgpy.tools.extractInfoFromckcsv
    "getROPFromCKCSV": "get_rop_from_ckcsv",
    "getConcentrationDictFromCKCSV": "get_concentration_dict_from_ckcsv",
    "getFluxGraphEdgesDict": "get_flux_graph_edges_dict",
    "getROPFlux": "get_rop_flux",
    # rmgpy.tools.observablesRegression
    "curvesSimilar": "curves_similar",
    # rmgpy.tools.regression
    "parseArguments": "parse_command_line_arguments",
}

# Class attributes
ATTRIBUTES1 = {
    # Arkane:
    "angleUnits": "angle_units",
    "energyUnits": "energy_units",
    "cosineRotor": "cosine_rotor",
    "fourierRotor": "fourier_rotor",
    "rotorIndex": "rotor_index",
    # rmgpy.data.base
    "shortDesc": "short_desc",
    "longDesc": "long_desc",
    "referenceType": "reference_type",
    "nodalDistance": "nodal_distance",
    # rmgpy.data.rmg
    "forbiddenStructures": "forbidden_structures",
    # rmgpy.data.thermo.ThermoDatabase
    "libraryOrder": "library_order",
    "deltaAtomicAdsorptionEnergy": "delta_atomic_adsorption_energy",
    "genericNodes": "generic_nodes",
    # rmgpy.data.kinetics.database.KineticsDatabase
    "recommendedFamilies": "recommended_families",
    # rmgpy.data.kinetics.depository.DepositoryReaction
    "specificCollider": "specific_collider",
    "transitionState": "transition_state",
    # rmgpy.data.kinetics.family.KineticsFamily
    "forwardTemplate": "forward_template",
    "forwardRecipe": "forward_recipe",
    "reverseTemplate": "reverse_template",
    "reverseRecipe": "reverse_recipe",
    "ownReverse": "own_reverse",
    "boundaryAtoms": "boundary_atoms",
    "treeDistance": "tree_distance",
    "reverseMap": "reverse_map",
    "reactantNum": "reactant_num",
    "productNum": "product_num",
    "autoGenerated": "auto_generated",
    # rmgpy.kinetics.diffusionLimited
    "solventData": "solvent_data",
    # rmgpy.molecule.atomtype
    "incrementBond": "increment_bond",
    "decrementBond": "decrement_bond",
    "formBond": "form_bond",
    "breakBond": "break_bond",
    "incrementRadical": "increment_radical",
    "decrementRadical": "decrement_radical",
    "incrementLonePair": "increment_lone_pair",
    "decrementLonePair": "decrement_lone_pair",
    "allDouble": "all_double",
    "rDouble": "r_double",
    "oDouble": "o_double",
    "sDouble": "s_double",
    "lonePairs": "lone_pairs",
    # rmgpy.molecule.element
    "chemkinName": "chemkin_name",
    "covRadius": "cov_radius",
    "elementList": "element_list",
    # rmgpy.molecule.graph
    "sortingLabel": "sorting_label",
    # rmgpy.molecule.molecule
    "radicalElectrons": "radical_electrons",
    "atomType": "atomtype",
    "symmetryNumber": "symmetry_number",
    # rmgpy.pdep.configuration
    "Elist": "e_list",
    "densStates": "dens_states",
    "sumStates": "sum_states",
    "activeJRotor": "active_j_rotor",
    "activeKRotor": "active_k_rotor",
    # rmgpy.pdep.network
    "pathReactions": "path_reactions",
    "bathGas": "bath_gas",
    "netReactions": "net_reactions",
    "Jlist": "j_list",
    "Nisom": "n_isom",
    "Nreac": "n_reac",
    "Nprod": "n_prod",
    "Ngrains": "n_grains",
    "NJ": "n_j",
    "grainSize": "grain_size",
    "grainCount": "grain_count",
    # rmgpy.qm.molecule
    "uniqueID": "unique_id",
    "uniqueIDlong": "unique_id_long",
    "outputFilePath": "output_file_path",
    "inputFilePath": "input_file_path",
    "scriptAttempts": "script_attempts",
    "maxAttempts": "max_attempts",
    "qmData": "qm_data",
    # rmgpy.qm.symmetry
    "pointGroup": "point_group",
    "attemptNumber": "attempt_number",
    "pointGroupFound": "point_group_found",
    # rmgpy.rmg.main
    "inputFile": "input_file",
    "outputDirectory": "output_directory",
    "modelSettingsList": "model_settings_list",
    "simulatorSettingsList": "simulator_settings_list",
    "databaseDirectory": "database_directory",
    "thermoLibraries": "thermo_libraries",
    "transportLibraries": "transport_libraries",
    "reactionLibraries": "reaction_libraries",
    "statmechLibraries": "statmech_libraries",
    "seedMechanisms": "seed_mechanisms",
    "kineticsFamilies": "kinetics_families",
    "kineticsDepositories": "kinetics_depositories",
    "kineticsEstimator": "kinetics_estimator",
    "diffusionLimiter": "diffusion_limiter",
    "bindingEnergies": "binding_energies",
    "reactionModel": "reaction_model",
    "reactionSystems": "reaction_systems",
    "balanceSpecies": "balance_species",
    "filterReactions": "filter_reactions",
    "unimolecularReact": "unimolecular_react",
    "bimolecularReact": "bimolecular_react",
    "trimolecularReact": "trimolecular_react",
    "generateOutputHTML": "generate_output_html",
    "generatePlots": "generate_plots",
    "saveSimulationProfiles": "save_simulation_profiles",
    "verboseComments": "verbose_comments",
    "saveEdgeSpecies": "save_edge_species",
    "keepIrreversible": "keep_irreversible",
    "trimolecularProductReversible": "trimolecular_product_reversible",
    "pressureDependence": "pressure_dependence",
    "quantumMechanics": "quantum_mechanics",
    "speciesConstraints": "species_constraints",
    "wallTime": "walltime",
    "initialSpecies": "initial_species",
    "initializationTime": "initialization_time",
    "kineticsdatastore": "kinetics_datastore",
    "coreSeedPath": "core_seed_path",
    "edgeSeedPath": "edge_seed_path",
    "filtersPath": "filters_path",
    "speciesMapPath": "species_map_path",
    "generateSeedEachIteration": "generate_seed_each_iteration",
    "saveSeedToDatabase": "save_seed_to_database",
    "thermoCentralDatabase": "thermo_central_database",
    "execTime": "exec_time",
    "reactionSystem": "reaction_system",
    "conditionList": "condition_list",
    "scaledConditionList": "scaled_condition_list",
    "randState": "rand_state",
    # rmgpy.rmg.model
    "networkDict": "network_dict",
    "networkList": "network_list",
    "networkCount": "network_count",
    "speciesDict": "species_dict",
    "reactionDict": "reaction_dict",
    "speciesCache": "species_cache",
    "speciesCounter": "species_counter",
    "reactionCounter": "reaction_counter",
    "newSpeciesList": "new_species_list",
    "newReactionList": "new_reaction_list",
    "outputSpeciesList": "output_species_list",
    "outputReactionList": "output_reaction_list",
    "indexSpeciesDict": "index_species_dict",
    "iterationNum": "iteration_num",
    "toleranceThermoKeepSpeciesInEdge": "thermo_tol_keep_spc_in_edge",
    "minCoreSizeForPrune": "min_core_size_for_prune",
    "maximumEdgeSpecies": "maximum_edge_species",
    "newSurfaceSpcsAdd": "new_surface_spcs_add",
    "newSurfaceRxnsAdd": "new_surface_rxns_add",
    "newSurfaceSpcsLoss": "new_surface_spcs_loss",
    "newSurfaceRxnsLoss": "new_surface_rxns_loss",
    "solventName": "solvent_name",
    # rmgpy.rmg.settings
    "fluxToleranceKeepInEdge": "tol_keep_in_edge",
    "fluxToleranceMoveToCore": "tol_move_to_core",
    "toleranceMoveEdgeReactionToCore": "tol_move_edge_rxn_to_core",
    "fluxToleranceInterrupt": "tol_interrupt_simulation",
    "minSpeciesExistIterationsForPrune": "min_species_exist_iterations_for_prune",
    "filterThreshold": "filter_threshold",
    "ignoreOverallFluxCriterion": "ignore_overall_flux_criterion",
    "toleranceMoveEdgeReactionToSurface": "tol_move_edge_rxn_to_surface",
    "toleranceMoveSurfaceSpeciesToCore": "tol_move_surface_spc_to_core",
    "toleranceMoveSurfaceReactionToCore": "tol_move_surface_rxn_to_core",
    "terminateAtMaxObjects": "terminate_at_max_objects",
    "dynamicsTimeScale": "dynamics_time_scale",
    "toleranceBranchReactionToCore": "tol_branch_rxn_to_core",
    "branchingIndex": "branching_index",
    "branchingRatioMax": "branching_ratio_max",
    "toleranceMoveEdgeReactionToSurfaceInterrupt": "tol_move_edge_rxn_to_surface_interrupt",
    "toleranceMoveEdgeReactionToCoreInterrupt": "tol_move_edge_rxn_to_core_interrupt",
    "maxNumSpecies": "max_num_species",
    "maxNumObjsPerIter": "max_num_objects_per_iter",
    # rmgpy.statmech.conformer
    "spinMultiplicity": "spin_multiplicity",
    "opticalIsomers": "optical_isomers",
    # rmgpy.tools.canteraModel
    "reactorType": "reactor_type",
    "reactionTime": "reaction_time",
    "molFrac": "mol_frac",
    "speciesList": "species_list",
    "reactionList": "reaction_list",
    "reactionMap": "reaction_map",
    # rmgpy.species
    "transportData": "transport_data",
    "molecularWeight": "molecular_weight",
    "energyTransferModel": "energy_transfer_model",
    "isSolvent": "is_solvent",
    "creationIteration": "creation_iteration",
    "explicitlyAllowed": "explicitly_allowed",
}

ATTRIBUTES2 = {
    # rmgpy.molecule.vf2
    "initialMapping": "initial_mapping",
    "findAll": "find_all",
    "isMatch": "is_match",
    "mappingList": "mapping_list",
    # rmgpy.molecule.molecule
    "InChI": "inchi",
    "SMILES": "smiles",
    # rmgpy.rmg.pdep
    "collFreq": "coll_freq",
    # rmgpy.solver.base
    "numCoreSpecies": "num_core_species",
    "numCoreReactions": "num_core_reactions",
    "numEdgeSpecies": "num_edge_species",
    "numEdgeReactions": "num_edge_reactions",
    "numPdepNetworks": "num_pdep_networks",
    "speciesIndex": "species_index",
    "reactionIndex": "reaction_index",
    "reactantIndices": "reactant_indices",
    "productIndices": "product_indices",
    "networkIndices": "network_indices",
    "networkLeakCoefficients": "network_leak_coefficients",
    "jacobianMatrix": "jacobian_matrix",
    "coreSpeciesConcentrations": "core_species_concentrations",
    "coreSpeciesRates": "core_species_rates",
    "coreReactionRates": "core_reaction_rates",
    "coreSpeciesProductionRates": "core_species_production_rates",
    "coreSpeciesConsumptionRates": "core_species_consumption_rates",
    "edgeSpeciesRates": "edge_species_rates",
    "edgeReactionRates": "edge_reaction_rates",
    "networkLeakRates": "network_leak_rates",
    "surfaceSpeciesIndices": "surface_species_indices",
    "surfaceReactionIndices": "surface_reaction_indices",
    "validLayeringIndices": "valid_layering_indices",
    "maxEdgeSpeciesRateRatios": "max_edge_species_rate_ratios",
    "maxNetworkLeakRateRatios": "max_network_leak_rate_ratios",
    "prunableSpecies": "prunable_species",
    "prunableNetworks": "prunable_networks",
    "prunableSpeciesIndices": "prunable_species_indices",
    "prunableNetworkIndices": "prunable_network_indices",
    "sensitivityCoefficients": "sensitivity_coefficients",
    "sensitiveSpecies": "sensitive_species",
    "sensitivityThreshold": "sensitivity_threshold",
    "unimolecularThreshold": "unimolecular_threshold",
    "bimolecularThreshold": "bimolecular_threshold",
    "trimolecularThreshold": "trimolecular_threshold",
    # rmgpy.solver.simple
    "constantVolume": "constant_volume",
    "initialMoleFractions": "initial_mole_fractions",
    "pdepColliderKinetics": "pdep_collider_kinetics",
    "colliderEfficiencies": "collider_efficiencies",
    "pdepColliderReactionIndices": "pdep_collision_reaction_indices",
    "pdepSpecificColliderKinetics": "pdep_specific_collider_kinetics",
    "specificColliderSpecies": "specific_collider_species",
    "pdepSpecificColliderReactionIndices": "pdep_specific_collider_reaction_indices",
    "sensConditions": "sens_conditions",
    "nSims": "n_sims",
    # rmgpy.solver.liquid
    "constSPCNames": "const_spc_names",
    "constSPCIndices": "const_spc_indices",
    "initialConcentrations": "initial_concentrations",
    # rmgpy.solver.surface
    "initialP": "P_initial",
    "initialGasMoleFractions": "initial_gas_mole_fractions",
    "initialSurfaceCoverages": "initial_surface_coverages",
    "surfaceVolumeRatio": "surface_volume_ratio",
    "surfaceSiteDensity": "surface_site_density",
    "reactionsOnSurface": "reactions_on_surface",
    "speciesOnSurface": "species_on_surface",
    # rmgpy.quantity
    "uncertaintyType": "uncertainty_type",
    "commonUnits": "common_units",
    "extraDimensionality": "extra_dimensionality",
    # rmgpy.statmech.ndTorsions
    "calcPath": "calc_path",
    "isLinear": "is_linear",
    "isTS": "is_ts",
    # rmgpy.tools.observablesRegression
    "oldDir": "old_dir",
    "newDir": "new_dir",
    "exptData": "expt_data",
    "oldSim": "old_sim",
    "newSim": "new_sim",
    # rmgpy.tools.plot
    "xVar": "x_var",
    "yVar": "y_var",
    "csvFile": "csv_file",
    "numSpecies": "num_species",
    "numReactions": "num_reactions",
    # rmgpy.tools.uncertainty
    "speciesSourcesDict": "species_sources_dict",
    "reactionSourcesDict": "reaction_sources_dict",
    "allThermoSources": "all_thermo_sources",
    "allKineticSources": "all_kinetic_sources",
    "thermoInputUncertainties": "thermo_input_uncertainties",
    "kineticInputUncertainties": "kinetic_input_uncertainties",
    "extraSpecies": "extra_species",
    # rmgpy.species
    "_molecularWeight": "_molecular_weight",
}

# Class methods
METHODS1 = {
    # Arkane:
    "loadInputFile": "load_input_file",
    "generateTemperatureList": "generate_T_list",
    "generatePressureList": "generate_P_list",
    "getNumberOfAtoms": "get_number_of_atoms",
    "loadForceConstantMatrix": "load_force_constant_matrix",
    "loadGeometry": "load_geometry",
    "loadConformer": "load_conformer",
    "loadEnergy": "load_energy",
    "loadZeroPointEnergy": "load_zero_point_energy",
    "loadScanEnergies": "load_scan_energies",
    "loadNegativeFrequency": "load_negative_frequency",
    "loadNecessaryDatabases": "load_necessary_databases",
    "getLibraries": "get_libraries",
    "visit_Call": "visit_call",
    "visit_List": "visit_list",
    "visit_Tuple": "visit_tuple",
    "visit_Dict": "visit_dict",
    "visit_Str": "visit_str",
    "visit_Num": "visit_num",
    "fitInterpolationModels": "fit_interpolation_models",
    "projectRotors": "project_rotors",
    # rmgpy.__init__
    "getPath": "get_path",
    # rmgpy.data.base
    "getAllDescendants": "get_all_descendants",
    "getEntriesToSave": "get_entries_to_save",
    "getSpecies": "get_species",
    "saveDictionary": "save_dictionary",
    "matchNodeToNode": "match_node_to_node",
    "matchNodeToChild": "match_node_to_child",
    "matchNodeToStructure": "match_node_to_structure",
    "descendTree": "descend_tree",
    "areSiblings": "are_siblings",
    "removeGroup": "remove_group",
    "matchToStructure": "match_to_structure",
    "matchLogicOr": "match_logic_or",
    "getPossibleStructures": "get_possible_structures",
    "isMoleculeForbidden": "is_molecule_forbidden",
    "loadEntry": "load_entry",
    "saveEntry": "save_entry",
    # rmgpy.data.rmg
    "loadThermo": "load_thermo",
    "loadTransport": "load_transport",
    "loadForbiddenStructures": "load_forbidden_structures",
    "loadKinetics": "load_kinetics",
    "loadSolvation": "load_solvation",
    "loadStatmech": "load_statmech",
    # rmgpy.data.solvation
    "getHAbsCorrection": "get_h_abs_correction",
    "getSolventViscosity": "get_solvent_viscosity",
    "getStokesDiffusivity": "get_stokes_diffusivity",
    "setMcGowanVolume": "set_mcgowan_volume",
    "getSolventData": "get_solvent_data",
    "getSolventStructure": "get_solvent_structure",
    "loadGroups": "load_groups",
    "saveLibraries": "save_libraries",
    "saveGroups": "save_groups",
    "getSoluteData": "get_solute_data",
    "getAllSoluteData": "get_all_solute_data",
    "getSoluteDataFromLibrary": "get_solute_data_from_library",
    "getSoluteDataFromGroups": "get_solute_data_from_groups",
    "transformLonePairs": "transform_lone_pairs",
    "removeHBonding": "remove_h_bonding",
    "estimateSoluteViaGroupAdditivity": "estimate_solute_via_group_additivity",
    "calcH": "calc_h",
    "calcG": "calc_g",
    "calcS": "calc_s",
    "getSolvationCorrection": "get_solvation_correction",
    "checkSolventinInitialSpecies": "check_solvent_in_initial_species",
    # rmgpy.data.statmech
    "getFrequencyGroups": "get_frequency_groups",
    "getStatmechData": "get_statmech_data",
    "loadDepository": "load_depository",
    "loadLibraries": "load_libraries",
    "saveDepository": "save_depository",
    "getStatmechDataFromDepository": "get_statmech_data_from_depository",
    "getStatmechDataFromLibrary": "get_statmech_data_from_library",
    "getStatmechDataFromGroups": "get_statmech_data_from_groups",
    "generateFrequencies": "generate_frequencies",
    # rmgpy.data.thermo
    "pruneHeteroatoms": "prune_heteroatoms",
    "recordPolycyclicGenericNodes": "record_polycyclic_generic_nodes",
    "recordRingGenericNodes": "record_ring_generic_nodes",
    "getThermoData": "get_thermo_data",
    "setDeltaAtomicAdsorptionEnergies": "set_delta_atomic_adsorption_energies",
    "correctBindingEnergy": "correct_binding_energy",
    "getThermoDataForSurfaceSpecies": "get_thermo_data_for_surface_species",
    "getThermoDataFromLibraries": "get_thermo_data_from_libraries",
    "getAllThermoData": "get_all_thermo_data",
    "getThermoDataFromDepository": "get_thermo_data_from_depository",
    "getThermoDataFromLibrary": "get_thermo_data_from_library",
    "getThermoDataFromGroups": "get_thermo_data_from_groups",
    "prioritizeThermo": "prioritize_thermo",
    "estimatRadicalThermoViaHBI": "estimate_radical_thermo_via_hbi",
    "estimateThermoViaGroupAdditivity": "estimate_thermo_via_group_additivity",
    "computeGroupAdditivityThermo": "compute_group_additivity_thermo",
    "getBicyclicCorrectionThermoDataFromHeuristic": "get_bicyclic_correction_thermo_data_from_heuristic",
    "getRingGroupsFromComments": "get_ring_groups_from_comments",
    "extractSourceFromComments": "extract_source_from_comments",
    # rmgpy.data.transport
    "getTransportProperties": "get_transport_properties",
    "getAllTransportProperties": "get_all_transport_properties",
    "getTransportPropertiesFromLibrary": "get_transport_properties_from_library",
    "getTransportPropertiesViaGroupEstimates": "get_transport_properties_via_group_estimates",
    "estimateCriticalPropertiesViaGroupAdditivity": "estimate_critical_properties_via_group_additivity",
    "getTransportPropertiesViaLennardJonesParameters": "get_transport_properties_via_lennard_jones_parameters",
    # rmgpy.data.kinetics.database
    "loadRecommendedFamiliesList": "load_recommended_families",
    "loadFamilies": "load_families",
    "saveRecommendedFamilies": "save_recommended_families",
    "saveFamilies": "save_families",
    "getForwardReactionForFamilyEntry": "get_forward_reaction_for_family_entry",
    "reconstructKineticsFromSource": "reconstruct_kinetics_from_source",
    # rmgpy.data.kinetics.family
    "applyForward": "apply_forward",
    "applyReverse": "apply_reverse",
    "loadTemplate": "load_template",
    "loadRecipe": "load_recipe",
    "loadForbidden": "load_forbidden",
    "saveTrainingReactions": "save_training_reactions",
    "generateProductTemplate": "generate_product_template",
    "addKineticsRulesFromTrainingSet": "add_rules_from_training",
    "getRootTemplate": "get_root_template",
    "fillKineticsRulesByAveragingUp": "fill_rules_by_averaging_up",
    "applyRecipe": "apply_recipe",
    "generateReactions": "generate_reactions",
    "addReverseAttribute": "add_reverse_attribute",
    "calculateDegeneracy": "calculate_degeneracy",
    "getReactionPairs": "get_reaction_pairs",
    "getReactionTemplate": "get_reaction_template",
    "getKineticsForTemplate": "get_kinetics_for_template",
    "getKineticsFromDepository": "get_kinetics_from_depository",
    "getKinetics": "get_kinetics",
    "estimateKineticsUsingGroupAdditivity": "estimate_kinetics_using_group_additivity",
    "estimateKineticsUsingRateRules": "estimate_kinetics_using_rate_rules",
    "getReactionTemplateLabels": "get_reaction_template_labels",
    "retrieveTemplate": "retrieve_template",
    "getLabeledReactantsAndProducts": "get_labeled_reactants_and_products",
    "addAtomLabelsForReaction": "add_atom_labels_for_reaction",
    "getTrainingDepository": "get_training_depository",
    "addEntry": "add_entry",
    # rmgpy.data.kinetics.library
    "getLibraryReactions": "get_library_reactions",
    "markValidDuplicates": "mark_valid_duplicates",
    "checkForDuplicates": "check_for_duplicates",
    "convertDuplicatesToMulti": "convert_duplicates_to_multi",
    # rmgpy.data.kinetics.rules
    "getEntries": "get_entries",
    "hasRule": "has_rule",
    "getRule": "get_rule",
    "getAllRules": "get_all_rules",
    "fillRulesByAveragingUp": "fill_rules_by_averaging_up",
    "estimateKinetics": "estimate_kinetics",
    # rmgpy.kinetics.model
    "isPressureDependent": "is_pressure_dependent",
    "isTemperatureValid": "is_temperature_valid",
    "getRateCoefficient": "get_rate_coefficient",
    "toHTML": "to_html",
    "isSimilarTo": "is_similar_to",
    "isIdenticalTo": "is_identical_to",
    "getCanteraEfficiencies": "get_cantera_efficiencies",
    "setCanteraKinetics": "set_cantera_kinetics",
    "isPressureValid": "is_pressure_valid",
    "getEffectivePressure": "get_effective_pressure",
    "getEffectiveColliderEfficiencies": "get_effective_collider_efficiencies",
    "calculateTunnelingFactor": "calculate_tunneling_factor",
    "calculateTunnelingFunction": "calculate_tunneling_function",
    # rmgpy.kinetics.arrhenius
    "changeT0": "change_t0",
    "fitToData": "fit_to_data",
    "changeRate": "change_rate",
    "toCanteraKinetics": "to_cantera_kinetics",
    "toArrheniusEP": "to_arrhenius_ep",
    "getActivationEnergy": "get_activation_energy",
    "toArrhenius": "to_arrhenius",
    "fitToReactions": "fit_to_reactions",
    "getAdjacentExpressions": "get_adjacent_expressions",
    # rmgpy.kinetics.chebyshev
    "getReducedTemperature": "get_reduced_temperature",
    "getReducedPressure": "get_reduced_pressure",
    # rmgpy.kinetics.diffusionLimited
    "getEffectiveRate": "get_effective_rate",
    "getDiffusionFactor": "get_diffusion_factor",
    "getDiffusionLimit": "get_diffusion_limit",
    # rmgpy.kinetics.surface
    "getStickingCoefficient": "get_sticking_coefficient",
    # rmgpy.kinetics.uncertainties
    "getExpectedLogUncertainty": "get_expected_log_uncertainty",
    # rmgpy.molecule.atomtype
    "setActions": "set_actions",
    "isSpecificCaseOf": "is_specific_case_of",
    # rmgpy.molecule.graph
    "resetConnectivityValues": "reset_connectivity_values",
    "getOtherVertex": "get_other_vertex",
    "addVertex": "add_vertex",
    "addEdge": "add_edge",
    "getAllEdges": "get_all_edges",
    "getEdges": "get_edges",
    "getEdge": "get_edge",
    "hasVertex": "has_vertex",
    "hasEdge": "has_edge",
    "removeVertex": "remove_vertex",
    "removeEdge": "remove_edge",
    "copyAndMap": "copy_and_map",
    "updateConnectivityValues": "update_connectivity_values",
    "sortVertices": "sort_vertices",
    "isIsomorphic": "is_isomorphic",
    "findIsomorphism": "find_isomorphism",
    "isSubgraphIsomorphic": "is_subgraph_isomorphic",
    "findSubgraphIsomorphisms": "find_subgraph_isomorphisms",
    "isCyclic": "is_cyclic",
    "isVertexInCycle": "is_vertex_in_cycle",
    "isEdgeInCycle": "is_edge_in_cycle",
    "getAllCyclicVertices": "get_all_cyclic_vertices",
    "getAllPolycyclicVertices": "get_all_polycyclic_vertices",
    "getPolycyclicRings": "get_polycycles",
    "getMonocyclicRings": "get_monocycles",
    "getDisparateRings": "get_disparate_cycles",
    "getAllCycles": "get_all_cycles",
    "getAllCyclesOfSize": "get_all_cycles_of_size",
    "getAllSimpleCyclesOfSize": "get_all_simple_cycles_of_size",
    "getSmallestSetOfSmallestRings": "get_smallest_set_of_smallest_rings",
    "getRelevantCycles": "get_relevant_cycles",
    "getMaxCycleOverlap": "get_max_cycle_overlap",
    "getLargestRing": "get_largest_ring",
    "isMappingValid": "is_mapping_valid",
    # rmgpy.molecule.molecule
    "isHydrogen": "is_hydrogen",
    "isNonHydrogen": "is_non_hydrogen",
    "isCarbon": "is_carbon",
    "isNitrogen": "is_nitrogen",
    "isOxygen": "is_oxygen",
    "isFluorine": "is_fluorine",
    "isSurfaceSite": "is_surface_site",
    "isSilicon": "is_silicon",
    "isSulfur": "is_sulfur",
    "isChlorine": "is_chlorine",
    "isIodine": "is_iodine",
    "isNOS": "is_nos",
    "setLonePairs": "set_lone_pairs",
    "incrementLonePairs": "increment_lone_pairs",
    "decrementLonePairs": "decrement_lone_pairs",
    "updateCharge": "update_charge",
    "applyAction": "apply_action",
    "getBondOrdersForAtom": "get_total_bond_order",
    "getBDE": "get_bde",
    "getOrderStr": "get_order_str",
    "setOrderStr": "set_order_str",
    "getOrderNum": "get_order_num",
    "setOrderNum": "set_order_num",
    "isVanDerWaals": "is_van_der_waals",
    "isOrder": "is_order",
    "incrementOrder": "increment_order",
    "decrementOrder": "decrement_order",
    "addAtom": "add_atom",
    "addBond": "add_bond",
    "getBonds": "get_bonds",
    "getBond": "get_bond",
    "hasAtom": "has_atom",
    "hasBond": "has_bond",
    "containsSurfaceSite": "contains_surface_site",
    "removeAtom": "remove_atom",
    "removeBond": "remove_bond",
    "removeVanDerWaalsBonds": "remove_van_der_waals_bonds",
    "sortAtoms": "sort_atoms",
    "getFormula": "get_formula",
    "getMolecularWeight": "get_molecular_weight",
    "getRadicalCount": "get_radical_count",
    "getSingletCarbeneCount": "get_singlet_carbene_count",
    "getNumAtoms": "get_num_atoms",
    "deleteHydrogens": "delete_hydrogens",
    "connectTheDots": "connect_the_dots",
    "updateAtomTypes": "update_atomtypes",
    "updateMultiplicity": "update_multiplicity",
    "clearLabeledAtoms": "clear_labeled_atoms",
    "containsLabeledAtom": "contains_labeled_atom",
    "getLabeledAtoms": "get_all_labeled_atoms",
    "getLabeledAtom": "get_labeled_atoms",
    "isAtomInCycle": "is_atom_in_cycle",
    "isBondInCycle": "is_bond_in_cycle",
    "fromInChI": "from_inchi",
    "fromAugmentedInChI": "from_augmented_inchi",
    "fromSMILES": "from_smiles",
    "fromSMARTS": "from_smarts",
    "fromXYZ": "from_xyz",
    "toSingleBonds": "to_single_bonds",
    "toInChI": "to_inchi",
    "toAugmentedInChI": "to_augmented_inchi",
    "toInChIKey": "to_inchi_key",
    "toAugmentedInChIKey": "to_augmented_inchi_key",
    "toSMARTS": "to_smarts",
    "toSMILES": "to_smiles",
    "find_H_bonds": "find_h_bonds",
    "generate_H_bonded_structures": "generate_h_bonded_structures",
    "remove_H_bonds": "remove_h_bonds",
    "isLinear": "is_linear",
    "isAromatic": "is_aromatic",
    "isHeterocyclic": "is_heterocyclic",
    "countInternalRotors": "count_internal_rotors",
    "calculateCp0": "calculate_cp0",
    "calculateCpInf": "calculate_cpinf",
    "getSymmetryNumber": "get_symmetry_number",
    "calculateSymmetryNumber": "calculate_symmetry_number",
    "isRadical": "is_radical",
    "isArylRadical": "is_aryl_radical",
    "getURL": "get_url",
    "getRadicalAtoms": "get_radical_atoms",
    "updateLonePairs": "update_lone_pairs",
    "getNetCharge": "get_net_charge",
    "getChargeSpan": "get_charge_span",
    "toGroup": "to_group",
    "getAromaticRings": "get_aromatic_rings",
    "assignAtomIDs": "assign_atom_ids",
    "atomIDValid": "atom_ids_valid",
    "isIdentical": "is_identical",
    "getNthNeighbor": "get_nth_neighbor",
    # rmgpy.molecule.group
    "hasWildcards": "has_wildcards",
    "countBonds": "count_bonds",
    "makeSampleAtom": "make_sample_atom",
    "makeBond": "make_bond",
    "sortByConnectivity": "sort_by_connectivity",
    "clearRegDims": "clear_reg_dims",
    "getExtensions": "get_extensions",
    "specifyAtomExtensions": "specify_atom_extensions",
    "specifyRingExtensions": "specify_ring_extensions",
    "specifyUnpairedExtensions": "specify_unpaired_extensions",
    "specifyInternalNewBondExtensions": "specify_internal_new_bond_extensions",
    "specifyExternalNewBondExtensions": "specify_external_new_bond_extensions",
    "specifyBondExtensions": "specify_bond_extensions",
    "updateFingerprint": "update_fingerprint",
    "standardizeAtomType": "standardize_atomtype",
    "createAndConnectAtom": "create_and_connect_atom",
    "addExplicitLigands": "add_explicit_ligands",
    "standardizeGroup": "standardize_group",
    "addImplicitAtomsFromAtomType": "add_implicit_atoms_from_atomtype",
    "classifyBenzeneCarbons": "classify_benzene_carbons",
    "addImplicitBenzene": "add_implicit_benzene",
    "pickWildcards": "pick_wildcards",
    "makeSampleMolecule": "make_sample_molecule",
    "isBenzeneExplicit": "is_benzene_explicit",
    "mergeGroups": "merge_groups",
    "resetRingMembership": "reset_ring_membership",
    # rmgpy.pdep.configuration
    "isUnimolecular": "is_unimolecular",
    "isBimolecular": "is_bimolecular",
    "isTermolecular": "is_termolecular",
    "isTransitionState": "is_transition_state",
    "calculateCollisionFrequency": "calculate_collision_frequency",
    "calculateDensityOfStates": "calculate_density_of_states",
    "mapDensityOfStates": "map_density_of_states",
    "mapSumOfStates": "map_sum_of_states",
    # rmgpy.pdep.network
    "getAllSpecies": "get_all_species",
    "calculateRateCoefficients": "calculate_rate_coefficients",
    "setConditions": "set_conditions",
    "selectEnergyGrains": "select_energy_grains",
    "calculateDensitiesOfStates": "calculate_densities_of_states",
    "mapDensitiesOfStates": "map_densities_of_states",
    "calculateMicrocanonicalRates": "calculate_microcanonical_rates",
    "calculateEquilibriumRatios": "calculate_equilibrium_ratios",
    "calculateCollisionModel": "calculate_collision_model",
    "solveFullME": "solve_full_me",
    "solveReducedME": "solve_reduced_me",
    # rmgpy.rmg.main
    "loadInput": "load_input",
    "loadThermoInput": "load_thermo_input",
    "checkInput": "check_input",
    "checkLibraries": "check_libraries",
    "saveInput": "save_input",
    "loadDatabase": "load_database",
    "makeSeedMech": "make_seed_mech",
    "makeSpeciesLabelsIndependent": "make_species_labels_independent",
    "processToSpeciesNetworks": "process_to_species_networks",
    "processPdepNetworks": "process_pdep_networks",
    "processReactionsToSpecies": "process_reactions_to_species",
    "generateCanteraFiles": "generate_cantera_files",
    "initializeReactionThresholdAndReactFlags": "initialize_reaction_threshold_and_react_flags",
    "updateReactionThresholdAndReactFlags": "update_reaction_threshold_and_react_flags",
    "saveEverything": "save_everything",
    "getGitCommit": "get_git_commit",
    "logHeader": "log_header",
    "loadRMGJavaInput": "load_rmg_java_input",
    "readMeaningfulLineJava": "read_meaningful_line_java",
    "determine_procnum_from_RAM": "determine_procnum_from_ram",
    # rmgpy.rmg.model
    "checkForExistingSpecies": "check_for_existing_species",
    "makeNewSpecies": "make_new_species",
    "checkForExistingReaction": "check_for_existing_reaction",
    "makeNewReaction": "make_new_reaction",
    "makeNewPDepReaction": "make_new_pdep_reaction",
    "addNewSurfaceObjects": "add_new_surface_objects",
    "adjustSurface": "adjust_surface",
    "clearSurfaceAdjustments": "clear_surface_adjustments",
    "processNewReactions": "process_new_reactions",
    "applyThermoToSpecies": "apply_thermo_to_species",
    "generateThermo": "generate_thermo",
    "applyKineticsToReaction": "apply_kinetics_to_reaction",
    "generateKinetics": "generate_kinetics",
    "printEnlargeSummary": "log_enlarge_summary",
    "addSpeciesToCore": "add_species_to_core",
    "addSpeciesToEdge": "add_species_to_edge",
    "setThermodynamicFilteringParameters": "set_thermodynamic_filtering_parameters",
    "thermoFilterSpecies": "thermo_filter_species",
    "thermoFilterDown": "thermo_filter_down",
    "removeEmptyPdepNetworks": "remove_empty_pdep_networks",
    "removeSpeciesFromEdge": "remove_species_from_edge",
    "addReactionToCore": "add_reaction_to_core",
    "addReactionToEdge": "add_reaction_to_edge",
    "getModelSize": "get_model_size",
    "getLists": "get_species_reaction_lists",
    "getStoichiometryMatrix": "get_stoichiometry_matrix",
    "addSeedMechanismToCore": "add_seed_mechanism_to_core",
    "addReactionLibraryToEdge": "add_reaction_library_to_edge",
    "addReactionLibraryToOutput": "add_reaction_library_to_output",
    "addReactionToUnimolecularNetworks": "add_reaction_to_unimolecular_networks",
    "updateUnimolecularReactionNetworks": "update_unimolecular_reaction_networks",
    "markChemkinDuplicates": "mark_chemkin_duplicates",
    "registerReaction": "register_reaction",
    "searchRetrieveReactions": "search_retrieve_reactions",
    "initializeIndexSpeciesDict": "initialize_index_species_dict",
    # rmgpy.rmg.pdep
    "getLeakCoefficient": "get_leak_coefficient",
    "getMaximumLeakSpecies": "get_maximum_leak_species",
    "getLeakBranchingRatios": "get_leak_branching_ratios",
    "exploreIsomer": "explore_isomer",
    "addPathReaction": "add_path_reaction",
    "solve_SS_network": "solve_ss_network",
    "updateConfigurations": "update_configurations",
    # rmgpy.solver.base
    "initializeModel": "initialize_model",
    "getLayeringIndices": "get_layering_indices",
    "addReactionsToSurface": "add_reactions_to_surface",
    "logRates": "log_rates",
    "logConversions": "log_conversions",
    "computeRateDerivative": "compute_rate_derivative",
    # rmgpy.solver.simple
    "convertInitialKeysToSpeciesObjects": "convert_initial_keys_to_species_objects",
    # rmgpy.solver.liquid
    "get_constSPCIndices": "get_const_spc_indices",
    # rmgpy.statmech.conformer
    "getPartitionFunction": "get_partition_function",
    "getHeatCapacity": "get_heat_capacity",
    "getEnthalpy": "get_enthalpy",
    "getEntropy": "get_entropy",
    "getFreeEnergy": "get_free_energy",
    "getSumOfStates": "get_sum_of_states",
    "getDensityOfStates": "get_density_of_states",
    "getTotalMass": "get_total_mass",
    "getCenterOfMass": "get_center_of_mass",
    "getNumberDegreesOfFreedom": "get_number_degrees_of_freedom",
    "getMomentOfInertiaTensor": "get_moment_of_inertia_tensor",
    "getPrincipalMomentsOfInertia": "get_principal_moments_of_inertia",
    "getInternalReducedMomentOfInertia": "get_internal_reduced_moment_of_inertia",
    "getSymmetricTopRotors": "get_symmetric_top_rotors",
    "getActiveModes": "get_active_modes",
    # rmgpy.thermo.nasa
    "changeBaseEnthalpy": "change_base_enthalpy",
    "changeBaseEntropy": "change_base_entropy",
    "selectPolynomial": "select_polynomial",
    "toThermoData": "to_thermo_data",
    "toWilhoit": "to_wilhoit",
    # rmgpy.thermo.thermodata
    "toNASA": "to_nasa",
    # rmgpy.reaction
    "toLabeledStr": "to_labeled_str",
    "isIsomerization": "is_isomerization",
    "isAssociation": "is_association",
    "isDissociation": "is_dissociation",
    "isSurfaceReaction": "is_surface_reaction",
    "hasTemplate": "has_template",
    "matchesSpecies": "matches_species",
    "getEnthalpyOfReaction": "get_enthalpy_of_reaction",
    "getEntropyOfReaction": "get_entropy_of_reaction",
    "getFreeEnergyOfReaction": "get_free_energy_of_reaction",
    "getEquilibriumConstant": "get_equilibrium_constant",
    "getEnthalpiesOfReaction": "get_enthalpies_of_reaction",
    "getEntropiesOfReaction": "get_entropies_of_reaction",
    "getFreeEnergiesOfReaction": "get_free_energies_of_reaction",
    "getEquilibriumConstants": "get_equilibrium_constants",
    "getStoichiometricCoefficient": "get_stoichiometric_coefficient",
    "getSurfaceRateCoefficient": "get_surface_rate_coefficient",
    "fixDiffusionLimitedA": "fix_diffusion_limited_a_factor",
    "fixBarrierHeight": "fix_barrier_height",
    "reverseThisArrheniusRate": "reverse_arrhenius_rate",
    "generateReverseRateCoefficient": "generate_reverse_rate_coefficient",
    "calculateTSTRateCoefficient": "calculate_tst_rate_coefficient",
    "calculateTSTRateCoefficients": "calculate_tst_rate_coefficients",
    "canTST": "can_tst",
    "calculateMicrocanonicalRateCoefficient": "calculate_microcanonical_rate_coefficient",
    "isBalanced": "is_balanced",
    "generatePairs": "generate_pairs",
    "generate3dTS": "generate_3d_ts",
    # rmgpy.species
    "toChemkin": "to_chemkin",
    "toCantera": "to_cantera",
    "hasStatMech": "has_statmech",
    "hasThermo": "has_thermo",
    "getResonanceHybrid": "get_resonance_hybrid",
    "getAugmentedInChI": "get_augmented_inchi",
    "generateTransportData": "generate_transport_data",
    "getTransportData": "get_transport_data",
    "generateStatMech": "generate_statmech",
    "setE0WithThermo": "set_e0_with_thermo",
    "generateEnergyTransferModel": "generate_energy_transfer_model",
    # rmgpy.transport
    "getCollisionFrequency": "get_collision_frequency",
    # rmgpy.quantity
    "getConversionFactorToSI": "get_conversion_factor_to_si",
    "getConversionFactorFromSI": "get_conversion_factor_from_si",
    "getConversionFactorFromSItoCmMolS": "get_conversion_factor_from_si_to_cm_mol_s",
    "isUncertaintyAdditive": "is_uncertainty_additive",
    "isUncertaintyMultiplicative": "is_uncertainty_multiplicative",
    # rmgpy.tools.canteraModel
    "generateConditions": "generate_conditions",
    "loadModel": "load_model",
    "refreshModel": "refresh_model",
    "loadChemkinModel": "load_chemkin_model",
    "modifyReactionKinetics": "modify_reaction_kinetics",
    "modifySpeciesThermo": "modify_species_thermo",
    # rmgpy.tools.plot
    "comparePlot": "compare_plot",
    "uncertaintyPlot": "uncertainty_plot",
}

METHODS2 = {
    # rmgpy.data.base
    "loadOld": "load_old",
    "loadOldDictionary": "load_old_dictionary",
    "__loadTree": "_load_tree",
    "loadOldTree": "load_old_tree",
    "loadOldLibrary": "load_old_library",
    "parseOldLibrary": "parse_old_library",
    "saveOld": "save_old",
    "saveOldDictionary": "save_old_dictionary",
    "generateOldTree": "generate_old_tree",
    "saveOldTree": "save_old_tree",
    "saveOldLibrary": "save_old_library",
    "__hashLabels": "_hash_labels",
    # rmgpy.data.reference
    "toPrettyRepr": "to_pretty_repr",
    "getAuthorString": "get_author_string",
    # rmgpy.data.solvation
    "__addGroupSoluteData": "_add_group_solute_data",
    # rmgpy.data.statmech
    "__countMatchesToNode": "_count_matches_to_node",
    "__getNode": "_get_node",
    # rmgpy.data.thermo
    "copyData": "copy_data",
    "__addPolycyclicCorrectionThermoData": "_add_polycyclic_correction_thermo_data",
    "__addPolyRingCorrectionThermoDataFromHeuristic": "_add_poly_ring_correction_thermo_data_from_heuristic",
    "__addRingCorrectionThermoDataFromTree": "_add_ring_correction_thermo_data_from_tree",
    "__averageChildrenThermo": "_average_children_thermo",
    "__addGroupThermoData": "_add_group_thermo_data",
    "__removeGroupThermoData": "_remove_group_thermo_data",
    "satisfyRegistrationRequirements": "satisfy_registration_requirements",
    "registerInCentralThermoDB": "register_in_central_thermo_db",
    # rmgpy.data.transport
    "__addCriticalPointContribution": "_add_critical_point_contribution",
    # rmgpy.data.kinetics.depository
    "getSource": "get_source",
    # rmgpy.data.kinetics.family
    "addAction": "add_action",
    "getReverse": "get_reverse",
    "__apply": "_apply",
    "loadOldTemplate": "load_old_template",
    "saveOldTemplate": "save_old_template",
    "distributeTreeDistances": "distribute_tree_distances",
    "__generateProductStructures": "_generate_product_structures",
    "__createReaction": "_create_reaction",
    "__matchReactantToTemplate": "_match_reactant_to_template",
    "__generateReactions": "_generate_reactions",
    "__selectBestKinetics": "_select_best_kinetics",
    "_splitReactions": "_split_reactions",
    "evalExt": "eval_ext",
    "getExtensionEdge": "get_extension_edge",
    "extendNode": "extend_node",
    "generateTree": "generate_tree",
    "getRxnBatches": "get_rxn_batches",
    "pruneTree": "prune_tree",
    "makeTreeNodes": "make_tree_nodes",
    "_absorbProcess": "_absorb_process",
    "makeBMRulesFromTemplateRxnMap": "make_bm_rules_from_template_rxn_map",
    "crossValidate": "cross_validate",
    "crossValidateOld": "cross_validate_old",
    "simpleRegularization": "simple_regularization",
    "checkTree": "check_tree",
    "makeTree": "make_tree",
    "cleanTreeRules": "clean_tree_rules",
    "cleanTreeGroups": "clean_tree_groups",
    "cleanTree": "clean_tree",
    "saveGeneratedTree": "save_generated_tree",
    "getTrainingSet": "get_training_set",
    "getReactionMatches": "get_reaction_matches",
    "isEntryMatch": "is_entry_match",
    "rxnsMatchNode": "rxns_match_node",
    "retrieveOriginalEntry": "retrieve_original_entry",
    "getSourcesForTemplate": "get_sources_for_template",
    "getBackboneRoots": "get_backbone_roots",
    "getEndRoots": "get_end_roots",
    "getTopLevelGroups": "get_top_level_groups",
    # rmgpy.data.kinetics.groups
    "__multipleKineticsData": "_multiple_kinetics_data",
    "generateGroupAdditivityValues": "generate_group_additivity_values",
    # rmgpy.data.kinetics.library
    "__loadOldReactions": "_load_old_reactions",
    # rmgpy.data.kinetics.rules
    "__loadOldComments": "_load_old_comments",
    "__getAverageKinetics": "_get_average_kinetics",
    # rmgpy.molecule.draw
    "createNewSurface": "create_new_surface",
    "__findRingGroups": "_find_ring_groups",
    "__generateCoordinates": "_generate_coordinates",
    "__findCyclicBackbone": "_find_cyclic_backbone",
    "__findStrightChainBackbone": "_find_straight_chain_backbone",
    "__findStraightChainPaths": "_find_straight_chain_paths",
    "__generateRingSystemCoordinates": "_generate_ring_system_coordinates",
    "__generateStraightChainCoordinates": "_generate_straight_chain_coordinates",
    "__generateNeighborCoordinates": "_generate_neighbor_coordinates",
    "__generateFunctionalGroupCoordinates": "_generate_functional_group_coordinates",
    "__generateAtomLabels": "_generate_atom_labels",
    "__drawLine": "_draw_line",
    "__renderBond": "_render_bond",
    "__renderAtom": "_render_atom",
    "__make_single_bonds": "_make_single_bonds",
    "__replace_bonds": "_replace_bonds",
    # rmgpy.molecule.graph
    "__isChainInCycle": "_is_chain_in_cycle",
    "__exploreCyclesRecursively": "_explore_cycles_recursively",
    "_sortCyclicVertices": "sort_cyclic_vertices",
    # rmgpy.molecule.vf2
    "addToMapping": "add_to_mapping",
    "removeFromMapping": "remove_from_mapping",
    # rmgpy.molecule.molecule
    "__changeBond": "_change_bond",
    "identifyRingMembership": "identify_ring_membership",
    "getDeterministicSmallestSetOfSmallestRings": "get_deterministic_sssr",
    # rmgpy.molecule.group
    "__formBond": "_form_bond",
    "__breakBond": "_break_bond",
    "__gainRadical": "_gain_radical",
    "__loseRadical": "_lose_radical",
    "__gainPair": "_gain_pair",
    "__losePair": "_lose_pair",
    # rmgpy.pdep.collision
    "getAlpha": "get_alpha",
    "generateCollisionMatrix": "generate_collision_matrix",
    "calculate_collision_efficiency": "calculate_collision_efficiency",
    # rmgpy.pdep.draw
    "__getEnergyRange": "_get_energy_range",
    "__useStructureForLabel": "_use_structure_for_label",
    "__getTextSize": "_get_text_size",
    "__drawText": "_draw_text",
    "__getLabelSize": "_get_label_size",
    "__drawLabel": "_draw_label",
    # rmgpy.pdep.network
    "__getEnergyGrains": "_get_energy_grains",
    "printSummary": "log_summary",
    # rmgpy.qm.gaussian
    "testReady": "test_ready",
    "verifyOutputFile": "verify_output_file",
    "inputFileKeywords": "input_file_keywords",
    "writeInputFile": "write_input_file",
    "generateQMData": "generate_qm_data",
    "getParser": "get_parser",
    # rmgpy.qm.main
    "checkAllSet": "check_all_set",
    "setDefaultOutputDirectory": "set_default_output_directory",
    "checkReady": "check_ready",
    "checkPaths": "check_paths",
    "runJobs": "run_jobs",
    # rmgpy.qm.molecule
    "getFilePath": "get_file_path",
    "getCrudeMolFilePath": "get_crude_mol_file_path",
    "getRefinedMolFilePath": "get_refined_mol_file_path",
    "generateRDKitGeometries": "generate_rdkit_geometries",
    "saveCoordinatesFromRDMol": "save_coordinates_from_rdmol",
    "saveCoordinatesFromQMData": "save_coordinates_from_qm_data",
    "getThermoFilePath": "get_thermo_file_path",
    "createGeometry": "create_geometry",
    "saveThermoData": "save_thermo_data",
    "loadThermoData": "load_thermo_data",
    "getInChiKeyAug": "get_augmented_inchi_key",
    "getMolFilePathForCalculation": "get_mol_file_path_for_calculation",
    "determinePointGroup": "determine_point_group",
    "calculateChiralityCorrection": "calculate_chirality_correction",
    "calculateThermoData": "calculate_thermo_data",
    # rmgpy.qm.qmdata
    "testValid": "test_valid",
    # rmgpy.qm.qmverifier
    "checkForInChIKeyCollision": "check_for_inchi_key_collision",
    "successfulJobExists": "successful_job_exists",
    # rmgpy.statmech.ndTorsions
    "getTorsions": "get_torsions",
    "readScan": "read_scan",
    "readGjf": "read_gjf",
    "writeXYZ": "write_xyz",
    "writePes": "write_pes",
    "writeInp": "write_inp",
    "getIcsFile": "get_ics_file",
    "fitFourier": "fit_fourier",
    "getSplistfile": "get_splist_file",
    "getEigvals": "get_eigvals",
    "readEigvals": "read_eigvals",
    "calcPartitionFunction": "calc_partition_function",
    "getFrequencies": "get_frequencies",
    # rmgpy.statmech.torsion
    "getRotationalConstantEnergy": "get_rotational_constant_energy",
    "getFrequency": "get_frequency",
    "getLevelEnergy": "get_level_energy",
    "getLevelDegeneracy": "get_level_degeneracy",
    "solveSchrodingerEquation": "solve_schrodinger_equation",
    "getHamiltonian": "get_hamiltonian",
    "getPotential": "get_potential",
    "fitFourierPotentialToData": "fit_fourier_potential_to_data",
    "fitCosinePotentialToData": "fit_cosine_potential_to_data",
    # rmgpy.thermo.wilhoit
    "__residual": "_residual",
    "fitToDataForConstantB": "fit_to_data_for_constant_b",
    # rmgpy.stats
    "saveExecutionStatistics": "save_execution_statistics",
    "generateExecutionPlots": "generate_execution_plots",
    # rmgpy.tools.observablesRegression
    "runSimulations": "run_simulations",
    # rmgpy.tools.uncertainty
    "getUncertaintyValue": "get_uncertainty_value",
    "getPartialUncertaintyValue": "get_partial_uncertainty_value",
    "getUncertaintyFactor": "get_uncertainty_factor",
    "retrieveSaturatedSpeciesFromList": "retrieve_saturated_species_from_list",
    "extractSourcesFromModel": "extract_sources_from_model",
    "compileAllSources": "compile_all_sources",
    "assignParameterUncertainties": "assign_parameter_uncertainties",
    "sensitivityAnalysis": "sensitivity_analysis",
    "localAnalysis": "local_analysis",
}

# Function and method arguments
ARGUMENTS1 = {
    # Arekane:
    "Vlist": "v_list",
    "maximumRadicalElectrons": "maximum_radical_electrons",
    "format": "file_format",
    "F": "hessian",
    "getProjectedOutFreqs": "get_projected_out_freqs",
    # rmgpy.data.base
    "shortDesc": "short_desc",
    "longDesc": "long_desc",
    "referenceType": "reference_type",
    "nodalDistance": "nodal_distance",
    "numLabels": "num_labels",
    # rmgpy.data.rmg
    "thermoLibraries": "thermo_libraries",
    "transportLibraries": "transport_libraries",
    "reactionLibraries": "reaction_libraries",
    "seedMechanisms": "seed_mechanisms",
    "kineticsFamilies": "kinetics_families",
    "kineticsDepositories": "kinetics_depositories",
    "statmechLibraries": "statmech_libraries",
    # rmgpy.data.statmech
    "thermoModel": "thermo_model",
    # rmgpy.data.thermo
    "groupAdditivity": "group_additivity",
    "thermoDataList": "thermo_data_list",
    "trainingSet": "training_set",
    "bindingEnergies": "binding_energies",
    # rmgpy.data.kinetics.database
    "thermoDatabase": "thermo_database",
    "fixBarrierHeight": "fix_barrier_height",
    "forcePositiveBarrier": "force_positive_barrier",
    # rmgpy.data.kinetics.family
    "forwardTemplate": "forward_template",
    "forwardRecipe": "forward_recipe",
    "reverseTemplate": "reverse_template",
    "reverseRecipe": "reverse_recipe",
    "boundaryAtoms": "boundary_atoms",
    "treeDistance": "tree_distance",
    "depositoryLabels": "depository_labels",
    "returnAllKinetics": "return_all_kinetics",
    # rmgpy.data.kinetics.library
    "markDuplicates": "mark_duplicates",
    # rmgpy.data.kinetics.rules
    "rootTemplate": "root_template",
    "alreadyDone": "already_done",
    "kList": "k_list",
    # rmgpy.kinetics.model
    "otherKinetics": "other_kinetics",
    "ctReaction": "ct_reaction",
    "speciesList": "species_list",
    # rmgpy.kinetics.arrhenius
    "threeParams": "three_params",
    # rmgpy.kinetics.diffusionLimited
    "solvationDatabase": "solvation_database",
    # rmgpy.molecule.adjlist
    "saturateH": "saturate_h",
    "removeH": "remove_h",
    "removeLonePairs": "remove_lone_pairs",
    "oldStyle": "old_style",
    # rmgpy.molecule.atomtype
    "incrementBond": "increment_bond",
    "decrementBond": "decrement_bond",
    "formBond": "form_bond",
    "breakBond": "break_bond",
    "incrementRadical": "increment_radical",
    "decrementRadical": "decrement_radical",
    "incrementLonePair": "increment_lone_pair",
    "decrementLonePair": "decrement_lone_pair",
    # rmgpy.molecule.converter
    "removeHs": "remove_h",
    "returnMapping": "return_mapping",
    # rmgpy.molecule.draw
    "heavyFirst": "heavy_first",
    "drawLonePairs": "draw_lone_pairs",
    # rmgpy.molecule.element
    "chemkinName": "chemkin_name",
    # rmgpy.molecule.graph
    "saveOrder": "save_order",
    "initialMap": "initial_map",
    "startingVertex": "starting_vertex",
    # rmgpy.molecule.molecule
    "lonePairs": "lone_pairs",
    "newOrder": "new_order",
    "otherOrder": "other_order",
    "isSingle": "is_single",
    "isDouble": "is_double",
    "isTriple": "is_triple",
    "isQuadruple": "is_quadruple",
    "isBenzene": "is_benzene",
    "isHydrogenBond": "is_hydrogen_bond",
    "radicalElectrons": "radical_electrons",
    "logSpecies": "log_species",
    "raiseException": "raise_exception",
    "generateInitialMap": "generate_initial_map",
    "atomicNums": "atomic_nums",
    "aromaticRings": "aromatic_rings",
    "startingAtoms": "starting_atoms",
    "distanceList": "distance_list",
    "ignoreList": "ignore_list",
    # rmgpy.molecule.group
    "connectingAtom": "connecting_atom",
    "bondOrders": "bond_orders",
    "keepIdenticalLabels": "keep_identical_labels",
    # rmgpy.pdep.collision
    "densStates": "dens_states",
    "Elist": "e_list",
    "Jlist": "j_list",
    "Ereac": "e_reac",
    # rmgpy.pdep.configuration
    "bathGas": "bath_gas",
    # rmgpy.pdep.cse
    "lumpingOrder": "lumping_order",
    # rmgpy.pdep.msc
    "efficiencyModel": "efficiency_model",
    # rmgpy.pdep.network
    "pathReactions": "path_reactions",
    "netReactions": "net_reactions",
    "Ngrains": "n_grains",
    "NJ": "n_j",
    "grainSize": "grain_size",
    "grainCount": "grain_count",
    "maximumGrainSize": "maximum_grain_size",
    "minimumGrainCount": "minimum_grain_count",
    "errorCheck": "error_check",
    # rmgpy.pdep.reaction
    "reacDensStates": "reac_dens_states",
    "prodDensStates": "prod_dens_states",
    # rmgpy.qm.symmetry
    "pointGroup": "point_group",
    "symmetryNumber": "symmetry_number",
    "uniqueID": "unique_id",
    "qmData": "qm_data",
    # rmgpy.rmg.main
    "inputFile": "input_file",
    "outputDirectory": "output_directory",
    "firstTime": "first_time",
    "rxnSysUnimolecularThreshold": "rxn_sys_unimol_threshold",
    "rxnSysBimolecularThreshold": "rxn_sys_bimol_threshold",
    "rxnSysTrimolecularThreshold": "rxn_sys_trimol_threshold",
    "skipUpdate": "skip_update",
    "modulePath": "module_path",
    "reactionSystem": "reaction_system",
    # rmgpy.rmg.model
    "checkForExisting": "check_existing",
    "checkExisting": "check_existing",
    "generateThermo": "generate_thermo",
    "newObject": "new_object",
    "reactEdge": "react_edge",
    "unimolecularReact": "unimolecular_react",
    "bimolecularReact": "bimolecular_react",
    "trimolecularReact": "trimolecular_react",
    "newSurfaceSpecies": "new_surface_species",
    "newReactions": "new_reactions",
    "newSpecies": "new_species",
    "pdepNetwork": "pdep_network",
    "newCoreSpecies": "new_core_species",
    "newCoreReactions": "new_core_reactions",
    "newEdgeSpecies": "new_edge_species",
    "newEdgeReactions": "new_edge_reactions",
    "reactionsMovedFromEdge": "reactions_moved_from_edge",
    "toleranceThermoKeepSpeciesInEdge": "thermo_tol_keep_spc_in_edge",
    "minCoreSizeForPrune": "min_core_size_for_prune",
    "maximumEdgeSpecies": "maximum_edge_species",
    "reactionSystems": "reaction_systems",
    "minSpeciesExistIterationsForPrune": "min_species_exist_iterations_for_prune",
    "seedMechanism": "seed_mechanism",
    "reactionLibrary": "reaction_library",
    # rmgpy.rmg.output
    "partCoreEdge": "part_core_edge",
    # rmgpy.rmg.pdep
    "pdepSettings": "pdep_settings",
    # rmgpy.rmg.settings
    "toleranceKeepInEdge": "tol_keep_in_edge",
    "toleranceMoveToCore": "tol_move_to_core",
    "toleranceInterruptSimulation": "tol_interrupt_simulation",
    "toleranceMoveEdgeReactionToCore": "tol_move_edge_rxn_to_core",
    "filterThreshold": "filter_threshold",
    "ignoreOverallFluxCriterion": "ignore_overall_flux_criterion",
    "toleranceMoveEdgeReactionToSurface": "tol_move_edge_rxn_to_surface",
    "toleranceMoveSurfaceSpeciesToCore": "tol_move_surface_spc_to_core",
    "toleranceMoveSurfaceReactionToCore": "tol_move_surface_rxn_to_core",
    "terminateAtMaxObjects": "terminate_at_max_objects",
    "dynamicsTimeScale": "dynamics_time_scale",
    "toleranceBranchReactionToCore": "tol_branch_rxn_to_core",
    "branchingIndex": "branching_index",
    "branchingRatioMax": "branching_ratio_max",
    "toleranceMoveEdgeReactionToSurfaceInterrupt": "tol_move_edge_rxn_to_surface_interrupt",
    "toleranceMoveEdgeReactionToCoreInterrupt": "tol_move_edge_rxn_to_core_interrupt",
    "maxNumSpecies": "max_num_species",
    "maxNumObjsPerIter": "max_num_objects_per_iter",
    # rmgpy.solver.base
    "sensitiveSpecies": "sensitive_species",
    "sensitivityThreshold": "sensitivity_threshold",
    "coreSpecies": "core_species",
    "coreReactions": "core_reactions",
    "edgeSpecies": "edge_species",
    "edgeReactions": "edge_reactions",
    "surfaceSpecies": "surface_species",
    "surfaceReactions": "surface_reactions",
    "pdepNetworks": "pdep_networks",
    "filterReactions": "filter_reactions",
    "newSurfaceReactions": "new_surface_reactions",
    "newSurfaceReactionInds": "new_surface_reaction_inds",
    "sensWorksheet": "sens_worksheet",
    "modelSettings": "model_settings",
    "simulatorSettings": "simulator_settings",
    "charRate": "char_rate",
    "speciesRate": "species_rate",
    "maxDifLnAccumNum": "max_dif_ln_accum_num",
    "networkRate": "network_rate",
    "speciesIndex": "species_index",
    # rmgpy.solver.simple
    "initialMoleFractions": "initial_mole_fractions",
    "nSims": "n_sims",
    "sensConditions": "sens_conditions",
    # rmgpy.solver.liquid
    "initialConcentrations": "initial_concentrations",
    "constSPCNames": "const_spc_names",
    # rmgpy.solver.surface
    "initialP": "P_initial",
    "initialGasMoleFractions": "initial_gas_mole_fractions",
    "initialSurfaceCoverages": "initial_surface_coverages",
    "surfaceVolumeRatio": "surface_volume_ratio",
    # rmgpy.statmech.conformer
    "spinMultiplicity": "spin_multiplicity",
    "opticalIsomers": "optical_isomers",
    # rmgpy.statmech.torsion
    "Nbasis": "n_basis",
    "sumStates0": "sum_states_0",
    "densStates0": "dens_states_0",
    # rmgpy.species
    "transportData": "transport_data",
    "molecularWeight": "molecular_weight",
    "energyTransferModel": "energy_transfer_model",
    "creationIteration": "creation_iteration",
    "explicitlyAllowed": "explicitly_allowed",
    "useChemkinIdentifier": "use_chemkin_identifier",
    "solventName": "solvent_name",
    # rmgpy.reaction
    "eitherDirection": "either_direction",
    "checkIdentical": "check_identical",
    "checkOnlyLabel": "check_only_label",
    "checkTemplateRxnProducts": "check_template_rxn_products",
    "surfaceSiteDensity": "surface_site_density",
    "forcePositive": "force_positive",
    "kForward": "k_forward",
    "reverseUnits": "reverse_units",
    # rmgpy.chemkin
    "speciesDict": "species_dict",
    "reactionList": "reaction_list",
    "javaLibrary": "java_library",
    "reactionModel": "reaction_model",
    "dictionaryPath": "dictionary_path",
    "transportPath": "transport_path",
    "thermoPath": "thermo_path",
    "saveEdgeSpecies": "save_edge_species",
    "checkForDuplicates": "check_for_duplicates",
    "elementCounts": "element_counts",
    "readComments": "read_comments",
    "useChemkinNames": "use_chemkin_names",
    "checkDuplicates": "check_duplicates",
}

ARGUMENTS2 = {
    # rmgpy.data.base
    "numParameters": "num_parameters",
    "nodeOther": "node_other",
    "parentNode": "parent_node",
    "childNode": "child_node",
    "groupToRemove": "group_to_remove",
    # rmgpy.data.solvation
    "solventViscosity": "solvent_viscosity",
    "saturatedStruct": "saturated_struct",
    "addedToRadicals": "added_to_radicals",
    "addedToPairs": "added_to_pairs",
    "soluteData": "solute_data",
    "solventData": "solvent_data",
    "solventStructure": "solvent_structure",
    # rmgpy.data.thermo
    "thermoData1": "thermo_data1",
    "thermoData2": "thermo_data2",
    "stableThermoEstimator": "stable_thermo_estimator",
    "thermoData": "thermo_data",
    # rmgpy.data.transport
    "groupData": "group_data",
    # rmgpy.data.kinetics.family
    "doForward": "forward",
    "reactantStructures": "reactant_structures",
    "templateReactant": "template_reactant",
    "kineticsList": "kinetics_list",
    "templateLabels": "template_labels",
    "templateRxnMap": "template_rxn_map",
    "minSplitableEntryNum": "min_splitable_entry_num",
    "minRxnsToSpawn": "min_rxns_to_spawn",
    "maxBatchSize": "max_batch_size",
    "outlierFraction": "outlier_fraction",
    "stratumNum": "stratum_num",
    "maxRxnsToReoptNode": "max_rxns_to_reopt_node",
    "fixLabels": "fix_labels",
    "exactMatchesOnly": "exact_matches_only",
    "getReverse": "get_reverse",
    "testRxnInds": "test_rxn_inds",
    "keepRoot": "keep_root",
    "removeDegeneracy": "remove_degeneracy",
    "estimateThermo": "estimate_thermo",
    "templateLabel": "template_label",
    "childConn": "child_conn",
    # rmgpy.data.kinetics.groups
    "referenceKinetics": "reference_kinetics",
    # rmgpy.molecule.vf2
    "initialMapping": "initial_mapping",
    "findAll": "find_all",
    "callDepth": "call_depth",
    # rmgpy.molecule.molecule
    "InChI": "inchi",
    "SMILES": "smiles",
    # rmgpy.molecule.group
    "atomList": "atom_list",
    "R": "r",
    "atmInd": "atm_ind",
    "atmInd2": "atm_ind2",
    "Nsplits": "n_splits",
    "Run": "r_un",
    "Rbonds": "r_bonds",
    # rmgpy.statmech.ndTorsions
    "calcPath": "calc_path",
    "isLinear": "is_linear",
    "isTS": "is_ts",
    # rmgpy.thermo.thermoengine
    "thermoClass": "thermo_class",
    # rmgpy.thermo.wilhoit
    "contCons": "cont_cons",
    # rmgpy.quantity
    "uncertaintyType": "uncertainty_type",
    "commonUnits": "common_units",
    "extraDimensionality": "extra_dimensionality",
    # rmgpy.tools.canteraModel
    "reactorType": "reactor_type",
    "reactionTime": "reaction_time",
    "molFrac": "mol_frac",
    "reactorTypeList": "reactor_type_list",
    "reactionTimeList": "reaction_time_list",
    "molFracList": "mol_frac_list",
    "reactionMap": "reaction_map",
    "chemkinFile": "chemkin_file",
    "transportFile": "transport_file",
    "rmgReactionIndex": "rmg_reaction_index",
    "rmgReaction": "rmg_reaction",
    "rmgSpeciesIndex": "rmg_species_index",
    "rmgSpecies": "rmg_species",
    "topSpecies": "top_species",
    "topSensitiveReactions": "top_sensitive_reactions",
    "userList": "user_list",
    "RMGList": "rmg_list",
    "ctSpec1": "ct_spec1",
    "ctSpec2": "ct_spec2",
    "ctRxn1": "ct_rxn1",
    "ctRxn2": "ct_rxn2",
    "checkID": "check_id",
    # rmgpy.tools.fluxdiagram
    "reactionRates": "reaction_rates",
    "centralSpeciesList": "central_species_list",
    "speciesDirectory": "species_directory",
    "outputFile": "output_file",
    "savePath": "save_path",
    "speciesPath": "species_path",
    "chemkinOutput": "chemkin_output",
    "saveStates": "save_states",
    "readStates": "read_states",
    # rmgpy.tools.isotopes
    "useOriginalReactions": "use_original_reactions",
    "kineticIsotopeEffect": "kinetic_isotope_effect",
    "maximumIsotopicAtoms": "maximum_isotopic_atoms",
    # rmgpy.tools.loader
    "generateImages": "generate_images",
    "useJava": "use_java",
    # rmgpy.tools.observablesRegression
    "oldDir": "old_dir",
    "newDir": "new_dir",
    "exptData": "expt_data",
    # rmgpy.tools.plot
    "reactionSystemIndex": "reaction_system_index",
    "sensitiveSpeciesList": "sensitive_species_list",
    "xArray": "x_array",
    "yArray": "y_array",
    "xValue": "x_value",
    "xVar": "x_var",
    "yVar": "y_var",
    "csvFile": "csv_file",
    "numSpecies": "num_species",
    "numReactions": "num_reactions",
    "totalVariance": "total_variance",
    # rmgpy.tools.simulate
    "diffusionLimited": "diffusion_limited",
    "dictFile": "dict_file",
    # rmgpy.tools.uncertainty
    "corrSourceType": "corr_source_type",
    "corrParam": "corr_param",
    "corrGroupType": "corr_group_type",
    "corrFamily": "corr_family",
    "gParamEngine": "g_param_engine",
    "kParamEngine": "k_param_engine",
    "chemkinPath": "chemkin_path",
    "terminationTime": "termination_time",
}

# Names which are risky to replace using regex
DANGEROUS = ["SMILES", "InChI", "R", "Run", "atomTypes", "format", "F"]


def main(path, write=False, nobackups=False):
    """
    Function to analyze and update files with potential name replacements.

    Args:
        path: path to the file to be analyzed
        write: save the fixed file to disk
        nobackups: disable saving a backup of the original file
    """

    if write and not nobackups:
        dirname, filename = os.path.split(path)
        backup_path = os.path.join(dirname, filename + ".backup")
        shutil.copyfile(path, backup_path)

    with open(path, "r") as f:
        original = f.read()

    fixed = original

    for pattern, newname in tqdm(replacements):
        result = find_name_in_line(pattern, fixed, replacement=newname)
        if result:
            fixed = result

    diff = difflib.unified_diff(
        original.splitlines(keepends=True),
        fixed.splitlines(keepends=True),
        fromfile=path,
        tofile=path,
    )

    if diff:
        sys.stdout.writelines(diff)
        if write:
            with open(path, "w") as f:
                f.writelines(fixed)
    else:
        print("No changes detected")


def compile_regex(replacements, args=False, attr=False, words=True, avoid_danger=True):
    """
    Compile regex expressions for the given replacements. By default, checks for
    the name to be replaced as a word, i.e. surrounded by non-word characters.

    Args:
        replacements: dictionary of oldname, newname pairs
        args: require equal sign following name and no period before name
        attr: require period preceding name
        words: require name to be a standalone word
        avoid_danger: do not replace names in the DANGEROUS list
    """
    patterns = []
    for oldname, newname in replacements.items():
        if avoid_danger and oldname in DANGEROUS:
            continue
        pattern = r""
        if words:
            if args:
                # Require no period or word characters preceding the name
                pattern += r"(?<![\w\.])"
            elif attr:
                # Require period preceding the name
                pattern += r"(?<=\.)"
            else:
                # Require no word characters preceding the name
                pattern += r"(?<!\w)"
        pattern += oldname
        if words:
            if args:
                # Require equal sign following the name
                pattern += r"(?=\s*=\s*)"
            else:
                # Require no word characters following the name
                pattern += r"(?!\w)"
        patterns.append((re.compile(pattern), newname))

    return patterns


def find_name_in_line(pattern, line, replacement=None):
    """
    Use regex to replace the name in the line.

    Args:
        pattern: regex pattern to search for
        line: string in which to search for pattern
        replacement: if provided, text to replace the pattern with
    """
    if re.search(pattern, line):
        if replacement is not None:
            return re.sub(pattern, replacement, line)
        else:
            return True
    else:
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("filename", type=str, nargs="+")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-0", "--both-stages", action="store_true", help="check all possible replacements")
    group.add_argument("-1", "--stage1", action="store_true", help="check common replacements")
    group.add_argument("-2", "--stage2", action="store_true", help="check uncommon replacements")

    parser.add_argument("-w", "--write", action="store_true", help="perform replacements and write to file")
    parser.add_argument("-n", "--nobackups", action="store_true", help="do not save a backup file")

    parser.add_argument("-x", "--ignore-danger", action="store_true", help="do not avoid dangerous replacements")
    parser.add_argument("-y", "--all-args", action="store_true", help="disable `=` requirement for arg replacements")
    parser.add_argument("-z", "--all-attr", action="store_true", help="disable `.` requirement for attr replacements")
    parser.add_argument("-a", "--all", action="store_true", help="look for all appearances of names, implies `-xyz`")

    arguments = parser.parse_args()

    # Compile regex expressions
    args = not (arguments.all or arguments.all_args)
    attr = not (arguments.all or arguments.all_attr)
    avoid_danger = not (arguments.all or arguments.ignore_danger)
    words = not arguments.all

    replacements = []
    if arguments.both_stages or arguments.stage1:
        replacements.extend(compile_regex(MODULES, args=False, attr=False, words=words, avoid_danger=avoid_danger))
        replacements.extend(compile_regex(GLOBALS1, args=False, attr=False, words=words, avoid_danger=avoid_danger))
        replacements.extend(compile_regex(METHODS1, args=False, attr=attr, words=words, avoid_danger=avoid_danger))
        replacements.extend(compile_regex(ATTRIBUTES1, args=False, attr=attr, words=words, avoid_danger=avoid_danger))
        replacements.extend(compile_regex(ARGUMENTS1, args=args, attr=False, words=words, avoid_danger=avoid_danger))
    if arguments.both_stages or arguments.stage2:
        replacements.extend(compile_regex(GLOBALS2, args=False, attr=False, words=words, avoid_danger=avoid_danger))
        replacements.extend(compile_regex(METHODS2, args=False, attr=attr, words=words, avoid_danger=avoid_danger))
        replacements.extend(compile_regex(ATTRIBUTES2, args=False, attr=attr, words=words, avoid_danger=avoid_danger))
        replacements.extend(compile_regex(ARGUMENTS2, args=args, attr=False, words=words, avoid_danger=avoid_danger))

    for fname in arguments.filename:
        filepath = os.path.abspath(fname)
        print("Processing {0}".format(filepath))
        main(
            filepath,
            write=arguments.write,
            nobackups=arguments.nobackups,
        )
