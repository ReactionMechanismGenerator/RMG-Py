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
This script can be used to compare two RMG-generated kinetics models. To use,
pass the chem.inp and species_dictionary.txt files to the script. The syntax
is as follows:

python diffModels.py CHEMKIN1 SPECIESDICT1 CHEMKIN2 SPECIESDICT2

Optionally, you may use the --thermo1 and/or --thermo2 flags to add separate
thermo chemkin files.

The optional --web flag is used for running this script through the RMG-website

With all the above options the syntax is as follows:

python diffModels.py CHEMKIN1 SPECIESDICT1 --thermo1 THERMO1 CHEMKIN2 SPECIESDICT2 --thermo2 THERMO2 --web

Further option flags:
======================= ====================================================================================
Flag                    Description
======================= ====================================================================================
--diffOnly              Only show species and reactions which are unique or have different values
--commonDiffOnly        Only show species and reactions present in BOTH models which have different values
======================= ====================================================================================
"""

import argparse
import logging
import math
import os

import matplotlib.pyplot as plt

from rmgpy.chemkin import load_chemkin_file
from rmgpy.rmg.model import ReactionModel
from rmgpy.rmg.output import save_diff_html
from rmgpy.kinetics.surface import StickingCoefficient, SurfaceArrhenius


################################################################################

def compare_model_kinetics(model1, model2):
    """
    Compare the kinetics of :class:`ReactionModel` objects `model1` and 
    `model2`, printing the results to stdout.
    """
    # Determine reactions that both models have in common
    common_reactions = {}
    for rxn1 in model1.reactions:
        for rxn2 in model2.reactions:
            if rxn1.is_isomorphic(rxn2):
                common_reactions[rxn1] = rxn2
                model2.reactions.remove(rxn2)
                break
    unique_reactions1 = [rxn for rxn in model1.reactions if rxn not in list(common_reactions.keys())]
    unique_reactions2 = model2.reactions

    logging.info('{0:d} reactions were found in both models:'.format(len(common_reactions)))
    for rxn in common_reactions:
        logging.info('    {0!s}'.format(rxn))
    logging.info('{0:d} reactions were only found in the first model:'.format(len(unique_reactions1)))
    for rxn in unique_reactions1:
        logging.info('    {0!s}'.format(rxn))
    logging.info('{0:d} reactions were only found in the second model:'.format(len(unique_reactions2)))
    for rxn in unique_reactions2:
        logging.info('    {0!s}'.format(rxn))

    T = 1000
    P = 1e5
    kinetics1 = []
    kinetics2 = []
    for rxn1, rxn2 in common_reactions.items():
        stick = isinstance(rxn1.kinetics, StickingCoefficient)
        if stick:
            rate1 = rxn1.get_sticking_coefficient(T)
            rate2 = rxn1.get_sticking_coefficient(T)
        else: 
            rate1 = rxn1.get_rate_coefficient(T, P)
            rate2 = rxn1.get_rate_coefficient(T, P)

        kinetics1.append(rate1)
        if rxn1.is_isomorphic(rxn2, either_direction=False):
            kinetics2.append(rate2)
        else:
            kinetics2.append(rate2 / rxn2.get_equilibrium_constant(T))
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    plt.loglog(kinetics1, kinetics2, 'o', picker=5)
    xlim = plt.xlim()
    ylim = plt.ylim()
    lim = (min(xlim[0], ylim[0]), max(xlim[1], ylim[1]))
    ax.loglog(lim, lim, '--k')
    plt.xlabel('Model 1 rate coefficient (SI units)')
    plt.ylabel('Model 2 rate coefficient (SI units)')
    plt.title('T = {0:g} K, P = {1:g} bar'.format(T, P / 1e5))
    plt.xlim(lim)
    plt.ylim(lim)

    def onpick(event):
        xdata = event.artist.get_xdata()
        ydata = event.artist.get_ydata()
        for ind in event.ind:
            logging.info(list(common_reactions.keys())[ind])
            logging.info('k(T,P) = {0:9.2e} from model 1'.format(xdata[ind]))
            logging.info('k(T,P) = {0:9.2e} from model 2'.format(ydata[ind]))
            logging.info('ratio = 10**{0:.2f}'.format(math.log10(xdata[ind] / ydata[ind])))

    connection_id = fig.canvas.mpl_connect('pick_event', onpick)

    plt.show()


def compare_model_species(model1, model2):
    """
    This function compares two RMG models and returns a list of common species (with a nested list containing
    both species objects as elements), as well as a list of unique species for each model.
    """

    common_species = []
    unique_species1 = model1.species[:]
    unique_species2 = []

    for spec2 in model2.species:
        for spec1 in unique_species1[:]:  # make a copy so you don't remove from the list you are iterating over
            if spec1.is_isomorphic(spec2):
                common_species.append([spec1, spec2])
                unique_species1.remove(spec1)
                break
        else:
            unique_species2.append(spec2)
    # Remove species in the mechanism that aren't identified (includes those called out as species
    # but not used)        
    for spec in unique_species1[:]:  # make a copy so you don't remove from the list you are iterating over
        if not len(spec.molecule):
            unique_species1.remove(spec)
            logging.warning("Removing species {!r} from model 1 because it has no molecule info".format(spec))
    for spec in unique_species2[:]:  # make a copy so you don't remove from the list you are iterating over
        if not spec.molecule:
            unique_species2.remove(spec)
            logging.warning("Removing species {!r} from model 2 because it has no molecule info".format(spec))
    return common_species, unique_species1, unique_species2


def compare_model_reactions(model1, model2):
    """
    This function compares two RMG models and returns a list of common reactions (with a nested list containing
    both reaction objects as elements), as well as a list of unique reactions for each model.
    """
    reaction_list1 = model1.reactions[:]
    reaction_list2 = model2.reactions[:]

    # remove reactions that have an unidentified species
    to_remove = []
    for reactionList in (reaction_list1, reaction_list2):
        for reaction in reactionList:
            for side in (reaction.products, reaction.reactants):
                for species in side:
                    if not species.molecule:
                        to_remove.append((reactionList, reaction))
                        logging.warning(
                            "Removing reaction {!r} that had unidentified species {!r}".format(reaction, species))
                        break
    for reactionList, reaction in to_remove:
        reactionList.remove(reaction)

    common_reactions = []
    unique_reactions1 = []
    unique_reactions2 = []
    for rxn1 in reaction_list1:
        for rxn2 in reaction_list2[:]:  # make a copy so you don't remove from the list you are iterating over
            if rxn1.is_isomorphic(rxn2):
                common_reactions.append([rxn1, rxn2])
                # Remove reaction 2 from being chosen a second time.
                # Let each reaction only appear only once in the diff comparison.
                # Otherwise this miscounts number of reactions in model 2.
                reaction_list2.remove(rxn2)
                break
    for rxn1 in reaction_list1:
        for r1, r2 in common_reactions:
            if rxn1 is r1:
                break
        else:
            unique_reactions1.append(rxn1)
    for rxn2 in reaction_list2:
        for r1, r2 in common_reactions:
            if rxn2 is r2:
                break
        else:
            unique_reactions2.append(rxn2)

    return common_reactions, unique_reactions1, unique_reactions2


def save_compare_html(outputDir, chemkin_path1, species_dict_path1, chemkin_path2, species_dict_path2,
                      read_comments1=True, read_comments2=True, surf_path1 = False, surf_path2=False):
    """
    Saves a model comparison HTML file based on two sets of chemkin and species dictionary
    files.
    """
    model1 = ReactionModel()
    model2 = ReactionModel()

    if surf_path1 and surf_path2:
        model1.species, model1.reactions = load_chemkin_file(
            chemkin_path1, species_dict_path1, read_comments=read_comments1,
            surface_path=surf_path1)
        model2.species, model2.reactions = load_chemkin_file(
            chemkin_path2, species_dict_path2, read_comments=read_comments2,
            surface_path=surf_path2)
    else: 
        if not surf_path1 or not surf_path2:
            logging.error("to compare gas+surface mechanism, both models need a surface chemkin file")
        model1.species, model1.reactions = load_chemkin_file(
            chemkin_path1, species_dict_path1, read_comments=read_comments1)
        model2.species, model2.reactions = load_chemkin_file(
            chemkin_path2, species_dict_path2, read_comments=read_comments2)
    
    common_reactions, unique_reactions1, unique_reactions2 = compare_model_reactions(model1, model2)
    common_species, unique_species1, unique_species2 = compare_model_species(model1, model2)

    output_path = outputDir + 'diff.html'
    save_diff_html(output_path, common_species, unique_species1, unique_species2, common_reactions, unique_reactions1,
                   unique_reactions2)


def enthalpy_diff(species):
    """
    Returns the enthalpy discrepancy between the same species in the two models
    """
    thermo0 = species[0].thermo
    thermo1 = species[1].thermo
    if thermo0 and thermo1:
        diff = species[0].thermo.discrepancy(species[1].thermo)
    else:
        diff = 99999999
    return -1 * diff


def kinetics_diff(reaction):
    """
    Returns some measure of the discrepancy between two reactions in a model
    """
    kinetics0 = reaction[0].kinetics
    kinetics1 = reaction[1].kinetics
    if kinetics0 and kinetics1:
        diff = reaction[0].kinetics.discrepancy(reaction[1].kinetics)
    else:
        diff = 9999999
    return -1 * diff


def identical_thermo(species_pair):
    return species_pair[0].thermo.is_identical_to(species_pair[1].thermo)


def identical_kinetics(reaction_pair):
    return reaction_pair[0].kinetics.is_identical_to(reaction_pair[1].kinetics)


################################################################################

def parse_command_line_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('chemkin1', metavar='CHEMKIN1', type=str, nargs=1,
                        help='the Chemkin file of the first model')
    parser.add_argument('speciesDict1', metavar='SPECIESDICT1', type=str, nargs=1,
                        help='the species dictionary file of the first model')
    parser.add_argument('--thermo1', metavar='THERMO1', type=str, nargs=1,
                        help='the thermo file of the first model')
    parser.add_argument('chemkin2', metavar='CHEMKIN2', type=str, nargs=1,
                        help='the Chemkin file of the second model')
    parser.add_argument('speciesDict2', metavar='SPECIESDICT2', type=str, nargs=1,
                        help='the species dictionary file of the second model')
    parser.add_argument('--thermo2', metavar='THERMO2', type=str, nargs=1,
                        help='the thermo file of the second model')
    parser.add_argument('--web', action='store_true', help='Running diff models through the RMG-website')
    parser.add_argument('--diffOnly', action='store_true', help='Do not show identical species thermo or reactions')
    parser.add_argument('--commonDiffOnly', action='store_true',
                        help='Only show species and reactions present in BOTH models which have different values')

    args = parser.parse_args()

    return args


def main():
    """
    Driver function that parses command line arguments and passes them to the execute function.
    """
    args = parse_command_line_arguments()

    chemkin1 = args.chemkin1[0]
    species_dict1 = args.speciesDict1[0]
    if args.thermo1:
        thermo1 = args.thermo1[0]
    else:
        thermo1 = None
    chemkin2 = args.chemkin2[0]
    species_dict2 = args.speciesDict2[0]
    if args.thermo2:
        thermo2 = args.thermo2[0]
    else:
        thermo2 = None

    kwargs = {
        'web': args.web,
        'wd': os.getcwd(),
        'diffOnly': args.diffOnly,
        'commonDiffOnly': args.commonDiffOnly,
    }

    execute(chemkin1, species_dict1, thermo1, chemkin2, species_dict2, thermo2, **kwargs)


def execute(chemkin1, species_dict1, thermo1, chemkin2, species_dict2, thermo2, **kwargs):
    model1 = ReactionModel()
    model1.species, model1.reactions = load_chemkin_file(chemkin1, species_dict1, thermo_path=thermo1)
    model2 = ReactionModel()
    model2.species, model2.reactions = load_chemkin_file(chemkin2, species_dict2, thermo_path=thermo2)

    common_species, unique_species1, unique_species2 = compare_model_species(model1, model2)
    common_reactions, unique_reactions1, unique_reactions2 = compare_model_reactions(model1, model2)

    try:
        diff_only = kwargs['diffOnly']
    except KeyError:
        diff_only = False

    try:
        common_diff_only = kwargs['commonDiffOnly']
    except KeyError:
        common_diff_only = False

    if diff_only or common_diff_only:
        common_species = [x for x in common_species if not identical_thermo(x)]
        common_reactions = [x for x in common_reactions if not identical_kinetics(x)]

    if common_diff_only:
        unique_species1 = []
        unique_species2 = []
        unique_reactions1 = []
        unique_reactions2 = []

    try:
        web = kwargs['web']
    except KeyError:
        web = False

    if not web:
        logging.info('{0:d} species were found in both models:'.format(len(common_species)))
        for spec1, spec2 in common_species:
            logging.info('    {0!s}'.format(spec1))
            if spec1.thermo and spec2.thermo:
                spec1.molecule[0].get_symmetry_number()
                logging.info(
                    '        {0:7.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f} {6:7.2f} {7:7.2f} {8:7.2f}'.format(
                        spec1.thermo.get_enthalpy(300) / 4184.,
                        spec1.thermo.get_entropy(300) / 4.184,
                        spec1.thermo.get_heat_capacity(300) / 4.184,
                        spec1.thermo.get_heat_capacity(400) / 4.184,
                        spec1.thermo.get_heat_capacity(500) / 4.184,
                        spec1.thermo.get_heat_capacity(600) / 4.184,
                        spec1.thermo.get_heat_capacity(800) / 4.184,
                        spec1.thermo.get_heat_capacity(1000) / 4.184,
                        spec1.thermo.get_heat_capacity(1500) / 4.184,
                    ))
                logging.info(
                    '        {0:7.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f} {6:7.2f} {7:7.2f} {8:7.2f}'.format(
                        spec2.thermo.get_enthalpy(300) / 4184.,
                        spec2.thermo.get_entropy(300) / 4.184,
                        spec2.thermo.get_heat_capacity(300) / 4.184,
                        spec2.thermo.get_heat_capacity(400) / 4.184,
                        spec2.thermo.get_heat_capacity(500) / 4.184,
                        spec2.thermo.get_heat_capacity(600) / 4.184,
                        spec2.thermo.get_heat_capacity(800) / 4.184,
                        spec2.thermo.get_heat_capacity(1000) / 4.184,
                        spec2.thermo.get_heat_capacity(1500) / 4.184,
                    ))
        logging.info('{0:d} species were only found in the first model:'.format(len(unique_species1)))
        for spec in unique_species1:
            logging.info('    {0!s}'.format(spec))
        logging.info('{0:d} species were only found in the second model:'.format(len(unique_species2)))
        for spec in unique_species2:
            logging.info('    {0!s}'.format(spec))

        logging.info('{0:d} reactions were found in both models:'.format(len(common_reactions)))
        for rxn1, rxn2 in common_reactions:
            logging.info('    {0!s}'.format(rxn1))
            stick = isinstance(rxn1.kinetics, StickingCoefficient)
            if rxn1.kinetics and rxn2.kinetics:
                if not stick:
                    logging.info('        {0:7.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f} {6:7.2f} {7:7.2f}'.format(
                        math.log10(rxn1.kinetics.get_rate_coefficient(300, 1e5)),
                        math.log10(rxn1.kinetics.get_rate_coefficient(400, 1e5)),
                        math.log10(rxn1.kinetics.get_rate_coefficient(500, 1e5)),
                        math.log10(rxn1.kinetics.get_rate_coefficient(600, 1e5)),
                        math.log10(rxn1.kinetics.get_rate_coefficient(800, 1e5)),
                        math.log10(rxn1.kinetics.get_rate_coefficient(1000, 1e5)),
                        math.log10(rxn1.kinetics.get_rate_coefficient(1500, 1e5)),
                        math.log10(rxn1.kinetics.get_rate_coefficient(2000, 1e5)),
                    ))
                    logging.info('        {0:7.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f} {6:7.2f} {7:7.2f}'.format(
                        math.log10(rxn2.kinetics.get_rate_coefficient(300, 1e5)),
                        math.log10(rxn2.kinetics.get_rate_coefficient(400, 1e5)),
                        math.log10(rxn2.kinetics.get_rate_coefficient(500, 1e5)),
                        math.log10(rxn2.kinetics.get_rate_coefficient(600, 1e5)),
                        math.log10(rxn2.kinetics.get_rate_coefficient(800, 1e5)),
                        math.log10(rxn2.kinetics.get_rate_coefficient(1000, 1e5)),
                        math.log10(rxn2.kinetics.get_rate_coefficient(1500, 1e5)),
                        math.log10(rxn2.kinetics.get_rate_coefficient(2000, 1e5)),
                    ))
                else: 
                    logging.info('        {0:7.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f} {6:7.2f} {7:7.2f}'.format(
                        math.log10(rxn1.kinetics.get_sticking_coefficient(300)),
                        math.log10(rxn1.kinetics.get_sticking_coefficient(400)),
                        math.log10(rxn1.kinetics.get_sticking_coefficient(500)),
                        math.log10(rxn1.kinetics.get_sticking_coefficient(600)),
                        math.log10(rxn1.kinetics.get_sticking_coefficient(800)),
                        math.log10(rxn1.kinetics.get_sticking_coefficient(1000)),
                        math.log10(rxn1.kinetics.get_sticking_coefficient(1500)),
                        math.log10(rxn1.kinetics.get_sticking_coefficient(2000)),
                    ))
                    logging.info('        {0:7.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f} {6:7.2f} {7:7.2f}'.format(
                        math.log10(rxn2.kinetics.get_sticking_coefficient(300)),
                        math.log10(rxn2.kinetics.get_sticking_coefficient(400)),
                        math.log10(rxn2.kinetics.get_sticking_coefficient(500)),
                        math.log10(rxn2.kinetics.get_sticking_coefficient(600)),
                        math.log10(rxn2.kinetics.get_sticking_coefficient(800)),
                        math.log10(rxn2.kinetics.get_sticking_coefficient(1000)),
                        math.log10(rxn2.kinetics.get_sticking_coefficient(1500)),
                        math.log10(rxn2.kinetics.get_sticking_coefficient(2000)),
                    ))
        logging.info('{0:d} reactions were only found in the first model:'.format(len(unique_reactions1)))
        for rxn in unique_reactions1:
            logging.info('    {0!s}'.format(rxn))
        logging.info('{0:d} reactions were only found in the second model:'.format(len(unique_reactions2)))
        for rxn in unique_reactions2:
            logging.info('    {0!s}'.format(rxn))

    logging.info("Saving output in diff.html")

    try:
        wd = kwargs['wd']
    except KeyError:
        wd = os.getcwd()

    output_path = os.path.join(wd, 'diff.html')
    save_diff_html(output_path, common_species, unique_species1, unique_species2,
                   common_reactions, unique_reactions1, unique_reactions2)
    logging.info("Finished!")

    return common_species, unique_species1, unique_species2, common_reactions, unique_reactions1, unique_reactions2
