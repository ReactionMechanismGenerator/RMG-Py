#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
Contains functions for generating reactions.
"""
import logging
from multiprocessing import Pool

from rmgpy.data.rmg import get_db


################################################################################


def react(spc_fam_tuples, procnum=1):
    """
    Generate reactions between the species in the list of species-family tuples
    for the optionally specified reaction families.

    Each item should be a tuple of a species tuple and an optional family list:
        [((spc1,), [family1, ...]), ((spc2, spc3), [family2, ...]), ...]
            OR
        [((spc1,),), ((spc2, spc3),), ...]

    If no family list is provided, all of the loaded families are considered.

    Args:
        spc_fam_tuples (list): list of tuples for reaction generation
        procnum (int, optional): number of processors used for reaction generation

    Returns:
        list of lists of reactions generated from each species tuple (note: empty lists are possible)
    """
    # Execute multiprocessing map. It blocks until the result is ready.
    # This method chops the iterable into a number of chunks which it
    # submits to the process pool as separate tasks.
    if procnum == 1:
        logging.info('For reaction generation {0} process is used.'.format(procnum))
        reactions = list(map(_react_species_star, spc_fam_tuples))
    else:
        logging.info('For reaction generation {0} processes are used.'.format(procnum))
        p = Pool(processes=procnum)
        reactions = p.map(_react_species_star, spc_fam_tuples)
        p.close()
        p.join()

    return reactions


def _react_species_star(args):
    """Wrapper to unpack zipped arguments for use with map"""
    return react_species(*args)


def react_species(species_tuple, only_families=None):
    """
    Given a tuple of Species objects, generates all possible reactions
    from the loaded reaction families and combines degenerate reactions.

    Args:
        species_tuple (tuple): tuple of 1-3 Species objects to react together
        only_families (list, optional): list of reaction families to consider

    Returns:
        list of generated reactions
    """
    reactions = get_db('kinetics').generate_reactions_from_families(species_tuple, only_families=only_families)

    return reactions


def react_all(core_spc_list, num_old_core_species, unimolecular_react, bimolecular_react, trimolecular_react=None, procnum=1):
    """
    Reacts the core species list via uni-, bi-, and trimolecular reactions.

    For parallel processing, reaction families are split per task for improved
    load balancing. This is currently hard-coded using reaction family labels.
    Species tuples are only generated if the associated family allows unimolecular or bimolecular reactions, respectively.

    Note:
        The bi- and trimolecular_react flags are only set for the upper part of the tensor.

    Args:
        core_spc_list (list): list of all core species
        num_old_core_species (int): current number of core species in the model
        unimolecular_react (np.ndarray): reaction filter flags indicating which species to react unimolecularly
        bimolecular_react (np.ndarray): reaction filter flags indicating which species to react bimolecularly
        trimolecular_react (np.ndarray, optional): reaction filter flags indicating which species to react trimolecularly
        procnum (int, optional): number of processors used for reaction generation

    Returns:
        a list of lists of reactions generated from each species tuple
        a list of species tuples corresponding to each list of reactions
    """
    # Retrieve the molecularity for each family
    fam_molecularity = {}
    for label, family in get_db('kinetics').families.items():
        fam_molecularity[label] = [
            family.reactant_num or len(family.forward_template.reactants),
            family.product_num or len(family.forward_template.products),
        ]

    # List of families that should not react together as they are likely to generate a lot of reactions and
    # therefore negatively impact load balancing for multiprocessing
    major_families = [
        'H_Abstraction', 'R_Recombination', 'Intra_Disproportionation', 'Intra_RH_Add_Endocyclic',
        'Singlet_Carbene_Intra_Disproportionation', 'Intra_ene_reaction', 'Disproportionation',
        '1,4_Linear_birad_scission', 'R_Addition_MultipleBond', '2+2_cycloaddition_Cd', 'Diels_alder_addition',
        'Intra_RH_Add_Exocyclic', 'Intra_Retro_Diels_alder_bicyclic', 'Intra_2+2_cycloaddition_Cd',
        'Birad_recombination', 'Intra_Diels_alder_monocyclic', '1,4_Cyclic_birad_scission', '1,2_Insertion_carbene',
    ]

    # Only employ family splitting for reactants that have a larger number than min_atoms
    min_atoms = 10

    # Select reactive species that can undergo unimolecular reactions:
    spc_tuples = []
    for i in range(num_old_core_species):
        family_list = []
        for k, (family, molecularity) in enumerate(fam_molecularity.items()):
            # Find reactions involving the species that are unimolecular and only generate tuple if family is not
            # bimolecular in forward and backward direction.
            if 1 in molecularity:
                if core_spc_list[i].reactive and unimolecular_react[i, k]:
                    if procnum > 1 and family in major_families and len(core_spc_list[i].molecule[0].atoms) > min_atoms:
                        spc_tuples.append(((core_spc_list[i], ), family))
                    else:
                        family_list.append(family)
        if family_list:
            spc_tuples.append(((core_spc_list[i], ), family_list))

    for i in range(num_old_core_species):
        for j in range(i, num_old_core_species):
            family_list = []
            for k, (family, molecularity) in enumerate(fam_molecularity.items()):
                # Find reactions involving the species that are bimolecular
                # This includes a species reacting with itself (if its own concentration is high enough)
                if 2 in molecularity:
                    if bimolecular_react[i, j, k] and core_spc_list[i].reactive and core_spc_list[j].reactive:
                        if (procnum > 1 and family in major_families and
                                (len(core_spc_list[i].molecule[0].atoms) > min_atoms or
                                 len(core_spc_list[j].molecule[0].atoms) > min_atoms)):
                            spc_tuples.append(((core_spc_list[i], core_spc_list[j]), family))
                        else:
                            family_list.append(family)
            if family_list:
                spc_tuples.append(((core_spc_list[i], core_spc_list[j]), family_list))

    if trimolecular_react is not None:
        for i in range(num_old_core_species):
            for j in range(i, num_old_core_species):
                for k in range(j, num_old_core_species):
                    family_list = []
                    for l, (family, molecularity) in enumerate(fam_molecularity.items()):
                        # Find reactions involving the species that are trimolecular
                        if 3 in molecularity:
                            if (trimolecular_react[i,j,k,l] and core_spc_list[i].reactive and
                                    core_spc_list[j].reactive and core_spc_list[k].reactive):
                                if (procnum > 1 and family in major_families and
                                        (len(core_spc_list[i].molecule[0].atoms) > min_atoms or
                                         len(core_spc_list[j].molecule[0].atoms) > min_atoms or
                                         len(core_spc_list[k].molecule[0].atoms) > min_atoms)):
                                    spc_tuples.append(((core_spc_list[i], core_spc_list[j], core_spc_list[k]),
                                                       family))
                                else:
                                    family_list.append(family)
                    if family_list:
                        spc_tuples.append(((core_spc_list[i], core_spc_list[j], core_spc_list[k]), family_list))

    return react(spc_tuples, procnum), [fam_tuple[0] for fam_tuple in spc_tuples]
