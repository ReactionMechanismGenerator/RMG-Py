#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
import itertools
import logging

from rmgpy.data.rmg import getDB
from multiprocessing import Pool

################################################################################
def react(spc_tuples, procnum=1):
    """
    Generate reactions between the species in the
    list of species tuples for all the reaction families available.

    For each tuple of one or more Species objects [(spc1,), (spc2, spc3), ...]
    the following is done:

    A list of tuples is created for each resonance isomer of the species.
    Each tuple consists of (Molecule, index) with the index the species index of the Species object.

    Possible combinations between the first spc in the tuple, and the second species in the tuple
    is obtained by taking the combinatorial product of the two generated [(Molecule, index)] lists.

    Returns a flat generator object containing the generated Reaction objects.
    """
    # Execute multiprocessing map. It blocks until the result is ready.
    # This method chops the iterable into a number of chunks which it
    # submits to the process pool as separate tasks.
    if procnum == 1:
        logging.info('For reaction generation {0} process is used.'.format(procnum))
        reactions = map(_react_species_star, spc_tuples)
    else:
        logging.info('For reaction generation {0} processes are used.'.format(procnum))
        p = Pool(processes=procnum)
        reactions = p.map(_react_species_star, spc_tuples)
        p.close()
        p.join()

    return itertools.chain.from_iterable(reactions)

def _react_species_star(args):
    """Wrapper to unpack zipped arguments for use with map"""
    return react_species(*args)


def react_species(species_tuple, only_families=None):
    """
    Given a tuple of Species objects, generates all possible reactions
    from the loaded reaction families and combines degenerate reactions.
    """

    species_tuple = tuple([spc.copy(deep=True) for spc in species_tuple])

    reactions = getDB('kinetics').generate_reactions_from_families(species_tuple, only_families=only_families)

    return reactions


def react_all(core_spc_list, numOldCoreSpecies, unimolecularReact, bimolecularReact, trimolecularReact=None, procnum=1):
    """
    Reacts the core species list via uni-, bi-, and trimolecular
    reactions and splits reaction families per task for improved load balancing in parallel runs.
    """

    families = getDB('kinetics').families.keys()

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
    for i in xrange(numOldCoreSpecies):
        fam_leftovers = []
        for k, family in enumerate(families):
            # Find reactions involving the species that are unimolecular
            if core_spc_list[i].reactive and unimolecularReact[i, k]:
                if family in major_families and len(core_spc_list[i].molecule[0].atoms) > min_atoms:
                    spc_tuples.append(((core_spc_list[i], ), family))
                else:
                    fam_leftovers.append(family)
        if fam_leftovers:
            spc_tuples.append(((core_spc_list[i], ), fam_leftovers))

    for i in xrange(numOldCoreSpecies):
        for j in xrange(i, numOldCoreSpecies):
            fam_leftovers = []
            for k, family in enumerate(families):
                # Find reactions involving the species that are bimolecular
                # This includes a species reacting with itself (if its own concentration is high enough)
                if bimolecularReact[i, j, k] and core_spc_list[i].reactive and core_spc_list[j].reactive:
                    if family in major_families and len(core_spc_list[i].molecule[0].atoms) > min_atoms or \
                            len(core_spc_list[j].molecule[0].atoms) > min_atoms:
                        spc_tuples.append(((core_spc_list[i], core_spc_list[j]), family))
                    else:
                        fam_leftovers.append(family)
            if fam_leftovers:
                spc_tuples.append(((core_spc_list[i], core_spc_list[j]), fam_leftovers))

    if trimolecularReact is not None:
        for i in xrange(numOldCoreSpecies):
            for j in xrange(i, numOldCoreSpecies):
                for k in xrange(j, numOldCoreSpecies):
                    fam_leftovers = []
                    for l, family in enumerate(families):
                        # Find reactions involving the species that are trimolecular
                        if trimolecularReact[i, j, k, l] and core_spc_list[i].reactive and \
                                core_spc_list[j].reactive and core_spc_list[k].reactive:
                            if family in major_families and len(core_spc_list[i].molecule[0].atoms) > min_atoms or \
                                    len(core_spc_list[j].molecule[0].atoms) > min_atoms or \
                                    len(core_spc_list[k].molecule[0].atoms) > min_atoms:
                                spc_tuples.append(((core_spc_list[i], core_spc_list[j],
                                                   core_spc_list[k]), family))
                            else:
                                fam_leftovers.append(family)
                    if fam_leftovers:
                        spc_tuples.append(((core_spc_list[i], core_spc_list[j], core_spc_list[k]), fam_leftovers))

    return list(react(spc_tuples, procnum))
