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
import resource
import psutil
import os
from sys import platform

from rmgpy.data.rmg import getDB
from multiprocessing import Pool

################################################################################

def determine_procnum_from_RAM():
    """
    Get available RAM (GB)and procnum dependent on OS.
    """
    from rmgpy.rmg.main import maxproc

    if platform.startswith('linux'):
        # linux
        memory_available = psutil.virtual_memory().free / (1000.0 ** 3)
        memory_use = psutil.Process(os.getpid()).memory_info()[0]/(1000.0 ** 3)
        tmp = divmod(memory_available, memory_use)
        tmp2 = min(maxproc, tmp[0])
        procnum = max(1, int(tmp2))
        if procnum == 1:
            logging.info('For reaction generation {0} process is used.'.format(procnum))
        else:
            logging.info('For reaction generation {0} processes are used.'.format(procnum))
    elif platform == "darwin":
        # OS X
        memory_available = psutil.virtual_memory().available/(1000.0 ** 3)
        memory_use = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1000.0 ** 3)
        tmp = divmod(memory_available, memory_use)
        tmp2 = min(maxproc, tmp[0])
        procnum = max(1, int(tmp2))
        if procnum == 1:
            logging.info('For reaction generation {0} process is used.'.format(procnum))
        else:
            logging.info('For reaction generation {0} processes are used.'.format(procnum))
    else:
        # Everything else
        procnum = 1
        logging.info('For reaction generation {0} process is used.'.format(procnum))

    # Return the maximal number of processes for multiprocessing
    return procnum

def react(*spc_tuples):
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

    procnum = determine_procnum_from_RAM()

    # Execute multiprocessing map. It blocks until the result is ready.
    # This method chops the iterable into a number of chunks which it
    # submits to the process pool as separate tasks.
    if procnum == 1:
        reactions = map(_react_species_star, spc_tuples)
    else:
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


def react_all(core_spc_list, numOldCoreSpecies, unimolecularReact, bimolecularReact, trimolecularReact=None):
    """
    Reacts the core species list via uni-, bi-, and trimolecular
    reactions and splits reaction families per task for improved load balancing in parallel runs.
    """

    procnum = determine_procnum_from_RAM()

    # Select reactive species that can undergo unimolecular reactions:
    spc_tuples = [(core_spc_list[i],)
                  for i in xrange(numOldCoreSpecies) if (unimolecularReact[i] and core_spc_list[i].reactive)]

    for i in xrange(numOldCoreSpecies):
        for j in xrange(i, numOldCoreSpecies):
            # Find reactions involving the species that are bimolecular.
            # This includes a species reacting with itself (if its own concentration is high enough).
            if bimolecularReact[i, j]:
                if core_spc_list[i].reactive and core_spc_list[j].reactive:
                    spc_tuples.append((core_spc_list[i], core_spc_list[j]))

    if trimolecularReact is not None:
        for i in xrange(numOldCoreSpecies):
            for j in xrange(i, numOldCoreSpecies):
                for k in xrange(j, numOldCoreSpecies):
                    # Find reactions involving the species that are trimolecular.
                    if trimolecularReact[i, j, k]:
                        if core_spc_list[i].reactive and core_spc_list[j].reactive and core_spc_list[k].reactive:
                            spc_tuples.append((core_spc_list[i], core_spc_list[j], core_spc_list[k]))

    if procnum == 1:
        # React all families like normal (provide empty argument for only_families)
        spc_fam_tuples = zip(spc_tuples)
    else:
        # Identify and split families that are prone to generate many reactions into sublists.
        family_list = getDB('kinetics').families.keys()
        major_families = [
            'H_Abstraction', 'R_Recombination', 'Intra_Disproportionation', 'Intra_RH_Add_Endocyclic',
            'Singlet_Carbene_Intra_Disproportionation', 'Intra_ene_reaction', 'Disproportionation',
            '1,4_Linear_birad_scission', 'R_Addition_MultipleBond', '2+2_cycloaddition_Cd', 'Diels_alder_addition',
            'Intra_RH_Add_Exocyclic', 'Intra_Retro_Diels_alder_bicyclic', 'Intra_2+2_cycloaddition_Cd',
            'Birad_recombination', 'Intra_Diels_alder_monocyclic', '1,4_Cyclic_birad_scission', '1,2_Insertion_carbene',
        ]

        split_list = []
        leftovers = []
        for fam in family_list:
            if fam in major_families:
                split_list.append([fam])
            else:
                leftovers.append(fam)
        split_list.append(leftovers)

        # Only employ family splitting for reactants that have a larger number than min_atoms
        min_atoms = 10
        spc_fam_tuples = []
        for i, spc_tuple in enumerate(spc_tuples):
            if any([len(spc.molecule[0].atoms) > min_atoms for spc in spc_tuple]):
                for item in split_list:
                    spc_fam_tuples.append((spc_tuple, item))
            else:
                spc_fam_tuples.append((spc_tuple,))

    return list(react(*spc_fam_tuples))

