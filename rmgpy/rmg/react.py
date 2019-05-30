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

def react(*spcTuples):
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

    from rmgpy.rmg.main import maxproc

    # Get available RAM (GB)and procnum dependent on OS.
    if platform.startswith('linux'):
        # linux
        memory_available = psutil.virtual_memory().free / (1000.0 ** 3)
        memory_use = psutil.Process(os.getpid()).memory_info()[0]/(1000.0 ** 3)
        tmp = divmod(memory_available, memory_use)
        tmp2 = min(maxproc, tmp[0])
        procnum = max(1, int(tmp2))
        if maxproc == 1:
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
        if maxproc == 1:
            logging.info('For reaction generation {0} process is used.'.format(procnum))
        else:
            logging.info('For reaction generation {0} processes are used.'.format(procnum))
    else:
        # Everything else
        procnum = 1
        logging.info('For reaction generation {0} process is used.'.format(procnum))

    # Execute multiprocessing map. It blocks until the result is ready.
    # This method chops the iterable into a number of chunks which it
    # submits to the process pool as separate tasks.
    if procnum == 1:
        reactions = map(
                    reactSpecies,
                    spc_tuples)
    else:
        p = Pool(processes=procnum)
        
        reactions = p.map(
                    reactSpecies,
	            spc_tuples)
	
        p.close()
        p.join()


    return itertools.chain.from_iterable(reactions)


def reactSpecies(speciesTupleTmp):
    """
    Given a tuple of Species objects, generates all possible reactions
    from the loaded reaction families and combines degenerate reactions.
    """

    speciesTuple = speciesTupleTmp[0:-1]
    own_families = speciesTupleTmp[-1]

    speciesTuple = tuple([spc.copy(deep=True) for spc in speciesTuple])

    reactions = getDB('kinetics').generate_reactions_from_families(speciesTuple, only_families=own_families)

    return reactions


def reactAll(coreSpcList, numOldCoreSpecies, unimolecularReact, bimolecularReact, trimolecularReact=None):
    """
    Reacts the core species list via uni-, bi-, and trimolecular
    reactions.
    """

    from rmgpy.rmg.main import maxproc

    # Load kineticsFamilies to be added to reactant tuple to allow for improved load balancing 
    # in parallel jobs.
    split_listOrig = []
    split_list_tmp = []
    for key in getDB('kinetics').families:
        split_listOrig.append(key)
        split_list_tmp.append(key)

    if maxproc == 1:
        # Select reactive species that can undergo unimolecular reactions:
        spc_tuplestmp = [(core_spc_list[i], split_listOrig)
         for i in xrange(numOldCoreSpecies) if (unimolecularReact[i] and core_spc_list[i].reactive)]
    
        for i in xrange(numOldCoreSpecies):
            for j in xrange(i, numOldCoreSpecies):
                # Find reactions involving the species that are bimolecular.
                # This includes a species reacting with itself (if its own concentration is high enough).
                if bimolecularReact[i,j]:
                    if core_spc_list[i].reactive and core_spc_list[j].reactive:
                        spc_tuplestmp.append((core_spc_list[i], core_spc_list[j], split_listOrig))
    
        if trimolecularReact is not None:
            for i in xrange(numOldCoreSpecies):
                for j in xrange(i, numOldCoreSpecies):
                    for k in xrange(j, numOldCoreSpecies):
                        # Find reactions involving the species that are trimolecular.
                        if trimolecularReact[i,j,k]:
                            if core_spc_list[i].reactive and core_spc_list[j].reactive and core_spc_list[k].reactive:
                                spc_tuplestmp.append((core_spc_list[i], core_spc_list[j], core_spc_list[k], split_listOrig))
    else:
        # Select reactive species that can undergo unimolecular reactions:
        spc_tuples = [(core_spc_list[i],)
         for i in xrange(numOldCoreSpecies) if (unimolecularReact[i] and core_spc_list[i].reactive)]
    
        for i in xrange(numOldCoreSpecies):
            for j in xrange(i, numOldCoreSpecies):
                # Find reactions involving the species that are bimolecular.
                # This includes a species reacting with itself (if its own concentration is high enough).
                if bimolecularReact[i,j]:
                    if core_spc_list[i].reactive and core_spc_list[j].reactive:
                        spc_tuples.append((core_spc_list[i], core_spc_list[j]))
    
        if trimolecularReact is not None:
            for i in xrange(numOldCoreSpecies):
                for j in xrange(i, numOldCoreSpecies):
                    for k in xrange(j, numOldCoreSpecies):
                        # Find reactions involving the species that are trimolecular.
                        if trimolecularReact[i,j,k]:
                            if core_spc_list[i].reactive and core_spc_list[j].reactive and core_spc_list[k].reactive:
                                spc_tuples.append((core_spc_list[i], core_spc_list[j], core_spc_list[k]))


        # Identify and split families that are prone to generate many reactions into sublists.
        split_list = []
        for i in range(len(split_list_tmp)):
            if split_list_tmp[i] == 'H_Abstraction':
                split_list_tmp[i] = []
                split_list.append(['H_Abstraction'])
            elif split_list_tmp[i] == 'R_Recombination':
                split_list_tmp[i] = []
                split_list.append(['R_Recombination'])
            elif split_list_tmp[i] == 'Intra_Disproportionation':
                split_list_tmp[i] = []
                split_list.append(['Intra_Disproportionation'])
            elif split_list_tmp[i] == 'Intra_RH_Add_Endocyclic':
                split_list_tmp[i] = []
                split_list.append(['Intra_RH_Add_Endocyclic'])
            elif split_list_tmp[i] == 'Singlet_Carbene_Intra_Disproportionation':
                split_list_tmp[i] = []
                split_list.append(['Singlet_Carbene_Intra_Disproportionation'])
            elif split_list_tmp[i] == 'Intra_ene_reaction':
                split_list_tmp[i] = []
                split_list.append(['Intra_ene_reaction'])
            elif split_list_tmp[i] == 'Disproportionation':
                split_list_tmp[i] = []
                split_list.append(['Disproportionation'])
            elif split_list_tmp[i] == '1,4_Linear_birad_scission':
                split_list_tmp[i] = []
                split_list.append(['1,4_Linear_birad_scission'])
            elif split_list_tmp[i] == 'R_Addition_MultipleBond':
                split_list_tmp[i] = []
                split_list.append(['R_Addition_MultipleBond'])
            elif split_list_tmp[i] == '2+2_cycloaddition_Cd':
                split_list_tmp[i] = []
                split_list.append(['2+2_cycloaddition_Cd'])
            elif split_list_tmp[i] == 'Diels_alder_addition':
                split_list_tmp[i] = []
                split_list.append(['Diels_alder_addition'])
            elif split_list_tmp[i] == 'Intra_RH_Add_Exocyclic':
                split_list_tmp[i] = []
                split_list.append(['Intra_RH_Add_Exocyclic'])
            elif split_list_tmp[i] == 'Intra_Retro_Diels_alder_bicyclic':
                split_list_tmp[i] = []
                split_list.append(['Intra_Retro_Diels_alder_bicyclic'])
            elif split_list_tmp[i] == 'Intra_2+2_cycloaddition_Cd':
                split_list_tmp[i] = []
                split_list.append(['Intra_2+2_cycloaddition_Cd'])
            elif split_list_tmp[i] == 'Birad_recombination':
                split_list_tmp[i] = []
                split_list.append(['Birad_recombination'])
            elif split_list_tmp[i] == 'Intra_Diels_alder_monocyclic':
                split_list_tmp[i] = []
                split_list.append(['Intra_Diels_alder_monocyclic'])
            elif split_list_tmp[i] == '1,4_Cyclic_birad_scission':
                split_list_tmp[i] = []
                split_list.append(['1,4_Cyclic_birad_scission'])
            elif split_list_tmp[i] == '1,2_Insertion_carbene':
                split_list_tmp[i] = []
                split_list.append(['1,2_Insertion_carbene'])
    
        # Remove empty lists from remaining split_list_tmp. It now contains only
        # families that are not mentioned above.
        split_list.append(filter(None, split_list_tmp))
    
        # Only employ family splitting for reactants that have a larger number than nAFS.
        nAFS = 10
    
        spc_tuplestmp = []
        # Append reaction families to reactant tuple.
        for tmpj in spc_tuples:
            if len(tmpj) == 1:
                if len(str(tmpj[0])) > nAFS:
                    for tmpl in split_list: 
                        tmpk = list(tmpj)
                        tmpk.append(tmpl)
                        spc_tuplestmp.append(tuple(tmpk))
                else:
                    tmpk = list(tmpj)
                    tmpk.append(split_listOrig)
                    spc_tuplestmp.append(tuple(tmpk))
            elif len(tmpj) == 2:
                if (len(str(tmpj[0])) > nAFS
                   ) or (len(str(tmpj[1])) > nAFS):
                    for tmpl in split_list:
                        tmpk = list(tmpj)
                        tmpk.append(tmpl)
                        spc_tuplestmp.append(tuple(tmpk))
                else:
                    tmpk = list(tmpj)
                    tmpk.append(split_listOrig)
                    spc_tuplestmp.append(tuple(tmpk))
            else:
                if (len(str(tmpj[0])) > nAFS
                   ) or (len(str(tmpj[1])) > nAFS
                        ) or (len(str(tmpj[2])) > nAFS):
                    for tmpl in split_list:
                        tmpk = list(tmpj)
                        tmpk.append(tmpl)
                        spc_tuplestmp.append(tuple(tmpk))
                else:
                    tmpk = list(tmpj)
                    tmpk.append(split_listOrig)
                    spc_tuplestmp.append(tuple(tmpk))

    rxns = list(react(*spc_tuplestmp))

    return rxns

def reactPdep(*spc_tuples):
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

    reactions = map(
                react_species_pdep,
                spc_tuples)

    return itertools.chain.from_iterable(reactions)


def react_species_pdep(species_tuple):
    """
    Given a tuple of Species objects, generates all possible reactions
    from the loaded reaction families and combines degenerate reactions.
    """

    species_tuple = tuple([spc.copy(deep=True) for spc in species_tuple])

    reactions = getDB('kinetics').generate_reactions_from_families(species_tuple)

    return reactions


