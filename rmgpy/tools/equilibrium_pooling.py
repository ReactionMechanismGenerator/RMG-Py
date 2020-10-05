#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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

import logging
import os
from collections import defaultdict
from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import rdFMCS
import numpy as np
import math
import cmath
import rmgpy.constants as constants
from scipy.optimize import fsolve
from rmgpy.chemkin import save_chemkin_file

def search_priority(rmg,reaction_system):
    # print(dir(reaction_system))
    # print(dir(rmg.reaction_model.core.reactions[1]))
    sim_T=reaction_system.T.value
    core_rxn_kf=reaction_system.kf
    core_rxn_kb=reaction_system.kb
    core_concentrations=reaction_system.core_species_concentrations
    core_reactants_idx=reaction_system.reactant_indices.tolist()
    core_products_idx=reaction_system.product_indices.tolist()
    #get species list
    core_species=rmg.reaction_model.core.species
    # print(core_species[5].index)
    reactive_species=[i for i in core_species if i.reactive]
    #List of all mols and pairs
    potential_reactants = []
    for i,spc0 in enumerate(reactive_species):
        rct0_idx=core_species.index(spc0)
        potential_reactants.append((rct0_idx,))
        for spc1 in reactive_species[i:]:
            rct1_idx=core_species.index(spc1)
            potential_reactants.append((rct0_idx,rct1_idx))
    #keyed elements for all mols and pairs
    for i,pr in enumerate(potential_reactants):
        element_dict=defaultdict(int)
        multiplicity_dict=defaultdict(int)
        for spc_idx in pr:
            for elmt in core_species[spc_idx].molecule[0].get_element_count():
                element_dict[elmt]+=core_species[spc_idx].molecule[0].get_element_count()[elmt]
            multiplicity_dict[core_species[spc_idx].molecule[0].multiplicity]+=1
        potential_reactants[i]=dict([('reactants',pr),('elements',element_dict),('multiplicities',multiplicity_dict)])
        # print(element_dict)
    # print(len(potential_reactants)) #629
    #assemble multiplicity ref list
    multiplicity_ref = []
    for i in [1,2,3]:multiplicity_ref.append(defaultdict(int,[(i,1)]))
    for i in [1,2,3]:
        for j in range(i,4):
            if j==i: multiplicity_ref.append(defaultdict(int,[(i,2)]))
            else: multiplicity_ref.append(defaultdict(int,[(i,1),(j,1)]))
    # print(multiplicity_ref)
    # [{1: 1}), {2: 1}), {3: 1}), {1: 2}), {1: 1, 2: 1}), {1: 1, 3: 1}), {2: 2}), {2: 1, 3: 1}), {3: 2})]
    # allowed_multiplicity=dict([(0,[0,2,3,6]),(1,[1,4,7]),(2,[2,0,6,8]),(3,[0,3,6]),(4,[1,4,7]),(5,[2,5,6]),(6,[6,0,2,3,5,8]),(7,[7,1,4]),(8,[8,2,6])])
    # disallow two singlets to two singlets
    allowed_multiplicity=dict([(0,[0,2,3,6]),(1,[1,4,7]),(2,[2,0,6,8]),(3,[0,6]),(4,[1,4,7]),(5,[2,5,6]),(6,[6,0,2,3,5,8]),(7,[7,1,4]),(8,[8,2,6])])
    #possible connections list
    possible_connections = [] #should this be bidirectional? Currently forward rxn only.
    for i,pr0 in enumerate(potential_reactants):
        pr0_elements=pr0['elements']
        pr0_mult=pr0['multiplicities']
        for pr1 in potential_reactants[i+1:]:
            pr1_elements=pr1['elements']
            pr1_mult=pr1['multiplicities']
            if pr0_elements!=pr1_elements: continue
            if multiplicity_ref.index(pr0_mult) not in allowed_multiplicity[multiplicity_ref.index(pr1_mult)]: continue
            if len(pr0['reactants'])==1 and len(pr1['reactants'])==1:
                possible_connections.append(Connection(core_species,
                                            core_concentrations,
                                            reactant_index=pr0['reactants'][0],
                                            product_index=pr1['reactants'][0],
                                            temperature=sim_T,
                                            rct_mult=pr0_mult,
                                            prd_mult=pr1_mult))
            elif len(pr0['reactants'])==2 and len(pr1['reactants'])==1:
                possible_connections.append(Connection(core_species,
                                            core_concentrations,
                                            reactant_index=pr0['reactants'][0],
                                            coreactant_index=pr0['reactants'][1],
                                            product_index=pr1['reactants'][0],
                                            temperature=sim_T,
                                            rct_mult=pr0_mult,
                                            prd_mult=pr1_mult))
                possible_connections.append(Connection(core_species,
                                            core_concentrations,
                                            reactant_index=pr0['reactants'][1],
                                            coreactant_index=pr0['reactants'][0],
                                            product_index=pr1['reactants'][0],
                                            temperature=sim_T,
                                            rct_mult=pr0_mult,
                                            prd_mult=pr1_mult))
            elif len(pr0['reactants'])==1 and len(pr1['reactants'])==2:
                possible_connections.append(Connection(core_species,
                                            core_concentrations,
                                            reactant_index=pr0['reactants'][0],
                                            coproduct_index=pr1['reactants'][1],
                                            product_index=pr1['reactants'][0],
                                            temperature=sim_T,
                                            rct_mult=pr0_mult,
                                            prd_mult=pr1_mult))
                possible_connections.append(Connection(core_species,
                                            core_concentrations,
                                            reactant_index=pr0['reactants'][0],
                                            coproduct_index=pr1['reactants'][0],
                                            product_index=pr1['reactants'][1],
                                            temperature=sim_T,
                                            rct_mult=pr0_mult,
                                            prd_mult=pr1_mult))
            elif len(pr0['reactants'])==2 and len(pr1['reactants'])==2:
                possible_connections.append(Connection(core_species,
                                            core_concentrations,
                                            reactant_index=pr0['reactants'][0],
                                            coreactant_index=pr0['reactants'][1],
                                            coproduct_index=pr1['reactants'][1],
                                            product_index=pr1['reactants'][0],
                                            temperature=sim_T,
                                            rct_mult=pr0_mult,
                                            prd_mult=pr1_mult))
                possible_connections.append(Connection(core_species,
                                            core_concentrations,
                                            reactant_index=pr0['reactants'][1],
                                            coreactant_index=pr0['reactants'][0],
                                            coproduct_index=pr1['reactants'][0],
                                            product_index=pr1['reactants'][1],
                                            temperature=sim_T,
                                            rct_mult=pr0_mult,
                                            prd_mult=pr1_mult))
                possible_connections.append(Connection(core_species,
                                            core_concentrations,
                                            reactant_index=pr0['reactants'][0],
                                            coreactant_index=pr0['reactants'][1],
                                            coproduct_index=pr1['reactants'][0],
                                            product_index=pr1['reactants'][1],
                                            temperature=sim_T,
                                            rct_mult=pr0_mult,
                                            prd_mult=pr1_mult))
                possible_connections.append(Connection(core_species,
                                            core_concentrations,
                                            reactant_index=pr0['reactants'][1],
                                            coreactant_index=pr0['reactants'][0],
                                            coproduct_index=pr1['reactants'][1],
                                            product_index=pr1['reactants'][0],
                                            temperature=sim_T,
                                            rct_mult=pr0_mult,
                                            prd_mult=pr1_mult)) 
            else: raise Exception('Invalid equilibrium pooling handling')
    # print(len(possible_connections)) #4656 vs 750k vs mult 335k vs elem 5010
    #direct connections list
    for i,ids in enumerate(core_reactants_idx):
        core_reactants_idx[i] = [j for j in ids if j != -1]
    for i,ids in enumerate(core_products_idx):
        core_products_idx[i] = [j for j in ids if j != -1]
    direct_connections_list=[]
    for i in range(len(core_reactants_idx)):
        r_ids=core_reactants_idx[i]
        p_ids=core_products_idx[i]
        for j in range(len(r_ids)):
            for k in range(len(p_ids)):
                if len(r_ids)==1 and len(p_ids)==1: #commented to make unidirectional
                    direct_connections_list.append(Connection(core_species,core_concentrations,
                                                            reactant_index=r_ids[j],
                                                            product_index=p_ids[k],
                                                            kf=core_rxn_kf[i],
                                                            kb=core_rxn_kb[i],
                                                            temperature=sim_T))
                    # direct_connections_list.append(Connection(core_species,core_concentrations,
                    #                                         reactant_index=p_ids[k],
                    #                                         product_index=r_ids[j],
                    #                                         kf=core_rxn_kb[i],
                    #                                         kb=core_rxn_kf[i],
                    #                                         temperature=sim_T))
                if len(r_ids)==2 and len(p_ids)==1:
                    direct_connections_list.append(Connection(core_species,core_concentrations,
                                                            reactant_index=r_ids[j],
                                                            coreactant_index=r_ids[1-j],
                                                            product_index=p_ids[k],
                                                            kf=core_rxn_kf[i],
                                                            kb=core_rxn_kb[i],
                                                            temperature=sim_T))
                    # direct_connections_list.append(Connection(core_species,core_concentrations,
                    #                                         reactant_index=p_ids[k],
                    #                                         product_index=r_ids[j],
                    #                                         coproduct_index=r_ids[1-j],
                    #                                         kf=core_rxn_kb[i],
                    #                                         kb=core_rxn_kf[i],
                    #                                         temperature=sim_T))
                if len(r_ids)==1 and len(p_ids)==2:
                    direct_connections_list.append(Connection(core_species,core_concentrations,
                                                            reactant_index=r_ids[j],
                                                            product_index=p_ids[k],
                                                            coproduct_index=p_ids[1-k],
                                                            kf=core_rxn_kf[i],
                                                            kb=core_rxn_kb[i],
                                                            temperature=sim_T))
                    # direct_connections_list.append(Connection(core_species,core_concentrations,
                    #                                         reactant_index=p_ids[k],
                    #                                         coreactant_index=p_ids[1-k],
                    #                                         product_index=r_ids[j],
                    #                                         kf=core_rxn_kb[i],
                    #                                         kb=core_rxn_kf[i],
                    #                                         temperature=sim_T))
                if len(r_ids)==2 and len(p_ids)==2:
                    direct_connections_list.append(Connection(core_species,core_concentrations,
                                                            reactant_index=r_ids[j],
                                                            coreactant_index=r_ids[1-j],
                                                            product_index=p_ids[k],
                                                            coproduct_index=p_ids[1-k],
                                                            kf=core_rxn_kf[i],
                                                            kb=core_rxn_kb[i],
                                                            temperature=sim_T))
                    # direct_connections_list.append(Connection(core_species,core_concentrations,
                    #                                         reactant_index=p_ids[k],
                    #                                         coreactant_index=p_ids[1-k],
                    #                                         product_index=r_ids[j],
                    #                                         coproduct_index=r_ids[1-j],
                    #                                         kf=core_rxn_kb[i],
                    #                                         kb=core_rxn_kf[i],
                    #                                         temperature=sim_T))
#remove already covered connections
    removed_connections_idx=set()
    hash_set=[]
    for conn0 in direct_connections_list:
        conn0_hash=conn0.hash
        for e,conn1 in enumerate(possible_connections):
            conn1_hash=conn1.hash
            if conn0_hash==conn1_hash:
                removed_connections_idx.add(e)
    for e,conn in enumerate(possible_connections):
        if conn.reactant_conc!=conn.min_reactant_conc or conn.product_conc!=conn.min_product_conc:
            removed_connections_idx.add(e)
        elif conn.hash in hash_set:
            removed_connections_idx.add(e)
        else:
            hash_set.append(conn.hash)
        
    possible_connections=[i for e,i in enumerate(possible_connections) if e not in removed_connections_idx]

    for conn in direct_connections_list:
        conn.calculate_equil_flux()
        conn.calculate_time_constant()
    pooling_times=pooling_time_cycle(direct_connections_list,core_species)
    for conn in possible_connections:
        conn.ratchet_pooling_time(pooling_times[conn.reactant,conn.product])
        # print(conn.pooling_time,conn.largest_substructure,conn.H_rxn,conn.reactant,conn.coreactant,conn.product,conn.coproduct)
    # chemkin_path=os.path.join(rmg.output_directory,'pooling_chem.inp')
    # save_chemkin_file(path=chemkin_path,species=core_species,reactions=rmg.reaction_model.core.reactions,check_for_duplicates=False)
    #classify connections by pooling time, delta H, subgroup similarity, unimolecular,distance from equil,multiplicity,same element
    ranked_conn_list=sort_possible_connections(possible_connections)
    # return ranked_conn_list

    ranked_core_reaction_indices=rank_core_reactions(rmg,pooling_times,core_species,core_concentrations,sim_T)
    return ranked_core_reaction_indices

    # for conn in ranked_conn_list:
    #     print(conn.pooling_time,conn.largest_substructure,conn.H_rxn,conn.reactant,conn.coreactant,conn.product,conn.coproduct)

def equilibrium_pooling(direct_connections_list,core_species,t_characteristic): #,previous_forward_pass=None):
    fast_connections=[]
    # strong_connections=[]
    for rct_idx in range(len(core_species)):
        fast_connections.append(set())
        # strong_connections.append(set())
    for conn in direct_connections_list:
        if conn.time_constant is None:
            continue
        elif conn.product_conc!=conn.min_product_conc or conn.reactant_conc!=conn.min_reactant_conc:
            continue
        elif conn.time_constant<t_characteristic:
            fast_connections[conn.reactant].add(conn.product)
            fast_connections[conn.product].add(conn.reactant)
            # if conn.min_reactant_conc>0.1*conn.min_product_conc:
                # strong_connections[conn.reactant].add(conn.product)
        # print("connections",conn.time_constant,conn.reactant,conn.product)
    # if previous_forward_pass is not None:
        # forward_pass=previous_forward_pass
    # else:
    # print('fast_connections')
    # for i in fast_connections: print(i)
    forward_pass=[]
    # strong_pass=[]
    for rct_idx in range(len(core_species)):
        forward_pass.append(set([rct_idx]))
        # strong_pass.append(set([rct_idx]))
    for iter in range(len(core_species)):
        for rct_idx in range(len(core_species)):
            for prd_idx in fast_connections[rct_idx]:
                forward_pass[rct_idx].update(forward_pass[prd_idx])
            # for prd_idx in strong_connections[rct_idx]:
                # strong_pass[rct_idx].update(strong_pass[prd_idx])
    # print('forward_pass',forward_pass[20],t_characteristic)
    # for e,i in enumerate(forward_pass):
    #     if 20 in i:
    #         print('forward pass to',e)
    # print('strong_pass',strong_pass[20])
    # for i in forward_pass: print(i)
    pool_groups=[]
    for rct_idx in range(len(core_species)):
        pool_groups.append(set())
    # broken_pools=set()
    for rct_idx in range(len(core_species)):
        for prd_idx in forward_pass[rct_idx]:
            # if rct_idx in forward_pass[prd_idx]:
                pool_groups[rct_idx].add(prd_idx)
            # else: 
            #     broken_pools.add(rct_idx)
            #     if prd_idx in strong_pass[rct_idx]:
            #         broken_pools.add(prd_idx)
    # print("broken_pools",broken_pools)
    # print('pools before')
    # for i,group in enumerate(pool_groups):print(i,group)
    # for rct_idx in broken_pools:
    #     for prd_idx in pool_groups[rct_idx]:
    #         pool_groups[prd_idx]=set()
    # print('pool after')
    # for i,group in enumerate(pool_groups):print(i,group)
    return pool_groups

def pooling_time_cycle(direct_connections_list,core_species):
    time_constants=[]
    for conn in direct_connections_list:
        if conn.time_constant is not None:
            time_constants.append(conn.time_constant)
    time_constants.sort()
    # print(time_constants)
    pooling_times=np.empty([len(core_species),len(core_species)])
    pooling_times[:]=np.inf
    # log10_minimum=math.floor(math.log10(time_constants[0]))
    # print(time_constants)
    log10_minimum=math.floor(math.log10([i for i in time_constants if i>0][0]))
    log10_maximum=math.ceil(math.log10(time_constants[-1]))+1
    # log10_minimum=math.floor(math.log2(time_constants[0]))
    # log10_maximum=math.ceil(math.log2(time_constants[-1]))+1
    # pool_groups=equilibrium_pooling(direct_connections_list,core_species,10**log10_maximum)
    for log10_time in list(range(log10_minimum,log10_maximum))[::-1]:
    # for log10_time in [i/100 for i in list(range(log10_minimum*100,log10_maximum*100))[::-1]]:
        # print(log10_time)
        pool_groups=equilibrium_pooling(direct_connections_list,core_species,10**log10_time)
        # pool_groups=equilibrium_pooling(direct_connections_list,core_species,2**log10_time)
        # print(pool_groups)
        for i,grp in enumerate(pool_groups):
            for j in grp:
                pooling_times[i,j]=log10_time
    # for i in pooling_times: print(i) 
    # print(pool_groups)
    # for i,spc in enumerate(core_species): print(i,spc)
    # for i in pool_groups[-1]:
        # print(core_species[i].label)
    # equilibrium_pooling(direct_connections_list,core_species,10**-15)
    # print(log10_minimum,log10_maximum)
    return pooling_times

def sort_possible_connections(possible_connections):
    for conn in possible_connections:
        conn.calculate_equil_flux()
    ranked_conn_list=possible_connections
    ranked_conn_list.sort(key=connection_sort_key,reverse=True)
    return ranked_conn_list

def connection_sort_key(conn): #reversed
    if conn.reactants_multiplicity=={1: 2} and conn.products_multiplicity=={1: 2}:
        mult_score=0.1
        print('mult')
    else:
        mult_score=1
    if conn.min_reactant_conc>conn.min_product_conc:
        H_score=-conn.H_rxn/1000
    else:
        H_score=conn.H_rxn/1000
    if abs(H_score)<10:
        H_score=1
    elif H_score<0:
        H_score=math.log10(-H_score)
    else:
        H_score=1/math.log10(H_score)
    if conn.species_count[1]==1:
        if conn.reactant==conn.product or conn.coreactant==conn.product:
            same_species_score=0.1
        else:
            same_species_score=1
    else:
        if conn.reactant in [conn.product,conn.coproduct] or conn.coreactant in [conn.product,conn.coproduct]:
            same_species_score=0.1
        else:
            same_species_score=1
    substructure_score=conn.largest_substructure+1
    composite=H_score*mult_score*same_species_score*substructure_score
    if None in (conn.pooling_time,conn.largest_substructure,composite): print(conn.pooling_time,conn.largest_substructure,composite)
    return (conn.pooling_time,composite)
    #classify connections by pooling time, delta H, subgroup similarity, unimolecular,distance from equil,multiplicity,same element

def rank_core_reactions(rmg,pooling_times,core_species,core_concentrations,temperature):
    core_reactions=rmg.reaction_model.core.reactions
    core_species_labels=[i.index for i in core_species]
    core_reaction_hash_list=[]
    # print(core_species)
    core_connections=[]
    for rxn in core_reactions:
        rxn_species=[[],[]]
        reactants_mult=defaultdict(int)
        products_mult=defaultdict(int)
        for rct in rxn.reactants:
            rxn_species[0].append(core_species_labels.index(rct.index))
        for prd in rxn.products:
            # print(prd.index)
            rxn_species[1].append(core_species_labels.index(prd.index))
        for rct in rxn_species[0]:
            reactants_mult[core_species[rct].molecule[0].multiplicity]+=1
        for prd in rxn_species[1]:
            # print(prd)
            products_mult[core_species[prd].molecule[0].multiplicity]+=1
        # reaction_elements.append((rxn_species,reactants_mult,products_mult))
        if len(rxn_species[0])==2:
            # print(rxn_species[0])
            # print(rxn_species[0][1])
            # print(rxn.products)
            if core_concentrations[rxn_species[0][0]]<core_concentrations[rxn_species[0][1]]:
                rct_idx=rxn_species[0][0]
                corct_idx=rxn_species[0][1]
            else:
                rct_idx=rxn_species[0][1]
                corct_idx=rxn_species[0][0]
        else:
            rct_idx=rxn_species[0][0]
            corct_idx=None
        if len(rxn_species[1])==2:
            if core_concentrations[rxn_species[1][0]]<core_concentrations[rxn_species[1][1]]:
                prd_idx=rxn_species[1][0]
                coprd_idx=rxn_species[1][1]
            else:
                prd_idx=rxn_species[1][1]
                coprd_idx=rxn_species[1][0]
        else:
            prd_idx=rxn_species[1][0]
            coprd_idx=None
        core_connections.append(Connection(core_species,
                core_concentrations,
                reactant_index=rct_idx,
                coreactant_index=corct_idx,
                coproduct_index=coprd_idx,
                product_index=prd_idx,
                temperature=temperature,
                rct_mult=reactants_mult,
                prd_mult=products_mult))

    for conn in core_connections:
        conn.ratchet_pooling_time(pooling_times[conn.reactant,conn.product])
        core_reaction_hash_list.append(conn.hash)
    ranked_conn_list=sort_possible_connections(core_connections)
    ranked_reaction_indices=[]
    for conn in ranked_conn_list:
        idx=core_reaction_hash_list.index(conn.hash)
        ranked_reaction_indices.append(idx)
        core_reaction_hash_list[idx]=None
    # ranked_reaction_indices=[core_reaction_hash_list.index(i.hash) for i in ranked_conn_list]
    # for conn in ranked_conn_list:
    #     ranked_reaction_indices.append([e for e,i in enumerate(core_reaction_hash_list) if i==conn.hash])
    # print(ranked_reaction_indices)
    return ranked_reaction_indices

    # for i,pr in enumerate(potential_reactants):
    # element_dict=defaultdict(int)
    # multiplicity_dict=defaultdict(int)
    # for spc_idx in pr:
    #     for elmt in core_species[spc_idx].molecule[0].get_element_count():
    #         element_dict[elmt]+=core_species[spc_idx].molecule[0].get_element_count()[elmt]
    #     multiplicity_dict[core_species[spc_idx].molecule[0].multiplicity]+=1
    # potential_reactants[i]=dict([('reactants',pr),('elements',element_dict),('multiplicities',multiplicity_dict)])


#################################################################################


class Connection:
    """
    ======================= ================ =========================================
    Attribute               Type             Description
    ======================= ================ =========================================
    `reactant_index`
    `reactant_conc`
    `coreactant_index`
    `coreactant_conc`
    `product_index`
    `product_conc`
    `coproduct_index`
    `coproduct_conc`       ``float``
    `kf`                  ``float``
    `kb`
    `temperature`           ``float``
    `pfo_k`              ``float``
    `time_constant`         ``float``
    `H_rxn`                 ``float``
    `largest_substructure`  ``int``
    `pooling_time`          ``float``
    `hash`                  ``tuple``
    `min_reactant_conc`
    `delta_reactant_conc`
    `min_product_conc`
    `delta_product_conc`
    `equil_flux`
   ======================== ================= =========================================

    """

    def __init__(self,core_species,core_concentrations,reactant_index,product_index,coreactant_index=None,coproduct_index=None,temperature=298,kf=None,kb=None,rct_mult=None,prd_mult=None):
        self.reactant=reactant_index
        self.reactant_conc=core_concentrations[reactant_index]
        self.product=product_index
        self.product_conc=core_concentrations[product_index]
        if coreactant_index is not None and coproduct_index is not None:
            self.species_count=(2,2)
        elif coproduct_index is not None:
            self.species_count=(1,2)
        elif coreactant_index is not None:
            self.species_count=(2,1)
        else:
            self.species_count=(1,1)
        if self.species_count[0]==2:
            self.coreactant=coreactant_index
            self.coreactant_conc=core_concentrations[coreactant_index]
            self.min_reactant_conc=min([self.coreactant_conc,self.reactant_conc])
            self.delta_reactant_conc=max([self.coreactant_conc,self.reactant_conc])-self.min_reactant_conc
        else:
            self.coreactant=None
            self.coreactant_conc=None
            self.min_reactant_conc=self.reactant_conc
            self.delta_reactant_conc=0
        if self.species_count[1]==2:
            self.coproduct=coproduct_index
            self.coproduct_conc=core_concentrations[coproduct_index]
            self.min_product_conc=min([self.coproduct_conc,self.product_conc])
            self.delta_product_conc=max([self.coproduct_conc,self.product_conc])-self.min_product_conc
        else:
            self.coproduct=None
            self.coproduct_conc=None
            self.min_product_conc=self.product_conc
            self.delta_product_conc=0
        self.kf=kf
        self.kb=kb
        self.T=temperature
        self.pooling_time=None
        self.hash=set(((reactant_index,coreactant_index if coreactant_index is not None else None),
            (product_index,coproduct_index if coproduct_index is not None else None)))
        self.reactants_multiplicity=rct_mult
        self.products_multiplicity=prd_mult
        #thermo rxn calcs
        self.H_rxn=core_species[self.product].get_enthalpy(temperature)-core_species[self.reactant].get_enthalpy(temperature)
        self.G_rxn=core_species[self.product].get_free_energy(temperature)-core_species[self.reactant].get_free_energy(temperature)
        if self.species_count[1]==2:
            self.H_rxn+=core_species[self.coproduct].get_enthalpy(temperature)
            self.G_rxn+=core_species[self.coproduct].get_free_energy(temperature)
        if self.species_count[0]==2:
            self.H_rxn-=core_species[self.coreactant].get_enthalpy(temperature)
            self.G_rxn-=core_species[self.coreactant].get_free_energy(temperature)
        #Keq calculation
        if self.kf is not None:
            if self.kb !=0:
                self.Keq=self.kf/self.kb
            else:
                self.Keq=float('inf')
        else:
            self.Keq=math.exp(-self.G_rxn/self.T/constants.R)*(100000/constants.R/self.T)**(self.species_count[1]-self.species_count[0])
        # self.calculate_equil_flux()
        # self.calculate_time_constant()

        lg=RDLogger.logger()
        lg.setLevel(RDLogger.ERROR)
        #substructure match
        if self.species_count==(2,2):
            rdkit_rct1=Chem.MolFromSmiles(core_species[reactant_index].molecule[0].smiles)
            rdkit_rct2=Chem.MolFromSmiles(core_species[coreactant_index].molecule[0].smiles)
            rdkit_prd1=Chem.MolFromSmiles(core_species[product_index].molecule[0].smiles)
            rdkit_prd2=Chem.MolFromSmiles(core_species[coproduct_index].molecule[0].smiles)
            rdkit_rct=Chem.CombineMols(rdkit_rct1,rdkit_rct2)
            rdkit_prd=Chem.CombineMols(rdkit_prd1,rdkit_prd2)
            self.largest_substructure=min(rdFMCS.FindMCS([rdkit_rct1,rdkit_prd]).numBonds+rdFMCS.FindMCS([rdkit_rct2,rdkit_prd]).numBonds,rdFMCS.FindMCS([rdkit_prd1,rdkit_rct]).numBonds+rdFMCS.FindMCS([rdkit_prd2,rdkit_rct]).numBonds)
        if self.species_count==(1,2):
            rdkit_rct=Chem.MolFromSmiles(core_species[reactant_index].molecule[0].smiles)
            rdkit_prd1=Chem.MolFromSmiles(core_species[product_index].molecule[0].smiles)
            rdkit_prd2=Chem.MolFromSmiles(core_species[coproduct_index].molecule[0].smiles)
            rdkit_prd=Chem.CombineMols(rdkit_prd1,rdkit_prd2)
            self.largest_substructure=min(rdFMCS.FindMCS([rdkit_rct,rdkit_prd]).numBonds,rdFMCS.FindMCS([rdkit_prd1,rdkit_rct]).numBonds+rdFMCS.FindMCS([rdkit_prd2,rdkit_rct]).numBonds)
        if self.species_count==(2,1):
            rdkit_rct1=Chem.MolFromSmiles(core_species[reactant_index].molecule[0].smiles)
            rdkit_rct2=Chem.MolFromSmiles(core_species[coreactant_index].molecule[0].smiles)
            rdkit_prd=Chem.MolFromSmiles(core_species[product_index].molecule[0].smiles)
            rdkit_rct=Chem.CombineMols(rdkit_rct1,rdkit_rct2)
            self.largest_substructure=min(rdFMCS.FindMCS([rdkit_rct1,rdkit_prd]).numBonds+rdFMCS.FindMCS([rdkit_rct2,rdkit_prd]).numBonds,rdFMCS.FindMCS([rdkit_prd,rdkit_rct]).numBonds)
        if self.species_count==(1,1):
            rdkit_rct=Chem.MolFromSmiles(core_species[reactant_index].molecule[0].smiles)
            rdkit_prd=Chem.MolFromSmiles(core_species[product_index].molecule[0].smiles)
            self.largest_substructure=rdFMCS.FindMCS([rdkit_rct,rdkit_prd]).numBonds

    
    def ratchet_pooling_time(self,new_time):
        if self.pooling_time is None:
            self.pooling_time=new_time
        elif self.pooling_time>new_time:
            self.pooling_time=new_time

    def is_same_connection(self,other_connection): #there's a bug here somehow with how it's filtering. self.reactant is a looser filter than self.reactant.label
        # if self.reactant == other_connection.reactant:
        #     if (self.product,self.coreactant,self.coproduct) == (other_connection.product,other_connection.coreactant,other_connection.coproduct):return True
        # if self.reactant == other_connection.product:
        #     if (self.product,self.coreactant,self.coproduct) == (other_connection.reactant,other_connection.coproduct,other_connection.coreactant): return True
        # return False
        return self.hash==other_connection.hash
    
    def calculate_equil_flux(self):
        if max(self.min_product_conc,self.min_reactant_conc)<1e-18:
            self.equil_flux=0
            return
        elif self.species_count==(2,2):
            if self.delta_reactant_conc!=0 and (self.product_conc+self.min_reactant_conc)*(self.coproduct_conc+self.min_reactant_conc)/(self.Keq)/self.delta_reactant_conc<1e-18:
                # print((self.product_conc+self.min_reactant_conc)*(self.coproduct_conc+self.min_reactant_conc)/(self.Keq)/self.delta_reactant_conc,'small remainder')
                self.equil_flux=0
                return
            if self.min_reactant_conc==0:
                inst_K=float('inf')
            else:
                inst_K=self.product_conc*self.coproduct_conc/self.reactant_conc/self.coreactant_conc
            self.flux_guess=((2*self.Keq*self.min_reactant_conc+self.Keq*self.delta_reactant_conc+2*self.min_product_conc+self.delta_product_conc-math.sqrt(self.Keq**2*self.delta_reactant_conc**2+4*self.Keq*self.min_reactant_conc**2+8*self.Keq*self.min_reactant_conc*self.min_product_conc+4*self.Keq*self.min_reactant_conc*self.delta_product_conc+4*self.Keq*self.min_reactant_conc*self.delta_reactant_conc+4*self.Keq*self.min_product_conc**2+4*self.Keq*self.min_product_conc*self.delta_product_conc+4*self.Keq*self.min_product_conc*self.delta_reactant_conc+2*self.Keq*self.delta_reactant_conc*self.delta_product_conc+self.delta_product_conc**2))/(2*(self.Keq-1)))
            # print(self.Keq,self.reactant_conc,self.coreactant_conc,self.product_conc,self.coproduct_conc)
            # print(self.flux_guess)
            if not self.pos_flux_check(self.flux_guess) or self.flux_guess==0:
                if inst_K/self.Keq<1:
                    self.flux_guess=self.min_reactant_conc/2
                else:
                    self.flux_guess=-self.min_product_conc/2
            # print(self.flux_guess)
            func= lambda x: (self.product_conc+x)*(self.coproduct_conc+x)-(self.reactant_conc-x)*(self.coreactant_conc-x)*self.Keq
            self.equil_flux=fsolve(func,self.flux_guess)[0]
            if not self.flux_K_check():
                # if inst_K/self.Keq<1 and abs((self.min_reactant_conc-self.equil_flux)/self.min_reactant_conc)<1e-5:
                #     func= lambda x: (self.product_conc+self.min_reactant_conc-x)*(self.coproduct_conc+self.min_reactant_conc-x)-x*(self.delta_reactant_conc+x)*self.Keq
                #     self.equil_flux=self.min_reactant_conc-fsolve(func,(self.product_conc+self.min_reactant_conc)*(self.coproduct_conc+self.min_reactant_conc)/self.delta_reactant_conc/self.Keq)[0]
                # if inst_K/self.Keq>1 and abs((self.min_product_conc+self.equil_flux)/self.min_product_conc)<1e-5:
                #     func= lambda x: (self.delta_product_conc-x)*x-(self.reactant_conc+self.min_product_conc+x)*(self.coreactant_conc+self.min_product_conc+x)*self.Keq
                #     self.equil_flux=-self.min_product_conc+fsolve(func,self.Keq*(self.reactant_conc+self.min_product_conc)*(self.coreactant_conc+self.min_product_conc)/self.delta_product_conc)[0]
                #     print('tried it')
                # if not self.flux_K_check():
                #     print(abs((self.min_product_conc+self.equil_flux)/self.min_product_conc))
                #     print(inst_K/self.Keq)
                #     print(-self.min_product_conc+self.Keq*(self.reactant_conc+self.min_product_conc)*(self.coreactant_conc+self.min_product_conc)/self.delta_product_conc)
                    # print(self.Keq,self.reactant_conc,self.coreactant_conc,self.product_conc)
                    # print(self.flux_guess,self.equil_flux)
                    raise Exception('bad equil')
            return
        # elif self.species_count==(2,1) and self.product_conc<1e-10:
        #     inst_K=self.product_conc/self.reactant_conc/self.coreactant_conc
        #     if self.delta_reactant_conc!=0 and inst_K/self.Keq<1 and (self.product_conc+self.min_reactant_conc)/(self.delta_reactant_conc)/self.Keq/self.min_reactant_conc<1e-5:
        #         func= lambda x: math.log((self.min_reactant_conc-x)/x/(self.delta_reactant_conc+x)/self.Keq)
        #         self.equil_flux=self.min_reactant_conc-fsolve(func,(self.min_reactant_conc)/self.delta_reactant_conc/self.Keq)[0]
        #         if not self.flux_K_check():
        #             raise Exception('bad equil')
        #         return
        #     self.flux_guess=(2*self.Keq*self.min_reactant_conc+self.Keq*self.delta_reactant_conc-math.sqrt(self.Keq**2*self.delta_reactant_conc**2+4*self.Keq*self.min_reactant_conc+4*self.Keq*self.product_conc+2*self.Keq*self.delta_reactant_conc+1)+1)/(2*self.Keq)
        #     if not self.pos_flux_check(self.flux_guess) or self.flux_guess==0:
        #         self.flux_guess=self.Keq*(self.min_reactant_conc)*(self.min_reactant_conc+self.delta_reactant_conc)/(self.Keq+1)
        #     print(self.Keq,self.reactant_conc,self.coreactant_conc,self.product_conc)
        #     print(self.flux_guess)
        #     func= lambda x: math.log((self.product_conc+x)/(self.reactant_conc-x)/(self.coreactant_conc-x)/self.Keq)
        #     self.equil_flux=fsolve(func,self.flux_guess)[0]
        #     if not self.flux_K_check():
        #         raise Exception('bad equil')
        #     return
        elif self.species_count==(2,1):
            if self.min_reactant_conc==0:
                inst_K=float('inf')
            else:
                inst_K=self.product_conc/self.reactant_conc/self.coreactant_conc
            self.flux_guess=(2*self.Keq*self.min_reactant_conc+self.Keq*self.delta_reactant_conc-math.sqrt(self.Keq**2*self.delta_reactant_conc**2+4*self.Keq*self.min_reactant_conc+4*self.Keq*self.product_conc+2*self.Keq*self.delta_reactant_conc+1)+1)/(2*self.Keq)
            if not self.pos_flux_check(self.flux_guess) or self.flux_guess==0:
                if inst_K/self.Keq<1:
                    self.flux_guess=self.min_reactant_conc/2
                else:
                    self.flux_guess=-self.min_product_conc/2
            # print(self.Keq,self.reactant_conc,self.coreactant_conc,self.product_conc)
            # print(self.flux_guess)
            func= lambda x: (self.product_conc+x)-(self.reactant_conc-x)*(self.coreactant_conc-x)*self.Keq
            self.equil_flux=fsolve(func,self.flux_guess)[0]
            if not self.flux_K_check():
                if self.delta_reactant_conc!=0 and inst_K/self.Keq<1 and (self.product_conc+self.min_reactant_conc)/(self.delta_reactant_conc)/self.Keq/self.min_reactant_conc<1e-5:
                    func= lambda x: (self.product_conc+self.min_reactant_conc-x)-x*(self.delta_reactant_conc+x)*self.Keq
                    self.equil_flux=self.min_reactant_conc-fsolve(func,(self.product_conc+self.min_reactant_conc)*(self.coproduct_conc+self.min_reactant_conc)/self.delta_reactant_conc/self.Keq)[0]
                    if not self.flux_K_check():
                        raise Exception('bad equil')
                    return
                raise Exception('bad equil')
            return
        elif self.species_count==(1,2):
            if self.min_reactant_conc==0:
                inst_K=float('inf')
            else:
                inst_K=self.product_conc*self.coproduct_conc/self.reactant_conc
            self.flux_guess=-self.Keq/2-self.min_product_conc-self.delta_product_conc/2+math.sqrt(self.Keq**2+4*self.Keq*self.reactant_conc+4*self.Keq*self.min_product_conc+2*self.Keq*self.delta_product_conc+self.delta_product_conc**2)/2
            if not self.pos_flux_check(self.flux_guess) or self.flux_guess==0:
                if inst_K/self.Keq<1:
                    self.flux_guess=self.min_reactant_conc/2
                else:
                    self.flux_guess=-self.min_product_conc/2
            # if self.Keq<1:
            #     func= lambda x: (self.product_conc+x)*(self.coproduct_conc+x)/(self.reactant_conc-x)-self.Keq
            # else:
            #     func= lambda x: (self.reactant_conc-x)/(self.product_conc+x)/(self.coproduct_conc+x)-1/self.Keq
            func= lambda x: (self.product_conc+x)*(self.coproduct_conc+x)-(self.reactant_conc-x)*self.Keq
            self.equil_flux=fsolve(func,self.flux_guess)[0]
            if not self.flux_K_check():
                raise Exception('bad equil')
            return
        else:
            if self.min_reactant_conc==0:
                inst_K=float('inf')
            else:
                inst_K=self.product_conc/self.reactant_conc
            self.flux_guess=(self.Keq*self.reactant_conc-self.product_conc)/(self.Keq+1)
            if not self.pos_flux_check(self.flux_guess) or self.flux_guess==0:
                if inst_K/self.Keq<1:
                    self.flux_guess=self.min_reactant_conc/2
                else:
                    self.flux_guess=-self.product_conc/2
            # func= lambda x: (self.product_conc+x)/(self.reactant_conc-x)-self.Keq
            func= lambda x: (self.product_conc+x)-(self.reactant_conc-x)*self.Keq
            self.equil_flux=fsolve(func,self.flux_guess)[0]
            if not self.flux_K_check():
                raise Exception('bad equil')
            return

    def calculate_time_constant(self):
        if self.kf is not None:
            if self.min_reactant_conc==0 and self.min_product_conc==0:
                self.time_constant=None
        #pseudo first order time constant
            # self.pfo_k=self.kf*self.coreactant_conc if self.species_count[0]==2 else self.kf
            # if self.pfo_k!=0:
            #     self.time_constant=1/self.pfo_k
            # else:
            #     self.time_constant = None
        #tangent rate to equilibrium flux
            # if self.species_count[0]==1:
            #     self.time_constant=self.equil_flux/(self.kf*self.reactant_conc)
            # else:
            #     self.time_constant=self.equil_flux/(self.kf*self.reactant_conc*self.coreactant_conc)
        #time to 0.9 equilibrium flux
            if abs(self.equil_flux)<1e-16 or max(self.min_reactant_conc,self.min_product_conc)<1e-18 or self.equil_flux<min(self.min_reactant_conc,self.min_product_conc)*1e-5:
            # if self.equil_flux<min(self.min_reactant_conc,self.min_product_conc)*1e-5:
                self.time_constant=0 #sufficiently equilibrated that equil_flux is in the range of numerical precision
            elif self.species_count==(1,1):
                # self.time_constant=math.log(-self.reactant_conc*self.kf+self.product_conc*self.kb)/(self.kf+self.kb)-math.log(-self.reactant_conc*self.kf+self.product_conc*self.kb+(self.kf+self.kb)*self.equil_flux*0.9)/(self.kf+self.kb)
                self.time_constant=math.log((-self.reactant_conc*self.kf+self.product_conc*self.kb)/(-self.reactant_conc*self.kf+self.product_conc*self.kb+(self.kf+self.kb)*self.equil_flux*0.9))/(self.kf+self.kb)
            elif self.species_count==(2,1):
                a1=math.sqrt(1/(4*self.min_reactant_conc*self.kf*self.kb + 4*self.product_conc*self.kf*self.kb + self.delta_reactant_conc**2*self.kf**2 + 2*self.delta_reactant_conc*self.kf*self.kb + self.kb**2))
                a2=(-4*self.min_reactant_conc*self.kf*self.kb*a1 - 2*self.min_reactant_conc*self.kf - 4*self.product_conc*self.kf*self.kb*a1 - self.delta_reactant_conc**2*self.kf**2*a1 - 2*self.delta_reactant_conc*self.kf*self.kb*a1 - self.delta_reactant_conc*self.kf - self.kb**2*a1 - self.kb)/(2*self.kf)
                a3=(4*self.min_reactant_conc*self.kf*self.kb*a1 - 2*self.min_reactant_conc*self.kf + 4*self.product_conc*self.kf*self.kb*a1 + self.delta_reactant_conc**2*self.kf**2*a1 + 2*self.delta_reactant_conc*self.kf*self.kb*a1 - self.delta_reactant_conc*self.kf + self.kb**2*a1 - self.kb)/(2*self.kf)
                # print(a1,a2,a3)
                if a3==0:
                    print(a1,a2,a3)
                    print(self.kf,self.kb)
                    print(self.reactant_conc,self.coreactant_conc,self.product_conc,self.coproduct_conc)
                    print(self.equil_flux,self.flux_guess,self.species_count)
                    print(self.Keq)
                # C1=-a1*math.log(a2)+a1*math.log(a3)
                self.time_constant=a1*(math.log((self.equil_flux*0.9+a2)*a3/a2/(self.equil_flux*0.9+a3)))
            elif self.species_count==(1,2):
                a1=math.sqrt(1/(4*self.min_product_conc*self.kb*self.kf + 4*self.min_reactant_conc*self.kb*self.kf + self.delta_product_conc**2*self.kb**2 + 2*self.delta_product_conc*self.kb*self.kf + self.kf**2))
                a2=(-4*self.min_product_conc*self.kb*self.kf*a1 - 2*self.min_product_conc*self.kb - 4*self.min_reactant_conc*self.kb*self.kf*a1 - self.delta_product_conc**2*self.kb**2*a1 - 2*self.delta_product_conc*self.kb*self.kf*a1 - self.delta_product_conc*self.kb - self.kf**2*a1 - self.kf)/(2*self.kb)
                a3=(4*self.min_product_conc*self.kb*self.kf*a1 - 2*self.min_product_conc*self.kb + 4*self.min_reactant_conc*self.kb*self.kf*a1 + self.delta_product_conc**2*self.kb**2*a1 + 2*self.delta_product_conc*self.kb*self.kf*a1 - self.delta_product_conc*self.kb + self.kf**2*a1 - self.kf)/(2*self.kb)
                # if a3==0:
                    # print(a1,a2,a3)
                    # print(self.kf,self.kb)
                    # print(self.reactant_conc,self.coreactant_conc,self.product_conc,self.coproduct_conc)
                    # print(self.equil_flux,self.species_count)                # C1=-a1*math.log(a2)+a1*math.log(a3)
                # self.time_constant=a1*math.log(-self.equil_flux*0.9+a2)-a1*math.log(-self.equil_flux*0.9+a3+C1)
                self.time_constant=a1*(math.log((a2-self.equil_flux*0.9)*a3/a2/(a3-self.equil_flux*0.9)))
            elif self.species_count==(2,2) and abs(self.kb/self.kf-1)>1e-3:
                a1=math.sqrt(1/(4*self.min_reactant_conc**2*self.kf*self.kb + 8*self.min_reactant_conc*self.min_product_conc*self.kf*self.kb + 4*self.min_reactant_conc*self.delta_product_conc*self.kf*self.kb + 4*self.min_reactant_conc*self.delta_reactant_conc*self.kf*self.kb + 4*self.min_product_conc**2*self.kf*self.kb + 4*self.min_product_conc*self.delta_product_conc*self.kf*self.kb + 4*self.min_product_conc*self.delta_reactant_conc*self.kf*self.kb + self.delta_product_conc**2*self.kb**2 + 2*self.delta_product_conc*self.delta_reactant_conc*self.kf*self.kb + self.delta_reactant_conc**2*self.kf**2))
                a2=(-4*self.min_reactant_conc**2*self.kf*self.kb*a1 - 8*self.min_reactant_conc*self.min_product_conc*self.kf*self.kb*a1 - 4*self.min_reactant_conc*self.delta_product_conc*self.kf*self.kb*a1 - 4*self.min_reactant_conc*self.delta_reactant_conc*self.kf*self.kb*a1 - 2*self.min_reactant_conc*self.kf - 4*self.min_product_conc**2*self.kf*self.kb*a1 - 4*self.min_product_conc*self.delta_product_conc*self.kf*self.kb*a1 - 4*self.min_product_conc*self.delta_reactant_conc*self.kf*self.kb*a1 - 2*self.min_product_conc*self.kb - self.delta_product_conc**2*self.kb**2*a1 - 2*self.delta_product_conc*self.delta_reactant_conc*self.kf*self.kb*a1 - self.delta_product_conc*self.kb - self.delta_reactant_conc**2*self.kf**2*a1 - self.delta_reactant_conc*self.kf)/(2*(self.kf - self.kb))
                a3=(4*self.min_reactant_conc**2*self.kf*self.kb*a1 + 8*self.min_reactant_conc*self.min_product_conc*self.kf*self.kb*a1 + 4*self.min_reactant_conc*self.delta_product_conc*self.kf*self.kb*a1 + 4*self.min_reactant_conc*self.delta_reactant_conc*self.kf*self.kb*a1 - 2*self.min_reactant_conc*self.kf + 4*self.min_product_conc**2*self.kf*self.kb*a1 + 4*self.min_product_conc*self.delta_product_conc*self.kf*self.kb*a1 + 4*self.min_product_conc*self.delta_reactant_conc*self.kf*self.kb*a1 - 2*self.min_product_conc*self.kb + self.delta_product_conc**2*self.kb**2*a1 + 2*self.delta_product_conc*self.delta_reactant_conc*self.kf*self.kb*a1 - self.delta_product_conc*self.kb + self.delta_reactant_conc**2*self.kf**2*a1 - self.delta_reactant_conc*self.kf)/(2*(self.kf - self.kb))
                # print(a1,a2,a3)
                self.time_constant=a1*(math.log((a2 + self.equil_flux*0.9)*a3/(a3 + self.equil_flux*0.9)/a2))
            elif self.species_count==(2,2) and abs(self.kb/self.kf-1)<=1e-3:
                self.time_constant=(math.log(-self.min_reactant_conc**2 - self.min_reactant_conc*self.delta_reactant_conc + self.min_product_conc**2 + self.min_product_conc*self.delta_product_conc) - math.log(-self.min_reactant_conc**2 - self.min_reactant_conc*self.delta_reactant_conc + 2*self.min_reactant_conc*self.equil_flux*0.9 + self.min_product_conc**2 + self.min_product_conc*self.delta_product_conc + 2*self.min_product_conc*self.equil_flux*0.9 + self.delta_product_conc*self.equil_flux*0.9 + self.delta_reactant_conc*self.equil_flux*0.9))/(self.kf*(2*self.min_reactant_conc + 2*self.min_product_conc + self.delta_product_conc + self.delta_reactant_conc))
            else: raise Exception('Unsupported reaction type')
            # print(self.time_constant)
            # if self.time_constant.imag!=0:
            #     print(self.time_constant)
            #     print(a1,a2,a3)
            #     print(self.kf,self.kb)
            #     print(self.reactant_conc,self.coreactant_conc,self.product_conc,self.coproduct_conc)
            #     print(self.equil_flux,self.species_count)
            #     raise Exception('complex number problem')
            # print(self.time_constant)
        else:
            self.pfo_k=None
            self.time_constant=None

    def pos_flux_check(self,flux):
        if self.species_count==(1,1):
            pos_check=all(i>0 for i in [self.reactant_conc-flux,self.product_conc+flux])
        if self.species_count==(1,2):
            pos_check=all(i>0 for i in [self.reactant_conc-flux,self.product_conc+flux,self.coproduct_conc+flux])
        if self.species_count==(2,1):
            pos_check=all(i>0 for i in [self.reactant_conc-flux,self.product_conc+flux,self.coreactant_conc-flux])
        if self.species_count==(2,2):
            pos_check=all(i>0 for i in [self.reactant_conc-flux,self.product_conc+flux,self.coreactant_conc-flux,self.coproduct_conc+flux])
        # if not pos_check:
        #     print("pos check fail")
        return pos_check

    def flux_K_check(self):
        if self.species_count==(1,1):
            flux_K=(self.product_conc+self.equil_flux)/(self.reactant_conc-self.equil_flux)
        if self.species_count==(1,2):
            flux_K=(self.product_conc+self.equil_flux)*(self.coproduct_conc+self.equil_flux)/(self.reactant_conc-self.equil_flux)
        if self.species_count==(2,1):
            flux_K=(self.product_conc+self.equil_flux)/(self.reactant_conc-self.equil_flux)/(self.coreactant_conc-self.equil_flux)
        if self.species_count==(2,2):
            flux_K=(self.product_conc+self.equil_flux)*(self.coproduct_conc+self.equil_flux)/(self.reactant_conc-self.equil_flux)/(self.coreactant_conc-self.equil_flux)
        res_K=abs(math.log10(self.Keq/flux_K))
        if res_K>1e-5: #revisit numerical precision of solution
            if abs((self.min_product_conc+self.equil_flux)/self.min_product_conc)<1e-5 or abs((self.min_reactant_conc-self.equil_flux)/self.min_reactant_conc)<1e-5:
                return True #depleted species
            print('res K',self.Keq,flux_K)
            print(self.flux_guess,self.equil_flux,self.species_count)
            print(self.min_reactant_conc,self.delta_reactant_conc,self.min_product_conc,self.delta_product_conc)
            return False
        else:
            return True
