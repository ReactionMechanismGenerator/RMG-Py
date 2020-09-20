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
from collections import defaultdict
from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import rdFMCS
import numpy as np
import math

def search_priority(rmg,reaction_system):
    # print(dir(reaction_system))
    sim_T=reaction_system.T.value
    core_rxn_kf=reaction_system.kf
    core_rxn_kb=reaction_system.kb
    core_concentrations=reaction_system.core_species_concentrations
    core_reactants_idx=reaction_system.reactant_indices.tolist()
    core_products_idx=reaction_system.product_indices.tolist()
    #get species list
    core_species=rmg.reaction_model.core.species
    # print(dir(core_species[5]))
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
    allowed_multiplicity=dict([(0,[0,2,3,6]),(1,[1,4,7]),(2,[2,0,6,8]),(3,[0,3,6]),(4,[1,4,7]),(5,[2,5,6]),(6,[6,0,2,3,5,8]),(7,[7,1,4]),(8,[8,2,6])])
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
                possible_connections.append(Connection(core_species,core_concentrations,reactant_index=pr0['reactants'][0],product_index=pr1['reactants'][0],temperature=sim_T))
            elif len(pr0['reactants'])==2 and len(pr1['reactants'])==1:
                possible_connections.append(Connection(core_species,core_concentrations,reactant_index=pr0['reactants'][0],coreactant_index=pr0['reactants'][1],product_index=pr1['reactants'][0],temperature=sim_T))
                possible_connections.append(Connection(core_species,core_concentrations,reactant_index=pr0['reactants'][1],coreactant_index=pr0['reactants'][0],product_index=pr1['reactants'][0],temperature=sim_T))
            elif len(pr0['reactants'])==1 and len(pr1['reactants'])==2:
                possible_connections.append(Connection(core_species,core_concentrations,reactant_index=pr0['reactants'][0],coproduct_index=pr1['reactants'][1],product_index=pr1['reactants'][0],temperature=sim_T))
                possible_connections.append(Connection(core_species,core_concentrations,reactant_index=pr0['reactants'][0],coproduct_index=pr1['reactants'][0],product_index=pr1['reactants'][1],temperature=sim_T))
            elif len(pr0['reactants'])==2 and len(pr1['reactants'])==2:
                possible_connections.append(Connection(core_species,core_concentrations,reactant_index=pr0['reactants'][0],coreactant_index=pr0['reactants'][1],coproduct_index=pr1['reactants'][1],product_index=pr1['reactants'][0],temperature=sim_T))
                possible_connections.append(Connection(core_species,core_concentrations,reactant_index=pr0['reactants'][1],coreactant_index=pr0['reactants'][0],coproduct_index=pr1['reactants'][0],product_index=pr1['reactants'][1],temperature=sim_T))
                possible_connections.append(Connection(core_species,core_concentrations,reactant_index=pr0['reactants'][0],coreactant_index=pr0['reactants'][1],coproduct_index=pr1['reactants'][0],product_index=pr1['reactants'][1],temperature=sim_T))
                possible_connections.append(Connection(core_species,core_concentrations,reactant_index=pr0['reactants'][1],coreactant_index=pr0['reactants'][0],coproduct_index=pr1['reactants'][1],product_index=pr1['reactants'][0],temperature=sim_T)) 
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
                if len(r_ids)==1 and len(p_ids)==1:
                    direct_connections_list.append(Connection(core_species,core_concentrations,
                                                            reactant_index=r_ids[j],
                                                            product_index=p_ids[k],
                                                            kf=core_rxn_kf[i],
                                                            kb=core_rxn_kb[i],
                                                            temperature=sim_T))
                    direct_connections_list.append(Connection(core_species,core_concentrations,
                                                            reactant_index=p_ids[k],
                                                            product_index=r_ids[j],
                                                            kf=core_rxn_kb[i],
                                                            kb=core_rxn_kf[i],
                                                            temperature=sim_T))
                if len(r_ids)==2 and len(p_ids)==1:
                    direct_connections_list.append(Connection(core_species,core_concentrations,
                                                            reactant_index=r_ids[j],
                                                            coreactant_index=r_ids[1-j],
                                                            product_index=p_ids[k],
                                                            kf=core_rxn_kf[i],
                                                            kb=core_rxn_kb[i],
                                                            temperature=sim_T))
                    direct_connections_list.append(Connection(core_species,core_concentrations,
                                                            reactant_index=p_ids[k],
                                                            product_index=r_ids[j],
                                                            coproduct_index=r_ids[1-j],
                                                            kf=core_rxn_kb[i],
                                                            kb=core_rxn_kf[i],
                                                            temperature=sim_T))
                if len(r_ids)==1 and len(p_ids)==2:
                    direct_connections_list.append(Connection(core_species,core_concentrations,
                                                            reactant_index=r_ids[j],
                                                            product_index=p_ids[k],
                                                            coproduct_index=p_ids[1-k],
                                                            kf=core_rxn_kf[i],
                                                            kb=core_rxn_kb[i],
                                                            temperature=sim_T))
                    direct_connections_list.append(Connection(core_species,core_concentrations,
                                                            reactant_index=p_ids[k],
                                                            coreactant_index=p_ids[1-k],
                                                            product_index=r_ids[j],
                                                            kf=core_rxn_kb[i],
                                                            kb=core_rxn_kf[i],
                                                            temperature=sim_T))
                if len(r_ids)==2 and len(p_ids)==2:
                    direct_connections_list.append(Connection(core_species,core_concentrations,
                                                            reactant_index=r_ids[j],
                                                            coreactant_index=r_ids[1-j],
                                                            product_index=p_ids[k],
                                                            coproduct_index=p_ids[1-k],
                                                            kf=core_rxn_kf[i],
                                                            kb=core_rxn_kb[i],
                                                            temperature=sim_T))
                    direct_connections_list.append(Connection(core_species,core_concentrations,
                                                            reactant_index=p_ids[k],
                                                            coreactant_index=p_ids[1-k],
                                                            product_index=r_ids[j],
                                                            coproduct_index=r_ids[1-j],
                                                            kf=core_rxn_kb[i],
                                                            kb=core_rxn_kf[i],
                                                            temperature=sim_T))
#remove already covered connections
    removed_connections=set()
    for conn0 in direct_connections_list:
        conn0_hash=conn0.hash
        for conn1 in possible_connections:
            conn1_hash=conn1.hash
            if conn0_hash==conn1_hash:
                possible_connections.remove(conn1)
                removed_connections.add(conn1)
    # print(len(possible_connections))
    # print(len(removed_connections))
    # for i in direct_connections_list: print(i.time_constant)
    pooling_time_cycle(direct_connections_list,core_species)

#forward pooling list at a given time
#bidirectional check on pooling at a given time
#cycle through different times
#classify connections by pooling time, delta H, subgroup similarity

def equilibrium_pooling(direct_connections_list,core_species,t_characteristic): #,previous_forward_pass=None):
    fast_connections=[]
    strong_connections=[]
    for rct_idx in range(len(core_species)):
        fast_connections.append(set())
        strong_connections.append(set())
    for conn in direct_connections_list:
        if conn.time_constant is None:
            continue
        elif conn.coproduct is not None and conn.coreactant is not None:
            if conn.coproduct_conc<conn.product_conc or conn.coreactant_conc<conn.reactant_conc:
                continue
            elif conn.time_constant<t_characteristic:
                fast_connections[conn.reactant].add(conn.product)
                if conn.min_reactant_conc>0.1*conn.min_product_conc:
                    strong_connections[conn.reactant].add(conn.product)
        elif conn.time_constant<t_characteristic:
            fast_connections[conn.reactant].add(conn.product)
            if conn.min_reactant_conc>0.1*conn.min_product_conc:
                strong_connections[conn.reactant].add(conn.product)
        # print("connections",conn.time_constant,conn.reactant,conn.product)
    # if previous_forward_pass is not None:
        # forward_pass=previous_forward_pass
    # else:
    # print('fast_connections')
    # for i in fast_connections: print(i)
    forward_pass=[]
    strong_pass=[]
    for rct_idx in range(len(core_species)):
        forward_pass.append(set([rct_idx]))
        strong_pass.append(set([rct_idx]))
    for iter in range(len(core_species)):
        for rct_idx in range(len(core_species)):
            for prd_idx in fast_connections[rct_idx]:
                forward_pass[rct_idx].update(forward_pass[prd_idx])
            for prd_idx in strong_connections[rct_idx]:
                strong_pass[rct_idx].update(strong_pass[prd_idx])
    # print('forward_pass',forward_pass)
    # print('strong_pass',strong_pass)
    # for i in forward_pass: print(i)
    pool_groups=[]
    for rct_idx in range(len(core_species)):
        pool_groups.append(set())
    broken_pools=set()
    for rct_idx in range(len(core_species)):
        for prd_idx in forward_pass[rct_idx]:
            if rct_idx in forward_pass[prd_idx]:
                pool_groups[rct_idx].add(prd_idx)
            else: 
                broken_pools.add(rct_idx)
                if prd_idx in strong_pass[rct_idx]:
                    broken_pools.add(prd_idx)
    # print("broken_pools",broken_pools)
    # print('pools before')
    # for i,group in enumerate(pool_groups):print(i,group)
    for rct_idx in broken_pools:
        for prd_idx in pool_groups[rct_idx]:
            pool_groups[prd_idx]=set()
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
    pooling_times[:]=np.nan
    log10_minimum=math.floor(math.log10(time_constants[0]))
    log10_maximum=math.ceil(math.log10(time_constants[-1]))+1
    # log10_minimum=math.floor(math.log2(time_constants[0]))
    # log10_maximum=math.ceil(math.log2(time_constants[-1]))+1
    # pool_groups=equilibrium_pooling(direct_connections_list,core_species,10**log10_maximum)
    for log10_time in list(range(log10_minimum,log10_maximum))[::-1]:
        pool_groups=equilibrium_pooling(direct_connections_list,core_species,10**log10_time)
        # pool_groups=equilibrium_pooling(direct_connections_list,core_species,2**log10_time)
        for i,grp in enumerate(pool_groups):
            for j in grp:
                pooling_times[i,j]=log10_time
    for i in pooling_times: print(i) 
    # print(pool_groups)
    for i,spc in enumerate(core_species): print(i,spc.label)
    # for i in pool_groups[-1]:
        # print(core_species[i].label)
    # equilibrium_pooling(direct_connections_list,core_species,10**-15)
    return pooling_times



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
    `rate`                  ``float``
    `temperature`           ``float``
    `pfo_rate`              ``float``
    `time_constant`         ``float``
    `H_rxn`                 ``float``
    `largest_substructure`  ``int``
    `pooling_time`          ``float``
    `hash`                  ``tuple``
    `min_reactant_conc`
    `min_product_conc`
   ======================== ================= =========================================

    """

    def __init__(self,core_species,core_concentrations,reactant_index,product_index,coreactant_index=None,coproduct_index=None,temperature=298,kf=None,kb=None):
        self.reactant=reactant_index
        self.reactant_conc=core_concentrations[reactant_index]
        self.product=product_index
        self.product_conc=core_concentrations[product_index]
        if coreactant_index is not None:
            self.coreactant=coreactant_index
            self.coreactant_conc=core_concentrations[coreactant_index]
            self.min_reactant_conc=min([self.coreactant_conc,self.reactant_conc])
        else:
            self.coreactant=None
            self.coreactant_conc=None
            self.min_reactant_conc=self.reactant_conc
        if coproduct_index is not None:
            self.coproduct=coproduct_index
            self.coproduct_conc=core_concentrations[coproduct_index]
            self.min_product_conc=min([self.coproduct_conc,self.product_conc])
        else:
            self.coproduct=None
            self.coproduct_conc=None
            self.min_product_conc=self.product_conc
        self.kf=kf
        self.kb=kb
        self.T=temperature
        self.pooling_time=None
        self.hash=set(((reactant_index,coreactant_index if coreactant_index is not None else None),
            (product_index,coproduct_index if coproduct_index is not None else None)))

        if kf is not None:
            self.pfo_k=kf*self.coreactant_conc if coreactant_index is not None else kf
            if self.pfo_k!=0:
                self.time_constant=1/self.pfo_k
            else:
                self.time_constant = None
        else:
            self.pfo_k=None
            self.time_constant=None
        self.H_rxn=core_species[self.product].get_enthalpy(temperature)-core_species[self.reactant].get_enthalpy(temperature)
        if coproduct_index is not None:
            self.H_rxn+=core_species[self.coproduct].get_enthalpy(temperature)
        if coreactant_index is not None:
            self.H_rxn-=core_species[self.coreactant].get_enthalpy(temperature)
        
        lg=RDLogger.logger()
        lg.setLevel(RDLogger.ERROR)
        rdkit_rct=Chem.MolFromSmiles(core_species[reactant_index].molecule[0].smiles)
        rdkit_prd=Chem.MolFromSmiles(core_species[product_index].molecule[0].smiles)
        substruct=rdFMCS.FindMCS([rdkit_rct,rdkit_prd])
        self.largest_substructure=substruct.numBonds
    
    def ratchet_pooling_time(self,new_time):
        if self.pooling_time is None:
            self.pooling_time=new_time
        elif self.pooling_time>new_time:
            self.pooling_time=new_time
    
    def set_pooling_time(self,new_time):
        self.pooling_time=new_time

    def is_same_connection(self,other_connection): #there's a bug here somehow with how it's filtering. self.reactant is a looser filter than self.reactant.label
        # if self.reactant == other_connection.reactant:
        #     if (self.product,self.coreactant,self.coproduct) == (other_connection.product,other_connection.coreactant,other_connection.coproduct):return True
        # if self.reactant == other_connection.product:
        #     if (self.product,self.coreactant,self.coproduct) == (other_connection.reactant,other_connection.coproduct,other_connection.coreactant): return True
        # return False
        return self.hash==other_connection.hash