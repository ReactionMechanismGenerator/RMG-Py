#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

#global imports
import copy
import os.path
import numpy as np
import re

#local imports
from rmgpy.chemkin import getSpeciesIdentifier
from rmgpy.scoop_framework.util import broadcast, get, WorkerWrapper, map_
from rmgpy.scoop_framework.util import logger as logging
from rmgpy.rmg.main import RMG

from model import ReductionReaction
from rates import isImportant


#global variables
reactions = None


def simulate_one(reactionModel, atol, rtol, reactionSystem):
    """

    Simulates one reaction system, listener registers results, 
    which are returned at the end.


    The returned data consists of a array of the species names, 
    and the concentration data.

    The concentration data consists of a number of elements for each timestep 
    the solver took to reach the end time of the batch reactor simulation.

    Each element consists of the time and the concentration data of the species at that 
    particular timestep in the order of the species names.

    """

    #register as a listener
    listener = ConcentrationListener()

    coreSpecies = reactionModel.core.species
    regex = r'\([0-9]+\)'#cut of '(one or more digits)'
    species_names = []
    for spc in coreSpecies:
        name = getSpeciesIdentifier(spc)
        name_cutoff = re.split(regex, name)[0]
        species_names.append(name_cutoff)

    listener.species_names = species_names

    reactionSystem.attach(listener)

    pdepNetworks = []
    for source, networks in reactionModel.networkDict.items():
        pdepNetworks.extend(networks)
    
    terminated, obj = reactionSystem.simulate(
        coreSpecies = reactionModel.core.species,
        coreReactions = reactionModel.core.reactions,
        edgeSpecies = reactionModel.edge.species,
        edgeReactions = reactionModel.edge.reactions,
        toleranceKeepInEdge = 0,
        toleranceMoveToCore = 1,
        toleranceInterruptSimulation = 1,
        pdepNetworks = pdepNetworks,
        absoluteTolerance = atol,
        relativeTolerance = rtol,
    ) 

    assert terminated

    #unregister as a listener
    reactionSystem.detach(listener) 

    return listener.species_names, listener.data

def simulate_all(rmg):
    """
    Simulate the RMG job, 
    for each of the simulated reaction systems.

    Each element i of the data corresponds to a reaction system.
    """
    reactionModel = rmg.reactionModel

    data = []

    atol, rtol = rmg.absoluteTolerance, rmg.relativeTolerance
    for reactionSystem in rmg.reactionSystems:
        data.append(simulate_one(reactionModel, atol, rtol, reactionSystem))

    return data
        

def initialize(wd, rxns):
    global working_dir, reactions
    working_dir = wd
    assert os.path.isdir(working_dir)
    
    #set global variable here such that functions executed in the root worker have access to it.
    
    reactions = [ReductionReaction(rxn) for rxn in rxns]
    broadcast(reactions, 'reactions')
    

def retrieve_reactions():
    """
    Reactions can be retrieved either through the global variable 'reactions' if parallel computing
    is not used.

    With the use of multiple workers, the reactions are retrieved from the previously broadcasted 
    constant.

    In any case, the references to the original reactions of the reaction model are assumed to be 
    broken.

    """
    global reactions    

    broadcasted_reactions = get('reactions')
    if broadcasted_reactions:
        reactions = broadcasted_reactions
    return reactions

def find_important_reactions(rmg, tolerance):
    """
    This function:

    - loops over all the species involved in a specific reaction
    - decides whether the specific reaction is important for the species.

    Whenever it is found that a reaction is important for a species, we break
    the species loop, and keep the reaction in the model.


    Returns:
        a list of rxns that can be removed.
    """
    
    # run the simulation, creating concentration profiles for each reaction system defined in input.
    simdata = simulate_all(rmg)


    reduce_reactions = retrieve_reactions()

    def chunks(l, n):
        """Yield successive n-sized chunks from l."""
        for i in xrange(0, len(l), n):
            yield l[i:i+n]

    CHUNKSIZE = 40
    boolean_array = []
    for chunk in chunks(reduce_reactions,CHUNKSIZE):
        N = len(chunk)
        partial_results = list(
            map_(
                WorkerWrapper(assess_reaction), chunk, [rmg.reactionSystems] * N, [tolerance] * N, [simdata] * N
                )
            )
        boolean_array.extend(partial_results)

    """
    Assuming that the order of the reduced reactions array and the core reactions of the reaction model
    are identical, iterate over the boolean array and retain those reactions of the reaction model
    that are deemed 'important'.
    """
    important_rxns = []
    for isImport, rxn in zip(boolean_array, rmg.reactionModel.core.reactions):
        logging.debug('Is rxn {rxn} important? {isImport}'.format(**locals()))
        if isImport:
            important_rxns.append(rxn)


    return important_rxns

def assess_reaction(rxn, reactionSystems, tolerance, data):
    """
    Returns whether the reaction is important or not in the reactions.

    It iterates over the reaction systems, and loads the concentration profile 
    of each reaction system.

    It iterates over a number of samples in profile and 
    evaluates the importance of the reaction at every sample.

    """


    logging.debug('Assessing reaction {}'.format(rxn))

    reactions = retrieve_reactions()    

    # read in the intermediate state variables

    for datum, reactionSystem in zip(data, reactionSystems):    
        T, P = reactionSystem.T.value_si, reactionSystem.P.value_si
        
        species_names, profile = datum

        # take N evenly spaced indices from the table with simulation results:

        """

        Number of time steps between start and end time of the batch reactor simulation at which the importance of 
        reactions should be evaluated.



        The more timesteps, the less chance we have to remove an important reactions, but the more simulations
        need to be carried out.
        """
        
        timesteps = len(profile) / 2
        logging.debug('Evaluating the importance of a reaction at {} time samples.'.format(timesteps))

        assert timesteps <= len(profile)
        indices = map(int, np.linspace(0, len(profile)-1, num = timesteps))
        for index in indices:
            assert profile[index] is not None
            timepoint, coreSpeciesConcentrations = profile[index]

            coreSpeciesConcentrations = {key: float(value) for (key, value) in zip(species_names, coreSpeciesConcentrations)}
            
            for species_i in rxn.reactants:
                if isImportant(rxn, species_i, reactions, 'reactant', tolerance, T, P, coreSpeciesConcentrations):
                    return True

            #only continue if the reaction is not important yet.
            for species_i in rxn.products:
                if isImportant(rxn, species_i, reactions, 'product', tolerance, T, P, coreSpeciesConcentrations):
                    return True

    return False

    
def search_target_index(target_label, reactionModel):
    """
    Searches for the Species object in the core species
    of the reaction that has the same label as the parameter string.
    """

    for i, spc in enumerate(reactionModel.core.species):
        if spc.label == target_label:
            return i

    raise Exception('{} could not be found...'.format(target_label))


def compute_observables(targets, reactionModel, reactionSystem, atol, rtol):
    """
    Computes the observables of the targets, provided in the function signature.

    Currently, the species mole fractions at the end time of the
    batch reactor simulation are the only observables that can be computed.

    - resetting the reaction system, initialing with empty variables
    - running the simulation at the conditions stored in the reaction system
    """
    reactionSystem.initializeModel(\
        reactionModel.core.species, reactionModel.core.reactions,\
        reactionModel.edge.species, reactionModel.edge.reactions, \
        [], atol, rtol)
        
    #reset reaction system variables:
    logging.info('No. of rxns in core reactions: {}'.format(len(reactionModel.core.reactions)))

    #run the simulation:
    simulate_one(reactionModel, atol, rtol, reactionSystem)

    observables = compute_mole_fractions(targets, reactionModel, reactionSystem)

    return observables

def compute_mole_fractions(targets, reactionModel, reactionSystem):
    """
    Computes the mole fractions of the targets, identified by the list 
    of species names in the function signature.

    Returns a numpy array with the mole fractions at the end time of the reactor
    simulation.

    - searching the index of the target species in the core species
    of the global reduction variable
    - fetching the computed moles variable y

    """
    mole_fractions = np.zeros(len(targets), np.float64)

    for i, label in enumerate(targets):
        target_index = search_target_index(label, reactionModel)

        mole_fractions[i] = reactionSystem.y[target_index]

    return mole_fractions

def compute_conversion(target_label, reactionModel, reactionSystem, atol, rtol):
    """
    Computes the conversion of a target molecule by

    - searching the index of the target species in the core species
    of the global reduction variable
    - resetting the reaction system, initialing with empty variables
    - fetching the initial moles variable y0
    - running the simulation at the conditions stored in the reaction system
    - fetching the computed moles variable y
    - computing conversion
    """

    target_index = search_target_index(target_label, reactionModel)

    #reset reaction system variables:
    logging.info('No. of rxns in core reactions: {}'.format(len(reactionModel.core.reactions)))

    reactionSystem.initializeModel(\
        reactionModel.core.species, reactionModel.core.reactions,\
        reactionModel.edge.species, reactionModel.edge.reactions, \
        [], atol, rtol)

    #get the initial moles:
    y0 = reactionSystem.y.copy()

    #run the simulation:
    simulate_one(reactionModel, atol, rtol, reactionSystem)

    #compute conversion:
    conv = 1 - (reactionSystem.y[target_index] / y0[target_index])
    return conv

def reduce_model(tolerance, targets, reactionModel, rmg, reaction_system_index):
    """
    Reduces the model for the given tolerance and evaluates the 
    target observables.
    """

    # reduce model with the tolerance specified earlier:
    important_reactions = find_important_reactions(rmg, tolerance)

    original_size = len(reactionModel.core.reactions)

    no_important_reactions = len(important_reactions)
    logging.info('Number of important reactions: {}'.format(no_important_reactions))

    #set the core reactions to the reduced reaction set:
    original_reactions = reactionModel.core.reactions
    rmg.reactionModel.core.reactions = important_reactions

    #re-compute observables: 
    observables = compute_observables(targets, rmg.reactionModel,\
     rmg.reactionSystems[reaction_system_index],\
     rmg.absoluteTolerance, rmg.relativeTolerance)

    #reset the reaction model to its original state:
    rmg.reactionModel.core.reactions = original_reactions

    logging.info('Observables of reduced model ({} rxns):'.format(no_important_reactions))
    for target, observable in zip(targets, observables):
        logging.info('{}: {:.2f}%'.format(target, observable * 100))

    return observables, important_reactions

class ConcentrationListener(object):
    """Returns the species concentration profiles at each time step."""

    def __init__(self):
        self.species_names = []
        self.data = []

    def update(self, subject):
        """
        Register the time (t) and the species mole fractions at the
        given time.

        The snapshots variable stores time and Volume as the first two 
        elements in the array.
        """
        data = subject.snapshots
        self.data = process(data)

def process(data):
    """
    The data is structured as a list of lists.

    Each list contains [time, Volume, [species mole fractions]]

    The volume is cut out of each list, the remaining part is stored as a tuple.
    """
    processed = []

    for d in data:
        processed.append((d[0], d[2:]))

    return processed
