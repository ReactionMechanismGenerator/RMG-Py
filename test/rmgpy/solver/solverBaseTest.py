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

import os.path
import pickle


import rmgpy
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.tools.loader import load_rmg_py_job


class ConcentrationPrinter(object):
    def __init__(self):
        self.species_names = []
        self.data = []

    def update(self, subject):
        self.data.append((subject.t, subject.core_species_concentrations))


class ReactionSystemTest:
    def setup_class(self):
        self.listener = ConcentrationPrinter()

        folder = os.path.join(os.path.dirname(rmgpy.__file__), "solver/files/listener/")
        input_file = os.path.join(folder, "input.py")
        chemkin_file = os.path.join(folder, "chemkin/chem.inp")
        spc_dict = os.path.join(folder, "chemkin/species_dictionary.txt")

        self.rmg = load_rmg_py_job(
            input_file,
            chemkin_file,
            spc_dict,
            generate_images=False,
            check_duplicates=False,
        )

    def test_surface_initialization(self):
        """
        test that initialize_surface is correctly removing species and reactions when
        they are no longer consistent with the surface (due to other species/reactions moving to the
        bulk core)
        """
        reaction_system = self.rmg.reaction_systems[0]
        reaction_system.attach(self.listener)
        reaction_model = self.rmg.reaction_model

        core_species = reaction_model.core.species
        core_reactions = reaction_model.core.reactions
        surface_species = [core_species[7], core_species[6]]
        surface_reactions = [core_reactions[0], core_reactions[2], core_reactions[3]]

        reaction_system.initialize_model(
            core_species,
            core_reactions,
            reaction_model.edge.species,
            reaction_model.edge.reactions,
            surface_species,
            surface_reactions,
        )

        assert len(surface_species) == 1  # only H should be left
        assert len(surface_reactions) == 2  # all the reactions with H should stay

    def test_surface_layering_constraint(self):
        """
        test that the correct maximum under the surface layering constraint is being
        found
        """
        reaction_system = self.rmg.reaction_systems[0]
        reaction_system.attach(self.listener)
        reaction_model = self.rmg.reaction_model
        core_species = reaction_model.core.species
        core_reactions = reaction_model.core.reactions

        edge_species = [core_species[6], core_species[7]]
        edge_reactions = core_reactions[1:]
        surface_species = [core_species[5]]
        surface_reactions = [core_reactions[0]]
        core_species = core_species[0:6] + [core_species[8]]
        core_reactions = surface_reactions[:]
        reaction_system.num_core_reactions = 1
        reaction_system.num_core_species = 7

        reaction_system.initialize_model(
            core_species,
            core_reactions,
            edge_species,
            edge_reactions,
            surface_species,
            surface_reactions,
        )

        assert len(reaction_system.surface_species_indices) == 1  # surface_species_indices calculated correctly
        assert reaction_system.surface_species_indices[0] == 5  # surface_species_indices calculated correctly

        inds = reaction_system.get_layering_indices()

        assert inds[0] == 1  # worked correctly
        assert inds[1] == 2

    def test_add_reactions_to_surface(self):
        """
        Test that add_reactions_to_surface gives the correct surface_species and surface_reactions lists after being called
        """
        reaction_system = self.rmg.reaction_systems[0]
        reaction_system.attach(self.listener)
        reaction_model = self.rmg.reaction_model
        species = reaction_model.core.species
        reactions = reaction_model.core.reactions

        core_species = species[0:6]
        core_reactions = [reactions[0]]
        surface_species = []
        surface_reactions = []
        edge_species = species[6:]
        edge_reactions = reactions[1:]

        reaction_system.initialize_model(
            core_species,
            core_reactions,
            edge_species,
            edge_reactions,
            surface_species,
            surface_reactions,
        )

        new_surface_reactions = edge_reactions
        new_surface_reaction_inds = [edge_reactions.index(i) for i in new_surface_reactions]

        surface_species, surface_reactions = reaction_system.add_reactions_to_surface(
            new_surface_reactions,
            new_surface_reaction_inds,
            surface_species,
            surface_reactions,
            edge_species,
        )

        assert set(surface_species) == set(edge_species)  # all edge species should now be in the surface
        assert set(surface_reactions) == set(edge_reactions)  # all edge reactions should now be in the surface

    def test_attach_detach(self):
        """
        Test that a ReactionSystem listener can be attached/detached.
        """
        # create observable

        reaction_system = self.rmg.reaction_systems[0]
        reaction_system.attach(self.listener)
        assert reaction_system._observers != []

        reaction_system.detach(self.listener)
        assert reaction_system._observers == []

    def test_listen(self):
        """
        Test that data can be retrieved from an attached ReactionSystem listener.
        """
        # create observable
        reaction_system = self.rmg.reaction_systems[0]
        reaction_system.attach(self.listener)

        reaction_model = self.rmg.reaction_model

        assert self.listener.data == []

        model_settings = ModelSettings(tol_move_to_core=1, tol_keep_in_edge=0, tol_interrupt_simulation=1)
        simulator_settings = SimulatorSettings()

        # run simulation:
        terminated, resurrected, obj, sspcs, srxns, t, conv = reaction_system.simulate(
            core_species=reaction_model.core.species,
            core_reactions=reaction_model.core.reactions,
            edge_species=reaction_model.edge.species,
            edge_reactions=reaction_model.edge.reactions,
            surface_species=[],
            surface_reactions=[],
            model_settings=model_settings,
            simulator_settings=simulator_settings,
        )

        assert self.listener.data != []

    def test_pickle(self):
        """
        Test that a ReactionSystem object can be un/pickled.
        """
        rxn_sys1 = self.rmg.reaction_systems[0]
        rxn_sys = pickle.loads(pickle.dumps(rxn_sys1))

        assert rxn_sys is not None
        assert isinstance(rxn_sys, rmgpy.solver.simple.SimpleReactor)
        assert rxn_sys.T.value_si == rxn_sys1.T.value_si
        assert rxn_sys.P.value_si == rxn_sys1.P.value_si
        assert rxn_sys.termination[0].conversion == rxn_sys1.termination[0].conversion
        assert rxn_sys.termination[1].time.value_si == rxn_sys1.termination[1].time.value_si
