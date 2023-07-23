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

import os


import numpy as np

import rmgpy
from rmgpy.quantity import Quantity
from rmgpy.tools.canteramodel import find_ignition_delay, CanteraCondition, Cantera


class CanteraTest:
    def test_ignition_delay(self):
        """
        Test that find_ignition_delay() works.
        """

        t = np.arange(0, 5, 0.5)
        P = np.array([0, 0.33, 0.5, 0.9, 2, 4, 15, 16, 16.1, 16.2])
        OH = np.array([0, 0.33, 0.5, 0.9, 2, 4, 15, 16, 7, 2])
        CO = OH * 0.9

        t_ign = find_ignition_delay(t, P)
        assert t_ign == 2.75

        t_ign = find_ignition_delay(t, OH, "maxHalfConcentration")
        assert t_ign == 3

        t_ign = find_ignition_delay(t, [OH, CO], "maxSpeciesConcentrations")
        assert t_ign == 3.5

    def test_repr(self):
        """
        Test that the repr function for a CanteraCondition object can reconstitute
        the same object
        """
        reactor_type = "IdealGasReactor"
        mol_frac = {"CC": 0.05, "[Ar]": 0.95}
        P = (3, "atm")
        T = (1500, "K")
        termination_time = (5e-5, "s")
        condition = CanteraCondition(reactor_type, termination_time, mol_frac, T0=T, P0=P)
        repr_condition = eval(condition.__repr__())
        assert repr_condition.T0.value_si == Quantity(T).value_si
        assert repr_condition.P0.value_si == Quantity(P).value_si
        assert repr_condition.V0 == None
        assert repr_condition.mol_frac == mol_frac


class RMGToCanteraTest:
    """
    Contains unit tests for the conversion of RMG species and reaction objects to Cantera objects.
    """

    def setup_class(self):
        from rmgpy.chemkin import load_chemkin_file

        folder = os.path.join(os.path.dirname(rmgpy.__file__), "tools/data/various_kinetics")

        chemkin_path = os.path.join(folder, "chem_annotated.inp")
        dictionary_path = os.path.join(folder, "species_dictionary.txt")
        transport_path = os.path.join(folder, "tran.dat")

        species, reactions = load_chemkin_file(chemkin_path, dictionary_path, transport_path)

        self.rmg_ctSpecies = [spec.to_cantera(use_chemkin_identifier=True) for spec in species]
        self.rmg_ctReactions = []
        for rxn in reactions:
            converted_reactions = rxn.to_cantera(species, use_chemkin_identifier=True)
            if isinstance(converted_reactions, list):
                self.rmg_ctReactions.extend(converted_reactions)
            else:
                self.rmg_ctReactions.append(converted_reactions)
        job = Cantera()
        job.load_chemkin_model(chemkin_path, transport_file=transport_path, quiet=True)
        self.ctSpecies = job.model.species()
        self.ctReactions = job.model.reactions()

    def test_species_conversion(self):
        """
        Test that species objects convert properly
        """
        from rmgpy.tools.canteramodel import check_equivalent_cantera_species

        for i in range(len(self.ctSpecies)):
            assert check_equivalent_cantera_species(self.ctSpecies[i], self.rmg_ctSpecies[i])

    def test_reaction_conversion(self):
        """
        Test that species objects convert properly
        """
        from rmgpy.tools.canteramodel import check_equivalent_cantera_reaction

        for i in range(len(self.ctReactions)):
            assert check_equivalent_cantera_reaction(self.ctReactions[i], self.rmg_ctReactions[i])
