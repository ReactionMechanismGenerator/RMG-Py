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

import itertools
import os


import numpy as np

from rmgpy import settings
from rmgpy.data.kinetics import TemplateReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.rmg.main import RMG
from rmgpy.rmg.react import react, react_all
from rmgpy.species import Species

TESTFAMILIES = [
    "H_Abstraction",
    "R_Recombination",
    "Disproportionation",
    "R_Addition_MultipleBond",
]


class TestReact:
    def setup_class(self):
        """
        A method that is run before each unit test in this class.
        """
        # set-up RMG object
        self.rmg = RMG()

        # load kinetic database and forbidden structures
        self.rmg.database = RMGDatabase()
        path = os.path.join(settings["database.directory"])

        # forbidden structure loading
        self.rmg.database.load_forbidden_structures(os.path.join(path, "forbiddenStructures.py"))
        # kinetics family loading
        self.rmg.database.load_kinetics(
            os.path.join(path, "kinetics"),
            kinetics_families=TESTFAMILIES,
            reaction_libraries=[],
        )

    def test_react(self):
        """
        Test that the ``react`` function works in serial
        """
        procnum = 1

        spc_a = Species().from_smiles("[OH]")
        spcs = [Species().from_smiles("CC"), Species().from_smiles("[CH3]")]
        spc_tuples = [((spc_a, spc), ["H_Abstraction"]) for spc in spcs]

        reaction_list = list(itertools.chain.from_iterable(react(spc_tuples, procnum)))
        assert reaction_list is not None
        assert len(reaction_list) == 3
        assert all([isinstance(rxn, TemplateReaction) for rxn in reaction_list])

    def test_react_parallel(self):
        """
        Test that the ``react`` function works in parallel using Python multiprocessing
        """
        import rmgpy.rmg.main

        rmgpy.rmg.main.maxproc = 2
        procnum = 2

        spc_a = Species().from_smiles("[OH]")
        spcs = [Species().from_smiles("CC"), Species().from_smiles("[CH3]")]
        spc_tuples = [((spc_a, spc), ["H_Abstraction"]) for spc in spcs]

        reaction_list = list(itertools.chain.from_iterable(react(spc_tuples, procnum)))
        assert reaction_list is not None
        assert len(reaction_list) == 3
        assert all([isinstance(rxn, TemplateReaction) for rxn in reaction_list])

        # Reset module level maxproc back to default
        rmgpy.rmg.main.maxproc = 1

    def test_react_all(self):
        """
        Test that the ``react_all`` function works in serial
        """
        procnum = 1

        spcs = [
            Species().from_smiles("C=C"),
            Species().from_smiles("[CH3]"),
            Species().from_smiles("[OH]"),
            Species().from_smiles("CCCCCCCCCCC"),
        ]

        n = len(spcs)
        reaction_list, spc_tuples = react_all(spcs, n, np.ones(n), np.ones([n, n]), np.ones([n, n, n]), procnum)
        assert reaction_list is not None
        assert len(reaction_list) == 34
        assert len(spc_tuples) == 34

        flat_rxn_list = list(itertools.chain.from_iterable(reaction_list))
        assert len(flat_rxn_list) == 44
        assert all([isinstance(rxn, TemplateReaction) for rxn in flat_rxn_list])

    def test_react_all_parallel(self):
        """
        Test that the ``react_all`` function works in parallel using Python multiprocessing
        """
        import rmgpy.rmg.main

        rmgpy.rmg.main.maxproc = 2
        procnum = 2

        spcs = [
            Species().from_smiles("C=C"),
            Species().from_smiles("[CH3]"),
            Species().from_smiles("[OH]"),
            Species().from_smiles("CCCCCCCCCCC"),
        ]

        n = len(spcs)
        reaction_list, spc_tuples = react_all(spcs, n, np.ones(n), np.ones([n, n]), np.ones([n, n, n]), procnum)
        assert reaction_list is not None
        assert len(reaction_list) == 94
        assert len(spc_tuples) == 94

        flat_rxn_list = list(itertools.chain.from_iterable(reaction_list))
        assert len(flat_rxn_list) == 44
        assert all([isinstance(rxn, TemplateReaction) for rxn in flat_rxn_list])

        # Reset module level maxproc back to default
        rmgpy.rmg.main.maxproc = 1

    def teardown_class(self):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None
