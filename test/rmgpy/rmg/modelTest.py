#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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
import pytest

from rmgpy import settings
from rmgpy.data.base import ForbiddenStructures
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.thermo import NASA, NASAPolynomial
from rmgpy.molecule import Molecule
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.rmg.react import react
from rmgpy.species import Species

THERMO_DICT = {
    "[OH]": NASA(
        polynomials=[
            NASAPolynomial(
                coeffs=[
                    3.51457,
                    2.92787e-05,
                    -5.32168e-07,
                    1.0195e-09,
                    -3.85947e-13,
                    3414.25,
                    2.10435,
                ],
                Tmin=(100, "K"),
                Tmax=(1145.75, "K"),
            ),
            NASAPolynomial(
                coeffs=[
                    3.07194,
                    0.000604014,
                    -1.39775e-08,
                    -2.13448e-11,
                    2.48067e-15,
                    3579.39,
                    4.578,
                ],
                Tmin=(1145.75, "K"),
                Tmax=(5000, "K"),
            ),
        ],
        Tmin=(100, "K"),
        Tmax=(5000, "K"),
        E0=(28.3945, "kJ/mol"),
        Cp0=(29.1007, "J/(mol*K)"),
        CpInf=(37.4151, "J/(mol*K)"),
        label="""OH(D)""",
        comment="""Thermo library: primaryThermoLibrary""",
    ),
    "C": NASA(
        polynomials=[
            NASAPolynomial(
                coeffs=[
                    4.20541,
                    -0.00535551,
                    2.51121e-05,
                    -2.1376e-08,
                    5.97513e-12,
                    -10161.9,
                    -0.921259,
                ],
                Tmin=(100, "K"),
                Tmax=(1084.13, "K"),
            ),
            NASAPolynomial(
                coeffs=[
                    0.908298,
                    0.011454,
                    -4.57171e-06,
                    8.29185e-10,
                    -5.66309e-14,
                    -9719.99,
                    13.9929,
                ],
                Tmin=(1084.13, "K"),
                Tmax=(5000, "K"),
            ),
        ],
        Tmin=(100, "K"),
        Tmax=(5000, "K"),
        E0=(-84.435, "kJ/mol"),
        Cp0=(33.2579, "J/(mol*K)"),
        CpInf=(108.088, "J/(mol*K)"),
        label="""CH4""",
        comment="""Thermo library: primaryThermoLibrary""",
    ),
    "[CH3]": NASA(
        polynomials=[
            NASAPolynomial(
                coeffs=[
                    3.67359,
                    0.00201095,
                    5.73022e-06,
                    -6.87117e-09,
                    2.54386e-12,
                    16445,
                    1.60456,
                ],
                Tmin=(200, "K"),
                Tmax=(1000, "K"),
            ),
            NASAPolynomial(
                coeffs=[
                    2.28572,
                    0.0072399,
                    -2.98714e-06,
                    5.95685e-10,
                    -4.67154e-14,
                    16775.6,
                    8.48007,
                ],
                Tmin=(1000, "K"),
                Tmax=(3500, "K"),
            ),
        ],
        Tmin=(200, "K"),
        Tmax=(3500, "K"),
        E0=(136.42, "kJ/mol"),
        Cp0=(33.2579, "J/(mol*K)"),
        CpInf=(83.1447, "J/(mol*K)"),
        label="""CH3""",
        comment="""Thermo library: GRI-Mech3.0""",
    ),
    "[CH2]": NASA(
        polynomials=[
            NASAPolynomial(
                coeffs=[
                    4.01192,
                    -0.000154978,
                    3.26298e-06,
                    -2.40422e-09,
                    5.69497e-13,
                    45867.7,
                    0.533201,
                ],
                Tmin=(100, "K"),
                Tmax=(1104.62, "K"),
            ),
            NASAPolynomial(
                coeffs=[
                    3.14983,
                    0.00296674,
                    -9.76056e-07,
                    1.54115e-10,
                    -9.50338e-15,
                    46058.1,
                    4.77808,
                ],
                Tmin=(1104.62, "K"),
                Tmax=(5000, "K"),
            ),
        ],
        Tmin=(100, "K"),
        Tmax=(5000, "K"),
        E0=(381.37, "kJ/mol"),
        Cp0=(33.2579, "J/(mol*K)"),
        CpInf=(58.2013, "J/(mol*K)"),
        label="""CH2(T)""",
        comment="""Thermo library: primaryThermoLibrary""",
    ),
    "O": NASA(
        polynomials=[
            NASAPolynomial(
                coeffs=[
                    4.05764,
                    -0.000787936,
                    2.90877e-06,
                    -1.47519e-09,
                    2.12842e-13,
                    -30281.6,
                    -0.311364,
                ],
                Tmin=(100, "K"),
                Tmax=(1130.24, "K"),
            ),
            NASAPolynomial(
                coeffs=[
                    2.84325,
                    0.00275109,
                    -7.81031e-07,
                    1.07244e-10,
                    -5.79392e-15,
                    -29958.6,
                    5.91042,
                ],
                Tmin=(1130.24, "K"),
                Tmax=(5000, "K"),
            ),
        ],
        Tmin=(100, "K"),
        Tmax=(5000, "K"),
        E0=(-251.755, "kJ/mol"),
        Cp0=(33.2579, "J/(mol*K)"),
        CpInf=(58.2013, "J/(mol*K)"),
        label="""H2O""",
        comment="""Thermo library: primaryThermoLibrary""",
    ),
}
# these species are used later, but their thermo doesn't matter so we use
# [H] as a placeholder
for key in ("[H]", "C=C[CH2]C", "C=C=CC", "C=CCC", "[H][H]", "[CC]", "C[CH2]"):
    THERMO_DICT[key] = THERMO_DICT["O"]


class TestSpecies:
    """
    Contains unit tests of the Species class.
    """

    @classmethod
    def setup_class(cls):
        """
        A method that is run before each unit test in this class.
        """
        # set-up RMG object
        cls.rmg = RMG()

        # load kinetic database and forbidden structures
        cls.rmg.database = RMGDatabase()
        path = os.path.join(settings["database.directory"])

        # forbidden structure loading
        cls.rmg.database.load_thermo(os.path.join(path, "thermo"))

    def test_get_thermo_data(self):
        """
        Test that get_thermo_data method of Species works.
        """
        spc = Species().from_smiles("CCC")

        assert not spc.thermo
        spc.get_thermo_data()
        assert spc.thermo
        thermo = spc.thermo
        spc.get_thermo_data()

        assert id(thermo) == id(spc.thermo)

        spc.thermo = None
        spc.get_thermo_data()
        assert id(thermo) != id(spc.thermo)

    @classmethod
    def tear_down_class(cls):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None


class TestCoreEdgeReactionModel:
    """
    Contains unit tests of the CoreEdgeReactionModel class.
    """

    @classmethod
    def setup_class(cls):
        """
        A method that is run before each unit test in this class.
        """
        test_family = "H_Abstraction"

        # set-up RMG object
        rmg = RMG()

        # load kinetic database and forbidden structures
        rmg.database = RMGDatabase()
        path = os.path.join(settings["test_data.directory"], "testing_database")

        # kinetics family loading
        rmg.database.load_kinetics(
            os.path.join(path, "kinetics"),
            kinetics_families=[test_family],
            reaction_libraries=[],
        )
        # load empty forbidden structures to avoid any dependence on forbidden structures
        # for these tests
        for family in rmg.database.kinetics.families.values():
            family.forbidden = ForbiddenStructures()
        rmg.database.forbidden_structures = ForbiddenStructures()

    def test_add_new_surface_objects(self):
        """
        basic test that surface movement object management works properly
        """

        # create object with ReactionSystem behavior
        class rsys:
            pass

        class item:
            pass

        T = item()
        P = item()
        T.value_si = 1000.0
        P.value_si = 101000.0
        rsys.T = T
        rsys.P = P
        procnum = 2

        cerm = CoreEdgeReactionModel()

        spcA = Species().from_smiles("[OH]")
        spcs = [Species().from_smiles("CC"), Species().from_smiles("[CH3]")]
        spc_tuples = [((spcA, spc), ["H_Abstraction"]) for spc in spcs]

        rxns = list(itertools.chain.from_iterable(react(spc_tuples, procnum)))
        rxns += list(itertools.chain.from_iterable(react([((spcs[0], spcs[1]), ["H_Abstraction"])], procnum)))

        for rxn in rxns:
            cerm.make_new_reaction(rxn, generate_thermo=False, generate_kinetics=False)

        cerm.core.species = [spcA] + spcs

        corerxns = []
        edgerxns = []
        edgespcs = set()
        for rxn in rxns:
            if set(rxn.reactants + rxn.products) <= set(cerm.core.species):
                corerxns.append(rxn)
            else:
                edgespcs |= set(cerm.core.species) - set(rxn.reactants + rxn.products)
                edgerxns.append(rxn)

        cerm.edge.species += list(edgespcs)

        cerm.core.reactions = corerxns
        cerm.edge.reactions = edgerxns

        cerm.surface.species = []
        cerm.surface.reactions = []

        new_surface_reactions = [cerm.edge.reactions[0]]
        new_surface_species = []
        obj = new_surface_reactions

        cerm.add_new_surface_objects(obj, new_surface_species, new_surface_reactions, rsys)

        empty = set()

        assert cerm.new_surface_spcs_add == empty
        assert cerm.new_surface_spcs_loss == empty
        assert cerm.new_surface_rxns_loss == empty
        assert cerm.new_surface_rxns_add == set([cerm.edge.reactions[0]])

    def test_make_new_species(self):
        """
        Test that CoreEdgeReactionModel.make_new_species method correctly stores the unique species.
        """

        # adding 3 unique species:
        cerm = CoreEdgeReactionModel()

        spcs = [
            Species().from_smiles("[OH]"),
            Species().from_smiles("CC"),
            Species().from_smiles("[CH3]"),
        ]

        for spc in spcs:
            cerm.make_new_species(spc)

        assert len(cerm.species_dict) == len(spcs)
        assert len(cerm.index_species_dict) == len(spcs)

        # adding 3 unique, and 1 already existing species:
        cerm = CoreEdgeReactionModel()

        spcs = [
            Species().from_smiles("[OH]"),
            Species().from_smiles("CC"),
            Species().from_smiles("[CH3]"),
            Species().from_smiles("CC"),
        ]  # duplicate species

        for spc in spcs:
            cerm.make_new_species(spc)

        assert len(cerm.species_dict) == len(spcs) - 1
        assert len(cerm.index_species_dict) == len(spcs) - 1

    def test_make_new_species_transfers_gas_volatile_veto(self):
        """
        The durable gas-volatile veto (props key set by the polymer handshake
        on a discrete gas product) must reach the solver-visible Species.
        make_new_species builds/returns that Species, so it must transfer the
        veto (1) onto a NEWLY created Species and (2) onto an already-existing
        Species when a later duplicate carries the veto -- dedup short-circuits
        before the proxy logic, so without an explicit transfer the verdict is
        dropped (the trap that defeated the earlier clear-based fixes).

        RED before the fix: neither the new nor the deduped Species carries the
        veto key.
        """
        from rmgpy.polymer import POLYMER_REFERENCE_STATE_GAS_VETO_KEY as VETO

        # (1) NEW species: incoming molecule carries the veto in props
        cerm = CoreEdgeReactionModel()
        spc = Species().from_smiles("C=C(C)c1ccccc1")
        spc.props[VETO] = True
        made, _ = cerm.make_new_species(spc)
        assert made.props.get(VETO) is True, (
            "make_new_species must transfer the gas veto onto a new Species"
        )

        # (2) EXISTING species: first register WITHOUT veto, then a duplicate
        # WITH veto must stamp the already-existing (solver-visible) Species.
        cerm = CoreEdgeReactionModel()
        first = Species().from_smiles("C=C(C)c1ccccc1")
        made1, is_new1 = cerm.make_new_species(first)
        assert is_new1
        assert made1.props.get(VETO) in (False, None)
        dup = Species().from_smiles("C=C(C)c1ccccc1")
        dup.props[VETO] = True
        made2, is_new2 = cerm.make_new_species(dup)
        assert is_new2 is False, "duplicate should dedup to the existing species"
        assert made2 is made1, "dedup must return the existing Species object"
        assert made1.props.get(VETO) is True, (
            "make_new_species must transfer the veto onto the existing "
            "(deduped) Species -- else the verdict is lost at first-write-wins"
        )

    def test_polymer_never_acquires_gas_veto_via_make_new_species(self):
        """
        Load-bearing invariant guard (code-review IMPORTANT #2): a Polymer
        (a melt chain) must NEVER acquire the durable gas-volatile veto through
        make_new_species. Correctness of the whole veto scheme rests on melt
        chains being Polymers routed to _register_polymer (fingerprint dedup)
        BEFORE the Species veto-transfer logic runs -- so a genuine chain can
        never be false-vetoed and silently dropped from the melt sum. This test
        pins that guarantee: even with a stray veto contaminating the chain's
        constituent molecule, the registered Polymer's props stay veto-free.

        RED if a future refactor removes the Polymer early-return and routes
        chains through the molecule-reading veto-transfer path.
        """
        from rmgpy.polymer import (Polymer,
                                   POLYMER_REFERENCE_STATE_GAS_VETO_KEY as VETO)

        cerm = CoreEdgeReactionModel()
        poly = Polymer(
            label="PS", monomer="[CH2][CH]c1ccccc1",
            end_groups=["[CH3]", "[H]"], cutoff=3,
            Mn=5000.0, Mw=6000.0, initial_mass=1.0,
        )
        # Contaminate the chain's constituent molecule(s) with a stray veto,
        # the exact hazard a molecule-dedup refactor would expose.
        for m in getattr(poly, "molecule", None) or []:
            if not isinstance(getattr(m, "props", None), dict):
                m.props = {}
            m.props[VETO] = True

        made, _ = cerm._register_polymer(poly, generate_thermo=False)
        # Route through the public entry too (Polymers early-return to
        # _register_polymer before any veto logic).
        made2, _ = cerm.make_new_species(poly, generate_thermo=False)
        for obj in (made, made2):
            assert isinstance(obj, Polymer)
            assert obj.props.get(VETO) in (False, None), (
                "a Polymer melt chain must never acquire the gas-volatile veto "
                "through make_new_species / _register_polymer"
            )

    def test_append_unreactive_structure(self):
        """
        Test that CERM.make_new_species correctly recognizes a non-representative resonance structure
        """

        cerm = CoreEdgeReactionModel()

        spcs = [
            Species().from_smiles("CCO"),  # a control species
            Species().from_smiles("[N]=O"),
            Species().from_adjacency_list(
                """1 O u1 p2 c0 {2,S}
                                               2 N u0 p2 c0 {1,S}"""
            ),  # a non-representative structure of '[N]=O'
        ]

        for spc in spcs:
            cerm.make_new_species(spc)

        assert len(cerm.species_dict) == 2
        assert len(cerm.index_species_dict) == 2
        assert len(cerm.index_species_dict[1].molecule) == 1
        assert cerm.index_species_dict[1].molecule[0].reactive
        assert len(cerm.index_species_dict[2].molecule) == 1
        assert cerm.index_species_dict[2].molecule[0].reactive

    def test_make_new_reaction(self):
        """
        Test that CoreEdgeReactionModel.make_new_reaction method correctly works.
        """

        procnum = 2
        spcA = Species().from_smiles("[OH]")
        spcs = [Species().from_smiles("CC"), Species().from_smiles("[CH3]")]
        spc_tuples = [((spcA, spc), ["H_Abstraction"]) for spc in spcs]

        rxns = list(itertools.chain.from_iterable(react(spc_tuples, procnum)))

        cerm = CoreEdgeReactionModel()

        for rxn in rxns:
            cerm.make_new_reaction(rxn, generate_thermo=False, generate_kinetics=False)

        """
        3 expected H-abstraction reactions:
            OH + CC = H2O + C[CH2]
            OH + [CH3] = H2O + [CH2]
            OH + [CH3] = [O] + C
        """

        # count no. of entries in reactionDict:
        counter = 0
        for fam, v1 in cerm.reaction_dict.items():
            for key2, v2 in v1.items():
                for key3, rxnList in v2.items():
                    counter += len(rxnList)

        assert counter == 3

    def test_thermo_filter_species(self):
        """
        test that thermo_filter_species leaves species alone if if toleranceThermoKeepInEdge
        is high and removes them if if toleranceThermoKeepInEdge is low
        """

        cerm = CoreEdgeReactionModel()

        spcs = [
            Species().from_smiles("[OH]"),
            Species().from_smiles("C"),
            Species().from_smiles("[CH3]"),
            Species().from_smiles("[CH2]"),
            Species().from_smiles("O"),
        ]

        for spc in spcs:
            spc.thermo = THERMO_DICT[spc.molecule[0].to_smiles()]

        for spc in spcs:
            cerm.make_new_species(spc, label=spc.molecule[0].to_smiles())
            spc.label = spc.molecule[0].to_smiles()

        for spc in spcs[:3]:
            cerm.add_species_to_core(spc)

        reaction = TemplateReaction(
            reactants=[spcs[0], spcs[2]],
            products=[spcs[-1], spcs[-2]],
            degeneracy=1,
            reversible=True,
            family="H_Abstraction",
        )

        cerm.process_new_reactions(
            new_reactions=[reaction],
            new_species=[],
            generate_kinetics=False,
            generate_thermo=False,
        )  # adds CH2 and O to edge

        for spc in cerm.core.species + cerm.edge.species:
            spc.thermo = THERMO_DICT[spc.molecule[0].to_smiles()]  # assign thermo

        cerm.set_thermodynamic_filtering_parameters(
            Tmax=300.0,
            thermo_tol_keep_spc_in_edge=1000.0,
            min_core_size_for_prune=0,
            maximum_edge_species=1,
            reaction_systems=[],
        )

        cerm.thermo_filter_species(cerm.edge.species)  # should not do anythinb because toleranceThermoKeepSpeciesInEdge is high

        difset = set([x.molecule[0].to_smiles() for x in cerm.edge.species]) - set([x.molecule[0].to_smiles() for x in cerm.core.species])

        assert len(difset) == 2  # no change in edge

        cerm.set_thermodynamic_filtering_parameters(
            Tmax=300.0,
            thermo_tol_keep_spc_in_edge=0.0,
            min_core_size_for_prune=0,
            maximum_edge_species=1,
            reaction_systems=[],
        )

        cerm.thermo_filter_species(cerm.edge.species)  # should remove stuff since CH2 and O have high thermo

        difset = set([x.molecule[0].to_smiles() for x in cerm.edge.species]) - set([x.molecule[0].to_smiles() for x in cerm.core.species])

        assert len(difset) < 2  # edge is smaller

    def test_thermo_filter_down(self):
        """
        test that thermo_filter_down with maximum_edge_species = 1 reduces
        the edge to one species
        """
        cerm = CoreEdgeReactionModel()

        spcs = [
            Species().from_smiles("[OH]"),
            Species().from_smiles("C"),
            Species().from_smiles("[CH3]"),
            Species().from_smiles("[CH2]"),
            Species().from_smiles("O"),
        ]

        for spc in spcs:
            spc.thermo = THERMO_DICT[spc.molecule[0].to_smiles()]

        for spc in spcs:
            cerm.make_new_species(spc, label=spc.molecule[0].to_smiles())
            spc.label = spc.molecule[0].to_smiles()

        for spc in spcs[:3]:
            cerm.add_species_to_core(spc)

        reaction = TemplateReaction(
            reactants=[spcs[0], spcs[2]],
            products=[spcs[-1], spcs[-2]],
            degeneracy=1,
            reversible=True,
            family="H_Abstraction",
        )

        cerm.process_new_reactions(
            new_reactions=[reaction],
            new_species=[],
            generate_thermo=False,
            generate_kinetics=False,
        )  # add CH2 and O to edge

        for spc in cerm.core.species + cerm.edge.species:
            spc.thermo = THERMO_DICT[spc.molecule[0].to_smiles()]  # assign thermo

        cerm.set_thermodynamic_filtering_parameters(
            Tmax=300.0,
            thermo_tol_keep_spc_in_edge=1000.0,
            min_core_size_for_prune=0,
            maximum_edge_species=1,
            reaction_systems=[],
        )

        difset = set([x.molecule[0].to_smiles() for x in cerm.edge.species]) - set([x.molecule[0].to_smiles() for x in cerm.core.species])

        assert len(difset) == 2  # no change because toleranceThermoKeepSpeciesInEdge is high

        cerm.thermo_filter_down(maximum_edge_species=1)

        difset = set([x.molecule[0].to_smiles() for x in cerm.edge.species]) - set([x.molecule[0].to_smiles() for x in cerm.core.species])

        assert len(difset) == 1  # should be one because we thermo filtered down to one edge species

    def test_check_for_existing_reaction_eliminates_identical_reactions(self):
        """
        Test that check_for_existing_reaction catches identical reactions.
        """
        cerm = CoreEdgeReactionModel()

        # make species' objects
        spcA = Species().from_smiles("[H]")
        spcB = Species().from_smiles("C=C[CH2]C")
        spcC = Species().from_smiles("C=C=CC")
        spcD = Species().from_smiles("[H][H]")
        spcA.label = "[H]"
        spcB.label = "C=C[CH2]C"
        spcC.label = "C=C=CC"
        spcD.label = "[H][H]"
        spcB.generate_resonance_structures()

        for spc in (spcA, spcB, spcC, spcD):
            spc.thermo = THERMO_DICT[spc.molecule[0].to_smiles()]

        cerm.add_species_to_core(spcA)
        cerm.add_species_to_core(spcB)
        cerm.add_species_to_core(spcC)
        cerm.add_species_to_core(spcD)

        reaction_in_model = TemplateReaction(
            reactants=[spcA, spcB],
            products=[spcC, spcD],
            family="H_Abstraction",
            template=["Csd", "H"],
        )
        reaction_in_model.reactants.sort()
        reaction_in_model.products.sort()

        reaction_to_add = TemplateReaction(
            reactants=[spcA, spcB],
            products=[spcC, spcD],
            family="H_Abstraction",
            template=["Csd", "H"],
        )
        cerm.add_reaction_to_core(reaction_in_model)
        cerm.register_reaction(reaction_in_model)

        found, rxn = cerm.check_for_existing_reaction(reaction_to_add)

        assert found, "check_for_existing_reaction failed to identify existing reaction"

    def test_check_for_existing_reaction_keeps_identical_reactions_with_duplicate_flag(
        self,
    ):
        """
        Test that check_for_existing_reaction keeps reactions with different templates and duplicate=True.
        """
        cerm = CoreEdgeReactionModel()

        # make species' objects
        spcA = Species().from_smiles("[H]")
        spcB = Species().from_smiles("C=C[CH2]C")
        spcC = Species().from_smiles("C=C=CC")
        spcD = Species().from_smiles("[H][H]")

        for spc in (spcA, spcB, spcC, spcD):
            spc.thermo = THERMO_DICT[spc.molecule[0].to_smiles()]

        spcA.label = "[H]"
        spcB.label = "C=C[CH2]C"
        spcC.label = "C=C=CC"
        spcD.label = "[H][H]"
        spcB.generate_resonance_structures()

        cerm.add_species_to_core(spcA)
        cerm.add_species_to_core(spcB)
        cerm.add_species_to_core(spcC)
        cerm.add_species_to_core(spcD)

        reaction_in_model = TemplateReaction(
            reactants=[spcA, spcB],
            products=[spcC, spcD],
            family="H_Abstraction",
            template=["Csd", "H"],
            duplicate=True,
        )
        reaction_in_model.reactants.sort()
        reaction_in_model.products.sort()

        reaction_to_add = TemplateReaction(
            reactants=[spcA, spcB],
            products=[spcC, spcD],
            family="H_Abstraction",
            template=["Cs12345", "H"],
            duplicate=True,
        )
        cerm.add_reaction_to_core(reaction_in_model)
        cerm.register_reaction(reaction_in_model)

        found, rxn = cerm.check_for_existing_reaction(reaction_to_add)

        assert not found, "check_for_existing_reaction failed to identify duplicate template reactions"

    def test_check_for_existing_reaction_eliminates_identical_reactions_without_duplicate_flag(
        self,
    ):
        """
        Test that check_for_existing_reaction eliminates reactions with different templates and duplicate=false
        """
        cerm = CoreEdgeReactionModel()

        # make species' objects
        spcA = Species().from_smiles("[H]")
        spcB = Species().from_smiles("C=C[CH2]C")
        spcC = Species().from_smiles("C=C=CC")
        spcD = Species().from_smiles("[H][H]")

        for spc in (spcA, spcB, spcC, spcD):
            spc.thermo = THERMO_DICT[spc.molecule[0].to_smiles()]
        spcA.label = "[H]"
        spcB.label = "C=C[CH2]C"
        spcC.label = "C=C=CC"
        spcD.label = "[H][H]"
        spcB.generate_resonance_structures()

        cerm.add_species_to_core(spcA)
        cerm.add_species_to_core(spcB)
        cerm.add_species_to_core(spcC)
        cerm.add_species_to_core(spcD)

        reaction_in_model = TemplateReaction(
            reactants=[spcA, spcB],
            products=[spcC, spcD],
            family="H_Abstraction",
            template=["Csd", "H"],
            duplicate=False,
        )
        reaction_in_model.reactants.sort()
        reaction_in_model.products.sort()

        reaction_to_add = TemplateReaction(
            reactants=[spcA, spcB],
            products=[spcC, spcD],
            family="H_Abstraction",
            template=["Cs12345", "H"],
            duplicate=False,
        )
        cerm.add_reaction_to_core(reaction_in_model)
        cerm.register_reaction(reaction_in_model)

        found, rxn = cerm.check_for_existing_reaction(reaction_to_add)

        assert found, "check_for_existing_reaction failed to eliminate reactions without duplicate tag"

    def test_check_for_existing_reaction_removes_duplicates_in_opposite_directions(
        self,
    ):
        """
        Test that check_for_existing_reaction removes duplicate reverse reactions
        """
        cerm = CoreEdgeReactionModel()

        # make species' objects
        s1 = Species().from_smiles("[H]")
        s2 = Species().from_smiles("CC")
        s3 = Species().from_smiles("[H][H]")
        s4 = Species().from_smiles("C[CH2]")
        s1.label = "H"
        s2.label = "CC"
        s3.label = "HH"
        s4.label = "C[CH2]"

        rxn_f = TemplateReaction(
            reactants=[s1, s2],
            products=[s3, s4],
            template=["C/H3/Cs/H3", "H_rad"],
            degeneracy=6,
            family="H_Abstraction",
            reverse=TemplateReaction(
                reactants=[s3, s4],
                products=[s1, s2],
                template=["H2", "C_rad/H2/Cs/H3"],
                degeneracy=2,
                family="H_Abstraction",
            ),
        )

        rxn_r = TemplateReaction(
            reactants=[s3, s4],
            products=[s1, s2],
            template=["H2", "C_rad/H2/Cs/H3"],
            degeneracy=2,
            family="H_Abstraction",
            reverse=TemplateReaction(
                reactants=[s1, s2],
                products=[s3, s4],
                template=["C/H3/Cs/H3", "H_rad"],
                degeneracy=6,
                family="H_Abstraction",
            ),
        )

        rxn_f.reactants.sort()
        rxn_f.products.sort()

        cerm.add_reaction_to_core(rxn_f)
        cerm.register_reaction(rxn_f)

        reactions = cerm.search_retrieve_reactions(rxn_r)
        assert 1 == len(reactions), "cerm.search_retrieve_reactions could not identify reverse reaction"

        found, rxn = cerm.check_for_existing_reaction(rxn_r)

        assert found, "check_for_existing_reaction failed to identify existing reaction in the reverse direction"
        assert rxn == rxn_f

    @classmethod
    def teardown_class(cls):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None


@pytest.mark.functional
class TestEnlarge:
    """
    Contains unit tests for CoreEdgeReactionModel.enlarge.
    """

    @classmethod
    def setup_class(cls):
        """
        A method that is run ONCE before all unit tests in this class.
        """
        cls.dirname = os.path.abspath(os.path.join(os.path.dirname(__file__), "temp"))
        os.makedirs(os.path.join(cls.dirname, "pdep"))

        test_family = "R_Recombination"

        cls.rmg = RMG()

        from rmgpy.rmg.input import set_global_rmg, pressure_dependence

        set_global_rmg(cls.rmg)

        pressure_dependence(
            method="modified strong collision",
            maximumGrainSize=(0.5, "kcal/mol"),
            minimumNumberOfGrains=250,
            temperatures=(300, 2100, "K", 8),
            pressures=(0.1, 100, "bar", 5),
            interpolation=("Chebyshev", 6, 4),
            maximumAtoms=10,
        )

        cls.rmg.output_directory = cls.rmg.pressure_dependence.output_file = cls.dirname

        cls.rmg.database = RMGDatabase()
        cls.rmg.database.load(
            path=settings["database.directory"],
            thermo_libraries=["primaryThermoLibrary"],
            kinetics_families=[test_family],
            reaction_libraries=[],
        )

        cls.rmg.reaction_model = CoreEdgeReactionModel()
        cls.rmg.reaction_model.pressure_dependence = cls.rmg.pressure_dependence

    def test_enlarge_1_add_nonreactive_species(self):
        """Test that we can add a nonreactive species to CERM"""
        m0 = Molecule(smiles="[He]")
        spc0 = self.rmg.reaction_model.make_new_species(m0, label="He", reactive=False)[0]
        self.rmg.reaction_model.enlarge(spc0)

        assert len(self.rmg.reaction_model.core.species) == 1
        assert not self.rmg.reaction_model.core.species[0].reactive

    def test_enlarge_2_add_reactive_species(self):
        """Test that we can add reactive species to CERM"""
        m1 = Molecule(smiles="CC")
        spc1 = self.rmg.reaction_model.make_new_species(m1, label="C2H4")[0]
        self.rmg.reaction_model.enlarge(spc1)

        assert len(self.rmg.reaction_model.core.species) == 2
        assert self.rmg.reaction_model.core.species[1].reactive

        m2 = Molecule(smiles="[CH3]")
        spc2 = self.rmg.reaction_model.make_new_species(m2, label="CH3")[0]
        self.rmg.reaction_model.enlarge(spc2)

        assert len(self.rmg.reaction_model.core.species) == 3
        assert self.rmg.reaction_model.core.species[2].reactive

    def test_enlarge_3_react_edge(self):
        """Test that enlarge properly generated reactions"""
        self.rmg.reaction_model.enlarge(
            react_edge=True,
            unimolecular_react=np.array([0, 1, 0], bool),
            bimolecular_react=np.zeros((3, 3), bool),
        )

        assert len(self.rmg.reaction_model.edge.species) == 2
        smiles = set([spc.smiles for spc in self.rmg.reaction_model.edge.species])
        assert smiles == {"[H]", "C[CH2]"}

        # We expect either C-C bond scission to be in the core and C-H bond scission to be in the edge
        assert len(self.rmg.reaction_model.core.reactions) == 1
        rxn = self.rmg.reaction_model.core.reactions[0]
        smiles = set([spc.smiles for spc in rxn.reactants + rxn.products])
        assert smiles == {"CC", "[CH3]"}

        assert len(self.rmg.reaction_model.edge.reactions) == 1
        rxn = self.rmg.reaction_model.edge.reactions[0]
        smiles = set([spc.smiles for spc in rxn.reactants + rxn.products])
        assert smiles == {"CC", "C[CH2]", "[H]"}

    def test_enlarge_4_create_pdep_network(self):
        """Test that enlarge properly creates a pdep network"""
        assert len(self.rmg.reaction_model.network_list) == 1
        assert len(self.rmg.reaction_model.network_list[0].source) == 1
        assert self.rmg.reaction_model.network_list[0].source[0].label == "C2H4"

        assert len(self.rmg.reaction_model.network_dict) == 1
        assert len(list(self.rmg.reaction_model.network_dict.keys())[0]) == 1
        assert list(self.rmg.reaction_model.network_dict.keys())[0][0].label == "C2H4"

    @classmethod
    def teardown_class(cls):
        """
        A method that is run ONCE after all unit tests in this class.

        Clear global variables and clean up files.
        """
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None
        import rmgpy.rmg.input

        rmgpy.rmg.input.rmg = None
        import shutil

        shutil.rmtree(cls.dirname)


# ---------------------------------------------------------------------------
# Helper and test class for Task 2: _convert_bm_kinetics_with_dHrxn
# ---------------------------------------------------------------------------

import copy

import numpy as np
from rmgpy.thermo import ThermoData
from rmgpy.kinetics import ArrheniusBM
from rmgpy.reaction import Reaction


def _spc_with_h298(smiles, h298_kJ):
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.thermo = ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
        Cpdata=([50, 60, 70, 80, 95, 105, 120], "J/(mol*K)"),
        H298=(h298_kJ, "kJ/mol"), S298=(300.0, "J/(mol*K)"),
        Cp0=(33.2578, "J/(mol*K)"), CpInf=(232.805, "J/(mol*K)"))
    return s


class TestPolymerPoolRegistryDedup:
    def test_spawn_pass_registry_dedups_freshly_promoted_daughter(self, monkeypatch):
        """Regression: a freshly-promoted daughter Polymer sits in BOTH
        core.species and new_species_list until the next enlarge clears it,
        so the pool_registry built by _apply_multipool_spawn_pass (core +
        edge + new_species_list) contained the SAME object twice (observed
        live as duplicated sidecar pools [PS, tail, tail_2, tail, tail_2]).
        The registry must be identity-deduped, order-preserving.

        RED before the fix: the captured registry has 2 entries."""
        import rmgpy.polymer as rmgpy_polymer
        from rmgpy.polymer import Polymer

        cerm = CoreEdgeReactionModel()
        poly = Polymer(
            label="PS", monomer="[CH2][CH]c1ccccc1",
            end_groups=["[CH3]", "[H]"], cutoff=3,
            Mn=5000.0, Mw=6000.0, initial_mass=1.0,
        )
        # The both-lists condition of a freshly-promoted daughter:
        cerm.core.species.append(poly)
        cerm.new_species_list.append(poly)

        captured = {}

        def fake_multipool(*args, **kwargs):
            captured["registry"] = list(kwargs["pool_registry"])
            return [], []  # (processed, no spawn intents)

        monkeypatch.setattr(rmgpy_polymer,
                            "process_polymer_candidates_multipool",
                            fake_multipool)

        proxy = Species().from_smiles("C=Cc1ccccc1")
        proxy.is_polymer_proxy = True
        cerm._apply_multipool_spawn_pass([proxy])

        registry = captured["registry"]
        assert len(registry) == 1, (
            f"pool_registry must dedup the same Polymer object appearing in "
            f"both core.species and new_species_list; got {len(registry)} "
            f"entries: {[p.label for p in registry]}"
        )
        assert registry[0] is poly


class TestPolymerBMRealDHrxn:
    def test_convert_bm_matches_fix_barrier_height(self):
        cerm = CoreEdgeReactionModel()
        r = _spc_with_h298("CCO", 50.0)
        p1 = _spc_with_h298("C=CO", 120.0)
        p2 = _spc_with_h298("[H][H]", 0.0)
        bm = ArrheniusBM(A=(1.293332e12, "s^-1"), n=0.0,
                         w0=(968.0, "kJ/mol"), E0=(182.946, "kJ/mol"))
        rxn_ours = Reaction(reactants=[r], products=[p1, p2], kinetics=copy.deepcopy(bm))
        rxn_ref = Reaction(reactants=[r], products=[p1, p2], kinetics=copy.deepcopy(bm))
        dH = rxn_ref.get_enthalpy_of_reaction(298)
        cerm._convert_bm_kinetics_with_dHrxn(rxn_ours, dH)
        rxn_ref.fix_barrier_height()
        from rmgpy.kinetics import Arrhenius
        assert isinstance(rxn_ours.kinetics, Arrhenius)
        assert np.isclose(rxn_ours.kinetics.A.value_si, rxn_ref.kinetics.A.value_si, rtol=1e-12)
        assert np.isclose(rxn_ours.kinetics.n.value_si, rxn_ref.kinetics.n.value_si, rtol=1e-12)
        assert np.isclose(rxn_ours.kinetics.Ea.value_si, rxn_ref.kinetics.Ea.value_si, rtol=1e-10)

    def test_real_dHrxn_reuses_existing_species_thermo(self):
        """When snapshot products are Species carrying thermo, _polymer_real_dHrxn
        reuses it (no estimation) -> deterministic ΔH, DB-free."""
        cerm = CoreEdgeReactionModel()
        react = _spc_with_h298("CCO", 50.0)            # reactant H298 = 50 kJ/mol
        prod1 = _spc_with_h298("C=CO", 120.0)          # real product H298 = 120 kJ/mol
        prod2 = _spc_with_h298("[H][H]", 0.0)          # real product H298 = 0 kJ/mol
        dH = cerm._polymer_real_dHrxn([react], [prod1, prod2])
        # ΔH = (120 + 0) - 50 = 70 kJ/mol = 70000 J/mol
        assert np.isclose(dH, 70000.0, rtol=1e-9)


# ---------------------------------------------------------------------------
# Task 3: _thermo_for_snapshot_product estimation path (DB-backed)
# ---------------------------------------------------------------------------


class TestPolymerRealDHrxnEstimation:
    @classmethod
    def setup_class(cls):
        import os
        from rmgpy import settings
        from rmgpy.data.rmg import RMGDatabase
        from rmgpy.rmg.main import RMG
        rmg = RMG()
        rmg.database = RMGDatabase()
        rmg.database.load_thermo(os.path.join(settings["database.directory"], "thermo"))
        cls.cerm = CoreEdgeReactionModel()

    def test_estimation_path_matches_independent_estimate(self):
        from rmgpy.molecule import Molecule
        # real product as a bare Molecule (no thermo) -> forces estimation
        styrene = Molecule().from_smiles("C=Cc1ccccc1")
        snap = [styrene.copy(deep=True)]
        spc = self.cerm._thermo_for_snapshot_product(snap[0])
        assert spc.has_thermo()
        # independent estimate of the same molecule via the model thermo path
        ref = Species(molecule=[Molecule().from_smiles("C=Cc1ccccc1")])
        ref.generate_resonance_structures()
        self.cerm.generate_thermo(ref)
        assert np.isclose(spc.get_enthalpy(298), ref.get_enthalpy(298), rtol=1e-6)


class TestMakeNewReactionBMGate:
    """
    Guard/skip tests: the BM pre-conversion gate must NOT over-fire.

    (c) Non-polymer reaction with ArrheniusBM: polymer_reactants is empty,
        so real_products_snapshot / relabeled stay at defaults; our block
        never runs.  Spy count == 0.

    (d1) Polymer reaction where _handshake_structures returns False (small
         gas products CO2 that cannot become polymer objects): relabeled
         stays False; gate skips.  Spy count == 0.
    """

    @classmethod
    def setup_class(cls):
        import os
        from rmgpy import settings
        from rmgpy.data.rmg import RMGDatabase
        from rmgpy.data.thermo import NASA, NASAPolynomial
        from rmgpy.rmg.main import RMG
        from rmgpy.polymer import Polymer

        rmg = RMG()
        rmg.database = RMGDatabase()
        rmg.database.load_thermo(os.path.join(settings["database.directory"], "thermo"))

        cls.cerm = CoreEdgeReactionModel()
        cls.ps = Polymer(
            label='PSGate',
            monomer='[CH2][CH]c1ccccc1',
            end_groups=['[CH3]', '[H]'],
            cutoff=3,
            Mn=5000.0,
            Mw=6000.0,
            initial_mass=1.0,
        )
        cls.cerm._register_polymer(cls.ps, generate_thermo=False)
        # Inject fake thermo so fix_barrier_height can call get_enthalpy on the
        # polymer reactant (Polymer objects have no real thermo by default).
        fake_thermo = NASA(
            polynomials=[
                NASAPolynomial(
                    coeffs=[3.5, 0.0, 0.0, 0.0, 0.0, -1000.0, 5.0],
                    Tmin=(200.0, 'K'), Tmax=(1000.0, 'K')
                ),
                NASAPolynomial(
                    coeffs=[3.5, 0.0, 0.0, 0.0, 0.0, -1000.0, 5.0],
                    Tmin=(1000.0, 'K'), Tmax=(6000.0, 'K')
                ),
            ],
            Tmin=(200.0, 'K'),
            Tmax=(6000.0, 'K'),
        )
        cls.ps.thermo = fake_thermo

    def test_non_polymer_bm_reaction_untouched(self):
        """(c) Non-polymer ArrheniusBM reaction: pre-conversion must not fire."""
        from rmgpy.kinetics import ArrheniusBM

        # Install spy
        called = {"n": 0}
        orig = self.cerm._convert_bm_kinetics_with_dHrxn
        def spy(*args, **kwargs):
            called["n"] += 1
            return orig(*args, **kwargs)
        self.cerm._convert_bm_kinetics_with_dHrxn = spy

        try:
            rxn = TemplateReaction(
                reactants=[Molecule().from_smiles('CCO')],
                products=[Molecule().from_smiles('C=CO'), Molecule().from_smiles('[H][H]')],
                family='Retroene',
                is_forward=True,
                kinetics=ArrheniusBM(
                    A=(1.0e12, 's^-1'), n=0.0,
                    w0=(968.0, 'kJ/mol'), E0=(150.0, 'kJ/mol'),
                ),
            )
            result, is_new = self.cerm.make_new_reaction(
                rxn,
                check_existing=False,
                generate_thermo=True,
                generate_kinetics=True,
            )
            assert result is not None, "make_new_reaction returned None unexpectedly"
            assert called["n"] == 0, (
                f"_convert_bm_kinetics_with_dHrxn called {called['n']} time(s); "
                "expected 0 for non-polymer reaction"
            )
        finally:
            self.cerm._convert_bm_kinetics_with_dHrxn = orig

    def test_unrelabeled_polymer_bm_stays_on_generic_path(self):
        """(d1) Polymer reaction where handshake does NOT relabel: pre-conversion skipped."""
        from rmgpy.kinetics import ArrheniusBM
        from rmgpy.polymer import Polymer

        # Fresh model so state from (c) doesn't interfere
        cerm2 = CoreEdgeReactionModel()
        ps2 = Polymer(
            label='PSGate2',
            monomer='[CH2][CH]c1ccccc1',
            end_groups=['[CH3]', '[H]'],
            cutoff=3,
            Mn=5000.0,
            Mw=6000.0,
            initial_mass=1.0,
        )
        cerm2._register_polymer(ps2, generate_thermo=False)
        # Inject fake thermo so fix_barrier_height can get enthalpy
        from rmgpy.data.thermo import NASA, NASAPolynomial
        fake_thermo = NASA(
            polynomials=[
                NASAPolynomial(
                    coeffs=[3.5, 0.0, 0.0, 0.0, 0.0, -1000.0, 5.0],
                    Tmin=(200.0, 'K'), Tmax=(1000.0, 'K')
                ),
                NASAPolynomial(
                    coeffs=[3.5, 0.0, 0.0, 0.0, 0.0, -1000.0, 5.0],
                    Tmin=(1000.0, 'K'), Tmax=(6000.0, 'K')
                ),
            ],
            Tmin=(200.0, 'K'),
            Tmax=(6000.0, 'K'),
        )
        ps2.thermo = fake_thermo

        # Install spy on cerm2
        called = {"n": 0}
        orig = cerm2._convert_bm_kinetics_with_dHrxn
        def spy(*args, **kwargs):
            called["n"] += 1
            return orig(*args, **kwargs)
        cerm2._convert_bm_kinetics_with_dHrxn = spy

        proxy_mol = ps2.baseline_proxy.molecule[0].copy(deep=True)
        # CO2: _handshake_structures returns False (cannot become a polymer fragment)
        co2 = Molecule().from_smiles('O=C=O')
        rxn = TemplateReaction(
            reactants=[proxy_mol],
            products=[co2],
            family='Retroene',
            is_forward=True,
            kinetics=ArrheniusBM(
                A=(1.0e12, 's^-1'), n=0.0,
                w0=(968.0, 'kJ/mol'), E0=(150.0, 'kJ/mol'),
            ),
        )
        result, is_new = cerm2.make_new_reaction(
            rxn,
            check_existing=False,
            generate_thermo=True,
            generate_kinetics=True,
        )
        assert result is not None, "make_new_reaction returned None unexpectedly"
        assert called["n"] == 0, (
            f"_convert_bm_kinetics_with_dHrxn called {called['n']} time(s); "
            "expected 0 for unrelabeled polymer reaction"
        )
