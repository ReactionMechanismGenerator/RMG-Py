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
import unittest

import numpy as np

from rmgpy import settings
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.kinetics.arrhenius import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.tools.isotopes import correct_entropy, apply_kinetic_isotope_effect_simple, \
    generate_isotope_reactions, get_reduced_mass, get_labeled_reactants, \
    is_enriched, generate_isotopomers, cluster, remove_isotope, \
    redo_isotope, ensure_reaction_direction, compare_isotopomers

database = None


def setUpModule():
    """A function that is run ONCE before all unit tests in this module."""
    global database
    database = RMGDatabase()
    database.load(
        path=os.path.join(settings['test_data.directory'], 'testing_database'),
        kinetics_families=[
            'H_Abstraction', 'intra_H_migration'
        ],
        testing=True,
        depository=False,
        solvation=False,
        surface=False,
    )
    database.load_forbidden_structures()

    # Prepare the database by loading training reactions and averaging the rate rules
    for family in database.kinetics.families.values():
        if not family.auto_generated:
            family.add_rules_from_training(thermo_database=database.thermo)
            family.fill_rules_by_averaging_up(verbose=True)


def tearDownModule():
    """A function that is run ONCE after all unit tests in this module."""
    from rmgpy.data import rmg
    rmg.database = None


class IsotopesTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        global database
        cls.database = database
        cls.family = database.kinetics.families['H_Abstraction']

    def test_cluster_with_species(self):
        """
        Test that isotope partitioning algorithm work with Reaction Objects.
        """

        eth = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")

        ethi = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")

        meth = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
""")

        spc_list = [eth, ethi]

        clusters = cluster(spc_list)
        self.assertEquals(len(clusters), 1)
        self.assertEquals(len(clusters[0]), 2)

        spc_list = [meth, ethi]

        clusters = cluster(spc_list)
        self.assertEquals(len(clusters), 2)
        self.assertEquals(len(clusters[0]), 1)

    def test_cluster_with_reactions(self):
        """
        Test that isotope partitioning algorithm works with Reaction objects
        """

        eth = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")

        ethi = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")

        meth = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
""")

        rxn0 = Reaction(reactants=[ethi, ethi], products=[ethi, eth])
        rxn1 = Reaction(reactants=[eth, ethi], products=[eth, eth])
        rxn2 = Reaction(reactants=[ethi, meth], products=[meth, ethi])
        rxn3 = Reaction(reactants=[eth, meth], products=[eth, meth])
        rxn4 = Reaction(reactants=[meth], products=[meth])
        rxn5 = Reaction(reactants=[ethi], products=[eth])

        same_cluster_list = [rxn0, rxn1]

        clusters = cluster(same_cluster_list)
        self.assertEquals(len(clusters), 1)
        self.assertEquals(len(clusters[0]), 2)

        same_cluster_list = [rxn2, rxn3]

        clusters = cluster(same_cluster_list)
        self.assertEquals(len(clusters), 1)
        self.assertEquals(len(clusters[0]), 2)

        multi_cluster_list = [rxn0, rxn1, rxn2, rxn3, rxn4, rxn5]

        clusters = cluster(multi_cluster_list)
        self.assertEquals(len(clusters), 4)
        self.assertEquals(len(clusters[0]), 1)

    def test_remove_isotope_for_reactions(self):
        """
        Test that remove isotope algorithm works with Reaction objects.
        """

        eth = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")

        ethi = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")
        unlabeled_rxn = Reaction(reactants=[eth], products=[eth])
        labeled_rxn = Reaction(reactants=[ethi], products=[ethi])
        stripped = remove_isotope(labeled_rxn)

        self.assertTrue(unlabeled_rxn.is_isomorphic(stripped))

    def test_remove_isotope_for_species(self):
        """
        Test that remove isotope algorithm works with Species.
        """

        eth = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")

        ethi = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")

        stripped = remove_isotope(ethi)

        self.assertTrue(eth.is_isomorphic(stripped))

    def test_inplace_remove_isotope_for_reactions(self):
        """
        Test that removeIsotope and redoIsotope works with reactions
        """

        eth = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")

        ethi = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")
        unlabeled_rxn = Reaction(reactants=[eth], products=[eth])
        labeled_rxn = Reaction(reactants=[ethi], products=[ethi])
        stored_labeled_rxn = labeled_rxn.copy()
        modified_atoms = remove_isotope(labeled_rxn, inplace=True)

        self.assertTrue(unlabeled_rxn.is_isomorphic(labeled_rxn))

        redo_isotope(modified_atoms)

        self.assertTrue(stored_labeled_rxn.is_isomorphic(labeled_rxn))

    def test_ensure_reaction_direction(self):
        """
        Tests that the direction of the reaction is constant for every isotopomer
        """

        # get reactions

        methyl = Species().from_smiles('[CH3]')
        methyli = Species().from_adjacency_list("""
multiplicity 2
1 C u1 p0 c0 i13 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""")
        methane = Species().from_smiles('C')
        methanei = Species().from_adjacency_list("""
1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
""")

        dipropyl = Species().from_smiles('[CH2]C[CH2]')
        dipropyli = Species().from_adjacency_list("""
multiplicity 3
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 C u0 p0 c0 {1,S} {5,S} {8,S} {9,S}
5 C u1 p0 c0 i13 {4,S} {6,S} {7,S}
6 H u0 p0 c0 {5,S}
7 H u0 p0 c0 {5,S}
8 H u0 p0 c0 {4,S}
9 H u0 p0 c0 {4,S}
""")

        propyl = Species().from_smiles('CC[CH2]')
        propyli = Species().from_adjacency_list("""
multiplicity 2
1  C u1 p0 c0 i13 {2,S} {3,S} {4,S}
2  H u0 p0 c0 {1,S}
3  H u0 p0 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
5  C u0 p0 c0 {4,S} {8,S} {9,S} {10,S}
6  H u0 p0 c0 {4,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {5,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
""")
        propyl.label = 'propyl'
        propyli.label = 'propyli'
        dipropyl.label = 'dipropyl'
        dipropyli.label = 'dipropyli'
        methyl.label = 'methyl'
        methyli.label = 'methyli'
        methane.label = 'methane'
        methanei.label = 'methanei'
        # make reactions
        rxn1 = TemplateReaction(reactants=[dipropyl, methane],
                                products=[methyl, propyl],
                                kinetics=Arrhenius(A=(1., 'cm^3/(mol*s)'), Ea=(2., 'kJ/mol'), n=0.),
                                family='H_Abstraction',
                                template=['a', 'c'],
                                degeneracy=8,
                                pairs=[[methane, methyl], [dipropyl, propyl]])
        rxn2 = TemplateReaction(reactants=[methyli, propyl],
                                products=[methanei, dipropyl],
                                kinetics=Arrhenius(A=(1e-20, 'cm^3/(mol*s)'), Ea=(2., 'kJ/mol'), n=0.),
                                family='H_Abstraction',
                                template=['b', 'd'],
                                degeneracy=3,
                                pairs=[[methyli, methanei], [propyl, dipropyl]])
        rxn3 = TemplateReaction(reactants=[methane, dipropyli],
                                products=[methyl, propyli],
                                kinetics=Arrhenius(A=(0.5, 'cm^3/(mol*s)'), Ea=(2., 'kJ/mol'), n=0.),
                                family='H_Abstraction',
                                template=['a', 'c'],
                                degeneracy=4,
                                pairs=[[methane, methyl], [dipropyli, propyli]])
        rxn4 = TemplateReaction(reactants=[methyli, propyli],
                                products=[methanei, dipropyli],
                                kinetics=Arrhenius(A=(1e-20, 'cm^3/(mol*s)'), Ea=(2., 'kJ/mol'), n=0.),
                                family='H_Abstraction',
                                template=['d', 'b'],
                                degeneracy=3,
                                pairs=[[methyli, methanei], [propyli, dipropyli]])

        rxns = [rxn1, rxn2, rxn3, rxn4]

        # call method
        ensure_reaction_direction(rxns)

        for rxn in rxns:
            # ensure there is a methane in reactants for each reaction
            self.assertTrue(any([compare_isotopomers(methane, reactant) for reactant in rxn.reactants]),
                            msg='ensureReactionDirection didnt flip the proper reactants and products')

            # ensure kinetics is correct
            if any([dipropyli.is_isomorphic(reactant) for reactant in rxn.reactants]):
                self.assertAlmostEqual(rxn.kinetics.A.value, 0.5,
                                       msg='The A value returned, {0}, is incorrect. '
                                           'Check the reactions degeneracy and how A.value is obtained. '
                                           'The reaction is:{1}'.format(rxn.kinetics.A.value, rxn))
            else:
                self.assertAlmostEqual(rxn.kinetics.A.value, 1.,
                                       msg='The A value returned, {0}, is incorrect. '
                                           'Check the reactions degeneracy and how A.value is obtained. '
                                           'The reaction is:{1}'.format(rxn.kinetics.A.value, rxn))

    def test_ensure_reaction_direction_with_multiple_ts(self):
        """Tests that ensure reaction direction can handle multiple transition states"""
        family = self.database.kinetics.families['intra_H_migration']
        r = Molecule().from_smiles("[CH2]CCC")
        p = Molecule().from_smiles("C[CH]CC")
        rxn = TemplateReaction(reactants=[r], products=[p])
        family.add_atom_labels_for_reaction(reaction=rxn)
        rxn.template = family.get_reaction_template_labels(reaction=rxn)
        rxn.degeneracy = family.calculate_degeneracy(rxn)
        rxn.family = 'intra_H_migration'
        rxn.kinetics = Arrhenius(A=(1, 's^-1'))
        ri = Molecule().from_adjacency_list("""
multiplicity 2
1  C u1 p0 c0 {2,S} {3,S} {4,S}
2  H u0 p0 c0 {1,S}
3  H u0 p0 c0 {1,S}
4  C u0 p0 c0 i13 {1,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 i13 {4,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
""")
        pi = Molecule().from_adjacency_list("""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u1 p0 c0 i13 {1,S} {3,S} {4,S}
3  H u0 p0 c0 {2,S}
4  C u0 p0 c0 i13 {2,S} {5,S} {9,S} {10,S}
5  C u0 p0 c0 {4,S} {11,S} {12,S} {13,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {5,S}
13 H u0 p0 c0 {5,S}
""")
        rxni = TemplateReaction(reactants=[pi], products=[ri], pairs=[[pi, ri]])
        family.add_atom_labels_for_reaction(reaction=rxni)
        rxni.template = family.get_reaction_template_labels(reaction=rxni)
        rxn_cluster = [rxn, rxni]
        ensure_reaction_direction(rxn_cluster)

        self.assertEqual(rxn_cluster[0].degeneracy, 2)
        self.assertEqual(rxn_cluster[1].degeneracy, 2)

        self.assertIn('R2Hall', rxn_cluster[0].template)
        self.assertIn('R2Hall', rxn_cluster[1].template)

        self.assertAlmostEqual(rxn_cluster[0].kinetics.get_rate_coefficient(298),
                               rxn_cluster[1].kinetics.get_rate_coefficient(298))

    def test_compare_isotopomers_works_on_species(self):
        """
        Test that compareIsotomers works on species objects
        """
        ethii = Species().from_adjacency_list("""
1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")
        ethi = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")
        self.assertTrue(compare_isotopomers(ethii, ethi))

    def test_compare_isotopomers_does_not_alter_species(self):
        """
        Test that compareIsotomers works on species objects
        """
        ethii = Species().from_adjacency_list("""
1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")
        ethi = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")
        compare_isotopomers(ethii, ethi)

        # ensure species still have labels
        for atom in ethii.molecule[0].atoms:
            if atom.element.symbol == 'C':
                self.assertEqual(atom.element.isotope, 13, 'compareIsotopomer removed the isotope of a species.')

    def test_compare_isotopomers_fails_on_species(self):
        """
        Test that compareIsotomers fails on species objects
        """
        ethane = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")
        ethenei = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 C u0 p0 c0 i13 {1,D} {6,S} {5,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
""")
        self.assertFalse(compare_isotopomers(ethane, ethenei))

    def test_compare_isotopomers_works_on_reactions(self):
        """
        Test that compareIsotomers works on different reaction objects
        """
        h = Species().from_adjacency_list("""
multiplicity 2
1 H u1 p0 c0
""")
        h2 = Species().from_adjacency_list("""
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""")
        propanei = Species().from_adjacency_list("""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 i13 {2,S} {9,S} {10,S} {11,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
""")
        propane = Species().from_adjacency_list("""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
""")
        npropyli = Species().from_adjacency_list("""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u1 p0 c0 i13 {2,S} {4,S} {5,S}
4  H u0 p0 c0 {3,S}
5  H u0 p0 c0 {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
""")
        npropyl = Species().from_adjacency_list("""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u1 p0 c0 {2,S} {4,S} {5,S}
4  H u0 p0 c0 {3,S}
5  H u0 p0 c0 {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
""")

        reaction2 = TemplateReaction(reactants=[propanei, h],
                                     products=[npropyli, h2],
                                     family='H_Abstraction')
        reaction3 = TemplateReaction(reactants=[propane, h],
                                     products=[h2, npropyl],
                                     family='H_Abstraction')
        self.assertTrue(compare_isotopomers(reaction2, reaction3))

    def test_compare_isotopomers_fails_on_reactions(self):
        """
        Test that compareIsotomers fails on different reaction objects
        """
        h = Species().from_adjacency_list("""
multiplicity 2
1 H u1 p0 c0
""")
        h2 = Species().from_adjacency_list("""
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""")
        propanei = Species().from_adjacency_list("""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 i13 {2,S} {9,S} {10,S} {11,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
""")
        propane = Species().from_adjacency_list("""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
""")
        npropyli = Species().from_adjacency_list("""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u1 p0 c0 i13 {2,S} {4,S} {5,S}
4  H u0 p0 c0 {3,S}
5  H u0 p0 c0 {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
""")

        reaction2 = TemplateReaction(reactants=[propanei, h],
                                     products=[npropyli, h2],
                                     family='H_Abstraction')

        magic_reaction = TemplateReaction(reactants=[propane, h],
                                         products=[propanei, h],
                                         family='H_Abstraction')
        self.assertFalse(compare_isotopomers(reaction2, magic_reaction))

    def test_correct_entropy(self):
        """
        Tests that correctEntropy effectively makes the isotopomer have the same
        thermo as isotopeless with the change in entropy corresponding to the
        symmetry difference.

        This example compares propane to a asymmetrically labeled propane.
        """
        from rmgpy.thermo.nasa import NASA, NASAPolynomial
        from copy import deepcopy
        from rmgpy import constants
        propanei = Species().from_adjacency_list("""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 i13 {2,S} {9,S} {10,S} {11,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
""")
        propane = Species().from_adjacency_list("""
1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2  C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3  C u0 p0 c0 {2,S} {9,S} {10,S} {11,S}
4  H u0 p0 c0 {1,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {2,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}
""")

        # get thermo of propane
        original_thermo = NASA(polynomials=[NASAPolynomial(
            coeffs=[1.10564, 0.0315339, -1.48274e-05, 3.39621e-09, -2.97953e-13, -14351.9, 18.775],
            Tmin=(100, 'K'), Tmax=(3370.6, 'K')), NASAPolynomial(
            coeffs=[-1.45473, 0.0241202, -6.87667e-06, 9.03634e-10, -4.48389e-14, -6688.59, 43.0459],
            Tmin=(3370.6, 'K'), Tmax=(5000, 'K'))], Tmin=(100, 'K'), Tmax=(5000, 'K'),
            comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + gauche(Cs(CsCsRR)) + other(R) + group(Cs-CsHHH) + gauche(Cs(Cs(CsRR)RRR)) + other(R) + group(Cs-CsHHH) + gauche(Cs(Cs(CsRR)RRR)) + other(R)""")
        propane.thermo = deepcopy(original_thermo)
        propanei.thermo = deepcopy(original_thermo)

        correct_entropy(propanei, propane)

        self.assertAlmostEqual(propane.get_enthalpy(298), propanei.get_enthalpy(298))
        self.assertAlmostEqual(propanei.get_entropy(298) - propane.get_entropy(298), constants.R * np.log(2))

    def test_generate_isotopomers(self):
        """
        Test that the generation of isotopomers with N isotopes works.
        """
        from rmgpy.thermo.nasa import NASAPolynomial, NASA

        spc = Species().from_smiles('CC')

        polynomial = NASAPolynomial(coeffs=[1., 1., 1., 1., 1., 1., 1.],
                                    Tmin=(200, 'K'), Tmax=(1600, 'K'), E0=(1., 'kJ/mol'),
                                    comment='made up thermo')

        spc.thermo = NASA(polynomials=[polynomial], Tmin=(200, 'K'),
                          Tmax=(1600, 'K'), E0=(1., 'kJ/mol'),
                          comment='made up thermo')

        spcs = generate_isotopomers(spc, 0)
        self.assertEquals(len(spcs), 0)

        spcs = generate_isotopomers(spc)
        self.assertEquals(len(spcs), 1)

        spcs = generate_isotopomers(spc, 2)
        self.assertEquals(len(spcs), 2)

        spcs = generate_isotopomers(spc, 3)
        self.assertEquals(len(spcs), 2)

    def test_is_enriched(self):
        """
        ensures the Enriched method functions
        """
        npropyl = Species().from_adjacency_list("""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3  C u1 p0 c0 {2,S} {4,S} {5,S}
4  H u0 p0 c0 {3,S}
5  H u0 p0 c0 {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
""")
        npropyli = Species().from_adjacency_list("""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 i13 {1,S} {3,S} {9,S} {10,S}
3  C u1 p0 c0 {2,S} {4,S} {5,S}
4  H u0 p0 c0 {3,S}
5  H u0 p0 c0 {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}
""")

        self.assertTrue(is_enriched(npropyli))
        self.assertFalse(is_enriched(npropyl))

        enriched_reaction = TemplateReaction(reactants=[npropyl],
                                            products=[npropyli],
                                            family='H_Abstraction')
        self.assertTrue(is_enriched(enriched_reaction))

        bare_reaction = TemplateReaction(reactants=[npropyl],
                                        products=[npropyl],
                                        family='H_Abstraction')

        self.assertFalse(is_enriched(bare_reaction))

    def test_get_labeled_reactants(self):
        """
        tests to ensure that get_labeled_reactants returns labeled reactants
        """
        reactant_pair = [Species().from_smiles("C"), Species().from_smiles("[H]")]
        product_pair = [Species().from_smiles("[H][H]"), Species().from_smiles("[CH3]")]
        rxn = TemplateReaction(reactants=reactant_pair,
                               products=product_pair,
                               family='H_Abstraction')
        labeled_reactants = get_labeled_reactants(rxn, self.family)
        r1_labels = labeled_reactants[0].get_all_labeled_atoms()
        self.assertIn("*1", list(r1_labels.keys()))
        self.assertIn("*2", list(r1_labels.keys()))
        r2_labels = labeled_reactants[1].get_all_labeled_atoms()
        self.assertIn("*3", list(r2_labels.keys()))

    def test_get_reduced_mass(self):
        """
        tests that get_reduced_mass returns the proper value for H_abstraction
        """
        labels = ['*1', '*3']
        reactants = [
            Molecule().from_adjacency_list("""
1 *1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 *2 H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
"""),
            Molecule().from_adjacency_list("""
multiplicity 2
1 *3 H u1 p0 c0
""")]
        reduced_mass = get_reduced_mass(reactants, labels, True)
        self.assertAlmostEqual(reduced_mass, 1 / (1 / 1.008 + 1 / (1.008 + 12.01)) / 1000, places=6)

    def test_get_reduced_mass2(self):
        """
        tests that get_reduced_mass returns proper value when isotopes are labeled
        """
        labels = ['*1', '*3']
        reactants = [
            Molecule().from_adjacency_list("""
1 *1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
2 *2 H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
"""),
            Molecule().from_adjacency_list("""
multiplicity 2
1 *3 H u1 p0 c0
""")]
        reduced_mass = get_reduced_mass(reactants, labels, True)
        self.assertAlmostEqual(reduced_mass, 1 / (1 / 1.008 + 1 / (1.008 + 13.01)) / 1000, places=6)

    def test_get_kinetic_isotope_effect_simple(self):
        reactant_pair = [Species().from_smiles("C"), Species().from_smiles("[H]")]
        product_pair = [Species().from_smiles("[H][H]"), Species().from_smiles("[CH3]")]
        rxn_unlabeled = TemplateReaction(reactants=reactant_pair,
                                         products=product_pair,
                                         family='H_Abstraction',
                                         kinetics=Arrhenius(A=(1e5, 'cm^3/(mol*s)'), Ea=(0, 'J/mol')))
        rxn_labeled = TemplateReaction(reactants=[
            Species().from_adjacency_list("""
1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
"""),
            Species().from_adjacency_list("""
multiplicity 2
1 H u1 p0 c0
""")],
                                       products=[
            Species().from_adjacency_list("""
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
            Species().from_adjacency_list("""
multiplicity 2
1 C u1 p0 c0 i13 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""")],
                                       family='H_Abstraction',
                                       kinetics=Arrhenius(A=(1e5, 'cm^3/(mol*s)'), Ea=(0, 'J/mol')))
        rxn_cluster = [[rxn_labeled, rxn_unlabeled]]
        apply_kinetic_isotope_effect_simple(rxn_cluster, self.database.kinetics)
        expected_kie = ((1 / 1.008 + 1 / (13.01 + 1.008)) / (1 / 1.008 + 1 / (12.01 + 1.008))) ** 0.5
        self.assertAlmostEqual(rxn_cluster[0][0].kinetics.A.value, 1e5 * expected_kie, places=-1)

    def test_generate_isotope_reactions(self):
        """
        shows that all isotope reactions are created with generateIsotopeReactions
        """
        methyl = Species().from_smiles('[CH3]')
        methyli = Species().from_adjacency_list("""
multiplicity 2
1 C u1 p0 c0 i13 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""")
        methane = Species().from_adjacency_list("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
""")
        methanei = Species().from_adjacency_list("""
1 C u0 p0 c0 i13 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
""")
        craigie = Species().from_adjacency_list("""
multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""")
        craigiei = Species().from_adjacency_list("""
multiplicity 3
1 C u2 p0 c0 i13 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""")
        reaction = TemplateReaction(reactants=[methyl, methyl],
                                    products=[methane, craigie],
                                    family='H_Abstraction',
                                    template=['Xrad_H', 'Y_rad'],
                                    degeneracy=3)
        reaction.kinetics = Arrhenius(A=(1e5, 'cm^3/(mol*s)'), Ea=(0, 'J/mol'))
        isotope_list = [[methyl, methyli],
                        [methane, methanei],
                        [craigie, craigiei]]

        new_reactions = generate_isotope_reactions([reaction], isotope_list)

        self.assertEqual(len(new_reactions), 4)

        degeneracies_found = set()
        for rxn in new_reactions:
            self.assertEqual(rxn.template, reaction.template)
            degeneracies_found.add(rxn.degeneracy)
            self.assertIsNotNone(rxn.kinetics, 'kinetics not obtained for reaction {}.'.format(rxn))
            self.assertAlmostEqual(reaction.kinetics.get_rate_coefficient(298),
                                   rxn.kinetics.get_rate_coefficient(298) * reaction.degeneracy / rxn.degeneracy)

        self.assertEqual(degeneracies_found, set([3]))

    def test_generate_isotope_reactions_limited_labeling(self):
        """
        shows that all isotope reactions are created with generateIsotopeReactions
        with limits of two isotopes per molecule
        """
        max_number_labels = 1
        methyl = Species().from_smiles('[CH3]')
        methyl_isotopologues = [methyl] + generate_isotopomers(methyl, max_number_labels)
        methane = Species().from_smiles('C')
        methane_isotopologues = [methane] + generate_isotopomers(methane, max_number_labels)
        ethyl = Species().from_smiles('C[CH2]')
        ethyl_isotopologues = [ethyl] + generate_isotopomers(ethyl, max_number_labels)
        ethane = Species().from_smiles('CC')
        ethane_isotopologues = [ethane] + generate_isotopomers(ethane, max_number_labels)

        self.assertEqual(len(methyl_isotopologues), 2)
        self.assertEqual(len(methane_isotopologues), 2)
        self.assertEqual(len(ethane_isotopologues), 2)
        self.assertEqual(len(ethyl_isotopologues), 3)

        reaction = TemplateReaction(reactants=[ethyl, methane],
                                    products=[ethane, methyl],
                                    family='H_Abstraction',
                                    template=['C/H4', 'Y_rad'],
                                    degeneracy=4)
        reaction.kinetics = Arrhenius(A=(1e5, 'cm^3/(mol*s)'), Ea=(0, 'J/mol'))
        isotope_list = [methyl_isotopologues,
                        methane_isotopologues,
                        ethyl_isotopologues,
                        ethane_isotopologues]

        new_reactions = generate_isotope_reactions([reaction], isotope_list)

        self.assertEqual(len(new_reactions), 6)

        degeneracies_found = set()
        for rxn in new_reactions:
            self.assertEqual(rxn.template, reaction.template)
            degeneracies_found.add(rxn.degeneracy)
            self.assertIsNotNone(rxn.kinetics, 'kinetics not obtained for reaction {}.'.format(rxn))
            self.assertAlmostEqual(reaction.kinetics.get_rate_coefficient(298),
                                   rxn.kinetics.get_rate_coefficient(298) * reaction.degeneracy / rxn.degeneracy)

        self.assertEqual(degeneracies_found, set([4]))
