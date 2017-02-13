import unittest
import os
import shutil
import numpy as np

from rmgpy.tools.loader import loadRMGJob
import rmgpy
from rmgpy.species import Species
from rmgpy.tools.isotopes import *

from rmgpy.reaction import Reaction
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.kinetics.arrhenius import Arrhenius

class IsotopesTest(unittest.TestCase):

    def testClusterWithSpecies(self):
        """
        Test that isotope partitioning algorithm work with Reaction Objects.
        """

        eth = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        ethi = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        meth = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
        """
    )

        spcList = [eth, ethi]

        clusters = cluster(spcList)
        self.assertEquals(len(clusters), 1)
        self.assertEquals(len(clusters[0]), 2)

        spcList = [meth, ethi]

        clusters = cluster(spcList)
        self.assertEquals(len(clusters), 2)
        self.assertEquals(len(clusters[0]), 1)

    def testClusterWithReactions(self):
        """
        Test that isotope partitioning algorithm works with Reaction objects
        """

        eth = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        ethi = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        meth = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
        """
    )

        rxn0 = Reaction(reactants=[ethi,ethi], products=[ethi,eth])
        rxn1 = Reaction(reactants=[eth,ethi], products=[eth,eth])
        rxn2 = Reaction(reactants=[ethi,meth], products=[meth,ethi])
        rxn3 = Reaction(reactants=[eth,meth], products=[eth,meth])
        rxn4 = Reaction(reactants=[meth], products=[meth])
        rxn5 = Reaction(reactants=[ethi], products=[eth])
        
        
        sameClusterList = [rxn0,rxn1]

        clusters = cluster(sameClusterList)
        self.assertEquals(len(clusters), 1)
        self.assertEquals(len(clusters[0]), 2)
        
        sameClusterList = [rxn2,rxn3]

        clusters = cluster(sameClusterList)
        self.assertEquals(len(clusters), 1)
        self.assertEquals(len(clusters[0]), 2)

        multiClusterList = [rxn0, rxn1, rxn2, rxn3, rxn4, rxn5]

        clusters = cluster(multiClusterList)
        self.assertEquals(len(clusters), 4)
        self.assertEquals(len(clusters[0]), 1)

        
    def testRemoveIsotopeForReactions(self):
        """
        Test that remove isotope algorithm works with Reaction objects.
        """

        eth = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        ethi = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )
        unlabeledRxn = Reaction(reactants=[eth], products = [eth])
        labeledRxn = Reaction(reactants=[ethi], products = [ethi])
        stripped = removeIsotope(labeledRxn)
        
        self.assertTrue(unlabeledRxn.isIsomorphic(stripped))

    def testRemoveIsotopeForSpecies(self):
        """
        Test that remove isotope algorithm works with Species.
        """

        eth = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        ethi = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )

        stripped = removeIsotope(ethi)
        
        self.assertTrue(eth.isIsomorphic(stripped))

    def testCorrectAFactorsForMethylRecombination(self):
        """
        Test that correct A factors are used for isotopomers with various symmetries - methyl recombination
        """   
        correctAdjustment = 1.
        methyl = Species().fromAdjacencyList(
        """
multiplicity 2
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
        """)
        
        methyli = Species().fromAdjacencyList(
        """
multiplicity 2
1 C u1 p0 c0 i13 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
        """)
        
        ethanei = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 i13 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
        """
    )
        reaction = TemplateReaction(reactants = [methyl, methyli],
                            products = [ethanei],
                            kinetics = Arrhenius(A=(1.,'cm^3/(mol*s)'), Ea = (2., 'kJ/mol'), n=0.),
                            family = 'R_Recombination')
        reactionRateChange = self.__AFactorRatio(reaction)
        self.assertAlmostEqual(reactionRateChange, correctAdjustment)
        
    def testCorrectAFactorsForHadditionAllyl(self):
        """
        Test that correct A factors are used for isotopomers with various symmetries - H addition allyl
        """   
        correctAdjustment = .5
        h = Species().fromAdjacencyList(
        """
multiplicity 2
1 H u1 p0 c0
        """)

        allyli = Species().fromAdjacencyList(
        """
multiplicity 2
1 C u0 p0 c0 i13 {2,D} {6,S} {7,S}
2 C u0 p0 c0 {1,D} {3,S} {8,S}
3 C u1 p0 c0 {2,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}
8 H u0 p0 c0 {2,S}

        """)
        
        propenei = Species().fromAdjacencyList(
        """
1 C u0 p0 c0 {2,D} {4,S} {5,S}
2 C u0 p0 c0 {1,D} {3,S} {6,S}
3 C u0 p0 c0 i13 {2,S} {7,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}
        """
    )
        reaction = TemplateReaction(reactants = [allyli,h],
                            products = [propenei],
                            kinetics = Arrhenius(A=(1.,'cm^3/(mol*s)'), Ea = (2., 'kJ/mol'), n=0.),
                            family = 'R_Recombination')
        reactionRateChange = self.__AFactorRatio(reaction)
        self.assertAlmostEqual(reactionRateChange, correctAdjustment)
        
    def testCorrectAFactorsForHabsPropane(self):
        """
        Test that correct A factors are used for isotopomers with various symmetries - H abstraction from center propane
        """   
        correctAdjustment = 1
        h = Species().fromAdjacencyList(
        """
multiplicity 2
1 H u1 p0 c0
        """)

        h2 = Species().fromAdjacencyList(
        """
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
        """)
        propanei = Species().fromAdjacencyList(
        """
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
        
        ipropyli = Species().fromAdjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 i13 {2,S} {5,S} {6,S} {7,S}
2  C u1 p0 c0 {1,S} {3,S} {4,S}
3  H u0 p0 c0 {2,S}
4  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
        """
    )
        reaction = TemplateReaction(reactants = [propanei,h],
                            products = [ipropyli,h2],
                            kinetics = Arrhenius(A=(1.,'cm^3/(mol*s)'), Ea = (2., 'kJ/mol'), n=0.),
                            family = 'H_Abstraction')
        reactionRateChange = self.__AFactorRatio(reaction)
        self.assertAlmostEqual(reactionRateChange, correctAdjustment)
        
    def testCorrectAFactorsForHabsPropane2(self):
        """
        Test that correct A factors are used for isotopomers with various symmetries - H abstraction from edge propane
        """   
        correctAdjustment = 0.5
        h = Species().fromAdjacencyList(
        """
multiplicity 2
1 H u1 p0 c0
        """)

        h2 = Species().fromAdjacencyList(
        """
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
        """)
        propanei = Species().fromAdjacencyList(
        """
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
        
        npropyli = Species().fromAdjacencyList(
        """
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

        """
    )
        reaction = TemplateReaction(reactants = [propanei,h],
                            products = [npropyli,h2],
                            kinetics = Arrhenius(A=(1.,'cm^3/(mol*s)'), Ea = (2., 'kJ/mol'), n=0.),
                            family = 'H_Abstraction')
        reactionRateChange = self.__AFactorRatio(reaction)
        self.assertAlmostEqual(reactionRateChange, correctAdjustment)
    
    def testCorrectAFactorsForRAddMultBond(self):
        """
        Test that correct A factors are used for isotopomers with identical reactants - R addition multiple bond
        
        The rate of reactions with identical reactants should already 
        be accounted for in the degeneracy term, so this method just ensures
        the correct rate is used.
        """   
        correctAdjustment = 1
        butenyl = Species().fromAdjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,D} {7,S} {8,S}
2  C u0 p0 c0 {1,D} {3,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
4  C u1 p0 c0 {3,S} {5,S} {6,S}
5  H u0 p0 c0 {4,S}
6  H u0 p0 c0 {4,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}

        """)

        butenyli = Species().fromAdjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,D} {7,S} {8,S}
2  C u0 p0 c0 {1,D} {3,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
4  C u1 p0 c0 i13 {3,S} {5,S} {6,S}
5  H u0 p0 c0 {4,S}
6  H u0 p0 c0 {4,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}

        """)
        product = Species().fromAdjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,D} {13,S} {14,S}
2  C u0 p0 c0 {1,D} {3,S} {15,S}
3  C u0 p0 c0 {2,S} {4,S} {16,S} {17,S}
4  C u0 p0 c0 {3,S} {5,S} {18,S} {19,S}
5  C u0 p0 c0 {4,S} {6,S} {9,S} {20,S}
6  C u1 p0 c0 {5,S} {7,S} {8,S}
7  H u0 p0 c0 {6,S}
8  H u0 p0 c0 {6,S}
9  C u0 p0 c0 {5,S} {10,S} {21,S} {22,S}
10 C u1 p0 c0 i13 {9,S} {11,S} {12,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
        """)

        reaction = TemplateReaction(reactants = [butenyl,butenyli],
                            products = [product],
                            kinetics = Arrhenius(A=(1.,'cm^3/(mol*s)'), Ea = (2., 'kJ/mol'), n=0.),
                            family = 'R_Addition_MultipleBond')
        reactionRateChange = self.__AFactorRatio(reaction)
        self.assertAlmostEqual(reactionRateChange, correctAdjustment)
    
    def __AFactorRatio(self, reaction):
        """
        This helper function uses `correctAFactorsFromIsotopomers` to modify the A factor on
        the reaction. It then returns the ratio of 'correct A factor/original A factor'
        """
        
        originalAFactor = reaction.kinetics.A.value
        correctAFactorsOfIsotopomers([reaction])
        newAFactor = reaction.kinetics.A.value
        return newAFactor / originalAFactor
        
    def testEntireAFactorProcessFunctions(self):
        """
        Test that correctAFactors method effectively can cluster and modify A factor rates.
        """
        from rmgpy.quantity import ScalarQuantity
        # reaction of one type
        butenyl = Species().fromAdjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,D} {7,S} {8,S}
2  C u0 p0 c0 {1,D} {3,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
4  C u1 p0 c0 {3,S} {5,S} {6,S}
5  H u0 p0 c0 {4,S}
6  H u0 p0 c0 {4,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}

        """)
        butenyli = Species().fromAdjacencyList(
        """
multiplicity 2
1  C u0 p0 c0 {2,D} {7,S} {8,S}
2  C u0 p0 c0 {1,D} {3,S} {9,S}
3  C u0 p0 c0 {2,S} {4,S} {10,S} {11,S}
4  C u1 p0 c0 i13 {3,S} {5,S} {6,S}
5  H u0 p0 c0 {4,S}
6  H u0 p0 c0 {4,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {3,S}

        """)
        producti = Species().fromAdjacencyList(
        """
multiplicity 3
1  C u0 p0 c0 {2,D} {13,S} {14,S}
2  C u0 p0 c0 {1,D} {3,S} {15,S}
3  C u0 p0 c0 {2,S} {4,S} {16,S} {17,S}
4  C u0 p0 c0 {3,S} {5,S} {18,S} {19,S}
5  C u0 p0 c0 {4,S} {6,S} {9,S} {20,S}
6  C u1 p0 c0 {5,S} {7,S} {8,S}
7  H u0 p0 c0 {6,S}
8  H u0 p0 c0 {6,S}
9  C u0 p0 c0 {5,S} {10,S} {21,S} {22,S}
10 C u1 p0 c0 i13 {9,S} {11,S} {12,S}
11 H u0 p0 c0 {10,S}
12 H u0 p0 c0 {10,S}
13 H u0 p0 c0 {1,S}
14 H u0 p0 c0 {1,S}
15 H u0 p0 c0 {2,S}
16 H u0 p0 c0 {3,S}
17 H u0 p0 c0 {3,S}
18 H u0 p0 c0 {4,S}
19 H u0 p0 c0 {4,S}
20 H u0 p0 c0 {5,S}
21 H u0 p0 c0 {9,S}
22 H u0 p0 c0 {9,S}
        """)

        reaction1 = TemplateReaction(reactants = [butenyl,butenyli],
                            products = [producti],
                            kinetics = Arrhenius(A=(1.,'cm^3/(mol*s)'), Ea = (2., 'kJ/mol'), n=0.),
                            family = 'R_Addition_MultipleBond')

        # adding two more reactions

        h = Species().fromAdjacencyList(
        """
multiplicity 2
1 H u1 p0 c0
        """)
        h2 = Species().fromAdjacencyList(
        """
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
        """)
        propanei = Species().fromAdjacencyList(
        """
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
        propane = Species().fromAdjacencyList(
        """
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
        npropyli = Species().fromAdjacencyList(
        """
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
        npropyl = Species().fromAdjacencyList(
        """
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

        reaction2 = TemplateReaction(reactants = [propanei,h],
                            products = [npropyli,h2],
                            kinetics = Arrhenius(A=(1.,'cm^3/(mol*s)'), Ea = (2., 'kJ/mol'), n=0.),
                            family = 'H_Abstraction')
        reaction3 = TemplateReaction(reactants = [propane,h],
                            products = [npropyl,h2],
                            kinetics = Arrhenius(A=(1.,'cm^3/(mol*s)'), Ea = (2., 'kJ/mol'), n=0.),
                            family = 'H_Abstraction')

        reactions = [reaction1,reaction2,reaction3]

        correctAFactors(reactions)

        self.assertAlmostEqual(reaction1.kinetics.A.value, ScalarQuantity(1,'cm^3/(mol*s)').value)
        self.assertAlmostEqual(reaction2.kinetics.A.value, ScalarQuantity(0.5,'cm^3/(mol*s)').value)
        self.assertAlmostEqual(reaction3.kinetics.A.value, ScalarQuantity(1,'cm^3/(mol*s)').value)

    def testComputeProbabilities(self):
        """
        Test that the retrieval of isotopomer concentrations works.
        """
        folder = os.path.join(os.path.dirname(rmgpy.__file__),'tools/data/isotopes')
        try:
            shutil.rmtree(os.path.join(folder,'solver'))
        except OSError, e:
            pass
        

        inputFile = os.path.join(folder, 'input_isotope.py')
        chemkinFile = os.path.join(folder, 'chemkin', 'chem_annotated.inp')
        dictFile = os.path.join(folder, 'chemkin', 'species_dictionary.txt')

        os.mkdir(os.path.join(folder, 'solver'))

        rmg = loadRMGJob(inputFile, chemkinFile, dictFile, generateImages=False, useChemkinNames=True)

        spcdata = solve(rmg)

        clusters = cluster(rmg.reactionModel.core.species)
        concs = retrieveConcentrations(spcdata, clusters)

        self.assertEquals(len(clusters), len(concs))
        for clust, df in zip(clusters, concs):
            self.assertEquals(len(clust), len(df.columns))

        # compute species probabilities
        probs = []
        for df in concs:
            df = computeProbabilities(df.copy())
            probs.append(df)

        for df in probs:
            sumConcs = df.sum(axis=1)
            sampleSize = 2
            for sample in sumConcs.sample(sampleSize):
                if not np.isnan(sample):
                    self.assertAlmostEqual(sample, 1.)            

        shutil.rmtree(os.path.join(folder,'solver'))

    def testGenerateIsotopomers(self):
        """
        Test that the generation of isotopomers with N isotopes works.
        """
        from rmgpy.thermo.nasa import NASAPolynomial, NASA

        spc = Species().fromSMILES('CC')

        polynomial = NASAPolynomial(coeffs=[1.,1.,1.,1.,1.,1.,1.],
                         Tmin=(200,'K'),Tmax=(1600,'K'),E0=(1.,'kJ/mol'),
                         comment='made up thermo')
        
        spc.thermo = NASA(polynomials=[polynomial],Tmin=(200,'K'),
                        Tmax=(1600,'K'),E0=(1.,'kJ/mol'),
                        comment='made up thermo')

        spcs = generateIsotopomers(spc, 0)
        self.assertEquals(len(spcs), 0)

        spcs = generateIsotopomers(spc)
        self.assertEquals(len(spcs), 1)

        spcs = generateIsotopomers(spc, 2)
        self.assertEquals(len(spcs), 2)

        spcs = generateIsotopomers(spc, 3)
        self.assertEquals(len(spcs), 2)
