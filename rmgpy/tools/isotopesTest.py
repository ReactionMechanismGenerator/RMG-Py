import unittest
import os
import shutil
import numpy as np

from rmgpy.species import Species
from rmgpy.tools.loader import loadRMGJob

from rmgpy.species import Species
from rmgpy.tools.isotopes import *

class IsotopesTest(unittest.TestCase):

    def testCluster(self):
        """
        Test that isotope partitioning algorithm works.
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

    def testRemoveIsotope(self):
        """
        Test that remove isotope algorithm works.
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

        spc = Species().fromSMILES('CC')

        spcs = generateIsotopomers(spc, 0)
        self.assertEquals(len(spcs), 0)

        spcs = generateIsotopomers(spc)
        self.assertEquals(len(spcs), 1)

        spcs = generateIsotopomers(spc, 2)
        self.assertEquals(len(spcs), 2)

        spcs = generateIsotopomers(spc, 3)
        self.assertEquals(len(spcs), 2)