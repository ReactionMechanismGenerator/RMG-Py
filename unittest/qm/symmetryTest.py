"""
Created on May 18, 2012

@author: nmvdewie
"""
import unittest
import os
import rmgpy.molecule as mol
import rmgpy.qm.calculator as calc
from rmgpy.thermo import ThermoData
import re
import rmgpy.qm.qmtp as qm
import rmgpy.qm.parsers as pars
import rmgpy.qm.symmetry as symm
class Test(unittest.TestCase):

    def testPGDictionary(self):
        """
        Check the dictionary class that collects information on symmetry groups, and their
        chirality flags
        """
        pgd = symm.PointGroupDictionary()
        
        self.assertFalse(len(pgd.library.keys()) == 0)

    def testPointGroupCalculator(self):
        """
        Tests PointGroupCalculator
        """
        name = 'AAAOFKFEDKWQNN-UHFFFAOYAY'
        InChIaug = 'InChI=1/C9H14O2/c1-6(2)9-5-8(11-10)4-7(9)3/h4-6,8,10H,1-3H3'
        molecule = mol.Molecule().fromInChI(InChIaug)
        inputFileExtension = '.log'
        directory = os.path.join(os.path.dirname(__file__),'data/QMfiles/G03')
        driver = qm.QMTP('gaussian03', 'pm3')
        parsingTool = pars.CCLibParser(os.path.join(directory, name+inputFileExtension), driver)
        
        data = parsingTool.parse(molecule)
        mf = qm.molFile(molecule, name, directory)
        pgc = symm.PointGroupCalculator(mf, data)
        pg = pgc.calculate()
        
        self.assertTrue(isinstance(pg, symm.PointGroup))
        self.assertFalse(pg.linear)
        self.assertTrue(pg.chiral)#the chiral flag is set to True for C1 symmetry groups!
        self.assertEqual(pg.symmetryNumber, 1)
        
if __name__ == "__main__":
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )