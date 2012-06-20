"""
Created on May 17, 2012

@author: nmvdewie
"""
import unittest
import rmgpy.qm.qmtp as qm
import rmgpy.qm.parsers as pars
import os
import rmgpy.molecule as mol
import rmgpy.qm.calculator as calc
from rmgpy.thermo import ThermoData
import re

class Test(unittest.TestCase):

    def testMOPAC_PM3_Parser(self):
        driver = qm.QMTP('mopac')
        name = 'GRWFGVWFFZKLTI-UHFFFAOYAF'
        InChIaug = 'InChI=1/C10H16/c1-7-4-5-8-6-9(7)10(8,2)3/h4,8-9H,5-6H2,1-3H3'
        molecule = mol.Molecule().fromInChI(InChIaug)
        directory = os.path.join(os.path.dirname(__file__),'data','QMfiles','MOPAC')
        mf = qm.molFile(molecule, name, directory)
        parser = pars.MOPACPM3Parser(mf, driver)
        result = parser.parse()
        assert isinstance(result, ThermoData)
        
    def testCCLibParser(self):
        
        #Tests to check whether the CCLibParser wrapper around CCLib works fine
        
        name = 'AAAOFKFEDKWQNN-UHFFFAOYAY'
        InChIaug = 'InChI=1/C9H14O2/c1-6(2)9-5-8(11-10)4-7(9)3/h4-6,8,10H,1-3H3'
        molecule = mol.Molecule().fromInChI(InChIaug)
        inputFileExtension = '.log'
        driver = qm.QMTP('gaussian03', 'pm3')
        directory = os.path.join(os.path.dirname(__file__),'data','QMfiles','G03')
        parsingTool = pars.CCLibParser(os.path.join(directory, name+inputFileExtension), driver)
        
        data = parsingTool.parse(molecule)
        
        self.assertEqual(data.groundStateDegeneracy, 1)
        self.assertEqual(data.cclib_data.natom, 25)
        self.assertEqual(data.cclib_data.molmass, 154.09938)
        self.assertEqual(len(data.cclib_data.atomcoords[-1]), data.cclib_data.natom)
        self.assertEqual(len(data.atomCoords), data.numberOfAtoms)
        self.assertEqual(len(data.cclib_data.rotcons[-1]), 3)
        self.assertEqual(len(data.cclib_data.atomnos), data.cclib_data.natom)
        
        print 'Z-coord Atom coord of 1st atom: '+str(data.cclib_data.atomcoords[-1][0][2])
        
     
    
    
if __name__ == "__main__":
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )