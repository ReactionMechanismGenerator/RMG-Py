"""
Created on May 19, 2012

@author: nmvdewie
"""
import unittest
import rmgpy.qm.qmtp as qm
import rmgpy.qm.parsers as pars
import os
import rmgpy.molecule as mol
import rmgpy.qm.calculator as calc
from rmgpy.thermo import ThermoData
import rmgpy.qm.jobs as job
from rmgpy.qm.jobs import MOPACJob
import rmgpy.qm.qmtp as qm
class Test(unittest.TestCase):

    def testMOPACJob(self):
        molecule = mol.Molecule()
        name = 'UMRZSTCPUPJPOJ-UHFFFAOYAR'
        directory = os.path.join(os.path.dirname(__file__),'data','QMfiles','MOPAC')
        InChIAug = ''
        
        molfile = qm.molFile(molecule, name, directory, InChIAug)
        mj = MOPACJob(molfile)
        success = mj.run()
        
        self.assertTrue(success)
        
        #remove generated output files
        os.remove(os.path.join(directory, name+'.out'))
        os.remove(os.path.join(directory, name+'.arc'))
        
    def testSymmetryJob(self):
        """
        Check whether external call to symmetry tool works fine. 
        """
        name = 'AAAOFKFEDKWQNN-UHFFFAOYAY'
        directory = os.path.join(os.path.dirname(__file__),'data','QMfiles','G03')
        InChIaug = 'InChI=1/C9H14O2/c1-6(2)9-5-8(11-10)4-7(9)3/h4-6,8,10H,1-3H3'
        molecule = mol.Molecule().fromInChI(InChIaug)
        molfile = qm.molFile(molecule, name, directory, InChIaug)
        
        inputFileExtension = '.log'
        driver = qm.QMTP('gaussian03', 'pm3')

        parsingTool = pars.CCLibParser(os.path.join(directory,name+inputFileExtension), driver)

        iqmdata = parsingTool.parse(molecule);
        symm_job = job.SymmetryJob(molfile, iqmdata)
        pointGroup = symm_job.calculate()
        self.assertTrue(os.path.exists(os.path.join(directory, 'AAAOFKFEDKWQNN-UHFFFAOYAY.symm')))
        
        self.assertEqual(pointGroup.symmetryNumber, 1)
        self.assertTrue(pointGroup.chiral)#the chiral flag is set to True for C1 symmetry groups!
        self.assertFalse(pointGroup.linear)
        

if __name__ == "__main__":
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )