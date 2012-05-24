'''
Created on May 19, 2012

@author: nmvdewie
'''
import unittest
import qmtp_package.qmtp as qm
import qmtp_package.rmg_qm_parsers as pars
import os
import rmgpy.molecule as mol
import qmtp_package.calculator as calc
from rmgpy.thermo import ThermoData
import qmtp_package.qmtpjobs as job
from qmtp_package.qmtpjobs import MOPACJob

class Test(unittest.TestCase):

    def testMOPACJob(self):
        name = 'UMRZSTCPUPJPOJ-UHFFFAOYAR'
        dir = os.path.join(os.getcwd(),'data/QMfiles/MOPAC')
        mj = MOPACJob(name, dir)
        success = mj.run()
        
        self.assertTrue(success)
        
        #remove generated output files
        os.remove(os.path.join(dir,name+'.out'))
        os.remove(os.path.join(dir,name+'.arc'))
        
        

    def testSymmetryJob(self):
        '''
        Check whether external call to symmetry tool works fine. 
        '''
        name = 'AAAOFKFEDKWQNN-UHFFFAOYAY'
        InChIaug = 'InChI=1/C9H14O2/c1-6(2)9-5-8(11-10)4-7(9)3/h4-6,8,10H,1-3H3'
        molecule = mol.Molecule().fromInChI(InChIaug)
        inputFileExtension = '.log'
        driver = qm.QMTP('gaussian03', 'pm3')
        dir = os.path.join(os.getcwd(),'data/QMfiles/G03')
        parsingTool = pars.CCLibParser(os.path.join(dir,name+inputFileExtension), driver)

        iqmdata = parsingTool.parse(molecule);
        path = os.environ.get('RMG_workingDirectory')
        symm_job = job.SymmetryJob(name, dir, iqmdata)
        pointGroup = symm_job.calculate()
        self.assertTrue(os.path.exists(os.path.join(dir, 'AAAOFKFEDKWQNN-UHFFFAOYAY.symm')))
        
        self.assertEqual(pointGroup.symmetryNumber, 1)
        self.assertTrue(pointGroup.chiral)#the chiral flag is set to True for C1 symmetry groups!
        self.assertFalse(pointGroup.linear)
        

if __name__ == "__main__":
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )