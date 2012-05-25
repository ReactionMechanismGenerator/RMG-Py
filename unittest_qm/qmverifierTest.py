'''
Created on May 17, 2012

@author: nmvdewie
'''
import unittest
import qmtp_package.qmtp as qm
import os
import qmtp_package.qmverifier as verif

class Test(unittest.TestCase):

    def testVerifierDoesNotExist(self):
        driver = qm.QMTP('both', 'pm3')
        driver.qmMethod = 'pm3'
        
        name = 'UMRZSTCPUPJPOJ-UHFFFAOYSA'
        InChIaug = 'InChI=1S/C7H12/c1-2-7-4-3-6(1)5-7/h6-7H,1-5H2'
        directory = os.path.join(os.path.dirname(__file__),'data','QMfiles')
        verifier = verif.QMVerifier(name, InChIaug, directory, QMTP = driver)
        verifier.verify()
        
        self.assertFalse(verifier.succesfulJobExists())
     
    def testVerifierMOPACResultExists(self):
        driver = qm.QMTP('mopac', 'pm3')

        
        name = 'GRWFGVWFFZKLTI-UHFFFAOYAF'
        InChIaug = 'InChI=1/C10H16/c1-7-4-5-8-6-9(7)10(8,2)3/h4,8-9H,5-6H2,1-3H3'
        directory = os.path.join(os.path.dirname(__file__),'data','QMfiles','MOPAC')
        verifier = verif.QMVerifier(name, InChIaug, directory, driver)
        verifier.verify()
        self.assertTrue(verifier.succesfulJobExists())
        self.assertTrue(verifier.mopacResultExists)
        self.assertFalse(verifier.gaussianResultExists)
    
    
    
if __name__ == "__main__":
        unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )