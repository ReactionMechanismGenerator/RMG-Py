"""
Created on May 17, 2012

@author: nmvdewie
"""
import unittest
import rmgpy.qm.qmtp as qm
import os
import rmgpy.qm.qmverifier as verif
import rmgpy.molecule as mol
class Test(unittest.TestCase):

    def testVerifierDoesNotExist(self):
        molecule = mol.Molecule()
        name = 'UMRZSTCPUPJPOJ-UHFFFAOYSA'
        directory = os.path.join(os.path.dirname(__file__),'data','QMfiles')
        InChIaug = 'InChI=1S/C7H12/c1-2-7-4-3-6(1)5-7/h6-7H,1-5H2'
        molfile = qm.molFile(molecule, name, directory, InChIaug)
        
        verifier = verif.QMVerifier(molfile)
        verifier.verify()
        
        self.assertFalse(verifier.succesfulJobExists())
     
    def testVerifierMOPACResultExists(self):
        molecule = mol.Molecule()
        name = 'GRWFGVWFFZKLTI-UHFFFAOYAF'
        directory = os.path.join(os.path.dirname(__file__),'data','QMfiles','MOPAC')
        InChIaug = 'InChI=1/C10H16/c1-7-4-5-8-6-9(7)10(8,2)3/h4,8-9H,5-6H2,1-3H3'
        molfile = qm.molFile(molecule, name, directory, InChIaug)
        
        verifier = verif.QMVerifier(molfile)
        verifier.verify()
        self.assertTrue(verifier.succesfulJobExists())
        self.assertTrue(verifier.mopacResultExists)
        self.assertFalse(verifier.gaussianResultExists)
    
    
    
if __name__ == "__main__":
        unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )