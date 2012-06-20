"""
Created on May 23, 2012

@author: nmvdewie
"""
import time
import unittest
import rmgpy.qm.inputwriters as writers
import os
import rmgpy.qm.qmtp as qm
import rmgpy.molecule as mol

class Test(unittest.TestCase):

    def testMOPACInputWriter(self):
        """
        Checks whether the .mop output file has been written based on the 3D coords file (.mol)
        """
        
        name = 'WTARULDDTDQWMU-UHFFFAOYAW'
        inchi = 'InChI=1/C10H16/c1-7-4-5-8-6-9(7)10(8,2)3/h8-9H,1,4-6H2,2-3H3'
        directory = os.path.join(os.path.dirname(__file__), 'data','QMfiles','3DMolfiles')
        target_file = os.path.join(directory, name+'.mop')
        if os.path.exists(target_file): os.remove(target_file)
        molecule = mol.Molecule().fromInChI(inchi)
        mf = qm.molFile(molecule, name, directory)
        
        writer = writers.MOPACPM3InputWriter(mf, attemptNumber=1, multiplicity=1)
        writer.write()
        
        time.sleep(3)#otherwise assertion fails before the file is written!
        self.assertTrue(os.path.exists(target_file))
    
    def testG03InputWriter(self):
        """
        Checks whether the .gjf output file has been written based on the 3D coords file (.mol)
        """
        name = 'WTARULDDTDQWMU-UHFFFAOYAW'
        inchi = 'InChI=1/C10H16/c1-7-4-5-8-6-9(7)10(8,2)3/h8-9H,1,4-6H2,2-3H3'
        directory = os.path.join(os.path.dirname(__file__), 'data','QMfiles','3DMolfiles')
        target_file = os.path.join(directory, name+'.gjf')
        if os.path.exists(target_file): os.remove(target_file)
        molecule = mol.Molecule().fromInChI(inchi)
        mf = qm.molFile(molecule, name, directory)
        
        writer = writers.GaussianPM3InputWriter(mf, attemptNumber=1, multiplicity=1)
        writer.write()
        
        time.sleep(3)#otherwise assertion fails before the file is written!
        self.assertTrue(os.path.exists(target_file))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()