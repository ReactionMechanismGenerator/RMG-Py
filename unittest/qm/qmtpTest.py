"""
Created on May 17, 2012

@author: nmvdewie
"""
import unittest

from rmgpy.molecule import *
import rmgpy.qm.qmtp as qm
from rmgpy.thermo import ThermoData

class QMTest(unittest.TestCase):

    def testGenerateIdentifiers(self):
        molecule = Molecule().fromSMILES('C12CC(CC1)CC2')#norbornane
        qmprogram = ''
        driver = qm.QMTP(qmprogram)
        name, inchi = driver.generateIdentifiers(molecule)
        self.assertEqual(name,'UMRZSTCPUPJPOJ-UHFFFAOYSA')
        self.assertEqual(inchi , 'InChI=1S/C7H12/c1-2-7-4-3-6(1)5-7/h6-7H,1-5H2')

    def test3DMolFileCreator(self):

        molecule = Molecule().fromSMILES('C12CC(CC1)CC2')#norbornane
        qmprogram = ''
        driver = qm.QMTP(qmprogram)
        name, inchi = driver.generateIdentifiers(molecule)#UMRZSTCPUPJPOJ-UHFFFAOYSA
        
        directory = os.path.join(os.path.dirname(__file__), 'data','ThreeDMolFileCreatorTest')
        print os.path.exists(directory)
        creator = qm.ThreeDMolFileCreator(name, directory, molecule)
        
        threedmolfile = creator.create()
        
        self.assertTrue(isinstance(threedmolfile, qm.molFile))
        

    def testGenerateQMThermoData(self):
        """
        Checks whether the QMTP.generateQMThermoData() finds the output file of the molecule
        for which already MOPAC PM3 results exist.
        """
        driver = qm.QMTP('mopac', 'pm3')
        #corresponds to GRWFGVWFFZKLTI-UHFFFAOYAF.out
        molecule = Molecule().fromInChI('InChI=1/C10H16/c1-7-4-5-8-6-9(7)10(8,2)3/h4,8-9H,5-6H2,1-3H3')
        qm.QMTP.qmfolder = os.path.join(os.path.dirname(__file__), 'data','QMfiles')
        result = driver.generateQMThermoData(molecule)
        
        self.assertTrue(isinstance(result, ThermoData))
        from glob import glob
        for f in glob (os.path.join(qm.QMTP.qmfolder,'GRWFGVWFFZKLTI-UHFFFAOYSA.*')):
            os.unlink (f)
        for f in glob (os.path.join(qm.QMTP.qmfolder,'3DMolfiles','GRWFGVWFFZKLTI-UHFFFAOYSA.*')):
            os.unlink (f)
        
if __name__ == "__main__":
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )