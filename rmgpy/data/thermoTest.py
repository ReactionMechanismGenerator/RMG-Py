#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import unittest
from external.wip import work_in_progress

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase, database
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import Species
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.molecule.molecule import Molecule
import rmgpy

################################################################################

class TestThermoDatabase(unittest.TestCase):
    """
    Contains unit tests of the ThermoDatabase class.
    """
    # Only load these once to save time
    database = ThermoDatabase()
    database.load(os.path.join(settings['database.directory'], 'thermo'))
#    oldDatabase = ThermoDatabase()
#    oldDatabase.loadOld(os.path.join(settings['database.directory'], '../output/RMG_database'))
    
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        
        self.database = self.__class__.database
#        self.oldDatabase = self.__class__.oldDatabase

        self.Tlist = [300, 400, 500, 600, 800, 1000, 1500]
        
        self.testCases = [
            # SMILES            symm  H298     S298     Cp300  Cp400  Cp500  Cp600  Cp800  Cp1000 Cp1500
            
            # 1,3-hexadiene decomposition products
            ['C=CC=CCC',        3,    13.45, 86.37, 29.49, 37.67, 44.54, 50.12, 58.66, 64.95, 74.71],
            ['[CH]=CC=CCC',     3,    72.55, 87.76, 29.30, 36.92, 43.18, 48.20, 55.84, 61.46, 70.18],
            ['C=[C]C=CCC',      3,    61.15, 87.08, 29.68, 36.91, 43.03, 48.11, 55.96, 61.78, 71.54],
            ['C=C[C]=CCC',      3,    61.15, 87.08, 29.68, 36.91, 43.03, 48.11, 55.96, 61.78, 71.54],
            ['C=CC=[C]CC',      3,    70.35, 88.18, 29.15, 36.46, 42.6, 47.6, 55.32, 61.04, 69.95],
            ['C=CC=C[CH]C',     6,    38.24, 84.41, 27.79, 35.46, 41.94, 47.43, 55.74, 61.92, 71.86],
            ['C=CC=CC[CH2]',    2,    62.45, 89.78, 28.72, 36.31, 42.63, 47.72, 55.50, 61.21, 70.05],
            ['[CH3]',           6,    34.81, 46.37, 9.14, 10.18, 10.81, 11.34, 12.57, 13.71, 15.2],
            ['C=CC=C[CH2]',     2,    46.11, 75.82, 22.54, 28.95, 34.24, 38.64, 45.14, 49.97, 57.85],
            ['[CH2]C',          6,    28.6, 59.87, 11.73, 14.47, 17.05, 19.34, 23.02, 25.91, 31.53],
            ['C=CC=[CH]',       1,    85.18, 69.37, 18.93, 23.55, 27.16, 29.92, 34.02, 37.03, 41.81],
            ['C=[CH]',          1,    71.62, 56.61, 10.01, 11.97, 13.66, 15.08, 17.32, 19.05, 21.85],
            ['[CH]=CCC',        3,    58.99, 75.0, 20.38, 25.34, 29.68, 33.36, 39.14, 43.48, 50.22],
            
            # Cyclic Structures
            ['C1CCCCC1',        12,   -29.45, 69.71, 27.20, 37.60, 46.60, 54.80, 67.50, 76.20, 88.50],
            ['C1CCC1',          8,     6.51, 63.35, 17.39, 23.91, 29.86, 34.76, 42.40, 47.98, 56.33],
            ['C1C=CC=C1',       2,    32.5, 65.5, 18.16, 24.71, 30.25, 34.7, 41.25, 45.83, 52.61],
        ]

    @work_in_progress
    def testNewThermoGeneration(self):
        """
        Test that the new ThermoDatabase generates appropriate thermo data.
        """
        
        for smiles, symm, H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500 in self.testCases:
            Cplist = [Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500]
            molecule=Molecule(SMILES=smiles)
            species = Species(molecule=molecule)
            species.generateResonanceIsomers()
            species.molecule[0]
            thermoData = self.database.getThermoDataFromGroups(species)
            molecule = species.molecule[0]
            for mol in species.molecule[1:]:
                thermoData0 = self.database.getAllThermoData(Species(molecule=[mol]))[0][0]
                for data in self.database.getAllThermoData(Species(molecule=[mol]))[1:]:
                    if data.getEnthalpy(298) < thermoData0.getEnthalpy(298):
                        thermoData0 = data
                if thermoData0.getEnthalpy(298) < thermoData.getEnthalpy(298):
                    thermoData = thermoData0
                    molecule = mol
            self.assertAlmostEqual(H298, thermoData.getEnthalpy(298) / 4184, places=1, msg="H298 error for {0}".format(smiles))
            self.assertAlmostEqual(S298, thermoData.getEntropy(298) / 4.184, places=1, msg="S298 error for {0}".format(smiles))
            for T, Cp in zip(self.Tlist, Cplist):
                self.assertAlmostEqual(Cp, thermoData.getHeatCapacity(T) / 4.184, places=1, msg="Cp{1} error for {0}".format(smiles,T))

    def testSymmetryContributionRadicals(self):
        """
        Test that the symmetry contribution is correctly added for radicals
        estimated via the HBI method. 
        """
        spc = Species(molecule=[Molecule().fromSMILES('[CH3]')])
        
        thermoData_lib = self.database.getThermoDataFromLibraries(spc)[0]
        
        thermoData_ga = self.database.getThermoDataFromGroups(spc)
        
        self.assertAlmostEqual(thermoData_lib.getEntropy(298.), thermoData_ga.getEntropy(298.), 0)

        
    @work_in_progress
    def testSymmetryNumberGeneration(self):
        """
        Test we generate symmetry numbers correctly.
        
        This uses the new thermo database to generate the H298, used 
        to select the stablest resonance isomer.
        """
        for smiles, symm, H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500 in self.testCases:
            molecule=Molecule(SMILES=smiles)
            species = Species(molecule=molecule)
            species.generateResonanceIsomers()
            thermoData = self.database.getThermoDataFromGroups(Species(molecule=[species.molecule[0]]))
            # pick the molecule with lowest H298
            molecule = species.molecule[0]
            for mol in species.molecule[1:]:
                thermoData0 = self.database.getAllThermoData(Species(molecule=[mol]))[0][0]
                for data in self.database.getAllThermoData(Species(molecule=[mol]))[1:]:
                    if data.getEnthalpy(298) < thermoData0.getEnthalpy(298):
                        thermoData0 = data
                if thermoData0.getEnthalpy(298) < thermoData.getEnthalpy(298):
                    thermoData = thermoData0
                    molecule = mol
            self.assertEqual(molecule.calculateSymmetryNumber(), symm, msg="Symmetry number error for {0}".format(smiles))

#    @work_in_progress
#    def testOldThermoGeneration(self):
#        """
#        Test that the old ThermoDatabase generates relatively accurate thermo data.
#        """
#        for smiles, symm, H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500 in self.testCases:
#            Cplist = [Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500]
#            species = Species(molecule=[Molecule(SMILES=smiles)])
#            species.generateResonanceIsomers()
#            thermoData = self.oldDatabase.getThermoData(Species(molecule=[species.molecule[0]]))
#            molecule = species.molecule[0]
#            for mol in species.molecule[1:]:
#                thermoData0 = self.oldDatabase.getAllThermoData(Species(molecule=[mol]))[0][0]
#                for data in self.oldDatabase.getAllThermoData(Species(molecule=[mol]))[1:]:
#                    if data.getEnthalpy(298) < thermoData0.getEnthalpy(298):
#                        thermoData0 = data
#                if thermoData0.getEnthalpy(298) < thermoData.getEnthalpy(298):
#                    thermoData = thermoData0
#                    molecule = mol
#            
#            self.assertAlmostEqual(H298, thermoData.getEnthalpy(298) / 4184, places=1, msg="H298 error for {0}".format(smiles))
#            self.assertAlmostEqual(S298, thermoData.getEntropy(298) / 4.184, places=1, msg="S298 error for {0}".format(smiles))
#            for T, Cp in zip(self.Tlist, Cplist):
#                self.assertAlmostEqual(Cp, thermoData.getHeatCapacity(T) / 4.184, places=1, msg="Cp{1} error for {0}".format(smiles, T))


class TestThermoDatabaseAromatics(TestThermoDatabase):
    """
    Test only Aromatic species.
    
    A copy of the above class, but with different test compounds
    """
    def setUp(self):
        TestThermoDatabase.setUp(self)
        self.testCases = [
            # SMILES            symm  H298     S298     Cp300  Cp400  Cp500  Cp600  Cp800  Cp1000 Cp1500
            ['c1ccccc1', 12, 19.80, 64.24, 19.44, 26.64, 32.76, 37.80, 45.24, 50.46, 58.38],
            ['c1ccc2ccccc2c1', 4, 36.0, 79.49, 31.94, 42.88, 52.08, 59.62, 70.72, 78.68, 90.24],
        ]
    def __init__(self, *args, **kwargs):
        super(TestThermoDatabaseAromatics, self).__init__(*args, **kwargs)
        self._testMethodDoc = self._testMethodDoc.strip().split('\n')[0] + " for Aromatics.\n"


class TestCyclicThermo(unittest.TestCase):
    """
    Contains unit tests of the ThermoDatabase class.
    """
    database = ThermoDatabase()
    database.load(os.path.join(settings['database.directory'], 'thermo'))
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        
        self.database = self.__class__.database
    
    def testComputeGroupAdditivityThermoForTwoRingMolecule(self):
        """
        The molecule being tested has two rings, one is 13cyclohexadiene5methylene
        the other is benzene ring. This method is to test thermo estimation will
        give two different corrections accordingly. 
        """
        spec = Species().fromSMILES('CCCCCCCCCCCC(CC=C1C=CC=CC1)c1ccccc1')
        spec.generateResonanceIsomers()
        thermo = self.database.getThermoDataFromGroups(spec)

        ringGroups, polycyclicGroups = self.database.getRingGroupsFromComments(thermo)
        self.assertEqual(len(ringGroups),2)
        self.assertEqual(len(polycyclicGroups),0)

        expected_matchedRingsLabels = ['13cyclohexadiene5methylene', 'Benzene']
        expected_matchedRings = [self.database.groups['ring'].entries[label] for label in expected_matchedRingsLabels]

        self.assertEqual(set(ringGroups), set(expected_matchedRings))
    
    def testThermoForMonocyclicAndPolycyclicSameMolecule(self):
        """
        Test a molecule that has both a polycyclic and a monocyclic ring in the same molecule
        """
        spec = Species().fromSMILES('C(CCC1C2CCC1CC2)CC1CCC1')
        spec.generateResonanceIsomers()
        thermo = self.database.getThermoDataFromGroups(spec)
        ringGroups, polycyclicGroups = self.database.getRingGroupsFromComments(thermo)
        self.assertEqual(len(ringGroups),1)
        self.assertEqual(len(polycyclicGroups),1)
        
        expected_matchedRingsLabels = ['Cyclobutane']
        expected_matchedRings = [self.database.groups['ring'].entries[label] for label in expected_matchedRingsLabels]
        self.assertEqual(set(ringGroups), set(expected_matchedRings))
        
        expected_matchedPolyringsLabels = ['norbornane']
        expected_matchedPolyrings = [self.database.groups['polycyclic'].entries[label] for label in expected_matchedPolyringsLabels]

        self.assertEqual(set(polycyclicGroups), set(expected_matchedPolyrings))
    
    def testPolycyclicPicksBestThermo(self):
        """
        Test that RMG prioritizes thermo correctly and chooses the thermo from the isomer which
        has a non generic polycyclic ring correction
        """
        
        spec = Species().fromSMILES('C1=C[C]2CCC=C2C1')
        spec.generateResonanceIsomers()
        
        thermoDataList = []
        for molecule in spec.molecule:
            thermo = self.database.estimateRadicalThermoViaHBI(molecule, self.database.computeGroupAdditivityThermo)
            thermoDataList.append(thermo)
            
        thermoDataList.sort(key=lambda x: x.getEnthalpy(298))
        most_stable_thermo = thermoDataList[0]
        ringGroups, polycyclicGroups = self.database.getRingGroupsFromComments(most_stable_thermo)
        
        selected_thermo = self.database.getThermoDataFromGroups(spec)
        
        self.assertNotEqual(selected_thermo, thermoDataList)
        
        selected_ringGroups, selected_polycyclicGroups = self.database.getRingGroupsFromComments(selected_thermo)
        
        # The group used to estimate the most stable thermo is the generic polycyclic group and
        # therefore is not selected.  Note that this unit test will have to change if the correction is fixed later.
        self.assertEqual(polycyclicGroups[0].label, 'PolycyclicRing')
        self.assertEqual(selected_polycyclicGroups[0].label, 'C12CCC=C1CC=C2')
    

    def testGetRingGroupsFromComments(self):
        """
        Test that getRingGroupsFromComments method works for fused polycyclics.
        """
        from rmgpy.thermo.thermoengine import generateThermoData
        
        # set-up RMG object
        rmg = RMG()

        # load kinetic database and forbidden structures
        rmg.database = RMGDatabase()
        path = os.path.join(settings['database.directory'])

        # forbidden structure loading
        rmg.database.loadThermo(os.path.join(path, 'thermo'))

        smi = 'C12C(C3CCC2C3)C4CCC1C4'#two norbornane rings fused together
        spc = Species().fromSMILES(smi)

        spc.thermo = generateThermoData(spc)

        thermodb = rmg.database.thermo
        thermodb.getRingGroupsFromComments(spc.thermo)

        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None

    def testRemoveGroup(self):
        """
        Test that removing groups using nodes near the root of radical.py
        """
        #load up test data designed for this test
        database2 = ThermoDatabase()
        path = os.path.join(os.path.dirname(rmgpy.__file__),'data/test_data/')
        database2.load(os.path.join(path, 'thermo'), depository = False)

        #load up the thermo radical database as a test
        radGroup = database2.groups['radical']

        #use root as removed groups parent, which should be an LogicOr node
        root = radGroup.top[0]
        #use group to remove as
        groupToRemove = radGroup.entries['RJ']
        children = groupToRemove.children

        #remove the group
        radGroup.removeGroup(groupToRemove)

        #afterwards groupToRemove should not be in the database or root's children
        self.assertFalse(groupToRemove in radGroup.entries.values())
        self.assertFalse(groupToRemove in root.children)

        for child in children:
            #groupToRemove children should all be in roots item.component and children attribuetes
            self.assertTrue(child.label in root.item.components)
            self.assertTrue(child in root.children)
            #the children should all have root a their parent now
            self.assertTrue(child.parent is root)

        #Specific to ThermoDatabase, (above test apply to all base class Database)
        #we check that unicode entry.data pointers are correctly reassigned

        #if groupToRemove is a pointer and another node pointed to it, we copy
        #groupToRemove pointer
        self.assertTrue(radGroup.entries['OJ'].data is groupToRemove.data)

        #Remove an entry with a ThermoData object
        groupToRemove2 = radGroup.entries['CsJ']
        radGroup.removeGroup(groupToRemove2)
        #If groupToRemove was a data object, we point toward parent instead
        self.assertTrue(radGroup.entries['RJ2_triplet'].data == groupToRemove2.parent.label)
        #If the parent pointed toward groupToRemove, we need should have copied data object
        Tlist=[300, 400, 500, 600, 800, 1000, 1500]
        self.assertFalse(isinstance(groupToRemove2.parent.data, basestring))
        self.assertTrue(groupToRemove2.parent.data.getEnthalpy(298) == groupToRemove2.data.getEnthalpy(298))
        self.assertTrue(groupToRemove2.parent.data.getEntropy(298) == groupToRemove2.data.getEntropy(298))
        self.assertFalse(False in [groupToRemove2.parent.data.getHeatCapacity(x) == groupToRemove2.data.getHeatCapacity(x) for x in Tlist])

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))

