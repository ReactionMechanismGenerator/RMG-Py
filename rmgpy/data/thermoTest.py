#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import unittest
from external.wip import work_in_progress

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase, database
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import Species
from rmgpy.data.thermo import *
from rmgpy.molecule.molecule import Molecule
import rmgpy

################################################################################

class TestThermoDatabaseLoading(unittest.TestCase):

    def testFailingLoadsThermoLibraries(self):

        database = ThermoDatabase()
        libraries = ['primaryThermoLibrary', 'GRI-Mech3.0', 'I am a library not existing in official RMG']
        path = os.path.join(settings['database.directory'], 'thermo')
        
        with self.assertRaises(Exception):
            database.loadLibraries(os.path.join(path, 'libraries'), libraries)

class TestThermoDatabase(unittest.TestCase):
    """
    Contains unit tests of the ThermoDatabase class.
    """
    # Only load these once to save time
    database = ThermoDatabase()
    database.load(os.path.join(settings['database.directory'], 'thermo'))

    
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



    def testParseThermoComments(self):
        """
        Test that the ThermoDatabase.extractSourceFromComments function works properly
        on various thermo comments.
        """
        from rmgpy.thermo import NASA, NASAPolynomial
        # Pure group additivity thermo
        propane = Species(index=3, label="Propane", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[3.05257,0.0125099,3.79386e-05,-5.12022e-08,1.87065e-11,-14454.2,10.0672], Tmin=(100,'K'), Tmax=(986.57,'K')),
         NASAPolynomial(coeffs=[5.91316,0.0218763,-8.17661e-06,1.49855e-09,-1.05991e-13,-16038.9,-8.86555], Tmin=(986.57,'K'), Tmax=(5000,'K'))],
         Tmin=(100,'K'), Tmax=(5000,'K'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + gauche(Cs(CsCsRR)) + other(R) + group(Cs-CsHHH) + gauche(Cs(Cs(CsRR)RRR)) + other(R) + group(Cs-CsHHH) + gauche(Cs(Cs(CsRR)RRR)) + other(R)"""), molecule=[Molecule(SMILES="CCC")])
        
        source = self.database.extractSourceFromComments(propane)
        self.assertTrue('GAV' in source, 'Should have found that propane thermo source is GAV.')
        self.assertEqual(len(source['GAV']['group']), 2)
        self.assertEqual(len(source['GAV']['other']), 1)
        self.assertEqual(len(source['GAV']['gauche']), 2)

        # Pure library thermo
        dipk = Species(index=1, label="DIPK", thermo=
        NASA(polynomials=[NASAPolynomial(coeffs=[3.35002,0.017618,-2.46235e-05,1.7684e-08,-4.87962e-12,35555.7,5.75335], Tmin=(100,'K'), Tmax=(888.28,'K')), 
        NASAPolynomial(coeffs=[6.36001,0.00406378,-1.73509e-06,5.05949e-10,-4.49975e-14,35021,-8.41155], Tmin=(888.28,'K'), Tmax=(5000,'K'))],
         Tmin=(100,'K'), Tmax=(5000,'K'), comment="""Thermo library: DIPK"""), molecule=[Molecule(SMILES="CC(C)C(=O)C(C)C")])
        
        source = self.database.extractSourceFromComments(dipk)
        self.assertTrue('Library' in source)
        
        # Mixed library and HBI thermo
        dipk_rad = Species(index=4, label="R_tert", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[2.90061,0.0298018,-7.06268e-05,6.9636e-08,-2.42414e-11,54431,5.44492], Tmin=(100,'K'), Tmax=(882.19,'K')), 
        NASAPolynomial(coeffs=[6.70999,0.000201027,6.65617e-07,-7.99543e-11,4.08721e-15,54238.6,-9.73662], Tmin=(882.19,'K'), Tmax=(5000,'K'))],
         Tmin=(100,'K'), Tmax=(5000,'K'), comment="""Thermo library: DIPK + radical(C2CJCHO)"""), molecule=[Molecule(SMILES="C[C](C)C(=O)C(C)C"), Molecule(SMILES="CC(C)=C([O])C(C)C")])
        
        source = self.database.extractSourceFromComments(dipk_rad)
        self.assertTrue('Library' in source)
        self.assertTrue('GAV' in source)
        self.assertEqual(len(source['GAV']['radical']),1)

        # Pure QM thermo
        cineole = Species(index=6, label="Cineole", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-0.324129,0.0619667,9.71008e-05,-1.60598e-07,6.28285e-11,-38699.9,29.3686], Tmin=(100,'K'), Tmax=(985.52,'K')),
         NASAPolynomial(coeffs=[20.6043,0.0586913,-2.22152e-05,4.19949e-09,-3.06046e-13,-46791,-91.4152], Tmin=(985.52,'K'), Tmax=(5000,'K'))],
         Tmin=(100,'K'), Tmax=(5000,'K'), comment="""QM MopacMolPM3 calculation attempt 1"""), molecule=[Molecule(SMILES="CC12CCC(CC1)C(C)(C)O2")])
        
        source = self.database.extractSourceFromComments(cineole)
        self.assertTrue('QM' in source)
        
        # Mixed QM and HBI thermo
        cineole_rad = Species(index=7, label="CineoleRad", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-0.2897,0.0627717,8.63299e-05,-1.47868e-07,5.81665e-11,-14017.6,31.0266], Tmin=(100,'K'), Tmax=(988.76,'K')),
         NASAPolynomial(coeffs=[20.4836,0.0562555,-2.13903e-05,4.05725e-09,-2.96023e-13,-21915,-88.1205], Tmin=(988.76,'K'), Tmax=(5000,'K'))],
         Tmin=(100,'K'), Tmax=(5000,'K'), comment="""QM MopacMolPM3 calculation attempt 1 + radical(Cs_P)"""), molecule=[Molecule(SMILES="[CH2]C12CCC(CC1)C(C)(C)O2")])
        
        
        source = self.database.extractSourceFromComments(cineole_rad)
        self.assertTrue('QM' in source)
        self.assertTrue('GAV' in source)
        self.assertEqual(len(source['GAV']['radical']),1)
        
        
        # No thermo comments
        other = Species(index=7, label="CineoleRad", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-0.2897,0.0627717,8.63299e-05,-1.47868e-07,5.81665e-11,-14017.6,31.0266], Tmin=(100,'K'), Tmax=(988.76,'K')),
         NASAPolynomial(coeffs=[20.4836,0.0562555,-2.13903e-05,4.05725e-09,-2.96023e-13,-21915,-88.1205], Tmin=(988.76,'K'), Tmax=(5000,'K'))],
         Tmin=(100,'K'), Tmax=(5000,'K'), ), molecule=[Molecule(SMILES="[CH2]C12CCC(CC1)C(C)(C)O2")])
        
        # Function should complain if there's no thermo comments
        self.assertRaises(self.database.extractSourceFromComments(cineole_rad)) 
        
        
        # Check a dummy species that has plus and minus thermo group contributions
        polycyclic = Species(index=7, label="dummy", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-0.2897,0.0627717,8.63299e-05,-1.47868e-07,5.81665e-11,-14017.6,31.0266], Tmin=(100,'K'), Tmax=(988.76,'K')),
         NASAPolynomial(coeffs=[20.4836,0.0562555,-2.13903e-05,4.05725e-09,-2.96023e-13,-21915,-88.1205], Tmin=(988.76,'K'), Tmax=(5000,'K'))],
         Tmin=(100,'K'), Tmax=(5000,'K'),  comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) - ring(Benzene)"""), molecule=[Molecule(SMILES="[CH2]C12CCC(CC1)C(C)(C)O2")])
   
        source = self.database.extractSourceFromComments(polycyclic)
        self.assertTrue('GAV' in source)
        self.assertEqual(source['GAV']['ring'][0][1],-1)  # the weight of benzene contribution should be -1
        self.assertEqual(source['GAV']['group'][0][1],2)  # weight of the group(Cs-CsCsHH) conbtribution should be 2
        
    def testSpeciesThermoGenerationHBILibrary(self):
        """Test thermo generation for species objects.

        Ensure that molecule list is only reordered, and not changed after matching library value"""
        spec = Species().fromSMILES('C[CH]c1ccccc1')
        spec.generateResonanceIsomers()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.getThermoData(spec)

        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertTrue('library' in thermo.comment, 'Thermo not found from library, test purpose not fulfilled.')

    def testSpeciesThermoGenerationHBIGAV(self):
        """Test thermo generation for species objects.

        Ensure that molecule list is only reordered, and not changed after group additivity"""
        spec = Species().fromSMILES('CCC[CH]c1ccccc1')
        spec.generateResonanceIsomers()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.getThermoData(spec)

        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertTrue('group additivity' in thermo.comment, 'Thermo not found from GAV, test purpose not fulfilled.')


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
        
        expected_matchedPolyringsLabels = ['s3_5_5_ane']
        expected_matchedPolyrings = [self.database.groups['polycyclic'].entries[label] for label in expected_matchedPolyringsLabels]

        self.assertEqual(set(polycyclicGroups), set(expected_matchedPolyrings))

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

    def testIsRingPartialMatched(self):
        
        # create testing molecule
        smiles = 'C1CC2CCCC3CCCC(C1)C23'
        mol = Molecule().fromSMILES(smiles)
        polyring = [atom for atom in mol.atoms if atom.isNonHydrogen()]

        # create matched group
        matched_group = self.database.groups['polycyclic'].entries['PolycyclicRing'].item
        
        # test
        self.assertTrue(isRingPartialMatched(polyring, matched_group))

    def testAddRingCorrectionThermoDataFromTreeForExistingTricyclic(self):

        # create testing molecule: C1CC2C3CCC(C3)C2C1
        # this tricyclic molecule is already in polycyclic database
        # so algorithm should give complete match: s2-3_5_5_5_ane
        smiles = 'C1CC2C3CCC(C3)C2C1'
        mol = Molecule().fromSMILES(smiles)
        polyring = mol.getDisparateRings()[1][0]

        poly_groups = self.database.groups['polycyclic']
        _, matched_entry, _ = self.database._ThermoDatabase__addRingCorrectionThermoDataFromTree(None, poly_groups, mol, polyring)

        self.assertEqual(matched_entry.label, 's2-3_5_5_5_ane')

    def testAddPolyRingCorrectionThermoDataFromHeuristicUsingPyrene(self):

        # create testing molecule: Pyrene with two ring of aromatic version
        # the other two ring of kekulized version
        #
        # creating it seems not natural in RMG, that's because
        # RMG cannot parse the adjacencyList of that isomer correctly
        # so here we start with pyrene radical and get the two aromatic ring isomer
        # then saturate it.
        smiles = 'C1C=C2C=CC=C3C=CC4=CC=CC=1C4=C23'
        spe = Species().fromSMILES(smiles)
        spe.generateResonanceIsomers()
        mols = []
        for mol in spe.molecule:
            sssr0 = mol.getSmallestSetOfSmallestRings()
            aromaticRingNum = 0
            for sr0 in sssr0:
                sr0mol = Molecule(atoms=sr0)
                if isAromaticRing(sr0mol):
                    aromaticRingNum += 1
            if aromaticRingNum == 2:
                mols.append(mol)
        
        ringGroupLabels = []
        polycyclicGroupLabels = []
        for mol in mols:
            polyring = mol.getDisparateRings()[1][0]

            thermoData = ThermoData(
                Tdata = ([300,400,500,600,800,1000,1500],"K"),
                Cpdata = ([0.0,0.0,0.0,0.0,0.0,0.0,0.0],"J/(mol*K)"),
                H298 = (0.0,"kJ/mol"),
                S298 = (0.0,"J/(mol*K)"),
            )

            self.database._ThermoDatabase__addPolyRingCorrectionThermoDataFromHeuristic(
                thermoData, polyring)

            ringGroups, polycyclicGroups = self.database.getRingGroupsFromComments(thermoData)

            ringGroupLabels += [ringGroup.label for ringGroup in ringGroups]
            polycyclicGroupLabels += [polycyclicGroup.label for polycyclicGroup in polycyclicGroups]

        self.assertIn('Benzene', ringGroupLabels)
        self.assertIn('Cyclohexene', ringGroupLabels)
        self.assertIn('s2_6_6_ben_ene_1', polycyclicGroupLabels)
        self.assertIn('s2_6_6_diene_2_7', polycyclicGroupLabels)

    def testAddPolyRingCorrectionThermoDataFromHeuristicUsingAromaticTricyclic(self):

        # create testing molecule
        #
        # creating it seems not natural in RMG, that's because
        # RMG cannot parse the adjacencyList of that isomer correctly
        # so here we start with kekulized version and generateResonanceIsomers
        # and pick the one with two aromatic rings
        smiles = 'C1=CC2C=CC=C3C=CC(=C1)C=23'
        spe = Species().fromSMILES(smiles)
        spe.generateResonanceIsomers()
        for mol in spe.molecule:
            sssr0 = mol.getSmallestSetOfSmallestRings()
            aromaticRingNum = 0
            for sr0 in sssr0:
                sr0mol = Molecule(atoms=sr0)
                if isAromaticRing(sr0mol):
                    aromaticRingNum += 1
            if aromaticRingNum == 2:
                break
        
        # extract polyring from the molecule
        polyring = mol.getDisparateRings()[1][0]

        thermoData = ThermoData(
            Tdata = ([300,400,500,600,800,1000,1500],"K"),
            Cpdata = ([0.0,0.0,0.0,0.0,0.0,0.0,0.0],"J/(mol*K)"),
            H298 = (0.0,"kJ/mol"),
            S298 = (0.0,"J/(mol*K)"),
        )

        self.database._ThermoDatabase__addPolyRingCorrectionThermoDataFromHeuristic(
            thermoData, polyring)

        ringGroups, polycyclicGroups = self.database.getRingGroupsFromComments(thermoData)

        ringGroupLabels = [ringGroup.label for ringGroup in ringGroups]
        polycyclicGroupLabels = [polycyclicGroup.label for polycyclicGroup in polycyclicGroups]

        self.assertIn('Benzene', ringGroupLabels)
        self.assertIn('Cyclopentene', ringGroupLabels)
        self.assertIn('s2_5_6_indene', polycyclicGroupLabels)
        self.assertIn('s2_6_6_naphthalene', polycyclicGroupLabels)

    def testAddPolyRingCorrectionThermoDataFromHeuristicUsingAlkaneTricyclic(self):

        # create testing molecule
        smiles = 'C1CC2CCCC3C(C1)C23'
        mol = Molecule().fromSMILES(smiles)
        
        # extract polyring from the molecule
        polyring = mol.getDisparateRings()[1][0]

        thermoData = ThermoData(
            Tdata = ([300,400,500,600,800,1000,1500],"K"),
            Cpdata = ([0.0,0.0,0.0,0.0,0.0,0.0,0.0],"J/(mol*K)"),
            H298 = (0.0,"kJ/mol"),
            S298 = (0.0,"J/(mol*K)"),
        )

        self.database._ThermoDatabase__addPolyRingCorrectionThermoDataFromHeuristic(
            thermoData, polyring)

        ringGroups, polycyclicGroups = self.database.getRingGroupsFromComments(thermoData)

        ringGroupLabels = [ringGroup.label for ringGroup in ringGroups]
        polycyclicGroupLabels = [polycyclicGroup.label for polycyclicGroup in polycyclicGroups]

        self.assertIn('Cyclohexane', ringGroupLabels)
        self.assertIn('Cyclopropane', ringGroupLabels)
        self.assertIn('s2_6_6_ane', polycyclicGroupLabels)
        self.assertIn('s2_3_6_ane', polycyclicGroupLabels)

    def testGetBicyclicCorrectionThermoDataFromHeuristic1(self):
        """
        Test bicyclic correction estimated properly from heuristic formula
        """
        smiles = 'C1=CCC2C1=C2'
        mol = Molecule().fromSMILES(smiles)

        # extract polyring from the molecule
        polyring = mol.getDisparateRings()[1][0]

        thermoData = self.database.getBicyclicCorrectionThermoDataFromHeuristic(polyring)

        ringGroups, polycyclicGroups = self.database.getRingGroupsFromComments(thermoData)

        ringGroupLabels = [ringGroup.label for ringGroup in ringGroups]
        polycyclicGroupLabels = [polycyclicGroup.label for polycyclicGroup in polycyclicGroups]

        self.assertIn('Cyclopentane', ringGroupLabels)
        self.assertIn('Cyclopropane', ringGroupLabels)
        self.assertIn('Cyclopentene', ringGroupLabels)
        self.assertIn('Cyclopropene', ringGroupLabels)
        self.assertIn('s2_3_5_ane', polycyclicGroupLabels)

    def testGetBicyclicCorrectionThermoDataFromHeuristic2(self):
        """
        Test bicyclic correction estimated properly from heuristic formula
        """
        smiles = 'C1=CCC2=C1C2'
        mol = Molecule().fromSMILES(smiles)

        # extract polyring from the molecule
        polyring = mol.getDisparateRings()[1][0]

        thermoData = self.database.getBicyclicCorrectionThermoDataFromHeuristic(polyring)

        ringGroups, polycyclicGroups = self.database.getRingGroupsFromComments(thermoData)

        ringGroupLabels = [ringGroup.label for ringGroup in ringGroups]
        polycyclicGroupLabels = [polycyclicGroup.label for polycyclicGroup in polycyclicGroups]

        self.assertIn('Cyclopentane', ringGroupLabels)
        self.assertIn('Cyclopropane', ringGroupLabels)
        self.assertIn('Cyclopentadiene', ringGroupLabels)
        self.assertIn('Cyclopropene', ringGroupLabels)
        self.assertIn('s2_3_5_ane', polycyclicGroupLabels)
        
class TestMolecularManipulationInvolvedInThermoEstimation(unittest.TestCase):
    """
    Contains unit tests for methods of molecular manipulations for thermo estimation
    """

    def testConvertRingToSubMolecule(self):

        # list out testing moleculess
        smiles1 = 'C1CCCCC1'
        smiles2 = 'C1CCC2CCCCC2C1'
        smiles3 = 'C1CC2CCCC3CCCC(C1)C23'
        mol1 = Molecule().fromSMILES(smiles1)
        mol2 = Molecule().fromSMILES(smiles2)
        mol3 = Molecule().fromSMILES(smiles3)

        # get ring structure by only extracting non-hydrogens
        ring1 = [atom for atom in mol1.atoms if atom.isNonHydrogen()]
        ring2 = [atom for atom in mol2.atoms if atom.isNonHydrogen()]
        ring3 = [atom for atom in mol3.atoms if atom.isNonHydrogen()]

        # convert to submolecules
        submol1, _ = convertRingToSubMolecule(ring1)
        submol2, _ = convertRingToSubMolecule(ring2)
        submol3, _ = convertRingToSubMolecule(ring3)

        # test against expected submolecules
        self.assertEqual(len(submol1.atoms), 6)
        self.assertEqual(len(submol2.atoms), 10)
        self.assertEqual(len(submol3.atoms), 13)

        bonds1 = []
        for atom in submol1.atoms:
            for bondAtom, bond in atom.edges.iteritems():
                if bond not in bonds1:
                    bonds1.append(bond)

        bonds2 = []
        for atom in submol2.atoms:
            for bondAtom, bond in atom.edges.iteritems():
                if bond not in bonds2:
                    bonds2.append(bond)

        bonds3 = []
        for atom in submol3.atoms:
            for bondAtom, bond in atom.edges.iteritems():
                if bond not in bonds3:
                    bonds3.append(bond)

        self.assertEqual(len(bonds1), 6)
        self.assertEqual(len(bonds2), 11)
        self.assertEqual(len(bonds3), 15)

    def testToFailCombineTwoRingsIntoSubMolecule(self):
        """
        Test that if two non-overlapped rings lead to AssertionError
        """

        smiles1 = 'C1CCCCC1'
        smiles2 = 'C1CCCCC1'
        mol1 = Molecule().fromSMILES(smiles1)
        mol2 = Molecule().fromSMILES(smiles2)

        ring1 = [atom for atom in mol1.atoms if atom.isNonHydrogen()]
        ring2 = [atom for atom in mol2.atoms if atom.isNonHydrogen()]
        with self.assertRaises(AssertionError):
            combined = combineTwoRingsIntoSubMolecule(ring1, ring2)

    def testCombineTwoRingsIntoSubMolecule(self):

        # create testing molecule
        smiles1 = 'C1CCC2CCCCC2C1'
        mol1 = Molecule().fromSMILES(smiles1)

        # get two SSSRs
        SSSR = mol1.getSmallestSetOfSmallestRings()
        ring1 = SSSR[0]
        ring2 = SSSR[1]

        # combine two rings into submolecule
        submol, _ = combineTwoRingsIntoSubMolecule(ring1, ring2)

        self.assertEqual(len(submol.atoms), 10)

        bonds = []
        for atom in submol.atoms:
            for bondAtom, bond in atom.edges.iteritems():
                if bond not in bonds:
                    bonds.append(bond)

        self.assertEqual(len(bonds), 11)

    def testIsAromaticRing(self):

        # create testing rings
        smiles1 = 'C1CCC1'
        smiles2 = 'C1CCCCC1'
        adj3 = """1  C u0 p0 c0 {2,B} {6,B} {7,S}
2  C u0 p0 c0 {1,B} {3,B} {8,S}
3  C u0 p0 c0 {2,B} {4,B} {9,S}
4  C u0 p0 c0 {3,B} {5,B} {10,S}
5  C u0 p0 c0 {4,B} {6,B} {11,S}
6  C u0 p0 c0 {1,B} {5,B} {12,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {2,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {4,S}
11 H u0 p0 c0 {5,S}
12 H u0 p0 c0 {6,S}
        """
        mol1 = Molecule().fromSMILES(smiles1)
        mol2 = Molecule().fromSMILES(smiles2)
        mol3 = Molecule().fromAdjacencyList(adj3)
        ring1mol = Molecule(atoms=[atom for atom in mol1.atoms if atom.isNonHydrogen()])
        ring2mol = Molecule(atoms=[atom for atom in mol2.atoms if atom.isNonHydrogen()])
        ring3mol = Molecule(atoms=[atom for atom in mol3.atoms if atom.isNonHydrogen()])

        # check with expected results
        self.assertEqual(isAromaticRing(ring1mol), False)
        self.assertEqual(isAromaticRing(ring2mol), False)
        self.assertEqual(isAromaticRing(ring3mol), True)

    def testFindAromaticBondsFromSubMolecule(self):

        smiles = "C1=CC=C2C=CC=CC2=C1"
        spe = Species().fromSMILES(smiles)
        spe.generateResonanceIsomers()
        mol = spe.molecule[1]

        # get two SSSRs
        SSSR = mol.getSmallestSetOfSmallestRings()
        ring1 = SSSR[0]
        ring2 = SSSR[1]

        # create two testing submols
        submol1 = Molecule(atoms=ring1)
        submol2 = Molecule(atoms=ring2)

        # check with expected results
        self.assertEqual(len(findAromaticBondsFromSubMolecule(submol1)), 6)
        self.assertEqual(len(findAromaticBondsFromSubMolecule(submol2)), 6)

    def testBicyclicDecompositionForPolyringUsingPyrene(self):

        # create testing molecule: Pyrene with two ring of aromatic version
        # the other two ring of kekulized version
        #
        # creating it seems not natural in RMG, that's because
        # RMG cannot parse the adjacencyList of that isomer correctly
        # so here we start with pyrene radical and get the two aromatic ring isomer
        # then saturate it.
        smiles = 'C1C=C2C=CC=C3C=CC4=CC=CC=1C4=C23'
        spe = Species().fromSMILES(smiles)
        spe.generateResonanceIsomers()
        for mol in spe.molecule:
            sssr0 = mol.getSmallestSetOfSmallestRings()
            aromaticRingNum = 0
            for sr0 in sssr0:
                sr0mol = Molecule(atoms=sr0)
                if isAromaticRing(sr0mol):
                    aromaticRingNum += 1
            if aromaticRingNum == 2:
                break

        # extract polyring from the molecule
        polyring = mol.getDisparateRings()[1][0]

        bicyclicList, ringOccurancesDict = bicyclicDecompositionForPolyring(polyring)

        # 1st test: number of cores
        self.assertEqual(len(bicyclicList), 5)

        # 2nd test: ringOccurancesDict
        ringInCoreOccurances = sorted(ringOccurancesDict.values())
        expectedRingInCoreOccurances = [2, 2, 3, 3]
        self.assertEqual(ringInCoreOccurances, expectedRingInCoreOccurances)

        # 3rd test: size of each bicyclic core
        bicyclicSizes = sorted([len(bicyclic.atoms) for bicyclic in bicyclicList])
        expectedBicyclicSizes = [10, 10, 10, 10, 10]
        self.assertEqual(bicyclicSizes, expectedBicyclicSizes)

        # 4th test: bond info for members of each core
        aromaticBondNumInBicyclics = []
        for bicyclic in bicyclicList:
            aromaticBondNum = len(findAromaticBondsFromSubMolecule(bicyclic))
            aromaticBondNumInBicyclics.append(aromaticBondNum)
        aromaticBondNumInBicyclics = sorted(aromaticBondNumInBicyclics)
        expectedAromaticBondNumInBicyclics = [0, 6, 6, 6, 6]
        self.assertEqual(aromaticBondNumInBicyclics, expectedAromaticBondNumInBicyclics)

    def testBicyclicDecompositionForPolyringUsingAromaticTricyclic(self):

        # create testing molecule
        #
        # creating it seems not natural in RMG, that's because
        # RMG cannot parse the adjacencyList of that isomer correctly
        # so here we start with kekulized version and generateResonanceIsomers
        # and pick the one with two aromatic rings
        smiles = 'C1=CC2C=CC=C3C=CC(=C1)C=23'
        spe = Species().fromSMILES(smiles)
        spe.generateResonanceIsomers()
        for mol in spe.molecule:
            sssr0 = mol.getSmallestSetOfSmallestRings()
            aromaticRingNum = 0
            for sr0 in sssr0:
                sr0mol = Molecule(atoms=sr0)
                if isAromaticRing(sr0mol):
                    aromaticRingNum += 1
            if aromaticRingNum == 2:
                break
        
        # extract polyring from the molecule
        polyring = mol.getDisparateRings()[1][0]

        bicyclicList, ringOccurancesDict = bicyclicDecompositionForPolyring(polyring)

        # 1st test: number of cores
        self.assertEqual(len(bicyclicList), 3)

        # 2nd test: ringOccurancesDict
        ringInCoreOccurances = sorted(ringOccurancesDict.values())
        expectedRingInCoreOccurances = [2, 2, 2]
        self.assertEqual(ringInCoreOccurances, expectedRingInCoreOccurances)

        # 3rd test: size of each bicyclic core
        bicyclicSizes = sorted([len(bicyclic.atoms) for bicyclic in bicyclicList])
        expectedBicyclicSizes = [9, 9, 10]
        self.assertEqual(bicyclicSizes, expectedBicyclicSizes)

        # 4th test: bond info for members of each core
        aromaticBondNumInBicyclics = []
        for bicyclic in bicyclicList:
            aromaticBondNum = len(findAromaticBondsFromSubMolecule(bicyclic))
            aromaticBondNumInBicyclics.append(aromaticBondNum)
        aromaticBondNumInBicyclics = sorted(aromaticBondNumInBicyclics)
        expectedAromaticBondNumInBicyclics = [6, 6, 11]
        self.assertEqual(aromaticBondNumInBicyclics, expectedAromaticBondNumInBicyclics)


    def testBicyclicDecompositionForPolyringUsingAlkaneTricyclic(self):

        # create testing molecule
        smiles = 'C1CC2CCCC3C(C1)C23'
        mol = Molecule().fromSMILES(smiles)
        
        # extract polyring from the molecule
        polyring = mol.getDisparateRings()[1][0]

        bicyclicList, ringOccurancesDict = bicyclicDecompositionForPolyring(polyring)

        # 1st test: number of cores
        self.assertEqual(len(bicyclicList), 3)

        # 2nd test: ringOccurancesDict
        ringInCoreOccurances = sorted(ringOccurancesDict.values())
        expectedRingInCoreOccurances = [2, 2, 2]
        self.assertEqual(ringInCoreOccurances, expectedRingInCoreOccurances)

        # 3rd test: size of each bicyclic core
        bicyclicSizes = sorted([len(bicyclic.atoms) for bicyclic in bicyclicList])
        expectedBicyclicSizes = [7, 7, 10]
        self.assertEqual(bicyclicSizes, expectedBicyclicSizes)

        # 4th test: bond info for members of each core
        aromaticBondNumInBicyclics = []
        for bicyclic in bicyclicList:
            aromaticBondNum = len(findAromaticBondsFromSubMolecule(bicyclic))
            aromaticBondNumInBicyclics.append(aromaticBondNum)
        aromaticBondNumInBicyclics = sorted(aromaticBondNumInBicyclics)
        expectedAromaticBondNumInBicyclics = [0, 0, 0]
        self.assertEqual(aromaticBondNumInBicyclics, expectedAromaticBondNumInBicyclics)

    def testSplitBicyclicIntoSingleRings(self):
        """
        Test bicyclic molecule can be divided into individual rings properly
        """
        smiles1 = 'C1=CCC2C1=C2'
        mol1 = Molecule().fromSMILES(smiles1)
        bicyclic1 = mol1.getDisparateRings()[1][0]

        bicyclic1_submol = convertRingToSubMolecule(bicyclic1)[0]
        single_ring_submols1 = splitBicyclicIntoSingleRings(bicyclic1_submol)
        self.assertEqual(len(single_ring_submols1), 2)

        single_ring_submol_a, single_ring_submol_b = sorted(single_ring_submols1, \
                                key=lambda submol: len(submol.atoms))

        expected_submol_a = Molecule().fromSMILES('C1=CC1')
        atomsToRemove = []
        for atom in expected_submol_a.atoms:
            if len(atom.bonds) == 1: 
                atomsToRemove.append(atom)

        for atom in atomsToRemove:
            expected_submol_a.removeAtom(atom)
        expected_submol_a.updateConnectivityValues()

        expected_submol_b = Molecule().fromSMILES('C1=CCCC1')
        atomsToRemove = []
        for atom in expected_submol_b.atoms:
            if len(atom.bonds) == 1: 
                atomsToRemove.append(atom)

        for atom in atomsToRemove:
            expected_submol_b.removeAtom(atom)
        expected_submol_b.updateConnectivityValues()


        self.assertTrue(single_ring_submol_a.isIsomorphic(expected_submol_a))
        self.assertTrue(single_ring_submol_b.isIsomorphic(expected_submol_b))


        smiles2 = 'C1=CCC2=C1C2'
        mol2 = Molecule().fromSMILES(smiles2)
        bicyclic2 = mol2.getDisparateRings()[1][0]

        bicyclic2_submol = convertRingToSubMolecule(bicyclic2)[0]
        single_ring_submols2 = splitBicyclicIntoSingleRings(bicyclic2_submol)
        self.assertEqual(len(single_ring_submols2), 2)

        single_ring_submol_a, single_ring_submol_b = sorted(single_ring_submols2, \
                                key=lambda submol: len(submol.atoms))

        expected_submol_a = Molecule().fromSMILES('C1=CC1')
        # remove hydrogen
        atomsToRemove = []
        for atom in expected_submol_a.atoms:
            if atom.isHydrogen(): 
                atomsToRemove.append(atom)

        for atom in atomsToRemove:
            expected_submol_a.removeAtom(atom)
        expected_submol_a.updateConnectivityValues()

        expected_submol_b = Molecule().fromSMILES('C1=CC=CC1')
        # remove hydrogen
        atomsToRemove = []
        for atom in expected_submol_b.atoms:
            if atom.isHydrogen(): 
                atomsToRemove.append(atom)

        for atom in atomsToRemove:
            expected_submol_b.removeAtom(atom)
        expected_submol_b.updateConnectivityValues()

        self.assertTrue(single_ring_submol_a.isIsomorphic(expected_submol_a))
        self.assertTrue(single_ring_submol_b.isIsomorphic(expected_submol_b))

    def testSaturateRingBonds1(self):
        """
        Test unsaturated bonds can be saturated properly
        """
        smiles = 'C1=CCC2=C1C2'
        mol = Molecule().fromSMILES(smiles)
        ring_submol = convertRingToSubMolecule(mol.getDisparateRings()[1][0])[0]

        saturated_ring_submol, alreadySaturated = saturateRingBonds(ring_submol)

        expected_saturated_ring_submol = Molecule().fromSMILES('C1CCC2C1C2')
        # remove hydrogen
        atomsToRemove = []
        for atom in expected_saturated_ring_submol.atoms:
            if atom.isHydrogen(): 
                atomsToRemove.append(atom)

        for atom in atomsToRemove:
            expected_saturated_ring_submol.removeAtom(atom)
        
        expected_saturated_ring_submol.updateConnectivityValues()

        self.assertFalse(alreadySaturated)
        self.assertEqual(saturated_ring_submol.multiplicity, \
                            expected_saturated_ring_submol.multiplicity)
        self.assertTrue(saturated_ring_submol.isIsomorphic(expected_saturated_ring_submol))

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))

