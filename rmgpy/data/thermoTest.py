#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import os
import unittest
from external.wip import work_in_progress
from nose.plugins.attrib import attr

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase, database
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import Species
from rmgpy.data.thermo import *
from rmgpy.molecule.molecule import Molecule
import rmgpy

################################################################################

def setUpModule():
    """A function that is run ONCE before all unit tests in this module."""
    global database
    database = RMGDatabase()
    database.loadThermo(os.path.join(settings['database.directory'], 'thermo'))

def tearDownModule():
    """A function that is run ONCE after all unit tests in this module."""
    from rmgpy.data import rmg
    rmg.database = None

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
    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        self.database = database.thermo

        self.databaseWithoutLibraries = ThermoDatabase()
        self.databaseWithoutLibraries.load(os.path.join(settings['database.directory'], 'thermo'),libraries = [])

    def testPickle(self):
        """
        Test that a ThermoDatabase object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        thermodb0 = cPickle.loads(cPickle.dumps(self.database))
        
        self.assertEqual(thermodb0.libraryOrder, self.database.libraryOrder)
        self.assertEqual(sorted(thermodb0.depository.keys()),
                         sorted(self.database.depository.keys()))

        self.assertEqual(sorted(thermodb0.libraries.keys()),
                         sorted(self.database.libraries.keys()))
        self.assertEqual(sorted(thermodb0.groups.keys()),
                         sorted(self.database.groups.keys()))

        for key, depository0 in thermodb0.depository.iteritems():
            depository = self.database.depository[key]
            self.assertTrue(type(depository0), type(depository))
            self.assertEqual(sorted(depository0.entries.keys()), sorted(depository.entries.keys()))

        for key, library0 in thermodb0.libraries.iteritems():
            library = self.database.libraries[key]
            self.assertTrue(type(library0), type(library))
            self.assertEqual(sorted(library0.entries.keys()), sorted(library.entries.keys()))

        for key, group0 in thermodb0.groups.iteritems():
            group = self.database.groups[key]
            self.assertTrue(type(group0), type(group))
            self.assertEqual(sorted(group0.entries.keys()), sorted(group.entries.keys()))

    def testSymmetryAddedByGetThermoData(self):
        """
        Test that `getThermoData` properly accounts for symmetry in thermo
        by comping with the method `estimateThermoViaGroupAdditivity`
        """
        
        spc = Species(molecule=[Molecule().fromSMILES('C[CH]C=CC')])
        
        thermoWithSym = self.databaseWithoutLibraries.getThermoData(spc)
        thermoWithoutSym = self.databaseWithoutLibraries.estimateThermoViaGroupAdditivity(spc.molecule[0])
        
        symmetryNumber = spc.getSymmetryNumber()
        self.assertNotEqual(symmetryNumber, spc.molecule[0].getSymmetryNumber(),
                            'For this test to be robust, species symmetry ({}) and molecule symmetry ({}) must be different'.format(symmetryNumber, spc.molecule[0].getSymmetryNumber()))
        
        symmetryContributionToEntropy = - constants.R * math.log(symmetryNumber)
        
        self.assertAlmostEqual(thermoWithSym.getEntropy(298.), 
                               thermoWithoutSym.getEntropy(298.) + symmetryContributionToEntropy, 
                               'The symmetry contribution is wrong {:.3f} /= {:.3f} + {:.3f}'.format(thermoWithSym.getEntropy(298.), thermoWithoutSym.getEntropy(298.), symmetryContributionToEntropy))
        
    def testSymmetryContributionRadicals(self):
        """
        Test that the symmetry contribution is correctly added for radicals
        estimated via the HBI method.
        
        This is done by testing thermoData from a database and from group 
        additivity and ensuring they give the correct value. 
        """
        spc = Species(molecule=[Molecule().fromSMILES('[CH3]')])
        
        thermoData_lib = self.database.getThermoData(spc)
        
        
        
        thermoData_ga = self.databaseWithoutLibraries.getThermoData(spc)
        
        self.assertAlmostEqual(thermoData_lib.getEntropy(298.), thermoData_ga.getEntropy(298.), 0)

    def testParseThermoComments(self):
        """
        Test that the ThermoDatabase.extractSourceFromComments function works properly
        on various thermo comments.
        """
        from rmgpy.thermo import NASA, NASAPolynomial
        # Pure group additivity thermo.
        GAVspecies = Species(index=3, label="c1c(O)c(O)c(CC(C)CC)cc1", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-1.18833,0.11272,-4.26393e-05,-2.12017e-08,1.441e-11,-51642.9,38.8904], Tmin=(100,'K'), Tmax=(1078.35,'K')),
         NASAPolynomial(coeffs=[26.6057,0.0538434,-2.22538e-05,4.22393e-09,-3.00808e-13,-60208.4,-109.218], Tmin=(1078.35,'K'), Tmax=(5000,'K'))],
         Tmin=(100,'K'), Tmax=(5000,'K'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) +
         group(Cs-CbCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cb-Cs) + group(Cb-O2s) + group(Cb-O2s) + group(Cb-H) +
         group(Cb-H) + group(Cb-H) + group(O2s-CbH) + group(O2s-CbH) + longDistanceInteraction_cyclic(o_OH_OH) +
         longDistanceInteraction_cyclic(o_OH_OH) + ring(Benzene)"""), molecule=[Molecule(SMILES="c1c(O)c(O)c(CC(C)CC)cc1")])

        source = self.database.extractSourceFromComments(GAVspecies)
        self.assertTrue('GAV' in source, 'Should have found that the thermo source is GAV.')
        self.assertEqual(len(source['GAV']['group']), 8)
        self.assertEqual(len(source['GAV']['longDistanceInteraction_noncyclic']), 1)
        self.assertEqual(len(source['GAV']['longDistanceInteraction_cyclic']), 1)
        self.assertEqual(len(source['GAV']['ring']), 1)

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
        """Test thermo generation for species objects for HBI correction on library value.

        Ensure that molecule list is only reordered, and not changed after matching library value"""
        spec = Species().fromSMILES('C[CH]c1ccccc1')
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.getThermoData(spec)

        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertTrue('library' in thermo.comment, 'Thermo not found from library, test purpose not fulfilled.')

    def testSpeciesThermoGenerationHBIGAV(self):
        """Test thermo generation for species objects for HBI correction on group additivity value.

        Ensure that molecule list is only reordered, and not changed after group additivity"""
        spec = Species().fromSMILES('C[CH]c1ccccc1')
        spec.generate_resonance_structures()
        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.databaseWithoutLibraries.getThermoData(spec)

        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertTrue('group additivity' in thermo.comment, 'Thermo not found from GAV, test purpose not fulfilled.')

    def testSpeciesThermoGenerationLibrary(self):
        """Test thermo generation for species objects for library value.

        Ensure that the matched molecule is placed at the beginning of the list."""
        spec = Species().fromSMILES('c12ccccc1c(C=[CH])ccc2')
        arom = Molecule().fromAdjacencyList("""
multiplicity 2
1  C u0 p0 c0 {2,B} {3,B} {5,B}
2  C u0 p0 c0 {1,B} {4,B} {7,B}
3  C u0 p0 c0 {1,B} {6,B} {11,S}
4  C u0 p0 c0 {2,B} {8,B} {13,S}
5  C u0 p0 c0 {1,B} {9,B} {16,S}
6  C u0 p0 c0 {3,B} {10,B} {17,S}
7  C u0 p0 c0 {2,B} {10,B} {19,S}
8  C u0 p0 c0 {4,B} {9,B} {14,S}
9  C u0 p0 c0 {5,B} {8,B} {15,S}
10 C u0 p0 c0 {6,B} {7,B} {18,S}
11 C u0 p0 c0 {3,S} {12,D} {20,S}
12 C u1 p0 c0 {11,D} {21,S}
13 H u0 p0 c0 {4,S}
14 H u0 p0 c0 {8,S}
15 H u0 p0 c0 {9,S}
16 H u0 p0 c0 {5,S}
17 H u0 p0 c0 {6,S}
18 H u0 p0 c0 {10,S}
19 H u0 p0 c0 {7,S}
20 H u0 p0 c0 {11,S}
21 H u0 p0 c0 {12,S}
""")
        spec.generate_resonance_structures()

        self.assertTrue(arom.isIsomorphic(spec.molecule[1]))  # The aromatic structure should be the second one

        initial = list(spec.molecule)  # Make a copy of the list
        thermo = self.database.getThermoData(spec)

        self.assertEqual(len(initial), len(spec.molecule))
        self.assertEqual(set(initial), set(spec.molecule))
        self.assertTrue(arom.isIsomorphic(spec.molecule[0]))  # The aromatic structure should now be the first one
        self.assertTrue('library' in thermo.comment, 'Thermo not found from library, test purpose not fulfilled.')

    def testThermoEstimationNotAffectDatabase(self):

        poly_root = self.database.groups['polycyclic'].entries['PolycyclicRing']
        previous_enthalpy = poly_root.data.getEnthalpy(298)/4184.0
        smiles = 'C1C2CC1C=CC=C2'
        spec = Species().fromSMILES(smiles)
        spec.generate_resonance_structures()

        thermo_gav = self.database.getThermoDataFromGroups(spec)
        _, polycyclicGroups = self.database.getRingGroupsFromComments(thermo_gav)

        polycyclicGroupLabels = [polycyclicGroup.label for polycyclicGroup in polycyclicGroups]

        self.assertIn('PolycyclicRing', polycyclicGroupLabels)

        latter_enthalpy = poly_root.data.getEnthalpy(298)/4184.0

        self.assertAlmostEqual(previous_enthalpy, latter_enthalpy, 2)

    def test_getAllThermoData_fails_quietly(self):
        """Test that getAllThermoData doesn't break when GAV fails."""
        spec = Species().fromSMILES('[Ne]')

        # Check that GAV fails
        with self.assertRaises(DatabaseError):
            self.database.getThermoDataFromGroups(spec)

        # Check that getAllThermoData doesn't break
        thermo = self.database.getAllThermoData(spec)
        self.assertEqual(len(thermo), 1)

    def test_lowest_h298_for_resonance_structures(self):
        """Test that the thermo entry with the lowest H298 is selected for a species with resonance structurers"""

        smiles = '[C]#C[O]'  # has H298 ~= 640 kJ/mol; has resonance structure `[C]=C=O` with H298 ~= 380 kJ/mol
        spec = Species().fromSMILES(smiles)
        thermo_gav1 = self.database.getThermoDataFromGroups(spec)
        spec.generate_resonance_structures()
        thermo_gav2 = self.database.getThermoDataFromGroups(spec)
        self.assertTrue(thermo_gav2.getEnthalpy(298) < thermo_gav1.getEnthalpy(298),
                        msg="Did not select the molecule with the lowest H298 as a the thermo entry for [C]#C[O] / [C]=C=O")

        smiles = 'C=C[CH][O]'  # has H298 ~= 209 kJ/mol; has (a reactive) resonance structure `C=CC=O` with H298 ~= -67 kJ/mol
        spec = Species().fromSMILES(smiles)
        thermo_gav1 = self.database.getThermoDataFromGroups(spec)
        spec.generate_resonance_structures()
        thermo_gav2 = self.database.getThermoDataFromGroups(spec)
        self.assertTrue(thermo_gav2.getEnthalpy(298) < thermo_gav1.getEnthalpy(298),
                        msg="Did not select the molecule with the lowest H298 as a the thermo entry for C=C[CH][O] / C=CC=O")

class TestThermoAccuracy(unittest.TestCase):
    """
    Contains tests for accuracy of thermo estimates and symmetry calculations.
    """

    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        self.database = database.thermo

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.Tlist = [300, 400, 500, 600, 800, 1000, 1500]

        self.testCases = [
            # SMILES         symm    H298   S298  Cp300  Cp400  Cp500  Cp600  Cp800 Cp1000 Cp1500

            # 1,3-hexadiene decomposition products
            ['C=CC=CCC',        3,  13.45, 86.37, 29.49, 37.67, 44.54, 50.12, 58.66, 64.95, 74.71],
            ['[CH]=CC=CCC',     3,  72.55, 87.76, 29.30, 36.92, 43.18, 48.20, 55.84, 61.46, 70.18],
            ['C=[C]C=CCC',      3,  61.15, 87.08, 29.68, 36.91, 43.03, 48.11, 55.96, 61.78, 71.54],
            ['C=C[C]=CCC',      3,  61.15, 87.08, 29.68, 36.91, 43.03, 48.11, 55.96, 61.78, 71.54],
            ['C=CC=[C]CC',      3,  70.35, 88.18, 29.15, 36.46, 42.6,  47.6,  55.32, 61.04, 69.95],
            ['C=CC=C[CH]C',     6,  38.24, 84.41, 27.79, 35.46, 41.94, 47.43, 55.74, 61.92, 71.86],
            ['C=CC=CC[CH2]',    2,  62.45, 89.78, 28.72, 36.31, 42.63, 47.72, 55.50, 61.21, 70.05],
            ['[CH3]',           6,  34.81, 46.37,  9.14, 10.18, 10.81, 11.34, 12.57, 13.71, 15.2],
            ['C=CC=C[CH2]',     2,  46.11, 75.82, 22.54, 28.95, 34.24, 38.64, 45.14, 49.97, 57.85],
            ['[CH2]C',          6,  28.6,  59.87, 11.73, 14.47, 17.05, 19.34, 23.02, 25.91, 31.53],
            ['C=CC=[CH]',       1,  85.18, 69.37, 18.93, 23.55, 27.16, 29.92, 34.02, 37.03, 41.81],
            ['C=[CH]',          1,  71.62, 56.61, 10.01, 11.97, 13.66, 15.08, 17.32, 19.05, 21.85],
            ['[CH]=CCC',        3,  58.99, 75.0,  20.38, 25.34, 29.68, 33.36, 39.14, 43.48, 50.22],

            # Cyclic Structures
            ['C1CCCCC1',       12, -29.45, 69.71, 27.20, 37.60, 46.60, 54.80, 67.50, 76.20, 88.50],
            ['C1CCC1',          8,   6.51, 63.35, 17.39, 23.91, 29.86, 34.76, 42.40, 47.98, 56.33],
            ['C1C=CC=C1',       2,  32.5,  65.5,  18.16, 24.71, 30.25, 34.7,  41.25, 45.83, 52.61],
        ]

    @work_in_progress
    def testNewThermoGeneration(self):
        """
        Test that the new ThermoDatabase generates appropriate thermo data.
        """
        for smiles, symm, H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500 in self.testCases:
            Cplist = [Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500]
            species = Species().fromSMILES(smiles)
            species.generate_resonance_structures()
            thermoData = self.database.getThermoDataFromGroups(species)
            molecule = species.molecule[0]
            for mol in species.molecule[1:]:
                thermoData0 = self.database.getAllThermoData(Species(molecule=[mol]))[0][0]
                for data in self.database.getAllThermoData(Species(molecule=[mol]))[1:]:
                    if data[0].getEnthalpy(298) < thermoData0.getEnthalpy(298):
                        thermoData0 = data[0]
                if thermoData0.getEnthalpy(298) < thermoData.getEnthalpy(298):
                    thermoData = thermoData0
                    molecule = mol
            self.assertAlmostEqual(H298, thermoData.getEnthalpy(298) / 4184, places=1,
                                   msg="H298 error for {0}. Expected {1}, but calculated {2}.".format(smiles, H298, thermoData.getEnthalpy(298) / 4184))
            self.assertAlmostEqual(S298, thermoData.getEntropy(298) / 4.184, places=1,
                                   msg="S298 error for {0}. Expected {1}, but calculated {2}.".format(smiles, S298, thermoData.getEntropy(298) / 4.184))
            for T, Cp in zip(self.Tlist, Cplist):
                self.assertAlmostEqual(Cp, thermoData.getHeatCapacity(T) / 4.184, places=1,
                                       msg="Cp{3} error for {0}. Expected {1} but calculated {2}.".format(smiles, Cp, thermoData.getHeatCapacity(T) / 4.184, T))

    @work_in_progress
    def testSymmetryNumberGeneration(self):
        """
        Test we generate symmetry numbers correctly.

        This uses the new thermo database to generate the H298, used
        to select the stablest resonance isomer.
        """
        for smiles, symm, H298, S298, Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500 in self.testCases:
            species = Species().fromSMILES(smiles)
            species.generate_resonance_structures()
            thermoData = self.database.getThermoDataFromGroups(species)
            # pick the molecule with lowest H298
            molecule = species.molecule[0]
            for mol in species.molecule[1:]:
                thermoData0 = self.database.getAllThermoData(Species(molecule=[mol]))[0][0]
                for data in self.database.getAllThermoData(Species(molecule=[mol]))[1:]:
                    if data[0].getEnthalpy(298) < thermoData0.getEnthalpy(298):
                        thermoData0 = data[0]
                if thermoData0.getEnthalpy(298) < thermoData.getEnthalpy(298):
                    thermoData = thermoData0
                    molecule = mol
            self.assertEqual(symm, molecule.calculateSymmetryNumber(),
                             msg="Symmetry number error for {0}. Expected {1} but calculated {2}.".format(smiles, symm, molecule.calculateSymmetryNumber()))


class TestThermoAccuracyAromatics(TestThermoAccuracy):
    """
    Contains tests for accuracy of thermo estimates and symmetry calculations for aromatics only.
    
    A copy of the above class, but with different test compounds.
    """
    def setUp(self):
        self.Tlist = [300, 400, 500, 600, 800, 1000, 1500]
        self.testCases = [
            # SMILES         symm    H298   S298  Cp300  Cp400  Cp500  Cp600  Cp800 Cp1000 Cp1500
            ['c1ccccc1',       12,  19.80, 64.24, 19.44, 26.64, 32.76, 37.80, 45.24, 50.46, 58.38],
            ['c1ccc2ccccc2c1',  4,   36.0, 79.49, 31.94, 42.88, 52.08, 59.62, 70.72, 78.68, 90.24],
        ]

    def __init__(self, *args, **kwargs):
        super(TestThermoAccuracyAromatics, self).__init__(*args, **kwargs)
        self._testMethodDoc = self._testMethodDoc.strip().split('\n')[0] + " for Aromatics.\n"

    def testLongDistanceInteractionInAromaticMolecule(self):
        """
        Test long distance interaction is properly caculated for aromatic molecule.
        """
        spec = Species().fromSMILES('c(O)1c(O)c(C=O)c(C=O)c(O)c(C=O)1')
        spec.generate_resonance_structures()
        thermo = self.database.getThermoDataFromGroups(spec)

        self.assertIn('o_OH_OH', thermo.comment)
        self.assertIn('o_OH_CHO', thermo.comment)
        self.assertIn('o_CHO_CHO', thermo.comment)
        self.assertIn('m_CHO_CHO', thermo.comment)
        self.assertIn('p_OH_OH', thermo.comment)
        self.assertIn('p_OH_CHO', thermo.comment)
        self.assertIn('p_CHO_CHO', thermo.comment)

    def testLongDistanceInteractionInAromaticRadical(self):
        """
        Test long distance interaction is properly caculated for aromatic radical.
        """
        spec = Species().fromSMILES('c([O])1c(C=O)c(C=O)c(OC)cc1')
        spec.generate_resonance_structures()
        thermo = self.database.getThermoDataFromGroups(spec)

        self.assertNotIn('o_OH_CHO', thermo.comment)
        self.assertNotIn('p_OH_MeO', thermo.comment)
        self.assertIn('o_Oj_CHO', thermo.comment)
        self.assertIn('m_Oj_CHO', thermo.comment)
        self.assertIn('p_Oj_OCH3', thermo.comment)
        self.assertIn('o_CHO_CHO', thermo.comment)
        self.assertIn('o_CHO_MeO', thermo.comment)

    def testLongDistanceInteractionInAromaticBiradical(self):
        """
        Test long distance interaction is properly caculated for aromatic biradical.
        """
        spec = Species().fromSMILES('c([O])1c([C]=O)cc(C=O)cc1')
        spec.generate_resonance_structures()
        thermo = self.database.getThermoDataFromGroups(spec)

        thermo = self.database.getThermoDataFromGroups(spec)
        self.assertNotIn('o_OH_CHO', thermo.comment)
        self.assertNotIn('m_CHO_CHO', thermo.comment)
        self.assertNotIn('p_OH_CHO', thermo.comment)
        self.assertNotIn('o_Oj_CHO', thermo.comment)
        self.assertIn('m_Cj=O_CHO', thermo.comment)

class TestCyclicThermo(unittest.TestCase):
    """
    Contains unit tests of the ThermoDatabase class.
    """
    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        self.database = database.thermo

    def testComputeGroupAdditivityThermoForTwoRingMolecule(self):
        """
        The molecule being tested has two rings, one is 13cyclohexadiene5methylene
        the other is benzene ring. This method is to test thermo estimation will
        give two different corrections accordingly. 
        """
        spec = Species().fromSMILES('CCCCCCCCCCCC(CC=C1C=CC=CC1)c1ccccc1')
        spec.generate_resonance_structures()
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
        spec.generate_resonance_structures()
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
        
        smi = 'C12C(C3CCC2C3)C4CCC1C4'#two norbornane rings fused together
        spc = Species().fromSMILES(smi)

        spc.thermo = generateThermoData(spc)

        self.database.getRingGroupsFromComments(spc.thermo)

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
        spe.generate_resonance_structures()
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
        # so here we start with kekulized version and generate_resonance_structures
        # and pick the one with two aromatic rings
        smiles = 'C1=CC2C=CC=C3C=CC(=C1)C=23'
        spe = Species().fromSMILES(smiles)
        spe.generate_resonance_structures()
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

    def testAddPolyRingCorrectionThermoDataFromHeuristicUsingHighlyUnsaturatedPolycyclics1(self):
        """
        Test proper thermo estimation for highly unsaturated polycyclic whose decomposed 
        bicyclics are not stored in database. Those bicyclics thermo will be estimated through
        a heuristic formula.

        In the future, the test assertion may be updated if some of the decomposed bicyclics
        have been added to database.
        """
        # create testing molecule
        smiles = '[CH]=C1C2=C=C3C=CC1C=C32'
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

        self.assertIn('1,4-Cyclohexadiene', ringGroupLabels)
        self.assertIn('Cyclopentene', ringGroupLabels)
        self.assertIn('cyclobutadiene_13', ringGroupLabels)
        self.assertIn('s3_5_6_ane', polycyclicGroupLabels)
        self.assertIn('s2_4_6_ane', polycyclicGroupLabels)
        self.assertIn('s2_4_5_ane', polycyclicGroupLabels)

    def testAddPolyRingCorrectionThermoDataFromHeuristicUsingHighlyUnsaturatedPolycyclics2(self):
        """
        Test proper thermo estimation for highly unsaturated polycyclic whose decomposed 
        bicyclics are not stored in database. Those bicyclics thermo will be estimated through
        a heuristic formula.
        
        In the future, the test assertion may be updated if some of the decomposed bicyclics
        have been added to database.
        """
        # create testing molecule
        smiles = 'C1=C2C#CC3C=CC1C=C23'
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

        self.assertIn('1,4-Cyclohexadiene', ringGroupLabels)
        self.assertIn('Cyclopentyne', ringGroupLabels)
        self.assertIn('Cyclopentadiene', ringGroupLabels)
        self.assertIn('s3_5_6_ane', polycyclicGroupLabels)
        self.assertIn('s2_5_6_ane', polycyclicGroupLabels)
        self.assertIn('s2_5_5_ane', polycyclicGroupLabels)

    def testGetBicyclicCorrectionThermoDataFromHeuristic1(self):
        """
        Test bicyclic correction estimated properly from heuristic formula
        The test molecule "C1=CCC2C1=C2" has a shared atom with Cd atomtype, 
        but in the correction estimation we stil expect the five-member ring 
        part to match Cyclopentene
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
        The test molecule "C1=CCC2=C1C2" doesn't have controversial shared 
        atomtypes in correction estimation, which is regarded as a simple case.
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

    def testGetCopyForOneRing(self):
        """
        This method tests the getCopyForOneRing method, which returns
        an atom object list that contains deep copies of the atoms
        """

        testAtomList=Molecule(SMILES='C1CCCCC1').atoms
        copiedAtomList=getCopyForOneRing(testAtomList)

        testMolecule=Molecule(atoms=testAtomList)
        copiedMolecule=Molecule(atoms=copiedAtomList)

        self.assertTrue(testAtomList!=copiedAtomList)
        self.assertTrue(len(testAtomList)==len(copiedAtomList))
        self.assertTrue(testMolecule.is_equal(copiedMolecule))

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

    def testIsBicyclic1(self):
        """
        Test isBicyclic identifies bicyclic correctly
        The test molecule is bicyclic, we expect isBicyclic()
        returns True.
        """
        smiles = 'C1=CCC2C1=C2'
        mol = Molecule().fromSMILES(smiles)
        polyring = mol.getDisparateRings()[1][0]

        self.assertTrue(isBicyclic(polyring))

    def testIsBicyclic2(self):
        """
        Test isBicyclic identifies bicyclic correctly
        The test molecule is tetracyclic, we expect 
        isBicyclic() returns False
        """
        smiles = 'C1C=C2C=CC=C3C=CC4=CC=CC=1C4=C23'
        mol = Molecule().fromSMILES(smiles)
        polyring = mol.getDisparateRings()[1][0]

        self.assertFalse(isBicyclic(polyring))

    def testFindAromaticBondsFromSubMolecule(self):

        smiles = "C1=CC=C2C=CC=CC2=C1"
        spe = Species().fromSMILES(smiles)
        spe.generate_resonance_structures()
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
        spe.generate_resonance_structures()
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
        for bicyclic in bicyclicList:
            bicyclic.deleteHydrogens()

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
        # so here we start with kekulized version and generate_resonance_structures
        # and pick the one with two aromatic rings
        smiles = 'C1=CC2C=CC=C3C=CC(=C1)C=23'
        spe = Species().fromSMILES(smiles)
        spe.generate_resonance_structures()
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
        for bicyclic in bicyclicList:
            bicyclic.deleteHydrogens()

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
        for bicyclic in bicyclicList:
            bicyclic.deleteHydrogens()

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

    def testCombineCycles(self):
        """
        This method tests the combineCycles method, which simply joins two lists
        together without duplication.
        """
        mainCycle=Molecule(SMILES='C1CCC2CCCCC2C1').atoms
        testCycle1=mainCycle[0:8]
        testCycle2=mainCycle[6:]
        joinedCycle=combineCycles(testCycle1,testCycle2)
        self.assertTrue(sorted(mainCycle)==sorted(joinedCycle))

    def testSplitBicyclicIntoSingleRings1(self):
        """
        Test bicyclic molecule "C1=CCC2C1=C2" can be divided into 
        individual rings properly
        """
        smiles = 'C1=CCC2C1=C2'
        mol = Molecule().fromSMILES(smiles)
        bicyclic = mol.getDisparateRings()[1][0]

        bicyclic_submol = convertRingToSubMolecule(bicyclic)[0]
        single_ring_submols = splitBicyclicIntoSingleRings(bicyclic_submol)
        self.assertEqual(len(single_ring_submols), 2)

        single_ring_submol_a, single_ring_submol_b = sorted(single_ring_submols,
                                key=lambda submol: len(submol.atoms))

        single_ring_submol_a.saturate_unfilled_valence()
        single_ring_submol_b.saturate_unfilled_valence()

        expected_submol_a = Molecule().fromSMILES('C1=CC1')
        expected_submol_a.updateConnectivityValues()

        expected_submol_b = Molecule().fromSMILES('C1=CCCC1')
        expected_submol_b.updateConnectivityValues()


        self.assertTrue(single_ring_submol_a.isIsomorphic(expected_submol_a))
        self.assertTrue(single_ring_submol_b.isIsomorphic(expected_submol_b))

    def testSplitBicyclicIntoSingleRings2(self):
        """
        Test bicyclic molecule "C1=CCC2=C1C2" can be divided into 
        individual rings properly
        """

        smiles = 'C1=CCC2=C1C2'
        mol = Molecule().fromSMILES(smiles)
        bicyclic = mol.getDisparateRings()[1][0]

        bicyclic_submol = convertRingToSubMolecule(bicyclic)[0]
        single_ring_submols = splitBicyclicIntoSingleRings(bicyclic_submol)
        self.assertEqual(len(single_ring_submols), 2)

        single_ring_submol_a, single_ring_submol_b = sorted(single_ring_submols,
                                key=lambda submol: len(submol.atoms))

        single_ring_submol_a.saturate_unfilled_valence()
        single_ring_submol_b.saturate_unfilled_valence()

        expected_submol_a = Molecule().fromSMILES('C1=CC1')
        expected_submol_a.updateConnectivityValues()

        expected_submol_b = Molecule().fromSMILES('C1=CC=CC1')
        expected_submol_b.updateConnectivityValues()

        self.assertTrue(single_ring_submol_a.isIsomorphic(expected_submol_a))
        self.assertTrue(single_ring_submol_b.isIsomorphic(expected_submol_b))

    def testSaturateRingBonds1(self):
        """
        Test unsaturated bonds of "C1=CCC2=C1C2" to be saturated properly
        """
        smiles = 'C1=CCC2=C1C2'
        mol = Molecule().fromSMILES(smiles)
        ring_submol = convertRingToSubMolecule(mol.getDisparateRings()[1][0])[0]

        saturated_ring_submol, alreadySaturated = saturate_ring_bonds(ring_submol)

        expected_saturated_ring_submol = Molecule().fromSMILES('C1CCC2C1C2')
        
        expected_saturated_ring_submol.updateConnectivityValues()

        self.assertFalse(alreadySaturated)
        self.assertEqual(saturated_ring_submol.multiplicity,
                            expected_saturated_ring_submol.multiplicity)
        self.assertTrue(saturated_ring_submol.isIsomorphic(expected_saturated_ring_submol))

    def testSaturateRingBonds2(self):
        """
        Test unsaturated bonds of "C1=CC=C2CCCCC2=C1" to be saturated properly
        """
        smiles = 'C1=CC=C2CCCCC2=C1'
        spe = Species().fromSMILES(smiles)
        spe.generate_resonance_structures()
        mol = spe.molecule[1]
        ring_submol = convertRingToSubMolecule(mol.getDisparateRings()[1][0])[0]

        saturated_ring_submol, alreadySaturated = saturate_ring_bonds(ring_submol)

        expected_spe = Species().fromSMILES('C1=CC=C2CCCCC2=C1')
        expected_spe.generate_resonance_structures()
        expected_saturated_ring_submol = expected_spe.molecule[1]

        expected_saturated_ring_submol.updateConnectivityValues()

        self.assertTrue(alreadySaturated)
        self.assertEqual(saturated_ring_submol.multiplicity,
                            expected_saturated_ring_submol.multiplicity)
        self.assertTrue(saturated_ring_submol.isIsomorphic(expected_saturated_ring_submol))

    def testSaturateRingBonds3(self):
        """
        Test unsaturated bonds of "C1=CC=C2CC=CCC2=C1" to be saturated properly
        """
        smiles = 'C1=CC=C2CC=CCC2=C1'
        spe = Species().fromSMILES(smiles)
        spe.generate_resonance_structures()
        mol = spe.molecule[1]
        ring_submol = convertRingToSubMolecule(mol.getDisparateRings()[1][0])[0]

        saturated_ring_submol, alreadySaturated = saturate_ring_bonds(ring_submol)

        expected_spe = Species().fromSMILES('C1=CC=C2CCCCC2=C1')
        expected_spe.generate_resonance_structures()
        expected_saturated_ring_submol = expected_spe.molecule[1]
        
        expected_saturated_ring_submol.updateConnectivityValues()

        self.assertFalse(alreadySaturated)
        self.assertEqual(saturated_ring_submol.multiplicity,
                            expected_saturated_ring_submol.multiplicity)
        self.assertTrue(saturated_ring_submol.isIsomorphic(expected_saturated_ring_submol))

@attr('auth')
class TestThermoCentralDatabaseInterface(unittest.TestCase):
    """
    Contains unit tests for methods of ThermoCentralDatabaseInterface
    """
    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        self.database = database.thermo

    def connectToTestCentralDatabase(self):

        host, port, username, password = getTestingTCDAuthenticationInfo()
        application = 'test'

        tcdi = ThermoCentralDatabaseInterface(host, port, username, password, application)
        return tcdi

    def testConnectFailure(self):

        host = 'somehost'
        port = 27017
        username = 'me'
        password = 'pswd'
        application = 'test'

        tcdi = ThermoCentralDatabaseInterface(host, port, username, password, application)

        self.assertTrue(tcdi.client is None)

    def testConnectSuccess(self):

        tcdi = self.connectToTestCentralDatabase()

        self.assertTrue(tcdi.client is not None)

    def testSatisfyRegistrationRequirements1(self):
        """
        the species is non-cyclic, currently regarded no need to 
        register in thermo central database
        """
        tcdi = self.connectToTestCentralDatabase()

        species = Species().fromSMILES('C[CH2]')

        thermoData = self.database.getThermoDataFromGroups(species)

        self.assertFalse(tcdi.satisfyRegistrationRequirements(species, thermoData, self.database))

    def testSatisfyRegistrationRequirements2(self):
        """
        the species is for non-cyclic, so no need to register in 
        thermo central database
        """
        tcdi = self.connectToTestCentralDatabase()

        species = Species().fromSMILES('CC')

        thermoData = self.database.getThermoDataFromGroups(species)

        self.assertFalse(tcdi.satisfyRegistrationRequirements(species, thermoData, self.database))


    def testSatisfyRegistrationRequirements3(self):
        """
        the thermo is exact match, so no need to register in 
        thermo central database
        """
        tcdi = self.connectToTestCentralDatabase()

        species = Species().fromSMILES('C1CC1')

        thermoData = self.database.getThermoDataFromGroups(species)

        self.assertFalse(tcdi.satisfyRegistrationRequirements(species, thermoData, self.database))

    def testSatisfyRegistrationRequirements4(self):
        """
        the thermo is from library, so no need to register in 
        thermo central database
        """
        tcdi = self.connectToTestCentralDatabase()

        species = Species().fromSMILES('[H][H]')

        thermoData = self.database.getThermoDataFromLibraries(species)

        self.assertFalse(tcdi.satisfyRegistrationRequirements(species, thermoData, self.database))

    def testSatisfyRegistrationRequirements5(self):
        """
        the thermo is matching generic node, so it needs to register in 
        thermo central database

        In the future, if RMG-database includes corresponding exact match
        this test should be modified.
        """
        tcdi = self.connectToTestCentralDatabase()

        species = Species().fromSMILES('C1C=CC2C=CC2=C1')

        thermoData = self.database.getThermoDataFromGroups(species)

        self.assertTrue(tcdi.satisfyRegistrationRequirements(species, thermoData, self.database))

    def testSatisfyRegistrationRequirements6(self):
        """
        the thermo is matching generic node, so it needs to register in 
        thermo central database

        In the future, if RMG-database includes corresponding exact match
        this test should be modified.
        """
        tcdi = self.connectToTestCentralDatabase()

        species = Species().fromSMILES('C1=C=C2CC23C=CC=1C=C3')

        thermoData = self.database.getThermoDataFromGroups(species)

        self.assertTrue(tcdi.satisfyRegistrationRequirements(species, thermoData, self.database))

    def testRegisterInCentralThermoDB1(self):
        """
        Test situation where both registration_table and results_table have no
        species as the one going to be registered
        """
        # connect to thermo central database
        host, port, username, password = getTestingTCDAuthenticationInfo()
        application = 'test'
        tcdi = ThermoCentralDatabaseInterface(host, port, username, password, application)

        # prepare species to register
        species = Species().fromSMILES('C1=C=C2CC23C=CC=1C=C3')
        expected_aug_inchi = "InChI=1S/C10H6/c1-2-9-7-10(9)5-3-8(1)4-6-10/h3-6H,7H2"

        # select registration table
        # and clean previous data
        db =  getattr(tcdi.client, 'thermoCentralDB')
        registration_table = getattr(db, 'registration_table')
        results_table = getattr(db, 'results_table')
        registration_table.delete_many({"aug_inchi": expected_aug_inchi})
        results_table.delete_many({"aug_inchi": expected_aug_inchi})

        tcdi.registerInCentralThermoDB(species)
        registered_species_entries = list(registration_table.find({"aug_inchi": expected_aug_inchi}))

        # should expect only one registered such species
        self.assertEqual(len(registered_species_entries), 1)
        registered_species_entry = registered_species_entries[0]

        # check all the columns are expected
        registered_species = Species().fromSMILES(str(registered_species_entry['SMILES_input']))
        self.assertEqual(registered_species_entry['aug_inchi'], expected_aug_inchi)
        self.assertTrue(registered_species.isIsomorphic(species))
        self.assertIn(registered_species_entry['status'], ['pending', 'submitted'])

        # clean up the table
        registration_table.delete_many({"aug_inchi": expected_aug_inchi})

    def testRegisterInCentralThermoDB2(self):
        """
        Test situation where registration_table has species as the one going 
        to be registered
        """

        # connect to thermo central database
        host, port, username, password = getTestingTCDAuthenticationInfo()
        application = 'test'

        tcdi = ThermoCentralDatabaseInterface(host, port, username, password, application)

        # prepare species to register
        species = Species().fromSMILES('C1=C=C2CC23C=CC=1C=C3')
        expected_aug_inchi = "InChI=1S/C10H6/c1-2-9-7-10(9)5-3-8(1)4-6-10/h3-6H,7H2"

        # select registration table
        # and clean previous data
        db =  getattr(tcdi.client, 'thermoCentralDB')
        registration_table = getattr(db, 'registration_table')
        results_table = getattr(db, 'results_table')
        registration_table.delete_many({"aug_inchi": expected_aug_inchi})
        registration_table.insert_one({"aug_inchi": expected_aug_inchi})
        results_table.delete_many({"aug_inchi": expected_aug_inchi})

        tcdi.registerInCentralThermoDB(species)
        registered_species_entries = list(registration_table.find({"aug_inchi": expected_aug_inchi}))

        # should expect only one registered such species
        self.assertEqual(len(registered_species_entries), 1)
        registered_species_entry = registered_species_entries[0]

        # check all the columns are expected
        self.assertEqual(registered_species_entry['aug_inchi'], expected_aug_inchi)
        self.assertTrue(len(registered_species_entry), 2)

        # clean up the table
        registration_table.delete_many({"aug_inchi": expected_aug_inchi})

    def testRegisterInCentralThermoDB3(self):
        """
        Test situation where results_table has species as the one going 
        to be registered
        """

        # connect to thermo central database
        host, port, username, password = getTestingTCDAuthenticationInfo()
        application = 'test'
        
        tcdi = ThermoCentralDatabaseInterface(host, port, username, password, application)

        # prepare species to register
        species = Species().fromSMILES('C1=C=C2CC23C=CC=1C=C3')
        expected_aug_inchi = "InChI=1S/C10H6/c1-2-9-7-10(9)5-3-8(1)4-6-10/h3-6H,7H2"

        # select registration table
        # and clean previous data
        db =  getattr(tcdi.client, 'thermoCentralDB')
        registration_table = getattr(db, 'registration_table')
        results_table = getattr(db, 'results_table')
        registration_table.delete_many({"aug_inchi": expected_aug_inchi})
        results_table.delete_many({"aug_inchi": expected_aug_inchi})
        results_table.insert_one({"aug_inchi": expected_aug_inchi})

        tcdi.registerInCentralThermoDB(species)
        registered_species_entries = list(registration_table.find({"aug_inchi": expected_aug_inchi}))

        # should expect only one registered such species
        self.assertEqual(len(registered_species_entries), 0)

        # clean up the table
        results_table.delete_many({"aug_inchi": expected_aug_inchi})

def getTestingTCDAuthenticationInfo():

    try:
        host = os.environ['TCD_HOST']
        port = int(os.environ['TCD_PORT'])
        username = os.environ['TCD_USER']
        password = os.environ['TCD_PW']
    except KeyError:
        print('Thermo Central Database Authentication Environment Variables Not Completely Set!')
        return 'None', 0, 'None', 'None'

    return host, port, username, password

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))

