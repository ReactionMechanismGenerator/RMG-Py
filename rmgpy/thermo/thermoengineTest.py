#!/usr/bin/env python
# encoding: utf-8 -*-

"""
This module contains unit tests of the rmgpy.parallel module.
"""

import os
import sys
import unittest
import random
from external.wip import work_in_progress

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.rmg.main import RMG
from rmgpy.scoop_framework.framework import TestScoopCommon

from rmgpy.species import Species
from rmgpy.thermo.thermoengine import submit, satisfyRegistrationRequirements

try:
    from scoop import futures, _control, shared
except ImportError, e:
    import logging as logging
    logging.debug("Could not properly import SCOOP.")

def setUpModule():
    """A function that is run ONCE before all unit tests in this module."""
    global database
    database = RMGDatabase()
    database.loadThermo(os.path.join(settings['database.directory'], 'thermo'))

def tearDownModule():
    """A function that is run ONCE after all unit tests in this module."""
    from rmgpy.data import rmg
    rmg.database = None

def load():
    tearDown()
    rmg = RMG()#for solvent
    database = RMGDatabase()
    database.loadThermo(os.path.join(settings['database.directory'], 'thermo'))
    database.loadTransport(os.path.join(settings['database.directory'], 'transport'))
    database.loadSolvation(os.path.join(settings['database.directory'], 'solvation'))

def tearDown():
    """
    Reset the loaded database
    """
    import rmgpy.data.rmg
    rmgpy.data.rmg.database = None

def funcSubmit():
    """
    Test that we can submit a number of species.
    """
    load()

    spcs = [
            Species().fromSMILES('C'),\
            Species().fromSMILES('CC'), \
            Species().fromSMILES('CCC')
            ]
    
    for spc in spcs:
        submit(spc)

    return True

def funcGet():
    """
    Test if we can retrieve thermo of species even before we have submitted them explicitly.
    """
    load()

    spcs = [
            Species().fromSMILES('C'),
            Species().fromSMILES('CC'), \
            Species().fromSMILES('CCC')
            ]
    
    output = []
    for spc in spcs:
        data = spc.getThermoData()
        output.append((spc, data))

    for spc, data in output:
        if not data:
            return False

    return True

def funcSubmitGet():
    """
    Test if we can retrieve thermo of species after submitting some of them.
    """
    load()

    spcs = [
            Species().fromSMILES('C'),\
            Species().fromSMILES('CC'), \
            Species().fromSMILES('CCC')
            ]
    
    for spc in spcs:
        submit(spc)

    absent = Species().fromSMILES('[CH3]')
    data = absent.getThermoData()
    if not data: return False

    present = Species().fromSMILES('CC')
    data = present.getThermoData()
    if not data: return False

    random.shuffle(spcs)
    for spc in spcs:
        data = spc.getThermoData()
        if not data: return False        

    return True

@work_in_progress
class AsyncThermoTest(TestScoopCommon):

    def __init__(self, *args, **kwargs):
        # Parent initialization
        super(self.__class__, self).__init__(*args, **kwargs)
        
        # Only setup the scoop framework once, and not in every test method:
        super(self.__class__, self).setUp()

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "test currently only runs on linux")
    def testSubmit(self):
        """
        Test that we can submit a request to generate
        thermo/transport for a number of species.
        """
        result = futures._startup(funcSubmit)
        self.assertEquals(result, True)

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "test currently only runs on linux")
    def testGet(self):
        """
        Test that we can get the data of a number of species.
        """
        result = futures._startup(funcGet)
        self.assertEquals(result, True)

class TestThermoCentralDatabase(unittest.TestCase):
    """
    Contains unit tests for methods related to thermo central database
    """
    @classmethod
    def setUpClass(self):
        """A function that is run ONCE before all unit tests in this class."""
        global database
        self.database = database.thermo

    def testSatisfyRegistrationRequirements1(self):
        """
        the species is radical, currently not allowed to register 
        in thermo central database
        """

        species = Species().fromSMILES('C[CH2]')

        thermoData = self.database.getThermoDataFromGroups(species)

        self.assertFalse(satisfyRegistrationRequirements(species, thermoData, self.database))

    def testSatisfyRegistrationRequirements2(self):
        """
        the species is for non-cyclic, so no need to register in 
        thermo central database
        """

        species = Species().fromSMILES('CC')

        thermoData = self.database.getThermoDataFromGroups(species)

        self.assertFalse(satisfyRegistrationRequirements(species, thermoData, self.database))


    def testSatisfyRegistrationRequirements3(self):
        """
        the thermo is exact match, so no need to register in 
        thermo central database
        """

        species = Species().fromSMILES('C1CC1')

        thermoData = self.database.getThermoDataFromGroups(species)

        self.assertFalse(satisfyRegistrationRequirements(species, thermoData, self.database))

    def testSatisfyRegistrationRequirements4(self):
        """
        the thermo is from library, so no need to register in 
        thermo central database
        """

        species = Species().fromSMILES('[H][H]')

        thermoData = self.database.getThermoDataFromLibraries(species)

        self.assertFalse(satisfyRegistrationRequirements(species, thermoData, self.database))

    def testSatisfyRegistrationRequirements5(self):
        """
        the thermo is matching generic node, so it needs to register in 
        thermo central database

        In the future, if RMG-database includes corresponding exact match
        this test should be modified.
        """

        species = Species().fromSMILES('C1C=CC2C=CC2=C1')

        thermoData = self.database.getThermoDataFromGroups(species)

        self.assertTrue(satisfyRegistrationRequirements(species, thermoData, self.database))

    def testSatisfyRegistrationRequirements6(self):
        """
        the thermo is matching generic node, so it needs to register in 
        thermo central database

        In the future, if RMG-database includes corresponding exact match
        this test should be modified.
        """

        species = Species().fromSMILES('C1=C=C2CC23C=CC=1C=C3')

        thermoData = self.database.getThermoDataFromGroups(species)

        self.assertTrue(satisfyRegistrationRequirements(species, thermoData, self.database))

if __name__ == '__main__' and os.environ.get('IS_ORIGIN', "1") == "1":
    unittest.main()
