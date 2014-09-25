#!/usr/bin/env python
# encoding: utf-8 -*-

import numpy
import unittest
import os

from rmgpy.cantherm.molepro import MoleProLog
import rmgpy.constants as constants

################################################################################

class MoleProTest(unittest.TestCase):
    """
    Contains unit tests for the chempy.io.gaussian module, used for reading
    and writing Molepro files.
    """
    
    def testLoadDzFromMoleProLog_F12(self):
        """
        Uses a Molepro log file for ethylene_dz (C2H4) to test that F12a
        energy can be properly read.
        """
        
        log=MoleProLog(os.path.join(os.path.dirname(__file__),'files','ethylene_f12_dz.out'))
        E0=log.loadCCSDEnergy()
        
        self.assertAlmostEqual(E0 / constants.Na / constants.E_h, -78.474353559604, 5)
    
    def testLoadQzFromMoleProLog_F12(self):
        """
        Uses a Molepro log file for ethylene_qz (C2H4) to test that F12b
        energy can be properly read.
        """
        
        log=MoleProLog(os.path.join(os.path.dirname(__file__),'files','ethylene_f12_qz.out'))
        E0=log.loadCCSDEnergy()
        
        self.assertAlmostEqual(E0 / constants.Na / constants.E_h, -78.472682755635, 5)

    def testLoadRadFromMoleProLog_F12(self):
        """
        Uses a Molepro log file for OH (C2H4) to test that radical
        energy can be properly read.
        """
        
        log=MoleProLog(os.path.join(os.path.dirname(__file__),'files','OH_f12.out'))
        E0=log.loadCCSDEnergy()
        
        self.assertAlmostEqual(E0 / constants.Na / constants.E_h, -75.663696424380, 5)