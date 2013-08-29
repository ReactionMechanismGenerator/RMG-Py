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
    
    def testLoadEthyleneFromMoleProLog_F12(self):
        """
        Uses a Molepro log file for ethylene (C2H4) to test that its
        energy can be properly read.
        """
        
        log=MoleProLog(os.path.join(os.path.dirname(__file__),'test','ethylene_f12.out'))
        E0=log.loadCCSDEnergy()
        
        self.assertAlmostEqual(E0 / constants.Na / constants.E_h, -78.474353559604, 5)