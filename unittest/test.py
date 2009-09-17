#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

from atomtest import *
from bondtest import *
from datatest import *
from structuretest import *
from thermotest import *
from speciestest import *
from reactiontest import *
from simulationtest import *
from canteraloadertest import *

################################################################################

if __name__ == '__main__':

	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )

