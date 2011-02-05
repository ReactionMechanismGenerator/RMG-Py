#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

from atomtest import *
from bondtest import *
from datatest import *
from graphtest import *
from structuretest import *
from thermotest import *
from speciestest import *
from reactiontest import *
## Now that we use cantera, this no longer works 
## (it used to compare our simulation solver with a cantera one)
# from simulationtest import *
## Also, this is broken, since now you can't make a species without a structure,
## and a chemkin or cantera mechanism does not include species structures
# from canteraloadertest import *

################################################################################

if __name__ == '__main__':

	unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )

