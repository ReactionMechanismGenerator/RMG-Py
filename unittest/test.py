#!/usr/bin/env python
# encoding: utf-8 -*-

"""
This module can be used to invoke all of the unit tests in a single command.
"""

import unittest

from cantherm.gaussianTest import *
from cantherm.geometryTest import *

from molecule.atomtypeTest import *
from molecule.elementTest import *
from molecule.graphTest import *
from molecule.groupTest import *
from molecule.moleculeTest import *

from kineticsTest import *
from quantityTest import *
from reactionTest import *
from speciesTest import *
from statmechTest import *
from thermoTest import *

unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
