#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

from gaussianTest import *
from geometryTest import *
from graphTest import *
from moleculeTest import *
from reactionTest import *
from statesTest import *
from thermoTest import *

unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
