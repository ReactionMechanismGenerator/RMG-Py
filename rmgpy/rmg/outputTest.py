#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import os
import unittest
import shutil 

from model import CoreEdgeReactionModel, ReactionModel 
from rmgpy.chemkin import loadChemkinFile

from output import *

###################################################

class TestOutput(unittest.TestCase):

	def testSaveOutputHTML(self):
		"""
		This example is to test if an HTML file can be generated
		for the provided chemkin model.
		"""
		folder = os.path.join(os.getcwd(),'rmgpy/rmg/test_data/saveOutputHTML/')
		
		chemkinPath = os.path.join(folder, 'eg6', 'chem_annotated.inp')
		dictionaryPath = os.path.join(folder,'eg6', 'species_dictionary.txt')

		# loadChemkinFile
		species, reactions = loadChemkinFile(chemkinPath, dictionaryPath) 

		# convert it into a reaction model:
		core = ReactionModel(species, reactions)
		cerm = CoreEdgeReactionModel(core)

		out = os.path.join(folder, 'output.html')
		saveOutputHTML(out, cerm)

		self.assertTrue(os.path.isfile(out))
		os.remove(out)
		shutil.rmtree(os.path.join(folder,'species'))
