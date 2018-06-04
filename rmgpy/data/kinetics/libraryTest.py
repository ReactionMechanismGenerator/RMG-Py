################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
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

import unittest
import os.path
import shutil

from rmgpy import settings
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.kinetics.model import PDepKineticsModel
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
###################################################

class TestLibrary(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        A function run ONCE before all unit tests in this class.
        """
        # Set up a dummy database
        cls.database = KineticsDatabase()
        cls.database.loadLibraries(os.path.join(settings['test_data.directory'], 'testing_database','kinetics','libraries'),
                                  libraries=None) #this loads all of them: ['GRI-Mech3.0', 'ethane-oxidation'])
        cls.libraries = cls.database.libraries

    def testGetLibraryReactions(self):
        """
        test that getLibraryReactions loads reactions correctly 
        """
        libRxns = self.libraries['GRI-Mech3.0'].getLibraryReactions()
        for rxn in libRxns:
            self.assertIsInstance(rxn,LibraryReaction)
        
        libRxns = self.libraries['ethane-oxidation'].getLibraryReactions() #should have no direct library reactions
        for rxn in libRxns:
            if isinstance(rxn.kinetics,PDepKineticsModel):
                self.assertIsInstance(rxn,LibraryReaction) #can load pdep as networks yet so load as libraries
            else:
                self.assertIsInstance(rxn,TemplateReaction) #all reactions are template based
                
    def testSaveLibrary(self):
        """
        This tests the the library.save method by writing a new temporary file and
        loading it and comparing the original and copied reactions
        """
        os.makedirs(os.path.join(settings['test_data.directory'], 'testing_database','kinetics','libraries','eth-oxcopy'))
        try:
            self.libraries['ethane-oxidation'].save(os.path.join(settings['test_data.directory'], 'testing_database','kinetics','libraries','eth-oxcopy','reactions.py'))
            self.database.loadLibraries(os.path.join(settings['test_data.directory'], 'testing_database','kinetics','libraries'),
                                      libraries=None) #this loads all of them: ['GRI-Mech3.0', 'ethane-oxidation','eth-oxcopy'])
            oriRxns = self.database.libraries['ethane-oxidation'].getLibraryReactions()
            copyRxns = self.database.libraries['eth-oxcopy'].getLibraryReactions()
            self.assertTrue(all([repr(oriRxns[i])==repr(copyRxns[i]) for i in xrange(len(oriRxns))]))
        finally:
            shutil.rmtree(os.path.join(settings['test_data.directory'], 'testing_database','kinetics','libraries','eth-oxcopy'))