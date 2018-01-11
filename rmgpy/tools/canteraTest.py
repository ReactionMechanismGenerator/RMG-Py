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
import os
import numpy
from rmgpy.tools.canteraModel import findIgnitionDelay, CanteraCondition, Cantera
from rmgpy.quantity import Quantity
import rmgpy
class CanteraTest(unittest.TestCase):

    def testIgnitionDelay(self):
        """
        Test that findIgnitionDelay() works.
        """

        t = numpy.arange(0,5,0.5)
        P = numpy.array([0,0.33,0.5,0.9,2,4,15,16,16.1,16.2])
        OH = numpy.array([0,0.33,0.5,0.9,2,4,15,16,7,2])
        CO = OH*0.9

        t_ign = findIgnitionDelay(t,P)
        self.assertEqual(t_ign,2.75)

        t_ign = findIgnitionDelay(t,OH,'maxHalfConcentration')
        self.assertEqual(t_ign,3)

        t_ign = findIgnitionDelay(t,[OH,CO], 'maxSpeciesConcentrations')
        self.assertEqual(t_ign,3.5)

    def testRepr(self):
        """
        Test that the repr function for a CanteraCondition object can reconstitute
        the same object
        """
        reactorType='IdealGasReactor'
        molFrac={'CC': 0.05, '[Ar]': 0.95}
        P=(3,'atm')
        T=(1500,'K')
        terminationTime=(5e-5,'s')
        condition = CanteraCondition(reactorType, 
                        terminationTime,
                        molFrac, 
                        T0=T,
                        P0=P)
        reprCondition=eval(condition.__repr__())
        self.assertEqual(reprCondition.T0.value_si,Quantity(T).value_si)
        self.assertEqual(reprCondition.P0.value_si,Quantity(P).value_si)
        self.assertEqual(reprCondition.V0,None)
        self.assertEqual(reprCondition.molFrac,molFrac)
        

class RMGToCanteraTest(unittest.TestCase):
    """
    Contains unit tests for the conversion of RMG species and reaction objects to Cantera objects.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        from rmgpy.chemkin import loadChemkinFile
        folder = os.path.join(os.path.dirname(rmgpy.__file__),'tools/data/various_kinetics')
        
        chemkinPath = os.path.join(folder, 'chem_annotated.inp')
        dictionaryPath = os.path.join(folder, 'species_dictionary.txt')
        transportPath = os.path.join(folder, 'tran.dat')
        
        species, reactions = loadChemkinFile(chemkinPath, dictionaryPath,transportPath) 
        
        self.rmg_ctSpecies = [spec.toCantera(useChemkinIdentifier = True) for spec in species]
        self.rmg_ctReactions = []
        for rxn in reactions:
            convertedReactions = rxn.toCantera(species, useChemkinIdentifier = True)
            if isinstance(convertedReactions,list):
                self.rmg_ctReactions.extend(convertedReactions)
            else:
                self.rmg_ctReactions.append(convertedReactions)
        job = Cantera()
        job.loadChemkinModel(chemkinPath, transportFile=transportPath,quiet=True)
        self.ctSpecies = job.model.species()
        self.ctReactions = job.model.reactions()
    
    def testSpeciesConversion(self):
        """
        Test that species objects convert properly
        """
        from rmgpy.tools.canteraModel import checkEquivalentCanteraSpecies
        for i in range(len(self.ctSpecies)):
            self.assertTrue(checkEquivalentCanteraSpecies(self.ctSpecies[i],self.rmg_ctSpecies[i]))
            
            
    def testReactionConversion(self):
        """
        Test that species objects convert properly
        """
        from rmgpy.tools.canteraModel import checkEquivalentCanteraReaction
        for i in range(len(self.ctReactions)):
            self.assertTrue(checkEquivalentCanteraReaction(self.ctReactions[i],self.rmg_ctReactions[i]))
        
