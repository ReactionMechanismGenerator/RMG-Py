import unittest
import os
import numpy
from rmgpy.tools.canteraModel import *

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
