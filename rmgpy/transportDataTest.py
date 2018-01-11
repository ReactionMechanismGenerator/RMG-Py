#!/usr/bin/env python
# encoding: utf-8

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

"""
This script contains unit test of the :mod: 'rmgpy.transport' module and :mod: 'rmgpy.data.transport' module
"""

import unittest

import rmgpy.constants as constants
from rmgpy.quantity import DipoleMoment, Length, Volume, Energy
from rmgpy.transport import TransportData

#################################################################################

class TestTransportData(unittest.TestCase):
    """
    Contains unit test of the :class: 'transportData' class
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.shapeIndex = 1
        self.epsilon = Energy(2.104, 'kJ/mol')
        self.sigma = Length(3.402, 'angstroms')
        self.dipoleMoment = DipoleMoment(1.000, 'C*m')
        self.polarizability = Volume(0.134, 'angstroms^3')
        self.rotrelaxcollnum = 0.000
        self.comment = 'test'

        self.transport = TransportData(
            shapeIndex=self.shapeIndex,
            epsilon=self.epsilon,
            sigma=self.sigma,
            dipoleMoment=self.dipoleMoment,
            polarizability=self.polarizability,
            rotrelaxcollnum=self.rotrelaxcollnum,
            comment=self.comment,
        )

    def test_shapeIndex(self):
        """
        Test that the TransportData shapeIndex property was properly set.
        """
        self.assertAlmostEqual(self.transport.shapeIndex, self.shapeIndex, 6)

    def test_epsilon(self):
        """
        Test that the TransportData epsilon property was properly set.
        """
        self.assertAlmostEqual(self.transport.epsilon.value_si, self.epsilon.value_si, 6)

    def test_sigma(self):
        """
        Test that the TransportData sigma property was properly set.
        """
        self.assertAlmostEqual(self.transport.sigma.value_si * 1e10, self.sigma.value_si * 1e10, 6)

    def test_dipoleMoment(self):
        """
        Test that the TransportData dipoleMoment property was properly set.
        """
        self.assertAlmostEqual(self.transport.dipoleMoment.value_si, self.dipoleMoment.value_si, 6)

    def test_polarizability(self):
        """
        Test that the TransportData polarizability property was properly set.
        """
        self.assertAlmostEqual(self.transport.polarizability.value_si, self.polarizability.value_si, 6)

    def test_rotrelaxcollnum(self):
        """
        Test that the TransportData rotrelaxcollnum property was properly set.
        """
        self.assertAlmostEqual(self.transport.rotrelaxcollnum, self.rotrelaxcollnum, 6)

    def test_comment(self):
        """
        Test that the TransportData comment property was properly set.
        """
        self.assertEqual(self.transport.comment, self.comment)

    def test_getCollisionFrequency(self):
        """
        Test the LennardJones.getCollisionFrequency() method.
        """
        T = 1000; P = 1.0e5
        M = P / constants.R / T
        mu = 1.0
        omega = self.transport.getCollisionFrequency(T, M, mu)
        self.assertAlmostEqual(omega / 1.17737e10, 1.0, 4)
 
    def test_pickle(self):
        """
        Test that a TransportData object can be pickled and unpickled with no loss of information.
        """
        import cPickle
        transport = cPickle.loads(cPickle.dumps(self.transport,-1))
        self.assertAlmostEqual(self.transport.shapeIndex, transport.shapeIndex, 4)
        self.assertAlmostEqual(self.transport.epsilon.value_si, transport.epsilon.value_si, 4)
        self.assertAlmostEqual(self.transport.sigma.value_si, transport.sigma.value_si, 4)
        self.assertAlmostEqual(self.transport.dipoleMoment.value_si, transport.dipoleMoment.value_si, 4)
        self.assertAlmostEqual(self.transport.polarizability.value_si, transport.polarizability.value_si, 4)
        self.assertAlmostEqual(self.transport.rotrelaxcollnum, transport.rotrelaxcollnum, 4)
        self.assertEqual(self.transport.comment, transport.comment)

    def test_repr(self):
        """
        Test that a TransportData object can be reconstructed from its repr() output with no loss of information
        """
        transport = None
        exec('transport = {0!r}'.format(self.transport))
        self.assertAlmostEqual(self.transport.shapeIndex, transport.shapeIndex, 4)
        self.assertAlmostEqual(self.transport.epsilon.value_si, transport.epsilon.value_si, 4)
        self.assertAlmostEqual(self.transport.sigma.value_si, transport.sigma.value_si, 4)
        self.assertAlmostEqual(self.transport.dipoleMoment.value_si, transport.dipoleMoment.value_si, 4)
        self.assertAlmostEqual(self.transport.polarizability.value_si, transport.polarizability.value_si, 4)
        self.assertAlmostEqual(self.transport.rotrelaxcollnum, transport.rotrelaxcollnum, 4)
        self.assertEqual(self.transport.comment, transport.comment)
        
    def test_toCantera(self):
        """
        Test that the Cantera GasTransportData creation is successful.
        """
        transport = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstrom'), dipoleMoment=(2,'De'), polarizability=(1,'angstrom^3'), rotrelaxcollnum=15.0, comment="""GRI-Mech""")
        rmg_ctTransport = transport.toCantera()
        import cantera as ct
        ctSpecies = ct.Species.fromCti("""species(name=u'Ar',
        atoms='Ar:1',
        transport=gas_transport(geom='atom',
                                diam=3.33,
                                well_depth=136.501,
                                dipole=2.0,
                                polar=1.0,
                                rot_relax=15.0))""")
        
        ctTransport = ctSpecies.transport
        
        self.assertAlmostEqual(rmg_ctTransport.geometry, ctTransport.geometry)
        self.assertAlmostEqual(rmg_ctTransport.acentric_factor, ctTransport.acentric_factor)
        self.assertAlmostEqual(rmg_ctTransport.diameter, ctTransport.diameter)
        self.assertAlmostEqual(rmg_ctTransport.dipole, ctTransport.dipole)
        self.assertAlmostEqual(rmg_ctTransport.polarizability, ctTransport.polarizability)
        self.assertAlmostEqual(rmg_ctTransport.rotational_relaxation, ctTransport.rotational_relaxation)
        self.assertAlmostEqual(rmg_ctTransport.well_depth, ctTransport.well_depth)


#################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
