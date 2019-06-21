#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
This script contains unit tests of the :mod:`rmgobject` module.
"""

import unittest

import numpy as np

from rmgpy.rmgobject import RMGObject

################################################################################

class PseudoRMGObject(RMGObject):
    """
    Child class of RMG Object with more attributes for testing.
    """

    def __init__(self, a=None, b=None, c=None, d=None):
        self.a = a
        self.b = b
        self.c = c
        self.d = d


class TestRMGObject(unittest.TestCase):
    """
    Contains unit tests for the RMGObject class
    """
    
    def test_save_int(self):
        """Test saving ints"""
        obj = PseudoRMGObject(a=1, b=5)
        result = obj.as_dict()

        expected = {'class': 'PseudoRMGObject', 'a': 1, 'b': 5}

        self.assertEqual(result, expected)

    def test_read_int(self):
        """Test reading ints"""
        data = {'a': 1, 'b': 5}
        obj = PseudoRMGObject()
        obj.make_object(data, class_dict={'PseudoRMGObject': PseudoRMGObject})

        self.assertIsInstance(obj.a, int)
        self.assertEqual(obj.a, 1)
        self.assertIsInstance(obj.b, int)
        self.assertEqual(obj.b, 5)

    def test_save_float(self):
        """Test saving floats"""
        obj = PseudoRMGObject(a=1.0, b=5.0)
        result = obj.as_dict()

        expected = {'class': 'PseudoRMGObject', 'a': 1.0, 'b': 5.0}

        self.assertEqual(result, expected)

    def test_read_float(self):
        """Test reading floats"""
        data = {'a': 1.0, 'b': 5.0}
        obj = PseudoRMGObject()
        obj.make_object(data, class_dict={'PseudoRMGObject': PseudoRMGObject})

        self.assertIsInstance(obj.a, float)
        self.assertEqual(obj.a, 1.0)
        self.assertIsInstance(obj.a, float)
        self.assertEqual(obj.b, 5.0)

    def test_save_str(self):
        """Test saving strings"""
        obj = PseudoRMGObject(a='foo', b='bar')
        result = obj.as_dict()

        expected = {'class': 'PseudoRMGObject', 'a': 'foo', 'b': 'bar'}

        self.assertEqual(result, expected)

    def test_read_str(self):
        """Test reading strings"""
        data = {'a': 'foo', 'b': 'bar'}
        obj = PseudoRMGObject()
        obj.make_object(data, class_dict={'PseudoRMGObject': PseudoRMGObject})

        self.assertEqual(obj.a, 'foo')
        self.assertEqual(obj.b, 'bar')

    def test_save_empty_str(self):
        """Test saving empty strings"""
        obj = PseudoRMGObject(a='', b='bar')
        result = obj.as_dict()

        expected = {'class': 'PseudoRMGObject', 'b': 'bar'}

        self.assertEqual(result, expected)

    def test_read_empty_str(self):
        """Test reading empty strings"""
        data = {'a': '', 'b': 'bar'}
        obj = PseudoRMGObject()
        obj.make_object(data, class_dict={'PseudoRMGObject': PseudoRMGObject})

        self.assertEqual(obj.a, '')
        self.assertEqual(obj.b, 'bar')

    def test_save_dict(self):
        """Test saving dictionaries"""
        obj = PseudoRMGObject(a={'foo': 'bar'})
        result = obj.as_dict()

        expected = {'class': 'PseudoRMGObject', 'a': {'foo': 'bar'}}

        self.assertEqual(result, expected)

    def test_read_dict(self):
        """Test reading dictionaries"""
        data = {'a': {'foo': 'bar'}}
        obj = PseudoRMGObject()
        obj.make_object(data, class_dict={'PseudoRMGObject': PseudoRMGObject})

        self.assertEqual(obj.a, {'foo': 'bar'})

    def test_save_mix(self):
        """Test saving mix of builtin types"""
        obj = PseudoRMGObject(a=1, b=5.0, c='foobar', d={'foo': 'bar'})
        result = obj.as_dict()

        expected = {'class': 'PseudoRMGObject', 'a': 1, 'b': 5.0, 'c': 'foobar', 'd': {'foo': 'bar'}}

        self.assertEqual(result, expected)

    def test_read_mix(self):
        """Test reading mix of builtin types"""
        data = {'a': 1, 'b': 5.0, 'c': 'foobar', 'd': {'foo': 'bar'}}
        obj = PseudoRMGObject()
        obj.make_object(data, class_dict={'PseudoRMGObject': PseudoRMGObject})

        self.assertIsInstance(obj.a, int)
        self.assertEqual(obj.a, 1)
        self.assertIsInstance(obj.b, float)
        self.assertEqual(obj.b, 5.0)
        self.assertEqual(obj.c, 'foobar')
        self.assertEqual(obj.d, {'foo': 'bar'})

    def test_save_numpy(self):
        """Test saving numpy array"""
        obj = PseudoRMGObject(a=np.array([1.0, 2.0, 3.0]))
        result = obj.as_dict()

        expected = {'class': 'PseudoRMGObject', 'a': [1.0, 2.0, 3.0]}

        self.assertEqual(result, expected)

    def test_save_object(self):
        """Test saving another object"""
        obj = PseudoRMGObject(a=PseudoRMGObject(b='foobar'))
        result = obj.as_dict()

        expected = {'class': 'PseudoRMGObject', 'a': {'class': 'PseudoRMGObject', 'b': 'foobar'}}

        self.assertEqual(result, expected)

    def test_read_object(self):
        """Test reading an object"""
        data = {'a': {'class': 'PseudoRMGObject', 'b': 'foobar'}}
        obj = PseudoRMGObject()
        obj.make_object(data, class_dict={'PseudoRMGObject': PseudoRMGObject})

        self.assertIsInstance(obj.a, PseudoRMGObject)
        self.assertEqual(obj.a.b, 'foobar')

    def test_save_object_list(self):
        """Test saving a list of objects"""
        obj = PseudoRMGObject(a=[PseudoRMGObject(b='foobar'), PseudoRMGObject(c=5.0)])
        result = obj.as_dict()

        expected = {'class': 'PseudoRMGObject',
                    'a': [{'class': 'PseudoRMGObject', 'b': 'foobar'},
                          {'class': 'PseudoRMGObject', 'c': 5.0}]}

        self.assertEqual(result, expected)

    def test_read_object_list(self):
        """Test reading a list of objects"""
        data = {'a': [{'class': 'PseudoRMGObject', 'b': 'foobar'},
                      {'class': 'PseudoRMGObject', 'c': 5.0}]}
        obj = PseudoRMGObject()
        obj.make_object(data, class_dict={'PseudoRMGObject': PseudoRMGObject})

        self.assertIsInstance(obj.a, list)
        self.assertEqual(len(obj.a), 2)
        self.assertIsInstance(obj.a[0], PseudoRMGObject)
        self.assertEqual(obj.a[0].b, 'foobar')
        self.assertIsInstance(obj.a[1], PseudoRMGObject)
        self.assertEqual(obj.a[1].c, 5.0)

    def test_save_empty_list(self):
        """Test saving an empty list"""
        obj = PseudoRMGObject(a=[PseudoRMGObject(b='foobar'), PseudoRMGObject(c=5.0)], b=[])
        result = obj.as_dict()

        expected = {'class': 'PseudoRMGObject',
                    'a': [{'class': 'PseudoRMGObject', 'b': 'foobar'},
                          {'class': 'PseudoRMGObject', 'c': 5.0}],
                    'b': []}

        self.assertEqual(result, expected)

    def test_read_empty_list(self):
        """Test reading an smpty list"""
        data = {'a': [{'class': 'PseudoRMGObject', 'b': 'foobar'},
                      {'class': 'PseudoRMGObject', 'c': 5.0}],
                'b': []}
        obj = PseudoRMGObject()
        obj.make_object(data, class_dict={'PseudoRMGObject': PseudoRMGObject})

        self.assertIsInstance(obj.a, list)
        self.assertEqual(len(obj.a), 2)
        self.assertIsInstance(obj.a[0], PseudoRMGObject)
        self.assertEqual(obj.a[0].b, 'foobar')
        self.assertIsInstance(obj.a[1], PseudoRMGObject)
        self.assertEqual(obj.a[1].c, 5.0)
        self.assertEqual(obj.b, [])

    def test_save_1_nested_dict(self):
        """Test saving dictionaries"""
        obj = PseudoRMGObject(a={'foo': {'bar': 'chocolate', 'chang': 'burger'}})
        result = obj.as_dict()

        expected = {'class': 'PseudoRMGObject', 'a': {'foo': {'bar': 'chocolate', 'chang': 'burger'}}}

        self.assertEqual(result, expected)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
