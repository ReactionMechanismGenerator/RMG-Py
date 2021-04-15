#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
from dataclasses import dataclass

import numpy as np

from rmgpy.quantity import ScalarQuantity, ArrayQuantity
from rmgpy.rmgobject import RMGObject, expand_to_dict, recursive_make_object


################################################################################


class PseudoRMGObject(RMGObject):
    """
    Child class of RMG Object with more attributes for testing.
    """

    def __init__(self, a=None, b=None, c=None, d=None):
        super(PseudoRMGObject, self).__init__()
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

        expected = {'class': 'PseudoRMGObject', 'a': {'object': [1.0, 2.0, 3.0], 'class': 'np_array'}}

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


class TestExpandAndMakeFromDictionaries(unittest.TestCase):
    """
    Contains unit tests for the expand_to_dict and recursive_make_object function utilized by RMGObjects
    """

    def setUp(self):
        self.np_array = np.array([[1, 2], [3, 4]])
        self.np_dict = {'class': 'np_array', 'object': [[1, 2], [3, 4]]}

        self.array_quantity = ArrayQuantity(value=self.np_array, units='kJ/mol')
        self.array_dict = {'class': 'ArrayQuantity', 'value': self.np_dict, 'units': 'kJ/mol'}

        self.scalar_quantity = ScalarQuantity(value=500.0, units='K')
        self.scalar_dict = {'class': 'ScalarQuantity', 'value': 500.0, 'units': 'K'}

        self.highly_nested_object = PseudoRMGObject(a=PseudoRMGObject(a=PseudoRMGObject(b=self.np_array,
                                                                                        c=PseudoRMGObject(
                                                                                            c=self.array_quantity,
                                                                                            d=PseudoRMGObject(
                                                                                                a=self.scalar_quantity,
                                                                                                b=PseudoRMGObject()
                                                                                                )
                                                                                            )
                                                                                        )
                                                                      ),
                                                    b=6
                                                    )
        self.highly_nest_dictionary = {'class': 'PseudoRMGObject',
                                       'a': {'class': 'PseudoRMGObject',
                                             'a': {'class': 'PseudoRMGObject',
                                                   'b': self.np_dict,
                                                   'c': {'class': 'PseudoRMGObject',
                                                         'c': self.array_dict,
                                                         'd': {'class': 'PseudoRMGObject',
                                                               'a': self.scalar_dict,
                                                               'b': {'class': 'PseudoRMGObject'}
                                                               }
                                                         }
                                                   }
                                             },
                                       'b': 6
                                       }

        self.list_of_objects = [1, 2.0, 'abc', self.np_array, self.array_quantity, self.scalar_quantity,
                                self.highly_nested_object]
        self.list_dict = [1, 2.0, 'abc', self.np_dict, self.array_dict, self.scalar_dict, self.highly_nest_dictionary]

        self.dictionary_of_objects = {'test_int': 1,
                                      'test_float': 2.0,
                                      'test_np': self.np_array,
                                      'test_array': self.array_quantity,
                                      'test_scalar': self.scalar_quantity,
                                      'test_nested': self.highly_nested_object}
        self.objects_dict = {'test_int': 1,
                             'test_float': 2.0,
                             'test_np': self.np_dict,
                             'test_array': self.array_dict,
                             'test_scalar': self.scalar_dict,
                             'test_nested': self.highly_nest_dictionary}

        self.input_dict = {'class': 'PseudoRMGObject', 'a': {'class': 'PseudoRMGObject', 'b': self.np_dict}}
        self.final_obj_dict = {'a': PseudoRMGObject(b=self.np_array)}

        self.class_dictionary = {'np_array': np.array,
                                 'ScalarQuantity': ScalarQuantity,
                                 'ArrayQuantity': ArrayQuantity,
                                 'PseudoRMGObject': PseudoRMGObject}

    def test_expanding_list_to_dict(self):
        """
        Test that objects nested inside of lists can be expanded
        """
        self.assertEqual(expand_to_dict(self.list_of_objects), self.list_dict)

    def test_expanding_objects_in_dictionary(self):
        """
        Test that objects nested inside of dictionaries can be expanded
        """
        self.assertEqual(expand_to_dict(self.dictionary_of_objects), self.objects_dict)

    def test_expanding_np_arrays(self):
        """
        Test that np_arrays are expanded properly
        """
        self.assertEqual(expand_to_dict(self.np_array), self.np_dict)

    def test_expanding_rmg_objects(self):
        """
        Test that RMGObjects (even when nested) can be expanded using the as_dict method
        """
        self.assertEqual(expand_to_dict(self.highly_nested_object), self.highly_nest_dictionary)
        self.assertEqual(self.highly_nested_object.as_dict(), self.highly_nest_dictionary)

    def test_make_object_from_dict(self):
        """
        Test that RMGObjects can be recreated from their dictionary representation
        """
        created_from_function = recursive_make_object(self.highly_nest_dictionary, self.class_dictionary)
        created_from_object = PseudoRMGObject.__new__(PseudoRMGObject)
        created_from_object.make_object(self.highly_nest_dictionary, self.class_dictionary)
        orig_obj = self.highly_nested_object

        for obj in (created_from_function, created_from_object):
            self.assertEqual(orig_obj.b, obj.b)
            self.assertEqual(type(orig_obj.a), type(obj.a))
            self.assertEqual(type(orig_obj.a.a), type(obj.a.a))
            self.assertTrue(np.array_equal(orig_obj.a.a.b, obj.a.a.b))
            self.assertEqual(type(orig_obj.a.a.c), type(obj.a.a.c))
            self.assertEqual(orig_obj.a.a.c.c.units, obj.a.a.c.c.units)
            self.assertTrue(np.array_equal(orig_obj.a.a.c.c.value, obj.a.a.c.c.value))
            self.assertEqual(type(orig_obj.a.a.c.d), type(obj.a.a.c.d))
            self.assertEqual(orig_obj.a.a.c.d.a.units, obj.a.a.c.d.a.units)
            self.assertTrue(orig_obj.a.a.c.d.a.value, obj.a.a.c.d.a.value)
            self.assertEqual(type(orig_obj.a.a.c.d.b), type(orig_obj.a.a.c.d.b))

    def test_make_all_but_final_object_from_dict(self):
        """
        Test the `make_final_object=False` option for the recursive_make_object function
        """
        final_obj_dict = recursive_make_object(self.input_dict, self.class_dictionary, make_final_object=False)
        self.assertTrue(np.array_equal(final_obj_dict['a'].b, self.final_obj_dict['a'].b))

    def test_float_creation(self):
        """
        Test that strings of floats are recreated as floats
        """
        obj = recursive_make_object('5.0', self.class_dictionary)
        self.assertEqual(obj, 5.0)
        self.assertEqual(type(obj), float)

    def test_int_creation(self):
        """
        Test that strings of ints are recreated as ints
        """
        obj = recursive_make_object('5', self.class_dictionary)
        self.assertEqual(obj, 5)
        self.assertEqual(type(obj), int)

    def test_np_array_creation(self):
        """
        Test that numpy arrays can be recreated from their dictionary representation
        """
        self.assertTrue(np.array_equal(recursive_make_object(self.np_dict, self.class_dictionary), self.np_array))

    def test_hashable_class_key_creation(self):
        """
        Test that dictionaries with hashable class instances as keys can be recreated
        """
        @dataclass(frozen=True)
        class Data:
            arg: str

        input_dict = {"Data(arg='test')": 'test'}
        class_dictionary = {'Data': Data}
        key, val = list(recursive_make_object(input_dict, class_dictionary).items())[0]

        self.assertIsInstance(key, Data)
        self.assertEqual(key, Data('test'))
        self.assertEqual(val, 'test')


################################################################################


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
