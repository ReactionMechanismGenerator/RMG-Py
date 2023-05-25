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

import collections

import numpy as np

################################################################################


cdef class RMGObject(object):
    """
    This class provides a general `as_dict` method to help with yml file construction.
    Used as a parent class for other classes in RMG to inherit from.
    """
    def __init__(self):
        pass

    cpdef dict as_dict(self):
        """
        A helper function for dumping objects as dictionaries for YAML files
        
        Returns: 
            dict: A dictionary representation of the object
        """
        cdef dict output_dict
        cdef list all_attributes
        cdef str attr

        output_dict = dict()
        output_dict['class'] = self.__class__.__name__
        all_attributes = [attr for attr in dir(self) if not attr.startswith('_')]
        for attr in all_attributes:
            val = getattr(self, attr)

            # Check if val is a numpy array first to prevent elementwise comparisons with ''
            if isinstance(val, np.ndarray) or (val is not None and not callable(val) and val != ''):
                output_dict[attr] = expand_to_dict(val)

        return output_dict

    cpdef make_object(self, dict data, dict class_dict):
        """
        A helper function for constructing objects from a dictionary (used when loading YAML files)
        
        Args:
            data (dict): The dictionary representation of the object
            class_dict (dict): A mapping of class names to the classes themselves

        Returns: 
            None
        """
        kwargs = recursive_make_object(data, class_dict, make_final_object=False)
        for key in ['aux', 'mol']:
            if key in kwargs.keys():
                del kwargs[key]
        self.__init__(**kwargs)

cpdef expand_to_dict(obj):
    """
    Takes an object of any type (list, dict, str, int, float, RMGObject, etc.) and returns a dictionary representation 
    of the object (useful for __repr__ methods and outputting to YAML). The function works recursively to all objects
    nested in the original object.
    
    Args:
        obj (Any): Object to be represented by a dictionary

    Returns: 
        Any: dictionary representation of the object (dict, unless str, int, or float, which are returned as themselves)

    """
    if isinstance(obj, (list, tuple)):
        return [expand_to_dict(x) for x in obj]

    elif isinstance(obj, dict):
        # Create a new dictionary to store the expanded objects
        new_obj = dict()

        for key, value in obj.items():
            new_key = expand_to_dict(key)
            new_value = expand_to_dict(value)
            try:
                new_obj[new_key] = new_value
            except TypeError:
                # Check if the key is a hashable object and use its string representation if so
                if isinstance(key, collections.Hashable):
                    new_obj[repr(key)] = new_value
                else:
                    raise NotImplementedError(f'Cannot expand objects that are serving as dictionary keys ({key} is'
                                              'serving as a key). The returned dictionary would not be hashable')
        return new_obj

    elif isinstance(obj, np.ndarray):
        # Output as a list, but add a class entry so that it can be recreated
        return {'class': 'np_array', 'object':obj.tolist()}

    else:
        if hasattr(obj, 'as_dict'):
            return expand_to_dict(obj.as_dict())
        else:
            return obj

cpdef recursive_make_object(obj, class_dictionary, make_final_object=True):
    """
    Takes a dictionary representation of an object and recreates the object from the mapping of class strings to classes
    provided in class_dictionary. This operates recursively to recreate objects that were nested inside other objects.
    The function can either return the arguments for the make_object method to make the final (topmost) object or can 
    return the recreated final object.

    Args:
        obj (Any): A dictionary representation of an object to be recreated
        class_dictionary (dict): A dictionary mapping of class strings to classes
        make_final_object (bool): If True (default) the topmost object will be created and returned. Else, all nested
                                  objects will be recreated but only the keyword arguments needed to recreate the
                                  topmost object will be returned.

    Returns: 
        Any: recreated object (default) or dictionary of keyword arguments to recreate the final (topmost) object
    """
    if isinstance(obj, dict):

        # Create a new dictionary object to store recreated keys and values
        new_obj = dict()

        for key, value in obj.items():
            if key != 'class':
                new_key = recursive_make_object(key, class_dictionary)
                new_value = recursive_make_object(value, class_dictionary)
                new_obj[new_key] = new_value

        if 'class' in obj:  # This is a dictionary of an object to be created
            try:
                class_to_make = class_dictionary[obj['class']]
            except KeyError:
                raise KeyError('Class {0} must be provided in the class_dictionary: {1}'.format(obj['class'],
                                                                                                class_dictionary))

            args = {key:new_obj[key] for key in new_obj if key != 'class'}
            if make_final_object:
                if hasattr(class_to_make, 'make_object'):
                    created_obj = class_to_make.__new__(class_to_make)
                    created_obj.make_object(args, class_dictionary)
                else:
                    created_obj = class_to_make(**args)
                return created_obj
            else:
                return args

        else:  # Not a dictionary of an object that needs reconstructing, so just return the new dictionary
            return new_obj

    elif isinstance(obj, list):
        return [recursive_make_object(x, class_dictionary) for x in obj]

    elif isinstance(obj, str):
        # First check if the string is a boolean, as these otherwise convert to floats
        if obj == 'True':
            return True

        elif obj == 'False':
            return False

        # Next, check to see if the string needs to be converted to a float or int
        else:
            try:
                if '.' in obj:  # This is could be a float, but can't be an integer
                    return float(obj)
                else:
                    return int(obj)

            except (ValueError, TypeError):
                # Next, check if the string is a representation of a class instance,
                # which can occur if hashable classes are used as keys
                for class_name in class_dictionary.keys():
                    if class_name in obj:
                        try:
                            return eval(obj, {'__builtins__': None}, class_dictionary)
                        except NameError:  # Probably just included the class name as a comment
                            pass
                return obj  # If we made it here then obj must be just a string

    else:
        return obj
