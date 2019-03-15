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

import logging
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
        """
        cdef dict output_dict
        cdef list all_attributes
        cdef str attr

        output_dict = dict()
        output_dict['class'] = self.__class__.__name__
        all_attributes = [attr for attr in dir(self) if not attr.startswith('_')]
        for attr in all_attributes:
            val = getattr(self, attr)
            if val is not None and not callable(val) and val != '':
                output_dict[attr] = val
        for key, val in output_dict.iteritems():
            if isinstance(val, list) and val:
                if isinstance(val[0], RMGObject):
                    output_dict[key] = [v.as_dict() for v in val]
                else:
                    output_dict[key] = val
            elif isinstance(val, np.ndarray):
                output_dict[key] = val.tolist()
            elif not isinstance(val, (int, float, str, dict)) and val:
                # this is an object, call as_dict() again
                output_dict[key] = val.as_dict()
        return output_dict

    cpdef make_object(self, dict data, dict class_dict):
        """
        A helper function for constructing objects from a dictionary (used when loading YAML files)
        """
        for key, val in data.iteritems():
            if isinstance(val, dict) and 'class' in val:
                # Call make_object to make another object within the parent object
                class_name = val['class']
                del val['class']
                try:
                    class_to_make = class_dict[class_name]
                except KeyError:
                    raise KeyError("Class {0} must be provided in the 'class_dict' parameter "
                                   "to make the object.".format(class_name))
                obj = class_to_make()
                obj.make_object(val, class_dict)
                logging.debug("made object {0}".format(class_name))
                data[key] = obj
            elif isinstance(val, list) and val:
                if isinstance(val[0], dict) and 'class' in val[0]:
                    # Call make_object to make a list of objects within the parent object (as in Conformer.Modes)
                    data[key] = list()
                    for entry in val:
                        class_name = entry['class']
                        del entry['class']
                        try:
                            class_to_make = class_dict[class_name]
                        except KeyError:
                            raise KeyError("Class {0} must be provided in the 'class_dict' parameter "
                                           "to make the object.".format(class_name))
                        obj = class_to_make()
                        if class_name in ['LinearRotor', 'NonlinearRotor', 'KRotor', 'SphericalTopRotor', 'HinderedRotor',
                                          'FreeRotor'] and 'rotationalConstant' in entry and 'inertia' in entry:
                                # Either `rotationalConstant` or `inertia` should be specified for a rotor.
                                # Here both are specified, so we delete `inertia`.
                                del entry['inertia']
                        obj.make_object(entry, class_dict)
                        logging.debug("made object {0}".format(class_name))
                        data[key].append(obj)
                else:
                    # print 'not a class dict', val
                    data[key] = val
            elif isinstance(val, str):
                try:
                    float(val)
                except ValueError:
                    pass
        self.__init__(**data)
