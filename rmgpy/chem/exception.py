#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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
This module contains exception classes for ChemPy-related exceptions. All such
exceptions should be placed within this module rather than scattered amongst
the others; this allows any ChemPy module that imports this one to see all of
the available ChemPy exceptions. Also, since this module contains only 
exception objecets, it is not among those that are compiled via Cython for 
speed.

All ChemPy exceptions derive from the base class :class:`ChemPyError`. This
base class can also be used as a generic exception, although this is generally
discouraged.
"""

################################################################################

class ChemPyError(Exception):
    """
    A generic ChemPy exception, and a base class for more detailed ChemPy
    exceptions. Contains a single attribute `msg` that should be used to
    provide information about the details of the exception.
    """

    def __init__(self, msg):
        self.msg = msg
    
    def __str__(self):
        return self.msg	

################################################################################

class InvalidThermoModelError(ChemPyError):
    """
    An exception used when working with a thermodynamics model to indicate that
    something went wrong while doing so.
    """
    pass

class InvalidKineticsModelError(ChemPyError):
    """
    An exception used when working with a kinetics model to indicate that
    something went wrong while doing so.
    """
    pass

class InvalidStatesModelError(ChemPyError):
    """
    An exception used when working with a states model to indicate that
    something went wrong while doing so.
    """
    pass
