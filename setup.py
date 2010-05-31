#!/usr/bin/env python
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
	
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import Cython.Compiler

# Create annotated HTML files for each of the Cython modules
Cython.Compiler.Options.annotate = True

# The Cython modules to setup
ext_modules = [
	Extension('chempy.constants', ['chempy/constants.py']),
	Extension('chempy.kinetics', ['chempy/kinetics.py']),
	Extension('chempy.reaction', ['chempy/reaction.py']),
	Extension('chempy.species', ['chempy/species.py']),
	Extension('chempy.thermo', ['chempy/thermo.py']),
]

setup(name='ChemPy',
	version='0.1.0',
	description='A chemistry toolkit for Python',
	author='Joshua W. Allen',
	author_email='jwallen@mit.edu',
	url='',
	packages=['chempy'],
	cmdclass = {'build_ext': build_ext},
	ext_modules = ext_modules,
)

