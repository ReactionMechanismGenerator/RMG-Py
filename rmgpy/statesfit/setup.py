#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
# RMG - Reaction Mechanism Generator
#
# Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
# RMG Team (rmg_dev@mit.edu)
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the 'Software'),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
#
################################################################################

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('',parent_package,top_path)
    config.add_library('libstatesfit', sources=['math.f90','cases.f90','modes.f90','states.f90','dqed.f90'], extra_compiler_args=['-fbounds-check'] )
    config.add_extension('_statesfit', sources=['fit.f90'], libraries=['libstatesfit'], extra_compile_args=['-fbounds-check'])
    return config

################################################################################

from numpy.distutils.core import setup

setup(
    configuration = configuration
)

