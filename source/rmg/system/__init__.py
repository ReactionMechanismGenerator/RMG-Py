#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
The rmg.system package is where reaction systems, such as chemical reactors,
should be defined. Each reaction system should be its own class derived from 
a base class ReactionSystem. Each module in this package can contain either a
single reaction system or a set of related reaction systems.
"""

# Import all of the modules in the package
import base
import batch

def getAvailableReactionSystems():
	"""
	Return a list of the available reaction systems. This is done by inspecting
	all of the modules in this package for classes that derive from
	:class:`base.ReactionSystem`.
	"""
	import inspect
	import sys

	# Get a reference to the current package
	package = sys.modules[__name__]

	# Get a list of all of the modules in the current package
	modules = []
	for name in dir(package):
		obj = getattr(package, name)
		if inspect.ismodule(obj):
			modules.append(obj)

	# In each module, get a list of the available reaction systems
	availableSystems = {}
	for module in modules:
		for name in dir(module):
			obj = getattr(module, name)
			if inspect.isclass(obj):
				availableSystems[obj.__name__] = obj

	# Manually remove the abstract base class base.ReactionSystem
	del availableSystems['ReactionSystem']

	return availableSystems

################################################################################
