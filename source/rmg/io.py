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

import xml.dom.minidom

"""
Contains functions for manipulation of RMG input and output files.
"""

################################################################################

class InvalidInputFileException(Exception):
	"""
	An exception used when parsing an RMG input file to indicate that the input
	file is invalid. The msg parameter is used to specify what about the file
	caused the exception to be raised.
	"""	

	def __init__(self, msg):
		self.msg = msg
	
	def __str__(self):
		return 'Invalid XML for RMG input file: ' + self.msg

################################################################################

def readInputFile(fstr):
	"""
	Parse an RMG input file.
	"""

	try:
		
		# Parse the RMG input XML file into a DOM tree
		dom = xml.dom.minidom.parse(fstr)

		# Process root element (must be a rmginput element)
		root = dom.documentElement
		if root.tagName != 'rmginput':
			raise InvalidInputFileException('Incorrect root element.')
	
	except InvalidInputFileException, e:
		print 'ERROR: ' + str(e)
	except IOError, e:
		print 'ERROR: Input file "' + e.filename + '" not found.'
		
