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
Contains classes relating to thermochemistry.
"""

import quantities as pq
import logging

import data

################################################################################




################################################################################

if __name__ == '__main__':
	
	datapath = '../data/RMG_database/'

	int15Database, gaucheDatabase, groupDatabase, otherDatabase, \
	        radicalDatabase, ringDatabase, primaryDatabase = \
			loadOldThermoDatabases(datapath)
			
	#saveNewThermoDatabases(datapath, int15Database, gaucheDatabase, \
	        #groupDatabase, otherDatabase, radicalDatabase, ringDatabase, \
			#primaryDatabase)
	
	
	#abrahamDatabase = loadDatabase(datapath + 'thermo/Abraham_Dictionary.txt', \
		#datapath + 'thermo/Abraham_Tree.txt', \
		#datapath + 'thermo/Abraham_Library.txt', True)
	
	#unifacDatabase = loadDatabase(datapath + 'thermo/Unifac_Dictionary.txt', \
		#datapath + 'thermo/Unifac_Tree.txt', \
		#datapath + 'thermo/Unifac_Library.txt', True)
	
	# Thermo data for H2
	#data = [0.0, 31.233, 6.895, 6.975, 6.994, 7.009, 7.081, 7.219, 7.720, 0, 0.0007, 0]
	#thermo = ThermoGAData()
	#thermo.fromDatabase(data, '')
	#Cp = thermo.heatCapacity(pq.Quantity(700, 'K')); Cp.units = 'kcal/(mol*K)'
	#print Cp

	