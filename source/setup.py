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

"""
This is the setup file for RMG.
"""

def setupCythonModules(name, version, description, author, author_email, 
	url, packages):
	"""
	Compile Cython modules into shared binary libraries. These appear to the
	end-user as normal Python modules, but execute much more rapidly because
	they are compiled rather than interpreted. Furthermore, the Cython modules
	add type information in certain places, which really gives Cython modules
	a speed boost when compiled.
	"""
	
	from distutils.core import setup
	from distutils.extension import Extension
	from Cython.Distutils import build_ext
	
	# Stop wasting my time compiling PowerPC-compatible C extensions on my intel Mac
	import distutils.sysconfig
	config = distutils.sysconfig.get_config_vars()
	for key,value in config.iteritems():
		location = str(value).find('-arch ppc')
		if location>=0:
			print "removing '-arch ppc' from %s"%(key)
			config[key] = value.replace('-arch ppc ','')
	# Tip from Lisandro Dalcin:
	# FYI, in Py>=2.6 the ARCHFLAGS environ var is honored by distutils, so
	# you should be able to
	#  export ARCHFLAGS='-arch i386' # or '-arch x86_64'
	# and get your job done.
	
	# Create annotated .html files for each of the cython modules
	# These show which bits of the code are not very C-friendly and help you optimize
	import Cython.Compiler
	Cython.Compiler.Options.annotate = True
	
	# The Cython modules to setup
	ext_modules = [
		Extension('rmg.chem', ['rmg/chem.py']),
		Extension('rmg.graph', ['rmg/graph.py']),
		Extension('rmg.thermo', ['rmg/thermo.py'])
	]

	setup(name=name,
		version=version,
		description=description,
		author=author,
		author_email=author_email,
		url=url,
		packages=packages,
		cmdclass = {'build_ext': build_ext},
		ext_modules = ext_modules
	)

################################################################################

def configuration(parent_package='',top_path=None):
	from numpy.distutils.misc_util import Configuration
	
	config = Configuration('rmg',parent_package,top_path)
	config.add_library('rmg_spectral_math', sources=['rmg/spectral/math.f90'], extra_compiler_args=['-fbounds-check'] )
	config.add_library('rmg_spectral_cases', sources=['rmg/spectral/cases.f90'], extra_compiler_args=['-fbounds-check'])
	config.add_library('rmg_spectral_dqed', sources=['rmg/spectral/dqed.f90'], extra_compiler_args=['-fbounds-check'])
	config.add_library('rmg_spectral_modes', sources=['rmg/spectral/_modes.f90'], extra_compiler_args=['-fbounds-check'])
	config.add_extension('spectral._modes', sources=['rmg/spectral/_modes.f90'], libraries=['rmg_spectral_math'], extra_compile_args=['-fbounds-check'])
	config.add_extension('spectral._fit', sources=['rmg/spectral/_fit.f90'], libraries=['rmg_spectral_cases', 'rmg_spectral_dqed', 'rmg_spectral_modes', 'rmg_spectral_math'], extra_compile_args=['-fbounds-check'])
	config.add_extension('spectral.states', sources=['rmg/spectral/states.f90'], libraries=['rmg_spectral_modes','rmg_spectral_math'], extra_compile_args=['-fbounds-check'])

	config.add_library('rmg_unirxn_mastereqn', sources=['rmg/unirxn/mastereqn.f90'], extra_compiler_args=['-fbounds-check'])
	config.add_extension('unirxn.mastereqn', sources=['rmg/unirxn/mastereqn.f90'], libraries=['blas', 'lapack'], extra_compile_args=['-fbounds-check'])
	config.add_extension('unirxn.msc', sources=['rmg/unirxn/msc.f90'], libraries=['blas', 'lapack'], extra_compile_args=['-fbounds-check'])
	config.add_extension('unirxn.rs', sources=['rmg/unirxn/rs.f90'], libraries=['blas', 'lapack'], extra_compile_args=['-fbounds-check'])
	config.add_extension('unirxn.cse', sources=['rmg/unirxn/cse.f90'], libraries=['rmg_unirxn_mastereqn', 'blas', 'lapack'], extra_compile_args=['-fbounds-check'])

	return config

def setupFortranModules(name, version, description, author, author_email, 
	url, packages):
	"""
	Compile Fortran files into shared binary libraries. These appear to the
	end-user as normal Python modules, but execute much more rapidly because
	they are compiled rather than interpreted.
	"""
	
	from numpy.distutils.core import setup
	
	setup(name=name,
		version=version,
		description=description,
		author=author,
		author_email=author_email,
		url=url,
		packages=packages,
		configuration = configuration
	)

################################################################################

if __name__ == '__main__':
	
	# Package details (passed to both setup() calls)
	name = 'RMG'
	version = '0.0.1'
	description = 'Reaction Mechanism Generator'
	author = 'Prof. William H. Green and the RMG Team'
	author_email = 'whgreen@mit.edu, rmg_dev@mit.edu'
	url = 'http://rmg.sourceforge.net/'
	packages = ['rmg']
	
	# Use f2py (in numpy) to compile the Fortran modules
	setupFortranModules(name, version, description, author, author_email, 
		url, packages)
	
	# Use Cython to compile the C modules
	setupCythonModules(name, version, description, author, author_email, 
		url, packages)
	