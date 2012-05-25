.. _introduction:

************
Introduction
************

**Reaction Mechanism Generator (RMG)** is an automatic chemical reaction mechanism generator that constructs kinetic models composed of elementary chemical reaction steps using a general understanding of how molecules react.

Dependencies
============

.. NOTE::
	RMG has been tested on the Python 2.5, 2.6, and 2.7 releases; dependency issues render it incompatible with Python 3.x releases


Briefly, RMG depends on the following packages:

* **NumPy:** fast matrix operations
* **SciPy:** fast mathematical toolkit
* **matplotlib:** generating plots
* **guppy:** memory profiling tools
* **OpenBabel:** species format conversion
* **Cython:** compiling Python modules to C
* **quantities:** unit conversion
* **nose:** advanced unit test controls
* **Sphinx:** documentation generation
* **pydot:** interface to Dot graph language
* **cairo:** molecular diagram generation
* **psutil:** system utilization diagnostic tool
* **xlwt:** generating Excel output files
* **Graphviz:** generating flux diagrams
* **PyDAS:** differential algebraic system solver
* **PyDQED:** constrained nonlinear optimization
* **RMG-database:** thermodynamic and kinetic libraries

Refer to platform-specific instructions for details on the best ways to install these packages before attempting to build RMG-Py:

* `Windows <installation/windows.html>`_

License
=======

RMG is an open source program, available to the general public free of charge. The primary RMG code is distributed under the terms of the `MIT/X11 License <http://www.opensource.org/licenses/mit-license.php>`_. However, RMG has a number of dependencies of various licenses, some of which may be more restrictive. **It is the user's responsibility to ensure these licenses have been obtained.** ::

	RMG - Reaction Mechanism Generator

	Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
	RMG Team (rmg_dev@mit.edu)
	
	Permission is hereby granted, free of charge, to any person obtaining a
	copy of this software and associated documentation files (the 'Software'),
	to deal in the Software without restriction, including without limitation
	the rights to use, copy, modify, merge, publish, distribute, sublicense,
	and/or sell copies of the Software, and to permit persons to whom the
	Software is furnished to do so, subject to the following conditions:
	
	The above copyright notice and this permission notice shall be included in
	all copies or substantial portions of the Software.
	
	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
	DEALINGS IN THE SOFTWARE.