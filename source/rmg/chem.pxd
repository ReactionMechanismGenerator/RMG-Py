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

# these class definitions are required by graph.pyx
#

cdef class Atom:
	"""
	Represent an atom in a chemical species or functional group. The `atomType`
	and `electronState` attributes contain lists of the allowed atom types and
	electron states, respectively, for this atom. The `charge` attribute stores
	the resulting formal charge. The `label` attribute can be used to tag
	individual atoms, e.g. center atoms or reactive sites in functional groups.
	"""
	# should only be declared public if access is needed from Python (not just C)
	cdef  public list _atomType
	cdef  public list _electronState
	cdef  public int charge
	cdef  public str label
	cpdef bint equivalent(Atom, Atom)


cdef class Bond:
	"""
	Represent a bond in a chemical species. Each bond has a list `atoms` of
	length two containing the two atoms in the bond and a `bondType` object,
	stored internally as a :class:`BondType` object.
	"""
	cdef public list atoms
	cdef public list _bondType
	cpdef bint equivalent(Bond, Bond)


