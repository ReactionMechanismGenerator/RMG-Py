************
Introduction
************

**MEASURE** - short for **M**\ aster **E**\ quation **A**\ utomatic **S**\ olver 
for **U**\ nimolecular **RE**\ actions - is a tool for estimating 
pressure-dependent phenomenological rate coefficients :math:`k(T,P)` for  
unimolecular reaction networks of arbitrary complexity. The approach is to
first generate a detailed model of the reaction network using the 
one-dimensional master equation, then apply one of several available model 
reduction methods of varying accuracy, speed, and robustness to simplify the 
detailed model into a set of phenomenological rate coefficients. The result 
is a set of :math:`k(T,P)` functions suitable for use in chemical reaction
mechanisms.

About MEASURE
=============

MEASURE was originally developed by Joshua W. Allen under the tutelage of
`Prof. William H. Green <http://web.mit.edu/greengp/>`_ in the 
`Department of Chemical Engineering <http://web.mit.edu/cheme/>`_ at the 
`Massachusetts Institute of Technology <http://web.mit.edu/>`_. It is a
component of the automatic reaction mechanism generation software 
`RMG <http://rmg.sourceforge.net/>`_, and provides functionality for estimating
pressure-dependent rate coefficients to that software. It is also useful as
a tool for exploring individual reaction networks in detail.

MEASURE is written in the `Python <http://www.python.org/>`_ programming
language to facilitate ease of development, installation, and use. Certain
computationally-intensive procedures have been adapted to the
`Cython <http://www.cython.org/>`_ programming language in order to improve the
execution time.

License
=======

MEASURE is provided as free, open source code under the terms of the 
`MIT/X11 License <http://www.opensource.org/licenses/mit-license.php>`_. The 
full, official license is reproduced below::

    Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu).

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
