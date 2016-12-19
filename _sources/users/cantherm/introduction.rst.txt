************
Introduction
************

**CanTherm** is a tool for computing the thermodynamic properties of chemical
species and high-pressure-limit rate coefficients for chemical reactions using
the results of a quantum chemistry calculation. Thermodynamic properties are
computed using the rigid rotor-harmonic oscillator approximation with optional
corrections for hindered internal rotors. Kinetic parameters are computed using
canonical transition state theory with optional tunneling correction.

CanTherm can also estimate 
pressure-dependent phenomenological rate coefficients :math:`k(T,P)` for  
unimolecular reaction networks of arbitrary complexity. The approach is to
first generate a detailed model of the reaction network using the 
one-dimensional master equation, then apply one of several available model 
reduction methods of varying accuracy, speed, and robustness to simplify the 
detailed model into a set of phenomenological rate coefficients. The result 
is a set of :math:`k(T,P)` functions suitable for use in chemical reaction
mechanisms.


About CanTherm
==============

CanTherm is written in the `Python <http://www.python.org/>`_ programming
language to facilitate ease of development, installation, and use. 

License
=======

CanTherm is provided as free, open source code under the terms of the 
`MIT/X11 License <http://www.opensource.org/licenses/mit-license.php>`_. The 
full, official license is reproduced below



.. literalinclude:: ../../../../COPYING.txt