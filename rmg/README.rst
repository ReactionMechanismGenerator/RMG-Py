RMG-Py - the Reaction Mechanism Generator in Python
===================================================

This is the Python version of RMG - Reaction Mechanism Generator - a tool for 
automatically generating kinetic models of chemical reaction mechanisms,
developed primarily by researchers in Prof. Green's research group at the 
Massachusetts Institute of Technology. Details can be found at 
http://rmg.sourceforge.net/ or by emailing rmg_dev@mit.edu
 
Installation
------------
To install, there are various dependencies:

 * Python (versions 2.5 and above are known to work)
 * Cantera (1.8 works)
 * ...a few more


To install, you have to compile the Fortran and C++ extensions 
(actually, building the C++ is optional but greatly speed things up). 
To build the extensions do::

    rmg/source$ python setup.py build_ext --inplace

or use the Makefile::

    rmg/source$ make

You may find that the GNU Fortran compiler gfortran helps:
http://gcc.gnu.org/wiki/GFortranBinaries

