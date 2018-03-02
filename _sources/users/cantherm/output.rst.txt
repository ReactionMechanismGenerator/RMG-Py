********************
Parsing Output Files
********************

Thermodynamic and High-pressure Limit Kinetics Calculations
===========================================================

The syntax of CanTherm output files closely mirrors that of the input files.
For each ``thermo()`` function in the input file, there will be a corresponding
``thermo()`` function in the output file containing the computed thermodynamic
model. Similarly, For each ``kinetics()`` function in the input file, there will 
be a corresponding ``kinetics()`` function in the output file containing the
computed kinetics model.


Pressure-Dependent Calculations
===============================
The output file contains the entire contents of the input file. In
addition, the output file contains a block of ``pdepreaction()`` calls. The 
parameters of each ``pdepreaction()`` block match those of the ``reaction()`` 
block from the input file, except that no transition state data is given and 
the ``kinetics`` are by definition pressure-dependent.

A ``pdepreaction()`` item is printed for each reaction pathway possible in the
network. Each reaction is reversible. Reactions in the opposite direction are
provided as commented out, so a user can choose to use them if she/he desires.


Chemkin Output File
===================

In addition to the ``output.py`` which contains the thermodynamic,
kinetic, and pressure dependent results from a cantherm run, a Chemkin 
input file, ``chem.inp`` is also returned. This file contains species and their 
thermodynamic parameters for each species that has the ``thermo()`` in the 
input file. The file also contains kinetics, both pressure dependent and high 
pressure limit, which have the ``kinetics()`` or ``pressureDependence()`` module 
called.

For the output file to function, all the names of species should be in valid
chemkin format. The butanol and ethyl examples both show how to obtain a valid 
chemkin file.

The ``chem.inp`` file can be used in Chemkin software package or converted to 
a Cantera input file for use in Cantera software.


Log File
========

A log file containing similar information to that displayed on the console
during CanTherm execution is also automatically saved. This file has the name
``cantherm.log`` and is found in the same directory as the output file. The
log file accepts logging messages at an equal or greater level of detail than
the console; thus, it is often useful (and recommended) to examine both if
something unexpected has occurred.

The ``examples/cantherm`` directory contains both CanTherm input files and the resulting
output files.

Species Dictionary
==================

Any species that had the ``thermo()`` method called and had the structure defined in the cantherm
input file will also have an RMG style adjacency list representation in ``species_dictionary.txt``.
This allows the user to input the corresponding thermo and kinetics into RMG in various ways
described in the RMG user guide.