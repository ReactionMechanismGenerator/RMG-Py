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

A ``pdepreaction()`` item is printed for the forward and reverse direction of 
every reaction involving isomers and reactant channels only. For reactions 
involving a product channel, there is only a ``pdepreaction()`` item for the 
direction in which the product channel is the product of the reaction. To use
this output, you must either keep all of the reactions and treat them as
irreversible, or discard the duplicate reverse directions and treat the
remaining reactions as reversible. This decision is left to the end user.

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
