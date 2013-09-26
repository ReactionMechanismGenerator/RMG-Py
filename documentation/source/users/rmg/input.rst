.. _input:

********************
Creating Input Files
********************

Syntax
======

The format of RMG-Py input files is based on Python syntax. 

Each section is made up of one or more function calls, where parameters are 
specified as text strings, numbers, or objects. Text strings must be wrapped in
either single or double quotes.

Datasources
===========

Species
=======
The following is an example of a typical species item, based on methane::

	species(
		label='CH4',
		reactive=True,
		structure=SMILES("C"),
	)

Reaction System
===============

On the fly Quantum Calculations
===============================

Pressure Dependence
===================

Miscellaneous Options
===================== 

Examples
========

Perhaps the best way to learn the input file syntax is by example. To that end,
a number of example input files and their corresponding output have been given
in the ``examples`` directory.
