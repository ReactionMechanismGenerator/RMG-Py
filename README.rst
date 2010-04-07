The RMG Database
================

This project should contain the data, and the tools for processing the data.

The original RMG_database (with all its history) is in currently in 'input'

There's a mostly ready kinetics project. Simply run::

	$ python kinetics_data.py

and it will process the kinetics_groups data from `input` and place it in 
`output`. It reads in the library.py files, and creates rateLibrary.txt files
(it currently also creates new library.py files, just to check that it can)


There's a half-finished thermo project.
