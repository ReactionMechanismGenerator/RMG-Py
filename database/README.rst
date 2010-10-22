The RMG Database
================

This project should contain the data, and the tools for processing the data.

If you just want the latest data in an RMG-Java compatible format, grab it from
the 'output' folder.


Details 
-------

The original RMG_database (with all its history) for the kinetics_groups it is in 'input'.
For the other data, which is currently not being processed, it is all in 'output'.

There's a mostly ready kinetics project. Simply run::

	$ python kinetics_data.py

and it will process the kinetics_groups data from `input` and place it in 
`output`. It reads in the library.py files, and creates rateLibrary.txt files
(it currently also creates new library.py files, just to check that it can).

When done, you'll have a complete RMG_database in output which can be used
directly with RMG-Java.

The results of this are stored in the repository, so that you don't have to
be able to run the script in order to get at the data (the script requires a working
installation of RMG-Py. See http://github.com/GreenGroup/RMG-Py ). 

There's also a half-finished thermo project.
