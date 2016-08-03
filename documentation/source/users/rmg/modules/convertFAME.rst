.. _convertFAME:

***********************************
Convert FAME to CanTherm Input File
***********************************

This module is utilized to convert FAME file types (used in RMG-Java) to CanTherm objects (used in RMG-Py) for pressure
dependent calculations.  

FAME is an early version of the pdep code in CanTherm written in Fortran and used by RMG-Java. This script enables importing FAME input files into CanTherm. Note that it is mostly designed to load the FAME input files generated automatically by RMG-Java, and may not load hand-crafted FAME input files. If you specify a `moleculeDict`, then this script will use it to associate the species with their structures. ::

    python convertFAME.py fame_object
    
where ``fame_object`` is the FAME file used to be converted into the CanTherm object.

Some additional options involve adding an RMG dictionary to process with the file.  The syntax for this is ::

	python convertFAME.py -d RMG_dictionary.txt fame_object
  
where ``RMG_dictionary.txt`` is the dictionary to process with the file.  

A max energy cuttoff is also possible when converting the file formats. ::

	python convertFAME.py -d RMG_dictionary.txt -x value units value units fame_object

where ``value`` represents the max energy amount and ``units`` represents its units
