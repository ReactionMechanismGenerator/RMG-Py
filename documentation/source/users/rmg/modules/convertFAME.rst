.. _convertFAME:

******************
convert FAME to MEASURE
******************

This module is utilized to convert FAME file types (used in RMG-Java) to MEASURE objects (used in RMG-Py).  

FAME is an early version of MEASURE written in Fortran and used by RMG-Java. This script enables importing FAME input files into MEASURE so we can use the additional functionality that MEASURE provides. Note that it is mostly designed to load the FAME input files generated automatically by RMG-Java, and may not load hand-crafted FAME input files. If you specify a `moleculeDict`, then this script will use it to associate the species with their structures.

    python $RMGPy/convertFAME.py fame_object
    
where ``fame_object`` is the FAME file used to be converted into the MEASURE object.

Some additional options involve adding an RMG dictionary to process with the file.  The syntax for this is

	python $RMGPy/convertFAME.py -d RMG_dictionary.txt fame_object
  
where ``RMG_dictionary.txt`` is the dictionary to process with the file.  

A max energy cuttoff is also possible when converting the file formats.

	python $RMGPy/convertFAME.py -d RMG_dictionary.txt -x value units value units fame_object

where ``value`` represents the max energy amount and ``units`` represents its units