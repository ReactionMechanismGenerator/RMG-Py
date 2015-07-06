.. _chemkin:

**************************
Working With Chemkin Files
**************************
This section provides information on how to import your own chemkin II-formatted mechanism files
into RMG so that they can be used as reaction libraries or seed mechanisms for further mechanism generation. 

Converting a Chemkin II file for use in RMG-Py
----------------------------------------------
This is useful when you have a mechanism (e.g., a 'seed' or 'foundation' mechanism) that you want RMG to build from.

In this case, :ref:`importChemkinLibrary.py` should be used to load a mechanism and its attending thermo and possibly transport data
to a user specified path (this script is located in `RMG-database`. Once loaded, it is convered into a python format compatible with RMG. The user must also provide a species dictionary in the correct and current adjacency list format. To run the script, type:

::

  python importChemkinLibrary.py chemkin.txt rmg_dictionary.txt library_name

chemkin.txt:
  the chemkin formated mechanism file containting ELEMENTS, SPECIES, REACTIONS, AND THERMO information
rmg_dictionary.txt:
  the RMG formatted dictionary of all species and their corresponding adjacency lists
library_name:
  the name you want to give your newly created RMG library (e.g., 'new_fuel')
