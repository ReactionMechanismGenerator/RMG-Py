.. _diffModels:

****************
Model Comparison
****************

The script ``diffModels`` compares two RMG generated models to determine their differences.  
To use this method you will need the chemkin and species dictionary outputs from RMG. These can be found in the 
``chemkin`` folder from the directory of the ``input.py`` file used for the RMG run.  The syntax is as follows::

	python diffModels.py CHEMKIN1 SPECIESDICT1 THERMO1 CHEMKIN2 SPECIESDICT2 THERMO2

where ``CHEMKIN`` represents the chemkin input file (``chem00XX.inp``), ``SPECIESDICT``
is the species diectionary from RMG (``species_dictionary.txt``) and ``THERMO`` 
is input as the chemkin file again (``chem00XX.inp``).  The numbers (``1`` and ``2``) represent 
which model to each file is from.  
 
Output of each comparison is printed, and the method then produces a html file (``diff.html``)
for easy viewing of the comparison.  

This method is also available to use with a web browser from the RMG website: `Model Comparison Tool <http://rmg.mit.edu/simulate/compare>`_.

