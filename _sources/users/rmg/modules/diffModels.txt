.. _diffModels:

****************
Model Comparison
****************

The script ``diffModels`` compares two RMG generated models to determine their differences.  
To use this method you will need the chemkin and species dictionary outputs from RMG. These can be found in the 
``chemkin`` folder from the directory of the ``input.py`` file used for the RMG run.  The syntax is as follows::

	python diffModels.py CHEMKIN1 SPECIESDICT1 --thermo1 THERMO1 CHEMKIN2 SPECIESDICT2 --thermo2 THERMO2 --web

where ``CHEMKIN`` represents the chemkin input file (``chem00XX.inp``), ``SPECIESDICT``
is the species diectionary from RMG (``species_dictionary.txt``) and  the optional ``--thermo`` flag can be used
to add separate thermo CHEMKIN files ``THERMO``.  The numbers (``1`` and ``2``) represent 
which model to each file is from.  The optional ``--web`` flag is used for running this script through the
RMG-website.  

Running the script without any optional flags looks like::

    python diffModels.py CHEMKIN1 SPECIESDICT1 CHEMKIN2 SPECIESDICT2

 
Output of each comparison is printed, and the method then produces a html file (``diff.html``)
for easy viewing of the comparison.  

This method is also available to use with a web browser from the RMG website: `Model Comparison Tool <http://rmg.mit.edu/simulate/compare>`_.

