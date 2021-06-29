.. _mergeModels:

**************
Merging Models
**************

This script combines up to 5 RMG models together.  The thermo and kinetics from common species and reactions is taken
from the first model with the commonality.  To better understand the difference in two models, use diffModels.py.  
To use this method type::

	python mergeModels.py --model1 chemkin1 speciesdict1 --model2 chemkin2 speciesdict2

where ``chemkin`` specifies the chemkin input file from the RMG run and ``speciesdict`` represents the 
species dictionary from the RMG run.  These can be found in the 
``chemkin`` folder from the directory of the ``input.py`` file used for the RMG run.  
The numbers are for different models that you want to merge.  To merge more than two files, 
you can add ``--model3 chemkin3 speciesdict3``. Up to 5 models can be merged together this way

Running this method will create a new species dictionary (species_dictionary.txt) 
and chemkin input file (chem.inp) in the parent directory of the terminal.


This method is also available to use with a web browser from the RMG website: `Model Merge Tool <http://rmg.mit.edu/tools/merge_models>`_.

