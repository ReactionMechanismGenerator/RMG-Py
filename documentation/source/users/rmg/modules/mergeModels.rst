.. _mergeModels:

**************
Merging Models
**************

This script combines up to 5 RMG models together. The thermo and kinetics of common species and reactions are taken
from the first model with the commonality. To better understand the difference in two models, see :ref:`diffModels`.

To use this script type::

    python mergeModels.py --model1 chemkin1 speciesdict1 --model2 chemkin2 speciesdict2

where ``chemkin`` specifies the Chemkin input file from the RMG run and ``speciesdict`` represents the
species dictionary from an RMG run. These can be found in the ``chemkin`` folder from the directory of the
``input.py`` file used for the RMG run. The numbers are for different models that you want to merge. To merge more
than two files, you can add ``--model3 chemkin3 speciesdict3``, etc.. Up to 5 models can be merged together this way.

Running this method will create a new species dictionary (species_dictionary_merged.txt)
and Chemkin input file (chem_merged.inp) in the parent directory of the terminal.

To get a transport file for the merged mechanism, provide the transport file for any of the input mechanisms::

    python mergeModels.py --model1 chemkin1 speciesdict1 transport1 --model2 chemkin2 speciesdict2

This feature is also available to use directly via the RMG website:
`Model Merge Tool <http://rmg.mit.edu/tools/merge_models>`_.

An additional feature available only when running the script locally is to generate cross reactions between the
different models. In other words, reactions will be generated between species which are exclusive to one model and
species which are exclusive to another model. To use this feature, an additional input file must be provided,
similar to an RMG input file. The only required component is the ``database`` block (see the :ref:`input` section).
The input file is provided via the ``--react`` argument::

    python mergeModels.py -m1 chem1 dict1 -m2 chem2 dict2 --react input.py

By default, this will only add new reactions if all of the species already exist in the model. To change this behavior
and force all new species generated during this process to be included in the model, add the ``--add-all`` flag::

    python mergeModels.py -m1 chem1 dict1 -m2 chem2 dict2 --react input.py --add-all

Similar to normal RMG operation, the ``--quiet``, ``--verbose``, ``--debug``,``--maxproc``, and ``--walltime``
options are also supported.
