.. _generateReactions:

******************
Generate Reactions
******************

The script generateReactions.py generates reactions between all species mentioned in an input file.
To call this method type::

	python generateReactions.py Input_File

where ``Input_File`` is a file similar to a general RMG input file which contains all the species
for RMG to generate reactions between.  An example file is placed in ``$RMGPy/examples/generateReactions/input.py``

.. literalinclude:: ../../../../../examples/generateReactions/input.py

This method will produce an ``output.html`` file in the directory of ``input.py`` which contains the all the reactions produced between the species.  


This method is also available to use with a web browser from the RMG website: `Populate Reactions <http://rmg.mit.edu/tools/populate_reactions>`_.

