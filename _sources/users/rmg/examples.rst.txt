.. _examples:

*******************
Example Input Files
*******************

Perhaps the best way to learn the input file syntax is by example. To that end,
a number of example input files and their corresponding output have been given
in the ``examples`` directory. Two of the RMG jobs are shown below.


Ethane pyrolysis (Minimal)
==========================

This is the minimal example file characterizing a very basic system for ethane pyrolysis and should run quickly if RMG is set up properly. It does not include any 
calculation of pressure-dependent reaction rates.

.. literalinclude:: ../../../../examples/rmg/minimal/input.py


1,3-hexadiene pyrolysis
=======================

This example models the pyrolysis of 1,3-hexadiene and demonstrates the effect of turning on the pressure-dependence module within RMG.  

.. literalinclude:: ../../../../examples/rmg/1,3-hexadiene/input.py


Commented input file
=======================

This is a fully commented input file with all optional blocks
for new users to better understand the options of rmg input files

.. literalinclude:: ../../../../examples/rmg/commented/input.py