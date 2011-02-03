********************
Parsing Output Files
********************

The syntax of MEASURE output files closely mirrors that of the input files.
The output files contain two sections: 

* The first section contains all of the species involved in the reaction 
  network, including the bath gas(es), as ``species()`` blocks. This is largely
  a reproduction of the same section of the input file. 

* The second section contains all of the path reactions - reactions between any 
  pair of molecular configurations on the potential energy surface, not just 
  those directly adjacent - as ``pdepreaction()`` blocks. The parameters of each 
  ``pdepreaction()`` block match those of the ``reaction()`` block from the  
  input file, except that no transition state data is given and the ``kinetics``  
  are by definition pressure-dependent.

The ``examples`` directory contains both MEASURE input files and the resulting
output files.
