.. _output:

**************************
Analyzing the Output Files
**************************

You will see that RMG has created multiple output files and folders: ::
	:file:'/chemkin  input.py output.html  /pdep  /plot  restart.pkl  RMG.log  /solver  /species'
 
The :file:'/chemkin' folder will likely have a large number of chemkin formatted files. In general, these can be disregarded, as you will be mainly interested in :file:'chem.inp', the chemkin formatted input file with a species list, thermochemical database, and a kinetic mechanism. The file :file:'chem_annotated.inp' is provided as a means to help make sense of species syntax and information sources. In addition, a species dictionary,:file:'species_dictionary.txt', is generated. Either chemkin file, in addition to the dictionary, may be used as inputs in the tools section of this website to better visualize the species and reactions: http://rmg.mit.edu/simulate/chemkin  
(alternatively, you can open :file:'output.html')

The :file:'/pdep' folder will contain files associated with the pressure-dependent reactions that RMG has generated. These files are formatted as input files for CanTherm. 

RMG currently includes a solver for isothermal batch reactors. This is in fact a critical part of the model enlargement algorithm. If you have included simulations in your input file, the solutions will be located in :file:'/solver'. You will probably only be interested in the files with the largest number tags.  