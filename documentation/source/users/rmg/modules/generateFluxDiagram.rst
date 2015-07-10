.. _generateFluxDiagram:

************************
Generating Flux Diagrams
************************

The script, ``generateFluxDiagrams.py``, will create a movie out of a completed RMG model
that shows interconnected arrows between species that represent fluxes.  

To use this method, you just need a completed RMG run.  THe syntax is as follows::

	python $RMGPy/generateFluxDiagram.py input.py
	
where ``input.py`` is the input file for the completed RMG run.  The program will use the automatically
generated file structure to find the other necessary files to create the movie.

