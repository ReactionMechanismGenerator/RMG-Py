.. _generateFluxDiagram:

************************
Generating Flux Diagrams
************************

The script, ``generateFluxDiagrams.py``, will create a movie out of a completed RMG model
that shows interconnected arrows between species that represent fluxes.  

To use this method, you just need a completed RMG run.  THe syntax is as follows::

	python generateFluxDiagram.py input.py
	
where ``input.py`` is the input file for the completed RMG run.  The program will use the automatically
generated file structure to find the other necessary files to create the movie.

This method is also available to use with a web browser from the RMG website: `Generate Flux Diagram <http://rmg.mit.edu/tools/flux>`_.

