.. _generateFluxDiagram:

************************
Generating Flux Diagrams
************************

The script, ``generateFluxDiagrams.py``, will create a movie out of a completed RMG model
that shows interconnected arrows between species that represent fluxes.  

To use this method, you just need a Chemkin input file and an RMG species dictionary.
The syntax is as follows::

    python generateFluxDiagram.py [-h] [--no-dlim] [-s SPECIES] [-f]
                                  [-n N] [-e N] [-c TOL] [-r TOL] [-t S]
                                  INPUT CHEMKIN DICTIONARY [CHEMKIN_OUTPUT]

Positional arguments::

    INPUT                  RMG input file
    CHEMKIN                Chemkin file
    DICTIONARY             RMG dictionary file
    CHEMKIN_OUTPUT         Chemkin output file

Optional arguments::

    -h, --help             show this help message and exit
    --no-dlim              Turn off diffusion-limited rates
    -s DIR, --species DIR  Path to folder containing species images
    -f, --foreign          Not an RMG generated Chemkin file (will be checked for duplicates)
    -n N, --maxnode N      Maximum number of nodes to show in diagram
    -e N, --maxedge N      Maximum number of edges to show in diagram
    -c TOL, --conctol TOL  Lowest fractional concentration to show
    -r TOL, --ratetol TOL  Lowest fractional species rate to show
    -t S, --tstep S        Multiplicative factor to use between consecutive time points

This method is also available to use with a web browser from the RMG website: `Generate Flux Diagram <https://rmg.mit.edu/tools/flux>`_.

