.. _reduction:

***********************************
Reducing an RMG output
***********************************

RMG's method of generating reactions between all species in a core mechanism and
including them in the resulting model 
is a robust process to obtain all chemistry. However, the huge number of cross reactions
lead to a non-sparse matrix, which can increase computational time when using
the resulting models in other simulations.

To help reduce the complexity of RMG produced mechanisms, a mechanism reduction
script was written that eliminates unimportant reactions up to a set threshold.
Though this method will reduce number of reactions and guarantee target species
concentrations at the given conditions are minimally affected, no guarantee is given
that it will result in optimally reduced mechanism.

To reduce an RMG job, you will need an additional file ``reduction_input.py``. 
This file contains two terms that tell the reduction algorithm what to do. The
example file located in ``rmgpy/reduction/test_data/minimal/chemkin`` is written
as followed. ::

    targets = ['ethane', 'C']
    tolerance = .05

``targets`` is a list of species labels whose concentration change should be minimized, and ``tolerance``
is the percent change the user can tolerate at the end of simulation. In the above
example, this would be 0.05%. 
Higher values of ``tolerance`` lead to fewer final reactions with more error in
output rates.

To run a simulation, type ::

    python $RMG/rmgpy/reduction/main.py input.py reduction_input.py chem_annotated.inp species_dictionary.txt

where ``input.py`` is the RMG input file, ``reduction_input.py`` contains the two
parameters above, ``chem_annotated.inp`` contains the reaction mechanism output
by RMG and ``species_dictionary.inp`` contains the reaction mechanism.

The algorithm will reduce the number of reactions until the tolerance is no 
longer met.