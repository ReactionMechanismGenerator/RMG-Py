.. _reduction:

***********************************
Reaction Reduction in an RMG Job
***********************************

This script is located at ``RMG-Py/rmgpy/reduction/main.py`` instead of the usual
``RMG-Py/scripts`` folder.

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
example, this would be 5%. 
Higher values of ``tolerance`` lead to fewer final reactions with more error in
output rates.

To run a simulation, type ::

    python $RMG/rmgpy/reduction/main.py input.py reduction_input.py chem_annotated.inp species_dictionary.txt

A command line interface to the reduction driver script is contained in
``rmgpy/reduction/main.py``. It accepts four files: 

* ``input.py``: RMG-Py input file containing the settings to evaluate state variables.
* ``reduction_input.py``: Reduction input file containing the target variables and associated error tolerances to allow in the reduced model
* ``chem_annotated.inp``: the reaction mechanism to reduce.
* ``species_dictionary.txt``: the species dictionary associated with the reaction mechanism to reduce.

The algorithm will reduce the number of reactions until the tolerance is no 
longer met. If everything goes as planned, a ``chem_reduced.inp`` is generated
containing the reduced mechanism. In addition, a number of files
``chem_reduced_{i}.inp`` are created and correspond to the intermediate
reduced mechanisms. They can be used in place of the final reduced model, in case
the reduction algorithm does not terminate normally.

You can go to ``$RMG/examples/reduction`` to try this module.

Background
----------

The reduction algorithm computes the ratio of species reaction rate
(:math:`r_{ij}`) to the total rate of formation/consumption (:math:`R_i`) of all species i,
and compares this ratio to a tolerance (:math:`\epsilon`), with values of epsilon
between 0 and 1. If the ratio of a reaction is greater than epsilon it
is deemed *important* for the species in question. When a reaction is
not important for a single species, at any given time between t=0 and
the user-defined end time, then it is deemed unimportant for the given
system. As a result, the reaction is removed from the mechanism.

The value of epsilon is determined by an optimization algorithm that
attempts to reduce the model as much as possible given the constraints
of the user-defined target variables. A logarithmic bisection
optimization algorithm is used to provide guesses for the value of
epsilon based on the two previous guesses that undershoot and overshoot
the user-defined relative deviation of the target variables

A value of 5% for the relative deviation of the target variable implies
that the mole fraction of the target variable at the end time of the
batch reactor simulation as computed by the reduced mechanism may
deviate up to 5% w.r.t. to the mole fraction of the target variable at
the end time of the batch reactor simulation as computed by the full
mechanism.
