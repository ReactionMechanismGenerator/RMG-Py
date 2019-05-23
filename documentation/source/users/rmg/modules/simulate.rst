.. _simulate:

***********************************************
Simulation and Sensitivity/Uncertainty Analysis
***********************************************

For sensitivity analysis, RMG-Py must be compiled with the DASPK solver, which is done by default but has
some dependency restrictions. (See :ref:`License Restrictions on Dependencies <dependenciesRestrictions>` for more details.) 
Sensitivity analysis or a simulation (without sensitivity) can be conducted in a standalone system for an existing
kinetics model in Chemkin format.

To run a simulation and/or sensitivity analysis, use the simulate module in RMG-Py/scripts::

    python simulate.py input.py chem.inp species_dictionary.txt
    
where ``chem.inp`` is the CHEMKIN file and the ``species_dictionary.txt`` contains the dictionary of
species associated with the CHEMKIN file.  ``input.py`` is an input file similar to one used for an RMG job but
does not generate a RMG job.  See the following ``input.py`` example file found under the
``$RMGPy/examples/sensitivity/input.py`` folder


.. literalinclude:: ../../../../../examples/rmg/minimal_sensitivity/input.py


The names of species named in the input file must coincide with the names specified in the CHEMKIN file.

Other options that can be specified for the ``simulate.py`` scripts are::

    --no-dlim      Turn off diffusion-limited rates for LiquidReactor
    -f, --foreign  Not an RMG generated Chemkin file (will be checked for duplicates)

Sensitivity analysis is conducted for the list of species given for the ``sensitivity`` argument in the input file.
The normalized concentration sensitivities with respect to the reaction rate coefficients dln(C_i)/dln(k_j) are saved to a csv file 
with the file name ``sensitivity_1_SPC_1.csv`` with the first index value indicating the reactor system and the second naming the index of the species 
the sensitivity analysis is conducted for.  Sensitivities to thermo of individual species is also saved as semi normalized sensitivities
dln(C_i)/d(G_j) where the units are given in 1/(kcal mol-1). The sensitivityThreshold is set to some value so that only
sensitivities for dln(C_i)/dln(k_j) > sensitivityThreshold  or dlnC_i/d(G_j) > sensitivityThreshold are saved to this file.  

Uncertainty analysis can also be requested via input file options. For more details, see :ref:`uncertaintyanalysis`.
The results of the analysis will be printed in the ``simulate.log`` file which is generated.
