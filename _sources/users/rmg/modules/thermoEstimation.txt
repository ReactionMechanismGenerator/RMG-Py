.. _thermoModule:

************************
Thermo Estimation Module
************************

The thermo estimation module can be run stand-alone. An example input file for 
this module is shown below:

::

	database(
	    thermoLibraries = ['primaryThermoLibrary', 'GRI-Mech3.0']   
	)
	
	species(
	    label='Cineole',
	    structure=SMILES('CC12CCC(CC1)C(C)(C)O2'),
	)
	
	quantumMechanics(
	    software='mopac',#mopac or gaussian
	    method='pm3',#pm3, pm6, pm7
	    fileStore='QMfiles', # defaults to inside the output folder.
	    onlyCyclics = True,#True, False
	    maxRadicalNumber = 0, # 0, 1
	)

The ``database`` block is used to specify species thermochemistry libraries.
Multiple libraries may be created, if so desired.
The order in which the thermo libraries are specified is important: 
If a species appears in multiple thermo libraries, the first instance will
be used.

Please see Section :ref:`thermoDatabase` for details on editing the
thermo library. In general, it is best to leave the ThermoLibrary
set to its default value.  In particular, the thermodynamic properties for H and H2
must be specified in one of the primary thermo libraries as they cannot be estimated
by Benson's method.

For example, if you wish to use the GRI-Mech 3.0 mechanism [GRIMech3.0]_ as a ThermoLibrary in your model, the syntax will be::

	thermoLibraries = ['primaryThermoLibrary','GRI-Mech3.0']
 

This library is located in the 
:file:`$RMG/RMG-database/input/thermo/libraries` directory.  All "Locations" for the
ThermoLibrary field must be with respect to the :file:`$RMG/RMG-database/input/thermo/libraries`
directory.


The optional ``quantumMechanics`` block is used when quantum mechanical calculations are desired to determine thermodynamic parameters.
These calculations are only run if the molecule is not included in a specified thermo library.	
The ``software`` option accepts either the ``mopac`` or ``gaussian`` string.
The ``method`` option refers to the level-of-theory, which can either be ``pm3``,``pm6``, or ``pm7``.
A folder can be specified to store the files used in these calculations,
however if not specified this defaults to a `QMfiles` folder in the output folder.
The ``onlyCyclics`` option, if ``True``, only runs these calculations for cyclic species.
In this case, group contribution estimates are used for all other species.
The calculations are also only run on species with a maximum radical number set by the user.
If a molecule has a higher radical number, the molecule is saturated with hydrogen atoms, then 
quantum mechanical calculations with subsequent hydrogen bond incrementation is used to determine the
thermodynamic parameters.

Submitting a job is easy::

	python thermoEstimator.py input.py

We recommend you make a job-specific directory for each thermoEstimator simulation.