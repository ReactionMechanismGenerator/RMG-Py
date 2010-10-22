********************
Creating Input Files
********************

Syntax
======

The format of CanTherm input files is based on Python syntax. In fact, CanTherm
input files are valid Python source code, and this is used to facilitate 
reading of the file. 

Each section is made up of one or more function calls, where parameters are 
specified as text strings, numbers, or objects. Text strings must be wrapped in
either single or double quotes. 

Model Chemistry
===============

The first item in the input file should be a ``modelChemistry()`` function,
which accepts a string describing the model chemistry. Currently the only
allowed model chemistries are ``'CBS-QB3'`` and ``'G3'``. CanTherm uses this
information to adjust the computed energies to the usual gas-phase reference
states. For example, below demonstrates how to specify CBS-QB3 as a model 
chemistry::

    modelChemistry("CBS-QB3")

Species
=======

Each species of interest must be specified using a ``species()`` function,
which accepts the following parameters:

====================== =========================================================
Parameter              Description
====================== =========================================================
``label``              A unique string label used as an identifier
``geomLog``            The path to the Gaussian log file containing the optimized geometry
``statesLog``          The path to the Gaussian log file containing the computed frequencies
``extSymmetry``        The external symmetry number for rotation
``freqScaleFactor``    A factor by which to scale all frequencies
``linear``             ``True`` if the molecule is linear, ``False`` if not
``rotors``             A list of :class:`HinderedRotor()` objects describing the hindered rotors
``atoms``              A dict associating atom symbols with the number of each atom in the molecule
``bonds``              A dict associating bond types with the number of each bond in the molecule      
====================== =========================================================

The geometry and states log files can be identical if you computed them in the
same Gaussian output. Allowed atom symbols for the ``atoms`` parameter are 
``'C'``, ``'H'``, ``'N'``, ``'O'``, and ``'P'``. Allowed bond types for the
``bonds`` parameter are ``'C-H'``, ``'C-C'``, ``'C=C'``, ``'C#C'``, ``'O-H'``,
``'C-O'``, ``'C-O'``, ``'C=O'``, ``'N#N'``, ``'O=O'``, ``'H-H'``, and
``'C#N'``. In both cases you can omit atoms and bonds not present in your
species, and their counts will be automatically set to zero.

Each :class:`HinderedRotor()` object requires the following parameters:

====================== =========================================================
Parameter              Description
====================== =========================================================
``scanLog``            The path to the Gaussian log file containing the scan
``pivots``             The indices of the atoms in the hindered rotor torsional bond
``top``                The indices of all atoms on one side of the torsional bond (including the pivot atom)
``symmetry``           The symmetry number for the torsional rotation
====================== =========================================================

The following is an example of a typical species item, based on ethane::

    species(
        label = 'ethane',
        geomLog = 'ethane_cbs.log',
        statesLog = 'ethane_cbs.log',
        extSymmetry = 2,
        freqScaleFactor = 0.99,
        linear = False,
        rotors = [
            HinderedRotor(scanLog='ethane_scan_1.log', pivots=[0,4], top=[0,1,2,3], symmetry=3),
        ]
        atoms = {'C': 2, 'H': 6},
        bonds = {'C-C': 1, 'C-H': 6},
    )

Transition State
================

Each transition state of interest must be specified using a 
``transitionState()`` function, which accepts exactly the same parameters as
the ``species()`` function described above. This is only required if you wish
to perform a kinetics computation.

The following is an example of a typical transition state item::

    transitionState(
        label = 'TS', 
        geomLog = 'H+C2H4.log', 
        statesLog = 'H+C2H4.log', 
        extSymmetry = 2,
        freqScaleFactor = 0.99,
        linear = False, 
        rotors = [],
        atoms = {'C': 2, 'H': 5},
        bonds = {'C-C': 1, 'C-H': 5},
    )

Reaction
========

Each reaction of interest must be specified using a ``reaction()`` function,
which accepts the following parameters: 

====================== =========================================================
Parameter              Description
====================== =========================================================
``label``              A unique string label used as an identifier
``reactants``          A list of strings indicating the labels of the reactant species
``products``           A list of strings indicating the labels of the product species
``transitionState``    The string label of the transition state
====================== =========================================================

This is only required if you wish to perform a kinetics computation. The
following is an example of a typical reaction item::

    reaction(
        label = 'H + C2H4 <=> C2H5',
        reactants = ['H', 'C2H4'],
        products = ['C2H5'],
        transitionState = 'TS',
    )

Thermodynamics Computations
===========================

Use a ``thermo()`` function to compute the thermodynamic parameters for a
species. Pass the string label of the species you wish to compute the 
thermodynamic parameters for and the type of thermodynamics model to
generate (either ``'Wilhoit'`` or ''`NASA`'' for a Wilhoit polynomial
model or NASA polynomial model). If you would like to see a plot of the
fitted thermodynamics, set the `plot` parameter to ``True``.

Below is a typical ``thermo()`` function::

    thermo('ethane', model='Wilhoit', plot=True)

Kinetics Computations
=====================

Use a ``kinetics()`` function to compute the kinetic parameters for a
reaction. Pass the string label of the reaction you wish to compute the 
reaction parameters for and the type of tunneling to use (``'Wigner'``,
``'Eckart'``, or ``''`` for no tunneling). A modified Arrhenius model will
automatically be fit to the generated rate coefficients. If you would like to 
see a plot of the fitted kinetics, set the `plot` parameter to ``True``.

Below is a typical ``kinetics()`` function::

    kinetics('H + C2H4 -> C2H5', tunneling='', plot=True)

Examples
========

Perhaps the best way to learn the input file syntax is by example. To that end,
a number of example input files and their corresponding output have been given
in the ``examples`` directory.
