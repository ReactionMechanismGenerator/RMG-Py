*************************************************************************************
Creating Input Files for Thermodynamics and High-Pressure Limit Kinetics Computations
*************************************************************************************

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
which accepts a string describing the model chemistry. Currently the 
allowed model chemistries are:
``'CBS-QB3'``
``'G3'``
``'M08SO/MG3S*'`` * indicates that the grid size used in the [QChem] electronic structure calculation utilized 75 radial points and 434 angular points
``'CCSD(T)-F12/cc-pVnZ-F12'``  n = D, T, Q
``'CCSD(T)-F12/aug-cc-pVnZ-F12'``  n = D, T, Q
``'MP2_rmp2_pVnZ'``  n = D, T, Q
``'FCI/cc-pVnZ'``  n = D, T, Q
``'DFT_G03_b3lyp'``  a B3LYP calculation with a moderately large basis set
``'BMK/cbsb7'`` or  ``'BMK/6-311G(2d,d,p)'``

CanTherm uses this information to adjust the computed energies to the usual gas-phase reference
states by applying atom, bond and spin-orbit coupling energy corrections. This is particularly important for ``thermo()`` calculations (see below). The example below 
demonstrates how to specify CBS-QB3 as a model chemistry::

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
``scanLog``            The path to the Gaussian/Qchem log file containing the scan
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

Note that the atoms identified within the rotor section should correspond to the geometry indicated by
``geomLog``. 

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
        tunneling='Eckart'        
    )

Note that in the above example, 'Wigner' is also an acceptable method of estimating the 
quantum tunneling factor. 

Thermodynamics Computations
===========================

Use a ``thermo()`` function to compute the thermodynamic parameters for a
species. Pass the string label of the species you wish to compute the 
thermodynamic parameters for and the type of thermodynamics model to
generate (either ``'Wilhoit'`` or ''`NASA`'' for a Wilhoit polynomial
model or NASA polynomial model). A table of thermodynamic parameters will
also be displayed in the output file. 

Below is a typical ``thermo()`` function::

    thermo('ethane', model='Wilhoit')

Kinetics Computations
=====================

Use a ``kinetics()`` function to compute the high-pressure limit kinetic parameters for a
reaction.  If desired, define a desired temperature range and number of temperatures 
at which the high-pressure rate coefficient will be tabulated and saved to 
the outupt file. 3-parameter modified Arrhenius coefficients will automatically be fit 
to the computed rate coefficients. The quantum tunneling factor will also be displayed

Below is a typical ``kinetics()`` function::

    kinetics(    
    label = 'H + C2H4 <=> C2H5',
    Tmin = (400,'K'), Tmax = (1200,'K'), Tcount = 6, 
    Tlist = ([400,500,700,900,1100,1200],'K'),
    )

This is also acceptable::

    kinetics('H + C2H4 <=> C2H5')

Examples
========

Perhaps the best way to learn the input file syntax is by example. To that end,
a number of example input files and their corresponding output have been given
in the ``examples`` directory.

Troubleshooting and FAQs
========================

1) The network that CanTherm generated and the resulting pdf file show abnormally large
absolute values. What's going on?

    This can happen if the number of atoms and atom types is not properly defined or consistent in your input file(s).

Cantherm User Checklist
========================

Using cantherm, or any rate theory package for that matter, requires careful consideration and management of a large amount of data, files, and input parameters. As a result, it is easy to make a mistake somewhere. This checklist was made to minimize such mistakes for users:

- Do correct paths exist for pointing to the files containing the electronic energies, molecular geometries and vibrational frequencies?

For calculations involving pressure dependence:

- Does the network pdf look reasonable? That is, are the relative energies what you expect based on the input?

For calculations using internal hindered rotors:

- Did you check to make sure the rotor has a reasonable potential (e.g., visually inspect the automatically generated rotor pdf files)?
- Within your input files, do all specified rotors point to the correct files?
- Do all of the atom label indices correspond to those in the file that is read by the logger (GaussianLog, QchemLog, etc.)?
- Why do the fourier fits look so much different than the results of the ab initio potential energy scan calculations? This is likely because the initial scan energy is not at a minimum. One solution is to simply shift the potential with respect to angle so that it starts at zero and, instead of having CanTherm read a Qchem or Gaussian output file, have CanTherm point to a 'ScanLog' file. Another problem can arise when the potential at 2*pi is also not [close] to zero.       
