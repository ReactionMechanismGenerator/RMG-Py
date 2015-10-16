********************************************************
Creating Input Files for Pressure Dependent Calculations
********************************************************

Syntax
======

There are four parts to a pressure-dependent calculation input file, giving the species, transition states,
path reactions, reaciton network, and algorithm
parameters. The species section must come before the reaction section. Before
discussing each of these sections, a brief word on the general input file
syntax will be given.

Species Parameters
==================

Each species in the network must be specified using a ``species()`` block.
This includes all unimolecular isomers, bimolecular reactants and products,
and the bath gas(es). A species that appears in multiple bimolecular channels
need only be specified with a single ``species()`` block.

There are a number of required and optional parameters associated with a species
block:

======================= =========================== ====================================
Parameter               Required?                   Description
======================= =========================== ====================================
``label``               all species                 A unique string label used as an identifier
``structure``           all species except bath gas A chemical structure for the species defined using either SMILES or InChI
``E0``                  all species                 The ground-state energy (including zero-point energy)
``modes``               all species                 The molecular degrees of freedom (see below)
``spinMultiplicity``    all species                 The ground-state spin multiplicity (degeneracy), sets to 1 by default if not used
``opticalIsomers``      all species                 The number of optical isomers of the species, sets to 1 by default if not used
``molecularWeight``     all species                 The molecular weight, if not given it is calculated based on the structure
``collisionModel``      optional                    Transport data for the species, if available
``energyTransferModel`` optional                    Assigned with ``SingleExponentialDown`` model if available
``thermo``              optional                    Thermo data for the species
======================= =========================== ====================================

The ``label`` parameter should be set to a string with the desired user name for the species. ::

    label = 'nButanol'

The ``structure`` parameter is defined by either SMILES or InChI.  For instance, either representation is
acceptable for the acetone molecule: ::

    structure = SMILES('CC(C)=O')

    structure = InChI('InChI=1S/C3H6O/c1-3(2)4/h1-2H3')

The ``E0`` ground state energy should be given in the quantity format ``(value, 'units')``, using units of either ``kJ/mol``, ``kcal/mol``, ``J/mol``, or ``cal/mol``: ::

    E0 = (-34.6,'kcal/mol')

The `modes` parameter is required for all unimolecular isomers and all
bimolecular reactant channels. When specifying the ``modes`` parameter, define a list
with the following types of degrees of freedom.  To understand how to define these
degrees of freedom, please click on the links below:

**Translational degrees of freedom**

.. currentmodule:: rmgpy.statmech

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`IdealGasTranslation`    A model of three-dimensional translation of an ideal gas
=============================== ================================================



**Rotational degrees of freedom**

.. currentmodule:: rmgpy.statmech

=========================== ====================================================
Class                       Description
=========================== ====================================================
:class:`LinearRotor`        A model of two-dimensional rigid rotation of a linear molecule
:class:`NonlinearRotor`     A model of three-dimensional rigid rotation of a nonlinear molecule
:class:`KRotor`             A model of one-dimensional rigid rotation of a K-rotor
:class:`SphericalTopRotor`  A model of three-dimensional rigid rotation of a spherical top molecule
=========================== ====================================================


**Vibrational degrees of freedom**

.. currentmodule:: rmgpy.statmech

=========================== ====================================================
Class                       Description
=========================== ====================================================
:class:`HarmonicOscillator` A model of a set of one-dimensional harmonic oscillators
=========================== ====================================================


**Torsional degrees of freedom**

.. currentmodule:: rmgpy.statmech

=========================== ====================================================
Class                       Description
=========================== ====================================================
:class:`HinderedRotor`      A model of a one-dimensional hindered rotation
=========================== ====================================================

The ``spinMultiplicity`` is defined using an integer, and is set to 1 if not indicated 
in the ``species`` block. ::

    spinMultiplicity = 2
    
Similarly, the ``opticalIsomers`` is also defined using an integer, and is set to 1
if not used in the ``species`` block. ::

    opticalIsomers = 6
    
The ``molecularWeight`` parameter should be defined in the quantity format ``(value, 'units')``
, for example: ::

    molecularWeight = (44.04, 'g/mol')

If the ``molecularWeight`` parameter is not given, it is calculated by CanTherm based
off the chemical structure.

The ``collisionModel`` is defined with the transport data, if available, using a 
``TransportData`` object: ::

    collisionModel = TransportData(sigma=(3.70,'angstrom'), epsilon=(94.9,'K'))
    
   
The ``energyTransferModel`` model available is a ``SingleExponentialDown``.

* ``SingleExponentialDown`` - Specify ``alpha0``, ``T0`` and ``n`` for the
  average energy transferred in a deactiving collision

  .. math :: \left< \Delta E_\mathrm{down} \right> = \alpha_0 \left( \frac{T}{T_0} \right)^n

An example of a typical ``energyTransferModel`` block is: ::

    energyTransferModel = SingleExponentialDown(
            alpha0 = (0.5718,'kcal/mol'),
            T0 = (300,'K'),
            n = 0.85,
        )
        
        

The following is an example of a typical species item, based on the acetylperoxy
radical :math:`\ce{CH3C(=O)OO.}`::

    species(
        label = 'acetylperoxy',
        structure = SMILES('CC(=O)O[O]'),
        E0 = (-34.6,'kcal/mol'),
        modes = [
            IdealGasTranslation(mass=(75.04,"g/mol")),
            NonlinearRotor(inertia=([54.2977,104.836,156.05],"amu*angstrom^2"), symmetry=1),
            HarmonicOscillator(frequencies=([319.695,500.474,536.674,543.894,727.156,973.365,1037.77,1119.72,1181.55,1391.11,1449.53,1454.72,1870.51,3037.12,3096.93,3136.39],"cm^-1")),
            HinderedRotor(inertia=(7.38359,"amu*angstrom^2"), symmetry=1, fourier=([[-1.95191,-11.8215,0.740041,-0.049118,-0.464522],[0.000227764,0.00410782,-0.000805364,-0.000548218,-0.000266277]],"kJ/mol")),
            HinderedRotor(inertia=(2.94723,"amu*angstrom^2"), symmetry=3, fourier=([[0.130647,0.0401507,-2.54582,-0.0436065,-0.120982],[-0.000701659,-0.000989654,0.00783349,-0.00140978,-0.00145843]],"kJ/mol")),
        ],
        spinMultiplicity = 2,
        opticalIsomers = 1,
        molecularWeight = (75.04,"g/mol"),
        collisionModel = TransportData(sigma=(5.09,'angstrom'), epsilon=(473,'K')),
        energyTransferModel = SingleExponentialDown(
            alpha0 = (0.5718,'kcal/mol'),
            T0 = (300,'K'),
            n = 0.85,
        ),
    )


Transition States
=================

Transition states for reactions in the pressure dependent network should be defined very similarly to ``species``
using a ``transitionState`` block, however it has less parameters:


====================== =================================================================================
Parameter              Description 
====================== =================================================================================
``label``              A unique string label used as an identifier
``E0``                 The ground-state energy (including zero-point energy)
``modes``              The molecular degrees of freedom (same as for ``species``, see above)
``spinMultiplicity``   The ground-state spin multiplicity (degeneracy), sets to 1 by default if not used
``opticalIsomers``     The number of optical isomers of the species, sets to 1 by default if not used
``frequency``          The negative frequency of the first-order saddle point
====================== =================================================================================

An example of a ``transitionState`` block is shown below. ::

    transitionState(
        label = 'isom1',
        E0 = (-5.8,'kcal/mol'),
        modes = [
            IdealGasTranslation(mass=(75.04,"g/mol")),
            NonlinearRotor(inertia=([49.3418,103.697,149.682],"u*angstrom**2"), symmetry=1, quantum=False),
            HarmonicOscillator(frequencies=([148.551,306.791,484.573,536.709,599.366,675.538,832.594,918.413,1022.28,1031.45,1101.01,1130.05,1401.51,1701.26,1844.17,3078.6,3163.07],"cm^-1"), quantum=True),
        ],
        spinMultiplicity = 2,
        opticalIsomers = 1,
        frequency = (-1679.04,'cm^-1'),
    )


Path Reactions
==============

Each path reaction - a reaction directly connecting two molecular configurations
in the network - is specified using a ``reaction()`` block. The following
parameters are available:

====================== ==================== ============================================================================================================
Parameter              Required?            Description
====================== ==================== ============================================================================================================
``label``              All reactions        A name for the reaction 
``reactants``          All reactions        A list of reactant species
``products``           All reactions        A list of product species
``transitionState``    All reactions        The transition state
``kinetics``           Optional             The high pressure-limit kinetics for the reaction
``tunneling``          Optional             The type of tunneling model (either 'Eckhart' or 'Wigner') to use for tunneling through the reaction barrier
====================== ==================== ============================================================================================================

A typical reaction block might look like this. ::

    reaction(
        label = 'isom1',
        reactants = ['acetylperoxy'],
        products = ['hydroperoxylvinoxy'],
        transitionState = 'isom1',
        kinetics = Arrhenius(A=(2.65e6,'m^3/(mol*s)'), n=0.0, Ea=(0.0,'kcal/mol'), T0=(1,"K")),
        tunneling = 'Eckart',
    )

Note that the reactants and products must have been previously declared using a ``species`` block,
using the same name labels.  Transition states must also be previously declared using a
``transitionState`` block.


Network
=======

A declaration for the overall network must be given using the ``network`` block.

This includes setting the following paramters:

====================== ================================================================================
Parameter              Description
====================== ================================================================================
``label``              A name for the network
``isomers``            A list of species participating in unimolecular reaction channels
``reactants``          A list of the species that participate in bimolecular reactant channels
``bathGas``            A dictionary of bath gases and their respective mole fractions, adding up to 1.0
====================== ================================================================================

CanTherm is largely able to determine the molecular configurations that define
the potential energy surface for your reaction network simply by inspecting the
path reactions. However, you must indicate which unimolecular and bimolecular
configurations you wish to include in the master equation formulation; all
others will be treated as irreversible sinks.

Note that all species and bath gases used in the ``network`` block must have been 
previously declared with the same name labels in a previous ``species`` block in the
input file.

You do not need to specify the product channels (infinite sinks) in this
manner, as any configuration not marked as an isomer or reactant channel will
be treated as a product channel. An example of the ``network`` block is shown below. ::


    network(
        label = 'acetyl + O2',
        isomers = [
            'acetylperoxy',
            'hydroperoxylvinoxy',
        ],
        reactants = [
            ('acetyl', 'oxygen'),
        ],
        bathGas = {
            'nitrogen': 0.4,
            'argon': 0.6,
        }
    )

Algorithm Parameters
====================

The overall parameters for the pressure-dependence calculation must be defined in a
``pressureDependence`` block at the end of the input file.   The following parameters are necessary:



====================== ====================================================================================================================================================
Parameter              Description
====================== ====================================================================================================================================================
``label``              Use the name for the ``network`` declared previously
``method``             Method to use for calculating the pdep network. Use either 'modified strong collision', 'reservoir state', or 'chemically-significant eigenvalues'
``interpolationModel`` Select the output type for the pdep kinetics, either in 'chebyshev' or 'pdeparrhenius' (plog) format
``activeKRotor``       A flag indicating whether to treat the K-rotor as active or adiabatic
``activeJRotor``       A flag indicating whether to treat the J-rotor as active or adiabatic
====================== ====================================================================================================================================================

Additionally, temperature/pressure ranges and energy grain sizes must be given.

**Temperature and Pressure Ranges**

CanTherm will compute the :math:`k(T,P)` values on a grid of temperature and
pressure points. ``Tmin``, ``Tmax``, and ``Tcount`` values, as well as ``Pmin``, ``Pmax``, and ``Pcount`` parameter values must be provided.  
CanTherm will automatically choose the intermediate temperatures based on the
interpolation model you wish to fit. This is the recommended approach.


**Energy Grains**
Determine the fineness of the energy grains to be used in the master equation calculations.  Dictate
the ``maximumGrainSize``, and the ``minimumGrainCount``.


An example of the algorithm parameters block for the acetyl + O2 network is shown below. ::

    pressureDependence(
        label='acetyl + O2',
        Tmin=(300.0,'K'), Tmax=(2000.0,'K'), Tcount=8,
        Pmin=(0.01,'bar'), Pmax=(100.0,'bar'), Pcount=5,
        maximumGrainSize = (1.0,'kcal/mol'),
        minimumGrainCount = 250,
        method = 'modified strong collision',
        #method = 'reservoir state',
        #method = 'chemically-significant eigenvalues',
        interpolationModel = ('chebyshev', 6, 4),
        #interpolationModel = ('pdeparrhenius'),
        #activeKRotor = True, 
        activeJRotor = True,
    )


Examples
========

Perhaps the best way to learn the input file syntax is by example. To that end,
a number of example input files and their corresponding output have been given
in the ``examples/cantherm/networks`` directory, which includes both an `acetyl+O2`
and `n-butanol` example.




