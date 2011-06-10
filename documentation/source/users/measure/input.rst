********************
Creating Input Files
********************

There are four parts to a MEASURE input file, giving the species, important
unimolecular and bimolecular configurations, path reactions, and algorithm 
parameters. The species section must come before the reaction section. Before 
discussing each of these sections, a brief word on the general input file 
syntax will be given.

Syntax
======

The format of MEASURE input files is based on Python syntax. In fact, MEASURE
input files are valid Python source code, and this is used to facilitate 
reading of the file. 

Each section is made up of one or more function calls, where parameters are 
specified as text strings, numbers, or objects. Text strings must be wrapped in
either single or double quotes. There are two ways to specify numeric 
quantities:

* Dimensionless quantities can be specified as a simple literal, such as 
  ``0.99``.

* Dimensional quantities are specified using a tuple of the form 
  ``(value, units)``, where ``value`` is either a single number or a list of
  numbers, and ``units`` is a string, e.g. ``(-10.0,'kcal/mol')`` or
  ``([1.0, 2.0, 3.0],"m")``.

.. note::

    MEASURE uses the ``quantities`` package to convert your numeric parameters
    into SI units, and therefore inherits all of the idiosyncracies from that 
    package. In particular, the ``quantities`` package does *not*
    follow the SI convention that all units after the forward slash are in the
    denominator. For example, ``J/mol*K`` would be interpreted as ``(J/mol)*K``
    rather than ``J/(mol*K)``. Thus we recommend using parentheses where
    necessary to make your intentions explicit.

Species Parameters
==================

Each species in the network must be specified using a ``species()`` block.
This includes all unimolecular isomers, bimolecular reactants and products,
and the bath gas(es). A species that appears in multiple bimolecular channels 
need only be specified with a single ``species()`` block.

There are a number of required and optional parameters associated with a species
block:

====================== ==================== ====================================
Parameter              Required?            Description
====================== ==================== ====================================
``label``              all species          A unique string label used as an identifier
``E0``                 all species          The ground-state energy (including zero-point energy)
``states``             isomers, reactants   The molecular degrees of freedom (see below)
``lennardJones``       isomers, bath gas    The Lennard-Jones parameters, using a ``LennardJones`` call
``molecularWeight``    isomers, bath gas    The molecular weight
``thermo``                                  The macroscopic thermodynamic parameters
``SMILES``                                  The `SMILES <http://en.wikipedia.org/wiki/SMILES>`_ string describing the chemical structure
``InChI``                                   The `InChI <http://en.wikipedia.org/wiki/InChI>`_ string describing the chemical structure
====================== ==================== ====================================

If you specify the molecular structure via SMILES or InChI strings and omit the
molecular weight, the code will compute the molecular weight for you.

The `states` parameter is required for all unimolecular isomers and all
bimolecular reactant channels. When specifying the ``states`` parameter, use a 
``States()`` function with the following parameters:

=========================== ====================================================
Parameter                   Description
=========================== ====================================================
``rotations``               Parameters describing the external rotational motion, as a ``RigidRotor()`` object
``vibrations``              Parameters describing the internal vibrational motion, as a ``HarmonicOscillator()`` object
``torsions``                Parameters describing the internal torsional motion, as a list of ``HinderedRotor()`` objects
``frequencyScaleFactor``    The frequency scale factor to use (1.0 if not specified)
``spinMultiplicity``        The ground-state spin multiplicity (degeneracy)
=========================== ====================================================

The ``RigidRotor()``, ``HarmonicOscillator()``, and ``HinderedRotor()``
constructors match the corresponding classes in the :mod:`rmgpy.statmech`
module. The parameters for each are also summarized below:

=========================== ====================================================
Parameter                   Description
=========================== ====================================================
``RigidRotor()``
--------------------------- ----------------------------------------------------
`linear`                    ``True`` if the associated molecule is linear, ``False`` if nonlinear
`inertia`                   A list of the moment(s) of inertia of the molecule (1 if linear, 3 if nonlinear)
`symmetry`                  The total external rotational symmetry number
--------------------------- ----------------------------------------------------
``HarmonicOscillator()``
--------------------------- ----------------------------------------------------
`frequencies`               The set of vibrational frequencies
--------------------------- ----------------------------------------------------
``HinderedRotor()``
--------------------------- ----------------------------------------------------
`inertia`                   The reduced moment of inertia of the hindered rotor
`symmetry`                  The symmetry number for the hindered rotation
`barrier`                   The barrier height of the cosine potential
`fourier`                   The :math:`2 \times C` array of Fourier coefficients for the Fourier series potential
=========================== ====================================================

For each ``HinderedRotor()``, you need only specify one of the barrier height
or Fourier series coefficients.

If ``states`` is specified and ``thermo`` is not, then the thermodynamic
parameters will be automatically computed. This is recommended unless you have
thermodynamic data that you believe to be more accurate than the molecular
degrees of freedom data. You can use any of the thermodynamics models in the
``rmgpy.thermo`` module; see that package for more information on the available
models and their syntax.

The following is an example of a typical species item, based on the acetylperoxy
radical :math:`\ce{CH3C(=O)OO.}`::

    species(
        label='acetylperoxy',
        SMILES='CC(=O)O[O]',
        E0=(-34.6,'kcal/mol'),
        states=States(
            rotations=RigidRotor(
                linear=False,
                inertia=([54.2978, 104.8364, 156.0495],"amu*angstrom^2"),
                symmetry=1,
            ),
            vibrations=HarmonicOscillator(
                frequencies=([321.607, 503.468, 539.885, 547.148, 731.506, 979.187, 1043.981, 1126.416, 1188.619, 1399.432, 1458.200, 1463.423, 1881.701, 3055.285, 3115.447, 3155.144], 'cm^-1'),
            ),
            torsions=[
                HinderedRotor(inertia=(7.38359,"amu*angstrom^2"), barrier=(6.11665,"kcal/mol"), symmetry=1),
                HinderedRotor(inertia=(2.94725,"amu*angstrom^2"), barrier=(1.22157,"kcal/mol"), symmetry=3),
            ],
            frequencyScaleFactor=0.99,
            spinMultiplicity=2,
        ),
        lennardJones=LennardJones(sigma=(5.09,'angstrom'), epsilon=(473,'K')),
    )

Molecular Configurations
========================

MEASURE is largely able to determine the molecular configurations that define
the potential energy surface for your reaction network simply by inspecting the
path reactions. However, you must indicate which unimolecular and bimolecular
configurations you wish to include in the master equation formulation; all 
others will be treated as irreversible sinks. 

* For a unimolecular configuration, use the ``isomer()`` method, passing as the
  only parameter a string containing the label of the species to treat as a 
  unimolecular isomer.

* For a bimolecular configuration, use the ``reactants()`` method, passing as
  the two parameters a pair of strings containing the labels of the species to 
  treat as a bimolecular reactant channel. (A reactant channel is allowed to be
  both a source and a sink, i.e. both association and dissociation pathways are
  kept.)

For example, the following input specifies acetylperoxy as a unimolecular
isomer and acetyl + oxygen as a bimolecular reactant channel.
  
    isomer('acetylperoxy')

    reactants('acetyl', 'oxygen')

You do not need to specify the product channels (infinite sinks) in this 
manner, as any configuration not marked as an isomer or reactant channel will
be treated as a product channel.

Path Reaction Parameters
========================

Each path reaction - a reaction directly connecting two molecular configurations
in the network - is specified using a ``reaction()`` block. The following
parameters are available:

====================== ==================== ====================================
Parameter              Required?            Description
====================== ==================== ====================================
``reactants``          All reactions        A list of strings indicating the labels of the reactant species
``products``           All reactions        A list of strings indicating the labels of the product species
``transitionState``    All reactions        Information about the transition state, using a ``TransitionState()`` block; see below
``kinetics``                                The high pressure-limit kinetics for the reaction
====================== ==================== ====================================

The type of information specified with each path reaction determines how the
microcanonical rate coefficient is computed:

* If detailed information is known about the transition state, pass both the
  ground-state energy `E0` and the molecular degrees of freedom information
  `states` as parameters to the ``TransitionState()`` block. MEASURE will then
  use RRKM theory to compute the :math:`k(E)` values. (The molecular degrees of 
  freedom information is given in the same way as for species.)

* If only the high pressure-limit kinetics are known, pass only `E0` to the
  ``TransitionState()`` block, and provide the `kinetics` to the ``reaction()``
  block using an ``Arrhenius()`` block, where you specify the Arrhenius
  parameters ``A``, ``n``, ``Ea``, and optionally ``T0`` (set to 1 K if not
  explicity given). MEASURE will then use the inverse Laplace transform (ILT)
  method to compute the :math:`k(E)` values.
 
MEASURE will automatically use the best method that it can, so if you provide
both the molecular degrees of freedom and the high pressure-limit kinetics -
as in the example below - RRKM theory will be used.

The following is an example of a typical reaction item, based on the reaction
:math:`\ce{CH3C(=O)OO. -> CH2C=O + HO2}`::

    reaction(
        reactants=['acetylperoxy'],
        products=['ketene', 'hydroperoxyl'],
        kinetics=Arrhenius(
            A=(2.62e9,'s^-1'),
            n=1.24,
            Ea=(34.06,'kcal/mol')
        ),
        transitionState=TransitionState(
            E0=(0.6,'kcal/mol'),
            states=States(
                rotations=RigidRotor(
                    linear=False,
                    inertia=([55.4256, 136.1886, 188.2442],"amu*angstrom^2"),
                    symmetry=1,
                ),
                vibrations=HarmonicOscillator(
                    frequencies=([59.306,  205.421,  354.483,  468.861,  482.875,  545.574,  657.825,  891.898, 1023.947, 1085.617, 1257.494, 1316.937, 1378.552, 1688.566, 2175.346, 3079.822, 3154.325], 'cm^-1'),
                ),
                frequencyScaleFactor=0.99,
                spinMultiplicity=2,
            ),
            frequency=(-1048.9950,'cm^-1'),
        )
    )

Algorithm Parameters
====================

Collision Model
---------------

The collision model to use when constructing the master equation is specified
using a ``collisionModel()`` block, where you must specify both the ``type``
of model, the value of any required ``parameters`` as a list of quantities,
and a dictionary specifying the labels of the species in the bath gas and
their relative mole fractions. The collision models available are:

* ``single exponential down`` - Specify ``alpha0``, ``T0`` and ``n`` for the 
  average energy transferred in a deactiving collision

  .. math :: \left< \Delta E_\mathrm{down} \right> = \alpha_0 \left( \frac{T}{T_0} \right)^n

An example of a typical ``collisionModel()`` block is given below::

    collisionModel(
        type='single exponential down',
        parameters={
            'alpha0': (0.5718,'kcal/mol'),
            'T0': (300,'K'),
            'n': 0.85,
        },
        bathGas={
            'nitrogen': 1.0,
        }
    )

Temperature and Pressure Ranges
-------------------------------

MEASURE will compute the :math:`k(T,P)` values on a grid of temperature and
pressure points. The discussion below is for temperatures, but applies 
identically to pressures as well.

There are two ways to specify the temperature range using a ``temperatures()``
block:

* Give an explicit list of temperature points using the ``Tlist`` parameter::

    temperatures(Tlist=([300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0],'K'))

* Give the minimum temperature ``Tmin``, maximum temperature ``Tmax``, and 
  number of temperatures ``count`` to use::

    temperatures(Tmin=(300.0,'K'), Tmax=(2000.0,'K'), count=8)

  MEASURE will automatically choose the intermediate temperatures based on the 
  interpolation model you wish to fit. This is the recommended approach.

An example of typical ``temperatures()`` and ``pressures()`` blocks is given
below::

    temperatures(Tmin=(300.0,'K'), Tmax=(2000.0,'K'), count=8)
    pressures(Pmin=(0.01,'bar'), Pmax=(100.0,'bar'), count=5)

Energy Grains
-------------

Use an ``energies()`` block to specify information about the energies to use.
The required parameters are the minimum grain size ``dE`` and/or the minimum
number of grains ``count``. MEASURE will use whichever of these results in a
more accurate calculation. 

.. note::

    You do not need to specify the minimum and maximum energies, as MEASURE can
    determine these automatically.

A typical ``energies()`` block is given below::

    energies(dE=(0.25,'kcal/mol'), count=250)

Method
------

Use a ``method()`` block to specify the approximate method to use when computing
:math:`k(T,P)` values from the full master equation. There are currently three
methods available: ``modified strong collision`` (MSC), ``reservoir state`` 
(RS), and ``chemically-significant eigenvalues`` (CSE). The details of each of 
these methods is provided in the :doc:`../../theory/measure/index`. In short: 
MSC is the fastest but least accurate, RS usually provides a good balance of 
speed and accuracy (except at very high temperatures), while CSE is the most 
accurate but is slow and not robust. Currently we recommend using MSC during 
initial explorations, then RS when more accurate numbers are needed.

An example ``method()`` block is given below::

    method('reservoir state')

Interpolation Model
-------------------

Finally, use a ``model()`` block to specify the interpolation model to fit to
the computed :math:`k(T,P)` values. Currently there are two such models 
available:

* **Chebyshev polynomials.** You must also provide the number of terms to use in
  the temperature and pressure domains, respectively. The following example uses
  six Chebyshev terms in temperature and four in pressure::

    model('chebyshev', 6, 4)

  You should use fewer terms than the number of grid points in each direction,
  and should allow MEASURE to choose the intermediate temperature and pressure
  grid points.

* **Pressure-dependent Arrhenius.** A modified Arrhenius expression is fit at
  each pressure; no additional parameters are required::

    model('pdeparrhenius')

Of the two, Chebyshev polynomials are more flexible and therefore more likely
to fit the complex behavior of large networks. However, support for these models
varies in other chemical kinetics packages.

Examples
========

Perhaps the best way to learn the input file syntax is by example. To that end,
a number of example input files and their corresponding output have been given
in the ``examples`` directory.
