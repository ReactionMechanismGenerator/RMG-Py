.. _thermo:

**************************
Thermochemistry Estimation
**************************


This section gives in-depth descriptions of the methods used for determining thermochemistry of 
species.

Thermochemistry of species is obtained via three possible ways:

#. Species thermochemistry libraries
#. Group contribution methods (GC)
#. On-the-fly Quantum-chemical calculation of Thermochemical Properties (QMTP)

1) Species thermochemistry libraries
====================================
These databases contain thermochemical parameters for species. 
In these databases each entry contains an unambiguous definition of the species 
(through the adjacency list representation), 
along with a values for the thermochemistry in a format that allows 
the evaluation of each thermodynamic variable as a function of temperature.

RMG is shipped with a number of species thermochemistry libraries, located in the
'libraries' folder of RMG-database. More information on these species thermochemistry libraries can be found in
:ref:`thermoDatabase`.


2) Group contribution methods
=============================
When the thermochemistry of a species is not present in one of the available
species thermochemistry libraries, RMG needs to estimate thermochemistry. One way
to do so, is by using group contribution methods that estimate the thermochemistry 
of a molecule based on the sub-molecular fragments present in the molecule. 
The Benson group additivity framework is such an example of a group contribution method 
that has proven to provide accurate estimates of the ideal gas thermochemistry 
for a large range of molecules.

Benson's Group Additivity approach ([Benson]_), divides a molecule into functional 
groups, and the contribution of each functional group to the overall thermochemistry is 
included. For example, the molecule 2-methylnonane consists of three types of groups:

* 1 tertiary carbon atom
* 6 secondary carbon atoms
* 3 primary carbon atoms


.. image:: images/BensonGA.PNG
	:align: center
	
Thermochemistry for the molecule X is calculated by summing up the values for each
of the contributions. E.g.:

.. math::

	\Delta_fH_{298}^{o}(X) = \sum_{i}GAV(C_{i})
 
The term 'group additive value' (GAV) denotes a polyvalent (ligancy > 1) monoatomic central atom 
:math:`C_{i}` surrounded by its nearest-neighbor ligands.
 
Values for each central atomtype (e.g. "tertiary carbon atom") and its surrounding ligands can be found in the thermo
group database, named group.py, of RMG. More information can be found here: :ref:`thermoDatabase`.

NNIs
----
Besides the main group-centered (GAV) contributions, non-next-nearest
neighbor interactions (NNI) may also be important to take into account. NNIs are interactions between
atoms separated by at least 2 atoms, such as alkane 1,4-gauche, alkane 1,5 (cf. figure), alkene 1,4-gauche, 
alkene single and double cis, ene-yne cis and ortho interactions. 

.. image:: images/15Alkanae.png
	:align: center

As a result, thermochemistry of the molecule X is determined as :

.. math::

	\Delta_fH_{298}^{o}(X) = \sum_{i}GAV(C_{i}) + \sum_{j}NNI_{j}

RMG contains a database with NNIs, named gauche.py and int15.py. More information on the nature on the available NNIs, and corresponding values 
can be found here: :ref:`thermoDatabase`.

Ring Strain
-----------
To account for ring strain, ring strain corrections (RSC) were
introduced. Because there is no obvious relation between the
RSC and the ring structure, a specific RSC is required for every
type of ring. For example, due to the significant ring strained induced in norbornane (cf. figure), a
ring correction (RSC) needs to be added to the sum of the GAVs of the individual carbon atoms:

.. image:: images/norbornane.png
	:align: center

As a result, thermochemistry of the molecule X is determined as :

.. math::

	\Delta_fH_{298}^{o}(X) = \sum_{i}GAV(C_{i}) + RSC

RMG contains a database with single-ring corrections, 'ring.py' and polycyclic ring corrections,
'polycyclic.py'. More information on the nature on the available NNIs, and corresponding values 
can be found here: :ref:`thermoDatabase`.


Hydrogen Bond Increment (HBI) method
------------------------------------
Lay et al. [Lay]_ introduced the hydrogen bond increment (HBI) method to
predict thermochemical properties of radicals. In contrast to Benson’s method,
the HBI method does not use the group-additivity concept. The HBI
enthalpy of formation of a radical (R*) is calculated from the enthalpy
of formation of the corresponding parent molecule (R-H) by adding a
HBI to account for the loss of a hydrogen atom. Hence, for
standard enthalpies of formation the HBI is defined as

.. math::

	HBI = \Delta_fH_{298}^{o}(R^*) - \Delta_fH_{298}^{o}(R-H)
    	= BDE(R-H) - \Delta_fH_{298}^{o}(H)

Similar expressions are valid for the entropy and heat capacity.

The HBI method is the default method use to estimate thermochemistry of radicals. Thus, 
the effect of resonance stabilization on the enthalpy of the radical will be accounted for
through the corresponding HBI. For example, the HBI labeled as "C=CC=CCJ" will account
for the resonance present in 1,4-Pentadien-3-yl radical.

RMG contains a database for with HBIs, named radical.py. More information on the nature on the available HBIs, and corresponding values 
can be found here: :ref:`thermoDatabase`.


On-the-fly Quantum-chemical calculation of Thermochemical Properties (QMTP)
===========================================================================
An interface for performing on-the-fly quantum and force field calculations 
has been developed and integrated into RMG to complement the species thermochemistry databases and
group contribution methods [Magoon and Green]_. This interface is particularly interesting for the estimation of
thermochemistry of molecules that are not present in one of the species thermochemistry databases,
and which cannot be estimated with sufficient accuracy using the Benson group additivity framework. This
pertains specifically to polycyclic fused ring containing species, whose ring strain cannot be modeled using
the available ring corrections in RMG's ring strain correction databases.

The QMTP interface involves a number of steps, summarized in the figure below.

.. image:: images/QMTP.jpg
	:align: center

In a first step the graph representation is converted into a three-dimensional representation of the molecule
through the generation of 3D coordinates for the atoms in the molecule. This is accomplished using the UFF force field available in
RDKit. Next, the 3D atomic coordinates are sent to a computational chemistry package, either OpenMopac or Gaussian,
that calculates the thermochemistry of the given molecule "on-the-fly". Finally, the calculated thermochemistry data is sent back to RMG.


Supported QM packages, and levels of theory
-------------------------------------------

Currently, the following table shows an overview of the computational chemistry packages
and levels of theory supported in the QMTP interface of RMG.

.. table::

    ======================================= ========================================
    QM Package                              Supported Levels of Theory
    ======================================= ========================================
    OpenMopac				                semi-empirical (PM3, PM6, PM7)
    Gaussian03		              			semi-empirical (PM3)
    MM4	[Allinger]_			                molecular mechanics (MM4)
    ======================================= ========================================
  	
  	
4) References
============================================

.. [Benson] Benson, Sidney William. "Thermochemical kinetics." (1976)

.. [Lay] Lay, T.; Bozzelli, J.; Dean, A.; Ritter, E. J. Phys. Chem. 1995, 99,14514-14527

.. [Magoon and Green] Magoon, Gregory R., and William H. Green. "Design and implementation of a next-generation software interface for on-the-fly quantum and force field calculations in automated reaction mechanism generation." Computers & Chemical Engineering 52 (2013): 35-45.

.. [Allinger] Allinger, N. L., & Lii, J.-H. (2008). MM4(2008) and MM4(2003).