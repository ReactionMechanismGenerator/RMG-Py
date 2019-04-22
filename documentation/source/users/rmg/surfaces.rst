.. _surfaces:

*****************************************************
Heterogeneous Catalysis Systems and Surface Reactions
*****************************************************

RMG can now be used to study heterogenous catalysis and surface reactions.
Initially developed in a fork of RMG called RMG-Cat ([Goldsmith2017]_),
this is now a part of the main RMG software.
Several surface specific features need to be considered when setting up an input file 
for surface reaction mechanism generation, as they cause certain aspects of it to 
deviate from the standard gas-phase RMG input file.


Reactor specifications
========================
For surface chemistry, RMG can model constant temperature and volume systems. 
The temperature, initial pressure, initial mole fractions of the reactant species, 
initial surface coverages, 
catalytic surface area to volume ratio in the reactor, 
and surface site density are defined for each individual reaction system in the input file.
As for the simple reactor model in RMG, the initial mole fractions are defined 
in a dictionary with the keys being names of species corresponding to molecules given
in the species block. 
The ``initialGasMoleFractions`` dictionary should contain gas-phase species,
the ``initialSurfaceCoverages`` dictionary should contain adsorbates and vacant sites.
Both will be normalized if the values given do not sum to 1.00.
The surface site density (``surfaceSiteDensity``) is the amount of active catalytic sites per unit surface area, 
and varies depending on the catalyst in question, but is held constant across simulations.
The ratio of catalyst surface area to gas phase volume (`surfaceVolumeRatio`) is determined by reactor geometry,
and so may be different in each ``surfaceReactor`` simulation. 

The following is an example of a surface reactor system for catalytic combustion of methane over a Pt catalyst::

    surfaceReactor(
        temperature=(800,'K'),
        initialPressure=(1.0, 'bar'),
        initialGasMoleFractions={
            "CH4": 0.0500,
            "O2": 0.1995,
            "N2": 0.7505,
        },
        initialSurfaceCoverages={
            "X": 1.0,
        },
        surfaceVolumeRatio=(1.0e4, 'm^-1'),
        surfaceSiteDensity=(2.72e-9, 'mol/cm^2'),
        terminationConversion = { "CH4":0.9 },
        terminationRateRatio=0.01
    )

                       
Adsorbate representation
-------------------------
Adsorbates are represented using chemical graph theory (ChemGraph), 
such that atoms are represented by nodes and bonds between them by edges. 
This way each adsorbate is represented by an :ref:`adjacency list <rmgpy.molecule.adjlist>`.
Specific to heterogeneous chemistry and surface reactions is the metal-adsorbate bond for adsorbates. 
Firstly, there needs to be a species representation for the binding site in the RMG-Cat input file. 
This can be added as follows::

    species(
        label='X',
        reactive=True,
        structure=adjacencyList("1 X u0"),
    )

The surface binding sites have an element ``X`` which 
has two atom-types:  either vacant, ``Xv``, 
or occupied, ``Xo`` in which case it has a molecule adsorbed on it.

Adsorbates in RMG-Cat can currently have one or two surface binding sites, 
and can be bound to the sites with single, double, triple, or quadruple bonds. 
The following is an example for the adjacency list of adsorbed methoxy with one binding site::

    1 X  u0 p0 c0 {3,S}
    2 C  u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
    3 O  u0 p2 c0 {1,S} {2,S}
    4 H  u0 p0 c0 {2,S}
    5 H  u0 p0 c0 {2,S}
    6 H  u0 p0 c0 {2,S}
    
Additionally, the adsorbates in RMG-Cat can be physisorbed (have a van der Waals bond to the surface). 
The adjacency lists for physisorbed species are structured the same way as for other adsorbates, 
except that there is no bond edge for the binding. 
The atom representing the surface site (``X``) must still be included.
Following is an adjacency list for physisorbed methane::

    1 X  u0 p0 c0
    2 C  u0 p0 c0 {3,S} {4,S} {5,S} {6,S}
    3 H  u0 p0 c0 {2,S}
    4 H  u0 p0 c0 {2,S}
    5 H  u0 p0 c0 {2,S}
    6 H  u0 p0 c0 {2,S}

This implementation currently does not distinguish between different binding sites on a given facet. 
Instead, the lowest energy binding site for a given adsorbate is assumed, 
consistent with the mean-field approach of the kinetics ([Goldsmith2017]_).

Thermochemistry
================
RMG will first check thermochemistry libraries for adsorbates.
Failing that, the gas phase thermochemistry will be used and an adsortion correction added.
The gas phase thermochemistry will be estimated using the methods specified for regular
gas phase species (libraries, automated quantum mechanics, machine learning, group additivity, etc.)
and the adsorption correction estimated as described below.
Finally, the adsorbed species will have its energy changed using linear scaling relationships,
to allow for metals other than Platinum(111) to be simulated.
These methods are all described below.

Use of thermo libraries for surface systems
---------------------------------------------
For surface species, thermo libraries provided in the input files are checked first 
for an exact match of a given adsorbate, and those thermodynamic properties are used. 

In order to predict the thermodynamic properties for species that are not in the database, 
RMG-Cat uses a precompiled adsorption correction with the thermodynamic properties 
of the gas-phase precursor ([Goldsmith2017]_).

Following is an example for how a thermo library for species adsorbed on platinum 
is provided in the input file database module::

    thermoLibraries=['surfaceThermoPt']
    
This can be added along with other gas-phase reaction libraries for coupling 
of gas-phase and surface reactions.  
For a full list of libraries check https://rmg.mit.edu/database/thermo/libraries/
or the folder ``RMG-database/input/thermo/libraries/`` in your RMG database.


Adsorption correction estimation
--------------------------------
The folder ``RMG-database/input/thermo/groups/`` contains the adsorption corrections 
for the change in thermodynamic properties of a species upon adsorption from the gas-phase. 
Currently the database has adsorption corrections for nickel 
(adsorptionNi.py) and platinum (adsorptionPt.py).

An example of an adsorption correction entry is shown below::

    entry(
        index = 40,
        label = "C-*R3",
        group =
    """
    1 X  u0 p0 c0 {2,S}
    2 C  u0 p0 c0 {1,S} {3,[S,D]} {4,[S,D]}
    3 R  u0 px c0 {2,[S,D]}
    4 R  u0 px c0 {2,[S,D]}
    """,
        thermo=ThermoData(
            Tdata=([300, 400, 500, 600, 800, 1000, 1500], 'K'),
            Cpdata=([-0.45, 0.61, 1.42, 2.02, 2.81, 3.26, 3.73], 'cal/(mol*K)'),
            H298=(-41.64, 'kcal/mol'),
            S298=(-32.73, 'cal/(mol*K)'),
        ),
        shortDesc=u"""Came from CH3 single-bonded on Pt(111)""",
        longDesc=u"""Calculated by Katrin Blondal at Brown University using statistical mechanics 
                (files: compute_NASA_for_Pt-adsorbates.ipynb and compute_NASA_for_Pt-gas_phase.ipynb). 
                Based on DFT calculations by Jelena Jelic at KIT.
                DFT binding energy: -1.770 eV.
                Linear scaling parameters: ref_adatom_C = -6.750 eV, psi = -0.08242 eV, gamma_C(X) = 0.250.

       CR3
       |
    ***********
    """
    )

Here, R in the label ``C-*R3`` represents any atom, where the H atoms in methyl have been replaced by wild cards.
This enables RMG-Cat to determine which species in the thermo database is the closest 
match for the adsorbate in question, using a hierachical tree for functional groups. 
This is defined at the bottom of the adsorption corrections files, e.g.::

    tree(
    """
    L1: R*
        L2: R*single_chemisorbed
            L3: C*
                L4: Cq*
                L4: C#*R
                    L5: C#*CR3
                    L5: C#*NR2
                    L5: C#*OR
                L4: C=*R2
                    L5: C=*RCR3
                    L5: C=*RNR2
                    L5: C=*ROR
                L4: C=*(=R)
                    L5: C=*(=C)
                    L5: C=*(=NR)
                L4: C-*R3
                    L5: C-*R2CR3
                    ...

When RMG-Cat has found the closest match, it reads the corresponding adsorption correction 
from the database and uses it with the thermo of the original adsorbate's 
gas-phase precursor to estimate its enthalpy, entropy and heat capacity.



Linear scaling relations
--------------------------
In surface reaction mechanism generation with RMG-Cat, 
linear scaling relations are used to investigate surface reaction systems occurring on surfaces 
that do not have DFT-based values in the database ([Mazeau2019]_). 
This is especially useful for alloy catalysts as conducting DFT calculations for such systems is impractical. 
Linear scaling relations for heterogeneous catalysis are based on the finding of 
Abild-Pedersen et al. ([Abild2007]_) that the adsorption energies of hydrogen-containing molecules 
of carbon, oxygen, sulfur, and nitrogen on transition metal surfaces scale linearly with the adsorption energy 
of the surface-bonded atom. 
Using this linear relationship, the energy of a species (`AH`\ :sub:`x`\ ) on any metal M2 can be estimated 
from the known energy on metal M1, :math:`\Delta E_{M1}^{AH_x}`, and the adsorption energies of 
atom A on the two metals M1 and M2 as follows:

.. math:: \Delta E_{M2}^{AH_x}=\Delta E_{M1}^{AH_x}+\gamma(x)(\Delta E_{M2}^A - \Delta E_{M1}^A),    
    :label: LSReqn
    
where

.. math:: \gamma (x)=(x_{max}-x)/x_{max}, 
    :label: gammaeqn

is the is the slope of the linear relationship between (`AH`\ :sub:`x`\ ) and A, and (`x`\ :sub:`max`\ ) 
is the maximum number of hydrogen atoms that can bond to the central atom A.
Since the adsorption energy of (`AH`\ :sub:`x`\ ) is proportional to adsorption energies on different metals, 
full calculations for every reaction intermediate on every metal are not necessary. 
Therefore, having this implemented in RMG-Cat allows for independent model generation for any metal surface. 
By effect it enables the expedient, high-throughput screening of catalysts for any surface 
catalyzed reaction of interest ([Mazeau2019]_).


Because of this feature,
it is required to provide the adsorption energies of C, N, O and H on the surface 
being investigated in the input file for RMG-Cat to generate a mechanism. 
The following is an example using the default binding energies of the four atoms on Pt(111).
Deviating from these values will result in adsorption energies being modified, 
even for species taken from the thermochemistry libraries::

    bindingEnergies = { # default values for Pt(111)
                       'H':(-2.479, 'eV/molecule'),
                       'O':(-3.586, 'eV/molecule'),
                       'C':(-6.750, 'eV/molecule'),
                       'N':(-4.352, 'eV/molecule'),
                       }


Reactions and kinetics
========================

Reaction families and libraries for surface reaction systems
------------------------------------------------------------
In the latest version of the RMG database, surface reaction families have been added. 
These include adsorption/desorption, bond fission and H abstraction ([Goldsmith2017]_). 
For surface reaction families to be considered in the mechanism generation, 
the 'surface' kinetics family keyword needs to be included in the database 
section of the input file as follows::

    kineticsFamilies=['surface', 'default']
    
This allows for RMG-Cat to consider both surface and gas reaction families. 
If inlcuding only surface reactions is desired, that can be attained by removing the 'default' keyword.

For surface reactions proposed by reaction families that do not have an exact 
match in the internal database of reactions, Arrhenius parameters are estimated 
according to a set of rules specific to that reaction family. 
The estimation rules are derived automatically from the database of known rate 
coefficients and formulated as Brønsted-Evans-Polanyi relationships ([Goldsmith2017]_).

The user can provide a surface reaction library containing a set of preferred rate coefficients for the mechanism. 
Just like for gas-phase reaction libraries, values in the provided reaction library are 
automatically used for the respective proposed reactions. 
The reactions in the reaction library are not required to be a part of the 
predefined reaction families ([Goldsmith2017]_).

Following is an example where a mechanism for catalytic partial oxidation of methane 
on platinum by Quiceno et al. ([Deutschmann2006]_) is provided as a reaction library 
in the database section of the input file::

    reactionLibraries = [('Surface/CPOX_Pt/Deutschmann2006', False)]   

Gas-phase reaction libraries should be included there as well 
for accurate coupled gas-phase/surface mechanism generation.

The following is a list of the current pre-packaged surface reaction libraries in RMG-Cat:

+-------------------------------------------------------------+------------------------------------------------------------------------------------------+
|Library                                                      |Description                                                                               |
+=============================================================+==========================================================================================+
|Surface/Deutschmann_Ni                                       |Steam- and CO2-Reforming as well as Oxidation of Methane over Nickel-Based Catalysts      |
+-------------------------------------------------------------+------------------------------------------------------------------------------------------+
|Surface/CPOX_Pt/Deutschmann2006                              |High-temperature catalytic partial oxidation of methane over platinum                     |
+-------------------------------------------------------------+------------------------------------------------------------------------------------------+




Examples
=========

Example input file: methane steam reforming
-------------------------------------------------------

This is a simple input file steam reforming of methane 

.. literalinclude:: ../../../../examples/rmg/catalysis/methane_steam/input.py

Example input file: methane oxidation
--------------------------------------------------------

This is an input file for catalytic partial oxidation (CPOX) of methane

.. literalinclude:: ../../../../examples/rmg/catalysis/ch4_o2/input.py


Additional Notes
-----------------
Other things to update:
* table of atom types in users/rmg/database/introduction.rst 
* table of atom types in reference/molecule/atomtype.rst


.. [Goldsmith2017] \ C.F. Goldsmith and R.H. West. "Automatic Generation of Microkinetic Mechanisms for Heterogeneous Catalysis." *J. Phys. Chem. C.* **121(18)**, p. 9970–9981 (2017).

.. [Deutschmann2006] \ R. Quiceno, J. Pérez-Ramírez, J. Warnatz and O. Deutschmann. "Modeling the high-temperature catalytic partion oxidation of methane over platinum gauze: Detailed gas-phase and surface chemistries coupled with 3D flow field simulations." *Appl. Catal., A* **303(2)**, p. 166-176 (2006).

.. [Mazeau2019] \ E.J. Mazeau, P. Satupte, K. Blondal, C.F. Goldsmith and R.H. West. "Linear Scaling Relationships and Sensitivity Analyses in RMG-Cat." *Unpublished*.

.. [Abild2007] \ F. Abild-Pedersen, J. Greeley, F. Studt, J. Rossmeisl, T.R. Munter, P.G. Moses, E. Skúlason, T. Bligaard, and J.K. Nørskov. "Scaling Properties of Adsorption Energies for Hydrogen-Containing Molecules on Transition-Metal Surfaces." *Phys. Rev. Lett.* **99(1)**, p. 4-7 (2007).

