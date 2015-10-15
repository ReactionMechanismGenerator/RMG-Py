.. _liquids:

********************
Liquid Phase Systems
********************

To simulate liquids in RMG requires a module in your input file for liquid-phase::


    solvation(
        solvent='octane'
    )

Your reaction system will also be different (liquidReactor rather than simpleReactor)::

    liquidReactor(
        temperature=(500,'K'),
        initialConcentrations={
            "octane": (6.154e-3,'mol/cm^3'),
            "oxygen": (4.953e-6,'mol/cm^3')
        },
        terminationTime=(5,'s'),
        sensitivity=['octane','oxygen'],
        sensitivityThreshold=0.001,

    )


For sensitivity analysis, RMG-Py must be compiled with the DASPK solver. 
(See :ref:`Compiling RMG-Py with Sensitivity Analysis  <compile_sensitivity>` for more details.)
Like for the simpleReactor, the sensitivity and sensitivityThrehold are optional arguments for when the
user would like to conduct sensitivity analysis with respect to the reaction rate
coefficients for the list of species given for ``sensitivity``.  

Sensitivity analysis is conducted for the list of species given for ``sensitivity`` argument in the input file.  
The normalized concentration sensitivities with respect to the reaction rate coefficients dln(C_i)/dln(k_j) are saved to a csv file 
with the file name ``sensitivity_1_SPC_1.csv`` with the first index value indicating the reactor system and the second naming the index of the species 
the sensitivity analysis is conducted for.  Sensitivities to thermo of individual species is also saved as semi normalized sensitivities
dln(C_i)/d(G_j) where the units are given in 1/(kcal mol-1). The sensitivityThreshold is set to some value so that only
sensitivities for dln(C_i)/dln(k_j) > sensitivityThreshold  or dlnC_i/d(G_j) > sensitivityThreshold are saved to this file.  

Note that in the RMG job, after the model has been generated to completion, sensitivity analysis will be conducted
in one final simulation (sensitivity is not performed in intermediate iterations of the job).

Equation of state
=================
Specifying a liquidReactor will have two effects:

#. disable the ideal gas law renormalization and instead rely on the concentrations you specified in the input file to initialize the system.
#. prevent the volume from changing when there is a net stoichiometry change due to a chemical reaction (A = B + C).


Solvation thermochemistry
=========================

The next correction for liquids is solvation effects on the thermochemistry. By specifying a solvent in the input file, we load the solvent parameters to use.

The free energy change associated with the process of transferring a
molecule from the gas phase to the solvent phase is defined as the free
energy of solvation (ΔG). Many different methods have been developed for
computing solvation energies among which continuum dielectric and force
field based methods are popular. Not all of these methods are easy to
automate, and many are not robust i.e. they either fail or give
unreasonable results for certain solute-solvent pairs. CPU time and
memory (RAM) requirements are also important considerations. A fairly
accurate and fast method for computing ΔG, which is used in RMG, is the
LSER approach described below.

.. _useofthermolibrariesliquid:

Use of thermo libraries in liquid phase system
----------------------------------------------

As it is for gas phase simulation, thermo libraries listed in the input files are checked first to find thermo for a given species and return the first match.
As it exists two types of thermo libraries, (more details on :ref:`thermo libraries <thermoDatabase>`),
thermo of species matching a library in a liquid phase simulation is obtained following those two cases:

If library is a "liquid thermo library", thermo data are directly used without applying solvation on it.

If library is a "gas thermo library", thermo data are extracted and then corrections are applied on it using the :ref:`LSER method <lserToEstimateThermo>`
for this specific species-solvent system. 

.. note::
	Gas phase libraries can be declared first, liquid thermo libraries will still be tested first but the order will be respected if several liquid libraries are provided. 

.. _lserToEstimateThermo:

Use of Abraham LSER to estimate thermochemistry
-----------------------------------------------

The Abraham LSER provides an estimate of the the partition coefficient (more specifically, the 
log (base 10) of the partition coefficient) of a solute between the vapor phase and a particular solvent 
(`K`\ :sub:`vs`\ ) (also known as gas-solvent partition coefficient) at 298 K:

.. math:: \log K_{vs} = c + eE + sS + aA + bB + lL
	:label: AbModelEqn

The Abraham model is used in RMG to estimate ΔG which is related to the `K`\ :sub:`vs`\  of a solute according to the following expression:

.. math:: ΔG = -RT \ln K_{vs} \\
	= -2.303RT \log K_{vs}
	:label: partition

The variables in the Abraham model represent solute (`E, S, A, B, V, L`) and solvent descriptors (`c, e, s, a, b, v, l`) 
for different interactions. The `sS` term is attributed to electrostatic interactions between the 
solute and the solvent (dipole-dipole interactions related to solvent dipolarity and the dipole-induced 
dipole interactions related to the polarizability of the solvent) [Vitha2006]_, [Abraham1999]_, [Jalan2010]_. The 
`lL` term accounts for the contribution from cavity formation and dispersion (dispersion interactions are 
known to scale with solute volume [Vitha2006]_, [Abraham1999]_. The `eE` term, like the `sS` term, 
accounts for residual contributions from dipolarity/polarizability related interactions for solutes 
whose blend of dipolarity/polarizability differs from that implicitly built into the `S` parameter [Vitha2006]_, [Abraham1999]_, [Jalan2010]_. 
The `aA` and `bB` terms account for the contribution of hydrogen bonding between the solute and 
the surrounding solvent molecules. H-bonding interactions require two terms as the solute (or solvent) 
can act as acceptor (donor) and vice versa. The descriptor `A` is a measure of the solute's ability 
to donate a hydrogen bond (acidity) and the solvent descriptor `a` is a measure of the solvent's ability 
to accept a hydrogen bond. A similar explanation applies to the `bB` term [Vitha2006]_, [Abraham1999]_, [Poole2009].


The solvent descriptors (`c, e, s, a, b, l`) are largely treated as regressed empirical coefficients. Parameters are provided in RMG's database for the following solvents:

#. acetonitrile
#. benzene
#. butanol
#. carbontet
#. chloroform
#. cyclohexane
#. decane
#. dibutylether
#. dichloroethane
#. dimethylformamide
#. dimethylsulfoxide
#. dodecane
#. ethanol
#. ethylacetate
#. heptane
#. hexadecane
#. hexane
#. isooctane
#. nonane
#. octane
#. octanol
#. pentane
#. toluene
#. undecane
#. water

Group additivity method for solute descriptor estimation
--------------------------------------------------------

Group additivity is a convenient way of estimating the thermochemistry for thousands of species sampled 
in a typical mechanism generation job. Use of the Abraham Model in RMG requires a similar approach 
to estimate the solute descriptors (`A, B, E, L,` and `S`). Platts et al. ([Platts1999]_) proposed such a scheme employing a set of 81 molecular fragments for estimating `B, E, L, V` and `S` and another set of 51 fragments for the estimation of `A`. Only those fragments containing C, H and O are implemented in order to match RMG's existing capabilities. The value of a given descriptor for a molecule is obtained by summing the contributions from each fragment found in the molecule and the intercept associated with that descriptor.

Mintz model for enthalpy of solvation
-------------------------------------

For estimating ΔG at temperatures other than 298 K, the enthalpy change associated with solvation, ΔH must be calculated separately and, along with ΔS, assumed to be independent of temperature. Recently, Mintz et al. ([Mintz2007]_, [Mintz2007a]_, [Mintz2007b]_, [Mintz2007c]_, [Mintz2007d]_, [Mintz2008]_, [Mintz2008a]_, [Mintz2009]_) have developed linear correlations similar to the Abraham model for estimating ΔH:

.. math:: ΔH(298 K) = c' + a'A+ b'B+ e'E+ s'S+ l'L
	:label: mintz

where `A, B, E, S` and `L` are the same solute descriptors used in the Abraham model for the estimation of ΔG. The lowercase coefficients `c', a', b', e', s'` and `l'` depend only on the solvent and were obtained by fitting to experimental data. In RMG, this equation is implemented and together with ΔG(298 K) can be used to find ΔS(298 K). From this data, ΔG at other temperatures is found by extrapolation.


Diffusion-limited kinetics
==========================
The next correction for liquid-phase reactions is to ensure that bimolecular reactions do not exceed their diffusion limits. The theory behind diffusive limits in solution phase reactions is well established ([Rice1985]_) and the effective rate constant of a bimolecular reaction is given as:

.. math::   k_{\textrm{eff}} = \frac {4\pi R\mathcal{D} k_{\textrm{int}}}{4\pi R\mathcal{D} + k_{\textrm{int}}}
   :label: diffusive_limit

where `k`\ :sub:`int` is the intrinsic reaction rate, `R` is the sum of radii of the reactants and 
`D` is the sum of the diffusivities of the reacting species. RMG uses the McGowan method for estimating 
radii, and diffusivities are estimated with the Stokes-Einstein equation using experimental solvent 
viscosities (`\eta` (T)).  In a unimolecular to bimolecular reaction, for example, the forward rate 
constant (`k`\ :sub:`f`\ ) can be slowed down if the reverse rate (`k`\ :sub:`r, eff`\ ) is diffusion 
limited since the equilibrium constant (`K`\ :sub:`eq`\ ) is not affected by diffusion limitations. In cases 
where both the forward and the reverse reaction rates are bimolecular, both diffusive limits are 
estimated and RMG uses the direction with the larger magnitude.

The viscosity of the solvent is calculated Pa.s using the solvent specified in the command line 
and a correlation for the viscosity using parameters `A, B, C, D, E`:

.. math:: \ln \eta = A + \frac{B}{T} + C\log T + DT^E
    :label: viscosity
       
To build accurate models of liquid phase chemical reactions you will also want to modify your kinetics libraries or correct gas-phase rates for intrinsic barrier solvation corrections (coming soon).

.. _exampleLiquidPhase:

Example liquid-phase input file
===============================
This is an example of an input file for a liquid-phase system::

    # Data sources
    database(
        thermoLibraries = ['primaryThermoLibrary'],
        reactionLibraries = [],
        seedMechanisms = [],
        kineticsDepositories = ['training'],
        kineticsFamilies = 'default',
        kineticsEstimator = 'rate rules',
    )

    # List of species
    species(
        label='octane',
        reactive=True,
        structure=SMILES("C(CCCCC)CC"),
    )

    species(
        label='oxygen',
        reactive=True,
        structure=SMILES("[O][O]"),
    )

    # Reaction systems
    liquidReactor(
        temperature=(500,'K'),
        initialConcentrations={
            "octane": (6.154e-3,'mol/cm^3'),
            "oxygen": (4.953e-6,'mol/cm^3')
        },
        terminationTime=(5,'s'),
    )

    solvation(
        solvent='octane'
    )

    simulator(
        atol=1e-16,
        rtol=1e-8,
    )

    model(
        toleranceKeepInEdge=1E-9,
        toleranceMoveToCore=0.001,
        toleranceInterruptSimulation=0.1,
        maximumEdgeSpecies=100000
    )

    options(
        units='si',
        saveRestartPeriod=None,
        generateOutputHTML=False,
        generatePlots=False,
        saveSimulationProfiles=True,
    )

.. [Vitha2006] \ M. Vitha and P.W. Carr. "The chemical interpretation and practice of linear solvation energy relationships in chromatography." *J. Chromatogr. A.* **1126(1-2)**, p. 143-194 (2006).

.. [Abraham1999] \ M.H. Abraham et al. "Correlation and estimation of gas-chloroform and water-chloroformpartition coefficients by a linear free energy relationship method." *J. Pharm. Sci.* **88(7)**, p. 670-679 (1999).

.. [Jalan2010] \ A. Jalan et al. "Predicting solvation energies for kinetic modeling." *Annu. Rep.Prog. Chem., Sect. C* **106**, p. 211-258 (2010).

.. [Poole2009] \ C.F. Poole et al. "Determination of solute descriptors by chromatographic methods." *Anal. Chim. Acta* **652(1-2)** p. 32-53 (2009).

.. [Platts1999] \ J. Platts and D. Butina. "Estimation of molecular linear free energy relation descriptorsusing a group contribution approach." *J. Chem. Inf. Comput. Sci.* **39**, p. 835-845 (1999).

.. [Mintz2007] \ C. Mintz et al. "Enthalpy of solvation correlations for gaseous solutes dissolved inwater and in 1-octanol based on the Abraham model." *J. Chem. Inf. Model.* **47(1)**, p. 115-121 (2007).

.. [Mintz2007a] \ C. Mintz et al. "Enthalpy of solvation corrections for gaseous solutes dissolved in benzene and in alkane solvents based on the Abraham model." *QSAR Comb. Sci.* **26(8)**, p. 881-888 (2007).

.. [Mintz2007b] \ C. Mintz et al. "Enthalpy of solvation correlations for gaseous solutes dissolved in toluene and carbon tetrachloride based on the Abraham model." *J. Sol. Chem.* **36(8)**, p. 947-966 (2007).

.. [Mintz2007c] \ C. Mintz et al. "Enthalpy of solvation correlations for gaseous solutes dissolved indimethyl sulfoxide and propylene carbonate based on the Abraham model." *Thermochim. Acta* **459(1-2)**, p, 17-25 (2007).

.. [Mintz2007d] \ C. Mintz et al. "Enthalpy of solvation correlations for gaseous solutes dissolved inchloroform and 1,2-dichloroethane based on the Abraham model." *Fluid Phase Equilib.* **258(2)**, p. 191-198 (2007).

.. [Mintz2008] \ C. Mintz et al. "Enthalpy of solvation correlations for gaseous solutes dissolved inlinear alkanes (C5-C16) based on the Abraham model." *QSAR Comb. Sci.* **27(2)**, p. 179-186 (2008).

.. [Mintz2008a] \ C. Mintz et al. "Enthalpy of solvation correlations for gaseous solutes dissolved inalcohol solvents based on the Abraham model." *QSAR Comb. Sci.* **27(5)**, p. 627-635 (2008).

.. [Mintz2009] \ C. Mintz et al. "Enthalpy of solvation correlations for organic solutes and gasesdissolved in acetonitrile and acetone." *Thermochim. Acta* **484(1-2)**, p. 65-69 (2009).

.. [Rice1985] \ S.A. Rice. "Diffusion-limited reactions". In *Comprehensive Chemical Kinetics*, EditorsC.H. Bamford, C.F.H. Tipper and R.G. Compton. **25**, (1985).
