.. _releaseNotes:

*************
Release Notes
*************

RMG-Py Version 2.1.9
====================
Date: May 1, 2018

- Cantherm:
    - Atom counts are no longer necessary in input files and are automatically determined from geometries
    - Custom atom energies can now be specified in input files
    - Removed atom energies for a few ambiguous model chemistries
    - Add atom energies for B3LYP/6-311+g(3df,2p)

- Changes:
    - Refactored molecule.parser and molecule.generator modules into molecule.converter and molecule.translator to improve code organization
    - SMILES generation now outputs canonical SMILES
    - Molecule.sortAtoms method restored for deterministic atom order
    - PDep reactions which match an existing library reaction are no longer added to the model

- Fixes:
    - Fix issue with reaction filter initiation when using seed mechanisms

RMG-database Version 2.1.9
==========================
Date: May 1, 2018

- Chlorine:
    - New Chlorinated_Hydrocarbons thermo library
    - Added group additivity values and long distance corrections for chlorinated species
    - Added chlorine groups and training reactions to H_Abstraction

- Additions:
    - New NOx2018 kinetics, thermo, and transport libraries
    - New N-S_interactions kinetics library
    - New SulfurHaynes thermo library
    - Added species to SOxNOx thermo library from quantum calculations

- Other changes:
    - Renamed NOx and SOx kinetics libraries to PrimaryNitrogenLibrary and PrimarySulfurLibrary
    - S2O2, SOO2, SO2O2, and N2SH were globally forbidden due to inability to optimize geometries

- Fixes:
    - Corrected some A-factor units in Nitrogen_Dean_and_Bozzelli kinetics library


RMG-Py Version 2.1.8
====================
Date: March 22, 2018

- New features:
    - Chlorine and iodine atom types have been added, bringing support for these elements to RMG-database
    - Forbidden structures now support Molecule and Species definitions in addition to Group definitions

- Changes:
    - Reaction pair generation will now fall back to generic method instead of raising an exception
    - Removed sensitivity.py script since it was effectively a duplicate of simulate.py
    - Thermo jobs in Cantherm now output a species dictionary
    - Fitted atom energy corrections added for B3LYP/6-31g**
    - Initial framework added for hydrogen bonding
    - Renamed molepro module and associated classes to molpro (MolPro) to match actual spelling of the program
    - Chemkin module is now cythonized to improve performance

- Fixes:
    - Allow delocalization of triradicals to prevent hysteresis in resonance structure generation
    - Fix reaction comment parsing issue with uncertainty analysis
    - Fix numerical issue causing a number of pressure dependent RMG jobs to crash
    - Template reactions from seed mechanisms are now loaded as library reactions if the original family is not loaded
    - Fix issues with degeneracy calculation for identical reactants

RMG-database Version 2.1.8
==========================
Date: March 22, 2018

- Changes:
    - Corrected name of JetSurf2.0 kinetics and thermo libraries to JetSurf1.0
    - Added actual JetSurf2.0 kinetics and thermo libraries
    - Updated thermo groups for near-aromatic radicals, including radical and polycyclic corrections


RMG-Py Version 2.1.7
====================
Date: February 12, 2018

- Charged atom types:
    - Atom types now have a charge attribute to cover a wider range of species
    - New atom types added for nitrogen and sulfur groups
    - Carbon and oxygen atom types renamed following new valence based naming scheme

- Ring perception:
    - Ring perception methods in the Graph class now use RingDecomposerLib
    - This includes the getSmallestSetOfSmallestRings methods and a newly added getRelevantCycles method
    - The set of relevant cycles is unique and generally more useful for chemical graphs
    - This also fixes inaccuracies with the original SSSR method

- Other changes:
    - Automatically load reaction libraries when using a seed mechanism
    - Default kinetics estimator has been changed to rate rules instead of group additivity
    - Kinetics families can now be set to be irreversible
    - Model enlargement now occurs after each reactor simulation rather than after all of them
    - Updated bond additivity corrections for CBS-QB3 in Cantherm

- Fixes:
    - Do not print SMILES when raising AtomTypeError to avoid further exceptions
    - Do not recalculate thermo if a species already has it
    - Fixes to parsing of family names in seed mechanisms


RMG-database Version 2.1.7
==========================
Date: February 12, 2018

- Charged atom types:
    - Update adjlists with new atom types across the entire database
    - Added sulfur groups to all relevant kinetics families
    - New thermo group additivity values for sulfur/oxygen species

- Additions:
    - Benzene bonds can now react in in R_Addition_MultipleBond
    - Many new training reactions and groups added in R_Addition_MultipleBond
    - New Singlet_Val6_to_triplet kinetics family
    - New Sulfur GlarborgBozzelli kinetics and thermo libraries
    - New Sulfur GlarborgMarshall kinetics and thermo libraries
    - New Sulfur GlarborgH2S kinetics and thermo libraries
    - New Sulfur GlarborgNS kinetics and thermo libraries
    - New NOx and NOx/LowT kinetics libraries
    - New SOx kinetics library
    - New BurcatNS thermo library
    - New SOxNOx thermo library
    - New 2+2_cycloaddition_CS kinetics family
    - New Cyclic_Thioether_Formation kinetics family
    - New Lai_Hexylbenzene kinetics and thermo libraries

- Changes:
    - 1,2-Birad_to_alkene family is now irreversible
    - OxygenSingTrip kinetics library removed (replaced by Singlet_Val6_to_triplet family)
    - Ozone is no longer forbidden

- Fixes:
    - Corrected adjlist for phenyl radical in JetSurf2.0 and USC-Mech-ii
    - Some singlet thermo groups relocated from radical.py to group.py


RMG-Py Version 2.1.6
====================
Date: December 21, 2017

- Model resurrection:
    - Automatically attempts to save simulation after encountering a DASPK error
    - Adds species and reactions in order to modify model dynamics and fix the error

- New features:
    - Add functionality to read RCCSD(T)-F12 energies from MolPro log files
    - Add liquidReactor support to flux diagram generation

- Other changes:
    - Removed rmgpy.rmg.model.Species class and merged functionality into main rmgpy.species.Species class
    - Refactored parsing of RMG-generated kinetics comments from Chemkin files and fixed related issues
    - Refactored framework for generating reactions to reduce code duplication
    - Resonance methods renamed from generateResonanceIsomers to generate_resonance_structures across all modules
    - Raise CpInf to Cphigh for entropy calculations to prevent invalid results

- Fixes:
    - Update sensitivity analysis to use ModelSettings and SimulatorSettings classes introduced in v2.1.5
    - Fixed generate_reactions methods in KineticsDatabase to be directly usable again
    - Fixed issues with aromaticity perception and generation of aromatic resonance structures

RMG-database Version 2.1.6
==========================
Date: December 21, 2017

- Additions:
    - New training reactions added for [NH2] related H_Abstractions
    - 14 new kinetics libraries related to aromatics formation (see RMG-database #222 for details)

- Other changes:
    - Removed some global forbidden groups which are no longer needed
    - Forbid CO and CS biradicals
    - Updated lone_electron_pair_bond family and removed from recommended list

- Fixes:
    - Fixed unit errors in some H_Abstraction and R_Addition_MultipleBond depositories


RMG-Py Version 2.1.5
====================
Date: October 18, 2017

- New bicyclic formula:
    - Estimates polycyclic corrections for unsaturated bicyclics by adjusting the correction for the saturated version
    - Can provide a decent estimate in many cases where there is not an exact match

- Other changes:
    - Refactored simulation algorithm to properly add multiple objects per iteration
    - Print equilibrium constant and reverse rate coefficient values when using Cantherm to calculate kinetics
    - Speed up degeneracy calculation by reducing unnecessary operations

- Fixes:
    - Loosen tolerance for bond order identification to account for floating point error
    - Fixed uncertainty analysis to allow floats as bond orders
    - Fixed some comment parsing issues in uncertainty analysis
    - Added product structure atom relabeling for families added in RMG-database v2.1.5
    - Fixed issue with automatic debugging of kinetics errors due to forbidden structures

RMG-database Version 2.1.5
==========================
Date: October 18, 2017

- Additions:
    - New thermo groups added for species relevant in cyclopentadiene and natural gas pyrolysis
    - Added C2H4+O_Klipp2017 kinetics library

- Fixes:
    - Prevent charged carbenes from reacting in Singlet_Carbene_Intra_Disproportionation
    - Updated H_Abstraction rates in ethylamine library and corresponding training reactions


RMG-Py Version 2.1.4
====================
Date: September 08, 2017

- Accelerator tools:
    - Dynamics criterion provides another method to expand the mechanism by adding reactions to the core
    - Surface algorithm enables better control of species movement to the core when using the dynamics criterion
    - Multiple sets of model parameters can now be specified in a input file to allow different stages of model generation
    - A species number termination criterion can now be set to limit model size
    - Multiple items can now be added per iteration to speed up model construction
    - New ModelSettings and SimulatorSettings classes for storing input parameters

- New features:
    - Kinetics libraries can now be automatically generated during RMG runs to be used as seeds for subsequent runs
    - Loading automatically generated seed mechanisms recreates the original template reaction objects to allow restarting runs from the seed mechanism
    - Carbene constraints can now be set in the species constraint block using maxSingletCarbenes and maxCarbeneRadicals
    - Chirality is now considered for determining symmetry numbers
    - Thermodynamic pruning has been added to allow removal of edge species with unfavorable free energy (beta)

- Other changes:
    - RMG-Py exception classes have been consolidated in the rmgpy.exceptions module
    - Species labels will now inherit the label from a matched thermo library entry
    - Sensitivity analysis is now available for LiquidReactor

- Fixes:
    - Fixed sensitivity analysis following changes to the simulate method
    - Add memory handling when generating collision matrix for pressure dependence
    - Improved error checking for MOPAC
    - Prevent infinite loops when retrieving thermo groups

- Known issues:
    - Seed mechanisms cannot be loaded if the database settings are different from the original ones used to generate the seed

RMG-database Version 2.1.4
==========================
Date: September 08, 2017

- New kinetics families for propargyl recombination route to benzene:
    - Singlet_Carbene_Intra_Disproportionation
    - Intra_5_membered_conjugated_C=C_C=C_addition
    - Intra_Diels_alder_monocyclic
    - Concerted_Intra_Diels_alder_monocyclic_1,2_shift
    - Intra_2+2_cycloaddition_Cd
    - Cyclopentadiene_scission
    - 6_membered_central_C-C_shift

- Renamed kinetics families:
    - Intra_Diels_Alder --> Intra_Retro_Diels_alder_bicyclic
    - H_shift_cyclopentadiene --> Intra_ene_reaction

- Other additions:
    - Klippenstein_Glarborg2016 kinetics and thermo libraries
    - Group additivity values added for singlet carbenes, which are no longer forbidden


RMG-Py Version 2.1.3
====================
Date: July 27, 2017

- Thermo central database:
    - Framework for tracking and submitting species to a central database have been added
    - Following species submission, the central database will queue and submit quantum chemistry jobs for thermochemistry calculation
    - This is an initial step towards self-improving thermochemistry prediction

- Rotor handling in Cantherm:
    - Free rotors can now be specified
    - Limit number of terms used when fitting hinder rotor scans
    - Fixed bug with ZPE calculation when using hindered rotors

- New reaction degeneracy algorithm:
    - Use atom ID's to distinguish degenerate reactions from duplicates due to other factors
    - Degeneracy calculation now operates across all families rather than within each separately
    - Multiple transition states are now identified based on template comparisons and kept as duplicate reactions

- Nodal distances:
    - Distances can now be assigned to trees in reaction families
    - This enables better rate averages with multiple trees
    - Fixed bug with finding the closest rate rule in the tree

- New features:
    - Added methods for automatically writing RMG-database files
    - New symmetry algorithm improves symmetry number calculations for resonant and cyclic species
    - Group additivity algorithm updated to apply new long distance corrections
    - Specific colliders can now be specified for pressure-dependent rates
    - Very short superminimal example added (hydrogen oxidation) for checking basic RMG operation
    - Cantera now outputs a Chemkin file which can be directly imported into Chemkin

- Fixes:
    - Fixed bug with negative activation energies when using Evans-Polanyi rates
    - Fixed walltime specification from command line when running RMG
    - Fixes and unit tests added for diffusionLimited module

- Known issues:
    - The multiple transition state algorithm can result in undesired duplicate reactions for reactants with multiple resonance structures

RMG-database Version 2.1.3
==========================
Date: July 27, 2017

- Long-distance interaction thermo corrections:
    - The gauche and int15 group files have been replaced by longDistanceInteraction_noncyclic
    - New corrections for cyclic ortho/meta/para interactions are now available in longDistanceInteraction_cyclic

- Changes:
    - Oa_R_Recombination family renamed to Birad_R_Recombination
    - More training reactions added for sulfur species in H_Abstraction
    - RMG-database tests have been moved to RMG-Py


RMG-Py Version 2.1.2
====================
Date: May 18, 2017

- Improvements:
    - New nitrogen atom types
    - Kinetics libraries can now be specified as a list of strings in the input file
    - New script to generate output HTML locally: generateChemkinHTML.py
    - New kekulization module replaces RDKit for generating Kekule structures
    - Benzene bonds can now be reacted in reaction families
    - Removed cantherm.geometry module due to redundancy with statmech.conformer

- Fixes:
    - Reaction direction is now more deterministic after accounting for floating point error
    - Multiple bugs with resonance structure generation for aromatics have been addressed


RMG-database Version 2.1.2
==========================
Date: May 18, 2017

- Nitrogen improvements:
    - Added ethylamine kinetics library
    - Updated group additivity values for nitrogen species
    - Added rate rules and training reactions for nitrogen species

- Additions:
    - New CO_Disproportionation family
    - Added CurranPentane kinetics and thermo libraries

- Fixes:
    - Corrected some rates in FFCM1(-) to use MultiArrhenius kinetics
    - Corrected a few adjlists in FFCM1(-)


RMG-Py Version 2.1.1
====================
Date: April 07, 2017

- Uncertainty analysis:
    - Local and global uncertainty analysis now available for RMG-generated models
    - Global uncertainty analysis uses MIT Uncertainty Quantification library, currently only supported on Linux systems
    - Examples for each module are available in localUncertainty.ipynb and globalUncertainty.ipynb

- Fixes:
    - Clar structure generation no longer intercepts signals
    - Fixes to SMILES generation
    - Fix default spin state of [CH]

RMG-database Version 2.1.1
==========================
Date: April 07, 2017

- Additions:
    - More species added to FFCM1(-) thermo library

- Changes:
    - Improved handling of excited species in FFCM1(-) kinetics library
    - Replaced Klippenstein H2O2 kinetics and thermo libraries with BurkeH2O2inN2 and BurkeH2O2inArHe

- Fixes:
    - Corrected adjlists for some species in JetSurf2.0 kinetics and thermo libraries (also renamed from JetSurf0.2)
    - Correct multiplicities for [C] and [CH] in multiple libraries ([C] from 5 to 3, [CH] from 4 to 2)


RMG-Py Version 2.1.0
====================
Date: March 07, 2017

- Clar structure generation
    - optimizes the aromatic isomer representations in RMG
    - lays the foundations for future development of poly-aromatic kinetics reaction families

- Flux pathway analysis
    - introduces an ipython notebook for post-generatation pathway analysis (``ipython.mechanism_analyzer.ipynb``)
    - visualizes reactions and provides flux statistics in a more transparent way

- Cantera mechanism
    - automatically writes cantera version of RMG-generated mechanism at the end of RMG jobs

- Fixes bugs
    - upgrades ``pruning`` to fix new memory leaks introduced by recent functionalities
    - fixes the bug of duplicated species creation caused by ``getThermoData`` removing isomers unexpectedly
    - fixes restart file generation and parsing problems and users can choose restart mode again
    - upgrades bicyclic decomposition method such that more deterministic behaviors are ensured
    - change bond order type to float from string to improve RMG's symmetry calculation for species with multiple resonance structures

RMG-database Version 2.1.0
==========================
Date: March 07, 2017

- Several new kinetics libraries added
    - FFCM-1
    - JetSurF 0.2
    - Chernov_aromatic_only
    - Narayanaswamy_aromatic_only
    - 1989_Stewart_2CH3_to_C2H5_H
    - 2005_Senosiain_OH_C2H2
    - 2006_Joshi_OH_CO
    - C6H5_C4H4_Mebel
    - c-C5H5_CH3_Sharma

- Several new thermochemistry libraries added
    - FFCM-1
    - JetSurF 0.2
    - Chernov_aromatic_only
    - Narayanaswamy_aromatic_only

- Improved kinetics tree accessibility
    - adds database tests ensuring groups in the tree to be accessible
    - improves definitions of group structures in the kinetics trees to ensure accessibility

- New oxygenates thermo groups are added based Paraskeva et al.

- Improved database tools
    - ``convertKineticsLibraryToTrainingReactions.ipynb`` now can visualize groups of matched rate rules that training reactions hit 
    - ``exportKineticsLibrarytoChemkin.py`` and ``importChemkinLibrary.py`` add more logging information on reaction sources


RMG-Py Version 2.0.0
====================
Date: September 16, 2016

This release includes several milestones of RMG project:

- Parallelization finally introduced in RMG:
    - Generates reactions during ``enlarge`` step in parallel fashion (``rmgpy.rmg.react``)
    - Enables concurrent computing for QMTP thermochemistry calculations (``rmgpy.thermo.thermoengine``)
    - Instructions of running RMG parallel mode can be found `here for SLURM scheduler <https://github.com/ReactionMechanismGenerator/RMG-Py/wiki/Running-RMG-in-parallel-with-a-SLURM-scheduler>`_ and `here for SGE scheduler <https://github.com/ReactionMechanismGenerator/RMG-Py/wiki/Running-RMG-in-parallel-with-a-SGE-scheduler>`_.

- Polycyclic thermochemistry estimation improved:
    - Extends group additivity method for polycyclics and estimates polycyclics of any large sizes by a heuristic method (bicyclics decomposition)

- New tree averaging for kinetics:
    - Fixes previous issue of imcomplete generation of cross-level rate rules
    - Implements Euclidean distance algorithm for the selection of the best rate rules to use in ``estimateKinetics``
    - Streamlines storage of kinetics comments for averaged rules, which can be analyzed by ``extractSourceFromComments``

- Database entry accessibility tests: 
    - Adds entry accessibility tests for future entries (``testing.databaseTest``)

- Fixes bugs
    - fluxdiagram generation is now fixed, one can use it to generate short video of fluxdigram evolution
    - mac environment yml file is introduced to make sure smooth RMG-Py installation and jobs on mac
    - fixes failure of ``checkForExistingSpecies`` for polyaromatics species
    - fixes execution failure when both pruning and pDep are turned on
    - fixes pDep irreversible reactions
    - fixes issue of valency of ``Cbf`` atom by dynamic benzene bond order assignment


RMG-database Version 2.0.0
==========================
Date: September 16, 2016

In conjunction with the release of RMG-Py v2.0.0, an updated package for the RMG-database has also been released.
This release brings some new additions and fixes:

- Polycyclic thermochemistry estimation improved:
    - polycyclic database reorganized and more entries added in systematic way (``input.thermo.groups.polycyclic``)

- Database entry accessibility tests:
    - Fixes existing inaccessible entries in solvation/statmech/thermo of RMG-database 


RMG-Py Version 1.0.4
====================
Date: March 28, 2016

- Cantera support in RMG (``rmgpy.tools.canteraModel``):
    - Provides functions to help simulate RMG models using Cantera.
    - Has capability to generate cantera conditions and convert CHEMKIN files to cantera models, or use RMG to directly convert species and reactions objects to Cantera objects.
    - Demonstrative example found in ``ipython/canteraSimulation.ipynb``

- Module for regression testing of models generated by RMG (``rmgpy.tools.observableRegression``):
    - Helps identify differences between two versions of models generated by RMG, using the "observables" that the user cares about.

- Automatic plotting of simulations and sensitivities when generating models (``rmgpy.tools.plot``):
    - Contains plotting classes useful for plotting simulations, sensitivities, and other data
    - Automatic plotting of simulations in the job's ``solver`` folder when ``saveSimulationProfiles`` is set to ``True`` in the input file. 
    - Sensitivities for top 10 most sensitivie reactions and thermo now plotted automatically and stored in the ``solver`` folder.

- Improved thermochemistry estimation (mostly for cyclics and polycyclics)
    - Add rank as an additional attribute in thermo database entries to determine trustworthiness

- Bug fixes:
    - Training reactions now load successfully regardless of ``generateSpeciesConstraints`` parameters
    - Transport data is now saved correctly to CHEMKIN ``tran.dat`` file and also imports successfully
    - Fixes appending of reactions to CHEMKIN file when reaction libraries are desired to be appended to output
    - Fixes writing of csv files for simulation and sensitivity results in Windows
    - Fixes ``Reaction.draw()`` function to draw the entire reaction rather than a single species


RMG-Py Version 1.0.3
====================
Date: February 4, 2016

This mini release contains the following updates:

- Pdep convergence issues in RMG-Py v1.0.2 are now fixed.
- RMG-database version information and anaconda binary version information is now recorded in RMG log file.


RMG-Py Version 1.0.2
====================
Date: January 29, 2016

This new release adds several new features and bug fixes. 

- Windows users can rejoice: RMG is now available in binary format on the Anaconda platform.  Building by source is also
  much easier now through the Anaconda managed python environment for dependencies. See the updated :ref:`Installation Page<installation>`
  for more details
- Reaction filtering for speeding up model generation has now been added.  It has been shown to speed up model convergence by
  7-10x.  See more details about how to use it in your RMG job :ref:`here <filterReactions>`.  Learn more about the theory 
  and algorithm on the :ref:`Rate-based Model Enlarging Algorithm <ratebasedmodelenlarger>` page.
- The RMG :ref:`native scripts <modules>` are now organized under the ``rmgpy.tools`` submodule for
  developer ease and better extensibility in external scripts.
- InChI conversion is now more robust for singlets and triplets, 
  and augmented InChIs and InChI keys are now possible with new radical electron, lone pair, and multiplicity flags.  
- Output HTML for visualizing models are now cleaned up and also more functional, including features to display thermo comments,
  display enthalpy, entropy, and free energy of reaction, as well as filter reactions by species.  You can use this new visualization format
  either by running a job in RMG v1.0.2 or revisualizing your CHEMKIN file and species dictionary using
  the `visualization web tool <http://rmg.mit.edu/simulate/chemkin>`_.
  
  
  
RMG-database Version 1.0.2
==========================
Date: January 29, 2016

In conjunction with the release of RMG-Py v1.0.2, an updated package for the RMG-database has also been released.
This release brings some new additions and fixes:

- New group additivity values for oxitene, oxerene, oexpane, and furan ring groups
- Improvements to sulfur chemistry:
    - Restructuring of radical trees in the kinetics families ``SubstitutionS`` and ``intra_substitutionCS_cyclization``
    - A reaction library for di-tert-butyl sulfide
- Improvements for the ``R_Addition_Multiple_Bond`` kinetics family through new rate rules
  for the addition of allyl radical to double bonds in ethene, propene, and butene-like
  compounds, based on CBS-QB3 estimates from K. Wang, S.M. Villano, A.M. Dean, 
  "Reactions of allylic radicals that impact molecular weight growth kinetics", *PCCP*,
  6255-6273 (2015).
- Several new thermodynamic and kinetics libraries for molecules associated with the
  pyrolysis of cyclopentadiene in the presence of ethene, based off of calculations from
  the paper A.G. Vandeputte, S.S. Merchant, M.R. Djokic, K.M. Van Geem, 
  G.B. Marin, W. H. Green, "Detailed study of cyclopentadiene pyrolysis in the 
  presence of ethene: realistic pathways from C5H5 to naphthalene" (2016)