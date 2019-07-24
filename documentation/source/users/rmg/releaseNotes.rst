.. _releaseNotes:

*************
Release Notes
*************

RMG-Py Version 2.4.1
====================
Date: July 23, 2019

- Bugfixes
    - Improve error handling in NASA as_dict method #1630
    - Fixes to Fluorine atomtypes #1656
    - Fix pressure dependent network generation #1658
    - Add support for reversing PDepArrhenius with MultiArrhenius rates #1659

- Arkane
    - Implement ZPE scaling factor #1619
    - Refactor infrastructure for bond additivity corrections #1605
    - Add frequency scale factors for wb97xd/def2tzvp and apfd/def2tzvpp #1653
    - Fix frequency scale factors in example files #1657
    - Get element counts from conformers #1651

- Miscellaneous
    - Update conda environment files #1623, #1644
    - Output RMS (Reaction Mechanism Simulator) format mechanism files #1629
    - Properly clean up files after running tests #1645
    - Documentation fixes #1650
    - Improve as_dict and make_object by making them recursive #1643


RMG-Py Version 2.4.0
====================
Date: June 14, 2019

- Heterogeneous catalysis!
    - RMG-cat fork has been merged #1573
        - Introduce SurfaceReactor
        - Thermo estimation for adsorbed species
        - Surface reaction generation and kinetics estimation
    - Introduce Van der Waals bonds (order 0) and quadruple bonds (order 4) #1542
- Arkane
    - Automatically detect rotor symmetry #1526
    - Introduce new YAML files for storing and loading species statmech data #1402, #1551
    - Don't create species dictionary file if there are no structures #1528
    - Improvements to network explorer tool #1545
    - Improved class inheritance for quantum log file classes #1571
    - Automatic determination of optical isomers and symmetry using ``symmetry`` package #1571
    - Parse CCSD(T) energies from Molpro output #1592
    - Automatically determine molecule linearity #1601
    - Determine frequency scaling factor based on geom/freq method rather than sp method #1612
    - Improve logging related to energy barriers #1575
    - Ensure that translational mode is calculated for atoms #1620
- Miscellaneous features
    - New ``enumerate_bonds`` method of Molecule to generate dictionary of bond types #1525
    - Introduce ``RMGObject`` parent class to support YAML dumping and loading #1402, #1540
    - Add support for fluorine atomtypes #1543
    - Introduce ``ArrheniusBM`` class for Blower-Masel kinetics #1461
    - Allow defining and using co-solvents for solvent libraries #1558
    - Introduce ``strict`` option to perform isomorphism between species/molecules while ignoring electrons and bond orders #1329
    - Molecule and Species objects can be instantiated by providing ``SMILES`` or ``InChI`` argument directly, and the identifiers can be accessed via the ``SMILES`` and ``InChI`` attributes #1329
    - Parallelization has been completely refactored using Python multiprocessing module in replacement of scoop, currently supports parallel reaction generation and QMTP #1459
    - Improvements to usability of uncertainty analysis functionality #1593
- Bug fixes
    - Various fixes for supporting mono-atomic molecules in Arkane #1513, #1521
    - Ensure ``keras_backend`` is set consistently #1535
    - Fix handling of disconnected graphs in VF2 isomorphism algorithm #1538
    - Ignore hydrogen bonds when converting to RDKit molecule #1552
    - Other miscellaneous bugs #1546, #1556, #1593, #1600, #1622
- Backward incompatible changes
    - Hydrogen bonds are now order 0.1 (instead of 0) #1542
- New dependencies
    - pyyaml (required) #1402
    - scikit-learn (required) #1461
    - textgenrnn (optional) #1573
- Other
    - Windows binaries are no longer officially supported. The new recommended way to use RMG on Windows computers is via a virtual machine or through the Linux subsystem. See documentation for updated installation instructions. #1531, #1534
    - Documentation updates #1544, #1567
    - Logging/exception improvements #1538, #1562
    - PEP-8 improvements #1566, #1592, #1596
    - Solver output files (png/csv) now report moles instead of mole fractions #1542
    - Replace global RMGDatabase object if the database is reloaded #1565
    - Print ML generated quote upon completion of RMG jobs #1573
    - Infrastructure for automatically generated reaction rate trees #1461
    - Testing related changes #1597, #1599
    - Updates to example Jupyter notebooks #1541, #1593

RMG-database Version 2.4.0
==========================
Date: June 14, 2019

- Heterogeneous catalysis!
    - RMG-cat fork has been merged #309
    - New kinetics families
        - Surface_Adsorption_Single
        - Surface_Adsorption_vdW
        - Surface_Adsorption_Dissociative
        - Surface_Dissociation
        - Surface_Abstraction
        - Surface_Adsorption_Double
        - Surface_Dissociation_vdW
        - Surface_Adsorption_Bidentate
        - Surface_Bidentate_Dissociation
        - Surface_Recombination (deprecated, use Surface_Dissociation instead)
    - New thermo group types
        - adsorptionNi
        - adsorptionPt
    - New thermo libraries
        - surfaceThermoNi
        - surfaceThermoPt
- New kinetics families
    - 1,2_NH3_elimination #326
    - 1,3_NH3_elimination #326
- New kinetics libraries
    - HydrazinePDep #326
- New transport libraries
    - OneDMinN2 #326
- Kinetics training reaction additions
    - 1,2_shiftC #306
    - Intra_R_Add_Endocyclic #306, #258
    - Intra_R_Add_Exocyclic #306, #258, #331
    - Intra_ene_reaction #306
    - R_Addition_COm #306
    - R_Addition_MultipleBond #306, #258
    - R_Recombination #306,  #326
    - Intra_H_migration #306
    - H_Abstraction #326
- Kinetics library additions
    - primaryNitrogenLibrary #326
    - Lai_Hexylbenzene #258
- Thermo library additions
    - CBS_QB3_1dHR, thermo_DFT_CCSDTF12_BAC #319
    - primaryNS #326
    - Lai_Hexylbenzene #258
- Thermo group additions
    - ring, polycyclic, radical #258
- Changes
    - [adjlist] kinetics/libraries/Klippenstein_Glarborg2016 #308
    - [labels] thermo/libraries/CBS_QB3_1dHR, Narayanaswamy #306
    - [units] kinetics/libraries/Sulfur/GlarborgMarhsall, Nitrogen_Dean_and_Bozzelli, primaryNitrogenLibrary, primarySulfurLibrary #311
    - [units] R_Addition_MultipleBond/training, R_Recombination/training #312
    - [adjlist] kinetics/libraries/GRI-Mech3.0-N #313
    - [adjlist] thermo/libraries/GRI-Mech3.0-N, GRI-Mech3.0 #313
    - [rates] Disproportionation/training, R_Addition_MultipleBond/training #326
    - [labels] kinetics/libraries/NOx2018 #326
    - [labels, attributes] kinetics/libraries/Nitrogen_Dean_and_Bozelli #326
    - [labels] kinetics/librariesNitrogen_Glarbog_Gimenez_et_al, Nitrogen_Glarborg_Zhang_et_al  #326
    - [labels, adjlist] thermo/libraries/BurcatNS #326
    - [labels] thermo/libraries/NOx2018, NitrogenCurran #326
    - [labels] transport/libraries/NOx2018 #326
    - [adjlist] Intra_R_Add_Endocyclic/training #332
    - [value] thermo/groups/ring/12dioxetane #327
    - [adjlist] thermo/libraries/GRI-Mech3.0 #336
    - [value] thermo/libraries/primaryThermoLibrary #338


RMG-Py Version 2.3.0
====================
Date: Dec 20, 2018

- Arkane (formerly CanTherm):
    - CanTherm had been renamed to Arkane (Automated Reaction Kinetics And Network Exploration)
    - New network exploration functionality using RMG-database
    - Support for all elements has been added for reading quantum output files
    - New supporting information output file with rotational constants and frequencies
    - Known thermo and kinetics can be provided in addition to quantum information
    - Improve general user experience and error handling

- New machine learning thermo estimator
    - Estimate species thermochemistry using a graph convolutional neural network
    - Estimator trained on quantum calculations at B3LYP and CCSD(T)-F12 levels
    - Currently supports C/H/O/N, with an emphasis on cyclic molecules

- Resonance:
    - New pathways added for lone-pair multiple-bond resonance, replacing
      two pathways which were more specific
    - New pathways added for aryne resonance
    - Aromatic resonance pathways simplified and refactored to use filtration
    - Kekule structures are now considered unreactive structures

- Miscellaneous changes:
    - Isotope support added for reading and writing InChI strings
    - New branching algorithm for picking up feedback loops implemented (beta)
    - Global forbidden structure checking is now only done for core species for
      efficiency, which may lead to forbidden species existing in the edge
    - Minor improvements to symmetry algorithm to fix a few incorrect cases

- Bug fixes:
    - Fixed issue where react flags were being reset when filterReactions was
      used with multiple reactors, resulting in no reactions generated
    - File paths for collision violators log changed to output directory
    - Fixed bug in local uncertainty introduced by ranged reactor changes
    - Fixed bug with diffusion limitation calculations for multi-molecular reactions
    - Various other minor fixes

RMG-database Version 2.3.0
==========================
Date: Dec 20, 2018

- Kinetics rules to training reactions
    - All kinetics rules have been converted into training reactions by converting
      each group to the smallest molecule that matches it
    - Training reactions are preferred over rules because they correspond to a
      specific reaction and are therefore easier to update
    - This conversion is in anticipation of upcoming changes to trees in kinetics families

- Additions:
    - R_Addition_MultipleBond training reactions
    - intra_NO2_ONO_conversion training reactions
    - SABIC_aromatics thermo library (CBS-QB3, RRHO)
    - McGowan volumes for noble gases
    - More entries added to Lai_Hexylbenzene libraries
    - Architecture and weights for neural network thermo estimator


RMG-Py Version 2.2.1
====================
Date July 23, 2018

This release is minor patch which fixes a number of issues discovered after 2.2.0.

- Collision limit checking:
    - RMG will now output a list of collision limit violations for the generated model

- Fixes:
    - Ambiguous chemical formulas in SMILES lookup leading to incorrect SMILES generation
    - Fixed issue with reading geometries from QChem output files
    - React flags for reaction filter were not properly updated on each iteration
    - Fixed issue with inconsistent symmetry number calculation


RMG-Py Version 2.2.0
====================
Date: July 5, 2018

- New features:
    - New ring membership attribute added to atoms. Can be specified in group adjacency lists in order to enforce
      ring membership of atoms during subgraph matching.
    - Reactors now support specification of T, P, X ranges. Different conditions are sampled on each iteration to
      optimally capture the full parameter space.
    - New termination type! Termination rate ratio stops the simulation when the characteristic rate falls to the
      specified fraction of the maximum characteristic rate. Currently not recommended for systems with two-stage ignition.
    - New resonance transitions implemented for species with lone pairs (particularly N and S containing species).
      A filtration algorithm was also added to select only the most representative structures.
    - Formal support for trimolecular reaction families.
    - New isotopes module allows post-processing of RMG mechanisms to generate a mechanism with isotopic labeling.

- Changes:
    - Library reactions can now be integrated into RMG pdep networks if the new elementary_high_p attribute is True
    - Library reactions may be duplicated by pdep reactions if the new allow_pdep_route attribute is True
    - Jupyter notebook for adding new training reactions has been revamped and is now located at ipython/kinetics_library_to_training.ipynb
    - Syntax for recommended families has changed to set notation instead of dictionaries, old style still compatible
    - Ranking system for database entries expanded to new 0-11 system from the old 0-5 system
    - Collision limit checking has been added for database entries

- Cantherm:
    - Improved support for MolPro output files
    - Added iodine support
    - Automatically read spin multiplicity from quantum output
    - Automatically assign frequency scale factor for supported model chemistries
    - Plot calculated rates and thermo by default
    - New sensitivity analysis feature analyzes sensitivity of reaction rates to isomer/TS energies in pdep networks

- Fixes:
    - Properly update charges when creating product templates in reaction families
    - Excessive duplicate reactions from different resonance structures has been fixed (bug introduced in 2.1.3)
    - Fixed rate calculation for MultiPdepArrhenius objects when member rates have different plists

- A more formal deprecation process is now being trialed. Deprecation warnings have been added to functions to be removed in version 2.3.0:
    - All methods related to saving or reading RMG-Java databases and old-style adjacency lists
    - The group additivity method for kinetics estimation (unrelated to thermo group additivity)
    - The saveRestartPeriod option and the old method of saving restart files

RMG-database Version 2.2.0
==========================
Date: July 5, 2018

- Additions:
    - New Intra_R_Add_Exo_Scission reaction family
    - New 1,2_ShiftC reaction family
    - New reaction families for peroxide chemistry in liquid systems
        - Korcek_step1_cat
        - Bimolec_Hydroperoxide_Decomposition
        - Peroxyl_Termination
        - Peroxyl_Disproportionation
        - Baeyer-Villiger_step1_cat
        - Baeyer-Villiger_step2
        - Baeyer-Villiger_step2_cat
    - Numerous new training reactions added to many families

- Changes:
    - New tree structure for Intra_R_Add_Endocyclic with consideration for cyclic species
    - Multiple bond on ring is no longer allowed in Intra_R_Add_Exocyclic and should react in Intra_R_Add_Endocyclic instead
    - Entry ranks rescaled to new 0-11 ranking system
    - Global forbidden structures has been cleaned up, leading to significant performance improvement

- Fixes:
    - Corrected shape indices in NOx2018 transport library
    - Removed or corrected some kinetics entries based on collision limit check


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