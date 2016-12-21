.. _releaseNotes:

*************
Release Notes
*************

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