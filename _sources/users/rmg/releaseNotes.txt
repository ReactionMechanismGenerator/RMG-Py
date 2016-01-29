.. _releaseNotes:

*************
Release Notes
*************


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