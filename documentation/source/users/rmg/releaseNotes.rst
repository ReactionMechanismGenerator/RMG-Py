.. _releaseNotes:

*************
Release Notes
*************


Version 1.0.2
=============
Date: February 1, 2016

This new release brings up several new features and bug fixes.

- Windows users can rejoice: RMG is now available in binary format on the Anaconda platform.  Building by source is also
much easier now through the Anaconda managed python environment for dependencies. See the updated :ref:`Installation Page<installation>`
for more details
- Reaction filtering for speeding up model generation has now been added.  It has been shown to speed up model convergence by
7-10x.  See more details about how to use it in your RMG job :ref:`here <filterReactions>`.  Learn more about the theory 
and algorithm on the :ref:`Rate-based Model Enlarging Algorithm <ratebasedmodelenlarger>` page.
- The RMG :ref:`native scripts <modules>` are now organized under the `rmgpy.tools` submodule for
  developer ease and better extensibility in external scripts.
- InChI conversion is now more robust for singlets and triplets, 
  and augmented InChIs and InChI keys are now possible with new radical electron, lone pair, and multiplicity flags.  
- Output HTML for visualizing models are now cleaned up and also more functional, including features to display thermo comments,
  display enthalpy, entropy, and free energy of reaction, as well as filter reactions by species.  You can use this new html format
  either by running a job in RMG v1.0.2 or revisualizing your CHEMKIN file and species dictionary using
  the :ref:`visualization web tool <http://rmg.mit.edu/simulate/chemkin>`_.

