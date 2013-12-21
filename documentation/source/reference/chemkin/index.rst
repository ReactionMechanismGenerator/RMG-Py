************************************
Chemkin files (:mod:`rmgpy.chemkin`)
************************************

.. module:: rmgpy.chemkin

The :mod:`rmgpy.chemkin` module contains functions for reading and writing of
Chemkin and Chemkin-like files.



Reading Chemkin files
=====================

.. currentmodule:: rmgpy.chemkin

=============================== ================================================
Function                        Description
=============================== ================================================
:func:`loadChemkinFile`         Load a reaction mechanism from a Chemkin file
:func:`loadSpeciesDictionary`   Load a species dictionary from a file
:func:`loadTransportFile`       Load a Chemkin transport properties file
------------------------------- ------------------------------------------------
:func:`readKineticsEntry`       Read a single reaction entry from a Chemkin file
:func:`readReactionComments`    Read the comments associated with a reaction entry  
:func:`readReactionsBlock`      Read the reactions block of a Chemkin file
:func:`readThermoEntry`         Read a single thermodynamics entry from a Chemkin file
:func:`removeCommentFromLine`   Remove comment text from a line of a Chemkin file or species dictionary
=============================== ================================================



Writing Chemkin files
=====================

.. currentmodule:: rmgpy.chemkin

=============================== ================================================
Function                        Description
=============================== ================================================
:func:`saveChemkinFile`         Save a reaction mechanism to a Chemkin file
:func:`saveSpeciesDictionary`   Save a species dictionary to a file
:func:`saveTransportFile`       Save a Chemkin transport properties file
:func:`saveHTMLFile`            Save an HTML file representing a Chemkin mechanism
:func:`saveJavaKineticsLibrary` Save a mechanism to a (Chemkin-like) kinetics library for RMG-Java
------------------------------- ------------------------------------------------
:func:`getSpeciesIdentifier`    Return the Chemkin-valid identifier for a given species
:func:`markDuplicateReactions`  Find and mark all duplicate reactions in a mechanism
:func:`writeKineticsEntry`      Write a single reaction entry to a Chemkin file
:func:`writeThermoEntry`        Write a single thermodynamics entry to a Chemkin file
=============================== ================================================



Exceptions
==========

.. currentmodule:: rmgpy.chemkin

=================================== ============================================
Exception                           Description
=================================== ============================================
:exc:`ChemkinError`                 Raised when an error occurs while working with a Chemkin file
=================================== ============================================



.. toctree::
    :hidden:
    
    reader
    writer
    
