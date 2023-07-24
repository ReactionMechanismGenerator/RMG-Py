************************************
Chemkin files (:mod:`rmgpy.chemkin`)
************************************

.. module:: rmgpy.chemkin

The :mod:`rmgpy.chemkin` module contains functions for reading and writing of
Chemkin and Chemkin-like files.



Reading Chemkin files
=====================

.. currentmodule:: rmgpy.chemkin

================================ ================================================
Function                         Description
================================ ================================================
:func:`load_chemkin_file`        Load a reaction mechanism from a Chemkin file
:func:`load_species_dictionary`  Load a species dictionary from a file
:func:`load_transport_file`      Load a Chemkin transport properties file
-------------------------------- ------------------------------------------------
:func:`read_kinetics_entry`      Read a single reaction entry from a Chemkin file
:func:`read_reaction_comments`   Read the comments associated with a reaction entry
:func:`read_reactions_block`     Read the reactions block of a Chemkin file
:func:`read_thermo_entry`        Read a single thermodynamics entry from a Chemkin file
:func:`remove_comment_from_line` Remove comment text from a line of a Chemkin file or species dictionary
================================ ================================================



Writing Chemkin files
=====================

.. currentmodule:: rmgpy.chemkin

================================== ================================================
Function                           Description
================================== ================================================
:func:`save_chemkin_file`          Save a reaction mechanism to a Chemkin file
:func:`save_species_dictionary`    Save a species dictionary to a file
:func:`save_transport_file`        Save a Chemkin transport properties file
---------------------------------- ------------------------------------------------
:func:`get_species_identifier`     Return the Chemkin-valid identifier for a given species
:func:`mark_duplicate_reactions`   Find and mark all duplicate reactions in a mechanism
:func:`write_kinetics_entry`       Write a single reaction entry to a Chemkin file
:func:`write_thermo_entry`         Write a single thermodynamics entry to a Chemkin file
================================== ================================================



.. toctree::
    :hidden:
    
    reader
    writer
    
