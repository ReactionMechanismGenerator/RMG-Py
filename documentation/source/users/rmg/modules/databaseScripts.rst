.. _databaseScripts:

****************
Database Scripts
****************

This section details usage for scripts available in ``RMG-database/scripts`` folder.


evansPolanyi.py
---------------
This script will generate an Evans-Polanyi plot for a single kinetics
depository.

Usage::

    python evansPolanyi.py [-h] <family> <kinetics_depository> [<kinetics_depository> ...]

Positional arguments::
  <family>              the family to use
  <kinetics_depository> the kineticsDepository to use, e.g., training, NIST

Optional arguments::

    -h, --help    show help message and exit


exportKineticsLibraryToChemkin.py
---------------------------------
This script exports an individual RMG-Py kinetics library to a chemkin
and dictionary file.  Thermo is taken from RMG's estimates and libraries.  
In order to use more specific thermo, you must tweak the thermoLibraries and 
estimators in use when loading the database. The script will save the
chem.inp and species_dictionary.txt files in the local directory.

Usage::

    python exportKineticsLibrarytoChemkin.py [-h] LIBRARYNAME

Positional arguments::

    LIBRARYNAME    the libraryname of the RMG-Py format kinetics library

Optional arguments::

    -h, --help     show help message and exit


importChemkinLibrary.py
-----------------------
This script imports a chemkin file (along with RMG dictionary) from a local directory and saves a set of
RMG-Py kinetics library and thermo library files.  These py files are automatically added to the 
input/kinetics/libraries and input/thermo/libraries folder under the user-specified `name` for the chemkin library.

Usage::

    python importChemkinLibrary.py [-h] CHEMKIN DICTIONARY NAME

Positional arguments::

    CHEMKIN     The path of the chemkin file
    DICTIONARY  The path of the RMG dictionary file
    NAME        Name of the chemkin library to be saved

Optional arguments::

    -h, --help  show help message and exit


process_family_images.py
------------------------
This script processes reaction family template images (saved as eps files) into user friendly files (pdf and pngs).
This should typically be run whenever a new family is added or an existing family is updated.

Notes:
  - Make sure you have a working LaTeX installation with pdflatex
  - Make sure you have a working GhostScript installation for epstopdf
  - Make sure you have ImageMagick installed for png generation
  - ImageMagick may have security limitations in place which prevent reading
    eps files. To circumvent these, edit the ``/etc/ImageMagick-6/policy.xml``
    file by changing ``<policy domain="coder" rights="none" pattern="EPS" />``
    to ``<policy domain="coder" rights="read|write" pattern="EPS" />``

Usage::

    python process_family_images.py


