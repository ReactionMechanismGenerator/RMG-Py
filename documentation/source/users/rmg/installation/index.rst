.. _installation:

************
Installation
************



.. NOTE::
    It is recommended that RMG be installed with Python 2.7, although it has been previously tested that Python 2.5 and 2.6 may also work. Dependency issues render it incompatible with Python 3.x releases.


For any questions related to RMG and its usage and installation, please
post an issue at https://github.com/ReactionMechanismGenerator/RMG-Py/issues and the RMG
developers will get back to you as soon as we can.  You can also search for your problem on the issues
page to see if there are already solutions in development.  Alternatively, you can email us at
rmg_dev@mit.edu

For Basic Users: Binary Installation Using Anaconda
===================================================

It is highly recommended to use the Python platform Anaconda to perform the installation of RMG-Py.
A binary installation is recommended for users who want to use RMG out of the box, and are
not interested in changing or recompiling the RMG code or making many additions to 
RMG's thermodynamic and kinetics databases. 
    
.. toctree::
    :maxdepth: 1
    
    anacondaUser
    anacondaUserWindows
    windowsEnvironment


For Developers: Installation by Source Using Anaconda Environment
=================================================================

RMG-Py can now be built by source using the Anaconda Python Platform to assist in installing
all necessary dependencies. This is recommended for a developer who may be altering the RMG source code
or someone who expects to manipulate the databases extensively.  You will also be able to access the latest
source code updates and patches through Github.

.. toctree::
    :maxdepth: 1
    
    anacondaDeveloper
    anacondaDeveloperWindows
    windowsEnvironment
    updatingSourceCode

For Developers: Direct Installation by Source without Anaconda
=================================================================

The installation approach in this section is not recommended and also not maintained by RMG developer team. This is only a record for people who don't want to use Anaconda.

.. toctree::
    :maxdepth: 1
    
    linux
    macos

Dependencies
============

Please visit the page below for detailed information on all of RMG's dependencies and their license restrictions

.. toctree::
    :maxdepth: 1
    
    dependencies

Installation FAQ
================

This section collects frequently asked questions on installation of RMG.

.. toctree::
    :maxdepth: 1
    
    faq