.. _installation:

************
Installation
************

.. NOTE::
    It is recommended that RMG be installed with Python 3.7. Legacy releases compatible with Python 2.7 are also available.

For any questions related to RMG and its usage and installation, please
post an issue at https://github.com/ReactionMechanismGenerator/RMG-Py/issues and the RMG
developers will get back to you as soon as we can.  You can also search for your problem on the issues
page to see if there are already solutions in development.  Alternatively, you can email us at
rmg_dev@mit.edu

Installation on a Windows Platform
====================================

Due to difficulties with dependencies, installation on Windows directly is no longer supported. Instead, it is
recommended to run a Linux virtual machine from Windows and follow either the instructions for basic users
(binary installation using Anaconda) or the instructions for developers. Alternatively, it is also possible to install
RMG in the Ubuntu subsystem now available on Windows 10.

.. toctree::
    :maxdepth: 1

    virtualMachineSetup
    linuxSubsystem

For users unfamiliar with bash or Linux, we recommend looking at
`online Linux tutorials <https://www.guru99.com/unix-linux-tutorial.html>`_. To start out
with, we recommend looking at the following tutorials: `Linux vs. Windows <https://www.guru99.com/linux-differences.html>`_,
`Terminal vs File Manager <https://www.guru99.com/terminal-file-manager.html>`_, and
`Must Know Linux/Unix Commands <https://www.guru99.com/must-know-linux-commands.html>`_.


For Basic Users: Binary Installation Using Anaconda
===================================================

It is highly recommended to use the Python platform Anaconda to perform the installation of RMG-Py.
A binary installation is recommended for users who want to use RMG out of the box, and are
not interested in changing or recompiling the RMG code or making many additions to 
RMG's thermodynamic and kinetics databases. 
    
.. toctree::
    :maxdepth: 1
    
    anacondaUser


For Developers: Installation by Source Using Anaconda Environment
=================================================================

RMG-Py can now be built by source using the Anaconda Python Platform to assist in installing
all necessary dependencies. This is recommended for a developer who may be altering the RMG source code
or someone who expects to manipulate the databases extensively.  You will also be able to access the latest
source code updates and patches through Github.

.. toctree::
    :maxdepth: 1
    
    anacondaDeveloper
    updatingSourceCode

Archive of Unsupported Installation Methods
===========================================

Below are old installation techniques that are no longer supported, including instructions for installation without
using Anaconda and the old installation instructions for Windows. These instructions are no longer maintained, and are
not recommended for use.

.. toctree::
    :maxdepth: 1
    
    linux
    macos
    anacondaUserWindows
    anacondaDeveloperWindows
    windowsEnvironment

Dependencies
============

Please visit the page below for detailed information on all of RMG's dependencies and their license restrictions

.. toctree::
    :maxdepth: 1
    
    dependencies
