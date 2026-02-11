.. _installation:

************
Installation
************

.. NOTE::
    For questions related to installing or using RMG please post an issue on the `RMG-Py GitHub repository <https://github.com/ReactionMechanismGenerator/RMG-Py/issues>`_.
    You can also search for your problem on the issues page to see if there are already solutions in development.  Alternatively, you can email us at
    rmg_dev@mit.edu

Recommended Install: Docker
===========================

RMG is primarily distributed using Docker, a software package for delivering applications.

#. Download and install `Docker <https://docs.docker.com/get-started/get-docker/>`_.

#. Open a terminal, powershell, or command prompt and run ``docker pull reactionmechanismgenerator/rmg:3.3.0``.

   This step may take some time as the image is downloaded.

#. Run ``docker run --name rmgcontainer -v "C:\Users\rmguser\myrmgfiles:/rmg/RMG-Py/myrmgfiles" -it reactionmechanismgenerator/rmg:3.3.0``

   This command will make the folder ``C:\Users\rmguser\myrmgfiles`` on your computer accessible from inside the container to easily edit and transfer input and output files.
   Change the path to match your individual computer.
   If the folder does not exist when the command is run, it will be created.

   If you want to use jupyter notebook inside the docker container, run ``docker run --name rmgcontainer -v "C:\Users\rmguser\myrmgfiles:/rmg/RMG-Py/myrmgfiles" -it -p 8888:8888 reactionmechanismgenerator/rmg:3.3.0`` instead.
   And you can start the jupyter notebook by running ``jupyter notebook --ip 0.0.0.0 --no-browser --allow-root`` inside the container.
   Then you can access the jupyter notebook from your browser by going to ``http://localhost:8888``.
   You may need to copy and paste the token from the terminal into the browser to access the notebook.

You are now operating inside an Ubuntu operating system (a container called "rmgcontainer") with a working installation of RMG-Py.
To leave this container run ``exit``, and to reconnect run ``docker start rmgcontainer --attach --interactive``.

For users unfamiliar with bash or Linux, we recommend looking at
`online Linux tutorials <https://www.guru99.com/unix-linux-tutorial.html>`_ particularly `Linux vs. Windows <https://www.guru99.com/linux-differences.html>`_,
`Terminal vs File Manager <https://www.guru99.com/terminal-file-manager.html>`_, and
`Must Know Linux/Unix Commands <https://www.guru99.com/must-know-linux-commands.html>`_.

.. NOTE::
    The docker instructions above which specify ``rmg:3.3.0`` are for the version 3.3.0 release of RMG-Py (July 2025), the user guide for which is archived in a `4mb PDF here <https://github.com/ReactionMechanismGenerator/RMG-Py/blob/3.3.0/documentation/RMG-Py_and_Arkane_Documentation.pdf?raw=true>`_.
    If you want to use the latest development version, which corresponds to the user guide you are reading online `on the RMG website <https://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/>`_, you can replace the version number in the docker pull command with "latest" or omit it entirely.
    For example, you can run ``docker pull reactionmechanismgenerator/rmg:latest`` or ``docker pull reactionmechanismgenerator/rmg`` to get the latest development version corresponding to `the main branch <https://github.com/ReactionMechanismGenerator/RMG-Py/commits/main/>`_.
    Be sure to also change the version number in the docker run commands to match the version you pulled.
    A full list of available docker images can be found at `Docker Hub <https://hub.docker.com/r/reactionmechanismgenerator/rmg/tags>`_.

Alternative Install: Binary Installation Using Anaconda
===========================================================

If you are accustomed to using the Anaconda package manager or cannot tolerate the storage overhead of Docker, installation from Conda is also available.
This is recommended for users who want to use RMG out of the box and are not interested in changing the RMG code or making many additions to RMG's thermodynamic and kinetics databases. 
Installation by Conda without compiling from source is limited to `specific versions of RMG-Py <https://anaconda.org/RMG/rmg>`_ .
If you would like to use the latest version or plan to modify the code or databases, please refer to the Developer Installation instructions below.

.. toctree::
    :maxdepth: 1

    anacondaUser

Developer Install: Installation from Source
===========================================================

RMG-Py can be built from source using the Anaconda package manager to assist in installing all necessary dependencies. 
This is recommended for a developer who may be altering the RMG source code or someone who expects to manipulate the databases. 
You will also be able to access the latest source code updates and patches through GitHub.

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
    virtualMachineSetup
    linuxSubsystem

Dependencies
============

Please visit the page below for detailed information on all of RMG's dependencies and their license restrictions

.. toctree::
    :maxdepth: 1
    
    dependencies
