.. _anacondaDeveloper:

*******************************************************************************************
Installation by Source Using Anaconda Environment for Unix-based Systems: Linux and Mac OSX
*******************************************************************************************

#. Download and install the `Anaconda Python Platform <https://www.anaconda.com/download/>`_ for Python 3.7.

   The download will be a .sh file with a name like ``Anaconda3-2019.07-Linux-x86_64.sh``. Open a terminal in the same
   directory as this file, and type the following to install Anaconda (replace the name of your .sh file below). ::

    bash Anaconda3-2019.07-Linux-x86_64.sh

   **When prompted to append Anaconda to your PATH, select or type Yes**.  Install the Anaconda folder inside your home
   directory (typically ``/home/YourUsername/`` in Linux and ``/Users/YourUsername`` in Mac). When prompted, you do not
   need to install Microsoft VSCode (but feel free to if you are looking for a lightweight IDE).

   Note that you should restart your terminal in order for the changes to take effect, as the installer will tell you.

#. There are a few system-level dependencies which are required and should not be installed via Anaconda. These include
   `Git <https://git-scm.com/>`_ for version control and GNU ``make``, ``gcc``, and ``g++`` for compiling RMG.

   For Linux users, you can check whether these are already installed by simply calling them via the command line, which
   will let you know if they are missing. To install any missing packages, you should use the appropriate package manager
   for your system. For example, on Ubuntu you would use the ``apt`` package manager. For Ubuntu 16 and newer, the
   necessary commands would be ::

    sudo apt install git
    sudo apt install gcc
    sudo apt install g++
    sudo apt install make

   For MacOS users, these packages will not come preinstalled, but can be easily obtained by installing the XCode Command Line Tools.
   These are a set of packages relevant for software development which have been bundled together by Apple. The easiest way
   to install this is to simply run one of the commands in the terminal, e.g. ``git``. The terminal will then prompt you on
   whether or not you would like to install the Command Line Tools.

#. Install the latest versions of RMG and RMG-database through cloning the source code via Git. Make sure to start in an
   appropriate local directory where you want both RMG-Py and RMG-database folders to exist. ::

    git clone https://github.com/ReactionMechanismGenerator/RMG-Py.git
    git clone https://github.com/ReactionMechanismGenerator/RMG-database.git

#. Now create the conda environment for RMG-Py ::

    cd RMG-Py
    conda env create -f environment.yml

   If the command errors due to being unable to find the `conda` command, try closing and re-opening your terminal
   window in order for the Anaconda settings to take effect.

#. Compile RMG-Py after activating the conda environment ::

    conda activate rmg_env
    make

   Note regarding differences between conda versions: Prior to Anaconda 4.4, the command to activate an environment was
   ``source activate rmg_env``. It has since been changed to ``conda activate rmg_env`` due to underlying changes to
   standardize operation across different operating systems. However, a prerequisite to using the new syntax is having
   run the ``conda init`` setup routine, which can be done at the end of the install procedure if the user requests.

#. Modify environment variables. Add RMG-Py to the PYTHONPATH to ensure that you can access RMG modules from any folder.
   Also, add your RMG-Py folder to PATH to launch ``rmg.py`` from any folder.

   In general, these commands should be placed in the appropriate shell initialization file. For Linux users using
   bash (the default on Ubuntu), these should be placed in ``~/.bashrc``. For MacOS users using bash (default before MacOS Catalina),
   these should be placed in ``~/.bash_profile``, which you should create if it doesn't exist. For MacOS users using zsh
   (default beginning in MacOS Catalina), these should be placed in ``~/.zshrc``. ::

    export PYTHONPATH=YourFolder/RMG-Py/:$PYTHONPATH
    export PATH=YourFolder/RMG-Py/:$PATH

   NOTE: Make sure to change ``YourFolder`` to the path leading to the ``RMG-Py`` code. Not doing so will lead to an error stating that python cannot find the module ``rmgpy``.

   Be sure to either close and reopen your terminal to refresh your environment variables, or type the following command ::

    source ~/.bashrc

#. Finally, you can run RMG from any location by typing the following (given that you have prepared the input file as ``input.py`` in the current folder). ::

    rmg.py input.py
   
#. Optional: If you wish to use the :ref:`QMTP interface <qm>` with `MOPAC <http://openmopac.net/>`_ to run quantum mechanical calculations for improved thermochemistry estimates of cyclic species, please obtain a legal license through the `MOPAC License Request Form <http://openmopac.net/form.php>`_.  Once you have it, type the following into your terminal ::
    
    mopac password_string_here    

You may now use RMG-Py, Arkane, as well as any of the :ref:`Standalone Modules <modules>` included in the RMG-Py package.




Test Suite
==========

There are a number of basic tests you can run on the newly installed RMG.  It is recommended to run them regularly to ensure the code and databases are behaving normally.  

#. **Unit test suite**: this will run all the unit tests in the ``rmgpy`` and ``arkane`` packages ::

    cd RMG-Py
    make test
    
#. **Functional test suite**: this will run all the functional tests in the ``rmgpy`` and ``arkane`` packages ::

    cd RMG-Py
    make test-functional


#. **Database test suite**: this will run the database unit tests to ensure that groups, rate rules, and libraries are well formed ::

    cd RMG-Py
    make test-database
    

Running Examples
================

A number of basic examples can be run immediately.  Additional example input files can be found in the ``RMG-Py/examples`` folder.  Please read more on :ref:`Example Input Files <examples>` in the documentation.
    
#. **Minimal Example**: this will run an Ethane pyrolysis model.  It should take less than a minute to complete. The results will be in the ``RMG-Py/testing/minimal`` folder::

    cd RMG-Py
    make eg1
    
#. **Hexadiene Example**: this will run a Hexadiene model with pressure dependence and QMTP.  Note that you must have MOPAC installed for this to run. The results will be in the ``RMG-Py/testing/hexadiene`` folder::

    cd RMG-Py
    make eg2
    
#. **Liquid Phase Example**: this will run a liquid phase RMG model.  The results will be in the ``RMG-Py/testing/liquid_phase`` folder ::

    cd RMG-Py
    make eg3
    
#. **ThermoEstimator Example**: this will run the :ref:`Thermo Estimation Module <thermoModule>` on a few molecules. Note that you must have MOPAC installed for this to run completely. The results will be in the ``RMG-Py/testing/thermoEstimator`` folder ::

    cd RMG-Py
    make eg4
