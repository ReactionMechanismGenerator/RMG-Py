.. _anacondaDeveloper:

*******************************************************************************************
Installation by Source Using Anaconda Environment for Unix-based Systems: Linux and Mac OSX
*******************************************************************************************

#. Download and install the `Anaconda Python Platform <https://www.anaconda.com/download/>`_ for Python 2.7 (make sure not to install Python 3.0+, which is incompatible with RMG).

   The download will be a .sh file with a name like ``Anaconda2-2018.12-Linux-x86_64.sh``. Open a terminal in the same
   directory as this file, and type the following to install Anaconda (replace the name of your .sh file below). ::

    bash Anaconda2-2018.12-Linux-x86_64.sh

   **When prompted to append Anaconda to your PATH, select or type Yes**.  Install the Anaconda folder inside your home directory (typically ``/home/YourUsername/`` in Linux and ``/Users/YourUsername`` in Mac). When prompted, you do NOT need to install Microsoft VSCode (but feel free to if you are looking for a lightweight IDE).

#. Install `Git <https://git-scm.com/>`_, the open source version control package through the Terminal. **For Mac OS X**: Git is already packages with OS X 10.9 or later, but requires installation of Xcode's Command Line Tools. Skip the git installation and run it through the terminal, where you will be prompted to install the Command Line Tools if they are not already installed. ::

    sudo apt-get install git
    

#. Make sure that you also have gcc and g++, and make installed (run the lines below if you are uncertain). ::

    sudo apt install gcc
    sudo apt install g++
    sudo apt install make

#. Install the latest versions of RMG and RMG-database through cloning the source code via Git. Make sure to start in an appropriate local directory where you want both RMG-Py and RMG-database folders to exist. ::

    git clone https://github.com/ReactionMechanismGenerator/RMG-Py.git
    git clone https://github.com/ReactionMechanismGenerator/RMG-database.git

#. Now create the anaconda environment for RMG-Py

   For Linux users: ::
    
    cd RMG-Py
    source ~/.bashrc
    conda env create -f environment_linux.yml
    
   For Mac users: ::
         
    cd RMG-Py
    source ~/.bash_profile
    conda env create -f environment_mac.yml

#. Compile RMG-Py after activating the anaconda environment ::

    source activate rmg_env
    make
    
#. Modify environment variables. Add RMG-Py to the PYTHONPATH to ensure that you can access RMG modules from any folder. Also, add your RMG-Py folder to PATH to launch ``rmg.py`` from any folder. **Modify your** ``~/.bashrc`` **file by adding the following lines**: ::

    export PYTHONPATH=$PYTHONPATH:YourFolder/RMG-Py/
    export PATH=$PATH:YourFolder/RMG-Py/

   NOTE: Make sure to change ``YourFolder`` to the path leading to the ``RMG-Py`` code. Not doing so will lead to an error stating that python cannot find the module ``rmgpy``.

   be sure to either close and reopen your terminal to refresh your environment variables, or type the following command ::
 
    source ~/.bashrc

#. Finally, you can run RMG from any location by typing the following (given that you have prepared the input file as ``input.py`` in the current folder). ::

    rmg.py input.py
   
#. Optional: If you wish to use the :ref:`QMTP interface <qm>` with `MOPAC <http://openmopac.net/>`_ to run quantum mechanical calculations for improved thermochemistry estimates of cyclic species, please obtain a legal license through the `MOPAC License Request Form <http://openmopac.net/form.php>`_.  Once you have it, type the following into your Terminal ::
    
    mopac password_string_here    

You may now use RMG-Py, Arkane, as well as any of the :ref:`Standalone Modules <modules>` included in the RMG-Py package.




Test Suite
==========

There are a number of basic tests you can run on the newly installed RMG.  It is recommended to run them regularly to ensure the code and databases are behaving normally.  

#. **Unit test suite**: this will run all the unit tests in the ``rmgpy`` package ::

    cd RMG-Py
    make test
    
    
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
