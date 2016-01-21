.. _anacondaDeveloper:

*******************************************************************************************
Installation by Source Using Anaconda Environment for Unix-based Systems: Linux and Mac OSX
*******************************************************************************************

RMG-Py can now be built by source using the Anaconda Python Platform to assist in installing
all necessary dependencies. This is recommended for a developer who may be altering the RMG source code
or someone who expects to manipulate the databases extensively.  You will also be able to access the latest
source code updates and patches through Github.

* Download and install the `Anaconda Python Platform <http://continuum.io/downloads>`_ for Python 2.7 (make sure not to install Python 3.0+, which is incompatible with RMG). When prompted to append Anaconda to your PATH, select or type Yes.

* Install `Git <https://git-scm.com/>`_, the open source version control package through the Terminal. **For Mac OS X**: Git is already packages with OS X 10.9 or later, but requires installation of Xcode's Command Line Tools. Skip the git installation and run it through the terminal, where you will be prompted to install the Command Line Tools if they are not already installed. ::

    sudo apt-get install git
    
* Install the latest versions of RMG and RMG-database through cloning the source code via Git. Make sure to start in an appropriate local directory where you want both RMG-Py and RMG-database folders to exist. ::

    git clone https://github.com/ReactionMechanismGenerator/RMG-Py.git
    git clone https://github.com/ReactionMechanismGenerator/RMG-database.git
    
* Compile RMG-Py ::
    
    cd RMG-Py
    conda env create
    source activate rmg_env
    make
    
* Modify environment variables. Add RMG-Py to the PYTHONPATH to ensure that you can access RMG modules from any python prompt.  Modify your ``~/.bashrc`` file by adding the following line ::

   export PYTHONPATH=$PYTHONPATH:YourFolder/RMG-Py/
   
* If you wish to always be able to run RMG-Py, you can modify the anaconda path to point to the RMG environment. Modify the following line in your ``~/.bashrc`` file ::

   export PATH=~/anaconda/bin:$PATH
   
  by changing it to the following line :: 

   export PATH=~/anaconda/envs/rmg/bin:$PATH
   
* Optional: If you wish to use the :ref:`QMTP interface <qm>` with `MOPAC <http://openmopac.net/>`_ to run quantum mechanical calculations for improved thermochemistry estimates of cyclic species, please obtain a legal license through the `MOPAC License Request Form <http://openmopac.net/form.php>`_.  Once you have it, type the following into your Terminal ::
    
    mopac password_string_here    

You may now use RMG-Py, CanTherm, as well as any of the :ref:`Standalone Modules <modules>` included in the RMG-Py package.




Test Suite
==========

There are a number of basic tests you can run on the newly installed RMG.  It is recommended to run them regularly to ensure the code and databases are behaving normally.  

* **Unit test suite**: this will run all the unit tests in the ``rmgpy`` package ::

    cd RMG-Py
    make test
    
    
* **Database test suite**: this will run the database unit tests to ensure that groups, rate rules, and libraries are well formed ::

    cd RMG-Py
    make test-database
    

Running Examples
================

A number of basic examples can be run immediately.  Additional example input files can be found in the ``RMG-Py/examples`` folder.  Please read more on :ref:`Example Input Files <examples>` in the documentation.
    
* **Minimal Example**: this will run an Ethane pyrolysis model.  It should take less than a minute to complete. The results will be in the ``RMG-Py/testing/minimal`` folder::

    cd RMG-Py
    make eg1
    
* **Hexadiene Example**: this will run a Hexadiene model with pressure dependence and QMTP.  Note that you must have MOPAC installed for this to run. The results will be in the ``RMG-Py/testing/hexadiene`` folder::

    cd RMG-Py
    make eg2
    
* **Liquid Phase Example**: this will run a liquid phase RMG model.  The results will be in the ``RMG-Py/testing/liquid_phase`` folder ::

    cd RMG-Py
    make eg3
    
* **ThermoEstimator Example**: this will run the :ref:`Thermo Estimation Module <thermoModule>` on a few molecules. Note that you must have MOPAC installed for this to run completely. The results will be in the ``RMG-Py/testing/thermoEstimator`` folder ::

    cd RMG-Py
    make eg4
