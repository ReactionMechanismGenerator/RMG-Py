.. _anacondaDeveloperWindows:

*************************************************************
Installation by Source Using Anaconda Environment for Windows
*************************************************************


* Download and install the `Anaconda Python Platform <http://continuum.io/downloads>`_ for Python 2.7 (make sure not to install Python 3.0+, which is incompatible with RMG). We recommend changing the default install path to ``C:\Anaconda\`` in order to avoid spaces in the install path and be easily accessible. It is recommended to append Anaconda to your PATH as well as setting it as your default Python executable.  All other settings can remain as their defaults.

.. image:: images/AnacondaInstallWindows.png
    :align: center

* Install `Git <http://git-scm.com/download/win>`_, the open source version control package. When asked, append Git tools to your Command Prompt. It is also recommended to commit Unix-style line endings:

.. image:: images/InstallGit.png
    :align: center
    
* Open Git CMD or a command prompt (either by finding it in your Program Files or by searching for ``cmd.exe``).  Install the latest versions of RMG and RMG-database through cloning the source code via Git. Make sure to start in an appropriate local directory where you want both RMG-Py and RMG-database folders to exist. We recommend creating a folder such as ``C:\Code\`` ::

    git clone https://github.com/ReactionMechanismGenerator/RMG-Py.git
    git clone https://github.com/ReactionMechanismGenerator/RMG-database.git
    
* Create and activate the RMG Anaconda environment ::
    
    cd RMG-Py
    conda env create -f environment_windows.yml
    activate rmg_env
    
  Every time you open a new command prompt and want to use RMG, you must reactivate this environment by typing ``activate rmg_env``.

* Now you can compile RMG-Py ::
    
    cd RMG-Py
    mingw32-make
    
* Now it is recommended to modify your system's environment variables.  Please see :ref:`Setting the RMG environment variable in Windows <windowsEnvironment>` for more information.  

  Additionally, set the ``PYTHONPATH`` environment variable to the path of your RMG-Py source folder to ensure that you can access RMG modules from any python prompt.  The prompt might look like this: 

    .. image:: images/AnacondaInstallWindows.png
        :align: center

* If you set any new environment variables, you must now close and reopen the command prompt so that those environment variables can be refreshed and used.
   
* Optional: If you wish to use the :ref:`QMTP interface <qm>` with `MOPAC <http://openmopac.net/>`_ to run quantum mechanical calculations for improved thermochemistry estimates of cyclic species, please obtain a legal license through the `MOPAC License Request Form <http://openmopac.net/form.php>`_.  Once you have it, type the following into your command prompt ::
    
     mopac password_string_here    

You may now use RMG-Py, CanTherm, as well as any of the :ref:`Standalone Modules <modules>` included in the RMG-Py package.



Test Suite
==========

There are a number of basic tests you can run on the newly installed RMG.  It is recommended to run them regularly to ensure the code and databases are behaving normally.  

* **Unit test suite**: this will run all the unit tests in the ``rmgpy`` package ::

    cd RMG-Py
    mingw32-make test
    
    
* **Database test suite**: this will run the database unit tests to ensure that groups, rate rules, and libraries are well formed ::

    cd RMG-Py
    mingw32-make test-database
    

Running Examples
================

A number of basic examples can be run immediately.  Additional example input files can be found in the ``RMG-Py\examples`` folder.  Please read more on :ref:`Example Input Files <examples>` in the documentation.
    
* **Minimal Example**: this will run an Ethane pyrolysis model.  It should take less than a minute to complete. The results will be in the ``RMG-Py\testing\minimal`` folder::

    cd RMG-Py
    mingw32-make eg1
    
* **Hexadiene Example**: this will run a Hexadiene model with pressure dependence and QMTP.  Note that you must have MOPAC installed for this to run. The results will be in the ``RMG-Py\testing\hexadiene`` folder::

    cd RMG-Py
    mingw32-make eg2
    
* **Liquid Phase Example**: this will run a liquid phase RMG model.  The results will be in the ``RMG-Py\testing\liquid_phase`` folder ::

    cd RMG-Py
    mingw32-make eg3
    
* **ThermoEstimator Example**: this will run the :ref:`Thermo Estimation Module <thermoModule>` on a few molecules. Note that you must have MOPAC installed for this to run completely. The results will be in the ``RMG-Py\testing\thermoEstimator`` folder ::

    cd RMG-Py
    mingw32-make eg4
