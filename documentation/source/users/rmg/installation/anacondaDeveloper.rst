.. _anacondaDeveloper:

*******************************************************************************************
Installation by Source Using Anaconda Environment for Unix-based Systems: Linux and Mac OSX
*******************************************************************************************

#. Install the `conda` package manager, if you do not already have it (or Anaconda).
   Select one of the following options:

   a. Users of Fedora Linux and Red Hat derivatives (RHEL, CentOS Stream) may install from the official repositories and EPEL, respectively, with the command ::

       sudo dnf install conda

   b. All other users, download and install `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_.

      The download will be a .sh file with a name like ``Miniconda3-latest-Linux-x86_64.sh``.
      Open a terminal in the same directory as this file, and type the following to install Conda
      (replace the name of your .sh file below). ::

       bash Miniconda3-latest-Linux-x86_64.sh

      **When prompted to append Anaconda to your PATH, select or type Yes**. 
      Install the Conda folder inside your home directory 
      (typically ``/home/YourUsername/`` in Linux and ``/Users/YourUsername`` in Mac).

      Note that you should reinitialize or restart your terminal in order for the changes to take effect, as the installer will tell you.

#. There are a few system-level dependencies which are required and should not be installed via Conda. These include
   `Git <https://git-scm.com/>`_ for version control, `GNU Make <https://www.gnu.org/software/make/>`_, and the C and C++ compilers from the `GNU Compiler Collection (GCC) <https://gcc.gnu.org/>`_ for compiling RMG.

   For Linux users, you can check whether these are already installed by simply calling them via the command line, which
   will let you know if they are missing. To install any missing packages, you should use the appropriate package manager
   for your system.

   a. On Ubuntu and Debian the package manager is ``apt`` ::

       sudo apt install git gcc g++ make

   b. On Fedora and Red Hat derivatives (RHEL 8+) the package manager is ``dnf`` ::

       sudo dnf install git gcc gcc-c++ make

   c. For Red Hat 7 and lower, replace ``dnf`` with ``yum`` in the preceding.

   d. On openSUSE the package manager is ``zypper``::

       sudo zypper install git gcc gcc-c++ make

   e. On Manjaro or Arch Linux the package manager is ``pacman`` ::

       sudo pacman -S git gcc make

   f. For MacOS users, the above packages can be easily obtained by installing the XCode Command Line Tools.
      These are a set of packages relevant for software development which have been bundled together by Apple.
      The easiest way to install this is to simply run one of the commands in the terminal, e.g. ``git``.
      The terminal will then prompt you to install the Command Line Tools.

#. Install the latest versions of RMG and RMG-database through cloning the source code via Git. Make sure to start in an
   appropriate local directory where you want both RMG-Py and RMG-database folders to exist.
   Github has deprecated password authentication from the command line, so it
   is preferred to clone the repositories using ``ssh``::

    git clone git@github.com:ReactionMechanismGenerator/RMG-Py.git
    git clone git@github.com:ReactionMechanismGenerator/RMG-database.git

   It is still possible to clone the repositories using ``https`` if you are
   unfamiliar with ``ssh``::
   
    git clone https://github.com/ReactionMechanismGenerator/RMG-Py.git
    git clone https://github.com/ReactionMechanismGenerator/RMG-database.git

   For information on using ``ssh`` with GitHub see the `Connecting to GitHub with SSH <https://docs.github.com/en/authentication/connecting-to-github-with-ssh>`_

#. Switch the conda solver backend to speed up creation of the RMG environment ::

    conda install -n base conda-libmamba-solver
    conda config --set solver libmamba

#. Navigate to the RMG-Py directory ::

    cd RMG-Py

#. Apple silicon (M1+) users only: execute the following commands
   **instead of** the following `conda env create -f environment.yml` step.
   (This will tell conda that we want to the environment to use x86 
   architecture rather than the native ARM64 architecture) ::

    conda create -n rmg_env
    conda activate rmg_env
    conda config --env --set subdir osx-64
    conda env update -f environment.yml

#. Create the conda environment for RMG-Py::

    conda env create -f environment.yml

   If either of these commands return an error due to being unable to find the ``conda`` command,
   try to either close and reopen your terminal to refresh your environment variables
   or type the following command.

   If on Linux or pre-Catalina MacOS (or if you have a bash shell)::

    source ~/.bashrc

   If on MacOS Catalina or later (or if you have a Z shell)::

    source ~/.zshrc

#. Activate conda environment ::

    conda activate rmg_env

#. Switch the conda solver to libmamba again, to accelerate any changes you might make to this conda environment in the future::

    conda config --set solver libmamba

#. Compile RMG-Py after activating the conda environment ::

    make

#. Modify environment variables. Add RMG-Py to the PYTHONPATH to ensure that you can access RMG modules from any folder.
   *This is important before the next step in which julia dependencies are installed.*
   Also, add your RMG-Py folder to PATH to launch ``rmg.py`` from any folder.

   In general, these commands should be placed in the appropriate shell initialization file.
   For Linux users using bash (the default on distributions mentioned here), these should be placed in ``~/.bashrc``.
   For MacOS users using bash (default before MacOS Catalina), these should be placed in ``~/.bash_profile``, which you should create if it doesn't exist.
   For MacOS users using zsh (default beginning in MacOS Catalina), these should be placed in ``~/.zshrc``. ::

    export PYTHONPATH=YourFolder/RMG-Py/:$PYTHONPATH
    export PATH=YourFolder/RMG-Py/:$PATH

   NOTE: Make sure to change ``YourFolder`` to the path leading to the ``RMG-Py`` code. Not doing so will lead to an error stating that python cannot find the module ``rmgpy``.

   Be sure to either close and reopen your terminal to refresh your environment variables (``source ~/.bashrc`` or ``source ~/.zshrc``).

#. Install and Link Julia dependencies: ::

     julia -e 'using Pkg; Pkg.add("PyCall");Pkg.build("PyCall");Pkg.add(PackageSpec(name="ReactionMechanismSimulator",rev="main")); using ReactionMechanismSimulator;'

     python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()"


#. Finally, you can run RMG from any location by typing the following (given that you have prepared the input file as ``input.py`` in the current folder). ::

    python replace/with/path/to/rmg.py input.py

You may now use RMG-Py, Arkane, as well as any of the :ref:`Standalone Modules <modules>` included in the RMG-Py package.


Debugging
=========

If you wish to debug using the (very helpful) debugger in `VSCode <https://code.visualstudio.com>`_,
here is an example launch configuration to put in your ``launch.json`` file,
which can be found in the ``.vscode`` folder.
You might have to edit them slightly to match your exact paths. Specifically, 
you will need ``/opt/miniconda3/envs/rmg_env`` to point to where your conda environment is located.

This configuration will allow you to debug the rms_constant_V example, running through
python. ::

        {
            "name": "Python: rmg.py rms_constant_V",
            "type": "python",
            "request": "launch",
            "cwd": "${workspaceFolder}/",
            "program": "rmg.py",
            "python": "/opt/miniconda3/envs/rmg_env/bin/python",
            "args": [
                "examples/rmg/rms_constant_V/input.py",
            ],
            "console": "integratedTerminal",
            "env": {
                "PATH": "/opt/miniconda3/envs/rmg_env/bin:${env:PATH}",
                "PYTHONPATH": "${workspaceFolder}/",
            }
        },

This configuration will allow you to debug a subset of the unit tests.
Open one of the many test files named ``*Test.py`` in ``test/rmgpy`` before you launch it::

        {
            "name": "Python: pytest Current File",
            "type": "python",
            "request": "launch",
            "program": "/opt/miniconda3/envs/rmg_env/bin/pytest",
            "python": "/opt/miniconda3/envs/rmg_env/bin/python",
            "args": [
                "--capture=no",
                "--verbose",
                "${file}"
            ],
            "console": "integratedTerminal",
            "env": {
                "PATH": "/opt/miniconda3/envs/rmg_env/bin:${env:PATH}",
                "PYTHONPATH": "${workspaceFolder}/",
            },
        },

This configuration will allow you to debug running all the database tests.::

        {
            "name": "Test RMG-database",
            "type": "python",
            "request": "launch",
            "program": "/opt/miniconda3/envs/rmg_env/bin/pytest",
            "python": "/opt/miniconda3/envs/rmg_env/bin/python",
            "args": [
                "--capture=no",
                "--verbose",
                "${workspaceFolder}/test/database/databaseTest.py"
            ],
            "console": "integratedTerminal",
            "env": {
                "PATH": "/opt/miniconda3/envs/rmg_env/bin:${env:PATH}",
                "PYTHONPATH": "${workspaceFolder}/",
            },
        },

This configuration will allow you to use the debugger breakpoints inside unit tests being run by the pytest framework::

        {
            "name": "Python: Debug Tests",
            "type": "python",
            "request": "launch",
            "program": "${file}",
            "purpose": ["debug-test"],
            "python": "/opt/miniconda3/envs/rmg_env/bin/python",
            "console": "integratedTerminal",
            "justMyCode": false,
            "env": {"PYTEST_ADDOPTS": "--no-cov",} // without disabling coverage VS Code doesn't stop at breakpoints while debugging because pytest-cov is using the same technique to access the source code being run
          }

See more about testing in VSCode in the :ref:`Testing in VSCode <vscode_testing>` section below.

Test Suite
==========

There are a number of basic tests you can run on the newly installed RMG.  It is recommended to run them regularly to ensure the code and databases are behaving normally.
Make sure that the environment is active before running the tests: ``conda activate rmg_env``.

#. **Unit test suite**: this will run all the unit tests in the ``rmgpy`` and ``arkane`` packages ::

    cd RMG-Py
    make test
    
#. **Functional test suite**: this will run all the functional tests in the ``rmgpy`` and ``arkane`` packages ::

    cd RMG-Py
    make test-functional


#. **Database test suite**: this will run the database unit tests to ensure that groups, rate rules, and libraries are well-formed ::

    cd RMG-Py
    make test-database
    

.. _vscode_testing:

Testing in VSCode
=================

Once you have the Python extension installed and a Python file open within the editor, 
a test beaker icon will be displayed on the VS Code Activity bar. 
The beaker icon is for the Test Explorer view. When opening the Test Explorer, 
you will see a Configure Tests button if you don't have a test framework enabled.
Once you select Configure Tests, you will be prompted to select a test framework 
(**select `pytest`**)
and a folder containing the tests
(**select `test`**).
To configure the rest of the settings, find the ``settings.json`` file in your ``.vscode`` folder.
You can use the following settings to configure the pytest framework::

    "python.testing.pytestEnabled": true,
    "python.testing.pytestPath": "python -m pytest",
    "python.testing.pytestArgs": [
        "-p", "julia.pytestplugin",
        "--julia-compiled-modules=no",
        "--ignore", "test/regression",
        "-m", "not functional",
        // "-n", "auto", // number of parallel processes, if you install pytest-xdist
        "test"
    ],

To run the tests, you can click the Run All Tests button in the Test Explorer view.
Learn more at the `Python testing in Visual Studio Code <https://code.visualstudio.com/docs/python/testing>`_ documentation.

Given the time taken for Julia to compile things every time it launches,
you might find this to be painfully slow even for a simple test.
It may be possible to use ``--julia-sysimage=JULIA_SYSIMAGE`` instead of ``--julia-compiled-modules=no``,
or disable PyJulia entirely.
If you find a better way to do this, or clearer instructions, 
please `update this section <https://github.com/ReactionMechanismGenerator/RMG-Py/edit/main/documentation/source/users/rmg/installation/anacondaDeveloper.rst>`_.


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


Building Documentation
======================
To build the documentation (to test that you have it right before pushing to GitHub) you will need to install sphinx::

    conda activate rmg_env
    conda install sphinx

Then you can build the documentation::

    make documentation
