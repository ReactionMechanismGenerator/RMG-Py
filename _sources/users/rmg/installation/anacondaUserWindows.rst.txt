.. _anacondaUserWindows:

**********************************************
Binary Installation Using Anaconda for Windows
**********************************************

* Download and install the `Anaconda Python Platform <http://continuum.io/downloads>`_ for Python 2.7 (make sure not to install Python 3.0+, which is incompatible with RMG). We recommend changing the default install path to ``C:\Anaconda\`` in order to avoid spaces in the install path and be easily accessible. It is recommended to append Anaconda to your PATH as well as setting it as your default Python executable.  All other settings can remain as their defaults.

.. image:: images/AnacondaInstallWindows.png
    :align: center

* Now we want to install both RMG and the RMG-database binaries via the command prompt.  Dependencies will be installed automatically.  It is safest to make a new Anaconda environment for RMG and all its dependencies. Open a command prompt (either by finding it in your Program Files or by searching for ``cmd.exe``) and type the following to create the new environment named 'rmg_env' containing the latest stable version of the RMG program and its database. ::

    conda create -c rmg --name rmg_env rmg rmgdatabase
    
* Whenever you wish to use it you must first activate the environment in the command prompt by typing::
    
    activate rmg_env
    
* Optional: If you wish to use the :ref:`QMTP interface <qm>` with `MOPAC <http://openmopac.net/>`_ to run quantum mechanical calculations for improved thermochemistry estimates of cyclic species, please obtain a legal license through the `MOPAC License Request Form <http://openmopac.net/form.php>`_.Once you have it, type the following into your command prompt (while the environment is activated) ::
    
    mopac password_string_here

* Now you must :ref:`set the RMG environment variable in Windows <windowsEnvironment>` to allow your system to find the RMG python files more easily.  

* If you set any new environment variables, you must now close and reopen the command prompt so that those environment variables can be refreshed and used.

* You may now run an RMG test job. Save the `Minimal Example Input File <https://raw.githubusercontent.com/ReactionMechanismGenerator/RMG-Py/master/examples/rmg/minimal/input.py>`_ to a local directory.  Use the command prompt to run your RMG job inside that folder by using the following command ::

    activate rmg_env
    python %RMGPy%\rmg.py input.py

You may now use RMG-Py, CanTherm, as well as any of the :ref:`Standalone Modules <modules>` included in the RMG-Py package.


Updating your binary installation of RMG for Windows
====================================================

If you had previously installed a binary version of the RMG package, you may
check and update your installation to the latest stable version available on Anaconda by typing the following command in a Command Prompt ::

    source activate rmg_env
    conda update rmg rmgdatabase -c rmg 
