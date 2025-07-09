.. _anacondaUser:

****************************************************************************
Binary Installation Using Anaconda for Unix-Based Systems: Linux and Mac OSX
****************************************************************************


#. Install the `conda` package manager via `miniforge`, if you do not already have it (or Anaconda), by following the `Miniforge installation instructions <https://github.com/conda-forge/miniforge?tab=readme-ov-file#install>`_.

#. Install both RMG and the RMG-database binaries through the terminal. Dependencies will be installed automatically. It is safest to make a new conda environment for RMG and its dependencies. Type the following command into the terminal to create the new environment named 'rmg_env' containing the latest stable version of the RMG program and its database. ::

    conda create --name rmg_env 'rmg::rmg'

   Whenever you wish to use it you must first activate the environment::

    conda activate rmg_env

   For more information about using conda, please check out the `conda user guide <https://conda.io/projects/conda/en/latest/user-guide/getting-started.html>`_.
  
   To install a specific version of RMG, add the version to the install command::

    conda create --name rmg_33_env 'rmg::rmg==3.3.0'
  
   Not all versions of RMG are available via conda for all platforms. Check the `official RMG conda channel <https://anaconda.org/RMG/rmg/files>`_ to see which are available for download.

#. You may now run an RMG test job. Save the `Minimal Example Input File <https://raw.githubusercontent.com/ReactionMechanismGenerator/RMG-Py/master/examples/rmg/minimal/input.py>`_
   to a local directory.  Use the terminal to run your RMG job inside that folder using the following command ::

    rmg.py input.py

   You will see a line saying ``MODEL GENERATION COMPLETED`` on your terminal if your RMG test job ran successfully.


Updating your binary installation of RMG in Linux or Mac OSX
============================================================

If you had previously installed a binary version of the RMG package and wish to update to a newer version, we suggest creating a new conda environment and installing the updated version there :: 

    conda create --name rmg_xyz_env 'rmg::rmg==x.y.z'
    
It is also possible, though not advisable, to update your existing installation to the latest stable version available on Anaconda by typing the following command on the terminal ::

    source activate rmg_env
    conda update 'rmg::rmg'

Doing this may break any other code present in the conda environment and RMG may not function correctly. If you attempt this update method and face issues, please attempt to install the new version of RMG in a new conda environment before reaching out for assistance.
