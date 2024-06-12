.. _anacondaUser:

****************************************************************************
Binary Installation Using Anaconda for Unix-Based Systems: Linux and Mac OSX
****************************************************************************


#. Install the ``conda`` package manager via `Miniconda <https://docs.anaconda.com/free/miniconda/miniconda-install/>`_.

#. Create and activate a new ``conda`` environment::

    conda create --name rmg_env python=3.7
    conda activate rmg_env

#. Install the latest RMG Python package::

     conda install rmg --channel rmg
    
    or a specific version of ``rmg`` and/or ``rmgdatabase`` by specifying the version number::

     conda install rmg=3.1.0 rmgdatabase=3.1.1 --channel rmg

#. Install RMS::

    conda install julia=1.9.1 pyjulia>=0.6 --channel conda-forge
    conda install pyrms diffeqpy --channel rmg
    python -c "import julia; julia.install(); import diffeqpy; diffeqpy.install()"
    julia -e 'using Pkg; Pkg.add(PackageSpec(name="ReactionMechanismSimulator",rev="main")); using ReactionMechanismSimulator'

.. NOTE::
    Arkane and RMG can both be run without installing RMS, though with a reduced feature set.
    Arkane will use the slower ``scipy`` solvers instead of ``diffeqpy`` when solving the master equation
    and RMG will lose some advanced reactor types and the sensitivity analysis feature.

#. You may now run an RMG test job. Save the `Minimal Example Input File <https://raw.githubusercontent.com/ReactionMechanismGenerator/RMG-Py/master/examples/rmg/minimal/input.py>`_
   to a local directory.  Use the terminal to run your RMG job inside that folder using the following command ::

    rmg.py input.py

   To run models requiring RMS, prefix the command with ``python-jl``::

    python-jl rmg.py input.py

   Whenever you wish to use RMG you must first activate the environment::

    conda activate rmg_env

   If you encounter the ImportError related to ``libmkl_rt.so.2``, refer to the :ref:`Fixing the ImportError related to libmkl_rt.so.2 <fixImportError>`
   section below to fix the error and re-run the RMG test job.

You may now use RMG-Py, Arkane, as well as any of the :ref:`Standalone Modules <modules>` included in the RMG-Py package.

.. _fixImportError:

Fixing the ImportError related to ``libmkl_rt.so.2``
============================================================

You may encounter the following ImportError when you try to run a RMG test job ::

    Traceback (most recent call last):
      File "/PATH-TO-YOUR-ANACONDA/envs/rmg_env/bin/rmg.py", line 48, in <module>
        from rmgpy.rmg.main import RMG, initialize_log, process_profile_stats, make_profile_graph
      File "/PATH-TO-YOUR-ANACONDA/envs/rmg_env/lib/python3.7/site-packages/rmgpy/rmg/main.py", line 51, in <module>
        from cantera import ck2cti
      File "/PATH-TO-YOUR-ANACONDA/envs/rmg_env/lib/python3.7/site-packages/cantera/__init__.py", line 4, in <module>
        from ._cantera import *
    ImportError: libmkl_rt.so.2: cannot open shared object file: No such file or directory

where ``PATH-TO-YOUR-ANACONDA`` is the path to the ``Anaconda3`` directory installed on your machine.
The default install location of Anaconda on Linux is ``/home/<your_username>/Anaconda3``.

To fix this issue, simply copy ``libmkl_rt.so.1`` as ``libmkl_rt.so.2`` under the rmg environment library by typing the following
line on your terminal::

    cp /PATH-TO-YOUR-ANACONDA/envs/rmg_env/lib/libmkl_rt.so.1 /PATH-TO-YOUR-ANACONDA/envs/rmg_env/lib/libmkl_rt.so.2

Note that ``PATH-TO-YOUR-ANACONDA`` needs to be the path to your ``Anaconda3`` directory, which mostly likely looks
like ``/home/<your_username>/Anaconda3``. If you are unable to locate the ``libmkl_rt.so.1`` file on your computer, you can find its location with the following command::

    locate libmkl_rt.so.1

After copying the file as ``libmkl_rt.so.2``, try running the RMG test job again::

    rmg.py input.py

You will see a line saying ``MODEL GENERATION COMPLETED`` on your terminal if your RMG test job ran successfully.


Updating your binary installation of RMG in Linux or Mac OSX
============================================================

If you had previously installed a binary version of the RMG package you may
check and update your installation to the latest stable version available on Anaconda by typing the following command on the terminal ::

    conda activate rmg_env
    conda update rmg --channel rmg
