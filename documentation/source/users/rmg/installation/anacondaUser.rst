.. _anacondaUser:

****************************************************************************
Binary Installation Using Anaconda for Unix-Based Systems: Linux and Mac OSX
****************************************************************************

* Download and install the `Anaconda Python Platform <http://continuum.io/downloads>`_ for Python 2.7 (make sure not to install Python 3.0+, which is incompatible with RMG). When prompted to append Anaconda to your PATH, select or type Yes.

* Install both RMG and the RMG-database binaries.  Dependencies will be installed automatically.  It is safest to make a new Anaconda environment for RMG. For example, type the following into your terminal to create a new environment named 'rmg_env' containing all that you need::

    conda create -c rmg --name rmg_env rmg rmgdatabase
    
  Whenever you wish to use it you must first activate the environment::
    
    source activate rmg_env
    
* Optional: If you wish to use the :ref:`QMTP interface <qm>` with `MOPAC <http://openmopac.net/>`_ to run quantum mechanical calculations for improved thermochemistry estimates of cyclic species, please obtain a legal license through the `MOPAC License Request Form <http://openmopac.net/form.php>`_.  Once you have it, type the following into your Terminal ::
    
    mopac password_string_here

* You may now run an RMG test job. Save the `Minimal Example Input File <https://raw.githubusercontent.com/ReactionMechanismGenerator/RMG-Py/master/examples/rmg/minimal/input.py>`_  
  to a local directory.  Use the Terminal to run your RMG job inside that folder using the following command ::

    rmg.py input.py

You may now use RMG-Py, CanTherm, as well as any of the 
:ref:`Standalone Modules <modules>` included in the RMG-Py package.
