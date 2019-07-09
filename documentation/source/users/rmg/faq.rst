**************************
Frequently Asked Questions
**************************

Introduction
============

We have compiled some common questions about installing and using RMG below.
For any other questions related to RMG and its usage and installation, please
post an issue on our `GitHub issues page <https://github.com/ReactionMechanismGenerator/RMG-Py/issues>`_,
where you can also search for any previous reports of your issue.
Alternatively, you can also ask questions via the `RMG-Py chat room <https://gitter.im/ReactionMechanismGenerator/RMG-Py>`_
or by contacting us directly at rmg_dev@mit.edu.


Installing RMG
==============

#. **How can I install RMG-Py without Anaconda?**

   Usually we don't recommend installing RMG-Py without Anaconda because it takes longer and is easier to get trouble
   with package management. But one still can try direct installation on Linux or MacOS by following
   :ref:`Linux instruction<linux>` or :ref:`MacOS instruction<macos>`. The RMG team does not use this install approach
   internally any more, so these instructions are not actively maintained.

#. **Why does RMG-Py not work natively on Windows?**

   One major challenge with supporting Windows is ensuring that all of our dependencies support Windows. This becomes
   non-trivial as we add more dependencies to support increasing RMG functionality. Ensuring that code within RMG is
   platform-agnostic is also challenging, since it is rarely the first priority for new development because our main
   focus is on research.

#. **What is the recommended way to run RMG-Py on Windows?**

   The currently recommended way to run RMG on Windows is to set up a Linux environment. There are multiple ways you
   can approach this. Windows 10 supports a Linux subsystem which allows one to set up a Linux environment within
   Windows without using virtualization. You can find instructions on setting up RMG within the Linux subsystem
   :ref:`here<linuxSubsystem>`.

   Another option would be to set up a full Linux virtual machine using something like VirtualBox or VMWare Workstation.
   The benefit of this option is being able to run in a full Linux environment. However, running two operating systems
   simultaneously does result in excess resource overhead, so it may not be suitable for running extended RMG jobs.
   Instructions for setting up a virtual machine can be found :ref:`here<virtualMachineSetup>`.

   A third option that we are currently beginning to explore using `Docker <https://www.docker.com/>`_, which is a
   container-based infrastructure which shares the same benefits as a virtual machine but with less overhead. There
   are some test images of RMG-Py which can be found on `Docker Hub <https://hub.docker.com/>`_ if you would like to
   give this a try. More detailed instructions will be made available once we officially support this approach.

#. **Windows binary installation gives ``WindowsError: [Error 5]``?**

   Error 5 is access is denied, so this is either a permissions error, or an issue with the Windows file lock.
   `These posts <https://github.com/conda/conda/issues/708>`_ suggest rebooting the computer (in case it's a file lock),
   and running the anaconda prompt, from which you run ``conda create -c rmg --name rmg_env rmg rmgdatabase``,
   as an administrator (in case it's a permissions error). Please checkout one example from a user having
   `Windows binary installation issue <https://github.com/ReactionMechanismGenerator/RMG-Py/issues/779>`_.


Running RMG
===========

#. **How do I run a basic RMG job?**

   Please see step-by-step instructions in the either the :ref:`binary<anacondaUser>` or :ref:`source <anacondaDeveloper>`
   installation instructions. In general, the syntax is ::

    rmg.py input.py

   if the RMG-Py directory has been properly added to PATH. For the binary installation, this is done automatically,
   but you will need to do this yourself if installing from source.

   For information on writing an RMG input file, please see the :ref:`documentation<input>`.

#. **Why did I get** ``ImportError: No module named graph`` **when trying to run RMG?**

   This error is commonly seen when RMG is run before compiling. The ``graph`` module just happens to be one of the
   first Cython modules to be imported. To resolve this, compile RMG using the ``make`` command while in the main
   RMG-Py directory.

#. **Why did I get** ``ImportError: No module named cantera`` **when trying to run RMG?**

   This error is commonly seen when RMG is run without properly setting up the Anaconda environment. The ``cantera``
   module happens to be one of the first dependencies to be imported. To resolve this, activate the RMG environment
   using ``conda activate rmg_env``. If the environment has not already been created, create the environment according
   to the :ref:`Linux installation instructions<anacondaDeveloper>`.

#. **What is an** ``UndeterminableKineticsError``?

   This is a common cause of crashed RMG jobs. It is often due to inconsistencies in how reaction templates are
   defined in RMG-database, and occasionally due to inconsistencies in resonance structure generation. Unfortunately,
   this type of error can be very difficult to debug. You can post such issues to our
   `GitHub issues page <https://github.com/ReactionMechanismGenerator/RMG-Py/issues>`_.

#. **What is an** ``InvalidMicrocanonicalRateError``?

   This is another common cause of crashed RMG jobs when using the pressure dependence module. It is due to a failure
   to converge the microcanonical rate calculation for a pressure dependent network. It can be due to a variety of
   factors, such as poor thermochemistry or rate constants. Unfortunately, there is currently no good way to debug and
   fix these types of errors.

#. **Why did I get** ``Segmentation fault:11`` **after installing RMG on my machine?**

   **Segmentation fault** is a typical error in C code, caused by a program trying to read or write an illegal memory
   location, i.e. one it is not allowed to access. The most common cause in RMG is a conflict between two different
   versions of a shared library. RMG has some dependencies which are written in C++, e.g. rdkit, openbabel. If you
   compile one of these with a different version of some compiler library, or you compile RMG using one version and
   run it with another, you will often get a Segmentation fault. Chances are those packages are not up to date, or
   maybe your environmental variable ``PATH`` is messed up so that the wrong version of something is being found.
   Please see one example from a user having same
   `Segmentation fault issue <https://github.com/ReactionMechanismGenerator/RMG-website/issues/125>`_.

#. **Why did I get** ``IOError: [Errno 13] Permission denied: 'C:\\RMG.log'``

   You do not have permission to write to the log file. Try running the RMG from a different folder that you do have
   write permission to, such as within your user's documents directory, or else try running the command prompt as an
   Administrator (so that you have write permission everywhere). See for example
   `issue #817 <https://github.com/ReactionMechanismGenerator/RMG-Py/issues/817>`_.


Miscellaneous
=============

#. **Why can't my adjacency lists be read any more?**

   The adjacency list syntax changed in July 2014.
   The minimal requirement for most translations is to prefix the number
   of unpaired electrons with the letter `u`.

   Example old syntax::

    HXD13
    1    C 0       {2,D}
    2    C 0 {1,D} {3,S}
    3    C 0 {2,S} {4,D}
    4    C 0 {3,D} {5,S}
    5 *1 C 0 {4,S} {6,S}
    6 *2 C 0 {5,S}

   Example new syntax::

    HXD13
    1    C u0       {2,D}
    2    C u0 {1,D} {3,S}
    3    C u0 {2,S} {4,D}
    4    C u0 {3,D} {5,S}
    5 *1 C u0 {4,S} {6,S}
    6 *2 C u0 {5,S}
    
   The new syntax, however, allows much greater flexibility, including definition of lone pairs, partial charges,
   wildcards, and molecule multiplicities, and was necessary to allow us to add Nitrogen chemistry.
   See :ref:`rmgpy.molecule.adjlist` for details of the new syntax.
