**************************
Frequently Asked Questions
**************************

Introduction
============

We have compiled some common questions about installing and using RMG below.
For any other questions related to RMG and its usage and installation, please
post an issue on our `GitHub issues page <https://github.com/ReactionMechanismGenerator/RMG-Py/issues>`_,
where you can also search for any previous reports of your issue.
Alternatively, you can contact us directly at rmg_dev@mit.edu.


Installing RMG
==============

#. **Why does RMG-Py not work natively from source on Windows?**

   One major challenge with supporting Windows is ensuring that all of our dependencies support Windows. This becomes
   non-trivial as we add more dependencies to support increasing RMG functionality. Ensuring that code within RMG is
   platform-agnostic is also challenging, since it is rarely the first priority for new development because our main
   focus is on research.


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
