.. _solverapireference:

************************************************
Reaction system simulation (:mod:`rmgpy.solver`)
************************************************

.. module:: rmgpy.solver

The :mod:`rmgpy.solver` module contains classes used to represent and simulate
reaction systems.



Reaction systems
================

.. currentmodule:: rmgpy.solver

=========================== ====================================================
Class                       Description
=========================== ====================================================
:class:`ReactionSystem`     Base class for all reaction systems
:class:`SimpleReactor`      A simple isothermal, isobaric, well-mixed batch reactor
=========================== ====================================================



Termination criteria
====================

.. currentmodule:: rmgpy.solver

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`TerminationTime`        Represent a time at which the simulation should be terminated
:class:`TerminationConversion`  Represent a species conversion at which the simulation should be terminated
=============================== ================================================

.. toctree::
    :hidden:
    
    reactionsystem
    simplereactor
    termination

