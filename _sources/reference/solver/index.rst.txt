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
:class:`LiquidReactor`      A homogeneous, isothermal, isobaric liquid batch reactor
:class:`SurfaceReactor`     A heterogeneous, isothermal, isochoric batch reactor
:class:`MBSampledReactor`   :class:`SimpleReactor` with sampling delay for simulating molecular beam experiments
=========================== ====================================================



Termination criteria
====================

.. currentmodule:: rmgpy.solver

=============================== ================================================
Class                           Description
=============================== ================================================
:class:`TerminationTime`        Represent a time at which the simulation should be terminated
:class:`TerminationConversion`  Represent a species conversion at which the simulation should be terminated
:class:`TerminationRateRatio`   Represent a fraction of the maximum characteristic rate at which the simulation should be terminated
=============================== ================================================

.. toctree::
    :hidden:
    
    reactionsystem
    simplereactor
    liquidreactor
    surfacereactor
    mbsampledreactor
    termination

