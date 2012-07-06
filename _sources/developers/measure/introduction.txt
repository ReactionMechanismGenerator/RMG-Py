************
Introduction
************

The following is a reference to the application programming interface (API) 
used within MEASURE. If you are a developer wishing to know the full 
implementation details of MEASURE - e.g. if you are trying to modify the 
codebase - this document is for you. If you are simply trying to learn how to
generate phenomenological rate coefficients with MEASURE, you should switch to
the :ref:`measureusersguide`. If you want to know the mathematical basis for the
master equation methods used by measure, you should read the 
:ref:`measuretheoryguide`. (In fact, if you haven't already, you may wish to 
read these documents first before continuing.)

The source for MEASURE is split across a number of Python and Cython modules. 
Each section of the developers' guide references a different module:

=============================== ================================================
Module                          Description
=============================== ================================================
:mod:`rmgpy.measure.input`      Reading network data from input files
:mod:`rmgpy.measure.network`    Representing network data in memory
:mod:`rmgpy.measure.collision`  Working with collision models and parameters
:mod:`rmgpy.measure.reaction`   Working with microcanonical and phenomenological rate coefficients
:mod:`rmgpy.measure.msc`        Determining :math:`k(T,P)` via the modified strong collision method
:mod:`rmgpy.measure.rs`         Determining :math:`k(T,P)` via the reservoir state method
:mod:`rmgpy.measure.output`     Writing network data to input and output files
:mod:`rmgpy.measure.settings`   Setting application-wide options
=============================== ================================================
