********************************
Kinetics (:mod:`rmgpy.kinetics`)
********************************

.. module:: rmgpy.kinetics

The :mod:`rmgpy.kinetics` subpackage contains classes that represent various
kinetics models of chemical reaction rates and models of quantum mechanical
tunneling through an activation barrier.



Pressure-independent kinetics models
====================================

.. currentmodule:: rmgpy.kinetics

======================= ========================================================
Class                   Description
======================= ========================================================
:class:`KineticsData`   A kinetics model based on a set of discrete rate coefficient points in temperature
:class:`Arrhenius`      A kinetics model based on the (modified) Arrhenius expression
:class:`MultiArrhenius` A kinetics model based on a sum of :class:`Arrhenius` expressions
======================= ========================================================



Pressure-dependent kinetics models
==================================

.. currentmodule:: rmgpy.kinetics

=========================== ====================================================
Class                       Description
=========================== ====================================================
:class:`PDepKineticsData`   A kinetics model based on a set of discrete rate coefficient points in temperature and pressure
:class:`PDepArrhenius`      A kinetics model based on a set of Arrhenius expressions for a range of pressures
:class:`MultiPDepArrhenius` A kinetics model based on a sum of :class:`PDepArrhenius` expressions
:class:`Chebyshev`          A kinetics model based on a Chebyshev polynomial representation
:class:`ThirdBody`          A low pressure-limit kinetics model based on the (modified) Arrhenius expression, with a third body 
:class:`Lindemann`          A kinetics model of pressure-dependent falloff based on the Lindemann model
:class:`Troe`               A kinetics model of pressure-dependent falloff based on the Lindemann model with the Troe falloff factor
=========================== ====================================================



Tunneling models
================

.. currentmodule:: rmgpy.kinetics

======================= ========================================================
Class                   Description
======================= ========================================================
:class:`Wigner`         A one-dimensional tunneling model based on the Wigner expression
:class:`Eckart`         A one-dimensional tunneling model based on the (asymmetric) Eckart expression
======================= ========================================================


.. toctree::
    :hidden:
    
    kineticsdata
    arrhenius
    multiarrhenius
    pdepkineticsdata
    pdeparrhenius
    multipdeparrhenius
    chebyshev
    thirdbody
    lindemann
    troe
    wigner
    eckart
