****************************
Statistical mechanics in RMG
****************************

.. module:: rmgpy.statmech

The :mod:`rmgpy.statmech` subpackage contains classes that represent various
statistical mechanical models of molecular degrees of freedom. These models
enable the computation of macroscopic parameters (e.g. thermodynamics, kinetics,
etc.) from microscopic parameters.

A molecular system consisting of :math:`N` atoms is described by :math:`3N`
molecular degrees of freedom. Three of these modes involve translation of the
system as a whole. Another three of these modes involve rotation of the system
as a whole, unless the system is linear (e.g. diatomics), for which there are
only two rotational modes. The remaining :math:`3N-6` (or :math:`3N-5` if 
linear) modes involve internal motion of the atoms within the system. Many of
these modes are well-described as harmonic oscillations, while others are
better modeled as torsional rotations around a bond within the system.

Molecular degrees of freedom are mathematically represented using the
Schrodinger equation :math:`\hat{H} \Psi = E \Psi`. By solving the
Schrodinger equation, we can determine the available energy states of the
molecular system, which enables computation of macroscopic parameters. 
Depending on the temperature of interest, some modes (e.g. vibrations) require 
a quantum mechanical treatment, while others (e.g. translation, rotation) can 
be described using a classical solution.

.. toctree::
    :maxdepth: 2
    
    translation
    rotation
    vibration
    torsion
    schrodinger
