************
Introduction
************

**Arkane** (Automated Reaction Kinetics and Network Exploration) is a tool for computing the thermodynamic properties
of chemical species and high-pressure-limit rate coefficients for chemical reactions using
the results of a quantum chemistry calculation. Thermodynamic properties are
computed using the rigid rotor-harmonic oscillator approximation with optional
corrections for hindered internal rotors. Kinetic parameters are computed using
canonical transition state theory with optional tunneling correction.

Arkane can also estimate
pressure-dependent phenomenological rate coefficients :math:`k(T,P)` for
unimolecular reaction networks of arbitrary complexity. The approach is to
first generate a detailed model of the reaction network using the
one-dimensional master equation, then apply one of several available model
reduction methods of varying accuracy, speed, and robustness to simplify the
detailed model into a set of phenomenological rate coefficients. The result
is a set of :math:`k(T,P)` functions suitable for use in chemical reaction
mechanisms. More information is available at `Allen et al. <http://dx.doi.org/10.1039/c1cp22765c>`_.

Arkane is developed and distributed as part of `RMG-Py <http://rmg.mit.edu/>`_, but can be used as a stand-alone
application for Thermochemistry, Transition State Theory, and Master Equation chemical kinetics calculations.

Arkane is written in the `Python <http://www.python.org/>`_ programming
language to facilitate ease of development, installation, and use.

Additional theoretical background can be found at `RMG's Theory Guide <http://reactionmechanismgenerator.github.io/RMG-Py/theory/index.html>`_
and `Arkane's Manual <manual.pdf>`_ as well as the `manual's supplement information <manual_supplement-Solving1DSchrodingerEquation.pdf>`_.

License
=======

Arkane is provided as free, open source code under the terms of the
`MIT/X11 License <http://www.opensource.org/licenses/mit-license.php>`_. The
full, official license is reproduced below



.. literalinclude:: ../../../../LICENSE.txt
