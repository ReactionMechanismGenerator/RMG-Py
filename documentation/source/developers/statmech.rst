******************************************************
:mod:`rmgpy.statmech` --- Molecular Degrees of Freedom
******************************************************

.. automodule:: rmgpy.statmech

Working With Molecular Degrees of Freedom
=========================================

Each atom in a molecular configuration has three spatial dimensions in which it
can move. Thus, a molecular configuration consisting of :math:`N` atoms has
:math:`3N` degrees of freedom. We can distinguish between those modes that
involve movement of atoms relative to the molecular center of mass (called
*internal* modes) and those that do not (called *external* modes). Of the
external degrees of freedom, three involve translation of the entire molecular
configuration, while either three (for a nonlinear molecule) or two (for a
linear molecule) involve rotation of the entire molecular configuration
around the center of mass. The remaining :math:`3N-6` (nonlinear) or
:math:`3N-5` (linear) degrees of freedom are the internal modes, and can be
divided into those that involve vibrational motions (symmetric and asymmetric
stretches, bends, etc.) and those that involve torsional rotation around single
bonds between nonterminal heavy atoms.

The mathematical description of these degrees of freedom falls under the purview
of quantum chemistry, and involves the solution of the time-independent
Schrodinger equation:

    .. math:: \hat{H} \psi = E \psi

where :math:`\hat{H}` is the Hamiltonian, :math:`\hat{H}` is the wavefunction,
and :math:`E` is the energy. The exact form of the Hamiltonian varies depending
on the degree of freedom you are modeling. Since this is a quantum system, the
energy can only take on discrete values. Once the allowed energy levels are
known, the partition function :math:`Q(\beta)` can be computed using the
summation

    .. math:: Q(\beta) = \sum_i g_i e^{-\beta E_i}

where :math:`g_i` is the degeneracy of energy level :math:`i` (i.e. the number
of energy states at that energy level) and
:math:`\beta \equiv (k_\mathrm{B} T)^{-1}`.

The partition function is an immensely useful quantity, as all sorts of
thermodynamic parameters can be evaluated using the partition function:

    .. math:: A = - k_\mathrm{B} T \ln Q

    .. math:: U = - \frac{\partial \ln Q}{\partial \beta}

    .. math:: S = \frac{\partial}{\partial T} \left( k_\mathrm{B} T \ln Q \right)

    .. math:: C_\mathrm{v} = \frac{1}{k_\mathrm{B} T} \frac{\partial^2 \ln Q}{\partial \beta^2}

Above, :math:`A`, :math:`U`, :math:`S`, and :math:`C_\mathrm{v}` are the
Helmholtz free energy, internal energy, entropy, and constant-volume heat
capacity, respectively.

The partition function for a molecular configuration is the product of the
partition functions for each invidual degree of freedom:

    .. math:: Q = Q_\mathrm{trans} Q_\mathrm{rot} Q_\mathrm{vib} Q_\mathrm{tors} Q_\mathrm{elec}

This means that the contributions to each thermodynamic quantity from each
molecular degree of freedom are additive.

Representing Molecular Degrees of Freedom
=========================================

.. autoclass:: rmgpy.statmech.StatesModel
    :members:

.. autoclass:: rmgpy.statmech.Mode
    :members:

External Degrees of Freedom
===========================

Translation
-----------

.. autoclass:: rmgpy.statmech.Translation
    :members:

Rotation
--------

.. autoclass:: rmgpy.statmech.RigidRotor
    :members:

Internal Degrees of Freedom
===========================

Vibration
---------

.. autoclass:: rmgpy.statmech.HarmonicOscillator
    :members:

Torsion
-------

.. autoclass:: rmgpy.statmech.HinderedRotor
    :members:
