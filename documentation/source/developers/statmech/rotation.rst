*****************************
Rotational degrees of freedom
*****************************

Rotational motion describes external rotation of the system as a whole. For
most systems it represents three degrees of freedom (one per dimension of
motion). However, linear systems have one less degree of rotational freedom,
as rotation about the linear axis of the system does not result in any movement.

Linear rotors
=============

A linear rigid rotor is modeled as a pair of point masses :math:`m_1` and 
:math:`m_2` separated by a distance :math:`R`. Since we are modeling the
rotation of this system, we choose to work in spherical coordinates. Following 
the physics convention -- where :math:`0 \le \theta \le \pi` is the zenith
angle and :math:`0 \le \phi \le 2\pi` is the azimuth -- the Schrodinger
equation for the rotor is given by

.. math:: -\frac{\hbar^2}{2I} \left[ \frac{1}{\sin \theta} \frac{\partial}{\partial \theta} \left( \sin \theta \frac{\partial}{\partial \theta} \right) + \frac{1}{\sin^2 \theta} \frac{\partial^2}{\partial \phi^2} \right] \Psi(\theta, \phi) &= E \Psi(\theta, \phi)

where :math:`I \equiv \mu R^2` is the moment of inertia of the rotating body,
and :math:`\mu \equiv m_1 m_2 / (m_1 + m_2)` is the reduced mass. Note that
there is no potential term in the above expression; for this reason, a rigid
rotor is often referred to as a *free* rotor. Solving the Schrodinger equation
gives the energy levels :math:`E_J` and corresponding degeneracies :math:`g_J`
for the linear rigid rotor as

.. math::

    E_J &= B J (J + 1) \hspace{2em} J = 0, 1, 2, \ldots \\
    g_J &= 2J + 1
    
where :math:`J` is the quantum number for the rotor -- sometimes called the
total angular momentum quantum number -- and :math:`B \equiv \hbar^2/2I` is the
rotational constant.

Using these expressions for the energy levels and corresponding degeneracies,
we can evaluate the partition function for the linear rigid rotor:

.. math::

    Q_\mathrm{rot}(T) = \frac{1}{\sigma} \sum_{J=0}^\infty (2J + 1) e^{- B J (J + 1) / k_\mathrm{B} T}

In many cases the temperature of interest is large relative to the energy
spacing; in this limit we can obtain a closed-form analytical expression
for the linear rotor partition function in the classical limit:

.. math:: Q_\mathrm{rot}^\mathrm{cl}(T) = \frac{1}{\sigma} \frac{8 \pi^2 I k_\mathrm{B} T}{h^2}

Above we have also introduced :math:`\sigma` as the symmetry number of the
rigid rotor.

In RMG, the rotational motion of a linear system represented using the 
:class:`LinearRotor <rmgpy.statmech.LinearRotor>` class.

.. autoclass:: rmgpy.statmech.LinearRotor
    :members:

Nonlinear rotors
================

A nonlinear rigid rotor is the generalization of the linear rotor to a 
nonlinear polyatomic system. Such a system is characterized by three moments
of inertia :math:`I_\mathrm{A}`, :math:`I_\mathrm{B}`, and :math:`I_\mathrm{C}` instead of just one.
The solution to the Schrodinger equation for the quantum nonlinear rotor is not 
well defined, so we will simply show the classical result instead:

.. math:: Q_\mathrm{rot}^\mathrm{cl}(T) = \frac{\pi^{1/2}}{\sigma} \left( \frac{8 k_\mathrm{B} T}{h^2} \right)^{3/2} \sqrt{I_\mathrm{A} I_\mathrm{B} I_\mathrm{C}}

In RMG, the rotational motion of a nonlinear system represented using the 
:class:`NonlinearRotor <rmgpy.statmech.NonlinearRotor>` class.

.. autoclass:: rmgpy.statmech.NonlinearRotor
    :members:

K-rotors
========

A K-rotor is simply the one-dimensional equivalent of a linear rigid rotor.
The energy levels :math:`E_K` of the K-rotor are given by

.. math::

    E_K &= B K^2 \hspace{2em} K = 0, \pm 1, \pm 2, \ldots
    
where :math:`K` is the quantum number for the rotor and 
:math:`B \equiv \hbar^2/2I` is the rotational constant.

Using these expressions for the energy levels and corresponding degeneracies,
we can evaluate the partition function for the K-rotor:

.. math::

    Q_\mathrm{rot}(T) = \frac{1}{\sigma} \left( 1 + \sum_{K=1}^\infty 2 e^{- B K^2 / k_\mathrm{B} T} \right)

In many cases the temperature of interest is large relative to the energy
spacing; in this limit we can obtain a closed-form analytical expression
for the linear rotor partition function in the classical limit:

.. math:: Q_\mathrm{rot}^\mathrm{cl}(T) = \frac{1}{\sigma} \left( \frac{8 \pi^2 I k_\mathrm{B} T}{h^2} \right)^{1/2}

where :math:`\sigma` is the symmetry number of the K-rotor.

In RMG, the rotational motion of a K-rotor is represented using the 
:class:`KRotor <rmgpy.statmech.KRotor>` class.

.. autoclass:: rmgpy.statmech.KRotor
    :members:

Spherical top rotors
====================

A spherical top rotor is simply the three-dimensional equivalent of a linear 
rigid rotor. Unlike the nonlinear rotor, all three moments of inertia of a
spherical top are equal, i.e. :math:`I_\mathrm{A} = I_\mathrm{B} = I_\mathrm{C} = I`.
The energy levels :math:`E_J` and corresponding degeneracies :math:`g_J` of the
spherial top rotor are given by

.. math::

    E_J &= B J (J + 1) \hspace{2em} J = 0, 1, 2, \ldots \\
    g_J &= (2J+1)^2
    
where :math:`J` is the quantum number for the rotor and 
:math:`B \equiv \hbar^2/2I` is the rotational constant.

Using these expressions for the energy levels and corresponding degeneracies,
we can evaluate the partition function for the spherical top rotor:

.. math::

    Q_\mathrm{rot}(T) = \frac{1}{\sigma} \sum_{J=0}^\infty (2J+1)^2 e^{- B J (J+1) / k_\mathrm{B} T}

In many cases the temperature of interest is large relative to the energy
spacing; in this limit we can obtain a closed-form analytical expression
for the linear rotor partition function in the classical limit:

.. math:: Q_\mathrm{rot}^\mathrm{cl}(T) = \frac{1}{\sigma} \left( \frac{8 \pi^2 I k_\mathrm{B} T}{h^2} \right)^{3/2}

where :math:`\sigma` is the symmetry number of the spherical top. Note that
the above differs from the nonlinear rotor partition function by a factor of
:math:`\pi`.

In RMG, the rotational motion of a spherical top rotor is represented using
the :class:`SphericalTopRotor <rmgpy.statmech.SphericalTopRotor>` class.

.. autoclass:: rmgpy.statmech.SphericalTopRotor
    :members:
