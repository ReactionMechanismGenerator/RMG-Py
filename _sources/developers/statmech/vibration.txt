******************************
Vibrational degrees of freedom
******************************

Vibrational motion describes internal motion of atoms in the system with
respect to one another.

Harmonic oscillators
====================

Many vibrational motions are well-described as one-dimensional quantum harmonic
oscillators. The time-independent Schrodinger equation for such an oscillator
is given by

.. math:: -\frac{\hbar^2}{2m} \frac{\partial^2}{\partial x^2} \Psi(x) + \frac{1}{2} m \omega^2 x^2 \Psi(x) = E \Psi(x)

where :math:`m` is the total mass of the particle. The harmonic potential
results in quantized solutions to the above with the following energy levels:

.. math:: E_n = \left( n + \frac{1}{2} \right) \hbar \omega \hspace{2em} n = 0, 1, 2, \ldots

Above we have introduced :math:`n` as the quantum number. Note that, even in
the ground state (:math:`n = 0`), the harmonic oscillator has an energy that is
not zero; this energy is called the *zero-point energy*.

The harmonic oscillator partition function is obtained by summing over
the above energy levels:

.. math:: Q_\mathrm{vib}(T) = \sum_{n=0}^\infty \exp \left( -\frac{\left( n + \frac{1}{2} \right) \hbar \omega}{k_\mathrm{B} T} \right)

This summation can be evaluated explicitly to give a closed-form analytical
expression for the vibrational partition function of a quantum harmonic
oscillator:

.. math:: Q_\mathrm{vib}(T) = \frac{e^{-\hbar \omega / 2 k_\mathrm{B} T}}{1 - e^{-\hbar \omega / k_\mathrm{B} T}}

In RMG the convention is to place the zero-point energy in with the
ground-state energy of the system instead of the numerator of the vibrational
partition function, which gives

.. math:: Q_\mathrm{vib}(T) = \frac{1}{1 - e^{-\hbar \omega / k_\mathrm{B} T}}

The energy levels of the harmonic oscillator in chemical systems are often
significant compared to the temperature of interest, so we usually use the
quantum result. However, the classical limit is provided here for completeness:

.. math:: Q_\mathrm{vib}^\mathrm{cl}(T) = \frac{k_\mathrm{B} T}{\hbar \omega}

In RMG, the vibrational motion of one or more harmonic oscillators is
represented using the 
:class:`HarmonicOscillator <rmgpy.statmech.HarmonicOscillator>` class.

.. autoclass:: rmgpy.statmech.HarmonicOscillator
    :members:
