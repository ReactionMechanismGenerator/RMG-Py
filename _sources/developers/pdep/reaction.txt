***************
Reaction events
***************

.. automodule:: rmgpy.pdep.reaction

Microcanonical rate coefficients
================================

.. autofunction:: rmgpy.pdep.calculateMicrocanonicalRateCoefficient

RRKM theory
===========

RRKM (Rice-Ramsperger-Kassel-Marcus) theory is the microcanonical analogue of
transition state theory. The microcanonical rate coefficient as a function of
total energy :math:`E` and total angular momentum quantum number :math:`J` is
given by

.. math:: k(E,J) = \frac{N^\ddagger(E, J)}{h \rho(E, J)}

where :math:`N^\ddagger(E, J)` is the sum of states of the transition state and
:math:`\rho(E, J)` is the density of states of the reactant. If the J-rotor is
treated as active, the J-dependence can be averaged in the above expression to
give

.. math:: k(E) = \frac{N^\ddagger(E)}{h \rho(E)}

as a function of total energy alone. This is reasonable at high temperatures,
but less accurate at low temperatures.

Use of RRKM theory requires detailed information about the statistical mechanics
of the reactant *and* transition state. However, it is generally more accurate
than the inverse Laplace transform method.

.. autofunction:: rmgpy.pdep.applyRRKMTheory

Inverse Laplace transform method
================================

The inverse Laplace transform method exploits the following relationship to
determine the microcanonical rate coefficient:

.. math:: \mathcal{L} \left[k(E) \rho(E) \right] = \int_0^\infty k(E) \rho(E) e^{-E/k_\mathrm{B} T} \ dE = k_\infty(T) Q(T)

Given a high-pressure limit rate coefficient :math:`k_\infty(T)` represented as
an Arrhenius expression with positive :math:`n` and :math:`E_\mathrm{a}`,
the microcanonical rate coefficient :math:`k(E)` can be determined via an
inverse Laplace transform. For :math:`n = 0` the transform can be defined
analytically:

.. math:: k(E) = A \frac{\rho(E - E_\mathrm{a})}{\rho(E)} \ \ \ \ (n = 0)

For :math:`n > 0` the transform is defined numerically. For :math:`n < 0` or
:math:`E_\mathrm{a} < 0` the transform is not defined; in this case we 
approximate by simply lumping the :math:`T^n` or :math:`e^{-E_\mathrm{a}/RT}`
terms into the preexponential factor, and use a different :math:`k(E)` at each
temperature.

The ILT method does not required detailed transition state information, but only
the high-pressure limit kinetics. However, it assumes that (1) 
:math:`k_\infty(T)` is valid over the temperature range from zero to infinity 
and (2) the activation energy :math:`E_\mathrm{a}` is physically identical to 
the reaction barrier :math:`E_0^\ddagger - E_0`.

.. autofunction:: rmgpy.pdep.applyInverseLaplaceTransformMethod
