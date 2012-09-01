********************************
Translational degrees of freedom
********************************

Translational motion describes external motion of the system as a whole, and
represents three degrees of freedom (one per dimension of motion).

Ideal gas translation
=====================

The translation of an *ideal gas* -- a gas composed of randomly-moving,
noninteracting particles of negligible size -- in three dimensions can be 
modeled using the particle-in-a-box model. In this model, a gas particle is
confined to a three-dimensional box of size :math:`L_x L_y L_z = V` with the
following potential:

.. math:: 

    V(x,y,z) = \begin{cases}
    0 & 0 \le x \le L_x, 0 \le y \le L_y, 0 \le z \le L_z \\
    \infty & \mathrm{otherwise}
    \end{cases}

The time-independent Schrodinger equation for this system (within the box) is
given by

.. math:: -\frac{\hbar^2}{2M} \left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} + \frac{\partial^2}{\partial z^2} \right) \Psi(x,y,z) = E \Psi(x,y,z)

where :math:`M` is the total mass of the particle. Because the box is finite in
all dimensions, the solution of the above is quantized with the following
energy levels:

.. math:: E_{n_x,n_y,n_z} = \frac{\hbar^2}{2M} \left[ \left( \frac{n_x \pi}{L_x} \right)^2 + \left( \frac{n_y \pi}{L_y} \right)^2 + \left( \frac{n_z \pi}{L_z} \right)^2 \right] \hspace{2em} n_x, n_y, n_z = 1, 2, \ldots

Above we have introduced :math:`n_x`, :math:`n_y`, and :math:`n_z` as quantum
numbers. The quantum mechanical partition function is obtained by summing over
the above energy levels:

.. math:: Q_\mathrm{trans}(T) = \sum_{n_x=1}^\infty \sum_{n_y=1}^\infty \sum_{n_z=1}^\infty \exp \left( -\frac{E_{n_x,n_y,n_z}}{k_\mathrm{B} T} \right)

In almost all cases the temperature of interest is large relative to the
energy spacing; in this limit we can obtain a closed-form analytical expression
for the translational partition function in the classical limit:

.. math:: Q_\mathrm{trans}^\mathrm{cl}(T) = \left( \frac{2 \pi M k_\mathrm{B} T}{h^2} \right)^{3/2} V

For a constant-pressure problem we can use the ideal gas law to replace 
:math:`V` with :math:`k_\mathrm{B} T/P`. This gives the partition function a
temperature dependence of :math:`T^{5/2}`.

In RMG, the three-dimensional translational motion of an ideal gas is
represented using the :class:`IdealGasTranslation <rmgpy.statmech.IdealGasTranslation>` class.

.. autoclass:: rmgpy.statmech.IdealGasTranslation
    :members:
