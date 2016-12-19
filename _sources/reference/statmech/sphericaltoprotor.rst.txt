********************************
rmgpy.statmech.SphericalTopRotor
********************************

.. autoclass:: rmgpy.statmech.SphericalTopRotor

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
