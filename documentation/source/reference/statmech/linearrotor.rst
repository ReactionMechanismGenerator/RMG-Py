**************************
rmgpy.statmech.LinearRotor
**************************

.. autoclass:: rmgpy.statmech.LinearRotor

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
