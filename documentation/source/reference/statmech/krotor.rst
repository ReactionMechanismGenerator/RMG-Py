*********************
rmgpy.statmech.KRotor
*********************

.. autoclass:: rmgpy.statmech.KRotor

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
