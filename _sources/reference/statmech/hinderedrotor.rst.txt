****************************
Torsional degrees of freedom
****************************

.. autoclass:: rmgpy.statmech.HinderedRotor

    The Schrodinger equation for a one-dimensional hindered rotor is given by

    .. math:: -\frac{\hbar^2}{2I} \frac{d^2}{d \phi^2} \Psi(\phi) + V(\phi) \Psi(\phi) &= E \Psi(\phi)

    where :math:`I` is the reduced moment of inertia of the torsion and 
    :math:`V(\phi)` describes the potential of the torsion. There are two common
    forms for the potential: a simple cosine of the form

    .. math:: V(\phi) = \frac{1}{2} V_0 \left( 1 - \cos \sigma \phi \right)

    where :math:`V_0` is the barrier height and :math:`\sigma` is the symmetry
    number, or a more general Fourier series of the form

    .. math:: V(\phi) = A + \sum_{k=1}^C \left( a_k \cos k \phi + b_k \sin k \phi \right)

    where :math:`A`, :math:`a_k` and :math:`b_k` are fitted coefficients. Both
    potentials are typically defined such that the minimum of the potential is zero
    and is found at :math:`\phi = 0`.

    For either the cosine or Fourier series potentials, the energy levels of the 
    quantum hindered rotor must be determined numerically. The cosine potential 
    does permit a closed-form representation of the classical partition function,
    however:

    .. math:: Q_\mathrm{hind}^\mathrm{cl}(T) = \left( \frac{2 \pi I k_\mathrm{B} T}{h^2} \right)^{1/2} \frac{2 \pi}{\sigma} \exp \left( -\frac{V_0}{2 k_\mathrm{B} T} \right) I_0 \left( \frac{V_0}{2 k_\mathrm{B} T} \right)

    A semiclassical correction to the above is usually required to provide a
    reasonable estiamate of the partition function:

    .. math:: 

        Q_\mathrm{hind}^\mathrm{semi}(T) &= \frac{Q_\mathrm{vib}^\mathrm{quant}(T)}{Q_\mathrm{vib}^\mathrm{cl}(T)} Q_\mathrm{hind}^\mathrm{cl}(T) \\
                                        &= \frac{h \nu}{k_\mathrm{B} T} \frac{1}{1 - \exp \left(- h \nu / k_\mathrm{B} T \right)} \left( \frac{2 \pi I k_\mathrm{B} T}{h^2} \right)^{1/2} \frac{2 \pi}{\sigma} \exp \left( -\frac{V_0}{2 k_\mathrm{B} T} \right) I_0 \left( \frac{V_0}{2 k_\mathrm{B} T} \right)

    Above we have defined :math:`\nu` as the vibrational frequency of the hindered
    rotor:

    .. math:: \nu \equiv \frac{\sigma}{2 \pi} \sqrt{\frac{V_0}{2 I}}

