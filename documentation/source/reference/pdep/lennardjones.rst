***********************
rmgpy.pdep.LennardJones
***********************

.. autoclass:: rmgpy.pdep.LennardJones

    The Lennard-Jones collision model is based on the empirical potential function

    .. math:: V(r) = 4 \epsilon \left[ \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^{6} \right]

    where :math:`\epsilon` is the depth of the potential well and :math:`\sigma` is
    related to the distance at which the potential is minimized by 
    :math:`r_\mathrm{m} = 2^{1/6} \sigma`. The collision frequency in the
    Lennard-Jones model is given by

    .. math:: \omega(T,[\mathrm{M}]) = \Omega_i^{(2,2)\ast} \sqrt{\frac{8 k_\mathrm{B} T}{\pi \mu}} \pi \sigma^2 [\mathrm{M}]

    where :math:`T` is temperature, :math:`[\mathrm{M}]` is the collider
    concentration, and :math:`\mu` is the reduced mass. The configurational 
    integral :math:`\Omega^{(2,2)*}` is well-approximated by a simple algebraic
    equation in temperature

    .. math:: \Omega_i^{(2,2)\ast} = 1.16145 \tilde{T}^{-0.14874} + 0.52487 e^{-0.7732 \tilde{T}} + 2.16178 e^{-2.437887 \tilde{T}}
            
    where :math:`\tilde{T} \equiv k_\mathrm{B} T/\epsilon`. Often the Lennard-Jones 
    parameters of the two colliders are averaged together when determining the
    collision frequency.
