****************
Collision events
****************

.. automodule:: rmgpy.pdep.collision

Collision frequency models
==========================

Lennard-Jones
-------------

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

.. autoclass:: rmgpy.pdep.LennardJones
    :members:

Collisional energy transfer models
==================================
    
Collisional energy transfer models seek to represent the probability 
:math:`P(E,E')` of a collision causing a molecule to move from energy 
:math:`E'` to energy :math:`E`. (In principle, we should also allow for a
change of total angular momentum, but very few collision models are
sophisticated enough to consider this.) All expressions for :math:`P(E,E')`
are constrained by two conditions. The first is detailed balance, the
requirement that at equilibrium the rate of energy transfer betweeen all pairs
of states is the same. This is expressed mathematically as

.. math:: P(E, E') \rho(E') e^{-E'/k_\mathrm{B}T} = P(E', E) \rho(E) e^{-E/k_\mathrm{B}T}

where :math:`\rho(E)` is the density of states. The second condition is
normalization, where the total probability of a transfer from :math:`E'` to
*any* :math:`E` must be unity:

.. math:: \int_0^\infty P(E, E') \ dE' = 1
    
Single exponential down
-----------------------

The single exponential down collision model is based around the collisional 
energy transfer probability function for a deactivating collision
    
.. math:: P(E, E') = C(E') \exp \left( - \frac{E' - E}{\alpha} \right) \hspace{40pt} E < E'
    
where the parameter :math:`\alpha = \left< \Delta E_\mathrm{d} \right>`
represents the average energy transferred in a deactivating collision. The
probability function for activating collisions is determined from detailed
balance, and the proportionality constant :math:`C(E')` is determined by the 
normalization condition.

This is the most commonly-used collision model, simply because it only has one
parameter to determine. The parameter :math:`\alpha` is often specified as a
power law in temperature

.. math:: \alpha = \alpha_0 \left( \frac{T}{T_0} \right)^n

where :math:`\alpha_0` is the value of :math:`\alpha` at temperature
:math:`T_0` in K. Set the exponent :math:`n` to zero to obtain a
temperature-independent value for :math:`\alpha`.

.. autoclass:: rmgpy.pdep.SingleExponentialDown
    :members:
