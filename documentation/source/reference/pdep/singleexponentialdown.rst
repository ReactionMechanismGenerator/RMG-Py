********************************
rmgpy.pdep.SingleExponentialDown
********************************

.. autoclass:: rmgpy.pdep.SingleExponentialDown

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
