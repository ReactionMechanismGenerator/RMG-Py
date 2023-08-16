*******************************
rmgpy.kinetics.SurfaceArrhenius
*******************************

.. autoclass:: rmgpy.kinetics.SurfaceArrhenius

    The SurfaceArrhenius class uses a rate expression that is identical to the 
    Arrhenius class:

    .. math:: k(T) = A \left( \frac{T}{T_0} \right)^n \exp \left( -\frac{E_\mathrm{a}}{RT} \right)

    Above, :math:`A` is the preexponential factor, :math:`T_0` is the reference
    temperature, :math:`n` is the temperature exponent, and :math:`E_\mathrm{a}`
    is the activation energy.

    The only difference is the units for :math:`A` are in terms of surface concentration 
    (:math:`\mathrm{mol/m^2}`) instead of volume concentration (:math:`\mathrm{mol/m^3}`).
