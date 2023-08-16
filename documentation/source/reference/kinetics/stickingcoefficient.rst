**********************************
rmgpy.kinetics.StickingCoefficient
**********************************

.. autoclass:: rmgpy.kinetics.StickingCoefficient

    The Sticking coefficient class uses an Arrhenius-like expression to determine
    the sticking coefficient :math:`\gamma` of a species on a surface.

    .. math:: \gamma = A (T/T_0)^{n} \exp\left(-\frac{E_\mathrm{a}}{RT}\right)

    Above, :math:`A` is the preexponential factor, :math:`T_0` is the reference
    temperature, :math:`n` is the temperature exponent, and :math:`E_\mathrm{a}`
    is the activation energy.

    The sticking coefficient is then used to calculate the rate coefficient :math:`k`
    using the following expression:

    .. math:: k = \frac{\gamma}{\Gamma} \sqrt{\frac{RT}{2 \pi W}}

    where :math:`\Gamma` is the surface site density and :math:`W` is the molecular 
    weight of the gas phase species. 
