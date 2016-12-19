************************
rmgpy.kinetics.ThirdBody
************************

.. autoclass:: rmgpy.kinetics.ThirdBody

    Third-body kinetics simply introduce an inert third body to the rate expression:

    .. math:: k(T,P) = k_0(T) [\mathrm{M}]

    Above, :math:`[\mathrm{M}] \approx P/RT` is the concentration of the bath gas.
    This formulation is equivalent to stating that the kinetics are always in the
    low-pressure limit.
