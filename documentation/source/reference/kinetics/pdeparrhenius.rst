****************************
rmgpy.kinetics.PDepArrhenius
****************************

.. autoclass:: rmgpy.kinetics.PDepArrhenius

    The pressure-dependent Arrhenius formulation is sometimes used to extend the
    Arrhenius expression to handle pressure-dependent kinetics. The formulation
    simply parameterizes :math:`A`, :math:`n`, and :math:`E_\mathrm{a}` to be
    dependent on pressure:

    .. math:: k(T,P) = A(P) \left( \frac{T}{T_0} \right)^{n(P)} \exp \left( -\frac{E_\mathrm{a}(P)}{RT} \right)

    Although this suggests some physical insight, the :math:`k(T,P)` data is often
    highly complex and non-Arrhenius, limiting the usefulness of this formulation
    to simple systems.
