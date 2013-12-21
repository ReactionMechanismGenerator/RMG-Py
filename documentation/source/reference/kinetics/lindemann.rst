************************
rmgpy.kinetics.Lindemann
************************

.. autoclass:: rmgpy.kinetics.Lindemann

    The Lindemann model qualitatively predicts the falloff of some simple
    pressure-dependent reaction kinetics. The formulation is as follows:

    .. math:: k(T,P) = k_\infty(T) \left[ \frac{P_\mathrm{r}}{1 + P_\mathrm{r}} \right]

    where

    .. math::

        P_\mathrm{r} &= \frac{k_0(T)}{k_\infty(T)} [\mathrm{M}]

        k_0(T) &= A_0 T^{n_0} \exp \left( - \frac{E_0}{RT} \right)

        k_\infty(T) &= A_\infty T^{n_\infty} \exp \left( - \frac{E_\infty}{RT} \right)

    and :math:`[\mathrm{M}] \approx P/RT` is the concentration of the bath gas. The
    Arrhenius expressions :math:`k_0(T)` and :math:`k_\infty(T)` represent the
    low-pressure and high-pressure limit kinetics, respectively.
