*******************
rmgpy.kinetics.Troe
*******************

.. autoclass:: rmgpy.kinetics.Troe

    The Troe model attempts to make the Lindemann model quantitative by introducing
    a broadening factor :math:`F`. The formulation is as follows:

    .. math:: k(T,P) = k_\infty(T) \left[ \frac{P_\mathrm{r}}{1 + P_\mathrm{r}} \right] F

    where

    .. math::

        P_\mathrm{r} &= \frac{k_0(T)}{k_\infty(T)} [\mathrm{M}]

        k_0(T) &= A_0 T^{n_0} \exp \left( - \frac{E_0}{RT} \right)

        k_\infty(T) &= A_\infty T^{n_\infty} \exp \left( - \frac{E_\infty}{RT} \right)

    and :math:`[\mathrm{M}] \approx P/RT` is the concentration of the bath gas. The
    Arrhenius expressions :math:`k_0(T)` and :math:`k_\infty(T)` represent the
    low-pressure and high-pressure limit kinetics, respectively. The broadening
    factor :math:`F` is computed via

    .. math::

        \log F &= \left\{1 + \left[ \frac{\log P_\mathrm{r} + c}{n - d (\log P_\mathrm{r} + c)} \right]^2 \right\}^{-1} \log F_\mathrm{cent}

        c &= -0.4 - 0.67 \log F_\mathrm{cent}

        n &= 0.75 - 1.27 \log F_\mathrm{cent}

        d &= 0.14

        F_\mathrm{cent} &= (1 - \alpha) \exp \left( -T/T_3 \right) + \alpha \exp \left( -T/T_1 \right) + \exp \left( -T_2/T \right)

