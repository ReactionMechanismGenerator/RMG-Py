*********************
Chebyshev polynomials
*********************

The Chebyshev polynomial formulation is a means of fitting a wide range of
complex :math:`k(T,P)` behavior. However, there is no meaningful physical
interpretation of the polynomial-based fit, and one must take care to minimize
the magnitude of Runge's phenomenon. The formulation is as follows:

.. math:: \log k(T,P) = \sum_{t=1}^{N_T} \sum_{p=1}^{N_P} \alpha_{tp} \phi_t(\tilde{T}) \phi_p(\tilde{P})
    
Above, :math:`\alpha_{tp}` is a constant, :math:`\phi_n(x)` is the
Chebyshev polynomial of degree :math:`n` evaluated at :math:`x`, and

.. math:: \tilde{T} \equiv \frac{2T^{-1} - T_\mathrm{min}^{-1} - T_\mathrm{max}^{-1}}{T_\mathrm{max}^{-1} - T_\mathrm{min}^{-1}}

.. math:: \tilde{P} \equiv \frac{2 \log P - \log P_\mathrm{min} - \log P_\mathrm{max}}{\log P_\mathrm{max} - \log P_\mathrm{min}}

are reduced temperature and reduced pressure designed to map the ranges
:math:`(T_\mathrm{min}, T_\mathrm{max})` and
:math:`(P_\mathrm{min}, P_\mathrm{max})` to :math:`(-1, 1)`.

.. autoclass:: rmgpy.kinetics.Chebyshev
    :members:
