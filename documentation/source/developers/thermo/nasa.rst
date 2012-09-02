***************
NASA polynomial
***************

The NASA polynomial is another representation of the heat capacity, enthalpy,
and entropy using seven or nine coefficients
:math:`\mathbf{a} = \left[a_{-2}\ a_{-1}\ a_0\ a_1\ a_2\ a_3\ a_4\ a_5\ a_6 \right]`.
The relevant thermodynamic parameters are evaluated via the expressions
    
.. math:: \frac{C_\mathrm{p}(T)}{R} = a_{-2} T^{-2} + a_{-1} T^{-1} + a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4

.. math:: \frac{H(T)}{RT} = - a_{-2} T^{-2} + a_{-1} T^{-1} \ln T + a_0 + \frac{1}{2} a_1 T + \frac{1}{3} a_2 T^2 + \frac{1}{4} a_3 T^3 + \frac{1}{5} a_4 T^4 + \frac{a_5}{T}

.. math:: \frac{S(T)}{R} = -\frac{1}{2} a_{-2} T^{-2} - a_{-1} T^{-1} + a_0 \ln T + a_1 T + \frac{1}{2} a_2 T^2 + \frac{1}{3} a_3 T^3 + \frac{1}{4} a_4 T^4 + a_6

In the seven-coefficient version, :math:`a_{-2} = a_{-1} = 0`.

As simple polynomial expressions, the NASA polynomial is faster to evaluate
when compared to the Wilhoit model; however, it does not have the nice
physical behavior of the Wilhoit representation. Often multiple NASA polynomials
are used to accurately represent the thermodynamics of a system over a wide
temperature range.

.. autoclass:: rmgpy.thermo.NASA
    :members:
