******************
Arrhenius kinetics
******************

The Arrhenius equation, given below, accurately reproduces the kinetics of many
reaction families:

.. math:: k(T) = A \left( \frac{T}{T_0} \right)^n \exp \left( -\frac{E_\mathrm{a}}{RT} \right)

Above, :math:`A` is the preexponential factor, :math:`T_0` is the reference
temperature, :math:`n` is the temperature exponent, and :math:`E_\mathrm{a}`
is the activation energy.

.. autoclass:: rmgpy.kinetics.Arrhenius
    :members:

.. autoclass:: rmgpy.kinetics.MultiArrhenius
    :members:
    
.. autoclass:: rmgpy.kinetics.ArrheniusEP
    :members:

