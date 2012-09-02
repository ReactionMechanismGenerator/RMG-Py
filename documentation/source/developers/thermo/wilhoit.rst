******************
Wilhoit polynomial
******************

The Wilhoit polynomial is an expression for heat capacity that is guaranteed to
give the correct limits at zero and infinite temperature, and gives a very
reasonable shape to the heat capacity profile in between:

.. math::

    C_\mathrm{p}(T) = C_\mathrm{p}(0) + \left[ C_\mathrm{p}(\infty) -
    C_\mathrm{p}(0) \right] y^2 \left[ 1 + (y - 1) \sum_{i=0}^3 a_i y^i \right]

Above, :math:`y \equiv T/(T + B)` is a scaled temperature that ranges
from zero to one based on the value of the coefficient :math:`B`, and
:math:`a_0`, :math:`a_1`, :math:`a_2`, and :math:`a_3` are the Wilhoit
polynomial coefficients.

The enthalpy is given by

.. math::

    H(T) & = H_0 +
    C_\mathrm{p}(0) T + \left[ C_\mathrm{p}(\infty) - C_\mathrm{p}(0) \right] T \\
    & \left\{ \left[ 2 + \sum_{i=0}^3 a_i \right]
    \left[ \frac{1}{2}y - 1 + \left( \frac{1}{y} - 1 \right) \ln \frac{T}{y} \right]
    + y^2 \sum_{i=0}^3 \frac{y^i}{(i+2)(i+3)} \sum_{j=0}^3 f_{ij} a_j
    \right\}

where :math:`f_{ij} = 3 + j` if :math:`i = j`, :math:`f_{ij} = 1` if
:math:`i > j`, and :math:`f_{ij} = 0` if :math:`i < j`.

The entropy is given by

.. math::

    S(T) = S_0 +
    C_\mathrm{p}(\infty) \ln T - \left[ C_\mathrm{p}(\infty) - C_\mathrm{p}(0) \right]
    \left[ \ln y + \left( 1 + y \sum_{i=0}^3 \frac{a_i y^i}{2+i} \right) y
    \right]


The low-temperature limit :math:`C_\mathrm{p}(0)` is :math:`3.5R` for linear
molecules and :math:`4R` for nonlinear molecules. The high-temperature limit 
:math:`C_\mathrm{p}(\infty)` is taken to be 
:math:`\left[ 3 N_\mathrm{atoms} - 1.5 \right] R` for linear molecules and
:math:`\left[ 3 N_\mathrm{atoms} - (2 + 0.5 N_\mathrm{rotors}) \right] R`
for nonlinear molecules, for a molecule composed of :math:`N_\mathrm{atoms}`
atoms and :math:`N_\mathrm{rotors}` internal rotors.

.. autoclass:: rmgpy.thermo.Wilhoit
    :members:
