**********
Unit types
**********

Units can be classified into categories based on the associated dimensionality.
For example, miles and kilometers are both units of length; seconds and hours
are both units of time, etc. Clearly, quantities of different unit types are
fundamentally different.

RMG provides functions that create physical quantities (scalar or array) and
validate the units for a variety of unit types. This prevents the user from
inadvertently mixing up their units - e.g. by setting an enthalpy with entropy
units - which should reduce errors. RMG recognizes the following unit types:

    ======================== ====================== ========================
    Function                 Unit type              SI unit
    ======================== ====================== ========================
    :func:`Acceleration`     acceleration           :math:`\mathrm{m/s^2}`
    :func:`Area`             area                   :math:`\mathrm{m^2}`
    :func:`Concentration`    concentration          :math:`\mathrm{mol/cm^3}`
    :func:`Dimensionless`    dimensionless          
    :func:`Energy`           energy                 :math:`\mathrm{J/mol}`
    :func:`Entropy`          entropy                :math:`\mathrm{J/mol \cdot K}`
    :func:`Flux`             flux                   :math:`\mathrm{mol/cm^2 \cdot s}`
    :func:`Frequency`        frequency              :math:`\mathrm{cm^{-1}}`
    :func:`Force`            force                  :math:`\mathrm{N}`
    :func:`Inertia`          inertia                :math:`\mathrm{kg \cdot m^2}`
    :func:`Length`           length                 :math:`\mathrm{m}`
    :func:`Mass`             mass                   :math:`\mathrm{kg}`
    :func:`Momentum`         momentum               :math:`\mathrm{kg \cdot m/s^2}`
    :func:`Power`            power                  :math:`\mathrm{W}`
    :func:`Pressure`         pressure               :math:`\mathrm{Pa}`
    :func:`RateCoefficient`  rate coefficient       :math:`\mathrm{s^{-1}}`, :math:`\mathrm{m^3/mol \cdot s}`, :math:`\mathrm{m^6/mol^2 \cdot s}`, :math:`\mathrm{m^9/mol^3 \cdot s}`
    :func:`Temperature`      temperature            :math:`\mathrm{K}`
    :func:`Time`             time                   :math:`\mathrm{s}`
    :func:`Velocity`         velocity               :math:`\mathrm{m/s}`
    :func:`Volume`           volume                 :math:`\mathrm{m^3}`
    ======================== ====================== ========================

In RMG, all energies, heat capacities, concentrations, fluxes, and rate
coefficients are treated as intensive; this means that these quantities are
always expressed "per mole" or "per molecule". All other unit types are
extensive. A special exception is added for mass so as to allow for coercion
of g/mol to amu.

RMG also handles rate coefficient units as a special case, as there are multiple
allowed dimensionalities based on the reaction order. Note that RMG generally
does not attempt to verify that the rate coefficient units match the reaction
order, but only that it matches one of the possibilities.

The table above gives the SI unit that RMG uses internally to work with
physical quantities. This does not necessarily correspond with the units used
when outputting values. For example, pressures are often output in units of
:math:`\mathrm{bar}` instead of :math:`\mathrm{Pa}`, and moments of inertia in 
:math:`\mathrm{amu*angstrom^2}` instead of :math:`\mathrm{kg*m^2}`. The
recommended rule of thumb is to use prefixed SI units (or aliases thereof) in 
the output; for example, use :math:`\mathrm{kJ/mol}` instead of 
:math:`\mathrm{kcal/mol}` for energy values.
