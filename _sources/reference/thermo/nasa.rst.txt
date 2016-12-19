*****************
rmgpy.thermo.NASA
*****************
.. _thermoNASA:

.. currentmodule:: rmgpy.thermo

.. autoclass:: rmgpy.thermo.NASA

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

Here is an example of a NASA entry::

	entry(
    index = 2,
    label = "octane",
    molecule = 
	"""
	1 C 0 {2,S}
	2 C 0 {1,S} {3,S}
	3 C 0 {2,S} {4,S}
	4 C 0 {3,S} {5,S}
	5 C 0 {4,S} {6,S}
	6 C 0 {5,S} {7,S}
	7 C 0 {6,S} {8,S}
	8 C 0 {7,S}
	""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.25245480E+01,-1.01018826E-02,2.21992610E-04,-2.84863722E-07,1.12410138E-10,-2.98434398E+04,-1.97109989E+01], Tmin=(200,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.09430708E+01,4.41691018E-02,-1.53261633E-05,2.30544803E-09,-1.29765727E-13,-3.55755088E+04,-8.10637726E+01], Tmin=(1000,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (200,'K'),
        Tmax = (6000,'K'),
    ),
    reference = Reference(authors=["check on burcat"], title='burcat', year="1999", url="http://www.me.berkeley.edu/gri-mech/version30/text30.html"),
    referenceType = "review",
    shortDesc = u"""""",
    longDesc = 
	u"""
	
	""",
	)
