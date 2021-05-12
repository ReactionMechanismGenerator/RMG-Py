# cython: embedsignature=True, cdivision=True

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import numpy as np
cimport numpy as np
from libc.math cimport exp, sqrt

cimport rmgpy.constants as constants
import rmgpy.quantity as quantity
from rmgpy.exceptions import KineticsError

# Prior to numpy 1.14, `numpy.linalg.lstsq` does not accept None as a value
RCOND = -1 if int(np.__version__.split('.')[1]) < 14 else None

################################################################################

cdef class StickingCoefficient(KineticsModel):
    """
    A kinetics model to give Sticking Coefficients for surface adsorption,
    following Arrhenius form. 
    Similar to :class:`Arrhenius` but with different units for `A`.
    The attributes are:

    ======================= =============================================================
    Attribute               Description
    ======================= =============================================================
    `A`                     The preexponential factor
    `T0`                    The reference temperature
    `n`                     The temperature exponent
    `Ea`                    The activation energy
    `Tmin`                  The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`                  The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`                  The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`                  The maximum pressure at which the model is valid, or zero if unknown or undefined
    `coverage_dependence`   A dictionary of coverage dependent parameters to a certain surface species with:
                             `a`, the coefficient for exponential dependence on the coverage,
                             `m`, the power-law exponent of coverage dependence, and
                             `E`, the activation energy dependence on coverage.
    `comment`               Information about the model (e.g. its source)
    ======================= =============================================================
    
    """

    def __init__(self, A=None, n=0.0, Ea=None, T0=(1.0, "K"), Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 coverage_dependence=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.A = A
        self.n = n
        self.Ea = Ea
        self.T0 = T0
        self.coverage_dependence = coverage_dependence

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        StickingCoefficient object.
        """
        string = 'StickingCoefficient(A={0!r}, n={1!r}, Ea={2!r}, T0={3!r}'.format(self.A, self.n, self.Ea, self.T0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.coverage_dependence:
            string += ", coverage_dependence={"
            for species, parameters in self.coverage_dependence.items():
                string += f"{species.to_chemkin()!r}: {{'a':{parameters['a']}, 'm':{parameters['m']}, 'E':({parameters['E'].value}, '{parameters['E'].units}')}},"
            string += "}"
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a StickingCoefficient object.
        """
        return (StickingCoefficient, (self.A, self.n, self.Ea, self.T0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                                      self.coverage_dependence, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.Dimensionless(value)

    property n:
        """The temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property Ea:
        """The activation energy."""
        def __get__(self):
            return self._Ea
        def __set__(self, value):
            self._Ea = quantity.Energy(value)

    property T0:
        """The reference temperature."""
        def __get__(self):
            return self._T0
        def __set__(self, value):
            self._T0 = quantity.Temperature(value)

    property coverage_dependence:
        """The coverage dependence parameters."""
        def __get__(self):
            return self._coverage_dependence
        def __set__(self, value):
             self._coverage_dependence = {}
             if value:
                 for species, parameters in value.items():
                     processed_parameters = {'E': quantity.Energy(parameters['E']),
                                             'm': quantity.Dimensionless(parameters['m']),
                                             'a': quantity.Dimensionless(parameters['a'])}
                     self._coverage_dependence[species] = processed_parameters

    cpdef double get_sticking_coefficient(self, double T) except -1:
        """
        Return the sticking coefficient (dimensionless) at temperature `T` in K. 
        """
        cdef double A, n, Ea, T0, stickingCoefficient
        A = self._A.value_si
        n = self._n.value_si
        Ea = self._Ea.value_si
        T0 = self._T0.value_si

        stickingCoefficient = A * (T / T0) ** n * exp(-Ea / (constants.R * T))
        if stickingCoefficient < 0:
            raise ValueError("Sticking coefficients cannot be negative, check your preexponential factor.")
        return min(stickingCoefficient, 1.0)

    cpdef fit_to_data(self, np.ndarray Tlist, np.ndarray klist, str kunits, double T0=1,
                    np.ndarray weights=None, bint three_params=True):
        """
        Fit Arrhenius parameters to a set of sticking coefficient data `klist`
        in units of `kunits` corresponding to a set of temperatures `Tlist` in
        K. A linear least-squares fit is used, which guarantees that the
        resulting parameters provide the best possible approximation to the
        data.
        """
        import scipy.stats

        assert len(Tlist) == len(klist), "length of temperatures and rates must be the same"
        if len(Tlist) < 3 + three_params:
            raise KineticsError('Not enough degrees of freedom to fit this Arrhenius expression')
        if three_params:
            A = np.zeros((len(Tlist), 3), np.float64)
            A[:, 0] = np.ones_like(Tlist)
            A[:, 1] = np.log(Tlist / T0)
            A[:, 2] = -1.0 / constants.R / Tlist
        else:
            A = np.zeros((len(Tlist), 2), np.float64)
            A[:, 0] = np.ones_like(Tlist)
            A[:, 1] = -1.0 / constants.R / Tlist
        b = np.log(klist)
        if weights is not None:
            for n in range(b.size):
                A[n, :] *= weights[n]
                b[n] *= weights[n]
        x, residues, rank, s = np.linalg.lstsq(A, b, rcond=RCOND)

        # Determine covarianace matrix to obtain parameter uncertainties
        count = klist.size
        cov = residues[0] / (count - 3) * np.linalg.inv(np.dot(A.T, A))
        t = scipy.stats.t.ppf(0.975, count - 3)

        if not three_params:
            x = np.array([x[0], 0, x[1]])
            cov = np.array([[cov[0, 0], 0, cov[0, 1]], [0, 0, 0], [cov[1, 0], 0, cov[1, 1]]])

        self.A = (exp(x[0]), kunits)
        self.n = x[1]
        self.Ea = (x[2] * 0.001, "kJ/mol")
        self.T0 = (T0, "K")
        self.Tmin = (np.min(Tlist), "K")
        self.Tmax = (np.max(Tlist), "K")
        self.comment = 'Fitted to {0:d} data points; dA = *|/ {1:g}, dn = +|- {2:g}, dEa = +|- {3:g} kJ/mol'.format(
            len(Tlist),
            exp(sqrt(cov[0, 0])),
            sqrt(cov[1, 1]),
            sqrt(cov[2, 2]) * 0.001,
        )

        return self

    cpdef change_t0(self, double T0):
        """
        Changes the reference temperature used in the exponent to `T0` in K, 
        and adjusts the preexponential factor accordingly.
        """
        self._A.value_si /= (self._T0.value_si / T0) ** self._n.value_si
        self._T0.value_si = T0

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's Arrhenius.) Otherwise returns ``False``.
        """
        if not isinstance(other_kinetics, StickingCoefficient):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if (not self.A.equals(other_kinetics.A) or not self.n.equals(other_kinetics.n)
                or not self.Ea.equals(other_kinetics.Ea) or not self.T0.equals(other_kinetics.T0)):
            return False

        return True

    cpdef change_rate(self, double factor):
        """
        Changes A factor in Arrhenius expression by multiplying it by a ``factor``.
        """
        self._A.value_si *= factor

    cpdef to_html(self):
        """
        Return an HTML rendering.
        """
        cdef double T
        cdef str string
        cdef list Tdata

        Tdata = [500, 1000, 1500, 2000]

        string = '<table class="KineticsData">\n<tr class="KineticsData_Tdata"><th>T/[K]</th>\n'
        try:
            for T in Tdata:
                string += '<td>{0:.0f}</td>'.format(T)

            string += '\n</tr><tr class="KineticsData_kdata"><th>Sticking Coefficient\n    '

            for T in Tdata:
                string += '<td>{0:+.2f}</td>'.format(self.get_sticking_coefficient(T))
        except:
            string += '<td>An error occurred in processing kinetics</td>'
        string += '\n</tr></table>'

        string += "<span class='KineticsData_repr'>{0!r}</span>".format(self)

        return string

################################################################################
cdef class StickingCoefficientBEP(KineticsModel):
    """
    A kinetics model based on the Arrhenius expression, to give 
    Sticking Coefficient for surface adsorption, using the
    Bronsted-Evans-Polanyi equation to determine the activation energy. 
    Similar to :class:`ArrheniusEP`, but with different units for `A`.
    Sticking Coefficients are between 0 and 1.
    The attributes are:

    ======================= =============================================================
    Attribute               Description
    ======================= =============================================================
    `A`                     The preexponential factor
    `n`                     The temperature exponent
    `alpha`                 The Evans-Polanyi slope
    `E0`                    The activation energy for a thermoneutral reaction
    `Tmin`                  The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`                  The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`                  The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`                  The maximum pressure at which the model is valid, or zero if unknown or undefined
    `coverage_dependence`   A dictionary of coverage dependent parameters to a certain surface species with:
                             `a`, the coefficient for exponential dependence on the coverage,
                             `m`, the power-law exponent of coverage dependence, and
                             `E`, the activation energy dependence on coverage.
    `comment`               Information about the model (e.g. its source)
    ======================= =============================================================
    
    """

    def __init__(self, A=None, n=0.0, alpha=0.0, E0=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 coverage_dependence=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.A = A
        self.n = n
        self.alpha = alpha
        self.E0 = E0
        self.coverage_dependence = coverage_dependence

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        StickingCoefficientBEP object.
        """
        string = 'StickingCoefficientBEP(A={0!r}, n={1!r}, alpha={2!r}, E0={3!r}'.format(self.A, self.n, self.alpha,
                                                                                         self.E0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.coverage_dependence:
            string += ", coverage_dependence={"
            for species, parameters in self.coverage_dependence.items():
                string += f"{species.to_chemkin()!r}: {{'a':{parameters['a']}, 'm':{parameters['m']}, 'E':({parameters['E'].value}, '{parameters['E'].units}')}},"
            string += "}"
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an StickingCoefficientBEP object.
        """
        return (StickingCoefficientBEP, (self.A, self.n, self.alpha, self.E0, self.Tmin, self.Tmax,
                                         self.Pmin, self.Pmax, self.coverage_dependence, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.Dimensionless(value)

    property n:
        """The temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property alpha:
        """The Bronsted-Evans-Polanyi slope."""
        def __get__(self):
            return self._alpha
        def __set__(self, value):
            self._alpha = quantity.Dimensionless(value)

    property E0:
        """The activation energy for a thermoneutral reaction."""
        def __get__(self):
            return self._E0
        def __set__(self, value):
            self._E0 = quantity.Energy(value)

    property coverage_dependence:
        """The coverage dependence parameters."""
        def __get__(self):
            return self._coverage_dependence
        def __set__(self, value):
             self._coverage_dependence = {}
             if value:
                 for species, parameters in value.items():
                     processed_parameters = {'E': quantity.Energy(parameters['E']),
                                             'm': quantity.Dimensionless(parameters['m']),
                                             'a': quantity.Dimensionless(parameters['a'])}
                     self._coverage_dependence[species] = processed_parameters

    cpdef double get_sticking_coefficient(self, double T, double dHrxn=0.0) except -1:
        """
        Return the sticking coefficient (dimensionless) at
        temperature `T` in K and enthalpy of reaction `dHrxn` in J/mol. 
        """
        cdef double A, n, Ea, stickingCoefficient
        Ea = self.get_activation_energy(dHrxn)
        A = self._A.value_si
        n = self._n.value_si

        stickingCoefficient = A * T ** n * exp(-Ea / (constants.R * T))
        assert 0 <= stickingCoefficient
        return min(stickingCoefficient, 1.0)

    cpdef double get_activation_energy(self, double dHrxn) except -1:
        """
        Return the activation energy in J/mol corresponding to the given
        enthalpy of reaction `dHrxn` in J/mol.
        """
        cdef double Ea
        Ea = self._alpha.value_si * dHrxn + self._E0.value_si
        if self._E0.value_si > 0:
            if dHrxn < 0.0 and Ea < 0.0:
                Ea = 0.0
            elif dHrxn > 0.0 and Ea < dHrxn:
                Ea = dHrxn
        return Ea

    cpdef StickingCoefficient to_arrhenius(self, double dHrxn):
        """
        Return an :class:`StickingCoefficient` instance of the kinetics model using the
        given enthalpy of reaction `dHrxn` to determine the activation energy.
        
        Note that despite its name it does not return a :class:`Arrhenius` object.
        """
        return StickingCoefficient(
            A=self.A,
            n=self.n,
            Ea=(self.get_activation_energy(dHrxn) * 0.001, "kJ/mol"),
            T0=(1, "K"),
            Tmin=self.Tmin,
            Tmax=self.Tmax,
            coverage_dependence=self.coverage_dependence,
            comment=self.comment,
        )

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match type, temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's StickingCoefficient.) Otherwise returns ``False``.
        """
        if not isinstance(other_kinetics, StickingCoefficientBEP):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if (not self.A.equals(other_kinetics.A) or not self.n.equals(other_kinetics.n)
                or not self.alpha.equals(other_kinetics.alpha) or not self.E0.equals(other_kinetics.E0)):
            return False

        return True

    cpdef change_rate(self, double factor):
        """
        Changes A factor by multiplying it by a ``factor``.
        """
        self._A.value_si *= factor

    def set_cantera_kinetics(self, ct_reaction, species_list=[]):
        """
        Sets a cantera ElementaryReaction() object in an Arrhenius form.
        """
        raise NotImplementedError('set_cantera_kinetics() is not implemented for StickingCoefficientBEP class kinetics.')

################################################################################

cdef class SurfaceArrhenius(Arrhenius):
    """
    A kinetics model based on (modified) Arrhenius for surface reactions.
    
    It is very similar to the gas phase :class:`Arrhenius`
    
    The attributes are:

    ======================= =============================================================
    Attribute               Description
    ======================= =============================================================
    `A`                     The preexponential factor
    `T0`                    The reference temperature
    `n`                     The temperature exponent
    `Ea`                    The activation energy
    `Tmin`                  The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`                  The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`                  The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`                  The maximum pressure at which the model is valid, or zero if unknown or undefined
    `coverage_dependence`   A dictionary of coverage dependent parameters to a certain surface species with:
                             `a`, the coefficient for exponential dependence on the coverage,
                             `m`, the power-law exponent of coverage dependence, and
                             `E`, the activation energy dependence on coverage.
    `uncertainty`           Uncertainty information
    `comment`               Information about the model (e.g. its source)
    ======================= =============================================================
    """
    def __init__(self, A=None, n=0.0, Ea=None, T0=(1.0, "K"), Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 coverage_dependence=None, uncertainty=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, uncertainty=uncertainty,
                               comment=comment)
        self.A = A
        self.n = n
        self.Ea = Ea
        self.T0 = T0
        self.coverage_dependence = coverage_dependence

    property A:
        """The preexponential factor. 
    
        This is the only thing different from a normal Arrhenius class."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.SurfaceRateCoefficient(value)

    property coverage_dependence:
        """The coverage dependence parameters."""
        def __get__(self):
            return self._coverage_dependence
        def __set__(self, value):
             self._coverage_dependence = {}
             if value:
                 for species, parameters in value.items():
                     processed_parameters = {'E': quantity.Energy(parameters['E']),
                                             'm': quantity.Dimensionless(parameters['m']),
                                             'a': quantity.Dimensionless(parameters['a'])}
                     self._coverage_dependence[species] = processed_parameters

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        SurfaceArrhenius object.
        """
        string = 'SurfaceArrhenius(A={0!r}, n={1!r}, Ea={2!r}, T0={3!r}'.format(self.A, self.n, self.Ea, self.T0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.coverage_dependence:
            string += ", coverage_dependence={"
            for species, parameters in self.coverage_dependence.items():
                string += f"{species.to_chemkin()!r}: {{'a':{parameters['a']}, 'm':{parameters['m']}, 'E':({parameters['E'].value}, '{parameters['E'].units}')}},"
            string += "}"
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.uncertainty is not None: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a SurfaceArrhenius object.
        """
        return (SurfaceArrhenius, (self.A, self.n, self.Ea, self.T0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                                   self.coverage_dependence, self.uncertainty, self.comment))

    cpdef SurfaceChargeTransfer to_surface_charge_transfer(self, double V0, double electrons=-1):
        """
        Return an :class:`SurfaceChargeTransfer` instance of the kinetics model with reversible
        potential `V0` in Volts and electron stochiometric coeff ` electrons`
        """
        return SurfaceChargeTransfer(
            A=self.A,
            n=self.n,
            electrons= electrons,
            Ea=self.Ea,
            V0=(V0,'V'),
            T0=(1, "K"),
            Tmin=self.Tmin,
            Tmax=self.Tmax,
            uncertainty = self.uncertainty,
            comment=self.comment,
        )


################################################################################

cdef class SurfaceArrheniusBEP(ArrheniusEP):
    """
    A kinetics model based on the (modified) Arrhenius equation, using the
    Bronsted-Evans-Polanyi equation to determine the activation energy. 
    
    It is very similar to the gas-phase :class:`ArrheniusEP`.
    The only differences being the A factor has different units,
    (and the catalysis community prefers to call it BEP rather than EP!)
    and has a coverage_dependence parameter for coverage dependence
    
    The attributes are:

    ======================= =============================================================
    Attribute               Description
    ======================= =============================================================
    `A`                     The preexponential factor
    `n`                     The temperature exponent
    `alpha`                 The Evans-Polanyi slope
    `E0`                    The activation energy for a thermoneutral reaction
    `Tmin`                  The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`                  The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`                  The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`                  The maximum pressure at which the model is valid, or zero if unknown or undefined
    `coverage_dependence`   A dictionary of coverage dependent parameters to a certain surface species with:
                             `a`, the coefficient for exponential dependence on the coverage,
                             `m`, the power-law exponent of coverage dependence, and
                             `E`, the activation energy dependence on coverage.
    `uncertainty`           Uncertainty information
    `comment`               Information about the model (e.g. its source)
    ======================= =============================================================
    
    """
    def __init__(self, A=None, n=0.0, alpha=0.0, E0=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, uncertainty=None,
                 coverage_dependence=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, uncertainty=uncertainty,
                               comment=comment,)
        self.A = A
        self.n = n
        self.alpha = alpha
        self.E0 = E0
        self.coverage_dependence = coverage_dependence

    property A:
        """The preexponential factor. 
    
        This is the only thing different from a normal ArrheniusEP class."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.SurfaceRateCoefficient(value)

    property coverage_dependence:
        """The coverage dependence parameters."""
        def __get__(self):
            return self._coverage_dependence
        def __set__(self, value):
             self._coverage_dependence = {}
             if value:
                 for species, parameters in value.items():
                     processed_parameters = {'E': quantity.Energy(parameters['E']),
                                             'm': quantity.Dimensionless(parameters['m']),
                                             'a': quantity.Dimensionless(parameters['a'])}
                     self._coverage_dependence[species] = processed_parameters

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        SurfaceArrheniusBEP object.
        """
        string = 'SurfaceArrheniusBEP(A={0!r}, n={1!r}, alpha={2!r}, E0={3!r}'.format(self.A, self.n, self.alpha,
                                                                                      self.E0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.coverage_dependence:
            string += ", coverage_dependence={"
            for species, parameters in self.coverage_dependence.items():
                string += f"{species.to_chemkin()!r}: {{'a':{parameters['a']}, 'm':{parameters['m']}, 'E':({parameters['E'].value}, '{parameters['E'].units}')}},"
            string += "}"
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.uncertainty is not None: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an SurfaceArrheniusBEP object.
        """
        return (SurfaceArrheniusBEP, (self.A, self.n, self.alpha, self.E0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                                      self.uncertainty, self.coverage_dependence, self.comment))

    cpdef SurfaceArrhenius to_arrhenius(self, double dHrxn):
        """
        Return an :class:`SurfaceArrhenius` instance of the kinetics model using the
        given enthalpy of reaction `dHrxn` to determine the activation energy.
        
        Note that despite its name it does not return a :class:`Arrhenius` object
        (although :class:`SurfaceArrhenius` is a subclass of :class:`Arrhenius` 
        so in a way, it does).
        """
        return SurfaceArrhenius(
            A=self.A,
            n=self.n,
            Ea=(self.get_activation_energy(dHrxn) * 0.001, "kJ/mol"),
            T0=(1, "K"),
            Tmin=self.Tmin,
            Tmax=self.Tmax,
            uncertainty=self.uncertainty,
            coverage_dependence=self.coverage_dependence,
            comment=self.comment,
        )

################################################################################

cdef class SurfaceChargeTransfer(KineticsModel):

    """
    A kinetics model for surface charge transfer reactions

    It is very similar to the :class:`SurfaceArrhenius`, but the Ea is potential-dependent


    The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `T0`            The reference temperature
    `n`             The temperature exponent
    `Ea`            The activation energy
    `electrons`     The stochiometry coeff for electrons (negative if reactant, positive if product)
    `V0`            The reference potential
    `alpha`         The charge transfer coefficient
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, A=None, n=0.0, Ea=None, V0=None, alpha=0.5, electrons=-1, T0=(1.0, "K"), Tmin=None, Tmax=None, 
                Pmin=None, Pmax=None, uncertainty=None, comment=''):

        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, uncertainty=uncertainty,
                comment=comment)

        self.alpha = alpha
        self.A = A
        self.n = n
        self.Ea = Ea
        self.T0 = T0
        self.electrons = electrons
        self.V0 = V0

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Arrhenius object.
        """
        string = 'SurfaceChargeTransfer(A={0!r}, n={1!r}, Ea={2!r}, V0={3!r}, alpha={4!r}, electrons={5!r}, T0={6!r}'.format(
            self.A, self.n, self.Ea, self.V0, self.alpha, self.electrons, self.T0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.uncertainty: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a SurfaceChargeTransfer object.
        """
        return (SurfaceChargeTransfer, (self.A, self.n, self.Ea, self.V0, self.alpha, self.electrons, self.T0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                            self.uncertainty, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.SurfaceRateCoefficient(value)

    property n:
        """The temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property Ea:
        """The activation energy."""
        def __get__(self):
            return self._Ea
        def __set__(self, value):
            self._Ea = quantity.Energy(value)

    property T0:
        """The reference temperature."""
        def __get__(self):
            return self._T0
        def __set__(self, value):
            self._T0 = quantity.Temperature(value)

    property V0:
        """The reference potential."""
        def __get__(self):
            return self._V0
        def __set__(self, value):
            self._V0 = quantity.Potential(value)

    property electrons:
        """The number of electrons transferred."""
        def __get__(self):
            return self._electrons
        def __set__(self, value):
            self._electrons = quantity.Dimensionless(value)

    property alpha:
        """The charge transfer coefficient."""
        def __get__(self):
            return self._alpha
        def __set__(self, value):
            self._alpha = quantity.Dimensionless(value)

    cpdef double get_activation_energy_from_potential(self, double V=0.0, bint non_negative=True):
        """
        Return the effective activation energy (in J/mol) at specificed potential (in Volts).
        """
        cdef double  electrons, alpha, Ea, V0
        
        electrons = self._electrons.value_si
        alpha = self._alpha.value_si
        Ea = self._Ea.value_si
        V0 = self._V0.value_si

        Ea -= alpha * electrons * constants.F * (V-V0)

        if non_negative is True:
            if Ea < 0:
                Ea = 0.0

        return Ea

    cpdef double get_rate_coefficient(self, double T, double V=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^2,
        mol, and s at temperature `T` in K.
        """
        cdef double A, n, V0, T0, Ea

        A = self._A.value_si
        n = self._n.value_si
        V0 = self._V0.value_si
        T0 = self._T0.value_si

        if V != V0:
            Ea = self.get_activation_energy_from_potential(V)
        else:
            Ea = self._Ea.value_si

        return A * (T / T0) ** n * exp(-Ea / (constants.R * T)) 

    cpdef change_t0(self, double T0):
        """
        Changes the reference temperature used in the exponent to `T0` in K,
        and adjusts the preexponential factor accordingly.
        """
        self._A.value_si /= (self._T0.value_si / T0) ** self._n.value_si
        self._T0.value_si = T0

    cpdef change_v0(self, double V0):
        """
        Changes the reference potential to `V0` in volts, and adjusts the
        activation energy `Ea` accordingly.
        """

        self._Ea.value_si = self.get_activation_energy_from_potential(V0)
        self._V0.value_si = V0

    cpdef fit_to_data(self, np.ndarray Tlist, np.ndarray klist, str kunits, double T0=1,
                      np.ndarray weights=None, bint three_params=False):
        """
        Fit the Arrhenius parameters to a set of rate coefficient data `klist`
        in units of `kunits` corresponding to a set of temperatures `Tlist` in
        K. A linear least-squares fit is used, which guarantees that the
        resulting parameters provide the best possible approximation to the
        data.
        """
        import scipy.stats
        if not all(np.isfinite(klist)):
            raise  ValueError("Rates must all be finite, not inf or NaN")
        if any(klist<0):
            if not all(klist<0):
                raise ValueError("Rates must all be positive or all be negative.")
            rate_sign_multiplier = -1
            klist = -1 * klist
        else:
            rate_sign_multiplier = 1

        assert len(Tlist) == len(klist), "length of temperatures and rates must be the same"
        if len(Tlist) < 3 + three_params:
            raise KineticsError('Not enough degrees of freedom to fit this Arrhenius expression')
        if three_params:
            A = np.zeros((len(Tlist), 3), np.float64)
            A[:, 0] = np.ones_like(Tlist)
            A[:, 1] = np.log(Tlist / T0)
            A[:, 2] = -1.0 / constants.R / Tlist
        else:
            A = np.zeros((len(Tlist), 2), np.float64)
            A[:, 0] = np.ones_like(Tlist)
            A[:, 1] = -1.0 / constants.R / Tlist
        b = np.log(klist)
        if weights is not None:
            for n in range(b.size):
                A[n, :] *= weights[n]
                b[n] *= weights[n]
        x, residues, rank, s = np.linalg.lstsq(A, b, rcond=RCOND)

        # Determine covarianace matrix to obtain parameter uncertainties
        count = klist.size
        cov = residues[0] / (count - 3) * np.linalg.inv(np.dot(A.T, A))
        t = scipy.stats.t.ppf(0.975, count - 3)

        if not three_params:
            x = np.array([x[0], 0, x[1]])
            cov = np.array([[cov[0, 0], 0, cov[0, 1]], [0, 0, 0], [cov[1, 0], 0, cov[1, 1]]])

        self.A = (rate_sign_multiplier * exp(x[0]), kunits)
        self.n = x[1]
        self.Ea = (x[2] * 0.001, "kJ/mol")
        self.T0 = (T0, "K")
        self.Tmin = (np.min(Tlist), "K")
        self.Tmax = (np.max(Tlist), "K")
        self.comment = 'Fitted to {0:d} data points; dA = *|/ {1:g}, dn = +|- {2:g}, dEa = +|- {3:g} kJ/mol'.format(
            len(Tlist),
            exp(sqrt(cov[0, 0])),
            sqrt(cov[1, 1]),
            sqrt(cov[2, 2]) * 0.001,
        )

        return self

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's Arrhenius.) Otherwise returns ``False``.
        """
        if not isinstance(other_kinetics, SurfaceChargeTransfer):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if (not self.A.equals(other_kinetics.A) or not self.n.equals(other_kinetics.n)
                or not self.Ea.equals(other_kinetics.Ea) or not self.T0.equals(other_kinetics.T0)
                or not self.alpha.equals(other_kinetics.alpha) or not self.electrons.equals(other_kinetics.electrons)
                or not self.V0.equals(other_kinetics.V0)):
            return False

        return True

    cpdef change_rate(self, double factor):
        """
        Changes A factor in Arrhenius expression by multiplying it by a ``factor``.
        """
        self._A.value_si *= factor

    cpdef SurfaceArrhenius to_surface_arrhenius(self):
        """
        Return an :class:`SurfaceArrhenius` instance of the kinetics model
        """
        return SurfaceArrhenius(
            A=self.A,
            n=self.n,
            Ea=self.Ea,
            T0=(1, "K"),
            Tmin=self.Tmin,
            Tmax=self.Tmax,
            uncertainty = self.uncertainty,
            comment=self.comment,
        )

    cpdef SurfaceChargeTransferBEP to_surface_charge_transfer_bep(self, double dGrxn, double V0=0.0):
        """
        Converts an SurfaceChargeTransfer object to SurfaceChargeTransferBEP
        """
        cdef double E0
    
        self.change_t0(1)
        self.change_v0(V0)

        E0 = self.Ea.value_si - self._alpha.value_si * dGrxn
        if E0 < 0:
            E0 = 0.0

        aep = SurfaceChargeTransferBEP(
                        A=self.A,
                        electrons=self.electrons,
                        n=self.n,
                        alpha=self.alpha,
                        V0=self.V0,
                        E0=(E0, 'J/mol'),
                        Tmin=self.Tmin,
                        Tmax=self.Tmax,
                        Pmin=self.Pmin,
                        Pmax=self.Pmax,
                        uncertainty=self.uncertainty,
                        comment=self.comment)
        return aep

cdef class SurfaceChargeTransferBEP(KineticsModel):
    """
    A kinetics model based on the (modified) Arrhenius equation, using the
    Evans-Polanyi equation to determine the activation energy. The attributes
    are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `n`             The temperature exponent
    `E0`            The activation energy at equilibiurm
    ` electrons`            The stochiometry coeff for electrons (negative if reactant, positive if product)
    `V0`            The reference potential
    `alpha`         The charge transfer coefficient
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, A=None, n=0.0, E0=None, V0=(0.0,'V'), alpha=0.5, electrons=-1, Tmin=None, Tmax=None,
                Pmin=None, Pmax=None, uncertainty=None, comment=''):

        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, uncertainty=uncertainty,
                comment=comment)

        self.alpha = alpha
        self.A = A
        self.n = n
        self.E0 = E0
        self.electrons =  electrons
        self.V0 = V0

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Arrhenius object.
        """
        string = 'SurfaceChargeTransferBEP(A={0!r}, n={1!r}, E0={2!r}, V0={3!r}, alpha={4!r}, electrons={5!r}'.format(
            self.A, self.n, self.E0, self.V0, self.alpha, self.electrons)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.uncertainty: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a SurfaceChargeTransfer object.
        """
        return (SurfaceChargeTransferBEP, (self.A, self.n, self.E0, self.V0, self.alpha, self.electrons, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                            self.uncertainty, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.SurfaceRateCoefficient(value)

    property n:
        """The temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property E0:
        """The activation energy."""
        def __get__(self):
            return self._E0
        def __set__(self, value):
            self._E0 = quantity.Energy(value)

    property V0:
        """The reference potential."""
        def __get__(self):
            return self._V0
        def __set__(self, value):
            self._V0 = quantity.Potential(value)

    property  electrons:
        """The number of electrons transferred."""
        def __get__(self):
            return self._electrons
        def __set__(self, value):
            self._electrons = quantity.Dimensionless(value)

    property alpha:
        """The charge transfer coefficient."""
        def __get__(self):
            return self._alpha
        def __set__(self, value):
            self._alpha = quantity.Dimensionless(value)

    cpdef change_v0(self, double V0):
        """
        Changes the reference potential to `V0` in volts, and adjusts the
        activation energy `E0` accordingly.
        """

        self._E0.value_si = self.get_activation_energy_from_potential(V0,0.0)
        self._V0.value_si = V0

    cpdef double get_activation_energy(self, double dGrxn) except -1:
        """
        Return the activation energy in J/mol corresponding to the given
        free energy of reaction `dGrxn` in J/mol at the reference potential.
        """
        cdef double Ea
        Ea = self._alpha.value_si * dGrxn + self._E0.value_si

        if Ea < 0.0:
            Ea = 0.0
        elif dGrxn > 0.0 and Ea < dGrxn:
            Ea = dGrxn

        return Ea

    cpdef double get_activation_energy_from_potential(self, double V, double dGrxn) except -1:
        """
        Return the activation energy in J/mol corresponding to the given
        free energy of reaction `dGrxn` in J/mol.
        """
        cdef double Ea
        Ea = self.get_activation_energy(dGrxn)
        Ea -= self._alpha.value_si * self._electrons.value_si * constants.F * (V-self._V0.value_si)

        return Ea

    cpdef double get_rate_coefficient_from_potential(self, double T, double V, double dGrxn) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K, potential `V` in volts, and 
        free of reaction `dGrxn` in J/mol.
        """
        cdef double A, n, Ea
        Ea = self.get_activation_energy_from_potential(V,dGrxn)
        A = self._A.value_si
        n = self._n.value_si
        return A * T ** n * exp(-Ea / (constants.R * T))

    cpdef SurfaceChargeTransfer to_surface_charge_transfer(self, double dGrxn):
        """
        Return an :class:`SurfaceChargeTransfer` instance of the kinetics model using the
        given free energy of reaction `dGrxn` to determine the activation energy.
        """
        return SurfaceChargeTransfer(
            A=self.A,
            n=self.n,
            electrons=self.electrons,
            Ea=(self.get_activation_energy(dGrxn) * 0.001, "kJ/mol"),
            V0=self.V0,
            T0=(1, "K"),
            Tmin=self.Tmin,
            Tmax=self.Tmax,
            Pmin=self.Pmin,
            Pmax=self.Pmax,
            uncertainty=self.uncertainty,
            comment=self.comment,
        )

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's Arrhenius.) Otherwise returns ``False``.
        """
        if not isinstance(other_kinetics, SurfaceChargeTransferBEP):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if (not self.A.equals(other_kinetics.A) or not self.n.equals(other_kinetics.n)
                or not self.E0.equals(other_kinetics.E0) or not self.alpha.equals(other_kinetics.alpha)
                or not self.electrons.equals(other_kinetics.electrons) or not self.V0.equals(other_kinetics.V0)):
            return False

        return True

    cpdef change_rate(self, double factor):
        """
        Changes A factor by multiplying it by a ``factor``.
        """
        self._A.value_si *= factor

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Sets a cantera ElementaryReaction() object with the modified Arrhenius object
        converted to an Arrhenius form.
        """
        raise NotImplementedError('set_cantera_kinetics() is not implemented for ArrheniusEP class kinetics.')
