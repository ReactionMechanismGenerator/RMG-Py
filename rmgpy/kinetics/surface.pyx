# cython: embedsignature=True, cdivision=True

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `T0`            The reference temperature
    `n`             The temperature exponent
    `Ea`            The activation energy
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `eps_ki`        The activation energy dependence on coverage
    `mu_ki`         The power-law exponent of coverage dependence
    `nu_ki`         The coefficient for exponential dependence on the coverage
    `species`       The species, i, that the coverage dependence parameters are dependent on
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================
    
    """

    def __init__(self, A=None, n=0.0, Ea=None, T0=(1.0, "K"), Tmin=None, Tmax=None, Pmin=None, Pmax=None, eps_ki=None,
                 mu_ki=None, nu_ki=None, species=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, eps_ki=eps_ki, mu_ki=mu_ki,
                               nu_ki=nu_ki, species=species, comment=comment)
        self.A = A
        self.n = n
        self.Ea = Ea
        self.T0 = T0
        self.eps_ki = eps_ki
        self.mu_ki = mu_ki
        self.nu_ki = nu_ki
        self.species = species

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        StickingCoefficient object.
        """
        string = 'StickingCoefficient(A={0!r}, n={1!r}, Ea={2!r}, T0={3!r}'.format(self.A, self.n, self.Ea, self.T0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.species is not None: string += f', species={self.species}'
        if self.eps_ki is not None: string += f', eps_ki={self.eps_ki}'
        if self.mu_ki is not None: string += f',  mu_ki={self.mu_ki}'
        if self.nu_ki is not None: string += f', nu_ki={self.nu_ki}'
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a StickingCoefficient object.
        """
        return (StickingCoefficient, (self.A, self.n, self.Ea, self.T0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                                      self.eps_ki, self.mu_ki, self.nu_ki, self.species, self.comment))

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

    property eps_ki:
        """The activation energy dependence on coverage."""
        def __get__(self):
            return self._eps_ki
        def __set__(self, value):
            self._eps_ki = quantity.Energy(value)

    property mu_ki:
        """The power-law exponent of coverage dependence."""
        def __get__(self):
            return self._mu_ki
        def __set__(self, value):
            self._mu_ki = quantity.Dimensionless(value)

    property nu_ki:
        """The coefficient for exponential dependence on the coverage."""
        def __get__(self):
            return self._nu_ki
        def __set__(self, value):
            self._nu_ki = quantity.Dimensionless(value)

    property species:
        """The species, i, that the coverage dependence parameters are dependent on."""
        def __get__(self):
            return self._species
        def __set__(self, value):
            # todo: not sure how to do this as species should be a rmg.species or a string?
            self._species = quantity.SurfaceConcentration

    cpdef double get_sticking_coefficient(self, double T, double species=1) except -1:
        """
        Return the sticking coefficient (dimensionless) at temperature `T` in K. 
        """
        cdef double A, n, Ea, T0, eps_ki, mu_ki, nu_ki, species, stickingCoefficient
        A = self._A.value_si
        n = self._n.value_si
        Ea = self._Ea.value_si
        T0 = self._T0.value_si
        eps_ki = self._eps_ki
        mu_ki = self._mu_ki
        nu_ki = self._nu_ki
        species = self._species  # todo: need to get the concentration/surface site fraction/surface coverage of said species
        if species:
            stickingCoefficient = A * (T / T0) ** n * exp(-Ea / (constants.R * T)) *
        else:
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
        )  # todo: add in coverage params here?

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

################################################################################
cdef class StickingCoefficientBEP(KineticsModel):
    """
    A kinetics model based on the Arrhenius expression, to give 
    Sticking Coefficient for surface adsorption, using the
    Bronsted-Evans-Polanyi equation to determine the activation energy. 
    Similar to :class:`ArrheniusEP`, but with different units for `A`.
    Sticking Coefficients are between 0 and 1.
    The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `n`             The temperature exponent
    `alpha`         The Evans-Polanyi slope
    `E0`            The activation energy for a thermoneutral reaction
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `eps_ki`        The activation energy dependence on coverage
    `mu_ki`         The power-law exponent of coverage dependence
    `nu_ki`         The coefficient for exponential dependence on the coverage
    `species`       The species, i, that the coverage dependence parameters are dependent on
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================
    
    """

    def __init__(self, A=None, n=0.0, alpha=0.0, E0=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, eps_ki=None,
                 mu_ki=None, nu_ki=None, species=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, eps_ki=eps_ki, mu_ki=mu_ki,
                               nu_ki=nu_ki, species=species, comment=comment)
        self.A = A
        self.n = n
        self.alpha = alpha
        self.E0 = E0
        self.eps_ki = eps_ki
        self.mu_ki = mu_ki
        self.nu_ki = nu_ki
        self.species = species

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        StickingCoefficientBEP object.
        """
        string = 'StickingCoefficientBEP(A={0!r}, n={1!r}, alpha={2!r}, E0={3!r}'.format(self.A, self.n, self.alpha,
                                                                                         self.E0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.species is not None: string += f', species={self.species}'
        if self.eps_ki is not None: string += f', eps_ki={self.eps_ki}'
        if self.mu_ki is not None: string += f',  mu_ki={self.mu_ki}'
        if self.nu_ki is not None: string += f', nu_ki={self.nu_ki}'
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an StickingCoefficientBEP object.
        """
        return (StickingCoefficientBEP, (self.A, self.n, self.Ea, self.T0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                                         self.eps_ki, self.mu_ki, self.nu_ki, self.species, self.comment))

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

    property eps_ki:
        """The activation energy dependence on coverage."""
        def __get__(self):
            return self._eps_ki
        def __set__(self, value):
            self._eps_ki = quantity.Energy(value)

    property mu_ki:
        """The power-law exponent of coverage dependence."""
        def __get__(self):
            return self._mu_ki
        def __set__(self, value):
            self._mu_ki = quantity.Dimensionless(value)

    property nu_ki:
        """The coefficient for exponential dependence on the coverage."""
        def __get__(self):
            return self._nu_ki
        def __set__(self, value):
            self._nu_ki = quantity.Dimensionless(value)

    property species:
        """The species, i, that the coverage dependence parameters are dependent on."""
        def __get__(self):
            return self._species
        def __set__(self, value):
            # todo: not sure how to do this as species should be a rmg.species or a string?
            self._species = quantity.SurfaceConcentration


    cpdef double get_sticking_coefficient(self, double T, double dHrxn=0.0) except -1:
        """
        Return the sticking coefficient (dimensionless) at
        temperature `T` in K and enthalpy of reaction `dHrxn` in J/mol. 
        """
        cdef double A, n, Ea, stickingCoefficient  # todo: update this function
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
            eps_ki=self.eps_ki,
            mu_ki=self.mu_ki,
            nu_ki=self.nu_ki,
            species=self.species,
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

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `T0`            The reference temperature
    `n`             The temperature exponent
    `Ea`            The activation energy
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `eps_ki`        The activation energy dependence on coverage
    `mu_ki`         The power-law exponent of coverage dependence
    `nu_ki`         The coefficient for exponential dependence on the coverage
    `species`       The species, i, that the coverage dependence parameters are dependent on
    `uncertainty`   Uncertainty information
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================
    """
    property A:  # todo: double check this is the only thing that changes even with coverage parameters
        """The preexponential factor. 
    
        This is the only thing different from a normal Arrhenius class."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.SurfaceRateCoefficient(value)

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        SurfaceArrhenius object.
        """
        string = 'SurfaceArrhenius(A={0!r}, n={1!r}, Ea={2!r}, T0={3!r}'.format(self.A, self.n, self.Ea, self.T0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.species is not None: string += f', species={self.species}'
        if self.eps_ki is not None: string += f', eps_ki={self.eps_ki}'
        if self.mu_ki is not None: string += f',  mu_ki={self.mu_ki}'
        if self.nu_ki is not None: string += f', nu_ki={self.nu_ki}'
        if self.uncertainty is not None: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a SurfaceArrhenius object.
        """
        return (SurfaceArrhenius, (self.A, self.n, self.Ea, self.T0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                                   self.eps_ki, self.mu_ki, self.nu_ki, self.species, self.uncertainty, self.comment))


################################################################################

cdef class SurfaceArrheniusBEP(ArrheniusEP):
    """
    A kinetics model based on the (modified) Arrhenius equation, using the
    Bronsted-Evans-Polanyi equation to determine the activation energy. 
    
    It is very similar to the gas-phase :class:`ArrheniusEP`.
    The only differences being the A factor has different units,
    (and the catalysis community prefers to call it BEP rather than EP!)
    
    The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `n`             The temperature exponent
    `alpha`         The Evans-Polanyi slope
    `E0`            The activation energy for a thermoneutral reaction
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `eps_ki`        The activation energy dependence on coverage
    `mu_ki`         The power-law exponent of coverage dependence
    `nu_ki`         The coefficient for exponential dependence on the coverage
    `species`       The species, i, that the coverage dependence parameters are dependent on
    `uncertainty`   Uncertainty information
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================
    
    """
    property A:  # todo: double check this is the only thing that changes even with coverage parameters
        """The preexponential factor. 
    
        This is the only thing different from a normal ArrheniusEP class."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.SurfaceRateCoefficient(value)

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        SurfaceArrheniusBEP object.
        """
        string = 'SurfaceArrheniusBEP(A={0!r}, n={1!r}, alpha={2!r}, E0={3!r}'.format(self.A, self.n, self.alpha,
                                                                                      self.E0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.species is not None: string += f', species={self.species}'
        if self.eps_ki is not None: string += f', eps_ki={self.eps_ki}'
        if self.mu_ki is not None: string += f',  mu_ki={self.mu_ki}'
        if self.nu_ki is not None: string += f', nu_ki={self.nu_ki}'
        if self.uncertainty is not None: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an SurfaceArrheniusBEP object.
        """
        return (SurfaceArrheniusBEP, (self.A, self.n, self.alpha, self.E0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                                   self.eps_ki, self.mu_ki, self.nu_ki, self.species, self.uncertainty, self.comment))

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
            eps_ki=self.eps_ki,
            mu_ki=self.mu_ki,
            nu_ki=self.nu_ki,
            species=self.species,
            uncertainty = self.uncertainty,
            comment=self.comment,
        )
