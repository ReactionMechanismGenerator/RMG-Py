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
from libc.math cimport exp, sqrt, log10
from scipy.optimize import curve_fit, fsolve

cimport rmgpy.constants as constants
import rmgpy.quantity as quantity
from rmgpy.exceptions import KineticsError
from rmgpy.kinetics.uncertainties import rank_accuracy_map
from rmgpy.molecule.molecule import Bond
import logging

# Prior to numpy 1.14, `numpy.linalg.lstsq` does not accept None as a value
RCOND = -1 if int(np.__version__.split('.')[1]) < 14 else None
################################################################################

cdef class Arrhenius(KineticsModel):
    """
    A kinetics model based on the (modified) Arrhenius equation. The attributes
    are:

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
    `solute`        Transition state solute data
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, A=None, n=0.0, Ea=None, T0=(1.0, "K"), Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 uncertainty=None, solute=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, uncertainty=uncertainty,
                               solute=solute, comment=comment)
        self.A = A
        self.n = n
        self.Ea = Ea
        self.T0 = T0

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Arrhenius object.
        """
        string = 'Arrhenius(A={0!r}, n={1!r}, Ea={2!r}, T0={3!r}'.format(self.A, self.n, self.Ea, self.T0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.uncertainty: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.solute: string += ', solute={0!r}'.format(self.solute)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an Arrhenius object.
        """
        return (Arrhenius, (self.A, self.n, self.Ea, self.T0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                            self.uncertainty, self.solute, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.RateCoefficient(value)

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

    cpdef double get_rate_coefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K.
        """
        cdef double A, n, Ea, T0
        A = self._A.value_si
        n = self._n.value_si
        Ea = self._Ea.value_si
        T0 = self._T0.value_si
        return A * (T / T0) ** n * exp(-Ea / (constants.R * T))

    cpdef change_t0(self, double T0):
        """
        Changes the reference temperature used in the exponent to `T0` in K,
        and adjusts the preexponential factor accordingly.
        """
        self._A.value_si /= (self._T0.value_si / T0) ** self._n.value_si
        self._T0.value_si = T0

    cpdef fit_to_data(self, np.ndarray Tlist, np.ndarray klist, str kunits, double T0=1, np.ndarray weights=None,
                      bint three_params=True):
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
        self.solute = None
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
        if not isinstance(other_kinetics, Arrhenius):
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

    def to_cantera_kinetics(self):
        """
        Converts the Arrhenius object to a cantera Arrhenius object

        Arrhenius(A,b,E) where A is in units of m^3/kmol/s, b is dimensionless, and E is in J/kmol
        """

        import cantera as ct

        rate_units_dimensionality = {'1/s': 0,
                                     's^-1': 0,
                                     'm^3/(mol*s)': 1,
                                     'm^6/(mol^2*s)': 2,
                                     'cm^3/(mol*s)': 1,
                                     'cm^6/(mol^2*s)': 2,
                                     'm^3/(molecule*s)': 1,
                                     'm^6/(molecule^2*s)': 2,
                                     'cm^3/(molecule*s)': 1,
                                     'cm^6/(molecule^2*s)': 2,
                                     }

        if self._T0.value_si != 1:
            A = self._A.value_si / (self._T0.value_si) ** self._n.value_si
        else:
            A = self._A.value_si

        try:
            A *= 1000 ** rate_units_dimensionality[self._A.units]
        except KeyError:
            raise Exception('Arrhenius A-factor units {0} not found among accepted units for converting to '
                            'Cantera Arrhenius object.'.format(self._A.units))

        b = self._n.value_si
        E = self._Ea.value_si * 1000  # convert from J/mol to J/kmol
        return ct.Arrhenius(A, b, E)

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Passes in a cantera ElementaryReaction() object and sets its
        rate to a Cantera Arrhenius() object.
        """
        import cantera as ct
        assert isinstance(ct_reaction, ct.ElementaryReaction), "Must be a Cantera ElementaryReaction object"

        # Set the rate parameter to a cantera Arrhenius object
        ct_reaction.rate = self.to_cantera_kinetics()

    cpdef ArrheniusEP to_arrhenius_ep(self, double alpha=0.0, double dHrxn=0.0):
        """
        Converts an Arrhenius object to ArrheniusEP

        If setting alpha, you need to also input dHrxn, which must be given
        in J/mol (and vise versa).
        """

        if bool(alpha) ^ bool(dHrxn):
            raise Exception('If you set alpha or dHrxn in to_arrhenius_ep, '
                            'you need to set the other value to non-zero.')
        self.change_t0(1)
        aep = ArrheniusEP(A=self.A,
                          n=self.n,
                          alpha=alpha,
                          E0=(self.Ea.value_si - alpha * dHrxn, 'J/mol'),
                          Tmin=self.Tmin,
                          Tmax=self.Tmax,
                          Pmin=self.Pmin,
                          Pmax=self.Pmax,
                          uncertainty=self.uncertainty,
                          solute=self.solute,
                          comment=self.comment)
        return aep
################################################################################


cdef class ArrheniusEP(KineticsModel):
    """
    A kinetics model based on the (modified) Arrhenius equation, using the
    Evans-Polanyi equation to determine the activation energy. The attributes
    are:

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
    `solute`        Transition state solute data
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, A=None, n=0.0, alpha=0.0, E0=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, uncertainty=None,
                 solute=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, uncertainty=uncertainty,
                               solute=solute, comment=comment)
        self.A = A
        self.n = n
        self.alpha = alpha
        self.E0 = E0

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        ArrheniusEP object.
        """
        string = 'ArrheniusEP(A={0!r}, n={1!r}, alpha={2!r}, E0={3!r}'.format(self.A, self.n, self.alpha, self.E0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.uncertainty is not None: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.solute is not None: string += ', solute={0!r}'.format(self.solute)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an ArrheniusEP object.
        """
        return (ArrheniusEP, (self.A, self.n, self.alpha, self.E0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                              self.uncertainty, self.solute, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.RateCoefficient(value)

    property n:
        """The temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property alpha:
        """The Evans-Polanyi slope."""
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

    cpdef double get_rate_coefficient(self, double T, double dHrxn=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K and enthalpy of reaction `dHrxn`
        in J/mol.
        """
        cdef double A, n, Ea
        Ea = self.get_activation_energy(dHrxn)
        A = self._A.value_si
        n = self._n.value_si
        return A * T ** n * exp(-Ea / (constants.R * T))

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

    cpdef Arrhenius to_arrhenius(self, double dHrxn):
        """
        Return an :class:`Arrhenius` instance of the kinetics model using the
        given enthalpy of reaction `dHrxn` to determine the activation energy.
        """
        return Arrhenius(
            A=self.A,
            n=self.n,
            Ea=(self.get_activation_energy(dHrxn) * 0.001, "kJ/mol"),
            T0=(1, "K"),
            Tmin=self.Tmin,
            Tmax=self.Tmax,
            Pmin=self.Pmin,
            Pmax=self.Pmax,
            uncertainty=self.uncertainty,
            solute=self.solute,
            comment=self.comment,
        )

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's Arrhenius.) Otherwise returns ``False``.
        """
        if not isinstance(other_kinetics, ArrheniusEP):
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

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Sets a cantera ElementaryReaction() object with the modified Arrhenius object
        converted to an Arrhenius form.
        """
        raise NotImplementedError('set_cantera_kinetics() is not implemented for ArrheniusEP class kinetics.')

################################################################################

cdef class ArrheniusBM(KineticsModel):
    """
    A kinetics model based on the (modified) Arrhenius equation, using the
    Blowers-Masel equation to determine the activation energy.
    Based on Blowers and Masel's 2000 paper Engineering Approximations for Activation
    Energies in Hydrogen Transfer Reactions.
    The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `n`             The temperature exponent
    `w0`            The average of the bond dissociation energies of the bond formed and the bond broken
    `E0`            The activation energy for a thermoneutral reaction
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `solute`        Transition state solute data
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, A=None, n=0.0, w0=(0.0, 'J/mol'), E0=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 uncertainty=None, solute=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, uncertainty=uncertainty,
                               solute=solute, comment=comment)
        self.A = A
        self.n = n
        self.w0 = w0
        self.E0 = E0

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        ArrheniusBM object.
        """
        string = 'ArrheniusBM(A={0!r}, n={1!r}, w0={2!r}, E0={3!r}'.format(self.A, self.n, self.w0, self.E0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.uncertainty is not None: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.solute is not None: string += ', solute={0!r}'.format(self.solute)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an ArrheniusEP object.
        """
        return (ArrheniusBM, (self.A, self.n, self.w0, self.E0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                              self.uncertainty, self.solute, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            return self._A
        def __set__(self, value):
            self._A = quantity.RateCoefficient(value)

    property n:
        """The temperature exponent."""
        def __get__(self):
            return self._n
        def __set__(self, value):
            self._n = quantity.Dimensionless(value)

    property w0:
        """The average of the bond dissociation energies of the bond formed and the bond broken."""
        def __get__(self):
            return self._w0
        def __set__(self, value):
            self._w0 = quantity.Energy(value)

    property E0:
        """The activation energy for a thermoneutral reaction."""
        def __get__(self):
            return self._E0
        def __set__(self, value):
            self._E0 = quantity.Energy(value)

    cpdef double get_rate_coefficient(self, double T, double dHrxn=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K and enthalpy of reaction `dHrxn`
        in J/mol.
        """
        cdef double A, n, Ea
        Ea = self.get_activation_energy(dHrxn)
        A = self._A.value_si
        n = self._n.value_si
        return A * T ** n * exp(-Ea / (constants.R * T))

    cpdef double get_activation_energy(self, double dHrxn) except -1:
        """
        Return the activation energy in J/mol corresponding to the given
        enthalpy of reaction `dHrxn` in J/mol.
        """
        cdef double w0, E0
        E0 = self._E0.value_si
        if dHrxn < -4 * self._E0.value_si:
            return 0.0
        elif dHrxn > 4 * self._E0.value_si:
            return dHrxn
        else:
            w0 = self._w0.value_si
            Vp = 2 * w0 * (2 * w0 + 2 * E0) / (2 * w0 - 2 * E0)
            return (w0 + dHrxn / 2.0) * (Vp - 2 * w0 + dHrxn) ** 2 / (Vp ** 2 - (2 * w0) ** 2 + dHrxn ** 2)

    cpdef Arrhenius to_arrhenius(self, double dHrxn):
        """
        Return an :class:`Arrhenius` instance of the kinetics model using the
        given enthalpy of reaction `dHrxn` to determine the activation energy.
        """
        return Arrhenius(
            A=self.A,
            n=self.n,
            Ea=(self.get_activation_energy(dHrxn) * 0.001, "kJ/mol"),
            T0=(1, "K"),
            Tmin=self.Tmin,
            Tmax=self.Tmax,
            uncertainty=self.uncertainty,
            solute=self.solute,
            comment=self.comment,
        )

    def fit_to_reactions(self, rxns, w0=None, recipe=None, Ts=None):
        """
        Fit an ArrheniusBM model to a list of reactions at the given temperatures,
        w0 must be either given or estimated using the family object
        """
        assert w0 is not None or recipe is not None, 'either w0 or recipe must be specified'

        if Ts is None:
            Ts = [300.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1500.0, 2000.0]
        if w0 is None:
            #estimate w0
            w0s = get_w0s(recipe, rxns)
            w0 = sum(w0s) / len(w0s)

        if len(rxns) == 1:
            T = 1000.0
            rxn = rxns[0]
            dHrxn = rxn.get_enthalpy_of_reaction(298.0)
            A = rxn.kinetics.A.value_si
            n = rxn.kinetics.n.value_si
            Ea = rxn.kinetics.Ea.value_si

            def kfcn(E0):
                Vp = 2 * w0 * (2 * w0 + 2 * E0) / (2 * w0 - 2 * E0)
                out = Ea - (w0 + dHrxn / 2.0) * (Vp - 2 * w0 + dHrxn) * (Vp - 2 * w0 + dHrxn) / (Vp * Vp - (2 * w0) * (2 * w0) + dHrxn * dHrxn)
                return out

            E0 = fsolve(kfcn, w0 / 10.0)[0]

            self.Tmin = rxn.kinetics.Tmin
            self.Tmax = rxn.kinetics.Tmax
            self.solute = None
            self.comment = 'Fitted to {0} reaction at temperature: {1} K'.format(len(rxns), T)
        else:
            # define optimization function
            def kfcn(xs, lnA, n, E0):
                T = xs[:,0]
                dHrxn = xs[:,1]
                Vp = 2 * w0 * (2 * w0 + 2 * E0) / (2 * w0 - 2 * E0)
                Ea = (w0 + dHrxn / 2.0) * (Vp - 2 * w0 + dHrxn) * (Vp - 2 * w0 + dHrxn) / (Vp * Vp - (2 * w0) * (2 * w0) + dHrxn * dHrxn)
                Ea = np.where(dHrxn< -4.0*E0, 0.0, Ea)
                Ea = np.where(dHrxn > 4.0*E0, dHrxn, Ea)
                return lnA + np.log(T ** n * np.exp(-Ea / (8.314 * T)))

            # get (T,dHrxn(T)) -> (Ln(k) mappings
            xdata = []
            ydata = []
            sigmas = []
            for rxn in rxns:
                # approximately correct the overall uncertainties to std deviations
                s = rank_accuracy_map[rxn.rank].value_si/2.0
                for T in Ts:
                    xdata.append([T, rxn.get_enthalpy_of_reaction(298.0)])
                    ydata.append(np.log(rxn.get_rate_coefficient(T)))

                    sigmas.append(s / (8.314 * T))

            xdata = np.array(xdata)
            ydata = np.array(ydata)

            # fit parameters
            boo = True
            xtol = 1e-8
            ftol = 1e-8
            while boo:
                boo = False
                try:
                    params = curve_fit(kfcn, xdata, ydata, sigma=sigmas, p0=[1.0, 1.0, w0 / 10.0], xtol=xtol, ftol=ftol)
                except RuntimeError:
                    if xtol < 1.0:
                        boo = True
                        xtol *= 10.0
                        ftol *= 10.0
                    else:
                        raise ValueError("Could not fit BM arrhenius to reactions with xtol<1.0")

            lnA, n, E0 = params[0].tolist()
            A = np.exp(lnA)

            self.Tmin = (np.min(Ts), "K")
            self.Tmax = (np.max(Ts), "K")
            self.solute = None
            self.comment = 'Fitted to {0} reactions at temperatures: {1}'.format(len(rxns), Ts)

        # fill in parameters
        A_units = ['', 's^-1', 'm^3/(mol*s)', 'm^6/(mol^2*s)']
        order = len(rxns[0].reactants)
        self.A = (A, A_units[order])

        self.n = n
        self.w0 = (w0, 'J/mol')
        self.E0 = (E0, 'J/mol')

        return self

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's Arrhenius.) Otherwise returns ``False``.
        """
        if not isinstance(other_kinetics, ArrheniusBM):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if (not self.A.equals(other_kinetics.A) or not self.n.equals(other_kinetics.n)
                or not self.w0.equals(other_kinetics.w0) or not self.E0.equals(other_kinetics.E0)):
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
        raise NotImplementedError('set_cantera_kinetics() is not implemented for ArrheniusBM class kinetics.')

################################################################################

cdef class PDepArrhenius(PDepKineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient :math:`k(T,P)` where
    a set of Arrhenius kinetics are stored at a variety of pressures and
    interpolated between on a logarithmic scale. The attributes are:

    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `pressures`     The list of pressures
    `arrhenius`     The list of :class:`Arrhenius` objects at each pressure
    `Tmin`          The minimum temperature in K at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature in K at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure in bar at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure in bar at which the model is valid, or zero if unknown or undefined
    `efficiencies`  A dict associating chemical species with associated efficiencies
    `order`         The reaction order (1 = first, 2 = second, etc.)
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================

    """

    def __init__(self, pressures=None, arrhenius=None, highPlimit=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 comment=''):
        PDepKineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, highPlimit=highPlimit,
                                   comment=comment)
        self.pressures = pressures
        self.arrhenius = arrhenius or []

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        PDepArrhenius object.
        """
        string = 'PDepArrhenius(pressures={0!r}, arrhenius={1!r}'.format(self.pressures, self.arrhenius)
        if self.highPlimit is not None: string += ', highPlimit={0!r}'.format(self.highPlimit)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a PDepArrhenius object.
        """
        return (PDepArrhenius, (self.pressures, self.arrhenius, self.highPlimit, self.Tmin, self.Tmax,
                                self.Pmin, self.Pmax, self.comment))

    property pressures:
        """The list of pressures."""
        def __get__(self):
            return self._pressures
        def __set__(self, value):
            self._pressures = quantity.Pressure(value)

    cdef get_adjacent_expressions(self, double P):
        """
        Returns the pressures and Arrhenius expressions for the pressures that
        most closely bound the specified pressure `P` in Pa.
        """
        cdef np.ndarray[np.float64_t, ndim=1] pressures
        cdef int i, ilow, ihigh

        pressures = self._pressures.value_si

        ilow = 0
        ihigh = -1
        for i in range(pressures.shape[0]):
            if pressures[i] <= P:
                ilow = i
            if pressures[i] >= P and ihigh == -1:
                ihigh = i
        return pressures[ilow], pressures[ihigh], self.arrhenius[ilow], self.arrhenius[ihigh]

    cpdef double get_rate_coefficient(self, double T, double P=0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K and pressure `P` in Pa.
        """
        cdef double Plow, Phigh, klow, khigh, k
        cdef KineticsModel alow, ahigh
        cdef int j

        if P == 0:
            raise ValueError('No pressure specified to pressure-dependent PDepArrhenius.get_rate_coefficient().')

        k = 0.0
        Plow, Phigh, alow, ahigh = self.get_adjacent_expressions(P)
        if Plow == Phigh:
            k = alow.get_rate_coefficient(T)
        else:
            klow = alow.get_rate_coefficient(T)
            khigh = ahigh.get_rate_coefficient(T)
            if klow == khigh == 0.0: return 0.0
            k = klow * 10 ** (log10(P / Plow) / log10(Phigh / Plow) * log10(khigh / klow))
        return k

    cpdef fit_to_data(self, np.ndarray Tlist, np.ndarray Plist, np.ndarray K, str kunits, double T0=1):
        """
        Fit the pressure-dependent Arrhenius model to a matrix of rate
        coefficient data `K` with units of `kunits` corresponding to a set of
        temperatures `Tlist` in K and pressures `Plist` in Pa. An Arrhenius
        model is fit cpdef change_rate(self, double factor)at each pressure.
        """
        cdef int i
        self.pressures = (Plist * 1e-5, "bar")
        self.arrhenius = []
        for i in range(len(Plist)):
            arrhenius = Arrhenius().fit_to_data(Tlist, K[:, i], kunits, T0)
            self.arrhenius.append(arrhenius)
        return self

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Each duplicate
        reaction must be matched and equal to that in the other PDepArrhenius model
        in the same order.  Otherwise returns ``False``
        """
        if not isinstance(other_kinetics, PDepArrhenius):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if len(self.arrhenius) != len(other_kinetics.arrhenius):
            return False
        if not self.pressures.equals(other_kinetics.pressures):
            return False
        for index in range(len(self.arrhenius)):
            if not self.arrhenius[index].is_identical_to(other_kinetics.arrhenius[index]):
                return False
        if self.highPlimit and not self.highPlimit.equals(other_kinetics.highPlimit):
            return False

        return True

    cpdef change_rate(self, double factor):
        """
        Changes kinetics rate by a multiple ``factor``.
        """
        for kin in self.arrhenius:
            kin.change_rate(factor)
        if self.highPlimit is not None:
            self.highPlimit.change_rate(factor)

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Sets a Cantera PlogReaction()'s `rates` attribute with
        A list of tuples containing [(pressure in Pa, cantera arrhenius object), (..)]
        """
        import cantera as ct
        import copy
        assert isinstance(ct_reaction, ct.PlogReaction), "Must be a Cantera PlogReaction object"

        pressures = copy.deepcopy(self._pressures.value_si)
        ctArrhenius = [arr.to_cantera_kinetics() for arr in self.arrhenius]

        ct_reaction.rates = list(zip(pressures, ctArrhenius))

################################################################################

cdef class MultiArrhenius(KineticsModel):
    """
    A kinetics model based on a set of (modified) Arrhenius equations, which
    are summed to obtain the overall rate. The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `arrhenius`     A list of the :class:`Arrhenius` kinetics
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, arrhenius=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.arrhenius = arrhenius

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        MultiArrhenius object.
        """
        string = 'MultiArrhenius(arrhenius={0!r}'.format(self.arrhenius)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an MultiArrhenius object.
        """
        return (MultiArrhenius, (self.arrhenius, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    cpdef double get_rate_coefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K.
        """
        cdef double k
        cdef Arrhenius arrh
        k = 0.0
        for arrh in self.arrhenius:
            k += arrh.get_rate_coefficient(T)
        return k

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Each duplicate
        reaction must be matched and equal to that in the other MultiArrhenius model
        in the same order.  Otherwise returns ``False``
        """
        if not isinstance(other_kinetics, MultiArrhenius):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if len(self.arrhenius) != len(other_kinetics.arrhenius):
            return False

        for index in range(len(self.arrhenius)):
            if not self.arrhenius[index].is_identical_to(other_kinetics.arrhenius[index]):
                return False

        return True

    cpdef Arrhenius to_arrhenius(self, double Tmin=-1, double Tmax=-1):
        """
        Return an :class:`Arrhenius` instance of the kinetics model

        Fit the Arrhenius parameters to a set of rate coefficient data generated
        from the MultiArrhenius kinetics, over the temperature range
        Tmin to Tmax, in Kelvin. If Tmin or Tmax are unspecified (or -1)
        then the MultiArrhenius's Tmin and Tmax are used.
        A linear least-squares fit is used, which guarantees that the
        resulting parameters provide the best possible approximation to the
        data.
        """
        cdef Arrhenius arrh
        cdef np.ndarray Tlist, klist
        cdef str kunits
        if Tmin == -1: Tmin = self.Tmin.value_si
        if Tmax == -1: Tmax = self.Tmax.value_si
        kunits = str(quantity.pq.Quantity(1.0, self.arrhenius[0].A.units).simplified).split()[-1]  # is this the best way to get the units returned by k??
        Tlist = np.logspace(log10(Tmin), log10(Tmax), num=25)
        klist = np.array(list(map(self.get_rate_coefficient, Tlist)), np.float64)
        arrh = Arrhenius().fit_to_data(Tlist, klist, kunits)
        arrh.comment = "Fitted to Multiple Arrhenius kinetics over range {Tmin}-{Tmax} K. {comment}".format(
            Tmin=Tmin, Tmax=Tmax, comment=self.comment)
        return arrh

    cpdef change_rate(self, double factor):
        """
        Change kinetics rate by a multiple ``factor``.
        """
        for kin in self.arrhenius:
            kin.change_rate(factor)

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Sets the kinetic rates for a list of cantera `Reaction` objects
        Here, ct_reaction must be a list rather than a single cantera reaction.
        """
        if len(ct_reaction) != len(self.arrhenius):
            raise Exception('The number of Cantera Reaction objects does not match the number of Arrhenius objects')

        for i, arr in enumerate(self.arrhenius):
            arr.set_cantera_kinetics(ct_reaction[i], species_list)

################################################################################

cdef class MultiPDepArrhenius(PDepKineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient :math:`k(T,P)` where
    sets of Arrhenius kinetics are stored at a variety of pressures and
    interpolated between on a logarithmic scale. The attributes are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `arrhenius`     A list of the :class:`PDepArrhenius` kinetics at each temperature
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, arrhenius=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        PDepKineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.arrhenius = arrhenius

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        MultiPDepArrhenius object.
        """
        string = 'MultiPDepArrhenius(arrhenius={0!r}'.format(self.arrhenius)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an MultiPDepArrhenius object.
        """
        return (MultiPDepArrhenius, (self.arrhenius, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    cpdef double get_rate_coefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K and pressure `P` in Pa.
        """
        cdef double k, klow, khigh, Plow, Phigh
        cdef PDepArrhenius arrh
        cdef KineticsModel arrh_low, arrh_high
        cdef np.ndarray Plist1, Plist2
        cdef int i

        if P == 0:
            raise ValueError('No pressure specified to pressure-dependent MultiPDepArrhenius.get_rate_coefficient().')

        k = 0
        for arrh in self.arrhenius:
            k += arrh.get_rate_coefficient(T, P)

        return k

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Each duplicate
        reaction must be matched and equal to that in the other MultiArrhenius model
        in the same order.  Otherwise returns ``False``
        """
        if not isinstance(other_kinetics, MultiPDepArrhenius):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if len(self.arrhenius) != len(other_kinetics.arrhenius):
            return False

        for index in range(len(self.arrhenius)):
            if not self.arrhenius[index].is_identical_to(other_kinetics.arrhenius[index]):
                return False

        return True

    cpdef change_rate(self, double factor):
        """
        Change kinetic rate by a multiple ``factor``.
        """
        for kin in self.arrhenius:
            kin.change_rate(factor)

    def set_cantera_kinetics(self, ct_reaction, species_list):
        """
        Sets the PLOG kinetics for multiple cantera `Reaction` objects, provided in a list.
        ct_reaction is a list of cantera reaction objects.
        """
        if len(ct_reaction) != len(self.arrhenius):
            raise Exception('The number of Cantera Reaction objects does not match the number of PdepArrhenius objects')

        for i, arr in enumerate(self.arrhenius):
            arr.set_cantera_kinetics(ct_reaction[i], species_list)

################################################################################

cdef class ArrheniusChargeTransfer(KineticsModel):

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
    `solute`        The transition state solute data
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, A=None, n=0.0, Ea=None, V0=None, alpha=0.5, electrons=-1, T0=(1.0, "K"), Tmin=None, Tmax=None,
                Pmin=None, Pmax=None, solute=None, uncertainty=None, comment=''):

        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, solute=solute, uncertainty=uncertainty,
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
        string = 'ArrheniusChargeTransfer(A={0!r}, n={1!r}, Ea={2!r}, V0={3!r}, alpha={4!r}, electrons={5!r}, T0={6!r}'.format(
            self.A, self.n, self.Ea, self.V0, self.alpha, self.electrons, self.T0)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.solute: string += ', solute={0!r}'.format(self.solute)
        if self.uncertainty: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a ArrheniusChargeTransfer object.
        """
        return (ArrheniusChargeTransfer, (self.A, self.n, self.Ea, self.V0, self.alpha, self.electrons, self.T0, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                            self.solute, self.uncertainty, self.comment))

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
        self.solute = None,
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
        if not isinstance(other_kinetics, ArrheniusChargeTransfer):
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

cdef class ArrheniusChargeTransferBM(KineticsModel):
    """
    A kinetics model based on the (modified) Arrhenius equation, using the
    Evans-Polanyi equation to determine the activation energy. The attributes
    are:

    =============== =============================================================
    Attribute       Description
    =============== =============================================================
    `A`             The preexponential factor
    `n`             The temperature exponent
    `w0`            The average of the bond dissociation energies of the bond formed and the bond broken
    `E0`            The activation energy for a thermoneutral reaction
    `electrons`     The stochiometry coeff for electrons (negative if reactant, positive if product)
    `V0`            The reference potential
    `alpha`         The charge transfer coefficient
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `solute`        Transition state solute data
    `comment`       Information about the model (e.g. its source)
    =============== =============================================================

    """

    def __init__(self, A=None, n=0.0, w0=(0.0, 'J/mol'), E0=None, V0=(0.0,'V'), alpha=0.5, electrons=-1, Tmin=None, Tmax=None,
                Pmin=None, Pmax=None, solute=None, uncertainty=None, comment=''):

        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, solute=solute, uncertainty=uncertainty,
                 comment=comment)

        self.alpha = alpha
        self.A = A
        self.n = n
        self.w0 = w0
        self.E0 = E0
        self.electrons =  electrons
        self.V0 = V0

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Arrhenius object.
        """
        string = 'ArrheniusChargeTransferBM(A={0!r}, n={1!r}, w0={2!r}, E0={3!r}, V0={4!r}, alpha={5!r}, electrons={6!r}'.format(
            self.A, self.n, self.w0, self.E0, self.V0, self.alpha, self.electrons)
        if self.Tmin is not None: string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None: string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None: string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None: string += ', Pmax={0!r}'.format(self.Pmax)
        if self.solute: string += ', solute={0!r}'.format(self.solute)
        if self.uncertainty: string += ', uncertainty={0!r}'.format(self.uncertainty)
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a ArrheniusChargeTransfer object.
        """
        return (ArrheniusChargeTransferBM, (self.A, self.n, self.w0, self.E0, self.V0, self.alpha, self.electrons, self.Tmin, self.Tmax, self.Pmin, self.Pmax,
                            self.solute, self.uncertainty, self.comment))

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

    property w0:
        """The average of the bond dissociation energies of the bond formed and the bond broken."""
        def __get__(self):
            return self._w0
        def __set__(self, value):
            self._w0 = quantity.Energy(value)

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

    cpdef double get_rate_coefficient(self, double T, double dHrxn=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K and enthalpy of reaction `dHrxn`
        in J/mol.
        """
        cdef double A, n, Ea
        Ea = self.get_activation_energy(dHrxn)
        A = self._A.value_si
        n = self._n.value_si
        return A * T ** n * exp(-Ea / (constants.R * T))

    cpdef double get_activation_energy(self, double dHrxn) except -1:
        """
        Return the activation energy in J/mol corresponding to the given
        enthalpy of reaction `dHrxn` in J/mol.
        """
        cdef double w0, E0
        E0 = self._E0.value_si
        if dHrxn < -4 * self._E0.value_si:
            return 0.0
        elif dHrxn > 4 * self._E0.value_si:
            return dHrxn
        else:
            w0 = self._w0.value_si
            Vp = 2 * w0 * (2 * w0 + 2 * E0) / (2 * w0 - 2 * E0)
            return (w0 + dHrxn / 2.0) * (Vp - 2 * w0 + dHrxn) ** 2 / (Vp ** 2 - (2 * w0) ** 2 + dHrxn ** 2)

    cpdef double get_rate_coefficient_from_potential(self, double T, double V, double dHrxn) except -1:
        """
        Return the rate coefficient in the appropriate combination of m^3,
        mol, and s at temperature `T` in K, potential `V` in volts, and
        heat of reaction `dHrxn` in J/mol.
        """
        cdef double A, n, Ea
        Ea = self.get_activation_energy_from_potential(V,dHrxn)
        Ea -= self._alpha.value_si * self._electrons.value_si * constants.F * (V-self._V0.value_si)
        A = self._A.value_si
        n = self._n.value_si
        return A * T ** n * exp(-Ea / (constants.R * T))

    def fit_to_reactions(self, rxns, w0=None, recipe=None, Ts=None):
        """
        Fit an ArrheniusBM model to a list of reactions at the given temperatures,
        w0 must be either given or estimated using the family object
        """
        assert w0 is not None or recipe is not None, 'either w0 or recipe must be specified'

        for rxn in rxns:
            if rxn.kinetics._V0.value_si != 0.0:
                rxn.kinetics.change_v0(0.0)

        if Ts is None:
            Ts = [300.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1500.0]
        if w0 is None:
            #estimate w0
            w0s = get_w0s(recipe, rxns)
            w0 = sum(w0s) / len(w0s)

        if len(rxns) == 1:
            rxn = rxns[0]
            dHrxn = rxn.get_enthalpy_of_reaction(298.0)
            A = rxn.kinetics.A.value_si
            n = rxn.kinetics.n.value_si
            Ea = rxn.kinetics.Ea.value_si

            def kfcn(E0):
                Vp = 2 * w0 * (2 * w0 + 2 * E0) / (2 * w0 - 2 * E0)
                out = Ea - (w0 + dHrxn / 2.0) * (Vp - 2 * w0 + dHrxn) * (Vp - 2 * w0 + dHrxn) / (Vp * Vp - (2 * w0) * (2 * w0) + dHrxn * dHrxn)
                return out

            E0 = fsolve(kfcn, w0 / 10.0)[0]

            self.Tmin = rxn.kinetics.Tmin
            self.Tmax = rxn.kinetics.Tmax
            self.comment = 'Fitted to 1 reaction'
        else:
            # define optimization function
            def kfcn(xs, lnA, n, E0):
                T = xs[:,0]
                dHrxn = xs[:,1]
                Vp = 2 * w0 * (2 * w0 + 2 * E0) / (2 * w0 - 2 * E0)
                Ea = (w0 + dHrxn / 2.0) * (Vp - 2 * w0 + dHrxn) * (Vp - 2 * w0 + dHrxn) / (Vp * Vp - (2 * w0) * (2 * w0) + dHrxn * dHrxn)
                Ea = np.where(dHrxn< -4.0*E0, 0.0, Ea)
                Ea = np.where(dHrxn > 4.0*E0, dHrxn, Ea)
                return lnA + np.log(T ** n * np.exp(-Ea / (8.314 * T)))

            # get (T,dHrxn(T)) -> (Ln(k) mappings
            xdata = []
            ydata = []
            sigmas = []
            for rxn in rxns:
                # approximately correct the overall uncertainties to std deviations
                s = rank_accuracy_map[rxn.rank].value_si/2.0
                for T in Ts:
                    xdata.append([T, rxn.get_enthalpy_of_reaction(298.0)])
                    ydata.append(np.log(rxn.get_rate_coefficient(T)))

                    sigmas.append(s / (8.314 * T))

            xdata = np.array(xdata)
            ydata = np.array(ydata)

            # fit parameters
            boo = True
            xtol = 1e-8
            ftol = 1e-8
            while boo:
                boo = False
                try:
                    params = curve_fit(kfcn, xdata, ydata, sigma=sigmas, p0=[1.0, 1.0, w0 / 10.0], xtol=xtol, ftol=ftol)
                except RuntimeError:
                    if xtol < 1.0:
                        boo = True
                        xtol *= 10.0
                        ftol *= 10.0
                    else:
                        raise ValueError("Could not fit BM arrhenius to reactions with xtol<1.0")

            lnA, n, E0 = params[0].tolist()
            A = np.exp(lnA)

            self.Tmin = (np.min(Ts), "K")
            self.Tmax = (np.max(Ts), "K")
            self.comment = 'Fitted to {0} reactions at temperatures: {1}'.format(len(rxns), Ts)

        # fill in parameters
        A_units = ['', 's^-1', 'm^3/(mol*s)', 'm^6/(mol^2*s)']
        order = len(rxns[0].reactants)
        self.A = (A, A_units[order])

        self.n = n
        self.w0 = (w0, 'J/mol')
        self.E0 = (E0, 'J/mol')
        self._V0.value_si = 0.0
        self.electrons = rxns[0].electrons

        return self

    cpdef ArrheniusChargeTransfer to_arrhenius_charge_transfer(self, double dHrxn):
        """
        Return an :class:`ArrheniusChargeTransfer` instance of the kinetics model using the
        given heat of reaction `dHrxn` to determine the activation energy.
        """
        return ArrheniusChargeTransfer(
            A=self.A,
            n=self.n,
            electrons=self.electrons,
            Ea=(self.get_activation_energy(dHrxn) * 0.001, "kJ/mol"),
            V0=self.V0,
            T0=(1, "K"),
            Tmin=self.Tmin,
            Tmax=self.Tmax,
            Pmin=self.Pmin,
            Pmax=self.Pmax,
            uncertainty=self.uncertainty,
            solute=self.solute,
            comment=self.comment,
        )

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """
        Returns ``True`` if kinetics matches that of another kinetics model.  Must match temperature
        and pressure range of kinetics model, as well as parameters: A, n, Ea, T0. (Shouldn't have pressure
        range if it's Arrhenius.) Otherwise returns ``False``.
        """
        if not isinstance(other_kinetics, ArrheniusChargeTransferBM):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        if (not self.A.equals(other_kinetics.A) or not self.n.equals(other_kinetics.n)
                or not self.E0.equals(other_kinetics.E0) or not self.w0.equals(other_kinetics.w0)
                or not self.alpha.equals(other_kinetics.alpha)
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

def get_w0(actions, rxn):
    """
    calculates the w0 for Blower Masel kinetics by calculating wf (total bond energy of bonds formed)
    and wb (total bond energy of bonds broken) with w0 = (wf+wb)/2
    """
    mol = None
    for r in rxn.reactants:
        m = r.molecule[0]
        if mol:
            mol = mol.merge(m)
        else:
            mol = m.copy(deep=True)
    a_dict = mol.get_all_labeled_atoms()

    recipe = actions

    wb = 0.0
    wf = 0.0
    for act in recipe:

        if act[0] in ['BREAK_BOND','FORM_BOND','CHANGE_BOND']:

            if act[1] == act[3]: # the labels are the same
                atom1 = a_dict[act[1]][0]
                atom2 = a_dict[act[3]][1]
            else:
                atom1 = a_dict[act[1]]
                atom2 = a_dict[act[3]]

        if act[0] == 'BREAK_BOND':
            bd = mol.get_bond(atom1, atom2)
            wb += bd.get_bde()
        elif act[0] == 'FORM_BOND':
            bd = Bond(atom1, atom2, act[2])
            wf += bd.get_bde()
        elif act[0] == 'CHANGE_BOND':
            bd1 = mol.get_bond(atom1, atom2)

            if act[2] + bd1.order == 0.5:
                mol2 = None
                for r in rxn.products:
                    m = r.molecule[0]
                    if mol2:
                        mol2 = mol2.merge(m)
                    else:
                        mol2 = m.copy(deep=True)
                a_dict_mol2 = mol2.get_all_labeled_atoms()
                if act[1] == act[3]: # the labels are the same
                    atom1_mol2 = a_dict_mol2[act[1]][0]
                    atom2_mol2 = a_dict_mol2[act[3]][1]
                else:
                    atom1_mol2 = a_dict_mol2[act[1]]
                    atom2_mol2 = a_dict_mol2[act[3]]

                bd2 = mol2.get_bond(atom1_mol2, atom2_mol2)
            else:
                bd2 = Bond(atom1, atom2, bd1.order + act[2])

            if bd2.order == 0:
                bd2_bde = 0.0
            else:
                bd2_bde = bd2.get_bde()
            bde_diff = bd2_bde - bd1.get_bde()
            if bde_diff > 0:
                wf += abs(bde_diff)
            else:
                wb += abs(bde_diff)
    return (wf + wb) / 2.0

def get_w0s(actions, rxns):
    return [get_w0(actions, rxn) for rxn in rxns]
