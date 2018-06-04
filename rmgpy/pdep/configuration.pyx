###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

"""
This module contains classes and functions for working with molecular
configurations on a potential energy surface. This includes local minima
(unimolecular and bimolecular) and first-order saddle points (transition
states).
"""

import math
import numpy
import logging
import cython

from libc.math cimport log, exp, sqrt

import rmgpy.constants as constants

from rmgpy.pdep.collision import *
from rmgpy.statmech import *
from rmgpy.statmech.conformer import getDensityOfStatesForst
from rmgpy.transport import TransportData

from rmgpy.species import Species, TransitionState
from rmgpy.reaction import Reaction

################################################################################

cdef class Configuration:
    """
    A representation of a molecular configuration on a potential energy
    surface.
    """
    
    def __init__(self, *species):
        self.species = list(species)
        self.Elist = None
        self.densStates = None
        self.sumStates = None
        self.activeJRotor = False
        self.activeKRotor = False
        
    def __str__(self):
        return ' + '.join([str(spec) for spec in self.species])
        
    def __repr__(self):
        string = 'Configuration('
        string += 'species="{0!r}", '.format(self.species)
        if self.Elist is not None: string += 'Elist={0}, '.format(self.Elist)
        if self.densStates is not None: string += 'densStates={0}, '.format(self.densStates)
        if self.sumStates is not None: string += 'sumStates={0}, '.format(self.sumStates)
        string += 'activeKRotor={0}, '.format(self.activeKRotor)
        string += 'activeJRotor={0}, '.format(self.activeJRotor)
        string += ')'
        return string

    property E0:
        """The ground-state energy of the configuration in J/mol."""
        def __get__(self):
            return sum([float(spec.conformer.E0.value_si) for spec in self.species])
    
    cpdef cleanup(self):
        """
        Delete intermediate arrays used in computing k(T,P) values.
        """
        self.Elist = None
        self.densStates = None
        self.sumStates = None
    
    cpdef bint isUnimolecular(self) except -2:
        """
        Return ``True`` if the configuration represents a unimolecular isomer,
        or ``False`` otherwise.
        """
        return len(self.species) == 1 and isinstance(self.species[0], Species)
    
    cpdef bint isBimolecular(self) except -2:
        """
        Return ``True`` if the configuration represents a bimolecular reactant
        or product channel, or ``False`` otherwise.
        """
        return len(self.species) == 2
    
    cpdef bint isTransitionState(self) except -2:
        """
        Return ``True`` if the configuration represents a transition state,
        or ``False`` otherwise.
        """
        return len(self.species) == 1 and isinstance(self.species[0], TransitionState)
    
    cpdef bint hasStatMech(self) except -2:
        """
        Return ``True`` if all species in the configuration have statistical
        mechanics parameters, or ``False`` otherwise.
        """
        return all([spec.hasStatMech() for spec in self.species])
    
    cpdef bint hasThermo(self) except -2:
        """
        Return ``True`` if all species in the configuration have thermodynamics
        parameters, or ``False`` otherwise.
        """
        return all([spec.hasThermo() for spec in self.species])
    
    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the constant-pressure heat capacity in J/mol*K at the
        specified temperature `T` in K.
        """
        cdef double Cp = 0.0
        for spec in self.species:
            Cp += spec.getHeatCapacity(T)
        return Cp

    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the enthalpy in kJ/mol at the specified temperature `T` in K.
        """
        cdef double H = 0.0
        for spec in self.species:
            H += spec.getEnthalpy(T)
        return H

    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the entropy in J/mol*K at the specified temperature `T` in K.
        """
        cdef double S = 0.0
        for spec in self.species:
            S += spec.getEntropy(T)
        return S

    cpdef double getFreeEnergy(self, double T) except 100000000:
        """
        Return the Gibbs free energy in kJ/mol at the specified temperature
        `T` in K.
        """
        cdef double G = 0.0
        for spec in self.species:
            G += spec.getFreeEnergy(T)
        return G
    
    cpdef double calculateCollisionFrequency(self, double T, double P, dict bathGas) except -1:
        """
        Return the value of the collision frequency in Hz at the given 
        temperature `T` in K and pressure `P` in Pa. If a dictionary `bathGas`
        of bath gas species and corresponding mole fractions is given, the
        collision parameters of the bas gas species will be averaged with those
        of the species before computing the collision frequency. 
        
        Only the Lennard-Jones collision model is currently supported.
        """
        cdef double bathGasSigma, bathGasEpsilon, bathGasMW
        cdef double sigma, epsilon, mu, gasConc, frac, Tred, omega22
        
        assert self.isUnimolecular()
        assert isinstance(self.species[0].getTransportData(), TransportData)
        for spec, frac in bathGas.items():
            assert isinstance(spec.getTransportData(), TransportData)
        
        bathGasSigma = 0.0; bathGasEpsilon = 1.0; bathGasMW = 0.0
        for spec, frac in bathGas.iteritems():
            bathGasSigma += spec.getTransportData().sigma.value_si * frac
            bathGasEpsilon *= spec.getTransportData().epsilon.value_si ** frac
            bathGasMW += spec._molecularWeight.value_si * frac
        
        sigma = 0.5 * (self.species[0].getTransportData().sigma.value_si + bathGasSigma)
        epsilon = sqrt((self.species[0].getTransportData().epsilon.value_si * bathGasEpsilon))
        mu = 1.0 / (1.0/self.species[0]._molecularWeight.value_si + 1.0/bathGasMW)
        gasConc = P / constants.kB / T
        
        # Evaluate configuration integral
        Tred = constants.R * T / epsilon
        omega22 = 1.16145 * Tred**(-0.14874) + 0.52487 * exp(-0.77320 * Tred) + 2.16178 * exp(-2.43787 * Tred)
    
        # Evaluate collision frequency
        return omega22 * sqrt(8 * constants.kB * T / constants.pi / mu) * constants.pi * sigma*sigma * gasConc
        
    cpdef numpy.ndarray generateCollisionMatrix(self, double T, numpy.ndarray densStates, numpy.ndarray Elist, numpy.ndarray Jlist=None):
        """
        Return the collisional energy transfer probabilities matrix for the
        configuration at the given temperature `T` in K using the given
        energies `Elist` in kJ/mol and total angular momentum quantum numbers
        `Jlist`. The density of states of the configuration `densStates` in
        mol/kJ is also required.
        """
        assert self.isUnimolecular()
        assert self.species[0].energyTransferModel is not None
        return self.species[0].energyTransferModel.generateCollisionMatrix(T, densStates, Elist, Jlist)
    
    cpdef calculateDensityOfStates(self, numpy.ndarray Elist, bint activeJRotor=True, bint activeKRotor=True, bint rmgmode=False):
        """
        Calculate the density (and sum) of states for the configuration at the
        given energies above the ground state `Elist` in J/mol. The 
        `activeJRotor` and `activeKRotor` flags control whether the J-rotor
        and/or K-rotor are treated as active (and therefore included in the
        density and sum of states). The computed density and sum of states
        arrays are stored on the object for future use.
        """
        cdef list modes
        cdef int i
        
        self.Elist = Elist
        self.activeJRotor = activeJRotor
        self.activeKRotor = activeKRotor
        
        # Get the active rovibrational modes for each species in the configuration
        modes = []
        for i, species in enumerate(self.species):
            modes.extend(species.conformer.getActiveModes(activeKRotor=activeKRotor, activeJRotor=activeJRotor))
        
        if rmgmode:
            # Include an arbitrary active rigid rotor if needed
            # The moments of inertia cancel in all subsequent calculations
            for mode in modes:
                if isinstance(mode, (LinearRotor,NonlinearRotor)):
                    break
            else:
                linear = False
                for species in self.species:
                    for molecule in species.molecule:
                        if molecule.isLinear(): 
                            linear = True
                            break
                if linear:
                    modes.insert(0, LinearRotor(inertia=(1.0,"amu*angstrom^2"), symmetry=1))
                else:
                    modes.insert(0, NonlinearRotor(inertia=([1.0,1.0,1.0],"amu*angstrom^2"), symmetry=1))
        
        if len(modes) == 0:
            self.densStates = None
            self.sumStates = None
        else:        
            # If the configuration is bimolecular, also include the relative
            # translational motion of the two molecules
            if self.isBimolecular():
                mass = []
                for species in self.species:
                    for mode in species.conformer.modes:
                        if isinstance(mode, IdealGasTranslation):
                            mass.append(mode.mass.value_si)
                            break
                    else:
                        if species.molecularWeight is not None:
                            mass.append(species.molecularWeight.value_si)
                        else:
                            m = 0
                            for atom in species.molecule[0].atoms:
                                m += atom.element.mass
                            mass.append(m * constants.amu * 1000)
                assert len(mass) == 2
                mu = 1.0/(1.0/mass[0] + 1.0/mass[1])
                modes.insert(0, IdealGasTranslation(mass=(mu/constants.amu,"amu")))
                
            if rmgmode:
                # Compute the density of states by direct count
                # This is currently faster than the method of steepest descents,
                # but requires classical hindered rotors
                densStates = None
                for mode in modes:
                    if not isinstance(mode,HarmonicOscillator):
                        densStates = mode.getDensityOfStates(Elist, densStates)
                    # Fix a numerical artifact that occurs when two modes have
                    # density of states expressions that are zero at the
                    # ground state
                    # Convoluting these modes gives a zero at the first excited
                    # energy grain as well, which is unphysical
                    # Instead, fill in an approximate value by extrapolation
                    # This should only occur in systems with IdealGasTranslation
                    # and NonlinearRotor modes
                    if densStates[1] == 0:
                        densStates[1] = densStates[2] * densStates[2] / densStates[3]
                for mode in modes:
                    if isinstance(mode,HarmonicOscillator):
                        densStates = mode.getDensityOfStates(Elist, densStates)
                self.densStates = densStates
                for spec in self.species:
                    self.densStates *= spec.conformer.spinMultiplicity * spec.conformer.opticalIsomers
                
            else:
                # Since the evaluation of quantum hindered rotors is slow, it is
                # actually faster (and probably negligibly less accurate) to use
                # interpolation in the steepest descents algorithm
                import scipy.interpolate
                
                logTdata = numpy.linspace(log(10.), log(10000.), 250.)
                Tdata = numpy.exp(logTdata)
                Qdata = numpy.ones_like(Tdata)
                for i in range(Tdata.shape[0]):
                    T = Tdata[i]
                    for mode in modes:
                        Qdata[i] = Qdata[i] * mode.getPartitionFunction(T)
                logQ = scipy.interpolate.InterpolatedUnivariateSpline(Tdata, numpy.log(Qdata))
                #logQ = LinearInterpolator(Tdata, numpy.log(Qdata))          
    
                self.densStates, self.sumStates = getDensityOfStatesForst(Elist, logQ)
            
                for spec in self.species:
                    self.densStates *= spec.conformer.spinMultiplicity * spec.conformer.opticalIsomers
                    self.sumStates *= spec.conformer.spinMultiplicity * spec.conformer.opticalIsomers
        if self.densStates is None:
            raise ValueError("Species {} has no active modes".format(species.label))
            
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def mapDensityOfStates(self, numpy.ndarray[numpy.float64_t,ndim=1] Elist, numpy.ndarray[numpy.int_t,ndim=1] Jlist=None):
        """
        Return a mapping of the density of states for the configuration to the
        given energies `Elist` in J/mol and, if the J-rotor is not active, the
        total angular momentum quantum numbers `Jlist`.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=2] densStates
        cdef double E0, B1, B2, E, dJ
        cdef int r0, r, s, t, Ngrains, NJ, J, J1, J2
        cdef list Blist
        
        import scipy.interpolate
        for r in range(self.Elist.shape[0]):
            if self.densStates[r] > 0:
                break
        f = scipy.interpolate.InterpolatedUnivariateSpline(self.Elist[r:], numpy.log(self.densStates[r:]))
        
        E0 = self.E0
        Ngrains = Elist.shape[0]
        dE = Elist[1] - Elist[0]
        dE0 = self.Elist[1] - self.Elist[0]
        
        if self.activeJRotor:
            densStates = numpy.zeros((Ngrains,1))
            for r0 in range(Ngrains):
                if Elist[r0] >= E0: break
            for r in range(r0, Ngrains):
                densStates[r,0] = f(Elist[r] - E0)
            densStates[r0:,0] = numpy.exp(densStates[r0:,0])
        else:
            assert Jlist is not None
            NJ = Jlist.shape[0]
            dJ = Jlist[1] - Jlist[0]
            densStates = numpy.zeros((Ngrains, NJ))
            
            Blist = []
            for spec in self.species:
                Jrotor, Krotor = spec.conformer.getSymmetricTopRotors()
                Blist.append(float(Jrotor.rotationalConstant.value_si))
        
            for r0 in range(Ngrains):
                if Elist[r0] >= E0: break
                
            if len(Blist) == 1:
                B1 = Blist[0] * 11.962        # cm^-1 to J/mol
                for r in range(r0, Ngrains):
                    for s in range(NJ):
                        J1 = Jlist[s]
                        E = Elist[r] - E0 - B1 * J1 * (J1 + 1)
                        if E < 0: break
                        densStates[r,s] = (2 * J1 + 1) * exp(f(E)) * dJ
                        
            elif len(Blist) == 2:
                B1 = Blist[0] * 11.962; B2 = Blist[1] * 11.962      # cm^-1 to J/mol
                for r in range(r0, Ngrains):
                    for s in range(NJ):
                        J = Jlist[s]
                        for t in range(s+1):
                            J1 = Jlist[t]; J2 = J - J1
                            E = Elist[r] - E0 - B1 * J1 * (J1 + 1) - B2 * J2 * (J2 + 1)
                            if E > 0:
                                densStates[r,s] += (2 * J1 + 1) * (2 * J2 + 2) * exp(f(E)) * dJ * dJ
        
        return densStates * dE / dE0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def mapSumOfStates(self, numpy.ndarray[numpy.float64_t,ndim=1] Elist, numpy.ndarray[numpy.int_t,ndim=1] Jlist=None):
        """
        Return a mapping of the density of states for the configuration to the
        given energies `Elist` in J/mol and, if the J-rotor is not active, the
        total angular momentum quantum numbers `Jlist`.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=2] sumStates
        cdef double E0, B1, B2, dJ
        cdef int r0, r, s, Ngrains, NJ, J1, J2
        
        import scipy.interpolate
        for r in range(self.Elist.shape[0]):
            if self.sumStates[r] > 0:
                break
        f = scipy.interpolate.InterpolatedUnivariateSpline(self.Elist[r:], numpy.log(self.sumStates[r:]))
        
        E0 = self.E0
        Ngrains = len(Elist)
        
        if self.activeJRotor:
            sumStates = numpy.zeros((Ngrains,1))
            for r0 in range(Ngrains):
                if Elist[r0] >= E0: break
            for r in range(r0, Ngrains):
                sumStates[r,0] = f(Elist[r] - E0)
            sumStates[r0:,0] = numpy.exp(sumStates[r0:,0])
        else:
            assert Jlist is not None
            NJ = len(Jlist)
            dJ = Jlist[1] - Jlist[0]
            sumStates = numpy.zeros((Ngrains, NJ))
            
            Blist = []
            for spec in self.species:
                Jrotor, Krotor = spec.conformer.getSymmetricTopRotors()
                Blist.append(float(Jrotor.rotationalConstant.value_si))
        
            for r0 in range(Ngrains):
                if Elist[r0] >= E0: break
                
            if len(Blist) == 1:
                B1 = Blist[0] * 11.962        # cm^-1 to J/mol
                for r in range(r0, Ngrains):
                    for s in range(NJ):
                        J1 = Jlist[s]
                        E = Elist[r] - E0 - B1 * J1 * (J1 + 1)
                        if E < 0: break
                        sumStates[r,s] = (2 * J1 + 1) * exp(f(E)) * dJ
                        
            elif len(Blist) == 2:
                B1 = Blist[0] * 11.962; B2 = Blist[1] * 11.962      # cm^-1 to J/mol
                for r in range(r0, Ngrains):
                    for s in range(NJ):
                        J = Jlist[s]
                        for t in range(s+1):
                            J1 = Jlist[t]; J2 = J - J1
                            E = Elist[r] - E0 - B1 * J1 * (J1 + 1) - B2 * J2 * (J2 + 1)
                            if E > 0:
                                sumStates[r,s] += (2 * J1 + 1) * (2 * J2 + 1) * exp(f(E)) * dJ * dJ
        
        return sumStates
