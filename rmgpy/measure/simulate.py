#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   MEASURE - Master Equation Automatic Solver for Unimolecular REactions
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
Contains functions for directly solving the full or reduced master equation
matrix to generate concentration profiles and population distributions. This
information is particularly useful in comparing and evaluating the methods for
reducing the master equation.
"""

import numpy
import scipy.integrate
#from pydas import DASSL

import rmgpy.constants as constants

################################################################################

def residual(t, y, K):
    return numpy.dot(K, y)

def jacobian(t, y, K):
    return K

################################################################################

def solveFullME(T, P, Elist, tlist, x0, M, indices, densStates, Nisom, Nreac, Nprod):
    """
    Directly solve the full master equation using a stiff ODE solver. Pass the
    reaction `network` to solve, the temperature `T` in K and pressure `P` in
    Pa to solve at, the energies `Elist` in J/mol to use, the output time
    points `tlist` in s, the initial total populations `x0`, the full master
    equation matrix `M`, the accounting matrix `indices` relating isomer and
    energy grain indices to indices of the master equation matrix, and the
    densities of states `densStates` in mol/J of each isomer.
    Returns the times in s, population distributions for each isomer, and total
    population profiles for each configuration.
    """

    Ngrains = len(Elist)
    Ntime = len(tlist)

    # Get equilibrium distributions
    eqDist = numpy.zeros_like(densStates)
    for i in range(Nisom):
        eqDist[i,:] = densStates[i,:] * numpy.exp(-Elist / constants.R / T)
        eqDist[i,:] /= sum(eqDist[i,:])

    # Set initial conditions
    p0 = numpy.zeros([M.shape[0]], float)
    for i in range(Nisom):
        for r in range(Ngrains):
            if indices[r,i] > 0:
                p0[indices[r,i]] = x0[i] * eqDist[i,r]
    for i in range(Nreac+Nprod):
        p0[-Nreac-Nprod + i] = x0[i+Nisom]


#    # Set up ODEs
#    me = MasterEquation(M)
#    me.initialize(t0=0, y0=p0, atol=1e-16, rtol=1e-8)
#
#    # Generate solution
#    t = numpy.zeros([Ntime], float)
#    p = numpy.zeros([Ntime, Nisom, Ngrains], float)
#    x = numpy.zeros([Ntime, Nisom+Nreac+Nprod], float)
#    for s in range(Ntime):
#        me.advance(tlist[s])
#        print me.t
#        t[s] = me.t
#        for r in range(Ngrains):
#            for i in range(0, Nisom):
#                if indices[r,i] > 0:
#                    p[s,i,r] += me.y[indices[r,i]]
#                    x[s,i] += me.y[indices[r,i]]
#        for n in range(Nisom, Nisom+Nreac+Nprod):
#            x[s,n] = me.y[-(Nisom+Nreac+Nprod)+n]

    # Set up ODEs
    ode = scipy.integrate.ode(residual, jacobian).set_integrator('vode', method='bdf', with_jacobian=True, atol=1e-16, rtol=1e-8)
    ode.set_initial_value(p0, 0.0).set_f_params(M).set_jac_params(M)

    # Generate solution
    t = numpy.zeros([Ntime], float)
    p = numpy.zeros([Ntime, Nisom, Ngrains], float)
    x = numpy.zeros([Ntime, Nisom+Nreac+Nprod], float)
    for s in range(Ntime):
        ode.integrate(tlist[s])
        t[s] = ode.t
        for r in range(Ngrains):
            for i in range(0, Nisom):
                if indices[r,i] > 0:
                    p[s,i,r] += ode.y[indices[r,i]]
                    x[s,i] += ode.y[indices[r,i]]
        for n in range(Nisom, Nisom+Nreac+Nprod):
            x[s,n] = ode.y[-(Nisom+Nreac+Nprod)+n]

    #import pylab
    #pylab.loglog(t,x)
    #pylab.show()

    return t, p, x

################################################################################

def solveReducedME(T, P, Elist, tlist, x0, K, p0, Nisom, Nreac, Nprod):
    """
    Directly solve a reduced master equation (set of phenomenological rate
    coefficients) using a stiff ODE solver. Pass the reaction `network` to
    solve, the temperature `T` in K and pressure `P` in Pa to solve at, the
    energies `Elist` in J/mol to use, the output time points `tlist` in s,
    the initial total populations `x0`, the matrix of phenomenological rate
    coefficients `K`, and the pseudo-steady population distributions `p0`.
    Returns the times in s, approximate population distributions for each
    isomer, and total population profiles for each configuration.
    """
    
    Ngrains = len(Elist)
    Ntime = len(tlist)

    # Set up ODEs
    ode = scipy.integrate.ode(residual, jacobian)
    ode.set_integrator('vode', method='bdf', with_jacobian=True, atol=1e-16, rtol=1e-8)
    ode.set_initial_value(x0, 0.0).set_f_params(K).set_jac_params(K)

    # Generate solution
    t = numpy.zeros([len(tlist)], float)
    x = numpy.zeros([len(tlist), len(x0)], float)
    for n in range(Ntime):
        ode.integrate(tlist[n])
        t[n] = ode.t
        x[n,:] = ode.y

    # Construct p from x
    p = numpy.zeros((len(t),Nisom,Ngrains), numpy.float64)
    for s in range(Ntime):
        for i in range(Nisom):
            for n in range(Nisom):
                p[s,i,:] += x[s,n] * p0[:,i,n]
            for n in range(Nisom, Nisom+Nreac):
                p[s,i,:] += x[s,n] * p0[:,i,n] * (0.21 * P / constants.R / T)

    return t, p, x

#################################################################################
#
#class MasterEquation(DASSL):
#    """
#    """
#
#    def __init__(self, M):
#        DASSL.__init__(self)
#        self.M = M
#
#    def residual(self, t, y, dydt):
#        return numpy.dot(self.M, y) - dydt, 0
#
#    def jacobian(self, t, y, dydt, cj):
#        return self.M - cj * numpy.identity(y.shape[0], numpy.float64)
