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

"""
Contains functionality for directly simulating the master equation
and implementing the SLS master equation reduction method
"""

from diffeqpy import de
from julia import Main
import scipy.sparse as sparse
import numpy as np
import scipy.linalg
import mpmath
import scipy.optimize as opt

import rmgpy.constants as constants
from rmgpy.pdep.reaction import calculate_microcanonical_rate_coefficient
from rmgpy.pdep.me import generate_full_me_matrix
from rmgpy.statmech.translation import IdealGasTranslation

def get_initial_condition(network,x0,indices):
    """
    generates a master equation initial condition
    from total species concentrations x0
    distributes population boltzmann within each well
    """
    e_list = network.e_list
    j_list = network.j_list
    dens_states = network.dens_states

    n_isom = network.n_isom
    n_reac = network.n_reac
    n_prod = network.n_prod
    n_grains = len(e_list)
    n_j = len(j_list)

    # Get equilibrium distributions
    eq_dist = np.zeros_like(dens_states)
    for i in range(n_isom):
        for s in range(n_j):
            eq_dist[i, :, s] = dens_states[i, :, s] * (2 * j_list[s] + 1) * np.exp(-e_list / constants.R / network.T)
        eq_dist[i, :, :] /= sum(eq_dist[i, :, :])

    # Set initial conditions
    p0 = np.zeros(np.sum(indices >= 0.0)+n_reac+n_prod, float)
    for i in range(n_isom):
        for r in range(n_grains):
            for s in range(n_j):
                index = indices[i, r, s]
                if indices[i, r, s] > 0:
                    p0[index] = x0[i] * eq_dist[i, r, s]
    for i in range(n_reac + n_prod):
        p0[-n_reac - n_prod + i] = x0[i + n_isom]

    return p0

def solve_me(M,p0,t):
    f = Main.eval("""
function f(u,M,t)
    return M*u
end""")
    jac = Main.eval("""
function jac(u,M,t)
    return M
end""")
    tspan = (0.0,t)
    fcn = de.ODEFunction(f,jac=jac)
    prob = de.ODEProblem(fcn,p0,tspan,M)
    sol = de.solve(prob,solver=de.CVODE_BDF(),abstol=1e-16,reltol=1e-6)
    return sol

def solve_me_fcns(f,jac,M,p0,t):
    tspan = (0.0,t)
    fcn = de.ODEFunction(f,jac=jac)
    prob = de.ODEProblem(fcn,p0,tspan,M)
    sol = de.solve(prob,solver=de.CVODE_BDF(),abstol=1e-16,reltol=1e-6)
    return sol

