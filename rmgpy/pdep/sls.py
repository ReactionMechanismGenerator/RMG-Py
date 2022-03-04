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
from rmgpy.pdep.me import generate_full_me_matrix, states_to_configurations
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
    p0 = np.zeros(np.sum(indices >= 0.0)+n_reac, float)
    for i in range(n_isom):
        for r in range(n_grains):
            for s in range(n_j):
                index = indices[i, r, s]
                if indices[i, r, s] > 0:
                    p0[index] = x0[i] * eq_dist[i, r, s]
    for i in range(n_reac):
        p0[-n_reac + i] = x0[i + n_isom]

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

def get_rate_coefficients_SLS(network,T,P,method="mexp",neglect_high_energy_collisions=False,high_energy_rate_tol=0.01):
    """
    Generate phenomenological rate coefficients using the simulation least squares method
    """
    if network.T != T or network.P != P:
        network.set_conditions(T,P)
        network.calculate_equilibrium_ratios()

    M, indices = generate_full_me_matrix(network, products=False,
                                                neglect_high_energy_collisions=neglect_high_energy_collisions,
                                                 high_energy_rate_tol=high_energy_rate_tol)

    if method == "eigen":
        eigvals,eigvecs = np.linalg.eig(M)
        omega = np.real(eigvals)
    else:
        omega = np.real(scipy.linalg.eigvals(M))

    n_chem = network.n_isom + network.n_reac

    # We can't assume that eigh returns them in sorted order
    ind = omega.argsort()

    #Find timescale
    n_cse = 0
    for i in range(n_chem):
        if abs(omega[ind[-n_chem - 1]] / omega[ind[-1 - i]]) > 3.0:
            n_cse += 1

    slowest_energy = omega[ind[-n_chem-1]]
    fastest_reaction = omega[ind[-n_chem]]
    #tau = np.sqrt(1.0/slowest_energy*1.0/fastest_reaction)
    tau = np.abs(1.0/fastest_reaction)

    if method == "ode":
        f = Main.eval("""
function f(u,M,t)
    return M*u
end""")
        jac = Main.eval("""
function jac(u,M,t)
    return M
end""")
        def run_single_source(network,channel):
            i = (network.isomers+network.reactants+network.products).index(channel)
            x0 = np.zeros(len(network.isomers+network.reactants+network.products))
            x0[i] = 1.0
            p0 = get_initial_condition(network,x0,indices)
            sol = solve_me_fcns(f,jac,M,p0,tau)
            out =  sol.u[-1]
            xsout = states_to_configurations(network,indices,out)
            ddt = states_to_configurations(network,indices,np.dot(M,out))
            return xsout,ddt

        def run_equilibrium(network,channel):
            i = (network.isomers+network.reactants+network.products).index(channel)
            x0 = network.eq_ratios/np.sum(network.eq_ratios)
            x0[i] = 0.0
            p0 = get_initial_condition(network,x0,indices)
            sol = solve_me_fcns(f,jac,M,p0,tau)
            out =  sol.u[-1]
            xsout = states_to_configurations(network,indices,out)
            ddt = states_to_configurations(network,indices,np.dot(M,out))
            return xsout,ddt

    elif method == "mexp":
        etM = scipy.linalg.expm(M*tau)
        def run_single_source(network,channel):
            i = (network.isomers+network.reactants+network.products).index(channel)
            x0 = np.zeros(len(network.isomers+network.reactants+network.products))
            x0[i] = 1.0
            p0 = get_initial_condition(network,x0,indices)
            out = etM @ p0
            xsout = states_to_configurations(network,indices,out)
            ddt = states_to_configurations(network,indices,np.dot(M,out))
            return xsout,ddt

        def run_equilibrium(network,channel):
            i = (network.isomers+network.reactants+network.products).index(channel)
            x0 = network.eq_ratios/np.sum(network.eq_ratios)
            x0[i] = 0.0
            p0 = get_initial_condition(network,x0,indices)
            out = etM @ p0
            xsout = states_to_configurations(network,indices,out)
            ddt = states_to_configurations(network,indices,np.dot(M,out))
            return xsout,ddt

    elif method == "eigen":
        def run_single_source(network,channel):
            i = (network.isomers+network.reactants+network.products).index(channel)
            x0 = np.zeros(len(network.isomers+network.reactants+network.products))
            x0[i] = 1.0
            p0 = get_initial_condition(network,x0,indices)
            c0 = np.linalg.solve(eigvecs,p0)
            cf = np.exp(eigvals*tau) * c0
            out = eigvecs @ cf
            xsout = states_to_configurations(network,indices,out)
            ddt = states_to_configurations(network,indices,np.dot(M,out))
            return xsout,ddt

        def run_equilibrium(network,channel):
            i = (network.isomers+network.reactants+network.products).index(channel)
            x0 = network.eq_ratios/np.sum(network.eq_ratios)
            x0[i] = 0.0
            p0 = get_initial_condition(network,x0,indices)
            c0 = np.linalg.solve(eigvecs,p0)
            cf = np.exp(eigvals*tau) * c0
            out = eigvecs @ cf
            xsout = states_to_configurations(network,indices,out)
            ddt = states_to_configurations(network,indices,np.dot(M,out))
            return xsout,ddt
    else:
        raise ValueError

    isomers = network.isomers[:]
    reactants = network.reactants[:]
    products = network.products[:]
    xsorder = isomers+reactants+products
    n_isomreac = len(isomers+reactants)
    xssource = []
    dxdtssource = []
    xseq = []
    dxdtseq = []

    #Single domainant source simulations
    for i,isomer in enumerate(isomers):
        xsout,dxdtout = run_single_source(network,isomer)
        for j in range(len(xsout)):
            if j >= n_isomreac:
                xsout[j] = 0.0
        xssource.append(xsout)
        dxdtssource.append(dxdtout)
    for i,prod in enumerate(reactants):
        xsout,dxdtout = run_single_source(network,prod)
        for j in range(len(xsout)):
            if j >= n_isomreac:
                xsout[j] = 0.0
        xssource.append(xsout)
        dxdtssource.append(dxdtout)

    #Equillibrium single loss simulations
    for i,isomer in enumerate(isomers):
        xsout,dxdtout = run_equilibrium(network,isomer)
        for j in range(len(xsout)):
            if j >= n_isomreac:
                xsout[j] = 0.0
        xseq.append(xsout)
        dxdtseq.append(dxdtout)
    for i,prod in enumerate(reactants+products):
        xsout,dxdtout = run_equilibrium(network,prod)
        for j in range(len(xsout)):
            if j >= n_isomreac:
                xsout[j] = 0.0
        xseq.append(xsout)
        dxdtseq.append(dxdtout)

    keqs = np.outer(network.eq_ratios,1.0/network.eq_ratios)

    #Generate initial guess
    kmat = np.zeros((len(xsorder),len(xsorder)))
    ks = np.zeros((len(xsorder)**2-len(xsorder)-len(products)**2+len(products))//2)

    for i,x in enumerate(xssource):
        dxdt = dxdtssource[i]
        kmat[:,i] = dxdt/x[i]
        kmat[i,i] = 0.0

    for i in range(kmat.shape[0]):
        for j in range(kmat.shape[1]):
            if i != j and kmat[i,j] == 0.0 and kmat[j,i] != 0.0:
                kmat[i,j] = kmat[j,i]/keqs[j,i]

    ks = ravel_kmat(kmat,n_isomreac)

    ks = np.abs(ks)

    #Do optimization
    def fcn(ks):
        kmat = unravel_ks(np.exp(ks),keqs,n_isomreac)
        s = np.zeros(len(xssource)*(len(isomers)+len(reactants)+len(products)-1)+len(xseq))
        ind = 0
        for i,x in enumerate(xssource):
            flux = calcfluxes(kmat,x)
            dxdt = dxdtssource[i]
            res = [(dxdt[j]-flux[j])/dxdt[j] for j in range(len(flux)) if i != j]
            s[ind:ind+len(res)] = res[:]
            ind = ind + len(res)
        for i,x in enumerate(xseq):
            flux = calcfluxes(kmat,x)
            dxdt = dxdtseq[i]
            s[ind] = (dxdt[i] - flux[i])/dxdt[i]
            ind += 1
        return s

    out = opt.least_squares(fcn,np.log(ks),ftol=1e-15,xtol=1e-15,gtol=1e-15)

    kmat = unravel_ks(np.exp(out.x),keqs,len(network.isomers+network.reactants))

    u = get_uncertainties(kmat,xssource,dxdtssource)

    return kmat,u

def get_uncertainties(kmat,xssource,dxdtssource):
    """
    approximate the independent factor uncertainty in individual
    rate coefficients
    """
    u = np.zeros(kmat.shape)
    for i,x in enumerate(xssource):
        flux = calcfluxes(kmat,x)
        dxdt = dxdtssource[i]
        for j in range(len(xssource)):
            u[j,i] = np.exp(np.abs(np.log(dxdt[j]/flux[j])))
    return u

def ravel_kmat(kmat,n_isomreac):
    """
    reduce rate coefficient matrix to a vector
    """
    ks = []
    ind = 0
    for i in range(kmat.shape[0]):
        if i < n_isomreac:
            indval = i
        else:
            indval = n_isomreac
        ks.extend(kmat[i,:indval].tolist())
        ind += indval
    return np.array(ks)

def ravel_kmat_mult(kmat,n_isomreac,keqs,bits):
    """
    reduce rate coefficient matrix to a vector in different ways
    """
    ks = []
    indk = 0
    for i in range(kmat.shape[0]):
        for j in range(i):
            if i > j and i < n_isomreac:
                if bits[indk]:
                    ks.append(kmat[i,j])
                else:
                    ks.append(kmat[j,i]/keqs[j,i])
                indk += 1
    return np.array(ks)

def unravel_ks(ks,keqs,n_isomreac):
    """
    convert SLS rate coefficient vector to rate coefficient matrix
    """
    kmat = np.zeros((keqs.shape[0],keqs.shape[0]))
    ind = 0
    for i in range(keqs.shape[0]):
        if i < n_isomreac:
            indval = i
        else:
            indval = n_isomreac
        kmat[i,:indval] = ks[ind:ind+indval]
        ind += i
    for i in range(keqs.shape[0]):
        for j in range(keqs.shape[0]):
            if i != j and kmat[j,i] == 0.0:
                kmat[j,i] = kmat[i,j]/keqs[i,j]
    return kmat

def get_names(channel):
    return [x.label for x in channel.species]

def calcfluxes(kmat,xs):
    """
    calculate fluxes from phenomenological rate coefficients
    """
    fluxes = np.zeros(len(xs))
    for i in range(kmat.shape[0]):
        for j in range(i):
            if i != j:
                flux = kmat[j,i]*xs[i]-kmat[i,j]*xs[j]
                fluxes[j] += flux
                fluxes[i] -= flux
    return fluxes

def apply_simulation_least_squares_method(network, method='mexp', neglect_high_energy_collisions=False, high_energy_rate_tol=0.01):
    return get_rate_coefficients_SLS(network, network.T, network.P, method=method,
                                     neglect_high_energy_collisions=neglect_high_energy_collisions, high_energy_rate_tol=high_energy_rate_tol)
