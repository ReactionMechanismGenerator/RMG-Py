#!/usr/bin/env python3

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
Contains classes for building RMS simulations
Note since rms currently can't be imported inside nosetests
This most of this can't be tested using nosetests
"""
import numpy as np
import sys
import logging
import itertools
try:
    from pyrms import rms
    from diffeqpy import de
except:
    pass

from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.thermo.nasa import NASAPolynomial, NASA
from rmgpy.thermo.wilhoit import Wilhoit
from rmgpy.thermo.thermodata import ThermoData
from rmgpy.kinetics.arrhenius import Arrhenius, ArrheniusEP, ArrheniusBM, PDepArrhenius, MultiArrhenius, MultiPDepArrhenius
from rmgpy.kinetics.kineticsdata import KineticsData
from rmgpy.kinetics.falloff import Troe, ThirdBody, Lindemann
from rmgpy.kinetics.chebyshev import Chebyshev
from rmgpy.data.solvation import SolventData
from rmgpy.kinetics.surface import StickingCoefficient
from rmgpy.solver.termination import TerminationTime, TerminationConversion, TerminationRateRatio
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.depository import DepositoryReaction

def to_rms(obj, species_list=None, rms_species_list=None):
    """
    Generate corresponding rms object
    """
    if isinstance(obj, Arrhenius):
        if obj._T0.value_si != 1:
            A = obj._A.value_si / (obj._T0.value_si) ** obj._n.value_si
        else:
            A = obj._A.value_si
        n = obj._n.value_si
        Ea = obj._Ea.value_si
        return rms.Arrhenius(A, n, Ea, rms.EmptyRateUncertainty())
    elif isinstance(obj, PDepArrhenius):
        Ps = obj._pressures.value_si
        arrs = [to_rms(arr) for arr in obj.arrhenius]
        return rms.PDepArrhenius(Ps, arrs)
    elif isinstance(obj, MultiArrhenius):
        arrs = [to_rms(arr) for arr in obj.arrhenius]
        return rms.MultiArrhenius(arrs)
    elif isinstance(obj, MultiPDepArrhenius):
        parrs = [to_rms(parr) for parr in obj.arrhenius]
        return rms.MultiPdepArrhenius(parrs)
    elif isinstance(obj, Chebyshev):
        Tmin = obj.Tmin.value_si
        Tmax = obj.Tmax.value_si
        Pmin = obj.Pmin.value_si
        Pmax = obj.Pmax.value_si
        coeffs = obj.coeffs.value_si.tolist()
        return rms.Chebyshev(coeffs, Tmin, Tmax, Pmin, Pmax)
    elif isinstance(obj, ThirdBody):
        arr = to_rms(obj.arrheniusLow)
        efficiencies = {spc.label: float(val) for spc, val in obj.efficiencies.items() if val != 1}
        return rms.ThirdBody(arr, nameefficiencies=efficiencies)
    elif isinstance(obj, Lindemann):
        arrlow = to_rms(obj.arrheniusLow)
        arrhigh = to_rms(obj.arrheniusHigh)
        efficiencies = {spc.label: float(val) for spc, val in obj.efficiencies.items() if val != 1}
        return rms.Lindemann(arrhigh, arrlow, nameefficiencies=efficiencies)
    elif isinstance(obj, Troe):
        arrlow = to_rms(obj.arrheniusLow)
        arrhigh = to_rms(obj.arrheniusHigh)
        efficiencies = {spc.label: float(val) for spc, val in obj.efficiencies.items() if val != 1}
        alpha = obj.alpha
        T1 = obj._T1.value_si if obj._T1 is not None else 0.0
        T2 = obj._T2.value_si if obj._T2 is not None else 0.0
        T3 = obj._T3.value_si if obj._T3 is not None else 0.0
        return rms.Troe(arrhigh, arrlow, alpha, T3, T1, T2, nameefficiencies=efficiencies)
    elif isinstance(obj, StickingCoefficient):
        if obj._T0.value_si != 1:
            A = obj._A.value_si / (obj._T0.value_si) ** obj._n.value_si
        else:
            A = obj._A.value_si
        n = obj._n.value_si
        Ea = obj._Ea.value_si
        return rms.StickingCoefficient(A, n, Ea)
    elif isinstance(obj, NASAPolynomial):
        return rms.NASApolynomial(obj.coeffs, obj.Tmin.value_si, obj.Tmax.value_si)
    elif isinstance(obj, NASA):
        return rms.NASA([to_rms(poly) for poly in obj.polynomials], rms.EmptyThermoUncertainty())
    elif isinstance(obj, Species):
        atomnums = dict()
        for atm in obj.molecule[0].atoms:
            if atomnums.get(atm.element.symbol):
                atomnums[atm.element.symbol] += 1
            else:
                atomnums[atm.element.symbol] = 1
        bondnum = len(obj.molecule[0].get_all_edges())
        rad = rms.getspeciesradius(atomnums, bondnum)
        diff = rms.StokesDiffusivity(rad)
        th = obj.get_thermo_data()
        thermo = to_rms(th)
        return rms.Species(obj.label, obj.index, "", "", "", thermo, atomnums, bondnum, diff, rad, obj.molecule[0].multiplicity-1, obj.molecular_weight.value_si)
    elif isinstance(obj, Reaction):
        reactantinds = [species_list.index(spc) for spc in obj.reactants]
        productinds = [species_list.index(spc) for spc in obj.products]
        reactants = [rms_species_list[i] for i in reactantinds]
        products = [rms_species_list[i] for i in productinds]
        kinetics = to_rms(obj.kinetics)
        radchange = sum([spc.molecule[0].multiplicity-1 for spc in obj.products]) - sum([spc.molecule[0].multiplicity-1 for spc in obj.reactants])
        electronchange = 0 #for now
        return rms.ElementaryReaction(obj.index, reactants, reactantinds, products, productinds, kinetics, electronchange, radchange, obj.reversible, [])
    elif isinstance(obj, SolventData):
        return rms.Solvent("solvent", rms.RiedelViscosity(float(obj.A), float(obj.B), float(obj.C), float(obj.D), float(obj.E)))
    elif isinstance(obj, TerminationTime):
        return rms.TerminationTime(obj.time.value_si)
    elif isinstance(obj, TerminationConversion):
        return rms.TerminationConversion(to_rms(obj.species), obj.conversion)
    elif isinstance(obj, TerminationRateRatio):
        return rms.TerminationRateRatio(obj.ratio)
    elif isinstance(obj, tuple): #Handle TerminationConversion when the obj doesn't have thermo yet
        return obj
    else:
        errortype = type(obj)
        raise ValueError(f"Couldn't convert object of type {errortype} to an RMS object")
