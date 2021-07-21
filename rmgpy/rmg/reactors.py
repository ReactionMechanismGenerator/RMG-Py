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


class PhaseSystem:
    """
    Class for tracking and managing all the phases and interfaces of species/reactions
    needed to run a simulation
    will typically have a core and edge PhaseSystem where the edge PhaseSystem
    contains all species/reactions in the core PhaseSystem
    """
    def __init__(self, phases, interfaces):
        self.phases = phases
        self.interfaces = interfaces
        self.species_dict = dict()
        if 'nose' in sys.modules.keys(): #pyrms and nosetests aren't compatible
            self.in_nose = True
        else:
            self.in_nose = False

    def remove_species(self, spc):
        """
        If the species is found
        removes that species and any reactions it is involved in
        from the PhaseSystem and its Phases/Interfaces
        """
        for key, phase in self.phases.items():
            rms_spc = phase.remove_species(spc.label)

        if rms_spc:
            for key, interface in self.interfaces.items():
                interface.remove_species(rms_spc)

            del self.species_dict[spc.label]

    def add_reaction(self, rxn, species_list):
        """
        adds reaction to the appropriate phase/interface
        within PhaseSystem
        """
        spclabels = [r.label for r in rxn.reactants+rxn.products]
        phaseinv = []
        for plabel, phase in self.phases.items():
            out = [label for label in spclabels if label in phase.names]
            if len(out) > 0:
                phaseinv.append(phase)

        if len(phaseinv) == 1:
            phaseinv[0].add_reaction(rxn, species_list)
        else:
            phases = set(phaseinv)
            for interface in self.interfaces:
                if interface.phaseset == phases:
                    interface.add_reaction(rxn, species_list)
                    break

    def pass_species(self, label, phasesys):
        """
        Adds a species from self to the input phase system phasesys
        (usually adding an edge species to the core PhaseSystem)
        also adds any reactions whose participating species are now all
        in the core to the core PhaseSystem
        reorders edge species ordering to match that of the core
        """
        spc = None
        for plabel, phase in self.phases.items():
            phase_label = plabel
            spc = phase.get_species(label)
            if spc is not None:
                break

        assert spc.name not in phasesys.phases[phase_label].names, spc.name

        phasesys.phases[phase_label].species.append(spc)
        phasesys.phases[phase_label].names.append(label)
        phasesys.species_dict[label] = self.species_dict[label]
        self.phases[phase_label].species.remove(spc)
        self.phases[phase_label].names.remove(label)
        self.phases[phase_label].species.insert(len(phasesys.phases[phase_label].species)-1, spc)
        self.phases[phase_label].names.insert(len(phasesys.phases[phase_label].species)-1, label)

        rxnlist = []
        for i, rxn in enumerate(self.phases[phase_label].reactions):
            if (spc in rxn.reactants or spc in rxn.products) and all([spec in phasesys.phases[phase_label].species for spec in rxn.reactants]) and all([spec in phasesys.phases[phase_label].species for spec in rxn.products]):
                rxnlist.append(rxn)

        phasesys.phases[phase_label].reactions.extend(rxnlist)

        for key, interface in self.interfaces.items():
            rxnlist = []
            for i, rxn in enumerate(interface.reactions):
                if (spc in rxn.reactants or spc in rxn.products) and all([spec in phasesys.interfaces[key].species for spec in rxn.reactants]) and all([spec in phasesys.interfaces[key].species for spec in rxn.products]):
                    rxnlist.append(rxn)

            phasesys.interfaces[key].reactions.extend(rxnlist)

        return

    def get_species(self, label):
        """
        Retrieve rms species associated with the input label
        """
        for plabel, phase in self.phases.items():
            spc = phase.get_species(label)
            if spc is not None:
                return spc
        else:
            return None

class Phase:
    """
    Class containing all species, reactions and properties necessary to describe
    kinetics within a specific phase of a simulation
    """
    def __init__(self, label="", solvent=None, site_density=None):
        self.species = []
        self.reactions = []
        self.names = []
        self.species_dict = dict()
        if solvent:
            self.solvent = to_rms(solvent)
        if site_density:
            self.site_density = site_density

    def get_species(self, label):
        """
        Retrieve rms species associated with the input label
        """
        try:
            ind = self.names.index(label)
            return self.species[ind]
        except IndexError:
            return None

    def set_solvent(self, solvent):
        """
        Set the solvent of the phase
        """
        self.solvent = to_rms(solvent)

    def add_reaction(self, rxn, species_list):
        """
        add a reaction to the phase
        """
        self.reactions.append(to_rms(rxn, species_list=species_list, rms_species_list=self.species))

    def add_species(self, spc, edge_phase=None):
        """
        add a species to the phase
        if an edge_phase is given the
        edge phase species order is reordered to match
        that of self
        """
        if spc.label in self.names: #already exists
            label = spc.label
            logging.debug(f"species {label} was already in phase skipping...")
            return

        label = spc.label
        spec = to_rms(spc)
        self.species.append(spec)
        self.names.append(spc.label)

        if edge_phase is not None:
            edge_phase.species.insert(len(self.species)-1, spec)
            edge_phase.names.insert(len(self.species)-1, label)

    def remove_species(self, label):
        """
        Remove species and associated reactions from the phase
        """
        try:
            ind = self.names.index(label)
        except ValueError:
            return
        spc = self.species[ind]
        rxninds = []
        for i, rxn in enumerate(self.reactions):
            for spc2 in itertools.chain(rxn.reactants, rxn.products):
                if label == spc2.name:
                    rxninds.append(i)
                    break

        del self.species[ind]
        del self.names[ind]

        for ind in reversed(rxninds):
            del self.reactions[ind]

        return spc

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
