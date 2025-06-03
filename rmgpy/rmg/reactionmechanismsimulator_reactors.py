#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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
"""
import itertools
import logging
import os

import numpy as np

from rmgpy.data.solvation import SolventData
from rmgpy.kinetics.arrhenius import (
    Arrhenius,
    ArrheniusChargeTransfer,
    Marcus,
    MultiArrhenius,
    MultiPDepArrhenius,
    PDepArrhenius,
)
from rmgpy.kinetics.chebyshev import Chebyshev
from rmgpy.kinetics.falloff import Lindemann, ThirdBody, Troe
from rmgpy.kinetics.surface import StickingCoefficient, SurfaceChargeTransfer
from rmgpy.molecule.fragment import Fragment
from rmgpy.reaction import Reaction
from rmgpy.solver.termination import (
    TerminationConversion,
    TerminationRateRatio,
    TerminationTime,
)
from rmgpy.species import Species
from rmgpy.thermo.nasa import NASA, NASAPolynomial
from rmgpy.thermo.thermodata import ThermoData


NO_JULIA = True

# Check if Julia should be disabled via module-level variable or environment variable
import rmgpy
if hasattr(rmgpy, 'DISABLE_JULIA') and rmgpy.DISABLE_JULIA:
    logging.info("Julia imports disabled via rmgpy.DISABLE_JULIA flag")
elif os.environ.get('RMG_DISABLE_JULIA', '').lower() in ('1', 'true', 'yes'):
    logging.info("Julia imports disabled via RMG_DISABLE_JULIA environment variable")
else:
    try:
        from juliacall import Main
        Main.seval("using PythonCall")
        Main.seval("using ReactionMechanismSimulator")
        Main.seval("using ReactionMechanismSimulator.Sundials")
        NO_JULIA = False
    except ImportError as e:
        logging.warning("Julia import failed, RMS reactors not available. "
                       "Exception: %s", str(e))


def to_julia(obj):
    """
    Convert native Python object to julia object. If the object is a Python dict, it will be converted to a Julia Dict. If the object is a Python list, it will be converted to a Julia Vector. If the object is a 1-d numpy array, it will be converted to a Julia Array. If the object is a n-d (n > 1) numpy array, it will be converted to a Julia Matrix.
    Otherwise, the object will be returned as is. Doesn't support RMG objects.

    Parameters
    ----------
    obj : dict | list | np.ndarray | object
        The native Python object to convert

    Returns
    -------
    object : Main.Dict | Main.Vector | Main.Matrix | object
        The Julia object
    """
    if isinstance(obj, dict):
        return Main.PythonCall.pyconvert(Main.Dict, obj)
    elif isinstance(obj, (list, np.ndarray)):
        if getattr(obj, "shape", False) and len(obj.shape) > 1:
            return Main.PythonCall.pyconvert(Main.Matrix, obj)
        return Main.PythonCall.pyconvert(Main.Vector, obj)
    else: # Other native Python project does not need special conversion.
        return obj


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

        if spc.label in self.species_dict:
            del self.species_dict[spc.label]

    def add_reaction(self, rxn, species_list):
        """
        adds reaction to the appropriate phase/interface
        within PhaseSystem
        """
        spclabels = [r.label for r in rxn.reactants + rxn.products]
        phaseinv = []
        for plabel, phase in self.phases.items():
            out = [label for label in spclabels if label in phase.names]
            if len(out) > 0:
                phaseinv.append(phase)

        rms_species_list = self.get_rms_species_list()

        if len(phaseinv) == 1:
            phaseinv[0].add_reaction(rxn, species_list, rms_species_list)
        else:
            phases = set(phaseinv)
            for interface in self.interfaces:
                if interface.phaseset == phases:
                    interface.add_reaction(rxn, species_list, rms_species_list)
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
        self.phases[phase_label].species.insert(len(phasesys.phases[phase_label].species) - 1, spc)
        self.phases[phase_label].names.insert(len(phasesys.phases[phase_label].species) - 1, label)

        rxnlist = []
        for i, rxn in enumerate(self.phases[phase_label].reactions):
            reactants = Main.pylist(rxn.reactants)
            products = Main.pylist(rxn.products)
            reacs_and_prods = reactants + products
            if (spc.name in [spec.name for spec in reacs_and_prods]) and all(
                [spec.name in phasesys.species_dict for spec in reacs_and_prods]
            ):
                rxnlist.append(rxn)

        for i, rxn in enumerate(rxnlist):
            phasesys.phases[phase_label].reactions.append(rxn)
            self.phases[phase_label].reactions.remove(rxn)
            self.phases[phase_label].reactions.insert(len(phasesys.phases[phase_label].reactions) - 1, rxn)

        for key, interface in self.interfaces.items():
            rxnlist = []
            for i, rxn in enumerate(interface.reactions):
                if (
                    (spc in rxn.reactants or spc in rxn.products)
                    and all(spec.name in phasesys.species_dict for spec in rxn.reactants)
                    and all(spec.name in phasesys.species_dict for spec in rxn.products)
                ):
                    rxnlist.append(rxn)

            for i, rxn in enumerate(rxnlist):
                phasesys.interfaces[key].reactions.append(rxn)
                self.interfaces[key].reactions.remove(rxn)
                self.interfaces[key].reactions.insert(len(phasesys.interfaces[key].reactions) - 1, rxn)

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

    def get_rms_species_list(self):
        rms_species_list = []
        for phase in self.phases.values():
            rms_species_list += phase.species

        return rms_species_list

    def get_species_names(self):
        names = []
        for phase in self.phases.values():
            names += phase.names
        return names


class Phase:
    """
    Class containing all species, reactions and properties necessary to describe
    kinetics within a specific phase of a simulation
    """

    def __init__(self, label="", solvent=None, site_density=None):
        self.label = label
        self.species = []
        self.reactions = []
        self.names = []
        self.species_dict = dict()
        self.add_later_reactions = []
        self.rmg_species = []
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
        except (IndexError, ValueError):
            return None

    def set_solvent(self, solvent):
        """
        Set the solvent of the phase
        """
        self.solvent = to_rms(solvent)

    def add_reaction(self, rxn, edge_phase=None):
        """
        add a reaction to the phase
        """
        reactions_to_remove = []
        if self.add_later_reactions != []:
            for reaction in self.add_later_reactions:
                added = False
                try:
                    rms_rxn = to_rms(reaction, species_names=self.names, rms_species_list=self.species, rmg_species=self.rmg_species)
                    self.reactions.append(rms_rxn)
                    reactions_to_remove.append(reaction)
                    added = True
                except ValueError:
                    pass

                if added and edge_phase is not None:
                    # If edge phase is provided, then self is the core phase, and this rxn was directly added to the core so not yet exist in edge
                    # Only add the corresponding rms_rxn to edge if it's successfully added to the core, otherwise this edge phase reaction maybe passed from the edge phase to core phase in pass_species and double added here again
                    edge_phase.reactions.insert(len(self.reactions) - 1, rms_rxn)

            for reaction in reactions_to_remove:
                self.add_later_reactions.remove(reaction)

        added = False
        try:
            rms_rxn = to_rms(rxn, species_names=self.names, rms_species_list=self.species, rmg_species=self.rmg_species)
            self.reactions.append(rms_rxn)
            added = True
        except ValueError:  # often reactions with efficiencies from seed mechanisms can't be fully constructed until input species are added
            if rxn not in self.add_later_reactions:
                self.add_later_reactions.append(rxn)

        if added and edge_phase is not None:
            edge_phase.reactions.insert(len(self.reactions) - 1, rms_rxn)

        return

    def add_species(self, spc, edge_phase=None):
        """
        add a species to the phase
        if an edge_phase is given the
        edge phase species order is reordered to match
        that of self
        """
        if spc.label in self.names:  # already exists
            label = spc.label
            logging.debug(f"species {label} was already in phase skipping...")
            return

        self.rmg_species.append(spc)

        label = spc.label
        spec = to_rms(spc)
        self.species.append(spec)
        self.names.append(spc.label)

        if edge_phase is not None:
            edge_phase.species.insert(len(self.species) - 1, spec)
            edge_phase.names.insert(len(self.species) - 1, label)

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


class Interface:
    """
    Class containing all reactions and properties necessary to describe
    kinetics within an interface between two phases in a simulation
    """

    def __init__(self, phases, reactions=[]):
        self.reactions = reactions or []
        self.phaseset = set(phases)

    def add_reaction(self, rxn, species_names, rms_species_list):
        """
        add a reaction to the interface
        """
        self.reactions.append(to_rms(rxn, species_names=species_names, rms_species_list=rms_species_list))

    def remove_species(self, spc):
        """
        Remove reactions associated with the input species from the interface
        """
        rxninds = []
        for i, rxn in enumerate(self.reactions):
            for spc2 in itertools.chain(rxn.reactants, rxn.products):
                if spc.name == spc2.name:
                    rxninds.append(i)
                    break
        for ind in reversed(rxninds):
            del self.reactions[ind]
        return


class Reactor:
    """
    Super class for running edge analysis
    Subclasses implement generate_reactor methods for
    generating the necessary RMS phase/domain/reactor objects
    """

    def __init__(self, core_phase_system, edge_phase_system, initial_conditions, terminations, constant_species=[]):
        self.core_phase_system = core_phase_system
        self.edge_phase_system = edge_phase_system
        self.initial_conditions = initial_conditions
        self.const_spc_names = constant_species
        self.n_sims = 1
        self.tf = 1.0e6
        for term in terminations:
            if isinstance(term, TerminationTime):
                self.tf = term.time.value_si

        self.terminations = [to_rms(term) for term in terminations]

    def get_const_spc_indices(self, spcs=None):
        rms_species_names = [spc.name for spc in self.core_phase_system.get_rms_species_list()]
        return [rms_species_names.index(name) for name in self.const_spc_names]

    def finish_termination_criteria(self):
        """
        Convert tuples into TerminationConversion objects
        this is necessary because in input.py some species objects
        are created before they can be converted to rms Species objects
        so we construct the rms TerminationConversion objects later
        """
        for i, term in enumerate(self.terminations):
            if isinstance(term, tuple):
                self.terminations[i] = to_rms(TerminationConversion(term[0], term[1]))

    def reset_max_edge_species_rate_ratios(self):
        """
        This function sets max_edge_species_rate_ratios back to zero
        for pruning of ranged reactors it is important to avoid doing this
        every initialization
        """
        self.max_edge_species_rate_ratios = np.zeros((len(self.prunable_species)), float)

    def simulate(self, model_settings, simulator_settings, conditions):
        """
        Run edge analysis of the reactor system
        """

        # We have to build the edge phase first, because the efficiencies dictionary in any kinetics that has
        # efficiencies gets re-written every time we build a phase. If we build the edge phase later, the
        # nameefficiencies will contain edge species, and the efficienies dictionary will have edge species index,
        # which ultimately cause the core simulation to crash.
        edge_react, edge_domains, edge_interfaces, edge_p = self.generate_reactor(self.edge_phase_system)
        core_react, core_domains, core_interfaces, core_p = self.generate_reactor(self.core_phase_system)

        (
            terminated,
            resurrected,
            invalid_objects,
            unimolecular_threshold,
            bimolecular_threshold,
            trimolecular_threshold,
            max_edge_species_rate_ratios,
            t,
            x,
        ) = Main.selectobjects(
            core_react,
            edge_react,
            edge_domains,
            edge_interfaces,
            core_domains,
            core_interfaces,
            core_p,
            edge_p,
            model_settings.tol_move_to_core,
            model_settings.tol_rad_move_to_core,
            model_settings.tol_interrupt_simulation,
            model_settings.ignore_overall_flux_criterion,
            model_settings.filter_reactions,
            model_settings.max_num_objects_per_iter,
            model_settings.tol_branch_rxn_to_core,
            model_settings.branching_ratio_max,
            model_settings.branching_index,
            model_settings.terminate_at_max_objects,
            self.terminations,
            model_settings.filter_threshold,
            model_settings.transitory_tol_dict,
            model_settings.transitory_step_period,
            model_settings.tol_rxn_to_core_deadend_radical,
            atol=simulator_settings.atol,
            rtol=simulator_settings.rtol,
            solver=Main.Sundials.CVODE_BDF(linear_solver=Main.Symbol("GMRES")),
        )

        return (
            terminated,
            resurrected,
            invalid_objects,
            unimolecular_threshold,
            bimolecular_threshold,
            trimolecular_threshold,
            max_edge_species_rate_ratios,
            t,
            x,
        )


class ConstantVIdealGasReactor(Reactor):
    def __init__(self, core_phase_system, edge_phase_system, initial_conditions, terminations, constant_species=[]):
        super().__init__(core_phase_system, edge_phase_system, initial_conditions, terminations, constant_species=[])

    def generate_reactor(self, phase_system):
        """
        Setup an RMS simulation for EdgeAnalysis
        """
        phase = phase_system.phases["Default"]
        ig = Main.IdealGas(phase.species, phase.reactions)
        domain, y0, p = Main.ConstantVDomain(phase=ig, initialconds=to_julia(self.initial_conditions))
        react = Main.Reactor(domain, y0, (0.0, self.tf), p=p)
        return react, domain, [], p


class ConstantTLiquidSurfaceReactor(Reactor):
    def __init__(self, core_phase_system, edge_phase_system, initial_conditions, terminations, constant_species):
        super().__init__(core_phase_system, edge_phase_system, initial_conditions, terminations, constant_species)

    def generate_reactor(self, phase_system):
        """
        Setup an RMS simulation for EdgeAnalysis
        """
        liq = phase_system.phases["Default"]
        surf = phase_system.phases["Surface"]
        interface = list(phase_system.interfaces.values())[0]
        if "mu" in self.initial_conditions["liquid"].keys():
            solv = Main.Solvent("solvent", Main.ConstantViscosity(self.initial_conditions["liquid"]["mu"]))
            liq_initial_cond = self.initial_conditions["liquid"].copy()
            del liq_initial_cond["mu"]
        else:
            solv = liq.solvent
            liq_initial_cond = self.initial_conditions["liquid"]
        liq = Main.IdealDiluteSolution(liq.species, liq.reactions, solv, name="liquid",diffusionlimited=True)
        surf = Main.IdealSurface(surf.species, surf.reactions, surf.site_density, name="surface")
        liq_constant_species = [cspc for cspc in self.const_spc_names if cspc in [spc.name for spc in liq.species]]
        cat_constant_species = [cspc for cspc in self.const_spc_names if cspc in [spc.name for spc in surf.species]]
        domainliq, y0liq, pliq = Main.ConstantTVDomain(phase=liq, initialconds=to_julia(liq_initial_cond), constantspecies=to_julia(liq_constant_species))
        domaincat, y0cat, pcat = Main.ConstantTAPhiDomain(
            phase=surf, initialconds=to_julia(self.initial_conditions["surface"]), constantspecies=to_julia(cat_constant_species),
        )
        if interface.reactions == []:
            inter, pinter = Main.ReactiveInternalInterfaceConstantTPhi(
                domainliq,
                domaincat,
                Main.seval("using ReactionMechanismSimulator; Vector{ElementaryReaction}()"),
                self.initial_conditions["liquid"]["T"],
                self.initial_conditions["surface"]["A"],
            )
        else:
            inter, pinter = Main.ReactiveInternalInterfaceConstantTPhi(
                domainliq, domaincat, interface.reactions, self.initial_conditions["liquid"]["T"], self.initial_conditions["surface"]["A"]
            )
        react, y0, p = Main.Reactor((domainliq, domaincat), (y0liq, y0cat), (0.0, self.tf), (inter,), (pliq, pcat, pinter))
        return react, (domainliq, domaincat), (inter,), p


class ConstantTVLiquidReactor(Reactor):
    def __init__(
        self,
        core_phase_system,
        edge_phase_system,
        initial_conditions,
        terminations,
        constant_species=[],
        inlet_conditions=dict(),
        outlet_conditions=dict(),
        evap_cond_conditions=dict(),
    ):
        super().__init__(core_phase_system, edge_phase_system, initial_conditions, terminations, constant_species=constant_species)

        self.inlet_conditions = inlet_conditions
        self.outlet_conditions = outlet_conditions
        self.evap_cond_conditions = evap_cond_conditions

    def generate_reactor(self, phase_system):
        """
        Setup an RMS simulation for EdgeAnalysis
        """
        phase = phase_system.phases["Default"]
        liq = Main.IdealDiluteSolution(phase.species, phase.reactions, phase.solvent)
        domain, y0, p = Main.ConstantTVDomain(phase=liq, initialconds=to_julia(self.initial_conditions), constantspecies=to_julia(self.const_spc_names))

        interfaces = []

        if self.inlet_conditions:
            inlet_conditions = {key: value for (key, value) in self.inlet_conditions.items() if key != "F"}
            total_molar_flow_rate = self.inlet_conditions["F"]
            inlet = Main.Inlet(domain, to_julia(inlet_conditions), Main.seval("x->" + str(total_molar_flow_rate)))
            interfaces.append(inlet)

        if self.outlet_conditions:
            total_volumetric_flow_rate = self.outlet_conditions["Vout"]
            outlet = Main.VolumetricFlowRateOutlet(domain, Main.seval("x->" + str(total_volumetric_flow_rate)))
            interfaces.append(outlet)

        if self.evap_cond_conditions:
            kLA_kH_evap_cond = Main.kLAkHCondensationEvaporationWithReservoir(domain, self.evap_cond_conditions)
            interfaces.append(kLA_kH_evap_cond)

        react = Main.Reactor(domain, y0, (0.0, self.tf), interfaces, p=p)
        return react, domain, interfaces, p


class ConstantTPIdealGasReactor(Reactor):
    def __init__(self, core_phase_system, edge_phase_system, initial_conditions, terminations, constant_species=None):
        if constant_species is None:
            constant_species = []
        super().__init__(core_phase_system, edge_phase_system, initial_conditions, terminations, constant_species=[])

    def generate_reactor(self, phase_system):
        """
        Setup an RMS simulation for EdgeAnalysis
        """
        phase = phase_system.phases["Default"]
        ig = Main.IdealGas(phase.species, phase.reactions)
        domain, y0, p = Main.ConstantTPDomain(phase=ig, initialconds=to_julia(self.initial_conditions))
        react = Main.Reactor(domain, y0, (0.0, self.tf), p=p)
        return react, domain, [], p


def to_rms(obj, species_names=None, rms_species_list=None, rmg_species=None):
    """
    Generate corresponding rms object
    """
    if isinstance(obj, ThermoData):
        obj = obj.to_nasa(Tmin=298, Tmax=2500, Tint=1000)
    if isinstance(obj, Arrhenius):
        if obj._T0.value_si != 1:
            A = obj._A.value_si / (obj._T0.value_si) ** obj._n.value_si
        else:
            A = obj._A.value_si
        n = obj._n.value_si
        Ea = obj._Ea.value_si
        return Main.Arrhenius(A, n, Ea, Main.EmptyRateUncertainty())
    elif isinstance(obj, ArrheniusChargeTransfer):
        A = obj._A.value_si
        if obj._T0.value_si != 1.0:
            A /= ((obj._T0.value_si) ** obj._n.value_si)
        n = obj._n.value_si
        Ea = obj._Ea.value_si
        q = obj._alpha.value_si*obj._electrons.value_si
        V0 = obj._V0.value_si
        return Main.Arrheniusq(A, n, Ea, q, V0, Main.EmptyRateUncertainty())
    elif isinstance(obj, SurfaceChargeTransfer):
        A = obj._A.value_si
        if obj._T0.value_si != 1.0:
            A /= ((obj._T0.value_si) ** obj._n.value_si)
        n = obj._n.value_si
        Ea = obj._Ea.value_si
        q = obj._alpha.value_si*obj._electrons.value_si
        V0 = obj._V0.value_si
        return Main.Arrheniusq(A, n, Ea, q, V0, Main.EmptyRateUncertainty())
    elif isinstance(obj, Marcus):
        A = obj._A.value_si
        n = obj._n.value_si
        return Main.Marcus(A,n,obj._lmbd_i_coefs.value_si,obj._lmbd_o.value_si, obj._wr.value_si, obj._wp.value_si, obj._beta.value_si, Main.EmptyRateUncertainty())
    elif isinstance(obj, PDepArrhenius):
        Ps = to_julia(obj._pressures.value_si)
        arrs = to_julia([to_rms(arr) for arr in obj.arrhenius])
        return Main.PdepArrhenius(Ps, arrs, Main.EmptyRateUncertainty())
    elif isinstance(obj, MultiArrhenius):
        arrs = to_julia([to_rms(arr) for arr in obj.arrhenius])
        return Main.MultiArrhenius(arrs, Main.EmptyRateUncertainty())
    elif isinstance(obj, MultiPDepArrhenius):
        parrs = to_julia([to_rms(parr) for parr in obj.arrhenius])
        return Main.MultiPdepArrhenius(parrs, Main.EmptyRateUncertainty())
    elif isinstance(obj, Chebyshev):
        Tmin = obj.Tmin.value_si
        Tmax = obj.Tmax.value_si
        Pmin = obj.Pmin.value_si
        Pmax = obj.Pmax.value_si
        coeffs = to_julia(obj.coeffs.value_si)
        return Main.Chebyshev(coeffs, Tmin, Tmax, Pmin, Pmax)
    elif isinstance(obj, ThirdBody):
        arrstr = arrhenius_to_julia_string(obj.arrheniusLow)
        efficiencies = {rmg_species[i].label: float(val) for i, val in enumerate(obj.get_effective_collider_efficiencies(rmg_species)) if val != 1}
        dstr = "Dict{String,Float64}(["
        for key, value in efficiencies.items():
            dstr += '"' + key + '"' "=>" + str(value) + ","
        dstr += "])"
        return Main.seval(
            "using ReactionMechanismSimulator; ThirdBody(" + arrstr + ", Dict{Int64,Float64}([]), " + dstr + "," + "EmptyRateUncertainty())"
        )
    elif isinstance(obj, Lindemann):
        arrlow = arrhenius_to_julia_string(obj.arrheniusLow)
        arrhigh = arrhenius_to_julia_string(obj.arrheniusHigh)
        efficiencies = {rmg_species[i].label: float(val) for i, val in enumerate(obj.get_effective_collider_efficiencies(rmg_species)) if val != 1}
        dstr = "Dict{String,Float64}(["
        for key, value in efficiencies.items():
            dstr += '"' + key + '"' "=>" + str(value) + ","
        dstr += "])"
        return Main.seval(
            "using ReactionMechanismSimulator; Lindemann("
            + arrhigh
            + ","
            + arrlow
            + ","
            + "Dict{Int64,Float64}([]),"
            + dstr
            + ","
            + "EmptyRateUncertainty())"
        )
    elif isinstance(obj, Troe):
        arrlow = arrhenius_to_julia_string(obj.arrheniusLow)
        arrhigh = arrhenius_to_julia_string(obj.arrheniusHigh)
        efficiencies = {rmg_species[i].label: float(val) for i, val in enumerate(obj.get_effective_collider_efficiencies(rmg_species)) if val != 1}
        alpha = obj.alpha
        T1 = obj._T1.value_si if obj._T1 is not None else 0.0
        T2 = obj._T2.value_si if obj._T2 is not None else 0.0
        T3 = obj._T3.value_si if obj._T3 is not None else 0.0
        dstr = "Dict{String,Float64}(["
        for key, value in efficiencies.items():
            dstr += '"' + key + '"' "=>" + str(value) + ","
        dstr += "])"
        return Main.seval(
            "using ReactionMechanismSimulator; Troe("
            + arrhigh
            + ","
            + arrlow
            + ","
            + str(alpha)
            + ","
            + str(T3)
            + ","
            + str(T1)
            + ","
            + str(T2)
            + ","
            + "Dict{Int64,Float64}([]),"
            + dstr
            + ","
            + "EmptyRateUncertainty())"
        )
    elif isinstance(obj, StickingCoefficient):
        if obj._T0.value_si != 1:
            A = obj._A.value_si / (obj._T0.value_si) ** obj._n.value_si
        else:
            A = obj._A.value_si
        n = obj._n.value_si
        Ea = obj._Ea.value_si
        return Main.StickingCoefficient(A, n, Ea, Main.EmptyRateUncertainty())
    elif isinstance(obj, NASAPolynomial):
        return Main.NASApolynomial(obj.coeffs, obj.Tmin.value_si, obj.Tmax.value_si)
    elif isinstance(obj, NASA):
        return Main.NASA([to_rms(poly) for poly in obj.polynomials], Main.EmptyThermoUncertainty())
    elif isinstance(obj, Species):

        if isinstance(obj.molecule[0], Fragment):
            obj.molecule[0].assign_representative_molecule()
            mol = obj.molecule[0].mol_repr
        else:
            mol = obj.molecule[0]

        atomnums = dict()
        for atm in mol.atoms:
            if atomnums.get(atm.element.symbol):
                atomnums[atm.element.symbol] += 1
            else:
                atomnums[atm.element.symbol] = 1
        atomnums = to_julia(atomnums)
        bondnum = len(mol.get_all_edges())

        if not obj.molecule[0].contains_surface_site():
            rad = Main.getspeciesradius(atomnums, bondnum)
            diff = Main.StokesDiffusivity(rad)
            th = obj.get_thermo_data()
            thermo = to_rms(th)
            if obj.henry_law_constant_data:
                kH = Main.TemperatureDependentHenryLawConstant(Ts=to_julia(obj.henry_law_constant_data.Ts), kHs=to_julia(obj.henry_law_constant_data.kHs))
            else:
                kH = Main.EmptyHenryLawConstant()
            if obj.liquid_volumetric_mass_transfer_coefficient_data:
                kLA = Main.TemperatureDependentLiquidVolumetricMassTransferCoefficient(
                    Ts=to_julia(obj.liquid_volumetric_mass_transfer_coefficient_data.Ts), kLAs=to_julia(obj.liquid_volumetric_mass_transfer_coefficient_data.kLAs)
                )
            else:
                kLA = Main.EmptyLiquidVolumetricMassTransferCoefficient()
            return Main.Species(
                name=obj.label,
                index=obj.index,
                inchi="",
                smiles="",
                adjlist="",
                thermo=thermo,
                atomnums=atomnums,
                bondnum=bondnum,
                diffusion=diff,
                radius=rad,
                radicalelectrons=obj.molecule[0].multiplicity - 1,
                molecularweight=obj.molecular_weight.value_si,
                henrylawconstant=kH,
                liquidvolumetricmasstransfercoefficient=kLA,
                comment=obj.thermo.comment,
            )
        else:
            th = obj.get_thermo_data()
            thermo = to_rms(th)
            return Main.Species(
                name=obj.label,
                index=obj.index,
                inchi="",
                smiles="",
                adjlist="",
                thermo=thermo,
                atomnums=atomnums,
                bondnum=bondnum,
                diffusion=Main.EmptyDiffusivity(),
                radius=0.0,
                radicalelectrons=obj.molecule[0].multiplicity - 1,
                molecularweight=0.0,
                henrylawconstant=Main.EmptyHenryLawConstant(),
                liquidvolumetricmasstransfercoefficient=Main.EmptyLiquidVolumetricMassTransferCoefficient(),
                comment=obj.thermo.comment,
            )
    elif isinstance(obj, Reaction):
        reactantinds = to_julia([species_names.index(spc.label) for spc in obj.reactants])
        productinds = to_julia([species_names.index(spc.label) for spc in obj.products])
        reactants = to_julia([rms_species_list[i] for i in reactantinds])
        products = to_julia([rms_species_list[i] for i in productinds])
        if isinstance(obj.kinetics, SurfaceChargeTransfer):
            obj.set_reference_potential(300)
        kinetics = to_rms(obj.kinetics, species_names=species_names, rms_species_list=rms_species_list, rmg_species=rmg_species)
        radchange = sum([spc.molecule[0].multiplicity-1 for spc in obj.products]) - sum([spc.molecule[0].multiplicity-1 for spc in obj.reactants])
        electronchange = -sum([spc.molecule[0].get_net_charge() for spc in obj.products]) + sum([spc.molecule[0].get_net_charge() for spc in obj.reactants])
        return Main.ElementaryReaction(index=obj.index, reactants=reactants, reactantinds=reactantinds, products=products, productinds=productinds, kinetics=kinetics, electronchange=electronchange, radicalchange=radchange, reversible=obj.reversible, pairs=to_julia([]), comment=obj.kinetics.comment)
    elif isinstance(obj, SolventData):
        return Main.Solvent("solvent", Main.RiedelViscosity(float(obj.A), float(obj.B), float(obj.C), float(obj.D), float(obj.E)))
    elif isinstance(obj, TerminationTime):
        return Main.TerminationTime(obj.time.value_si)
    elif isinstance(obj, TerminationConversion):
        return Main.TerminationConversion(to_rms(obj.species), obj.conversion)
    elif isinstance(obj, TerminationRateRatio):
        return Main.TerminationRateRatio(obj.ratio)
    elif isinstance(obj, tuple):  # Handle TerminationConversion when the obj doesn't have thermo yet
        return obj
    else:
        errortype = type(obj)
        raise ValueError(f"Couldn't convert object of type {errortype} to an RMS object")


def arrhenius_to_julia_string(obj):
    return "Arrhenius(" + str(obj.A.value_si) + "," + str(obj.n.value_si) + "," + str(obj.Ea.value_si) + ", EmptyRateUncertainty())"
