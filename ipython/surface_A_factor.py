#!/usr/bin/env python3

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.constants import kB,h,R
from rmgpy.thermo.thermoengine import submit
from rmgpy.kinetics import get_rate_coefficient_units_from_reaction_order

import logging
import numpy as np

def get_A_units(rxn):
    
    try:
        surf_reacts = [spcs for spcs in rxn.reactants if spcs.contains_surface_site()]
    except IndexError:
        surf_prods = []
        logging.warning(f"Species do not have an rmgpy.molecule.Molecule "  
                        "Cannot determine phases of species. We will assume gas"
                        )
    n_surf = len(surf_reacts)
    n_gas = len(rxn.reactants) - len(surf_reacts)
    return get_rate_coefficient_units_from_reaction_order(n_gas, n_surf)

def get_surface_reaction_stoich(rxn):
    
    n_gas = 0
    n_ads = 0
    n_x = 0
    
    for s in rxn.reactants:
        if s.is_surface_site():
            n_x -= 1
        elif s.contains_surface_site():
            n_ads -= 1
        else:
            n_gas -= 1
    for s in rxn.products:
        if s.is_surface_site():
            n_x += 1
        elif s.contains_surface_site():
            n_ads += 1
        else:
            n_gas += 1
            
    return (n_gas,n_ads,n_x)

def get_surface_reaction_type(rxn):
    
    n_gas, n_ads, n_x = get_surface_reaction_stoich(rxn)
    
    if (n_gas, n_ads, n_x) == (1,-1,1):
        return "desorption"
    elif (n_gas, n_ads, n_x) == (-1,1,-1):
        return "adsorption"
    elif (n_gas, n_ads, n_x) == (0,1,-1):
        return "dissociation"
    elif (n_gas, n_ads, n_x) == (0,-1,1):
        return "association"
    else:
        return "unknown"

def estimate_surface_desorption_A(adsorbate: Species):

    Ar_mw = 39.87750372310267 # amu
    dummy_species = get_dummy_species(adsorbate)
    S = dummy_species.get_entropy(298.)
    mw = dummy_species.molecular_weight.value

    A = kB*298./h
    A *= np.exp(0.3*S/R+3.3-(18.6+np.log((mw/Ar_mw)**(5./2.)))/3)
    commment = f"A factor estimated from gas-phase smiles {dummy_species.smiles} from "\
    f"{dummy_species.thermo.comment} and S298={S/4.184:.2f} cal/mol/K"
    return A,commment

def get_dummy_species(adsorbate: Species):

    dummy_molecules = adsorbate.molecule[0].get_desorbed_molecules()
    for mol in dummy_molecules:
        mol.clear_labeled_atoms()
        
    gas_phase_species_from_libraries = []
    gas_phase_species_estimates = []
    for dummy_molecule in dummy_molecules:
        dummy_species = Species()
        dummy_species.molecule = [dummy_molecule]
        dummy_species.generate_resonance_structures()
        submit(dummy_species)
        if dummy_species.thermo.label:
            gas_phase_species_from_libraries.append(dummy_species)
        else:
            gas_phase_species_estimates.append(dummy_species)

    # define the comparison function to find the lowest energy
    def lowest_energy(species):
        if hasattr(species.thermo, 'H298'):
            return species.thermo.H298.value_si
        else:
            return species.thermo.get_enthalpy(298.0)

    if gas_phase_species_from_libraries:
        species = min(gas_phase_species_from_libraries, key=lowest_energy)
    else:
        species = min(gas_phase_species_estimates, key=lowest_energy)
        
    return species

def estimate_surface_dissociation_A(adsorbate: Species):
    
    A, comment = estimate_surface_desorption_A(adsorbate)
    return 0.001 * A, comment

def get_adsorbate(rxn, on='reactants'):
    
    if on == 'reactants':
        for s in rxn.reactants:
            if s.contains_surface_site() and not s.is_surface_site():
                return s
    elif on == 'products':
        for s in rxn.products:
            if s.contains_surface_site() and not s.is_surface_site():
                return s

def estimate_surface_A(rxn):
    
    reaction_type = get_surface_reaction_type(rxn)
    comment = "A factor estimation:\n"
    
    if reaction_type == 'desorption':
        comment += f"A factor estimate for {reaction_type}\n"
        adsorbate = get_adsorbate(rxn,'reactants')
        A, _comment = estimate_surface_desorption_A(adsorbate)
        comment += _comment
        
    elif reaction_type == 'dissociation':
        comment += f"A factor estimate for {reaction_type}\n"
        adsorbate = get_adsorbate(rxn,'reactants')
        A, _comment = estimate_surface_dissociation_A(adsorbate)
        comment += _comment
    
    elif reaction_type == 'association':
        comment += f"A factor estimate for {reaction_type}\n"
        adsorbate = get_adsorbate(rxn,'products')
        A, _comment = estimate_surface_dissociation_A(adsorbate)
        comment += _comment
        if not adsorbate.thermo:
            submit(adsorbate)
        deltaS = adsorbate.thermo.get_entropy(298.)
        assert len(rxn.reactants) == 2
        for s in rxn.reactants:
            if s.contains_surface_site():
                if not s.thermo:
                    submit(s)
                deltaS -= s.thermo.get_entropy(298.)
        A *= np.exp(deltaS/R)
        
    else:
        A = kB*298./h
        comment = "Could not determine reaction type "\
        f"estimating A = kb/298/h = {A:.2e}"
        
    units = get_A_units(rxn)
    surface_site_density=2.483e-05 # Pt111 in metal DB
    if units in ('m^2/(mol*s)','m^5/(mol^2*s)'):
        A /= surface_site_density
        comment += '\nA/=2.483e-5 mol/m^2 (Pt111 site density)'
    elif units == 'm^4/(mol^2*s)':
        A /= surface_site_density
        A /= surface_site_density
        comment += '\nA/=(2.483e-5 mol/m^2)^2 (Pt111 site density)'

    return (A,units), comment
