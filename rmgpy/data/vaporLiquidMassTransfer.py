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
This module contains classes related to properties that are used to calculate vapor liquid mass transfer coefficients
"""

import logging
import numpy as np

from rmgpy.rmgobject import RMGObject

class VaporLiquidMassTransfer(RMGObject):
    """
    Class for working with coefficients related to vapor liquid mass transfer
    
    The attributes are:
    =========================================================  ===================================================================
    Attribute                                                  Description
    =========================================================  ===================================================================
    `enabled`                                                  Boolean for whether the vapor liquid mass transfer coefficients should be calculated
    `solvent_data`                                             Solvent data for calculating vapor liquid mass transfer coefficients
    `database`                                                 Solvation database to obtain the solute data
    `liquid_volumetric_mass_transfer_coefficient_power_law`    Power law coefficients for kLA calculations
    =========================================================  ===================================================================
    """

    def __init__(self):
        # default is false, enabled if there is a solvent and vapor liquid mass transfer power law are provided
        self.enabled = False

    def enable(self, solvent_data, solvation_database, liquid_volumetric_mass_transfer_coefficient_power_law):
        # vapor_liquid_mass_transfer is enabled if a solvent and vapor liquid mass transfer power law have been added to the RMG object.
        logging.info("Enabling vapor liquid mass transfer...")
        vapor_liquid_mass_transfer.enabled = True
        vapor_liquid_mass_transfer.solvent_data = solvent_data
        vapor_liquid_mass_transfer.database = solvation_database
        vapor_liquid_mass_transfer.liquid_volumetric_mass_transfer_coefficient_power_law = liquid_volumetric_mass_transfer_coefficient_power_law

    def get_liquid_volumetric_mass_transfer_coefficient_data(self, spec, Ts=[]):
        solute_data = self.database.get_solute_data(spec)
        solvent_data = self.solvent_data
        prefactor = self.liquid_volumetric_mass_transfer_coefficient_power_law.prefactor
        diffusion_coefficient_power = self.liquid_volumetric_mass_transfer_coefficient_power_law.diffusion_coefficient_power
        solvent_viscosity_power = self.liquid_volumetric_mass_transfer_coefficient_power_law.solvent_viscosity_power
        solvent_density_power = self.liquid_volumetric_mass_transfer_coefficient_power_law.solvent_density_power 

        if not Ts:
            Tmin = solvent_data.get_solvent_coolprop_Tmin()
            Tcrit = solvent_data.get_solvent_coolprop_Tcrit()
            Ts = [float(T) for T in np.linspace(Tmin, Tcrit-0.01, 50)]
        kLAs = [prefactor * solute_data.get_stokes_diffusivity(T, solvent_data.get_solvent_viscosity(T))**diffusion_coefficient_power * solvent_data.get_solvent_viscosity(T)**solvent_viscosity_power * solvent_data.get_solvent_density(T)**solvent_density_power for T in Ts]

        return LiquidVolumetricMassTransferCoefficientData(Ts=Ts,kLAs=kLAs)

    def get_henry_law_constant_data(self, spec, Ts=[]):
        solute_data = self.database.get_solute_data(spec)
        solvent_data = self.solvent_data

        if not Ts:
            Tmin = solvent_data.get_solvent_coolprop_Tmin()
            Tcrit = solvent_data.get_solvent_coolprop_Tcrit()
            Ts = [float(T) for T in np.linspace(Tmin, Tcrit-0.01, 50)]

        # The function `get_T_dep_solvation_energy_from_LSER_298`` returns (delG, Kfactor(T,Psat), kH(T))
        kHs = [self.database.get_T_dep_solvation_energy_from_LSER_298(solute_data, solvent_data, T)[2] for T in Ts]

        return HenryLawConstantData(Ts=Ts,kHs=kHs)

class LiquidVolumetricMassTransferCoefficientData(RMGObject):
    """
    Sampled liquid volumetric mass transfer coefficients.
    
    The attributes are:
    =================  ============================================================
    Attribute          Description
    =================  ============================================================
    `Ts`        		Sampled temperatures
    `kLAs`           	Liquid volumetric mass transfer coefficients corresponding to Ts
    =================  ============================================================
    """

    def __init__(self, Ts=[], kLAs=[]) -> None:
        super().__init__()
        self.Ts = Ts
        self.kLAs = kLAs

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        LiquidVolumetricMassTransferCoefficientData object.
        """
        attributes = []
        if self.Ts:
            attributes.append(f'Ts={self.Ts}')
        if self.kLAs:
            attributes.append(f'kLAs={self.kLAs}')
        string = 'LiquidVolumetricMassTransferCoefficientData({0!s})'.format(', '.join(attributes))
        return string

class HenryLawConstantData(RMGObject):
    """
    Sampled Henry's law constants.
    
    The attributes are:
    =================  ============================================================
    Attribute          Description
    =================  ============================================================
    `Ts`        		Sampled temperatures
    `kHs`           	Henry's law constants corresponding to Ts
    =================  ============================================================
    """

    def __init__(self, Ts=[], kHs=[]):
        self.Ts = Ts
        self.kHs = kHs

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        HenryLawConstantData object.
        """
        attributes = []
        if self.Ts:
            attributes.append(f'Ts={self.Ts}')
        if self.kHs:
            attributes.append(f'kHs={self.kHs}')
        string = 'HenryLawConstantData({0!s})'.format(', '.join(attributes))
        return string

class liquidVolumetricMassTransferCoefficientPowerLaw(RMGObject):
    """
    Power law coefficients used to calculation liquid volumetric mass transfer coefficient (kLA).
    kLA for species i with solvent solv is calculated by
    kLA_i = prefactor * D_i ^ diffusion_coefficient_power * mu_solv ^ solvent_viscosity_power * rho_solv ^ solvent_density_power
    
    The attributes are:
    =============================  ============================================================
    Attribute                      Description
    =============================  ============================================================
    `prefactor`        	           prefactor used in kLA calculation in SI unit
    `diffusion_coefficient_power`  Power law coefficient for species diffusion coefficient
    `solvent_viscosity_power`      Power law coefficient for solvent viscosity
    `solvent_density_power`        Power law coefficient for solvent density
    =============================  ============================================================
    
    """
    def __init__(self, prefactor=0, diffusion_coefficient_power=0, solvent_viscosity_power=0, solvent_density_power=0):
        self.prefactor = prefactor
        self.diffusion_coefficient_power = diffusion_coefficient_power
        self.solvent_viscosity_power = solvent_viscosity_power
        self.solvent_density_power = solvent_density_power

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        liquidVolumetricMassTransferCoefficientPowerLaw object.
        """
        attributes = []
        if self.prefactor:
            attributes.append(f'prefactor={self.prefactor}')
        if self.diffusion_coefficient_power:
            attributes.append(f'diffusion_coefficient_power={self.diffusion_coefficient_power}')
        if self.solvent_viscosity_power:
            attributes.append(f'solvent_viscosity_power={self.solvent_viscosity_power}')
        if self.solvent_density_power:
            attributes.append(f'solvent_density_power={self.solvent_density_power}')
        string = 'liquidVolumetricMassTransferCoefficientPowerLaw({0!s})'.format(', '.join(attributes))
        return string

# module level variable. There should only ever be one. It starts off disabled
vapor_liquid_mass_transfer = VaporLiquidMassTransfer()