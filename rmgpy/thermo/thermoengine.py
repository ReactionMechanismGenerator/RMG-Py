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

import logging as logging
import math

import numpy as np

import rmgpy.constants as constants
from rmgpy.data.rmg import get_db
from rmgpy.statmech import Conformer
from rmgpy.thermo import Wilhoit, NASA, ThermoData
from rmgpy.molecule import Molecule
from rmgpy.molecule.fragment import Fragment


def process_thermo_data(spc, thermo0, thermo_class=NASA, solvent_name=''):
    """
    Converts via Wilhoit into required `thermo_class` and sets `E0`.
    
    Resulting thermo is returned.
    """
    thermo = None

    # Always convert to Wilhoit so we can compute E0
    if isinstance(thermo0, Wilhoit):
        wilhoit = thermo0
    elif isinstance(thermo0, ThermoData):
        wilhoit = thermo0.to_wilhoit(B=1000.)
    else:
        wilhoit = thermo0.to_wilhoit()

    # Add on solvation correction
    solvation_database = get_db('solvation')
    if not solvent_name or solvation_database is None:
        logging.debug('Solvent database or solvent_name not found. Solvent effect was not utilized')
        solvent_data = None
    else:
        solvent_data = solvation_database.get_solvent_data(solvent_name)
    if solvent_data and not "Liquid thermo library" in thermo0.comment:
        solvation_database = get_db('solvation')
        solute_data = solvation_database.get_solute_data(spc)
        solvation_correction = solvation_database.get_solvation_correction(solute_data, solvent_data)
        # correction is added to the entropy and enthalpy
        wilhoit.S0.value_si = (wilhoit.S0.value_si + solvation_correction.entropy)
        wilhoit.H0.value_si = (wilhoit.H0.value_si + solvation_correction.enthalpy)
        wilhoit.comment += ' + Solvation correction with {} as solvent and solute estimated using {}'.format(solvent_name, solute_data.comment)

    # Compute E0 by extrapolation to 0 K
    if spc.conformer is None:
        spc.conformer = Conformer()
    spc.conformer.E0 = wilhoit.E0

    # Convert to desired thermo class
    if thermo_class is Wilhoit:
        thermo = wilhoit
    elif thermo_class is NASA:
        if solvent_data:
            # If liquid phase simulation keep the nasa polynomial if it comes from a liquid phase thermoLibrary.
            # Otherwise convert wilhoit to NASA
            if "Liquid thermo library" in thermo0.comment and isinstance(thermo0, NASA):
                thermo = thermo0
                if thermo.E0 is None:
                    thermo.E0 = wilhoit.E0
            else:
                thermo = wilhoit.to_nasa(Tmin=100.0, Tmax=5000.0, Tint=1000.0)
        else:
            # gas phase with species matching thermo library keep the NASA from library or convert if group additivity
            if "Thermo library" in thermo0.comment and isinstance(thermo0, NASA):
                thermo = thermo0
                if thermo.E0 is None:
                    thermo.E0 = wilhoit.E0
            else:
                thermo = wilhoit.to_nasa(Tmin=100.0, Tmax=5000.0, Tint=1000.0)
    else:
        raise Exception('thermo_class neither NASA nor Wilhoit.  Cannot process thermo data.')

    return thermo


def generate_thermo_data(spc, thermo_class=NASA, solvent_name=''):
    """
    Generates thermo data, first checking Libraries, then using either QM or Database.
    
    The database generates the thermo data for each structure (resonance isomer),
    picks that with lowest H298 value.
    
    It then calls :meth:`process_thermo_data`, to convert (via Wilhoit) to NASA
    and set the E0.
    
    Result stored in `spc.thermo` and returned.
    """

    try:
        thermodb = get_db('thermo')
        if not thermodb: raise Exception
    except Exception:
        logging.debug('Could not obtain the thermo database. Not generating thermo...')
        return None

    thermo0 = thermodb.get_thermo_data(spc)

    # 1. maybe only submit cyclic core
    # 2. to help radical prediction, HBI should also
    #    look up centrailThermoDB for its saturated version
    #    currently it only looks up libraries or estimates via GAV 
    from rmgpy.rmg.input import get_input

    try:
        thermo_central_database = get_input('thermo_central_database')
    except Exception:
        logging.debug('thermoCentralDatabase could not be found.')
        thermo_central_database = None

    if thermo_central_database and thermo_central_database.client \
            and thermo_central_database.satisfy_registration_requirements(spc, thermo0, thermodb):
        thermo_central_database.register_in_central_thermo_db(spc)

    return process_thermo_data(spc, thermo0, thermo_class, solvent_name)


def evaluator(spc, solvent_name=''):
    """
    Module-level function passed to workers.

    generates a thermodata object for the 
    identifier and stores it in the 
    thermo database.

    Next, thermo is generated of this species.

    """
    logging.debug("Evaluating spc %s ", spc)

    if not isinstance(spc.molecule[0], Fragment):
        spc.generate_resonance_structures()
        thermo = generate_thermo_data(spc, solvent_name=solvent_name)
    else:
        # assume it's a species for Fragment
        spc.molecule[0].assign_representative_species()
        spc_repr = spc.molecule[0].species_repr
        spc_repr.generate_resonance_structures()
        thermo = generate_thermo_data(spc_repr, solvent_name=solvent_name)

    return thermo


def submit(spc, solvent_name=''):
    """
    Submits a request to calculate chemical data for the Species object.

    In a parallel run, the thermo attribute will
    store the future object, until the get method
    is called, which replaces the future object with 
    the result.

    """
    spc.thermo = evaluator(spc, solvent_name=solvent_name)
