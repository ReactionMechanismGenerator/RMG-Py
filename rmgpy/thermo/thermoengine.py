#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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

from __future__ import division

import logging as logging
import math

import numpy as np

import rmgpy.constants as constants
from rmgpy.data.rmg import get_db
from rmgpy.statmech import Conformer
from rmgpy.thermo import Wilhoit, NASA, ThermoData


def processThermoData(spc, thermo0, thermoClass=NASA, solventName=''):
    """
    Converts via Wilhoit into required `thermoClass` and sets `E0`.
    
    Resulting thermo is returned.
    """
    thermo = None

    # Always convert to Wilhoit so we can compute E0
    if isinstance(thermo0, Wilhoit):
        wilhoit = thermo0
    elif isinstance(thermo0, ThermoData):
        wilhoit = thermo0.toWilhoit(B=1000.)
    else:
        wilhoit = thermo0.toWilhoit()

    # Add on solvation correction
    solvation_database = get_db('solvation')
    if not solventName or solvation_database is None:
        logging.debug('Solvent database or solventName not found. Solvent effect was not utilized')
        solvent_data = None
    else:
        solvent_data = solvation_database.get_solvent_data(solventName)
    if solvent_data and not "Liquid thermo library" in thermo0.comment:
        solvation_database = get_db('solvation')
        soluteData = solvation_database.get_solute_data(spc)
        solvation_correction = solvation_database.get_solvation_correction(soluteData, solvent_data)
        # correction is added to the entropy and enthalpy
        wilhoit.S0.value_si = (wilhoit.S0.value_si + solvation_correction.entropy)
        wilhoit.H0.value_si = (wilhoit.H0.value_si + solvation_correction.enthalpy)

    # Compute E0 by extrapolation to 0 K
    if spc.conformer is None:
        spc.conformer = Conformer()
    spc.conformer.e0 = wilhoit.E0

    # Convert to desired thermo class
    if thermoClass is Wilhoit:
        thermo = wilhoit
    elif thermoClass is NASA:
        if solvent_data:
            # If liquid phase simulation keep the nasa polynomial if it comes from a liquid phase thermoLibrary.
            # Otherwise convert wilhoit to NASA
            if "Liquid thermo library" in thermo0.comment and isinstance(thermo0, NASA):
                thermo = thermo0
                if thermo.E0 is None:
                    thermo.E0 = wilhoit.E0
            else:
                thermo = wilhoit.toNASA(Tmin=100.0, Tmax=5000.0, Tint=1000.0)
        else:
            # gas phase with species matching thermo library keep the NASA from library or convert if group additivity
            if "Thermo library" in thermo0.comment and isinstance(thermo0, NASA):
                thermo = thermo0
                if thermo.E0 is None:
                    thermo.E0 = wilhoit.E0
            else:
                thermo = wilhoit.toNASA(Tmin=100.0, Tmax=5000.0, Tint=1000.0)
    else:
        raise Exception('thermoClass neither NASA nor Wilhoit.  Cannot process thermo data.')

    if thermo.__class__ != thermo0.__class__:
        # Compute RMS error of overall transformation
        Tlist = np.array([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], np.float64)
        err = 0.0
        for T in Tlist:
            err += (thermo.get_heat_capacity(T) - thermo0.get_heat_capacity(T)) ** 2
        err = math.sqrt(err / len(Tlist)) / constants.R
        # logging.log(logging.WARNING if err > 0.2 else 0, 'Average RMS error in heat capacity fit to {0} = {1:g}*R'.format(spc, err))

    return thermo


def generateThermoData(spc, thermoClass=NASA, solventName=''):
    """
    Generates thermo data, first checking Libraries, then using either QM or Database.
    
    The database generates the thermo data for each structure (resonance isomer),
    picks that with lowest H298 value.
    
    It then calls :meth:`processThermoData`, to convert (via Wilhoit) to NASA
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
    from rmgpy.rmg.input import getInput

    try:
        thermo_central_database = getInput('thermoCentralDatabase')
    except Exception:
        logging.debug('thermoCentralDatabase could not be found.')
        thermo_central_database = None

    if thermo_central_database and thermo_central_database.client \
            and thermo_central_database.satisfyRegistrationRequirements(spc, thermo0, thermodb):
        thermo_central_database.registerInCentralThermoDB(spc)

    return processThermoData(spc, thermo0, thermoClass, solventName)


def evaluator(spc, solventName=''):
    """
    Module-level function passed to workers.

    generates a thermodata object for the 
    identifier and stores it in the 
    thermo database.

    Next, thermo is generated of this species.

    """
    logging.debug("Evaluating spc %s ", spc)

    spc.generate_resonance_structures()
    thermo = generateThermoData(spc, solventName=solventName)

    return thermo


def submit(spc, solventName=''):
    """
    Submits a request to calculate chemical data for the Species object.

    In a parallel run, the thermo attribute will
    store the future object, until the get method
    is called, which replaces the future object with 
    the result.

    """
    spc.thermo = evaluator(spc, solventName=solventName)
