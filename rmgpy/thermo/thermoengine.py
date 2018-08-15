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

import numpy
import math

import logging as logging
from rmgpy.data.rmg import getDB
import rmgpy.constants as constants
from rmgpy.statmech import Conformer
from rmgpy.thermo import Wilhoit, NASA, ThermoData
import rmgpy.data.rmg

def processThermoData(spc, thermo0, thermoClass=NASA, solventName = ''):
    """
    Converts via Wilhoit into required `thermoClass` and sets `E0`.
    
    Resulting thermo is returned.
    """
    # TODO moving this as a global import leads to circular imports.
    from rmgpy.rmg.model import Species

    thermo = None

    # Always convert to Wilhoit so we can compute E0
    if isinstance(thermo0, Wilhoit):
        wilhoit = thermo0
    elif isinstance(thermo0, ThermoData):
        wilhoit = thermo0.toWilhoit(B=1000.)
    else:
        wilhoit = thermo0.toWilhoit()

    # Add on solvation correction
    solvationdatabase = getDB('solvation')
    if not solventName or solvationdatabase is None:
        logging.debug('Solvent database or solventName not found. Solvent effect was not utilized')
        solventData = None
    else:
        solventData = solvationdatabase.getSolventData(solventName)
    if solventData and not "Liquid thermo library" in thermo0.comment:
        solvationdatabase = getDB('solvation')
        #logging.info("Making solvent correction for {0}".format(Species.solventName))
        soluteData = solvationdatabase.getSoluteData(spc)
        solvation_correction = solvationdatabase.getSolvationCorrection(soluteData, solventData)
        # correction is added to the entropy and enthalpy
        wilhoit.S0.value_si = (wilhoit.S0.value_si + solvation_correction.entropy)
        wilhoit.H0.value_si = (wilhoit.H0.value_si + solvation_correction.enthalpy)
        
    # Compute E0 by extrapolation to 0 K
    if spc.conformer is None:
        spc.conformer = Conformer()
    spc.conformer.E0 = wilhoit.E0
    
    # Convert to desired thermo class
    if thermoClass is Wilhoit:
        thermo = wilhoit
    elif thermoClass is NASA:
        if solventData:
            #if liquid phase simulation keep the nasa polynomial if it comes from a liquid phase thermoLibrary. Otherwise convert wilhoit to NASA
            if "Liquid thermo library" in thermo0.comment and isinstance(thermo0, NASA):
                thermo = thermo0
                if thermo.E0 is None:
                    thermo.E0 = wilhoit.E0
            else:
                thermo = wilhoit.toNASA(Tmin=100.0, Tmax=5000.0, Tint=1000.0)
        else: 
            #gas phase with species matching thermo library keep the NASA from library or convert if group additivity
            if "Thermo library" in thermo0.comment and isinstance(thermo0,NASA):
                thermo=thermo0
                if thermo.E0 is None:
                    thermo.E0 = wilhoit.E0
            else:
                thermo = wilhoit.toNASA(Tmin=100.0, Tmax=5000.0, Tint=1000.0)
    else:
        raise Exception('thermoClass neither NASA nor Wilhoit.  Cannot process thermo data.')
    
    if thermo.__class__ != thermo0.__class__:
        # Compute RMS error of overall transformation
        Tlist = numpy.array([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], numpy.float64)
        err = 0.0
        for T in Tlist:
            err += (thermo.getHeatCapacity(T) - thermo0.getHeatCapacity(T))**2
        err = math.sqrt(err/len(Tlist))/constants.R
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
        thermodb = getDB('thermo')
        if not thermodb: raise Exception
    except Exception:
        logging.debug('Could not obtain the thermo database. Not generating thermo...')
        return None

    thermo0 = thermodb.getThermoData(spc)

    # 1. maybe only submit cyclic core
    # 2. to help radical prediction, HBI should also
    #    look up centrailThermoDB for its saturated version
    #    currently it only looks up libraries or estimates via GAV 
    from rmgpy.rmg.input import getInput
    
    try:
        thermoCentralDatabase = getInput('thermoCentralDatabase')
    except Exception:
        logging.debug('thermoCentralDatabase could not be found.')
        thermoCentralDatabase = None
    
    if thermoCentralDatabase and thermoCentralDatabase.client \
        and thermoCentralDatabase.satisfyRegistrationRequirements(spc, thermo0, thermodb):
        
        thermoCentralDatabase.registerInCentralThermoDB(spc)
        
    return processThermoData(spc, thermo0, thermoClass, solventName)


def evaluator(spc, solventName = ''):
    """
    Module-level function passed to workers.

    generates a thermodata object for the 
    identifier and stores it in the 
    thermo database.

    Next, thermo is generated of this species.

    """
    logging.debug("Evaluating spc %s ", spc)

    if isinstance(spc.molecule[0], Molecule):
        spc.generate_resonance_structures()
        thermo = generateThermoData(spc, solventName=solventName)
    else:
        # assume it's a species for Fragment
        spc.molecule[0].assign_representative_species()
        spc_repr = spc.molecule[0].species_repr
        spc_repr.generate_resonance_structures()
        thermo = generateThermoData(spc_repr, solventName=solventName)

    return thermo

def submit(spc, solventName = ''):
    """
    Submits a request to calculate chemical data for the Species object.

    In a parallel run, the thermo attribute will
    store the future object, until the get method
    is called, which replaces the future object with 
    the result.

    """
    spc.thermo = evaluator(spc, solventName= solventName)

