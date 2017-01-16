
import numpy
import math

import logging as logging
from rmgpy.scoop_framework.util import submit_
from rmgpy.data.rmg import getDB
import rmgpy.constants as constants
from rmgpy.statmech import Conformer
from rmgpy.thermo import Wilhoit, NASA, ThermoData

def processThermoData(spc, thermo0, thermoClass=NASA):
    """
    Converts via Wilhoit into required `thermoClass` and sets `E0`.
    
    Resulting thermo is returned.
    """
    # TODO moving this as a global import leads to circular imports.
    from rmgpy.rmg.main import solvent

    thermo = None

    # Always convert to Wilhoit so we can compute E0
    if isinstance(thermo0, Wilhoit):
        wilhoit = thermo0
    elif isinstance(thermo0, ThermoData):
        wilhoit = thermo0.toWilhoit(B=1000.)
    else:
        wilhoit = thermo0.toWilhoit()

    # Add on solvation correction
    if solvent and not "Liquid thermo library" in thermo0.comment:
        solvationDatabase = getDB('solvation')
        #logging.info("Making solvent correction for {0}".format(solvent.solventName))
        soluteData = solvationDatabase.getSoluteData(spc)
        solvationCorrection = solvationDatabase.getSolvationCorrection(soluteData, solvent.solventData)
        # correction is added to the entropy and enthalpy
        wilhoit.S0.value_si = (wilhoit.S0.value_si + solvationCorrection.entropy)
        wilhoit.H0.value_si = (wilhoit.H0.value_si + solvationCorrection.enthalpy)
        
    # Compute E0 by extrapolation to 0 K
    if spc.conformer is None:
        spc.conformer = Conformer()
    spc.conformer.E0 = wilhoit.E0
    
    # Convert to desired thermo class
    if thermoClass is Wilhoit:
        thermo = wilhoit
    elif thermoClass is NASA:
        if solvent:
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
    

def generateThermoData(spc, thermoClass=NASA):
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
    except Exception, e:
        logging.debug('Could not obtain the thermo database. Not generating thermo...')
        return None
    
    thermo0 = thermodb.getThermoData(spc) 
        
    return processThermoData(spc, thermo0, thermoClass)    


def evaluator(spc):
    """
    Module-level function passed to workers.

    generates a thermodata object for the 
    identifier and stores it in the 
    thermo database.

    Next, thermo is generated of this species.

    """
    logging.debug("Evaluating spc %s ", spc)

    spc.generateResonanceIsomers()
    thermo = generateThermoData(spc)

    return thermo

def submit(spc):
    """
    Submits a request to calculate chemical data for the Species object.

    In a parallel run, the thermo attribute will
    store the future object, until the get method
    is called, which replaces the future object with 
    the result.

    """
    spc.thermo = submit_(evaluator, spc)