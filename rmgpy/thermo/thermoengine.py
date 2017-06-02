
import numpy
import math

import logging as logging
from rmgpy.scoop_framework.util import submit_
from rmgpy.data.rmg import getDB
import rmgpy.constants as constants
from rmgpy.statmech import Conformer
from rmgpy.thermo import Wilhoit, NASA, ThermoData

def processThermoData(spc, thermo0, solventThermoEstimator, thermoClass=NASA):
    """
    Converts via Wilhoit into required `thermoClass` and sets `E0`.
    
    Resulting thermo is returned.
    """
    # TODO moving this as a global import leads to circular imports.
    # Because thermoEstimator.py does not use rmgpy.rmg.main, solvent data must be directly passed on
    # from thermoEstimator input file. If this is normal rmg.py run, it will get the solvent data from
    # rmgpy.rmg.main.
    if solventThermoEstimator is None:
        from rmgpy.rmg.main import solvent
    else:
        solvent = solventThermoEstimator

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
        if solvent.solventData.inCoolProp:
            if spc.label == solvent.solventName: # the species is the solvent
                wilhoit, nasa = solvationDatabase.getSolventThermo(solvent.solventData, wilhoit)
            else: # the species is the solute
                wilhoit, nasa = solvationDatabase.getSolvationThermo(soluteData, solvent.solventData, wilhoit)
        else:
            solvationCorrection = solvationDatabase.getSolvationCorrection298(soluteData, solvent.solventData)
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
                if solvent.solventData.inCoolProp:
                    # Tmax is set high to prevent value error when Cp is calculated for chemkin output file
                    # Tmax should actually be the critical temperature of the solvent, Tc.
                    # RMG requires NASA to have two polynomials while liquid phase thermo is fitted to one NASA polynomial
                    # To prevent error, two NASA polynomials with the same coefficients are used
                    # Liquid phase thermo should actually have one NASA polynomial with Tmin = 298K and Tmax = Tc
                    thermo = wilhoit.toNASA(Tmin=298., Tmax=5000.0, Tint=1000.0)
                    thermo.poly1.coeffs = nasa.poly1.coeffs
                    thermo.poly1.coeffs = nasa.poly1.coeffs
                else:
                    thermo = wilhoit.toNASA(Tmin=298., Tmax=5000.0, Tint=1000.0)
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
    

def generateThermoData(spc, solventThermoEstimator, thermoClass=NASA):
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

    # 1. maybe only submit cyclic core
    # 2. to help radical prediction, HBI should also
    #    look up centrailThermoDB for its saturated version
    #    currently it only looks up libraries or estimates via GAV 
    from rmgpy.rmg.input import getInput
    
    try:
        thermoCentralDatabase = getInput('thermoCentralDatabase')
    except Exception, e:
        logging.debug('thermoCentralDatabase could not be found.')
        thermoCentralDatabase = None
    
    if thermoCentralDatabase and thermoCentralDatabase.client \
        and thermoCentralDatabase.satisfyRegistrationRequirements(spc, thermo0, thermodb):
        
        thermoCentralDatabase.registerInCentralThermoDB(spc)
        
    return processThermoData(spc, thermo0, solventThermoEstimator, thermoClass)


def evaluator(spc, solventThermoEstimator):
    """
    Module-level function passed to workers.

    generates a thermodata object for the 
    identifier and stores it in the 
    thermo database.

    Next, thermo is generated of this species.

    """
    logging.debug("Evaluating spc %s ", spc)

    spc.generateResonanceIsomers()
    thermo = generateThermoData(spc, solventThermoEstimator)

    return thermo

def submit(spc, solventThermoEstimator=None):
    """
    Submits a request to calculate chemical data for the Species object.

    In a parallel run, the thermo attribute will
    store the future object, until the get method
    is called, which replaces the future object with 
    the result.

    """
    spc.thermo = submit_(evaluator, spc, solventThermoEstimator)