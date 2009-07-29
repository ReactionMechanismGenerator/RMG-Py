#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

import rmg
import math
import logging
import sys

# get quantities package from http://pypi.python.org/pypi/quantities
# read about it at http://packages.python.org/quantities/index.html
import quantities as pq
pq.UnitQuantity('kilocalories', pq.cal*1e3, symbol='kcal')
pq.UnitQuantity('kilojoules', pq.J*1e3, symbol='kJ')
pq.UnitQuantity('kilomoles', pq.mol*1e3, symbol='kmol')


# First we import some Cantera ctml_writer classes 
# (which do some checks and allow you to write 'ctml' XML)
# then we can modify the classes we want to do more
sys.path.append('../external')
import ctml_writer as cti
from ctml_writer import units, OneAtm

_uConc = None # set these once base units have been set by 
_uTime = None # importing .cti file inside loadCanteraFile()

if cti._species:
    logging.info("Resetting cantera memory")
    reload(cti)
_species=cti._species
_reactions=cti._reactions
_speciesnames=cti._speciesnames
_speciesByName={} # a dictionary to look up species objects by their cantera name

# Classes should usually be named with CamelCase but need to stick with 
# Cantera .cti names, hence lowercase names here:

class ideal_gas(cti.ideal_gas):
    pass # for now, don't modify this class, but must define it

class state(cti.state):
    pass # for now, don't modify this class, but must define it 
    
class species(cti.species):
    """A species derived from ctml_writer.species.
    
    species(name,atoms,note,thermo,transport,charge)
    """
    def __init__(self, *args, **kwargs):
        # first call the parent __init__
        cti.species.__init__(self, *args, **kwargs)
        # now put it in the _speciesByName dictionary
        name=self._name
        global _speciesByName
        _speciesByName[name]=self
        # if overloading this __init__ to fill in the _speciesByName dictionary
        # does not work well, then replace uses of _speciesByName dictionary
        # with calls to the function getSpeciesByName()
        
    def getRmgSpecies(self):
        """Fetch the RMG Species object associated with this cantera species"""
        if not hasattr(self,'_RmgSpecies'):
            self.makeRmgSpecies()
        return self._RmgSpecies
        
    def makeRmgSpecies(self):
        """Create an RMG Species object based on this cantera species"""
        name=self._name
        id=cti._speciesnames.index(name)
        self._RmgSpecies = rmg.species.Species(id=id, label=name)
        # initialise thermo now
        self._RmgSpecies.thermoData = rmg.species.ThermoNASAData(
            comment="Imported from Cantera file" )
        for canterapoly in self._thermo:
            poly = canterapoly.getRmgPolynomial()
            self._RmgSpecies.thermoData.addPolynomial(poly)
            
    def getCorrectNasaPoly(self,T):
        """Select the first NASA polynomial that is valid at temperature T"""
        for nasapoly in self._thermo:
            if nasapoly.ValidTemperature(T): return nasapoly
        return None
                    
class NASA(cti.NASA):
    """NASA polynomial representation of thermo. For importing Cantera files."""
    def getRmgPolynomial(self):
        """Fetch the RMG ThermoNASAPolynomial associated with this Cantera NASA"""
        if not hasattr(self,'_RmgThermoNASAPolynomial'):
            self.makeRmgPolynomial()
        return self._RmgThermoNASAPolynomial
    def makeRmgPolynomial(self):
        """Make the RMG ThermoNASAPolynomial associated with this Cantera NASA"""
        if self._pref > 0.0: raise Exception(
            "Not sure what to do with customised standard state pressure")
        poly = rmg.species.ThermoNASAPolynomial( T_range = self._t, 
                                  coeffs = self._coeffs,
                                  comment = "Imported from Cantera file")
        self._RmgThermoNASAPolynomial = poly
        
class reaction(cti.reaction):
    """A chemical reaction. For importing Cantera files"""
    def getRmgReaction(self):
        """Fetch the RMG Reaction object associated with this cantera reaction"""
        if not hasattr(self,'_RmgReaction'):
            self.makeRmgReaction()
        return self._RmgReaction
    def makeRmgReaction(self):
        """Create an RMG Reaction object based on this cantera reaction"""        
        reactants=[]
        for speciesName,order in self._r.items(): #  reactants
            assert order==int(order), "RMG currently needs elementary reactions"
            order = int(order)
            species=_speciesByName[speciesName].getRmgSpecies()
            reactants.extend([species]*order)
        products=[]
        for speciesName,order in self._p.items(): #  reactants
            assert order==int(order), "RMG currently needs elementary reactions"
            order = int(order)
            species=_speciesByName[speciesName].getRmgSpecies()
            products.extend([species]*order)
        rxn = rmg.reaction.Reaction(reactants=reactants, products=products)
        rxn.kinetics = [self.makeRmgKinetics()]
        self._RmgReaction = rxn
    def makeRmgKinetics(self):
        """Create and return an RMG ArrheniusKinetics object based on self._kf"""
        A = self._kf[0] / _uConc**(self.getReactantNu()-1) / _uTime
        A = float(A.simplified)
        n = self._kf[1]
        #assert cti._ue=='kcal/mol', 'unit conversion not implemented'
        Ea=float(pq.Quantity(self._kf[2], cti._ue).simplified)
        kinetics=rmg.reaction.ArrheniusKinetics(A=A, Ea=Ea, n=n)
        kinetics.comment="Generated from Cantera expression %s"%(self._kf)
#        kinetics.Trange = [0.0, 100000.0]
        return kinetics
    def getReactantNu(self):
        """Returns the stoichiometry in the forwards direction.
        
        i.e. the number of reactant molecules.
        uses self._reactantNu to cache the answer"""
        if hasattr(self,'_reactantNu'): return self._reactantNu
        reactantNu=0
        for speciesName,order in self._r.items(): # add reactants
            reactantNu += order
        self._reactantNu=reactantNu
        return reactantNu
    def getProductNu(self):
        """Returns the stoichiometry in the reverse direction. 
        
        i.e. the number of product molecules.
        uses self._productNu to cache the answer"""
        if hasattr(self,'_productNu'): return self._productNu
        productNu=0
        for speciesName,order in self._p.items(): # add products
            productNu += order
        self._productNu=productNu
        return productNu
    def getDeltaNu(self):
        """Returns the change in stoichiometry of the reaction, delta Nu.
        
        deltaNu= productNu - reactantNu
        uses self._deltaNu to cache the answer"""
        if hasattr(self,'_deltaNu'): return self._deltaNu
        self._deltaNu=self.getProductNu()-self.getReactantNu()
        return self._deltaNu        

def getSpeciesByName(name):
    """Select a species by its name."""
    for s in _species:
        if s._name==name:
            return s
    return None

def loadCanteraFile(filepath):
    """Load a Cantera input (.cti) file. Returns a CoreEdgeReactionModel"""
    import os
    logging.info("Loading Cantera file %s"%filepath)
    logging.debug("which is located at %s"%os.path.abspath(filepath))
    base = os.path.basename(filepath)
    root, ext = os.path.splitext(base)
    cti.dataset(root)
    assert not _species, "Oops! There's already cantera species info in memory"
    try:
        execfile(filepath)
    except NameError, error: # python 2.5 syntax
        logging.exception("RMG's cantera interpreter can't read this cti file."
            + error.message.lstrip('name '))
        raise    
    # the cti file will have called units() which may have changed _umol etc
    global _uConc
    global _uTime
    _uConc=pq.Quantity(1,cti._umol) / pq.Quantity(1,cti._ulen)**3
    _uTime=pq.Quantity(1,cti._utime)

    model = rmg.model.CoreEdgeReactionModel()
    for s in _species:
        model.addSpeciesToCore( s.getRmgSpecies() )
    for r in _reactions:
        model.addReactionToCore( r.getRmgReaction() )
    return model

def loadChemkinFile(filepath='chem.inp', thermodb='',trandb=''):
    """Loads a Chemkin file. (Converts to Cantera first)"""
    from Cantera import ck2cti
    import shutil
    import os
    oldpath = os.path.abspath(filepath)
    oldwd   = os.getcwd()
    filename = os.path.basename(filepath)
    newpath = os.path.join(oldwd,rmg.constants.scratchDir, filename)
    
    if not os.path.exists(newpath) or not os.path.samefile(oldpath,newpath):
        logging.debug("Copying %s to %s"%(oldpath,newpath))
        shutil.copy2(oldpath, newpath)
    ctifile = os.path.splitext(newpath)[0]+'.cti'
    if os.path.exists(ctifile):
        #move aside old file
        os.rename(ctifile,'%s-%d.cti'%(os.path.splitext(ctifile)[0],
                                         os.lstat(ctifile).st_mtime) )
    nm='chem'
    # convert from chemkin to cti
    try:
        os.chdir(rmg.constants.scratchDir)
        ck2cti.ck2cti(infile = newpath, thermodb = thermodb,
                  trandb = trandb, idtag = nm, debug=0, validate=1)
    except: 
        logging.error("Conversion from Chemkin to Cantera failed:")
        for line in open('ck2cti.log'):
             logging.error(line.rstrip())
        os.chdir(oldwd) # change back to where you were
        raise
    logging.info('Converting %s from Chemkin to Cantera:'%filepath)
    for line in open('ck2cti.log'):
        logging.debug(line.rstrip())
    os.chdir(oldwd) # change back to where you were
    model = loadCanteraFile(ctifile)
    return model


###### END OF MODULE


if __name__ == "__main__":
    #reload(cti)
    import sys, os

    rmg.initializeLog(verbose=10)
    
    if len(sys.argv)>1:
        filename = sys.argv[1]
    else:
        logging.debug("No file specificed; using default, chem.cti")
        filename = 'chem.cti' # default input file
    if os.path.splitext(filename)[1]=='.cti':
        model = loadCanteraFile(filename)
    else:
        model = loadChemkinFile(filename)
    

    
