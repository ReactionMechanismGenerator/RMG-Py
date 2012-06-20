"""
Created on Apr 29, 2012

@author: nmvdewie
"""
import os
import math

import rmgpy.quantity
from rmgpy.thermo import ThermoData
from rmgpy.statmech import RigidRotor, HarmonicOscillator, Translation, StatesModel

import symmetry
 
class TDPropertiesCalculator:
    
    def __init__(self, molfile, qmdata, environ = os.environ.get("RMG.workingDirectory"), pointGroup = None):
        self.qmdata = qmdata
               
        self.molfile = molfile
        
        self.environ = environ   
        
        self.determinePointGroup()
        
        trans = Translation(mass=(qmdata.molecularMass,"amu"))
        rot = RigidRotor(linear=self.pointGroup.linear, inertia=self.qmdata.rotationalConstants, symmetry=self.pointGroup.symmetryNumber)
        vib = HarmonicOscillator(frequencies=qmdata.frequencies)
        self.statesmodel = StatesModel(modes=[trans, rot, vib], spinMultiplicity=self.qmdata.groundStateDegeneracy)
     
        
    def determinePointGroup(self):
        #determine point group using the SYMMETRY Program

        pgc = symmetry.PointGroupCalculator(self.molfile, self.qmdata, self.environ);
        self.pointGroup = pgc.calculate();
    
    def calculateChiralityCorrection(self):
        if self.pointGroup.chiral:
            return rmgpy.quantity.constants.R * math.log(2);
        else:
            return 0.

    def calculate(self):
        #we will use number of atoms from above (alternatively, we could use the chemGraph); this is needed to test whether the species is monoatomic
        Hartree_to_kcal = 627.5095#conversion from Hartree to kcal/mol taken from Gaussian thermo white paper
        Hf298 = self.qmdata.energy * Hartree_to_kcal;

        S298 = self.statesmodel.getEntropy(298.0)
        
        Tdata = [300.0,400.0,500.0,600.0,800.0,1000.0,1500.0]
        Cp = []
        for T in Tdata:
            Cp.append(self.statesmodel.getHeatCapacity(T))

        S298 = S298 + self.calculateChiralityCorrection()

        Tmin = 300.0
        
        Tmax = 2000.0
        
        comment = "PM3 or MM4 calculation"
        
        return ThermoData(Tdata, Cp, Hf298, S298, Tmin, Tmax, comment)
