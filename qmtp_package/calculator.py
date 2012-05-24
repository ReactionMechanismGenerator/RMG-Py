'''
Created on Apr 29, 2012

@author: nmvdewie
'''


import symmetry as symm
from rmgpy.quantity import constants as cnts, Quantity
import os
import math
from rmgpy.thermo import ThermoData

class TDPropertiesCalculator:
    
    def __init__(self, name, dir, qmdata, environ = os.environ.get("RMG.workingDirectory"), pointGroup = None):
        self.qmdata = qmdata
        
        self.name = name
        
        self.dir = dir
        
        self.environ = environ
        

        self.determinePointGroup()
        
        
    def determinePointGroup(self):
        #determine point group using the SYMMETRY Program

        pgc = symm.PointGroupCalculator(self.name, self.dir, self.qmdata, self.environ);
        self.pointGroup = pgc.calculate();
    
   
    def calculateSymmetryCorrection(self):
        '''
     * gets the statistical correction for S in dimensionless units (divided by R)

        //determine statistical correction factor for 1. external rotational symmetry (affects rotational partition function) and 2. chirality (will add R*ln2 to entropy) based on point group
        //ref: http://cccbdb.nist.gov/thermo.asp
        //assumptions below for Sn, T, Th, O, I seem to be in line with expectations based on order reported at: http://en.wikipedia.org/w/index.php?title=List_of_character_tables_for_chemically_important_3D_point_groups&oldid=287261611 (assuming order = symmetry number * 2 (/2 if chiral))...this appears to be true for all point groups I "know" to be correct
        //minor concern: does SYMMETRY appropriately calculate all Sn groups considering 2007 discovery of previous errors in character tables (cf. Wikipedia article above)
        '''
        return cnts.R * math.log(self.pointGroup.symmetryNumber);   
    
    def calculateChiralityCorrection(self):
        if self.pointGroup.chiral:
            return cnts.R * math.log(2);
        
        else:
            return 0.
    
    
    def calcElecS(self):
        '''
        * electronic partition function
        '''
        return cnts.R * math.log(self.qmdata.groundStateDegeneracy)
    
    def calcVibS(self,temperature):
        '''
        //gmagoon 6/8/09
    //calculate the vibrational contribution (divided by R, dimensionless) at temperature, T, in Kelvin to entropy
    //p_freqs in cm^-1; c in cm/s; k in J/K; h in J-s
    //ref.: http://cccbdb.nist.gov/thermo.asp
        '''
        Scontrib = 0.;
        dr = 0.
        for freq in self.qmdata.frequencies:
            dr = cnts.h * cnts.c * freq / (cnts.kB * temperature)#frequently used dimensionless ratio
            Scontrib = Scontrib - math.log(1.-math.exp(-dr)) + dr * math.exp(-dr) / (1.-math.exp(-dr))
        

        return cnts.R*Scontrib;
    
    
    
    def calcVibH(self,temperature):
        '''
        //gmagoon 6/23/10
    //calculate the vibrational contribution (divided by R, units of K) at temperature, T, in Kelvin to Hthermal (ZPE not included)
    //p_freqs in cm^-1; c in cm/s; k in J/K; h in J-s
    //we need to ignore zero frequencies, as MM4 does, to avoid dividing by zero; based on equation for Ev in http://www.gaussian.com/g_whitepap/thermo.htm; however, we ignore zero point contribution to be consistent with the input that CanTherm currently takes; note, however, that this could introduce some small inaccuracies, as the frequencies may be slightly different in CanTherm vs. MM4, particularly for a frequency < 7.7 cm^-1 (treated as zero in MM4)
        '''
        Hcontrib = 0.;
        dr = 0.;
        for freq in self.qmdata.frequencies:
            if freq > 0.0:#ignore zero frequencies as MM4 does
                dr = cnts.h*cnts.c * freq / (cnts.kB * temperature)#frequently used dimensionless ratio
                Hcontrib = Hcontrib + dr*temperature/(math.exp(dr) - 1.);
        return Hcontrib;

    def calcTransS(self,temperature):
        '''
        Translational contribution to the molar entropy at standard conditions (1 bar)
        '''
        return cnts.R * (3./2.*math.log(2.*cnts.pi * self.qmdata.molecularMass
                /(1000.*cnts.Na*math.pow(cnts.h,2.)))
                +5./2.*math.log(cnts.kB*temperature)
                -math.log(100000.)+5./2.)
       
    def calcTransCp(self):
        '''
        translational contribution to standard heat capacity at constant pressure: 
        '''
        return 5./2.*cnts.R;
    

    def calcRotCp(self):
        if self.pointGroup.linear:
            return cnts.R
        
        else:
            return 3./2. * cnts.R;
      
    def calcVibCp(self, temperature):
        '''
          //gmagoon 6/8/09
    //calculate the vibrational contribution (divided by R, dimensionless) at temperature, T, in Kelvin to heat capacity, Cp
    //p_freqs in cm^-1; c in cm/s; k in J/K; h in J-s
    //ref.: http://cccbdb.nist.gov/thermo.asp
        '''
        Cpcontrib = 0
        for freq in self.qmdata.frequencies:
            dr = cnts.h * cnts.c * freq / (cnts.kB * temperature)#frequently used dimensionless ratio
            Cpcontrib = Cpcontrib + math.pow(dr, 2.) * math.exp(-dr) / math.pow(1.-math.exp(-dr),2.)
        
            
        return Cpcontrib;
    
    def calcRotS(self, temperature):
        '''
        * rotational contribution to entropy
     * 
     * for linear molecules, one rotational constant will be zero, 
     * the other two will be identical
        '''

        if self.pointGroup.linear:
            if self.qmdata.rotationalConstants[0] > 0.0001:
                rotCons = self.qmdata.rotationalConstants[0]
            else:
                rotCons = self.qmdata.rotationalConstants[1]

            return cnts.R*(math.log(cnts.kB*temperature/(cnts.h*rotCons))+1);
        
        else:
            return cnts.R*(3./2.*math.log(cnts.kB*temperature/cnts.h)-1./2.*math.log(self.qmdata.rotationalConstants[0]*self.qmdata.rotationalConstants[1]*self.qmdata.rotationalConstants[2]/cnts.pi)+3./2.)+cnts.R*self.calcVibS(temperature);
        

    def calculate(self):
        '''
        //calculate thermo quantities using stat. mech. equations
        '''

        #we will use number of atoms from above (alternatively, we could use the chemGraph); this is needed to test whether the species is monoatomic
        Hartree_to_kcal = 627.5095#conversion from Hartree to kcal/mol taken from Gaussian thermo white paper
        Hf298 = self.qmdata.energy * Hartree_to_kcal;

        '''
         * electronic + translation; note use of 10^5 Pa for standard pressure; 
         * also note that molecular mass needs to be divided by 1000 for kg units
        '''
        S298 = self.calcElecS() + self.calcTransS(298.0);
        Tdata = [300.0,400.0,500.0,600.0,800.0,1000.0,1500.0]
        Cp = [0,0,0,0,0,0,0,]
        Cp = [t + self.calcTransCp() for t in Cp]

        '''
         * //include statistical correction and rotational
         *  (without symmetry number, vibrational contributions if species is polyatomic
        '''
        if self.qmdata.numberOfAtoms > 1:
            Cp = [t + self.calcRotCp() for t in Cp]
            for i in range(len(Tdata)):
                Cp[i] = Cp[i] +  self.calcVibCp(Tdata[i])

            S298 = S298 
            -self.calculateSymmetryCorrection()
            +self.calculateChiralityCorrection()
            +self.calcRotS(298.15)
            +self.calcVibS(298.15);
        Cpdata = Cp
        S298data = S298
        Hf298data = Hf298
        Tmin = 300.0
        Tmax = 2000.0
        comment = "PM3 or MM4 calculation"
        
        return ThermoData(Tdata, Cpdata, Hf298data, S298data,Tmin, Tmax, comment)
    