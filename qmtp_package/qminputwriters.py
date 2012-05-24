'''
Module that writes input files for the various QM packages.

Input files contain information on the 3D coordinates of the molecule, but also 
command strings to specify the level of theory, and other flags.



'''



import qmtp
import platform
from subprocess import Popen
from symmetry import *
import logging
import os
import math



class QMInputWriter:
    '''
     Supertype for all input writers for quantum chemistry methods
    '''
    
    #the number of keyword permutations available update as additional options are added
    scriptAttempts = 0
    
    #we will try a second time with crude coordinates if the UFF refined coordinates do not work
    maxAttemptNumber = 0

        
    def __init__(self, name, directory, molfile ='', attemptNumber = 0, multiplicity = -1):
        self.name = name
        
        self.directory = directory
        
        self.molfile = molfile
        
        self.attemptNumber = attemptNumber
        
        self.multiplicity = multiplicity
        
        self.keywords = ''
         
class MOPACKEYWORDS:
    BOTH = 'Both'
    BOTTOM = 'Bottom'
    TOP = 'Top'
    POLAR = 'Polar'
    
class MOPACPM3InputWriter(QMInputWriter):
    
    scriptAttempts = 5
    maxAttemptNumber = 2* scriptAttempts
    
    def fillMultiplicityKeywords(self):
        self.multiplicityKeywords[1] = ''
        self.multiplicityKeywords[2] = 'uhf doublet'
        self.multiplicityKeywords[3] = 'uhf triplet'
        self.multiplicityKeywords[4] = 'uhf quartet'
        self.multiplicityKeywords[5] = 'uhf quintet'
        self.multiplicityKeywords[6] = 'uhf sextet'
        self.multiplicityKeywords[7] = 'uhf septet'
        self.multiplicityKeywords[8] = 'uhf octet'
        self.multiplicityKeywords[9] = 'uhf nonet'
        
    def __init__(self, name, directory, p_molfile, attemptNumber, multiplicity):
        QMInputWriter.__init__(self, name = name, directory = directory, molfile = p_molfile, attemptNumber = attemptNumber, multiplicity = multiplicity)
        self.multiplicityKeywords = {}
        self.fillMultiplicityKeywords()
        self.keywords = {}
        
        self.inputExtension = '.mop'
   
    def getMopacRadicalString(self, multipl):
        '''
        
        returns the extra Mopac keywords to use for radical species, given the spin multiplicity 
        (radical number + 1)
     
        '''
        try:
            keyword = self.multiplicityKeywords.get(multipl)
        except Exception as e:
            logging.critical('Invalid multiplicity encountered: "+multipl'+str(multipl))
        return keyword  
    
    
     
    def createKeywords(self):
        radicalString = self.getMopacRadicalString(self.multiplicity)
        try:
            if self.attemptNumber % MOPACPM3InputWriter.scriptAttempts == 1:
                inpKeyStrBoth="pm3 "+radicalString
                inpKeyStrTop=" precise nosym"
                '''
                 * 7/10/09: based on a quick review of recent results, 
                 * keyword combo #1 rarely works, and when it did (CJAINEUZFLXGFA-UHFFFAOYAUmult3 
                 * (InChI=1/C8H16O5Si/c1-4-11-14(9,12-5-2)13-8-6-10-7(8)3/h7-8H,3-6H2,1-2H3/mult3)), 
                 * the grad. norm on the force step was about 1.7 (too large) 
                 * I manually removed this result and re-ran...the entropy was increased 
                 * by nearly 20 cal/mol-K...perhaps we should add a check for the "WARNING" 
                 * that MOPAC prints out when the gradient is high 7/22/09: for the case of 
                 * FUGDBSHZYPTWLG-UHFFFAOYADmult3 (InChI=1/C5H8/c1-4-3-5(4)2/h4-5H,1-3H2/mult3), 
                 * adding nosym seemed to resolve 1. large discrepancies from Gaussian and 2. 
                 * negative frequencies in mass-weighted coordinates and possibly related issue 
                 * in discrepancies between regular and mass-weighted coordinate frequencies
                 '''
                inpKeyStrBottom="oldgeo thermo nosym precise "
            elif self.attemptNumber%self.scriptAttempts == 2:
                '''
                 * #7/9/09: used for VCSJVABXVCFDRA-UHFFFAOYAI 
                 * (InChI=1/C8H19O5Si/c1-5-10-8(4)13-14(9,11-6-2)12-7-3/h8H,5-7H2,1-4H3) 
                 * all existing Gaussian keywords also failed the Gaussian result 
                 * was also rectified, but the resulting molecule was over 70 kcal/mol less stable,
                 *  probably due to a large amount of spin contamination (~1.75 in fixed Gaussian 
                 *  result vs. 0.754 for MOPAC)
                '''
                inpKeyStrBoth="pm3 "+radicalString
                inpKeyStrTop=" precise nosym gnorm=0.0 nonr"
                inpKeyStrBottom="oldgeo thermo nosym precise "
            elif self.attemptNumber%self.scriptAttempts == 3:
                '''
                 * #7/8/09: used for ADMPQLGIEMRGAT-UHFFFAOYAUmult3 
                 * (InChI=1/C6H14O5Si/c1-4-9-12(8,10-5-2)11-6(3)7/h6-7H,3-5H2,1-2H3/mult3) 
                 * all existing Gaussian keywords also failed however, the Gaussian result
                 *  was also rectified, and the resulting conformation was about 1.0 kcal/mol 
                 *  more stable than the one resulting from this, so fixed Gaussian result was 
                 *  manually copied to QMFiles folder
                '''
                inpKeyStrBoth="pm3 "+radicalString
                inpKeyStrTop=" precise nosym gnorm=0.0"
                '''
                 * #precise appeared to be necessary for the problematic case (to avoid negative frequencies)
                '''
                inpKeyStrBottom="oldgeo thermo nosym precise " 
            elif self.attemptNumber%self.scriptAttempts == 4:
                '''
                 * #7/8/09: used for GYFVJYRUZAKGFA-UHFFFAOYALmult3 
                 * (InChI=1/C6H14O6Si/c1-3-10-13(8,11-4-2)12-6-5-9-7/h6-7H,3-5H2,1-2H3/mult3) 
                 * case (negative frequency issues in MOPAC) (also, none of the existing Gaussian
                 *  combinations worked with it) note that the Gaussian result appears to be a 
                 *  different conformation as it is about 0.85 kcal/mol more stable, so the Gaussian 
                 *  result was manually copied to QMFiles directory note that the MOPAC output included 
                 *  a very low frequency (4-5 cm^-1)
                '''
                inpKeyStrBoth="pm3 "+radicalString
                inpKeyStrTop=" precise nosym gnorm=0.0 bfgs"
                '''
                 * precise appeared to be necessary for the problematic case (to avoid negative frequencies)
                '''
                inpKeyStrBottom="oldgeo thermo nosym precise " 
            elif self.attemptNumber%self.scriptAttempts == 0:
                '''
                 * used for troublesome HGRZRPHFLAXXBT-UHFFFAOYAVmult3 
                 * (InChI=1/C3H2O4/c4-2(5)1-3(6)7/h1H2/mult3) case 
                 * (negative frequency and large gradient issues)                
                '''
                inpKeyStrBoth="pm3 "+radicalString
                inpKeyStrTop=" precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000"
                inpKeyStrBottom="oldgeo thermo nosym precise "
        except Exception as e:
            logging.error('Error in writing inputkeywords.txt \n'+str(e))
            
            
        if qmtp.QMTP.usePolar:
            if self.multiplicity == 1:
               polarString = "\n" + "\n" + "\n"+ "oldgeo polar nosym precise " + inpKeyStrBoth
            else:
               polarString = "\n" + "\n" + "\n"+ "oldgeo static nosym precise " + inpKeyStrBoth
            self.keywords[MOPACKEYWORDS.POLAR] = polarString
        else:
            self.keywords[MOPACKEYWORDS.POLAR] = ''
            
        self.keywords[MOPACKEYWORDS.BOTH] = inpKeyStrBoth
        self.keywords[MOPACKEYWORDS.TOP] = inpKeyStrTop
        self.keywords[MOPACKEYWORDS.BOTTOM] = inpKeyStrBottom

        
        return self.keywords
    
    
    def write(self):
        self.createKeywords()

        return self.createInputFile()
    
    def createInputFile(self):
        '''
        TODO instead calling the external executable, import obabel directly!!!
        '''
        inpKeyStrTopCombined = self.keywords[MOPACKEYWORDS.BOTH] + self.keywords[MOPACKEYWORDS.TOP]
        if self.attemptNumber <= self.scriptAttempts: #use UFF-refined coordinates
                command = ["babel", "-imol", self.molfile.path, "-xk", inpKeyStrTopCombined,"--title", self.molfile.molecule.toAugmentedInChI(),"-omop", os.path.join(self.directory,self.name + self.inputExtension) ]
        else:
                command = [ "babel", "-imol", self.molfile.crudepath, "-xk",  inpKeyStrTopCombined,"--title", self.molfile.molecule.toAugmentedInChI(),"-omop",os.path.join(self.directory,self.name + self.inputExtension) ]        
        try:   
            process = Popen(command)
            process.communicate()#necessary to wait for completion of process!
        except Exception as e:
            err = 'Error in running OpenBabel MOL to MOP process \n' + str(e)
            logging.error(err)   

        with open(os.path.join(self.directory,self.name + self.inputExtension), 'a') as mop:#append 'a' instead of overwrite 'w'
            mop.write('\n'+self.keywords[MOPACKEYWORDS.BOTTOM]+self.keywords[MOPACKEYWORDS.BOTH]+self.keywords[MOPACKEYWORDS.POLAR])
        
        return self.name + '.mop'
