"""
Module that writes input files for the various QM packages.

Input files contain information on the 3D coordinates of the molecule, but also 
command strings to specify the level of theory, and other flags.


"""

import qmtp
import platform
from subprocess import Popen
from symmetry import *
import logging
import os
import math
import openbabel

class QMInputWriter:
    """
     Supertype for all input writers for quantum chemistry methods
    """
    
    #the number of keyword permutations available update as additional options are added
    scriptAttempts = 0
    
    #we will try a second time with crude coordinates if the UFF refined coordinates do not work
    maxAttemptNumber = 0

        
    def __init__(self, molfile ='', attemptNumber = 0, multiplicity = -1):
        
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
        
    def __init__(self, p_molfile, attemptNumber, multiplicity):
        QMInputWriter.__init__(self, p_molfile, attemptNumber, multiplicity)
        self.multiplicityKeywords = {}
        self.fillMultiplicityKeywords()
        self.keywords = {}
        
        self.keywordsTop = {}#keywords that will be added at the top of the qm input file
        self.keywordsTop[1] = " precise nosym"
        self.keywordsTop[2] = " precise nosym gnorm=0.0 nonr"
        self.keywordsTop[3] = " precise nosym gnorm=0.0"
        self.keywordsTop[4] = " precise nosym gnorm=0.0 bfgs"
        self.keywordsTop[5] = " precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000"
        
        self.keywordsBottom = {}#keywords that will be added at the bottom of the qm input file
        self.keywordsBottom[1] = "oldgeo thermo nosym precise "
        self.keywordsBottom[2] = "oldgeo thermo nosym precise "
        self.keywordsBottom[3] = "oldgeo thermo nosym precise "
        self.keywordsBottom[4] = "oldgeo thermo nosym precise "
        self.keywordsBottom[5] = "oldgeo thermo nosym precise "
        
        self.inputExtension = '.mop'

    def createKeywords(self):
        """
        Based on the attempt number keywords will be added to an QM input file.
        """
        inpKeyStrBoth = "pm3 "+self.multiplicityKeywords[self.multiplicity]
        inpKeyStrTop = self.keywordsTop[self.attemptNumber]
        inpKeyStrBottom = self.keywordsBottom[self.attemptNumber]
           
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

        inpKeyStrTopCombined = self.keywords[MOPACKEYWORDS.BOTH] + self.keywords[MOPACKEYWORDS.TOP]
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "mop")
        mol = openbabel.OBMol()
        
        if self.attemptNumber <= self.scriptAttempts: #use UFF-refined coordinates
            obConversion.ReadFile(mol, self.molfile.path)
        else:
            obConversion.ReadFile(mol, self.molfile.crudepath)
    
        mol.SetTitle(self.molfile.molecule.toAugmentedInChI()) 
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
        'TODO still dont know how to write keywords, therefore they will be prepended afterwards...'
        obConversion.WriteFile(mol, os.path.join(self.molfile.directory,self.molfile.name + self.inputExtension))

        #pre-pend keywords:
        with open(os.path.join(self.molfile.directory,self.molfile.name + self.inputExtension), 'r+') as mop:
            old = mop.read() # read everything in the file
            mop.seek(0) # rewind
            mop.write(inpKeyStrTopCombined + old) # write the new line before

        with open(os.path.join(self.molfile.directory,self.molfile.name + self.inputExtension), 'a') as mop:#append 'a' instead of overwrite 'w'
            mop.write('\n'+self.keywords[MOPACKEYWORDS.BOTTOM]+self.keywords[MOPACKEYWORDS.BOTH]+self.keywords[MOPACKEYWORDS.POLAR])
        
        return self.molfile.name + '.mop'

class G03PM3KEYWORDS:
    INPUT = 'Input'
    
    
class GaussianPM3InputWriter(QMInputWriter):
    
    """static fields"""
    scriptAttempts = 18
    
    maxAttemptNumber = 2* scriptAttempts
    
    def __init__(self, p_molfile, attemptNumber, multiplicity):
        QMInputWriter.__init__(self, p_molfile, attemptNumber, multiplicity)
        
        self.keywords = {}
        
        self.inputExtension = '.gjf'
        
        self.keywordsTop = {}#keywords that will be added to the qm input file based on the attempt number
        self.keywordsTop[1] = "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3)"
        self.keywordsTop[2] = "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)"
        self.keywordsTop[3] = "# pm3 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm" 
        self.keywordsTop[4] = "# pm3 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm"
        self.keywordsTop[5] = "# pm3 opt=(verytight,gdiis,small) freq IOP(2/16=3)"
        self.keywordsTop[6] = "# pm3 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)"
        self.keywordsTop[7] = "# pm3 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)"
        self.keywordsTop[8] = "# pm3 opt=tight freq IOP(2/16=3)"
        self.keywordsTop[9] = "# pm3 opt=tight freq=numerical IOP(2/16=3)"
        self.keywordsTop[10] = "# pm3 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)"
        self.keywordsTop[11] = "# pm3 opt freq IOP(2/16=3)"
        self.keywordsTop[12] = "# pm3 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)"
        self.keywordsTop[13] = "# pm3 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm"
        self.keywordsTop[14] = "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm"
        self.keywordsTop[15] = "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm"
        self.keywordsTop[16] = "# pm3 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm"
        self.keywordsTop[17] = "# pm3 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm"
        self.keywordsTop[18] = "# pm3 opt=(calcall,small,maxcyc=100) IOP(2/16=3)"
        
    def createInputFile(self):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "gjf")
        mol = openbabel.OBMol()
        
        if self.attemptNumber <= GaussianPM3InputWriter.scriptAttempts: #use UFF-refined coordinates
            obConversion.ReadFile(mol, self.molfile.path)
        else:
            obConversion.ReadFile(mol, self.molfile.crudepath)
    
        mol.SetTitle(self.molfile.molecule.toAugmentedInChI()) 
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
        'TODO still dont know how to write keywords, therefore they will be prepended afterwards...'
        obConversion.WriteFile(mol, os.path.join(self.molfile.directory,self.molfile.name + self.inputExtension))

        #pre-pend keywords:
        with open(os.path.join(self.molfile.directory,self.molfile.name + self.inputExtension), 'r+') as gjf:
            old = gjf.read() # read everything in the file
            gjf.seek(0) # rewind
            gjf.write(self.keywords[G03PM3KEYWORDS.INPUT] + old) # write the new line before
 
        return self.molfile.name+".gjf"
    
    def write(self):
        self.createKeywords()
        inputFile = self.createInputFile()
        return inputFile
    
    def createKeywords(self):
        inpKeyStr ="%chk=" + self.molfile.directory + "/RMGrunCHKfile.chk\n"
        inpKeyStr =inpKeyStr + "%mem=6MW\n"
        inpKeyStr =inpKeyStr + "%nproc=1\n"
        inpKeyStr =inpKeyStr + self.keywordsTop[self.attemptNumber]
        
        if qmtp.QMTP.usePolar:
            inpKeyStr = inpKeyStr+" polar"

        self.keywords[G03PM3KEYWORDS.INPUT] = inpKeyStr
        
        return self.keywords
