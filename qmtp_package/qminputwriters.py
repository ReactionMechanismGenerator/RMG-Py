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
import openbabel


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
        '''
        Based on the attempt number keywords will be added to an QM input file.
        '''
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
        obConversion.WriteFile(mol, os.path.join(self.directory,self.name + self.inputExtension))

        #pre-pend keywords:
        with open(os.path.join(self.directory,self.name + self.inputExtension), 'r+') as mop:
            old = mop.read() # read everything in the file
            mop.seek(0) # rewind
            mop.write(inpKeyStrTopCombined + old) # write the new line before

        with open(os.path.join(self.directory,self.name + self.inputExtension), 'a') as mop:#append 'a' instead of overwrite 'w'
            mop.write('\n'+self.keywords[MOPACKEYWORDS.BOTTOM]+self.keywords[MOPACKEYWORDS.BOTH]+self.keywords[MOPACKEYWORDS.POLAR])
        
        return self.name + '.mop'
