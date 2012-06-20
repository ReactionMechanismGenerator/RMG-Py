"""
Module that writes input files for the various QM packages.

Input files contain information on the 3D coordinates of the molecule, but also 
command strings to specify the level of theory, and other flags.


"""

import os

import openbabel

import qmtp

class QMInputWriter:
    """
    Supertype for all input writers for quantum chemistry methods
    """
    
    "the number of keyword permutations available. update as additional options are added"
    scriptAttempts = 0
    
    #we will try a second time with crude coordinates if the UFF refined coordinates do not work
    "Maximum number of attempts"
    maxAttemptNumber = 0

    def __init__(self, molfile='', attemptNumber=0, multiplicity=-1):
        
        self.molfile = molfile
        self.attemptNumber = attemptNumber
        self.multiplicity = multiplicity

class MOPACPM3InputWriter(QMInputWriter):
    
    scriptAttempts = 5
    maxAttemptNumber = 2 * scriptAttempts
    
    inputExtension = '.mop'
    
    "Keywords for the multiplicity"
    multiplicityKeywords = {}
    multiplicityKeywords[1] = ''
    multiplicityKeywords[2] = 'uhf doublet'
    multiplicityKeywords[3] = 'uhf triplet'
    multiplicityKeywords[4] = 'uhf quartet'
    multiplicityKeywords[5] = 'uhf quintet'
    multiplicityKeywords[6] = 'uhf sextet'
    multiplicityKeywords[7] = 'uhf septet'
    multiplicityKeywords[8] = 'uhf octet'
    multiplicityKeywords[9] = 'uhf nonet'
    
    "Keywords that will be added at the top of the qm input file"
    keywordsTop = {}
    keywordsTop[1] = "precise nosym"
    keywordsTop[2] = "precise nosym gnorm=0.0 nonr"
    keywordsTop[3] = "precise nosym gnorm=0.0"
    keywordsTop[4] = "precise nosym gnorm=0.0 bfgs"
    keywordsTop[5] = "precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000"
    
    "Keywords that will be added at the bottom of the qm input file"
    keywordsBottom = {}
    keywordsBottom[1] = "oldgeo thermo nosym precise "
    keywordsBottom[2] = "oldgeo thermo nosym precise "
    keywordsBottom[3] = "oldgeo thermo nosym precise "
    keywordsBottom[4] = "oldgeo thermo nosym precise "
    keywordsBottom[5] = "oldgeo thermo nosym precise "
        
    def __init__(self, molfile, attemptNumber, multiplicity):
        QMInputWriter.__init__(self, molfile, attemptNumber, multiplicity)


    def write(self):
        return self.createInputFile()
    
    def createInputFile(self):
        multiplicity_keywords = self.multiplicityKeywords[self.multiplicity]
        top_keywords = "pm3 {0} {1}".format(
                multiplicity_keywords,
                self.keywordsTop[self.attemptNumber],
                )
        bottom_keywords = "{0} pm3 {1}".format(
                self.keywordsBottom[self.attemptNumber],
                multiplicity_keywords,
                )
        polar_keywords = "oldgeo {0} nosym precise pm3 {1}".format(
                'polar' if self.multiplicity == 1 else 'static',
                multiplicity_keywords,
                )
        
        inputFilePath = os.path.join(self.molfile.directory, self.molfile.name + self.inputExtension)
        
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "mop")
        mol = openbabel.OBMol()
        
        if self.attemptNumber <= self.scriptAttempts: #use UFF-refined coordinates
            obConversion.ReadFile(mol, self.molfile.path)
        else:
            obConversion.ReadFile(mol, self.molfile.crudepath)
        
        mol.SetTitle(self.molfile.molecule.toAugmentedInChI()) 
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)

        input_string = obConversion.WriteString(mol)
        
        with open(inputFilePath, 'w') as mopacFile:
            mopacFile.write(top_keywords)
            mopacFile.write(input_string)
            mopacFile.write('\n')
            mopacFile.write(bottom_keywords)
            if qmtp.QMTP.usePolar:
                mopacFile.write('\n\n\n')
                mopacFile.write(polar_keywords)
        
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
