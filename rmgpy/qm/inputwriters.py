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

        
    def __init__(self, molfile, attemptNumber, multiplicity):
        QMInputWriter.__init__(self, molfile, attemptNumber, multiplicity)


    def write(self):
        return self.createInputFile()
    
    def createInputFile(self):
        # moved
        
class G03PM3KEYWORDS:
    INPUT = 'Input'
    
    
class GaussianPM3InputWriter(QMInputWriter):
    
    """static fields"""
    scriptAttempts = 18
    
    maxAttemptNumber = 2* scriptAttempts
    
    keywords = {}
    
    inputExtension = '.gjf'
    
    keywordsTop = {}#keywords that will be added to the qm input file based on the attempt number
    keywordsTop[1] = "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3)"
    keywordsTop[2] = "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)"
    keywordsTop[3] = "# pm3 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm" 
    keywordsTop[4] = "# pm3 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm"
    keywordsTop[5] = "# pm3 opt=(verytight,gdiis,small) freq IOP(2/16=3)"
    keywordsTop[6] = "# pm3 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)"
    keywordsTop[7] = "# pm3 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)"
    keywordsTop[8] = "# pm3 opt=tight freq IOP(2/16=3)"
    keywordsTop[9] = "# pm3 opt=tight freq=numerical IOP(2/16=3)"
    keywordsTop[10] = "# pm3 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)"
    keywordsTop[11] = "# pm3 opt freq IOP(2/16=3)"
    keywordsTop[12] = "# pm3 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)"
    keywordsTop[13] = "# pm3 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm"
    keywordsTop[14] = "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm"
    keywordsTop[15] = "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm"
    keywordsTop[16] = "# pm3 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm"
    keywordsTop[17] = "# pm3 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm"
    keywordsTop[18] = "# pm3 opt=(calcall,small,maxcyc=100) IOP(2/16=3)"
    
    def __init__(self, molfile, attemptNumber, multiplicity):
        QMInputWriter.__init__(self, molfile, attemptNumber, multiplicity)

        
    def createInputFile(self):
        
        keywords = "\n".join((
            "%chk={0}/RMGrunCHKfile.chk".format(self.molfile.directory),
            "%mem=6MW",
            "%nproc=1",
            self.keywordsTop[self.attemptNumber] + ' polar' if qmtp.QMTP.usePolar else ''
            ))
        inputFilePath = os.path.join(self.molfile.directory,self.molfile.name + self.inputExtension)
        
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "gjf")
        mol = openbabel.OBMol()
        
        if self.attemptNumber <= GaussianPM3InputWriter.scriptAttempts: #use UFF-refined coordinates
            obConversion.ReadFile(mol, self.molfile.path)
        else:
            obConversion.ReadFile(mol, self.molfile.crudepath)
    
        mol.SetTitle(self.molfile.molecule.toAugmentedInChI())
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
        
        input_string = obConversion.WriteString(mol)
        
        with open(inputFilePath, 'w') as gjfFile:
            gjfFile.write(keywords)
            gjfFile.write(input_string)

        return self.molfile.name + self.inputExtension
    
    def write(self):
        inputFile = self.createInputFile()
        return inputFile
    
