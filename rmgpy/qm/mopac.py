import os
from subprocess import Popen, PIPE

import openbabel

class Mopac:
    
    directory = 'QMfiles'
    inputFileExtension = '.mop'
    outputFileExtension = '.out'
    executablePath = os.path.join(os.getenv('MOPAC_DIR', default="/opt/mopac") , 'MOPAC2009.exe')
    
    scriptAttempts = 5
    maxAttempts = 2 * scriptAttempts
    
    usePolar = False#use polar keyword in MOPAC
    
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
    
    failureKeys = ['IMAGINARY FREQUENCIES', 'EXCESS NUMBER OF OPTIMIZATION CYCLES', 'NOT ENOUGH TIME FOR ANOTHER CYCLE']
    
    def writeInputFile(self, top_keys, bottom_keys, polar_keys, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        inputFilePath = os.path.join(self.directory , self.geometry.uniqueID + self.inputFileExtension)
        
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "mop")
        mol = openbabel.OBMol()
    
        if attempt <= self.scriptAttempts: #use UFF-refined coordinates
            obConversion.ReadFile(mol, self.geometry.getRefinedMolFilePath() )
        else:
            obConversion.ReadFile(mol, self.geometry.getCrudeMolFilePath() )
    
        mol.SetTitle(self.geometry.uniqueID) 
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    
        input_string = obConversion.WriteString(mol)
        
        with open(inputFilePath, 'w') as mopacFile:
            mopacFile.write(top_keys)
            mopacFile.write(input_string)
            mopacFile.write('\n')
            mopacFile.write(bottom_keys)
            if self.usePolar:
                mopacFile.write('\n\n\n')
                mopacFile.write(polar_keys)
        
        return self.geometry.uniqueID + self.inputFileExtension
       
    def run(self, inputFileName):
        # submits the input file to mopac
        command = os.path.join(self.directory, self.geometry.uniqueID + self.inputFileExtension)
        process = Popen([self.executablePath, command])
        process.communicate()# necessary to wait for executable termination!
    
        return self.checkNoFailure()
        
    def checkNoFailure(self):
        """
        checks whether the output file contains any of the 
        failure keywords
        """
        file = os.path.join(self.directory,self.geometry.uniqueID+self.outputFileExtension)
        with open(file) as qmfile:    
            for each_line in qmfile:
                each_line = each_line.rstrip().strip()
                for element in self.failureKeys:#search for failure keywords
                    if element in each_line:
                        logging.error("MOPAC output file contains the following error %s")%element
                        return False
                    
        return True
        
    def read():
        # reads the output file
        import ipdb; ipdb.set_trace()
        pass
    
    def calculate():
        # calculators for the parsing
        import ipdb; ipdb.set_trace()
        pass
        
    def parse():
        # parses the output file to generate the TPs
        import ipdb; ipdb.set_trace()
        pass