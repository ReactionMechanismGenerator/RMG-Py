class Mopac:
    
    inputFileExtension = '.mop'
    outputFileExtension = '.out'
    executablePath = os.path.join(os.getenv('MOPAC_DIR', default="/opt/mopac") , 'MOPAC2009.exe')
    
    script_attempts = 5
    max_attempts = 2 * script_attempts
    
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
    
    # **may need to have some check to verify the mopac directory    
    
    def __init__(self):
        self.failureKeys = ['IMAGINARY FREQUENCIES', 'EXCESS NUMBER OF OPTIMIZATION CYCLES', 'NOT ENOUGH TIME FOR ANOTHER CYCLE']
        
    def writeInputFile(self, top_keywords, bottom_keywords, polar_keywords, geometry, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        multiplicity_keyword = this.multiplicityKeywords[geometry.multiplicity]
    
        inputFilePath = os.path.join( where_stuff_lives , geometry.uniqueID + self.inputExtension)
    
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "mop")
        mol = openbabel.OBMol()
    
        if attempt <= self.scriptAttempts: #use UFF-refined coordinates
            obConversion.ReadFile(mol, geometry.getRefinedMolFilePath() )
        else:
            obConversion.ReadFile(mol, geometry.getCrudeMolFilePath() )
    
        mol.SetTitle( geometry.uniqueID ) 
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
       
    def run(self):
        # submits the input file to mopac
        process = Popen([self.executablePath, self.command])
        process.communicate()# necessary to wait for executable termination!
    
        return self.check()
        
    def checkNoFailure(self):
        """
        checks whether the output file contains any of the 
        failure keywords
        """
        file = os.path.join(self.molfile.directory,self.molfile.name+self.outputExtension)
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
    
    def calculate():
        # calculators for the parsing
        
    def parse():
        # parses the output file to generate the TPs