"""
Module that collects all parsers that read in a output file of a successfully finished 
QM run (G03, MOPAC, MM4, etc...) and parses this to a workable object using an external parser library.

Next this data is converted in to useful thermodynamic quantities such as H, S, Cp by calling a calculator module.

"""

import os

import cclib.parser

import logging
from calculator import TDPropertiesCalculator
from qmdata import CCLibData

class QMParser:
    """
    All sub-classes of QMParser have the following attributes:
        *inputFileExtension
        *executable #the keywords denoting the executable
        *command ##the command line command
        *scriptFile##parsing script file
        *parsingTool
        *qmdata
    """
    
    def __init__(self, molfile, qmtp):
        self.molfile = molfile
        
        self.qmtp = qmtp

    def read(self):
        self.qmdata = self.parsingTool.parse(self.molfile.molecule)
        return self.qmdata

    def parse(self):
        self.read()
        
        calculator = TDPropertiesCalculator(self.molfile, self.qmdata)
        
        return calculator.calculate()

class CCLibParser:
    def __init__(self, path, qmtp, qmdata = ''):
        self.path = path #path to QM output file being parsed
        self.qmdata = qmdata
        self.qmtp = qmtp
        
    def parse(self, molecule):
        try:
            #parse the Mopac file using cclib
            """
            energy = 0#PM3 energy (Hf298) in Hartree (***note: in the case of MOPAC, the MOPAC file will contain in units of kcal/mol, but modified ccLib will return in Hartree)
            
             * calculate ground state degeneracy from the number of radicals 
             * this should give the same result as spin multiplicity in Gaussian 
             * input file (and output file), but we do not explicitly check this 
             * (we could use "mult" which cclib reads in if we wanted to do so) also,
             *  note that this is not always correct, as there can apparently be additional 
             *  spatial degeneracy for non-symmetric linear molecules like OH radical 
             *  (cf. http:#cccbdb.nist.gov/thermo.asp)            
            """
            if self.qmtp.qmprogram == 'gaussian03':
                myfile=cclib.parser.ccopen(self.path)
            elif self.qmtp.qmprogram == 'mm4':
                myfile=cclib.parser.MM4(self.path)
            elif self.qmtp.qmprogram == 'mopac':
                myfile=cclib.parser.Mopac(self.path)
                
            myfile.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
            
            cclib_data=myfile.parse()
            
            radicalNumber = sum([i.radicalElectrons for i in molecule.atoms])
            
            self.qmdata = CCLibData(cclib_data, radicalNumber + 1)
            
        except Exception as e:
            logging.error('Error in reading/parsing ccLib Python process \n')
            logging.error(str(e))
            raise

        return self.qmdata
        
class MOPACPM3Parser(QMParser):
    def __init__(self, molfile, qmtp):
        QMParser.__init__(self, molfile, qmtp)
        
        self.inputFileExtension = ".out"
        
        path = os.path.join(self.molfile.directory,self.molfile.name+self.inputFileExtension)
        
        self.parsingTool = CCLibParser(path, qmtp)

