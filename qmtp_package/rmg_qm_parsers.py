
'''
Module that collects all parsers that read in a output file of a successfully finished 
QM run (G03, MOPAC, MM4, etc...) and parses this to a workable object using an external parser library.

Next this data is converted in to useful thermodynamic quantities such as H, S, Cp by calling a calculator module.

 
'''

from rmgpy.molecule import Molecule
from rmgpy.thermo import ThermoData
import calculator as calc
import logging
import re
from subprocess import Popen
import platform
import os
from qmdata import CCLibData
from cclib.parser import ccopen, MM4, Mopac

class QMParser:
    '''
    
    All sub-classes of QMParser have the following attributes:
        *inputFileExtension
        *executable #the keywords denoting the executable
        *command ##the command line command
        *scriptFile##parsing script file
        *parsingTool
        *qmdata
    '''
    
    def __init__(self, name, directory, molecule, qmtp, environ = os.environ.get("RMG_workingDirectory")):
        self.name = name
        self.directory = directory
        self.molecule = molecule
        
        self.environ = environ#working directory
        
        self.qmtp = qmtp
        
        self.scripts = "scripts/"##directory of parsing scripts

        
    def read(self):
            self.qmdata = self.parsingTool.parse(self.molecule)    
            return self.qmdata

    def parse(self) :
        self.read()
        
        calculator = calc.TDPropertiesCalculator(self.name, self.directory, self.qmdata, environ = self.environ)
        
        return calculator.calculate()

class CCLibParser:
    def __init__(self, path, qmtp, qmdata = ''):
        self.path = path #path to QM output file being parsed
        self.qmdata = qmdata
        self.qmtp = qmtp
        
    def parse(self, molecule):
        try:
            #parse the Mopac file using cclib
            '''
            
            energy = 0#PM3 energy (Hf298) in Hartree (***note: in the case of MOPAC, the MOPAC file will contain in units of kcal/mol, but modified ccLib will return in Hartree)
            
             * calculate ground state degeneracy from the number of radicals 
             * this should give the same result as spin multiplicity in Gaussian 
             * input file (and output file), but we do not explicitly check this 
             * (we could use "mult" which cclib reads in if we wanted to do so) also,
             *  note that this is not always correct, as there can apparently be additional 
             *  spatial degeneracy for non-symmetric linear molecules like OH radical 
             *  (cf. http:#cccbdb.nist.gov/thermo.asp)            
            '''
            if self.qmtp.qmprogram == 'mopac':
                myfile=ccopen(self.path)
            elif self.qmtp.qmprogram == 'mm4':
                myfile=MM4(self.path)
            elif self.qmtp.qmprogram == 'gaussian03':
                myfile=Mopac(self.path)
                
            myfile.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
            
            cclib_data=myfile.parse()
            
            radicalNumber = sum([i.radicalElectrons for i in molecule.atoms])
            
            self.qmdata = CCLibData(cclib_data, radicalNumber + 1)
            
        except Exception as e:
            logging.error('Error in reading/parsing ccLib Python process \n')
            logging.error(str(e))

        return self.qmdata
        
class MOPACPM3Parser(QMParser):
    def __init__(self, name, directory, molecule, qmtp, environ = os.environ.get("RMG_workingDirectory")):
        QMParser.__init__(self, name, directory, molecule, qmtp, environ)
        
        self.inputFileExtension = ".out"
        
        path = os.path.join(self.directory,self.name+self.inputFileExtension)
        
        self.parsingTool = CCLibParser(path, qmtp)
            