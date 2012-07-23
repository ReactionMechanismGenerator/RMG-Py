import os

import openbabel
import cclib.parser
import logging
from subprocess import Popen, PIPE

from qmdata import CCLibData
import symmetry

class Mopac:
    
    directory = 'QMfiles'
    inputFileExtension = '.mop'
    outputFileExtension = '.out'
    executablePath = os.path.join(os.getenv('MOPAC_DIR', default="/opt/mopac") , 'MOPAC2009.exe')
    
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
    
    import ipdb; ipdb.set_trace()
    scriptAttempts = len(keywordsTop)
    maxAttempts = 2 * scriptAttempts
    
    failureKeys = ['IMAGINARY FREQUENCIES', 'EXCESS NUMBER OF OPTIMIZATION CYCLES', 'NOT ENOUGH TIME FOR ANOTHER CYCLE']
    
    def writeInputFile(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        from mopacpm3mol import MopacPM3
        
        inputFilePath = os.path.join(self.directory , self.geometry.uniqueID + self.inputFileExtension)
        
        method = MopacPM3(self)
        top_keys, bottom_keys, polar_keys = method.inputFileKeys(attempt)
        
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
        
    def parse(self):
        # parses the output file to generate the TPs
        path = os.path.join(self.directory,self.geometry.uniqueID+'.out')
        try:
            myfile = cclib.parser.Mopac(path)
            myfile.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
            
            cclibData = myfile.parse()
            radicalNumber = sum([i.radicalElectrons for i in self.molecule.atoms])
            
            qmData = CCLibData(cclibData, radicalNumber+1)
                
        except Exception as e:
            logging.error('Error in reading/parsing ccLib Python process.')
            logging.error(str(e))
            raise
        
        return qmData
        
class SymmetryJob:
    """
    Determine the point group using the SYMMETRY program 
    (http://www.cobalt.chem.ucalgary.ca/ps/symmetry/).


    Required input is a line with number of atoms followed by lines for each atom 
    including:
    1) atom number
    2) x,y,z coordinates

    finalTol determines how loose the point group criteria are;
    values are comparable to those specified in the GaussView point group interface

    """
    directory = 'QMfiles'
    
    RMG_path = os.environ.get("RMGpy")
    if RMG_path is None:
        RMG_path = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..'))
        logging.info("Setting RMG_path to {0}".format(RMG_path))
    executable_path = os.path.join(RMG_path, 'bin', 'symmetry')
    if not os.path.exists(executable_path):
        raise Exception("Symmetry program not found at {0}.".format(executable_path))

    'Argumenst that will be passed as an argument for the consecutive attempts'
    argumentsList = [
         ['-final', '0.02'],
         ['-final', '0.1'],
         ['-primary', '0.2', '-final' ,'0.1'],
         ['-final', '0.0'],
        ]

    inputFileExtension = '.symm'

    def __init__(self, molfile, iqmdata):
        self.molfile = molfile
        
        'the command line command'
        self.command = []

        "IQMData is the object that holds information from a previous QM Job on 3D coords, molecule etc..."
        self.qmdata = iqmdata

        self.inputFile = self.molfile + self.inputFileExtension

        self.attemptNumber = 1
        self.pointGroupFound = False

    def check(self, output):
        """
        Check the `output` string and extract the resulting point group, which is returned.
        """
        for line in output.split('\n'):
            if line.startswith("It seems to be the "): # "It seems to be the [x] point group" indicates point group.
                result = line.split(" ")[5]
                break
        else:
            logging.exception("Couldn't find point group from symmetry output:\n{0}".format(output))
            raise RuntimeError("Couldn't find point group in symmetry output.")

        logging.info("Point group: "+ result)
        return result;

    def run(self):
        """
        Run the command and wait for it to finish.
        """
        pp = Popen(self.command, stdout=PIPE, stderr=PIPE)
        stdout, stderr = pp.communicate()
        return self.check(stdout)    

    def writeInputFile(self):
        """
        Write the input file for the SYMMETRY program.
        """
        geom = str(self.qmdata.numberOfAtoms) + "\n"
        for i in range(self.qmdata.numberOfAtoms):
            geom = geom + " ".join((str(self.qmdata.atomicNumbers[i]),
                                    str(self.qmdata.atomCoords[i][0]),
                                    str(self.qmdata.atomCoords[i][1]),
                                    str(self.qmdata.atomCoords[i][2])
                                   )) + "\n"
        with open(os.path.join(self.directory, self.inputFile), 'w') as input_file:
            input_file.write(geom)
        input_file.close()
        logging.info("Symmetry input file written to %s"%os.path.join(self.directory, self.inputFile))
        return input_file

    def calculate(self):
        """
        Do the entire point group calculation.

        This writes the input file, then tries several times to run 'symmetry'
        with different parameters, until a point group is found and returned.
        """
        self.writeInputFile();

        result = "";

        #continue trying to generate symmetry group until too many no. of attempts or until a point group is found: 
        for attempt, arguments in enumerate(self.argumentsList):
            """
            TODO only *nix case works!
            """
            self.command = [self.executable_path]
            self.command.extend(arguments)
            self.command.append(os.path.join(self.directory, self.inputFile))

            #call the program and read the result
            result = self.run();

            #check for a recognized point group
            symmPGC = symmetry.PointGroupDictionary()

            if symmPGC.contains(result):
                self.pointGroupFound = True;
                return symmetry.PointGroup(result)
            else:
                logging.info("Attempt number {0} did not identify a recognized point group ({1}).".format(attempt,result))
        logging.critical("Final attempt did not identify a recognized point group. Exiting.")
        return None