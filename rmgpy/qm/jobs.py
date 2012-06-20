
"""
QMTP Jobs is a module that collects all classes that execute external programs in the 
framework of on-the-fly generation of thermodynamic properties via quantum chemical software packages.

A superclass, QMJob, is created with two methods: .run() and .check(output):
.run() is the method that calls the external program
.check(output) checks the output produced by the .run() method for specific lines that indicate 
that a succesful termination of the external tool has occurred.

The class has several attributes:
-name:  the filename of the file containing all the settings (QM method, 3D coords for molecule, etc...)
-directory
-inputFileExtension
-outputFileExtension: the extension of the file produced by the external tool

The following implementations of QMJob are created:

-MOPACJob
    Calls MOPAC2009 from the command line

"""

import os
from subprocess import Popen, PIPE

import logging
import symmetry
from qmverifier import QMVerifier

class QMJob:
    """
    supertype for all wrapper classes of 
    third-party quantum chemistry packages such as 
    OpenMopac, Gaussian 03, MM4.
    
    New subclasses of QMJob need to specify
    -the name of the executable
    -the command that is given to the executable
    -the name of the molecule (and file)
    -the directory in which the input file can be found
    -the input file extension which will be appended to the name and used in the command
    -the output file extension which will be used to check the output file for the correct termination.
    """
    
    "The path to the executable."
    executablePath = ''
    inputFileExtension = ''
    outputFileExtension = ''
    
    def __init__(self, molfile):
        "The molfile associated with this job."
        # does everything have to be done via molfiles? do we even mean molfile?
        self.molfile = molfile

        "The command line command."
        self.command = ''
 
    def run(self):
        """
        run the quantum chemistry package by initiating 
        a Process and pointing to the external executable
        """ 
        return False
    

    def check(self, output):
        """
        Check whether the Process has terminated succesfully
        
        Returns boolean flag
        """
        return False
    
class MOPACJob(QMJob):
    """
    MOPACJob is the wrapper class for the Open Mopac executable.
    
    The executable is in a directory stored in an environment variable $MOPAC_DIR and is named
    MOPAC2009.exe (.exe also for *nix!)
    
    The input file is of extension .mop and comprises of :
    -a top section with instructions
    -a section with 3D coordinates
    -a bottom section with instructions
    
    The output file is of extension .out.
    
    """
    
    inputFileExtension = '.mop'
    outputFileExtension = '.out'
    executablePath = os.path.join(os.getenv('MOPAC_DIR') , 'MOPAC2009.exe')
    
    def __init__(self, molfile):
        QMJob.__init__(self, molfile)
        assert os.path.exists(self.executablePath), "Please set the environment variable MOPAC_DIR to the directory containing MOPAC2009.exe"

        'specify the input file'
        self.command = os.path.join(self.molfile.directory, self.molfile.name + self.inputFileExtension)
        
    def check(self):
        verifier = QMVerifier(self.molfile)
        return verifier.verifyNoFailure()
    
    def run(self):
        process = Popen([self.executablePath, self.command])
        process.communicate()# necessary to wait for executable termination!

        return self.check()

class SymmetryJob(QMJob):
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
        QMJob.__init__(self, molfile)

        'the command line command'
        self.command = []
        
        "IQMData is the object that holds information from a previous QM Job on 3D coords, molecule etc..."
        self.qmdata = iqmdata
        
        self.inputFile = self.molfile.name + self.inputFileExtension
        
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
        with open(os.path.join(self.molfile.directory, self.inputFile), 'w') as input_file:
            input_file.write(geom)
        input_file.close()
        logging.info("Symmetry input file written to %s"%os.path.join(self.molfile.directory, self.inputFile))
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
            self.command.append(os.path.join(self.molfile.directory, self.inputFile))
            
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

