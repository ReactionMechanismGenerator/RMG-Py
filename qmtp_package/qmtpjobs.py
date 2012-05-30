
'''
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

    
'''

import subprocess
import logging
import platform
import os
import symmetry as sym
import qmverifier as verif
import qmtp as qm

class QMJob:
    '''
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
    '''
    
    def __init__(self, molfile):
        self.molfile = molfile
        self.inputFileExtension = ''
        self.outputFileExtension = ''
        
        'the keywords denoting the executable'
        self.executable = ''
        
        'the command line command'
        self.command = ''
 
    def run(self):
        '''
        run the quantum chemistry package by initiating 
        a Process and pointing to the external executable
        ''' 
        return -1
    

    def check(self, output):
        '''
        Check whether the Process has terminated succesfully
        
        Returns boolean flag
        '''
        return -1
    
class MOPACJob(QMJob):
    '''
    MOPACJob is the wrapper class for the Open Mopac executable.
    
    The executable is in a directory stored in an environment variable $MOPAC_DIR and is named
    MOPAC2009.exe (.exe also for *nix!)
    
    The input file is of extension .mop and comprises of :
    -a top section with instructions
    -a section with 3D coordinates
    -a bottom section with instructions
    
    The output file is of extension .out.
    
    
    ''' 
    def __init__(self, molfile):
        QMJob.__init__(self, molfile)
        
        self.inputFileExtension = '.mop'
        self.outputFileExtension = '.out'
        
        '''
        TODO maybe it's better to call the alias 'mopac'. However, this did not work yet... 
        '''
        self.executable = os.path.join(os.getenv('MOPAC_DIR') , 'MOPAC2009.exe')#assumes this env var is pointing to install directory of mopac!
                
        'specify the input file'
        self.command = os.path.join(self.molfile.directory ,self.molfile.name+ self.inputFileExtension)
        
    def check(self):
        verifier = verif.QMVerifier(self.molfile)
        return verifier.verifyNoFailure()
         
    
    def run(self):
        process = subprocess.Popen([self.executable, self.command])
        process.communicate()# necessary to wait for executable termination!


        return self.check()

class SymmetryJob(QMJob):
     '''
     determine the point group using the SYMMETRY program (http://www.cobalt.chem.ucalgary.ca/ps/symmetry/)
     required input is a line with number of atoms followed by lines for each atom including atom number and x,y,z coordinates
     finalTol determines how loose the point group criteria are; values are comparable to those specifed in the GaussView point group interface
     public String determinePointGroupUsingSYMMETRYProgram(String geom, double finalTol){ 
     '''   
    
     maxAttemptNumber = 4;
         
     def __init__(self, molfile, iqmdata, environ = os.environ.get("RMG_workingDirectory")):
        QMJob.__init__(self, molfile)
        
        'the keywords denoting the executable'
        self.executable = ''
        
        'the command line command'
        self.command = []
        
        
        '''
        IQMData is the object that holds information from a previous QM Job on 3D coords, molecule etc...
        '''
        self.qmdata = iqmdata
        
        self.inputFile = self.molfile.name + '.symm'
        
        self.environ = environ
        
        self.attemptNumber = 1
        
        self.pointGroupFound = False
     
     def check(self, output):
        result = "";
        output = output.split('\n')
        try:
            #check for errors and display the error if there is one
            for line in output:
                if line.startswith("It seems to be the "):#last line, ("It seems to be the [x] point group") indicates point group
                    lineArray = line.split(" ")#split the line around spaces
                    result = lineArray[5]#point group string should be the 6th word
            
        except Exception as e:
            err = "Error in veryfying point group calculation process using SYMMETRY \n";
            err += e.toString();
            logging.error(err+str(e))

        logging.info("Point group: "+ result)#print result, at least for debugging purposes
        return result;
           
     def run(self):
        try:
            pp = subprocess.Popen(self.command, stdout=subprocess.PIPE)
            stdout, stderr = pp.communicate()
        except Exception as e:
            logging.error('Error in running Symmetry process \n')
            logging.error(str(e))
        
        return self.check(stdout)    
    
     def writeInputFile(self):
          geom = str(self.qmdata.numberOfAtoms) + "\n"
          for i in range(self.qmdata.numberOfAtoms):
            geom = geom +\
            str(self.qmdata.atomicNumbers[i]) +\
            " "+ str(self.qmdata.atomCoords[i][0]) +\
            " " +str(self.qmdata.atomCoords[i][1]) +\
            " " +str(self.qmdata.atomCoords[i][2]) +"\n"
            
          '''
          Write the input file for the SYMMETRY program based on the passed-in string.
          '''
          with open(os.path.join(self.molfile.directory, self.inputFile), 'w') as input_file:
               input_file.write(geom)
          input_file.close()
          logging.info('Symmetry input file written to %s'%os.path.join(self.molfile.directory, self.inputFile))
          return input_file     
        
     def calculate(self):
            
        self.writeInputFile();

        result = "";

        #continue trying to generate symmetry group until too many no. of attempts or until a point group is found: 
        while self.attemptNumber <= SymmetryJob.maxAttemptNumber and not self.pointGroupFound:
            result = "";
            '''
            TODO only *nix case works!
            '''
            path_symm = os.path.join(self.environ, 'bin/symmetry')
            print 'Symmetry exec: '+str(os.path.exists(path_symm))
            #path to the symmetry input file 'InChIKey.symm'
            path_mol = os.path.join(self.molfile.directory, self.inputFile)
            print 'mol symm file: '+str(os.path.exists(path_mol))
            if self.attemptNumber == 1:
                self.command = [path_symm  , '-final', '0.02',path_mol]
            elif self.attemptNumber == 2:
                self.command = [path_symm  , '-final', '0.1',path_mol]
            elif self.attemptNumber == 3:
                self.command = [path_symm  , '-primary', '0.2', '-final' ,'0.1', path_mol]
            elif self.attemptNumber == 4:
                logging.warning("*****WARNING****: Using last-resort symmetry estimation options; symmetry may be underestimated")
                self.command = [path_symm  , '-final', '0.0',path_mol]
            else:
                logging.critical("Invalid attemptNumber: "+ self.attemptNumber)
            
            #call the program and read the result
            result = self.run();

            #check for a recognized point group
            symmPGC = sym.PointGroupDictionary()
            #symmPGC.initiate()
                        
            if symmPGC.contains(result):
                self.pointGroupFound = True;
            else:
                if self.attemptNumber < self.maxAttemptNumber:
                    logging.info("Attempt number "+self.attemptNumber+" did not identify a recognized point group (" +result+"). Will retry with looser point group criteria.")
                else:
                    logging.critical("Final attempt number "+self.attemptNumber+" did not identify a recognized point group (" +result+"). Exiting.")

                self.attemptNumber = self.attemptNumber + 1
                 
            return sym.PointGroup(result)
        