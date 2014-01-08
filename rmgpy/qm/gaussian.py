import os

import re
import external.cclib as cclib
import logging
import time
import math
import numpy
from subprocess import Popen, PIPE

from rmgpy.molecule import Molecule, Atom, getElement
from rmgpy.reaction import Reaction
from qmdata import CCLibData
from molecule import QMMolecule
from reaction import QMReaction
from rmgpy.data.base import Entry
from rmgpy.data.kinetics import saveEntry
from rmgpy.data.kinetics.transitionstates import DistanceData

class Gaussian:
    """
    A base class for all QM calculations that use Gaussian.
    
    Classes such as :class:`GaussianMol` will inherit from this class.
    """
    
    inputFileExtension = '.gjf'
    outputFileExtension = '.log'
    
    gaussEnv = os.getenv('GAUSS_EXEDIR') or os.getenv('g09root') or os.getenv('g03root') or ""
    if os.path.exists(os.path.join(gaussEnv , 'g09')):
        executablePath = os.path.join(gaussEnv , 'g09')
    elif os.path.exists(os.path.join(gaussEnv , 'g03')):
        executablePath = os.path.join(gaussEnv , 'g03')
    else:
        executablePath = os.path.join(gaussEnv , '(g03 or g09)')

    usePolar = False
    
    #: List of phrases that indicate failure
    #: NONE of these must be present in a succesful job.
    failureKeys = [
                   'ERROR TERMINATION',
                   'IMAGINARY FREQUENCIES'
                   ]
    
    #: List of phrases to indicate success.
    #: ALL of these must be present in a successful job.
    successKeys = [
                   'Normal termination of Gaussian'
                  ]
 
    def testReady(self):
        if not os.path.exists(self.executablePath):
            raise Exception("Couldn't find Gaussian executable at {0}. Try setting your GAUSS_EXEDIR environment variable.".format(self.executablePath))

   
    def run(self):
        self.testReady()
        # submits the input file to Gaussian
        process = Popen([self.executablePath, self.inputFilePath, self.outputFilePath])
        process.communicate()# necessary to wait for executable termination!
        
        return self.verifyOutputFile()
        
    def verifyOutputFile(self):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful GAUSSIAN simulation already exists for the molecule with the 
        given (augmented) InChI Key.
        
        The definition of finding a successful simulation is based on these criteria:
        1) finding an output file with the file name equal to the InChI Key
        2) NOT finding any of the keywords that are denote a calculation failure
        3) finding all the keywords that denote a calculation success.
        4) finding a match between the InChI of the given molecule and the InchI found in the calculation files
        
        If any of the above criteria is not matched, False will be returned and the procedures to start a new calculation 
        will be initiated.
        """
        if not os.path.exists(self.outputFilePath):
            logging.info("Output file {0} does not exist.".format(self.outputFilePath))
            return False
    
        InChIMatch=False #flag (1 or 0) indicating whether the InChI in the file matches InChIaug this can only be 1 if InChIFound is also 1
        InChIFound=False #flag (1 or 0) indicating whether an InChI was found in the log file
        
        # Initialize dictionary with "False"s 
        successKeysFound = dict([(key, False) for key in self.successKeys])
        
        with open(self.outputFilePath) as outputFile:
            for line in outputFile:
                line = line.strip()
                
                for element in self.failureKeys: #search for failure keywords
                    if element in line:
                        logging.error("Gaussian output file contains the following error: {0}".format(element) )
                        return False
                    
                for element in self.successKeys: #search for success keywords
                    if element in line:
                        successKeysFound[element] = True
               
                if line.startswith("InChI="):
                    logFileInChI = line #output files should take up to 240 characters of the name in the input file
                    InChIFound = True
                    if logFileInChI == self.uniqueIDlong:
                        InChIMatch = True
                    elif self.uniqueIDlong.startswith(logFileInChI):
                        logging.info("InChI too long to check, but beginning matches so assuming OK.")
                        InChIMatch = True
                    else:
                        logging.info("InChI in log file ({0}) didn't match that in geometry ({1}).".format(logFileInChI, self.geometry.uniqueIDlong))
        
        # Check that ALL 'success' keywords were found in the file.
        if not all( successKeysFound.values() ):
            logging.error('Not all of the required keywords for sucess were found in the output file!')
            return False
        
        if not InChIFound:
            logging.error("No InChI was found in the Gaussian output file {0}".format(self.outputFilePath))
            return False
        
        if InChIMatch:
            logging.info("Successful Gaussian quantum result found in {0}".format(self.outputFilePath))
            # " + self.molfile.name + " ("+self.molfile.InChIAug+") has been found. This log file will be used.")
            return True
        else:
            return False # until the next line works
        
        #InChIs do not match (most likely due to limited name length mirrored in log file (240 characters), but possibly due to a collision)
        return self.checkForInChiKeyCollision(logFileInChI) # Not yet implemented!
    
    
class GaussianMol(QMMolecule, Gaussian):
    """
    A base Class for calculations of molecules using Gaussian. 
    
    Inherits from both :class:`QMMolecule` and :class:`Gaussian`.
    """
    
    def writeInputFile(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        molfile = self.getMolFilePathForCalculation(attempt) 
        atomline = re.compile('\s*([\- ][0-9.]+\s+[\-0-9.]+\s+[\-0-9.]+)\s+([A-Za-z]+)')
        
        output = ['', self.geometry.uniqueIDlong, '' ]
        output.append("{charge}   {mult}".format(charge=0, mult=(self.molecule.getRadicalCount() + 1) ))
        
        atomCount = 0
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    output.append("{0:8s} {1}".format(match.group(2), match.group(1)))
                    atomCount += 1
        assert atomCount == len(self.molecule.atoms)
    
        output.append('')
        input_string = '\n'.join(output)
        
        top_keys = self.inputFileKeywords(attempt)
        with open(self.inputFilePath, 'w') as gaussianFile:
            gaussianFile.write(top_keys)
            gaussianFile.write('\n')
            gaussianFile.write(input_string)
            gaussianFile.write('\n')
            if self.usePolar:
                gaussianFile.write('\n\n\n')
                gaussianFile.write(polar_keys)
    
    def inputFileKeywords(self, attempt):
        """
        Return the top keywords.
        """
        raise NotImplementedError("Should be defined by subclass, eg. GaussianMolPM3")
    
    def generateQMData(self):
        """
        Calculate the QM data and return a QMData object, or None if it fails.
        """
        for atom in self.molecule.vertices:
            if atom.atomType.label == 'N5s' or atom.atomType.label == 'N5d' or atom.atomType.label =='N5dd' or atom.atomType.label == 'N5t' or atom.atomType.label == 'N5b':
                return None
                
        if self.verifyOutputFile():
            logging.info("Found a successful output file already; using that.")
            source = "QM Gaussian result file found from previous run."
        else:
            self.createGeometry()
            success = False
            for attempt in range(1, self.maxAttempts+1):
                self.writeInputFile(attempt)
                logging.info('Trying {3} attempt {0} of {1} on molecule {2}.'.format(attempt, self.maxAttempts, self.molecule.toSMILES(), self.__class__.__name__))
                success = self.run()
                if success:
                    logging.info('Attempt {0} of {1} on species {2} succeeded.'.format(attempt, self.maxAttempts, self.molecule.toAugmentedInChI()))
                    source = "QM Gaussian result created during this run."
                    break
            else:
                logging.error('QM thermo calculation failed for {0}.'.format(self.molecule.toAugmentedInChI()))
                return None
        result = self.parse() # parsed in cclib
        result.source = source
        return result # a CCLibData object
    
    def parse(self):
        """
        Parses the results of the Gaussian calculation, and returns a CCLibData object.
        """
        parser = cclib.parser.Gaussian(self.outputFilePath)
        parser.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
        cclibData = parser.parse()
        radicalNumber = sum([i.radicalElectrons for i in self.molecule.atoms])
        qmData = CCLibData(cclibData, radicalNumber+1)
        return qmData


class GaussianMolPM3(GaussianMol):

    #: Keywords that will be added at the top of the qm input file
    keywords = [
               "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3)",
               "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)",
               "# pm3 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm" ,
               "# pm3 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm",
               "# pm3 opt=(verytight,gdiis,small) freq IOP(2/16=3)",
               "# pm3 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)",
               "# pm3 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)",
               "# pm3 opt=tight freq IOP(2/16=3)",
               "# pm3 opt=tight freq=numerical IOP(2/16=3)",
               "# pm3 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)",
               "# pm3 opt freq IOP(2/16=3)",
               "# pm3 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)",
               "# pm3 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm",
               "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
               "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
               "# pm3 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm",
               "# pm3 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm",
               "# pm3 opt=(calcall,small,maxcyc=100) IOP(2/16=3)",
               ]

    @property
    def scriptAttempts(self):
        "The number of attempts with different script keywords"
        return len(self.keywords)

    @property
    def maxAttempts(self):
        "The total number of attempts to try"
        return 2 * len(self.keywords)

    def inputFileKeywords(self, attempt):
        """
        Return the top keywords for attempt number `attempt`.

        NB. `attempt`s begin at 1, not 0.
        """
        assert attempt <= self.maxAttempts
        if attempt > self.scriptAttempts:
            attempt -= self.scriptAttempts
        return self.keywords[attempt-1]



class GaussianMolPM6(GaussianMol):

    #: Keywords that will be added at the top of the qm input file
    keywords = [
               "# pm6 opt=(verytight,gdiis) freq IOP(2/16=3)",
               "# pm6 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)",
               "# pm6 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm" ,
               "# pm6 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm",
               "# pm6 opt=(verytight,gdiis,small) freq IOP(2/16=3)",
               "# pm6 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)",
               "# pm6 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)",
               "# pm6 opt=tight freq IOP(2/16=3)",
               "# pm6 opt=tight freq=numerical IOP(2/16=3)",
               "# pm6 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)",
               "# pm6 opt freq IOP(2/16=3)",
               "# pm6 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)",
               "# pm6 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm",
               "# pm6 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
               "# pm6 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
               "# pm6 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm",
               "# pm6 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm",
               "# pm6 opt=(calcall,small,maxcyc=100) IOP(2/16=3)",
               ]

    @property
    def scriptAttempts(self):
        "The number of attempts with different script keywords"
        return len(self.keywords)

    @property
    def maxAttempts(self):
        "The total number of attempts to try"
        return 2 * len(self.keywords)

    def inputFileKeywords(self, attempt):
        """
        Return the top keywords for attempt number `attempt`.

        NB. `attempt`s begin at 1, not 0.
        """
        assert attempt <= self.maxAttempts
        if attempt > self.scriptAttempts:
            attempt -= self.scriptAttempts
        return self.keywords[attempt-1]
        
class GaussianMolMP2(GaussianMol):

    #: Keywords that will be added at the top of the qm input file
    keywords = [
               "# mp2 opt=(verytight,gdiis) freq IOP(2/16=3)",
               "# mp2 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)",
               "# mp2 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm" ,
               "# mp2 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm",
               "# mp2 opt=(verytight,gdiis,small) freq IOP(2/16=3)",
               "# mp2 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)",
               "# mp2 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)",
               "# mp2 opt=tight freq IOP(2/16=3)",
               "# mp2 opt=tight freq=numerical IOP(2/16=3)",
               "# mp2 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)",
               "# mp2 opt freq IOP(2/16=3)",
               "# mp2 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)",
               "# mp2 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm",
               "# mp2 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
               "# mp2 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
               "# mp2 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm",
               "# mp2 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm",
               "# mp2 opt=(calcall,small,maxcyc=100) IOP(2/16=3)",
               ]

    @property
    def scriptAttempts(self):
        "The number of attempts with different script keywords"
        return len(self.keywords)

    @property
    def maxAttempts(self):
        "The total number of attempts to try"
        return 2 * len(self.keywords)

    def inputFileKeywords(self, attempt):
        """
        Return the top keywords for attempt number `attempt`.

        NB. `attempt`s begin at 1, not 0.
        """
        assert attempt <= self.maxAttempts
        if attempt > self.scriptAttempts:
            attempt -= self.scriptAttempts
        return self.keywords[attempt-1]

class GaussianMolB3LYP(GaussianMol):
    """
    For use with the automated transition state estimator. This will find the
    stable species geomtries when required for TST rate calculation.
    """

    #: Keywords that will be added at the top of the qm input file
    # removed 'gdiis' from the keywords; http://www.gaussian.com/g_tech/g_ur/d_obsolete.htm
    keywords = [
               "# b3lyp/6-31+g(d,p) opt=(verytight) freq IOP(2/16=3)",
               "# b3lyp/6-31+g(d,p) opt=(verytight) freq IOP(2/16=3) IOP(4/21=2)",
               "# b3lyp/6-31+g(d,p) opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm" ,
               "# b3lyp/6-31+g(d,p) opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm",
               "# b3lyp/6-31+g(d,p) opt=(verytight,small) freq IOP(2/16=3)",
               "# b3lyp/6-31+g(d,p) opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)",
               "# b3lyp/6-31+g(d,p) opt=(verytight,maxcyc=200) freq=numerical IOP(2/16=3)",
               "# b3lyp/6-31+g(d,p) opt=tight freq IOP(2/16=3)",
               "# b3lyp/6-31+g(d,p) opt=tight freq=numerical IOP(2/16=3)",
               "# b3lyp/6-31+g(d,p) opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)",
               "# b3lyp/6-31+g(d,p) opt freq IOP(2/16=3)",
               "# b3lyp/6-31+g(d,p) opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)",
               "# b3lyp/6-31+g(d,p) opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm",
               "# b3lyp/6-31+g(d,p) opt=(tight,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
               "# b3lyp/6-31+g(d,p) opt=(tight,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
               "# b3lyp/6-31+g(d,p) opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm",
               "# b3lyp/6-31+g(d,p) opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm",
               "# b3lyp/6-31+g(d,p) opt=(calcall,small,maxcyc=100) IOP(2/16=3)",
               ]

    @property
    def scriptAttempts(self):
        "The number of attempts with different script keywords"
        return len(self.keywords)

    @property
    def maxAttempts(self):
        "The total number of attempts to try"
        return 2 * len(self.keywords)

    def inputFileKeywords(self, attempt):
        """
        Return the top keywords for attempt number `attempt`.

        NB. `attempt`s begin at 1, not 0.
        """
        assert attempt <= self.maxAttempts
        if attempt > self.scriptAttempts:
            attempt -= self.scriptAttempts
        return self.keywords[attempt-1]

##########################################################################################

class GaussianTS(QMReaction, Gaussian):
    """
    A base Class for calculations of transition states using Gaussian. 

    Inherits from both :class:`QMReaction` and :class:`Gaussian`.
    """
    
    #: List of phrases to indicate success.
    #: ALL of these must be present in a successful job.
    successKeys = [
                   'Normal termination of Gaussian',
                   '******    1 imaginary frequencies (negative Signs) ******',
                  ]
    
    #: List of phrases that indicate failure
    #: NONE of these must be present in a succesful job.
    failureKeys = [
                   'ERROR TERMINATION',
                   'Error in internal coordinate system.',
                   ]
    
    def writeInputFile(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        numProc = '%nprocshared=' + '11' + '\n' # could be something that is set in the qmSettings
        mem = '%mem=' + '1GB' + '\n' # could be something that is set in the qmSettings
        chk_file = '%chk=' + os.path.join(self.settings.fileStore, self.uniqueID) + '\n'
        
        molfile = self.geometry.getRefinedMolFilePath()
        atomline = re.compile('\s*([\- ][0-9.]+\s+[\-0-9.]+\s+[\-0-9.]+)\s+([A-Za-z]+)')
        
        output = ['', self.geometry.uniqueID, '' ]
        output.append("{charge}   {mult}".format(charge=0, mult=(self.geometry.molecule.getRadicalCount() + 1) ))
        
        atomCount = 0
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    output.append("{0:8s} {1}".format(match.group(2), match.group(1)))
                    atomCount += 1
        assert atomCount == len(self.geometry.molecule.atoms)
        
        output.append('')
        input_string = '\n'.join(output) + '\n'
        top_keys = self.keywords[attempt - 1] + '\n'
        
        atomTypes = []
        for atom in self.geometry.molecule.atoms:
            if not atom.element.symbol in atomTypes:
                atomTypes.append(atom.element.symbol)
        
        with open(self.inputFilePath, 'w') as gaussianFile:
            gaussianFile.write(numProc)
            gaussianFile.write(mem)
            gaussianFile.write(chk_file)
            gaussianFile.write(top_keys)
            if attempt == 1:
                gaussianFile.write(input_string)
            else:
                gaussianFile.write('\n')
            # for atom in atomTypes:
            #     gaussianFile.write(self.mg3s[atom])
            gaussianFile.write('\n')
    
    def writeIRCFile(self):
        """
        Using the :class:`Geometry` object, write the input file for the 
        IRC calculation on the transition state. The geometry is taken 
        from the checkpoint file created during the geometry search.
        """
        
        numProc = '%nprocshared=' + '11' + '\n' # could be something that is set in the qmSettings
        mem = '%mem=' + '1GB' + '\n' # could be something that is set in the qmSettings
        chk_file = '%chk=' + os.path.join(self.settings.fileStore, self.uniqueID) + '\n'
        top_keys = self.keywords[4] + '\n\n'
        output = "{charge}   {mult}".format(charge=0, mult=(self.geometry.molecule.getRadicalCount() + 1) )
        
        atomTypes = []
        for atom in self.geometry.molecule.atoms:
            if not atom.element.symbol in atomTypes:
                atomTypes.append(atom.element.symbol)
        
        with open(self.ircInputFilePath, 'w') as gaussianFile:
            gaussianFile.write(numProc)
            gaussianFile.write(mem)
            gaussianFile.write(chk_file)
            gaussianFile.write(top_keys)
            gaussianFile.write(output)
            gaussianFile.write('\n')
            # for atom in atomTypes:
            #     gaussianFile.write(self.mg3s[atom]) 
            gaussianFile.write('\n')
            
    def inputFileKeywords(self, attempt):
        """
        Return the top keywords.
        """
        raise NotImplementedError("Should be defined by subclass, eg. GaussianTSM062X")

    def generateQMKinetics(self):
        """
        Calculate the QM data and return a QMData object.
        """
        self.createGeometry()
        if self.verifyOutputFile():
            logging.info("Found a successful output file already; using that.")
        else:
            success = False
            for attempt in range(1, self.maxAttempts+1):
                self.writeInputFile(attempt)
                self.writeIRCFile()
                success = self.run()
                if success:
                    logging.info('Attempt {0} of {1} on species {2} succeeded.'.format(attempt, self.maxAttempts, self.molecule.toAugmentedInChI()))
                    break
            else:
                logging.error('QM thermo calculation failed for {0}.'.format(self.molecule.toAugmentedInChI()))
                return None
        result = self.parse() # parsed in cclib
        return result
    
    def runIRC(self):
        self.testReady()
        # submits the input file to Gaussian
        process = Popen([self.executablePath, self.ircInputFilePath, self.ircOutputFilePath])
        process.communicate()# necessary to wait for executable termination!
        
        return self.verifyIRCOutputFile()
    
    def verifyOutputFile(self):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful GAUSSIAN simulation already exists for the molecule with the 
        given file name.
        
        The definition of finding a successful simulation is based on these criteria:
        1) finding an output file with the file name equal to the the reaction unique ID
        2) NOT finding any of the keywords that are denote a calculation failure
        3) finding all the keywords that denote a calculation success.
        
        If any of the above criteria is not matched, False will be returned and the procedures to start a new calculation 
        will be initiated. The second boolean flag indicates if there was a failure in the internal coordinate system.
        This will initiate a subsequent calculation in cartesian coordinates.
        """
        if not os.path.exists(self.outputFilePath):
            logging.info("Output file {0} does not exist.".format(self.outputFilePath))
            return False
        
        # Initialize dictionary with "False"s 
        successKeysFound = dict([(key, False) for key in self.successKeys])
        failureKeysFound = dict([(key, False) for key in self.failureKeys])
        
        with open(self.outputFilePath) as outputFile:
            for line in outputFile:
                line = line.strip()
                
                for element in self.failureKeys: #search for failure keywords
                    if element in line:
                        logging.error("Gaussian output file contains the following error: {0}".format(element) )
                        failureKeysFound[element] = True
                    
                for element in self.successKeys: #search for success keywords
                    if element in line:
                        successKeysFound[element] = True
        
        if any(failureKeysFound.values()):
            if failureKeysFound['Error in internal coordinate system.']:
                return False, True
            else:
                return False, False
        
        # Check that ALL 'success' keywords were found in the file.
        if not all( successKeysFound.values() ):
            logging.error('Not all of the required keywords for success were found in the output file!')
            return False, False
        else:
            return True, False
    
    def verifyIRCOutputFile(self):
        """
        Check's that the resulting geometries of the path analysis match the reaction.
        """
        
        """
        Compares IRC geometries to input geometries.
        """
        if not os.path.exists(self.ircOutputFilePath):
            logging.info("Output file {0} does not exist.".format(self.ircOutputFilePath))
            return False
        
        # Initialize dictionary with "False"s 
        successKeysFound = dict([(key, False) for key in self.successKeys])
        
        pth1 = list()
        steps = list()
        with open(self.ircOutputFilePath) as outputFile:
            for line in outputFile:
                line = line.strip()
                
                for element in self.failureKeys: #search for failure keywords
                    if element in line:
                        logging.error("Gaussian IRC output file contains the following error: {0}".format(element) )
                        return False
                    
                for element in self.successKeys: #search for success keywords
                    if element in line:
                        successKeysFound[element] = True
                
                if line.startswith('Point Number:'):
                    if int(line.split()[2]) > 0:
                        if int(line.split()[-1]) == 1:
                            ptNum = int(line.split()[2])
                            pth1.append(ptNum)
                        else:
                            pass
                elif line.startswith('# OF STEPS ='):
                    numStp = int(line.split()[-1])
                    steps.append(numStp)
        
        # Check that ALL 'success' keywords were found in the file.
        if not successKeysFound['Normal termination of Gaussian']:
            logging.error('Not all of the required keywords for success were found in the IRC output file!')
            return False
        # This indexes the coordinate to be used from the parsing
        elif steps == []:
            logging.error('No steps taken in the IRC calculation!')
            return False
        else:
            pth1End = sum(steps[:pth1[-1]])		
            # Compare the reactants and products
            ircParse = cclib.parser.Gaussian(self.ircOutputFilePath)
            ircParse.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
            ircParse = ircParse.parse()
        
            atomnos = ircParse.atomnos
            atomcoords = ircParse.atomcoords
        
            # Convert the IRC geometries into RMG molecules
            # We don't know which is reactant or product, so take the two at the end of the
            # paths and compare to the reactants and products
            mol1 = Molecule()
            mol1.fromXYZ(atomnos, atomcoords[pth1End])
            mol2 = Molecule()
            mol2.fromXYZ(atomnos, atomcoords[-1])
            
            targetReaction = Reaction(
                                    reactants = [reactant.toSingleBonds() for reactant in self.reaction.reactants],
                                    products = [product.toSingleBonds() for product in self.reaction.products],
                                    )
            testReaction = Reaction(
                                    reactants = mol1.split(),
                                    products = mol2.split(),                     
                                    )
                                                            
            if targetReaction.isIsomorphic(testReaction):
                return True
            else:
                return False
    
    def parse(self):
        """
        Parses the results of the Gaussian calculation, and returns a CCLibData object.
        """
        parser = cclib.parser.Gaussian(self.outputFilePath)
        parser.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
        cclibData = parser.parse()
        radicalNumber = 0
        for molecule in self.reaction.reactants:
            radicalNumber += sum([i.radicalElectrons for i in molecule.atoms])
        self.qmData = CCLibData(cclibData, radicalNumber+1)
        return self.qmData
    
    def parseTS(self, labels):
        
        def getDistance(coordinates1, coordinates2):
            """
            Return the square of the distance (in Angstrom) between the two atoms.
            """
            diff = (coordinates1.coords - coordinates2.coords)
            return math.sqrt(sum(diff * diff))
        
        self.parse()
        atomCoords = self.qmData.atomCoords.getValue()
        atomNums = self.qmData.atomicNumbers
        
        atom1 = Atom(element=getElement(int(atomNums[labels[0]])), coords=atomCoords[labels[0]])
        atom2 = Atom(element=getElement(int(atomNums[labels[1]])), coords=atomCoords[labels[1]])
        atom3 = Atom(element=getElement(int(atomNums[labels[2]])), coords=atomCoords[labels[2]])
        
        at12 = getDistance(atom1, atom2)
        at23 = getDistance(atom2, atom3)
        at13 = getDistance(atom1, atom3)
        
        atomDist = [at12, at23, at13]
        
        return atomDist
    
    def writeRxnOutputFile(self, labels):
        
        product = self.reaction.products[0].merge(self.reaction.products[1])
        star3 = product.getLabeledAtom('*1').sortingLabel
        star1 = product.getLabeledAtom('*3').sortingLabel
        product.atoms[star1].label = '*1'
        product.atoms[star3].label = '*3'
        
        atomDist = self.parseTS(labels)
        
        distances = {'d12':float(atomDist[0]), 'd23':float(atomDist[1]), 'd13':float(atomDist[2])}
        user = "Pierre Bhoorasingh <bhoorasingh.p@husky.neu.edu>"
        description = "Found via group estimation strategy using automatic transition state generator"
        entry = Entry(
            index = 1,
            item = self.reaction,
            data = DistanceData(distances=distances, method='B3LYP/6-31+G(d,p)'),
            shortDesc = "B3LYP/6-31+G(d,p) calculation via group additive TS generator.",
            history = [(time.asctime(), user, 'action', description)]
        )
        
        outputDataFile = os.path.join(self.file_store_path, self.uniqueID + '.data')
        with open(outputDataFile, 'w') as parseFile:
            saveEntry(parseFile, entry)
        
class GaussianTSM062X(GaussianTS):
    """
    M06-2X requires a minimally augmented basis set for good prediction of interatomic distances.
    These required basis sets are not available in the standard G09 kit, so the MG3S has been implemented here for use.
    Recommended by Xu et al. J. Chem. Theory Comput. (2011)
    """
  

    #: Keywords that will be added at the top of the qm input file
    keywords = [
                "# m062x/gen opt=(ts,calcall,tight,noeigentest)  int=ultrafine nosymm",
                "# m062x/gen opt=(ts,calcall,tight,noeigentest,cartesian) int=ultrafine geom=allcheck guess=check nosymm",
                "# m062x/gen opt=(ts,calcall,noeigentest,cartesian) nosymm geom=allcheck guess=check nosymm",
                "# m062x/gen opt=(ts,calcall,noeigentest) nosymm",
                "# m062x/gen irc=(calcall,report=read) geom=allcheck guess=check nosymm",
               ]
    """
    This needs some work, to determine options that are best used. Commented out the
    methods for now.
    """
    
    mg3s = {
            'H'  : """H  0 \nS  3 1.00\n  3.386500000000E+01 2.549380000000E-02\n  5.094790000000E+00 1.903730000000E-01\n  1.158790000000E+00 8.521610000000E-01\nS  1 1.00\n  3.258400000000E-01 1.000000000000E+00\nS  1 1.00\n  1.027410000000E-01 1.000000000000E+00\nP  1 1.00\n  1.500000000000E+00 1.000000000000E+00\nP  1 1.00\n  3.750000000000E-01 1.000000000000E+00\n****\n""",
            'He' : """He 0 \nS  3 1.00\n  9.812430000000E+01 2.874520000000E-02\n  1.476890000000E+01 2.080610000000E-01\n  3.318830000000E+00 8.376350000000E-01\nS  1 1.00\n  8.740470000000E-01 1.000000000000E+00\nS  1 1.00\n  2.445640000000E-01 1.000000000000E+00\nS  1 1.00\n  8.600000000000E-02 1.000000000000E+00\nP  1 1.00\n  1.500000000000E+00 1.000000000000E+00\nP  1 1.00\n  3.750000000000E-01 1.000000000000E+00\n****\n""",
            'Li' : """Li 0 \nS  6 1.00\n  9.004600000000E+02 2.287040000000E-03\n  1.344330000000E+02 1.763500000000E-02\n  3.043650000000E+01 8.734340000000E-02\n  8.626390000000E+00 2.809770000000E-01\n  2.483320000000E+00 6.587410000000E-01\n  3.031790000000E-01 1.187120000000E-01\nS  3 1.00\n  4.868900000000E+00 9.332930000000E-02\n  8.569240000000E-01 9.430450000000E-01\n  2.432270000000E-01 -2.798270000000E-01\nS  1 1.00\n  6.350700000000E-02 1.000000000000E+00\nS  1 1.00\n  2.436830000000E-02 1.000000000000E+00\nS  1 1.00\n  7.400000000000E-03 1.000000000000E+00\nP  3 1.00\n  4.868900000000E+00 3.276610000000E-02\n  8.569240000000E-01 1.597920000000E-01\n  2.432270000000E-01 8.856670000000E-01\nP  1 1.00\n  6.350700000000E-02 1.000000000000E+00\nP  1 1.00\n  2.436830000000E-02 1.000000000000E+00\nP  1 1.00\n  7.400000000000E-03 1.000000000000E+00\nD  1 1.00\n  4.000000000000E-01 1.000000000000E+00\nD  1 1.00\n  1.000000000000E-01 1.000000000000E+00\nF  1 1.00\n  1.500000000000E-01 1.000000000000E+00\n****\n""",
            'Be' : """Be 0 \nS  6 1.00\n  1.682800000000E+03 2.285740000000E-03\n  2.517150000000E+02 1.759380000000E-02\n  5.741160000000E+01 8.633150000000E-02\n  1.651710000000E+01 2.818350000000E-01\n  4.853640000000E+00 6.405940000000E-01\n  6.268630000000E-01 1.444670000000E-01\nS  3 1.00\n  8.093800000000E+00 1.086210000000E-01\n  1.740750000000E+00 9.273010000000E-01\n  4.858160000000E-01 -2.971690000000E-01\nS  1 1.00\n  1.636130000000E-01 1.000000000000E+00\nS  1 1.00\n  5.672850000000E-02 1.000000000000E+00\nS  1 1.00\n  2.070000000000E-02 1.000000000000E+00\nP  3 1.00\n  8.093800000000E+00 3.613440000000E-02\n  1.740750000000E+00 2.169580000000E-01\n  4.858160000000E-01 8.418390000000E-01\nP  1 1.00\n  1.636130000000E-01 1.000000000000E+00\nP  1 1.00\n  5.672850000000E-02 1.000000000000E+00\nP  1 1.00\n  2.070000000000E-02 1.000000000000E+00\nD  1 1.00\n  5.100000000000E-01 1.000000000000E+00\nD  1 1.00\n  1.275000000000E-01 1.000000000000E+00\nF  1 1.00\n  2.600000000000E-01 1.000000000000E+00\n****\n""",
            'B'  : """B  0 \nS  6 1.00\n  2.858890000000E+03 2.153750000000E-03\n  4.281400000000E+02 1.658230000000E-02\n  9.752820000000E+01 8.218700000000E-02\n  2.796930000000E+01 2.766180000000E-01\n  8.215770000000E+00 6.293160000000E-01\n  1.112780000000E+00 1.737700000000E-01\nS  3 1.00\n  1.324150000000E+01 1.174430000000E-01\n  3.001660000000E+00 9.180020000000E-01\n  9.128560000000E-01 -2.651050000000E-03\nS  1 1.00\n  3.154540000000E-01 1.000000000000E+00\nS  1 1.00\n  9.885630000000E-02 1.000000000000E+00\nS  1 1.00\n  3.150000000000E-02 1.000000000000E+00\nP  3 1.00\n  1.324150000000E+01 4.181000000000E-02\n  3.001660000000E+00 2.365750000000E-01\n  9.128560000000E-01 8.162140000000E-01\nP  1 1.00\n  3.154540000000E-01 1.000000000000E+00\nP  1 1.00\n  9.885630000000E-02 1.000000000000E+00\nP  1 1.00\n  3.150000000000E-02 1.000000000000E+00\nD  1 1.00\n  8.020000000000E-01 1.000000000000E+00\nD  1 1.00\n  2.005000000000E-01 1.000000000000E+00\nF  1 1.00\n  5.000000000000E-01 1.000000000000E+00\n****\n""",
            'C'  : """C  0 \nS  6 1.00\n  4.563240000000E+03 1.966650000000E-03\n  6.820240000000E+02 1.523060000000E-02\n  1.549730000000E+02 7.612690000000E-02\n  4.445530000000E+01 2.608010000000E-01\n  1.302900000000E+01 6.164620000000E-01\n  1.827730000000E+00 2.210060000000E-01\nS  3 1.00\n  2.096420000000E+01 1.146600000000E-01\n  4.803310000000E+00 9.199990000000E-01\n  1.459330000000E+00 -3.030680000000E-03\nS  1 1.00\n  4.834560000000E-01 1.000000000000E+00\nS  1 1.00\n  1.455850000000E-01 1.000000000000E+00\nS  1 1.00\n  4.380000000000E-02 1.000000000000E+00\nP  3 1.00\n  2.096420000000E+01 4.024870000000E-02\n  4.803310000000E+00 2.375940000000E-01\n  1.459330000000E+00 8.158540000000E-01\nP  1 1.00\n  4.834560000000E-01 1.000000000000E+00\nP  1 1.00\n  1.455850000000E-01 1.000000000000E+00\nP  1 1.00\n  4.380000000000E-02 1.000000000000E+00\nD  1 1.00\n  1.252000000000E+00 1.000000000000E+00\nD  1 1.00\n  3.130000000000E-01 1.000000000000E+00\nF  1 1.00\n  8.000000000000E-01 1.000000000000E+00\n****\n""",
            'N'  : """N  0 \nS  6 1.00\n  6.293480000000E+03 1.969790000000E-03\n  9.490440000000E+02 1.496130000000E-02\n  2.187760000000E+02 7.350060000000E-02\n  6.369160000000E+01 2.489370000000E-01\n  1.882820000000E+01 6.024600000000E-01\n  2.720230000000E+00 2.562020000000E-01\nS  3 1.00\n  3.063310000000E+01 1.119060000000E-01\n  7.026140000000E+00 9.216660000000E-01\n  2.112050000000E+00 -2.569190000000E-03\nS  1 1.00\n  6.840090000000E-01 1.000000000000E+00\nS  1 1.00\n  2.008780000000E-01 1.000000000000E+00\nS  1 1.00\n  6.390000000000E-02 1.000000000000E+00\nP  3 1.00\n  3.063310000000E+01 3.831190000000E-02\n  7.026140000000E+00 2.374030000000E-01\n  2.112050000000E+00 8.175920000000E-01\nP  1 1.00\n  6.840090000000E-01 1.000000000000E+00\nP  1 1.00\n  2.008780000000E-01 1.000000000000E+00\nP  1 1.00\n  6.390000000000E-02 1.000000000000E+00\nD  1 1.00\n  1.826000000000E+00 1.000000000000E+00\nD  1 1.00\n  4.565000000000E-01 1.000000000000E+00\nF  1 1.00\n  1.000000000000E+00 1.000000000000E+00\n****\n""",
            'O'  : """O  0 \nS  6 1.00\n  8.588500000000E+03 1.895150000000E-03\n  1.297230000000E+03 1.438590000000E-02\n  2.992960000000E+02 7.073200000000E-02\n  8.737710000000E+01 2.400010000000E-01\n  2.567890000000E+01 5.947970000000E-01\n  3.740040000000E+00 2.808020000000E-01\nS  3 1.00\n  4.211750000000E+01 1.138890000000E-01\n  9.628370000000E+00 9.208110000000E-01\n  2.853320000000E+00 -3.274470000000E-03\nS  1 1.00\n  9.056610000000E-01 1.000000000000E+00\nS  1 1.00\n  2.556110000000E-01 1.000000000000E+00\nS  1 1.00\n  8.450000000000E-02 1.000000000000E+00\nP  3 1.00\n  4.211750000000E+01 3.651140000000E-02\n  9.628370000000E+00 2.371530000000E-01\n  2.853320000000E+00 8.197020000000E-01\nP  1 1.00\n  9.056610000000E-01 1.000000000000E+00\nP  1 1.00\n  2.556110000000E-01 1.000000000000E+00\nP  1 1.00\n  8.450000000000E-02 1.000000000000E+00\nD  1 1.00\n  2.584000000000E+00 1.000000000000E+00\nD  1 1.00\n  6.460000000000E-01 1.000000000000E+00\nF  1 1.00\n  1.400000000000E+00 1.000000000000E+00\n****\n""",
            'F'  : """F  0 \nS  6 1.00\n  1.142710000000E+04 1.800930000000E-03\n  1.722350000000E+03 1.374190000000E-02\n  3.957460000000E+02 6.813340000000E-02\n  1.151390000000E+02 2.333250000000E-01\n  3.360260000000E+01 5.890860000000E-01\n  4.919010000000E+00 2.995050000000E-01\nS  3 1.00\n  5.544410000000E+01 1.145360000000E-01\n  1.263230000000E+01 9.205120000000E-01\n  3.717560000000E+00 -3.378040000000E-03\nS  1 1.00\n  1.165450000000E+00 1.000000000000E+00\nS  1 1.00\n  3.218920000000E-01 1.000000000000E+00\nS  1 1.00\n  1.076000000000E-01 1.000000000000E+00\nP  3 1.00\n  5.544410000000E+01 3.546090000000E-02\n  1.263230000000E+01 2.374510000000E-01\n  3.717560000000E+00 8.204580000000E-01\nP  1 1.00\n  1.165450000000E+00 1.000000000000E+00\nP  1 1.00\n  3.218920000000E-01 1.000000000000E+00\nP  1 1.00\n  1.076000000000E-01 1.000000000000E+00\nD  1 1.00\n  3.500000000000E+00 1.000000000000E+00\nD  1 1.00\n  8.750000000000E-01 1.000000000000E+00\nF  1 1.00\n  1.850000000000E+00 1.000000000000E+00\n****\n""",
            'Ne' : """Ne 0 \nS  6 1.00\n  1.399570000000E+04 1.832760000000E-03\n  2.117100000000E+03 1.388270000000E-02\n  4.904250000000E+02 6.806870000000E-02\n  1.438330000000E+02 2.313280000000E-01\n  4.192650000000E+01 5.858900000000E-01\n  6.156840000000E+00 3.058830000000E-01\nS  3 1.00\n  6.912110000000E+01 1.191490000000E-01\n  1.583500000000E+01 9.173750000000E-01\n  4.673260000000E+00 -4.058390000000E-03\nS  1 1.00\n  1.457560000000E+00 1.000000000000E+00\nS  1 1.00\n  3.970570000000E-01 1.000000000000E+00\nS  1 1.00\n  1.300000000000E-01 1.000000000000E+00\nP  3 1.00\n  6.912110000000E+01 3.565740000000E-02\n  1.583500000000E+01 2.394770000000E-01\n  4.673260000000E+00 8.184610000000E-01\nP  1 1.00\n  1.457560000000E+00 1.000000000000E+00\nP  1 1.00\n  3.970570000000E-01 1.000000000000E+00\nP  1 1.00\n  1.300000000000E-01 1.000000000000E+00\nD  1 1.00\n  4.608000000000E+00 1.000000000000E+00\nD  1 1.00\n  1.152000000000E+00 1.000000000000E+00\nF  1 1.00\n  2.500000000000E+00 1.000000000000E+00\n****\n""",
            'Na' : """Na 0 \nS  6 1.00\n  3.616640000000E+04 1.032000000000E-03\n  5.372580000000E+03 8.071000000000E-03\n  1.213210000000E+03 4.212900000000E-02\n  3.396230000000E+02 1.697890000000E-01\n  1.095530000000E+02 5.146210000000E-01\n  3.877730000000E+01 3.798170000000E-01\nS  3 1.00\n  3.877730000000E+01 3.747620000000E-01\n  1.457590000000E+01 5.757690000000E-01\n  5.269930000000E+00 1.129330000000E-01\nS  1 1.00\n  1.827770000000E+00 1.000000000000E+00\nS  1 1.00\n  6.199480000000E-01 1.000000000000E+00\nS  1 1.00\n  5.724000000000E-02 1.000000000000E+00\nS  1 1.00\n  2.404800000000E-02 1.000000000000E+00\nS  1 1.00\n  7.600000000000E-03 1.000000000000E+00\nP  4 1.00\n  1.446450000000E+02 1.148500000000E-02\n  3.390740000000E+01 8.238300000000E-02\n  1.062850000000E+01 3.196580000000E-01\n  3.823890000000E+00 7.012950000000E-01\nP  2 1.00\n  1.444290000000E+00 6.385060000000E-01\n  5.526210000000E-01 4.253650000000E-01\nP  1 1.00\n  1.887200000000E-01 1.000000000000E+00\nP  1 1.00\n  4.650100000000E-02 1.000000000000E+00\nP  1 1.00\n  1.628500000000E-02 1.000000000000E+00\nP  1 1.00\n  7.600000000000E-03 1.000000000000E+00\nD  1 1.00\n  7.000000000000E-01 1.000000000000E+00\nD  1 1.00\n  1.750000000000E-01 1.000000000000E+00\nD  1 1.00\n  4.375000000000E-02 1.000000000000E+00\nF  1 1.00\n  3.000000000000E-01 1.000000000000E+00\nF  1 1.00\n  7.500000000000E-02 1.000000000000E+00\n****\n""",
            'Mg' : """Mg 0 \nS  6 1.00\n  4.386650000000E+04 9.180000000000E-04\n  6.605370000000E+03 7.047000000000E-03\n  1.513260000000E+03 3.594100000000E-02\n  4.323170000000E+02 1.414610000000E-01\n  1.421490000000E+02 4.267640000000E-01\n  5.139830000000E+01 4.979750000000E-01\nS  3 1.00\n  5.139830000000E+01 2.513550000000E-01\n  1.991960000000E+01 6.186710000000E-01\n  8.024740000000E+00 1.884170000000E-01\nS  1 1.00\n  2.508170000000E+00 1.000000000000E+00\nS  1 1.00\n  8.715310000000E-01 1.000000000000E+00\nS  1 1.00\n  1.081880000000E-01 1.000000000000E+00\nS  1 1.00\n  4.013000000000E-02 1.000000000000E+00\nS  1 1.00\n  1.460000000000E-02 1.000000000000E+00\nP  4 1.00\n  1.938540000000E+02 1.018800000000E-02\n  4.544200000000E+01 7.536000000000E-02\n  1.418640000000E+01 3.074190000000E-01\n  5.057510000000E+00 7.175750000000E-01\nP  2 1.00\n  1.888610000000E+00 6.673390000000E-01\n  7.226520000000E-01 3.946490000000E-01\nP  1 1.00\n  2.364170000000E-01 1.000000000000E+00\nP  1 1.00\n  9.335800000000E-02 1.000000000000E+00\nP  1 1.00\n  3.480900000000E-02 1.000000000000E+00\nP  1 1.00\n  1.460000000000E-02 1.000000000000E+00\nD  1 1.00\n  7.000000000000E-01 1.000000000000E+00\nD  1 1.00\n  1.750000000000E-01 1.000000000000E+00\nD  1 1.00\n  4.375000000000E-02 1.000000000000E+00\nF  1 1.00\n  4.000000000000E-01 1.000000000000E+00\nF  1 1.00\n  1.000000000000E-01 1.000000000000E+00\n****\n""",
            'Al' : """Al 0 \nS  6 1.00\n  5.486648900000E+04 8.390000000000E-04\n  8.211766500000E+03 6.527000000000E-03\n  1.866176100000E+03 3.366600000000E-02\n  5.311293400000E+02 1.329020000000E-01\n  1.751179700000E+02 4.012660000000E-01\n  6.400550000000E+01 5.313380000000E-01\nS  3 1.00\n  6.400550000000E+01 2.023050000000E-01\n  2.529250700000E+01 6.247900000000E-01\n  1.053491000000E+01 2.274390000000E-01\nS  1 1.00\n  3.206711000000E+00 1.000000000000E+00\nS  1 1.00\n  1.152555000000E+00 1.000000000000E+00\nS  1 1.00\n  1.766780000000E-01 1.000000000000E+00\nS  1 1.00\n  6.523700000000E-02 1.000000000000E+00\nS  1 1.00\n  3.180000000000E-02 1.000000000000E+00\nP  4 1.00\n  2.592836200000E+02 9.448000000000E-03\n  6.107687000000E+01 7.097400000000E-02\n  1.930323700000E+01 2.956360000000E-01\n  7.010882000000E+00 7.282190000000E-01\nP  2 1.00\n  2.673865000000E+00 6.444670000000E-01\n  1.036596000000E+00 4.174130000000E-01\nP  1 1.00\n  3.168190000000E-01 1.000000000000E+00\nP  1 1.00\n  1.142570000000E-01 1.000000000000E+00\nP  1 1.00\n  4.139700000000E-02 1.000000000000E+00\nP  1 1.00\n  3.180000000000E-02 1.000000000000E+00\nD  1 1.00\n  1.300000000000E+00 1.000000000000E+00\nD  1 1.00\n  3.250000000000E-01 1.000000000000E+00\nD  1 1.00\n  8.125000000000E-02 1.000000000000E+00\nF  1 1.00\n  5.000000000000E-01 1.000000000000E+00\nF  1 1.00\n  1.250000000000E-01 1.000000000000E+00\n****\n""",
            'Si' : """Si 0 \nS  6 1.00\n  6.937923000000E+04 7.570000000000E-04\n  1.035494000000E+04 5.932000000000E-03\n  2.333879600000E+03 3.108800000000E-02\n  6.571429500000E+02 1.249670000000E-01\n  2.143011300000E+02 3.868970000000E-01\n  7.762916800000E+01 5.548880000000E-01\nS  3 1.00\n  7.762916800000E+01 1.778810000000E-01\n  3.063080700000E+01 6.277650000000E-01\n  1.280129500000E+01 2.476230000000E-01\nS  1 1.00\n  3.926866000000E+00 1.000000000000E+00\nS  1 1.00\n  1.452343000000E+00 1.000000000000E+00\nS  1 1.00\n  2.562340000000E-01 1.000000000000E+00\nS  1 1.00\n  9.427900000000E-02 1.000000000000E+00\nS  1 1.00\n  3.310000000000E-02 1.000000000000E+00\nP  4 1.00\n  3.354831900000E+02 8.866000000000E-03\n  7.890036600000E+01 6.829900000000E-02\n  2.498815000000E+01 2.909580000000E-01\n  9.219711000000E+00 7.321170000000E-01\nP  2 1.00\n  3.621140000000E+00 6.198790000000E-01\n  1.451310000000E+00 4.391480000000E-01\nP  1 1.00\n  5.049770000000E-01 1.000000000000E+00\nP  1 1.00\n  1.863170000000E-01 1.000000000000E+00\nP  1 1.00\n  6.543200000000E-02 1.000000000000E+00\nP  1 1.00\n  3.310000000000E-02 1.000000000000E+00\nD  1 1.00\n  1.800000000000E+00 1.000000000000E+00\nD  1 1.00\n  4.500000000000E-01 1.000000000000E+00\nD  1 1.00\n  1.125000000000E-01 1.000000000000E+00\nF  1 1.00\n  6.400000000000E-01 1.000000000000E+00\nF  1 1.00\n  1.600000000000E-01 1.000000000000E+00\n****\n""",
            'P'  : """P  0 \nS  6 1.00\n  7.749240000000E+04 7.869212000000E-04\n  1.160580000000E+04 6.108245000000E-03\n  2.645960000000E+03 3.139689000000E-02\n  7.549700000000E+02 1.242379000000E-01\n  2.487500000000E+02 3.811538000000E-01\n  9.115000000000E+01 5.595372000000E-01\nS  3 1.00\n  9.115000000000E+01 1.641617000000E-01\n  3.622000000000E+01 6.259097000000E-01\n  1.521000000000E+01 2.620744000000E-01\nS  1 1.00\n  4.710000000000E+00 1.000000000000E+00\nS  1 1.00\n  1.780000000000E+00 1.000000000000E+00\nS  1 1.00\n  3.400000000000E-01 1.000000000000E+00\nS  1 1.00\n  1.200000000000E-01 1.000000000000E+00\nS  1 1.00\n  3.480000000000E-02 1.000000000000E+00\nP  4 1.00\n  3.848400000000E+02 8.967875000000E-03\n  9.055000000000E+01 6.904902000000E-02\n  2.880000000000E+01 2.928770000000E-01\n  1.068000000000E+01 7.292494000000E-01\nP  2 1.00\n  4.250000000000E+00 6.325822000000E-01\n  1.740000000000E+00 4.232996000000E-01\nP  1 1.00\n  5.900000000000E-01 1.000000000000E+00\nP  1 1.00\n  2.200000000000E-01 1.000000000000E+00\nP  1 1.00\n  8.000000000000E-02 1.000000000000E+00\nP  1 1.00\n  3.480000000000E-02 1.000000000000E+00\nD  1 1.00\n  2.200000000000E+00 1.000000000000E+00\nD  1 1.00\n  5.500000000000E-01 1.000000000000E+00\nD  1 1.00\n  1.375000000000E-01 1.000000000000E+00\nF  1 1.00\n  9.000000000000E-01 1.000000000000E+00\nF  1 1.00\n  2.250000000000E-01 1.000000000000E+00\n****\n""",
            'S'  : """S  0 \nS  6 1.00\n  9.341340000000E+04 7.420791000000E-04\n  1.396170000000E+04 5.787658000000E-03\n  3.169910000000E+03 2.994067000000E-02\n  9.024500000000E+02 1.189282000000E-01\n  2.971500000000E+02 3.681822000000E-01\n  1.087000000000E+02 5.776336000000E-01\nS  3 1.00\n  1.087000000000E+02 1.427905000000E-01\n  4.315000000000E+01 6.246934000000E-01\n  1.810000000000E+01 2.834835000000E-01\nS  1 1.00\n  5.570000000000E+00 1.000000000000E+00\nS  1 1.00\n  2.140000000000E+00 1.000000000000E+00\nS  1 1.00\n  4.300000000000E-01 1.000000000000E+00\nS  1 1.00\n  1.500000000000E-01 1.000000000000E+00\nS  1 1.00\n  4.050000000000E-02 1.000000000000E+00\nP  4 1.00\n  4.950400000000E+02 8.196253000000E-03\n  1.172200000000E+02 6.364204000000E-02\n  3.750000000000E+01 2.788060000000E-01\n  1.391000000000E+01 7.447404000000E-01\nP  2 1.00\n  5.500000000000E+00 6.168248000000E-01\n  2.240000000000E+00 4.402946000000E-01\nP  1 1.00\n  7.700000000000E-01 1.000000000000E+00\nP  1 1.00\n  2.900000000000E-01 1.000000000000E+00\nP  1 1.00\n  1.000000000000E-01 1.000000000000E+00\nP  1 1.00\n  4.050000000000E-02 1.000000000000E+00\nD  1 1.00\n  2.600000000000E+00 1.000000000000E+00\nD  1 1.00\n  6.500000000000E-01 1.000000000000E+00\nD  1 1.00\n  1.625000000000E-01 1.000000000000E+00\nF  1 1.00\n  1.100000000000E+00 1.000000000000E+00\nF  1 1.00\n  2.750000000000E-01 1.000000000000E+00\n****\n""",
            'Cl' : """Cl 0 \nS  6 1.00\n  1.058190000000E+05 7.423627000000E-04\n  1.587200000000E+04 5.747318000000E-03\n  3.619650000000E+03 2.964876000000E-02\n  1.030800000000E+03 1.178998000000E-01\n  3.399000000000E+02 3.648532000000E-01\n  1.245300000000E+02 5.816968000000E-01\nS  3 1.00\n  1.245300000000E+02 1.370443000000E-01\n  4.951000000000E+01 6.231380000000E-01\n  2.080000000000E+01 2.903279000000E-01\nS  1 1.00\n  6.460000000000E+00 1.000000000000E+00\nS  1 1.00\n  2.520000000000E+00 1.000000000000E+00\nS  1 1.00\n  5.300000000000E-01 1.000000000000E+00\nS  1 1.00\n  1.900000000000E-01 1.000000000000E+00\nS  1 1.00\n  4.830000000000E-02 1.000000000000E+00\nP  4 1.00\n  5.897800000000E+02 7.873332000000E-03\n  1.398500000000E+02 6.155460000000E-02\n  4.479000000000E+01 2.742514000000E-01\n  1.661000000000E+01 7.498994000000E-01\nP  2 1.00\n  6.590000000000E+00 6.147640000000E-01\n  2.710000000000E+00 4.413416000000E-01\nP  1 1.00\n  9.500000000000E-01 1.000000000000E+00\nP  1 1.00\n  3.500000000000E-01 1.000000000000E+00\nP  1 1.00\n  1.200000000000E-01 1.000000000000E+00\nP  1 1.00\n  4.830000000000E-02 1.000000000000E+00\nD  1 1.00\n  3.000000000000E+00 1.000000000000E+00\nD  1 1.00\n  7.500000000000E-01 1.000000000000E+00\nD  1 1.00\n  1.875000000000E-01 1.000000000000E+00\nF  1 1.00\n  1.400000000000E+00 1.000000000000E+00\nF  1 1.00\n  3.500000000000E-01 1.000000000000E+00\n****\n""",
            'Ar' : """Ar 0 \nS  6 1.00\n  1.180220000000E+05 7.461902000000E-04\n  1.768350000000E+04 5.786362000000E-03\n  4.027770000000E+03 2.990098000000E-02\n  1.145400000000E+03 1.191287000000E-01\n  3.771600000000E+02 3.687839000000E-01\n  1.381600000000E+02 5.767726000000E-01\nS  3 1.00\n  1.381600000000E+02 1.435931000000E-01\n  5.498000000000E+01 6.231142000000E-01\n  2.317000000000E+01 2.840810000000E-01\nS  1 1.00\n  7.370000000000E+00 1.000000000000E+00\nS  1 1.00\n  2.920000000000E+00 1.000000000000E+00\nS  1 1.00\n  6.500000000000E-01 1.000000000000E+00\nS  1 1.00\n  2.300000000000E-01 1.000000000000E+00\nS  1 1.00\n  6.000000000000E-02 1.000000000000E+00\nP  4 1.00\n  6.630600000000E+02 7.820021000000E-03\n  1.570900000000E+02 6.148333000000E-02\n  5.023000000000E+01 2.754731000000E-01\n  1.863000000000E+01 7.488402000000E-01\nP  2 1.00\n  7.440000000000E+00 -6.282210000000E-01\n  3.090000000000E+00 -4.260202000000E-01\nP  1 1.00\n  1.100000000000E+00 1.000000000000E+00\nP  1 1.00\n  4.100000000000E-01 1.000000000000E+00\nP  1 1.00\n  1.400000000000E-01 1.000000000000E+00\nP  1 1.00\n  6.000000000000E-02 1.000000000000E+00\nD  1 1.00\n  3.400000000000E+00 1.000000000000E+00\nD  1 1.00\n  8.500000000000E-01 1.000000000000E+00\nD  1 1.00\n  2.125000000000E-01 1.000000000000E+00\nF  1 1.00\n  1.700000000000E+00 1.000000000000E+00\nF  1 1.00\n  4.250000000000E-01 1.000000000000E+00\n****""",
            }


class GaussianTSB3LYP(GaussianTS):

    #: Keywords that will be added at the top of the qm input file
    keywords = [
               "# b3lyp/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest) int=ultrafine nosymm",
               "# b3lyp/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest,cartesian) int=ultrafine geom=allcheck guess=check nosymm",
               "# b3lyp/6-31+g(d,p) opt=(ts,calcall,noeigentest) nosymm",
               "# b3lyp/6-31+g(d,p) opt=(ts,calcall,noeigentest,cartesian) nosymm geom=allcheck guess=check nosymm",
               "# b3lyp/6-31+g(d,p) irc=(calcall,report=read) geom=allcheck guess=check nosymm",
               ]
    """
    This needs some work, to determine options that are best used. Commented out the
    methods for now.
    """
    
    def getQMMolecule(self, molecule):
        """
        The TST calculation must use the same electronic structure and basis set for the
        reactant species as the transition state. This method will ensure this by creating
        and returning the corresponding QMMolecule from the child class GaussianMolB3LYP.
        """
        
        qmMolecule = GaussianMolB3LYP(molecule, self.settings)
        return qmMolecule
