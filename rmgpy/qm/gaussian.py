import os

import openbabel
import external.cclib as cclib
import logging
from subprocess import Popen, PIPE

from qmdata import CCLibData
from molecule import QMMolecule
from reaction import QMReaction

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
    
    
    
class GaussianMol(QMMolecule, Gaussian):
    """
    A base Class for calculations of molecules using Gaussian. 
    
    Inherits from both :class:`QMMolecule` and :class:`Gaussian`.
    """
    
    #: List of phrases that indicate failure
    #: NONE of these must be present in a succesful job.
    failureKeys = [
                   'ERROR TERMINATION',
                   'IMAGINARY FREQUENCIES'
                   ]
    
    def writeInputFile(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
    
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "gjf")
        mol = openbabel.OBMol()
    
        obConversion.ReadFile(mol, self.getMolFilePathForCalculation(attempt) )
    
        mol.SetTitle(self.geometry.uniqueIDlong)
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
        input_string = obConversion.WriteString(mol)
        top_keys = self.inputFileKeywords(attempt)
        with open(self.inputFilePath, 'w') as gaussianFile:
            gaussianFile.write(top_keys)
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
        Calculate the QM data and return a QMData object.
        """
        self.createGeometry()
        if self.verifyOutputFile():
            logging.info("Found a successful output file already; using that.")
        else:
            success = False
            for attempt in range(1, self.maxAttempts+1):
                self.writeInputFile(attempt)
                success = self.run()
                if success:
                    logging.info('Attempt {0} of {1} on species {2} succeeded.'.format(attempt, self.maxAttempts, self.molecule.toAugmentedInChI()))
                    break
            else:
                raise Exception('QM thermo calculation failed for {0}.'.format(self.molecule.toAugmentedInChI()))
        result = self.parse() # parsed in cclib
        return result
    
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
                    if logFileInChI == self.geometry.uniqueIDlong:
                        InChIMatch = True
                    elif self.geometry.uniqueIDlong.startswith(logFileInChI):
                        logging.info("InChI too long to check, but beginning matches so assuming OK.")
                        InChIMatch = True
                    else:
                        logging.info("InChI in log file didn't match that in geometry.")
                        logging.info(self.geometry.uniqueIDlong)
                        logging.info(logFileInChI)
        
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

class GaussianMolPM3(GaussianMol):

    #: Keywords that will be added at the top of the qm input files
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

##########################################################################################

class GaussianTS(QMReaction, Gaussian):
    """
    A base Class for calculations of transition states using Gaussian. 

    Inherits from both :class:`QMReaction` and :class:`Gaussian`.
    """
    
    
    #: List of phrases that indicate failure
    #: NONE of these must be present in a succesful job.
    failureKeys = [
                   'ERROR TERMINATION',
                   ]
    
    def writeInputFile(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        chk_file = '%chk=' + self.uniqueID
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "gjf")
        mol = openbabel.OBMol()
        
        obConversion.ReadFile(mol, self.geometry.getRefinedMolFilePath() )
        
        mol.SetTitle(self.uniqueID)
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
        input_string = obConversion.WriteString(mol)
        numProc = '%nprocshared=' + '4' # could be something that could be set in the qmSettings
        top_keys = self.keywords[0]
        title = ' ' + self.uniqueID
        with open(self.inputFilePath, 'w') as gaussianFile:
            gaussianFile.write(numProc)
            gaussianFile.write('\n')
            gaussianFile.write(chk_file)
            gaussianFile.write('\n')
            gaussianFile.write(top_keys)
            gaussianFile.write(input_string)
    
    def writeIRCFile(self):
        """
        Using the :class:`Geometry` object, write the input file for the 
        IRC calculation on the transition state. The geometry is taken 
        from the checkpoint file created during the geometry search.
        """
        chk_file = '%chk=' + self.uniqueID
        numProc = '%nprocshared=' + '4' # could be something that could be set in the qmSettings
        top_keys = self.keywords[2]
        chrgMult = '1 2'
        with open(self.inputFilePath, 'w') as gaussianFile:
            gaussianFile.write(numProc)
            gaussianFile.write('\n')
            gaussianFile.write(chk_file)
            gaussianFile.write('\n')
            gaussianFile.write(top_keys)
            gaussianFile.write('\n\n')
            gaussianFile.write(chrgMult)
            gaussianFile.write('\n\n')
            
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
                raise Exception('QM thermo calculation failed for {0}.'.format(self.molecule.toAugmentedInChI()))
        result = self.parse() # parsed in cclib
        return result
    
    def verifyOutputFile(self):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful GAUSSIAN simulation already exists for the molecule with the 
        given (augmented) InChI Key.
        
        The definition of finding a successful simulation is based on these criteria:
        1) finding an output file with the file name equal to the the reaction unique ID
        2) NOT finding any of the keywords that are denote a calculation failure
        3) finding all the keywords that denote a calculation success.
        
        If any of the above criteria is not matched, False will be returned and the procedures to start a new calculation 
        will be initiated.
        """
        if not os.path.exists(self.outputFilePath):
            logging.info("Output file {0} does not exist.".format(self.outputFilePath))
            return False
        
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
        
        # Check that ALL 'success' keywords were found in the file.
        if not all( successKeysFound.values() ):
            logging.error('Not all of the required keywords for sucess were found in the output file!')
            return False
        else:
            return True
    
    def verifyTSGeometry(self):
        """
        Check's that the resulting geometries of the path analysis match the reaction.
        """
        
        """
        Compares IRC geometries to input geometries.
        """
        # Search IRC output for steps on each side of the path
        readFile = file(ircOutput)
        pth1 = list()
        steps = list()
        for line in readFile.readlines():
            if line.startswith(' Point Number:'):
                if int(line.split()[2]) > 0:
                    if int(line.split()[-1]) == 1:
                        ptNum = int(line.split()[2])
                        pth1.append(ptNum)
                    else:
                        pass
            elif line.startswith('  # OF STEPS ='):
                numStp = int(line.split()[-1])
                steps.append(numStp)
        
        # This indexes the coordinate to be used from the parsing
        if steps == []:
            notes = ' IRC failed '
            return 0, notes
        else:
            pth1End = sum(steps[:pth1[-1]])		
            # Compare the reactants and products
            ircParse = external.cclib.parser.Gaussian(ircOutput)
            ircParse = ircParse.parse()
        
            atomnos = ircParse.atomnos
            atomcoords = ircParse.atomcoords
        
            # Convert the IRC geometries into RMG molecules
            # We don't know which is reactant or product, so take the two at the end of the
            # paths and compare to the reactants and products
        
            mol1 = fromXYZ(atomcoords[pth1End], atomnos)
            mol2 = fromXYZ(atomcoords[-1], atomnos)
            
            # Had trouble with isIsomorphic, but resetting the connectivity seems to fix it (WHY??).
            reactant.resetConnectivityValues()
            product.resetConnectivityValues()
            mol1.resetConnectivityValues()
            mol2.resetConnectivityValues()
            
            if reactant.isIsomorphic(mol1) and product.isIsomorphic(mol2):
                    notes = 'Verified TS'
                    return 1, notes
            elif reactant.isIsomorphic(mol2) and product.isIsomorphic(mol1):
                    notes = 'Verified TS'
                    return 1, notes
            else:
                notes = 'Saddle found, but wrong one'
                return 0, notes
        
class GaussianTSM062X(GaussianTS):

    #: Keywords that will be added at the top of the qm input file
    keywords = [
               "# m062x/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest)  int=ultrafine nosymm",
               "# m062x/6-31+g(d,p) opt=(ts,calcall,noeigentest) nosymm",
               "# m062x/6-31+g(d,p) irc=(calcall,report=read) geom=allcheck guess=check nosymm",
               ]
    """
    This needs some work, to determine options that are best used. Commented out the
    methods for now.
    """
    
    # @property
    # def scriptAttempts(self):
    #     "The number of attempts with different script keywords"
    #     return len(self.keywords)
    # 
    # @property
    # def maxAttempts(self):
    #     "The total number of attempts to try"
    #     return 2 * len(self.keywords)
    # 
    # def inputFileKeywords(self, attempt):
    #     """
    #     Return the top keywords for attempt number `attempt`.
    # 
    #     NB. `attempt`s begin at 1, not 0.
    #     """
    #     assert attempt <= self.maxAttempts
    #     if attempt > self.scriptAttempts:
    #         attempt -= self.scriptAttempts
    #     return self.keywords[attempt-1]