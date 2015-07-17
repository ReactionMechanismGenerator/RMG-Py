import os

import external.cclib as cclib
import logging
from subprocess import Popen
import re

from rmgpy.molecule import Molecule
from qmdata import CCLibData
from molecule import QMMolecule

class Gaussian:
    """
    A base class for all QM calculations that use Gaussian.
    
    Classes such as :class:`GaussianMol` will inherit from this class.
    """
    
    inputFileExtension = '.gjf'
    outputFileExtension = '.log'
    
    gaussEnv = os.getenv('GAUSS_EXEDIR') or os.getenv('g09root') or os.getenv('g03root') or ""
    
    # GAUSS_EXEDIR may be a list like "path1:path2:path3"
    for possibleDir in gaussEnv.split(':'):
        if os.path.exists(os.path.join(possibleDir , 'g09')):
            executablePath = os.path.join(possibleDir , 'g09')
            break
        elif os.path.exists(os.path.join(possibleDir , 'g03')):
            executablePath = os.path.join(possibleDir , 'g03')
            break
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
        5) checking that the optimized geometry, when connected by single bonds, is isomorphic with self.molecule (converted to single bonds)

        If any of the above criteria is not matched, False will be returned.
        If all are satisfied, it will return True.
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
                    if self.geometry.uniqueIDlong in logFileInChI:
                        InChIMatch = True
                    elif self.geometry.uniqueIDlong.startswith(logFileInChI):
                        logging.info("InChI too long to check, but beginning matches so assuming OK.")
                        InChIMatch = True
                    else:
                        logging.warning("InChI in log file ({0}) didn't match that in geometry ({1}).".format(logFileInChI, self.geometry.uniqueIDlong))                    
                        if self.geometry.uniqueIDlong.startswith(logFileInChI):
                            logging.warning("but the beginning matches so it's probably just a truncation problem.")
                            InChIMatch = True
        # Check that ALL 'success' keywords were found in the file.
        if not all( successKeysFound.values() ):
            logging.error('Not all of the required keywords for success were found in the output file!')
            return False
        
        if not InChIFound:
            logging.error("No InChI was found in the Gaussian output file {0}".format(self.outputFilePath))
            return False
        
        if not InChIMatch:
            #InChIs do not match (most likely due to limited name length mirrored in log file (240 characters), but possibly due to a collision)
            return self.checkForInChiKeyCollision(logFileInChI) # Not yet implemented!

        # Compare the optimized geometry to the original molecule
        qmData = self.parse()
        cclibMol = Molecule()
        cclibMol.fromXYZ(qmData.atomicNumbers, qmData.atomCoords.value)
        testMol = self.molecule.toSingleBonds()
        if not cclibMol.isIsomorphic(testMol):
            logging.info("Incorrect connectivity for optimized geometry in file {0}".format(self.outputFilePath))
            return False

        logging.info("Successful {1} quantum result in {0}".format(self.outputFilePath, self.__class__.__name__))
        return True
        
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
    
    def inputFileKeywords(self, attempt):
        """
        Return the top keywords for attempt number `attempt`.
    
        NB. `attempt` begins at 1, not 0.
        """
        assert attempt <= self.maxAttempts
        if attempt > self.scriptAttempts:
            attempt -= self.scriptAttempts
        return self.keywords[attempt-1]
    
    def writeInputFile(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attempt`.
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
                raise NotImplementedError("Not sure what should be here, if anything.")
                #gaussianFile.write(polar_keys)
    
    def generateQMData(self):
        """
        Calculate the QM data and return a QMData object.
        """
        for atom in self.molecule.vertices:
            if atom.atomType.label in ('N5s', 'N5d', 'N5dd', 'N5t', 'N5b'):
                return None
                
        if self.verifyOutputFile():
            logging.info("Found a successful output file already; using that.")
            source = "QM {0} result file found from previous run.".format(self.__class__.__name__)
        else:
            self.createGeometry()
            success = False
            for attempt in range(1, self.maxAttempts+1):
                self.writeInputFile(attempt)
                logging.info('Trying {3} attempt {0} of {1} on molecule {2}.'.format(attempt, self.maxAttempts, self.molecule.toSMILES(), self.__class__.__name__))
                success = self.run()
                if success:
                    logging.info('Attempt {0} of {1} on species {2} succeeded.'.format(attempt, self.maxAttempts, self.molecule.toAugmentedInChI()))
                    source = "QM {0} calculation attempt {1}".format(self.__class__.__name__, attempt )
                    break
            else:
                logging.error('QM thermo calculation failed for {0}.'.format(self.molecule.toAugmentedInChI()))
                return None
        result = self.parse() # parsed in cclib
        result.source = source
        return result # a CCLibData object
        
    def getParser(self, outputFile):
        """
        Returns the appropriate cclib parser.
        """
        return cclib.parser.Gaussian(outputFile)

class GaussianMolPM3(GaussianMol):
    """
    Gaussian PM3 calculations for molecules
    
    This is a class of its own in case you wish to do anything differently,
    but for now it's only the 'pm3' in the keywords that differs.
    """
    #: Keywords that will be added at the top of the qm input file
    keywords = [
                # The combinations of keywords were derived by Greg Magoon for pm3 in Gaussian. His comments are attached to each combination.
               "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3)", # added IOP option to avoid aborting when symmetry changes; 3 is supposed to be default according to documentation, but it seems that 0 (the default) is the only option that doesn't work from 0-4; also, it is interesting to note that all 4 options seem to work for test case with z-matrix input rather than xyz coords; cf. http://www.ccl.net/cgi-bin/ccl/message-new?2006+10+17+005 for original idea for solution
               "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)", # use different SCF method; this addresses at least one case of failure for a C4H7J species
               "# pm3 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm" , # try multiple different options (no gdiis, use calcfc, nosymm); 7/21/09: added maxcyc option to fix case of MPTBUKVAJYJXDE-UHFFFAOYAPmult3 (InChI=1/C4H10O5Si/c1-3-7-9-10(5,6)8-4-2/h4-5H,3H2,1-2H3/mult3) (file manually copied to speed things along)
               "# pm3 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm", # numerical frequency keyword version of keyword #3; used to address GYFVJYRUZAKGFA-UHFFFAOYALmult3 (InChI=1/C6H14O6Si/c1-3-10-13(8,11-4-2)12-6-5-9-7/h6-7H,3-5H2,1-2H3/mult3) case; (none of the existing Gaussian or MOPAC combinations worked with it)
               "# pm3 opt=(verytight,gdiis,small) freq IOP(2/16=3)", #  somehow, this worked for problematic case of ZGAWAHRALACNPM-UHFFFAOYAF (InChI=1/C8H17O5Si/c1-3-11-14(10,12-4-2)13-8-5-7(9)6-8/h7-9H,3-6H2,1-2H3); (was otherwise giving l402 errors); even though I had a keyword that worked for this case, I manually copied the fixed log file to QMfiles folder to speed things along; note that there are a couple of very low frequencies (~5-6 cm^-1 for this case)
               "# pm3 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)", # used for troublesome C5H7J2 case (similar error to C5H7J below); calcfc is not necessary for this particular species, but it speeds convergence and probably makes it more robust for other species
               "# pm3 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)", # use numerical frequencies; this takes a relatively long time, so should only be used as one of the last resorts; this seemed to address at least one case of failure for a C6H10JJ species; 7/15/09: maxcyc=200 added to address GVCMURUDAUQXEY-UHFFFAOYAVmult3 (InChI=1/C3H4O7Si/c1-2(9-6)10-11(7,8)3(4)5/h6-7H,1H2/mult3)...however, result was manually pasted in QMfiles folder to speed things along
               "# pm3 opt=tight freq IOP(2/16=3)", # this worked for problematic case of SZSSHFMXPBKYPR-UHFFFAOYAF (InChI=1/C7H15O5Si/c1-3-10-13(8,11-4-2)12-7-5-6-9-7/h7H,3-6H2,1-2H3) (otherwise, it had l402.exe errors); corrected log file was manually copied to QMfiles to speed things along; we could also add a freq=numerical version of this keyword combination for added robustness; UPDATE: see below
               "# pm3 opt=tight freq=numerical IOP(2/16=3)", # used for problematic case of CIKDVMUGTARZCK-UHFFFAOYAImult4 (InChI=1/C8H15O6Si/c1-4-12-15(10,13-5-2)14-7-6-11-8(7,3)9/h7H,3-6H2,1-2H3/mult4 (most other cases had l402.exe errors); corrected log file was manually copied to QMfiles to speed things along
               "# pm3 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)", # similar to existing #5, but uses tight rather than verytight; used for ADMPQLGIEMRGAT-UHFFFAOYAUmult3 (InChI=1/C6H14O5Si/c1-4-9-12(8,10-5-2)11-6(3)7/h6-7H,3-5H2,1-2H3/mult3)
               "# pm3 opt freq IOP(2/16=3)", # use default (not verytight) convergence criteria; use this as last resort
               "# pm3 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)", # to address problematic C10H14JJ case
               "# pm3 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm", # for very troublesome RRMZRNPRCUANER-UHFFFAOYAQ (InChI=1/C5H7/c1-3-5-4-2/h3H,1-2H3) case...there were troubles with negative frequencies, where I don't think they should have been; step size of numerical frequency was adjusted to give positive result; accuracy of result is questionable; it is possible that not all of these keywords are needed; note that for this and other nearly free rotor cases, I think heat capacity will be overestimated by R/2 (R vs. R/2) (but this is a separate issue)
               "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm", #  for troublesome QDERTVAGQZYPHT-UHFFFAOYAHmult3(InChI=1/C6H14O4Si/c1-4-8-11(7,9-5-2)10-6-3/h4H,5-6H2,1-3H3/mult3); key aspects appear to be tight (rather than verytight) convergence criteria, no calculation of frequencies during optimization, use of numerical frequencies, and probably also the use of opt=small
               "# pm3 opt=(verytight,gdiis,calcall) IOP(2/16=3)", # used for troublesome C5H7J case; note that before fixing, I got errors like the following: "Incomplete coordinate system.  Try restarting with Geom=Check Guess=Read Opt=(ReadFC,NewRedundant) Incomplete coordinate system. Error termination via Lnk1e in l103.exe"; we could try to restart, but it is probably preferrable to have each keyword combination standalone; another keyword that may be helpful if additional problematic cases are encountered is opt=small; 6/9/09 note: originally, this had # pm3 opt=(verytight,gdiis,calcall) freq IOP(2/16=3)" (with freq keyword), but I discovered that in this case, there are two thermochemistry sections and cclib parses frequencies twice, giving twice the number of desired frequencies and hence produces incorrect thermo; this turned up on C5H6JJ isomer 
               "# pm3 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm", # worked for troublesome ketene case: CCGKOQOJPYTBIH-UHFFFAOYAO (InChI=1/C2H2O/c1-2-3/h1H2) (could just increase number of iterations for similar keyword combination above (#6 at the time of this writing), allowing symmetry, but nosymm seemed to reduce # of iterations; I think one of nosymm or higher number of iterations would allow the similar keyword combination to converge; both are included here for robustness)
               "# pm3 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm", # added for case of ZWMVZWMBTVHPBS-UHFFFAOYAEmult3 (InChI=1/C4H4O2/c1-3-5-6-4-2/h1-2H2/mult3)
               "# pm3 opt=(calcall,small,maxcyc=100) IOP(2/16=3)", # used to address troublesome FILUFGAZMJGNEN-UHFFFAOYAImult3 case (InChI=1/C5H6/c1-3-5-4-2/h3H,1H2,2H3/mult3)
               ]

class GaussianMolPM6(GaussianMol):
    """
    Gaussian PM6 calculations for molecules
    
    This is a class of its own in case you wish to do anything differently,
    but for now it's only the 'pm6' in the keywords that differs.
    """
    #: Keywords that will be added at the top of the qm input file
    keywords = [
                # The combinations of keywords were derived by Greg Magoon for pm3. For now, we assume similar ones will work for pm6:
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
               "# pm6 opt=(verytight,gdiis,calcall) IOP(2/16=3)",
               "# pm6 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm",
               "# pm6 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm",
               "# pm6 opt=(calcall,small,maxcyc=100) IOP(2/16=3)",
               ]
