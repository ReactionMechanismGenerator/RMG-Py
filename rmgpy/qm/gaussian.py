import os

import re
import external.cclib as cclib
import logging
import time
import math
import numpy
import itertools
from subprocess import Popen, PIPE

from rmgpy.molecule import Molecule, Atom, getElement
from rmgpy.species import Species
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
        
        with open(self.inputFilePath) as infile:
            print "Running GAUSSIAN input file {0!s}:".format(self.inputFilePath)
            for line in infile:
                print line.rstrip()

        # submits the input file to Gaussian
        process = Popen([self.executablePath, self.inputFilePath, self.outputFilePath])
        process.communicate()# necessary to wait for executable termination!
            
        with open(self.outputFilePath) as outfile:
            print "Gaussian output file {0!s}:".format(self.outputFilePath)
            for line in outfile:
                print line.rstrip()

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
                    if logFileInChI == self.geometry.uniqueIDlong:
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

        logging.info("Successful Gaussian quantum result found in {0}".format(self.outputFilePath))
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
    
    def writeInputFile(self, output, attempt=None, top_keys=None, numProcShared=None, memory=None, checkPoint=False, bottomKeys=None):
        """
        Takes the output from the createInputFile method and prints the
        file. Options provided allow the 
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        if not top_keys:
            top_keys = self.inputFileKeywords(attempt)
        output = [top_keys] + output
        
        if checkPoint:
            chk_file = '%chk=' + os.path.join(self.settings.fileStore, self.uniqueID)
            output = [chk_file] + output
        if memory:
            mem = '%mem={0}'.format(memory)
            output = [mem] + output
        if numProc:
            numProc = '%nprocshared={0}'.format(numProcShared)
            output = [numProc] + output
        if bottomKeys:
            output = output + [bottomKeys]
        
        input_string = '\n'.join(output)
        
        with open(self.inputFilePath, 'w') as gaussianFile:
            gaussianFile.write(input_string)
            gaussianFile.write('\n')                
    
class GaussianMol(QMMolecule, Gaussian):
    """
    A base Class for calculations of molecules using Gaussian. 
    
    Inherits from both :class:`QMMolecule` and :class:`Gaussian`.
    """
    
    def inputFileKeywords(self, attempt):
        """
        Return the top keywords for attempt number `attempt`.
    
        NB. `attempt`s begin at 1, not 0.
        """
        assert attempt <= self.maxAttempts
        if attempt > self.scriptAttempts:
            attempt -= self.scriptAttempts
        return self.keywords[attempt-1]
    
    def createInputFile(self, attempt):
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
        self.writeInputFile(output, attempt)
    
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
                self.createInputFile(attempt)
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
                   'Error termination',
                   'Error in internal coordinate system.',
                   ]
    
    def inputFileKeywords(self, attempt):
        """
        Return the top keywords for attempt number `attempt`.
    
        NB. `attempt`s begin at 1, not 0.
        """
        
        return self.keywords[attempt-1]
    
    def createInputFile(self, attempt, fromSddl=False, fromQST2=False, fromInt=False, fromNEB=False):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        
        if fromSddl:
            molfile = self.getFilePath('.arc')
            atomline = re.compile('\s*([A-Za-z]+)\s+([\- ][0-9.]+)\s+([\+ ][0-9.]+)\s+([\- ][0-9.]+)\s+([\+ ][0-9.]+)\s+([\- ][0-9.]+)')
            
            output = ['', self.geometry.uniqueID, '' ]
            output.append("{charge}   {mult}".format(charge=0, mult=(fromSddl + 1) ))
            
            atomCount = 0
            with open(molfile) as molinput:
                for line in molinput:
                    match = atomline.match(line)
                    if match:
                        output.append("{0:8s} {1}  {2}  {3}".format(match.group(1), match.group(2), match.group(4), match.group(6)))
                        atomCount += 1
        elif fromNEB:
            output = ['', self.geometry.uniqueID, '' ]
            output.append("{charge}   {mult}".format(charge=0, mult=self.geometry.multiplicity ))
            
            filePath = self.getFilePath('peak.xyz')
            assert os.path.exists(filePath)
            atomsymbols, atomcoords = self.geometry.parseXYZ(filePath)
            atomCount = 0
            for symbol, coordinates in zip(atomsymbols, atomcoords):
                output.append("{0:8s} {1:+.6f}  {2:+.6f}  {3:+.6f}".format(symbol, coordinates[0], coordinates[1], coordinates[2]).replace('+',' '))
                atomCount += 1
        elif fromQST2:
            # output = []
            # atomCount = len(self.geometry.molecule.atoms)
            # Really should be using above as we take advantage of gaussian's checkpoint files
            # However with the seg faults on the Discovery cluster, we use below as a temporary workaround.
            # Also see the keywords that we have changed to use the workaround.
            
            output = ['', self.geometry.uniqueID, '' ]
            output.append("{charge}   {mult}".format(charge=0, mult=(self.geometry.molecule.getRadicalCount() + 1) ))
            
            
            molfile = self.geometry.getRefinedMolFilePath()
            atomline = re.compile('\s*([\- ][0-9.]+\s+[\-0-9.]+\s+[\-0-9.]+)\s+([A-Za-z]+)')
            
            atomCount = 0
            atomnos = []
            with open(molfile) as molinput:
                for line in molinput:
                    match = atomline.match(line)
                    if match:
                        atomnos.append(match.group(2))
                        atomCount += 1
            
            assert atomCount == len(self.geometry.molecule.atoms)
            
            parser = cclib.parser.Gaussian(self.outputFilePath)
            parser.logger.setLevel(logging.ERROR)
            parser = parser.parse()
            
            lines = []
            for line in parser.atomcoords[-1]:
                lineList = ''
                for item in line:
                    if item > 0:
                        lineList = lineList + ' ' + str(format(item, '.6f')) + ' '
                    else:
                        lineList = lineList + str(format(item, '.6f')) + ' '
                lines.append(lineList)
                
            atomCount = 0
            for i, line in enumerate(lines):
                output.append("{0:8s} {1}".format(atomnos[i], line))
                atomCount += 1
        elif fromInt:
            output = ['', self.geometry.uniqueID, '' ]
            output.append("{charge}   {mult}".format(charge=0, mult=self.geometry.multiplicity ))
            
            atomsymbols, atomcoords = self.geometry.parseLOG(self.outputFilePath)
            atomCount = 0
            for symbol, coordinates in zip(atomsymbols, atomcoords):
                output.append("{0:8s} {1:+.6f}  {2:+.6f}  {3:+.6f}".format(symbol, coordinates[0], coordinates[1], coordinates[2]).replace('+',' '))
                atomCount += 1
        elif attempt > 2:
            # Until checkpointing is fixed, do the following
            output = ['', self.geometry.uniqueID, '' ]
            output.append("{charge}   {mult}".format(charge=0, mult=(self.geometry.molecule.getRadicalCount() + 1) ))
            
            parser = cclib.parser.Gaussian(self.outputFilePath)
            parser.logger.setLevel(logging.ERROR)
            parser = parser.parse()
            
            lines = []
            for line in parser.atomcoords[-1]:
                lineList = ''
                for item in line:
                    if item > 0:
                        lineList = lineList + ' ' + str(format(item, '.6f')) + ' '
                    else:
                        lineList = lineList + str(format(item, '.6f')) + ' '
                lines.append(lineList)
                
            atomCount = 0
            atomnos = parser.atomnos
            for i, line in enumerate(lines):
                output.append("{0:8s} {1}".format(getElement(int(atomnos[i])).symbol, line))
                atomCount += 1
        else:
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
        self.writeInputFile(output, attempt, numProcShared=20, memory='800MB', checkPoint=True)
    
    def createIRCFile(self):
        """
        Using the :class:`Geometry` object, write the input file for the 
        IRC calculation on the transition state. The geometry is taken 
        from the checkpoint file created during the geometry search.
        """
        
        top_keys = self.keywords[4]
        output = "{charge}   {mult}".format(charge=0, mult=(self.geometry.molecule.getRadicalCount() + 1) )
        
        # atomTypes = []
        # for atom in self.geometry.molecule.atoms:
        #     if not atom.element.symbol in atomTypes:
        #         atomTypes.append(atom.element.symbol)
        self.writeInputFile(output, top_keys=top_keys, numProcShared=20, memory='800MB', checkPoint=True)
        
    def setImages(self, pGeom):
        """
        Set and return the initial and final ase images for the NEB calculation
        """
        import ase
        from ase import io, Atoms
        
        # Give ase the atom positions for each side of the reaction path
        atomsymbols, atomcoords = self.geometry.parseLOG(self.outputFilePath)
        
        newImage = Atoms([getElement(i).number for i in atomsymbols])
        newImage.set_positions(atomcoords)
        initial = newImage.copy()
        
        atomsymbols, atomcoords = self.geometry.parseLOG(pGeom.getFilePath(self.outputFileExtension))
        
        newImage = Atoms([getElement(i).number for i in atomsymbols])
        newImage.set_positions(atomcoords)
        final = newImage.copy()
        
        return initial, final
            
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
                self.createInputFile(attempt)
                self.createIRCFile()
                success = self.run()
                if success:
                    logging.info('Attempt {0} of {1} on species {2} succeeded.'.format(attempt, self.maxAttempts, self.molecule.toAugmentedInChI()))
                    break
            else:
                logging.error('QM thermo calculation failed for {0}.'.format(self.molecule.toAugmentedInChI()))
                return None
        result = self.parse() # parsed in cclib
        return result

    def runDouble(self, inputFilePath):
        self.testReady()
        with open(inputFilePath) as infile:
            print "Running GAUSSIAN input file {0!s}:".format(inputFilePath)
            for line in infile:
                print line.rstrip()
        # submits the input file to Gaussian
        process = Popen([self.executablePath, inputFilePath])
        process.communicate()# necessary to wait for executable termination!
        
        logFilePath = os.path.splitext(inputFilePath)[0]+self.outputFileExtension
        with open(logFilePath) as outfile:
            print "Gaussian output file {0!s}:".format(logFilePath)
            for line in outfile:
                print line.rstrip()
        return logFilePath
        
    def runQST2(self):
        self.testReady()
        with open(self.inputFilePath) as infile:
            print "Running GAUSSIAN input file {0!s}:".format(self.inputFilePath)
            for line in infile:
                print line.rstrip()
        # submits the input file to Gaussian
        process = Popen([self.executablePath, self.inputFilePath])
        process.communicate()# necessary to wait for executable termination!
        
        logFilePath = self.outputFilePath
        with open(logFilePath) as outfile:
            print "Gaussian output file {0!s}:".format(logFilePath)
            for line in outfile:
                print line.rstrip()
        return self.verifyQST2OutputFile(), logFilePath
    

        
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
            
    def verifyQST2OutputFile(self):
        """
        Check's that a qst2 output file exists and was successful.
        
        Returns a boolean flag that states whether a QST2 GAUSSIAN simulation that hasn't failed already exists for the molecule with the 
        given file name.
        
        This checks that the calculation was attempted, only checking that a QST3 calculation is not required.
        
        If QST3 is required, False will be returned and the double-ended procedure has failed for the reaction.
        """
        
        failureKeys = [
                       '***** Convergence failure in GTrans *****',
                       'Try using 3 structures as input for',
                       ]
        
        if not os.path.exists(self.outputFilePath):
            logging.info("Output file {0} does not exist.".format(self.outputFilePath))
            return False
        
        # Initialize dictionary with "False"s 
        failureKeysFound = dict([(key, False) for key in failureKeys])
        
        with open(self.outputFilePath) as outputFile:
            for line in outputFile:
                line = line.strip()
                
                for element in failureKeys: #search for failure keywords
                    if element in line:
                        logging.error("Gaussian output file contains the following error: {0}".format(element) )
                        failureKeysFound[element] = True
        
        if any(failureKeysFound.values()):
            return False
        else:
            return True
    
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
            molfile = self.getFilePath('.crude.mol')
            atomline = re.compile('\s*([\- ][0-9.]+\s+[\-0-9.]+\s+[\-0-9.]+)\s+([A-Za-z]+)')
            
            atomCount = 0
            atomnosPrep = []
            with open(molfile) as molinput:
                for line in molinput:
                    match = atomline.match(line)
                    if match:
                        atomnosPrep.append(match.group(2))
                        atomCount += 1
            
            atomnos = []
            for atom in atomnosPrep:
                if atom == 'H':
                    atomnos.append(1)
                elif atom == 'C':
                    atomnos.append(6)
                elif atom == 'O':
                    atomnos.append(8)
            
            atomcoords = ircParse.atomcoords
            atomnos = numpy.array(atomnos)
            # Convert the IRC geometries into RMG molecules
            # We don't know which is reactant or product, so take the two at the end of the
            # paths and compare to the reactants and products
            mol1 = Molecule()
            mol1.fromXYZ(atomnos, atomcoords[pth1End])
            mol2 = Molecule()
            mol2.fromXYZ(atomnos, atomcoords[-1])
            
            testReaction = Reaction(
                                    reactants = mol1.split(),
                                    products = mol2.split(),                     
                                    )
            
            if isinstance(self.reaction.reactants[0], Molecule):
                targetReaction = Reaction(
                                        reactants = [reactant.toSingleBonds() for reactant in self.reaction.reactants],
                                        products = [product.toSingleBonds() for product in self.reaction.products],
                                        )
            elif isinstance(self.reaction.reactants[0], Species):
                targetReaction = Reaction(
                                        reactants = [reactant.molecule[0].toSingleBonds() for reactant in self.reaction.reactants],
                                        products = [product.molecule[0].toSingleBonds() for product in self.reaction.products],
                                        )
                                                                      
            if targetReaction.isIsomorphic(testReaction):
                return True
            else:
                return False
    
    def testTSGeometry(self, reactant):
        """
        Generates a dictionary of distances between all labeled atoms pairs. This may help identify
        transition state geometries with a small imaginary frequency that is NOT a reaction saddle
        point. This would be valueable as it would prevent the need for an IRC in these cases, saving
        on the computational cost.
        """
        labeledList = reactant.getLabeledAtoms()
        labeledAtoms = labeledList.values()
        labels = labeledList.keys()
        atomSymbols, atomCoords = self.geometry.parseLOG(self.outputFilePath)
        
        distances = {}
        dist_combo_it = itertools.combinations(labels, 2)
        dist_combo_l = list(dist_combo_it)
        for combo in dist_combo_l[:len(labels)]:
            atom1 = reactant.getLabeledAtom(combo[0])
            atom2 = reactant.getLabeledAtom(combo[1])
            
            coords1 = atomCoords[atom1.sortingLabel]
            coords2 = atomCoords[atom2.sortingLabel]
            
            distances[combo] = self.geometry.getDistance(coords1, coords2)
    
    def checkGeometry(self, outputFilePath, molecule):
        """
        Takes an output file, and extracts the geometry from it to ensure it has the same structure as the molecule provided.
        """
        parser = cclib.parser.Gaussian(outputFilePath)
        parser.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
        cclibData = parser.parse()
        
        mol = Molecule()
        mol.fromXYZ(cclibData.atomnos, cclibData.atomcoords[-1])
                                                        
        if mol.isIsomorphic(molecule.toSingleBonds()):
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
        if isinstance(self.reaction.reactants[0], Molecule):
            for molecule in self.reaction.reactants:
                radicalNumber += molecule.getRadicalCount()
        elif isinstance(self.reaction.reactants[0], Species):
            for molecule in self.reaction.reactants:
                radicalNumber += molecule.molecule[0].getRadicalCount() 
        
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
    
    def writeRxnOutputFile(self, labels, doubleEnd=False):                
        atomDist = self.parseTS(labels)
        
        distances = {'d12':float(atomDist[0]), 'd23':float(atomDist[1]), 'd13':float(atomDist[2])}
        user = "Pierre Bhoorasingh <bhoorasingh.p@husky.neu.edu>"
        if doubleEnd:
            description = "Found via double-ended search by the automatic transition state generator"
        else:
            description = "Found via group additive estimation by the automatic transition state generator"
        entry = Entry(
            index = 1,
            item = self.reaction,
            data = DistanceData(distances=distances, method='B3LYP/6-31+G(d,p)'),
            shortDesc = "B3LYP/6-31+G(d,p) calculation via group additive TS generator.",
        )
        
        outputDataFile = os.path.join(self.settings.fileStore, self.uniqueID + '.data')
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
                "#  m062x/gen opt=(ts,calcall,tight,noeigentest)  int=ultrafine nosymm",
                "#  m062x/gen opt=(ts,calcall,tight,noeigentest,cartesian) int=ultrafine geom=allcheck guess=check nosymm",
                "#  m062x/gen opt=(ts,calcall,noeigentest,cartesian) nosymm geom=allcheck guess=check nosymm",
                "#  m062x/gen opt=(ts,calcall,noeigentest) nosymm",
                "#  m062x/gen irc=(calcall,report=read) geom=allcheck guess=check nosymm",
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
                "# b3lyp/6-31+g(d,p) opt=(ts,calcfc,noeigentest) freq", # nosymm
                "# b3lyp/6-31+g(d,p) opt=(ts,calcfc,noeigentest,cartesian) freq", # nosymm geom=allcheck guess=check
                "# b3lyp/6-31+g(d,p) opt=(ts,calcfc,noeigentest) freq nosymm geom=allcheck guess=read",
                "# b3lyp/6-31+g(d,p) opt=(ts,calcfc,noeigentest,cartesian) freq nosymm geom=allcheck guess=check",
                "# b3lyp/6-31+g(d,p) irc=(calcall,report=read) freq geom=allcheck guess=check nosymm",
                "# b3lyp/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest) freq int=ultrafine nosymm",
                "# b3lyp/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest,cartesian) freq int=ultrafine geom=allcheck guess=check nosymm",
                
                ]
               # "# b3lyp/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest) int=ultrafine nosymm",
               # "# b3lyp/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest,cartesian) int=ultrafine geom=allcheck guess=check nosymm",
               # "# b3lyp/6-31+g(d,p) opt=(ts,calcall,noeigentest) nosymm",
               # "# b3lyp/6-31+g(d,p) opt=(ts,calcall,noeigentest,cartesian) nosymm geom=allcheck guess=check nosymm",
               # "# b3lyp/6-31+g(d,p) irc=(calcall,report=read) geom=allcheck guess=check nosymm",
               # ]
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
    
    def setCalculator(self, images):
        """
        Set up the Gaussian calculator for the Atomic Simulation Environment
        """
        import ase
        from ase.calculators.gaussian import Gaussian
        
        label=os.path.join(os.path.abspath(self.settings.fileStore), 'g09')
        for image in images[1:len(images)-1]:
            image.set_calculator(ase.calculators.gaussian.Gaussian(command=self.executablePath + '< ' + label + '.com' + ' > ' + label + '.log', label=label))
            image.get_calculator().set(multiplicity=self.geometry.molecule.getRadicalCount() + 1, method='b3lyp', basis='6-31+g(d,p)')
    
    def createGeomInputFile(self, freezeAtoms, otherGeom=None):
        
        output = [ '', self.geometry.uniqueIDlong, '', "{charge}   {mult}".format(charge=0, mult=(self.geometry.molecule.getRadicalCount() + 1) ) ]
        
        atomline = re.compile('\s*([\- ][0-9.]+\s+[\-0-9.]+\s+[\-0-9.]+)\s+([A-Za-z]+)')
        
        if otherGeom:
            molfile = otherGeom.getRefinedMolFilePath() # Now get the product geometry
            
            atomCount = 0
            with open(molfile) as molinput:
                for line in molinput:
                    match = atomline.match(line)
                    if match:
                        output.append("{0:8s} {1}".format(match.group(2), match.group(1)))
                        atomCount += 1
            inputFilePath = otherGeom.getFilePath(self.inputFileExtension)
            bottom_keys = "{atom1} {atom3} F\n{atom1} {atom2} F\n".format(atom1=freezeAtoms[0] + 1, atom2=freezeAtoms[1] + 1, atom3=freezeAtoms[2] + 1)
        else:
            molfile = self.geometry.getRefinedMolFilePath() # Get the reactant geometry
            
            atomCount = 0
            with open(molfile) as molinput:
                for line in molinput:
                    match = atomline.match(line)
                    if match:
                        output.append("{0:8s} {1}".format(match.group(2), match.group(1)))
                        atomCount += 1
            inputFilePath = self.inputFilePath
            bottom_keys = "{atom1} {atom3} F\n{atom2} {atom3} F\n".format(atom1=freezeAtoms[0] + 1, atom2=freezeAtoms[1] + 1, atom3=freezeAtoms[2] + 1)
        
        assert atomCount == len(self.geometry.molecule.atoms)
        
        output.append('')
        top_keys = "#  b3lyp/6-31+g(d,p) opt=(modredundant,MaxCycles={N}) nosymm\n".format(N=max(100,atomCount*10))
        self.writeInputFile(output, top_keys=top_keys, numProcShared=20, memory='800MB', bottomKeys=bottom_keys)
    
    def createQST2InputFile(self, pGeom):
        # For now we don't do this, until seg faults are fixed on Discovery.
        # chk_file = '%chk=' + os.path.join(self.settings.fileStore, self.uniqueID) + '\n'
        output = ['', self.geometry.uniqueID, '' ]
        output.append("{charge}   {mult}".format(charge=0, mult=(self.geometry.molecule.getRadicalCount() + 1) ))
        
        molfile = self.geometry.getRefinedMolFilePath()
        atomline = re.compile('\s*([\- ][0-9.]+\s+[\-0-9.]+\s+[\-0-9.]+)\s+([A-Za-z]+)')
        
        atomCount = 0
        atomnos = []
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    atomnos.append(match.group(2))
                    atomCount += 1
        
        assert atomCount == len(self.geometry.molecule.atoms)
        
        parser = cclib.parser.Gaussian(self.outputFilePath)
        parser.logger.setLevel(logging.ERROR)
        parser = parser.parse()
        
        lines = []
        for line in parser.atomcoords[-1]:
            lineList = ''
            for item in line:
                if item > 0:
                    lineList = lineList + ' ' + str(format(item, '.6f')) + ' '
                else:
                    lineList = lineList + str(format(item, '.6f')) + ' '
            lines.append(lineList)
            
        atomCount = 0
        
        for i, line in enumerate(lines):
            output.append("{0:8s} {1}".format(atomnos[i], line))
            atomCount += 1
        
        assert atomCount == len(self.geometry.molecule.atoms)
        
        output.append('')
        output.append(pGeom.uniqueIDlong)
        output.append('')
        output.append("{charge}   {mult}".format(charge=0, mult=(pGeom.molecule.getRadicalCount() + 1) ))
        
        parser = cclib.parser.Gaussian(pGeom.getFilePath(self.outputFileExtension))
        parser.logger.setLevel(logging.ERROR)
        parser = parser.parse()
        
        lines = []
        for line in parser.atomcoords[-1]:
            lineList = ''
            for item in line:
                if item > 0:
                    lineList = lineList + ' ' + str(format(item, '.6f')) + ' '
                else:
                    lineList = lineList + str(format(item, '.6f')) + ' '
            lines.append(lineList)
            
        atomCount = 0
        for i, line in enumerate(lines):
            output.append("{0:8s} {1}".format(atomnos[i], line))
            atomCount += 1
        
        assert atomCount == len(self.geometry.molecule.atoms)
        # 
        # #####################
        # 
        # atomCount = 0
        # with open(molfile) as molinput:
        #     for line in molinput:
        #         match = atomline.match(line)
        #         if match:
        #             output.append("{0:8s} {1}".format(match.group(2), match.group(1)))
        #             atomCount += 1
        # 
        # assert atomCount == len(self.geometry.molecule.atoms)
        # 
        # output.append('')
        # output.append(pGeom.uniqueIDlong)
        # output.append('')
        # output.append("{charge}   {mult}".format(charge=0, mult=(pGeom.molecule.getRadicalCount() + 1) ))
        # 
        # molfile = pGeom.getRefinedMolFilePath() # Now get the product geometry
        # 
        # atomCount = 0
        # with open(molfile) as molinput:
        #     for line in molinput:
        #         match = atomline.match(line)
        #         if match:
        #             output.append("{0:8s} {1}".format(match.group(2), match.group(1)))
        #             atomCount += 1
        # 
        # assert atomCount == len(self.geometry.molecule.atoms)
        
        output.append('')
        top_keys = "#  b3lyp/6-31+g(d,p) opt=(qst2,calcall,noeigentest,MaxCycles={N}) nosymm\n".format(N=max(100,atomCount*10))
        self.writeInputFile(output, top_keys=top_keys, numProcShared=20, memory='800MB')

class GaussianTSPM6(GaussianTS):

    #: Keywords that will be added at the top of the qm input file
    keywords = [
                "#  pm6 opt=(ts,calcfc,noeigentest) freq", # nosymm
                "#  pm6 opt=(ts,calcfc,noeigentest,cartesian) freq", # nosymm geom=allcheck guess=check
                "#  pm6 opt=(ts,calcfc,noeigentest) freq nosymm geom=allcheck guess=read",
                "#  pm6 opt=(ts,calcfc,noeigentest,cartesian) freq nosymm geom=allcheck guess=check",
                "#  pm6 irc=(calcall,report=read) geom=allcheck guess=check nosymm",
                "#  pm6 opt=(ts,calcall,tight,noeigentest) freq int=ultrafine nosymm",
                "#  pm6 opt=(ts,calcall,tight,noeigentest,cartesian) freq int=ultrafine geom=allcheck guess=check nosymm",

                ]
               # "# b3lyp/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest) int=ultrafine nosymm",
               # "# b3lyp/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest,cartesian) int=ultrafine geom=allcheck guess=check nosymm",
               # "# b3lyp/6-31+g(d,p) opt=(ts,calcall,noeigentest) nosymm",
               # "# b3lyp/6-31+g(d,p) opt=(ts,calcall,noeigentest,cartesian) nosymm geom=allcheck guess=check nosymm",
               # "# b3lyp/6-31+g(d,p) irc=(calcall,report=read) geom=allcheck guess=check nosymm",
               # ]
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
    
    def setCalculator(self, images):
        """
        Set up the Gaussian calculator for the Atomic Simulation Environment
        """
        import ase
        from ase.calculators.gaussian import Gaussian
        
        for image in images[1:len(images)+1]:
            image.set_calculator(ase.calculators.gaussian.Gaussian(command=self.executablePath + '< g09.com > g09.log'))
            image.get_calculator().set(multiplicity=self.geometry.molecule.getRadicalCount() + 1, method='pm6', basis='')
    
    def createGeomInputFile(self, freezeAtoms, otherGeom=None):
        
        output = [ '', self.geometry.uniqueIDlong, '', "{charge}   {mult}".format(charge=0, mult=(self.geometry.molecule.getRadicalCount() + 1) ) ]
        
        atomline = re.compile('\s*([\- ][0-9.]+\s+[\-0-9.]+\s+[\-0-9.]+)\s+([A-Za-z]+)')
        
        if otherGeom:
            molfile = otherGeom.getRefinedMolFilePath() # Now get the product geometry
            
            atomCount = 0
            with open(molfile) as molinput:
                for line in molinput:
                    match = atomline.match(line)
                    if match:
                        output.append("{0:8s} {1}".format(match.group(2), match.group(1)))
                        atomCount += 1
            inputFilePath = otherGeom.getFilePath(self.inputFileExtension)
            bottom_keys = "{atom1} {atom3} F\n{atom1} {atom2} F\n".format(atom1=freezeAtoms[0] + 1, atom2=freezeAtoms[1] + 1, atom3=freezeAtoms[2] + 1)
        else:
            molfile = self.geometry.getRefinedMolFilePath() # Get the reactant geometry
            
            atomCount = 0
            with open(molfile) as molinput:
                for line in molinput:
                    match = atomline.match(line)
                    if match:
                        output.append("{0:8s} {1}".format(match.group(2), match.group(1)))
                        atomCount += 1
            inputFilePath = self.inputFilePath
            bottom_keys = "{atom1} {atom3} F\n{atom2} {atom3} F\n".format(atom1=freezeAtoms[0] + 1, atom2=freezeAtoms[1] + 1, atom3=freezeAtoms[2] + 1)
        
        assert atomCount == len(self.geometry.molecule.atoms)
        
        output.append('')
        input_string = '\n'.join(output) + '\n'
        top_keys = "#  pm6 opt=(modredundant,MaxCycles={N}) nosymm\n".format(N=max(100,atomCount*10))
        
        with open(inputFilePath, 'w') as gaussianFile:
            # gaussianFile.write(numProc)
            # gaussianFile.write(mem)
            # gaussianFile.write(chk_file)
            gaussianFile.write(top_keys)
            gaussianFile.write(input_string)
            gaussianFile.write(bottom_keys)
            gaussianFile.write('\n')
    
    def createQST2InputFile(self, pGeom):
        # numProc = '%nprocshared=' + '20' + '\n' # could be something that is set in the qmSettings
        # mem = '%mem=' + '800MB' + '\n' # could be something that is set in the qmSettings
        # For now we don't do this, until seg faults are fixed on Discovery.
        # chk_file = '%chk=' + os.path.join(self.settings.fileStore, self.uniqueID) + '\n'
        output = ['', self.geometry.uniqueID, '' ]
        output.append("{charge}   {mult}".format(charge=0, mult=(self.geometry.molecule.getRadicalCount() + 1) ))
        
        molfile = self.geometry.getRefinedMolFilePath()
        atomline = re.compile('\s*([\- ][0-9.]+\s+[\-0-9.]+\s+[\-0-9.]+)\s+([A-Za-z]+)')
        
        atomCount = 0
        atomnos = []
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    atomnos.append(match.group(2))
                    atomCount += 1
        
        assert atomCount == len(self.geometry.molecule.atoms)
        
        parser = cclib.parser.Gaussian(self.outputFilePath)
        parser.logger.setLevel(logging.ERROR)
        parser = parser.parse()
        
        lines = []
        for line in parser.atomcoords[-1]:
            lineList = ''
            for item in line:
                if item > 0:
                    lineList = lineList + ' ' + str(format(item, '.6f')) + ' '
                else:
                    lineList = lineList + str(format(item, '.6f')) + ' '
            lines.append(lineList)
            
        atomCount = 0
        
        for i, line in enumerate(lines):
            output.append("{0:8s} {1}".format(atomnos[i], line))
            atomCount += 1
        
        assert atomCount == len(self.geometry.molecule.atoms)
        
        output.append('')
        output.append(pGeom.uniqueIDlong)
        output.append('')
        output.append("{charge}   {mult}".format(charge=0, mult=(pGeom.molecule.getRadicalCount() + 1) ))
        
        parser = cclib.parser.Gaussian(pGeom.getFilePath(self.outputFileExtension))
        parser.logger.setLevel(logging.ERROR)
        parser = parser.parse()
        
        lines = []
        for line in parser.atomcoords[-1]:
            lineList = ''
            for item in line:
                if item > 0:
                    lineList = lineList + ' ' + str(format(item, '.6f')) + ' '
                else:
                    lineList = lineList + str(format(item, '.6f')) + ' '
            lines.append(lineList)
            
        atomCount = 0
        for i, line in enumerate(lines):
            output.append("{0:8s} {1}".format(atomnos[i], line))
            atomCount += 1
        
        assert atomCount == len(self.geometry.molecule.atoms)
        # 
        # #####################
        # 
        # atomCount = 0
        # with open(molfile) as molinput:
        #     for line in molinput:
        #         match = atomline.match(line)
        #         if match:
        #             output.append("{0:8s} {1}".format(match.group(2), match.group(1)))
        #             atomCount += 1
        # 
        # assert atomCount == len(self.geometry.molecule.atoms)
        # 
        # output.append('')
        # output.append(pGeom.uniqueIDlong)
        # output.append('')
        # output.append("{charge}   {mult}".format(charge=0, mult=(pGeom.molecule.getRadicalCount() + 1) ))
        # 
        # molfile = pGeom.getRefinedMolFilePath() # Now get the product geometry
        # 
        # atomCount = 0
        # with open(molfile) as molinput:
        #     for line in molinput:
        #         match = atomline.match(line)
        #         if match:
        #             output.append("{0:8s} {1}".format(match.group(2), match.group(1)))
        #             atomCount += 1
        # 
        # assert atomCount == len(self.geometry.molecule.atoms)
        
        output.append('')
        input_string = '\n'.join(output) + '\n'
        top_keys = "#  pm6 opt=(qst2,calcall,noeigentest,MaxCycles={N}) nosymm\n".format(N=max(100,atomCount*10))
        
        with open(self.inputFilePath, 'w') as gaussianFile:
            # gaussianFile.write(numProc)
            # gaussianFile.write(mem)
            # gaussianFile.write(chk_file)
            gaussianFile.write(top_keys)
            gaussianFile.write(input_string)
            gaussianFile.write('\n')
