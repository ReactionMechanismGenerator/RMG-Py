import os
import re
import external.cclib as cclib
import logging
from subprocess import Popen, PIPE

from qmdata import CCLibData
from molecule import QMMolecule

class Mopac:
    """
    A base class for all QM calculations that use MOPAC.
    
    Classes such as :class:`MopacMol` will inherit from this class.
    """

    inputFileExtension = '.mop'
    outputFileExtension = '.out'
    
    mopacEnv = os.getenv('MOPAC_DIR', default="/opt/mopac")
    if os.path.exists(os.path.join(mopacEnv , 'MOPAC2012.exe')):
        executablePath = os.path.join(mopacEnv , 'MOPAC2012.exe')
        logging.debug("{0} is found.".format(executablePath))
    elif os.path.exists(os.path.join(mopacEnv , 'MOPAC2009.exe')):
        executablePath = os.path.join(mopacEnv , 'MOPAC2009.exe')
        logging.debug("{0} is found.".format(executablePath))
    else:
        executablePath = os.path.join(mopacEnv , '(MOPAC 2009 or 2012)')
        logging.debug("{0} is found.".format(executablePath))
    
    usePolar = False #use polar keyword in MOPAC
    
    "Keywords for the multiplicity"
    multiplicityKeywords = {
                            1: '',
                            2: 'uhf doublet',
                            3: 'uhf triplet',
                            4: 'uhf quartet',
                            5: 'uhf quintet',
                            6: 'uhf sextet',
                            7: 'uhf septet',
                            8: 'uhf octet',
                            9: 'uhf nonet',
                           }
    
    #: List of phrases that indicate failure
    #: NONE of these must be present in a succesful job.
    failureKeys = [
                   'IMAGINARY FREQUENCIES',
                   'EXCESS NUMBER OF OPTIMIZATION CYCLES',
                   'NOT ENOUGH TIME FOR ANOTHER CYCLE',
                   ]
    #: List of phrases to indicate success.
    #: ALL of these must be present in a successful job.
    successKeys = [
                   'DESCRIPTION OF VIBRATIONS',
                   'MOPAC DONE'
                  ]

    def testReady(self):
        if not os.path.exists(self.executablePath):
            logging.debug("{0} is not found.").format(self.executablePath)
            raise Exception("Couldn't find MOPAC executable at {0}. Try setting your MOPAC_DIR environment variable.".format(self.executablePath))

    def run(self):
        self.testReady()
        # submits the input file to mopac
        process = Popen([self.executablePath, self.inputFilePath])
        process.communicate()# necessary to wait for executable termination!
    
        return self.verifyOutputFile()
        
    def verifyOutputFile(self):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful MOPAC simulation already exists for the molecule with the 
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
            logging.debug("Output file {0} does not (yet) exist.".format(self.outputFilePath))
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
                        logging.error("MOPAC output file contains the following error: {0}".format(element) )
                        return False
                    
                for element in self.successKeys: #search for success keywords
                    if element in line:
                        successKeysFound[element] = True
               
                if "InChI=" in line:
                    logFileInChI = line #output files should take up to 240 characters of the name in the input file
                    InChIFound = True
                    if logFileInChI == self.uniqueIDlong:
                        InChIMatch = True
                    else:
                        logging.warning("InChI in log file ({0}) didn't match that in geometry ({1}).".format(logFileInChI, self.uniqueIDlong))                    
                        # Use only up to first 80 characters to match due to MOPAC bug which deletes 81st character of InChI string
                        if self.uniqueIDlong.startswith(logFileInChI[:80]):
                            logging.warning("but the beginning matches so it's probably just a truncation problem.")
                            InChIMatch = True
        # Check that ALL 'success' keywords were found in the file.
        if not all( successKeysFound.values() ):
            logging.error('Not all of the required keywords for success were found in the output file!')
            return False
        
        if not InChIFound:
            logging.error("No InChI was found in the MOPAC output file {0}".format(self.outputFilePath))
            return False
        
        if InChIMatch:
            logging.info("Successful MOPAC quantum result found in {0}".format(self.outputFilePath))
            # " + self.molfile.name + " ("+self.molfile.InChIAug+") has been found. This log file will be used.")
            return True
        
        #InChIs do not match (most likely due to limited name length mirrored in log file (240 characters), but possibly due to a collision)
        return self.checkForInChiKeyCollision(logFileInChI) # Not yet implemented!


    def parse(self):
        """
        Parses the results of the Mopac calculation, and returns a CCLibData object.
        """
        parser = cclib.parser.Mopac(self.outputFilePath)
        parser.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
        cclibData = parser.parse()
        radicalNumber = sum([i.radicalElectrons for i in self.molecule.atoms])
        qmData = CCLibData(cclibData, radicalNumber+1)
        return qmData

class MopacMol(QMMolecule, Mopac):
    """
    A base Class for calculations of molecules using MOPAC. 
    
    Inherits from both :class:`QMMolecule` and :class:`Mopac`.
    """

    def writeInputFile(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        
        molfile = self.getMolFilePathForCalculation(attempt) 
        atomline = re.compile('\s*([\- ][0-9.]+)\s+([\- ][0-9.]+)+\s+([\- ][0-9.]+)\s+([A-Za-z]+)')
        
        output = [ self.geometry.uniqueIDlong, '' ]
 
        atomCount = 0
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(4), match.group(1), match.group(2), match.group(3)))
                    atomCount += 1
        assert atomCount == len(self.molecule.atoms)
    
        output.append('')
        input_string = '\n'.join(output)
        
        top_keys, bottom_keys, polar_keys = self.inputFileKeywords(attempt)
        with open(self.inputFilePath, 'w') as mopacFile:
            mopacFile.write(top_keys)
            mopacFile.write('\n')
            mopacFile.write(input_string)
            mopacFile.write('\n')
            mopacFile.write(bottom_keys)
            if self.usePolar:
                mopacFile.write('\n\n\n')
                mopacFile.write(polar_keys)
                
    def inputFileKeywords(self, attempt):
        """
        Return the top, bottom, and polar keywords.
        """
        raise NotImplementedError("Should be defined by subclass, eg. MopacMolPM3")
        
class MopacMolPM3(MopacMol):

    #: Keywords that will be added at the top and bottom of the qm input file
    keywords = [
                {'top':"precise nosym", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0 nonr", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0 bfgs", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000", 'bottom':"oldgeo thermo nosym precise "},
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
        Return the top, bottom, and polar keywords for attempt number `attempt`.
        
        NB. `attempt`s begin at 1, not 0.
        """
        assert attempt <= self.maxAttempts
        
        if attempt > self.scriptAttempts:
            attempt -= self.scriptAttempts
        
        multiplicity_keys = self.multiplicityKeywords[self.geometry.multiplicity]

        top_keys = "pm3 {0} {1}".format(
                multiplicity_keys,
                self.keywords[attempt-1]['top'],
                )
        bottom_keys = "{0} pm3 {1}".format(
                self.keywords[attempt-1]['bottom'],
                multiplicity_keys,
                )
        polar_keys = "oldgeo {0} nosym precise pm3 {1}".format(
                'polar' if self.geometry.multiplicity == 1 else 'static',
                multiplicity_keys,
                )

        return top_keys, bottom_keys, polar_keys
    
class MopacMolPM6(MopacMol):

    #: Keywords that will be added at the top and bottom of the qm input file
    keywords = [
                {'top':"precise nosym", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0 nonr", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0 bfgs", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000", 'bottom':"oldgeo thermo nosym precise "},
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
        Return the top, bottom, and polar keywords for attempt number `attempt`.
        
        NB. `attempt`s begin at 1, not 0.
        """
        assert attempt <= self.maxAttempts
        
        if attempt > self.scriptAttempts:
            attempt -= self.scriptAttempts
        
        multiplicity_keys = self.multiplicityKeywords[self.geometry.multiplicity]

        top_keys = "pm6 {0} {1}".format(
                multiplicity_keys,
                self.keywords[attempt-1]['top'],
                )
        bottom_keys = "{0} pm6 {1}".format(
                self.keywords[attempt-1]['bottom'],
                multiplicity_keys,
                )
        polar_keys = "oldgeo {0} nosym precise pm6 {1}".format(
                'polar' if self.geometry.multiplicity == 1 else 'static',
                multiplicity_keys,
                )

        return top_keys, bottom_keys, polar_keys

class MopacMolPM7(MopacMol):

    #: Keywords that will be added at the top and bottom of the qm input file
    keywords = [
                {'top':"precise nosym", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0 nonr", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0 bfgs", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000", 'bottom':"oldgeo thermo nosym precise "},
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
        Return the top, bottom, and polar keywords for attempt number `attempt`.
        
        NB. `attempt`s begin at 1, not 0.
        """
        assert attempt <= self.maxAttempts
        
        if attempt > self.scriptAttempts:
            attempt -= self.scriptAttempts
        
        multiplicity_keys = self.multiplicityKeywords[self.geometry.multiplicity]

        top_keys = "pm7 {0} {1}".format(
                multiplicity_keys,
                self.keywords[attempt-1]['top'],
                )
        bottom_keys = "{0} pm7 {1}".format(
                self.keywords[attempt-1]['bottom'],
                multiplicity_keys,
                )
        polar_keys = "oldgeo {0} nosym precise pm7 {1}".format(
                'polar' if self.geometry.multiplicity == 1 else 'static',
                multiplicity_keys,
                )

        return top_keys, bottom_keys, polar_keys    
