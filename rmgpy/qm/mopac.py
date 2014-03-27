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
    elif os.path.exists(os.path.join(mopacEnv , 'MOPAC2009.exe')):
        executablePath = os.path.join(mopacEnv , 'MOPAC2009.exe')
    else:
        executablePath = os.path.join(mopacEnv , '(MOPAC 2009 or 2012)')
    
    usePolar = False#use polar keyword in MOPAC
    
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
            raise Exception("Couldn't find MOPAC executable at {0}. Try setting your MOPAC_DIR environment variable.".format(self.executablePath))

    def run(self):
        self.testReady()
        # submits the input file to mopac
        process = Popen([self.executablePath, self.inputFilePath])
        process.communicate()# necessary to wait for executable termination!
    
        return self.verifyOutputFile()
        
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
                
    def inputFileKeywords(self, attempt):
        """
        Return the top, bottom, and polar keywords.
        """
        raise NotImplementedError("Should be defined by subclass, eg. MopacMolPM3")
        
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
    
    def verifyOutputFile(self):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful MOPAC simulation already exists for the molecule with the 
        given (augmented) InChI Key.
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
            # Compare the optimized geometry to the original molecule
            parser = cclib.parser.Mopac(self.outputFilePath)
            parser.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
            cclibData = parser.parse()
            cclibMol = Molecule()
            cclibMol.fromXYZ(cclibData.atomnos, cclibData.atomcoords[-1])
            testMol = self.molecule.toSingleBonds()
            
            if cclibMol.isIsomorphic(testMol):
                logging.info("Successful MOPAC quantum result found in {0}".format(self.outputFilePath))
                # " + self.molfile.name + " ("+self.molfile.InChIAug+") has been found. This log file will be used.")
                return True
            else:
                logging.info("Incorrect connectivity for optimized geometry in file {0}".format(self.outputFilePath))
                # " + self.molfile.name + " ("+self.molfile.InChIAug+") has been found. This log file will be used.")
                return False
        
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

    def generateQMData(self):
        """
        Calculate the QM data and return a QMData object, or None if it fails.
        """
        for atom in self.molecule.vertices:
            if atom.atomType.label == 'N5s' or atom.atomType.label == 'N5d' or atom.atomType.label =='N5dd' or atom.atomType.label == 'N5t' or atom.atomType.label == 'N5b':
                return None

        if self.verifyOutputFile():
            logging.info("Found a successful output file already; using that.")
            source = "QM MOPAC result file found from previous run."
        else:
            self.createGeometry()
            success = False
            for attempt in range(1, self.maxAttempts+1):
                self.writeInputFile(attempt)
                logging.info('Trying {3} attempt {0} of {1} on molecule {2}.'.format(attempt, self.maxAttempts, self.molecule.toSMILES(), self.__class__.__name__))
                success = self.run()
                if success:
                    source = "QM {0} calculation attempt {1}".format(self.__class__.__name__, attempt )
                    break
            else:
                logging.error('QM thermo calculation failed for {0}.'.format(self.molecule.toAugmentedInChI()))
                return None
        result = self.parse() # parsed in cclib
        result.source = source
        return result # a CCLibData object


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

class MopacTS(QMReaction, Mopac):
    
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

    "Keywords that will be added at the top of the qm input file"
    keywords = [
                {'top':"ts", 'bottom':"oldgeo force let"},
                {'top':"ts", 'bottom':"oldgeo force esp "},
                {'top':"ts", 'bottom':"oldgeo force vectors "},
                {'top':"ts", 'bottom':"oldgeo force vectors esp "},
                ]
                
    scriptAttempts = len(keywords)

    failureKeys = ['GRADIENT IS TOO LARGE', 
                   'EXCESS NUMBER OF OPTIMIZATION CYCLES', 
                   'NOT ENOUGH TIME FOR ANOTHER CYCLE',
                   # '6 IMAGINARY FREQUENCIES',
                   # '5 IMAGINARY FREQUENCIES',
                   # '4 IMAGINARY FREQUENCIES',
                   # '3 IMAGINARY FREQUENCIES',
                   # '2 IMAGINARY FREQUENCIES'
                   ]
    
    def runDouble(self, inputFilePath):
        self.testReady()
        with open(inputFilePath) as infile:
            print "Running MOPAC input file {0!s}:".format(inputFilePath)
            for line in infile:
                print line.rstrip()
        # submits the input file to mopac
        process = Popen([self.executablePath, inputFilePath])
        process.communicate()# necessary to wait for executable termination!
        
        logFilePath = os.path.splitext(inputFilePath)[0]+self.outputFileExtension
        with open(logFilePath) as outfile:
            print "MOPAC output file {0!s}:".format(logFilePath)
            for line in outfile:
                print line.rstrip()
        return logFilePath
        
    def runIRC(self):
        self.testReady()
        # submits the input file to mopac
        process = Popen([self.executablePath, self.inputFilePath])
        process.communicate()# necessary to wait for executable termination!
    
        return self.verifyIRCOutputFile()
    
    def writeInputFile(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        
        molfile = self.geometry.getRefinedMolFilePath()
        atomline = re.compile('\s*([\- ][0-9.]+)\s+([\- ][0-9.]+)+\s+([\- ][0-9.]+)\s+([A-Za-z]+)')
        
        output = [ self.geometry.uniqueID, '' ]
        
        atomCount = 0
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(4), match.group(1), match.group(2), match.group(3)))
                    atomCount += 1
        assert atomCount == len(self.geometry.molecule.atoms)
        
        output.append('')
        input_string = '\n'.join(output)
        
        top_keys, bottom_keys, polar_keys = self.inputFileKeywords(attempt, self.geometry.molecule.getRadicalCount() + 1)
        with open(self.inputFilePath, 'w') as mopacFile:
            mopacFile.write(top_keys)
            mopacFile.write('\n')
            mopacFile.write(input_string)
            mopacFile.write('\n')
            mopacFile.write(bottom_keys)
            if self.usePolar:
                mopacFile.write('\n\n\n')
                mopacFile.write(polar_keys)
    
    def writeGeomInputFile(self, freezeAtoms, otherGeom=None):
        
        atomline = re.compile('\s*([\- ][0-9.]+)\s+([\- ][0-9.]+)+\s+([\- ][0-9.]+)\s+([A-Za-z]+)')
        
        output = [ self.geometry.uniqueID ]
        
        if otherGeom:
            freezeAtoms = [freezeAtoms[0], freezeAtoms[1]]
            molfile = otherGeom.getRefinedMolFilePath() # Now get the product geometry
            inputFilePath = otherGeom.getFilePath(self.inputFileExtension)
        else:
            freezeAtoms = [freezeAtoms[1], freezeAtoms[2]]
            molfile = self.geometry.getRefinedMolFilePath() # Get the reactant geometry
            inputFilePath = self.inputFilePath
        
        atomCount = 0
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    if atomCount in freezeAtoms:
                        output.append("{0:4s} {1} 0 {2} 0 {3} 0".format(match.group(4), match.group(1), match.group(2), match.group(3)))
                    else:
                        output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(4), match.group(1), match.group(2), match.group(3)))
                    atomCount += 1
        
        assert atomCount == len(self.geometry.molecule.atoms)
        
        output.append('')
        input_string = '\n'.join(output)
        
        top_keys = "precise nosym {spin}\n".format(spin=self.multiplicityKeywords[self.geometry.multiplicity])
        
        with open(inputFilePath, 'w') as mopacFile:
            mopacFile.write(top_keys)
            mopacFile.write('\n')
            mopacFile.write(input_string)
            mopacFile.write('\n')
                        
    def writeIRCFile(self):
        output = ['irc=1* let', self.geometry.uniqueID, '' ]
        atomCount = 0
        
        molfile = self.getFilePath('.arc')
        atomline = re.compile('\s*([A-Za-z]+)\s+([\- ][0-9.]+)\s+([\+ ][0-9.]+)\s+([\- ][0-9.]+)\s+([\+ ][0-9.]+)\s+([\- ][0-9.]+)')
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(1), match.group(2), match.group(4), match.group(6)))
                    atomCount += 1
        
        assert atomCount == len(self.geometry.molecule.atoms)
        
        output.append('')
        input_string = '\n'.join(output)
        
        with open(self.inputFilePath, 'w') as mopacFile:
            mopacFile.write(input_string)
        
    def writeReferenceFile(self, freezeAtoms, otherGeom=None):#, inputFilePath, molFilePathForCalc, geometry, attempt, outputFile=None):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        
        output = [ '', self.geometry.uniqueIDlong, '' ]
        atomline = re.compile('\s*([A-Za-z]+)\s+([\- ][0-9.]+)\s+([\+ ][0-9.]+)\s+([\- ][0-9.]+)\s+([\+ ][0-9.]+)\s+([\- ][0-9.]+)')
        
        if otherGeom:
            inputFilePath = otherGeom.getFilePath(self.inputFileExtension)
            molfile = otherGeom.getFilePath('.arc')
        else:
            inputFilePath = self.inputFilePath
            molfile = self.getFilePath('.arc')
            
        # freezeAtoms = [freezeAtoms[0], freezeAtoms[2]]
        atomCount = 0
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(1), match.group(2), match.group(4), match.group(6)))
                    # if atomCount in freezeAtoms:
                    #     output.append("{0:4s} {1} 0 {2} 0 {3} 0".format(match.group(1), match.group(2), match.group(4), match.group(6)))
                    # else:
                    #     output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(1), match.group(2), match.group(4), match.group(6)))
                    atomCount += 1
        
        assert atomCount == len(self.geometry.molecule.atoms)
        
        output.append('')
        input_string = '\n'.join(output)
        
        with open(inputFilePath, 'w') as mopacFile:
            mopacFile.write(input_string)
            mopacFile.write('\n')
    
    def writeGeoRefInputFile(self, otherGeom, freezeAtoms, otherSide=False):
        if otherSide:
            molfile = otherGeom.getFilePath('.arc')
            refFile = self.inputFilePath
            inputFilePath = otherGeom.getFilePath(self.inputFileExtension)
        else:
            molfile = self.getFilePath('.arc')
            refFile = otherGeom.getFilePath(self.inputFileExtension)
            inputFilePath = self.inputFilePath
            
        freezeAtoms = [freezeAtoms[0], freezeAtoms[2]]
        atomline = re.compile('\s*([A-Za-z]+)\s+([\- ][0-9.]+)\s+([\+ ][0-9.]+)\s+([\- ][0-9.]+)\s+([\+ ][0-9.]+)\s+([\- ][0-9.]+)')
        
        output = [ 'geo_ref="{0}"'.format(refFile), self.geometry.uniqueIDlong, '' ]
        
        atomCount = 0
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(1), match.group(2), match.group(4), match.group(6)))
                    # if atomCount in freezeAtoms:
                    #     output.append("{0:4s} {1} 0 {2} 0 {3} 0".format(match.group(1), match.group(2), match.group(4), match.group(6)))
                    # else:
                    #     output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(1), match.group(2), match.group(4), match.group(6)))
                    atomCount += 1
        assert atomCount == len(self.geometry.molecule.atoms)
        
        output.append('')
        input_string = '\n'.join(output)
        with open(inputFilePath, 'w') as mopacFile:
            mopacFile.write(input_string)
            mopacFile.write('\n')
    
    def writeSaddleInputFile(self, otherGeom):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        output = [ 'saddle', self.geometry.uniqueIDlong, '' ]
        atomline = re.compile('\s*([A-Za-z]+)\s+([\- ][0-9.]+)\s+([\+ ][0-9.]+)\s+([\- ][0-9.]+)\s+([\+ ][0-9.]+)\s+([\- ][0-9.]+)')
        
        # Reactant side
        molfile = self.getFilePath('.arc')
        atomCount = 0
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(1), match.group(2), match.group(4), match.group(6)))
                    atomCount += 1
        
        assert atomCount == len(self.geometry.molecule.atoms)
        atomCount = 0
        
        output.append('')
        
        # Product side
        molfile = otherGeom.getFilePath('.arc') 
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(1), match.group(2), match.group(4), match.group(6)))
                    atomCount += 1
        
        # atomline = re.compile('\s*([A-Za-z]+)\s+([\- ][0-9.]+)\s+([\+ ][0-9.]+)\s+([\- ][0-9.]+)\s+([\+ ][0-9.]+)\s+([\- ][0-9.]+)')
        # output = [ 'saddle', self.geometry.uniqueIDlong, '' ]
        # 
        # atomCount = 0
        # molfile = self.getFilePath('.arc')
        # with open(molfile) as molinput:
        #     for line in molinput:
        #         match = atomline.match(line)
        #         if match:
        #             output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(1), match.group(2), match.group(4), match.group(6)))
        #             atomCount += 1
        #             
        # assert atomCount == len(self.geometry.molecule.atoms)
        # atomCount = 0
        # 
        # output.append('')
        # 
        # molfile = otherGeom.getFilePath('.arc')
        # with open(molfile) as molinput:
        #     for line in molinput:
        #         match = atomline.match(line)
        #         if match:
        #             output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(1), match.group(2), match.group(4), match.group(6)))
        #             atomCount += 1
        
        assert atomCount == len(self.geometry.molecule.atoms)
        
        output.append('')
        input_string = '\n'.join(output)
    
        with open(self.inputFilePath, 'w') as mopacFile:
            mopacFile.write(input_string)
            mopacFile.write('\n')
    
    def parseArc(outputFile):
        """
        This can be reworked with regex to be more general
        """
        arcFile = outputFile.split('.')[0] + '.arc'
        geomLines = list()
        readFile = file(arcFile)
        for line in reversed(readFile.readlines()):
            if line.startswith('geo_ref'):
                break
            else:
                geomLines.append(line)
    
        geomLines.reverse()
    
        return geomLines
    
    def writeTSInputFile(self, doubleEnd=False):
        
        output = [ self.geometry.uniqueID, '' ]
        atomCount = 0
        
        if doubleEnd:
            molfile = self.outputFilePath
            atomline = re.compile('\s*([0-9]+)\s+([A-Za-z]+)\s+([\- ][0-9.]+)\s+([*])\s+([\- ][0-9.]+)\s+([*])\s+([\- ][0-9.]+)')
            with open(molfile) as molinput:
                total = len(self.geometry.molecule.atoms)
                matchLines = []
                for line in molinput:
                    match = atomline.match(line)
                    if match:
                        matchLines.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(2), match.group(3), match.group(5), match.group(7)))
                if len(matchLines) == total:
                    for line in matchLines:
                        output.append(line)
                        atomCount += 1
                else:
                    for line in matchLines[total*-1:]:
                        output.append(line)
                        atomCount += 1
        else:
            molfile = self.geometry.getRefinedMolFilePath()
            atomline = re.compile('\s*([\- ][0-9.]+)\s+([\- ][0-9.]+)+\s+([\- ][0-9.]+)\s+([A-Za-z]+)')
            with open(molfile) as molinput:
                for line in molinput:
                    match = atomline.match(line)
                    if match:
                        output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(4), match.group(1), match.group(2), match.group(3)))
                        atomCount += 1
        
        assert atomCount == len(self.geometry.molecule.atoms)
        
        output.append('')
        input_string = '\n'.join(output)
    
        top_keys = 'ts recalc=5\n'
        bottom_keys = 'oldgeo force let\n'
        with open(self.inputFilePath, 'w') as mopacFile:
            mopacFile.write(top_keys)
            mopacFile.write(input_string)
            mopacFile.write('\n')
            mopacFile.write(bottom_keys)
    
    def verifyOutputFile(self):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful MOPAC simulation already exists for the molecule with the 
        given (augmented) InChI Key.
        """
        if not os.path.exists(self.outputFilePath):
            logging.info("Output file {0} does not exist.".format(self.outputFilePath))
            return False, False
        
        # Initialize dictionary with "False"s 
        successKeysFound = dict([(key, False) for key in self.successKeys])
        
        with open(self.outputFilePath) as outputFile:
            for line in outputFile:
                line = line.strip()
                
                for element in self.failureKeys: #search for failure keywords
                    if element in line:
                        logging.error("MOPAC output file contains the following error: {0}".format(element) )
                        return False, False
                    
                for element in self.successKeys: #search for success keywords
                    if element in line:
                        successKeysFound[element] = True
                      
        # Check that ALL 'success' keywords were found in the file.
        if not all( successKeysFound.values() ):
            logging.error('Not all of the required keywords for success were found in the output file!')
            return False, False
        else:
            logging.info("Successful MOPAC quantum result found in {0}".format(self.outputFilePath))
            return True, False
        
        #InChIs do not match (most likely due to limited name length mirrored in log file (240 characters), but possibly due to a collision)
        return self.checkForInChiKeyCollision(logFileInChI) # Not yet implemented!
    
    def convertMol(self, geomLines):
        
        atomcoords = []
        atomnos = []
        for line in geomLines:
            atType, x, y, z = line.split()
            atomnos.append(getElement(atType).number)
            atomcoords.append([float(x),float(y),float(z)])
        atomnos = numpy.array(atomnos, dtype=int)
        atomcoords = numpy.array(atomcoords)
        mol = Molecule()
        mol.fromXYZ(atomnos, atomcoords)
        
        return mol
    
    def verifyIRCOutputFile(self):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful MOPAC simulation already exists for the molecule with the 
        given (augmented) InChI Key.
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
                        logging.error("MOPAC output file contains the following error: {0}".format(element) )
                        return False
                    
                for element in self.successKeys: #search for success keywords
                    if element in line:
                        successKeysFound[element] = True
        
        if not successKeysFound['MOPAC DONE']:
            logging.error('Not all of the required keywords for success were found in the IRC output file!')
            return False
        
        with open(self.outputFilePath.split('.')[0] + '.xyz') as geomFile:
            geomFile = geomFile.readlines()
            geomFile.pop(0)
            geomFile.pop(0)
            geom1 = []
            for line in geomFile:
                if not line.startswith('  reversed'):
                    geom1.append(line)
                else:
                    break
            geom1.pop()
            
            geom2 = []
            for line in reversed(geomFile):
                if not line.startswith(' DRC'):
                    geom2.append(line)
                else:
                    break
        
        mol1 = self.convertMol(geom1)
        mol2 = self.convertMol(geom2)
        
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
    
    def parseTS(self, labels):
        
        def getDistance(coordinates1, coordinates2):
            """
            Return the square of the distance (in Angstrom) between the two atoms.
            """
            diff = (coordinates1.coords - coordinates2.coords)
            return math.sqrt(sum(diff * diff))
        
        tsParse = cclib.parser.Mopac(os.path.join(self.file_store_path, self.uniqueID + self.outputFileExtension))
        tsParse = tsParse.parse()
        geom = tsParse.atomcoords[-1]
        atomNums = tsParse.atomnos
        
        atom1 = Atom(element=getElement(int(atomNums[labels[0]])), coords=geom[labels[0]])
        atom2 = Atom(element=getElement(int(atomNums[labels[1]])), coords=geom[labels[1]])
        atom3 = Atom(element=getElement(int(atomNums[labels[2]])), coords=geom[labels[2]])
        
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
            shortDesc = "B3LYP/6-31+G(d,p) calculation via group estimated TS generator.",
            history = [(time.asctime(), user, 'action', description)]
        )
        
        outputDataFile = os.path.join(self.file_store_path, self.uniqueID + '.data')
        with open(outputDataFile, 'w') as parseFile:
            saveEntry(parseFile, entry)
    
class MopacTSPM3(MopacTS):
    def inputFileKeywords(self, attempt, multiplicity):
        """
        Inherits the writeInputFile methods from mopac.py
        """
        multiplicity_keys = self.multiplicityKeywords[multiplicity]

        top_keys = "pm3 {0} {1}".format(
                multiplicity_keys,
                self.keywords[attempt-1]['top'],
                )
        bottom_keys = "{0} pm3 {1}".format(
                self.keywords[attempt-1]['bottom'],
                multiplicity_keys,
                )
        polar_keys = "oldgeo {0} nosym precise pm3 {1}".format(
                'polar' if multiplicity == 1 else 'static',
                multiplicity_keys,
                )

        return top_keys, bottom_keys, polar_keys

class MopacTSPM7(MopacTS):
    def inputFileKeywords(self, attempt, multiplicity):
        """
        Inherits the writeInputFile methods from mopac.py
        """
        multiplicity_keys = self.multiplicityKeywords[multiplicity]

        top_keys = "pm7 {0} {1}".format(
                multiplicity_keys,
                self.keywords[attempt-1]['top'],
                )
        bottom_keys = "{0} pm7 {1}".format(
                self.keywords[attempt-1]['bottom'],
                multiplicity_keys,
                )
        polar_keys = "oldgeo {0} nosym precise pm7 {1}".format(
                'polar' if multiplicity == 1 else 'static',
                multiplicity_keys,
                )

        return top_keys, bottom_keys, polar_keys
