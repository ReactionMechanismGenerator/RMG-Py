import os
import re
import external.cclib as cclib
import logging
from subprocess import Popen
import time
import math
import numpy

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
    
    def geomToString(self, atomsymbols, atomcoords, outputString='', freezeAtoms=[]):
        atomCount = 0
        for atomsymbol, atomcoord in zip(atomsymbols, atomcoords):
            if atomCount in freezeAtoms:
                inputline = "{0:4s} {1: .6f}  0  {2: .6f}  0  {3: .6f}  0".format(atomsymbol, atomcoord[0], atomcoord[1], atomcoord[2])
            else:
                inputline = "{0:4s} {1: .6f}  1  {2: .6f}  1  {3: .6f}  1".format(atomsymbol, atomcoord[0], atomcoord[1], atomcoord[2])
            outputString.append(inputline)
            atomCount += 1
        
        return outputString, atomCount
    
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
        5) checking that the optimized geometry, when connected by single bonds, is isomorphic with self.molecule (converted to single bonds)
        
        If any of the above criteria is not matched, False will be returned.
        If all succeed, then it will return True.
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
                    if self.uniqueIDlong in logFileInChI:
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
    
    def getParser(self, outputFile):
             """
             Returns the appropriate cclib parser.
             """
             return cclib.parser.Mopac(outputFile)
    
    def parse(self):
        """
        Returns the appropriate cclib parser.
        """
        return cclib.parser.Mopac(outputFile)
    
    def writeInputFile(self, output, attempt=None, top_keys=None, bottom_keys=None, inputFilePath=None, refFile=False):
        """
        Takes the output from the writeInputFile method and prints the
        file. Options provided allow the 
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        
        if not refFile:
            if not top_keys:
                top_keys, bottom_keys, polar_keys = self.inputFileKeywords(attempt)
            output = [top_keys] + output
        
        if bottom_keys:
            output = output + [bottom_keys]
        
        if self.usePolar:
            polar_keys = '\n\n\n' + polar_keys
            output = output + [polar_keys]
        
        input_string = '\n'.join(output)
        
        if not inputFilePath:
            inputFilePath = self.inputFilePath
        
        with open(inputFilePath, 'w') as mopacFile:
            mopacFile.write(input_string)

class MopacMol(QMMolecule, Mopac):
    """
    A base Class for calculations of molecules using MOPAC. 
    
    Inherits from both :class:`QMMolecule` and :class:`Mopac`.
    """

    #: Keywords that will be added at the top and bottom of the qm input file
    keywords = [
                {'top':"precise nosym", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0 nonr", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0 bfgs", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000", 'bottom':"oldgeo thermo nosym precise "},
                ]

    def createInputFile(self, attempt):
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
        self.writeInputFile(output, attempt)
                
    def inputFileKeywords(self, attempt):
        """
        Return the top, bottom, and polar keywords.
        """
        raise NotImplementedError("Should be defined by subclass, eg. MopacMolPM3")  
                  
    def generateQMData(self):
        """
        Calculate the QM data and return a QMData object, or None if it fails.
        """
        for atom in self.molecule.vertices:
            if atom.atomType.label in ('N5s', 'N5d', 'N5dd', 'N5t', 'N5b'):
                return None

        if self.verifyOutputFile():
            logging.info("Found a successful output file already; using that.")
            source = "QM {0} calculation found from previous run.".format(self.__class__.__name__)
        else:
            self.createGeometry()
            success = False
            for attempt in range(1, self.maxAttempts+1):
                self.createInputFile(attempt)
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


class MopacMolPMn(MopacMol):
    """
    Mopac PMn calculations for molecules (n undefined here)
    
    This is a parent class for MOPAC PMn calculations.
    Inherit it, and define the pm_method, then redefine 
    anything you wish to do differently.
    """
    pm_method = '(should be defined by sub class)'
    def inputFileKeywords(self, attempt):
        """
        Return the top, bottom, and polar keywords for attempt number `attempt`.
        
        NB. `attempt`s begin at 1, not 0.
        """
        assert attempt <= self.maxAttempts
        
        if attempt > self.scriptAttempts:
            attempt -= self.scriptAttempts
        
        multiplicity_keys = self.multiplicityKeywords[self.geometry.molecule.multiplicity]

        top_keys = "{method} {mult} {top}".format(
                method = self.pm_method,
                mult = multiplicity_keys,
                top = self.keywords[attempt-1]['top'],
                )
        bottom_keys = "{bottom} {method} {mult}".format(
                method = self.pm_method,
                bottom = self.keywords[attempt-1]['bottom'],
                mult = multiplicity_keys,
                )
        polar_keys = "oldgeo {polar} nosym precise {method} {mult}".format(
                method = self.pm_method,
                polar = ('polar' if self.geometry.molecule.multiplicity == 1 else 'static'),
                mult = multiplicity_keys,
                )

        return top_keys, bottom_keys, polar_keys

class MopacMolPM3(MopacMolPMn):
    """
    Mopac PM3 calculations for molecules
    
    This is a class of its own in case you wish to do anything differently,
    but for now it's the same as all the MOPAC PMn calculations, only pm3
    """
    pm_method = 'pm3'

class MopacMolPM6(MopacMolPMn):
    """
    Mopac PM6 calculations for molecules
    
    This is a class of its own in case you wish to do anything differently,
    but for now it's the same as all the MOPAC PMn calculations, only pm6
    """
    pm_method = 'pm6'

class MopacMolPM7(MopacMolPMn):
    """
    Mopac PM7 calculations for molecules
    
    This is a class of its own in case you wish to do anything differently,
    but for now it's the same as all the MOPAC PMn calculations, only pm7
    """
    pm_method = 'pm7'

    
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

    def setImages(self, pGeom):
        """
        Set and return the initial and final ase images for the NEB calculation
        """
        import ase
        from ase import Atoms

        # ASE doesn't keep the atoms in the same order as it's positions (weird!),
        # so get the correct atom list and recreate the images
        molfileR = self.getFilePath('.arc')
        molfileP = pGeom.getFilePath('.arc')
        atomline = re.compile('\s*([A-Za-z]+)\s+([\- ][0-9.]+)\s+([\+ ][0-9.]+)\s+([\- ][0-9.]+)\s+([\+ ][0-9.]+)\s+([\- ][0-9.]+)')

        atomCount = 0
        atomsymbols = []
        atomcoords = []
        with open(molfileR) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    atomsymbols.append(match.group(1))
                    atomcoords.append([float(match.group(2)), float(match.group(4)), float(match.group(6))])
                    atomCount += 1

        assert atomCount == len(self.geometry.molecule.atoms)

        newImage = Atoms(atomsymbols)
        newImage.set_positions(atomcoords)
        initial = newImage.copy()

        atomCount = 0
        atomsymbols = []
        atomcoords = []
        with open(molfileP) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    atomsymbols.append(match.group(1))
                    atomcoords.append([float(match.group(2)), float(match.group(4)), float(match.group(6))])
                    atomCount += 1

        assert atomCount == len(self.geometry.molecule.atoms)

        newImage = Atoms(atomsymbols)
        newImage.set_positions(atomcoords)
        final = newImage.copy()

        return initial, final

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

    def createInputFile(self, attempt):
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
        self.writeInputFile(output, attempt)

    def createGeomInputFile(self, freezeAtoms, product=False):

        if product:
            output = [ self.productGeom.uniqueID, '' ]
            freezeAtoms = [freezeAtoms[0], freezeAtoms[1]]
            atomsymbols, atomcoords = self.productGeom.parseMOL(self.productGeom.getRefinedMolFilePath())
            inputFilePath = self.productGeom.getFilePath(self.inputFileExtension)
            multiplicity = self.productGeom.molecule.multiplicity
        else:
            output = [ self.reactantGeom.uniqueID, '' ]
            freezeAtoms = [freezeAtoms[1], freezeAtoms[2]]
            atomsymbols, atomcoords = self.reactantGeom.parseMOL(self.reactantGeom.getRefinedMolFilePath())
            inputFilePath = self.reactantGeom.getFilePath(self.inputFileExtension)
            multiplicity = self.reactantGeom.molecule.multiplicity
            
        output, atomCount = self.geomToString(atomsymbols, atomcoords, outputString=output, freezeAtoms=freezeAtoms)

        assert atomCount == len(self.reactantGeom.molecule.atoms)
        
        output.append('')
        
        top_keys = "precise nosym {spin}".format(spin=self.multiplicityKeywords[multiplicity])
        self.writeInputFile(output, top_keys=top_keys, inputFilePath=inputFilePath)

    def createIRCFile(self):
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

    def createReferenceFile(self, productSide=False):#, inputFilePath, molFilePathForCalc, geometry, attempt, outputFile=None):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        if productSide:
            inputFilePath = self.productGeom.getFilePath(self.inputFileExtension)
            atomsymbols, atomcoords = self.productGeom.parseOUT(self.productGeom.getFilePath(self.outputFileExtension))
            output = [ '', self.productGeom.uniqueIDlong, '' ]
        else:
            inputFilePath = self.reactantGeom.getFilePath(self.inputFileExtension)
            atomsymbols, atomcoords = self.reactantGeom.parseARC(self.reactantGeom.getFilePath('.arc'))
            output = [ '', self.reactantGeom.uniqueIDlong, '' ]
        
        output, atomCount = self.geomToString(atomsymbols, atomcoords, outputString=output)

        assert atomCount == len(self.reactantGeom.molecule.atoms)

        output.append('')
        self.writeInputFile(output, inputFilePath=inputFilePath, refFile=True)

    def createGeoRefInputFile(self, productSide=False):
        if productSide:
            atomsymbols, atomcoords = self.productGeom.parseARC(self.productGeom.getFilePath('.arc'))
            refFile = self.reactantGeom.getFilePath(self.inputFileExtension)
            inputFilePath = self.productGeom.getFilePath(self.inputFileExtension)
            output = [ 'geo_ref="{0}"'.format(refFile), self.productGeom.uniqueIDlong, '' ]
        else:
            atomsymbols, atomcoords = self.reactantGeom.parseARC(self.reactantGeom.getFilePath('.arc'))
            refFile = self.productGeom.getFilePath(self.inputFileExtension)
            inputFilePath = self.reactantGeom.getFilePath(self.inputFileExtension)
            output = [ 'geo_ref="{0}"'.format(refFile), self.reactantGeom.uniqueIDlong, '' ]
        
        output, atomCount = self.geomToString(atomsymbols, atomcoords, outputString=output)

        assert atomCount == len(self.reactantGeom.molecule.atoms)

        output.append('')
        self.writeInputFile(output, inputFilePath=inputFilePath, refFile=True)

    def createSaddleInputFile(self):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        output = [ 'saddle', self.uniqueID, '' ]

        # Reactant side
        atomsymbols, atomcoords = self.reactantGeom.parseARC(self.reactantGeom.getFilePath('.arc'))
        
        output, atomCount = self.geomToString(atomsymbols, atomcoords, outputString=output)

        assert atomCount == len(self.reactantGeom.molecule.atoms)

        output.append('')

        # Product side
        atomsymbols, atomcoords = self.productGeom.parseARC(self.productGeom.getFilePath('.arc'))
        
        output, atomCount = self.geomToString(atomsymbols, atomcoords, outputString=output)

        assert atomCount == len(self.productGeom.molecule.atoms)

        output.append('')
        self.writeInputFile(output, refFile=True)

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
    
    def prepDoubleEnded(self, labels, notes):
        self.createGeomInputFile(freezeAtoms=labels)
        logFilePath = self.runDouble(self.reactantGeom.getFilePath(self.inputFileExtension))
        
        self.createGeomInputFile(freezeAtoms=labels, product=True)
        logFilePath = self.runDouble(self.productGeom.getFilePath(self.inputFileExtension))
        
        # A check is needed to ensure the geometry that comes out has not been altered
        return True, notes
    
    def conductDoubleEnded(self, NEB=False):
        if NEB:
            self.runNEB()
        else:
            # MOPAC Saddle Calculation
            self.createReferenceFile()
            self.createGeoRefInputFile(productSide=True)
            logFilePath = self.runDouble(self.productGeom.getFilePath(self.inputFileExtension))
            # shutil.copy(logFilePath, logFilePath+'.ref1.out')
                
            if not os.path.exists(self.productGeom.getFilePath('.arc')):
                notes = notes + 'product .arc file does not exits\n'
                return False, notes
            
            # Reactant that references the product geometry
            self.createReferenceFile(productSide=True)
            self.createGeoRefInputFile()
            logFilePath = self.runDouble(self.reactantGeom.getFilePath(self.inputFileExtension))
            # shutil.copy(logFilePath, logFilePath+'.ref2.out')
            
            if not os.path.exists(self.reactantGeom.getFilePath('.arc')):
                notes = notes + 'reactant .arc file does not exits\n'
                return False, notes
            
            # Write saddle calculation file using the outputs of the reference calculations
            self.createSaddleInputFile()
            self.runDouble(self.inputFilePath)
            
            atomSymbols, atomCoords = self.reactantGeom.parseOUT(self.outputFilePath)
            self.writeXYZ(atomSymbols, atomCoords)
            
            return True, notes # True, self.geometry, labels, notes
            
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

    def setCalculator(self, images):
        """
        Set up the Mopac calculator for the Atomic Simulation Environment
        """
        import ase
        from ase.calculators.mopac import Mopac

        label=os.path.join(os.path.abspath(self.settings.fileStore), 'ase')
        for image in images[1:len(images)-1]:
            image.set_calculator(ase.calculators.mopac.Mopac(command=self.executablePath, label=label, functional='PM7'))
            image.get_calculator().set(spin=self.geometry.molecule.getRadicalCount())

