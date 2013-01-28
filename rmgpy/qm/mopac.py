import os

import openbabel
import cclib.parser
import logging
from subprocess import Popen, PIPE

from qmdata import CCLibData
from molecule import QMMolecule
from reaction import QMReaction

class Mopac:
    """
    A base class for all QM calculations that use MOPAC.
    
    Classes such as :class:`MopacMol` will inherit from this class.
    """

    inputFileExtension = '.mop'
    outputFileExtension = '.out'
    executablePath = os.path.join(os.getenv('MOPAC_DIR', default="/opt/mopac") , 'MOPAC2012.exe')
    assert os.path.exists(executablePath), "Couldn't find MOPAC 2009 executable at {0}. Try setting your MOPAC_DIR environment variable.".format(executablePath)
    
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
                  ]


    def run(self):
        # submits the input file to mopac
        process = Popen([self.executablePath, self.inputFilePath])
        process.communicate()# necessary to wait for executable termination!
    
        return self.verifyOutputFile()
        
    def verifyOutputFile(self):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful MOPAC simulation already exists for the molecule with the 
        given (augmented) InChI Key.
        """
        
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
                        logging.info("InChI in log file didn't match that in geometry.")
        
        # Check that ALL 'success' keywords were found in the file.
        if not all( successKeysFound.values() ):
            logging.error('Not all of the required keywords for sucess were found in the output file!')
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
        
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "mop")
        mol = openbabel.OBMol()

        obConversion.ReadFile(mol, self.getMolFilePathForCalculation(attempt) )
        
        mol.SetTitle(self.geometry.uniqueIDlong)
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
        input_string = obConversion.WriteString(mol)
        top_keys, bottom_keys, polar_keys = self.inputFileKeywords(attempt)
        with open(self.inputFilePath, 'w') as mopacFile:
            mopacFile.write(top_keys)
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
        
    def generateQMData(self):
        """
        Calculate the QM data and return a QMData object, or None if it fails.
        """
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

#TS
class MopacTS(QMReaction, Mopac):
    #*****change this for TS
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
    keywordsTop[1] = "ts"
    keywordsTop[2] = "ts recalc=5"
    keywordsTop[3] = "ts ddmin=0.0001"
    keywordsTop[4] = "ts recalc=5 ddmin=0.0001"

    "Keywords that will be added at the bottom of the qm input file"
    keywordsBottom = {}
    keywordsBottom[1] = "oldgeo force vectors esp"
    keywordsBottom[2] = "oldgeo force vectors esp"
    keywordsBottom[3] = "oldgeo force vectors esp"
    keywordsBottom[4] = "oldgeo force vectors esp"

    scriptAttempts = len(keywordsTop)

    failureKeys = ['GRADIENT IS TOO LARGE', 
                'EXCESS NUMBER OF OPTIMIZATION CYCLES', 
                'NOT ENOUGH TIME FOR ANOTHER CYCLE',
                '6 IMAGINARY FREQUENCIES',
                '5 IMAGINARY FREQUENCIES',
                '4 IMAGINARY FREQUENCIES',
                '3 IMAGINARY FREQUENCIES',
                '2 IMAGINARY FREQUENCIES'
                ]

    def __init__(self, reaction):
        self.reaction = reaction
        self.reactants = reaction.reactants
        self.products = reaction.products
        self.family = reaction.family
        self.rdmol = None

    def generateTransitionState(self):
        """
        make TS geometry
        """
        if not os.path.exists(self.reaction.family.name):
            logging.info("Creating directory %s for mol files."%os.path.abspath(self.reaction.family.name))
            os.makedirs(self.reaction.family.name)
        inputFilePath = os.path.join(self.reaction.family.name, self.reactants[0].toAugmentedInChIKey())
        if os.path.exists(inputFilePath):
            inputFilePath = os.path.join(self.reaction.family.name, self.products[0].toAugmentedInChIKey())
            if os.path.exists(inputFilePath):
                inputFilePath = os.path.join(self.reaction.family.name, self.reactants[0].toAugmentedInChIKey() + self.products[0].toAugmentedInChIKey())
        with open(inputFilePath, 'w') as mopacFile:
            for reactant in self.reactants:
                mopacFile.write(reactant.toSMILES())
                mopacFile.write('\n')
                mopacFile.write(reactant.toAdjacencyList())
                mopacFile.write('\n')
            for product in self.products:
                mopacFile.write(product.toSMILES())
                mopacFile.write('\n')
                mopacFile.write(product.toAdjacencyList())
                mopacFile.write('\n')
        # if self.reaction.family.name.lower() == 'intra_r_add_exocyclic' or self.reaction.family.name.lower() == 'intra_r_add_endocyclic':
        #     rdMol, tsBM, mult, lbl, other = self.getTSBMatrix()
        #     self.geometry.uniqueID = self.reactants[0].toSMILES() + '_' + self.products[0].toSMILES()
        #     import copy
        #     initialID = copy.deepcopy(self.geometry.uniqueID)
        #     success = False
        #     check = 0
        #     self.geometry.rd_embed(rdMol, 1, tsBM)
        #     inputString = self.convertMolFile('mopin', 1, self.scriptAttempts)
        #     while not success and check <= 5:
        #         inputString = self.fixBond(inputString, lbl)
        #         check += 1
        #         attempt = 0
        #         while not success and attempt < self.scriptAttempts:
        #             attempt += 1
        #             self.geometry.uniqueID = initialID + str(check) + str(attempt)
        #             top_keys, bottom_keys, polar_keys = self.inputFileKeys(attempt, mult)
        #             inputFileName = self.writeInputFile(attempt, top_keys, bottom_keys, polar_keys, self.scriptAttempts, input_string=inputString)
        #             success = self.run(inputFileName)
        #     import ipdb; ipdb.set_trace()
        # else:
        #     pass

class MopacPM3(MopacTS):
    def inputFileKeys(self, attempt, multiplicity):
        """
        Inherits the writeInputFile methods from mopac.py
        """
        multiplicity_keys = self.multiplicityKeywords[multiplicity]

        top_keys = "pm3 {0} {1}".format(
                multiplicity_keys,
                self.keywordsTop[attempt],
                )
        bottom_keys = "{0} pm3 {1}".format(
                self.keywordsBottom[attempt],
                multiplicity_keys,
                )
        polar_keys = "oldgeo {0} nosym precise pm3 {1}".format(
                'polar' if multiplicity == 1 else 'static',
                multiplicity_keys,
                )

        return top_keys, bottom_keys, polar_keys
