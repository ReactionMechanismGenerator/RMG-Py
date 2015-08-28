import os
import logging
import math

import numpy

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    logging.debug("To use QM calculations you must correctly install rdkit.")
    
import rmgpy.quantity
from rmgpy.thermo import ThermoData
import rmgpy.statmech
import rmgpy.thermo
import rmgpy.molecule
import symmetry
import qmdata
from qmdata import CCLibData
from rmgpy.molecule import parser

class Geometry:
    """
    A geometry, used for quantum calculations.
    
    Created from a molecule. Geometry estimated by RDKit.
    
    The attributes are:
    
    =================== ======================= ====================================
    Attribute           Type                    Description
    =================== ======================= ====================================
    `settings`          :class:`QMSettings`     Settings for QM calculations
    `uniqueID`          ``str``                 A short ID such as an augmented InChI Key
    `molecule`          :class:`Molecule`       RMG Molecule object
    `uniqueIDlong`      ``str``                 A long, truly unique ID such as an augmented InChI
    =================== ======================= ====================================
    
    """
    def __init__(self, settings, uniqueID, molecule, uniqueIDlong=None):
        self.settings = settings
        #: A short unique ID such as an augmented InChI Key.
        self.uniqueID = uniqueID
        self.molecule = molecule
        if uniqueIDlong is None:
            self.uniqueIDlong = uniqueID
        else:
            #: Long, truly unique, ID, such as the augmented InChI.
            self.uniqueIDlong = uniqueIDlong
        
        # ToDo: why do we copy self.settings.fileStore into self.fileStore ?
        # (and same for .scratchDirectory)
        if self.settings:
            self.fileStore = self.settings.fileStore
            self.scratchDirectory = self.settings.scratchDirectory
        else:
            self.fileStore = None
            self.scratchDirectory = None
        
        if self.fileStore and not os.path.exists(self.fileStore):
            logging.info("Creating permanent directory %s for qm files."%os.path.abspath(self.fileStore))
            os.makedirs(self.fileStore)
            
        if self.scratchDirectory and not os.path.exists(self.scratchDirectory):
            logging.info("Creating scratch directory %s for qm files."%os.path.abspath(self.scratchDirectory))
            os.makedirs(self.scratchDirectory)

    def getFilePath(self, extension, scratch=True):
        """
        Returns the path to the file with the given extension.
        
        The provided extension should include the leading dot.
        If called with `scratch=False` then it will be in the `fileStore` directory,
        else `scratch=True` is assumed and it will be in the `scratchDirectory` directory.
        """
        return os.path.join(
            self.settings.scratchDirectory if scratch else self.settings.fileStore,
            self.uniqueID + extension
            )
        
    def getCrudeMolFilePath(self):
        "Returns the path of the crude mol file."
        return self.getFilePath('.crude.mol')

    def getRefinedMolFilePath(self):
        "Returns the path the the refined mol file."
        return self.getFilePath('.refined.mol')

    def generateRDKitGeometries(self):
        """
        Use RDKit to guess geometry.

        Save mol files of both crude and refined.
        Saves coordinates on atoms.
        """
        rdmol, rdAtIdx = self.rd_build()
        
        atoms = len(self.molecule.atoms)
        distGeomAttempts=1
        if atoms > 3:#this check prevents the number of attempts from being negative
            distGeomAttempts = 5*(atoms-3) #number of conformer attempts is just a linear scaling with molecule size, due to time considerations in practice, it is probably more like 3^(n-3) or something like that
        
        rdmol, minEid = self.rd_embed(rdmol, distGeomAttempts)
        self.saveCoordinatesFromRDMol(rdmol, minEid, rdAtIdx)
        
    def rd_build(self):
        """
        Import rmg molecule and create rdkit molecule with the same atom labeling.
        """
        return self.molecule.toRDKitMol(removeHs=False, returnMapping=True)


    def rd_embed(self, rdmol, numConfAttempts):
        """
        Embed the RDKit molecule and create the crude molecule file.
        """
        AllChem.EmbedMultipleConfs(rdmol, numConfAttempts,randomSeed=1)
        
        energy=0.0
        minEid=0;
        lowestE=9.999999e99;#start with a very high number, which would never be reached

        crude = Chem.Mol(rdmol.ToBinary())
        
        for i in range(rdmol.GetNumConformers()):
            AllChem.UFFOptimizeMolecule(rdmol,confId=i)
            energy=AllChem.UFFGetMoleculeForceField(rdmol,confId=i).CalcEnergy()
            if energy < lowestE:
                minEid = i
                lowestE = energy 
        
        with open(self.getCrudeMolFilePath(), 'w') as out3Dcrude:
            out3Dcrude.write(Chem.MolToMolBlock(crude,confId=minEid))
        
        with open(self.getRefinedMolFilePath(), 'w') as out3D:
            out3D.write(Chem.MolToMolBlock(rdmol,confId=minEid))

        return rdmol, minEid

    def saveCoordinatesFromRDMol(self, rdmol, minEid, rdAtIdx):
        # Save xyz coordinates on each atom in molecule ****
        for atom in self.molecule.atoms:
            point = rdmol.GetConformer(minEid).GetAtomPosition(atom.sortingLabel)
            atom.coords = numpy.array([point.x, point.y, point.z])
    
    def saveCoordinatesFromQMData(self, qmdata):
        """
        Save geometry info from QMData (eg CCLibData)
        """
        raise NotImplementedError

def loadThermoDataFile(filePath):
    """
    Load the specified thermo data file and return the dictionary of its contents.
    
    Returns `None` if the file is invalid or missing.
    
    Checks that the returned dictionary contains at least InChI, adjacencyList, thermoData.
    """
    if not os.path.exists(filePath):
        return None
    try:
        with open(filePath) as resultFile:
            logging.info('Reading existing thermo file {0}'.format(filePath))
            global_context = { '__builtins__': None }
            local_context = {
                '__builtins__': None,
                'True': True,
                'False': False,
                'ThermoData': rmgpy.thermo.ThermoData,
                'PointGroup': symmetry.PointGroup,
                'QMData': qmdata.QMData,
                'array': numpy.array,
                'int32': numpy.int32,
            }
            exec resultFile in global_context, local_context
    except IOError, e:
        logging.info("Couldn't read thermo file {0}".format(filePath))
        return None
    except (NameError, TypeError, SyntaxError), e:
        logging.error('The thermo file "{0}" was invalid:'.format(filePath))
        logging.exception(e)
        return None
    if not 'InChI' in local_context:
        logging.error('The thermo file "{0}" did not contain an InChI.'.format(filePath))
        return None
    if not 'adjacencyList' in local_context:
        logging.error('The thermo file "{0}" did not contain adjacencyList.'.format(filePath))
        return None
    if not 'thermoData' in local_context:
        logging.error('The thermo file "{0}" did not contain thermoData.'.format(filePath))
        return None
    return local_context

class QMMolecule:
    """ 
    A base class for QM Molecule calculations.
    
    Specific programs and methods should inherit from this and define some
    extra attributes and methods:
    
     * outputFileExtension
     * inputFileExtension
     * generateQMData() ...and whatever else is needed to make this method work.
     
    The attributes are:
    
    =================== ======================= ====================================
    Attribute           Type                    Description
    =================== ======================= ====================================
    `molecule`          :class:`Molecule`       RMG Molecule object
    `settings`          :class:`QMSettings`     Settings for QM calculations
    `uniqueID`          ``str``                 A short ID such as an augmented InChI Key
    `uniqueIDlong`      ``str``                 A long, truly unique ID such as an augmented InChI
    =================== ======================= ====================================
    
    """
    
    def __init__(self, molecule, settings):
        self.molecule = molecule
        self.settings = settings
        
        self.uniqueID = self.molecule.toAugmentedInChIKey()
        self.uniqueIDlong = self.molecule.toAugmentedInChI()
        
    def getFilePath(self, extension, scratch=True):
        """
        Returns the path to the file with the given extension.
        
        The provided extension should include the leading dot.
        If called with `scratch=False` then it will be in the `fileStore` directory,
        else `scratch=True` is assumed and it will be in the `scratchDirectory` directory.
        """
        #ToDo: this is duplicated in Geometry class. Should be refactored.
        return os.path.join(
            self.settings.scratchDirectory if scratch else self.settings.fileStore,
            self.uniqueID + extension
            )
        
    @property
    def outputFilePath(self):
        """Get the output file name."""
        return self.getFilePath(self.outputFileExtension)
    
    @property
    def inputFilePath(self):
        """Get the input file name."""
        return self.getFilePath(self.inputFileExtension)
    
    def getThermoFilePath(self):
        "Returns the path the thermo data file."
        return self.getFilePath('.thermo', scratch=False)
    
    @property
    def scriptAttempts(self):
        "The number of attempts with different script keywords"
        return len(self.keywords)
    
    @property
    def maxAttempts(self):
        "The total number of attempts to try"
        return 2 * len(self.keywords)
    
    def initialize(self):
        """
        Do any startup tasks.
        """
        self.checkReady()
    
    def checkReady(self):
        """
        Check that it's ready to run calculations.
        """
        self.settings.checkAllSet()
        self.checkPaths()
    
    def checkPaths(self):
        """
        Check the paths in the settings are OK. Make folders as necessary.
        """        
        self.settings.fileStore = os.path.expandvars(self.settings.fileStore) # to allow things like $HOME or $RMGpy
        self.settings.scratchDirectory = os.path.expandvars(self.settings.scratchDirectory)
        for path in [self.settings.fileStore, self.settings.scratchDirectory]:
            if not os.path.exists(path):
                logging.info("Creating directory %s for QM files."%os.path.abspath(path))
                os.makedirs(path)
        
    def createGeometry(self):
        """
        Creates self.geometry with RDKit geometries
        """
        self.geometry = Geometry(self.settings, self.uniqueID, self.molecule, uniqueIDlong=self.uniqueIDlong)
        self.geometry.generateRDKitGeometries()
        return self.geometry
        
    def parse(self):
        """
        Parses the results of the Mopac calculation, and returns a CCLibData object.
        """
        parser = self.getParser(self.outputFilePath)
        parser.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
        cclibData = parser.parse()
        radicalNumber = self.molecule.getRadicalCount()
        qmData = CCLibData(cclibData, radicalNumber+1) # Should `radicalNumber+1` be `self.molecule.multiplicity` in the next line of code? It's the electronic ground state degeneracy.
        return qmData
    
    def generateQMData(self):
        """
        Calculate the QM data somehow and return a CCLibData object, or None if it fails.
        """
        raise NotImplementedError("This should be defined in a subclass that inherits from QMMolecule")
        return qmdata.QMData() or None
    
    def generateThermoData(self):
        """
        Generate Thermo Data via a QM calc. 
        
        Returns None if it fails.
        """
        self.initialize()
        
        # First, see if we already have it.
        if self.loadThermoData():
            return self.thermo
        
        # If not, generate the QM data
        self.qmData = self.generateQMData()
        
        # If that fails, give up and return None.
        if self.qmData  is None:
            return None
            
        self.determinePointGroup()
        
        # If that fails, give up and return None.
        if self.pointGroup is None:
            return None
            
        self.calculateThermoData()
        Cp0 = self.molecule.calculateCp0()
        CpInf = self.molecule.calculateCpInf()
        self.thermo.Cp0 = (Cp0,"J/(mol*K)")
        self.thermo.CpInf = (CpInf,"J/(mol*K)")
        self.saveThermoData()
        return self.thermo
        
    def saveThermoData(self):
        """
        Save the generated thermo data.
        """
        self.thermo.H298.units = 'kcal/mol'
        self.thermo.S298.units = 'cal/mol/K'
        self.thermo.Cpdata.units = 'cal/mol/K'
        with open(self.getThermoFilePath(), 'w') as resultFile:
            resultFile.write('InChI = "{0!s}"\n'.format(self.uniqueIDlong))
            resultFile.write("thermoData = {0!r}\n".format(self.thermo))
            resultFile.write("pointGroup = {0!r}\n".format(self.pointGroup))
            resultFile.write("qmData = {0!r}\n".format(self.qmData))
            resultFile.write('adjacencyList = """\n{0!s}"""\n'.format(self.molecule.toAdjacencyList(removeH=False)))

    def loadThermoData(self):
        """
        Try loading a thermo data from a previous run.
        """
        filePath = self.getThermoFilePath()
        local_context = loadThermoDataFile(filePath)
        if local_context is None:
            # file does not exist or is invalid
            return None
        if local_context['InChI'] != self.uniqueIDlong:
            logging.error('The InChI in the thermo file {0} did not match the current molecule {1}'.format(filePath,self.uniqueIDlong))
            return None
        loadedMolecule = rmgpy.molecule.Molecule().fromAdjacencyList(local_context['adjacencyList'])
        if not loadedMolecule.isIsomorphic(self.molecule):
            logging.error('The adjacencyList in thermo file {0} did not match the current molecule {1}'.format(filePath,self.uniqueIDlong))
            return None
        thermo = local_context['thermoData']
        assert isinstance(thermo, rmgpy.thermo.ThermoData)
        self.thermo = thermo
        
        self.pointGroup = symmetry.pointGroupDictionary[local_context['pointGroup'].pointGroup] # point to the one in the module level dictionary with the same name
        self.qmData = local_context['qmData']
        return thermo


    def getInChiKeyAug(self):
        """
        Returns the augmented InChI from self.molecule 
        """        
        return self.molecule.toAugmentedInChIKey()

    def getMolFilePathForCalculation(self, attempt):
        """
        Get the path to the MOL file of the geometry to use for calculation `attempt`.
        
        If attempt <= self.scriptAttempts then we use the refined coordinates,
        then we start to use the crude coordinates.
        """
        if attempt <= self.scriptAttempts: # use UFF-refined coordinates
            return self.geometry.getRefinedMolFilePath()
        else:
            return self.geometry.getCrudeMolFilePath()
    
    def determinePointGroup(self):
        """
        Determine point group using the SYMMETRY Program
        
        Stores the resulting :class:`PointGroup` in self.pointGroup
        """
        assert self.qmData, "Need QM Data first in order to calculate point group."
        pgc = symmetry.PointGroupCalculator(self.settings, self.uniqueID, self.qmData)
        self.pointGroup = pgc.calculate()

    def calculateChiralityCorrection(self):
        """
        Returns the chirality correction to entropy (R*ln(2) if chiral) in J/mol/K.
        """
        if self.pointGroup.chiral:
            return rmgpy.quantity.constants.R * math.log(2)
        else:
            return 0.

    def calculateThermoData(self):
        """
        Calculate the thermodynamic properties.
        
        Stores and returns a ThermoData object as self.thermo.
        self.qmData and self.pointGroup need to be generated before this method is called.
        """
        assert self.qmData, "Need QM Data first in order to calculate thermo."
        assert self.pointGroup, "Need Point Group first in order to calculate thermo."
        
        trans = rmgpy.statmech.IdealGasTranslation( mass=self.qmData.molecularMass )
        if self.pointGroup.linear:
            # there should only be one rotational constant for a linear rotor
            rotationalConstant = rmgpy.quantity.Frequency(max(self.qmData.rotationalConstants.value), self.qmData.rotationalConstants.units)
            rot = rmgpy.statmech.LinearRotor(
                                         rotationalConstant = rotationalConstant,
                                         symmetry = self.pointGroup.symmetryNumber,
                                        )
        else:
            rot = rmgpy.statmech.NonlinearRotor(
                                         rotationalConstant = self.qmData.rotationalConstants,
                                         symmetry = self.pointGroup.symmetryNumber,
                                        )
        # @todo: should we worry about spherical top rotors?
        vib = rmgpy.statmech.HarmonicOscillator( frequencies=self.qmData.frequencies )
        
        # @todo: We need to extract or calculate E0 somehow from the qmdata
        E0 = (0, "kJ/mol")
        self.statesmodel = rmgpy.statmech.Conformer(E0=E0,
                                                    modes=[trans, rot, vib],
                                spinMultiplicity = self.qmData.groundStateDegeneracy )
        
        #we will use number of atoms from above (alternatively, we could use the chemGraph); this is needed to test whether the species is monoatomic
        # SI units are J/mol, but converted to kJ/mol for generating the thermo.
        Hf298 = self.qmData.energy.value_si / 1000
        
        S298 = self.statesmodel.getEntropy(298.0)
        Tdata = [300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0]
        Cp = [self.statesmodel.getHeatCapacity(T) for T in Tdata]
        S298 = S298 + self.calculateChiralityCorrection()
        comment = self.qmData.source or "QM calculation of some sort."
        
        thermo = ThermoData( 
                           Tdata = (Tdata,"K"),
                           Cpdata = (Cp,"J/(mol*K)"),
                           H298 = (Hf298,"kJ/mol"),
                           S298 = (S298,"J/(mol*K)"),
                           Tmin = (300.0,"K"),
                           Tmax = (2000.0,"K"),
                           comment = comment
                          )
        self.thermo = thermo
        return thermo
