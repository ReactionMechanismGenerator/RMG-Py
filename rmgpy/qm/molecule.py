import os
import logging
import math

import numpy
from rdkit import Chem
from rdkit.Chem import AllChem, Pharm3D

import rmgpy.quantity
from rmgpy.thermo import ThermoData
import rmgpy.statmech
import symmetry
import qmdata

class RDKitFailedError(Exception):
    """For when RDkit failed. try the next reaction """
    pass
    
class Geometry:
    """
    A geometry, used for quantum calculations.
    
    Created from a molecule. Geometry estimated by RDKit.
    """
    def __init__(self, settings, uniqueID, molecule, multiplicity, uniqueIDlong=None):
        self.settings = settings
        #: A short unique ID such as an augmented InChI Key.
        self.uniqueID = uniqueID
        self.molecule = molecule
        #: The multiplicity, eg. the number of free radicals plus one.
        self.multiplicity = multiplicity
        if uniqueIDlong is None:
            self.uniqueIDlong = uniqueID
        else:
            #: Long, truly unique, ID, such as the augmented InChI.
            self.uniqueIDlong = uniqueIDlong

    def getFilePath(self, extension):
        """
        Returns the path to the file with the given extension.
        
        The provided extension should include the leading dot.
        """
        return os.path.join(self.settings.fileStore, self.uniqueID  + extension)
        
    def getCrudeMolFilePath(self):
        "Returns the path of the crude mol file."
        return self.getFilePath('.crude.mol')

    def getRefinedMolFilePath(self):
        "Returns the path the the refined mol file."
        return self.getFilePath('.refined.mol')

    def generateRDKitGeometries(self, boundsMatrix=None, atomMatch=None):
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
        
        rdmol, minEid = self.rd_embed(rdmol, distGeomAttempts, bm=boundsMatrix, match=atomMatch)
        self.save_coordinates(rdmol, minEid, rdAtIdx)
        
    def rd_build(self):
        """
        Import rmg molecule and create rdkit molecule with the same atom labeling.
        """
        rdAtomIdx = {} # dictionary of rdkit atom indices
        # Initialize a blank Editable molecule and add all the atoms from RMG molecule
        rdmol = AllChem.rdchem.EditableMol(AllChem.rdchem.Mol())
        for index, atom in enumerate(self.molecule.vertices):
            rdAtom = AllChem.rdchem.Atom(atom.element.symbol)
            rdAtom.SetNumRadicalElectrons(atom.radicalElectrons)
            rdmol.AddAtom(rdAtom)
            rdAtomIdx[atom] = index

        # Add the bonds
        for atom1 in self.molecule.vertices:
            for atom2, bond in atom1.edges.items():
                index1 = rdAtomIdx[atom1]
                index2 = rdAtomIdx[atom2]
                if index1 > index2:
                    # Check the RMG bond order and add the appropriate rdkit bond.
                    if bond.order == 'S':
                        rdBond = AllChem.rdchem.BondType.SINGLE
                    elif bond.order == 'D':
                        rdBond = AllChem.rdchem.BondType.DOUBLE
                    elif bond.order == 'T':
                        rdBond = AllChem.rdchem.BondType.TRIPLE
                    elif bond.order == 'B':
                        rdBond = AllChem.rdchem.BondType.AROMATIC
                    else:
                        logging.error('Unknown bond order')
                    rdmol.AddBond(index1, index2, rdBond)

        # Make editable mol into a mol and rectify the molecule
        rdmol = rdmol.GetMol()
        Chem.SanitizeMol(rdmol)

        return rdmol, rdAtomIdx

    def rd_embed(self, rdmol, numConfAttempts, bm=None, match=None):
        """
        Embed the RDKit molecule and create the crude molecule file.
        """
        if bm == None:
            AllChem.EmbedMultipleConfs(rdmol, numConfAttempts,randomSeed=1)
            crude = Chem.Mol(rdmol.ToBinary())
            rdmol, minEid = self.optimize(rdmol)
        else:
            for i in range(0,numConfAttempts):
                """
                Embed the molecule according to the bounds matrix. Built to handle possible failures
                of some of the embedding attempts.
                """
                try:
                    Pharm3D.EmbedLib.EmbedMol(rdmol, bm, atomMatch=match)
                except ValueError:
                    pass
                except RuntimeError:
                    raise RDKitFailedError()
            """
            RDKit currently embeds the conformers and sets the id as 0, so even though multiple
            conformers have been generated, only 1 can be called. Below the id's are resolved.
            """
            for i in range(len(rdmol.GetConformers())):
                rdmol.GetConformers()[i].SetId(i)
            crude = Chem.Mol(rdmol.ToBinary())
            rdmol, minEid = self.optimize(rdmol, boundsMatrix=bm, atomMatch=match)
        
        self.writeMolFile(crude, self.getCrudeMolFilePath(), minEid)
        self.writeMolFile(rdmol, self.getRefinedMolFilePath(), minEid)
        
        return rdmol, minEid
    
    def optimize(self, rdmol, boundsMatrix=None, atomMatch=None):
        
        energy=0.0
        minEid=0;
        lowestE=9.999999e99;#start with a very high number, which would never be reached
        
        for i in range(rdmol.GetNumConformers()):
            if boundsMatrix == None:    
                AllChem.UFFOptimizeMolecule(rdmol,confId=i)
                energy=AllChem.UFFGetMoleculeForceField(rdmol,confId=i).CalcEnergy()
            else:
                eBefore, energy = Pharm3D.EmbedLib.OptimizeMol(rdmol, boundsMatrix)
            
            if energy < lowestE:
                minEid = i
                lowestE = energy
                
        return rdmol, minEid
        
    def writeMolFile(self, mol, path, minEid):
        with open(path, 'w') as out3Dcrude:
            out3Dcrude.write(Chem.MolToMolBlock(mol,confId=minEid))
    
    def save_coordinates(self, rdmol, minEid, rdAtIdx):
        # Save xyz coordinates on each atom in molecule ****
        for atom in self.molecule.atoms:
            point = rdmol.GetConformer(minEid).GetAtomPosition(atom.sortingLabel)
            atom.coords = [point.x, point.y, point.z]

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
    """
    
    def __init__(self, molecule, settings):
        self.molecule = molecule
        self.settings = settings
        
        self.uniqueID = self.molecule.toAugmentedInChIKey()
        self.uniqueIDlong = self.molecule.toAugmentedInChI()
        

    def getFilePath(self, extension):
        """
        Returns the path to the file with the given extension.
        
        The provided extension should include the leading dot.
        """
        return os.path.join(self.settings.fileStore, self.uniqueID  + extension)
        
    @property
    def outputFilePath(self):
        """Get the output file name."""
        return self.getFilePath(self.outputFileExtension)
    
    @property
    def inputFilePath(self):
        """Get the input file name."""
        return self.getFilePath(self.inputFileExtension)
        
    def createGeometry(self, boundsMatrix=None, atomMatch=None):
        """
        Creates self.geometry with RDKit geometries
        """
        multiplicity = sum([i.radicalElectrons for i in self.molecule.atoms]) + 1
        self.geometry = Geometry(self.settings, self.uniqueID, self.molecule, multiplicity, uniqueIDlong=self.uniqueIDlong)
        self.geometry.generateRDKitGeometries(boundsMatrix, atomMatch)
        return self.geometry
    
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
        
        # First, see if we already have it.
        if self.loadThermoData():
            return self.thermo
        
        # If not, generate the QM data
        self.qmData = self.generateQMData()
        
        # If that fails, give up and return None.
        if self.qmData  is None:
            return None
            
        self.determinePointGroup()
        self.calculateThermoData()
        self.saveThermoData()
        return self.thermo
        
    def saveThermoData(self):
        """
        Save the generated thermo data.
        """
        self.thermo.H298.units = 'kcal/mol'
        self.thermo.S298.units = 'cal/mol/K'
        self.thermo.Cpdata.units = 'cal/mol/K'
        with open(self.getFilePath('.thermo'), 'w') as resultFile:
            resultFile.write('InChI = "{0!s}"\n'.format(self.uniqueIDlong))
            resultFile.write("thermoData = {0!r}\n".format(self.thermo))
            resultFile.write("pointGroup = {0!r}\n".format(self.pointGroup))
            resultFile.write("qmData = {0!r}\n".format(self.qmData))
            resultFile.write('adjacencyList = """\n{0!s}"""\n'.format(self.molecule.toAdjacencyList(removeH=True)))

    def loadThermoData(self):
        """
        Try loading a thermo data from a previous run.
        """
        filePath = self.getFilePath('.thermo')
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
        
        # self.qmData.rotationalConstants is in Hz.
        h_over_8_pi_squared = 8.39201316e-36 # kg m^2 / s
        inertias = ([h_over_8_pi_squared / frequency for frequency in self.qmData.rotationalConstants], "kg*m^2")
        
        trans = rmgpy.statmech.Translation( mass=(self.qmData.molecularMass,"g/mol") )
        rot = rmgpy.statmech.RigidRotor( linear = self.pointGroup.linear,
                                         inertia = inertias,
                                         symmetry = self.pointGroup.symmetryNumber,
                                        )
        vib = rmgpy.statmech.HarmonicOscillator( frequencies=self.qmData.frequencies )
        self.statesmodel = rmgpy.statmech.StatesModel( modes=[trans, rot, vib],
                                spinMultiplicity = self.qmData.groundStateDegeneracy )
        
        #we will use number of atoms from above (alternatively, we could use the chemGraph); this is needed to test whether the species is monoatomic
        # Hartree_to_kcal = 627.5095 # from Gaussian thermo white paper
        Hartree_to_kJmol = 2625.49962 # from Wikipedia, which cites CODATA2010
        Hf298 = self.qmData.energy * Hartree_to_kJmol
        
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
