import os
import logging
import external.cclib.parser
import openbabel
import time
from subprocess import Popen
from copy import deepcopy
import numpy

from rmgpy.molecule import Molecule
from rmgpy.species import TransitionState
from molecule import QMMolecule, Geometry
from collections import defaultdict
from rmgpy.data.base import Entry
from rmgpy.data.kinetics import KineticsFamily, ReactionRecipe, saveEntry
from rmgpy.data.kinetics.transitionstates import TransitionStates, DistanceData

try:
    import rdkit
    from rdkit import DistanceGeometry
    from rdkit.Chem.Pharm3D import EmbedLib
except ImportError:
    logging.info("To use transition state searches, you must correctly install rdkit")

transitionStates = TransitionStates()
transitionStates.load(os.path.join(os.getenv('HOME'), 'Code/RMG-database/input/kinetics/families/H_Abstraction'), None, None)

class QMReaction:
    
    file_store_path = 'QMfiles'
    if not os.path.exists(file_store_path):
        logging.info("Creating directory %s for mol files."%os.path.abspath(file_store_path))
        os.makedirs(file_store_path)
    
    def __init__(self, reaction, settings):
        self.reaction = reaction
        self.settings = settings
        
        reactants = sorted([s.toSMILES() for s in self.reaction.reactants])
        products = sorted([s.toSMILES() for s in self.reaction.products])
        stringID = "+".join(reactants) + "_" + "+".join(products)
        
        self.uniqueID = stringID
        
        self.geometry = None
        self.transitionState = None
    
    def getFilePath(self, extension):
        """
        Should return the path to the file with the given extension.
        
        The provided extension should include the leading dot.
        
        Need to define some reaction line notation.
        Possibly '<Reaction_Family>/<reactant1SMILES>+<reactant2SMILES>--<product1SMILES>+<product2SMILES>' ???
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
        
    def fixSortLabel(self, molecule):
        """
        This may not be required anymore. Was needed as when molecules were created, the
        rmg sorting labels would be set after where we tried to generate the TS.
        """
        sortLbl = 0
        for vertex in molecule.vertices:
            vertex.sortingLabel = sortLbl
            sortLbl += 1
        return molecule
    
    def getGeometry(self, molecule, settings):
        
        multiplicity = sum([i.radicalElectrons for i in molecule.atoms]) + 1
        geom = Geometry(settings, molecule.toAugmentedInChIKey(), molecule, multiplicity)
        
        return geom, multiplicity
        
    def getRDKitMol(self, geometry):
        """
        Check there is no RDKit mol file already made. If so, use rdkit to make a rdmol from
        a mol file. If not, make rdmol from geometry.
        """ 
        if not os.path.exists(geometry.getCrudeMolFilePath()):
            geometry.generateRDKitGeometries()
        rdKitMol = rdkit.Chem.MolFromMolFile(geometry.getCrudeMolFilePath(), removeHs=False)      
        
        return rdKitMol
        
    def generateBoundsMatrix(self, molecule):
        """
        Uses rdkit to generate the bounds matrix of a rdkit molecule.
        """
        self.geometry, multiplicity = self.getGeometry(molecule, self.settings)
        rdKitMol = self.getRDKitMol(self.geometry)
        boundsMatrix = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(rdKitMol)
        
        return rdKitMol, boundsMatrix, multiplicity
    
    def setLimits(self, bm, lbl1, lbl2, num, diff):
        if lbl1 > lbl2:
            bm[lbl2][lbl1] = num + diff
            bm[lbl1][lbl2] = num
        else:
            bm[lbl2][lbl1] = num
            bm[lbl1][lbl2] = num + diff
    
        return bm
    
    def editMatrix(self, reactant, bm):
        
        # edit bounds distances to align reacting atoms
        if self.reaction.label.lower() == 'h_abstraction':
            
            """
            Reduce the minimum distance between atoms on the two reactants
            to 1.8 Angstrom. 
            """
            sect = len(reactant.split()[1].atoms)        
            bm[sect:,:sect] = 1.8
        
            lbl1 = reactant.getLabeledAtom('*1').sortingLabel
            lbl2 = reactant.getLabeledAtom('*2').sortingLabel
            lbl3 = reactant.getLabeledAtom('*3').sortingLabel
            
            labels = [lbl1, lbl2, lbl3]
            atomMatch = ((lbl1,),(lbl2,),(lbl3,))
            
            distanceData = transitionStates.estimateDistances(self.reaction)
            
            bm = self.setLimits(bm, lbl1, lbl2, distanceData.distances['d12'], 0.1)
            bm = self.setLimits(bm, lbl2, lbl3, distanceData.distances['d23'], 0.001)
            bm = self.setLimits(bm, lbl1, lbl3, distanceData.distances['d13'], 0.001)
        # elif .....:
        
        return bm, labels, atomMatch
        
    def generateGeometry(self):
        """
        
        """
        if len(self.reaction.reactants)==2:
            reactant = self.reaction.reactants[0].merge(self.reaction.reactants[1])
        if len(self.reaction.products)==2:
            product = self.reaction.products[0].merge(self.reaction.products[1])
        
        reactant = self.fixSortLabel(reactant)
        product = self.fixSortLabel(product)
        
        tsRDMol, tsBM, tsMult = self.generateBoundsMatrix(reactant)
        tsBM, labels, atomMatch = self.editMatrix(reactant, tsBM)
        atoms = len(reactant.atoms)
        distGeomAttempts = 5*(atoms-3) # number of conformers embedded from the bounds matrix
        
        setBM = rdkit.DistanceGeometry.DoTriangleSmoothing(tsBM)
        if setBM:
            self.geometry.rd_embed(tsRDMol, distGeomAttempts, bm=tsBM, match=atomMatch)
            
            self.writeInputFile(1)
            import ipdb; ipdb.set_trace()
            self.run()
            self.writeIRCFile()
            self.run()
            self.verifyTSGeometry()
            