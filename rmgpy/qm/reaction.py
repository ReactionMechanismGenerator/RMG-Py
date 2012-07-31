import os
import logging

from rmgpy.molecule import Molecule
from rmgpy.species import TransitionState
from molecule import QMMolecule, Geometry

import rdkit
from rdkit.Chem.Pharm3D import EmbedLib

class QMReaction:
    
    file_store_path = 'QMfiles'
    if not os.path.exists(file_store_path):
        logging.info("Creating directory %s for mol files."%os.path.abspath(file_store_path))
        os.makedirs(file_store_path)
    
    def __init__(self):
        self.geometry = None
        self.transitionState = None
        
    def fixSortLabel(molecule):
        """
        This may not be required anymore. Was needed as when molecules were created, the
        rmg sorting labels would be set after where we tried to generate the TS.
        """
        sortLbl = 0
        for atom in molecule.atoms:
            atom.sortingLabel = sortLbl
            sortLbl += 1
        return molecule
    
    def getGeometry(self, molecule):
        
        multiplicity = sum([i.radicalElectrons for i in molecule.atoms]) + 1
        geom = Geometry(molecule.toAugmentedInChIKey(), molecule, multiplicity)
        
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
        
    def write(self):
        pass
    
    def chooseMol(self):
        if len(self.reactants) == 1:
            if self.reactants[0].atoms[0].sortingLabel == -1:
                self.reactants[0] = self.fixSortLabel(self.reactants[0])
            buildTS = self.reactants[0].copy()
            actionList = self.family.forwardRecipe.actions
        else:
            if self.products[0].atoms[0].sortingLabel == -1:
                self.products[0] = self.fixSortLabel(reaction.products[0])
            buildTS = self.products[0].copy()
            actionList = self.family.reverseRecipe.actions
        
        return buildTS, actionList
        
    def generateBoundsMatrix(self, molecule):
        """
        Uses rdkit to generate the bounds matrix of a rdkit molecule.
        """
        self.geometry, multiplicity = self.getGeometry(molecule)
        rdKitMol = self.getRDKitMol(self.geometry)
        boundsMatrix = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(rdKitMol)
        
        return rdKitMol, boundsMatrix, multiplicity
        
    def editBoundsMatrix(self, molecule, boundsMatrix, actionList):
        
        for action in actionList:
            lbl1 = action[1]
            atom1 = molecule.getLabeledAtom(lbl1)
            idx1 = atom1.sortingLabel
            if len(action) ==4:
                lbl2 = action[3]
                atom2 = molecule.getLabeledAtom(lbl2)
                idx2 = atom2.sortingLabel
            if action[0].lower() == 'change_bond':
                if action[2] == '1':
                    boundsMatrix[idx1][idx2] -= 0.15
                    boundsMatrix[idx2][idx1] -= 0.15
                elif action[2] == '-1':
                    boundsMatrix[idx1][idx2] += 0.15
                    boundsMatrix[idx2][idx1] += 0.15
            elif action[0].lower() == 'form_bond':
                boundsMatrix[idx1][idx2] -= 0.25
                boundsMatrix[idx2][idx1] -= 0.25
            elif action[0].lower() == 'break_bond':
                boundsMatrix[idx1][idx2] += 0.25
                boundsMatrix[idx2][idx1] += 0.25
            elif action[0].lower() == 'gain_radical':
                # may not be as simple as the others.
                pass
            elif action[0].lower() == 'lose_radical':
                pass
        
        rdkit.DistanceGeometry.DistGeom.DoTriangleSmoothing(boundsMatrix)
        return boundsMatrix
    
    def combineBoundsMatrices(bM1, bM2):
        # creates a bounds matrix for a TS by merging the matrices for its 2 reactants
        pass