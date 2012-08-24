"""
authors: P Bhoorasingh, S Troiano
"""
import os
import logging
import numpy

from rmgpy.molecule import Molecule
from rmgpy.species import TransitionState
from molecule import QMMolecule, Geometry
from collections import defaultdict

import rdkit
from rdkit import DistanceGeometry
from rdkit.Chem.Pharm3D import EmbedLib

class QMReaction:
    
    file_store_path = 'QMfiles'
    if not os.path.exists(file_store_path):
        logging.info("Creating directory %s for mol files."%os.path.abspath(file_store_path))
        os.makedirs(file_store_path)
    
    def __init__(self):
        self.geometry = None
        self.transitionState = None
        
    def fixSortLabel(self, molecule):
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
    
    def getDeepCopy(self, molecule):
        """
        Gets a deep copy of the molecule so that any edits to the molecule does not affect the
        original rmg molecule.
        """
        copyMol = molecule.copy(deep=True)
        
        return copyMol
        
    def generateBoundsMatrix(self, molecule):
        """
        Uses rdkit to generate the bounds matrix of a rdkit molecule.
        """
        self.geometry, multiplicity = self.getGeometry(molecule)
        rdKitMol = self.getRDKitMol(self.geometry)
        boundsMatrix = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(rdKitMol)
        
        return rdKitMol, boundsMatrix, multiplicity
    
    def getBMParameters(self, reactant, product):
        for action in self.reaction.family.forwardRecipe.actions:
            if action[0].lower() == 'form_bond' or action[0].lower() == 'break_bond':
                atlbl1 = action[1]
                atlbl2 = action[3]
        
        if reactant.isCyclic():
            rdMol, bMatrix, mult = self.generateBoundsMatrix(reactant)
            lbl = (reactant.getLabeledAtom(atlbl1).sortingLabel, reactant.getLabeledAtom(atlbl2).sortingLabel)
            oRDMol, oBM, oMult = self.generateBoundsMatrix(product)
            other = (product.getLabeledAtom(atlbl1).sortingLabel, product.getLabeledAtom(atlbl2).sortingLabel)
            other = (oBM[other[0]][other[1]], oBM[other[1]][other[0]])
        elif product.isCyclic():
            rdMol, bMatrix, mult = self.generateBoundsMatrix(product)
            lbl = (product.getLabeledAtom(atlbl1).sortingLabel, product.getLabeledAtom(atlbl2).sortingLabel)
            oRDMol, oBM, oMult = self.generateBoundsMatrix(reactant)
            other = (reactant.getLabeledAtom(atlbl1).sortingLabel, reactant.getLabeledAtom(atlbl2).sortingLabel)
            other = (oBM[other[0]][other[1]], oBM[other[1]][other[0]])
        
        return rdMol, bMatrix, mult, lbl, other
    
            
        
    
    def getTSBMatrix(self):
        
        reactant = self.getDeepCopy(self.reactants[0])
        product = self.getDeepCopy(self.products[0])
        if reactant.atoms[0].sortingLabel == -1:
            reactant = self.fixSortLabel(reactant)
        if product.atoms[0].sortingLabel == -1:
            product = self.fixSortLabel(product)
            
        rdMol, tsBM, mult, lbl, other = self.getBMParameters(reactant, product)
        
        DistanceGeometry.DoTriangleSmoothing(tsBM)
        
        return rdMol, tsBM, mult, lbl, other