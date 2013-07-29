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
    
    def __init__(self, reaction, settings):
        self.reaction = reaction
        self.settings = settings
        
        self.geometry = None
        self.transitionState = None
        
    def getFilePath(self, extension):
        """
        Should return the path to the file with the given extension.
        
        The provided extension should include the leading dot.
        
        Need to define some reaction line notation.
        Possibly '<Reaction_Family>/<reactant1SMILES>+<reactant2SMILES>--<product1SMILES>+<product2SMILES>' ???
        """
        raise NotImplementedError("This should be a unique string representing the reaction")
        return os.path.join(self.settings.fileStore, self.uniqueID  + extension)
        
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
    
    def getDeepCopy(self, molecule):
        """
        Gets a deep copy of the molecule so that any edits to the molecule does not affect the
        original rmg molecule.
        """
        copyMol = molecule.copy(deep=True)
        
        return copyMol
        
    def generateBoundsMatrix(self, molecule, settings):
        """
        Uses rdkit to generate the bounds matrix of a rdkit molecule.
        """
        self.geometry, multiplicity = self.getGeometry(molecule, settings)
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
    
    def insertString(self, string, insert):
        num = len(string)
        for i in range(1, len(string)):
            j = num - i
            string.insert(j, insert)
        return string
    
    def stretchBond(self, editLine, line):
        bondLen = float(editLine.pop(1))
        bondLen += 0.1
        bondLen = str(round(bondLen, 6))
        editLine.insert(1, bondLen)
        difSplit = line.split('    ')
        oSplit = difSplit[1].split()
        oSplit.pop(0)
        oSplit.insert(0, bondLen + '  ')
        difSplit[1] = ''.join(oSplit)
        
        return difSplit
    
    def editAngle(self, editLine, line):
        angle = float(editLine.pop(3))
        angle -= 10
        angle = str(round(angle, 6))
        editLine.insert(3, angle)
        difSplit = line.split('    ')
        oSplit = difSplit[1].split()
        oSplit.pop(0)
        oSplit.insert(0, bondLen + '  ')
        difSplit[1] = ''.join(oSplit)
        
        return difSplit
        
    def nonEditableVar(self, splitString):
        for i in range(1, 4):
            splitString[i] = splitString[i][:-1] + '0'
            
        return splitString
    
    def fixBond(self, inputString, lbl):
        splits = inputString.splitlines()
        if lbl[0] > lbl[1]:
            line = splits[lbl[0] + 3]
            otherLine = splits[lbl[1] + 3]
            editLine = line.split()
            if line.split()[-3] == str(lbl[1] + 1):
                difSplit = self.stretchBond(editLine, line)
            elif line.split()[-2] == str(lbl[1] + 1):
                disSplit = self.editAngle(editLine, line)
                import ipdb; ipdb.set_trace()
            otherSplit = otherLine.split('    ')
            difSplit = self.nonEditableVar(difSplit)
            otherSplit = self.nonEditableVar(otherSplit)
            
            difSplit = self.insertString(difSplit, '    ')
            otherSplit = self.insertString(otherSplit, '    ')
            
            splits[lbl[0] + 3] = ''.join(difSplit)
            splits[lbl[1] + 3] = ''.join(otherSplit)
            splits = self.insertString(splits, '\n')
            
            inputString = ''.join(splits) + '\n'
        elif lbl[1] > lbl[0]:
            line = splits[lbl[1] + 3]
            otherLine = splits[lbl[0] + 3]
            editLine = line.split()
            if line.split()[-3] == str(lbl[0] + 1):
                difSplit = self.stretchBond(editLine, line)
            elif line.split()[-2] == str(lbl[0] + 1):
                import ipdb; ipdb.set_trace()
                disSplit = self.editAngle(editLine, line)
            otherSplit = otherLine.split('    ')
            difSplit = self.nonEditableVar(difSplit)
            otherSplit = self.nonEditableVar(otherSplit)
            
            difSplit = self.insertString(difSplit, '    ')
            otherSplit = self.insertString(otherSplit, '    ')
            
            splits[lbl[1] + 3] = ''.join(difSplit)
            splits[lbl[0] + 3] = ''.join(otherSplit)
            splits = self.insertString(splits, '\n')
            
            inputString = ''.join(splits) + '\n'
        
        return inputString
    
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