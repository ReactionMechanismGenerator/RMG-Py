import os
import logging
import numpy

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
        
        # Smooth the bounds matrix to speed up the optimization
        rdkit.DistanceGeometry.DistGeom.DoTriangleSmoothing(boundsMatrix)
        
        return boundsMatrix
    
    def combineBoundsMatrices(self, boundsMat1, boundsMat2, lblAt1, lblAt2):
        # creates a bounds matrix for a TS by merging the matrices for its 2 reactants
        # Add bounds matrix 1 to corresponding place of the TS bounds matrix
        totSize = len(boundsMat1) + len(boundsMat2) - 1
        boundsMat = numpy.ones((totSize, totSize)) * 1000
        
        boundsMat[:len(boundsMat1),:len(boundsMat1)] = boundsMat1
        
        # Fill the bottom left of the bounds matrix with minima
        boundsMat[len(boundsMat1):, :len(boundsMat1)] = numpy.ones((len(boundsMat)-len(boundsMat1), len(boundsMat1))) * 1.07
        
        # Add bounds matrix 2, but it has to shift to the end of bounds matrix 1, and shift 
        # numbers for the reacting atom which has already been included from above
        boundsMat[len(boundsMat1):len(boundsMat1)+lblAt2, lblAt1] = boundsMat2[lblAt2, :lblAt2]
        boundsMat[lblAt1, len(boundsMat1):len(boundsMat1)+lblAt2] = boundsMat2[:lblAt2, lblAt2]
        boundsMat[lblAt1, len(boundsMat1)+lblAt2+1:] = boundsMat2[lblAt2, lblAt2+1:]
        boundsMat[len(boundsMat1)+lblAt2+1:, lblAt1] = boundsMat2[lblAt2+1:, lblAt2]
        
        # Remove all the parts of the transfered atom from the second bounds matrix
        # Incorporate the rest into the TS bounds matrix
        boundsMat2 = numpy.delete(numpy.delete(boundsMat2, lblAt2, 1), lblAt2, 0)
        boundsMat[-len(boundsMat2):, -len(boundsMat2):] = boundsMat2
        
        return boundsMat
        
    def generateGeometry(self):
        # A --> B or A + B --> C + D
        if len(self.reactants) == len(self.products):
            # 1 reactant
            if len(self.reactants) == 1:
                pass
            
            # 2 reactants
            else:
                actionList = self.family.forwardRecipe.actions
                for action in actionList:
                    if action[0].lower() == 'form_bond':
                        lbl1 = action[1]
                        lbl2 = action[3]
                    elif action[0].lower() == 'break_bond':
                        lbl3 = action[1]
                        lbl4 = action[3]
                
                # Find the atom being transferred in the reaction
                if lbl1 == lbl3 or lbl1 == lbl4:
                    lblAt = lbl1
                else:
                    lblAt = lbl2
                
                # Derive the bounds matrix from the reactants and products
                try:
                    self.reactants[0].getLabeledAtom(lblAt)
                    reactant = self.getDeepCopy(self.reactants[0])
                    reactant2 = self.getDeepCopy(self.reactants[1])
                except ValueError:
                    reactant = self.getDeepCopy(self.reactants[1])
                    reactant2 = self.getDeepCopy(self.reactants[0])
                    
                try:
                    self.products[0].getLabeledAtom(lblAt)
                    product = self.getDeepCopy(self.products[0])
                except ValueError:
                    product = self.getDeepCopy(self.products[1])
                
                # Merge the reactants to generate the TS template
                buildTS = reactant.merge(reactant2)
                
                # Check for sorting labels
                if reactant.atoms[0].sortingLabel != 0:
                    reactant = self.fixSortLabel(reactant)
                if product.atoms[0].sortingLabel != 0:
                    product = self.fixSortLabel(product)
                
                # Generate the bounds matrices for the reactant and product with the transfered atom
                reactant.rdmol, boundsMatR, multiplicityR = self.generateBoundsMatrix(reactant)
                product.rdmol, boundsMatP, multiplicityP = self.generateBoundsMatrix(product)
                
                # Calculate the multiplicity
                multiplicity = sum([i.radicalElectrons for i in buildTS.atoms]) +1
                
                rAtLbl = reactant.getLabeledAtom(lblAt).sortingLabel
                pAtLbl = product.getLabeledAtom(lblAt).sortingLabel
                
                # Get the total size of the TS bounds matrix and initialize it
                boundsMat = self.combineBoundsMatrices(boundsMatR, boundsMatP, rAtLbl, pAtLbl)
                
                self.fixSortLabel(buildTS)
                
                boundsMat = self.editBoundsMatrix(buildTS, boundsMat, actionList)
                
                try:
                    rdkit.Chem.Pharm3D.EmbedLib.OptimizeMol(self.rdmol, boundsMat, maxPasses = 10)
                except RuntimeError:
                    pass
                
                multiplicity = sum([i.radicalElectrons for i in buildTS.atoms]) +1
                
        # A --> B + C or A + B --> C
        else:
            # The single species is used as the base for the transition state 
            buildTS, actionList = self.chooseMol()
            
            # Generate the RDKit::Mol from the RMG molecule and get the bounds matrix
            self.rdmol, boundsMat, multiplicity = self.generateBoundsMatrix(buildTS)
            
            # Alter the bounds matrix based on the reaction recipe
            boundsMat = self.editBoundsMatrix(buildTS, boundsMat, actionList)
            
            # Optimize the TS geometry in place, outputing the initial and final energies
            try:
                rdkit.Chem.Pharm3D.EmbedLib.OptimizeMol(self.rdmol, boundsMat, maxPasses = 10)
            except RuntimeError:
                pass
        
        return buildTS, boundsMat, multiplicity