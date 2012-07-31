import os

from reaction import QMReaction
from mopac import Mopac

import rdkit
from rdkit.Chem.Pharm3D import EmbedLib

class MopacTS(QMReaction, Mopac):
    #*****change this stuff for TS
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
    keywordsBottom[1] = ""
    keywordsBottom[2] = ""
    keywordsBottom[3] = ""
    keywordsBottom[4] = ""
    
    scriptAttempts = len(keywordsTop)
    
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
        tsMol, boundsMat, multiplicity = self.generateGeometry()
        self.geometry.uniqueID = self.geometry.uniqueID + 'ts'
        self.geometry.writeMolFile(self.rdmol, self.geometry.getRefinedMolFilePath(), 0)
        method = MopacPM3(self)
        for attempt in range(1, self.scriptAttempts+1):    
            top_keys, bottom_keys, polar_keys = method.inputFileKeys(attempt, multiplicity)
            inputFileName = self.writeInputFile(attempt, top_keys, bottom_keys, polar_keys, self.scriptAttempts)
            import ipdb; ipdb.set_trace()
            success = self.run(inputFileName)
            
    
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
                    reactant = self.reactants[0]
                    reactant2 = self.reactants[1]
                except ValueError:
                    reactant = self.reactants[1]
                    reactant2 = self.reactants[0]
                    
                try:
                    self.products[0].getLabeledAtom(lblAt)
                    product = self.products[0]
                except ValueError:
                    product = self.products[1]
                
                # Merge the reactants to generate the TS template
                import ipdb; ipdb.set_trace()
                buildTS = self.reactant.merge(reactant2)
                
                # Check for rdkit molecules
                # done differently, search for rdkit mol file in QMfiles
                if reactant.rdMol == None:
                    Molecule.generate3dGeometry(reactant)
                if product.rdMol == None:
                    Molecule.generate3dGeometry(product)
                if self.rdmol == None:
                    Molecule.generate3dGeometry(buildTS)
                    
                # Check for sorting labels
                if reactant.atoms[0].sortingLabel != 0:
                    reactant = fixSortLabel(reactant)
                if product.atoms[0].sortingLabel != 0:
                    product = fixSortLabel(product)
                
                # Generate the bounds matrices for the reactant and product with the transfered atom
                boundsMat1 = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(reactant.rdMol)
                boundsMat2 = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(product.rdMol)
                
                # Get the total size of the TS bounds matrix and initialize it
                totSize = len(boundsMat1) + len(boundsMat2) - 1
                boundsMat = numpy.ones((totSize, totSize)) * 1000
                
                # Add bounds matrix 1 to corresponding place of the TS bounds matrix
                boundsMat[:len(boundsMat1),:len(boundsMat1)] = boundsMat1
                
                # Fill the bottom left of the bounds matrix with minima
                boundsMat[len(boundsMat1):, :len(boundsMat1)] = numpy.ones((len(boundsMat)-len(boundsMat1), len(boundsMat1))) * 1.07
                
                # Add bounds matrix 2, but it has to shift to the end of bounds matrix 1, and shift 
                # numbers for the reacting atom which has already been included from above
                rAtLbl = reactant.getLabeledAtom(lblAt).sortingLabel
                pAtLbl = product.getLabeledAtom(lblAt).sortingLabel
                boundsMat[len(boundsMat1):len(boundsMat1)+pAtLbl, rAtLbl] = boundsMat2[pAtLbl, :pAtLbl]
                boundsMat[rAtLbl, len(boundsMat1):len(boundsMat1)+pAtLbl] = boundsMat2[:pAtLbl, pAtLbl]
                boundsMat[rAtLbl, len(boundsMat1)+pAtLbl+1:] = boundsMat2[pAtLbl, pAtLbl+1:]
                boundsMat[len(boundsMat1)+pAtLbl+1:, rAtLbl] = boundsMat2[pAtLbl+1:, pAtLbl]
                
                # Remove all the parts of the transfered atom from the second bounds matrix
                # Incorporate the rest into the TS bounds matrix
                boundsMat2 = numpy.delete(numpy.delete(boundsMat2, pAtLbl, 1), pAtLbl, 0)
                boundsMat[-len(boundsMat2):, -len(boundsMat2):] = boundsMat2
                
        # A --> B + C or A + B --> C
        else:
            # The single species is used as the base for the transition state 
            buildTS, actionList = self.chooseMol()
            
            # Generate the RDKit::Mol from the RMG molecule and get the bounds matrix
            self.rdmol, boundsMat, multiplicity = self.generateBoundsMatrix(buildTS)
            
            # Alter the bounds matrix based on the reaction recipe
            boundsMat = self.editBoundsMatrix(buildTS, boundsMat, actionList)
            
            # Smooth the bounds matrix to speed up the optimization
            # Optimize the TS geometry in place, outputing the initial and final energies
            try:
                rdkit.Chem.Pharm3D.EmbedLib.OptimizeMol(self.rdmol, boundsMat, maxPasses = 10)
            except RuntimeError:
                pass
        
        return buildTS, boundsMat, multiplicity
        

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