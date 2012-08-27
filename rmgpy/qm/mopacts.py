import os

from reaction import QMReaction
from mopac import Mopac

import rdkit
from rdkit.Chem.Pharm3D import EmbedLib
import logging

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