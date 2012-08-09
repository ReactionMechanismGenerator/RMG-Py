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
    keywordsTop[1] = "ts vectors xyz"
    keywordsTop[2] = "ts vectors xyz recalc=5"
    keywordsTop[3] = "ts vectors xyz ddmin=0.0001"
    keywordsTop[4] = "ts vectors xyz recalc=5 ddmin=0.0001"
    
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
        reactants, products = self.twoEnded()
        self.geometry.uniqueID = self.geometry.uniqueID + 'ts'
        self.geometry.writeMolFile(self.rdmol, self.geometry.getRefinedMolFilePath(), 0)
        method = MopacPM3(self)
        for attempt in range(1, self.scriptAttempts+1):    
            top_keys, bottom_keys, polar_keys = method.inputFileKeys(attempt, multiplicity)
            inputFileName = self.writeInputFile(attempt, top_keys, bottom_keys, polar_keys, self.scriptAttempts)
            success = self.run(inputFileName)
    
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