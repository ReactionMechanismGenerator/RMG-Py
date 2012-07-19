import os

import logging

from mopacmol import MopacMol
from mopac import Mopac

class MopacPM3(MopacMol):
    
    def __init__(self, molecule):
        self.molecule = molecule
        
    def createInputFile(self, geometry, attempt):
        """
        Inherits the writeInputFile methods from mopac.py
        """
        multiplicity_keys = self.multiplicityKeywords[geometry.multiplicity]
        
        top_keys = "pm3 {0} {1}".format(
                multiplicity_keys,
                self.keywordsTop[attempt],
                )
        bottom_keys = "{0} pm3 {1}".format(
                self.keywordsBottom[attempt],
                multiplicity_keys,
                )
        polar_keys = "oldgeo {0} nosym precise pm3 {1}".format(
                'polar' if geometry.multiplicity == 1 else 'static',
                multiplicity_keys,
                )
        
        inputFileName = self.writeInputFile(top_keys, bottom_keys, polar_keys, attempt)
        
        return inputFileName
        
    def generateQMThermoData(self):
        # call the methods to generate the thermodata
        self.createGeometry()
        
        success = False
        for attempt in range(1, self.maxAttempts+1):
            inputFileName = self.createInputFile(self.geometry, attempt)
            success = self.run(inputFileName)
            if success:
                logging.info('Attempt {0} of {1} on species {2} succeeded.'.format(attempt, self.maxAttempts, self.molecule.toAugmentedInChI()))
                result = self.parse() # parsed in cclib
                break
        else:
            raise Exception('QM thermo calculation failed for {0}.'.format(InChIaug))
        import ipdb; ipdb.set_trace()
        # so far, result has our calculated values from cclib
        # need to be converted to thermo        
        return thermo