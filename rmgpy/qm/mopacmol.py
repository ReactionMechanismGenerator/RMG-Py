import os

import logging

from molecule import QMMolecule
from mopac import Mopac

class MopacMol(QMMolecule, Mopac):
    def __init__(self, molecule):
        self.molecule = molecule
    
    def generateQMThermoData(self):
        # call the methods to generate the thermodata
        self.createGeometry()
        
        success = False
        for attempt in range(1, self.maxAttempts+1):
            inputFileName = self.writeInputFile(attempt)
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