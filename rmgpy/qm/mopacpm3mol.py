import os

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
        
        attemptNumber = 1
        success = False
        while not success and (attemptNumber <= self.maxAttempts):
            inputFileName = self.createInputFile(self.geometry, attemptNumber)
            success = self.run(inputFileName)
            if success:
                logging.info('Attempt {0} on species {1} succeeded.'.format(attemptNumber, InChIaug))
                """
                TODO Rotor Scan not yet implemented here.
                """
            else:
                if attemptNumber == maxAttemptNumber:
                    logging.info('Last attempt on species {0} failed.'.format(InChIaug))
        for attempt in range(method.max_attempts):
            writer.writeInputFile(self.geometry, attempt)
            success = method.runJob()
            method.parseResult()
            if success: break
        else:
            raise Exception("Couldn't generate thermo data.")
        thermoData = method.processResult()
        thermo = MopacMol.generateQMThermoData(molecule)
        
        return thermo