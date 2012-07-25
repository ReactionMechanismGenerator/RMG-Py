import os

import logging

from molecule import QMMolecule, TDPropertiesCalculator
from mopac import Mopac

class MopacMol(QMMolecule, Mopac):
    def generateQMThermoData(self):
        # call the methods to generate the thermodata
        self.createGeometry()
        
        success = False
        method = MopacPM3(self)
        for attempt in range(1, self.maxAttempts+1):
            top_keys, bottom_keys, polar_keys = method.inputFileKeys(attempt)
            inputFileName = self.writeInputFile(attempt, top_keys, bottom_keys, polar_keys)
            success = self.run(inputFileName)
            if success:
                logging.info('Attempt {0} of {1} on species {2} succeeded.'.format(attempt, self.maxAttempts, self.molecule.toAugmentedInChI()))
                result = self.parse() # parsed in cclib
                break
        else:
            raise Exception('QM thermo calculation failed for {0}.'.format(InChIaug))
        
        thermo = TDPropertiesCalculator(result, self.getInChiAug())

        return thermo.calculate()


        
class MopacPM3(MopacMol):
    def inputFileKeys(self, attempt):
        """
        Inherits the writeInputFile methods from mopac.py
        """
        multiplicity_keys = self.multiplicityKeywords[self.molecule.geometry.multiplicity]

        top_keys = "pm3 {0} {1}".format(
                multiplicity_keys,
                self.keywordsTop[attempt],
                )
        bottom_keys = "{0} pm3 {1}".format(
                self.keywordsBottom[attempt],
                multiplicity_keys,
                )
        polar_keys = "oldgeo {0} nosym precise pm3 {1}".format(
                'polar' if self.molecule.geometry.multiplicity == 1 else 'static',
                multiplicity_keys,
                )

        return top_keys, bottom_keys, polar_keys