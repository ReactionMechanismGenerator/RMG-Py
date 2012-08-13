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
        
        if self.verifyOutputFile():
            logging.info("Found a successful output file already; using that.")
        else:
            for attempt in range(1, self.maxAttempts+1):
                top_keys, bottom_keys, polar_keys = method.inputFileKeys(attempt)
                self.writeInputFile(attempt, top_keys, bottom_keys, polar_keys)
                success = self.run()
                if success:
                    logging.info('Attempt {0} of {1} on species {2} succeeded.'.format(attempt, self.maxAttempts, self.molecule.toAugmentedInChI()))
                    break
            else:
                raise Exception('QM thermo calculation failed for {0}.'.format(self.molecule.toAugmentedInChI()))
        result = self.parse() # parsed in cclib
        thermo = TDPropertiesCalculator(result, self.getInChiKeyAug())

        return thermo.calculate()


        
class MopacPM3(MopacMol):
    
    "Keywords that will be added at the top of the qm input file"
    keywordsTop = {}
    keywordsTop[1] = "precise nosym"
    keywordsTop[2] = "precise nosym gnorm=0.0 nonr"
    keywordsTop[3] = "precise nosym gnorm=0.0"
    keywordsTop[4] = "precise nosym gnorm=0.0 bfgs"
    keywordsTop[5] = "precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000"
    
    "Keywords that will be added at the bottom of the qm input file"
    keywordsBottom = {}
    keywordsBottom[1] = "oldgeo thermo nosym precise "
    keywordsBottom[2] = "oldgeo thermo nosym precise "
    keywordsBottom[3] = "oldgeo thermo nosym precise "
    keywordsBottom[4] = "oldgeo thermo nosym precise "
    keywordsBottom[5] = "oldgeo thermo nosym precise "
    
    scriptAttempts = len(keywordsTop)
    maxAttempts = 2 * scriptAttempts
    
    
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