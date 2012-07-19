import os

import logging

from mopacmol import MopacMol

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