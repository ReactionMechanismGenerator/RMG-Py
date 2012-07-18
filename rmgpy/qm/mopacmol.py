import os

from molecule import QMMolecule
from mopac import Mopac

class MopacMol(QMMolecule, Mopac):
    def generateQMThermoData(molecule):
        # call the methods to generate the thermodynamic data
        thermo = Mopac.generateThermoData()
        return thermo