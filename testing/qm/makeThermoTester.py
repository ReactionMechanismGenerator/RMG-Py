"""
A script to generate an input file to the RMG-Java QM Data estimator,
to aid comparison.
"""

from rmgpy.molecule import Molecule
from rmgpy.species import *
from rmgpy.thermo import ThermoData
from rmgpy.statmech import *
from rmgpy.quantity import constants

import rmgpy.qm
import rmgpy.qm.main
import re

folder = 'QMfiles'
settings = rmgpy.qm.main.QMSettings()
settings.software = 'mopac'
settings.fileStore = folder
settings.scratchDirectory = folder
settings.onlyCyclics = False

import os


with open(os.path.expandvars('$RMG/thermotest/input.txt'),'w') as output:
    output.write(
"""
Database: RMG_database

//QM? true/false
true
//method: both/gaussian03/mopac/mm4/mm4hr
mopac
//ForCyclicsOnly? true/false
false
//maxradnumforQM?
0
//CheckConnectivity? off/check/confirm
check

MaxCarbonNumberPerSpecies:     20
MaxOxygenNumberPerSpecies:     20
MaxRadicalNumberPerSpecies:    20
MaxSulfurNumberPerSpecies:     20
MaxSiliconNumberPerSpecies:    20
MaxHeavyAtomNumberPerSpecies: 100
MaxCycleNumberPerSpecies:      20
END

PrimaryThermoLibrary:
Name: RMG-minimal
Location: primaryThermoLibrary
END

""")


    for f in os.listdir(folder):
        f = os.path.join(folder,f)
        stem, ext = os.path.splitext(f)
        if ext != '.thermo':
            continue
        
        data = rmgpy.qm.molecule.loadThermoDataFile(f)
        mol = Molecule().fromAdjacencyList(data['adjacencyList'])
    
        output.write('// {0!s}\n'.format(data['InChI']))
        output.write('// {0!r}\n'.format(data['thermoData']).replace('),','),\n//           '))
        output.write(mol.toSMILES())
        output.write(data['adjacencyList']+'\n')
            