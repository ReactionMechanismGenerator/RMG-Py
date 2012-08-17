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

mols = []

for f in os.listdir(folder):
    f = os.path.join(folder,f)
    stem, ext = os.path.splitext(f)
    if ext != '.thermo':
        continue
    
    with open(f) as thermofile:
        line = thermofile.readline()
        match = re.match('InChI = "(.*?)(\/mult\d)*"',line)
        if not match: continue
        inchi = match.group(1)
        mult = match.group(2)
        if mult: continue
        mol = Molecule()
        mol.fromInChI(inchi)
        mmol = rmgpy.qm.mopac.MopacMolPM3(mol, settings)
        print mmol.uniqueID, f
        thermo = mmol.loadThermoData()
        if not thermo: continue
        print "Trying ",mol
        mols.append(mmol)


with open('/Users/rwest/Code/RMG-Java/thermotest/input.txt','w') as f:
    f.write(
"""
Database: RMG_database

//QM? true/false
true
//method: both/gaussian03/mopac/mm4/mm4hr 
mopac
//ForCyclicsOnly? true/false
false
//maxradnumforQM?
3
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
    for mmol in mols:
        f.write(mmol.molecule.toSMILES()+'\n')
        f.write(mmol.molecule.toAdjacencyList()+'\n')
        