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

inchidict = {}
thermodict = {}
keydict = {}

with open('/Users/rwest/Code/RMG-Java/thermotest/RMG.log') as f:
    for line in f:
        match = re.match('Attempt #(\d+) on species ([A-Z]{14})-([A-Z]{10}) \(InChI=1/(.*?)\) succeeded.',line)
        if match:
            attempt = match.group(1)
            keystart = match.group(2)
            keyend = match.group(3)
            inchi = match.group(4)
            inchidict[keystart] = inchi
            keydict[inchi] = keystart
        match = re.match('Thermo for ([A-Z]{14})-([A-Z]{10}):\s+(.*)',line)
        if match:
            keystart = match.group(1)
            keyend = match.group(2)
            thermo = match.group(3).split()
            thermo = [float(t) for t in thermo]
            
            td = ThermoData( 
             Tdata = ([300,400,500,600,800,1000,1500],"K"),
             Cpdata = (thermo[2:9],"cal/(mol*K)"),
             H298 = (thermo[0],"kcal/mol"),
             S298 = (thermo[1],"cal/(mol*K)"),
             Tmin = (300.0,"K"),
             Tmax = (2000.0,"K"),
             comment = "from RMG-Java TDE"
            )
            thermodict[keystart] = td

            
                    
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
        
        assert inchi.startswith('InChI=1S/')
        inchimain = inchi[9:]
        keystart = keydict[inchimain]
        jthermo = thermodict[keystart]
        print mol.toSMILES()
        print inchi
        print keystart
        print thermo
        print jthermo
        
        print 'H', thermo.H298.value/jthermo.H298.value, thermo.H298.value-jthermo.H298.value 
        print 'S', thermo.S298.value/jthermo.S298.value, thermo.S298.value-jthermo.S298.value
        print 'Cp', thermo.Cpdata.values/jthermo.Cpdata.values
        print 
        
        
                