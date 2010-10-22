#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os

rmgpy_dir = os.environ.get('RMGpy')
if rmgpy_dir is None: raise RuntimeError("Please set RMGpy environment variable")
sys.path.append(os.path.join(rmgpy_dir,'source'))

import rmg.thermo
from rmg.structure import *
from rmg.species import *

import pybel
sys.path.append(os.path.join(rmgpy_dir,'source','external','cinfony'))
import webel


for line in file('input.dat'):
	if line.startswith("GEOM File"):
		logfile = line.replace("GEOM File","").strip()

pymol = pybel.readfile("g03",logfile).next()

s=Structure()
s.fromOBMol(pymol.OBMol)
smiles = s.toSMILES()
print "The input species is probably %s"%smiles

print "The chemgraph for that is"
print s.toAdjacencyList(strip_hydrogens=True)

outfilename='chemgraph.txt'

outfile = file(outfilename,'w')
outfile.write(s.toAdjacencyList(strip_hydrogens=True))
outfile.close()

webmol = webel.readstring("smi", smiles)
names = webmol.write("names")
print "Some names for this include:"
for n in names: print " ",n

print "This has been written to %s"%outfilename
print "Please check that it is correct, especially for radicals and singlet/triplet issues"




################################################################################

if __name__ == '__main__':
	pass