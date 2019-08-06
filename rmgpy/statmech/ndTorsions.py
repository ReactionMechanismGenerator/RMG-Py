#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import os.path
import numpy as np
import logging
import subprocess
import os
from scipy import interpolate
from scipy import integrate as inte
import numdifftools as nd
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model

import rmgpy.constants as constants
from rmgpy.statmech.conformer import Conformer
from rmgpy.statmech import schrodinger
from rmgpy.statmech.mode import Mode
from rmgpy.molecule.element import elementList

class HinderedRotor2D(Mode):
   """
   A statistical mechanical model of a 2D-dimensional hindered rotor.
   Computation of eigenvalues and fourier fitting outsourced to
   Q2DTor software

   The attributes are:

   ======================== ===================================================
   Attribute                Description
   ======================== ===================================================
   `name`                   The Q2DTor name of the rotor
   `torsion1`               The 1-indexed atom indices of the atoms involved in the first rotor
   `torsion2`               The 1-indexed atom indices of the atoms involved in the second rotor
   `torsigma1`              The symmetry number of the first rotor
   `torsigma2`              The symmetry number of the second rotor
   `calcPath`               The directory containing all of the rotor point calculations formated:  name_angle1_angle2
   `symmetry`               Q2DTor symmetry identifier ('none','a','b','c','ab','bc','ac','abc')
   `evals`                  Array of energy levels for the 2D-HR
   `energy`                 Function mapping quantum number to energy level
   ======================== ===================================================
   """
   
   q2dtor_path = os.path.join(os.path.split(os.path.split(os.path.abspath(os.path.dirname(__file__)))[0])[0], 'external', 'Q2DTor', 'src', 'Q2DTor.py')
   q2dtor_message = """\nUsing Q2DTor...
   Q2DTor is a software for calculating the partition functions and themodynamic properties of molecular systems with two or more
   torsional modes developed by David Ferro Costas (david.ferro@usc.es) and Antonio Fernandez Ramos (qf.ramos@usc.es) at
   the Universidade de Santiago de Compostela. Arkane can integrate Q2DTor to compute the quantum mechanical partition function 
   of 2D rotors.  

   For use of HinderedRotor2D within Arkane please cite:  
   D. Ferro-Costas, M. N. D. S. Cordeiro, D. G. Truhlar, A. Fernández-Ramos, Comput. Phys. Commun. 232, 190-205, 2018.
   """
   q2dtor_message_used = False
   
   def __init__(self, name, torsigma1, torsigma2, calcPath, symmetry='none', pivots1=None, pivots2=None, top1=None, top2=None):
       Mode.__init__(self, True)
       self.dof = 2

       self.name = name
       self.calcPath = calcPath

       self.torsion1 = None
       self.torsion2 = None
       self.torsigma1 = torsigma1
       self.torsigma2 = torsigma2
       self.symmetry = symmetry
       self.charge = None
       self.multiplicity = None

       self.pivots1 = pivots1
       self.pivots2 = pivots2
       self.top1 = top1
       self.top2 = top2

       self.xyzs = []
       self.phi1s = []
       self.phi2s = []
       self.Es = []
       self.atnums = []
       self.element_names = []
       
       if not self.q2dtor_message_used:
           logging.info(self.q2dtor_message)
           self.q2dtor_message_used = True

       # prepare a directory to run q2dtor in
       self.q2dtor_dir = os.path.join(os.path.split(calcPath)[0],name)
       try:
           os.mkdir(self.q2dtor_dir)
           os.mkdir(os.path.join(self.q2dtor_dir,'IOfiles'))
       except OSError:
           pass

   def getTorsions(self):
       """
       determine torsions, not entirely necessary for 2D-NS (2D rotor not separated), 
       but important for E2DT (couples frequencies with a 2D-NS rotor)
       """
       if not self.torsion1:
           self.readGjf()  # check if there is a gaussian format file

       Natoms = len(self.xyzs[0])  # define a feasible torsion from pivots and tops
       aset = set(list(xrange(1,Natoms+1)))
       if not self.torsion1 and self.pivots1 and self.top1:
           if self.pivots1[0] in self.top1:
               self.torsion1 = [list(aset-set(self.top1))[0],self.pivots1[0],self.pivots1[1],self.top1[0]]
           else:
               self.torsion1 = [list(aset-set(self.top1))[0],self.pivots1[1],self.pivots1[0],self.top1[0]]

       if not self.torsion2 and self.pivots2 and self.top2:
           if self.pivots2[0] in self.top2:
               self.torsion2 = [list(aset-set(self.top2))[0],self.pivots2[0],self.pivots2[1],self.top2[0]]
           else:
               self.torsion2 = [list(aset-set(self.top2))[0],self.pivots2[1],self.pivots2[0],self.top2[0]]

   def readScan(self):
       """
       Read quantum optimization job files at self.calcPath to determine
       vectors of angles (self.phi1s, self.phi2s), xyz coordinates (self.xyzs)
       energies (self.Es) and atom numbers (self.atnums) for each point
       """
       from arkane.util import determine_qm_software
       from arkane.common import symbol_by_number
       phi1s = []
       phi2s = []
       xyzs = []
       Es = []
       atnums = []
       for f in os.listdir(self.calcPath):
           if len(f.split('_')) != 4:
               continue
           s, name, phi1, phi2 = f.split('_') # scangeom_r0_0.0_360.0.log
           phi2,identifier = '.'.join(phi2.split('.')[:-1]),phi2.split('.')[-1]
           if identifier != 'out':
               continue
           phi1s.append(float(phi1))
           phi2s.append(float(phi2.split(".")[0]))

           fpath = os.path.join(self.calcPath,f)
           lg = determine_qm_software(fpath)
           
           Es.append(lg.loadEnergy())
           xyz,atnums,_ = lg.loadGeometry()
           xyzs.append(xyz)

       self.xyzs = xyzs
       self.phi1s = phi1s
       self.phi2s = phi2s
       self.Es = Es
       self.atnums = atnums
       self.element_names = [symbol_by_number[k] for k in self.atnums]
       
   def readGjf(self):
       """
       read gaussian input file to determine torsions, charge and multiplicity
       unnecessary for 2D-NS
       """
       for f in os.listdir(self.calcPath):
           if len(f.split('_')) != 4:
               continue
           s, name, phi1, phi2 = f.split('_') # scangeom_r0_0.0_360.0.log
           phi2,identifier = phi2.split('.')
           if identifier == 'gjf' and float(phi1) == 0.0 and float(phi2) == 0.0:
               with open(os.path.join(self.calcPath,f), 'r') as fop:
                   lines = fop.readlines()
                   for i,line in enumerate(lines):
                       split_line = line.split()
                       if len(split_line) < 2:
                           continue
                       elif split_line[-1] == 'F' and len(split_line) == 5: # F 1 2 11 14
                           if not self.torsion1:
                               self.torsion1 = [int(split_line[i]) for i in xrange(4)]
                           elif not self.torsion2:
                               self.torsion2 = [int(split_line[i]) for i in xrange(4)]
                       elif (not self.charge or not self.multiplicity) and len(split_line) == 2 and \
                             len(lines[i+1].split()) == 4 and len(lines[i+1].split()[0]) <= 2:
                           self.charge = int(split_line[0])
                           self.multiplicity = int(split_line[1])

   def writeXYZ(self):
       """
       write an .xyz file for Q2DTor
       done based on the angle coordinates (0.0,0.0)
       """
       for i in xrange(len(self.phi1s)):
           if self.phi1s[i] == 0.0 and self.phi2s[i] == 0.0:
               with open(os.path.join(self.q2dtor_dir,self.name+".xyz"),'w') as f:
                   f.write(str(len(self.atnums))+'\n')
                   f.write("reference geometry for {0}\n".format(self.name))
                   for j,ename in enumerate(self.element_names):
                       f.write('{ename}    {x}   {y}   {z}\n'.format(ename=ename,
                               x=self.xyzs[i][j][0],y=self.xyzs[i][j][1],z=self.xyzs[i][j][2]))
                   break

   def writePes(self):
       """
       write a .pes file for Q2DTor based on the
       read in scans
       """
       if not len(self.Es) > 0:
           raise ValueError("Cannot write PES file with no scan information")
       with open(os.path.join(self.q2dtor_dir,'IOfiles',self.name+".pes"),'w') as f:
           for i in xrange(len(self.phi1s)):
               f.write(str(len(self.atnums))+'\n')
               f.write("Geometry   {E}   {phi1}   {phi2}  {name}_{phi1}_{phi2}  YES\n".format(E=self.Es[i]/(constants.E_h*constants.Na),
                       phi1=self.phi1s[i],phi2=self.phi2s[i],name=self.name))
               for j,ename in enumerate(self.element_names):
                   f.write('{ename}    {x}   {y}   {z}\n'.format(ename=ename,
                               x=self.xyzs[i][j][0],y=self.xyzs[i][j][1],z=self.xyzs[i][j][2]))

   def writeInp(self):
       """
       Write an input file for Q2DTor based on object
       information
       """
       inp = """#----------------------------------#
# Torsional information            #
#----------------------------------#
start_torsions
   torsion1         {tor1}         # atoms involved in torsion 1
   torsion2         {tor2}         # atoms involved in torsion 2
   tsigma1          {torsigma1}               # torsional symmetry number of hindered rotor 1
   tsigma2          {torsigma2}               # torsional symmetry number of hindered rotor 2
end_torsions
#----------------------------------#
# Calculations                     #
#----------------------------------#
start_calcs
   level            hf 3-21g        # the calculation level
   charge           0               # charge
   multiplicity     0               # spin multiplicity, shouldn't matter in context of Arkane
end_calcs
#----------------------------------#
# Torsional PES                    #
#----------------------------------#
start_pes
   t1step           10.0            # step in phi1 for scan calculation [degrees]
   t2step           10.0            # step in phi2 for scan calculation [degrees]
   symmetry         {symmetry}               # Symmetry condition for PES: [a,b,c,ab,ac,bc,abc] or none
end_pes
#----------------------------------#
# Fitting details                  #
#----------------------------------#
start_fourier
   weight           0.9             #
   ignore           0.0             # Set to zero coefficients with smaller absolute value (in cm^-1)
   # Fourier Terms (Even)           #
   cos1             1-9             # i values in cos(i*Phi_1)
   cos2             1-9             # j values in cos(j*Phi_2)
   cos1cos2         1-7 , 1-7       # i,j values in cos(i*Phi_1) * cos(j*Phi_2)
   sin1sin2         1-7 , 1-7       # i,j values in sin(i*Phi_1) * sin(j*Phi_2)
   # Fourier Terms (Odd)            #
   sin1             1-9            # i values in sin(i*Phi_1)
   sin2             1-9            # j values in sin(j*Phi_2)
   cos1sin2         1-7 , 1-7            # i,j values in cos(i*Phi_1) * sin(j*Phi_2)
   sin1cos2         1-7 , 1-7            # i,j values in sin(i*Phi_1) * cos(j*Phi_2)
end_fourier
#----------------------------------#
# Search and Opt stationary points #
#----------------------------------#
start_statpoint
   tolerance        1.0             # step (in degrees) to explore torsional PES when looking for CPs
   freqscal         1.000           # scaling factor for frequencies
end_statpoint
#----------------------------------#
# 2D-NS Hamiltonian                #
#----------------------------------#
start_tor2dns
   dijvar           yes             # yes (dij not constant) or no (dij constant)
   kmax             100             # check 2013-JChemPhys_138_134112, eq (14)
   maxeigen         1e4             # threshold for H eigenvalues (in cm^-1)
end_tor2dns
#----------------------------------#
# Partition functions              #
#----------------------------------#
start_rovibpf
   interpolation    fourier         # fourier or spline order (1,3,5)
   integrationstep  1.0             # integration dphi
end_rovibpf
#----------------------------------#
# Working temperatures             #
#----------------------------------#
start_temperatures                 #
   100.0   150.0   200.0          #
   250.0   300.0   400.0          #
   500.0   700.0  1000.0          #
   1500.0  2000.0  2500.0         #
end_temperatures                   #
#----------------------------------#""".format(tor1="-".join([str(x) for x in self.torsion1]),
                  tor2="-".join([str(x) for x in self.torsion2]),
           torsigma1=self.torsigma1,torsigma2=self.torsigma2,
           charge=self.charge,multiplicity=self.multiplicity,
           symmetry=self.symmetry)
       f = open(os.path.join(self.q2dtor_dir,self.name+".inp"),'w')
       f.write(inp)
       f.close()

   def getIcsFile(self):
       """
       use Q2DTor to generate a .ics file the Q2DTor file that
       has torsional information
       """
       out = subprocess.check_call(['python2',self.q2dtor_path,self.name,'--init'],
                       cwd=self.q2dtor_dir)

   def fitFourier(self):
       """
       use Q2DTor to fit fourier coefficients
       to the potential
       """
       out = subprocess.check_call(['python2',self.q2dtor_path,self.name,'--fourier'],
                       cwd=self.q2dtor_dir)

   def getSplistfile(self):
       """
       use Q2DTor to generate a .splist file
       """
       out = subprocess.check_call(['python2',self.q2dtor_path,self.name,'--findsp'],
                       cwd=self.q2dtor_dir)

   def getEigvals(self):
       """
       use Q2DTor to determine the QM energy levels for the 2D-NS
       rotors
       writes a .evals file and reads it to fill self.evals and self.energy
       """
       out = subprocess.check_call(['python2',self.q2dtor_path,self.name,'--tor2dns'],
                       cwd=self.q2dtor_dir)
       self.readEigvals()

   def readEigvals(self):
       """
       reads an available .evals file to get the QM energy levels
       for the 2D-NS rotors
       """
       with open(os.path.join(self.q2dtor_dir,'IOfiles',self.name+'.evals'),'r') as f:
           out = f.readlines()
           evals = [float(x.split()[1]) for x in out[2:]]  # cm^-1
           self.evals = np.array(evals)*10**2*constants.c*constants.h*constants.Na  # J/mol
           self.energy = lambda x: self.evals[x]

   def getPartitionFunction(self, T):
       return schrodinger.getPartitionFunction(T,self.energy,nmax=len(self.evals))

   def getHeatCapacity(self, T):
       return schrodinger.getHeatCapacity(T,self.energy,nmax=len(self.evals))

   def getEnthalpy(self, T):
       return schrodinger.getEnthalpy(T,self.energy,nmax=len(self.evals))

   def getEntropy(self, T):
       return schrodinger.getEntropy(T,self.energy,nmax=len(self.evals))

   def getSumOfStates(self, Elist, sumStates0=None):
       if sumStates0:
           return schrodinger.getSumOfStates(Elist,self.energy,sumStates0=sumStates0,nmax=len(self.evals))
       else:
           return schrodinger.getSumOfStates(Elist,self.energy,nmax=len(self.evals))

   def getDensityOfStates(self, Elist, densStates0=None):
       if densStates0:
           return schrodinger.getDensityOfStates(Elist,self.energy,densStates0=densStates0,nmax=len(self.evals))
       else:
           return schrodinger.getDensityOfStates(Elist,self.energy,nmax=len(self.evals))
