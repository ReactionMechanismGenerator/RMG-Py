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
from rmgpy.quantity import Quantity
from rmgpy.molecule.element import elementList

from arkane.util import determine_qm_software
from arkane.output import prettify
from arkane.log import Log
from arkane.common import symbol_by_number


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
   
   q2dtor_path = os.path.join(os.path.split(os.path.abspath(os.path.dirname(__file__)))[0], 'external', 'Q2DTor', 'src', 'Q2DTor.py')
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
       determine torsions, not entirely necessary for 2D-NS (2D rotor not separated), but important for E2DT (couples frequencies with a 2D-NS rotor)
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
