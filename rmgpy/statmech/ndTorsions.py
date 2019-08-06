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

 