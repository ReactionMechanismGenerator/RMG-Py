#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import os.path

import rmgpy.constants as constants

class MoleProLog:
    """
    Represents a MolePro log file. The attribute `path` refers to the
    location on disk of the MolePro log file of interest. Methods are provided
    to extract a variety of information into CanTherm classes and/or NumPy
    arrays. 
    """
    
    def __init__(self, path):
        self.path = path
            
    def loadCCSDEnergy(self):
        """
        Return the f12 energy in J/mol from a MolePro Logfile of a CCSD(T)-f12 job. 
        This function determines which energy (f12a or f12b) to use based on the basis set,
        which it will parse out of the molepro file. For the vtz and dtz basis sets f12a is
        better approximation, but for higher basis sets f12b is a better approximation
        """
        f = open(self.path, 'r')
        line=f.readline()
        
        #search for basisSet
        while line!='':
            if 'basis' in line.lower():
                if 'vtz' in line.lower() or'vdz' in line.lower():
                    f12a=True
                else: f12a=False
                break
            line=f.readline()
        else: raise Exception('Could not find basis set in MolePro File')
        #search for energy
        E0=None
        if f12a:
            while line!='':
                #first one is for radicals second is for non radicals
                if 'RHF-UCCSD(T)-F12a energy' in line or 'CCSD(T)-F12a total energy  ' in line:
                    E0=float(line.split()[-1])
                    break
                line=f.readline()
        else:
            while line!='':
                if 'RHF-UCCSD(T)-F12b energy' in line or 'CCSD(T)-F12b total energy  ' in line:
                    E0=float(line.split()[-1])
                    break
                line=f.readline()
        
        f.close()
        
        #multiply E0 by correct constants
        if E0 is not None:
            E0 = E0 * constants.E_h * constants.Na
            return E0
        else: raise Exception('Unable to find energy in MolePro log file.')