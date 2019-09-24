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

"""
This module contains the QMVerifier class for confirming quantum job success.
"""

import logging
import os


class QMVerifier(object):
    """
    Verifies whether a QM job (externalized) was succesfully completed by 
      * searching for specific keywords in the output files, 
      * located in a specific directory (e.g. "QMFiles")
    """

    def __init__(self, molfile):
        self.molfile = molfile
        self.gaussianResultExists = False
        self.mopacResultExists = False
        self.mm4ResultExists = False

        self.outputExtension = '.out'
        self.inputExtension = '.mop'

    def check_for_inchi_key_collision(self, log_file_inchi):
        """
        This method is designed in the case a MOPAC output file was found but the InChI found in the file did not
        correspond to the InChI of the given molecule.
        
        This could mean two things:
        1) that the InChI Key hash does not correspond to the InChI it is hashed from. This is the rarest case of them all
        2) the complete InChI did not fit onto just one line in the MOPAC output file. Therefore it was continued on the
        second line and only a part of the InChI was actually taken as the 'whole' InChI.
        
        This method reads in the MOPAC input file and compares the found InChI in there to the InChI of the given molecule.
        """
        # if InChIPartialMatch == 1:#case where the InChI in memory begins with the InChI in the log file we will continue and check the input file, pring a warning if there is no match
        # look in the input file if the InChI doesn't match (apparently, certain characters can be deleted in MOPAC output file for long InChIs)
        input_file = os.path.join(self.molfile.directory, self.molfile.name + self.inputExtension)

        assert os.path.exists(input_file)

        with open(input_file) as input_file:  # read the MOPAC input_file
            lineI = input_file.readline()
            for line in input_file:
                if line.startswith("InChI="):
                    input_file_inchi = lineI.rstrip()
                    break

            if input_file_inchi == self.molfile.InChIAug:
                logging.info('An InChI match could be found in the input file, but not in the output file. Anywho, a\
                pre-existing successful MOPAC result exists.')
                return True

            else:
                logging.info("*****Warning: potential InChIKey collision: InChIKey(augmented) = " +
                             self.molfile.name + " RMG Augmented InChI = " + self.molfile.InChIAug +
                             " Log file Augmented InChI = " + log_file_inchi +
                             " . InChI could not be found in the MOPAC input file. You should manually check that the output file contains the ended species.")
                return False

        # returns True if an MM4 output file for the given name and directory (.mm4out suffix) exists and indicates successful completion (same criteria as used after calculation runs) terminates if the InChI doesn't match the InChI in the file or if there is no InChI in the file returns False otherwise

    def succesful_job_exists(self):
        """
        checks whether one of the flags is true.
        If so, it returns true.
        """
        return self.gaussianResultExists or self.mopacResultExists or self.mm4ResultExists
