#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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

import distutils.spawn
import logging
import os
import re
import shutil
import tempfile
from subprocess import Popen, PIPE

import cclib

from rmgpy.exceptions import DependencyError
from rmgpy.molecule.molecule import Molecule
from rmgpy.qm.molecule import QMMolecule


class Mopac(object):
    """
    A base class for all QM calculations that use MOPAC.
    
    Classes such as :class:`MopacMol` will inherit from this class.
    """

    inputFileExtension = '.mop'
    outputFileExtension = '.out'

    executablesToTry = ('MOPAC2016.exe', 'MOPAC2012.exe', 'MOPAC2009.exe', 'mopac')

    for exe in executablesToTry:
        try:
            executablePath = distutils.spawn.find_executable(exe)
        except:
            executablePath = None
        if executablePath is not None:
            break
    else:  # didn't break
        logging.debug("Did not find MOPAC on path, checking if it exists in a declared MOPAC_DIR...")
        mopacEnv = os.getenv('MOPAC_DIR', default="/opt/mopac")
        for exe in executablesToTry:
            executablePath = os.path.join(mopacEnv, exe)
            if os.path.exists(executablePath):
                break
        else:  # didn't break
            executablePath = os.path.join(mopacEnv, '(MOPAC 2009 or 2012 or 2016)')

    usePolar = False  # use polar keyword in MOPAC

    "Keywords for the multiplicity"
    multiplicityKeywords = {
        1: '',
        2: 'uhf doublet',
        3: 'uhf triplet',
        4: 'uhf quartet',
        5: 'uhf quintet',
        6: 'uhf sextet',
        7: 'uhf septet',
        8: 'uhf octet',
        9: 'uhf nonet',
    }

    #: List of phrases that indicate failure
    #: NONE of these must be present in a succesful job.
    failureKeys = [
        'IMAGINARY FREQUENCIES',
        'EXCESS NUMBER OF OPTIMIZATION CYCLES',
        'NOT ENOUGH TIME FOR ANOTHER CYCLE',
    ]
    #: List of phrases to indicate success.
    #: ALL of these must be present in a successful job.
    successKeys = [
        'DESCRIPTION OF VIBRATIONS',
        'MOPAC DONE'
    ]

    def test_ready(self):
        if not os.path.exists(self.executablePath):
            raise DependencyError("Couldn't find MOPAC executable at {0}. Try setting your MOPAC_DIR "
                                  "environment variable.".format(self.executablePath))

        # Check if MOPAC executable works properly
        process = Popen(self.executablePath,
                        stdin=PIPE,
                        stdout=PIPE,
                        stderr=PIPE)
        stdout, stderr = process.communicate()

        self.expired = False
        stderr = stderr.decode('utf-8')
        if 'has expired' in stderr:
            # The MOPAC executable is expired
            logging.warning('\n'.join(stderr.split('\n')[2:7]))
            self.expired = True
        elif 'To install the MOPAC license' in stderr:
            # The MOPAC executable exists, but the license has not been installed
            raise DependencyError('\n'.join(stderr.split('\n')[0:9]))
        elif 'MOPAC_LICENSE' in stderr:
            # The MOPAC executable is in the wrong location on Windows; MOPAC_LICENSE must be set
            raise DependencyError('\n'.join(stderr.split('\n')[0:11]))

    def run(self):
        self.test_ready()
        # submits the input file to mopac

        dirpath = tempfile.mkdtemp()
        # copy input file to temp dir:
        tempInpFile = os.path.join(dirpath, os.path.basename(self.input_file_path))
        shutil.copy(self.input_file_path, dirpath)

        process = Popen([self.executablePath, tempInpFile], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        command = b'\n' if self.expired else None  # press enter to pass expiration notice
        stdout, stderr = process.communicate(input=command)  # necessary to wait for executable termination!
        if b"ended normally" not in stderr.strip():
            logging.warning("Mopac error message:" + stderr.decode('utf-8'))

        # copy output file from temp dir to output dir:
        tempOutFile = os.path.join(dirpath, os.path.basename(self.output_file_path))
        shutil.copy(tempOutFile, self.output_file_path)

        # delete temp folder:
        shutil.rmtree(dirpath)
        return self.verify_output_file()

    def verify_output_file(self):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful MOPAC simulation already exists for the molecule with the 
        given (augmented) InChI Key.
        
        The definition of finding a successful simulation is based on these criteria:
        1) finding an output file with the file name equal to the InChI Key
        2) NOT finding any of the keywords that are denote a calculation failure
        3) finding all the keywords that denote a calculation success.
        4) finding a match between the InChI of the given molecule and the InchI found in the calculation files
        5) checking that the optimized geometry, when connected by single bonds, is isomorphic with self.molecule (converted to single bonds)
        
        If any of the above criteria is not matched, False will be returned.
        If all succeed, then it will return True.
        """
        if not os.path.exists(self.output_file_path):
            logging.debug("Output file {0} does not (yet) exist.".format(self.output_file_path))
            return False

        inchi_found = False  # flag (1 or 0) indicating whether an InChI was found in the log file

        # Initialize dictionary with "False"s 
        success_keys_found = dict([(key, False) for key in self.successKeys])

        with open(self.output_file_path) as outputFile:
            for line in outputFile:
                line = line.strip()

                for element in self.failureKeys:  # search for failure keywords
                    if element in line:
                        logging.error("MOPAC output file contains the following error: {0}".format(element))
                        return False

                for element in self.successKeys:  # search for success keywords
                    if element in line:
                        success_keys_found[element] = True

                if "InChI=" in line:
                    log_file_inchi = line  # output files should take up to 240 characters of the name in the input file
                    inchi_found = True
                    if self.unique_id_long in log_file_inchi:
                        pass
                    elif self.unique_id_long.startswith(log_file_inchi):
                        logging.info("InChI too long to check, but beginning matches so assuming OK.")

                    else:
                        logging.warning("InChI in log file ({0}) didn't match that in geometry "
                                        "({1}).".format(log_file_inchi, self.unique_id_long))
                        # Use only up to first 80 characters to match due to MOPAC bug which deletes 81st character of InChI string
                        if self.unique_id_long.startswith(log_file_inchi[:80]):
                            logging.warning("but the beginning matches so it's probably just a truncation problem.")

        # Check that ALL 'success' keywords were found in the file.
        if not all(success_keys_found.values()):
            logging.error('Not all of the required keywords for success were found in the output file!')
            return False

        if not inchi_found:
            logging.error("No InChI was found in the MOPAC output file {0}".format(self.output_file_path))
            return False

        # Compare the optimized geometry to the original molecule
        qm_data = self.parse()
        cclib_mol = Molecule()
        cclib_mol.from_xyz(qm_data.atomicNumbers, qm_data.atomCoords.value)
        test_mol = self.molecule.to_single_bonds()
        if not cclib_mol.is_isomorphic(test_mol):
            logging.info("Incorrect connectivity for optimized geometry in file {0}".format(self.output_file_path))
            return False

        logging.info("Successful {1} quantum result in {0}".format(self.output_file_path, self.__class__.__name__))
        return True

    def get_parser(self, output_file):
        """
        Returns the appropriate cclib parser.
        """
        return cclib.parser.MOPAC(output_file)


class MopacMol(QMMolecule, Mopac):
    """
    A base Class for calculations of molecules using MOPAC. 
    
    Inherits from both :class:`QMMolecule` and :class:`Mopac`.
    """

    #: Keywords that will be added at the top and bottom of the qm input file
    keywords = [
        {'top': "precise nosym THREADS=1", 'bottom': "oldgeo thermo nosym precise THREADS=1 "},
        {'top': "precise nosym gnorm=0.0 nonr THREADS=1", 'bottom': "oldgeo thermo nosym precise THREADS=1 "},
        {'top': "precise nosym gnorm=0.0 THREADS=1", 'bottom': "oldgeo thermo nosym precise THREADS=1 "},
        {'top': "precise nosym gnorm=0.0 bfgs THREADS=1", 'bottom': "oldgeo thermo nosym precise THREADS=1 "},
        {'top': "precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000 THREADS=1", 'bottom': "oldgeo thermo nosym precise THREADS=1 "},
    ]

    def write_input_file(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attempt`.
        """

        molfile = self.get_mol_file_path_for_calculation(attempt)
        atomline = re.compile('\s*([\- ][0-9.]+)\s+([\- ][0-9.]+)+\s+([\- ][0-9.]+)\s+([A-Za-z]+)')

        output = [self.geometry.unique_id_long, '']

        atom_count = 0
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(4), match.group(1),
                                                                    match.group(2), match.group(3)))
                    atom_count += 1
        assert atom_count == len(self.molecule.atoms)

        output.append('')
        input_string = '\n'.join(output)

        top_keys, bottom_keys, polar_keys = self.input_file_keywords(attempt)
        with open(self.input_file_path, 'w') as mopac_file:
            mopac_file.write(top_keys)
            mopac_file.write('\n')
            mopac_file.write(input_string)
            mopac_file.write('\n')
            mopac_file.write(bottom_keys)
            if self.usePolar:
                mopac_file.write('\n\n\n')
                mopac_file.write(polar_keys)

    def input_file_keywords(self, attempt):
        """
        Return the top, bottom, and polar keywords.
        """
        raise NotImplementedError("Should be defined by subclass, eg. MopacMolPM3")

    def generate_qm_data(self):
        """
        Calculate the QM data and return a QMData object, or None if it fails.
        """
        for atom in self.molecule.vertices:
            if atom.charge != 0:
                return None

        if self.verify_output_file():
            logging.info("Found a successful output file already; using that.")
            source = "QM {0} calculation found from previous run.".format(self.__class__.__name__)
        else:
            self.create_geometry()
            success = False
            for attempt in range(1, self.max_attempts + 1):
                self.write_input_file(attempt)
                logging.info('Trying {3} attempt {0} of {1} on molecule {2}.'.format(attempt, self.max_attempts,
                                                                                     self.molecule.to_smiles(),
                                                                                     self.__class__.__name__))
                success = self.run()
                if success:
                    logging.info('Attempt {0} of {1} on species {2} succeeded.'.format(attempt, self.max_attempts,
                                                                                       self.molecule.to_augmented_inchi()))
                    source = "QM {0} calculation attempt {1}".format(self.__class__.__name__, attempt)
                    break
            else:
                logging.error('QM thermo calculation failed for {0}.'.format(self.molecule.to_augmented_inchi()))
                return None
        result = self.parse()  # parsed in cclib
        result.source = source
        return result


class MopacMolPMn(MopacMol):
    """
    Mopac PMn calculations for molecules (n undefined here)
    
    This is a parent class for MOPAC PMn calculations.
    Inherit it, and define the pm_method, then redefine 
    anything you wish to do differently.
    """
    pm_method = '(should be defined by sub class)'

    def input_file_keywords(self, attempt):
        """
        Return the top, bottom, and polar keywords for attempt number `attempt`.
        
        NB. `attempt` begins at 1, not 0.
        """
        assert attempt <= self.max_attempts

        if attempt > self.script_attempts:
            attempt -= self.script_attempts

        multiplicity_keys = self.multiplicityKeywords[self.geometry.molecule.multiplicity]

        top_keys = "{method} {mult} {top}".format(
            method=self.pm_method,
            mult=multiplicity_keys,
            top=self.keywords[attempt - 1]['top'],
        )
        bottom_keys = "{bottom} {method} {mult}".format(
            method=self.pm_method,
            bottom=self.keywords[attempt - 1]['bottom'],
            mult=multiplicity_keys,
        )
        polar_keys = "oldgeo {polar} nosym precise {method} {mult}".format(
            method=self.pm_method,
            polar=('polar' if self.geometry.molecule.multiplicity == 1 else 'static'),
            mult=multiplicity_keys,
        )

        return top_keys, bottom_keys, polar_keys


class MopacMolPM3(MopacMolPMn):
    """
    Mopac PM3 calculations for molecules
    
    This is a class of its own in case you wish to do anything differently,
    but for now it's the same as all the MOPAC PMn calculations, only pm3
    """
    pm_method = 'pm3'


class MopacMolPM6(MopacMolPMn):
    """
    Mopac PM6 calculations for molecules
    
    This is a class of its own in case you wish to do anything differently,
    but for now it's the same as all the MOPAC PMn calculations, only pm6
    """
    pm_method = 'pm6'


class MopacMolPM7(MopacMolPMn):
    """
    Mopac PM7 calculations for molecules
    
    This is a class of its own in case you wish to do anything differently,
    but for now it's the same as all the MOPAC PMn calculations, only pm7
    """
    pm_method = 'pm7'
