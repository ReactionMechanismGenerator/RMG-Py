#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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
import itertools
import logging
import os
import re
from subprocess import Popen

import cclib

from rmgpy.molecule.molecule import Molecule
from rmgpy.qm.molecule import QMMolecule
from rmgpy.qm.qmdata import parse_cclib_data


class Gaussian:
    """
    A base class for all QM calculations that use Gaussian.

    Classes such as :class:`GaussianMol` will inherit from this class.
    """

    input_file_extension = ".gjf"
    output_file_extension = ".log"

    executablesToTry = ("g16", "g09", "g03")

    for exe in executablesToTry:
        try:
            executable_path = distutils.spawn.find_executable(exe)
        except:
            executable_path = None
        if executable_path is not None:
            break
    else:  # didn't break
        logging.debug(
            "Did not find Gaussian on path, checking if it exists in a declared GAUSS_EXEDIR, g09root or g03root..."
        )
        gaussEnv = (
            os.getenv("GAUSS_EXEDIR")
            or os.getenv("g09root")
            or os.getenv("g03root")
            or ""
        )
        possibleDirs = gaussEnv.split(
            ":"
        )  # GAUSS_EXEDIR may be a list like "path1:path2:path3"
        for exe, possibleDir in itertools.product(executablesToTry, possibleDirs):
            executable_path = os.path.join(possibleDir, exe)
            if os.path.exists(executable_path):
                break
        else:  # didn't break
            executable_path = os.path.join(gaussEnv, "(Gaussian 2003 or 2009)")

    use_polar = False

    #: List of phrases that indicate failure
    #: NONE of these must be present in a succesful job.
    failure_keys = ["ERROR TERMINATION", "IMAGINARY FREQUENCIES"]

    #: List of phrases to indicate success.
    #: ALL of these must be present in a successful job.
    success_keys = ["Normal termination of Gaussian"]

    def test_ready(self):
        if not os.path.exists(self.executable_path):
            raise Exception(
                "Couldn't find Gaussian executable at {0}. "
                "Try setting your GAUSS_EXEDIR environment variable.".format(
                    self.executable_path
                )
            )

    def run(self):
        self.test_ready()
        # submits the input file to Gaussian
        process = Popen(
            [self.executable_path, self.input_file_path, self.output_file_path]
        )
        process.communicate()  # necessary to wait for executable termination!

        return self.verify_output_file()

    def verify_output_file(self):
        """
        Check's that an output file exists and was successful.

        Returns a boolean flag that states whether a successful GAUSSIAN simulation already exists for the molecule with the
        given (augmented) InChI Key.

        The definition of finding a successful simulation is based on these criteria:
        1) finding an output file with the file name equal to the InChI Key
        2) NOT finding any of the keywords that are denote a calculation failure
        3) finding all the keywords that denote a calculation success.
        4) finding a match between the InChI of the given molecule and the InchI found in the calculation files
        5) checking that the optimized geometry, when connected by single bonds, is isomorphic with self.molecule (converted to single bonds)

        If any of the above criteria is not matched, False will be returned.
        If all are satisfied, it will return True.
        """
        if not os.path.exists(self.output_file_path):
            logging.info(
                "Output file {0} does not exist.".format(self.output_file_path)
            )
            return False

        inchi_match = False  # flag (1 or 0) indicating whether the InChI in the file matches InChIaug this can only be 1 if inchi_found is also 1
        inchi_found = (
            False  # flag (1 or 0) indicating whether an InChI was found in the log file
        )

        # Initialize dictionary with "False"s
        success_keys_found = dict([(key, False) for key in self.success_keys])

        with open(self.output_file_path) as outputFile:
            for line in outputFile:
                line = line.strip()

                for element in self.failure_keys:  # search for failure keywords
                    if element in line:
                        logging.error(
                            "Gaussian output file contains the following error: {0}".format(
                                element
                            )
                        )
                        return False

                for element in self.success_keys:  # search for success keywords
                    if element in line:
                        success_keys_found[element] = True

                if line.startswith("InChI="):
                    log_file_inchi = line  # output files should take up to 240 characters of the name in the input file
                    inchi_found = True
                    if self.unique_id_long in log_file_inchi:
                        inchi_match = True
                    elif self.unique_id_long.startswith(log_file_inchi):
                        logging.info(
                            "InChI too long to check, but beginning matches so assuming OK."
                        )
                        inchi_match = True
                    else:
                        logging.warning(
                            "InChI in log file ({0}) didn't match that in geometry "
                            "({1}).".format(
                                log_file_inchi, self.geometry.unique_id_long
                            )
                        )
                        if self.geometry.unique_id_long.startswith(log_file_inchi):
                            logging.warning(
                                "but the beginning matches so it's probably just a truncation problem."
                            )
                            inchi_match = True
        # Check that ALL 'success' keywords were found in the file.
        if not all(success_keys_found.values()):
            logging.error(
                "Not all of the required keywords for success were found in the output file!"
            )
            return False

        if not inchi_found:
            logging.error(
                "No InChI was found in the Gaussian output file {0}".format(
                    self.output_file_path
                )
            )
            return False

        if not inchi_match:
            # InChIs do not match (most likely due to limited name length mirrored in log file (240 characters), but possibly due to a collision)
            return self.checkForInChiKeyCollision(
                log_file_inchi
            )  # Not yet implemented!

        # Compare the optimized geometry to the original molecule
        qm_data = self.parse()
        cclib_mol = Molecule()
        cclib_mol.from_xyz(qm_data.atomicNumbers, qm_data.atomCoords.value)
        test_mol = self.molecule.to_single_bonds()
        if not cclib_mol.is_isomorphic(test_mol):
            logging.info(
                "Incorrect connectivity for optimized geometry in file {0}".format(
                    self.output_file_path
                )
            )
            return False

        logging.info(
            "Successful {1} quantum result in {0}".format(
                self.output_file_path, self.__class__.__name__
            )
        )
        return True

    def parse(self):
        """
        Parses the results of the Gaussian calculation, and returns a QMData object.
        """
        parser = cclib.parser.Gaussian(self.output_file_path)
        parser.logger.setLevel(
            logging.ERROR
        )  # cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
        cclib_data = parser.parse()
        radical_number = sum([i.radical_electrons for i in self.molecule.atoms])
        qm_data = parse_cclib_data(cclib_data, radical_number + 1)
        return qm_data


class GaussianMol(QMMolecule, Gaussian):
    """
    A base Class for calculations of molecules using Gaussian.

    Inherits from both :class:`QMMolecule` and :class:`Gaussian`.
    """

    def input_file_keywords(self, attempt):
        """
        Return the top keywords for attempt number `attempt`.

        NB. `attempt` begins at 1, not 0.
        """
        assert attempt <= self.max_attempts
        if attempt > self.script_attempts:
            attempt -= self.script_attempts
        return self.keywords[attempt - 1]

    def write_input_file(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attempt`.
        """
        molfile = self.get_mol_file_path_for_calculation(attempt)
        atomline = re.compile(
            r"\s*([\- ][0-9.]+\s+[\-0-9.]+\s+[\-0-9.]+)\s+([A-Za-z]+)"
        )

        output = ["", self.geometry.unique_id_long, ""]
        output.append(
            "{charge}   {mult}".format(
                charge=0, mult=(self.molecule.get_radical_count() + 1)
            )
        )

        atom_count = 0
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    output.append("{0:8s} {1}".format(match.group(2), match.group(1)))
                    atom_count += 1
        assert atom_count == len(self.molecule.atoms)

        output.append("")
        input_string = "\n".join(output)

        top_keys = self.input_file_keywords(attempt)
        with open(self.input_file_path, "w") as gaussian_file:
            gaussian_file.write(top_keys)
            gaussian_file.write("\n")
            gaussian_file.write(input_string)
            gaussian_file.write("\n")
            if self.use_polar:
                gaussian_file.write("\n\n\n")
                raise NotImplementedError("Not sure what should be here, if anything.")
                # gaussian_file.write(polar_keys)

    def generate_qm_data(self):
        """
        Calculate the QM data and return a QMData object.
        """
        # still can't handle charged atoms for QMData
        for atom in self.molecule.vertices:
            if atom.charge != 0:
                return None

        if self.verify_output_file():
            logging.info("Found a successful output file already; using that.")
            source = "QM {0} calculation found from previous run.".format(
                self.__class__.__name__
            )
        else:
            self.create_geometry()
            success = False
            for attempt in range(1, self.max_attempts + 1):
                self.write_input_file(attempt)
                logging.info(
                    "Trying {3} attempt {0} of {1} on molecule {2}.".format(
                        attempt,
                        self.max_attempts,
                        self.molecule.to_smiles(),
                        self.__class__.__name__,
                    )
                )
                success = self.run()
                if success:
                    logging.info(
                        "Attempt {0} of {1} on species {2} succeeded.".format(
                            attempt,
                            self.max_attempts,
                            self.molecule.to_augmented_inchi(),
                        )
                    )
                    source = "QM {0} calculation attempt {1}".format(
                        self.__class__.__name__, attempt
                    )
                    break
            else:
                logging.error(
                    "QM thermo calculation failed for {0}.".format(
                        self.molecule.to_augmented_inchi()
                    )
                )
                return None
        result = self.parse()  # parsed in cclib
        result.source = source
        return result  # a CCLibData object

    def get_parser(self, output_file):
        """
        Returns the appropriate cclib parser.
        """
        return cclib.parser.Gaussian(output_file)


class GaussianMolPM3(GaussianMol):
    """
    Gaussian PM3 calculations for molecules

    This is a class of its own in case you wish to do anything differently,
    but for now it's only the 'pm3' in the keywords that differs.
    """

    #: Keywords that will be added at the top of the qm input file
    keywords = [
        # The combinations of keywords were derived by Greg Magoon for pm3 in Gaussian. His comments are attached to each combination.
        "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3)",  # added IOP option to avoid aborting when symmetry changes; 3 is supposed to be default according to documentation, but it seems that 0 (the default) is the only option that doesn't work from 0-4; also, it is interesting to note that all 4 options seem to work for test case with z-matrix input rather than xyz coords; cf. http://www.ccl.net/cgi-bin/ccl/message-new?2006+10+17+005 for original idea for solution
        "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)",  # use different SCF method; this addresses at least one case of failure for a C4H7J species
        "# pm3 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm",  # try multiple different options (no gdiis, use calcfc, nosymm); 7/21/09: added maxcyc option to fix case of MPTBUKVAJYJXDE-UHFFFAOYAPmult3 (InChI=1/C4H10O5Si/c1-3-7-9-10(5,6)8-4-2/h4-5H,3H2,1-2H3/mult3) (file manually copied to speed things along)
        "# pm3 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm",  # numerical frequency keyword version of keyword #3; used to address GYFVJYRUZAKGFA-UHFFFAOYALmult3 (InChI=1/C6H14O6Si/c1-3-10-13(8,11-4-2)12-6-5-9-7/h6-7H,3-5H2,1-2H3/mult3) case; (none of the existing Gaussian or MOPAC combinations worked with it)
        "# pm3 opt=(verytight,gdiis,small) freq IOP(2/16=3)",  # somehow, this worked for problematic case of ZGAWAHRALACNPM-UHFFFAOYAF (InChI=1/C8H17O5Si/c1-3-11-14(10,12-4-2)13-8-5-7(9)6-8/h7-9H,3-6H2,1-2H3); (was otherwise giving l402 errors); even though I had a keyword that worked for this case, I manually copied the fixed log file to QMfiles folder to speed things along; note that there are a couple of very low frequencies (~5-6 cm^-1 for this case)
        "# pm3 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)",  # used for troublesome C5H7J2 case (similar error to C5H7J below); calcfc is not necessary for this particular species, but it speeds convergence and probably makes it more robust for other species
        "# pm3 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)",  # use numerical frequencies; this takes a relatively long time, so should only be used as one of the last resorts; this seemed to address at least one case of failure for a C6H10JJ species; 7/15/09: maxcyc=200 added to address GVCMURUDAUQXEY-UHFFFAOYAVmult3 (InChI=1/C3H4O7Si/c1-2(9-6)10-11(7,8)3(4)5/h6-7H,1H2/mult3)...however, result was manually pasted in QMfiles folder to speed things along
        "# pm3 opt=tight freq IOP(2/16=3)",  # this worked for problematic case of SZSSHFMXPBKYPR-UHFFFAOYAF (InChI=1/C7H15O5Si/c1-3-10-13(8,11-4-2)12-7-5-6-9-7/h7H,3-6H2,1-2H3) (otherwise, it had l402.exe errors); corrected log file was manually copied to QMfiles to speed things along; we could also add a freq=numerical version of this keyword combination for added robustness; UPDATE: see below
        "# pm3 opt=tight freq=numerical IOP(2/16=3)",  # used for problematic case of CIKDVMUGTARZCK-UHFFFAOYAImult4 (InChI=1/C8H15O6Si/c1-4-12-15(10,13-5-2)14-7-6-11-8(7,3)9/h7H,3-6H2,1-2H3/mult4 (most other cases had l402.exe errors); corrected log file was manually copied to QMfiles to speed things along
        "# pm3 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)",  # similar to existing #5, but uses tight rather than verytight; used for ADMPQLGIEMRGAT-UHFFFAOYAUmult3 (InChI=1/C6H14O5Si/c1-4-9-12(8,10-5-2)11-6(3)7/h6-7H,3-5H2,1-2H3/mult3)
        "# pm3 opt freq IOP(2/16=3)",  # use default (not verytight) convergence criteria; use this as last resort
        "# pm3 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)",  # to address problematic C10H14JJ case
        "# pm3 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm",  # for very troublesome RRMZRNPRCUANER-UHFFFAOYAQ (InChI=1/C5H7/c1-3-5-4-2/h3H,1-2H3) case...there were troubles with negative frequencies, where I don't think they should have been; step size of numerical frequency was adjusted to give positive result; accuracy of result is questionable; it is possible that not all of these keywords are needed; note that for this and other nearly free rotor cases, I think heat capacity will be overestimated by R/2 (R vs. R/2) (but this is a separate issue)
        "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",  #  for troublesome QDERTVAGQZYPHT-UHFFFAOYAHmult3(InChI=1/C6H14O4Si/c1-4-8-11(7,9-5-2)10-6-3/h4H,5-6H2,1-3H3/mult3); key aspects appear to be tight (rather than verytight) convergence criteria, no calculation of frequencies during optimization, use of numerical frequencies, and probably also the use of opt=small
        "# pm3 opt=(verytight,gdiis,calcall) IOP(2/16=3)",  # used for troublesome C5H7J case; note that before fixing, I got errors like the following: "Incomplete coordinate system.  Try restarting with Geom=Check Guess=Read Opt=(ReadFC,NewRedundant) Incomplete coordinate system. Error termination via Lnk1e in l103.exe"; we could try to restart, but it is probably preferrable to have each keyword combination standalone; another keyword that may be helpful if additional problematic cases are encountered is opt=small; 6/9/09 note: originally, this had # pm3 opt=(verytight,gdiis,calcall) freq IOP(2/16=3)" (with freq keyword), but I discovered that in this case, there are two thermochemistry sections and cclib parses frequencies twice, giving twice the number of desired frequencies and hence produces incorrect thermo; this turned up on C5H6JJ isomer
        "# pm3 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm",  # worked for troublesome ketene case: CCGKOQOJPYTBIH-UHFFFAOYAO (InChI=1/C2H2O/c1-2-3/h1H2) (could just increase number of iterations for similar keyword combination above (#6 at the time of this writing), allowing symmetry, but nosymm seemed to reduce # of iterations; I think one of nosymm or higher number of iterations would allow the similar keyword combination to converge; both are included here for robustness)
        "# pm3 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm",  # added for case of ZWMVZWMBTVHPBS-UHFFFAOYAEmult3 (InChI=1/C4H4O2/c1-3-5-6-4-2/h1-2H2/mult3)
        "# pm3 opt=(calcall,small,maxcyc=100) IOP(2/16=3)",  # used to address troublesome FILUFGAZMJGNEN-UHFFFAOYAImult3 case (InChI=1/C5H6/c1-3-5-4-2/h3H,1H2,2H3/mult3)
    ]


class GaussianMolPM6(GaussianMol):
    """
    Gaussian PM6 calculations for molecules

    This is a class of its own in case you wish to do anything differently,
    but for now it's only the 'pm6' in the keywords that differs.
    """

    #: Keywords that will be added at the top of the qm input file
    keywords = [
        # The combinations of keywords were derived by Greg Magoon for pm3. For now, we assume similar ones will work for pm6:
        "# pm6 opt=(verytight,gdiis) freq IOP(2/16=3)",
        "# pm6 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)",
        "# pm6 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm",
        "# pm6 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm",
        "# pm6 opt=(verytight,gdiis,small) freq IOP(2/16=3)",
        "# pm6 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)",
        "# pm6 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)",
        "# pm6 opt=tight freq IOP(2/16=3)",
        "# pm6 opt=tight freq=numerical IOP(2/16=3)",
        "# pm6 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)",
        "# pm6 opt freq IOP(2/16=3)",
        "# pm6 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)",
        "# pm6 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm",
        "# pm6 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
        "# pm6 opt=(verytight,gdiis,calcall) IOP(2/16=3)",
        "# pm6 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm",
        "# pm6 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm",
        "# pm6 opt=(calcall,small,maxcyc=100) IOP(2/16=3)",
    ]
