#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
Module that collects all classes related to symmetry of molecules
"""

import logging
import os
from subprocess import Popen, PIPE


class PointGroup(object):
    """
    A symmetry Point Group.

    Attributes are:

     * point_group
     * symmetry_number
     * chiral
     * linear
    """

    def __init__(self, point_group, symmetry_number, chiral):
        self.point_group = point_group
        self.symmetry_number = symmetry_number
        self.chiral = chiral

        # determine linearity from 3D-geometry; changed to correctly consider linear ketene radical case
        if self.point_group in ["Cinfv", "Dinfh"]:
            self.linear = True
        else:
            self.linear = False

    def __repr__(self):
        return 'PointGroup("{0}", symmetry_number={1}, chiral={2})'.format(
            self.point_group, self.symmetry_number, self.chiral
        )


def _make_point_group_dictionary():
    """
    A function to make and fill the point group dictionary.

    This will be stored once as a module level variable.
    """
    point_group_dictionary = dict()
    for point_group, symmetry_number, chiral in [
        ("C1", 1, True),
        ("Cs", 1, False),
        ("Ci", 1, False),
        ("C2", 2, True),
        ("C3", 3, True),
        ("C4", 4, True),
        ("C5", 5, True),
        ("C6", 6, True),
        ("C7", 7, True),
        ("C8", 8, True),
        ("D2", 4, True),
        ("D3", 6, True),
        ("D4", 8, True),
        ("D5", 10, True),
        ("D6", 12, True),
        ("D7", 14, True),
        ("D8", 16, True),
        ("C2v", 2, False),
        ("C3v", 3, False),
        ("C4v", 4, False),
        ("C5v", 5, False),
        ("C6v", 6, False),
        ("C7v", 7, False),
        ("C8v", 8, False),
        ("C2h", 2, False),
        ("C3h", 3, False),
        ("C4h", 4, False),
        ("C5h", 5, False),
        ("C6h", 6, False),
        ("C8h", 8, False),
        ("D2h", 4, False),
        ("D3h", 6, False),
        ("D4h", 8, False),
        ("D5h", 10, False),
        ("D6h", 12, False),
        ("D7h", 14, False),
        ("D8h", 16, False),
        ("D2d", 4, False),
        ("D3d", 6, False),
        ("D4d", 8, False),
        ("D5d", 10, False),
        ("D6d", 12, False),
        ("D7d", 14, False),
        ("D8d", 16, False),
        ("S4", 2, True),
        ("S6", 3, True),
        ("S8", 4, True),
        ("T", 12, True),
        ("Th", 12, False),
        ("Td", 12, False),
        ("O", 24, True),
        ("Oh", 24, False),
        ("Cinfv", 1, False),
        ("Dinfh", 2, False),
        ("I", 60, True),
        ("Ih", 60, False),
        ("Kh", 1, False),
    ]:
        point_group_dictionary[point_group] = PointGroup(
            point_group, symmetry_number, chiral
        )
    return point_group_dictionary


#: A dictionary of PointGroup objects, stored as a module level variable.
POINT_GROUP_DICTIONARY = _make_point_group_dictionary()


class PointGroupCalculator(object):
    """
    Wrapper type to determine molecular symmetry point groups based on 3D coords information.

    Will point to a specific algorithm, like SYMMETRY that is able to do this.
    """

    def __init__(self, settings, unique_id, qm_data):
        self.unique_id = unique_id
        self.qm_data = qm_data  # QMDdata object that contains 3D coords of molecule used in symmetry calculation
        self.calculator = SymmetryJob(settings, unique_id, qm_data)

    def calculate(self):
        return self.calculator.calculate()


class SymmetryJob(object):
    """
    Determine the point group using the SYMMETRY program

    (Originally ``http://www.cobalt.chem.ucalgary.ca/ps/symmetry/``
     now mirrored at https://github.com/nquesada/symmetry).

    Required input is a line with number of atoms followed by lines for each atom
    including:
    1) atom number
    2) x,y,z coordinates

    finalTol determines how loose the point group criteria are;
    values are comparable to those specified in the GaussView point group interface

    """

    "Arguments that will be passed as an argument for the consecutive attempts"
    argumentsList = [
        ["-final", "0.02"],
        ["-final", "0.1"],
        ["-primary", "0.2", "-final", "0.1"],
        ["-final", "0.0"],
    ]

    input_file_extension = ".symm"

    def __init__(self, settings, unique_id, qm_data):
        self.settings = settings
        self.unique_id = unique_id

        "The object that holds information from a previous QM Job on 3D coords, molecule etc..."
        self.qm_data = qm_data
        self.attempt_number = 1
        self.point_group_found = False

    @property
    def input_file_path(self):
        """The input file's path"""
        return os.path.join(
            self.settings.scratchDirectory, self.unique_id + self.input_file_extension
        )

    def parse(self, output):
        """
        Check the `output` string and extract the resulting point group, which is returned.
        """
        for line in output.split("\n"):
            if line.startswith(
                "It seems to be the "
            ):  # "It seems to be the [x] point group" indicates point group.
                result = line.split(" ")[5]
                break
        else:
            logging.exception(
                "Couldn't find point group from symmetry output:\n{0}".format(output)
            )
            return "Not found"

        logging.info("Point group: " + result)
        return result

    def run(self, command):
        """
        Run the command, wait for it to finish, and return the stdout.
        """
        try:
            pp = Popen(command, stdout=PIPE, stderr=PIPE)
        except OSError as e:
            logging.error(e)
            raise Exception(
                "Running symmetry on the point group calculation has failed. Please check if symmetry "
                "program is installed on your system in RMG-Py/bin/symmetry or on your path."
            )
        stdout, stderr = pp.communicate()
        if stderr:
            logging.error("Error message from SYMMETRY calculation:")
            logging.error(stderr)
        return stdout.decode("utf-8")

    def write_input_file(self):
        """
        Write the input file for the SYMMETRY program.
        """
        geom = str(self.qm_data.numberOfAtoms) + "\n"
        coords_in_angstrom = self.qm_data.atomCoords.value_si * 1e10
        for i in range(self.qm_data.numberOfAtoms):
            geom = (
                geom
                + " ".join(
                    (
                        str(self.qm_data.atomicNumbers[i]),
                        str(coords_in_angstrom[i][0]),
                        str(coords_in_angstrom[i][1]),
                        str(coords_in_angstrom[i][2]),
                    )
                )
                + "\n"
            )
        with open(self.input_file_path, "w") as input_file:
            input_file.write(geom)
        input_file.close()
        logging.info("Symmetry input file written to {0}".format(self.input_file_path))
        return input_file

    def calculate(self):
        """
        Do the entire point group calculation.

        This writes the input file, then tries several times to run 'symmetry'
        with different parameters, until a point group is found and returned.
        """
        self.write_input_file()

        # continue trying to generate symmetry group until too many no. of attempts or until a point group is found:
        for attempt, arguments in enumerate(self.argumentsList):
            """
            TODO only *nix case works!
            """
            command = [self.settings.symmetryPath]
            command.extend(arguments)
            command.append(self.input_file_path)

            # call the program and read the result
            output = self.run(command)
            # parse the output to get a point group name
            point_group_name = self.parse(output)

            if point_group_name in POINT_GROUP_DICTIONARY:
                self.point_group_found = True
                return POINT_GROUP_DICTIONARY[point_group_name]
            else:
                logging.info(
                    "Attempt number {0} did not identify a recognized "
                    "point group ({1}).".format(attempt, point_group_name)
                )
                if attempt + 2 == len(self.argumentsList):
                    logging.warning(
                        "Using last-resort symmetry estimation options; symmetry may be underestimated."
                    )
        logging.critical(
            "Final attempt did not identify a recognized point group. Exiting."
        )
        return None
