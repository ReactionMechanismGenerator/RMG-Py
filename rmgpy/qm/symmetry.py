"""
Created on Apr 29, 2012
@author: nmvdewie and rwest

Module that collects all classes related to symmetry of molecules
"""

import os

import logging
from subprocess import Popen, PIPE

class PointGroup:
    """
    A symmetry Point Group.
    
    Attributes are:
    
     * pointGroup
     * symmetryNumber
     * chiral
     * linear 
    """
    def __init__(self, pointGroup, symmetryNumber, chiral):

        self.pointGroup = pointGroup
        self.symmetryNumber = symmetryNumber
        self.chiral = chiral

        #determine linearity from 3D-geometry; changed to correctly consider linear ketene radical case
        if self.pointGroup in ["Cinfv", "Dinfh"]:
            self.linear = True;
        else:
            self.linear = False;

    def __repr__(self):
        return 'PointGroup("{0}", symmetryNumber={1}, chiral={2})'.format(self.pointGroup, self.symmetryNumber, self.chiral)

def makePointGroupDictionary():
    """
    A function to make and fill the point group dictionary.
    
    This will be stored once as a module level variable.
    """
    pointGroupDictionary = dict()
    for (pointGroup, symmetryNumber, chiral) in [
        ("C1"   ,  1,  True ),
        ("Cs"   ,  1,  False),
        ("Ci"   ,  1,  False),
        ("C2"   ,  2,  True ),
        ("C3"   ,  3,  True ),
        ("C4"   ,  4,  True ),
        ("C5"   ,  5,  True ),
        ("C6"   ,  6,  True ),
        ("C7"   ,  7,  True ),
        ("C8"   ,  8,  True ),
    
        ("D2"   ,  4,  True ),
        ("D3"   ,  6,  True ),
        ("D4"   ,  8,  True ),
        ("D5"   , 10,  True ),
        ("D6"   , 12,  True ),
        ("D7"   , 14,  True ),
        ("D8"   , 16,  True ),
    
        ("C2v"  ,  2,  False),
        ("C3v"  ,  3,  False),
        ("C4v"  ,  4,  False),
        ("C5v"  ,  5,  False),
        ("C6v"  ,  6,  False),
        ("C7v"  ,  7,  False),
        ("C8v"  ,  8,  False),
    
        ("C2h"  ,  2,  False),
        ("C3h"  ,  3,  False),
        ("C4h"  ,  4,  False),
        ("C5h"  ,  5,  False),
        ("C6h"  ,  6,  False),
        ("C8h"  ,  8,  False),
    
        ("D2h"  ,  4,  False),
        ("D3h"  ,  6,  False),
        ("D4h"  ,  8,  False),
        ("D5h"  , 10,  False),
        ("D6h"  , 12,  False),
        ("D7h"  , 14,  False),
        ("D8h"  , 16,  False),
    
        ("D2d"  ,  4,  False),
        ("D3d"  ,  6,  False),
        ("D4d"  ,  8,  False),
        ("D5d"  , 10,  False),
        ("D6d"  , 12,  False),
        ("D7d"  , 14,  False),
        ("D8d"  , 16,  False),
    
        ("S4"   ,  2,  True ),
        ("S6"   ,  3,  True ),
        ("S8"   ,  4,  True ),
    
        ("T"    , 12,  True ),
        ("Th"   , 12,  False),
        ("Td"   , 12,  False),
    
        ("O"    , 24,  True ),
        ("Oh"   , 24,  False),
    
        ("Cinfv",  1,  False),
        ("Dinfh",  2,  False),
        ("I"    , 60,  True ),
        ("Ih"   , 60,  False),
        ("Kh"   ,  1,  False),
    ]:
        pointGroupDictionary[pointGroup] = PointGroup(pointGroup, symmetryNumber, chiral)
    return pointGroupDictionary
    
#: A dictionary of PointGroup objects, stored as a module level variable.
pointGroupDictionary = makePointGroupDictionary()


class PointGroupCalculator:
    """
    Wrapper type to determine molecular symmetry point groups based on 3D coords information.

    Will point to a specific algorithm, like SYMMETRY that is able to do this.
    """
    def __init__(self, settings, uniqueID, qmData):
        self.uniqueID = uniqueID
        self.qmData = qmData # QMDdata object that contains 3D coords of molecule used in symmetry calculation
        self.calculator = SymmetryJob(settings, uniqueID, qmData)

    def calculate(self):
        return self.calculator.calculate();


class SymmetryJob:
    """
    Determine the point group using the SYMMETRY program 
    
    (http://www.cobalt.chem.ucalgary.ca/ps/symmetry/).

    Required input is a line with number of atoms followed by lines for each atom 
    including:
    1) atom number
    2) x,y,z coordinates

    finalTol determines how loose the point group criteria are;
    values are comparable to those specified in the GaussView point group interface

    """

    'Arguments that will be passed as an argument for the consecutive attempts'
    argumentsList = [
         ['-final', '0.02'],
         ['-final', '0.1'],
         ['-primary', '0.2', '-final' ,'0.1'],
         ['-final', '0.0'],
        ]

    inputFileExtension = '.symm'

    def __init__(self, settings, uniqueID, qmData):
        self.settings = settings
        self.uniqueID = uniqueID

        "The object that holds information from a previous QM Job on 3D coords, molecule etc..."
        self.qmData = qmData
        self.attemptNumber = 1
        self.pointGroupFound = False

        if os.sys.platform == 'win32':
            self.executable_path = os.path.join(settings.RMG_bin_path, 'symmetry.exe')
        else:
            self.executable_path = os.path.join(settings.RMG_bin_path, 'symmetry')
        if not os.path.exists(self.executable_path):
            raise Exception("Symmetry program not found at {0}.".format(self.executable_path))
        

    @property
    def inputFilePath(self):
        "The input file's path"
        return os.path.join(self.settings.fileStore, self.uniqueID + self.inputFileExtension)
        
    def parse(self, output):
        """
        Check the `output` string and extract the resulting point group, which is returned.
        """
        for line in output.split('\n'):
            if line.startswith("It seems to be the "): # "It seems to be the [x] point group" indicates point group.
                result = line.split(" ")[5]
                break
        else:
            logging.exception("Couldn't find point group from symmetry output:\n{0}".format(output))
            return "Not found"

        logging.info("Point group: "+ result)
        return result;

    def run(self, command):
        """
        Run the command, wait for it to finish, and return the stdout.
        """
        pp = Popen(command, stdout=PIPE, stderr=PIPE)
        stdout, stderr = pp.communicate()
        if stderr:
            logging.error("Error message from SYMMETRY calculation:")
            logging.error(stderr)
        return stdout

    def writeInputFile(self):
        """
        Write the input file for the SYMMETRY program.
        """
        geom = str(self.qmData.numberOfAtoms) + "\n"
        coords_in_angstrom = self.qmData.atomCoords.value_si * 1e10
        for i in range(self.qmData.numberOfAtoms):
            geom = geom + " ".join((str(self.qmData.atomicNumbers[i]),
                                    str(coords_in_angstrom[i][0]),
                                    str(coords_in_angstrom[i][1]),
                                    str(coords_in_angstrom[i][2])
                                   )) + "\n"
        with open(self.inputFilePath, 'w') as input_file:
            input_file.write(geom)
        input_file.close()
        logging.info("Symmetry input file written to {0}".format(self.inputFilePath))
        return input_file

    def calculate(self):
        """
        Do the entire point group calculation.

        This writes the input file, then tries several times to run 'symmetry'
        with different parameters, until a point group is found and returned.
        """
        self.writeInputFile();

        #continue trying to generate symmetry group until too many no. of attempts or until a point group is found: 
        for attempt, arguments in enumerate(self.argumentsList):
            """
            TODO only *nix case works!
            """
            command = [self.executable_path]
            command.extend(arguments)
            command.append(self.inputFilePath)

            #call the program and read the result
            output = self.run(command)
            # parse the output to get a point group name
            pointGroupName = self.parse(output)
            
            if pointGroupDictionary.has_key(pointGroupName):
                self.pointGroupFound = True;
                return pointGroupDictionary[pointGroupName]
            else:
                logging.info("Attempt number {0} did not identify a recognized point group ({1}).".format(attempt,pointGroupName))
                if attempt+2 == len(self.argumentsList):
                    logging.warning("Using last-resort symmetry estimation options; symmetry may be underestimated.")
        logging.critical("Final attempt did not identify a recognized point group. Exiting.")
        return None


