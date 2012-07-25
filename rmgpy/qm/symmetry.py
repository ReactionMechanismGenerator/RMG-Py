"""
Created on Apr 29, 2012
@author: nmvdewie

Module that collects all classes related to symmetry of molecules
"""

import os
import mopac

import logging
from subprocess import Popen, PIPE


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
    directory = 'QMfiles'

    RMG_path = os.environ.get("RMGpy")
    if RMG_path is None:
        RMG_path = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..'))
        logging.info("Setting RMG_path to {0}".format(RMG_path))
    executable_path = os.path.join(RMG_path, 'bin', 'symmetry')
    if not os.path.exists(executable_path):
        raise Exception("Symmetry program not found at {0}.".format(executable_path))

    'Argumenst that will be passed as an argument for the consecutive attempts'
    argumentsList = [
         ['-final', '0.02'],
         ['-final', '0.1'],
         ['-primary', '0.2', '-final' ,'0.1'],
         ['-final', '0.0'],
        ]

    inputFileExtension = '.symm'

    def __init__(self, molfile, iqmdata):
        self.molfile = molfile

        'the command line command'
        self.command = []

        "IQMData is the object that holds information from a previous QM Job on 3D coords, molecule etc..."
        self.qmdata = iqmdata

        self.inputFile = self.molfile + self.inputFileExtension

        self.attemptNumber = 1
        self.pointGroupFound = False

    def check(self, output):
        """
        Check the `output` string and extract the resulting point group, which is returned.
        """
        for line in output.split('\n'):
            if line.startswith("It seems to be the "): # "It seems to be the [x] point group" indicates point group.
                result = line.split(" ")[5]
                break
        else:
            logging.exception("Couldn't find point group from symmetry output:\n{0}".format(output))
            raise RuntimeError("Couldn't find point group in symmetry output.")

        logging.info("Point group: "+ result)
        return result;

    def run(self):
        """
        Run the command and wait for it to finish.
        """
        pp = Popen(self.command, stdout=PIPE, stderr=PIPE)
        stdout, stderr = pp.communicate()
        return self.check(stdout)    

    def writeInputFile(self):
        """
        Write the input file for the SYMMETRY program.
        """
        geom = str(self.qmdata.numberOfAtoms) + "\n"
        for i in range(self.qmdata.numberOfAtoms):
            geom = geom + " ".join((str(self.qmdata.atomicNumbers[i]),
                                    str(self.qmdata.atomCoords[i][0]),
                                    str(self.qmdata.atomCoords[i][1]),
                                    str(self.qmdata.atomCoords[i][2])
                                   )) + "\n"
        with open(os.path.join(self.directory, self.inputFile), 'w') as input_file:
            input_file.write(geom)
        input_file.close()
        logging.info("Symmetry input file written to %s"%os.path.join(self.directory, self.inputFile))
        return input_file

    def calculate(self):
        """
        Do the entire point group calculation.

        This writes the input file, then tries several times to run 'symmetry'
        with different parameters, until a point group is found and returned.
        """
        self.writeInputFile();

        result = "";

        #continue trying to generate symmetry group until too many no. of attempts or until a point group is found: 
        for attempt, arguments in enumerate(self.argumentsList):
            """
            TODO only *nix case works!
            """
            self.command = [self.executable_path]
            self.command.extend(arguments)
            self.command.append(os.path.join(self.directory, self.inputFile))

            #call the program and read the result
            result = self.run();

            #check for a recognized point group
            symmPGC = symmetry.PointGroupDictionary()

            if symmPGC.contains(result):
                self.pointGroupFound = True;
                return symmetry.PointGroup(result)
            else:
                logging.info("Attempt number {0} did not identify a recognized point group ({1}).".format(attempt,result))
        logging.critical("Final attempt did not identify a recognized point group. Exiting.")
        return None
        

class PointGroupDictionary:
    """
     Dictionary type with all the  symmetry groups supported in RMG
    """
    def __init__(self):
        #Map with the symmetry groups and their associated symmmetry numbers
        self.library = {}
        
        #Map with the symmetry groups and a flag whether it is chiral or not
        self.chiralLibrary = {}
        
        self.initiate()
    
    def populateGroups(self):
        """
        Could be externalized, so that users could add/remove symmetry groups
        """
        self.library["C1"] = 1
        self.library["Cs"] = 1
        self.library["Ci"] = 1
        self.library["C2"] = 2
        self.library["C3"] = 3
        self.library["C4"] = 4
        self.library["C5"] = 5
        self.library["C6"] = 6
        self.library["C7"] = 7
        self.library["C8"] = 8

        self.library["D2"] = 4
        self.library["D3"] = 6
        self.library["D4"] = 8
        self.library["D5"] = 10
        self.library["D6"] = 12
        self.library["D7"] = 14
        self.library["D8"] = 16

        self.library["C2v"] = 2
        self.library["C3v"] = 3
        self.library["C4v"] = 4
        self.library["C5v"] = 5
        self.library["C6v"] = 6
        self.library["C7v"] = 7
        self.library["C8v"] = 8

        self.library["C2h"] = 2
        self.library["C3h"] = 3
        self.library["C4h"] = 4
        self.library["C5h"] = 5
        self.library["C6h"] = 6
        self.library["C8h"] = 8

        self.library["D2h"] = 4
        self.library["D3h"] = 6
        self.library["D4h"] = 8
        self.library["D5h"] = 10
        self.library["D6h"] = 12
        self.library["D7h"] = 14
        self.library["D8h"] = 16

        self.library["D2d"] = 4
        self.library["D3d"] = 6
        self.library["D4d"] = 8
        self.library["D5d"] = 10
        self.library["D6d"] = 12
        self.library["D7d"] = 14
        self.library["D8d"] = 16

        self.library["S4"] = 2
        self.library["S6"] = 3
        self.library["S8"] = 4

        self.library["T"] = 12
        self.library["Th"] = 12
        self.library["Td"] = 12

        self.library["O"] = 24
        self.library["Oh"] = 24

        self.library["Cinfv"] = 1
        self.library["Dinfh"] = 2
        self.library["I"] = 60
        self.library["Ih"] = 60
        self.library["Kh"] = 1
        
    def populateChiralityFlags(self):
        """
        Could be externalized, so that users could add/remove symmetry groups
        """
        self.chiralLibrary["C1"] = True
        self.chiralLibrary["Cs"] = False
        self.chiralLibrary["Ci"] = False
        self.chiralLibrary["C2"] = True
        self.chiralLibrary["C3"] = True
        self.chiralLibrary["C4"] = True
        self.chiralLibrary["C5"] = True
        self.chiralLibrary["C6"] = True
        self.chiralLibrary["C7"] = True
        self.chiralLibrary["C8"] = True

        self.chiralLibrary["D2"] = True
        self.chiralLibrary["D3"] = True
        self.chiralLibrary["D4"] = True
        self.chiralLibrary["D5"] = True
        self.chiralLibrary["D6"] = True
        self.chiralLibrary["D7"] = True
        self.chiralLibrary["D8"] = True

        self.chiralLibrary["C2v"] = False
        self.chiralLibrary["C3v"] = False
        self.chiralLibrary["C4v"] = False
        self.chiralLibrary["C5v"] = False
        self.chiralLibrary["C6v"] = False
        self.chiralLibrary["C7v"] = False
        self.chiralLibrary["C8v"] = False

        self.chiralLibrary["C2h"] = False
        self.chiralLibrary["C3h"] = False
        self.chiralLibrary["C4h"] = False
        self.chiralLibrary["C5h"] = False
        self.chiralLibrary["C6h"] = False
        self.chiralLibrary["C8h"] = False

        self.chiralLibrary["D2h"] = False
        self.chiralLibrary["D3h"] = False
        self.chiralLibrary["D4h"] = False
        self.chiralLibrary["D5h"] = False
        self.chiralLibrary["D6h"] = False
        self.chiralLibrary["D7h"] = False
        self.chiralLibrary["D8h"] = False

        self.chiralLibrary["D2d"] = False
        self.chiralLibrary["D3d"] = False
        self.chiralLibrary["D4d"] = False
        self.chiralLibrary["D5d"] = False
        self.chiralLibrary["D6d"] = False
        self.chiralLibrary["D7d"] = False
        self.chiralLibrary["D8d"] = False

        self.chiralLibrary["S4"] = True
        self.chiralLibrary["S6"] = True
        self.chiralLibrary["S8"] = True

        self.chiralLibrary["T"] = True
        self.chiralLibrary["Th"] = False
        self.chiralLibrary["Td"] = False

        self.chiralLibrary["O"] = True
        self.chiralLibrary["Oh"] = False

        self.chiralLibrary["Cinfv"] = False
        self.chiralLibrary["Dinfh"] = False
        self.chiralLibrary["I"] = True
        self.chiralLibrary["Ih"] = False
        self.chiralLibrary["Kh"] = False
    
    def initiate(self):
        self.populateGroups()
        self.populateChiralityFlags()
        
    def contains(self,pointGroup):
        return pointGroup in self.library
    
    def get(self, pointGroup):
        symmNumber = self.library.get(pointGroup);

        if symmNumber:
            return symmNumber;
        else:
            return None;
    
    def isChiral(self, pointGroup):
        chiral = self.chiralLibrary.get(pointGroup);

        if chiral:
            return chiral;
        else:
            return None;
        
        

class PointGroupCalculator:
    """
    Wrapper type to determine molecular symmetry point groups based on 3D coords information.
     
    Will point to a specific algorithm, like SYMMETRY that is able to do this.
    """
    def __init__(self, molfile, iqmdata):
        self.molfile = molfile
        self.qmdata = iqmdata#data object that contains 3D coords of molecule used in symmetry calculation
        self.calculator = SymmetryJob(molfile, iqmdata)

    def calculate(self):
        return self.calculator.calculate();

class PointGroup:
    
    def __init__(self, pointGroup):
        self.pointGroup = pointGroup
        
        dic = PointGroupDictionary()
        
        self.symmetryNumber = dic.get(pointGroup)#integer
        self.chiral = dic.isChiral(pointGroup)#boolean
        
        #determine linearity from 3D-geometry; changed to correctly consider linear ketene radical case
        if self.pointGroup == "Cinfv" or self.pointGroup == "Dinfh":
            self.linear = True;
        else:
            self.linear = False;
        
    def equals(self,group):
        return self.pointGroup == group
    
    def isLinear(self):
        return self.linear;
    
