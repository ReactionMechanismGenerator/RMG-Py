"""
Created on Apr 29, 2012

@author: nmvdewie
"""

import os

import logging

class QMVerifier:
    """
    Verifies whether a QM job (externalized) was succesfully completed by 
      * searching for specific keywords in the output files, 
      * located in a specific directory (e.g. "QMFiles")
    """
    def __init__(self,molfile):
        self.molfile = molfile
        self.gaussianResultExists = False
        self.mopacResultExists = False
        self.mm4ResultExists = False
        
        self.outputExtension = '.out'
        self.inputExtension = '.mop'
    
    def checkForInChiKeyCollision(self,logFileInChI):
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
        #look in the input file if the InChI doesn't match (apparently, certain characters can be deleted in MOPAC output file for long InChIs)
        inputFile = os.path.join(self.molfile.directory,self.molfile.name+self.inputExtension)
            
        assert os.path.exists(inputFile)

        with open(inputFile) as inputFile:#read the MOPAC inputFile
            lineI = inputFile.readline()
            for line in inputFile:
                if line.startswith("InChI="):
                    inputFileInChI = lineI.rstrip()
                    break
            
            if inputFileInChI == self.molfile.InChIAug:
                logging.info('An InChI match could be found in the input file, but not in the output file. Anywho, a\
                pre-existing successful MOPAC result exists.')
                return True
            
            else:
                logging.info("*****Warning: potential InChIKey collision: InChIKey(augmented) = " +
                             self.molfile.name + " RMG Augmented InChI = "+ self.molfile.InChIAug +
                             " Log file Augmented InChI = "+logFileInChI + 
                             " . InChI could not be found in the MOPAC input file. You should manually check that the output file contains the ended species.")
                return False
    

        #returns True if an MM4 output file for the given name and directory (.mm4out suffix) exists and indicates successful completion (same criteria as used after calculation runs) terminates if the InChI doesn't match the InChI in the file or if there is no InChI in the file returns False otherwise

    def succesfulJobExists(self):
        """
        checks whether one of the flags is true.
        If so, it returns true.
        """
        return self.gaussianResultExists or self.mopacResultExists or self.mm4ResultExists
