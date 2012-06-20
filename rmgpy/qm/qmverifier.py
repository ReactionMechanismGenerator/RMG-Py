"""
Created on Apr 29, 2012

@author: nmvdewie
"""

import os
import qmtp as qm
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
        
        self.failureKeys = ['IMAGINARY FREQUENCIES', 'EXCESS NUMBER OF OPTIMIZATION CYCLES', 'NOT ENOUGH TIME FOR ANOTHER CYCLE']
        
        self.successKeys = {}#maps whether a particular success keyword has been found in the output file.
        self.successKeys['DESCRIPTION OF VIBRATIONS'] = False
    
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
                    logging.info("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + self.molfile.name + " RMG Augmented InChI = "+ self.molfile.InChIAug + " Log file Augmented InChI = "+logFileInChI + " . InChI could not be found in the MOPAC input file. You should manually check that the output file contains the ended species.")
                    return False
        
    def successfulMopacResultExists(self):
        """
        Returns a boolean flag that states whether a successful MOPAC simulation already exists for the molecule with the 
        given (augmented) InChI Key.
        
        The definition of finding a successful simulation is based on these criteria:
        1) finding an output file with the file name equal to the InChI Key
        2) NOT finding any of the keywords that are denote a calculation failure
        3) finding all the keywords that denote a calculation success.
        4) finding a match between the InChI of the given molecule and the InchI found in the calculation files

        If any of the above criteria is not matched, False will be returned and the procedures to start a new calculation 
        will be initiated.
        
        """
        file = os.path.join(self.molfile.directory,self.molfile.name+self.outputExtension)
        
        if os.path.exists(file):#if the file exists, do further checks otherwise, we will skip to final statement and return False
            InChIMatch=False#flag (1 or 0) indicating whether the InChI in the file matches InChIaug this can only be 1 if InChIFound is also 1
            InChIFound=False#flag (1 or 0) indicating whether an InChI was found in the log file

            with open(file) as qmfile:    
                   for each_line in qmfile:
                       each_line = each_line.rstrip().strip()
                       
                       for element in self.failureKeys:#search for failure keywords
                           if element in each_line:
                               logging.error("MOPAC output file contains the following error %s")%element
                               return False
                           
                       for element in self.successKeys:#search for success keywords
                           if element in each_line:
                               self.successKeys[element] = True
                      
                       if "InChI=" in each_line:
                            logFileInChI = each_line#output files should take up to 240 characters of the name in the input file
                            InChIFound=True
                            if logFileInChI == self.molfile.InChIAug:
                                InChIMatch = True
                                
                   #check if ALL 'success' keywords were found in the file.    
                   success = True     
                   success = [success and t for t in self.successKeys.values()]
                   if not success:
                       logging.info('Not all of the required keywords for sucess were found in the output file!')
                       return False
            
            #if the failure flag is still 0, the process should have been successful
            if InChIMatch:
                logging.info("Pre-existing successful MOPAC quantum result for " + self.molfile.name + " ("+self.molfile.InChIAug+") has been found. This log file will be used.")
                return True

            elif InChIFound:#InChIs do not match (most likely due to limited name length mirrored in log file (240 characters), but possibly due to a collision)
                return self.checkForInChiKeyCollision(logFileInChI)          
            else:
                logging.critical("An InChI was not found in the MOPAC output file: " +self.molfile.name+".out nor in the input file.")
                return False  
        
        return False
    
        #returns True if an MM4 output file for the given name and directory (.mm4out suffix) exists and indicates successful completion (same criteria as used after calculation runs) terminates if the InChI doesn't match the InChI in the file or if there is no InChI in the file returns False otherwise
       
    def verify(self):
       
       self.mopacResultExists = self.successfulMopacResultExists()
        
    def verifyNoFailure(self):
        """
        checks whether the output file contains any of the 
        failure keywords
        """
        file = os.path.join(self.molfile.directory,self.molfile.name+self.outputExtension)
        with open(file) as qmfile:    
                   for each_line in qmfile:
                       each_line = each_line.rstrip().strip()
                       for element in self.failureKeys:#search for failure keywords
                           if element in each_line:
                               logging.error("MOPAC output file contains the following error %s")%element
                               return False
                           
        return True
    
    def succesfulJobExists(self):
        """
        checks whether one of the flags is true.
        If so, it returns true.
        """
     return self.gaussianResultExists or self.mopacResultExists or self.mm4ResultExists
