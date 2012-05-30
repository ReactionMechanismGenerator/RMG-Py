''''
Created on Apr 29, 2012

@author: nmvdewie
'''

import os
import qmtp as qm
import logging

class QMVerifier:
    '''
     * Verifies whether a QM job (externalized) was succesfully completed by 
 * searching for specific keywords in the output files, 
 * located in a specific directory (e.g. "QMFiles")
    '''
    def __init__(self,name,InChIaug, directory, QMTP):
        self.name = name
        self.InChIaug = InChIaug
        self.directory = directory
        self.QMTP = QMTP
        self.gaussianResultExists = False
        self.mopacResultExists = False
        self.mm4ResultExists = False
    
    def successfulMopacResultExistsQ(self):
        '''
        #returns True if a MOPAC output file for the given name and directory (.out suffix) exists and indicates
        successful completion (same criteria as used after calculation runs)
        terminates if the InChI doesn't match the InChI in the file or if there is no InChI in the file returns False otherwise
                                                    
        
        #look in the output file to check for the successful termination of the calculation 
        (assumed to be successful if "description of vibrations appears)
        '''
        file = os.path.join(self.directory,self.name+".out")
        
        if os.path.exists(file):#if the file exists, do further checks otherwise, we will skip to final statement and return False
            failureFlag=True#flag (1 or 0) indicating whether the MOPAC job failed
            failureOverrideFlag=False#flag (1 or 0) to override success as measured by failureFlag
            InChIMatch=False#flag (1 or 0) indicating whether the InChI in the file matches InChIaug this can only be 1 if InChIFound is also 1
            InChIFound=False#flag (1 or 0) indicating whether an InChI was found in the log file
            InChIPartialMatch=False#flag (1 or 0) indicating whether the InChI in the log file is a substring of the InChI in RMG's memory
            logFileInChI=""
            try:
               with open(file) as qmfile:
                   failure_expr = []
                       #negative frequencies notice example:
                        #         NOTE: SYSTEM IS NOT A GROUND STATE, THEREFORE ZERO POINT
                        #         ENERGY IS NOT MEANINGFULL. ZERO POINT ENERGY PRINTED
                        #         DOES NOT INCLUDE THE  2 IMAGINARY FREQUENCIES
                   failure_expr.append("IMAGINARY FREQUENCIES")
                   failure_expr.append("EXCESS NUMBER OF OPTIMIZATION CYCLES")#exceeding max cycles error
                   failure_expr.append("NOT ENOUGH TIME FOR ANOTHER CYCLE")#timeout error
                       
                   for each_line in qmfile:
                       each_line = each_line.rstrip().strip()
                       
                       if "DESCRIPTION OF VIBRATIONS" in each_line:
                            # if !MopacFileContainsNegativeFreqsQ(name, directory)) failureFlag=0#check for this line if it is here, check for negative frequencies
                            failureFlag = False
                        
                       elif "InChI=" in each_line:
                            logFileInChI = each_line#output files should take up to 240 characters of the name in the input file
                            InChIFound=1
                            if logFileInChI == self.InChIaug:
                                InChIMatch=1
                            elif self.InChIaug.startswith(logFileInChI):
                                InChIPartialMatch=1
                                 
                       for element in failure_expr:
                           if element in each_line:
                               failureOverrideFlag=True
                               failureFlag=True #job will be considered a failure if there are imaginary frequencies or if job terminates to to excess time/cycles
                               logging.error("MOPAC output file contains the following error %s")%element
                               break;
                       
                        
                       
            
            except Exception as e:
                err = "Error in reading preexisting MOPAC output file \n"
                err = err+str(e)
                logging.error(err)
               
            
            #if the failure flag is still 0, the process should have been successful
            if not failureFlag and InChIMatch:
                logging.info("Pre-existing successful MOPAC quantum result for " + self.name + " ("+self.InChIaug+") has been found. This log file will be used.")
                return True
            
            elif InChIFound and not InChIMatch:#InChIs do not match (most likely due to limited name length mirrored in log file (240 characters), but possibly due to a collision)
                # if InChIPartialMatch == 1:#case where the InChI in memory begins with the InChI in the log file we will continue and check the input file, pring a warning if there is no match
                #look in the input file if the InChI doesn't match (apparently, certain characters can be deleted in MOPAC output file for long InChIs)
                inputFile = os.path.join(self.directory,self.name+".mop")
                
                if os.path.exists(inputFile):#read the MOPAC inputFile
                    inputFileInChI=""
                    try:
                       with open(inputFile) as inputFile:
                            lineI = inputFile.readline()

                            while not lineI:
                                if lineI.startswith("InChI="):
                                    inputFileInChI = lineI.rstrip()
                                
                                lineI = inputFile.readline()

                    
                    except Exception as e:
                        err = "Error in reading preexisting MOPAC input file \n"
                        err = err+str(e)
                        logging.error(err)
                       
                    
                    if inputFileInChI == self.InChIaug:
                        if not failureFlag:
                            logging.info("Pre-existing successful MOPAC quantum result for " + self.name + " ("+self.InChIaug+") has been found. This log file will be used. *Note that input file was read to confirm lack of InChIKey collision (InChI probably more than 240 characters or characters probably deleted from InChI in .out file)")
                            return True
                        
                        else:#otherwise, failureFlag==1
                            logging.info("Pre-existing MOPAC quantum result for " + self.name + " ("+self.InChIaug+") has been found, but the result was apparently unsuccessful. The file will be overwritten with a new calculation or Gaussian result (if available) will be used. *Note that input file was read to confirm lack of InChIKey collision (InChI probably more than 240 characters or characters probably deleted from InChI in .out file)")
                            return False
                        
                    
                    else:
                        if inputFileInChI == "":#InChI was not found in input file
                            logging.info("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + self.name + " RMG Augmented InChI = "+ self.InChIaug + " Log file Augmented InChI = "+logFileInChI + " . InChI could not be found in the MOPAC input file. You should manually check that the output file contains the ended species.")
                            return True
                        
                        else:#InChI was found but doesn't match
                            logging.critical("Congratulations! You appear to have discovered the first recorded instance of an InChIKey collision: InChIKey(augmented) = " + self.name + " RMG Augmented InChI = "+ self.InChIaug + " MOPAC input file Augmented InChI = " + inputFileInChI + " Log file Augmented InChI = "+logFileInChI)
                           
                        
                    
                
                else:
                    logging.info("*****Warning: potential InChIKey collision: InChIKey(augmented) = " + self.name + " RMG Augmented InChI = "+ self.InChIaug + " Log file Augmented InChI = "+logFileInChI + " . MOPAC input file could not be found to check full InChI. You should manually check that the log file contains the ended species.")
                    return True
                
                #  
                #                else:
                #                    logging.critical("Congratulations! You appear to have discovered the first recorded instance of an InChIKey collision: InChIKey(augmented) = " + self.name + " RMG Augmented InChI = "+ self.InChIaug + " MOPAC output file Augmented InChI = "+logFileInChI)
                #                   
                #                
            
            elif not InChIFound:
                logging.critical("An InChI was not found in file: " +self.name+".out")
               
            
            elif failureFlag:#note these should cover all possible results for this block, and if the file.exists block is entered, it should return from within the block and should not reach the return statement below
                logging.info("Pre-existing MOPAC quantum result for " + self.name + " ("+self.InChIaug+") has been found, but the result was apparently unsuccessful. The file will be overwritten with a new calculation or Gaussian result (if available) will be used.")
                return False
            
        
        #we could pr a line here for cases where the file doesn't exist, but this would probably be too verbose
        return False
    
        #returns True if an MM4 output file for the given name and directory (.mm4out suffix) exists and indicates successful completion (same criteria as used after calculation runs) terminates if the InChI doesn't match the InChI in the file or if there is no InChI in the file returns False otherwise
       
    def verify(self):
        '''
        /**
     * Depending on the QM method used, a specific verification method will be called,
     * and assigns the correct flag stored as an attribute.
     */
        ''' 
        if self.QMTP.qmMethod == "pm3":
            #first, check to see if the result already exists and the job terminated successfully
            self.mopacResultExists = self.successfulMopacResultExistsQ()
        
       
    def succesfulJobExists(self):
     '''
      /**
     * checks whether one of the flags is true.
     * If so, it returns true.
     * @return
     */
     '''
     return self.gaussianResultExists or self.mopacResultExists or self.mm4ResultExists
     
    