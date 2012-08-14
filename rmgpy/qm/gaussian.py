import os

import openbabel
import cclib.parser
import logging
from subprocess import Popen, PIPE

from qmdata import CCLibData

class Gaussian:

    directory = 'QMfiles'
    inputFileExtension = '.gif'
    outputFileExtension = ''
    lavaQueueFileExtension = '.lav'
    # edit for gaussian
    #executablePath = os.path.join(os.getenv('MOPAC_DIR', default="/opt/mopac") , 'MOPAC2009.exe')

    usePolar = False
    
    scriptAttempts = len(keywordsTop)
    maxAttempts = 2 * scriptAttempts

    failureKeys = ['IMAGINARY FREQUENCIES', 'EXCESS NUMBER OF OPTIMIZATION CYCLES', 'NOT ENOUGH TIME FOR ANOTHER CYCLE']

    def writeLavaSub(self, inputFileName):
        """
        Write the submission script for gaussian files to a cluster using the
        lava queueing system.
        """
        inputFilePath = os.path.join(directory + inputFileName)
        with open(inputFilePath, 'w') as lava:
            lava.write('#!/bin/sh\n')
            lava.write('#BSUB -q normal\n')
            lava.write('#BSUB -o ./output/' + title + '.out\n')
            lava.write('#BSUB -J reactionTS\n\n')
            lava.write('export GAUSS_EXEDIR=/share/apps/g09\n')
            lava.write('export PATH=$GAUSS_EXEDIR:$PATH\n\n')
            lava.write('g09 < ' + filename + ' > ./log/' + title + '.log\n\n')
            lava.close()
        
    def run(self):
        # submits the input file to gaussian
        command = os.path.join(self.directory, self.geometry.uniqueID + self.inputFileExtension)
        process = Popen([self.executablePath, command])
        process.communicate()# necessary to wait for executable termination!

        return self.checkNoFailure()

    def checkNoFailure(self):
        """
        checks whether the output file contains any of the 
        failure keywords
        """
        file = os.path.join(self.directory,self.geometry.uniqueID+self.outputFileExtension)
        with open(file) as qmfile:    
            for each_line in qmfile:
                each_line = each_line.rstrip().strip()
                for element in self.failureKeys:#search for failure keywords
                    if element in each_line:
                        logging.error("G09 output file contains the following error %s")%element
                        return False

        return True

    def read():
        # reads the output file
        import ipdb; ipdb.set_trace()
        pass

    def calculate():
        # calculators for the parsing
        import ipdb; ipdb.set_trace()
        pass

    def parse(self):
        # parses the output file to generate the TPs
        path = os.path.join(self.directory,self.geometry.uniqueID+'.out')
        try:
            myfile = cclib.parser.ccopen(path)
            myfile.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information

            cclibData = myfile.parse()
            radicalNumber = sum([i.radicalElectrons for i in self.molecule.atoms])

            qmData = CCLibData(cclibData, radicalNumber+1)

        except Exception as e:
            logging.error('Error in reading/parsing ccLib Python process.')
            logging.error(str(e))
            raise

        return qmData