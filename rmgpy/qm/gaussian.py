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
    # edit for gaussian
    #executablePath = os.path.join(os.getenv('MOPAC_DIR', default="/opt/mopac") , 'MOPAC2009.exe')

    usePolar = False#use polar keyword in MOPAC

    "Keywords that will be added at the top of the qm input file"
    # these are just for pm3. need to edit and put a keyword editor in specific gausspm3 class
    keywordsTop = {}#keywords that will be added to the qm input file based on the attempt number
    keywordsTop[1] = "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3)"
    keywordsTop[2] = "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)"
    keywordsTop[3] = "# pm3 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm" 
    keywordsTop[4] = "# pm3 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm"
    keywordsTop[5] = "# pm3 opt=(verytight,gdiis,small) freq IOP(2/16=3)"
    keywordsTop[6] = "# pm3 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)"
    keywordsTop[7] = "# pm3 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)"
    keywordsTop[8] = "# pm3 opt=tight freq IOP(2/16=3)"
    keywordsTop[9] = "# pm3 opt=tight freq=numerical IOP(2/16=3)"
    keywordsTop[10] = "# pm3 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)"
    keywordsTop[11] = "# pm3 opt freq IOP(2/16=3)"
    keywordsTop[12] = "# pm3 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)"
    keywordsTop[13] = "# pm3 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm"
    keywordsTop[14] = "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm"
    keywordsTop[15] = "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm"
    keywordsTop[16] = "# pm3 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm"
    keywordsTop[17] = "# pm3 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm"
    keywordsTop[18] = "# pm3 opt=(calcall,small,maxcyc=100) IOP(2/16=3)"
    
    scriptAttempts = len(keywordsTop)
    maxAttempts = 2 * scriptAttempts

    failureKeys = ['IMAGINARY FREQUENCIES', 'EXCESS NUMBER OF OPTIMIZATION CYCLES', 'NOT ENOUGH TIME FOR ANOTHER CYCLE']

    def writeInputFile(self, attempt, top_keys, bottom_keys, polar_keys):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        keywords = "\n".join((
            "%chk={0}/RMGrunCHKfile.chk".format(self.molfile.directory),
            "%mem=6MW",
            "%nproc=1",
            self.keywordsTop[self.attemptNumber] + ' polar' if qmtp.QMTP.usePolar else ''
            ))
        inputFilePath = os.path.join(self.directory,self.geometry.uniqueID + self.inputFileExtension)
        
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "gjf")
        mol = openbabel.OBMol()
        
        if self.attemptNumber <= self.scriptAttempts: #use UFF-refined coordinates
            obConversion.ReadFile(mol, self.geometry.getRefinedMolFilePath())
        else:
            obConversion.ReadFile(mol, self.geometry.getCrudeMolFilePath())
        
        mol.SetTitle(self.geometry.uniqueID)
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
        
        input_string = obConversion.WriteString(mol)
        
        with open(inputFilePath, 'w') as gjfFile:
            gjfFile.write(keywords)
            gjfFile.write(input_string)
        
        return self.geometry.uniqueID + self.inputFileExtension
        
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
                        logging.error("MOPAC output file contains the following error %s")%element
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