import os

import logging

from molecule import QMMolecule, TDPropertiesCalculator
from gaussian import Gaussian

class G09Mol(QMMolecule, Gaussian):
    
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
    
    def writeInputFile(self, attempt, keywordsTop):
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
    
    def generateQMThermoData(self):
        # call the methods to generate the thermodata
        self.createGeometry()

        success = False
        method = G09PM3(self)
        for attempt in range(1, self.maxAttempts+1):
            top_keys= method.inputFileKeys(attempt)
            inputFileName = self.writeInputFile(attempt, top_keys)
            success = self.run(inputFileName)
            if success:
                logging.info('Attempt {0} of {1} on species {2} succeeded.'.format(attempt, self.maxAttempts, self.molecule.toAugmentedInChI()))
                result = self.parse() # parsed in cclib
                break
        else:
            raise Exception('QM thermo calculation failed for {0}.'.format(InChIaug))

        thermo = TDPropertiesCalculator(result, self.getInChiAug())

        return thermo.calculate()



class G09PM3(G09Mol):
    def inputFileKeys(self, attempt):
        """
        Inherits the writeInputFile methods from mopac.py
        """
        multiplicity_keys = self.multiplicityKeywords[self.molecule.geometry.multiplicity]

        top_keys = "pm3 {0} {1}".format(
                multiplicity_keys,
                self.keywordsTop[attempt],
                )
        return top_keys