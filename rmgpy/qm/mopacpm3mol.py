import os

from mopacmol import MopacMol

class MopacPM3(MopacMol):
    def __init__(self, molecule):
        self.molecule = molecule
        
    def createInputFile():
        """
        Inherits the writeInputFile methods from mopac.py
        """
        top_keywords = "pm3 {0} {1}".format(
                multiplicity_keyword,
                self.keywordsTop[attempt],
                )
        bottom_keywords = "{0} pm3 {1}".format(
                self.keywordsBottom[attempt],
                multiplicity_keyword,
                )
        polar_keywords = "oldgeo {0} nosym precise pm3 {1}".format(
                'polar' if geometry.multiplicity == 1 else 'static',
                multiplicity_keyword,
                )
        
        self.writeInputFile(top_keywords, bottom_keywords, polar_keywords, geometry, attempt)
        
    def run():
        # Send the file to run
        self.MopacMol.run()
            
    def check():
        # might have more or less keywords to search for
        pass
    
    def parse():
        # might have different parsing
        pass
    
    def generateQMThermoData(self):
        # call the methods to generate the thermodata
        self.createGeometry()
        
        import ipdb; ipdb.set_trace()
        writer = method.InputWriter()
        verifier = method.Verifier()
        
        success = False
        for attempt in range(method.max_attempts):
            writer.writeInputFile(self.geometry, attempt)
            success = method.runJob()
            method.parseResult()
            if success: break
        else:
            raise Exception("Couldn't generate thermo data.")
        thermoData = method.processResult()
        thermo = MopacMol.generateQMThermoData(molecule)
        
        return thermo