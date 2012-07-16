from mopacmol import MopacMol

class MopacPM3(MopacMol):
    def __init__():
        # what it needs to be initialized
    
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
        MopacMol.run()
            
    def check():
        # might have more or less keywords to search for
    
    def parse():
        # might have different parsing
    
    def generateQMThermoData():
        # call the methods to generate the thermodata
        