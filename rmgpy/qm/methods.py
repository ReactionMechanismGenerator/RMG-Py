
class QMMethod:
    max_attempts = 0
    
    def writeInputFile(self, geometry, attempt):
    
    def checkJob(self):
        # verify successful run
        # return 'False' by default
        return False
    
    def runJob(self):
        # send the job to be run
    
    def parseResult(self):
        # convert the output to useful parameters
        
    def processResult(self):
        # derive thermo properties

class MopacPM3(QMMethod):
    writeInputFile = qm.inputwriters.MOPACPM3InputWriter.write
    max_attempts = qm.inputwriters.MOPACPM3InputWriter.max_attempts
    
    def checkJob():
        # Check the job was completed
        
    def runJob(self):
        # Run the job. Returns 'True' if the job successfully ran, 'False' if not.
        job = qm.jobs.MOPACJob(self)
        
        return self.checkJob()
        
    def parseResult(self):
        # cclib, Quixote?
        qmData = qm.parsers.CCLibParser.parse
        
        return qmData
    
    def processResult(self):
        # package thermo, calculator
        
        return self.thermo

class Gaussian03PM3(QMMethod):
    writeInputFile = qm.inputwriters.GaussianPM3InputWriter.write
    max_attempts = qm.inputwriters.GaussianPM3InputWriter.max_attempts
    
    def checkJob():
        # Check the job was completed
        
    def runJob(self):
        # Run the job. Returns 'True' if the job successfully ran, 'False' if not.
        # not current 'GaussianJob' method
        job = qm.jobs.GaussianJob(self)
        
        return self.checkJob()
        
    def parseResult(self):
        # cclib, Quixote?
        qmData = qm.parsers.CCLibParser.parse
        
        return qmData
    
    def processResult(self):
        # package thermo, calculator
        
        return self.thermo