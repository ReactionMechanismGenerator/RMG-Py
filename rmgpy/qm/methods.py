
class QMMethod:
    max_attempts = 0
    
    def writeInputFile(self, geometry, attempt):
    
    def checkJob(self):
        # verify successful run
    
    def runJob(self):
        # send the job to be run
    
    def parseResult(self):
        # convert the output to useful parameters
        
    def processResult(self):
        # derive thermo properties

class MopacPM3(QMMethod):
    writeInputFile = qm.inputwriters.MOPACPM3InputWriter.write
    max_attempts = qm.inputwriters.MOPACPM3InputWriter.max_attempts
    
    def checkJob(self):
        # check that the QM job successfully ran
        
    def runJob(self):
        qm.jobs.MOPACJob.run
        
        return self.checkJob()
        
    def parseResult(self):
        qmData = qm.parsers.CCLibParser.parse
        
        return qmData
    
    def processResult(self):
        # package thermo
        
        return self.thermo