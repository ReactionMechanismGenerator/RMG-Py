class Geometry:
    file_store_path = 'QMfiles'
    
    def __init__( uniqueID, rmg_molecule ):
        self.uniqueID = uniqueID
        self.rmg_molecule = rmg_molecule
    
    def generateRDKitGeometries(boundsMatrix=None):
        """
        Use RDKit to guess geometry.
        
        Save mol files of both crude and refined.
        Saves coordinates on atoms.
        """
        rdmol = self.rd_embed(boundsMatrix)
        rdmol = self.rd_optimize(rdmol, boundsMatrix)
        self.save_coordinates(rdmol)

    def rd_embed(boundsMatrix=None):
        "Create rdmol from rmg_molecule"
        "Embed in 3D using rdkit"
        "Save crude geom in mol file, named  uniqueID.crude.mol"
        return rdmol

    def rd_optimize(rdmol, boundsMatrix=None):
        "Optimize using UFF in RDKit"
        "Save optimized geom in mol file, named after uniqueID.refined.mol"
        
        return rdmol
    
    def save_coordinates(rdmol):
        "Save xyz coordinates on each atom in rmg_molecule"


class QM_Molecule:
    def __init__( rmg_molecule ):
        self.rmg_molecule = rmg_molecule
        
    def createGeometry(self):
        """
        Creates self.geometry with RDKit geometries
        """
        uniqueID = self.generateInChiAug()
        self.geometry = Geometry( uniqueID, rmg_molecule )
        self.geometry.generateRDKitGeometries()
        

    def generateThermoData(self, option):
        self.createGeometry()
        
        method = qm.methods.MopacPM3 # for example
        
        for attempt in range(method.max_attempts):
            method.writeInputFile(self.geometry, attempt)
            method.runJob()
            method.parseResult()
            success = method.validate()
            if success: break
        else:
            raise Exception("Couldn't generate thermo data")
        thermoData = method.processResult()
        return thermoData
        
    def getInChiAug(self):
        """
        Returns the augmented InChI from self.rmg_molecule 
        """
        return augmented_InChI
        
        
class QM_Reaction:
    
    def createGeometry(self):
        """
        Creates self.geometry with RDKit geometries
        """
        
        "create merged rmg_molecule"
        "create boundsMatrix"
        "create uniqueID for reaction"
        self.geometry = Geometry( uniqueID, rmg_molecule )
        self.geometry.generateRDKitGeometries(boundsMatrix)
    
    def generateKineticData():
        self.createGeometry()
        "write input file"
        "run job"
        "parse result"
        "process result"
        return kineticData
        
## This goes in methods:

class QMMethod:
    max_attempts = 0
    def writeInputFile(self, geometry, attempt):
    def runJob(self):
    def parseResult(self):


class MopacPM3(QMMethod):
    writeInputFile = qm.inputwriters.MOPACPM3InputWriter.write
    max_attempts = qm.inputwriters.MOPACPM3InputWriter.max_attempts
    
