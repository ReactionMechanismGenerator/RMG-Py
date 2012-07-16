class Geometry:
    file_store_path = 'QMfiles'
    if not os.path.exists(file_store_path):
        logging.info("Creating directory %s for mol files."%os.path.abspath(file_store_path))
        os.makedirs(file_store_path)
    
    def __init__(self, uniqueID, rmg_molecule ):
        self.uniqueID = uniqueID
        self.rmg_molecule = rmg_molecule
        self.multiplicity = # get it from rmg_molecule
        
    def getCrudeMolFilePath(self):
        # os.join, uniqueID, suffix
        return self.file_store_path + self.uniqueID + '.crude.mol'
        
    def getRefinedMolFilePath(self):
        "Returns the path the the refined mol file"
        return self.file_store_path + self.uniqueID + '.refined.mol'
        
    def generateRDKitGeometries(boundsMatrix=None):
        """
        Use RDKit to guess geometry.
        
        Save mol files of both crude and refined.
        Saves coordinates on atoms.
        """
        rdmol = self.rd_build()
        rdmol, rdAtIdx = self.rd_embed(rdmol, boundsMatrix)
        rdmol = self.rd_optimize(rdmol, boundsMatrix)
        self.save_coordinates(rdmol, rdAtIdx)

    def rd_build():
        """
        Import rmg molecule and create rdkit molecule with the same atom labeling.
        """
        rdAtomIdx = {} # dictionary of rdkit atom indices
        # Initialize a blank Editable molecule and add all the atoms from RMG molecule
        rdmol = AllChem.rdchem.EditableMol(AllChem.rdchem.Mol())
        for index, atom in enumerate(self.molecule.vertices):
            rdAtom = AllChem.rdchem.Atom(atom.element.symbol)
            rdAtom.SetNumRadicalElectrons(atom.radicalElectrons)
            rdmol.AddAtom(rdAtom)
            rdAtomIdx[atom] = index
        
        # Add the bonds
        for atom1 in self.molecule.vertices:
            for atom2, bond in atom1.edges.items():
                index1 = rdAtomIdx[atom1] # atom1.sortingLabel
                index2 = rdAtomIdx[atom2] # atom2.sortingLabel
                if index1 > index2:
                    # Check the RMG bond order and add the appropriate rdkit bond.
                    if bond.order == 'S':
                        rdBond = AllChem.rdchem.BondType.SINGLE
                    elif bond.order == 'D':
                        rdBond = AllChem.rdchem.BondType.DOUBLE
                    elif bond.order == 'T':
                        rdBond = AllChem.rdchem.BondType.TRIPLE
                    elif bond.order == 'B':
                        rdBond = AllChem.rdchem.BondType.AROMATIC
                    else:
                        print "Unknown bond order"
                    rdmol.AddBond(index1, index2, rdBond)
        
        # Make editable mol into a mol and rectify the molecule
        rdmol = rdmol.GetMol()
        Chem.SanitizeMol(rdmol)
        
        return rdmol
        
    def rd_embed(rdmol, boundsMatrix=None):
        """
        Embed the RDKit molecule and create the crude molecule file.
        """
        AllChem.EmbedMultipleConfs(rdMol, numConfAttempts,randomSeed=1)
        crude = Chem.Mol(rdMol.ToBinary()) 
        
        with open(uniqueID + '.crude.mol', 'w') as out3Dcrude:
            out3Dcrude.write(Chem.MolToMolBlock(crude,confId=minEid))
            
        return rdmol, rdAtIdx
        
    def rd_optimize(rdmol, boundsMatrix=None):
        """
        Optimize the RDKit molecule and create the refined molecule file from this
        UFF-optimized molecular structure.
        """
        energy=0.0
        minEid=0;
        lowestE=9.999999e99;#start with a very high number, which would never be reached
        
        for i in range(rdmol.GetNumConformers()):
            AllChem.UFFOptimizeMolecule(rdmol,confId=i)
            energy=AllChem.UFFGetMoleculeForceField(rdmol,confId=i).CalcEnergy()
            if energy < lowestE:
                minEid = i
                lowestE = energy
        with open(uniqueID + '.refined.mol', 'w') as out3D:
            out3D.write(Chem.MolToMolBlock(rdmol,confId=minEid))
        
        return rdmol
    
    def save_coordinates(rdmol, rdAtIdx):
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
        
        if option.program == 'mopac':
            method = qm.methods.mopac_pm3 # for example
        elif option.program == 'gaussian03':
            method = qm.methods.Gaussian03PM3
        else:
            raise Exception('Unknown QM Method')
            
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