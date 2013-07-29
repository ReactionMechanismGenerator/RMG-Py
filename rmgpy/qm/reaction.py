import os
import logging
import numpy

from rmgpy.molecule import Molecule
from rmgpy.species import TransitionState
from molecule import QMMolecule, Geometry
from collections import defaultdict

try:
    import rdkit
    from rdkit import DistanceGeometry
    from rdkit.Chem.Pharm3D import EmbedLib
except ImportError:
    logging.info("To use transition state searches, you must correctly install rdkit")

class QMReaction:
    
    file_store_path = 'QMfiles'
    if not os.path.exists(file_store_path):
        logging.info("Creating directory %s for mol files."%os.path.abspath(file_store_path))
        os.makedirs(file_store_path)
    
    def __init__(self, reaction, settings):
        self.reaction = reaction
        self.settings = settings
        
        self.geometry = None
        self.transitionState = None
        
    def getFilePath(self, extension):
        """
        Should return the path to the file with the given extension.
        
        The provided extension should include the leading dot.
        
        Need to define some reaction line notation.
        Possibly '<Reaction_Family>/<reactant1SMILES>+<reactant2SMILES>--<product1SMILES>+<product2SMILES>' ???
        """
        raise NotImplementedError("This should be a unique string representing the reaction")
        return os.path.join(self.settings.fileStore, self.uniqueID  + extension)
        
    def fixSortLabel(self, molecule):
        """
        This may not be required anymore. Was needed as when molecules were created, the
        rmg sorting labels would be set after where we tried to generate the TS.
        """
        sortLbl = 0
        for vertex in molecule.vertices:
            vertex.sortingLabel = sortLbl
            sortLbl += 1
        return molecule
    
    def getGeometry(self, molecule, settings):
        
        multiplicity = sum([i.radicalElectrons for i in molecule.atoms]) + 1
        geom = Geometry(settings, molecule.toAugmentedInChIKey(), molecule, multiplicity)
        
        return geom, multiplicity
        
    def getRDKitMol(self, geometry):
        """
        Check there is no RDKit mol file already made. If so, use rdkit to make a rdmol from
        a mol file. If not, make rdmol from geometry.
        """ 
        if not os.path.exists(geometry.getCrudeMolFilePath()):
            geometry.generateRDKitGeometries()
        rdKitMol = rdkit.Chem.MolFromMolFile(geometry.getCrudeMolFilePath(), removeHs=False)      
        
        return rdKitMol
        
    def generateBoundsMatrix(self, molecule, settings):
        """
        Uses rdkit to generate the bounds matrix of a rdkit molecule.
        """
        self.geometry, multiplicity = getGeometry(molecule, settings)
        rdKitMol = getRDKitMol(self.geometry)
        boundsMatrix = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(rdKitMol)
        
        return rdKitMol, boundsMatrix, multiplicity
        
    def generateGeometry(self):
        quantumMechanics = QMCalculator()
        quantumMechanics.settings.software = 'gaussian'
        quantumMechanics.settings.fileStore = 'QMfiles'
        quantumMechanics.settings.scratchDirectory = 'scratch'
        quantumMechanics.settings.onlyCyclics = False
        quantumMechanics.settings.maxRadicalNumber = 0
        
        
        notes = ''
        
        reactant = fixSortLabel(TS[0])
        product = fixSortLabel(TS[1])
        
        rRDMol, tsBM, tsMult = generateBoundsMatrix(reactant, quantumMechanics.settings)
        
        # edit bounds distances to align reacting atoms
        if family.lower() == 'h_abstraction':
            sect = len(reactant.split()[1].atoms)
        
            tsBM[sect:,:sect] = 1.8
        
            lbl1 = reactant.getLabeledAtom('*1').sortingLabel
            lbl2 = reactant.getLabeledAtom('*2').sortingLabel
            lbl3 = reactant.getLabeledAtom('*3').sortingLabel
            
            labels = [lbl1, lbl2, lbl3]
            atomMatch = ((lbl1,),(lbl2,),(lbl3,))
            
            reaction = Reaction(reactants=reactant.split(), products=product.split())
            distanceData = transitionStates.estimateDistances(reaction)
            
            tsBM = editMatrix(tsBM, lbl1, lbl2, distanceData.distances['d12'], 0.1)
            tsBM = editMatrix(tsBM, lbl2, lbl3, distanceData.distances['d23'], 0.001)
            tsBM = editMatrix(tsBM, lbl1, lbl3, distanceData.distances['d13'], 0.001)
        
        setBM = rdkit.DistanceGeometry.DoTriangleSmoothing(tsBM)
        
        if setBM:
            tssorted_atom_list = reactant.vertices[:]
            qmcalcTS = rmgpy.qm.gaussian.GaussianMolPM3(reactant, quantumMechanics.settings)
            reactant.vertices = tssorted_atom_list
        
            qmcalcTS.createGeometry(tsBM,atomMatch)
        
            geometryTS = qmcalcTS.geometry
            tsmolFilePathForCalc = qmcalcTS.getMolFilePathForCalculation(attempt)
        
            # Create filenames
            tsFilePath = os.path.join(quantumMechanics.settings.fileStore, 'transitionState' + str(count) + '.gjf')
            tsOutPath = os.path.join(quantumMechanics.settings.fileStore, 'transitionState' + str(count) + '.log')
            ircInPath = os.path.join(quantumMechanics.settings.fileStore, 'transitionStateIRC' + str(count) + '.gjf')
            ircOutPath = os.path.join(quantumMechanics.settings.fileStore, 'transitionStateIRC' + str(count) + '.log')
            tsOutputDataFile = os.path.join(quantumMechanics.settings.fileStore, 'data' + str(count) + outputFileExtension)
        
            # QM saddle search
            # Write and run the TS optimization
            tsConverge = 0
            if os.path.exists(tsOutPath):
                tsConverge = checkOutput(tsOutPath)
            else:
                writeTSInputFile(tsFilePath, tsmolFilePathForCalc, geometryTS, family)
                run(executablePath, tsFilePath, tsOutPath)
                tsConverge = checkOutput(tsOutPath)
                # writeSubFile(tsFilePath.split('.')[0])
                # os.system('bsub < ' + tsFilePath.split('.')[0] + '.sh')
        
            # Validation
            if tsConverge == 2:
                # Error in internal coodinate system, continue calculation in cartesian
                writeTSCartInput(tsFilePath, count)
                run(executablePath, tsFilePath, tsOutPath)
                tsConverge = checkOutput(tsOutPath)
                # os.system('bsub < ' + tsFilePath.split('.')[0] + '.sh')
            
            # If saddle found, write and run the IRC calculation to check the TS geometry
            if tsConverge == 1:
                writeIRCInput(ircInPath, count)
                run(executablePath, ircInPath, ircOutPath)
                # writeSubFile(ircInPath.split('.')[0])
                # os.system('bsub < ' + ircInPath.split('.')[0] + '.sh')
                ircCheck, notes = testGeometries(reactant, product, ircOutPath, notes)
                if ircCheck == 1:
                    vibFreq, activeAts, atomDist = parse(tsOutPath, tsOutputDataFile, labels)
                    writeRxnOutputFile(tsOutputDataFile, reactant, product, vibFreq, activeAts, atomDist, notes)
                    