import os
import logging
import external.cclib.parser
import time
from subprocess import Popen
from copy import deepcopy
import numpy

from rmgpy.molecule import Molecule
from rmgpy.species import Species, TransitionState
from rmgpy.kinetics import Wigner
from molecule import QMMolecule, Geometry
from rmgpy.cantherm.gaussian import GaussianLog
from rmgpy.cantherm.kinetics import KineticsJob
from rmgpy.data.kinetics.transitionstates import TransitionStates
import symmetry

try:
    import rdkit
    from rdkit import DistanceGeometry
    from rdkit.Chem.Pharm3D import EmbedLib
except ImportError:
    logging.info("To use transition state searches, you must correctly install rdkit")

transitionStates = TransitionStates()
transitionStates.load(os.path.join(os.getenv('HOME'), 'Code/RMG-database/input/kinetics/families/H_Abstraction'), None, None)

class QMReaction:
    
    file_store_path = 'QMfiles'
    if not os.path.exists(file_store_path):
        logging.info("Creating directory %s for mol files."%os.path.abspath(file_store_path))
        os.makedirs(file_store_path)
    
    def __init__(self, reaction, settings):
        self.reaction = reaction
        self.settings = settings
        
        reactants = sorted([s.toSMILES() for s in self.reaction.reactants])
        products = sorted([s.toSMILES() for s in self.reaction.products])
        stringID = "+".join(reactants) + "_" + "+".join(products)
        
        self.uniqueID = stringID
        
        self.geometry = None
        self.transitionState = None
    
    def getFilePath(self, extension):
        """
        Should return the path to the file with the given extension.
        
        The provided extension should include the leading dot.
        
        Need to define some reaction line notation.
        Possibly '<Reaction_Family>/<reactant1SMILES>+<reactant2SMILES>--<product1SMILES>+<product2SMILES>' ???
        """
        return os.path.join(self.settings.fileStore, self.uniqueID  + extension)
    
    @property
    def outputFilePath(self):
        """Get the output file name."""
        return self.getFilePath(self.outputFileExtension)
    
    @property
    def inputFilePath(self):
        """Get the input file name."""
        return self.getFilePath(self.inputFileExtension)
    
    @property
    def ircOutputFilePath(self):
        """Get the irc output file name."""
        return self.getFilePath('IRC' + self.outputFileExtension)
    
    @property
    def ircInputFilePath(self):
        """Get the irc input file name."""
        return self.getFilePath('IRC' + self.inputFileExtension)
        
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
        
    def generateBoundsMatrix(self, molecule):
        """
        Uses rdkit to generate the bounds matrix of a rdkit molecule.
        """
        self.geometry, multiplicity = self.getGeometry(molecule, self.settings)
        rdKitMol = self.getRDKitMol(self.geometry)
        boundsMatrix = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(rdKitMol)
        
        return rdKitMol, boundsMatrix, multiplicity
    
    def setLimits(self, bm, lbl1, lbl2, value, uncertainty):
        if lbl1 > lbl2:
            bm[lbl2][lbl1] = value + uncertainty//2
            bm[lbl1][lbl2] = max(0,value - uncertainty//2)
        else:
            bm[lbl2][lbl1] = max(0,value - uncertainty//2)
            bm[lbl1][lbl2] = value + uncertainty//2
    
        return bm
    
    def editMatrix(self, reactant, bm):
        
        """
        For bimolecular reactions, reduce the minimum distance between atoms
        of the two reactants. 
        """
        if self.reaction.label.lower() == 'h_abstraction':
            
            lbl1 = reactant.getLabeledAtom('*1').sortingLabel
            lbl2 = reactant.getLabeledAtom('*2').sortingLabel
            lbl3 = reactant.getLabeledAtom('*3').sortingLabel
        
        elif self.reaction.label.lower() == 'disproportionation':
            
            lbl1 = reactant.getLabeledAtom('*2').sortingLabel
            lbl2 = reactant.getLabeledAtom('*4').sortingLabel
            lbl3 = reactant.getLabeledAtom('*1').sortingLabel
            
        labels = [lbl1, lbl2, lbl3]
        atomMatch = ((lbl1,),(lbl2,),(lbl3,))
        
        distanceData = transitionStates.estimateDistances(self.reaction)
        
        sect = len(reactant.split()[1].atoms)
        
        uncertainties = distanceData.uncertainties or {'d12':0.1, 'd13':0.1, 'd23':0.1 } # default if uncertainty is None
        
        bm = self.setLimits(bm, lbl1, lbl2, distanceData.distances['d12'], uncertainties['d12'])
        bm = self.setLimits(bm, lbl2, lbl3, distanceData.distances['d23'], uncertainties['d23'])
        bm = self.setLimits(bm, lbl1, lbl3, distanceData.distances['d13'], uncertainties['d13'])
        
        for i in range(sect,len(bm)):
            for j in range(0,sect):
                for k in range(len(bm)):
                    if k==i or k==j: continue
                    Uik = bm[i,k] if k>i else bm[k,i]
                    Ukj = bm[j,k] if k>j else bm[k,j]
                    
                    maxLij = Uik + Ukj - 0.15
                    if bm[i,j] >  maxLij:
                        print "CHANGING {0} to {1}".format(bm[i,j], maxLij)
                        bm[i,j] = maxLij
            
        return bm, labels, atomMatch
        
    def generateTSGeometry(self):
        """
        
        """
        if os.path.exists(os.path.join(self.file_store_path, self.uniqueID + '.data')):
            return True
        else:
            if len(self.reaction.reactants)==2:
                reactant = self.reaction.reactants[0].merge(self.reaction.reactants[1])
            if len(self.reaction.products)==2:
                product = self.reaction.products[0].merge(self.reaction.products[1])
            
            reactant = self.fixSortLabel(reactant)
            product = self.fixSortLabel(product)
            tsRDMol, tsBM, tsMult = self.generateBoundsMatrix(reactant)
            
            self.geometry.uniqueID = self.uniqueID
            
            tsBM, labels, atomMatch = self.editMatrix(reactant, tsBM)
            atoms = len(reactant.atoms)
            distGeomAttempts = 15*(atoms-3) # number of conformers embedded from the bounds matrix
            
            setBM = rdkit.DistanceGeometry.DoTriangleSmoothing(tsBM)
            
            if setBM:
                for i in range(len(tsBM)):
                    for j in range(i,len(tsBM)):
                        if tsBM[j,i] > tsBM[i,j]:
                                print "BOUNDS MATRIX FLAWED {0}>{1}".format(tsBM[j,i], tsBM[i,j])
            
                self.geometry.rd_embed(tsRDMol, distGeomAttempts, bm=tsBM, match=atomMatch)
                
                if not os.path.exists(self.outputFilePath):
                    self.writeInputFile(1)
                    converged, internalCoord = self.run()
                else:
                    converged, internalCoord = self.verifyOutputFile()
                
                if internalCoord and not converged:
                    self.writeInputFile(2)
                    converged = self.run()
                
                if converged:
                    if not os.path.exists(self.ircOutputFilePath):
                        self.writeIRCFile()
                        rightTS = self.runIRC()
                    else:
                        rightTS = self.verifyIRCOutputFile()
                    if rightTS:
                        self.writeRxnOutputFile(labels)
                        return True
                    else:
                        return False
                else:
                    return False
    
    def calculateQMData(self, moleculeList):
        """
        If the transition state is found, optimize reactant and product geometries for use in
        TST calculations.
        """
        molecules = []
        
        for molecule in moleculeList:
            molecule = self.fixSortLabel(molecule)
            qmMolecule = self.getQMMolecule(molecule)
            result = qmMolecule.generateQMData()
            if result:
                log = GaussianLog(qmMolecule.outputFilePath)
                species = Species(label=qmMolecule.molecule.toSMILES(), conformer=log.loadConformer(), molecule=[molecule])
                molecules.append(species)
        return molecules
    
    def calculateKinetics(self):
        # provides transitionstate geometry
        tsFound = self.generateTSGeometry()
        
        if tsFound:
            reactants = self.calculateQMData(self.reaction.reactants)
            products = self.calculateQMData(self.reaction.products)
            
            if len(reactants)==len(self.reaction.reactants) and len(products)==len(self.reaction.products):
                #self.determinePointGroup()
                tsLog = GaussianLog(self.outputFilePath)
                self.reaction.transitionState = TransitionState(label=self.uniqueID + 'TS', conformer=tsLog.loadConformer(), frequency=(tsLog.loadNegativeFrequency(), 'cm^-1'), tunneling=Wigner(frequency=None))
                                
                self.reaction.reactants = reactants
                self.reaction.products = products

                
                kineticsJob = KineticsJob(self.reaction)
                kineticsJob.generateKinetics()
                
                """
                What do I do with it? For now just save it.
                Various parameters are not considered in the calculations so far e.g. symmetry.
                This is just a crude calculation, calculating the partition functions
                from the molecular properties and plugging them through the equation. 
                """
                kineticsJob.save(self.getFilePath('.kinetics'))
                # return self.reaction.kinetics
        else:
            # fall back on group additivity
            return None
    
    def determinePointGroup(self):
        """
        Determine point group using the SYMMETRY Program
        
        Stores the resulting :class:`PointGroup` in self.pointGroup
        """
        assert self.qmData, "Need QM Data first in order to calculate point group."
        pgc = symmetry.PointGroupCalculator(self.settings, self.uniqueID, self.qmData)
        self.pointGroup = pgc.calculate()
