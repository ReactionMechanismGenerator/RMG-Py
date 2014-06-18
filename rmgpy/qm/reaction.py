import os
import re
import logging
import external.cclib.parser
import time
from subprocess import Popen
from copy import deepcopy
import numpy
import shutil
from rmgpy.molecule import Molecule
from rmgpy.species import Species, TransitionState
from rmgpy.kinetics import Wigner
import mopac
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
# transitionStates.load(os.path.join(os.getenv('HOME'), 'Code/RMG-database/input/kinetics/families/intra_H_migration'), None, None)
# transitionStates.load(os.path.join(os.getenv('HOME'), 'Code/RMG-database/input/kinetics/families/R_Addition_MultipleBond'), None, None)

def matrixToString(matrix):
    """Returns a string representation of a matrix, for printing to the console"""
    text = '\n'.join([ ' '.join([str(round(item, 1)) for item in line]) for line in matrix ])
    return text.replace('1000.0', '1e3')

                
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
        geometry.generateRDKitGeometries()
        rdKitMol = rdkit.Chem.MolFromMolFile(geometry.getCrudeMolFilePath(), removeHs=False)      
        
        return rdKitMol
        
    def generateBoundsMatrix(self, molecule):
        """
        Uses rdkit to generate the bounds matrix of a rdkit molecule.
        """
        geometry, multiplicity = self.getGeometry(molecule, self.settings)
        rdKitMol = self.getRDKitMol(geometry)
        boundsMatrix = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(rdKitMol)
        
        return rdKitMol, boundsMatrix, multiplicity, geometry
    
    def setLimits(self, bm, lbl1, lbl2, value, uncertainty):
        if lbl1 > lbl2:
            bm[lbl2][lbl1] = value + uncertainty/2
            bm[lbl1][lbl2] = max(0,value - uncertainty/2)
        else:
            bm[lbl2][lbl1] = max(0,value - uncertainty/2)
            bm[lbl1][lbl2] = value + uncertainty/2
    
        return bm
    
    def bmPreEdit(self, bm, sect):
        """
        Clean up some of the atom distance limits before attempting triangle smoothing.
        This ensures any edits made do not lead to unsolvable scenarios for the molecular
        embedding algorithm.
        
        sect is the list of atom indices belonging to one species.
        """
        others = range(len(bm))
        for idx in sect: others.remove(idx)
            
        for i in range(len(bm)):#sect:
            for j in range(i):#others:
                if i<j: continue
                for k in range(len(bm)):
                    if k==i or k==j or i==j: continue
                    Uik = bm[i,k] if k>i else bm[k,i]
                    Ukj = bm[j,k] if k>j else bm[k,j]
                    
                    maxLij = Uik + Ukj - 0.1
                    if bm[i,j] >  maxLij:
                        print "CHANGING Lower limit {0} to {1}".format(bm[i,j], maxLij)
                        bm[i,j] = maxLij
        
        return bm
    
    def bmTest(self, bm):
        """
        Test to show why a bounds matrix is not valid
        """
        for k in range(len(bm)):
            for i in range(len(bm)-1):
                if i==k: continue
                Uik = bm[i,k] if k>i else bm[k,i]
                Lik = bm[i,k] if i>k else bm[k,i]
                for j in range(i+1, len(bm)):
                    if j==k: continue
                    Ujk = bm[j,k] if k>j else bm[k,j]
                    Ljk = bm[j,k] if j>k else bm[j,k]
                    Uij = bm[i,j] if j>i else bm[j,i]
                    Lij = bm[i,j] if i>j else bm[j,i]
                    sumUikUjk = Uik + Ujk
                    if Uij > sumUikUjk:
                        print "Upper limit for {i} and {j} is too high".format(i=i, j=j)
                    
                    diffLikUjk = Lik - Ujk
                    diffLjkUik = Ljk - Uik
                    if Uij < diffLikUjk or Uik < diffLjkUik:
                        print "Lower limit for {i} and {j} is too low".format(i=i, j=j)
    
    def editDoubMatrix(self, reactant, product, bm1, bm2):
        """
        For bimolecular reactions, reduce the minimum distance between atoms
        of the two reactanting species, in preparation for a double-ended search.
        `bm1` is typically the bounds matrix for the reactant side, and `bm2` for
        the products. 
        """
        def fixMatrix(bm, lbl1, lbl2, lbl3, num, diff):
            if lbl2 > lbl1:
                dnDiff = bm[lbl1][lbl2]
                upDiff = bm[lbl2][lbl1]
            else:
                dnDiff = bm[lbl2][lbl1]
                upDiff = bm[lbl1][lbl2]
            
            if lbl1 > lbl3:
                bm[lbl3][lbl1] = num + diff/2.
                bm[lbl1][lbl3] = num - diff/2.
                if lbl2 > lbl3:
                    bm[lbl3][lbl2] = bm[lbl3][lbl1] - upDiff
                    bm[lbl2][lbl3] = bm[lbl1][lbl3] - dnDiff
                else:
                    bm[lbl2][lbl3] = bm[lbl3][lbl1] - upDiff
                    bm[lbl3][lbl2] = bm[lbl1][lbl3] - dnDiff
            else:
                bm[lbl1][lbl3] = num + diff/2.
                bm[lbl3][lbl1] = num - diff/2.
                if lbl2 > lbl3:
                    bm[lbl3][lbl2] = bm[lbl1][lbl3] - upDiff
                    bm[lbl2][lbl3] = bm[lbl3][lbl1] - dnDiff
                else:
                    bm[lbl2][lbl3] = bm[lbl1][lbl3] - upDiff
                    bm[lbl3][lbl2] = bm[lbl3][lbl1] - dnDiff
            return bm
                
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
            
        if (reactant.atoms[lbl1].symbol == 'H' and reactant.atoms[lbl3].symbol == 'C') or (reactant.atoms[lbl1].symbol == 'C' and reactant.atoms[lbl3].symbol == 'H'):
            bm1 = fixMatrix(bm1, lbl1, lbl2, lbl3, 2.3, 0.1)
            bm2 = fixMatrix(bm2, lbl3, lbl2, lbl1, 2.3, 0.1)
        # elif (reactant.atoms[lbl1].symbol == 'H' and reactant.atoms[lbl3].symbol == 'O') or (reactant.atoms[lbl1].symbol == 'O' and reactant.atoms[lbl3].symbol == 'H'):
        #     bm1 = fixMatrix(bm1, lbl1, lbl2, lbl3, 2.2, 0.1)
        #     bm2 = fixMatrix(bm2, lbl3, lbl2, lbl1, 2.2, 0.1)
        # elif reactant.atoms[lbl1].symbol == 'O' and reactant.atoms[lbl3].symbol == 'O':
        #     bm1 = fixMatrix(bm1, lbl1, lbl2, lbl3, 2.2, 0.1)
        #     bm2 = fixMatrix(bm2, lbl3, lbl2, lbl1, 2.2, 0.1)
        else:
            bm1 = fixMatrix(bm1, lbl1, lbl2, lbl3, 2.7, 0.1)
            bm2 = fixMatrix(bm2, lbl3, lbl2, lbl1, 2.7, 0.1)
        
        # sect = len(reactant.split()[1].atoms)
        rSect = []
        for atom in reactant.split()[0].atoms: rSect.append(atom.sortingLabel)
        
        pSect = []
        for atom in product.split()[0].atoms: pSect.append(atom.sortingLabel)
            
        bm1 = self.bmPreEdit(bm1, rSect)
        bm2 = self.bmPreEdit(bm2, pSect)
        
        return bm1, bm2, labels, atomMatch
    
    def editMatrix(self, reactant, bm):
        
        """
        For bimolecular reactions, reduce the minimum distance between atoms
        of the two reactants. 
        """
        if self.reaction.label.lower() in ['h_abstraction', 'r_addition_multiplebond', 'intra_h_migration']:
            
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
        
        sect = []
        for atom in reactant.split()[0].atoms: sect.append(atom.sortingLabel)
        
        uncertainties = {'d12':0.1, 'd13':0.1, 'd23':0.1 }#distanceData.uncertainties or {'d12':0.1, 'd13':0.1, 'd23':0.1 } # default if uncertainty is None
        bm = self.setLimits(bm, lbl1, lbl2, distanceData.distances['d12'], uncertainties['d12'])
        bm = self.setLimits(bm, lbl2, lbl3, distanceData.distances['d23'], uncertainties['d23'])
        bm = self.setLimits(bm, lbl1, lbl3, distanceData.distances['d13'], uncertainties['d13'])
        
        bm = self.bmPreEdit(bm, sect)
            
        return bm, labels, atomMatch
        
    def generateTSGeometryDoubleEnded(self, doubleEnd=None):
        """
        Generate a Transition State geometry using the double-ended search method
        
        Returns (mopac, fromDbl, labels, notes) where mopac and fromDbl are 
        booleans (fromDbl is always True), and notes is a string of comments on what happened.
        """
        assert doubleEnd is not None and len(doubleEnd)==2, "You must provide the two ends of the search using 'doubleEnd' argument."
        notes = ''
        if os.path.exists(os.path.join(self.file_store_path, self.uniqueID + '.data')):
            logging.info("Not generating TS geometry because it's already done.")
            return True, None, None, "Already done!"

        reactant = doubleEnd[0]
        product = doubleEnd[1]

        rRDMol, rBM, rMult, self.geometry = self.generateBoundsMatrix(reactant)
        pRDMol, pBM, pMult, pGeom = self.generateBoundsMatrix(product)
        
        # # Smooth the inital matrix derived in rdkit
        # reactantSmoothingSuccessful = rdkit.DistanceGeometry.DoTriangleSmoothing(rBM)
        # productSmoothingSuccessful = rdkit.DistanceGeometry.DoTriangleSmoothing(pBM)
        
        print "Reactant original matrix (smoothed)"
        print matrixToString(rBM)
        print "Product original matrix (smoothed)"
        print matrixToString(pBM)
        
        self.geometry.uniqueID = self.uniqueID
        rBM, pBM, labels, atomMatch = self.editDoubMatrix(reactant, product, rBM, pBM)
        
        print "Reactant edited matrix"
        print matrixToString(rBM)
        print "Product edited matrix"
        print matrixToString(pBM)
        
        reactantSmoothingSuccessful = rdkit.DistanceGeometry.DoTriangleSmoothing(rBM)
        productSmoothingSuccessful  = rdkit.DistanceGeometry.DoTriangleSmoothing(pBM)
        
        if reactantSmoothingSuccessful:
            print "Reactant matrix is embeddable"
            print "Smoothed reactant matrix"
            print matrixToString(rBM)
        else:
            print "Reactant matrix is NOT embeddable"
        if productSmoothingSuccessful:
            print "Product matrix is embeddable"
            print "Smoothed product matrix"
            print matrixToString(pBM)
        else:
            print "Product matrix is NOT embeddable"
            
        if not (reactantSmoothingSuccessful and productSmoothingSuccessful):
            notes = 'Bounds matrix editing failed\n'
            return False, None, None, notes
        
        atoms = len(reactant.atoms)
        distGeomAttempts = 15*(atoms-3) # number of conformers embedded from the bounds matrix
         
        rdmol, minEid = self.geometry.rd_embed(rRDMol, distGeomAttempts, bm=rBM, match=atomMatch)
        if not rdmol:
            print "RDKit failed all attempts to embed"
            notes = notes + "RDKit failed all attempts to embed"
            return False, None, None, notes
        rRDMol = rdkit.Chem.MolFromMolFile(self.geometry.getCrudeMolFilePath(), removeHs=False)
        # Make product pRDMol a copy of the reactant rRDMol geometry
        for atom in reactant.atoms:
            i = atom.sortingLabel
            pRDMol.GetConformer(0).SetAtomPosition(i, rRDMol.GetConformer(0).GetAtomPosition(i))

        # don't re-embed the product, just optimize at UFF, constrained with the correct bounds matrix
        pRDMol, minEid = pGeom.optimize(pRDMol, boundsMatrix=pBM, atomMatch=atomMatch)
        pGeom.writeMolFile(pRDMol, pGeom.getRefinedMolFilePath(), minEid)
             
        if os.path.exists(self.outputFilePath):
            logging.info("File {0} already exists.".format(self.outputFilePath))
            # I'm not sure why that should be a problem, but we used to do nothin in this case
            notes = notes + 'Already have an output, check the IRC\n'
            rightTS = self.verifyIRCOutputFile()
            if rightTS:
                self.writeRxnOutputFile(labels)
                return True, self.geometry, labels, notes
            else:
                return False, None, None, notes

        if self.settings.software.lower() == 'mopac':
            # all below needs to change
            print "Optimizing reactant geometry"
            self.writeGeomInputFile(freezeAtoms=labels)
            logFilePath = self.runDouble(self.inputFilePath)
            shutil.copy(logFilePath, logFilePath+'.reactant.out')
            print "Optimizing product geometry"
            self.writeGeomInputFile(freezeAtoms=labels, otherGeom=pGeom)
            logFilePath = self.runDouble(pGeom.getFilePath(self.inputFileExtension))
            shutil.copy(logFilePath, logFilePath+'.product.out')
                
            print "Product geometry referencing reactant"
            self.writeReferenceFile(freezeAtoms=labels)#inputFilePath, molFilePathForCalc, geometry, attempt, outputFile=None)
            self.writeGeoRefInputFile(pGeom, freezeAtoms=labels, otherSide=True)#inputFilePath, molFilePathForCalc, refFilePath, geometry)
            logFilePath = self.runDouble(pGeom.getFilePath(self.inputFileExtension))
            shutil.copy(logFilePath, logFilePath+'.ref1.out')
                
            if not os.path.exists(pGeom.getFilePath('.arc')):
                notes = notes + 'product .arc file does not exits\n'
                return False, None, None, notes
            
            # Reactant that references the product geometry
            print "Reactant referencing product on slope"
            self.writeReferenceFile(freezeAtoms=labels, otherGeom=pGeom)
            self.writeGeoRefInputFile(pGeom, freezeAtoms=labels)
            logFilePath = self.runDouble(self.inputFilePath)
            shutil.copy(logFilePath, logFilePath+'.ref2.out')
            
            if not os.path.exists(self.getFilePath('.arc')):
                notes = notes + 'reactant .arc file does not exits\n'
                return False, None, None, notes
            
            # Write saddle calculation file using the outputs of the reference calculations
            print "Running Saddle from optimized geometries"
            self.writeSaddleInputFile(pGeom)
            self.runDouble(self.inputFilePath)
            return True, self.geometry, labels, notes
            # # Optimize the transition state using the TS protocol
            # self.writeInputFile(1, fromQST2=True)
            # converged, cartesian = self.run()
            # 
            # if converged:
            #     notes = notes + 'Transition state converged\n'
            #     self.writeIRCFile()
            #     rightTS = self.runIRC()
            #     if rightTS:
            #         notes = notes + 'Correct geometry found\n'
            #         return True, self.geometry, labels, notes
            #     else:
            #         notes = notes + 'Failure at IRC\n'
            #         return False, None, None, notes
            # else:
            #     notes = notes + 'Transition state not converged\n'
            #     return False, None, None, notes
        elif self.settings.software.lower() == 'gaussian':
            # all below needs to change
            print "Optimizing reactant geometry"
            self.writeGeomInputFile(freezeAtoms=labels)
            logFilePath = self.runDouble(self.inputFilePath)
            rightReactant = self.checkGeometry(logFilePath, self.geometry.molecule)
            shutil.copy(logFilePath, logFilePath+'.reactant.log')
            
            print "Optimizing product geometry"
            self.writeGeomInputFile(freezeAtoms=labels, otherGeom=pGeom)
            logFilePath = self.runDouble(pGeom.getFilePath(self.inputFileExtension))
            rightProduct = self.checkGeometry(logFilePath, pGeom.molecule)
            shutil.copy(logFilePath, logFilePath+'.product.log')
            
            if not (rightReactant and rightProduct):
                if not rightReactant:
                    print "Reactant geometry failure, see:" + self.settings.fileStore
                    notes = notes + 'Reactant geometry failure\n'
                else:
                    print "Reactant geometry success"
                
                if not rightProduct:
                    print "Product geometry failure, see:" + self.settings.fileStore
                    notes = notes + 'Product geometry failure\n'
                else:
                    print "Product geometry success"
                # Don't run if the geometries have optimized to another geometry
                return False, None, None, notes
                
            print "Running QST2 from optimized geometries"
            self.writeQST2InputFile(pGeom)
            qst2, logFilePath = self.runQST2()
            shutil.copy(logFilePath, logFilePath+'.QST2.log')
            
            if not qst2:
                print "QST3 needed, see:" + self.settings.fileStore
                notes = notes + 'QST3 needed\n'
                return False, None, None, notes
                
            print "Optimizing TS once"
            self.writeInputFile(1, fromQST2=True)
            converged, internalCoord = self.run()
            shutil.copy(self.outputFilePath, self.outputFilePath+'.TS1.log')
            
            if internalCoord and not converged:
                print "Internal coordinate error, trying in cartesian"
                self.writeInputFile(2, fromQST2=True)
                converged, internalCoord = self.run()
            
            if not converged:
                notes = notes + 'Transition state failed\n'
                return False, None, None, notes
            
            if os.path.exists(self.ircOutputFilePath):
                rightTS = self.verifyIRCOutputFile()
            else:
                self.writeIRCFile()
                rightTS = self.runIRC()
            
            if not rightTS:
                notes = notes + 'IRC failed\n'
                return False, None, None, notes
            
            self.writeRxnOutputFile(labels, doubleEnd=True)
            return True, None, None, notes
        else:
            raise NotImplementedError("self.settings.software.lower() should be gaussian or mopac")
            return False, None, None, notes
    
    def runInterplolation(self, pGeom):
        """
        Takes the reactant geometry (in `self`) and the product geometry (`pGeom`)
        and does an interpolation. This can run nudged-elastic band calculations using
        the Atomic Simulation Environment (`ASE <https://wiki.fysik.dtu.dk/ase/>`).
        """
        import ase
        from ase import io, Atoms
        from ase.neb import NEB, SingleCalculatorNEB
        import ase.calculators
        from ase.calculators.emt import EMT
        from ase.calculators.mopac import Mopac
        from ase.calculators.gaussian import Gaussian
        from ase.optimize import BFGS, FIRE
        
        # Give ase the atom positions for each side of the reaction path
        #initial = ase.io.read(self.outputFilePath, format='gaussian_out')
        # final = ase.io.read(pGeom.getFilePath(self.outputFileExtension), format='gaussian_out')
        
        # ASE doesn't keep the atoms in the same order as it's positions (weird!),
        # so get the correct atom list and recreate the images
        molfileR = self.geometry.getRefinedMolFilePath()
        molfileP = pGeom.getRefinedMolFilePath()
        atomline = re.compile('\s*([\- ][0-9.]+\s+[\-0-9.]+\s+[\-0-9.]+)\s+([A-Za-z]+)')
        
        atomCount = 0
        atomnos = []
        atomcoords = []
        with open(molfileR) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    atomnos.append(match.group(2))
                    position = numpy.array([float(i) for i in match.group(1).split()])
                    atomcoords.append(position)
                    atomCount += 1
        atomcoords = numpy.array(atomcoords)
        
        
        newImage = Atoms(atomnos)
        newImage.set_positions(atomcoords)
        initial = newImage.copy()
        
        atomCount = 0
        atomnos = []
        atomcoords = []
        with open(molfileP) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    atomnos.append(match.group(2))
                    position = numpy.array([float(i) for i in match.group(1).split()])
                    atomcoords.append(position)
                    atomCount += 1
        atomcoords = numpy.array(atomcoords)
                    
        
        newImage = Atoms(atomnos)
        newImage.set_positions(atomcoords)
        final = newImage.copy()
        
        # Now make a band of x + 2 images (x plus the initial and final geometries)
        x = 5
        images = [initial]
        images += [initial.copy() for i in range(x)]
        images += [final]
        
        # We use the linear NEB, but we can use the 'climbing image' variation by adding `climb=True`
        neb = ase.neb.NEB(images, climb=False)
        
        # Interpolate the positions of the middle images linearly, then set calculators
        neb.interpolate()
        
        # Set up the calculator
        #calc = ase.calculators.emt.EMT()
        #calc.set(multiplicity=self.geometry.molecule.getRadicalCount() + 1)

        calc = ase.calculators.mopac.Mopac(command=mopac.Mopac.executablePath, functional='PM7')
        calc.set(spin=self.geometry.molecule.getRadicalCount() )
        
        for image in images[1:x+1]:
            image.set_calculator(calc)
        
        if os.path.exists('NEB.log'): os.unlink('NEB.log')
        optimizer = BFGS(neb, trajectory='trajNEB.traj', logfile='NEB.log')
        optimizer.run(fmax=0.05)
         
        for j, image in enumerate(neb.images):
            image.write('optimized' + str(j+1), format='xyz')
        
        
    def generateTSGeometryNEB(self, doubleEnd=None):
        """
        Generate a Transition State geometry using the double-ended search method
        
        Returns (mopac, fromDbl, labels, notes) where mopac and fromDbl are 
        booleans (fromDbl is always True), and notes is a string of comments on what happened.
        """
        assert doubleEnd is not None and len(doubleEnd)==2, "You must provide the two ends of the search using 'doubleEnd' argument."
        notes = ''
        if os.path.exists(os.path.join(self.file_store_path, self.uniqueID + '.data')):
            logging.info("Not generating TS geometry because it's already done.")
            return True, None, None, "Already done!"
    
        reactant = doubleEnd[0]
        product = doubleEnd[1]
    
        rRDMol, rBM, rMult, self.geometry = self.generateBoundsMatrix(reactant)
        pRDMol, pBM, pMult, pGeom = self.generateBoundsMatrix(product)
        
        print "Reactant original matrix (smoothed)"
        print matrixToString(rBM)
        print "Product original matrix (smoothed)"
        print matrixToString(pBM)
        
        self.geometry.uniqueID = self.uniqueID
        rBM, pBM, labels, atomMatch = self.editDoubMatrix(reactant, product, rBM, pBM)
        
        print "Reactant edited matrix"
        print matrixToString(rBM)
        print "Product edited matrix"
        print matrixToString(pBM)
        
        reactantSmoothingSuccessful = rdkit.DistanceGeometry.DoTriangleSmoothing(rBM)
        productSmoothingSuccessful  = rdkit.DistanceGeometry.DoTriangleSmoothing(pBM)
        
        if reactantSmoothingSuccessful:
            print "Reactant matrix is embeddable"
            print "Smoothed reactant matrix"
            print matrixToString(rBM)
        else:
            print "Reactant matrix is NOT embeddable"
        if productSmoothingSuccessful:
            print "Product matrix is embeddable"
            print "Smoothed product matrix"
            print matrixToString(pBM)
        else:
            print "Product matrix is NOT embeddable"
            
        if not (reactantSmoothingSuccessful and productSmoothingSuccessful):
            notes = 'Bounds matrix editing failed\n'
            return False, None, None, notes
        
        atoms = len(reactant.atoms)
        distGeomAttempts = 15*(atoms-3) # number of conformers embedded from the bounds matrix
         
        rdmol, minEid = self.geometry.rd_embed(rRDMol, distGeomAttempts, bm=rBM, match=atomMatch)
        if not rdmol:
            print "RDKit failed all attempts to embed"
            notes = notes + "RDKit failed all attempts to embed"
            return False, None, None, notes
        rRDMol = rdkit.Chem.MolFromMolFile(self.geometry.getCrudeMolFilePath(), removeHs=False)
        # Make product pRDMol a copy of the reactant rRDMol geometry
        for atom in reactant.atoms:
            i = atom.sortingLabel
            pRDMol.GetConformer(0).SetAtomPosition(i, rRDMol.GetConformer(0).GetAtomPosition(i))
    
        # don't re-embed the product, just optimize at UFF, constrained with the correct bounds matrix
        pRDMol, minEid = pGeom.optimize(pRDMol, boundsMatrix=pBM, atomMatch=atomMatch)
        pGeom.writeMolFile(pRDMol, pGeom.getRefinedMolFilePath(), minEid)
             
        # if os.path.exists(self.outputFilePath):
        #     logging.info("File {0} already exists.".format(self.outputFilePath))
        #     # I'm not sure why that should be a problem, but we used to do nothin in this case
        #     notes = notes + 'Already have an output, check the IRC\n'
        #     rightTS = self.verifyIRCOutputFile()
        #     if rightTS:
        #         self.writeRxnOutputFile(labels)
        #         return True, self.geometry, labels, notes
        #     else:
        #         return False, None, None, notes
    
#         if self.settings.software.lower() == 'gaussian' or True:
#             # all below needs to change
#             if os.path.exists(self.getFilePath('.log.reactant.log')):
#                 print "Already have reactant"
#                 rightReactant = self.checkGeometry(self.getFilePath('.log.reactant.log'), self.geometry.molecule)
#             else:
#                 print "Optimizing reactant geometry"
#                 self.writeGeomInputFile(freezeAtoms=labels)
#                 logFilePath = self.runDouble(self.inputFilePath)
#                 rightReactant = self.checkGeometry(logFilePath, self.geometry.molecule)
#                 shutil.copy(logFilePath, logFilePath+'.reactant.log')
#             
#             if os.path.exists(pGeom.getFilePath('.log.product.log')):
#                 print "Already have product"
#                 rightProduct = self.checkGeometry(pGeom.getFilePath('.log.product.log'), pGeom.molecule)
#             else:
#                 print "Optimizing product geometry"
#                 self.writeGeomInputFile(freezeAtoms=labels, otherGeom=pGeom)
#                 logFilePath = self.runDouble(pGeom.getFilePath(self.inputFileExtension))
#                 rightProduct = self.checkGeometry(logFilePath, pGeom.molecule)
#                 shutil.copy(logFilePath, logFilePath+'.product.log')
#             
#             if not (rightReactant and rightProduct):
#                 if not rightReactant:
#                     print "Reactant geometry failure, see:" + self.settings.fileStore
#                     notes = notes + 'Reactant geometry failure\n'
#                 else:
#                     print "Reactant geometry success"
#                 
#                 if not rightProduct:
#                     print "Product geometry failure, see:" + self.settings.fileStore
#                     notes = notes + 'Product geometry failure\n'
#                 else:
#                     print "Product geometry success"
#                 # Don't run if the geometries have optimized to another geometry
#                 return False, None, None, notes
                
        print "Running NEB from optimized geometries"
        # Atomic Simulation Environment can take the two geometries and
        # do the calculation on its own
        self.runInterplolation(pGeom)
            
        print "Optimizing TS once"
        self.writeInputFile(1, fromNEB=True)
        converged, internalCoord = self.run()
        shutil.copy(self.outputFilePath, self.outputFilePath+'.TS1.log')
        
        if internalCoord and not converged:
            print "Internal coordinate error, trying in cartesian"
            self.writeInputFile(2, fromInt=True)
            converged, internalCoord = self.run()
        
        if not converged:
            notes = notes + 'Transition state failed\n'
            return False, None, None, notes
        
        if os.path.exists(self.ircOutputFilePath):
            rightTS = self.verifyIRCOutputFile()
        else:
            self.writeIRCFile()
            rightTS = self.runIRC()
        
        if not rightTS:
            notes = notes + 'IRC failed\n'
            return False, None, None, notes
        
        self.writeRxnOutputFile(labels, doubleEnd=True)
        return True, None, None, notes


    def generateTSGeometryDirectGuess(self):
        """
        Generate a transition state geometry, using the direct guess (group additive) method.
        
        Returns (success, notes) where success is a True if it worked, else False,
        and notes is a string describing what happened.
        """
        notes = ''
        if os.path.exists(os.path.join(self.file_store_path, self.uniqueID + '.data')):
            logging.info("Not generating TS geometry because it's already done.")
            return True, "Already done!"
        
        if len(self.reaction.reactants)==2:
            reactant = self.reaction.reactants[0].merge(self.reaction.reactants[1])
        else:
            reactant = self.reaction.reactants[0]
        if len(self.reaction.products)==2:
            product = self.reaction.products[0].merge(self.reaction.products[1])
        else:
            product = self.reaction.products[0]
        
        reactant = self.fixSortLabel(reactant)
        product = self.fixSortLabel(product)
        
        tsRDMol, tsBM, tsMult, self.geometry = self.generateBoundsMatrix(reactant)
        
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
                print "Optimizing TS once"
                self.writeInputFile(1)
                converged, internalCoord = self.run()
            else:
                converged, internalCoord = self.verifyOutputFile()
            
            if internalCoord and not converged:
                notes = 'Internal coordinate error, trying cartesian\n'
                print "Optimizing TS in cartesian"
                shutil.copy(self.outputFilePath, self.outputFilePath+'.TS1.log')
                self.writeInputFile(2)
                converged = self.run()
                
            if converged:
                notes = 'TS converged, now for IRC\n'
                print "IRC calculation"
                if not os.path.exists(self.ircOutputFilePath):
                    self.writeIRCFile()
                    rightTS = self.runIRC()
                else:
                    rightTS = self.verifyIRCOutputFile()
                if rightTS:
                    print "Found a transition state"
                    notes = 'Success\n'
                    self.writeRxnOutputFile(labels)
                    return True, notes
                else:
                    print "Graph matching failed"
                    notes = 'IRC failed\n'
                    return False, notes
            else:
                print "TS failed"
                notes = 'TS not converged\n'
                return False, notes
        else:
            notes = 'Bounds matrix editing failed\n'
            return False, notes
    
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
        tsFound = self.generateTSGeometryDirectGuess()
        
        if not tsFound:
            # fall back on group additivity
            return None
            
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
    
    def determinePointGroup(self):
        """
        Determine point group using the SYMMETRY Program
        
        Stores the resulting :class:`PointGroup` in self.pointGroup
        """
        assert self.qmData, "Need QM Data first in order to calculate point group."
        pgc = symmetry.PointGroupCalculator(self.settings, self.uniqueID, self.qmData)
        self.pointGroup = pgc.calculate()
