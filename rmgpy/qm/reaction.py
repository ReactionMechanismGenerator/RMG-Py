#!/usr/bin/env python
# encoding: utf-8

import os
import re
import logging
import external.cclib as cclib
import time
from subprocess import Popen
from copy import deepcopy
import numpy
import shutil
import math
import sqlite3 as lite

from rmgpy.data.kinetics.transitionstates import TransitionStates
from rmgpy.molecule import Molecule, Atom, getElement
from rmgpy.species import Species, TransitionState
from rmgpy.kinetics import Wigner
from molecule import QMMolecule, Geometry
from rmgpy.cantherm.main import CanTherm
from rmgpy.cantherm.kinetics import KineticsJob
import symmetry

try:
    import rdkit
    from rdkit import DistanceGeometry
    from rdkit.Chem.Pharm3D import EmbedLib
except ImportError:
    logging.info("To use transition state searches, you must correctly install rdkit")

def matrixToString(matrix):
    """Returns a string representation of a matrix, for printing to the console"""
    text = '\n'.join([ ' '.join([str(round(item, 1)) for item in line]) for line in matrix ])
    return text.replace('1000.0', '1e3')


class QMReaction:
    """
    A base class for QM Reaction calculations.

    Specific programs and methods should inherit from this and define some
    extra attributes and methods:

    The attributes are:

    =================== =========================== ====================================
    Attribute           Type                        Description
    =================== =========================== ====================================
    `reaction`          :class:`Reaction`           RMG Reaction object
    `settings`          :class:`QMSettings`         Settings for QM calculations
    `database`                                      The RMG Database
    `uniqueID`          ``str``                     A reaction ID with canonical smiles of the reactants and products
    `geometry`          :class:`Geometry`           Geometry object for handling 3-dimensional structures
    `transitionState`   :class:`TransitionState`    Storage of the transition state object
    =================== =========================== ====================================

    """

    def __init__(self, reaction, settings, tsDatabase):
        self.reaction = reaction
        self.settings = settings
        self.tsDatabase = tsDatabase
        # self.tsDatabase = TransitionStates()
        # self.tsDatabase.load(
        #         path=os.path.abspath(os.path.join(os.getenv('RMGpy'), '..', 'RMG-database', 'input', 'kinetics', 'families', self.reaction.family.label)),
        #         local_context={},
        #         global_context={},
        #         )

        if isinstance(self.reaction.reactants[0], Molecule):
            reactants = sorted([s.toSMILES() for s in self.reaction.reactants])
            products = sorted([s.toSMILES() for s in self.reaction.products])
        elif isinstance(self.reaction.reactants[0], Species):
            reactants = sorted([s.molecule[0].toSMILES() for s in self.reaction.reactants])
            products = sorted([s.molecule[0].toSMILES() for s in self.reaction.products])
        stringID = "+".join(reactants) + "_" + "+".join(products)
        revID = "+".join(products) + "_" + "+".join(reactants)

        self.uniqueID = stringID
        self.revID = revID

        # self.geometry = None
        self.reactantGeom = None # Geometry(settings, molecule.toAugmentedInChIKey(), molecule)
        self.productGeom = None
        self.tsGeom = None

        self.fileStore = self.settings.fileStore
        self.scratchDirectory = self.settings.scratchDirectory

        if not os.path.exists(self.fileStore):
            logging.info("Creating permanent directory %s for qm files."%os.path.abspath(self.fileStore))
            os.makedirs(self.fileStore)

        if not os.path.exists(self.scratchDirectory):
            logging.info("Creating scratch directory %s for qm files."%os.path.abspath(self.scratchDirectory))
            os.makedirs(self.scratchDirectory)

    def getFilePath(self, extension):
        """
        Should return the path to the file with the given extension.

        The provided extension should include the leading dot.

        Need to define some reaction line notation.
        Possibly '<Reaction_Family>/<reactant1SMILES>+<reactant2SMILES>--<product1SMILES>+<product2SMILES>' ???
        """
        return os.path.join(self.fileStore, self.uniqueID  + extension)

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

    def setOutputDirectory(self):
        """
        Set up the fileStore and scratchDirectory if not already done.
        """
        subPath = os.path.join('Reactions', self.reaction.family.label, self.uniqueID, self.settings.method)

        setFileStore = True
        setScratch = True
        if self.fileStore:
            if not self.fileStore.endswith(subPath):
                self.fileStore = os.path.join(self.settings.fileStore, subPath)
                logging.info("Setting the qm kinetics fileStore to {0}".format(self.settings.fileStore))
            setFileStore = False

        if self.scratchDirectory:
            if not self.scratchDirectory.endswith(subPath):
                self.scratchDirectory = os.path.join(self.settings.scratchDirectory, subPath)
                logging.info("Setting the qm kinetics scratchDirectory fileStore to {0}".format(self.settings.scratchDirectory))
            setScratch = False

        if setFileStore:
            self.fileStore = os.path.join('QMfiles', subPath)
            logging.info("Set the qm kinetics fileStore to {0}".format(self.fileStore))
        if setScratch:
            self.scratchDirectory = os.path.join('QMscratch', subPath)
            logging.info("Set the qm kinetics scratchDirectory to {0}".format(self.scratchDirectory))

    def initialize(self):
        """
        Do any startup tasks.
        """
        self.checkReady()

    def checkReady(self):
        """
        Check that it's ready to run calculations.
        """
        self.settings.checkAllSet()
        self.checkPaths()

    def checkPaths(self):
        """
        Check the paths in the settings are OK. Make folders as necessary.
        """
        if not os.path.exists(self.settings.RMG_bin_path):
            raise Exception("RMG-Py 'bin' directory {0} does not exist.".format(self.settings.RMG_bin_path))
        if not os.path.isdir(self.settings.RMG_bin_path):
            raise Exception("RMG-Py 'bin' directory {0} is not a directory.".format(self.settings.RMG_bin_path))

        self.setOutputDirectory()
        self.fileStore = os.path.expandvars(self.fileStore) # to allow things like $HOME or $RMGpy
        self.scratchDirectory = os.path.expandvars(self.scratchDirectory)
        for path in [self.fileStore, self.scratchDirectory]:
            if not os.path.exists(path):
                logging.info("Creating directory %s for QM files."%os.path.abspath(path))
                os.makedirs(path)

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

        geom = Geometry(settings, molecule.toAugmentedInChIKey(), molecule)

        return geom

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
        geometry = self.getGeometry(molecule, self.settings)
        rdKitMol = self.getRDKitMol(geometry)
        boundsMatrix = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(rdKitMol)

        return rdKitMol, boundsMatrix, geometry

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

    def bmLenComp(self, bm, i, j):

        listDist = []
        for k in range(len(bm)):
            if k==i or k==j or i==j: continue
            Uik = bm[i,k] if k>i else bm[k,i]
            Ukj = bm[j,k] if k>j else bm[k,j]
            if i>j:
                dist = Uik + Ukj - bm[i,j]
            else:
                dist = Uik + Ukj - bm[j,i]
            listDist.append(dist)

        return listDist
            # maxLij = Uik + Ukj - 0.1
            # if bm[i,j] >  maxLij:
            #     print "CHANGING Lower limit {0} to {1}".format(bm[i,j], maxLij)
            #     bm[i,j] = maxLij

    def getLabels(self, reactant):
        """
        Creates the list of sorting labels for the reacting atoms. These labels are
        also put intoa tuple for the atomMatch weighting for RDKit. The weighting tells
        RDKit to place greater importance in maintaining these distance limits when
        generating conformers.
        """
        if self.reaction.family.label.lower() in ['h_abstraction', 'r_addition_multiplebond', 'intra_h_migration']:
            lbl1 = reactant.getLabeledAtom('*1').sortingLabel
            lbl2 = reactant.getLabeledAtom('*2').sortingLabel
            lbl3 = reactant.getLabeledAtom('*3').sortingLabel
        elif self.reaction.family.label.lower() in ['disproportionation']:
            lbl1 = reactant.getLabeledAtom('*2').sortingLabel
            lbl2 = reactant.getLabeledAtom('*4').sortingLabel
            lbl3 = reactant.getLabeledAtom('*1').sortingLabel

        labels = [lbl1, lbl2, lbl3]
        atomMatch = ((lbl1,),(lbl2,),(lbl3,))

        return labels, atomMatch

    def editDoubMatrix(self, reactant, product, bm1, bm2, labels):
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

        lbl1, lbl2, lbl3 = labels

        if reactant.atoms[lbl1].symbol == 'H' or reactant.atoms[lbl3].symbol == 'H':
            bm1 = fixMatrix(bm1, lbl1, lbl2, lbl3, 2.3, 0.1)
            bm2 = fixMatrix(bm2, lbl3, lbl2, lbl1, 2.3, 0.1)
        else:
            bm1 = fixMatrix(bm1, lbl1, lbl2, lbl3, 2.7, 0.1)
            bm2 = fixMatrix(bm2, lbl3, lbl2, lbl1, 2.7, 0.1)

        rSect = []
        for atom in reactant.split()[0].atoms: rSect.append(atom.sortingLabel)

        pSect = []
        for atom in product.split()[0].atoms: pSect.append(atom.sortingLabel)

        bm1 = self.bmPreEdit(bm1, rSect)
        bm2 = self.bmPreEdit(bm2, pSect)

        return bm1, bm2

    def editMatrix(self, reactant, bm, labels):

        """
        For bimolecular reactions, reduce the minimum distance between atoms
        of the two reactants.
        """
        lbl1, lbl2, lbl3 = labels

        distanceData = self.tsDatabase.estimateDistances(self.reaction)

        sect = []
        for atom in reactant.split()[0].atoms: sect.append(atom.sortingLabel)

        uncertainties = {'d12':0.02, 'd13':0.02, 'd23':0.02 } # distanceData.uncertainties or {'d12':0.1, 'd13':0.1, 'd23':0.1 } # default if uncertainty is None
        bm = self.setLimits(bm, lbl1, lbl2, distanceData.distances['d12'], uncertainties['d12'])
        bm = self.setLimits(bm, lbl2, lbl3, distanceData.distances['d23'], uncertainties['d23'])
        bm = self.setLimits(bm, lbl1, lbl3, distanceData.distances['d13'], uncertainties['d13'])

        with open(os.path.join(self.fileStore, 'estDists.txt'), 'w') as distFile:
            distFile.write('d12: {0:.6f}, d13: {1:.6f}, d23: {2:.6f}'.format(distanceData.distances['d12'], distanceData.distances['d13'], distanceData.distances['d23']))

        bm = self.bmPreEdit(bm, sect)

        return bm

    def runNEB(self):
        """
        Takes the reactant and product geometries and does an interpolation. This
        can run nudged-elastic band calculations using the Atomic Simulation
        Environment (`ASE <https://wiki.fysik.dtu.dk/ase/>`).
        """
        import ase
        from ase.neb import NEB
        from ase.optimize import BFGS, FIRE

        initial, final = self.setImages()

        # Now make a band of x + 2 images (x plus the initial and final geometries)
        x = 11
        images = [initial]
        images += [initial.copy() for i in range(x)]
        images += [final]

        # We use the linear NEB, but we can use the 'climbing image' variation by adding `climb=True`
        # Options recommended in Gonzalez-Garcia et al., doi: 10.1021/ct060032y (They do recommend climbing image, but I'll test without if first)
        neb = ase.neb.NEB(images, k=0.01, climb=False)

        # Interpolate the positions of the middle images linearly, then set calculators
        neb.interpolate()

        self.setCalculator(images)

        nebLog = os.path.join(self.settings.fileStore, 'NEB.log')
        if os.path.exists(nebLog): os.unlink(nebLog)
        optimizer = BFGS(neb, logfile=nebLog)
        optimized = True
        try:
            optimizer.run(steps=200)
        except Exception, e:
            print str(e)
            optimized = False
            pass

        if optimized:
            energies = numpy.empty(neb.nimages - 2)
            for i in range(1, neb.nimages - 1):
                energies[i - 1] = neb.images[i].get_potential_energy()
            imax = 1 + numpy.argsort(energies)[-1]
            image = neb.images[imax]
            image.write(self.getFilePath('peak.xyz'), format='xyz')
        else:
            print "Not optimized"

    def setupMolecules(self, doubleEnded=False):
        """
        Setup the reactant and product molecules for the transition state calculations.
        If there are 2 reactants and/or products, they are merged. This also handles
        species as well as molecules, but returns the reactant and product as merged molecules.
        """
        if doubleEnded:
            kineticsFamily = self.reaction.family
            if isinstance(self.reaction.reactants[0], Species):
                reactants = []
                for reactant in self.reaction.reactants:
                    reactants.append(reactant.molecule[0])
            else:
                reactants = self.reaction.reactants
            prodStruct, tsStructures = kineticsFamily.applyRecipe(reactants, getTS=True)

            reactant = tsStructures[0]
            product = tsStructures[1]
        else:
            if len(self.reaction.reactants)==2:
                if isinstance(self.reaction.reactants[0], Molecule):
                    reactant = self.reaction.reactants[0].merge(self.reaction.reactants[1])
                elif isinstance(self.reaction.reactants[0], Species):
                    reactant = self.reaction.reactants[0].molecule[0].merge(self.reaction.reactants[1].molecule[0])
            else:
                if isinstance(self.reaction.reactants[0], Molecule):
                    reactant = self.reaction.reactants[0]
                elif isinstance(self.reaction.reactants[0], Species):
                    reactant = self.reaction.reactants[0].molecule[0]

            if len(self.reaction.products)==2:
                if isinstance(self.reaction.reactants[0], Molecule):
                    product = self.reaction.products[0].merge(self.reaction.products[1])
                elif isinstance(self.reaction.reactants[0], Species):
                    product = self.reaction.products[0].molecule[0].merge(self.reaction.products[1].molecule[0])
            else:
                if isinstance(self.reaction.reactants[0], Molecule):
                    product = self.reaction.products[0]
                elif isinstance(self.reaction.reactants[0], Species):
                    product = self.reaction.products[0].molecule[0]

        reactant.multiplicity = reactant.getRadicalCount() + 1
        product.multiplicity = product.getRadicalCount() + 1

        reactant = self.fixSortLabel(reactant)
        product = self.fixSortLabel(product)

        return reactant, product

    def optimizeTS(self, labels, fromDoubleEnded=False):
        """
        Conduct the optimization step of the transition state search.
        """
        # if os.path.exists(self.outputFilePath):
        #     os.remove(checkpointFile)
            # complete = self.checkComplete(self.outputFilePath)
            # 
            # if complete:
            #     converged, internalCoord = self.verifyOutputFile()
            # else:
            #     # Delete the output and checkpoint files so we redo the calc
            #     os.remove(self.outputFilePath)
            #     checkpointFile = os.path.join(self.settings.fileStore, self.uniqueID + ".chk")
            #     assert os.path.exists(checkpointFile)
            #     os.remove(checkpointFile) # Checkpoint file path
        if os.path.exists(self.outputFilePath):
            converged, internalCoord = self.verifyOutputFile()
        else:
            optEst = self.optEstimate(labels)
            optRC = self.optRxnCenter(labels)
            print "Optimizing TS once"
            self.createInputFile(1, fromDoubleEnded=fromDoubleEnded, optEst=optRC)
            converged, internalCoord = self.run()
            shutil.copy(self.outputFilePath, self.outputFilePath+'.TS1.log')
        
            if os.path.exists(self.ircOutputFilePath):# Remove it
                os.remove(self.ircOutputFilePath)

        if internalCoord and not converged:
            notes = 'Internal coordinate error, trying cartesian\n'
            print "Optimizing TS in cartesian"
            self.createInputFile(2)
            converged = self.run()
            shutil.copy(self.outputFilePath, self.outputFilePath+'.TS2.log')
            
        # if not os.path.exists(self.outputFilePath):
        #     optEst = self.getFilePath('Est{0}'.format(self.outputFileExtension))
        #     optRC = self.getFilePath('RxnC{0}'.format(self.outputFileExtension))
        #     if os.path.exists(optEst):
        #         optEst = self.optEstimate(labels)
        #     if os.path.exists(optRC):
        #         optRC = self.optRxnCenter(labels)
        #     print "Optimizing TS once"
        #     self.createInputFile(1, fromDoubleEnded=fromDoubleEnded, optEst=optRC)
        #     converged, internalCoord = self.run()
        #     shutil.copy(self.outputFilePath, self.outputFilePath+'.TS1.log')
        # 
        #     if internalCoord and not converged:
        #         notes = 'Internal coordinate error, trying cartesian\n'
        #         print "Optimizing TS in cartesian"
        #         self.createInputFile(2)
        #         converged = self.run()
        #         shutil.copy(self.outputFilePath, self.outputFilePath+'.TS2.log')
        # else:
        #     converged = self.verifyOutputFile()
        
        if not converged:
            # Check for convergence failures
            complete, convergenceFailure = self.checkComplete(self.outputFilePath)
            if convergenceFailure:
                # Rerun the calculation with `scf=qc`
                self.createInputFile(1, fromDoubleEnded=fromDoubleEnded, optEst=optRC, scf=True)
                converged, internalCoord = self.run()
                if internalCoord:
                    self.createInputFile(2, scf=True)
                    converged = self.run()
                shutil.copy(self.outputFilePath, self.outputFilePath+'.QC.log')

        return converged

    def validateTS(self):
        """
        Conduct the path analysis calculations and validate the transition state.
        """
        if not os.path.exists(self.ircOutputFilePath):
            self.createIRCFile()
            rightTS = self.runIRC()
        else:
            complete, convergenceFailure = self.checkComplete(self.ircOutputFilePath)
            if convergenceFailure:
                rightTS = False
            else:
                rightTS = self.verifyIRCOutputFile()
        
        if not rightTS:
            # Check for convergence failures
            complete, convergenceFailure = self.checkComplete(self.ircOutputFilePath)
            if convergenceFailure:
                # Rerun the calculation with `scf=qc`
                self.createIRCFile(scf=True)
                rightTS = self.runIRC()
        
        return rightTS
    
    def checkSQL(self, path_to_db):
        """
        Check the SQL database for a transition state. This should be run before
        the automated TS geometry search.
        """
        con = lite.connect(os.path.join(path_to_db, 'ts_data.db'))
        with con:
            cur = con.cursor()
            cur.execute("SELECT * FROM TSData")
            
            rows = cur.fetchall()
            
            for row in rows:
                print row

    def tsSearch(self, notes, labels, fromDoubleEnded=False):
        """
        Once the transition state estimate is made, this runs the optimization and the
        path analysis calculation. The ts estimate can be from the group additive or
        double-ended search methods.
        """
        # Check SQL database for transition state
        # self.checkSQL('/scratch/westgroup')

        successfulTS = self.optimizeTS(labels, fromDoubleEnded=fromDoubleEnded)
        if not successfulTS:
            notes = 'TS not converged\n'
            return successfulTS, notes

        validTS = self.validateTS()
        if not validTS:
            notes = 'IRC failed\n'
        else:
            self.writeRxnOutputFile(labels)
            notes = 'Success\n'

        return validTS, notes
            
    def generateTSGeometryDirectGuess(self):
        """
        Generate a transition state geometry, using the direct guess (group additive) method.

        Returns (success, notes) where success is a True if it worked, else False,
        and notes is a string describing what happened.
        """
        def getDistance(coordinates1, coordinates2):
            """
            Return the square of the distance (in Angstrom) between the two atoms.
            """
            diff = (coordinates1.coords - coordinates2.coords)
            return math.sqrt(sum(diff * diff))
            
        self.settings.fileStore = self.fileStore
        self.settings.scratchDirectory = self.scratchDirectory
        split_fileStore = self.fileStore.split(self.uniqueID)
        rev_fileStore = self.revID.join(split_fileStore)
        split_scratch = self.scratchDirectory.split(self.uniqueID)
        rev_scratch = self.revID.join(split_scratch)
        notes = ''

        reactant, product = self.setupMolecules()

        tsRDMol, tsBM, self.reactantGeom = self.generateBoundsMatrix(reactant)

        self.reactantGeom.uniqueID = self.uniqueID

        labels, atomMatch = self.getLabels(reactant)
        
        if os.path.exists(os.path.join(self.fileStore, self.uniqueID + '.data')):
            estFilePath = self.getFilePath('Est{0}'.format(self.outputFileExtension))
            rcFilePath = self.getFilePath('RxnC{0}'.format(self.outputFileExtension))
            if os.path.exists(estFilePath):
                parser = cclib.parser.Gaussian(estFilePath)
                parser.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
                cclib_data = parser.parse()
                atomNums = cclib_data.atomnos
                atomCoords = cclib_data.atomcoords[-1]
                
                atom1 = Atom(element=getElement(int(atomNums[labels[0]])), coords=atomCoords[labels[0]])
                atom2 = Atom(element=getElement(int(atomNums[labels[1]])), coords=atomCoords[labels[1]])
                atom3 = Atom(element=getElement(int(atomNums[labels[2]])), coords=atomCoords[labels[2]])

                at12 = getDistance(atom1, atom2)
                at23 = getDistance(atom2, atom3)
                at13 = getDistance(atom1, atom3)
                
                with open(os.path.join(self.fileStore, 'fzEstDists.txt'), 'w') as distFile:
                    distFile.write('d12: {0:.6f}, d13: {1:.6f}, d23: {2:.6f}'.format(at12, at13, at23))
            if os.path.exists(rcFilePath):
                parser = cclib.parser.Gaussian(rcFilePath)
                parser.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
                cclib_data = parser.parse()
                atomNums = cclib_data.atomnos
                atomCoords = cclib_data.atomcoords[-1]
                
                atom1 = Atom(element=getElement(int(atomNums[labels[0]])), coords=atomCoords[labels[0]])
                atom2 = Atom(element=getElement(int(atomNums[labels[1]])), coords=atomCoords[labels[1]])
                atom3 = Atom(element=getElement(int(atomNums[labels[2]])), coords=atomCoords[labels[2]])

                at12 = getDistance(atom1, atom2)
                at23 = getDistance(atom2, atom3)
                at13 = getDistance(atom1, atom3)
                
                with open(os.path.join(self.fileStore, 'rcDists.txt'), 'w') as distFile:
                    distFile.write('d12: {0:.6f}, d13: {1:.6f}, d23: {2:.6f}'.format(at12, at13, at23))
            parser = cclib.parser.Gaussian(self.outputFilePath)
            parser.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
            cclib_data = parser.parse()
            atomNums = cclib_data.atomnos
            atomCoords = cclib_data.atomcoords[-1]
            
            atom1 = Atom(element=getElement(int(atomNums[labels[0]])), coords=atomCoords[labels[0]])
            atom2 = Atom(element=getElement(int(atomNums[labels[1]])), coords=atomCoords[labels[1]])
            atom3 = Atom(element=getElement(int(atomNums[labels[2]])), coords=atomCoords[labels[2]])

            at12 = getDistance(atom1, atom2)
            at23 = getDistance(atom2, atom3)
            at13 = getDistance(atom1, atom3)
            
            with open(os.path.join(self.fileStore, 'optDists.txt'), 'w') as distFile:
                distFile.write('d12: {0:.6f}, d13: {1:.6f}, d23: {2:.6f}'.format(at12, at13, at23))
            return True, "Already done!"
        
        tsBM = self.editMatrix(reactant, tsBM, labels)
        atoms = len(reactant.atoms)
        if atoms>3:
            distGeomAttempts = 15*(atoms-3) # number of conformers embedded from the bounds matrix
        else:
            distGeomAttempts = 15

        setBM = rdkit.DistanceGeometry.DoTriangleSmoothing(tsBM)

        if not setBM:
            notes = 'Bounds matrix editing failed\n'
            return False, notes

        for i in range(len(tsBM)):
            for j in range(i,len(tsBM)):
                if tsBM[j,i] > tsBM[i,j]:
                        print "BOUNDS MATRIX FLAWED {0}>{1}".format(tsBM[j,i], tsBM[i,j])

        self.reactantGeom.rd_embed(tsRDMol, distGeomAttempts, bm=tsBM, match=atomMatch)
        atomSymbols, atomCoords = self.reactantGeom.parseMOL(self.reactantGeom.getRefinedMolFilePath())

        d12sq = (atomCoords[labels[0]]-atomCoords[labels[1]])**2
        d23sq = (atomCoords[labels[1]]-atomCoords[labels[2]])**2
        d13sq = (atomCoords[labels[0]]-atomCoords[labels[2]])**2
        d12 = math.sqrt(d12sq.sum())
        d23 = math.sqrt(d23sq.sum())
        d13 = math.sqrt(d13sq.sum())

        with open(os.path.join(self.fileStore, 'rdkitDists.txt'), 'w') as distFile:
            distFile.write('d12: {0:.6f}, d13: {1:.6f}, d23: {2:.6f}'.format(d12, d13, d23))

        check, notes =  self.tsSearch(notes, labels)
        
        if check:
            estFilePath = self.getFilePath('Est{0}'.format(self.outputFileExtension))
            rcFilePath = self.getFilePath('RxnC{0}'.format(self.outputFileExtension))
            if os.path.exists(estFilePath):
                parser = cclib.parser.Gaussian(estFilePath)
                parser.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
                cclib_data = parser.parse()
                atomNums = cclib_data.atomnos
                atomCoords = cclib_data.atomcoords[-1]
                
                atom1 = Atom(element=getElement(int(atomNums[labels[0]])), coords=atomCoords[labels[0]])
                atom2 = Atom(element=getElement(int(atomNums[labels[1]])), coords=atomCoords[labels[1]])
                atom3 = Atom(element=getElement(int(atomNums[labels[2]])), coords=atomCoords[labels[2]])

                at12 = getDistance(atom1, atom2)
                at23 = getDistance(atom2, atom3)
                at13 = getDistance(atom1, atom3)
                
                with open(os.path.join(self.fileStore, 'fzEstDists.txt'), 'w') as distFile:
                    distFile.write('d12: {0:.6f}, d13: {1:.6f}, d23: {2:.6f}'.format(at12, at13, at23))
            if os.path.exists(rcFilePath):
                parser = cclib.parser.Gaussian(rcFilePath)
                parser.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
                cclib_data = parser.parse()
                atomNums = cclib_data.atomnos
                atomCoords = cclib_data.atomcoords[-1]
                
                atom1 = Atom(element=getElement(int(atomNums[labels[0]])), coords=atomCoords[labels[0]])
                atom2 = Atom(element=getElement(int(atomNums[labels[1]])), coords=atomCoords[labels[1]])
                atom3 = Atom(element=getElement(int(atomNums[labels[2]])), coords=atomCoords[labels[2]])

                at12 = getDistance(atom1, atom2)
                at23 = getDistance(atom2, atom3)
                at13 = getDistance(atom1, atom3)
                
                with open(os.path.join(self.fileStore, 'rcDists.txt'), 'w') as distFile:
                    distFile.write('d12: {0:.6f}, d13: {1:.6f}, d23: {2:.6f}'.format(at12, at13, at23))
            parser = cclib.parser.Gaussian(self.outputFilePath)
            parser.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
            cclib_data = parser.parse()
            atomNums = cclib_data.atomnos
            atomCoords = cclib_data.atomcoords[-1]
            
            atom1 = Atom(element=getElement(int(atomNums[labels[0]])), coords=atomCoords[labels[0]])
            atom2 = Atom(element=getElement(int(atomNums[labels[1]])), coords=atomCoords[labels[1]])
            atom3 = Atom(element=getElement(int(atomNums[labels[2]])), coords=atomCoords[labels[2]])

            at12 = getDistance(atom1, atom2)
            at23 = getDistance(atom2, atom3)
            at13 = getDistance(atom1, atom3)
            
            with open(os.path.join(self.fileStore, 'optDists.txt'), 'w') as distFile:
                distFile.write('d12: {0:.6f}, d13: {1:.6f}, d23: {2:.6f}'.format(at12, at13, at23))

        return check, notes

    def generateTSGeometryDoubleEnded(self, neb=False):
        """
        Generate a Transition State geometry using the double-ended search method

        Returns (mopac, fromDbl, labels, notes) where mopac and fromDbl are
        booleans (fromDbl is always True), and notes is a string of comments on what happened.
        """

        notes = ''
        if os.path.exists(os.path.join(self.settings.fileStore, self.uniqueID + '.data')):
            logging.info("Not generating TS geometry because it's already done.")
            return True, "Output used from a previous run."

        reactant, product = self.setupMolecules(doubleEnded=True)

        rRDMol, rBM, self.reactantGeom = self.generateBoundsMatrix(reactant)
        pRDMol, pBM, self.productGeom = self.generateBoundsMatrix(product)

        self.reactantGeom.uniqueID = 'reactant'
        self.productGeom.uniqueID = 'product'

        labels, atomMatch = self.getLabels(reactant)
        rBM, pBM = self.editDoubMatrix(reactant, product, rBM, pBM, labels)

        reactantSmoothingSuccessful = rdkit.DistanceGeometry.DoTriangleSmoothing(rBM)
        productSmoothingSuccessful  = rdkit.DistanceGeometry.DoTriangleSmoothing(pBM)

        if not (reactantSmoothingSuccessful and productSmoothingSuccessful):
            notes = 'Bounds matrix editing failed\n'
            return False, None, None, notes

        atoms = len(reactant.atoms)
        distGeomAttempts = 15*(atoms-3) # number of conformers embedded from the bounds matrix

        rdmol, minEid = self.reactantGeom.rd_embed(rRDMol, distGeomAttempts, bm=rBM, match=atomMatch)
        if not rdmol:
            notes = notes + "RDKit failed all attempts to embed"
            return False, None, None, notes
        rRDMol = rdkit.Chem.MolFromMolFile(self.reactantGeom.getCrudeMolFilePath(), removeHs=False)
        # Make product pRDMol a copy of the reactant rRDMol geometry
        for atom in reactant.atoms:
            i = atom.sortingLabel
            pRDMol.GetConformer(0).SetAtomPosition(i, rRDMol.GetConformer(0).GetAtomPosition(i))

        # don't re-embed the product, just optimize at UFF, constrained with the correct bounds matrix
        pRDMol, minEid = self.productGeom.optimize(pRDMol, boundsMatrix=pBM, atomMatch=atomMatch)
        self.productGeom.writeMolFile(pRDMol, self.productGeom.getRefinedMolFilePath(), minEid)

        if os.path.exists(self.outputFilePath):
            logging.info("File {0} already exists.".format(self.outputFilePath))
            # I'm not sure why that should be a problem, but we used to do nothin in this case
            notes = notes + 'Already have an output, checking the IRC\n'
            rightTS = self.verifyIRCOutputFile()
            if rightTS:
                self.writeRxnOutputFile(labels)
                return True, notes#True, self.geometry, labels, notes
            else:
                return False, notes#False, None, None, notes

        check, notes = self.conductDoubleEnded(notes, NEB=neb, labels=labels)
        if check:
            # Optimize the TS
            check, notes =  self.tsSearch(notes, labels, fromDoubleEnded=True)

        return check, notes

    def writeXYZ(self, atomSymbols, atomCoords):
        """
        Takes a list of atom symbols stings, and an array of coordiantes and
        writes a file to store the geometry.
        """
        xyzOutput = []
        for atomsymbol, atomcoord in zip(atomsymbols, atomcoords):
            inputline = "{0:4s}  {1:4f}  {2:4f}  {3:4f}".format(atomsymbol, atomcoord[0], atomcoord[1], atomcoord[2])
            inputline = inputline.replace(' -', '-')
            xyzOutput.append(inputline)

        input_string = '\n'.join(xyzOutput)

        with open(self.getFilePath('.xyz'), 'w') as xyzFile:
            xyzFile.write(input_string)

    def getBonds(self, qmMolecule):
        bondList = []
        for atom in qmMolecule.molecule.atoms:
            for bond in atom.bonds.values():
                bondList.append(bond)
        bonds = list(set(bondList))
        bondDict = {}
        for bond in bonds:
            if bond.isSingle():
                if bond.atom1.symbol=='C' and bond.atom2.symbol=='C':
                    bondType = 'C-C'
                elif (bond.atom1.symbol=='H' and bond.atom2.symbol=='H'):
                    bondType = 'H-H'
                elif (bond.atom1.symbol=='C' and bond.atom2.symbol=='H') or (bond.atom1.symbol=='H' and bond.atom2.symbol=='C'):
                    bondType = 'C-H'
                elif (bond.atom1.symbol=='O' and bond.atom2.symbol=='O'):
                    bondType = 'O-O'
                elif (bond.atom1.symbol=='C' and bond.atom2.symbol=='O') or (bond.atom1.symbol=='O' and bond.atom2.symbol=='C'):
                    bondType = 'C-O'
                elif (bond.atom1.symbol=='H' and bond.atom2.symbol=='O') or (bond.atom1.symbol=='O' and bond.atom2.symbol=='H'):
                    bondType = 'O-H'
                elif bond.atom1.symbol=='N' and bond.atom2.symbol=='N':
                    bondType = 'N-N'
                elif (bond.atom1.symbol=='C' and bond.atom2.symbol=='N') or (bond.atom1.symbol=='N' and bond.atom2.symbol=='C'):
                    bondType = 'N-C'
                elif (bond.atom1.symbol=='O' and bond.atom2.symbol=='N') or (bond.atom1.symbol=='N' and bond.atom2.symbol=='O'):
                    bondType = 'N-O'
                elif (bond.atom1.symbol=='H' and bond.atom2.symbol=='N') or (bond.atom1.symbol=='N' and bond.atom2.symbol=='H'):
                    bondType = 'N-H'
                elif bond.atom1.symbol=='S' and bond.atom2.symbol=='S':
                    bondType = 'S-S'
                elif (bond.atom1.symbol=='H' and bond.atom2.symbol=='S') or (bond.atom1.symbol=='S' and bond.atom2.symbol=='H'):
                    bondType = 'S-H'
            elif bond.isDouble:
                if bond.atom1.symbol=='C' and bond.atom2.symbol=='C':
                    bondType = 'C=C'
                elif (bond.atom1.symbol=='O' and bond.atom2.symbol=='O'):
                    bondType = 'O=O'
                elif (bond.atom1.symbol=='C' and bond.atom2.symbol=='O') or (bond.atom1.symbol=='O' and bond.atom2.symbol=='C'):
                    bondType = 'C=O'
                elif bond.atom1.symbol=='N' and bond.atom2.symbol=='N':
                    bondType = 'N=N'
                elif (bond.atom1.symbol=='C' and bond.atom2.symbol=='N') or (bond.atom1.symbol=='N' and bond.atom2.symbol=='C'):
                    bondType = 'N=C'
                elif (bond.atom1.symbol=='O' and bond.atom2.symbol=='N') or (bond.atom1.symbol=='N' and bond.atom2.symbol=='O'):
                    bondType = 'N=O'
                elif (bond.atom1.symbol=='O' and bond.atom2.symbol=='S') or (bond.atom1.symbol=='S' and bond.atom2.symbol=='O'):
                    bondType = 'S=O'
            elif bond.isTriple:
                if bond.atom1.symbol=='C' and bond.atom2.symbol=='C':
                    bondType = 'C#C'
                elif bond.atom1.symbol=='N' and bond.atom2.symbol=='N':
                    bondType = 'N#N'
                elif (bond.atom1.symbol=='C' and bond.atom2.symbol=='N') or (bond.atom1.symbol=='N' and bond.atom2.symbol=='C'):
                    bondType = 'N#C'
            try:
                bondDict[bondType] += 1
            except KeyError:
                bondDict[bondType] = 1

        return bondDict

    def writeCanThermStatMech(self, allAtomTypes, bondDict, multiplicity, path, symmetry=1):
        """
        Write the file required by cantherm to do the statmech calculation for the molecule
        """
        output = ['#!/usr/bin/env python', '# -*- coding: utf-8 -*-', '', 'atoms = {']

        # add the atoms
        atomTypes = set(list(allAtomTypes))
        for atom in atomTypes:
            output.append("    '{0}': {1},".format(atom, allAtomTypes.count(atom)))
        output = output + ['}', '']

        if bondDict!={}:
            output.append('bonds = {')
            for bondType, num in bondDict.iteritems():
                output.append("    '{0}': {1},".format(bondType, num))
            output.append('}')
        else:
            output.append('bonds = {}')

        # add the bonds
        output = output + ['', 'linear = False', '', 'externalSymmetry = {0}'.format(symmetry), '', 'spinMultiplicity = {0}'.format(multiplicity), '', 'opticalIsomers = 1', '']

        # Energy ~ get it from the QM output file
        if self.method=='b3lyp':
            output = output + ['energy = {', "    'DFT_G03_b3lyp': GaussianLog('{0}'),".format(path.split('/')[-1]), '}', '']
        elif self.method=='m062x':
            output = output + ['energy = {', "    'M062X/MG3S': GaussianLog('{0}'),".format(path.split('/')[-1]), '}', '']

        # Geometry - get it from the QM output file
        output = output + ["geometry = GaussianLog('{0}')".format(path.split('/')[-1]), '']

        # Frequencies - get them from the QM output file
        output = output + ["frequencies = GaussianLog('{0}')".format(path.split('/')[-1]), '']

        # Rotors -  no rotors yet
        output = output + ['rotors = []', '']

        input_string = '\n'.join(output)
        with open(path.rsplit('.',1)[0]+'.py', 'w') as statMechFile:
            statMechFile.write(input_string)

    def writeCanThermInput(self, reactants, products, filePath, fileStore, scratchDirectory):
        """
        Write the CanTherm input file
        """
        output = ['#!/usr/bin/env python', '# -*- coding: utf-8 -*-', '']

        if self.method=='b3lyp':
            output.append('modelChemistry = "DFT_G03_b3lyp"')
            output.append('frequencyScaleFactor = 0.964')
        elif self.method=='m062x':
            output.append('modelChemistry = "M062X/MG3S"')
            output.append('frequencyScaleFactor = 0.982') # source - http://comp.chem.umn.edu/freqscale/version3b1.htm
        output.append('useHinderedRotors = False')
        output.append('useBondCorrections = False\n')

        if not fileStore.startswith('/'):
            fileStore = os.path.abspath(fileStore)
            
        speciesList = [] # Check if species path already specified. If so, DON'T put it again.
        for reactant in reactants:
            if reactant.uniqueID not in speciesList:
                reactant.settings.fileStore = fileStore
                reactant.settings.scratchDirectory = scratchDirectory
                reactant.checkPaths()
                output.append("species('{0}', '{1}')".format(reactant.uniqueID, reactant.getFilePath('.py')))
                speciesList.append(reactant.uniqueID)
        for product in products:
            if product.uniqueID not in speciesList:
                product.settings.fileStore = fileStore
                product.settings.scratchDirectory = scratchDirectory
                product.checkPaths()
                output.append("species('{0}', '{1}')".format(product.uniqueID, product.getFilePath('.py')))
                speciesList.append(product.uniqueID)

        output.append("transitionState('TS', '{0}')".format(self.getFilePath('.py')))
        output.append('')
        output.append('reaction(')
        output.append("    label = '{0}',".format(self.uniqueID))
        if len(reactants)==2:
            output.append("    reactants = ['{0}', '{1}'],".format(reactants[0].uniqueID, reactants[1].uniqueID))
        else:
            output.append("    reactants = ['{0}'],".format(reactants[0].uniqueID))

        if len(products)==2:
            output.append("    products = ['{0}', '{1}'],".format(products[0].uniqueID, products[1].uniqueID))
        else:
            output.append("    products = ['{0}'],".format(products[0].uniqueID))
        output.append("    transitionState = 'TS',")
        output.append("    tunneling = 'Eckart',\n)\n")
        output.append("statmech('TS')\nkinetics('{0}')\n".format(self.uniqueID))

        input_string = '\n'.join(output)
        with open(filePath, 'w') as canThermInp:
            canThermInp.write(input_string)

    def calculateQMData(self, moleculeList, fileStore, scratchDirectory):
        """
        If the transition state is found, optimize reactant and product geometries for use in
        TST calculations.
        """
        molecules = []
        for molecule in moleculeList:
            if isinstance(molecule, Species):
                molecule = molecule.molecule[0]
            qmMolecule = self.getQMMolecule(molecule)
            qmMolecule.settings.fileStore = fileStore
            qmMolecule.settings.scratchDirectory = scratchDirectory
            qmMolecule.checkPaths()
            qmMolecule.qmData = qmMolecule.generateQMData()
            if qmMolecule.qmData:
                qmMolecule.determinePointGroup()
                allAtoms = []
                for atom in qmMolecule.molecule.atoms:
                    allAtoms.append(atom.symbol)
                bondDict = self.getBonds(qmMolecule)
                self.writeCanThermStatMech(allAtoms, bondDict, qmMolecule.molecule.multiplicity, qmMolecule.outputFilePath, symmetry=1)#qmMolecule.pointGroup.symmetryNumber)
                molecules.append(qmMolecule)
                # log = GaussianLog(qmMolecule.outputFilePath)
                # species = Species(label=qmMolecule.molecule.toSMILES(), conformer=log.loadConformer(), molecule=[molecule])
                # molecules.append(species)
            else:
                raise Exception('The reactant or product geometry did not optimize. Cannot calculate the kinetics.')
        return molecules

    def generateKineticData(self):
        self.initialize()
        # provides transitionstate geometry
        fileStore = self.settings.fileStore #  To ensure all files are found in the same base directory
        scratchDirectory = self.settings.scratchDirectory #  To ensure all files are found in the same base directory
        tsFound, notes = self.generateTSGeometryDirectGuess()

        self.settings.fileStore = fileStore
        with open(os.path.join(self.fileStore, 'error.txt'), 'w') as errorFile:
            errorFile.write(notes)

        if not tsFound:
            # Return the reaction without the kinetics included. Fall back on group additivity.
            return self.reaction

        # cantopt = ['CC(C)C([O])(OO)C(C)C', 'C[C](CC(=O)C(C)(C)OO)OO', 'C[C](OO)C(=O)C(C)(C)OO', 'CC(C)(C)[CH]OO', 'OO[C]1CCCCC1', '[CH2]OO', '[O]Cc1ccccc1'] # A local minimum for the following geometries has not been found at B3LYP/6-31+G(d,p)
        cantopt = [] # They have been found at m062x
        for mol in self.reaction.reactants:
            if isinstance(mol, Species):
                mol = mol.molecule[0]
            if mol.toSMILES() in cantopt:
                print 'Cannot optimize geometry for {0}'.format(mol.toSMILES())
                return self.reaction
        for mol in self.reaction.products:
            if isinstance(mol, Species):
                mol = mol.molecule[0]
            if mol.toSMILES() in cantopt:
                print 'Cannot optimize geometry for {0}'.format(mol.toSMILES())
                return self.reaction

        reactants = self.calculateQMData(self.reaction.reactants, fileStore, scratchDirectory)
        products = self.calculateQMData(self.reaction.products, fileStore, scratchDirectory)

        allAtoms = []
        multiplicity=1
        for molecule in self.reaction.reactants:
            if isinstance(molecule, Species):
                molecule = molecule.molecule[0]
            multiplicity += molecule.getRadicalCount()
            for atom in molecule.atoms:
                allAtoms.append(atom.symbol)

        self.writeCanThermStatMech(allAtoms, {}, multiplicity, self.outputFilePath)
        canThermFilePath = os.path.join(self.fileStore, 'input.py')
        self.writeCanThermInput(reactants, products, canThermFilePath, fileStore, scratchDirectory)
        canThermJob = CanTherm()
        canThermJob.outputDirectory = self.fileStore
        canThermJob.inputFile = canThermFilePath
        canThermJob.plot = False
        try:
            canThermJob.execute()
    	    jobResult = "Successful kinetics calculation in CanTherm."
    	    for job in canThermJob.jobList:
    	        if isinstance(job, KineticsJob):
    		    # Return the reaction with the kinetics.
    	            return job.reaction
    	except ValueError, e:
    	    jobResult = e

    	logging.info(jobResult)
        return self.reaction
        #
        # if len(reactants)==len(self.reaction.reactants) and len(products)==len(self.reaction.products):
        #     #self.determinePointGroup()
        #     tsLog = GaussianLog(self.outputFilePath)
        #     self.reaction.transitionState = TransitionState(label=self.uniqueID + 'TS', conformer=tsLog.loadConformer(), frequency=(tsLog.loadNegativeFrequency(), 'cm^-1'), tunneling=Wigner(frequency=None))
        #
        #     self.reaction.reactants = reactants
        #     self.reaction.products = products
        #
        #
        #     kineticsJob = KineticsJob(self.reaction)
        #     kineticsJob.generateKinetics()
        #
        #     """
        #     What do I do with it? For now just save it.
        #     Various parameters are not considered in the calculations so far e.g. symmetry.
        #     This is just a crude calculation, calculating the partition functions
        #     from the molecular properties and plugging them through the equation.
        #     """
        #     kineticsJob.save(self.getFilePath('.kinetics'))
        #     # return self.reaction.kinetics

    def determinePointGroup(self):
        """
        Determine point group using the SYMMETRY Program

        Stores the resulting :class:`PointGroup` in self.pointGroup
        """
        assert self.qmData, "Need QM Data first in order to calculate point group."
        pgc = symmetry.PointGroupCalculator(self.settings, self.uniqueID, self.qmData)
        self.pointGroup = pgc.calculate()
