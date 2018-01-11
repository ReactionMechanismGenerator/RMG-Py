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

import rmgpy
from rmgpy.data.kinetics.transitionstates import TransitionStates
from rmgpy.molecule import Molecule, Atom, getElement
from rmgpy.species import Species, TransitionState
from rmgpy.reaction import Reaction
from rmgpy.kinetics import Arrhenius, Eckart
from rmgpy.statmech import Conformer, IdealGasTranslation, NonlinearRotor, HarmonicOscillator, LinearRotor
from molecule import QMMolecule, Geometry
from rmgpy.cantherm.main import CanTherm
from rmgpy.cantherm.kinetics import KineticsJob
import qmdata
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

def loadTSDataFile(filePath):
    """
    Load the specified transition state data file and return the dictionary of its contents.

    Returns `None` if the file is invalid or missing.

    Checks that the returned dictionary contains at least rxnLabel, method, qmData.
    """
    
    if not os.path.exists(filePath):
        return None
    try:
        with open(filePath) as resultFile:
            logging.info('Reading existing ts file {0}'.format(filePath))
            global_context = { '__builtins__': None }
            local_context = {
                '__builtins__': None,
                'True': True,
                'False': False,
                'QMData': qmdata.QMData,
                'array': numpy.array,
                'int32': numpy.int32,
            }
            exec resultFile in global_context, local_context
    except IOError, e:
        logging.info("Couldn't read ts file {0}".format(filePath))
        return None
    except (NameError, TypeError, SyntaxError), e:
        logging.error('The ts file "{0}" was invalid:'.format(filePath))
        logging.exception(e)
        return None
    if not 'rxnLabel' in local_context:
        logging.error('The ts file "{0}" did not contain a rxnLabel.'.format(filePath))
        return None
    if not 'method' in local_context:
        logging.error('The ts file "{0}" did not contain a method.'.format(filePath))
        return None
    if not 'qmData' in local_context:
        logging.error('The ts file "{0}" did not contain thermoData.'.format(filePath))
        return None
    return local_context

def loadKineticsDataFile(filePath):
    """
    Load the specified kinetic data file and return the dictionary of its contents.

    Returns `None` if the file is invalid or missing.

    Checks that the returned dictionary contains at least method, Reaction.
    """
    
    if not os.path.exists(filePath):
        return None
    try:
        with open(filePath) as resultFile:
            logging.info('Reading existing kinetics file {0}'.format(filePath))
            global_context = { '__builtins__': None }
            local_context = {
                '__builtins__': None,
                'True': True,
                'False': False,
                'Reaction': Reaction,
                'Species': Species,
                'TransitionState': TransitionState,
                'Arrhenius': Arrhenius,
                'Eckart': Eckart,
                'Conformer': Conformer,
                'IdealGasTranslation': IdealGasTranslation,
                'NonlinearRotor': NonlinearRotor,
                'HarmonicOscillator': HarmonicOscillator,
                'LinearRotor': LinearRotor,
                'array': numpy.array,
                'int32': numpy.int32,
            }
            exec resultFile in global_context, local_context
    except IOError, e:
        logging.error("Couldn't read kinetics file {0}".format(filePath))
        return None
    except (NameError, TypeError, SyntaxError), e:
        logging.error('The kinetics file "{0}" was invalid:'.format(filePath))
        logging.exception(e)
        return None
    if not 'method' in local_context:
        logging.error('The kinetics file "{0}" did not contain a method.'.format(filePath))
        return None
    if not 'Reaction' in local_context:
        logging.error('The kinetics file "{0}" did not contain a reaction.'.format(filePath))
        return None
    return local_context

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

    def getFilePath(self, extension, scratch=True):
        """
        Returns the path to the file with the given extension.

        The provided extension should include the leading dot.

        If called with `scratch=False` then it will be in the `fileStore` directory,
        else `scratch=True` is assumed and it will be in the `scratchDirectory` directory.
        """
        return os.path.join(
            self.settings.scratchDirectory if scratch else self.settings.fileStore,
            self.uniqueID + extension
            )

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
        
    @property
    def duplicateFam(self):
        """Get reacton families that should are it's own reverse."""
        duplicateFam = {
            'H_Abstraction': True,
            'R_Addition_MultipleBond': False,
            'intra_H_migration': True,
            'Disproportionation': False,
        }
        return duplicateFam
    
    @property
    def getTSFilePath(self):
        "Returns the path the transition state data file."
        return self.getFilePath('.ts', scratch=False)
    
    @property
    def getCanThermFilePath(self):
        "Returns the path the cantherm file."
        return self.getFilePath('.canth.py')
        
    @property
    def getKineticsFilePath(self):
        "Returns the path the kinetics data file."
        return self.getFilePath('.kinetics', scratch=False)

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
        logging.info("Generating bounds matrix for {}".format(molecule.toSMILES()))
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
                        logging.info("Changing lower limit {0} to {1}".format(bm[i, j], maxLij))
                        bm[i,j] = maxLij

        return bm

    def getLabels(self, reactant):
        """
        Creates the list of sorting labels for the reacting atoms. These labels are
        also put intoa tuple for the atomMatch weighting for RDKit. The weighting tells
        RDKit to place greater importance in maintaining these distance limits when
        generating conformers.
        """
        if self.reaction.family.lower() in ['h_abstraction', 'r_addition_multiplebond', 'intra_h_migration']:
            lbl1 = reactant.getLabeledAtom('*1').sortingLabel
            lbl2 = reactant.getLabeledAtom('*2').sortingLabel
            lbl3 = reactant.getLabeledAtom('*3').sortingLabel
        elif self.reaction.family.lower() in ['disproportionation']:
            lbl1 = reactant.getLabeledAtom('*2').sortingLabel
            lbl2 = reactant.getLabeledAtom('*4').sortingLabel
            lbl3 = reactant.getLabeledAtom('*1').sortingLabel

        labels = [lbl1, lbl2, lbl3]
        atomMatch = ((lbl1,),(lbl2,),(lbl3,))

        return labels, atomMatch

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

        bm = self.bmPreEdit(bm, sect)

        return bm

    def setupMolecules(self):
        """
        Setup the reactant and product molecules for the transition state calculations.
        If there are 2 reactants and/or products, they are merged. This also handles
        species as well as molecules, but returns the reactant and product as merged molecules.
        """
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

    def optimizeTS(self, labels):
        """
        Conduct the optimization step of the transition state search.
        """
        logging.info("Starting optimization steps of the TS search.")
        if os.path.exists(self.outputFilePath):
            completeOutputFileExists, _ = self.checkComplete(self.outputFilePath)
            if not completeOutputFileExists:
                # Looks like calculation was interrupted, rather than failed.
                logging.info("Output file {} exists but looks incomplete, so deleting it to try again.".format(self.outputFilePath))
                os.remove(self.outputFilePath)
        else:
            completeOutputFileExists = False

        if completeOutputFileExists:
            logging.info("Output file {} exists and looks complete. Trying that.".format(self.outputFilePath))
            converged, internalCoord = self.verifyOutputFile()
        else:
            optEst = self.optEstimate(labels)
            optRC = self.optRxnCenter(labels)

            logging.info("Optimizing TS attempt 1")
            self.createInputFile(1, optEst=optRC)
            converged, internalCoord = self.run()
            shutil.copy(self.outputFilePath, self.outputFilePath+'.TS1.log')
        
            if os.path.exists(self.ircOutputFilePath):# Remove it
                os.remove(self.ircOutputFilePath)

        if internalCoord and not converged:
            logging.info("Optimizing TS attempt 2")
            self.createInputFile(2)
            converged = self.run()
            shutil.copy(self.outputFilePath, self.outputFilePath+'.TS2.log')
            
        if not converged:
            # Check for convergence failures
            complete, convergenceFailure = self.checkComplete(self.outputFilePath)
            if convergenceFailure:
                # Rerun the calculation with `scf=qc`
                logging.info("Optimizing TS attempt 1 but with quadratic scf")
                self.createInputFile(1, optEst=optRC, scf=True)
                converged, internalCoord = self.run()
                if internalCoord:
                    logging.info("Optimizing TS attempt 2 but with quadratic scf")
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

    def tsSearch(self, labels):
        """
        Once the transition state estimate is made, this runs the optimization and the
        path analysis calculation. The ts estimate is derived from the group additive prediction.
        """
        successfulTS = self.optimizeTS(labels)
        if not successfulTS:
            return successfulTS

        validTS = self.validateTS()
        if validTS:
            self.writeRxnOutputFile(labels)
            self.saveTSData()

        return validTS
            
    def generateTSGeometryDirectGuess(self):
        """
        Generate a transition state geometry, using the direct guess (group additive) method.

        Returns (success, notes) where success is a True if it worked, else False,
        and notes is a string describing what happened.
        """
        logging.info("Generating a TS geometry via the direct guess method")
        def getDistance(coordinates1, coordinates2):
            """
            Return the square of the distance (in Angstrom) between the two atoms.
            """
            diff = (coordinates1.coords - coordinates2.coords)
            return math.sqrt(sum(diff * diff))
            
        reactant, product = self.setupMolecules()

        tsRDMol, tsBM, self.reactantGeom = self.generateBoundsMatrix(reactant)

        self.reactantGeom.uniqueID = self.uniqueID

        labels, atomMatch = self.getLabels(reactant)
        
        if os.path.exists(self.getTSFilePath):
            logging.info("TS file {} exists. Loading it to save time.".format(self.getTSFilePath))
            self.loadTSData()
            return True
        
        tsBM = self.editMatrix(reactant, tsBM, labels)
        atoms = len(reactant.atoms)
        if atoms>3:
            distGeomAttempts = 15*(atoms-3) # number of conformers embedded from the bounds matrix
        else:
            distGeomAttempts = 15

        setBM = rdkit.DistanceGeometry.DoTriangleSmoothing(tsBM)

        if not setBM:
            logging.error("Couldn't create smoothed bounds matrix")
            return False

        for i in range(len(tsBM)):
            for j in range(i,len(tsBM)):
                if tsBM[j,i] > tsBM[i,j]:
                        logging.error("BOUNDS MATRIX FLAWED {0}>{1}".format(tsBM[j, i], tsBM[i, j]))

        self.reactantGeom.rd_embed(tsRDMol, distGeomAttempts, bm=tsBM, match=atomMatch)
        atomSymbols, atomCoords = self.reactantGeom.parseMOL(self.reactantGeom.getRefinedMolFilePath())
        logging.info("TS estimate made. About to try the search...")
        check =  self.tsSearch(labels)
        
        return check
    
    def saveTSData(self):
        """
        Save the generated TS data.
        """
        logging.info("Saving TS result file {}".format(self.getTSFilePath))
        with open(self.getTSFilePath, 'w') as resultFile:
            resultFile.write('rxnLabel = "{0!s}"\n'.format(self.uniqueID))
            resultFile.write('method = "{0!s}"\n'.format(self.method))
            resultFile.write("qmData = {0!r}\n".format(self.parse()))
    
    def loadTSData(self):
        """
        Try loading TS data from a previous run.
        """
        filePath = self.getTSFilePath
        
        local_context = loadTSDataFile(filePath)
        if local_context is None:
            # file does not exist or is invalid
            logging.warning("TS data file {} does not exist or is invalid".format(filePath))
            return None
        
        if local_context['rxnLabel'] != self.uniqueID:
            if local_context['rxnLabel'] != self.revID:
                logging.error('The rxnLabel in the ts file {0} did not match the current reaction {1}'.format(filePath,self.uniqueID))
                return None
        if local_context['method'] != self.method:
            logging.error('The method in the ts file {0} did not match the current method {1}'.format(filePath,self.method))
            return None

        self.qmData = local_context['qmData']
        
        return self.qmData
    
    def saveKineticsData(self, reaction):
        """
        Save the calculated kinetics. `reaction` is a CanTherm reaction object that
        should include the molecular parameters.
        """
        logging.info("Saving kinetics data file {}".format(self.getKineticsFilePath))
        with open(self.getKineticsFilePath, 'w') as resultFile:
            resultFile.write('method = "{0!s}"\n'.format(self.method))
            resultFile.write('reaction = {0!r}\n'.format(reaction))
    
    def loadKineticsData(self):
        """
        Try loading kinetics from a previous run.
        """
        filePath = self.getKineticsFilePath
        
        local_context = loadKineticsDataFile(filePath)
        if local_context is None:
            # file does not exist or is invalid
            logging.info("Kinetics results file {} does not exist or is invalid".format(filePath))
            return None
        
        if local_context['method'] != self.method:
            logging.error('The method in the ts file {0} did not match the current method {1}'.format(filePath,self.method))
            return None
        
        if local_context['reaction'].label != self.uniqueID:
            if local_context['reaction'].label != self.revID:
                logging.error('The reaction in the ts file {0} did not match the current reaction {1}'.format(filePath,self.uniqueID))
                return None

        self.reaction = local_context['reaction']
        
        return self.reaction
        
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

    def writeCanThermInput(self, reactants, products):
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

        if not self.fileStore.startswith('/'):
            fileStore = os.path.abspath(self.fileStore)
            
        speciesList = [] # Check if species path already specified. If so, DON'T put it again.
        for reactant in reactants:
            if reactant.uniqueID not in speciesList:
                output.append("species('{0}', '{1}')".format(reactant.uniqueID, reactant.getFilePath('.py')))
                speciesList.append(reactant.uniqueID)
        for product in products:
            if product.uniqueID not in speciesList:
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
        with open(self.getCanThermFilePath, 'w') as canThermInp:
            canThermInp.write(input_string)

    def calculateQMData(self, moleculeList):
        """
        If the transition state is found, optimize reactant and product geometries for use in
        TST calculations.
        """
        molecules = []
        for molecule in moleculeList:
            if isinstance(molecule, Species):
                molecule = molecule.molecule[0]
            qmMolecule = self.getQMMolecule(molecule)
            qmMolecule.qmData = qmMolecule.generateQMData()
            if qmMolecule.qmData:
                qmMolecule.determinePointGroup()
                allAtoms = []
                for atom in qmMolecule.molecule.atoms:
                    allAtoms.append(atom.symbol)
                bondDict = self.getBonds(qmMolecule)
                self.writeCanThermStatMech(allAtoms, bondDict, qmMolecule.molecule.multiplicity, qmMolecule.outputFilePath, symmetry=1)#qmMolecule.pointGroup.symmetryNumber)
                molecules.append(qmMolecule)
            else:
                raise Exception('The reactant or product geometry did not optimize. Cannot calculate the kinetics.')
        return molecules

    def generateKineticData(self):
        """
        Generate Kinetic Data via a QM calc.

        Returns :class:`Reaction`. :class:`Reaction` has a kinetics attribute, which is None if the calculation fails.
        """
        if self.loadKineticsData():
            return self.reaction
        
        self.initialize()
        
        tsFound = self.generateTSGeometryDirectGuess()

        if not tsFound:
            # Return the reaction without the kinetics included. Fall back on group additivity.
            logging.warning("Couldn't find transition state. Not using TST")
            return self.reaction

        reactants = self.calculateQMData(self.reaction.reactants)
        products = self.calculateQMData(self.reaction.products)

        allAtoms = []
        multiplicity=1
        for molecule in self.reaction.reactants:
            if isinstance(molecule, Species):
                molecule = molecule.molecule[0]
            multiplicity += molecule.getRadicalCount()
            for atom in molecule.atoms:
                allAtoms.append(atom.symbol)

        self.writeCanThermStatMech(allAtoms, {}, multiplicity, self.outputFilePath)
        self.writeCanThermInput(reactants, products)
        canThermJob = CanTherm()
        canThermJob.outputDirectory = self.getFilePath('')
        if not os.path.exists(canThermJob.outputDirectory):
            os.makedirs(canThermJob.outputDirectory)
        canThermJob.inputFile = self.getCanThermFilePath
        canThermJob.plot = False
        
        try:
            canThermJob.execute()
            jobResult = "Successful kinetics calculation in CanTherm."
            for job in canThermJob.jobList:
                if isinstance(job, KineticsJob):
                    # Return the reaction with the kinetics.
                    self.saveKineticsData(job.reaction)
                    return job.reaction
        except ValueError, e:
            jobResult = e

        logging.info(jobResult)
        return self.reaction

    def determinePointGroup(self):
        """
        Determine point group using the SYMMETRY Program

        Stores the resulting :class:`PointGroup` in self.pointGroup
        """
        assert self.qmData, "Need QM Data first in order to calculate point group."
        pgc = symmetry.PointGroupCalculator(self.settings, self.uniqueID, self.qmData)
        self.pointGroup = pgc.calculate()
