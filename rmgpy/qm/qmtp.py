"""
Main entry point to access the on-the-fly estimation of 
thermochemical properties of species using quantum chemical packages such 
as G03, (Open)Mopac and MM4.

This module contains the driver class called QMTP with its principal method
 .generateQMThermoData that generates a ThermoData object.
 
A molfile class is a rather artificial and probably redundant class from RMG-Java
that collects information on paths to the 3D mol files.

ThreeDMolFileCreator is a wrapper class that:
-converts a RMG-Py molecule into a 2D mol file using OpenBabel
-generates 3D coordinates using RDKit 
"""

import os

from rdkit import Chem
from rdkit.Chem import AllChem

import logging
from rmgpy.molecule import Molecule
import qmverifier
import parsers
import jobs
import inputwriters


class molFile:
    """
    molFile class to store and perform operations on mol files
    may eventually want to have a super class threeDGeom and say "public class molFile extends threeDGeom {"
    """
    
    def __init__(self, molecule, name='', directory='', InChIAug = ''):
        #construct a molFile object while writing a file with text determined by chemGraph
        self.name = name
        self.directory = directory
        self.molecule = molecule
        self.InChIAug = InChIAug#Augmented InChI
        self.path = os.path.join(self.directory,self.name +'.mol')
        self.crudepath = os.path.join(self.directory, self.name+'.cmol')
     
    
class ThreeDMolFileCreator:
    
    def __init__(self, name, directory, molecule):
        self.name = name
        self.directory = directory
        self.molecule = molecule
        
        self.mol2D = name + '.mol2D'
        self.mol3D = name + '.mol'
        self.mol3Dcrude = name + '.cmol'
        
        self.dir2D = os.path.join(directory,'2DMolfiles/')
        self.dir3D = os.path.join(directory,'3DMolfiles/')
        
        for directory in (self.dir2D, self.dir3D):
            if not os.path.exists(directory):
                logging.info("Creating directory %s for mol files."%os.path.abspath(directory))
                os.makedirs(directory)

    def createRDMol(self):
        """
        Import rmg molecule and create rdkit molecule with the same atom labeling.
        """
        rdAtomIdx = {} # dictionary of rdkit atom indices
        # Initialize a blank Editable molecule and add all the atoms from RMG molecule
        rdMol = AllChem.rdchem.EditableMol(AllChem.rdchem.Mol())
        for index, atom in enumerate(self.molecule.vertices):
            rdAtom = AllChem.rdchem.Atom(atom.element.symbol)
            rdAtom.SetNumRadicalElectrons(atom.radicalElectrons)
            rdMol.AddAtom(rdAtom)
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
                    rdMol.AddBond(index1, index2, rdBond)

        # Make editable mol into a mol and rectify the molecule
        rdMol = rdMol.GetMol()
        Chem.SanitizeMol(rdMol)

        return rdMol, rdAtomIdx
        
    def embed3D(self, numConfAttempts):
        """
        #embed a molecule in 3D, using RDKit
        """
        m, mIdx = self.createRDMol()

        AllChem.EmbedMultipleConfs(m, numConfAttempts,randomSeed=1)
        
        m2crude = Chem.Mol(m.ToBinary()) #make a copy of the (crude) coordinates via ToBinary
        energy=0.0
        minEid=0;
        lowestE=9.999999e99;#start with a very high number, which would never be reached
        
        for i in range(m.GetNumConformers()):
            AllChem.UFFOptimizeMolecule(m,confId=i)
            energy=AllChem.UFFGetMoleculeForceField(m,confId=i).CalcEnergy()
            if energy < lowestE:
                minEid = i
                lowestE = energy
        with open(os.path.join(self.dir3D, self.mol3D), 'w') as out3D:
            out3D.write(Chem.MolToMolBlock(m,confId=minEid))
        
        with open(os.path.join(self.dir3D, self.mol3Dcrude), 'w') as out3Dcrude:
            out3Dcrude.write(Chem.MolToMolBlock(m2crude,confId=minEid))        

        #construct molFile pointer to new file (name will be same as 2D mol file
        return molFile(self.molecule, self.name, self.dir3D)
    
    def create(self):
        #2. convert from 2D to 3D using RDKit if the 2D molfile is for a molecule with 2 or more atoms
        atoms = len(self.molecule.atoms)
        distGeomAttempts=1
        if atoms > 3:#this check prevents the number of attempts from being negative
            distGeomAttempts = 5*(atoms-3) #number of conformer attempts is just a linear scaling with molecule size, due to time considerations in practice, it is probably more like 3^(n-3) or something like that
            
        molfile = self.embed3D(distGeomAttempts)
        
        return molfile

class QMTP:
    """
    Quantum mechanics thermo property estimator analog of GATP
    """
    
    #static fields
    qmfolder= "QMfiles"
    
    
    """
         * the qmprogram can be 
     * "mopac", 
     * "gaussian03",
     *  "both" (MOPAC and Gaussian), 
     *  "mm4",
     *  "mm4hr"
     """

    
    usePolar = False#use polar keyword in MOPAC
    useCanTherm = True#whether to use CanTherm in MM4 cases for interpreting output via force-constant matrix this will hopefully avoid zero frequency issues
    
    def __init__(self, qmprogram, qmMethod = '', useHindRot = False):
        self.qmprogram = qmprogram
        #set global attribute qmMethod
        if self.qmprogram == "mm4" or self.qmprogram == "mm4hr":
            self.qmMethod = "mm4"
            if self.qmprogram == "mm4hr":#whether to use HinderedRotor scans with MM4 (requires useCanTherm=true)
                self.qmMethod = "mm4hr"
                self.useHindRot=True
            
        else:
            #may eventually want to pass this to various functions to choose which "sub-function" to call
            self.qmMethod="pm3"
        

        self.mapMaxAttemptNumber = { }
        self.mapMaxAttemptNumber["gaussian03"]= 36
        self.mapMaxAttemptNumber["mopac"] = 10
        self.mapMaxAttemptNumber["both"] = 10
        self.mapMaxAttemptNumber["mm4"] = 4
        self.mapMaxAttemptNumber["mm4hr"] = 4

    def generateIdentifiers(self, molecule):
        """
            #determine the QM filename (element 0) and augmented InChI (element 1) for a ChemGraph
            #QM filename is InChIKey appended with mult3, mult4, mult5, or mult6 for multiplicities of 3 or higher
            #augmented InChI is InChI appended with /mult3, /mult4, /mult5, or /mult6 for multiplicities of 3 or higher
        """ 
        #inchikey_mod = molecule.getModifiedInChIKeyAnew()
        inchikey_mod = molecule.toAugmentedInChIKey()#need to generate InChI and key anew because ChemGraph may have changed (in particular, adding/removing hydrogens in HBI process)
        inchi_mod = molecule.toAugmentedInChI()#
        return inchikey_mod, inchi_mod
 
    def parseOutput(self, molfile):
        """
        wrapper method for parser types
        """

        if self.qmMethod == "pm3" :
                        
            if  self.qmprogram == "mopac" or self.qmprogram == "both" :
                parser = parsers.MOPACPM3Parser(molfile, self)
                result = parser.parse()
                result.comment = result.comment +'MOPAC PM3 calculation'
                logging.info("Thermo for " + molfile.name + ": "+ result.__repr__())#print result, at least for debugging purposes
                return result
            
        else:
            logging.critical("Unexpected situation in QMTP thermo estimation"   )    

        return result

    
    def generateQMThermoData(self, molecule):
        """
              /**
        * #if there is no data in the libraries, calculate the result based on QM or MM calculations
        *  the below steps will be generalized later to allow for other quantum mechanics packages, etc.
        *  
        *  Several steps can be identified:
        *  <LI> InChI generation of ChemGraph
        *  <LI> Verification whether this species has already been processed (and parse it if so)
        *  <LI> Generating 3D coords
        *  <LI> Creating QM input file and Feeding species to QM Program
        *  <LI> QM Program Output File Parsing
        * @param molecule
        * @return
        */
        """
        name, InChIaug = self.generateIdentifiers(molecule)#determine the filename (InChIKey) and InChI with appended info for triplets, etc.
        
        #check for existing hold file before starting calculations (to avoid the possibility of interference with other jobs using the same QMfiles folder)
        while os.path.exists(os.path.join(QMTP.qmfolder,name+'.hold')):#can lead to race condition!
            logging.info("Existence of hold file for "+name+" suggests that another RMG process is currently running calculations on this molecule; waiting for other RMG process to finish; will check again in 60 seconds...")
            import time
            time.sleep(60) # delays for 60 seconds
       
        #verify whether a succesful QM results exists for this particular species:
        molfile = molFile(Molecule(), name, InChIaug, QMTP.qmfolder)
        verifier = qmverifier.QMVerifier(molfile)
        verifier.verify()   
         
        #if a succesful job exists (by one of the QM Programs), you can readily parse it.
        if verifier.succesfulJobExists():
            result = self.parseOutput(name, molecule)
            return result
        else:#no successful result exists, we have to calculate from zero
            molfile = self.generate3DCoords(molecule, name)
            molfile.InChIAug = InChIaug
            multiplicity = sum([i.radicalElectrons for i in molecule.atoms]) + 1
            attemptNumber = 1
            success = False
            maxAttemptNumber = self.mapMaxAttemptNumber[self.qmprogram]
            while not success and (attemptNumber <= maxAttemptNumber):
                self.createQMInput(molfile, attemptNumber, multiplicity)
                success = self.runQM(molfile)
                if success:
                    logging.info('Attempt {0} on species {1} succeeded.'.format(attemptNumber, InChIaug))
                    """
                    TODO Rotor Scan not yet implemented here.
                    """
                else:
                    if attemptNumber == maxAttemptNumber:
                        logging.info('Last attempt on species {0} failed.'.format(InChIaug))
            result = self.parseOutput(molfile)
            return result

    def runQM(self, molfile):
        if self.qmprogram == "mopac"  or  self.qmprogram == "both":
            """
             * name and directory are the name and directory for the input (and output) file
             * input is assumed to be preexisting and have the .mop suffix
             * returns an integer indicating success or failure of the MOPAC calculation: 1 for success, 0 for failure
             * this function is based on the Gaussian analogue
            """
            jobMOPAC = jobs.MOPACJob(molfile) 
            return jobMOPAC.run()
         
        else :
            logging.critical("Unsupported quantum chemistry program")
         
        return -1 
    
    def createQMInput(self, molfile, attemptNumber, multiplicity):
        #3. create the Gaussian or MOPAC input file
        if self.qmprogram == "mopac"  or  self.qmprogram == "both":
            #write a file with the input keywords
            writer = inputwriters.MOPACPM3InputWriter(molfile, attemptNumber, multiplicity)
            inputFile = writer.write()
            
         

    def generate3DCoords(self, molecule, name):
        #steps 1 and 2: create 2D and 3D mol files
        creator = ThreeDMolFileCreator(name, self.qmfolder, molecule)
        molfile = creator.create()
        return molfile
     
    
