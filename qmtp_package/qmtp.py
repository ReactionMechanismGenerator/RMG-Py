'''
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
'''

import os
import platform
from subprocess import Popen
from rmgpy.molecule import Molecule
import qmverifier as verif
import rmg_qm_parsers as pars
import qmtpjobs as job
import qminputwriters as writers
import symmetry as symm
import logging
import openbabel
from rdkit import Chem
from rdkit.Chem import AllChem

class molFile:
    '''
    molFile class to store and perform operations on mol files
    may eventually want to have a super class threeDGeom and say "public class molFile extends threeDGeom {"
    '''
    
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
        
    def create2D(self):
        '''
        #1. create a 2D file
        #use the absolute path for directory, so we can easily reference from other directories in command-line paths
        #can't use RMG_workingDirectory, since this basically holds the RMG environment variable, not the workingDirectory
        String directory = "2Dmolfiles/"
        '''
        obmol = self.molecule.toOBMol()
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat('mol')
        output = obConversion.WriteString(obmol)
        with open(os.path.join(self.dir2D, self.mol2D),'w') as out:
            out.write(output)
        
    def embed3D(self, numConfAttempts):
        '''
        #embed a molecule in 3D, using RDKit
        '''

        m = Chem.MolFromMolFile(os.path.join(self.dir2D, self.mol2D), removeHs=False)

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
        self.create2D()

        #2. convert from 2D to 3D using RDKit if the 2D molfile is for a molecule with 2 or more atoms
        atoms = len(self.molecule.atoms)
        distGeomAttempts=1
        if atoms > 3:#this check prevents the number of attempts from being negative
            distGeomAttempts = 5*(atoms-3) #number of conformer attempts is just a linear scaling with molecule size, due to time considerations in practice, it is probably more like 3^(n-3) or something like that
            
        molfile = self.embed3D(distGeomAttempts)
        
        return molfile

class QMTP:
    '''
    Quantum mechanics thermo property estimator analog of GATP
    '''
    
    #static fields
    qmfolder= "QMfiles"
    
    
    '''
         * the qmprogram can be 
     * "mopac", 
     * "gaussian03",
     *  "both" (MOPAC and Gaussian), 
     *  "mm4",
     *  "mm4hr"
     '''

    
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
        '''
            #determine the QM filename (element 0) and augmented InChI (element 1) for a ChemGraph
            #QM filename is InChIKey appended with mult3, mult4, mult5, or mult6 for multiplicities of 3 or higher
            #augmented InChI is InChI appended with /mult3, /mult4, /mult5, or /mult6 for multiplicities of 3 or higher
        ''' 
        #inchikey_mod = molecule.getModifiedInChIKeyAnew()
        inchikey_mod = molecule.toAugmentedInChIKey()#need to generate InChI and key anew because ChemGraph may have changed (in particular, adding/removing hydrogens in HBI process)
        inchi_mod = molecule.toAugmentedInChI()#
        return inchikey_mod, inchi_mod
 
    def parseOutput(self, name, molecule):
        '''
        wrapper method for parser types
        '''

        if self.qmMethod == "pm3" :
                        
            if  self.qmprogram == "mopac" or self.qmprogram == "both" :
                parser = pars.MOPACPM3Parser(name, QMTP.qmfolder, molecule, self)
                result = parser.parse()
                result.comment = result.comment +'MOPAC PM3 calculation'
                logging.info("Thermo for " + name + ": "+ result.__repr__())#print result, at least for debugging purposes
                return result
            
        else:
                logging.critical("Unexpected situation in QMTP thermo estimation"   )    

        return result

    
    def generateQMThermoData(self, molecule):
         '''
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
         '''
         name, InChIaug = self.generateIdentifiers(molecule)#determine the filename (InChIKey) and InChI with appended info for triplets, etc.
         
         #check for existing hold file before starting calculations (to avoid the possibility of interference with other jobs using the same QMfiles folder)
         while os.path.exists(os.path.join(QMTP.qmfolder,name+'.hold')):#can lead to race condition!
             logging.info("Existence of hold file for "+name+" suggests that another RMG process is currently running calculations on this molecule; waiting for other RMG process to finish; will check again in 60 seconds...")
             import time
             time.sleep(60) # delays for 60 seconds
        
         #verify whether a succesful QM results exists for this particular species:
         verifier = verif.QMVerifier(name, InChIaug, QMTP.qmfolder, self)
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
                    self.createQMInput(name, molecule, molfile, attemptNumber, InChIaug, multiplicity)
                    success = self.runQM(molfile)
                    if success:
                        logging.info('Attempt {0} on species {1} succeeded.'.format(attemptNumber, InChIaug))
                        '''
                        TODO Rotor Scan not yet implemented here.
                        '''
                    else:
                        if attemptNumber == maxAttemptNumber:
                            logging.info('Last attempt on species {0} failed.'.format(InChIaug))
                result = self.parseOutput(name, molecule)
                return result
    
    def runQM(self, molfile):
        if self.qmprogram == "mopac"  or  self.qmprogram == "both":
            '''
             * name and directory are the name and directory for the input (and output) file
             * input is assumed to be preexisting and have the .mop suffix
             * returns an integer indicating success or failure of the MOPAC calculation: 1 for success, 0 for failure
             * this function is based on the Gaussian analogue
            '''
            jobMOPAC = job.MOPACJob(molfile, QMTP.qmfolder) 
            return jobMOPAC.run()
         
        else :
            logging.critical("Unsupported quantum chemistry program")
         
        return -1 
    
    def createQMInput(self, name, molecule, molfile, attemptNumber, inChIaug, multiplicity):
        #3. create the Gaussian or MOPAC input file
        if self.qmprogram == "mopac"  or  self.qmprogram == "both":
            #write a file with the input keywords
            writer = writers.MOPACPM3InputWriter(name, QMTP.qmfolder, molfile, attemptNumber, multiplicity)
            inputFile = writer.write()
            
         

    def generate3DCoords(self, molecule, name):
        #steps 1 and 2: create 2D and 3D mol files
        creator = ThreeDMolFileCreator(name, self.qmfolder, molecule)
        molfile = creator.create()
        return molfile
     
    
     