import os
import logging
import math

from rdkit import Chem
from rdkit.Chem import AllChem

import rmgpy.quantity
from rmgpy.thermo import ThermoData
from rmgpy.statmech import RigidRotor, HarmonicOscillator, Translation, StatesModel
import symmetry
import qmdata

class Geometry:
    file_store_path = 'QMfiles/'
    if not os.path.exists(file_store_path):
        logging.info("Creating directory %s for mol files."%os.path.abspath(file_store_path))
        os.makedirs(file_store_path)

    def __init__(self, uniqueID, molecule, multiplicity, uniqueIDlong=None):
        #: A short unique ID such as an augmented InChI Key.
        self.uniqueID = uniqueID
        self.molecule = molecule
        #: The multiplicity, eg. the number of free radicals plus one.
        self.multiplicity = multiplicity
        if uniqueIDlong is None:
            self.uniqueIDlong = uniqueID
        else:
            #: Long, truly unique, ID, such as the augmented InChI.
            self.uniqueIDlong = uniqueIDlong

    def getCrudeMolFilePath(self):
        # os.join, uniqueID, suffix
        return self.file_store_path + self.uniqueID + '.crude.mol'

    def getRefinedMolFilePath(self):
        "Returns the path the the refined mol file"
        return self.file_store_path + self.uniqueID + '.refined.mol'

    def generateRDKitGeometries(self, boundsMatrix=None):
        """
        Use RDKit to guess geometry.

        Save mol files of both crude and refined.
        Saves coordinates on atoms.
        """
        rdmol, rdAtIdx = self.rd_build()
        
        atoms = len(self.molecule.atoms)
        distGeomAttempts=1
        if atoms > 3:#this check prevents the number of attempts from being negative
            distGeomAttempts = 5*(atoms-3) #number of conformer attempts is just a linear scaling with molecule size, due to time considerations in practice, it is probably more like 3^(n-3) or something like that
        
        rdmol, minEid = self.rd_embed(rdmol, distGeomAttempts, boundsMatrix)
        self.save_coordinates(rdmol, minEid, rdAtIdx)
        
    def rd_build(self):
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

        return rdmol, rdAtomIdx

    def rd_embed(self, rdmol, numConfAttempts, boundsMatrix=None):
        """
        Embed the RDKit molecule and create the crude molecule file.
        """
        AllChem.EmbedMultipleConfs(rdmol, numConfAttempts,randomSeed=1)
        
        energy=0.0
        minEid=0;
        lowestE=9.999999e99;#start with a very high number, which would never be reached

        crude = Chem.Mol(rdmol.ToBinary())
        
        for i in range(rdmol.GetNumConformers()):
            AllChem.UFFOptimizeMolecule(rdmol,confId=i)
            energy=AllChem.UFFGetMoleculeForceField(rdmol,confId=i).CalcEnergy()
            if energy < lowestE:
                minEid = i
                lowestE = energy 
        
        with open(self.getCrudeMolFilePath(), 'w') as out3Dcrude:
            out3Dcrude.write(Chem.MolToMolBlock(crude,confId=minEid))
        
        with open(self.getRefinedMolFilePath(), 'w') as out3D:
            out3D.write(Chem.MolToMolBlock(rdmol,confId=minEid))

        return rdmol, minEid

    def save_coordinates(self, rdmol, minEid, rdAtIdx):
        # Save xyz coordinates on each atom in molecule ****
        for atom in self.molecule.atoms:
            point = rdmol.GetConformer(minEid).GetAtomPosition(atom.sortingLabel)
            atom.coords = [point.x, point.y, point.z]
            


class QMMolecule:
    def __init__(self, molecule):
        self.molecule = molecule

    def createGeometry(self):
        """
        Creates self.geometry with RDKit geometries
        """
        uniqueID = self.molecule.toAugmentedInChIKey()
        uniqueIDlong = self.molecule.toAugmentedInChI()
        multiplicity = sum([i.radicalElectrons for i in self.molecule.atoms]) + 1
        self.geometry = Geometry(uniqueID, self.molecule, multiplicity, uniqueIDlong=uniqueIDlong)
        
        self.geometry.generateRDKitGeometries()
        
        return self.geometry
    
    def generateQMData(self):
        """
        Calculate the QM data somehow and return a CCLibData object.
        """
        raise NotImplementedError("This should be defined in a subclass that inherits from QMMolecule")
        return qmdata.QMData()
        
    
    def generateThermoData(self):
        """
        Generate Thermo Data via a QM calc
        """
        
        # First generate the QM data
        result = self.generateQMData()
        
        
        thermo = TDPropertiesCalculator(result, self.getInChiKeyAug())
        
        
        return thermo.calculate()

    def getInChiKeyAug(self):
        """
        Returns the augmented InChI from self.molecule 
        """        
        return self.molecule.toAugmentedInChIKey()


class TDPropertiesCalculator:
    """
    A thermodynamic properties calculator.
    """
    def __init__(self, qmdata, molfile, pointGroup = None):
        self.qmdata = qmdata

        self.molfile = molfile

        self.determinePointGroup()

        trans = Translation(mass=(qmdata.molecularMass,"amu"))
        rot = RigidRotor(linear=self.pointGroup.linear, inertia=self.qmdata.rotationalConstants, symmetry=self.pointGroup.symmetryNumber)
        vib = HarmonicOscillator(frequencies=qmdata.frequencies)
        self.statesmodel = StatesModel(modes=[trans, rot, vib], spinMultiplicity=self.qmdata.groundStateDegeneracy)


    def determinePointGroup(self):
        #determine point group using the SYMMETRY Program

        pgc = symmetry.PointGroupCalculator(self.molfile, self.qmdata);
        self.pointGroup = pgc.calculate();

    def calculateChiralityCorrection(self):
        if self.pointGroup.chiral:
            return rmgpy.quantity.constants.R * math.log(2);
        else:
            return 0.

    def calculate(self):
        #we will use number of atoms from above (alternatively, we could use the chemGraph); this is needed to test whether the species is monoatomic
        Hartree_to_kcal = 627.5095#conversion from Hartree to kcal/mol taken from Gaussian thermo white paper
        Hf298 = self.qmdata.energy * Hartree_to_kcal;

        S298 = self.statesmodel.getEntropy(298.0)

        Tdata = [300.0,400.0,500.0,600.0,800.0,1000.0,1500.0]
        Cp = []
        for T in Tdata:
            Cp.append(self.statesmodel.getHeatCapacity(T))

        S298 = S298 + self.calculateChiralityCorrection()

        Tmin = 300.0

        Tmax = 2000.0

        comment = "PM3 or MM4 calculation"

        return ThermoData(Tdata, Cp, Hf298, S298, Tmin, Tmax, comment)