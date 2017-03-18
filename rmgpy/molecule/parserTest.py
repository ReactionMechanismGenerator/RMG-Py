import unittest

from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.parser import *
from rmgpy.molecule.atomtype import atomTypes
from external.wip import work_in_progress


class ParserTest(unittest.TestCase):

    def setUp(self):

        self.methane = Molecule().fromAdjacencyList("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
""")
        self.methylamine = Molecule().fromAdjacencyList("""
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 C u0 p0 c0 {1,S} {5,S} {6,S} {7,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
""")

    def test_fromAugmentedInChI(self):
        aug_inchi = 'InChI=1S/CH4/h1H4'
        mol = fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(not mol.InChI == '')
        self.assertTrue(mol.isIsomorphic(self.methane))

        aug_inchi = 'InChI=1/CH4/h1H4'
        mol = fromAugmentedInChI(Molecule(), aug_inchi)
        self.assertTrue(not mol.InChI == '')
        self.assertTrue(mol.isIsomorphic(self.methane))

    def test_fromSMILES(self):
        smiles = 'C'
        mol = fromSMILES(Molecule(), smiles)
        self.assertTrue(mol.isIsomorphic(self.methane))

        #Test that atomtypes that rely on lone pairs for identity are typed correctly
        smiles = 'CN'
        mol = fromSMILES(Molecule(), smiles)
        self.assertEquals(mol.atoms[1].atomType, atomTypes['N3s'] )

    def test_fromInChI(self):
        inchi = 'InChI=1S/CH4/h1H4'
        mol = fromInChI(Molecule(), inchi)
        self.assertTrue(mol.isIsomorphic(self.methane))
        #Test that atomtypes that rely on lone pairs for identity are typed correctly
        inchi = "InChI=1S/CH5N/c1-2/h2H2,1H3"
        mol = fromInChI(Molecule(), inchi)
        self.assertEquals(mol.atoms[1].atomType, atomTypes['N3s'] )

    #current implementation of SMARTS is broken
    def test_fromSMARTS(self):
        smarts = '[CH4]'
        mol = fromSMARTS(Molecule(), smarts)
        self.assertTrue(mol.isIsomorphic(self.methane))

    def test_toRDKitMol(self):
        """
        Test that toRDKitMol returns correct indices and atom mappings.
        """
        bondOrderDict = {'SINGLE':1,'DOUBLE':2,'TRIPLE':3,'AROMATIC':1.5}
        mol = fromSMILES(Molecule(), 'C1CCC=C1C=O')
        rdkitmol, rdAtomIndices = mol.toRDKitMol(removeHs=False, returnMapping=True, sanitize=True)
        for atom in mol.atoms:
            # Check that all atoms are found in mapping
            self.assertTrue(atom in rdAtomIndices)
            # Check that all bonds are in rdkitmol with correct mapping and order
            for connectedAtom, bond in atom.bonds.iteritems():
                bondType = str(rdkitmol.GetBondBetweenAtoms(rdAtomIndices[atom],rdAtomIndices[connectedAtom]).GetBondType())
                rdkitBondOrder = bondOrderDict[bondType]
                self.assertEqual(bond.order, rdkitBondOrder)
        
        # Test for removeHs = True        
        rdkitmol2, rdAtomIndices2 = mol.toRDKitMol(removeHs=True, returnMapping=True, sanitize=True)
        
        for atom in mol.atoms:
            # Check that all non-hydrogen atoms are found in mapping
            if atom.symbol != 'H':
                self.assertTrue(atom in rdAtomIndices)
                # Check that all bonds connected to non-hydrogen have the correct mapping and order
                for connectedAtom, bond in atom.bonds.iteritems():
                    if connectedAtom.symbol != 'H':
                        bondType = str(rdkitmol.GetBondBetweenAtoms(rdAtomIndices[atom],rdAtomIndices[connectedAtom]).GetBondType())
                        rdkitBondOrder = bondOrderDict[bondType]
                        self.assertEqual(bond.order, rdkitBondOrder)   

class ResetLonePairsTest(unittest.TestCase):

    def test_Methane(self):
        smi = 'C'
        mol = Molecule().fromSMILES(smi)
        p_indices = []

        reset_lone_pairs(mol, p_indices)

        for at in mol.atoms:
            self.assertEquals(at.lonePairs, 0)

    def test_SingletMethylene(self):
        adjlist = """
multiplicity 1
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""
        mol = Molecule().fromAdjacencyList(adjlist)
        p_indices = [1]

        reset_lone_pairs(mol, p_indices)

        for at in mol.atoms:
            if at.symbol == 'C':
                self.assertEquals(at.lonePairs, 1)
            else:
                self.assertEquals(at.lonePairs, 0)