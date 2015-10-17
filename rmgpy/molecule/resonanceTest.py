import unittest

from rmgpy.molecule.molecule import Molecule

class ResonanceTest(unittest.TestCase):

    def test_C9H9_aro(self):
        """CyclopropylBenzene-radical, aromatic bonds"""
        mol = Molecule(SMILES="[CH]1CC1c1ccccc1")
        mol.generateResonanceIsomers()
    
    def test_C9H9_kek(self):
        """CyclopropylBenzene-radical, kekulized bonds"""
        mol = Molecule(SMILES="[CH]1CC1C1C=CC=CC=1")
        mol.generateResonanceIsomers()

    def test_Benzene_aro(self):
        mol = Molecule(SMILES="c1ccccc1")
        mol.generateResonanceIsomers()
    
    def test_Benzene_kek(self):
        mol = Molecule(SMILES="C1C=CC=CC=1")
        mol.generateResonanceIsomers()

    def test_C9H11_aro(self):
        """PropylBenzene-radical"""
        mol = Molecule(SMILES="[CH2]CCc1ccccc1")
        mol.generateResonanceIsomers()

    def test_C10H11_aro(self):
        """CyclobutylBenzene-radical"""
        mol = Molecule(SMILES="[CH]1CCC1c1ccccc1")
        mol.generateResonanceIsomers()

    def test_C9H10_aro(self):
        """CyclopropylBenzene, aromatic bonds"""
        mol = Molecule(SMILES="C1CC1c1ccccc1")
        mol.generateResonanceIsomers()

    def test_C10H12_aro(self):
        """CyclopropylMethylBenzene"""
        mol = Molecule(SMILES="C1CC1c1c(C)cccc1")
        mol.generateResonanceIsomers()

    def test_C9H10_aro(self):
        """CyclopropylBenzene, generate aro resonance isomers"""
        mol = Molecule(SMILES="C1CC1c1ccccc1")
        mol.getAromaticResonanceIsomers()