#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

# global imports

import cython

# local imports
try:
    import openbabel
except:
    pass
from rdkit import Chem

from .molecule import Atom
from rmgpy.molecule.converter import toOBMol, toRDKitMol

import rmgpy.molecule.inchi as inchiutil

# global variables:

#: This dictionary is used to shortcut lookups of a molecule's SMILES string from its chemical formula.
_known_smiles_molecules = {
                 'N2': 'N#N',
                 'CH4': 'C',
                 'CH2O': 'C=O',
                 'H2O': 'O',
                 'C2H6': 'CC',
                 'H2': '[H][H]',
                 'H2O2': 'OO',
                 'C3H8': 'CCC',
                 'Ar': '[Ar]',
                 'He': '[He]',
                 'CH4O': 'CO',
                 'CO2': 'O=C=O',
                 'CO': '[C-]#[O+]',
                 'O2': 'O=O',
                 'C': '[C]',  # for this to be in the "molecule" list it must be singlet with 2 lone pairs
                 'H2S': 'S',
                 'N2O': 'N#[N+][O-]',
                 'NH3': 'N',
                 'O3': '[O-][O+]=O',
                 'Cl2': '[Cl][Cl]',
                 'ClH': 'Cl',
                 'I2': '[I][I]',
                 'HI': 'I',
             }

_known_smiles_radicals = {
                 'CH3': '[CH3]',
                 'HO': '[OH]',
                 'C2H5': 'C[CH2]',
                 'O': '[O]',
                 'S': '[S]',
                 'N': '[N]',
                 'HO2': '[O]O',
                 'CH': '[CH]',
                 'CH2': '[CH2]',
                 'H': '[H]',
                 'C': '[C]',  # this, in the radical list, could be triplet or quintet.
                 'O2': '[O][O]',
                 'S2': '[S][S]',
                 'OS': '[S][O]',
                 'HS': '[SH]',
                 'H2N': '[NH2]',
                 'HN': '[NH]',
                 'NO': '[N]=O',
                 'NO2': 'N(=O)[O]',
                 'Cl': '[Cl]',
                 'I': '[I]',
             }

def toInChI(mol):
    """
    Convert a molecular structure to an InChI string. Uses
    `RDKit <http://rdkit.org/>`_ to perform the conversion.
    Perceives aromaticity.
    
    or
    
    Convert a molecular structure to an InChI string. Uses
    `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
    """
    try:
        if not Chem.inchi.INCHI_AVAILABLE:
            return "RDKitInstalledWithoutInChI"
        rdkitmol = toRDKitMol(mol)
        return Chem.inchi.MolToInchi(rdkitmol, options='-SNon')
    except:
        pass

    obmol = toOBMol(mol)
    obConversion = openbabel.OBConversion()
    obConversion.SetOutFormat('inchi')
    obConversion.SetOptions('w', openbabel.OBConversion.OUTOPTIONS)
    return obConversion.WriteString(obmol).strip()


def toAugmentedInChI(mol):
    """
    This function generates the augmented InChI canonical identifier, and that allows for the differentiation
    between structures with spin states and multiple unpaired electrons.

    Two additional layers are added to the InChI:
    - unpaired electrons layer: the position of the unpaired electrons in the molecule

    """

    cython.declare(
                inchi=str,
                ulayer=str,
                aug_inchi=str,
               )
    inchi = toInChI(mol)

    ulayer, player = inchiutil.create_augmented_layers(mol)

    aug_inchi = inchiutil.compose_aug_inchi(inchi, ulayer, player)

    return aug_inchi

def toInChIKey(mol):
    """
    Convert a molecular structure to an InChI Key string. Uses
    `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
    
    or 
    
    Convert a molecular structure to an InChI Key string. Uses
    `RDKit <http://rdkit.org/>`_ to perform the conversion.
    
    Removes check-sum dash (-) and character so that only 
    the 14 + 9 characters remain.
    """
    try:
        if not Chem.inchi.INCHI_AVAILABLE:
            return "RDKitInstalledWithoutInChI"
        inchi = toInChI(mol)
        return Chem.inchi.InchiToInchiKey(inchi)[:-2]
    except:
        pass
    

#        for atom in mol.vertices:
#           if atom.isNitrogen():
    obmol = toOBMol(mol)
    obConversion = openbabel.OBConversion()
    obConversion.SetOutFormat('inchi')
    obConversion.SetOptions('w', openbabel.OBConversion.OUTOPTIONS)
    obConversion.SetOptions('K', openbabel.OBConversion.OUTOPTIONS)
    return obConversion.WriteString(obmol).strip()[:-2]

def toAugmentedInChIKey(mol):
    """
    Adds additional layers to the InChIKey,
    generating the "augmented" InChIKey.
    """
    
    cython.declare(
            key=str,
            ulayer=str
        )

    key = toInChIKey(mol)

    ulayer, player = inchiutil.create_augmented_layers(mol)

    return inchiutil.compose_aug_inchi_key(key, ulayer, player)

def toSMARTS(mol):
    """
    Convert a molecular structure to an SMARTS string. Uses
    `RDKit <http://rdkit.org/>`_ to perform the conversion.
    Perceives aromaticity and removes Hydrogen atoms.
    """
    rdkitmol = toRDKitMol(mol)
    
    return Chem.MolToSmarts(rdkitmol)


def toSMILES(mol):
    """
    Convert a molecular structure to an SMILES string. 
    
    If there is a Nitrogen atom present it uses
    `OpenBabel <http://openbabel.org/>`_ to perform the conversion,
    and the SMILES may or may not be canonical.
    
    Otherwise, it uses `RDKit <http://rdkit.org/>`_ to perform the 
    conversion, so it will be canonical SMILES.
    While converting to an RDMolecule it will perceive aromaticity
    and removes Hydrogen atoms.
    """
    
    # If we're going to have to check the formula anyway,
    # we may as well shortcut a few small known molecules.
    # Dictionary lookups are O(1) so this should be fast:
    # The dictionary is defined at the top of this file.

    cython.declare(
            atom=Atom,
            # obmol=,
            # rdkitmol=,
        )

    try:
        if mol.isRadical():
            return _known_smiles_radicals[mol.getFormula()]
        else:
            return _known_smiles_molecules[mol.getFormula()]
    except KeyError:
        # It wasn't in the above list.
        pass
    for atom in mol.vertices:
        if atom.isNitrogen():
            obmol = toOBMol(mol)
            try:
                SMILEwriter = openbabel.OBConversion()
                SMILEwriter.SetOutFormat('smi')
                SMILEwriter.SetOptions("i",SMILEwriter.OUTOPTIONS) # turn off isomer and stereochemistry information (the @ signs!)
            except:
                pass
            return SMILEwriter.WriteString(obmol).strip()

    rdkitmol = toRDKitMol(mol, sanitize=False)
    if not mol.isAromatic():
        return Chem.MolToSmiles(rdkitmol, kekuleSmiles=True)
    return Chem.MolToSmiles(rdkitmol)


