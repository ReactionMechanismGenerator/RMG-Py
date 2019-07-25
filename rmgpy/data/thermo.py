#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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

"""

"""

import os.path
import re
import math
import logging
import numpy
import time
import itertools
from copy import deepcopy

from base import Database, Entry, makeLogicNode, DatabaseError

import rmgpy.constants as constants
from rmgpy.thermo import NASAPolynomial, NASA, ThermoData, Wilhoit
from rmgpy.molecule import Molecule, Bond, Group
import rmgpy.molecule
from rmgpy.species import Species
import rmgpy.quantity
from rmgpy.ml.estimator import MLEstimator

#: This dictionary is used to add multiplicity to species label
_multiplicity_labels = {1:'S',2:'D',3:'T',4:'Q',5:'V',}

################################################################################

def saveEntry(f, entry):
    """
    Write a Pythonic string representation of the given `entry` in the thermo
    database to the file object `f`.
    """
    
    f.write('entry(\n')
    f.write('    index = {0:d},\n'.format(entry.index))
    f.write('    label = "{0}",\n'.format(entry.label))

    if isinstance(entry.item, Molecule):
        f.write('    molecule = \n')
        f.write('"""\n')
        f.write(entry.item.toAdjacencyList(removeH=False))
        f.write('""",\n')
    elif isinstance(entry.item, Group):
        f.write('    group = \n')
        f.write('"""\n')
        f.write(entry.item.toAdjacencyList())
        f.write('""",\n')
    else:
        f.write('    group = "{0}",\n'.format(entry.item))

    if isinstance(entry.data, ThermoData):
        f.write('    thermo = ThermoData(\n')
        f.write('        Tdata = {0!r},\n'.format(entry.data.Tdata))
        f.write('        Cpdata = {0!r},\n'.format(entry.data.Cpdata))
        f.write('        H298 = {0!r},\n'.format(entry.data.H298))
        f.write('        S298 = {0!r},\n'.format(entry.data.S298))
        if entry.data.Tmin is not None: f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None: f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        f.write('    ),\n')
    elif isinstance(entry.data, Wilhoit):
        f.write('    thermo = Wilhoit(\n')
        f.write('        cp0 = {0!r},\n'.format(entry.data.cp0))
        f.write('        cpInf = {0!r},\n'.format(entry.data.cpInf))
        f.write('        a0 = {0:g},\n'.format(entry.data.a0))
        f.write('        a1 = {0:g},\n'.format(entry.data.a1))
        f.write('        a2 = {0:g},\n'.format(entry.data.a2))
        f.write('        a3 = {0:g},\n'.format(entry.data.a3))
        f.write('        B = {0!r},\n'.format(entry.data.B))
        f.write('        H0 = {0!r},\n'.format(entry.data.H0))
        f.write('        S0 = {0!r},\n'.format(entry.data.S0))
        if entry.data.Tmin is not None: f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None: f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        f.write('    ),\n')
    elif isinstance(entry.data, NASA):
        f.write('    thermo = NASA(\n')
        f.write('        polynomials = [\n')
        for poly in entry.data.polynomials:
            f.write('            {0!r},\n'.format(poly))
        f.write('        ],\n')
        if entry.data.Tmin is not None: f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None: f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        if entry.data.E0 is not None: f.write('        E0 = {0!r},\n'.format(entry.data.E0))
        if entry.data.Cp0 is not None: f.write('        Cp0 = {0!r},\n'.format(entry.data.Cp0))
        if entry.data.CpInf is not None: f.write('        CpInf = {0!r},\n'.format(entry.data.CpInf))
        f.write('    ),\n')
    else:
        f.write('    thermo = {0!r},\n'.format(entry.data))

    if entry.reference is not None: f.write('    reference = {0!r},\n'.format(entry.reference))
    if entry.referenceType != "": f.write('    referenceType = "{0}",\n'.format(entry.referenceType))
    f.write('    shortDesc = u"""')
    try:
        f.write(entry.shortDesc.encode('utf-8'))
    except (UnicodeEncodeError, UnicodeDecodeError):
        f.write(entry.shortDesc.strip().encode('ascii', 'replace'))
    f.write('""",\n')
    f.write('    longDesc = \n')
    f.write('u"""\n')
    try:
        f.write(entry.longDesc.strip().encode('utf-8') + "\n")    
    except (UnicodeEncodeError, UnicodeDecodeError):
        f.write(entry.longDesc.strip().encode('ascii', 'replace') + "\n")
    f.write('""",\n')
    if entry.rank:
        f.write("    rank = {0},\n".format(entry.rank))

    f.write(')\n\n')

def generateOldLibraryEntry(data):
    """
    Return a list of values used to save entries to the old-style RMG
    thermo database based on the thermodynamics object `data`.
    """
    if isinstance(data, ThermoData):
        return '{0:9g} {1:9g} {2:9g} {3:9g} {4:9g} {5:9g} {6:9g} {7:9g} {8:9g} {9:9g} {10:9g} {11:9g}'.format(
            data.H298.value_si/4184.,
            data.S298.value_si/4.184,
            data.Cpdata.value_si[0]/4.184,
            data.Cpdata.value_si[1]/4.184,
            data.Cpdata.value_si[2]/4.184,
            data.Cpdata.value_si[3]/4.184,
            data.Cpdata.value_si[4]/4.184,
            data.Cpdata.value_si[5]/4.184,
            data.Cpdata.value_si[6]/4.184,
            data.H298.uncertainty/4184.,
            data.S298.uncertainty/4.184,
            max(data.Cpdata.uncertainty)/4.184,
        )
    elif isinstance(data, basestring):
        return data
    else:
        return '{0:9g} {1:9g} {2:9g} {3:9g} {4:9g} {5:9g} {6:9g} {7:9g} {8:9g} {9:9g} {10:9g} {11:9g}'.format(
            data.getEnthalpy(298)/4184.,
            data.getEntropy(298)/4.184,
            data.getHeatCapacity(300)/4.184,
            data.getHeatCapacity(400)/4.184,
            data.getHeatCapacity(500)/4.184,
            data.getHeatCapacity(600)/4.184,
            data.getHeatCapacity(800)/4.184,
            data.getHeatCapacity(1000)/4.184,
            data.getHeatCapacity(1500)/4.184,
            0,
            0,
            0,
        )

def processOldLibraryEntry(data):
    """
    Process a list of parameters `data` as read from an old-style RMG
    thermo database, returning the corresponding thermodynamics object.
    """
    return ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],"K"),
        Cpdata = ([float(d) for d in data[2:9]],"cal/(mol*K)","+|-",float(data[11])),
        H298 = (float(data[0]),"kcal/mol","+|-",float(data[9])),
        S298 = (float(data[1]),"cal/(mol*K)","+|-",float(data[10])),
    )


def addThermoData(thermoData1, thermoData2, groupAdditivity=False, verbose=False):
        """
        Add the thermodynamic data `thermoData2` to the data `thermoData1`,
        and return `thermoData1`.
        
        If `groupAdditivity` is True, append comments related to group additivity estimation
        If `verbose` is False, omit the comments from a "zero entry", whose H298, S298, and Cp are all 0.
        If `verbose` is True, or thermoData2 is not a zero entry, add thermoData2.comment to thermoData1.comment.
        """
        if len(thermoData1.Tdata.value_si) != len(thermoData2.Tdata.value_si) or any([T1 != T2 for T1, T2 in zip(thermoData1.Tdata.value_si, thermoData2.Tdata.value_si)]):
            raise Exception('Cannot add these ThermoData objects due to their having different temperature points.')
        
        for i in range(thermoData1.Tdata.value_si.shape[0]):
            thermoData1.Cpdata.value_si[i] += thermoData2.Cpdata.value_si[i]
        thermoData1.H298.value_si += thermoData2.H298.value_si
        thermoData1.S298.value_si += thermoData2.S298.value_si

        testZero = sum(abs(value) for value in [thermoData2.H298.value_si, thermoData2.S298.value_si]+thermoData2.Cpdata.value_si.tolist())
        # Used to check if all of the entries in thermoData2 are zero

        if groupAdditivity:
            if verbose or testZero !=0:
            # If verbose==True or testZero!=0, add thermoData2.comment to thermoData1.comment.
                if thermoData1.comment:
                    thermoData1.comment += ' + {0}'.format(thermoData2.comment)
                else:
                    thermoData1.comment = 'Thermo group additivity estimation: ' + thermoData2.comment
            
        return thermoData1
    
def removeThermoData(thermoData1, thermoData2, groupAdditivity=False, verbose=False):
    """
    Remove the thermodynamic data `thermoData2` from the data `thermoData1`,
    and return `thermoData1`.
    If `verbose` is True, append ' - thermoData2.comment' to the thermoData1.comment.
    If `verbose` is False, remove the thermoData2.comment from the thermoData1.comment.
    """
    if len(thermoData1.Tdata.value_si) != len(thermoData2.Tdata.value_si) or any([T1 != T2 for T1, T2 in zip(thermoData1.Tdata.value_si, thermoData2.Tdata.value_si)]):
        raise Exception('Cannot take the difference between these ThermoData objects due to their having different temperature points.')

    for i in range(thermoData1.Tdata.value_si.shape[0]):
        thermoData1.Cpdata.value_si[i] -= thermoData2.Cpdata.value_si[i]
    thermoData1.H298.value_si -= thermoData2.H298.value_si
    thermoData1.S298.value_si -= thermoData2.S298.value_si

    if groupAdditivity:
        if verbose:
            thermoData1.comment += ' - {0}'.format(thermoData2.comment)
        else:
            thermoData1.comment = re.sub(re.escape(' + '+thermoData2.comment),'',thermoData1.comment, 1)
    return thermoData1

def averageThermoData(thermoDataList=None):
    """
    Average a list of thermoData values together.
    Sets uncertainty values to be the approximately the 95% confidence interval, equivalent to
    2 standard deviations calculated using the sample standard variance:
    
    Uncertainty = 2s
    s = sqrt( sum(abs(x - x.mean())^2) / N - 1) where N is the number of values averaged
    
    Note that uncertainties are only computed when number of values is greater than 1.
    """
    if thermoDataList is None:
        thermoDataList = []

    numValues = len(thermoDataList)
        
    if numValues == 0:
        raise Exception('No thermo data values were inputted to be averaged.')
    else:
        logging.debug('Averaging thermo data over {0} value(s).'.format(numValues))
        
        if numValues == 1:
            return deepcopy(thermoDataList[0])
        
        else:
            averagedThermoData = deepcopy(thermoDataList[0])
            for thermoData in thermoDataList[1:]:
                averagedThermoData = addThermoData(averagedThermoData, thermoData)


            for i in range(averagedThermoData.Tdata.value_si.shape[0]):
                averagedThermoData.Cpdata.value_si[i] /= numValues

                cpData = [thermoData.Cpdata.value_si[i] for thermoData in thermoDataList]
                averagedThermoData.Cpdata.uncertainty[i] = 2*numpy.std(cpData, ddof=1)

            HData = [thermoData.H298.value_si for thermoData in thermoDataList]
            averagedThermoData.H298.value_si /= numValues
            averagedThermoData.H298.uncertainty_si = 2*numpy.std(HData, ddof=1)

            SData = [thermoData.S298.value_si for thermoData in thermoDataList]
            averagedThermoData.S298.value_si /= numValues
            averagedThermoData.S298.uncertainty_si = 2*numpy.std(SData, ddof=1)
            return averagedThermoData

def commonAtoms(cycle1, cycle2):
    """
    INPUT: two cycles with type: list of atoms
    OUTPUT: a set of common atoms
    """
    set1 = set(cycle1)
    set2 = set(cycle2)
    return set1.intersection(set2)

def combineCycles(cycle1, cycle2):
    """
    INPUT: two cycles with type: list of atoms
    OUTPUT: a combined cycle with type: list of atoms
    """
    set1 = set(cycle1)
    set2 = set(cycle2)
    return list(set1.union(set2))

def isAromaticRing(submol):
    """
    This method takes a monoring submol (Molecule initialized with a list of atoms containing just 
    the ring), and check if it is a aromatic ring.
    """
    ring_size = len(submol.atoms)
    if ring_size not in [5, 6]:
        return False
    for ringAtom in submol.atoms:
        for bondedAtom, bond in ringAtom.edges.iteritems():
            if bondedAtom in submol.atoms:
                if not bond.isBenzene():
                    return False
    return True

def isBicyclic(polyring):
    """
    Given a polyring (a list of `Atom`s)
    returns True if it's a bicyclic, False otherwise
    """
    submol, _ = convertRingToSubMolecule(polyring)
    sssr = submol.getSmallestSetOfSmallestRings()

    return len(sssr) == 2

def findAromaticBondsFromSubMolecule(submol):
    """
    This method finds all the aromatic bonds within a input submolecule and 
    returns a set of unique aromatic bonds
    """

    aromaticBonds = []
    for atom in submol.atoms:
        bonds = submol.getBonds(atom)
        for atom_j in bonds:
            if atom_j in submol.atoms:
                bond = bonds[atom_j]
                if bond.isBenzene():
                    aromaticBonds.append(bond)
    return set(aromaticBonds)

def convertRingToSubMolecule(ring):
    """
    This function takes a ring structure (can either be monoring or polyring) to create a new 
    submolecule with newly deep copied atoms

    Outputted submolecules may have incomplete valence and may cause errors with some Molecule.methods(), such
    as updateAtomTypes() or update(). In the future we may consider using groups for the sub-molecules.
    """
    
    atomsMapping = {}
    for atom in ring:
        atomsMapping[atom] = atom.copy() # this copy is deep copy of origin atom with empty edges

    mol0 = Molecule(atoms=atomsMapping.values())

    for atom in ring:
        for bondedAtom, bond in atom.edges.iteritems():
            if bondedAtom in ring:
                if not mol0.hasBond(atomsMapping[atom],atomsMapping[bondedAtom]):
                    mol0.addBond(Bond(atomsMapping[atom],atomsMapping[bondedAtom],order=bond.order))

    mol0.updateMultiplicity()
    mol0.updateConnectivityValues()
    return mol0, atomsMapping

def combineTwoRingsIntoSubMolecule(ring1, ring2):
    """
    This function combines 2 rings (with common atoms) to create a new 
    submolecule with newly deep copied atoms
    """

    assert len(commonAtoms(ring1, ring2))>0, "The two input rings don't have common atoms."

    atomsMapping = {}
    for atom in ring1 + ring2:
        if atom not in atomsMapping:
            atomsMapping[atom] = atom.copy()

    mol0 = Molecule(atoms=atomsMapping.values())

    for atom in ring1:
        for bondedAtom, bond in atom.edges.iteritems():
            if bondedAtom in ring1:
                if not mol0.hasBond(atomsMapping[atom],atomsMapping[bondedAtom]):
                    mol0.addBond(Bond(atomsMapping[atom],atomsMapping[bondedAtom],order=bond.order))
    
    for atom in ring2:
        for bondedAtom, bond in atom.edges.iteritems():
            if bondedAtom in ring2:
                if not mol0.hasBond(atomsMapping[atom],atomsMapping[bondedAtom]):
                    mol0.addBond(Bond(atomsMapping[atom],atomsMapping[bondedAtom],order=bond.order))
    
    mol0.updateMultiplicity()
    mol0.updateConnectivityValues()

    return mol0, atomsMapping

def getCopyForOneRing(ring):

    _, atomsMapping = convertRingToSubMolecule(ring)

    ringCopy = [atomsMapping[atom] for atom in ring]
    
    return ringCopy


def getCopyFromTwoRingsWithCommonAtoms(ring1, ring2):

    mergedRing, atomsMapping = combineTwoRingsIntoSubMolecule(ring1, ring2)

    ring1Copy = [atomsMapping[atom] for atom in ring1]
    ring2Copy = [atomsMapping[atom] for atom in ring2]
    
    return ring1Copy, ring2Copy, mergedRing

def isRingPartialMatched(ring, matched_group):
    """
    An example of ring partial match is tricyclic ring is matched by a bicyclic group
    usually because of not enough data in polycyclic tree. The method takes a matched group 
    returned from descendTree and the ring (a list of non-hydrogen atoms in the ring)
    """
    # if matched group has less atoms than the target ring
    # it's surely a partial match
    if len(ring) > len(matched_group.atoms):
        return True
    else:
        submol_ring, _ = convertRingToSubMolecule(ring)
        sssr = submol_ring.getSmallestSetOfSmallestRings()
        sssr_grp = matched_group.getSmallestSetOfSmallestRings()
        if sorted([len(sr) for sr in sssr]) == sorted([len(sr_grp) for sr_grp in sssr_grp]):
            return False
        else:
            return True

def bicyclicDecompositionForPolyring(polyring):
    """
    Decompose a polycyclic ring into all possible bicyclic combinations: `bicyclicsMergedFromRingPair`
    and return a `ringOccurancesDict` that contains all single ring tuples as keys and the number of times
    they appear each bicyclic submolecule.  These bicyclic and single rings are used 
    later in the heuristic polycyclic thermo algorithm.
    """

    submol, _ = convertRingToSubMolecule(polyring)
    SSSR = submol.getDeterministicSmallestSetOfSmallestRings()

    ringPairWithCommonAtomsList = []
    ringOccurancesDict = {}
    
    # Initialize ringOccuranceDict
    for ring in SSSR:
        ringOccurancesDict[tuple(ring)] = 0

    ringNum = len(SSSR)
    for i in range(ringNum):
        for j in range(i+1,ringNum):
            if commonAtoms(SSSR[i], SSSR[j]):
                # Copy the SSSR's again because these ones are going to be merged into bicyclics
                # and manipulated (aromatic bonds have to be screened and changed to single if needed)
                SSSRi, SSSRj, mergedRing = getCopyFromTwoRingsWithCommonAtoms(SSSR[i], SSSR[j])
                ringPairWithCommonAtomsList.append([SSSRi, SSSRj, mergedRing])
                # Save the single ring SSSRs that appear in bicyclics using the original copy
                # because they will be manipulated (differently) in __addPolyRingCorrectionThermoDataFromHeuristic
                ringOccurancesDict[tuple(SSSR[i])] += 1
                ringOccurancesDict[tuple(SSSR[j])] += 1

    bicyclicsMergedFromRingPair = []
    # pre-process 2-ring cores
    for ringA, ringB, mergedRing in ringPairWithCommonAtomsList:
        submolA = Molecule(atoms=ringA)
        submolB = Molecule(atoms=ringB)
        isA_aromatic = isAromaticRing(submolA)
        isB_aromatic = isAromaticRing(submolB)
        # if ringA and ringB are both aromatic or not aromatic
        # don't need to do anything extra
        if (isA_aromatic and isB_aromatic):
            pass
        elif (not isA_aromatic and not isB_aromatic):
            aromaticBonds_inA = findAromaticBondsFromSubMolecule(submolA)
            for aromaticBond_inA in aromaticBonds_inA:
                aromaticBond_inA.setOrderNum(1)

            aromaticBonds_inB = findAromaticBondsFromSubMolecule(submolB)
            for aromaticBond_inB in aromaticBonds_inB:
                aromaticBond_inB.setOrderNum(1)
        elif isA_aromatic:
            aromaticBonds_inB = findAromaticBondsFromSubMolecule(submolB)
            for aromaticBond_inB in aromaticBonds_inB:
                # Make sure the aromatic bond in ringB is in ringA, and both ringB atoms are in ringA 
                # If so, preserve the B bond status, otherwise change to single bond order
                if (aromaticBond_inB.atom1 in submolA.atoms) and (aromaticBond_inB.atom2 in submolA.atoms) and (submolA.hasBond(aromaticBond_inB.atom1, aromaticBond_inB.atom2)):
                    pass
                else:
                    aromaticBond_inB.setOrderNum(1)
        else:
            aromaticBonds_inA = findAromaticBondsFromSubMolecule(submolA)
            for aromaticBond_inA in aromaticBonds_inA:
                if (aromaticBond_inA.atom1 in submolB.atoms) and (aromaticBond_inA.atom2 in submolB.atoms) and (submolB.hasBond(aromaticBond_inA.atom1, aromaticBond_inA.atom2)):
                    pass
                else:
                    aromaticBond_inA.setOrderNum(1)
        mergedRing.saturate_unfilled_valence(update = True)
        bicyclicsMergedFromRingPair.append(mergedRing)

    return bicyclicsMergedFromRingPair, ringOccurancesDict

def splitBicyclicIntoSingleRings(bicyclic_submol):
    """
    Splits a given bicyclic submolecule into two individual single 
    ring submolecules (a list of `Molecule`s ).
    """
    SSSR = bicyclic_submol.getDeterministicSmallestSetOfSmallestRings()

    return [convertRingToSubMolecule(SSSR[0])[0],
                convertRingToSubMolecule(SSSR[1])[0]]

def saturate_ring_bonds(ring_submol):
    """
    Given a ring submolelcule (`Molecule`), makes a deep copy and converts non-single bonds 
    into single bonds, returns a new saturated submolecule (`Molecule`)
    """
    atomsMapping = {}
    for atom in ring_submol.atoms:
        if atom not in atomsMapping:
            atomsMapping[atom] = atom.copy()

    mol0 = Molecule(atoms=atomsMapping.values())

    alreadySaturated = True
    for atom in ring_submol.atoms:
        for bondedAtom, bond in atom.edges.iteritems():
            if bondedAtom in ring_submol.atoms:
                if bond.order > 1.0 and not bond.isBenzene(): alreadySaturated = False
                if not mol0.hasBond(atomsMapping[atom],atomsMapping[bondedAtom]):
                    bond_order = 1.0
                    if bond.isBenzene():
                        bond_order = 1.5
                    mol0.addBond(Bond(atomsMapping[atom],atomsMapping[bondedAtom],order=bond_order))

    mol0.saturate_unfilled_valence()
    mol0.updateAtomTypes()
    mol0.updateMultiplicity()
    mol0.updateConnectivityValues()
    return mol0, alreadySaturated

################################################################################

class ThermoDepository(Database):
    """
    A class for working with the RMG thermodynamics depository.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self, index, label, molecule, thermo, reference=None, referenceType='', shortDesc='', longDesc='', rank=None):
        entry = Entry(
            index = index,
            label = label,
            item = Molecule().fromAdjacencyList(molecule),
            data = thermo,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            rank = rank,
        )
        self.entries[label] = entry
        return entry

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return saveEntry(f, entry)

################################################################################

class ThermoLibrary(Database):
    """
    A class for working with a RMG thermodynamics library.
    """

    def __init__(self, label='', name='',solvent=None, shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self,
                  index,
                  label,
                  molecule,
                  thermo,
                  reference=None,
                  referenceType='',
                  shortDesc='',
                  longDesc='',
                  rank=None,
                  ):
        
        molecule = Molecule().fromAdjacencyList(molecule)
        
        # Internal checks for adding entry to the thermo library
        if label in self.entries.keys():
            raise DatabaseError('Found a duplicate molecule with label {0} in the thermo library {1}.  Please correct your library.'.format(label, self.name))
        
        for entry in self.entries.values():
            if molecule.isIsomorphic(entry.item):
                if molecule.multiplicity == entry.item.multiplicity:
                    raise DatabaseError('Adjacency list and multiplicity of {0} matches that of existing molecule {1} in thermo library {2}.  Please correct your library.'.format(label, entry.label, self.name))
        
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = molecule,
            data = thermo,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            rank = rank,
        )

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return saveEntry(f, entry)

    def generateOldLibraryEntry(self, data):
        """
        Return a list of values used to save entries to the old-style RMG
        thermo database based on the thermodynamics object `data`.
        """
        return generateOldLibraryEntry(data)

    def processOldLibraryEntry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        thermo database, returning the corresponding thermodynamics object.
        """
        return processOldLibraryEntry(data)

################################################################################

class ThermoGroups(Database):
    """
    A class for working with an RMG thermodynamics group additivity database.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self,
                  index,
                  label,
                  group,
                  thermo,
                  reference=None,
                  referenceType='',
                  shortDesc='',
                  longDesc='',
                  rank=None,
                  ):

        if group[0:3].upper() == 'OR{' or group[0:4].upper() == 'AND{' or group[0:7].upper() == 'NOT OR{' or group[0:8].upper() == 'NOT AND{':
            item = makeLogicNode(group)
        else:
            item = Group().fromAdjacencyList(group)
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = item,
            data = thermo,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            rank = rank,
        )
    
    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return saveEntry(f, entry)

    def generateOldLibraryEntry(self, data):
        """
        Return a list of values used to save entries to the old-style RMG
        thermo database based on the thermodynamics object `data`.
        """
        
        return generateOldLibraryEntry(data)

    def processOldLibraryEntry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        thermo database, returning the corresponding thermodynamics object.
        """
        return processOldLibraryEntry(data)

    def copyData(self, source, destination):
        """
        This method copys the ThermoData object and all meta data
        from source to destination
        Args:
            source: The entry for which data is being copied
            destination: The entry for which data is being overwritten

        """
        destination.data = source.data
        destination.reference = source.reference
        destination._longDesc = source._longDesc
        destination._shortDesc = source._shortDesc
        destination._longDesc = source.longDesc
        destination._shortDesc = source.shortDesc
        destination.rank = source.rank
        destination.referenceType = source.referenceType    


    def removeGroup(self, groupToRemove):
        """
        Removes a group that is in a tree from the database. For thermo
        groups we also, need to re-point any unicode thermoData that may
        have pointed to the entry.

        Returns the removed group
        """

        #First call base class method
        Database.removeGroup(self, groupToRemove)

        parentR = groupToRemove.parent

        #look for other pointers that point toward entry
        for entry in self.entries.itervalues():
            if isinstance(entry.data, basestring):
                if entry.data == groupToRemove.label:
                    #if the entryToRemove.data is also a pointer, then copy
                    if isinstance(groupToRemove.data, basestring):
                        entry.data = groupToRemove.data
                    #if the parent points toward entry and the data is
                    #not a base string, we need to copy the data to the parent
                    elif entry is parentR:
                        self.copyData(groupToRemove, parentR)
                    #otherwise, point toward entryToRemove's parent
                    else:
                        entry.data = unicode(parentR.label)

        return groupToRemove
################################################################################

class ThermoDatabase(object):
    """
    A class for working with the RMG thermodynamics database.
    """

    def __init__(self):
        self.depository = {}
        self.libraries = {}
        self.groups = {}
        self.libraryOrder = []
        self.local_context = {
            'ThermoData': ThermoData,
            'Wilhoit': Wilhoit,
            'NASAPolynomial': NASAPolynomial,
            'NASA': NASA,
        }
        self.global_context = {}

        # Catalyst properties
        self.setDeltaAtomicAdsorptionEnergies()

    def __reduce__(self):
        """
        A helper function used when pickling a ThermoDatabase object.
        """
        d = {
            'depository': self.depository,
            'libraries': self.libraries,
            'groups': self.groups,
            'libraryOrder': self.libraryOrder,
        }
        return (ThermoDatabase, (), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling a ThermoDatabase object.
        """
        self.depository = d['depository']
        self.libraries = d['libraries']
        self.groups = d['groups']
        self.libraryOrder = d['libraryOrder']

    def load(self, path, libraries=None, depository=True):
        """
        Load the thermo database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        if depository:
            self.loadDepository(os.path.join(path, 'depository'))
        else:
            self.depository = {}
        self.loadLibraries(os.path.join(path, 'libraries'), libraries)
        self.loadGroups(os.path.join(path, 'groups'))

    def loadDepository(self, path):
        """
        Load the thermo database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        self.depository = {}
        self.depository['stable']  = ThermoDepository().load(os.path.join(path, 'stable.py'), self.local_context, self.global_context)
        self.depository['radical'] = ThermoDepository().load(os.path.join(path, 'radical.py'), self.local_context, self.global_context)

    def loadLibraries(self, path, libraries=None):
        """
        Load the thermo database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        
        If no libraries are given, all are loaded.
        """
        self.libraries = {}; self.libraryOrder = []
        if libraries is None:
            for (root, dirs, files) in os.walk(os.path.join(path)):
                for f in files:
                    name, ext = os.path.splitext(f)
                    if ext.lower() == '.py':
                        logging.info('Loading thermodynamics library from {0} in {1}...'.format(f, root))
                        library = ThermoLibrary()
                        library.load(os.path.join(root, f), self.local_context, self.global_context)
                        library.label = os.path.splitext(f)[0]
                        self.libraries[library.label] = library
                        self.libraryOrder.append(library.label)

        else:
            for libraryName in libraries:
                f = libraryName + '.py'
                if os.path.exists(os.path.join(path, f)):
                    logging.info('Loading thermodynamics library from {0} in {1}...'.format(f, path))
                    library = ThermoLibrary()
                    library.load(os.path.join(path, f), self.local_context, self.global_context)
                    library.label = os.path.splitext(f)[0]
                    self.libraries[library.label] = library
                    self.libraryOrder.append(library.label)
                else:
                    if libraryName == "KlippensteinH2O2":
                        logging.info("""\n** Note: The thermo library KlippensteinH2O2 was replaced and is no longer available in RMG.
                        For H2 combustion chemistry consider using the BurkeH2O2 library instead\n""")
                    raise Exception('Library {} not found in {}...Please check if your library is correctly placed'.format(libraryName, path))

    def loadGroups(self, path):
        """
        Load the thermo database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        logging.info('Loading thermodynamics group database from {0}...'.format(path))
        categories = [
            'group',
            'ring',
            'radical',
            'polycyclic',
            'other',
            'longDistanceInteraction_cyclic',
            'longDistanceInteraction_noncyclic',
            'adsorptionPt',
        ]
        self.groups = {
            category: ThermoGroups(label=category).load(os.path.join(path, category + '.py'), self.local_context, self.global_context)
            for category in categories
        }

        self.recordRingGenericNodes()
        self.recordPolycylicGenericNodes()

    def save(self, path):
        """
        Save the thermo database to the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        path = os.path.abspath(path)
        if not os.path.exists(path): os.mkdir(path)
        self.saveDepository(os.path.join(path, 'depository'))
        self.saveLibraries(os.path.join(path, 'libraries'))
        self.saveGroups(os.path.join(path, 'groups'))

    def saveDepository(self, path):
        """
        Save the thermo depository to the given `path` on disk, where `path`
        points to the top-level folder of the thermo depository.
        """
        if not os.path.exists(path): os.mkdir(path)
        for depo in self.depository.keys():
            self.depository[depo].save(os.path.join(path, depo+'.py'))

    def saveLibraries(self, path):
        """
        Save the thermo libraries to the given `path` on disk, where `path`
        points to the top-level folder of the thermo libraries.
        """
        if not os.path.exists(path): os.mkdir(path)
        for library in self.libraries.values():
            library.save(os.path.join(path, '{0}.py'.format(library.label)))

    def saveGroups(self, path):
        """
        Save the thermo groups to the given `path` on disk, where `path`
        points to the top-level folder of the thermo groups.
        """
        if not os.path.exists(path): os.mkdir(path)
        for group in self.groups.keys():
            self.groups[group].save(os.path.join(path, group+'.py'))
            
    def loadOld(self, path):
        """
        Load the old RMG thermo database from the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        # The old database does not have a depository, so create an empty one
        self.depository = {}
        self.depository['stable']  = ThermoDepository(label='stable', name='Stable Molecules')
        self.depository['radical'] = ThermoDepository(label='radical', name='Radical Molecules')
        
        for (root, dirs, files) in os.walk(os.path.join(path, 'thermo_libraries')):
            if os.path.exists(os.path.join(root, 'Dictionary.txt')) and os.path.exists(os.path.join(root, 'Library.txt')):
                library = ThermoLibrary(label=os.path.basename(root), name=os.path.basename(root))
                library.loadOld(
                    dictstr = os.path.join(root, 'Dictionary.txt'),
                    treestr = '',
                    libstr = os.path.join(root, 'Library.txt'),
                    numParameters = 12,
                    numLabels = 1,
                    pattern = False,
                )
                library.label = os.path.basename(root)
                self.libraries[library.label] = library

        self.groups = {}
        self.groups['group'] = ThermoGroups(label='group', name='Functional Group Additivity Values').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', 'Group_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', 'Group_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', 'Group_Library.txt'),
            numParameters = 12,
            numLabels = 1,
            pattern = True,
        )
        self.groups['gauche'] = ThermoGroups(label='gauche', name='Gauche Interaction Corrections').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', 'Gauche_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', 'Gauche_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', 'Gauche_Library.txt'),
            numParameters = 12,
            numLabels = 1,
            pattern = True,
        )
        self.groups['int15'] = ThermoGroups(label='int15', name='1,5-Interaction Corrections').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', '15_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', '15_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', '15_Library.txt'),
            numParameters = 12,
            numLabels = 1,
            pattern = True,
        )
        self.groups['radical'] = ThermoGroups(label='radical', name='Radical Corrections').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', 'Radical_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', 'Radical_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', 'Radical_Library.txt'),
            numParameters = 12,
            numLabels = 1,
            pattern = True,
        )
        self.groups['ring'] = ThermoGroups(label='ring', name='Ring Corrections').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', 'Ring_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', 'Ring_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', 'Ring_Library.txt'),
            numParameters = 12,
            numLabels = 1,
            pattern = True,
        )
        self.groups['polycyclic'] = ThermoGroups(label='other', name='Polycyclic Ring Corrections').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', 'Polycyclic_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', 'Polycyclic_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', 'Polycyclic_Library.txt'),
            numParameters = 12,
            numLabels = 1,
            pattern = True,
        )
        self.groups['other'] = ThermoGroups(label='other', name='Other Corrections').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', 'Other_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', 'Other_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', 'Other_Library.txt'),
            numParameters = 12,
            numLabels = 1,
            pattern = True,
        )
        
    def pruneHeteroatoms(self, allowed=['C','H','O','S']):
        """
        Remove all species from thermo libraries that contain atoms other than those allowed.
        
        This is useful before saving the database for use in RMG-Java
        """
        allowedElements = [rmgpy.molecule.element.getElement(label) for label in allowed]
        for library in self.libraries.values():
            logging.info("Removing hetoroatoms from thermo library '{0}'".format(library.name))
            toDelete = []
            for entry in library.entries.values():
                for atom in entry.item.atoms:
                    if atom.element not in allowedElements:
                        toDelete.append(entry.label)
                        break
            for label in toDelete:
                logging.info(" {0}".format(label))
                library.entries.pop(label)

    def saveOld(self, path):
        """
        Save the old RMG thermo database to the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """

        # Depository not used in old database, so it is not saved

        librariesPath = os.path.join(path, 'thermo_libraries')
        if not os.path.exists(librariesPath): os.mkdir(librariesPath)
        for library in self.libraries.values():
            libraryPath = os.path.join(librariesPath, library.label)
            if not os.path.exists(libraryPath): os.mkdir(libraryPath)
            library.saveOld(
                dictstr = os.path.join(libraryPath, 'Dictionary.txt'),
                treestr = '',
                libstr = os.path.join(libraryPath, 'Library.txt'),
            )

        groupsPath = os.path.join(path, 'thermo_groups')
        if not os.path.exists(groupsPath): os.mkdir(groupsPath)
        self.groups['group'].saveOld(
            dictstr = os.path.join(groupsPath, 'Group_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Group_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Group_Library.txt'),
        )
        self.groups['gauche'].saveOld(
            dictstr = os.path.join(groupsPath, 'Gauche_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Gauche_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Gauche_Library.txt'),
        )
        self.groups['int15'].saveOld(
            dictstr = os.path.join(groupsPath, '15_Dictionary.txt'),
            treestr = os.path.join(groupsPath, '15_Tree.txt'),
            libstr = os.path.join(groupsPath, '15_Library.txt'),
        )
        self.groups['radical'].saveOld(
            dictstr = os.path.join(groupsPath, 'Radical_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Radical_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Radical_Library.txt'),
        )
        self.groups['ring'].saveOld(
            dictstr = os.path.join(groupsPath, 'Ring_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Ring_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Ring_Library.txt'),
        )
        self.groups['polycyclic'].saveOld(
            dictstr = os.path.join(groupsPath, 'Polycyclic_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Polycyclic_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Polycyclic_Library.txt'),
        )
        self.groups['other'].saveOld(
            dictstr = os.path.join(groupsPath, 'Other_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Other_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Other_Library.txt'),
        )

    def recordPolycylicGenericNodes(self):
        """
        Identify generic nodes in tree for polycyclic groups.
        Saves them as a list in the `genericNodes` attribute
        in the polycyclic :class:`ThermoGroups` object, which
        must be pre-loaded.

        Necessary for polycyclic heuristic.
        """
        self.groups['polycyclic'].genericNodes = ['PolycyclicRing']
        for label, entry in self.groups['polycyclic'].entries.iteritems():
            if isinstance(entry.data, ThermoData): 
                continue
            self.groups['polycyclic'].genericNodes.append(label)

    def recordRingGenericNodes(self):
        """
        Identify generic nodes in tree for ring groups.
        Saves them as a list in the `genericNodes` attribute
        in the ring :class:`ThermoGroups` object, which
        must be pre-loaded.

        Necessary for polycyclic heuristic.
        """
        self.groups['ring'].genericNodes = ['Ring']
        for label, entry in self.groups['ring'].entries.iteritems():
            if isinstance(entry.data, ThermoData): 
                continue
            self.groups['ring'].genericNodes.append(label)

    def getThermoData(self, species, trainingSet=None):
        """
        Return the thermodynamic parameters for a given :class:`Species`
        object `species`. This function first searches the loaded libraries
        in order, returning the first match found, before falling back to
        estimation via machine learning and then group additivity.
        
        The method corrects for symmetry when the molecule uses machine
        learning or group additivity. Libraries and direct QM calculations
        are already corrected.
        
        Returns: ThermoData
        """
        from rmgpy.rmg.input import getInput
        
        thermo0 = self.getThermoDataFromLibraries(species)
        
        if thermo0 is not None:
            logging.debug("Found thermo for {0} in {1}".format(species.label,thermo0[0].comment.lower()))
            assert len(thermo0) == 3, "thermo0 should be a tuple at this point: (thermoData, library, entry)"
            thermo0 = thermo0[0]

            if species.containsSurfaceSite():
                thermo0 = self.correctBindingEnergy(thermo0, species)
            return thermo0

        if species.containsSurfaceSite():
            thermo0 = self.getThermoDataForSurfaceSpecies(species)
            thermo0 = self.correctBindingEnergy(thermo0, species)
            return thermo0

        try:
            quantumMechanics = getInput('quantumMechanics')
        except Exception:
            logging.debug('Quantum Mechanics DB could not be found.')
            quantumMechanics = None

        try:
            ml_estimator, ml_settings = getInput('MLEstimator')
        except Exception:
            logging.debug('ML estimator could not be found.')
            ml_estimator, ml_settings = None, None
            
        if quantumMechanics:
            original_molecule = species.molecule[0]
            if quantumMechanics.settings.onlyCyclics and not original_molecule.isCyclic():
                pass
            else: # try a QM calculation
                if original_molecule.getRadicalCount() > quantumMechanics.settings.maxRadicalNumber:
                    # Too many radicals for direct calculation: use HBI.
                    logging.info("{0} radicals on {1} exceeds limit of {2}. Using HBI method.".format(
                        original_molecule.getRadicalCount(),
                        species.label,
                        quantumMechanics.settings.maxRadicalNumber,
                        ))
                    
                    # Need to estimate thermo via each resonance isomer
                    thermo = []
                    for molecule in species.molecule:
                        molecule.clearLabeledAtoms()
                        # Try to see if the saturated molecule can be found in the libraries
                        tdata = self.estimateRadicalThermoViaHBI(molecule, self.getThermoDataFromLibraries)
                        priority = 1
                        if tdata is None:
                            # Then attempt quantum mechanics job on the saturated molecule
                            tdata = self.estimateRadicalThermoViaHBI(molecule, quantumMechanics.getThermoData)
                            priority = 2
                        if tdata is None:
                            # Fall back to group additivity
                            tdata = self.estimateThermoViaGroupAdditivity(molecule)
                            priority = 3
                        
                        thermo.append((priority, tdata.getEnthalpy(298.), molecule, tdata))
                    
                    if len(thermo) > 1:
                        # Sort thermo first by the priority, then by the most stable H298 value
                        thermo = sorted(thermo, key=lambda x: (x[0], x[1])) 
                        for i, therm in enumerate(thermo):
                            logging.debug("Resonance isomer {0} {1} gives H298={2:.0f} J/mol".format(i+1,
                                            therm[2].toSMILES(), therm[1]))
                        # Save resonance isomers reordered by their thermo
                        species.molecule = [item[2] for item in thermo]
                        original_molecule = species.molecule[0]
                    thermo0 = thermo[0][3] 

                    # update entropy by symmetry correction
                    thermo0.S298.value_si -= constants.R * math.log(species.getSymmetryNumber())
                    
                else: # Not too many radicals: do a direct calculation.
                    thermo0 = quantumMechanics.getThermoData(original_molecule) # returns None if it fails
                
        if thermo0 is None:
            # First try finding stable species in libraries and using HBI
            for mol in species.molecule:
                if mol.reactive:
                    original_molecule = mol
                    break
            else:
                for mol in species.molecule:
                    logging.info(mol.toAdjacencyList())
                    logging.info('reactive = {0}'.format(mol.reactive))
                    logging.info('\n')
                raise ValueError('Could not process a species with no reactive structures')
            if original_molecule.getRadicalCount() > 0:
                # If the molecule is a radical, check if any of the saturated forms are in the libraries
                # first and perform an HBI correction on them
                thermo = []
                for molecule in species.molecule:
                    if molecule.reactive:
                        molecule.clearLabeledAtoms()
                        # First see if the saturated molecule is in the libaries
                        tdata = self.estimateRadicalThermoViaHBI(molecule, self.getThermoDataFromLibraries)
                        if tdata:
                            thermo.append((tdata.getEnthalpy(298.), molecule, tdata))
                
                if thermo:
                    # Sort thermo by the most stable H298 value when choosing between thermoLibrary values
                    thermo = sorted(thermo, key=lambda x: x[0])
                    # Sort thermo by the structure reactive attribute, with `reactive=True` structures first
                    thermo.sort(key=lambda x: x[1].reactive, reverse=True)
                    for i, therm in enumerate(thermo):
                        if therm[1].reactive:
                            logging.debug("Resonance isomer {0} {1} gives H298={2:.0f} J/mol".format(i+1,
                                            therm[1].toSMILES(), therm[0]))
                        else:
                            logging.debug("Non-reactive resonance isomer {0} {1} gives H298={2:.0f} J/mol".format(i+1,
                                            therm[1].toSMILES(), therm[0]))
                    # Save resonance isomers reordered by their thermo
                    newMolList = [item[1] for item in thermo]
                    if len(newMolList) < len(species.molecule):
                        newMolList.extend([mol for mol in species.molecule if mol not in newMolList])
                    species.molecule = newMolList
                    thermo0 = thermo[0][2]

            if thermo0 is None:
                # If we still don't have thermo, use ML to estimate it, but
                # only if the molecule is made up of H, C, N, and O atoms and
                # is not a singlet carbene. ML settings are checked in
                # `self.get_thermo_data_from_ml`.
                if (ml_estimator is not None
                        and all(a.element.number in {1, 6, 7, 8} for a in species.molecule[0].atoms)
                        and species.molecule[0].getSingletCarbeneCount() == 0):

                    thermo0 = self.get_thermo_data_from_ml(species,
                                                           ml_estimator,
                                                           ml_settings)

            if thermo0 is None:
                # And lastly, resort back to group additivity to determine thermo for molecule
                thermo0 = self.getThermoDataFromGroups(species)

            # Update entropy by symmetry correction (not included in trained ML model)
            thermo0.S298.value_si -= constants.R * math.log(species.getSymmetryNumber())

        # Make sure to calculate Cp0 and CpInf if it wasn't done already
        findCp0andCpInf(species, thermo0)

        # Return the resulting thermo parameters
        return thermo0

    def setDeltaAtomicAdsorptionEnergies(self, bindingEnergies=None):
        """
        Sets and stores the change in atomic binding energy between
        the desired and the Pt(111) default.

        This depends on the two metal surfaces: the reference one used in
        the database of adsorption energies, and the desired surface.

        If bindingEnergies are not provided, resets the values to those
        of the Pt(111) default.

        Args:
            bindingEnergies (dict, optional): the desired binding energies with
                elements as keys and binding energy/unit tuples as values

        Returns:
            None, stores result in self.deltaAtomicAdsorptionEnergy
        """
        referenceBindingEnergies = {
            'C': rmgpy.quantity.Energy(-6.750, 'eV/molecule'),
            'H': rmgpy.quantity.Energy(-2.479, 'eV/molecule'),
            'O': rmgpy.quantity.Energy(-3.586, 'eV/molecule'),
            'N': rmgpy.quantity.Energy(-4.352, 'eV/molecule'),
        }

        # Use Pt(111) reference if no binding energies are provided
        if bindingEnergies is None:
            bindingEnergies = referenceBindingEnergies

        self.deltaAtomicAdsorptionEnergy = {
            'C': rmgpy.quantity.Energy(0.0, 'eV/molecule'),
            'H': rmgpy.quantity.Energy(0.0, 'eV/molecule'),
            'O': rmgpy.quantity.Energy(0.0, 'eV/molecule'),
            'N': rmgpy.quantity.Energy(0.0, 'eV/molecule'),
        }

        for element, deltaEnergy in self.deltaAtomicAdsorptionEnergy.iteritems():
            deltaEnergy.value_si = bindingEnergies[element].value_si - referenceBindingEnergies[element].value_si

    def correctBindingEnergy(self, thermo, species):
        """
        Changes the provided thermo, by applying a linear scaling relation
        to correct the adsorption energy.

        :param thermo: starting thermo data
        :param species: the species (which is an adsorbate)
        :return: corrected thermo
        """
        molecule = species.molecule[0]
        # only want/need to do one resonance structure
        surfaceSites = []
        for atom in molecule.atoms:
            if atom.isSurfaceSite():
                surfaceSites.append(atom)
        normalizedBonds = {'C':0., 'O':0., 'N':0., 'H':0.}
        maxBondOrder = {'C':4., 'O':2., 'N':3., 'H':1.}
        for site in surfaceSites:
            numbonds = len(site.bonds)
            if numbonds == 0:
                #vanDerWaals
                pass
            else:
                assert len(site.bonds) == 1, "Each surface site can only be bonded to 1 atom"
                bondedAtom = site.bonds.keys()[0]
                bond = site.bonds[bondedAtom]
                if bond.isSingle():
                    bondOrder = 1.
                elif bond.isDouble():
                    bondOrder = 2.
                elif bond.isTriple():
                    bondOrder = 3.
                elif bond.isQuadruple():
                    bondOrder = 4.
                else:
                    raise NotImplementedError("Unsupported bond order {0} for binding energy correction.".format(bond.order))

                normalizedBonds[bondedAtom.symbol] += bondOrder / maxBondOrder[bondedAtom.symbol]

        if not isinstance(thermo, ThermoData):
            thermo = thermo.toThermoData()
            findCp0andCpInf(species, thermo)

        ## now edit the adsorptionThermo using LSR
        for element in 'CHON':
            changeInBindingEnergy = self.deltaAtomicAdsorptionEnergy[element].value_si * normalizedBonds[element]
            thermo.H298.value_si += changeInBindingEnergy
        thermo.comment += " Binding energy corrected by LSR."
        return thermo

    def getThermoDataForSurfaceSpecies(self, species):
        """
        Get the thermo data for an adsorbed species,
        by desorbing it, finding the thermo of the gas-phase
        species, then adding an adsorption correction that
        is found from the groups/adsorption tree.
        Does not apply linear scaling relationship.
        
        Returns a :class:`ThermoData` object, with no Cp0 or CpInf
        """

        if species.isSurfaceSite():
            raise DatabaseError("Can't estimate thermo of vacant site. Should be in library (and should be 0).")

        logging.debug(("Trying to generate thermo for surface species"
                        " with these {} resonance isomer(s):").format(len(species.molecule)))
        molecule = species.molecule[0]
        # only want/need to do one resonance structure,
        # because will need to regenerate others in gas phase
        dummyMolecule = molecule.copy(deep=True)
        sitesToRemove = []
        for atom in dummyMolecule.atoms:
            if atom.isSurfaceSite():
                sitesToRemove.append(atom)
        for site in sitesToRemove:
            numbonds = len(site.bonds)
            if numbonds == 0:
                #vanDerWaals
                pass
            else:
                assert len(site.bonds) == 1, "Each surface site can only be bonded to 1 atom"
                bondedAtom = site.bonds.keys()[0]
                bond = site.bonds[bondedAtom]
                dummyMolecule.removeBond(bond)
                if bond.isSingle():
                    bondedAtom.incrementRadical()
                elif bond.isDouble():
                    bondedAtom.incrementRadical()
                    bondedAtom.incrementRadical()
                elif bond.isTriple():
                    bondedAtom.incrementRadical()
                    bondedAtom.incrementLonePairs()
                elif bond.isQuadruple():
                    bondedAtom.incrementRadical()
                    bondedAtom.incrementRadical()
                    bondedAtom.incrementLonePairs()
                else:
                    raise NotImplementedError("Can't remove surface bond of type {}".format(bond.order))

            dummyMolecule.removeAtom(site)
        dummyMolecule.update()

        logging.debug("Before removing from surface:\n" + molecule.toAdjacencyList())
        logging.debug("After removing from surface:\n" + dummyMolecule.toAdjacencyList())

        dummySpecies = Species()
        dummySpecies.molecule.append(dummyMolecule)
        dummySpecies.generate_resonance_structures()
        thermo = self.getThermoData(dummySpecies)

        thermo.comment = "Gas phase thermo from {0}. Adsorption correction:".format(thermo.comment)
        logging.debug("Using thermo from gas phase for species {}\n".format(species.label) + repr(thermo))

        if not isinstance(thermo, ThermoData):
            thermo = thermo.toThermoData()
            findCp0andCpInf(species, thermo)

        ## Get the adsorption energy
        # Create the ThermoData object
        adsorptionThermo = ThermoData(
            Tdata = ([300,400,500,600,800,1000,1500],"K"),
            Cpdata = ([0.0,0.0,0.0,0.0,0.0,0.0,0.0],"J/(mol*K)"),
            H298 = (0.0,"kJ/mol"),
            S298 = (0.0,"J/(mol*K)"),
        )
        try:
            self.__addGroupThermoData(adsorptionThermo, self.groups['adsorptionPt'], molecule, {})
        except KeyError:
            logging.error("Couldn't find in adsorption thermo database:")
            logging.error(molecule)
            logging.error(molecule.toAdjacencyList())
            raise

        # (groupAdditivity=True means it appends the comments)
        addThermoData(thermo, adsorptionThermo, groupAdditivity=True)

        if thermo.label:
            thermo.label += 'X'

        findCp0andCpInf(species, thermo)
        return thermo
        
    def getThermoDataFromLibraries(self, species, trainingSet=None):
        """
        Return the thermodynamic parameters for a given :class:`Species`
        object `species`. This function first searches the loaded libraries
        in order, returning the first match found, before failing and returning None.
        `trainingSet` is used to identify if function is called during training set or not.
        During training set calculation we want to use gas phase thermo to not affect reverse
        rate calculation.
        
        Returns: ThermoData or None
        """
        import rmgpy.rmg.main
        thermoData = None
        
        #chatelak 11/15/14: modification to introduce liquid phase thermo libraries
        libraryList=deepcopy(self.libraryOrder) #copy the value to not affect initial object

        if rmgpy.rmg.main.solvent is not None:
            liqLibraries=[]
            #Liquid phase simulation part: 
            #This bloc "for": Identify liquid phase libraries and store them in liqLibraries
            for iterLib in libraryList:
                if self.libraries[iterLib].solvent:
                    liqLibraries.append(iterLib)
            #Check in liqLibraries if thermo for species exists and return the first match. Only if function not called by trainingSet
            if liqLibraries and trainingSet is None:
                for label in liqLibraries:
                    thermoData = self.getThermoDataFromLibrary(species, self.libraries[label])
                    if thermoData is not None:
                        assert len(thermoData) == 3, "thermoData should be a tuple at this point"
                        #Watch out comments changed: this is used later to apply solvation or not on species matching thermo. If required, Modify this carefully.
                        thermoData[0].comment += 'Liquid thermo library: ' + label
                        return thermoData
            #Remove liqLibraries from libraryList if: called by training set (trainingSet=True) or if no thermo found in liqLibrairies
            #if no liquid library found this does nothing.   
            for libIter in liqLibraries:
                libraryList.remove(libIter)

        # Condition to execute this part: gas phase simulation or training set or liquid phase simulation with : noliquid libraries found or no matching species found in liquid libraries       
        # If gas phase simulation libraryList = self.libraryOrder (just like before modifications) and they are all gas phase, already checked by checkLibrairies function in database.load()
        # Check the libraries in order; return the first successful match
        for label in libraryList:
            thermoData = self.getThermoDataFromLibrary(species, self.libraries[label])
            if thermoData is not None:
                assert len(thermoData) == 3, "thermoData should be a tuple at this point"
                if rmgpy.rmg.main.solvent is not None and trainingSet is None:
                    thermoData[0].comment += 'Thermo library corrected for liquid phase: ' + label
                else:
                    thermoData[0].comment += 'Thermo library: ' + label
                return thermoData

        return None                
                
    def getAllThermoData(self, species):
        """
        Return all possible sets of thermodynamic parameters for a given
        :class:`Species` object `species`. The hits from the depository come
        first, then the libraries (in order), and then the group additivity
        estimate. This method is useful for a generic search job.
        
        Returns: a list of tuples (ThermoData, source, entry) 
        (Source is a library or depository, or None)
        """
        thermoDataList = []
        # Data from depository comes first
        thermoDataList.extend(self.getThermoDataFromDepository(species))
        # Data from libraries comes second
        for label in self.libraryOrder:
            data = self.getThermoDataFromLibrary(species, self.libraries[label])
            if data: 
                assert len(data) == 3, "thermoData should be a tuple at this point"
                data[0].comment += label
                thermoDataList.append(data)

        # Last entry is always the estimate from group additivity
        # Make it a tuple
        # Distinguish surface species, as orignial getThermoDataFromGroups does
        # not work for surface sites or surface species
        if species.isSurfaceSite():
            # Cannot estimate thermo of vacant site. Thermo stores in library
            pass
        elif species.containsSurfaceSite():
            try:
                # Estimate thermo of surface species based on modfied GA method
                data = (self.getThermoDataForSurfaceSpecies(species), None, None)
            except DatabaseError:
                 # We don't have a GAV estimate, e.g. unsupported element
                pass
            else:
                thermoDataList.append(data)
        else:
            try:
                data = (self.getThermoDataFromGroups(species), None, None)
            except DatabaseError:
                # We don't have a GAV estimate, e.g. unsupported element
                pass
            else:
                # update group activity for symmetry
                data[0].S298.value_si -= constants.R * math.log(species.getSymmetryNumber())
                thermoDataList.append(data)
        # Return all of the resulting thermo parameters
        return thermoDataList

    def getThermoDataFromDepository(self, species):
        """
        Return all possible sets of thermodynamic parameters for a given
        :class:`Species` object `species` from the depository. If no
        depository is loaded, a :class:`DatabaseError` is raised.
        
        Returns: a list of tuples (thermoData, depository, entry) without any Cp0 or CpInf data.
        """
        items = []
        for entry in self.depository['stable'].entries.itervalues():
            for molecule in species.molecule:
                if molecule.isIsomorphic(entry.item):
                    items.append((deepcopy(entry.data), self.depository['stable'], entry))
                    break
        for entry in self.depository['radical'].entries.itervalues():
            for molecule in species.molecule:
                if molecule.isIsomorphic(entry.item):
                    items.append((deepcopy(entry.data), self.depository['radical'], entry))
                    break
        return items

    def getThermoDataFromLibrary(self, species, library):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Species` object `species` from the specified thermodynamics
        `library`. If `library` is a string, the list of libraries is searched
        for a library with that name. If no match is found in that library,
        ``None`` is returned. If no corresponding library is found, a
        :class:`DatabaseError` is raised.
        
        Returns a tuple: (ThermoData, library, entry)  or None.
        """
        match = None
        for entry in library.entries.itervalues():
            for molecule in species.molecule:
                if molecule.isIsomorphic(entry.item) and entry.data is not None:
                    thermoData = deepcopy(entry.data)
                    thermoData.label = entry.label
                    findCp0andCpInf(species, thermoData)
                    match = (thermoData, library, entry)
                    break
            if match is not None:
                break
        if match is not None:
            # Move the matched molecule to the first position in the list
            species.molecule.remove(molecule)
            species.molecule.insert(0, molecule)
        return match

    def getThermoDataFromGroups(self, species):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Species` object `species` by estimation using the group
        additivity values. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        
        The resonance isomer (molecule) with the lowest H298 is used, and as a side-effect
        the resonance isomers (items in `species.molecule` list) are sorted in ascending order.
        
        This does not account for symmetry. The method calling this sould correct for it.
        
        Returns: ThermoData
        """       
        thermo = []
        for molecule in species.molecule:
            molecule.clearLabeledAtoms()
            molecule.updateAtomTypes()
            tdata = self.estimateThermoViaGroupAdditivity(molecule)
            thermo.append(tdata)

        indices = self.prioritizeThermo(species, thermo)
        
        species.molecule = [species.molecule[ind] for ind in indices]
        
        thermoData = thermo[indices[0]]
        findCp0andCpInf(species, thermoData)
        return thermoData

    def get_thermo_data_from_ml(self, species, ml_estimator, ml_settings):
        """
        Return the set of thermodynamic parameters corresponding to a
        given :class:`Species` object `species` by estimation using the
        ML estimator. Also compare the estimated uncertainties to the
        user-defined cutoffs. If any of the uncertainties are larger
        than their corresponding cutoffs, return None. Also check all
        other options in `ml_settings`.

        For HBI, the resonance isomer with the lowest H298 is used and
        the resonance isomers in species are sorted in ascending order.

        The entropy is not corrected for the symmetry of the molecule.
        This should be done later by the calling function.
        """
        molecule = species.molecule[0]

        min_heavy = ml_settings['min_heavy_atoms'] or 1
        max_heavy = ml_settings['max_heavy_atoms'] or numpy.inf
        min_carbon = ml_settings['min_carbon_atoms'] or 0
        max_carbon = ml_settings['max_carbon_atoms'] or numpy.inf
        min_oxygen = ml_settings['min_oxygen_atoms'] or 0
        max_oxygen = ml_settings['max_oxygen_atoms'] or numpy.inf
        min_nitrogen = ml_settings['min_nitrogen_atoms'] or 0
        max_nitrogen = ml_settings['max_nitrogen_atoms'] or numpy.inf

        element_count = molecule.get_element_count()
        n_heavy = sum(count for element, count in element_count.iteritems() if element != 'H')

        if not(min_heavy <= n_heavy <= max_heavy):
            return None
        if not(min_carbon <= element_count.get('C', 0) <= max_carbon):
            return None
        if not(min_oxygen <= element_count.get('O', 0) <= max_oxygen):
            return None
        if not(min_nitrogen <= element_count.get('N', 0) <= max_nitrogen):
            return None
        if ml_settings['only_heterocyclics'] and not molecule.isHeterocyclic():
            return None
        if ml_settings['only_cyclics'] and not molecule.isCyclic():
            return None
        min_cycle_overlap = ml_settings['min_cycle_overlap']
        if min_cycle_overlap > 0 and molecule.getMaxCycleOverlap() < min_cycle_overlap:
            return None

        if molecule.isRadical():
            thermo = [self.estimateRadicalThermoViaHBI(mol, ml_estimator.get_thermo_data) for mol in species.molecule]
            H298 = numpy.array([tdata.H298 for tdata in thermo])
            indices = H298.argsort()
            species.molecule = [species.molecule[ind] for ind in indices]
            thermo0 = thermo[indices[0]]
        else:
            thermo0 = ml_estimator.get_thermo_data_for_species(species)

        # The keys for this dictionary should match the keys in
        # `RMG.ml_uncertainty_cutoffs`. Use a temperature-weighted
        # average to estimate uncertainty for Cp.
        uncertainties = dict(
            H298=thermo0.H298.uncertainty_si,
            S298=thermo0.S298.uncertainty_si,
            Cp=numpy.average(thermo0.Cpdata.uncertainty_si, weights=thermo0.Tdata.value_si)
        )

        ml_uncertainty_cutoffs = ml_settings['uncertainty_cutoffs']
        if any(uncertainties[p] > ml_uncertainty_cutoffs[p].value_si for p in ml_uncertainty_cutoffs):
            return None
        else:
            return thermo0


    def prioritizeThermo(self, species, thermoDataList):
        """
        Use some metrics to reorder a list of thermo data from best to worst.
        Return a list of indices with the desired order associated with the index of thermo from the data list.
        """
        if len(species.molecule) > 1:
            # Go further only if there is more than one isomer
            if species.molecule[0].isCyclic():
                # Special treatment for cyclic compounds
                entries = []
                for thermo in thermoDataList:
                    ringGroups, polycyclicGroups = self.getRingGroupsFromComments(thermo)
                    
                    # Use rank as a metric for prioritizing thermo. 
                    # The smaller the rank, the better.
                    sumRank = numpy.sum([3 if entry.rank is None else entry.rank for entry in ringGroups + polycyclicGroups])
                    entries.append((thermo, sumRank))
                
                # Sort first by rank, then by enthalpy at 298 K
                entries = sorted(entries, key=lambda entry: (entry[1], entry[0].getEnthalpy(298.)))
                indices = [thermoDataList.index(entry[0]) for entry in entries]
            else:
                # For noncyclics, default to original algorithm of ordering thermo based on the most stable enthalpy
                H298 = numpy.array([t.getEnthalpy(298.) for t in thermoDataList])
                indices = H298.argsort()
            indices = numpy.array([i for i in indices if species.molecule[i].reactive] +\
                        [i for i in indices if not species.molecule[i].reactive])
        else:
            indices = [0]
        return indices

    def estimateRadicalThermoViaHBI(self, molecule, stableThermoEstimator):
        """
        Estimate the thermodynamics of a radical by saturating it,
        applying the provided stableThermoEstimator method on the saturated species,
        then applying hydrogen bond increment corrections for the radical
        site(s) and correcting for the symmetry.
        
        No entropy is included in the returning term.
        This should be done later by the calling function.
        """
        
        assert molecule.isRadical(), "Method only valid for radicals."
        saturatedStruct = molecule.copy(deep=True)
        added = saturatedStruct.saturate_radicals()
        saturatedStruct.props['saturated'] = True
        
        # Get thermo estimate for saturated form of structure
        if stableThermoEstimator == self.getThermoDataFromLibraries:
            # Get data from libraries
            saturatedSpec = Species(molecule=[saturatedStruct])
            thermoData_sat = stableThermoEstimator(saturatedSpec)
            if thermoData_sat:
                assert len(thermoData_sat) == 3, "thermoData should be a tuple at this point: (thermoData, library, entry)"
                thermoData_sat = thermoData_sat[0]
        else:
            thermoData_sat = stableThermoEstimator(saturatedStruct)
        if thermoData_sat is None:
            # logging.info("Thermo data of saturated {0} of molecule {1} is None.".format(saturatedStruct, molecule))
            return None
        assert thermoData_sat is not None, "Thermo data of saturated {0} of molecule {1} is None!".format(saturatedStruct, molecule)
        
        # Convert to ThermoData object if necessary in order to add and subtract from enthalpy and entropy values
        if not isinstance(thermoData_sat, ThermoData):
            thermoData_sat = thermoData_sat.toThermoData()
        
        if not (stableThermoEstimator == self.computeGroupAdditivityThermo
                or isinstance(stableThermoEstimator.__self__, MLEstimator)):
            #remove the symmetry contribution to the entropy of the saturated molecule
            ##assumes that the thermo data comes from QMTP or from a thermolibrary
            thermoData_sat.S298.value_si += constants.R * math.log(saturatedStruct.getSymmetryNumber())
        
        thermoData = thermoData_sat
        
        # For each radical site, get radical correction
        # Only one radical site should be considered at a time; all others
        # should be saturated with hydrogen atoms
        for atom in added:
            # Remove the added hydrogen atoms and bond and restore the radical
            for H, bond in added[atom]:
                saturatedStruct.removeBond(bond)
                saturatedStruct.removeAtom(H)
                atom.incrementRadical()
            saturatedStruct.update()
            try:
                self.__addGroupThermoData(thermoData, self.groups['radical'], saturatedStruct, {'*':atom})
            except KeyError:
                logging.error("Couldn't find in radical thermo database:")
                logging.error(molecule)
                logging.error(molecule.toAdjacencyList())
                raise
            # Re-saturate
            for H, bond in added[atom]:
                saturatedStruct.addAtom(H)
                saturatedStruct.addBond(bond)
                atom.decrementRadical()
            # Subtract the enthalpy of the added hydrogens
            for H, bond in added[atom]:
                thermoData.H298.value_si -= 52.103 * 4184

        # Remove all of the interactions of the saturated structure. Then add the interactions of the radical.
        # Take C1=CC=C([O])C(O)=C1 as an example, we need to remove the interation of OH-OH, then add the interaction of Oj-OH.
        # For now, we only apply this part to cyclic structure because we only have radical interaction data for aromatic radical.
        if saturatedStruct.isCyclic():
            SSSR = saturatedStruct.getSmallestSetOfSmallestRings()
            for ring in SSSR:
                for atomPair in itertools.permutations(ring, 2):
                    try:
                        self.__removeGroupThermoData(thermoData,self.groups['longDistanceInteraction_cyclic'], saturatedStruct, {'*1':atomPair[0], '*2':atomPair[1]})
                    except KeyError: pass
            SSSR = molecule.getSmallestSetOfSmallestRings()
            for ring in SSSR:
                for atomPair in itertools.permutations(ring, 2):
                    try:
                        self.__addGroupThermoData(thermoData,self.groups['longDistanceInteraction_cyclic'], molecule, {'*1':atomPair[0], '*2':atomPair[1]})
                    except KeyError: pass
        
        thermoData.label = '' #prevents the original thermo species name being used for the HBI corrected radical in species generation
        
        return thermoData
        
        
    def estimateThermoViaGroupAdditivity(self, molecule):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Molecule` object `molecule` by estimation using the group
        additivity values. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        
        The entropy is not corrected for the symmetry of the molecule.
        This should be done later by the calling function.
        """
        # For thermo estimation we need the atoms to already be sorted because we
        # iterate over them; if the order changes during the iteration then we
        # will probably not visit the right atoms, and so will get the thermo wrong
        molecule.sortAtoms()

        if molecule.isRadical():
            thermoData = self.estimateRadicalThermoViaHBI(molecule, self.computeGroupAdditivityThermo)
        else:
            thermoData = self.computeGroupAdditivityThermo(molecule)
        return thermoData


    def computeGroupAdditivityThermo(self, molecule):
        """
        Return the set of thermodynamic parameters corresponding to a given
        :class:`Molecule` object `molecule` by estimation using the group
        additivity values. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        
        The entropy is not corrected for the symmetry of the molecule.
        This should be done later by the calling function.
        """

        assert not molecule.isRadical(), "This method is only for saturated non-radical species."
        # For thermo estimation we need the atoms to already be sorted because we
        # iterate over them; if the order changes during the iteration then we
        # will probably not visit the right atoms, and so will get the thermo wrong
        molecule.sortAtoms()

        # Create the ThermoData object
        thermoData = ThermoData(
            Tdata = ([300,400,500,600,800,1000,1500],"K"),
            Cpdata = ([0.0,0.0,0.0,0.0,0.0,0.0,0.0],"J/(mol*K)"),
            H298 = (0.0,"kJ/mol"),
            S298 = (0.0,"J/(mol*K)"),
        )

        cyclic = molecule.isCyclic()
        # Generate estimate of thermodynamics
        for atom in molecule.atoms:
            # Iterate over heavy (non-hydrogen) atoms
            if atom.isNonHydrogen():
                # Get initial thermo estimate from main group database
                try:
                    self.__addGroupThermoData(thermoData, self.groups['group'], molecule, {'*':atom})
                except KeyError:
                    logging.error("Couldn't find in main thermo database:")
                    logging.error(molecule)
                    logging.error(molecule.toAdjacencyList())
                    raise
                # Correct for gauche and 1,5- interactions
                # Pair atom with its 1st and 2nd nonHydrogen neighbors, 
                # Then match the pair with the entries in the database longDistanceInteraction_noncyclic.py
                # Currently we only have gauche(1,4) and 1,5 interactions in that file. 
                # If you want to add more corrections for longer distance, please call getNthNeighbor() method accordingly.
                # Potentially we could include other.py in this database, but it's a little confusing how to label atoms for the entries in other.py
                if not molecule.isAtomInCycle(atom):
                    for atom_2 in molecule.getNthNeighbor([atom],[1,2]):
                        if not molecule.isAtomInCycle(atom_2):
                        # This is the correction for noncyclic structure. If `atom` or `atom_2` is in a cycle, do not apply this correction.
                        # Note that previously we do not do gauche for cyclic molecule, which is unreasonable for cyclic molecule with a long tail.
                            try:
                                self.__addGroupThermoData(thermoData, self.groups['longDistanceInteraction_noncyclic'], molecule, {'*1':atom, '*2': atom_2})
                            except KeyError: pass
                try:
                    self.__addGroupThermoData(thermoData, self.groups['other'], molecule, {'*':atom})
                except KeyError: pass
        
        # Do long distance interaction correction for cyclic molecule. 
        # First get smallest set of smallest rings. 
        # Then for every single ring, generate the atom pairs by itertools.permutation.
        # Finally match the atom pair with the database.
        # WIPWIPWIPWIPWIPWIPWIP         #########################################         WIPWIPWIPWIPWIPWIPWIP
        # WIP: For now, in the database, if an entry describes the interaction between same groups, 
        # it will be halved because it will be counted twice here. 
        # Alternatively we could keep all the entries as their full values by using combinations instead of permutations here.
        # In that case, we need to add more lines to match from reverse side when we didn't hit the most specific level from the forward side.
        # PS: by saying 'forward side', I mean {'*1':atomPair[0], '*2':atomPair[1]}. So the following is the reverse side '{'*1':atomPair[1], '*2':atomPair[0]}'
        # In my opinion, it's cleaner to do it in the current way.
        # WIPWIPWIPWIPWIPWIPWIP         #########################################         WIPWIPWIPWIPWIPWIPWIP
        if cyclic:
            SSSR = molecule.getSmallestSetOfSmallestRings()
            for ring in SSSR:
                for atomPair in itertools.permutations(ring, 2):
                    try:
                        self.__addGroupThermoData(thermoData,self.groups['longDistanceInteraction_cyclic'], molecule, {'*1':atomPair[0], '*2':atomPair[1]})
                    except KeyError: pass

        # Do ring corrections separately because we only want to match
        # each ring one time
        
        if cyclic:                
            monorings, polyrings = molecule.getDisparateRings()
            for ring in monorings:
                # Make a temporary structure containing only the atoms in the ring
                # NB. if any of the ring corrections depend on ligands not in the ring, they will not be found!
                try:
                    self.__addRingCorrectionThermoDataFromTree(thermoData, self.groups['ring'], molecule, ring)
                except KeyError:
                    logging.error("Couldn't find a match in the monocyclic ring database even though monocyclic rings were found.")
                    logging.error(molecule)
                    logging.error(molecule.toAdjacencyList())
                    raise
            for polyring in polyrings:
                # Make a temporary structure containing only the atoms in the ring
                # NB. if any of the ring corrections depend on ligands not in the ring, they will not be found!
                try:
                    self.__addPolycyclicCorrectionThermoData(thermoData, molecule, polyring)
                except KeyError:
                    logging.error("Couldn't find a match in the polycyclic ring database even though polycyclic rings were found.")
                    logging.error(molecule)
                    logging.error(molecule.toAdjacencyList())
                    raise

        return thermoData

    def __addPolycyclicCorrectionThermoData(self, thermoData, molecule, polyring):
        """
        INPUT: `polyring` as a list of `Atom` forming a polycyclic ring
        OUTPUT: if the input `polyring` can be fully matched in polycyclic database, the correction
        will be directly added to `thermoData`; otherwise, a heuristic approach will 
        be applied.
        """
        # look up polycylic tree directly
        matched_group_thermodata, matched_group, isPartialMatch = self.__addRingCorrectionThermoDataFromTree(None, self.groups['polycyclic'], molecule, polyring)
        
        # if partial match (non-H atoms number same between 
        # polycylic ring in molecule and match group)
        # otherwise, apply heuristic algorithm
        if not isPartialMatch:
            if isBicyclic(polyring) and matched_group.label in self.groups['polycyclic'].genericNodes:
                # apply secondary decompostion formula
                # to get a estimated_group_thermodata
                estimated_bicyclic_thermodata = self.getBicyclicCorrectionThermoDataFromHeuristic(polyring)
                if not estimated_bicyclic_thermodata:
                    estimated_bicyclic_thermodata = matched_group_thermodata
                thermoData = addThermoData(thermoData, estimated_bicyclic_thermodata, groupAdditivity=True, verbose=True)
            else:
                # keep matched_group_thermodata as is
                thermoData = addThermoData(thermoData, matched_group_thermodata, groupAdditivity=True, verbose=True)
                # By setting verbose=True, we turn on the comments of polycyclic correction to pass the unittest.
                # Typically this comment is very short and also very helpful to check if the ring correction is calculated correctly.
        else:
            self.__addPolyRingCorrectionThermoDataFromHeuristic(thermoData, polyring)
            

    def __addPolyRingCorrectionThermoDataFromHeuristic(self, thermoData, polyring):
        """
        INPUT: `polyring` as a list of `Atom` forming a polycyclic ring, which can 
        only be partially matched.
        OUTPUT: `polyring` will be decomposed into a combination of 2-ring polycyclics
        and each one will be looked up from polycyclic database. The heuristic formula 
        is "polyring thermo correction = sum of correction of all 2-ring sub-polycyclics - 
        overlapped single-ring correction"; the calculated polyring thermo correction 
        will be finally added to input `thermoData`.
        """

        # pring decomposition
        bicyclicsMergedFromRingPair, ringOccurancesDict = bicyclicDecompositionForPolyring(polyring)
        
        # loop over 2-ring cores
        for bicyclic in bicyclicsMergedFromRingPair:
            matched_group_thermodata, matched_group, _ = self.__addRingCorrectionThermoDataFromTree(None, 
                self.groups['polycyclic'], bicyclic, bicyclic.atoms)

            if matched_group.label in self.groups['polycyclic'].genericNodes:
                # apply secondary decompostion formula
                # to get a estimated_group_thermodata
                estimated_bicyclic_thermodata = self.getBicyclicCorrectionThermoDataFromHeuristic(bicyclic.atoms)
                if not estimated_bicyclic_thermodata:
                    estimated_bicyclic_thermodata = matched_group_thermodata 
                thermoData = addThermoData(thermoData, estimated_bicyclic_thermodata, groupAdditivity=True, verbose=True)
            else:
                # keep matched_group_thermodata as is
                thermoData = addThermoData(thermoData, matched_group_thermodata, groupAdditivity=True, verbose=True)

        # loop over 1-ring 
        for singleRingTuple, occurance in ringOccurancesDict.iteritems():
            singleRing = list(singleRingTuple)

            if occurance >= 2:
                submol, _ = convertRingToSubMolecule(singleRing)
                
                if not isAromaticRing(submol):
                    aromaticBonds = findAromaticBondsFromSubMolecule(submol)
                    for aromaticBond in aromaticBonds:
                        aromaticBond.setOrderNum(1)

                    submol.saturate_unfilled_valence()
                    singleRingThermodata = self.__addRingCorrectionThermoDataFromTree(None, \
                                                self.groups['ring'], submol, submol.atoms)[0]
                    
                else:
                    submol.update()
                    singleRingThermodata = self.__addRingCorrectionThermoDataFromTree(None, \
                                                    self.groups['ring'], submol, submol.atoms)[0]
            for _ in range(occurance-1):
                thermoData = removeThermoData(thermoData, singleRingThermodata, True, True)
                # By setting verbose=True, we turn on the comments of polycyclic correction to pass the unittest.
                # Typically this comment is very short and also very helpful to check if the ring correction is calculated correctly.


    def getBicyclicCorrectionThermoDataFromHeuristic(self, bicyclic):

        # saturate if the bicyclic has unsaturated bonds
        # otherwise return None
        bicyclic_submol = convertRingToSubMolecule(bicyclic)[0]
        saturated_bicyclic_submol, alreadySaturated = saturate_ring_bonds(bicyclic_submol)

        if alreadySaturated:
            return None
        # split bicyclic into two single ring submols
        single_ring_submols = splitBicyclicIntoSingleRings(bicyclic_submol)

        # split saturated bicyclic into two single ring submols
        saturated_single_ring_submols = splitBicyclicIntoSingleRings(saturated_bicyclic_submol)

        # apply formula: 
        # bicyclic correction ~= saturated bicyclic correction - 
        # saturated single ring corrections + single ring corrections

        estimated_bicyclic_thermodata = ThermoData(
            Tdata = ([300,400,500,600,800,1000,1500],"K"),
            Cpdata = ([0.0,0.0,0.0,0.0,0.0,0.0,0.0],"J/(mol*K)"),
            H298 = (0.0,"kJ/mol"),
            S298 = (0.0,"J/(mol*K)")
        )
        
        saturated_bicyclic_thermoData = self.__addRingCorrectionThermoDataFromTree(None, 
                self.groups['polycyclic'], saturated_bicyclic_submol, saturated_bicyclic_submol.atoms)[0]
        
        estimated_bicyclic_thermodata = addThermoData(estimated_bicyclic_thermodata,
                                                    saturated_bicyclic_thermoData , 
                                                    groupAdditivity=True)

        estimated_bicyclic_thermodata.comment = "Estimated bicyclic component: " + saturated_bicyclic_thermoData.comment

        for submol in saturated_single_ring_submols:
            
            if not isAromaticRing(submol):
                aromaticBonds = findAromaticBondsFromSubMolecule(submol)
                for aromaticBond in aromaticBonds:
                    aromaticBond.setOrderNum(1)

                submol.saturate_unfilled_valence()
                single_ring_thermoData = self.__addRingCorrectionThermoDataFromTree(None,
                                            self.groups['ring'], submol, submol.atoms)[0]
                
            else:
                submol.update()
                single_ring_thermoData = self.__addRingCorrectionThermoDataFromTree(None,
                                                self.groups['ring'], submol, submol.atoms)[0]
            estimated_bicyclic_thermodata = removeThermoData(estimated_bicyclic_thermodata,
                                                    single_ring_thermoData, groupAdditivity=True, verbose=True)

        for submol in single_ring_submols:

            if not isAromaticRing(submol):
                aromaticBonds = findAromaticBondsFromSubMolecule(submol)
                for aromaticBond in aromaticBonds:
                    aromaticBond.setOrderNum(1)

                submol.saturate_unfilled_valence()
                single_ring_thermoData = self.__addRingCorrectionThermoDataFromTree(None,
                                            self.groups['ring'], submol, submol.atoms)[0]
                
            else:
                submol.update()
                single_ring_thermoData = self.__addRingCorrectionThermoDataFromTree(None,
                                                self.groups['ring'], submol, submol.atoms)[0]
                
            estimated_bicyclic_thermodata = addThermoData(estimated_bicyclic_thermodata,
                                                    single_ring_thermoData, groupAdditivity=True, verbose=True)

        return estimated_bicyclic_thermodata

    def __addRingCorrectionThermoDataFromTree(self, thermoData, ring_database, molecule, ring):
        """
        Determine the ring correction group additivity thermodynamic data for the given
         `ring` in the `molecule`, and add it to the existing thermo data
        `thermoData`.
        Also returns the matched ring group from the database from which the data originated.
        """
        matchedRingEntries = []
        # label each atom in the ring individually to try to match the group
        # for each ring, save only the ring that is matches the most specific leaf in the tree.
        for atom in ring:
            atoms = {'*':atom}
            entry = ring_database.descendTree(molecule, atoms)
            matchedRingEntries.append(entry)
        
        if matchedRingEntries is []:
            raise KeyError('Node not found in database.')
        # Decide which group to keep
        isPartialMatch = True
        completeMatchedGroups = [entry for entry in matchedRingEntries if not isRingPartialMatched(ring, entry.item)]

        if completeMatchedGroups:
            isPartialMatch = False
            matchedRingEntries = completeMatchedGroups

        depthList = [len(ring_database.ancestors(entry)) for entry in matchedRingEntries]
        mostSpecificMatchIndices = [i for i, x in enumerate(depthList) if x == max(depthList)]
        
        mostSpecificMatchedEntries = [matchedRingEntries[idx] for idx in mostSpecificMatchIndices]
        if len(set(mostSpecificMatchedEntries)) != 1:
            logging.debug('More than one type of node was found to be most specific for this ring.')
            logging.debug('This is either due to a database error in the ring or polycyclic groups, or a partial match between the group and the full ring.')
            logging.debug(mostSpecificMatchedEntries)
            
        # Condense the number of most specific groups down to one
        mostSpecificMatchedEntry = matchedRingEntries[mostSpecificMatchIndices[0]]
        
        node = mostSpecificMatchedEntry
        
        if node is None:
            raise DatabaseError('Unable to determine thermo parameters for {0}: no data for {1} or any of its ancestors.'.format(molecule, mostSpecificGroup) )

        while node is not None and node.data is None:
            # do average of its children
            success, averagedThermoData = self.__averageChildrenThermo(node)
            if success:
                node.data = averagedThermoData
            else:
                node = node.parent

        data = node.data; comment = node.label
        while isinstance(data, basestring) and data is not None:
            for entry in ring_database.entries.values():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    node = entry
                    break
        data.comment = '{0}({1})'.format(ring_database.label, comment)
        
        if thermoData is None:
            return data, node, isPartialMatch
        else:
            return addThermoData(thermoData, data, groupAdditivity=True, verbose=True), node, isPartialMatch
            # By setting verbose=True, we turn on the comments of ring correction to pass the unittest.
            # Typically this comment is very short and also very helpful to check if the ring correction is calculated correctly.

    def __averageChildrenThermo(self, node):
        """
        Use children's thermo data to guess thermo data of parent `node` 
        that doesn't have thermo data built-in in tree yet. 
        For `node` has children that have thermo data, return success flag 
        `True` and the average thermo data.
        For `node` whose children that all have no thermo data, return flag
        `False` and None for the thermo data.
        """
        if not node.children:
            if node.data is None:
                return (False, None)
            else:
                return (True, node.data)
        else:
            childrenThermoDataList = []
            for child in node.children:
                if child.data is None:
                    success, childThermoData_average = self.__averageChildrenThermo(child)
                    if success:
                        childrenThermoDataList.append(childThermoData_average)
                else:
                    childrenThermoDataList.append(child.data)
            if childrenThermoDataList:
                return (True, averageThermoData(childrenThermoDataList))
            else:
                return (False, None)

    def __addGroupThermoData(self, thermoData, database, molecule, atom):
        """
        Determine the group additivity thermodynamic data for the atom `atom`
        in the structure `structure`, and add it to the existing thermo data
        `thermoData`.
        The parameter `atom` is a dictionary of label-atom pairs like {'*',atom}
        """
        node0 = database.descendTree(molecule, atom, None)
        if node0 is None:
            raise KeyError('Node not found in thermo database for atom {0} in molecule {1}.'.format(atom, molecule))

        # It's possible (and allowed) that items in the tree may not be in the
        # library, in which case we need to fall up the tree until we find an
        # ancestor that has an entry in the library
        node = node0
        while node is not None and node.data is None:
            node = node.parent
        if node is None:
            raise DatabaseError('Unable to determine thermo parameters for {0}: no data for node {1} or any of its ancestors.'.format(molecule, node0) )

        data = node.data
        comment = node.label
        loop_count = 0
        while isinstance(data, basestring):
            loop_count += 1
            if loop_count > 100:
                raise DatabaseError("Maximum iterations reached while following thermo group data pointers. A circular"
                                    " reference may exist. Last node was {0} pointing to group called {1} in "
                                    "database {2}".format(node.label, data, database.label))

            for entry in database.entries.itervalues():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    break
            else:
                raise DatabaseError("Node {0} points to a non-existant group called {1} in database: {2}".format(node.label, data, database.label))
        data.comment = '{0}({1})'.format(database.label, comment)

        # This code prints the hierarchy of the found node; useful for debugging
#        result = ''
#        while node is not None:
#           result = ' -> ' + node.label + result
#           node = node.parent
#        print result[4:]
        
        if thermoData is None:
            return data
        else:
            return addThermoData(thermoData, data, groupAdditivity=True)

    def __removeGroupThermoData(self, thermoData, database, molecule, atom):
        """
        Based on the __addGroupThermoData method. Just replace the last line with 'return removeThermoData()'.
        Determine the group additivity thermodynamic data for the atom `atom` in the structure `structure`,
        and REMOVE it from the existing thermo data `thermoData`.
        """
        node0 = database.descendTree(molecule, atom, None)
        if node0 is None:
            raise KeyError('Node not found in database.')

        # It's possible (and allowed) that items in the tree may not be in the
        # library, in which case we need to fall up the tree until we find an
        # ancestor that has an entry in the library
        node = node0
        while node.data is None and node is not None:
            node = node.parent
        if node is None:
            raise DatabaseError('Unable to determine thermo parameters for {0}: no data for node {1} or any of'
                                ' its ancestors.'.format(molecule, node0))

        data = node.data; comment = node.label
        loop_count = 0
        while isinstance(data, basestring):
            loop_count += 1
            if loop_count > 100:
                raise DatabaseError(
                    "Maximum iterations reached while following thermo group data pointers. A circular"
                    " reference may exist. Last node was {0} pointing to group called {1} in"
                    " database {2}".format(node.label, data, database.label))
            for entry in database.entries.itervalues():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    break
            else: raise DatabaseError("Node {0} points to a non-existant group called {1} in database: {2}".format(
                node.label, data, database.label))
        data.comment = '{0}({1})'.format(database.label, comment)

        # This code prints the hierarchy of the found node; useful for debugging
#        result = ''
#        while node is not None:
#           result = ' -> ' + node.label + result
#           node = node.parent
#        print result[4:]
        
        if thermoData is None:
            return data
        else:
            return removeThermoData(thermoData, data, True)

    def getRingGroupsFromComments(self, thermoData):
        """
        Takes a string of comments from group additivity estimation, and extracts the ring and polycyclic ring groups
        from them, returning them as lists.
        """
        tokens = thermoData.comment.split()
        ringGroups = []
        polycyclicGroups = []
        regex = "\((.*)\)" #only hit outermost parentheses
        for token in tokens:
            if token.startswith('ring'):
                splitTokens = re.split(regex, token)
                assert len(splitTokens) == 3, 'token: {}'.format(token)
                groupLabel = splitTokens[1]
                ringGroups.append(self.groups['ring'].entries[groupLabel])
            if token.startswith('polycyclic'):
                splitTokens = re.split(regex, token)
                assert len(splitTokens) == 3, 'token: {}'.format(token)
                groupLabel = splitTokens[1]
                polycyclicGroups.append(self.groups['polycyclic'].entries[groupLabel])    
        
        return ringGroups, polycyclicGroups
    
    def extractSourceFromComments(self, species):
        """
        `species`: A species object containing thermo data and thermo data comments
        
        Parses the verbose string of comments from the thermo data of the species object,
        and extracts the thermo sources.

        Returns a dictionary with keys of either 'Library', 'QM', and/or 'GAV'.
        Commonly, species thermo are estimated using only one of these sources.
        However, a radical can be estimated with more than one type of source, for 
        instance a saturated library value and a GAV HBI correction, or a QM saturated value
        and a GAV HBI correction.  
        
        source = {'Library': String_Name_of_Library_Used,
                  'QM': String_of_Method_Used,
                  'GAV': Dictionary_of_Groups_Used 
                  }
                  
        The Dictionary_of_Groups_Used looks like 
        {'groupType':[List of tuples containing (Entry, Weight)]
        """
        comment = species.thermo.comment
        tokens = comment.split()
        
        
        source = {}
        
        if comment.startswith('Thermo library'):
            # Store name of the library source, which is the 3rd token in the comments
            source['Library'] = tokens[2]
            
        elif comment.startswith('QM'):
            # Store the level of the calculation, which is the 2nd token in the comments
            source['QM'] = tokens[1]
        
        
        # Check for group additivity contributions to the thermo in this species            
        
        # The contribution of the groups can be either additive or substracting
        # after changes to the polycyclic algorithm
        
        comment = comment.replace(' + ',' +')
        comment = comment.replace(' - ', ' -')
        tokens = comment.split()
        
        groups = {}
        groupTypes = self.groups.keys()
        
        
        regex = "\((.*)\)" #only hit outermost parentheses
        for token in tokens:
            weight = 1  # default contribution is additive
            if token.startswith('+'):
                token = token[1:]
            elif token.startswith('-'):
                weight = -1
                token = token[1:]
            for groupType in groupTypes:
                if token.startswith(groupType+'(') and token.endswith(')'):
                    splitTokens = re.split(regex, token)
                    groupLabel = splitTokens[1]
                    groupEntry = self.groups[groupType].entries[groupLabel]
                    # Use dictionary to combine into weights when necessary
                    if not groupType in groups:
                        groups[groupType] = {groupEntry:weight}
                    else:
                        if groupEntry in groups[groupType]:
                            groups[groupType][groupEntry] += weight
                        else:
                            groups[groupType][groupEntry] = weight
                    break
            

        if groups:
            # Indicate that group additivity is used when it is either an HBI correction
            # onto a  thermo library or QM value, or if the entire molecule is estimated using group additivity
            # Save the groups into the source dictionary
            
            # Convert groups back into tuples 
            for groupType, groupDict in groups.iteritems():
                groups[groupType] = groupDict.items()
            
            source['GAV'] = groups
            
        # Perform a sanity check that this molecule is estimated by at least one method
        if not source.keys():
            raise Exception('Species {0} thermo appears to not be estimated using any methods.'.format(species))
        
        return source

class ThermoCentralDatabaseInterface(object):
    """
    A class for interfacing with RMG online thermo central database.
    """

    def __init__(self, host, port, username, password, application):
        self.host = host
        self.port = port
        self.username = username
        self.password = password
        self.application = application
        self.client = self.connect()

    def connect(self):
        
        import pymongo

        remote_address = 'mongodb://{0}:{1}@{2}/thermoCentralDB'.format(self.username, 
                                                            self.password,
                                                            self.host)
        client = pymongo.MongoClient(remote_address, 
                                    self.port, 
                                    serverSelectionTimeoutMS=2000)
        try:
            client.server_info()
            logging.info("\nConnection success to RMG Thermo Central Database!\n")
            return client
        
        except (pymongo.errors.ServerSelectionTimeoutError,
                pymongo.errors.OperationFailure):
            logging.info("\nConnection failure to RMG Thermo Central Database...")
            logging.info("This RMG job still can run but cannot utilize data from central database.\n")
            return None

    def satisfyRegistrationRequirements(self, species, thermo, thermodb):
        """
        Given a species, check if it's allowed to register in 
        central thermo database.

        Requirements for now: 
        cyclic, 
        its thermo is estimated by GAV and no exact match/use heuristics
        """
        if not species.molecule[0].isCyclic():
            return False

        GAV_keywords = 'Thermo group additivity estimation'
        if isinstance(thermo, ThermoData) and thermo.comment.startswith(GAV_keywords):
            ringGroups, polycyclicGroups = thermodb.getRingGroupsFromComments(thermo)
            
            # use GAV generic node to estimate thermo
            for group in ringGroups + polycyclicGroups:
                if group.label in thermodb.groups['ring'].genericNodes + thermodb.groups['polycyclic'].genericNodes:
                    return True
            
            # used some heuristic way to estimate thermo
            if ") - ring(" in thermo.comment:
                return True
            else:
                return False
        else:
            return False

    def registerInCentralThermoDB(self, species):

        # choose registration table
        db =  getattr(self.client, 'thermoCentralDB')
        registration_table = getattr(db, 'registration_table')
        results_table = getattr(db, 'results_table')
        
        # prepare registration entry
        try:
            aug_inchi = species.getAugmentedInChI()

            # check if it's registered before or
            # already have available data in results_table
            registered_entries = list(registration_table.find({"aug_inchi": aug_inchi}))
            finished_entries = list(results_table.find({"aug_inchi": aug_inchi}))
            
            if len(registered_entries) + len(finished_entries) > 0 :
                return

            SMILES_input = species.molecule[0].toSMILES()
            status = 'pending'
            species_registration_entry = {'aug_inchi': aug_inchi,
                                        'SMILES_input': SMILES_input,
                                        'radical_number': species.molecule[0].getRadicalCount(),
                                        'status': status,
                                        'user': self.username,
                                        'application': self.application,
                                        'timestamp': time.time()
                                        }

            registration_table.insert(species_registration_entry)

        except ValueError:
            logging.info('Fail to generate inchi/smiles for species below:\n{0}'.format(species.toAdjacencyList()))

def findCp0andCpInf(species, heatCap):
    """
    Calculate the Cp0 and CpInf values, and add them to the HeatCapacityModel object.
    """
    if heatCap.Cp0 is None:
        Cp0 = species.calculateCp0()
        heatCap.Cp0 = (Cp0,"J/(mol*K)")
    if heatCap.CpInf is None:
        CpInf = species.calculateCpInf()  
        heatCap.CpInf = (CpInf,"J/(mol*K)")
