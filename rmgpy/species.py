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

"""
This module contains classes and functions for working with chemical species.

From the `IUPAC Compendium of Chemical Terminology 
<http://dx.doi.org/10.1351/goldbook>`_, a chemical species is "an 
ensemble of chemically identical molecular entities that can explore the same 
set of molecular energy levels on the time scale of the experiment". This
definition is purposefully vague to allow the user flexibility in application.

In RMG Py, a chemical species -- a local minimum on a potential energy surface
-- is represented in memory as a :class:`Species` object. This module also
contains the :class:`TransitionState` class for representing chemical reaction
transition states (first-order saddle points on a potential energy surface).
"""

import numpy
import cython
import logging
from operator import itemgetter

import rmgpy.quantity as quantity

from rmgpy.molecule.molecule import Atom, Bond, Molecule
from rmgpy.pdep import SingleExponentialDown
from rmgpy.statmech.conformer import Conformer
from rmgpy.thermo import Wilhoit, NASA, ThermoData
from rmgpy.exceptions import SpeciesError
#: This dictionary is used to add multiplicity to species label
_multiplicity_labels = {1:'S',2:'D',3:'T',4:'Q',5:'V',}

################################################################################

class Species(object):
    """
    A chemical species, representing a local minimum on a potential energy
    surface. The attributes are:

    ======================= ====================================================
    Attribute               Description
    ======================= ====================================================
    `index`                 A unique nonnegative integer index
    `label`                 A descriptive string label
    `thermo`                The heat capacity model for the species
    `conformer`             The molecular conformer for the species
    `molecule`              A list of the :class:`Molecule` objects describing the molecular structure
    `transportData`          A set of transport collision parameters
    `molecularWeight`       The molecular weight of the species
    `energyTransferModel`   The collisional energy transfer model to use
    `reactive`              ``True`` if the species participates in reaction families, ``False`` if not
                                Reaction libraries and seed mechanisms that include the species are
                                always considered regardless of this variable
    `props`                 A generic 'properties' dictionary to store user-defined flags
    `aug_inchi`             Unique augmented inchi
    `isSolvent`             Boolean describing whether this species is the solvent
    `creationIteration`     Iteration which the species is created within the reaction mechanism generation algorithm
    ======================= ====================================================

    note: :class:`rmg.model.Species` inherits from this class, and adds some extra methods.
    """

    # these are class level attributes?


    def __init__(self, index=-1, label='', thermo=None, conformer=None, 
                 molecule=None, transportData=None, molecularWeight=None, 
                 energyTransferModel=None, reactive=True, props=None, aug_inchi=None,
                 symmetryNumber = -1, creationIteration = 0):
        self.index = index
        self.label = label
        self.thermo = thermo
        self.conformer = conformer
        self.molecule = molecule or []
        self.transportData = transportData
        self.reactive = reactive
        self.molecularWeight = molecularWeight
        self.energyTransferModel = energyTransferModel        
        self.props = props or {}
        self.aug_inchi = aug_inchi
        self.symmetryNumber = symmetryNumber
        self.isSolvent = False
        self.creationIteration = creationIteration
        # Check multiplicity of each molecule is the same
        if molecule is not None and len(molecule)>1:
            mult = molecule[0].multiplicity
            for m in molecule[1:]:
                if mult != m.multiplicity:
                    raise SpeciesError('Multiplicities of molecules in species {species} do not match.'.format(species=label))
        


    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'Species('
        if self.index != -1: string += 'index={0:d}, '.format(self.index)
        if self.label != -1: string += 'label="{0}", '.format(self.label)
        if self.thermo is not None: string += 'thermo={0!r}, '.format(self.thermo)
        if self.conformer is not None: string += 'conformer={0!r}, '.format(self.conformer)
        if len(self.molecule) > 0: string += 'molecule={0!r}, '.format(self.molecule)
        if self.transportData is not None: string += 'transportData={0!r}, '.format(self.transportData)
        if not self.reactive: string += 'reactive={0}, '.format(self.reactive)
        if self.molecularWeight is not None: string += 'molecularWeight={0!r}, '.format(self.molecularWeight)
        if self.energyTransferModel is not None: string += 'energyTransferModel={0!r}, '.format(self.energyTransferModel)
        string = string[:-2] + ')'
        return string
    
    def _repr_png_(self):
        if len(self.molecule) > 0:
            return self.molecule[0]._repr_png_()
        else:
            return None

    def __str__(self):
        """
        Return a string representation of the species, in the form 'label(id)'.
        """
        if not self.label:
            self.label = self.molecule[0].toSMILES()
        if self.index == -1: return self.label
        else: return '{0}({1:d})'.format(self.label, self.index)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Species, (self.index, self.label, self.thermo, self.conformer, self.molecule, self.transportData, self.molecularWeight, self.energyTransferModel, self.reactive, self.props))

    def getMolecularWeight(self):
        return self._molecularWeight
    def setMolecularWeight(self, value):
        self._molecularWeight = quantity.Mass(value)
    molecularWeight = property(getMolecularWeight, setMolecularWeight, """The molecular weight of the species. (Note: value_si is in kg/molecule not kg/mole)""")

    def generate_resonance_structures(self, keep_isomorphic=True):
        """
        Generate all of the resonance structures of this species. The isomers are
        stored as a list in the `molecule` attribute. If the length of
        `molecule` is already greater than one, it is assumed that all of the
        resonance structures have already been generated.
        """
        if len(self.molecule) == 1:
            if not self.molecule[0].atomIDValid():
                self.molecule[0].assignAtomIDs()
            self.molecule = self.molecule[0].generate_resonance_structures(keep_isomorphic)
    
    def isIsomorphic(self, other, generate_res=False):
        """
        Return ``True`` if the species is isomorphic to `other`, which can be
        either a :class:`Molecule` object or a :class:`Species` object.
        If generate_res is ``True`` and other is a :class:`Species` object, the resonance structures of other will
        be generated and isomorphically compared against self. This is useful for situations where a
        "non-representative" resonance structure of self is generated, and it should be identified as the same Species,
        and be assigned a reactive=False flag.
        """
        if isinstance(other, Molecule):
            for molecule in self.molecule:
                if molecule.isIsomorphic(other):
                    return True
        elif isinstance(other, Species):
            for molecule1 in self.molecule:
                for molecule2 in other.molecule:
                    if molecule1.isIsomorphic(molecule2):
                        return True
            if generate_res:
                other_copy = other.copy(deep=True)
                other_copy.generate_resonance_structures(keepIsomorphic=False)
                for molecule1 in self.molecule:
                    for molecule2 in other_copy.molecule:
                        if molecule1.isIsomorphic(molecule2):
                            # If they are isomorphic and this was found only by generating resonance structures, append
                            # the structure in other to self.molecule as unreactive, since it is a non-representative
                            # resonance structure of it, and return `True`.
                            other_copy.molecule[0].reactive = False
                            self.molecule.append(other_copy.molecule[0])
                            return True
        else:
            raise ValueError('Unexpected value "{0!r}" for other parameter; should be a Molecule or Species object.'.format(other))
        return False

    def isIdentical(self, other):
        """
        Return ``True`` if at least one molecule of the species is identical to `other`,
        which can be either a :class:`Molecule` object or a :class:`Species` object.
        """
        if isinstance(other, Molecule):
            for molecule in self.molecule:
                if molecule.isIdentical(other):
                    return True
        elif isinstance(other, Species):
            for molecule1 in self.molecule:
                for molecule2 in other.molecule:
                    if molecule1.isIdentical(molecule2):
                        return True
        else:
            raise ValueError('Unexpected value "{0!r}" for other parameter;'
                             ' should be a Molecule or Species object.'.format(other))
        return False

    
    def fromAdjacencyList(self, adjlist):
        """
        Load the structure of a species as a :class:`Molecule` object from the
        given adjacency list `adjlist` and store it as the first entry of a 
        list in the `molecule` attribute. Does not generate resonance isomers
        of the loaded molecule.
        """
        self.molecule = [Molecule().fromAdjacencyList(adjlist)]
        # If the first line is a label, then save it to the label attribute
        for label in adjlist.splitlines():
            if label.strip():
                break
        else:
            label = ''
        if len(label.split()) > 0 and not label.split()[0].isdigit():
            self.label = label.strip()
        # Return a reference to itself so we can use e.g. Species().fromAdjacencyList()
        return self
        
    def fromSMILES(self, smiles):
        """
        Load the structure of a species as a :class:`Molecule` object from the
        given SMILES string `smiles` and store it as the first entry of a 
        list in the `molecule` attribute. Does not generate resonance isomers
        of the loaded molecule.
        """
        self.molecule = [Molecule().fromSMILES(smiles)]
        # Return a reference to itself so we can use e.g. Species().fromAdjacencyList()
        return self
    
    def toAdjacencyList(self):
        """
        Return a string containing each of the molecules' adjacency lists.
        """
        output = '\n\n'.join([m.toAdjacencyList(label=self.label, removeH=False) for m in self.molecule])
        return output
    
    def toChemkin(self):
        """
        Return the chemkin-formatted string for this species.
        """
        from rmgpy.chemkin import getSpeciesIdentifier
        return getSpeciesIdentifier(self)
        
    def toCantera(self, useChemkinIdentifier = False):
        """
        Converts the RMG Species object to a Cantera Species object
        with the appropriate thermo data.

        If useChemkinIdentifier is set to False, the species label is used
        instead. Be sure that species' labels are unique when setting it False.
        """
        import cantera as ct
        
        # Determine the number of each type of element in the molecule
        elementDict = {} # elementCounts = [0,0,0,0]
        for atom in self.molecule[0].atoms:
            # The atom itself
            symbol = atom.element.symbol
            if symbol not in elementDict:
                elementDict[symbol] = 1
            else:
                elementDict[symbol] += 1
        if useChemkinIdentifier:
            ctSpecies = ct.Species(self.toChemkin(), elementDict)
        else:
            ctSpecies = ct.Species(self.label, elementDict)
        if self.thermo:
            try:
                ctSpecies.thermo = self.thermo.toCantera()
            except Exception:
                logging.error('Could not convert thermo to create Cantera Species object. Check that thermo is a NASA polynomial.')
                raise
        
        if self.transportData:
            ctSpecies.transport = self.transportData.toCantera()
            
        return ctSpecies

    def hasStatMech(self):
        """
        Return ``True`` if the species has statistical mechanical parameters,
        or ``False`` otherwise.
        """
        return self.conformer is not None and (len(self.conformer.modes) > 0 or (len(self.molecule) > 0 and len(self.molecule[0].atoms) == 1))

    def hasThermo(self):
        """
        Return ``True`` if the species has thermodynamic parameters, or 
        ``False`` otherwise.
        """
        return self.thermo is not None

    def getPartitionFunction(self, T):
        """
        Return the partition function for the species at the specified
        temperature `T` in K.
        """
        cython.declare(Q=cython.double)
        if self.hasStatMech():
            Q = self.conformer.getPartitionFunction(T)
        else:
            raise Exception('Unable to calculate partition function for species {0!r}: no statmech data available.'.format(self.label))
        return Q
        
    def getHeatCapacity(self, T):
        """
        Return the heat capacity in J/mol*K for the species at the specified
        temperature `T` in K.
        """
        cython.declare(Cp=cython.double)
        Cp = 0.0
        if self.hasThermo():
            Cp = self.getThermoData().getHeatCapacity(T)
        elif self.hasStatMech():
            Cp = self.conformer.getHeatCapacity(T)
        else:
            raise Exception('Unable to calculate heat capacity for species {0!r}: no thermo or statmech data available.'.format(self.label))
        return Cp
    
    def getEnthalpy(self, T):
        """
        Return the enthalpy in J/mol for the species at the specified
        temperature `T` in K.
        """
        cython.declare(H=cython.double)
        H = 0.0
        if self.hasThermo():
            H = self.getThermoData().getEnthalpy(T)
        elif self.hasStatMech():
            H = self.conformer.getEnthalpy(T) + self.conformer.E0.value_si
        else:
            raise Exception('Unable to calculate enthalpy for species {0!r}: no thermo or statmech data available.'.format(self.label))
        return H
    
    def getEntropy(self, T):
        """
        Return the entropy in J/mol*K for the species at the specified
        temperature `T` in K.
        """
        cython.declare(S=cython.double)
        S = 0.0
        if self.hasThermo():
            S = self.getThermoData().getEntropy(T)
        elif self.hasStatMech():
            S = self.conformer.getEntropy(T)
        else:
            raise Exception('Unable to calculate entropy for species {0!r}: no thermo or statmech data available.'.format(self.label))
        return S

    def getFreeEnergy(self, T):
        """
        Return the Gibbs free energy in J/mol for the species at the specified
        temperature `T` in K.
        """
        cython.declare(G=cython.double)
        G = 0.0
        if self.hasThermo():
            G = self.getThermoData().getFreeEnergy(T)
        elif self.hasStatMech():
            G = self.conformer.getFreeEnergy(T) + self.conformer.E0.value_si
        else:
            raise Exception('Unable to calculate free energy for species {0!r}: no thermo or statmech data available.'.format(self.label))
        return G
        
    def getSumOfStates(self, Elist):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in J/mol.
        """
        if self.hasStatMech():
            return self.conformer.getSumOfStates(Elist)
        else:
            raise Exception('Unable to calculate sum of states for species {0!r}: no statmech data available.'.format(self.label))
        
    def getDensityOfStates(self, Elist):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` in J/mol above the ground state.
        """
        if self.hasStatMech():
            return self.conformer.getDensityOfStates(Elist)
        else:
            raise Exception('Unable to calculate density of states for species {0!r}: no statmech data available.'.format(self.label))

    def getSymmetryNumber(self):
        """
        Get the symmetry number for the species, which is the highest symmetry number amongst
        its resonance isomers and the resonance hybrid.  
        This function is currently used for website purposes and testing only as it
        requires additional calculateSymmetryNumber calls.
        """
        if self.symmetryNumber < 1:
            cython.declare(resonanceHybrid = Molecule, maxSymmetryNum = cython.short)
            resonanceHybrid = self.getResonanceHybrid()
            try:
                self.symmetryNumber = resonanceHybrid.getSymmetryNumber()
            except KeyError:
                logging.error('Wrong bond order generated by resonance hybrid.')
                logging.error('Resonance Hybrid: {}'.format(resonanceHybrid.toAdjacencyList()))
                for index, mol in enumerate(self.molecule):
                    logging.error("Resonance Structure {}: {}".format(index,mol.toAdjacencyList()))
                raise
        return self.symmetryNumber
        
    def getResonanceHybrid(self):
        """
        Returns a molecule object with bond orders that are the average 
        of all the resonance structures.
        """
        # get labeled resonance isomers
        self.generate_resonance_structures(keep_isomorphic=True)

        # return if no resonance
        if len(self.molecule) == 1:
            return self.molecule[0]

        # create a sorted list of atom objects for each resonance structure
        cython.declare(atomsFromStructures = list, oldAtoms = list, newAtoms = list,
                       numResonanceStructures=cython.short, structureNum = cython.short,
                       oldBondOrder = cython.float,
                       index1 = cython.short, index2 = cython.short,
                      newMol=Molecule, oldMol = Molecule,
                      atom1=Atom, atom2=Atom, 
                      bond=Bond,  
                      atoms=list,)

        atomsFromStructures = []
        for newMol in self.molecule:
            newMol.atoms.sort(key=lambda atom: atom.id)
            atomsFromStructures.append(newMol.atoms)

        numResonanceStructures = len(self.molecule)

        # make original structure with no bonds
        newMol = Molecule()
        originalAtoms = atomsFromStructures[0]
        for atom1 in originalAtoms:
            atom = newMol.addAtom(Atom(atom1.element))
            atom.id = atom1.id

        newAtoms = newMol.atoms

        # initialize bonds to zero order
        for index1, atom1 in enumerate(originalAtoms):
            for atom2 in atom1.bonds:
                index2 = originalAtoms.index(atom2)
                bond = Bond(newAtoms[index1],newAtoms[index2], 0)
                newMol.addBond(bond)

        # set bonds to the proper value
        for structureNum, oldMol in enumerate(self.molecule):
            oldAtoms = atomsFromStructures[structureNum]

            for index1, atom1 in enumerate(oldAtoms):
                for atom2 in atom1.bonds:
                    index2 = oldAtoms.index(atom2)

                    newBond = newMol.getBond(newAtoms[index1], newAtoms[index2])
                    oldBondOrder = oldMol.getBond(oldAtoms[index1], oldAtoms[index2]).getOrderNum()
                    newBond.applyAction(('CHANGE_BOND',None,oldBondOrder / numResonanceStructures / 2))

        newMol.updateAtomTypes(logSpecies = False, raiseException=False)
        return newMol

    def calculateCp0(self):
        """
        Return the value of the heat capacity at zero temperature in J/mol*K.
        """
        return self.molecule[0].calculateCp0()

    def calculateCpInf(self):
        """
        Return the value of the heat capacity at infinite temperature in J/mol*K.
        """
        return self.molecule[0].calculateCpInf()

    def copy(self, deep=False):
        """
        Create a copy of the current species. If the 
        kw argument 'deep' is True, then a deep copy will be made of the 
        Molecule objects in self.molecule.

        For other complex attributes, a deep copy will always be made.
        """
        from copy import deepcopy

        cython.declare(other=Species)

        other = Species.__new__(Species)

        other.index = self.index

        other.label = self.label

        other.thermo = deepcopy(self.thermo)

        other.molecule = []
        for mol in self.molecule:
            other.molecule.append(mol.copy(deep=deep))

        other.conformer = deepcopy(self.conformer)

        other.transportData = deepcopy(self.transportData)

        other.molecularWeight = deepcopy(self.molecularWeight)
        other.energyTransferModel = deepcopy(self.energyTransferModel)
        other.reactive = self.reactive        
        other.props = deepcopy(self.props)

        return other

    def getAugmentedInChI(self):
        if self.aug_inchi is None:
            self.aug_inchi = self.generate_aug_inchi()
        return self.aug_inchi

    def generate_aug_inchi(self):
        candidates = []
        self.generate_resonance_structures()
        for mol in self.molecule:
            try:
                cand = [mol.toAugmentedInChI(),mol]
            except ValueError:
                pass  # not all resonance structures can be parsed into InChI (e.g. if containing a hypervalance atom)
            else:
                candidates.append(cand)
        candidates = sorted(candidates, key=itemgetter(0))
        for cand in candidates:
            if all(atom.charge == 0 for atom in cand[1].vertices):
                return cand[0]
        return candidates[0][0]

    def getThermoData(self, solventName = ''):
        """
        Returns a `thermoData` object of the current Species object.

        If the thermo object already exists, it is either of the (Wilhoit, ThermoData)
        type, or it is a Future.

        If the type of the thermo attribute is Wilhoit, or ThermoData,
        then it is converted into a NASA format.

        If it is a Future, then a blocking call is made to retrieve the NASA object.
        If the thermo object did not exist yet, the thermo object is generated.        
        """

        from rmgpy.thermo.thermoengine import submit
        
        if self.thermo:
            if not isinstance(self.thermo, (NASA, Wilhoit, ThermoData)):
                self.thermo = self.thermo.result()
        else:
            submit(self, solventName)
            if not isinstance(self.thermo, (NASA, Wilhoit, ThermoData)):
                self.thermo = self.thermo.result()

        return self.thermo       
            
    def generateTransportData(self):
        """
        Generate the transportData parameters for the species.
        """
        from rmgpy.data.rmg import getDB
        try:
            transportDB = getDB('transport')        
            if not transportDB: raise Exception
        except Exception, e:
            logging.debug('Could not obtain the transport database. Not generating transport...')
            raise e

        #count = sum([1 for atom in self.molecule[0].vertices if atom.isNonHydrogen()])
        self.transportData = transportDB.getTransportProperties(self)[0]


    def getTransportData(self):
        """
        Returns the transport data associated with this species, and
        calculates it if it is not yet available.
        """

        if not self.transportData:
            self.generateTransportData()

        return self.transportData

    def generateStatMech(self):
        """
        Generate molecular degree of freedom data for the species. You must
        have already provided a thermodynamics model using e.g.
        :meth:`generateThermoData()`.
        """

        from rmgpy.data.rmg import getDB
        try:
            statmechDB = getDB('statmech')        
            if not statmechDB: raise Exception
        except Exception, e:
            logging.debug('Could not obtain the stat. mech database. Not generating stat. mech...')
            raise e

        molecule = self.molecule[0]
        conformer = statmechDB.getStatmechData(molecule, self.getThermoData())

        if self.conformer is None:
            self.conformer = Conformer()
        self.conformer.E0 = self.getThermoData().E0
        self.conformer.modes = conformer.modes
        self.conformer.spinMultiplicity = conformer.spinMultiplicity
        
    def generateEnergyTransferModel(self):
        """
        Generate the collisional energy transfer model parameters for the
        species. This "algorithm" is *very* much in need of improvement.
        """
        self.energyTransferModel = SingleExponentialDown(
            alpha0 = (300*0.011962,"kJ/mol"),
            T0 = (300,"K"),
            n = 0.85,
        ) 
################################################################################

class TransitionState():
    """
    A chemical transition state, representing a first-order saddle point on a
    potential energy surface. The attributes are:

    =============== ============================================================
    Attribute       TDescription
    =============== ============================================================
    `label`         A descriptive string label
    `conformer`     The molecular degrees of freedom model for the species
    `frequency`     The negative frequency of the first-order saddle point
    `tunneling`     The type of tunneling model to use for tunneling through the reaction barrier
    `degeneracy`    The reaction path degeneracy
    =============== ============================================================
    """

    def __init__(self, label='', conformer=None, frequency=None, tunneling=None, degeneracy=1):
        self.label = label
        self.conformer = conformer
        self.frequency = frequency
        self.tunneling = tunneling
        self.degeneracy = degeneracy

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'TransitionState('
        if self.label != '': string += 'label="{0}", '.format(self.label)
        if self.conformer is not None: string += 'conformer={0!r}, '.format(self.conformer)
        if self.frequency is not None: string += 'frequency={0!r}, '.format(self.frequency)
        if self.tunneling is not None: string += 'tunneling={0!r}, '.format(self.tunneling)
        if self.degeneracy != 1: string += 'degeneracy={0}, '.format(self.degeneracy)
        string = string[:-2] + ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (TransitionState, (self.label, self.conformer, self.frequency, self.tunneling, self.degeneracy))

    def getFrequency(self):
        return self._frequency
    def setFrequency(self, value):
        self._frequency = quantity.Frequency(value)
    frequency = property(getFrequency, setFrequency, """The negative frequency of the first-order saddle point.""")

    def getPartitionFunction(self, T):
        """
        Return the partition function for the transition state at the
        specified temperature `T` in K.
        """
        cython.declare(Q=cython.double)
        if self.conformer is not None and len(self.conformer.modes) > 0:
            Q = self.conformer.getPartitionFunction(T)
        else:
            raise Exception('Unable to calculate partition function for transition state {0!r}: no statmech data available.'.format(self.label))
        return Q
        
    def getHeatCapacity(self, T):
        """
        Return the heat capacity in J/mol*K for the transition state at the
        specified temperature `T` in K.
        """
        cython.declare(Cp=cython.double)
        Cp = 0.0

        if self.getThermoData() is not None:
            Cp = self.getThermoData().getHeatCapacity(T)
        elif self.conformer is not None and len(self.conformer.modes) > 0:
            Cp = self.conformer.getHeatCapacity(T)
        else:
            raise Exception('Unable to calculate heat capacity for transition state {0!r}: no thermo or statmech data available.'.format(self.label))
        return Cp
    
    def getEnthalpy(self, T):
        """
        Return the enthalpy in J/mol for the transition state at the
        specified temperature `T` in K.
        """
        cython.declare(H=cython.double)
        H = 0.0

        if self.getThermoData() is not None:
            H = self.getThermoData().getEnthalpy(T)
        elif self.conformer is not None and len(self.conformer.modes) > 0:
            H = self.conformer.getEnthalpy(T)
        else:
            raise Exception('Unable to calculate enthalpy for transition state {0!r}: no thermo or statmech data available.'.format(self.label))
        return H
    
    def getEntropy(self, T):
        """
        Return the entropy in J/mol*K for the transition state at the
        specified temperature `T` in K.
        """
        cython.declare(S=cython.double)
        S = 0.0

        if self.getThermoData() is not None:
            S = self.getThermoData().getEntropy(T)
        elif self.conformer is not None and len(self.conformer.modes) > 0:
            S = self.conformer.getEntropy(T)
        else:
            raise Exception('Unable to calculate entropy for transition state {0!r}: no thermo or statmech data available.'.format(self.label))
        return S

    def getFreeEnergy(self, T):
        """
        Return the Gibbs free energy in J/mol for the transition state at the
        specified temperature `T` in K.
        """
        cython.declare(G=cython.double)
        G = 0.0

        if self.getThermoData() is not None:
            G = self.getThermoData().getFreeEnergy(T)
        elif self.conformer is not None and len(self.conformer.modes) > 0:
            G = self.conformer.getFreeEnergy(T)
        else:
            raise Exception('Unable to calculate free energy for transition state {0!r}: no thermo or statmech data available.'.format(self.label))
        return G
        
    def getSumOfStates(self, Elist):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in J/mol.
        """
        if self.conformer is not None and len(self.conformer.modes) > 0:
            return self.conformer.getSumOfStates(Elist)
        else:
            raise Exception('Unable to calculate sum of states for transition state {0!r}: no statmech data available.'.format(self.label))
        
    def getDensityOfStates(self, Elist):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` in J/mol above the ground state.
        """
        if self.conformer is not None and len(self.conformer.modes) > 0:
            return self.conformer.getDensityOfStates(Elist)
        else:
            raise Exception('Unable to calculate density of states for transition state {0!r}: no statmech data available.'.format(self.label))

    def calculateTunnelingFactor(self, T):
        """
        Calculate and return the value of the canonical tunneling correction 
        factor for the reaction at the given temperature `T` in K.
        """
        if self.tunneling is not None:
            return self.tunneling.calculateTunnelingFactor(T)
        else:
            # Return unity
            return 1.0
    
    def calculateTunnelingFunction(self, Elist):
        """
        Calculate and return the value of the microcanonical tunneling 
        correction for the reaction at the given energies `Elist` in J/mol.
        """
        if self.tunneling is not None:
            return self.tunneling.calculateTunnelingFunction(Elist)
        else:
            # Return step function
            kappa = numpy.ones_like(Elist)
            E0 = float(self.conformer.E0.value_si)
            for r in range(Elist.shape[0]):
                if Elist[r] >= E0:
                    break
                kappa[r] = 0.0
            return kappa
