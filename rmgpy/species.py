#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

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

import rmgpy.quantity as quantity
from rmgpy.molecule import Molecule

#: This dictionary is used to add multiplicity to species label
_multiplicity_labels = {1:'S',2:'D',3:'T',4:'Q',5:'V',}
                           
################################################################################

class SpeciesError(Exception):
    """
    An exception class for exceptional behavior that occurs while working with
    chemical species. Pass a string describing the circumstances that caused the
    exceptional behavior.
    """
    pass

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
    `dipoleMoment`          The molecular dipole moment
    `polarizability`        The polarizability alpha
    `Zrot`                  The rotational relaxation collision number
    `energyTransferModel`   The collisional energy transfer model to use
    `reactive`              ``True`` if the species participates in reactions, ``False`` if not
    'props'                 A generic 'properties' dictionary to store user-defined flags
    ======================= ====================================================

    note::
    
        :class:`rmg.model.Species` inherits from this class, and adds some extra methods.
    """

    def __init__(self, index=-1, label='', thermo=None, conformer=None, 
                 molecule=None, transportData=None, molecularWeight=None, 
                 dipoleMoment=None, polarizability=None, Zrot=None, 
                 energyTransferModel=None, reactive=True, props=None):
        self.index = index
        self.label = label
        self.thermo = thermo
        self.conformer = conformer
        self.molecule = molecule or []
        self.transportData = transportData
        self.reactive = reactive
        self.molecularWeight = molecularWeight
        self.dipoleMoment = dipoleMoment
        self.polarizability = polarizability
        self.Zrot = Zrot
        self.energyTransferModel = energyTransferModel        
        self.props = props or {}
        
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
        if len(self.molecule) > 0: string += 'molecule=[{0!r}], '.format(self.molecule[0])
        if self.transportData is not None: string += 'transportData={0!r}, '.format(self.transportData)
        if not self.reactive: string += 'reactive={0}, '.format(self.reactive)
        if self.molecularWeight is not None: string += 'molecularWeight={0!r}, '.format(self.molecularWeight)
        if self.dipoleMoment is not None: string += 'dipoleMoment={0!r}, '.format(self.dipoleMoment)
        if self.polarizability is not None: string += 'polarizability={0!r}, '.format(self.polarizability)
        if self.Zrot is not None: string += 'Zrot={0!r}, '.format(self.Zrot)
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
        if self.index == -1: return self.label
        else: return '{0}({1:d})'.format(self.label, self.index)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Species, (self.index, self.label, self.thermo, self.conformer, self.molecule, self.transportData, self.molecularWeight, self.dipoleMoment, self.polarizability, self.Zrot, self.energyTransferModel, self.reactive, self.props))

    def getMolecularWeight(self):
        return self._molecularWeight
    def setMolecularWeight(self, value):
        self._molecularWeight = quantity.Mass(value)
    molecularWeight = property(getMolecularWeight, setMolecularWeight, """The molecular weight of the species.""")

    def getDipoleMoment(self):
        return self._dipoleMoment
    def setDipoleMoment(self, value):
        self._dipoleMoment = quantity.DipoleMoment(value)
    dipoleMoment = property(getDipoleMoment, setDipoleMoment, """The molecular dipole moment.""")

    def getPolarizability(self):
        return self._polarizability
    def setPolarizability(self, value):
        self._polarizability = quantity.Volume(value)
    polarizability = property(getPolarizability, setPolarizability, """The polarizability alpha.""")

    def getZrot(self):
        return self._Zrot
    def setZrot(self, value):
        self._Zrot = quantity.Dimensionless(value)
    Zrot = property(getZrot, setZrot, """The rotational relaxation collision number.""")

    def generateResonanceIsomers(self):
        """
        Generate all of the resonance isomers of this species. The isomers are
        stored as a list in the `molecule` attribute. If the length of
        `molecule` is already greater than one, it is assumed that all of the
        resonance isomers have already been generated.
        """
        if len(self.molecule) == 1:
            self.molecule = self.molecule[0].generateResonanceIsomers()
    
    def isIsomorphic(self, other):
        """
        Return ``True`` if the species is isomorphic to `other`, which can be
        either a :class:`Molecule` object or a :class:`Species` object.
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
        else:
            raise ValueError('Unexpected value "{0!r}" for other parameter; should be a Molecule or Species object.'.format(other))
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
            Cp = self.thermo.getHeatCapacity(T)
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
            H = self.thermo.getEnthalpy(T)
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
            S = self.thermo.getEntropy(T)
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
            G = self.thermo.getFreeEnergy(T)
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
        its resonance isomers.  This function is currently used for website purposes and testing only as it
        requires additional calculateSymmetryNumber calls.
        """
        cython.declare(symmetryNumber=cython.int)
        symmetryNumber = numpy.max([mol.getSymmetryNumber() for mol in self.molecule])
        return symmetryNumber
        
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
        if self.thermo is not None:
            Cp = self.thermo.getHeatCapacity(T)
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
        if self.thermo is not None:
            H = self.thermo.getEnthalpy(T)
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
        if self.thermo is not None:
            S = self.thermo.getEntropy(T)
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
        if self.thermo is not None:
            G = self.thermo.getFreeEnergy(T)
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
