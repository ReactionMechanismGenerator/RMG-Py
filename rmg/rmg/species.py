#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
Contains classes describing chemical species.
"""

import log as logging
import pybel
import os
import xml.sax.saxutils
import quantities as pq

import settings
import structure
import thermo.model
import thermo.data
import spectral.modes
import spectral.data
import chempy.ctml_writer as ctml_writer

import chempy.species
from chempy.species import LennardJones

################################################################################

class ThermoSnapshot:
    """
    A set of thermodynamic state variables for a given state under a single
    set of conditions. During dynamic simulations, this class enables the
    simulator to update the thermodynamics of each species one per iteration.
    The attributes are:

    ============== ========================================================
    Attribute      Description
    ============== ========================================================
    `temperature`  the temperature at which the snapshot was taken in K
    `heatCapacity` the heat capacity at the current conditions in J/mol*K
    `enthalpy`     the enthalpy at the current conditions in J/mol
    `entropy`      the entropy at the current conditions in J/mol*K
    `freeEnergy`   the Gibbs free energy at the current conditions in J/mol
    ============== ========================================================

    """

    def __init__(self, temperature=0.0, heatCapacity=0.0, enthalpy=0.0, entropy=0.0):
        self.temperature = temperature
        self.heatCapacity = heatCapacity
        self.enthalpy = enthalpy
        self.entropy = entropy
        self.freeEnergy = enthalpy - temperature * entropy

    def isValid(self, temperature):
        """
        Return :data:`True` if the current snapshot is still valid at the new
        `temperature` (in K), or :data:`False` if it needs to be recalculated.
        The snapshot is considered valid if the	temperatures are equal to within
        five signficiant figures. The primary benefit is a speedup of simulations
        that are isothermal or steady-state in temperature.
        """
        if self.temperature == 0: return False
        if temperature==self.temperature: return True

        ratio =  temperature / self.temperature
        return (ratio > 0.99999 and ratio < 1.00001)


    def update(self, temperature, thermoData):
        """
        Update the thermodynamics snapshot to a new `temperature` (in K) using
        the information in `thermoData`. The latter can be any object that
        provides the necessary getHeatCapacity(), getEnthalpy(), getEntropy(),
        and getFreeEnergy() methods.
        """
        self.temperature = temperature
        self.heatCapacity = thermoData.getHeatCapacity(temperature)
        self.enthalpy = thermoData.getEnthalpy(temperature)
        self.entropy = thermoData.getEntropy(temperature)
        self.freeEnergy = self.enthalpy - self.temperature * self.entropy

################################################################################

class Species:
    """
    Represent a chemical species (including all of its resonance forms). The
    attributes are:

    ================  ==========================================================
    Attributes        Description
    ================  ==========================================================
    `index`           A unique integer identifier
    `label`           A more descriptive (but not necessarily unique) string
                      label
    `lennardJones`    The Lennard-Jones parameters for the species
    `reactive`        :data:`True` if the species is reactive, :data:`False` if
                      inert
    `spectralData`    The spectral data (degrees of freedom) for the species
    `E0`              The ground-state energy of the species in J/mol
    `structure`       A list of :class:`structure.Structure` objects
                      representing the set of resonance isomers
    `thermoData`      The thermodynamic parameters for the species, always a
                      derived class of `thermo.ThermoModel`
    `thermoSnapshot`  The thermodynamic parameters at the current point in the
                      simulation
    ================  ==========================================================

    """

    def toCantera(self):
        """Return a Cantera ctml_writer instance"""
        # contrivedly get a list like ['C', '3', 'H', '9', 'Si', '1']
        atoms = self.structure[0].toOBMol().GetSpacedFormula().split()
        # magically turn that lst into a string like 'C:3 H:9 Si:1'
        atoms = ' '.join([i+':'+j for i,j in zip(*[iter(atoms)]*2)])
        return ctml_writer.species(name = str(self),
            atoms = " %s "%atoms,
            thermo = self.thermoData.toCantera(),
            # this escaping should really be done by ctml_writer, but it doesn't do it
            note = xml.sax.saxutils.escape("%s (%s)"%(self.label,self.thermoData.comment))
               )

    def getFormula(self):
        """
        Return the chemical formula for the species.
        """
        return self.structure[0].getFormula()

    def fromAdjacencyList(self, adjstr):
        """
        Convert an adjacency list string `adjstr` to a Species object.
        """
        self.structure = [structure.Structure()]
        self.structure[0].fromAdjacencyList(adjstr)

    def fromCML(self, cmlstr):
        """
        Convert a string of CML `cmlstr` to a Species object.
        """
        self.structure = [structure.Structure()]
        self.structure[0].fromCML(cmlstr)

    def fromInChI(self, inchistr):
        """
        Convert an InChI string `inchistr` to a Species object.
        """
        self.structure = [structure.Structure()]
        self.structure[0].fromInChI(inchistr)

    def fromSMILES(self, smilesstr):
        """
        Convert a SMILES string `smilesstr` to a Species object.
        """
        self.structure = [structure.Structure()]
        self.structure[0].fromSMILES(smilesstr)

    def fromOBMol(self, obmol):
        """
        Convert an OpenBabel OBMol object `obmol` to a Species object.
        """
        self.structure = [structure.Structure()]
        self.structure[0].fromOBMol(obmol)

    def toCML(self):
        """
        Convert a Species object to CML.
        """
        return self.structure[0].toCML()

    def toInChI(self):
        """
        Convert a Species object to an InChI string.
        """
        return self.structure[0].toInChI()

    def toOBMol(self):
        """
        Convert a Species object to an OBMol object.
        """
        return self.structure[0].toOBMol()

    def toSMILES(self):
        """
        Convert a Species object to an InChI string.
        """
        return self.structure[0].toSMILES()

    def toAdjacencyList(self, strip_hydrogens=False):
        """
        Convert a Species object to an adjacency list.

        If strip_hydrogens=True then Hydrogen atoms are not reported
        (this is a valid shorthand: they will be replaced on importing such an
        adjacency list, provided that the free electron numbers are accurate)
        """
        return str(self) + '\n' + self.structure[0].toAdjacencyList(strip_hydrogens=strip_hydrogens)

    def getThermoData(self):
        """
        Get thermodynamic data for the species, generating it if necessary.

        The specific type of thermo data returned is unknown, but it will be
        a subclass of :class:`thermo.ThermoModel`.  If `self.thermoData` does
        not yet exist, it is created using :func:`generateThermoData`.
        """
        if not self.thermoData:
            self.generateThermoData()
        return self.thermoData

    def generateThermoData(self, thermoClass=thermo.model.NASAModel):
        """
        Generate thermodynamic data for the species using the thermo database.

        Generates the thermo data for each structure (resonance isomer),
        picks that with lowest H298 value, and saves it to `self.thermoData`.
        """

        thermoData = []
        for structure in self.structure:
            structure.updateAtomTypes()
            thermoData.append(thermo.data.generateThermoData(structure, thermoClass))

        # If multiple resonance isomers are present, use the thermo data of
        # the most stable isomer (i.e. one with lowest enthalpy of formation)
        # as the thermo data of the species
        self.thermoData = thermoData[0]
        lowestH298 = self.thermoData.getEnthalpy(298)
        for tdata in thermoData[1:]:
            thisH298 = tdata.getEnthalpy(298)
            if thisH298 < lowestH298:
                self.thermoData = tdata
                lowestH298 = thisH298

        # Put the most stable structure first in the list of structures
        i = thermoData.index(self.thermoData)
        s = self.structure.pop(i)
        self.structure.insert(0, s)

        return self.thermoData

    def generateSpectralData(self):
        """
        Generate spectral data for species using the frequency database. The
        spectral data comprises the molecular degrees of freedom of the
        molecule.
        """
        # Calculate the thermo data if necessary
        if not self.thermoData:
            self.generateThermoData()
        # Generate the spectral data
        return spectral.data.generateSpectralData(self.structure[0], self.thermoData)

    def getHeatCapacity(self, T):
        """
        Return the heat capacity of the species at temperature `T`.
        """
        if self.thermoSnapshot.isValid(temperature=T):
            return self.thermoSnapshot.heatCapacity
        #if self.thermoData is None: self.getThermoData()
        return self.thermoData.getHeatCapacity(T)

    def getEnthalpy(self, T):
        """
        Return the enthalpy of the species at temperature `T`.
        """
        if self.thermoSnapshot.isValid(temperature=T):
            return self.thermoSnapshot.enthalpy
        #if self.thermoData is None: self.getThermoData()
        return self.thermoData.getEnthalpy(T)

    def getEntropy(self, T):
        """
        Return the entropy of the species at temperature `T`.
        """
        if self.thermoSnapshot.isValid(temperature=T):
            return self.thermoSnapshot.entropy
        #if self.thermoData is None: self.getThermoData()
        return self.thermoData.getEntropy(T)

    def getFreeEnergy(self, T):
        """
        Return the Gibbs free energy of the species at temperature `T`.
        """
        if self.thermoSnapshot.isValid(temperature=T):
            return self.thermoSnapshot.freeEnergy
        #if self.thermoData is None: self.getThermoData()
        return self.thermoData.getFreeEnergy(T)

    def getMolecularWeight(self):
        """
        Return the molecular weight of the species in kg/mol.
        """
        # If we have structure data, use that to calculate the molecular weight
        if self.structure and len(self.structure) > 0:
            return self.structure[0].getMolecularWeight()
        # If not, try to return the 'molWt' attribute; this will raise an
        # AttributeError if it is undefined
        else:
            return self.molWt

    def isIsomorphic(self, other):
        """
        Returns :data:`True` if the two species are isomorphic and data:`False`
        otherwise.
        """
        if other.__class__ == Species:
            for struct1 in self.structure:
                for struct2 in other.structure:
                    if struct1.isIsomorphic(struct2):
                        return True
        elif other.__class__ == structure.Structure:
            for struct1 in self.structure:
                if struct1.isIsomorphic(other):
                    return True
        return False

    def isSubgraphIsomorphic(self, other):
        """
        Returns :data:`True` if the species is isomorphic with the other
        functional group and data:`False` otherwise.
        """
        for struct1 in self.structure:
            if struct1.isSubgraphIsomorphic(other):
                return True
        return False

    def findSubgraphIsomorphisms(self, other):
        """
        Returns a list of the subgraph matches between the species and the
        functional group `other`.
        """
        maps12 = []; maps21 = []
        for struct1 in self.structure:
            ismatch, map21, map12 = struct1.findSubgraphIsomorphisms(other)
            maps12.extend(map12)
            maps21.extend(map21)
        return (len(maps12) > 0), maps21, maps12

    def calculateLennardJonesParameters(self):
        """
        Calculate the Lennard-Jones collision parameters for the species. The
        parameters are not returned, but are instead stored in the
        `lennardJones` attribute.
        """
        sigma, epsilon = self.structure[0].calculateLennardJonesParameters()
        self.lennardJones = LennardJones(sigma, epsilon)

################################################################################

# The global list of species created at any point during RMG execution
# The list is stored in reverse of the order in which the species are created;
# when searching the list, it is more likely to match a recently created species
# than an older species
#: Global list of species currently in memory.
speciesList = []

# A cache of recently visited species
speciesCache = []
speciesCacheMaxSize = 4

global speciesCounter 
#: Used to label species uniquely. Incremented each time a new species is made.
speciesCounter = 0 

def checkForExistingSpecies(structure):
    """
    Check to see if an existing species contains the same
    :class:`structure.Structure` as `structure`. Returns :data:`True` or
    :data:`False`, the matched species (if found), structure (if found), and
    mapping (if found).
    """

    # First check cache and return if species is found
    for i, spec in enumerate(speciesCache):
        for struct in spec.structure:
            found, map12, map21 = structure.findIsomorphism(struct)
            if found:
                speciesCache.pop(i)
                speciesCache.insert(0, spec)
                return True, spec, struct, map21

    # Return an existing species if a match is found
    for spec in speciesList:
        for struct in spec.structure:
            found, map12, map21 = structure.findIsomorphism(struct)
            if found:
                speciesCache.pop(i)
                speciesCache.insert(0, spec)
                return True, spec, struct, map21

    # At this point we can conclude that the structure does not exist
    return False, None, None, None

def makeNewSpecies(structure, label='', reactive=True, checkExisting=True):
    """
    Attempt to make a new species based on a chemical `structure`, which is a
    :class:`Structure` object.

    The proposed species is checked against the list of existing species; if the
    species already exists, this function returns
    the existing species. If the species does not exist, a :class:`Species`
    object is created and returned after being appended to the global species
    list.
    """
    global speciesCounter
#	# Recalculate atom types for proposed structure (hopefully not necessary)
#	structure.simplifyAtomTypes()
#	structure.updateAtomTypes()

    # Check to ensure that the species is new; return the existing species if
    # not new
    if checkExisting:
        found, spec, struct, map = checkForExistingSpecies(structure)
        if found: return spec, False

    # Return None if the species has a forbidden structure
    if thermo.data.isStructureForbidden(structure): return None, False

    # Otherwise make a new species
    if label == '':
#		label = structure.getFormula()
#		for atom in structure.atoms():
#			if atom.hasFreeElectron(): label += 'J'
        label = structure.toSMILES()

    # Note in the log
    spec = Species(speciesCounter+1, label, structure, reactive)
    logging.verbose('Creating new species %s' % str(spec))
    return processNewSpecies(spec), True

def processNewSpecies(spec):
    """
    Once a species `spec` has been created (e.g. via :meth:`makeNewSpecies`),
    this function handles other aspects	of preparing it for RMG.
    """
    global speciesCounter

    speciesList.insert(0, spec)
    speciesCounter += 1
    spec.index = speciesCounter

    spec.getResonanceIsomers()
    spec.generateThermoData()

    # Generate spectral data
    if settings.spectralDataEstimation and spec.thermoData and spec.reactive:
        import spectral.data
        spec.spectralData = spectral.data.generateSpectralData(spec.structure[0], spec.thermoData)

    # Generate Lennard-Jones parameters
    spec.calculateLennardJonesParameters()

    # Draw species
    if settings.drawMolecules:
        mol = pybel.Molecule(spec.toOBMol())
        mol.draw(False, os.path.join(settings.outputDirectory, 'species/' + str(spec) + '.png'))

    # Return the newly created species
    speciesCache.insert(0, spec)
    if len(speciesCache) > speciesCacheMaxSize: speciesCache.pop()
    return spec

################################################################################
