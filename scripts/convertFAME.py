#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Convert a FAME input file to a MEASURE input file.
"""

import argparse
import logging
import numpy
import os.path

from rmgpy.molecule import Molecule
import rmgpy.constants as constants
from rmgpy.quantity import Quantity, Energy

from rmgpy.cantherm.main import CanTherm
from rmgpy.cantherm.pdep import PressureDependenceJob

from rmgpy.pdep import Network, Configuration, SingleExponentialDown
from rmgpy.species import Species, TransitionState
from rmgpy.reaction import Reaction
from rmgpy.transport import TransportData
from rmgpy.statmech import HarmonicOscillator, HinderedRotor, Conformer
from rmgpy.thermo import ThermoData
from rmgpy.kinetics import Arrhenius

################################################################################

def parseCommandLineArguments():
    """
    Parse the command-line arguments being passed to MEASURE. These are
    described in the module docstring.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('file', metavar='FILE', type=str, nargs='+',
        help='a file to convert')
    parser.add_argument('-d', '--dictionary', metavar='DICTFILE', type=str, nargs=1,
        help='the RMG dictionary corresponding to these files')
    parser.add_argument('-x', '--max-energy', metavar='VALUE UNITS', type=str, nargs=2, 
        help='A maximum energy to crop at')

    return parser.parse_args()

################################################################################
    
def loadFAMEInput(path, moleculeDict=None):
    """
    Load the contents of a FAME input file into the MEASURE object. FAME
    is an early version of MEASURE written in Fortran and used by RMG-Java.
    This script enables importing FAME input files into MEASURE so we can
    use the additional functionality that MEASURE provides. Note that it
    is mostly designed to load the FAME input files generated automatically
    by RMG-Java, and may not load hand-crafted FAME input files. If you
    specify a `moleculeDict`, then this script will use it to associate
    the species with their structures.
    """
    
    def readMeaningfulLine(f):
        line = f.readline()
        while line != '':
            line = line.strip()
            if len(line) > 0 and line[0] != '#':
                return line
            else:
                line = f.readline()
        return ''

    moleculeDict = moleculeDict or {}

    logging.info('Loading file "{0}"...'.format(path))
    f = open(path)

    job = PressureDependenceJob(network=None)
    
    # Read method
    method = readMeaningfulLine(f).lower()
    if method == 'modifiedstrongcollision': 
        job.method = 'modified strong collision'
    elif method == 'reservoirstate': 
        job.method = 'reservoir state'

    # Read temperatures
    Tcount, Tunits, Tmin, Tmax = readMeaningfulLine(f).split()
    job.Tmin = Quantity(float(Tmin), Tunits) 
    job.Tmax = Quantity(float(Tmax), Tunits)
    job.Tcount = int(Tcount)
    Tlist = []
    for i in range(int(Tcount)):
        Tlist.append(float(readMeaningfulLine(f)))
    job.Tlist = Quantity(Tlist, Tunits)
    
    # Read pressures
    Pcount, Punits, Pmin, Pmax = readMeaningfulLine(f).split()
    job.Pmin = Quantity(float(Pmin), Punits) 
    job.Pmax = Quantity(float(Pmax), Punits)
    job.Pcount = int(Pcount)
    Plist = []
    for i in range(int(Pcount)):
        Plist.append(float(readMeaningfulLine(f)))
    job.Plist = Quantity(Plist, Punits)
    
    # Read interpolation model
    model = readMeaningfulLine(f).split()
    if model[0].lower() == 'chebyshev':
        job.interpolationModel = ('chebyshev', int(model[1]), int(model[2]))
    elif model[0].lower() == 'pdeparrhenius':
        job.interpolationModel = ('pdeparrhenius',)
    
    # Read grain size or number of grains
    job.minimumGrainCount = 0
    job.maximumGrainSize = None
    for i in range(2):
        data = readMeaningfulLine(f).split()
        if data[0].lower() == 'numgrains':
            job.minimumGrainCount = int(data[1])
        elif data[0].lower() == 'grainsize':
            job.maximumGrainSize = (float(data[2]), data[1])

    # A FAME file is almost certainly created during an RMG job, so use RMG mode
    job.rmgmode = True

    # Create the Network
    job.network = Network()

    # Read collision model
    data = readMeaningfulLine(f)
    assert data.lower() == 'singleexpdown'
    alpha0units, alpha0 = readMeaningfulLine(f).split()
    T0units, T0 = readMeaningfulLine(f).split()
    n = readMeaningfulLine(f)
    energyTransferModel = SingleExponentialDown(
        alpha0 = Quantity(float(alpha0), alpha0units),
        T0 = Quantity(float(T0), T0units),
        n = float(n),
    )
    
    speciesDict = {}

    # Read bath gas parameters
    bathGas = Species(label='bath_gas', energyTransferModel=energyTransferModel)
    molWtunits, molWt = readMeaningfulLine(f).split()
    if molWtunits == 'u': molWtunits = 'amu'
    bathGas.molecularWeight = Quantity(float(molWt), molWtunits)
    sigmaLJunits, sigmaLJ = readMeaningfulLine(f).split()
    epsilonLJunits, epsilonLJ = readMeaningfulLine(f).split()
    assert epsilonLJunits == 'J'
    bathGas.transportData = TransportData(
        sigma = Quantity(float(sigmaLJ), sigmaLJunits),
        epsilon = Quantity(float(epsilonLJ) / constants.kB, 'K'),
    )
    job.network.bathGas = {bathGas: 1.0}
    
    # Read species data
    Nspec = int(readMeaningfulLine(f))
    for i in range(Nspec):
        species = Species()
        species.conformer = Conformer()
        species.energyTransferModel = energyTransferModel
        
        # Read species label
        species.label = readMeaningfulLine(f)
        speciesDict[species.label] = species
        if species.label in moleculeDict:
            species.molecule = [moleculeDict[species.label]]
        
        # Read species E0
        E0units, E0 = readMeaningfulLine(f).split()
        species.conformer.E0 = Quantity(float(E0), E0units)
        species.conformer.E0.units = 'kJ/mol'
        
        # Read species thermo data
        H298units, H298 = readMeaningfulLine(f).split()
        S298units, S298 = readMeaningfulLine(f).split()
        Cpcount, Cpunits = readMeaningfulLine(f).split()
        Cpdata = []
        for i in range(int(Cpcount)):
            Cpdata.append(float(readMeaningfulLine(f)))
        if S298units == 'J/mol*K': S298units = 'J/(mol*K)'
        if Cpunits == 'J/mol*K': Cpunits = 'J/(mol*K)'
        species.thermo = ThermoData(
            H298 = Quantity(float(H298), H298units),
            S298 = Quantity(float(S298), S298units),
            Tdata = Quantity([300,400,500,600,800,1000,1500], "K"),
            Cpdata = Quantity(Cpdata, Cpunits),
            Cp0 = (Cpdata[0], Cpunits),
            CpInf = (Cpdata[-1], Cpunits),
        )
        
        # Read species collision parameters
        molWtunits, molWt = readMeaningfulLine(f).split()
        if molWtunits == 'u': molWtunits = 'amu'
        species.molecularWeight = Quantity(float(molWt), molWtunits)
        sigmaLJunits, sigmaLJ = readMeaningfulLine(f).split()
        epsilonLJunits, epsilonLJ = readMeaningfulLine(f).split()
        assert epsilonLJunits == 'J'
        species.transportData = TransportData(
            sigma = Quantity(float(sigmaLJ), sigmaLJunits),
            epsilon = Quantity(float(epsilonLJ) / constants.kB, 'K'),
        )
        
        # Read species vibrational frequencies
        freqCount, freqUnits = readMeaningfulLine(f).split()
        frequencies = []
        for j in range(int(freqCount)):
            frequencies.append(float(readMeaningfulLine(f)))
        species.conformer.modes.append(HarmonicOscillator(
            frequencies = Quantity(frequencies, freqUnits),
        ))
        
        # Read species external rotors
        rotCount, rotUnits = readMeaningfulLine(f).split()
        if int(rotCount) > 0:
            raise NotImplementedError('Cannot handle external rotational modes in FAME input.')
        
        # Read species internal rotors
        freqCount, freqUnits = readMeaningfulLine(f).split()
        frequencies = []
        for j in range(int(freqCount)):
            frequencies.append(float(readMeaningfulLine(f)))
        barrCount, barrUnits = readMeaningfulLine(f).split()
        barriers = []
        for j in range(int(barrCount)):
            barriers.append(float(readMeaningfulLine(f)))
        if barrUnits == 'cm^-1':
            barrUnits = 'J/mol'
            barriers = [barr * constants.h * constants.c * constants.Na * 100. for barr in barriers]
        elif barrUnits in ['Hz', 's^-1']:
            barrUnits = 'J/mol'
            barriers = [barr * constants.h * constants.Na for barr in barriers]
        elif barrUnits != 'J/mol':
            raise Exception('Unexpected units "{0}" for hindered rotor barrier height.'.format(barrUnits))
        inertia = [V0 / 2.0 / (nu * constants.c * 100.)**2 / constants.Na for nu, V0 in zip(frequencies, barriers)]
        for I, V0 in zip(inertia, barriers):
            species.conformer.modes.append(HinderedRotor(
                inertia = Quantity(I,"kg*m^2"), 
                barrier = Quantity(V0,barrUnits), 
                symmetry = 1,
                semiclassical = False,
            ))
            
        # Read overall symmetry number
        species.conformer.spinMultiplicity = int(readMeaningfulLine(f))
        
    # Read isomer, reactant channel, and product channel data
    Nisom = int(readMeaningfulLine(f))
    Nreac = int(readMeaningfulLine(f))
    Nprod = int(readMeaningfulLine(f))
    for i in range(Nisom):
        data = readMeaningfulLine(f).split()
        assert data[0] == '1'
        job.network.isomers.append(speciesDict[data[1]])
    for i in range(Nreac):
        data = readMeaningfulLine(f).split()
        assert data[0] == '2'
        job.network.reactants.append([speciesDict[data[1]], speciesDict[data[2]]])
    for i in range(Nprod):
        data = readMeaningfulLine(f).split()
        if data[0] == '1':
            job.network.products.append([speciesDict[data[1]]])
        elif data[0] == '2':
            job.network.products.append([speciesDict[data[1]], speciesDict[data[2]]])

    # Read path reactions
    Nrxn = int(readMeaningfulLine(f))
    for i in range(Nrxn):
        
        # Read and ignore reaction equation
        equation = readMeaningfulLine(f)
        reaction = Reaction(transitionState=TransitionState(), reversible=True)
        job.network.pathReactions.append(reaction)
        reaction.transitionState.conformer = Conformer()
        
        # Read reactant and product indices
        data = readMeaningfulLine(f).split()
        reac = int(data[0]) - 1
        prod = int(data[1]) - 1
        if reac < Nisom:
            reaction.reactants = [job.network.isomers[reac]]
        elif reac < Nisom+Nreac:
            reaction.reactants = job.network.reactants[reac-Nisom]
        else:
            reaction.reactants = job.network.products[reac-Nisom-Nreac]
        if prod < Nisom:
            reaction.products = [job.network.isomers[prod]]
        elif prod < Nisom+Nreac:
            reaction.products = job.network.reactants[prod-Nisom]
        else:
            reaction.products = job.network.products[prod-Nisom-Nreac]
        
        # Read reaction E0
        E0units, E0 = readMeaningfulLine(f).split()
        reaction.transitionState.conformer.E0 = Quantity(float(E0), E0units)
        reaction.transitionState.conformer.E0.units = 'kJ/mol'
        
        # Read high-pressure limit kinetics
        data = readMeaningfulLine(f)
        assert data.lower() == 'arrhenius'
        Aunits, A = readMeaningfulLine(f).split()
        if '/' in Aunits:
            index = Aunits.find('/')
            Aunits = '{0}/({1})'.format(Aunits[0:index], Aunits[index+1:])
        Eaunits, Ea = readMeaningfulLine(f).split()
        n = readMeaningfulLine(f)
        reaction.kinetics = Arrhenius(
            A = Quantity(float(A), Aunits),
            Ea = Quantity(float(Ea), Eaunits),
            n = Quantity(float(n)),
        )
        reaction.kinetics.Ea.units = 'kJ/mol'

    f.close()
    
    job.network.isomers = [Configuration(isomer) for isomer in job.network.isomers]
    job.network.reactants = [Configuration(*reactants) for reactants in job.network.reactants]
    job.network.products = [Configuration(*products) for products in job.network.products]

    return job

def pruneNetwork(network, Emax):
    """
    Prune the network by removing any configurations with ground-state energy
    above `Emax` in J/mol and any reactions with transition state energy above
    `Emax` from the network. All reactions involving removed configurations
    are also removed. Any configurations that have zero reactions as a result
    of this process are also removed.
    """
    
    # Remove configurations with ground-state energies above the given Emax
    isomersToRemove = []
    for isomer in network.isomers:
        if isomer.E0 > Emax:
            isomersToRemove.append(isomer)
    for isomer in isomersToRemove:
        network.isomers.remove(isomer)
    
    reactantsToRemove = []
    for reactant in network.reactants:
        if reactant.E0 > Emax:
            reactantsToRemove.append(reactant)
    for reactant in reactantsToRemove:
        network.reactants.remove(reactant)
    
    productsToRemove = []
    for product in network.products:
        if product.E0 > Emax:
            productsToRemove.append(product)
    for product in productsToRemove:
        network.products.remove(product)
        
    # Remove path reactions involving the removed configurations
    removedConfigurations = []
    removedConfigurations.extend([isomer.species for isomer in isomersToRemove])
    removedConfigurations.extend([reactant.species for reactant in reactantsToRemove])
    removedConfigurations .extend([product.species for product in productsToRemove])
    reactionsToRemove = []
    for rxn in network.pathReactions:
        if rxn.reactants in removedConfigurations or rxn.products in removedConfigurations:
            reactionsToRemove.append(rxn)
    for rxn in reactionsToRemove:
        network.pathReactions.remove(rxn)
        
    # Remove path reactions with barrier heights above the given Emax
    reactionsToRemove = []
    for rxn in network.pathReactions:
        if rxn.transitionState.conformer.E0.value_si > Emax:
            reactionsToRemove.append(rxn)
    for rxn in reactionsToRemove:
        network.pathReactions.remove(rxn)

    def ismatch(speciesList1, speciesList2):
        if len(speciesList1) == len(speciesList2) == 1:
            return (speciesList1[0] is speciesList2[0])
        elif len(speciesList1) == len(speciesList2) == 2:
            return ((speciesList1[0] is speciesList2[0] and speciesList1[1] is speciesList2[1]) or
                    (speciesList1[0] is speciesList2[1] and speciesList1[1] is speciesList2[0]))
        elif len(speciesList1) == len(speciesList2) == 3:
            return ((speciesList1[0] is speciesList2[0] and speciesList1[1] is speciesList2[1] and speciesList1[2] is speciesList2[2]) or
                    (speciesList1[0] is speciesList2[0] and speciesList1[1] is speciesList2[2] and speciesList1[2] is speciesList2[1]) or
                    (speciesList1[0] is speciesList2[1] and speciesList1[1] is speciesList2[0] and speciesList1[2] is speciesList2[2]) or
                    (speciesList1[0] is speciesList2[1] and speciesList1[1] is speciesList2[2] and speciesList1[2] is speciesList2[0]) or
                    (speciesList1[0] is speciesList2[2] and speciesList1[1] is speciesList2[0] and speciesList1[2] is speciesList2[1]) or
                    (speciesList1[0] is speciesList2[2] and speciesList1[1] is speciesList2[1] and speciesList1[2] is speciesList2[0]))
        else:
            return False
        
    # Remove orphaned configurations (those with zero path reactions involving them)
    isomersToRemove = []
    for isomer in network.isomers:
        for rxn in network.pathReactions:
            if ismatch(rxn.reactants, isomer.species) or ismatch(rxn.products, isomer.species):
                break
        else:
            isomersToRemove.append(isomer)
    for isomer in isomersToRemove:
        network.isomers.remove(isomer)
    
    reactantsToRemove = []
    for reactant in network.reactants:
        for rxn in network.pathReactions:
            if ismatch(rxn.reactants, reactant.species) or ismatch(rxn.products, reactant.species):
                break
        else:
            reactantsToRemove.append(reactant)
    for reactant in reactantsToRemove:
        network.reactants.remove(reactant)
    
    productsToRemove = []
    for product in network.products:
        for rxn in network.pathReactions:
            if ismatch(rxn.reactants, product.species) or ismatch(rxn.products, product.species):
                break
        else:
            productsToRemove.append(product)
    for product in productsToRemove:
        network.products.remove(product)

################################################################################

if __name__ == '__main__':
    
    # Parse the command-line arguments
    args = parseCommandLineArguments()
    
    if args.max_energy:
        Emax = float(args.max_energy[0])
        Eunits = str(args.max_energy[1])
        Emax = Energy(Emax, Eunits).value_si
    else:
        Emax = None
    
    # Load RMG dictionary if specified
    moleculeDict = {}
    if args.dictionary is not None:
        f = open(args.dictionary[0])
        adjlist = ''; label = ''
        for line in f:
            if len(line.strip()) == 0:
                if len(adjlist.strip()) > 0:
                    molecule = Molecule()
                    molecule.fromAdjacencyList(adjlist)
                    moleculeDict[label] = molecule
                adjlist = ''; label = ''
            else:
                if len(adjlist.strip()) == 0:
                    label = line.strip()
                adjlist += line
                    
        f.close()
    
    method = None

    for fstr in args.file:

        # Construct CanTherm job from FAME input
        job = loadFAMEInput(fstr, moleculeDict)
        
        if Emax is not None:
            pruneNetwork(job.network, Emax)

        # Save MEASURE input file based on the above
        dirname, basename = os.path.split(os.path.abspath(fstr))
        basename, ext = os.path.splitext(basename)
        path = os.path.join(dirname, basename + '.py')
        job.saveInputFile(path)
