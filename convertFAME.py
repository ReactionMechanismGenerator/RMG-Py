#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Convert a FAME input file to a MEASURE input file.
"""

import argparse
import numpy
import os.path

from rmgpy.quantity import Quantity
from rmgpy.molecule import Molecule
from rmgpy.species import Species, TransitionState
from rmgpy.reaction import Reaction
from rmgpy.species import LennardJones
from rmgpy.statmech import *
from rmgpy.thermo import ThermoData
from rmgpy.kinetics import Arrhenius

from rmgpy.measure.main import MEASURE
from rmgpy.measure.network import Network
from rmgpy.measure.collision import SingleExponentialDownModel

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

    return parser.parse_args()

################################################################################

def readMeaningfulLine(f):

    line = f.readline()
    while line != '':
        line = line.strip()
        if len(line) > 0 and line[0] != '#':
            return line
        else:
            line = f.readline()
    return ''

if __name__ == '__main__':
    
    # Parse the command-line arguments
    args = parseCommandLineArguments()
    
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

        print 'Loading file "%s"...' % fstr
        
        f = open(fstr)

        # Construct MEASURE job
        measure = MEASURE()
        
        # Read method
        method = readMeaningfulLine(f).lower()
        if method == 'modifiedstrongcollision': 
            measure.method = 'modified strong collision'
        elif method == 'reservoirstate': 
            measure.method = 'reservoir state'

        # Read temperatures
        Tcount, Tunits, Tmin, Tmax = readMeaningfulLine(f).split()
        measure.Tmin = Quantity(float(Tmin), Tunits) 
        measure.Tmax = Quantity(float(Tmax), Tunits)
        measure.Tcount = int(Tcount)
        Tlist = []
        for i in range(int(Tcount)):
            Tlist.append(float(readMeaningfulLine(f)))
        measure.Tlist = Quantity(Tlist, Tunits)
        
        # Read pressures
        Pcount, Punits, Pmin, Pmax = readMeaningfulLine(f).split()
        measure.Pmin = Quantity(float(Pmin), Punits) 
        measure.Pmax = Quantity(float(Pmax), Punits)
        measure.Pcount = int(Pcount)
        Plist = []
        for i in range(int(Pcount)):
            Plist.append(float(readMeaningfulLine(f)))
        measure.Plist = Quantity(Plist, Punits)
        
        # Read interpolation model
        model = readMeaningfulLine(f).split()
        if model[0].lower() == 'chebyshev':
            measure.model = ['chebyshev', int(model[1]), int(model[2])]
        elif model[0].lower() == 'pdeparrhenius':
            measure.model = ['pdeparrhenius']
        
        # Read grain size or number of grains
        data = readMeaningfulLine(f).split()
        if data[0].lower() == 'numgrains':
            measure.grainCount = int(data[1])
            measure.grainSize = Quantity(0.0, "J/mol")
        elif data[0].lower() == 'grainsize':
            measure.grainCount = 0
            measure.grainSize = Quantity(float(data[2]), data[1])

        # Create the Network
        measure.network = Network()

        # Read collision model
        data = readMeaningfulLine(f)
        assert data.lower() == 'singleexpdown'
        alpha0units, alpha0 = readMeaningfulLine(f).split()
        T0units, T0 = readMeaningfulLine(f).split()
        n = readMeaningfulLine(f)
        measure.network.collisionModel = SingleExponentialDownModel(
            alpha0 = Quantity(float(alpha0), alpha0units),
            T0 = Quantity(float(T0), T0units),
            n = Quantity(float(n)),
        )
        
        speciesDict = {}

        # Read bath gas parameters
        bathGas = Species(label='bath_gas')
        molWtunits, molWt = readMeaningfulLine(f).split()
        if molWtunits == 'u': molWtunits = 'g/mol'
        bathGas.molecularWeight = Quantity(float(molWt), molWtunits)
        sigmaLJunits, sigmaLJ = readMeaningfulLine(f).split()
        epsilonLJunits, epsilonLJ = readMeaningfulLine(f).split()
        bathGas.lennardJones = LennardJones(
            sigma = Quantity(float(sigmaLJ), sigmaLJunits),
            epsilon = Quantity(float(epsilonLJ), epsilonLJunits),
        )
        measure.network.bathGas = {bathGas: 1.0}
        
        # Read species data
        Nspec = int(readMeaningfulLine(f))
        for i in range(Nspec):
            species = Species()
            
            # Read species label
            species.label = readMeaningfulLine(f)
            speciesDict[species.label] = species
            if species.label in moleculeDict:
                species.molecule = [moleculeDict[species.label]]
            
            # Read species E0
            E0units, E0 = readMeaningfulLine(f).split()
            species.E0 = Quantity(float(E0), E0units)
            
            # Read species thermo data
            H298units, H298 = readMeaningfulLine(f).split()
            S298units, S298 = readMeaningfulLine(f).split()
            Cpcount, Cpunits = readMeaningfulLine(f).split()
            Cpdata = []
            for i in range(int(Cpcount)):
                Cpdata.append(float(readMeaningfulLine(f)))
            species.thermo = ThermoData(
                H298 = Quantity(float(H298), H298units),
                S298 = Quantity(float(S298), S298units),
                Tdata = Quantity([300,400,500,600,800,1000,1500], "K"),
                Cpdata = Quantity(Cpdata, Cpunits),
            )
            
            # Read species collision parameters
            molWtunits, molWt = readMeaningfulLine(f).split()
            if molWtunits == 'u': molWtunits = 'g/mol'
            species.molecularWeight = Quantity(float(molWt), molWtunits)
            sigmaLJunits, sigmaLJ = readMeaningfulLine(f).split()
            epsilonLJunits, epsilonLJ = readMeaningfulLine(f).split()
            species.lennardJones = LennardJones(
                sigma = Quantity(float(sigmaLJ), sigmaLJunits),
                epsilon = Quantity(float(epsilonLJ), epsilonLJunits),
            )
            
            species.states = StatesModel()
            
            # Read species vibrational frequencies
            freqCount, freqUnits = readMeaningfulLine(f).split()
            frequencies = []
            for j in range(int(freqCount)):
                frequencies.append(float(readMeaningfulLine(f)))
            species.states.modes.append(HarmonicOscillator(
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
                species.states.modes.append(HinderedRotor(
                    inertia = Quantity(I,"kg*m^2"), 
                    barrier = Quantity(V0,barrUnits), 
                    symmetry = 1,
                ))
                
            # Read overall symmetry number
            species.states.spinMultiplicity = int(readMeaningfulLine(f))
            
        # Read isomer, reactant channel, and product channel data
        Nisom = int(readMeaningfulLine(f))
        Nreac = int(readMeaningfulLine(f))
        Nprod = int(readMeaningfulLine(f))
        for i in range(Nisom):
            data = readMeaningfulLine(f).split()
            assert data[0] == '1'
            measure.network.isomers.append(speciesDict[data[1]])
        for i in range(Nreac):
            data = readMeaningfulLine(f).split()
            assert data[0] == '2'
            measure.network.reactants.append([speciesDict[data[1]], speciesDict[data[2]]])
        for i in range(Nprod):
            data = readMeaningfulLine(f).split()
            if data[0] == '1':
                measure.network.products.append([speciesDict[data[1]]])
            elif data[0] == '2':
                measure.network.products.append([speciesDict[data[1]], speciesDict[data[2]]])

        # Read path reactions
        Nrxn = int(readMeaningfulLine(f))
        for i in range(Nrxn):
            
            # Read and ignore reaction equation
            equation = readMeaningfulLine(f)
            reaction = Reaction(transitionState=TransitionState(), reversible=True)
            measure.network.pathReactions.append(reaction)
            
            # Read reactant and product indices
            data = readMeaningfulLine(f).split()
            reac = int(data[0]) - 1
            prod = int(data[1]) - 1
            if reac < Nisom:
                reaction.reactants = [measure.network.isomers[reac]]
            elif reac < Nisom+Nreac:
                reaction.reactants = measure.network.reactants[reac-Nisom]
            else:
                reaction.reactants = measure.network.products[reac-Nisom-Nreac]
            if prod < Nisom:
                reaction.products = [measure.network.isomers[prod]]
            elif prod < Nisom+Nreac:
                reaction.products = measure.network.reactants[prod-Nisom]
            else:
                reaction.products = measure.network.products[prod-Nisom-Nreac]
            
            # Read reaction E0
            E0units, E0 = readMeaningfulLine(f).split()
            reaction.transitionState.E0 = Quantity(float(E0), E0units)
            
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
            
        # Close file
        f.close()

        # Save MEASURE input file based on the above
        dirname, basename = os.path.split(os.path.abspath(fstr))
        basename, ext = os.path.splitext(basename)
        path = os.path.join(dirname, basename + '.py')
        measure.saveInput(path)
