#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Convert a FAME input file to a MEASURE input file.
"""

import argparse
import numpy
import os.path

from rmgpy.chem.molecule import Molecule
from rmgpy.chem.species import Species, TransitionState
from rmgpy.chem.reaction import Reaction
from rmgpy.chem.species import LennardJones
from rmgpy.chem.states import *
from rmgpy.chem.kinetics import ArrheniusModel

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

        # Read method
        method = readMeaningfulLine(f).lower()
        if method == 'modifiedstrongcollision': method = 'modified strong collision'
        elif method == 'reservoirstate': method = 'reservoir state'

        # Read temperatures
        data = readMeaningfulLine(f).split()
        assert data[1] == 'K'
        Tmin = float(data[2]); Tmax = float(data[3])
        Tlist = numpy.zeros(int(data[0]), numpy.float64)
        for i in range(int(data[0])):
            Tlist[i] = float(readMeaningfulLine(f))
        
        # Read pressures
        data = readMeaningfulLine(f).split()
        assert data[1] == 'Pa'
        Pmin = float(data[2]); Pmax = float(data[3])
        Plist = numpy.zeros(int(data[0]), numpy.float64)
        for i in range(int(data[0])):
            Plist[i] = float(readMeaningfulLine(f))

        # Read interpolation model
        model = readMeaningfulLine(f).split()

        # Read grain size or number of grains
        data = readMeaningfulLine(f).split()
        if data[0].lower() == 'numgrains':
            Ngrains = int(data[1])
            grainSize = 0.0
        elif data[0].lower() == 'grainsize':
            assert data[2] == 'J/mol'
            Ngrains = 0
            grainSize = float(data[2])

        network = Network()

        # Read collision model
        data = readMeaningfulLine(f).split()
        assert data[0].lower() == 'singleexpdown'
        assert data[1] == 'J/mol'
        network.collisionModel = SingleExponentialDownModel(alpha0=float(data[2]), T0=298, n=0.0)
        
        # Read bath gas parameters
        bathGas = Species()
        bathGas.molecularWeight = float(readMeaningfulLine(f).split()[1]) / 1000.0
        bathGas.lennardJones = LennardJones(sigma=float(readMeaningfulLine(f).split()[1]), epsilon=float(readMeaningfulLine(f).split()[1]))
        
        # Read species data
        Nspec = int(readMeaningfulLine(f))
        speciesDict = {}
        for i in range(Nspec):
            spec = Species()
            # Read species label
            spec.label = readMeaningfulLine(f)
            speciesDict[spec.label] = spec
            if spec.label in moleculeDict:
                spec.molecule = [moleculeDict[spec.label]]
            # Read species E0
            data = readMeaningfulLine(f).split()
            assert data[0] == 'J/mol'
            spec.E0 = float(data[1])
            # Read and ignore species thermo data
            for j in range(10):
                data = readMeaningfulLine(f)
            # Read species collision parameters
            spec.molecularWeight = float(readMeaningfulLine(f).split()[1]) / 1000.0
            spec.lennardJones = LennardJones(sigma=float(readMeaningfulLine(f).split()[1]), epsilon=float(readMeaningfulLine(f).split()[1]))
            # Read species frequencies
            spec.states = StatesModel()
            data = readMeaningfulLine(f).split()
            assert data[1] == 'cm^-1'
            frequencies = []
            for j in range(int(data[0])):
                frequencies.append(float(readMeaningfulLine(f)))
            spec.states.modes.append(HarmonicOscillator(frequencies))
            # Read species external rotors
            data = readMeaningfulLine(f).split()
            assert data[0] == '0'
            assert data[1] == 'cm^-1'
            # Read species internal rotors
            data = readMeaningfulLine(f).split()
            assert data[1] == 'cm^-1'
            frequencies = []
            for j in range(int(data[0])):
                frequencies.append(float(readMeaningfulLine(f)) * 2.9979e10)
            data = readMeaningfulLine(f).split()
            assert data[1] == 'cm^-1'
            barriers = []
            for j in range(int(data[0])):
                barriers.append(float(readMeaningfulLine(f)) * 11.96)
            inertia = [V0 / 2.0 / nu**2 / 6.022e23 for nu, V0 in zip(frequencies, barriers)]
            for I, V0 in zip(inertia, barriers):
                spec.states.modes.append(HinderedRotor(inertia=I, barrier=V0, symmetry=1))
            # Read overall symmetry number
            symm = int(readMeaningfulLine(f))
            
        # Read isomer, reactant channel, and product channel data
        Nisom = int(readMeaningfulLine(f))
        Nreac = int(readMeaningfulLine(f))
        Nprod = int(readMeaningfulLine(f))
        for i in range(Nisom):
            data = readMeaningfulLine(f).split()
            assert data[0] == '1'
            network.isomers.append(speciesDict[data[1]])
        for i in range(Nreac):
            data = readMeaningfulLine(f).split()
            assert data[0] == '2'
            network.reactants.append([speciesDict[data[1]], speciesDict[data[2]]])
        for i in range(Nprod):
            data = readMeaningfulLine(f).split()
            if data[0] == '1':
                network.products.append([speciesDict[data[1]]])
            elif data[0] == '2':
                network.products.append([speciesDict[data[1]], speciesDict[data[2]]])

        # Read path reactions
        Nrxn = int(readMeaningfulLine(f))
        for i in range(Nrxn):
            # Read and ignore reaction equation
            equation = readMeaningfulLine(f)
            rxn = Reaction(transitionState=TransitionState(), reversible=True)
            network.pathReactions.append(rxn)
            # Read reactant and product indices
            data = readMeaningfulLine(f).split()
            reac = int(data[0]) - 1
            prod = int(data[1]) - 1
            if reac < Nisom:
                rxn.reactants = [network.isomers[reac]]
            elif reac < Nisom+Nreac:
                rxn.reactants = network.reactants[reac-Nisom]
            else:
                rxn.reactants = network.products[reac-Nisom-Nreac]
            if prod < Nisom:
                rxn.products = [network.isomers[prod]]
            elif prod < Nisom+Nreac:
                rxn.products = network.reactants[prod-Nisom]
            else:
                rxn.products = network.products[prod-Nisom-Nreac]
            # Read reaction E0
            data = readMeaningfulLine(f).split()
            assert data[0] == 'J/mol'
            rxn.transitionState.E0 = float(data[1])
            # Read high-pressure limit kinetics
            data = readMeaningfulLine(f)
            assert data.lower() == 'arrhenius'
            rxn.kinetics = ArrheniusModel(
                A=float(readMeaningfulLine(f).split()[1]),
                Ea=float(readMeaningfulLine(f).split()[1]),
                n=float(readMeaningfulLine(f).split()[0]),
            )
            
        # Close file
        f.close()
        
        dirname, basename = os.path.split(os.path.abspath(fstr))
        basename, ext = os.path.splitext(basename)
        output = os.path.join(dirname, basename + '.svg')
        
        network.drawPotentialEnergySurface(output, Eunits='kcal/mol')
    