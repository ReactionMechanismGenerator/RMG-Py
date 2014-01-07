#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   MEASURE - Master Equation Automatic Solver for Unimolecular REactions
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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
This module contains the main :meth:`execute()` function for MEASURE.
"""

import logging
import time
import os.path
import numpy

import rmgpy.constants as constants
from rmgpy.quantity import Quantity
from rmgpy.kinetics import Chebyshev, PDepArrhenius
from rmgpy.reaction import Reaction

################################################################################

class MEASURE:
    """
    A representation of a Master Equation Automatic Solver for Unimolecular
    REactions (MEASURE) job. The attributes are:
    
    =================== ======================= ================================
    Attribute           Type                    Description
    =================== ======================= ================================
    `inputFile`         ``str``                 The path to the input file
    `logFile`           ``str``                 The path to the log file
    `outputFile`        ``str``                 The path to the output file
    `drawFile`          ``str``                 The path to the PES drawing file (PNG, SVG, PDF, or PS)
    ------------------- ----------------------- --------------------------------
    `Tmin`              :class:`Quantity`       The minimum temperature at which to compute :math:`k(T,P)` values
    `Tmax`              :class:`Quantity`       The maximum temperature at which to compute :math:`k(T,P)` values
    `Tcount`            ``int``                 The number of temperatures at which to compute :math:`k(T,P)` values
    `Pmin`              :class:`Quantity`       The minimum pressure at which to compute :math:`k(T,P)` values
    `Pmax`              :class:`Quantity`       The maximum pressure at which to compute :math:`k(T,P)` values
    `Pcount`            ``int``                 The number of pressures at which to compute :math:`k(T,P)` values
    `Emin`              :class:`Quantity`       The minimum energy to use to compute :math:`k(T,P)` values
    `Emax`              :class:`Quantity`       The maximum energy to use to compute :math:`k(T,P)` values
    `grainSize`         :class:`Quantity`       The maximum energy grain size to use to compute :math:`k(T,P)` values
    `grainCount`        ``int``                 The minimum number of energy grains to use to compute :math:`k(T,P)` values
    `method`            ``str``                 The method to use to reduce the master equation to :math:`k(T,P)` values
    `model`             ``str``                 The interpolation model to fit to the computed :math:`k(T,P)` values
    ------------------- ----------------------- --------------------------------
    `network`           :class:`Network`        The unimolecular reaction network
    `Tlist`             :class:`Quantity`       An array of temperatures at which to compute :math:`k(T,P)` values
    `Plist`             :class:`Quantity`       An array of pressures at which to compute :math:`k(T,P)` values
    `Elist`             :class:`Quantity`       An array of energies to use to compute :math:`k(T,P)` values
    =================== ======================= ================================

    """
    
    def __init__(self, inputFile=None, outputFile=None, logFile=None, drawFile=None):
        self.inputFile = inputFile
        self.logFile = logFile
        self.outputFile = outputFile
        self.drawFile = drawFile
        self.clear()
    
    def clear(self):
        """
        Clear all loaded information about the job (except the file paths).
        """
        self.Tmin = None
        self.Tmax = None
        self.Tcount = None
        self.Pmin = None
        self.Pmax = None
        self.Pcount = None
        self.Emin = None
        self.Emax = None
        self.grainSize = None
        self.grainCount = None
        
        self.method = None
        self.model = None
        
        self.network = None
        self.Tlist = None
        self.Plist = None
        self.Elist = None
    
    def copy(self):
        """
        Return a copy of the current MEASURE job.
        """
        measure = MEASURE()
        
        measure.inputFile = self.inputFile
        measure.logFile = self.logFile
        measure.outputFile = self.outputFile
        measure.drawFile = self.drawFile
        
        if self.Tmin is not None: measure.Tmin = Quantity(self.Tmin)
        if self.Tmax is not None: measure.Tmax = Quantity(self.Tmax)
        measure.Tcount = self.Tcount
        if self.Pmin is not None: measure.Pmin = Quantity(self.Pmin)
        if self.Pmax is not None: measure.Pmax = Quantity(self.Pmax)
        measure.Pcount = self.Pcount
        if self.Emin is not None: measure.Emin = Quantity(self.Emin)
        if self.Emax is not None: measure.Emax = Quantity(self.Emax)
        if self.grainSize is not None: measure.grainSize = Quantity(self.grainSize)
        measure.grainCount = self.grainCount
        
        measure.method = self.method
        measure.model = self.model
        
        measure.network = self.network
        measure.Tlist = self.Tlist
        measure.Plist = self.Plist
        measure.Elist = self.Elist
        
        return measure
    
    def loadInput(self, inputFile=None):
        """
        Load a MEASURE job from the input file located at `inputFile`, or
        from the `inputFile` attribute if not given as a parameter.
        """
        from input import readFile
        
        # If an input file is specified, then it overrides the inputFile attribute
        if inputFile is not None:
            self.inputFile = inputFile
        # No matter where we got the input filename, make sure that it exists
        if not os.path.exists(self.inputFile):
            raise PDepError('Input file "{0}" does not exist.'.format(self.inputFile))
        
        # Set locations of log and output files to be in same folder as input file
        # (unless already set previously)
        inputDirectory = os.path.dirname(os.path.relpath(self.inputFile))
        if not self.outputFile:
            self.outputFile = os.path.join(inputDirectory, 'output.py')
        if not self.logFile:
            self.logFile = os.path.join(inputDirectory, 'MEASURE.log')
        
        # Load the data from the input file
        readFile(self.inputFile, self)
        
    def loadOutput(self, outputFile=None):
        """
        Load a MEASURE job from the output file located at `outputFile`, or
        from the `outputFile` attribute if not given as a parameter.
        """
        from input import readFile
        
        # If an input file is specified, then it overrides the inputFile attribute
        if outputFile is not None:
            self.outputFile = outputFile
        # No matter where we got the input filename, make sure that it exists
        if not os.path.exists(self.outputFile):
            raise PDepError('Output file "{0}" does not exist.'.format(self.outputFile))
        
        # Load the data from the output file
        readFile(self.outputFile, self)
        
    def saveInput(self, inputFile=None):
        """
        Save a MEASURE job to the output file located at `outputFile`, or
        from the `outputFile` attribute if not given as a parameter.
        """
        from output import writeFile
        
        # If an input file is specified, then it overrides the inputFile attribute
        if inputFile is not None:
            self.inputFile = inputFile
        
        writeFile(self.inputFile, self)

    def saveOutput(self, outputFile=None):
        """
        Save a MEASURE job to the output file located at `outputFile`, or
        from the `outputFile` attribute if not given as a parameter.
        """
        from output import writeFile
        
        # If an output file is specified, then it overrides the outputFile attribute
        if outputFile is not None:
            self.outputFile = outputFile
        
        writeFile(self.outputFile, self)

    def draw(self):
        """
        Draw the potential energy surface corresponding to the loaded MEASURE 
        calculation.
        """
        logging.info('Drawing potential energy surface...')
        self.network.drawPotentialEnergySurface(self.drawFile)
    
    def compute(self):
        """
        Compute the pressure-dependent rate coefficients :math:`k(T,P)` for
        the loaded MEASURE calculation.
        """
        
        # Only proceed if the input network is valid
        if self.network is None or self.network.errorString != '':
            raise PDepError('Attempted to run MEASURE calculation with invalid input.')
    
        Nisom = len(self.network.isomers)
        Nreac = len(self.network.reactants)
        Nprod = len(self.network.products)

        network = self.network   
        Tmin = self.Tmin.value_si
        Tmax = self.Tmax.value_si
        Tlist = self.Tlist.value_si
        Pmin = self.Pmin.value_si
        Pmax = self.Pmax.value_si
        Plist = self.Plist.value_si
        method = self.method
        model = self.model
        
        # Calculate the rate coefficients
        K = network.calculateRateCoefficients(Tlist, Plist, method, grainCount=self.grainCount, grainSize=self.grainSize.value_si)

        # Fit interpolation model
        from rmgpy.reaction import Reaction
        from rmgpy.measure.reaction import fitInterpolationModel
        if model[0] != '':
            logging.info('Fitting {0} interpolation models...'.format(model[0]))
        configurations = []
        configurations.extend([[isom] for isom in network.isomers])
        configurations.extend([reactants for reactants in network.reactants])
        configurations.extend([products for products in network.products])
        for i in range(Nisom+Nreac+Nprod):
            for j in range(Nisom+Nreac):
                if i != j:
                    # Check that we have nonzero k(T,P) values
                    if (numpy.any(K[:,:,i,j]) and not numpy.all(K[:,:,i,j])):
                        raise NetworkError('Zero rate coefficient encountered while updating network {0}.'.format(network))

                    # Make a new net reaction
                    netReaction = Reaction(
                        reactants=configurations[j],
                        products=configurations[i],
                        kinetics=None,
                        reversible=(i<Nisom+Nreac),
                    )
                    network.netReactions.append(netReaction)
                    
                    # Set/update the net reaction kinetics using interpolation model
                    netReaction.kinetics = fitInterpolationModel(netReaction, Tlist, Plist,
                        K[:,:,i,j],
                        model, Tmin, Tmax, Pmin, Pmax, errorCheck=True)
        logging.info('')
    
    def loadFAMEInput(self, path, moleculeDict=None):
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
        
        from network import Network
        from collision import SingleExponentialDown
        from rmgpy.species import Species, TransitionState
        from rmgpy.reaction import Reaction
        from rmgpy.transport import TransportData
        from rmgpy.statmech import HarmonicOscillator, HinderedRotor, StatesModel
        from rmgpy.thermo import ThermoData
        from rmgpy.kinetics import Arrhenius

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

        # Read method
        method = readMeaningfulLine(f).lower()
        if method == 'modifiedstrongcollision': 
            self.method = 'modified strong collision'
        elif method == 'reservoirstate': 
            self.method = 'reservoir state'

        # Read temperatures
        Tcount, Tunits, Tmin, Tmax = readMeaningfulLine(f).split()
        self.Tmin = Quantity(float(Tmin), Tunits) 
        self.Tmax = Quantity(float(Tmax), Tunits)
        self.Tcount = int(Tcount)
        Tlist = []
        for i in range(int(Tcount)):
            Tlist.append(float(readMeaningfulLine(f)))
        self.Tlist = Quantity(Tlist, Tunits)
        
        # Read pressures
        Pcount, Punits, Pmin, Pmax = readMeaningfulLine(f).split()
        self.Pmin = Quantity(float(Pmin), Punits) 
        self.Pmax = Quantity(float(Pmax), Punits)
        self.Pcount = int(Pcount)
        Plist = []
        for i in range(int(Pcount)):
            Plist.append(float(readMeaningfulLine(f)))
        self.Plist = Quantity(Plist, Punits)
        
        # Read interpolation model
        model = readMeaningfulLine(f).split()
        if model[0].lower() == 'chebyshev':
            self.model = ['chebyshev', int(model[1]), int(model[2])]
        elif model[0].lower() == 'pdeparrhenius':
            self.model = ['pdeparrhenius']
        
        # Read grain size or number of grains
        self.grainCount = 0
        self.grainSize = Quantity(0.0, "J/mol")
        for i in range(2):
            data = readMeaningfulLine(f).split()
            if data[0].lower() == 'numgrains':
                self.grainCount = int(data[1])
            elif data[0].lower() == 'grainsize':
                self.grainSize = Quantity(float(data[2]), data[1])

        # Create the Network
        self.network = Network()

        # Read collision model
        data = readMeaningfulLine(f)
        assert data.lower() == 'singleexpdown'
        alpha0units, alpha0 = readMeaningfulLine(f).split()
        T0units, T0 = readMeaningfulLine(f).split()
        n = readMeaningfulLine(f)
        collisionModel = SingleExponentialDown(
            alpha0 = Quantity(float(alpha0), alpha0units),
            T0 = Quantity(float(T0), T0units),
            n = Quantity(float(n)),
        )
        
        speciesDict = {}

        # Read bath gas parameters
        bathGas = Species(label='bath_gas', collisionModel=collisionModel)
        molWtunits, molWt = readMeaningfulLine(f).split()
        if molWtunits == 'u': molWtunits = 'g/mol'
        bathGas.molecularWeight = Quantity(float(molWt), molWtunits)
        sigmaLJunits, sigmaLJ = readMeaningfulLine(f).split()
        epsilonLJunits, epsilonLJ = readMeaningfulLine(f).split()
        bathGas.transportData = TransportData(
            sigma = Quantity(float(sigmaLJ), sigmaLJunits),
            epsilon = Quantity(float(epsilonLJ), epsilonLJunits),
        )
        self.network.bathGas = {bathGas: 1.0}
        
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
            species.transportData = TransportData(
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
            self.network.isomers.append(speciesDict[data[1]])
        for i in range(Nreac):
            data = readMeaningfulLine(f).split()
            assert data[0] == '2'
            self.network.reactants.append([speciesDict[data[1]], speciesDict[data[2]]])
        for i in range(Nprod):
            data = readMeaningfulLine(f).split()
            if data[0] == '1':
                self.network.products.append([speciesDict[data[1]]])
            elif data[0] == '2':
                self.network.products.append([speciesDict[data[1]], speciesDict[data[2]]])

        # Read path reactions
        Nrxn = int(readMeaningfulLine(f))
        for i in range(Nrxn):
            
            # Read and ignore reaction equation
            equation = readMeaningfulLine(f)
            reaction = Reaction(transitionState=TransitionState(), reversible=True)
            self.network.pathReactions.append(reaction)
            
            # Read reactant and product indices
            data = readMeaningfulLine(f).split()
            reac = int(data[0]) - 1
            prod = int(data[1]) - 1
            if reac < Nisom:
                reaction.reactants = [self.network.isomers[reac]]
            elif reac < Nisom+Nreac:
                reaction.reactants = self.network.reactants[reac-Nisom]
            else:
                reaction.reactants = self.network.products[reac-Nisom-Nreac]
            if prod < Nisom:
                reaction.products = [self.network.isomers[prod]]
            elif prod < Nisom+Nreac:
                reaction.products = self.network.reactants[prod-Nisom]
            else:
                reaction.products = self.network.products[prod-Nisom-Nreac]
            
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
    
        f.close()
    
    def loadFAMEOutput(self, path):
        """
        Load the contents of a FAME ourput file into the MEASURE object. This
        method assumes that you have already loaded the corresponding input
        file via :meth:`loadFAMEInput()`.
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
            
        with open(path, 'r') as f:
        
            method = readMeaningfulLine(f).strip()
            Tlist = numpy.array([float(d) for d in readMeaningfulLine(f).strip().split()[2:]])
            Plist = numpy.array([float(d) for d in readMeaningfulLine(f).strip().split()[2:]])
            model = readMeaningfulLine(f).strip().split()
            
            Nspec = int(readMeaningfulLine(f).strip())
            Nisom = int(readMeaningfulLine(f).strip())
            Nreac = int(readMeaningfulLine(f).strip())
            Nprod = int(readMeaningfulLine(f).strip())
            Npath = int(readMeaningfulLine(f).strip())
            Nnet = int(readMeaningfulLine(f).strip())
            
            assert Nisom == len(self.network.isomers)
            assert Nreac == len(self.network.reactants)
            assert Nprod == len(self.network.products)
            assert Npath == len(self.network.pathReactions)
            
            for n in range(Nnet):
                reac, prod = readMeaningfulLine(f).strip().split()
                reac = int(reac) - 1; prod = int(prod) - 1
                
                if reac < Nisom:
                    reactants = [self.network.isomers[reac]]
                elif reac < Nisom + Nreac:
                    reactants = self.network.reactants[reac-Nisom]
                elif reac < Nisom + Nreac + Nprod:
                    reactants = self.network.products[reac-Nisom-Nreac]
                else:
                    reactants = []
                
                if prod < Nisom:
                    products = [self.network.isomers[prod]]
                elif prod < Nisom + Nreac:
                    products = self.network.reactants[prod-Nisom]
                elif prod < Nisom + Nreac + Nprod:
                    products = self.network.products[prod-Nisom-Nreac]
                else:
                    products = []
                
                readMeaningfulLine(f)
                
                K = numpy.zeros((len(Tlist), len(Plist)), numpy.float64)
                for t in range(len(Tlist)):
                    K[t,:] = [float(d) for d in readMeaningfulLine(f).strip().split()[1:]]
                
                if len(reactants) > 1:
                    # FAME returns k(T,P) values in cm^3/mol*s and s^-1, when
                    # we want m^3/mol*s and s^-1
                    K /= 1e6
                    kunits = 'm^3/(mol*s)'
                else:
                    kunits = 's^-1'
                    
                if model[0].lower() == 'chebyshev':
                    degreeT = int(model[1]); degreeP = int(model[2])
                    coeffs = numpy.zeros((degreeT, degreeP), numpy.float64)
                    for t in range(degreeT):
                        coeffs[t,:] = [float(d) for d in readMeaningfulLine(f).strip().split()]
                    if kunits == 'm^3/(mol*s)':
                        coeffs[0,0] -= 6.0
                    kinetics = Chebyshev(coeffs=coeffs, kunits=kunits, Tmin=self.Tmin, Tmax=self.Tmax, Pmin=self.Pmin, Pmax=self.Pmax)
                elif model[0].lower() == 'pdeparrhenius':
                    pressures = []
                    arrhenius = []
                    for p in range(len(Plist)):
                        P, A, n, Ea = [float(d) for d in readMeaningfulLine(f).strip().split()]
                        if kunits == 'm^3/(mol*s)':
                            A /= 1e6
                        pressures.append(P)
                        arrhenius.append(Arrhenius(
                            A = (A,kunits), n = n, Ea = (Ea,"J/mol"), T0=(1,"K"), 
                            Tmin=self.Tmin, Tmax=self.Tmax
                        ))
                    kinetics = PDepArrhenius(pressures=(pressures,"Pa"), arrhenius=arrhenius, Tmin=self.Tmin, Tmax=self.Tmax, Pmin=self.Pmin, Pmax=self.Pmax)
                    
                netReaction = Reaction(
                    reactants = reactants,
                    products = products,
                    kinetics = kinetics
                )
                
                self.network.netReactions.append(netReaction)
        
################################################################################

def initializeLogging(level, logFile=None):
    """
    Initialize the logging system. The level of information printed is 
    determined by looking at the ``args.quiet`` and ``args.verbose`` attributes
    to see if either of the corresponding flags were set. The parameter `args`
    is an object returned by the ``argparse`` module.
    """

    # Reassign the level names so that they look better on printing
    logging.addLevelName(logging.CRITICAL, 'CRITICAL: ')
    logging.addLevelName(logging.ERROR, 'ERROR: ')
    logging.addLevelName(logging.WARNING, 'Warning: ')
    logging.addLevelName(logging.INFO, '')
    logging.addLevelName(logging.DEBUG, '')

    # Create logger
    logger = logging.getLogger()
    logger.setLevel(level)

    # Remove any old handlers that might exist
    while logger.handlers:
        logger.removeHandler(logger.handlers[0])

    # Create formatter
    formatter = logging.Formatter('%(levelname)s%(message)s')
    
    # Create console handler and set level to debug
    # Also send everything to stdout rather than stderr
    import sys
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    # create file handler
    if logFile is not None:
        fh = logging.FileHandler(filename=logFile, mode='w')
        fh.setLevel(min(logging.DEBUG,level))
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    
################################################################################

def logHeader(level=logging.INFO):
    """
    Output a header containing identifying information about RMG to the log.
    """

    logging.log(level, '###############################################################')
    logging.log(level, '# Master Equation Automatic Solver for Unimolecular REactions #')
    logging.log(level, '# (MEASURE)                                                   #')
    logging.log(level, '# Release: 0.1.0 (7 July 2010)                                #')
    logging.log(level, '# Author: Joshua W. Allen (jwallen@mit.edu)                   #')
    logging.log(level, '# Website: http://jwallen.github.com/MEASURE                  #')
    logging.log(level, '###############################################################\n')

################################################################################

def execute(inputFile, outputFile=None, drawFile=None, logFile=None, quiet=False, verbose=False):
    """
    Execute a MEASURE job using the file located at `inputFile` as the input
    file.
    """
    # We will save our output files to the directory containing the input file,
    # NOT the current working directory
    outputDirectory = os.path.dirname(os.path.relpath(inputFile))

    # Determine output level for logging system
    if quiet: 
        level = logging.WARNING
    elif verbose: 
        level = logging.DEBUG
    else:
        level = logging.INFO
        
    # Initialize the logging system
    if logFile is not None:
        logFile = os.path.abspath(logFile)
    else:
        logFile = os.path.join(outputDirectory, 'MEASURE.log')  
    initializeLogging(level, logFile)
    
    # Log start timestamp
    logging.info('MEASURE execution initiated at ' + time.asctime() + '\n')
    
    # Log header
    logHeader()
    
    # Initialize the MEASURE job
    measure = MEASURE(inputFile=inputFile, outputFile=outputFile, logFile=logFile, drawFile=drawFile)
        
    # Load input file
    measure.loadInput()
    
    # Proceed with the desired job
    if measure.network is not None and measure.network.errorString == '':
        if drawFile is not None:
            # Draw the potential energy surface
            measure.draw()
        else:
            # Compute the k(T,P) values
            measure.compute()
            # Save results to output file
            measure.saveOutput()

    # Log end timestamp
    logging.info('')
    logging.info('MEASURE execution terminated at ' + time.asctime())
