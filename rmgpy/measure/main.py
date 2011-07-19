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

from rmgpy.quantity import Quantity

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
        self.network, self.Tlist, self.Plist, self.Elist, self.method, self.model, self.Tmin, self.Tmax, self.Pmin, self.Pmax = readFile(self.inputFile)
        
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
        self.network, self.Tlist, self.Plist, self.Elist, self.method, self.model, self.Tmin, self.Tmax, self.Pmin, self.Pmax = readFile(self.outputFile)
        
    def saveInput(self, inputFile=None):
        """
        Save a MEASURE job to the output file located at `outputFile`, or
        from the `outputFile` attribute if not given as a parameter.
        """
        from output import writeFile
        
        # If an input file is specified, then it overrides the inputFile attribute
        if inputFile is not None:
            self.inputFile = inputFile
        
        writeFile(self.inputFile, self.network, self.Tlist, self.Plist, self.Elist, self.method, self.model, self.Tmin, self.Tmax, self.Pmin, self.Pmax)

    def saveOutput(self, outputFile=None):
        """
        Save a MEASURE job to the output file located at `outputFile`, or
        from the `outputFile` attribute if not given as a parameter.
        """
        from output import writeFile
        
        # If an output file is specified, then it overrides the outputFile attribute
        if outputFile is not None:
            self.outputFile = outputFile
        
        writeFile(self.outputFile, self.network, self.Tlist, self.Plist, self.Elist, self.method, self.model, self.Tmin, self.Tmax, self.Pmin, self.Pmax)

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
    
    # Only proceed if the input network is valid
    if measure.network is not None and measure.network.errorString == '':
    
        network, Tlist, Plist, Elist, method, model, Tmin, Tmax, Pmin, Pmax = measure.network, measure.Tlist, measure.Plist, measure.Elist, measure.method, measure.model, measure.Tmin, measure.Tmax, measure.Pmin, measure.Pmax

        Nisom = len(network.isomers)
        Nreac = len(network.reactants)
        Nprod = len(network.products)

        # Draw potential energy surface
        if drawFile is not None:
            logging.info('Drawing potential energy surface...')
            network.drawPotentialEnergySurface(drawFile)

        else:
            # Automatically choose a suitable set of energy grains if they were not
            # explicitly specified in the input file
            if len(Elist) == 2:
                logging.info('Automatically determining energy grains...')
                grainSize, Ngrains = Elist
                Elist = network.autoGenerateEnergyGrains(Tmax=Tmax, grainSize=grainSize, Ngrains=Ngrains)
                
            # Calculate the rate coefficients
            K, p0 = network.calculateRateCoefficients(Tlist, Plist, Elist, method)

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
            
            # Save results to file
            measure.saveOutput()

    # Log end timestamp
    logging.info('')
    logging.info('MEASURE execution terminated at ' + time.asctime())
