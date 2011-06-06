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
Contains the :meth:`writeOutput()` method for saving output files to disk and
the :meth:`writeInput()` method for saving input files to disk. A number of
helper functions are used by both functions due to the high degree of similarity
between the input and output file syntax.
"""

import logging
import os.path
import time

from rmgpy.thermo import *
from rmgpy.kinetics import *
from rmgpy.statmech import *

from collision import *

################################################################################

def writeNetworkMetadata(f, network):
    """
    Write the metadata for the `network` to the given file object `f`.
    """
    if network.title != '' and network.title != 'Untitled':
        f.write('title = "{0}"\n'.format(network.title))
    if network.description != '':
        f.write('description = \\\n')
        f.write('"""\n')
        f.write('{0}\n'.format(network.description))
        f.write('"""\n')
    f.write('\n')
    
################################################################################

def writeStates(f, states, prefix=''):
    """
    Write the :class:`StatesModel` data `states` containing molecular degree of
    freedom data to the file object `f`. The optional parameter `prefix`
    is prepended to each line of the output file, which provides an easy way to
    adjust the indentation.
    """
    f.write(prefix + 'states=States(\n')
    
    # Sort modes into rotational, vibrational, and torsional categories
    rotations = []; vibrations = []; torsions = []
    for mode in states.modes:
        if isinstance(mode, RigidRotor):
            rotations.append(mode)
        elif isinstance(mode, HarmonicOscillator):
            vibrations.append(mode)
        elif isinstance(mode, HinderedRotor):
            torsions.append(mode)
    
    # Write rotational modes
    if len(rotations) == 1:
        f.write(prefix + '    rotations={0!r},\n'.format(rotations[0]))
    elif len(rotations) > 1:
        f.write(prefix + '    rotations=[\n')
        for rotation in rotations:
            f.write(prefix + '        {0!r},\n'.format(rotation))
        f.write(prefix + '    ],\n')
    
    # Write vibrational modes
    if len(vibrations) == 1:
        f.write(prefix + '    vibrations={0!r},\n'.format(vibrations[0]))
    elif len(rotations) > 1:
        f.write(prefix + '    vibrations=[\n')
        for vibration in vibrations:
            f.write(prefix + '        {0!r},\n'.format(vibration))
        f.write(prefix + '    ],\n')
    
    # Write torsional modes
    f.write(prefix + '    torsions=[\n')
    for torsion in torsions:
        f.write(prefix + '        {0!r},\n'.format(torsion))
    f.write(prefix + '    ],\n')
    
    # Write other parameters
    f.write(prefix + '    frequencyScaleFactor=1.0,\n')
    f.write(prefix + '    spinMultiplicity={0:d},\n'.format(states.spinMultiplicity))
    f.write(prefix + '),\n')

################################################################################

def writeNetworkConfigurations(f, network):
    """
    Write all configurations in the given unimolecular reaction `network` to a
    file object `f`.
    """
    # Write isomer configurations
    for isomer in network.isomers:
        f.write('isomer("{0}")\n\n'.format(isomer.label))
    # Write reactant configurations
    for reactants in network.reactants:
        f.write('reactants("{0}", "{1}")\n\n'.format(reactants[0].label, reactants[1].label))
    # No need to write product configurations, as these are assumed

################################################################################

def writeNetworkSpecies(f, network):
    """
    Write all species in the given unimolecular reaction `network` to a file
    object`f`. All isomer, reactant, product, and bath gas species are
    automatically written one time each.
    """
    for spec in network.getAllSpecies():
        writeSpecies(f, spec)

def writeSpecies(f, spec):
    """
    Write a :class:`Species` object `spec` to a file object `f`.
    """
    f.write('species(\n')
    f.write('    label="{0}",\n'.format(spec))
    if len(spec.molecule) > 0:
        f.write('    SMILES="{0}",\n'.format(spec.molecule[0].toSMILES()))
    f.write('    E0={0!r},\n'.format(spec.E0))
    if spec.states is not None:
        writeStates(f, spec.states, prefix='    ')
    if spec.thermo is not None:
        f.write('    thermo={0!r},\n'.format(spec.thermo))
    if spec.lennardJones is not None:
        f.write('    lennardJones={0!r},\n'.format(spec.lennardJones))
    if spec.molecularWeight is not None:
        f.write('    molecularWeight={0!r},\n'.format(spec.molecularWeight))
    f.write(')\n\n')

################################################################################

def writeNetworkPathReactions(f, network):
    """
    Write all path reactions in the given unimolecular reaction `network` to a
    file object`f`. The path reactions are those reactions that directly connect
    adjacent molecular configurations; these are the reactions that remain in
    the high-pressure limit.
    """
    for rxn in network.pathReactions:
        writeReaction(f, rxn)

def writeReaction(f, rxn):
    """
    Write a :class:`Reaction` object `rxn` to a file object `f`.
    """
    f.write('reaction(\n')
    f.write('    reactants=[{0}],\n'.format(', '.join([('"{0}"'.format(spec)) for spec in rxn.reactants])))
    f.write('    products=[{0}],\n'.format(', '.join([('"{0}"'.format(spec)) for spec in rxn.products])))
    f.write('    reversible={0!r},\n'.format(rxn.reversible))
    if rxn.kinetics is not None:
        f.write('    kinetics={0!r},\n'.format(rxn.kinetics))
    if rxn.transitionState is not None:
        f.write('    transitionState=TransitionState(\n')
        f.write('        E0={0!r},\n'.format(rxn.transitionState.E0))
        if rxn.transitionState.states is not None:
            writeStates(f, rxn.transitionState.states, prefix='        ')
        f.write('    ),\n')
    f.write(')\n\n')

################################################################################

def writeNetworkNetReactions(f, network):
    """
    Write all net reactions in the given unimolecular reaction `network` to a
    file object`f`. The net reactions are those reactions that can connect
    any pair of molecular configurations, not just those directly adjacent.
    These are the reactions that have pressure-dependent rate coefficients.
    """
    for rxn in network.netReactions:
        writePDepReaction(f, rxn)

def writePDepReaction(f, rxn):
    """
    Write a :class:`Reaction` object `rxn` that has pressure-dependent kinetics
    to a file object `f`.
    """

    f.write('pdepreaction(\n')
    f.write('    reactants=[{0}],\n'.format(', '.join([('"{0}"'.format(spec)) for spec in rxn.reactants])))
    f.write('    products=[{0}],\n'.format(', '.join([('"{0}"'.format(spec)) for spec in rxn.products])))
    f.write('    reversible={0!r},\n'.format(rxn.reversible))
    if rxn.kinetics is not None:
        if isinstance(rxn.kinetics, Chebyshev):
            f.write('    kinetics=Chebyshev(\n')
            f.write('        Tmin={0!r},\n'.format(rxn.kinetics.Tmin))
            f.write('        Tmax={0!r},\n'.format(rxn.kinetics.Tmax))
            f.write('        Pmin={0!r},\n'.format(rxn.kinetics.Pmin))
            f.write('        Pmax={0!r},\n'.format(rxn.kinetics.Pmax))
            f.write('        coeffs=[\n')
            for t in range(rxn.kinetics.degreeT):
                f.write('            [{0}],\n'.format(', '.join([('{0:g}'.format(rxn.kinetics.coeffs[t,p])) for p in range(rxn.kinetics.degreeP)])))
            f.write('        ],\n')
            f.write('        kunits="{0}",\n'.format(rxn.kinetics.kunits))
            f.write('    ),\n')
        elif isinstance(rxn.kinetics, PDepArrhenius):
            f.write('    kinetics=PDepArrhenius(\n')
            f.write('        pressures={0!r},\n'.format(rxn.kinetics.pressures))
            f.write('        arrhenius=[\n')
            for arrh in rxn.kinetics.arrhenius:
                f.write('            {0!r},\n'.format(arrh))
            f.write('        ],\n')
            f.write('    ),\n')

    f.write(')\n\n')

################################################################################

def writeOutput(path, network, Tlist, Plist, Elist, method, model):
    """
    Write a MEASURE output file to `path` on disk. The parameters needed mirror
    those returned by :meth:`readInput()`:

    * The :class:`Network` object `network` representing the unimolecular
      reaction network

    * The list of temperatures `Tlist` in K to be used in the master equation
      calculation

    * The list of pressures `Plist` in Pa to be used in the master equation
      calculation

    * A tuple `Elist` containing the maximum energy grain size in J/mol and the
      minimum number of energy grains to use in the master equation calculation;
      whichever of these results in more energy grains

    * The approximate `method` to use to estimate the phenomenological rate
      coefficients :math:`k(T,P)`

    * The interpolation `model` to fit the estimated :math:`k(T,P)` values to

    If successful, the file created on disk will contain all of the species
    and net reaction data, including all phenomenological rate coefficients
    :math:`k(T,P)`.
    """

    logging.info('Saving output to "{0}"...'.format(path))

    f = open(os.path.relpath(path), 'w')

    f.write('################################################################################\n')
    f.write('#\n')
    f.write('#   MEASURE output file for {0}\n'.format(network))
    f.write('#\n')
    f.write('#   Generated on {0}\n'.format(time.asctime()))
    f.write('#\n')
    f.write('################################################################################\n\n')
    
    # Write metadata
    writeNetworkMetadata(f, network)
    
    # Write each species to file
    writeNetworkSpecies(f, network)
    f.write('################################################################################\n\n')
    # Write each configuration to file
    writeNetworkConfigurations(f, network)
    f.write('################################################################################\n\n')
    
    # Write each net reaction to file
    writeNetworkNetReactions(f, network)
    
    f.close()

################################################################################

def writeInput(path, network, Tlist, Plist, Elist, method, model):
    """
    Write a MEASURE input file to `path` on disk. The parameters needed mirror
    those returned by :meth:`readInput()`:

    * The :class:`Network` object `network` representing the unimolecular
      reaction network

    * The list of temperatures `Tlist` in K to be used in the master equation
      calculation

    * The list of pressures `Plist` in Pa to be used in the master equation
      calculation

    * A tuple `Elist` containing the maximum energy grain size in J/mol and the
      minimum number of energy grains to use in the master equation calculation;
      whichever of these results in more energy grains

    * The approximate `method` to use to estimate the phenomenological rate
      coefficients :math:`k(T,P)`

    * The interpolation `model` to fit the estimated :math:`k(T,P)` values to

    If successful, the file created on disk should be able to be read in by
    :meth:`readInput()` with (hopefully) no loss of fidelity.
    """

    f = open(path, 'w')

    f.write('################################################################################\n')
    f.write('#\n')
    f.write('#   MEASURE input file for {0}\n'.format(network))
    f.write('#\n')
    f.write('#   Generated on {0}\n'.format(time.asctime()))
    f.write('#\n')
    f.write('################################################################################\n\n')
    
    # Write metadata
    writeNetworkMetadata(f, network)
    
    # Write each species to file
    writeNetworkSpecies(f, network)
    f.write('################################################################################\n\n')
    # Write each configuration to file
    writeNetworkConfigurations(f, network)
    f.write('################################################################################\n\n')
    # Write each path reaction to file
    writeNetworkPathReactions(f, network)
    f.write('################################################################################\n\n')

    f.write('collisionModel(\n')
    if isinstance(network.collisionModel, SingleExponentialDownModel):
        f.write('    type="single exponential down",\n')
        f.write('    parameters={\n')
        f.write('        "alpha0": ({0:g},"kJ/mol"),\n'.format(network.collisionModel.alpha0/1000.0))
        f.write('        "T0": ({0:g},"K"),\n'.format(network.collisionModel.T0))
        f.write('        "n": {0:g},\n'.format(network.collisionModel.n))
        f.write('    },\n')
        f.write('    bathGas={\n')
        for spec, frac in network.bathGas.iteritems():
            f.write('        "{0}": {1:g},\n'.format(spec, frac))
        f.write('    },\n')
    f.write(')\n\n')
    
    f.write('temperatures(([{0}],"K"))\n'.format(', '.join([('{0:g}'.format(T)) for T in Tlist])))
    f.write('pressures(([{0}],"bar"))\n'.format(', '.join([('{0:g}'.format(P/1e5)) for P in Plist])))
    dE, count = Elist
    f.write('energies(dE=({0:g},"kJ/mol"), count={1:d})\n'.format(dE/1000.0, count))
    f.write('method("{0}")\n'.format(method))
    if model[0].lower() == 'chebyshev':
        f.write('interpolationModel("chebyshev", {0:d}, {1:d})\n'.format(model[1], model[2]))
    else:
        f.write('interpolationModel("{0}")\n'.format(model[0].lower()))

    f.close()
