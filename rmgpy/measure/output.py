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

from rmgpy.chem.thermo import *
from rmgpy.chem.kinetics import *
from rmgpy.chem.states import *

from collision import *

################################################################################

def writeStates(f, states, prefix=''):
    """
    Write the :class:`StatesModel` data `states` containing molecular degree of
    freedom data to the file object `f`. The optional parameter `prefix`
    is prepended to each line of the output file, which provides an easy way to
    adjust the indentation.
    """
    f.write(prefix + 'states=States(\n')
    for mode in states.modes:
        if isinstance(mode, RigidRotor):
            f.write(prefix + '    rotationalConstants=([%s], "amu/angstrom^2"),\n' % (', '.join(['%g' % (inertia*constants.Na*1e23) for inertia in mode.inertia])))
            f.write(prefix + '    symmetry=%g,\n' % (mode.symmetry))
        elif isinstance(mode, HarmonicOscillator):
            f.write(prefix + '    frequencies=([%s], "cm^-1"),\n' % (', '.join([('%g' % freq) for freq in mode.frequencies])))
            f.write(prefix + '    frequencyScaleFactor=1.0,\n')
    if any([isinstance(mode, HinderedRotor) for mode in states.modes]):
        f.write(prefix + '    hinderedRotors=[\n')
        for mode in states.modes:
            if isinstance(mode, HinderedRotor):
                f.write(prefix + '        ((%g,"amu*angstrom^2"), (%g,"kJ/mol"), %i),\n' % (mode.inertia*constants.Na*1e23, mode.barrier/1000.0, mode.symmetry))
        f.write(prefix + '    ],\n')
    f.write(prefix + '    spinMultiplicity=%i,\n' % (states.spinMultiplicity))
    f.write(prefix + '),\n')

################################################################################

def writeNetworkSpecies(f, network):
    """
    Write all species in the given unimolecular reaction `network` to a file
    object`f`. All isomer, reactant, product, and bath gas species are
    automatically written one time each.
    """

    # Get list of all species in network
    speciesList = []
    for isomer in network.isomers:
        if isomer not in speciesList: speciesList.append(isomer)
    for reactants in network.reactants:
        for spec in reactants:
            if spec not in speciesList: speciesList.append(spec)
    for products in network.products:
        for spec in products:
            if spec not in speciesList: speciesList.append(spec)
    for spec in network.bathGas:
        if spec not in speciesList: speciesList.append(spec)

    for spec in speciesList:
        writeSpecies(f, spec)

def writeSpecies(f, spec):
    """
    Write a :class:`Species` object `spec` to a file object `f`.
    """
    f.write('species(\n')
    f.write('    label="%s",\n' % (spec))
    if len(spec.molecule) > 0:
        f.write('    SMILES="%s",\n' % (spec.molecule[0].toSMILES()))
    f.write('    E0=(%g,"kJ/mol"),\n' % (spec.E0 / 1000))
    if spec.states is not None:
        writeStates(f, spec.states, prefix='    ')
    if spec.thermo is not None:
        if isinstance(spec.thermo, ThermoData):
            f.write('    thermo=ThermoData(\n')
            f.write('        Tdata=([%s], "K"),\n' % (', '.join([('%g' % T) for T in spec.thermo.Tdata])))
            f.write('        Cpdata=([%s], "J/(mol*K)"),\n' % (', '.join([('%g' % Cp) for Cp in spec.thermo.Cpdata])))
            f.write('        H298=(%g, "kJ/mol"),\n' % (spec.thermo.H298/1000.0))
            f.write('        S298=(%g, "J/(mol*K)"),\n' % (spec.thermo.S298))
            f.write('    ),\n')
        else:
            f.write('    thermo=%r,\n' % spec.thermo)
    if spec.lennardJones is not None:
        f.write('    lennardJones=LennardJones(sigma=(%g,"m"), epsilon=(%g,"J")),\n' % (spec.lennardJones.sigma, spec.lennardJones.epsilon))
    if spec.molecularWeight != 0.0:
        f.write('    molecularWeight=(%g,"g/mol"),\n' % (spec.molecularWeight*1000))
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
    f.write('    reactants=[%s],\n' % (', '.join([('"%s"' % spec) for spec in rxn.reactants])))
    f.write('    products=[%s],\n' % (', '.join([('"%s"' % spec) for spec in rxn.products])))
    f.write('    reversible=%s,\n' % ('True' if rxn.reversible else 'False'))
    if rxn.kinetics is not None:
        if isinstance(rxn.kinetics, ArrheniusModel):
            f.write('    kinetics=Arrhenius(\n')
            if len(rxn.reactants) == 1: units = 's^-1'
            else: units = 'm^%g/(mol^%g*s)' % (3*(len(rxn.reactants)-1), len(rxn.reactants)-1)
            f.write('        A=(%g,"%s"),\n' % (rxn.kinetics.A, units))
            f.write('        n=%g,\n' % (rxn.kinetics.n))
            f.write('        Ea=(%g,"kJ/mol"),\n' % (rxn.kinetics.Ea / 1000.0))
            f.write('        T0=(%g,"K"),\n' % (rxn.kinetics.T0))
            f.write('    ),\n')
        else:
            f.write('    kinetics=%r,\n' % rxn.kinetics)
    if rxn.transitionState is not None:
        f.write('    transitionState=TransitionState(\n')
        f.write('        E0=(%g,"kJ/mol"),\n' % (rxn.transitionState.E0/1000))
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
    f.write('    reactants=[%s],\n' % (', '.join([('"%s"' % spec) for spec in rxn.reactants])))
    f.write('    products=[%s],\n' % (', '.join([('"%s"' % spec) for spec in rxn.products])))
    f.write('    reversible=%s,\n' % ('True' if rxn.reversible else 'False'))
    if rxn.kinetics is not None:
        if len(rxn.reactants) == 1: Aunits = 's^-1'
        else: Aunits = 'm^%g/(mol^%g*s)' % (3*(len(rxn.reactants)-1), len(rxn.reactants)-1)
        if isinstance(rxn.kinetics, ChebyshevModel):
            f.write('    kinetics=Chebyshev(\n')
            f.write('        Tmin=(%g,"K"),\n' % (rxn.kinetics.Tmin))
            f.write('        Tmax=(%g,"K"),\n' % (rxn.kinetics.Tmax))
            f.write('        Pmin=(%g,"bar"),\n' % (rxn.kinetics.Pmin/1e5))
            f.write('        Pmax=(%g,"bar"),\n' % (rxn.kinetics.Pmax/1e5))
            f.write('        coeffs=[\n')
            for t in range(rxn.kinetics.degreeT):
                f.write('            [%s],\n' % (', '.join([('%g' % rxn.kinetics.coeffs[t,p]) for p in range(rxn.kinetics.degreeP)])))
            f.write('        ],\n')
            f.write('    ),\n')
        elif isinstance(rxn.kinetics, PDepArrheniusModel):
            f.write('    kinetics=PDepArrhenius(\n')
            f.write('        pressures=([%s],"bar"),\n' % (', '.join([('%g' % (P / 1e5)) for P in rxn.kinetics.pressures])))
            f.write('        arrhenius=[\n')
            for arrh in rxn.kinetics.arrhenius:
                f.write('            Arrhenius(A=(%g,"%s"), n=(%g), Ea=(%g,"kJ/mol"), T0=(%g,"K")),\n' % (arrh.A, Aunits, arrh.n, arrh.Ea/1000.0, arrh.T0))
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

    logging.info('Saving output to "%s"...' % path)

    f = open(os.path.relpath(path), 'w')

    # Write each species to file
    writeNetworkSpecies(f, network)
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

    # Write each species to file
    writeNetworkSpecies(f, network)
    f.write('################################################################################\n\n')
    # Write each path reaction to file
    writeNetworkPathReactions(f, network)
    f.write('################################################################################\n\n')

    f.write('collisionModel(\n')
    if isinstance(network.collisionModel, SingleExponentialDownModel):
        f.write('    type="single exponential down",\n')
        f.write('    parameters={\n')
        f.write('        "alpha0": (%g,"kJ/mol"),\n' % (network.collisionModel.alpha0/1000.0))
        f.write('        "T0": (%g,"K"),\n' % (network.collisionModel.T0))
        f.write('        "n": %g,\n' % (network.collisionModel.n))
        f.write('    },\n')
        f.write('    bathGas={\n')
        for spec, frac in network.bathGas.iteritems():
            f.write('        "%s": %g,\n' % (spec, frac))
        f.write('    },\n')
    f.write(')\n\n')
    
    f.write('temperatures(([%s],"K"))\n' % (', '.join([('%g' % T) for T in Tlist])))
    f.write('pressures(([%s],"bar"))\n' % (', '.join([('%g' % (P/1e5)) for P in Plist])))
    dE, count = Elist
    f.write('energies(dE=(%g,"kJ/mol"), count=%i)\n' % (dE/1000.0, count))
    f.write('method("%s")\n' % (method))
    if model[0].lower() == 'chebyshev':
        f.write('interpolationModel("chebyshev", %i, %i)\n' % (model[1], model[2]))
    else:
        f.write('interpolationModel("%s")\n' % (model[0].lower()))

    f.close()
