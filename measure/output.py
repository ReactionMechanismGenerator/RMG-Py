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

import logging
import os.path

from chempy.thermo import *
from chempy.kinetics import *
from chempy.states import *

################################################################################

def writeStates(f, states, prefix=''):
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
    f.write(prefix + '),\n')

################################################################################

def writeNetworkSpecies(f, network):

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
    if network.bathGas not in speciesList: speciesList.append(network.bathGas)

    for spec in speciesList:
        writeSpecies(f, spec)

def writeSpecies(f, spec):
    f.write('species(\n')
    f.write('    label="%s",\n' % (spec.label))
    if len(spec.molecule) > 0:
        f.write('    SMILES="%s",\n' % (spec.molecule[0].toSMILES()))
    f.write('    E0=(%g,"kJ/mol"),\n' % (spec.E0 / 1000))
    if spec.states is not None:
        writeStates(f, spec.states, prefix='    ')
    if spec.thermo is not None:
        if isinstance(spec.thermo, ThermoGAModel):
            f.write('    thermo=ThermoGAModel(\n')
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
    for rxn in network.pathReactions:
        writeReaction(f, rxn)

def writeReaction(f, rxn):
    f.write('reaction(\n')
    f.write('    reactants=[%s],\n' % (', '.join([('"%s"' % spec.label) for spec in rxn.reactants])))
    f.write('    products=[%s],\n' % (', '.join([('"%s"' % spec.label) for spec in rxn.products])))
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
    for rxn in network.netReactions:
        writePDepReaction(f, rxn)

def writePDepReaction(f, rxn):
    
    f.write('pdepreaction(\n')
    f.write('    reactants=[%s],\n' % (', '.join([('"%s"' % spec.label) for spec in rxn.reactants])))
    f.write('    products=[%s],\n' % (', '.join([('"%s"' % spec.label) for spec in rxn.products])))
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
    Write a MEASURE input file to `path` on disk.
    """

    from output import writeNetworkSpecies, writeNetworkPathReactions

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
        f.write('    parameters=[\n')
        f.write('        (%g,"kJ/mol"),\n' % (network.collisionModel.alpha/1000.0))
        f.write('    ],\n')
    f.write(')\n\n')

    f.write('bathGas("%s")\n\n' % (network.bathGas.label))

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
