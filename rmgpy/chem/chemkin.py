#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
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
This module contains functions for writing of Chemkin input files.
"""

import re

from thermo import MultiNASA
from kinetics import *

################################################################################

class ChemkinError(Exception):
    """
    An exception class for exceptional behavior involving Chemkin files. Pass a
    string describing the circumstances that caused the exceptional behavior.
    """
    pass

################################################################################

def getSpeciesIdentifier(species):
    """
    Return a string identifier for the provided `species` that can be used in a
    Chemkin file. Although the Chemkin format allows up to 16 characters for a
    species identifier, this function uses a maximum of 10 to ensure that all
    reaction equations fit in the maximum limit of 52 characters.
    """

    # First try to use the label and index
    # The label can only contain alphanumeric characters, hyphens, and underscores
    if len(species.label) > 0 and species.index >= 0 and not re.search('[^A-Za-z0-9\-_]+', species.label):
        name = '%s(%i)' % (species.label, species.index)
        if len(name) <= 10:
            return name

    # Next try the chemical formula
    if len(species.molecule) > 0:
        # Try the chemical formula
        name = '%s(%i)' % (species.molecule[0].getFormula(), species.index)
        if len(name) <= 10:
            return name

    # As a last resort, just use the index
    if species.index >= 0:
        name = 'S(%i)' % (species.index)
        if len(name) <= 10:
            return name

    # If we're here then we just can't come up with a valid Chemkin name
    # for this species, so raise an exception
    raise ChemkinError("Unable to determine valid Chemkin identifier for species %s." % species)

################################################################################

def writeThermoEntry(species):
    """
    Return a string representation of the NASA model readable by Chemkin.
    To use this method you must have exactly two NASA polynomials in your
    model, and you must use the seven-coefficient forms for each.
    """

    thermo = species.thermo
    if not isinstance(thermo, MultiNASA):
        return ''
        raise ChemkinError('Cannot generate Chemkin string for species "%s": Thermodynamics data must be a MultiNASA object.' % species)

    assert len(thermo.polynomials) == 2
    assert thermo.polynomials[0].Tmin < thermo.polynomials[1].Tmin
    assert thermo.polynomials[0].Tmax == thermo.polynomials[1].Tmin
    assert thermo.polynomials[0].cm2 == 0 and thermo.polynomials[0].cm1 == 0
    assert thermo.polynomials[1].cm2 == 0 and thermo.polynomials[1].cm1 == 0

    # Determine the number of each type of element in the molecule
    elements = ['C','H','N','O']; elementCounts = [0,0,0,0]
    for atom in species.molecule[0].atoms:
        # The atom itself
        symbol = atom.element.symbol
        if symbol not in elements:
            elements.append(symbol)
            elementCounts.append(1)
        else:
            elementCounts[elements.index(symbol)] += 1
        # Also handle implicit hydrogen atoms
        symbol = 'H'
        if symbol not in elements:
            elements.append(symbol)
            elementCounts.append(atom.implicitHydrogens)
        else:
            elementCounts[elements.index(symbol)] += atom.implicitHydrogens
    # Remove elements with zero count
    index = 0
    while index < len(elementCounts):
        if elementCounts[index] == 0:
            del elements[index]
            del elementCounts[index]
        else:
            index += 1

    # Line 1
    string = '%-16s        ' % (getSpeciesIdentifier(species))
    if len(elements) <= 4:
        # Use the original Chemkin syntax for the element counts
        for symbol, count in zip(elements, elementCounts):
            string += '%-2s%3i' % (symbol, count)
        string += '     ' * (4 - len(elements))
    else:
        string += '     ' * 4
    string += 'G%10.3f%10.3f%8.2f      1' % (thermo.polynomials[0].Tmin, thermo.polynomials[1].Tmax, thermo.polynomials[0].Tmax)
    if len(elements) > 4:
        string += '&\n'
        # Use the new-style Chemkin syntax for the element counts
        # This will only be recognized by Chemkin 4 or later
        for symbol, count in zip(elements, elementCounts):
            string += '%-2s%3i' % (symbol, count)
    string += '\n'

    # Line 2
    string += '%15.8E%15.8E%15.8E%15.8E%15.8E    2\n' % (thermo.polynomials[0].c0, thermo.polynomials[0].c1, thermo.polynomials[0].c2, thermo.polynomials[0].c3, thermo.polynomials[0].c4)

    # Line 3
    string += '%15.8E%15.8E%15.8E%15.8E%15.8E    3\n' % (thermo.polynomials[0].c5, thermo.polynomials[0].c6, thermo.polynomials[1].c0, thermo.polynomials[1].c1, thermo.polynomials[1].c2)

    # Line 4
    string += '%15.8E%15.8E%15.8E%15.8E                   4\n' % (thermo.polynomials[1].c3, thermo.polynomials[1].c4, thermo.polynomials[1].c5, thermo.polynomials[1].c6)

    return string

################################################################################

def writeKineticsEntry(reaction):
    """
    Return a string representation of the reaction as used in a Chemkin
    file.
    """

    kinetics = reaction.kinetics
    numReactants = len(reaction.reactants)
    
    thirdBody = ''
    if kinetics.isPressureDependent():
        if isinstance(kinetics, ThirdBody) and not isinstance(kinetics, Lindemann) and not isinstance(kinetics, Troe):
            thirdBody = '+M'
        else:
            thirdBody = '(+M)'
    
    string = '+'.join([getSpeciesIdentifier(reactant) for reactant in reaction.reactants])
    string += thirdBody
    string += '=>' if not reaction.reversible else '='
    string += '+'.join([getSpeciesIdentifier(product) for product in reaction.products])
    string += thirdBody

    string = '%-52s' % string

    if isinstance(kinetics, Arrhenius):
        string += '%9.3e %9.3f %9.3f' % (
            kinetics.A / (kinetics.T0 ** kinetics.n) * 1.0e6 ** (numReactants - 1),
            kinetics.n,
            kinetics.Ea / 4184.
        )
    elif isinstance(kinetics, ThirdBody):
        arrhenius = kinetics.arrheniusHigh
        string += '%9.3e %9.3f %9.3f' % (
            arrhenius.A / (arrhenius.T0 ** arrhenius.n) * 1.0e6 ** (numReactants - 1),
            arrhenius.n,
            arrhenius.Ea / 4184.
        )
    else:
        # Print dummy values that Chemkin parses but ignores
        string += '%9.3e %9.3f %9.3f' % (1, 0, 0)

    string += '\n'

    if isinstance(kinetics, ThirdBody):
        # Write collider efficiencies
        for collider, efficiency in kinetics.efficiencies.iteritems():
            string += '%s/%4.2f/ ' % (getSpeciesIdentifier(collider), efficiency)
        string += '\n'
        
        if isinstance(kinetics, Lindemann):
            # Write low-P kinetics
            arrhenius = kinetics.arrheniusLow
            string += '    LOW/ %9.3e %9.3f %9.3f/\n' % (
                arrhenius.A / (arrhenius.T0 ** arrhenius.n) * 1.0e6 ** (numReactants - 1),
                arrhenius.n,
                arrhenius.Ea / 4184.
            )
            if isinstance(kinetics, Troe):
                # Write Troe parameters
                if kinetics.T2 == 1e100:
                    string += '    TROE/ %9.3e %9.3f %9.3f/\n' % (kinetics.alpha, kinetics.T3, kinetics.T1)
                else:
                    string += '    TROE/ %9.3e %9.3f %9.3f %9.3f/\n' % (kinetics.alpha, kinetics.T3, kinetics.T1, kinetics.T2)
    elif isinstance(kinetics, PDepArrhenius):
        for P, arrhenius in zip(kinetics.pressures, kinetics.arrhenius):
            string += '    PLOG/ %9.3f %9.3f %9.3f %9.3f/\n' % (P / 101325.,
                arrhenius.A / (arrhenius.T0 ** arrhenius.n) * 1.0e6 ** (numReactants - 1),
                arrhenius.n,
                arrhenius.Ea / 4184.
            )
    elif isinstance(kinetics, Chebyshev):
        string += '    TCHEB/ %9.3f %9.3f/\n' % (kinetics.Tmin, kinetics.Tmax)
        string += '    PCHEB/ %9.3f %9.3f/\n' % (kinetics.Pmin / 101325., kinetics.Pmax / 101325.)
        string += '    CHEB/ %i %i/\n' % (kinetics.degreeT, kinetics.degreeP)
        if kinetics.degreeP < 6:
            coeffs = kinetics.coeffs.copy()
            coeffs[0,0] += 6 * (numReactants - 1)
            for i in range(kinetics.degreeT):
                string += '    CHEB/'
                for j in range(kinetics.degreeP):
                    string += ' %12.3e' % (coeffs[i,j])
                string += '/\n'
        else:
            coeffs = []
            for i in range(kinetics.degreeT):
                for j in range(kinetics.degreeP):
                    coeffs.append(kinetics.coeffs[i,j])
            coeffs[0] += 6 * (numReactants - 1)
            for i in range(len(coeffs)):
                if i % 5 == 0: string += '    CHEB/'
                string += ' %12.3e' % (kinetics.coeffs[i,j])
                if i % 5 == 4: string += '/\n'

    return string

################################################################################

def saveChemkinFile(path, species, reactions):
    """
    Save a Chemkin input file to `path` on disk containing the provided lists
    of `species` and `reactions`.
    """
    f = open(path, 'w')

    # Elements section
    f.write('ELEMENTS H C O N Ne Ar He Si S END\n\n')

    # Species section
    f.write('SPECIES\n')
    for spec in species:
        label = getSpeciesIdentifier(spec)
        f.write('    %-16s    ! %s\n' % (label, str(spec)))
    f.write('END\n\n\n\n')

    # Thermodynamics section
    f.write('THERM ALL\n')
    f.write('    300.000  1000.000  5000.000\n\n')
    for spec in species:
        f.write(writeThermoEntry(spec))
        f.write('\n')
    f.write('END\n\n\n\n')

    ## Transport section would go here
    #f.write('TRANSPORT\n')
    #f.write('END\n\n')

    # Reactions section
    f.write('REACTIONS    KCAL/MOLE   MOLES\n\n')
    for rxn in reactions:
        f.write(writeKineticsEntry(rxn))
        # Don't forget to mark duplicates!
        f.write('\n')
    f.write('END\n\n')

    f.close()
