#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
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
from reaction import Reaction

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
        name = '{0}({1:d})'.format(species.label, species.index)
        if len(name) <= 10:
            return name

    # Next try the chemical formula
    if len(species.molecule) > 0:
        # Try the chemical formula
        name = '{0}({1:d})'.format(species.molecule[0].getFormula(), species.index)
        if len(name) <= 10:
            return name

    # As a last resort, just use the index
    if species.index >= 0:
        name = 'S({0:d})'.format(species.index)
        if len(name) <= 10:
            return name

    # If we're here then we just can't come up with a valid Chemkin name
    # for this species, so raise an exception
    raise ChemkinError("Unable to determine valid Chemkin identifier for species {0}.".format(species))

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
        raise ChemkinError('Cannot generate Chemkin string for species "{0}": Thermodynamics data must be a MultiNASA object.'.format(species))

    assert len(thermo.polynomials) == 2
    assert thermo.polynomials[0].Tmin.value < thermo.polynomials[1].Tmin.value
    assert thermo.polynomials[0].Tmax.value == thermo.polynomials[1].Tmin.value
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
    string = '{0:<16}        '.format(getSpeciesIdentifier(species))
    if len(elements) <= 4:
        # Use the original Chemkin syntax for the element counts
        for symbol, count in zip(elements, elementCounts):
            string += '{0!s:<2}{1:<3d}'.format(symbol, count)
        string += '     ' * (4 - len(elements))
    else:
        string += '     ' * 4
    string += 'G{0:<10.3f}{1:<10.3f}{2:<8.2f}      1'.format(thermo.polynomials[0].Tmin.value, thermo.polynomials[1].Tmax.value, thermo.polynomials[0].Tmax.value)
    if len(elements) > 4:
        string += '&\n'
        # Use the new-style Chemkin syntax for the element counts
        # This will only be recognized by Chemkin 4 or later
        for symbol, count in zip(elements, elementCounts):
            string += '{0!s:<2}{1:<3d}'.format(symbol, count)
    string += '\n'

    # Line 2
    string += '{0:<15.8E}{1:<15.8E}{2:<15.8E}{3:<15.8E}{4:<15.8E}    2\n'.format(thermo.polynomials[1].c0, thermo.polynomials[1].c1, thermo.polynomials[1].c2, thermo.polynomials[1].c3, thermo.polynomials[1].c4)

    # Line 3
    string += '{0:<15.8E}{1:<15.8E}{2:<15.8E}{3:<15.8E}{4:<15.8E}    3\n'.format(thermo.polynomials[1].c5, thermo.polynomials[1].c6, thermo.polynomials[0].c0, thermo.polynomials[0].c1, thermo.polynomials[0].c2)

    # Line 4
    string += '{0:<15.8E}{1:<15.8E}{2:<15.8E}{3:<15.8E}                   4\n'.format(thermo.polynomials[0].c3, thermo.polynomials[0].c4, thermo.polynomials[0].c5, thermo.polynomials[0].c6)

    return string

################################################################################

def writeKineticsEntry(reaction, speciesList):
    """
    Return a string representation of the reaction as used in a Chemkin
    file.
    """
    string = ""
    
    if isinstance(reaction.kinetics, MultiKinetics):
        if reaction.kinetics.comment:
            for line in reaction.kinetics.comment.split("\n"):
                string += "! {0}\n".format(line) 
        for kinetics in reaction.kinetics.kineticsList:
            new_reaction = Reaction(reactants=reaction.reactants,
                     products=reaction.products,
                     reversible=reaction.reversible,
                     kinetics=kinetics)
            string += writeKineticsEntry(new_reaction, speciesList)
            string += "DUPLICATE\n"

    if reaction.kinetics.comment:
        for line in reaction.kinetics.comment.split("\n"):
            string += "! {0}\n".format(line) 
    kinetics = reaction.kinetics
    numReactants = len(reaction.reactants)
    
    thirdBody = ''
    if kinetics.isPressureDependent():
        if isinstance(kinetics, ThirdBody) and not isinstance(kinetics, Lindemann) and not isinstance(kinetics, Troe):
            thirdBody = '+M'
        elif isinstance(kinetics, PDepArrhenius):
            thirdBody = ''
        else:
            thirdBody = '(+M)'
    
    reaction_string = '+'.join([getSpeciesIdentifier(reactant) for reactant in reaction.reactants])
    reaction_string += thirdBody
    reaction_string += '=>' if not reaction.reversible else '='
    reaction_string += '+'.join([getSpeciesIdentifier(product) for product in reaction.products])
    reaction_string += thirdBody
    
    string += '{0!s:<52}'.format(reaction_string)

    if isinstance(kinetics, Arrhenius):
        string += '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(
            kinetics.A.value/ (kinetics.T0.value ** kinetics.n.value) * 1.0e6 ** (numReactants - 1),
            kinetics.n.value,
            kinetics.Ea.value / 4184.
        )
    elif isinstance(kinetics, ThirdBody):
        arrhenius = kinetics.arrheniusHigh
        string += '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(
            arrhenius.A.value / (arrhenius.T0.value ** arrhenius.n.value) * 1.0e6 ** (numReactants - 1),
            arrhenius.n.value,
            arrhenius.Ea.value / 4184.
        )
    elif hasattr(kinetics,'highPlimit'):
        arrhenius = kinetics.highPlimit
        string += '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(
            arrhenius.A.value / (arrhenius.T0.value ** arrhenius.n.value) * 1.0e6 ** (numReactants - 1),
            arrhenius.n.value,
            arrhenius.Ea.value / 4184.
            )
    else:
        # Print dummy values that Chemkin parses but ignores
        string += '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(1, 0, 0)

    string += '\n'

    if isinstance(kinetics, ThirdBody):
        # Write collider efficiencies
        for collider, efficiency in kinetics.efficiencies.iteritems():
            for species in speciesList:
                if any([collider.isIsomorphic(molecule) for molecule in species.molecule]):
                    string += '{0!s}/{1:<4.2f}/ '.format(getSpeciesIdentifier(species), efficiency)
                    break
        string += '\n'
        
        if isinstance(kinetics, Lindemann):
            # Write low-P kinetics
            arrhenius = kinetics.arrheniusLow
            string += '    LOW/ {0:<9.3e} {1:<9.3f} {2:<9.3f}/\n'.format(
                arrhenius.A.value / (arrhenius.T0.value ** arrhenius.n.value) * 1.0e6 ** (numReactants - 1),
                arrhenius.n.value,
                arrhenius.Ea.value / 4184.
            )
            if isinstance(kinetics, Troe):
                # Write Troe parameters
                if kinetics.T2 is None:
                    string += '    TROE/ {0:<9.3e} {1:<9.3g} {2:<9.3g}/\n'.format(kinetics.alpha.value, kinetics.T3.value, kinetics.T1.value)
                else:
                    string += '    TROE/ {0:<9.3e} {1:<9.3g} {2:<9.3g} {3:<9.3g}/\n'.format(kinetics.alpha.value, kinetics.T3.value, kinetics.T1.value, kinetics.T2.value)
    elif isinstance(kinetics, PDepArrhenius):
        for P, arrhenius in zip(kinetics.pressures.values, kinetics.arrhenius):
            string += '    PLOG/ {0:<9.3f} {1:<9.3e} {2:<9.3f} {3:<9.3f}/\n'.format(P / 101325.,
                arrhenius.A.value / (arrhenius.T0.value ** arrhenius.n.value) * 1.0e6 ** (numReactants - 1),
                arrhenius.n.value,
                arrhenius.Ea.value / 4184.
            )
    elif isinstance(kinetics, Chebyshev):
        string += '    TCHEB/ {0:<9.3f} {1:<9.3f}/\n'.format(kinetics.Tmin.value, kinetics.Tmax.value)
        string += '    PCHEB/ {0:<9.3f} {1:<9.3f}/\n'.format(kinetics.Pmin.value / 101325., kinetics.Pmax.value / 101325.)
        string += '    CHEB/ {0:d} {1:d}/\n'.format(kinetics.degreeT, kinetics.degreeP)
        if kinetics.degreeP < 6:
            coeffs = kinetics.coeffs.copy()
            coeffs[0,0] += 6 * (numReactants - 1)
            for i in range(kinetics.degreeT):
                string += '    CHEB/'
                for j in range(kinetics.degreeP):
                    string += ' {0:<12.3e}'.format(coeffs[i,j])
                string += '/\n'
        else:
            coeffs = []
            for i in range(kinetics.degreeT):
                for j in range(kinetics.degreeP):
                    coeffs.append(kinetics.coeffs[i,j])
            coeffs[0] += 6 * (numReactants - 1)
            for i in range(len(coeffs)):
                if i % 5 == 0: string += '    CHEB/'
                string += ' {0:<12.3e}'.format(kinetics.coeffs[i,j])
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
        f.write('    {0!s:<16}    ! {1}\n'.format(label, str(spec)))
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
        f.write(writeKineticsEntry(rxn, speciesList=species))
        # Don't forget to mark duplicates!
        f.write('\n')
    f.write('END\n\n')

    f.close()
