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
import logging
import os.path
from thermo import MultiNASA
from kinetics import *
from reaction import Reaction
from species import Species
from thermo import NASA, MultiNASA
from quantity import constants, Quantity

__chemkin_reaction_count = None
    
################################################################################

class ChemkinError(Exception):
    """
    An exception class for exceptional behavior involving Chemkin files. Pass a
    string describing the circumstances that caused the exceptional behavior.
    """
    pass

################################################################################

def readThermoEntry(entry):
    """
    Read a thermodynamics `entry` for one species in a Chemkin file. Returns
    the label of the species and the thermodynamics model as a 
    :class:`MultiNASA` object.
    """
    lines = entry.splitlines()
    species = str(lines[0][0:24].split()[0].strip())
        
    # Extract the NASA polynomial coefficients
    # Remember that the high-T polynomial comes first!
    try:
        Tmin = float(lines[0][45:55].strip())
        Tmax = float(lines[0][55:65].strip())
        Tint = float(lines[0][65:75].strip())
    
        a0_high = float(lines[1][0:15].strip())
        a1_high = float(lines[1][15:30].strip())
        a2_high = float(lines[1][30:45].strip())
        a3_high = float(lines[1][45:60].strip())
        a4_high = float(lines[1][60:75].strip())
    
        a5_high = float(lines[2][0:15].strip())
        a6_high = float(lines[2][15:30].strip())
        a0_low = float(lines[2][30:45].strip())
        a1_low = float(lines[2][45:60].strip())
        a2_low = float(lines[2][60:75].strip())
    
        a3_low = float(lines[3][0:15].strip())
        a4_low = float(lines[3][15:30].strip())
        a5_low = float(lines[3][30:45].strip())
        a6_low = float(lines[3][45:60].strip())
    except (IndexError, ValueError):
        raise ChemkinError('Error while reading thermo entry for species {0}'.format(species))
    
    # Construct and return the thermodynamics model
    thermo = MultiNASA(
        polynomials = [
            NASA(Tmin=(Tmin,"K"), Tmax=(Tint,"K"), coeffs=[a0_low, a1_low, a2_low, a3_low, a4_low, a5_low, a6_low]),
            NASA(Tmin=(Tint,"K"), Tmax=(Tmax,"K"), coeffs=[a0_high, a1_high, a2_high, a3_high, a4_high, a5_high, a6_high])
        ],
        Tmin = (Tmin,"K"),
        Tmax = (Tmax,"K"),
    )

    return species, thermo

################################################################################

def readKineticsEntry(entry, speciesDict, energyUnits, moleculeUnits):
    """
    Read a kinetics `entry` for a single reaction as loaded from a Chemkin
    file. The associated mapping of labels to species `speciesDict` should also
    be provided. Returns a :class:`Reaction` object with the reaction and its
    associated kinetics.
    """
    
    if energyUnits.lower() in ['kcal/mole', 'kcal/mol']:
        energyFactor = 1.0 
    elif energyUnits.lower() in ['cal/mole', 'cal/mol']:
        energyFactor = 0.001 
    else:
        raise ChemkinError('Unexpected energy units "{0}" in reaction block.'.format(energyUnits))
    if moleculeUnits.lower() in ['moles']:
        moleculeFactor = 1.0 
    else:
        raise ChemkinError('Unexpected molecule units "{0}" in reaction block.'.format(energyUnits))
    
    lines = entry.strip().splitlines()
    
    # Extract the reaction equation
    reaction = str(lines[0][0:52].strip())
    thirdBody = False
    duplicate = False
    
    # Split the reaction equation into reactants and products
    reversible = True
    reactants, products = reaction.split('=')
    if '=>' in reaction:
        products = products[1:]
        reversible = False
    if '(+M)' in reactants: reactants = reactants.replace('(+M)','')
    if '(+m)' in reactants: reactants = reactants.replace('(+m)','')
    if '(+M)' in products:  products = products.replace('(+M)','')
    if '(+m)' in products:  products = products.replace('(+m)','')
    
    # Create a new Reaction object for this reaction
    reaction = Reaction(reactants=[], products=[], reversible=reversible)
    
    # Convert the reactants and products to Species objects using the speciesDict
    for reactant in reactants.split('+'):
        reactant = reactant.strip()
        if reactant == 'M':
            thirdBody = True
        elif reactant not in speciesDict:
            raise ChemkinError('Unexpected reactant "{0}" in reaction {1}.'.format(reactant, reaction))
        else:
            reaction.reactants.append(speciesDict[reactant])
    for product in products.split('+'):
        product = product.strip()
        if product.upper() == 'M':
            pass
        elif product not in speciesDict:
            raise ChemkinError('Unexpected product "{0}" in reaction {1}.'.format(product, reaction))
        else:
            reaction.products.append(speciesDict[product])
    
    # Determine the appropriate units for k(T) and k(T,P) based on the number of reactants
    # This assumes elementary kinetics for all reactions
    if len(reaction.reactants) + (1 if thirdBody else 0) == 3:
        kunits = "cm^6/(mol^2*s)"
        klow_units = "cm^9/(mol^3*s)"
    elif len(reaction.reactants) + (1 if thirdBody else 0) == 2:
        kunits = "cm^3/(mol*s)" 
        klow_units = "cm^6/(mol^2*s)"
    elif len(reaction.reactants) + (1 if thirdBody else 0) == 1:
        kunits = "s^-1" 
        klow_units = "cm^3/(mol*s)"
    else:
        raise ChemkinError('Invalid number of reactant species for reaction {0}.'.format(reaction))
    
    # The rest of the first line contains the high-P limit Arrhenius parameters (if available)
    tokens = lines[0][52:].split()
    arrheniusHigh = Arrhenius(
        A = (float(tokens[0].strip()),kunits),
        n = float(tokens[1].strip()),
        Ea = (float(tokens[2].strip()) * energyFactor,"kcal/mol"),
        T0 = (1,"K"),
    )
    
    if len(lines) == 1:
        # If there's only one line then we know to use the high-P limit kinetics as-is
        reaction.kinetics = arrheniusHigh
    else:
        # There's more kinetics information to be read
        arrheniusLow = None
        troe = None
        lindemann = None
        chebyshev = None
        pdepArrhenius = None
        efficiencies = {}
        chebyshevCoeffs = []
    
        # Note that the subsequent lines could be in any order
        for line in lines[1:]:
            tokens = line.split('/')
            if 'DUP' in line:
                # Duplicate reaction
                duplicate = True
            
            elif 'LOW' in line:
                # Low-pressure-limit Arrhenius parameters
                tokens = tokens[1].split()
                arrheniusLow = Arrhenius(
                    A = (float(tokens[0].strip()),klow_units),
                    n = float(tokens[1].strip()),
                    Ea = (float(tokens[2].strip()) * energyFactor,"kcal/mol"),
                    T0 = (1,"K"),
                )
            
            elif 'TROE' in line:
                # Troe falloff parameters
                tokens = tokens[1].split()
                alpha = float(tokens[0].strip())
                T3 = float(tokens[1].strip())
                T1 = float(tokens[2].strip())
                try:
                    T2 = float(tokens[3].strip())
                except (IndexError, ValueError):
                    T2 = None
                
                troe = Troe(
                    alpha = (alpha,''),
                    T3 = (T3,"K"),
                    T1 = (T1,"K"),
                    T2 = (T2,"K") if T2 is not None else None,
                )
            
            elif 'CHEB' in line:
                # Chebyshev parameters
                if chebyshev is None:
                    chebyshev = Chebyshev()
                tokens = [t.strip() for t in tokens]
                if 'TCHEB' in line:
                    index = tokens.index('TCHEB')
                    tokens2 = tokens[index+1].split()
                    chebyshev.Tmin = Quantity(float(tokens2[0].strip()),"K")
                    chebyshev.Tmax = Quantity(float(tokens2[1].strip()),"K")
                if 'PCHEB' in line:
                    index = tokens.index('PCHEB')
                    tokens2 = tokens[index+1].split()
                    chebyshev.Pmin = Quantity(float(tokens2[0].strip()),"atm")
                    chebyshev.Pmax = Quantity(float(tokens2[1].strip()),"atm")
                if 'TCHEB' in line or 'PCHEB' in line:
                    pass
                elif chebyshev.degreeT == 0 or chebyshev.degreeP == 0:
                    tokens2 = tokens[1].split()
                    chebyshev.degreeT = int(tokens2[0].strip())
                    chebyshev.degreeP = int(tokens2[1].strip())
                    chebyshev.coeffs = numpy.zeros((chebyshev.degreeT,chebyshev.degreeP), numpy.float64)
                else:
                    tokens2 = tokens[1].split()
                    chebyshevCoeffs.extend([float(t.strip()) for t in tokens2])
                    
            elif 'PLOG' in line:
                # Pressure-dependent Arrhenius parameters
                if pdepArrhenius is None:
                    pdepArrhenius = []
                tokens = tokens[1].split()
                pdepArrhenius = [float(tokens[0].strip()), Arrhenius(
                    A = (float(tokens[1].strip()),kunits),
                    n = float(tokens[2].strip()),
                    Ea = (float(tokens[3].strip()) * energyFactor,"kcal/mol"),
                    T0 = (1,"K"),
                )]

            else:
                # Assume a list of collider efficiencies
                for collider, efficiency in zip(tokens[0::2], tokens[1::2]):
                    efficiencies[speciesDict[collider.strip()]] = float(efficiency.strip())
    
        # Decide which kinetics to keep and store them on the reaction object
        # Only one of these should be true at a time!
        if chebyshev is not None:
            if chebyshev.Tmin is None or chebyshev.Tmax is None:
                raise ChemkinError('Missing TCHEB line for reaction {0}'.format(reaction))
            if chebyshev.Pmin is None or chebyshev.Pmax is None:
                raise ChemkinError('Missing PCHEB line for reaction {0}'.format(reaction))
            index = 0
            for t in range(chebyshev.degreeT):
                for p in range(chebyshev.degreeP):
                    chebyshev.coeffs[t,p] = chebyshevCoeffs[index]
                    index += 1
            reaction.kinetics = chebyshev
        elif pdepArrhenius is not None:
            reaction.kinetics = PDepArrhenius(
                pressures = ([P for P, arrh in pdepArrhenius],"atm"),
                arrhenius = [arrh in pdepArrhenius],
            )
        elif troe is not None:
            troe.arrheniusHigh = arrheniusHigh
            troe.arrheniusLow = arrheniusLow
            troe.efficiencies = efficiencies
            reaction.kinetics = troe
        elif arrheniusLow is not None:
            reaction.kinetics = Lindemann(arrheniusHigh=arrheniusHigh, arrheniusLow=arrheniusLow)
            reaction.kinetics.efficiencies = efficiencies
        elif thirdBody:
            reaction.kinetics = ThirdBody(arrheniusHigh=arrheniusHigh)
            reaction.kinetics.efficiencies = efficiencies
        elif duplicate:
            reaction.kinetics = arrheniusHigh
        else:
            raise ChemkinError('Unable to determine pressure-dependent kinetics for reaction {0}.'.format(reaction))
   
    return reaction

################################################################################

def loadChemkinFile(path, dictionaryPath=None):
    """
    Load a Chemkin input file to `path` on disk, returning lists of the species
    and reactions in the Chemkin file.
    """
    
    speciesList = []; speciesDict = {}
    reactionList = []

    # If the dictionary path is given, the read it and generate Molecule objects
    # You need to append an additional adjacency list for nonreactive species, such
    # as N2, or else the species objects will not store any structures for the final
    # HTML output.
    if dictionaryPath:
        with open(dictionaryPath, 'r') as f:
            adjlist = ''
            for line in f:
                if line.strip() == '' and adjlist.strip() != '':
                    # Finish this adjacency list
                    species = Species().fromAdjacencyList(adjlist)
                    species.generateResonanceIsomers()
                    speciesDict[species.label] = species
                    adjlist = ''
                else:
                    if '//' in line:
                        index = line.index('//')
                        line = line[0:index]
                    adjlist += line
    
    def removeCommentFromLine(line):
        if '!' in line:
            index = line.index('!')
            comment = line[index+1:-1]
            line = line[0:index] + '\n'
            return line, comment
        else:
            comment = ''
            return line, comment

    def checkDuplicateKinetics(reaction, kinetics,comments,dupReactionList,reactionList):
        if 'DUP' in kinetics:
            kinetics = kinetics.replace('\nDUP','')
            reaction = readKineticsEntry(kinetics,speciesDict,energyUnits,moleculeUnits)
            reaction.kinetics.comment = comments
            if dupReactionList:
                if not reaction.hasTemplate(dupReactionList[-1].reactants,dupReactionList[-1].products):
                    # It's not the same kind of duplicate reaction
                    oldReactionKinetics = MultiKinetics()
                    for item in dupReactionList:
                        oldReactionKinetics.kineticsList.append(item.kinetics)
                    oldReaction = dupReactionList[0]
                    oldReaction.kinetics = oldReactionKinetics
                    reactionList.append(oldReaction)
                    dupReactionList = []
            dupReactionList.append(reaction)
            kinetics = ''
            comments = ''
            return reaction, kinetics, comments, dupReactionList, reactionList

        else:
            # No more duplicate reactions
            if dupReactionList:
                # add previous reaction if they were duplicate reactions
                oldReactionKinetics = MultiKinetics()
                for item in dupReactionList:
                    oldReactionKinetics.kineticsList.append(item.kinetics)
                oldReaction = dupReactionList[0]
                oldReaction.kinetics = oldReactionKinetics
                reactionList.append(oldReaction)
                dupReactionList = []
            # add this new, nonduplicate reaction
            reaction = readKineticsEntry(kinetics,speciesDict,energyUnits,moleculeUnits)
            reaction.kinetics.comment = comments
            print reaction.kinetics
            reactionList.append(reaction)
            kinetics = ''
            comments = ''
            return reaction, kinetics, comments, dupReactionList, reactionList

    with open(path, 'r') as f:
    
        line = f.readline()
        while line != '':        
            line = removeCommentFromLine(line)[0]
            line = line.strip()
            tokens = line.split()
            
            if 'SPECIES' in line:
                # List of species identifiers
                index = tokens.index('SPECIES')
                tokens = tokens[index+1:]
                while 'END' not in tokens:
                    line = f.readline()
                    line = removeCommentFromLine(line)[0]
                    line = line.strip()
                    tokens.extend(line.split())
                
                for token in tokens:
                    if token == 'END':
                        break
                    if token in speciesDict:
                        species = speciesDict[token]
                    else:
                        species = Species(label=token)
                        speciesDict[token] = species
                    speciesList.append(species)         
                
            elif 'THERM' in line:
                # List of thermodynamics (hopefully one per species!)
                line = f.readline()
                thermo = ''
                while line != '' and 'END' not in line:
                    line = removeCommentFromLine(line)[0]
                    if len(line) >= 80:
                        if line[79] in ['1', '2', '3', '4']:
                            thermo += line
                            if line[79] == '4':
                                label, thermo = readThermoEntry(thermo)
                                try:
                                    speciesDict[label].thermo = thermo
                                except KeyError:
                                    if label in ['Ar', 'N2', 'He', 'Ne']:
                                        pass
                                    else:
                                        raise ChemkinError('Unexpected species "{0}" while reading thermodynamics entry.'.format(label))
                                thermo = ''
                    line = f.readline()
                
            elif 'REACTIONS' in line:
                # Reactions section
                energyUnits, moleculeUnits = tokens[1:3]

                line = f.readline()
                kinetics = ''
                comments = ''
                dupReactionList = []
                reaction = None
                while line != '' and 'END' not in line:
                    line, comment = removeCommentFromLine(line)
                    if '=' in line and kinetics.strip() != '':
                        reaction, kinetics,comments,dupReactionList,reactionList = checkDuplicateKinetics(reaction, kinetics, comments, dupReactionList, reactionList)

                    kinetics += line
                    comments += comment
                    
                    line = f.readline()
                # Don't forget the last reaction!
                if kinetics.strip() != '':
                    reaction, kinetics,comments,dupReactionList,reactionList = checkDuplicateKinetics(reaction, kinetics, comments, dupReactionList, reactionList)
                    if dupReactionList:
                    # add previous reaction if they were duplicate reactions
                        oldReactionKinetics = MultiKinetics()
                        for item in dupReactionList:
                            oldReactionKinetics.kineticsList.append(item.kinetics)
                        oldReaction = dupReactionList[0]
                        oldReaction.kinetics = oldReactionKinetics
                        reactionList.append(oldReaction)
                        dupReactionList = []

            line = f.readline()

    from rmgpy.data.kinetics import LibraryReaction,TemplateReaction, KineticsLibrary, KineticsFamily
    from rmgpy.rmg.pdep import PDepReaction, PDepNetwork
    newReactionList = []
    index = 0
    # Create sources for each reaction to include the additional info that it came from RMG-Java
    # and so it will generate family categories with the saveOutputHTML function
    for reaction in reactionList:
        index += 1
        
        #duplicate reactions only occur in libraries
        if reaction.kinetics.comment:
       
            comment = reaction.kinetics.comment
            if comment.find('PDepNetwork') > 0:
                number = comment.split(' ')[3][1:]
                network = PDepNetwork(index = int(number))
                newReaction = PDepReaction(reactants = reaction.reactants, products = reaction.products, kinetics = reaction.kinetics, network = network, index = index)
            else:
                comment = comment.split(' ')
                comment0 = comment[0]
                if comment0.find('Library:') > 0:
                    library = KineticsLibrary(label= comment[1])
                    newReaction = LibraryReaction(reactants = reaction.reactants, products = reaction.products, kinetics = reaction.kinetics, library = library, index = index)
                else:
                    family = KineticsFamily(label = comment[0])
                    newReaction = TemplateReaction(reactants = reaction.reactants, products = reaction.products, kinetics = reaction.kinetics, family = family, index = index)
        else:
            print reaction.kinetics.kineticsList[0].comment.split(' ')[1]
            library = KineticsLibrary(label = reaction.kinetics.kineticsList[0].comment.split(' ')[1])
            newReaction = LibraryReaction(reactants = reaction.reactants, products = reaction.products, kinetics = reaction.kinetics, library=library,index=index)

        newReactionList.append(newReaction)

    return speciesList, newReactionList
    
################################################################################

def saveHTMLFile(path):
    """
    Save an output HTML file from the contents of a RMG-Java output folder
    """
    from rmgpy.rmg.model import CoreEdgeReactionModel
    from rmgpy.rmg.output import saveOutputHTML
    chemkinPath= path + '/chemkin/chem.inp'
    dictionaryPath = path + 'RMG_Dictionary.txt'
    model = CoreEdgeReactionModel()
    model.core.species, model.core.reactions = loadChemkinFile(chemkinPath,dictionaryPath)
    outputPath = path + 'output.html'
    speciesPath = path + '/species/'
    if not os.path.isdir(speciesPath):
        os.makedirs(speciesPath)
    saveOutputHTML(outputPath, model)

################################################################################
def getSpeciesIdentifier(species):
    """
    Return a string identifier for the provided `species` that can be used in a
    Chemkin file. Although the Chemkin format allows up to 16 characters for a
    species identifier, this function uses a maximum of 10 to ensure that all
    reaction equations fit in the maximum limit of 52 characters.
    """

    # Special case for inert colliders - just use the label if possible
    if not species.reactive and 0 < len(species.label) < 10:
        return species.label

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

    if species.index == -1:
        # this didn't come from an RMG-job.  It came from a preexisting chemkin file.
        return species.label

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
    string += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}{4:< 15.8E}    2\n'.format(thermo.polynomials[1].c0, thermo.polynomials[1].c1, thermo.polynomials[1].c2, thermo.polynomials[1].c3, thermo.polynomials[1].c4)

    # Line 3
    string += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}{4:< 15.8E}    3\n'.format(thermo.polynomials[1].c5, thermo.polynomials[1].c6, thermo.polynomials[0].c0, thermo.polynomials[0].c1, thermo.polynomials[0].c2)

    # Line 4
    string += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}                   4\n'.format(thermo.polynomials[0].c3, thermo.polynomials[0].c4, thermo.polynomials[0].c5, thermo.polynomials[0].c6)

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
            new_reaction = Reaction( index=reaction.index,
                     reactants=reaction.reactants,
                     products=reaction.products,
                     reversible=reaction.reversible,
                     kinetics=kinetics)
            string += writeKineticsEntry(new_reaction, speciesList)
            string += "DUPLICATE\n"
        return string + "\n"
    
    global __chemkin_reaction_count
    if __chemkin_reaction_count is not None:
        __chemkin_reaction_count += 1
        string += "! Chemkin # {0}. RMG # {1}.\n".format(__chemkin_reaction_count, reaction.index)
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
    elif hasattr(kinetics,'highPlimit') and kinetics.highPlimit is not None:
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

def saveSpeciesDictionary(path, species):
    """
    Save the given list of `species` as adjacency lists in a text file `path` 
    on disk.
    """
    with open(path, 'w') as f:
        for spec in species:
            f.write(spec.molecule[0].toAdjacencyList(label=getSpeciesIdentifier(spec), removeH=True))
            f.write('\n')

def saveChemkinFile(path, species, reactions):
    """
    Save a Chemkin input file to `path` on disk containing the provided lists
    of `species` and `reactions`.
    """
    f = open(path, 'w')
    
    sorted_species = sorted(species, key=lambda species: species.index)

    # Elements section
    f.write('ELEMENTS H C O N Ne Ar He Si S END\n\n')

    # Species section
    f.write('SPECIES\n')
    for spec in sorted_species:
        label = getSpeciesIdentifier(spec)
        f.write('    {0!s:<16}    ! {1}\n'.format(label, str(spec)))
    f.write('END\n\n\n\n')

    # Thermodynamics section
    f.write('THERM ALL\n')
    f.write('    300.000  1000.000  5000.000\n\n')
    for spec in sorted_species:
        f.write(writeThermoEntry(spec))
        f.write('\n')
    f.write('END\n\n\n\n')

    ## Transport section would go here
    #f.write('TRANSPORT\n')
    #f.write('END\n\n')

    # Reactions section
    f.write('REACTIONS    KCAL/MOLE   MOLES\n\n')
    global __chemkin_reaction_count
    __chemkin_reaction_count = 0
    for rxn in reactions:
        f.write(writeKineticsEntry(rxn, speciesList=species))
        # Don't forget to mark duplicates!
        f.write('\n')
    f.write('END\n\n')
    f.close()
    logging.info("Chemkin file contains {0} reactions.".format(__chemkin_reaction_count))
    __chemkin_reaction_count = None
