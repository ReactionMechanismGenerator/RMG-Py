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
import shutil
import math
import re
import logging
import textwrap
import os.path
import numpy

import rmgpy.kinetics as _kinetics
from rmgpy.reaction import Reaction
#from species import Species
from rmgpy.rmg.model import Species
from rmgpy.rmg.pdep import PDepReaction
from rmgpy.thermo import NASAPolynomial, NASA
import rmgpy.constants as constants
from rmgpy.quantity import Quantity
from rmgpy.data.base import Entry 
from rmgpy.data.kinetics.library import LibraryReaction, KineticsLibrary
from rmgpy.data.kinetics.family import TemplateReaction, KineticsFamily
from rmgpy.rmg.pdep import PDepNetwork
from rmgpy.molecule import Molecule
from rmgpy.transport import TransportData

__chemkin_reaction_count = None
    
################################################################################

class ChemkinError(Exception):
    """
    An exception class for exceptional behavior involving Chemkin files. Pass a
    string describing the circumstances that caused the exceptional behavior.
    """
    pass

################################################################################

def Ffloat(string):
    """
    Parse a Fortran-ish string into a float, like "1.00D 03" or "1.00d+03"
    """
    return float(string.replace("e", "E").replace("d", "D").replace("D", "E").replace("E ", "E+"))
    
def readThermoEntry(entry, Tmin=0, Tint=0, Tmax=0):
    """
    Read a thermodynamics `entry` for one species in a Chemkin file. Returns
    the label of the species and the thermodynamics model as a :class:`NASA`
    object.
    
    Format specification at http://www2.galcit.caltech.edu/EDL/public/formats/chemkin.html
    """
    lines = entry.splitlines()
    species = str(lines[0][0:18].split()[0].strip())
    
    comment = lines[0][len(species):24].strip()
    formula = {}
    for i in [24, 29, 34, 39, 74]:
        element,count = lines[0][i:i+2].strip(), lines[0][i+2:i+5].strip()
        if element:
            try:
                count = int(count)
            except ValueError:
                # Chemkin allows float values for the number of atoms, so try this next.
                try:
                    count = int(float(count))
                except ValueError:
                    logging.info("Trouble reading line '{0}' element segment '{1}'".format(lines[0].strip(),lines[0][i:i+5]))
                    if count == '' and re.match('\.?0*', element):
                        if i == 74: logging.info("Assuming it's spillover from Tint, and ignoring.")
                        else: logging.warning("Assuming it's not meant to be there, although it would be good to fix the chemkin file.")
                        count = 0
                    else:
                        raise
            if count != 0: # Some people put garbage elements in, with zero count. Ignore these. Allow negative counts though (eg. negative one electron)
                formula[element] = count
    phase = lines[0][44]
    if phase.upper() != 'G':
        logging.warning("Was expecting gas phase thermo data for {0}. Skipping thermo data.".format(species))
        return species, None, None
    
    # Extract the NASA polynomial coefficients
    # Remember that the high-T polynomial comes first!
    try:
        try:
            Tmin = float(lines[0][45:55].strip())
        except ValueError:
            pass
        try:
            Tmax = float(lines[0][55:65].strip())
        except ValueError:
            pass
        try:
            Tint = float(lines[0][65:73].strip())
        except ValueError:
            pass
        a0_high = Ffloat(lines[1][0:15].strip())
        a1_high = Ffloat(lines[1][15:30].strip())
        a2_high = Ffloat(lines[1][30:45].strip())
        a3_high = Ffloat(lines[1][45:60].strip())
        a4_high = Ffloat(lines[1][60:75].strip())
    
        a5_high = Ffloat(lines[2][0:15].strip())
        a6_high = Ffloat(lines[2][15:30].strip())
        a0_low = Ffloat(lines[2][30:45].strip())
        a1_low = Ffloat(lines[2][45:60].strip())
        a2_low = Ffloat(lines[2][60:75].strip())
    
        a3_low = Ffloat(lines[3][0:15].strip())
        a4_low = Ffloat(lines[3][15:30].strip())
        a5_low = Ffloat(lines[3][30:45].strip())
        a6_low = Ffloat(lines[3][45:60].strip())
    except (IndexError, ValueError), e:
        logging.warning('Error while reading thermo entry for species {0}'.format(species))
        logging.warning(e.message)
        return species, None, None
    
    # Construct and return the thermodynamics model
    thermo = NASA(
        polynomials = [
            NASAPolynomial(Tmin=(Tmin,"K"), Tmax=(Tint,"K"), coeffs=[a0_low, a1_low, a2_low, a3_low, a4_low, a5_low, a6_low]),
            NASAPolynomial(Tmin=(Tint,"K"), Tmax=(Tmax,"K"), coeffs=[a0_high, a1_high, a2_high, a3_high, a4_high, a5_high, a6_high])
        ],
        Tmin = (Tmin,"K"),
        Tmax = (Tmax,"K"),
    )
    if comment:
        thermo.comment = comment

    return species, thermo, formula

################################################################################

def readKineticsEntry(entry, speciesDict, Aunits, Eunits):
    """
    Read a kinetics `entry` for a single reaction as loaded from a Chemkin
    file. The associated mapping of labels to species `speciesDict` should also
    be provided. Returns a :class:`Reaction` object with the reaction and its
    associated kinetics.
    """
    Afactor = {
        'cm^3/(mol*s)': 1.0e6,
        'cm^3/(molecule*s)': 1.0e6,
        'm^3/(mol*s)': 1.0,
        'm^3/(molecule*s)': 1.0,
    }[Aunits[2]]
    
    lines = entry.strip().splitlines()
    
    # The first line contains the reaction equation and a set of
    # modified Arrhenius parameters
    reaction, thirdBody, kinetics, kunits, klow_units = _readKineticsReaction(
        line=lines[0], speciesDict=speciesDict, Aunits=Aunits, Eunits=Eunits)

    if len(lines) == 1 and not thirdBody:
        # If there's only one line then we know to use the high-P limit kinetics as-is
        reaction.kinetics = kinetics['arrhenius high']
    else:
        # There's more kinetics information to be read
        kinetics.update({
            'chebyshev coefficients': [],
            'efficiencies': {},
        })

        # Note that the subsequent lines could be in any order
        for line in lines[1:]:
            kinetics = _readKineticsLine(
                line=line, reaction=reaction, speciesDict=speciesDict, Eunits=Eunits,
                kunits=kunits, klow_units=klow_units,
                kinetics=kinetics)

        # Decide which kinetics to keep and store them on the reaction object
        # Only one of these should be true at a time!
        if 'chebyshev' in kinetics:
            chebyshev = kinetics['chebyshev']
            if chebyshev.Tmin is None or chebyshev.Tmax is None:
                raise ChemkinError('Missing TCHEB line for reaction {0}'.format(reaction))
            if chebyshev.Pmin is None or chebyshev.Pmax is None:
                raise ChemkinError('Missing PCHEB line for reaction {0}'.format(reaction))
            index = 0
            for t in range(chebyshev.degreeT):
                for p in range(chebyshev.degreeP):
                    chebyshev.coeffs.value_si[t,p] = kinetics[
                        'chebyshev coefficients'][index]
                    index += 1
            # Don't forget to convert the Chebyshev coefficients to SI units!
            # This assumes that s^-1, cm^3/mol*s, etc. are compulsory
            chebyshev.coeffs.value_si[0,0] -= (len(reaction.reactants) - 1) * math.log10(Afactor)
            reaction.kinetics = chebyshev
        elif 'pressure-dependent arrhenius' in kinetics:
            pdepArrhenius = kinetics['pressure-dependent arrhenius']
            # Check for duplicates and combine them to MultiArrhenius objects
            duplicatesToRemove = []
            duplicatesToAdd = []
            for index1 in range(len(pdepArrhenius)):
                reaction1 = pdepArrhenius[index1]
                P1, kinetics1 = reaction1
                if reaction1 in duplicatesToRemove:
                    continue
                for index2 in range(index1+1, len(pdepArrhenius)):
                    reaction2 = pdepArrhenius[index2]
                    P2, kinetics2 = reaction2
                    if P1 == P2:
                        if reaction1 not in duplicatesToRemove:
                            new_kinetics = _kinetics.MultiArrhenius()
                            duplicatesToAdd.append([P1,new_kinetics])
                            new_kinetics.arrhenius = [kinetics1]
                            duplicatesToRemove.append(reaction1)
                        new_kinetics.arrhenius.append(kinetics2)
                        duplicatesToRemove.append(reaction2)
            for item in duplicatesToRemove:
                pdepArrhenius.remove(item)
            pdepArrhenius.extend(duplicatesToAdd)
            
            pdepArrhenius = sorted(pdepArrhenius, key=lambda reaction: reaction[0])  # sort by ascending pressures

            reaction.kinetics = _kinetics.PDepArrhenius(
                pressures = ([P for P, arrh in pdepArrhenius],"atm"),
                arrhenius = [arrh for P, arrh in pdepArrhenius],
            )
        elif 'troe' in kinetics:
            troe = kinetics['troe']
            troe.arrheniusHigh = kinetics['arrhenius high']
            troe.arrheniusLow = kinetics['arrhenius low']
            troe.efficiencies = kinetics['efficiencies']
            reaction.kinetics = troe
        elif thirdBody:
            reaction.kinetics = _kinetics.ThirdBody(
                arrheniusLow=kinetics['arrhenius low'])
            reaction.kinetics.efficiencies = kinetics['efficiencies']
        elif 'arrhenius low' in kinetics:
            reaction.kinetics = _kinetics.Lindemann(
                arrheniusHigh=kinetics['arrhenius high'],
                arrheniusLow=kinetics['arrhenius low'])
            reaction.kinetics.efficiencies = kinetics['efficiencies']
        elif 'explicit reverse' in kinetics or reaction.duplicate:
            # it's a normal high-P reaction - the extra lines were only either REV (explicit reverse) or DUP (duplicate)
            reaction.kinetics = kinetics['arrhenius high']
        else:
            raise ChemkinError(
                'Unable to understand all additional information lines for reaction {0}.'.format(entry))
        
        # These things may *also* be true
        if 'sri' in kinetics:
            reaction.kinetics.comment += "Warning: SRI parameters from chemkin file ignored on import. "

        if 'explicit reverse' in kinetics:
            reaction.kinetics.comment += (
                "Chemkin file stated explicit reverse rate: {0}"
                ).format(kinetics['explicit reverse'])

    return reaction


def _readKineticsReaction(line, speciesDict, Aunits, Eunits):
    """
    Parse the first line of of a Chemkin reaction entry.
    """
    tokens = line.split()
    
    rmg = True
    try:
        float(tokens[-6])
    except (ValueError, IndexError):
        rmg = False
    AuncertaintyType = '+|-'
    if rmg:
        A = float(tokens[-6])
        n = float(tokens[-5])
        Ea = float(tokens[-4])
        try:
            dA = float(tokens[-3])
        except ValueError:
            AuncertaintyType = '*|/'
            dA = float(tokens[-3][1:])
        dn = float(tokens[-2])
        dEa = float(tokens[-1])
        reaction = ''.join(tokens[:-6])
    else:
        A = Ffloat(tokens[-3])
        n = Ffloat(tokens[-2])
        Ea = Ffloat(tokens[-1])
        dA = 0.0
        dn = 0.0
        dEa = 0.0
        reaction = ''.join(tokens[:-3])
    thirdBody = False
    
    # Split the reaction equation into reactants and products
    reversible = True
    reactants, products = reaction.split('=')
    if '<=>' in reaction:
        reactants = reactants[:-1]
        products = products[1:]
    elif '=>' in reaction:
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
        stoichiometry = 1
        if reactant[0].isdigit():
            # This allows for reactions to be of the form 2A=B+C instead of A+A=B+C
            # The implementation below assumes an integer between 0 and 9, inclusive
            stoichiometry = int(reactant[0])
            reactant = reactant[1:]               
        if reactant.upper() == 'M':
            thirdBody = True
        elif reactant not in speciesDict:
            raise ChemkinError('Unexpected reactant "{0}" in reaction {1}.'.format(reactant, reaction))
        else:
            for i in range(stoichiometry):
                reaction.reactants.append(speciesDict[reactant])
    for product in products.split('+'):
        product = product.strip()
        stoichiometry = 1
        if product[0].isdigit():
            # This allows for reactions to be of the form A+B=2C instead of A+B=C+C
            # The implementation below assumes an integer between 0 and 9, inclusive
            stoichiometry = int(product[0])
            product = product[1:]
        if product.upper() == 'M':
            pass
        elif product not in speciesDict:
            if re.match('[0-9.]+',product):
                logging.warning("Looks like reaction {0!r} has fractional stoichiometry, which RMG cannot handle. Ignoring".format(line))
                raise ChemkinError('Skip reaction!')
            raise ChemkinError('Unexpected product "{0}" in reaction {1}.'.format(product, reaction))
        else:
            for i in range(stoichiometry):
                reaction.products.append(speciesDict[product])
    
    # Determine the appropriate units for k(T) and k(T,P) based on the number of reactants
    # This assumes elementary kinetics for all reactions
    try:
        Nreac = len(reaction.reactants) + (1 if thirdBody else 0)
        kunits = Aunits[Nreac]
        klow_units = Aunits[Nreac+1]
    except IndexError:
        raise ChemkinError('Invalid number of reactant species for reaction {0}.'.format(reaction))
    
    key = 'arrhenius low' if thirdBody else 'arrhenius high'
    
    kinetics = {
        key: _kinetics.Arrhenius(
            A = (A,kunits,AuncertaintyType,dA),
            n = (n,'','+|-',dn),
            Ea = (Ea,Eunits,'+|-',dEa),
            T0 = (1,"K"),
        ),
    }
    return (reaction, thirdBody, kinetics, kunits, klow_units)


def _readKineticsLine(line, reaction, speciesDict, Eunits, kunits, klow_units, kinetics):
    """
    Parse the subsequent lines of of a Chemkin reaction entry.
    """
    case_preserved_tokens = line.split('/')
    line = line.upper()
    tokens = line.split('/')

    if 'DUP' in line:
        # Duplicate reaction
        reaction.duplicate = True

    elif 'LOW' in line:
        # Low-pressure-limit Arrhenius parameters
        tokens = tokens[1].split()
        kinetics['arrhenius low'] = _kinetics.Arrhenius(
            A = (float(tokens[0].strip()),klow_units),
            n = float(tokens[1].strip()),
            Ea = (float(tokens[2].strip()),Eunits),
            T0 = (1,"K"),
        )

    elif 'HIGH' in line:
        # What we thought was high, was in fact low-pressure
        kinetics['arrhenius low'] = kinetics['arrhenius high']
        kinetics['arrhenius low'].A = (
            kinetics['arrhenius low'].A.value, klow_units)
        # High-pressure-limit Arrhenius parameters
        tokens = tokens[1].split()
        kinetics['arrhenius high'] = _kinetics.Arrhenius(
            A = (float(tokens[0].strip()),kunits),
            n = float(tokens[1].strip()),
            Ea = (float(tokens[2].strip()),Eunits),
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

        kinetics['troe'] = _kinetics.Troe(
            alpha = alpha,
            T3 = (T3,"K"),
            T1 = (T1,"K"),
            T2 = (T2,"K") if T2 is not None else None,
        )
    elif line.strip().startswith('SRI'):
        kinetics['sri'] = True
        """To define an SRI pressure-dependent reaction, in addition to the LOW or HIGH parameters, the
        keyword SRI followed by three or five parameters must be included in the following order: a, b,
        c, d, and e [Eq. (74)]. The fourth and fifth parameters are options. If only the first three are
        stated, then by default d = 1 and e = 0.
        """
        # see eg. http://www.dipic.unipd.it/faculty/canu/files/Comb/Docs/chemkinCK.pdf

    elif 'CHEB' in line:
        # Chebyshev parameters
        chebyshev = kinetics.get(
            'chebyshev', _kinetics.Chebyshev(kunits=kunits))
        kinetics['chebyshev'] = chebyshev
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
            chebyshev.degreeT = int(float(tokens2[0].strip()))
            chebyshev.degreeP = int(float(tokens2[1].strip()))
            chebyshev.coeffs = numpy.zeros((chebyshev.degreeT,chebyshev.degreeP), numpy.float64)
        else:
            tokens2 = tokens[1].split()
            kinetics['chebyshev coefficients'].extend(
                [float(t.strip()) for t in tokens2])

    elif 'PLOG' in line:
        pdepArrhenius = kinetics.get('pressure-dependent arrhenius', [])
        kinetics['pressure-dependent arrhenius'] = pdepArrhenius
        tokens = tokens[1].split()
        pdepArrhenius.append([float(tokens[0].strip()),
            _kinetics.Arrhenius(
                A = (float(tokens[1].strip()),kunits),
                n = float(tokens[2].strip()),
                Ea = (float(tokens[3].strip()),Eunits),
                T0 = (1,"K"),
        )])
    elif tokens[0].startswith('REV'):
        reverseA = float(tokens[1].split()[0])
        kinetics['explicit reverse'] = line.strip()
        if reverseA == 0:
            logging.info("Reverse rate is 0 so making irreversible for reaction {0}".format(reaction))
            reaction.reversible = False
        else:
            logging.info("Ignoring explicit reverse rate for reaction {0}".format(reaction))

    else:
        # Assume a list of collider efficiencies
        try:
            for collider, efficiency in zip(case_preserved_tokens[0::2], case_preserved_tokens[1::2]):
                try:
                    efficiency = float(efficiency.strip())
                except ValueError:
                    error_msg = "{0!r} doesn't look like a collision efficiency for species {1} in line {2!r}".format(efficiency,collider.strip(),line)
                    logging.error(error_msg)
                    raise ChemkinError(error_msg)
                if collider.strip() in speciesDict:
                    kinetics['efficiencies'][speciesDict[collider.strip()].molecule[0]] = efficiency
                else: # try it with capital letters? Not sure whose malformed chemkin files this is needed for.
                    kinetics['efficiencies'][speciesDict[collider.strip().upper()].molecule[0]] = efficiency
        except IndexError:
            error_msg = 'Could not read collider efficiencies for reaction: {0}.\n'.format(reaction)
            error_msg += 'The following line was parsed incorrectly:\n{0}'.format(line)
            error_msg += "\n(Case-preserved tokens: {0!r} )".format(case_preserved_tokens)
            raise ChemkinError(error_msg)
    return kinetics


def readReactionComments(reaction, comments, read = True):
    """
    Parse the `comments` associated with a given `reaction`. If the comments
    come from RMG (Py or Java), parse them and extract the useful information.
    Return the reaction object based on the information parsed from these
    comments. If `read` if False, the reaction is returned as an "Unclassified"
    LibraryReaction.
    """
    
    if read == False:
        # The chemkin file was not generated by either RMG-Py or RMG-Java, thus, there should be no parsing
        # of the comments.  Instead, return as an unclassified LibraryReaction.
        reaction = LibraryReaction(
            index = reaction.index,
            reactants = reaction.reactants, 
            products = reaction.products, 
            kinetics = reaction.kinetics,
            reversible = reaction.reversible,
            duplicate = reaction.duplicate,
            library = KineticsLibrary(label='Unclassified'),
        )        
        
        return reaction  
    
    atKineticsComments = False
    lines = comments.strip().splitlines()
        
    for line in lines:
        
        
        tokens = line.split()
        if 'Reaction index:' in line:
            # Don't store the reaction indices
            pass
        
        elif 'Template reaction:' in line:
            label = str(tokens[-2])
            template = tokens[-1][1:-1].split(',')
            reaction = TemplateReaction(
                index = reaction.index,
                reactants = reaction.reactants, 
                products = reaction.products, 
                kinetics = reaction.kinetics,
                reversible = reaction.reversible,
                duplicate = reaction.duplicate,
                family = KineticsFamily(label=label),
                template = [Entry(label=g) for g in template],
            )
            
        elif 'Library reaction:' in line or 'Seed mechanism:' in line:
            label = str(tokens[-1])
            reaction = LibraryReaction(
                index = reaction.index,
                reactants = reaction.reactants, 
                products = reaction.products, 
                kinetics = reaction.kinetics,
                reversible = reaction.reversible,
                duplicate = reaction.duplicate,
                library = KineticsLibrary(label=label),
            )   
            
        elif 'PDep reaction:' in line:
            networkIndex = int(tokens[-1][1:])
            reaction = PDepReaction(
                index = reaction.index,
                reactants = reaction.reactants, 
                products = reaction.products, 
                kinetics = reaction.kinetics, 
                reversible = reaction.reversible,
                duplicate = reaction.duplicate,
                network = PDepNetwork(index=networkIndex), 
            )
            
        elif 'Flux pairs:' in line:
            reaction.pairs = []
            for reacStr, prodStr in zip(tokens[2::2], tokens[3::2]):
                if reacStr[-1] == ',': reacStr = reacStr[:-1]
                for reactant in reaction.reactants:
                    if reactant.label == reacStr:
                        break
                else:
                    import pdb; pdb.set_trace()
                    raise ChemkinError('Unexpected species identifier {0} encountered in flux pairs for reaction {1}.'.format(reacStr, reaction))
                if prodStr[-1] == ';': prodStr = prodStr[:-1]
                for product in reaction.products:
                    if product.label == prodStr:
                        break
                else:
                    import pdb; pdb.set_trace()
                    raise ChemkinError('Unexpected species identifier {0} encountered in flux pairs for reaction {1}.'.format(prodStr, reaction))
                reaction.pairs.append((reactant, product))
            assert len(reaction.pairs) == max(len(reaction.reactants), len(reaction.products))

        elif 'Kinetics comments:' in line:
            atKineticsComments = True

        elif atKineticsComments:
            reaction.kinetics.comment += line.strip() + "\n"


        # Comment parsing from old RMG-Java chemkin files
        elif 'PDepNetwork' in line:
            networkIndex = int(tokens[3][1:])
            reaction = PDepReaction(
                index = reaction.index,
                reactants = reaction.reactants, 
                products = reaction.products,
                kinetics = reaction.kinetics,
                reversible = reaction.reversible,
                duplicate = reaction.duplicate,
                network = PDepNetwork(index=networkIndex)
                )
            reaction.kinetics.comment = line

        elif 'ReactionLibrary:' in line or 'Seed Mechanism:' in line:
            label = str(tokens[-1])
            reaction = LibraryReaction(
                index = reaction.index,
                reactants = reaction.reactants, 
                products = reaction.products, 
                kinetics = reaction.kinetics,
                reversible = reaction.reversible,
                duplicate = reaction.duplicate,
                library = KineticsLibrary(label=label),
            )
            reaction.kinetics.comment = line
            
        elif 'exact' in line or 'estimate' in line:
            index1 = line.find('[')
            index2 = line.find(']')
            template = [s.strip() for s in line[index1:index2].split(',')]
            label = str(tokens[0])
            reaction = TemplateReaction(
                index = reaction.index,
                reactants = reaction.reactants, 
                products = reaction.products, 
                kinetics = reaction.kinetics,
                reversible = reaction.reversible,
                duplicate = reaction.duplicate,
                family = KineticsFamily(label=label),
                template = [Entry(label=g) for g in template],
            )
            reaction.kinetics.comment = line

    if not isinstance(reaction, LibraryReaction) and not isinstance(reaction, TemplateReaction) and not isinstance(reaction,PDepReaction):
        reaction = LibraryReaction(
            index = reaction.index,
            reactants = reaction.reactants, 
            products = reaction.products, 
            kinetics = reaction.kinetics,
            reversible = reaction.reversible,
            duplicate = reaction.duplicate,
            library = KineticsLibrary(label='Unclassified'),
        )  
            
    return reaction

################################################################################

def loadSpeciesDictionary(path):
    """
    Load an RMG dictionary - containing species identifiers and the associated
    adjacency lists - from the file located at `path` on disk. Returns a dict
    mapping the species identifiers to the loaded species. Resonance isomers
    for each species are automatically generated.
    """
    speciesDict = {}
    
    with open(path, 'r') as f:
        adjlist = ''
        for line in f:
            if line.strip() == '' and adjlist.strip() != '':
                # Finish this adjacency list
                species = Species().fromAdjacencyList(adjlist)
                species.generateResonanceIsomers()
                label = species.label
                speciesDict[label] = species
                adjlist = ''
            else:
                if "InChI" in line:
                    line = line.split()[0] + '\n'
                if '//' in line:
                    index = line.index('//')
                    line = line[0:index]
                adjlist += line

    return speciesDict

def removeCommentFromLine(line):
    """
    Remove a comment from a line of a Chemkin file or species dictionary file.
    
    Returns the line and the comment.
    If the comment is encoded with latin-1, it is converted to utf-8.
    """
    try:
        index1 = line.index('!')
    except ValueError:
        index1 = len(line)
    try:
        index2 = line.index('//')
    except ValueError:
        index2 = len(line)
    
    index = min(index1, index2)
    comment = line[index+1:-1]
    if index < len(line):
        line = line[0:index] + '\n'

    try:
        ucomment = comment.decode('utf-8')
    except UnicodeDecodeError:
        try:
            ucomment = comment.decode('latin-1')
        except UnicodeDecodeError:
            ucomment = comment.decode('windows-1252', errors='replace')
    # Convert back to utf-8.
    # Other parts of RMG-Py expect a string object not a unicode object,
    # but at least this way we know what the encoding is.
    comment = ucomment.encode('utf-8', 'replace')
    return line, comment

def loadTransportFile(path, speciesDict):
    """
    Load a Chemkin transport properties file located at `path` and store the
    properties on the species in `speciesDict`.
    """
    with open(path, 'r') as f:
        for line0 in f:
            line = removeCommentFromLine(line0)[0]
            line = line.strip()
            if line != '':
                # This line contains an entry, so parse it
                label = line[0:16].strip()
                data = line[16:].split()
                species = speciesDict[label]
                species.transportData = TransportData(
                    sigma = (float(data[2]),'angstrom'),
                    epsilon = (float(data[1]),'K'),
                )
                species.dipoleMoment = (float(data[3]),'De')
                species.polarizability = (float(data[4]),'angstrom^3')
                species.Zrot = (float(data[5]),'')

def loadChemkinFile(path, dictionaryPath=None, transportPath=None, readComments = True, thermoPath = None):
    """
    Load a Chemkin input file located at `path` on disk to `path`, returning lists of the species
    and reactions in the Chemkin file. The 'thermoPath' point to a separate thermo file, or, if 'None' is 
    specified, the function will look for the thermo database within the chemkin mechanism file
    """
    
    speciesList = []; speciesDict = {}; speciesAliases = {}
    reactionList = []

    # If the dictionary path is given, then read it and generate Molecule objects
    # You need to append an additional adjacency list for nonreactive species, such
    # as N2, or else the species objects will not store any structures for the final
    # HTML output.
    if dictionaryPath:
        speciesDict = loadSpeciesDictionary(dictionaryPath)
    
    with open(path, 'r+b') as f:
    
        line0 = f.readline()
        while line0 != '':        
            line = removeCommentFromLine(line0)[0]
            line = line.strip()
            tokens = line.split()
            tokens_upper = line.upper().split()
            
            if 'SPECIES' in line.upper():
                # Unread the line (we'll re-read it in readReactionBlock())
                f.seek(-len(line0), 1)
                readSpeciesBlock(f, speciesDict, speciesAliases, speciesList)
                
                # Also always add in a few bath gases (since RMG-Java does)
                for label, smiles in [('Ar','[Ar]'), ('He','[He]'), ('Ne','[Ne]'), ('N2','N#N')]:
                    molecule = Molecule().fromSMILES(smiles)
                    for species in speciesList:
                        if species.label == label:
                            if len(species.molecule) == 0:
                                species.molecule = [molecule]
                            break
                        if species.isIsomorphic(molecule):
                            break
                    else:
                        species = Species(label=label, molecule=[molecule])
                        speciesList.append(species)
                        speciesDict[label.upper()] = species                            
                
            elif 'THERM' in line.upper() and thermoPath is None:
                # Skip this if a thermo file is specified
                # Unread the line (we'll re-read it in readThermoBlock())
                f.seek(-len(line0), 1)
                readThermoBlock(f, speciesDict)
                
            elif 'REACTIONS' in line.upper():
                # Reactions section
                # Unread the line (we'll re-read it in readReactionBlock())
                f.seek(-len(line0), 1)
                reactionList = readReactionsBlock(f, speciesDict, readComments = readComments)
                    
            line0 = f.readline()
            
    # Read in the thermo data from the thermo file        
    if thermoPath:
        with open(thermoPath, 'r') as f:
            line0 = f.readline()
            while line0 != '':
                line = removeCommentFromLine(line0)[0]
                line = line.strip()
                if 'THERM' in line.upper():
                    f.seek(-len(line0), 1)
                    readThermoBlock(f, speciesDict)  
                    break
                line0 = f.readline()
    # Index the reactions now to have identical numbering as in Chemkin 
    index = 0
    for reaction in reactionList:
        index += 1
        reaction.index = index

    # Check for marked (and unmarked!) duplicate reactions
    # Combine marked duplicate reactions into a single reaction using MultiKinetics
    # Raise exception for unmarked duplicate reactions
    duplicateReactionsToRemove = []
    duplicateReactionsToAdd = []
    for index1 in range(len(reactionList)):
        reaction1 = reactionList[index1]
        if reaction1 in duplicateReactionsToRemove:
            continue

        for index2 in range(index1+1, len(reactionList)):
            reaction2 = reactionList[index2]
            if reaction1.reactants == reaction2.reactants and reaction1.products == reaction2.products:
                if reaction1.duplicate and reaction2.duplicate:
                    
                    if isinstance(reaction1, LibraryReaction) and isinstance(reaction2, LibraryReaction):
                        assert reaction1.library.label == reaction2.library.label
                        if reaction1 not in duplicateReactionsToRemove:
                            # already created duplicate reaction, move on to appending any additional duplicate kinetics
                            if isinstance(reaction1.kinetics,
                                          _kinetics.PDepArrhenius):
                                kinetics = _kinetics.MultiPDepArrhenius()
                            elif isinstance(reaction1.kinetics,
                                            _kinetics.Arrhenius):
                                kinetics = _kinetics.MultiArrhenius()
                            else:
                                logging.warning('Unexpected kinetics type {0} for duplicate reaction {1}. Not combining reactions.'.format(reaction1.kinetics.__class__, reaction1))
                                continue
                            reaction = LibraryReaction(
                                index = reaction1.index,
                                reactants = reaction1.reactants,
                                products = reaction1.products,
                                kinetics = kinetics,
                                library = reaction1.library,
                                duplicate = False,
                            )
                            duplicateReactionsToAdd.append(reaction)
                            kinetics.arrhenius = [reaction1.kinetics]
                            duplicateReactionsToRemove.append(reaction1)

                    else:
                        # Do not use as duplicate reactions if it's not a library reaction
                        # Template reactions should be kept separate
                        continue
                    
                    if (isinstance(reaction.kinetics,
                                   _kinetics.MultiPDepArrhenius) and
                            isinstance(reaction2.kinetics,
                                       _kinetics.PDepArrhenius)):
                        reaction.kinetics.arrhenius.append(reaction2.kinetics)
                    elif (isinstance(reaction.kinetics,
                                     _kinetics.MultiArrhenius) and
                              isinstance(reaction2.kinetics,
                                         _kinetics.Arrhenius)):
                        reaction.kinetics.arrhenius.append(reaction2.kinetics)
                    else:
                        raise ChemkinError('Mixed kinetics for duplicate reaction {0}.'.format(reaction))
                    
                    duplicateReactionsToRemove.append(reaction2)
                elif reaction1.kinetics.isPressureDependent() == reaction2.kinetics.isPressureDependent():
                    # If both reactions are pressure-independent or both are pressure-dependent, then they need duplicate tags
                    # Chemkin treates pdep and non-pdep reactions as different, so those are okay
                    raise ChemkinError('Encountered unmarked duplicate reaction {0}.'.format(reaction1))
                    
    for reaction in duplicateReactionsToRemove:
        reactionList.remove(reaction)
    reactionList.extend(duplicateReactionsToAdd)

    # If the transport path is given, then read it to obtain the transport
    # properties
    if transportPath:
        loadTransportFile(transportPath, speciesDict)
    
    # Apply species aliases if known
    for spec in speciesList:
        try:
            spec.label = speciesAliases[spec.label]
        except KeyError:
            pass
    
    # Attempt to extract index from species label
    indexPattern = re.compile(r'\(\d+\)$')
    for spec in speciesList:
        if indexPattern.search(spec.label):
            label, sep, index = spec.label[:-1].rpartition('(')
            spec.label = label
            spec.index = int(index)

    reactionList.sort(key=lambda reaction: reaction.index)
    return speciesList, reactionList


def readSpeciesBlock(f, speciesDict, speciesAliases, speciesList):
    """
    Read a Species block from a chemkin file.
    
    f is a file-like object that is just before the 'SPECIES' statement. When finished, it will have just passed the 'END' statement.
    speciesDict is a dictionary of species that will be updated.
    speciesAliases is a dictionary of species aliases that will be updated.
    speciesList is a list of species that will be extended.
    """
    line = f.readline()
    line = removeCommentFromLine(line)[0]
    line = line.strip()
    tokens = line.split()
    tokens_upper = line.upper().split() 
    firstToken = tokens.pop(0)
    firstToken = tokens_upper.pop(0) # pop from both lists
    assert firstToken in ['SPECIES', 'SPEC'] # should be first token in first line
    # Build list of species identifiers
    while 'END' not in tokens_upper:
        line = f.readline()
        # If the line contains only one species, and also contains
        # a comment with only one token, assume that token is 
        # intended to be the true identifier for the species, but
        # was not used e.g. due to a length limitation
        if '!' in line and len(line.split('!')) == 2:
            label, alias = line.split('!')
            label = label.strip()
            alias = alias.strip()
            if len(label.split()) == 1 and len(alias.split()) == 1:
                speciesAliases[label] = alias
        line = removeCommentFromLine(line)[0]
        line = line.strip()
        tokens.extend(line.split())
        tokens_upper.extend(line.upper().split())
    # Now process list of tokens
    processed_tokens = []
    for token in tokens:
        if token in processed_tokens:
            continue # ignore species declared twice
        token_upper = token.upper()
        if token_upper in ['SPECIES', 'SPEC']:
            continue # there may be more than one SPECIES statement
        if token_upper == 'END':
            break
        processed_tokens.append(token)
        if token in speciesDict:
            logging.debug("Re-using species {0} already in speciesDict".format(token))
            species = speciesDict[token]
        else:
            species = Species(label=token)
            speciesDict[token] = species
        speciesList.append(species)
        
def readThermoBlock(f, speciesDict):
    """
    Read a thermochemistry block from a chemkin file.
    
    f is a file-like object that is just before the 'THERM' statement.
    When finished, it will have just passed the 'END' statement.
    speciesDict is a dictionary of species that will be updated with the given thermodynamics.
    
    Returns a dictionary of molecular formulae for each species, in the form
    `{'methane': {'C':1, 'H':4}}
    
    If duplicate entries are found, the FIRST is used, and a warning is printed.
    """
    # List of thermodynamics (hopefully one per species!)
    formulaDict = {}
    line = f.readline()
    assert line.upper().strip().startswith('THER'), "'{0}' doesn't begin with THERM statement.".format(line)
    line = f.readline()
    
    # In case there are commented lines immediately after THER
    meaningfulline, comment = removeCommentFromLine(line)
    while not meaningfulline.strip():
        line = f.readline()
        meaningfulline, comment = removeCommentFromLine(line)
    Tmin = Tint = Tmax = None
    try:
        Tmin = float(meaningfulline[0:9].strip())
        Tint = float(meaningfulline[10:19].strip())
        Tmax = float(meaningfulline[20:29].strip())
        if [Tmin, Tint, Tmax] != [float(i) for i in meaningfulline.split()[0:3]]:
            logging.warning("Default temperature range line {0!r} may be badly formatted.".format(line))
            logging.warning("It should have Tmin in columns 1-10, Tmid in columns 11-20, and Tmax in columns 21-30")
        logging.info("Thermo file has default temperature range {0} to {1} and {1} to {2}".format(Tmin, Tint, Tmax))
        line = f.readline()
    except:
        logging.info("Thermo file has no default temperature ranges")
        logging.info("(The line it would be on is {0!r} but that is not formatted as such)".format(line))
        logging.info("(It should have Tmin in columns 1-10, Tmid in columns 11-20, and Tmax in columns 21-30)")
    thermoBlock = ''
    comments = ''
    while line != '' and not line.upper().strip().startswith('END'):
        line, comment = removeCommentFromLine(line)
        if comment: comments += comment.strip().replace('\t',', ') + '\n'

        if len(line) < 80:
            if line.strip():
                logging.info("Ignoring short but non-empty line: {0!r}".format(line))
            line = f.readline()
            continue

        if line[79] not in ['1', '2', '3', '4']:
            logging.warning("Ignoring line without 1,2,3 or 4 in 80th column: {0!r}".format(line))
            line = f.readline()
            continue

        thermoBlock += line
        if line[79] == '4':
            try:
                label, thermo, formula = readThermoEntry(thermoBlock, Tmin=Tmin, Tint=Tint, Tmax=Tmax)
            except:
                logging.error("Error reading thermo block:\n" + thermoBlock)
                raise
            if label not in speciesDict:
                logging.info("Ignoring thermo data for {0} because it's not in the requested list of species.".format(label))
                thermoBlock = ''
                line = f.readline()
                continue
            elif speciesDict[label].thermo:
                logging.warning('Skipping duplicate thermo for the species {0}'.format(label))
                thermoBlock = ''
                line = f.readline()
                continue
            else:
                if thermo is None:
                    logging.error("Problematic thermo block:\n{0}".format(thermoBlock))
                    raise ChemkinError('Error while reading thermo entry for required species {0}'.format(label))
            try:
                formulaDict[label] = formula
                speciesDict[label].thermo = thermo
                speciesDict[label].thermo.comment = getattr(speciesDict[label].thermo,'comment','') 
                if comments:
                    speciesDict[label].thermo.comment += '\n{0}'.format(comments)
                comments = ''
            except KeyError:
                if label.upper() in ['AR', 'N2', 'HE', 'NE']:
                    logging.warning('Skipping species"{0}" while reading thermodynamics entry.'.format(label))
                else:
                    logging.warning('Skipping unexpected species "{0}" while reading thermodynamics entry.'.format(label))
            thermoBlock = ''
        assert len(thermoBlock.split('/n'))<=4, "Should only have 4 lines in a thermo block:\n{0}".format(thermoBlock)
        line = f.readline()
    return formulaDict

        
def readReactionsBlock(f, speciesDict, readComments = True):
    """
    Read a reactions block from a Chemkin file stream.
    
    This function can also read the ``reactions.txt`` and ``pdepreactions.txt``
    files from RMG-Java kinetics libraries, which have a similar syntax.
    """    
    energyUnits = 'cal/mol'
    moleculeUnits = 'moles'
    volumeUnits = 'cm3'
    timeUnits = 's'
    
    line = f.readline()
    found = False
    while line != '' and not found:
    
        line = removeCommentFromLine(line)[0]
        line = line.strip()
        tokens = line.split()
        
        if len(tokens) > 0 and tokens[0].upper() == 'REACTIONS':
            # Regular Chemkin file
            found = True
            for token in tokens[1:]:
                unit = token.lower()
                if unit in ['molecules', 'moles', 'mole', 'mol', 'molecule']:
                    moleculeUnits = unit
                elif unit in ['kcal/mole', 'kcal/mol', 'cal/mole', 'cal/mol', 'kj/mole', 'kj/mol', 'j/mole', 'j/mol', 'kelvins']:
                    energyUnits = unit
                else:
                    raise ChemkinError('Unknown unit type "{0}"'.format(unit))

        elif len(tokens) > 0 and tokens[0].lower() == 'unit:':
            # RMG-Java kinetics library file
            found = True
            while 'reactions:' not in line.lower():
                line = f.readline()
                line = removeCommentFromLine(line)[0]
                line = line.strip()
                
                if 'A:' in line or 'E:' in line:
                    units = line.split()[1]
                    if 'A:' in line:
                        moleculeUnits, volumeUnits, timeUnits = units.lower().split('/') # Assume this is a 3-tuple: moles or molecules, volume, time
                    elif 'E:' in line:
                        energyUnits = units.lower()
        else:
            line = f.readline()
            
    if not found:
        raise ChemkinError('Invalid reaction block.')
    
    # Check that the units are valid
    assert moleculeUnits in ['molecules', 'moles', 'mole', 'mol', 'molecule']
    assert volumeUnits in ['cm3', 'm3']
    assert timeUnits in ['s']
    assert energyUnits in ['kcal/mole', 'kcal/mol', 'cal/mole', 'cal/mol', 'kj/mole', 'kj/mol', 'j/mole', 'j/mol', 'kelvins']
    
    # Homogenize units
    if moleculeUnits == 'molecules':
        moleculeUnits = 'molecule'
    elif moleculeUnits == 'moles' or moleculeUnits == 'mole':
        moleculeUnits = 'mol'
    volumeUnits = {'cm3': 'cm', 'm3': 'm'}[volumeUnits]
    if energyUnits == 'kcal/mole':
        energyUnits = 'kcal/mol'
    elif energyUnits == 'cal/mole':
        energyUnits = 'cal/mol'
    elif energyUnits == 'kj/mole':
        energyUnits = 'kj/mol'
    elif energyUnits == 'j/mole':
        energyUnits = 'j/mol'
    elif energyUnits == 'kelvins':
        energyUnits = 'K'
    energyUnits = energyUnits.replace('j/mol', 'J/mol')
    
    # Set up kinetics units
    Aunits = [
        '',                                                                 # Zeroth-order
        's^-1'.format(timeUnits),                                           # First-order
        '{0}^3/({1}*{2})'.format(volumeUnits, moleculeUnits, timeUnits),    # Second-order
        '{0}^6/({1}^2*{2})'.format(volumeUnits, moleculeUnits, timeUnits),  # Third-order
        '{0}^9/({1}^3*{2})'.format(volumeUnits, moleculeUnits, timeUnits),  # Fourth-order
    ]
    Eunits = energyUnits
    
    kineticsList = []
    commentsList = []
    kinetics = ''
    comments = ''
    
    line = f.readline()
    while line != '':

        lineStartsWithComment = line.lstrip().startswith('!') or line.lstrip().startswith('//')
        line, comment = removeCommentFromLine(line)
        line = line.strip(); comment = comment.strip()
    
        if 'end' in line or 'END' in line:
            break

        # if 'rev' in line or 'REV' in line:
            # can no longer name reactants rev...
        #    line = f.readline()
        #    continue  # need to re-do the comment stripping!

        if '=' in line and not lineStartsWithComment:
            # Finish previous record
            kineticsList.append(kinetics)
            commentsList.append(comments)
            kinetics = ''
            comments = ''
            
        if line: kinetics += line + '\n'
        if comment: comments += comment + '\n'
        
        line = f.readline()
        
    # Don't forget the last reaction!
    if kinetics.strip() != '':
        kineticsList.append(kinetics)
        commentsList.append(comments)
    
    if len(kineticsList) == 0 and len(commentsList) == 0:
        # No reactions found
        pass
    elif kineticsList[0] == '' and commentsList[-1] == '':
        # True for Chemkin files generated from RMG-Py
        kineticsList.pop(0)
        commentsList.pop(-1)
    elif kineticsList[0] == '' and commentsList[0] == '':
        # True for Chemkin files generated from RMG-Java
        kineticsList.pop(0)
        commentsList.pop(0)
    else:
        # In reality, comments can occur anywhere in the Chemkin
        # file (e.g. either or both of before and after the
        # reaction equation)
        # If we can't tell what semantics we are using, then just
        # throw the comments away
        # (This is better than failing to load the Chemkin file at
        # all, which would likely occur otherwise)
        if kineticsList[0] == '':
            kineticsList.pop(0)
        if len(kineticsList) != len(commentsList):
            logging.warning("Discarding comments from Chemkin file because not sure which reaction they apply to")
            commentsList = ['' for kinetics in kineticsList]
        
    reactionList = []
    for kinetics, comments in zip(kineticsList, commentsList):
        try:
            reaction = readKineticsEntry(kinetics, speciesDict, Aunits, Eunits)
            reaction = readReactionComments(reaction, comments, read = readComments)
        except ChemkinError, e:
            if e.message == "Skip reaction!":
                logging.warning("Skipping the reaction {0!r}".format(kinetics))
                continue
            else:
                raise e
        reactionList.append(reaction)
        
    return reactionList

################################################################################

def saveHTMLFile(path, readComments = True):
    """
    Save an output HTML file from the contents of a RMG-Java output folder
    """
    from rmgpy.rmg.model import CoreEdgeReactionModel
    from rmgpy.rmg.output import saveOutputHTML
    chemkinPath= path + '/chemkin/chem.inp'
    dictionaryPath = path + 'RMG_Dictionary.txt'
    model = CoreEdgeReactionModel()
    model.core.species, model.core.reactions = loadChemkinFile(chemkinPath,dictionaryPath, readComments = readComments)
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
    label = species.label
    # Special case for inert colliders - just use the label if possible
    if not species.reactive and 0 < len(label) <= 10:
        return label

    # The algorithm is slightly different depending on whether or not the
    # species has an index
    # If so, we want to include the index in the identifier
    if species.index == -1:
        # No index present -- probably not in RMG job
        # In this case just return the label (if the right size)
        if len(label) > 0 and not re.search(r'[^A-Za-z0-9\-_,\(\)\*#]+', label):
            if len(label) <= 10:
                return label
            elif len(label) <= 15:
                #logging.warning('Species label {0} is longer than 10 characters and may exceed chemkin string limit'.format(label))
                return label            
            else:
                logging.warning('Species label is longer than 15 characters and will break CHEMKIN 2.0')
                return label
        else:
            # try the chemical formula if the species label is not present
            if len(species.molecule) > 0:
                # Try the chemical formula
                return '{0}'.format(species.molecule[0].getFormula())
    else:
        
        # Index present - the index will be included in the identifier
        # (at the expense of the current label or formula if need be)

        # First try to use the label and index
        # The label can only contain alphanumeric characters, and -()*#_,
        if len(label) > 0 and species.index >= 0 and not re.search(r'[^A-Za-z0-9\-_,\(\)\*#]+', label):
            name = '{0}({1:d})'.format(label, species.index)
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

def writeThermoEntry(species, verbose = True):
    """
    Return a string representation of the NASA model readable by Chemkin.
    To use this method you must have exactly two NASA polynomials in your
    model, and you must use the seven-coefficient forms for each.
    """

    thermo = species.thermo
    if not isinstance(thermo, NASA):
        return ''
        raise ChemkinError('Cannot generate Chemkin string for species "{0}": Thermodynamics data must be a NASA object.'.format(species))

    assert len(thermo.polynomials) == 2
    assert thermo.polynomials[0].Tmin.value_si < thermo.polynomials[1].Tmin.value_si
    assert thermo.polynomials[0].Tmax.value_si == thermo.polynomials[1].Tmin.value_si
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
    # Remove elements with zero count
    index = 0
    while index < len(elementCounts):
        if elementCounts[index] == 0:
            del elements[index]
            del elementCounts[index]
        else:
            index += 1
    
    string = ''
    # Write thermo comments
    if verbose:
        if thermo.comment:
            for line in thermo.comment.split("\n"):
                if len(line) > 150:
                    short_lines = textwrap.fill(line,150).split("\n")
                    for short_line in short_lines:
                        string += "! {0}\n".format(short_line) 
                else:
                    string += "! {0}\n".format(line) 

    # Line 1
    string += '{0:<16}        '.format(getSpeciesIdentifier(species))
    if len(elements) <= 4:
        # Use the original Chemkin syntax for the element counts
        for symbol, count in zip(elements, elementCounts):
            string += '{0!s:<2}{1:<3d}'.format(symbol, count)
        string += '     ' * (4 - len(elements))
    else:
        string += '     ' * 4
    string += 'G{0:<10.3f}{1:<10.3f}{2:<8.2f}      1'.format(thermo.polynomials[0].Tmin.value_si, thermo.polynomials[1].Tmax.value_si, thermo.polynomials[0].Tmax.value_si)
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

def writeReactionString(reaction, javaLibrary = False):
    """
    Return a reaction string in chemkin format.
    """
    kinetics = reaction.kinetics
    
    if kinetics is None:
        reaction_string = ' + '.join([getSpeciesIdentifier(reactant) for reactant in reaction.reactants])
        reaction_string += ' => ' if not reaction.reversible else ' = '
        reaction_string += ' + '.join([getSpeciesIdentifier(product) for product in reaction.products])
        return reaction_string
    
    if javaLibrary:
        thirdBody = ''
        if kinetics.isPressureDependent():
            if (isinstance(kinetics, _kinetics.ThirdBody) and
                    not isinstance(kinetics, _kinetics.Lindemann) and
                    not isinstance(kinetics, _kinetics.Troe)):
                thirdBody = ' + M'
            elif isinstance(kinetics, _kinetics.PDepArrhenius):
                thirdBody = ''
            elif isinstance(kinetics, _kinetics.Chebyshev):
                thirdBody = ''
            else:
                thirdBody = ' (+M)'
        
        reaction_string = ' + '.join([getSpeciesIdentifier(reactant) for reactant in reaction.reactants])
        reaction_string += thirdBody
        reaction_string += ' => ' if not reaction.reversible else ' = '
        reaction_string += ' + '.join([getSpeciesIdentifier(product) for product in reaction.products])
        reaction_string += thirdBody
    
    else:
        thirdBody = ''
        if kinetics.isPressureDependent():
            if (isinstance(kinetics, _kinetics.ThirdBody) and
                    not isinstance(kinetics,
                                   (_kinetics.Lindemann, _kinetics.Troe))):
                thirdBody = '+M'
            elif isinstance(kinetics, _kinetics.PDepArrhenius):
                thirdBody = ''
            else:
                thirdBody = '(+M)'
        
        reaction_string = '+'.join([getSpeciesIdentifier(reactant) for reactant in reaction.reactants])
        reaction_string += thirdBody
        reaction_string += '=>' if not reaction.reversible else '='
        reaction_string += '+'.join([getSpeciesIdentifier(product) for product in reaction.products])
        reaction_string += thirdBody

    if len(reaction_string) > 52:
        logging.warning("Chemkin reaction string {0!r} is too long for Chemkin 2!".format(reaction_string))
    return reaction_string

################################################################################

def writeTransportEntry(species, verbose = True):
    """
    Return a string representation of the reaction as used in a Chemkin file. Lists the 
    """
    
################################################################################

def writeKineticsEntry(reaction, speciesList, verbose = True, javaLibrary = False):
    """
    Return a string representation of the reaction as used in a Chemkin
    file. Use verbose = True to turn on comments.  Use javaLibrary = True in order to 
    generate a kinetics entry suitable for an RMG-Java kinetics library.  
    """
    string = ""
    
    if isinstance(reaction.kinetics,
                  (_kinetics.MultiArrhenius, _kinetics.MultiPDepArrhenius)):
        if verbose:
            if reaction.kinetics.comment:
                for line in reaction.kinetics.comment.split("\n"):
                    string += "! {0}\n".format(line) 
        for kinetics in reaction.kinetics.arrhenius:
            if isinstance(reaction,LibraryReaction):
                new_reaction = LibraryReaction( index=reaction.index,
                     reactants=reaction.reactants,
                     products=reaction.products,
                     reversible=reaction.reversible,
                     kinetics=kinetics,
                     library=reaction.library
                     )
            else:
                new_reaction = Reaction( index=reaction.index,
                         reactants=reaction.reactants,
                         products=reaction.products,
                         reversible=reaction.reversible,
                         kinetics=kinetics)
            string += writeKineticsEntry(new_reaction, speciesList, verbose, javaLibrary)
            string += "DUPLICATE\n"
        return string + "\n"
    
    # Add to global chemkin reaction count if the kinetics is not a duplicate
    global __chemkin_reaction_count
    if __chemkin_reaction_count is not None:
        __chemkin_reaction_count += 1
            
    if verbose:        
        # Next line of comment contains Chemkin and RMG indices
        if __chemkin_reaction_count is not None:
            string += "! Reaction index: Chemkin #{0:d}; RMG #{1:d}\n".format(__chemkin_reaction_count, reaction.index)
        
        # Next line of comment contains information about the type of reaction
        if isinstance(reaction, TemplateReaction):
            string += '! Template reaction: {0!s}\n'.format(reaction.family.label)
        elif isinstance(reaction, LibraryReaction):
            string += '! Library reaction: {0!s}\n'.format(reaction.library.label)
        elif isinstance(reaction, PDepReaction):
            string += '! PDep reaction: {0!s}\n'.format(reaction.network)          
            if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
                # Print additional information about the pdep network's high-P limit reactions if in debug mode.
                for rxn in reaction.network.pathReactions:
                    if isinstance(rxn, LibraryReaction):
                        string += '! High-P limit: {0} (Library reaction: {1!s})\n'.format(rxn, rxn.library.label)
                    else:
                        string += '! High-P limit: {0} (Template reaction: {1!s})\n'.format(rxn, rxn.family.label)   
    
        # Remaining lines of comments taken from reaction kinetics
        if reaction.kinetics.comment:
            for line in reaction.kinetics.comment.split("\n"):
                if len(line) > 150:
                    short_lines = textwrap.fill(line,150).split("\n")
                    for short_line in short_lines:
                        string += "! {0}\n".format(short_line) 
                else:
                    string += "! {0}\n".format(line)                               
    
    kinetics = reaction.kinetics
    numReactants = len(reaction.reactants)
    reaction_string = writeReactionString(reaction, javaLibrary)    
    
    string += '{0!s:<51} '.format(reaction_string)

    if isinstance(kinetics, _kinetics.Arrhenius):
        string += '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(
            kinetics.A.value_si/ (kinetics.T0.value_si ** kinetics.n.value_si) * 1.0e6 ** (numReactants - 1),
            kinetics.n.value_si,
            kinetics.Ea.value_si / 4184.
        )
    elif isinstance(kinetics, (_kinetics.Lindemann, _kinetics.Troe)):
        arrhenius = kinetics.arrheniusHigh
        string += '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(
            arrhenius.A.value_si / (arrhenius.T0.value_si ** arrhenius.n.value_si) * 1.0e6 ** (numReactants - 1),
            arrhenius.n.value_si,
            arrhenius.Ea.value_si / 4184.
        )
    elif isinstance(kinetics, _kinetics.ThirdBody):
        arrhenius = kinetics.arrheniusLow
        string += '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(
            arrhenius.A.value_si / (arrhenius.T0.value_si ** arrhenius.n.value_si) * 1.0e6 ** (numReactants),
            arrhenius.n.value_si,
            arrhenius.Ea.value_si / 4184.
        )
    elif hasattr(kinetics,'highPlimit') and kinetics.highPlimit is not None:
        arrhenius = kinetics.highPlimit
        string += '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(
            arrhenius.A.value_si / (arrhenius.T0.value_si ** arrhenius.n.value_si) * 1.0e6 ** (numReactants - 1),
            arrhenius.n.value_si,
            arrhenius.Ea.value_si / 4184.
            )
    else:
        # Print dummy values that Chemkin parses but ignores
        string += '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(1, 0, 0)
        
    if javaLibrary:
        # Assume uncertainties are zero (when parsing from chemkin), may need to adapt later
        string += '{0:<9.1f} {1:<9.1f} {2:<9.1f}'.format(0, 0, 0)

    string += '\n'

    if isinstance(kinetics,
                  (_kinetics.ThirdBody, _kinetics.Lindemann, _kinetics.Troe)):
        # Write collider efficiencies
        for collider, efficiency in sorted(kinetics.efficiencies.items(), key=lambda item: id(item[0])):
            for species in speciesList:
                if any([collider.isIsomorphic(molecule) for molecule in species.molecule]):
                    string += '{0!s}/{1:<4.2f}/ '.format(getSpeciesIdentifier(species), efficiency)
                    break
        string += '\n'
        
        if isinstance(kinetics, (_kinetics.Lindemann, _kinetics.Troe)):
            # Write low-P kinetics
            arrhenius = kinetics.arrheniusLow
            string += '    LOW/ {0:<9.3e} {1:<9.3f} {2:<9.3f}/\n'.format(
                arrhenius.A.value_si / (arrhenius.T0.value_si ** arrhenius.n.value_si) * 1.0e6 ** (numReactants),
                arrhenius.n.value_si,
                arrhenius.Ea.value_si / 4184.
            )
            if isinstance(kinetics, _kinetics.Troe):
                # Write Troe parameters
                if kinetics.T2 is None:
                    string += '    TROE/ {0:<9.3e} {1:<9.3g} {2:<9.3g}/\n'.format(kinetics.alpha, kinetics.T3.value_si, kinetics.T1.value_si)
                else:
                    string += '    TROE/ {0:<9.3e} {1:<9.3g} {2:<9.3g} {3:<9.3g}/\n'.format(kinetics.alpha, kinetics.T3.value_si, kinetics.T1.value_si, kinetics.T2.value_si)
    elif isinstance(kinetics, _kinetics.PDepArrhenius):
        for P, arrhenius in zip(kinetics.pressures.value_si, kinetics.arrhenius):
            if isinstance(arrhenius, _kinetics.MultiArrhenius):
                for arrh in arrhenius.arrhenius:
                    string += '    PLOG/ {0:<9.3f} {1:<9.3e} {2:<9.3f} {3:<9.3f}/\n'.format(P / 101325.,
                    arrh.A.value_si / (arrh.T0.value_si ** arrh.n.value_si) * 1.0e6 ** (numReactants - 1),
                    arrh.n.value_si,
                    arrh.Ea.value_si / 4184.
                    )
            else:
                string += '    PLOG/ {0:<9.3f} {1:<9.3e} {2:<9.3f} {3:<9.3f}/\n'.format(P / 101325.,
                    arrhenius.A.value_si / (arrhenius.T0.value_si ** arrhenius.n.value_si) * 1.0e6 ** (numReactants - 1),
                    arrhenius.n.value_si,
                    arrhenius.Ea.value_si / 4184.
                )
    elif isinstance(kinetics, _kinetics.Chebyshev):
        string += '    TCHEB/ {0:<9.3f} {1:<9.3f}/\n'.format(kinetics.Tmin.value_si, kinetics.Tmax.value_si)
        string += '    PCHEB/ {0:<9.3f} {1:<9.3f}/\n'.format(kinetics.Pmin.value_si / 101325., kinetics.Pmax.value_si / 101325.)
        string += '    CHEB/ {0:d} {1:d}/\n'.format(kinetics.degreeT, kinetics.degreeP)
        if kinetics.degreeP < 6:
            coeffs = kinetics.coeffs.value_si.copy()
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
                    coeffs.append(kinetics.coeffs.value_si[i,j])
            coeffs[0] += 6 * (numReactants - 1)
            for i in range(len(coeffs)):
                if i % 5 == 0: string += '    CHEB/'
                string += ' {0:<12.3e}'.format(coeffs[i])
                if i % 5 == 4: string += '/\n'

    if reaction.duplicate:
        string += 'DUPLICATE\n'

    return string

################################################################################

def markDuplicateReaction(test_reaction, reaction_list):
    """
    If the test_reaction is a duplicate (in Chemkin terms) of one in reaction_list, then set `duplicate=True` on both instances.
    `reaction_list` can be any iterator.
    It does not add the testReaction to the reactionList - you probably want to do this yourself afterwards.
    """
    reaction1 = test_reaction
    for reaction2 in reaction_list:
        if reaction1.__class__ != reaction2.__class__:
            # TemplateReaction, LibraryReaction, and PDepReaction cannot be
            # duplicates of one another.
            # RHW question: why can't TemplateReaction be duplicate of LibraryReaction, in Chemkin terms? I guess it shouldn't happen in RMG.
            continue
        if (reaction1.reactants == reaction2.reactants and reaction1.products == reaction2.products) \
        or (reaction1.products == reaction2.reactants and reaction1.reactants == reaction2.products):
            if reaction1.duplicate and reaction2.duplicate:                
                if reaction1.kinetics.isPressureDependent() != reaction2.kinetics.isPressureDependent():
                    logging.warning('Marked reaction {0} as not duplicate because of mixed pressure dependence for saving to Chemkin file.'.format(reaction1))
                    reaction1.duplicate = False
                    reaction2.duplicate = False
            else:
                if reaction1.kinetics.isPressureDependent() == reaction2.kinetics.isPressureDependent():
                    # Only mark as duplicate if both reactions are pressure dependent or both are
                    # not pressure dependent.  Do not mark as duplicates otherwise.
                    logging.warning('Marked reaction {0} as duplicate for saving to Chemkin file.'.format(reaction1))
                    reaction1.duplicate = True
                    reaction2.duplicate = True

def markDuplicateReactions(reactions):
    """
    For a given list of `reactions`, mark all of the duplicate reactions as
    understood by Chemkin.
    
    This is pretty slow (quadratic in size of reactions list) so only call it if you're really worried
    you may have undetected duplicate reactions.
    """
    for index1 in range(len(reactions)):
        reaction1 = reactions[index1]
        remainingList = reactions[index1+1:]
        markDuplicateReaction(reaction1, remainingList)
 

def saveSpeciesDictionary(path, species, oldStyle=False):
    """
    Save the given list of `species` as adjacency lists in a text file `path` 
    on disk.
    
    If `oldStyle==True` then it saves it in the old RMG-Java syntax.
    """
    with open(path, 'w') as f:
        for spec in species:
            if oldStyle:
                try:
                    f.write(spec.molecule[0].toAdjacencyList(label=getSpeciesIdentifier(spec), removeH=True, oldStyle=True))
                except:
                    newAdjList = spec.molecule[0].toAdjacencyList(label=getSpeciesIdentifier(spec), removeH=False)
                    f.write("// Couldn't save {0} in old RMG-Java syntax, but here it is in newer RMG-Py syntax:".format(getSpeciesIdentifier(spec)))
                    f.write("\n// " + "\n// ".join(newAdjList.splitlines()) + '\n')
            else:
                try:
                    f.write(spec.molecule[0].toAdjacencyList(label=getSpeciesIdentifier(spec), removeH=False))
                except:
                    raise ChemkinError('Ran into error saving dictionary for species {0}. Please check your files.'.format(getSpeciesIdentifier(spec)))
            f.write('\n')

def saveTransportFile(path, species):
    r"""
    Save a Chemkin transport properties file to `path` on disk containing the
    transport properties of the given list of `species`.
    
    The syntax is from the Chemkin TRANSPORT manual.
    The first 16 columns in each line of the database are reserved for the species name
    (Presently CHEMKIN is programmed to allow no more than 16-character names.) 
    Columns 17 through 80 are free-format, and they contain the molecular parameters for each species. They are, in order:
    
    1. An index indicating whether the molecule has a monatomic, linear or nonlinear geometrical configuration.
       If the index is 0, the molecule is a single atom. 
       If the index is 1 the molecule is linear, and 
       if it is 2, the molecule is nonlinear.
    2. The Lennard-Jones potential well depth  :math:`\epsilon / k_B` in Kelvins.
    3. The Lennard-Jones collision diameter :math:`\sigma` in Angstroms.
    4. The dipole moment :math:`\mu` in Debye. Note: a Debye is :math:`10^{-18} cm^{3/2}erg^{1/2}`.
    5. The polarizability :math:`\alpha` in cubic Angstroms.
    6. The rotational relaxation collision number :math:`Z_rot` at 298K.
    7. After the last number, a comment field can be enclosed in parenthesis.

    """
    with open(path, 'w') as f:
        f.write("! {0:15} {1:8} {2:9} {3:9} {4:9} {5:9} {6:9} {7:9}\n".format('Species','Shape', 'LJ-depth', 'LJ-diam', 'DiplMom', 'Polzblty', 'RotRelaxNum','Data'))
        f.write("! {0:15} {1:8} {2:9} {3:9} {4:9} {5:9} {6:9} {7:9}\n".format('Name','Index', 'epsilon/k_B', 'sigma', 'mu', 'alpha', 'Zrot','Source'))
        for spec in species:            
            if (not spec.transportData or
                len(spec.molecule) == 0):
                missingData = True
            else:
                missingData = False
            
            label = getSpeciesIdentifier(spec)
            
            molecule = spec.molecule[0]
            if len(molecule.atoms) == 1:
                shapeIndex = 0
            elif molecule.isLinear():
                shapeIndex = 1
            else:
                shapeIndex = 2
            
            if missingData:
                f.write('! {0:19s} {1!r}\n'.format(label, spec.transportData))
            else:
                f.write('{0:19} {1:d}   {2:9.3f} {3:9.3f} {4:9.3f} {5:9.3f} {6:9.3f}    ! {7:s}\n'.format(
                    label,
                    shapeIndex,
                    spec.transportData.epsilon.value_si / constants.R,
                    spec.transportData.sigma.value_si * 1e10,
                    (spec.transportData.dipoleMoment.value_si * constants.c * 1e21 if spec.transportData.dipoleMoment else 0),
                    (spec.transportData.polarizability.value_si * 1e30 if spec.transportData.polarizability else 0),
                    (spec.Zrot.value_si if spec.Zrot else 0),
                    spec.transportData.comment,
                ))

def saveChemkinFile(path, species, reactions, verbose = True, checkForDuplicates=True):
    """
    Save a Chemkin input file to `path` on disk containing the provided lists
    of `species` and `reactions`.
    If checkForDuplicates is False then we don't check for unlabeled duplicate reactions,
    thus saving time (eg. if you are sure you've already labeled them as duplicate).
    """
    # Check for duplicate
    if checkForDuplicates:
        markDuplicateReactions(reactions)
    
    f = open(path, 'w')
    
    sorted_species = sorted(species, key=lambda species: species.index)

    # Elements section
    f.write('ELEMENTS H C O N Ne Ar He Si S Cl END\n\n')

    # Species section
    f.write('SPECIES\n')
    for spec in sorted_species:
        label = getSpeciesIdentifier(spec)
        if verbose:
            f.write('    {0!s:<16}    ! {1}\n'.format(label, str(spec)))
        else:
            f.write('    {0!s:<16}\n'.format(label))
    f.write('END\n\n\n\n')

    # Thermodynamics section
    f.write('THERM ALL\n')
    f.write('    300.000  1000.000  5000.000\n\n')
    for spec in sorted_species:
        f.write(writeThermoEntry(spec, verbose=verbose))
        f.write('\n')
    f.write('END\n\n\n\n')

    ## Transport section would go here
    #f.write('TRANSPORT\n')
    #for spec in sorted_species:
        #f.write(writeTransportEntry(spec)
        #f.write('\n')
    #f.write('END\n\n')

    # Reactions section
    f.write('REACTIONS    KCAL/MOLE   MOLES\n\n')
    global __chemkin_reaction_count
    __chemkin_reaction_count = 0
    for rxn in reactions:
        f.write(writeKineticsEntry(rxn, speciesList=species, verbose=verbose))
        # Don't forget to mark duplicates!
        f.write('\n')
    f.write('END\n\n')
    f.close()
    logging.info("Chemkin file contains {0} reactions.".format(__chemkin_reaction_count))
    __chemkin_reaction_count = None

def saveJavaKineticsLibrary(path, species, reactions):
    """
    Save the reaction files for a RMG-Java kinetics library: pdepreactions.txt
    and reactions.txt given a list of reactions, with species.txt containing the
    RMG-Java formatted dictionary.
    """
    # Check for duplicate
    markDuplicateReactions(reactions)
    
    f = open(os.path.join(os.path.dirname(path), 'reactions.txt'), 'w')
    f2 = open(os.path.join(os.path.dirname(path), 'pdepreactions.txt'), 'w')

    # Headers
    f.write('Unit:\n')
    f.write('A: mol/cm3/s\n')
    f.write('E: kcal/mol\n')
    f.write('\n')
    f.write('Reactions:\n')
    f.write('\n')
    
    f2.write('Unit:\n')
    f2.write('A: mol/cm3/s\n')
    f2.write('E: kcal/mol\n')
    f2.write('\n')
    f2.write('Reactions:\n')
    f2.write('\n')
    

    for rxn in reactions:
        if rxn.kinetics.isPressureDependent():
            f2.write(writeKineticsEntry(rxn, speciesList=species, verbose = False, javaLibrary = True))
            f2.write('\n')
        else:  
            f.write(writeKineticsEntry(rxn, speciesList=species, verbose = False, javaLibrary = True))
            f.write('\n')
    f.close()
    f2.close()
    
    saveSpeciesDictionary(os.path.join(os.path.dirname(path), 'species.txt'), species, oldStyle=True)

def saveChemkin(reactionModel, path, verbose_path, dictionaryPath=None, transportPath=None, saveEdgeSpecies=False):
    """
    Save a Chemkin file for the current model as well as any desired output
    species and reactions to `path`. If `saveEdgeSpecies` is True, then 
    a chemkin file and dictionary file for the core and edge species and reactions
    will be saved.  
    """
    
    if saveEdgeSpecies == False:
        speciesList = reactionModel.core.species + reactionModel.outputSpeciesList
        rxnList = reactionModel.core.reactions + reactionModel.outputReactionList
        saveChemkinFile(path, speciesList, rxnList, verbose = False, checkForDuplicates=False) # We should already have marked everything as duplicates by now        
        logging.info('Saving current model to verbose Chemkin file...')
        saveChemkinFile(verbose_path, speciesList, rxnList, verbose = True, checkForDuplicates=False)
        if dictionaryPath:
            saveSpeciesDictionary(dictionaryPath, speciesList)
        if transportPath:
            saveTransportFile(transportPath, speciesList)
        
    else:
        speciesList = reactionModel.core.species + reactionModel.edge.species + reactionModel.outputSpeciesList
        rxnList = reactionModel.core.reactions + reactionModel.edge.reactions + reactionModel.outputReactionList
        saveChemkinFile(path, speciesList, rxnList, verbose = False, checkForDuplicates=False)        
        logging.info('Saving current core and edge to verbose Chemkin file...')
        saveChemkinFile(verbose_path, speciesList, rxnList, verbose = True, checkForDuplicates=False)
        if dictionaryPath:
            saveSpeciesDictionary(dictionaryPath, speciesList)
        if transportPath:
            saveTransportFile(transportPath, speciesList)

def saveChemkinFiles(rmg):
        """
        Save the current reaction model to a set of Chemkin files.
        """        
        logging.info('Saving current model core to Chemkin file...')
        this_chemkin_path = os.path.join(rmg.outputDirectory, 'chemkin', 'chem{0:04d}.inp'.format(len(rmg.reactionModel.core.species)))
        latest_chemkin_path = os.path.join(rmg.outputDirectory, 'chemkin','chem.inp')
        latest_chemkin_verbose_path = os.path.join(rmg.outputDirectory, 'chemkin', 'chem_annotated.inp')
        latest_dictionary_path = os.path.join(rmg.outputDirectory, 'chemkin','species_dictionary.txt')
        latest_transport_path = os.path.join(rmg.outputDirectory, 'chemkin', 'tran.dat')
        saveChemkin(rmg.reactionModel, this_chemkin_path, latest_chemkin_verbose_path, latest_dictionary_path, latest_transport_path, False)
        if os.path.exists(latest_chemkin_path):
            os.unlink(latest_chemkin_path)
        shutil.copy2(this_chemkin_path,latest_chemkin_path)
        
        if rmg.saveEdgeSpecies == True:
            logging.info('Saving current model core and edge to Chemkin file...')
            this_chemkin_path = os.path.join(rmg.outputDirectory, 'chemkin', 'chem_edge%04i.inp' % len(rmg.reactionModel.core.species)) # len() needs to be core to have unambiguous index
            latest_chemkin_path = os.path.join(rmg.outputDirectory, 'chemkin','chem_edge.inp')
            latest_chemkin_verbose_path = os.path.join(rmg.outputDirectory, 'chemkin', 'chem_edge_annotated.inp')
            latest_dictionary_path = os.path.join(rmg.outputDirectory, 'chemkin','species_edge_dictionary.txt')
            latest_transport_path = None
            saveChemkin(rmg.reactionModel, this_chemkin_path, latest_chemkin_verbose_path, latest_dictionary_path, latest_transport_path, rmg.saveEdgeSpecies)
            if os.path.exists(latest_chemkin_path):
                os.unlink(latest_chemkin_path)
            shutil.copy2(this_chemkin_path,latest_chemkin_path)

class ChemkinListener(object):
    """docstring for ChemkinListener"""
    def __init__(self):
        super(ChemkinListener, self).__init__()
    
    def update(self, rmg):
        saveChemkinFiles(rmg)

        
    