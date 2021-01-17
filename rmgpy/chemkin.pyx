###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
This module contains functions for writing of Chemkin input files.
"""

import logging
import math
import os.path
import re
import shutil
import textwrap
import warnings

import numpy as np

import rmgpy.constants as constants
import rmgpy.kinetics as _kinetics
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.exceptions import ChemkinError
from rmgpy.molecule.element import get_element
from rmgpy.molecule.util import get_element_count
from rmgpy.quantity import Quantity
from rmgpy.reaction import Reaction
from rmgpy.rmg.pdep import PDepNetwork, PDepReaction
from rmgpy.species import Species
from rmgpy.thermo import NASAPolynomial, NASA
from rmgpy.transport import TransportData
from rmgpy.util import make_output_subdirectory

_chemkin_reaction_count = None

################################################################################


def fortran_float(string):
    """
    Parse a Fortran-ish string into a float, like "1.00D 03" or "1.00d+03"
    """
    return float(string.replace("e", "E").replace("d", "D").replace("D", "E").replace("E ", "E+"))


def read_thermo_entry(entry, Tmin=0, Tint=0, Tmax=0):
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
    for i in [24, 29, 34, 39, 73]:
        element, count = lines[0][i:i + 2].strip(), lines[0][i + 2:i + 5].strip()
        if element:
            try:
                count = int(count)
            except ValueError:
                # Chemkin allows float values for the number of atoms, so try this next.
                try:
                    count = int(float(count))
                except ValueError:
                    logging.info("Trouble reading line '{0}' element segment "
                                 "'{1}'".format(lines[0].strip(),lines[0][i:i+5]))
                    if count == '' and re.match('\.?0*', element):
                        if i == 74:
                            logging.info("Assuming it's spillover from Tint, and ignoring.")
                        else:
                            logging.warning("Assuming it's not meant to be there, although it would be "
                                            "good to fix the chemkin file.")
                        count = 0
                    else:
                        raise
            if count != 0:  # Some people put garbage elements in, with zero count. Ignore these. Allow negative counts though (eg. negative one electron)
                formula[element] = count

    # Parsing for extended elemental composition syntax, adapted from Cantera ck2cti.py
    if lines[0].rstrip().endswith('&'):
        complines = []
        for i in range(len(lines)-1):
            if lines[i].rstrip().endswith('&'):
                complines.append(lines[i+1])
            else:
                break
        lines = [lines[0]] + lines[i+1:]
        elements = ' '.join(line.rstrip('&\n') for line in complines).split()
        formula = {}
        for i in range(0, len(elements), 2):
            formula[elements[i].capitalize()] = int(elements[i+1])

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
        a0_high = fortran_float(lines[1][0:15].strip())
        a1_high = fortran_float(lines[1][15:30].strip())
        a2_high = fortran_float(lines[1][30:45].strip())
        a3_high = fortran_float(lines[1][45:60].strip())
        a4_high = fortran_float(lines[1][60:75].strip())

        a5_high = fortran_float(lines[2][0:15].strip())
        a6_high = fortran_float(lines[2][15:30].strip())
        a0_low = fortran_float(lines[2][30:45].strip())
        a1_low = fortran_float(lines[2][45:60].strip())
        a2_low = fortran_float(lines[2][60:75].strip())

        a3_low = fortran_float(lines[3][0:15].strip())
        a4_low = fortran_float(lines[3][15:30].strip())
        a5_low = fortran_float(lines[3][30:45].strip())
        a6_low = fortran_float(lines[3][45:60].strip())
    except (IndexError, ValueError) as e:
        logging.warning('Error while reading thermo entry for species {0}'.format(species))
        logging.warning(str(e))
        return species, None, None

    # Construct and return the thermodynamics model
    thermo = NASA(
        polynomials=[
            NASAPolynomial(Tmin=(Tmin, "K"), Tmax=(Tint, "K"),
                           coeffs=[a0_low, a1_low, a2_low, a3_low, a4_low, a5_low, a6_low]),
            NASAPolynomial(Tmin=(Tint, "K"), Tmax=(Tmax, "K"),
                           coeffs=[a0_high, a1_high, a2_high, a3_high, a4_high, a5_high, a6_high])
        ],
        Tmin=(Tmin, "K"),
        Tmax=(Tmax, "K"),
    )
    if comment:
        thermo.comment = comment.strip()

    return species, thermo, formula

################################################################################


def read_kinetics_entry(entry, species_dict, Aunits, Eunits):
    """
    Read a kinetics `entry` for a single reaction as loaded from a Chemkin
    file. The associated mapping of labels to species `species_dict` should also
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
    reaction, third_body, kinetics, k_units, k_low_units = _read_kinetics_reaction(
        line=lines[0], species_dict=species_dict, Aunits=Aunits, Eunits=Eunits)

    if len(lines) == 1 and not third_body:
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
            kinetics = _read_kinetics_line(
                line=line, reaction=reaction, species_dict=species_dict, Eunits=Eunits,
                kunits=k_units, klow_units=k_low_units,
                kinetics=kinetics)

        # Decide which kinetics to keep and store them on the reaction object
        # Only one of these should be true at a time!
        if 'chebyshev' in kinetics:
            chebyshev = kinetics['chebyshev']
            if chebyshev.Tmin is None or chebyshev.Tmax is None:
                raise ChemkinError('Missing TCHEB line for reaction {0}'.format(reaction))
            if chebyshev.Pmin is None or chebyshev.Pmax is None:
                raise ChemkinError('Missing PCHEB line for reaction {0}'.format(reaction))
            if len(kinetics['chebyshev coefficients']) != (chebyshev.degreeT * chebyshev.degreeP):
                raise ChemkinError('Wrong number of Chebyshev coefficients '
                    'for reaction {0}'.format(reaction))
            index = 0
            for t in range(chebyshev.degreeT):
                for p in range(chebyshev.degreeP):
                    chebyshev.coeffs.value_si[t, p] = kinetics[
                        'chebyshev coefficients'][index]
                    index += 1
            # Don't forget to convert the Chebyshev coefficients to SI units!
            # This assumes that s^-1, cm^3/mol*s, etc. are compulsory
            chebyshev.coeffs.value_si[0, 0] -= (len(reaction.reactants) - 1) * math.log10(Afactor)
            reaction.kinetics = chebyshev
        elif 'pressure-dependent arrhenius' in kinetics:
            pdep_arrhenius = kinetics['pressure-dependent arrhenius']
            # Check for duplicates and combine them to MultiArrhenius objects
            duplicates_to_remove = []
            duplicates_to_add = []
            for index1 in range(len(pdep_arrhenius)):
                reaction1 = pdep_arrhenius[index1]
                p1, kinetics1 = reaction1
                if reaction1 in duplicates_to_remove:
                    continue
                for index2 in range(index1 + 1, len(pdep_arrhenius)):
                    reaction2 = pdep_arrhenius[index2]
                    p2, kinetics2 = reaction2
                    if p1 == p2:
                        if reaction1 not in duplicates_to_remove:
                            new_kinetics = _kinetics.MultiArrhenius()
                            duplicates_to_add.append([p1, new_kinetics])
                            new_kinetics.arrhenius = [kinetics1]
                            duplicates_to_remove.append(reaction1)
                        new_kinetics.arrhenius.append(kinetics2)
                        duplicates_to_remove.append(reaction2)
            for item in duplicates_to_remove:
                pdep_arrhenius.remove(item)
            pdep_arrhenius.extend(duplicates_to_add)

            pdep_arrhenius = sorted(pdep_arrhenius, key=lambda reaction: reaction[0])  # sort by ascending pressures

            reaction.kinetics = _kinetics.PDepArrhenius(
                pressures=([P for P, arrh in pdep_arrhenius], "atm"),
                arrhenius=[arrh for P, arrh in pdep_arrhenius],
            )
        elif 'troe' in kinetics:
            troe = kinetics['troe']
            troe.arrheniusHigh = kinetics['arrhenius high']
            troe.arrheniusLow = kinetics['arrhenius low']
            troe.efficiencies = kinetics['efficiencies']
            reaction.kinetics = troe
        elif third_body:
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
        elif 'sticking coefficient' in kinetics:
            reaction.kinetics = kinetics['sticking coefficient']
        elif 'surface arrhenius' in kinetics:
            reaction.kinetics = kinetics['surface arrhenius']
        else:
            raise ChemkinError(
                'Unable to understand all additional information lines for reaction {0}.'.format(entry))

        # These things may *also* be true
        if 'sri' in kinetics:
            reaction.kinetics.comment += "Warning: SRI parameters from chemkin file ignored on import. "

        if 'explicit reverse' in kinetics:
            reaction.kinetics.comment += \
                "Chemkin file stated explicit reverse rate: {0}".format(kinetics['explicit reverse'])

    return reaction


def _read_kinetics_reaction(line, species_dict, Aunits, Eunits):
    """
    Parse the first line of of a Chemkin reaction entry.
    """
    tokens = line.split()

    rmg = True
    try:
        float(tokens[-6])
    except (ValueError, IndexError):
        rmg = False
    A_uncertainty_type = '+|-'
    if rmg:
        A = float(tokens[-6])
        n = float(tokens[-5])
        Ea = float(tokens[-4])
        try:
            dA = float(tokens[-3])
        except ValueError:
            A_uncertainty_type = '*|/'
            dA = float(tokens[-3][1:])
        dn = float(tokens[-2])
        dEa = float(tokens[-1])
        reaction = ''.join(tokens[:-6])
    else:
        A = fortran_float(tokens[-3])
        n = fortran_float(tokens[-2])
        Ea = fortran_float(tokens[-1])
        dA = 0.0
        dn = 0.0
        dEa = 0.0
        reaction = ''.join(tokens[:-3])
    third_body = False

    # Split the reaction equation into reactants and products
    reversible = True
    reactants, products = reaction.split('=')
    if '<=>' in reaction:
        reactants = reactants[:-1]
        products = products[1:]
    elif '=>' in reaction:
        products = products[1:]
        reversible = False
    specific_collider = None
    # search for a third body collider, e.g., '(+M)', '(+m)', or a specific species like '(+N2)',
    #     matching `(+anything_other_than_ending_parenthesis)`:
    collider = re.search(r'\(\+[^)]+\)', reactants)
    if collider is not None:
        collider = collider.group(0)  # save string value rather than the object
        if collider != re.search(r'\(\+[^)]+\)', products).group(0):
            raise ChemkinError(
                'Third body colliders in reactants and products of reaction {0} are not identical!'.format(reaction))
        extra_parenthesis = collider.count('(') - 1
        for i in range(extra_parenthesis):
            # allow for species like N2(5) or CH2(T)(15) to be read as specific colliders,
            #     although currently not implemented in Chemkin. See RMG-Py #1070
            collider += ')'
        reactants = reactants.replace(collider, '')
        products = products.replace(collider, '')
        if collider.upper().strip() != "(+M)":  # the collider is a specific species, not (+M) or (+m)
            if collider.strip()[2:-1] not in species_dict:  # stripping spaces, '(+' and ')'
                raise ChemkinError(
                    'Unexpected third body collider "{0}" in reaction {1}.'.format(collider.strip()[2:-1], reaction))
            specific_collider = species_dict[collider.strip()[2:-1]]

    # Create a new Reaction object for this reaction
    reaction = Reaction(reactants=[], products=[], specific_collider=specific_collider, reversible=reversible)

    # Convert the reactants and products to Species objects using the species_dict
    for reactant in reactants.split('+'):
        reactant = reactant.strip()
        stoichiometry = 1
        if reactant not in species_dict and reactant[0].isdigit():
            # This allows for reactions to be of the form 2A=B+C instead of A+A=B+C
            # The implementation below assumes an integer between 0 and 9, inclusive
            stoichiometry = int(reactant[0])
            reactant = reactant[1:]
        if reactant.upper() == 'M':
            # this identifies reactions like 'H+H+M=H2+M' as opposed to 'H+H(+M)=H2(+M)' as identified above
            third_body = True
        elif reactant not in species_dict:
            raise ChemkinError('Unexpected reactant "{0}" in reaction {1}.'.format(reactant, reaction))
        else:
            reactant_species = species_dict[reactant]
            if not reactant_species.reactive:
                reactant_species.reactive = True
            for i in range(stoichiometry):
                reaction.reactants.append(reactant_species)
    for product in products.split('+'):
        product = product.strip()
        stoichiometry = 1
        if product not in species_dict and product[0].isdigit():
            # This allows for reactions to be of the form A+B=2C instead of A+B=C+C
            # The implementation below assumes an integer between 0 and 9, inclusive
            stoichiometry = int(product[0])
            product = product[1:]
        if product.upper() == 'M':
            pass
        elif product not in species_dict:
            if re.match('[0-9.]+', product):
                logging.warning("Looks like reaction {0!r} has fractional stoichiometry, which RMG cannot handle. "
                                "Ignoring".format(line))
                raise ChemkinError('Skip reaction!')
            raise ChemkinError(
                'Unexpected product "{0}" in reaction {1} from line {2}.'.format(product, reaction, line))
        else:
            product_species = species_dict[product]
            if not product_species.reactive:
                product_species.reactive = True
            for i in range(stoichiometry):
                reaction.products.append(product_species)

    # Determine the appropriate units for k(T) and k(T,P) based on the number of reactants
    # This assumes elementary kinetics for all reactions
    try:
        n_react = len(reaction.reactants) + (1 if third_body else 0)
        k_units = Aunits[n_react]
        k_low_units = Aunits[n_react + 1]
    except IndexError:
        raise ChemkinError('Invalid number of reactant species for reaction {0}.'.format(reaction))

    key = 'arrhenius low' if third_body else 'arrhenius high'

    kinetics = {
        key: _kinetics.Arrhenius(
            A=(A, k_units, A_uncertainty_type, dA),
            n=(n, '', '+|-', dn),
            Ea=(Ea, Eunits, '+|-', dEa),
            T0=(1, "K"),
        ),
    }
    return reaction, third_body, kinetics, k_units, k_low_units


def _read_kinetics_line(line, reaction, species_dict, Eunits, kunits, klow_units, kinetics):
    """
    Parse the subsequent lines of of a Chemkin reaction entry.
    """
    case_preserved_tokens = line.split('/')
    line = line.upper()
    tokens = line.split('/')

    if 'DUP' in line:
        # Duplicate reaction
        reaction.duplicate = True

    elif 'COV' in line:
        try:
            k = kinetics['sticking coefficient']
        except KeyError:
            k = kinetics['arrhenius high']
            k = _kinetics.SurfaceArrhenius(
                A=(k.A.value, kunits),
                n=k.n,
                Ea=k.Ea,
                T0=k.T0,
            )
            kinetics['surface arrhenius'] = k
            del kinetics['arrhenius high']

        tokens = case_preserved_tokens[1].split()
        cov_dep_species = species_dict[tokens[0].strip()]
        k.coverage_dependence[cov_dep_species] = {'a':float(tokens[1]), 'm':float(tokens[2]), 'E':(float(tokens[3]), Eunits)}

    elif 'LOW' in line:
        # Low-pressure-limit Arrhenius parameters
        tokens = tokens[1].split()
        kinetics['arrhenius low'] = _kinetics.Arrhenius(
            A=(float(tokens[0].strip()), klow_units),
            n=float(tokens[1].strip()),
            Ea=(float(tokens[2].strip()), Eunits),
            T0=(1, "K"),
        )

    elif 'HIGH' in line:
        # What we thought was high, was in fact low-pressure
        kinetics['arrhenius low'] = kinetics['arrhenius high']
        kinetics['arrhenius low'].A = (
            kinetics['arrhenius low'].A.value, klow_units)
        # High-pressure-limit Arrhenius parameters
        tokens = tokens[1].split()
        kinetics['arrhenius high'] = _kinetics.Arrhenius(
            A=(float(tokens[0].strip()), kunits),
            n=float(tokens[1].strip()),
            Ea=(float(tokens[2].strip()), Eunits),
            T0=(1, "K"),
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
            alpha=alpha,
            T3=(T3, "K"),
            T1=(T1, "K"),
            T2=(T2, "K") if T2 is not None else None,
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
            tokens2 = tokens[index + 1].split()
            chebyshev.Tmin = Quantity(float(tokens2[0].strip()), "K")
            chebyshev.Tmax = Quantity(float(tokens2[1].strip()), "K")
        if 'PCHEB' in line:
            index = tokens.index('PCHEB')
            tokens2 = tokens[index + 1].split()
            chebyshev.Pmin = Quantity(float(tokens2[0].strip()), "atm")
            chebyshev.Pmax = Quantity(float(tokens2[1].strip()), "atm")
        if 'TCHEB' in line or 'PCHEB' in line:
            pass
        elif chebyshev.degreeT == 0 or chebyshev.degreeP == 0:
            tokens2 = tokens[1].split()
            chebyshev.degreeT = int(float(tokens2[0].strip()))
            chebyshev.degreeP = int(float(tokens2[1].strip()))
            chebyshev.coeffs = np.zeros((chebyshev.degreeT, chebyshev.degreeP), np.float64)
            # There may be some coefficients on this first line
            kinetics['chebyshev coefficients'].extend(
                [float(t.strip()) for t in tokens2[2:]])
        else:
            tokens2 = tokens[1].split()
            kinetics['chebyshev coefficients'].extend(
                [float(t.strip()) for t in tokens2])

    elif 'PLOG' in line:
        pdep_arrhenius = kinetics.get('pressure-dependent arrhenius', [])
        kinetics['pressure-dependent arrhenius'] = pdep_arrhenius
        tokens = tokens[1].split()
        pdep_arrhenius.append([float(tokens[0].strip()),
                               _kinetics.Arrhenius(
                                   A=(float(tokens[1].strip()), kunits),
                                   n=float(tokens[2].strip()),
                                   Ea=(float(tokens[3].strip()), Eunits),
                                   T0=(1, "K"),
                               )])

    elif tokens[0].startswith('REV'):
        reverse_A = float(tokens[1].split()[0])
        kinetics['explicit reverse'] = line.strip()
        if reverse_A == 0:
            logging.info("Reverse rate is 0 so making irreversible for reaction {0}".format(reaction))
            reaction.reversible = False
        else:
            logging.info("Ignoring explicit reverse rate for reaction {0}".format(reaction))

    elif line.strip() == 'STICK':
        # Convert what we thought was Arrhenius into StickingCoefficient
        k = kinetics['arrhenius high']
        kinetics['sticking coefficient'] = _kinetics.StickingCoefficient(
            A=k.A.value,
            n=k.n,
            Ea=k.Ea,
            T0=k.T0,
        )
        del kinetics['arrhenius high']

    else:
        # Assume a list of collider efficiencies
        try:
            for collider, efficiency in zip(case_preserved_tokens[0::2], case_preserved_tokens[1::2]):
                try:
                    efficiency = float(efficiency.strip())
                except ValueError:
                    raise ChemkinError("{0!r} doesn't look like a collision efficiency for species {1} in "
                                       "line {2!r}".format(efficiency, collider.strip(), line))
                if collider.strip() in species_dict:
                    kinetics['efficiencies'][species_dict[collider.strip()].molecule[0]] = efficiency
                else:  # try it with capital letters? Not sure whose malformed chemkin files this is needed for.
                    kinetics['efficiencies'][species_dict[collider.strip().upper()].molecule[0]] = efficiency
        except IndexError:
            error_msg = 'Could not read collider efficiencies for reaction: {0}.\n'.format(reaction)
            error_msg += 'The following line was parsed incorrectly:\n{0}'.format(line)
            error_msg += "\n(Case-preserved tokens: {0!r} )".format(case_preserved_tokens)
            raise ChemkinError(error_msg)
    return kinetics


def _remove_line_breaks(comments):
    """
    This method removes any extra line breaks in reaction comments, so they
    can be parsed by read_reaction_comments.
    """
    comments = comments.replace('\n', ' ')
    new_statement_indicators = ['Reaction index', 'Template reaction', 'Library reaction',
                                'PDep reaction', 'Flux pairs', 'BM rule fitted to',
                                'Uncertainty in Total Std:',
                                'Estimated using', 'Exact match found', 'Average of ',
                                'Euclidian distance', 'Matched node ', 'Matched reaction ',
                                'Multiplied by reaction path degeneracy ',
                                'Kinetics were estimated in this direction',
                                'dGrxn(298 K) = ', 'Both directions are estimates',
                                'Other direction matched ', 'Both directions matched ',
                                'This direction matched an entry in ', 'From training reaction',
                                'This reaction matched rate rule', 'family: ', 'Warning:',
                                'Chemkin file stated explicit reverse rate:', 'Ea raised from',
                                'Fitted to', 'Reaction library', 'Estimated from node', 'Matched node',
                                ]
    for indicator in new_statement_indicators:
        comments = comments.replace(' ' + indicator, '\n' + indicator, 1)
    return comments


def read_reaction_comments(reaction, comments, read=True):
    """
    Parse the `comments` associated with a given `reaction`. If the comments
    come from RMG (Py or Java), parse them and extract the useful information.
    Return the reaction object based on the information parsed from these
    comments. If `read` if False, the reaction is returned as an "Unclassified"
    LibraryReaction.
    """

    if not read:
        # The chemkin file was not generated by either RMG-Py or RMG-Java, thus, there should be no parsing
        # of the comments.  Instead, return as an unclassified LibraryReaction.
        reaction = LibraryReaction(
            index=reaction.index,
            reactants=reaction.reactants,
            products=reaction.products,
            specific_collider=reaction.specific_collider,
            kinetics=reaction.kinetics,
            reversible=reaction.reversible,
            duplicate=reaction.duplicate,
            library='Unclassified',
        )

        return reaction

    # the comments could have line breaks that will mess up reading
    # we will now combine the lines and separate them based on statements
    comments = _remove_line_breaks(comments)
    lines = comments.strip().splitlines()

    for line in lines:

        tokens = line.split()
        if 'Reaction index:' in line:
            # Don't store the reaction indices
            pass

        elif 'Template reaction:' in line:
            label = str(tokens[-1])
            reaction = TemplateReaction(
                index=reaction.index,
                reactants=reaction.reactants,
                products=reaction.products,
                specific_collider=reaction.specific_collider,
                kinetics=reaction.kinetics,
                reversible=reaction.reversible,
                duplicate=reaction.duplicate,
                family=label,
            )

        elif 'Library reaction:' in line or 'Seed mechanism:' in line:
            label = str(tokens[2])
            reaction = LibraryReaction(
                index=reaction.index,
                reactants=reaction.reactants,
                products=reaction.products,
                specific_collider=reaction.specific_collider,
                kinetics=reaction.kinetics,
                reversible=reaction.reversible,
                duplicate=reaction.duplicate,
                library=label,
            )

        elif 'PDep reaction:' in line:
            network_index = int(tokens[-1][1:])
            reaction = PDepReaction(
                index=reaction.index,
                reactants=reaction.reactants,
                products=reaction.products,
                specific_collider=reaction.specific_collider,
                kinetics=reaction.kinetics,
                reversible=reaction.reversible,
                duplicate=reaction.duplicate,
                network=PDepNetwork(index=network_index),
            )

        elif 'Flux pairs:' in line:
            reaction.pairs = []
            for reac_str, prod_str in zip(tokens[2::2], tokens[3::2]):
                if reac_str[-1] == ',':
                    reac_str = reac_str[:-1]
                for reactant in reaction.reactants:
                    if reactant.label == reac_str:
                        break
                else:
                    raise ChemkinError('Unexpected species identifier {0} encountered in flux pairs '
                                       'for reaction {1}.'.format(reac_str, reaction))
                if prod_str[-1] == ';':
                    prod_str = prod_str[:-1]
                for product in reaction.products:
                    if product.label == prod_str:
                        break
                else:
                    raise ChemkinError('Unexpected species identifier {0} encountered in flux pairs '
                                       'for reaction {1}.'.format(prod_str, reaction))
                reaction.pairs.append((reactant, product))
            assert len(reaction.pairs) == max(len(reaction.reactants), len(reaction.products))

        elif isinstance(reaction, TemplateReaction) and 'rate rule ' in line:
            bracketed_rule = tokens[-1]
            templates = bracketed_rule[1:-1].split(';')
            reaction.template = templates
            # still add kinetic comment
            reaction.kinetics.comment += line.strip() + "\n"

        elif isinstance(reaction, TemplateReaction) and \
                'Multiplied by reaction path degeneracy ' in line:
            degen = float(tokens[-1])
            reaction.degeneracy = degen
            # undo the kinetic manipulation caused by setting degneracy
            if reaction.kinetics:
                reaction.kinetics.change_rate(1. / degen)
            # do not add comment because setting degeneracy does so already
            reaction.kinetics.comment += "\n"

        elif 'BM rule fitted to' in line:
            reaction.kinetics.comment += line.strip() + "\n"

        elif 'Uncertainty in Total Std:' in line:
            reaction.kinetics.comment += line.strip() + "\n"

        elif line.strip() != '':
            # Any lines which are commented out but don't have any specific flag are simply kinetics comments
            reaction.kinetics.comment += line.strip() + "\n"


        # Comment parsing from old RMG-Java chemkin files
        elif 'PDepNetwork' in line:
            network_index = int(tokens[3][1:])
            reaction = PDepReaction(
                index=reaction.index,
                reactants=reaction.reactants,
                products=reaction.products,
                specific_collider=reaction.specific_collider,
                kinetics=reaction.kinetics,
                reversible=reaction.reversible,
                duplicate=reaction.duplicate,
                network=PDepNetwork(index=network_index)
            )
            reaction.kinetics.comment = line

        elif 'ReactionLibrary:' in line or 'Seed Mechanism:' in line:
            label = str(tokens[-1])
            reaction = LibraryReaction(
                index=reaction.index,
                reactants=reaction.reactants,
                products=reaction.products,
                specific_collider=reaction.specific_collider,
                kinetics=reaction.kinetics,
                reversible=reaction.reversible,
                duplicate=reaction.duplicate,
                library=label,
            )
            reaction.kinetics.comment = line

        elif 'exact:' in line or 'estimate:' in line:
            index1 = line.find('[')
            index2 = line.find(']')
            template = [s.strip() for s in line[index1:index2].split(',')]
            label = str(tokens[0])
            reaction = TemplateReaction(
                index=reaction.index,
                reactants=reaction.reactants,
                products=reaction.products,
                specific_collider=reaction.specific_collider,
                kinetics=reaction.kinetics,
                reversible=reaction.reversible,
                duplicate=reaction.duplicate,
                family=label,
                template=[Entry(label=g) for g in template],
            )
            reaction.kinetics.comment = line

    if (not isinstance(reaction, LibraryReaction)
            and not isinstance(reaction, TemplateReaction)
            and not isinstance(reaction, PDepReaction)):
        reaction = LibraryReaction(
            index=reaction.index,
            reactants=reaction.reactants,
            products=reaction.products,
            specific_collider=reaction.specific_collider,
            kinetics=reaction.kinetics,
            reversible=reaction.reversible,
            duplicate=reaction.duplicate,
            library='Unclassified',
        )

    # Clean up line endings on comments so that there aren't any blank commented lines
    reaction.kinetics.comment = reaction.kinetics.comment.strip()
    return reaction

################################################################################


def load_species_dictionary(path, generate_resonance_structures=True):
    """
    Load an RMG dictionary - containing species identifiers and the associated
    adjacency lists - from the file located at `path` on disk. Returns a dict
    mapping the species identifiers to the loaded species.
    If `generate_resonance_structures` is True (default if omitted)
    then resonance isomers for each species are generated.
    """
    from rmgpy.molecule.fragment import Fragment
    import re
    species_dict = {}

    inerts = [Species().from_smiles(inert) for inert in ('[He]', '[Ne]', 'N#N', '[Ar]')]
    with open(path, 'r') as f:
        adjlist = ''
        for line in f:
            if line.strip() == '' and adjlist.strip() != '':
                # Finish this adjacency list
                if len(re.findall(r'([LR]\d?)', adjlist)) != 0:
                    frag = Fragment().from_adjacency_list(adjlist)
                    species = Species(molecule = [frag])
                    for label in adjlist.splitlines():
                        if label.strip():
                            break
                    else:
                        label = ''
                    if len(label.split()) > 0 and not label.split()[0].isdigit():
                        species.label = label.strip()
                else:
                    species = Species().from_adjacency_list(adjlist)
                if generate_resonance_structures:
                    species.generate_resonance_structures()
                label = species.label
                for inert in inerts:
                    if inert.is_isomorphic(species):
                        species.reactive = False
                        break
                species_dict[label] = species
                adjlist = ''
            else:
                if "InChI" in line:
                    line = line.split()[0] + '\n'
                if '//' in line:
                    index = line.index('//')
                    line = line[0:index]
                adjlist += line
        else:  #reach end of file
            if adjlist.strip() != '':
                if len(re.findall(r'([LR]\d?)', adjlist)) != 0:
                    frag = Fragment().from_adjacency_list(adjlist)
                    species = Species(molecule = [frag])
                    for label in adjlist.splitlines():
                        if label.strip():
                            break
                    else:
                        label = ''
                    if len(label.split()) > 0 and not label.split()[0].isdigit():
                        species.label = label.strip()
                else:
                    species = Species().from_adjacency_list(adjlist)
                if generate_resonance_structures:
                    species.generate_resonance_structures()
                label = species.label
                for inert in inerts:
                    if inert.is_isomorphic(species):
                        species.reactive = False
                        break
                species_dict[label] = species

    return species_dict


def remove_comment_from_line(line):
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
    comment = line[index + 1:-1]
    if index < len(line):
        line = line[0:index] + '\n'

    return line, comment


def load_transport_file(path, species_dict):
    """
    Load a Chemkin transport properties file located at `path` and store the
    properties on the species in `species_dict`.
    """
    with open(path, 'r') as f:
        for line0 in f:
            line, comment = remove_comment_from_line(line0)
            line = line.strip()
            if line != '':
                # This line contains an entry, so parse it
                label = line[0:16].strip()
                data = line[16:].split()
                species = species_dict[label]
                species.transport_data = TransportData(
                    shapeIndex=int(data[0]),
                    sigma=(float(data[2]), 'angstrom'),
                    epsilon=(float(data[1]), 'K'),
                    dipoleMoment=(float(data[3]), 'De'),
                    polarizability=(float(data[4]), 'angstrom^3'),
                    rotrelaxcollnum=float(data[5]),
                    comment=comment.strip(),
                )


def load_chemkin_file(path, dictionary_path=None, transport_path=None, read_comments=True, thermo_path=None,
                      use_chemkin_names=False, check_duplicates=True, generate_resonance_structures=True):
    """
    Load a Chemkin input file located at `path` on disk to `path`, returning lists of the species
    and reactions in the Chemkin file. The 'thermo_path' point to a separate thermo file, or, if 'None' is
    specified, the function will look for the thermo database within the chemkin mechanism file.
    If `generate_resonance_structures` is True (default if omitted) then resonance isomers for
    each species are generated.
    """
    species_list = []
    species_dict = {}
    species_aliases = {}
    reaction_list = []

    # If the dictionary path is given, then read it and generate Molecule objects
    # You need to append an additional adjacency list for nonreactive species, such
    # as N2, or else the species objects will not store any structures for the final
    # HTML output.
    if dictionary_path:
        species_dict = load_species_dictionary(dictionary_path, generate_resonance_structures=generate_resonance_structures)

    with open(path, 'r') as f:
        previous_line = f.tell()
        line0 = f.readline()
        while line0 != '':
            line = remove_comment_from_line(line0)[0]
            line = line.strip()

            if 'SPECIES' in line.upper():
                # Unread the line (we'll re-read it in readReactionBlock())
                f.seek(previous_line)
                read_species_block(f, species_dict, species_aliases, species_list)
            
            elif 'SITE' in line.upper():
                # Unread the line (we'll re-read it in readReactionBlock())
                f.seek(previous_line)
                read_species_block(f, species_dict, species_aliases, species_list)

            elif 'THERM' in line.upper() and thermo_path is None:
                # Skip this if a thermo file is specified
                # Unread the line (we'll re-read it in read_thermo_block())
                f.seek(previous_line)
                read_thermo_block(f, species_dict)

            elif 'REACTIONS' in line.upper():
                # Reactions section
                # Unread the line (we'll re-read it in readReactionBlock())
                f.seek(previous_line)
                reaction_list = read_reactions_block(f, species_dict, read_comments=read_comments)

            previous_line = f.tell()
            line0 = f.readline()

    # Read in the thermo data from the thermo file        
    if thermo_path:
        with open(thermo_path, 'r') as f:
            line0 = f.readline()
            while line0 != '':
                line = remove_comment_from_line(line0)[0]
                line = line.strip()
                if 'THERM' in line.upper():
                    f.seek(-len(line0), 1)
                    read_thermo_block(f, species_dict)
                    break
                line0 = f.readline()
    # Index the reactions now to have identical numbering as in Chemkin
    index = 0
    for reaction in reaction_list:
        index += 1
        reaction.index = index

    # Process duplicate reactions
    if check_duplicates:
        _process_duplicate_reactions(reaction_list)

    # If the transport path is given, then read it to obtain the transport
    # properties
    if transport_path:
        load_transport_file(transport_path, species_dict)

    if not use_chemkin_names:
        # Apply species aliases if known
        for spec in species_list:
            try:
                spec.label = species_aliases[spec.label]
            except KeyError:
                pass

    # Attempt to extract index from species label
    indexPattern = re.compile(r'\(\d+\)$')
    for spec in species_list:
        if indexPattern.search(spec.label):
            label, sep, index = spec.label[:-1].rpartition('(')
            spec.label = label
            spec.index = int(index)

    reaction_list.sort(key=lambda reaction: reaction.index)
    return species_list, reaction_list


cpdef _process_duplicate_reactions(list reaction_list):
    """
    Check for marked (and unmarked!) duplicate reactions
    Combine marked duplicate reactions into a single reaction using MultiKinetics
    Raise exception for unmarked duplicate reactions
    """
    cdef list duplicate_reactions_to_remove = []
    cdef list duplicate_reactions_to_add = []
    cdef int index1, index2
    cdef Reaction reaction, reaction1, reaction2
    cdef KineticsModel kinetics

    for index1 in range(len(reaction_list)):
        reaction1 = reaction_list[index1]
        if reaction1 in duplicate_reactions_to_remove:
            continue

        for index2 in range(index1 + 1, len(reaction_list)):
            reaction2 = reaction_list[index2]
            if (reaction1.reactants == reaction2.reactants
                    and reaction1.products == reaction2.products
                    and reaction1.specific_collider == reaction2.specific_collider):
                if reaction1.duplicate and reaction2.duplicate:

                    if isinstance(reaction1, LibraryReaction) and isinstance(reaction2, LibraryReaction):
                        if reaction1.library != reaction2.library:
                            raise ChemkinError("Identical reactions {0} and {1} taken from different libraries: {2}, "
                                               "{3}".format(reaction1, reaction2, reaction1.library, reaction2.library))
                        if reaction1 not in duplicate_reactions_to_remove:
                            # already created duplicate reaction, move on to appending any additional duplicate kinetics
                            if isinstance(reaction1.kinetics,
                                          _kinetics.PDepArrhenius):
                                kinetics = _kinetics.MultiPDepArrhenius()
                            elif isinstance(reaction1.kinetics,
                                            _kinetics.Arrhenius):
                                kinetics = _kinetics.MultiArrhenius()
                            else:
                                logging.warning(
                                    'Unexpected kinetics type {0} for duplicate reaction {1}. '
                                    'Not combining reactions.'.format(reaction1.kinetics.__class__, reaction1)
                                )
                                continue
                            reaction = LibraryReaction(
                                index=reaction1.index,
                                reactants=reaction1.reactants,
                                products=reaction1.products,
                                specific_collider=reaction1.specific_collider,
                                kinetics=kinetics,
                                library=reaction1.library,
                                duplicate=False,
                            )
                            duplicate_reactions_to_add.append(reaction)
                            kinetics.arrhenius = [reaction1.kinetics]
                            duplicate_reactions_to_remove.append(reaction1)

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

                    duplicate_reactions_to_remove.append(reaction2)
                elif reaction1.kinetics.is_pressure_dependent() == reaction2.kinetics.is_pressure_dependent():
                    # If both reactions are pressure-independent or both are pressure-dependent, then they need
                    # duplicate tags. Chemkin treates pdep and non-pdep reactions as different, so those are okay
                    raise ChemkinError('Encountered unmarked duplicate reaction {0}.'.format(reaction1))

    for reaction in duplicate_reactions_to_remove:
        reaction_list.remove(reaction)
    reaction_list.extend(duplicate_reactions_to_add)


def read_species_block(f, species_dict, species_aliases, species_list):
    """
    Read a Species block from a chemkin file.
    
    f is a file-like object that is just before the 'SPECIES' statement. When finished, it will have just passed the 'END' statement.
    species_dict is a dictionary of species that will be updated.
    species_aliases is a dictionary of species aliases that will be updated.
    species_list is a list of species that will be extended.
    """
    line = f.readline()
    line = remove_comment_from_line(line)[0]
    line = line.strip()
    tokens = line.split()
    tokens_upper = line.upper().split()
    first_token = tokens.pop(0)
    first_token = tokens_upper.pop(0)  # pop from both lists
    assert first_token in ['SPECIES', 'SPEC', 'SITE']  # should be first token in first line
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
                species_aliases[label] = alias
        line = remove_comment_from_line(line)[0]
        line = line.strip()
        tokens.extend(line.split())
        tokens_upper.extend(line.upper().split())
    # Now process list of tokens
    processed_tokens = []
    for token in tokens:
        if token in processed_tokens:
            continue  # ignore species declared twice
        token_upper = token.upper()
        if token_upper in ['SPECIES', 'SPEC', 'SITE']:
            continue  # there may be more than one SPECIES statement
        if token_upper == 'END':
            break
        processed_tokens.append(token)
        if token in species_dict:
            logging.debug("Re-using species {0} already in species_dict".format(token))
            species = species_dict[token]
        else:
            species = Species(label=token)
            species_dict[token] = species
        species_list.append(species)


def read_thermo_block(f, species_dict):
    """
    Read a thermochemistry block from a chemkin file.
    
    f is a file-like object that is just before the 'THERM' statement.
    When finished, it will have just passed the 'END' statement.
    species_dict is a dictionary of species that will be updated with the given thermodynamics.
    
    Returns a dictionary of molecular formulae for each species, in the form
    `{'methane': {'C':1, 'H':4}}
    
    If duplicate entries are found, the FIRST is used, and a warning is printed.
    """
    # List of thermodynamics (hopefully one per species!)
    formula_dict = {}
    line = f.readline()
    assert line.upper().strip().startswith('THER'), "'{0}' doesn't begin with THERM statement.".format(line)
    line = f.readline()

    # In case there are commented lines immediately after THER
    meaningfulline, comment = remove_comment_from_line(line)
    while not meaningfulline.strip():
        line = f.readline()
        meaningfulline, comment = remove_comment_from_line(line)
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
    thermo_block = ''
    comments = ''
    while line != '' and not line.upper().strip().startswith('END'):
        line, comment = remove_comment_from_line(line)
        if comment: comments += comment.strip().replace('\t', ', ') + '\n'

        if len(line) < 80:
            if line.strip():
                logging.info("Ignoring short but non-empty line: {0!r}".format(line))
            line = f.readline()
            continue

        if line[79] not in ['1', '2', '3', '4']:
            logging.warning("Ignoring line without 1,2,3 or 4 in 80th column: {0!r}".format(line))
            line = f.readline()
            continue

        thermo_block += line

        # check for extended elemental composition line
        if line.rstrip().endswith('1&'):
            # this thermo entry has extended elemental composition line
            # read the elements line, append it to thermo_block, and continue 
            line = f.readline()
            thermo_block += line
            line = f.readline()
            continue

        if line[79] == '4':
            try:
                label, thermo, formula = read_thermo_entry(thermo_block, Tmin=Tmin, Tint=Tint, Tmax=Tmax)
            except:
                logging.error("Error reading thermo block:\n" + thermo_block)
                raise
            if label not in species_dict:
                logging.info("Ignoring thermo data for {0} because it's not in the requested "
                             "list of species.".format(label))
                thermo_block = ''
                line = f.readline()
                continue
            elif species_dict[label].has_thermo():
                logging.warning('Skipping duplicate thermo for the species {0}'.format(label))
                thermo_block = ''
                line = f.readline()
                continue
            else:
                if thermo is None:
                    logging.error("Problematic thermo block:\n{0}".format(thermo_block))
                    raise ChemkinError('Error while reading thermo entry for required species {0}'.format(label))
            try:
                formula_dict[label] = formula
                species_dict[label].thermo = thermo
                species_dict[label].thermo.comment = getattr(species_dict[label].thermo, 'comment', '')
                if comments:
                    species_dict[label].thermo.comment += '\n{0}'.format(comments)
                # Make sure to strip whitespace
                species_dict[label].thermo.comment = species_dict[label].thermo.comment.strip()
                comments = ''
            except KeyError:
                if label.upper() in ['AR', 'N2', 'HE', 'NE']:
                    logging.warning('Skipping species"{0}" while reading thermodynamics entry.'.format(label))
                else:
                    logging.warning('Skipping unexpected species "{0}" while reading '
                                    'thermodynamics entry.'.format(label))
            thermo_block = ''
        if len(thermo_block.split('/n')) > 5:
            raise ChemkinError('Should only have a maximum of 5 lines in a thermo block:\n{0}'.format(thermo_block))
        line = f.readline()
    return formula_dict


def read_reactions_block(f, species_dict, read_comments=True):
    """
    Read a reactions block from a Chemkin file stream.
    
    This function can also read the ``reactions.txt`` and ``pdepreactions.txt``
    files from RMG-Java kinetics libraries, which have a similar syntax.
    """
    energy_units = 'cal/mol'
    molecule_units = 'moles'
    volume_units = 'cm3'
    time_units = 's'

    line = f.readline()
    found = False
    while line != '' and not found:

        line = remove_comment_from_line(line)[0]
        line = line.strip()
        tokens = line.split()

        if len(tokens) > 0 and tokens[0].upper() == 'REACTIONS':
            # Regular Chemkin file
            found = True
            for token in tokens[1:]:
                unit = token.lower()
                if unit in ['molecules', 'moles', 'mole', 'mol', 'molecule']:
                    molecule_units = unit
                elif unit in ['kcal/mole', 'kcal/mol', 'cal/mole', 'cal/mol', 'kj/mole', 'kj/mol', 'j/mole', 'j/mol',
                              'kelvins']:
                    energy_units = unit
                else:
                    raise ChemkinError('Unknown unit type "{0}"'.format(unit))

        elif len(tokens) > 0 and tokens[0].lower() == 'unit:':
            # RMG-Java kinetics library file
            warnings.warn("The RMG-Java kinetic library files are"
                          " no longer supported and may be"
                          " removed in version 2.3.", DeprecationWarning)
            found = True
            while 'reactions:' not in line.lower():
                line = f.readline()
                line = remove_comment_from_line(line)[0]
                line = line.strip()

                if 'A:' in line or 'E:' in line:
                    units = line.split()[1]
                    if 'A:' in line:
                        molecule_units, volume_units, time_units = units.lower().split(
                            '/')  # Assume this is a 3-tuple: moles or molecules, volume, time
                    elif 'E:' in line:
                        energy_units = units.lower()
        else:
            line = f.readline()

    if not found:
        raise ChemkinError('Invalid reaction block.')

    # Check that the units are valid
    assert molecule_units in ['molecules', 'moles', 'mole', 'mol', 'molecule']
    assert volume_units in ['cm3', 'm3']
    assert time_units in ['s']
    assert energy_units in ['kcal/mole', 'kcal/mol', 'cal/mole', 'cal/mol', 'kj/mole', 'kj/mol', 'j/mole', 'j/mol',
                           'kelvins']

    # Homogenize units
    if molecule_units == 'molecules':
        molecule_units = 'molecule'
    elif molecule_units == 'moles' or molecule_units == 'mole':
        molecule_units = 'mol'
    volume_units = {'cm3': 'cm', 'm3': 'm'}[volume_units]
    if energy_units == 'kcal/mole':
        energy_units = 'kcal/mol'
    elif energy_units == 'cal/mole':
        energy_units = 'cal/mol'
    elif energy_units == 'kj/mole':
        energy_units = 'kj/mol'
    elif energy_units == 'j/mole':
        energy_units = 'j/mol'
    elif energy_units == 'kelvins':
        energy_units = 'K'
    energy_units = energy_units.replace('j/mol', 'J/mol')

    # Set up kinetics units
    Aunits = [
        '',  # Zeroth-order
        's^-1'.format(time_units),  # First-order
        '{0}^3/({1}*{2})'.format(volume_units, molecule_units, time_units),  # Second-order
        '{0}^6/({1}^2*{2})'.format(volume_units, molecule_units, time_units),  # Third-order
        '{0}^9/({1}^3*{2})'.format(volume_units, molecule_units, time_units),  # Fourth-order
    ]
    Eunits = energy_units

    kinetics_list = []
    comments_list = []
    kinetics = ''
    comments = ''

    line = f.readline()
    while line != '':

        line_starts_with_comment = line.lstrip().startswith('!') or line.lstrip().startswith('//')
        line, comment = remove_comment_from_line(line)
        line = line.strip()
        comment = comment.strip()

        if 'end' in line or 'END' in line:
            break

        # if 'rev' in line or 'REV' in line:
        # can no longer name reactants rev...
        #    line = f.readline()
        #    continue  # need to re-do the comment stripping!

        if '=' in line and not line_starts_with_comment:
            # Finish previous record
            kinetics_list.append(kinetics)
            comments_list.append(comments)
            kinetics = ''
            comments = ''

        if line: kinetics += line + '\n'
        if comment: comments += comment + '\n'

        line = f.readline()

    # Don't forget the last reaction!
    if kinetics.strip() != '':
        kinetics_list.append(kinetics)
        comments_list.append(comments)

    if len(kinetics_list) == 0 and len(comments_list) == 0:
        # No reactions found
        pass
    elif kinetics_list[0] == '' and comments_list[-1] == '':
        # True for Chemkin files generated from RMG-Py
        kinetics_list.pop(0)
        comments_list.pop(-1)
    elif kinetics_list[0] == '' and comments_list[0] == '':
        # True for Chemkin files generated from RMG-Java
        warnings.warn("RMG-Java loading is no longer supported and may be"
                      " removed in version 2.3.", DeprecationWarning)
        kinetics_list.pop(0)
        comments_list.pop(0)
    else:
        # In reality, comments can occur anywhere in the Chemkin
        # file (e.g. either or both of before and after the
        # reaction equation)
        # If we can't tell what semantics we are using, then just
        # throw the comments away
        # (This is better than failing to load the Chemkin file at
        # all, which would likely occur otherwise)
        if kinetics_list[0] == '':
            kinetics_list.pop(0)
        if len(kinetics_list) != len(comments_list):
            logging.warning("Discarding comments from Chemkin file because not sure which reaction they apply to")
            comments_list = ['' for kinetics in kinetics_list]

    reaction_list = []
    for kinetics, comments in zip(kinetics_list, comments_list):
        try:
            reaction = read_kinetics_entry(kinetics, species_dict, Aunits, Eunits)
            reaction = read_reaction_comments(reaction, comments, read=read_comments)
        except ChemkinError as e:
            if "Skip reaction!" in str(e):
                logging.warning("Skipping the reaction {0!r}".format(kinetics))
                continue
            else:
                raise
        reaction_list.append(reaction)

    return reaction_list

################################################################################


def save_html_file(path, read_comments=True):
    """
    Save an output HTML file from the contents of a RMG-Java output folder
    """
    warnings.warn("RMG-Java loading is no longer supported and may be"
                  " removed in version 2.3.", DeprecationWarning)
    from rmgpy.rmg.model import CoreEdgeReactionModel
    from rmgpy.rmg.output import save_output_html
    chemkin_path = os.path.join(path, 'chemkin', 'chem.inp')
    dictionary_path = os.path.join(path, 'RMG_Dictionary.txt')
    model = CoreEdgeReactionModel()
    model.core.species, model.core.reactions = load_chemkin_file(chemkin_path, dictionary_path,
                                                                 read_comments=read_comments)
    output_path = os.path.join(path, 'output.html')
    species_path = os.path.join(path, 'species')
    if not os.path.isdir(species_path):
        os.makedirs(species_path)
    save_output_html(output_path, model)

################################################################################


def get_species_identifier(species):
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
                return '{0}'.format(species.molecule[0].get_formula())
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
            name = '{0}({1:d})'.format(species.molecule[0].get_formula(), species.index)
            if len(name) <= 10:
                if 'obs' in label:
                    # For MBSampledReactor, keep observed species tag
                    return name + '_obs'
                else:
                    return name

        # As a last resort, just use the index
        if species.index >= 0:
            if 'X' in name:
                # Helpful to keep X in the names of all surface species.
                name = 'SX({0:d})'.format(species.index)
            else:
                name = 'S({0:d})'.format(species.index)
            if len(name) <= 10:
                if 'obs' in label:
                    # For MBSampledReactor, keep observed species tag
                    return name + '_obs'
                else:
                    return name

    # If we're here then we just can't come up with a valid Chemkin name
    # for this species, so raise an exception
    raise ChemkinError("Unable to determine valid Chemkin identifier for species {0}.".format(species))

################################################################################


def write_thermo_entry(species, element_counts=None, verbose=True):
    """
    Return a string representation of the NASA model readable by Chemkin.
    To use this method you must have exactly two NASA polynomials in your
    model, and you must use the seven-coefficient forms for each.
    """

    thermo = species.get_thermo_data()

    if not isinstance(thermo, NASA):
        raise ChemkinError('Cannot generate Chemkin string for species "{0}": '
                           'Thermodynamics data must be a NASA object.'.format(species))

    assert len(thermo.polynomials) == 2
    assert thermo.polynomials[0].Tmin.value_si < thermo.polynomials[1].Tmin.value_si
    assert thermo.polynomials[0].Tmax.value_si == thermo.polynomials[1].Tmin.value_si
    assert thermo.polynomials[0].cm2 == 0 and thermo.polynomials[0].cm1 == 0
    assert thermo.polynomials[1].cm2 == 0 and thermo.polynomials[1].cm1 == 0

    # Determine the number of each type of element in the molecule
    if element_counts is None:
        element_counts = get_element_count(species.molecule[0])

    # Sort the element_counts dictionary so that it's C, H, Al, B, Cl, D, etc.
    # if there's any C, else Al, B, Cl, D, H, if not. This is the "Hill" system
    # done by Molecule.get_formula
    if 'C' in element_counts:
        sorted_elements = sorted(element_counts, key = lambda e: {'C':'0','H':'1'}.get(e, e))
    else:
        sorted_elements = sorted(element_counts)
    element_counts = {e: element_counts[e] for e in sorted_elements}

    string = ''
    # Write thermo comments
    if verbose:
        if thermo.comment:
            for line in thermo.comment.strip().split("\n"):
                if len(line) > 150:
                    short_lines = textwrap.fill(line, 150).split("\n")
                    for short_line in short_lines:
                        string += "! {0}\n".format(short_line)
                else:
                    string += "! {0}\n".format(line)

    # Compile element count string
    extended_syntax = len(element_counts) > 4  # If there are more than 4 elements, use extended syntax
    elements = []
    for key, count in element_counts.items():
        if isinstance(key, tuple):
            symbol, isotope = key
            chemkin_name = get_element(symbol, isotope=isotope).chemkin_name
        else:
            chemkin_name = key
        if extended_syntax:
            # Create a list of alternating elements and counts
            elements.extend([chemkin_name, str(count)])
        else:
            # Create a list of 5-column wide formatted element counts, e.g. 'C  10'
            elements.append('{0!s:<2}{1:>3d}'.format(chemkin_name, count))
    if extended_syntax:
        # Use the new-style Chemkin syntax for the element counts
        # Place all elements in space delimited format on new line
        # This will only be recognized by Chemkin 4 or later
        elem_1 = ' ' * 20
        elem_2 = '&\n' + ' '.join(elements)
    else:
        # Use the original Chemkin syntax for the element counts
        # Place up to 4 elements in columns 24-43 of the first line
        # (don't use the space in columns 74-78 for the 5th element
        #  because nobody else does and Cantera can't read it)
        elem_1 = ''.join(elements)
        elem_2 = ''

    # Line 1
    string += '{ident:<16}        {elem_1:<20}G{Tmin:>10.3f}{Tint:>10.3f}{Tmax:>8.2f}      1{elem_2}\n'.format(
        ident=get_species_identifier(species),
        elem_1=elem_1,
        Tmin=thermo.polynomials[0].Tmin.value_si,
        Tint=thermo.polynomials[1].Tmax.value_si,
        Tmax=thermo.polynomials[0].Tmax.value_si,
        elem_2=elem_2,
    )

    # Line 2
    string += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}{4:< 15.8E}    2\n'.format(thermo.polynomials[1].c0,
                                                                                      thermo.polynomials[1].c1,
                                                                                      thermo.polynomials[1].c2,
                                                                                      thermo.polynomials[1].c3,
                                                                                      thermo.polynomials[1].c4)

    # Line 3
    string += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}{4:< 15.8E}    3\n'.format(thermo.polynomials[1].c5,
                                                                                      thermo.polynomials[1].c6,
                                                                                      thermo.polynomials[0].c0,
                                                                                      thermo.polynomials[0].c1,
                                                                                      thermo.polynomials[0].c2)

    # Line 4
    string += '{0:< 15.8E}{1:< 15.8E}{2:< 15.8E}{3:< 15.8E}                   4\n'.format(thermo.polynomials[0].c3,
                                                                                          thermo.polynomials[0].c4,
                                                                                          thermo.polynomials[0].c5,
                                                                                          thermo.polynomials[0].c6)

    return string

################################################################################


def write_reaction_string(reaction, java_library=False):
    """
    Return a reaction string in chemkin format.
    """
    kinetics = reaction.kinetics

    if kinetics is None:
        reaction_string = ' + '.join([get_species_identifier(reactant) for reactant in reaction.reactants])
        reaction_string += ' => ' if not reaction.reversible else ' = '
        reaction_string += ' + '.join([get_species_identifier(product) for product in reaction.products])
        return reaction_string

    if reaction.specific_collider is not None and not isinstance(kinetics, (_kinetics.Lindemann, _kinetics.Troe)):
        raise ChemkinError('A third body collider `(+species)` is only allowed in either the Troe or Lindemann formats '
                           'that support different reaction orders for the Low and High pressures limits. '
                           'You should revise reaction {0}'.format(reaction.label))

    if java_library:
        warnings.warn("Writing RMG-Java format is no longer supported and may be"
                      " removed in version 2.3.", DeprecationWarning)
        third_body = ''
        if kinetics.is_pressure_dependent():
            if (isinstance(kinetics, _kinetics.ThirdBody) and
                    not isinstance(kinetics, (_kinetics.Lindemann, _kinetics.Troe))):
                third_body = ' + M'
            elif isinstance(kinetics, _kinetics.PDepArrhenius):
                third_body = ''
            elif isinstance(kinetics, _kinetics.Chebyshev):
                third_body = ''
            else:
                third_body = ' (+{0})'.format(
                    get_species_identifier(reaction.specific_collider)) if reaction.specific_collider else ' (+M)'

        reaction_string = ' + '.join([get_species_identifier(reactant) for reactant in reaction.reactants])
        reaction_string += third_body
        reaction_string += ' = ' if reaction.reversible else ' => '
        reaction_string += ' + '.join([get_species_identifier(product) for product in reaction.products])
        reaction_string += third_body

    else:
        third_body = ''
        if kinetics.is_pressure_dependent():
            if (isinstance(kinetics, _kinetics.ThirdBody) and
                    not isinstance(kinetics, (_kinetics.Lindemann, _kinetics.Troe))):
                third_body = '+M'
            elif isinstance(kinetics, (_kinetics.PDepArrhenius, _kinetics.MultiPDepArrhenius)):
                third_body = ''
            else:
                third_body = '(+{0})'.format(
                    get_species_identifier(reaction.specific_collider)) if reaction.specific_collider else '(+M)'

        reaction_string = '+'.join([get_species_identifier(reactant) for reactant in reaction.reactants])
        reaction_string += third_body
        reaction_string += '=' if reaction.reversible else '=>'
        reaction_string += '+'.join([get_species_identifier(product) for product in reaction.products])
        reaction_string += third_body

    if len(reaction_string) > 52:
        logging.warning("Chemkin reaction string {0!r} is too long for Chemkin 2!".format(reaction_string))
    return reaction_string

################################################################################


def write_kinetics_entry(reaction, species_list, verbose=True, java_library=False, commented=False):
    """
    Return a string representation of the reaction as used in a Chemkin
    file. Use `verbose = True` to turn on kinetics comments.
    Use `commented = True` to comment out the entire reaction.
    Use java_library = True in order to generate a kinetics entry suitable
    for an RMG-Java kinetics library.
    """
    string = ""

    if isinstance(reaction.kinetics, (_kinetics.MultiArrhenius, _kinetics.MultiPDepArrhenius)):
        if verbose:
            if reaction.kinetics.comment:
                for line in reaction.kinetics.comment.split("\n"):
                    string += "! {0}\n".format(line)
        for kinetics in reaction.kinetics.arrhenius:
            if isinstance(reaction, LibraryReaction):
                new_reaction = LibraryReaction(index=reaction.index,
                                               reactants=reaction.reactants,
                                               products=reaction.products,
                                               specific_collider=reaction.specific_collider,
                                               reversible=reaction.reversible,
                                               kinetics=kinetics,
                                               library=reaction.library
                                               )
            else:
                new_reaction = Reaction(index=reaction.index,
                                        reactants=reaction.reactants,
                                        products=reaction.products,
                                        specific_collider=reaction.specific_collider,
                                        reversible=reaction.reversible,
                                        kinetics=kinetics)
            string += write_kinetics_entry(new_reaction, species_list, verbose, java_library, commented)
            string += "DUPLICATE\n"

        if commented:
            # add comments to the start of each line
            string = '! ' + string.replace('\n', '\n! ')
        return string + "\n"

    # Add to global chemkin reaction count if the kinetics is not a duplicate
    global _chemkin_reaction_count
    if _chemkin_reaction_count is not None:
        _chemkin_reaction_count += 1

    if verbose:
        # Next line of comment contains Chemkin and RMG indices
        if _chemkin_reaction_count is not None:
            string += "! Reaction index: Chemkin #{0:d}; RMG #{1:d}\n".format(_chemkin_reaction_count, reaction.index)

        # Next line of comment contains information about the type of reaction
        if isinstance(reaction, TemplateReaction):
            string += '! Template reaction: {0!s}\n'.format(reaction.family)
        elif isinstance(reaction, LibraryReaction):
            string += '! Library reaction: {0!s}\n'.format(reaction.library)
        elif isinstance(reaction, PDepReaction):
            string += '! PDep reaction: {0!s}\n'.format(reaction.network)
            if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
                # Print additional information about the pdep network's high-P limit reactions if in debug mode.
                for rxn in reaction.network.path_reactions:
                    if isinstance(rxn, LibraryReaction):
                        string += '! High-P limit: {0} (Library reaction: {1!s})\n'.format(rxn, rxn.library)
                    else:
                        string += '! High-P limit: {0} (Template reaction: {1!s})\n'.format(rxn, rxn.family)

        if reaction.specific_collider is not None:
            string += "! Specific third body collider: {0}\n".format(reaction.specific_collider.label)
        # Next line of comment contains information about the pairs of reaction
        pairs = []
        if reaction.pairs:
            for pair in reaction.pairs:
                pairs.append([get_species_identifier(spec) for spec in pair])
            string += "! Flux pairs: "

            for pair in pairs:
                string += pair[0] + ", " + pair[1] + "; "
            string += "\n"

        # Remaining lines of comments taken from reaction kinetics
        if reaction.kinetics.comment:
            for line in reaction.kinetics.comment.split("\n"):
                if len(line) > 150:
                    short_lines = textwrap.fill(line, 150).split("\n")
                    for short_line in short_lines:
                        string += "! {0}\n".format(short_line)
                else:
                    string += "! {0}\n".format(line)

    kinetics = reaction.kinetics
    num_reactants = len(reaction.reactants)
    reaction_string = write_reaction_string(reaction, java_library)

    string += '{0!s:<51} '.format(reaction_string)

    if isinstance(kinetics, _kinetics.StickingCoefficient):
        string += '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(
            kinetics.A.value_si / (kinetics.T0.value_si ** kinetics.n.value_si),
            kinetics.n.value_si,
            kinetics.Ea.value_si / 4184.
        )
        string += '\n    STICK'
    elif isinstance(kinetics, _kinetics.Arrhenius):
        conversion_factor = kinetics.A.get_conversion_factor_from_si_to_cm_mol_s()
        if not isinstance(kinetics, _kinetics.SurfaceArrhenius):
            assert 0.999 < conversion_factor / 1.0e6 ** (num_reactants - 1) < 1.001, """
              Gas phase reaction \n{}
              with kinetics \n{}
              with {} reactants was expected to have
              kinetics.A.get_conversion_factor_from_si_to_cm_mol_s() = {}
              but instead it is {}
              """.format(string, repr(kinetics), num_reactants, 1.0e6 ** (num_reactants - 1), conversion_factor)
            # debugging; for gas phase only
        string += '{0:<9.6e} {1:<9.3f} {2:<9.3f}'.format(
            kinetics.A.value_si / (kinetics.T0.value_si ** kinetics.n.value_si) * conversion_factor,
            kinetics.n.value_si,
            kinetics.Ea.value_si / 4184.
        )
    elif isinstance(kinetics, (_kinetics.Lindemann, _kinetics.Troe)):
        arrhenius = kinetics.arrheniusHigh
        conversion_factor = arrhenius.A.get_conversion_factor_from_si_to_cm_mol_s()
        assert 0.999 < conversion_factor / 1.0e6 ** (num_reactants - 1) < 1.001  # for gas phase only
        string += '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(
            arrhenius.A.value_si / (arrhenius.T0.value_si ** arrhenius.n.value_si) * conversion_factor,
            arrhenius.n.value_si,
            arrhenius.Ea.value_si / 4184.
        )
    elif isinstance(kinetics, _kinetics.ThirdBody):
        arrhenius = kinetics.arrheniusLow
        conversion_factor = arrhenius.A.get_conversion_factor_from_si_to_cm_mol_s()
        assert 0.999 < conversion_factor / 1.0e6 ** num_reactants < 1.001  # for gas phase only
        string += '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(
            arrhenius.A.value_si / (arrhenius.T0.value_si ** arrhenius.n.value_si) * conversion_factor,
            arrhenius.n.value_si,
            arrhenius.Ea.value_si / 4184.
        )
    elif hasattr(kinetics, 'highPlimit') and kinetics.highPlimit is not None:
        arrhenius = kinetics.highPlimit
        conversion_factor = arrhenius.A.get_conversion_factor_from_si_to_cm_mol_s()
        assert 0.999 < conversion_factor / 1.0e6 ** (num_reactants - 1) < 1.001  # for gas phase only
        string += '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(
            arrhenius.A.value_si / (arrhenius.T0.value_si ** arrhenius.n.value_si) * conversion_factor,
            arrhenius.n.value_si,
            arrhenius.Ea.value_si / 4184.
        )
    else:
        # Print dummy values that Chemkin parses but ignores
        string += '{0:<9.3e} {1:<9.3f} {2:<9.3f}'.format(1, 0, 0)

    if java_library:
        warnings.warn("RMG-Java libraries are no longer supported and may be"
                      " removed in version 2.3.", DeprecationWarning)
        # Assume uncertainties are zero (when parsing from chemkin), may need to adapt later
        string += '{0:<9.1f} {1:<9.1f} {2:<9.1f}'.format(0, 0, 0)

    string += '\n'

    if getattr(kinetics, 'coverage_dependence', None):
        # Write coverage dependence parameters for surface reactions
        for species, cov_params in kinetics.coverage_dependence.items():
            label = get_species_identifier(species)
            string += f'    COV / {label:<41} '
            string += f"{cov_params['a'].value:<9.3g} {cov_params['m'].value:<9.3g} {cov_params['E'].value_si/4184.:<9.3f} /\n"

    if isinstance(kinetics, (_kinetics.ThirdBody, _kinetics.Lindemann, _kinetics.Troe)):
        # Write collider efficiencies
        for collider, efficiency in sorted(list(kinetics.efficiencies.items()), key=lambda item: id(item[0])):
            for species in species_list:
                if any([collider.is_isomorphic(molecule) for molecule in species.molecule]):
                    string += '{0!s}/{1:<4.2f}/ '.format(get_species_identifier(species), efficiency)
                    break
        string += '\n'

        if isinstance(kinetics, (_kinetics.Lindemann, _kinetics.Troe)):
            # Write low-P kinetics
            arrhenius = kinetics.arrheniusLow
            conversion_factor = arrhenius.A.get_conversion_factor_from_si_to_cm_mol_s()
            assert 0.999 < conversion_factor / 1.0e6 ** num_reactants < 1.001  # for gas phase only
            string += '    LOW/ {0:<9.3e} {1:<9.3f} {2:<9.3f}/\n'.format(
                arrhenius.A.value_si / (arrhenius.T0.value_si ** arrhenius.n.value_si) * conversion_factor,
                arrhenius.n.value_si,
                arrhenius.Ea.value_si / 4184.
            )
            if isinstance(kinetics, _kinetics.Troe):
                # Write Troe parameters
                if kinetics.T2 is None:
                    string += '    TROE/ {0:<9.3e} {1:<9.3g} {2:<9.3g}/\n'.format(kinetics.alpha, kinetics.T3.value_si,
                                                                                  kinetics.T1.value_si)
                else:
                    string += '    TROE/ {0:<9.3e} {1:<9.3g} {2:<9.3g} {3:<9.3g}/\n'.format(kinetics.alpha,
                                                                                            kinetics.T3.value_si,
                                                                                            kinetics.T1.value_si,
                                                                                            kinetics.T2.value_si)
    elif isinstance(kinetics, _kinetics.PDepArrhenius):
        for P, arrhenius in zip(kinetics.pressures.value_si, kinetics.arrhenius):
            if isinstance(arrhenius, _kinetics.MultiArrhenius):
                for arrh in arrhenius.arrhenius:
                    conversion_factor = arrh.A.get_conversion_factor_from_si_to_cm_mol_s()
                    assert 0.999 < conversion_factor / 1.0e6 ** (num_reactants - 1) < 1.001  # for gas phase only
                    string += '    PLOG/ {0:<9.6f} {1:<9.3e} {2:<9.3f} {3:<9.3f}/\n'.format(
                        P / 101325.,
                        arrh.A.value_si / (arrh.T0.value_si ** arrh.n.value_si) * conversion_factor,
                        arrh.n.value_si,
                        arrh.Ea.value_si / 4184.
                    )
            else:
                conversion_factor = arrhenius.A.get_conversion_factor_from_si_to_cm_mol_s()
                assert 0.999 < conversion_factor / 1.0e6 ** (num_reactants - 1) < 1.001  # for gas phase only
                string += '    PLOG/ {0:<9.6f} {1:<9.3e} {2:<9.3f} {3:<9.3f}/\n'.format(
                    P / 101325.,
                    arrhenius.A.value_si / (arrhenius.T0.value_si ** arrhenius.n.value_si) * conversion_factor,
                    arrhenius.n.value_si,
                    arrhenius.Ea.value_si / 4184.
                )
    elif isinstance(kinetics, _kinetics.Chebyshev):
        string += '    TCHEB/ {0:<9.3f} {1:<9.3f}/\n'.format(kinetics.Tmin.value_si, kinetics.Tmax.value_si)
        string += '    PCHEB/ {0:<9.3f} {1:<9.3f}/\n'.format(kinetics.Pmin.value_si / 101325.,
                                                             kinetics.Pmax.value_si / 101325.)
        string += '    CHEB/ {0:d} {1:d}/\n'.format(kinetics.degreeT, kinetics.degreeP)
        # rwest: bypassing the Units.get_conversion_factor_from_si_to_cm_mol_s() because it's in log10 space?
        if kinetics.degreeP < 6:
            coeffs = kinetics.coeffs.value_si.copy()
            coeffs[0, 0] += 6 * (num_reactants - 1)
            for i in range(kinetics.degreeT):
                string += '    CHEB/'
                for j in range(kinetics.degreeP):
                    string += ' {0:<12.3e}'.format(coeffs[i, j])
                string += '/\n'
        else:
            coeffs = []
            for i in range(kinetics.degreeT):
                for j in range(kinetics.degreeP):
                    coeffs.append(kinetics.coeffs.value_si[i, j])
            coeffs[0] += 6 * (num_reactants - 1)
            for i in range(len(coeffs)):
                if i % 5 == 0: string += '    CHEB/'
                string += ' {0:<12.3e}'.format(coeffs[i])
                if i % 5 == 4: string += '/\n'

    if reaction.duplicate:
        string += 'DUPLICATE\n'

    if commented:
        # add comments to the start of each line
        string = '! ' + string.replace('\n', '\n! ')
    return string

################################################################################


def mark_duplicate_reaction(test_reaction, reaction_list):
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
        same_dir_match = (reaction1.reactants == reaction2.reactants and reaction1.products == reaction2.products)
        opposite_dir_match = (reaction1.products == reaction2.reactants and reaction1.reactants == reaction2.products)
        if (same_dir_match or opposite_dir_match) and (reaction1.specific_collider == reaction2.specific_collider):
            if reaction1.duplicate and reaction2.duplicate:
                if reaction1.kinetics.is_pressure_dependent() != reaction2.kinetics.is_pressure_dependent():
                    # Reactions with mixed pressure dependence do not need to be marked duplicate in Chemkin
                    logging.warning('Marked reaction {0} as not duplicate because of mixed pressure dependence '
                                    'for saving to Chemkin file.'.format(reaction1))
                    reaction1.duplicate = False
                    reaction2.duplicate = False
                elif opposite_dir_match and not reaction1.reversible and not reaction2.reversible:
                    # Irreversible reactions in opposite directions do not need to be marked duplicate in Chemkin
                    logging.warning('Marked reaction {0} as not duplicate because they are irreversible '
                                    'in opposite directions for saving to Chemkin file.'.format(reaction1))
                    reaction1.duplicate = False
                    reaction2.duplicate = False
            else:
                if (reaction1.kinetics.is_pressure_dependent() == reaction2.kinetics.is_pressure_dependent()
                        and ((reaction1.reversible and reaction2.reversible)
                             or (same_dir_match and not reaction1.reversible and not reaction2.reversible))):
                    # Only mark as duplicate if both reactions are pressure dependent or both are
                    # not pressure dependent. Also, they need to both be reversible or both be
                    # irreversible in the same direction.  Do not mark as duplicates otherwise.
                    logging.warning('Marked reaction {0} as duplicate of {1} for saving '
                                    'to Chemkin file.'.format(reaction1, reaction2))
                    reaction1.duplicate = True
                    reaction2.duplicate = True


def mark_duplicate_reactions(reactions):
    """
    For a given list of `reactions`, mark all of the duplicate reactions as
    understood by Chemkin.
    
    This is pretty slow (quadratic in size of reactions list) so only call it if you're really worried
    you may have undetected duplicate reactions.
    """
    for index1 in range(len(reactions)):
        reaction1 = reactions[index1]
        remaining_list = reactions[index1 + 1:]
        mark_duplicate_reaction(reaction1, remaining_list)


def save_species_dictionary(path, species, old_style=False):
    """
    Save the given list of `species` as adjacency lists in a text file `path` 
    on disk.
    
    If `old_style==True` then it saves it in the old RMG-Java syntax.
    """
    with open(path, 'w') as f:
        for spec in species:
            if old_style:
                try:
                    f.write(spec.molecule[0].to_adjacency_list(label=get_species_identifier(spec),
                                                               remove_h=True, old_style=True))
                except:
                    new_adjlist = spec.molecule[0].to_adjacency_list(label=get_species_identifier(spec), remove_h=False)
                    f.write("// Couldn't save {0} in old RMG-Java syntax, but here it is in "
                            "newer RMG-Py syntax:".format(get_species_identifier(spec)))
                    f.write("\n// " + "\n// ".join(new_adjlist.splitlines()) + '\n')
            else:
                try:
                    for mol in spec.molecule:
                        if mol.reactive:
                            f.write(mol.to_adjacency_list(label=get_species_identifier(spec), remove_h=False))
                            break
                    else:
                        raise ValueError('No reactive structures were found for species '
                                         '{0}.'.format(get_species_identifier(spec)))
                except:
                    raise ChemkinError('Ran into error saving dictionary for species {0}. '
                                       'Please check your files.'.format(get_species_identifier(spec)))
            f.write('\n')


def save_transport_file(path, species):
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
        f.write("! {0:15} {1:8} {2:9} {3:9} {4:9} {5:9} {6:9} {7:9}\n".format(
            'Species', 'Shape', 'LJ-depth', 'LJ-diam', 'DiplMom', 'Polzblty', 'RotRelaxNum', 'Data'))
        f.write("! {0:15} {1:8} {2:9} {3:9} {4:9} {5:9} {6:9} {7:9}\n".format(
            'Name', 'Index', 'epsilon/k_B', 'sigma', 'mu', 'alpha', 'Zrot', 'Source'))
        for spec in species:
            transport_data = spec.get_transport_data()
            if not transport_data:
                missing_data = True
            else:
                missing_data = False

            label = get_species_identifier(spec)

            if missing_data:
                f.write('! {0:19s} {1!r}\n'.format(label, transport_data))
            else:
                f.write('{0:19} {1:d}   {2:9.3f} {3:9.3f} {4:9.3f} {5:9.3f} {6:9.3f}    ! {7:s}\n'.format(
                    label,
                    transport_data.shapeIndex,
                    transport_data.epsilon.value_si / constants.R,
                    transport_data.sigma.value_si * 1e10,
                    (transport_data.dipoleMoment.value_si * constants.c * 1e21 if transport_data.dipoleMoment else 0),
                    (transport_data.polarizability.value_si * 1e30 if transport_data.polarizability else 0),
                    (transport_data.rotrelaxcollnum if transport_data.rotrelaxcollnum else 0),
                    transport_data.comment,
                ))


def save_chemkin_file(path, species, reactions, verbose=True, check_for_duplicates=True):
    """
    Save a Chemkin input file to `path` on disk containing the provided lists
    of `species` and `reactions`.
    If check_for_duplicates is False then we don't check for unlabeled duplicate reactions,
    thus saving time (eg. if you are sure you've already labeled them as duplicate).
    """
    # Check for duplicate
    if check_for_duplicates:
        mark_duplicate_reactions(reactions)

    f = open(path, 'w')

    sorted_species = sorted(species, key=lambda species: species.index)

    # Elements section
    write_elements_section(f)

    # Species section
    f.write('SPECIES\n')
    for spec in sorted_species:
        label = get_species_identifier(spec)
        if verbose:
            f.write('    {0!s:<16}    ! {1}\n'.format(label, str(spec)))
        else:
            f.write('    {0!s:<16}\n'.format(label))
    f.write('END\n\n\n\n')

    # Thermodynamics section
    f.write('THERM ALL\n')
    f.write('   300.000  1000.000  5000.000\n\n')
    for spec in sorted_species:
        f.write(write_thermo_entry(spec, verbose=verbose))
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
    global _chemkin_reaction_count
    _chemkin_reaction_count = 0
    for rxn in reactions:
        f.write(write_kinetics_entry(rxn, species_list=species, verbose=verbose))
        # Don't forget to mark duplicates!
        f.write('\n')
    f.write('END\n\n')
    f.close()
    logging.info("Chemkin file contains {0} reactions.".format(_chemkin_reaction_count))
    _chemkin_reaction_count = None


def save_chemkin_surface_file(path, species, reactions, verbose=True, check_for_duplicates=True,
                              surface_site_density=None):
    """
    Save a Chemkin *surface* input file to `path` on disk containing the provided lists
    of `species` and `reactions`.
    If check_for_duplicates is False then we don't check for unlabeled duplicate reactions,
    thus saving time (eg. if you are sure you've already labeled them as duplicate).
    """
    # Check for duplicate
    if check_for_duplicates:
        mark_duplicate_reactions(reactions)

    f = open(path, 'w')

    sorted_species = sorted(species, key=lambda species: species.index)

    # Species section
    surface_name = None
    if surface_name:
        f.write('SITE/{}/'.format(surface_name))
    else:
        f.write('SITE ')
    if surface_site_density:
        f.write('  SDEN/{0:.4E}/ ! mol/cm^2\n'.format(surface_site_density.value_si * 1e-4))
    else:
        f.write('  SDEN/2.72E-9/ ! mol/cm^2 DEFAULT!\n')
    # todo: add surface site density from reactor simulation
    for spec in sorted_species:
        label = get_species_identifier(spec)
        number_of_sites = spec.molecule[0].get_num_atoms('X')
        if  number_of_sites >= 2:
            label +=  f"/{number_of_sites}/"
        if verbose:
            f.write('    {0!s:<16}    ! {1}\n'.format(label, str(spec)))
        else:
            f.write('    {0!s}\n'.format(label))
    f.write('END\n\n\n\n')

    # Thermodynamics section
    f.write('THERM ALL\n')
    f.write('    300.000  1000.000  5000.000\n\n')
    for spec in sorted_species:
        f.write(write_thermo_entry(spec, verbose=verbose))
        f.write('\n')
    f.write('END\n\n\n\n')

    # Reactions section
    f.write('REACTIONS    KCAL/MOLE   MOLES\n\n')
    global _chemkin_reaction_count
    _chemkin_reaction_count = 0
    for rxn in reactions:
        f.write(write_kinetics_entry(rxn, species_list=species, verbose=verbose))
        f.write('\n')
    f.write('END\n\n')
    f.close()
    logging.info("Chemkin file contains {0} reactions.".format(_chemkin_reaction_count))
    _chemkin_reaction_count = None


def save_java_kinetics_library(path, species, reactions):
    """
    Save the reaction files for a RMG-Java kinetics library: pdepreactions.txt
    and reactions.txt given a list of reactions, with species.txt containing the
    RMG-Java formatted dictionary.
    """
    warnings.warn("Java kinetics libararies are no longer supported and may be" \
                  "removed in version 2.3.", DeprecationWarning)
    # Check for duplicate
    mark_duplicate_reactions(reactions)

    f = open(os.path.join(path, 'reactions.txt'), 'w')
    f2 = open(os.path.join(path, 'pdepreactions.txt'), 'w')

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
        if rxn.kinetics.is_pressure_dependent():
            f2.write(write_kinetics_entry(rxn, species_list=species, verbose=False, java_library=True))
            f2.write('\n')
        else:
            f.write(write_kinetics_entry(rxn, species_list=species, verbose=False, java_library=True))
            f.write('\n')
    f.close()
    f2.close()

    save_species_dictionary(os.path.join(path, 'species.txt'), species, old_style=True)


def save_chemkin(reaction_model, path, verbose_path, dictionary_path=None, transport_path=None, 
                 save_edge_species=False):
    """
    Save a Chemkin file for the current model as well as any desired output
    species and reactions to `path`. If `save_edge_species` is True, then 
    a chemkin file and dictionary file for the core AND edge species and reactions
    will be saved.  It also saves verbose versions of each file.
    """
    if save_edge_species:
        species_list = reaction_model.core.species + reaction_model.edge.species
        rxn_list = reaction_model.core.reactions + reaction_model.edge.reactions
    else:
        species_list = reaction_model.core.species + reaction_model.output_species_list
        rxn_list = reaction_model.core.reactions + reaction_model.output_reaction_list

    if any([s.contains_surface_site() for s in reaction_model.core.species]):
        # it's a surface model
        root, ext = os.path.splitext(path)
        gas_path = root + '-gas' + ext
        surface_path = root + '-surface' + ext
        root, ext = os.path.splitext(verbose_path)
        gas_verbose_path = root + '-gas' + ext
        surface_verbose_path = root + '-surface' + ext

        surface_species_list = []
        gas_species_list = []
        surface_rxn_list = []
        gas_rxn_list = []

        for s in species_list:
            if s.contains_surface_site():
                surface_species_list.append(s)
            else:
                gas_species_list.append(s)
        for r in rxn_list:
            if r.is_surface_reaction():
                surface_rxn_list.append(r)
            else:
                gas_rxn_list.append(r)

        # We should already have marked everything as duplicates by now so use check_for_duplicates=False
        save_chemkin_file(gas_path, gas_species_list, gas_rxn_list, verbose=False, check_for_duplicates=False)
        save_chemkin_surface_file(surface_path, surface_species_list, surface_rxn_list, verbose=False,
                                  check_for_duplicates=False, surface_site_density=reaction_model.surface_site_density)
        logging.info('Saving annotated version of Chemkin files...')
        save_chemkin_file(gas_verbose_path, gas_species_list, gas_rxn_list, verbose=True, check_for_duplicates=False)
        save_chemkin_surface_file(surface_verbose_path, surface_species_list, surface_rxn_list, verbose=True,
                                  check_for_duplicates=False, surface_site_density=reaction_model.surface_site_density)

    else:
        # Gas phase only
        save_chemkin_file(path, species_list, rxn_list, verbose=False, check_for_duplicates=False)
        logging.info('Saving annotated version of Chemkin file...')
        save_chemkin_file(verbose_path, species_list, rxn_list, verbose=True, check_for_duplicates=False)
    if dictionary_path:
        save_species_dictionary(dictionary_path, species_list)
    if transport_path:
        save_transport_file(transport_path, species_list)


def save_chemkin_files(rmg):
    """
    Save the current reaction model to a set of Chemkin files.
    """

    # todo: make this an attribute or method of reactionModel
    is_surface_model = any([s.contains_surface_site() for s in rmg.reaction_model.core.species])

    logging.info('Saving current model core to Chemkin file...')
    this_chemkin_path = os.path.join(rmg.output_directory, 'chemkin',
                                     'chem{0:04d}.inp'.format(len(rmg.reaction_model.core.species)))
    latest_chemkin_path = os.path.join(rmg.output_directory, 'chemkin', 'chem.inp')
    latest_chemkin_verbose_path = os.path.join(rmg.output_directory, 'chemkin', 'chem_annotated.inp')
    latest_dictionary_path = os.path.join(rmg.output_directory, 'chemkin', 'species_dictionary.txt')
    latest_transport_path = os.path.join(rmg.output_directory, 'chemkin', 'tran.dat')
    save_chemkin(rmg.reaction_model,
                 this_chemkin_path,
                 latest_chemkin_verbose_path,
                 latest_dictionary_path,
                 latest_transport_path,
                 save_edge_species=False)

    if is_surface_model:
        paths = []
        for phase in ['surface', 'gas']:
            root, ext = os.path.splitext(this_chemkin_path)
            path1 = root + '-' + phase + ext
            root, ext = os.path.splitext(latest_chemkin_path)
            path2 = root + '-' + phase + ext
            paths.append((path1, path2))
    else:
        paths = [(this_chemkin_path, latest_chemkin_path)]

    for this_chemkin_path, latest_chemkin_path in paths:
        if os.path.exists(latest_chemkin_path):
            os.unlink(latest_chemkin_path)
        shutil.copy2(this_chemkin_path, latest_chemkin_path)

    if rmg.save_edge_species:
        logging.info('Saving current model core and edge to Chemkin file...')
        this_chemkin_path = os.path.join(rmg.output_directory, 'chemkin',
                                         'chem_edge{0:04d}.inp'.format(len(rmg.reaction_model.core.species)))
        latest_chemkin_path = os.path.join(rmg.output_directory, 'chemkin', 'chem_edge.inp')
        latest_chemkin_verbose_path = os.path.join(rmg.output_directory, 'chemkin', 'chem_edge_annotated.inp')
        latest_dictionary_path = os.path.join(rmg.output_directory, 'chemkin', 'species_edge_dictionary.txt')
        latest_transport_path = None
        save_chemkin(rmg.reaction_model, this_chemkin_path, latest_chemkin_verbose_path, latest_dictionary_path,
                     latest_transport_path, rmg.save_edge_species)

        if is_surface_model:
            paths = []
            for phase in ['surface', 'gas']:
                root, ext = os.path.splitext(this_chemkin_path)
                path1 = root + '-' + phase + ext
                root, ext = os.path.splitext(latest_chemkin_path)
                path2 = root + '-' + phase + ext
                paths.append((path1, path2))
        else:
            paths = [(this_chemkin_path, latest_chemkin_path)]

        for this_chemkin_path, latest_chemkin_path in paths:
            if os.path.exists(latest_chemkin_path):
                os.unlink(latest_chemkin_path)
            shutil.copy2(this_chemkin_path, latest_chemkin_path)


def write_elements_section(f):
    """
    Write the ELEMENTS section of the chemkin file.  This file currently lists
    all elements and isotopes available in RMG. It may become useful in the future
    to only include elements/isotopes present in the current RMG run. 
    """

    s = 'ELEMENTS\n'

    # map of isotope elements with chemkin-compatible element representation:
    elements = ('H', ('H', 2), ('H', 3), 'C', ('C', 13), 'O', ('O', 18), 'N', 'Ne', 'Ar', 'He', 'Si', 'S',
                'F', 'Cl', 'Br', 'I')
    for el in elements:
        if isinstance(el, tuple):
            symbol, isotope = el
            chemkin_name = get_element(symbol, isotope=isotope).chemkin_name
            mass = 1000 * get_element(symbol, isotope=isotope).mass
            s += '\t{0} /{1:.3f}/\n'.format(chemkin_name, mass)
        else:
            s += '\t' + el + '\n'
    s += '\tX /195.083/\n'
    s += 'END\n\n'

    f.write(s)


class ChemkinWriter(object):
    """
    This class listens to a RMG subject
    and writes a chemkin file with the current state of the RMG model,
    to a chemkin subfolder.


    A new instance of the class can be appended to a subject as follows:
    
    rmg = ...
    listener = ChemkinWriter(outputDirectory)
    rmg.attach(listener)

    Whenever the subject calls the .notify() method, the
    .update() method of the listener will be called.

    To stop listening to the subject, the class can be detached
    from its subject:

    rmg.detach(listener)
    
    """
    def __init__(self, output_directory=''):
        super(ChemkinWriter, self).__init__()
        make_output_subdirectory(output_directory, 'chemkin')

    def update(self, rmg):
        save_chemkin_files(rmg)
