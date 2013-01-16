#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
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
This module contains functionality for working with kinetics families.
"""

import os
import os.path
import re
import logging
import codecs
from copy import copy, deepcopy

from rmgpy.data.base import *

from rmgpy.quantity import Quantity
from rmgpy.reaction import Reaction, ReactionError
from rmgpy.kinetics import Arrhenius, ArrheniusEP, ThirdBody, Lindemann, Troe, \
                           PDepArrhenius, MultiArrhenius, MultiPDepArrhenius, \
                           Chebyshev, KineticsData, PDepKineticsModel
from rmgpy.molecule import Bond, GroupBond, Group
from rmgpy.species import Species
from .common import KineticsError, saveEntry

################################################################################

class LibraryReaction(Reaction):
    """
    A Reaction object generated from a reaction library. In addition to the
    usual attributes, this class includes `library` and `entry` attributes to
    store the library and the entry in that library that it was created from.
    """
    
    def __init__(self,
                 index=-1,
                 reactants=None,
                 products=None,
                 kinetics=None,
                 reversible=True,
                 transitionState=None,
                 duplicate=False,
                 degeneracy=1,
                 pairs=None,
                 library=None,
                 entry=None
                 ):
        Reaction.__init__(self,
                          index=index,
                          reactants=reactants,
                          products=products,
                          kinetics=kinetics,
                          reversible=reversible,
                          transitionState=transitionState,
                          duplicate=duplicate,
                          degeneracy=degeneracy,
                          pairs=pairs
                          )
        self.library = library
        self.family = library
        self.entry = entry

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (LibraryReaction, (self.index,
                                  self.reactants,
                                  self.products,
                                  self.kinetics,
                                  self.reversible,
                                  self.transitionState,
                                  self.duplicate,
                                  self.degeneracy,
                                  self.pairs,
                                  self.library,
                                  self.entry
                                  ))

    def getSource(self):
        """
        Return the database that was the source of this reaction. For a
        LibraryReaction this should be a KineticsLibrary object.
        """
        return self.library

################################################################################

class KineticsLibrary(Database):
    """
    A class for working with an RMG kinetics library.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def __repr__(self):
        return '<KineticsLibrary "{0}">'.format(self.label)

    def getSpecies(self):
        """
        Return a dictionary containing all of the species in this kinetics
        library.
        """
        speciesDict = {}
        
        def speciesMatch(speciesA, speciesB):
            for moleculeA in speciesA.molecule:
                for moleculeB in speciesB.molecule:
                    if moleculeA.isIsomorphic(moleculeB):
                        return True
            return False
        
        entries = self.entries.values()
        for entry in entries:
            for reactant in entry.item.reactants:
                if reactant.label not in speciesDict:
                    speciesDict[reactant.label] = reactant
                elif not speciesMatch(reactant, speciesDict[reactant.label]):
                    print reactant.molecule[0].toAdjacencyList()
                    print speciesDict[reactant.label].molecule[0].toAdjacencyList()
                    raise DatabaseError('Species label "{0}" used for multiple species in kinetics library {1}.'.format(reactant.label, self.label))
            for product in entry.item.products:
                if product.label not in speciesDict:
                    speciesDict[product.label] = product
                elif not speciesMatch(product, speciesDict[product.label]):
                    import pdb; pdb.set_trace()
                    print product.molecule[0].toAdjacencyList()
                    print speciesDict[product.label].molecule[0].toAdjacencyList()
                    print product.molecule[0].isIsomorphic(speciesDict[product.label].molecule[0])
                    raise DatabaseError('Species label "{0}" used for multiple species in kinetics library {1}.'.format(product.label, self.label))
        
        return speciesDict

    def markValidDuplicates(self, reactions1, reactions2):
        """
        Check for reactions that appear in both lists,
        and mark them as (valid) duplicates.
        """
        for r1 in reactions1:
            for r2 in reactions2:
                if (r1.reactants == r2.reactants and
                    r1.products == r2.products and
                    r1.reversible == r2.reversible
                    ):
                    r1.duplicate = True
                    r2.duplicate = True

    def checkForDuplicates(self):
        """
        Check that all duplicate reactions in the kinetics library are
        properly marked (i.e. with their ``duplicate`` attribute set to 
        ``True``).
        """
        for entry0 in self.entries.values():
            reaction0 = entry0.item
            if not reaction0.duplicate:
                # This reaction is not marked as a duplicate reaction
                # This means that if we find any duplicate reactions, it is an error
                for entry in self.entries.values():
                    reaction = entry.item
                    if (reaction0 is not reaction and
                        reaction0.reactants == reaction.reactants and
                        reaction0.products == reaction.products and
                        reaction0.reversible == reaction.reversible
                        ):
                        # We found a duplicate reaction that wasn't marked!
                        raise DatabaseError('Unexpected duplicate reaction {0} in kinetics library {1}.'.format(reaction0, self.label))                   

    def convertDuplicatesToMulti(self):
        """
        Merge all marked duplicate reactions in the kinetics library
        into single reactions with multiple kinetics.
        """
        print "trying to find duplicates"
        entries_to_remove = []
        for entry0 in self.entries.values():
            if entry0 in entries_to_remove:
                continue
            reaction0 = entry0.item
            if not reaction0.duplicate:
                continue
            print "Found a duplicate reaction: {0}".format(reaction0)
            duplicates = [entry0]
            for entry in self.entries.values():
                reaction = entry.item
                if reaction0 is reaction:
                    continue
                if reaction0.isIsomorphic(reaction, eitherDirection=False):
                    if reaction0.reversible != reaction.reversible:
                        print "Reactions isomorphic but with different reversibilities"
                        continue
                    duplicates.append(entry)
            
            assert len(duplicates)>1
            kineticsList = []
            longDesc = ''
            
            for entry in duplicates:
                kinetics = entry.data
                kineticsList.append(kinetics)
                Tmin = kinetics.Tmin
                Tmax = kinetics.Tmax
                if kinetics.isPressureDependent():
                    Pmin = kinetics.Pmin
                    Pmax = kinetics.Pmax
                else:
                    Pmin = None
                    Pmax = None
                longDesc += entry.longDesc+'\n'
            
            if len(kineticsList) == 2 and isinstance(kineticsList[0],ThirdBody) and not isinstance(kineticsList[1],ThirdBody):
                continue
            elif len(kineticsList) == 2 and isinstance(kineticsList[1],ThirdBody) and not isinstance(kineticsList[0],ThirdBody):
                continue
                    
            
            if all([isinstance(k, Arrhenius) for k in kineticsList]):
                entry0.data = MultiArrhenius(arrhenius=kineticsList, Tmin=Tmin, Tmax=Tmax)
            elif all([isinstance(k, PDepArrhenius) for k in kineticsList]):
                entry0.data = MultiPDepArrhenius(arrhenius=kineticsList, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax)
            else:
                raise Exception('Only Arrhenius and PDepArrhenius kinetics supported for duplicate reactions.')
            entry0.longDesc = longDesc
            entries_to_remove.extend(duplicates[1:])
        for entry in entries_to_remove:
            print "removing duplicate reaction with index {0}.".format(entry.index)
            del(self.entries[entry.index])
        print "NB. the entries have not been renumbered, so these indices are missing."
        
        
    def load(self, path, local_context=None, global_context=None):
        Database.load(self, path, local_context, global_context)
        
        # Generate a unique set of the species in the kinetics library
        speciesDict = self.getSpecies()
        # Make sure all of the reactions draw from only this set
        entries = self.entries.values()
        for entry in entries:
            entry.item.reactants = [speciesDict[spec.label] for spec in entry.item.reactants]
            entry.item.products = [speciesDict[spec.label] for spec in entry.item.products]
            
        self.checkForDuplicates()
        
    def loadEntry(self,
                  index,
                  reactant1,
                  product1,
                  kinetics,
                  reactant2=None,
                  reactant3=None,
                  product2=None,
                  product3=None,
                  degeneracy=1,
                  label='',
                  duplicate=False,
                  reversible=True,
                  reference=None,
                  referenceType='',
                  shortDesc='',
                  longDesc='',
                  history=None
                  ):
        
        reactants = [Species(label=reactant1.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(reactant1)])]
        if reactant2 is not None: reactants.append(Species(label=reactant2.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(reactant2)]))
        if reactant3 is not None: reactants.append(Species(label=reactant3.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(reactant3)]))

        products = [Species(label=product1.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(product1)])]
        if product2 is not None: products.append(Species(label=product2.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(product2)]))
        if product3 is not None: products.append(Species(label=product3.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(product3)]))
        
        comment = "Reaction and kinetics from {0}.".format(self.label)
        if shortDesc.strip(): 
            comment += "{0!s}\n".format(shortDesc.strip())
        if longDesc.strip():
            comment += str(re.sub('\s*\n\s*','\n',longDesc))
        kinetics.comment = comment.strip()
        
        self.entries['{0:d}:{1}'.format(index,label)] = Entry(
            index = index,
            label = label,
            item = Reaction(reactants=reactants, products=products, degeneracy=degeneracy, duplicate=duplicate, reversible=reversible),
            data = kinetics,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            history = history or [],
        )
        
        # Convert SMILES to Molecule objects in collision efficiencies
        if isinstance(kinetics, PDepKineticsModel):
            efficiencies = {}
            for smiles, eff in kinetics.efficiencies.items():
                if isinstance(smiles, str):
                    efficiencies[Molecule().fromSMILES(smiles)] = eff
            kinetics.efficiencies = efficiencies

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the kinetics library to the file object `f`.
        """
        return saveEntry(f, entry)

    def loadOld(self, path):
        """
        Load an old-style RMG kinetics library from the location `path`.
        """
        path = os.path.abspath(path)

        self.loadOldDictionary(os.path.join(path,'species.txt'), pattern=False)
        species = dict([(label, Species(label=label, molecule=[entry.item])) for label, entry in self.entries.iteritems()])
        
        reactions = []
        reactions.extend(self.__loadOldReactions(os.path.join(path,'reactions.txt'), species))
        if os.path.exists(os.path.join(path,'pdepreactions.txt')):
            pdep_reactions = self.__loadOldReactions(os.path.join(path,'pdepreactions.txt'), species)
            # RMG-Py likes otherwise equivalent PDep and non-pdep reactions to be marked as duplicates
            self.markValidDuplicates(reactions, pdep_reactions)
            reactions.extend(pdep_reactions)

        self.entries = {}
        for index, reaction in enumerate(reactions):
            entry = Entry(
                index = index+1,
                item = reaction,
                data = reaction.kinetics,
            )
            entry.longDesc = reaction.kinetics.comment
            reaction.kinetics.comment = ''
            self.entries[index+1] = entry
            reaction.kinetics = None
        
        self.checkForDuplicates()
        self.convertDuplicatesToMulti()

    def __loadOldReactions(self, path, species):
        """
        Load an old-style reaction library from `path`. This algorithm can
        handle both the pressure-independent and pressure-dependent reaction
        files. If the pressure-dependent file is read, the extra pressure-
        dependent kinetics information is ignored unless the kinetics database
        is a seed mechanism.
        """
        reactions = []
        
        # Process the reactions or pdepreactions file
        try:
            inUnitSection = False; inReactionSection = False
            Aunits = []; Eunits = ''
            reaction = None; kinetics = None
            next_reaction_comment = ''
            factorSI = 1.0
            
            fdict = open(path, 'r')
            for line in fdict:
                line, comment = splitLineAndComment(line)
                line = line.strip()
                if len(line) == 0:
                    comment = comment.strip()
                    # collect all comment lines and assume they're for the following reaction 
                    next_reaction_comment += comment + '\n'
                    continue
                else: # line is not empty
                    if inUnitSection:
                        if 'A:' in line or 'E:' in line:
                            units = line.split()[1]
                            if 'A:' in line:
                                Aunits0 = units.split('/') # Assume this is a 3-tuple: moles or molecules, volume, time
                                Aunits0[1] = Aunits0[1][0:-1] # Remove '3' from e.g. 'm3' or 'cm3'; this is assumed
                                Aunits = [
                                    '',                                                         # Zeroth-order
                                    '{0}^-1'.format(Aunits0[2]),                                     # First-order
                                    '{0}^3/({1}*{2})'.format(Aunits0[1], Aunits0[0], Aunits0[2]),      # Second-order
                                    '{0}^6/({1}^2*{2})'.format(Aunits0[1], Aunits0[0], Aunits0[2]),    # Third-order
                                ]
                                assert Aunits0[1] in ['cm', 'm']
                                if Aunits0[1] == 'cm':
                                    factorSI = 1e6
                            elif 'E:' in line:
                                Eunits = units
                    elif inReactionSection:
                        if '=' in line:
                            # This line contains a reaction equation, (high-P) Arrhenius parameters, and uncertainties

                            # Strip out "(+M)" from line
                            line = line.replace("(+M)", "")
                            line = line.replace("(+m)", "")

                            items = line.split()

                            # Find the reaction arrow
                            for arrow in ['<=>', '=>', '=', '->']:
                                if arrow in items:
                                    arrowIndex = items.index(arrow)
                                    break

                            # Find the start of the data
                            try:
                                temp = float(items[-6])
                                dataIndex = -6
                            except ValueError:
                                dataIndex = -3

                            # Find the reactant and product items
                            hasThirdBody = False
                            reactantItems = []
                            for item in items[0:arrowIndex]:
                                if item != '+':
                                    for i in item.split('+'):
                                        if i != '' and i != 'M' and i != 'm': reactantItems.append(i)
                                        elif i != '' and (i == 'M' or i == 'm'): hasThirdBody = True
                            productItems = []
                            for item in items[arrowIndex+1:dataIndex]:
                                if item != '+':
                                    for i in item.split('+'):
                                        if i != '' and i != 'M' and i != 'm': productItems.append(i)
                                        elif i != '' and (i == 'M' or i == 'm'): hasThirdBody = True
                            
                            reactants = []; products = []
                            for item in reactantItems:
                                try:
                                    reactants.append(species[item])
                                except KeyError:
                                    raise DatabaseError('Reactant {0} not found in species dictionary.'.format(item))
                            for item in productItems:
                                try:
                                    products.append(species[item])
                                except KeyError:
                                    raise DatabaseError('Product {0} not found in species dictionary.'.format(item))

                            if dataIndex == -6:
                                A, n, Ea, dA, dn, dEa = items[-6:]
                                A = float(A)
                            else:
                                A, n, Ea = items[-3:]
                                dA = '0'; dn = '0'; dEa = '0'
                            
                            A = float(A)
                            kunits = Aunits[len(reactants)+1] if hasThirdBody else Aunits[len(reactants)]
                            if dA[0] == '*':
                                A = Quantity(A,kunits,'*|/',float(dA[1:]))
                            else:
                                dA = float(dA)
                                if dA != 0:
                                    A = Quantity(A,kunits,'+|-',dA)
                                else:
                                    A = Quantity(A,kunits)

                            n = float(n); dn = float(dn)
                            if dn != 0:
                                n = Quantity(n,'','+|-',dn)
                            else:
                                n = Quantity(n,'')

                            Ea = float(Ea); dEa = float(dEa)
                            if dEa != 0:
                                Ea = Quantity(Ea,Eunits,'+|-',dEa)
                            else:
                                Ea = Quantity(Ea,Eunits)

                            kinetics = Arrhenius(A=A, n=n, Ea=Ea, T0=(1.0,"K"))
                            if hasThirdBody:
                                kinetics = ThirdBody(arrheniusLow=kinetics)
                            
                            reaction = Reaction(
                                reactants=reactants,
                                products=products,
                                kinetics=kinetics,
                                reversible=(arrow in ['<=>', '=']),
                            )
                            reaction.kinetics.comment = next_reaction_comment
                            next_reaction_comment = ""
                            reactions.append(reaction)

                        elif 'PLOG' in line:
                            # This line contains pressure-dependent Arrhenius parameters in Chemkin format
                            items = line.split('/')
                            P, A, n, Ea = items[1].split()
                            P = float(P)
                            arrhenius = Arrhenius(
                                A = (float(A), Aunits[len(reactants)]), 
                                n = float(n), 
                                Ea = (float(Ea), Eunits), 
                                T0 = (1.0,"K"),
                            )
                            
                            if not isinstance(kinetics, PDepArrhenius):
                                old_kinetics = kinetics
                                comment = old_kinetics.comment
                                old_kinetics.comment = ''
                                assert isinstance(old_kinetics, Arrhenius)
                                kinetics = PDepArrhenius(pressures=([P],"atm"), arrhenius=[arrhenius], highPlimit=old_kinetics, comment=comment)
                            else:
                                pressures = list(kinetics.pressures.value_si)
                                pressures.append(P*101325.)
                                kinetics.pressures.value_si = numpy.array(pressures, numpy.float)
                                kinetics.arrhenius.append(arrhenius)
                            reaction.kinetics = kinetics
                            
                        elif 'CHEB' in line or 'cheb' in line:
                            # Chebyshev parameters
                            if not isinstance(kinetics, Chebyshev):
                                kinetics = Chebyshev()
                                kinetics.kunits = Aunits[len(reactants)]
                                reaction.kinetics = kinetics
                                chebyshevCoeffs = []
                            tokens = [t.strip() for t in line.split('/')]
                            if 'TCHEB' in line:
                                index = tokens.index('TCHEB')
                                tokens2 = tokens[index+1].split()
                                kinetics.Tmin = Quantity(float(tokens2[0].strip()),"K")
                                kinetics.Tmax = Quantity(float(tokens2[1].strip()),"K")
                            if 'PCHEB' in line:
                                index = tokens.index('PCHEB')
                                tokens2 = tokens[index+1].split()
                                kinetics.Pmin = Quantity(float(tokens2[0].strip()),"atm")
                                kinetics.Pmax = Quantity(float(tokens2[1].strip()),"atm")
                            if 'TCHEB' in line or 'PCHEB' in line:
                                pass
                            elif kinetics.degreeT == 0 or kinetics.degreeP == 0:
                                tokens2 = tokens[1].split()
                                kinetics.degreeT = int(float(tokens2[0].strip()))
                                kinetics.degreeP = int(float(tokens2[1].strip()))
                                kinetics.coeffs = numpy.zeros((kinetics.degreeT,kinetics.degreeP), numpy.float64)
                            else:
                                tokens2 = tokens[1].split()
                                coeffs = [float(t.strip()) for t in tokens2]
                                for index, C in enumerate(coeffs):
                                    i, j = divmod(index + len(chebyshevCoeffs), kinetics.degreeP)
                                    if i == 0 and j == 0:
                                        C -= math.log10(factorSI) * (len(reactants) - 1)
                                    kinetics.coeffs.value_si[i,j] = C
                                chebyshevCoeffs.extend(coeffs)

                        elif 'LOW' in line:
                            # This line contains low-pressure-limit Arrhenius parameters in Chemkin format

                            # Upgrade the kinetics to a Lindemann if not already done
                            if isinstance(kinetics, Lindemann):
                                pass
                            elif isinstance(kinetics, ThirdBody):
                                kinetics = Lindemann(arrheniusHigh=kinetics.arrheniusLow,
                                                     efficiencies=kinetics.efficiencies,
                                                     comment=kinetics.comment)
                                reaction.kinetics = kinetics
                            elif isinstance(kinetics, Arrhenius):
                                kinetics = Lindemann(arrheniusHigh=kinetics, comment=kinetics.comment)
                                kinetics.arrheniusHigh.comment = ''
                                reaction.kinetics = kinetics

                            items = line.split('/')
                            A, n, Ea = items[1].split()
                            kinetics.arrheniusLow = Arrhenius(
                                A = (float(A), Aunits[len(reactants)+1]), 
                                n = float(n), 
                                Ea = (float(Ea), Eunits), 
                                T0 = (1.0,"K"),
                            )

                        elif 'TROE' in line:
                            # This line contains Troe falloff parameters in Chemkin format

                            # Upgrade the kinetics to a Troe if not already done
                            if isinstance(kinetics, Lindemann):
                                kinetics = Troe(arrheniusLow=kinetics.arrheniusLow,
                                                arrheniusHigh=kinetics.arrheniusHigh,
                                                efficiencies=kinetics.efficiencies,
                                                comment=kinetics.comment)
                                reaction.kinetics = kinetics
                            elif isinstance(kinetics, ThirdBody):
                                kinetics = Troe(arrheniusHigh=kinetics.arrheniusLow,
                                                efficiencies=kinetics.efficiencies,
                                                comment=kinetics.comment)
                                reaction.kinetics = kinetics
                            elif isinstance(kinetics, Arrhenius):
                                kinetics = Troe(arrheniusHigh=kinetics, comment=kinetics.comment)
                                kinetics.arrheniusHigh.comment = ''
                                reaction.kinetics = kinetics

                            items = line.split('/')
                            items = items[1].split()
                            if len(items) == 3:
                                alpha, T3, T1 = items; T2 = None
                            else:
                                alpha, T3, T1, T2 = items

                            kinetics.alpha = float(alpha)
                            kinetics.T1 = (float(T1),"K")
                            if T2 is not None:
                                kinetics.T2 = (float(T2),"K")
                            else:
                                kinetics.T2 = None
                            kinetics.T3 = (float(T3),"K")

                        elif 'DUPLICATE' in line or 'DUP' in line:
                            reaction.duplicate = True
                        
                        else:
                            # This line contains collider efficiencies

                            # Upgrade the kinetics to a Lindemann if not already done
                            if isinstance(kinetics, Arrhenius):
                                kinetics = Lindemann(arrheniusHigh=kinetics, comment=kinetics.comment)
                                kinetics.arrheniusHigh.comment = ''
                                reaction.kinetics = kinetics

                            items = line.split('/')
                            for spec, eff in zip(items[0::2], items[1::2]):
                                spec = str(spec).strip()

                                # In old database, N2, He, Ne, and Ar were treated as special "bath gas" species
                                # These bath gas species were not required to be in the species dictionary
                                # The new database removes this special case, and requires all colliders to be explicitly defined
                                # This is hardcoding to handle these special colliders
                                if spec.upper() in ['N2', 'HE', 'AR', 'NE'] and spec not in species:
                                    if spec.upper() == 'N2':
                                        species[spec] = Species(label='N2', molecule=[Molecule().fromSMILES('N#N')])
                                    elif spec.upper() == 'HE':
                                        species[spec] = Species(label='He', molecule=[Molecule().fromAdjacencyList('1 He 0')])
                                    elif spec.upper() == 'AR':
                                        species[spec] = Species(label='Ar', molecule=[Molecule().fromAdjacencyList('1 Ar 0')])
                                    elif spec.upper() == 'NE':
                                        species[spec] = Species(label='Ne', molecule=[Molecule().fromAdjacencyList('1 Ne 0')])
                                
                                if spec not in species:
                                    logging.warning('Collider {0} for reaction {1} not found in species dictionary.'.format(spec, reaction))
                                else:
                                    kinetics.efficiencies[species[spec].molecule[0]] = float(eff)

                    if 'Unit:' in line:
                        inUnitSection = True; inReactionSection = False
                    elif 'Reactions:' in line:
                        inUnitSection = False; inReactionSection = True
                        
        except (DatabaseError, InvalidAdjacencyListError), e:
            logging.exception('Error while reading old reactions file {0}.'.format(path))
            logging.exception(str(e))
            raise
        except IOError, e:
            logging.exception('Database dictionary file "' + e.filename + '" not found.')
            raise
        finally:
            if fdict: fdict.close()

        return reactions
    
    def saveOld(self, path):
        """
        Save an old-style reaction library to `path`. This creates files named
        ``species.txt``, ``reactions.txt``, and ``pdepreactions.txt`` in the
        given directory; these contain the species dictionary, high-pressure
        limit reactions and kinetics, and pressure-dependent reactions and
        kinetics, respectively.
        """
        try:
            os.makedirs(path)
        except OSError:
            pass
        
        def writeArrhenius(f, arrhenius):
            f.write(' {0:<12.3E} {1:>7.3f} {2:>11.2f}    {3}{4:g} {5:g} {6:g}\n'.format(
                arrhenius.A.value_si,
                arrhenius.n.value_si,
                arrhenius.Ea.value_si / 4.184,
                '*' if arrhenius.A.isUncertaintyMultiplicative() else '',
                arrhenius.A.uncertainty,
                arrhenius.n.uncertainty,
                arrhenius.Ea.uncertainty / 4.184,
            ))
        
        # Gather all of the species used in this kinetics library
        speciesDict = self.getSpecies()
        # Also include colliders in the above
        for entry in self.entries.values():
            if isinstance(entry.data, ThirdBody):
                for molecule in entry.data.efficiencies:
                    formula = molecule.getFormula()
                    if formula in ['He', 'Ar', 'N2', 'Ne']:
                        pass
                    else:
                        found = False
                        for species in speciesDict.values():
                            for mol in species.molecule:
                                if mol.isIsomorphic(molecule):
                                    found = True
                                    break
                        if not found:
                            speciesDict[formula] = Species(label=formula, molecule=[molecule])
        
        entries = self.entries.values()
        entries.sort(key=lambda x: x.index)
        
        # Save the species dictionary
        speciesList = speciesDict.values()
        speciesList.sort(key=lambda x: x.label)
        f = open(os.path.join(path, 'species.txt'), 'w')
        for species in speciesList:
            f.write(species.molecule[0].toAdjacencyList(label=species.label, removeH=True) + "\n")
        f.close()
        
        # Save the high-pressure limit reactions
        # Currently only Arrhenius kinetics are allowed
        f = open(os.path.join(path, 'reactions.txt'), 'w')
        f.write('Unit:\n')
        f.write('A: mol/m3/s\n')
        f.write('E: cal/mol\n\n')
        f.write('Reactions:\n')
        for entry in entries:
            kinetics = entry.data
            rateList = []
            if isinstance(kinetics, MultiArrhenius):
                entry.item.duplicate = True
                rateList = kinetics.arrhenius[:]
            else:
                if not kinetics.isPressureDependent():
                    rateList.append(kinetics)
            for rate in rateList:
                # Write reaction equation
                f.write('{0:<59}'.format(entry.item))
                # Write kinetics
                if isinstance(rate, Arrhenius):
                    writeArrhenius(f, rate)
                else:
                    raise DatabaseError('Unexpected kinetics type "{0}" encountered while saving old kinetics library (reactions.txt).'.format(rate.__class__))
                # Mark as duplicate if needed
                if entry.item.duplicate:
                    f.write(' DUPLICATE\n')
        f.close()
        
        # Save the pressure-dependent reactions
        # Currently only ThirdBody, Lindemann, Troe, and PDepArrhenius kinetics are allowed
        f = open(os.path.join(path, 'pdepreactions.txt'), 'w')
        f.write('Unit:\n')
        f.write('A: mol/m3/s\n')
        f.write('E: cal/mol\n\n')
        f.write('Reactions:\n')
        for entry in entries:
            kinetics = entry.data
            if not kinetics.isPressureDependent():
                continue
            rateList = []
            if isinstance(kinetics, MultiPDepArrhenius):
                entry.item.duplicate = True
                rateList = kinetics.arrhenius[:]
            else:
                rateList.append(kinetics)
            for rate in rateList:
                # Write reaction equation
                equation = str(entry.item)
                if entry.item.reversible:
                    index = equation.find('<=>')
                else:
                    index = equation.find('=>')
                if isinstance(rate, ThirdBody) and not isinstance(rate, Lindemann):
                    equation = '{0}+ M {1} + M'.format(equation[0:index], equation[index:])
                elif isinstance(rate, PDepArrhenius):
                    pass
                else:
                    equation = '{0}(+M) {1} (+M)'.format(equation[0:index], equation[index:]) 
                f.write('{0:<59}'.format(equation))
                # Write kinetics
                if isinstance(rate, (ThirdBody, Lindemann, Troe)):
                    if isinstance(rate, Lindemann):
                        # Lindemann (and Troe) fall-off have the High-P as default, and Low-P labeled LOW
                        writeArrhenius(f, rate.arrheniusHigh)
                    else:
                        # Non-falloff ThirdBody reactions are always in the Low-P limit
                        writeArrhenius(f, rate.arrheniusLow)
                    if len(rate.efficiencies) > 0:
                        eff_line = ''
                        for molecule, efficiency in rate.efficiencies.iteritems():
                            for spec in speciesDict.values():
                                if molecule in spec.molecule:
                                    mol_label = spec.label
                                    break
                            else:
                                mol_label = molecule.getFormula().upper()
                            eff_line += '{0}/{1:g}/  '.format(mol_label, efficiency)
                        f.write(eff_line.strip() + '\n')
                    if isinstance(rate, Lindemann):
                        f.write('     LOW  /  {0:10.3e} {1:9.3f} {2:10.2f}/\n'.format(
                            rate.arrheniusLow.A.value_si,
                            rate.arrheniusLow.n.value_si,
                            rate.arrheniusLow.Ea.value_si / 4.184,
                        ))
                    if isinstance(rate, Troe):
                        if rate.T2 is not None:
                            f.write('     TROE /  {0:10.4f} {1:10.2g} {2:10.2g} {3:10.2g}/\n'.format(
                                rate.alpha,
                                rate.T3.value_si,
                                rate.T1.value_si,
                                rate.T2.value_si,
                            ))
                        else:
                            f.write('     TROE /  {0:10.4f} {1:10.2g} {2:10.2g}/\n'.format(
                                rate.alpha,
                                rate.T3.value_si,
                                rate.T1.value_si,
                            ))
                        
                elif isinstance(rate, PDepArrhenius):
                    writeArrhenius(f, rate.arrhenius[-1])
                    for pressure, arrhenius in zip(rate.pressures.value_si, rate.arrhenius):
                        f.write('     PLOG /  {0:10g} {1:10.3e} {2:9.3f} {3:10.2f} /\n'.format(
                            pressure / 1e5,
                            arrhenius.A.value_si,
                            arrhenius.n.value_si,
                            arrhenius.Ea.value_si / 4.184,
                        ))
                else:
                    raise DatabaseError('Unexpected kinetics type "{0}" encountered while saving old kinetics library (reactions.txt).'.format(rate.__class__))
                # Mark as duplicate if needed
                if entry.item.duplicate:
                    f.write(' DUPLICATE\n')
                f.write('\n')
        f.close()
    