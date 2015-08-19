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

import os.path
import logging

from rmgpy.data.base import DatabaseError, Database, Entry

from rmgpy.reaction import Reaction
from rmgpy.kinetics import Arrhenius, ThirdBody, Lindemann, Troe, \
                           PDepArrhenius, MultiArrhenius, MultiPDepArrhenius, \
                           PDepKineticsModel
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from .common import saveEntry

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

    def __init__(self, label='', name='', solvent=None, shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)
        
    def __str__(self):
        return 'Kinetics Library {0}'.format(self.label)
    
    def __repr__(self):
        return '<KineticsLibrary "{0}">'.format(self.label)

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
                    
    def checkForDuplicates(self, markDuplicates=False):
        """
        Check that all duplicate reactions in the kinetics library are
        properly marked (i.e. with their ``duplicate`` attribute set to 
        ``True``).
        If ``markDuplicates`` is set to ``True``, then ignore and
        mark all duplicate reactions as duplicate.
        """
        for entry0 in self.entries.values():
            reaction0 = entry0.item
            if not reaction0.duplicate:
                # This reaction is not marked as a duplicate reaction
                # This means that if we find any duplicate reactions, it is an error
                for entry in self.entries.values():
                    reaction = entry.item
                    if reaction0 is not reaction and reaction0.isIsomorphic(reaction): 
                        # We found a duplicate reaction that wasn't marked!
                        # RMG requires all duplicate reactions to be marked, unlike CHEMKIN
                        if markDuplicates:
                            reaction0.duplicate = reaction.duplicate = True
                            logging.warning('Reaction indices {0} and {1} were marked as duplicate.'.format(entry0.index, entry.index))
                            continue
                        
                        raise DatabaseError('Unexpected duplicate reaction {0} in kinetics library {1}. Reaction index {2} matches index {3}.'.format(reaction0, self.label, entry.index, entry0.index))        

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
            if len(duplicates)<=1:
                continue
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
                logging.warning('Only Arrhenius and PDepArrhenius kinetics supported for duplicate reactions.')
                continue
            entry0.longDesc = longDesc
            entries_to_remove.extend(duplicates[1:])
        for entry in entries_to_remove:
            print "removing duplicate reaction with index {0}.".format(entry.index)
            del(self.entries[entry.index])
        print "NB. the entries have not been renumbered, so these indices are missing."
        
        
    def load(self, path, local_context=None, global_context=None):
        Database.load(self, path, local_context, global_context)
        
        # Load a unique set of the species in the kinetics library
        speciesDict = self.getSpecies(os.path.join(os.path.dirname(path),'dictionary.txt'))
        # Make sure all of the reactions draw from only this set
        entries = self.entries.values()
        for entry in entries:
            # Create a new reaction per entry
            rxn = entry.item
            rxn_string = entry.label
            # Convert the reactants and products to Species objects using the speciesDict
            reactants, products = rxn_string.split('=')
            reversible = True
            if '<=>' in rxn_string:
                reactants = reactants[:-1]
                products = products[1:]
            elif '=>' in rxn_string:
                products = products[1:]
                reversible = False
            assert reversible == rxn.reversible, "Reaction string reversibility (=>) and entry attribute `reversible` (set to `False`) must agree if reaction is irreversible."
            for reactant in reactants.split('+'):
                reactant = reactant.strip()
                if reactant not in speciesDict:
                    raise DatabaseError('Species {0} in kinetics library {1} is missing from its dictionary.'.format(reactant, self.label))
                rxn.reactants.append(speciesDict[reactant])
            for product in products.split('+'):
                product = product.strip()
                if product not in speciesDict:
                    raise DatabaseError('Species {0} in kinetics library {1} is missing from its dictionary.'.format(product, self.label))
                rxn.products.append(speciesDict[product])
                
            if not rxn.isBalanced():
                raise DatabaseError('Reaction {0} in kinetics library {1} was not balanced! Please reformulate.'.format(rxn, self.label))    
            
        self.checkForDuplicates()
        
    def loadEntry(self,
                  index,
                  label,
                  kinetics,
                  degeneracy=1,
                  duplicate=False,
                  reversible=True,
                  reference=None,
                  referenceType='',
                  shortDesc='',
                  longDesc='',
                  ):
        
#        reactants = [Species(label=reactant1.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(reactant1)])]
#        if reactant2 is not None: reactants.append(Species(label=reactant2.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(reactant2)]))
#        if reactant3 is not None: reactants.append(Species(label=reactant3.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(reactant3)]))
#
#        products = [Species(label=product1.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(product1)])]
#        if product2 is not None: products.append(Species(label=product2.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(product2)]))
#        if product3 is not None: products.append(Species(label=product3.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(product3)]))
#        
        # Make a blank reaction
        rxn = Reaction(reactants=[], products=[], degeneracy=degeneracy, duplicate=duplicate, reversible=reversible)
#        if not rxn.isBalanced():
#            raise DatabaseError('Reaction {0} in kinetics library {1} was not balanced! Please reformulate.'.format(rxn, self.label))        
#        label = str(rxn)
        assert index not in self.entries, "Reaction {0} already present!".format(label)
        self.entries[index] = Entry(
            index = index,
            label = label,
            item = rxn,
            data = kinetics,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
        )

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
        
        # Add common bath gases (Ar, Ne, He, N2) if not already present
        for label, smiles in [('Ar','[Ar]'), ('He','[He]'), ('Ne','[Ne]'), ('N2','N#N')]:
            if label not in species:
                molecule = Molecule().fromSMILES(smiles)
                spec = Species(label=label, molecule=[molecule])
                species[label] = spec                         

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
                label = str(reaction)
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
        from rmgpy.chemkin import readReactionsBlock
        f = open(path, 'r')
        reactionList = readReactionsBlock(f, speciesDict=species)
        f.close()
        return reactionList
    
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
            f.write(species.molecule[0].toAdjacencyList(label=species.label, removeH=False) + "\n")
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
    
