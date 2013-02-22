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

import os
import os.path
import logging
from copy import copy, deepcopy
import numpy

import rmgpy.constants as constants
from rmgpy.kinetics import Arrhenius, ArrheniusEP, ThirdBody, Lindemann, Troe, \
                           PDepArrhenius, MultiArrhenius, MultiPDepArrhenius, \
                           Chebyshev, KineticsData, PDepKineticsModel
from rmgpy.molecule import Molecule, Group
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.data.base import LogicNode

from .common import KineticsError, saveEntry
from .depository import DepositoryReaction, KineticsDepository
from .family import TemplateReaction, KineticsFamily, KineticsGroups, \
    ReactionRecipe, InvalidActionError, ReactionPairsError, \
    UndeterminableKineticsError
from .library import LibraryReaction, KineticsLibrary

################################################################################

class KineticsDatabase(object):
    """
    A class for working with the RMG kinetics database.
    """

    def __init__(self):
        self.families = {}
        self.libraries = {}
        self.libraryOrder = []
        self.local_context = {
            'KineticsData': KineticsData,
            'Arrhenius': Arrhenius,
            'ArrheniusEP': ArrheniusEP,
            'MultiArrhenius': MultiArrhenius,
            'MultiPDepArrhenius': MultiPDepArrhenius,
            'PDepArrhenius': PDepArrhenius,
            'Chebyshev': Chebyshev,
            'ThirdBody': ThirdBody,
            'Lindemann': Lindemann,
            'Troe': Troe,
            'R': constants.R,
        }
        self.global_context = {}

    def __reduce__(self):
        """
        A helper function used when pickling a KineticsDatabase object.
        """
        d = {
            'families': self.families,
            'libraries': self.libraries,
            'libraryOrder': self.libraryOrder,
        }
        return (KineticsDatabase, (), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling a KineticsDatabase object.
        """
        self.families = d['families']
        self.libraries = d['libraries']
        self.libraryOrder = d['libraryOrder']

    def load(self, path, families=None, libraries=None, depositories=None):
        """
        Load the kinetics database from the given `path` on disk, where `path`
        points to the top-level folder of the families database.
        """
        self.loadFamilies(os.path.join(path, 'families'), families, depositories)
        self.loadLibraries(os.path.join(path, 'libraries'), libraries)
        
    def loadFamilies(self, path, families=None, depositories=None):
        """
        Load the kinetics families from the given `path` on disk, where `path`
        points to the top-level folder of the kinetics families.
        """
        logging.info('Loading kinetics families from {0}'.format(path))
        
        familiesToLoad = []
        for (root, dirs, files) in os.walk(os.path.join(path)):
            if root == path:
                if families is None or families == 'all':
                    # All families are loaded by default
                    for d in dirs:
                        familiesToLoad.append(d)
                elif isinstance(families, list) or isinstance(families, tuple):
                    # If all items in the list start with !, all families will be loaded except these
                    if all([label.startswith('!') for label in families]):
                        for d in dirs:
                            if '!{0}'.format(d) not in families:
                                familiesToLoad.append(d)
                    # Otherwise only the families given will be loaded
                    else:
                        for d in dirs:
                            if d in families:
                                familiesToLoad.append(d)
        
        # Now we know what families to load, so let's load them
        self.families = {}
        for label in familiesToLoad:
            familyPath = os.path.join(path, label)
            family = KineticsFamily(label=label)
            family.load(familyPath, self.local_context, self.global_context, depositoryLabels=depositories)
            self.families[label] = family

    def loadLibraries(self, path, libraries=None):
        """
        Load the listed kinetics libraries from the given `path` on disk.
        
        Loads them all if `libraries` list is not specified or `None`.
        The `path` points to the folder of kinetics libraries in the database,
        and the libraries should be in files like :file:`<path>/<library>.py`.
        """
        self.libraries = {}; self.libraryOrder = []
        
        if libraries is not None:
            for library_name in libraries:
                library_file = os.path.join(path, library_name+'.py')
                if os.path.exists(library_file):
                    logging.info('Loading kinetics library {0} from {1}...'.format(library_name, library_file))
                    library = KineticsLibrary(label=library_name)
                    library.load(library_file, self.local_context, self.global_context)
                    self.libraries[library.label] = library
                    self.libraryOrder.append(library.label)
                else:
                    raise IOError("Couldn't find kinetics library {0}".format(library_file))
            assert (self.libraryOrder == libraries)
        else:# load all the libraries you can find
            for (root, dirs, files) in os.walk(os.path.join(path)):
                for f in files:
                    name, ext = os.path.splitext(f)
                    if ext.lower() == '.py':
                        library_file = os.path.join(root, f)
                        label=library_file[len(path)+1:-3]
                        logging.info('Loading kinetics library {0} from {1}...'.format(label, library_file))
                        library = KineticsLibrary(label=label)
                        library.load(library_file, self.local_context, self.global_context)
                        self.libraries[library.label] = library
                        self.libraryOrder.append(library.label)

    def save(self, path):
        """
        Save the kinetics database to the given `path` on disk, where `path`
        points to the top-level folder of the kinetics database.
        """
        path = os.path.abspath(path)
        if not os.path.exists(path): os.mkdir(path)
        self.saveFamilies(os.path.join(path, 'families'))
        self.saveLibraries(os.path.join(path, 'libraries'))

    def saveFamilies(self, path):
        """
        Save the kinetics families to the given `path` on disk, where `path`
        points to the top-level folder of the kinetics families.
        """
        if not os.path.exists(path): os.mkdir(path)
        for label, family in self.families.iteritems():
            familyPath = os.path.join(path, label)
            if not os.path.exists(familyPath): os.mkdir(familyPath)
            family.save(familyPath)

    def saveLibraries(self, path):
        """
        Save the kinetics libraries to the given `path` on disk, where `path`
        points to the top-level folder of the kinetics libraries.
        """
        for label, library in self.libraries.iteritems():
            folders = label.split(os.sep)
            try:
                os.makedirs(os.path.join(path, *folders[:-1]))
            except OSError:
                pass
            library.save(os.path.join(path, '{0}.py'.format(label)))

    def loadOld(self, path):
        """
        Load the old RMG kinetics database from the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        self.families = {}
        self.libraries = {}
        
        librariesPath = os.path.join(path, 'kinetics_libraries')
        for (root, dirs, files) in os.walk(os.path.join(path, 'kinetics_libraries')):
            if os.path.exists(os.path.join(root, 'species.txt')) and os.path.exists(os.path.join(root, 'reactions.txt')):
                library = KineticsLibrary(label=root[len(librariesPath)+1:], name=root[len(librariesPath)+1:])
                library.loadOld(root)
                self.libraries[library.label] = library
                
        for (root, dirs, files) in os.walk(os.path.join(path, 'kinetics_groups')):
            if os.path.exists(os.path.join(root, 'dictionary.txt')) and os.path.exists(os.path.join(root, 'rateLibrary.txt')):
                label = os.path.split(root)[1]
                family = KineticsFamily(label=label)
                family.loadOld(root)
                self.families[family.label] = family

        return self

    def saveOld(self, path):
        """
        Save the old RMG kinetics database to the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        librariesPath = os.path.join(path, 'kinetics_libraries')
        if not os.path.exists(librariesPath): os.mkdir(librariesPath)
        for library in self.libraries.values():
            libraryPath = os.path.join(librariesPath, library.label)
            library.saveOld(libraryPath)

        groupsPath = os.path.join(path, 'kinetics_groups')
        if not os.path.exists(groupsPath): os.mkdir(groupsPath)
        for label, family in self.families.iteritems():
            groupPath = os.path.join(groupsPath, label)
            family.saveOld(groupPath)
    
    def generateReactions(self, reactants, products=None, **options):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be :class:`Molecule` objects. This method
        searches the depository, libraries, and groups, in that order.
        """
        reactionList = []
        reactionList.extend(self.generateReactionsFromLibraries(reactants, products, **options))
        reactionList.extend(self.generateReactionsFromFamilies(reactants, products, **options))
        return reactionList

    def generateReactionsFromLibraries(self, reactants, products, **options):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be :class:`Molecule` objects. This method
        searches the depository.
        """
        reactionList = []
        for label in self.libraryOrder:
            reactionList.extend(self.generateReactionsFromLibrary(reactants, products, self.libraries[label], **options))
        return reactionList

    def generateReactionsFromLibrary(self, reactants, products, library, **options):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be :class:`Molecule` objects. This method
        searches the depository.
        """
        reactionList = []
        for entry in library.entries.values():
            if entry.item.matchesMolecules(reactants):
                reaction = LibraryReaction(
                    reactants = entry.item.reactants[:],
                    products = entry.item.products[:],
                    degeneracy = entry.item.degeneracy,
                    reversible = entry.item.reversible,
                    kinetics = deepcopy(entry.data),
                    library = library,
                    entry = entry,
                )
                reactionList.append(reaction)
        if products:
            reactionList = filterReactions(reactants, products, reactionList)
        return reactionList

    def generateReactionsFromFamilies(self, reactants, products, only_families=None, **options):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be :class:`Molecule` objects. This method
        applies the reaction family.
        If `only_families` is a list of strings, only families with those labels
        are used.
        """
        # If there are two structures and they are the same, then make a copy
        # of the second one so we can independently manipulate both of them 
        # This is for the case where A + A --> products
        if len(reactants) == 2 and reactants[0] == reactants[1]:
            reactants[1] = reactants[1].copy(deep=True)
        
        reactionList = []
        for label, family in self.families.iteritems():
            if only_families is None or label in only_families:
                reactionList.extend(family.generateReactions(reactants, **options))
        if products:
            reactionList = filterReactions(reactants, products, reactionList)
        return reactionList

    def getForwardReactionForFamilyEntry(self, entry, family, thermoDatabase):
        """
        For a given `entry` for a reaction of the given reaction `family` (the
        string label of the family), return the reaction with kinetics and
        degeneracy for the "forward" direction as defined by the reaction 
        family. For families that are their own reverse, the direction the
        kinetics is given in will be preserved. If the entry contains 
        functional groups for the reactants, assume that it is given in the 
        forward direction and do nothing. Returns the reaction in the direction
        consistent with the reaction family template, and the matching template.
        Note that the returned reaction will have its kinetics and degeneracy
        set appropriately.
        
        In order to reverse the reactions that are given in the reverse of the
        direction the family is defined, we need to compute the thermodynamics
        of the reactants and products. For this reason you must also pass
        the `thermoDatabase` to use to generate the thermo data.
        """

        def generateThermoData(species, thermoDatabase):
            thermoData = [thermoDatabase.getThermoData(species)]
            thermoData.sort(key=lambda x: x.getEnthalpy(298))
            return thermoData[0]
        
        def matchSpeciesToMolecules(species, molecules):
            if len(species) == len(molecules) == 1:
                return species[0].isIsomorphic(molecules[0])
            elif len(species) == len(molecules) == 2:
                if species[0].isIsomorphic(molecules[0]) and species[1].isIsomorphic(molecules[1]):
                    return True
                elif species[0].isIsomorphic(molecules[1]) and species[1].isIsomorphic(molecules[0]):
                    return True
            return False

        reaction = None; template = None

        # Get the indicated reaction family
        try:
            groups = self.families[family].groups
        except KeyError:
            raise ValueError('Invalid value "{0}" for family parameter.'.format(family))

        if all([(isinstance(reactant, Group) or isinstance(reactant, LogicNode)) for reactant in entry.item.reactants]):
            # The entry is a rate rule, containing functional groups only
            # By convention, these are always given in the forward direction and
            # have kinetics defined on a per-site basis
            reaction = Reaction(
                reactants = entry.item.reactants[:],
                products = [],
                kinetics = entry.data,
                degeneracy = 1,
            )
            template = [groups.entries[label] for label in entry.label.split(';')]

        elif (all([isinstance(reactant, Molecule) for reactant in entry.item.reactants]) and
            all([isinstance(product, Molecule) for product in entry.item.products])):
            # The entry is a real reaction, containing molecules
            # These could be defined for either the forward or reverse direction
            # and could have a reaction-path degeneracy

            reaction = Reaction(reactants=[], products=[])
            for molecule in entry.item.reactants:
                reactant = Species(molecule=[molecule])
                reactant.generateResonanceIsomers()
                reactant.thermo = generateThermoData(reactant, thermoDatabase)
                reaction.reactants.append(reactant)
            for molecule in entry.item.products:
                product = Species(molecule=[molecule])
                product.generateResonanceIsomers()
                product.thermo = generateThermoData(product, thermoDatabase)
                reaction.products.append(product)

            # Generate all possible reactions involving the reactant species
            generatedReactions = self.generateReactionsFromFamilies([reactant.molecule for reactant in reaction.reactants], [], only_families=[family])

            # Remove from that set any reactions that don't produce the desired reactants and products
            forward = []; reverse = []
            for rxn in generatedReactions:
                if matchSpeciesToMolecules(reaction.reactants, rxn.reactants) and matchSpeciesToMolecules(reaction.products, rxn.products):
                    forward.append(rxn)
                if matchSpeciesToMolecules(reaction.reactants, rxn.products) and matchSpeciesToMolecules(reaction.products, rxn.reactants):
                    reverse.append(rxn)

            # We should now know whether the reaction is given in the forward or
            # reverse direction
            if len(forward) == 1 and len(reverse) == 0:
                # The reaction is in the forward direction, so use as-is
                reaction = forward[0]
                template = reaction.template
                # Don't forget to overwrite the estimated kinetics from the database with the kinetics for this entry
                reaction.kinetics = entry.data
            elif len(reverse) == 1 and len(forward) == 0:
                # The reaction is in the reverse direction
                # First fit Arrhenius kinetics in that direction
                Tdata = 1000.0 / numpy.arange(0.5, 3.301, 0.1, numpy.float64)
                kdata = numpy.zeros_like(Tdata)
                for i in range(Tdata.shape[0]):
                    kdata[i] = entry.data.getRateCoefficient(Tdata[i]) / reaction.getEquilibriumConstant(Tdata[i])
                kunits = 'm^3/(mol*s)' if len(reverse[0].reactants) == 2 else 's^-1'
                kinetics = Arrhenius().fitToData(Tdata, kdata, kunits, T0=1.0)
                kinetics.Tmin = entry.data.Tmin
                kinetics.Tmax = entry.data.Tmax
                kinetics.Pmin = entry.data.Pmin
                kinetics.Pmax = entry.data.Pmax
                # Now flip the direction
                reaction = reverse[0]
                reaction.kinetics = kinetics
                template = reaction.template
            elif len(reverse) > 0 and len(forward) > 0:
                print 'FAIL: Multiple reactions found for {0!r}.'.format(entry.label)
            elif len(reverse) == 0 and len(forward) == 0:
                print 'FAIL: No reactions found for "%s".' % (entry.label)
            else:
                print 'FAIL: Unable to estimate kinetics for {0!r}.'.format(entry.label)

        assert reaction is not None
        assert template is not None
        return reaction, template

################################################################################

def filterReactions(reactants, products, reactionList):
    """
    Remove any reactions from the given `reactionList` whose reactants do
    not involve all the given `reactants` or whose products do not involve 
    all the given `products`. This method checks both forward and reverse
    directions, and only filters out reactions that don't match either.
    """
    
    # Convert from molecules to species and generate resonance isomers.
    reactant_species = []
    for mol in reactants:
        s = Species(molecule=[mol])
        s.generateResonanceIsomers()
        reactant_species.append(s)
    reactants = reactant_species
    product_species = []
    for mol in products:
        s = Species(molecule=[mol])
        s.generateResonanceIsomers()
        product_species.append(s)
    products = product_species
    
    reactions = reactionList[:]
    
    for reaction in reactionList:
        # Forward direction
        reactants0 = [r for r in reaction.reactants]
        for reactant in reactants:
            for reactant0 in reactants0:
                if reactant.isIsomorphic(reactant0):
                    reactants0.remove(reactant0)
                    break
        products0 = [p for p in reaction.products]
        for product in products:
            for product0 in products0:
                if product.isIsomorphic(product0):
                    products0.remove(product0)
                    break
        forward = not (len(reactants0) != 0 or len(products0) != 0)
        # Reverse direction
        reactants0 = [r for r in reaction.products]
        for reactant in reactants:
            for reactant0 in reactants0:
                if reactant.isIsomorphic(reactant0):
                    reactants0.remove(reactant0)
                    break
        products0 = [p for p in reaction.reactants]
        for product in products:
            for product0 in products0:
                if product.isIsomorphic(product0):
                    products0.remove(product0)
                    break
        reverse = not (len(reactants0) != 0 or len(products0) != 0)
        if not forward and not reverse:
            reactions.remove(reaction)
    return reactions
