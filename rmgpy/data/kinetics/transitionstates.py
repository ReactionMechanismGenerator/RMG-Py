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
This module contains functionality for working with transition state family functional
groups, including support for using group additivity to estimate TS geometries.
"""

import os
import os.path
import re
import logging
import codecs
import math
import numpy
from copy import copy, deepcopy

from rmgpy.data.base import *

from rmgpy.quantity import Quantity, constants
from rmgpy.reaction import Reaction, ReactionError
from rmgpy.molecule import Bond, GroupBond, Group
from rmgpy.species import Species

from .common import KineticsError, UndeterminableKineticsError, saveEntry

################################################################################

class DistanceData():
    """
    A class for storing distance matrix data, for geometry estimation
    """
    def __init__(self, distances={}, uncertainties=None, method=None):
        self.distances = distances
        self.uncertainties = uncertainties
        self.method = method
        self.comment = u''
        assert isinstance(distances,dict), "distances should be a dict"
        if method: assert isinstance(method,str), "method should be a string"
        
    def __repr__(self):
        strings = ["DistanceData("]

        strings.append("distances={")
        for key in sorted(self.distances.keys()):
            strings.append("{0!r}: {1:.6f},".format(key, self.distances[key]))
        strings.append("}")

        if self.uncertainties is not None:
            strings.append(", uncertainties={")
            for key in sorted(self.uncertainties.keys()):
                strings.append("{0!r}: {1:.6f},".format(key, self.uncertainties[key]))
            strings.append("}")

        if self.method:
            strings.append(", method={0!r}".format(self.method))
        if self.comment:
            strings.append(", comment={0!r}".format(self.comment))
        strings.append(")")
        return ''.join(strings)
    
    def add(self, other):
        """Adds the `other` distances to these."""
        assert len(self.distances)==len(other.distances), "self and other must have the same size dictionary of distances, but self={0!r} and other={1!r}".format(self,other)
        for key, value in other.distances.iteritems():
            self.distances[key] += value
        if self.uncertainties and other.uncertainties:
            for key, value in other.uncertainties.iteritems():
                self.uncertainties[key] += value
        else:
            self.uncertainties = None
        
    def __copy__(self):
        return DistanceData(distances=self.distances.copy(), uncertainties=self.uncertainties.copy(), method=self.method)
        

class TransitionStates(Database):
    """
    loads, and contains, both Depository and Groups
    """
    def __init__(self):
        self.groups = None
        self.depository = None
        self.family = None

    def load(self, path, local_context, global_context):
        """
        Load the TS database
        """
        if local_context is None: local_context = {}
        local_context['DistanceData'] = DistanceData
        
        fpath = os.path.join(path,'TS_training', 'reactions.py')
        logging.debug("Loading transitions state family training set from {0}".format(fpath))
        depository = TransitionStateDepository(label='{0}/TS_training'.format(path.split('/')[-1]))#'intra_H_migration/TS_training')
        depository.load(fpath, local_context, global_context )
        self.depository = depository
        
        fpath = os.path.join(path,'TS_groups.py')
        logging.debug("Loading transitions state family groups from {0}".format(fpath))
        groups = TSGroups(label='{0}/TS_groups'.format(path.split('/')[-1]))#'intra_H_migration/TS_groups')
        groups.load(fpath , local_context, global_context )
        
        self.family.forwardTemplate.reactants = [groups.entries[entry.label] for entry in self.family.forwardTemplate.reactants]
        # self.family.forwardTemplate.products = [groups.entries[entry.label] for entry in self.family.forwardTemplate.products]
        self.family.entries = groups.entries
        self.family.groups = groups
        groups.numReactants = len(self.family.forwardTemplate.reactants)
        self.groups = groups
    
    def estimateDistances(self, reaction):
        """
        Return estimated DistanceData for the given reaction
        """
        # Should check depository first, but for now just go straight to group additive estimate:
        return self.groups.estimateDistancesUsingGroupAdditivity(reaction)
    
    def saveTransitionStateGroups(self, path, entryName='entry'):
        """
        Save the current database to the file at location `path` on disk. The
        optional `entryName` parameter specifies the identifier used for each
        data entry.
        """
        entries = self.groups.getEntriesToSave()
                
        # Write the header
        f = codecs.open(path, 'w', 'utf-8')
        f.write('#!/usr/bin/env python\n')
        f.write('# encoding: utf-8\n\n')
        f.write('name = "{0}"\n'.format(self.groups.name))
        f.write('shortDesc = u"{0}"\n'.format(self.groups.shortDesc))
        f.write('longDesc = u"""\n')
        f.write(self.groups.longDesc)
        f.write('\n"""\n\n')

        # Save the entries
        for entry in entries:
            saveEntry(f, entry)

        # Write the tree
        if len(self.groups.top) > 0:
            f.write('tree(\n')
            f.write('"""\n')
            f.write(self.generateOldTree(self.groups.top, 1))
            f.write('"""\n')
            f.write(')\n\n')

        f.close()
    
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
        
    def generateReactionsFromFamilies(self, reactants, products, only_families=None, families=None, **options):
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
        reactionList.extend(families.generateReactions(reactants, **options))
        
        if products:
            reactionList = filterReactions(reactants, products, reactionList)
        
        return reactionList
        
    def getForwardReactionForFamilyEntry(self, entry, family, groups, rxnFamily):
        """
        For a given `entry` for a reaction of the given reaction `family` (the
        string label of the family), return the reaction with transition state
        distances for the "forward" direction as defined by the reaction 
        family. For families that are their own reverse, the direction the
        kinetics is given in will be preserved. If the entry contains 
        functional groups for the reactants, assume that it is given in the 
        forward direction and do nothing. Returns the reaction in the direction
        consistent with the reaction family template, and the matching template.
        """
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
        if groups == None:
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
    
        elif (all([isinstance(reactant, (Molecule, Species)) for reactant in entry.item.reactants]) and
            all([isinstance(product, (Molecule, Species)) for product in entry.item.products])):
            # The entry is a real reaction, containing molecules
            # These could be defined for either the forward or reverse direction
            # and could have a reaction-path degeneracy
    
            reaction = Reaction(reactants=[], products=[])
            for molecule in entry.item.reactants:
                if isinstance(molecule, Molecule):
                    reactant = Species(molecule=[molecule])
                else:
                    reactant = molecule
                reactant.generateResonanceIsomers()
                reaction.reactants.append(reactant)
            for molecule in entry.item.products:
                if isinstance(molecule, Molecule):
                    product = Species(molecule=[molecule])
                else:
                    product = molecule
                product.generateResonanceIsomers()
                reaction.products.append(product)
            
            # Generate all possible reactions involving the reactant species
            generatedReactions = self.generateReactionsFromFamilies([reactant.molecule[0] for reactant in reaction.reactants], [], only_families=[family], families=rxnFamily)
            
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
                # Don't forget to overwrite the estimated distances from the database with the distances for this entry
                reaction.distances = entry.data
            elif len(reverse) == 1 and len(forward) == 0:
                # The reaction is in the reverse direction
                # The reaction is in the forward direction, so use as-is
                reaction = forward[0]
                template = reaction.template
                # The distances to the H atom are reversed
                reaction.distances = entry.data
                reaction.distances['d12'] = entry.data['d23']
                reaction.distances['d23'] = entry.data['d12']
            elif len(reverse) > 0 and len(forward) > 0:
                print 'FAIL: Multiple reactions found for {0!r}.'.format(entry.label)
            elif len(reverse) == 0 and len(forward) == 0:
                print 'FAIL: No reactions found for "%s".' % (entry.label)
            else:
                print 'FAIL: Unable to estimate distances for {0!r}.'.format(entry.label)
                
        assert reaction is not None
        assert template is not None
        return reaction, template

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
        s = Species(molecule=mol)
        s.generateResonanceIsomers()
        reactant_species.append(s)
    reactants = reactant_species
    product_species = []
    for mol in products:
        s = Species(molecule=mol)
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
    
################################################################################

class TransitionStateDepository(Database):
    """
    A class for working with an RMG transition state depository. Each depository 
    corresponds to a reaction family (a :class:`KineticsFamily` object). Each
    entry in a transition state depository involves a reaction defined either by a
    real reactant and product species.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def __repr__(self):
        return '<TransitionStateDepository "{0}">'.format(self.label)
    
    def load(self, path, local_context=None, global_context=None):
        
        Database.load(self, path, local_context, global_context)
        
        # Load the species in the kinetics library
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
            assert reversible == rxn.reversible
            for reactant in reactants.split('+'):
                reactant = reactant.strip()
                if reactant not in speciesDict:
                    raise DatabaseError('Species {0} in kinetics depository {1} is missing from its dictionary.'.format(reactant, self.label))
                # For some reason we need molecule objects in the depository rather than species objects
                rxn.reactants.append(speciesDict[reactant])
            for product in products.split('+'):
                product = product.strip()
                if product not in speciesDict:
                    raise DatabaseError('Species {0} in kinetics depository {1} is missing from its dictionary.'.format(product, self.label))
                # For some reason we need molecule objects in the depository rather than species objects
                rxn.products.append(speciesDict[product])
                
            if not rxn.isBalanced():
                raise DatabaseError('Reaction {0} in kinetics depository {1} was not balanced! Please reformulate.'.format(rxn, self.label))
                
    def loadEntry(self,
                  index,
                  reactant1=None,
                  reactant2=None,
                  reactant3=None,
                  product1=None,
                  product2=None,
                  product3=None,
                  distances=None,
                  degeneracy=1,
                  label='',
                  duplicate=False,
                  reversible=True,
                  reference=None,
                  referenceType='',
                  shortDesc='',
                  longDesc='',
                  rank=None,
                  ):

        reaction = Reaction(reactants=[], products=[], degeneracy=degeneracy, duplicate=duplicate, reversible=reversible)

        entry = Entry(
            index = index,
            label = label,
            item = reaction,
            data = distances,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            rank = rank,
        )
        self.entries['{0:d}:{1}'.format(index,label)] = entry
        return entry

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the database to the file object `f`.
        """
        return saveEntry(f, entry)


################################################################################

class TSGroups(Database):
    """
    A class for working with group additivity values for transition state distances.
    """

    def __init__(self,
                 entries=None,
                 top=None,
                 label='',
                 name='',
                 shortDesc='',
                 longDesc='',
                 forwardTemplate=None,
                 forwardRecipe=None,
                 reverseTemplate=None,
                 reverseRecipe=None,
                 forbidden=None
                 ):
        Database.__init__(self, entries, top, label, name, shortDesc, longDesc)
        self.numReactants = 0
        
    def __repr__(self):
        return '<TSGroups "{0}">'.format(self.label)

    def loadEntry(self, index, label, group, distances, reference=None, referenceType='', shortDesc='', longDesc=''):
        if group[0:3].upper() == 'OR{' or group[0:4].upper() == 'AND{' or group[0:7].upper() == 'NOT OR{' or group[0:8].upper() == 'NOT AND{':
            item = makeLogicNode(group)
        else:
            item = Group().fromAdjacencyList(group)
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = item,
            data = distances,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
        )

    def getReactionTemplate(self, reaction):
        """
        For a given `reaction` with properly-labeled :class:`Molecule` objects
        as the reactants, determine the most specific nodes in the tree that
        describe the reaction.
        """
        # from .family import TemplateReaction
        #assert isinstance(reaction, TemplateReaction), "Can only match TemplateReactions"
        # Get forward reaction template and remove any duplicates
        forwardTemplate = self.top[:]
        temporary = []
        symmetricTree = False
        for entry in forwardTemplate:
            if entry not in temporary:
                temporary.append(entry)
            else:
                # duplicate node found at top of tree
                # eg. R_recombination: ['Y_rad', 'Y_rad']
                assert len(forwardTemplate)==2 , 'Can currently only do symmetric trees with nothing else in them'
                symmetricTree = True
        forwardTemplate = temporary

        # Descend reactant trees as far as possible
        template = []
        for entry in forwardTemplate:
            # entry is a top-level node that should be matched
            group = entry.item

            # To sort out "union" groups, descend to the first child that's not a logical node
            # ...but this child may not match the structure.
            # eg. an R3 ring node will not match an R4 ring structure.
            # (but at least the first such child will contain fewest labels - we hope)
            if isinstance(entry.item, LogicNode):
                group = entry.item.getPossibleStructures(self.entries)[0]

            atomList = group.getLabeledAtoms() # list of atom labels in highest non-union node
            
            for reactant in reaction.reactants:
                if isinstance(reactant, Species):
                    reactant = reactant.molecule[0]
                # Match labeled atoms
                # Check this reactant has each of the atom labels in this group
                if not all([reactant.containsLabeledAtom(label) for label in atomList]):
                    continue # don't try to match this structure - the atoms aren't there!
                # Match structures
                atoms = reactant.getLabeledAtoms()
                
                matched_node = self.descendTree(reactant, atoms, root=entry)
                if matched_node is not None:
                    template.append(matched_node)
                #else:
                #    logging.warning("Couldn't find match for {0} in {1}".format(entry,atomList))
                #    logging.warning(reactant.toAdjacencyList())

        # Get fresh templates (with duplicate nodes back in)
        forwardTemplate = self.top[:]
        if self.label.lower().startswith('r_recombination'):
            forwardTemplate.append(forwardTemplate[0])

        # Check that we were able to match the template.
        # template is a list of the actual matched nodes
        # forwardTemplate is a list of the top level nodes that should be matched
        if len(template) != len(forwardTemplate):
            logging.warning('Unable to find matching template for reaction {0} in reaction family {1}'.format(str(reaction), str(self)) )
            logging.warning(" Trying to match " + str(forwardTemplate))
            logging.warning(" Matched "+str(template))
            print str(self), template, forwardTemplate
            for n,reactant in enumerate(reaction.reactants):
                print "Reactant", n
                print reactant.toAdjacencyList() + '\n'
            for n,product in enumerate(reaction.products):
                print "Product", n
                print product.toAdjacencyList() + '\n'
            raise UndeterminableKineticsError(reaction)
        
        for reactant in reaction.reactants:
            if isinstance(reactant, Species):
                reactant = reactant.molecule[0]
            #reactant.clearLabeledAtoms()

        return template

    def estimateDistancesUsingGroupAdditivity(self, reaction):
        """
        Determine the appropriate transition state distances for a reaction 
        with the given `template` using group additivity.
        """
        template = self.getReactionTemplate(reaction)
        referenceDistances = self.top[0].data # or something like that

        # Start with the generic distances of the top-level nodes
        # Make a copy so we don't modify the original
        tsDistances = deepcopy(referenceDistances)
        
        # Now add in more specific corrections if possible
        for node in template:
            entry = node
            comment_line = "Matched node "
            while not entry.data.distances and entry not in self.top:
                # Keep climbing tree until you find a (non-top) node with distances.
                comment_line += "{0} >> ".format(entry.label)
                entry = entry.parent
            if entry.data.distances and entry not in self.top:
                tsDistances.add(entry.data)
                comment_line += "{0} ({1})".format(entry.label, entry.longDesc.split('\n')[0])
            elif entry in self.top:
                comment_line += "{0} (Top node)".format(entry.label)
            tsDistances.comment += comment_line + '\n'
        
        return tsDistances


    def generateGroupAdditivityValues(self, trainingSet):
        """
        Generate the group additivity values using the given `trainingSet`,
        a list of 2-tuples of the form ``(template, kinetics)``. You must also
        specify the `kunits` for the family and the `method` to use when
        generating the group values. Returns ``True`` if the group values have
        changed significantly since the last time they were fitted, or ``False``
        otherwise.
        """        
        # keep track of previous values so we can detect if they change
        old_entries = dict()
        for label,entry in self.entries.items():
            if entry.data is not None:
                old_entries[label] = entry.data
        
        # Determine a complete list of the entries in the database, sorted as in the tree
        groupEntries = self.top[:]
        
        for entry in self.top:
            groupEntries.extend(self.descendants(entry)) # Entries in the TS_group.py tree
        
        # Determine a unique list of the groups we will be able to fit parameters for
        groupList = []
        for template, distances in trainingSet:
            for group in template:
                if isinstance(group, str): group = self.entries[group]
                if group not in self.top:
                    groupList.append(group)
                    groupList.extend(self.ancestors(group)[:-1])
        groupList = list(set(groupList))
        groupList.sort(key=lambda x: x.index)
        
        if True: # should remove this IF block, as we only have one method.
            # Initialize dictionaries of fitted group values and uncertainties
            groupValues = {}; groupUncertainties = {}; groupCounts = {}; groupComments = {}
            for entry in groupEntries:
                groupValues[entry] = []
                groupUncertainties[entry] = []
                groupCounts[entry] = []
                groupComments[entry] = set()
            
            # Generate least-squares matrix and vector
            A = []; b = []
            
            distance_keys = sorted(trainingSet[0][1].distances.keys())  # ['d12', 'd13', 'd23']
            distance_data = []
            for template, distanceData in trainingSet:
                d = [distanceData.distances[key] for key in distance_keys]
                distance_data.append(d)
                    
                # Create every combination of each group and its ancestors with each other
                combinations = []
                for group in template:
                    groups = [group]; groups.extend(self.ancestors(group)) # Groups from the group.py tree
                    combinations.append(groups)
                combinations = getAllCombinations(combinations)
                # Add a row to the matrix for each combination
                for groups in combinations:
                    Arow = [1 if group in groups else 0 for group in groupList]
                    Arow.append(1)
                    brow = d
                    A.append(Arow); b.append(brow)
                    
                    for group in groups:
                        if isinstance(group, str): group = self.entries[group]
                        groupComments[group].add("{0!s}".format(template))
            
            if len(A) == 0:
                logging.warning('Unable to fit kinetics groups for family "{0}"; no valid data found.'.format(self.label))
                return
            A = numpy.array(A)
            b = numpy.array(b)
            distance_data = numpy.array(distance_data)
            
            x, residues, rank, s = numpy.linalg.lstsq(A, b)
            
            for t, distance_key in enumerate(distance_keys):
                
                # Determine error in each group
                stdev = numpy.zeros(len(groupList)+1, numpy.float64)
                count = numpy.zeros(len(groupList)+1, numpy.int)
                
                for index in range(len(trainingSet)):
                    template, distances = trainingSet[index]
                    d = numpy.float64(distance_data[index,t])
                    dm = x[-1,t] + sum([x[groupList.index(group),t] for group in template if group in groupList])
                    variance = (dm - d)**2
                    for group in template:
                        groups = [group]
                        groups.extend(self.ancestors(group))
                        for g in groups:
                            if g.label not in [top.label for top in self.top]:
                                ind = groupList.index(g)
                                stdev[ind] += variance
                                count[ind] += 1
                    stdev[-1] += variance
                    count[-1] += 1
                
                import scipy.stats
                ci = numpy.zeros(len(count))
                for i in range(len(count)):
                    if count[i] > 1:
                        stdev[i] = numpy.sqrt(stdev[i] / (count[i] - 1))
                        ci[i] = scipy.stats.t.ppf(0.975, count[i] - 1) * stdev[i]
                    else:
                        stdev[i] = None
                        ci[i] = None
                # Update dictionaries of fitted group values and uncertainties
                for entry in groupEntries:
                    if entry == self.top[0]:
                        groupValues[entry].append(x[-1,t])
                        groupUncertainties[entry].append(ci[-1])
                        groupCounts[entry].append(count[-1])
                    elif entry.label in [group.label for group in groupList]:
                        index = [group.label for group in groupList].index(entry.label)
                        groupValues[entry].append(x[index,t])
                        groupUncertainties[entry].append(ci[index])
                        groupCounts[entry].append(count[index])
                    else:
                        groupValues[entry] = None
                        groupUncertainties[entry] = None
                        groupCounts[entry] = None
            
            # Store the fitted group values and uncertainties on the associated entries
            for entry in groupEntries:
                if groupValues[entry] is not None:
                    if not any(numpy.isnan(numpy.array(groupUncertainties[entry]))):
                        # should be entry.data.* (e.g. entry.data.uncertainties)
                        uncertainties = numpy.array(groupUncertainties[entry])
                        uncertaintyType = '+|-'
                    else:
                        uncertainties = {}
                    # should be entry.*
                    shortDesc = "Fitted to {0} distances.\n".format(groupCounts[entry][0])
                    longDesc = "\n".join(groupComments[entry.label])
                    distances_dict = {key:distance for key, distance in zip(distance_keys, groupValues[entry])}
                    uncertainties_dict = {key:distance for key, distance in zip(distance_keys, uncertainties)}
                    entry.data = DistanceData(distances=distances_dict, uncertainties=uncertainties_dict)
                    entry.shortDesc = shortDesc
                    entry.longDesc = longDesc
                else:
                    entry.data = DistanceData()
        
        changed = False
        for label, entry in self.entries.items():
            if entry.data is not None:
                continue # because this is broken:
                if old_entries.has_key(label):
                    old_entry = old_entries[label][label][0]
                    for key, distance in entry.data.iteritems():
                        diff = 0
                        for k in range(3):
                            diff += abs(distance[0][k]/old_entry[k] - 1)
                        if diff > 0.01:
                            changed = True
                            entry.history.append(event)
            else:
                changed = True
                entry.history.append(event)
        return True # because the thing above is broken
        return changed
        # below is what has been updated
        # # Add a note to the history of each changed item indicating that we've generated new group values
        # import time
        # changed = False
        # for label, entry in self.entries.items():
        #     if entry.data is not None and old_entries.has_key(label):
        #         if (isinstance(entry.data, KineticsData) and 
        #             isinstance(old_entries[label], KineticsData) and
        #             len(entry.data.kdata.value_si) == len(old_entries[label].kdata.value_si) and
        #             all(abs(entry.data.kdata.value_si / old_entries[label].kdata.value_si - 1) < 0.01)):
        #             #print "New group values within 1% of old."
        #             pass
        #         elif (isinstance(entry.data, Arrhenius) and 
        #             isinstance(old_entries[label], Arrhenius) and
        #             abs(entry.data.A.value_si / old_entries[label].A.value_si - 1) < 0.01 and
        #             abs(entry.data.n.value_si / old_entries[label].n.value_si - 1) < 0.01 and
        #             abs(entry.data.Ea.value_si / old_entries[label].Ea.value_si - 1) < 0.01 and
        #             abs(entry.data.T0.value_si / old_entries[label].T0.value_si - 1) < 0.01):
        #             #print "New group values within 1% of old."
        #             pass
        #         else:
        #             changed = True
        #             break
        #     else:
        #         changed = True
        #         break
        # 
        # return changed

        
