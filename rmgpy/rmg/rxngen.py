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
Contains functions for generating reactions based on the loaded templates in
the kinetics database.
"""

import logging

from rmgpy.chem.pattern import MoleculePattern

from rmgpy.data.base import LogicNode
import rmgpy.data.kinetics

################################################################################

def matchReactantToTemplate(reactant, templateReactant, family):
    """
    Return ``True`` if the provided reactant matches the provided
    template reactant and ``False`` if not, along with a complete list of the
    mappings.
    """

    if templateReactant.__class__ == list: templateReactant = templateReactant[0]
    struct = family.dictionary[templateReactant]

    if isinstance(struct, LogicNode):
        mappings = []
        for child_structure in struct.getPossibleStructures(family.dictionary):
            ismatch, map = reactant.findSubgraphIsomorphisms(child_structure)
            if ismatch: mappings.extend(map)
        return len(mappings) > 0, mappings
    elif isinstance(struct, MoleculePattern):
        return reactant.findSubgraphIsomorphisms(struct)

def generateProductStructures(family, reactantStructures, maps, forward):
    """
    For a given set of `reactantStructures` and a given set of `maps`,
    generate and return the corresponding product structures. The
    `reactantStructures` parameter should be given in the order the
    reactants are stored in the reaction family template. The `maps`
    parameter is a list of mappings of the top-level tree node of each
    *template* reactant to the corresponding *structure*. This function
    returns the product structures, species, and a boolean that tells
    whether any species are new.
    """

    if not forward: template = family.reverseTemplate
    else:           template = family.forwardTemplate

    # Clear any previous atom labeling from all reactant structures
    for struct in reactantStructures: struct.clearLabeledAtoms()

    # If there are two structures and they are the same, then make a copy
    # of the second one and adjust the second map to point to its atoms
    # This is for the case where A + A --> products
    if len(reactantStructures) == 2 and reactantStructures[0] == reactantStructures[1]:
        reactantStructures[1] = reactantStructures[1].copy(deep=True)
        newMap = {}
        for reactantAtom, templateAtom in maps[1].iteritems():
            index = reactantStructures[0].atoms.index(reactantAtom)
            newMap[reactantStructures[1].atoms[index]] = templateAtom
        maps[1] = newMap
        
    # Tag atoms with labels
    for map in maps:
        for reactantAtom, templateAtom in map.iteritems():
            reactantAtom.label = templateAtom.label

    # Generate the product structures by applying the forward reaction recipe
    try:
        productStructures = family.applyRecipe(reactantStructures, forward=forward)
        if not productStructures: return None
    except rmgpy.data.kinetics.InvalidActionError, e:
        print 'Unable to apply reaction recipe!'
        print 'Reaction family is %s in %s direction' % (family.label, 'forward' if forward else 'reverse')
        print 'Reactant structures are:'
        for struct in reactantStructures:
            print struct.toAdjacencyList()
        raise

    # If there are two product structures, place the one containing '*1' first
    if len(productStructures) == 2:
        if not productStructures[0].containsLabeledAtom('*1') and \
            productStructures[1].containsLabeledAtom('*1'):
            productStructures.reverse()

    # Check that reactant and product structures are allowed in this family
    # If not, then stop
    if family.forbidden is not None:
        for label, struct2 in family.forbidden.iteritems():
            for struct in reactantStructures:
                if struct.isSubgraphIsomorphic(struct2): return None
            for struct in productStructures:
                if struct.isSubgraphIsomorphic(struct2): return None

    return productStructures

def createReaction(family, reactants, reactantStructures, productStructures, isForward):
    """
    Create a :class:`Reaction` object representing the reaction that
    converts the structures in `reactantStructures` corresponding to the
    species in `reactants` to the structures in `productStructures`. The
    atom labels for the reactants should already be known, and they are
    passed in the `reactantAtomLabels` parameter.
    """

    from model import Species, Reaction

    # Make sure the products are in fact different than the reactants
    if len(reactants) == len(productStructures) == 1:
        if any([molecule.isIsomorphic(productStructures[0]) for molecule in reactants[0].molecule]):
            return None
    elif len(reactants) == len(productStructures) == 2:
        if ((any([molecule.isIsomorphic(productStructures[0]) for molecule in reactants[0].molecule]) and
            any([molecule.isIsomorphic(productStructures[1]) for molecule in reactants[1].molecule])) or
            (any([molecule.isIsomorphic(productStructures[1]) for molecule in reactants[0].molecule]) and
            any([molecule.isIsomorphic(productStructures[0]) for molecule in reactants[1].molecule]))):
            return None

    # Convert product structures to product species
    # This does *not* check that the products do not already exist
    products = []
    for molecule in productStructures:
        product = Species(molecule=[molecule])
        products.append(product)

    # Create reaction object
    forward = Reaction(reactants=reactants, products=products, family=family, isForward=isForward)
    reverse = Reaction(reactants=products, products=reactants, family=family, isForward=not isForward)
    forward.reverse = reverse
    reverse.reverse = forward

    # We need to save the reactant and product structures with atom labels so
    # we can generate the kinetics
    # We make copies so the structures aren't trampled on by later actions
    # Once we have the kinetics we can delete these to recover the memory
    forward.reactantMolecules = [s.copy(deep=True) for s in reactantStructures]
    reverse.reactantMolecules = [s.copy(deep=True) for s in productStructures]
    
    # Return the created reaction (forward direction only)
    return forward

def generateReactionsForFamily(reactants, family, model, forward=True):
    """
    Generate all possible reactions involving `reactants`, a list of either one
    or two species to react, that can occur for a single reaction `family`.
    """

    rxnList = []; speciesList = []

    if forward: template = family.forwardTemplate
    elif not forward and family.reverseTemplate is not None: template = family.reverseTemplate
    else: return rxnList, speciesList

    # Unimolecular reactants: A --> products
    if len(reactants) == 1 and len(template.reactants) == 1:

        # Iterate over all resonance isomers of the reactant
        for molecule in reactants[0].molecule:

            ismatch, mappings = matchReactantToTemplate(molecule, template.reactants[0], family)
            if ismatch:
                for map in mappings:
                    reactantAtomLabels = [{}]
                    for atom1, atom2 in map.iteritems():
                        reactantAtomLabels[0][atom1.label] = atom2

                    reactantStructures = [molecule]
                    productStructures = generateProductStructures(family, reactantStructures, [map], forward)
                    if productStructures:
                        rxn = createReaction(family, reactants, reactantStructures, productStructures, forward)
                        if rxn: rxnList.append(rxn)
        
    # Bimolecular reactants: A + B --> products
    elif len(reactants) == 2 and len(template.reactants) == 2:

        moleculesA = reactants[0].molecule
        moleculesB = reactants[1].molecule

        # Iterate over all resonance isomers of the reactant
        for moleculeA in moleculesA:
            for moleculeB in moleculesB:

                # Reactants stored as A + B
                ismatchA, mappingsA = matchReactantToTemplate(moleculeA, template.reactants[0], family)
                ismatchB, mappingsB = matchReactantToTemplate(moleculeB, template.reactants[1], family)

                # Iterate over each pair of matches (A, B)
                if ismatchA and ismatchB:
                    for mapA in mappingsA:
                        for mapB in mappingsB:

                            reactantAtomLabels = [{},{}]
                            for atom1, atom2 in mapA.iteritems():
                                reactantAtomLabels[0][atom1.label] = atom2
                            for atom1, atom2 in mapB.iteritems():
                                reactantAtomLabels[1][atom1.label] = atom2

                            reactantStructures = [moleculeA, moleculeB]
                            productStructures = generateProductStructures(family, reactantStructures, [mapA, mapB], forward)
                            if productStructures:
                                rxn = createReaction(family, reactants, reactantStructures, productStructures, forward)
                                if rxn: rxnList.append(rxn)

                # Only check for swapped reactants if they are different
                if reactants[0] is not reactants[1]:

                    # Reactants stored as B + A
                    ismatchA, mappingsA = matchReactantToTemplate(moleculeA, template.reactants[1], family)
                    ismatchB, mappingsB = matchReactantToTemplate(moleculeB, template.reactants[0], family)

                    # Iterate over each pair of matches (A, B)
                    if ismatchA and ismatchB:
                        for mapA in mappingsA:
                            for mapB in mappingsB:

                                reactantAtomLabels = [{},{}]
                                for atom1, atom2 in mapA.iteritems():
                                    reactantAtomLabels[0][atom1.label] = atom2
                                for atom1, atom2 in mapB.iteritems():
                                    reactantAtomLabels[1][atom1.label] = atom2

                                reactantStructures = [moleculeA, moleculeB]
                                productStructures = generateProductStructures(family, reactantStructures, [mapA, mapB], forward)
                                if productStructures:
                                    rxn = createReaction(family, reactants, reactantStructures, productStructures, forward)
                                    if rxn: rxnList.append(rxn)

    # Check the product species to see if any of them have been encountered before
    # If so, replace them with the existing species
    # If not, formally create a new species
    for rxn in rxnList:
        for i, product in enumerate(rxn.products):
            rxn.products[i], isNew = model.makeNewSpecies(product.molecule[0])
            if isNew: speciesList.append(rxn.products[i])
        # Sort reactants and products
        rxn.reactants.sort()
        rxn.products.sort()

    # Merge duplicate reactions and increment multiplier
    # In this context we already know that the family and the reactants
    # match, so we only need to check the products
    reactionsToRemove = []
    for index1, rxn1 in enumerate(rxnList):
        for index2, rxn2, in enumerate(rxnList[index1+1:]):
            if rxn1.products == rxn2.products:
                rxn1.degeneracy += 1
                if rxn2 not in reactionsToRemove:
                    reactionsToRemove.append(rxn2)
    for rxn2 in reactionsToRemove:
        rxnList.remove(rxn2)

    # For R_Recombination reactions, the multiplier is twice what it should
    # be, so divide those by two
    # This is hardcoding of reaction families!
    if family.label.lower() == 'unimolecular homolysis':
        for rxn in rxnList:
            rxn.degeneracy /= 2

    # This reaction list has only checked for duplicates within itself, not
    # with the global list of reactions
    return rxnList, speciesList

################################################################################

def generateReactionsForPrimary(kineticsDatabase, species, model):
    """
    Generate all possible reactions involving `species`, a list of either one
    or two species to react, that are present within a single primary kinetics
    database `kineticsDatabase`.
    """

    from model import Reaction

    rxnList0 = kineticsDatabase.getReactionList(species)

    # Check the product species to see if any of them have been encountered before
    # If so, replace them with the existing species
    # If not, formally create a new species
    rxnList = []; speciesList = []
    dictionary = kineticsDatabase.database.dictionary
    for rxn in rxnList0:
        forward = Reaction(reactants=rxn.reactants[:], products=rxn.products[:], family=kineticsDatabase, kinetics=rxn.kinetics, isForward=True)
        for i, product in enumerate(forward.products):
            label = dictionary.keys()[dictionary.values().index(product)]
            forward.products[i], isNew = model.makeNewSpecies(product, label=label)
            if isNew: speciesList.append(forward.products[i])
        # Sort reactants and products
        rxn.reactants.sort()
        rxn.products.sort()

        reverse = Reaction(reactants=rxn.products, products=rxn.reactants, family=kineticsDatabase, isForward=False)
        forward.reverse = reverse
        reverse.reverse = forward

    return rxnList, speciesList

def generateReactions(species, model):
    """
    Generate all possible reactions involving `species`, a list of either one
    or two species to react.
    """

    reactionList = []; speciesList = []

    # Don't bother if any or all of the species are marked as nonreactive
    if not all([spec.reactive for spec in species]):
        return reactionList

    log_text = ' + '.join([str(spec) for spec in species])

    logging.info('Generating reactions for %s...'%(log_text))

    # Must use explicit hydrogens for reaction generation
    implicitH = [spec.molecule[0].implicitHydrogens for spec in species]
    for spec in species:
        for molecule in spec.molecule: molecule.makeHydrogensExplicit()

    for kineticsDatabase in rmgpy.data.kinetics.kineticsDatabases:
        if isinstance(kineticsDatabase, rmgpy.data.kinetics.KineticsPrimaryDatabase) and kineticsDatabase.reactionLibrary and not kineticsDatabase.seedMechanism:
            rxnList, specList = generateReactionsForPrimary(kineticsDatabase, species, model)
            speciesList.extend(specList)
            # Formally make the new reactions
            for rxn in rxnList:
                forward = species[0] in rxn.reactants
                r, isNew = model.makeNewReaction(rxn if forward else rxn.reverse)
                reactionList.append(r)
        elif isinstance(kineticsDatabase, rmgpy.data.kinetics.KineticsGroupDatabase):
            for key, family in rmgpy.data.kinetics.kineticsDatabases[-1].families.iteritems():
                for forward in [True, False]:
                    rxnList, specList = generateReactionsForFamily(species, family, model, forward=forward)
                    speciesList.extend(specList)
                    # Formally make the new reactions
                    for rxn in rxnList:
                        r, isNew = model.makeNewReaction(rxn if forward else rxn.reverse)
                        reactionList.append(r)

    # Restore implicit hydrogens if necessary
    for implicit, spec in zip(implicitH, species):
        if implicit:
            for molecule in spec.molecule: molecule.makeHydrogensImplicit()

    # This list contains all generated reactions, not just those that are new
    return reactionList, speciesList
