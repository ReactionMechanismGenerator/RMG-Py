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
