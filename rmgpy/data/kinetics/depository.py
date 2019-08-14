#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains functionality for working with kinetics depositories.
"""

from rmgpy.data.base import Database, Entry, DatabaseError
import re

from rmgpy.reaction import Reaction
from .common import saveEntry

################################################################################

class DepositoryReaction(Reaction):
    """
    A Reaction object generated from a reaction depository. In addition to the
    usual attributes, this class includes `depository` and `entry` attributes to
    store the library and the entry in that depository that it was created from.
    """

    def __init__(self,
                 index=-1,
                 reactants=None,
                 products=None,
                 specificCollider=None,
                 kinetics=None,
                 reversible=True,
                 transitionState=None,
                 duplicate=False,
                 degeneracy=1,
                 pairs=None,
                 depository=None,
                 family=None,
                 entry=None
                 ):
        Reaction.__init__(self,
                          index=index,
                          reactants=reactants,
                          products=products,
                          specificCollider=specificCollider,
                          kinetics=kinetics,
                          reversible=reversible,
                          transitionState=transitionState,
                          duplicate=duplicate,
                          degeneracy=degeneracy,
                          pairs=pairs
                          )
        self.depository = depository
        self.family = family
        self.entry = entry

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (DepositoryReaction, (self.index,
                                     self.reactants,
                                     self.products,
                                     self.specificCollider,
                                     self.kinetics,
                                     self.reversible,
                                     self.transitionState,
                                     self.duplicate,
                                     self.degeneracy,
                                     self.pairs,
                                     self.depository,
                                     self.family,
                                     self.entry
                                     ))

    def getSource(self):
        """
        Return the database that was the source of this reaction. For a
        DepositoryReaction this should be a KineticsDepository object.
        """
        return self.depository.label

################################################################################

class KineticsDepository(Database):
    """
    A class for working with an RMG kinetics depository. Each depository 
    corresponds to a reaction family (a :class:`KineticsFamily` object). Each
    entry in a kinetics depository involves a reaction defined either by a
    real reactant and product species (as in a kinetics library).
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)
        
    def __str__(self):
        return 'Kinetics Depository {0}'.format(self.label)

    def __repr__(self):
        return '<KineticsDepository "{0}">'.format(self.label)
    
    def load(self, path, local_context=None, global_context=None):
        import os
        Database.load(self, path, local_context, global_context)
        
        # Load the species in the kinetics library
        # Do not generate resonance structures, since training reactions may be written for a specific resonance form
        speciesDict = self.getSpecies(os.path.join(os.path.dirname(path),'dictionary.txt'), resonance=False)
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

            specificCollider = None
            collider = re.search('\(\+[^\)]+\)',reactants)
            if collider is not None:
                collider = collider.group(0) # save string value rather than the object
                assert collider == re.search('\(\+[^\)]+\)',products).group(0), "Third body colliders in reaction {0} in kinetics library {1} are not identical!".format(rxn_string, self.label)
                extraParenthesis = collider.count('(') -1
                for i in range(extraParenthesis):
                    collider += ')' # allow for species like N2(5) or CH2(T)(15) to be read as specific colliders, although currently not implemented in Chemkin. See RMG-Py #1070
                reactants = reactants.replace(collider,'')
                products = products.replace(collider,'')
                if collider.upper().strip() != "(+M)": # the collider is a specific species, not (+M) or (+m)
                    if collider.strip()[2:-1] not in speciesDict: # stripping spaces, '(+' and ')'
                        raise DatabaseError('Collider species {0} in kinetics library {1} is missing from its dictionary.'.format(collider.strip()[2:-1], self.label))
                    specificCollider = speciesDict[collider.strip()[2:-1]]

            for reactant in reactants.split('+'):
                reactant = reactant.strip()
                if reactant not in speciesDict:
                    raise DatabaseError('Species {0} in kinetics depository {1} is missing from its dictionary.'.format(reactant, self.label))
                # Depository reactions should have molecule objects because they are needed in order to descend the
                # tree using `getReactionTemplate()` later, but species objects work because `getReactionTemplate()`
                # will simply pick the first molecule object in `Species().molecule`.
                rxn.reactants.append(speciesDict[reactant])
            for product in products.split('+'):
                product = product.strip()
                if product not in speciesDict:
                    raise DatabaseError('Species {0} in kinetics depository {1} is missing from its dictionary.'.format(product, self.label))
                # Same comment about molecule vs species objects as above.
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
                  specificCollider=None,
                  kinetics=None,
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
        
#        reactants = [Species().fromAdjacencyList(reactant1)]
#        if reactant2 is not None: reactants.append(Species().fromAdjacencyList(reactant2))
#        if reactant3 is not None: reactants.append(Species().fromAdjacencyList(reactant3))
#
#        products = [Species().fromAdjacencyList(product1)]
#        if product2 is not None: products.append(Species().fromAdjacencyList(product2))
#        if product3 is not None: products.append(Species().fromAdjacencyList(product3))
        
        reaction = Reaction(reactants=[], products=[], specificCollider=specificCollider,
          degeneracy=degeneracy, duplicate=duplicate, reversible=reversible)
        
        entry = Entry(
            index = index,
            label = label,
            item = reaction,
            data = kinetics,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            rank = rank,
        )
        assert index not in self.entries
        self.entries[index] = entry
        return entry

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the kinetics database to the file object `f`.
        """
        return saveEntry(f, entry)
