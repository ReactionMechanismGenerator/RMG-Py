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
This module contains functionality for working with kinetics depositories.
"""

from rmgpy.data.base import Database, Entry
from rmgpy.molecule import Molecule

from rmgpy.reaction import Reaction, ReactionError
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
        return self.depository

################################################################################

class KineticsDepository(Database):
    """
    A class for working with an RMG kinetics depository. Each depository 
    corresponds to a reaction family (a :class:`KineticsFamily` object). Each
    entry in a kinetics depository involves a reaction defined either by a
    real reactant and product species (as in a kinetics library).
    """

    def __init__(self, label='', name='', shortDesc='', longDesc='', recommended=False):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc, recommended=recommended)

    def __repr__(self):
        return '<KineticsDepository "{0}">'.format(self.label)

    def loadEntry(self,
                  index,
                  reactant1=None,
                  reactant2=None,
                  reactant3=None,
                  product1=None,
                  product2=None,
                  product3=None,
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
                  history=None
                  ):
        
        reactants = [Molecule().fromAdjacencyList(reactant1)]
        if reactant2 is not None: reactants.append(Molecule().fromAdjacencyList(reactant2))
        if reactant3 is not None: reactants.append(Molecule().fromAdjacencyList(reactant3))

        products = [Molecule().fromAdjacencyList(product1)]
        if product2 is not None: products.append(Molecule().fromAdjacencyList(product2))
        if product3 is not None: products.append(Molecule().fromAdjacencyList(product3))
        
        reaction = Reaction(reactants=reactants, products=products, degeneracy=degeneracy, duplicate=duplicate, reversible=reversible)
        
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
            history = history or [],
        )
        self.entries['{0:d}:{1}'.format(index,label)] = entry
        return entry

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the kinetics database to the file object `f`.
        """
        return saveEntry(f, entry)
