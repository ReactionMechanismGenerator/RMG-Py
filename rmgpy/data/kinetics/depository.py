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
    corresponds to a reaction family (a :class:`KineticsGroups` object). Each
    entry in a kinetics depository involves a reaction defined either by
    real reactant and product species (as in a kinetics library) or a set of
    functional groups (as in a reaction family).
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
                  group1=None,
                  group2=None,
                  group3=None,
                  group4=None,
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
        
        if reactant1 is not None and product1 is not None:
            # The reaction involves real reactants and products
            assert group1 is None and group2 is None and group3 is None and group4 is None
            
            reactants = [Molecule().fromAdjacencyList(reactant1)]
            if reactant2 is not None: reactants.append(Molecule().fromAdjacencyList(reactant2))
            if reactant3 is not None: reactants.append(Molecule().fromAdjacencyList(reactant3))

            products = [Molecule().fromAdjacencyList(product1)]
            if product2 is not None: products.append(Molecule().fromAdjacencyList(product2))
            if product3 is not None: products.append(Molecule().fromAdjacencyList(product3))
            
            reaction = Reaction(reactants=reactants, products=products, degeneracy=degeneracy, duplicate=duplicate, reversible=reversible)
        
        elif group1 is not None:
            # The reaction involves functional groups
            assert reactant1 is None and reactant2 is None and reactant3 is None
            assert product1 is None and product2 is None and product3 is None
            
            reactants = []
            
            if group1[0:3].upper() == 'OR{' or group1[0:4].upper() == 'AND{' or group1[0:7].upper() == 'NOT OR{' or group1[0:8].upper() == 'NOT AND{':
                reactants.append(makeLogicNode(group1))
            else:
                reactants.append(Group().fromAdjacencyList(group1))
            if group2 is not None: 
                if group2[0:3].upper() == 'OR{' or group2[0:4].upper() == 'AND{' or group2[0:7].upper() == 'NOT OR{' or group2[0:8].upper() == 'NOT AND{':
                    reactants.append(makeLogicNode(group2))
                else:
                    reactants.append(Group().fromAdjacencyList(group2))
            if group3 is not None: 
                if group3[0:3].upper() == 'OR{' or group3[0:4].upper() == 'AND{' or group3[0:7].upper() == 'NOT OR{' or group3[0:8].upper() == 'NOT AND{':
                    reactants.append(makeLogicNode(group3))
                else:
                    reactants.append(Group().fromAdjacencyList(group3))
            if group4 is not None: 
                if group4[0:3].upper() == 'OR{' or group4[0:4].upper() == 'AND{' or group4[0:7].upper() == 'NOT OR{' or group4[0:8].upper() == 'NOT AND{':
                    reactants.append(makeLogicNode(group4))
                else:
                    reactants.append(Group().fromAdjacencyList(group4))
            
            reaction = Reaction(reactants=reactants, products=[])
        else:
            raise ValueError("Entry doesn't seem to involve reactants and products, or groups.")
            
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
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return saveEntry(f, entry)

    def processOldLibraryEntry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        thermo database, returning the corresponding kinetics object.
        """
        # This is hardcoding of reaction families!
        label = os.path.split(self.label)[-2]
        if label in ['H_Abstraction',
                     'R_Addition_MultipleBond',
                     'R_Recombination',
                     'Disproportionation',
                     '1+2_Cycloaddition',
                     '2+2_cycloaddition_Cd',
                     '2+2_cycloaddition_CO',
                     '2+2_cycloaddition_CCO',
                     'Diels_alder_addition',
                     '1,2_Insertion',
                     '1,3_Insertion_CO2',
                     '1,3_Insertion_ROR',
                     'R_Addition_COm',
                     'Oa_R_Recombination',
                     'Substitution_O',
                     'SubstitutionS',
                     'R_Addition_CSm',
                     '1,3_Insertion_RSR',
                     ]:
            Aunits = 'cm^3/(mol*s)'
        elif label in ['intra_H_migration',
                      'Birad_recombination',
                      'intra_OH_migration',
                      'HO2_Elimination_from_PeroxyRadical',
                      'Cyclic_Ether_Formation',
                      'Enol_keto_tautomerism',
                      'Intra_R_Add_Exocyclic',
                      'Intra_R_Add_Endocyclic',
                      '1,2-Birad_to_alkene',
                      'Intra_Disproportionation',
                      'Korcek_step1',
                      'Korcek_step2',
                      '1,2_shiftS',
                      'intra_substitutionCS_cyclization',
                      'intra_substitutionCS_isomerization',
                      'intra_substitutionS_cyclization',
                      'intra_substitutionS_isomerization',
                      ]:
            Aunits = 's^-1'
        else:
            raise Exception('Unable to determine preexponential units for old reaction family "{0}".'.format(self.label))

        try:
            Tmin, Tmax = data[0].split('-')
            Tmin = (float(Tmin),"K")
            Tmax = (float(Tmax),"K")
        except ValueError:
            Tmin = (float(data[0]),"K")
            Tmax = None

        A, n, alpha, E0, dA, dn, dalpha, dE0 = data[1:9]
        
        A = float(A)
        if dA[0] == '*':
            A = Quantity(A,Aunits,'*|/',float(dA[1:]))
        else:
            dA = float(dA)
            if dA != 0:
                A = Quantity(A,Aunits,'+|-',dA)
            else:
                A = Quantity(A,Aunits)
        
        n = float(n); dn = float(dn)
        if dn != 0:
            n = Quantity(n,'','+|-',dn)
        else:
            n = Quantity(n,'')
                
        alpha = float(alpha); dalpha = float(dalpha)
        if dalpha != 0:
            alpha = Quantity(alpha,'','+|-',dalpha)
        else:
            alpha = Quantity(alpha,'')
        
        E0 = float(E0); dE0 = float(dE0)
        if dE0 != 0:
            E0 = Quantity(E0,'kcal/mol','+|-',dE0)
        else:
            E0 = Quantity(E0,'kcal/mol')
        
        rank = int(data[9])
        
        return ArrheniusEP(A=A, n=n, alpha=alpha, E0=E0, Tmin=Tmin, Tmax=Tmax), rank

    def loadOldRateRules(self, path, groups, numLabels):
        """
        Load a set of old rate rules for kinetics groups into this depository.
        """
        # Parse the old library
        entries = self.parseOldLibrary(os.path.join(path, 'rateLibrary.txt'), numParameters=10, numLabels=numLabels)
        
        self.entries = {}
        for entry in entries:
            index, label, data, shortDesc = entry
            if isinstance(data, (str,unicode)):
                kinetics = data
                rank = 0
            elif isinstance(data, tuple) and len(data) == 2:
                kinetics, rank = data
            else:
                raise DatabaseError('Unexpected data {0!r} for entry {1!s}.'.format(data, entry))
            reactants = [groups.entries[l].item for l in label.split(';')]
            item = Reaction(reactants=reactants, products=[])
            self.entries['{0:d}:{1}'.format(index,label)] = Entry(
                index = index,
                label = label,
                item = item,
                data = kinetics,
                rank = rank,
                shortDesc = shortDesc
            )
        self.__loadOldComments(path)
    
    def __loadOldComments(self, path):
        """
        Load a set of old comments from the ``comments.txt`` file for the old
        kinetics groups. This function assumes that the groups have already
        been loaded.
        """
        index = 'General' #mops up comments before the first rate ID
        
        re_underline = re.compile('^\-+')
        
        comments = {}
        comments[index] = ''
        
        # Load the comments into a temporary dictionary for now
        # If no comments file then do nothing
        try:
            f = codecs.open(os.path.join(path, 'comments.rst'), 'r', 'utf-8')
        except IOError:
            return
        for line in f:
            match = re_underline.match(line)
            if match:
                index = f.next().strip()
                assert line.rstrip() == f.next().rstrip(), "Overline didn't match underline"
                if not comments.has_key(index):
                    comments[index] = ''
                line = f.next()
            comments[index] += line
        f.close()
        
        # Transfer the comments to the longDesc attribute of the associated entry
        unused = []
        for index, longDesc in comments.iteritems():
            try:
                index = int(index)
            except ValueError:
                unused.append(index)
                
            if isinstance(index, int):
                for entry in self.entries.values():
                    if entry.index == index:
                        entry.longDesc = longDesc
                        break
                #else:
                #    unused.append(str(index))
            
        # Any unused comments are placed in the longDesc attribute of the depository
        self.longDesc = comments['General'] + '\n'
        unused.remove('General')
        for index in unused:
            try:
                self.longDesc += comments[index] + '\n'
            except KeyError:
                import pdb; pdb.set_trace()
                
    def saveOldRateRules(self, path, groups):
        """
        Save a set of old rate rules for kinetics groups from this depository.
        """
        
        # This is hardcoding of reaction families!
        label = os.path.split(self.label)[-2]
        if label in ['H_Abstraction',
                     'R_Addition_MultipleBond',
                     'R_Recombination',
                     'Disproportionation',
                     '1+2_Cycloaddition',
                     '2+2_cycloaddition_Cd',
                     '2+2_cycloaddition_CO',
                     '2+2_cycloaddition_CCO',
                     'Diels_alder_addition',
                     '1,2_Insertion',
                     '1,3_Insertion_CO2',
                     '1,3_Insertion_ROR',
                     'R_Addition_COm',
                     'Oa_R_Recombination',
                     'Substitution_O',
                     'SubstitutionS',
                     'R_Addition_CSm',
                     '1,3_Insertion_RSR',
                     ]:
            factor = 1.0e6
        elif label in ['intra_H_migration',
                      'Birad_recombination',
                      'intra_OH_migration',
                      'HO2_Elimination_from_PeroxyRadical',
                      'Cyclic_Ether_Formation',
                      'Enol_keto_tautomerism',
                      'Intra_R_Add_Exocyclic',
                      'Intra_R_Add_Endocyclic',
                      '1,2-Birad_to_alkene',
                      'Intra_Disproportionation',
                      'Korcek_step1',
                      'Korcek_step2',
                      '1,2_shiftS',
                      'intra_substitutionCS_cyclization',
                      'intra_substitutionCS_isomerization',
                      'intra_substitutionS_cyclization',
                      'intra_substitutionS_isomerization',
                      ]:
            factor = 1.0
        else:
            raise ValueError('Unable to determine preexponential units for old reaction family "{0}".'.format(self.label))

        entries = self.entries.values()
        entries.sort(key=lambda x: x.index)
        
        flib = codecs.open(os.path.join(path, 'rateLibrary.txt'), 'w', 'utf-8')
        flib.write('// The format for the data in this rate library\n')
        flib.write('Arrhenius_EP\n\n')
        
        fcom = codecs.open(os.path.join(path, 'comments.rst'), 'w', 'utf-8')
        fcom.write('-------\n')
        fcom.write('General\n')
        fcom.write('-------\n')
        fcom.write(self.longDesc.strip() + '\n\n')
        
        for entry in entries:
            flib.write('{0:<5d} '.format(entry.index))
            for label in entry.label.split(';'):
                flib.write('{0:<23} '.format(label))
            if entry.data.Tmax is None:
                # Tmin contains string of Trange
                Trange = '{0}    '.format(entry.data.Tmin)
            else:
                Trange = '{0:g}-{1:g}    '.format(entry.data.Tmin.value_si, entry.data.Tmax.value_si)
            flib.write('{0:<12}'.format(Trange))
            flib.write('{0:11.2e} {1:9.2f} {2:9.2f} {3:11.2f} '.format(
                            entry.data.A.value_si * factor,
                            entry.data.n.value_si,
                            entry.data.alpha.value_si,
                            entry.data.E0.value_si / 4184.
                            ))
            if entry.data.A.isUncertaintyMultiplicative():
                flib.write('*{0:<6g} '.format(entry.data.A.uncertainty))
            else:
                flib.write('{0:<7g} '.format(entry.data.A.uncertainty * factor))
            flib.write('{0:6g} {1:6g} {2:6g} '.format(
                            entry.data.n.uncertainty,
                            entry.data.alpha.uncertainty,
                            entry.data.E0.uncertainty / 4184.
                            ))

            if not entry.rank:
                entry.rank = 0
            flib.write(u'    {0:<4d}     {1}\n'.format(entry.rank, entry.shortDesc))
            
            fcom.write('------\n')
            fcom.write('{0}\n'.format(entry.index))
            fcom.write('------\n')
            fcom.write(entry.longDesc.strip() + '\n\n')
            
        flib.close()
        fcom.close()
    
