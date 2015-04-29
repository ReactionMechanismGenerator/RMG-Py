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
This module contains functionality for working with kinetics family functional
groups, including support for using group additivity to estimate rate
coefficients.
"""
import logging
import math
import numpy
from copy import deepcopy

from rmgpy.data.base import Database, DatabaseError, Entry, Group, LogicNode, getAllCombinations, makeLogicNode

from rmgpy.kinetics import Arrhenius, ArrheniusEP, KineticsData
from rmgpy.species import Species
from rmgpy.quantity import constants
from .common import KineticsError, UndeterminableKineticsError

################################################################################

class KineticsGroups(Database):
    """
    A class for working with an RMG kinetics family group additivity values. 
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
        return '<KineticsGroups "{0}">'.format(self.label)

    def loadEntry(self, index, label, group, kinetics, reference=None, referenceType='', shortDesc='', longDesc=''):
        if group[0:3].upper() == 'OR{' or group[0:4].upper() == 'AND{' or group[0:7].upper() == 'NOT OR{' or group[0:8].upper() == 'NOT AND{':
            item = makeLogicNode(group)
        else:
            item = Group().fromAdjacencyList(group)
            
        if label in self.entries:
            raise DatabaseError("Duplicate group name {label} found in kinetics groups for {family} family.".format(label=label,family=self.label))
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = item,
            data = kinetics,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
        )

    def saveGroups(self, path, entryName='entry'):
        """
        Save the current database to the file at location `path` on disk. The
        optional `entryName` parameter specifies the identifier used for each
        data entry.
        """

        entries = self.getEntriesToSave()
                
        # Write the header
        f = codecs.open(path, 'w', 'utf-8')
        f.write('#!/usr/bin/env python\n')
        f.write('# encoding: utf-8\n\n')
        f.write('name = "{0}/groups"\n'.format(self.name))
        f.write('shortDesc = u"{0}"\n'.format(self.shortDesc))
        f.write('longDesc = u"""\n')
        f.write(self.longDesc)
        f.write('\n"""\n\n')

        # Write the template
        f.write('template(reactants=[{0}], products=[{1}], ownReverse={2})\n\n'.format(
            ', '.join(['"{0}"'.format(entry.label) for entry in self.forwardTemplate.reactants]),
            ', '.join(['"{0}"'.format(entry.label) for entry in self.forwardTemplate.products]),
            self.ownReverse))

        # Write reverse name
        if not self.ownReverse:
            f.write('reverse = "{0}"\n\n'.format(self.reverse))

        # Write the recipe
        f.write('recipe(actions=[\n')
        for action in self.forwardRecipe.actions:
            f.write('    {0!r},\n'.format(action))
        f.write('])\n\n')

        # Save the entries
        for entry in entries:
            self.saveEntry(f, entry)

        # Write the tree
        if len(self.groups.top) > 0:
            f.write('tree(\n')
            f.write('"""\n')
            f.write(self.generateOldTree(self.groups.top, 1))
            f.write('"""\n')
            f.write(')\n\n')

        # Save forbidden structures, if present
        if self.forbidden is not None:
            entries = self.forbidden.entries.values()
            entries.sort(key=lambda x: x.label)
            for entry in entries:
                self.forbidden.saveEntry(f, entry, name='forbidden')
    
        f.close()
        
    def getReactionTemplate(self, reaction):
        """
        For a given `reaction` with properly-labeled :class:`Molecule` objects
        as the reactants, determine the most specific nodes in the tree that
        describe the reaction.
        """

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
#            print 'len(template):', len(template)
#            print 'len(forwardTemplate):', len(forwardTemplate)
            msg = 'Unable to find matching template for reaction {0} in reaction family {1}.'.format(str(reaction), str(self)) 
            msg += 'Trying to match {0} but matched {1}'.format(str(forwardTemplate),str(template))
#            print 'reactants'
#            for reactant in reaction.reactants:
#                print reactant.toAdjacencyList() + '\n'
#            print 'products'
#            for product in reaction.products:
#                print product.toAdjacencyList() + '\n'
            raise UndeterminableKineticsError(reaction, message=msg)

        return template

    def estimateKineticsUsingGroupAdditivity(self, template, referenceKinetics, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template` using group additivity.
        """

        # Start with the generic kinetics of the top-level nodes
        # Make a copy so we don't modify the original
        kinetics = deepcopy(referenceKinetics)
        
        # Now add in more specific corrections if possible
        for node in template:
            entry = node
            comment_line = "Matched node "
            while entry.data is None and entry not in self.top:
                # Keep climbing tree until you find a (non-top) node with data.
                comment_line += "{0} >> ".format(entry.label)
                entry = entry.parent
            if entry.data is not None and entry not in self.top:
                kinetics = self.__multiplyKineticsData(kinetics, entry.data)
                comment_line += "{0} ({1})".format(entry.label, entry.longDesc.split('\n')[0])
            elif entry in self.top:
                comment_line += "{0} (Top node)".format(entry.label)
            kinetics.comment += comment_line + '\n'

        # Also include reaction-path degeneracy
        if isinstance(kinetics, KineticsData):
            kinetics.kdata.value_si *= degeneracy
        elif isinstance(kinetics, Arrhenius):
            kinetics.A.value_si *= degeneracy
        elif kinetics is not None:
            raise KineticsError('Unexpected kinetics type "{0}" encountered while generating kinetics from group values.'.format(kinetics.__class__))
        kinetics.comment += "Multiplied by reaction path degeneracy {0}".format(degeneracy)
        
        return kinetics

    def __multiplyKineticsData(self, kinetics1, kinetics2):
        """
        Multiply two kinetics objects `kinetics1` and `kinetics2` of the same
        class together, returning their product as a new kinetics object of 
        that class. Currently this only works for :class:`KineticsData` or
        :class:`Arrhenius` objects.
        """
        if isinstance(kinetics1, KineticsData) and isinstance(kinetics2, KineticsData):
            if len(kinetics1.Tdata.value_si) != len(kinetics2.Tdata.value_si) or any([T1 != T2 for T1, T2 in zip(kinetics1.Tdata.value_si, kinetics2.Tdata.value_si)]):
                raise KineticsError('Cannot add these KineticsData objects due to their having different temperature points.')
            kinetics = KineticsData(
                Tdata = (kinetics1.Tdata.value, kinetics2.Tdata.units),
                kdata = (kinetics1.kdata.value * kinetics2.kdata.value, kinetics1.kdata.units),
            )
        elif isinstance(kinetics1, Arrhenius) and isinstance(kinetics2, Arrhenius):
            assert kinetics1.A.units == kinetics2.A.units
            assert kinetics1.Ea.units == kinetics2.Ea.units
            assert kinetics1.T0.units == kinetics2.T0.units
            assert kinetics1.T0.value == kinetics2.T0.value
            kinetics = Arrhenius(
                A = (kinetics1.A.value * kinetics2.A.value, kinetics1.A.units),
                n = (kinetics1.n.value + kinetics2.n.value, kinetics1.n.units),
                Ea = (kinetics1.Ea.value + kinetics2.Ea.value, kinetics1.Ea.units),
                T0 = (kinetics1.T0.value, kinetics1.T0.units),
            )
        else:
            raise KineticsError('Unable to multiply kinetics types "{0}" and "{1}".'.format(kinetics1.__class__, kinetics2.__class__))
        
        if kinetics1.Tmin is not None and kinetics2.Tmin is not None:
            kinetics.Tmin = kinetics1.Tmin if kinetics1.Tmin.value_si > kinetics2.Tmin.value_si else kinetics2.Tmin
        elif kinetics1.Tmin is not None and kinetics2.Tmin is None:
            kinetics.Tmin = kinetics1.Tmin
        elif kinetics1.Tmin is None and kinetics2.Tmin is not None:
            kinetics.Tmin = kinetics2.Tmin
        
        if kinetics1.Tmax is not None and kinetics2.Tmax is not None:
            kinetics.Tmax = kinetics1.Tmax if kinetics1.Tmax.value_si < kinetics2.Tmax.value_si else kinetics2.Tmax
        elif kinetics1.Tmax is not None and kinetics2.Tmax is None:
            kinetics.Tmax = kinetics1.Tmax
        elif kinetics1.Tmax is None and kinetics2.Tmax is not None:
            kinetics.Tmax = kinetics2.Tmax
        
        if kinetics1.Pmin is not None and kinetics2.Pmin is not None:
            kinetics.Pmin = kinetics1.Pmin if kinetics1.Pmin.value_si > kinetics2.Pmin.value_si else kinetics2.Pmin
        elif kinetics1.Pmin is not None and kinetics2.Pmin is None:
            kinetics.Pmin = kinetics1.Pmin
        elif kinetics1.Pmin is None and kinetics2.Pmin is not None:
            kinetics.Pmin = kinetics2.Pmin
        
        if kinetics1.Pmax is not None and kinetics2.Pmax is not None:
            kinetics.Pmax = kinetics1.Pmax if kinetics1.Pmax.value_si < kinetics2.Pmax.value_si else kinetics2.Pmax
        elif kinetics1.Pmax is not None and kinetics2.Pmax is None:
            kinetics.Pmax = kinetics1.Pmax
        elif kinetics1.Pmax is None and kinetics2.Pmax is not None:
            kinetics.Pmax = kinetics2.Pmax
        
        if kinetics1.comment == '': kinetics.comment = kinetics2.comment
        elif kinetics2.comment == '': kinetics.comment = kinetics1.comment
        else: kinetics.comment = kinetics1.comment + ' + ' + kinetics2.comment
        return kinetics

    def generateGroupAdditivityValues(self, trainingSet, kunits, method='Arrhenius'):
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
            groupEntries.extend(self.descendants(entry))
        
        # Determine a unique list of the groups we will be able to fit parameters for
        groupList = []
        for template, kinetics in trainingSet:
            for group in template:
                if group not in self.top:
                    groupList.append(group)
                    groupList.extend(self.ancestors(group)[:-1])
        groupList = list(set(groupList))
        groupList.sort(key=lambda x: x.index)

        if method == 'KineticsData':
            # Fit a discrete set of k(T) data points by training against k(T) data
            
            Tdata = numpy.array([300,400,500,600,800,1000,1500,2000])
            
            # Initialize dictionaries of fitted group values and uncertainties
            groupValues = {}; groupUncertainties = {}; groupCounts = {}; groupComments = {}
            for entry in groupEntries:
                groupValues[entry] = []
                groupUncertainties[entry] = []
                groupCounts[entry] = []
                groupComments[entry] = set()
            
            # Generate least-squares matrix and vector
            A = []; b = []
            
            kdata = []
            for template, kinetics in trainingSet:
                
                if isinstance(kinetics, (Arrhenius, KineticsData)):
                    kd = [kinetics.getRateCoefficient(T) for T in Tdata]
                elif isinstance(kinetics, ArrheniusEP):
                    kd = [kinetics.getRateCoefficient(T, 0) for T in Tdata]
                else:
                    raise Exception('Unexpected kinetics model of type {0} for template {1}.'.format(kinetics.__class__, template))
                kdata.append(kd)
                    
                # Create every combination of each group and its ancestors with each other
                combinations = []
                for group in template:
                    groups = [group]; groups.extend(self.ancestors(group))
                    combinations.append(groups)
                combinations = getAllCombinations(combinations)
                # Add a row to the matrix for each combination
                for groups in combinations:
                    Arow = [1 if group in groups else 0 for group in groupList]
                    Arow.append(1)
                    brow = [math.log10(k) for k in kd]
                    A.append(Arow); b.append(brow)
                    
                    for group in groups:
                        groupComments[group].add("{0!s}".format(template))
                
            if len(A) == 0:
                logging.warning('Unable to fit kinetics groups for family "{0}"; no valid data found.'.format(self.label))
                return
            A = numpy.array(A)
            b = numpy.array(b)
            kdata = numpy.array(kdata)
            
            x, residues, rank, s = numpy.linalg.lstsq(A, b)
            
            for t, T in enumerate(Tdata):
                
                # Determine error in each group (on log scale)
                stdev = numpy.zeros(len(groupList)+1, numpy.float64)
                count = numpy.zeros(len(groupList)+1, numpy.int)
                
                for index in range(len(trainingSet)):
                    template, kinetics = trainingSet[index]
                    kd = math.log10(kdata[index,t])
                    km = x[-1,t] + sum([x[groupList.index(group),t] for group in template if group in groupList])
                    variance = (km - kd)**2
                    for group in template:
                        groups = [group]; groups.extend(self.ancestors(group))
                        for g in groups:
                            if g not in self.top:
                                ind = groupList.index(g)
                                stdev[ind] += variance
                                count[ind] += 1
                    stdev[-1] += variance
                    count[-1] += 1
                stdev = numpy.sqrt(stdev / (count - 1))
                import scipy.stats
                ci = scipy.stats.t.ppf(0.975, count - 1) * stdev
                
                # Update dictionaries of fitted group values and uncertainties
                for entry in groupEntries:
                    if entry == self.top[0]:
                        groupValues[entry].append(10**x[-1,t])
                        groupUncertainties[entry].append(10**ci[-1])
                        groupCounts[entry].append(count[-1])
                    elif entry in groupList:
                        index = groupList.index(entry)
                        groupValues[entry].append(10**x[index,t])
                        groupUncertainties[entry].append(10**ci[index])
                        groupCounts[entry].append(count[index])
                    else:
                        groupValues[entry] = None
                        groupUncertainties[entry] = None
                        groupCounts[entry] = None
            
            # Store the fitted group values and uncertainties on the associated entries
            for entry in groupEntries:
                if groupValues[entry] is not None:
                    entry.data = KineticsData(Tdata=(Tdata,"K"), kdata=(groupValues[entry],kunits))
                    if not any(numpy.isnan(numpy.array(groupUncertainties[entry]))):
                        entry.data.kdata.uncertainties = numpy.array(groupUncertainties[entry])
                        entry.data.kdata.uncertaintyType = '*|/'
                    entry.shortDesc = "Group additive kinetics."
                    entry.longDesc = "Fitted to {0} rates.\n".format(groupCounts[entry])
                    entry.longDesc += "\n".join(groupComments[entry])
                else:
                    entry.data = None
        
        elif method == 'Arrhenius':
            # Fit Arrhenius parameters (A, n, Ea) by training against k(T) data
            
            Tdata = numpy.array([300,400,500,600,800,1000,1500,2000])
            logTdata = numpy.log(Tdata)
            Tinvdata = 1000. / (constants.R * Tdata)
            
            A = []; b = []
            
            kdata = []
            for template, kinetics in trainingSet:
                
                if isinstance(kinetics, (Arrhenius, KineticsData)):
                    kd = [kinetics.getRateCoefficient(T) for T in Tdata]
                elif isinstance(kinetics, ArrheniusEP):
                    kd = [kinetics.getRateCoefficient(T, 0) for T in Tdata]
                else:
                    raise Exception('Unexpected kinetics model of type {0} for template {1}.'.format(kinetics.__class__, template))
                kdata.append(kd)
                
                # Create every combination of each group and its ancestors with each other
                combinations = []
                for group in template:
                    groups = [group]; groups.extend(self.ancestors(group))
                    combinations.append(groups)
                combinations = getAllCombinations(combinations)
                
                # Add a row to the matrix for each combination at each temperature
                for t, T in enumerate(Tdata):
                    logT = logTdata[t]
                    Tinv = Tinvdata[t]
                    for groups in combinations:
                        Arow = []
                        for group in groupList:
                            if group in groups:
                                Arow.extend([1,logT,-Tinv])
                            else:
                                Arow.extend([0,0,0])
                        Arow.extend([1,logT,-Tinv])
                        brow = math.log(kd[t])
                        A.append(Arow); b.append(brow)
            
            if len(A) == 0:
                logging.warning('Unable to fit kinetics groups for family "{0}"; no valid data found.'.format(self.label))
                return
            A = numpy.array(A)
            b = numpy.array(b)
            kdata = numpy.array(kdata)
            
            x, residues, rank, s = numpy.linalg.lstsq(A, b)
            
            # Store the results
            self.top[0].data = Arrhenius(
                A = (math.exp(x[-3]),kunits),
                n = x[-2],
                Ea = (x[-1],"kJ/mol"),
                T0 = (1,"K"),
            )
            for i, group in enumerate(groupList):
                group.data = Arrhenius(
                    A = (math.exp(x[3*i]),kunits),
                    n = x[3*i+1],
                    Ea = (x[3*i+2],"kJ/mol"),
                    T0 = (1,"K"),
                )
        
        elif method == 'Arrhenius2':
            # Fit Arrhenius parameters (A, n, Ea) by training against (A, n, Ea) values
            
            A = []; b = []
            
            for template, kinetics in trainingSet:
                
                # Create every combination of each group and its ancestors with each other
                combinations = []
                for group in template:
                    groups = [group]; groups.extend(self.ancestors(group))
                    combinations.append(groups)
                combinations = getAllCombinations(combinations)
                        
                # Add a row to the matrix for each parameter
                if isinstance(kinetics, Arrhenius) or (isinstance(kinetics, ArrheniusEP) and kinetics.alpha.value_si == 0):
                    for groups in combinations:
                        Arow = []
                        for group in groupList:
                            if group in groups:
                                Arow.append(1)
                            else:
                                Arow.append(0)
                        Arow.append(1)
                        Ea = kinetics.E0.value_si if isinstance(kinetics, ArrheniusEP) else kinetics.Ea.value_si
                        brow = [math.log(kinetics.A.value_si), kinetics.n.value_si, Ea / 1000.]
                        A.append(Arow); b.append(brow)
            
            if len(A) == 0:
                logging.warning('Unable to fit kinetics groups for family "{0}"; no valid data found.'.format(self.label))
                return
            A = numpy.array(A)
            b = numpy.array(b)
            
            x, residues, rank, s = numpy.linalg.lstsq(A, b)
            
            # Store the results
            self.top[0].data = Arrhenius(
                A = (math.exp(x[-1,0]),kunits),
                n = x[-1,1],
                Ea = (x[-1,2],"kJ/mol"),
                T0 = (1,"K"),
            )
            for i, group in enumerate(groupList):
                group.data = Arrhenius(
                    A = (math.exp(x[i,0]),kunits),
                    n = x[i,1],
                    Ea = (x[i,2],"kJ/mol"),
                    T0 = (1,"K"),
                )
        
        # Add a note to the history of each changed item indicating that we've generated new group values
        changed = False
        for label, entry in self.entries.items():
            if entry.data is not None and old_entries.has_key(label):
                if (isinstance(entry.data, KineticsData) and 
                    isinstance(old_entries[label], KineticsData) and
                    len(entry.data.kdata.value_si) == len(old_entries[label].kdata.value_si) and
                    all(abs(entry.data.kdata.value_si / old_entries[label].kdata.value_si - 1) < 0.01)):
                    #print "New group values within 1% of old."
                    pass
                elif (isinstance(entry.data, Arrhenius) and 
                    isinstance(old_entries[label], Arrhenius) and
                    abs(entry.data.A.value_si / old_entries[label].A.value_si - 1) < 0.01 and
                    abs(entry.data.n.value_si / old_entries[label].n.value_si - 1) < 0.01 and
                    abs(entry.data.Ea.value_si / old_entries[label].Ea.value_si - 1) < 0.01 and
                    abs(entry.data.T0.value_si / old_entries[label].T0.value_si - 1) < 0.01):
                    #print "New group values within 1% of old."
                    pass
                else:
                    changed = True
                    break
            else:
                changed = True
                break
        
        return changed
