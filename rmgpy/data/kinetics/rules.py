#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
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
This module contains functionality for working with kinetics "rate rules",
which provide rate coefficient parameters for various combinations of 
functional groups.
"""

import os.path
import re
import codecs
import math
import numpy
from copy import  deepcopy

from rmgpy.data.base import Database, Entry, getAllCombinations

from rmgpy.quantity import Quantity, ScalarQuantity
from rmgpy.reaction import Reaction
from rmgpy.kinetics import ArrheniusEP, Arrhenius
from .common import saveEntry
from rmgpy.exceptions import KineticsError, DatabaseError

################################################################################

class KineticsRules(Database):
    """
    A class for working with a set of "rate rules" for a RMG kinetics family. 
    """
    
    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def __repr__(self):
        return '<KineticsRules "{0}">'.format(self.label)

    def loadEntry(self,
                  index,
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
                  nodalDistance=None,
                  treeDistances=None
                  ):
            
        if isinstance(kinetics,Arrhenius):
            kinetics = kinetics.toArrheniusEP()
        entry = Entry(
            index = index,
            label = label,
            # item = reaction,
            data = kinetics,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            rank = rank,
            nodalDistance=nodalDistance,
        )
        try:
            self.entries[label].append(entry)
        except KeyError:
            self.entries[label] = [entry]
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
        
        # The names of all of the RMG reaction families that are bimolecular
        BIMOLECULAR_KINETICS_FAMILIES = [
            'H_Abstraction',
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
            'lone_electron_pair_bond',
        ]
        
        # The names of all of the RMG reaction families that are unimolecular
        UNIMOLECULAR_KINETICS_FAMILIES = [
            'intra_H_migration',
            'Birad_recombination',
            'intra_OH_migration',
            'HO2_Elimination_from_PeroxyRadical',
            'H_shift_cyclopentadiene',
            'Cyclic_Ether_Formation',
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
            'intra_NO2_ONO_conversion',
            '1,4_Cyclic_birad_scission',
            '1,4_Linear_birad_scission',
            'Intra_Diels_alder',
            'ketoenol',
            'Retroen'
        ]
        # This is hardcoding of reaction families!
        label = os.path.split(self.label)[-2]
        if label in BIMOLECULAR_KINETICS_FAMILIES:
            Aunits = 'cm^3/(mol*s)'
        elif label in UNIMOLECULAR_KINETICS_FAMILIES:
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

    def loadOld(self, path, groups, numLabels):
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
            entry = Entry(
                index = index,
                label = label,
                item = item,
                data = kinetics,
                rank = rank,
                shortDesc = shortDesc
            )
            try:
                self.entries[label].append(entry)
            except KeyError:
                self.entries[label] = [entry]
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
        entries = self.getEntries()
        unused = []
        for index, longDesc in comments.iteritems():
            try:
                index = int(index)
            except ValueError:
                unused.append(index)
                
            if isinstance(index, int):
                for entry in entries:
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
                
    def saveOld(self, path, groups):
        """
        Save a set of old rate rules for kinetics groups from this depository.
        """
        
        # This is hardcoding of reaction families!
        label = os.path.split(self.label)[-2]
        reactionOrder = groups.groups.numReactants
        if reactionOrder == 2:
            factor = 1.0e6
        elif reactionOrder == 1:
            factor = 1.0
        else:
            raise ValueError('Unable to determine preexponential units for old reaction family "{0}".'.format(self.label))

        entries = self.getEntries()
        
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
            line = ''
            for label in entry.label.split(';'):
                line = line + '{0:<23} '.format(label)
            flib.write(line)
            if len(line)>48: # make long lines line up in 10-space columns
                flib.write(' '*(10-len(line)%10))
            if entry.data.Tmax is None:
                if re.match('\d+\-\d+',str(entry.data.Tmin).strip()):
                    # Tmin contains string of Trange
                    Trange = '{0} '.format(entry.data.Tmin)
                elif isinstance(entry.data.Tmin, ScalarQuantity):
                    # Tmin is a temperature. Make range 1 degree either side!
                    Trange = '{0:4g}-{1:g} '.format(entry.data.Tmin.value_si-1, entry.data.Tmin.value_si+1)
                else:
                    # Range is missing, but we have to put something:
                    Trange = '   1-9999 '
            else:
                Trange = '{0:4g}-{1:g} '.format(entry.data.Tmin.value_si, entry.data.Tmax.value_si)
            flib.write('{0:<12}'.format(Trange))
            flib.write('{0:11.2e} {1:9.2f} {2:9.2f} {3:11.2f} '.format(
                            entry.data.A.value_si * factor,
                            entry.data.n.value_si,
                            entry.data.alpha.value_si,
                            entry.data.E0.value_si / 4184.
                            ))
            if entry.data.A.isUncertaintyMultiplicative():
                flib.write('*{0:<6g} '.format(entry.data.A.uncertainty_si))
            else:
                flib.write('{0:<7g} '.format(entry.data.A.uncertainty_si * factor))
            flib.write('{0:6g} {1:6g} {2:6g} '.format(
                            entry.data.n.uncertainty_si,
                            entry.data.alpha.uncertainty_si,
                            entry.data.E0.uncertainty_si / 4184.
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

    def getEntries(self):
        """
        Return a list of all of the entries in the rate rules database,
        sorted by index.
        """
        entries = []
        for e in self.entries.values(): 
            if isinstance(e,list):
                entries.extend(e)
            else:
                entries.append(e)
        entries.sort(key=lambda x: x.index)
        return entries

    def getEntriesToSave(self):
        """
        Return a sorted list of all of the entries in the rate rules database
        to save.
        """
        return self.getEntries()

    def hasRule(self, template):
        """
        Return ``True`` if a rate rule with the given `template` currently 
        exists, or ``False`` otherwise.
        """
        return self.getRule(template) is not None

    def getRule(self, template):
        """
        Return the exact rate rule with the given `template`, or ``None`` if no
        corresponding entry exists.
        """
        entries = self.getAllRules(template)
            
        if len(entries) == 1:
            return entries[0]
        elif len(entries) > 1:
            if any([entry.rank > 0 for entry in entries]):
                entries = [entry for entry in entries if entry.rank > 0]
                entries.sort(key=lambda x: (x.rank, x.index))
                return entries[0]
            else:
                entries.sort(key=lambda x: x.index)
                return entries[0]
        else:
            return None

    def getAllRules(self, template):
        """
        Return all of the exact rate rules with the given `template`. Raises a 
        :class:`ValueError` if no corresponding entry exists.
        """
        entries = []
        templateLabels = ';'.join([group.label for group in template])
        try:
            entries.extend(self.entries[templateLabels])
        except KeyError:
            pass
        
        family = os.path.split(self.label)[0]   # i.e. self.label = 'R_Recombination/rules'
        if family.lower() == 'r_recombination':
            template.reverse()
            templateLabels = ';'.join([group.label for group in template])
            try:
                entries.extend(self.entries[templateLabels])
            except KeyError:
                pass
            template.reverse()
        
        return entries

    def fillRulesByAveragingUp(self, rootTemplate, alreadyDone, verbose=False):
        """
        Fill in gaps in the kinetics rate rules by averaging child nodes.
        If verbose is set to True, then exact sources of kinetics are saved in the kinetics comments
        (warning: this uses up a lot of memory due to the extensively long comments)
        """
        rootLabel = ';'.join([g.label for g in rootTemplate])
        
        if rootLabel in alreadyDone:
            return alreadyDone[rootLabel]

        # Generate the distance 1 pairings which must be averaged for this root template.
        # The distance 1 template is created by taking the parent node from one or more trees
        # and creating the combinations with children from a single remaining tree.  
        # i.e. for some node (A,B), we want to fetch all combinations for the pairing of (A,B's children) and
        # (A's children, B).  For node (A,B,C), we would retrieve all combinations of (A,B,C's children) 
        # (A,B's children,C) etc...  
        # If a particular node has no children, it is skipped from the children expansion altogether.

        childrenList = []
        distanceList = []
        for i, parent in enumerate(rootTemplate):
            # Start with the root template, and replace the ith member with its children
            if parent.children:
                childrenSet = [[group] for group in rootTemplate]
                childrenSet[i] = parent.children
                childrenList.extend(getAllCombinations(childrenSet))
                distanceList.extend([k.nodalDistance for k in parent.children])
                
        if distanceList != []: #average the minimum distance neighbors
            minDist = min(distanceList) 
            closeChildrenList = [childrenList[i] for i in xrange(len(childrenList)) if distanceList[i]==minDist]
        else:
            closeChildrenList = []
            
        kineticsList = []
        for template in childrenList:
            label = ';'.join([g.label for g in template])
            
            if label in alreadyDone:
                kinetics = alreadyDone[label]
            else:
                kinetics = self.fillRulesByAveragingUp(template, alreadyDone, verbose)
            
            if template in closeChildrenList and kinetics is not None:
                kineticsList.append([kinetics, template])
        
        # See if we already have a rate rule for this exact template instead
        # and return it now that we have finished searching its children
        entry = self.getRule(rootTemplate)
        
        if entry is not None and entry.rank > 0:
            # We already have a rate rule for this exact template
            # If the entry has rank of zero, then we have so little faith
            # in it that we'd rather use an averaged value if possible
            # Since this entry does not have a rank of zero, we keep its
            # value
            alreadyDone[rootLabel] = entry.data
            return entry.data
        
        if len(kineticsList) > 0:
            
            if len(kineticsList) > 1:
                # We found one or more results! Let's average them together
                kinetics = self.__getAverageKinetics([k for k, t in kineticsList])
                
                if verbose:
                    kinetics.comment = 'Average of [{0}]'.format(
                         ' + '.join(k.comment if k.comment != '' else ';'.join(g.label for g in t) for k, t in kineticsList))
                
                else:
                    kinetics.comment = 'Average of [{0}]'.format(
                     ' + '.join(';'.join(g.label for g in t) for k, t in kineticsList))

            else:
                k,t = kineticsList[0]
                kinetics = deepcopy(k)
                # Even though we are using just a single set of kinetics, it's still considered
                # an average.  It just happens that the other distance 1 children had no data.
                
                if verbose:
                    kinetics.comment = 'Average of [{0}]'.format(k.comment if k.comment != '' else ';'.join(g.label for g in t))
                else:
                    kinetics.comment = 'Average of [{0}]'.format(';'.join(g.label for g in t))
                

            
            entry = Entry(
                index = 0,
                label = rootLabel,
                item = rootTemplate,
                data = kinetics,
                rank = 10, # Indicates this is an averaged estimate
            )
            self.entries[entry.label] = [entry]
            alreadyDone[rootLabel] = entry.data
            return entry.data
            
        alreadyDone[rootLabel] = None
        return None

    def __getAverageKinetics(self, kineticsList):
        """
        Based on averaging log k. For most complex case:
        k = AT^n * exp(-Ea+alpha*H)
        log k = log(A) * nlog(T) * (-Ea + alpha*H)
        
        Hence we average n, Ea, and alpha arithmetically, but we
        average log A (geometric average) 
        """
        logA = 0.0; n = 0.0; E0 = 0.0; alpha = 0.0
        count = len(kineticsList)
        for kinetics in kineticsList:
            logA += math.log10(kinetics.A.value_si)
            n += kinetics.n.value_si
            alpha += kinetics.alpha.value_si
            E0 += kinetics.E0.value_si
        logA /= count
        n /= count
        alpha /= count
        E0 /= count
        Aunits = kineticsList[0].A.units
        if Aunits == 'cm^3/(mol*s)' or Aunits == 'cm^3/(molecule*s)' or Aunits == 'm^3/(molecule*s)':
            Aunits = 'm^3/(mol*s)'
        elif Aunits == 'cm^6/(mol^2*s)' or Aunits == 'cm^6/(molecule^2*s)' or Aunits == 'm^6/(molecule^2*s)':
            Aunits = 'm^6/(mol^2*s)'
        elif Aunits == 's^-1' or Aunits == 'm^3/(mol*s)' or Aunits == 'm^6/(mol^2*s)':
            pass
        else:
            raise Exception('Invalid units {0} for averaging kinetics.'.format(Aunits))
        averagedKinetics = ArrheniusEP(
            A = (10**logA,Aunits),
            n = n,
            alpha = alpha,
            E0 = (E0*0.001,"kJ/mol"),
        )
        return averagedKinetics
    
    def estimateKinetics(self, template, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template` using rate rules.
        
        Returns a tuple (kinetics, entry) where `entry` is the database
        entry used to determine the kinetics only if it is an exact match,
        and is None if some averaging or use of a parent node took place.
        """
        entry = self.getRule(template)
        
        originalLeaves = getTemplateLabel(template)
        templateList = [template]
        distanceList = [numpy.zeros(len(template))]
        minNorm = numpy.inf
        savedKinetics = []
        
        if entry is not None and entry.data:
            savedKinetics = [[deepcopy(entry.data),template]]
            templateList = []
            minNorm = 0
            
        while len(templateList) > 0:
            
            kineticsList = []
            distances = []
            for i,t in enumerate(templateList):
                entry = self.getRule(t)
                if entry is None: 
                    continue
                kinetics = deepcopy(entry.data)
                kineticsList.append([kinetics, t])
                distances.append(distanceList[i])
            
            
            if len(kineticsList) > 0:                 
                # Filter the kinetics to use templates with the lowest minimum euclidean distance 
                # from the specified template
                norms = [numpy.linalg.norm(d) for d in distances]
                newMinNorm = min(norms)
                if newMinNorm == minNorm:
                    savedKinetics.extend([pair for pair, norm in zip(kineticsList,norms) if norm == min(norms)])
                elif newMinNorm < minNorm:
                    minNorm = newMinNorm
                    savedKinetics = [pair for pair, norm in zip(kineticsList,norms) if norm == min(norms)]
                
            templateList0 = templateList #keep the old template list
            distanceList0 = distanceList #keep thge old distance list
            distanceList = []
            templateList = []
            
            if minNorm > 0:  #filter out stuff too large to be used
                toDelete = []
                norms = [numpy.linalg.norm(d) for d in distanceList0]
                for i in range(len(templateList0)):
                    if norms[i] > minNorm:
                        toDelete.append(i)
                toDelete.reverse()
                for k in toDelete:
                    del templateList0[k]
                    del distanceList0[k]
                        
            
            for i,template0 in enumerate(templateList0):
                for index in xrange(len(template0)):
                    if not template0[index].parent: # We're at the top-level node in this subtreee
                        continue
                    dist = deepcopy(distanceList0[i])
                    t = template0[:]
                    dist[index] += t[index].nodalDistance
                    t[index] = t[index].parent
                     
                    if t not in templateList:
                        templateList.append(t)
                        distanceList.append(dist)

            if templateList != [] and minNorm != 0:
                continue
            
        kineticsList = removeIdenticalKinetics(savedKinetics)
        
        if len(kineticsList) == 0:
            raise KineticsError('Unable to determine kinetics for reaction with template {0} in family {1}.'.format(template, self.label))
            
        elif len(kineticsList) == 1:
            kinetics, t = kineticsList[0]
            # Check whether the exact rate rule for the original template (most specific
            # leaves) were found or not.
            matchedLeaves = getTemplateLabel(t)
            if kinetics.comment:
                kinetics.comment += '\n'
            if matchedLeaves == originalLeaves:
                if 'Average' in kinetics.comment:
                    kinetics.comment += 'Estimated using an average'
                else:
                    kinetics.comment += 'Exact match found' 
            else:
                #Using a more general node to estimate original template
                kinetics.comment +='Estimated using template ' + matchedLeaves
                            
        else:
            # We found one or more results! Let's average them together
            kinetics = self.__getAverageKinetics([k for k, t in kineticsList])
            # Unlike in the case of a single rule, the verbose comments for averaging are lost unless they are 
            # appended in the following lines.  Verbose comments are filtered out in 
            # rmgpy.rmg.model.CoreEdgeReactionModel.generateKinetics
            kinetics.comment = 'Average of [{0}]'.format(
                    ' + '.join(k.comment if k.comment != '' else ';'.join(g.label for g in t) for k, t in kineticsList))
            kinetics.comment +='\n'
            # Append standard portion of kinetics comments that appear in non-verbose mode.
            kinetics.comment += 'Estimated using average of templates {0}'.format(
                        ' + '.join([getTemplateLabel(t) for k, t in kineticsList]),
                    )
                
        kinetics.comment += ' for rate rule ' + originalLeaves
        kinetics.comment += '\nEuclidian distance = {}'.format(minNorm)
        kinetics.A.value_si *= degeneracy
        if degeneracy > 1:
            kinetics.comment += "\n"
            kinetics.comment += "Multiplied by reaction path degeneracy {0}".format(degeneracy)
        
        kinetics.comment += "\n"
        kinetics.comment += "family: {0}".format(self.label.replace('/rules',''))
        
        return kinetics, (entry if 'Exact' in kinetics.comment else None)

def removeIdenticalKinetics(kList):
    """
    removes all identical kinetics entries in kList
    takes in a list of kinetics entries
    returns the list with the identical kinetics entries removed
    
    does this based on strings, which should be fine for this specifically, since we shouldn't have any
    identical kinetics entries in the families and all of the identical kinetics should look exactly the same
    """
    outSet = set()
    outList = []
    for k in kList:
        sk = str(k)
        if sk in outSet:
            continue
        else:
            outSet.add(sk)
            outList.append(k)
            
    return outList

def getTemplateLabel(template):
    # Get string format of the template in the form "(leaf1,leaf2)"
    return '[{0}]'.format(';'.join([g.label for g in template]))
