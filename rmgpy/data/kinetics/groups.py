#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains functionality for working with kinetics family functional
groups, including support for using group additivity to estimate rate
coefficients.
"""
import logging
import math
import warnings
from copy import deepcopy

import numpy as np

import rmgpy.constants as constants
from rmgpy.data.base import Database, Entry, Group, LogicNode, get_all_combinations, make_logic_node
from rmgpy.exceptions import KineticsError, UndeterminableKineticsError, DatabaseError
from rmgpy.kinetics import Arrhenius, ArrheniusEP, KineticsData
from rmgpy.species import Species
from rmgpy.molecule import Molecule
from rmgpy.molecule.fragment import Fragment

# Prior to np 1.14, `np.linalg.lstsq` does not accept None as a value
RCOND = -1 if int(np.__version__.split('.')[1]) < 14 else None


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
                 short_desc='',
                 long_desc='',
                 forwardTemplate=None,
                 forwardRecipe=None,
                 reverseTemplate=None,
                 reverseRecipe=None,
                 forbidden=None
                 ):
        Database.__init__(self, entries, top, label, name, short_desc, long_desc)
        self.reactant_num = 0

    def __repr__(self):
        return '<KineticsGroups "{0}">'.format(self.label)

    def load_entry(self, index, label, group, kinetics, reference=None, referenceType='', shortDesc='', longDesc='',
                   nodalDistance=None):
        """
        Method for parsing entries in database files.
        Note that these argument names are retained for backward compatibility.

        nodal_distance is the distance between a given entry and its parent specified by a float
        """
        if (group[0:3].upper() == 'OR{' or
                group[0:4].upper() == 'AND{' or
                group[0:7].upper() == 'NOT OR{' or
                group[0:8].upper() == 'NOT AND{'):
            item = make_logic_node(group)
        else:
            item = Group().from_adjacency_list(group)

        if label in self.entries:
            raise DatabaseError("Duplicate group name {label} found in kinetics groups for {family} "
                                "family.".format(label=label, family=self.label))
        self.entries[label] = Entry(
            index=index,
            label=label,
            item=item,
            data=kinetics,
            reference=reference,
            reference_type=referenceType,
            short_desc=shortDesc,
            long_desc=longDesc.strip(),
            nodal_distance=nodalDistance
        )

    def get_reaction_template(self, reaction):
        """
        For a given `reaction` with properly-labeled :class:`Molecule` objects
        as the reactants, determine the most specific nodes in the tree that
        describe the reaction.
        """

        # Get forward reaction template and remove any duplicates
        forward_template = self.top[:]

        temporary = []
        symmetric_tree = False
        for entry in forward_template:
            if entry not in temporary:
                temporary.append(entry)
            else:
                # duplicate node found at top of tree
                # eg. R_recombination: ['Y_rad', 'Y_rad']
                if len(forward_template) != 2:
                    raise DatabaseError('Can currently only do symmetric trees with nothing else in them')
                symmetric_tree = True
        forward_template = temporary

        # Descend reactant trees as far as possible
        template = []
        if (len(forward_template) == 1 and len(reaction.reactants) > len(forward_template) and
                self.label.lower().split('/')[0]):
            entry = forward_template[0]
            group = entry.item

            r = None
            for react in reaction.reactants:
                if isinstance(react, Species):
                    react = react.molecule[0]
                if r:
                    if isinstance(r, Molecule) and isinstance(react,Fragment):
                        r = react.merge(r)
                    else:
                        r = r.merge(react)
                else:
                    r = deepcopy(react)

            atoms = r.get_all_labeled_atoms()

            matched_node = self.descend_tree(r, atoms, root=entry, strict=True)

            if matched_node is not None:
                template.append(matched_node)

        else:
            for entry in forward_template:
                # entry is a top-level node that should be matched
                group = entry.item

                # Identify the atom labels in a group if it is not a logical node
                atom_list = []
                if not isinstance(entry.item, LogicNode):
                    atom_list = group.get_all_labeled_atoms()

                for reactant in reaction.reactants:
                    if isinstance(reactant, Species):
                        reactant = reactant.molecule[0]
                    # Match labeled atoms
                    # Check that this reactant has each of the atom labels in this group.
                    # If it is a LogicNode, the atom_list is empty and
                    # it will proceed directly to the descend_tree step.
                    if not all([reactant.contains_labeled_atom(label) for label in atom_list]):
                        continue  # don't try to match this structure - the atoms aren't there!
                    # Match structures
                    atoms = reactant.get_all_labeled_atoms()
                    # Descend the tree, making sure to match atomlabels exactly using strict = True
                    matched_node = self.descend_tree(reactant, atoms, root=entry, strict=True)
                    if matched_node is not None:
                        template.append(matched_node)
                    # else:
                    #    logging.warning("Couldn't find match for {0} in {1}".format(entry,atom_list))
                    #    logging.warning(reactant.to_adjacency_list())

            # Get fresh templates (with duplicate nodes back in)
            forward_template = self.top[:]

        # Check that we were able to match the template.
        # template is a list of the actual matched nodes
        # forward_template is a list of the top level nodes that should be matched
        if len(template) != len(forward_template):
            msg = 'Unable to find matching template for reaction {0} in reaction family {1}.'.format(str(reaction),
                                                                                                     str(self))
            msg += 'Trying to match {0} but matched {1}'.format(str(forward_template), str(template))
            raise UndeterminableKineticsError(reaction, message=msg)

        return template

    def estimate_kinetics_using_group_additivity(self, template, reference_kinetics, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template` using group additivity.
        
        Returns just the kinetics.
        """
        warnings.warn("Group additivity is no longer supported and may be"
                      " removed in version 2.3.", DeprecationWarning)
        # Start with the generic kinetics of the top-level nodes
        # Make a copy so we don't modify the original
        kinetics = deepcopy(reference_kinetics)

        # Now add in more specific corrections if possible
        for node in template:
            entry = node
            comment_line = "Matched node "
            while entry.data is None and entry not in self.top:
                # Keep climbing tree until you find a (non-top) node with data.
                comment_line += "{0} >> ".format(entry.label)
                entry = entry.parent
            if entry.data is not None and entry not in self.top:
                kinetics = self._multiply_kinetics_data(kinetics, entry.data)
                comment_line += "{0} ({1})".format(entry.label, entry.long_desc.split('\n')[0])
            elif entry in self.top:
                comment_line += "{0} (Top node)".format(entry.label)
            kinetics.comment += comment_line + '\n'

        # Also include reaction-path degeneracy

        kinetics.change_rate(degeneracy)

        kinetics.comment += "Multiplied by reaction path degeneracy {0}".format(degeneracy)

        return kinetics

    def _multiply_kinetics_data(self, kinetics1, kinetics2):
        """
        Multiply two kinetics objects `kinetics1` and `kinetics2` of the same
        class together, returning their product as a new kinetics object of 
        that class. Currently this only works for :class:`KineticsData`, :class:`ArrheniusEP` or
        :class:`Arrhenius` objects.
        """
        if isinstance(kinetics1, KineticsData) and isinstance(kinetics2, KineticsData):
            if (len(kinetics1.Tdata.value_si) != len(kinetics2.Tdata.value_si) or
                    any([T1 != T2 for T1, T2 in zip(kinetics1.Tdata.value_si, kinetics2.Tdata.value_si)])):
                raise KineticsError('Cannot add these KineticsData objects due to '
                                    'their having different temperature points.')
            kinetics = KineticsData(
                Tdata=(kinetics1.Tdata.value, kinetics2.Tdata.units),
                kdata=(kinetics1.kdata.value * kinetics2.kdata.value, kinetics1.kdata.units),
            )
        elif isinstance(kinetics1, Arrhenius) and isinstance(kinetics2, Arrhenius):
            assert kinetics1.A.units == kinetics2.A.units
            assert kinetics1.T0.units == kinetics2.T0.units
            assert kinetics1.T0.value == kinetics2.T0.value
            kinetics = Arrhenius(
                A=(kinetics1.A.value * kinetics2.A.value, kinetics1.A.units),
                n=(kinetics1.n.value + kinetics2.n.value, kinetics1.n.units),
                Ea=(kinetics1.Ea.value_si + kinetics2.Ea.value_si, 'J/mol'),
                T0=(kinetics1.T0.value, kinetics1.T0.units),
            )
        elif isinstance(kinetics1, ArrheniusEP) and isinstance(kinetics2, ArrheniusEP):
            assert kinetics1.A.units == kinetics2.A.units
            kinetics = ArrheniusEP(
                A=(kinetics1.A.value * kinetics2.A.value, kinetics1.A.units),
                n=(kinetics1.n.value + kinetics2.n.value, kinetics1.n.units),
                alpha=kinetics1.alpha + kinetics2.alpha,
                E0=(kinetics1.E0.value_si + kinetics2.E0.value_si, 'J/mol'),
            )
        elif isinstance(kinetics1, Arrhenius) and isinstance(kinetics2, ArrheniusEP):
            assert kinetics1.A.units == kinetics2.A.units
            assert kinetics1.T0.units == 'K'
            assert kinetics1.T0.value == 1.0
            kinetics = ArrheniusEP(
                A=(kinetics1.A.value * kinetics2.A.value, kinetics1.A.units),
                n=(kinetics1.n.value + kinetics2.n.value, kinetics1.n.units),
                alpha=kinetics2.alpha,
                E0=(kinetics1.Ea.value_si + kinetics2.E0.value_si, 'J/mol'),
            )
        elif isinstance(kinetics1, ArrheniusEP) and isinstance(kinetics2, Arrhenius):
            assert kinetics1.A.units == kinetics2.A.units
            assert 'K' == kinetics2.T0.units
            assert 1.0 == kinetics2.T0.value
            kinetics = ArrheniusEP(
                A=(kinetics1.A.value * kinetics2.A.value, kinetics1.A.units),
                n=(kinetics1.n.value + kinetics2.n.value, kinetics1.n.units),
                alpha=kinetics1.alpha,
                E0=(kinetics1.E0.value_si + kinetics2.Ea.value_si, 'J/mol'),
            )
        else:
            raise KineticsError('Unable to multiply kinetics types "{0}" and '
                                '"{1}".'.format(kinetics1.__class__, kinetics2.__class__))

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

        if kinetics1.comment == '':
            kinetics.comment = kinetics2.comment
        elif kinetics2.comment == '':
            kinetics.comment = kinetics1.comment
        else:
            kinetics.comment = kinetics1.comment + ' + ' + kinetics2.comment
        return kinetics

    def generate_group_additivity_values(self, training_set, kunits, method='Arrhenius'):
        """
        Generate the group additivity values using the given `training_set`,
        a list of 2-tuples of the form ``(template, kinetics)``. You must also
        specify the `kunits` for the family and the `method` to use when
        generating the group values. Returns ``True`` if the group values have
        changed significantly since the last time they were fitted, or ``False``
        otherwise.
        """
        warnings.warn("Group additivity is no longer supported and may be"
                      " removed in version 2.3.", DeprecationWarning)
        # keep track of previous values so we can detect if they change
        old_entries = dict()
        for label, entry in self.entries.items():
            if entry.data is not None:
                old_entries[label] = entry.data

        # Determine a complete list of the entries in the database, sorted as in the tree
        group_entries = self.top[:]
        for entry in self.top:
            group_entries.extend(self.descendants(entry))

        # Determine a unique list of the groups we will be able to fit parameters for
        group_list = []
        for template, kinetics in training_set:
            for group in template:
                if group not in self.top:
                    group_list.append(group)
                    group_list.extend(self.ancestors(group)[:-1])
        group_list = list(set(group_list))
        group_list.sort(key=lambda x: x.index)

        if method == 'KineticsData':
            # Fit a discrete set of k(T) data points by training against k(T) data

            Tdata = np.array([300, 400, 500, 600, 800, 1000, 1500, 2000])

            # Initialize dictionaries of fitted group values and uncertainties
            group_values = {}
            group_uncertainties = {}
            group_counts = {}
            group_comments = {}
            for entry in group_entries:
                group_values[entry] = []
                group_uncertainties[entry] = []
                group_counts[entry] = []
                group_comments[entry] = set()

            # Generate least-squares matrix and vector
            A = []
            b = []

            kdata = []
            for template, kinetics in training_set:

                if isinstance(kinetics, (Arrhenius, KineticsData)):
                    kd = [kinetics.get_rate_coefficient(T) for T in Tdata]
                elif isinstance(kinetics, ArrheniusEP):
                    kd = [kinetics.get_rate_coefficient(T, 0) for T in Tdata]
                else:
                    raise TypeError('Unexpected kinetics model of type {0} for template '
                                    '{1}.'.format(kinetics.__class__, template))
                kdata.append(kd)

                # Create every combination of each group and its ancestors with each other
                combinations = []
                for group in template:
                    groups = [group]
                    groups.extend(self.ancestors(group))
                    combinations.append(groups)
                combinations = get_all_combinations(combinations)
                # Add a row to the matrix for each combination
                for groups in combinations:
                    Arow = [1 if group in groups else 0 for group in group_list]
                    Arow.append(1)
                    brow = [math.log10(k) for k in kd]
                    A.append(Arow)
                    b.append(brow)

                    for group in groups:
                        group_comments[group].add("{0!s}".format(template))

            if len(A) == 0:
                logging.warning('Unable to fit kinetics groups for family "{0}"; '
                                'no valid data found.'.format(self.label))
                return
            A = np.array(A)
            b = np.array(b)
            kdata = np.array(kdata)

            x, residues, rank, s = np.linalg.lstsq(A, b, rcond=RCOND)

            for t, T in enumerate(Tdata):

                # Determine error in each group (on log scale)
                stdev = np.zeros(len(group_list) + 1, np.float64)
                count = np.zeros(len(group_list) + 1, np.int)

                for index in range(len(training_set)):
                    template, kinetics = training_set[index]
                    kd = math.log10(kdata[index, t])
                    km = x[-1, t] + sum([x[group_list.index(group), t] for group in template if group in group_list])
                    variance = (km - kd) ** 2
                    for group in template:
                        groups = [group]
                        groups.extend(self.ancestors(group))
                        for g in groups:
                            if g not in self.top:
                                ind = group_list.index(g)
                                stdev[ind] += variance
                                count[ind] += 1
                    stdev[-1] += variance
                    count[-1] += 1
                stdev = np.sqrt(stdev / (count - 1))
                import scipy.stats
                ci = scipy.stats.t.ppf(0.975, count - 1) * stdev

                # Update dictionaries of fitted group values and uncertainties
                for entry in group_entries:
                    if entry == self.top[0]:
                        group_values[entry].append(10 ** x[-1, t])
                        group_uncertainties[entry].append(10 ** ci[-1])
                        group_counts[entry].append(count[-1])
                    elif entry in group_list:
                        index = group_list.index(entry)
                        group_values[entry].append(10 ** x[index, t])
                        group_uncertainties[entry].append(10 ** ci[index])
                        group_counts[entry].append(count[index])
                    else:
                        group_values[entry] = None
                        group_uncertainties[entry] = None
                        group_counts[entry] = None

            # Store the fitted group values and uncertainties on the associated entries
            for entry in group_entries:
                if group_values[entry] is not None:
                    entry.data = KineticsData(Tdata=(Tdata, "K"), kdata=(group_values[entry], kunits))
                    if not any(np.isnan(np.array(group_uncertainties[entry]))):
                        entry.data.kdata.uncertainties = np.array(group_uncertainties[entry])
                        entry.data.kdata.uncertainty_type = '*|/'
                    entry.short_desc = "Group additive kinetics."
                    entry.long_desc = "Fitted to {0} rates.\n".format(group_counts[entry])
                    entry.long_desc += "\n".join(group_comments[entry])
                else:
                    entry.data = None

        elif method == 'Arrhenius':
            # Fit Arrhenius parameters (A, n, Ea) by training against k(T) data

            Tdata = np.array([300, 400, 500, 600, 800, 1000, 1500, 2000])
            logTdata = np.log(Tdata)
            Tinvdata = 1000. / (constants.R * Tdata)

            A = []
            b = []

            kdata = []
            for template, kinetics in training_set:

                if isinstance(kinetics, (Arrhenius, KineticsData)):
                    kd = [kinetics.get_rate_coefficient(T) for T in Tdata]
                elif isinstance(kinetics, ArrheniusEP):
                    kd = [kinetics.get_rate_coefficient(T, 0) for T in Tdata]
                else:
                    raise TypeError('Unexpected kinetics model of type {0} for template '
                                    '{1}.'.format(kinetics.__class__, template))
                kdata.append(kd)

                # Create every combination of each group and its ancestors with each other
                combinations = []
                for group in template:
                    groups = [group]
                    groups.extend(self.ancestors(group))
                    combinations.append(groups)
                combinations = get_all_combinations(combinations)

                # Add a row to the matrix for each combination at each temperature
                for t, T in enumerate(Tdata):
                    logT = logTdata[t]
                    Tinv = Tinvdata[t]
                    for groups in combinations:
                        Arow = []
                        for group in group_list:
                            if group in groups:
                                Arow.extend([1, logT, -Tinv])
                            else:
                                Arow.extend([0, 0, 0])
                        Arow.extend([1, logT, -Tinv])
                        brow = math.log(kd[t])
                        A.append(Arow)
                        b.append(brow)

            if len(A) == 0:
                logging.warning('Unable to fit kinetics groups for family "{0}"; '
                                'no valid data found.'.format(self.label))
                return
            A = np.array(A)
            b = np.array(b)
            kdata = np.array(kdata)

            x, residues, rank, s = np.linalg.lstsq(A, b, rcond=RCOND)

            # Store the results
            self.top[0].data = Arrhenius(
                A=(math.exp(x[-3]), kunits),
                n=x[-2],
                Ea=(x[-1], "kJ/mol"),
                T0=(1, "K"),
            )
            for i, group in enumerate(group_list):
                group.data = Arrhenius(
                    A=(math.exp(x[3 * i]), kunits),
                    n=x[3 * i + 1],
                    Ea=(x[3 * i + 2], "kJ/mol"),
                    T0=(1, "K"),
                )

        elif method == 'Arrhenius2':
            # Fit Arrhenius parameters (A, n, Ea) by training against (A, n, Ea) values

            A = []
            b = []

            for template, kinetics in training_set:

                # Create every combination of each group and its ancestors with each other
                combinations = []
                for group in template:
                    groups = [group]
                    groups.extend(self.ancestors(group))
                    combinations.append(groups)
                combinations = get_all_combinations(combinations)

                # Add a row to the matrix for each parameter
                if (isinstance(kinetics, Arrhenius) or
                        (isinstance(kinetics, ArrheniusEP) and kinetics.alpha.value_si == 0)):
                    for groups in combinations:
                        Arow = []
                        for group in group_list:
                            if group in groups:
                                Arow.append(1)
                            else:
                                Arow.append(0)
                        Arow.append(1)
                        Ea = kinetics.E0.value_si if isinstance(kinetics, ArrheniusEP) else kinetics.Ea.value_si
                        brow = [math.log(kinetics.A.value_si), kinetics.n.value_si, Ea / 1000.]
                        A.append(Arow)
                        b.append(brow)

            if len(A) == 0:
                logging.warning('Unable to fit kinetics groups for family "{0}"; '
                                'no valid data found.'.format(self.label))
                return
            A = np.array(A)
            b = np.array(b)

            x, residues, rank, s = np.linalg.lstsq(A, b, rcond=RCOND)

            # Store the results
            self.top[0].data = Arrhenius(
                A=(math.exp(x[-1, 0]), kunits),
                n=x[-1, 1],
                Ea=(x[-1, 2], "kJ/mol"),
                T0=(1, "K"),
            )
            for i, group in enumerate(group_list):
                group.data = Arrhenius(
                    A=(math.exp(x[i, 0]), kunits),
                    n=x[i, 1],
                    Ea=(x[i, 2], "kJ/mol"),
                    T0=(1, "K"),
                )

        # Add a note to the history of each changed item indicating that we've generated new group values
        changed = False
        for label, entry in self.entries.items():
            if entry.data is not None and label in old_entries:
                if (isinstance(entry.data, KineticsData) and
                        isinstance(old_entries[label], KineticsData) and
                        len(entry.data.kdata.value_si) == len(old_entries[label].kdata.value_si) and
                        all(abs(entry.data.kdata.value_si / old_entries[label].kdata.value_si - 1) < 0.01)):
                    # New group values within 1% of old
                    pass
                elif (isinstance(entry.data, Arrhenius) and
                        isinstance(old_entries[label], Arrhenius) and
                        abs(entry.data.A.value_si / old_entries[label].A.value_si - 1) < 0.01 and
                        abs(entry.data.n.value_si / old_entries[label].n.value_si - 1) < 0.01 and
                        abs(entry.data.Ea.value_si / old_entries[label].Ea.value_si - 1) < 0.01 and
                        abs(entry.data.T0.value_si / old_entries[label].T0.value_si - 1) < 0.01):
                    # New group values within 1% of old
                    pass
                else:
                    changed = True
                    break
            else:
                changed = True
                break

        return changed
