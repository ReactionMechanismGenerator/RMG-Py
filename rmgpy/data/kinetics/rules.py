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
This module contains functionality for working with kinetics "rate rules",
which provide rate coefficient parameters for various combinations of 
functional groups.
"""
import codecs
import math
import os.path
import re
import warnings
from copy import deepcopy

import numpy as np

from rmgpy.data.base import Database, Entry, get_all_combinations
from rmgpy.data.kinetics.common import save_entry
from rmgpy.exceptions import KineticsError, DatabaseError
from rmgpy.kinetics import ArrheniusEP, Arrhenius, StickingCoefficientBEP, SurfaceArrheniusBEP
from rmgpy.quantity import Quantity, ScalarQuantity
from rmgpy.reaction import Reaction


################################################################################


class KineticsRules(Database):
    """
    A class for working with a set of "rate rules" for a RMG kinetics family.
    """

    def __init__(self, label="", name="", short_desc="", long_desc="", auto_generated=False):
        Database.__init__(self, label=label, name=name, short_desc=short_desc, long_desc=long_desc)
        self.auto_generated = auto_generated

    def __repr__(self):
        return '<KineticsRules "{0}">'.format(self.label)

    def load_entry(
        self,
        index,
        kinetics=None,
        degeneracy=1,
        label="",
        duplicate=False,
        reversible=True,
        reference=None,
        referenceType="",
        shortDesc="",
        longDesc="",
        rank=None,
        nodalDistance=None,
        treeDistances=None,
    ):
        """
        Method for parsing entries in database files.
        Note that these argument names are retained for backward compatibility.
        """

        if isinstance(kinetics, Arrhenius):
            kinetics = kinetics.to_arrhenius_ep()
        entry = Entry(
            index=index,
            label=label,
            # item = reaction,
            data=kinetics,
            reference=reference,
            reference_type=referenceType,
            short_desc=shortDesc,
            long_desc=longDesc.strip(),
            rank=rank,
            nodal_distance=nodalDistance,
        )
        try:
            self.entries[label].append(entry)
        except KeyError:
            self.entries[label] = [entry]
        return entry

    def save_entry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return save_entry(f, entry)

    def get_entries(self):
        """
        Return a list of all of the entries in the rate rules database,
        sorted by index.
        """
        entries = []
        for entry in self.entries.values():
            if isinstance(entry, list):
                entries.extend(entry)
            else:
                entries.append(entry)
        entries.sort(key=lambda x: x.index)
        return entries

    def get_entries_to_save(self):
        """
        Return a sorted list of all of the entries in the rate rules database
        to save.
        """
        return self.get_entries()

    def has_rule(self, template):
        """
        Return ``True`` if a rate rule with the given `template` currently
        exists, or ``False`` otherwise.
        """
        return self.get_rule(template) is not None

    def get_rule(self, template):
        """
        Return the exact rate rule with the given `template`, or ``None`` if no
        corresponding entry exists.
        """
        entries = self.get_all_rules(template)

        if len(entries) == 1:
            return entries[0]
        elif len(entries) > 1:
            # Take the entry with the highest rank (smaller numbers are higher) and smallest index
            # If an entry has rank 0 or None, give it an effective rank of 1000 for sorting
            entries.sort(key=lambda x: (1000 if not x.rank else x.rank, x.index))
            return entries[0]
        else:
            return None

    def get_all_rules(self, template):
        """
        Return all of the exact rate rules with the given `template`. Raises a
        :class:`ValueError` if no corresponding entry exists.
        """
        entries = []
        template_labels = ";".join([group.label for group in template])
        try:
            entries.extend(self.entries[template_labels])
        except KeyError:
            pass

        return entries

    def fill_rules_by_averaging_up(self, root_template, already_done, verbose=False):
        """
        Fill in gaps in the kinetics rate rules by averaging child nodes.
        If verbose is set to True, then exact sources of kinetics are saved in the kinetics comments
        (warning: this uses up a lot of memory due to the extensively long comments)
        """
        root_label = ";".join([g.label for g in root_template])

        if root_label in already_done:
            return already_done[root_label]

        # Generate the distance 1 pairings which must be averaged for this root template.
        # The distance 1 template is created by taking the parent node from one or more trees
        # and creating the combinations with children from a single remaining tree.
        # i.e. for some node (A,B), we want to fetch all combinations for the pairing of (A,B's children) and
        # (A's children, B).  For node (A,B,C), we would retrieve all combinations of (A,B,C's children)
        # (A,B's children,C) etc...
        # If a particular node has no children, it is skipped from the children expansion altogether.

        children_list = []
        distance_list = []
        for i, parent in enumerate(root_template):
            # Start with the root template, and replace the ith member with its children
            if parent.children:
                children_set = [[group] for group in root_template]
                children_set[i] = parent.children
                children_list.extend(get_all_combinations(children_set))
                distance_list.extend([k.nodal_distance for k in parent.children])

        if distance_list != []:  # average the minimum distance neighbors
            min_dist = min(distance_list)
            close_children_list = [children_list[i] for i in range(len(children_list)) if distance_list[i] == min_dist]
        else:
            close_children_list = []

        kinetics_list = []
        for template in children_list:
            label = ";".join([g.label for g in template])

            if label in already_done:
                kinetics = already_done[label]
            else:
                kinetics = self.fill_rules_by_averaging_up(template, already_done, verbose)

            if template in close_children_list and kinetics is not None:
                kinetics_list.append([kinetics, template])

        # See if we already have a rate rule for this exact template instead
        # and return it now that we have finished searching its children
        entry = self.get_rule(root_template)

        if entry is not None and entry.rank > 0:
            # We already have a rate rule for this exact template
            # If the entry has rank of zero, then we have so little faith
            # in it that we'd rather use an averaged value if possible
            # Since this entry does not have a rank of zero, we keep its
            # value
            already_done[root_label] = entry.data
            return entry.data

        if len(kinetics_list) > 0:
            if len(kinetics_list) > 1:
                # We found one or more results! Let's average them together
                kinetics = self._get_average_kinetics([k for k, t in kinetics_list])

                if verbose:
                    kinetics.comment = "Average of [{0}]".format(
                        " + ".join(k.comment if k.comment != "" else ";".join(g.label for g in t) for k, t in kinetics_list)
                    )

                else:
                    kinetics.comment = "Average of [{0}]".format(" + ".join(";".join(g.label for g in t) for k, t in kinetics_list))

            else:
                k, t = kinetics_list[0]
                kinetics = deepcopy(k)
                # Even though we are using just a single set of kinetics, it's still considered
                # an average.  It just happens that the other distance 1 children had no data.

                if verbose:
                    kinetics.comment = "Average of [{0}]".format(k.comment if k.comment != "" else ";".join(g.label for g in t))
                else:
                    kinetics.comment = "Average of [{0}]".format(";".join(g.label for g in t))

            entry = Entry(
                index=0,
                label=root_label,
                item=root_template,
                data=kinetics,
                rank=11,  # Indicates this is an averaged estimate
            )
            self.entries[entry.label] = [entry]
            already_done[root_label] = entry.data
            return entry.data

        already_done[root_label] = None
        return None

    def _get_average_kinetics(self, kinetics_list):
        """
        Based on averaging log k. For most complex case:
        k = AT^n * exp(-Ea+alpha*H)
        log k = log(A) * nlog(T) * (-Ea + alpha*H)

        Hence we average n, Ea, and alpha arithmetically, but we
        average log A (geometric average)
        """
        logA = 0.0
        n = 0.0
        E0 = 0.0
        alpha = 0.0
        count = len(kinetics_list)
        for kinetics in kinetics_list:
            logA += math.log10(kinetics.A.value_si)
            n += kinetics.n.value_si
            alpha += kinetics.alpha.value_si
            E0 += kinetics.E0.value_si
        logA /= count
        n /= count
        alpha /= count
        E0 /= count
        Aunits = kinetics_list[0].A.units
        if Aunits == "cm^3/(mol*s)" or Aunits == "cm^3/(molecule*s)" or Aunits == "m^3/(molecule*s)":
            Aunits = "m^3/(mol*s)"
        elif Aunits == "cm^6/(mol^2*s)" or Aunits == "cm^6/(molecule^2*s)" or Aunits == "m^6/(molecule^2*s)":
            Aunits = "m^6/(mol^2*s)"
        elif Aunits == "s^-1" or Aunits == "m^3/(mol*s)" or Aunits == "m^6/(mol^2*s)":
            # they were already in SI
            pass
        elif Aunits in ["m^2/(mol*s)", "cm^2/(mol*s)", "m^2/(molecule*s)", "cm^2/(molecule*s)"]:
            # surface: bimolecular (Langmuir-Hinshelwood)
            Aunits = "m^2/(mol*s)"
        elif Aunits in ["m^5/(mol^2*s)", "cm^5/(mol^2*s)", "m^5/(molecule^2*s)", "cm^5/(molecule^2*s)"]:
            # surface: dissociative adsorption
            Aunits = "m^5/(mol^2*s)"
        elif Aunits == "":
            # surface: sticking coefficient
            pass
        else:
            raise Exception("Invalid units {0} for averaging kinetics.".format(Aunits))

        if type(kinetics) not in [ArrheniusEP, SurfaceArrheniusBEP, StickingCoefficientBEP]:
            raise Exception("Invalid kinetics type {0!r} for {1!r}.".format(type(kinetics), self))

        averaged_kinetics = type(kinetics)(
            A=(10**logA, Aunits),
            n=n,
            alpha=alpha,
            E0=(E0 * 0.001, "kJ/mol"),
        )
        return averaged_kinetics

    def estimate_kinetics(self, template, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template` using rate rules.

        Returns a tuple (kinetics, entry) where `entry` is the database
        entry used to determine the kinetics only if it is an exact match,
        and is None if some averaging or use of a parent node took place.
        """
        entry = self.get_rule(template)

        if self.auto_generated:
            entry0 = entry
            while entry.parent is not None:
                parent = entry.parent
                err_parent = abs(parent.data.uncertainty.data_mean + parent.data.uncertainty.mu - entry.data.uncertainty.data_mean) + sqrt(
                    2.0 * parent.data.uncertainty.var / pi
                )
                err_entry = abs(entry.data.uncertainty.mu) + sqrt(2.0 * entry.data.uncertainty.var / pi)
                if err_entry > err_parent:
                    entry = entry.parent

            kinetics = deepcopy(entry.data)
            if entry0 == entry:
                kinetics.comment = "Estimated from node {}".format(entry.label)
                kinetics.A.value_si *= degeneracy
                if degeneracy > 1:
                    kinetics.comment += "\n"
                    kinetics.comment += "Multiplied by reaction path degeneracy {0}".format(degeneracy)
                return kinetics, entry
            else:
                kinetics.comment = "Matched node {}\n".format(entry0.label)
                kinetics.comment += "Estimated from node {}".format(entry.label)
                kinetics.A.value_si *= degeneracy
                if degeneracy > 1:
                    kinetics.comment += "\n"
                    kinetics.comment += "Multiplied by reaction path degeneracy {0}".format(degeneracy)
                return kinetics, None

        original_leaves = get_template_label(template)
        template_list = [template]
        distance_list = [np.zeros(len(template))]
        min_norm = np.inf
        saved_kinetics = []

        if entry is not None and entry.data:
            saved_kinetics = [[deepcopy(entry.data), template]]
            template_list = []
            min_norm = 0

        while len(template_list) > 0:
            kinetics_list = []
            distances = []
            for i, t in enumerate(template_list):
                entry = self.get_rule(t)
                if entry is None:
                    continue
                kinetics = deepcopy(entry.data)
                kinetics_list.append([kinetics, t])
                distances.append(distance_list[i])

            if len(kinetics_list) > 0:
                # Filter the kinetics to use templates with the lowest minimum euclidean distance
                # from the specified template
                norms = [np.linalg.norm(d) for d in distances]
                new_min_norm = min(norms)
                if new_min_norm == min_norm:
                    saved_kinetics.extend([pair for pair, norm in zip(kinetics_list, norms) if norm == min(norms)])
                elif new_min_norm < min_norm:
                    min_norm = new_min_norm
                    saved_kinetics = [pair for pair, norm in zip(kinetics_list, norms) if norm == min(norms)]

            template_list0 = template_list  # keep the old template list
            distance_list0 = distance_list  # keep thge old distance list
            distance_list = []
            template_list = []

            if min_norm > 0:  # filter out stuff too large to be used
                to_delete = []
                norms = [np.linalg.norm(d) for d in distance_list0]
                for i in range(len(template_list0)):
                    if norms[i] > min_norm:
                        to_delete.append(i)
                to_delete.reverse()
                for k in to_delete:
                    del template_list0[k]
                    del distance_list0[k]

            for i, template0 in enumerate(template_list0):
                for index in range(len(template0)):
                    if not template0[index].parent:  # We're at the top-level node in this subtreee
                        continue
                    dist = deepcopy(distance_list0[i])
                    t = template0[:]
                    dist[index] += t[index].nodal_distance
                    t[index] = t[index].parent

                    if t not in template_list:
                        template_list.append(t)
                        distance_list.append(dist)

            if template_list != [] and min_norm != 0:
                continue

        kinetics_list = remove_identical_kinetics(saved_kinetics)

        if len(kinetics_list) == 0:
            raise KineticsError("Unable to determine kinetics for reaction with template {0} in family " "{1}.".format(template, self.label))

        elif len(kinetics_list) == 1:
            kinetics, t = kinetics_list[0]
            # Check whether the exact rate rule for the original template (most specific
            # leaves) were found or not.
            matched_leaves = get_template_label(t)
            if kinetics.comment:
                kinetics.comment += "\n"
            if matched_leaves == original_leaves:
                if "Average" in kinetics.comment:
                    kinetics.comment += "Estimated using an average"
                else:
                    kinetics.comment += "Exact match found"
            else:
                # Using a more general node to estimate original template
                kinetics.comment += "Estimated using template " + matched_leaves

        else:
            # We found one or more results! Let's average them together
            kinetics = self._get_average_kinetics([k for k, t in kinetics_list])
            # Unlike in the case of a single rule, the verbose comments for averaging are lost unless they are
            # appended in the following lines.  Verbose comments are filtered out in
            # rmgpy.rmg.model.CoreEdgeReactionModel.generate_kinetics
            kinetics.comment = "Average of [{0}]".format(
                " + ".join(k.comment if k.comment != "" else ";".join(g.label for g in t) for k, t in kinetics_list)
            )
            kinetics.comment += "\n"
            # Append standard portion of kinetics comments that appear in non-verbose mode.
            kinetics.comment += "Estimated using average of templates {0}".format(
                " + ".join([get_template_label(t) for k, t in kinetics_list]),
            )

        kinetics.comment += " for rate rule " + original_leaves
        kinetics.comment += "\nEuclidian distance = {}".format(min_norm)
        kinetics.A.value_si *= degeneracy
        if degeneracy > 1:
            kinetics.comment += "\n"
            kinetics.comment += "Multiplied by reaction path degeneracy {0}".format(degeneracy)

        kinetics.comment += "\n"
        kinetics.comment += "family: {0}".format(self.label.replace("/rules", ""))

        return kinetics, (entry if "Exact" in kinetics.comment else None)


def remove_identical_kinetics(k_list):
    """
    removes all identical kinetics entries in k_list
    takes in a list of kinetics entries
    returns the list with the identical kinetics entries removed

    does this based on strings, which should be fine for this specifically, since we shouldn't have any
    identical kinetics entries in the families and all of the identical kinetics should look exactly the same
    """
    out_set = set()
    out_list = []
    for k in k_list:
        sk = str(k)
        if sk in out_set:
            continue
        else:
            out_set.add(sk)
            out_list.append(k)

    return out_list


def get_template_label(template):
    # Get string format of the template in the form "(leaf1,leaf2)"
    return "[{0}]".format(";".join([g.label for g in template]))
