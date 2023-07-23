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
RCOND = -1 if int(np.__version__.split(".")[1]) < 14 else None


################################################################################


class KineticsGroups(Database):
    """
    A class for working with an RMG kinetics family group additivity values.
    """

    def __init__(
        self,
        entries=None,
        top=None,
        label="",
        name="",
        short_desc="",
        long_desc="",
        forwardTemplate=None,
        forwardRecipe=None,
        reverseTemplate=None,
        reverseRecipe=None,
        forbidden=None,
    ):
        Database.__init__(self, entries, top, label, name, short_desc, long_desc)
        self.reactant_num = 0

    def __repr__(self):
        return '<KineticsGroups "{0}">'.format(self.label)

    def load_entry(self, index, label, group, kinetics, reference=None, referenceType="", shortDesc="", longDesc="", nodalDistance=None):
        """
        Method for parsing entries in database files.
        Note that these argument names are retained for backward compatibility.

        nodal_distance is the distance between a given entry and its parent specified by a float
        """
        if group[0:3].upper() == "OR{" or group[0:4].upper() == "AND{" or group[0:7].upper() == "NOT OR{" or group[0:8].upper() == "NOT AND{":
            item = make_logic_node(group)
        else:
            item = Group().from_adjacency_list(group)

        if label in self.entries:
            raise DatabaseError(
                "Duplicate group name {label} found in kinetics groups for {family} " "family.".format(label=label, family=self.label)
            )
        self.entries[label] = Entry(
            index=index,
            label=label,
            item=item,
            data=kinetics,
            reference=reference,
            reference_type=referenceType,
            short_desc=shortDesc,
            long_desc=longDesc.strip(),
            nodal_distance=nodalDistance,
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
                    raise DatabaseError("Can currently only do symmetric trees with nothing else in them")
                symmetric_tree = True
        forward_template = temporary

        # Descend reactant trees as far as possible
        template = []
        if len(forward_template) == 1 and len(reaction.reactants) > len(forward_template) and self.label.lower().split("/")[0]:
            entry = forward_template[0]
            group = entry.item

            r = None
            for react in reaction.reactants:
                if isinstance(react, Species):
                    react = react.molecule[0]
                if r:
                    if isinstance(r, Molecule) and isinstance(react, Fragment):
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
            msg = "Unable to find matching template for reaction {0} in reaction family {1}.".format(str(reaction), str(self))
            msg += "Trying to match {0} but matched {1}".format(str(forward_template), str(template))
            raise UndeterminableKineticsError(reaction, message=msg)

        return template

    def _multiply_kinetics_data(self, kinetics1, kinetics2):
        """
        Multiply two kinetics objects `kinetics1` and `kinetics2` of the same
        class together, returning their product as a new kinetics object of
        that class. Currently this only works for :class:`KineticsData`, :class:`ArrheniusEP` or
        :class:`Arrhenius` objects.
        """
        if isinstance(kinetics1, KineticsData) and isinstance(kinetics2, KineticsData):
            if len(kinetics1.Tdata.value_si) != len(kinetics2.Tdata.value_si) or any(
                [T1 != T2 for T1, T2 in zip(kinetics1.Tdata.value_si, kinetics2.Tdata.value_si)]
            ):
                raise KineticsError("Cannot add these KineticsData objects due to " "their having different temperature points.")
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
                Ea=(kinetics1.Ea.value_si + kinetics2.Ea.value_si, "J/mol"),
                T0=(kinetics1.T0.value, kinetics1.T0.units),
            )
        elif isinstance(kinetics1, ArrheniusEP) and isinstance(kinetics2, ArrheniusEP):
            assert kinetics1.A.units == kinetics2.A.units
            kinetics = ArrheniusEP(
                A=(kinetics1.A.value * kinetics2.A.value, kinetics1.A.units),
                n=(kinetics1.n.value + kinetics2.n.value, kinetics1.n.units),
                alpha=kinetics1.alpha + kinetics2.alpha,
                E0=(kinetics1.E0.value_si + kinetics2.E0.value_si, "J/mol"),
            )
        elif isinstance(kinetics1, Arrhenius) and isinstance(kinetics2, ArrheniusEP):
            assert kinetics1.A.units == kinetics2.A.units
            assert kinetics1.T0.units == "K"
            assert kinetics1.T0.value == 1.0
            kinetics = ArrheniusEP(
                A=(kinetics1.A.value * kinetics2.A.value, kinetics1.A.units),
                n=(kinetics1.n.value + kinetics2.n.value, kinetics1.n.units),
                alpha=kinetics2.alpha,
                E0=(kinetics1.Ea.value_si + kinetics2.E0.value_si, "J/mol"),
            )
        elif isinstance(kinetics1, ArrheniusEP) and isinstance(kinetics2, Arrhenius):
            assert kinetics1.A.units == kinetics2.A.units
            assert "K" == kinetics2.T0.units
            assert 1.0 == kinetics2.T0.value
            kinetics = ArrheniusEP(
                A=(kinetics1.A.value * kinetics2.A.value, kinetics1.A.units),
                n=(kinetics1.n.value + kinetics2.n.value, kinetics1.n.units),
                alpha=kinetics1.alpha,
                E0=(kinetics1.E0.value_si + kinetics2.Ea.value_si, "J/mol"),
            )
        else:
            raise KineticsError('Unable to multiply kinetics types "{0}" and ' '"{1}".'.format(kinetics1.__class__, kinetics2.__class__))

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

        if kinetics1.comment == "":
            kinetics.comment = kinetics2.comment
        elif kinetics2.comment == "":
            kinetics.comment = kinetics1.comment
        else:
            kinetics.comment = kinetics1.comment + " + " + kinetics2.comment
        return kinetics
