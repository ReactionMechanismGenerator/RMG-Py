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
A class for returning and estimating the transport properties of a species

"""
from __future__ import division

import logging
import os.path
from copy import deepcopy

import rmgpy.constants as constants
from rmgpy.data.base import Database, Entry, make_logic_node, DatabaseError
from rmgpy.molecule import Molecule, Group
from rmgpy.transport import TransportData


def save_entry(f, entry):
    """
    Write a Pythonic string representation of the given `entry` in the transport
    database to the file object `f`.
    """
    f.write('entry(\n')
    f.write('    index = {0:d},\n'.format(entry.index))
    f.write('    label = "{0}",\n'.format(entry.label))

    if isinstance(entry.item, Molecule):
        f.write('    molecule = \n')
        f.write('"""\n')
        f.write(entry.item.to_adjacency_list(remove_h=False))
        f.write('""",\n')
    elif isinstance(entry.item, Group):
        f.write('    group = \n')
        f.write('"""\n')
        f.write(entry.item.to_adjacency_list())
        f.write('""",\n')
    else:
        f.write('    group = "{0}",\n'.format(entry.item))

    if isinstance(entry.data, CriticalPointGroupContribution):
        f.write('    transportGroup = CriticalPointGroupContribution(\n')
        f.write('        Tc = {0!r},\n'.format(entry.data.Tc))
        f.write('        Pc = {0!r},\n'.format(entry.data.Pc))
        f.write('        Vc = {0!r},\n'.format(entry.data.Vc))
        f.write('        Tb = {0!r},\n'.format(entry.data.Tb))
        f.write('        structureIndex = {0!r},\n'.format(entry.data.structureIndex))
        f.write('    ),\n')
    elif entry.data is None:
        f.write('    transportGroup = None,\n')
    elif isinstance(entry.data, TransportData):
        f.write('    transport = TransportData(\n')
        f.write('        shapeIndex = {0!r},\n'.format(entry.data.shapeIndex))
        f.write('        epsilon = {0!r},\n'.format(entry.data.epsilon))
        f.write('        sigma = {0!r},\n'.format(entry.data.sigma))
        f.write('        dipoleMoment = {0!r},\n'.format(entry.data.dipoleMoment))
        f.write('        polarizability = {0!r},\n'.format(entry.data.polarizability))
        f.write('        rotrelaxcollnum = {0!r},\n'.format(entry.data.rotrelaxcollnum))
        f.write('    ),\n')
    else:
        raise DatabaseError("Not sure how to save {0!r}".format(entry.data))

    f.write(f'    shortDesc = """{entry.shortDesc.strip()}""",\n')
    f.write(f'    longDesc = \n"""\n{entry.longDesc.strip()}\n""",\n')

    f.write(')\n\n')


def generate_old_library_entry(data):
    """
    Return a list of values used to save entries to the old-style RMG
    transport database based on the transport object `data`.
    """
    raise NotImplementedError


def process_old_library_entry(data):
    """
    Process a list of parameters `data` as read from an old-style RMG
    transport database, returning the corresponding transport object.
    """
    raise NotImplementedError


class TransportLibrary(Database):
    """
    A class for working with a RMG transport library.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def load_entry(self,
                   index,
                   label,
                   molecule,
                   transport,
                   reference=None,
                   referenceType='',
                   shortDesc='',
                   longDesc='',
                   ):
        item = Molecule().from_adjacency_list(molecule)

        self.entries[label] = Entry(
            index=index,
            label=label,
            item=item,
            data=transport,
            reference=reference,
            referenceType=referenceType,
            shortDesc=shortDesc,
            longDesc=longDesc.strip(),
        )

    def save_entry(self, f, entry):
        """
        Write the given `entry` in the transport database to the file object `f`.
        """
        return save_entry(f, entry)

    def generate_old_library_entry(self, data):
        """
        Return a list of values used to save entries to the old-style RMG
        transport database based on the transport object `data`.
        """
        return generate_old_library_entry(data)

    def process_old_library_entry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        transport database, returning the corresponding transport object.
        """
        return process_old_library_entry(data)


class TransportGroups(Database):
    """
    A class for working with an RMG transport group additivity database.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def load_entry(self,
                   index,
                   label,
                   group,
                   transportGroup,
                   reference=None,
                   referenceType='',
                   shortDesc='',
                   longDesc='',
                   ):
        if (group[0:3].upper() == 'OR{' or
                group[0:4].upper() == 'AND{' or
                group[0:7].upper() == 'NOT OR{' or
                group[0:8].upper() == 'NOT AND{'):
            item = make_logic_node(group)
        else:
            item = Group().from_adjacency_list(group)
        self.entries[label] = Entry(
            index=index,
            label=label,
            item=item,
            data=transportGroup,
            reference=reference,
            referenceType=referenceType,
            shortDesc=shortDesc,
            longDesc=longDesc.strip(),
        )

    def save_entry(self, f, entry):
        """
        Write the given `entry` in the transport database to the file object `f`.
        """
        return save_entry(f, entry)

    def generate_old_library_entry(self, data):
        """
        Return a list of values used to save entries to the old-style RMG
        transport database based on the transport object `data`.
        """
        return generate_old_library_entry(data)

    def process_old_library_entry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        transport database, returning the corresponding transport object.
        """
        return process_old_library_entry(data)


class TransportDatabase(object):
    """
    A class for working with the RMG transport database.
    """

    def __init__(self):
        self.libraries = {}
        self.groups = {}
        self.library_order = []
        self.local_context = {
            'CriticalPointGroupContribution': CriticalPointGroupContribution,
            'TransportData': TransportData,
        }
        self.global_context = {}

    def __reduce__(self):
        """
        A helper function used when pickling a TransportDatabase object.
        """
        d = {
            'libraries': self.libraries,
            'groups': self.groups,
            'library_order': self.library_order,
        }
        return (TransportDatabase, (), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling a TransportDatabase object.
        """
        self.libraries = d['libraries']
        self.groups = d['groups']
        self.library_order = d['library_order']

    def load(self, path, libraries=None):
        """
        Load the transport database from the given `path` on disk, where `path`
        points to the top-level folder of the transport database.
        """
        self.load_libraries(os.path.join(path, 'libraries'), libraries)
        self.load_groups(os.path.join(path, 'groups'))

    def load_libraries(self, path, libraries=None):
        """
        Load the transport libraries from the given `path` on disk, where `path`
        points to the libraries folder of the transport database.
        """
        self.libraries = {}
        self.library_order = []
        for (root, dirs, files) in os.walk(os.path.join(path)):
            for f in files:
                name, ext = os.path.splitext(f)
                if ext.lower() == '.py' and (libraries is None or name in libraries):
                    logging.info('Loading transport library from {0} in {1}...'.format(f, root))
                    library = TransportLibrary()
                    library.load(os.path.join(root, f), self.local_context, self.global_context)
                    library.label = os.path.splitext(f)[0]
                    self.libraries[library.label] = library
                    self.library_order.append(library.label)
        if libraries is not None:
            self.library_order = libraries

    def load_groups(self, path):
        """
        Load the transport groups from the given `path` on disk, where `path`
        points to the groups folder of the transport database.
        """
        logging.info('Loading transport group database from {0}...'.format(path))
        self.groups = {'ring': TransportGroups(label='ring').load(os.path.join(path, 'ring.py'),
                                                                  self.local_context, self.global_context),
                       'nonring': TransportGroups(label='nonring').load(os.path.join(path, 'nonring.py'),
                                                                        self.local_context, self.global_context)}

    def save(self, path):
        """
        Save the transport database to the given `path` on disk, where `path`
        points to the top-level folder of the transport database.
        """
        path = os.path.abspath(path)
        if not os.path.exists(path): os.mkdir(path)
        self.save_libraries(os.path.join(path, 'libraries'))
        self.save_groups(os.path.join(path, 'groups'))

    def save_libraries(self, path):
        """
        Save the trasnport libraries to the given `path` on disk, where `path`
        points to the top-level folder of the transport libraries.
        """
        if not os.path.exists(path): os.mkdir(path)
        for library in self.libraries.values():
            library.save(os.path.join(path, '{0}.py'.format(library.label)))

    def save_groups(self, path):
        """
        Save the transport groups to the given `path` on disk, where `path`
        points to the top-level folder of the transport groups.
        """
        if not os.path.exists(path): os.mkdir(path)
        for group in self.groups.keys():
            self.groups[group].save(os.path.join(path, group + '.py'))

    def get_transport_properties(self, species):
        """
        Return the transport properties for a given :class:`Species`
        object `species`. This function first searches the loaded libraries
        in order, returning the first match found, before falling back to
        estimation via group additivity.
        """
        transport = (None, None, None)

        if species.contains_surface_site():
            return transport

        for label in self.library_order:
            transport = self.get_transport_properties_from_library(species, self.libraries[label])
            if transport is not None:
                transport[0].comment = label
                break
        else:
            try:
                # Transport not found in any loaded libraries, so estimate
                transport = self.get_transport_properties_via_group_estimates(species)
            except (KeyError, AssertionError):
                transport = self.get_transport_properties_via_lennard_jones_parameters(species)

        return transport

    def get_all_transport_properties(self, species):
        """
        Return all possible sets of transport parameters for a given
        :class:`Species` object `species`. The hits from the libraries (in order) come first, and then the group additivity
        estimate. This method is useful for a generic search job.
        """
        transport = []

        # Data from libraries comes first
        for label in self.library_order:
            data = self.get_transport_properties_from_library(species, self.libraries[label])
            if data:
                data[0].comment = label
                transport.append(data)
        # Last entry is always the estimate from group additivity
        try:
            transport.append(self.get_transport_properties_via_group_estimates(species))
        except (KeyError, AssertionError):
            transport.append(self.get_transport_properties_via_lennard_jones_parameters(species))

        # Return all of the resulting transport parameters
        return transport

    def get_transport_properties_from_library(self, species, library):
        """
        Return the set of transport properties corresponding to a given
        :class:`Species` object `species` from the specified transport
        `library`. If `library` is a string, the list of libraries is searched
        for a library with that name. If no match is found in that library,
        ``None`` is returned. If no corresponding library is found, a
        :class:`DatabaseError` is raised.
        """
        for entry in library.entries.values():
            if species.is_isomorphic(entry.item) and entry.data is not None:
                return deepcopy(entry.data), library, entry
        return None

    def get_transport_properties_via_group_estimates(self, species):
        """
        Return the set of transport parameters corresponding to a given
        :class:`Species` object `species` by estimation using the group
        additivity values. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        """

        # assume that the stablest resonance isomer has already been put as the first
        # and that we want the transport properties of this isomer
        molecule = species.molecule[0]
        molecule.clear_labeled_atoms()
        molecule.update_atomtypes()
        critical_point = self.estimate_critical_properties_via_group_additivity(molecule)
        Tc = critical_point.Tc
        Pc = critical_point.Pc
        Vc = critical_point.Vc
        Tb = critical_point.Tb
        if critical_point.linear != molecule.is_linear():
            logging.warning("Group-based structure index and is_linear() function disagree about "
                            "linearity of {mol!r}".format(mol=molecule))

        if len(molecule.atoms) == 1:
            shape_index = 0
        elif molecule.is_linear():
            shape_index = 1
        else:
            shape_index = 2

        # Acetone values from Joback thesis: Tc = 511.455  (based on experimental Tb)  Pc = 47.808    Vc = 209.000    Tb = 322.082
        # print "Tc={Tc:.2f} K, Pc={Pc:.4g} bar, Vc={Vc:.4g} cm3/mol, Tb={Tb:.4g} K, average of {isomers} isomers".format(Tc=Tc,Pc=Pc,Vc=Vc,Tb=Tb,isomers=counter)
        # print 'Estimated with Tc={Tc:.2f} K, Pc={Pc:.4g} bar (from Joback method)'.format(Tc=Tc,Pc=Pc)
        transport = TransportData(
            shapeIndex=shape_index,
            epsilon=(.77 * Tc * constants.R, 'J/mol'),
            sigma=(2.44 * (Tc / Pc) ** (1. / 3), 'angstroms'),
            dipoleMoment=(0, 'C*m'),
            polarizability=(0, 'angstroms^3'),
            rotrelaxcollnum=0,  # rotational relaxation collision number at 298 K
            comment='Epsilon & sigma estimated with Tc={Tc:.2f} K, Pc={Pc:.4g} bar (from Joback method)'.format(Tc=Tc, Pc=Pc),
        )
        return transport, None, None
        # Things calling this expect a tuple with the library and entry that it came from.

    def estimate_critical_properties_via_group_additivity(self, molecule):
        """
        Return the set of transport parameters corresponding to a given
        :class:`Molecule` object `molecule` by estimation using Joback's group
        additivity values. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        
        Radicals are saturated with H atoms and the parent molecule properties
        are returned.
        """
        # For transport estimation we need the atoms to already be sorted because we
        # iterate over them; if the order changes during the iteration then we
        # will probably not visit the right atoms, and so will get the transport wrong

        if molecule.is_radical():  # radical species
            # Make a copy of the structure so we don't change the original
            saturated_struct = molecule.copy(deep=True)

            # Saturate structure by replacing all radicals with bonds to
            # hydrogen atoms
            saturated_struct.saturate_radicals()

            # Get critical point contribution estimates for saturated form of structure
            critical_point = self.estimate_critical_properties_via_group_additivity(saturated_struct)
            if critical_point is None:
                raise ValueError('Critical point estimate of saturated {0} of molecule {1} '
                                 'is None!'.format(saturated_struct, molecule))

            # We have no radical corrections for critical point estimates, so this is commented out:
            # # For each radical site, get radical correction
            # # Only one radical site should be considered at a time; all others
            # # should be saturated with hydrogen atoms
            # for atom in added:
            #    # Remove the added hydrogen atoms and bond and restore the radical
            #    for H, bond in added[atom]:
            #        saturated_struct.remove_bond(bond)
            #        saturated_struct.remove_atom(H)
            #        atom.increment_radical()
            #    saturated_struct.update_connectivity_values()
            return critical_point

        # non-radical species
        num_atoms = 0
        group_data = CriticalPointGroupContribution(
            Tc=0,
            Pc=0,
            Vc=0,
            Tb=0,
            structureIndex=0,
        )

        # Generate estimate of critical point contribution data
        for atom in molecule.atoms:
            num_atoms += 1
            # Iterate over heavy (non-hydrogen) atoms
            if atom.is_non_hydrogen():
                try:
                    if molecule.is_vertex_in_cycle(atom):
                        self._add_critical_point_contribution(group_data, self.groups['ring'], molecule, {'*': atom})
                    else:
                        self._add_critical_point_contribution(group_data, self.groups['nonring'], molecule, {'*': atom})
                except KeyError:
                    raise

        Tb = 198.18 + group_data.Tb
        Vc = 17.5 + group_data.Vc
        Tc = Tb / (0.584 + 0.965 * group_data.Tc - (group_data.Tc * group_data.Tc))
        Pc = 1 / (0.113 + 0.0032 * num_atoms + group_data.Pc) ** 2
        is_linear = (group_data.structureIndex == 0)

        critical_point = CriticalPoint(
            Tc=Tc,
            Pc=Pc,
            Vc=Vc,
            Tb=Tb,
            linear=is_linear,
        )
        return critical_point

    def _add_critical_point_contribution(self, group_data, database, molecule, atom):
        """
        Determine the critical point contribution values for the atom `atom`
        in the structure `structure`, and add it to the existing criticalPointContribution
        `criticalPointContribution`.
        """

        node0 = database.descend_tree(molecule, atom, None)

        if node0 is None:
            raise KeyError('Node not found in database.')

        # It's possible (and allowed) that items in the tree may not be in the
        # library, in which case we need to fall up the tree until we find an
        # ancestor that has an entry in the library
        node = node0

        while node is not None and node.data is None:
            node = node.parent
        if node is None:
            raise KeyError('Node {!r} has no parent with data in the transport database.'.format(node0))
        data = node.data
        comment = node.label
        while isinstance(data, str) and data is not None:
            for entry in database.entries.values():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    break
        comment = '{0}({1})'.format(database.label, comment)

        group_data.Tc += data.Tc
        group_data.Pc += data.Pc
        group_data.Vc += data.Vc
        group_data.Tb += data.Tb
        group_data.structureIndex += data.structureIndex

        return group_data

    def get_transport_properties_via_lennard_jones_parameters(self, species):
        """
        Serves as last resort if every other method to estimate Transport Properties fails.
        
        Generate the Lennard-Jones parameters for the species.
        """
        count = sum([1 for atom in species.molecule[0].vertices if atom.is_non_hydrogen()])

        if count == 1:
            sigma = (3.758e-10, "m")
            epsilon = (148.6, "K")
        elif count == 2:
            sigma = (4.443e-10, "m")
            epsilon = (110.7, "K")
        elif count == 3:
            sigma = (5.118e-10, "m")
            epsilon = (237.1, "K")
        elif count == 4:
            sigma = (4.687e-10, "m")
            epsilon = (531.4, "K")
        elif count == 5:
            sigma = (5.784e-10, "m")
            epsilon = (341.1, "K")
        else:
            sigma = (5.949e-10, "m")
            epsilon = (399.3, "K")

        if len(species.molecule[0].atoms) == 1:
            shape_index = 0
        elif species.molecule[0].is_linear():
            shape_index = 1
        else:
            shape_index = 2

        transport = TransportData(
            shapeIndex=shape_index,
            epsilon=epsilon,
            sigma=sigma,
            dipoleMoment=(0, 'C*m'),
            polarizability=(0, 'angstroms^3'),
            rotrelaxcollnum=0,  # rotational relaxation collision number at 298 K
            comment='Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!'
        )

        return transport, None, None


class CriticalPoint(object):
    """
    The critical properties of the species (and structureIndex)
    """

    def __init__(self, Tc=None, Pc=None, Vc=None, Tb=None, linear=None):
        self.Tc = Tc
        self.Pc = Pc
        self.Vc = Vc
        self.Tb = Tb
        self.linear = linear

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        CriticalPoint object
        """
        string = 'CriticalPoint(Tc={0!r}, Pc={1!r}, Vc={2!r}, Tb={3!r}, linear={4!r}'.format(
            self.Tc, self.Pc, self.Vc, self.Tb, self.linear)
        string += ')'
        return string


class CriticalPointGroupContribution(object):
    """Joback group contribution to estimate critical properties"""

    def __init__(self, Tc=None, Pc=None, Vc=None, Tb=None, structureIndex=None):
        self.Tc = Tc
        self.Pc = Pc
        self.Vc = Vc
        self.Tb = Tb
        self.structureIndex = structureIndex  # 0 if linear, 1 if makes molecule nonlinear

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        CriticalPointGroupContribution object
        """
        string = 'CriticalPointGroupContribution(Tc={0!r}, Pc={1!r}, Vc={2!r}, Tb={3!r}, structureIndex={4!r}'.format(
            self.Tc, self.Pc, self.Vc, self.Tb, self.structureIndex)
        string += ')'
        return string
