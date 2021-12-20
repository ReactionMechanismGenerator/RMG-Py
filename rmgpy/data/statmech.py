#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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

import logging
import os.path

import numpy as np

import rmgpy.constants as constants
from rmgpy.data.base import Database, Entry, LogicOr, make_logic_node
from rmgpy.data.statmechfit import fit_statmech_to_heat_capacity
from rmgpy.molecule import Molecule, Group
from rmgpy.statmech import Conformer, HarmonicOscillator, LinearRotor, NonlinearRotor, HinderedRotor, \
                           IdealGasTranslation


################################################################################

def save_entry(f, entry):
    """
    Write a Pythonic string representation of the given `entry` in the thermo
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

    if isinstance(entry.data, GroupFrequencies):
        f.write('    statmech = GroupFrequencies(\n')
        f.write('        frequencies = [\n')
        for lower, upper, degeneracy in entry.data.frequencies:
            f.write('            ({0:g}, {1:g}, {2:d}),\n'.format(lower, upper, degeneracy))
        f.write('        ],\n')
        f.write('        symmetry = {0:d},\n'.format(entry.data.symmetry))
        f.write('    ),\n')
    else:
        f.write('    statmech = {0!r},\n'.format(entry.data))

    if entry.reference is not None:
        f.write('    reference = {0!r},\n'.format(entry.reference))
    if entry.reference_type != "":
        f.write('    referenceType = "{0}",\n'.format(entry.reference_type))
    f.write(f'    shortDesc = """{entry.short_desc.strip()}""",\n')
    f.write(f'    longDesc = \n"""\n{entry.long_desc.strip()}\n""",\n')

    f.write(')\n\n')


def generate_old_library_entry(data):
    """
    Return a list of values used to save entries to the old-style RMG
    thermo database based on the thermodynamics object `data`.
    """
    items = '{0:3d}'.format(data.symmetry)
    for lower, upper, degeneracy in data.frequencies:
        items += '     {0:3d} {1:9.1f} {2:9.1f}'.format(degeneracy, lower, upper)
    return items


def process_old_library_entry(data, format):
    """
    Process a list of parameters `data` as read from an old-style RMG
    thermo database, returning the corresponding thermodynamics object.
    """
    if format == 'library':
        pass
    elif format == 'groups':
        frequencies = []
        for degeneracy, lower, upper in zip(data[1::3], data[2::3], data[3::3]):
            frequencies.append((float(lower), float(upper), int(degeneracy)))
        return GroupFrequencies(frequencies=frequencies, symmetry=int(data[0]))
    else:
        raise ValueError('format parameter must be "library" or "groups"; got "{0}" instead'.format(format))


################################################################################

class StatmechDepository(Database):
    """
    A class for working with the RMG statistical mechanics (frequencies)
    depository.
    """

    def __init__(self, label='', name='', short_desc='', long_desc=''):
        Database.__init__(self, label=label, name=name, short_desc=short_desc, long_desc=long_desc)

    def load_entry(self,
                   index,
                   label,
                   molecule,
                   statmech,
                   reference=None,
                   referenceType='',
                   shortDesc='',
                   longDesc='',
                   ):
        """
        Method for parsing entries in database files.
        Note that these argument names are retained for backward compatibility.
        """
        item = Molecule().from_adjacency_list(molecule)

        self.entries[label] = Entry(
            index=index,
            label=label,
            item=item,
            data=statmech,
            reference=reference,
            reference_type=referenceType,
            short_desc=shortDesc,
            long_desc=longDesc.strip(),
        )

    def save_entry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return save_entry(f, entry)


################################################################################

class StatmechLibrary(Database):
    """
    A class for working with a RMG statistical mechanics (frequencies) library.
    """

    def __init__(self, label='', name='', short_desc='', long_desc=''):
        Database.__init__(self, label=label, name=name, short_desc=short_desc, long_desc=long_desc)

    def load_entry(self,
                   index,
                   label,
                   molecule,
                   statmech,
                   reference=None,
                   referenceType='',
                   shortDesc='',
                   longDesc='',
                   ):
        """
        Method for parsing entries in database files.
        Note that these argument names are retained for backward compatibility.
        """
        self.entries[label] = Entry(
            index=index,
            label=label,
            item=Molecule().from_adjacency_list(molecule),
            data=statmech,
            reference=reference,
            reference_type=referenceType,
            short_desc=shortDesc,
            long_desc=longDesc.strip(),
        )

    def save_entry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return save_entry(f, entry)

    def generate_old_library_entry(self, data):
        """
        Return a list of values used to save entries to the old-style RMG
        thermo database based on the thermodynamics object `data`.
        """
        return generate_old_library_entry(data)

    def process_old_library_entry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        thermo database, returning the corresponding thermodynamics object.
        """
        return process_old_library_entry(data, "library")


################################################################################

class StatmechGroups(Database):
    """
    A class for working with an RMG statistical mechanics (frequencies) group
    database.
    """

    def __init__(self, label='', name='', short_desc='', long_desc=''):
        Database.__init__(self, label=label, name=name, short_desc=short_desc, long_desc=long_desc)

    def load_entry(self,
                   index,
                   label,
                   group,
                   statmech,
                   reference=None,
                   referenceType='',
                   shortDesc='',
                   longDesc='',
                   ):
        """
        Method for parsing entries in database files.
        Note that these argument names are retained for backward compatibility.
        """
        if (group[0:3].upper() == 'OR{' or
                group[0:4].upper() == 'AND{' or
                group[0:7].upper() == 'NOT OR{' or
                group[0:8].upper() == 'NOT AND{'
        ):
            item = make_logic_node(group)
        else:
            item = Group().from_adjacency_list(group)
        self.entries[label] = Entry(
            index=index,
            label=label,
            item=item,
            data=statmech,
            reference=reference,
            reference_type=referenceType,
            short_desc=shortDesc,
            long_desc=longDesc.strip(),
        )

    def save_entry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return save_entry(f, entry)

    def generate_old_library_entry(self, data):
        """
        Return a list of values used to save entries to the old-style RMG
        thermo database based on the thermodynamics object `data`.
        """
        return generate_old_library_entry(data)

    def process_old_library_entry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        statmech database, returning the corresponding thermodynamics object.
        """
        return process_old_library_entry(data, "groups")

    def _count_matches_to_node(self, molecule, node):
        """
        Count the number of matches in the given :class:`Molecule` object
        `molecule` for the given `node` in the group database.
        """
        count = 0
        if isinstance(self.dictionary[node], Group):
            ismatch, mappings = molecule.find_subgraph_isomorphisms(self.dictionary[node])
            count = len(mappings) if ismatch else 0
        elif isinstance(self.dictionary[node], LogicOr):
            for child in self.dictionary[node].components:
                count += self._count_matches_to_node(molecule, child)
        return count

    def _get_node(self, molecule, atom):
        """
        For a given :class:`Molecule` object `molecule` with central `atom`,
        determine the most specific functional group that describes that atom
        center and has characteristic frequencies associated with it.
        """

        node0 = self.descend_tree(molecule, atom, None)

        if node0 is None:
            raise KeyError('Node not found in database.')

        # It's possible (and allowed) that items in the tree may not be in the
        # library, in which case we need to fall up the tree until we find an
        # ancestor that has an entry in the library
        node = node0
        while node.data is None and node.parent is not None:
            node = node.parent
        if node.data is None:
            logging.warning('Statmech node {0!r} and all its parents have data = None'.format(node0.label))
            return None
            raise KeyError('Statmech node {0!r} and all its parents have data = None'.format(node0.label))
        return node

    def get_frequency_groups(self, molecule):
        """
        Return the set of characteristic group frequencies corresponding to the
        specified `molecule`. This is done by searching the molecule for
        certain functional groups for which characteristic frequencies are
        known, and using those frequencies.
        """

        frequencies = []
        group_count = {}

        # This is an additional hardcoded functional group for C-H with C in a ring
        # It is hardcoded because the adjacency list format isn't well-equipped
        # to handle this sort of functional group
        ring_ch = Entry(
            label='ringCH',
            item=None,
            data=GroupFrequencies([(2750., 3150., 1), (900., 1100., 1)]),
        )

        # Generate estimate of thermodynamics
        for atom in molecule.atoms:
            # Iterate over heavy (non-hydrogen) atoms
            if atom.is_hydrogen() or atom.is_halogen(): continue
            if molecule.is_atom_in_cycle(atom) and not atom.is_bonded_to_halogen():
                # Atom is in cycle
                # Add each C-H bond to the ringCH group
                # This is hardcoding of functional groups!
                if atom.is_carbon():
                    for atom2 in atom.edges:
                        if atom2.is_hydrogen():
                            try:
                                group_count[ring_ch] += 1
                            except KeyError:
                                group_count[ring_ch] = 1
            else:
                # Atom is not in cycle, so find a group for it
                node = self._get_node(molecule, {'*': atom})
                if node is not None:
                    try:
                        group_count[node] += 1
                    except KeyError:
                        group_count[node] = 1

        return group_count

    def get_statmech_data(self, molecule, thermo_model):
        """
        Use the previously-loaded frequency database to generate a set of
        characteristic group frequencies corresponding to the speficied
        `molecule`. The provided thermo data in `thermo_model` is used to fit
        some frequencies and all hindered rotors to heat capacity data.
        """
        conformer = Conformer()

        # Compute spin multiplicity
        # For closed-shell molecule the spin multiplicity is 1
        # For monoradicals the spin multiplicity is 2
        # For higher-order radicals the highest allowed spin multiplicity is assumed
        conformer.spin_multiplicity = molecule.get_radical_count() + 1

        # No need to determine rotational and vibrational modes for single atoms
        if len(molecule.atoms) < 2:
            return conformer, None, None

        linear = molecule.is_linear()
        num_rotors = molecule.count_internal_rotors()
        num_vibrations = 3 * len(molecule.atoms) - (5 if linear else 6) - num_rotors

        # Get characteristic frequency groups and the associated frequencies
        group_count = self.get_frequency_groups(molecule)
        logging.debug('Found frequencies from groups {}'.format(group_count))
        frequencies = []
        for entry, count in group_count.items():
            if count != 0 and entry.data is not None:
                frequencies.extend(entry.data.generate_frequencies(count))

        # Check that we have the right number of degrees of freedom specified
        if len(frequencies) > num_vibrations:
            # We have too many vibrational modes
            difference = len(frequencies) - num_vibrations
            # First try to remove hindered rotor modes until the proper number of modes remain
            if num_rotors > difference:
                num_rotors -= difference
                num_vibrations = len(frequencies)
                logging.warning('For {0}, more characteristic frequencies were generated than '
                                'vibrational modes allowed. Removed {1:d} internal rotors to '
                                'compensate.'.format(molecule.to_smiles(), difference))
            # If that won't work, turn off functional groups until the problem is underspecified again
            else:
                groups_removed = 0
                freqs_removed = 0
                total_freqs_removed = 0
                freq_count = len(frequencies)
                while freq_count > num_vibrations:
                    min_entry = min((entry for entry in group_count if group_count[entry] > 0),
                                    key=lambda x: x.data.symmetry)
                    if group_count[min_entry] > 1:
                        group_count[min_entry] -= 1
                    else:
                        del group_count[min_entry]
                    groups_removed += 1
                    freqs_removed = len(min_entry.data.generate_frequencies())
                    total_freqs_removed += freqs_removed
                    freq_count -= freqs_removed
                # Log warning
                logging.warning('For {0}, more characteristic frequencies were generated than '
                                'vibrational modes allowed. Removed {1:d} groups ({2:d} frequencies) to '
                                'compensate.'.format(molecule.to_smiles(), groups_removed, total_freqs_removed))
                # Regenerate characteristic frequencies
                frequencies = []
                for entry, count in group_count.items():
                    if count != 0:
                        frequencies.extend(entry.data.generate_frequencies(count))

        # Subtract out contributions to heat capacity from the group frequencies
        Tlist = np.arange(300.0, 1501.0, 100.0, float)
        Cv = np.array([thermo_model.get_heat_capacity(T) / constants.R for T in Tlist], float)
        logging.debug('Fitting statmech with heat capacities {0}'.format(Cv))
        ho = HarmonicOscillator(frequencies=(frequencies, "cm^-1"))
        for i in range(Tlist.shape[0]):
            Cv[i] -= ho.get_heat_capacity(Tlist[i]) / constants.R
        logging.debug('After removing found frequencies, the heat capacities are {0}'.format(Cv))
        # Subtract out translational modes
        Cv -= 1.5
        # Subtract out external rotational modes
        Cv -= (1.5 if not linear else 1.0)
        # Subtract out PV term (Cp -> Cv)
        Cv -= 1.0
        logging.debug('After removing translation, rotation, and Cp->Cv, the heat capacities are {0}'.format(Cv))
        # Fit remaining frequencies and hindered rotors to the heat capacity data
        modes = fit_statmech_to_heat_capacity(Tlist, Cv, num_vibrations - len(frequencies), num_rotors, molecule)
        for mode in modes:
            if isinstance(mode, HarmonicOscillator):
                uncertainties = [0 for f in frequencies]  # probably shouldn't be zero
                frequencies.extend(mode.frequencies.value_si)
                uncertainties.extend(mode.frequencies.uncertainty)
                mode.frequencies.value_si = np.array(frequencies, float)
                mode.frequencies.uncertainty = np.array(uncertainties, float)
                break
        else:
            modes.insert(0, HarmonicOscillator(frequencies=(frequencies, "cm^-1")))

        conformer.modes = modes

        return conformer, None, None


################################################################################

class StatmechDatabase(object):
    """
    A class for working with the RMG statistical mechanics (frequencies) database.
    """

    def __init__(self):
        self.depository = {}
        self.libraries = {}
        self.groups = {}
        self.library_order = []
        self.local_context = {
            'HarmonicOscillator': HarmonicOscillator,
            'LinearRotor': LinearRotor,
            'NonlinearRotor': NonlinearRotor,
            'HinderedRotor': HinderedRotor,
            'IdealGasTranslation': IdealGasTranslation,
            'GroupFrequencies': GroupFrequencies,
            'Conformer' : Conformer,
        }
        self.global_context = {}

    def __reduce__(self):
        """
        A helper function used when pickling a StatmechDatabase object.
        """
        d = {
            'depository': self.depository,
            'libraries': self.libraries,
            'groups': self.groups,
            'library_order': self.library_order,
        }
        return StatmechDatabase, (), d

    def __setstate__(self, d):
        """
        A helper function used when unpickling a StatmechDatabase object.
        """
        self.depository = d['depository']
        self.libraries = d['libraries']
        self.groups = d['groups']
        self.library_order = d['library_order']

    def load(self, path, libraries=None, depository=True):
        """
        Load the statmech database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        if depository:
            self.load_depository(os.path.join(path, 'depository'))
        else:
            self.depository = {}
        self.load_libraries(os.path.join(path, 'libraries'), libraries)
        self.load_groups(os.path.join(path, 'groups'))

    def load_depository(self, path):
        """
        Load the statmech database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        self.depository = {'depository': StatmechDepository().load(os.path.join(path, 'depository.py'),
                                                                   self.local_context, self.global_context)}

    def load_libraries(self, path, libraries=None):
        """
        Load the statmech database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        self.libraries = {}
        self.library_order = []
        for (root, dirs, files) in os.walk(os.path.join(path)):
            for f in files:
                name, ext = os.path.splitext(f)
                if ext.lower() == '.py' and (libraries is None or name in libraries):
                    logging.info('Loading frequencies library from {0} in {1}...'.format(f, root))
                    library = StatmechLibrary()
                    library.load(os.path.join(root, f), self.local_context, self.global_context)
                    library.label = os.path.splitext(f)[0]
                    self.libraries[library.label] = library
                    self.library_order.append(library.label)
        if libraries is not None:
            self.library_order = libraries

    def load_groups(self, path):
        """
        Load the statmech database from the given `path` on disk, where `path`
        points to the top-level folder of the thermo database.
        """
        logging.info('Loading frequencies group database from {0}...'.format(path))
        self.groups = {'groups': StatmechGroups().load(os.path.join(path, 'groups.py'),
                                                       self.local_context, self.global_context)}

    def save(self, path):
        """
        Save the statmech database to the given `path` on disk, where `path`
        points to the top-level folder of the statmech database.
        """
        path = os.path.abspath(path)
        if not os.path.exists(path):
            os.mkdir(path)
        self.save_depository(os.path.join(path, 'depository'))
        self.save_libraries(os.path.join(path, 'libraries'))
        self.save_groups(os.path.join(path, 'groups'))

    def save_depository(self, path):
        """
        Save the statmech depository to the given `path` on disk, where `path`
        points to the top-level folder of the statmech depository.
        """
        if not os.path.exists(path):
            os.mkdir(path)
        for name, depository in self.depository.items():
            depository.save(os.path.join(path, name + '.py'))

    def save_libraries(self, path):
        """
        Save the statmech libraries to the given `path` on disk, where `path`
        points to the top-level folder of the statmech libraries.
        """
        if not os.path.exists(path):
            os.mkdir(path)
        for library in self.libraries.values():
            library.save(os.path.join(path, '{0}.py'.format(library.label)))

    def save_groups(self, path):
        """
        Save the statmech groups to the given `path` on disk, where `path`
        points to the top-level folder of the statmech groups.
        """
        if not os.path.exists(path):
            os.mkdir(path)
        for name, groups in self.groups.items():
            groups.save(os.path.join(path, name + '.py'))

    def load_old(self, path):
        """
        Load the old RMG thermo database from the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        # The old database does not have a depository, so create an empty one
        self.depository = {'depository': StatmechDepository(label='depository', name='Statmech Depository')}

        for (root, dirs, files) in os.walk(os.path.join(path, 'frequencies_libraries')):
            if (os.path.exists(os.path.join(root, 'Dictionary.txt')) and
                    os.path.exists(os.path.join(root, 'Library.txt'))):
                library = StatmechLibrary(label=os.path.basename(root), name=os.path.basename(root))
                library.load_old(
                    dictstr=os.path.join(root, 'Dictionary.txt'),
                    treestr='',
                    libstr=os.path.join(root, 'Library.txt'),
                    num_parameters=-1,
                    num_labels=1,
                    pattern=False,
                )
                library.label = os.path.basename(root)
                self.libraries[library.label] = library

        self.groups['groups'] = StatmechGroups(label='group', name='Functional Group Values').load_old(
            dictstr=os.path.join(path, 'frequencies_groups', 'Dictionary.txt'),
            treestr=os.path.join(path, 'frequencies_groups', 'Tree.txt'),
            libstr=os.path.join(path, 'frequencies_groups', 'Library.txt'),
            num_parameters=-1,
            num_labels=1,
            pattern=True,
        )

    def save_old(self, path):
        """
        Save the old RMG thermo database to the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """

        # Depository not used in old database, so it is not saved

        libraries_path = os.path.join(path, 'frequencies_libraries')
        for library in self.libraries.values():
            if not os.path.exists(libraries_path):
                os.mkdir(libraries_path)
            library_path = os.path.join(libraries_path, library.label)
            if not os.path.exists(library_path):
                os.mkdir(library_path)
            library.save_old(
                dictstr=os.path.join(library_path, 'Dictionary.txt'),
                treestr='',
                libstr=os.path.join(library_path, 'Library.txt'),
            )

        groups_path = os.path.join(path, 'frequencies_groups')
        if not os.path.exists(groups_path):
            os.mkdir(groups_path)
        self.groups.save_old(
            dictstr=os.path.join(groups_path, 'Dictionary.txt'),
            treestr=os.path.join(groups_path, 'Tree.txt'),
            libstr=os.path.join(groups_path, 'Library.txt'),
        )

    def get_statmech_data(self, molecule, thermo_model=None):
        """
        Return the thermodynamic parameters for a given :class:`Molecule`
        object `molecule`. This function first searches the loaded libraries
        in order, returning the first match found, before falling back to
        estimation via group additivity.
        """
        logging.debug('Retrieving stat mech data for {}.'.format(molecule.to_smiles()))
        statmech_model = None
        # Check the libraries in order first; return the first successful match
        for label in self.library_order:
            statmech_model = self.get_statmech_data_from_library(molecule, self.libraries[label])
            if statmech_model:
                break
        else:
            # Thermo not found in any loaded libraries, so estimate
            statmech_model = self.get_statmech_data_from_groups(molecule, thermo_model)
        return statmech_model[0]

    def get_statmech_data_from_depository(self, molecule):
        """
        Return statmech data for the given :class:`Molecule` object `molecule`
        by searching the entries in the depository.
        Returns a list of tuples  (statmechData, depository, entry).
        """
        items = []
        for name, depository in self.depository.items():
            for label, entry in depository.entries.items():
                if molecule.is_isomorphic(entry.item):
                    items.append((entry.data, self.depository[name], entry))
        return items

    def get_statmech_data_from_library(self, molecule, library):
        """
        Return statmech data for the given :class:`Molecule` object `molecule`
        by searching the entries in the specified :class:`StatmechLibrary` object
        `library`. Returns ``None`` if no data was found.
        """
        for label, entry in library.entries.items():
            if molecule.is_isomorphic(entry.item):
                return entry.data, library, entry
        return None

    def get_statmech_data_from_groups(self, molecule, thermo_model):
        """
        Return statmech data for the given :class:`Molecule` object `molecule`
        by estimating using characteristic group frequencies and fitting the
        remaining internal modes to heat capacity data from the given thermo
        model `thermo_model`. This always returns valid degrees of freedom data.
        """
        return self.groups['groups'].get_statmech_data(molecule, thermo_model)


################################################################################

class GroupFrequencies(object):
    """
    Represent a set of characteristic frequencies for a group in the frequency
    database. These frequencies are stored in the `frequencies` attribute, which
    is a ``list`` of ``tuples``, where each ``tuple`` defines a lower bound,
    upper bound, and degeneracy. Each group also has a `symmetry` correction.
    """

    def __init__(self, frequencies=None, symmetry=1):
        self.frequencies = frequencies or []
        self.symmetry = symmetry

    def generate_frequencies(self, count=1):
        """
        Generate a set of frequencies. For each characteristic frequency group,
        the number of frequencies returned is degeneracy * count, and these are
        distributed linearly between the lower and upper bounds.
        """
        frequencies = []
        for lower, upper, degeneracy in self.frequencies:
            number = degeneracy * count
            if number == 1:
                frequencies.append((lower + upper) / 2.0)
            else:
                frequencies.extend(list(np.linspace(lower, upper, number, endpoint=True)))
        return frequencies
