#!/usr/bin/env python3

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

"""
import logging
import math
import os.path
from copy import deepcopy

import rmgpy.constants as constants
from rmgpy.data.base import Database, Entry, make_logic_node, DatabaseError
from rmgpy.molecule import Molecule, Group, ATOMTYPES
from rmgpy.species import Species


################################################################################

def save_entry(f, entry):
    """
    Write a Pythonic string representation of the given `entry` in the solvation
    database to the file object `f`.
    """
    f.write('entry(\n')
    f.write('    index = {0:d},\n'.format(entry.index))
    f.write('    label = "{0}",\n'.format(entry.label))

    if isinstance(entry.item, Species):
        if Molecule(smiles=entry.item.molecule[0].to_smiles()).is_isomorphic(entry.item.molecule[0]):
            # The SMILES representation accurately describes the molecule, so we can save it that way.
            f.write('    molecule = "{0}",\n'.format(entry.item.molecule[0].to_smiles()))
        else:
            f.write('    molecule = \n')
            f.write('"""\n')
            f.write(entry.item.to_adjacency_list(remove_h=False))
            f.write('""",\n')
    elif isinstance(entry.item, Group):
        f.write('    group = \n')
        f.write('"""\n')
        f.write(entry.item.to_adjacency_list())
        f.write('""",\n')
    elif entry.item is not None:
        f.write('    group = "{0}",\n'.format(entry.item))

    if isinstance(entry.data, SoluteData):
        f.write('    solute = SoluteData(\n')
        f.write('        S = {0!r},\n'.format(entry.data.S))
        f.write('        B = {0!r},\n'.format(entry.data.B))
        f.write('        E = {0!r},\n'.format(entry.data.E))
        f.write('        L = {0!r},\n'.format(entry.data.L))
        f.write('        A = {0!r},\n'.format(entry.data.A))
        if entry.data.V is not None: f.write('        V = {0!r},\n'.format(entry.data.V))
        f.write('    ),\n')
    elif isinstance(entry.data, SolventData):
        f.write('    solvent = SolventData(\n')
        f.write('        s_g = {0!r},\n'.format(entry.data.s_g))
        f.write('        b_g = {0!r},\n'.format(entry.data.b_g))
        f.write('        e_g = {0!r},\n'.format(entry.data.e_g))
        f.write('        l_g = {0!r},\n'.format(entry.data.l_g))
        f.write('        a_g = {0!r},\n'.format(entry.data.a_g))
        f.write('        c_g = {0!r},\n'.format(entry.data.c_g))
        f.write('        s_h = {0!r},\n'.format(entry.data.s_h))
        f.write('        b_h = {0!r},\n'.format(entry.data.b_h))
        f.write('        e_h = {0!r},\n'.format(entry.data.e_h))
        f.write('        l_h = {0!r},\n'.format(entry.data.l_h))
        f.write('        a_h = {0!r},\n'.format(entry.data.a_h))
        f.write('        c_h = {0!r},\n'.format(entry.data.c_h))
        f.write('        A = {0!r},\n'.format(entry.data.A))
        f.write('        B = {0!r},\n'.format(entry.data.B))
        f.write('        C = {0!r},\n'.format(entry.data.C))
        f.write('        D = {0!r},\n'.format(entry.data.D))
        f.write('        E = {0!r},\n'.format(entry.data.E))
        f.write('        alpha = {0!r},\n'.format(entry.data.alpha))
        f.write('        beta = {0!r},\n'.format(entry.data.beta))
        f.write('        eps = {0!r},\n'.format(entry.data.eps))
        f.write('    ),\n')
    elif entry.data is None:
        f.write('    solute = None,\n')
    else:
        raise DatabaseError("Not sure how to save {0!r}".format(entry.data))

    f.write(f'    shortDesc = """{entry.short_desc.strip()}""",\n')
    f.write(f'    longDesc = \n"""\n{entry.long_desc.strip()}\n""",\n')

    f.write(')\n\n')


def generate_old_library_entry(data):
    """
    Return a list of values used to save entries to the old-style RMG
    thermo database based on the thermodynamics object `data`.
    """
    raise NotImplementedError()


def process_old_library_entry(data):
    """
    Process a list of parameters `data` as read from an old-style RMG
    thermo database, returning the corresponding thermodynamics object.
    """
    raise NotImplementedError()


class SolventData(object):
    """
    Stores Abraham/Mintz parameters for characterizing a solvent.
    """

    def __init__(self, s_h=None, b_h=None, e_h=None, l_h=None, a_h=None,
                 c_h=None, s_g=None, b_g=None, e_g=None, l_g=None, a_g=None, c_g=None, A=None, B=None,
                 C=None, D=None, E=None, alpha=None, beta=None, eps=None):
        self.s_h = s_h
        self.b_h = b_h
        self.e_h = e_h
        self.l_h = l_h
        self.a_h = a_h
        self.c_h = c_h
        self.s_g = s_g
        self.b_g = b_g
        self.e_g = e_g
        self.l_g = l_g
        self.a_g = a_g
        self.c_g = c_g
        # These are parameters for calculating viscosity
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E
        # These are SOLUTE parameters used for intrinsic rate correction in H-abstraction rxns
        self.alpha = alpha
        self.beta = beta
        # This is the dielectric constant
        self.eps = eps

    def get_h_abs_correction(self):
        """
        If solvation is on, this will give the log10 of the ratio of the intrinsic rate
        constants log10(k_sol/k_gas) for H-abstraction rxns
        """
        return -8.3 * self.alpha * self.beta

    def get_solvent_viscosity(self, T):
        """
        Returns the viscosity in Pa s, according to correlation in Perry's Handbook
        and coefficients in DIPPR
        """
        return math.exp(self.A + (self.B / T) + (self.C * math.log(T)) + (self.D * (T ** self.E)))


class SolvationCorrection(object):
    """
    Stores corrections for enthalpy, entropy, and Gibbs free energy when a species is solvated.
    Enthalpy and Gibbs free energy is in J/mol; entropy is in J/mol/K
    """

    def __init__(self, enthalpy=None, gibbs=None, entropy=None):
        self.enthalpy = enthalpy
        self.entropy = entropy
        self.gibbs = gibbs


class SoluteData(object):
    """
    Stores Abraham parameters to characterize a solute
    """
    # Set class variable with McGowan volumes
    mcgowan_volumes = {
        1: 8.71, 2: 6.75,
        6: 16.35, 7: 14.39, 8: 12.43, 9: 10.47, 10: 8.51,
        14: 26.83, 15: 24.87, 16: 22.91, 17: 20.95, 18: 18.99,
        35: 26.21,
    }

    def __init__(self, S=None, B=None, E=None, L=None, A=None, V=None, comment=""):
        self.S = S
        self.B = B
        self.E = E
        self.L = L
        self.A = A
        self.V = V
        self.comment = comment

    def __repr__(self):
        return "SoluteData(S={0},B={1},E={2},L={3},A={4},comment={5!r})".format(
            self.S, self.B, self.E, self.L, self.A, self.comment)

    def get_stokes_diffusivity(self, T, solvent_viscosity):
        """
        Get diffusivity of solute using the Stokes-Einstein sphere relation. 
        Radius is found from the McGowan volume.
        solvent_viscosity should be given in  kg/s/m which equals Pa.s
        (water is about 9e-4 Pa.s at 25C, propanol is 2e-3 Pa.s)
        Returns D in m2/s
        """
        radius = math.pow((75 * self.V / constants.pi / constants.Na),
                          (1.0 / 3.0)) / 100  # in meters, V is in MgGowan volume in cm3/mol/100
        D = constants.kB * T / 6 / constants.pi / solvent_viscosity / radius  # m2/s
        return D  # m2/s

    def set_mcgowan_volume(self, species):
        """
        Find and store the McGowan's Volume
        Returned volumes are in cm^3/mol/100 (see note below)
        See Table 2 in Abraham & McGowan, Chromatographia Vol. 23, No. 4, p. 243. April 1987
        doi: 10.1007/BF02311772
        Also see Table 1 in Zhao et al., J. Chem. Inf. Comput. Sci. Vol. 43, p.1848. 2003
        doi: 10.1021/ci0341114
        
        "V is scaled to have similar values to the other
        descriptors by division by 100 and has units of (cm3mol−1/100)."
        the contibutions in this function are in cm3/mol, and the division by 100 is done at the very end.
        """
        molecule = species.molecule[0]  # any will do, use the first.
        Vtot = 0.0

        for atom in molecule.atoms:
            try:
                Vtot += self.mcgowan_volumes[atom.element.number]
            except KeyError:
                raise DatabaseError('McGowan volume not available for element {}'.format(atom.element.nubmer))

            # divide contribution in half since all bonds would be counted twice this way
            Vtot -= len(molecule.get_bonds(atom)) * 6.56 / 2

        self.V = Vtot / 100  # division by 100 to get units correct.


################################################################################


################################################################################

class SolventLibrary(Database):
    """
    A class for working with a RMG solvent library.
    """

    def __init__(self, label='', name='', short_desc='', long_desc=''):
        Database.__init__(self, label=label, name=name, short_desc=short_desc, long_desc=long_desc)

    def load_entry(self,
                   index,
                   label,
                   solvent,
                   molecule=None,
                   reference=None,
                   referenceType='',
                   shortDesc='',
                   longDesc='',
                   ):
        """
        Method for parsing entries in database files.
        Note that these argument names are retained for backward compatibility.
        """
        if molecule is not None:
            if not isinstance(molecule, list):
                molecule = [molecule]
            spc_list = []
            for mol in molecule:
                spc0 = Species(label=label)
                spc0.set_structure(mol)
                spc_list.append(spc0)
        else:
            spc_list = None

        self.entries[label] = Entry(
            index=index,
            label=label,
            item=spc_list,
            data=solvent,
            reference=reference,
            reference_type=referenceType,
            short_desc=shortDesc,
            long_desc=longDesc.strip(),
        )

    def load(self, path):
        """
        Load the solvent library from the given path
        """
        Database.load(self, path, local_context={'SolventData': SolventData}, global_context={})

    def save_entry(self, f, entry):
        """
        Write the given `entry` in the solute database to the file object `f`.
        """
        return save_entry(f, entry)

    def get_solvent_data(self, label):
        """
        Get a solvent's data from its name
        """
        return self.entries[label].data

    def get_solvent_structure(self, label):
        """
        Get a solvent's molecular structure as SMILES or adjacency list from its name
        """
        return self.entries[label].item


class SoluteLibrary(Database):
    """
    A class for working with a RMG solute library. Not currently used.
    """

    def __init__(self, label='', name='', short_desc='', long_desc=''):
        Database.__init__(self, label=label, name=name, short_desc=short_desc, long_desc=long_desc)

    def load_entry(self,
                   index,
                   label,
                   molecule,
                   solute,
                   reference=None,
                   referenceType='',
                   shortDesc='',
                   longDesc='',
                   ):
        """
        Method for parsing entries in database files.
        Note that these argument names are retained for backward compatibility.
        """
        try:
            spc = Species().from_smiles(molecule)
        except:
            logging.debug("Solute '{0}' does not have a valid SMILES '{1}'".format(label, molecule))
            try:
                spc = Species().from_adjacency_list(molecule)
            except:
                logging.error("Can't understand '{0}' in solute library '{1}'".format(molecule, self.name))
                raise

        self.entries[label] = Entry(
            index=index,
            label=label,
            item=spc,
            data=solute,
            reference=reference,
            reference_type=referenceType,
            short_desc=shortDesc,
            long_desc=longDesc.strip(),
        )

    def load(self, path):
        """
        Load the solute library from the given path
        """
        Database.load(self, path, local_context={'SoluteData': SoluteData}, global_context={})

    def save_entry(self, f, entry):
        """
        Write the given `entry` in the solute database to the file object `f`.
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
        return process_old_library_entry(data)


################################################################################

class SoluteGroups(Database):
    """
    A class for working with an RMG solute group additivity database.
    """

    def __init__(self, label='', name='', short_desc='', long_desc=''):
        Database.__init__(self, label=label, name=name, short_desc=short_desc, long_desc=long_desc)

    def load_entry(self,
                   index,
                   label,
                   group,
                   solute,
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
                group[0:8].upper() == 'NOT AND{'):
            item = make_logic_node(group)
        else:
            item = Group().from_adjacency_list(group)
        self.entries[label] = Entry(
            index=index,
            label=label,
            item=item,
            data=solute,
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
        return process_old_library_entry(data)


################################################################################

class SolvationDatabase(object):
    """
    A class for working with the RMG solvation database.
    """

    def __init__(self):
        self.libraries = {}
        self.libraries['solvent'] = SolventLibrary()
        self.libraries['solute'] = SoluteLibrary()
        self.groups = {}
        self.local_context = {
            'SoluteData': SoluteData,
            'SolventData': SolventData
        }
        self.global_context = {}

    def __reduce__(self):
        """
        A helper function used when pickling a SolvationDatabase object.
        """
        d = {
            'libraries': self.libraries,
            'groups': self.groups,
        }
        return (SolvationDatabase, (), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling a SolvationDatabase object.
        """
        self.libraries = d['libraries']
        self.groups = d['groups']

    def load(self, path, libraries=None, depository=True):
        """
        Load the solvation database from the given `path` on disk, where `path`
        points to the top-level folder of the solvation database.
        
        Load the solvent and solute libraries, then the solute groups.
        """

        self.libraries['solvent'].load(os.path.join(path, 'libraries', 'solvent.py'))
        self.libraries['solute'].load(os.path.join(path, 'libraries', 'solute.py'))

        self.load_groups(os.path.join(path, 'groups'))

    def get_solvent_data(self, solvent_name):
        try:
            solvent_data = self.libraries['solvent'].get_solvent_data(solvent_name)
        except:
            raise DatabaseError('Solvent {0!r} not found in database'.format(solvent_name))
        return solvent_data

    def get_solvent_structure(self, solvent_name):
        try:
            solvent_structure = self.libraries['solvent'].get_solvent_structure(solvent_name)
        except:
            raise DatabaseError('Solvent {0!r} not found in database'.format(solvent_name))
        return solvent_structure

    def load_groups(self, path):
        """
        Load the solute database from the given `path` on disk, where `path`
        points to the top-level folder of the solute database.
        
        Three sets of groups for additivity, atom-centered ('abraham'), non atom-centered 
        ('nonacentered'), and radical corrections ('radical')
        """
        logging.info('Loading Platts additivity group database from {0}...'.format(path))
        self.groups = {
            'abraham': SoluteGroups(label='abraham').load(os.path.join(path, 'abraham.py'),
                                                          self.local_context, self.global_context),
            'nonacentered': SoluteGroups(label='nonacentered').load(os.path.join(path, 'nonacentered.py'),
                                                                    self.local_context, self.global_context),
            'radical': SoluteGroups(label='radical').load(os.path.join(path, 'radical.py'),
                                                          self.local_context, self.global_context)
        }

    def save(self, path):
        """
        Save the solvation database to the given `path` on disk, where `path`
        points to the top-level folder of the solvation database.
        """
        path = os.path.abspath(path)
        if not os.path.exists(path):
            os.mkdir(path)
        self.save_libraries(os.path.join(path, 'libraries'))
        self.save_groups(os.path.join(path, 'groups'))

    def save_libraries(self, path):
        """
        Save the solute libraries to the given `path` on disk, where `path`
        points to the top-level folder of the solute libraries.
        """
        if not os.path.exists(path):
            os.mkdir(path)
        for library in self.libraries.keys():
            self.libraries[library].save(os.path.join(path, library + '.py'))

    def save_groups(self, path):
        """
        Save the solute groups to the given `path` on disk, where `path`
        points to the top-level folder of the solute groups.
        """
        if not os.path.exists(path):
            os.mkdir(path)
        for group in self.groups.keys():
            self.groups[group].save(os.path.join(path, group + '.py'))

    def load_old(self, path):
        """
        Load the old RMG solute database from the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """

        for (root, dirs, files) in os.walk(os.path.join(path, 'thermo_libraries')):
            if (os.path.exists(os.path.join(root, 'Dictionary.txt')) and
                    os.path.exists(os.path.join(root, 'Library.txt'))):
                library = SoluteLibrary(label=os.path.basename(root), name=os.path.basename(root))
                library.load_old(
                    dictstr=os.path.join(root, 'Dictionary.txt'),
                    treestr='',
                    libstr=os.path.join(root, 'Library.txt'),
                    num_parameters=5,
                    num_labels=1,
                    pattern=False,
                )
                library.label = os.path.basename(root)
                self.libraries[library.label] = library

        self.groups = {}
        self.groups['abraham'] = SoluteGroups(
            label='abraham',
            name='Platts Group Additivity Values for Abraham Solute Descriptors'
        ).load_old(
            dictstr=os.path.join(path, 'thermo_groups', 'Abraham_Dictionary.txt'),
            treestr=os.path.join(path, 'thermo_groups', 'Abraham_Tree.txt'),
            libstr=os.path.join(path, 'thermo_groups', 'Abraham_Library.txt'),
            num_parameters=5,
            num_labels=1,
            pattern=True,
        )

    def save_old(self, path):
        """
        Save the old RMG Abraham database to the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        # Depository not used in old database, so it is not saved

        libraries_path = os.path.join(path, 'thermo_libraries')
        if not os.path.exists(libraries_path):
            os.mkdir(libraries_path)
        for library in self.libraries.values():
            library_path = os.path.join(libraries_path, library.label)
            if not os.path.exists(library_path):
                os.mkdir(library_path)
            library.save_old(
                dictstr=os.path.join(library_path, 'Dictionary.txt'),
                treestr='',
                libstr=os.path.join(library_path, 'Library.txt'),
            )

        groups_path = os.path.join(path, 'thermo_groups')
        if not os.path.exists(groups_path):
            os.mkdir(groups_path)
        self.groups['abraham'].save_old(
            dictstr=os.path.join(groups_path, 'Abraham_Dictionary.txt'),
            treestr=os.path.join(groups_path, 'Abraham_Tree.txt'),
            libstr=os.path.join(groups_path, 'Abraham_Library.txt'),
        )

    def get_solute_data(self, species):
        """
        Return the solute descriptors for a given :class:`Species`
        object `species`. This function first searches the loaded libraries
        in order, returning the first match found, before falling back to
        estimation via Platts group additivity.
        """
        solute_data = None

        # Check the library first
        solute_data = self.get_solute_data_from_library(species, self.libraries['solute'])
        if solute_data is not None:
            assert len(solute_data) == 3, "solute_data should be a tuple (solute_data, library, entry)"
            solute_data[0].comment += "Data from solute library"
            solute_data = solute_data[0]
        else:
            # Solute not found in any loaded libraries, so estimate
            solute_data = self.get_solute_data_from_groups(species)
            # No Platts group additivity for V, so set using atom sizes
            solute_data.set_mcgowan_volume(species)
        # Return the resulting solute parameters S, B, E, L, A
        return solute_data

    def get_all_solute_data(self, species):
        """
        Return all possible sets of Abraham solute descriptors for a given
        :class:`Species` object `species`. The hits from the library come
        first, then the group additivity  estimate. This method is useful
        for a generic search job. Right now, there should either be 1 or 
        2 sets of descriptors, depending on whether or not we have a 
        library entry.
        """
        solute_data_list = []

        # Data from solute library
        data = self.get_solute_data_from_library(species, self.libraries['solute'])
        if data is not None:
            assert len(data) == 3, "solute_data should be a tuple (solute_data, library, entry)"
            data[0].comment += "Data from solute library"
            solute_data_list.append(data)
        # Estimate from group additivity
        # Make it a tuple
        data = (self.get_solute_data_from_groups(species), None, None)
        solute_data_list.append(data)
        return solute_data_list

    def get_solute_data_from_library(self, species, library):
        """
        Return the set of Abraham solute descriptors corresponding to a given
        :class:`Species` object `species` from the specified solute
        `library`. If `library` is a string, the list of libraries is searched
        for a library with that name. If no match is found in that library,
        ``None`` is returned. If no corresponding library is found, a
        :class:`DatabaseError` is raised.
        """
        for label, entry in library.entries.items():
            if species.is_isomorphic(entry.item) and entry.data is not None:
                return deepcopy(entry.data), library, entry
        return None

    def get_solute_data_from_groups(self, species):
        """
        Return the set of Abraham solute parameters corresponding to a given
        :class:`Species` object `species` by estimation using the Platts group
        additivity method. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        
        It averages (linearly) over the desciptors for each Molecule (resonance isomer)
        in the Species.
        """
        solute_data = SoluteData(0.0, 0.0, 0.0, 0.0, 0.0)
        count = 0
        comments = []
        for molecule in species.molecule:
            molecule.clear_labeled_atoms()
            molecule.update_atomtypes()
            sdata = self.estimate_solute_via_group_additivity(molecule)
            solute_data.S += sdata.S
            solute_data.B += sdata.B
            solute_data.E += sdata.E
            solute_data.L += sdata.L
            solute_data.A += sdata.A
            count += 1
            comments.append(sdata.comment)

        solute_data.S /= count
        solute_data.B /= count
        solute_data.E /= count
        solute_data.L /= count
        solute_data.A /= count

        # Print groups that are used for debugging purposes
        solute_data.comment = "Average of {0}".format(" and ".join(comments))

        return solute_data

    def transform_lone_pairs(self, molecule):
        """
        Changes lone pairs in a molecule to two radicals for purposes of finding
        solute data via group additivity. Transformed for each atom based on valency.
        """
        saturated_struct = molecule.copy(deep=True)
        added_to_pairs = {}

        for atom in saturated_struct.atoms:
            added_to_pairs[atom] = 0
            if atom.lone_pairs > 0:
                charge = atom.charge  # Record this so we can conserve it when checking
                bonds = saturated_struct.get_bonds(atom)
                sum_bond_orders = 0
                for key, bond in bonds.items():
                    sum_bond_orders += bond.order  # We should always have 2 'B' bonds (but what about Cbf?)
                if ATOMTYPES['Val4'] in atom.atomtype.generic:  # Carbon, Silicon
                    while atom.radical_electrons + charge + sum_bond_orders < 4:
                        atom.decrement_lone_pairs()
                        atom.increment_radical()
                        atom.increment_radical()
                        added_to_pairs[atom] += 1
                if ATOMTYPES['Val5'] in atom.atomtype.generic:  # Nitrogen
                    while atom.radical_electrons + charge + sum_bond_orders < 3:
                        atom.decrement_lone_pairs()
                        atom.increment_radical()
                        atom.increment_radical()
                        added_to_pairs[atom] += 1
                if ATOMTYPES['Val6'] in atom.atomtype.generic:  # Oxygen, sulfur
                    while atom.radical_electrons + charge + sum_bond_orders < 2:
                        atom.decrement_lone_pairs()
                        atom.increment_radical()
                        atom.increment_radical()
                        added_to_pairs[atom] += 1
                if ATOMTYPES['Val7'] in atom.atomtype.generic:  # Chlorine
                    while atom.radical_electrons + charge + sum_bond_orders < 1:
                        atom.decrement_lone_pairs()
                        atom.increment_radical()
                        atom.increment_radical()
                        added_to_pairs[atom] += 1

        saturated_struct.update()
        saturated_struct.update_lone_pairs()

        return saturated_struct, added_to_pairs

    def remove_h_bonding(self, saturated_struct, added_to_radicals, added_to_pairs, solute_data):

        # Remove hydrogen bonds and restore the radical
        for atom in added_to_radicals:
            for H, bond in added_to_radicals[atom]:
                saturated_struct.remove_bond(bond)
                saturated_struct.remove_atom(H)
                atom.increment_radical()

        # Change transformed lone pairs back
        for atom in added_to_pairs:
            if added_to_pairs[atom] > 0:
                for pair in range(1, added_to_pairs[atom]):
                    saturated_struct.decrement_radical()
                    saturated_struct.decrement_radical()
                    saturated_struct.increment_lone_pairs()

        # Update Abraham 'A' H-bonding parameter for unsaturated struct
        for atom in saturated_struct.atoms:
            # Iterate over heavy (non-hydrogen) atoms
            if atom.is_non_hydrogen() and atom.radical_electrons > 0:
                for electron in range(atom.radical_electrons):
                    # Get solute data for radical group
                    try:
                        self._add_group_solute_data(solute_data, self.groups['radical'], saturated_struct, {'*': atom})
                    except KeyError:
                        pass

        return solute_data

    def estimate_solute_via_group_additivity(self, molecule):
        """
        Return the set of Abraham solute parameters corresponding to a given
        :class:`Molecule` object `molecule` by estimation using the Platts' group
        additivity method. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        """
        # For thermo estimation we need the atoms to already be sorted because we
        # iterate over them; if the order changes during the iteration then we
        # will probably not visit the right atoms, and so will get the thermo wrong
        molecule.sort_atoms()

        # Create the SoluteData object with the intercepts from the Platts groups
        solute_data = SoluteData(
            S=0.277,
            B=0.071,
            E=0.248,
            L=0.13,
            A=0.003
        )

        added_to_radicals = {}  # Dictionary of key = atom, value = dictionary of {H atom: bond}
        added_to_pairs = {}  # Dictionary of key = atom, value = # lone pairs changed
        saturated_struct = molecule.copy(deep=True)

        # Convert lone pairs to radicals, then saturate with H.

        # Change lone pairs to radicals based on valency
        if (sum([atom.lone_pairs for atom in saturated_struct.atoms]) > 0 and  # molecule contains lone pairs
                not any([atom.atomtype.label == 'C2tc' for atom in saturated_struct.atoms])):  # and is not [C-]#[O+]
            saturated_struct, added_to_pairs = self.transform_lone_pairs(saturated_struct)

        # Now saturate radicals with H
        if sum([atom.radical_electrons for atom in saturated_struct.atoms]) > 0:  # radical species
            added_to_radicals = saturated_struct.saturate_radicals()

        # Saturated structure should now have no unpaired electrons, and only "expected" lone pairs
        # based on the valency
        for atom in saturated_struct.atoms:
            # Iterate over heavy (non-hydrogen) atoms
            if atom.is_non_hydrogen():
                # Get initial solute data from main group database. Every atom must
                # be found in the main abraham database
                try:
                    self._add_group_solute_data(solute_data, self.groups['abraham'], saturated_struct, {'*': atom})
                except KeyError:
                    logging.error("Couldn't find in main abraham database:")
                    logging.error(saturated_struct)
                    logging.error(saturated_struct.to_adjacency_list())
                    raise
                # Get solute data for non-atom centered groups (being found in this group
                # database is optional)    
                try:
                    self._add_group_solute_data(solute_data, self.groups['nonacentered'], saturated_struct, {'*': atom})
                except KeyError:
                    pass

        solute_data = self.remove_h_bonding(saturated_struct, added_to_radicals, added_to_pairs, solute_data)

        return solute_data

    def _add_group_solute_data(self, solute_data, database, molecule, atom):
        """
        Determine the Platts group additivity solute data for the atom `atom`
        in the structure `structure`, and add it to the existing solute data
        `solute_data`.
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
            raise KeyError('Node has no parent with data in database.')
        data = node.data
        comment = node.label
        while isinstance(data, str) and data is not None:
            for entry in database.entries.values():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    break
        comment = '{0}({1})'.format(database.label, comment)

        # This code prints the hierarchy of the found node; useful for debugging
        # result = ''
        # while node is not None:
        #   result = ' -> ' + node + result
        #   node = database.tree.parent[node]
        # print result[4:]

        # Add solute data for each atom to the overall solute data for the molecule.
        solute_data.S += data.S
        solute_data.B += data.B
        solute_data.E += data.E
        solute_data.L += data.L
        solute_data.A += data.A
        solute_data.comment += comment + "+"

        return solute_data

    def calc_h(self, solute_data, solvent_data):
        """
        Returns the enthalpy of solvation, at 298K, in J/mol
        """
        # Use Mintz parameters for solvents. Multiply by 1000 to go from kJ->J to maintain consistency
        delH = 1000 * ((solute_data.S * solvent_data.s_h) +
                       (solute_data.B * solvent_data.b_h) +
                       (solute_data.E * solvent_data.e_h) +
                       (solute_data.L * solvent_data.l_h) +
                       (solute_data.A * solvent_data.a_h) + solvent_data.c_h)
        return delH

    def calc_g(self, solute_data, solvent_data):
        """
        Returns the Gibbs free energy of solvation, at 298K, in J/mol
        """
        # Use Abraham parameters for solvents to get log K
        logK = ((solute_data.S * solvent_data.s_g) +
                (solute_data.B * solvent_data.b_g) +
                (solute_data.E * solvent_data.e_g) +
                (solute_data.L * solvent_data.l_g) +
                (solute_data.A * solvent_data.a_g) + solvent_data.c_g)
        # Convert to delG with units of J/mol
        delG = -8.314 * 298 * 2.303 * logK
        return delG

    def calc_s(self, delG, delH):
        """
        Returns the entropy of solvation, at 298K, in J/mol/K
        """
        delS = (delH - delG) / 298
        return delS

    def get_solvation_correction(self, solute_data, solvent_data):
        """ 
        Given a solute_data and solvent_data object, calculates the enthalpy, entropy,
        and Gibbs free energy of solvation at 298 K. Returns a SolvationCorrection
        object
        """
        correction = SolvationCorrection(0.0, 0.0, 0.0)
        correction.enthalpy = self.calc_h(solute_data, solvent_data)
        correction.gibbs = self.calc_g(solute_data, solvent_data)
        correction.entropy = self.calc_s(correction.gibbs, correction.enthalpy)
        return correction

    def check_solvent_in_initial_species(self, rmg, solvent_structure):
        """
        Given the instance of RMG class and the solvent_structure, it checks whether the solvent is listed as one
        of the initial species.
        If the SMILES / adjacency list for all the solvents exist in the solvent library, it uses the solvent's
        molecular structure to determine whether the species is the solvent or not.
        If the solvent library does not have SMILES / adjacency list, then it uses the solvent's string name
        to determine whether the species is the solvent or not
        """
        for spec in rmg.initial_species:
            if solvent_structure is not None:
                spec.is_solvent = spec.is_isomorphic(solvent_structure)
            else:
                spec.is_solvent = rmg.solvent == spec.label
        if not any([spec.is_solvent for spec in rmg.initial_species]):
            if solvent_structure is not None:
                logging.info('One of the initial species must be the solvent')
                raise ValueError('One of the initial species must be the solvent')
            else:
                logging.info('One of the initial species must be the solvent with the same string name')
                raise ValueError('One of the initial species must be the solvent with the same string name')
