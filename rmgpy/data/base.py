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
Contains classes and functions for working with the various RMG databases. In
particular, this module is devoted to functionality that is common across all
components of the RMG database.
"""

import codecs
import logging
import os
import re
from collections import OrderedDict

from rmgpy.data.reference import Reference, Article, Book, Thesis
from rmgpy.exceptions import DatabaseError, InvalidAdjacencyListError
from rmgpy.kinetics.uncertainties import RateUncertainty
from rmgpy.kinetics.arrhenius import ArrheniusChargeTransfer, ArrheniusChargeTransferBM
from rmgpy.molecule import Molecule, Group


################################################################################

class Entry(object):
    """
    A class for representing individual records in an RMG database. Each entry
    in the database associates a chemical item (generally a species, functional
    group, or reaction) with a piece of data corresponding to that item. A
    significant amount of metadata can also be stored with each entry.

    The attributes are:

    ====================== ========================================================
    Attribute              Description
    ====================== ========================================================
    `index`                A unique nonnegative integer index for the entry
    `label`                A unique string identifier for the entry (or '' if not used)
    `item`                 The item that this entry represents
    `parent`               The parent of the entry in the hierarchy (or ``None`` if not used)
    `children`             A list of the children of the entry in the hierarchy (or ``None`` if not used)
    `data`                 The data to associate with the item
    `data_count`           The number of data used to fit the group values in the group additivity method
    `reference`            A :class:`Reference` object containing bibliographic reference information to the source of the data
    `reference_type`       The way the data was determined: ``'theoretical'``, ``'experimental'``, or ``'review'``
    `short_desc`           A brief (one-line) description of the data
    `long_desc`            A long, verbose description of the data
    `rank`                 An integer indicating the degree of confidence in the entry data, or ``None`` if not used
    `nodal_distance`       A float representing the distance of a given entry from it's parent entry
     --                    For surface species thermo calculations:
    `metal`                Which metal the thermo calculation was done on (``None`` if not used)
    `facet`                Which facet the thermo calculation was done on (``None`` if not used)
    `site`                 Which surface site the molecule prefers (``None`` if not used)
    `binding_energies`     The surface binding energies for C,H,O, and N
    `surface_site_density` The surface site density
    ====================== ========================================================
    """

    def __init__(self,
                 index=-1,
                 label='',
                 item=None,
                 parent=None,
                 children=None,
                 data=None,
                 data_count=None,
                 reference=None,
                 reference_type='',
                 short_desc='',
                 long_desc='',
                 rank=None,
                 nodal_distance=None,
                 metal=None,
                 facet=None,
                 site=None,
                 binding_energies=None,
                 surface_site_density=None,
                 ):
        self.index = index
        self.label = label
        self.item = item
        self.parent = parent
        self.children = children or []
        self.data = data
        self.data_count = data_count
        self.reference = reference
        self.reference_type = reference_type
        self.short_desc = short_desc
        self.long_desc = long_desc
        self.rank = rank
        self.nodal_distance = nodal_distance
        self.metal = metal
        self.facet = facet
        self.site = site
        self.binding_energies = binding_energies
        self.surface_site_density = surface_site_density

    def __str__(self):
        return self.label

    def __repr__(self):
        return '<Entry index={0:d} label="{1}">'.format(self.index, self.label)

    def get_all_descendants(self):
        """
        retrieve all the descendants of entry
        """
        new_nodes = [self]
        tot_nodes = []
        temp_nodes = []
        while new_nodes:
            for entry in new_nodes:
                temp_nodes.extend(entry.children)
            tot_nodes.extend(new_nodes)
            new_nodes = temp_nodes
            temp_nodes = []

        tot_nodes.remove(self)
        return tot_nodes


################################################################################

class Database(object):
    """
    An RMG-style database, consisting of a dictionary of entries (associating
    items with data), and an optional tree for assigning a hierarchy to the
    entries. The use of the tree enables the database to be easily extensible
    as more parameters are available.

    In constructing the tree, it is important to develop a hierarchy such that
    siblings are mutually exclusive, to ensure that there is a unique path of
    descent down a tree for each structure. If non-mutually exclusive siblings
    are encountered, a warning is raised and the parent of the siblings is
    returned.

    There is no requirement that the children of a node span the range of
    more specific permutations of the parent. As the database gets more complex,
    attempting to maintain complete sets of children for each parent in each
    database rapidly becomes untenable, and is against the spirit of
    extensibility behind the database development.

    You must derive from this class and implement the :meth:`load_entry`,
    :meth:`save_entry`, :meth:`process_old_library_entry`, and
    :meth:`generate_old_library_entry` methods in order to load and save from the
    new and old database formats.
    """

    local_context = {}
    local_context['Reference'] = Reference
    local_context['Article'] = Article
    local_context['Book'] = Book
    local_context['Thesis'] = Thesis

    def __init__(self,
                 entries=None,
                 top=None,
                 label='',
                 name='',
                 solvent=None,
                 short_desc='',
                 long_desc='',
                 metal=None,
                 site=None,
                 facet=None,
                 ):
        self.entries = OrderedDict(entries or {})
        self.top = top or []
        self.label = label
        self.name = name
        self.solvent = solvent
        self.short_desc = short_desc
        self.long_desc = long_desc
        self.metal = metal
        self.site = site
        self.facet = facet

    def load(self, path, local_context=None, global_context=None):
        """
        Load an RMG-style database from the file at location `path` on disk.
        The parameters `local_context` and `global_context` are used to
        provide specialized mapping of identifiers in the input file to
        corresponding functions to evaluate. This method will automatically add
        a few identifiers required by all data entries, so you don't need to
        provide these.
        """

        # Clear any previously-loaded data
        self.entries = OrderedDict()
        self.top = []

        # Set up global and local context
        if global_context is None: global_context = {}
        global_context['__builtins__'] = None
        global_context['True'] = True
        global_context['False'] = False
        if local_context is None: local_context = {}
        local_context['__builtins__'] = None
        local_context['entry'] = self.load_entry
        local_context['tree'] = self._load_tree
        local_context['name'] = self.name
        local_context['solvent'] = self.solvent
        local_context['shortDesc'] = self.short_desc
        local_context['longDesc'] = self.long_desc
        local_context['RateUncertainty'] = RateUncertainty
        local_context['ArrheniusChargeTransfer'] = ArrheniusChargeTransfer
        local_context['ArrheniusChargeTransferBM'] = ArrheniusChargeTransferBM
        local_context['metal'] = self.metal
        local_context['site'] = self.site
        local_context['facet'] = self.facet
        # add in anything from the Class level dictionary.
        for key, value in Database.local_context.items():
            local_context[key] = value

        # Process the file
        f = open(path, 'r')
        try:
            exec(f.read(), global_context, local_context)
        except Exception:
            logging.error('Error while reading database {0!r}.'.format(path))
            raise
        f.close()

        # Extract the database metadata
        self.name = local_context['name']
        self.solvent = local_context['solvent']
        self.short_desc = local_context['shortDesc']
        self.long_desc = local_context['longDesc'].strip()
        self.metal = local_context['metal']
        self.site = local_context['site']
        self.facet = local_context['facet']

        # Return the loaded database (to allow for Database().load() syntax)
        return self

    def get_entries_to_save(self):
        """
        Return a sorted list of the entries in this database that should be
        saved to the output file.

        Then renumber the entry indexes so that we never have any duplicate indexes.
        """
        entries = self.top[:]
        if len(self.top) > 0:
            # Save the entries in the same order as the tree (so that it saves
            # in the same order each time)
            for entry in self.top:
                entries.extend(self.descendants(entry))
            # It may be that a logical or is defined such that its children
            # are not in the tree; this ensures that they still get saved
            index = 0
            while index < len(entries):
                entry = entries[index]
                if isinstance(entry.item, LogicOr):
                    descendants = self.descendants(entry)
                    for child in entry.item.components:
                        if self.entries[child] not in descendants:
                            entries.append(self.entries[child])
                index += 1
        else:
            # Otherwise save the entries sorted by index, if defined
            entries = list(self.entries.values())
            entries.sort(key=lambda x: x.index)

        for index, entry in enumerate(entries):
            entry.index = index

        return entries

    def get_species(self, path, resonance=True):
        """
        Load the dictionary containing all of the species in a kinetics library or depository.
        """
        from rmgpy.species import Species
        species_dict = OrderedDict()
        with open(path, 'r') as f:
            adjlist = ''
            for line in f:
                if line.strip() == '' and adjlist.strip() != '':
                    # Finish this adjacency list
                    species = Species().from_adjacency_list(adjlist)
                    if resonance:
                        species.generate_resonance_structures()
                    label = species.label
                    if label in species_dict:
                        raise DatabaseError('Species label "{0}" used for multiple species in {1}.'.format(label,
                                                                                                           str(self)))
                    species_dict[label] = species
                    adjlist = ''
                else:
                    adjlist += line
            else:  # reached end of file
                if adjlist.strip() != '':
                    # Finish this adjacency list
                    species = Species().from_adjacency_list(adjlist)
                    if resonance:
                        species.generate_resonance_structures()
                    label = species.label
                    if label in species_dict:
                        raise DatabaseError('Species label "{0}" used for multiple species in {1}.'.format(label,
                                                                                                           str(self)))
                    species_dict[label] = species

        return species_dict

    def save_dictionary(self, path):
        """
        Extract species from all entries associated with a kinetics library or depository and save them 
        to the path given.
        """
        try:
            os.makedirs(os.path.dirname(path))
        except OSError:
            pass
        # Extract species from all the entries
        species_dict = {}
        entries = self.entries.values()
        for entry in entries:
            for reactant in entry.item.reactants:
                if reactant.label not in species_dict:
                    species_dict[reactant.label] = reactant

            for product in entry.item.products:
                if product.label not in species_dict:
                    species_dict[product.label] = product

        with open(path, 'w') as f:
            for label in species_dict.keys():
                f.write(species_dict[label].molecule[0].to_adjacency_list(label=label, remove_h=False))
                f.write('\n')

    def save(self, path, reindex=True):
        """
        Save the current database to the file at location `path` on disk. 
        """
        try:
            os.makedirs(os.path.dirname(path))
        except OSError:
            pass
        
        if reindex:
            entries = self.get_entries_to_save()
        else: 
            entries = self.entries.values()

        f = codecs.open(path, 'w', 'utf-8')
        f.write('#!/usr/bin/env python\n')
        f.write('# encoding: utf-8\n\n')
        f.write('name = "{0}"\n'.format(self.name))
        f.write('shortDesc = "{0}"\n'.format(self.short_desc))
        f.write('longDesc = """\n')
        f.write(self.long_desc.strip() + '\n')
        f.write('"""\n')

        for entry in entries:
            self.save_entry(f, entry)

        # Write the tree
        if len(self.top) > 0:
            f.write('tree(\n')
            f.write('"""\n')
            f.write(self.generate_old_tree(self.top, 1))
            f.write('"""\n')
            f.write(')\n\n')

        f.close()

    def load_old(self, dictstr, treestr, libstr, num_parameters, num_labels=1, pattern=True):
        """
        Load a dictionary-tree-library based database. The database is stored
        in three files: `dictstr` is the path to the dictionary, `treestr` to
        the tree, and `libstr` to the library. The tree is optional, and should
        be set to '' if not desired.
        """

        # Load dictionary, library, and (optionally) tree
        try:
            self.load_old_dictionary(dictstr, pattern)
        except Exception:
            logging.error('Error while reading database {0!r}.'.format(os.path.dirname(dictstr)))
            raise

        try:
            if treestr != '': self.load_old_tree(treestr)
        except Exception:
            logging.error('Error while reading database {0!r}.'.format(os.path.dirname(treestr)))
            raise

        try:
            self.load_old_library(libstr, num_parameters, num_labels)
        except Exception:
            logging.error('Error while reading database {0!r}.'.format(os.path.dirname(libstr)))
            raise

        return self

    def load_old_dictionary(self, path, pattern):
        """
        Parse an old-style RMG database dictionary located at `path`. An RMG
        dictionary is a list of key-value pairs of a one-line string key and a
        multi-line string value. Each record is separated by at least one empty
        line. Returns a ``dict`` object with the values converted to
        :class:`Molecule` or :class:`Group` objects depending on the
        value of `pattern`.
        """

        # The dictionary being loaded
        self.entries = {}
        # The current record
        record = ''

        f_dict = None
        # Process the dictionary file
        try:
            f_dict = open(path, 'r')
            for line in f_dict:
                line = line.strip()
                # If at blank line, end of record has been found
                if len(line) == 0 and len(record) > 0:
                    # Label is first line of record
                    lines = record.splitlines()
                    label = lines[0]
                    # Add record to dictionary
                    self.entries[label] = Entry(label=label, item=record)
                    # Clear record in preparation for next iteration
                    record = ''
                # Otherwise append line to record (if not empty and not a comment line)
                else:
                    line = remove_comment_from_line(line).strip()
                    if len(line) > 0:
                        record += line + '\n'
            # process the last record! (after end of for loop)
            # Label is first line of record
            if record:
                label = record.splitlines()[0]
                # Add record to dictionary
                self.entries[label] = Entry(label=label, item=record)
        except DatabaseError as e:
            logging.exception(str(e))
            raise
        except IOError as e:
            logging.exception('Database dictionary file "' + e.filename + '" not found.')
            raise
        finally:
            if f_dict:
                f_dict.close()

        # Convert the records in the dictionary to Molecule, Group, or
        # logical objects
        try:
            for label in self.entries:
                record = self.entries[label].item
                lines = record.splitlines()
                # If record is a logical node, make it into one.
                if re.match("(?i)\s*(NOT\s)?\s*(OR|AND|UNION)\s*(\{.*\})", lines[1]):
                    self.entries[label].item = make_logic_node(' '.join(lines[1:]))
                # Otherwise convert adjacency list to molecule or pattern
                elif pattern:
                    self.entries[label].item = Group().from_adjacency_list(record)
                else:
                    self.entries[label].item = Molecule().from_adjacency_list(record, saturate_h=True)
        except InvalidAdjacencyListError:
            logging.error('Error while loading old-style dictionary "{0}"'.format(path))
            logging.error('Error occurred while parsing adjacency list "{0}"'.format(label))
            raise

    def _load_tree(self, tree):
        """
        Parse an group tree located at `tree`. An RMG tree is an n-ary
        tree representing the hierarchy of items in the dictionary.
        """

        if len(self.entries) == 0:
            raise DatabaseError("Load the dictionary before you load the tree.")

        # should match '  L3 : foo_bar '  and 'L3:foo_bar'
        parser = re.compile('^\s*L(?P<level>\d+)\s*:\s*(?P<label>\S+)')

        parents = [None]
        for line in tree.splitlines():
            line = remove_comment_from_line(line).strip()
            if len(line) > 0:
                # Extract level
                match = parser.match(line)
                if not match:
                    raise DatabaseError("Couldn't parse line '{0}'".format(line.strip()))
                level = int(match.group('level'))
                label = match.group('label')

                # Find immediate parent of the new node
                parent = None
                if len(parents) < level:
                    raise DatabaseError("Invalid level specified in line '{0}'".format(line.strip()))
                else:
                    while len(parents) > level:
                        parents.remove(parents[-1])
                    if len(parents) > 0:
                        parent = parents[level - 1]

                if parent is not None: parent = self.entries[parent]
                try:
                    entry = self.entries[label]
                except KeyError:
                    raise DatabaseError('Unable to find entry "{0}" from tree in dictionary.'.format(label))

                if isinstance(parent, str):
                    raise DatabaseError('Unable to find parent entry "{0}" of entry "{1}" in tree.'.format(parent,
                                                                                                           label))

                # Update the parent and children of the nodes accordingly
                if parent is not None:
                    entry.parent = parent
                    parent.children.append(entry)
                else:
                    entry.parent = None
                    self.top.append(entry)

                # Save the level of the tree into the entry
                entry.level = level

                # Add node to list of parents for subsequent iteration
                parents.append(label)

    def load_old_tree(self, path):
        """
        Parse an old-style RMG database tree located at `path`. An RMG
        tree is an n-ary tree representing the hierarchy of items in the
        dictionary.
        """
        tree = []
        try:
            ftree = open(path, 'r')
            tree = ftree.read()

        except IOError:
            logging.exception('Database tree file "' + e.filename + '" not found.')
        finally:
            ftree.close()

        self._load_tree(tree)

    def load_old_library(self, path, num_parameters, num_labels=1):
        """
        Parse an RMG database library located at `path`.
        """

        if len(self.entries) == 0:
            raise DatabaseError("Load the dictionary before you load the library.")

        entries = self.parse_old_library(path, num_parameters, num_labels)

        # Load the parsed entries into the database, skipping duplicate entries
        skipped_count = 0
        for index, label, parameters, comment in entries:
            if label not in self.entries:
                raise DatabaseError('Entry {0!r} in library was not found in dictionary.'.format(label))
            if self.entries[label].index != -1:
                # The entry is a duplicate, so skip it
                logging.debug("There was already something labeled {0} in the {1} library. "
                              "Ignoring '{2}' ({3})".format(label, self.label, index, parameters))
                skipped_count += 1
            else:
                # The entry is not a duplicate
                self.entries[label].index = index
                self.entries[label].data = parameters
                self.entries[label].short_desc = comment
        if skipped_count > 0:
            logging.warning("Skipped {0:d} duplicate entries in {1} library.".format(skipped_count, self.label))

        # Make sure each entry with data has a nonnegative index
        entries2 = list(self.entries.values())
        entries2.sort(key=lambda entry: entry.index)
        index = entries2[-1].index + 1
        if index < 1: index = 1
        for index0, label, parameters, comment in entries:
            if self.entries[label].index < 0:
                self.entries[label].index = index
                index += 1

    def parse_old_library(self, path, num_parameters, num_labels=1):
        """
        Parse an RMG database library located at `path`, returning the loaded
        entries (rather than storing them in the database). This method does
        not discard duplicate entries.
        """

        entries = []

        flib = None
        try:
            flib = codecs.open(path, 'r', 'utf-8', errors='replace')
            for line in flib:
                line = remove_comment_from_line(line).strip()
                if len(line) > 0:

                    info = line.split()

                    # Skip if the number of items on the line is invalid
                    if len(info) < 2:
                        continue

                    # Determine if the first item is an index
                    # This index is optional in the old library format
                    index = -1
                    offset = 0
                    try:
                        index = int(float(info[0]))
                        offset = 1
                    except ValueError:
                        pass
                    # Extract label(s)
                    label = self._hash_labels(info[offset:offset + num_labels])
                    offset += num_labels
                    # Extract numeric parameter(s) or label of node with data to use
                    if num_parameters < 0:
                        parameters = self.process_old_library_entry(info[offset:])
                        comment = ''
                    else:
                        try:
                            parameters = self.process_old_library_entry(info[offset:offset + num_parameters])
                            offset += num_parameters
                        except (IndexError, ValueError):
                            parameters = info[offset]
                            offset += 1
                        # Remaining part of string is comment
                        comment = ' '.join(info[offset:])
                        comment = comment.strip('"')

                    entries.append((index, label, parameters, comment))

        except DatabaseError as e:
            logging.exception(str(e))
            logging.exception("path = '{0}'".format(path))
            logging.exception("line = '{0}'".format(line))
            raise
        except IOError as e:
            logging.exception('Database library file "' + e.filename + '" not found.')
            raise
        finally:
            if flib: flib.close()

        return entries

    def save_old(self, dictstr, treestr, libstr):
        """
        Save the current database to a set of text files using the old-style
        syntax.
        """
        self.save_old_dictionary(dictstr)
        if treestr != '':
            self.save_old_tree(treestr)

        # RMG-Java does not require a frequencies_groups/Library.txt file to
        # operate, but errors are raised upon importing to Py if this file is
        # not found. This check prevents the placeholder from being discarded.
        if 'StatesGroups' not in self.__class__.__name__:
            self.save_old_library(libstr)

    def save_old_dictionary(self, path):
        """
        Save the current database dictionary to a text file using the old-style
        syntax.
        """

        entries = []
        entries_not_in_tree = []

        # If we have tree information, save the dictionary in the same order as
        # the tree (so that it saves in the same order each time)
        def get_logic_node_components(entry_or_item):
            """
            If we want to save an entry, but that is a logic node, we also want
            to save its components, recursively. This is a horribly complicated way
            to *not* save in the dictionary any things which are not accessed from
            (or needed to define things that are accessed from) the tree.
            """
            if isinstance(entry_or_item, Entry):
                entry = entry_or_item
                item = entry.item
                nodes = [entry]
            else:
                entry = None
                item = entry_or_item
                nodes = []
            if isinstance(item, LogicNode):
                for child in item.components:
                    if isinstance(child, LogicNode):
                        nodes.extend(get_logic_node_components(child))
                    else:
                        nodes.extend(get_logic_node_components(self.entries[child]))
                return nodes
            else:
                return [entry]

        if len(self.top) > 0:
            for entry in self.top:
                entries.extend(get_logic_node_components(entry))
                for descendant in self.descendants(entry):
                    for entry2 in get_logic_node_components(descendant):
                        if entry2 not in entries:
                            entries.append(entry2)

            # Don't forget entries that aren't in the tree
            for entry in self.entries.values():
                if entry not in entries:
                    entries_not_in_tree.append(entry)
            entries_not_in_tree.sort(key=lambda x: (x.index, x.label))
        # Otherwise save the dictionary in any order
        else:
            # Save the library in order by index
            entries = list(self.entries.values())
            entries.sort(key=lambda x: (x.index, x.label))

        def comment(s):
            """Return the string, with each line prefixed with '// '"""
            return '\n'.join('// ' + line if line else '' for line in s.split('\n'))

        try:
            f = open(path, 'w')
            f.write('////////////////////////////////////////////////////////////////////////////////\n')
            f.write('//\n')
            f.write('//  {0} dictionary\n'.format(self.name))
            f.write('//\n')
            f.write('////////////////////////////////////////////////////////////////////////////////\n')
            f.write('\n')
            for entry in entries:
                f.write(entry.label + '\n')
                if isinstance(entry.item, Molecule):
                    try:
                        f.write(entry.item.to_adjacency_list(remove_h=True, old_style=True) + '\n')
                    except InvalidAdjacencyListError:
                        f.write("// Couldn't save in old syntax adjacency list. Here it is in new syntax:\n")
                        f.write(comment(entry.item.to_adjacency_list(remove_h=False, old_style=False) + '\n'))
                elif isinstance(entry.item, Group):
                    f.write(entry.item.to_adjacency_list(old_style=True).replace('{2S,2T}', '2') + '\n')
                elif isinstance(entry.item, LogicOr):
                    f.write('{0}\n\n'.format(entry.item).replace('OR{', 'Union {'))
                elif entry.label[0:7] == 'Others-':
                    assert isinstance(entry.item, LogicNode)
                    f.write('{0}\n\n'.format(entry.item))
                else:
                    raise DatabaseError('Unexpected item with label {0} encountered in dictionary while '
                                        'attempting to save.'.format(entry.label))

            if entries_not_in_tree:
                f.write(comment("These entries do not appear in the tree:\n\n"))
            for entry in entries_not_in_tree:
                f.write(comment(entry.label + '\n'))
                if isinstance(entry.item, Molecule):
                    f.write(comment(entry.item.to_adjacency_list(remove_h=False) + '\n'))
                elif isinstance(entry.item, Group):
                    f.write(comment(entry.item.to_adjacency_list().replace('{2S,2T}', '2') + '\n'))
                elif isinstance(entry.item, LogicOr):
                    f.write(comment('{0}\n\n'.format(entry.item).replace('OR{', 'Union {')))
                elif entry.label[0:7] == 'Others-':
                    assert isinstance(entry.item, LogicNode)
                    f.write(comment('{0}\n\n'.format(entry.item)))
                else:
                    raise DatabaseError('Unexpected item with label {0} encountered in dictionary while '
                                        'attempting to save.'.format(entry.label))

            f.close()
        except IOError:
            logging.exception('Unable to save old-style dictionary to "{0}".'.format(os.path.abspath(path)))
            raise

    def generate_old_tree(self, entries, level):
        """
        Generate a multi-line string representation of the current tree using
        the old-style syntax.
        """
        string = ''
        for entry in entries:
            # Write current node
            string += '{0}L{1:d}: {2}\n'.format('    ' * (level - 1), level, entry.label)
            # Recursively descend children (depth-first)
            string += self.generate_old_tree(entry.children, level + 1)
        return string

    def save_old_tree(self, path):
        """
        Save the current database tree to a text file using the old-style
        syntax.
        """
        try:
            f = open(path, 'w')
            f.write('////////////////////////////////////////////////////////////////////////////////\n')
            f.write('//\n')
            f.write('//  {0} tree\n'.format(self.name))
            f.write('//\n')
            f.write('////////////////////////////////////////////////////////////////////////////////\n')
            f.write('\n')
            f.write(self.generate_old_tree(self.top, 1))
            f.close()
        except IOError:
            logging.exception('Unable to save old-style tree to "{0}".'.format(os.path.abspath(path)))
            raise

    def save_old_library(self, path):
        """
        Save the current database library to a text file using the old-style
        syntax.
        """
        try:
            # Save the library in order by index
            entries = list(self.entries.values())
            entries.sort(key=lambda x: x.index)

            f = codecs.open(path, 'w', 'utf-8')
            records = []
            for entry in entries:
                if entry.data is not None:
                    data = entry.data
                    if not isinstance(data, str):
                        data = self.generate_old_library_entry(data)
                    records.append((entry.index, [entry.label], data, entry.short_desc))

            records.sort()

            f.write('////////////////////////////////////////////////////////////////////////////////\n')
            f.write('//\n')
            f.write('//  {0} library\n'.format(self.name))
            f.write('//\n')
            f.write('////////////////////////////////////////////////////////////////////////////////\n')
            f.write('\n')
            for index, labels, data, comment in records:
                f.write('{:<6d} '.format(index))
                for label in labels:
                    f.write('{:<32s} '.format(label))
                if isinstance(data, str):
                    f.write('{:s} '.format(data))
                else:
                    f.write('{:s} '.format(' '.join(['{:<10g}'.format(d) for d in data])))
                f.write(u'    {:s}\n'.format(comment))
            f.close()
        except IOError:
            logging.exception('Unable to save old-style library to "{0}".'.format(os.path.abspath(path)))
            raise

    def _hash_labels(self, labels):
        """
        Convert a list of string `labels` to a list of single strings that
        represent permutations of the individual strings in the `labels` list::

            >>> hashLabels(['a','b'])
            ['a;b', 'b;a']
        """
        return ';'.join(labels)

    def ancestors(self, node):
        """
        Returns all the ancestors of a node, climbing up the tree to the top.
        """
        if isinstance(node, str):
            node = self.entries[node]
        ancestors = []
        parent = node.parent
        if parent is not None:
            ancestors = [parent]
            ancestors.extend(self.ancestors(parent))
        return ancestors

    def descendants(self, node):
        """
        Returns all the descendants of a node, climbing down the tree to the bottom.
        """
        if isinstance(node, str):
            node = self.entries[node]
        descendants = []
        for child in node.children:
            descendants.append(child)
            descendants.extend(self.descendants(child))
        return descendants

    def match_node_to_node(self, node, node_other):
        """ 
        Return `True` if `node` and `node_other` are identical.  Otherwise, return `False`.
        Both `node` and `node_other` must be Entry types with items containing Group or LogicNode types.
        """
        if isinstance(node.item, Group) and isinstance(node_other.item, Group):
            return (self.match_node_to_structure(node, node_other.item, atoms=node_other.item.get_all_labeled_atoms()) and
                    self.match_node_to_structure(node_other, node.item, atoms=node.item.get_all_labeled_atoms()))
        elif isinstance(node.item, LogicOr) and isinstance(node_other.item, LogicOr):
            return node.item.match_logic_or(node_other.item)
        else:
            # Assume nonmatching
            return False

    def match_node_to_child(self, parent_node, child_node):
        """ 
        Return `True` if `parent_node` is a parent of `child_node`.  Otherwise, return `False`.
        Both `parent_node` and `child_node` must be Entry types with items containing Group or LogicNode types.
        If `parent_node` and `child_node` are identical, the function will also return `False`.
        """

        if isinstance(parent_node.item, Group) and isinstance(child_node.item, Group):
            if (self.match_node_to_structure(parent_node, child_node.item,
                                             atoms=child_node.item.get_all_labeled_atoms(), strict=True) and
                    not self.match_node_to_structure(child_node, parent_node.item,
                                                     atoms=parent_node.item.get_all_labeled_atoms(), strict=True)):
                return True
            return False

        # If the parent_node is a Group and the child_node is a LogicOr there is nothing to check,
        # except that the parent is listed in the attributes. However, we do need to check that everything down this
        # family line is consistent, which is done in the databaseTest unitTest
        elif isinstance(parent_node.item, Group) and isinstance(child_node.item, LogicOr):
            ancestor_node = child_node.parent
            while ancestor_node:
                if ancestor_node is parent_node:
                    return True
                else:
                    ancestor_node = ancestor_node.parent
            else:
                return False

        elif isinstance(parent_node.item, LogicOr):
            return child_node.label in parent_node.item.components

    def match_node_to_structure(self, node, structure, atoms, strict=False):
        """
        Return :data:`True` if the `structure` centered at `atom` matches the
        structure at `node` in the dictionary. The structure at `node` should
        have atoms with the appropriate labels because they are set on loading
        and never change. However, the atoms in `structure` may not have the
        correct labels, hence the `atoms` parameter. The `atoms` parameter may
        include extra labels, and so we only require that every labeled atom in
        the functional group represented by `node` has an equivalent labeled
        atom in `structure`.
        
        Matching to structure is more strict than to node.  All labels in structure must 
        be found in node.  However the reverse is not true, unless `strict` is set to True.
        
        =================== ========================================================
        Attribute           Description
        =================== ========================================================
        `node`              Either an Entry or a key in the self.entries dictionary which has a Group or LogicNode as its Entry.item
        `structure`         A Group or a Molecule
        `atoms`             Dictionary of {label: atom} in the structure.  A possible dictionary is the one produced by structure.get_all_labeled_atoms()
        `strict`            If set to ``True``, ensures that all the node's atomLabels are matched by in the structure
        =================== ========================================================
        """
        if isinstance(node, str):
            node = self.entries[node]
        group = node.item
        if isinstance(group, LogicNode):
            return group.match_to_structure(self, structure, atoms, strict)
        else:
            # try to pair up labeled atoms
            centers = group.get_all_labeled_atoms()
            initial_map = {}
            for label in centers.keys():
                # Make sure the labels are in both group and structure.
                if label not in atoms:
                    logging.log(0, "Label {0} is in group {1} but not in structure".format(label, node))
                    if strict:
                        # structure must match all labeled atoms in node if strict is set to True
                        return False
                    continue  # with the next label - ring structures might not have all labeled atoms
                center = centers[label]
                atom = atoms[label]
                # Make sure labels actually point to atoms.
                if center is None or atom is None:
                    return False
                # Semantic check #1: atoms with same label are equivalent
                if isinstance(center, list) or isinstance(atom, list):
                    pass
                else:
                    if not atom.is_specific_case_of(center):
                        return False
                # Semantic check #2: labeled atoms that share bond in the group (node)
                # also share equivalent (or more specific) bond in the structure
                for atom2, atom1 in initial_map.items():
                    if group.has_bond(center, atom1) and structure.has_bond(atom, atom2):
                        bond1 = group.get_bond(center, atom1)  # bond1 is group
                        bond2 = structure.get_bond(atom, atom2)  # bond2 is structure
                        if not bond2.is_specific_case_of(bond1):
                            return False
                    elif group.has_bond(center, atom1):  # but structure doesn't
                        return False
                    elif structure.has_bond(atom, atom2):  # but group doesn't
                        # We don't mind that the structure has bond but the group doesn't
                        pass
                # Passed semantic checks, so add to maps of already-matched atoms
                if not (isinstance(center, list) or isinstance(atom, list)):
                    initial_map[atom] = center
            # Labeled atoms in the structure that are not in the group should
            # not be considered in the isomorphism check, so flag them temporarily
            # Without this we would hit a lot of nodes that are ambiguous
            flagged_atoms = [atom for label, atom in structure.get_all_labeled_atoms().items() if label not in centers]
            for atom in flagged_atoms:
                atom.ignore = True

            # use mapped (labeled) atoms to try to match subgraph
            result = structure.is_subgraph_isomorphic(group, initial_map)

            # Restore atoms flagged in previous step
            for atom in flagged_atoms:
                atom.ignore = False

            return result

    def descend_tree(self, structure, atoms, root=None, strict=False):
        """
        Descend the tree in search of the functional group node that best
        matches the local structure around `atoms` in `structure`.

        If root=None then uses the first matching top node.

        Returns None if there is no matching root.
        
        Set strict to ``True`` if all labels in final matched node must match that of the
        structure.  This is used in kinetics groups to find the correct reaction template, but
        not generally used in other GAVs due to species generally not being prelabeled.
        """

        if root is None:
            for root in self.top:
                if self.match_node_to_structure(root, structure, atoms, strict):
                    break  # We've found a matching root
            else:  # didn't break - matched no top nodes
                return None
        elif not self.match_node_to_structure(root, structure, atoms, strict):
            return None

        next_node = []
        for child in root.children:
            if self.match_node_to_structure(child, structure, atoms, strict):
                next_node.append(child)

        if len(next_node) == 1:
            return self.descend_tree(structure, atoms, next_node[0], strict)
        elif len(next_node) == 0:
            if len(root.children) > 0 and root.children[-1].label.startswith('Others-'):
                return root.children[-1]
            else:
                return root
        else:
            # logging.warning('For {0}, a node {1} with overlapping children {2} was encountered '
            #                 'in tree with top level nodes {3}. Assuming the first match is the '
            #                 'better one.'.format(structure, root, next, self.top))
            return self.descend_tree(structure, atoms, next_node[0], strict)

    def are_siblings(self, node, node_other):
        """
        Return `True` if `node` and `node_other` have the same parent node.  Otherwise, return `False`.
        Both `node` and `node_other` must be Entry types with items containing Group or LogicNode types.
        """
        if node.parent is node_other.parent:
            return True
        else:
            return False

    def remove_group(self, group_to_remove):
        """
        Removes a group that is in a tree from the database. In addition to deleting from self.entries,
        it must also update the parent/child relationships

        Returns the removed group
        """
        # Don't remove top nodes or LogicOrs as this will cause lots of problems
        if group_to_remove in self.top:
            raise ValueError("Cannot remove top node: {0} from {1} because it is a top node".format(group_to_remove, self))
        elif isinstance(group_to_remove.item, LogicOr):
            raise ValueError("Cannot remove top node: {0} from {1} because it is a LogicOr".format(group_to_remove, self))
        # Remove from entryToRemove from entries
        self.entries.pop(group_to_remove.label)

        # If there is a parent, then the group exists in a tree and we should edit relatives
        parent_r = group_to_remove.parent
        if parent_r is not None:
            # Remove from parent's children attribute
            parent_r.children.remove(group_to_remove)

            # change children's parent attribute to former grandparent
            for child in group_to_remove.children:
                child.parent = parent_r

            # extend parent's children attribute with new children
            parent_r.children.extend(group_to_remove.children)

            # A few additional changes needed if parent_r is a LogicOr node
            if isinstance(parent_r.item, LogicOr):
                parent_r.item.components.remove(group_to_remove.label)
                parent_r.item.components.extend([child.label for child in group_to_remove.children])

        return group_to_remove


class LogicNode(object):
    """
    A base class for AND and OR logic nodes.
    """

    symbol = "<TBD>"  # To be redefined by subclass

    def __init__(self, items, invert):
        self.components = []
        for item in items:
            if re.match("(?i)\s*(NOT\s)?\s*(OR|AND|UNION)\s*(\{.*\})", item):
                component = make_logic_node(item)
            else:
                component = item
            self.components.append(component)
        self.invert = bool(invert)

    def __str__(self):
        result = ''
        if self.invert: result += 'NOT '
        result += self.symbol
        result += "{{{0}}}".format(', '.join([str(c) for c in self.components]))
        return result


class LogicOr(LogicNode):
    """
    A logical OR node. Structure can match any component.

    Initialize with a list of component items and a boolean instruction to invert the answer.
    """

    symbol = "OR"

    def match_to_structure(self, database, structure, atoms, strict=False):
        """
        Does this node in the given database match the given structure with the labeled atoms?
        
        Setting `strict` to True makes enforces matching of atomLabels in the structure to every
        atomLabel in the node.
        """
        for node in self.components:
            if isinstance(node, LogicNode):
                match = node.match_to_structure(database, structure, atoms, strict)
            else:
                match = database.match_node_to_structure(node, structure, atoms, strict)
            if match:
                return True != self.invert
        return False != self.invert

    def match_logic_or(self, other):
        """
        Is other the same LogicOr group as self?
        """
        if len(self.components) != len(other.components):
            return False
        else:
            for node in self.components:
                if node not in other.components:
                    return False
        return True

    def get_possible_structures(self, entries):
        """
        Return a list of the possible structures below this node.
        """
        if self.invert:
            raise NotImplementedError("Finding possible structures of NOT OR nodes not implemented.")
        structures = []
        for item in self.components:
            struct = entries[item].item
            if isinstance(struct, LogicNode):
                structures.extend(struct.get_possible_structures(entries))
            else:
                structures.append(struct)
        for struct in structures:  # check this worked
            assert isinstance(struct, Group)
        return structures


class LogicAnd(LogicNode):
    """A logical AND node. Structure must match all components."""

    symbol = "AND"

    def match_to_structure(self, database, structure, atoms, strict=False):
        """
        Does this node in the given database match the given structure with the labeled atoms?
        
        Setting `strict` to True makes enforces matching of atomLabels in the structure to every
        atomLabel in the node.
        """
        for node in self.components:
            if isinstance(node, LogicNode):
                match = node.match_to_structure(database, structure, atoms, strict)
            else:
                match = database.match_node_to_structure(node, structure, atoms, strict)
            if not match:
                return False != self.invert
        return True != self.invert


def make_logic_node(string):
    """
    Creates and returns a node in the tree which is a logic node.

    String should be of the form:

    * OR{}
    * AND{}
    * NOT OR{}
    * NOT AND{}

    And the returned object will be of class LogicOr or LogicAnd
    """

    match = re.match(r"(?i)\s*(NOT\s)?\s*(OR|AND|UNION)\s*(\{.*\})", string)  # the (?i) makes it case-insensitive
    if not match:
        raise ValueError("Unexpected string for Logic Node: {0}".format(string))

    if match.group(1):
        invert = True
    else:
        invert = False

    logic = match.group(2)  # OR or AND (or Union)

    contents = match.group(3).strip()
    while contents.startswith('{'):
        if not contents.endswith('}'):
            raise ValueError("Unbalanced braces in Logic Node: {0}".format(string))
        contents = contents[1:-1]

    items = []
    chars = []
    brace_depth = 0
    for character in contents:
        if character == '{':
            brace_depth += 1
        if character == '}':
            brace_depth -= 1
        if character == ',' and brace_depth == 0:
            items.append(''.join(chars).lstrip().rstrip())
            chars = []
        else:
            chars.append(character)
    if chars:  # add last item
        items.append(''.join(chars).lstrip().rstrip())
    if brace_depth != 0:
        raise ValueError("Unbalanced braces in Logic Node: {0}".format(string))

    if logic.upper() in ['OR', 'UNION']:
        return LogicOr(items, invert)
    if logic == 'AND':
        return LogicAnd(items, invert)

    raise ValueError("Could not create Logic Node from {0}".format(string))


################################################################################

def remove_comment_from_line(line):
    """
    Remove a C++/Java style comment from a line of text. This refers
    particularly to comments that begin with a double-slash '//' and continue
    to the end of the line.
    """

    index = line.find('//')
    if index >= 0:
        line = line[0:index]
    return line


def split_line_and_comment(line):
    """
    Returns a tuple(line, comment) based on a '//' comment delimiter.
    
    Either `line` or `comment` may be ''.
    Does not strip whitespace, nor remove more than two slashes.
    """
    split = line.split('//', 1)
    if len(split) == 1:
        return split[0], ''
    else:
        return tuple(split)


def get_all_combinations(node_lists):
    """
    Generate a list of all possible combinations of items in the list of
    lists `node_lists`. Each combination takes one item from each list
    contained within `node_lists`. The order of items in the returned lists
    reflects the order of lists in `node_lists`. For example, if `node_lists` was
    [[A, B, C], [N], [X, Y]], the returned combinations would be
    [[A, N, X], [A, N, Y], [B, N, X], [B, N, Y], [C, N, X], [C, N, Y]].
    """

    items = [[]]
    for node_list in node_lists:
        items = [item + [node] for node in node_list for item in items]

    return items


################################################################################

class ForbiddenStructures(Database):
    """
    A database consisting solely of structures that are forbidden
    from occurring.
    """

    def is_molecule_forbidden(self, molecule):
        """
        Return ``True`` if the given :class:`Molecule` object `molecule`
        contains forbidden functionality, or ``False`` if not. Labeled atoms
        on the forbidden structures and the molecule are honored.
        """
        from rmgpy.species import Species

        for entry in self.entries.values():
            if isinstance(entry.item, Molecule) or isinstance(entry.item, Species):
                # Perform an isomorphism check
                if entry.item.is_isomorphic(molecule):
                    return True
            elif isinstance(entry.item, Group):
                # We need to do subgraph isomorphism
                entry_labeled_atoms = entry.item.get_all_labeled_atoms()
                molecule_labeled_atoms = molecule.get_all_labeled_atoms()
                for label in entry_labeled_atoms:
                    # all group labels must be present in the molecule
                    if label not in molecule_labeled_atoms: break
                else:
                    if molecule.is_subgraph_isomorphic(entry.item, generate_initial_map=True):
                        return True
            else:
                raise NotImplementedError('Checking is only implemented for forbidden Groups, Molecule, and Species.')

        # Until we have more thermodynamic data of molecular ions we will forbid them
        # if molecule.get_net_charge() != 0:
        #     return True

        return False

    def load_old(self, path):
        """
        Load an old forbidden structures file from the location `path` on disk.
        """
        self.load_old_dictionary(path, pattern=True)
        return self

    def save_old(self, path):
        """
        Save an old forbidden structures file to the location `path` on disk.
        """
        self.save_old_dictionary(path)

    def load_entry(self, label, group=None, molecule=None, species=None, shortDesc='', longDesc='',
                   metal=None, facet=None, site=None):
        """
        Load an entry from the forbidden structures database. This method is
        automatically called during loading of the forbidden structures 
        database.
        """
        from rmgpy.species import Species

        if sum([bool(molecule), bool(group), bool(species)]) != 1:
            raise DatabaseError('A forbidden group should be defined with exactly one item from '
                                'the following options: group, molecule, or species.')
        if molecule is not None:
            item = Molecule().from_adjacency_list(molecule)
        elif species is not None:
            item = Species().from_adjacency_list(species)
            item.generate_resonance_structures()
        elif group is not None:
            if (group[0:3].upper() == 'OR{' or
                    group[0:4].upper() == 'AND{' or
                    group[0:7].upper() == 'NOT OR{' or
                    group[0:8].upper() == 'NOT AND{'):
                item = make_logic_node(group)
            else:
                item = Group().from_adjacency_list(group)
        self.entries[label] = Entry(
            label=label,
            item=item,
            short_desc=shortDesc,
            long_desc=longDesc.strip(),
            metal=metal,
            facet=facet,
            site=site,
        )

    def save_entry(self, f, entry, name='entry'):
        """
        Save an `entry` from the forbidden structures database. This method is
        automatically called during saving of the forbidden structures 
        database.
        """

        f.write('{0}(\n'.format(name))
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

        if entry.metal:
            f.write('    metal = "{0}",\n'.format(entry.metal))
        if entry.facet:
            f.write('    facet = "{0}",\n'.format(entry.facet))
        if entry.site:
            f.write('    site = "{0}",\n'.format(entry.site))

        f.write(f'    shortDesc = """{entry.short_desc.strip()}""",\n')
        f.write(f'    longDesc = \n"""\n{entry.long_desc.strip()}\n""",\n')

        f.write(')\n\n')
