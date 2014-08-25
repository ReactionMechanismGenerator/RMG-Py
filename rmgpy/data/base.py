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
Contains classes and functions for working with the various RMG databases. In
particular, this module is devoted to functionality that is common across all
components of the RMG database.
"""

import os
import logging
import re
import codecs
try:
    from collections import OrderedDict
except ImportError:
    logging.warning("Upgrade to Python 2.7 or later to ensure your database entries are read and written in the same order each time!")
    OrderedDict = dict
from rmgpy.molecule import Molecule, Group, InvalidAdjacencyListError

from reference import Reference, Article, Book, Thesis

################################################################################

class DatabaseError(Exception):
    """
    A exception that occurs when working with an RMG database. Pass a string
    giving specifics about the exceptional behavior.
    """
    pass

################################################################################

class Entry:
    """
    A class for representing individual records in an RMG database. Each entry
    in the database associates a chemical item (generally a species, functional
    group, or reaction) with a piece of data corresponding to that item. A
    significant amount of metadata can also be stored with each entry.

    The attributes are:

    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `index`             A unique nonnegative integer index for the entry
    `label`             A unique string identifier for the entry (or '' if not used)
    `item`              The item that this entry represents
    `parent`            The parent of the entry in the hierarchy (or ``None`` if not used)
    `children`          A list of the children of the entry in the hierarchy (or ``None`` if not used)
    `data`              The data to associate with the item
    `reference`         A :class:`Reference` object containing bibliographic reference information to the source of the data
    `referenceType`     The way the data was determined: ``'theoretical'``, ``'experimental'``, or ``'review'``
    `shortDesc`         A brief (one-line) description of the data
    `longDesc`          A long, verbose description of the data
    `rank`              An integer indicating the degree of confidence in the entry data, or ``None`` if not used
    =================== ========================================================

    """

    def __init__(self,
                 index=-1,
                 label='',
                 item=None,
                 parent=None,
                 children=None,
                 data=None,
                 reference=None,
                 referenceType='',
                 shortDesc='',
                 longDesc='',
                 rank=None,
                 ):
        self.index = index
        self.label = label
        self.item = item
        self.parent = parent
        self.children = children or []
        self.data = data
        self.reference = reference
        self.referenceType = referenceType
        self.shortDesc = shortDesc
        self.longDesc = longDesc
        self.rank = rank

    def __str__(self):
        return self.label

    def __repr__(self):
        return '<Entry index={0:d} label="{1}">'.format(self.index, self.label)

################################################################################

class Database:
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

    You must derive from this class and implement the :meth:`loadEntry`,
    :meth:`saveEntry`, :meth:`processOldLibraryEntry`, and
    :meth:`generateOldLibraryEntry` methods in order to load and save from the
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
                 shortDesc='',
                 longDesc='',
                 ):
        self.entries = OrderedDict(entries or {})
        self.top = top or []
        self.label = label
        self.name = name
        self.shortDesc = shortDesc
        self.longDesc = longDesc

    def load(self, path, local_context=None, global_context=None):
        """
        Load an RMG-style database from the file at location `path` on disk.
        The `entryName` parameter specifies the identifier used for each data
        entry. The parameters `local_context` and `global_context` are used to
        provide specialized mapping of identifiers in the input file to
        corresponding functions to evaluate. This method will automatically add
        a few identifiers required by all data entries, so you don't need to
        provide these.
        """

        # Collision efficiencies are in SMILES format, so we'll need RDKit
        # to convert them to Molecule objects
        # Do the import here to ensure it is imported from a pure Python
        # environment (as opposed to a Cythonized environment, which is not
        # allowed during an exec() call)
        from rdkit import Chem

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
        local_context['entry'] = self.loadEntry
        local_context['tree'] = self.__loadTree
        local_context['name'] = self.name
        local_context['shortDesc'] = self.shortDesc
        local_context['longDesc'] = self.longDesc
        # add in anything from the Class level dictionary.
        for key, value in Database.local_context.iteritems():
            local_context[key]=value
        
        # Process the file
        f = open(path, 'r')
        try:
            exec f in global_context, local_context
        except Exception, e:
            logging.error('Error while reading database {0!r}.'.format(path))
            raise
        f.close()

        # Extract the database metadata
        self.name = local_context['name']
        self.shortDesc = local_context['shortDesc']
        self.longDesc = local_context['longDesc'].strip()
        
        # Return the loaded database (to allow for Database().load() syntax)
        return self

    def getEntriesToSave(self):
        """
        Return a sorted list of the entries in this database that should be
        saved to the output file.
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
            entries = self.entries.values()
            entries.sort(key=lambda x: (x.index))
        return entries
    
    def getSpecies(self, path):
        """
        Load the dictionary containing all of the species in a kinetics library or depository.
        """
        from rmgpy.species import Species
        speciesDict = {}
        with open(path, 'r') as f:
            adjlist = ''
            for line in f:
                if line.strip() == '' and adjlist.strip() != '':
                    # Finish this adjacency list
                    species = Species().fromAdjacencyList(adjlist)
                    species.generateResonanceIsomers()
                    label = species.label
                    if label in speciesDict:
                        raise DatabaseError('Species label "{0}" used for multiple species in {1}.'.format(label, str(self)))
                    speciesDict[label] = species
                    adjlist = ''
                else:
                    adjlist += line
        
        return speciesDict
    
    def saveDictionary(self, path):
        """
        Extract species from all entries associated with a kinetics library or depository and save them 
        to the path given.
        """
        try:
            os.makedirs(os.path.dirname(path))
        except OSError:
            pass
        # Extract species from all the entries
        speciesDict = {}
        entries = self.entries.values()
        for entry in entries:
            for reactant in entry.item.reactants:
                if reactant.label not in speciesDict:
                    speciesDict[reactant.label] = reactant
                elif not reactant.isIsomorphic(speciesDict[reactant.label]):
                    print reactant.molecule[0].toAdjacencyList()
                    print speciesDict[reactant.label].molecule[0].toAdjacencyList()
                    speciesDict[reactant.label] = reactant
                    raise DatabaseError('Species label "{0}" used for multiple species in {1}.'.format(reactant.label, str(self)))
            for product in entry.item.products:
                if product.label not in speciesDict:
                    speciesDict[product.label] = product
                elif not product.isIsomorphic(speciesDict[product.label]):
                    print product.molecule[0].toAdjacencyList()
                    print speciesDict[product.label].molecule[0].toAdjacencyList()
                    raise DatabaseError('Species label "{0}" used for multiple species in {1}.'.format(product.label, str(self)))
            
        with open(path, 'w') as f:
            for label in speciesDict.keys():
                f.write(speciesDict[label].molecule[0].toAdjacencyList(label=label, removeH=False))
                f.write('\n')

    def save(self, path):
        """
        Save the current database to the file at location `path` on disk. The
        optional `entryName` parameter specifies the identifier used for each
        data entry.
        """
        try:
            os.makedirs(os.path.dirname(path))
        except OSError:
            pass
        entries = self.getEntriesToSave()

        f = codecs.open(path, 'w', 'utf-8')
        f.write('#!/usr/bin/env python\n')
        f.write('# encoding: utf-8\n\n')
        f.write('name = "{0}"\n'.format(self.name))
        f.write('shortDesc = u"{0}"\n'.format(self.shortDesc))
        f.write('longDesc = u"""\n')
        f.write(self.longDesc.strip() + '\n')
        f.write('"""\n')
        
        for entry in entries:
            self.saveEntry(f, entry)

        # Write the tree
        if len(self.top) > 0:
            f.write('tree(\n')
            f.write('"""\n')
            f.write(self.generateOldTree(self.top, 1))
            f.write('"""\n')
            f.write(')\n\n')

        f.close()

    def loadOld(self, dictstr, treestr, libstr, numParameters, numLabels=1, pattern=True):
        """
        Load a dictionary-tree-library based database. The database is stored
        in three files: `dictstr` is the path to the dictionary, `treestr` to
        the tree, and `libstr` to the library. The tree is optional, and should
        be set to '' if not desired.
        """

        # Load dictionary, library, and (optionally) tree
        try:
            self.loadOldDictionary(dictstr, pattern)
        except Exception, e:
            logging.error('Error while reading database {0!r}.'.format(os.path.dirname(dictstr)))
            raise
        
        try:
            if treestr != '': self.loadOldTree(treestr)
        except Exception, e:
            logging.error('Error while reading database {0!r}.'.format(os.path.dirname(treestr)))
            raise
        
        try:
            self.loadOldLibrary(libstr, numParameters, numLabels)
        except Exception, e:
            logging.error('Error while reading database {0!r}.'.format(os.path.dirname(libstr)))
            raise
          
        return self

    def loadOldDictionary(self, path, pattern):
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

        fdict=None
        # Process the dictionary file
        try:
            fdict = open(path, 'r')
            for line in fdict:
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
                    line = removeCommentFromLine(line).strip()
                    if len(line) > 0:
                        record += line + '\n'
            # process the last record! (after end of for loop)
            # Label is first line of record
            if record:
                label = record.splitlines()[0]
                # Add record to dictionary
                self.entries[label] = Entry(label=label, item=record)
        except DatabaseError, e:
            logging.exception(str(e))
            raise
        except IOError, e:
            logging.exception('Database dictionary file "' + e.filename + '" not found.')
            raise
        finally:
            if fdict: fdict.close()

        # Convert the records in the dictionary to Molecule, Group, or
        # logical objects
        try:
            for label in self.entries:
                record = self.entries[label].item
                lines = record.splitlines()
                # If record is a logical node, make it into one.
                if re.match("(?i)\s*(NOT\s)?\s*(OR|AND|UNION)\s*(\{.*\})", lines[1]):
                    self.entries[label].item = makeLogicNode(' '.join(lines[1:]) )
                # Otherwise convert adjacency list to molecule or pattern
                elif pattern:
                    self.entries[label].item = Group().fromAdjacencyList(record)
                else:
                    self.entries[label].item = Molecule().fromAdjacencyList(record,saturateH=True)
        except InvalidAdjacencyListError, e:
            logging.error('Error while loading old-style dictionary "{0}"'.format(path))
            logging.error('Error occurred while parsing adjacency list "{0}"'.format(label))
            raise

    def __loadTree(self, tree):
        """
        Parse an old-style RMG tree located at `tree`. An RMG tree is an n-ary
        tree representing the hierarchy of items in the dictionary.
        """

        if len(self.entries) == 0:
            raise DatabaseError("Load the dictionary before you load the tree.")

        # should match '  L3 : foo_bar '  and 'L3:foo_bar'
        parser = re.compile('^\s*L(?P<level>\d+)\s*:\s*(?P<label>\S+)')

        parents = [None]
        for line in tree.splitlines():
            line = removeCommentFromLine(line).strip()
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
                        parent = parents[level-1]

                if parent is not None: parent = self.entries[parent]
                try:
                    entry = self.entries[label]
                except KeyError:
                    raise DatabaseError('Unable to find entry "{0}" from tree in dictionary.'.format(label))
                
                if isinstance(parent, str):
                    raise DatabaseError('Unable to find parent entry "{0}" of entry "{1}" in tree.'.format(parent, label))

                # Update the parent and children of the nodes accordingly
                if parent is not None:
                    entry.parent = parent
                    parent.children.append(entry)
                else:
                    entry.parent = None
                    self.top.append(entry)

                # Add node to list of parents for subsequent iteration
                parents.append(label)


    def loadOldTree(self, path):
        """
        Parse an old-style RMG database tree located at `path`. An RMG
        tree is an n-ary tree representing the hierarchy of items in the
        dictionary.
        """
        tree = []
        try:
            ftree = open(path, 'r')
            tree = ftree.read()

        except IOError, e:
            logging.exception('Database tree file "' + e.filename + '" not found.')
        finally:
            ftree.close()

        self.__loadTree(tree)

    def loadOldLibrary(self, path, numParameters, numLabels=1):
        """
        Parse an RMG database library located at `path`.
        """

        if len(self.entries) == 0:
            raise DatabaseError("Load the dictionary before you load the library.")

        entries = self.parseOldLibrary(path, numParameters, numLabels)

        # Load the parsed entries into the database, skipping duplicate entries
        skippedCount = 0
        for index, label, parameters, comment in entries:
            if label not in self.entries:
                raise DatabaseError('Entry {0!r} in library was not found in dictionary.'.format(label))
            if self.entries[label].index != -1:
                # The entry is a duplicate, so skip it
                logging.debug("There was already something labeled {0} in the {1} library. Ignoring '{2}' ({3})".format(label, self.label, index, parameters))
                skippedCount += 1
            else:
                # The entry is not a duplicate
                self.entries[label].index = index
                self.entries[label].data = parameters
                self.entries[label].shortDesc = comment
        if skippedCount > 0:
            logging.warning("Skipped {0:d} duplicate entries in {1} library.".format(skippedCount, self.label))

        # Make sure each entry with data has a nonnegative index
        entries2 = self.entries.values()
        entries2.sort(key=lambda entry: entry.index)
        index = entries2[-1].index + 1
        if index < 1: index = 1
        for index0, label, parameters, comment in entries:
            if self.entries[label].index < 0:
                self.entries[label].index = index
                index += 1

    def parseOldLibrary(self, path, numParameters, numLabels=1):
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
                line = removeCommentFromLine(line).strip()
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
                    label = self.__hashLabels(info[offset:offset+numLabels])
                    offset += numLabels
                    # Extract numeric parameter(s) or label of node with data to use
                    if numParameters < 0:
                        parameters = self.processOldLibraryEntry(info[offset:])
                        comment = ''
                    else:
                        try:
                            parameters = self.processOldLibraryEntry(info[offset:offset+numParameters])
                            offset += numParameters
                        except (IndexError, ValueError), e:
                            parameters = info[offset]
                            offset += 1
                        # Remaining part of string is comment
                        comment = ' '.join(info[offset:])
                        comment = comment.strip('"')

                    entries.append((index, label, parameters, comment))

        except DatabaseError, e:
            logging.exception(str(e))
            logging.exception("path = '{0}'".format(path))
            logging.exception("line = '{0}'".format(line))
            raise
        except IOError, e:
            logging.exception('Database library file "' + e.filename + '" not found.')
            raise
        finally:
            if flib: flib.close()

        return entries

    def saveOld(self, dictstr, treestr, libstr):
        """
        Save the current database to a set of text files using the old-style
        syntax.
        """
        self.saveOldDictionary(dictstr)
        if treestr != '':
            self.saveOldTree(treestr)

        # RMG-Java does not require a frequencies_groups/Library.txt file to
        # operate, but errors are raised upon importing to Py if this file is
        # not found. This check prevents the placeholder from being discarded.
        if 'StatesGroups' not in self.__class__.__name__:
            self.saveOldLibrary(libstr)

    def saveOldDictionary(self, path):
        """
        Save the current database dictionary to a text file using the old-style
        syntax.
        """

        entries = []
        entriesNotInTree = []
        # If we have tree information, save the dictionary in the same order as
        # the tree (so that it saves in the same order each time)
        def getLogicNodeComponents(entry_or_item):
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
                        nodes.extend(getLogicNodeComponents(child))
                    else:
                        nodes.extend(getLogicNodeComponents(self.entries[child]))
                return nodes
            else:
                return [entry]

        if len(self.top) > 0:
            for entry in self.top:
                entries.extend(getLogicNodeComponents(entry))
                for descendant in self.descendants(entry):
                    for entry2 in getLogicNodeComponents(descendant):
                        if entry2 not in entries:
                            entries.append(entry2)
                
            # Don't forget entries that aren't in the tree
            for entry in self.entries.values():
                if entry not in entries:
                    entriesNotInTree.append(entry)
            entriesNotInTree.sort(key=lambda x: (x.index, x.label))
        # Otherwise save the dictionary in any order
        else:
            # Save the library in order by index
            entries = self.entries.values()
            entries.sort(key=lambda x: (x.index, x.label))

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
                    f.write(entry.item.toAdjacencyList(removeH=False) + '\n')
                elif isinstance(entry.item, Group):
                    f.write(entry.item.toAdjacencyList().replace('{2S,2T}','2') + '\n')
                elif isinstance(entry.item, LogicOr):
                    f.write('{0}\n\n'.format(entry.item).replace('OR{', 'Union {'))
                elif entry.label[0:7] == 'Others-':
                    assert isinstance(entry.item, LogicNode)
                    f.write('{0}\n\n'.format(entry.item))
                else:
                    raise DatabaseError('Unexpected item with label {0} encountered in dictionary while attempting to save.'.format(entry.label))
            
            def comment(s):
                "Return the string, with each line prefixed with '// '"
                return '\n'.join('// '+line if line else '' for line in s.split('\n'))
            if entriesNotInTree:
                f.write(comment("These entries do not appear in the tree:\n\n"))
            for entry in entriesNotInTree:
                f.write(comment(entry.label + '\n'))
                if isinstance(entry.item, Molecule):
                    f.write(comment(entry.item.toAdjacencyList(removeH=False) + '\n'))
                elif isinstance(entry.item, Group):
                    f.write(comment(entry.item.toAdjacencyList().replace('{2S,2T}','2') + '\n'))
                elif isinstance(entry.item, LogicOr):
                    f.write(comment('{0}\n\n'.format(entry.item).replace('OR{', 'Union {')))
                elif entry.label[0:7] == 'Others-':
                    assert isinstance(entry.item, LogicNode)
                    f.write(comment('{0}\n\n'.format(entry.item)))
                else:
                    raise DatabaseError('Unexpected item with label {0} encountered in dictionary while attempting to save.'.format(entry.label))
           
           
            f.close()
        except IOError, e:
            logging.exception('Unable to save old-style dictionary to "{0}".'.format(os.path.abspath(path)))
            raise

    def generateOldTree(self, entries, level):
        """
        Generate a multi-line string representation of the current tree using
        the old-style syntax.
        """
        string = ''
        for entry in entries:
            # Write current node
            string += '{0}L{1:d}: {2}\n'.format('    ' * (level-1), level, entry.label)
            # Recursively descend children (depth-first)
            string += self.generateOldTree(entry.children, level+1)
        return string

    def saveOldTree(self, path):
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
            f.write(self.generateOldTree(self.top, 1))
            f.close()
        except IOError, e:
            logging.exception('Unable to save old-style tree to "{0}".'.format(os.path.abspath(path)))
            raise

    def saveOldLibrary(self, path):
        """
        Save the current database library to a text file using the old-style
        syntax.
        """
        try:
            # Save the library in order by index
            entries = self.entries.values()
            entries.sort(key=lambda x: x.index)
            
            f = codecs.open(path, 'w',  'utf-8')
            records = []
            for entry in entries:
                if entry.data is not None:
                    data = entry.data
                    if not isinstance(data, str):
                        data = self.generateOldLibraryEntry(data)
                    records.append((entry.index, [entry.label], data, entry.shortDesc))

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
                if isinstance(data, basestring):
                    f.write('{:s} '.format(data))
                else:
                    f.write('{:s} '.format(' '.join(['{:<10g}'.format(d) for d in data])))
                f.write(u'    {:s}\n'.format(comment))
            f.close()
        except IOError, e:
            logging.exception('Unable to save old-style library to "{0}".'.format(os.path.abspath(path)))
            raise

    def __hashLabels(self, labels):
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
        if isinstance(node, str): node = self.entries[node]
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
        if isinstance(node, str): node = self.entries[node]
        descendants = []
        for child in node.children:
            descendants.append(child)
            descendants.extend(self.descendants(child))
        return descendants
    
    def matchNodeToNode(self, node, nodeOther):
        """ 
        Return `True` if `node` and `nodeOther` are identical.  Otherwise, return `False`.
        Both `node` and `nodeOther` must be Entry types with items containing Group or LogicNode types.
        """
        if isinstance(node.item, Group) and isinstance(nodeOther.item, Group):
            return self.matchNodeToStructure(node,nodeOther.item, atoms=nodeOther.item.getLabeledAtoms()) and self.matchNodeToStructure(nodeOther,node.item,atoms=node.item.getLabeledAtoms())
        elif isinstance(node.item,LogicOr) and isinstance(nodeOther.item,LogicOr):
            return node.item.matchLogicOr(nodeOther.item)
        else:
            # Assume nonmatching
            return False
        
    def matchNodeToChild(self, parentNode, childNode):        
        """ 
        Return `True` if `parentNode` is a parent of `childNode`.  Otherwise, return `False`.
        Both `parentNode` and `childNode` must be Entry types with items containing Group or LogicNode types.
        If `parentNode` and `childNode` are identical, the function will also return `False`.
        """
        
        if isinstance(parentNode.item, Group) and isinstance(childNode.item, Group):
            if self.matchNodeToStructure(parentNode,childNode.item, atoms=childNode.item.getLabeledAtoms()) is True:
                if self.matchNodeToStructure(childNode,parentNode.item, atoms=parentNode.item.getLabeledAtoms()) is False:
                    return True                
            return False
        
        #If the parentNode is a Group and the childNode is a LogicOr there is nothing to check,
        #so it gets an automatic pass. However, we do need to check that everything down this
        #family line is consistent, which is done in the databaseTest unitTest
        elif isinstance(parentNode.item, Group) and isinstance(childNode.item, LogicOr):
            return True
        
        elif isinstance(parentNode.item,LogicOr):
            return childNode.label in parentNode.item.components

    def matchNodeToStructure(self, node, structure, atoms):
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
        be found in node.  However the reverse is not true.
        
        Usage: node = either an Entry or a key in the self.entries dictionary which has
                      a Group or LogicNode as its Entry.item
               structure = a Group or a Molecule
               atoms = dictionary of {label: atom} in the structure.  A possible dictionary
                       is the one produced by structure.getLabeledAtoms()
        """
        if isinstance(node, str): node = self.entries[node]
        group = node.item
        if isinstance(group, LogicNode):
            return group.matchToStructure(self, structure, atoms)
        else:
            # try to pair up labeled atoms
            centers = group.getLabeledAtoms()
            initialMap = {}
            for label in centers.keys():
                # Make sure the labels are in both group and structure.
                if label not in atoms:
                    logging.log(0, "Label {0} is in group {1} but not in structure".format(label, node))
                    continue # with the next label - ring structures might not have all labeled atoms
                    # return False # force it to have all the labeled atoms
                center = centers[label]
                atom = atoms[label]
                # Make sure labels actually point to atoms.
                if center is None or atom is None:
                    return False
                if isinstance(center, list):
                    center = center[0]
                # Semantic check #1: atoms with same label are equivalent
                elif not atom.isSpecificCaseOf(center):
                    return False
                # Semantic check #2: labeled atoms that share bond in the group (node)
                # also share equivalent (or more specific) bond in the structure
                for atom2, atom1 in initialMap.iteritems():
                    if group.hasBond(center, atom1) and structure.hasBond(atom, atom2):
                        bond1 = group.getBond(center, atom1)   # bond1 is group
                        bond2 = structure.getBond(atom, atom2) # bond2 is structure
                        if not bond2.isSpecificCaseOf(bond1):
                            return False
                    elif group.hasBond(center, atom1): # but structure doesn't
                        return False
                    # but we don't mind if...
                    elif structure.hasBond(atom, atom2): # but group doesn't
                        logging.debug("We don't mind that structure "+ str(structure) +
                            " has bond but group {0} doesn't".format(node))
                # Passed semantic checks, so add to maps of already-matched atoms
                initialMap[atom] = center
            # Labeled atoms in the structure that are not in the group should
            # not be considered in the isomorphism check, so remove them temporarily
            # Without this we would hit a lot of nodes that are ambiguous
            removedAtoms = []
            for label, atom in structure.getLabeledAtoms().iteritems():
                if label not in centers:
                    removedAtoms.append(atom)
                    structure.atoms.remove(atom)
            # use mapped (labeled) atoms to try to match subgraph
            result = structure.isSubgraphIsomorphic(group, initialMap)
            # Restore atoms removed in previous step
            for atom in removedAtoms:
                structure.atoms.append(atom)
            return result

    def descendTree(self, structure, atoms, root=None):
        """
        Descend the tree in search of the functional group node that best
        matches the local structure around `atoms` in `structure`.

        If root=None then uses the first matching top node.

        Returns None if there is no matching root.
        """

        if root is None:
            for root in self.top:
                if self.matchNodeToStructure(root, structure, atoms):
                    break # We've found a matching root
            else: # didn't break - matched no top nodes
                return None
        elif not self.matchNodeToStructure(root, structure, atoms):
            return None
        
        next = []
        for child in root.children:
            if self.matchNodeToStructure(child, structure, atoms):
                next.append(child)

        if len(next) == 1:
            return self.descendTree(structure, atoms, next[0])
        elif len(next) == 0:
            if len(root.children) > 0 and root.children[-1].label.startswith('Others-'):
                return root.children[-1]
            else:
                return root
        else:
            #logging.warning('For {0}, a node {1} with overlapping children {2} was encountered in tree with top level nodes {3}. Assuming the first match is the better one.'.format(structure, root, next, self.top))
            return self.descendTree(structure, atoms, next[0])

################################################################################

class LogicNode:
    """
    A base class for AND and OR logic nodes.
    """

    symbol="<TBD>" # To be redefined by subclass

    def __init__(self,items,invert):
        self.components = []
        for item in items:
            if re.match("(?i)\s*(NOT\s)?\s*(OR|AND|UNION)\s*(\{.*\})",item):
                component = makeLogicNode(item)
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

    def matchToStructure(self,database,structure,atoms):
        """
        Does this node in the given database match the given structure with the labeled atoms?
        """
        for node in self.components:
            if isinstance(node,LogicNode):
                match = node.matchToStructure(database,structure,atoms)
            else:
                match = database.matchNodeToStructure(node, structure, atoms)
            if match:
                return True != self.invert
        return False != self.invert

    def matchLogicOr(self, other):
        """
        Is other the same LogicOr group as self?
        """
        if len(self.components)!=len(other.components):
            return False
        else:
            for node in self.components:
                if node not in other.components:
                    return False
        return True
        
    def getPossibleStructures(self, entries):
        """
        Return a list of the possible structures below this node.
        """
        if self.invert: raise NotImplementedError("Finding possible structures of NOT OR nodes not implemented.")
        structures = []
        for item in self.components:
            struct = entries[item].item
            if isinstance(struct, LogicNode):
                structures.extend(struct.getPossibleStructures(entries))
            else:
                structures.append(struct)
        for struct in structures: # check this worked
            assert isinstance(struct,Group)
        return structures

class LogicAnd(LogicNode):
    """A logical AND node. Structure must match all components."""

    symbol = "AND"

    def matchToStructure(self,database,structure,atoms):
        """
        Does this node in the given database match the given structure with the labeled atoms?
        """
        for node in self.components:
            if isinstance(node,LogicNode):
                match = node.matchToStructure(database,structure,atoms)
            else:
                match = database.matchNodeToStructure(node, structure, atoms)
            if not match:
                return False != self.invert
        return True != self.invert

def makeLogicNode(string):
    """
    Creates and returns a node in the tree which is a logic node.

    String should be of the form:

    * OR{}
    * AND{}
    * NOT OR{}
    * NOT AND{}

    And the returned object will be of class LogicOr or LogicAnd
    """

    match = re.match("(?i)\s*(NOT\s)?\s*(OR|AND|UNION)\s*(\{.*\})",string)  # the (?i) makes it case-insensitive
    if not match:
        raise Exception("Unexpected string for Logic Node: {0}".format(string))

    if match.group(1): invert = True
    else: invert = False

    logic = match.group(2)  # OR or AND (or Union)

    contents = match.group(3).strip()
    while contents.startswith('{'):
        if not contents.endswith('}'):
            raise Exception("Unbalanced braces in Logic Node: {0}".format(string))
        contents = contents[1:-1]

    items=[]
    chars=[]
    brace_depth = 0
    for character in contents:
        if character == '{':
            brace_depth += 1
        if character == '}':
            brace_depth -= 1
        if character == ',' and brace_depth == 0:
            items.append(''.join(chars).lstrip().rstrip() )
            chars = []
        else:
            chars.append(character)
    if chars: # add last item
        items.append(''.join(chars).lstrip().rstrip() )
    if brace_depth != 0: raise Exception("Unbalanced braces in Logic Node: {0}".format(string))

    if logic.upper() in ['OR', 'UNION']:
        return LogicOr(items, invert)
    if logic == 'AND':
        return LogicAnd(items, invert)

    raise Exception("Could not create Logic Node from {0}".format(string))

################################################################################

def removeCommentFromLine(line):
    """
    Remove a C++/Java style comment from a line of text. This refers
    particularly to comments that begin with a double-slash '//' and continue
    to the end of the line.
    """

    index = line.find('//')
    if index >= 0:
        line = line[0:index]
    return line

def splitLineAndComment(line):
    """
    Returns a tuple(line, comment) based on a '//' comment delimiter.
    
    Either `line` or `comment` may be ''.
    Does not strip whitespace, nor remove more than two slashes.
    """
    split = line.split('//',1)
    if len(split) == 1:
        return (split[0],'')
    else:
        return tuple(split)

def getAllCombinations(nodeLists):
    """
    Generate a list of all possible combinations of items in the list of
    lists `nodeLists`. Each combination takes one item from each list
    contained within `nodeLists`. The order of items in the returned lists
    reflects the order of lists in `nodeLists`. For example, if `nodeLists` was
    [[A, B, C], [N], [X, Y]], the returned combinations would be
    [[A, N, X], [A, N, Y], [B, N, X], [B, N, Y], [C, N, X], [C, N, Y]].
    """

    items = [[]]
    for nodeList in nodeLists:
        items = [ item + [node] for node in nodeList for item in items ]

    return items

################################################################################

class ForbiddenStructureException(Exception):
    """
    Made a forbidden structure.
    """
    pass

class ForbiddenStructures(Database):
    """
    A database consisting solely of structures that are forbidden
    from occurring.
    """
    
    def isMoleculeForbidden(self, molecule):
        """
        Return ``True`` if the given :class:`Molecule` object `molecule`
        contains forbidden functionality, or ``False`` if not. Labeled atoms
        on the forbidden structures and the molecule are honored.
        """
        for entry in self.entries.values():
            entryLabeledAtoms = entry.item.getLabeledAtoms()
            moleculeLabeledAtoms = molecule.getLabeledAtoms()
            initialMap = {}
            for label in entryLabeledAtoms:
                if label not in moleculeLabeledAtoms: continue
                initialMap[moleculeLabeledAtoms[label]] = entryLabeledAtoms[label]
            if molecule.isMappingValid(entry.item, initialMap) and molecule.isSubgraphIsomorphic(entry.item, initialMap):
                return True
            
        # Until we have more thermodynamic data of molecular ions we will forbid them
        molecule_charge = 0
        for atom in molecule.atoms:
            molecule_charge += atom.charge
        if molecule_charge != 0:
            return True
        
        return False
    
    def loadOld(self, path):
        """
        Load an old forbidden structures file from the location `path` on disk.
        """
        self.loadOldDictionary(path, pattern=True)
        return self

    def saveOld(self, path):
        """
        Save an old forbidden structures file to the location `path` on disk.
        """
        self.saveOldDictionary(path)

    def loadEntry(self, label, molecule=None, group=None, shortDesc='', longDesc=''):
        """
        Load an entry from the forbidden structures database. This method is
        automatically called during loading of the forbidden structures 
        database.
        """
        assert molecule is not None or group is not None
        assert not (molecule is not None and group is not None)
        if molecule is not None:
            item = Molecule.fromAdjacencyList(molecule)
        elif group is not None:
            if ( group[0:3].upper() == 'OR{' or
                 group[0:4].upper() == 'AND{' or
                 group[0:7].upper() == 'NOT OR{' or
                 group[0:8].upper() == 'NOT AND{'
                ):
                item = makeLogicNode(group)
            else:
                item = Group().fromAdjacencyList(group)
        self.entries[label] = Entry(
            label = label,
            item = item,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
        )
    
    def saveEntry(self, f, entry, name='entry'):
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
            f.write(entry.item.toAdjacencyList(removeH=False))
            f.write('""",\n')
        elif isinstance(entry.item, Group):
            f.write('    group = \n')
            f.write('"""\n')
            f.write(entry.item.toAdjacencyList())
            f.write('""",\n')
        else:
            f.write('    group = "{0}",\n'.format(entry.item))

        f.write('    shortDesc = u"""{0}""",\n'.format(entry.shortDesc))
        f.write('    longDesc = \n')
        f.write('u"""\n')
        f.write(entry.longDesc.strip() + "\n")
        f.write('""",\n')

        f.write(')\n\n')
