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


import os.path
import logging
from copy import deepcopy
import numpy

import rmgpy.constants as constants
from rmgpy.kinetics import Arrhenius, ArrheniusEP, ThirdBody, Lindemann, Troe, \
                           PDepArrhenius, MultiArrhenius, MultiPDepArrhenius, \
                           Chebyshev, KineticsData
from rmgpy.molecule import Molecule, Group
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.data.base import LogicNode, DatabaseError

from .family import  KineticsFamily
from .library import LibraryReaction, KineticsLibrary
from .common import filterReactions

################################################################################

class KineticsDatabase(object):
    """
    A class for working with the RMG kinetics database.
    """

    def __init__(self):
        self.recommendedFamilies = {}
        self.families = {}
        self.libraries = {}
        self.libraryOrder = []     # a list of tuples in the format ('library_label', LibraryType),
                                   # where LibraryType is set to either 'Reaction Library' or 'Seed'.  
        self.local_context = {
            'KineticsData': KineticsData,
            'Arrhenius': Arrhenius,
            'ArrheniusEP': ArrheniusEP,
            'MultiArrhenius': MultiArrhenius,
            'MultiPDepArrhenius': MultiPDepArrhenius,
            'PDepArrhenius': PDepArrhenius,
            'Chebyshev': Chebyshev,
            'ThirdBody': ThirdBody,
            'Lindemann': Lindemann,
            'Troe': Troe,
            'R': constants.R,
        }
        self.global_context = {}

    def __reduce__(self):
        """
        A helper function used when pickling a KineticsDatabase object.
        """
        d = {
            'families': self.families,
            'libraries': self.libraries,
            'libraryOrder': self.libraryOrder,
        }
        return (KineticsDatabase, (), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling a KineticsDatabase object.
        """
        self.families = d['families']
        self.libraries = d['libraries']
        self.libraryOrder = d['libraryOrder']

    def load(self, path, families=None, libraries=None, depositories=None):
        """
        Load the kinetics database from the given `path` on disk, where `path`
        points to the top-level folder of the families database.
        """
        self.loadRecommendedFamiliesList(os.path.join(path, 'families', 'recommended.py')),
        self.loadFamilies(os.path.join(path, 'families'), families, depositories)
        self.loadLibraries(os.path.join(path, 'libraries'), libraries)

    def loadRecommendedFamiliesList(self, filepath):
        """
        Load the list of recommended families from the given file
        
        The file is usually 'kinetics/families/recommended.py'.
        This is stored as a dictionary of True or False values (checked here),
        and should contain entries for every available family (checked in loadFamilies).
        """
        try:
            global_context = {}
            global_context['__builtins__'] = None
            global_context['True'] = True
            global_context['False'] = False
            local_context = {}
            local_context['__builtins__'] = None
            f = open(filepath, 'r')
            exec f in global_context, local_context
            f.close()
            self.recommendedFamilies = local_context['recommendedFamilies']
        except Exception, e:
            raise DatabaseError('Error while reading list of recommended families from {0}/recommended.py.\n{1}'.format(filepath,e))
        for recommended in self.recommendedFamilies.values():
            if not isinstance(recommended, bool):
                raise DatabaseError("recommendedFamilies dictionary should contain only True or False values")

    def loadFamilies(self, path, families=None, depositories=None):
        """
        Load the kinetics families from the given `path` on disk, where `path`
        points to the top-level folder of the kinetics families.
        """
        
        familiesToLoad = []
        for (root, dirs, files) in os.walk(os.path.join(path)):
            if root == path:
                break  # all we wanted was the list of dirs in the base path

        if families == 'default':
            logging.info('Loading default kinetics families from {0}'.format(path))
            for d in dirs:  # load them in directory listing order, like other methods (better than a random dict order)
                try:
                    recommended = self.recommendedFamilies[d]
                except KeyError:
                    raise DatabaseError('Family {0} not found in recommendation list (probably at {1}/recommended.py)'.format(d, path))
                if recommended:
                    familiesToLoad.append(d)
            for label, value in self.recommendedFamilies.iteritems():
                if label not in dirs:
                    raise DatabaseError('Family {0} found (in {1}/recommended.py) not found on disk.'.format(label, path))

        elif families == 'all':
            # All families are loaded
            logging.info('Loading all of the kinetics families from {0}'.format(path))
            for d in dirs:
                familiesToLoad.append(d)
        elif families == 'none':
            logging.info('Not loading any of the kinetics families from {0}'.format(path))
            # Don't load any of the families
            familiesToLoad = []
        elif isinstance(families, list):
            logging.info('Loading the user-specified kinetics families from {0}'.format(path))
            # If all items in the list start with !, all families will be loaded except these
            if len(families) == 0:
                raise DatabaseError('Kinetics families should be a non-empty list, or set to `default`, `all`, or `none`.')
            elif all([label.startswith('!') for label in families]):
                for d in dirs:
                    if '!{0}'.format(d) not in families:
                        familiesToLoad.append(d)
            elif any([label.startswith('!') for label in families]):
                raise DatabaseError('Families list must either all or none have prefix "!", but not a mix.')
            else:  # only the families given will be loaded
                for d in dirs:
                    if d in families:
                        familiesToLoad.append(d)
                for label in families:
                    if label not in dirs:
                        raise DatabaseError('Family {0} not found on disk.'.format(label))
        else:
            raise DatabaseError('Kinetics families was not specified properly.  Should be set to `default`,`all`,`none`, or a list.')
        
        # Now we know what families to load, so let's load them
        self.families = {}
        for label in familiesToLoad:
            familyPath = os.path.join(path, label)
            family = KineticsFamily(label=label)
            family.load(familyPath, self.local_context, self.global_context, depositoryLabels=depositories)
            self.families[label] = family

    def loadLibraries(self, path, libraries=None):
        """
        Load the listed kinetics libraries from the given `path` on disk.
        
        Loads them all if `libraries` list is not specified or `None`.
        The `path` points to the folder of kinetics libraries in the database,
        and the libraries should be in files like :file:`<path>/<library>.py`.
        """
        self.libraries = {}
        
        if libraries is not None:
            for library_name in libraries:
                library_file = os.path.join(path, library_name,'reactions.py')
                if os.path.exists(library_file):
                    logging.info('Loading kinetics library {0} from {1}...'.format(library_name, library_file))
                    library = KineticsLibrary(label=library_name)
                    library.load(library_file, self.local_context, self.global_context)
                    self.libraries[library.label] = library
                else:
                    raise IOError("Couldn't find kinetics library {0}".format(library_file))
            # library order should've been set prior to this, with the given seed mechs and reaction libraries
            assert (len(self.libraryOrder) == len(libraries))
        else:# load all the libraries you can find (this cannot be activated in a normal RMG job.  Only activated when loading the database for other purposes)
            self.libraryOrder = []
            for (root, dirs, files) in os.walk(os.path.join(path)):
                for f in files:
                    name, ext = os.path.splitext(f)
                    if ext.lower() == '.py':
                        library_file = os.path.join(root, f)
                        label=os.path.dirname(library_file)[len(path)+1:]
                        logging.info('Loading kinetics library {0} from {1}...'.format(label, library_file))
                        library = KineticsLibrary(label=label)
                        library.load(library_file, self.local_context, self.global_context)
                        self.libraries[library.label] = library
                        self.libraryOrder.append((library.label,'Reaction Library'))

    def save(self, path):
        """
        Save the kinetics database to the given `path` on disk, where `path`
        points to the top-level folder of the kinetics database.
        """
        path = os.path.abspath(path)
        if not os.path.exists(path): os.mkdir(path)
        self.saveRecommendedFamilies(os.path.join(path, 'families'))
        self.saveFamilies(os.path.join(path, 'families'))
        self.saveLibraries(os.path.join(path, 'libraries'))
        
    def saveRecommendedFamilies(self, path):
        """ 
        Save the list of recommended families in a dictionary stored at 
        `path`/recommended.py
        """
        import codecs
        
        if not os.path.exists(path): os.mkdir(path)
        f = codecs.open(os.path.join(path,'recommended.py'), 'w', 'utf-8')
        f.write('''# This file contains a dictionary of kinetics families.  The families
# set to `True` are recommended by RMG and turned on by default by setting
# kineticsFamilies = 'default' in the RMG input file. Families set to `False` 
# are not turned on by default because the family is severely lacking in data.
# These families should only be turned on with caution.''')
        f.write('\n\n')
        f.write('recommendedFamilies = {\n')
        for label in sorted(self.recommendedFamilies.keys()):
            f.write("'{label}':{value},\n".format(label=label,value=self.recommendedFamilies[label]))
        f.write('}')
        f.close()
        
    def saveFamilies(self, path):
        """
        Save the kinetics families to the given `path` on disk, where `path`
        points to the top-level folder of the kinetics families.
        """
        if not os.path.exists(path): os.mkdir(path)
        for label, family in self.families.iteritems():
            familyPath = os.path.join(path, label)
            if not os.path.exists(familyPath): os.mkdir(familyPath)
            family.save(familyPath)

    def saveLibraries(self, path):
        """
        Save the kinetics libraries to the given `path` on disk, where `path`
        points to the top-level folder of the kinetics libraries.
        """
        for label, library in self.libraries.iteritems():
            folders = label.split(os.sep)
            try:
                os.makedirs(os.path.join(path, *folders))
            except OSError:
                pass
            library.save(os.path.join(path, label, 'reactions.py'))
            library.saveDictionary(os.path.join(path, label, 'dictionary.txt'))

    def loadOld(self, path):
        """
        Load the old RMG kinetics database from the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        self.families = {}
        self.libraries = {}
        
        librariesPath = os.path.join(path, 'kinetics_libraries')
        for (root, dirs, files) in os.walk(os.path.join(path, 'kinetics_libraries')):
            if os.path.exists(os.path.join(root, 'species.txt')) and os.path.exists(os.path.join(root, 'reactions.txt')):
                library = KineticsLibrary(label=root[len(librariesPath)+1:], name=root[len(librariesPath)+1:])
                logging.warning("Loading {0}".format(root))
                library.loadOld(root)
                self.libraries[library.label] = library
                
        for (root, dirs, files) in os.walk(os.path.join(path, 'kinetics_groups')):
            if os.path.exists(os.path.join(root, 'dictionary.txt')) and os.path.exists(os.path.join(root, 'rateLibrary.txt')):
                label = os.path.split(root)[1]
                family = KineticsFamily(label=label)
                logging.warning("Loading {0}".format(root))
                family.loadOld(root)
                self.families[family.label] = family

        return self

    def saveOld(self, path):
        """
        Save the old RMG kinetics database to the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        librariesPath = os.path.join(path, 'kinetics_libraries')
        if not os.path.exists(librariesPath): os.mkdir(librariesPath)
        for library in self.libraries.values():
            libraryPath = os.path.join(librariesPath, library.label)
            library.saveOld(libraryPath)

        groupsPath = os.path.join(path, 'kinetics_groups')
        if not os.path.exists(groupsPath): os.mkdir(groupsPath)
        for label, family in self.families.iteritems():
            groupPath = os.path.join(groupsPath, label)
            family.saveOld(groupPath)
            
        with open(os.path.join(path,'kinetics_groups','families.txt'),'w') as f:
            f.write("""
////////////////////////////////////////////////////////////////////////////////
//
// REACTION FAMILIES USED BY RMG 
//
// Notes:
//
// The name of each forward reaction should match up with a folder
// in this directory containing the dictionary, tree, library, and adjusts
// for that family.
//
// Families can be deactivated by simply changing the "on/off" column to off.
//
// This file was generated by exporting the entirety of RMG-Py,
// so the on/off decisions come from the defaults there.
//
////////////////////////////////////////////////////////////////////////////////

// No.  on/off  Forward reaction
""")
            for number, label in enumerate(sorted(self.families.keys())):
                onoff = 'on ' if self.recommendedFamilies[label] else 'off'
                f.write("{num:<2d}    {onoff}     {label}\n".format(num=number, label=label, onoff=onoff))
    
    def generateReactions(self, reactants, products=None):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be :class:`Molecule` objects. This method
        searches the depository, libraries, and groups, in that order.
        """
        reactionList = []
        reactionList.extend(self.generateReactionsFromLibraries(reactants, products))
        reactionList.extend(self.generateReactionsFromFamilies(reactants, products))
        return reactionList

    def generateReactionsFromLibraries(self, reactants, products):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be :class:`Molecule` objects. This method
        searches the depository.
        """
        reactionList = []
        for label, libraryType in self.libraryOrder:
            # Generate reactions from reaction libraries (no need to generate them from seeds)
            if libraryType == "Reaction Library":
                reactionList.extend(self.generateReactionsFromLibrary(reactants, products, self.libraries[label]))
        return reactionList

    def generateReactionsFromLibrary(self, reactants, products, library):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be :class:`Molecule` objects. This method
        searches the depository.
        """
        reactionList = []
        for entry in library.entries.values():
            if entry.item.matchesMolecules(reactants):
                reaction = LibraryReaction(
                    reactants = entry.item.reactants[:],
                    products = entry.item.products[:],
                    degeneracy = entry.item.degeneracy,
                    reversible = entry.item.reversible,
                    duplicate = entry.item.duplicate,
                    kinetics = deepcopy(entry.data),
                    library = library,
                    entry = entry,
                )
                reactionList.append(reaction)
        if products:
            reactionList = filterReactions(reactants, products, reactionList)
        return reactionList

    def generateReactionsFromFamilies(self, reactants, products, only_families=None):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be :class:`Molecule` objects. This method
        applies the reaction family.
        If `only_families` is a list of strings, only families with those labels
        are used.
        """
        # If there are two structures and they are the same, then make a copy
        # of the second one so we can independently manipulate both of them 
        # This is for the case where A + A --> products
        if len(reactants) == 2 and reactants[0] == reactants[1]:
            reactants[1] = reactants[1].copy(deep=True)
        
        reactionList = []
        for label, family in self.families.iteritems():
            if only_families is None or label in only_families:
                reactionList.extend(family.generateReactions(reactants))
        if products:
            reactionList = filterReactions(reactants, products, reactionList)
        return reactionList

    def getForwardReactionForFamilyEntry(self, entry, family, thermoDatabase):
        """
        For a given `entry` for a reaction of the given reaction `family` (the
        string label of the family), return the reaction with kinetics and
        degeneracy for the "forward" direction as defined by the reaction 
        family. For families that are their own reverse, the direction the
        kinetics is given in will be preserved. If the entry contains 
        functional groups for the reactants, assume that it is given in the 
        forward direction and do nothing. Returns the reaction in the direction
        consistent with the reaction family template, and the matching template.
        Note that the returned reaction will have its kinetics and degeneracy
        set appropriately.
        
        In order to reverse the reactions that are given in the reverse of the
        direction the family is defined, we need to compute the thermodynamics
        of the reactants and products. For this reason you must also pass
        the `thermoDatabase` to use to generate the thermo data.
        """
        
        def matchSpeciesToMolecules(species, molecules):
            if len(species) == len(molecules) == 1:
                return species[0].isIsomorphic(molecules[0])
            elif len(species) == len(molecules) == 2:
                if species[0].isIsomorphic(molecules[0]) and species[1].isIsomorphic(molecules[1]):
                    return True
                elif species[0].isIsomorphic(molecules[1]) and species[1].isIsomorphic(molecules[0]):
                    return True
            return False

        reaction = None; template = None

        # Get the indicated reaction family
        try:
            groups = self.families[family].groups
        except KeyError:
            raise ValueError('Invalid value "{0}" for family parameter.'.format(family))

        if all([(isinstance(reactant, Group) or isinstance(reactant, LogicNode)) for reactant in entry.item.reactants]):
            # The entry is a rate rule, containing functional groups only
            # By convention, these are always given in the forward direction and
            # have kinetics defined on a per-site basis
            reaction = Reaction(
                reactants = entry.item.reactants[:],
                products = [],
                kinetics = entry.data,
                degeneracy = 1,
            )
            template = [groups.entries[label] for label in entry.label.split(';')]

        elif (all([isinstance(reactant, Molecule) for reactant in entry.item.reactants]) and
            all([isinstance(product, Molecule) for product in entry.item.products])):
            # The entry is a real reaction, containing molecules
            # These could be defined for either the forward or reverse direction
            # and could have a reaction-path degeneracy

            reaction = Reaction(reactants=[], products=[])
            for molecule in entry.item.reactants:
                reactant = Species(molecule=[molecule])
                reactant.generateResonanceIsomers()
                reactant.thermo = thermoDatabase.getThermoData(reactant)
                reaction.reactants.append(reactant)
            for molecule in entry.item.products:
                product = Species(molecule=[molecule])
                product.generateResonanceIsomers()
                product.thermo = thermoDatabase.getThermoData(product)
                reaction.products.append(product)

            # Generate all possible reactions involving the reactant species
            generatedReactions = self.generateReactionsFromFamilies([reactant.molecule for reactant in reaction.reactants], [], only_families=[family])

            # Remove from that set any reactions that don't produce the desired reactants and products
            forward = []; reverse = []
            for rxn in generatedReactions:
                if matchSpeciesToMolecules(reaction.reactants, rxn.reactants) and matchSpeciesToMolecules(reaction.products, rxn.products):
                    forward.append(rxn)
                if matchSpeciesToMolecules(reaction.reactants, rxn.products) and matchSpeciesToMolecules(reaction.products, rxn.reactants):
                    reverse.append(rxn)

            # We should now know whether the reaction is given in the forward or
            # reverse direction
            if len(forward) == 1 and len(reverse) == 0:
                # The reaction is in the forward direction, so use as-is
                reaction = forward[0]
                template = reaction.template
                # Don't forget to overwrite the estimated kinetics from the database with the kinetics for this entry
                reaction.kinetics = entry.data
            elif len(reverse) == 1 and len(forward) == 0:
                # The reaction is in the reverse direction
                # First fit Arrhenius kinetics in that direction
                Tdata = 1000.0 / numpy.arange(0.5, 3.301, 0.1, numpy.float64)
                kdata = numpy.zeros_like(Tdata)
                for i in range(Tdata.shape[0]):
                    kdata[i] = entry.data.getRateCoefficient(Tdata[i]) / reaction.getEquilibriumConstant(Tdata[i])
                kunits = 'm^3/(mol*s)' if len(reverse[0].reactants) == 2 else 's^-1'
                kinetics = Arrhenius().fitToData(Tdata, kdata, kunits, T0=1.0)
                kinetics.Tmin = entry.data.Tmin
                kinetics.Tmax = entry.data.Tmax
                kinetics.Pmin = entry.data.Pmin
                kinetics.Pmax = entry.data.Pmax
                # Now flip the direction
                reaction = reverse[0]
                reaction.kinetics = kinetics
                template = reaction.template
            elif len(reverse) > 0 and len(forward) > 0:
                print 'FAIL: Multiple reactions found for {0!r}.'.format(entry.label)
            elif len(reverse) == 0 and len(forward) == 0:
                print 'FAIL: No reactions found for "%s".' % (entry.label)
            else:
                print 'FAIL: Unable to estimate kinetics for {0!r}.'.format(entry.label)

        assert reaction is not None
        assert template is not None
        return reaction, template

