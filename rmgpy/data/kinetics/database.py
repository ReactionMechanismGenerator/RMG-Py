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
from rmgpy.data.base import LogicNode

from .family import  KineticsFamily
from .library import LibraryReaction, KineticsLibrary
from .common import ensure_species, generate_molecule_combos, \
                    find_degenerate_reactions, ensure_independent_atom_ids
from rmgpy.exceptions import DatabaseError

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
                    if library_name == "KlippensteinH2O2":
                        logging.info("""\n** Note: The KlippensteinH2O2 library was replaced and is no longer available in RMG.
For H2 combustion chemistry consider using either the BurkeH2inN2 or BurkeH2inArHe
library instead, depending on the main bath gas (N2 or Ar/He, respectively)\n""")
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
    
    def generate_reactions(self, reactants, products=None, only_families=None, resonance=True):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be :class:`Molecule` objects. This method
        searches the depository, libraries, and groups, in that order.
        """
        reactionList = []
        if only_families is None:
            reactionList.extend(self.generate_reactions_from_libraries(reactants, products))
        reactionList.extend(self.generate_reactions_from_families(reactants, products, only_families=None, resonance=True))
        return reactionList

    def generate_reactions_from_libraries(self, reactants, products=None):
        """
        Find all reactions from all loaded kinetics library involving the
        provided `reactants`, which can be either :class:`Molecule` objects or
        :class:`Species` objects.
        """
        reaction_list = []
        for label, library_type in self.libraryOrder:
            # Generate reactions from reaction libraries (no need to generate them from seeds)
            if library_type == "Reaction Library":
                reaction_list.extend(self.generate_reactions_from_library(self.libraries[label], reactants, products=products))
        return reaction_list

    def generate_reactions_from_library(self, library, reactants, products=None):
        """
        Find all reactions from the specified kinetics library involving the
        provided `reactants`, which can be either :class:`Molecule` objects or
        :class:`Species` objects.
        """
        reactants = ensure_species(reactants)

        reaction_list = []
        for entry in library.entries.values():
            if entry.item.matchesSpecies(reactants, products=products):
                reaction = LibraryReaction(
                    reactants = entry.item.reactants[:],
                    products = entry.item.products[:],
                    specificCollider = entry.item.specificCollider,
                    degeneracy = entry.item.degeneracy,
                    reversible = entry.item.reversible,
                    duplicate = entry.item.duplicate,
                    kinetics = deepcopy(entry.data),
                    library = library,
                    entry = entry,
                )
                reaction_list.append(reaction)

        return reaction_list

    def generate_reactions_from_families(self, reactants, products=None, only_families=None, resonance=True):
        """
        Generate all reactions between the provided list or tuple of one or two
        `reactants`, which can be either :class:`Molecule` objects or :class:`Species`
        objects. This method can apply all kinetics families or a selected subset.

        Args:
            reactants:      Molecules or Species to react
            products:       List of Molecules or Species of desired product structures (optional)
            only_families:  List of family labels to generate reactions from (optional)
                            Default is to generate reactions from all families
            resonance:      Flag to generate resonance structures for reactants and products (optional)
                            Default is True, resonance structures will be generated

        Returns:
            List of reactions containing Species objects with the specified reactants and products.
        """
        # Check if the reactants are the same
        # If they refer to the same memory address, then make a deep copy so
        # they can be manipulated independently
        same_reactants = False
        if len(reactants) == 2:
            if reactants[0] is reactants[1]:
                reactants[1] = reactants[1].copy(deep=True)
                same_reactants = True
            elif reactants[0].isIsomorphic(reactants[1]):
                same_reactants = True

        # Convert to Species objects if necessary
        reactants = ensure_species(reactants)

        # Label reactant atoms for proper degeneracy calculation
        ensure_independent_atom_ids(reactants, resonance=resonance)

        combos = generate_molecule_combos(reactants)

        reaction_list = []
        for combo in combos:
            reaction_list.extend(self.react_molecules(combo, products=products, only_families=only_families, prod_resonance=resonance))

        # Calculate reaction degeneracy
        reaction_list = find_degenerate_reactions(reaction_list, same_reactants, kinetics_database=self)
        # Add reverse attribute to families with ownReverse
        to_delete = []
        for i, rxn in enumerate(reaction_list):
            family = self.families[rxn.family]
            if family.ownReverse:
                successful = family.addReverseAttribute(rxn)
                if not successful:
                    to_delete.append(i)
        # Delete reactions which we could not find a reverse reaction for
        for i in reversed(to_delete):
            del reaction_list[i]

        return reaction_list

    def react_molecules(self, molecules, products=None, only_families=None, prod_resonance=True):
        """
        Generate reactions from all families for the input molecules.
        """
        reaction_list = []
        for label, family in self.families.iteritems():
            if only_families is None or label in only_families:
                reaction_list.extend(family.generateReactions(molecules, products=products, prod_resonance=prod_resonance))

        for reactant in molecules:
            reactant.clearLabeledAtoms()

        return reaction_list

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
                specificCollider = entry.item.specificCollider,
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
                reactant.generate_resonance_structures()
                reactant.thermo = thermoDatabase.getThermoData(reactant)
                reaction.reactants.append(reactant)
            for molecule in entry.item.products:
                product = Species(molecule=[molecule])
                product.generate_resonance_structures()
                product.thermo = thermoDatabase.getThermoData(product)
                reaction.products.append(product)

            # Generate all possible reactions involving the reactant species
            generatedReactions = self.generate_reactions_from_families([reactant.molecule for reactant in reaction.reactants], [], only_families=[family])

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
    
    def extractSourceFromComments(self, reaction):
        """
        `reaction`: A reaction object containing kinetics data and kinetics data comments.  
            Should be either a PDepReaction, LibraryReaction, or TemplateReaction object
            as loaded from the rmgpy.chemkin.loadChemkinFile function
        
        Parses the verbose string of comments from the thermo data of the species object,
        and extracts the thermo sources.

        Returns a dictionary with keys of either 'Rate Rules', 'Training', 'Library', or 'PDep'.
        A reaction can only be estimated using one of these methods.
        
        source = {'RateRules': (Family_Label, OriginalTemplate, RateRules),
                  'Library': String_Name_of_Library_Used,
                  'PDep': Network_Index,
                  'Training':  (Family_Label, Training_Reaction_Entry),
                  }
        """
        from rmgpy.rmg.pdep import PDepReaction
        from rmgpy.data.kinetics.library import LibraryReaction
        from rmgpy.data.kinetics.family import TemplateReaction
        
        source = {}
        
        if isinstance(reaction, TemplateReaction):
            # This reaction comes from rate rules
            training, dataSource = self.families[reaction.family].extractSourceFromComments(reaction)
            if training:
                source['Training'] = dataSource
            else:
                source['Rate Rules'] = dataSource
        elif isinstance(reaction, LibraryReaction):
            # This reaction comes from a reaction library or seed mechanism
            source['Library'] = reaction.library
            
        elif isinstance(reaction, PDepReaction):
            # This reaction is a pressure-dependent reaction
            source['PDep'] = reaction.network.index
        
        else:
            raise Exception('Reaction {} must be either a TemplateReaction, LibraryReaction, or PDepReaction object for source data to be extracted.'.format(reaction))
            
        return source
    
    def reconstructKineticsFromSource(self, reaction, source, fixBarrierHeight=False, forcePositiveBarrier=False):
        """
        Reaction is the original reaction with original kinetics.
        Note that for Library and PDep reactions this function does not do anything other than return the original kinetics...
        
        You must enter source data in the appropriate format such as returned from returned from self.extractSourceFromComments, 
        self-constructed.  
        fixBarrierHeight and forcePositiveBarrier will change the kinetics based on the Reaction.fixBarrierHeight function.
        Return Arrhenius form kinetics if the source is from training reaction or rate rules.
        """
        from rmgpy.data.thermo import findCp0andCpInf
        from rmgpy.thermo import Wilhoit
        if 'Library' in source:
            return reaction.kinetics
        elif 'PDep' in source:
            return reaction.kinetics
        else:
            rxnCopy = deepcopy(reaction)
            if 'Training' in source:
                trainingEntry = source['Training'][1]
                reverse = source['Training'][2]
                if reverse:
                    reverseKinetics = trainingEntry.data
                    rxnCopy.kinetics = reverseKinetics
                    forwardKinetics = rxnCopy.generateReverseRateCoefficient()
                    kinetics = forwardKinetics
                else:
                    kinetics = trainingEntry.data
            elif 'Rate Rules' in source:
    
                sourceDict = source['Rate Rules'][1]
                rules = sourceDict['rules']
                training = sourceDict['training']
                degeneracy = sourceDict['degeneracy']
    
                logA = 0
                n = 0
                alpha = 0
                E0 = 0
                for ruleEntry, weight in rules:
                    logA += numpy.log10(ruleEntry.data.A.value_si)*weight
                    n += ruleEntry.data.n.value_si*weight
                    alpha +=ruleEntry.data.alpha.value_si*weight
                    E0 +=ruleEntry.data.E0.value_si*weight
                for ruleEntry, trainingEntry, weight in training:
                    logA += numpy.log10(ruleEntry.data.A.value_si)*weight
                    n += ruleEntry.data.n.value_si*weight
                    alpha +=ruleEntry.data.alpha.value_si*weight
                    E0 +=ruleEntry.data.E0.value_si*weight
                
                Aunits = ruleEntry.data.A.units 
                if Aunits == 'cm^3/(mol*s)' or Aunits == 'cm^3/(molecule*s)' or Aunits == 'm^3/(molecule*s)':
                    Aunits = 'm^3/(mol*s)'
                elif Aunits == 'cm^6/(mol^2*s)' or Aunits == 'cm^6/(molecule^2*s)' or Aunits == 'm^6/(molecule^2*s)':
                    Aunits = 'm^6/(mol^2*s)'
                elif Aunits == 's^-1' or Aunits == 'm^3/(mol*s)' or Aunits == 'm^6/(mol^2*s)':
                    pass
                else:
                    raise Exception('Invalid units {0} for averaging kinetics.'.format(Aunits))
                kinetics = ArrheniusEP(
                    A = (degeneracy*10**logA, Aunits),
                    n = n,
                    alpha = alpha,
                    E0 = (E0*0.001,"kJ/mol"),
                )
            else:
                raise Exception("Source data must be either 'Library', 'PDep','Training', or 'Rate Rules'.")
                
            
            # Convert ArrheniusEP to Arrhenius
            if fixBarrierHeight:
                for spc in rxnCopy.reactants + rxnCopy.products:
                    # Need wilhoit to do this
                    if not isinstance(spc.thermo, Wilhoit):
                        findCp0andCpInf(spc, spc.thermo)
                        wilhoit = spc.thermo.toWilhoit()
                        spc.thermo = wilhoit
                        
                rxnCopy.kinetics = kinetics
                rxnCopy.fixBarrierHeight(forcePositive=forcePositiveBarrier)
                
                return rxnCopy.kinetics
            else:
                
                H298 = rxnCopy.getEnthalpyOfReaction(298)
                if isinstance(kinetics, ArrheniusEP):
                    kinetics = kinetics.toArrhenius(H298)
                return kinetics
