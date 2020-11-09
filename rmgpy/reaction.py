#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains classes and functions for working with chemical reactions.

From the `IUPAC Compendium of Chemical Terminology 
<http://dx.doi.org/10.1351/goldbook>`_, a chemical reaction is "a process that 
results in the interconversion of chemical species".

In RMG Py, a chemical reaction is represented in memory as a :class:`Reaction`
object. This module also provides the :class:`ReactionModel` class for
representing a set of chemical reactions and the species involved.
"""

import logging
import math
import os.path
from copy import deepcopy
from functools import reduce
from urllib.parse import quote

import cython
import numpy as np

import rmgpy.constants as constants
from rmgpy.exceptions import ReactionError, KineticsError
from rmgpy.kinetics import KineticsData, ArrheniusBM, ArrheniusEP, ThirdBody, Lindemann, Troe, Chebyshev, \
    PDepArrhenius, MultiArrhenius, MultiPDepArrhenius, get_rate_coefficient_units_from_reaction_order, \
    StickingCoefficient, SurfaceArrheniusBEP, StickingCoefficientBEP
from rmgpy.kinetics.arrhenius import Arrhenius  # Separate because we cimport from rmgpy.kinetics.arrhenius
from rmgpy.kinetics.surface import SurfaceArrhenius  # Separate because we cimport from rmgpy.kinetics.surface
from rmgpy.kinetics.diffusionLimited import diffusion_limiter
from rmgpy.molecule.element import Element, element_list
from rmgpy.molecule.molecule import Molecule, Atom
from rmgpy.pdep.reaction import calculate_microcanonical_rate_coefficient
from rmgpy.species import Species

################################################################################


class Reaction:
    """
    A chemical reaction. The attributes are:
    
    =================== =========================== ============================
    Attribute           Type                        Description
    =================== =========================== ============================
    `index`             :class:`int`                A unique nonnegative integer index
    `label`             ``str``                     A descriptive string label
    `reactants`         :class:`list`               The reactant species (as :class:`Species` objects)
    `products`          :class:`list`               The product species (as :class:`Species` objects)
    'specific_collider'  :class:`Species`            The collider species (as a :class:`Species` object)
    `kinetics`          :class:`KineticsModel`      The kinetics model to use for the reaction
    `network_kinetics`  :class:`Arrhenius`          The kinetics model to use for PDep network exploration if the `kinetics` attribute is :class:PDepKineticsModel:
    `reversible`        ``bool``                    ``True`` if the reaction is reversible, ``False`` if not
    `transition_state`   :class:`TransitionState`    The transition state
    `duplicate`         ``bool``                    ``True`` if the reaction is known to be a duplicate, ``False`` if not
    `degeneracy`        :class:`double`             The reaction path degeneracy for the reaction
    `pairs`             ``list``                    Reactant-product pairings to use in converting reaction flux to species flux
    `allow_pdep_route`  ``bool``                    ``True`` if the reaction has an additional PDep pathway, ``False`` if not (by default), used for LibraryReactions
    `elementary_high_p` ``bool``                    If ``True``, pressure dependent kinetics will be generated (relevant only for unimolecular library reactions)
                                                    If ``False`` (by default), this library reaction will not be explored.
                                                    Only unimolecular library reactions with high pressure limit kinetics should be flagged (not if the kinetics were measured at some relatively low pressure)
    `comment`           ``str``                     A description of the reaction source (optional)
    `is_forward`        ``bool``                    Indicates if the reaction was generated in the forward (true) or reverse (false)
    `rank`              ``int``                     Integer indicating the accuracy of the kinetics for this reaction
    =================== =========================== ============================
    
    """

    def __init__(self,
                 index=-1,
                 label='',
                 reactants=None,
                 products=None,
                 specific_collider=None,
                 kinetics=None,
                 network_kinetics=None,
                 reversible=True,
                 transition_state=None,
                 duplicate=False,
                 degeneracy=1,
                 pairs=None,
                 allow_pdep_route=False,
                 elementary_high_p=False,
                 allow_max_rate_violation=False,
                 rank=None,
                 comment='',
                 is_forward=None,
                 ):
        self.index = index
        self.label = label
        self.reactants = reactants
        self.products = products
        self.specific_collider = specific_collider
        self._degeneracy = degeneracy
        self.kinetics = kinetics
        self.network_kinetics = network_kinetics
        self.reversible = reversible
        self.transition_state = transition_state
        self.duplicate = duplicate
        self.pairs = pairs
        self.allow_pdep_route = allow_pdep_route
        self.elementary_high_p = elementary_high_p
        self.comment = comment
        self.k_effective_cache = {}
        self.is_forward = is_forward
        self.allow_max_rate_violation = allow_max_rate_violation
        self.rank = rank

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'Reaction('
        if self.index != -1: string += 'index={0:d}, '.format(self.index)
        if self.label != '': string += 'label={0!r}, '.format(self.label)
        if self.reactants is not None: string += 'reactants={0!r}, '.format(self.reactants)
        if self.products is not None: string += 'products={0!r}, '.format(self.products)
        if self.specific_collider is not None: string += 'specific_collider={0!r}, '.format(self.specific_collider)
        if self.kinetics is not None: string += 'kinetics={0!r}, '.format(self.kinetics)
        if self.network_kinetics is not None: string += 'network_kinetics={0!r}, '.format(self.network_kinetics)
        if not self.reversible: string += 'reversible={0}, '.format(self.reversible)
        if self.transition_state is not None: string += 'transition_state={0!r}, '.format(self.transition_state)
        if self.duplicate: string += 'duplicate={0}, '.format(self.duplicate)
        if self.degeneracy != 1: string += 'degeneracy={0:.1f}, '.format(self.degeneracy)
        if self.pairs is not None: string += 'pairs={0}, '.format(self.pairs)
        if self.allow_pdep_route: string += 'allow_pdep_route={0}, '.format(self.allow_pdep_route)
        if self.elementary_high_p: string += 'elementary_high_p={0}, '.format(self.elementary_high_p)
        if self.comment != '': string += 'comment={0!r}, '.format(self.comment)
        if self.rank is not None: string += 'rank={0!r},'.format(self.rank)
        string = string[:-2] + ')'
        return string

    def __str__(self):
        """
        Return a string representation of the reaction, in the form 'A + B <=> C + D'.
        If a specific_collider exists, the srting representation is 'A + B (+S) <=> C + D (+S)'.
        """
        return self.to_labeled_str(use_index=True)

    def to_labeled_str(self, use_index=False):
        """
        the same as __str__ except that the labels are assumed to exist and used for reactant and products rather than 
        the labels plus the index in parentheses
        """
        arrow = ' <=> ' if self.reversible else ' => '
        reactants = ' + '.join([str(s) if use_index else s.label for s in self.reactants])
        products = ' + '.join([str(s) if use_index else s.label for s in self.products])
        collider = ' (+{0!s})'.format(self.specific_collider) if self.specific_collider else ''

        return reactants + collider + arrow + products + collider

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Reaction, (self.index,
                           self.label,
                           self.reactants,
                           self.products,
                           self.specific_collider,
                           self.kinetics,
                           self.network_kinetics,
                           self.reversible,
                           self.transition_state,
                           self.duplicate,
                           self.degeneracy,
                           self.pairs,
                           self.allow_pdep_route,
                           self.elementary_high_p,
                           self.rank,
                           self.comment
                           ))

    @property
    def degeneracy(self):
        """
        The reaction path degeneracy for this reaction.

        If the reaction has kinetics, changing the degeneracy
        will adjust the reaction rate by a ratio of the new
        degeneracy to the old degeneracy.
        """
        return self._degeneracy

    @degeneracy.setter
    def degeneracy(self, new):
        # modify rate if kinetics exists
        if self.kinetics is not None:
            if self._degeneracy < 2:
                degeneracy_ratio = new
            else:
                degeneracy_ratio = (new * 1.0) / self._degeneracy
            # fix kinetics comment with new degeneracy
            if 'Multiplied by reaction path degeneracy {}'.format(self._degeneracy) in self.kinetics.comment:
                self.kinetics.comment = self.kinetics.comment.replace(
                    'Multiplied by reaction path degeneracy {}'.format(self._degeneracy),
                    'Multiplied by reaction path degeneracy {}'.format(float(new)))
            elif self.kinetics.comment:
                self.kinetics.comment += 'Multiplied by reaction path degeneracy {}'.format(float(new))
            self.kinetics.change_rate(degeneracy_ratio)
        # set new degeneracy
        self._degeneracy = new

    def to_chemkin(self, species_list=None, kinetics=True):
        """
        Return the chemkin-formatted string for this reaction.
        
        If `kinetics` is set to True, the chemkin format kinetics will also
        be returned (requires the `species_list` to figure out third body colliders.)
        Otherwise, only the reaction string will be returned.
        """
        import rmgpy.chemkin
        if kinetics:
            return rmgpy.chemkin.write_kinetics_entry(self, species_list)
        else:
            return rmgpy.chemkin.write_reaction_string(self)

    def to_cantera(self, species_list=None, use_chemkin_identifier=False):
        """
        Converts the RMG Reaction object to a Cantera Reaction object
        with the appropriate reaction class.

        If use_chemkin_identifier is set to False, the species label is used
        instead. Be sure that species' labels are unique when setting it False.
        """
        import cantera as ct

        if species_list is None:
            species_list = []

        # Create the dictionaries containing species strings and their stoichiometries
        # for initializing the cantera reaction object
        ct_reactants = {}
        ct_collider = {}
        for reactant in self.reactants:
            if use_chemkin_identifier:
                reactant_name = reactant.to_chemkin()
            else:
                reactant_name = reactant.label
            if reactant_name in ct_reactants:
                ct_reactants[reactant_name] += 1
            else:
                ct_reactants[reactant_name] = 1
        ct_products = {}
        for product in self.products:
            if use_chemkin_identifier:
                product_name = product.to_chemkin()
            else:
                product_name = product.label
            if product_name in ct_products:
                ct_products[product_name] += 1
            else:
                ct_products[product_name] = 1
        if self.specific_collider:  # add a specific collider if exists
            ct_collider[self.specific_collider.to_chemkin() if use_chemkin_identifier else self.specific_collider.label] = 1

        if self.kinetics:
            if isinstance(self.kinetics, Arrhenius):
                # Create an Elementary Reaction
                ct_reaction = ct.ElementaryReaction(reactants=ct_reactants, products=ct_products)
            elif isinstance(self.kinetics, MultiArrhenius):
                # Return a list of elementary reactions which are duplicates
                ct_reaction = [ct.ElementaryReaction(reactants=ct_reactants, products=ct_products)
                               for arr in self.kinetics.arrhenius]

            elif isinstance(self.kinetics, PDepArrhenius):
                ct_reaction = ct.PlogReaction(reactants=ct_reactants, products=ct_products)

            elif isinstance(self.kinetics, MultiPDepArrhenius):
                ct_reaction = [ct.PlogReaction(reactants=ct_reactants, products=ct_products)
                               for arr in self.kinetics.arrhenius]

            elif isinstance(self.kinetics, Chebyshev):
                ct_reaction = ct.ChebyshevReaction(reactants=ct_reactants, products=ct_products)

            elif isinstance(self.kinetics, ThirdBody):
                if ct_collider is not None:
                    ct_reaction = ct.ThreeBodyReaction(reactants=ct_reactants, products=ct_products, tbody=ct_collider)
                else:
                    ct_reaction = ct.ThreeBodyReaction(reactants=ct_reactants, products=ct_products)

            elif isinstance(self.kinetics, Lindemann) or isinstance(self.kinetics, Troe):
                if ct_collider is not None:
                    ct_reaction = ct.FalloffReaction(reactants=ct_reactants, products=ct_products, tbody=ct_collider)
                else:
                    ct_reaction = ct.FalloffReaction(reactants=ct_reactants, products=ct_products)
            else:
                raise NotImplementedError('Unable to set cantera kinetics for {0}'.format(self.kinetics))

            # Set reversibility, duplicate, and ID attributes
            if isinstance(ct_reaction, list):
                for rxn in ct_reaction:
                    rxn.reversible = self.reversible
                    # Set the duplicate flag to true since this reaction comes from multiarrhenius or multipdeparrhenius 
                    rxn.duplicate = True
                    # Set the ID flag to the original rmg index 
                    rxn.ID = str(self.index)
            else:
                ct_reaction.reversible = self.reversible
                ct_reaction.duplicate = self.duplicate
                ct_reaction.ID = str(self.index)

            self.kinetics.set_cantera_kinetics(ct_reaction, species_list)

            return ct_reaction

        else:
            raise Exception('Cantera reaction cannot be created because there was no kinetics.')

    def get_url(self):
        """
        Get a URL to search for this reaction in the rmg website.
        """
        # eg. http://dev.rmg.mit.edu/database/kinetics/reaction/reactant1=1%20C%200%20%7B2,S%7D;2%20O%200%20%7B1,S%7D;__reactant2=1%20C%202T;__product1=1%20C%201;__product2=1%20C%200%20%7B2,S%7D;2%20O%201%20%7B1,S%7D;

        base_url = "http://rmg.mit.edu/database/kinetics/reaction/"

        rxn_string = ''
        for i, species in enumerate(self.reactants):
            adjlist = species.molecule[0].to_adjacency_list(remove_h=False)
            rxn_string += "reactant{0}={1}__".format(i + 1, adjlist)
        for i, species in enumerate(self.products):
            adjlist = species.molecule[0].to_adjacency_list(remove_h=False)
            rxn_string += "product{0}={1}__".format(i + 1, adjlist)

        url = base_url + quote(rxn_string)
        return url.strip('_')

    def is_isomerization(self):
        """
        Return ``True`` if the reaction represents an isomerization reaction
        :math:`\\ce{A <=> B}` or ``False`` if not.
        """
        return len(self.reactants) == 1 and len(self.products) == 1

    def is_association(self):
        """
        Return ``True`` if the reaction represents an association reaction
        :math:`\\ce{A + B <=> C}` or ``False`` if not.
        """
        return len(self.reactants) > 1 and len(self.products) == 1

    def is_dissociation(self):
        """
        Return ``True`` if the reaction represents a dissociation reaction
        :math:`\\ce{A <=> B + C}` or ``False`` if not.
        """
        return len(self.reactants) == 1 and len(self.products) > 1

    def is_unimolecular(self):
        """
        Return ``True`` if the reaction has a single molecule as either reactant or product (or both)
        :math:`\\ce{A <=> B + C}` or :math:`\\ce{A + B <=> C}` or :math:`\\ce{A <=> B}`,
        or ``False`` if not.
        """
        return len(self.reactants) == 1 or len(self.products) == 1

    def is_surface_reaction(self):
        """
        Return ``True`` if one or more reactants or products are surface species (or surface sites)
        """
        for spec in self.reactants:
            if spec.contains_surface_site():
                return True
        for spec in self.products:
            if spec.contains_surface_site():
                return True
        return False

    def has_template(self, reactants, products):
        """
        Return ``True`` if the reaction matches the template of `reactants`
        and `products`, which are both lists of :class:`Species` objects, or
        ``False`` if not.
        """
        return ((all([spec in self.reactants for spec in reactants]) and
                 all([spec in self.products for spec in products])) or
                (all([spec in self.products for spec in reactants]) and
                 all([spec in self.reactants for spec in products])))

    def matches_species(self, reactants, products=None):
        """
        Compares the provided reactants and products against the reactants
        and products of this reaction. Both directions are checked.

        Args:
            reactants (list): Species required on one side of the reaction
            products (list, optional): Species required on the other side
        """
        # Check forward direction
        if same_species_lists(self.reactants, reactants):
            if products is None or same_species_lists(self.products, products):
                return True
            else:
                return False
        elif same_species_lists(self.products, reactants):
            if products is None or same_species_lists(self.reactants, products):
                return True
            else:
                return False
        else:
            return False

    def is_isomorphic(self, other, either_direction=True, check_identical=False, check_only_label=False,
                      check_template_rxn_products=False, generate_initial_map=False, strict=True, save_order=False):
        """
        Return ``True`` if this reaction is the same as the `other` reaction,
        or ``False`` if they are different. The comparison involves comparing
        isomorphism of reactants and products, and doesn't use any kinetic
        information.

        Args:
            either_direction (bool, optional):            if ``False``,then the reaction direction must match.
            check_identical (bool, optional):             if ``True``, check that atom ID's match (used for checking degeneracy)
            check_only_label (bool, optional):            if ``True``, only check the string representation,
                                                          ignoring molecular structure comparisons
            check_template_rxn_products (bool, optional): if ``True``, only check isomorphism of reaction products
                                                          (used when we know the reactants are identical, i.e. in generating reactions)
            generate_initial_map (bool, optional):        if ``True``, initialize map by pairing atoms with same labels
            strict (bool, optional):                      if ``False``, perform isomorphism ignoring electrons
            save_order (bool, optional):                  if ``True``, perform isomorphism saving atom order
        """
        if check_template_rxn_products:
            try:
                species1 = self.products if self.is_forward else self.reactants
                species2 = other.products if other.is_forward else other.reactants
            except AttributeError:
                raise TypeError('Only use check_template_rxn_products flag for TemplateReactions.')

            return same_species_lists(species1, species2,
                                      check_identical=check_identical,
                                      only_check_label=check_only_label,
                                      generate_initial_map=generate_initial_map,
                                      strict=strict,
                                      save_order=save_order)

        # Compare reactants to reactants
        forward_reactants_match = same_species_lists(self.reactants, other.reactants,
                                                     check_identical=check_identical,
                                                     only_check_label=check_only_label,
                                                     generate_initial_map=generate_initial_map,
                                                     strict=strict,
                                                     save_order=save_order)

        # Compare products to products
        forward_products_match = same_species_lists(self.products, other.products,
                                                    check_identical=check_identical,
                                                    only_check_label=check_only_label,
                                                    generate_initial_map=generate_initial_map,
                                                    strict=strict,
                                                    save_order=save_order)

        # Compare specific_collider to specific_collider
        collider_match = (self.specific_collider == other.specific_collider)

        # Return now, if we can
        if forward_reactants_match and forward_products_match and collider_match:
            return True
        if not either_direction:
            return False

        # Compare reactants to products
        reverse_reactants_match = same_species_lists(self.reactants, other.products,
                                                     check_identical=check_identical,
                                                     only_check_label=check_only_label,
                                                     generate_initial_map=generate_initial_map,
                                                     strict=strict,
                                                     save_order=save_order)

        # Compare products to reactants
        reverse_products_match = same_species_lists(self.products, other.reactants,
                                                    check_identical=check_identical,
                                                    only_check_label=check_only_label,
                                                    generate_initial_map=generate_initial_map,
                                                    strict=strict,
                                                    save_order=save_order)

        # should have already returned if it matches forwards, or we're not allowed to match backwards
        return reverse_reactants_match and reverse_products_match and collider_match

    def get_enthalpy_of_reaction(self, T):
        """
        Return the enthalpy of reaction in J/mol evaluated at temperature
        `T` in K.
        """
        cython.declare(dHrxn=cython.double, reactant=Species, product=Species)
        dHrxn = 0.0
        for reactant in self.reactants:
            dHrxn -= reactant.get_enthalpy(T)
        for product in self.products:
            dHrxn += product.get_enthalpy(T)
        return dHrxn

    def get_entropy_of_reaction(self, T):
        """
        Return the entropy of reaction in J/mol*K evaluated at temperature `T`
        in K.
        """
        cython.declare(dSrxn=cython.double, reactant=Species, product=Species)
        dSrxn = 0.0
        for reactant in self.reactants:
            dSrxn -= reactant.get_entropy(T)
        for product in self.products:
            dSrxn += product.get_entropy(T)
        return dSrxn

    def get_free_energy_of_reaction(self, T):
        """
        Return the Gibbs free energy of reaction in J/mol evaluated at
        temperature `T` in K.
        """
        cython.declare(dGrxn=cython.double, reactant=Species, product=Species)
        dGrxn = 0.0
        for reactant in self.reactants:
            try:
                dGrxn -= reactant.get_free_energy(T)
            except Exception:
                logging.error("Problem with reactant {!r} in reaction {!s}".format(reactant, self))
                raise
        for product in self.products:
            try:
                dGrxn += product.get_free_energy(T)
            except Exception:
                logging.error("Problem with product {!r} in reaction {!s}".format(reactant, self))
                raise
        return dGrxn

    def get_equilibrium_constant(self, T, type='Kc', surface_site_density=2.5e-05):
        """
        Return the equilibrium constant for the reaction at the specified
        temperature `T` in K and reference `surface_site_density`
        in mol/m^2 (2.5e-05 default) The `type` parameter lets you specify
        the quantities used in the equilibrium constant: ``Ka`` for activities,
        ``Kc`` for concentrations (default), or ``Kp`` for pressures.  This
        function assumes a reference pressure of 1e5 Pa for gas phases species
        and uses the ideal gas law to determine reference concentrations. For
        surface species, the `surface_site_density` is the assumed reference.
        """
        cython.declare(dGrxn=cython.double, K=cython.double, C0=cython.double, P0=cython.double)
        # Use free energy of reaction to calculate Ka
        dGrxn = self.get_free_energy_of_reaction(T)
        K = np.exp(-dGrxn / constants.R / T)
        # Convert Ka to Kc or Kp if specified
        # Assume a pressure of 1e5 Pa for gas phase species
        P0 = 1e5
        # Determine the number of gas phase reactants and products. For gas species,
        # we will use 1e5 Pa and ideal gas law to determine reference concentration.
        try:
            number_of_gas_reactants = len([spcs for spcs in self.reactants if not spcs.contains_surface_site()])
            number_of_gas_products = len([spcs for spcs in self.products if not spcs.contains_surface_site()])
        except IndexError:
            logging.warning("Species do not have an rmgpy.molecule.Molecule "
                            "Cannot determine phases of species. We will assume "
                            "ideal gas mixture when calculating Kc and Kp.")
            number_of_gas_reactants = len(self.reactants)
            number_of_gas_products = len(self.products)

        # Determine the number of surface reactants and products.  For surface species,
        # we will use the provided `surface_site_density` as the reference
        number_of_surface_reactants = len(self.reactants) - number_of_gas_reactants
        number_of_surface_products = len(self.products) - number_of_gas_products

        # Determine the change in the number of mols of gas and surface species in the reaction
        dN_surf = number_of_surface_products - number_of_surface_reactants # change in mols of surface spcs
        dN_gas = number_of_gas_products - number_of_gas_reactants # change in mols of gas spcs

        if type == 'Kc':
            # Convert from Ka to Kc; C0 is the reference concentration
            if dN_gas:
                C0 = P0 / constants.R / T
                K *= C0 ** dN_gas
            if dN_surf:
                K *= surface_site_density ** dN_surf
        elif type == 'Kp':
            # Convert from Ka to Kp; P0 is the reference pressure
            K *= P0 ** dN_gas
        elif type != 'Ka' and type != '':
            raise ReactionError('Invalid type "{0}" passed to Reaction.get_equilibrium_constant(); '
                                'should be "Ka", "Kc", or "Kp".'.format(type))
        if K == 0:
            raise ReactionError('Got equilibrium constant of 0')
        return K

    def get_enthalpies_of_reaction(self, Tlist):
        """
        Return the enthalpies of reaction in J/mol evaluated at temperatures
        `Tlist` in K.
        """
        return np.array([self.get_enthalpy_of_reaction(T) for T in Tlist], np.float64)

    def get_entropies_of_reaction(self, Tlist):
        """
        Return the entropies of reaction in J/mol*K evaluated at temperatures
        `Tlist` in K.
        """
        return np.array([self.get_entropy_of_reaction(T) for T in Tlist], np.float64)

    def get_free_energies_of_reaction(self, Tlist):
        """
        Return the Gibbs free energies of reaction in J/mol evaluated at
        temperatures `Tlist` in K.
        """
        return np.array([self.get_free_energy_of_reaction(T) for T in Tlist], np.float64)

    def get_equilibrium_constants(self, Tlist, type='Kc'):
        """
        Return the equilibrium constants for the reaction at the specified
        temperatures `Tlist` in K. The `type` parameter lets you specify the
        quantities used in the equilibrium constant: ``Ka`` for	activities,
        ``Kc`` for concentrations (default), or ``Kp`` for pressures. Note that
        this function currently assumes an ideal gas mixture.
        """
        return np.array([self.get_equilibrium_constant(T, type) for T in Tlist], np.float64)

    def get_stoichiometric_coefficient(self, spec):
        """
        Return the stoichiometric coefficient of species `spec` in the reaction.
        The stoichiometric coefficient is increased by one for each time `spec`
        appears as a product and decreased by one for each time `spec` appears
        as a reactant.
        """
        cython.declare(stoich=cython.int, reactant=Species, product=Species)
        stoich = 0
        for reactant in self.reactants:
            if reactant is spec: stoich -= 1
        for product in self.products:
            if product is spec: stoich += 1
        return stoich

    def get_rate_coefficient(self, T, P=0):
        """
        Return the overall rate coefficient for the forward reaction at
        temperature `T` in K and pressure `P` in Pa, including any reaction
        path degeneracies.
        
        If diffusion_limiter is enabled, the reaction is in the liquid phase and we use
        a diffusion limitation to correct the rate. If not, then use the intrinsic rate
        coefficient.
        """
        if diffusion_limiter.enabled:
            try:
                k = self.k_effective_cache[T]
            except KeyError:
                k = diffusion_limiter.get_effective_rate(self, T)
                self.k_effective_cache[T] = k
            return k
        else:
            return self.kinetics.get_rate_coefficient(T, P)

    def get_surface_rate_coefficient(self, T, surface_site_density):
        """
        Return the overall surface rate coefficient for the forward reaction at
        temperature `T` in K with surface site density `surface_site_density` in mol/m2.
        Value is returned in combination of [m,mol,s]
        """
        cython.declare(rateCoefficient=cython.double,
                       molecularWeight_kg=cython.double, )

        if diffusion_limiter.enabled:
            raise NotImplementedError()
        if not self.is_surface_reaction():
            raise ReactionError("This is not a surface reaction!")

        if isinstance(self.kinetics, StickingCoefficient):
            rate_coefficient = self.kinetics.get_sticking_coefficient(T)
            adsorbate = None
            for r in self.reactants:
                if r.contains_surface_site():
                    rate_coefficient /= surface_site_density
                else:
                    if adsorbate is None:
                        adsorbate = r
                    else:
                        logging.error("Error in kinetics for reaction {0!s}: "
                                      "more than one adsorbate detected".format(self))
                        raise ReactionError("More than one adsorbate detected")

            if adsorbate is None or adsorbate.contains_surface_site():
                logging.error("Problem reaction: {0!s}".format(self))
                raise ReactionError("Couldn't find the adsorbate!")
            molecular_weight_kg = adsorbate.molecular_weight.value_si
            # molecular_weight_kg in kg per molecule
            rate_coefficient *= math.sqrt(constants.kB * T / (2 * math.pi * molecular_weight_kg))

            # ToDo: missing the sigma terms for bidentate species. only works for single site adsorption
            return rate_coefficient

        if isinstance(self.kinetics, SurfaceArrhenius):
            return self.kinetics.get_rate_coefficient(T, P=0)

        raise NotImplementedError("Can't get_surface_rate_coefficient for kinetics type {!r}".format(type(self.kinetics)))

    def fix_diffusion_limited_a_factor(self, T):
        """
        Decrease the pre-exponential factor (A) by the diffusion factor
        to account for the diffusion limit at the specified temperature.
        """
        if not diffusion_limiter.enabled:
            return
        # Obtain effective rate
        try:
            k = self.k_effective_cache[T]
        except KeyError:
            k = diffusion_limiter.get_effective_rate(self, T)
            self.k_effective_cache[T] = k

        # calculate diffusion factor
        diffusion_factor = k / self.kinetics.get_rate_coefficient(T, P=0)
        # update preexponential factor
        self.kinetics.A = self.kinetics.A * diffusion_factor
        # Add a comment to self.kinetics.comment
        self.kinetics.comment.append(
            ("Pre-exponential factor A has been decreased by the "
             "diffusion factor {0.2g} evaluated at {1} K.").format(
                diffusion_factor, T))

    def fix_barrier_height(self, force_positive=False):
        """
        Turns the kinetics into Arrhenius (if they were ArrheniusEP)
        and ensures the activation energy is at least the endothermicity
        for endothermic reactions, and is not negative only as a result 
        of using Evans Polanyi with an exothermic reaction.
        If `force_positive` is True, then all reactions
        are forced to have a non-negative barrier.
        """
        cython.declare(H0=cython.double, H298=cython.double, Ea=cython.double)

        if self.kinetics is None:
            raise KineticsError("Cannot fix barrier height for reactions with no kinetics attribute")

        H298 = self.get_enthalpy_of_reaction(298)
        H0 = sum([spec.get_thermo_data().E0.value_si for spec in self.products]) \
             - sum([spec.get_thermo_data().E0.value_si for spec in self.reactants])
        if isinstance(self.kinetics, (ArrheniusEP, SurfaceArrheniusBEP, StickingCoefficientBEP, ArrheniusBM)):
            Ea = self.kinetics.E0.value_si  # temporarily using Ea to store the intrinsic barrier height E0
            self.kinetics = self.kinetics.to_arrhenius(H298)
            if self.kinetics.Ea.value_si < 0.0 and self.kinetics.Ea.value_si < Ea:
                # Calculated Ea (from Evans-Polanyi) is negative AND below than the intrinsic E0
                Ea = min(0.0, Ea)  # (the lowest we want it to be)
                self.kinetics.comment += "\nEa raised from {0:.1f} to {1:.1f} kJ/mol.".format(
                    self.kinetics.Ea.value_si / 1000., Ea / 1000.)
                logging.info("For reaction {0!s} Ea raised from {1:.1f} to {2:.1f} kJ/mol.".format(
                    self, self.kinetics.Ea.value_si / 1000., Ea / 1000.))
                self.kinetics.Ea.value_si = Ea
        if isinstance(self.kinetics, (Arrhenius, StickingCoefficient)):  # SurfaceArrhenius is a subclass of Arrhenius
            Ea = self.kinetics.Ea.value_si
            if H0 >= 0 and Ea < H0:
                self.kinetics.Ea.value_si = H0
                self.kinetics.comment += "\nEa raised from {0:.1f} to {1:.1f} kJ/mol to match endothermicity of " \
                                         "reaction.".format( Ea / 1000., H0 / 1000.)
                logging.info("For reaction {2!s}, Ea raised from {0:.1f} to {1:.1f} kJ/mol to match "
                             "endothermicity of reaction.".format( Ea / 1000., H0 / 1000., self))
        if force_positive and isinstance(self.kinetics, (Arrhenius, StickingCoefficient)) and self.kinetics.Ea.value_si < 0:
            self.kinetics.comment += "\nEa raised from {0:.1f} to 0 kJ/mol.".format(self.kinetics.Ea.value_si / 1000.)
            logging.info("For reaction {1!s} Ea raised from {0:.1f} to 0 kJ/mol.".format(
                self.kinetics.Ea.value_si / 1000., self))
            self.kinetics.Ea.value_si = 0
        if self.kinetics.is_pressure_dependent() and self.network_kinetics is not None:
            Ea = self.network_kinetics.Ea.value_si
            if H0 >= 0 and Ea < H0:
                self.network_kinetics.Ea.value_si = H0
                self.network_kinetics.comment += "\nEa raised from {0:.1f} to {1:.1f} kJ/mol to match endothermicity of" \
                                                 " reaction.".format(Ea / 1000., H0 / 1000.)
                logging.info("For reaction {2!s}, Ea of the high pressure limit kinetics raised from {0:.1f} to {1:.1f}"
                             " kJ/mol to match endothermicity of reaction.".format(Ea / 1000., H0 / 1000., self))
            if force_positive and isinstance(self.kinetics, Arrhenius) and self.kinetics.Ea.value_si < 0:
                self.network_kinetics.comment += "\nEa raised from {0:.1f} to 0 kJ/mol.".format(
                    self.kinetics.Ea.value_si / 1000.)
                logging.info("For reaction {1!s} Ea of the high pressure limit kinetics raised from {0:.1f} to 0"
                             " kJ/mol.".format(self.kinetics.Ea.value_si / 1000., self))
                self.kinetics.Ea.value_si = 0

    def reverse_arrhenius_rate(self, k_forward, reverse_units, Tmin=None, Tmax=None):
        """
        Reverses the given k_forward, which must be an Arrhenius type.
        You must supply the correct units for the reverse rate.
        The equilibrium constant is evaluated from the current reaction instance (self).
        """
        cython.declare(kf=Arrhenius, kr=Arrhenius)
        cython.declare(Tlist=np.ndarray, klist=np.ndarray, i=cython.int)
        kf = k_forward
        assert isinstance(kf, Arrhenius), "Only reverses Arrhenius rates"
        if Tmin is not None and Tmax is not None:
            Tlist = 1.0 / np.linspace(1.0 / Tmax.value, 1.0 / Tmin.value, 50)
        else:
            Tlist = 1.0 / np.arange(0.0005, 0.0034, 0.0001)
        # Determine the values of the reverse rate coefficient k_r(T) at each temperature
        klist = np.zeros_like(Tlist)
        for i in range(len(Tlist)):
            klist[i] = kf.get_rate_coefficient(Tlist[i]) / self.get_equilibrium_constant(Tlist[i])
        kr = Arrhenius()
        kr.fit_to_data(Tlist, klist, reverse_units, kf.T0.value_si)
        return kr

    def reverse_surface_arrhenius_rate(self, k_forward, reverse_units, Tmin=None, Tmax=None):
        """
        Reverses the given k_forward, which must be a SurfaceArrhenius type.
        You must supply the correct units for the reverse rate.
        The equilibrium constant is evaluated from the current reaction instance (self).
        """
        cython.declare(kf=SurfaceArrhenius, kr=SurfaceArrhenius)
        cython.declare(Tlist=np.ndarray, klist=np.ndarray, i=cython.int)
        kf = k_forward
        if not isinstance(kf, SurfaceArrhenius): # Only reverse SurfaceArrhenius rates
            raise TypeError(f'Expected a SurfaceArrhenius object for k_forward but received {kf}')
        if Tmin is not None and Tmax is not None:
            Tlist = 1.0 / np.linspace(1.0 / Tmax.value, 1.0 / Tmin.value, 50)
        else:
            Tlist = 1.0 / np.arange(0.0005, 0.0034, 0.0001)
        # Determine the values of the reverse rate coefficient k_r(T) at each temperature
        klist = np.zeros_like(Tlist)
        for i in range(len(Tlist)):
            klist[i] = kf.get_rate_coefficient(Tlist[i]) / self.get_equilibrium_constant(Tlist[i])
        kr = SurfaceArrhenius()
        kr.fit_to_data(Tlist, klist, reverse_units, kf.T0.value_si)
        return kr

    def generate_reverse_rate_coefficient(self, network_kinetics=False, Tmin=None, Tmax=None):
        """
        Generate and return a rate coefficient model for the reverse reaction. 
        Currently this only works if the `kinetics` attribute is one of several
        (but not necessarily all) kinetics types.
        """
        cython.declare(Tlist=np.ndarray, Plist=np.ndarray, K=np.ndarray,
                       rxn=Reaction, klist=np.ndarray, i=cython.size_t,
                       Tindex=cython.size_t, Pindex=cython.size_t)

        supported_types = (
            KineticsData.__name__,
            Arrhenius.__name__,
            SurfaceArrhenius.__name__,
            MultiArrhenius.__name__,
            PDepArrhenius.__name__,
            MultiPDepArrhenius.__name__,
            Chebyshev.__name__,
            ThirdBody.__name__,
            Lindemann.__name__,
            Troe.__name__,
        )

        # Get the units for the reverse rate coefficient
        try:
            surf_prods = [spcs for spcs in self.products if spcs.contains_surface_site()]
        except IndexError:
            surf_prods = []
            logging.warning(f"Species do not have an rmgpy.molecule.Molecule "  
                            "Cannot determine phases of species. We will assume gas"
                            )
        n_surf = len(surf_prods)
        n_gas = len(self.products) - len(surf_prods)
        kunits = get_rate_coefficient_units_from_reaction_order(n_gas, n_surf)

        kf = self.kinetics
        if isinstance(kf, KineticsData):

            Tlist = kf.Tdata.value_si
            klist = np.zeros_like(Tlist)
            for i in range(len(Tlist)):
                klist[i] = kf.get_rate_coefficient(Tlist[i]) / self.get_equilibrium_constant(Tlist[i])

            kr = KineticsData(Tdata=(Tlist, "K"), kdata=(klist, kunits), Tmin=(np.min(Tlist), "K"),
                              Tmax=(np.max(Tlist), "K"))
            return kr

        elif isinstance(kf, Arrhenius):
            if isinstance(kf, SurfaceArrhenius):
                return self.reverse_surface_arrhenius_rate(kf, kunits, Tmin, Tmax)
            else:
                return self.reverse_arrhenius_rate(kf, kunits, Tmin, Tmax)

        elif network_kinetics and self.network_kinetics is not None:
            kf = self.network_kinetics
            return self.reverse_arrhenius_rate(kf, kunits)

        elif isinstance(kf, Chebyshev):
            Tlist = 1.0 / np.linspace(1.0 / kf.Tmax.value, 1.0 / kf.Tmin.value, 50)
            Plist = np.linspace(kf.Pmin.value, kf.Pmax.value, 20)
            K = np.zeros((len(Tlist), len(Plist)), np.float64)
            for Tindex, T in enumerate(Tlist):
                for Pindex, P in enumerate(Plist):
                    K[Tindex, Pindex] = kf.get_rate_coefficient(T, P) / self.get_equilibrium_constant(T)
            kr = Chebyshev()
            kr.fit_to_data(Tlist, Plist, K, kunits, kf.degreeT, kf.degreeP, kf.Tmin.value, kf.Tmax.value, kf.Pmin.value,
                         kf.Pmax.value)
            return kr

        elif isinstance(kf, PDepArrhenius):
            kr = PDepArrhenius()
            kr.pressures = kf.pressures
            kr.arrhenius = []
            rxn = Reaction(reactants=self.reactants, products=self.products)
            for kinetics in kf.arrhenius:
                rxn.kinetics = kinetics
                kr.arrhenius.append(rxn.generate_reverse_rate_coefficient(kf.Tmin, kf.Tmax))
            return kr

        elif isinstance(kf, MultiArrhenius):
            kr = MultiArrhenius()
            kr.arrhenius = []
            rxn = Reaction(reactants=self.reactants, products=self.products)
            for kinetics in kf.arrhenius:
                rxn.kinetics = kinetics
                kr.arrhenius.append(rxn.generate_reverse_rate_coefficient())
            return kr

        elif isinstance(kf, MultiPDepArrhenius):
            kr = MultiPDepArrhenius()
            kr.arrhenius = []
            rxn = Reaction(reactants=self.reactants, products=self.products)
            for kinetics in kf.arrhenius:
                rxn.kinetics = kinetics
                kr.arrhenius.append(rxn.generate_reverse_rate_coefficient())
            return kr

        elif isinstance(kf, ThirdBody):
            lowPkunits = get_rate_coefficient_units_from_reaction_order(n_gas + 1, n_surf)
            krLow = self.reverse_arrhenius_rate(kf.arrheniusLow, lowPkunits)
            parameters = kf.__reduce__()[1]  # use the pickle helper to get all the other things needed
            kr = ThirdBody(krLow, *parameters[1:])
            return kr

        elif isinstance(kf, Lindemann):
            krHigh = self.reverse_arrhenius_rate(kf.arrheniusHigh, kunits)
            lowPkunits = get_rate_coefficient_units_from_reaction_order(n_gas + 1, n_surf)
            krLow = self.reverse_arrhenius_rate(kf.arrheniusLow, lowPkunits)
            parameters = kf.__reduce__()[1]  # use the pickle helper to get all the other things needed
            kr = Lindemann(krHigh, krLow, *parameters[2:])
            return kr

        elif isinstance(kf, Troe):
            krHigh = self.reverse_arrhenius_rate(kf.arrheniusHigh, kunits)
            lowPkunits = get_rate_coefficient_units_from_reaction_order(n_gas + 1, n_surf)
            krLow = self.reverse_arrhenius_rate(kf.arrheniusLow, lowPkunits)
            parameters = kf.__reduce__()[1]  # use the pickle helper to get all the other things needed
            kr = Troe(krHigh, krLow, *parameters[2:])
            return kr
        else:
            raise ReactionError("Unexpected kinetics type {0}; "
                                "should be one of {1}".format(self.kinetics.__class__, supported_types))

    def calculate_tst_rate_coefficients(self, Tlist):
        return np.array([self.calculate_tst_rate_coefficient(T) for T in Tlist], np.float64)

    def calculate_tst_rate_coefficient(self, T):
        """
        Evaluate the forward rate coefficient for the reaction with
        corresponding transition state `TS` at temperature `T` in K using
        (canonical) transition state theory. The TST equation is

        .. math:: k(T) = \\kappa(T) \\frac{k_\\mathrm{B} T}{h} \\frac{Q^\\ddagger(T)}{Q^\\mathrm{A}(T) Q^\\mathrm{B}(T)} \\exp \\left( -\\frac{E_0}{k_\\mathrm{B} T} \\right)

        where :math:`Q^\\ddagger` is the partition function of the transition state,
        :math:`Q^\\mathrm{A}` and :math:`Q^\\mathrm{B}` are the partition function
        of the reactants, :math:`E_0` is the ground-state energy difference from
        the transition state to the reactants, :math:`T` is the absolute
        temperature, :math:`k_\\mathrm{B}` is the Boltzmann constant, and :math:`h`
        is the Planck constant. :math:`\\kappa(T)` is an optional tunneling
        correction.
        """
        # Determine TST rate constant at each temperature
        Qreac = 1.0
        E0 = 0.0
        for spec in self.reactants:
            logging.debug('    Calculating Partition function for ' + spec.label)
            Qreac *= spec.get_partition_function(T) / (constants.R * T / 101325.)
            E0 -= spec.conformer.E0.value_si
        logging.debug('    Calculating Partition function for ' + self.transition_state.label)
        Qts = self.transition_state.get_partition_function(T) / (constants.R * T / 101325.)
        E0 += self.transition_state.conformer.E0.value_si
        k = (constants.kB * T / constants.h * Qts / Qreac) * math.exp(-E0 / constants.R / T)

        # Apply tunneling correction
        k *= self.transition_state.calculate_tunneling_factor(T)

        return k

    def can_tst(self):
        """
        Return ``True`` if the necessary parameters are available for using
        transition state theory -- or the microcanonical equivalent, RRKM
        theory -- to compute the rate coefficient for this reaction, or
        ``False`` otherwise.
        """
        return len(self.transition_state.conformer.modes) > 0

    def calculate_microcanonical_rate_coefficient(self, e_list, j_list, reac_dens_states, prod_dens_states=None, T=0.0):
        """
        Calculate the microcanonical rate coefficient :math:`k(E)` for the reaction
        `reaction` at the energies `e_list` in J/mol. `reac_dens_states` and
        `prod_dens_states` are the densities of states of the reactant and product
        configurations for this reaction. If the reaction is irreversible, only the
        reactant density of states is required; if the reaction is reversible, then
        both are required. This function will try to use the best method that it
        can based on the input data available:
        
        * If detailed information has been provided for the transition state (i.e.
          the molecular degrees of freedom), then RRKM theory will be used.
        
        * If the above is not possible but high-pressure limit kinetics
          :math:`k_\\infty(T)` have been provided, then the inverse Laplace 
          transform method will be used.
    
        The density of states for the product `prod_dens_states` and the temperature
        of interest `T` in K can also be provided. For isomerization and association
        reactions `prod_dens_states` is required; for dissociation reactions it is
        optional. The temperature is used if provided in the detailed balance
        expression to determine the reverse kinetics, and in certain cases in the
        inverse Laplace transform method.
        """
        return calculate_microcanonical_rate_coefficient(self, e_list, j_list, reac_dens_states, prod_dens_states, T)

    def is_balanced(self):
        """
        Return ``True`` if the reaction has the same number of each atom on
        each side of the reaction equation, or ``False`` if not.
        """
        cython.declare(reactantElements=dict, productElements=dict, molecule=Molecule, atom=Atom, element=Element)

        reactant_elements = {}
        product_elements = {}
        for element in element_list:
            reactant_elements[element] = 0
            product_elements[element] = 0

        for reactant in self.reactants:
            if isinstance(reactant, Species):
                molecule = reactant.molecule[0]
            elif isinstance(reactant, Molecule):
                molecule = reactant
            for atom in molecule.atoms:
                reactant_elements[atom.element] += 1

        for product in self.products:
            if isinstance(product, Species):
                molecule = product.molecule[0]
            elif isinstance(product, Molecule):
                molecule = product
            for atom in molecule.atoms:
                product_elements[atom.element] += 1

        for element in element_list:
            if reactant_elements[element] != product_elements[element]:
                return False

        return True

    def generate_pairs(self):
        """
        Generate the reactant-product pairs to use for this reaction when
        performing flux analysis. The exact procedure for doing so depends on
        the reaction type:
        
        =================== =============== ========================================
        Reaction type       Template        Resulting pairs
        =================== =============== ========================================
        Isomerization       A     -> C      (A,C)
        Dissociation        A     -> C + D  (A,C), (A,D)
        Association         A + B -> C      (A,C), (B,C)
        Bimolecular         A + B -> C + D  (A,C), (B,D) *or* (A,D), (B,C)
        =================== =============== ========================================
        
        There are a number of ways of determining the correct pairing for 
        bimolecular reactions. Here we try a simple similarity analysis by comparing
        the number of heavy atoms. This should work most of the time, but a more
        rigorous algorithm may be needed for some cases.
        """
        self.pairs = []

        if len(self.reactants) == 1 or len(self.products) == 1:
            # Pair each reactant with each product
            for reactant in self.reactants:
                for product in self.products:
                    self.pairs.append((reactant, product))

        else:  # this is the bimolecular case
            reactants = self.reactants[:]
            products = self.products[:]

            def get_sorting_key(spc):
                # List of elements to sort by, order is intentional
                numbers = [6, 8, 7, 14, 16, 15, 17, 53, 9, 35]  # C, O, N, Si, S, P, Cl, I, F, Br
                return tuple(sum([1 for atom in spc.molecule[0].atoms if atom.element.number == n]) for n in numbers)

            # Sort the reactants and products by element counts
            reactants.sort(key=get_sorting_key)
            products.sort(key=get_sorting_key)

            while len(reactants) > 1 and len(products) > 1:
                self.pairs.append((reactants.pop(), products.pop()))

            for reactant in reactants:
                for product in products:
                    self.pairs.append((reactant, product))

    def draw(self, path):
        """
        Generate a pictorial representation of the chemical reaction using the
        :mod:`draw` module. Use `path` to specify the file to save
        the generated image to; the image type is automatically determined by
        extension. Valid extensions are ``.png``, ``.svg``, ``.pdf``, and
        ``.ps``; of these, the first is a raster format and the remainder are
        vector formats.
        """
        from rmgpy.molecule.draw import ReactionDrawer
        img_format = os.path.splitext(path)[1].lower()[1:]
        ReactionDrawer().draw(self, img_format, path)

    def _repr_png_(self):
        """
        Return a png picture of the reaction, useful for ipython-qtconsole.
        """
        from rmgpy.molecule.draw import ReactionDrawer
        temp_file_name = 'temp_reaction.png'
        ReactionDrawer().draw(self, 'png', temp_file_name)
        png = open(temp_file_name, 'rb').read()
        os.unlink(temp_file_name)
        return png

    # Build the transition state geometry
    def generate_3d_ts(self, reactants, products):
        """
        Generate the 3D structure of the transition state. Called from 
        model.generate_kinetics().
        
        self.reactants is a list of reactants
        self.products is a list of products
        """

        """
        Iterate through each reactant, then iterate through its atoms to find the
        atoms involved in the reaction. If a radical is involved, can find the atom
        with radical electrons. If a more reliable method can be found, would greatly
        improve the method.
        
        Repeat for the products
        """
        for i in range(0, len(reactants)):
            mol = reactants[i].molecule[0]
            for j in range(0, mol.rdMol.GetNumAtoms()):
                if mol.rdMol.GetAtomWithIdx(j).GetNumRadicalElectrons():
                    point = mol.rdMol.GetConformer(mol.rdMolConfId).GetAtomPosition(j)
                    neighbor = mol.rdMol.GetAtomWithIdx(j).GetNeighbors()
                    dir_vec = [{} for k in range(len(neighbor))]
                    len_vec = [None] * len(neighbor)
                    for k in range(0, len(neighbor)):
                        new_idx = neighbor[k].GetIdx()
                        new_pt = mol.rdMol.GetConformer(mol.rdMolConfId).GetAtomPosition(new_idx)
                        dir_vec[k] = point.DirectionVector(new_pt)
                        len_vec[k] = point.Distance(new_pt)
                    x_coord = [None] * len(neighbor)
                    y_coord = [None] * len(neighbor)
                    z_coord = [None] * len(neighbor)
                    for k in range(0, len(neighbor)):
                        x_coord[k] = dir_vec[k].x * len_vec[k]
                        y_coord[k] = dir_vec[k].y * len_vec[k]
                        z_coord[k] = dir_vec[k].z * len_vec[k]
            reaction_axis = [sum(x_coord), sum(y_coord), sum(z_coord)]
            reactants[i].reactionAxis = reaction_axis

        for i in range(0, len(products)):
            mol = products[i].molecule[0]
            for j in range(0, mol.rdMol.GetNumAtoms()):
                if mol.rdMol.GetAtomWithIdx(j).GetNumRadicalElectrons():
                    point = mol.rdMol.GetConformer(mol.rdMolConfId).GetAtomPosition(j)
                    neighbor = mol.rdMol.GetAtomWithIdx(j).GetNeighbors()
                    dir_vec = [{} for k in range(len(neighbor))]
                    len_vec = [None] * len(neighbor)
                    for k in range(0, len(neighbor)):
                        new_idx = neighbor[k].GetIdx()
                        new_pt = mol.rdMol.GetConformer(mol.rdMolConfId).GetAtomPosition(new_idx)
                        dir_vec[k] = point.DirectionVector(new_pt)
                        len_vec[k] = point.Distance(new_pt)
                    x_coord = [None] * len(neighbor)
                    y_coord = [None] * len(neighbor)
                    z_coord = [None] * len(neighbor)
                    for k in range(0, len(neighbor)):
                        x_coord[k] = dir_vec[k].x * len_vec[k]
                        y_coord[k] = dir_vec[k].y * len_vec[k]
                        z_coord[k] = dir_vec[k].z * len_vec[k]
            reaction_axis = [sum(x_coord), sum(y_coord), sum(z_coord)]
            products[i].reactionAxis = reaction_axis

    def copy(self):
        """
        Create a deep copy of the current reaction.
        """

        cython.declare(other=Reaction)

        other = Reaction.__new__(Reaction)
        other.index = self.index
        other.label = self.label
        other.reactants = []
        for reactant in self.reactants:
            other.reactants.append(reactant.copy(deep=True))
        other.products = []
        for product in self.products:
            other.products.append(product.copy(deep=True))
        other.degeneracy = self.degeneracy
        other.specific_collider = self.specific_collider
        other.kinetics = deepcopy(self.kinetics)
        other.network_kinetics = deepcopy(self.network_kinetics)
        other.reversible = self.reversible
        other.transition_state = deepcopy(self.transition_state)
        other.duplicate = self.duplicate
        other.pairs = deepcopy(self.pairs)
        other.allow_pdep_route = self.allow_pdep_route
        other.elementary_high_p = self.elementary_high_p
        other.comment = deepcopy(self.comment)

        return other

    def ensure_species(self, reactant_resonance=False, product_resonance=False):
        """
        Ensure the reaction contains species objects in its reactant and product
        attributes. If the reaction is found to hold molecule objects, it
        modifies the reactant, product and pairs to hold
        Species objects.

        Generates resonance structures for Molecules if the corresponding options,
        reactant_resonance and/or product_resonance, are True. Does not generate
        resonance for reactants or products that start as Species objects.
        """
        from rmgpy.data.kinetics.common import ensure_species
        # if already species' objects, return none
        if isinstance(self.reactants[0], Species):
            return None
        # obtain species with all resonance isomers
        if self.is_forward:
            ensure_species(self.reactants, resonance=reactant_resonance, keep_isomorphic=True)
            ensure_species(self.products, resonance=product_resonance, keep_isomorphic=True)
        else:
            ensure_species(self.reactants, resonance=product_resonance, keep_isomorphic=True)
            ensure_species(self.products, resonance=reactant_resonance, keep_isomorphic=True)

        # convert reaction.pairs object to species
        if self.pairs:
            new_pairs = []
            for reactant, product in self.pairs:
                new_pair = []
                for reactant0 in self.reactants:
                    if reactant0.is_isomorphic(reactant):
                        new_pair.append(reactant0)
                        break
                for product0 in self.products:
                    if product0.is_isomorphic(product):
                        new_pair.append(product0)
                        break
                new_pairs.append(new_pair)
            self.pairs = new_pairs

        try:
            self.reverse.ensure_species()
        except AttributeError:
            pass

    def check_collision_limit_violation(self, t_min, t_max, p_min, p_max):
        """
        Warn if a core reaction violates the collision limit rate in either the forward or reverse direction
        at the relevant extreme T/P conditions. Assuming a monotonic behaviour of the kinetics.
        Returns a list with the reaction object and the direction in which the violation was detected.
        """
        conditions = [[t_min, p_min]]
        if t_min != t_max:
            conditions.append([t_max, p_min])
        if self.kinetics.is_pressure_dependent() and p_max != p_min:
            conditions.append([t_min, p_max])
            if t_min != t_max:
                conditions.append([t_max, p_max])
        logging.debug("Checking whether reaction {0} violates the collision rate limit...".format(self))
        violator_list = []
        kf_list = []
        kr_list = []
        collision_limit_f = []
        collision_limit_r = []
        for condition in conditions:
            if len(self.reactants) >= 2:
                try:
                    collision_limit_f.append(self.calculate_coll_limit(temp=condition[0], reverse=False))
                except ValueError:
                    continue
                else:
                    kf_list.append(self.get_rate_coefficient(condition[0], condition[1]))
            if len(self.products) >= 2:
                try:
                    collision_limit_r.append(self.calculate_coll_limit(temp=condition[0], reverse=True))
                except ValueError:
                    continue
                else:
                    kr_list.append(self.generate_reverse_rate_coefficient().get_rate_coefficient(condition[0], condition[1]))
        if len(self.reactants) >= 2:
            for i, k in enumerate(kf_list):
                if k > collision_limit_f[i]:
                    ratio = k / collision_limit_f[i]
                    condition = '{0} K, {1:.1f} bar'.format(conditions[i][0], conditions[i][1] / 1e5)
                    violator_list.append([self, 'forward', ratio, condition])
        if len(self.products) >= 2:
            for i, k in enumerate(kr_list):
                if k > collision_limit_r[i]:
                    ratio = k / collision_limit_r[i]
                    condition = '{0} K, {1:.1f} bar'.format(conditions[i][0], conditions[i][1] / 1e5)
                    violator_list.append([self, 'reverse', ratio, condition])
        return violator_list

    def calculate_coll_limit(self, temp, reverse=False):
        """
        Calculate the collision limit rate in m3/mol-s for the given temperature
        implemented as recommended in Wang et al. doi 10.1016/j.combustflame.2017.08.005 (Eq. 1)
        """
        reduced_mass = self.get_reduced_mass(reverse)
        sigma, epsilon = self.get_mean_sigma_and_epsilon(reverse)
        Tr = temp * constants.kB * constants.Na / epsilon
        reduced_coll_integral = 1.16145 * Tr ** (-0.14874) + 0.52487 * math.exp(-0.7732 * Tr) + 2.16178 * math.exp(
            -2.437887 * Tr)
        k_coll = (math.sqrt(8 * math.pi * constants.kB * temp * constants.Na / reduced_mass) * sigma ** 2
                  * reduced_coll_integral * constants.Na)
        return k_coll

    def get_reduced_mass(self, reverse=False):
        """
        Returns the reduced mass of the reactants if reverse is ``False``
        Returns the reduced mass of the products if reverse is ``True``
        """
        if reverse:
            mass_list = [spc.molecule[0].get_molecular_weight() for spc in self.products]
        else:
            mass_list = [spc.molecule[0].get_molecular_weight() for spc in self.reactants]
        reduced_mass = reduce((lambda x, y: x * y), mass_list) / sum(mass_list)
        return reduced_mass

    def get_mean_sigma_and_epsilon(self, reverse=False):
        """
        Calculates the collision diameter (sigma) using an arithmetic mean
        Calculates the well depth (epsilon) using a geometric mean
        If reverse is ``False`` the above is calculated for the reactants, otherwise for the products
        """
        sigmas = []
        epsilons = []
        if reverse:
            for spc in self.products:
                trans = spc.get_transport_data()
                sigmas.append(trans.sigma.value_si)
                epsilons.append(trans.epsilon.value_si)
            num_of_spcs = len(self.products)
        else:
            for spc in self.reactants:
                trans = spc.get_transport_data()
                sigmas.append(trans.sigma.value_si)
                epsilons.append(trans.epsilon.value_si)
            num_of_spcs = len(self.reactants)
        if any([x == 0 for x in sigmas + epsilons]):
            raise ValueError
        mean_sigmas = sum(sigmas) / num_of_spcs
        mean_epsilons = reduce((lambda x, y: x * y), epsilons) ** (1 / len(epsilons))
        return mean_sigmas, mean_epsilons

    def generate_high_p_limit_kinetics(self):
        """
        Used for incorporating library reactions with pressure-dependent kinetics in PDep networks.
        Only implemented for LibraryReaction
        """
        raise NotImplementedError("generate_high_p_limit_kinetics is not implemented for all Reaction subclasses.")


def same_species_lists(list1, list2, check_identical=False, only_check_label=False, generate_initial_map=False,
                       strict=True, save_order=False):
    """
    This method compares whether two lists of species or molecules are the same,
    given the comparison options provided. It is used for the `is_same` method
    of :class:`Reaction`, but may also be useful in other situations.

    Args:
        list1 (list):                          list of :class:`Species` or :class:`Molecule` objects
        list2 (list):                          list of :class:`Species` or :class:`Molecule` objects
        check_identical (bool, optional):      if ``True``, use is_identical comparison and compare atom IDs
        only_check_label (bool, optional):     if ``True``, only compare the label attribute of each species
        generate_initial_map (bool, optional): if ``True``, initialize map by pairing atoms with same labels
        strict (bool, optional):               if ``False``, perform isomorphism ignoring electrons
        save_order (bool, optional):           if ``True``, perform isomorphism saving atom order

    Returns:
        ``True`` if the lists are the same and ``False`` otherwise
    """

    def same(object1, object2, _check_identical=check_identical, _only_check_label=only_check_label,
             _generate_initial_map=generate_initial_map, _strict=strict, save_order=save_order):
        if _only_check_label:
            return str(object1) == str(object2)
        elif _check_identical:
            return object1.is_identical(object2, strict=_strict)
        else:
            return object1.is_isomorphic(object2, generate_initial_map=_generate_initial_map,
                                         strict=_strict, save_order=save_order)

    if len(list1) == len(list2) == 1:
        if same(list1[0], list2[0]):
            return True
    elif len(list1) == len(list2) == 2:
        if same(list1[0], list2[0]) and same(list1[1], list2[1]):
            return True
        elif same(list1[0], list2[1]) and same(list1[1], list2[0]):
            return True
    elif len(list1) == len(list2) == 3:
        if same(list1[0], list2[0]):
            if same(list1[1], list2[1]):
                if same(list1[2], list2[2]):
                    return True
            elif same(list1[1], list2[2]):
                if same(list1[2], list2[1]):
                    return True
        elif same(list1[0], list2[1]):
            if same(list1[1], list2[0]):
                if same(list1[2], list2[2]):
                    return True
            elif same(list1[1], list2[2]):
                if same(list1[2], list2[0]):
                    return True
        elif same(list1[0], list2[2]):
            if same(list1[1], list2[0]):
                if same(list1[2], list2[1]):
                    return True
            elif same(list1[1], list2[1]):
                if same(list1[2], list2[0]):
                    return True
    elif len(list1) == len(list2):
        raise NotImplementedError("Can't check isomorphism of lists with {0} species/molecules".format(len(list1)))
    # nothing found
    return False
