#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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

import cython
import math
import numpy
import logging
import re
import os.path
from copy import copy, deepcopy
import urllib

import rmgpy.constants as constants
from rmgpy.molecule.molecule import Molecule, Atom
from rmgpy.molecule.element import Element
from rmgpy.species import Species
from rmgpy.kinetics.arrhenius import Arrhenius #PyDev: @UnresolvedImport
from rmgpy.kinetics import KineticsData, ArrheniusEP, ThirdBody, Lindemann, Troe, Chebyshev, PDepArrhenius, MultiArrhenius, MultiPDepArrhenius, getRateCoefficientUnitsFromReactionOrder  #PyDev: @UnresolvedImport
from rmgpy.pdep.reaction import calculateMicrocanonicalRateCoefficient
from rmgpy.exceptions import ReactionError
from rmgpy.kinetics.diffusionLimited import diffusionLimiter

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
    'specificCollider'  :class:`Species`            The collider species (as a :class:`Species` object)
    `kinetics`          :class:`KineticsModel`      The kinetics model to use for the reaction
    `network_kinetics`  :class:`Arrhenius`          The kinetics model to use for PDep network exploration if the `kinetics` attribute is :class:PDepKineticsModel:
    `reversible`        ``bool``                    ``True`` if the reaction is reversible, ``False`` if not
    `transitionState`   :class:`TransitionState`    The transition state
    `duplicate`         ``bool``                    ``True`` if the reaction is known to be a duplicate, ``False`` if not
    `degeneracy`        :class:`double`             The reaction path degeneracy for the reaction
    `pairs`             ``list``                    Reactant-product pairings to use in converting reaction flux to species flux
    `allow_pdep_route`  ``bool``                    ``True`` if the reaction has an additional PDep pathway, ``False`` if not (by default), used for LibraryReactions
    `elementary_high_p` ``bool``                    If ``True``, pressure dependent kinetics will be generated (relevant only for unimolecular library reactions)
                                                    If ``False`` (by default), this library reaction will not be explored.
                                                    Only unimolecular library reactions with high pressure limit kinetics should be flagged (not if the kinetics were measured at some relatively low pressure)
    `comment`           ``str``                     A description of the reaction source (optional)
    =================== =========================== ============================
    
    """
    
    def __init__(self,
                 index=-1,
                 label='',
                 reactants=None,
                 products=None,
                 specificCollider=None,
                 kinetics=None,
                 network_kinetics=None,
                 reversible=True,
                 transitionState=None,
                 duplicate=False,
                 degeneracy=1,
                 pairs=None,
                 allow_pdep_route=False,
                 elementary_high_p=False,
                 comment='',
                 ):
        self.index = index
        self.label = label
        self.reactants = reactants
        self.products = products
        self.specificCollider = specificCollider
        self._degeneracy = degeneracy
        self.kinetics = kinetics
        self.network_kinetics = network_kinetics
        self.reversible = reversible
        self.transitionState = transitionState
        self.duplicate = duplicate
        self.pairs = pairs
        self.allow_pdep_route = allow_pdep_route
        self.elementary_high_p = elementary_high_p
        self.comment = comment
        self.k_effective_cache = {}

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
        if self.specificCollider is not None: string += 'specificCollider={0!r}, '.format(self.specificCollider)
        if self.kinetics is not None: string += 'kinetics={0!r}, '.format(self.kinetics)
        if self.network_kinetics is not None: string += 'network_kinetics={0!r}, '.format(self.network_kinetics)
        if not self.reversible: string += 'reversible={0}, '.format(self.reversible)
        if self.transitionState is not None: string += 'transitionState={0!r}, '.format(self.transitionState)
        if self.duplicate: string += 'duplicate={0}, '.format(self.duplicate)
        if self.degeneracy != 1: string += 'degeneracy={0:.1f}, '.format(self.degeneracy)
        if self.pairs is not None: string += 'pairs={0}, '.format(self.pairs)
        if self.allow_pdep_route: string += 'allow_pdep_route={0}, '.format(self.allow_pdep_route)
        if self.elementary_high_p: string += 'elementary_high_p={0}, '.format(self.elementary_high_p)
        if self.comment != '': string += 'comment={0!r}, '.format(self.comment)
        string = string[:-2] + ')'
        return string

    def __str__(self):
        """
        Return a string representation of the reaction, in the form 'A + B <=> C + D'.
        If a specificCollider exists, the srting representation is 'A + B (+S) <=> C + D (+S)'.
        """
        return self.toLabeledStr(use_index=True)
    
    def toLabeledStr(self, use_index=False):
        """
        the same as __str__ except that the labels are assumed to exist and used for reactant and products rather than 
        the labels plus the index in parentheses
        """
        arrow = ' <=> '
        if not self.reversible: arrow = ' => '
        
        if self.specificCollider:
            return arrow.join([' + '.join([str(s) if use_index else s.label for s in self.reactants])+' (+'+str(self.specificCollider)+')', ' + '.join([str(s) if use_index else s.label for s in self.products])+' (+'+str(self.specificCollider)+')'])
        else:
            return arrow.join([' + '.join([str(s) if use_index else s.label for s in self.reactants]), ' + '.join([str(s) if use_index else s.label for s in self.products])])

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Reaction, (self.index,
                           self.label,
                           self.reactants,
                           self.products,
                           self.specificCollider,
                           self.kinetics,
                           self.network_kinetics,
                           self.reversible,
                           self.transitionState,
                           self.duplicate,
                           self.degeneracy,
                           self.pairs,
                           self.allow_pdep_route,
                           self.elementary_high_p,
                           self.comment
                           ))

    def __getDegneneracy(self):
        return self._degeneracy

    def __setDegeneracy(self, new):
        # modify rate if kinetics exists
        if self.kinetics is not None:
            if self._degeneracy < 2:
                degeneracyRatio = new
            else:
                degeneracyRatio = (new*1.0) / self._degeneracy
            # fix kinetics comment with new degeneracy
            if 'Multiplied by reaction path degeneracy {}'.format(self._degeneracy) in self.kinetics.comment:
                self.kinetics.comment = self.kinetics.comment.replace(
                                                  'Multiplied by reaction path degeneracy {}'.format(self._degeneracy),
                                                  'Multiplied by reaction path degeneracy {}'.format(float(new)))
            elif self.kinetics.comment:
                self.kinetics.comment += 'Multiplied by reaction path degeneracy {}'.format(float(new))
            self.kinetics.changeRate(degeneracyRatio)
        # set new degeneracy
        self._degeneracy = new
    degeneracy = property(__getDegneneracy, __setDegeneracy)

    def toChemkin(self, speciesList=None, kinetics=True):
        """
        Return the chemkin-formatted string for this reaction.
        
        If `kinetics` is set to True, the chemkin format kinetics will also
        be returned (requires the `speciesList` to figure out third body colliders.)
        Otherwise, only the reaction string will be returned.
        """
        import rmgpy.chemkin
        if kinetics:
            return rmgpy.chemkin.writeKineticsEntry(self, speciesList)
        else:
            return rmgpy.chemkin.writeReactionString(self)
    
    def toCantera(self, speciesList=None, useChemkinIdentifier = False):
        """
        Converts the RMG Reaction object to a Cantera Reaction object
        with the appropriate reaction class.

        If useChemkinIdentifier is set to False, the species label is used
        instead. Be sure that species' labels are unique when setting it False.
        """
        from rmgpy.kinetics import Arrhenius, ArrheniusEP, MultiArrhenius, PDepArrhenius, MultiPDepArrhenius, Chebyshev, ThirdBody, Lindemann, Troe
                    
        import cantera as ct
        
        if speciesList is None:
            speciesList = []
        
        # Create the dictionaries containing species strings and their stoichiometries
        # for initializing the cantera reaction object
        ctReactants = {}
        ctCollider = {}
        for reactant in self.reactants:
            if useChemkinIdentifier:
                reactantName = reactant.toChemkin()
            else:
                reactantName = reactant.label
            if reactantName in ctReactants:
                ctReactants[reactantName] += 1
            else:
                ctReactants[reactantName] = 1
        ctProducts = {}
        for product in self.products:
            if useChemkinIdentifier:
                productName = product.toChemkin()
            else:
                productName = product.label
            if productName in ctProducts:
                ctProducts[productName] += 1
            else:
                ctProducts[productName] = 1
        if self.specificCollider:              # add a specific collider if exists
            ctCollider[self.specificCollider.toChemkin() if useChemkinIdentifier else self.specificCollider.label] = 1
                
        if self.kinetics:
            if isinstance(self.kinetics, Arrhenius):
                # Create an Elementary Reaction
                ctReaction = ct.ElementaryReaction(reactants=ctReactants, products=ctProducts)
            elif isinstance(self.kinetics, MultiArrhenius):
                # Return a list of elementary reactions which are duplicates
                ctReaction = [ct.ElementaryReaction(reactants=ctReactants, products=ctProducts) for arr in self.kinetics.arrhenius]
                
            elif isinstance(self.kinetics, PDepArrhenius):
                ctReaction = ct.PlogReaction(reactants=ctReactants, products=ctProducts)
                
            elif isinstance(self.kinetics, MultiPDepArrhenius):
                ctReaction = [ct.PlogReaction(reactants=ctReactants, products=ctProducts) for arr in self.kinetics.arrhenius]
                
            
            elif isinstance(self.kinetics, Chebyshev):
                ctReaction = ct.ChebyshevReaction(reactants=ctReactants, products=ctProducts)
            
            elif isinstance(self.kinetics, ThirdBody):
                if ctCollider is not None:
                    ctReaction = ct.ThreeBodyReaction(reactants=ctReactants, products=ctProducts, tbody=ctCollider)
                else:
                    ctReaction = ct.ThreeBodyReaction(reactants=ctReactants, products=ctProducts)
                
            elif isinstance(self.kinetics, Lindemann) or isinstance(self.kinetics, Troe):
                if ctCollider is not None:
                    ctReaction = ct.FalloffReaction(reactants=ctReactants, products=ctProducts, tbody=ctCollider)
                else:
                    ctReaction = ct.FalloffReaction(reactants=ctReactants, products=ctProducts)
            else:
                raise NotImplementedError('Unable to set cantera kinetics for {0}'.format(self.kinetics))
            
            
            # Set reversibility, duplicate, and ID attributes
            if isinstance(ctReaction,list):
                for rxn in ctReaction:
                    rxn.reversible = self.reversible
                    # Set the duplicate flag to true since this reaction comes from multiarrhenius or multipdeparrhenius 
                    rxn.duplicate = True
                    # Set the ID flag to the original rmg index 
                    rxn.ID = str(self.index) 
            else:
                ctReaction.reversible = self.reversible
                ctReaction.duplicate = self.duplicate
                ctReaction.ID = str(self.index)
                
            
            self.kinetics.setCanteraKinetics(ctReaction, speciesList)
            
            return ctReaction
                
        else:
            raise Exception('Cantera reaction cannot be created because there was no kinetics.')
    
    def getURL(self):
        """
        Get a URL to search for this reaction in the rmg website.
        """
        # eg. http://dev.rmg.mit.edu/database/kinetics/reaction/reactant1=1%20C%200%20%7B2,S%7D;2%20O%200%20%7B1,S%7D;__reactant2=1%20C%202T;__product1=1%20C%201;__product2=1%20C%200%20%7B2,S%7D;2%20O%201%20%7B1,S%7D;

        base_url = "http://rmg.mit.edu/database/kinetics/reaction/"

        rxn_string = ''
        for i,species in enumerate(self.reactants):
            adjlist = species.molecule[0].toAdjacencyList(removeH=False)
            rxn_string += "reactant{0}={1}__".format(i+1, adjlist)
        for i,species in enumerate(self.products):
            adjlist = species.molecule[0].toAdjacencyList(removeH=False)
            rxn_string += "product{0}={1}__".format(i+1, adjlist)

        url = base_url + urllib.quote(rxn_string)
        return url.strip('_')
        
    def isIsomerization(self):
        """
        Return ``True`` if the reaction represents an isomerization reaction
        :math:`\\ce{A <=> B}` or ``False`` if not.
        """
        return len(self.reactants) == 1 and len(self.products) == 1

    def isAssociation(self):
        """
        Return ``True`` if the reaction represents an association reaction
        :math:`\\ce{A + B <=> C}` or ``False`` if not.
        """
        return len(self.reactants) > 1 and len(self.products) == 1

    def isDissociation(self):
        """
        Return ``True`` if the reaction represents a dissociation reaction
        :math:`\\ce{A <=> B + C}` or ``False`` if not.
        """
        return len(self.reactants) == 1 and len(self.products) > 1
    
    def isUnimolecular(self):
        """
        Return ``True`` if the reaction has a single molecule as either reactant or product (or both)
        :math:`\\ce{A <=> B + C}` or :math:`\\ce{A + B <=> C}` or :math:`\\ce{A <=> B}`,
        or ``False`` if not.
        """
        return len(self.reactants) == 1 or len(self.products) == 1

    def hasTemplate(self, reactants, products):
        """
        Return ``True`` if the reaction matches the template of `reactants`
        and `products`, which are both lists of :class:`Species` objects, or
        ``False`` if not.
        """
        return ((all([spec in self.reactants for spec in reactants]) and
            all([spec in self.products for spec in products])) or
            (all([spec in self.products for spec in reactants]) and
            all([spec in self.reactants for spec in products])))

    def matchesSpecies(self, reactants, products=None):
        """
        Compares the provided reactants and products against the reactants
        and products of this reaction. Both directions are checked.

        Args:
            reactants (list): Species required on one side of the reaction
            products (list, optional): Species required on the other side
        """
        # Check forward direction
        if _isomorphicSpeciesList(self.reactants, reactants):
            if products is None or _isomorphicSpeciesList(self.products, products):
                return True
            else:
                return False
        elif _isomorphicSpeciesList(self.products, reactants):
            if products is None or _isomorphicSpeciesList(self.reactants, products):
                return True
            else:
                return False
        else:
            return False

    def isIsomorphic(self, other, eitherDirection=True, checkIdentical = False,
                     checkOnlyLabel = False, checkTemplateRxnProducts=False):
        """
        Return ``True`` if this reaction is the same as the `other` reaction,
        or ``False`` if they are different. The comparison involves comparing
        isomorphism of reactants and products, and doesn't use any kinetic
        information.

        If `eitherDirection=False` then the directions must match.

        `checkIdentical` indicates that atom ID's must match and is used in
                        checking degeneracy
        `checkOnlyLabel` indicates that the string representation will be 
                        checked, ignoring the molecular structure comparisons
        `checkTemplateRxnProducts` indicates that only the products of the
                        reaction are checked for isomorphism. This is used when
                        we know the reactants are identical, i.e. in generating
                        reactions.
        """
        if checkTemplateRxnProducts:
            try:
                species1 = self.products if self.isForward else self.reactants
                species2 = other.products if other.isForward else other.reactants
            except AttributeError:
                raise TypeError('Only use checkTemplateRxnProducts flag for TemplateReactions.')

            return _isomorphicSpeciesList(species1, species2,
                                          checkIdentical=checkIdentical,
                                          checkOnlyLabel=checkOnlyLabel)

        # Compare reactants to reactants
        forwardReactantsMatch = _isomorphicSpeciesList(self.reactants, 
                                    other.reactants,checkIdentical = checkIdentical,
                                    checkOnlyLabel = checkOnlyLabel)
        
        # Compare products to products
        forwardProductsMatch = _isomorphicSpeciesList(self.products, 
                                    other.products,checkIdentical = checkIdentical,
                                    checkOnlyLabel = checkOnlyLabel)

        # Compare specificCollider to specificCollider
        ColliderMatch = (self.specificCollider == other.specificCollider)

        # Return now, if we can
        if (forwardReactantsMatch and forwardProductsMatch and ColliderMatch):
            return True
        if not eitherDirection:
            return False
        
        # Compare reactants to products
        reverseReactantsMatch = _isomorphicSpeciesList(self.reactants, 
                                    other.products,checkIdentical = checkIdentical,
                                    checkOnlyLabel = checkOnlyLabel)

        # Compare products to reactants
        reverseProductsMatch = _isomorphicSpeciesList(self.products, 
                                    other.reactants,checkIdentical = checkIdentical,
                                    checkOnlyLabel = checkOnlyLabel)

        # should have already returned if it matches forwards, or we're not allowed to match backwards
        return  (reverseReactantsMatch and reverseProductsMatch and ColliderMatch)

    def getEnthalpyOfReaction(self, T):
        """
        Return the enthalpy of reaction in J/mol evaluated at temperature
        `T` in K.
        """
        cython.declare(dHrxn=cython.double, reactant=Species, product=Species)
        dHrxn = 0.0
        for reactant in self.reactants:
            dHrxn -= reactant.getEnthalpy(T)
        for product in self.products:
            dHrxn += product.getEnthalpy(T)
        return dHrxn

    def getEntropyOfReaction(self, T):
        """
        Return the entropy of reaction in J/mol*K evaluated at temperature `T`
        in K.
        """
        cython.declare(dSrxn=cython.double, reactant=Species, product=Species)
        dSrxn = 0.0
        for reactant in self.reactants:
            dSrxn -= reactant.getEntropy(T)
        for product in self.products:
            dSrxn += product.getEntropy(T)
        return dSrxn

    def getFreeEnergyOfReaction(self, T):
        """
        Return the Gibbs free energy of reaction in J/mol evaluated at
        temperature `T` in K.
        """
        cython.declare(dGrxn=cython.double, reactant=Species, product=Species)
        dGrxn = 0.0
        for reactant in self.reactants:
            try:
                dGrxn -= reactant.getFreeEnergy(T)
            except Exception:
                logging.error("Problem with reactant {!r} in reaction {!s}".format(reactant, self))
                raise
        for product in self.products:
            try:
                dGrxn += product.getFreeEnergy(T)
            except Exception as e:
                logging.error("Problem with product {!r} in reaction {!s}".format(reactant, self))
                raise
        return dGrxn

    def getEquilibriumConstant(self, T, type='Kc'):
        """
        Return the equilibrium constant for the reaction at the specified
        temperature `T` in K. The `type` parameter lets	you specify the
        quantities used in the equilibrium constant: ``Ka`` for	activities,
        ``Kc`` for concentrations (default), or ``Kp`` for pressures. Note that
        this function currently assumes an ideal gas mixture.
        """
        cython.declare(dGrxn=cython.double, K=cython.double, C0=cython.double, P0=cython.double)
        # Use free energy of reaction to calculate Ka
        dGrxn = self.getFreeEnergyOfReaction(T)
        K = numpy.exp(-dGrxn / constants.R / T)
        # Convert Ka to Kc or Kp if specified
        P0 = 1e5
        if type == 'Kc':
            # Convert from Ka to Kc; C0 is the reference concentration
            C0 = P0 / constants.R / T
            K *= C0 ** (len(self.products) - len(self.reactants))
        elif type == 'Kp':
            # Convert from Ka to Kp; P0 is the reference pressure
            K *= P0 ** (len(self.products) - len(self.reactants))
        elif type != 'Ka' and type != '':
            raise ReactionError('Invalid type "%s" passed to Reaction.getEquilibriumConstant(); should be "Ka", "Kc", or "Kp".')
        if K == 0:
            raise ReactionError('Got equilibrium constant of 0')
        return K

    def getEnthalpiesOfReaction(self, Tlist):
        """
        Return the enthalpies of reaction in J/mol evaluated at temperatures
        `Tlist` in K.
        """
        return numpy.array([self.getEnthalpyOfReaction(T) for T in Tlist], numpy.float64)

    def getEntropiesOfReaction(self, Tlist):
        """
        Return the entropies of reaction in J/mol*K evaluated at temperatures
        `Tlist` in K.
        """
        return numpy.array([self.getEntropyOfReaction(T) for T in Tlist], numpy.float64)

    def getFreeEnergiesOfReaction(self, Tlist):
        """
        Return the Gibbs free energies of reaction in J/mol evaluated at
        temperatures `Tlist` in K.
        """
        return numpy.array([self.getFreeEnergyOfReaction(T) for T in Tlist], numpy.float64)

    def getEquilibriumConstants(self, Tlist, type='Kc'):
        """
        Return the equilibrium constants for the reaction at the specified
        temperatures `Tlist` in K. The `type` parameter lets you specify the
        quantities used in the equilibrium constant: ``Ka`` for	activities,
        ``Kc`` for concentrations (default), or ``Kp`` for pressures. Note that
        this function currently assumes an ideal gas mixture.
        """
        return numpy.array([self.getEquilibriumConstant(T, type) for T in Tlist], numpy.float64)

    def getStoichiometricCoefficient(self, spec):
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

    def getRateCoefficient(self, T, P=0):
        """
        Return the overall rate coefficient for the forward reaction at
        temperature `T` in K and pressure `P` in Pa, including any reaction
        path degeneracies.
        
        If diffusionLimiter is enabled, the reaction is in the liquid phase and we use
        a diffusion limitation to correct the rate. If not, then use the intrinsic rate
        coefficient.
        """
        if diffusionLimiter.enabled:
            try:
                k = self.k_effective_cache[T]
            except KeyError:
                k = diffusionLimiter.getEffectiveRate(self, T)
                self.k_effective_cache[T] = k
            return k
        else:
            return self.kinetics.getRateCoefficient(T, P)

    def fixDiffusionLimitedA(self, T):
        """
        Decrease the pre-exponential factor (A) by the diffusion factor
        to account for the diffusion limit at the specified temperature.
        """
        if not diffusionLimiter.enabled:
            return
        # Obtain effective rate
        try:
            k = self.k_effective_cache[T]
        except KeyError:
            k = diffusionLimiter.getEffectiveRate(self, T)
            self.k_effective_cache[T] = k

        # calculate diffusion factor
        diffusionFactor = k / self.kinetics.getRateCoefficient(T, P=0)
        # update preexponential factor
        self.kinetics.A = self.kinetics.A * diffusionFactor
        # Add a comment to self.kinetics.comment
        self.kinetics.comment.append(
            ("Pre-exponential factor A has been decreased by the "
             "diffusion factor {0.2g} evaluated at {1} K.").format(
                diffusionFactor, T))

    def fixBarrierHeight(self, forcePositive=False):
        """
        Turns the kinetics into Arrhenius (if they were ArrheniusEP)
        and ensures the activation energy is at least the endothermicity
        for endothermic reactions, and is not negative only as a result 
        of using Evans Polanyi with an exothermic reaction.
        If `forcePositive` is True, then all reactions
        are forced to have a non-negative barrier.
        """
        cython.declare(H0=cython.double, H298=cython.double, Ea=cython.double)

        H298 = self.getEnthalpyOfReaction(298)
        H0 = sum([spec.getThermoData().E0.value_si for spec in self.products]) \
            - sum([spec.getThermoData().E0.value_si for spec in self.reactants])
        if isinstance(self.kinetics, ArrheniusEP):
            Ea = self.kinetics.E0.value_si # temporarily using Ea to store the intrinsic barrier height E0
            self.kinetics = self.kinetics.toArrhenius(H298)
            if self.kinetics.Ea.value_si < 0.0 and self.kinetics.Ea.value_si < Ea:
                # Calculated Ea (from Evans-Polanyi) is negative AND below than the intrinsic E0
                Ea = min(0.0,Ea) # (the lowest we want it to be)
                self.kinetics.comment += "\nEa raised from {0:.1f} to {1:.1f} kJ/mol.".format(self.kinetics.Ea.value_si/1000., Ea/1000.)
                logging.info("For reaction {0!s} Ea raised from {1:.1f} to {2:.1f} kJ/mol.".format(self, self.kinetics.Ea.value_si/1000., Ea/1000.))
                self.kinetics.Ea.value_si = Ea
        if isinstance(self.kinetics, Arrhenius):
            Ea = self.kinetics.Ea.value_si
            if H0 >= 0 and Ea < H0:
                self.kinetics.Ea.value_si = H0
                self.kinetics.comment += "\nEa raised from {0:.1f} to {1:.1f} kJ/mol to match endothermicity of reaction.".format(Ea/1000.,H0/1000.)
                logging.info("For reaction {2!s}, Ea raised from {0:.1f} to {1:.1f} kJ/mol to match endothermicity of reaction.".format(Ea/1000., H0/1000., self))
        if forcePositive and isinstance(self.kinetics, Arrhenius) and self.kinetics.Ea.value_si < 0:
            self.kinetics.comment += "\nEa raised from {0:.1f} to 0 kJ/mol.".format(self.kinetics.Ea.value_si/1000.)
            logging.info("For reaction {1!s} Ea raised from {0:.1f} to 0 kJ/mol.".format(self.kinetics.Ea.value_si/1000., self))
            self.kinetics.Ea.value_si = 0
        if self.kinetics.isPressureDependent() and self.network_kinetics is not None:
            Ea = self.network_kinetics.Ea.value_si
            if H0 >= 0 and Ea < H0:
                self.network_kinetics.Ea.value_si = H0
                self.network_kinetics.comment += "\nEa raised from {0:.1f} to {1:.1f} kJ/mol to match endothermicity of" \
                                                 " reaction.".format(Ea / 1000., H0 / 1000.)
                logging.info("For reaction {2!s}, Ea of the high pressure limit kinetics raised from {0:.1f} to {1:.1f}"
                             " kJ/mol to match endothermicity of reaction.".format(Ea / 1000., H0 / 1000., self))
            if forcePositive and isinstance(self.kinetics, Arrhenius) and self.kinetics.Ea.value_si < 0:
                self.network_kinetics.comment += "\nEa raised from {0:.1f} to 0 kJ/mol.".format(self.kinetics.Ea.value_si / 1000.)
                logging.info("For reaction {1!s} Ea of the high pressure limit kinetics raised from {0:.1f} to 0"
                             " kJ/mol.".format(self.kinetics.Ea.value_si / 1000.,self))
                self.kinetics.Ea.value_si = 0

    def reverseThisArrheniusRate(self, kForward, reverseUnits):
        """
        Reverses the given kForward, which must be an Arrhenius type.
        You must supply the correct units for the reverse rate.
        The equilibrium constant is evaluated from the current reaction instance (self).
        """
        cython.declare(kf=Arrhenius, kr=Arrhenius)
        cython.declare(Tlist=numpy.ndarray, klist=numpy.ndarray, i=cython.int)
        kf = kForward
        assert isinstance(kf, Arrhenius), "Only reverses Arrhenius rates"
        Tlist = 1.0 / numpy.arange(0.0005, 0.0034, 0.0001)  # 294 K to 2000 K
        # Determine the values of the reverse rate coefficient k_r(T) at each temperature
        klist = numpy.zeros_like(Tlist)
        for i in range(len(Tlist)):
            klist[i] = kf.getRateCoefficient(Tlist[i]) / self.getEquilibriumConstant(Tlist[i])
        kr = Arrhenius()
        kr.fitToData(Tlist, klist, reverseUnits, kf.T0.value_si)
        return kr
        
    def generateReverseRateCoefficient(self, network_kinetics=False):
        """
        Generate and return a rate coefficient model for the reverse reaction. 
        Currently this only works if the `kinetics` attribute is one of several
        (but not necessarily all) kinetics types.
        """
        cython.declare(Tlist=numpy.ndarray, Plist=numpy.ndarray, K=numpy.ndarray,
                       rxn=Reaction, klist=numpy.ndarray, i=cython.size_t,
                       Tindex=cython.size_t, Pindex=cython.size_t)

        supported_types = (
                            KineticsData.__name__,
                            Arrhenius.__name__,
                            MultiArrhenius.__name__,
                            PDepArrhenius.__name__,
                            MultiPDepArrhenius.__name__,
                            Chebyshev.__name__,
                            ThirdBody.__name__,
                            Lindemann.__name__,
                            Troe.__name__,
                            )

        # Get the units for the reverse rate coefficient
        kunits = getRateCoefficientUnitsFromReactionOrder(len(self.products))
            
        kf = self.kinetics
        if isinstance(kf, KineticsData):
            
            Tlist = kf.Tdata.value_si
            klist = numpy.zeros_like(Tlist)
            for i in range(len(Tlist)):
                klist[i] = kf.getRateCoefficient(Tlist[i]) / self.getEquilibriumConstant(Tlist[i])
            
            kr = KineticsData(Tdata=(Tlist,"K"), kdata=(klist,kunits), Tmin=(numpy.min(Tlist),"K"), Tmax=(numpy.max(Tlist),"K"))
            return kr
            
        elif isinstance(kf, Arrhenius):
            return self.reverseThisArrheniusRate(kf, kunits)

        elif network_kinetics and self.network_kinetics is not None:
            kf = self.network_kinetics
            return self.reverseThisArrheniusRate(kf, kunits)
                    
        elif isinstance (kf, Chebyshev):
            Tlist = 1.0/numpy.linspace(1.0/kf.Tmax.value, 1.0/kf.Tmin.value, 50)
            Plist = numpy.linspace(kf.Pmin.value, kf.Pmax.value, 20)
            K = numpy.zeros((len(Tlist), len(Plist)), numpy.float64)
            for Tindex, T in enumerate(Tlist):
                for Pindex, P in enumerate(Plist):
                    K[Tindex, Pindex] = kf.getRateCoefficient(T, P) / self.getEquilibriumConstant(T)
            kr = Chebyshev()
            kr.fitToData(Tlist, Plist, K, kunits, kf.degreeT, kf.degreeP, kf.Tmin.value, kf.Tmax.value, kf.Pmin.value, kf.Pmax.value)
            return kr
        
        elif isinstance(kf, PDepArrhenius):  
            if kf.Tmin is not None and kf.Tmax is not None:
                Tlist = 1.0/numpy.linspace(1.0/kf.Tmax.value, 1.0/kf.Tmin.value, 50)
            else:
                Tlist = 1.0/numpy.arange(0.0005, 0.0035, 0.0001)
            Plist = kf.pressures.value_si
            K = numpy.zeros((len(Tlist), len(Plist)), numpy.float64)
            for Tindex, T in enumerate(Tlist):
                for Pindex, P in enumerate(Plist):
                    K[Tindex, Pindex] = kf.getRateCoefficient(T, P) / self.getEquilibriumConstant(T)
            kr = PDepArrhenius()
            kr.fitToData(Tlist, Plist, K, kunits, kf.arrhenius[0].T0.value)
            return kr       
        
        elif isinstance(kf, MultiArrhenius):
            kr = MultiArrhenius()
            kr.arrhenius = []            
            rxn = Reaction(reactants = self.reactants, products = self.products)            
            for kinetics in kf.arrhenius:
                rxn.kinetics = kinetics
                kr.arrhenius.append(rxn.generateReverseRateCoefficient())
            return kr
        
        elif isinstance(kf, MultiPDepArrhenius):
            kr = MultiPDepArrhenius()              
            kr.arrhenius = []                
            rxn = Reaction(reactants = self.reactants, products = self.products)            
            for kinetics in kf.arrhenius:
                rxn.kinetics = kinetics
                kr.arrhenius.append(rxn.generateReverseRateCoefficient())
            return kr

        elif isinstance(kf, ThirdBody):
            lowPkunits = getRateCoefficientUnitsFromReactionOrder(len(self.products) + 1)
            krLow = self.reverseThisArrheniusRate(kf.arrheniusLow, lowPkunits)
            parameters = kf.__reduce__()[1]  # use the pickle helper to get all the other things needed
            kr = ThirdBody(krLow, *parameters[1:])
            return kr

        elif isinstance(kf, Lindemann):
            krHigh = self.reverseThisArrheniusRate(kf.arrheniusHigh, kunits)
            lowPkunits = getRateCoefficientUnitsFromReactionOrder(len(self.products) + 1)
            krLow = self.reverseThisArrheniusRate(kf.arrheniusLow, lowPkunits)
            parameters = kf.__reduce__()[1]  # use the pickle helper to get all the other things needed
            kr = Lindemann(krHigh, krLow, *parameters[2:])
            return kr

        elif isinstance(kf, Troe):
            krHigh = self.reverseThisArrheniusRate(kf.arrheniusHigh, kunits)
            lowPkunits = getRateCoefficientUnitsFromReactionOrder(len(self.products) + 1)
            krLow = self.reverseThisArrheniusRate(kf.arrheniusLow, lowPkunits)
            parameters = kf.__reduce__()[1]  # use the pickle helper to get all the other things needed
            kr = Troe(krHigh, krLow, *parameters[2:])
            return kr
        else:
            raise ReactionError(("Unexpected kinetics type {0}; should be one of {1}").format(self.kinetics.__class__, supported_types))

    def calculateTSTRateCoefficients(self, Tlist):
        return numpy.array([self.calculateTSTRateCoefficient(T) for T in Tlist], numpy.float64)

    def calculateTSTRateCoefficient(self, T):
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
            Qreac *= spec.getPartitionFunction(T) / (constants.R * T / 101325.)
            E0 -= spec.conformer.E0.value_si
        logging.debug('    Calculating Partition function for ' + self.transitionState.label)
        Qts = self.transitionState.getPartitionFunction(T) / (constants.R * T / 101325.)
        E0 += self.transitionState.conformer.E0.value_si
        k = (constants.kB * T / constants.h * Qts / Qreac) * math.exp(-E0 / constants.R / T)
        
        # Apply tunneling correction
        k *= self.transitionState.calculateTunnelingFactor(T)
        
        return k
        
    def canTST(self):
        """
        Return ``True`` if the necessary parameters are available for using
        transition state theory -- or the microcanonical equivalent, RRKM
        theory -- to compute the rate coefficient for this reaction, or
        ``False`` otherwise.
        """
        return len(self.transitionState.conformer.modes) > 0

    def calculateMicrocanonicalRateCoefficient(self, Elist, Jlist, reacDensStates, prodDensStates=None, T=0.0):
        """
        Calculate the microcanonical rate coefficient :math:`k(E)` for the reaction
        `reaction` at the energies `Elist` in J/mol. `reacDensStates` and 
        `prodDensStates` are the densities of states of the reactant and product
        configurations for this reaction. If the reaction is irreversible, only the
        reactant density of states is required; if the reaction is reversible, then
        both are required. This function will try to use the best method that it
        can based on the input data available:
        
        * If detailed information has been provided for the transition state (i.e.
          the molecular degrees of freedom), then RRKM theory will be used.
        
        * If the above is not possible but high-pressure limit kinetics
          :math:`k_\\infty(T)` have been provided, then the inverse Laplace 
          transform method will be used.
    
        The density of states for the product `prodDensStates` and the temperature
        of interest `T` in K can also be provided. For isomerization and association
        reactions `prodDensStates` is required; for dissociation reactions it is
        optional. The temperature is used if provided in the detailed balance
        expression to determine the reverse kinetics, and in certain cases in the
        inverse Laplace transform method.
        """
        return calculateMicrocanonicalRateCoefficient(self, Elist, Jlist, reacDensStates, prodDensStates, T)
    
    def isBalanced(self):
        """
        Return ``True`` if the reaction has the same number of each atom on
        each side of the reaction equation, or ``False`` if not.
        """
        from rmgpy.molecule.element import elementList
        
        cython.declare(reactantElements=dict, productElements=dict, molecule=Molecule, atom=Atom, element=Element)
        
        reactantElements = {}; productElements = {}
        for element in elementList:
            reactantElements[element] = 0
            productElements[element] = 0
        
        for reactant in self.reactants:
            if isinstance(reactant, Species):
                molecule = reactant.molecule[0]
            elif isinstance(reactant, Molecule):
                molecule = reactant
            for atom in molecule.atoms:
                reactantElements[atom.element] += 1
        
        for product in self.products:
            if isinstance(product, Species):
                molecule = product.molecule[0]
            elif isinstance(product, Molecule):
                molecule = product
            for atom in molecule.atoms:
                productElements[atom.element] += 1
         
        for element in elementList:
            if reactantElements[element] != productElements[element]:
                return False
        
        return True
    
    def generatePairs(self):
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
        the number of heavy atoms (C/O/N/S at the moment). This should
        work most of the time, but a more rigorous algorithm may be needed for
        some cases.
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

            reactantCarbons   = [sum([1 for atom in reactant.molecule[0].atoms if atom.isCarbon()])   for reactant in reactants]
            productCarbons    = [sum([1 for atom in  product.molecule[0].atoms if atom.isCarbon()])   for product  in products ]
            reactantOxygens   = [sum([1 for atom in reactant.molecule[0].atoms if atom.isOxygen()])   for reactant in reactants]
            productOxygens    = [sum([1 for atom in  product.molecule[0].atoms if atom.isOxygen()])   for product  in products ]
            reactantNitrogens = [sum([1 for atom in reactant.molecule[0].atoms if atom.isNitrogen()]) for reactant in reactants]
            productNitrogens  = [sum([1 for atom in  product.molecule[0].atoms if atom.isNitrogen()]) for product  in products ]
            reactantSilicons  = [sum([1 for atom in reactant.molecule[0].atoms if atom.isSilicon()])  for reactant in reactants]
            productSilicons   = [sum([1 for atom in  product.molecule[0].atoms if atom.isSilicon()])  for product  in products ]
            reactantSulfurs   = [sum([1 for atom in reactant.molecule[0].atoms if atom.isSulfur()])   for reactant in reactants]
            productSulfurs    = [sum([1 for atom in  product.molecule[0].atoms if atom.isSulfur()])   for product  in products ]
            reactantChlorines = [sum([1 for atom in reactant.molecule[0].atoms if atom.isChlorine()]) for reactant in reactants]
            productChlorines  = [sum([1 for atom in  product.molecule[0].atoms if atom.isChlorine()]) for product  in products ]
            reactantIodines   = [sum([1 for atom in reactant.molecule[0].atoms if atom.isChlorine()]) for reactant in reactants]
            productIodines    = [sum([1 for atom in  product.molecule[0].atoms if atom.isChlorine()]) for product  in products ]
            
            # Sort the reactants and products by C/O/N/S numbers
            reactants = [(carbon, oxygen, nitrogen, silicon, sulfur, chlorine, iodine, reactant) for carbon, oxygen, nitrogen, silicon, sulfur, chlorine, iodine, reactant
                         in zip(reactantCarbons,reactantOxygens,reactantNitrogens,reactantSilicons,reactantSulfurs,reactantChlorines, reactantIodines, reactants)]
            reactants.sort()
            products = [(carbon, oxygen, nitrogen, silicon, sulfur, chlorine, iodine, product) for carbon, oxygen, nitrogen, silicon, sulfur, chlorine, iodine, product
                        in zip(productCarbons,productOxygens,productNitrogens,productSilicons,productSulfurs,productChlorines, productIodines, products)]
            products.sort()
            
            while len(reactants) > 1 and len(products) > 1:
                self.pairs.append((reactants[-1][-1], products[-1][-1]))
                reactants.pop()
                products.pop()
            for reactant in reactants:
                for product in products:
                    self.pairs.append((reactant[-1], product[-1]))
    
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
        format = os.path.splitext(path)[1].lower()[1:]
        ReactionDrawer().draw(self, format, path)
        
    def _repr_png_(self):
        """
        Return a png picture of the reaction, useful for ipython-qtconsole.
        """
        from rmgpy.molecule.draw import ReactionDrawer
        tempFileName = 'temp_reaction.png'
        ReactionDrawer().draw(self, 'png', tempFileName)
        png = open(tempFileName,'rb').read()
        os.unlink(tempFileName)
        return png
            
    # Build the transition state geometry
    def generate3dTS(self, reactants, products):
        """
        Generate the 3D structure of the transition state. Called from 
        model.generateKinetics().
        
        self.reactants is a list of reactants
        self.products is a list of products
        """
        import rdkit
        import rdkit.Chem
        import rdkit.Chem.AllChem
        import rdkit.Geometry
        
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
                    dirVec = [{} for k in range(len(neighbor))]
                    lenVec = [None]*len(neighbor)
                    for k in range(0, len(neighbor)):
                        newIdx = neighbor[k].GetIdx()
                        newPt = mol.rdMol.GetConformer(mol.rdMolConfId).GetAtomPosition(newIdx)
                        dirVec[k] = point.DirectionVector(newPt)
                        lenVec[k] = point.Distance(newPt)
                    xCoord = [None]*len(neighbor)
                    yCoord = [None]*len(neighbor)
                    zCoord = [None]*len(neighbor) 
                    for k in range(0, len(neighbor)):
                        xCoord[k] = dirVec[k].x*lenVec[k]
                        yCoord[k] = dirVec[k].y*lenVec[k]
                        zCoord[k] = dirVec[k].z*lenVec[k]
            reactionAxis = [sum(xCoord), sum(yCoord), sum(zCoord)]
            reactants[i].reactionAxis = reactionAxis
        
        for i in range(0, len(products)):
            mol = products[i].molecule[0]
            for j in range(0, mol.rdMol.GetNumAtoms()):
                if mol.rdMol.GetAtomWithIdx(j).GetNumRadicalElectrons():
                    point = mol.rdMol.GetConformer(mol.rdMolConfId).GetAtomPosition(j)
                    neighbor = mol.rdMol.GetAtomWithIdx(j).GetNeighbors()
                    dirVec = [{} for k in range(len(neighbor))]
                    lenVec = [None]*len(neighbor)
                    for k in range(0, len(neighbor)):
                        newIdx = neighbor[k].GetIdx()
                        newPt = mol.rdMol.GetConformer(mol.rdMolConfId).GetAtomPosition(newIdx)
                        dirVec[k] = point.DirectionVector(newPt)
                        lenVec[k] = point.Distance(newPt)
                    xCoord = [None]*len(neighbor)
                    yCoord = [None]*len(neighbor)
                    zCoord = [None]*len(neighbor) 
                    for k in range(0, len(neighbor)):
                        xCoord[k] = dirVec[k].x*lenVec[k]
                        yCoord[k] = dirVec[k].y*lenVec[k]
                        zCoord[k] = dirVec[k].z*lenVec[k]
            reactionAxis = [sum(xCoord), sum(yCoord), sum(zCoord)]
            products[i].reactionAxis = reactionAxis
            
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
        other.specificCollider = self.specificCollider
        other.kinetics = deepcopy(self.kinetics)
        other.network_kinetics = deepcopy(self.network_kinetics)
        other.reversible = self.reversible
        other.transitionState = deepcopy(self.transitionState)
        other.duplicate = self.duplicate
        other.pairs = deepcopy(self.pairs)
        other.allow_pdep_route = self.allow_pdep_route
        other.elementary_high_p = self.elementary_high_p
        other.comment = deepcopy(self.comment)
        
        return other

    def ensure_species(self, reactant_resonance=False, product_resonance=True):
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
        if self.isForward:
            ensure_species(self.reactants, resonance=reactant_resonance, keepIsomorphic=True)
            ensure_species(self.products, resonance=product_resonance, keepIsomorphic=True)
        else:
            ensure_species(self.reactants, resonance=product_resonance, keepIsomorphic=True)
            ensure_species(self.products, resonance=reactant_resonance, keepIsomorphic=True)

        # convert reaction.pairs object to species
        if self.pairs:
            new_pairs = []
            for reactant, product in self.pairs:
                new_pair = []
                for reactant0 in self.reactants:
                    if reactant0.isIsomorphic(reactant):
                        new_pair.append(reactant0)
                        break
                for product0 in self.products:
                    if product0.isIsomorphic(product):
                        new_pair.append(product0)
                        break
                new_pairs.append(new_pair)
            self.pairs = new_pairs

        try:
            self.reverse.ensure_species()
        except AttributeError:
            pass

def _isomorphicSpeciesList(list1, list2, checkIdentical=False, checkOnlyLabel = False):
    """
    This method compares whether lists of species or molecules are isomorphic
    or identical. It is used for the 'isIsomorphic' method of Reaction class.
    It likely can be useful elswehere as well:
        
        list1 - list of species/molecule objects of reaction1
        list2 - list of species/molecule objects of reaction2
        checkIdentical - if true, uses the 'isIdentical' comparison
                         if false, uses the 'isIsomorphic' comparison
        checkOnlyLabel - only look at species' labels, no isomorphism checks
                         
    Returns True if the lists are isomorphic/identical & false otherwise
    """

    def comparison_method(other1, other2, checkIdentical=checkIdentical, checkOnlyLabel=checkOnlyLabel):
        if checkOnlyLabel:
            return str(other1) == str(other2)
        elif checkIdentical:
            return other1.isIdentical(other2)
        else:
            return other1.isIsomorphic(other2)

    if len(list1) == len(list2) == 1:
        if comparison_method(list1[0], list2[0]):
            return True
    elif len(list1) == len(list2) == 2:
        if comparison_method(list1[0], list2[0]) \
                    and comparison_method(list1[1], list2[1]):
            return True
        elif comparison_method(list1[0], list2[1]) \
                    and comparison_method(list1[1], list2[0]):
            return True
    elif len(list1) == len(list2) == 3:
        if (    comparison_method(list1[0], list2[0]) and
                comparison_method(list1[1], list2[1]) and
                comparison_method(list1[2], list2[2]) ):
            return True
        elif (  comparison_method(list1[0], list2[0]) and
                comparison_method(list1[1], list2[2]) and
                comparison_method(list1[2], list2[1]) ):
            return True
        elif (  comparison_method(list1[0], list2[1]) and
                comparison_method(list1[1], list2[0]) and
                comparison_method(list1[2], list2[2]) ):
            return True
        elif (  comparison_method(list1[0], list2[2]) and
                comparison_method(list1[1], list2[0]) and
                comparison_method(list1[2], list2[1]) ):
            return True
        elif (  comparison_method(list1[0], list2[1]) and
                comparison_method(list1[1], list2[2]) and
                comparison_method(list1[2], list2[0]) ):
            return True
        elif (  comparison_method(list1[0], list2[2]) and
                comparison_method(list1[1], list2[1]) and
                comparison_method(list1[2], list2[0]) ):
            return True
    elif len(list1) == len(list2):
        raise NotImplementedError("Can't check isomorphism of lists with {0} species/molecules".format(len(list1)))
    # nothing found
    return False
