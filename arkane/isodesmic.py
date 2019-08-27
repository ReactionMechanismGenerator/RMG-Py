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
This module provides the :class:`ErrorCancelingScheme` and related classes for the automatic generation of error
canceling reactions (e.g. isodesmic reactions). This code is heavily based on algorithms and ideas found in the existing
literature, including the following:

Buerger, P., Akroyd, J., Mosbach, S., & Kraft, M. (2018). A systematic method to estimate and validate enthalpies of
formation using error-cancelling balanced reactions. Combustion and Flame (Vol. 187).
https://doi.org/10.1016/j.combustflame.2017.08.013

Dobek, F. J., Ranasinghe, D. S., Throssell, K., & Petersson, G. A. (2013). Evaluation of the heats of formation of
corannulene and C60 by means of inexpensive theoretical procedures. Journal of Physical Chemistry A, 117(22), 4726–4730.
https://doi.org/10.1021/jp404158v
"""

from __future__ import division

import signal
from collections import deque,defaultdict,OrderedDict

from lpsolve55 import lpsolve, EQ, LE
import numpy as np
import pyomo.environ as pyo

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.quantity import ScalarQuantity
from rmgpy.molecule.graph import getVertexConnectivityValue
import logging

class ErrorCancelingSpecies(object):
    """Class for target and known (benchmark) species participating in an error canceling reaction"""

    def __init__(self, molecule, low_level_hf298, model_chemistry, high_level_hf298=None, fod=None, source=None):
        """

        Args:
            molecule (rmgpy.molecule.Molecule): molecule object to represent the species
            low_level_hf298 (tuple): (Hf298, unit) evaluated using a lower level of theory (e.g. DFT)
            model_chemistry (str): Level of theory used to calculate the low level thermo
            high_level_hf298 (tuple, optional): (Hf298 , unit) evaluated using a high level of theory
                (e.g. expt. data) that is serving as the "reference" for the isodesmic calculation
            source (str): Literature source from which the high level data was taken
        """
        if isinstance(molecule, Molecule):
            self.molecule = molecule
        else:
            raise ValueError('ErrorCancelingSpecies molecule attribute must be an rmgpy Molecule object. Instead a '
                             '{0} object was given'.format(type(molecule)))

        if isinstance(model_chemistry, str):
            self.model_chemistry = model_chemistry
        else:
            raise ValueError('The model chemistry string used to calculate the low level Hf298 must be provided '
                             'consistency checks. Instead, a {0} object was given'.format(type(model_chemistry)))

        self.low_level_hf298 = ScalarQuantity(*low_level_hf298)
        self.fod = fod

        # If the species is a reference species, then the high level data is already known
        self.high_level_hf298 = ScalarQuantity(*high_level_hf298) if high_level_hf298 else None
        self.source = source

    def __repr__(self):
        return '<ErrorCancelingSpecies {0}>'.format(self.molecule.toSMILES())


class ErrorCancelingReaction(object):
    """Class for representing an error canceling reaction, with the target species being an implicit reactant"""

    def __init__(self, target, species):
        """
        Initialize an error canceling reaction from ErrorCancelingSpecies objects

        The reaction must be written with the target species participating as a reactant with stoichiometric coefficient
        v=-1

        The species dictionary should be given as species/stochiometric coefficient pairs in the following format:
        {ErrorCancelingSpecies_1: v_1, ErrorCancelingSpecies_2: v_2, ... }

        Args:
            target (ErrorCancelingSpecies): high level H_f(298 K) will be estimated for this species
            species (dict): species taking place in the reaction (excluding the target)
        """

        self.target = target
        self.model_chemistry = self.target.model_chemistry

        # Perform a consistency check that all species are using the same model chemistry
        for spcs in species.keys():
            if spcs.model_chemistry != self.model_chemistry:
                raise ValueError('Species {0} has model chemistry {1}, which does not match the model chemistry of the '
                                 'reaction of {2}'.format(spcs, spcs.model_chemistry, self.model_chemistry))

        # Does not include the target, which is handled separately.
        self.species = species

        self.fod = self.calculate_fod()
        self.species_dict = self.get_dict()

    def __repr__(self):
        reactant_string = '1.0*{0} + '.format(self.target.molecule.toSMILES())
        product_string = ''
        for spcs, coeff in self.species.items():
            if coeff > 0:
                product_string += '{0}*{1} + '.format(coeff, spcs.molecule.toSMILES())
            else:
                reactant_string += '{0}*{1} + '.format(-1*coeff, spcs.molecule.toSMILES())

        return '<ErrorCancelingReaction {0}== {1}>'.format(reactant_string[:-2], product_string[:-2])

    def get_dict(self):

        species_dict = {}
        for spcs, coeff in self.species.items():
            species_dict[spcs] = coeff

        return species_dict

    def calculate_target_thermo(self):
        """
        Estimate the high level thermochemistry for the target species using the error canceling scheme

        Returns:
            rmgpy.quantity.ScalarQuantity: Hf298 in 'J/mol' estimated for the target species

        """
        low_level_h_rxn = sum(map(lambda spec: spec[0].low_level_hf298.value_si*spec[1], self.species.items())) - \
            self.target.low_level_hf298.value_si

        target_thermo = sum(map(lambda spec: spec[0].high_level_hf298.value_si*spec[1], self.species.items())) - \
            low_level_h_rxn
        return ScalarQuantity(target_thermo, 'J/mol')

    def calculate_fod(self):
        """
        Returns the sum of the fod number for all the species in the reaction.
        If a species is in the reaction does not have an fod number, it will return None.
        """

        fods = []

        if self.target.fod:
            fod_sum = self.target.fod
            fods.append(-1.0*self.target.fod)
        else:
            return None

        for spcs,v in self.species.items():
            if spcs.fod:
                fod_sum += spcs.fod * abs(v)
                fods.append(v*spcs.fod)
            else:
                return None

        fod_diff = sum(fods)

        return fod_sum

class SpeciesConstraints(object):
    """
    A class for defining and enumerating constraints to ReferenceSpecies objects for error canceling reactions
    """

    def __init__(self, target=None, reference_list=None, constraint_classes=None, conserve_bonds=True, conserve_ring_size=True):
        """
        Define the constraints that will be enforced, and determine the mapping of indices in the constraint vector to
        the labels for these constraints.

        To reduce the size of the linear programming problem that will try to find error canceling reactions of the
        target and subsets of the reference species, the `reference_species` list is automatically pruned to remove
        species that have additional atom, bond, and/or ring attributes not found in the target molecule.

        Args:
            target (ErrorCancelingSpecies): Molecule object for the target of the error canceling reaction scheme
            reference_list(:obj:`list` of :obj:`ErrorCancelingSpecies`): A list of molecule objects for the reference
                species that can participate in the error canceling reaction scheme
            conserve_bonds (bool, optional): Enforce the number of each bond type be conserved
            conserve_ring_size (bool, optional): Enforce that the number of each ring size be conserved
        """

        self.target = target
        self.reference_species = reference_list

        all_constraint_classes = ['class_0', 'class_1', 'class_2', 'class_3', 'class_4', 'class_5']

        if constraint_classes:
            if isinstance(constraint_classes,str):
                constraint_classes = constraint_classes.lower()
                if constraint_classes.lower() == 'all':
                    self.constraint_classes = all_constraint_classes
                else:
                    assert(constraint_classes in all_constraint_classes)
                    self.constraint_classes = [constraint_classes]
            elif isinstance(constraint_classes,list):
                self.constraint_classes = []
                for c in constraint_classes:
                    c = c.lower()
                    assert(c in constraint_classes)
                    self.constraint_classes.append(c)
            else:
                raise ValueError('constraint class must be None, all, or a <str> or <list> of'
                'supported constraint classes {}'.format(all_constraint_claseses))
        else:
            self.constraint_classes = all_constraint_classes
        
        self.conserve_bonds = conserve_bonds
        self.conserve_ring_size = conserve_ring_size
        self.target_constraint_classes = None

        if self.reference_species:
            self.constraint_map, self.descriptor_weights, self.fod_vector = self.get_constraint_map()
        else:
            self.constraint_map = None

    def get_descriptors(self, species):
        """

        """
        if isinstance(species, ErrorCancelingSpecies):
            mol = species.molecule
        elif isinstance(species, Molecule):
            mol = species
        elif isinstance(species,str): # assume smiles string
            mol = Molecule(SMILES=species)
        else:
            raise ValueError('species must be a {},{},or a SMILES string, not {}'.format(type(ErrorCancelingSpecies),
            type(Molecule), type(species)))
        
        descriptors = defaultdict(list)

        for atom in mol.atoms:
            
            atom_general = (atom.number, atom.radicalElectrons)
            atom_specific = (atom.number, atom.lonePairs, atom.charge, atom.radicalElectrons)

            descriptors['class_0'].append(atom_general)
            descriptors['class_1'].append(atom_specific)
            descriptors['class_2'].append(atom_specific)
            descriptors['class_4'].append(atom_specific)

            bonds_general = [] # for constraint class 3
            bonds_specific = [] # for constraint class 5
            radical_count = atom.radicalElectrons
            for bonded_atom, bond in atom.bonds.items():
                radical_count += bonded_atom.radicalElectrons
                order_number = bond.getOrderNum()
                bonds_general.append((bonded_atom.number,order_number))
                bonds_specific.append((bonded_atom.number,bonded_atom.radicalElectrons,bonded_atom.lonePairs,bonded_atom.charge,order_number))
            bonds_general.sort()
            bonds_specific.sort()
            class_3 = atom_specific + tuple(bonds_general) + (radical_count,)
            class_5 = atom_specific + tuple(bonds_specific) + (radical_count,)
            descriptors['class_3'].append(class_3)
            descriptors['class_5'].append(class_5)

        for bond in mol.getAllEdges():
            atom1 = bond.atom1
            atom2 = bond.atom2
            radical_count = atom1.radicalElectrons + atom2.radicalElectrons
            class_2 = tuple(sorted([atom1.number, atom2.number]) + [bond.getOrderNum()]) + (radical_count,)
            class_4 = tuple(sorted((a.number, a.radicalElectrons, a.lonePairs, a.charge) for a in (atom1, atom2)) + [bond.getOrderNum()]) + (radical_count,)
            descriptors['class_2'].append(class_2)
            descriptors['class_4'].append(class_4)

        if self.conserve_ring_size:
            rings = mol.getSmallestSetOfSmallestRings()
            if len(rings) > 0:
                for ring in rings:
                    for key in descriptors.keys():
                        descriptors[key].append(len(ring))

        return descriptors
    

    def get_constraint_map(self):
        # Enumerate all of the constraints in the target molecule to initialize the constraint mapping

        # constraint_map = {label: i for i, label in enumerate(self.target.molecule.get_element_count().keys())}

        # if self.conserve_bonds:
        #     j = len(constraint_map)
        #     constraint_map.update(
        #         {label: j + i for i, label in enumerate(self.target.molecule.enumerate_bonds().keys())})
        # if self.conserve_ring_size:
        #     j = len(constraint_map)
        #     possible_rings_sizes = set(map(lambda x: '{0}_ring'.format(len(x)),
        #                                    self.target.molecule.getSmallestSetOfSmallestRings()))
        #     constraint_map.update({label: j + i for i, label in enumerate(possible_rings_sizes)})

        # constraint_map = []
        # all_species = self.all_reference_species + [self.target]
        # for spcs in all_species:
        #     for atom in spcs.molecule.atoms:
        #         descriptor = atom.get_descriptor()
        #         if descriptor not in constraint_map:
        #             constraint_map.append(descriptor)

        constraint_map = defaultdict(defaultdict)
        objective_vectors = defaultdict(np.array)
        fod_vector = np.zeros(len(self.reference_species))
        
        # constraint_map['class_0'] = defaultdict()
        # constraint_map['class_1'] = defaultdict()
        # constraint_map['class_2'] = defaultdict()
        # constraint_map['class_3'] = defaultdict()
        # constraint_map['class_4'] = defaultdict()
        # constraint_map['class_5'] = defaultdict()
        # constraint_map['class_connectivity'] = defaultdict()

        #all_species = self.all_reference_species + [self.target]

        for i,species in enumerate(self.reference_species):
            if species.fod is None:
                fod = 0.01
            else:
                fod = species.fod
            fod_vector[i] = fod
            descriptors = self.get_descriptors(species)
            for constraint_class, descriptor_list in descriptors.items():
                for d in descriptor_list:
                    if d not in constraint_map[constraint_class].keys():
                        constraint_map[constraint_class][d] = 1
                    else:
                        constraint_map[constraint_class][d] += 1


        for constraint_class, descriptors in constraint_map.items():
            objective_vector = np.zeros(shape=(1,len(descriptors.keys())))
            for i,d in enumerate(descriptors.keys()):
                if isinstance(d,int): # descriptor is ring
                    weight = 1.5
                else:
                    radical_electrons = d[-1]
                    weight = float(radical_electrons) + 1.0
                objective_vector[0][i] = weight
            objective_vectors[constraint_class] = objective_vector

        fod_vector.shape = (len(self.reference_species),1)

        return constraint_map, objective_vectors, fod_vector

    def filter_constraint_classes(self,target=None):

        if not target: # use self.target
            if self.target:
                target = self.target
    
        # if not self.constraint_map and len(self.reference_species) > 0:
        #     self.constraint_map, self.descriptor_weights, self.fod_vector = self.get_constraint_map()
        # else:
        #     raise ValueError('There are no reference species for this SpeciesConstraints object from which' 
        #     'to generate a constraint map')

        target_constraint_classes = []
        unique_descriptors = defaultdict(list)

        descriptors = self.get_descriptors(target)
        for constraint_class, descriptor_list in descriptors.items():
            if constraint_class not in self.constraint_classes:
                continue
            use_constraint_class = True
            for d in descriptor_list:
                if d not in self.constraint_map[constraint_class].keys():
                    use_constraint_class = False
                    unique_descriptors[constraint_class].append(d)
            if use_constraint_class:
                target_constraint_classes.append(constraint_class)

        # for constraint_class in self.constraint_classes:
        #     if constraint_class not in unique_descriptors.keys():
        #         target_constraint_classes.append(constraint_class)

        target_constraint_classes.sort()
        self.target_constraint_classes = target_constraint_classes

        return target_constraint_classes, unique_descriptors


    def _enumerate_constraints(self, species):
        """
        Determine the constraint vector for a molecule given the enforced constraints

        Args:
            species (ErrorCancelingSpecies): Species whose constraints are to be enumerated

        Returns:
            np.ndarray: vector of the number of instances of each constraining feature e.g. number of carbon atoms
        """
        # constraint_vector = np.zeros(len(self.constraint_map))
        # molecule = species.molecule

        # try:
        #     atoms = molecule.get_element_count()
        #     for atom_label, count in atoms.items():
        #         constraint_vector[self.constraint_map[atom_label]] += count

        #     if self.conserve_bonds:
        #         bonds = molecule.enumerate_bonds()
        #         for bond_label, count in bonds.items():
        #             constraint_vector[self.constraint_map[bond_label]] += count

        # if self.conserve_ring_size:
        #     rings = molecule.getSmallestSetOfSmallestRings()
        #     if len(rings) > 0:
        #         for ring in rings:
        #             constraint_vector[self.constraint_map[constraint_class].keys().index(descriptor)] += 1
        # except KeyError:  # This molecule has a feature not found in the target molecule. Return None to exclude this
        #     return None

        constraint_vectors = dict()

        if not self.target_constraint_classes:
            self.filter_constraint_classes()

        for constraint_class in self.target_constraint_classes:
            constraint_vector = np.zeros(len(self.constraint_map[constraint_class]))
            descriptors = self.get_descriptors(species)
            for descriptor in descriptors[constraint_class]:
                constraint_vector[self.constraint_map[constraint_class].keys().index(descriptor)] += 1
            constraint_vectors[constraint_class] = constraint_vector
            
        # for atom in molecule.atoms:
        #     descriptor = list(atom.get_descriptor())
        #     descriptor.pop(1)
        #     bonds = []
        #     for atom,bond in atom.bonds.items():
        #         descriptor2 = list(atom.get_descriptor())
        #         descriptor2.pop(1)
        #         descriptor2.append(bond.getOrderNum())
        #         #bonds.append((atom.number,bond.getOrderNum()))
        #         bonds.append(descriptor2)
        #     bonds.sort()
        #     descriptor.extend(bonds)
        #     constraint_vector[self.constraint_map.index(descriptor)] += 1

        # else:
        #     for atom in molecule.atoms:
        #         descriptor = list(atom.get_descriptor())
        #         descriptor.pop(1)
        #         constraint_vector[self.constraint_map.index(descriptor)] += 1
        #     for bond,count in molecule.enumerate_bonds().items():
        #         constraint_vector[self.constraint_map.index(bond)] += count

        return constraint_vectors

    def calculate_constraints(self, constraint_class= None):
        """
        Calculate the constraint vector for the target and the constraint matrix for all allowable reference species

        Returns:
            np.ndarray: target constraint vector (1xn_constraints)
            np.ndarray: constraint matrix for allowable reference species (len(self.reference_species)xn_constraints)
        """

        if not self.target_constraint_classes:
            self.filter_constraint_classes()

        if len(self.target_constraint_classes) == 0:
            logging.info('target species has descriptors that do not match with any of the reference species descriptors'
            'for all of the constraint classes')
            return None, None, None
  
        #constraint_class = self.constraint_classes[-1]

        #target_constraints = self._enumerate_constraints(self.target,constraint_class)
        target_constraints_dict = self._enumerate_constraints(self.target)
        #constraint_matrix = []
        constraint_matrix_dict = defaultdict(list)
        for spcs in self.reference_species:
            spcs_constraints_dict = self._enumerate_constraints(spcs)
            for constraint_class, constraint_vector in spcs_constraints_dict.items():
                constraint_matrix_dict[constraint_class].append(constraint_vector)
            #self.reference_species.append(spcs)
            #constraint_matrix.append(spcs_constraints)

        return target_constraints_dict, constraint_matrix_dict
        # return target_constraints, np.array(constraint_matrix, dtype=int), self.descriptor_weights[constraint_class]


class ErrorCancelingScheme(object):
    """
    A Base class for calculating target species thermochemistry using error canceling reactions
    """

    def __init__(self, target, reference_set,constraint_classes=None, conserve_bonds=True, conserve_ring_size=True):
        """

        Args:
            target (ErrorCancelingSpecies): Species whose Hf298 will be calculated using error canceling reactions
            reference_set (:obj:`list` of :obj:`ErrorCancelingSpecies`): list of reference species that can participate
                in error canceling reactions to help calculate the thermochemistry of the target
            conserve_bonds (bool): Flag to determine if the number and type of each bond must be conserved in each error
                canceling reaction
            conserve_ring_size (bool): Flag to determine if the number of each ring size must be conserved in each error
                canceling reaction
        """

        self.target = target
        self.reference_set = reference_set
        #self.constraint_classes = constraint_classes

        self.constraints = SpeciesConstraints(target, reference_set, constraint_classes = constraint_classes, conserve_bonds = conserve_bonds,
                                              conserve_ring_size=conserve_ring_size)

        self.reference_species = self.constraints.reference_species

        try:
            self.target_constraints, self.constraint_matrix = self.constraints.calculate_constraints()
        except:
            self.target_constraints = None
            self.constraint_matrix = None
        # except:
        #     self.target_constraint = self.constraint_matrix = None

    def _find_error_canceling_reaction(self, reference_subset, constraint_class, milp_software='lpsolve'):
        """
        Automatically find a valid error canceling reaction given a subset of the available benchmark species. This
        is done by solving a mixed integer linear programming (MILP) problem similiar to
        Buerger et al. (https://doi.org/10.1016/j.combustflame.2017.08.013)

        Args:
            reference_subset (list): A list of indices from self.reference_species that can participate in the reaction
            milp_software (str, optional): 'lpsolve' (default) or 'pyomo'. lpsolve is usually faster.

        Returns:
            ErrorCancelingReaction: reaction with the target species (if a valid reaction is found, else `None`)
            np.ndarray: indices (of the subset) for the species that participated in the return reaction
        """

        if constraint_class is None:
            constraint_class = self.constraints.target_constraint_classes[-1]
        
        #class_penalty = constraint_class_penalty[constraint_class]
        constraint_matrix = np.array(self.constraint_matrix[constraint_class], dtype=int)
        #weights = self.constraints.descriptor_weights[constraint_class]

        # Define the constraints based on the provided subset
      
        c_matrix = np.take(constraint_matrix, reference_subset, axis=0)
        c_matrix = np.tile(c_matrix, (2, 1))
        #fod_vector = self.constraints.fod_vector[reference_subset]
        #fod_vector = np.tile(fod_vector,(2,1))
        # objective_fn: vector, one per reference species (doubled) 
        # objective_fn = np.sum(c_matrix, 1, dtype=int) # indicating ~size of molecule
        # objective_fn = np.sum(c_matrix * fod_vector * weights * class_penalty, 1, dtype=float)
        objective_matrix = np.ones_like(c_matrix,dtype=int) # all species are equal
        #fod_penalty = 5.0 # arbitrary chosen penalty for high FOD
        # objective_fn = 1.0 + fod_penalty*fod_vector # penalize things with higher FOD
        #objective_fn = np.sum(objective_matrix * (2.0*fod_vector) * weights * class_penalty,1,dtype=float)
        objective_fn = np.sum(objective_matrix,1,dtype=float)
        targets = -1 * self.target_constraints[constraint_class]
        m, n  = c_matrix.shape
        split = int(m/2)
        solution = None
        
        if milp_software == 'pyomo':
            """This pyomo code not maintained (by David)"""
            # Setup the MILP problem using pyomo
            lp_model = pyo.ConcreteModel()
            lp_model.i = pyo.RangeSet(0, m - 1)
            lp_model.j = pyo.RangeSet(0, n - 1)
            lp_model.r = pyo.RangeSet(0, split-1)  # indices before the split correspond to reactants
            lp_model.p = pyo.RangeSet(split, m - 1)  # indices after the split correspond to products
            lp_model.v = pyo.Var(lp_model.i, domain=pyo.NonNegativeIntegers)  # The stoich. coef. we are solving for
            lp_model.c = pyo.Param(lp_model.i, lp_model.j, initialize=lambda _, i, j: c_matrix[i, j])
            lp_model.s = pyo.Param(lp_model.i, initialize=lambda _, i: sum_constraints[i])
            lp_model.t = pyo.Param(lp_model.j, initialize=lambda _, j: targets[j])

            def obj_expression(model):
                return pyo.summation(model.v, model.s, index=model.i)

            lp_model.obj = pyo.Objective(rule=obj_expression)

            def constraint_rule(model, j):
                return sum(model.v[i] * model.c[i, j] for i in model.r) - \
                       sum(model.v[i] * model.c[i, j] for i in model.p) == model.t[j]

            lp_model.constraints = pyo.Constraint(lp_model.j, rule=constraint_rule)

            # Solve the MILP problem using the CBC MILP solver (https://www.coin-or.org/Cbc/)
            opt = pyo.SolverFactory('glpk')
            results = opt.solve(lp_model)

            # Return None if a valid reaction is not found
            if results.solver.status != pyo.SolverStatus.ok:
                return None, None

            # Extract the solution and find the species with non-zero stoichiometric coefficients
            solution = lp_model.v.extract_values().values()

        elif milp_software == 'lpsolve':
            # Save the current signal handler
            sig = signal.getsignal(signal.SIGINT)

            # Setup the MILP problem using lpsolve
            lp = lpsolve('make_lp', 0, m)
            lpsolve('set_verbose', lp, 2)  # Reduce the logging from lpsolve
            #lpsolve('set_obj_fn', lp, sum_constraints)
            lpsolve('set_obj_fn', lp, objective_fn)
            lpsolve('set_minim', lp)


            for j in range(n):
                lpsolve('add_constraint', lp, np.concatenate((c_matrix[:split, j], -1*c_matrix[split:, j])), EQ,
                        targets[j])

            lpsolve('add_constraint', lp, np.ones(m), LE, 20)  # Use at most 20 species (including replicates)
            lpsolve('set_timeout', lp, 1)  # Move on if lpsolve can't find a solution quickly

            # Constrain v_i to be 4 or less
            for i in range(m):
                lpsolve('set_upbo', lp, i, 4)

            # All v_i must be integers
            lpsolve('set_int', lp, [True]*m)

            status = lpsolve('solve', lp)

            # Reset signal handling since lpsolve changed it
            try:
                signal.signal(signal.SIGINT, sig)
            except ValueError:
                # This is not being run in the main thread, so we cannot reset signal
                pass

            if status != 0:
                return None, None, None

            else:
                obj, solution = lpsolve('get_solution', lp)[:2]

        reaction = ErrorCancelingReaction(self.target, dict())
        subset_indices = []
        for index, v in enumerate(solution):
            if v > 0:
                subset_indices.append(index % split)
                if index < split:
                    reaction.species.update({self.reference_species[reference_subset[index]]: -v})
                else:
                    reaction.species.update({self.reference_species[reference_subset[index % split]]: v})

        return reaction, np.array(subset_indices), obj

    def multiple_error_canceling_reaction_search(self, reactions=None, rejected_reactions=None, n_reactions_max=20, reject_reactions = False, milp_software='lpsolve'):
        """
        Generate multiple error canceling reactions involving the target and a subset of the reference species.

        To do this, a rudimentary search is implemented whereby all possible combinations of the species participating
        in the previously found reaction are excluded from the reference species subset for the next generation process.
        This is implemented using a FIFO queue structure.

        Args:
            n_reactions_max (int, optional): The maximum number of found reactions that will returned, after which no
                further searching will occur even if there are possible subsets left in the queue.
            milp_software (str, optional): 'lpsolve' (default) or 'pyomo'. lpsolve is usually faster.

        Returns:
            :obj:list of :obj:ErrorCancelingReaction: A list of the found error canceling reactions
        """
        
        if (self.constraint_matrix is None) or (len(self.constraints.target_constraint_classes) == 0):
            return None,None
        
        if not reactions:
            reactions = OrderedDict()
        if not rejected_reactions:
            rejected_reactions = OrderedDict()

        current_constraint_class = self.constraints.target_constraint_classes[-1]
        full_set = np.arange(0, len(self.reference_species))
        if reject_reactions and len(self.constraints.target_constraint_classes) > 1:
            _, _, lower_class_obj = self._find_error_canceling_reaction(full_set,self.constraints.target_constraint_classes[-2])
        else:
            lower_class_obj = 1e6
        subset_queue = deque()
        subset_queue.append(full_set)
        max_attempts = 1000
        attempts = 0
        rejected = 0

        while (len(subset_queue) != 0) and (len(reactions) < n_reactions_max) and (attempts<max_attempts):
            subset = subset_queue.popleft()
            if len(subset) == 0:
                continue
            reaction, subset_indices, obj = self._find_error_canceling_reaction(subset,current_constraint_class, milp_software=milp_software)
            attempts += 1

            if reaction is None:
                continue

            if len(reactions) == 0:
                h298 = reaction.calculate_target_thermo()
                fod = reaction.calculate_fod()
                ref_obj = obj
                reactions[(reaction,current_constraint_class)] = (ref_obj,fod,h298,1.0)
                for index in subset_indices:
                    subset_queue.append((np.delete(subset, index)))
                continue
            
            # elif current_constraint_class not in [key[1] for key in reactions.keys()]:
            #     ref_obj = obj
            
            unique = True
            for rxn in [r[0] for r in reactions.keys() + rejected_reactions.keys()]:
                if rxn.species == reaction.species:
                    unique = False
                    for index in subset_indices:
                        subset_queue.append((np.delete(subset, index)))
                    continue
           
            if unique:
                h298 = reaction.calculate_target_thermo()
                fod = reaction.calculate_fod()
                
                weight = 1.0
                # if ref_obj <= obj:
                #     weight = 1-0.5*((obj-ref_obj)/(obj+ref_obj))
                # else:
                #     weight = 1.0
                #     ref_obj = obj

                #if (len(reactions) <= 10) or (not reject_reactions):
                
                reactions[(reaction,current_constraint_class)] = (obj,fod,h298,weight)
                for index in subset_indices:
                    subset_queue.append((np.delete(subset, index)))
                continue
                
                # if reject_reactions:
                #     data = np.array(reactions.values(),ndmin=2)
                #     sum_of_weights = np.sum(data[:,-1])

                #     h298_sum = np.sum(np.array([h.value_si for h in data[:,-2]]) * data[:,-1])
                #     h298_mean = h298_sum/sum_of_weights
                #     new_h298_mean = (h298.value_si*weight + h298_sum)/(sum_of_weights + weight)
                #     if (abs(new_h298_mean-h298_mean) >= abs(2 * np.std([h.value_si for h in data[:,-2]]))) or (obj>lower_class_obj):
                #         rejected += 1
                #         rejected_reactions[(reaction,current_constraint_class)] = (obj,fod,h298,weight)
                #         for index in subset_indices:
                #             subset_queue.append((np.delete(subset, index)))
                #     else:
                #         reactions[(reaction,current_constraint_class)] = (obj,fod,h298,weight)
                #         for index in subset_indices:
                #             subset_queue.append((np.delete(subset, index)))


        self.constraints.target_constraint_classes.pop()
        if len(reactions)<n_reactions_max and len(self.constraints.constraint_classes) > 0:
            #self.target_constraint, self.constraint_matrix, self.descriptor_weights = self.constraints.calculate_constraints()
            self.multiple_error_canceling_reaction_search(reactions,rejected_reactions)

        return reactions,rejected_reactions

    def calculate_target_enthalpy(self, n_reactions_max=20, milp_software='lpsolve'):
        """
        Perform a multiple error canceling reactions search and calculate hf298 for the target species by taking the
        median hf298 value from among the error canceling reactions found

        Args:
            n_reactions_max (int, optional): The maximum number of found reactions that will returned, after which no
                further searching will occur even if there are possible subsets left in the queue.
            milp_software (str, optional): 'lpsolve' (default) or 'pyomo'. lpsolve is usually faster.

        Returns:
            ScalarQuantity: Standard heat of formation at 298 K calculated for the target species
            list: reaction list containing all error canceling reactions found

        """
        if self.constraint_matrix is None:
            return None,None,None
        reactions, rejected_reactions = self.multiple_error_canceling_reaction_search(n_reactions_max=n_reactions_max, milp_software=milp_software)
        if reactions is None:
            return None, None, None
        data = np.array(reactions.values())
        if len(reactions) == 0:
            return None,None,None
        elif len(reactions) == 1:
            h298_mean = float(data[0][-2].value_si)
        else:
            try:
                sum_of_weights = np.sum(data[:,-1])
                h298_sum = np.sum(np.array([h.value_si for h in data[:,-2]]) * data[:,-1])
                h298_mean = h298_sum/sum_of_weights
            except:
                print "something went wrong when calculating H298 for {}".format(self.target)
                logging.info("something went wrong when calculating H298 for {}".format(self.target))
                logging.info(data)
                return None,None,None
        # h298_list = np.zeros(len(reaction_list))
        # h298_kcal_mol = np.zeros(len(reaction_list))
        # fod_dict = dict()

        # for i, rxn in enumerate(reaction_list):
        #     h298 = rxn.calculate_target_thermo()
        #     h298_list[i] = h298.value_si
        #     h298_kcal_mol[i] = h298.value_si/h298.conversionFactors['kcal/mol']
        #     fod_dict[rxn] = (h298_kcal_mol[i],rxn.target.fod) + tuple(s.fod for s in rxn.species)

        # return ScalarQuantity(np.median(h298_list), 'J/mol'), zip(reaction_list,h298_kcal_mol), fod_dict
        return ScalarQuantity(h298_mean, 'J/mol'), reactions, rejected_reactions

class IsodesmicScheme(ErrorCancelingScheme):
    """
    An error canceling reaction where the number and type of both atoms and bonds are conserved
    """
    def __init__(self, target, reference_set):
        super(IsodesmicScheme, self).__init__(target, reference_set, conserve_bonds=True, conserve_ring_size=False)


class IsodesmicRingScheme(ErrorCancelingScheme):
    """
    A stricter form of the traditional isodesmic reaction scheme where the number of each ring size is also conserved
    """
    def __init__(self, target, reference_set):
        super(IsodesmicRingScheme, self).__init__(target, reference_set, conserve_bonds=True, conserve_ring_size=True)


if __name__ == '__main__':
    pass
