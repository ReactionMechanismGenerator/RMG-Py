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

"""
import itertools
import logging
import math
import os.path
import re
import numpy as np
from copy import deepcopy
from CoolProp.CoolProp import PropsSI

import rmgpy.constants as constants
from rmgpy.data.base import Database, Entry, make_logic_node, DatabaseError
from rmgpy.molecule import Molecule, Group, ATOMTYPES
from rmgpy.species import Species
from rmgpy.exceptions import InputError
from rmgpy.data.thermo import is_aromatic_ring, is_bicyclic, find_aromatic_bonds_from_sub_molecule, \
    convert_ring_to_sub_molecule, is_ring_partial_matched, bicyclic_decomposition_for_polyring, \
    split_bicyclic_into_single_rings, saturate_ring_bonds


################################################################################

def save_entry(f, entry):
    """
    Write a Pythonic string representation of the given `entry` in the solvation
    database to the file object `f`.
    """
    f.write('entry(\n')
    f.write('    index = {0:d},\n'.format(entry.index))
    f.write('    label = "{0}",\n'.format(entry.label))

    if isinstance(entry.item, list):
        if len(entry.item) == 1:
            item = entry.item[0]
            if isinstance(item, Species):
                if Molecule(smiles=item.molecule[0].to_smiles()).is_isomorphic(item.molecule[0]):
                    # The SMILES representation accurately describes the molecule, so we can save it that way.
                    f.write('    molecule = "{0}",\n'.format(item.molecule[0].to_smiles()))
                else:
                    f.write('    molecule = \n')
                    f.write('"""\n')
                    f.write(item.to_adjacency_list(remove_h=False))
                    f.write('""",\n')
            else:
                raise DatabaseError("Not sure how to save {0!r}".format(entry.item))
        else:
            f.write('    molecule = [')
            for i in range(len(entry.item)):
                item = entry.item[i]
                if i > 0:
                    f.write(', ')
                if isinstance(item, Species):
                    if Molecule(smiles=item.molecule[0].to_smiles()).is_isomorphic(item.molecule[0]):
                        f.write('"{0}"'.format(item.molecule[0].to_smiles()))
                    else:
                        f.write('\n"""\n')
                        f.write(item.molecule[0].to_adjacency_list())
                        f.write('"""')
                else:
                    raise DatabaseError("Not sure how to save {0!r}".format(entry.item))
            f.write('],\n')
    elif isinstance(entry.item, Species):
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
        f.write('        name_in_coolprop = "{0!s}",\n'.format(entry.data.name_in_coolprop))
        f.write('    ),\n')
    elif entry.data is None:
        f.write('    solute = None,\n')
    else:
        raise DatabaseError("Not sure how to save {0!r}".format(entry.data))

    if isinstance(entry.data_count, DataCountGAV):
        f.write('    dataCount = DataCountGAV(\n')
        f.write('        S = {0!r},\n'.format(entry.data_count.S))
        f.write('        B = {0!r},\n'.format(entry.data_count.B))
        f.write('        E = {0!r},\n'.format(entry.data_count.E))
        f.write('        L = {0!r},\n'.format(entry.data_count.L))
        f.write('        A = {0!r},\n'.format(entry.data_count.A))
        f.write('    ),\n')
    elif isinstance(entry.data_count, DataCountSolvent):
        f.write('    dataCount = DataCountSolvent(\n')
        f.write('        dGsolvCount = {0!r},\n'.format(entry.data_count.dGsolvCount))
        f.write('        dGsolvMAE = {0!r},\n'.format(entry.data_count.dGsolvMAE))
        f.write('        dHsolvCount = {0!r},\n'.format(entry.data_count.dHsolvCount))
        f.write('        dHsolvMAE = {0!r},\n'.format(entry.data_count.dHsolvMAE))
        f.write('    ),\n')
    elif entry.data_count is None:
        # only write entry.data_count if this is not a solute library entry
        if isinstance(entry.item, Group) or isinstance(entry.data, SolventData):
            f.write('    dataCount = None,\n')
    else:
        raise DatabaseError("Not sure how to save {0!r}".format(entry.data_count))

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


def add_solute_data(solute_data1, solute_data2, group_additivity=False, verbose=False):
    """
    Add the solute data `solute_data2` to the data `solute_data1`,
    and return `solute_data1`.

    If `group_additivity` is True, append comments related to group additivity estimation
    If `verbose` is False, omit the comments from a "zero entry", whose E, S, A, B, amd L are all 0.
    If `verbose` is True, or solute_data2 is not a zero entry, add solute_data2.comment to solute_data1.comment.
    """
    solute_data1.A += solute_data2.A
    solute_data1.B += solute_data2.B
    solute_data1.L += solute_data2.L
    solute_data1.E += solute_data2.E
    solute_data1.S += solute_data2.S

    test_zero = sum(abs(value) for value in
                    [solute_data2.A, solute_data2.B, solute_data2.L, solute_data2.E, solute_data2.S])
    # Used to check if all of the entries in solute_data2 are zero

    if group_additivity:
        if verbose or test_zero != 0:
            # If verbose==True or test_zero!=0, add solute_data2.comment to solute_data1.comment.
            if solute_data1.comment:
                solute_data1.comment += ' + {0}'.format(solute_data2.comment)
            else:
                solute_data1.comment = solute_data2.comment

    return solute_data1


def remove_solute_data(solute_data1, solute_data2, group_additivity=False, verbose=False):
    """
    Remove the solute data `solute_data2` from the data `solute_data1`,
    and return `solute_data1`.
    If `verbose` is True, append ' - solute_data2.comment' to the solute_data1.comment.
    If `verbose` is False, remove the solute_data2.comment from the solute_data1.comment.
    """
    solute_data1.A -= solute_data2.A
    solute_data1.B -= solute_data2.B
    solute_data1.L -= solute_data2.L
    solute_data1.E -= solute_data2.E
    solute_data1.S -= solute_data2.S

    if group_additivity:
        if verbose:
            solute_data1.comment += ' - {0}'.format(solute_data2.comment)
        else:
            solute_data1.comment = re.sub(re.escape(' + ' + solute_data2.comment), '', solute_data1.comment, 1)
    return solute_data1


def average_solute_data(solute_data_list):
    """
    Average a list of solute_data values together.
    """

    num_values = len(solute_data_list)

    if num_values == 0:
        raise ValueError('No solute data values were inputted to be averaged.')
    else:
        logging.debug('Averaging solute data over {0} value(s).'.format(num_values))

        if num_values == 1:
            return deepcopy(solute_data_list[0])

        else:
            averaged_solute_data = deepcopy(solute_data_list[0])
            for solute_data in solute_data_list[1:]:
                averaged_solute_data = add_solute_data(averaged_solute_data, solute_data)
            averaged_solute_data.S /= num_values
            averaged_solute_data.B /= num_values
            averaged_solute_data.E /= num_values
            averaged_solute_data.L /= num_values
            averaged_solute_data.A /= num_values
            return averaged_solute_data


def get_critical_temperature(compound_name):
    """Returns the critical temperature of a compound in K for the given compound_name.

    The critical temperature is given by the CoolProp function, PropsSI

    Args:
        compound_name (str): a name of the compound used in CoolProp.

    Returns:
        Tc (float): critical temperature of the given compound in K.

    Raises:
        DatabaseError: If the given compound_name is not available in CoolProp.

    """

    try:
        Tc = PropsSI('T_critical', compound_name)
    except:
        raise DatabaseError(f"Critical temperature is not available for the given compound: {compound_name}")
    return Tc


def get_liquid_saturation_density(compound_name, temp):
    """Returns the liquid-phase saturation density of a compound in mol/m^3 at an input temperature.

    The density is calculated by the CoolProp function, PropsSI.

    Args:
        compound_name (str): a name of the compound used in CoolProp.
        temp (float): Temperature [K] at which the liquid-phase saturation density is calculated.

    Returns:
        rho_l (float): liquid-phase saturation density at the given temperature in mol/m^3.

    Raises:
        DatabaseError: If the given compound_name is not available in CoolProp or the given temperature is out of
        the calculable range.

    """

    try:
        rho_l = PropsSI('Dmolar', 'T', temp, 'Q', 0, compound_name)  # saturated liquid phase density in mol/m^3
    except:
        raise DatabaseError(f"Liquid-phase saturation density is not available for {compound_name} at {temp} K.")
    return rho_l


def get_gas_saturation_density(compound_name, temp):
    """Returns the gas-phase saturation density of a compound in mol/m^3 at an input temperature.

    The density is calculated by the CoolProp function, PropsSI.

    Args:
        compound_name (str): a name of the compound used in CoolProp.
        temp (float): Temperature [K] at which the liquid-phase saturation density is calculated.

    Returns:
        rho_g (float): gas-phase saturation density at the given temperature in mol/m^3.

    Raises:
        DatabaseError: If the given compound_name is not available in CoolProp or the given temperature is out of
        the calculable range.

    """

    try:
        rho_g = PropsSI('Dmolar', 'T', temp, 'Q', 1, compound_name)  # saturated gas phase density in mol/m^3
    except:
        raise DatabaseError(f"Gas-phase saturation density is not available for {compound_name} at {temp} K.")
    return rho_g

################################################################################

class SolventData(object):
    """
    Stores Abraham/Mintz parameters for characterizing a solvent.
    """

    def __init__(self, s_h=None, b_h=None, e_h=None, l_h=None, a_h=None,
                 c_h=None, s_g=None, b_g=None, e_g=None, l_g=None, a_g=None, c_g=None, A=None, B=None,
                 C=None, D=None, E=None, alpha=None, beta=None, eps=None, name_in_coolprop=None):
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
        # This corresponds to the solvent's name in CoolProp. CoolProp is an external package used for
        # fluid property calculation. If the solvent is not available in CoolProp, this is set to None
        self.name_in_coolprop = name_in_coolprop

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


    def get_solvent_saturation_pressure(self, T):
        """
        Returns the saturation pressure of the solvent in Pa if the solvent is available in CoolProp (i.e. name_in_coolprop is
        not None); raises DatabaseError if the solvent is not available in CoolProp (i.e. name_in_coolprop is None).
        """
        if self.name_in_coolprop is not None:
            return PropsSI('P', 'T', T, 'Q', 0, self.name_in_coolprop)
        else:
            raise DatabaseError("Saturation pressure is not available for the solvent whose `name_in_coolprop` is None")

    def get_solvent_density(self, T):
        """
        Returns the density of the solvent in Pa if the solvent is available in CoolProp (i.e. name_in_coolprop is
        not None); raises DatabaseError if the solvent is not available in CoolProp (i.e. name_in_coolprop is None).
        """
        if self.name_in_coolprop is not None:
            return PropsSI('D', 'T', T, 'Q', 0, self.name_in_coolprop)
        else:
            raise DatabaseError("Saturation pressure is not available for the solvent whose `name_in_coolprop` is None")

    def get_solvent_coolprop_Tmin(self):
        """
        Returns the minimum temperature range available for the solvent in  CoolProp
        """ 
        if self.name_in_coolprop is not None:
            return PropsSI('Tmin', self.name_in_coolprop)
        else:
            raise DatabaseError("Tmin is not available for the solvent whose `name_in_coolprop` is None")

    def get_solvent_coolprop_Tcrit(self):
        """
        Returns the critical temperature of solvent from CoolProp
        """ 
        if self.name_in_coolprop is not None:
            return PropsSI('Tcrit', self.name_in_coolprop)
        else:
            raise DatabaseError("Critical temperature is not available for the solvent whose `name_in_coolprop` is None")

class SolvationCorrection(object):
    """
    Stores corrections for enthalpy, entropy, and Gibbs free energy when a species is solvated.
    Enthalpy and Gibbs free energy is in J/mol; entropy is in J/mol/K
    """

    def __init__(self, enthalpy=None, gibbs=None, entropy=None):
        self.enthalpy = enthalpy
        self.entropy = entropy
        self.gibbs = gibbs

class KfactorParameters(object):
    """
    Stores 4 coefficients (A, B, C, D) in the following K-factor relationships and
    the transition temperature (T_transition) in K.
    1) T <= T_transition : Harvey's semi-empirical relationship
                    Tr*ln(K-factor) = A + B(1-Tr)^0.355 + Cexp(1-Tr)(Tr)^0.59
    2) T_transition <= T < T_c : Japas and Levelt Sengers' linear relationship
                    Tr*ln(K-factor) = D(rho_l / rho_c -1)
    Relevant definitions:
    rho_l = saturated liquid phase density of the solvent [=] mol / m^3
    rho_c = critical density of the solvent [=] mol / m^3
    Tr = reduced temperature = T / Tc [=] K
    Tc = critical temperature of the solvent [=] K
    K-factor = y_solute / x_solute.
    y_solute = mole fraction of the solute in a gas phase at equilibrium in a binary dilute mixture
    x_solute = mole fraction of the solute in a liquid phase at equilibrium in a binary dilute mixture
    """
    def __init__(self, A=None, B=None, C=None, D=None, T_transition=None):
        self.lower_T = [A, B, C]
        self.higher_T = D
        self.T_transition = T_transition # in K

class SoluteData(object):
    """
    Stores Abraham parameters to characterize a solute
    """
    # Set class variable with McGowan volumes
    mcgowan_volumes = {
        1: 8.71, 2: 6.75,
        6: 16.35, 7: 14.39, 8: 12.43, 9: 10.47, 10: 8.51,
        14: 26.83, 15: 24.87, 16: 22.91, 17: 20.95, 18: 18.99,
        35: 26.21, 53: 34.53,
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
        if molecule.contains_surface_site():
            molecule = molecule.get_desorbed_molecules()[0]
            molecule.saturate_unfilled_valence()
        Vtot = 0.0

        for atom in molecule.atoms:
            try:
                Vtot += self.mcgowan_volumes[atom.element.number]
            except KeyError:
                raise DatabaseError('McGowan volume not available for element {}'.format(atom.element.nubmer))

            # divide contribution in half since all bonds would be counted twice this way
            Vtot -= len(molecule.get_bonds(atom)) * 6.56 / 2

        self.V = Vtot / 100  # division by 100 to get units correct.

class DataCountGAV(object):
    """
    A class for storing the number of data used to fit each solute parameter group value in the solute group additivity.

    ...

    Attributes
    ----------
    S : int
        The number of data used to fit S value of the solute GAV group.
    B : int
        The number of data used to fit B value of the solute GAV group.
    E : int
        The number of data used to fit E value of the solute GAV group.
    L : int
        The number of data used to fit L value of the solute GAV group.
    A : int
        The number of data used to fit A value of the solute GAV group.
    comment: str
        Comment with extra information.
    """

    def __init__(self, S=None, B=None, E=None, L=None, A=None, comment=""):
        self.S = S
        self.B = B
        self.E = E
        self.L = L
        self.A = A
        self.comment = comment

    def __repr__(self):
        return "DataCountGAV(S={0},B={1},E={2},L={3},A={4},comment={5!r})".format(
            self.S, self.B, self.E, self.L, self.A, self.comment)

class DataCountSolvent(object):
        """
        A class for storing the number of data used to fit the solvent parameters and the solvation energy
        and enthalpy mean absolute error (MAE) associated with the fitted solvent parameters.

        ...

        Attributes
        ----------
        dGsolvCount : int
            The number of solvation free energy data used to fit Abraham solvent parameters (s_g, b_g, e_g, l_g, a_g, c_g).
        dGsolvMAE : (float, str)
            The solvation free energy mean absolute error associated with the fitted Abraham solvent parameters.
            The second element is the unit of the solvation free energy error.
        dHsolvCount : int
            The number of solvation enthalpy data used to fit Mintz solvent parameters (s_h, b_h, e_h, l_h, a_h, c_h).
        dHsolvMAE : (float, str)
            The solvation enthalpy mean absolute error associated with the fitted Mintz solvent parameters.
            The second element is the unit of the solvation enthalpy error.
        comment: str
            Comment with extra information.
        """

        def __init__(self, dGsolvCount=None, dGsolvMAE=None, dHsolvCount=None, dHsolvMAE=None, comment=""):
            self.dGsolvCount = dGsolvCount
            self.dGsolvMAE = dGsolvMAE
            self.dHsolvCount = dHsolvCount
            self.dHsolvMAE = dHsolvMAE
            self.comment = comment

        def __repr__(self):
            return "DataCountSolvent(dGsolvCount={0},dGsolvMAE={1!r},dHsolvCount={2},dHsolvMAE={3!r},comment={4!r})".format(
                self.dGsolvCount, self.dGsolvMAE, self.dHsolvCount, self.dHsolvMAE, self.comment)

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
                   dataCount=None,
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
            data_count=dataCount,
            reference=reference,
            reference_type=referenceType,
            short_desc=shortDesc,
            long_desc=longDesc.strip(),
        )

    def load(self, path):
        """
        Load the solvent library from the given path
        """
        Database.load(self, path, local_context={'SolventData': SolventData, 'DataCountSolvent': DataCountSolvent}, global_context={})

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

    def get_solvent_data_count(self, label):
        """
        Get a solvent's data count information from its name
        """
        return self.entries[label].data_count

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
                   dataCount=None,
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
            data_count=dataCount,
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
            'DataCountGAV': DataCountGAV,
            'SolventData': SolventData,
            'DataCountSolvent': DataCountSolvent
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

    def get_all_solvent_data(self, solvent_species):
        """Return all possible sets of solvent data for a given :class:`Species` object `species`.

        This searches for the solvent data from the solvent library. The mixture solvent data containing
        the given species are also searched.

        Args:
            solvent_species: :class:`Species` object

        Returns:
            solvent_data_list: A list of tuple. The tuple is (solvent_label, solvent_entry) for the matched solvent.
            `solvent_label` is the label of the solvent and `solvent_entry` is a :class:`Entry` object containing
            solvent information stored in the solvent library of RMG-database.

        """

        solvent_data_list = []

        # Searches for the solvent data from library
        for label, value in self.libraries['solvent'].entries.items():
            spc_list = value.item
            # This also returns a mixture solvent entry if any of the mixture solvents matches the given solvent species
            for spc in spc_list:
                if solvent_species.is_isomorphic(spc):
                    solvent_data_list.append((label, value))
                    break
        return solvent_data_list

    def find_solvent_from_smiles(self, smiles):
        """
        This function uses the given SMILES string to find matching solvents from the solvent library. It uses
        the isomorphism check, which does not take stereochemistry into account, so it's possible that more than
        one solvent is matched from the solvent library. For example, for a SMILES "Cl/C=C/Cl", both
        "cis-1,2-dichloroethene" and "trans-1,2-dichloroethene" will be found as matching solvents.
        If any solvent is matched, it returns a list of a tuple containing the solvent label and solvent entry.
        If no solvent is matched, it returns an empty list.
        This method only works for pure solvent search, and it will not work for mixture solvents.

        Parameters
        ----------
        smiles : str
            A SMILES string of a compound to search for

        Returns
        -------
        match_list: list of tuple
            A list of tuple. The tuple is (solvent_label, solvent_entry) for the matched solvent. `solvent_label` is
            the label of the solvent and `solvent_entry` is a :class:`Entry` object containing solvent information
            stored in the solvent library of RMG-database.

        """

        solvent_spc = Species().from_smiles(smiles)
        solvent_spc.generate_resonance_structures()
        return self.get_all_solvent_data(solvent_spc)

    def get_solvent_data(self, solvent_name):
        try:
            solvent_data = self.libraries['solvent'].get_solvent_data(solvent_name)
        except:
            raise DatabaseError('Solvent {0!r} not found in database'.format(solvent_name))
        return solvent_data

    def get_solvent_data_count(self, solvent_name):
        try:
            solvent_data_count = self.libraries['solvent'].get_solvent_data_count(solvent_name)
        except:
            raise DatabaseError('Solvent {0!r} not found in database'.format(solvent_name))
        return solvent_data_count

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
        """
        logging.info('Loading solvation thermodynamics group database from {0}...'.format(path))
        categories = [
            'radical',
            'group',
            'ring',
            'polycyclic',
            'longDistanceInteraction_cyclic',
            'longDistanceInteraction_noncyclic',
            'halogen'
        ]
        self.groups = {
            category: SoluteGroups(label=category).load(os.path.join(path, category + '.py'),
                                                        self.local_context, self.global_context)
            for category in categories
        }

        self.record_ring_generic_nodes()
        self.record_polycylic_generic_nodes()

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

    def record_polycylic_generic_nodes(self):
        """
        Identify generic nodes in tree for polycyclic groups.
        Saves them as a list in the `generic_nodes` attribute
        in the polycyclic :class:`SoluteGroups` object, which
        must be pre-loaded.

        Necessary for polycyclic heuristic.
        """
        self.groups['polycyclic'].generic_nodes = ['PolycyclicRing']
        for label, entry in self.groups['polycyclic'].entries.items():
            if isinstance(entry.data, SoluteData):
                continue
            self.groups['polycyclic'].generic_nodes.append(label)

    def record_ring_generic_nodes(self):
        """
        Identify generic nodes in tree for ring groups.
        Saves them as a list in the `generic_nodes` attribute
        in the ring :class:`SoluteGroups` object, which
        must be pre-loaded.

        Necessary for polycyclic heuristic.
        """
        self.groups['ring'].generic_nodes = ['Ring']
        for label, entry in self.groups['ring'].entries.items():
            if isinstance(entry.data, SoluteData):
                continue
            self.groups['ring'].generic_nodes.append(label)

    def get_solute_data(self, species, skip_library=False):
        """
        For a given :class:`Species` object `species`, this function first searches
        the loaded libraries in order, returning the first match found, before falling back to
        estimation via solute data group additivity.

        Args:
            species: :class:`Species` object for which the solute data is searched or estimated.
            skip_library: a Boolean to skip the library search and do the group additivity estimation.

        Returns:
            solute_data: :class:`SoluteData` object for the given species object.

        """

        # Check the library first
        if species.molecule[0].contains_surface_site(): #desorb and saturate s
            molecule = species.molecule[0].get_desorbed_molecules()[0]
            molecule.saturate_unfilled_valence()
            species = Species(molecule=[molecule])

        solute_data = None
        if not skip_library:
            solute_data = self.get_solute_data_from_library(species, self.libraries['solute'])
        if solute_data is not None:
            assert len(solute_data) == 3, "solute_data should be a tuple (solute_data, library, entry)"
            solute_data[0].comment += "Solute library: " + solute_data[2].label
            solute_data = solute_data[0]
        else:
            # Solute not found in any loaded libraries, so estimate
            # First try finding stable species in libraries and then use HBI / halogen solvation correction
            molecule = species.molecule[0]
            if molecule.is_radical():
                # If the molecule is a radical, check if any of the saturated forms are in the libraries and
                # perform a HBI correction on the original molecule.
                # If the saturated form is not found in the libraries and the molecule has any halogen atoms, the
                # halogen atoms in the molecule are replaced by the hydrogen atoms and the method will check if
                # the saturated, replaced form is in the libraries. Then halogen corrections will be applied on the
                # saturated form, and lastly a HBI correction will be applied on the original molecule.
                molecule.clear_labeled_atoms()
                # First see if the saturated molecule is in the libraries.
                solute_data = self.estimate_radical_solute_data_via_hbi(molecule, self.get_solute_data_from_library)
            elif molecule.has_halogen():
                # If the molecule is halogenated, check if any of the replaced forms are in the libraries
                # first and perform halogen correction on them
                solute_data = self.estimate_halogen_solute_data(molecule, self.get_solute_data_from_library)

        if solute_data is None:
            # Solute or its saturated structure is not found in libraries. Use group additivty to determine solute data
            solute_data = self.get_solute_data_from_groups(species)

        # No group additivity for V, so set using atom sizes
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
            data[0].comment += "Solute library: " + data[2].label
            solute_data_list.append(data)
        # Estimate from group additivity
        # Make it a tuple
        data = (self.get_solute_data(species, skip_library=True), None, None)
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
        :class:`Species` object `species` by estimation using the group
        additivity method. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        
        It estimates the solute data for the first item in the species's
        molecule list because it is the most stable resonance structure found
        by gas-phase thermo estimate.
        """
        molecule = species.molecule[0]
        molecule.clear_labeled_atoms()
        molecule.update_atomtypes()
        solute_data = self.estimate_solute_via_group_additivity(molecule)

        return solute_data

    def estimate_solute_via_group_additivity(self, molecule):
        """
        Return the set of Abraham solute parameters corresponding to a given
        :class:`Molecule` object `molecule` by estimation using the group
        additivity values.

        Parameters
        ----------
        molecule : :class:`Molecule`

        Returns
        -------
        solute_data: :class:`SoluteData` object
            Contains the Abraham solute parameters estimated from the group additivity.
        """
        # For solute data estimation we need the atoms to already be sorted because we
        # iterate over them; if the order changes during the iteration then we
        # will probably not visit the right atoms, and so will get the solute data wrong
        molecule.sort_atoms()

        if molecule.is_radical():
            # If the molecule is a radical, calculate the solute data for the saturated form using
            # group additivity and then perform an HBI correction on the original molecule.
            # If the molecule has any halogen atoms, estimate_halogen_solute_data method is used for the saturated
            # form within the estimate_radical_solute_data_via_hbi method.
            molecule.clear_labeled_atoms()
            solute_data = self.estimate_radical_solute_data_via_hbi(molecule, self.compute_group_additivity_solute)
        elif molecule.has_halogen():
            # If the molecule is halogenated, compute the solute data for the replaced form whose halogen atoms are
            # replaced by hydrogen atoms first and perform halogen correction on it.
            solute_data = self.estimate_halogen_solute_data(molecule, self.compute_group_additivity_solute)
        else:
            solute_data = self.compute_group_additivity_solute(molecule)
        return solute_data

    def estimate_radical_solute_data_via_hbi(self, molecule, stable_solute_data_estimator):
        """
        Estimate the solute data of a radical by saturating it,
        applying the provided stable_solute_data_estimator method on the saturated species,
        then applying hydrogen bond increment corrections for the radical site(s).
        If the radical molecule contains halogens and the solute data of the saturated species
        are not found in the libraries, the solute data of the saturated species are found
        using estimate_halogen_solute_data method.

        Parameters
        ----------
        molecule : :class:`Molecule`
        stable_solute_data_estimator : function
            A solute_data estimator method for the saturated / halogen-replaced form.
            It's either `get_solute_data_from_library` or `compute_group_additivity_solute`.

        Returns
        -------
        solute_data: :class:`SoluteData` object
            Contains the Abraham solute parameters estimated from the group additivity and HBI correction.
        """
        if not molecule.is_radical():
            raise ValueError("Method only valid for radicals.")

        saturated_struct = molecule.copy(deep=True)
        added = saturated_struct.saturate_radicals()
        saturated_struct.props['saturated'] = True

        # Get solute data estimate for saturated form of structure
        if stable_solute_data_estimator == self.get_solute_data_from_library:
            # Get data from libraries
            saturated_spec = Species(molecule=[saturated_struct])
            solute_data_sat = stable_solute_data_estimator(saturated_spec, library=self.libraries['solute'])
            if solute_data_sat:
                if len(solute_data_sat) != 3:
                    raise ValueError("solute_data should be a tuple (solute_data, library, entry), "
                                       "not {0}".format(solute_data_sat))
                sat_label = solute_data_sat[2].label
                solute_data_sat = solute_data_sat[0]
                solute_data_sat.comment += "Solute library: " + sat_label
            else:
                # If the molecule has any halogen atoms, check whether the solute data for the saturated, replaced form
                # can be found from libraries.
                if molecule.has_halogen():
                    solute_data_sat = self.estimate_halogen_solute_data(saturated_struct, stable_solute_data_estimator)
        else:
            if molecule.has_halogen():
                solute_data_sat = self.estimate_halogen_solute_data(saturated_struct, stable_solute_data_estimator)
            else:
                solute_data_sat = stable_solute_data_estimator(saturated_struct)

        if solute_data_sat is None:
            # We couldn't get solute data for the saturated species from libraries.
            # However, if we were trying group additivity, this could be a problem
            if stable_solute_data_estimator == self.compute_group_additivity_solute:
                logging.info("Solute data of saturated {0} of molecule {1} is None.".format(saturated_struct, molecule))
            return None

        solute_data = solute_data_sat

        # For each radical site, get radical correction
        # Only one radical site should be considered at a time; all others
        # should be saturated with hydrogen atoms
        for atom in added:
            # Remove the added hydrogen atoms and bond and restore the radical
            for H, bond in added[atom]:
                saturated_struct.remove_bond(bond)
                saturated_struct.remove_atom(H)
                atom.increment_radical()
            saturated_struct.update()
            try:
                self._add_group_solute_data(solute_data, self.groups['radical'], saturated_struct, {'*': atom})
            except KeyError:
                pass
            # Re-saturate
            for H, bond in added[atom]:
                saturated_struct.add_atom(H)
                saturated_struct.add_bond(bond)
                atom.decrement_radical()

        # Remove all of the long distance interactions of the saturated structure. Then add the long interactions of the radical.
        # Take C1=CC=C([O])C(O)=C1 as an example, we need to remove the interation of OH-OH, then add the interaction of Oj-OH.
        # For now, we only apply this part to cyclic structure because we only have radical interaction data for aromatic radical.
        if saturated_struct.is_cyclic():
            sssr = saturated_struct.get_smallest_set_of_smallest_rings()
            for ring in sssr:
                for atomPair in itertools.permutations(ring, 2):
                    try:
                        self._remove_group_solute_data(solute_data, self.groups['longDistanceInteraction_cyclic'],
                                                       saturated_struct, {'*1': atomPair[0], '*2': atomPair[1]})
                    except KeyError:
                        pass
            sssr = molecule.get_smallest_set_of_smallest_rings()
            for ring in sssr:
                for atomPair in itertools.permutations(ring, 2):
                    try:
                        self._add_group_solute_data(solute_data, self.groups['longDistanceInteraction_cyclic'],
                                                    molecule,
                                                    {'*1': atomPair[0], '*2': atomPair[1]})
                    except KeyError:
                        pass

        # prevents the original species name being used for the HBI corrected radical in species generation
        solute_data.label = ''

        return solute_data

    def estimate_halogen_solute_data(self, molecule, stable_solute_data_estimator):
        """
        Estimate the solute data of a halogenated molecule by replacing halogens with hydrogen atoms,
        applying the provided stable_solute_data_estimator method on the replaced_struct species,
        then applying halogen corrections for the halogenated site(s).

        Parameters
        ----------
        molecule : :class:`Molecule`
        stable_solute_data_estimator : function
            A solute_data estimator method for the halogen-replaced form.
            It's either `get_solute_data_from_library` or `compute_group_additivity_solute`.

        Returns
        -------
        solute_data: :class:`SoluteData` object
            Contains the Abraham solute parameters estimated from the group additivity and halogen correction.
        """
        if not molecule.has_halogen():
            raise ValueError("Method only valid for halogenated molecule.")

        if molecule.is_radical():
            raise ValueError("Method only valid for non-radical molecule.")

        replaced_struct = molecule.copy(deep=True)
        replaced_struct.replace_halogen_with_hydrogen()

        # Get solute data estimate for replaced form of structure
        if stable_solute_data_estimator == self.get_solute_data_from_library:
            # Get data from libraries
            replaced_spec = Species(molecule=[replaced_struct])
            solute_data_replaced = stable_solute_data_estimator(replaced_spec, library=self.libraries['solute'])
            if solute_data_replaced:
                if len(solute_data_replaced) != 3:
                    raise ValueError("solute_data should be a tuple (solute_data, library, entry), "
                                       "not {0}".format(solute_data_replaced))
                replaced_label = solute_data_replaced[2].label
                solute_data_replaced = solute_data_replaced[0]
                solute_data_replaced.comment += "Solute library: " + replaced_label
        else:
            solute_data_replaced = stable_solute_data_estimator(replaced_struct)

        if solute_data_replaced is None:
            # We couldn't get solute data for the replaced species from libraries.
            # However, if we were trying group additivity, this could be a problem
            if stable_solute_data_estimator == self.compute_group_additivity_solute:
                logging.info("Solute data of halogen-replaced form {0} of molecule {1} is None.".format(replaced_struct, molecule))
            return None

        solute_data = solute_data_replaced

        # For each halogenated site, get halogen correction
        for atom in molecule.atoms:
            if atom.is_halogen():
                try:
                    self._add_group_solute_data(solute_data, self.groups['halogen'], molecule, {'*': atom})
                except KeyError:
                    pass

        # Remove all of the long distance interactions of the replaced structure. Then add the long interactions of the halogenated molecule.
        if replaced_struct.is_cyclic():
            sssr = replaced_struct.get_smallest_set_of_smallest_rings()
            for ring in sssr:
                for atomPair in itertools.permutations(ring, 2):
                    try:
                        self._remove_group_solute_data(solute_data, self.groups['longDistanceInteraction_cyclic'],
                                                       replaced_struct, {'*1': atomPair[0], '*2': atomPair[1]})
                    except KeyError:
                        pass
            sssr = molecule.get_smallest_set_of_smallest_rings()
            for ring in sssr:
                for atomPair in itertools.permutations(ring, 2):
                    try:
                        self._add_group_solute_data(solute_data, self.groups['longDistanceInteraction_cyclic'],
                                                    molecule,
                                                    {'*1': atomPair[0], '*2': atomPair[1]})
                    except KeyError:
                        pass

        # prevents the original species name being used for the HBI corrected radical in species generation
        solute_data.label = ''

        return solute_data

    def compute_group_additivity_solute(self, molecule):

        """
        Return the set of Abraham solute parameters corresponding to a given
        :class:`Molecule` object `molecule` by estimation using the group
        additivity values. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        """

        assert not molecule.is_radical(), "This method is only for saturated non-radical species."
        assert not molecule.has_halogen(), "This method is only for non-halogenated species."
        # For solute data estimation we need the atoms to already be sorted because we
        # iterate over them; if the order changes during the iteration then we
        # will probably not visit the right atoms, and so will get the solute data wrong
        molecule.sort_atoms()

        # Create the SoluteData object
        solute_data = SoluteData(
            S=0.0,
            B=0.0,
            E=0.0,
            L=0.0,
            A=0.0,
        )
        cyclic = molecule.is_cyclic()
        # Generate estimate of thermodynamics
        for atom in molecule.atoms:
            # Iterate over heavy (non-hydrogen) atoms
            if atom.is_non_hydrogen():
                # Get initial thermo estimate from main group database
                try:
                    self._add_group_solute_data(solute_data, self.groups['group'], molecule, {'*': atom})
                except KeyError:
                    logging.error("Couldn't find in main solute database:")
                    logging.error(molecule)
                    logging.error(molecule.to_adjacency_list())
                    raise

                # Correct for gauche and 1,5- interactions
                # Pair atom with its 1st and 2nd nonHydrogen neighbors,
                # Then match the pair with the entries in the database longDistanceInteraction_noncyclic.py
                # Currently we only have gauche(1,4) and 1,5 interactions in that file.
                # If you want to add more corrections for longer distance, please call get_nth_neighbor() method accordingly.
                # Potentially we could include other.py in this database, but it's a little confusing how to label atoms for the entries in other.py
                if not molecule.is_atom_in_cycle(atom):
                    for atom_2 in molecule.get_nth_neighbor([atom], [1, 2]):
                        if not molecule.is_atom_in_cycle(atom_2):
                            # This is the correction for noncyclic structure. If `atom` or `atom_2` is in a cycle, do not apply this correction.
                            # Note that previously we do not do gauche for cyclic molecule, which is unreasonable for cyclic molecule with a long tail.
                            try:
                                self._add_group_solute_data(solute_data,
                                                            self.groups['longDistanceInteraction_noncyclic'],
                                                            molecule, {'*1': atom, '*2': atom_2})
                            except KeyError:
                                pass

        # Do long distance interaction correction for cyclic molecule.
        # First get smallest set of smallest rings.
        # Then for every single ring, generate the atom pairs by itertools.permutation.
        # Finally match the atom pair with the database.
        # WIPWIPWIPWIPWIPWIPWIP         #########################################         WIPWIPWIPWIPWIPWIPWIP
        # WIP: For now, in the database, if an entry describes the interaction between same groups,
        # it will be halved because it will be counted twice here.
        # Alternatively we could keep all the entries as their full values by using combinations instead of permutations here.
        # In that case, we need to add more lines to match from reverse side when we didn't hit the most specific level from the forward side.
        # PS: by saying 'forward side', I mean {'*1':atomPair[0], '*2':atomPair[1]}. So the following is the reverse side '{'*1':atomPair[1], '*2':atomPair[0]}'
        # In my opinion, it's cleaner to do it in the current way.
        # WIPWIPWIPWIPWIPWIPWIP         #########################################         WIPWIPWIPWIPWIPWIPWIP
        if cyclic:
            sssr = molecule.get_smallest_set_of_smallest_rings()
            for ring in sssr:
                for atomPair in itertools.permutations(ring, 2):
                    try:
                        self._add_group_solute_data(solute_data, self.groups['longDistanceInteraction_cyclic'],
                                                    molecule,
                                                    {'*1': atomPair[0], '*2': atomPair[1]})
                    except KeyError:
                        pass

        # Do ring corrections separately because we only want to match
        # each ring one time
        if cyclic:
            monorings, polyrings = molecule.get_disparate_cycles()
            for ring in monorings:
                # Make a temporary structure containing only the atoms in the ring
                # NB. if any of the ring corrections depend on ligands not in the ring, they will not be found!
                try:
                    self._add_ring_correction_solute_data_from_tree(solute_data, self.groups['ring'], molecule, ring)
                except KeyError:
                    logging.error("Couldn't find a match in the monocyclic ring database even though "
                                  "monocyclic rings were found.")
                    logging.error(molecule)
                    logging.error(molecule.to_adjacency_list())
                    raise
            for polyring in polyrings:
                # Make a temporary structure containing only the atoms in the ring
                # NB. if any of the ring corrections depend on ligands not in the ring, they will not be found!
                try:
                    self._add_polycyclic_correction_solute_data(solute_data, molecule, polyring)
                except KeyError:
                    logging.error("Couldn't find a match in the polycyclic ring database even though "
                                  "polycyclic rings were found.")
                    logging.error(molecule)
                    logging.error(molecule.to_adjacency_list())
                    raise

        return solute_data

    def _add_polycyclic_correction_solute_data(self, solute_data, molecule, polyring):
        """
        INPUT: `polyring` as a list of `Atom` forming a polycyclic ring
        OUTPUT: if the input `polyring` can be fully matched in polycyclic database, the correction
        will be directly added to `solute_data`; otherwise, a heuristic approach will
        be applied.
        """
        # look up polycylic tree directly
        matched_group_solutedata, matched_group, is_partial_match = self._add_ring_correction_solute_data_from_tree(
            None, self.groups['polycyclic'], molecule, polyring)

        # if partial match (non-H atoms number same between
        # polycylic ring in molecule and match group)
        # otherwise, apply heuristic algorithm
        if not is_partial_match:
            if is_bicyclic(polyring) and matched_group.label in self.groups['polycyclic'].generic_nodes:
                # apply secondary decompostion formula
                # to get a estimated_group_thermodata
                estimated_bicyclic_solutedata = self.get_bicyclic_correction_solute_data_from_heuristic(polyring)
                if not estimated_bicyclic_solutedata:
                    estimated_bicyclic_solutedata = matched_group_solutedata
                solute_data = add_solute_data(solute_data, estimated_bicyclic_solutedata, group_additivity=True,
                                              verbose=True)
            else:
                # keep matched_group_solutedata as is
                solute_data = add_solute_data(solute_data, matched_group_solutedata, group_additivity=True,
                                              verbose=True)
                # By setting verbose=True, we turn on the comments of polycyclic correction to pass the unittest.
                # Typically this comment is very short and also very helpful to check if the ring correction is calculated correctly.
        else:
            self._add_poly_ring_correction_solute_data_from_heuristic(solute_data, polyring)

    def _add_poly_ring_correction_solute_data_from_heuristic(self, solute_data, polyring):
        """
        INPUT: `polyring` as a list of `Atom` forming a polycyclic ring, which can
        only be partially matched.
        OUTPUT: `polyring` will be decomposed into a combination of 2-ring polycyclics
        and each one will be looked up from polycyclic database. The heuristic formula
        is "polyring solute correction = sum of correction of all 2-ring sub-polycyclics -
        overlapped single-ring correction"; the calculated polyring solute correction
        will be finally added to input `solute_data`.
        """

        # polyring decomposition
        bicyclics_merged_from_ring_pair, ring_occurrences_dict = bicyclic_decomposition_for_polyring(polyring)

        # loop over 2-ring cores
        for bicyclic in bicyclics_merged_from_ring_pair:
            matched_group_solutedata, matched_group, _ = self._add_ring_correction_solute_data_from_tree(
                None, self.groups['polycyclic'], bicyclic, bicyclic.atoms)

            if matched_group.label in self.groups['polycyclic'].generic_nodes:
                # apply secondary decompostion formula
                # to get a estimated_group_solutedata
                estimated_bicyclic_solutedata = self.get_bicyclic_correction_solute_data_from_heuristic(bicyclic.atoms)
                if not estimated_bicyclic_solutedata:
                    estimated_bicyclic_solutedata = matched_group_solutedata
                solute_data = add_solute_data(solute_data, estimated_bicyclic_solutedata, group_additivity=True,
                                              verbose=True)
            else:
                # keep matched_group_solutedata as is
                solute_data = add_solute_data(solute_data, matched_group_solutedata, group_additivity=True,
                                              verbose=True)

        # loop over 1-ring
        for singleRingTuple, occurrence in ring_occurrences_dict.items():
            single_ring = list(singleRingTuple)

            if occurrence >= 2:
                submol, _ = convert_ring_to_sub_molecule(single_ring)

                if not is_aromatic_ring(submol):
                    aromatic_bonds = find_aromatic_bonds_from_sub_molecule(submol)
                    for aromaticBond in aromatic_bonds:
                        aromaticBond.set_order_num(1)

                    submol.saturate_unfilled_valence()
                    single_ring_solutedata = self._add_ring_correction_solute_data_from_tree(
                        None, self.groups['ring'], submol, submol.atoms)[0]

                else:
                    submol.update()
                    single_ring_solutedata = self._add_ring_correction_solute_data_from_tree(
                        None, self.groups['ring'], submol, submol.atoms)[0]
            for _ in range(occurrence - 1):
                solute_data = remove_solute_data(solute_data, single_ring_solutedata, True, True)
                # By setting verbose=True, we turn on the comments of polycyclic correction to pass the unittest.
                # Typically this comment is very short and also very helpful to check if the ring correction is calculated correctly.

    def get_bicyclic_correction_solute_data_from_heuristic(self, bicyclic):
        # saturate if the bicyclic has unsaturated bonds
        # otherwise return None
        bicyclic_submol = convert_ring_to_sub_molecule(bicyclic)[0]
        saturated_bicyclic_submol, already_saturated = saturate_ring_bonds(bicyclic_submol)

        if already_saturated:
            return None
        # split bicyclic into two single ring submols
        single_ring_submols = split_bicyclic_into_single_rings(bicyclic_submol)

        # split saturated bicyclic into two single ring submols
        saturated_single_ring_submols = split_bicyclic_into_single_rings(saturated_bicyclic_submol)

        # apply formula:
        # bicyclic correction ~= saturated bicyclic correction -
        # saturated single ring corrections + single ring corrections

        estimated_bicyclic_solute_data = SoluteData(
            S=0.0,
            B=0.0,
            E=0.0,
            L=0.0,
            A=0.0,
        )

        saturated_bicyclic_solute_data = self._add_ring_correction_solute_data_from_tree(
            None, self.groups['polycyclic'], saturated_bicyclic_submol, saturated_bicyclic_submol.atoms)[0]

        estimated_bicyclic_solute_data = add_solute_data(estimated_bicyclic_solute_data,
                                                         saturated_bicyclic_solute_data,
                                                         group_additivity=True)

        estimated_bicyclic_solute_data.comment = "Estimated bicyclic component: " + \
                                                 saturated_bicyclic_solute_data.comment

        for submol in saturated_single_ring_submols:

            if not is_aromatic_ring(submol):
                aromatic_bonds = find_aromatic_bonds_from_sub_molecule(submol)
                for aromatic_bond in aromatic_bonds:
                    aromatic_bond.set_order_num(1)

                submol.saturate_unfilled_valence()
                single_ring_solute_data = self._add_ring_correction_solute_data_from_tree(
                    None, self.groups['ring'], submol, submol.atoms)[0]

            else:
                submol.update()
                single_ring_solute_data = self._add_ring_correction_solute_data_from_tree(
                    None, self.groups['ring'], submol, submol.atoms)[0]
            estimated_bicyclic_solute_data = remove_solute_data(estimated_bicyclic_solute_data,
                                                                single_ring_solute_data,
                                                                group_additivity=True, verbose=True)

        for submol in single_ring_submols:

            if not is_aromatic_ring(submol):
                aromatic_bonds = find_aromatic_bonds_from_sub_molecule(submol)
                for aromatic_bond in aromatic_bonds:
                    aromatic_bond.set_order_num(1)

                submol.saturate_unfilled_valence()
                single_ring_solute_data = self._add_ring_correction_solute_data_from_tree(
                    None, self.groups['ring'], submol, submol.atoms)[0]

            else:
                submol.update()
                single_ring_solute_data = self._add_ring_correction_solute_data_from_tree(
                    None, self.groups['ring'], submol, submol.atoms)[0]

            estimated_bicyclic_solute_data = add_solute_data(estimated_bicyclic_solute_data,
                                                             single_ring_solute_data, group_additivity=True,
                                                             verbose=True)

        return estimated_bicyclic_solute_data

    def _add_ring_correction_solute_data_from_tree(self, solute_data, ring_database, molecule, ring):
        """
        Determine the ring correction group additivity solute data for the given
         `ring` in the `molecule`, and add it to the existing solute data
        `solute_data`.
        Also returns the matched ring group from the database from which the data originated.
        """
        matched_ring_entries = []
        # label each atom in the ring individually to try to match the group
        # for each ring, save only the ring that matches the most specific leaf in the tree.
        for atom in ring:
            atoms = {'*': atom}
            entry = ring_database.descend_tree(molecule, atoms)
            matched_ring_entries.append(entry)

        if matched_ring_entries is []:
            raise KeyError('Node not found in database.')
        # Decide which group to keep
        is_partial_match = True
        complete_matched_groups = [entry for entry in matched_ring_entries
                                   if not is_ring_partial_matched(ring, entry.item)]

        if complete_matched_groups:
            is_partial_match = False
            matched_ring_entries = complete_matched_groups

        depth_list = [len(ring_database.ancestors(entry)) for entry in matched_ring_entries]
        most_specific_match_indices = [i for i, x in enumerate(depth_list) if x == max(depth_list)]

        most_specific_matched_entries = [matched_ring_entries[idx] for idx in most_specific_match_indices]
        if len(set(most_specific_matched_entries)) != 1:
            logging.debug('More than one type of node was found to be most specific for this ring.')
            logging.debug('This is either due to a database error in the ring or polycyclic groups, '
                          'or a partial match between the group and the full ring.')
            logging.debug(most_specific_matched_entries)

        # Condense the number of most specific groups down to one
        most_specific_matched_entry = matched_ring_entries[most_specific_match_indices[0]]

        node = most_specific_matched_entry

        if node is None:
            raise DatabaseError('Unable to determine thermo parameters for {0}: no data for {1} or '
                                'any of its ancestors.'.format(molecule, most_specific_match_indices[0]))

        while node is not None and node.data is None:
            # do average of its children
            success, averaged_solute_data = self._average_children_solute(node, ring_database)
            if success:
                node.data = averaged_solute_data
            else:
                node = node.parent

        data = node.data
        comment = node.label
        while isinstance(data, str) and data is not None:
            for entry in ring_database.entries.values():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    node = entry
                    break
        data.comment = '{0}({1})'.format(ring_database.label, comment)

        if solute_data is None:
            return data, node, is_partial_match
        else:
            return add_solute_data(solute_data, data, group_additivity=True, verbose=True), node, is_partial_match
            # By setting verbose=True, we turn on the comments of ring correction to pass the unittest.
            # Typically this comment is very short and also very helpful to check if the ring correction is calculated correctly.

    def _average_children_solute(self, node, database):
        """
        Use children's solute data to guess solute data of parent `node`
        that doesn't have solute data built-in in tree yet.
        For `node` has children that have solute data, return success flag
        `True` and the average solute data.
        For `node` whose children that all have no solute data, return flag
        `False` and None for the solute data.
        """
        if not node.children:
            if node.data is None:
                return False, None
            else:
                return True, node.data
        else:
            children_solute_data_list = []
            for child in node.children:
                if child.data is None:
                    success, child_solute_data_average = self._average_children_solute(child, database)
                    if success:
                        children_solute_data_list.append(child_solute_data_average)
                else:
                    data = child.data
                    while isinstance(data, str):
                        data = database.entries[data].data
                    children_solute_data_list.append(data)
            if children_solute_data_list:
                return True, average_solute_data(children_solute_data_list)
            else:
                return False, None

    def _add_group_solute_data(self, solute_data, database, molecule, atom):
        """
        Determine the group additivity solute data for the atom `atom`
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
        data.comment = '{0}({1})'.format(database.label, comment)

        # This code prints the hierarchy of the found node; useful for debugging
        # result = ''
        # while node is not None:
        #   result = ' -> ' + node + result
        #   node = database.tree.parent[node]
        # print result[4:]
        if solute_data is None:
            return data
        else:
            return add_solute_data(solute_data, data, group_additivity=True)

    def _remove_group_solute_data(self, solute_data, database, molecule, atom):
        """
        Based on the _add_group_solute_data method. Just replace the last line with 'return remove_solute_data()'.
        Determine the group additivity solute data for the atom `atom` in the structure `structure`,
        and REMOVE it from the existing solute data `solute_data`.
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
        data.comment = '{0}({1})'.format(database.label, comment)

        if solute_data is None:
            return data
        else:
            return remove_solute_data(solute_data, data, True)

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

    def get_Kfactor(self, delG298, delH298, delS298, solvent_name, T):
        """Returns a K-factor for the input temperature given the solvation properties of a solute in a solvent
        at 298 K.

        Args:
            delG298 (float): solvation free energy at 298 K in J/mol.
            delH298 (float): solvation enthalpy at 298 K in J/mol.
            delS298 (float): solvation entropy at 298 K in J/mol/K.
            solvent_name (str): name of the solvent that is used in CoolProp.
            T (float): input temperature in K.

        Returns:
            Kfactor (float): K-factor, which is a ratio of the mole fraction of a solute in a gas-phase to
            the mole fraction of a solute in a liquid-phase at equilibrium.

        Raises:
            InputError: if the input temperature is above the critical temperature of the solvent.
            DatabaseError: if the given solvent_name is not available in CoolProp.

        """

        if solvent_name is not None:
            Tc = get_critical_temperature(solvent_name)
            if T < Tc:
                kfactor_parameters = self.get_Kfactor_parameters(delG298, delH298, delS298, solvent_name)
                A = kfactor_parameters.lower_T[0]
                B = kfactor_parameters.lower_T[1]
                C = kfactor_parameters.lower_T[2]
                D = kfactor_parameters.higher_T
                T_transition = kfactor_parameters.T_transition
                rho_c = PropsSI('rhomolar_critical', solvent_name) # critical density of the solvent in mol/m^3
                rho_l = get_liquid_saturation_density(solvent_name, T)  # saturated liquid phase density of the solvent, in mol/m^3
                if T < T_transition:
                    Kfactor = math.exp((A + B * (1 - T / Tc) ** 0.355 + C * math.exp(1 - T / Tc) * (T / Tc) ** 0.59) / (T / Tc))
                else:
                    Kfactor = math.exp(D * (rho_l / rho_c -1) / (T / Tc))
            else:
                raise InputError("The input temperature {0} K cannot be greater than "
                                 "or equal to the critical temperature, {1} K".format(T, Tc))
        else:
            raise DatabaseError("K-factor calculation or temperature-dependent solvation free energy calculation "
                                f"is not available for the given solvent name: {solvent_name}")
        return Kfactor

    def get_T_dep_solvation_energy_from_LSER_298(self, solute_data, solvent_data, T):
        """Returns solvation free energy and K-factor for the input temperature based on the 298 K solvation properties
          calculated from the LSER method.

        Args:
            solute_data: :class:`SoluteData` object `solute_data`.
            solvent_data: :class:`SolventData` object `solvent_data`.
            T (float): input temperature in K.

        Returns:
            delG (float): solvation free energy at the input temperature in J/mol.
            Kfactor (float): K-factor at the input temperature. K-factor is defined as a ratio of the mole fraction
            of a solute in a gas-phase to the mole fraction of a solute in a liquid-phase at equilibrium.

        Raises:
            DatabaseError: if `solute_data.name_in_coolprop` is None or `solute_data` has any missing Abarham or
            Mintz solvent parameters.

        """

        solvent_name = solvent_data.name_in_coolprop
        # check whether all solvent parameters exist
        solvent_param_list = [solvent_data.s_g, solvent_data.b_g, solvent_data.e_g,
                              solvent_data.l_g, solvent_data.a_g, solvent_data.c_g,
                              solvent_data.s_h, solvent_data.b_h, solvent_data.e_h,
                              solvent_data.l_h, solvent_data.a_h, solvent_data.c_h,
                              ]
        param_found_list = [param is not None for param in solvent_param_list]

        if solvent_name is not None and all(param_found_list):
            correction = self.get_solvation_correction(solute_data, solvent_data)
            delG298 = correction.gibbs  # in J/mol
            delH298 = correction.enthalpy  # in J/mol
            delS298 = correction.entropy  # in J/mol/K
        else:
            raise DatabaseError("K-factor parameter calculation is not available for the solvent "
                                "whose `name_in_coolprop` is None or that is missing the Abraham and Mintz "
                                "solvent parameters.")

        return self.get_T_dep_solvation_energy_from_input_298(delG298, delH298, delS298, solvent_name, T)

    def get_T_dep_solvation_energy_from_input_298(self, delG298, delH298, delS298, solvent_name, T):
        """Returns solvation free energy and K-factor for the input temperature based on the given solvation properties
          values at 298 K.

        Args:
            delG298 (float): solvation free energy at 298 K in J/mol.
            delH298 (float): solvation enthalpy at 298 K in J/mol.
            delS298 (float): solvation entropy at 298 K in J/mol/K.
            solvent_name (str): name of the solvent that is used in CoolProp.
            T (float): input temperature in K.

        Returns:
            delG (float): solvation free energy at the input temperature in J/mol.
            Kfactor (float): K-factor at the input temperature. K-factor is defined as a ratio of the mole fraction
            of a solute in a gas-phase to the mole fraction of a solute in a liquid-phase at equilibrium

        """

        Kfactor = self.get_Kfactor(delG298, delH298, delS298, solvent_name, T)
        rho_g = get_gas_saturation_density(solvent_name, T)
        rho_l = get_liquid_saturation_density(solvent_name, T)
        delG = constants.R * T * math.log(Kfactor * rho_g / (rho_l))  # in J/mol
        return delG, Kfactor

    def get_Kfactor_parameters(self, delG298, delH298, delS298, solvent_name, T_trans_factor=0.75):
        """Returns fitted parameters for the K-factor piecewise functions given the solvation properties at 298 K.

        Given delG298 (J/mol), delH298 (J/mol), delS298 (J/mol/K), and solvent_name,
        it finds the fitted K-factor parameters for the solvent-solute pair.
        The parameters (A, B, C, D) are determined by enforcing the smooth continuity of the piecewise
        functions at transition temperature and by making K-factor match in value and temperature
        gradient at 298 K with those estimated from Abraham and Mintz LSERs.
        The four equations are listed here:
            1) Match in value of K-factor with the estimate from the Abraham LSER at 298 K
                A + B(1-Tr)^0.355 + Cexp(1-Tr)(Tr)^0.59 = Tr*ln(K-factor)
            2) Match in temperature gradient of K-factor with the estimate from the Mintz LSER at 298 K
                -0.355B / Tc (1-Tr)^-0.645 + Cexp(1-Tr) / Tc * (0.59(Tr)^-0.41 - (Tr)^0.59) = d(Tr*ln(K-factor))/dT
            3) Continuity of the piecewise function at T_transition:
                A + B(1-Tr)^0.355 + Cexp(1-Tr)(Tr)^0.59 = D(rho_l / rho_c - 1)
            4) Continuous temperature gradient of the piecewise function at T_transition
                -0.355B / Tc (1-Tr)^-0.645 + Cexp(1-Tr) / Tc * (0.59(Tr)^-0.41 - (Tr)^0.59) = D / rho_c * d(rho_l)/dT
        The conversion between dGsolv estimate from the Abraham and K-factor is shown below:
                dGsolv = RTln(K-factor * rho_g / rho_l)
        where rho_g is the saturated gas phase density of the solvent.
        See the RMG documentation "Liquid Phase Systems" section on a temperature-dependent model for more details.

        Args:
            delG298 (float): solvation free energy at 298 K in J/mol.
            delH298 (float): solvation enthalpy at 298 K in J/mol.
            delS298 (float): solvation entropy at 298 K in J/mol/K.
            solvent_name (str): name of the solvent that is used in CoolProp.
            T_trans_factor (float): a temperature [K] for transitioning from the first piecewise function to the
            second piecewise function.

        """

        Tc = get_critical_temperature(solvent_name)
        T_transition = Tc * T_trans_factor  # T_trans_factor is empirically set to 0.75 by default
        rho_c = PropsSI('rhomolar_critical', solvent_name)  # the critical density of the solvent, in mol/m^3

        # Generate Amatrix and bvector for Ax = b
        Amatrix = np.zeros((4, 4))
        bvec = np.zeros((4, 1))
        # 1. Tr*ln(K-factor) value at T = 298 K
        rho_g_298 = get_gas_saturation_density(solvent_name, 298)
        rho_l_298 = get_liquid_saturation_density(solvent_name, 298)
        K298 = math.exp(delG298 / (298 * constants.R)) / rho_g_298 * rho_l_298  # K-factor
        x298 = 298. / Tc * math.log(K298)  # Tr*ln(K-factor), in K
        Amatrix[0][0] = 1
        Amatrix[0][1] = (1 - 298 / Tc) ** 0.355
        Amatrix[0][2] = math.exp(1 - 298 / Tc) * (298 / Tc) ** 0.59
        Amatrix[0][3] = 0
        bvec[0] = x298
        # 2. d(Tr*ln(K-factor)) / dT at T = 298. Use finite difference method to get the temperature gradient from
        # delG, delH, and delS at 298 K
        T2 = 299
        delG_T2 = delH298 - delS298 * T2
        rho_g_T2 = get_gas_saturation_density(solvent_name, T2)
        rho_l_T2 = get_liquid_saturation_density(solvent_name, T2)
        K_T2 = math.exp(delG_T2 / (T2 * constants.R)) / rho_g_T2 * rho_l_T2
        x_T2 = T2 / Tc * math.log(K_T2)  # Tln(K-factor) at 299 K, in K
        slope298 = (x_T2 - x298) / (T2 - 298)
        Amatrix[1][0] = 0
        Amatrix[1][1] = -0.355 / Tc * ((1 - 298 / Tc) ** (-0.645))
        Amatrix[1][2] = 1 / Tc * math.exp(1 - 298 / Tc) * (0.59 * (298 / Tc) ** (-0.41) - (298 / Tc) ** 0.59)
        Amatrix[1][3] = 0
        bvec[1] = slope298
        # 3. Tln(K-factor) continuity at T = T_transition
        rho_l_Ttran = get_liquid_saturation_density(solvent_name, T_transition)
        Amatrix[2][0] = 1
        Amatrix[2][1] = (1 - T_transition / Tc) ** 0.355
        Amatrix[2][2] = math.exp(1 - T_transition / Tc) * (T_transition / Tc) ** 0.59
        Amatrix[2][3] = -(rho_l_Ttran - rho_c) / rho_c
        bvec[2] = 0
        # 4. d(Tln(K-factor)) / dT smooth transition at T = T_transition
        T3 = T_transition + 1
        rho_l_T3 = get_liquid_saturation_density(solvent_name, T3)
        Amatrix[3][0] = 0
        Amatrix[3][1] = -0.355 / Tc * ((1 - T_transition / Tc) ** (-0.645))
        Amatrix[3][2] = 1 / Tc * math.exp(1 - T_transition / Tc) * (
                    0.59 * (T_transition / Tc) ** (-0.41) - (T_transition / Tc) ** 0.59)
        Amatrix[3][3] = - ((rho_l_T3 - rho_l_Ttran) / rho_c / (T3 - T_transition))
        bvec[3] = 0
        # solve for the parameters
        param, residues, ranks, s = np.linalg.lstsq(Amatrix, bvec, rcond=None)
        # store the results in kfactor_parameters class
        kfactor_parameters = KfactorParameters()
        kfactor_parameters.lower_T = [float(param[0]), float(param[1]), float(param[2])]
        kfactor_parameters.higher_T = float(param[3])
        kfactor_parameters.T_transition = T_transition

        return kfactor_parameters
    
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
