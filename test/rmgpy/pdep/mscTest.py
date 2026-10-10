#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
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
Tests for the diagnostics around the modified strong collision solve.
"""

import numpy as np
import pytest

from rmgpy.exceptions import ModifiedStrongCollisionError
from rmgpy.pdep.msc import apply_modified_strong_collision_method


class _EnergyTransferModel:
    def get_alpha(self, T):
        return 1000.0


class _Species:
    def __init__(self):
        self.energy_transfer_model = _EnergyTransferModel()


class _Isomer:
    def __init__(self):
        self.species = [_Species()]


class _StubNetwork:
    """
    The smallest network that drives ``a_mat`` singular.

    Two isomers whose isomerization rate coefficients dwarf the collision frequency make the
    collisional terms on the diagonal negligible, so the matrix reduces to
    ``[[-k21, k12], [k21, -k12]]`` -- determinant zero, exactly rank-1. That is the shape the
    solve fails on in practice, when k(E) has been inflated by a network whose wells and
    transition states are described with different degrees of freedom.
    """

    label = "stub"

    def __init__(self, k_isomerization=1.0e27, coll_freq=1.0e07):
        n_isom, n_grains, n_j = 2, 4, 1
        self.T = 1000.0
        self.P = 1.0e05
        self.e_list = np.linspace(0.0, 3.0e04, n_grains)
        self.j_list = np.array([0], dtype=np.int_)
        self.dens_states = np.ones((n_isom, n_grains, n_j))
        self.coll_freq = np.array([coll_freq, coll_freq])
        self.Kij = np.zeros((n_isom, n_isom, n_grains, n_j))
        self.Kij[1, 0, :, :] = k_isomerization
        self.Kij[0, 1, :, :] = k_isomerization
        self.Fim = np.zeros((n_isom, 0, n_grains, n_j))
        self.Gnj = np.zeros((0, n_isom, n_grains, n_j))
        self.E0 = np.zeros(n_isom)
        self.n_isom = n_isom
        self.n_reac = 0
        self.n_prod = 0
        self.n_grains = n_grains
        self.n_j = n_j
        self.isomers = [_Isomer() for _ in range(n_isom)]
        self.reactants = []
        self.products = []


class TestModifiedStrongCollisionDiagnostics:
    def test_singular_matrix_reports_where_and_why(self):
        """
        A singular solve must say which grain, J index and temperature it failed at, and how far
        k(E) had outrun the collision frequency. Bare ``numpy.linalg.LinAlgError: Singular
        matrix`` carries none of that, and locating it otherwise means bisecting the grains by
        hand.
        """
        with pytest.raises(ModifiedStrongCollisionError) as exc_info:
            apply_modified_strong_collision_method(_StubNetwork(), efficiency_model="none")

        message = str(exc_info.value)
        assert "Singular matrix" in message
        assert "grain 0" in message
        assert "J index 0" in message
        assert "T = 1000 K" in message
        # The two scales, and the ratio between them -- the number that says the input was
        # unphysical rather than the algorithm unlucky.
        assert "1e+27" in message
        assert "1e+07" in message
        assert "1e+20" in message

    def test_diagnostic_survives_a_network_with_no_reactant_channels(self):
        """
        ``Gnj`` is zero-sized when a network has no reactant or product channels, and calling
        ``.max()`` on an empty array raises. The diagnostic must not trade one confusing error
        for another.
        """
        network = _StubNetwork()
        assert network.Gnj.size == 0

        with pytest.raises(ModifiedStrongCollisionError):
            apply_modified_strong_collision_method(network, efficiency_model="none")

    def test_well_conditioned_network_is_untouched(self):
        """
        The diagnostic must only fire on the path that was already failing: a network whose rate
        coefficients sit below the collision frequency solves as before.
        """
        network = _StubNetwork(k_isomerization=1.0e03, coll_freq=1.0e07)
        apply_modified_strong_collision_method(network, efficiency_model="none")
