#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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
Tests for rmgpy.tools.compare_cantera_yaml.
"""

import pytest

from rmgpy.tools.compare_cantera_yaml import _normalize_site_densities, compare_values


class TestNormalizeSiteDensities:
    def test_bare_number_uses_file_units(self):
        """A bare site density is read in the file's default units, as ck2yaml writes it."""
        metadata = {'units': {'length': 'cm', 'quantity': 'mol'},
                    'phases': [{'name': 'gas'}, {'name': 'SURF0', 'site-density': 2.483e-9}]}
        _normalize_site_densities(metadata)
        assert metadata['phases'][1]['site-density'] == pytest.approx(2.483e-8)
        assert 'site-density' not in metadata['phases'][0]

    def test_explicit_units_match_bare_number(self):
        """A site density with explicit units compares equal to the same value given as a bare number."""
        ck = {'units': {'length': 'cm', 'quantity': 'mol'},
              'phases': [{'name': 'SURF0', 'site-density': 2.483e-9}]}
        rmg = {'units': {'length': 'm', 'quantity': 'kmol'},
               'phases': [{'name': 'SURF0', 'site-density': '2.483000e-08 kmol/m^2'}]}
        for metadata in (ck, rmg):
            _normalize_site_densities(metadata)
        assert compare_values(ck['phases'], rmg['phases'], 'metadata.phases') == []

    def test_wrong_site_density_is_reported(self):
        """A site density off by a factor of ten is reported as a numerical difference."""
        ck = {'phases': [{'name': 'SURF0', 'site-density': '2.483e-9 mol/cm^2'}]}
        rmg = {'phases': [{'name': 'SURF0', 'site-density': '2.483e-9 kmol/m^2'}]}
        for metadata in (ck, rmg):
            _normalize_site_densities(metadata)
        differences = compare_values(ck, rmg, 'metadata')
        assert len(differences) == 1
        assert 'Numerical difference at metadata.phases[0].site-density' in differences[0]
