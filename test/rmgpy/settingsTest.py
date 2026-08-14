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
This script contains unit tests of the :mod:`rmgpy` settings object.
"""

import os

from rmgpy import Settings


class TestDefaultDatabaseDirectory:
    """
    Tests for how the default database directory is chosen when no rmgrc
    supplies one.
    """

    def test_paired_database_worktree_is_used(self, tmp_path):
        """
        A checkout named RMG-Py-<suffix> uses RMG-database-<suffix> beside it.
        """
        (tmp_path / "RMG-Py-feature").mkdir()
        (tmp_path / "RMG-database-feature" / "input").mkdir(parents=True)

        directory, source = Settings.default_database_directory(str(tmp_path / "RMG-Py-feature"))

        assert directory == os.path.realpath(str(tmp_path / "RMG-database-feature" / "input"))
        assert "paired" in source

    def test_falls_back_when_the_paired_worktree_is_absent(self, tmp_path):
        """
        The suffix alone is not enough: with no RMG-database-<suffix> on disk,
        the historical default is used, so an unpaired worktree is unaffected.
        """
        (tmp_path / "RMG-Py-feature").mkdir()
        (tmp_path / "RMG-database" / "input").mkdir(parents=True)

        directory, source = Settings.default_database_directory(str(tmp_path / "RMG-Py-feature"))

        assert directory == os.path.realpath(str(tmp_path / "RMG-database" / "input"))
        assert source == "Default, relative to RMG-Py source code"

    def test_plain_checkout_name_is_unaffected(self, tmp_path):
        """
        An ordinary RMG-Py checkout keeps the behaviour it has always had,
        even when a similarly named database directory happens to be present.
        """
        (tmp_path / "RMG-Py").mkdir()
        (tmp_path / "RMG-database" / "input").mkdir(parents=True)

        directory, source = Settings.default_database_directory(str(tmp_path / "RMG-Py"))

        assert directory == os.path.realpath(str(tmp_path / "RMG-database" / "input"))
        assert source == "Default, relative to RMG-Py source code"
