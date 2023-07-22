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
This module contains unit tests of the rmgpy.element module.
"""


import rmgpy.molecule.element
from rmgpy.molecule.element import Element


class TestElement:
    """
    Contains unit tests of the Element class.
    """

    def setup_class(self):
        """
        A function run before each unit test in this class.
        """
        self.element = rmgpy.molecule.element.C
        self.element_x = rmgpy.molecule.element.X

    def test_pickle(self):
        """
        Test that an Element object can be successfully pickled and
        unpickled with no loss of information.
        """
        import pickle

        element = pickle.loads(pickle.dumps(self.element))
        assert self.element.number == element.number
        assert self.element.symbol == element.symbol
        assert self.element.name == element.name
        assert self.element.mass == element.mass

    def test_output(self):
        """
        Test that we can reconstruct an Element object from its repr()
        output with no loss of information.
        """
        namespace = {}
        exec("element = {0!r}".format(self.element), globals(), namespace)
        assert "element" in namespace
        element = namespace["element"]
        assert self.element.number == element.number
        assert self.element.symbol == element.symbol
        assert self.element.name == element.name
        assert self.element.mass == element.mass

    def test_get_element(self):
        """
        Test the rmgpy.elements.get_element() method.
        """
        assert rmgpy.molecule.element.get_element(6) is self.element
        assert rmgpy.molecule.element.get_element("C") is self.element
        assert rmgpy.molecule.element.get_element(0) is self.element_x
        assert rmgpy.molecule.element.get_element("X") is self.element_x

    def test_get_element_isotope(self):
        """
        Test that the rmgpy.elements.get_element() method works for isotopes.
        """
        assert isinstance(rmgpy.molecule.element.get_element("C", isotope=13), Element)
        assert isinstance(rmgpy.molecule.element.get_element(6, isotope=13), Element)

    def test_chemkin_name(self):
        """
        Test that retrieving the chemkin name of an element works.
        """
        d = rmgpy.molecule.element.get_element("H", isotope=2)
        assert d.chemkin_name == "D"

        c13 = rmgpy.molecule.element.get_element("C", isotope=13)
        assert c13.chemkin_name == "CI"

        o18 = rmgpy.molecule.element.get_element("O", isotope=18)
        assert o18.chemkin_name == "OI"
