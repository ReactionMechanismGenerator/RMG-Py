#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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

"""
Contains classes that define an internal representation of a unimolecular
reaction network.
"""

import math
import numpy
import cython

import chempy.constants as constants
import chempy.states as states

################################################################################

class Network:
    """
    A representation of a unimolecular reaction network. The attributes are:

    =================== ======================= ================================
    Attribute           Type                    Description
    =================== ======================= ================================
    `isomers`           ``list``                A list of the unimolecular isomers in the network
    `reactants`         ``list``                A list of the bimolecular reactant channels in the network
    `products`          ``list``                A list of the bimolecular product channels in the network
    `pathReactions`     ``list``                A list of reaction objects that connect adjacent isomers (the high-pressure-limit)
    `bathGas`           :class:`Species`        The bath gas
    `collisionModel`    :class:`CollisionModel` The collision model to use
    `netReactions`      ``list``                A list of reaction objects that connect any pair of isomers
    =================== ======================= ================================

    """

    def __init__(self, isomers=None, reactants=None, products=None, pathReactions=None, bathGas=None):
        self.isomers = isomers or []
        self.reactants = reactants or []
        self.products = products or []
        self.pathReactions = pathReactions or []
        self.bathGas = bathGas
        self.netReactions = []

################################################################################
