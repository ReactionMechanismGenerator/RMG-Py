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

import unittest
import os
import os.path
import shutil
from nose.plugins.attrib import attr
import rmgpy
from rmgpy.tools.fluxdiagram import createFluxDiagram
@attr('functional')
class FluxDiagramTest(unittest.TestCase):

    def test_avi_simple(self):
        folder = os.path.join(os.path.dirname(rmgpy.__file__),'tools','data','flux')
        
        inputFile = os.path.join(folder,'input_simple.py')
        chemkinFile = os.path.join(folder, 'chemkin', 'chem.inp')
        dictFile = os.path.join(folder, 'chemkin', 'species_dictionary.txt')
        settings = {'maximumNodeCount': 50,
                    'maximumEdgeCount': 50,
                    'concentrationTolerance': 1e-6,
                    'speciesRateTolerance': 1e-6,
                    'maximumNodePenWidth': 10.0,
                    'maximumEdgePenWidth': 10.0,
                    'radius': 2,
                    'centralReactionCount': 2,
                    'timeStep': 10**0.1}
        createFluxDiagram(inputFile, chemkinFile, dictFile,
                          centralSpeciesList=[1], superimpose=True, settings=settings)

        outputdir = os.path.join(folder,'flux')
        simfile = os.path.join(outputdir,'1','flux_diagram.avi')

        speciesdir = os.path.join(folder,'species')

        self.assertTrue(os.path.isfile(simfile))
        
        shutil.rmtree(outputdir)
        shutil.rmtree(speciesdir)

    def test_avi_liquid(self):
        folder = os.path.join(os.path.dirname(rmgpy.__file__), 'tools', 'data', 'flux')

        inputFile = os.path.join(folder, 'input_liquid.py')
        chemkinFile = os.path.join(folder, 'chemkin', 'chem.inp')
        dictFile = os.path.join(folder, 'chemkin', 'species_dictionary.txt')
        createFluxDiagram(inputFile, chemkinFile, dictFile, diffusionLimited=False)

        outputdir = os.path.join(folder, 'flux')
        simfile = os.path.join(outputdir, '1', 'flux_diagram.avi')

        speciesdir = os.path.join(folder, 'species')

        self.assertTrue(os.path.isfile(simfile))

        shutil.rmtree(outputdir)
        shutil.rmtree(speciesdir)

    def tearDown(self):
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None
