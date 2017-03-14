#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2015 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
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

import sys
import os

try:
    from distutils.core import setup
    from distutils.extension import Extension
except ImportError:
    print 'The distutils package is required to build or install RMG Py.'
    
try:
    from Cython.Distutils import build_ext
    from Cython.Build import cythonize
    import Cython.Compiler.Options
except ImportError:
    print 'Cython (http://www.cython.org/) is required to build or install RMG Py.'
    
try:
    import numpy
except ImportError:
    print 'NumPy (http://numpy.scipy.org/) is required to build or install RMG Py.'

# Create annotated HTML files for each of the Cython modules
Cython.Compiler.Options.annotate = True

# Turn on profiling capacity for all Cython modules
#Cython.Compiler.Options.directive_defaults['profile'] = True
# turn on or off cython coverage
cython_coverage = 1

################################################################################

def getMainExtensionModules():
    return [
        # Kinetics
        Extension('rmgpy.kinetics.arrhenius', ['rmgpy/kinetics/arrhenius.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.kinetics.chebyshev', ['rmgpy/kinetics/chebyshev.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.kinetics.kineticsdata', ['rmgpy/kinetics/kineticsdata.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.kinetics.falloff', ['rmgpy/kinetics/falloff.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.kinetics.model', ['rmgpy/kinetics/model.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.kinetics.tunneling', ['rmgpy/kinetics/tunneling.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        # Molecules and molecular representations
        Extension('rmgpy.molecule.atomtype', ['rmgpy/molecule/atomtype.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.molecule.element', ['rmgpy/molecule/element.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.molecule.graph', ['rmgpy/molecule/graph.pyx'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.molecule.group', ['rmgpy/molecule/group.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.molecule.molecule', ['rmgpy/molecule/molecule.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.molecule.symmetry', ['rmgpy/molecule/symmetry.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.molecule.vf2', ['rmgpy/molecule/vf2.pyx'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.molecule.parser', ['rmgpy/molecule/parser.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.molecule.generator', ['rmgpy/molecule/generator.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.molecule.util', ['rmgpy/molecule/util.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.molecule.inchi', ['rmgpy/molecule/inchi.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.molecule.resonance', ['rmgpy/molecule/resonance.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.molecule.pathfinder', ['rmgpy/molecule/pathfinder.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        # Pressure dependence
        Extension('rmgpy.pdep.collision', ['rmgpy/pdep/collision.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.pdep.configuration', ['rmgpy/pdep/configuration.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.pdep.me', ['rmgpy/pdep/me.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.pdep.msc', ['rmgpy/pdep/msc.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.pdep.reaction', ['rmgpy/pdep/reaction.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.pdep.rs', ['rmgpy/pdep/rs.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.pdep.cse', ['rmgpy/pdep/cse.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        # Statistical mechanics
        Extension('rmgpy.statmech.conformer', ['rmgpy/statmech/conformer.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.statmech.mode', ['rmgpy/statmech/mode.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.statmech.rotation', ['rmgpy/statmech/rotation.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.statmech.schrodinger', ['rmgpy/statmech/schrodinger.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.statmech.torsion', ['rmgpy/statmech/torsion.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.statmech.translation', ['rmgpy/statmech/translation.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.statmech.vibration', ['rmgpy/statmech/vibration.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        # Thermodynamics
        Extension('rmgpy.thermo.thermodata', ['rmgpy/thermo/thermodata.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.thermo.model', ['rmgpy/thermo/model.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.thermo.nasa', ['rmgpy/thermo/nasa.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.thermo.wilhoit', ['rmgpy/thermo/wilhoit.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        # Miscellaneous
        Extension('rmgpy.constants', ['rmgpy/constants.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.quantity', ['rmgpy/quantity.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.reaction', ['rmgpy/reaction.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.species', ['rmgpy/species.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
    ]
    
def getSolverExtensionModules():
    return [
        Extension('rmgpy.solver.base', ['rmgpy/solver/base.pyx'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.solver.simple', ['rmgpy/solver/simple.pyx'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.solver.liquid', ['rmgpy/solver/liquid.pyx'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
    ]

def getCanthermExtensionModules():
    return [
        # Kinetics
        Extension('rmgpy.kinetics.arrhenius', ['rmgpy/kinetics/arrhenius.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.kinetics.chebyshev', ['rmgpy/kinetics/chebyshev.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.kinetics.kineticsdata', ['rmgpy/kinetics/kineticsdata.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.kinetics.falloff', ['rmgpy/kinetics/falloff.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.kinetics.model', ['rmgpy/kinetics/model.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.kinetics.tunneling', ['rmgpy/kinetics/tunneling.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        # Pressure dependence
        Extension('rmgpy.pdep.collision', ['rmgpy/pdep/collision.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.pdep.configuration', ['rmgpy/pdep/configuration.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.pdep.me', ['rmgpy/pdep/me.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.pdep.msc', ['rmgpy/pdep/msc.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.pdep.reaction', ['rmgpy/pdep/reaction.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.pdep.rs', ['rmgpy/pdep/rs.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.pdep.cse', ['rmgpy/pdep/cse.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        # Statistical mechanics
        Extension('rmgpy.statmech.conformer', ['rmgpy/statmech/conformer.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.statmech.mode', ['rmgpy/statmech/mode.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.statmech.rotation', ['rmgpy/statmech/rotation.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.statmech.schrodinger', ['rmgpy/statmech/schrodinger.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.statmech.torsion', ['rmgpy/statmech/torsion.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.statmech.translation', ['rmgpy/statmech/translation.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.statmech.vibration', ['rmgpy/statmech/vibration.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        # Thermodynamics
        Extension('rmgpy.thermo.thermodata', ['rmgpy/thermo/thermodata.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.thermo.model', ['rmgpy/thermo/model.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.thermo.nasa', ['rmgpy/thermo/nasa.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.thermo.wilhoit', ['rmgpy/thermo/wilhoit.pyx'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        # Miscellaneous
        Extension('rmgpy.constants', ['rmgpy/constants.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
        Extension('rmgpy.quantity', ['rmgpy/quantity.py'], include_dirs=['.'], define_macros=[('CYTHON_TRACE',cython_coverage)]),
    ]

################################################################################

ext_modules = []
if 'install' in sys.argv:
    # This is so users can still do simply `python setup.py install`
    ext_modules.extend(getMainExtensionModules())
    ext_modules.extend(getSolverExtensionModules())
if 'main' in sys.argv:
    # This is for `python setup.py build_ext main`
    sys.argv.remove('main')
    ext_modules.extend(getMainExtensionModules())
if 'solver' in sys.argv:
    # This is for `python setup.py build_ext solver`
    sys.argv.remove('solver')
    ext_modules.extend(getSolverExtensionModules())
if 'cantherm' in sys.argv:
    # This is for `python setup.py build_ext cantherm`
    sys.argv.remove('cantherm')
    ext_modules.extend(getMainExtensionModules())
    ext_modules.extend(getCanthermExtensionModules())
if 'minimal' in sys.argv:
    # This starts with the full install list, but removes anything that has a pure python mode
    # i.e. in only includes things whose source is .pyx
    sys.argv.remove('minimal')
    temporary_list = []
    temporary_list.extend(getMainExtensionModules())
    temporary_list.extend(getSolverExtensionModules())
    for module in temporary_list:
        for source in module.sources:
            if os.path.splitext(source)[1] == '.pyx':
                ext_modules.append(module)

scripts=['cantherm.py', 'rmg.py', 'scripts/diffModels.py', 'scripts/generateFluxDiagram.py',
         'scripts/generateReactions.py', 'scripts/mergeModels.py','scripts/sensitivity.py', 'scripts/thermoEstimator.py',
         'testing/databaseTest.py']

modules = []
for root, dirs, files in os.walk('rmgpy'):
    if 'test_data' in root:
        continue
    for file in files:
        if file.endswith('.py') or file.endswith('.pyx'):
            if 'Test' not in file and '__init__' not in file:
                if not root.endswith('rmgpy/cantherm/files'):
                    module = 'rmgpy' + root.partition('rmgpy')[-1].replace('/','.') + '.' + file.partition('.py')[0]
                    modules.append(module)       

# Initiate the build and/or installation

# Read the version number
exec(open('rmgpy/version.py').read())

setup(name='RMG-Py',
    version= __version__,
    description='Reaction Mechanism Generator',
    author='William H. Green and the RMG Team',
    author_email='rmg_dev@mit.edu',
    url='http://reactionmechanismgenerator.github.io',
    packages=['rmgpy'],
    py_modules = modules,
    scripts=scripts,
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize(ext_modules, compiler_directives={'linetrace': bool(cython_coverage)}),
    include_dirs=['.', numpy.get_include()],
)
