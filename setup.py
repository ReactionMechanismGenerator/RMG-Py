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

from setuptools import setup
    
try:
    from Cython.Build import cythonize
    from Cython.Compiler import Options
except ImportError:
    print('Cython (https://cython.org/) is required to build or install RMG Py.')
    raise
    
try:
    import numpy
except ImportError:
    print('NumPy (https://numpy.scipy.org/) is required to build or install RMG Py.')
    raise

from setuptools import find_packages

# Create annotated HTML files for each of the Cython modules
Options.annotate = False

directives = {
    # Set input language version to python 3
    'language_level': 3,
    # Turn on profiling capacity for all Cython modules
    # 'profile': True,
    # Embed call signatures in cythonized files - enable when building documentation
    # 'embedsignature': True,
}

ext_modules = [
    # RMG
    'rmgpy/rmgobject.pyx',
    # Kinetics
    'rmgpy/kinetics/arrhenius.pyx',
    'rmgpy/kinetics/chebyshev.pyx',
    'rmgpy/kinetics/kineticsdata.pyx',
    'rmgpy/kinetics/falloff.pyx',
    'rmgpy/kinetics/model.pyx',
    'rmgpy/kinetics/tunneling.pyx',
    'rmgpy/kinetics/surface.pyx',
    'rmgpy/kinetics/uncertainties.pyx',
    # Molecules and molecular representations
    'rmgpy/molecule/atomtype.py',
    'rmgpy/molecule/element.py',
    'rmgpy/molecule/graph.pyx',
    'rmgpy/molecule/group.py',
    'rmgpy/molecule/molecule.py',
    'rmgpy/molecule/symmetry.py',
    'rmgpy/molecule/vf2.pyx',
    'rmgpy/molecule/converter.py',
    'rmgpy/molecule/translator.py',
    'rmgpy/molecule/util.py',
    'rmgpy/molecule/inchi.py',
    'rmgpy/molecule/resonance.py',
    'rmgpy/molecule/pathfinder.py',
    'rmgpy/molecule/kekulize.pyx',
    # Pressure dependence
    'rmgpy/pdep/collision.pyx',
    'rmgpy/pdep/configuration.pyx',
    'rmgpy/pdep/me.pyx',
    'rmgpy/pdep/msc.pyx',
    'rmgpy/pdep/reaction.pyx',
    'rmgpy/pdep/rs.pyx',
    'rmgpy/pdep/cse.pyx',
    # Statistical mechanics
    'rmgpy/statmech/conformer.pyx',
    'rmgpy/statmech/mode.pyx',
    'rmgpy/statmech/rotation.pyx',
    'rmgpy/statmech/schrodinger.pyx',
    'rmgpy/statmech/torsion.pyx',
    'rmgpy/statmech/translation.pyx',
    'rmgpy/statmech/vibration.pyx',
    # Thermodynamics
    'rmgpy/thermo/thermodata.pyx',
    'rmgpy/thermo/model.pyx',
    'rmgpy/thermo/nasa.pyx',
    'rmgpy/thermo/wilhoit.pyx',
    # Miscellaneous
    'rmgpy/constants.py',
    'rmgpy/quantity.py',
    'rmgpy/reaction.py',
    'rmgpy/species.py',
    'rmgpy/chemkin.pyx',
    # solvers
    'rmgpy/solver/base.pyx',
    'rmgpy/solver/simple.pyx',
    'rmgpy/solver/liquid.pyx',
    'rmgpy/solver/mbSampled.pyx',
    'rmgpy/solver/surface.pyx',
]

scripts = [
    'rmg.py',
    'Arkane.py',
    'scripts/checkModels.py',
    'scripts/diffModels.py',
    'scripts/generateChemkinHTML.py',
    'scripts/generateFluxDiagram.py',
    'scripts/generateReactions.py',
    'scripts/machineWriteDatabase.py',
    'scripts/mergeModels.py',
    'scripts/rmg2to3.py',
    'scripts/simulate.py',
    'scripts/standardizeModelSpeciesNames.py',
    'scripts/thermoEstimator.py',
    'scripts/isotopes.py',
]

# Read the version number
exec(open('rmgpy/version.py').read())

# Initiate the build and/or installation
setup(
    name='reactionmechanismgenerator',
    version=__version__,
    description='Reaction Mechanism Generator',
    author='William H. Green and the RMG Team',
    author_email='rmg_dev@mit.edu',
    url='http://reactionmechanismgenerator.github.io',
    python_requires='>=3.9,<3.12',
    setup_requires=['numpy'],
    packages=find_packages(where='.', include=["rmgpy*"]) + find_packages(where='.', include=["arkane*"]),
    scripts=scripts,
    entry_points={
        'console_scripts': [
            'rmg.py = rmgpy.__main__:main',
            'Arkane.py = arkane.__main__:main',
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.pxd"],
    },
    ext_modules=cythonize(ext_modules, compiler_directives=directives),
    include_dirs=numpy.get_include(),
    install_requires=["Cython", "numpy"],
)
