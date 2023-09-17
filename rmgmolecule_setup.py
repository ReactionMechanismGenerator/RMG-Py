import sys
import os
from collections import OrderedDict

try:
    from distutils.core import setup
    from distutils.extension import Extension
except ImportError:
    print("The distutils package is required to build or install RMG Py.")
    raise

try:
    from Cython.Build import cythonize
    from Cython.Compiler import Options
except ImportError:
    print("Cython (http://www.cython.org/) is required to build or install RMG Py.")
    raise

try:
    import numpy
except ImportError:
    print("NumPy (http://numpy.scipy.org/) is required to build or install RMG Py.")
    raise

# Create annotated HTML files for each of the Cython modules
Options.annotate = True

directives = {
    # Set input language version to python 3
    "language_level": 3,
    # Turn on profiling capacity for all Cython modules
    # 'profile': True,
    # Embed call signatures in cythonized files - enable when building documentation
    # 'embedsignature': True,
}


main_ext_modules = [
    # Molecules and molecular representations
    Extension(
        "rmgpy.molecule.atomtype",
        ["rmgpy/molecule/atomtype.py"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.molecule.element",
        ["rmgpy/molecule/element.py"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.molecule.graph",
        ["rmgpy/molecule/graph.pyx"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.molecule.group",
        ["rmgpy/molecule/group.py"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.molecule.molecule",
        ["rmgpy/molecule/molecule.py"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.molecule.symmetry",
        ["rmgpy/molecule/symmetry.py"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.molecule.vf2",
        ["rmgpy/molecule/vf2.pyx"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.molecule.converter",
        ["rmgpy/molecule/converter.py"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.molecule.translator",
        ["rmgpy/molecule/translator.py"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.molecule.util",
        ["rmgpy/molecule/util.py"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.molecule.inchi",
        ["rmgpy/molecule/inchi.py"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.molecule.resonance",
        ["rmgpy/molecule/resonance.py"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.molecule.pathfinder",
        ["rmgpy/molecule/pathfinder.py"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.molecule.kekulize",
        ["rmgpy/molecule/kekulize.pyx"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.constants",
        ["rmgpy/constants.py"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.quantity",
        ["rmgpy/quantity.py"],
        include_dirs=["."],
    ),
    Extension(
        "rmgpy.rmgobject",
        ["rmgpy/rmgobject.pyx"],
    ),
]

ext_modules = []
if "install" in sys.argv:
    # This is so users can still do simply `python setup.py install`
    ext_modules.extend(main_ext_modules)
if "main" in sys.argv:
    # This is for `python setup.py build_ext main`
    sys.argv.remove("main")
    ext_modules.extend(main_ext_modules)
if "minimal" in sys.argv:
    # This starts with the full install list, but removes anything that has a pure python mode
    # i.e. in only includes things whose source is .pyx
    sys.argv.remove("minimal")
    temporary_list = []
    temporary_list.extend(main_ext_modules)
    for module in temporary_list:
        for source in module.sources:
            if os.path.splitext(source)[1] == ".pyx":
                ext_modules.append(module)

# Remove duplicates while preserving order:
ext_modules = list(OrderedDict.fromkeys(ext_modules))

scripts = []

modules = ["rmgpy.exceptions", "rmgpy.version"]

__version__ = "1.0.0"

import logging

logging.error(ext_modules)
# Initiate the build and/or installation
setup(
    name="RMG-Py",
    version=__version__,
    description="Reaction Mechanism Generator",
    author="William H. Green and the RMG Team",
    author_email="rmg_dev@mit.edu",
    url="http://reactionmechanismgenerator.github.io",
    packages=[
        "rmgpy.molecule",
    ],
    py_modules=modules,
    scripts=scripts,
    ext_modules=cythonize(
        ext_modules,
        build_dir="build",
        compiler_directives=directives,
    ),
    include_dirs=[".", numpy.get_include()],
)
