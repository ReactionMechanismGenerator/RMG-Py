#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script checks for the major dependencies that RMG-Py requires.
"""

import re
import subprocess

missing = False

print '\nChecking vital dependencies...\n'
print '{0:<15}{1:<15}{2}'.format('Package', 'Version', 'Location')

# Check for symmetry
try:
    result = subprocess.check_output('symmetry -h', stderr=subprocess.STDOUT, shell=True)
except subprocess.CalledProcessError:
    print '{0:<30}{1}'.format('symmetry', 'Not found. Please install in order to use QM.')
    missing = True
else:
    match = re.search(r'\$Revision: (\S*) \$', result)
    version = match.group(1)
    location = subprocess.check_output('which symmetry', shell=True)

    print '{0:<15}{1:<15}{2}'.format('symmetry', version, location.strip())

# Check for RDKit
try:
    import rdkit
    from rdkit import Chem
except ImportError:
    print '{0:<30}{1}'.format('RDKit', 'Not found. Please install RDKit version 2015.03.1 or later with InChI support.')
    missing = True
else:
    try:
        version = rdkit.__version__
    except AttributeError:
        version = False
    location = rdkit.__file__
    inchi = Chem.inchi.INCHI_AVAILABLE

    if version:
        print '{0:<15}{1:<15}{2}'.format('RDKit', version, location)
        if not inchi:
            print 'RDKit installed without InChI Support. Please install with InChI.'
            missing = True
    else:
        print 'RDKit version out of date, please install RDKit version 2015.03.1 or later with InChI support.'
        missing = True

# Check for OpenBabel
try:
    import openbabel
except ImportError:
    print '{0:<30}{1}'.format('OpenBabel', 'Not found. Necessary for SMILES/InChI functionality for nitrogen compounds.')
    missing = True
else:
    location = openbabel.__file__
    print '{0:<30}{1}'.format('OpenBabel', location)

# Check for lpsolve
try:
    import lpsolve55
except ImportError:
    print '{0:<30}{1}'.format('lpsolve55', 'Not found. Necessary for generating Clar structures for aromatic species.')
    missing = True
else:
    location = lpsolve55.__file__
    print '{0:<30}{1}'.format('lpsolve55', location)

# Check for pyrdl
try:
    import py_rdl
except ImportError:
    print '{0:<30}{1}'.format('pyrdl', 'Not found. Necessary for ring perception algorithms.')
    missing = True
else:
    location = py_rdl.__file__
    print '{0:<30}{1}'.format('pyrdl', location)

if missing:
    print """
There are missing dependencies as listed above. Please install them before proceeding.

Using Anaconda, these dependencies can be installed from the RMG channel as follows:

    conda install -c rmg <package name>

Be sure to activate your conda environment (rmg_env by default) before installing.
"""
else:
    print """
Everything was found :)
"""
