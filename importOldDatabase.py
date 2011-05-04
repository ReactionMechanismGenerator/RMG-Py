#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script converts the old RMG database into the new format. This includes
generating the kinetics group additivity values from the old rate rules. To
use, pass the root directories of the old database and the directory to save
the new database.
"""

import os
import time
import argparse

from rmgpy.data.rmg import RMGDatabase

################################################################################

def parseCommandLineArguments():
    """
    Parse the command-line arguments being passed to RMG Py. This uses the
    :mod:`argparse` module, which ensures that the command-line arguments are
    sensible, parses them, and returns them.
    """

    parser = argparse.ArgumentParser(description=
"""
Convert the old RMG-Java database into the new format, including generating
kinetics group additivity values from the rate rules.
""")
    parser.add_argument('-i', '--input', metavar='OLDDATABASE', type=str, nargs=1, required=True,
        help='The path to the old RMG database')
    parser.add_argument('-o', '--output', metavar='NEWDATABASE', type=str, nargs=1, required=True,
        help='The path to save the new RMG database')

    return parser.parse_args()

################################################################################

if __name__ == '__main__':

    # Parse the command-line arguments
    args = parseCommandLineArguments()
    inputDirectory = os.path.abspath(args.input[0])
    outputDirectory = os.path.abspath(args.output[0])
    
    # Load the old RMG-Java database
    print 'Loading old RMG-Java database...'
    database = RMGDatabase()
    database.loadOld(inputDirectory)

    # Generate kinetics group additivity values from old rate rules
    print 'Warning: Not yet generating kinetics group additivity values!'

    # Add history item to each entry in each database (saying that JWA added them all!)
    print 'Setting history of all entries in database...'
    event = [time.asctime(),'jwallen','action','jwallen imported this entry from the old RMG database.']
    for depository in database.thermo.depository.values():
        for label, entry in depository.entries.iteritems():
            entry.history.append(event)
    for library in database.thermo.libraries.values():
        for label, entry in library.entries.iteritems():
            entry.history.append(event)
    for groups in database.thermo.groups.values():
        for label, entry in groups.entries.iteritems():
            entry.history.append(event)
    for depository in database.kinetics.depository.values():
        for label, entry in depository.entries.iteritems():
            entry.history.append(event)
    for library in database.kinetics.libraries.values():
        for label, entry in library.entries.iteritems():
            entry.history.append(event)
    for groups in database.kinetics.groups.values():
        for label, entry in groups.entries.iteritems():
            entry.history.append(event)

    # Save the new database
    print 'Saving the new RMG-Py database...'
    database.save(outputDirectory)
