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
import math
import argparse

from rmgpy.chem.kinetics import KineticsData
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.kinetics import KineticsDatabase, KineticsGroups

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

def generateKineticsGroups(inputDirectory, family, Tdata):
    """
    Generate kinetics group additivity values for the given reaction `family`
    name at the specified temperatures `Tdata`. The `inputDirectory` is the
    path to the old-style RMG database to load the rate rules from.
    """

    database = KineticsDatabase()

    root = os.path.join(inputDirectory, 'kinetics_groups', family)
    if os.path.exists(os.path.join(root, 'dictionary.txt')) and os.path.exists(os.path.join(root, 'rateLibrary.txt')):

        # Load the old kinetics rate rules
        group = KineticsGroups(label=os.path.basename(root), name=os.path.basename(root))
        group.loadOld(root)
        database.groups[group.label] = group

        # Generate the new kinetics group additivity values
        groupValues0, groupUncertainties0, groupCounts, templates, kdata, kmodel = group.fitGroupValuesFromOldLibrary(root, Tdata)

        # Postprocess the fitted group values and uncertainies
        groupValues = {}
        for entry, data in groupValues0.iteritems():
            if data is not None:
                groupValues[entry.label] = data
        groupUncertainties = {}
        for entry, data in groupUncertainties0.iteritems():
            if data is not None:
                if not any([(math.isinf(d) or math.isnan(d)) for d in data]):
                    groupUncertainties[entry.label] = data

    return groupValues, groupUncertainties

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
    print 'Generating kinetics group additivity values...'
    Tdata = [300,400,500,600,800,1000,1500,2000]
    for familyLabel, family in database.kinetics.groups.iteritems():
        print '    %s...' % (familyLabel)

        # Determine units of generated rate coefficients
        numReactants = len(family.forwardTemplate.reactants)
        if numReactants == 1:
            kunits = 's^-1'
        elif numReactants == 2:
            kunits = 'm^3/(mol*s)'
        else:
            raise Exception("Unexpected number of reactants %i in forward template of reaction family %s." % (numReactants, familyLabel))

        # Generate and store the group values
        groupValues, groupUncertainties = generateKineticsGroups(inputDirectory, familyLabel, Tdata)
        for label, entry in family.entries.iteritems():
            if label in groupValues:
                if label in groupUncertainties:
                    kdata = (groupValues[label],kunits,"*|/",groupUncertainties[label])
                else:
                    kdata = (groupValues[label],kunits)
                entry.data = KineticsData(
                    Tdata=(Tdata,"K"),
                    kdata=kdata,
                    Tmin=(Tdata[0],"K"),
                    Tmax=(Tdata[-1],"K"),
                    comment='Fitted from RMG-Java rate rules',
                )
            else:
                entry.data = None

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
