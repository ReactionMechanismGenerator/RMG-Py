"""
This tool converts all thermo and kinetics library files in the 
toCompare folder.
"""
import re
import os
import argparse
import logging
import numpy


import rmgpy
from rmgpy.thermo import NASA, ThermoData, Wilhoit, NASAPolynomial
import rmgpy.constants as constants
from rmgpy.kinetics import Arrhenius, ArrheniusEP, ThirdBody, Lindemann, Troe, \
                           PDepArrhenius, MultiArrhenius, MultiPDepArrhenius, \
                           Chebyshev, KineticsData, PDepKineticsModel
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.data.thermo import ThermoLibrary
from __builtin__ import True


name_path_re = re.compile('\.*\/?(.+?)\/RMG-Py-.*-library.*')
def nameFromPath(library):
    match = name_path_re.match(library)
    if match:
        return match.group(1)
    else:
        return os.path.split(library)[0]


    
def readKineticsLibs(libraries, root=''):
    
    # Define local context to allow for loading of the library
    local_context = {
                'KineticsData': KineticsData,
                'Arrhenius': Arrhenius,
                'ArrheniusEP': ArrheniusEP,
                'MultiArrhenius': MultiArrhenius,
                'MultiPDepArrhenius': MultiPDepArrhenius,
                'PDepArrhenius': PDepArrhenius,
                'Chebyshev': Chebyshev,
                'ThirdBody': ThirdBody,
                'Lindemann': Lindemann,
                'Troe': Troe,
                'R': constants.R,
            }
      
    compareDict = {}
    for fileName in libraries:
        # Load the library
        logging.info("Loading reaction library {0}".format(fileName))
        library = KineticsLibrary(label=nameFromPath(fileName))
        library.ALLOW_UNMARKED_DUPLICATES = True
        try:
            # fileName = fixKineticsLibrary(fileName)
            library.load(fileName, local_context=local_context)
        except Exception, e:
            logging.error("Error reading {0}:".format(fileName), exc_info=True)
            logging.warning("Will continue without that model")
            continue

        library.convertDuplicatesToMulti()
        
        # Iterate through the species in the loaded kinetics library to find matches
        # in the comparison data structure
        reactionsToAdd = []
        for entry in library.entries:
            kinetics = library.entries[entry].data
            chemkinReaction = library.entries[entry].item
            
            #containsCarbon = False
            #for species in chemkinReaction.reactants:
            #    if 'C' in species.molecule[0].getFormula():
            #        containsCarbon = True
            #        break
            #if containsCarbon:
            #    continue  # to next library entry
            
            # If the comparison dictionary is empty, add the reaction and skip the rest of the loop
            if not compareDict:
                compareDict[chemkinReaction] = {library.label:kinetics}
                continue
            
            # Iterate through the reactions currently in the comparison dictionary to find matches
            for reaction in compareDict:
                if chemkinReaction.isIsomorphic(reaction, eitherDirection = False):
                    if library.label in compareDict[reaction]:
                        # If the matched reaction from this library is already in the comparison dictionary, skip it and generate a warning
                        logging.warning("Skipping duplicate of reaction {0}, {1} in library {2}".format(reaction, 
                                                                                                        entry, library.label))
                        compareDict[reaction][library.label].comment = "WARNING: ignored a duplicate reaction! " + compareDict[reaction][library.label].comment
                        break
                    # Add the matched reaction to the comparison dictionary
                    logging.debug('Found matching reaction {0} in kinetics library {1}'.format(reaction, library.label))
                    compareDict[reaction][library.label] = kinetics
                    break
            else:
                # If the chemkin reaction did not match any of the reactions in the dict, we will add it later
                logging.debug('Found NEW reaction {0} in kinetics library {1}'.format(reaction, library.label))
                reactionsToAdd.append((chemkinReaction, kinetics))
        
        # Add the reactions which did not match any of the reactions in the comparison file
        for reaction, rate in reactionsToAdd:
            for reaction2 in compareDict:
                if reaction.isIsomorphic(reaction2, eitherDirection = False):
                    logging.warning(("Skipping duplicate of reaction {0} in library {1}, "
                                    "likely a mismatch or different energy state").format(reaction,
                                                                                         library.label))
                    break
            else:
                compareDict[reaction] = {library.label:rate}
        
    return compareDict



def findLibraryFiles(path):
    thermoLibs = []
    kineticsLibs = []
    for root, dirs, files in os.walk(path):
        for name in files:
            path = os.path.join(root, name)
            if root.endswith('RMG-Py-thermo-library') and name == 'ThermoLibrary.py':
                logging.info("Found thermo library {0}".format(path))
                thermoLibs.append(path)
            elif root.endswith('RMG-Py-kinetics-library') and name == 'reactions.py':
                logging.info("Found kinetics file {0}".format(path))
                kineticsLibs.append(path)
            else:
                logging.debug('{0} unread because it is not named like a kinetics or thermo '
                              'library generated by the chemkin importer'.format(path))
    return thermoLibs, kineticsLibs

def main(args):

    thermoLibs = []
    kineticsLibs = []
    for path in args.paths:
        t, k = findLibraryFiles(path)
        thermoLibs.extend(t)
        kineticsLibs.extend(k)

    # This will force the "first" file to be first
    if args.names:
        search = os.path.normpath(args.names)
        for libraryList in (thermoLibs, kineticsLibs):
            for path in libraryList:
                if os.path.normpath(path.split('RMG-Py-')[0]) == search:
                    break  # now 'path' is the one we want first
            else:  # didn't break
                logging.error("Couldn't find {0} in the list {1} to use for names".format(search, libraryList))
                continue
            libraryList.remove(path)
            libraryList.insert(0, path)
            logging.warning('Putting {0} first so its names get used '.format(path))


    kineticsDict = readKineticsLibs(kineticsLibs)

    
    print 'Finished'


def parseCommandLineArguments():
    parser = argparse.ArgumentParser()
    # Options for controlling the amount of information printed to the console
    # By default a moderate level of information is printed; you can either
    # ask for less (quiet), more (verbose), or much more (debug)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-q', '--quiet', action='store_true', help='only print warnings and errors')
    group.add_argument('-v', '--verbose', action='store_true', help='print more verbose output')
    group.add_argument('-d', '--debug', action='store_true', help='print debug information')
    # The important argument(s)
    parser.add_argument('paths', metavar='path', type=str, nargs='+', help='the path(s) to search for mechanisms')
    parser.add_argument('--names', metavar='path', type=str, help='the path to the model which will be used for naming species')
    parser.add_argument('--limit', action='store_true', help='limit to only species in the --names file')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    # Parse the command-line arguments (requires the argparse module)
    args = parseCommandLineArguments()

    # Initialize the logging system
    level = logging.INFO
    if args.debug: level = 0
    elif args.verbose: level = logging.DEBUG
    elif args.quiet: level = logging.WARNING
    logging.basicConfig(level=level)

    main(args)

