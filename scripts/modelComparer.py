# This tool converts all thermo and kinetics library files in the 
# toCompare folder to the data structures:
#
# Kinetics:
# {reaction1:{model1:rate, model2:rate, ...}, 
#  reaction2:{model1:rate, model2:rate, ...}, ...}
#
# Thermo:
# {species1:{model1:thermo, model2:thermo, ...}, 
#  species2:{model1:thermo, model2:thermo, ...}, ...}
# 
# These dictionaries are then used to perform statistical analysis on 
# the group of models.
import re
import os
import argparse
import logging
import numpy
from math import log10
from xlwt import Workbook, Formula
from xlwt.Utils import rowcol_to_cell

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


def readThermoLibs(libraries, root=''):

    # Define local context to allow for loading of the library
    local_context = {
            'ThermoData': ThermoData,
            'Wilhoit': Wilhoit,
            'NASAPolynomial': NASAPolynomial,
            'NASA': NASA,
        }
      
    compareDict = {}
    namesDict = {}
    for fileName in libraries:
        # Load the library
        library = ThermoLibrary(label=nameFromPath(fileName))
        library.SKIP_DUPLICATES = True
        library.load(fileName, local_context = local_context)
        
        # Iterate through the species in the loaded thermo library to find matches
        # in the comparison data structure
        speciesToAdd = []
        for entry in library.entries:
            thermo = library.entries[entry].data
            chemkinMolecule = library.entries[entry].item
            ##todo: temporary skip!
            #if 'C' in chemkinMolecule.getFormula():
            #    continue

            name = library.entries[entry].label
            
            # Iterate through the species currently in the comparison dictionary
            # to find matches of the current chemkinMolecule
            for species in compareDict:
                if chemkinMolecule.isIsomorphic(species):
                    logging.debug('Found isomorphic species {0} in thermo library {1}'.format(species.toSMILES(), library.label))
                    if library.label in compareDict[species]:
                        # If the matched species is already in the comparison dictionary, skip it and generate a warning
                        old_name = namesDict[species][library.label]
                        logging.warning('Duplicate "{1}" of species "{old}" {0} in library {2}, likely a mismatch or different energy state'.format(species.toSMILES(), entry, library.label, old=old_name))
                        previous_thermo = compareDict[species][library.label]
                        if previous_thermo.getEnthalpy(300) < thermo.getEnthalpy(300) :
                            logging.warning('Using thermo of lower energy version {old} and ignoring subsequent {new}'.format(new=name, old=old_name))
                            break
                        else:
                            logging.warning('Using thermo of lower energy version {new} and dropping initial {old}'.format(new=name, old=old_name))
                    # Add the matched species to the comparison dictionary
                    compareDict[species][library.label] = thermo
                    namesDict[species][library.label] = name
                    break
            else:
                # If the chemkinMolecule did not match any of the species already in the dict,
                # we will add it later
                logging.debug('Found NEW species {0} in thermo library {1}'.format(chemkinMolecule.toSMILES(), library.label))
                speciesToAdd.append((chemkinMolecule, thermo, name))
        
        # Add the species which did not match any of the species in the comparison file
        for species, thermo, name in speciesToAdd:
            for species2 in compareDict:
                if species.isIsomorphic(species2):
                    old_name = namesDict[species2][library.label]
                    logging.warning('Duplicate "{1}" of species "{old}" {0} in library {2}, likely a mismatch or different energy state'.format(species.toSMILES(), name, library.label, old=old_name))
                    previous_thermo = compareDict[species2][library.label]
                    if previous_thermo.getEnthalpy(300) < thermo.getEnthalpy(300) :
                        logging.warning('Using thermo of lower energy version {old} and ignoring subsequent {new}'.format(new=name, old=old_name))
                        break
                    else:
                        logging.warning('Using thermo of lower energy version {new} and dropping initial {old}'.format(new=name, old=old_name))

            else:
                compareDict[species] = {library.label: thermo}
                namesDict[species] = {library.label: name}

    return compareDict, namesDict

def fixKineticsLibrary(fileName):
    """
    Fix malformed kinetics files from previous version of the importer script.
    Should now be deprecated.
    """
    stem,ext = os.path.splitext(fileName)
    fixedFileName = stem+'.fixed'+ext
    line=''
    with open(fileName) as incoming:
        with open(fixedFileName, 'w') as outgoing:
            entry = []
            line = ''
            opening = True
            for line in incoming:
                if line.startswith('entry('):
                    if entry[-2]==')\n' or opening:
                        outgoing.write(''.join(entry))
                        opening = False
                    else:
                        outgoing.write("#"+"#".join(entry))
                    entry = []
                entry.append(line)
    return fixedFileName
    
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

def statsThermo(compareDict, namesDict, book, libraries):
    """
    Outputs a text file with basic statistics on the thermodynamic data collected
    """
    
    # Adjust this to determine the minimum number of species a model must have to be
    # in the output file
    compareLimit = 2
    temp = 298.15 # Temperature in K
    
    def numberModels(species):
        """
        Returns the number of models with data for the given species
        """
        return len(compareDict[species])
    
    # Create a list of all the species sorted by the number of models they are in
    speciesList = list(compareDict)
    speciesList.sort(key = numberModels, reverse = True)
    
    
    # Output to the statistics file
    f = open('thermo.txt', 'w+')
    f.write('Reaction Model Species Comparison File\n')
    f.write('The total number of models is {0}\n'.format(len(libraries)))
    f.write('The temperature is {0} K\n\n'.format(temp))
    f2 = open('thermo_ranges.txt', 'w')
    f3 = open('entropy_ranges.txt', 'w')
    
    
    sheet = book.add_sheet('enthalpy')
    sheet.write(0, 0, 'Reaction Model Comparison File. enthalpy')
    sheet.write(1, 0, 'The total number of models is {0}'.format(len(libraries)))
    sheet.write(2, 0, 'The temperature is {0} K'.format(temp))
    sheet.write(4, 0, 'N')
    sheet.write(4, 1, 'Species')
    sheet.write(4, 2, 'Formula')
    
    Stemp = 1000
    entropy_sheet = book.add_sheet('entropy')
    entropy_sheet.write(0, 0, 'Reaction Model Comparison File. entropy')
    entropy_sheet.write(1, 0, 'The total number of models is {0}'.format(len(libraries)))
    entropy_sheet.write(2, 0, 'The temperature is {0} K'.format(Stemp))
    entropy_sheet.write(4, 0, 'N')
    entropy_sheet.write(4, 1, 'Species')
    entropy_sheet.write(4, 2, 'Formula')

    names_sheet = book.add_sheet('names')
    names_sheet.write(0, 0, 'Reaction Model Comparison File. names')
    names_sheet.write(1, 0, 'The total number of models is {0}'.format(len(libraries)))
    names_sheet.write(4, 0, 'N')
    names_sheet.write(4, 1, 'Species')
    names_sheet.write(4, 2, 'Formula')


    models = {}
    prefix_columns = 3
    stats_columns = len(libraries) + prefix_columns + 1
    for i, library in enumerate(libraries, prefix_columns):
        model = nameFromPath(library)
        models[model] = i
        for s in [sheet, entropy_sheet, names_sheet]:
            s.write(4, i, model)

    for s in [sheet, entropy_sheet]:
        s.write(4, stats_columns, 'Mean')
        s.write(4, stats_columns+1, 'Stdev')
        s.write(4, stats_columns+2, 'Range')
        s.write(4, stats_columns+3, 'Median')

    for i, species in enumerate(speciesList, 5):
        
        if len(compareDict[species]) < compareLimit: 
            break
        
        f.write(str(species.toSMILES())+'\n')
        for s in [sheet, names_sheet, entropy_sheet]:
            s.write(i, 0, numberModels(species))
            s.write(i, 1, species.toSMILES())
            s.write(i, 2, species.getFormula())
        
        stats = []
        entropy_stats = []
        for model in compareDict[species]:
            f.write(model.ljust(35))
            thermo = compareDict[species][model].getEnthalpy(temp)
            entropy = compareDict[species][model].getEntropy(Stemp)
            f.write(str(round(thermo/1000,3)).rjust(15)+'\n')
            entropy_sheet.write(i, models[model], entropy)
            sheet.write(i, models[model], thermo/1000)
            stats.append(thermo/1000)
            entropy_stats.append(entropy)
            names_sheet.write(i, models[model], namesDict[species][model])

        f.write('N{0} | Mean | Std  | Max  | Min\n'.format(' ' if len(compareDict[species])>9 else ''))
        f.write('{0} | {1:.2f} | {2:.2f} | {3:.2f} | {4:.2f}\n\n\n'.format(len(compareDict[species]),
                                                                           numpy.mean(stats), numpy.std(stats), 
                                                                           max(stats), min(stats)))
                
        # rowcol_to_cell(row, col, row_abs=False, col_abs=False):
        first =  rowcol_to_cell(i, 3) # the cell for the first model
        last =  rowcol_to_cell(i, stats_columns-2, col_abs=False) # the cell for the last model
        for s in [sheet, entropy_sheet]:
            s.write(i, stats_columns, Formula('AVERAGE({first}:{last})'.format(first=first, last=last)))
            s.write(i, stats_columns + 1, Formula('STDEV({first}:{last})'.format(first=first, last=last)))
            s.write(i, stats_columns + 2, Formula('MAX({first}:{last}) - MIN({first}:{last})'.format(first=first, last=last)))
            s.write(i, stats_columns + 3, Formula('MEDIAN({first}:{last})'.format(first=first, last=last)))
                
        f2.write("{0:.4f}\t{1}\n".format(max(stats)-min(stats),species.toSMILES()))
        f3.write("{0:.4f}\t{1}\n".format(max(entropy_stats)-min(entropy_stats),species.toSMILES()))

    f.close()
    f2.close()
    f3.close()

def statsKinetics(compareDict, book, libraries):
    """
    Outputs a text file with basic statistics on the kinetics data collected
    """
    
    # Adjust this to determine the minimum number of species a model must have to be
    # in the output file
    compareLimit = 2
    
    
    def numberModels(reaction):
        """
        Returns the number of models with data for the given species
        """
        return len(compareDict[reaction])
    
    # Create a list of all the species sorted by the number of models they are in
    reactionsList = list(compareDict)
    reactionsList.sort(key = numberModels, reverse = True)
    for press in [0.01, 100]:  # Pressure in bar
      for temp in [500, 1000, 1500, 'AnE']:  # Temperature in K
        filename = 'kinetics-{0}K-{1}bar.txt'.format(temp, press)
        sheetname = 'k-{0}K-{1}bar'.format(temp, press)
        if temp == 'AnE':
            if press == 100:
                filename = 'kinetics-expressions.txt'
                sheetname = 'k-expressions'
            else:
                continue
        rangefilename = filename[:-4]+'-range.txt'
        f2 = open(rangefilename,'w')
        # Output to the statistics file
        f = open(filename, 'w+')
        f.write('Reaction Model Comparison File\n')
        f.write('The total number of models is {0}\n'.format(len(libraries)))
        f.write('The temperature is {0} K'.format(temp))
        f.write('The pressure is {0} bar'.format(press))
        f.write('\n\n')
        
        sheet = book.add_sheet(sheetname)
        sheet.write(0, 0, 'Reaction Model Comparison File. {0}'.format(sheetname))
        sheet.write(1, 0, 'The total number of models is {0}'.format(len(libraries)))
        sheet.write(2, 0, 'The temperature is {0} K'.format(temp))
        sheet.write(3, 0, 'The pressure is {0} bar'.format(press))
        sheet.write(4, 0, 'N')
        sheet.write(4, 1, 'Reaction')
        
        models = {}
        for i, library in enumerate(libraries, 2):
            model = nameFromPath(library)
            sheet.write(4, i, model)
            models[model] = i
        stats_columns = len(models) + 3
        
        if temp != 'AnE':
            sheet.write(4, stats_columns, 'Mean')
            sheet.write(4, stats_columns + 1, 'Stdev')
            sheet.write(4, stats_columns + 2, 'Range')
            sheet.write(4, stats_columns + 3, 'Median')

        for i, reaction in enumerate(reactionsList, 5):

            if len(compareDict[reaction]) < compareLimit:
                break

            sheet.write(i, 0, numberModels(reaction))
            sheet.write(i, 1, str(reaction))
            f.write(str(reaction) + '\n')

            stats = []
            for model in compareDict[reaction]:
                f.write(model.ljust(35))
                
                if temp=='AnE':
                    rate = repr(compareDict[reaction][model])
                    rate = rate.decode(errors='replace').encode('ascii', errors='replace')
                    sheet.write(i, models[model], "'{0}".format(rate))
                    f.write(rate +'\n')
                else:
                    try:
                        rate = compareDict[reaction][model].getRateCoefficient(temp, press * 1e5)
                    except:
                        f.write(" Error evaluating rate\n")
                        sheet.write(i, models[model], '#ERROR')
                        continue
                    if rate <= 0:
                        f.write(" k={0:g} so can't take log\n".format(rate))
                        sheet.write(i, models[model], 'LOG10({0:g})'.format(rate))
                        continue
                    f.write(str(round(log10(rate), 3)).rjust(15) + '\n')
                    sheet.write(i, models[model], log10(rate))
                    stats.append(log10(rate))
                    
            if temp=="AnE":
                f.write('\n\n')
            else:
                f.write('N{0} | Mean | Std  | Max  | Min\n'.format(' ' if len(compareDict[reaction])>9 else ''))
                f.write('{0} | {1:.2f} | {2:.2f} | {3:.2f} | {4:.2f}\n\n\n'.format(len(compareDict[reaction]),
                                                                                   numpy.mean(stats), numpy.std(stats),
                                                                                   max(stats), min(stats)))
                # rowcol_to_cell(row, col, row_abs=False, col_abs=False):
                first =  rowcol_to_cell(i, 2) # the cell for the first model
                last =  rowcol_to_cell(i, stats_columns-2, col_abs=False) # the cell for the last model
                sheet.write(i, stats_columns, Formula('AVERAGE({first}:{last})'.format(first=first, last=last)))
                sheet.write(i, stats_columns + 1, Formula('STDEV({first}:{last})'.format(first=first, last=last)))
                sheet.write(i, stats_columns + 2, Formula('MAX({first}:{last}) - MIN({first}:{last})'.format(first=first, last=last)))
                sheet.write(i, stats_columns + 3, Formula('MEDIAN({first}:{last})'.format(first=first, last=last)))
            
            if temp!="AnE":
                f2.write("{0:.4f}\t{1}\n".format(max(stats)-min(stats),str(reaction)))
        f.close()
        f2.close()

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


    book = Workbook()
    thermoDict, namesDict = readThermoLibs(thermoLibs)
    kineticsDict = readKineticsLibs(kineticsLibs)
    statsThermo(thermoDict, namesDict, book, thermoLibs)
    statsKinetics(kineticsDict, book, kineticsLibs)
    
    book.save('output.xls')
    
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

