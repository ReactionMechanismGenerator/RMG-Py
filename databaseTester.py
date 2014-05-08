"""
This scripts runs tests on the database 
"""
import os.path
import logging

from rmgpy.data.rmg import RMGDatabase

def checkFamilies(FullDatabase):
    familyStatus={}
    for family in FullDatabase.kinetics.families:
        print 'Checking ' + family + "..."
        familyStatus[family]=FullDatabase.kinetics.families[family].checkWellFormed()
    
    with open(r'C:\RMG-database\DatabaseWellFormedSummary.txt', 'wb') as outputFile:
        for family, problems in familyStatus.iteritems():
            problemsExist=[]
            for problem in problems:
                problemsExist.append(not problem==[] and not problem=={})
            if True in problemsExist:
                outputFile.write(family + '\n')
                if problemsExist[0]:
                    outputFile.write('\n' + 'These groups exist in rules.py but not groups.py:' + '\n' + "A suggested match could be incorrect, but if 'No match' is written, it is true (and most unfortunate)" + '\n')
                    for group, matchedGroups in problems[0].iteritems():
                        outputFile.write(group + ', Suggested match from groups.py: ')
                        for matchedGroup in matchedGroups:
                            if matchedGroup==matchedGroups[-1]:
                                if len(matchedGroups)>1:
                                    outputFile.write('and ')
                                outputFile.write(matchedGroup + '\n')
                            else:
                                outputFile.write(matchedGroup +', ' )
                if problemsExist[1]:
                    outputFile.write('\n' + 'These groups do not match the definition in the rule' + '\n')
                    for rule, groups in problems[1].iteritems():
                        for group in groups:
                            if group==groups[-1]:
                                if len(groups)>1:
                                    outputFile.write('and ')
                                outputFile.write(group + ' ')
                            else:
                                outputFile.write(group +', ' )
                        outputFile.write('in ' + rule + '\n')
                if problemsExist[2]:
                    outputFile.write('\n' + 'These groups are not in the tree:' + '\n')
                    for group in problems[2]:
                        outputFile.write(group + '\n')
                if problemsExist[3]:
                    outputFile.write('\n' + 'These groups are not unique' + '\n')
                    for key, groups in problems[3].iteritems():
                        outputFile.write(key + ' matches ')
                        for group in groups:
                            if group==groups[-1]:
                                if len(groups)>1:
                                    outputFile.write('and ')
                                outputFile.write(group + '\n')
                            else:
                                outputFile.write(group +', ' )
                if problemsExist[4]:
                    outputFile.write('\n' + 'These groups are not actually subgroups of their parent' + '\n')
                    for group, parent in problems[4].iteritems():
                        outputFile.write('Child: ' + group + ', Parent: ' + parent + '\n') 
                if problemsExist[5]:
                    outputFile.write('\n' + 'These groups are probably products, but you should check them anyway' + '\n')
                    for group in problems[5]:
                        outputFile.write(group + '\n')
                outputFile.write('\n\n')

if __name__ == '__main__':
    #set up paths for database and logger
#     databaseProjectRootPath = os.path.dirname( os.path.abspath( __file__ ))
    databaseProjectRootPath = "C:\RMG-database"
    path = os.path.join(databaseProjectRootPath, 'input')
    logpath = os.path.join(databaseProjectRootPath, 'database.log')
    if os.path.exists('database.log'): os.remove('database.log')
    logging.basicConfig(filename=logpath, level=logging.ERROR)
    #load rmg database
    FullDatabase=RMGDatabase()
#     FullDatabase.load(thermoLibraries=)
    FullDatabase.load(path, kineticsFamilies='all')
    checkFamilies(FullDatabase)

