"""
This scripts runs tests on the database 
"""
import os.path
import logging

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase

def checkFamilies(FullDatabase):
    databaseLog=logging.getLogger('databaseLog')
    familyStatus={}
    for family in FullDatabase.kinetics.families:
        databaseLog.error('\nChecking ' + family + "...")
        print 'Checking ' + family + "..."
        familyStatus[family]=FullDatabase.kinetics.families[family].checkWellFormed()
    
if __name__ == '__main__':
    # Set up paths for database and logger
    databaseDirectory = settings['database.directory']    # RMG-database/input    
    logPath = os.path.join(databaseDirectory, '..', 'database.log')
    #clear logger if it exists
    if os.path.exists(logPath):
        with open(logPath, 'w'): 
            pass
    databaseLog=logging.getLogger('databaseLog')
    fh=logging.FileHandler(logPath)
    fh.setLevel(logging.DEBUG)
    databaseLog.addHandler(fh)
    databaseLog.propagate=False #prevents these logging messages to being sent to ancestors (so that it doesn't print on console)
#     logging.basicConfig(filename=logpath, level=logging.ERROR)

    FullDatabase=RMGDatabase()
    FullDatabase.load(databaseDirectory, kineticsFamilies='all')
    checkFamilies(FullDatabase)

