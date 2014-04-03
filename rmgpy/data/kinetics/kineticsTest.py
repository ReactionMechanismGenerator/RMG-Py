import os
import unittest 
from external.wip import work_in_progress

from rmgpy import settings
from rmgpy.data.kinetics import *
from rmgpy.data.base import DatabaseError
###################################################

class TestKineticsDatabase(unittest.TestCase):
    
    def testLoadFamilies(self):
        """
        Test that the loadFamilies function raises the correct exceptions
        """
        path = os.path.join(settings['database.directory'],'kinetics','families')
        database = KineticsDatabase()
        
        with self.assertRaises(DatabaseError):
            database.loadFamilies(path, families='random')
        with self.assertRaises(DatabaseError):
            database.loadFamilies(path, families=['!H_Abstraction','Disproportionation'])
        with self.assertRaises(DatabaseError):
            database.loadFamilies(path, families=['fake_family'])
        with self.assertRaises(DatabaseError):
            database.loadFamilies(path, families=[])