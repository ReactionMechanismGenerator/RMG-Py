"""
This script machine writes all the families and groups (no libraries) so that RMG-database looks cleaner
and does not have duplicate indexes
"""

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase

database = RMGDatabase()
database.load(settings['database.directory'], kineticsFamilies = 'all')

database.save(settings['database.directory'])