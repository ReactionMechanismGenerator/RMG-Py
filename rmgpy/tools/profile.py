import os
from functools import wraps
import cProfile

from rmgpy.rmg.main import processProfileStats, makeProfileGraph

def profilefn(fn):
    @wraps(fn)
    def profile(*args, **kwargs):
        command = """fn(*args, **kwargs)"""
        stats_file = os.path.join(os.getcwd(), fn.__name__+'.profile')
        cProfile.runctx(command, globals(), locals(), stats_file)

        # postprocess the stats
        
        processProfileStats(stats_file, os.path.join(os.getcwd(), fn.__name__+'.log'))
        makeProfileGraph(stats_file)
        return 
    return profile