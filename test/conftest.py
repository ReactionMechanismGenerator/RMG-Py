import multiprocessing

multiprocessing.set_start_method('fork')

from openbabel import pybel

pybel.ob.obErrorLog.SetOutputLevel(0)
