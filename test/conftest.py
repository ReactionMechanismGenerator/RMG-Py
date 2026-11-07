import multiprocessing

multiprocessing.set_start_method('fork')

# Load CoolProp before openbabel, or importing CoolProp later hangs test collection;
# see rmgpy/molecule/converter.py
import CoolProp  # noqa: F401
from openbabel import pybel

pybel.ob.obErrorLog.SetOutputLevel(0)
