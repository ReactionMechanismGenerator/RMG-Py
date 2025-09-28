"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
gmagoon 4/5/10-4/6/10 (this notice added 4/29/10): Gregory Magoon modified this file from cclib 1.0
"""

__revision__ = "$Revision: 863 $"

# These import statements are added for the convenience of users...

# Rather than having to type:
#         from cclib.parser.gaussianparser import Gaussian
# they can use:
#         from cclib.parser import Gaussian

from adfparser import ADF
from gamessparser import GAMESS
from gamessukparser import GAMESSUK
from gaussianparser import Gaussian
from jaguarparser import Jaguar
from molproparser import Molpro
from orcaparser import ORCA
from mopacparser import Mopac
from mm4parser import MM4

# This allow users to type:
#         from cclib.parser import ccopen

from ccopen import ccopen
