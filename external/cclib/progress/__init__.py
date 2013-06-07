"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 620 $"


import sys

if 'qt' in sys.modules.keys():
    from qtprogress import QtProgress
if 'PyQt4' in sys.modules.keys():
    from qt4progress import Qt4Progress

from textprogress import TextProgress
