"""
cclib (http://cclib.sf.net) is (c) 2006, the cclib development team
and licensed under the LGPL (http://www.gnu.org/copyleft/lgpl.html).
"""

__revision__ = "$Revision: 238 $"


from PyQt4 import QtGui,QtCore


class Qt4Progress(QtGui.QProgressDialog):

    def __init__(self, title, parent=None):

        QtGui.QProgressDialog.__init__(self, parent)

        self.nstep = 0
        self.text = None
        self.oldprogress = 0
        self.progress = 0
        self.calls = 0
        self.loop=QtCore.QEventLoop(self)
        self.setWindowTitle(title)

    def initialize(self, nstep, text=None):

        self.nstep = nstep
        self.text = text
        self.setRange(0,nstep)
        if text:
            self.setLabelText(text)
        self.setValue(1)
        #sys.stdout.write("\n")

    def update(self, step, text=None):

        if text:
            self.setLabelText(text)
        self.setValue(step)
        self.loop.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)

