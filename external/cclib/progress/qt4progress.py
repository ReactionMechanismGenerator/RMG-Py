# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

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

