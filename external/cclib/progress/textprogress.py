# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

from __future__ import print_function
import sys


class TextProgress:

    def __init__(self):

        self.nstep = 0
        self.text = None
        self.oldprogress = 0
        self.progress = 0
        self.calls = 0

    def initialize(self, nstep, text=None):

        self.nstep = float(nstep)
        self.text = text

        #sys.stdout.write("\n")

    def update(self, step, text=None):

        self.progress = int(step * 100 / self.nstep)

        if self.progress/2 >= self.oldprogress/2 + 1 or self.text != text:
        # just went through at least an interval of ten, ie. from 39 to 41,
        # so update

            mystr = "\r["
            prog = int(self.progress / 10)
            mystr += prog * "=" + (10-prog) * "-"
            mystr += "] %3i" % self.progress + "%"

            if text:
                mystr += "    "+text

            sys.stdout.write("\r" + 70 * " ")
            sys.stdout.flush()
            sys.stdout.write(mystr)
            sys.stdout.flush()
            self.oldprogress = self.progress

            if self.progress >= 100 and text == "Done":
                print(" ")

        return
