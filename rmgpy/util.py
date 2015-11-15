#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import os.path
import shutil


class Subject(object):
    """Subject in Observer Pattern"""
    def __init__(self):
        self._observers = []

    """
    Call this method when your (self-implemented)
    observer class should start listening to the Subject
    class.

    e.g.:

    listener = YourOwnListener()
    subject.attach(listener)
    """
    def attach(self, observer):
        if not observer in self._observers:
            self._observers.append(observer)


    """
    Call this method when your (self-implemented)
    observer class should stop listening to the Subject
    class.

    e.g.:
    listener = YourOwnListener()
    subject.attach(listener)

    ...<do some work>...

    subject.detach(listener)

    """
    def detach(self, observer):
        try:
            self._observers.remove(observer)
        except ValueError:
            pass

    """
    Call this method in classes that implement
    Subject, when the data that your interested in,
    is available.

    e.g.:
    class YourClass(Subject):
        ...
        def simulate(...)
            <stuff is being done>

            self.notify()

            <continue doing other stuff>            


    Make sure that your listener class implements the update(subject)
    method!

    e.g.:

    class YourOwnListener(object):
        def __init__(self):
            self.data = []

        def update(self, subject):
            self.data.append(subject.data)

    """
    def notify(self, modifier=None):
        for observer in self._observers:
            if modifier != observer:
                observer.update(self)

def makeOutputSubdirectory(outputDirectory, folder):
    """
    Create a subdirectory `folder` in the output directory. If the folder
    already exists (e.g. from a previous job) its contents are deleted.
    """
    dir = os.path.join(outputDirectory, folder)
    if os.path.exists(dir):
        # The directory already exists, so delete it (and all its content!)
        shutil.rmtree(dir)
    os.mkdir(dir)