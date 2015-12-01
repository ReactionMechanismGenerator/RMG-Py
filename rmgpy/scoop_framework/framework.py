#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
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

"""
This module contains functionality for the correct execution of unit tests 
that use the SCOOP framework.
"""

import unittest
import subprocess
import time
import os
import sys
import signal

try:

    import scoop
    scoop.DEBUG = False

    from scoop import futures, _control, utils, shared
    from scoop._types import FutureQueue
    from scoop.broker.structs import BrokerInfo

    from scoop import logger as logging

except ImportError, e:
    import logging as logging
    logging.debug("Could not properly import SCOOP.")


subprocesses = []
def cleanSubprocesses():
    [a.kill() for a in subprocesses]

try:
    signal.signal(signal.SIGQUIT, cleanSubprocesses)
except AttributeError:
    # SIGQUIT doesn't exist on Windows
    signal.signal(signal.SIGTERM, cleanSubprocesses)


def port_ready(port, socket):
    """Checks if a given port is already binded"""
    try:
        socket.connect(('127.0.0.1', port))
    except IOError:
        return False
    else:
        socket.shutdown(2)
        socket.close()
        return True


class TestScoopCommon(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        # Parent initialization
        super(TestScoopCommon, self).__init__(*args, **kwargs)

    @classmethod
    def setUpClass(cls):
        global subprocesses
        import socket, datetime, time

        # Start the server
        cls.server = subprocess.Popen([sys.executable, "-m", "scoop.broker.__main__",
        "--tPort", "5555", "--mPort", "5556"])
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        begin = datetime.datetime.now()
        while not port_ready(5555, s):
            if (datetime.datetime.now() - begin > datetime.timedelta(seconds=3)):
                raise Exception('Could not start server!')
        subprocesses.append(cls.server)

        # Setup worker environment
        scoop.IS_RUNNING = True
        scoop.IS_ORIGIN = True
        scoop.WORKER_NAME = 'origin'.encode()
        scoop.BROKER_NAME = 'broker'.encode()
        scoop.BROKER = BrokerInfo("127.0.0.1",
                                  5555,
                                  5556,
                                  "127.0.0.1")
        scoop.worker = (scoop.WORKER_NAME, scoop.BROKER_NAME)
        scoop.MAIN_MODULE = "tests.py"
        scoop.VALID = True
        scoop.DEBUG = False
        scoop.SIZE = 1
        _control.execQueue = FutureQueue()


    @classmethod
    def tearDownClass(cls):
        global subprocesses
        import socket, datetime, time
        _control.execQueue.shutdown()
        del _control.execQueue
        _control.futureDict.clear()
        try:
            cls.w.terminate()
            cls.w.wait()
        except:
            pass
        # Destroy the server
        if cls.server.poll() == None:
            try:
                cls.server.terminate()
                cls.server.wait()
            except:
                pass
        # Stabilise zmq after a deleted socket
        del subprocesses[:]

        # Wait for the previous server to be correctly terminated
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        begin = datetime.datetime.now()
        while port_ready(5555, s):
            if (datetime.datetime.now() - begin > datetime.timedelta(seconds=3)):
                raise Exception('Could not terminate server!')
        s.close()