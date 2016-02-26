import unittest
import os.path
import re

from rmgpy.tools.profile import *

class ProfileTest(unittest.TestCase):

    def testFunc(self):
        """
        Test if we decorate a function with @profilefn.
        """
        
        @profilefn
        def func():
            pass

        func()

        module = sys.modules[func.__module__]
        dirname = os.path.join(os.path.dirname(module.__file__))

        log = os.path.join(dirname, func.__name__ + '.log')
        prf = os.path.join(dirname, func.__name__ + '.profile')
        dot = os.path.join(dirname, func.__name__ + '.profile.dot')
        pdf = os.path.join(dirname, func.__name__ + '.profile.dot.pdf')
        ps2 = os.path.join(dirname, func.__name__ + '.profile.dot.ps2')

        for f in [log, prf, dot, pdf, ps2]:
            self.assertTrue(os.path.isfile(f))

        for f in os.listdir(dirname):
            if re.search(func.__name__, f):
                os.remove(os.path.join(dirname, f))