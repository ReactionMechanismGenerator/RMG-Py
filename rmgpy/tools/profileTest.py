import unittest
import os.path
import re

from rmgpy.tools.profile import *

class ProfileTest(unittest.TestCase):

    def testFunc(self):
        """
        Test if we decorate a function with @profilefn.
        """
        
        name = 'func'

        @profiler
        def func():
            pass

        func()

        module = sys.modules[func.__module__]
        dirname = os.path.join(os.path.dirname(module.__file__))

        prf = os.path.join(dirname, name + '.profile')
        dot = os.path.join(dirname, name + '.profile.dot')
        pdf = os.path.join(dirname, name + '.profile.dot.pdf')
        ps2 = os.path.join(dirname, name + '.profile.dot.ps2')

        for f in [prf, dot, pdf, ps2]:
            self.assertTrue(os.path.isfile(f), f)

        for f in os.listdir(dirname):
            if re.search(name, f):
                os.remove(os.path.join(dirname, f))

    def testCalledFunc(self):
        """
        Test if we decorate a function with @profilefn
        that is called by another function.
        """
        
        name = 'called'

        @profiler
        def called():
            pass

        def func():
            for i in range(3):
                called()

        func()

        module = sys.modules[func.__module__]
        dirname = os.path.join(os.path.dirname(module.__file__))

        prf = os.path.join(dirname, name+ '.profile')
        dot = os.path.join(dirname, name+ '.profile.dot')
        pdf = os.path.join(dirname, name+ '.profile.dot.pdf')
        ps2 = os.path.join(dirname, name+ '.profile.dot.ps2')

        for f in [prf, dot, pdf, ps2]:
            self.assertTrue(os.path.isfile(f), f)

        for f in os.listdir(dirname):
            if re.search(name, f):
                os.remove(os.path.join(dirname, f))