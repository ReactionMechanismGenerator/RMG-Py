#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

from rmgpy.chem.pattern import MoleculePattern

################################################################################

class PatternCheck(unittest.TestCase):

    def testPickle(self):
        """
        Test that a MoleculePattern object can be successfully pickled and
        unpickled with no loss of information.
        """
        pattern0 = MoleculePattern().fromAdjacencyList(
        """
        1 *2 {Cs,Cd} 0 {2,{S,D}} {3,S}
        2 *1 {Os,Od} 0 {1,{S,D}}
        3    R!H     0 {1,S}
        """)

        pattern0.updateConnectivityValues()
        import cPickle
        pattern = cPickle.loads(cPickle.dumps(pattern0))

#        print
#        print pattern0.toAdjacencyList()
#        print pattern.toAdjacencyList()

        self.assertEqual(len(pattern0.atoms), len(pattern.atoms))
        for atom0, atom in zip(pattern.atoms, pattern0.atoms):
            self.assertTrue(atom0.equivalent(atom))

        self.assertTrue(pattern0.isIsomorphic(pattern))
        self.assertTrue(pattern.isIsomorphic(pattern0))
        
################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )