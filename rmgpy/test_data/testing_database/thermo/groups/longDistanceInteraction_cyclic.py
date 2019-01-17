#!/usr/bin/env python
# encoding: utf-8

name = "Long distance interaction correction for cyclic molecule"
shortDesc = u""
longDesc = u"""
Designed to account for the long distance interaction for cyclic molecule.
Currently we only have the data for aromatic ortho/meta/para corrections.
Watch out:  if the groups on the two labeled atoms are identical, it's value should be halved because it'll be counted twice.
            For example, for the interaction of ortho OH and OH, which is labeled as 'o_OH_OH' in this database,
            it'll be counted in both {'*1': atom1, '*2'atom2} and {'*2': atom1, '*1'atom2}.
            It should be claimed in the 'longDesc' if a entry was halved.

Source: 
For aromatic molecule: [1] Ince et al., AIChE 2015, DOI 10.1002/aic.15008
For aromatic radical: [2] Ince et al., AIChE 2016, DOI 10.1002/aic.15588

Jan-23-2017 PZ
"""


