#!/usr/bin/env python
# -*- coding: utf-8 -*-

energy = {
    'CBS-QB3': Log('composite.out'),
}

geometry = Log('freq.out')

frequencies = Log('freq.out')

"""
2D rotors require the pivots and tops (much like regular hindered rotors) for the 2 involved rotors
they also require a scandir that has scans at regular angle differences in (phi1,phi2) angle space
the scans output file names must be formatted in  {}_{}_{phi1}_{phi2}.out where {} can be anything without a '_'
in it and phi1 and phi2 are in degrees.    
symmetry numbers for 2D rotors must be specified as they are not auto-detected
Q2DTor overall symmetry type specification ('symmetry') can speed up calculations a bit if they are known 
otherwise they default to 'none'.  
"""
rotors = [HinderedRotor2D(scandir='CH2CHOOHscans', pivots1=[6,7], pivots2=[1,6], top1=[7,8],
                          top2=[6,7,8], symmetry1=1, symmetry2=1, symmetry='none')]

