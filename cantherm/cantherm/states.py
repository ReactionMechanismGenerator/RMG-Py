#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   CanTherm
#    
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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

import math
import numpy
import logging

import chempy.constants as constants
from chempy.states import *

################################################################################

def applyEnergyCorrections(E0, modelChemistry, atoms, bonds):
    """
    Given an energy `E0` in J/mol as read from the output of a quantum chemistry
    calculation at a given `modelChemistry`, adjust the energy such that it
    is consistent with the normal gas-phase reference states. `atoms` is a
    dictionary associating element symbols with the number of that element in
    the molecule. `bonds` is a dictionary associating bond types with the number
    of that bond in the molecule.
    """
    
    # Step 1: Reference all energies to a model chemistry-independent basis
    # by subtracting out that model chemistry's atomic energies
    if modelChemistry == 'CBS-QB3':
        atomEnergies = {'H':-0.499818 , 'N':-54.520543, 'O':-74.987624, 'C':-37.785385, 'P':-340.817186}
    elif modelChemistry == 'G3':
        atomEnergies = {'H':-0.5010030, 'N':-54.564343, 'O':-75.030991, 'C':-37.827717, 'P':-341.116432}
    else:
        logging.warning('Unknown model chemistry "%s"; not applying energy corrections.' % modelChemistry)
        return E0
    for symbol, count in atoms.iteritems():
        if symbol in atomEnergies: E0 -= count * atomEnergies[symbol] * 4.35974394e-18 * constants.Na
        else:
            logging.warning('Ignored unknown atom type "%s".' % symbol)
    
    # Step 2: Atom energy corrections to reach gas-phase reference state
    # Experimental number for H includes H + TC + SOC (SOC = spin-orbit coupling)
    atomEnergies = {'H': 50.62, 'N': 111.49, 'O': 58.163, 'C': 169.8147}
    for symbol, count in atoms.iteritems():
        if symbol in atomEnergies: E0 += count * atomEnergies[symbol] * 4184
    
    # Step 3: Bond energy corrections
    bondEnergies = { 'C-H': -0.11, 'C-C': -0.3, 'C=C': -0.08, 'C#C': -0.64,
        'O-H': 0.02, 'C-O': 0.33, 'C=O': 0.55, 'N#N': -2.0, 'O=O': -0.2, 
        'H-H': 1.1, 'C#N': -0.89 }
    for symbol, count in bonds.iteritems():
        if symbol in bondEnergies: E0 += count * bondEnergies[symbol] * 4184
        else:
            logging.warning('Ignored unknown bond type "%s".' % symbol)
    
    return E0

################################################################################

def projectRotors(geom, F, rotors, linear, TS):
    """
    For a given geometry `geom` with associated force constant matrix `F`,
    lists of rotor information `rotors`, `pivots`, and `top1`, and the linearity
    of the molecule `linear`, project out the nonvibrational modes from the
    force constant matrix and use this to determine the vibrational frequencies.
    The list of vibrational frequencies is returned in cm^-1.
    """
    
    Nrotors = len(rotors)
    Natoms = len(geom.mass)
    Nvib = 3 * Natoms - (5 if linear else 6) - Nrotors - (1 if (TS) else 0)
    
    if linear:
        D = numpy.zeros((Natoms*3,5+Nrotors), numpy.float64)
    else:
        D = numpy.zeros((Natoms*3,6+Nrotors), numpy.float64)

    for i in range(Natoms):
        # Projection vectors for translation
        D[3*i+0,0] = 1.0
        D[3*i+1,1] = 1.0
        D[3*i+2,2] = 1.0
        # Projection vectors for [external] rotation
        D[3*i:3*i+3,3] = numpy.array([0, -geom.coordinates[i,2], geom.coordinates[i,1]], numpy.float64)
        D[3*i:3*i+3,4] = numpy.array([geom.coordinates[i,2], 0, -geom.coordinates[i,0]], numpy.float64)
        if not linear:
            D[3*i:3*i+3,5] = numpy.array([-geom.coordinates[i,1], geom.coordinates[i,0], 0], numpy.float64)
    for i, rotor in enumerate(rotors):
        scanLog, pivots, top, symmetry = rotor
        # Determine pivot atom
        if pivots[0] in top: pivot = pivots[0]
        elif pivots[1] in top: pivot = pivots[1]
        else: raise Exception('Could not determine pivot atom.')
        # Projection vectors for internal rotation
        e12 = geom.coordinates[pivots[0],:] - geom.coordinates[pivots[1],:]
        e12 /= numpy.linalg.norm(e12)
        for atom in top:
            e31 = geom.coordinates[atom,:] - geom.coordinates[pivot,:]
            D[3*atom:3*atom+3,-Nrotors+i] = numpy.cross(e31, e12)

    # Make sure projection matrix is orthonormal
    import scipy.linalg
    D = scipy.linalg.orth(D)

    # Project out the non-vibrational modes from the force constant matrix
    P = numpy.dot(D, D.transpose())
    I = numpy.identity(Natoms*3, numpy.float64)
    F = numpy.dot(I - P, numpy.dot(F, I - P))

    # Generate mass-weighted force constant matrix
    # This converts the axes to mass-weighted Cartesian axes
    # Units of Fm are J/m^2*kg = 1/s^2
    Fm = F.copy()
    for i in range(Natoms):
        for j in range(Natoms):
            for u in range(3):
                for v in range(3):
                    Fm[3*i+u,3*j+v] /= math.sqrt(geom.mass[i] * geom.mass[j]) / constants.Na

    # Get eigenvalues of mass-weighted force constant matrix
    eig, V = numpy.linalg.eigh(Fm)
    eig.sort()

    # Convert eigenvalues to vibrational frequencies in cm^-1
    # Only keep the modes that don't correspond to translation, rotation, or internal rotation
    return numpy.sqrt(eig[-Nvib:]) / (2 * math.pi * constants.c * 100)

################################################################################   

def saveStates(species, label, path):
    """
    Append the molecular degrees of freedom for `species` with associated
    string `label` to the file located at `path` on disk.
    """
    
    f = open(path, 'a')
    f.write('states(\n')
    f.write('    label = "%s",\n' % label)
    
    f.write('    modes = [\n')
    for mode in species.states.modes:
        if isinstance(mode, Translation):
            f.write('        Translation(mass=(%g,"g/mol")),\n' % (mode.mass*1000))
        elif isinstance(mode, RigidRotor):
            f.write('        RigidRotor(linear=%s, inertia=[%s], "amu*angstrom^2"), symmetry=%i),\n' % (mode.linear, ', '.join(['%g' % (I * 6.022e46) for I in mode.inertia]), mode.symmetry))
        elif isinstance(mode, HarmonicOscillator):
            f.write('        HarmonicOscillator(frequencies=([%s], "cm^-1")),\n' % (', '.join(['%g' % (freq) for freq in mode.frequencies])))
        elif isinstance(mode, HinderedRotor):
            f.write('        HinderedRotor(inertia=(%g, "amu*angstrom^2"), symmetry=%i, fourier=[[%g, %g, %g, %g, %g], [%g, %g, %g, %g, %g]]),\n' % (mode.inertia * 6.022e46, mode.symmetry, 
                mode.fourier[0,0], mode.fourier[0,1], mode.fourier[0,2], mode.fourier[0,3], mode.fourier[0,4],
                mode.fourier[1,0], mode.fourier[1,1], mode.fourier[1,2], mode.fourier[1,3], mode.fourier[1,4],
                ))
    f.write('    ],\n')
    try:
        f.write('    frequency=(%g,"cm^-1"),\n' % species.frequency)
    except AttributeError: pass
    f.write('    short_comment = "",\n')
    f.write('    long_comment = \n')
    f.write('"""\n')
    f.write('\n')
    f.write('""",\n')
    f.write(')\n\n')
    
    f.close()
