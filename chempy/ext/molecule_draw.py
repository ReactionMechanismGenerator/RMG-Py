#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
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

"""
This module provides functionality for automatic drawing of two-dimensional
molecules.
"""

import math
import numpy
import os.path

import sys
sys.path.append(os.path.abspath('.'))
from chempy.molecule import *

################################################################################

def render(atoms, bonds, coordinates, fstr):
    """
    Uses the Cairo graphics library to create the drawing of the molecule.
    """

    try:
        import cairo
    except ImportError:
        print 'Cairo not found; potential energy surface will not be drawn.'
        return

    # Initialize Cairo surface and context
    width = 640; height = 480
    ext = os.path.splitext(fstr)[1].lower()
    if ext == '.svg':
        surface = cairo.SVGSurface(fstr, width, height)
    elif ext == '.pdf':
        surface = cairo.PDFSurface(fstr, width, height)
    elif ext == '.ps':
        surface = cairo.PSSurface(fstr, width, height)
    cr = cairo.Context(surface)

    bondLength = 32
    coordinates[:,1] *= -1
    coordinates = coordinates * bondLength + 160

    # Some global settings
    cr.select_font_face("sans")
    cr.set_font_size(10)

    # Draw bond skeleton (for now)
    cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
    cr.set_line_width(1.0)
    for atom1 in bonds:
        for atom2, bond in bonds[atom1].iteritems():
            index1 = atoms.index(atom1)
            index2 = atoms.index(atom2)
            if index1 < index2:
                x1, y1 = coordinates[index1,:]
                x2, y2 = coordinates[index2,:]
                cr.move_to(x1, y1)
                cr.line_to(x2, y2)
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.stroke()
                # Text label (for now)
                labels = {1:'S', 2:'D', 3:'T'}
                label = labels[bond.order]
                extents = cr.text_extents(label)
                x0 = 0.5 * (x1 + x2); y0 = 0.5 * (y1 + y2)
                cr.move_to(x0 - extents[0] - extents[2] / 2.0, y0 - extents[1] - extents[3] / 2.0)
                cr.set_source_rgba(0.0, 0.0, 1.0, 1.0)
                cr.set_font_size(8)
                cr.show_text(label)
                cr.set_font_size(10)

    # Draw atoms (for now)
    for atom in atoms:
        symbol = atom.symbol
        index = atoms.index(atom)
        x0, y0 = coordinates[index,:]
        extents = cr.text_extents(symbol)
        width = extents[2]; height = extents[3]
        # Background
        cr.rectangle(x0 - width / 2.0 - 2.0, y0 - height / 2.0 - 2.0, width + 4.0, height + 4.0)
        cr.set_source_rgba(1.0, 1.0, 1.0, 1.0)
        cr.fill()
        # Text itself
        cr.move_to(x0 - extents[0] - width / 2.0, y0 - extents[1] - height / 2.0)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.show_text(symbol)


    # Finish Cairo drawing
    surface.finish()

################################################################################

def findLongestPath(chemGraph, atoms):
    """
    Finds the longest path containing the list of `atoms` in the `chemGraph`.
    The atoms are assumed to already be in a path, with ``atoms[0]`` being a
    terminal atom.
    """
    atom1 = atoms[-1]
    paths = [atoms]
    for atom2 in chemGraph.bonds[atom1]:
        if atom2 not in atoms:
            atoms.append(atom2)
            paths.append(findLongestPath(chemGraph, atoms))
            atoms = atoms[:-1]
    lengths = [len(path) for path in paths]
    index = lengths.index(max(lengths))
    return paths[index]

def getBackbone(chemGraph):
    """
    Return the atoms that make up the backbone of the molecule.
    """

    if chemGraph.isCyclic():
        raise NotImplementedError('Currently cannot find backbone of cyclic molecules!')
    else:
        # Make a shallow copy of the chemGraph so we don't modify the original
        chemGraph = chemGraph.copy()

        # Remove hydrogen atoms from consideration, as they cannot be part of
        # the backbone
        chemGraph.makeHydrogensImplicit()

        # If there are only one or two atoms remaining, these are the backbone
        if len(chemGraph.atoms) == 1 or len(chemGraph.atoms) == 2:
            return chemGraph.atoms[:]

        # Find the terminal atoms - those that only have one explicit bond
        terminalAtoms = []
        for atom in chemGraph.atoms:
            if len(chemGraph.bonds[atom]) == 1:
                terminalAtoms.append(atom)

        # Starting from each terminal atom, find the longest straight path to
        # another terminal; this defines the backbone
        backbone = []
        for atom in terminalAtoms:
            path = findLongestPath(chemGraph, [atom])
            if len(path) > len(backbone):
                backbone = path

        return backbone

################################################################################

def generateCoordinates(chemGraph, atoms, bonds):
    """
    Generate the 2D coordinates to be used when drawing the `chemGraph`, a
    :class:`ChemGraph` object. Use the `atoms` parameter to pass a list
    containing the atoms in the molecule for which coordinates are needed. If
    you don't specify this, all atoms in the molecule will be used.

    Because we are working in two dimensions, we call the horizontal direction
    :math:`u` and the vertical direction :math:`v`. The vertices are arranged
    based on a standard bond length of unity, and can be scaled later for
    longer bond lengths.

    This function ignores any previously-existing coordinate information.
    """

    # Initialize array of coordinates
    coordinates = numpy.zeros((len(atoms), 2), numpy.float64)

    # If there are only one or two atoms to draw, then determining the
    # coordinates is trivial
    if len(atoms) == 1:
        return coordinates
    elif len(atoms) == 2:
        coordinates[1,:] = [1.0, 0.0]
        return coordinates

    # Find the backbone of the molecule
    backbone = getBackbone(chemGraph)

    # Generate coordinates for atoms in backbone
    if chemGraph.isCyclic():
        raise NotImplementedError('Currently cannot find backbone of cyclic molecules!')
    else:
        # Straight chain backbone

        # First atom in backbone goes at origin
        index0 = atoms.index(backbone[0])
        coordinates[index0,:] = [0.0, 0.0]

        # Second atom in backbone goes on x-axis (for now; this could be improved!)
        index1 = atoms.index(backbone[1])
        vector = numpy.array([1.0, 0.0], numpy.float64)
        if bonds[backbone[0]][backbone[1]].isTriple():
            rotatePositive = False
        else:
            rotatePositive = True
            rot = numpy.array([[math.cos(-math.pi / 6), math.sin(-math.pi / 6)], [-math.sin(-math.pi / 6), math.cos(-math.pi / 6)]], numpy.float64)
            vector = numpy.array([1.0, 0.0], numpy.float64)
            vector = numpy.dot(rot, vector)
        coordinates[index1,:] = coordinates[index0,:] + vector

        # Other atoms in backbone
        for i in range(2, len(backbone)):
            atom1 = backbone[i-1]
            atom2 = backbone[i]
            index1 = atoms.index(atom1)
            index2 = atoms.index(atom2)
            bond0 = bonds[backbone[i-2]][atom1]
            bond = bonds[atom1][atom2]
            # Angle of next bond depends on the number of bonds to the start atom
            numBonds = len(bonds[atom1])
            if numBonds == 2:
                if (bond0.isTriple() or bond.isTriple()) or (bond0.isDouble() and bond.isDouble()):
                    # Rotate by 0 degrees towards horizontal axis (to get angle of 180)
                    angle = 0.0
                else:
                    # Rotate by 60 degrees towards horizontal axis (to get angle of 120)
                    angle = math.pi / 3
            elif numBonds == 3:
                # Rotate by 60 degrees towards horizontal axis (to get angle of 120)
                angle = math.pi / 3
            elif numBonds == 4:
                # Rotate by 90 degrees towards horizontal axis (to get angle of 90)
                angle = math.pi / 2
            elif numBonds == 5:
                # Rotate by 36 degrees towards horizontal axis (to get angle of 144)
                angle = math.pi / 5
            elif numBonds == 6:
                # Rotate by 0 degrees towards horizontal axis (to get angle of 180)
                angle = 0.0
            # Determine coordinates for atom
            if angle != 0:
                if not rotatePositive: angle = -angle
                rot = numpy.array([[math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]], numpy.float64)
                vector = numpy.dot(rot, vector)
                rotatePositive = not rotatePositive
            coordinates[index2,:] = coordinates[index1,:] + vector

    # Center backbone at origin
    origin = numpy.zeros(2, numpy.float64)
    for atom in backbone:
        index = atoms.index(atom)
        origin += coordinates[index,:]
    origin /= len(backbone)
    coordinates -= origin

    return coordinates

################################################################################

def drawMolecule(chemGraph, fstr=''):
    """
    Primary function for generating a drawing of a :class:`ChemGraph` object
    `chemGraph`. The parameter `fstr` is the name of a file to save the drawing
    to; valid file extensions are ``.pdf``, ``.svg``, and ``.ps``.
    """

    atoms = chemGraph.atoms[:]
    bonds = chemGraph.bonds.copy()

    # Special cases: H, H2, anything with one heavy atom

    # Remove all unlabeled hydrogen atoms from the molecule, as they are not drawn
    atomsToRemove = []
    for atom in atoms:
        if atom.isHydrogen() and atom.label == '': atomsToRemove.append(atom)
    for atom in atomsToRemove:
        atoms.remove(atom)
    for atom in bonds:
        if atom not in atoms: del bonds[atom]
        for atom2 in bonds[atom]:
            if atom2 not in atoms: del bonds[atom][atom2]

    # Generate the coordinates to use to draw the molecule
    coordinates = generateCoordinates(chemGraph, atoms, bonds)

    # Render using Cairo
    render(atoms, bonds, coordinates, fstr)

################################################################################

if __name__ == '__main__':

    molecule = Molecule()
    molecule.fromSMILES('C=CC=CCC')
    #molecule.fromSMILES('CCC(C)CCC(CCC)C')
    #molecule.fromSMILES('C=CC(C)=CCC')
    #molecule.fromSMILES('COC(C)(C)C(C)(C)N(C)C')
    #molecule.fromSMILES('CCC=C=CCCC')
    #molecule.fromSMILES('C1CCCCC1CCC2CCCC2')

    drawMolecule(molecule.resonanceForms[0], 'molecule.svg')

    