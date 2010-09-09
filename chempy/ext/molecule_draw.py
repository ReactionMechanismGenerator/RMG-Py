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
This module provides functionality for automatic two-dimensional drawing of the
`skeletal formulae <http://en.wikipedia.org/wiki/Skeletal_formula>`_ of a wide
variety of organic and inorganic molecules. The general method for creating
these drawings is to utilize the :meth:`draw()` method of the :class:`Molecule`
or :class:`ChemGraph` you wish to draw; this wraps a call to 
:meth:`drawMolecule()`, where the molecule drawing algorithm begins. Advanced
use may require calling of the :meth:`drawMolecule()` method directly.

The `Cairo <http://cairographics.org/>`_ 2D graphics library is used to create
the drawings. The :meth:`drawMolecule()` method module will fail gracefully if
Cairo is not installed.

The general procedure for creating drawings of skeletal formula is as follows:

1.  **Find the molecular backbone.** If the molecule contains no cycles, the
    longest straight chain of heavy atoms is used as the backbone. If the 
    molecule contains cycles, the largest independent cycle group is used as the
    backbone. The :meth:`findBackbone()` method is used for this purpose.

2.  **Generate coordinates for the backbone atoms.** Straight-chain backbones
    are laid out in a horizontal seesaw pattern. Cyclic backbones are laid out
    as regular polygons (or as close to this as is possible). The
    :meth:`generateStraightChainCoordinates()` and 
    :meth:`generateRingSystemCoordinates()` methods are used for this purpose.

3.  **Generate coordinates for immediate neighbors to backbone.** Each neighbor
    atom represents the start of a functional group attached to the backbone.
    Generating coordinates for these means that we have determined the bonds
    for all backbone atoms. The :meth:`generateNeighborCoordinates()` method is
    used for this purpose.

4.  **Continue generating coordinates for atoms in functional groups.** Moving
    away from the molecular backbone and its immediate neighbors, the
    coordinates for each atom in each functional group are determined such that
    the functional groups tend to radiate away from the center of the backbone
    (to reduce chances of overlap). If cycles are encountered in the functional
    groups, their coordinates are processed as a unit. This continues until
    the coordinates of all atoms in the molecule have been assigned. The
    :meth:`generateFunctionalGroupCoordinates()` recursive method is used for
    this.

5.  **Use the generated coordinates and the atom and bond types to render the
    skeletal formula.** The :meth:`render()`,  and :meth:`renderBond()`, and
    :meth:`renderAtom()` methods are used for this.

The developed procedure seems to be rather robust, but occasionally it will
encounter a molecule that it renders incorrectly. In particular, features which
have not yet been implemented by this drawing algorithm include:

* cis-trans isomerism

* stereoisomerism

* bridging atoms in fused rings

"""

import math
import numpy
import os.path
import re

from chempy.molecule import *

################################################################################

# Parameters that control the Cairo output
fontFamily = 'sans'
fontSizeNormal = 10
fontSizeSubscript = 6
bondLength = 24
    
################################################################################

class MoleculeRenderError(Exception): pass

################################################################################

def render(atoms, bonds, coordinates, symbols, cr, offset=(0,0)):
    """
    Uses the Cairo graphics library to create a skeletal formula drawing of a
    molecule containing the list of `atoms` and dict of `bonds` to be drawn.
    The 2D position of each atom in `atoms` is given in the `coordinates` array.
    The symbols to use at each atomic position are given by the list `symbols`.
    You must specify the Cairo context `cr` to render to.
    """

    import cairo

    # Adjust coordinates such that the top left corner is (0,0) and determine
    # the bounding rect for the molecule
    # Find the atoms on each edge of the bounding rect
    sorted = numpy.argsort(coordinates[:,0])
    left = sorted[0]; right = sorted[-1]
    sorted = numpy.argsort(coordinates[:,1])
    top = sorted[0]; bottom = sorted[-1]
    # Get rough estimate of bounding box size using atom coordinates
    left = coordinates[left,0] + offset[0]
    top = coordinates[top,1] + offset[1]
    right = coordinates[right,0] + offset[0]
    bottom = coordinates[bottom,1] + offset[1]
    # Shift coordinates by offset value
    coordinates[:,0] += offset[0]
    coordinates[:,1] += offset[1]
    
    # Draw bonds
    for atom1 in bonds:
        for atom2, bond in bonds[atom1].iteritems():
            index1 = atoms.index(atom1)
            index2 = atoms.index(atom2)
            if index1 < index2: # So we only draw each bond once
                renderBond(index1, index2, bond, coordinates, symbols, cr)

    # Draw atoms
    for i, atom in enumerate(atoms):
        symbol = symbols[i]
        index = atoms.index(atom)
        x0, y0 = coordinates[index,:]
        vector = numpy.zeros(2, numpy.float64)
        if atom in bonds:
            for atom2 in bonds[atom]:
                vector += coordinates[atoms.index(atom2),:] - coordinates[index,:]
        heavyFirst = vector[0] <= 0
        if len(atoms) == 1 and not atoms[0].symbol not in ['C', 'N'] and atoms[0].charge == 0 and atoms[0].radicalElectrons == 0:
            # This is so e.g. water is rendered as H2O rather than OH2
            heavyFirst = False
            cr.set_font_size(fontSizeNormal)
            x0 += cr.text_extents(symbols[0])[2] / 2.0
        atomBoundingRect = renderAtom(symbol, atom, coordinates, atoms, bonds, x0, y0, cr, heavyFirst)
        # Update bounding rect to ensure atoms are included
        if atomBoundingRect[0] < left:
            left = atomBoundingRect[0]
        if atomBoundingRect[1] < top:
            top = atomBoundingRect[1]
        if atomBoundingRect[2] > right:
            right = atomBoundingRect[2]
        if atomBoundingRect[3] > bottom:
            bottom = atomBoundingRect[3]
    
    # Add a small amount of whitespace on all sides
    padding = 2
    left -= padding; top -= padding; right += padding; bottom += padding

    # Return a tuple containing the bounding rectangle for the drawing
    return (left, top, right-left, bottom-top)

################################################################################

def renderBond(atom1, atom2, bond, coordinates, symbols, cr):
    """
    Render an individual `bond` between atoms with indices `atom1` and `atom2`
    on the Cairo context `cr`.
    """

    import cairo

    cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
    cr.set_line_width(1.0)
    cr.set_line_cap(cairo.LINE_CAP_ROUND)

    x1, y1 = coordinates[atom1,:]
    x2, y2 = coordinates[atom2,:]
    angle = math.atan2(y2 - y1, x2 - x1)

    dx = x2 - x1; dy = y2 - y1
    du = math.cos(angle + math.pi / 2)
    dv = math.sin(angle + math.pi / 2)
    if bond.isDouble() and (symbols[atom1] != '' or symbols[atom2] != ''):
        # Draw double bond centered on bond axis
        du *= 2; dv *= 2
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.move_to(x1 - du, y1 - dv); cr.line_to(x2 - du, y2 - dv)
        cr.stroke()
        cr.move_to(x1 + du, y1 + dv); cr.line_to(x2 + du, y2 + dv)
        cr.stroke()
    elif bond.isTriple() and (symbols[atom1] != '' or symbols[atom2] != ''):
        # Draw triple bond centered on bond axis
        du *= 3; dv *= 3
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.move_to(x1 - du, y1 - dv); cr.line_to(x2 - du, y2 - dv)
        cr.stroke()
        cr.move_to(x1, y1); cr.line_to(x2, y2)
        cr.stroke()
        cr.move_to(x1 + du, y1 + dv); cr.line_to(x2 + du, y2 + dv)
        cr.stroke()
    else:
        # Draw bond on skeleton
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.move_to(x1, y1); cr.line_to(x2, y2)
        cr.stroke()
        # Draw other bonds
        if bond.isDouble():
            du *= 4; dv *= 4; dx = 4 * dx / bondLength; dy = 4 * dy / bondLength
            cr.move_to(x1 + du + dx, y1 + dv + dy); cr.line_to(x2 + du - dx, y2 + dv - dy)
            cr.stroke()
        elif bond.isTriple():
            du *= 3; dv *= 3; dx = 3 * dx / bondLength; dy = 3 * dy / bondLength
            cr.move_to(x1 - du + dx, y1 - dv + dy); cr.line_to(x2 - du - dx, y2 - dv - dy)
            cr.stroke()
            cr.move_to(x1 + du + dx, y1 + dv + dy); cr.line_to(x2 + du - dx, y2 + dv - dy)
            cr.stroke()
    
################################################################################

def renderAtom(symbol, atom, coordinates0, atoms, bonds, x0, y0, cr, heavyFirst=True):
    """
    Render the `label` for an atom centered around the coordinates (`x0`, `y0`)
    onto the Cairo context `cr`. If `heavyFirst` is ``False``, then the order
    of the atoms will be reversed in the symbol. This method also causes
    radical electrons and charges to be drawn adjacent to the rendered symbol.
    """

    import cairo

    if symbol != '':
        heavyAtom = symbol[0]

        # Split label by atoms
        labels = re.findall('[A-Z][0-9]*', symbol)
        if not heavyFirst: labels.reverse()
        symbol = ''.join(labels)

        # Determine positions of each character in the symbol
        coordinates = []

        cr.set_font_size(fontSizeNormal)
        y0 += max([cr.text_extents(char)[3] for char in symbol if char.isalpha()]) / 2

        for i, label in enumerate(labels):
            for j, char in enumerate(label):
                cr.set_font_size(fontSizeSubscript if char.isdigit() else fontSizeNormal)
                xbearing, ybearing, width, height, xadvance, yadvance = cr.text_extents(char)
                if i == 0 and j == 0:
                    # Center heavy atom at (x0, y0)
                    x = x0 - width / 2.0 - xbearing
                    y = y0
                else:
                    # Left-justify other atoms (for now)
                    x = x0
                    y = y0
                if char.isdigit(): y += height / 2.0
                coordinates.append((x,y))
                x0 = x + xadvance

        x = 1000000; y = 1000000; width = 0; height = 0
        startWidth = 0; endWidth = 0
        for i, char in enumerate(symbol):
            cr.set_font_size(fontSizeSubscript if char.isdigit() else fontSizeNormal)
            extents = cr.text_extents(char)
            if coordinates[i][0] + extents[0] < x: x = coordinates[i][0] + extents[0]
            if coordinates[i][1] + extents[1] < y: y = coordinates[i][1] + extents[1]
            width += extents[4] if i < len(symbol) - 1 else extents[2]
            if extents[3] > height: height = extents[3]
            if i == 0: startWidth = extents[2]
            if i == len(symbol) - 1: endWidth = extents[2]

        if not heavyFirst:
            for i in range(len(coordinates)):
                coordinates[i] = (coordinates[i][0] - (width - startWidth / 2 - endWidth / 2), coordinates[i][1])
            x -= width - startWidth / 2 - endWidth / 2

        # Background
        x1 = x - 2; y1 = y - 2; x2 = x + width + 2; y2 = y + height + 2; r = 4
        cr.move_to(x1 + r, y1)
        cr.line_to(x2 - r, y1)
        cr.curve_to(x2 - r/2, y1, x2, y1 + r/2, x2, y1 + r)
        cr.line_to(x2, y2 - r)
        cr.curve_to(x2, y2 - r/2, x2 - r/2, y2, x2 - r, y2)
        cr.line_to(x1 + r, y2)
        cr.curve_to(x1 + r/2, y2, x1, y2 - r/2, x1, y2 - r)
        cr.line_to(x1, y1 + r)
        cr.curve_to(x1, y1 + r/2, x1 + r/2, y1, x1 + r, y1)
        cr.close_path()
        cr.set_operator(cairo.OPERATOR_CLEAR)
        cr.set_source_rgba(1.0, 1.0, 1.0, 1.0)
        cr.fill()
        cr.set_operator(cairo.OPERATOR_OVER)
        boundingRect = [x1, y1, x2, y2]

        # Set color for text
        if heavyAtom == 'C':    cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        elif heavyAtom == 'N':  cr.set_source_rgba(0.0, 0.0, 1.0, 1.0)
        elif heavyAtom == 'O':  cr.set_source_rgba(1.0, 0.0, 0.0, 1.0)
        elif heavyAtom == 'F':  cr.set_source_rgba(0.5, 0.75, 1.0, 1.0)
        elif heavyAtom == 'Si': cr.set_source_rgba(0.5, 0.5, 0.75, 1.0)
        elif heavyAtom == 'Al': cr.set_source_rgba(0.75, 0.5, 0.5, 1.0)
        elif heavyAtom == 'P':  cr.set_source_rgba(1.0, 0.5, 0.0, 1.0)
        elif heavyAtom == 'S':  cr.set_source_rgba(1.0, 0.75, 0.5, 1.0)
        elif heavyAtom == 'Cl': cr.set_source_rgba(0.0, 1.0, 0.0, 1.0)
        elif heavyAtom == 'Br': cr.set_source_rgba(0.6, 0.2, 0.2, 1.0)
        elif heavyAtom == 'I':  cr.set_source_rgba(0.5, 0.0, 0.5, 1.0)
        else:                   cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)

        # Text itself
        for i, char in enumerate(symbol):
            cr.set_font_size(fontSizeSubscript if char.isdigit() else fontSizeNormal)
            xbearing, ybearing, width, height, xadvance, yadvance = cr.text_extents(char)
            xi, yi = coordinates[i]
            cr.move_to(xi, yi)
            cr.show_text(char)

        x, y = coordinates[0] if heavyFirst else coordinates[-1]
            
    else:
        x = x0; y = y0; width = 0; height = 0
        boundingRect = [x0 - 0.5, y0 - 0.5, x0 + 0.5, y0 + 0.5]
        heavyAtom = ''

    # Draw radical electrons and charges
    # These will be placed either horizontally along the top or bottom of the
    # atom or vertically along the left or right of the atom
    orientation = ' '
    if atom not in bonds or len(bonds[atom]) == 0:
        if len(symbol) == 1:  orientation = 'r'
        else:                 orientation = 'l'
    elif len(bonds[atom]) == 1:
        # Terminal atom - we require a horizontal arrangement if there are
        # more than just the heavy atom
        atom1 = bonds[atom].keys()[0]
        vector = coordinates0[atoms.index(atom),:] - coordinates0[atoms.index(atom1),:]
        if len(symbol) <= 1:
            angle = math.atan2(vector[1], vector[0])
            if 3 * math.pi / 4 <= angle or angle < -3 * math.pi / 4:  orientation = 'l'
            elif -3 * math.pi / 4 <= angle < -1 * math.pi / 4:        orientation = 'b'
            elif -1 * math.pi / 4 <= angle <  1 * math.pi / 4:        orientation = 'r'
            else:                                                     orientation = 't'
        else:
            if vector[1] <= 0:
                orientation = 'b'
            else:
                orientation = 't'
    else:
        # Internal atom
        # First try to see if there is a "preferred" side on which to place the
        # radical/charge data, i.e. if the bonds are unbalanced
        vector = numpy.zeros(2, numpy.float64)
        for atom1 in bonds[atom]:
            vector += coordinates0[atoms.index(atom),:] - coordinates0[atoms.index(atom1),:]
        if numpy.linalg.norm(vector) < 1e-4:
            # All of the bonds are balanced, so we'll need to be more shrewd
            angles = []
            for atom1 in bonds[atom]:
                vector = coordinates0[atoms.index(atom1),:] - coordinates0[atoms.index(atom),:]
                angles.append(math.atan2(vector[1], vector[0]))
            # Try one more time to see if we can use one of the four sides
            # (due to there being no bonds in that quadrant)
            # We don't even need a full 90 degrees open (using 60 degrees instead)
            if   all([ 1 * math.pi / 3 >= angle or angle >=  2 * math.pi / 3 for angle in angles]):  orientation = 't'
            elif all([-2 * math.pi / 3 >= angle or angle >= -1 * math.pi / 3 for angle in angles]):  orientation = 'b'
            elif all([-1 * math.pi / 6 >= angle or angle >=  1 * math.pi / 6 for angle in angles]):  orientation = 'r'
            elif all([ 5 * math.pi / 6 >= angle or angle >= -5 * math.pi / 6 for angle in angles]):  orientation = 'l'
            else:
                # If we still don't have it (e.g. when there are 4+ equally-
                # spaced bonds), just put everything in the top right for now
                orientation = 'tr'
        else:
            # There is an unbalanced side, so let's put the radical/charge data there
            angle = math.atan2(vector[1], vector[0])
            if 3 * math.pi / 4 <= angle or angle < -3 * math.pi / 4:  orientation = 'l'
            elif -3 * math.pi / 4 <= angle < -1 * math.pi / 4:        orientation = 'b'
            elif -1 * math.pi / 4 <= angle <  1 * math.pi / 4:        orientation = 'r'
            else:                                                     orientation = 't'
        
    cr.set_font_size(fontSizeNormal)
    extents = cr.text_extents(heavyAtom)

    # (xi, yi) mark the center of the space in which to place the radicals and charges
    if orientation[0] == 'l':
        xi = x - 2
        yi = y - extents[3]/2
    elif orientation[0] == 'b':
        xi = x + extents[0] + extents[2]/2
        yi = y - extents[3] - 3
    elif orientation[0] == 'r':
        xi = x + extents[0] + extents[2] + 3
        yi = y - extents[3]/2
    elif orientation[0] == 't':
        xi = x + extents[0] + extents[2]/2
        yi = y + 3

    # If we couldn't use one of the four sides, then offset the radical/charges
    # horizontally by a few pixels, in hope that this avoids overlap with an
    # existing bond
    if len(orientation) > 1: xi += 4

    # Get width and height
    cr.set_font_size(fontSizeSubscript)
    width = 0.0; height = 0.0
    if orientation[0] == 'b' or orientation[0] == 't':
        if atom.radicalElectrons > 0:
            width += atom.radicalElectrons * 2 + (atom.radicalElectrons - 1)
            height = atom.radicalElectrons * 2
        text = ''
        if atom.radicalElectrons > 0 and atom.charge != 0: width += 1
        if atom.charge == 1:          text = '+'
        elif atom.charge > 1:         text = '%i+' % atom.charge
        elif atom.charge == -1:       text = u'\u2013'
        elif atom.charge < -1:        text = u'%i\u2013' % abs(atom.charge)
        if text != '':
            extents = cr.text_extents(text)
            width += extents[2] + 1
            height = extents[3]
    elif orientation[0] == 'l' or orientation[0] == 'r':
        if atom.radicalElectrons > 0:
            height += atom.radicalElectrons * 2 + (atom.radicalElectrons - 1)
            width = atom.radicalElectrons * 2
        text = ''
        if atom.radicalElectrons > 0 and atom.charge != 0: height += 1
        if atom.charge == 1:          text = '+'
        elif atom.charge > 1:         text = '%i+' % atom.charge
        elif atom.charge == -1:       text = u'\u2013'
        elif atom.charge < -1:        text = u'%i\u2013' % abs(atom.charge)
        if text != '':
            extents = cr.text_extents(text)
            height += extents[3] + 1
            width = extents[2]
    # Move (xi, yi) to top left corner of space in which to draw radicals and charges
    xi -= width / 2.0; yi -= height / 2.0

    # Update bounding rectangle if necessary
    if width > 0 and height > 0:
        if xi < boundingRect[0]:
            boundingRect[0] = xi
        if yi < boundingRect[1]:
            boundingRect[1] = yi
        if xi + width > boundingRect[2]:
            boundingRect[2] = xi + width
        if yi + height > boundingRect[3]:
            boundingRect[3] = yi + height
        
    if orientation[0] == 'b' or orientation[0] == 't':
        # Draw radical electrons first
        for i in range(atom.radicalElectrons):
            cr.new_sub_path()
            cr.arc(xi + 3 * i + 1, yi + height/2, 1, 0, 2 * math.pi)
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            cr.fill()
        if atom.radicalElectrons > 0: xi += atom.radicalElectrons * 2 + (atom.radicalElectrons - 1) + 1
        # Draw charges second
        text = ''
        if atom.charge == 1:       text = '+'
        elif atom.charge > 1:      text = '%i+' % atom.charge
        elif atom.charge == -1:    text = u'\u2013'
        elif atom.charge < -1:     text = u'%i\u2013' % abs(atom.charge)
        if text != '':
            extents = cr.text_extents(text)
            cr.move_to(xi, yi - extents[1])
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            cr.show_text(text)
    elif orientation[0] == 'l' or orientation[0] == 'r':
        # Draw charges first
        text = ''
        if atom.charge == 1:       text = '+'
        elif atom.charge > 1:      text = '%i+' % atom.charge
        elif atom.charge == -1:    text = u'\u2013'
        elif atom.charge < -1:     text = u'%i\u2013' % abs(atom.charge)
        if text != '':
            extents = cr.text_extents(text)
            cr.move_to(xi - extents[2]/2, yi - extents[1])
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            cr.show_text(text)
        if atom.charge != 0: yi += extents[3] + 1
        # Draw radical electrons second
        for i in range(atom.radicalElectrons):
            cr.new_sub_path()
            cr.arc(xi + width/2, yi + 3 * i + 1, 1, 0, 2 * math.pi)
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            cr.fill()

    return boundingRect

################################################################################

def findLongestPath(chemGraph, atoms0):
    """
    Finds the longest path containing the list of `atoms` in the `chemGraph`.
    The atoms are assumed to already be in a path, with ``atoms[0]`` being a
    terminal atom.
    """
    atom1 = atoms0[-1]
    paths = [atoms0]
    for atom2 in chemGraph.bonds[atom1]:
        if atom2 not in atoms0:
            atoms = atoms0[:]
            atoms.append(atom2)
            paths.append(findLongestPath(chemGraph, atoms))
    lengths = [len(path) for path in paths]
    index = lengths.index(max(lengths))
    return paths[index]

################################################################################

def findBackbone(chemGraph, ringSystems):
    """
    Return the atoms that make up the backbone of the molecule. For acyclic
    molecules, the longest straight chain of heavy atoms will be used. For
    cyclic molecules, the largest independent ring system will be used.
    """

    if chemGraph.isCyclic():
        # Find the largest ring system and use it as the backbone
        # Only count atoms in multiple cycles once
        count = [len(set([atom for ring in ringSystem for atom in ring])) for ringSystem in ringSystems]
        index = 0
        for i in range(1, len(ringSystems)):
            if count[i] > count[index]:
                index = i
        return ringSystems[index]
    
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
    you don't specify this, all atoms in the molecule will be used. The vertices
    are arranged based on a standard bond length of unity, and can be scaled
    later for longer bond lengths. This function ignores any previously-existing
    coordinate information.
    """

    # Initialize array of coordinates
    coordinates = numpy.zeros((len(atoms), 2), numpy.float64)

    # If there are only one or two atoms to draw, then determining the
    # coordinates is trivial
    if len(atoms) == 1:
        coordinates[0,:] = [0.0, 0.0]
        return coordinates
    elif len(atoms) == 2:
        coordinates[0,:] = [0.0, 0.0]
        coordinates[1,:] = [1.0, 0.0]
        return coordinates

    # If the molecule contains cycles, find them and group them
    if chemGraph.isCyclic():
        # This is not a robust method of identifying the ring systems, but will work as a starting point
        cycles = chemGraph.getSmallestSetOfSmallestRings()

        # Split the list of cycles into groups
        # Each atom in the molecule should belong to exactly zero or one such groups
        ringSystems = []
        for cycle in cycles:
            found = False
            for ringSystem in ringSystems:
                for ring in ringSystem:
                    if any([atom in ring for atom in cycle]) and not found:
                        ringSystem.append(cycle)
                        found = True
            if not found:
                ringSystems.append([cycle])
    else:
        ringSystems = []

    # Find the backbone of the molecule
    backbone = findBackbone(chemGraph, ringSystems)

    # Generate coordinates for atoms in backbone
    if chemGraph.isCyclic():
        # Cyclic backbone
        coordinates = generateRingSystemCoordinates(backbone, atoms)
        
        # Flatten backbone so that it contains a list of the atoms in the
        # backbone, rather than a list of the cycles in the backbone
        backbone = list(set([atom for cycle in backbone for atom in cycle]))

    else:
        # Straight chain backbone
        coordinates = generateStraightChainCoordinates(backbone, atoms, bonds)

    # If backbone is linear, then rotate so that the bond is parallel to the
    # horizontal axis
    vector0 = coordinates[atoms.index(backbone[1]),:] - coordinates[atoms.index(backbone[0]),:]
    linear = True
    for i in range(2, len(backbone)):
        vector = coordinates[atoms.index(backbone[i]),:] - coordinates[atoms.index(backbone[i-1]),:]
        if numpy.linalg.norm(vector - vector0) > 1e-4:
            linear = False
            break
    if linear:
        angle = math.atan2(vector0[0], vector0[1]) - math.pi / 2
        rot = numpy.array([[math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]], numpy.float64)
        coordinates = numpy.dot(coordinates, rot)

    # Center backbone at origin
    origin = numpy.zeros(2, numpy.float64)
    for atom in backbone:
        index = atoms.index(atom)
        origin += coordinates[index,:]
    origin /= len(backbone)
    for atom in backbone:
        index = atoms.index(atom)
        coordinates[index,:] -= origin
    
    # We now proceed by calculating the coordinates of the functional groups
    # attached to the backbone
    # Each functional group is independent, although they may contain further
    # branching and cycles
    # In general substituents should try to grow away from the origin to
    # minimize likelihood of overlap
    generateNeighborCoordinates(backbone, atoms, bonds, coordinates, ringSystems)

    return coordinates

################################################################################

def generateStraightChainCoordinates(backbone, atoms, bonds):
    """
    Generate the coordinates for a mutually-adjacent straight chain of atoms
    `backbone`, for which `atoms` and `bonds` are the list and dict of atoms
    and bonds to be rendered, respectively. The general approach is to work from
    one end of the chain to the other, using a horizontal seesaw pattern to lay
    out the coordinates.
    """

    coordinates = numpy.zeros((len(atoms), 2), numpy.float64)

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
            # Rotate by 0 degrees towards horizontal axis (to get angle of 90)
            angle = 0.0
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

    return coordinates

################################################################################

def generateNeighborCoordinates(backbone, atoms, bonds, coordinates, ringSystems):
    """
    Each atom in the backbone must be directly connected to another atom in the
    backbone.
    """

    for i in range(len(backbone)):
        atom0 = backbone[i]
        index0 = atoms.index(atom0)

        # Determine bond angles of all previously-determined bond locations for
        # this atom
        bondAngles = []
        for atom1 in bonds[atom0]:
            index1 = atoms.index(atom1)
            if atom1 in backbone:
                vector = coordinates[index1,:] - coordinates[index0,:]
                angle = math.atan2(vector[1], vector[0])
                bondAngles.append(angle)
        bondAngles.sort()
        
        bestAngle = 2 * math.pi / len(bonds[atom0])
        regular = True
        for angle1, angle2 in zip(bondAngles[0:-1], bondAngles[1:]):
            if all([abs(angle2 - angle1 - (i+1) * bestAngle) > 1e-4 for i in range(len(bonds[atom0]))]):
                regular = False

        if regular:
            # All the bonds around each atom are equally spaced
            # We just need to fill in the missing bond locations

            # Determine rotation angle and matrix
            rot = numpy.array([[math.cos(bestAngle), -math.sin(bestAngle)], [math.sin(bestAngle), math.cos(bestAngle)]], numpy.float64)
            # Determine the vector of any currently-existing bond from this atom
            vector = None
            for atom1 in bonds[atom0]:
                index1 = atoms.index(atom1)
                if atom1 in backbone or numpy.linalg.norm(coordinates[index1,:]) > 1e-4:
                    vector = coordinates[index1,:] - coordinates[index0,:]

            # Iterate through each neighboring atom to this backbone atom
            # If the neighbor is not in the backbone and does not yet have
            # coordinates, then we need to determine coordinates for it
            for atom1 in bonds[atom0]:
                if atom1 not in backbone and numpy.linalg.norm(coordinates[atoms.index(atom1),:]) < 1e-4:
                    occupied = True; count = 0
                    # Rotate vector until we find an unoccupied location
                    while occupied and count < len(bonds[atom0]):
                        count += 1; occupied = False
                        vector = numpy.dot(rot, vector)
                        for atom2 in bonds[atom0]:
                            index2 = atoms.index(atom2)
                            if numpy.linalg.norm(coordinates[index2,:] - coordinates[index0,:] - vector) < 1e-4:
                                occupied = True
                    coordinates[atoms.index(atom1),:] = coordinates[index0,:] + vector
                    generateFunctionalGroupCoordinates(atom0, atom1, atoms, bonds, coordinates, ringSystems)

        else:

            # The bonds are not evenly spaced (e.g. due to a ring)
            # We place all of the remaining bonds evenly over the reflex angle
            startAngle = max(bondAngles)
            endAngle = min(bondAngles)
            if 0.0 < endAngle - startAngle < math.pi: endAngle += 2 * math.pi
            elif 0.0 > endAngle - startAngle > -math.pi: startAngle -= 2 * math.pi
            dAngle = (endAngle - startAngle) / (len(bonds[atom0]) - len(bondAngles) + 1)
            
            index = 1
            for atom1 in bonds[atom0]:
                if atom1 not in backbone and numpy.linalg.norm(coordinates[atoms.index(atom1),:]) < 1e-4:
                    angle = startAngle + index * dAngle
                    index += 1
                    vector = numpy.array([math.cos(angle), math.sin(angle)], numpy.float64)
                    vector /= numpy.linalg.norm(vector)
                    coordinates[atoms.index(atom1),:] = coordinates[index0,:] + vector
                    generateFunctionalGroupCoordinates(atom0, atom1, atoms, bonds, coordinates, ringSystems)

################################################################################

def generateRingSystemCoordinates(ringSystem, atoms):
    """
    Generate the coordinates for all atoms in a mutually-adjacent set of rings
    `ringSystem`, where `atoms` is a list of all atoms to be rendered. The
    general procedure is to (1) find and map the coordinates of the largest
    ring in the system, then (2) iteratively map the coordinates of adjacent
    rings to those already mapped until all rings are processed. This approach
    works well for flat ring systems, but will probably not work when bridge
    atoms are needed.
    """

    coordinates = numpy.zeros((len(atoms), 2), numpy.float64)
    ringSystem = ringSystem[:]
    processed = []

    # Lay out largest cycle in ring system first
    cycle = ringSystem[0]
    for cycle0 in ringSystem[1:]:
        if len(cycle0) > len(cycle):
            cycle = cycle0
    angle = - 2 * math.pi / len(cycle)
    radius = 1.0 / (2 * math.sin(math.pi / len(cycle)))
    for i, atom in enumerate(cycle):
        index = atoms.index(atom)
        coordinates[index,:] = [math.cos(math.pi / 2 + i * angle), math.sin(math.pi / 2 + i * angle)]
        coordinates[index,:] *= radius
    ringSystem.remove(cycle)
    processed.append(cycle)

    # If there are other cycles, then try to lay them out as well
    while len(ringSystem) > 0:

        # Find the largest cycle that shares one or two atoms with a ring that's
        # already been processed
        cycle = None
        for cycle0 in ringSystem:
            for cycle1 in processed:
                count = sum([1 for atom in cycle0 if atom in cycle1])
                if (count == 1 or count == 2):
                    if cycle is None or len(cycle0) > len(cycle): cycle = cycle0
        cycle0 = cycle1
        ringSystem.remove(cycle)

        # Shuffle atoms in cycle such that the common atoms come first
        # Also find the average center of the processed cycles that touch the
        # current cycles
        found = False
        commonAtoms = []
        count = 0
        center0 = numpy.zeros(2, numpy.float64)
        for cycle1 in processed:
            found = False
            for atom in cycle1:
                if atom in cycle and atom not in commonAtoms:
                    commonAtoms.append(atom)
                    found = True
            if found:
                center1 = numpy.zeros(2, numpy.float64)
                for atom in cycle1:
                    center1 += coordinates[atoms.index(atom),:]
                center1 /= len(cycle1)
                center0 += center1
                count += 1
        center0 /= count

        if len(commonAtoms) > 1:
            index0 = cycle.index(commonAtoms[0])
            index1 = cycle.index(commonAtoms[1])
            if (index0 == 0 and index1 == len(cycle) - 1) or (index1 == 0 and index0 == len(cycle) - 1):
                cycle = cycle[-1:] + cycle[0:-1]
            if cycle.index(commonAtoms[1]) < cycle.index(commonAtoms[0]):
                cycle.reverse()
            index = cycle.index(commonAtoms[0])
            cycle = cycle[index:] + cycle[0:index]

        # Determine center of cycle based on already-assigned positions of
        # common atoms (which won't be changed)
        if len(commonAtoms) == 1 or len(commonAtoms) == 2:
            # Center of new cycle is reflection of center of adjacent cycle
            # across common atom or bond
            center = numpy.zeros(2, numpy.float64)
            for atom in commonAtoms:
                center += coordinates[atoms.index(atom),:]
            center /= len(commonAtoms)
            vector = center - center0
            center += vector
            radius = 1.0 / (2 * math.sin(math.pi / len(cycle)))
            
        else:
            # Use any three points to determine the point equidistant from these
            # three; this is the center
            index0 = atoms.index(commonAtoms[0])
            index1 = atoms.index(commonAtoms[1])
            index2 = atoms.index(commonAtoms[2])
            A = numpy.zeros((2,2), numpy.float64)
            b = numpy.zeros((2), numpy.float64)
            A[0,:] = 2 * (coordinates[index1,:] - coordinates[index0,:])
            A[1,:] = 2 * (coordinates[index2,:] - coordinates[index0,:])
            b[0] = coordinates[index1,0]**2 + coordinates[index1,1]**2 - coordinates[index0,0]**2 - coordinates[index0,1]**2
            b[1] = coordinates[index2,0]**2 + coordinates[index2,1]**2 - coordinates[index0,0]**2 - coordinates[index0,1]**2
            center = numpy.linalg.solve(A, b)
            radius = numpy.linalg.norm(center - coordinates[index0,:])
            
        startAngle = 0.0; endAngle = 0.0
        if len(commonAtoms) == 1:
            # We will use the full 360 degrees to place the other atoms in the cycle
            startAngle = math.atan2(-vector[1], vector[0])
            endAngle = startAngle + 2 * math.pi
        elif len(commonAtoms) >= 2:
            # Divide other atoms in cycle equally among unused angle
            vector = coordinates[atoms.index(commonAtoms[-1]),:] - center
            startAngle = math.atan2(vector[1], vector[0])
            vector = coordinates[atoms.index(commonAtoms[0]),:] - center
            endAngle = math.atan2(vector[1], vector[0])
        
        # Place remaining atoms in cycle
        if endAngle < startAngle:
            endAngle += 2 * math.pi
            dAngle = (endAngle - startAngle) / (len(cycle) - len(commonAtoms) + 1)
        else:
            endAngle -= 2 * math.pi
            dAngle = (endAngle - startAngle) / (len(cycle) - len(commonAtoms) + 1)
        
        count = 1
        for i in range(len(commonAtoms), len(cycle)):
            angle = startAngle + count * dAngle
            index = atoms.index(cycle[i])
            # Check that we aren't reassigning any atom positions
            # This version assumes that no atoms belong at the origin, which is
            # usually fine because the first ring is centered at the origin
            if numpy.linalg.norm(coordinates[index,:]) < 1e-4:
                vector = numpy.array([math.cos(angle), math.sin(angle)], numpy.float64)
                coordinates[index,:] = center + radius * vector
            count += 1

        # We're done assigning coordinates for this cycle, so mark it as processed
        processed.append(cycle)
    
    return coordinates

################################################################################

def generateFunctionalGroupCoordinates(atom0, atom1, atoms, bonds, coordinates, ringSystems):
    """
    For the functional group starting with the bond from `atom0` to `atom1`,
    generate the coordinates of the rest of the functional group. `atom0` is
    treated as if a terminal atom. `atom0` and `atom1` must already have their
    coordinates determined. `atoms` is a list of the atoms to be drawn, `bonds`
    is a dictionary of the bonds to draw, and `coordinates` is an array of the
    coordinates for each atom to be drawn. This function is designed to be
    recursive.
    """

    index0 = atoms.index(atom0)
    index1 = atoms.index(atom1)

    # Determine the vector of any currently-existing bond from this atom
    # (We use the bond to the previous atom here)
    vector = coordinates[index0,:] - coordinates[index1,:]

    # Check to see if atom1 is in any cycles in the molecule
    ringSystem = None
    for ringSys in ringSystems:
        if any([atom1 in ring for ring in ringSys]):
            ringSystem = ringSys

    if ringSystem is not None:
        # atom1 is part of a ring system, so we need to process the entire
        # ring system at once

        # Generate coordinates for all atoms in the ring system
        coordinates_cycle = generateRingSystemCoordinates(ringSystem, atoms)

        # Rotate the ring system coordinates so that the line connecting atom1
        # and the center of mass of the ring is parallel to that between
        # atom0 and atom1
        cycleAtoms = list(set([atom for ring in ringSystem for atom in ring]))
        center = numpy.zeros(2, numpy.float64)
        for atom in cycleAtoms:
            center += coordinates_cycle[atoms.index(atom),:]
        center /= len(cycleAtoms)
        vector0 = center - coordinates_cycle[atoms.index(atom1),:]
        angle = math.atan2(vector[1] - vector0[1], vector[0] - vector0[0])
        rot = numpy.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]], numpy.float64)
        coordinates_cycle = numpy.dot(coordinates_cycle, rot)
        
        # Translate the ring system coordinates to the position of atom1
        coordinates_cycle += coordinates[atoms.index(atom1),:] - coordinates_cycle[atoms.index(atom1),:]
        for atom in cycleAtoms:
            coordinates[atoms.index(atom),:] = coordinates_cycle[atoms.index(atom),:]

        # Generate coordinates for remaining neighbors of ring system,
        # continuing to recurse as needed
        generateNeighborCoordinates(cycleAtoms, atoms, bonds, coordinates, ringSystems)
        
    else:
        # atom1 is not in any rings, so we can continue as normal

        # Determine rotation angle and matrix
        numBonds = len(bonds[atom1])
        angle = 0.0
        if numBonds == 2:
            bond0, bond = bonds[atom1].values()
            if (bond0.isTriple() or bond.isTriple()) or (bond0.isDouble() and bond.isDouble()):
                angle = math.pi
            else:
                angle = 2 * math.pi / 3
                # Make sure we're rotating such that we move away from the origin,
                # to discourage overlap of functional groups
                rot1 = numpy.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]], numpy.float64)
                rot2 = numpy.array([[math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]], numpy.float64)
                vector1 = coordinates[index1,:] + numpy.dot(rot1, vector)
                vector2 = coordinates[index1,:] + numpy.dot(rot2, vector)
                if numpy.linalg.norm(vector1) < numpy.linalg.norm(vector2):
                    angle = -angle
        else:
            angle = 2 * math.pi / numBonds
        rot = numpy.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]], numpy.float64)

        # Iterate through each neighboring atom to this backbone atom
        # If the neighbor is not in the backbone, then we need to determine
        # coordinates for it
        for atom, bond in bonds[atom1].iteritems():
            if atom is not atom0:
                occupied = True; count = 0
                # Rotate vector until we find an unoccupied location
                while occupied and count < len(bonds[atom1]):
                    count += 1; occupied = False
                    vector = numpy.dot(rot, vector)
                    for atom2 in bonds[atom1]:
                        index2 = atoms.index(atom2)
                        if numpy.linalg.norm(coordinates[index2,:] - coordinates[index1,:] - vector) < 1e-4:
                            occupied = True
                coordinates[atoms.index(atom),:] = coordinates[index1,:] + vector

                # Recursively continue with functional group
                generateFunctionalGroupCoordinates(atom1, atom, atoms, bonds, coordinates, ringSystems)

################################################################################

def createNewSurface(type, path=None, width=1024, height=768):
    """
    Create a new surface of the specified `type`: "png" for
    :class:`ImageSurface`, "svg" for :class:`SVGSurface`, "pdf" for
    :class:`PDFSurface`, or "ps" for :class:`PSSurface`. If the surface is to
    be saved to a file, use the `path` parameter to give the path to the file.
    You can also optionally specify the `width` and `height` of the generated
    surface if you know what it is; otherwise a default size of 1024 by 768 is
    used.
    """
    import cairo
    type = type.lower()
    if type == 'png':
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, int(width), int(height))
    elif type == 'svg':
        surface = cairo.SVGSurface(path, width, height)
    elif type == 'pdf':
        surface = cairo.PDFSurface(path, width, height)
    elif type == 'ps':
        surface = cairo.PSSurface(path, width, height)
    else:
        raise ValueError('Invalid value "%s" for type parameter; valid values are "png", "svg", "pdf", and "ps".' % type)
    return surface

def drawMolecule(molecule, path=None, surface=''):
    """
    Primary function for generating a drawing of a :class:`Molecule` object
    `molecule`. You can specify the render target in a few ways:

    * If you wish to create an image file (PNG, SVG, PDF, or PS), use the `path`
      parameter to pass a string containing the location at which you wish to
      save the file; the extension will be used to identify the proper target
      type.

    * If you want to render the molecule onto a Cairo surface without saving it
      to a file (e.g. as part of another drawing you are constructing), use the
      `surface` paramter to pass the type of surface you wish to use: "png",
      "svg", "pdf", or "ps".

    This function returns the Cairo surface and context used to create the
    drawing, as well as a bounding box for the molecule being drawn as the
    tuple (`left`, `top`, `width`, `height`).
    """

    try:
        import cairo
    except ImportError:
        print 'Cairo not found; molecule will not be drawn.'
        return

    # This algorithm requires that the hydrogen atoms be implicit
    implicitH = molecule.implicitHydrogens
    molecule.makeHydrogensImplicit()

    atoms = molecule.atoms[:]
    bonds = molecule.bonds.copy()

    # Special cases: H, H2, anything with one heavy atom

    # Remove all unlabeled hydrogen atoms from the molecule, as they are not drawn
    # However, if this would remove all atoms, then don't remove any
    atomsToRemove = []
    for atom in atoms:
        if atom.isHydrogen() and atom.label == '': atomsToRemove.append(atom)
    if len(atomsToRemove) < len(atoms):
        for atom in atomsToRemove:
            atoms.remove(atom)
            for atom2 in bonds[atom]: del bonds[atom2][atom]
            del bonds[atom]

    # Generate the coordinates to use to draw the molecule
    coordinates = generateCoordinates(molecule, atoms, bonds)
    coordinates[:,1] *= -1
    coordinates = coordinates * bondLength

    # Generate labels to use
    symbols = [atom.symbol for atom in atoms]
    for i in range(len(symbols)):
        # Don't label carbon atoms, unless there is only one heavy atom
        if symbols[i] == 'C' and len(symbols) > 1:
            if len(bonds[atoms[i]]) > 1 or (atoms[i].radicalElectrons == 0 and atoms[i].charge == 0):
                symbols[i] = ''
    # Do label atoms that have only double bonds to one or more labeled atoms
    changed = True
    while changed:
        changed = False
        for i in range(len(symbols)):
            if symbols[i] == '' and all([(bond.isDouble() or bond.isTriple()) for bond in bonds[atoms[i]].values()]) and any([symbols[atoms.index(atom)] != '' for atom in bonds[atoms[i]]]):
                symbols[i] = atoms[i].symbol
                changed = True
    # Add implicit hydrogens
    for i in range(len(symbols)):
        if symbols[i] != '':
            if atoms[i].implicitHydrogens == 1: symbols[i] = symbols[i] + 'H'
            elif atoms[i].implicitHydrogens > 1: symbols[i] = symbols[i] + 'H%i' % (atoms[i].implicitHydrogens)

    # Create a dummy surface to draw to, since we don't know the bounding rect
    # We will copy this to another surface with the correct bounding rect
    if path is not None and surface == '':
        type = os.path.splitext(path)[1].lower()[1:]
    else:
        type = surface.lower()
    surface0 = createNewSurface(type=type, path=None)
    cr0 = cairo.Context(surface0)

    # Render using Cairo
    left, top, width, height = render(atoms, bonds, coordinates, symbols, cr0)
    
    # Create the real surface with the appropriate size
    surface = createNewSurface(type=type, path=path, width=width, height=height)
    cr = cairo.Context(surface)
    left, top, width, height = render(atoms, bonds, coordinates, symbols, cr, offset=(-left,-top))

    if path is not None:
        # Finish Cairo drawing
        if surface is not None:
            surface.finish()
        # Save PNG of drawing if appropriate
        ext = os.path.splitext(path)[1].lower()
        if ext == '.png':
            surface.write_to_png(path)

    if not implicitH: molecule.makeHydrogensExplicit()

    return surface, cr, (0, 0, width, height)

################################################################################

if __name__ == '__main__':

    molecule = Molecule()

    # Test #1: Straight chain backbone, no functional groups
    molecule.fromSMILES('C=CC=CCC') # 1,3-hexadiene

    # Test #2: Straight chain backbone, small functional groups
    #molecule.fromSMILES('OCC(O)C(O)C(O)C(O)C(=O)') # glucose

    # Test #3: Straight chain backbone, large functional groups
    #molecule.fromSMILES('CCCCCCCCC(CCCC(CCC)(CCC)CCC)CCCCCCCCC')

    # Test #4: For improved rendering
    # Double bond test #1
    #molecule.fromSMILES('C=CCC=CC(=C)C(=C)C(=O)CC')
    # Double bond test #2
    #molecule.fromSMILES('C=C=O')
    # Radicals
    #molecule.fromSMILES('[O][CH][C]([O])[C]([O])[CH][O]')
    
    # Test #5: Cyclic backbone, no functional groups
    #molecule.fromSMILES('C1=CC=CCC1') # 1,3-cyclohexadiene
    #molecule.fromSMILES('c1ccc2ccccc2c1') # naphthalene
    #molecule.fromSMILES('c1ccc2cc3ccccc3cc2c1') # anthracene
    #molecule.fromSMILES('c1ccc2c(c1)ccc3ccccc32') # phenanthrene
    #molecule.fromSMILES('C1CC2CCCC3C2C1CCC3')

    # Tests #6: Small molecules
    #molecule.fromSMILES('[O]C([O])([O])[O]')

    # Test #7: Cyclic backbone with functional groups
    molecule.fromSMILES('c1ccc(OCc2cc([CH]C)cc2)cc1')

    #molecule.fromSMILES('C=CC(C)(C)CCC')
    #molecule.fromSMILES('CCC(C)CCC(CCC)C')
    #molecule.fromSMILES('C=CC(C)=CCC')
    #molecule.fromSMILES('COC(C)(C)C(C)(C)N(C)C')
    #molecule.fromSMILES('CCC=C=CCCC')
    #molecule.fromSMILES('C1CCCCC1CCC2CCCC2')

    drawMolecule(molecule, 'molecule.pdf')
