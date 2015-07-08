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

"""
This module provides functionality for automatic two-dimensional drawing of the
`skeletal formulae <http://en.wikipedia.org/wiki/Skeletal_formula>`_ of a wide
variety of organic and inorganic molecules. The general method for creating
these drawings is to utilize the :meth:`draw()` method of the :class:`Molecule`
you wish to draw; this wraps a call to :meth:`MoleculeDrawer.draw()`, where the
molecule drawing algorithm begins. Advanced use may require use of the
:class:`MoleculeDrawer` class directly.

The `Cairo <http://cairographics.org/>`_ 2D graphics library is used to create
the drawings. The :class:`MoleculeDrawer` class module will fail gracefully if
Cairo is not installed.

The implementation uses the 2D coordinate generation of rdKit to find coordinates,
then uses Cairo to render the atom.

"""

import math
import numpy
import os.path
import re
import logging

from rmgpy.qm.molecule import Geometry
from rdkit.Chem import AllChem

from numpy.linalg import LinAlgError

################################################################################

def createNewSurface(format, path=None, width=1024, height=768):
    """
    Create a new surface of the specified `type`: "png" for
    :class:`ImageSurface`, "svg" for :class:`SVGSurface`, "pdf" for
    :class:`PDFSurface`, or "ps" for :class:`PSSurface`. If the surface is to
    be saved to a file, use the `path` parameter to give the path to the file.
    You can also optionally specify the `width` and `height` of the generated
    surface if you know what it is; otherwise a default size of 1024 by 768 is
    used.
    """
    try:
        import cairocffi as cairo
    except ImportError:
        import cairo
    format = format.lower()
    if format == 'png':
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, int(width), int(height))
    elif format == 'svg':
        surface = cairo.SVGSurface(path, width, height)
    elif format == 'pdf':
        surface = cairo.PDFSurface(path, width, height)
    elif format == 'ps':
        surface = cairo.PSSurface(path, width, height)
    else:
        raise ValueError('Invalid value "{0}" for type parameter; valid values are "png", "svg", "pdf", and "ps".'.format(type))
    return surface
    
################################################################################

class MoleculeDrawer:
    """
    This class provides functionality for drawing the skeletal formula of
    molecules using the Cairo 2D graphics engine. The most common use case is
    simply::
    
        MoleculeDrawer().draw(molecule, format='png', path='molecule.png')
    
    where ``molecule`` is the :class:`Molecule` object to draw. You can also
    pass a dict of options to the constructor to affect how the molecules are
    drawn.
    """
    
    def __init__(self, options=None):
        self.options = {
            'fontFamily': 'sans',
            'fontSizeNormal': 12,
            'fontSizeSubscript': 8,
            'bondLength': 24,
            'padding': 2,
        }
        if options: self.options.update(options)
        self.clear()
    
    def clear(self):
        self.molecule = None
        self.cycles = None
        self.ringSystems = None
        self.coordinates = None
        self.symbols = None
        self.implicitHydrogens = None
        self.left = 0.0
        self.top = 0.0
        self.right = 0.0
        self.bottom = 0.0
        self.surface = None
        self.cr = None
        
    def draw(self, molecule, format, path=None):
        """
        Draw the given `molecule` using the given image `format` - pdf, svg, ps, or
        png. If `path` is given, the drawing is saved to that location on disk. The
        `options` dict is an optional set of key-value pairs that can be used to
        control the generated drawing.
        
        This function returns the Cairo surface and context used to create the
        drawing, as well as a bounding box for the molecule being drawn as the
        tuple (`left`, `top`, `width`, `height`).
        """
        
        # The Cairo 2D graphics library (and its Python wrapper) is required for
        # the molecule drawing algorithm
        try:
            import cairocffi as cairo
        except ImportError:
            try:
                import cairo
            except ImportError:
                print 'Cairo not found; molecule will not be drawn.'
                return
        
        # Make a copy of the molecule so we don't modify the original
        self.molecule = molecule.copy(deep=True)
        
        # Remove all unlabeled hydrogen atoms from the copied atoms and bonds, as
        # they are not drawn
        # However, if this would remove all atoms, then don't remove any
        atomsToRemove = []
        self.implicitHydrogens = {}
        for atom in self.molecule.atoms:
            if atom.isHydrogen() and atom.label == '': atomsToRemove.append(atom)
        if len(atomsToRemove) < len(self.molecule.atoms):
            for atom in atomsToRemove:
                for atom2 in atom.bonds:
                    try:
                        self.implicitHydrogens[atom2] += 1
                    except KeyError:
                        self.implicitHydrogens[atom2] = 1
                self.molecule.removeAtom(atom)
    
        # Generate information about any cycles present in the molecule, as
        # they will need special attention
        self.__findRingGroups()
        # Handle carbon monoxide special case
        if self.molecule.getFormula() == 'CO' and len(atomsToRemove) == 0:
            # RDKit does not accept atom type Ot
            self.molecule.removeAtom(self.molecule.atoms[-1])
            self.symbols = ['CO']
            self.coordinates = numpy.array([[0,0]], numpy.float64)
        else:
            # Generate the coordinates to use to draw the molecule
            try:
                self.__generateCoordinates()
                
                # Generate labels to use
                self.__generateAtomLabels()
        
            except (ValueError, numpy.linalg.LinAlgError), e:
                logging.error('Error while drawing molecule {0}: {1}'.format(molecule.toSMILES(), e))
                import sys, traceback
                exc_type, exc_value, exc_traceback = sys.exc_info()
                traceback.print_exc()
                return None, None, None

        self.coordinates[:,1] *= -1
        self.coordinates *= self.options['bondLength']
        
        # Handle some special cases
        if self.symbols == ['H','H']:
            # Render as H2 instead of H-H
            self.molecule.removeAtom(self.molecule.atoms[-1])
            self.symbols = ['H2']
            self.coordinates = numpy.array([[0,0]], numpy.float64)
        elif self.symbols == ['O', 'O']:
            # Render as O2 instead of O-O
            self.molecule.removeAtom(self.molecule.atoms[-1])
            self.molecule.atoms[0].radicalElectrons = 0
            self.symbols = ['O2']
            self.coordinates = numpy.array([[0,0]], numpy.float64)
        elif self.symbols == ['OH', 'O'] or self.symbols == ['O', 'OH']:
            # Render as HO2 instead of HO-O or O-OH
            self.molecule.removeAtom(self.molecule.atoms[-1])
            self.symbols = ['O2H']
            self.coordinates = numpy.array([[0,0]], numpy.float64)
        elif self.symbols == ['OH', 'OH']:
            # Render as H2O2 instead of HO-OH or O-OH
            self.molecule.removeAtom(self.molecule.atoms[-1])
            self.symbols = ['O2H2']
            self.coordinates = numpy.array([[0,0]], numpy.float64)
        elif self.symbols == ['O', 'C', 'O']:
            # Render as CO2 instead of O=C=O
            self.molecule.removeAtom(self.molecule.atoms[0])
            self.molecule.removeAtom(self.molecule.atoms[-1])
            self.symbols = ['CO2']
            self.coordinates = numpy.array([[0,0]], numpy.float64)
  
        # Create a dummy surface to draw to, since we don't know the bounding rect
        # We will copy this to another surface with the correct bounding rect
        surface0 = createNewSurface(format=format, path=None)
        cr0 = cairo.Context(surface0)
    
        # Render using Cairo
        self.render(cr0)
        
        # Create the real surface with the appropriate size
        xoff = self.left
        yoff = self.top
        width = self.right - self.left
        height = self.bottom - self.top
        self.surface = createNewSurface(format=format, path=path, width=width, height=height)
        self.cr = cairo.Context(self.surface)

        # Draw white background
        self.cr.set_source_rgba(1.0, 1.0, 1.0, 1.0)
        self.cr.paint()        

        self.render(self.cr, offset=(-xoff,-yoff))

        if path is not None:
            # Finish Cairo drawing
            # Save PNG of drawing if appropriate
            ext = os.path.splitext(path)[1].lower()
            if ext == '.png':
                self.surface.write_to_png(path)
            else:
                self.surface.finish()
    
        return self.surface, self.cr, (xoff, yoff, width, height)

    def __findRingGroups(self):
        """
        Find all of the cycles in the current molecule, and group them into
        sets of adjacent cycles.
        """
                
        # Find all of the cycles in the molecule
        self.cycles = self.molecule.getSmallestSetOfSmallestRings()
        self.ringSystems = []
        
        # If the molecule contains cycles, find them and group them
        if len(self.cycles) > 0:
            # Split the list of cycles into groups
            # Each atom in the molecule should belong to exactly zero or one such groups
            for cycle in self.cycles:
                found = False
                for ringSystem in self.ringSystems:
                    for ring in ringSystem:
                        if any([atom in ring for atom in cycle]) and not found:
                            ringSystem.append(cycle)
                            found = True
                if not found:
                    self.ringSystems.append([cycle])
    
    def __generateCoordinates(self):
        """
        Generate the 2D coordinates to be used when drawing the current 
        molecule. The function uses rdKits 2D coordinate generation.
        """
        atoms = self.molecule.atoms
        Natoms = len(atoms)
        flag_charge = 0
        
        for atom in self.molecule.atoms:
            if atom.charge != 0: #atomType.label in ['N5s','N5d','N5dd','N5t','N5b']:
                 flag_charge = 1
                 break
        
        # Initialize array of coordinates
        self.coordinates = coordinates = numpy.zeros((Natoms, 2))
        
        if flag_charge == 1:
            # If there are only one or two atoms to draw, then determining the
            # coordinates is trivial
            if Natoms == 1:
                self.coordinates[0,:] = [0.0, 0.0]
                return self.coordinates
            elif Natoms == 2:
                self.coordinates[0,:] = [-0.5, 0.0]
                self.coordinates[1,:] = [0.5, 0.0]
                return self.coordinates
        
            if len(self.cycles) > 0:
                # Cyclic molecule
                backbone = self.__findCyclicBackbone()
                self.__generateRingSystemCoordinates(backbone)
                # Flatten backbone so that it contains a list of the atoms in the
                # backbone, rather than a list of the cycles in the backbone
                backbone = list(set([atom for cycle in backbone for atom in cycle]))
            else:
                # Straight chain molecule
                backbone = self.__findStraightChainBackbone()
                self.__generateStraightChainCoordinates(backbone)
                
                # If backbone is linear, then rotate so that the bond is parallel to the
                # horizontal axis
                vector0 = coordinates[atoms.index(backbone[1]),:] - coordinates[atoms.index(backbone[0]),:]
                for i in range(2, len(backbone)):
                    vector = coordinates[atoms.index(backbone[i]),:] - coordinates[atoms.index(backbone[i-1]),:]
                    if numpy.linalg.norm(vector - vector0) > 1e-4:
                        break
                else:
                    angle = math.atan2(vector0[0], vector0[1]) - math.pi / 2
                    rot = numpy.array([[math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]], numpy.float64)
                    coordinates = numpy.dot(coordinates, rot)
                
            # Center backbone at origin
            xmin = numpy.min(coordinates[:,0])
            xmax = numpy.max(coordinates[:,0])
            ymin = numpy.min(coordinates[:,1])
            ymax = numpy.max(coordinates[:,1])
            xmid = 0.5 * (xmax + xmin)
            ymid = 0.5 * (ymax + ymin)
            for atom in backbone:
                index = atoms.index(atom)
                coordinates[index,0] -= xmid
                coordinates[index,1] -= ymid
            
            # We now proceed by calculating the coordinates of the functional groups
            # attached to the backbone
            # Each functional group is independent, although they may contain further
            # branching and cycles
            # In general substituents should try to grow away from the origin to
            # minimize likelihood of overlap
            self.__generateNeighborCoordinates(backbone)
            
            return coordinates
            
        else:
            
            # Use rdkit 2D coordinate generation:
            
            # Generate the RDkit molecule from the RDkit molecule, use geometry
            # in order to match the atoms in the rdmol with the atoms in the
            # RMG molecule (which is required to extract coordinates).
            self.geometry = Geometry(None, None, self.molecule, None)
            
            rdmol, rdAtomIdx = self.geometry.rd_build()
            AllChem.Compute2DCoords(rdmol)
            
            # Extract the coordinates from each atom.
            for atom in atoms:
                index = rdAtomIdx[atom]
                point = rdmol.GetConformer(0).GetAtomPosition(index)
                coordinates[index,:]= [point.x*0.6, point.y*0.6]
            
            # RDKit generates some molecules more vertically than horizontally,
            # Especially linear ones. This will reflect any molecule taller than
            # it is wide across the line y=x
            ranges = numpy.ptp(coordinates, axis = 0)
            if ranges[1] > ranges[0]:
                temp = numpy.copy(coordinates)
                coordinates[:,0] = temp[:,1]
                coordinates[:,1] = temp[:,0]
            
            return coordinates
    
    def __findCyclicBackbone(self):
        """
        Return a set of atoms to use as the "backbone" of the molecule. For
        cyclics this is simply the largest ring system.
        """
        count = [len(set([atom for ring in ringSystem for atom in ring])) for ringSystem in self.ringSystems]
        index = 0
        for i in range(1, len(self.ringSystems)):
            if count[i] > count[index]:
                index = i
        return self.ringSystems[index]
    
    def __findStraightChainBackbone(self):
        """
        Return a set of atoms to use as the "backbone" of the molecule. For
        non-cyclics this is the largest straight chain between atoms. If carbon
        atoms are present, then we define the backbone only in terms of them.
        """
        # Find the terminal atoms - those that only have one explicit bond
        terminalAtoms = [atom for atom in self.molecule.atoms if len(atom.bonds) == 1]
        assert len(terminalAtoms) >= 2
        
        # Starting from each terminal atom, find the longest straight path to
        # another terminal
        # The longest found is the backbone
        backbone = []
        paths = []
        for atom in terminalAtoms:
            paths.extend(self.__findStraightChainPaths([atom]))
        
        # Remove any paths that don't end in a terminal atom
        # (I don't think this should remove any!)
        paths = [path for path in paths if path[-1] in terminalAtoms]
        
        # Remove all paths shorter than the maximum
        length = max([len(path) for path in paths])
        paths = [path for path in paths if len(path) == length]
        
        # Prefer the paths with the most carbon atoms
        carbons = [sum([1 for atom in path if atom.isCarbon()]) for path in paths]
        maxCarbons = max(carbons)
        paths = [path for path, carbon in zip(paths, carbons) if carbon == maxCarbons]
        
        # At this point we could choose any remaining path, so simply choose the first
        backbone = paths[0]

        assert len(backbone) > 1
        assert backbone[0] in terminalAtoms
        assert backbone[-1] in terminalAtoms
        
        return backbone
    
    def __findStraightChainPaths(self, atoms0):
        """
        Finds the paths containing the list of atoms `atoms0` in the
        current molecule. The atoms are assumed to already be in a path, with
         ``atoms0[0]`` being a terminal atom.
        """
        atom1 = atoms0[-1]
        paths = []
        for atom2 in atom1.bonds:
            if atom2 not in atoms0:
                atoms = atoms0[:]
                atoms.append(atom2)
                if not self.molecule.isAtomInCycle(atom2):
                    paths.extend(self.__findStraightChainPaths(atoms))
        if len(paths) == 0:
            paths.append(atoms0[:])
        return paths

    def __generateRingSystemCoordinates(self, atoms):
        """
        For a ring system composed of the given cycles of `atoms`, update the
        coordinates of each atom in the system.
        """
        coordinates = self.coordinates
        atoms = atoms[:]
        processed = []
    
        # Lay out largest cycle in ring system first
        cycle = atoms[0]
        for cycle0 in atoms[1:]:
            if len(cycle0) > len(cycle):
                cycle = cycle0
        angle = - 2 * math.pi / len(cycle)
        radius = 1.0 / (2 * math.sin(math.pi / len(cycle)))
        for i, atom in enumerate(cycle):
            index = self.molecule.atoms.index(atom)
            coordinates[index,:] = [math.cos(math.pi / 2 + i * angle), math.sin(math.pi / 2 + i * angle)]
            coordinates[index,:] *= radius
        atoms.remove(cycle)
        processed.append(cycle)
    
        # If there are other cycles, then try to lay them out as well
        while len(atoms) > 0:
    
            # Find the largest cycle that shares one or two atoms with a ring that's
            # already been processed
            cycle = None
            for cycle0 in atoms:
                for cycle1 in processed:
                    count = sum([1 for atom in cycle0 if atom in cycle1])
                    if (count == 1 or count == 2):
                        if cycle is None or len(cycle0) > len(cycle): cycle = cycle0
            cycle0 = cycle1
            atoms.remove(cycle)
    
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
                        center1 += coordinates[cycle1.index(atom),:]
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
                    center += coordinates[self.molecule.atoms.index(atom),:]
                center /= len(commonAtoms)
                vector = center - center0
                center += vector
                radius = 1.0 / (2 * math.sin(math.pi / len(cycle)))
                
            else:
                # Use any three points to determine the point equidistant from these
                # three; this is the center
                index0 = self.molecule.atoms.index(commonAtoms[0])
                index1 = self.molecule.atoms.index(commonAtoms[1])
                index2 = self.molecule.atoms.index(commonAtoms[2])
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
                vector = coordinates[cycle.index(commonAtoms[-1]),:] - center
                startAngle = math.atan2(vector[1], vector[0])
                vector = coordinates[cycle.index(commonAtoms[0]),:] - center
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
                index = self.molecule.atoms.index(cycle[i])
                # Check that we aren't reassigning any atom positions
                # This version assumes that no atoms belong at the origin, which is
                # usually fine because the first ring is centered at the origin
                if numpy.linalg.norm(coordinates[index,:]) < 1e-4:
                    vector = numpy.array([math.cos(angle), math.sin(angle)], numpy.float64)
                    coordinates[index,:] = center + radius * vector
                count += 1
    
            # We're done assigning coordinates for this cycle, so mark it as processed
            processed.append(cycle)
    
    def __generateStraightChainCoordinates(self, atoms):
        """
        Update the coordinates for the linear straight chain of `atoms` in
        the current molecule. 
        """
        coordinates = self.coordinates
    
        # First atom goes at origin
        index0 = self.molecule.atoms.index(atoms[0])
        coordinates[index0,:] = [0.0, 0.0]
    
        # Second atom goes on x-axis (for now; this could be improved!)
        index1 = self.molecule.atoms.index(atoms[1])
        vector = numpy.array([1.0, 0.0], numpy.float64)
        if atoms[0].bonds[atoms[1]].isTriple():
            rotatePositive = False
        else:
            rotatePositive = True
            rot = numpy.array([[math.cos(-math.pi / 6), math.sin(-math.pi / 6)], [-math.sin(-math.pi / 6), math.cos(-math.pi / 6)]], numpy.float64)
            vector = numpy.array([1.0, 0.0], numpy.float64)
            vector = numpy.dot(rot, vector)
        coordinates[index1,:] = coordinates[index0,:] + vector
    
        # Other atoms
        for i in range(2, len(atoms)):
            atom0 = atoms[i-2]
            atom1 = atoms[i-1]
            atom2 = atoms[i]
            index1 = self.molecule.atoms.index(atom1)
            index2 = self.molecule.atoms.index(atom2)
            bond0 = atom0.bonds[atom1]
            bond = atom1.bonds[atom2]
            # Angle of next bond depends on the number of bonds to the start atom
            numBonds = len(atom1.bonds)
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
    
    def __generateNeighborCoordinates(self, backbone):
        """
        Recursively update the coordinates for the atoms immediately adjacent
        to the atoms in the molecular `backbone`.
        """
        atoms = self.molecule.atoms
        coordinates = self.coordinates
        
        for i in range(len(backbone)):
            atom0 = backbone[i]
            index0 = atoms.index(atom0)
    
            # Determine bond angles of all previously-determined bond locations for
            # this atom
            bondAngles = []
            for atom1 in atom0.bonds:
                index1 = atoms.index(atom1)
                if atom1 in backbone:
                    vector = coordinates[index1,:] - coordinates[index0,:]
                    angle = math.atan2(vector[1], vector[0])
                    bondAngles.append(angle)
            bondAngles.sort()
            
            bestAngle = 2 * math.pi / len(atom0.bonds)
            regular = True
            for angle1, angle2 in zip(bondAngles[0:-1], bondAngles[1:]):
                if all([abs(angle2 - angle1 - (i+1) * bestAngle) > 1e-4 for i in range(len(atom0.bonds))]):
                    regular = False
    
            if regular:
                # All the bonds around each atom are equally spaced
                # We just need to fill in the missing bond locations
    
                # Determine rotation angle and matrix
                rot = numpy.array([[math.cos(bestAngle), -math.sin(bestAngle)], [math.sin(bestAngle), math.cos(bestAngle)]], numpy.float64)
                # Determine the vector of any currently-existing bond from this atom
                vector = None
                for atom1 in atom0.bonds:
                    index1 = atoms.index(atom1)
                    if atom1 in backbone or numpy.linalg.norm(coordinates[index1,:]) > 1e-4:
                        vector = coordinates[index1,:] - coordinates[index0,:]
    
                # Iterate through each neighboring atom to this backbone atom
                # If the neighbor is not in the backbone and does not yet have
                # coordinates, then we need to determine coordinates for it
                for atom1 in atom0.bonds:
                    if atom1 not in backbone and numpy.linalg.norm(coordinates[atoms.index(atom1),:]) < 1e-4:
                        occupied = True; count = 0
                        # Rotate vector until we find an unoccupied location
                        while occupied and count < len(atom0.bonds):
                            count += 1; occupied = False
                            vector = numpy.dot(rot, vector)
                            for atom2 in atom0.bonds:
                                index2 = atoms.index(atom2)
                                if numpy.linalg.norm(coordinates[index2,:] - coordinates[index0,:] - vector) < 1e-4:
                                    occupied = True
                        coordinates[atoms.index(atom1),:] = coordinates[index0,:] + vector
                        self.__generateFunctionalGroupCoordinates(atom0, atom1)
    
            else:
    
                # The bonds are not evenly spaced (e.g. due to a ring)
                # We place all of the remaining bonds evenly over the reflex angle
                startAngle = max(bondAngles)
                endAngle = min(bondAngles)
                if 0.0 < endAngle - startAngle < math.pi: endAngle += 2 * math.pi
                elif 0.0 > endAngle - startAngle > -math.pi: startAngle -= 2 * math.pi
                dAngle = (endAngle - startAngle) / (len(atom0.bonds) - len(bondAngles) + 1)
                
                index = 1
                for atom1 in atom0.bonds:
                    if atom1 not in backbone and numpy.linalg.norm(coordinates[atoms.index(atom1),:]) < 1e-4:
                        angle = startAngle + index * dAngle
                        index += 1
                        vector = numpy.array([math.cos(angle), math.sin(angle)], numpy.float64)
                        vector /= numpy.linalg.norm(vector)
                        coordinates[atoms.index(atom1),:] = coordinates[index0,:] + vector
                        self.__generateFunctionalGroupCoordinates(atom0, atom1)

    def __generateFunctionalGroupCoordinates(self, atom0, atom1):
        """
        For the functional group starting with the bond from `atom0` to `atom1`,
        generate the coordinates of the rest of the functional group. `atom0` is
        treated as if a terminal atom. `atom0` and `atom1` must already have their
        coordinates determined. `atoms` is a list of the atoms to be drawn, `bonds`
        is a dictionary of the bonds to draw, and `coordinates` is an array of the
        coordinates for each atom to be drawn. This function is designed to be
        recursive.
        """
    
        atoms = self.molecule.atoms
        coordinates = self.coordinates

        index0 = atoms.index(atom0)
        index1 = atoms.index(atom1)
    
        # Determine the vector of any currently-existing bond from this atom
        # (We use the bond to the previous atom here)
        vector = coordinates[index0,:] - coordinates[index1,:]
        bondAngle = math.atan2(vector[1], vector[0])
        
        # Check to see if atom1 is in any cycles in the molecule
        ringSystem = None
        for ringSys in self.ringSystems:
            if any([atom1 in ring for ring in ringSys]):
                ringSystem = ringSys
    
        if ringSystem is not None:
            # atom1 is part of a ring system, so we need to process the entire
            # ring system at once
    
            # Generate coordinates for all atoms in the ring system
            self.__generateRingSystemCoordinates(ringSystem)
    
            cycleAtoms = list(set([atom for ring in ringSystem for atom in ring]))
            
            coordinates_cycle = numpy.zeros_like(self.coordinates)
            for atom in cycleAtoms:
                coordinates_cycle[atoms.index(atom),:] = coordinates[atoms.index(atom),:]
    
            # Rotate the ring system coordinates so that the line connecting atom1
            # and the center of mass of the ring is parallel to that between
            # atom0 and atom1
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
            self.__generateNeighborCoordinates(cycleAtoms)
            
        else:
            # atom1 is not in any rings, so we can continue as normal
    
            # Determine rotation angle and matrix
            numBonds = len(atom1.bonds)
            angle = 0.0
            if numBonds == 2:
                bond0, bond = atom1.bonds.values()
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
                    if bondAngle < -0.5 * math.pi or bondAngle > 0.5 * math.pi:
                        angle = abs(angle)
                    else:
                        angle = -abs(angle)
            else:
                angle = 2 * math.pi / numBonds
            rot = numpy.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]], numpy.float64)
    
            # Iterate through each neighboring atom to this backbone atom
            # If the neighbor is not in the backbone, then we need to determine
            # coordinates for it
            for atom, bond in atom1.bonds.iteritems():
                if atom is not atom0:
                    occupied = True; count = 0
                    # Rotate vector until we find an unoccupied location
                    while occupied and count < len(atom1.bonds):
                        count += 1; occupied = False
                        vector = numpy.dot(rot, vector)
                        for atom2 in atom1.bonds:
                            index2 = atoms.index(atom2)
                            if numpy.linalg.norm(coordinates[index2,:] - coordinates[index1,:] - vector) < 1e-4:
                                occupied = True
                    coordinates[atoms.index(atom),:] = coordinates[index1,:] + vector
    
                    # Recursively continue with functional group
                    self.__generateFunctionalGroupCoordinates(atom1, atom)
    
    def __generateAtomLabels(self):
        """
        Generate the labels to use for each atom in the drawing. In general,
        all atoms are labeled with their symbols except carbon. Some carbon
        atoms are also labeled in certain circumstances. The labels also
        contain any implicit hydrogen atoms (i.e. those hydrogen atoms not
        explicitly drawn in the skeletal formula).
        """
        atoms = self.molecule.atoms
        
        self.symbols = symbols = [atom.symbol for atom in atoms]
        for i in range(len(symbols)):
            # Don't label carbon atoms, unless there are only one or two heavy atoms
            if symbols[i] == 'C' and len(symbols) > 2:
                if len(atoms[i].bonds) > 1 or (atoms[i].radicalElectrons == 0 and atoms[i].charge == 0):
                    symbols[i] = ''
        # Do label atoms that have only double bonds to one or more labeled atoms
        changed = True
        while changed:
            changed = False
            for i in range(len(symbols)):
                if symbols[i] == '' and all([(bond.isDouble() or bond.isTriple()) for bond in atoms[i].bonds.values()]) and any([symbols[atoms.index(atom)] != '' for atom in atoms[i].bonds]):
                    symbols[i] = atoms[i].symbol
                    changed = True
        # Add implicit hydrogens
        for i in range(len(symbols)):
            if symbols[i] != '':
                try:
                    Hcount = self.implicitHydrogens[atoms[i]]
                except KeyError:
                    continue
                if Hcount == 1: symbols[i] = symbols[i] + 'H'
                elif Hcount > 1: symbols[i] = symbols[i] + 'H{0:d}'.format(Hcount)
        
        return symbols

    def render(self, cr, offset=None):
        """
        Uses the Cairo graphics library to create a skeletal formula drawing of a
        molecule containing the list of `atoms` and dict of `bonds` to be drawn.
        The 2D position of each atom in `atoms` is given in the `coordinates` array.
        The symbols to use at each atomic position are given by the list `symbols`.
        You must specify the Cairo context `cr` to render to.
        """
    
        try:
            import cairocffi as cairo
        except ImportError:
            import cairo
        
        coordinates = self.coordinates
        atoms = self.molecule.atoms
        symbols = self.symbols
        
        drawLonePairs = False
        
        for atom in atoms:
            if atom.isNitrogen():
                drawLonePairs = True
    
        left = 0.0
        top = 0.0
        right = 0.0
        bottom = 0.0
        
        # Shift coordinates by offset value
        if offset is not None:
            coordinates[:,0] += offset[0]
            coordinates[:,1] += offset[1]
        
        # Draw bonds
        for atom1 in atoms:
            for atom2, bond in atom1.bonds.items():
                index1 = atoms.index(atom1)
                index2 = atoms.index(atom2)
                if index1 < index2: # So we only draw each bond once
                    self.__renderBond(index1, index2, bond, cr)
    
        # Draw aromatic bonds
        for cycle in self.cycles:
            cycleBonds = []
            for atom1, atom2 in zip(cycle[0:-1], cycle[1:]):
                cycleBonds.append(atom1.bonds[atom2])
            cycleBonds.append(cycle[0].bonds[cycle[-1]])
            if all([bond.isBenzene() for bond in cycleBonds]):
                # We've found an aromatic ring, so draw a circle in the center to represent the benzene bonds
                center = numpy.zeros(2, numpy.float64)
                for atom in cycle:
                    index = atoms.index(atom)
                    center += coordinates[index,:]
                center /= len(cycle)
                index1 = atoms.index(cycle[0])
                index2 = atoms.index(cycle[1])
                radius = math.sqrt(
                    (center[0] - (coordinates[index1,0] + coordinates[index2,0]) / 2)**2 +
                    (center[1] - (coordinates[index1,1] + coordinates[index2,1]) / 2)**2
                ) - 4
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.set_line_width(1.0)
                cr.set_line_cap(cairo.LINE_CAP_ROUND)
                cr.arc(center[0], center[1], radius, 0.0, 2 * math.pi)
                cr.stroke()
    
        # Draw atoms
        for i, atom in enumerate(atoms):
            symbol = symbols[i]
            index = atoms.index(atom)
            x0, y0 = coordinates[index,:]
            vector = numpy.zeros(2, numpy.float64)
            for atom2 in atom.bonds:
                vector += coordinates[atoms.index(atom2),:] - coordinates[index,:]
            heavyFirst = vector[0] <= 0
            if len(atoms) == 1 and atoms[0].symbol not in ['C', 'N'] and atoms[0].charge == 0 and atoms[0].radicalElectrons == 0:
                # This is so e.g. water is rendered as H2O rather than OH2
                heavyFirst = False
                cr.set_font_size(self.options['fontSizeNormal'])
                x0 += cr.text_extents(symbols[0])[2] / 2.0
            atomBoundingRect = self.__renderAtom(symbol, atom, x0, y0, cr, heavyFirst, drawLonePairs)
        
        # Add a small amount of whitespace on all sides
        padding = self.options['padding']
        self.left -= padding; self.top -= padding; self.right += padding; self.bottom += padding
    
    def __drawLine(self, cr, x1, y1, x2, y2):
        """
        Draw a line on the given Cairo context `cr` from (`x1`, `y1`) to
        (`x2`,`y2`), and update the bounding rectangle if necessary.
        """
        try:
            import cairocffi as cairo
        except ImportError:
            import cairo
        cairo
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.set_line_width(1.0)
        cr.set_line_cap(cairo.LINE_CAP_ROUND)
        cr.move_to(x1, y1); cr.line_to(x2, y2)
        cr.stroke()
        if x1 < self.left: self.left = x1
        if x1 > self.right: self.right = x1
        if y1 < self.top: self.top = y1
        if y1 > self.bottom: self.bottom = y1
        if x2 < self.left: self.left = x2
        if x2 > self.right: self.right = x2
        if y2 < self.top: self.top = y2
        if y2 > self.bottom: self.bottom = y2
    
    def __renderBond(self, atom1, atom2, bond, cr):
        """
        Render an individual `bond` between atoms with indices `atom1` and `atom2`
        on the Cairo context `cr`.
        """
    
        try:
            import cairocffi as cairo
        except ImportError:
            import cairo
    
        bondLength = self.options['bondLength']
    
        x1, y1 = self.coordinates[atom1,:]
        x2, y2 = self.coordinates[atom2,:]
        angle = math.atan2(y2 - y1, x2 - x1)
    
        dx = x2 - x1; dy = y2 - y1
        du = math.cos(angle + math.pi / 2)
        dv = math.sin(angle + math.pi / 2)
        if bond.isDouble() and (self.symbols[atom1] != '' or self.symbols[atom2] != ''):
            # Draw double bond centered on bond axis
            du *= 1.6; dv *= 1.6
            self.__drawLine(cr, x1 - du, y1 - dv, x2 - du, y2 - dv)
            self.__drawLine(cr, x1 + du, y1 + dv, x2 + du, y2 + dv)
        elif bond.isTriple() and (self.symbols[atom1] != '' or self.symbols[atom2] != ''):
            # Draw triple bond centered on bond axis
            du *= 3; dv *= 3
            self.__drawLine(cr, x1 - du, y1 - dv, x2 - du, y2 - dv)
            self.__drawLine(cr, x1     , y1     , x2     , y2     )
            self.__drawLine(cr, x1 + du, y1 + dv, x2 + du, y2 + dv)
        else:
            # Draw bond on skeleton
            self.__drawLine(cr, x1, y1, x2, y2)
            # Draw other bonds
            if bond.isDouble():
                du *= 3.2; dv *= 3.2; dx = 2 * dx / bondLength; dy = 2 * dy / bondLength
                self.__drawLine(cr, x1 + du + dx, y1 + dv + dy, x2 + du - dx, y2 + dv - dy)
            elif bond.isTriple():
                du *= 3; dv *= 3; dx = 2 * dx / bondLength; dy = 2 * dy / bondLength
                self.__drawLine(cr, x1 - du + dx, y1 - dv + dy, x2 - du - dx, y2 - dv - dy)
                self.__drawLine(cr, x1 + du + dx, y1 + dv + dy, x2 + du - dx, y2 + dv - dy)
        
    def __renderAtom(self, symbol, atom, x0, y0, cr, heavyFirst=True, drawLonePairs=False):
        """
        Render the `label` for an atom centered around the coordinates (`x0`, `y0`)
        onto the Cairo context `cr`. If `heavyFirst` is ``False``, then the order
        of the atoms will be reversed in the symbol. This method also causes
        radical electrons and charges to be drawn adjacent to the rendered symbol.
        """
    
        try:
            import cairocffi as cairo
        except ImportError:
            import cairo
    
        atoms = self.molecule.atoms
    
        if symbol != '':
            heavyAtom = symbol[0]
    
            # Split label by atoms
            labels = re.findall('[A-Z][a-z]*[0-9]*', symbol)
            if not heavyFirst: labels.reverse()
            if 'C' not in symbol and 'O' not in symbol and len(atoms) == 1: labels.sort()
            symbol = ''.join(labels)
    
            # Determine positions of each character in the symbol
            coordinates = []
    
            cr.set_font_size(self.options['fontSizeNormal'])
            y0 += max([cr.text_extents(char)[3] for char in symbol if char.isalpha()]) / 2
    
            for i, label in enumerate(labels):
                for j, char in enumerate(label):
                    cr.set_font_size(self.options['fontSizeSubscript' if char.isdigit() else 'fontSizeNormal'])
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
                cr.set_font_size(self.options['fontSizeSubscript' if char.isdigit() else 'fontSizeNormal'])
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
            
            cr.save()
            cr.set_operator(cairo.OPERATOR_SOURCE)
            cr.set_source_rgba(1.0, 1.0, 1.0, 1.0)
            cr.fill()
            cr.restore()
            
            boundingRect = [x1, y1, x2, y2]
    
            # Set color for text
            if   heavyAtom == 'C':  cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
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
                cr.set_font_size(self.options['fontSizeSubscript' if char.isdigit() else 'fontSizeNormal'])
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
        if len(atom.bonds) == 0:
            if len(symbol) == 1:  orientation = 'r'
            else:                 orientation = 'l'
        elif len(atom.bonds) == 1:
            # Terminal atom - we require a horizontal arrangement if there are
            # more than just the heavy atom
            atom1 = atom.bonds.keys()[0]
            vector = self.coordinates[atoms.index(atom),:] - self.coordinates[atoms.index(atom1),:]
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
            for atom1 in atom.bonds:
                vector += self.coordinates[atoms.index(atom),:] - self.coordinates[atoms.index(atom1),:]
            if numpy.linalg.norm(vector) < 1e-4:
                # All of the bonds are balanced, so we'll need to be more shrewd
                angles = []
                for atom1 in atom.bonds:
                    vector = self.coordinates[atoms.index(atom1),:] - self.coordinates[atoms.index(atom),:]
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
            
        cr.set_font_size(self.options['fontSizeNormal'])
        extents = cr.text_extents(heavyAtom)
    
        # (xi, yi) mark the center of the space in which to place the radicals and charges
        if orientation[0] == 'l':
            xi = x - 3
            yi = y - extents[3]/2
        elif orientation[0] == 'b':
            xi = x + extents[0] + extents[2]/2
            yi = y - extents[3] - 4
        elif orientation[0] == 'r':
            xi = x + extents[0] + extents[2] + 4
            yi = y - extents[3]/2
        elif orientation[0] == 't':
            xi = x + extents[0] + extents[2]/2
            yi = y + 4
    
        # If we couldn't use one of the four sides, then offset the radical/charges
        # horizontally by a few pixels, in hope that this avoids overlap with an
        # existing bond
        if len(orientation) > 1: xi += 4
    
        # Get width and height
        cr.set_font_size(self.options['fontSizeSubscript'])
        width = 0.0; height = 0.0
        if orientation[0] == 'b' or orientation[0] == 't':
            if atom.radicalElectrons > 0:
                width += atom.radicalElectrons * 2 + (atom.radicalElectrons - 1)
                height = atom.radicalElectrons * 2
            text = ''
            if atom.radicalElectrons > 0 and atom.charge != 0: width += 1
            if atom.charge == 1:          text = '+'
            elif atom.charge > 1:         text = '{0:d}+'.format(atom.charge)
            elif atom.charge == -1:       text = u'\u2013'
            elif atom.charge < -1:        text = u'{0:d}\u2013'.format(abs(atom.charge))
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
            elif atom.charge > 1:         text = '{0:d}+'.format(atom.charge)
            elif atom.charge == -1:       text = u'\u2013'
            elif atom.charge < -1:        text = u'{0:d}\u2013'.format(abs(atom.charge))
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
            elif atom.charge > 1:      text = '{0:d}+'.format(atom.charge)
            elif atom.charge == -1:    text = u'\u2013'
            elif atom.charge < -1:     text = u'{0:d}\u2013'.format(abs(atom.charge))
            if text != '':
                extents = cr.text_extents(text)
                cr.move_to(xi, yi - extents[1])
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.show_text(text)
                
            # Draw lone electron pairs            
            # Draw them for nitrogen containing molecules only
            if drawLonePairs:
                for i in range(atom.lonePairs):
                    cr.new_sub_path()
                    if i == 0:
                        x1lp = x-2
                        y1lp = y-8
                        x2lp = x+2
                        y2lp = y-12
                    elif i == 1:
                        x1lp = x+12
                        y1lp = y-8
                        x2lp = x+8
                        y2lp = y-12
                    elif i == 2:
                        x1lp = x-2
                        y1lp = y-1
                        x2lp = x+2
                        y2lp = y+3
                    self.__drawLine(cr, x1lp, y1lp, x2lp, y2lp)
                
        elif orientation[0] == 'l' or orientation[0] == 'r':
            # Draw charges first
            text = ''
            if atom.charge == 1:       text = '+'
            elif atom.charge > 1:      text = '{0:d}+'.format(atom.charge)
            elif atom.charge == -1:    text = u'\u2013'
            elif atom.charge < -1:     text = u'{0:d}\u2013'.format(abs(atom.charge))
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
            # Draw lone electron pairs
            # Draw them for nitrogen species only
            if drawLonePairs:
                for i in range (atom.lonePairs):
                    cr.new_sub_path()
                    if i == 0:
                        x1lp = x-2
                        y1lp = y-8
                        x2lp = x+2
                        y2lp = y-12
                    elif i == 1:
                        x1lp = x+12
                        y1lp = y-8
                        x2lp = x+8
                        y2lp = y-12
                    elif i == 2:
                        x1lp = x-2
                        y1lp = y-1
                        x2lp = x+2
                        y2lp = y+3
                    self.__drawLine(cr, x1lp, y1lp, x2lp, y2lp)
                
        # Update bounding rect to ensure atoms are included
        if boundingRect[0] < self.left:
            self.left = boundingRect[0]
        if boundingRect[1] < self.top:
            self.top = boundingRect[1]
        if boundingRect[2] > self.right:
            self.right = boundingRect[2]
        if boundingRect[3] > self.bottom:
            self.bottom = boundingRect[3]

################################################################################

class ReactionDrawer:
    """
    This class provides functionality for drawing chemical reactions using the
    skeletal formula of each reactant and product molecule via the Cairo 2D
    graphics engine. The most common use case is simply::
    
        ReactionDrawer().draw(reaction, format='png', path='reaction.png')
    
    where ``reaction`` is the :class:`Reaction` object to draw. You can also
    pass a dict of options to the constructor to affect how the molecules are
    drawn.
    """
    
    def __init__(self, options=None):
        self.options = MoleculeDrawer().options.copy()
        self.options.update({
            'arrowLength': 36,
        })
        if options: self.options.update(options)
    
    def draw(self, reaction, format, path=None):
        """
        Draw the given `reaction` using the given image `format` - pdf, svg, 
        ps, or png. If `path` is given, the drawing is saved to that location
        on disk.
        
        This function returns the Cairo surface and context used to create the
        drawing, as well as a bounding box for the molecule being drawn as the
        tuple (`left`, `top`, `width`, `height`).
        """
        # The Cairo 2D graphics library (and its Python wrapper) is required for
        # the reaction drawing algorithm
        try:
            import cairocffi as cairo
        except ImportError:
            try:
                import cairo
            except ImportError:
                print 'Cairo not found; molecule will not be drawn.'
                return

        from .molecule import Molecule
        from rmgpy.species import Species

        fontFamily = self.options['fontFamily']
        fontSizeNormal = self.options['fontSizeNormal']
        
        # First draw each of the reactants and products
        reactants = []; products = []
        for reactant in reaction.reactants:
            if isinstance(reactant, Species):
                molecule = reactant.molecule[0]
            elif isinstance(reactant, Molecule):
                molecule = reactant
            reactants.append(MoleculeDrawer().draw(molecule, format))
        for product in reaction.products:
            if isinstance(product, Species):
                molecule = product.molecule[0]
            elif isinstance(product, Molecule):
                molecule = product
            products.append(MoleculeDrawer().draw(molecule, format))
            
        # Next determine size required for surface
        rxn_width = 0; rxn_height = 0; rxn_top = 0
        for surface, cr, rect in reactants:
            left, top, width, height = rect
            rxn_width += width
            if height > rxn_height: rxn_height = height
            if height + top > rxn_top: rxn_top = height + top
        for surface, cr, rect in products:
            left, top, width, height = rect
            rxn_width += width
            if height > rxn_height: rxn_height = height
            if height + top > rxn_top: rxn_top = height + top
        
        rxn_top = 0.5 * rxn_height - rxn_top
        
        # Also include '+' and reaction arrow in width
        cr.set_font_size(fontSizeNormal)
        plus_extents = cr.text_extents(' + ')
        arrow_width = self.options['arrowLength']
        rxn_width += (len(reactants)-1) * plus_extents[4] + arrow_width + (len(products)-1) * plus_extents[4]
        
        # Now make the surface for the reaction and render each molecule on it
        rxn_surface = createNewSurface(format, path, width=rxn_width, height=rxn_height)
        rxn_cr = cairo.Context(rxn_surface)
        
        # Draw white background
        rxn_cr.set_source_rgba(1.0, 1.0, 1.0, 1.0)
        rxn_cr.paint()
    
        # Draw reactants
        rxn_x = 0.0; rxn_y = 0.0
        for index, reactant in enumerate(reactants):
            surface, cr, rect = reactant
            left, top, width, height = rect
            if index > 0:
                # Draw the "+" between the reactants
                rxn_cr.save()
                rxn_cr.set_font_size(fontSizeNormal)
                rxn_y = rxn_top + 0.5 * (rxn_height - plus_extents[3])
                rxn_cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                rxn_cr.move_to(rxn_x, rxn_y - plus_extents[1])
                rxn_cr.show_text(' + ')
                rxn_cr.restore()
                rxn_x += plus_extents[4]
            # Draw the reactant
            rxn_y = top + rxn_top + 0.5 * rxn_height
            rxn_cr.save()
            rxn_cr.set_source_surface(surface, rxn_x, rxn_y)
            rxn_cr.paint()
            rxn_cr.restore()
            rxn_x += width            
        
        # Draw reaction arrow
        # Unfortunately Cairo does not have arrow drawing built-in, so we must
        # draw the arrow head ourselves
        rxn_cr.save()
        rxn_cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        rxn_cr.set_line_width(1.0)
        rxn_cr.move_to(rxn_x + 8, rxn_top + 0.5 * rxn_height)
        rxn_cr.line_to(rxn_x + arrow_width - 8, rxn_top + 0.5 * rxn_height)
        rxn_cr.move_to(rxn_x + arrow_width - 14, rxn_top + 0.5 * rxn_height - 3.0)
        rxn_cr.line_to(rxn_x + arrow_width - 8, rxn_top + 0.5 * rxn_height)
        rxn_cr.line_to(rxn_x + arrow_width - 14, rxn_top + 0.5 * rxn_height + 3.0)
        rxn_cr.stroke()
        rxn_cr.restore()
        rxn_x += arrow_width
        
        # Draw products
        for index, product in enumerate(products):
            surface, cr, rect = product
            left, top, width, height = rect
            if index > 0:
                # Draw the "+" between the products
                rxn_cr.save()
                rxn_cr.set_font_size(fontSizeNormal)
                rxn_y = rxn_top + 0.5 * (rxn_height - plus_extents[3])
                rxn_cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                rxn_cr.move_to(rxn_x, rxn_y - plus_extents[1])
                rxn_cr.show_text(' + ')
                rxn_cr.restore()
                rxn_x += plus_extents[4]
            # Draw the product
            rxn_y = top + rxn_top + 0.5 * rxn_height
            rxn_cr.save()
            rxn_cr.set_source_surface(surface, rxn_x, rxn_y)
            rxn_cr.paint()
            rxn_cr.restore()
            rxn_x += width
        
        # Finish Cairo drawing
        if format == 'png':
            surface.write_to_png(path)
        else:
            surface.finish()
