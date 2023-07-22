#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2021 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

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

import logging
import math
import os.path
import re

try:
    import cairocffi as cairo
except ImportError:
    try:
        import cairo
    except ImportError:
        cairo = None
import numpy as np
from rdkit.Chem import AllChem

from rmgpy.molecule.molecule import Atom, Molecule
from rmgpy.qm.molecule import Geometry


################################################################################


def create_new_surface(file_format, target=None, width=1024, height=768):
    """
    Create a new surface of the specified `file_format`:
        "png" for :class:`ImageSurface`
        "svg" for :class:`SVGSurface`
        "pdf" for :class:`PDFSurface`
        "ps" for :class:`PSSurface`
    The surface will be written to the `target` parameter , which can be a
    path to save the surface to, or file-like object with a `write()` method.
    You can also optionally specify the `width` and `height` of the generated
    surface if you know what it is; otherwise a default size of 1024 by 768 is
    used.
    """
    file_format = file_format.lower()
    if file_format == "png":
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, int(width), int(height))
    elif file_format == "svg":
        surface = cairo.SVGSurface(target, width, height)
    elif file_format == "pdf":
        surface = cairo.PDFSurface(target, width, height)
    elif file_format == "ps":
        surface = cairo.PSSurface(target, width, height)
    else:
        raise ValueError('Invalid value "{0}" for type parameter; valid values are "png", "svg", "pdf", and "ps".'.format(type))
    return surface


################################################################################


class MoleculeDrawer(object):
    """
    This class provides functionality for drawing the skeletal formula of
    molecules using the Cairo 2D graphics engine. The most common use case is
    simply::

        MoleculeDrawer().draw(molecule, file_format='png', path='molecule.png')

    where ``molecule`` is the :class:`Molecule` object to draw. You can also
    pass a dict of options to the constructor to affect how the molecules are
    drawn.
    """

    def __init__(self, options=None):
        self.options = {
            "fontFamily": "sans",
            "fontSizeNormal": 12,
            "fontSizeSubscript": 8,
            "bondLength": 24,
            "padding": 2,
        }
        if options:
            self.options.update(options)

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

    def draw(self, molecule, file_format, target=None):
        """
        Draw the given `molecule` using the given image `file_format` - pdf, svg, ps, or
        png. If `path` is given, the drawing is saved to that location on disk. The
        `options` dict is an optional set of key-value pairs that can be used to
        control the generated drawing.

        This function returns the Cairo surface and context used to create the
        drawing, as well as a bounding box for the molecule being drawn as the
        tuple (`left`, `top`, `width`, `height`).
        """

        # The Cairo 2D graphics library (and its Python wrapper) is required for
        # the molecule drawing algorithm
        if cairo is None:
            logging.error("Cairo not found; molecule will not be drawn.")
            return

        # Make a copy of the molecule so we don't modify the original
        self.molecule = molecule.copy(deep=True)

        # Remove all unlabeled hydrogen atoms from the copied atoms and bonds, as
        # they are not drawn
        # However, if this would remove all atoms, then don't remove any
        atoms_to_remove = []
        self.implicitHydrogens = {}
        surface_sites = []
        for atom in self.molecule.atoms:
            if isinstance(atom, Atom) and atom.is_hydrogen() and atom.label == "":
                atoms_to_remove.append(atom)
            elif atom.is_surface_site():
                surface_sites.append(atom)
        if len(atoms_to_remove) < len(self.molecule.atoms) - len(surface_sites):
            for atom in atoms_to_remove:
                for atom2 in atom.bonds:
                    try:
                        self.implicitHydrogens[atom2] += 1
                    except KeyError:
                        self.implicitHydrogens[atom2] = 1
                self.molecule.remove_atom(atom)

        # Generate information about any cycles present in the molecule, as
        # they will need special attention
        self._find_ring_groups()
        # Handle carbon monoxide special case
        if self.molecule.get_formula() == "CO" and len(atoms_to_remove) == 0:
            # RDKit does not accept atom type O4tc
            for atom in self.molecule.atoms:
                if atom.symbol == "O":
                    self.molecule.remove_atom(atom)
            self.symbols = ["CO"]
            self.molecule.atoms[0].charge = 0  # don't label the C as - if you're not drawing the O with a +
            self.coordinates = np.array([[0, 0]], float)
        else:
            # Generate the coordinates to use to draw the molecule
            try:
                # before getting coordinates, make all bonds single and then
                # replace the bonds after generating coordinates. This avoids
                # bugs with RDKit
                old_bond_dictionary = self._make_single_bonds()
                self._generate_coordinates()
                self._replace_bonds(old_bond_dictionary)

                # Generate labels to use
                self._generate_atom_labels()

            except (ValueError, np.linalg.LinAlgError) as e:
                logging.error("Error while drawing molecule {0}: {1}".format(molecule.to_smiles(), e))
                import sys, traceback

                exc_type, exc_value, exc_traceback = sys.exc_info()
                traceback.print_exc()
                return None, None, None
            except KeyError:
                logging.error(
                    "KeyError occured when drawing molecule, likely because"
                    " the molecule contained non-standard bond orders in the"
                    " get_resonance_hybrid method. These cannot be drawn since"
                    " they cannot be sent to RDKit for coordinate placing."
                )
                raise

        self.coordinates[:, 1] *= -1
        self.coordinates *= self.options["bondLength"]

        # Handle some special cases
        if self.symbols == ["H", "H"]:
            # Render as H2 instead of H-H
            self.molecule.remove_atom(self.molecule.atoms[-1])
            self.symbols = ["H2"]
            self.coordinates = np.array([[0, 0]], float)
        elif molecule.is_isomorphic(Molecule(smiles="[O][O]")):
            # Render as O2 instead of O-O
            self.molecule.remove_atom(self.molecule.atoms[-1])
            self.molecule.atoms[0].radical_electrons = 0
            self.symbols = ["O2"]
            self.coordinates = np.array([[0, 0]], float)
        elif self.symbols == ["OH", "O"] or self.symbols == ["O", "OH"]:
            # Render as HO2 instead of HO-O or O-OH
            self.molecule.remove_atom(self.molecule.atoms[-1])
            self.symbols = ["O2H"]
            self.coordinates = np.array([[0, 0]], float)
        elif self.symbols == ["OH", "OH"]:
            # Render as H2O2 instead of HO-OH or O-OH
            self.molecule.remove_atom(self.molecule.atoms[-1])
            self.symbols = ["O2H2"]
            self.coordinates = np.array([[0, 0]], float)
        elif self.symbols == ["O", "C", "O"]:
            # Render as CO2 instead of O=C=O
            self.molecule.remove_atom(self.molecule.atoms[0])
            self.molecule.remove_atom(self.molecule.atoms[-1])
            self.symbols = ["CO2"]
            self.coordinates = np.array([[0, 0]], float)
        elif self.symbols == ["H", "H", "X"]:
            # Render as H2::X instead of crashing on H-H::X (vdW bond)
            self.molecule.remove_atom(self.molecule.atoms[0])
            self.symbols = ["H2", "X"]
            self.coordinates = np.array([[0, -0.5], [0, 0.5]], float) * self.options["bondLength"]

        # Create a dummy surface to draw to, since we don't know the bounding rect
        # We will copy this to another surface with the correct bounding rect
        surface0 = create_new_surface(file_format=file_format, target=None)
        cr0 = cairo.Context(surface0)

        # Render using Cairo
        self.render(cr0)

        # Create the real surface with the appropriate size
        xoff = self.left
        yoff = self.top
        width = self.right - self.left
        height = self.bottom - self.top
        self.surface = create_new_surface(file_format=file_format, target=target, width=width, height=height)
        self.cr = cairo.Context(self.surface)

        # Draw white background
        self.cr.set_source_rgba(1.0, 1.0, 1.0, 1.0)
        self.cr.paint()

        self.render(self.cr, offset=(-xoff, -yoff))

        if target is not None:
            # Finish Cairo drawing
            # Save PNG of drawing if appropriate
            if isinstance(target, str):
                ext = os.path.splitext(target)[1].lower()
                if ext == ".png":
                    self.surface.write_to_png(target)
                else:
                    self.surface.finish()
            else:
                self.surface.finish()

        return self.surface, self.cr, (xoff, yoff, width, height)

    def _find_ring_groups(self):
        """
        Find all of the cycles in the current molecule, and group them into
        sets of adjacent cycles.
        """

        # Find all of the cycles in the molecule
        self.cycles = self.molecule.get_smallest_set_of_smallest_rings()
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

    def _generate_coordinates(self):
        """
        Generate the 2D coordinates to be used when drawing the current
        molecule. The function uses rdKits 2D coordinate generation.
        Updates the self.coordinates Array in place.
        """
        atoms = self.molecule.atoms
        natoms = len(atoms)

        # Initialize array of coordinates
        self.coordinates = coordinates = np.zeros((natoms, 2))

        # If there are only one or two atoms to draw, then determining the
        # coordinates is trivial
        if natoms == 1:
            self.coordinates[0, :] = [0.0, 0.0]
            return self.coordinates
        elif natoms == 2:
            if atoms[0].is_surface_site():
                self.coordinates[0, :] = [0.0, -0.5]
                self.coordinates[1, :] = [0.0, 0.5]
            elif atoms[1].is_surface_site():
                self.coordinates[0, :] = [0.0, 0.5]
                self.coordinates[1, :] = [0.0, -0.5]
            else:
                self.coordinates[0, :] = [-0.5, 0.0]
                self.coordinates[1, :] = [0.5, 0.0]
            return self.coordinates

        # Decide whether we can use RDKit or have to generate coordinates ourselves
        for atom in self.molecule.atoms:
            if atom.charge != 0:
                use_rdkit = False
                break
        else:  # didn't break
            use_rdkit = True

        if not use_rdkit:
            if len(self.cycles) > 0:
                # Cyclic molecule
                backbone = self._find_cyclic_backbone()
                self._generate_ring_system_coordinates(backbone)
                # Flatten backbone so that it contains a list of the atoms in the
                # backbone, rather than a list of the cycles in the backbone
                backbone = list(set([atom for cycle in backbone for atom in cycle]))
            else:
                # Straight chain molecule
                backbone = self._find_straight_chain_backbone()
                self._generate_straight_chain_coordinates(backbone)

                # If backbone is linear, then rotate so that the bond is parallel to the
                # horizontal axis
                vector0 = coordinates[atoms.index(backbone[1]), :] - coordinates[atoms.index(backbone[0]), :]
                for i in range(2, len(backbone)):
                    vector = coordinates[atoms.index(backbone[i]), :] - coordinates[atoms.index(backbone[i - 1]), :]
                    if np.linalg.norm(vector - vector0) > 1e-4:
                        break
                else:
                    angle = math.atan2(vector0[0], vector0[1]) - math.pi / 2
                    rot = np.array([[math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]], float)
                    # need to keep self.coordinates and coordinates referring to the same object
                    self.coordinates = coordinates = np.dot(coordinates, rot)

            # If two atoms lie on top of each other, push them apart a bit
            # This is ugly, but at least the mess you end up with isn't as misleading
            # as leaving everything piled on top of each other at the origin
            import itertools

            for atom1, atom2 in itertools.combinations(backbone, 2):
                i1, i2 = atoms.index(atom1), atoms.index(atom2)
                if np.linalg.norm(coordinates[i1, :] - coordinates[i2, :]) < 0.5:
                    coordinates[i1, 0] -= 0.3
                    coordinates[i2, 0] += 0.3
                    coordinates[i1, 1] -= 0.2
                    coordinates[i2, 1] += 0.2

            # If two atoms lie on top of each other, push them apart a bit
            # This is ugly, but at least the mess you end up with isn't as misleading
            # as leaving everything piled on top of each other at the origin
            import itertools

            for atom1, atom2 in itertools.combinations(backbone, 2):
                i1, i2 = atoms.index(atom1), atoms.index(atom2)
                if np.linalg.norm(coordinates[i1, :] - coordinates[i2, :]) < 0.5:
                    coordinates[i1, 0] -= 0.3
                    coordinates[i2, 0] += 0.3
                    coordinates[i1, 1] -= 0.2
                    coordinates[i2, 1] += 0.2

            # Center backbone at origin
            xmin = np.min(coordinates[:, 0])
            xmax = np.max(coordinates[:, 0])
            ymin = np.min(coordinates[:, 1])
            ymax = np.max(coordinates[:, 1])
            xmid = 0.5 * (xmax + xmin)
            ymid = 0.5 * (ymax + ymin)
            for atom in backbone:
                index = atoms.index(atom)
                coordinates[index, 0] -= xmid
                coordinates[index, 1] -= ymid

            # We now proceed by calculating the coordinates of the functional groups
            # attached to the backbone
            # Each functional group is independent, although they may contain further
            # branching and cycles
            # In general substituents should try to grow away from the origin to
            # minimize likelihood of overlap
            self._generate_neighbor_coordinates(backbone)

        else:
            # Use RDKit 2D coordinate generation:

            # Generate the RDkit molecule from the RDkit molecule, use geometry
            # in order to match the atoms in the rdmol with the atoms in the
            # RMG molecule (which is required to extract coordinates).
            self.geometry = Geometry(None, None, self.molecule, None)

            rdmol, rd_atom_idx = self.geometry.rd_build()
            AllChem.Compute2DCoords(rdmol)

            # Extract the coordinates from each atom.
            for atom in atoms:
                index = rd_atom_idx[atom]
                point = rdmol.GetConformer(0).GetAtomPosition(index)
                coordinates[index, :] = [point.x * 0.6, point.y * 0.6]

            # RDKit generates some molecules more vertically than horizontally,
            # Especially linear ones. This will reflect any molecule taller than
            # it is wide across the line y=x
            ranges = np.ptp(coordinates, axis=0)
            if ranges[1] > ranges[0]:
                temp = np.copy(coordinates)
                coordinates[:, 0] = temp[:, 1]
                coordinates[:, 1] = temp[:, 0]

        # For surface species, rotate them so the site is at the bottom.
        if self.molecule.contains_surface_site():
            if len(self.molecule.atoms) == 1:
                return coordinates
            for site in self.molecule.atoms:
                if site.is_surface_site():
                    break
            else:
                raise Exception("Can't find surface site")
            if site.bonds:
                adsorbate = next(iter(site.bonds))
                vector0 = coordinates[atoms.index(site), :] - coordinates[atoms.index(adsorbate), :]
                angle = math.atan2(vector0[0], vector0[1]) - math.pi
                rot = np.array([[math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]], float)
                self.coordinates = coordinates = np.dot(coordinates, rot)
            else:
                # van der waals
                index = atoms.index(site)
                coordinates[index, 1] = min(coordinates[:, 1]) - 0.8  # just move the site down a bit
                coordinates[index, 0] = coordinates[:, 0].mean()  # and center it

    def _find_cyclic_backbone(self):
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

    def _find_straight_chain_backbone(self):
        """
        Return a set of atoms to use as the "backbone" of the molecule. For
        non-cyclics this is the largest straight chain between atoms. If carbon
        atoms are present, then we define the backbone only in terms of them.
        """
        # Find the terminal atoms - those that only have one explicit bond
        terminal_atoms = [atom for atom in self.molecule.atoms if len(atom.bonds) == 1]
        assert len(terminal_atoms) >= 2

        # Starting from each terminal atom, find the longest straight path to
        # another terminal
        # The longest found is the backbone
        backbone = []
        paths = []
        for atom in terminal_atoms:
            paths.extend(self._find_straight_chain_paths([atom]))

        # Remove any paths that don't end in a terminal atom
        # (I don't think this should remove any!)
        paths = [path for path in paths if path[-1] in terminal_atoms]

        # Remove all paths shorter than the maximum
        length = max([len(path) for path in paths])
        paths = [path for path in paths if len(path) == length]

        # Prefer the paths with the most carbon atoms
        carbons = [sum([1 for atom in path if atom.is_carbon()]) for path in paths]
        max_carbons = max(carbons)
        paths = [path for path, carbon in zip(paths, carbons) if carbon == max_carbons]

        # At this point we could choose any remaining path, so simply choose the first
        backbone = paths[0]

        assert len(backbone) > 1
        assert backbone[0] in terminal_atoms
        assert backbone[-1] in terminal_atoms

        return backbone

    def _find_straight_chain_paths(self, atoms0):
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
                if not self.molecule.is_atom_in_cycle(atom2):
                    paths.extend(self._find_straight_chain_paths(atoms))
        if len(paths) == 0:
            paths.append(atoms0[:])
        return paths

    def _generate_ring_system_coordinates(self, atoms):
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
        angle = -2 * math.pi / len(cycle)
        radius = 1.0 / (2 * math.sin(math.pi / len(cycle)))
        for i, atom in enumerate(cycle):
            index = self.molecule.atoms.index(atom)
            coordinates[index, :] = [math.cos(math.pi / 2 + i * angle), math.sin(math.pi / 2 + i * angle)]
            coordinates[index, :] *= radius
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
                    if count == 1 or count == 2:
                        if cycle is None or len(cycle0) > len(cycle):
                            cycle = cycle0
            cycle0 = cycle1
            if cycle is None:
                break
            atoms.remove(cycle)

            # Shuffle atoms in cycle such that the common atoms come first
            # Also find the average center of the processed cycles that touch the
            # current cycles
            found = False
            common_atoms = []
            count = 0
            center0 = np.zeros(2, float)
            for cycle1 in processed:
                found = False
                for atom in cycle1:
                    if atom in cycle and atom not in common_atoms:
                        common_atoms.append(atom)
                        found = True
                if found:
                    center1 = np.zeros(2, float)
                    for atom in cycle1:
                        center1 += coordinates[cycle1.index(atom), :]
                    center1 /= len(cycle1)
                    center0 += center1
                    count += 1
            center0 /= count

            if len(common_atoms) > 1:
                index0 = cycle.index(common_atoms[0])
                index1 = cycle.index(common_atoms[1])
                if (index0 == 0 and index1 == len(cycle) - 1) or (index1 == 0 and index0 == len(cycle) - 1):
                    cycle = cycle[-1:] + cycle[0:-1]
                if cycle.index(common_atoms[1]) < cycle.index(common_atoms[0]):
                    cycle.reverse()
                index = cycle.index(common_atoms[0])
                cycle = cycle[index:] + cycle[0:index]

            # Determine center of cycle based on already-assigned positions of
            # common atoms (which won't be changed)
            if len(common_atoms) == 1 or len(common_atoms) == 2:
                # Center of new cycle is reflection of center of adjacent cycle
                # across common atom or bond
                center = np.zeros(2, float)
                for atom in common_atoms:
                    center += coordinates[self.molecule.atoms.index(atom), :]
                center /= len(common_atoms)
                vector = center - center0
                center += vector
                radius = 1.0 / (2 * math.sin(math.pi / len(cycle)))

            else:
                # Use any three points to determine the point equidistant from these
                # three; this is the center
                index0 = self.molecule.atoms.index(common_atoms[0])
                index1 = self.molecule.atoms.index(common_atoms[1])
                index2 = self.molecule.atoms.index(common_atoms[2])
                A = np.zeros((2, 2), float)
                b = np.zeros((2), float)
                A[0, :] = 2 * (coordinates[index1, :] - coordinates[index0, :])
                A[1, :] = 2 * (coordinates[index2, :] - coordinates[index0, :])
                b[0] = coordinates[index1, 0] ** 2 + coordinates[index1, 1] ** 2 - coordinates[index0, 0] ** 2 - coordinates[index0, 1] ** 2
                b[1] = coordinates[index2, 0] ** 2 + coordinates[index2, 1] ** 2 - coordinates[index0, 0] ** 2 - coordinates[index0, 1] ** 2
                center = np.linalg.solve(A, b)
                radius = np.linalg.norm(center - coordinates[index0, :])

            start_angle = 0.0
            end_angle = 0.0
            if len(common_atoms) == 1:
                # We will use the full 360 degrees to place the other atoms in the cycle
                start_angle = math.atan2(-vector[1], vector[0])
                end_angle = start_angle + 2 * math.pi
            elif len(common_atoms) >= 2:
                # Divide other atoms in cycle equally among unused angle
                vector = coordinates[cycle.index(common_atoms[-1]), :] - center
                start_angle = math.atan2(vector[1], vector[0])
                vector = coordinates[cycle.index(common_atoms[0]), :] - center
                end_angle = math.atan2(vector[1], vector[0])

            # Place remaining atoms in cycle
            if end_angle < start_angle:
                end_angle += 2 * math.pi
                d_angle = (end_angle - start_angle) / (len(cycle) - len(common_atoms) + 1)
            else:
                end_angle -= 2 * math.pi
                d_angle = (end_angle - start_angle) / (len(cycle) - len(common_atoms) + 1)

            count = 1
            for i in range(len(common_atoms), len(cycle)):
                angle = start_angle + count * d_angle
                index = self.molecule.atoms.index(cycle[i])
                # Check that we aren't reassigning any atom positions
                # This version assumes that no atoms belong at the origin, which is
                # usually fine because the first ring is centered at the origin
                if np.linalg.norm(coordinates[index, :]) < 1e-4:
                    vector = np.array([math.cos(angle), math.sin(angle)], float)
                    coordinates[index, :] = center + radius * vector
                count += 1

            # We're done assigning coordinates for this cycle, so mark it as processed
            processed.append(cycle)

    def _generate_straight_chain_coordinates(self, atoms):
        """
        Update the coordinates for the linear straight chain of `atoms` in
        the current molecule.
        """
        coordinates = self.coordinates

        # First atom goes at origin
        index0 = self.molecule.atoms.index(atoms[0])
        coordinates[index0, :] = [0.0, 0.0]

        # Second atom goes on x-axis (for now; this could be improved!)
        index1 = self.molecule.atoms.index(atoms[1])
        vector = np.array([1.0, 0.0], float)
        if atoms[0].bonds[atoms[1]].is_triple():
            rotate_positive = False
        else:
            rotate_positive = True
            rot = np.array([[math.cos(-math.pi / 6), math.sin(-math.pi / 6)], [-math.sin(-math.pi / 6), math.cos(-math.pi / 6)]], float)
            vector = np.array([1.0, 0.0], float)
            vector = np.dot(rot, vector)
        coordinates[index1, :] = coordinates[index0, :] + vector

        # Other atoms
        for i in range(2, len(atoms)):
            atom0 = atoms[i - 2]
            atom1 = atoms[i - 1]
            atom2 = atoms[i]
            index1 = self.molecule.atoms.index(atom1)
            index2 = self.molecule.atoms.index(atom2)
            bond0 = atom0.bonds[atom1]
            bond = atom1.bonds[atom2]
            # Angle of next bond depends on the number of bonds to the start atom
            num_bonds = len(atom1.bonds)
            if num_bonds == 2:
                if (bond0.is_triple() or bond.is_triple()) or (bond0.is_double() and bond.is_double()):
                    # Rotate by 0 degrees towards horizontal axis (to get angle of 180)
                    angle = 0.0
                else:
                    # Rotate by 60 degrees towards horizontal axis (to get angle of 120)
                    angle = math.pi / 3
            elif num_bonds == 3:
                # Rotate by 60 degrees towards horizontal axis (to get angle of 120)
                angle = math.pi / 3
            elif num_bonds == 4:
                # Rotate by 0 degrees towards horizontal axis (to get angle of 90)
                angle = 0.0
            elif num_bonds == 5:
                # Rotate by 36 degrees towards horizontal axis (to get angle of 144)
                angle = math.pi / 5
            elif num_bonds == 6:
                # Rotate by 0 degrees towards horizontal axis (to get angle of 180)
                angle = 0.0
            # Determine coordinates for atom
            if angle != 0:
                if not rotate_positive:
                    angle = -angle
                rot = np.array([[math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]], float)
                vector = np.dot(rot, vector)
                rotate_positive = not rotate_positive
            coordinates[index2, :] = coordinates[index1, :] + vector

    def _generate_neighbor_coordinates(self, backbone):
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
            bond_angles = []
            for atom1 in atom0.bonds:
                index1 = atoms.index(atom1)
                if atom1 in backbone:
                    vector = coordinates[index1, :] - coordinates[index0, :]
                    angle = math.atan2(vector[1], vector[0])
                    bond_angles.append(angle)
            bond_angles.sort()

            best_angle = 2 * math.pi / len(atom0.bonds)
            regular = True
            for angle1, angle2 in zip(bond_angles[0:-1], bond_angles[1:]):
                if all([abs(angle2 - angle1 - (i + 1) * best_angle) > 1e-4 for i in range(len(atom0.bonds))]):
                    regular = False

            if regular:
                # All the bonds around each atom are equally spaced
                # We just need to fill in the missing bond locations

                # Determine rotation angle and matrix
                rot = np.array([[math.cos(best_angle), -math.sin(best_angle)], [math.sin(best_angle), math.cos(best_angle)]], float)
                # Determine the vector of any currently-existing bond from this atom
                vector = None
                for atom1 in atom0.bonds:
                    index1 = atoms.index(atom1)
                    if atom1 in backbone or np.linalg.norm(coordinates[index1, :]) > 1e-4:
                        vector = coordinates[index1, :] - coordinates[index0, :]

                # Iterate through each neighboring atom to this backbone atom
                # If the neighbor is not in the backbone and does not yet have
                # coordinates, then we need to determine coordinates for it
                for atom1 in atom0.bonds:
                    if atom1 not in backbone and np.linalg.norm(coordinates[atoms.index(atom1), :]) < 1e-4:
                        occupied = True
                        count = 0
                        # Rotate vector until we find an unoccupied location
                        while occupied and count < len(atom0.bonds):
                            count += 1
                            occupied = False
                            vector = np.dot(rot, vector)
                            for atom2 in atom0.bonds:
                                index2 = atoms.index(atom2)
                                if np.linalg.norm(coordinates[index2, :] - coordinates[index0, :] - vector) < 1e-4:
                                    occupied = True
                        coordinates[atoms.index(atom1), :] = coordinates[index0, :] + vector
                        self._generate_functional_group_coordinates(atom0, atom1)

            else:
                # The bonds are not evenly spaced (e.g. due to a ring)
                # We place all of the remaining bonds evenly over the reflex angle
                start_angle = max(bond_angles)
                end_angle = min(bond_angles)
                if 0.0 < end_angle - start_angle < math.pi:
                    end_angle += 2 * math.pi
                elif 0.0 > end_angle - start_angle > -math.pi:
                    start_angle -= 2 * math.pi
                d_angle = (end_angle - start_angle) / (len(atom0.bonds) - len(bond_angles) + 1)

                index = 1
                for atom1 in atom0.bonds:
                    if atom1 not in backbone and np.linalg.norm(coordinates[atoms.index(atom1), :]) < 1e-4:
                        angle = start_angle + index * d_angle
                        index += 1
                        vector = np.array([math.cos(angle), math.sin(angle)], float)
                        vector /= np.linalg.norm(vector)
                        coordinates[atoms.index(atom1), :] = coordinates[index0, :] + vector
                        self._generate_functional_group_coordinates(atom0, atom1)

    def _generate_functional_group_coordinates(self, atom0, atom1):
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
        vector = coordinates[index0, :] - coordinates[index1, :]
        bond_angle = math.atan2(vector[1], vector[0])

        # Check to see if atom1 is in any cycles in the molecule
        ring_system = None
        for ring_sys in self.ringSystems:
            if any([atom1 in ring for ring in ring_sys]):
                ring_system = ring_sys

        if ring_system is not None:
            # atom1 is part of a ring system, so we need to process the entire
            # ring system at once

            # Generate coordinates for all atoms in the ring system
            self._generate_ring_system_coordinates(ring_system)

            cycle_atoms = list(set([atom for ring in ring_system for atom in ring]))

            coordinates_cycle = np.zeros_like(self.coordinates)
            for atom in cycle_atoms:
                coordinates_cycle[atoms.index(atom), :] = coordinates[atoms.index(atom), :]

            # Rotate the ring system coordinates so that the line connecting atom1
            # and the center of mass of the ring is parallel to that between
            # atom0 and atom1
            center = np.zeros(2, float)
            for atom in cycle_atoms:
                center += coordinates_cycle[atoms.index(atom), :]
            center /= len(cycle_atoms)
            vector0 = center - coordinates_cycle[atoms.index(atom1), :]
            angle = math.atan2(vector[1] - vector0[1], vector[0] - vector0[0])
            rot = np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]], float)
            coordinates_cycle = np.dot(coordinates_cycle, rot)

            # Translate the ring system coordinates to the position of atom1
            coordinates_cycle += coordinates[atoms.index(atom1), :] - coordinates_cycle[atoms.index(atom1), :]
            for atom in cycle_atoms:
                coordinates[atoms.index(atom), :] = coordinates_cycle[atoms.index(atom), :]

            # Generate coordinates for remaining neighbors of ring system,
            # continuing to recurse as needed
            self._generate_neighbor_coordinates(cycle_atoms)

        else:
            # atom1 is not in any rings, so we can continue as normal

            # Determine rotation angle and matrix
            num_bonds = len(atom1.bonds)
            angle = 0.0
            if num_bonds == 2:
                bond0, bond = list(atom1.bonds.values())
                if (bond0.is_triple() or bond.is_triple()) or (bond0.is_double() and bond.is_double()):
                    angle = math.pi
                else:
                    angle = 2 * math.pi / 3
                    # Make sure we're rotating such that we move away from the origin,
                    # to discourage overlap of functional groups
                    rot1 = np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]], float)
                    rot2 = np.array([[math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]], float)
                    vector1 = coordinates[index1, :] + np.dot(rot1, vector)
                    vector2 = coordinates[index1, :] + np.dot(rot2, vector)
                    if bond_angle < -0.5 * math.pi or bond_angle > 0.5 * math.pi:
                        angle = abs(angle)
                    else:
                        angle = -abs(angle)
            else:
                angle = 2 * math.pi / num_bonds
            rot = np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]], float)

            # Iterate through each neighboring atom to this backbone atom
            # If the neighbor is not in the backbone, then we need to determine
            # coordinates for it
            for atom, bond in atom1.bonds.items():
                if atom is not atom0:
                    occupied = True
                    count = 0
                    # Rotate vector until we find an unoccupied location
                    while occupied and count < len(atom1.bonds):
                        count += 1
                        occupied = False
                        vector = np.dot(rot, vector)
                        for atom2 in atom1.bonds:
                            index2 = atoms.index(atom2)
                            if np.linalg.norm(coordinates[index2, :] - coordinates[index1, :] - vector) < 1e-4:
                                occupied = True
                    coordinates[atoms.index(atom), :] = coordinates[index1, :] + vector

                    # Recursively continue with functional group
                    self._generate_functional_group_coordinates(atom1, atom)

    def _generate_atom_labels(self):
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
            # or they are isotopically labeled
            if symbols[i] == "C" and len(symbols) > 2:
                if (len(atoms[i].bonds) > 1 or (atoms[i].radical_electrons == 0 and atoms[i].charge == 0)) and atoms[i].element.isotope == -1:
                    symbols[i] = ""
        # Do label atoms that have only double bonds to one or more labeled atoms
        changed = True
        while changed:
            changed = False
            for i in range(len(symbols)):
                if (
                    symbols[i] == ""
                    and all([(bond.is_double() or bond.is_triple()) for bond in atoms[i].bonds.values()])
                    and any([symbols[atoms.index(atom)] != "" for atom in atoms[i].bonds])
                ):
                    symbols[i] = atoms[i].symbol
                    changed = True
        # Add implicit hydrogens
        for i in range(len(symbols)):
            if symbols[i] != "":
                try:
                    h_count = self.implicitHydrogens[atoms[i]]
                except KeyError:
                    continue
                if h_count == 1:
                    symbols[i] = symbols[i] + "H"
                elif h_count > 1:
                    symbols[i] = symbols[i] + "H{0:d}".format(h_count)

        return symbols

    def render(self, cr, offset=None):
        """
        Uses the Cairo graphics library to create a skeletal formula drawing of a
        molecule containing the list of `atoms` and dict of `bonds` to be drawn.
        The 2D position of each atom in `atoms` is given in the `coordinates` array.
        The symbols to use at each atomic position are given by the list `symbols`.
        You must specify the Cairo context `cr` to render to.
        """

        coordinates = self.coordinates
        atoms = self.molecule.atoms
        symbols = self.symbols

        draw_lone_pairs = False

        for atom in atoms:
            if isinstance(atom, Atom) and atom.is_nitrogen():
                draw_lone_pairs = True

        left = 0.0
        top = 0.0
        right = 0.0
        bottom = 0.0

        # Shift coordinates by offset value
        if offset is not None:
            coordinates[:, 0] += offset[0]
            coordinates[:, 1] += offset[1]

        # Draw bonds
        for atom1 in atoms:
            for atom2, bond in atom1.bonds.items():
                index1 = atoms.index(atom1)
                index2 = atoms.index(atom2)
                if index1 < index2:  # So we only draw each bond once
                    self._render_bond(index1, index2, bond, cr)

        # Draw aromatic bonds
        for cycle in self.cycles:
            cycle_bonds = []
            for atom1, atom2 in zip(cycle[0:-1], cycle[1:]):
                cycle_bonds.append(atom1.bonds[atom2])
            cycle_bonds.append(cycle[0].bonds[cycle[-1]])
            if all([bond.is_benzene() for bond in cycle_bonds]):
                # We've found an aromatic ring, so draw a circle in the center to represent the benzene bonds
                center = np.zeros(2, float)
                for atom in cycle:
                    index = atoms.index(atom)
                    center += coordinates[index, :]
                center /= len(cycle)
                index1 = atoms.index(cycle[0])
                index2 = atoms.index(cycle[1])
                radius = (
                    math.sqrt(
                        (center[0] - (coordinates[index1, 0] + coordinates[index2, 0]) / 2) ** 2
                        + (center[1] - (coordinates[index1, 1] + coordinates[index2, 1]) / 2) ** 2
                    )
                    - 4
                )
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.set_line_width(1.0)
                cr.set_line_cap(cairo.LINE_CAP_ROUND)
                cr.arc(center[0], center[1], radius, 0.0, 2 * math.pi)
                cr.stroke()

        # Draw atoms
        for i, atom in enumerate(atoms):
            symbol = symbols[i]
            index = atoms.index(atom)
            x0, y0 = coordinates[index, :]
            vector = np.zeros(2, float)
            for atom2 in atom.bonds:
                vector += coordinates[atoms.index(atom2), :] - coordinates[index, :]
            heavy_first = vector[0] <= 0
            if len(atoms) == 1 and atoms[0].symbol not in ["C", "N"] and atoms[0].charge == 0 and atoms[0].radical_electrons == 0:
                # This is so e.g. water is rendered as H2O rather than OH2
                heavy_first = False
                cr.set_font_size(self.options["fontSizeNormal"])
                x0 += cr.text_extents(symbols[0])[2] / 2.0
            self._render_atom(symbol, atom, x0, y0, cr, heavy_first, draw_lone_pairs)

        # Add a small amount of whitespace on all sides
        padding = self.options["padding"]
        self.left -= padding
        self.top -= padding
        self.right += padding
        self.bottom += padding

    def _draw_line(self, cr, x1, y1, x2, y2, dashed=False, dash_sizes=None):
        """
        Draw a line on the given Cairo context `cr` from (`x1`, `y1`) to
        (`x2`,`y2`), and update the bounding rectangle if necessary.

        For a dashed line set ``dashed=True``.
        Then the optional `dash_sizes` can be a list of on/off segment lengths,
        which defaults to [3.5, 3.5] if not specified.
        """

        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.set_line_width(1.0)
        if dashed:
            if dash_sizes is None:
                dash_sizes = [3.5, 3.5]
            cr.set_dash(dash_sizes)
        cr.set_line_cap(cairo.LINE_CAP_ROUND)
        cr.move_to(x1, y1)
        cr.line_to(x2, y2)
        cr.stroke()
        # remove dashes for next method call
        if dashed:
            cr.set_dash([])
        if x1 < self.left:
            self.left = x1
        if x1 > self.right:
            self.right = x1
        if y1 < self.top:
            self.top = y1
        if y1 > self.bottom:
            self.bottom = y1
        if x2 < self.left:
            self.left = x2
        if x2 > self.right:
            self.right = x2
        if y2 < self.top:
            self.top = y2
        if y2 > self.bottom:
            self.bottom = y2

    def _render_bond(self, atom1, atom2, bond, cr):
        """
        Render an individual `bond` between atoms with indices `atom1` and `atom2`
        on the Cairo context `cr`.
        """

        bond_length = self.options["bondLength"]

        # determine if aromatic
        is_aromatic = False
        for cycle in self.cycles:
            if self.molecule.atoms[atom1] in cycle and self.molecule.atoms[atom2] in cycle:
                all_benzenes = True
                for index in range(len(cycle)):
                    if not cycle[index - 1].bonds[cycle[index]].is_benzene():
                        all_benzenes = False
                        break
                if all_benzenes:
                    is_aromatic = True
                    break

        x1, y1 = self.coordinates[atom1, :]
        x2, y2 = self.coordinates[atom2, :]
        angle = math.atan2(y2 - y1, x2 - x1)

        dx = x2 - x1
        dy = y2 - y1
        du = math.cos(angle + math.pi / 2)
        dv = math.sin(angle + math.pi / 2)
        if self.symbols[atom1] != "" or self.symbols[atom2] != "":
            if bond.is_quadruple():
                # Draw quadruple bond centered on bond axis
                du *= 1.5
                dv *= 1.5
                self._draw_line(cr, x1 - du, y1 - dv, x2 - du, y2 - dv)
                self._draw_line(cr, x1 + du, y1 + dv, x2 + du, y2 + dv)
                du *= 2.2
                dv *= 2.2
                self._draw_line(cr, x1 - du, y1 - dv, x2 - du, y2 - dv)
                self._draw_line(cr, x1 + du, y1 + dv, x2 + du, y2 + dv)
            elif bond.is_triple():
                # Draw triple bond centered on bond axis
                du *= 3
                dv *= 3
                self._draw_line(cr, x1 - du, y1 - dv, x2 - du, y2 - dv)
                self._draw_line(cr, x1, y1, x2, y2)
                self._draw_line(cr, x1 + du, y1 + dv, x2 + du, y2 + dv)
            elif 2 < bond.get_order_num() < 3:
                du *= 3
                dv *= 3
                self._draw_line(cr, x1 - du, y1 - dv, x2 - du, y2 - dv)
                self._draw_line(cr, x1, y1, x2, y2)
                self._draw_line(cr, x1 + du, y1 + dv, x2 + du, y2 + dv, dashed=True)
            elif bond.is_double():
                # Draw double bond centered on bond axis
                du *= 1.6
                dv *= 1.6
                self._draw_line(cr, x1 - du, y1 - dv, x2 - du, y2 - dv)
                self._draw_line(cr, x1 + du, y1 + dv, x2 + du, y2 + dv)
            elif 1 < bond.get_order_num() < 2 and not is_aromatic:
                # Draw dashed double bond centered on bond axis
                du *= 1.6
                dv *= 1.6
                self._draw_line(cr, x1 - du, y1 - dv, x2 - du, y2 - dv)
                self._draw_line(cr, x1 + du, y1 + dv, x2 + du, y2 + dv, dashed=True)
            elif bond.is_hydrogen_bond():
                # Draw a dashed line
                self._draw_line(cr, x1, y1, x2, y2, dashed=True, dash_sizes=[0.5, 3.5])
            else:
                self._draw_line(cr, x1, y1, x2, y2)
        else:
            # Draw bond on skeleton
            self._draw_line(cr, x1, y1, x2, y2)
            # Draw other bonds
            if bond.is_double():
                du *= 3.2
                dv *= 3.2
                dx = 2 * dx / bond_length
                dy = 2 * dy / bond_length
                self._draw_line(cr, x1 + du + dx, y1 + dv + dy, x2 + du - dx, y2 + dv - dy)
            elif bond.is_triple():
                du *= 3
                dv *= 3
                dx = 2 * dx / bond_length
                dy = 2 * dy / bond_length
                self._draw_line(cr, x1 - du + dx, y1 - dv + dy, x2 - du - dx, y2 - dv - dy)
                self._draw_line(cr, x1 + du + dx, y1 + dv + dy, x2 + du - dx, y2 + dv - dy)
            elif 1 < bond.get_order_num() < 2 and not is_aromatic:
                du *= 3.2
                dv *= 3.2
                dx = 2 * dx / bond_length
                dy = 2 * dy / bond_length
                self._draw_line(cr, x1 + du + dx, y1 + dv + dy, x2 + du - dx, y2 + dv - dy, dashed=True)
            elif 2 < bond.get_order_num() < 3:
                du *= 3
                dv *= 3
                dx = 2 * dx / bond_length
                dy = 2 * dy / bond_length
                self._draw_line(cr, x1 - du + dx, y1 - dv + dy, x2 - du - dx, y2 - dv - dy)
                self._draw_line(cr, x1 + du + dx, y1 + dv + dy, x2 + du - dx, y2 + dv - dy, dashed=True)
            elif bond.is_quadruple():
                du *= 3
                dv *= 3
                dx = 2 * dx / bond_length
                dy = 2 * dy / bond_length
                self._draw_line(cr, x1 - du + dx, y1 - dv + dy, x2 - du - dx, y2 - dv - dy)
                self._draw_line(cr, x1 + du + dx, y1 + dv + dy, x2 + du - dx, y2 + dv - dy)
                self._draw_line(cr, x1 + 2 * du + dx, y1 + 2 * dv + dy, x2 + 2 * du - dx, y2 + 2 * dv - dy)

    def _render_atom(self, symbol, atom, x0, y0, cr, heavy_first=True, draw_lone_pairs=False):
        """
        Render the `label` for an atom centered around the coordinates (`x0`, `y0`)
        onto the Cairo context `cr`. If `heavyFirst` is ``False``, then the order
        of the atoms will be reversed in the symbol. This method also causes
        radical electrons and charges to be drawn adjacent to the rendered symbol.
        """

        atoms = self.molecule.atoms

        if symbol != "":
            heavy_atom = symbol[0]

            # Split label by atoms
            labels = re.findall(r"[A-Z][a-z]*[0-9]*", symbol)
            if not heavy_first:
                labels.reverse()
            if "C" not in symbol and "O" not in symbol and len(atoms) == 1:
                labels.sort()
            symbol = "".join(labels)

            # Determine positions of each character in the symbol
            coordinates = []

            cr.set_font_size(self.options["fontSizeNormal"])
            y0 += max([cr.text_extents(char)[3] for char in symbol if char.isalpha()]) / 2

            for i, label in enumerate(labels):
                for j, char in enumerate(label):
                    cr.set_font_size(self.options["fontSizeSubscript" if char.isdigit() else "fontSizeNormal"])
                    xbearing, ybearing, width, height, xadvance, yadvance = cr.text_extents(char)
                    if i == 0 and j == 0:
                        # Center heavy atom at (x0, y0)
                        x = x0 - width / 2.0 - xbearing
                        y = y0
                    else:
                        # Left-justify other atoms (for now)
                        x = x0
                        y = y0
                    if char.isdigit():
                        y += height / 2.0
                    coordinates.append((x, y))
                    x0 = x + xadvance

            x = y = 1000000
            width = height = 0
            start_width = end_width = 0
            for i, char in enumerate(symbol):
                cr.set_font_size(self.options["fontSizeSubscript" if char.isdigit() else "fontSizeNormal"])
                extents = cr.text_extents(char)
                if coordinates[i][0] + extents[0] < x:
                    x = coordinates[i][0] + extents[0]
                if coordinates[i][1] + extents[1] < y:
                    y = coordinates[i][1] + extents[1]
                width += extents[4] if i < len(symbol) - 1 else extents[2]
                if extents[3] > height:
                    height = extents[3]
                if i == 0:
                    start_width = extents[2]
                if i == len(symbol) - 1:
                    end_width = extents[2]

            if not heavy_first:
                for i in range(len(coordinates)):
                    coordinates[i] = (coordinates[i][0] - (width - start_width / 2 - end_width / 2), coordinates[i][1])
                x -= width - start_width / 2 - end_width / 2

            # Background
            x1 = x - 2
            y1 = y - 2
            x2 = x + width + 2
            y2 = y + height + 2
            r = 4
            cr.move_to(x1 + r, y1)
            cr.line_to(x2 - r, y1)
            cr.curve_to(x2 - r / 2, y1, x2, y1 + r / 2, x2, y1 + r)
            cr.line_to(x2, y2 - r)
            cr.curve_to(x2, y2 - r / 2, x2 - r / 2, y2, x2 - r, y2)
            cr.line_to(x1 + r, y2)
            cr.curve_to(x1 + r / 2, y2, x1, y2 - r / 2, x1, y2 - r)
            cr.line_to(x1, y1 + r)
            cr.curve_to(x1, y1 + r / 2, x1 + r / 2, y1, x1 + r, y1)
            cr.close_path()

            cr.save()
            cr.set_operator(cairo.OPERATOR_SOURCE)
            cr.set_source_rgba(1.0, 1.0, 1.0, 1.0)
            cr.fill()
            cr.restore()

            bounding_rect = [x1, y1, x2, y2]

            # Set color for text
            if not isinstance(atom, Atom) or atom.element.isotope != -1:
                cr.set_source_rgba(0.0, 0.5, 0.0, 1.0)
            elif heavy_atom == "C":
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            elif heavy_atom == "N":
                cr.set_source_rgba(0.0, 0.0, 1.0, 1.0)
            elif heavy_atom == "O":
                cr.set_source_rgba(1.0, 0.0, 0.0, 1.0)
            elif heavy_atom == "F":
                cr.set_source_rgba(0.5, 0.75, 1.0, 1.0)
            elif heavy_atom == "Si":
                cr.set_source_rgba(0.5, 0.5, 0.75, 1.0)
            elif heavy_atom == "Al":
                cr.set_source_rgba(0.75, 0.5, 0.5, 1.0)
            elif heavy_atom == "P":
                cr.set_source_rgba(1.0, 0.5, 0.0, 1.0)
            elif heavy_atom == "S":
                cr.set_source_rgba(1.0, 0.75, 0.5, 1.0)
            elif heavy_atom == "Cl":
                cr.set_source_rgba(0.0, 1.0, 0.0, 1.0)
            elif heavy_atom == "Br":
                cr.set_source_rgba(0.6, 0.2, 0.2, 1.0)
            elif heavy_atom == "I":
                cr.set_source_rgba(0.5, 0.0, 0.5, 1.0)
            elif heavy_atom == "X":
                cr.set_source_rgba(0.5, 0.25, 0.5, 1.0)
            else:
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)

            # Text itself
            for i, char in enumerate(symbol):
                cr.set_font_size(self.options["fontSizeSubscript" if char.isdigit() else "fontSizeNormal"])
                xbearing, ybearing, width, height, xadvance, yadvance = cr.text_extents(char)
                xi, yi = coordinates[i]
                cr.move_to(xi, yi)
                cr.show_text(char)

            x, y = coordinates[0] if heavy_first else coordinates[-1]

        else:
            x, y = x0, y0
            width = height = 0
            bounding_rect = [x0 - 0.5, y0 - 0.5, x0 + 0.5, y0 + 0.5]
            heavy_atom = ""

        # Draw radical electrons and charges
        # These will be placed either horizontally along the top or bottom of the
        # atom or vertically along the left or right of the atom
        orientation = " "
        if len(atom.bonds) == 0:
            if len(symbol) == 1:
                orientation = "r"
            else:
                orientation = "l"
        elif len(atom.bonds) == 1:
            # Terminal atom - we require a horizontal arrangement if there are
            # more than just the heavy atom
            atom1 = next(iter(atom.bonds))
            vector = self.coordinates[atoms.index(atom), :] - self.coordinates[atoms.index(atom1), :]
            if len(symbol) <= 1:
                angle = math.atan2(vector[1], vector[0])
                if 3 * math.pi / 4 <= angle or angle < -3 * math.pi / 4:
                    orientation = "l"
                elif -3 * math.pi / 4 <= angle < -1 * math.pi / 4:
                    orientation = "b"
                elif -1 * math.pi / 4 <= angle < 1 * math.pi / 4:
                    orientation = "r"
                else:
                    orientation = "t"
            else:
                if vector[1] <= 0:
                    orientation = "b"
                else:
                    orientation = "t"
        else:
            # Internal atom
            # First try to see if there is a "preferred" side on which to place the
            # radical/charge data, i.e. if the bonds are unbalanced
            vector = np.zeros(2, float)
            for atom1 in atom.bonds:
                vector += self.coordinates[atoms.index(atom), :] - self.coordinates[atoms.index(atom1), :]
            if np.linalg.norm(vector) < 1e-4:
                # All of the bonds are balanced, so we'll need to be more shrewd
                angles = []
                for atom1 in atom.bonds:
                    vector = self.coordinates[atoms.index(atom1), :] - self.coordinates[atoms.index(atom), :]
                    angles.append(math.atan2(vector[1], vector[0]))
                # Try one more time to see if we can use one of the four sides
                # (due to there being no bonds in that quadrant)
                # We don't even need a full 90 degrees open (using 60 degrees instead)
                if all([1 * math.pi / 3 >= angle or angle >= 2 * math.pi / 3 for angle in angles]):
                    orientation = "t"
                elif all([-2 * math.pi / 3 >= angle or angle >= -1 * math.pi / 3 for angle in angles]):
                    orientation = "b"
                elif all([-1 * math.pi / 6 >= angle or angle >= 1 * math.pi / 6 for angle in angles]):
                    orientation = "r"
                elif all([5 * math.pi / 6 >= angle or angle >= -5 * math.pi / 6 for angle in angles]):
                    orientation = "l"
                else:
                    # If we still don't have it (e.g. when there are 4+ equally-
                    # spaced bonds), just put everything in the top right for now
                    orientation = "tr"
            else:
                # There is an unbalanced side, so let's put the radical/charge data there
                angle = math.atan2(vector[1], vector[0])
                if 3 * math.pi / 4 <= angle or angle < -3 * math.pi / 4:
                    orientation = "l"
                elif -3 * math.pi / 4 <= angle < -1 * math.pi / 4:
                    orientation = "b"
                elif -1 * math.pi / 4 <= angle < 1 * math.pi / 4:
                    orientation = "r"
                else:
                    orientation = "t"

        cr.set_font_size(self.options["fontSizeNormal"])
        extents = cr.text_extents(heavy_atom)

        # (xi, yi) mark the center of the space in which to place the radicals and charges
        if orientation[0] == "l":
            xi = x - 3
            yi = y - extents[3] / 2
        elif orientation[0] == "b":
            xi = x + extents[0] + extents[2] / 2
            yi = y - extents[3] - 4
        elif orientation[0] == "r":
            xi = x + extents[0] + extents[2] + 4
            yi = y - extents[3] / 2
        elif orientation[0] == "t":
            xi = x + extents[0] + extents[2] / 2
            yi = y + 4

        # If we couldn't use one of the four sides, then offset the radical/charges
        # horizontally by a few pixels, in hope that this avoids overlap with an
        # existing bond
        if len(orientation) > 1:
            xi += 4

        # Get width and height
        cr.set_font_size(self.options["fontSizeSubscript"])
        width = 0.0
        height = 0.0
        if orientation[0] == "b" or orientation[0] == "t":
            if atom.radical_electrons > 0:
                width += atom.radical_electrons * 2 + (atom.radical_electrons - 1)
                height = atom.radical_electrons * 2
            text = ""
            if atom.radical_electrons > 0 and atom.charge != 0:
                width += 1
            if atom.charge == 1:
                text = "+"
            elif atom.charge > 1:
                text = "{0:d}+".format(atom.charge)
            elif atom.charge == -1:
                text = "\u2013"
            elif atom.charge < -1:
                text = "{0:d}\u2013".format(abs(atom.charge))
            if text != "":
                extents = cr.text_extents(text)
                width += extents[2] + 1
                height = extents[3]
        elif orientation[0] == "l" or orientation[0] == "r":
            if atom.radical_electrons > 0:
                height += atom.radical_electrons * 2 + (atom.radical_electrons - 1)
                width = atom.radical_electrons * 2
            text = ""
            if atom.radical_electrons > 0 and atom.charge != 0:
                height += 1
            if atom.charge == 1:
                text = "+"
            elif atom.charge > 1:
                text = "{0:d}+".format(atom.charge)
            elif atom.charge == -1:
                text = "\u2013"
            elif atom.charge < -1:
                text = "{0:d}\u2013".format(abs(atom.charge))
            if text != "":
                extents = cr.text_extents(text)
                height += extents[3] + 1
                width = extents[2]
        # Move (xi, yi) to top left corner of space in which to draw radicals and charges
        xi -= width / 2.0
        yi -= height / 2.0

        # Update bounding rectangle if necessary
        if width > 0 and height > 0:
            if xi < bounding_rect[0]:
                bounding_rect[0] = xi
            if yi < bounding_rect[1]:
                bounding_rect[1] = yi
            if xi + width > bounding_rect[2]:
                bounding_rect[2] = xi + width
            if yi + height > bounding_rect[3]:
                bounding_rect[3] = yi + height

        if orientation[0] == "b" or orientation[0] == "t":
            # Draw radical electrons first
            for i in range(atom.radical_electrons):
                cr.new_sub_path()
                cr.arc(xi + 3 * i + 1, yi + height / 2, 1, 0, 2 * math.pi)
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.fill()
            if atom.radical_electrons > 0:
                xi += atom.radical_electrons * 2 + (atom.radical_electrons - 1) + 1
            # Draw charges second
            text = ""
            if atom.charge == 1:
                text = "+"
            elif atom.charge > 1:
                text = "{0:d}+".format(atom.charge)
            elif atom.charge == -1:
                text = "\u2013"
            elif atom.charge < -1:
                text = "{0:d}\u2013".format(abs(atom.charge))
            if text != "":
                extents = cr.text_extents(text)
                cr.move_to(xi, yi - extents[1])
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.show_text(text)

            # Draw lone electron pairs
            # Draw them for nitrogen containing molecules only
            if draw_lone_pairs:
                for i in range(atom.lone_pairs):
                    cr.new_sub_path()
                    if i == 0:
                        x1lp = x - 2
                        y1lp = y - 8
                        x2lp = x + 2
                        y2lp = y - 12
                    elif i == 1:
                        x1lp = x + 12
                        y1lp = y - 8
                        x2lp = x + 8
                        y2lp = y - 12
                    elif i == 2:
                        x1lp = x - 2
                        y1lp = y - 1
                        x2lp = x + 2
                        y2lp = y + 3
                    self._draw_line(cr, x1lp, y1lp, x2lp, y2lp)

        elif orientation[0] == "l" or orientation[0] == "r":
            # Draw charges first
            text = ""
            if atom.charge == 1:
                text = "+"
            elif atom.charge > 1:
                text = "{0:d}+".format(atom.charge)
            elif atom.charge == -1:
                text = "\u2013"
            elif atom.charge < -1:
                text = "{0:d}\u2013".format(abs(atom.charge))
            if text != "":
                extents = cr.text_extents(text)
                cr.move_to(xi - extents[2] / 2, yi - extents[1])
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.show_text(text)
            if atom.charge != 0:
                yi += extents[3] + 1
            # Draw radical electrons second
            for i in range(atom.radical_electrons):
                cr.new_sub_path()
                cr.arc(xi + width / 2, yi + 3 * i + 1, 1, 0, 2 * math.pi)
                cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                cr.fill()
            # Draw lone electron pairs
            # Draw them for nitrogen species only
            if draw_lone_pairs:
                for i in range(atom.lone_pairs):
                    cr.new_sub_path()
                    if i == 0:
                        x1lp = x - 2
                        y1lp = y - 8
                        x2lp = x + 2
                        y2lp = y - 12
                    elif i == 1:
                        x1lp = x + 12
                        y1lp = y - 8
                        x2lp = x + 8
                        y2lp = y - 12
                    elif i == 2:
                        x1lp = x - 2
                        y1lp = y - 1
                        x2lp = x + 2
                        y2lp = y + 3
                    self._draw_line(cr, x1lp, y1lp, x2lp, y2lp)

        # Update bounding rect to ensure atoms are included
        if bounding_rect[0] < self.left:
            self.left = bounding_rect[0]
        if bounding_rect[1] < self.top:
            self.top = bounding_rect[1]
        if bounding_rect[2] > self.right:
            self.right = bounding_rect[2]
        if bounding_rect[3] > self.bottom:
            self.bottom = bounding_rect[3]

    def _make_single_bonds(self):
        """This method converts all bonds to single bonds and then returns
        a dictionary of Bond object keys with the old bond order as a value"""
        dictionary = {}
        for atom1 in self.molecule.atoms:
            for atom2, bond in atom1.bonds.items():
                if not bond.is_single():
                    dictionary[bond] = bond.get_order_num()
                    bond.set_order_num(1)
        return dictionary

    def _replace_bonds(self, bond_order_dictionary):
        """
        Sets the bond order in self.molecule equal to the orders in bond_order_dictionary
        which is obtained from _make_single_bonds().
        """
        for bond, order in bond_order_dictionary.items():
            bond.set_order_num(order)


################################################################################


class ReactionDrawer(object):
    """
    This class provides functionality for drawing chemical reactions using the
    skeletal formula of each reactant and product molecule via the Cairo 2D
    graphics engine. The most common use case is simply::

        ReactionDrawer().draw(reaction, file_format='png', path='reaction.png')

    where ``reaction`` is the :class:`Reaction` object to draw. You can also
    pass a dict of options to the constructor to affect how the molecules are
    drawn.
    """

    def __init__(self, options=None):
        self.options = MoleculeDrawer().options.copy()
        self.options.update(
            {
                "arrowLength": 36,
            }
        )
        if options:
            self.options.update(options)

    def draw(self, reaction, file_format, path=None):
        """
        Draw the given `reaction` using the given image `file_format` - pdf, svg,
        ps, or png. If `path` is given, the drawing is saved to that location
        on disk.

        This function returns the Cairo surface and context used to create the
        drawing, as well as a bounding box for the molecule being drawn as the
        tuple (`left`, `top`, `width`, `height`).
        """
        # The Cairo 2D graphics library (and its Python wrapper) is required for
        # the reaction drawing algorithm
        if cairo is None:
            logging.error("Cairo not found; molecule will not be drawn.")
            return

        font_family = self.options["fontFamily"]
        font_size_normal = self.options["fontSizeNormal"]

        # First draw each of the reactants and products
        reactants, products = [], []
        for reactant in reaction.reactants:
            if isinstance(reactant, Molecule):
                molecule = reactant
            elif hasattr(reactant, "molecule"):
                molecule = reactant.molecule[0]
            else:
                raise TypeError("Expected Molecule or Species object, not {0}".format(reactant.__class__.__name__))
            reactants.append(MoleculeDrawer().draw(molecule, file_format))
        for product in reaction.products:
            if isinstance(product, Molecule):
                molecule = product
            elif hasattr(product, "molecule"):
                molecule = product.molecule[0]
            else:
                raise TypeError("Expected Molecule or Species object, not {0}".format(product.__class__.__name__))
            products.append(MoleculeDrawer().draw(molecule, file_format))

        # Next determine size required for surface
        rxn_width = rxn_height = rxn_top = 0
        for surface, cr, rect in reactants + products:
            left, top, width, height = rect
            rxn_width += width
            if height > rxn_height:
                rxn_height = height
            if height + top > rxn_top:
                rxn_top = height + top

        rxn_top = 0.5 * rxn_height - rxn_top

        # Also include '+' and reaction arrow in width
        cr.set_font_size(font_size_normal)
        plus_extents = cr.text_extents(" + ")
        arrow_width = self.options["arrowLength"]
        rxn_width += (len(reactants) - 1) * plus_extents[4] + arrow_width + (len(products) - 1) * plus_extents[4]

        # Now make the surface for the reaction and render each molecule on it
        rxn_surface = create_new_surface(file_format, path, width=rxn_width, height=rxn_height)
        rxn_cr = cairo.Context(rxn_surface)

        # Draw white background
        rxn_cr.set_source_rgba(1.0, 1.0, 1.0, 1.0)
        rxn_cr.paint()

        # Draw reactants
        rxn_x = 0.0
        rxn_y = 0.0
        for index, reactant in enumerate(reactants):
            surface, cr, rect = reactant
            left, top, width, height = rect
            if index > 0:
                # Draw the "+" between the reactants
                rxn_cr.save()
                rxn_cr.set_font_size(font_size_normal)
                rxn_y = rxn_top + 0.5 * (rxn_height - plus_extents[3])
                rxn_cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                rxn_cr.move_to(rxn_x, rxn_y - plus_extents[1])
                rxn_cr.show_text(" + ")
                rxn_cr.restore()
                rxn_x += plus_extents[4]
            # Draw the reactant
            rxn_y = top + rxn_top + 0.5 * rxn_height
            if rxn_y < 0:
                rxn_y = 0
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
                rxn_cr.set_font_size(font_size_normal)
                rxn_y = rxn_top + 0.5 * (rxn_height - plus_extents[3])
                rxn_cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
                rxn_cr.move_to(rxn_x, rxn_y - plus_extents[1])
                rxn_cr.show_text(" + ")
                rxn_cr.restore()
                rxn_x += plus_extents[4]
            # Draw the product
            rxn_y = top + rxn_top + 0.5 * rxn_height
            if rxn_y < 0:
                rxn_y = 0
            rxn_cr.save()
            rxn_cr.set_source_surface(surface, rxn_x, rxn_y)
            rxn_cr.paint()
            rxn_cr.restore()
            rxn_x += width

        # Finish Cairo drawing
        if file_format == "png":
            rxn_surface.write_to_png(path)
        else:
            rxn_surface.finish()
